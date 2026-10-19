# Prototype of the correction c = w − v_s on the new core, as method overrides
# evaluated into the package's modules. Planning evidence for
# docs/Mesh_Framework_Phase2_Plan.md (Appendix A); nothing here is package code.
# Run only in a dedicated Julia session (see README.md): the overrides persist
# in the session and must never share one with the package's tests.
ENV["PYROLYSIS_CORE"] = "new";
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
using Pyrolysis, Statistics, LinearAlgebra, SparseArrays, Printf
using OrdinaryDiffEq: DiscreteCallback
const Cl = Pyrolysis.Closures;
const K = Pyrolysis.Kernels;
const St = Pyrolysis.States
const Ig = Pyrolysis.Integration;
const CJ = Pyrolysis.ColoredJacobian;
const As = Pyrolysis.Assembly
const De = Pyrolysis.Deactivation;
const Me = Pyrolysis.Meshes
const REPO = normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(Pyrolysis.Closures, :CORR)
    @eval Pyrolysis.Closures begin
        mutable struct Corr
            on::Bool
            τ::Float64
            α::Float64
            β::Float64
            r_hi::Float64
            r_lo::Float64
            lagrangian::Vector{Int}
            h0::Vector{Float64}
            solid::Vector{Bool}
            ρ_ref::Float64
            mode::Symbol   # :or (smooth or of the two cells' activations) or :and (product)
        end
        const CORR =
            Corr(false, 1.0, 1.0, 0.0, 0.3, 0.1, Int[], Float64[], Bool[], 1.0, :or)
        @inline _sstep(x) = x <= 0 ? zero(x) : (x >= 1 ? one(x) : x * x * (3 - 2x))
        @inline function _act(h, h0)
            r = h / h0
            r <= 1e-9 && return one(r)
            return _sstep((log(CORR.r_hi) - log(r)) / (log(CORR.r_hi) - log(CORR.r_lo)))
        end
        @inline function _monitor(ds, i)
            CORR.β == 0 && return CORR.α + zero(eltype(ds.V_eval))
            ms = zero(eltype(ds.m))
            @inbounds for j in eachindex(CORR.solid)
                CORR.solid[j] && (ms += ds.m[j, i])
            end
            return CORR.α + CORR.β * (ms / ds.V_eval[i]) / CORR.ρ_ref
        end
        # The correction velocity at node `nd` of a 1D mesh (cell i between nodes i, i+1).
        @inline function corr_c(ds, topo, nd::Int)
            T = eltype(ds.x)
            CORR.on || return zero(T)
            N = n_cells(topo)
            (nd == 1 || nd == N + 1) && return zero(T)
            nd in CORR.lagrangian && return zero(T)
            ib, ia = nd - 1, nd
            (ds.active[ib] && ds.active[ia]) || return zero(T)
            hb = ds.x[nd] - ds.x[nd-1]
            ha = ds.x[nd+1] - ds.x[nd]
            ab = _act(hb, CORR.h0[ib])
            aa = _act(ha, CORR.h0[ia])
            a = CORR.mode === :or ? 1 - (1 - ab) * (1 - aa) : ab * aa
            return a * (_monitor(ds, ia) * ha - _monitor(ds, ib) * hb) / CORR.τ
        end
    end
end

@eval Pyrolysis.Closures begin
    @inline function _velocity_rows!(c, cl::BalanceRowClosure, ::Nothing, ds, topo, t)
        g = ds.geometry
        θ_A = ds.θ_A[]
        c.r_velocity[cl.pinned] = ds.w[cl.pinned]
        @inbounds for k in cl.order
            s = zero(eltype(c.r_velocity))
            for f in topo.cells[k].faces
                σ = topo.faces[f].cells[1] == k ? 1 : -1
                s +=
                    σ *
                    g.face_area[f] *
                    g.face_normal[f][1] *
                    corr_c(ds, topo, topo.faces[f].nodes[1])
            end
            c.r_velocity[cl.row_node[k]] =
                _sweep_rate(g, topo, k, ds.w) - (c.θ[k] - θ_A) * ds.V_eval[k] - s
        end
        return nothing
    end
    function _consistent_velocities!(w, cl::BalanceRowClosure, ::Nothing, c, ds, topo, t)
        g = ds.geometry
        θ_A = ds.θ_A[]
        w[cl.pinned] = 0
        for k in cl.order
            row = cl.row_node[k]
            if !ds.active[k]
                w[row] = 0
                continue
            end
            nodes = topo.cells[k].nodes
            other = nodes[1] == row ? nodes[2] : nodes[1]
            s = zero(eltype(w))
            for f in topo.cells[k].faces
                σ = topo.faces[f].cells[1] == k ? 1 : -1
                s +=
                    σ *
                    g.face_area[f] *
                    g.face_normal[f][1] *
                    corr_c(ds, topo, topo.faces[f].nodes[1])
            end
            rate =
                (c.θ[k] - θ_A) * ds.V_eval[k] + s -
                _coefficient(g, topo, k, other) * w[other]
            w[row] = rate / _coefficient(g, topo, k, row)
        end
        return w
    end
end

@eval Pyrolysis.Kernels function contribute!(
    c::Contributions{NC},
    kernel::RelativeFluxKernel{L,NC},
    ds::DerivedState{NC},
    topo::MeshTopology{1},
    material,
    t,
) where {L,NC}
    g = ds.geometry
    @inbounds for f = 1:n_faces(topo)
        u =
            -Main.Pyrolysis.Closures.corr_c(ds, topo, topo.faces[f].nodes[1]) *
            g.face_normal[f][1]
        s = _stencil(kernel, g, topo, f, u)
        U = s[1]
        T_f = _face_value(i -> ds.T[i], s..., kernel.ε[2, U])
        h = component_enthalpies(kernel.laws, T_f)
        q = zero(eltype(c.F_energy))
        F_solid = zero(q)
        for j = 1:NC
            F = u * _face_value(i -> ds.ξ[j, i], s..., kernel.ε[1, U])
            c.F_mass[j, f] += F
            q += F * h[j]
            kernel.solid[j] && (F_solid += F)
        end
        c.F_energy[f] += q
        kernel.history || continue
        ε = kernel.ε[3, U]
        g_f = _face_value(i -> ds.G[i] / _solid_mass(kernel, ds, i), s..., ε)
        μ_f = _face_value(i -> ds.M[i] / _solid_mass(kernel, ds, i), s..., ε)
        c.F_carried[1, f] += F_solid * g_f
        c.F_carried[2, f] += F_solid * μ_f
    end
    return c
end

@eval Pyrolysis.ColoredJacobian function jacobian_pattern(def::CoreDefinition)
    reconstruction = find_kernel(RelativeFluxKernel, def.kernels) !== nothing
    P = jacobian_pattern(
        def.layout,
        def.setup.topology;
        closure = def.closure,
        reconstruction,
    )
    Main.Pyrolysis.Closures.CORR.on || return P
    cl = def.closure
    cl isa BalanceRowClosure || return P
    cl.override === nothing || return P
    layout = def.layout
    N = n_cells(layout)
    rows = Int[]
    cols = Int[]
    for k in cl.order
        r = velocity_index(layout, cl.row_node[k])
        for nd = max(1, k-1):min(N+1, k+2)
            push!(rows, r)
            push!(cols, position_index(layout, nd))
        end
        for kk = max(1, k-1):min(N, k+1), c in cell_rows(layout, kk)
            push!(rows, r)
            push!(cols, c)
        end
    end
    return SparseMatrixCSC{Bool,Int}(
        P .| sparse(rows, cols, trues(length(rows)), size(P)...),
    )
end

# Block deactivation: the contiguous run of below-threshold cells ending at the
# top live cell merges, cell by cell, into the first live cell beneath it (or
# into the lowest cell of the column, which survives).
isdefined(Pyrolysis.Deactivation, :BLOCK) ||
    @eval Pyrolysis.Deactivation const BLOCK = Ref(true)
isdefined(Pyrolysis.Deactivation, :BAND) ||
    @eval Pyrolysis.Deactivation const BAND = Ref(1.0)
isdefined(Pyrolysis.Deactivation, :BAND_LOG) ||
    @eval Pyrolysis.Deactivation const BAND_LOG = Tuple{Float64,Int,Float64}[]
@eval Pyrolysis.Deactivation function deactivate!(
    u::AbstractVector,
    d::CellDeactivation,
    def::CoreDefinition,
    t::Real,
)
    layout = def.layout
    if !BLOCK[]
        while true
            top = top_live_cell(layout)
            (top == 0 || !below_threshold(d, layout, u, top)) && break
            if top == 1
                _lose!(u, d, def, top, Float64(t))
                d.outcome[] = :burned_out
                return :burned_out
            end
            _transfer!(u, d, d.policy, def, top, top - 1, Float64(t))
        end
        if any_below(d, layout, u)
            d.outcome[] = :interior
            return :interior
        end
        return :transferred
    end
    top = top_live_cell(layout)
    # the near-threshold band: the run from the top of cells below BAND·ε_m
    inband(i) =
        layout.mask[energy_index(layout, i)] &&
        d.reference[i] > 0 &&
        condensed_mass(d, layout, u, i) < BAND[] * d.threshold * d.reference[i]
    lo = top
    while lo >= 1 && inband(lo)
        lo -= 1
    end
    lo += 1
    if lo <= top && any(i -> below_threshold(d, layout, u, i), lo:top)
        fr = maximum(condensed_mass(d, layout, u, i) / d.reference[i] for i = lo:top)
        push!(BAND_LOG, (Float64(t), top - lo + 1, fr))
        if lo == 1 && top == 1
            _lose!(u, d, def, 1, Float64(t))
            d.outcome[] = :burned_out
            return :burned_out
        end
        for i = top:-1:max(lo, 2)
            _transfer!(u, d, d.policy, def, i, i - 1, Float64(t))
        end
    end
    if any_below(d, layout, u)
        d.outcome[] = :interior
        return :interior
    end
    return :transferred
end

include(joinpath(REPO, "test", "verification", "casedefs.jl"))
const CASES = load_verification_cases()
const CORR = Pyrolysis.Closures.CORR
const BLOCK = Pyrolysis.Deactivation.BLOCK
const BAND = Pyrolysis.Deactivation.BAND
const BAND_LOG = Pyrolysis.Deactivation.BAND_LOG

function build_vbo(
    n_pmma;
    τ = 1.0,
    on = true,
    α = 1.0,
    β = 0.0,
    r_hi = 0.3,
    r_lo = 0.1,
    mode = :or,
    reltol = 1e-6,
    abstol = 1e-8,
    threshold = 1e-2,
    builder = nothing,
    case_kw = (;),
)
    prob = builder === nothing ? Base.invokelatest(_vbo_problem, n_pmma) : builder(n_pmma)
    CORR.on = false
    cp = Ig.setup(
        prob;
        use_ale = true,
        deactivation = true,
        abstol,
        reltol,
        threshold,
        time_scale = 1e-3,
        case_kw...,
    )
    def = cp.def
    layout = def.layout
    setup = def.setup
    topo = setup.topology
    g0 = Me.MeshGeometry{Float64}(topo)
    Me.geometry!(g0, setup.reference.positions, topo.cross_section, topo)
    rk = K.RelativeFluxKernel(
        topo,
        g0,
        setup.material,
        setup.laws,
        setup.reference;
        history = false,
        temperature_scale = abstol,
        history_scale = reltol * threshold,
    )
    CORR.h0 = diff(setup.reference.positions)
    ξ1 = [
        prob.ξ_initial[1](z) for z in
        (setup.reference.positions[1:(end-1)] .+ setup.reference.positions[2:end]) ./ 2
    ]
    # layer interface: the face where the first component appears
    iface = findfirst(i -> (ξ1[i] > 0) != (ξ1[i+1] > 0), 1:(length(ξ1)-1))
    CORR.lagrangian = iface === nothing ? Int[] : [iface + 1]
    CORR.solid = [!Pyrolysis.is_gas(c) for c in setup.material.components]
    CORR.ρ_ref = maximum(
        sum(cp.m0[j, i] for j in eachindex(CORR.solid) if CORR.solid[j]) /
        g0.cell_volume[i] for i = 1:Me.n_cells(topo)
    )
    CORR.τ = τ
    CORR.α = α
    CORR.β = β
    CORR.r_hi = r_hi
    CORR.r_lo = r_lo
    CORR.mode = mode
    CORR.on = on
    kernels = on ? (def.kernels..., rk) : def.kernels
    def2 = As.CoreDefinition(layout, setup, kernels, def.closure, def.sources)
    ws2 = As.CoreWorkspace(def2, cp.T0)
    u0 = copy(cp.u0)
    Ig.consistent_state!(u0, def2, ws2, prob.tspan[1])
    cj2 = CJ.ColoredJacobianCache(def2, ws2, u0, prob.tspan[1])
    cp2 = Ig.CoreProblem(
        cp.problem,
        def2,
        ws2,
        u0,
        cp.reltol,
        cp.abstol,
        cj2,
        cp.T0,
        cp.m0,
        cp.deactivation,
    )
    return cp2, prob
end

# recorder: smallest live-cell thickness over the run, smallest step
function recorder(cp, hmin, dtmin, cmax)
    layout = cp.def.layout
    return DiscreteCallback(
        (u, t, i) -> true,
        integ -> begin
            u = integ.u
            hm = Inf
            for i = 1:St.n_cells(layout)
                layout.mask[St.energy_index(layout, i)] || continue
                hm = min(
                    hm,
                    u[St.position_index(layout, i + 1)] - u[St.position_index(layout, i)],
                )
            end
            hmin[] = min(hmin[], hm)
            dtmin[] = min(dtmin[], integ.t - integ.tprev)
            nothing
        end;
        save_positions = (false, false),
    )
end

function run_vbo(n_pmma; min_thickness = 1e-3, nrep = 2, kw...)
    cp2, prob = build_vbo(n_pmma; kw...)
    hmin = Ref(Inf)
    dtmin = Ref(Inf)
    cmax = Ref(0.0)
    sol = nothing
    times = Float64[]
    for r = 1:nrep
        hmin[] = Inf
        dtmin[] = Inf
        empty!(BAND_LOG)
        t = @elapsed sol = Ig.integrate(
            cp2;
            callback = recorder(cp2, hmin, dtmin, cmax),
            min_thickness,
        )
        push!(times, t)
        println(
            "      solve $(r): retcode=$(sol.retcode) outcome=$(sol.outcome) events=$(length(sol.events)) invocations=$(length(unique(e.t for e in sol.events))) t_end=$(round(sol.t[end], digits=2)) hmin=$(hmin[])",
        )
    end
    return (;
        sol,
        wall = minimum(times),
        times,
        hmin = hmin[],
        dtmin = dtmin[],
        colors = cp2.jacobian.n_colors,
        cp = cp2,
    )
end

function report(tag, r)
    sol = r.sol
    st = sol.statistics
    ev = sol.events
    line = @sprintf(
        "[%s] retcode=%s outcome=%s t_end=%.2f events=%d wall=%.3f s (%s) steps=%d rej=%d jac=%d nf=%d domain_rej=%d cres=%.2e hmin=%.3e m dtmin=%.2e colors=%d",
        tag,
        sol.retcode,
        sol.outcome,
        sol.t[end],
        length(ev),
        r.wall,
        join(round.(r.times, digits = 3), ","),
        st.accepted,
        st.rejected,
        st.jacobians,
        st.rhs_evaluations,
        st.domain_rejections,
        st.largest_constraint_residual,
        r.hmin,
        r.dtmin,
        r.colors
    )
    println(line)
    if !isempty(ev)
        strs = ["($(round(e.t, digits=2)), $(e.cell)→$(e.receiver))" for e in ev]
        println(
            "    invocations (distinct event times) = ",
            length(unique(e.t for e in ev)),
            "; events: ",
            join(first(strs, 12), " "),
            " … ",
            join(last(strs, 12), " "),
        )
    end
    isempty(BAND_LOG) || println(
        "    band merges (t, cells, largest mass fraction): ",
        join(
            [
                "($(round(b[1], digits=2)), $(b[2]), $(round(b[3], sigdigits=3)))" for
                b in BAND_LOG
            ],
            " ",
        ),
    )
    layout = r.cp.def.layout
    u = sol.raw_solution.u[end]
    d = r.cp.deactivation
    De.replay_events!(d, layout, sol.t[end])
    below = [i for i = 1:St.n_cells(layout) if De.below_threshold(d, layout, u, i)]
    println(
        "    live cells at the end = ",
        De.top_live_cell(layout),
        "; cells below threshold at the end: ",
        below,
    )
    return nothing
end

function vbo_qois(sol)
    m = Base.invokelatest(_vbo_measures, sol)
    @printf(
        "    peak=%.3f (ref %.2f, %.2f%%) t_peak=%.1f (ref %.0f) t_bo=%.1f (ref %.1f) mass_loss=%.3f (exact %.3f) ledger=%.2e L_end=%.6f\n",
        m.peak,
        VBO_PEAK_REF,
        100 * (m.peak / VBO_PEAK_REF - 1),
        m.t_peak,
        VBO_T_PEAK_REF,
        m.t_bo,
        VBO_T_BO_REF,
        m.mass_loss,
        VBO_MASS_LOSS_EXACT,
        m.ledger,
        m.L_end
    )
    return m
end
println("proto_common loaded")
