# Shared pieces of the Phase 0 Jacobian baselines (items 4 and 5 of the
# mesh-framework plan): deck construction copied from the gitignored
# performance-audit harness, the union of the production sparsity patterns,
# the lifted-pattern model, and a greedy distance-1 column coloring.
#
# Loaded by jacobian_colors.jl, ad_overhead.jl and lifted_containment.jl; run
# those from the repo root with
# `julia --project=test --startup-file=no --threads=1 <script>`.

using Pyrolysis, SparseArrays, LinearAlgebra, Statistics, Printf
using Pyrolysis.Jacobian:
    DirectSolve, ApproxSolve, JacobianSpec, Structured, compute_jacobian!, build_sparsity
using Pyrolysis.Problem: build_ad_workspace, chi_mode, radiation_mode, WithChi, BeerLambert
using Pyrolysis.Problem: mesh_mode, ALE
using Pyrolysis.Residual: residual!

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULTS = joinpath(ROOT, "scripts", "mesh_phase0", "results")
mkpath(RESULTS)

# ----------------------------------------------------------------------------
# Decks (copied from scripts/performance_audit/audit.jl `problem_for`)
# ----------------------------------------------------------------------------

function problem_for(
    kind,
    n;
    t_end = 80.0,
    outputs = 3,
    solid_λ = nothing,
    lateral = nothing,
    flux = 50e3,
)
    wood = kind in (:ale_darcy, :eulerian_darcy)
    mat = if wood
        Material(
            name = :audit_wood,
            components = (
                SolidComponent(
                    :wood,
                    ρ = 500.0,
                    ρ_intrinsic = 1500.0,
                    c = 1500.0,
                    k = 0.15,
                    ε = 0.9,
                    α = 5.0,
                    κ = 1e-14,
                    λ = something(solid_λ, 1e-5),
                ),
                GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5),
            ),
            reactions = (Reaction(:pyrolysis, 1=>2, A = 1e10, E = 150e3, h = 300e3),),
            lateral_shrinkage_law = lateral,
        )
    else
        Material(
            name = :audit_pmma,
            components = (
                SolidComponent(
                    :pmma,
                    ρ = 1190.0,
                    c = 1420.0,
                    k = 0.21,
                    ε = 0.86,
                    λ = something(solid_λ, 0.0),
                ),
                GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
            ),
            reactions = (Reaction(:pyrolysis, 1=>2, A = 8.5e12, E = 188e3, h = 870e3),),
        )
    end
    thermal = (
        top = RadiativeFluxBC(flux = flux) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    )
    mass =
        wood ?
        (
            top = DiffusiveMassBC(0.05)+DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ) : (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC())
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(; thermal, mass),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = outputs)),
    )
    return problem,
    (;
        use_ale = kind in (:ale, :ale_darcy),
        radiation_model = wood ? BEER_LAMBERT : SURFACE_ABSORPTION,
        progress = false,
    )
end

"""
	build_deck(name, n) -> (def, ws, u, NC, n_cells)

`name` is one of :eulerian, :ale, :ale_darcy, :ale_darcy_withchi, :vbo01.
The V-BO-01 deck ignores `n` and uses its CI cell count.
"""
function build_deck(name::Symbol, n::Int)
    if name === :vbo01
        case = Base.invokelatest(_verification_case_vbo01)
        prob = case.build(case.ci.n_cells)
        kw = case.solve_kwargs
        def, ws, u = Pyrolysis.initialize_state(
            prob;
            use_ale = kw.use_ale,
            radiation_model = kw.radiation_model,
            radiation_intensity = 50e3,
        )
    else
        kind = name === :ale_darcy_withchi ? :ale_darcy : name
        lateral = name === :ale_darcy_withchi ? (χ -> 1.0 - 0.2χ) : nothing
        prob, kw = problem_for(kind, n; lateral)
        def, ws, u = Pyrolysis.initialize_state(
            prob;
            use_ale = kw.use_ale,
            radiation_model = kw.radiation_model,
            radiation_intensity = 50e3,
        )
    end
    layout = ws.layout
    NC = Pyrolysis.Materials.n_components(layout)
    return def, ws, u, NC, layout.n_cells
end

# Synthetic mixed, reacting state (audit `kernels(...; hot = true)`); keeps the
# ALE node positions. Assumes NC = 2 with the gas as component 2.
function hot_state!(u, n)
    for i = 1:n
        u[i] = 450.0 + 200.0*i/n
        u[n+i] *= 0.65 + 0.3*i/n
        u[2n+i] = 0.02 + 0.08*i/n
    end
    return u
end

# ----------------------------------------------------------------------------
# Current production patterns
# ----------------------------------------------------------------------------

"""
	current_pattern(ws, profile) -> SparseMatrixCSC{Bool,Int}

`S_phys ∪ S_geom` for `profile ∈ (:full, :ale_reduced)`, exactly as `solve`
forms `jac_prototype`.
"""
function current_pattern(ws, profile::Symbol)
    S_phys, S_geom = build_sparsity(ws.layout, ws.layout.mode; profile)
    return S_phys .| S_geom
end

# ----------------------------------------------------------------------------
# Lifted-pattern model (plan §3.7): every nonlocal coupling is a local
# algebraic row. Interleaved cell-major ordering:
#   cell i: T_i, ξ_{1..NC,i}, [G_i, M_i], z_i, w_i, [I_i]
#   node n+1: z_{n+1}, w_{n+1}, [I_{n+1}]
#   [A, θ_A]  (global scalars; their rows are dense and are NOT part of the
#              colored pattern — they carry only a structural diagonal here;
#              their columns stay in the cell and balance rows)
#
# χ is carried as the extensive pair (G, M) of the design's §3.2, two cell
# unknowns coupled like every other cell unknown.
#
# A cell row (T_i, ξ_{j,i}, G_i, M_i) reads, with h = (s − 1)/2 for the
# species stencil s ∈ (3, 5):
#   * every unknown of cells i−h..i+h. With s = 3 (label c0: no relative
#     solid flux, c = 0) the two-point face fluxes read the neighbours'
#     states; with s = 5 (label c_active) the second-order limited
#     reconstruction of the relative solid flux reads the second neighbours
#     too. The energy row carries the same face mass fluxes (the enthalpy
#     flux Σ_j F_{j,f} h_j(T_f) uses the reconstructed state), so every cell
#     row takes the same stencil.
#   * node positions z_{i−h}..z_{i+h+1}: a face flux divides by the distance
#     between the neighbouring cell centres, each a function of two nodes,
#     and a neighbour's concentration ξ = m/V reads that neighbour's volume;
#     with s = 5 the limiter slopes read the second neighbours' centres too.
#   * node velocities w_i, w_{i+1} (relative flux through the two faces)
#     and, under Beer–Lambert, the face intensities I_i, I_{i+1};
#   * the global cross-section A under WithChi (volumes V = A·Δz).
# Radial geometry has identical locality (areas and volumes are local
# functions of the node positions), so no separate radial pattern is built.
# ----------------------------------------------------------------------------

const STENCILS = (3, 5)
stencil_label(s::Int) = s == 3 ? "c0" : "c_active"

struct LiftedLayout
    n::Int
    NC::Int
    chi::Bool
    bl::Bool
    stencil::Int      # cells read by a cell row (3 or 5)
    per_cell::Int
    N::Int
end

function LiftedLayout(n, NC; chi::Bool, bl::Bool, species_stencil::Int = 3)
    species_stencil in STENCILS ||
        throw(ArgumentError("species_stencil must be 3 or 5, got $species_stencil"))
    per_cell = 1 + NC + (chi ? 2 : 0) + 2 + (bl ? 1 : 0)
    N = n * per_cell + 2 + (bl ? 1 : 0) + (chi ? 2 : 0)
    return LiftedLayout(n, NC, chi, bl, species_stencil, per_cell, N)
end

n_cell_unknowns(L::LiftedLayout) = 1 + L.NC + (L.chi ? 2 : 0)
_cell0(L::LiftedLayout, i) = (i - 1) * L.per_cell
lT(L, i) = _cell0(L, i) + 1
lξ(L, j, i) = _cell0(L, i) + 1 + j
lG(L, i) = _cell0(L, i) + 2 + L.NC
lM(L, i) = _cell0(L, i) + 3 + L.NC
lz(L, k) = k <= L.n ? _cell0(L, k) + n_cell_unknowns(L) + 1 : L.n * L.per_cell + 1
lw(L, k) = lz(L, k) + 1
lI(L, f) = lz(L, f) + 2
lA(L) = L.n * L.per_cell + 2 + (L.bl ? 1 : 0) + 1
lθA(L) = lA(L) + 1
ncols_local(L::LiftedLayout) = L.N - (L.chi ? 2 : 0)

# Every unknown of cell i (T, ξ, G, M).
function cell_unknowns(L::LiftedLayout, i)
    cols = Int[lT(L, i)]
    for j = 1:L.NC
        push!(cols, lξ(L, j, i))
    end
    L.chi && append!(cols, (lG(L, i), lM(L, i)))
    return cols
end

"Cells whose unknowns a cell row of cell `i` reads under the layout's stencil."
function cell_row_cells(L::LiftedLayout, i)
    h = (L.stencil - 1) ÷ 2
    return max(1, i-h):min(L.n, i+h)
end

"Nodes whose positions a cell row of cell `i` reads under the layout's stencil."
function cell_row_nodes(L::LiftedLayout, i)
    h = (L.stencil - 1) ÷ 2
    return max(1, i-h):min(L.n+1, i+h+1)
end

"""
	lifted_color_formula(NC, s; chi, bl) -> Int

Degree of an interior cell row of the model, the clique lower bound of its
column coloring: `s·(NC + 1 + 2χ) + (s + 1) + 2 + 2·bl + χ` (cell unknowns
of `s` cells, `s + 1` node positions, two node velocities, two face
intensities, the cross-section column).
"""
lifted_color_formula(NC, s; chi::Bool, bl::Bool) =
    s * (NC + 1 + 2chi) + (s + 1) + 2 + 2bl + chi

"""
	lifted_pattern(n, NC; chi, bl, species_stencil = 3) -> (S, L)

Model of the lifted Jacobian pattern (planar; radial is identical).
"""
function lifted_pattern(n, NC; chi::Bool, bl::Bool, species_stencil::Int = 3)
    L = LiftedLayout(n, NC; chi, bl, species_stencil)
    I = Int[]
    J = Int[]
    add!(r, c) = (push!(I, r); push!(J, c))
    # Cell rows: T_i, ξ_{j,i}, G_i, M_i.
    for i = 1:n
        rows = cell_unknowns(L, i)
        cols = Int[]
        for nb in cell_row_cells(L, i)
            append!(cols, cell_unknowns(L, nb))
        end
        for k in cell_row_nodes(L, i)
            push!(cols, lz(L, k))
        end
        append!(cols, (lw(L, i), lw(L, i+1)))
        bl && append!(cols, (lI(L, i), lI(L, i+1)))
        chi && push!(cols, lA(L))
        for r in rows, c in cols
            add!(r, c)
        end
    end
    # Node rows: dz_k/dt = w_k.
    for k = 1:(n+1)
        add!(lz(L, k), lz(L, k))
        add!(lz(L, k), lw(L, k))
    end
    # Velocity rows: node 1 pinned (w_1 = 0); node k+1 carries the volume
    # balance of cell k.
    add!(lw(L, 1), lw(L, 1))
    for k = 1:n
        r = lw(L, k+1)
        for c in (lw(L, k), lw(L, k+1), lz(L, k), lz(L, k+1))
            add!(r, c)
        end
        for c in cell_unknowns(L, k)
            add!(r, c)
        end
        chi && add!(r, lθA(L))
    end
    # Intensity rows: I_{n+1} = a_eff(T_n, ξ_n)·I₀ at the exposed face;
    # I_f = I_{f+1}·exp(−α_f Δz_f) below it.
    if bl
        r = lI(L, n+1)
        add!(r, r)
        for c in cell_unknowns(L, n)
            add!(r, c)
        end
        for f = 1:n
            r = lI(L, f)
            for c in (lI(L, f), lI(L, f+1), lz(L, f), lz(L, f+1))
                add!(r, c)
            end
            for c in cell_unknowns(L, f)
                add!(r, c)
            end
        end
    end
    # Global scalars: structural diagonal only (dense rows are assembled by
    # hand after the colored pass).
    if chi
        add!(lA(L), lA(L))
        add!(lθA(L), lθA(L))
    end
    S = sparse(I, J, trues(length(I)), L.N, L.N, |)
    return S, L
end

"""
	model_row_dependencies(S, L, row) -> (; cells, nodes, A)

Cell indices, node indices and the cross-section column a model row reads,
from the stored entries of `S`.
"""
function model_row_dependencies(S::SparseMatrixCSC, L::LiftedLayout, row::Int)
    cells = Int[]
    nodes = Int[]
    for i = 1:L.n
        any(c -> S[row, c], cell_unknowns(L, i)) && push!(cells, i)
    end
    for k = 1:(L.n+1)
        S[row, lz(L, k)] && push!(nodes, k)
    end
    return (; cells, nodes, A = L.chi && S[row, lA(L)])
end

# ----------------------------------------------------------------------------
# Pattern statistics and greedy distance-1 column coloring
# ----------------------------------------------------------------------------

"Max |row − col| over the stored entries of `S`, restricted to rows and columns ≤ `nloc`."
function pattern_bandwidth(S::SparseMatrixCSC, nloc::Int = size(S, 2))
    bw = 0
    rv = rowvals(S)
    for j = 1:min(size(S, 2), nloc)
        for p in nzrange(S, j)
            r = rv[p]
            r <= nloc && (bw = max(bw, abs(r - j)))
        end
    end
    return bw
end

function max_row_degree(S::SparseMatrixCSC)
    deg = zeros(Int, size(S, 1))
    for r in rowvals(S)
        deg[r] += 1
    end
    return maximum(deg)
end

# Columns conflicting with each column (distinct count), via the transpose.
function column_degrees(S::SparseMatrixCSC, St::SparseMatrixCSC)
    N = size(S, 2)
    deg = zeros(Int, N)
    stamp = zeros(Int, N)
    rv, rvt = rowvals(S), rowvals(St)
    for j = 1:N
        d = 0
        for p in nzrange(S, j)
            r = rv[p]
            for q in nzrange(St, r)
                c = rvt[q]
                if c != j && stamp[c] != j
                    stamp[c] = j
                    d += 1
                end
            end
        end
        deg[j] = d
    end
    return deg
end

"""
	greedy_column_coloring(S, St, order) -> (color::Vector{Int}, ncolors)

Distance-1 column coloring: two columns conflict when they share a stored row.
Columns are visited in `order`; each takes the smallest color unused by its
already-colored neighbours.
"""
function greedy_column_coloring(S::SparseMatrixCSC, St::SparseMatrixCSC, order)
    N = size(S, 2)
    color = zeros(Int, N)
    forbidden = zeros(Int, N + 1)
    rv, rvt = rowvals(S), rowvals(St)
    ncolors = 0
    for j in order
        for p in nzrange(S, j)
            r = rv[p]
            for q in nzrange(St, r)
                col = color[rvt[q]]
                col > 0 && (forbidden[col] = j)
            end
        end
        k = 1
        while forbidden[k] == j
            k += 1
        end
        color[j] = k
        ncolors = max(ncolors, k)
    end
    return color, ncolors
end

"""
	color_pattern(S) -> (; colors, coloring, ordering, max_row_degree, nnz)

Colors `S` with the natural and the largest-degree-first ordering and keeps
the smaller result. `max_row_degree` is the clique lower bound.
"""
function color_pattern(S::SparseMatrixCSC)
    St = sparse(transpose(S))
    N = size(S, 2)
    c_nat, k_nat = greedy_column_coloring(S, St, 1:N)
    deg = column_degrees(S, St)
    order = sortperm(deg; rev = true)
    c_ldf, k_ldf = greedy_column_coloring(S, St, order)
    if k_ldf < k_nat
        return (;
            colors = k_ldf,
            coloring = c_ldf,
            ordering = "largest_degree_first",
            max_row_degree = max_row_degree(S),
            nnz = nnz(S),
        )
    else
        return (;
            colors = k_nat,
            coloring = c_nat,
            ordering = "natural",
            max_row_degree = max_row_degree(S),
            nnz = nnz(S),
        )
    end
end

# ----------------------------------------------------------------------------
# Timing: median and min of `repeats` warmed samples, each a loop of `nevals`
# calls sized so one sample takes about `target_s`.
# ----------------------------------------------------------------------------
function time_op(f; repeats = 7, target_s = 0.05)
    f()
    f()
    t1 = @elapsed f()
    nevals = clamp(ceil(Int, target_s / max(t1, 1e-9)), 1, 100_000)
    samples = Float64[]
    for _ = 1:repeats
        t = @elapsed for _ = 1:nevals
            f()
        end
        push!(samples, t / nevals)
    end
    return (; median = median(samples), min = minimum(samples), nevals, repeats)
end

# CSV output: RFC 4180 — a field holding a comma, a double quote or a line
# break is quoted (quotes doubled); everything else is written bare.
function csv_field(x::AbstractString)
    needs = any(c -> c == ',' || c == '"' || c == '\n' || c == '\r', x)
    return needs ? '"' * replace(x, '"' => "\"\"") * '"' : String(x)
end
csv_field(x) = csv_field(string(x))

function csv_row(io, xs...)
    println(io, join((csv_field(x) for x in xs), ','))
    flush(io)
end
