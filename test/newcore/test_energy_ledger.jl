using Test
using ForwardDiff
using LinearAlgebra
using Pyrolysis
using Pyrolysis.Properties: guarded_harmonic_mean
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.Kernels: GasTransportKernel
using Pyrolysis.Integration: setup
using OrdinaryDiffEq: ODEFunction, ODEProblem
import OrdinaryDiffEq
using Pyrolysis.Solver: KenCarp4

# The char-layer column of V-CPL-06 with a reaction heat: a reactive band low
# in a hot, gas-permeable matrix gasifies, and the gas diffuses up through the
# hot layer to leave through a film at the top (Fick), picking up a few
# hundred kelvin on the way. Under Darcy transport the band sits against an
# adiabatic, sealed back and the pores start above the ambient pressure,
# falling toward the pressure vent at the top, so the gas flows up through
# every face throughout: the first-order upwind switch of the advective flux
# is a kink at zero velocity, where fixed steps at a converged Newton
# tolerance fail on either core, and a face through which the flow reverses
# would test that instead of the ledger.
function _column(n; darcy::Bool)
    L = 0.01
    matrix = SolidComponent(
        :matrix,
        ρ = 2000.0,
        c = (800.0, 0.4),
        k = 1.0,
        λ = 1e-4,
        κ = darcy ? 1e-13 : 0.0,
        ρ_intrinsic = 2600.0,
    )
    fuel = SolidComponent(
        :fuel,
        ρ = 2000.0,
        c = 1000.0,
        k = 1.0,
        λ = 1e-4,
        ρ_intrinsic = 2600.0,
    )
    gas = GaseousComponent(:gas, M = 0.05, c = (900.0, 0.3), k = 0.02, λ = 1e-4)
    rxn = Reaction(:gasify, 2 => 3; A = 9.3e3, E = 8.0e4, h = 3e5)
    mat = Material(name = :column, components = (matrix, fuel, gas), reactions = (rxn,))
    vent = darcy ? DirichletPressureBC(101325.0) : ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0)
    z_band = darcy ? 0.1 : 0.3
    band(z) = 100.0 * exp(-((z / L - z_band) / 0.08)^2)
    T0(z) = 550.0 + 3e4 * z
    P0(z) = 101325.0 * (1.3 - 0.3 * z / L)
    gas0(z) =
        darcy ?
        background_gas_concentration(mat, 3, (400.0, band(z), 0.0), T0(z); P = P0(z)) : 0.05
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = ConvectiveBC(h = 1000.0, T_inf = 900.0),
                bottom = darcy ? AdiabaticBC() : ConvectiveBC(h = 1000.0, T_inf = 300.0),
            ),
            mass = (top = vent, bottom = ImpermeableBC()),
        ),
        T_initial = T0,
        ξ_initial = [z -> 400.0, band, gas0],
        tspan = (0.0, 20.0),
    )
end

# The total-enthalpy rate the inventories imply: what enters through the
# boundary faces (conduction or the boundary condition, and the enthalpy the
# boundary mass fluxes carry), the reaction heat and the formation terms.
function _ledger_rate(def, c, ds)
    topo = def.setup.topology
    g = ds.geometry
    L = zero(eltype(c.F_energy))
    for f = 1:n_faces(topo)
        topo.faces[f].cells[2] == 0 && (L -= c.F_energy[f] * g.face_area[f])
    end
    for i = 1:n_cells(topo), j in axes(c.S_mass, 1)
        L += c.S_mass[j, i] * ds.h[j, i]
    end
    return L + sum(c.S_energy)
end

# The total enthalpy the inventories hold: the H inventories themselves in the
# inventory form, Σᵢ Σⱼ mᵢⱼ hⱼ(Tᵢ) in the temperature form.
_enthalpy(def, ws, u) = _enthalpy(energy_state(def.layout), def, ws, u)
_enthalpy(::EnthalpyState, def, ws, u) =
    sum(u[energy_index(def.layout, i)] for i = 1:n_cells(def.layout))
function _enthalpy(::TemperatureState, def, ws, u)
    ds = ws.derived
    derive!(ds, u, def.layout, def.setup)
    return sum(ds.m[j, i] * ds.h[j, i] for j in axes(ds.m, 1), i in axes(ds.m, 2))
end

@testset "energy row term by term against the ledger: $(darcy ? "Darcy" : "Fick"), $energy form" for darcy in
                                                                                                     (
        false,
        true,
    ),
    energy in (:temperature, :enthalpy)

    cp = setup(_column(10; darcy); energy_state = energy)
    def, ws = cp.def, cp.ws
    layout, topo = def.layout, def.setup.topology
    laws = def.setup.laws
    kernel = only(k for k in def.kernels if k isa GasTransportKernel)
    u = copy(cp.u0)
    du = zeros(length(u))
    residual!(du, u, def, 3.0, ws)
    ds, c = buffers(ws, def, Float64)
    g = ds.geometry
    @test maximum(abs, c.F_mass) > 0

    # Every face flux is conduction (Fourier between cells, the condition on the
    # boundary) plus the enthalpy its mass flux carries at the face temperature.
    carried = zeros(n_faces(topo))
    conducted = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        if nb == 0
            T_f = c.T_face[f]
        else
            d_o = center_to_face_distance(g, topo, f, o)
            d_n = center_to_face_distance(g, topo, f, nb)
            T_f = (d_n * ds.T[o] + d_o * ds.T[nb]) / (d_o + d_n)
            k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
            @test c.F_energy[f] - sum(c.F_mass[:, f] .* component_enthalpies(laws, T_f)) ≈
                  -k_f * (ds.T[nb] - ds.T[o]) / (d_o + d_n) rtol = 1e-12
        end
        carried[f] = sum(c.F_mass[:, f] .* component_enthalpies(laws, T_f))
        conducted[f] = c.F_energy[f] - carried[f]
    end
    # The carried term is what the old core's S_conv approximated; it is not
    # small here, so a double count or an omission would show.
    @test maximum(abs, carried) > 1e-3 * maximum(abs, conducted)

    for i = 1:n_cells(layout)
        faces = topo.cells[i].faces
        σ(f) = topo.faces[f].cells[1] == i ? -1.0 : 1.0
        cond = sum(σ(f) * conducted[f] * g.face_area[f] for f in faces)
        carr = sum(σ(f) * carried[f] * g.face_area[f] for f in faces)
        formation = sum(c.S_mass[j, i] * ds.h[j, i] for j = 1:3)
        expected = cond + carr + c.S_energy[i] + formation
        dm = [du[mass_index(layout, i, j)] for j = 1:3]
        rate =
            energy === :enthalpy ? du[energy_index(layout, i)] :
            ds.cp_sum[i] * du[energy_index(layout, i)] + sum(dm .* ds.h[:, i])
        scale = abs(cond) + abs(carr) + abs(c.S_energy[i]) + abs(formation)
        @test abs(rate - expected) <= 1e-11 * scale
    end
end

# Converged solves: the Newton tolerance of the fixed steps, with the core's
# own per-row absolute tolerances at the same relative tolerance.
const _LEDGER_RELTOL = 1e-10

# The solve augmented with the ledger integral E(t), so the integrator's own
# quadrature accounts for it: Σ H(t) − E(t) is a linear invariant of the
# inventory form, exact to round-off at every step, and converges at the
# integrator's order in the temperature form.
function _ledger_defect(cp, dt, t_end)
    def, ws = cp.def, cp.ws
    n = length(cp.u0)
    function rhs!(dv, v, p, t)
        u = view(v, 1:n)
        du = view(dv, 1:n)
        residual!(du, u, def, t, ws)
        ds, c = buffers(ws, def, eltype(v))
        dv[n+1] = _ledger_rate(def, c, ds)
        return nothing
    end
    function jac!(J, v, p, t)
        ForwardDiff.jacobian!(J, (dv, x) -> rhs!(dv, x, p, t), zeros(n + 1), v)
        return nothing
    end
    v0 = vcat(cp.u0, 0.0)
    f = ODEFunction(rhs!; jac = jac!)
    sol = OrdinaryDiffEq.solve(
        ODEProblem(f, v0, (0.0, t_end)),
        KenCarp4(autodiff = OrdinaryDiffEq.AutoFiniteDiff());
        adaptive = false,
        dt,
        reltol = _LEDGER_RELTOL,
        abstol = vcat(cp.abstol, 1e-6),
    )
    H0 = _enthalpy(def, ws, cp.u0)
    defect = maximum(abs(_enthalpy(def, ws, v[1:n]) - H0 - v[n+1]) for v in sol.u)
    scale = maximum(abs(v[n+1]) for v in sol.u)
    return sol.retcode, defect, scale
end

@testset "energy ledger over a solve: $(darcy ? "Darcy" : "Fick")" for darcy in
                                                                       (false, true)
    t_end = 10.0
    steps = darcy ? (0.05, 0.025, 0.0125) : (0.1, 0.05, 0.025)
    tols = (reltol = _LEDGER_RELTOL, abstol = _LEDGER_RELTOL)
    retcode, defect, scale = _ledger_defect(
        setup(_column(10; darcy); energy_state = :enthalpy, tols...),
        steps[1],
        t_end,
    )
    @test retcode == OrdinaryDiffEq.ReturnCode.Success
    @test scale > 1e5
    @test defect <= 1e-13 * scale

    defects = Float64[]
    for dt in steps
        retcode, defect, _ = _ledger_defect(
            setup(_column(10; darcy); energy_state = :temperature, tols...),
            dt,
            t_end,
        )
        @test retcode == OrdinaryDiffEq.ReturnCode.Success
        push!(defects, defect)
    end
    orders = log2.(defects[1:2] ./ defects[2:3])
    @info "temperature-form ledger defect ($(darcy ? "Darcy" : "Fick"))" defects orders
    @test all(>(0), defects)
    @test issorted(defects; rev = true)
    @test all(>=(2), orders)
end
