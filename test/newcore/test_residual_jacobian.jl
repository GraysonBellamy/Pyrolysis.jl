using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# A reacting, conducting slab: a solid that decomposes to a char and a
# retained gas, temperature-dependent properties, a Dirichlet back face and a
# heated front face.
function _slab(n; top = HeatFluxBC(3e4), bottom = DirichletThermalBC(400.0))
    virgin = SolidComponent(:virgin, ρ = 1100.0, c = (900.0, 1.5), k = (0.2, 2e-4))
    char = SolidComponent(
        :char,
        ρ = 300.0,
        c = PolynomialProperty((700.0, 0.8, -1e-4)),
        k = 0.1,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :slab,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    mesh = create_stretched_mesh(0.01, n, Val(3), 1.15)
    bcs = BoundaryConditionSet(thermal = (top = top, bottom = bottom))
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bcs,
        T_initial = z -> 400.0 + 4e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.5],
        tspan = (0.0, 10.0),
    )
end

# States on and off the physical manifold.
function _states(cp)
    u0 = cp.u0
    layout = cp.def.layout
    hot = copy(u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 3
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    zero_state = zeros(length(u0))
    negative = -abs.(u0)
    return (; u0, hot, perturbed, zero_state, negative)
end

@testset "residual and colored Jacobian: $(energy) form, $(n) cells" for energy in (
        :temperature,
        :enthalpy,
    ),
    n in (6, 8)

    cp = setup(_slab(n); energy_state = energy)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    NC = n_components(layout)
    pattern = cj.pattern
    @test same_structure(jac_prototype(cj), pattern)
    @test all(pattern[r, r] for r in axes(pattern, 1))
    # Fixed mesh: a cell row reaches the three cells i−1…i+1.
    @test cj.n_colors == 3 * (NC + 1)

    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test all(isfinite, du)
    residual!(du, cp.u0, def, 0.0, ws)
    @test (@allocated residual!(du, cp.u0, def, 0.0, ws)) == 0

    J = jac_prototype(cj)
    for (name, u) in pairs(_states(cp))
        @testset "$name" begin
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, pattern))
            jacobian!(J, cj, u, 1.0)
            @test same_structure(J, pattern)
            cmp = compare_entries(Matrix(J), Jd; rtol = 1e-12)
            @test cmp.passed
        end
    end
end

@testset "row assembly conserves mass and scatters fluxes antisymmetrically" begin
    cp = setup(
        _slab(7; top = AdiabaticBC(), bottom = AdiabaticBC());
        energy_state = :enthalpy,
    )
    def, ws = cp.def, cp.ws
    layout = def.layout
    du = zeros(length(cp.u0))
    u = _states(cp).hot
    residual!(du, u, def, 0.0, ws)
    # Closed slab: the component mass rates sum to the reaction sources, whose
    # total over components vanishes (unit stoichiometry).
    total = sum(du[mass_index(layout, i, j)] for i = 1:n_cells(layout), j = 1:3)
    scale = maximum(abs(du[mass_index(layout, i, 1)]) for i = 1:n_cells(layout))
    @test abs(total) <= 1e-12 * scale
    # Inventory form, adiabatic faces, every reaction heat and formation term
    # included: Σ dHᵢ/dt equals the summed reaction heat plus formation terms,
    # with the conduction fluxes cancelling pairwise.
    ds, c = buffers(ws, def, Float64)
    expected = sum(
        c.S_energy[i] + sum(c.S_mass[j, i] * ds.h[j, i] for j = 1:3) for
        i = 1:n_cells(layout)
    )
    dH = sum(du[energy_index(layout, i)] for i = 1:n_cells(layout))
    @test dH ≈ expected rtol = 1e-12
end
