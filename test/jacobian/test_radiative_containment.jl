using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: RADIATIVE_TRANSFER
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Closures: EquidistributionCorrection
using Pyrolysis.Integration: setup, DomainCheck

# The radiative transfer rows in the colored Jacobian: the pattern contains
# the dense Jacobian at every state, the colored Jacobian equals it per entry,
# and one residual, one colored Jacobian and one domain test allocate nothing.

function _pmma()
    virgin = SolidComponent(
        :virgin,
        φ = 0.05,
        λ = 1e-5,
        ρ = 1190.0,
        c = 1420.0,
        k = 0.21,
        ε = 0.86,
        α = 2.0,
    )
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood(; α = 1000.0, κ = 1e-14, law = nothing)
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        ρ_intrinsic = 1500.0,
        c = (1500.0, 0.5),
        k = 0.15,
        ε = (0.9, -1e-4),
        α = α,
        κ = κ,
        λ = 1e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    return Material(
        name = :wood,
        components = (virgin, gas),
        reactions = (rxn,),
        lateral_shrinkage_law = law,
    )
end

const _top = RadiativeFluxBC(t -> 5e4 * (1 + 0.2 * sin(t))) + RadiativeBC(T_inf = 300.0)
const _open = RadiativeBC(T_inf = 350.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)
# A collimated heater at its own absorptivity, and a cone beside a diffuse flux.
const _collimated =
    RadiativeFluxBC(
        t -> 5e4 * (1 + 0.2 * sin(t));
        absorptivity = t -> 0.9 - 0.01 * sin(t),
        distribution = Collimated(),
    ) + RadiativeBC(T_inf = 300.0)
const _cone =
    RadiativeFluxBC(4e4; distribution = Cone(0.35, 1.27)) +
    RadiativeFluxBC(1e4) +
    RadiativeBC(T_inf = 300.0)

_column(mat, n; top = _top, bottom = AdiabaticBC(), stretch = nothing) = PyrolysisProblem(
    mesh = stretch === nothing ? uniform_mesh(0.01, n) : stretched_mesh(0.01, n, stretch),
    material = mat,
    bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
    tspan = (0.0, 100.0),
)

const _configurations = (
    eulerian = (problem = n -> _column(_pmma(), n), options = (;), wood = false),
    eulerian_open = (
        problem = n -> _column(_pmma(), n; bottom = _open),
        options = (;),
        wood = false,
    ),
    ale_darcy = (
        problem = n -> _column(_wood(), n),
        options = (use_ale = true,),
        wood = true,
    ),
    ale_darcy_moderate = (
        problem = n -> _column(_wood(α = 0.5), n; bottom = _open),
        options = (use_ale = true,),
        wood = true,
    ),
    ale_wall_moderate = (
        problem = n -> _column(_wood(α = 0.5, κ = 0.0), n),
        options = (use_ale = true,),
        wood = true,
    ),
    cross_section = (
        problem = n -> _column(_wood(α = 5.0, law = χ -> 1 - 0.3χ), n),
        options = (use_ale = true,),
        wood = true,
    ),
    equidistribution = (
        problem = n -> _column(_pmma(), n; stretch = 1.05),
        options = (
            use_ale = true,
            correction = EquidistributionCorrection(activation = (2.0, 0.5)),
        ),
        wood = false,
    ),
    eulerian_collimated = (
        problem = n -> _column(_pmma(), n; top = _collimated, bottom = _open),
        options = (;),
        wood = false,
        rows = 3,
    ),
    ale_darcy_cone = (
        problem = n -> _column(_wood(α = 0.5), n; top = _cone, bottom = _open),
        options = (use_ale = true,),
        wood = true,
        rows = 5,
    ),
    equidistribution_collimated = (
        problem = n -> _column(_pmma(), n; top = _collimated, stretch = 1.05),
        options = (
            use_ale = true,
            correction = EquidistributionCorrection(activation = (2.0, 0.5)),
        ),
        wood = false,
        rows = 3,
    ),
)

# The state of `cp` with cell temperatures `T(i)` and concentrations `ξ(i)`
# [kg/m³]; every other unknown as at the start.
function _cell_state(cp, T, ξ)
    u = copy(cp.u0)
    layout = cp.def.layout
    V = cp.def.setup.reference.volume
    for i = 1:n_cells(layout)
        m = Tuple(ξ(i) .* V[i])
        for j in eachindex(m)
            u[mass_index(layout, i, j)] = m[j]
        end
        u[energy_index(layout, i)] =
            energy_state_value(energy_state(layout), cp.def.setup, T(i), m, i)
    end
    return u
end

function _states(cp, wood::Bool)
    n = n_cells(cp.def.layout)
    u0 = cp.u0
    solid = wood ? 500.0 : 1190.0
    uniform(T, ξ) = _cell_state(cp, i -> T, i -> ξ)
    return (;
        u0,
        mixed = _cell_state(
            cp,
            i -> 450.0 + 200.0 * i / n,
            i -> (solid * (0.65 + 0.3 * i / n), 0.02 + 0.08 * i / n),
        ),
        hot = uniform(1000.0, wood ? (400.0, 0.0) : (1100.0, 0.0)),
        gas_heavy = wood ? uniform(1000.0, (400.0, 5.0)) : uniform(1200.0, (400.0, 8.0)),
        partial = uniform(800.0, wood ? (200.0, 2.0) : (700.0, 2.0)),
        perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0))),
        zero_state = zeros(length(u0)),
        negative = -abs.(u0),
    )
end

@testset "radiative transfer $(name): $(energy) form, $(n) cells" for name in
                                                                      keys(_configurations),
    energy in (:temperature, :enthalpy),
    n in (6, 8)

    config = _configurations[name]
    cp = setup(
        config.problem(n);
        energy_state = energy,
        radiation_model = RADIATIVE_TRANSFER,
        config.options...,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    @test radiation_rows(def.layout) == get(config, :rows, 2)
    du = zeros(length(cp.u0))
    for (state, u) in pairs(_states(cp, config.wood))
        @testset "$state" begin
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, cj.pattern))
            J = jacobian!(jac_prototype(cj), cj, u, 1.0)
            @test same_structure(J, cj.pattern)
            @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
        end
    end
end

# Allocations of one residual, one colored Jacobian and one domain test after
# warm-up, measured inside a function, at a hundred cells.
function _allocations(cp)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    residual!(du, u, cp.def, 1.0, cp.ws)
    jacobian!(J, cp.jacobian, u, 1.0)
    jacobian!(J, cp.jacobian, u, 1.0)
    check = DomainCheck(cp)
    check(u, nothing, 1.0)
    r = @allocated residual!(du, u, cp.def, 1.0, cp.ws)
    j = @allocated jacobian!(J, cp.jacobian, u, 1.0)
    d = @allocated check(u, nothing, 1.0)
    return r, j, d
end

@testset "radiative transfer allocations: $(name), $(energy) form" for name in (
        :eulerian_open,
        :ale_darcy,
        :cross_section,
        :equidistribution,
        :eulerian_collimated,
        :ale_darcy_cone,
        :equidistribution_collimated,
    ),
    energy in (:temperature, :enthalpy)

    config = _configurations[name]
    cp = setup(
        config.problem(100);
        energy_state = energy,
        radiation_model = RADIATIVE_TRANSFER,
        config.options...,
    )
    r, j, d = _allocations(cp)
    @test r == 0
    @test j == 0
    @test d == 0
end
