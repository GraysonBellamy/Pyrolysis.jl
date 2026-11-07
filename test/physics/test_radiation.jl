"""
Tests for the exposed face under radiative transfer. The face is a gray
interface whose transmissivity is the first live cell's effective emissivity,
at its temperature and composition, for the radiation entering and for the
radiation leaving (detailed balance). With no component absorptivity a beam
enters at that emissivity too. The cells absorb with the volumetric absorption
coefficient α_eff = Σ_k ξ_k α_k of their own composition.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: RADIATIVE_TRANSFER, STEFAN_BOLTZMANN
using Pyrolysis.Properties: mixture_emissivity, effective_absorption

# Two solids with different ε so ε_eff is composition-dependent, and non-zero
# absorption coefficients. No component sets an absorptivity.
function _radiation_material()
    virgin = SolidComponent(
        :virgin,
        φ = 0.3,
        ρ = 600.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.95,
        α = 2.0,
        λ = 1e-5,
    )
    char = SolidComponent(
        :char,
        φ = 0.3,
        ρ = 150.0,
        c = 1700.0,
        k = 0.08,
        ε = 0.85,
        α = 50.0,
        λ = 2e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03)
    rxn = Reaction(:pyro, 1 => (2, 3), yields = (0.3, 0.7), A = 1e10, E = 150e3, h = -300e3)
    return Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
end

# Graded hot state with a graded mixed composition (ε_eff strictly between
# the component values in every cell), with consistent radiation rows.
_radiation_composition(i) = (400.0 - 10.0 * i, 80.0 + 5.0 * i, 4.0)

function _radiation_state(; n = 8)
    mat = _radiation_material()
    problem = PyrolysisProblem(
        mesh = uniform_mesh(0.01, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(50e3) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
        ),
        tspan = (0.0, 100.0),
    )
    cp = Pyrolysis.Integration.setup(problem; radiation_model = RADIATIVE_TRANSFER)
    def, ws = cp.def, cp.ws
    layout = def.layout
    V = Pyrolysis.Problem.cell_volumes(problem.mesh)
    u = copy(cp.u0)
    T = [600.0 + 300.0 * (i - 1) / (n - 1) for i = 1:n]
    for i = 1:n
        u[Pyrolysis.States.energy_index(layout, i)] = T[i]
        for (j, ξ) in enumerate(_radiation_composition(i))
            u[Pyrolysis.States.mass_index(layout, i, j)] = ξ * V[i]
        end
    end
    Pyrolysis.Integration._consistent_radiation!(u, def, ws, 0.0)
    du = zeros(length(u))
    Pyrolysis.Assembly.residual!(du, u, def, 0.0, ws)
    ds, c = Pyrolysis.Assembly.buffers(ws, def, Float64)
    K = Pyrolysis.Kernels
    kernel = K.find_kernel(K.RadiativeTransferKernel, def.kernels)
    return (; mat, cp, def, layout, u, ds, c, kernel, T, n)
end

@testset "the absorption coefficient is Σ ξ_k α_k of each cell's composition" begin
    s = _radiation_state()
    one_component = Material(
        name = :one,
        components = (
            SolidComponent(:s, λ = 0.0, ρ = 1000.0, c = 1000.0, k = 0.2, α = 5.0),
        ),
        reactions = (),
    )
    @test effective_absorption((1000.0,), one_component) ≈ 1000.0 * 5.0 rtol = 1e-14
    # Two solids and a gas without an absorption coefficient.
    @test effective_absorption((350.0, 120.0, 4.0), s.mat) ≈ 350.0 * 2.0 + 120.0 * 50.0 rtol =
        1e-14
    # The cells' coefficients are those sums at their own concentrations.
    for i = 1:s.n
        ξ = _radiation_composition(i)
        @test s.ds.α[i] ≈ ξ[1] * 2.0 + ξ[2] * 50.0 rtol = 1e-13
    end
    @test allunique(s.ds.α)
end

# Under radiative transfer the exposed face is a gray interface whose
# transmissivity is the first live cell's emissivity, at its temperature and
# composition, for the radiation entering and for the radiation leaving: the
# exposed row reads q⁻ = ε Q + (1 − ε) q⁺, so ε of the diffuse Q enters and
# ε of the q⁺ arriving from inside leaves.
@testset "Radiative transfer: the exposed face transmits at the first live cell's emissivity, both ways" begin
    s = _radiation_state()
    mat, cp, def, layout, u, ds, c, kernel, T, n =
        s.mat, s.cp, s.def, s.layout, s.u, s.ds, s.c, s.kernel, s.T, s.n
    K = Pyrolysis.Kernels
    Q = 50e3 + STEFAN_BOLTZMANN * 300.0^4
    for i = 1:n
        ξ = _radiation_composition(i)
        @test ds.ε[i] ≈ mixture_emissivity(ξ, T[i], mat) rtol = 1e-13
        # With no component absorptivity a beam enters at the emissivity.
        @test ds.a[i] === ds.ε[i]
    end
    @test 0.85 < ds.ε[n] < 0.95
    e = n + 1
    ε = ds.ε[n]
    @test ds.q[1, e] ≈ ε * Q + (1 - ε) * ds.q[2, e] rtol = 1e-13
    # The colored Jacobian: the exposed row reflects 1 − ε of q⁺.
    J = Pyrolysis.ColoredJacobian.jac_prototype(cp.jacobian)
    Pyrolysis.ColoredJacobian.jacobian!(J, cp.jacobian, u, 0.0)
    r, col = Pyrolysis.States.radiation_index(layout, e, 1),
    Pyrolysis.States.radiation_index(layout, e, 2)
    @test J[r, col] ≈ -(1 - ε) rtol = 1e-13
    # With the exposed cell dead the next cell's emissivity takes over, at the
    # face below it.
    ds.active[n] = false
    K.reset!(c)
    K.contribute!(c, kernel, ds, def.setup.topology, mat, 0.0)
    ε′ = ds.ε[n-1]
    @test ε′ != ε
    @test c.r_radiation[1, n] ≈ ds.q[1, n] - (ε′ * Q + (1 - ε′) * ds.q[2, n]) rtol = 1e-13
    @test c.r_radiation[1, e] == 0 && c.r_radiation[2, e] == 0
end
