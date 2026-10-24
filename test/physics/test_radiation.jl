"""
Tests for Beer–Lambert surface absorptivity. Absorption and emission use
the same surface coefficient, so switching SURFACE_ABSORPTION →
BEER_LAMBERT does not change the absorbed power by 1/ε_eff (≈ 11 % at
ε = 0.9).

Convention:
- Kirchhoff path (no explicit BC absorptivity, no component absorptivity):
  the flux entering the medium is ε_eff·I, with ε_eff the exposed cell's
  effective emissivity — absorption and emission use the same coefficient
  (detailed balance).
- Explicit-absorptivity path: the condition's α(t) replaces it.
- The cells attenuate the beam with the volumetric absorption coefficient
  α_eff = Σ_k ξ_k α_k of their own composition.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Properties: mixture_emissivity, effective_absorption

# Two solids with different ε so ε_eff is composition-dependent; non-zero
# absorption coefficients so Beer–Lambert is non-trivial. No component sets
# an absorptivity: the surface absorbs at its emissivity.
function _radiation_material()
    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.95,
        α = 2.0,
        λ = 1e-5,
    )
    char =
        SolidComponent(:char, ρ = 150.0, c = 1700.0, k = 0.08, ε = 0.85, α = 50.0, λ = 2e-5)
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyro, 1 => (2, 3), yields = (0.3, 0.7), A = 1e10, E = 150e3, h = -300e3)
    return Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
end

# Graded hot state with a graded mixed composition (ε_eff strictly between
# the component values in every cell), with consistent intensities.
_radiation_composition(i) = (400.0 - 10.0 * i, 80.0 + 5.0 * i, 4.0)

function _radiation_state(; n = 8)
    mat = _radiation_material()
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50e3) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
        ),
        tspan = (0.0, 100.0),
    )
    cp = Pyrolysis.Integration.setup(problem; radiation_model = BEER_LAMBERT)
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
    Pyrolysis.Integration._consistent_intensities!(u, def, ws, 0.0)
    du = zeros(length(u))
    Pyrolysis.Assembly.residual!(du, u, def, 0.0, ws)
    ds, _ = Pyrolysis.Assembly.buffers(ws, def, Float64)
    kernel = Pyrolysis.Kernels.find_kernel(Pyrolysis.Kernels.BeerLambertKernel, def.kernels)
    return (; mat, def, ds, kernel, T, n)
end

@testset "Beer–Lambert surface absorptivity" begin
    s = _radiation_state()
    ds, kernel, mat, n = s.ds, s.kernel, s.mat, s.n
    i_top = kernel.order[1]
    @test i_top == n

    @testset "Kirchhoff path: the exposed cell absorbs at the emissivity it emits with" begin
        # No component absorptivity: the surface absorptivity is the
        # emissivity itself, between the components' (mixed composition).
        @test ds.a[i_top] === ds.ε[i_top]
        @test 0.85 < ds.ε[i_top] < 0.95
        @test ds.ε[i_top] ≈ mixture_emissivity(_radiation_composition(n), s.T[n], mat) rtol =
            1e-13
        # The beam entering the medium carries that absorptivity.
        @test Pyrolysis.Kernels.incident_intensity(kernel, ds, 0.0) ≈ ds.ε[i_top] * 50e3 rtol =
            1e-13
        @test ds.I[kernel.exposed] ≈ ds.ε[i_top] * 50e3 rtol = 1e-13
    end

    @testset "the absorption coefficient is Σ ξ_k α_k of each cell's composition" begin
        one_component = Material(
            name = :one,
            components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2, α = 5.0),),
            reactions = (),
        )
        @test effective_absorption((1000.0,), one_component) ≈ 1000.0 * 5.0 rtol = 1e-14
        # Two solids and a gas without an absorption coefficient.
        @test effective_absorption((350.0, 120.0, 4.0), mat) ≈ 350.0 * 2.0 + 120.0 * 50.0 rtol =
            1e-14
        # The cells' coefficients are those sums at their own concentrations.
        for i = 1:n
            ξ = _radiation_composition(i)
            @test ds.α[i] ≈ ξ[1] * 2.0 + ξ[2] * 50.0 rtol = 1e-13
        end
        @test allunique(ds.α)
    end
end
