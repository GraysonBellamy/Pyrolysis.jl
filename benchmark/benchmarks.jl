"""
PkgBenchmark.jl entry point for Pyrolysis.jl.

Defines a top-level `SUITE::BenchmarkGroup` covering three canonical
problems:

- PMMA cone (fixed mesh, surface absorption)
- White-pine-style charring (moving mesh, Beer–Lambert + lumped and Darcy)
- Generic charring polymer (moving mesh, lumped, no radiation)

For each problem it benchmarks the three hot operations on the problem
`setup` builds:

- One residual evaluation
- One colored Jacobian
- An end-to-end solve over a short time span

Run from the package root with:

```julia
using PkgBenchmark
results = benchmarkpkg("Pyrolysis")
export_markdown("benchmark/report.md", results)
```

For run-over-run regression detection, compare against a baseline:

```julia
judge = judge("Pyrolysis", "main")
```
"""

using BenchmarkTools
using Pyrolysis
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian: jac_prototype, jacobian!

const SUITE = BenchmarkGroup()

const _N_CELLS = 100

function _pmma()
    virgin = SolidComponent(
        :virgin,
        ρ = 1190.0,
        φ = 0.01,
        λ = 1e-6,
        c = 1420.0,
        k = 0.21,
        ε = 0.86,
    )
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3)
    return Material(
        name = :PMMA,
        gas_transport = LumpedTransfer(
            endpoint = LumpedEndpoint(λ = 1e-6, κ = 0.0, max_skeletal_occupancy = 0.001),
        ),
        components = (virgin, gas),
        reactions = (rxn,),
    )
end

function _wood()
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = (1500.0, 0.5),
        k = 0.15,
        ε = 0.9,
        α = 500.0,
        κ = 1e-14,
        λ = 1e-5,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = 300e3)
    return Material(
        name = :wood,
        gas_transport = LumpedTransfer(
            endpoint = LumpedEndpoint(λ = 1e-5, κ = 1e-14, max_skeletal_occupancy = 0.001),
        ),
        components = (virgin, gas),
        reactions = (rxn,),
    )
end

# A 1 cm column under the cone boundary conditions; a permeable material
# vents at the ambient pressure.
function _column(material; darcy::Bool = false, t_final = 100.0)
    top_mass =
        darcy ?
        AmbientGasBC(
            GasReservoir(P = 101325.0, T = 300.0, mole_fractions = (pyrgas = 1.0,)),
        ) : AmbientGasBC(GasReservoir(P = 0.0, T = 300.0))
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, _N_CELLS),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50e3) +
                      RadiativeBC(T_inf = 300.0) +
                      ConvectiveBC(h = 10.0, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = top_mass, bottom = ImpermeableBC()),
        ),
        tspan = (0.0, t_final),
        dt_initial = 1e-3,
        output_times = [0.0, t_final / 2, t_final],
    )
end

const _CONFIGURATIONS = (
    pmma_fixed_surface = (
        material = _pmma,
        darcy = false,
        options = (radiation_model = SURFACE_ABSORPTION,),
    ),
    wood_moving_beer_lambert_darcy = (
        material = _wood,
        darcy = true,
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
    ),
    polymer_moving_lumped = (material = _pmma, darcy = false, options = (use_ale = true,)),
)

for (name, config) in pairs(_CONFIGURATIONS)
    group = SUITE[string(name)] = BenchmarkGroup()
    cp = setup(_column(config.material(); config.darcy); config.options...)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    group["residual"] = @benchmarkable residual!($du, $u, $(cp.def), 0.0, $(cp.ws))
    group["jacobian"] = @benchmarkable jacobian!($J, $(cp.jacobian), $u, 0.0)
    problem = _column(config.material(); config.darcy, t_final = 5.0)
    group["solve"] = @benchmarkable solve($problem; $(config.options)...)
end
