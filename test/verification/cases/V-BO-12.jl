# =============================================================================
# V-BO-12 — what deactivation does to a burn at an order other than one
# =============================================================================
#
# Two decks that a burn without deactivation can finish, so each has its own
# reference. In `a` a column at order zero leaves 0.2 % of its mass as a char of
# its own density: 0.2 % of its volume, so its cells go, and what that costs must
# stay inside what the mesh already costs, the difference of the same deck
# between 40 and 80 cells without deactivation. In `b` a column at order 0.5 leaves
# a char that keeps its volume, a layer the front burns beneath: its cells are not
# tested, no event happens, and the burn is the burn without deactivation.

const VBO12_L, VBO12_FLUX = 10e-3, 4e4
const VBO12_T_END = 400.0

function _vbo12_problem(N; n, y, ρ_c, k_c, φ_c = 0.3)
    material = Material(
        name = :VBO12,
        components = (
            SolidComponent(:a, ρ = 500.0, φ = 0.3, c = 1500.0, k = 0.2, ε = 0.9, λ = 1e-6),
            SolidComponent(:char, ρ = ρ_c, φ = φ_c, c = 1500.0, k = k_c, ε = 0.9, λ = 1e-6),
            GaseousComponent(:gas, M = 0.03, c = 1500.0, k = 0.02),
        ),
        reactions = (
            Reaction(
                :burn,
                :a => (:char, :gas);
                A = 1e10,
                E = 1.5e5,
                h = 3e5,
                n,
                yields = (y, 1 - y),
            ),
        ),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(VBO12_L, N),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(VBO12_FLUX) + RadiativeBC(ε = 0.9, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 500.0, z -> 0.0, z -> 0.0],
        tspan = (0.0, VBO12_T_END),
        output_times = collect(0.0:0.5:VBO12_T_END),
    )
end

const VBO12_SOLVE_KWARGS = (;
    abstol = 1e-9,
    reltol = 1e-7,
    use_ale = true,
    min_thickness = 1e-5,
    dtmin = 1e-14,
    deactivation = true,
)

# Over 5 s to the end: the integrated difference of the gas generation of `test`
# and `ref` over the mass the reference lost, and the largest difference of their
# masses over the initial mass.
function _vbo12_difference(ref, test)
    window = findall(t -> t >= 5.0, ref.t)
    t = ref.t[window]
    d = abs.(
        test.extras.gas_generation_rate[window] .- ref.extras.gas_generation_rate[window],
    )
    M = ref.extras.total_mass
    lost = M[1] - M[window[end]]
    integral = sum((d[j] + d[j+1]) / 2 * (t[j+1] - t[j]) for j = 1:(length(t)-1))
    return integral / lost,
    maximum(abs.(test.extras.total_mass[window] .- M[window])) / M[1]
end

function _verification_case_vbo12a()
    deck = (; n = 0.0, y = 0.002, ρ_c = 500.0, k_c = 0.2)
    off = (; VBO12_SOLVE_KWARGS..., deactivation = false)
    function qois(sol)
        coarse = suite_solve(_vbo12_problem(40; deck...); off...)
        fine = suite_solve(_vbo12_problem(80; deck...); off...)
        rate, mass = _vbo12_difference(coarse, sol)
        mesh_rate, mesh_mass = _vbo12_difference(coarse, fine)
        return [
            QoIResult(
                name = "reaches the end",
                value = completed(sol) ? 1.0 : 0.0,
                exact = 1.0,
                error = completed(sol) ? 0.0 : 1.0,
                tol = 0.5,
                comment = "the solve ended with outcome :$(sol.outcome)",
            ),
            QoIResult(
                name = "integrated gas generation, deactivation over the mesh",
                value = rate / mesh_rate,
                error = rate / mesh_rate,
                tol = 1.0,
                comment = "$(rate) against $(mesh_rate)",
            ),
            QoIResult(
                name = "mass, deactivation over the mesh",
                value = mass / mesh_mass,
                error = mass / mesh_mass,
                tol = 1.0,
                comment = "$(mass) against $(mesh_mass)",
            ),
        ]
    end
    return VerificationCase(
        id = "V-BO-12a",
        title = "Burnout at order zero against the mesh's own difference",
        tags = [:ale, :burnout, :depletion, :coupled],
        references = ["the same deck without deactivation at 40 and 80 cells"],
        build = N -> _vbo12_problem(N; deck...),
        solve_kwargs = VBO12_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = 40, budget_s = 240.0),
    )
end

function _verification_case_vbo12b()
    deck = (; n = 0.5, y = 0.2, ρ_c = 100.0, k_c = 0.1, φ_c = 0.8)
    function qois(sol)
        off = suite_solve(
            _vbo12_problem(40; deck...);
            VBO12_SOLVE_KWARGS...,
            deactivation = false,
        )
        rate, mass = _vbo12_difference(off, sol)
        return [
            QoIResult(
                name = "events",
                value = Float64(length(sol.events)),
                exact = 0.0,
                error = Float64(length(sol.events)),
                tol = 0.5,
                comment = "the char keeps the cells' volume",
            ),
            QoIResult(
                name = "integrated gas generation against no deactivation",
                value = rate,
                error = rate,
                tol = 1e-3,
                comment = "over the mass the reference lost",
            ),
            QoIResult(
                name = "mass against no deactivation",
                value = mass,
                error = mass,
                tol = 1e-4,
                comment = "over the initial mass",
            ),
        ]
    end
    return VerificationCase(
        id = "V-BO-12b",
        title = "A char that keeps its volume keeps its cells",
        tags = [:ale, :burnout, :depletion, :coupled],
        references = ["the same deck without deactivation"],
        build = N -> _vbo12_problem(N; deck...),
        solve_kwargs = VBO12_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = 40, budget_s = 120.0),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo12a)
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo12b)
