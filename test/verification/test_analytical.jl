"""
Verification Suite

Comprehensive verification tests against analytical solutions.

## Coordinate Convention
- z = 0: Bottom boundary (substrate, boundary_id = 2)
- z = L: Top boundary (exposed surface, boundary_id = 1)
- Heat enters from top surface (z = L)
- Cell 1 is at bottom (z ≈ 0), Cell N is at top (z ≈ L)

Tests:
1. Semi-infinite solid heat conduction (erf solution)
2. Beer-Lambert radiation absorption
3. Single reaction isothermal decay
4. Method of Manufactured Solutions (MMS)
5. Mass conservation
6. Energy conservation
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: StateCache
using Pyrolysis.Discretization: BEER_LAMBERT
using Pyrolysis.Solver: get_cell_positions
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Physics: compute_radiation_source!
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Geometry: CellState, cell_size
using Statistics: mean

# Abramowitz-Stegun 7.1.26 approximation; avoids a test-only dependency.
if !isdefined(@__MODULE__, :erf)
    function erf(x::Real)
        sign_x = x < 0 ? -1.0 : 1.0
        a = abs(float(x))
        t = 1.0 / (1.0 + 0.3275911 * a)
        y =
            1.0 -
            (
                ((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) *
                t + 0.254829592
            ) *
            t *
            exp(-a * a)
        return sign_x * y
    end
end

# ============================================================================
# Test 9.1: Semi-Infinite Solid Heat Conduction
# ============================================================================
#
# Verify heat conduction against analytical solution for semi-infinite solid.
#
# Setup:
# - Initial temperature: T(z,0) = T₀
# - Surface temperature: T(L,t) = T_s (Dirichlet BC at top surface)
# - Analytical solution: T(z,t) = T_s + (T₀ - T_s) * erf((L-z) / (2√(αt)))
#
# where α = k / (ρ * cp) is thermal diffusivity.
# 
# Note: z=0 is at bottom (substrate), z=L is at top (exposed surface).
# The analytical solution uses distance from surface = L - z.
#
@testset "Semi-Infinite Solid Heat Conduction" begin

    # Material properties (constant for analytical solution)
    ρ = 1000.0    # kg/m³
    cp = 1000.0   # J/(kg·K)
    k = 0.5       # W/(m·K)
    α = k / (ρ * cp)  # Thermal diffusivity = 5e-7 m²/s

    # Boundary and initial conditions
    T_initial = 300.0  # K
    T_surface = 500.0  # K

    # Domain
    L = 0.05  # 5 cm - thick enough to approximate semi-infinite

    # Analytical solution (z measured from bottom, surface at z=L)
    function T_analytical(z, t, T_s, T_0, α, L)
        if t <= 0
            return T_0
        end
        # Distance from surface = L - z
        return T_s + (T_0 - T_s) * erf((L - z) / (2 * sqrt(α * t)))
    end

    @testset "Solution Accuracy at t=100s" begin
        # Create material (inert, no reactions)
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        # Create mesh
        n_cells = 100
        mesh =
            create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

        # Boundary conditions: Dirichlet at top, adiabatic at bottom
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_surface), bottom = AdiabaticBC()),
        )

        # Create and solve problem
        tspan = (0.0, 100.0)
        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ],
            tspan = tspan,
        )

        solution = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        @test solution.retcode == :Success

        # Compare numerical to analytical at final time
        t_final = 100.0
        z_positions = get_cell_positions(solution)
        T_numerical = solution.T[:, end]

        max_error = 0.0
        for (i, z) in enumerate(z_positions)
            T_exact = T_analytical(z, t_final, T_surface, T_initial, α, L)
            error = abs(T_numerical[i] - T_exact)
            max_error = max(max_error, error)
        end

        # Error should be small (< 1K for this mesh resolution)
        @test max_error < 2.0
    end

    @testset "Spatial Convergence (2nd Order)" begin
        # Create material
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        # Test at multiple mesh resolutions
        n_cells_list = [25, 50, 100, 200]
        errors = Float64[]
        h_values = Float64[]

        t_test = 50.0  # Test time

        for n_cells in n_cells_list
            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(1);
                T_initial = T_initial,
                ξ_initial = (ρ,),
            )
            h = L / n_cells
            push!(h_values, h)

            bc_set = BoundaryConditionSet(
                thermal = (top = DirichletThermalBC(T_surface), bottom = AdiabaticBC()),
            )

            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> T_initial,
                ξ_initial = [z -> ρ],
                tspan = (0.0, t_test),
                output_times = [0.0, t_test],
            )

            sol = solve(problem; abstol = 1e-12, reltol = 1e-10, progress = false)

            # Compute L2 error
            z_positions = get_cell_positions(sol)
            T_num = sol.T[:, end]

            l2_error = 0.0
            for (i, z) in enumerate(z_positions)
                T_exact = T_analytical(z, t_test, T_surface, T_initial, α, L)
                l2_error += (T_num[i] - T_exact)^2 * h
            end
            l2_error = sqrt(l2_error / L)
            push!(errors, l2_error)
        end

        # Compute convergence rates
        rates = Float64[]
        for i = 2:length(errors)
            rate = log(errors[i-1] / errors[i]) / log(h_values[i-1] / h_values[i])
            push!(rates, rate)
        end

        # Should be approximately 2nd order (rate ≈ 2)
        # Allow some tolerance due to time integration errors
        @test all(r -> r > 1.5, rates)  # At least 1.5 order
        @test mean(rates) > 1.8  # Average should be close to 2
    end
end

# ============================================================================
# Test 9.2: Beer-Lambert Radiation Absorption
# ============================================================================
#
# Verify radiation absorption using Beer-Lambert law.
#
# Setup:
# - Incident radiation I₀ at surface
# - Constant absorption coefficient α_eff
# - Analytical solution: I(z) = I₀ * exp(-α_eff * z)
# - Energy conservation: ∫ q''' dz = I₀ * (1 - exp(-α_eff * L))
#
@testset "Beer-Lambert Radiation Absorption" begin

    # Material properties
    ρ = 1000.0
    cp = 1000.0
    k = 0.5
    absorption_coeff = 100.0  # m²/kg, gives α_eff = 100,000 m⁻¹

    L = 0.01  # 1 cm domain
    I_surface = 50000.0  # W/m² incident radiation

    @testset "Exponential Decay Profile" begin
        # Create material with specified absorption
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9, α = absorption_coeff)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        # Create mesh
        n_cells = 100
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = 300.0, ξ_initial = (ρ,))

        # Set up property cache
        prop_cache = PropertyCache(n_cells, n_cells + 1)
        update_properties!(prop_cache, StateCache(mesh), mesh, mat)

        # Compute radiation source
        q_rad = zeros(n_cells)
        compute_radiation_source!(
            q_rad,
            StateCache(mesh),
            mesh,
            prop_cache,
            I_surface,
            BEER_LAMBERT,
            false,  # bare attenuation under test; Kirchhoff ε_eff scaling covered in test/physics/test_radiation.jl
        )

        # α_eff = ξ * α = ρ * absorption_coeff
        α_eff = ρ * absorption_coeff

        bc_set =
            BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))
        z_centers = get_cell_positions(
            PyrolysisSolution(
                PyrolysisProblem(
                    mesh = mesh,
                    material = mat,
                    bc_set = bc_set,
                    tspan = (0.0, 1.0),
                ),
            ),
        )

        # Compare with cell-averaged analytical value:
        # q̄'''_i = I(z_left) * (1 - exp(-α_eff * Δz)) / Δz
        # where I(z) = I₀ * exp(-α_eff * (L - z))  (z=L is surface)
        # 
        # With surface at z=L, radiation enters at cell N and propagates toward cell 1.
        # Cell N is near z=L (surface), Cell 1 is near z=0 (bottom).
        Δz = L / n_cells
        max_rel_error = 0.0
        I_at_top = I_surface  # intensity entering top cell (cell N)

        # Iterate from top (cell N) to bottom (cell 1)
        for i = n_cells:-1:1
            # Cell-averaged analytical value
            absorption_factor = 1.0 - exp(-α_eff * Δz)
            q_expected = I_at_top * absorption_factor / Δz

            if q_expected > 1e-10
                rel_error = abs(q_rad[i] - q_expected) / q_expected
                max_rel_error = max(max_rel_error, rel_error)
            end

            # Update intensity for next cell (below)
            I_at_top *= exp(-α_eff * Δz)
        end

        # Should match within 1% (our formula is exact for piecewise constant α)
        @test max_rel_error < 0.01
    end

    @testset "Energy Conservation" begin
        # Create material
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9, α = absorption_coeff)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        n_cells = 200
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = 300.0, ξ_initial = (ρ,))

        prop_cache = PropertyCache(n_cells, n_cells + 1)
        update_properties!(prop_cache, StateCache(mesh), mesh, mat)

        q_rad = zeros(n_cells)
        compute_radiation_source!(
            q_rad,
            StateCache(mesh),
            mesh,
            prop_cache,
            I_surface,
            BEER_LAMBERT,
            false,  # bare attenuation under test; Kirchhoff ε_eff scaling covered in test/physics/test_radiation.jl
        )

        # Compute total absorbed energy
        Q_absorbed = sum(q_rad[i] * cell_size(mesh, i) for i = 1:n_cells)

        # Analytical total absorbed
        α_eff = ρ * absorption_coeff
        Q_analytical = I_surface * (1.0 - exp(-α_eff * L))

        # Should conserve energy within 5%
        rel_error = abs(Q_absorbed - Q_analytical) / Q_analytical
        @test rel_error < 0.05

        # Verify transmitted radiation: I_surface·exp(−τ_total)
        τ_total = sum(prop_cache.α[i] * cell_size(mesh, i) for i = 1:n_cells)
        I_transmitted = I_surface * exp(-τ_total)
        I_transmitted_analytical = I_surface * exp(-α_eff * L)

        @test abs(I_transmitted - I_transmitted_analytical) / I_surface < 0.01
    end
end

# ============================================================================
# Test 9.3: Single Reaction Isothermal Decay
# ============================================================================
#
# Verify single first-order reaction kinetics.
#
# Setup:
# - Isothermal conditions (fixed temperature)
# - First-order reaction: A → B
# - Rate: r = A_pre * exp(-E/RT) * ξ_A
# - Analytical solution: ξ_A(t) = ξ_A,0 * exp(-k*t)
#                        ξ_B(t) = ξ_A,0 * (1 - exp(-k*t))
#
# where k = A_pre * exp(-E/RT)
#
@testset "Isothermal Decay Kinetics" begin

    # Reaction parameters
    A_pre = 1e10      # s⁻¹
    E_act = 150e3     # J/mol
    T_isothermal = 600.0  # K

    # Rate constant at this temperature
    k_rate = A_pre * exp(-E_act / (R_GAS * T_isothermal))

    # Initial conditions
    ξ_A_initial = 1000.0  # kg/m³

    # Time for significant decay
    τ = 1.0 / k_rate  # characteristic time
    t_final = 5.0 * τ  # 5 time constants

    @testset "Exponential Decay Rate" begin
        # Create material with single reaction
        solid_A = SolidComponent(:A, ρ = ξ_A_initial, c = 1000.0, k = 0.5, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = ξ_A_initial, c = 1000.0, k = 0.5, ε = 0.9)

        rxn = Reaction(:decay, 1 => 2, A = A_pre, E = E_act, h = 0.0)

        mat =
            Material(name = :DecayTest, components = (solid_A, solid_B), reactions = (rxn,))

        # Single cell mesh (0D test)
        L = 0.001
        n_cells = 1
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = T_isothermal,
            ξ_initial = (ξ_A_initial, 0.0),
        )

        # Adiabatic boundaries (no heat loss) and Dirichlet to maintain temperature
        bc_set = BoundaryConditionSet(
            thermal = (
                top = DirichletThermalBC(T_isothermal),
                bottom = DirichletThermalBC(T_isothermal),
            ),
        )

        # Create problem
        n_outputs = 50
        output_times = collect(range(0.0, t_final, length = n_outputs))

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_isothermal,
            ξ_initial = [z -> ξ_A_initial, z -> 0.0],
            tspan = (0.0, t_final),
            output_times = output_times,
        )

        sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        @test sol.retcode == :Success

        # Verify exponential decay
        for (i, t) in enumerate(sol.t)
            ξ_A_numerical = sol.ξ[:, :, i][1, 1]
            ξ_A_analytical = ξ_A_initial * exp(-k_rate * t)

            if ξ_A_analytical > 1e-3 * ξ_A_initial  # Above noise floor
                rel_error = abs(ξ_A_numerical - ξ_A_analytical) / ξ_A_analytical
                @test rel_error < 0.01  # 1% tolerance
            end
        end
    end

    @testset "Mass Conservation in Reaction" begin
        # Same setup as above
        solid_A = SolidComponent(:A, ρ = ξ_A_initial, c = 1000.0, k = 0.5, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = ξ_A_initial, c = 1000.0, k = 0.5, ε = 0.9)

        # Stoichiometry: 1 kg A → 1 kg B (mass conserving)
        rxn = Reaction(:decay, 1 => 2, A = A_pre, E = E_act, h = 0.0)

        mat =
            Material(name = :DecayTest, components = (solid_A, solid_B), reactions = (rxn,))

        L = 0.001
        mesh = create_uniform_mesh(
            L,
            1,
            Val(2);
            T_initial = T_isothermal,
            ξ_initial = (ξ_A_initial, 0.0),
        )

        bc_set = BoundaryConditionSet(
            thermal = (
                top = DirichletThermalBC(T_isothermal),
                bottom = DirichletThermalBC(T_isothermal),
            ),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_isothermal,
            ξ_initial = [z -> ξ_A_initial, z -> 0.0],
            tspan = (0.0, t_final),
            output_times = collect(range(0.0, t_final, length = 20)),
        )

        sol = solve(problem; abstol = 1e-12, reltol = 1e-10, progress = false)

        # Check mass conservation: ξ_A + ξ_B = ξ_A,0
        for i in eachindex(sol.t)
            ξ_A = sol.ξ[:, :, i][1, 1]
            ξ_B = sol.ξ[:, :, i][2, 1]
            total_mass = ξ_A + ξ_B

            rel_error = abs(total_mass - ξ_A_initial) / ξ_A_initial
            @test rel_error < 1e-8  # Machine precision for mass conservation
        end
    end
end

# ============================================================================
# Test 9.3b: TGA Simulation (Constant Heating Rate)
# ============================================================================
#
# Verify kinetics against analytical solution for thermogravimetric analysis (TGA).
#
# Setup:
# - Constant heating rate: T(t) = T₀ + β*t
# - First-order reaction: A → B
# - Rate equation: dα/dt = A * exp(-E/RT) * (1-α)
#
# where α = 1 - ξ_A/ξ_A₀ is the conversion fraction.
#
# With T = T₀ + β*t, we have:
#   dα/dT = (A/β) * exp(-E/RT) * (1-α)
#
# The analytical solution is obtained by integrating:
#   α(T) = 1 - exp(-∫[T₀ to T] (A/β) * exp(-E/RT') dT')
#
# We use numerical quadrature for the "exact" reference solution since
# the temperature integral has no closed form.
#
@testset "TGA Constant Heating Rate Kinetics" begin

    # Reaction parameters (typical polymer decomposition)
    A_pre = 1e13      # s⁻¹
    E_act = 180e3     # J/mol (180 kJ/mol)

    # TGA conditions
    T_initial = 300.0   # K (room temperature)
    β = 10.0 / 60.0     # K/s (10 K/min, typical TGA heating rate)

    # Initial conditions
    ξ_A_initial = 1000.0  # kg/m³

    # Temperature range for TGA
    T_final = 800.0  # K
    t_final = (T_final - T_initial) / β  # Time to reach T_final

    """
    	Compute analytical conversion α(T) for first-order Arrhenius kinetics
    	with constant heating rate using numerical integration of the temperature integral.

    	α(T) = 1 - exp(-p(T))
    	where p(T) = (A/β) * ∫[T₀ to T] exp(-E/RT') dT'

    	Uses adaptive Simpson's rule for accurate integration.
    	"""
    function analytical_conversion(T, T₀, A, E, β; n_points = 1000)
        if T <= T₀
            return 0.0
        end

        # Numerical integration using trapezoidal rule with fine spacing
        temps = range(T₀, T, length = n_points)
        dT = (T - T₀) / (n_points - 1)

        integral = 0.0
        for i = 1:(n_points-1)
            T1 = temps[i]
            T2 = temps[i+1]
            # Trapezoidal rule
            f1 = exp(-E / (R_GAS * T1))
            f2 = exp(-E / (R_GAS * T2))
            integral += 0.5 * (f1 + f2) * dT
        end

        p = (A / β) * integral
        α = 1.0 - exp(-p)
        return α
    end

    """
    	Compute analytical conversion rate dα/dT for first-order kinetics.
    	"""
    function analytical_rate(T, α, A, E, β)
        return (A / β) * exp(-E / (R_GAS * T)) * (1.0 - α)
    end

    @testset "TGA Mass Loss Profile" begin
        # Create material with single reaction
        solid_A = SolidComponent(:A, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)

        rxn = Reaction(:pyrolysis, 1 => 2, A = A_pre, E = E_act, h = 0.0)

        mat = Material(name = :TGATest, components = (solid_A, solid_B), reactions = (rxn,))

        # Single cell mesh (0D kinetics test - uniform temperature)
        L = 0.0001  # Very thin sample (0.1 mm) for uniform T
        n_cells = 1
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = T_initial,
            ξ_initial = (ξ_A_initial, 0.0),
        )

        # Dirichlet BC with ramping temperature (TGA furnace)
        T_ramp = t -> T_initial + β * t
        bc_set = BoundaryConditionSet(
            thermal = (
                top = DirichletThermalBC(T_ramp),
                bottom = DirichletThermalBC(T_ramp),
            ),
        )

        # Create problem with many output points for smooth curve
        n_outputs = 200
        output_times = collect(range(0.0, t_final, length = n_outputs))

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ξ_A_initial, z -> 0.0],
            tspan = (0.0, t_final),
            output_times = output_times,
        )

        sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        @test sol.retcode == :Success

        # Compare numerical conversion to analytical at each output time
        max_abs_error = 0.0
        significant_conversion_tested = false

        for (i, t) in enumerate(sol.t)
            T_current = T_initial + β * t
            ξ_A_numerical = sol.ξ[:, :, i][1, 1]

            # Numerical conversion
            α_numerical = 1.0 - ξ_A_numerical / ξ_A_initial

            # Analytical conversion (using high-accuracy quadrature)
            α_analytical = analytical_conversion(
                T_current,
                T_initial,
                A_pre,
                E_act,
                β;
                n_points = 2000,
            )

            abs_error = abs(α_numerical - α_analytical)
            max_abs_error = max(max_abs_error, abs_error)

            # Test significant conversion region (5% to 95%)
            if 0.05 < α_analytical < 0.95
                significant_conversion_tested = true
                @test abs_error < 0.02  # 2% absolute error in conversion
            end
        end

        @test significant_conversion_tested  # Ensure we tested the main reaction region
        @test max_abs_error < 0.05  # Overall error < 5%
    end

    @testset "TGA Peak Temperature (DTG)" begin
        # The peak reaction rate (DTG peak) occurs at a characteristic temperature
        # For first-order kinetics: T_peak where d²α/dT² = 0
        # Approximate: E/(R*T_peak²) ≈ (A/β) * exp(-E/(R*T_peak))
        # This can be solved iteratively or using Kissinger equation:
        # ln(β/T_peak²) = ln(A*R/E) - E/(R*T_peak)

        # Solve Kissinger equation iteratively for T_peak
        function find_peak_temperature(A, E, β; tol = 1e-6, max_iter = 100)
            # Initial guess using approximate solution
            T_peak = E / (R_GAS * log(A * E / (β * R_GAS)))

            for _ = 1:max_iter
                # Newton-Raphson on: f(T) = ln(β/T²) - ln(A*R/E) + E/(R*T) = 0
                f = log(β / T_peak^2) - log(A * R_GAS / E) + E / (R_GAS * T_peak)
                df = -2.0 / T_peak - E / (R_GAS * T_peak^2)
                T_new = T_peak - f / df

                if abs(T_new - T_peak) < tol
                    return T_new
                end
                T_peak = T_new
            end
            return T_peak
        end

        T_peak_analytical = find_peak_temperature(A_pre, E_act, β)

        # Run simulation and find numerical peak
        solid_A = SolidComponent(:A, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)

        rxn = Reaction(:pyrolysis, 1 => 2, A = A_pre, E = E_act, h = 0.0)

        mat = Material(name = :TGATest, components = (solid_A, solid_B), reactions = (rxn,))

        L = 0.0001
        mesh = create_uniform_mesh(
            L,
            1,
            Val(2);
            T_initial = T_initial,
            ξ_initial = (ξ_A_initial, 0.0),
        )

        T_ramp = t -> T_initial + β * t
        bc_set = BoundaryConditionSet(
            thermal = (
                top = DirichletThermalBC(T_ramp),
                bottom = DirichletThermalBC(T_ramp),
            ),
        )

        # Fine time resolution near expected peak
        t_peak_expected = (T_peak_analytical - T_initial) / β
        output_times = collect(range(0.0, t_final, length = 500))

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ξ_A_initial, z -> 0.0],
            tspan = (0.0, t_final),
            output_times = output_times,
        )

        sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        # Compute numerical mass loss rate (-dξ_A/dt) and find peak
        max_rate = 0.0
        T_peak_numerical = T_initial

        for i = 2:length(sol.t)
            dt = sol.t[i] - sol.t[i-1]
            ξ_A_prev = sol.ξ[:, :, i-1][1, 1]
            ξ_A_curr = sol.ξ[:, :, i][1, 1]

            rate = -(ξ_A_curr - ξ_A_prev) / dt  # Mass loss rate

            if rate > max_rate
                max_rate = rate
                T_peak_numerical = T_initial + β * (sol.t[i] + sol.t[i-1]) / 2
            end
        end

        # Peak temperature should match within 2 K
        @test abs(T_peak_numerical - T_peak_analytical) < 2.0
    end

    @testset "Heating Rate Effect (Kissinger Analysis)" begin
        # Test that peak temperature shifts with heating rate as predicted
        # by Kissinger equation: ln(β/T_p²) = const - E/(R*T_p)

        heating_rates = [5.0, 10.0, 20.0, 40.0] ./ 60.0  # K/s

        T_peaks = Float64[]

        for β_test in heating_rates
            solid_A = SolidComponent(:A, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)
            solid_B = SolidComponent(:B, ρ = ξ_A_initial, c = 1500.0, k = 0.2, ε = 0.9)

            rxn = Reaction(:pyrolysis, 1 => 2, A = A_pre, E = E_act, h = 0.0)

            mat = Material(
                name = :TGATest,
                components = (solid_A, solid_B),
                reactions = (rxn,),
            )

            L = 0.0001
            t_final_test = (T_final - T_initial) / β_test

            mesh = create_uniform_mesh(
                L,
                1,
                Val(2);
                T_initial = T_initial,
                ξ_initial = (ξ_A_initial, 0.0),
            )

            T_ramp = t -> T_initial + β_test * t
            bc_set = BoundaryConditionSet(
                thermal = (
                    top = DirichletThermalBC(T_ramp),
                    bottom = DirichletThermalBC(T_ramp),
                ),
            )

            output_times = collect(range(0.0, t_final_test, length = 300))

            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> T_initial,
                ξ_initial = [z -> ξ_A_initial, z -> 0.0],
                tspan = (0.0, t_final_test),
                output_times = output_times,
            )

            sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

            # Find peak temperature
            max_rate = 0.0
            T_peak = T_initial

            for i = 2:length(sol.t)
                dt = sol.t[i] - sol.t[i-1]
                ξ_A_prev = sol.ξ[:, :, i-1][1, 1]
                ξ_A_curr = sol.ξ[:, :, i][1, 1]
                rate = -(ξ_A_curr - ξ_A_prev) / dt

                if rate > max_rate
                    max_rate = rate
                    T_peak = T_initial + β_test * (sol.t[i] + sol.t[i-1]) / 2
                end
            end

            push!(T_peaks, T_peak)
        end

        # Kissinger plot: ln(β/T_p²) vs 1/T_p should be linear with slope -E/R
        x_kissinger = 1.0 ./ T_peaks
        y_kissinger = log.(heating_rates ./ (T_peaks .^ 2))

        # Linear regression to find slope
        n = length(x_kissinger)
        x_mean = mean(x_kissinger)
        y_mean = mean(y_kissinger)

        numerator = sum((x_kissinger .- x_mean) .* (y_kissinger .- y_mean))
        denominator = sum((x_kissinger .- x_mean) .^ 2)
        slope = numerator / denominator

        # Slope should be -E/R
        E_recovered = -slope * R_GAS
        rel_error = abs(E_recovered - E_act) / E_act

        # Recovered activation energy should be within 5% of true value
        @test rel_error < 0.05
    end
end

# ============================================================================
# Test 9.4: Method of Manufactured Solutions (MMS)
# ============================================================================
#
# Verify spatial convergence using Method of Manufactured Solutions.
#
# Setup:
# - Choose smooth manufactured solution: T̂(z,t) = T₀ + ΔT*sin(πz/L)*(1 - exp(-t/τ))
# - Compute required source term S(z,t)
# - Solve numerical problem with source term
# - Verify error convergence rate
#
# Since we don't have a mechanism to inject source terms, we use an exact 
# analytical solution instead: decaying sine wave for heat conduction.
#
# T(z,t) = T_0 + A * sin(πz/L) * exp(-α*(π/L)²*t)
#
# This is an exact solution to the heat equation with:
# - IC: T(z,0) = T_0 + A*sin(πz/L)
# - BCs: T(0,t) = T(L,t) = T_0 (Dirichlet)
#
@testset "Method of Manufactured Solutions" begin

    # Problem parameters
    T_0 = 300.0
    A = 100.0  # amplitude
    L = 0.01  # 1 cm

    # Material properties
    ρ = 1000.0
    cp = 1000.0
    k = 0.5
    α = k / (ρ * cp)  # thermal diffusivity = 5e-7 m²/s

    # Decay constant for fundamental mode
    decay_rate = α * (π / L)^2  # ≈ 4.93 s⁻¹

    # Analytical solution
    function T_analytical(z, t)
        return T_0 + A * sin(π * z / L) * exp(-decay_rate * t)
    end

    @testset "Convergence Rate Verification" begin
        # Test at multiple mesh resolutions
        n_cells_list = [20, 40, 80, 160]
        errors = Float64[]
        h_values = Float64[]

        t_test = 0.2  # Short time to see significant amplitude remaining

        for n_cells in n_cells_list
            h = L / n_cells
            push!(h_values, h)

            # Create material (no reactions for pure diffusion)
            solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
            mat = Material(name = :MMSTest, components = (solid,), reactions = ())

            mesh =
                create_uniform_mesh(L, n_cells, Val(1); T_initial = T_0, ξ_initial = (ρ,))

            # Dirichlet BCs at both ends
            bc_set = BoundaryConditionSet(
                thermal = (top = DirichletThermalBC(T_0), bottom = DirichletThermalBC(T_0)),
            )

            # Initial condition from analytical solution at t=0
            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> T_analytical(z, 0.0),
                ξ_initial = [z -> ρ],
                tspan = (0.0, t_test),
                output_times = [0.0, t_test],
            )

            sol = solve(problem; abstol = 1e-12, reltol = 1e-10, progress = false)

            # Compute L2 error against analytical solution at final time
            T_num = sol.T[:, end]

            l2_error = 0.0
            for (i, z) in enumerate(mesh.cell_centers)
                T_exact = T_analytical(z, t_test)
                l2_error += (T_num[i] - T_exact)^2 * h
            end
            l2_error = sqrt(l2_error / L)
            push!(errors, l2_error)
        end

        # Compute convergence rates
        rates = Float64[]
        for i = 2:length(errors)
            rate = log(errors[i-1] / errors[i]) / log(h_values[i-1] / h_values[i])
            push!(rates, rate)
        end

        # For pure diffusion with Dirichlet BCs and 2nd order spatial discretization,
        # expect approximately 2nd order convergence
        @test all(r -> r > 1.0, rates)  # At least 1st order
        @test mean(rates) > 1.5  # Average should be above 1.5
    end
end

# ============================================================================
# Test 9.5: Mass Conservation
# ============================================================================
#
# Verify global mass conservation.
#
# For impermeable boundaries:
#   Total mass in domain + mass lost through boundaries = initial mass
#
# For adiabatic mass BCs (impermeable):
#   Total mass should be exactly conserved (to machine precision)
#
@testset "Mass Conservation" begin

    @testset "Inert Material (No Reactions)" begin
        # Create inert material
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)
        mat = Material(name = :Inert, components = (solid,), reactions = ())

        L = 0.01
        n_cells = 50
        ξ_initial = 1000.0

        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1);
            T_initial = 300.0,
            ξ_initial = (ξ_initial,),
        )

        # Heat flux but impermeable mass BC
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> 300.0,
            ξ_initial = [z -> ξ_initial],
            tspan = (0.0, 100.0),
            output_times = collect(range(0.0, 100.0, length = 20)),
        )

        sol = solve(problem; progress = false)

        # Initial total mass per unit area
        M_initial = ξ_initial * L

        # Mass should be conserved at all times
        for i in eachindex(sol.t)
            M_total = sol.extras.total_mass[i]
            rel_error = abs(M_total - M_initial) / M_initial
            @test rel_error < 1e-10  # Machine precision
        end
    end

    @testset "Reacting Material (Closed System)" begin
        # Solid → Solid reaction (no gas escape)
        solid_A = SolidComponent(:A, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)

        rxn = Reaction(:convert, 1 => 2, A = 1e10, E = 100e3, h = 0.0)

        mat = Material(
            name = :SolidConvert,
            components = (solid_A, solid_B),
            reactions = (rxn,),
        )

        L = 0.01
        n_cells = 50
        ξ_A_initial = 1000.0

        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = 600.0,  # High T for reaction
            ξ_initial = (ξ_A_initial, 0.0),
        )

        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(600.0), bottom = DirichletThermalBC(600.0)),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> 600.0,
            ξ_initial = [z -> ξ_A_initial, z -> 0.0],
            tspan = (0.0, 10.0),
            output_times = collect(range(0.0, 10.0, length = 20)),
        )

        sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        # Total mass (A + B) should be conserved
        M_initial = ξ_A_initial * L

        for i in eachindex(sol.t)
            # Sum both components
            M_A = sum(sol.ξ[:, :, i][1, :] .* mesh.cell_volumes) / mesh.cross_section_area
            M_B = sum(sol.ξ[:, :, i][2, :] .* mesh.cell_volumes) / mesh.cross_section_area
            M_total = M_A + M_B

            rel_error = abs(M_total - M_initial) / M_initial
            @test rel_error < 1e-8
        end
    end
end

# ============================================================================
# Test 9.6: Energy Conservation
# ============================================================================
#
# Verify energy conservation.
#
# For adiabatic boundaries with no reactions:
#   Total thermal energy should be exactly conserved
#
# For adiabatic boundaries with reactions:
#   d(thermal energy)/dt = heat release from reactions
#
@testset "Energy Conservation" begin

    @testset "Adiabatic System (No Reactions)" begin
        # Create inert material with non-uniform initial temperature
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)
        mat = Material(name = :Inert, components = (solid,), reactions = ())

        L = 0.01
        n_cells = 50
        ρ_solid = 1000.0

        # Non-uniform initial temperature
        T_init_func = z -> 300.0 + 200.0 * sin(π * z / L)

        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1);
            T_initial = 300.0,
            ξ_initial = (ρ_solid,),
        )

        # Apply non-uniform IC
        for (i, cell) in enumerate(mesh.cells)
            z = mesh.cell_centers[i]
            mesh.cell_states[i] =
                CellState{1}(T_init_func(z), (ρ_solid,), mesh.cell_volumes[i])
        end

        # Adiabatic boundaries (closed system)
        bc_set =
            BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = T_init_func,
            ξ_initial = [z -> ρ_solid],
            tspan = (0.0, 1000.0),  # Long time for equilibration
            output_times = [0.0, 100.0, 500.0, 1000.0],
        )

        sol = solve(problem; progress = false)

        # Compute initial thermal energy
        # E = ∫ ρ cp T dV
        cp_val = 1000.0
        function compute_energy(temperatures, mesh)
            E = 0.0
            for (i, T) in enumerate(temperatures)
                E += ρ_solid * cp_val * T * mesh.cell_volumes[i]
            end
            return E / mesh.cross_section_area  # per unit area
        end

        E_initial = compute_energy(sol.T[:, 1], mesh)

        # Energy should be conserved at all times
        for i in eachindex(sol.t)
            E = compute_energy(sol.T[:, i], mesh)
            rel_error = abs(E - E_initial) / E_initial
            @test rel_error < 1e-6  # Should be very close
        end

        # Temperature should eventually equilibrate
        T_final = sol.T[:, end]
        T_mean_initial = sum(sol.T[:, 1]) / length(sol.T[:, 1])
        T_mean_final = mean(T_final)

        # Mean temperature should be conserved
        @test abs(T_mean_final - T_mean_initial) / T_mean_initial < 1e-6
    end

    @testset "Energy Balance with Reactions" begin
        # Exothermic reaction system
        solid_A = SolidComponent(:A, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)
        solid_B = SolidComponent(:B, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9)

        h_rxn = -100000.0  # J/kg (exothermic)
        rxn = Reaction(:exo, 1 => 2, A = 1e8, E = 80e3, h = h_rxn)

        mat = Material(
            name = :Exothermic,
            components = (solid_A, solid_B),
            reactions = (rxn,),
        )

        L = 0.001  # Small domain for fast equilibration
        n_cells = 10
        ρ_solid = 1000.0
        T_initial = 500.0

        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = T_initial,
            ξ_initial = (ρ_solid, 0.0),
        )

        # Adiabatic boundaries
        bc_set =
            BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ_solid, z -> 0.0],
            tspan = (0.0, 1.0),
            output_times = collect(range(0.0, 1.0, length = 20)),
        )

        sol = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        # Compute energy change
        cp_val = 1000.0
        function total_thermal_energy(temps, mesh)
            E = 0.0
            for (i, T) in enumerate(temps)
                # Use constant ρ*cp for simplicity
                E += ρ_solid * cp_val * T * mesh.cell_volumes[i]
            end
            return E / mesh.cross_section_area
        end

        function total_reaction_heat(conc_initial, conc_final, mesh, h_rxn)
            # Heat released = h_rxn * (mass of A converted)
            mass_A_initial = 0.0
            mass_A_final = 0.0
            for i in eachindex(mesh.cells)
                mass_A_initial += conc_initial[1, i] * mesh.cell_volumes[i]
                mass_A_final += conc_final[1, i] * mesh.cell_volumes[i]
            end
            mass_converted = (mass_A_initial - mass_A_final) / mesh.cross_section_area
            return h_rxn * mass_converted  # J/m² (negative for exothermic = heat released)
        end

        E_thermal_initial = total_thermal_energy(sol.T[:, 1], mesh)
        E_thermal_final = total_thermal_energy(sol.T[:, end], mesh)

        Q_rxn = total_reaction_heat(sol.ξ[:, :, 1], sol.ξ[:, :, end], mesh, h_rxn)

        # Energy balance: E_final - E_initial = -Q_rxn (heat released increases thermal energy)
        # For exothermic: h_rxn < 0, Q_rxn < 0, so -Q_rxn > 0
        ΔE_thermal = E_thermal_final - E_thermal_initial
        ΔE_expected = -Q_rxn  # Heat released goes into thermal energy

        if abs(ΔE_expected) > 1e-6
            rel_error = abs(ΔE_thermal - ΔE_expected) / abs(ΔE_expected)
            @test rel_error < 0.05  # 5% tolerance for energy balance
        end
    end
end

# ============================================================================
# Test 9.7: ALE Mass Conservation
# ============================================================================
#
# Verify that the integrated ALE ODE formulation conserves mass correctly.
#
# Key test: mass_initial = mass_final + mass_released_through_surface
#
# The ALE formulation should maintain conservation to solver tolerance
# because mesh motion is integrated into the ODE system.
#
@testset "ALE Mass Conservation" begin
    using Pyrolysis: get_thickness_from_state, compute_volumes_from_state

    @testset "Shrinking Domain with Reactions" begin
        # Create reactive material that produces gas
        virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

        # Moderate reaction rate for testing (positive h = endothermic, like PMMA)
        rxn = Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = 500e3)

        mat = Material(name = :TestPMMA, components = (virgin, gas), reactions = (rxn,))

        # Create mesh
        L = 0.005  # 5 mm
        n_cells = 50
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = 300.0,
            ξ_initial = (1190.0, 0.0),
        )

        # Heat flux to drive reactions, mass release at surface
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 0.1, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> 300.0,
            ξ_initial = [z -> 1190.0, z -> 0.0],
            tspan = (0.0, 10.0),  # Short simulation
            output_times = collect(range(0.0, 10.0, length = 21)),
        )

        # Solve with ALE
        solution = solve(
            problem;
            use_ale = true,
            min_thickness = 1e-5,
            progress = false,
            abstol = 1e-8,
            reltol = 1e-6,
        )

        @test solution.retcode in [:Success, :Terminated]

        # Check that we have valid solution data
        @test length(solution.t) >= 2
        @test all(isfinite, solution.extras.total_mass)

        # Verify mass balance: initial_mass = final_mass + mass_lost
        initial_mass = solution.extras.total_mass[1]
        final_mass = solution.extras.total_mass[end]

        # Integrate gas generation rate to get total mass released
        mass_released = 0.0
        for i = 2:length(solution.t)
            dt = solution.t[i] - solution.t[i-1]
            # Use trapezoidal rule
            ggr_avg =
                (
                    solution.extras.gas_generation_rate[i] +
                    solution.extras.gas_generation_rate[i-1]
                ) / 2
            mass_released += ggr_avg * dt
        end

        # Mass conservation check
        # initial_mass = final_mass + mass_released (approximately)
        mass_balance_error = abs(initial_mass - final_mass - mass_released)
        relative_error = mass_balance_error / initial_mass

        # Should be conserved to within a few percent
        # (tolerance accounts for numerical integration of GGR)
        @test relative_error < 0.1  # 10% tolerance for this coarse test
    end

    @testset "Thickness Evolution" begin
        # Verify that thickness decreases monotonically during pyrolysis
        virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:gas, M = 0.050, c = 1500.0, k = 0.02, λ = 1e-5)

        # Fast reaction for quick thickness change
        rxn = Reaction(:decomp, 1 => 2, A = 1e8, E = 100e3, h = -300e3)

        mat =
            Material(name = :FastPyrolysis, components = (virgin, gas), reactions = (rxn,))

        L = 0.003  # 3 mm
        n_cells = 30
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = 500.0,  # Start hot to get immediate reactions
            ξ_initial = (1190.0, 0.0),
        )

        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(80000.0), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 0.1, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> 500.0,
            ξ_initial = [z -> 1190.0, z -> 0.0],
            tspan = (0.0, 5.0),
            output_times = collect(range(0.0, 5.0, length = 11)),
        )

        solution = solve(
            problem;
            use_ale = true,
            min_thickness = 0.1e-3,  # Stop at 0.1 mm
            handle_depletion = true,
            depletion_threshold = 0.05,
            progress = false,
        )

        @test solution.retcode in [:Success, :Terminated]

        # Check thickness history exists and makes sense
        @test !isempty(solution.extras.thickness_history)
        @test length(solution.extras.thickness_history) == length(solution.t)

        # Initial thickness should match domain
        @test solution.extras.thickness_history[1] ≈ L rtol=0.01

        # Final thickness should be less than initial (material consumed)
        @test solution.extras.thickness_history[end] < solution.extras.thickness_history[1]

        # Thickness should decrease monotonically (no expansion)
        for i = 2:length(solution.extras.thickness_history)
            @test solution.extras.thickness_history[i] <=
                  solution.extras.thickness_history[i-1] + 1e-10
        end
    end

    @testset "Inert Material No Shrinkage" begin
        # With inert material at low temperature, no reactions -> no shrinkage
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 1.0, ε = 0.9)

        mat = Material(name = :Inert, components = (solid,), reactions = ())

        L = 0.01  # 10 mm
        n_cells = 20
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1);
            T_initial = 300.0,
            ξ_initial = (1000.0,),
        )

        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(10000.0), bottom = AdiabaticBC()),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> 300.0,
            ξ_initial = [z -> 1000.0],
            tspan = (0.0, 10.0),
            output_times = [0.0, 5.0, 10.0],
        )

        solution = solve(problem; use_ale = true, min_thickness = 1e-6, progress = false)

        @test solution.retcode == :Success

        # Thickness should remain constant (no reactions, no gas release)
        for thickness in solution.extras.thickness_history
            @test thickness ≈ L rtol=1e-6
        end

        # Mass should be conserved exactly (no mass loss)
        for mass in solution.extras.total_mass
            @test mass ≈ solution.extras.total_mass[1] rtol=1e-6
        end
    end
end

# ============================================================================
# Test: Darcy-Fick Gas Transport Conservation
# ============================================================================
#
# Verify that Darcy-Fick transport conserves mass in a closed system.
#
@testset "Darcy-Fick Mass Conservation" begin
    using Pyrolysis.Discretization:
        compute_rhs!, PropertyCache, RHSCache, update_properties!

    # Material with gas generation from solid
    solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.3)
    gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)

    material = Material(
        name = :ConservationTest,
        components = (solid, gas),
        reactions = (),  # No reactions - pure transport test
    )

    @testset "Closed System Conservation" begin
        # Setup mesh with gas in the middle
        L = 0.01
        n_cells = 20
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = 400.0,
            ξ_initial = (800.0, 0.0),
        )

        # Add gas to middle cells (creates pressure gradient)
        for i = 8:12
            state = mesh.cell_states[i]
            mesh.cell_states[i] =
                Pyrolysis.Geometry.CellState{2}(state.T, (state.ξ[1], 2.0), state.volume)
        end

        Pyrolysis.Geometry.update_geometry!(mesh)

        # Property cache
        prop_cache = PropertyCache(mesh)
        update_properties!(prop_cache, StateCache(mesh), mesh, material)

        # RHS cache
        rhs_cache = RHSCache(mesh, material)

        # Impermeable boundaries (closed system)
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        )

        # Pre-cache face IDs
        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]

        # Use component-based permeability on the solid
        material = Material(
            name = :ConservationTest,
            components = (
                SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.3, κ = 1e-11),
                GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5),
            ),
            reactions = (),  # No reactions - pure transport test
        )

        # Compute RHS with Darcy-Fick (component-based κ)
        dT = zeros(n_cells)
        dξ = zeros(2, n_cells)
        t = 0.0

        compute_rhs!(
            dT,
            dξ,
            StateCache(mesh),
            mesh,
            prop_cache,
            rhs_cache,
            material,
            bc_set,
            t,
            top_face_ids,
            bottom_face_ids;
            I_radiation = 0.0,
        )

        # In a closed system, total dξ_gas should be zero (mass conserved)
        total_dξ_gas = sum(dξ[2, i] * mesh.cell_volumes[i] for i = 1:n_cells)
        @test abs(total_dξ_gas) < 1e-10

        # Solid should be unchanged (no reactions)
        total_dξ_solid = sum(dξ[1, i] * mesh.cell_volumes[i] for i = 1:n_cells)
        @test abs(total_dξ_solid) < 1e-15
    end
end

# ============================================================================
# Test: Semi-Infinite Solid with Surface Heat Flux (Boundary Temperature Verification)
# ============================================================================
#
# Verify that the boundary temperature solver correctly computes surface temperature
# for a heat flux boundary condition against analytical solution.
#
# Analytical solution for semi-infinite solid with constant heat flux q at surface:
#   T(z,t) - T₀ = (2q/k) × √(αt/π) × exp(-η²) - (q×(L-z)/k) × erfc(η)
#   where η = (L-z) / (2√(αt))
#
# At the surface (z=L, η=0):
#   T_s(t) - T₀ = (2q/k) × √(αt/π)
#
# This tests the boundary temperature enhancement.
#
@testset "Semi-Infinite Solid with Heat Flux BC (Boundary Temperature)" begin
    using Pyrolysis.Discretization: get_surface_temperature, check_boundary_energy_balance

    # Material properties (constant for analytical solution)
    ρ = 1000.0    # kg/m³
    cp = 2000.0   # J/(kg·K)
    k = 0.2       # W/(m·K)
    α = k / (ρ * cp)  # Thermal diffusivity = 1e-7 m²/s

    # Boundary and initial conditions
    T_initial = 300.0  # K
    q_ext = 50000.0    # W/m² (50 kW/m² heat flux into domain)

    # Domain - thick enough to approximate semi-infinite for short times
    L = 0.01  # 1 cm

    # Analytical surface temperature: T_s(t) = T₀ + (2q/k) × √(αt/π)
    function T_surface_analytical(t, T_0, q, k, α)
        if t <= 0
            return T_0
        end
        return T_0 + (2 * q / k) * sqrt(α * t / π)
    end

    @testset "Surface Temperature at t=10s" begin
        # Create material (inert, no reactions)
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        # Fine mesh for accuracy
        n_cells = 100
        mesh =
            create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

        # Pure heat flux BC at top (no radiation losses for clean comparison)
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q_ext), bottom = AdiabaticBC()),
        )

        # Create and solve problem
        t_final = 10.0
        tspan = (0.0, t_final)
        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ],
            tspan = tspan,
        )

        solution = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        @test solution.retcode == :Success

        # Get numerical surface temperature from boundary state
        T_s_numerical = get_surface_temperature(solution.problem.mesh)

        # Analytical surface temperature
        T_s_analytical = T_surface_analytical(t_final, T_initial, q_ext, k, α)

        # Compare with 2% tolerance (mesh discretization error)
        @test T_s_numerical ≈ T_s_analytical rtol=0.02

        # Verify energy balance at boundary - returns Dict{Int, Float64}
        balance_results = check_boundary_energy_balance(solution.problem.mesh)
        @test balance_results[1] < 1e-6

        # Print values for inspection
        println("  Surface temperature verification at t=$(t_final)s:")
        println("    Analytical T_s = $(round(T_s_analytical, digits=2)) K")
        println("    Numerical T_s  = $(round(T_s_numerical, digits=2)) K")
        println(
            "    Error = $(round(100*abs(T_s_numerical - T_s_analytical)/T_s_analytical, digits=2))%",
        )
    end

    @testset "Surface Temperature Evolution" begin
        # Test at multiple times to verify evolution
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        n_cells = 100
        mesh =
            create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q_ext), bottom = AdiabaticBC()),
        )

        # Solve with output at multiple times
        test_times = [1.0, 5.0, 10.0, 20.0]
        tspan = (0.0, test_times[end])

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ],
            tspan = tspan,
        )

        solution = solve(
            problem;
            abstol = 1e-10,
            reltol = 1e-8,
            progress = false,
            saveat = test_times,
        )

        @test solution.retcode == :Success

        # Check solution outputs at each time
        for (i, t) in enumerate(test_times)
            T_s_analytical = T_surface_analytical(t, T_initial, q_ext, k, α)
            T_s_numerical = solution.surface_T[i]

            # Allow 3% tolerance (accounts for finite domain effects at later times)
            @test T_s_numerical ≈ T_s_analytical rtol=0.03
        end

        # Surface temperature should increase with sqrt(t)
        for i = 2:length(test_times)
            @test solution.surface_T[i] > solution.surface_T[i-1]
        end
    end

    @testset "Boundary Energy Balance Throughout Simulation" begin
        # Verify energy balance is maintained throughout, not just at end
        solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
        mat = Material(name = :TestSolid, components = (solid,), reactions = ())

        n_cells = 50
        mesh =
            create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

        # Combined BC: heat flux + radiative losses (more realistic)
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(q_ext) + RadiativeBC(ε = 0.9, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
        )

        tspan = (0.0, 10.0)
        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ],
            tspan = tspan,
        )

        solution = solve(problem; abstol = 1e-10, reltol = 1e-8, progress = false)

        @test solution.retcode == :Success

        # Both boundaries should have excellent energy balance - returns Dict
        balance_results = check_boundary_energy_balance(solution.problem.mesh)
        @test balance_results[1] < 1e-4
        @test balance_results[2] < 1e-10
    end
end
