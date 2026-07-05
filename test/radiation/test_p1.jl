"""
Unit tests for P1 volumetric radiation transport.

Tests the core P1 implementation:
- Radiation diffusion coefficient
- Source terms (emission, absorption)
- Boundary condition fluxes
- Radiation property updates
- Optical thickness diagnostics
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: StateCache
using Pyrolysis.Materials: STEFAN_BOLTZMANN, P1_QUASI_STEADY, AbstractRadiationBC
using Pyrolysis.Experimental:
    RadiationBCSet,
    MarshakRadiationBC,
    ReflectiveRadiationBC,
    DirichletRadiationBC,
    default_radiation_bcs
using Pyrolysis.Materials: evaluate_param, evaluate_param_T
using Pyrolysis.Properties:
    PropertyCache,
    RadiationCache,
    update_properties!,
    resize_radiation_cache!,
    initialize_radiation_intensity!,
    radiation_energy_density,
    net_radiation_source,
    total_radiation_flux
using Pyrolysis.Experimental:
    compute_radiation_diffusion_coefficient,
    compute_radiation_source_terms!,
    compute_radiation_fluxes!,
    radiation_bc_flux,
    optical_thickness,
    optical_thickness_cell,
    check_p1_validity,
    radiation_energy_balance
using Pyrolysis.Experimental: update_radiation_properties!

@testset "P1 Radiation" begin

    @testset "RadiationModel enum" begin
        @test Int(NO_RADIATION) == 0
        @test Int(SURFACE_ABSORPTION) == 1
        @test Int(BEER_LAMBERT) == 2
        @test Int(P1_QUASI_STEADY) == 3

        # Type check
        @test NO_RADIATION isa RadiationModel
        @test P1_QUASI_STEADY isa RadiationModel
    end

    @testset "Radiation diffusion coefficient" begin
        # D_rad = 1/(3α)

        # Typical virgin polymer (α ≈ 1500 m⁻¹)
        α = 1500.0
        D_rad = compute_radiation_diffusion_coefficient(α)
        @test D_rad ≈ 1.0 / (3.0 * 1500.0) atol=1e-12
        @test D_rad ≈ 2.222e-4 rtol=1e-3

        # Char layer (α ≈ 8000 m⁻¹)
        α = 8000.0
        D_rad = compute_radiation_diffusion_coefficient(α)
        @test D_rad ≈ 1.0 / (3.0 * 8000.0) atol=1e-12
        @test D_rad ≈ 4.167e-5 rtol=1e-3

        # Very optically thick (α = 50000 m⁻¹)
        α = 50000.0
        D_rad = compute_radiation_diffusion_coefficient(α)
        @test D_rad ≈ 1.0 / (3.0 * 50000.0) atol=1e-12

        # With scattering
        α = 5000.0
        σ_s = 1000.0
        D_rad = compute_radiation_diffusion_coefficient(α, σ_s)
        @test D_rad ≈ 1.0 / (3.0 * 6000.0) atol=1e-12

        # Transparent material (α → 0): clamped to a large finite value so the
        # quasi-steady Newton system stays regular (D = 0 would zero the G row).
        D_clamp = 1.0 / (3.0 * Pyrolysis.Experimental.P1_MIN_EXTINCTION)
        D_rad = compute_radiation_diffusion_coefficient(0.0)
        @test D_rad == D_clamp

        D_rad = compute_radiation_diffusion_coefficient(1e-15)
        @test D_rad == D_clamp  # Below threshold

        # Type stability
        @test @inferred(compute_radiation_diffusion_coefficient(1000.0)) isa Float64
        @test @inferred(compute_radiation_diffusion_coefficient(1000.0, 500.0)) isa Float64
    end

    @testset "Radiation boundary conditions" begin
        @testset "MarshakRadiationBC construction" begin
            # Simple constant BC
            bc = MarshakRadiationBC(I_external = 25e3, T_ambient = 300.0, emissivity = 0.9)
            @test bc.I_external == 25e3
            @test bc.T_ambient == 300.0
            @test bc.emissivity == 0.9
            @test bc isa AbstractRadiationBC
            @test bc isa MarshakRadiationBC

            # With defaults
            bc = MarshakRadiationBC(I_external = 50e3)
            @test bc.I_external == 50e3
            @test bc.T_ambient == 300.0  # default
            @test bc.emissivity == 0.9   # default

            # Time-dependent external flux
            bc = MarshakRadiationBC(
                I_external = t -> 25e3 * (1 - exp(-t/10)),
                T_ambient = 300.0,
            )
            @test bc.I_external isa Function
            @test bc.I_external(0.0) ≈ 0.0 atol=1e-10
            @test bc.I_external(100.0) ≈ 25e3 rtol=1e-3

            # Temperature-dependent emissivity
            bc = MarshakRadiationBC(
                I_external = 25e3,
                emissivity = T -> 0.8 + 0.0001*(T-300),
            )
            @test bc.emissivity isa Function
            @test bc.emissivity(300.0) ≈ 0.8
            @test bc.emissivity(500.0) ≈ 0.82

            # With validation
            @test_throws ArgumentError MarshakRadiationBC(
                I_external = -100.0,
                validate = true,
            )
            @test_throws ArgumentError MarshakRadiationBC(
                I_external = 25e3,
                T_ambient = -10.0,
                validate = true,
            )
            @test_throws ArgumentError MarshakRadiationBC(
                I_external = 25e3,
                emissivity = 1.5,
                validate = true,
            )
        end

        @testset "ReflectiveRadiationBC" begin
            bc = ReflectiveRadiationBC()
            @test bc isa AbstractRadiationBC
            @test bc isa ReflectiveRadiationBC
        end

        @testset "DirichletRadiationBC" begin
            # Constant G
            bc = DirichletRadiationBC(G = 1000.0)
            @test bc.G == 1000.0
            @test bc isa AbstractRadiationBC

            # Time-dependent G
            bc = DirichletRadiationBC(G = t -> 1000.0 + 100.0*t)
            @test bc.G isa Function
        end

        @testset "default_radiation_bcs" begin
            bcs = default_radiation_bcs(25e3, T_ambient = 300.0)
            @test bcs isa RadiationBCSet
            @test bcs.top isa MarshakRadiationBC
            @test bcs.bottom isa ReflectiveRadiationBC
            @test bcs.top.I_external == 25e3
        end
    end

    @testset "Parameter evaluation helpers" begin
        # Constant values
        @test evaluate_param(25e3, 0.0) == 25e3
        @test evaluate_param(300.0, 10.0) == 300.0

        # Functions of time
        f_t = t -> 100.0 * t
        @test evaluate_param(f_t, 5.0) == 500.0

        # Temperature-dependent
        @test evaluate_param_T(0.9, 500.0) == 0.9
        f_T = T -> 0.8 + 0.0001*T
        @test evaluate_param_T(f_T, 500.0) ≈ 0.85  # Use ≈ for floating point

        # Type stability
        @test @inferred(evaluate_param(100.0, 1.0)) isa Float64
        @test @inferred(evaluate_param_T(0.9, 500.0)) isa Float64
    end

    @testset "Marshak flux derivative coefficient" begin
        T_surf = 500.0
        G_surf = 12_345.0
        α_surf = 1000.0
        δG = 1e-3

        for ε in (0.5, 0.85, 0.9, 0.95, 1.0)
            bc = MarshakRadiationBC(I_external = 25e3, T_ambient = 300.0, emissivity = ε)
            coeff = Pyrolysis.Experimental._marshak_exchange_coefficient(bc, T_surf)
            fd =
                (
                    radiation_bc_flux(bc, G_surf + δG, T_surf, α_surf, ε, 0.0) -
                    radiation_bc_flux(bc, G_surf - δG, T_surf, α_surf, ε, 0.0)
                ) / (2δG)

            @test coeff ≈ ε / (2.0 * (2.0 - ε))
            @test fd ≈ -coeff rtol=1e-7 atol=1e-9
        end

        bc_T = MarshakRadiationBC(
            I_external = 25e3,
            T_ambient = 300.0,
            emissivity = T -> 0.8 + 0.0001 * (T - 300.0),
        )
        ε_T = bc_T.emissivity(T_surf)
        @test Pyrolysis.Experimental._marshak_exchange_coefficient(bc_T, T_surf) ≈
              ε_T / (2.0 * (2.0 - ε_T))
    end

    @testset "Half-cell effective Marshak coefficient" begin
        # c_eff = 1/(1/c + δ/D): the surface exchange c in series with the
        # half-cell diffusive conductance D/δ. Both the applied boundary flux
        # and the quasi-steady Newton Jacobian use this same coefficient.
        bc = MarshakRadiationBC(I_external = 0.0, T_ambient = 300.0, emissivity = 0.9)
        T = 500.0
        c = Pyrolysis.Experimental._marshak_exchange_coefficient(bc, T)

        D = 1.0 / (3.0 * 5000.0)
        δ = 5e-5
        c_eff = Pyrolysis.Experimental._marshak_effective_coefficient(bc, T, D, δ)
        @test c_eff ≈ 1.0 / (1.0 / c + δ / D)
        @test c_eff < c  # the half-cell resistance always reduces the exchange

        # Limits: vanishing half-cell distance or huge D recover the exact c
        @test Pyrolysis.Experimental._marshak_effective_coefficient(bc, T, D, 1e-300) ≈ c
        @test Pyrolysis.Experimental._marshak_effective_coefficient(bc, T, 1e10, δ) ≈ c
    end

    @testset "RadiationCache" begin
        @testset "Construction" begin
            # From cell count
            n_cells = 50
            cache = RadiationCache(n_cells)
            @test length(cache.D_rad) == 50
            @test length(cache.D_rad_face) == 51  # n_cells + 1
            @test length(cache.G) == 50
            @test length(cache.q_rad) == 51
            @test length(cache.Q_emission) == 50
            @test length(cache.Q_absorption) == 50

            # All initialized to zero
            @test all(cache.D_rad .== 0.0)
            @test all(cache.G .== 0.0)
        end

        @testset "Resize" begin
            cache = RadiationCache(20)
            @test length(cache.G) == 20

            resize_radiation_cache!(cache, 50)
            @test length(cache.G) == 50
            @test length(cache.q_rad) == 51
            @test length(cache.D_rad_face) == 51
        end

        @testset "Initialize radiation intensity" begin
            cache = RadiationCache(10)

            # From ambient temperature
            T_amb = 300.0
            initialize_radiation_intensity!(cache, T_amb)

            G_expected = 4.0 * STEFAN_BOLTZMANN * T_amb^4
            @test all(cache.G .≈ G_expected)
            @test all(cache.q_rad .== 0.0)
            @test all(cache.Q_emission .== 0.0)

            # Higher temperature
            T_high = 800.0
            initialize_radiation_intensity!(cache, T_high)
            G_high = 4.0 * STEFAN_BOLTZMANN * T_high^4
            @test all(cache.G .≈ G_high)
            @test G_high > G_expected  # Higher T means more radiation
        end

        @testset "Diagnostic functions" begin
            cache = RadiationCache(10)
            cache.G[5] = 1000.0
            cache.Q_emission[5] = 500.0
            cache.Q_absorption[5] = 800.0

            # Radiation energy density
            c_light = 2.998e8
            @test radiation_energy_density(cache, 5) ≈ 1000.0 / c_light

            # Net radiation source
            @test net_radiation_source(cache, 5) == 800.0 - 500.0  # absorption - emission
        end
    end

    @testset "Marshak BC flux computation" begin
        # Setup: surface at 500 K, ambient at 300 K, external flux 25 kW/m²
        bc = MarshakRadiationBC(I_external = 25e3, T_ambient = 300.0, emissivity = 0.9)

        T_surf = 500.0
        G_surf = 1000.0  # Some incident radiation
        α_surf = 5000.0  # absorption coefficient
        ε_surf = 0.9
        t = 0.0

        flux = radiation_bc_flux(bc, G_surf, T_surf, α_surf, ε_surf, t)

        # The Marshak relation derived from the P1 half-range fluxes is
        #   q = c · (G_ext + 4·I_ext − G_surf),  c = ε/(2(2−ε)),  G_ext = 4σT_amb⁴
        # The external flux carries the SAME half-range weight 4c = 2ε/(2−ε)
        # as the ambient term (2.0 at ε = 1, 1.64 at ε = 0.9) — a bare I_ext
        # (weight 1) under-injects the external flux by that factor.
        σ = STEFAN_BOLTZMANN
        G_ext = 4.0 * σ * 300.0^4
        marshak_coef(ε) = ε / (2.0 * (2.0 - ε))
        expected_flux = marshak_coef(0.9) * (G_ext + 4.0 * 25e3 - G_surf)

        @test flux ≈ expected_flux rtol=1e-10

        # Regression: the I_ext response is exactly 4c·I_ext = 2ε/(2−ε)·I_ext
        bc0 = MarshakRadiationBC(I_external = 0.0, T_ambient = 300.0, emissivity = 0.9)
        flux0 = radiation_bc_flux(bc0, G_surf, T_surf, α_surf, ε_surf, t)
        @test flux - flux0 ≈ 4.0 * marshak_coef(0.9) * 25e3 rtol=1e-12
        for ε in (0.5, 0.9, 1.0)
            bcA = MarshakRadiationBC(I_external = 1e4, T_ambient = 300.0, emissivity = ε)
            bcB = MarshakRadiationBC(I_external = 0.0, T_ambient = 300.0, emissivity = ε)
            ΔI =
                radiation_bc_flux(bcA, G_surf, T_surf, α_surf, ε, t) -
                radiation_bc_flux(bcB, G_surf, T_surf, α_surf, ε, t)
            @test ΔI ≈ 2.0 * ε / (2.0 - ε) * 1e4 rtol=1e-12
        end

        # Equilibrium: when G_surf = G_ext and I_ext = 0, flux should be zero
        # (independent of the coefficient form, as long as the (G_ext − G_surf)
        # factor is preserved).
        bc_eq = MarshakRadiationBC(I_external = 0.0, T_ambient = 300.0, emissivity = 1.0)
        G_eq = 4.0 * σ * 300.0^4
        flux_eq = radiation_bc_flux(bc_eq, G_eq, 300.0, 5000.0, 1.0, 0.0)
        @test abs(flux_eq) < 1e-6  # Zero at equilibrium

        # Type stability
        @test @inferred(radiation_bc_flux(bc, G_surf, T_surf, α_surf, ε_surf, t)) isa
              Float64
    end

    @testset "Reflective BC flux" begin
        bc = ReflectiveRadiationBC()
        flux = radiation_bc_flux(bc, 1000.0, 500.0, 5000.0, 0.9, 0.0)
        @test flux == 0.0

        @test @inferred(radiation_bc_flux(bc, 1000.0, 500.0, 5000.0, 0.9, 0.0)) isa Float64
    end

    @testset "Source terms computation" begin
        # Create simple test material with known properties
        # Using a simple solid component
        # Note: α in SolidComponent is MASS absorption coefficient [m²/kg]
        # Effective α = ξ * α_mass, so with ξ=1000 kg/m³ and α_mass=5.0 m²/kg
        # we get α_eff = 5000 m⁻¹
        α_mass = 5.0  # m²/kg
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = α_mass)
        material = Material(name = :TestMat, components = (solid,), reactions = ())

        # Create mesh
        n_cells = 10
        L = 0.01  # 1 cm
        ξ_solid = 1000.0  # kg/m³
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 500.0,
            ξ_initial = (ξ_solid,),
        )

        # Create caches
        prop_cache = PropertyCache(mesh)
        rad_cache = RadiationCache(n_cells)

        # Set radiation intensity
        fill!(rad_cache.G, 2000.0)  # W/m²

        # Update properties
        update_properties!(prop_cache, StateCache(mesh, material), mesh, material)

        # Check that α is correctly populated
        # α_eff = ξ * α_mass = 1000 * 5.0 = 5000 m⁻¹
        α_eff = ξ_solid * α_mass
        @test all(prop_cache.α .≈ α_eff)
        @test all(prop_cache.ε .≈ 0.9)

        # Compute source terms
        compute_radiation_source_terms!(
            rad_cache.Q_emission,
            rad_cache.Q_absorption,
            mesh,
            rad_cache,
            prop_cache,
        )

        # Verify emission: Q_emission = 4σαT⁴ (no emissivity factor - blackbody inside)
        # Note: In P1 approximation, volumetric emission is blackbody (ε=1).
        # Emissivity only appears at surface boundaries (Marshak BC).
        T = 500.0
        Q_em_expected = 4.0 * STEFAN_BOLTZMANN * α_eff * T^4
        @test all(rad_cache.Q_emission .≈ Q_em_expected)

        # Verify absorption: Q_absorption = αG
        G = 2000.0
        Q_abs_expected = α_eff * G
        @test all(rad_cache.Q_absorption .≈ Q_abs_expected)

        # Verify dimensions: W/m³
        # At T=500K, Q_emission ≈ 4 * 5.67e-8 * 0.9 * 5000 * 500^4 ≈ 6.4e7 W/m³
        @test rad_cache.Q_emission[1] > 1e6  # Should be large

        # At G=2000, Q_absorption = 5000 * 2000 = 1e7 W/m³
        @test rad_cache.Q_absorption[1] ≈ 1e7
    end

    @testset "Radiation property update" begin
        # Create simple test setup
        # α_mass = 5.0 m²/kg, ξ = 1000 kg/m³ → α_eff = 5000 m⁻¹
        α_mass = 5.0
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = α_mass)
        material = Material(name = :TestMat, components = (solid,), reactions = ())

        n_cells = 5
        L = 0.01
        ξ_solid = 1000.0
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )

        prop_cache = PropertyCache(mesh)
        rad_cache = RadiationCache(n_cells)

        # Update thermal properties first
        update_properties!(prop_cache, StateCache(mesh, material), mesh, material)

        # Update radiation properties
        update_radiation_properties!(rad_cache, mesh, prop_cache)

        # Check D_rad = 1/(3α)
        α_eff = ξ_solid * α_mass  # 5000 m⁻¹
        D_rad_expected = 1.0 / (3.0 * α_eff)
        @test all(rad_cache.D_rad .≈ D_rad_expected)

        # Check face values (interior faces should be harmonic mean, 
        # but with uniform α, same value)
        # For uniform mesh and properties, all face values should be same
        @test all(rad_cache.D_rad_face .≈ D_rad_expected)
    end

    @testset "Radiation flux computation" begin
        # Create test setup
        # α_mass = 5.0 m²/kg, ξ = 1000 kg/m³ → α_eff = 5000 m⁻¹
        α_mass = 5.0
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = α_mass)
        material = Material(name = :TestMat, components = (solid,), reactions = ())

        n_cells = 10
        L = 0.01
        ξ_solid = 1000.0
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )

        prop_cache = PropertyCache(mesh)
        rad_cache = RadiationCache(n_cells)

        # Set up a linear G profile: G increases from bottom to top
        for i = 1:n_cells
            z = mesh.cell_centers[i]
            rad_cache.G[i] = 1000.0 + 50000.0 * (z / L)  # Linear increase
        end

        # Update properties
        update_properties!(prop_cache, StateCache(mesh, material), mesh, material)
        update_radiation_properties!(rad_cache, mesh, prop_cache)

        # Compute fluxes
        radiation_bcs = RadiationBCSet(
            top = MarshakRadiationBC(I_external = 25e3, T_ambient = 300.0),
            bottom = ReflectiveRadiationBC(),
        )

        compute_radiation_fluxes!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

        # Check interior fluxes: q = -D_rad * dG/dx
        # With positive dG/dx, flux should be negative (toward low G)
        for f = 2:n_cells  # Interior faces
            face = mesh.faces[f]
            left_id, right_id = face.cell_ids
            if left_id > 0 && right_id > 0
                @test rad_cache.q_rad[f] < 0  # Negative flux (downward)
            end
        end

        # Reflective BC at bottom: flux should be zero
        @test rad_cache.q_rad[1] == 0.0
    end

    @testset "Optical thickness" begin
        # Create test setup with known α
        # α_mass = 5.0 m²/kg, ξ = 1000 kg/m³ → α_eff = 5000 m⁻¹
        # With L = 0.002 m (2 mm), τ = α_eff × L = 5000 × 0.002 = 10
        α_mass = 5.0
        ξ_solid = 1000.0
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = α_mass)
        material = Material(name = :TestMat, components = (solid,), reactions = ())

        n_cells = 10
        L = 0.002  # 2 mm
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )

        prop_cache = PropertyCache(mesh)
        update_properties!(prop_cache, StateCache(mesh, material), mesh, material)

        # α_eff = ξ × α_mass = 1000 × 5.0 = 5000 m⁻¹
        # τ = α_eff × L = 5000 × 0.002 = 10
        α_eff = ξ_solid * α_mass
        τ = optical_thickness(mesh, prop_cache)
        @test τ ≈ α_eff * L rtol=1e-6

        # Cell optical thickness
        Δx = L / n_cells
        τ_cell = optical_thickness_cell(mesh, prop_cache, 5)
        @test τ_cell ≈ α_eff * Δx rtol=1e-6

        # Validity check
        validity = check_p1_validity(mesh, prop_cache, warn = false)
        @test validity.τ_total ≈ α_eff * L rtol=1e-6
        @test validity.valid == true
        @test occursin("optically thick", validity.regime)
    end

    @testset "P1 validity regimes" begin
        # Use mass absorption coefficients that give desired optical thicknesses
        # with ξ = 1000 kg/m³ and L = 0.01 m (1 cm)
        # τ = ξ × α_mass × L = 1000 × α_mass × 0.01 = 10 × α_mass
        ξ_solid = 1000.0
        L = 0.01  # 1 cm
        n_cells = 10

        # Optically thick: τ = 100 → α_mass = 10 m²/kg
        solid_thick =
            SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 10.0)
        # Intermediate: τ = 5 → α_mass = 0.5 m²/kg
        solid_intermediate =
            SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.5)
        # Optically thin: τ = 0.5 → α_mass = 0.05 m²/kg
        solid_thin =
            SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.05)
        # Transparent: τ = 0.05 → α_mass = 0.005 m²/kg
        solid_transparent =
            SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.005)

        # Optically thick: τ = 10 × 10.0 = 100
        mat_thick = Material(name = :Thick, components = (solid_thick,), reactions = ())
        mesh_thick = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )
        prop_thick = PropertyCache(mesh_thick)
        update_properties!(
            prop_thick,
            StateCache(mesh_thick, mat_thick),
            mesh_thick,
            mat_thick,
        )
        validity = check_p1_validity(mesh_thick, prop_thick, warn = false)
        @test validity.τ_total ≈ 100.0 rtol=1e-3
        @test validity.valid == true
        @test occursin("optically thick", validity.regime)

        # Intermediate: τ = 10 × 0.5 = 5
        mat_int = Material(
            name = :Intermediate,
            components = (solid_intermediate,),
            reactions = (),
        )
        mesh_int = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )
        prop_int = PropertyCache(mesh_int)
        update_properties!(prop_int, StateCache(mesh_int, mat_int), mesh_int, mat_int)
        validity = check_p1_validity(mesh_int, prop_int, warn = false)
        @test 1.0 < validity.τ_total < 10.0
        @test validity.valid == true
        @test occursin("intermediate", validity.regime)

        # Optically thin: τ = 10 × 0.05 = 0.5
        mat_thin = Material(name = :Thin, components = (solid_thin,), reactions = ())
        mesh_thin = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )
        prop_thin = PropertyCache(mesh_thin)
        update_properties!(prop_thin, StateCache(mesh_thin, mat_thin), mesh_thin, mat_thin)
        validity = check_p1_validity(mesh_thin, prop_thin, warn = false)
        @test 0.1 < validity.τ_total < 1.0
        @test validity.valid == true
        @test occursin("optically thin", validity.regime)

        # Transparent: τ = 10 × 0.005 = 0.05
        mat_trans =
            Material(name = :Transparent, components = (solid_transparent,), reactions = ())
        mesh_trans = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 300.0,
            ξ_initial = (ξ_solid,),
        )
        prop_trans = PropertyCache(mesh_trans)
        update_properties!(
            prop_trans,
            StateCache(mesh_trans, mat_trans),
            mesh_trans,
            mat_trans,
        )
        validity = check_p1_validity(mesh_trans, prop_trans, warn = false)
        @test validity.τ_total < 0.1
        @test validity.valid == false
        @test occursin("transparent", validity.regime)
    end

    @testset "Radiation energy balance" begin
        # Create test setup
        # α_mass = 5.0 m²/kg, ξ = 1000 kg/m³ → α_eff = 5000 m⁻¹
        α_mass = 5.0
        ξ_solid = 1000.0
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = α_mass)
        material = Material(name = :TestMat, components = (solid,), reactions = ())

        n_cells = 20
        L = 0.01
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(1),
            T_initial = 500.0,
            ξ_initial = (ξ_solid,),
        )

        prop_cache = PropertyCache(mesh)
        rad_cache = RadiationCache(n_cells)

        # Initialize G and compute sources
        initialize_radiation_intensity!(rad_cache, 500.0)
        update_properties!(prop_cache, StateCache(mesh, material), mesh, material)

        compute_radiation_source_terms!(
            rad_cache.Q_emission,
            rad_cache.Q_absorption,
            mesh,
            rad_cache,
            prop_cache,
        )

        # Compute balance
        balance = radiation_energy_balance(rad_cache, mesh)

        @test haskey(balance, :total_emission)
        @test haskey(balance, :total_absorption)
        @test haskey(balance, :net_source)
        @test haskey(balance, :flux_top)
        @test haskey(balance, :flux_bottom)
        @test haskey(balance, :balance_error)

        # At equilibrium (G = 4σT⁴ with ε=1), emission ≈ absorption
        # With ε=0.9, there will be some imbalance
        @test balance.total_emission > 0
        @test balance.total_absorption > 0
    end

end  # @testset "P1 Radiation"
