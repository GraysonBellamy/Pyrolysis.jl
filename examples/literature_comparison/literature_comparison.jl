"""
Literature Property Comparison — Douglas Fir Cone Calorimeter

Runs 10 simulations (1 baseline + 9 variants) with literature-sourced property
values substituted one category at a time, using the sensitivity analysis
baseline as the reference configuration (240 cells, 12 mm slab, 50 kW/m²,
6% MC, 600 s).

Variants:
  0. Baseline — our measured values
  1. Virgin cp — Harada et al. (1998), low
  2. Virgin cp — Dupont et al. (2014), high
  3. Char cp — Hostikka & Matala (2017), low
  4. Char cp — Edrich (1980), high
  5. Heats of reaction — Jones et al. (2020), ΔH_P = +306.92 kJ/kg (uniform, f=2.654);
     Jones with-lid linear fit evaluated at Y_c = 0.17680 (our DF TGA model at
     Jones's 5 K/min, 823 K endpoint) — see literature_comparison_plan.md for derivation
  7. Heat of sorption — Free water only (2442 kJ/kg)
  8. Heats of reaction — Rinta-Paavola et al. (2023), ΔH_P = +112 kJ/kg (uniform, f=0.969)
  9. Heats of reaction — Rath et al. (2003), two-zone (f_P=1.234, f_S=0.638);
     delivers Rath's H_P + H_S = +179.46 kJ/kg at Y_c = 0.20741 (our DF TGA model at
     Rath's 10 K/min, 773 K endpoint) — see literature_comparison_plan.md for derivation

Outputs (saved to results/):
  - {label}_results.csv per run (time, mlr, surface_temp, back_temp, total_mass, thickness)
  - literature_comparison_qois.csv summary table

Usage:
  julia -p 10 --project=. examples/literature_comparison/literature_comparison.jl
  julia --project=. examples/literature_comparison/literature_comparison.jl   # serial
"""

using Distributed
using Printf

println("Running with $(nworkers()) worker(s)")

# Load Printf on all workers before the main @everywhere block so that
# @printf is available when the block is lowered (macro-expanded).
@everywhere using Printf

# ============================================================================
# Load packages and define all functions on every worker
# ============================================================================

@everywhere begin
    using Pyrolysis
    using Pyrolysis: compute_surface_mass_flux, effective_gas_transfer
    using OrdinaryDiffEq
    using OrdinaryDiffEqSDIRK: KenCarp4   # v7 split: KenCarp4 not re-exported by the umbrella pkg
    using LinearSolve
    using Sparspak
    using LinearSolve: SparspakFactorization
    using Printf

    # ========================================================================
    # Constants
    # ========================================================================

    const LC_FINAL_CHAR_YIELD = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
    const LC_ALPHA = [
        (1 - 0.95) / (1 - LC_FINAL_CHAR_YIELD),
        (1 - 0.95 * 0.80) / (1 - LC_FINAL_CHAR_YIELD),
        (1 - 0.95 * 0.80 * 0.37) / (1 - LC_FINAL_CHAR_YIELD),
        (1 - 0.95 * 0.80 * 0.37 * 0.79) / (1 - LC_FINAL_CHAR_YIELD),
    ]
    const LC_MOISTURE_IDX = 7
    const LC_ρ_DRY_BASIS = 530.0
    const LC_NOMINAL_HRXN = [
        -19200.66288032251,
        -5263.201466602739,
        264424.79037373835,
        -217151.69043762444,
        -326.40885469229025,
    ]

    # ========================================================================
    # Baseline property functions
    # ========================================================================

    const LC_T_HFM_MAX = 323.16
    const LC_T_STA_MIN = 425.02219679025956
    const LC_T_INTERP_WIDTH = LC_T_STA_MIN - LC_T_HFM_MAX

    function baseline_virgin_cp(T)
        if T <= LC_T_HFM_MAX
            -18.01630709338176 + 4.133836600273531 * T
        elseif T <= LC_T_STA_MIN
            cp_low = -18.01630709338176 + 4.133836600273531 * LC_T_HFM_MAX
            cp_high = -2139.3044552435705 + 9.136799235573225 * LC_T_STA_MIN
            cp_low + (cp_high - cp_low) * (T - LC_T_HFM_MAX) / LC_T_INTERP_WIDTH
        else
            -2139.3044552435705 + 9.136799235573225 * T
        end
    end

    baseline_char_cp(T) = 312.4273785167445 + 1.8372509265710314 * T

    function baseline_hos(T, ξ)
        MC = 100.0 * ξ[LC_MOISTURE_IDX] / LC_ρ_DRY_BASIS
        return (2450.420909987441 + 875.6441221735336 * exp(-0.1121684273316275 * MC)) * 1e3
    end

    # ========================================================================
    # Literature substitution functions
    # ========================================================================

    # Virgin cp — Harada et al. (1998): cp = 1.200 + 2.45e-3*(T-273.15) kJ/(kg·K)
    harada_virgin_cp(T) = 1200.0 + 2.45 * (T - 273.15)

    # Virgin cp — Dupont et al. (2014): cp = 5.2e-3*T - 0.2847 kJ/(kg·K)
    dupont_virgin_cp(T) = 5.2 * T - 284.7

    # Char cp — Hostikka & Matala (2017): cp ≈ -0.370 + 2.122e-3*T kJ/(kg·K)
    hostikka_char_cp(T) = -370.0 + 2.122 * T

    # Char cp — Edrich (1980): cp = -0.29227 + 4.1387e-3*T - 1.58261e-6*T² kJ/(kg·K)
    edrich_char_cp(T) = -292.27 + 4.1387 * T - 1.58261e-3 * T^2

    # HoS — Free water only: constant 2442 kJ/kg = 2.442e6 J/kg
    free_water_hos(T, ξ) = 2.442e6

    # ========================================================================
    # Symbol → function dispatch (avoids serializing closures across processes)
    # ========================================================================

    const VIRGIN_CP_DISPATCH = Dict{Symbol,Function}(
        :baseline => baseline_virgin_cp,
        :harada => harada_virgin_cp,
        :dupont => dupont_virgin_cp,
    )
    const CHAR_CP_DISPATCH = Dict{Symbol,Function}(
        :baseline => baseline_char_cp,
        :hostikka => hostikka_char_cp,
        :edrich => edrich_char_cp,
    )
    const HOS_DISPATCH =
        Dict{Symbol,Function}(:baseline => baseline_hos, :free_water => free_water_hos)

    # ========================================================================
    # Run specification — fully serializable (no closures)
    # ========================================================================

    struct LitCompRunSpec
        label::String
        virgin_cp_key::Symbol
        char_cp_key::Symbol
        hrxn_scale::Float64      # uniform scale: multiply all NOMINAL_HRXN by this
        hrxn_scale_secondary::Float64  # two-zone: secondary zone (rxns 4-5) scale; NaN = use hrxn_scale
        hos_key::Symbol
    end

    # ========================================================================
    # Simulation builder
    # ========================================================================

    function build_and_solve(spec::LitCompRunSpec)
        virgin_cp = VIRGIN_CP_DISPATCH[spec.virgin_cp_key]
        char_cp = CHAR_CP_DISPATCH[spec.char_cp_key]
        hos_func = HOS_DISPATCH[spec.hos_key]

        # Two-zone scaling: primary (rxns 1-3) vs secondary (rxns 4-5)
        f_p = spec.hrxn_scale
        f_s = isnan(spec.hrxn_scale_secondary) ? spec.hrxn_scale : spec.hrxn_scale_secondary
        hrxn = [
            LC_NOMINAL_HRXN[1] * f_p,
            LC_NOMINAL_HRXN[2] * f_p,
            LC_NOMINAL_HRXN[3] * f_p,
            LC_NOMINAL_HRXN[4] * f_s,
            LC_NOMINAL_HRXN[5] * f_s,
        ]

        # Intermediate cp's via cumulative-yield interpolation
        int1_cp = CallableProperty(
            T -> (1 - LC_ALPHA[1]) * virgin_cp(T) + LC_ALPHA[1] * char_cp(T),
        )
        int2_cp = CallableProperty(
            T -> (1 - LC_ALPHA[2]) * virgin_cp(T) + LC_ALPHA[2] * char_cp(T),
        )
        int3_cp = CallableProperty(
            T -> (1 - LC_ALPHA[3]) * virgin_cp(T) + LC_ALPHA[3] * char_cp(T),
        )
        int4_cp = CallableProperty(
            T -> (1 - LC_ALPHA[4]) * virgin_cp(T) + LC_ALPHA[4] * char_cp(T),
        )

        df_cp = CallableProperty(virgin_cp)
        char_cp_prop = CallableProperty(char_cp)

        # Moisture cp (NIST Shomate, fixed)
        moisture_cp = CallableProperty(T -> begin
            t = T / 1000.0
            M = 0.01801528
            A = -203.6060
            B = 1523.290
            C = -3196.413
            D = 2474.455
            E = 3.855326
            (A + B*t + C*t^2 + D*t^3 + E/t^2) / M
        end)

        # Water vapor cp (fixed)
        water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

        # Heat of sorption
        hos = StateDependentProperty(hos_func)

        # Gas transfer / permeability constants (physically grounded)
        λ_wood = 5.0e-7     # m²/s - virgin wood (ε~0.05, τ~5)
        λ_char = 1.0e-5     # m²/s - porous char (ε~0.7, τ~2)
        λ_gas = 3.0e-5     # m²/s - reference free diffusion (M~29, ~500K)
        κ_virgin = 1.0e-15    # m² - tight virgin wood
        κ_char = 1.0e-12    # m² - open char pore network

        material = Material(
            name = :DouglasFir,
            components = (
                SolidComponent(
                    :DF,
                    ρ = 530.0,
                    c = df_cp,
                    k = 0.145,
                    ε = 0.81,
                    α = 0.0,
                    λ = λ_wood,
                    κ = κ_virgin,
                ),
                SolidComponent(
                    :DF_int1,
                    ρ = 507.2,
                    c = int1_cp,
                    k = 0.292951,
                    ε = 0.790965,
                    α = 0.0,
                    λ = λ_wood * 2.0,
                    κ = κ_virgin * 3.0,
                ),
                SolidComponent(
                    :DF_int2,
                    ρ = 419.8,
                    c = int2_cp,
                    k = 0.329938,
                    ε = 0.786207,
                    α = 0.0,
                    λ = λ_wood * 5.0,
                    κ = κ_virgin * 10.0,
                ),
                SolidComponent(
                    :DF_int3,
                    ρ = 203.2,
                    c = int3_cp,
                    k = 0.644833,
                    ε = 0.745694,
                    α = 0.0,
                    λ = λ_char * 0.3,
                    κ = κ_char * 0.01,
                ),
                SolidComponent(
                    :DF_int4,
                    ρ = 176.6,
                    c = int4_cp,
                    k = 0.777700,
                    ε = 0.728600,
                    α = 0.0,
                    λ = λ_char * 0.7,
                    κ = κ_char * 0.1,
                ),
                SolidComponent(
                    :Char,
                    ρ = 150.0,
                    c = char_cp_prop,
                    k = 1.0,
                    ε = 0.70,
                    α = 0.0,
                    λ = λ_char,
                    κ = κ_char,
                ),
                LiquidComponent(
                    :MOISTURE,
                    ρ = 1e6,
                    c = moisture_cp,
                    k = 0.145,
                    ε = 0.81,
                    α = 0.0,
                    κ = 0.0,
                ),
                GaseousComponent(
                    :DF_gas1,
                    M = 0.060,
                    c = 1400.0,
                    k = 0.020,
                    λ = λ_gas * sqrt(0.029/0.060),
                ),
                GaseousComponent(
                    :DF_gas2,
                    M = 0.070,
                    c = 1400.0,
                    k = 0.020,
                    λ = λ_gas * sqrt(0.029/0.070),
                ),
                GaseousComponent(
                    :DF_gas3,
                    M = 0.110,
                    c = 1300.0,
                    k = 0.018,
                    λ = λ_gas * sqrt(0.029/0.110),
                ),
                GaseousComponent(
                    :DF_gas4,
                    M = 0.070,
                    c = 1500.0,
                    k = 0.020,
                    λ = λ_gas * sqrt(0.029/0.070),
                ),
                GaseousComponent(
                    :DF_gas5,
                    M = 0.050,
                    c = 1600.0,
                    k = 0.025,
                    λ = λ_gas * sqrt(0.029/0.050),
                ),
                GaseousComponent(
                    :WATER_VAPOR,
                    M = 0.018,
                    c = water_vapor_cp,
                    k = 0.025,
                    λ = λ_gas * sqrt(0.029/0.018),
                ),
            ),
            reactions = (
                Reaction(
                    :pyrolysis_step1,
                    1 => (2, 8),
                    A = 74816.95,
                    E = 71.59e3,
                    h = hrxn[1],
                    n = 1.0,
                    yields = (0.95, 0.05),
                ),
                Reaction(
                    :pyrolysis_step2,
                    2 => (3, 9),
                    A = 23550492839.0,
                    E = 140.9e3,
                    h = hrxn[2],
                    n = 1.0,
                    yields = (0.80, 0.20),
                ),
                Reaction(
                    :pyrolysis_step3,
                    3 => (4, 10),
                    A = 7.744618e12,
                    E = 181.0e3,
                    h = hrxn[3],
                    n = 1.0,
                    yields = (0.37, 0.63),
                ),
                Reaction(
                    :pyrolysis_step4,
                    4 => (5, 11),
                    A = 103276.14,
                    E = 98.94e3,
                    h = hrxn[4],
                    n = 1.0,
                    yields = (0.79, 0.21),
                ),
                Reaction(
                    :pyrolysis_step5,
                    5 => (6, 12),
                    A = 6.397348,
                    E = 56.01e3,
                    h = hrxn[5],
                    n = 1.0,
                    yields = (0.74, 0.26),
                ),
                Reaction(
                    :moisture_evaporation,
                    7 => 13,
                    A = 2.997e5,
                    E = 5.0e4,
                    h = hos,
                    n = 1.0,
                ),
            ),
            mixing = MaterialMixing(
                k = MixingSpec(WEIGHTED, 0.5),
                λ = MixingSpec(WEIGHTED, 0.5),
                κ = MixingSpec(WEIGHTED, 0.5),
            ),
        )

        # Problem setup
        T_ambient = 300.0
        ρ_DF_initial = 530.0
        ρ_moisture_initial = 0.06 * ρ_DF_initial

        mesh = create_uniform_mesh(
            12.0e-3,
            240,
            Val(13);
            T_initial = T_ambient,
            ξ_initial = (
                ρ_DF_initial,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                ρ_moisture_initial,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
        )

        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50000.0) +
                      RadiativeBC(T_inf = 300) +
                      ConvectiveBC(h = 10.0, T_inf = 300),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
                bottom = ImpermeableBC(),
            ),
        )

        problem = PyrolysisProblem(
            mesh = mesh,
            material = material,
            bc_set = bc_set,
            T_initial = z -> T_ambient,
            ξ_initial = [
                z -> ρ_DF_initial,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> ρ_moisture_initial,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
                z -> 0.0,
            ],
            tspan = (0.0, 600.0),
            dt_initial = 1e-3,
            output_times = collect(range(0.0, 600.0, length = 3601)),
        )

        solution = solve(
            problem;
            abstol = 1e-9,
            reltol = 1e-7,
            dtmin = 1e-14,
            progress = false,
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            min_thickness = 0.5e-3,
            handle_depletion = false,
            # Default Jacobian is `Structured(strategy = DirectSolve())`.
            integrator = KenCarp4(linsolve = SparspakFactorization()),
        )

        return solution
    end

    # ========================================================================
    # QoI extraction
    # ========================================================================

    const LC_MLR_IGNITION_THRESHOLD = 4e-3   # kg/(m²·s) — matches sensitivity analysis

    # --- New QoI constants ---
    const LC_FIRST_PEAK_T_START = 30.0       # [s] skip initial moisture transient
    const LC_PYROLYSIS_ONSET_THRESHOLD = 0.01 # 1% dry mass loss
    const LC_SOLID_COMPONENT_INDICES = 1:6    # DF, int1-4, Char
    # Combustible gas indices (DF_gas1-5, components 8-12) — excludes water vapor (13)
    const LC_COMBUSTIBLE_GAS_INDICES = 8:12
    # Critical combustible MLR threshold for ignition/extinction [kg/(m²·s)]
    const LC_COMBUSTIBLE_MLR_THRESHOLD = 1.5e-3

    """
    		time_to_threshold(t, curve, fraction)

    	Find the time at which `curve` first reaches `fraction` of its total change,
    	defined as T₀ + fraction*(T_max - T₀) where T₀ = curve[1].
    	Uses linear interpolation between time steps. Returns Inf if never reached.
    	"""
    function time_to_threshold(t, curve, fraction)
        T0 = curve[1]
        T_max = maximum(curve)
        threshold = T0 + fraction * (T_max - T0)
        idx = findfirst(x -> x >= threshold, curve)
        if idx === nothing
            return Inf
        elseif idx == 1
            return t[1]
        else
            f = (threshold - curve[idx-1]) / (curve[idx] - curve[idx-1])
            return t[idx-1] + f * (t[idx] - t[idx-1])
        end
    end

    """
    		parabolic_peak(t, y, i) -> (peak_value, peak_time)

    	Fit a parabola through (i-1, i, i+1) and return the interpolated peak.
    	Falls back to the grid point if the parabola opens upward or the vertex
    	falls outside the interval.
    	"""
    function parabolic_peak(t, y, i)
        t0, t1, t2 = t[i-1], t[i], t[i+1]
        y0, y1, y2 = y[i-1], y[i], y[i+1]
        denom = y0 - 2y1 + y2
        if denom >= 0
            return (y1, t1)
        end
        h = t1 - t0
        dt = 0.5 * h * (y0 - y2) / denom
        t_star = t1 + dt
        if t_star < t0 || t_star > t2
            return (y1, t1)
        end
        frac = dt / h
        y_star = y1 - 0.25 * (y0 - y2) * frac
        return (y_star, t_star)
    end

    function compute_dry_mass_series(sol)
        nt = length(sol.t)
        A = sol.problem.mesh.cross_section_area
        dry_mass = Vector{Float64}(undef, nt)
        for i = 1:nt
            ξ = sol.ξ[:, :, i]
            n_cells = size(ξ, 2)
            m = 0.0
            for j = 1:n_cells
                v = sol.problem.mesh.cell_volumes[j]
                for k in LC_SOLID_COMPONENT_INDICES
                    m += ξ[k, j] * v
                end
            end
            dry_mass[i] = m / A
        end
        return dry_mass
    end

    function find_first_pyrolysis_peak(t, mlr, global_peak_idx)
        i_start = findfirst(x -> x >= LC_FIRST_PEAK_T_START, t)
        if i_start === nothing || i_start >= global_peak_idx
            return (NaN, NaN)
        end
        for i = (i_start+1):(global_peak_idx-1)
            if mlr[i] >= mlr[i-1] && mlr[i] >= mlr[i+1]
                return parabolic_peak(t, mlr, i)
            end
        end
        i_mid = div(i_start + global_peak_idx, 2)
        if i_mid > i_start
            local_range = i_start:i_mid
            idx = local_range[argmax(mlr[local_range])]
            if idx > 1 && idx < length(t)
                return parabolic_peak(t, mlr, idx)
            end
            return (mlr[idx], t[idx])
        end
        return (NaN, NaN)
    end

    function find_pyrolysis_onset(t, dry_mass)
        m0 = dry_mass[1]
        threshold_mass = m0 * (1.0 - LC_PYROLYSIS_ONSET_THRESHOLD)
        idx = findfirst(m -> m <= threshold_mass, dry_mass)
        if idx === nothing
            return Inf
        elseif idx == 1
            return t[1]
        else
            f = (dry_mass[idx-1] - threshold_mass) / (dry_mass[idx-1] - dry_mass[idx])
            return t[idx-1] + f * (t[idx] - t[idx-1])
        end
    end

    function find_back_temp_at_mass_loss(t, total_mass, back_temp, fraction)
        m0 = total_mass[1]
        threshold = m0 * (1.0 - fraction)
        idx = findfirst(m -> m <= threshold, total_mass)
        if idx === nothing
            return NaN
        elseif idx == 1
            return back_temp[1]
        else
            f = (total_mass[idx-1] - threshold) / (total_mass[idx-1] - total_mass[idx])
            return back_temp[idx-1] + f * (back_temp[idx] - back_temp[idx-1])
        end
    end

    """
    		compute_combustible_mlr_series(sol) -> Vector{Float64}

    	Compute the combustible-only surface mass flux [kg/(m²·s)] at each output
    	time. Mirrors the solver's per-species flux computation but sums only over
    	combustible gas components (indices 8-12), excluding water vapor (13).
    	"""
    function compute_combustible_mlr_series(sol)
        nt = length(sol.t)
        comb_mlr = zeros(nt)

        mesh = sol.problem.mesh
        material = sol.problem.material
        bc_set = sol.problem.bc_set
        mass_bc_top = get_mass_bc(bc_set, Val(:top))
        n_cells = length(mesh.cells)

        d_surface = mesh.cell_volumes[n_cells] / mesh.cross_section_area / 2

        for i = 1:nt
            t = sol.t[i]
            ξ_profile = sol.ξ[:, :, i]
            T_surface = sol.surface_T[i]

            ξ_surface = ntuple(k -> ξ_profile[k, n_cells], Val(length(material.components)))

            flux = 0.0
            for k in LC_COMBUSTIBLE_GAS_INDICES
                ξ_k = ξ_surface[k]
                if ξ_k > 0.0
                    λρ =
                        effective_gas_transfer(ξ_surface, T_surface, material) *
                        sum(ξ_surface)
                    ρ_surface = sum(ξ_surface)
                    J_k = compute_surface_mass_flux(
                        mass_bc_top,
                        ξ_k,
                        λρ,
                        ρ_surface,
                        d_surface,
                        t,
                        k,
                    )
                    flux += J_k
                end
            end
            comb_mlr[i] = flux
        end

        return comb_mlr
    end

    """
    		find_threshold_crossing(t, y, threshold, direction) -> Float64

    	Find time at which `y` crosses `threshold`. `direction` is `:rising` or
    	`:falling`. Uses linear interpolation. Returns `Inf` if no crossing found.
    	"""
    function find_threshold_crossing(t, y, threshold, direction::Symbol)
        for i = 2:length(t)
            crossed = if direction == :rising
                y[i-1] <= threshold && y[i] > threshold
            else
                y[i-1] >= threshold && y[i] < threshold
            end
            if crossed
                f = (threshold - y[i-1]) / (y[i] - y[i-1])
                return t[i-1] + f * (t[i] - t[i-1])
            end
        end
        return Inf
    end

    function extract_qois(sol)
        t = sol.t
        mlr = sol.extras.surface_mass_flux
        T_surf = sol.surface_T
        T_back = sol.extras.back_temperature
        total_mass = sol.extras.total_mass

        peak_idx = argmax(mlr)
        peak_mlr_val, time_to_peak = parabolic_peak(t, mlr, peak_idx)

        tig_idx = findfirst(x -> x > LC_MLR_IGNITION_THRESHOLD, mlr)
        if tig_idx === nothing
            time_to_ignition = Inf
        elseif tig_idx == 1
            time_to_ignition = t[1]
        else
            f =
                (LC_MLR_IGNITION_THRESHOLD - mlr[tig_idx-1]) /
                (mlr[tig_idx] - mlr[tig_idx-1])
            time_to_ignition = t[tig_idx-1] + f * (t[tig_idx] - t[tig_idx-1])
        end

        total_mass_loss = (1.0 - total_mass[end] / total_mass[1]) * 100.0

        # Threshold-crossing times for temperature curves
        t_surf_70 = time_to_threshold(t, T_surf, 0.70)
        t_surf_95 = time_to_threshold(t, T_surf, 0.95)
        t_back_5 = time_to_threshold(t, T_back, 0.05)
        t_back_95 = time_to_threshold(t, T_back, 0.95)

        # --- New QoIs ---
        # 1. First pyrolysis peak MLR (skip moisture transient)
        peak_first, time_peak_first = find_first_pyrolysis_peak(t, mlr, peak_idx)

        # 2. Pyrolysis onset (1% dry mass loss)
        dry_mass = compute_dry_mass_series(sol)
        onset_time = find_pyrolysis_onset(t, dry_mass)

        # 3. Back temperature at 50% total mass loss
        T_back_50 = find_back_temp_at_mass_loss(t, total_mass, T_back, 0.50)

        # 4. Combustible-only MLR: ignition, extinction, and average burning rate
        comb_mlr = compute_combustible_mlr_series(sol)
        t_comb_ign =
            find_threshold_crossing(t, comb_mlr, LC_COMBUSTIBLE_MLR_THRESHOLD, :rising)
        t_comb_ext = if isfinite(t_comb_ign)
            ign_idx = findfirst(ti -> ti >= t_comb_ign, t)
            find_threshold_crossing(
                t[ign_idx:end],
                comb_mlr[ign_idx:end],
                LC_COMBUSTIBLE_MLR_THRESHOLD,
                :falling,
            )
        else
            Inf
        end

        avg_comb_mlr =
            if isfinite(t_comb_ign) && isfinite(t_comb_ext) && t_comb_ext > t_comb_ign
                burn_mask = t_comb_ign .<= t .<= t_comb_ext
                burn_idx = findall(burn_mask)
                if length(burn_idx) >= 2
                    total_integral = 0.0
                    for j = 2:length(burn_idx)
                        i1, i2 = burn_idx[j-1], burn_idx[j]
                        total_integral +=
                            0.5 * (comb_mlr[i1] + comb_mlr[i2]) * (t[i2] - t[i1])
                    end
                    total_integral / (t_comb_ext - t_comb_ign)
                else
                    NaN
                end
            else
                NaN
            end

        # 5. Peak surface temperature
        peak_T_surf = maximum(T_surf)

        return (
            peak_mlr = peak_mlr_val,
            time_to_peak_mlr = time_to_peak,
            time_to_ignition = time_to_ignition,
            total_mass_loss_pct = total_mass_loss,
            t_surf_70pct = t_surf_70,
            t_surf_95pct = t_surf_95,
            t_back_5pct = t_back_5,
            t_back_95pct = t_back_95,
            final_thickness = sol.extras.thickness_history[end],
            peak_mlr_first = peak_first,
            time_to_peak_mlr_first = time_peak_first,
            pyrolysis_onset_time = onset_time,
            T_back_50pct_mass_loss = T_back_50,
            time_to_comb_ignition = t_comb_ign,
            time_to_comb_extinction = t_comb_ext,
            avg_mlr_comb_burning = avg_comb_mlr,
            peak_surface_temp = peak_T_surf,
            dry_mass_curve = dry_mass,
            combustible_mlr_curve = comb_mlr,
        )
    end

    # ========================================================================
    # Worker entry point: solve one run, save CSV, return QoIs
    # ========================================================================

    function run_single(spec::LitCompRunSpec, results_dir::String)
        mkpath(results_dir)
        t_start = time()
        sol = build_and_solve(spec)
        elapsed = time() - t_start

        qois = extract_qois(sol)

        # Save time-series CSV (including dry mass and combustible MLR for
        # new QoI reproducibility)
        filepath = joinpath(results_dir, "$(spec.label)_results.csv")
        open(filepath, "w") do io
            println(
                io,
                "time,mlr,surface_temp,back_temp,total_mass,thickness,dry_mass,combustible_mlr",
            )
            for i in eachindex(sol.t)
                @printf(
                    io,
                    "%.6g,%.10g,%.6g,%.6g,%.10g,%.10g,%.10g,%.10g\n",
                    sol.t[i],
                    sol.extras.surface_mass_flux[i],
                    sol.surface_T[i],
                    sol.extras.back_temperature[i],
                    sol.extras.total_mass[i],
                    sol.extras.thickness_history[i],
                    qois.dry_mass_curve[i],
                    qois.combustible_mlr_curve[i]
                )
            end
        end

        return (label = spec.label, qois = qois, elapsed = elapsed, retcode = sol.retcode)
    end

end # @everywhere

# ============================================================================
# Run definitions
# ============================================================================

const RESULTS_DIR = joinpath(@__DIR__, "results")

const RUNS = [
    LitCompRunSpec("baseline", :baseline, :baseline, 1.0, NaN, :baseline),
    LitCompRunSpec("virgin_cp_harada", :harada, :baseline, 1.0, NaN, :baseline),
    LitCompRunSpec("virgin_cp_dupont", :dupont, :baseline, 1.0, NaN, :baseline),
    LitCompRunSpec("char_cp_hostikka", :baseline, :hostikka, 1.0, NaN, :baseline),
    LitCompRunSpec("char_cp_edrich", :baseline, :edrich, 1.0, NaN, :baseline),
    # Rath uniform removed — two-zone (run 9) is the rigorous Rath comparison
    LitCompRunSpec("hrxn_jones", :baseline, :baseline, 306.92 / 115.627, NaN, :baseline),     # Jones with-lid linear fit at Y_c=0.17680 (our DF at 5 K/min, 823 K) → +306.92 kJ/kg
    LitCompRunSpec("hos_free_water", :baseline, :baseline, 1.0, NaN, :free_water),
    LitCompRunSpec("hrxn_rp_spruce", :baseline, :baseline, 112.0 / 115.627, NaN, :baseline),  # Rinta-Paavola uniform (+112 kJ/kg)
    LitCompRunSpec(
        "hrxn_rath_2zone",
        :baseline,
        :baseline,
        218.36 / 177.0,
        38.9 / 61.0,
        :baseline,
    ),  # Rath two-zone (f_P=1.234, f_S=0.638)
]

# ============================================================================
# Main
# ============================================================================

function main()
    mkpath(RESULTS_DIR)

    println("\n", "="^60)
    println("  Literature Property Comparison — Douglas Fir")
    println(
        "  $(length(RUNS)) simulations (1 baseline + $(length(RUNS)-1) variants) on $(nworkers()) worker(s)",
    )
    println("="^60, "\n")

    # Run all simulations in parallel via pmap
    results = pmap(spec -> run_single(spec, RESULTS_DIR), RUNS)

    # Print results and save summary
    for r in results
        @printf("%-22s  %.1f s  (retcode: %s)\n", r.label, r.elapsed, r.retcode)
        q = r.qois
        @printf("  Peak MLR: %.4e kg/(m²·s) at t=%.1f s\n", q.peak_mlr, q.time_to_peak_mlr)
        @printf(
            "  1st peak MLR: %.4e kg/(m²·s) at t=%.1f s\n",
            q.peak_mlr_first,
            q.time_to_peak_mlr_first
        )
        @printf("  Time to ignition (4 g/m²·s): %.1f s\n", q.time_to_ignition)
        @printf("  Total mass loss: %.2f%%\n", q.total_mass_loss_pct)
        @printf("  Pyrolysis onset (1%% dry): %.1f s\n", q.pyrolysis_onset_time)
        @printf("  T_back at 50%% mass loss: %.1f K\n", q.T_back_50pct_mass_loss)
        @printf(
            "  t_surf 70%%: %.1f s, t_surf 95%%: %.1f s\n",
            q.t_surf_70pct,
            q.t_surf_95pct
        )
        @printf(
            "  t_back 5%%: %.1f s, t_back 95%%: %.1f s\n",
            q.t_back_5pct,
            q.t_back_95pct
        )
        @printf(
            "  Combustible ignition: %.1f s, extinction: %.1f s\n",
            q.time_to_comb_ignition,
            q.time_to_comb_extinction
        )
        @printf("  Avg combustible MLR (burning): %.4e kg/(m²·s)\n", q.avg_mlr_comb_burning)
        @printf("  Peak surface temp: %.1f K\n", q.peak_surface_temp)
        @printf("  Final thickness: %.4e m\n\n", q.final_thickness)
    end

    # Save summary QoI table
    qoi_file = joinpath(RESULTS_DIR, "literature_comparison_qois.csv")
    open(qoi_file, "w") do io
        println(
            io,
            "run_id,label,peak_mlr,time_to_peak_mlr,time_to_ignition," *
            "total_mass_loss_pct,t_surf_70pct,t_surf_95pct,t_back_5pct,t_back_95pct," *
            "final_thickness,peak_mlr_first,time_to_peak_mlr_first," *
            "pyrolysis_onset_time,T_back_50pct_mass_loss," *
            "time_to_comb_ignition,time_to_comb_extinction,avg_mlr_comb_burning," *
            "peak_surface_temp",
        )
        for (i, r) in enumerate(results)
            q = r.qois
            @printf(
                io,
                "%d,%s,%.10g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.10g,%.10g,%.6g,%.6g,%.6g,%.6g,%.6g,%.10g,%.6g\n",
                i,
                r.label,
                q.peak_mlr,
                q.time_to_peak_mlr,
                q.time_to_ignition,
                q.total_mass_loss_pct,
                q.t_surf_70pct,
                q.t_surf_95pct,
                q.t_back_5pct,
                q.t_back_95pct,
                q.final_thickness,
                q.peak_mlr_first,
                q.time_to_peak_mlr_first,
                q.pyrolysis_onset_time,
                q.T_back_50pct_mass_loss,
                q.time_to_comb_ignition,
                q.time_to_comb_extinction,
                q.avg_mlr_comb_burning,
                q.peak_surface_temp
            )
        end
    end

    println("="^60)
    println("  All runs complete. Results saved to: $RESULTS_DIR")
    println("="^60)
end

main()
