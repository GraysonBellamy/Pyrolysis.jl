"""
Literature Property Comparison — Douglas Fir Cone Calorimeter

Runs 9 simulations with literature-sourced property values substituted one
category at a time, using the sensitivity analysis baseline as the reference
configuration (240 cells, 12 mm slab, 50 kW/m², 6% MC, 600 s).

Variants:
  1. Virgin cp — Harada et al. (1998), low
  2. Virgin cp — Dupont et al. (2014), high
  3. Char cp — Hostikka & Matala (2017), low
  4. Char cp — Edrich (1980), high
  5. Heats of reaction — Rath et al. (2003), ΔH_P = +428 kJ/kg (uniform)
  6. Heats of reaction — Jones et al. (2020), ΔH_P = +382 kJ/kg (uniform)
  7. Heat of sorption — Free water only (2442 kJ/kg)
  8. Heats of reaction — Jones et al. (2020), two-zone (f_P=2.244, f_S=0.279)
  9. Heats of reaction — Rath et al. (2003), two-zone (f_P=2.468, f_S=0.185)

Outputs (saved to results/):
  - {label}_results.csv per run (time, mlr, surface_temp, back_temp, total_mass, thickness)
  - literature_comparison_qois.csv summary table

Usage:
  julia -p 9 --project=. examples/literature_comparison/literature_comparison.jl
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
	using OrdinaryDiffEq
	using LinearSolve
	using KLU
	using LinearSolve: KLUFactorization
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
	const LC_NOMINAL_HRXN = [-19153.3, -5573.79, 266699.7, -219638.7, -12381.6]

	# ========================================================================
	# Baseline property functions
	# ========================================================================

	function baseline_virgin_cp(T)
		if T <= 323.0
			-53.052 + 4.2453 * T
		elseif T <= 425.0
			cp_low  = -53.052 + 4.2453 * 323.0
			cp_high = -2131.33 + 9.11878 * 425.0
			cp_low + (cp_high - cp_low) * (T - 323.0) / 102.0
		else
			-2131.33 + 9.11878 * T
		end
	end

	baseline_char_cp(T) = 347.212 + 1.78423 * T

	function baseline_hos(T, ξ)
		MC = 100.0 * ξ[LC_MOISTURE_IDX] / LC_ρ_DRY_BASIS
		return (2450.4209 + 875.6441 * exp(-0.11217 * MC)) * 1e3
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

	const VIRGIN_CP_DISPATCH = Dict{Symbol, Function}(
		:baseline => baseline_virgin_cp,
		:harada   => harada_virgin_cp,
		:dupont   => dupont_virgin_cp,
	)
	const CHAR_CP_DISPATCH = Dict{Symbol, Function}(
		:baseline => baseline_char_cp,
		:hostikka => hostikka_char_cp,
		:edrich   => edrich_char_cp,
	)
	const HOS_DISPATCH = Dict{Symbol, Function}(
		:baseline   => baseline_hos,
		:free_water => free_water_hos,
	)

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
		char_cp   = CHAR_CP_DISPATCH[spec.char_cp_key]
		hos_func  = HOS_DISPATCH[spec.hos_key]

		# Two-zone scaling: primary (rxns 1-3) vs secondary (rxns 4-5)
		f_p = spec.hrxn_scale
		f_s = isnan(spec.hrxn_scale_secondary) ? spec.hrxn_scale : spec.hrxn_scale_secondary
		hrxn = [LC_NOMINAL_HRXN[1] * f_p, LC_NOMINAL_HRXN[2] * f_p, LC_NOMINAL_HRXN[3] * f_p,
			LC_NOMINAL_HRXN[4] * f_s, LC_NOMINAL_HRXN[5] * f_s]

		# Intermediate cp's via cumulative-yield interpolation
		int1_cp = CallableProperty(T -> (1 - LC_ALPHA[1]) * virgin_cp(T) + LC_ALPHA[1] * char_cp(T))
		int2_cp = CallableProperty(T -> (1 - LC_ALPHA[2]) * virgin_cp(T) + LC_ALPHA[2] * char_cp(T))
		int3_cp = CallableProperty(T -> (1 - LC_ALPHA[3]) * virgin_cp(T) + LC_ALPHA[3] * char_cp(T))
		int4_cp = CallableProperty(T -> (1 - LC_ALPHA[4]) * virgin_cp(T) + LC_ALPHA[4] * char_cp(T))

		df_cp        = CallableProperty(virgin_cp)
		char_cp_prop = CallableProperty(char_cp)

		# Moisture cp (NIST Shomate, fixed)
		moisture_cp = CallableProperty(T -> begin
			t = T / 1000.0
			M = 0.01801528
			A = -203.6060;
			B = 1523.290;
			C = -3196.413;
			D = 2474.455;
			E = 3.855326
			(A + B*t + C*t^2 + D*t^3 + E/t^2) / M
		end)

		# Water vapor cp (fixed)
		water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)

		# Heat of sorption
		hos = StateDependentProperty(hos_func)

		# Gas transfer / permeability constants
		λ_wood   = 1.0e-6
		λ_char   = 1.0e-5
		λ_gas    = 1.0e-5
		κ_virgin = 2.3e-15
		κ_char   = 1.0e-14

		material = Material(
			name = :DouglasFir,
			components = (
				SolidComponent(:DF,
					ρ = 530.0, c = df_cp, k = 0.145,
					ε = 0.81, α = 0.0, λ = λ_wood, κ = κ_virgin),
				SolidComponent(:DF_int1,
					ρ = 507.2, c = int1_cp, k = 0.292951,
					ε = 0.790965, α = 0.0, λ = λ_wood * 1.5, κ = 2.7e-15),
				SolidComponent(:DF_int2,
					ρ = 419.8, c = int2_cp, k = 0.329938,
					ε = 0.786207, α = 0.0, λ = λ_wood * 3.0, κ = 3.1e-15),
				SolidComponent(:DF_int3,
					ρ = 203.2, c = int3_cp, k = 0.644833,
					ε = 0.745694, α = 0.0, λ = λ_wood * 6.0, κ = 3.6e-15),
				SolidComponent(:DF_int4,
					ρ = 176.6, c = int4_cp, k = 0.777700,
					ε = 0.728600, α = 0.0, λ = λ_char * 0.8, κ = 8.3e-15),
				SolidComponent(:Char,
					ρ = 150.0, c = char_cp_prop, k = 1.0,
					ε = 0.70, α = 0.0, λ = λ_char, κ = κ_char),
				LiquidComponent(:MOISTURE,
					ρ = 1e6, c = moisture_cp, k = 0.145,
					ε = 0.81, α = 0.0, κ = 0.0),
				GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas * 0.85),
				GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas * 1.0),
				GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.6),
				GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas * 1.3),
				GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas * 1.4),
				GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.6),
			),
			reactions = (
				Reaction(:pyrolysis_step1, 1 => (2, 8),
					A = 74816.95, E = 71.59e3, h = hrxn[1], n = 1.0, yields = (0.95, 0.05)),
				Reaction(:pyrolysis_step2, 2 => (3, 9),
					A = 23550492839.0, E = 140.9e3, h = hrxn[2], n = 1.0, yields = (0.80, 0.20)),
				Reaction(:pyrolysis_step3, 3 => (4, 10),
					A = 7.744618e12, E = 181.0e3, h = hrxn[3], n = 1.0, yields = (0.37, 0.63)),
				Reaction(:pyrolysis_step4, 4 => (5, 11),
					A = 103276.14, E = 98.94e3, h = hrxn[4], n = 1.0, yields = (0.79, 0.21)),
				Reaction(:pyrolysis_step5, 5 => (6, 12),
					A = 6.397348, E = 56.01e3, h = hrxn[5], n = 1.0, yields = (0.74, 0.26)),
				Reaction(:moisture_evaporation, 7 => 13,
					A = 2.997e5, E = 5.0e4, h = hos, n = 1.0),
			),
			conductivity_mixing = WEIGHTED, β = 0.5,
			gas_transfer_mixing = WEIGHTED, β_λ = 0.5,
			permeability_mixing = WEIGHTED, β_κ = 0.5,
		)

		# Problem setup
		T_ambient = 300.0
		ρ_DF_initial = 530.0
		ρ_moisture_initial = 0.06 * ρ_DF_initial

		mesh = create_uniform_mesh(
			12.0e-3, 240, Val(13);
			T_initial = T_ambient,
			ξ_initial = (ρ_DF_initial, 0.0, 0.0, 0.0, 0.0, 0.0, ρ_moisture_initial,
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
		)

		bc_set = BoundaryConditionSet(
			thermal_top = RadiativeFluxBC(flux = 50e3) +
						  RadiativeBC(T_inf = 300) +
						  ConvectiveBC(h = 10.0, T_inf = 300),
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
			mass_bottom = ImpermeableBC(),
		)

		problem = PyrolysisProblem(
			mesh = mesh,
			material = material,
			bc_set = bc_set,
			T_initial = z -> T_ambient,
			ξ_initial = [
				z -> ρ_DF_initial,
				z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
				z -> ρ_moisture_initial,
				z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
			],
			external_flux = t -> 0.0,
			ambient_temperature = t -> T_ambient,
			tspan = (0.0, 600.0),
			dt_initial = 1e-3,
			output_times = collect(range(0.0, 600.0, length = 3601)),
		)

		solution = solve(problem;
			abstol = 1e-9,
			reltol = 1e-7,
			dtmin = 1e-14,
			progress = false,
			use_ale = true,
			radiation_model = SURFACE_ABSORPTION,
			min_thickness = 0.5e-3,
			handle_depletion = false,
			use_analytical_jacobian = true,
			solver = KenCarp4(linsolve = KLUFactorization()),
		)

		return solution
	end

	# ========================================================================
	# QoI extraction
	# ========================================================================

	const LC_MLR_IGNITION_THRESHOLD = 3e-3

	# --- New QoI constants ---
	const LC_FIRST_PEAK_T_START = 30.0       # [s] skip initial moisture transient
	const LC_PYROLYSIS_ONSET_THRESHOLD = 0.01 # 1% dry mass loss
	const LC_SOLID_COMPONENT_INDICES = 1:6    # DF, int1-4, Char

	function time_to_threshold(t, curve, fraction)
		T0 = curve[1]
		T_max = maximum(curve)
		threshold = T0 + fraction * (T_max - T0)
		idx = findfirst(x -> x >= threshold, curve)
		return idx !== nothing ? t[idx] : Inf
	end

	function compute_dry_mass_series(sol)
		nt = length(sol.times)
		A = sol.problem.mesh.cross_section_area
		dry_mass = Vector{Float64}(undef, nt)
		for i in 1:nt
			ξ = sol.concentrations[i]
			n_cells = size(ξ, 2)
			m = 0.0
			for j in 1:n_cells
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
		for i in (i_start+1):(global_peak_idx-1)
			if mlr[i] >= mlr[i-1] && mlr[i] >= mlr[i+1]
				return (mlr[i], t[i])
			end
		end
		i_mid = div(i_start + global_peak_idx, 2)
		if i_mid > i_start
			local_range = i_start:i_mid
			idx = local_range[argmax(mlr[local_range])]
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

	function extract_qois(sol)
		t = sol.times
		mlr = sol.surface_mass_flux
		T_surf = sol.surface_temperature
		T_back = sol.back_temperature
		total_mass = sol.total_mass

		peak_mlr = maximum(mlr)
		peak_idx = argmax(mlr)
		time_to_peak = t[peak_idx]

		tig_idx = findfirst(x -> x > LC_MLR_IGNITION_THRESHOLD, mlr)
		time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf

		total_mass_loss = (1.0 - total_mass[end] / total_mass[1]) * 100.0

		# Threshold-crossing times for temperature curves
		t_surf_70 = time_to_threshold(t, T_surf, 0.70)
		t_surf_95 = time_to_threshold(t, T_surf, 0.95)
		t_back_5  = time_to_threshold(t, T_back, 0.05)
		t_back_95 = time_to_threshold(t, T_back, 0.95)

		# --- New QoIs ---
		peak_first, time_peak_first = find_first_pyrolysis_peak(t, mlr, peak_idx)
		dry_mass = compute_dry_mass_series(sol)
		onset_time = find_pyrolysis_onset(t, dry_mass)
		T_back_50 = find_back_temp_at_mass_loss(t, total_mass, T_back, 0.50)

		return (
			peak_mlr = peak_mlr,
			time_to_peak_mlr = time_to_peak,
			time_to_ignition = time_to_ignition,
			total_mass_loss_pct = total_mass_loss,
			t_surf_70pct = t_surf_70,
			t_surf_95pct = t_surf_95,
			t_back_5pct = t_back_5,
			t_back_95pct = t_back_95,
			final_thickness = sol.thickness_history[end],
			peak_mlr_first = peak_first,
			time_to_peak_mlr_first = time_peak_first,
			pyrolysis_onset_time = onset_time,
			T_back_50pct_mass_loss = T_back_50,
			dry_mass_curve = dry_mass,
		)
	end

	# ========================================================================
	# Worker entry point: solve one run, save CSV, return QoIs
	# ========================================================================

	function run_single(spec::LitCompRunSpec, results_dir::String)
		t_start = time()
		sol = build_and_solve(spec)
		elapsed = time() - t_start

		qois = extract_qois(sol)

		# Save time-series CSV (including dry mass for new QoI reproducibility)
		filepath = joinpath(results_dir, "$(spec.label)_results.csv")
		open(filepath, "w") do io
			println(io, "time,mlr,surface_temp,back_temp,total_mass,thickness,dry_mass")
			for i in eachindex(sol.times)
				@printf(io, "%.6g,%.10g,%.6g,%.6g,%.10g,%.10g,%.10g\n",
					sol.times[i],
					sol.surface_mass_flux[i],
					sol.surface_temperature[i],
					sol.back_temperature[i],
					sol.total_mass[i],
					sol.thickness_history[i],
					qois.dry_mass_curve[i])
			end
		end

		return (label = spec.label, qois = qois, elapsed = elapsed,
			retcode = sol.retcode)
	end

end # @everywhere

# ============================================================================
# Run definitions
# ============================================================================

const RESULTS_DIR = joinpath(@__DIR__, "results")

const RUNS = [
	LitCompRunSpec("virgin_cp_harada", :harada, :baseline, 1.0, NaN, :baseline),
	LitCompRunSpec("virgin_cp_dupont", :dupont, :baseline, 1.0, NaN, :baseline),
	LitCompRunSpec("char_cp_hostikka", :baseline, :hostikka, 1.0, NaN, :baseline),
	LitCompRunSpec("char_cp_edrich", :baseline, :edrich, 1.0, NaN, :baseline),
	LitCompRunSpec("hrxn_rath", :baseline, :baseline, 428.0 / 113.74, NaN, :baseline),
	LitCompRunSpec("hrxn_jones", :baseline, :baseline, 382.0 / 113.74, NaN, :baseline),
	LitCompRunSpec("hos_free_water", :baseline, :baseline, 1.0, NaN, :free_water),
	LitCompRunSpec("hrxn_jones_2zone", :baseline, :baseline, 400.0 / 178.25, -18.0 / -64.51, :baseline),
	LitCompRunSpec("hrxn_rath_2zone", :baseline, :baseline, 440.0 / 178.25, -12.0 / -64.51, :baseline),
]

# ============================================================================
# Main
# ============================================================================

function main()
	mkpath(RESULTS_DIR)

	println("\n", "="^60)
	println("  Literature Property Comparison — Douglas Fir")
	println("  $(length(RUNS)) variant simulations (incl. 2 two-zone) on $(nworkers()) worker(s)")
	println("="^60, "\n")

	# Run all simulations in parallel via pmap
	results = pmap(spec -> run_single(spec, RESULTS_DIR), RUNS)

	# Print results and save summary
	for r in results
		@printf("%-22s  %.1f s  (retcode: %s)\n", r.label, r.elapsed, r.retcode)
		q = r.qois
		@printf("  Peak MLR: %.4e kg/(m²·s) at t=%.1f s\n", q.peak_mlr, q.time_to_peak_mlr)
		@printf("  1st peak MLR: %.4e kg/(m²·s) at t=%.1f s\n", q.peak_mlr_first, q.time_to_peak_mlr_first)
		@printf("  Total mass loss: %.2f%%\n", q.total_mass_loss_pct)
		@printf("  Pyrolysis onset (1%% dry): %.1f s\n", q.pyrolysis_onset_time)
		@printf("  T_back at 50%% mass loss: %.1f K\n", q.T_back_50pct_mass_loss)
		@printf("  t_surf 70%%: %.1f s, t_surf 95%%: %.1f s\n", q.t_surf_70pct, q.t_surf_95pct)
		@printf("  t_back 5%%: %.1f s, t_back 95%%: %.1f s\n", q.t_back_5pct, q.t_back_95pct)
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
			"pyrolysis_onset_time,T_back_50pct_mass_loss",
		)
		for (i, r) in enumerate(results)
			q = r.qois
			@printf(io, "%d,%s,%.10g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.10g,%.10g,%.6g,%.6g,%.6g\n",
				i, r.label,
				q.peak_mlr, q.time_to_peak_mlr, q.time_to_ignition,
				q.total_mass_loss_pct, q.t_surf_70pct, q.t_surf_95pct,
				q.t_back_5pct, q.t_back_95pct,
				q.final_thickness,
				q.peak_mlr_first, q.time_to_peak_mlr_first,
				q.pyrolysis_onset_time, q.T_back_50pct_mass_loss)
		end
	end

	println("="^60)
	println("  All runs complete. Results saved to: $RESULTS_DIR")
	println("="^60)
end

main()
