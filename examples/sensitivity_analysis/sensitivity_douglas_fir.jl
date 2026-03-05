"""
Sensitivity Analysis of Douglas Fir Cone Calorimeter Simulation

Propagates uncertainties in thermophysical input parameters through the
07_douglas_fir_cone.jl simulation using Sobol quasi-random sampling.

Uncertain parameters (17 total):
  - Systematic (Type B) instrument biases: ε_HFM, ε_STA, ε_water (3)
	  Multiplicative calibration biases that scale measurements from each
	  instrument uniformly. Following GUM (JCGM 100:2008):
	  ε_HFM ~ N(0, 0.00577²)  — HFM heat flow accuracy (1%/√3, manufacturer spec)
	  ε_STA ~ N(0, 0.0261²)   — DSC multi-point enthalpy calibration (7 standards)
	  ε_water ~ N(0, 0.01732²) — isothermal water-evaporation calibration (22 replicates)
	  ε_STA applies to DSC-derived Cp (virgin seg 3, char) and ΔH.
	  ε_water applies to HoS (independent calibration pathway, not ε_STA).
	  ε_HFM applies to HFM-derived Cp (virgin seg 1).
  - Heat capacity regression coefficients for DF virgin (4) and Char (2)
	(Segment 2 [323–425 K] is interpolated from segments 1 & 3, not independent)
	(Intermediate cp's are derived via cumulative-yield interpolation)
  - Heat of sorption equation parameters (3)
  - Heats of reaction for pyrolysis steps 1-5 (5)

Covariance matrices from: Bellamy, McKinnon, Stoliarov, "Thermochemical Properties
of Wood and Sorbed Moisture for Modeling of Wood Pyrolysis" (energetics.tex).
Type A parameters within each group are sampled jointly via MvNormal to preserve
cross-correlations (e.g., intercept–slope anti-correlation in cp regressions,
amplitude–decay correlation in the HoS fit, shared baseline uncertainty in ΔH).
Type B parameters are sampled independently per instrument (STA, HFM) and applied
as multiplicative factors (1 + ε) to all Cp values from that instrument.

Kinetic parameters (A, E) and all other properties are held fixed.

Outputs:
  - Scalar QoIs: peak MLR, time to ignition, total mass loss, peak temperatures
  - Time-series: MLR curve, surface temperature, back temperature per run
  - Sensitivity ranking via Standardized Rank Regression Coefficients (SRRC)
  - Uncertainty envelope plots (median, 50% CI, 90% CI)

Usage:
  julia -p 8 sensitivity_douglas_fir.jl          # 8 parallel workers
  julia sensitivity_douglas_fir.jl               # serial execution
"""

# ============================================================================
# Setup: Load packages on all workers
# ============================================================================

using Distributed

# Add workers if launched with -p flag (already added by Julia)
# If running interactively, uncomment:
# nprocs() == 1 && addprocs(8)

println("Running with $(nworkers()) worker(s)")

@everywhere begin
	using Pyrolysis
	using LinearSolve
	using KLU
	using OrdinaryDiffEq
	using LinearSolve: KLUFactorization
end

using Distributions
using QuasiMonteCarlo
using Statistics
using Printf
using Plots
using DelimitedFiles
using LinearAlgebra

# ============================================================================
# Parameter Space Definition — Correlated Multivariate Normal Groups
# ============================================================================
#
# Parameters are sampled jointly within correlated groups using MvNormal
# distributions. Covariance matrices from the regression JSON files in
# oxidative_pyrolysis/ (Bellamy, McKinnon, Stoliarov, energetics.tex).
#
# Eight correlated groups:
#   0a. ε_HFM (systematic bias)           — 1×1 N(0, 0.00577²)
#   0b. ε_STA (systematic bias)           — 1×1 N(0, 0.0261²)
#   0c. ε_water (systematic bias)          — 1×1 N(0, 0.01732²)
#   1. DF virgin cp segment 1 (a1, b1)    — 2×2 from HFM regression
#   2. DF virgin cp segment 3 (a3, b3)    — 2×2 from STA/DSC regression
#   3. Char cp (a, b)                      — 2×2 from independent char STA
#   4. Heat of sorption (a, b, c)          — 3×3 from WLS exponential fit
#   5. Heats of reaction (ΔH1–ΔH5)        — 5×5 conditional covariance (analytical delta method)
#
# Cp regressions were fit in J/(g·K); values and covariances are converted
# to J/(kg·K) below (×1e3 for values, ×1e6 for covariance entries).
# HoS parameters are in kJ/kg and 1/% (same units as build_and_solve).
# Heats of reaction conditional covariance is in kJ²/kg², converted to J²/kg².
#
# Complete uncertainty budget traceability (combined_uncertainty_budget.json):
#
# INCLUDED:
#   Type B (systematic):
#     ε_HFM  — HFM heat flow accuracy, 1%/√3 = 0.577% (manufacturer spec)     [Group 0a]
#     ε_STA  — DSC multi-point enthalpy calibration, 2.61% (7 standards)        [Group 0b]
#     ε_water — isothermal water calibration for HoS, 1.73% (22 replicates)     [Group 0c]
#   Type A (statistical/regression):
#     HFM cp regression (seg 1)              — 2×2 covariance                   [Group 1]
#     STA cp regression (seg 3)              — 2×2 covariance                   [Group 2]
#     Independent char STA cp regression     — 2×2 covariance                   [Group 3]
#     HoS WLS exponential fit                — 3×3 covariance                   [Group 4]
#     ΔH conditional (WLS + char cp + cross) — 5×5 covariance                   [Group 5]
#   Coupling:
#     Virgin cp → ΔH Jacobian correction     — J_HRXN_CP × δ(a₃,b₃) at runtime
#     ΔH × ε_STA multiplicative scaling      — applied after Jacobian correction
#
# OMITTED (justified):
#   STA temperature calibration on Cp        — 4.6% of STA Type B variance; negligible
#   STA temperature calibration on ΔH        — second-order: O(δT²/ΔT²_peak); negligible
#   STA mass measurement                     — u_rel = 0.03%; ~100× smaller than ε_STA
#   HFM temperature (Type K TC) on Cp        — additive intercept shift δCp = b·δT ≈ 0.004
#                                               J/(g·K), negligible vs regression SE of
#                                               intercept (0.218 J/(g·K))
#   HFM density systematic                   — u_rel = 0.06%; negligible vs heat flow
#   Kinetic parameters (A, E)                — held fixed (outside scope of this study)

"""
	CorrelatedGroup

Defines a group of parameters with a joint multivariate normal distribution.
`indices` maps into the 17-element θ vector passed to build_and_solve.
"""
struct CorrelatedGroup
	name::String
	indices::Vector{Int}
	μ::Vector{Float64}
	Σ::Matrix{Float64}
end

const PARAM_NAMES = [
	:epsilon_hfm, :epsilon_sta, :epsilon_water,                  # indices 1–3  (Type B biases)
	:df_cp_a1, :df_cp_b1,                                       # indices 4–5  (Type A, HFM)
	:df_cp_a3, :df_cp_b3,                                       # indices 6–7  (Type A, STA/DSC)
	:char_cp_a, :char_cp_b,                                     # indices 8–9  (Type A, char STA)
	:hos_baseline, :hos_amplitude, :hos_decay,                   # indices 10–12
	:h_rxn1, :h_rxn2, :h_rxn3, :h_rxn4, :h_rxn5,               # indices 13–17
]
const N_PARAMS = 17

const NOMINAL_VALUES = [
	0.0, 0.0, 0.0,                # ε_HFM, ε_STA, ε_water (zero-mean biases)
	-53.052, 4.2453,               # DF cp seg 1 [J/(kg·K)]  (HFM regression)
	-2131.33, 9.11878,             # DF cp seg 3 [J/(kg·K)]  (STA regression)
	347.212, 1.78423,              # Char cp [J/(kg·K)]       (independent char STA)
	2279.76, 786.962, 0.11436,    # HoS [kJ/kg, kJ/kg, 1/%] (WLS exponential fit)
	-19153.3, -5573.79, 266699.7, -219638.7, -12381.6,  # ΔH [J/kg]
]

# --- Group 0a: Type B systematic bias for HFM (ε_HFM) ---
# HFM heat flow accuracy ±1% (manufacturer spec), rectangular → u_rel = 1%/√3 = 0.577%.
# Temperature systematic (Type K TC, ±1.5°C → 0.87°C at 1σ) propagates through the
# Cp slope as an intercept shift (~0.29% relative at midrange) — small and partially
# degenerate with the Type A intercept uncertainty. Density systematic is negligible
# (<0.07%). Only the dominant multiplicative heat flow component is sampled here.
# Source: combined_uncertainty_budget.json → shared_systematics.HFM
const GROUP_EPSILON_HFM = CorrelatedGroup("eps_HFM", [1],
	[0.0],
	[0.005774^2;;])

# --- Group 0b: Type B systematic bias for STA/DSC (ε_STA) ---
# DSC multi-point enthalpy calibration (7 standards after RbNO₃ exclusion):
# u_rel = 2.61% (1σ). Temperature systematic (~0.3%) and mass (~0.03%) are
# negligible additions (RSS dominated by calibration at 95.4% of variance).
# Applied as (1 + ε_STA) to all DSC-derived Cp (virgin seg 3, char) and ΔH.
# Source: dsc_calibration_model.json, combined_uncertainty_budget.json
const GROUP_EPSILON_STA = CorrelatedGroup("eps_STA", [2],
	[0.0],
	[0.02613^2;;])

# --- Group 0c: Type B systematic bias for HoS water calibration (ε_water) ---
# HoS was measured isothermally on the STA, but calibrated independently via
# water-evaporation runs (S = |∫q dt| / (Δm·Lvap)), NOT the standard multi-point
# S(T) polynomial. This gives a separate, smaller calibration uncertainty.
# u_rel = 1.73% (1σ) from 22 replicates (11 DF + 11 RO), 21 dof.
# Applied as (1 + ε_water) to HoS parameters a (asymptote) and b (amplitude);
# c (decay constant) is a shape parameter unaffected by sensitivity scaling.
# Since (1+ε)(a + b·exp(-c·MC)) = (1+ε)·a + (1+ε)·b·exp(-c·MC), scaling a and b
# is equivalent to scaling the entire HoS expression.
# Source: isothermal_dsc_calibration_model.json, combined_uncertainty_budget.json
const GROUP_EPSILON_WATER = CorrelatedGroup("eps_water", [3],
	[0.0],
	[0.01732^2;;])

# --- Group 1: DF virgin cp segment 1 (HFM), T ≤ 323 K ---
# Source: douglas_fir_hfm/regression_results.json
# JSON units: J/(g·K) → script units: J/(kg·K), covariance ×1e6
# NOTE: This is the Type A (regression) covariance only. Type B is handled
# by ε_HFM above, applied as (1 + ε_HFM) × Cp in build_and_solve.
const GROUP_DF_CP_SEG1 = CorrelatedGroup("DF_cp_seg1", [4, 5],
	[-53.052, 4.2453],
	#         a1              b1          (intercept, slope)
	[0.21819^2 * 1e6   -1.5350e-4 * 1e6;
								   -1.5350e-4 * 1e6   7.038e-4^2 * 1e6])

# --- Group 2: DF virgin cp segment 3 (STA/DSC), T > 425 K ---
# Source: douglas_fir_sta/regression_results.json
# NOTE: Type A only. Type B handled by ε_STA.
const GROUP_DF_CP_SEG3 = CorrelatedGroup("DF_cp_seg3", [6, 7],
	[-2131.33, 9.11878],
	[0.056399^2 * 1e6    -7.1999e-6 * 1e6;
								   -7.1999e-6 * 1e6    1.2769e-4^2 * 1e6])

# --- Group 3: Char cp (independent STA measurement) ---
# Source: douglas_fir_independent_char_sta/regression_results.json
# NOTE: Type A only. Type B handled by ε_STA (same DSC instrument).
const GROUP_CHAR_CP = CorrelatedGroup("Char_cp", [8, 9],
	[347.212, 1.78423],
	[3.486e-3^2 * 1e6    -2.0047e-8 * 1e6;
								   -2.0047e-8 * 1e6    5.782e-6^2 * 1e6])

# --- Group 4: Heat of sorption (WLS exponential fit) ---
# Source: papers/energetics/figures/hos_fit_results.json
# Parameters: a=baseline [kJ/kg], b=amplitude [kJ/kg], c=decay [1/%]
# NOTE: Nominal values from WLS fit (hos_fit_results.json). Published as
#   a = 2280(13), b = 787(31), c = 0.114(10) in manuscript Eq. 18.
# The ~7% offset from known latent heat is a known bias (not random), not
# propagated here. Only the Type A fit covariance is sampled.
#
# Type B (calibration): HoS was measured isothermally on the STA using an
# independent water-evaporation calibration — NOT the standard multi-point
# DSC calibration. Therefore ε_STA does not apply to HoS. Instead, ε_water
# (Group 0c) provides the Type B calibration uncertainty from the isothermal
# water calibration reproducibility (u_rel = 1.73%, 22 replicates).
const GROUP_HOS = CorrelatedGroup("HoS", [10, 11, 12],
	[2279.76, 786.962, 0.11436],
	[173.957     22.025     0.09442;
									22.025    941.547     0.18431;
									 0.09442    0.18431   9.715e-5])

# --- Group 5: Heats of reaction (CONDITIONAL covariance from two-stage optimization) ---
# Source: scripts/energetics/results/douglas_fir_two_stage_optimized.json
#         → conditional_dh_uncertainty.covariance_matrix_kJ_per_kg_sq
#
# We use the *conditional* covariance from the analytical delta method (IFT +
# chain rule), which decomposes as:
#   Σ_ΔH|cp_virgin = Σ_WLS + J_char @ Σ_char @ J_charᵀ + cross terms
# where Σ_WLS is the Stage 2 weighted least-squares contribution and J_char
# propagates Stage 1 pseudo-char cp estimation uncertainty into ΔH.
# This covariance holds virgin cp at its nominal value, so Groups 2 (DF cp
# seg 3) and 5 (ΔH) are conditionally independent. The virgin cp → ΔH coupling
# is handled explicitly by the Jacobian correction J_HRXN_CP × δ(a₃,b₃) in
# build_and_solve. Using the full bootstrap covariance would double-count:
# virgin cp uncertainty would enter once through Group 2 (direct perturbation)
# and again through inflated ΔH variances (implicit perturbation baked in).
#
# Type B (DSC calibration): ΔH_fit = (1+ε)·ΔH_true (calibration scales, does not
# cancel). Applied as (1+ε_STA) multiplicative factor in build_and_solve.
#
# Units: kJ²/kg², converted to J²/kg² by ×1e6.
const GROUP_H_RXN = CorrelatedGroup("H_rxn", [13, 14, 15, 16, 17],
	[-19153.3, -5573.79, 266699.7, -219638.7, -12381.6],
	# Conditional covariance (analytical delta method): WLS + char cp + cross terms
	[3.42252    8.04561   11.3715    50.1601    12.3280;
								   8.04561   20.7989    28.7027   127.515     31.2934;
								  11.3715    28.7027    41.1541   178.306     45.5874;
								  50.1601   127.515    178.306    805.909    176.265;
								  12.3280    31.2934    45.5874   176.265    313.704] .* 1e6)

const CORR_GROUPS = [
	GROUP_EPSILON_HFM,
	GROUP_EPSILON_STA,
	GROUP_EPSILON_WATER,
	GROUP_DF_CP_SEG1,
	GROUP_DF_CP_SEG3,
	GROUP_CHAR_CP,
	GROUP_HOS,
	GROUP_H_RXN,
]

# ============================================================================
# Sampling Configuration
# ============================================================================

const N_SAMPLES = 2048      # Number of Sobol samples (power of 2 for optimal coverage)
const RESULTS_DIR = joinpath(@__DIR__, "results")

# ============================================================================
# Sobol Quasi-Random Sampling
# ============================================================================

"""
	generate_samples(groups, n_samples)

Generate Sobol quasi-random samples with correlated multivariate normal groups.

For each correlated group:
  1. Draw Sobol samples in [0,1]^k for the group's k dimensions
  2. Transform marginals to standard normal via Φ⁻¹
  3. Apply Cholesky factor of the group covariance to induce correlations:
	 x = μ + L·z,  where Σ = L·Lᵀ
This preserves low-discrepancy coverage in each marginal while reproducing the
target covariance structure.
"""
function generate_samples(groups::Vector{CorrelatedGroup}, n_samples::Int)
	n_params = sum(length(g.indices) for g in groups)

	# Sobol sequence in unit hypercube [0,1]^n_params
	sobol_unit = QuasiMonteCarlo.sample(n_samples, zeros(n_params), ones(n_params),
		SobolSample())

	samples = Matrix{Float64}(undef, n_params, n_samples)

	for group in groups
		idx = group.indices
		k = length(idx)
		L = cholesky(Symmetric(group.Σ)).L

		# Transform LHS uniform marginals → standard normal
		Z = Matrix{Float64}(undef, k, n_samples)
		for (gi, pi) in enumerate(idx)
			for j in 1:n_samples
				Z[gi, j] = quantile(Normal(), sobol_unit[pi, j])
			end
		end

		# Apply Cholesky to induce correlations: x = μ + L·z
		for j in 1:n_samples
			samples[idx, j] = group.μ + L * Z[:, j]
		end
	end

	return samples
end

# ============================================================================
# Simulation Function (runs on workers)
# ============================================================================

# Constants needed on all workers for the HoS closure
@everywhere const SA_MOISTURE_IDX = 7
@everywhere const SA_ρ_DRY_BASIS = 530.0

# Cumulative solid yield interpolation weights for intermediate cp's
# Reaction yields: 0.95, 0.80, 0.37, 0.79, 0.74
# Cumulative: int1=0.95, int2=0.76, int3=0.2812, int4=0.22215
# Final char yield: 0.95 * 0.80 * 0.37 * 0.79 * 0.74 = 0.164388...
@everywhere const FINAL_CHAR_YIELD = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
@everywhere const INTERP_PROGRESS = [
	(1 - 0.95) / (1 - FINAL_CHAR_YIELD),                   # int1: α ≈ 0.0598
	(1 - 0.95 * 0.80) / (1 - FINAL_CHAR_YIELD),            # int2: α ≈ 0.2872
	(1 - 0.95 * 0.80 * 0.37) / (1 - FINAL_CHAR_YIELD),    # int3: α ≈ 0.8602
	(1 - 0.95 * 0.80 * 0.37 * 0.79) / (1 - FINAL_CHAR_YIELD), # int4: α ≈ 0.9309
]

# Nominal virgin cp (high-T) regression coefficients — used as reference point
# for the Jacobian correction to ΔH.
@everywhere const NOM_CP_A3 = -2131.33
@everywhere const NOM_CP_B3 = 9.11878

# Jacobian J = ∂ΔH/∂(a₃, b₃): total sensitivity of fitted heats of reaction to
# the virgin cp (high-T) regression coefficients used as the optimization baseline.
# Source: conditional_dh_uncertainty.jacobian_dH_wrt_virgin_cp (analytical, IFT +
# chain rule through pseudo-char cp re-optimization). Includes both the direct
# baseline shift and the indirect path through Stage 1 pseudo-char cp.
# Rows = [ΔH₁, ΔH₂, ΔH₃, ΔH₄, ΔH₅], Cols = [a₃ (intercept), b₃ (slope)]
# Units col 1: (kJ/kg)/(J/(g·K)) = (J/kg)/(J/(kg·K)) → [K]
# Units col 2: (kJ/kg)/(J/(g·K²)) = (J/kg)/(J/(kg·K²)) → [K²]
# JSON column order was [slope, intercept]; swapped here to match [a₃, b₃].
@everywhere const J_HRXN_CP = [
	-101.8546   -50959.70;
	 -37.0298   -22263.28;
	 -12.3170    -7496.43;
	   1.17179    1237.13;
	  -3.57531   -2190.76]

@everywhere begin

	"""
		build_and_solve(θ)

	Build a Douglas Fir simulation with the given parameter vector θ (length 17),
	solve it, and return the PyrolysisSolution.
	"""
	function build_and_solve(θ::Vector{Float64})
		# Unpack the 17 uncertain parameters
		ε_hfm, ε_sta, ε_water,
		df_cp_a1, df_cp_b1, df_cp_a3, df_cp_b3,
		char_cp_a, char_cp_b,
		hos_baseline, hos_amplitude, hos_decay,
		h_rxn1, h_rxn2, h_rxn3, h_rxn4, h_rxn5 = θ

		# Type B multiplicative factors for each instrument
		scale_hfm   = 1.0 + ε_hfm    # HFM-derived Cp (segment 1)
		scale_sta   = 1.0 + ε_sta    # STA/DSC-derived Cp (segment 3, char) and ΔH
		scale_water = 1.0 + ε_water  # isothermal water calibration (HoS only)

		# --- Correct ΔH for virgin cp regression uncertainty (Type A coupling) ---
		# The two-stage optimization used nominal (a₃, b₃) as fixed input for the
		# reactant-side sensible heat baseline. When the regression coefficients
		# deviate from nominal (Type A), the fitted ΔH shifts by J × δ(a₃,b₃).
		δa3 = df_cp_a3 - NOM_CP_A3
		δb3 = df_cp_b3 - NOM_CP_B3
		h_rxn1_adj = h_rxn1 + J_HRXN_CP[1, 1] * δa3 + J_HRXN_CP[1, 2] * δb3
		h_rxn2_adj = h_rxn2 + J_HRXN_CP[2, 1] * δa3 + J_HRXN_CP[2, 2] * δb3
		h_rxn3_adj = h_rxn3 + J_HRXN_CP[3, 1] * δa3 + J_HRXN_CP[3, 2] * δb3
		h_rxn4_adj = h_rxn4 + J_HRXN_CP[4, 1] * δa3 + J_HRXN_CP[4, 2] * δb3
		h_rxn5_adj = h_rxn5 + J_HRXN_CP[5, 1] * δa3 + J_HRXN_CP[5, 2] * δb3

		# --- Correct ΔH for DSC calibration bias (Type B) ---
		# The optimization derived ΔH from DSC heat flow using DSC-measured cp as
		# baseline: q_meas = (1+ε)·q_true, cp_meas = (1+ε)·cp_true, so the residual
		# yields ΔH_fit = (1+ε)·ΔH_true. The calibration does NOT cancel — it scales
		# ΔH proportionally, just like cp. Both must be scaled by (1+ε_STA).
		h_rxn1_adj *= scale_sta
		h_rxn2_adj *= scale_sta
		h_rxn3_adj *= scale_sta
		h_rxn4_adj *= scale_sta
		h_rxn5_adj *= scale_sta

		# --- DF virgin heat capacity (piecewise linear, 3 segments) ---
		# Segment 1 (T ≤ 323 K): HFM instrument → scale by (1 + ε_HFM)
		# Segment 3 (T > 425 K): STA/DSC instrument → scale by (1 + ε_STA)
		# Segment 2 (323–425 K): interpolated from endpoints of seg 1 & 3
		#   so it inherits both biases through the endpoint values
		function virgin_cp(T)
			if T <= 323.0
				(df_cp_a1 + df_cp_b1 * T) * scale_hfm
			elseif T <= 425.0
				cp_low  = (df_cp_a1 + df_cp_b1 * 323.0) * scale_hfm
				cp_high = (df_cp_a3 + df_cp_b3 * 425.0) * scale_sta
				cp_low + (cp_high - cp_low) * (T - 323.0) / 102.0
			else
				(df_cp_a3 + df_cp_b3 * T) * scale_sta
			end
		end

		df_cp = CallableProperty(virgin_cp)

		# --- Char heat capacity (linear, STA instrument → scale by ε_STA) ---
		char_cp_func = CallableProperty(T -> (char_cp_a + char_cp_b * T) * scale_sta)

		# --- Intermediate cp's derived via cumulative-yield interpolation ---
		# char endpoint also scaled by scale_sta (same DSC instrument)
		int1_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[1]) * virgin_cp(T) +
										INTERP_PROGRESS[1] * (char_cp_a + char_cp_b * T) * scale_sta)

		int2_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[2]) * virgin_cp(T) +
										INTERP_PROGRESS[2] * (char_cp_a + char_cp_b * T) * scale_sta)

		int3_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[3]) * virgin_cp(T) +
										INTERP_PROGRESS[3] * (char_cp_a + char_cp_b * T) * scale_sta)

		int4_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[4]) * virgin_cp(T) +
										INTERP_PROGRESS[4] * (char_cp_a + char_cp_b * T) * scale_sta)

		# --- Moisture cp (FIXED — not uncertain) ---
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

		# --- Water vapor cp (FIXED) ---
		water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)

		# --- Heat of sorption (parametric) ---
		# HoS uses an independent isothermal water-evaporation calibration,
		# NOT the standard DSC multi-point calibration. ε_water (not ε_STA)
		# provides the Type B calibration uncertainty. Scaling a and b by
		# (1+ε_water) is equivalent to scaling the full expression since c is
		# a shape parameter (decay constant) unaffected by sensitivity scaling.
		hos = StateDependentProperty((T, ξ) -> begin
			MC = 100.0 * ξ[SA_MOISTURE_IDX] / SA_ρ_DRY_BASIS
			return (hos_baseline + hos_amplitude * exp(-hos_decay * MC)) * 1e3 * scale_water
		end)

		# --- Gas transfer coefficients and permeability (FIXED) ---
		λ_wood   = 1.0e-6
		λ_char   = 1.0e-5
		λ_gas    = 1.0e-5
		κ_virgin = 2.3e-15
		κ_char   = 1.0e-14

		# --- Build material ---
		material = Material(
			name = :DouglasFir,
			components = (
				# 1. Virgin Douglas Fir
				SolidComponent(:DF,
					ρ = 530.0, c = df_cp, k = 0.145,
					ε = 0.81, α = 0.0, λ = λ_wood, κ = κ_virgin),
				# 2. Intermediate 1
				SolidComponent(:DF_int1,
					ρ = 507.2, c = int1_cp, k = 0.292951,
					ε = 0.790965, α = 0.0, λ = λ_wood * 1.5, κ = 2.7e-15),
				# 3. Intermediate 2
				SolidComponent(:DF_int2,
					ρ = 419.8, c = int2_cp, k = 0.329938,
					ε = 0.786207, α = 0.0, λ = λ_wood * 3.0, κ = 3.1e-15),
				# 4. Intermediate 3
				SolidComponent(:DF_int3,
					ρ = 203.2, c = int3_cp, k = 0.644833,
					ε = 0.745694, α = 0.0, λ = λ_wood * 6.0, κ = 3.6e-15),
				# 5. Intermediate 4
				SolidComponent(:DF_int4,
					ρ = 176.6, c = int4_cp, k = 0.777700,
					ε = 0.728600, α = 0.0, λ = λ_char * 0.8, κ = 8.3e-15),
				# 6. Char
				SolidComponent(:Char,
					ρ = 150.0, c = char_cp_func, k = 1.0,
					ε = 0.70, α = 0.0, λ = λ_char, κ = κ_char),
				# 7. Moisture (liquid, pseudo-density)
				LiquidComponent(:MOISTURE,
					ρ = 1e6, c = moisture_cp, k = 0.145,
					ε = 0.81, α = 0.0, κ = 0.0),
				# 8-12. Gas products
				GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas * 0.85),
				GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas * 1.0),
				GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.6),
				GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas * 1.3),
				GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas * 1.4),
				# 13. Water vapor
				GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.6),
			),
			reactions = (
				Reaction(:pyrolysis_step1, 1 => (2, 8),
					A = 74816.95, E = 71.59e3, h = h_rxn1_adj, n = 1.0, yields = (0.95, 0.05)),
				Reaction(:pyrolysis_step2, 2 => (3, 9),
					A = 23550492839.0, E = 140.9e3, h = h_rxn2_adj, n = 1.0, yields = (0.80, 0.20)),
				Reaction(:pyrolysis_step3, 3 => (4, 10),
					A = 7.744618e12, E = 181.0e3, h = h_rxn3_adj, n = 1.0, yields = (0.37, 0.63)),
				Reaction(:pyrolysis_step4, 4 => (5, 11),
					A = 103276.14, E = 98.94e3, h = h_rxn4_adj, n = 1.0, yields = (0.79, 0.21)),
				Reaction(:pyrolysis_step5, 5 => (6, 12),
					A = 6.397348, E = 56.01e3, h = h_rxn5_adj, n = 1.0, yields = (0.74, 0.26)),
				Reaction(:moisture_evaporation, 7 => 13,
					A = 2.997e5, E = 5.0e4, h = hos, n = 1.0),
			),
			conductivity_mixing = WEIGHTED, β = 0.5,
			gas_transfer_mixing = WEIGHTED, β_λ = 0.5,
			permeability_mixing = WEIGHTED, β_κ = 0.5,
		)

		# --- Problem setup (all FIXED at nominal) ---
		sample_thickness = 12.0e-3
		n_cells = 240
		T_ambient = 300.0
		q_external = 50e3
		h_conv = 10.0
		moisture_content_db = 0.06
		ρ_DF_initial = 530.0
		ρ_moisture_initial = moisture_content_db * ρ_DF_initial

		mesh = create_uniform_mesh(
			sample_thickness, n_cells, Val(13);
			T_initial = T_ambient,
			ξ_initial = (ρ_DF_initial, 0.0, 0.0, 0.0, 0.0, 0.0, ρ_moisture_initial,
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
		)

		bc_set = BoundaryConditionSet(
			thermal_top = RadiativeFluxBC(flux = q_external) +
						  RadiativeBC(T_inf = 300) +
						  ConvectiveBC(h = h_conv, T_inf = 300),
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
			mass_bottom = ImpermeableBC(),
		)

		T_initial_func = z -> T_ambient
		ξ_initial_funcs = [
			z -> ρ_DF_initial,
			z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
			z -> ρ_moisture_initial,
			z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
		]

		problem = PyrolysisProblem(
			mesh = mesh,
			material = material,
			bc_set = bc_set,
			T_initial = T_initial_func,
			ξ_initial = ξ_initial_funcs,
			external_flux = t -> 0.0,
			ambient_temperature = t -> T_ambient,
			tspan = (0.0, 600.0),
			dt_initial = 1e-3,
			output_times = collect(range(0.0, 600.0, length = 3601)),
		)

		# --- Solve ---
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

end # @everywhere

# ============================================================================
# QoI Extraction (must be @everywhere — called inside pmap on workers)
# ============================================================================

@everywhere const MLR_IGNITION_THRESHOLD = 3e-3  # kg/(m²·s)

@everywhere struct SensitivityResult
	# Scalars
	peak_mlr::Float64
	time_to_peak_mlr::Float64
	time_to_ignition::Float64
	total_mass_loss_pct::Float64
	peak_surface_temp::Float64
	peak_back_temp::Float64
	final_thickness::Float64
	retcode::Symbol
	solve_time::Float64
	# Curves
	times::Vector{Float64}
	mlr_curve::Vector{Float64}
	surface_temp_curve::Vector{Float64}
	back_temp_curve::Vector{Float64}
end

const SCALAR_QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr", "time_to_ignition",
	"total_mass_loss_pct", "peak_surface_temp", "peak_back_temp",
	"final_thickness", "solve_time",
]

@everywhere function extract_qois(sol)
	t = sol.times
	mlr = sol.surface_mass_flux
	T_surf = sol.surface_temperature
	T_back = sol.back_temperature

	peak_mlr_val = maximum(mlr)
	peak_idx = argmax(mlr)
	time_to_peak = t[peak_idx]

	tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
	time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf

	total_mass_loss = (1.0 - sol.total_mass[end] / sol.total_mass[1]) * 100.0

	return SensitivityResult(
		peak_mlr_val, time_to_peak, time_to_ignition,
		total_mass_loss,
		maximum(T_surf), maximum(T_back),
		sol.thickness_history[end],
		sol.retcode, sol.solve_time,
		t, mlr, T_surf, T_back,
	)
end

function scalar_qoi_vector(r::SensitivityResult)
	return [r.peak_mlr, r.time_to_peak_mlr, r.time_to_ignition,
		r.total_mass_loss_pct, r.peak_surface_temp, r.peak_back_temp,
		r.final_thickness, r.solve_time]
end

# ============================================================================
# I/O Helpers
# ============================================================================

@everywhere function save_curve_csv(filepath, times, values)
	open(filepath, "w") do io
		println(io, "time,value")
		for (t, v) in zip(times, values)
			println(io, "$t,$v")
		end
	end
end

function save_parameters_csv(filepath, param_names, samples)
	open(filepath, "w") do io
		header = "sample_id," * join(string.(param_names), ",")
		println(io, header)
		for j in 1:size(samples, 2)
			row = "$j," * join([@sprintf("%.10g", samples[i, j]) for i in 1:size(samples, 1)], ",")
			println(io, row)
		end
	end
end

function save_scalar_qois_csv(filepath, sample_ids, results)
	open(filepath, "w") do io
		header = "sample_id,retcode," * join(SCALAR_QOI_NAMES, ",")
		println(io, header)
		for (id, r) in zip(sample_ids, results)
			if r !== nothing
				qois = scalar_qoi_vector(r)
				row = "$id,$(r.retcode)," * join([@sprintf("%.10g", q) for q in qois], ",")
				println(io, row)
			else
				row = "$id,FAILED," * join(["NaN" for _ in SCALAR_QOI_NAMES], ",")
				println(io, row)
			end
		end
	end
end

# ============================================================================
# Sensitivity Analysis: Standardized Rank Regression Coefficients
# ============================================================================

"""
	compute_srrc(X, y)

Compute Standardized Rank Regression Coefficients (SRRC).
X is [n_samples × n_params], y is [n_samples].
Returns (β, R²) where β is the vector of SRRC values (one per parameter)
and R² is the coefficient of determination of the rank regression (measures
how well a linear-in-ranks model fits; low R² indicates non-monotonic effects).
"""
function compute_srrc(X::Matrix{Float64}, y::Vector{Float64})
	n, p = size(X)

	# Rank-transform both X and y
	X_ranked = similar(X)
	for j in 1:p
		X_ranked[:, j] = Float64.(sortperm(sortperm(X[:, j])))
	end
	y_ranked = Float64.(sortperm(sortperm(y)))

	# Standardize ranks to zero mean, unit variance
	for j in 1:p
		μ = mean(X_ranked[:, j])
		σ = std(X_ranked[:, j])
		X_ranked[:, j] .= (X_ranked[:, j] .- μ) ./ σ
	end
	μ_y = mean(y_ranked)
	σ_y = std(y_ranked)
	y_ranked .= (y_ranked .- μ_y) ./ σ_y

	# Ordinary least squares on ranked data
	β = X_ranked \ y_ranked

	# R² of the rank regression
	y_pred = X_ranked * β
	ss_res = sum((y_ranked .- y_pred) .^ 2)
	ss_tot = sum(y_ranked .^ 2)  # already zero-mean
	R2 = 1.0 - ss_res / ss_tot

	return β, R2
end

# ============================================================================
# Main Execution
# ============================================================================

function main()
	println("\n" * "="^70)
	println("  ThermaKin3 Sensitivity Analysis — Douglas Fir Cone Calorimeter")
	println("="^70)
	println("  Parameters: $N_PARAMS uncertain")
	println("  Samples:    $N_SAMPLES (Sobol)")
	println("  Workers:    $(nworkers())")
	println("  Est. time:  ~$(round(N_SAMPLES * 7.0 / 60 / nworkers(), digits=1)) hours")
	println("="^70 * "\n")

	# --- Create output directories ---
	mkpath(RESULTS_DIR)
	mkpath(joinpath(RESULTS_DIR, "mlr_curves"))
	mkpath(joinpath(RESULTS_DIR, "surface_temp_curves"))
	mkpath(joinpath(RESULTS_DIR, "back_temp_curves"))

	# --- Generate or reload samples ---
	# NOTE: We always generate the full Sobol sequence of length N_SAMPLES to
	# preserve the low-discrepancy property. Appending a fresh sequence of
	# length (N_SAMPLES - n_existing) would break the space-filling guarantee.
	param_file = joinpath(RESULTS_DIR, "parameters.csv")
	if isfile(param_file)
		println("Reloading existing samples from $param_file...")
		param_csv = readdlm(param_file, ','; header = true)
		samples = Matrix{Float64}(param_csv[1][:, 2:end]')  # [n_params × n_samples]
		n_existing = size(samples, 2)
		println("  Loaded $n_existing samples with $(size(samples, 1)) parameters")

		if n_existing < N_SAMPLES
			println("  Need $N_SAMPLES samples but only have $n_existing.")
			println("  Regenerating full Sobol sequence of length $N_SAMPLES...")
			samples = generate_samples(CORR_GROUPS, N_SAMPLES)
			# Verify existing samples match the first n_existing columns of
			# the regenerated sequence (they should, since Sobol is deterministic)
			param_csv_check = Matrix{Float64}(param_csv[1][:, 2:end]')
			max_diff = maximum(abs.(samples[:, 1:n_existing] .- param_csv_check))
			if max_diff > 1e-10
				@warn "Regenerated Sobol sequence differs from saved samples " *
					  "(max |Δ| = $max_diff). Overwriting with fresh sequence."
			end
			save_parameters_csv(param_file, PARAM_NAMES, samples)
			println("  Updated $param_file ($n_existing → $N_SAMPLES samples)")
		end
	else
		println("Generating $N_SAMPLES Sobol samples for $N_PARAMS parameters...")
		samples = generate_samples(CORR_GROUPS, N_SAMPLES)
		save_parameters_csv(param_file, PARAM_NAMES, samples)
		println("  Parameter samples saved to $param_file")
	end

	# --- Check for existing results (resumability) ---
	completed_ids = Set{Int}()

	# Check scalar_qois.csv for completed runs
	scalar_file = joinpath(RESULTS_DIR, "scalar_qois.csv")
	if isfile(scalar_file)
		lines = readlines(scalar_file)
		for line in lines[2:end]  # skip header
			parts = split(line, ",")
			if length(parts) > 1 && parts[2] != "FAILED"
				push!(completed_ids, parse(Int, parts[1]))
			end
		end
	end

	# Also check for saved curve files (covers crash before scalar CSV was written)
	for i in 1:size(samples, 2)
		mlr_file = joinpath(RESULTS_DIR, "mlr_curves", @sprintf("run_%04d.csv", i))
		surf_file = joinpath(RESULTS_DIR, "surface_temp_curves", @sprintf("run_%04d.csv", i))
		back_file = joinpath(RESULTS_DIR, "back_temp_curves", @sprintf("run_%04d.csv", i))
		if isfile(mlr_file) && isfile(surf_file) && isfile(back_file)
			push!(completed_ids, i)
		end
	end

	if !isempty(completed_ids)
		println("  Found $(length(completed_ids)) completed runs — resuming...")
	end

	# Backfill scalar_qois.csv from curve files for runs not yet in the scalar CSV
	scalar_csv_ids = Set{Int}()
	if isfile(scalar_file)
		for line in readlines(scalar_file)[2:end]
			parts = split(line, ",")
			if length(parts) > 0
				push!(scalar_csv_ids, parse(Int, parts[1]))
			end
		end
	end
	backfill_ids = sort(collect(setdiff(completed_ids, scalar_csv_ids)))
	if !isempty(backfill_ids)
		println("  Backfilling $(length(backfill_ids)) scalar QoIs from saved curves...")
		need_header = !isfile(scalar_file)
		open(scalar_file, "a") do io
			if need_header
				println(io, "sample_id,retcode," * join(SCALAR_QOI_NAMES, ","))
			end
			for id in backfill_ids
				try
					mlr_data = readdlm(joinpath(RESULTS_DIR, "mlr_curves", @sprintf("run_%04d.csv", id)), ','; header = true)[1]
					surf_data = readdlm(joinpath(RESULTS_DIR, "surface_temp_curves", @sprintf("run_%04d.csv", id)), ','; header = true)[1]
					back_data = readdlm(joinpath(RESULTS_DIR, "back_temp_curves", @sprintf("run_%04d.csv", id)), ','; header = true)[1]

					t = Float64.(mlr_data[:, 1])
					mlr = Float64.(mlr_data[:, 2])
					T_surf = Float64.(surf_data[:, 2])
					T_back = Float64.(back_data[:, 2])

					peak_mlr_val = maximum(mlr)
					time_to_peak = t[argmax(mlr)]
					tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
					time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf
					total_mass_loss = NaN  # not recoverable from curves alone
					peak_surface_temp = maximum(T_surf)
					peak_back_temp = maximum(T_back)
					final_thickness = NaN  # not recoverable from curves alone
					solve_time = NaN

					qois = [peak_mlr_val, time_to_peak, time_to_ignition,
						total_mass_loss, peak_surface_temp, peak_back_temp,
						final_thickness, solve_time]
					row = "$id,Success," * join([@sprintf("%.10g", q) for q in qois], ",")
					println(io, row)
				catch e
					@warn "Could not backfill sample $id" exception=(e, catch_backtrace())
					row = "$id,FAILED," * join(["NaN" for _ in SCALAR_QOI_NAMES], ",")
					println(io, row)
				end
			end
		end
		println("  Backfill complete.")
	end

	remaining_ids = [i for i in 1:N_SAMPLES if i ∉ completed_ids]
	println("  Running $(length(remaining_ids)) remaining simulations...\n")

	if isempty(remaining_ids)
		println("All samples already completed. Skipping to analysis.\n")
	else
		# --- Run simulations in parallel ---
		t_start = time()
		n_done = Threads.Atomic{Int}(0)

		results = pmap(remaining_ids) do i
			θ = samples[:, i]

			try
				sol = build_and_solve(θ)
				result = extract_qois(sol)

				# Save curves immediately
				save_curve_csv(
					joinpath(RESULTS_DIR, "mlr_curves", @sprintf("run_%04d.csv", i)),
					result.times, result.mlr_curve)
				save_curve_csv(
					joinpath(RESULTS_DIR, "surface_temp_curves", @sprintf("run_%04d.csv", i)),
					result.times, result.surface_temp_curve)
				save_curve_csv(
					joinpath(RESULTS_DIR, "back_temp_curves", @sprintf("run_%04d.csv", i)),
					result.times, result.back_temp_curve)

				return (i, result)
			catch e
				@warn "Sample $i FAILED" exception=(e, catch_backtrace())
				return (i, nothing)
			end
		end

		elapsed = time() - t_start
		n_success = count(r -> r[2] !== nothing, results)
		n_fail = count(r -> r[2] === nothing, results)
		println("\nBatch complete:")
		println("  Succeeded: $n_success / $(length(remaining_ids))")
		println("  Failed:    $n_fail")
		println("  Elapsed:   $(round(elapsed / 60, digits=1)) minutes")
		println("  Avg/run:   $(round(elapsed / length(remaining_ids), digits=1)) seconds")

		# --- Save scalar QoIs (append mode for resumability) ---
		new_ids = [r[1] for r in results]
		new_results = [r[2] for r in results]

		if isfile(scalar_file)
			# Append new results
			open(scalar_file, "a") do io
				for (id, r) in zip(new_ids, new_results)
					if r !== nothing
						qois = scalar_qoi_vector(r)
						row = "$id,$(r.retcode)," * join([@sprintf("%.10g", q) for q in qois], ",")
						println(io, row)
					else
						row = "$id,FAILED," * join(["NaN" for _ in SCALAR_QOI_NAMES], ",")
						println(io, row)
					end
				end
			end
		else
			save_scalar_qois_csv(scalar_file, new_ids, new_results)
		end
	end

	# ========================================================================
	# Post-Processing: Load all results
	# ========================================================================

	println("\n" * "="^70)
	println("  Post-Processing & Sensitivity Analysis")
	println("="^70 * "\n")

	# Re-read all scalar QoIs
	scalar_data = readdlm(scalar_file, ','; header = true)
	data_matrix = scalar_data[1]
	col_names = vec(scalar_data[2])

	# Filter successful runs
	retcodes = data_matrix[:, 2]
	success_mask = retcodes .!= "FAILED"
	n_success_total = sum(success_mask)
	println("Total successful runs: $n_success_total / $N_SAMPLES")

	sample_ids = Int.(data_matrix[success_mask, 1])
	qoi_data = Float64.(data_matrix[success_mask, 3:end])  # skip id and retcode

	# Parameter values for successful runs
	param_data = Matrix(samples[:, sample_ids]')  # [n_success × n_params]

	# --- SRRC Sensitivity Analysis (screening diagnostic) ---
	# NOTE: SRRC is a rank-based screening tool that assumes monotonic
	# input-output relationships. It does NOT decompose variance. For
	# publication-quality variance attribution, use the PCE-based Sobol
	# indices in sensitivity_sobol.jl instead.
	println("\nSRRC Screening Analysis (Standardized Rank Regression Coefficients):")
	println("-"^70)

	srrc_results = Dict{String, Vector{Float64}}()
	srrc_r2 = Dict{String, Float64}()

	# Skip solve_time (not physical) and time_to_ignition (quantized/threshold)
	skip_qois = Set(["solve_time", "time_to_ignition"])

	for (qi, qname) in enumerate(SCALAR_QOI_NAMES)
		if qname ∈ skip_qois
			continue
		end

		y = qoi_data[:, qi]

		# Skip if all values are identical or contain Inf/NaN
		if all(isfinite, y) && std(y) > 0
			srrc, rank_r2 = compute_srrc(param_data, y)
			srrc_results[qname] = srrc
			srrc_r2[qname] = rank_r2

			println(@sprintf("\n  %s (rank R² = %.4f):", qname, rank_r2))
			if rank_r2 < 0.7
				println("    ⚠  Low rank R² — non-monotonic or higher-order effects likely")
			end
			sorted = sortperm(abs.(srrc), rev = true)
			for (rank, j) in enumerate(sorted)
				bar = repeat("█", round(Int, abs(srrc[j]) * 30))
				sign_char = srrc[j] >= 0 ? "+" : "-"
				println(@sprintf("    %2d. %-18s %s%.4f  %s",
					rank, PARAM_NAMES[j], sign_char, abs(srrc[j]), bar))
			end
		else
			println("\n  $qname: SKIPPED (constant or non-finite values)")
		end
	end

	# Save SRRC results
	srrc_file = joinpath(RESULTS_DIR, "srrc_sensitivity.csv")
	open(srrc_file, "w") do io
		qoi_keys = sort(collect(keys(srrc_results)))
		header = "parameter," * join(qoi_keys, ",")
		println(io, header)
		for (pi, pname) in enumerate(PARAM_NAMES)
			row = "$pname," * join(
				[haskey(srrc_results, qname) ? @sprintf("%.6f", srrc_results[qname][pi]) : "NaN"
				 for qname in qoi_keys], ",")
			println(io, row)
		end
		# Append rank R² as a summary row
		println(io, "rank_R2," * join(
			[@sprintf("%.6f", srrc_r2[q]) for q in qoi_keys], ","))
	end
	println("\n  SRRC results saved to $srrc_file")

	# ========================================================================
	# Visualization
	# ========================================================================

	println("\nGenerating plots...")

	# --- 1. MLR Uncertainty Envelope ---
	println("  Loading MLR curves...")
	all_mlr = Matrix{Float64}(undef, 3601, 0)
	common_times = nothing
	for id in sample_ids
		fpath = joinpath(RESULTS_DIR, "mlr_curves", @sprintf("run_%04d.csv", id))
		if isfile(fpath)
			data = readdlm(fpath, ','; header = true)[1]
			if common_times === nothing
				common_times = data[:, 1]
			end
			all_mlr = hcat(all_mlr, data[:, 2])
		end
	end

	if common_times !== nothing && size(all_mlr, 2) > 10
		n_t = length(common_times)
		p5  = [quantile(all_mlr[i, :], 0.05) for i in 1:n_t]
		p25 = [quantile(all_mlr[i, :], 0.25) for i in 1:n_t]
		p50 = [quantile(all_mlr[i, :], 0.50) for i in 1:n_t]
		p75 = [quantile(all_mlr[i, :], 0.75) for i in 1:n_t]
		p95 = [quantile(all_mlr[i, :], 0.95) for i in 1:n_t]

		plt_mlr = plot(common_times, p50 .* 1000, label = "Median", lw = 2.5, color = :blue,
			xlabel = "Time [s]", ylabel = "MLR [g/(m²·s)]",
			title = "Mass Loss Rate — Uncertainty Envelope ($n_success_total samples)")
		plot!(common_times, p5 .* 1000, fillrange = p95 .* 1000,
			alpha = 0.15, label = "90% CI", color = :blue, lw = 0)
		plot!(common_times, p25 .* 1000, fillrange = p75 .* 1000,
			alpha = 0.3, label = "50% CI", color = :blue, lw = 0)

		savefig(plt_mlr, joinpath(RESULTS_DIR, "mlr_uncertainty_envelope.png"))
		println("  Saved mlr_uncertainty_envelope.png")

		# --- 2. Surface Temperature Envelope ---
		println("  Loading surface temperature curves...")
		all_tsurf = Matrix{Float64}(undef, 3601, 0)
		for id in sample_ids
			fpath = joinpath(RESULTS_DIR, "surface_temp_curves", @sprintf("run_%04d.csv", id))
			if isfile(fpath)
				data = readdlm(fpath, ','; header = true)[1]
				all_tsurf = hcat(all_tsurf, data[:, 2])
			end
		end

		if size(all_tsurf, 2) > 10
			p50_ts = [quantile(all_tsurf[i, :], 0.50) for i in 1:n_t]
			p5_ts  = [quantile(all_tsurf[i, :], 0.05) for i in 1:n_t]
			p95_ts = [quantile(all_tsurf[i, :], 0.95) for i in 1:n_t]
			p25_ts = [quantile(all_tsurf[i, :], 0.25) for i in 1:n_t]
			p75_ts = [quantile(all_tsurf[i, :], 0.75) for i in 1:n_t]

			plt_ts = plot(common_times, p50_ts, label = "Median", lw = 2.5, color = :red,
				xlabel = "Time [s]", ylabel = "Temperature [K]",
				title = "Surface Temperature — Uncertainty Envelope")
			plot!(common_times, p5_ts, fillrange = p95_ts, alpha = 0.15, label = "90% CI", color = :red, lw = 0)
			plot!(common_times, p25_ts, fillrange = p75_ts, alpha = 0.3, label = "50% CI", color = :red, lw = 0)

			savefig(plt_ts, joinpath(RESULTS_DIR, "surface_temp_uncertainty_envelope.png"))
			println("  Saved surface_temp_uncertainty_envelope.png")
		end

		# --- 3. Back Temperature Envelope ---
		println("  Loading back temperature curves...")
		all_tback = Matrix{Float64}(undef, 3601, 0)
		for id in sample_ids
			fpath = joinpath(RESULTS_DIR, "back_temp_curves", @sprintf("run_%04d.csv", id))
			if isfile(fpath)
				data = readdlm(fpath, ','; header = true)[1]
				all_tback = hcat(all_tback, data[:, 2])
			end
		end

		if size(all_tback, 2) > 10
			p50_tb = [quantile(all_tback[i, :], 0.50) for i in 1:n_t]
			p5_tb  = [quantile(all_tback[i, :], 0.05) for i in 1:n_t]
			p95_tb = [quantile(all_tback[i, :], 0.95) for i in 1:n_t]
			p25_tb = [quantile(all_tback[i, :], 0.25) for i in 1:n_t]
			p75_tb = [quantile(all_tback[i, :], 0.75) for i in 1:n_t]

			plt_tb = plot(common_times, p50_tb, label = "Median", lw = 2.5, color = :green,
				xlabel = "Time [s]", ylabel = "Temperature [K]",
				title = "Back Temperature — Uncertainty Envelope")
			plot!(common_times, p5_tb, fillrange = p95_tb, alpha = 0.15, label = "90% CI", color = :green, lw = 0)
			plot!(common_times, p25_tb, fillrange = p75_tb, alpha = 0.3, label = "50% CI", color = :green, lw = 0)

			savefig(plt_tb, joinpath(RESULTS_DIR, "back_temp_uncertainty_envelope.png"))
			println("  Saved back_temp_uncertainty_envelope.png")
		end
	end

	# --- 4. Scalar QoI Histograms ---
	println("  Generating scalar QoI histograms...")
	hist_qois = ["peak_mlr", "time_to_ignition", "total_mass_loss_pct", "peak_surface_temp"]
	hist_labels = ["Peak MLR [kg/(m²·s)]", "Time to Ignition [s]",
		"Total Mass Loss [%]", "Peak Surface Temp [K]"]
	hist_indices = [findfirst(==(q), SCALAR_QOI_NAMES) for q in hist_qois]

	plots_hist = []
	for (qi, label) in zip(hist_indices, hist_labels)
		y = qoi_data[:, qi]
		valid = isfinite.(y)
		if sum(valid) > 10
			push!(plots_hist, histogram(y[valid], xlabel = label, ylabel = "Count",
				title = label, legend = false, bins = 30, alpha = 0.7, color = :steelblue))
		end
	end

	if !isempty(plots_hist)
		plt_hist = plot(plots_hist..., layout = (2, 2), size = (1000, 800))
		savefig(plt_hist, joinpath(RESULTS_DIR, "scalar_qoi_histograms.png"))
		println("  Saved scalar_qoi_histograms.png")
	end

	# --- 5. Tornado Charts (one per QoI) ---
	println("  Generating tornado charts...")
	for (qname, srrc_vals) in srrc_results
		sorted_idx = sortperm(abs.(srrc_vals), rev = true)

		xmax = max(1.0, 1.1 * maximum(abs.(srrc_vals)))
		plt_tornado = bar(1:N_PARAMS,
			srrc_vals[sorted_idx],
			orientation = :horizontal,
			yticks = (1:N_PARAMS, string.(PARAM_NAMES[sorted_idx])),
			xlims = (-xmax, xmax),
			xlabel = "SRRC", title = "Parameter Sensitivity — $qname",
			legend = false, color = :steelblue, alpha = 0.8,
			size = (800, 600), left_margin = 15Plots.mm)
		vline!([0], color = :black, lw = 1)

		fname = "tornado_$(qname).png"
		savefig(plt_tornado, joinpath(RESULTS_DIR, fname))
		println("  Saved $fname")
	end

	println("\n" * "="^70)
	println("  Sensitivity analysis complete!")
	println("  Results in: $RESULTS_DIR")
	println("="^70)
end

# Run
main()
