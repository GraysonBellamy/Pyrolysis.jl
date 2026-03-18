"""
Sensitivity Analysis of Douglas Fir Cone Calorimeter Simulation

Propagates uncertainties in thermophysical input parameters through the
07_douglas_fir_cone.jl simulation using Sobol quasi-random sampling.

Uncertain parameters (17 total):
  - Systematic (Type B) instrument biases: ε_HFM, ε_DSC, ε_HoS (3)
	  Multiplicative calibration biases that scale measurements from each
	  instrument uniformly. Following GUM (JCGM 100:2008):
	  ε_HFM ~ N(0, 0.00462²)  — HFM heat flow calibration (data-driven, SRM 1453/1453e)
	  ε_DSC ~ N(0, 0.01132²)  — DSC enthalpy calibration (certified-only RMS 1.13%, dof=20)
	  ε_HoS ~ N(0, 0.00597²)  — HoS asymptote calibration (u_rel = 0.60%)
	  ε_DSC applies to DSC-derived Cp (virgin high-T, char) and ΔH.
	  ε_HFM applies to HFM-derived Cp (virgin low-T).
	  ε_HoS applies as multiplicative bias on HoS parameters.
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
Type B parameters are sampled independently per instrument (HFM, DSC, HoS) and
applied as multiplicative factors (1 + ε) to all values from that source.
HoS calibration uncertainty is now a separate Type B group (ε_HoS), with the
CR2 covariance (conditional on f_cal exact) used for the fit parameters.

Kinetic parameters (A, E) and all other properties are held fixed.

Outputs:
  - Scalar QoIs: peak MLR, time to ignition, total mass loss, temperature rise times
  - Time-series: MLR curve, surface temperature, back temperature per run
  - Sensitivity ranking via Standardized Rank Regression Coefficients (SRRC)
  - Uncertainty envelope plots (median, 50% CI, 90% CI)
  - Uncertainty envelope CSVs: mlr_envelope.csv, surface_temp_envelope.csv,
	back_temp_envelope.csv (columns: time, p5, p25, p50, p75, p95)

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
#   1. ε_HFM (systematic bias)              — 1×1 N(0, 0.00462²)
#   2. ε_DSC (systematic bias)              — 1×1 N(0, 0.01132²)
#   3. ε_HoS (systematic bias)              — 1×1 N(0, 0.00597²)
#   4. DF virgin cp low-T (a1, b1)          — 2×2 from HFM regression
#   5. DF virgin cp high-T (a3, b3)         — 2×2 per-curve empirical marginal
#   6. Char cp (a, b)                        — 2×2 from independent char STA
#   7. Heat of sorption (a, b, c)            — 3×3 CR2 covariance (conditional on f_cal)
#   8. Heats of reaction (ΔH1–ΔH5)          — 5×5 conditional covariance
#
# Cp regressions were fit in kJ/(kg·K); values and covariances are converted
# to J/(kg·K) below (×1e3 for values, ×1e6 for covariance entries).
# HoS parameters are in kJ/kg and 1/% (same units as build_and_solve).
# Heats of reaction conditional covariance is in kJ²/kg², converted to J²/kg².
#
# Complete uncertainty budget traceability (combined_uncertainty_budget.json):
#
# INCLUDED:
#   Type B (systematic):
#     ε_HFM  — HFM heat flow calibration, 0.462% (data-driven, SRM cal.)      [Group 1]
#     ε_DSC  — DSC enthalpy calibration, 1.13% RMS (certified-only, dof=20)     [Group 2]
#     ε_HoS  — HoS asymptote calibration, u_rel = 0.60%                         [Group 3]
#   Type A (statistical/regression):
#     HFM cp regression (low-T)              — 2×2 covariance                   [Group 4]
#     STA cp per-curve (high-T)              — 2×2 covariance                   [Group 5]
#     Independent char STA cp regression     — 2×2 covariance                   [Group 6]
#     HoS CR2 (conditional on f_cal)         — 3×3 covariance                   [Group 7]
#     ΔH conditional (per-curve + char cp)   — 5×5 covariance                   [Group 8]
#   Coupling:
#     Virgin cp → ΔH Jacobian correction     — J_HRXN_CP × δ(slope,int) at runtime
#     ΔH × ε_DSC multiplicative scaling      — applied after Jacobian correction
#     HoS × ε_HoS multiplicative scaling     — applied to HoS output
#
# OMITTED (justified):
#   STA temperature calibration on Cp        — 4.6% of STA Type B variance; negligible
#   STA temperature calibration on ΔH        — second-order: O(δT²/ΔT²_peak); negligible
#   STA mass measurement                     — u_rel = 0.03%; ~100× smaller than ε_DSC
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
	:epsilon_hfm, :epsilon_dsc, :epsilon_hos,                   # indices 1–3  (Type B biases)
	:df_cp_a1, :df_cp_b1,                                       # indices 4–5  (Type A, HFM)
	:df_cp_a3, :df_cp_b3,                                       # indices 6–7  (Type A, STA per-curve)
	:char_cp_a, :char_cp_b,                                     # indices 8–9  (Type A, independent char STA)
	:hos_a, :hos_b, :hos_c,                                     # indices 10–12 (CR2 cov, conditional on f_cal)
	:h_rxn1, :h_rxn2, :h_rxn3, :h_rxn4, :h_rxn5,               # indices 13–17
]
const N_PARAMS = 17

const NOMINAL_VALUES = [
	0.0, 0.0, 0.0,                    # ε_HFM, ε_DSC, ε_HoS (zero-mean biases)
	4.133836600273531, -18.01630709338176,                       # DF cp low-T [J/(kg·K)]  (HFM regression, ×1e3)
	9.136799235573225, -2139.3044552435705,                      # DF cp high-T [J/(kg·K)] (STA per-curve, ×1e3)
	1.8372509265710314, 312.4273785167445,                       # Char cp [J/(kg·K)]       (independent char STA, ×1e3)
	2450.420909987441, 875.6441221735336, 0.1121684273316275,    # HoS [kJ/kg, kJ/kg, 1/%] (CR2 fit)
	-19200.66288032251, -5263.201466602739, 264424.79037373835, -217151.69043762444, -326.40885469229025,  # ΔH [J/kg] (×1e3)
]

# --- Group 1: Type B systematic bias for HFM (ε_HFM) ---
# HFM heat flow calibration: data-driven from NIST SRM 1453/1453e at 3 thicknesses,
# multiDX piecewise-linear interpolation to DF sample thickness (19.4 mm).
# u_rel = 0.462% (1σ, RSS of replicate scatter + NIST certified uncertainty).
# Type K TC offset absorbed into f_cal — not a separate term.
# Source: hfm_calibration_uncertainty.py → interpolate_cal_at_thickness(DF)
const GROUP_EPSILON_HFM = CorrelatedGroup("eps_HFM", [1],
	[0.0],
	[0.0046154243860687^2;;])

# --- Group 2: Type B systematic bias for DSC (ε_DSC) ---
# DSC enthalpy calibration (certified-only RMS):
# u_rel = 1.13% (1σ, pooled across 5 calibration epochs, dof = 20).
# Applied as (1 + ε_DSC) to all DSC-derived Cp (virgin high-T, char) and ΔH.
# Source: dsc_calibration_model.json, gsa_parameter_config_douglas_fir.json
const GROUP_EPSILON_DSC = CorrelatedGroup("eps_DSC", [2],
	[0.0],
	[0.011322211895623097^2;;])

# --- Group 3: Type B systematic bias for HoS calibration (ε_HoS) ---
# Multiplicative bias on HoS from asymptote calibration:
# f_cal = L_vap(T_iso) / â_uncal, u_rel = 0.60%.
# Previously embedded in the total HoS covariance; now separated so that the
# fit covariance (Group 7) uses CR2 (conditional on f_cal exact).
# Source: gsa_parameter_config_douglas_fir.json
const GROUP_EPSILON_HOS = CorrelatedGroup("eps_HoS", [3],
	[0.0],
	[0.005974653346763125^2;;])

# --- Group 4: DF virgin cp low-T (HFM), T ≤ 323 K ---
# Source: douglas_fir_hfm/regression_results.json
# JSON units: kJ/(kg·K) → script units: J/(kg·K), mean ×1e3, covariance ×1e6
# NOTE: This is the Type A (regression) covariance only. Type B is handled
# by ε_HFM above, applied as (1 + ε_HFM) × Cp in build_and_solve.
# dof = 11 (12 specimens, N-1)
const GROUP_DF_CP_SEG1 = CorrelatedGroup("DF_cp_seg1", [4, 5],
	[0.004133836600273531 * 1e3, -0.01801630709338176 * 1e3],
	#         slope                   intercept
	[7.183142841163425e-05^2 * 1e6    -1.1925585013168409e-06 * 1e6;
	 -1.1925585013168409e-06 * 1e6    0.01850672773953015^2 * 1e6])

# --- Group 5: DF virgin cp high-T (STA per-curve), T > 425 K ---
# Source: gsa_parameter_config_douglas_fir.json
# Per-curve empirical marginal covariance (NOT the old WLS regression).
# JSON units: kJ/(kg·K) → script units: J/(kg·K), mean ×1e3, covariance ×1e6
# NOTE: Type A only. Type B handled by ε_DSC.
# dof = 9 (10 replicates, N-1)
const GROUP_DF_CP_SEG3 = CorrelatedGroup("DF_cp_seg3", [6, 7],
	[0.009136799235573225 * 1e3, -2.1393044552435705 * 1e3],
	[2.7703324106603792e-08 * 1e6    -9.7846678915291e-06 * 1e6;
	 -9.7846678915291e-06 * 1e6       0.0038340012469450084 * 1e6])

# --- Group 6: Independent char cp (dedicated charring experiments) ---
# Source: gsa_parameter_config_douglas_fir.json
# NOT the pseudo char Cp from the optimization — independently measured.
# JSON units: kJ/(kg·K) → script units: J/(kg·K), mean ×1e3, covariance ×1e6
# NOTE: Type A only. Type B handled by ε_DSC (same DSC instrument).
# dof = 4 (5 replicates, N-1)
const GROUP_CHAR_CP = CorrelatedGroup("Char_cp", [8, 9],
	[0.0018372509265710314 * 1e3, 0.3124273785167445 * 1e3],
	[0.00040051156523772673^2 * 1e6    -9.472629865704193e-05 * 1e6;
	 -9.472629865704193e-05 * 1e6       0.2369805093355783^2 * 1e6])

# --- Group 7: Heat of sorption (CR2 covariance, conditional on f_cal exact) ---
# Source: papers/energetics/figures/hos_fit_results.json → covariance_matrix_cr2
# Parameters: a=asymptote [kJ/kg], b=amplitude [kJ/kg], c=decay [1/%]
#
# Now uses CR2 (conditional on f_cal exact), with the calibration uncertainty
# handled by the separate ε_HoS group (Group 3). The old code used the total
# covariance which baked in calibration.
#
# NOTE: ε_DSC does NOT apply to HoS. HoS calibration is independent of the
# standard DSC multi-point enthalpy calibration; it uses the asymptote
# constraint referencing L_vap(T_iso) from the NIST-Webbook Antoine equation.
# dof = 19 (22 clusters - 3 params)
const GROUP_HOS = CorrelatedGroup("HoS", [10, 11, 12],
	[2450.420909987441, 875.6441221735336, 0.1121684273316275],
	[255.54778632642243     826.4091656069307     0.24363365219954056;
	 826.4091656069309     4539.412144899307      1.054987399933727;
	 0.24363365219954064    1.0549873999337271    0.00028091975868012163])

# --- Group 8: Heats of reaction (CONDITIONAL covariance from per-curve decomposition) ---
# Source: gsa_parameter_config_douglas_fir.json
#         → conditional covariance from per-curve empirical decomposition
#
# This covariance holds virgin cp at its nominal value, so Groups 5 (DF cp
# high-T) and 8 (ΔH) are conditionally independent. The virgin cp → ΔH coupling
# is handled explicitly by the Jacobian correction J_HRXN_CP × δ(slope,int) in
# build_and_solve.
#
# Type B (DSC calibration): ΔH_fit = (1+ε)·ΔH_true (calibration scales, does not
# cancel). Applied as (1+ε_DSC) multiplicative factor in build_and_solve.
#
# Units: kJ²/kg², converted to J²/kg² by ×1e6.
const GROUP_H_RXN = CorrelatedGroup("H_rxn", [13, 14, 15, 16, 17],
	[-19.20066288032251 * 1e3, -5.263201466602739 * 1e3, 264.42479037373835 * 1e3, -217.15169043762444 * 1e3, -0.32640885469229025 * 1e3],
	# Conditional covariance (per-curve decomposition): kJ²/kg² → J²/kg²
	[12.188429198060906    24.748705356848376   29.002904740480126    88.4320951009737     -6.604044259616362;
	 24.748705356848376   50.700350890015784    59.67848721837838    184.93304070273382   -13.71270702182141;
	 29.002904740480123   59.67848721837838     80.78981196904051    246.88080651001275    -7.889710589815444;
	 88.4320951009737    184.93304070273382     246.88080651001275    837.274115396477       0.9278934423897951;
	 -6.604044259616367  -13.712707021821425    -7.889710589815465     0.9278934423897454  41.07609754085615] .* 1e6)

const CORR_GROUPS = [
	GROUP_EPSILON_HFM,
	GROUP_EPSILON_DSC,
	GROUP_EPSILON_HOS,
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
# for the Jacobian correction to ΔH. Per-curve grand mean values in J/(kg·K).
@everywhere const NOM_CP_SLOPE = 9.136799235573225
@everywhere const NOM_CP_INT = -2139.3044552435705

# Jacobian J = ∂ΔH/∂(slope, intercept): total sensitivity of fitted heats of
# reaction to the virgin cp (high-T) regression coefficients used as the
# optimization baseline.
# Source: gsa_parameter_config_douglas_fir.json → per-curve empirical J
# Rows = [ΔH₁, ΔH₂, ΔH₃, ΔH₄, ΔH₅], Cols = [slope, intercept]
# Units: J_sim = J_kJ (no conversion needed — see derivation in config spec)
@everywhere const J_HRXN_CP = [
	-36266.85748394125   -73.02438035269826;
	-76518.48774219683   -157.57351983306967;
	-108968.04417626567  -244.32381607518548;
	-284878.45699294383  -635.1246074430231;
	27455.93513146396    7.708707967010774]

@everywhere begin

	"""
		build_and_solve(θ)

	Build a Douglas Fir simulation with the given parameter vector θ (length 17),
	solve it, and return the PyrolysisSolution.
	"""
	function build_and_solve(θ::Vector{Float64})
		# Unpack the 17 uncertain parameters
		ε_hfm, ε_dsc, ε_hos,
		df_cp_a1, df_cp_b1, df_cp_a3, df_cp_b3,
		char_cp_a, char_cp_b,
		hos_a, hos_b, hos_c,
		h_rxn1, h_rxn2, h_rxn3, h_rxn4, h_rxn5 = θ

		# Type B multiplicative factors for each source
		scale_hfm = 1.0 + ε_hfm    # HFM-derived Cp (low-T)
		scale_dsc = 1.0 + ε_dsc    # DSC-derived Cp (high-T, char) and ΔH
		scale_hos = 1.0 + ε_hos    # HoS asymptote calibration

		# --- Correct ΔH for virgin cp regression uncertainty (Type A coupling) ---
		# The two-stage optimization used nominal (slope, intercept) as fixed input for the
		# reactant-side sensible heat baseline. When the regression coefficients
		# deviate from nominal (Type A), the fitted ΔH shifts by J × δ(slope, int).
		δslope = df_cp_a3 - NOM_CP_SLOPE
		δint = df_cp_b3 - NOM_CP_INT
		h_rxn1_adj = h_rxn1 + J_HRXN_CP[1, 1] * δslope + J_HRXN_CP[1, 2] * δint
		h_rxn2_adj = h_rxn2 + J_HRXN_CP[2, 1] * δslope + J_HRXN_CP[2, 2] * δint
		h_rxn3_adj = h_rxn3 + J_HRXN_CP[3, 1] * δslope + J_HRXN_CP[3, 2] * δint
		h_rxn4_adj = h_rxn4 + J_HRXN_CP[4, 1] * δslope + J_HRXN_CP[4, 2] * δint
		h_rxn5_adj = h_rxn5 + J_HRXN_CP[5, 1] * δslope + J_HRXN_CP[5, 2] * δint

		# --- Correct ΔH for DSC calibration bias (Type B) ---
		# The optimization derived ΔH from DSC heat flow using DSC-measured cp as
		# baseline: q_meas = (1+ε)·q_true, cp_meas = (1+ε)·cp_true, so the residual
		# yields ΔH_fit = (1+ε)·ΔH_true. The calibration does NOT cancel — it scales
		# ΔH proportionally, just like cp. Both must be scaled by (1+ε_DSC).
		h_rxn1_adj *= scale_dsc
		h_rxn2_adj *= scale_dsc
		h_rxn3_adj *= scale_dsc
		h_rxn4_adj *= scale_dsc
		h_rxn5_adj *= scale_dsc

		# --- DF virgin heat capacity (piecewise linear, 3 segments) ---
		# Segment 1 (T ≤ 323 K): HFM instrument → scale by (1 + ε_HFM)
		# Segment 3 (T > 425 K): STA/DSC instrument → scale by (1 + ε_DSC)
		# Segment 2 (323–425 K): interpolated from endpoints of seg 1 & 3
		#   so it inherits both biases through the endpoint values
		function virgin_cp(T)
			if T <= 323.0
				(df_cp_b1 + df_cp_a1 * T) * scale_hfm
			elseif T <= 425.0
				cp_low  = (df_cp_b1 + df_cp_a1 * 323.0) * scale_hfm
				cp_high = (df_cp_b3 + df_cp_a3 * 425.0) * scale_dsc
				cp_low + (cp_high - cp_low) * (T - 323.0) / 102.0
			else
				(df_cp_b3 + df_cp_a3 * T) * scale_dsc
			end
		end

		df_cp = CallableProperty(virgin_cp)

		# --- Char heat capacity (linear, STA instrument → scale by ε_DSC) ---
		char_cp_func = CallableProperty(T -> (char_cp_b + char_cp_a * T) * scale_dsc)

		# --- Intermediate cp's derived via cumulative-yield interpolation ---
		# char endpoint also scaled by scale_dsc (same DSC instrument)
		int1_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[1]) * virgin_cp(T) +
										INTERP_PROGRESS[1] * (char_cp_b + char_cp_a * T) * scale_dsc)

		int2_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[2]) * virgin_cp(T) +
										INTERP_PROGRESS[2] * (char_cp_b + char_cp_a * T) * scale_dsc)

		int3_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[3]) * virgin_cp(T) +
										INTERP_PROGRESS[3] * (char_cp_b + char_cp_a * T) * scale_dsc)

		int4_cp = CallableProperty(T -> (1 - INTERP_PROGRESS[4]) * virgin_cp(T) +
										INTERP_PROGRESS[4] * (char_cp_b + char_cp_a * T) * scale_dsc)

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
		water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

		# --- Heat of sorption (parametric) ---
		# HoS CR2 covariance is conditional on f_cal exact; calibration uncertainty
		# is handled by the separate ε_HoS multiplicative factor (Group 3).
		# ε_DSC does not apply to HoS — the calibration pathway is independent
		# of the standard DSC multi-point calibration.
		hos = StateDependentProperty((T, ξ) -> begin
			MC = 100.0 * ξ[SA_MOISTURE_IDX] / SA_ρ_DRY_BASIS
			return (hos_a + hos_b * exp(-hos_c * MC)) * 1e3 * scale_hos
		end)

		# --- Gas transfer coefficients and permeability (FIXED) ---
		# λ_solid: D_free × ε/τ; λ_gas: free molecular diffusion at ~500K
		# κ: intrinsic matrix permeability for Darcy flow
		λ_wood   = 5.0e-7     # m²/s - virgin wood (ε~0.05, τ~5)
		λ_char   = 1.0e-5     # m²/s - porous char (ε~0.7, τ~2)
		λ_gas    = 3.0e-5     # m²/s - reference free diffusion (M~29, ~500K)
		κ_virgin = 1.0e-15    # m² - tight virgin wood
		κ_char   = 1.0e-12    # m² - open char pore network

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
					ε = 0.790965, α = 0.0, λ = λ_wood * 2.0, κ = κ_virgin * 3.0),
				# 3. Intermediate 2
				SolidComponent(:DF_int2,
					ρ = 419.8, c = int2_cp, k = 0.329938,
					ε = 0.786207, α = 0.0, λ = λ_wood * 5.0, κ = κ_virgin * 10.0),
				# 4. Intermediate 3
				SolidComponent(:DF_int3,
					ρ = 203.2, c = int3_cp, k = 0.644833,
					ε = 0.745694, α = 0.0, λ = λ_char * 0.3, κ = κ_char * 0.01),
				# 5. Intermediate 4
				SolidComponent(:DF_int4,
					ρ = 176.6, c = int4_cp, k = 0.777700,
					ε = 0.728600, α = 0.0, λ = λ_char * 0.7, κ = κ_char * 0.1),
				# 6. Char
				SolidComponent(:Char,
					ρ = 150.0, c = char_cp_func, k = 1.0,
					ε = 0.70, α = 0.0, λ = λ_char, κ = κ_char),
				# 7. Moisture (liquid, pseudo-density)
				LiquidComponent(:MOISTURE,
					ρ = 1e6, c = moisture_cp, k = 0.145,
					ε = 0.81, α = 0.0, κ = 0.0),
				# 8-12. Gas products (λ scaled by √(M_ref/M) from kinetic theory)
				GaseousComponent(:DF_gas1, M = 0.060, c = 1400.0, k = 0.020, λ = λ_gas * sqrt(0.029/0.060)),
				GaseousComponent(:DF_gas2, M = 0.070, c = 1400.0, k = 0.020, λ = λ_gas * sqrt(0.029/0.070)),
				GaseousComponent(:DF_gas3, M = 0.110, c = 1300.0, k = 0.018, λ = λ_gas * sqrt(0.029/0.110)),
				GaseousComponent(:DF_gas4, M = 0.070, c = 1500.0, k = 0.020, λ = λ_gas * sqrt(0.029/0.070)),
				GaseousComponent(:DF_gas5, M = 0.050, c = 1600.0, k = 0.025, λ = λ_gas * sqrt(0.029/0.050)),
				# 13. Water vapor
				GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.025, λ = λ_gas * sqrt(0.029/0.018)),
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

# --- New QoI constants ---
# First pyrolysis peak: skip the initial moisture transient (complete by ~20 s
# at 50 kW/m²). 30 s is conservative across all parameter realizations.
@everywhere const FIRST_PEAK_T_START = 30.0  # [s]
# Pyrolysis onset: 1% dry mass loss, analogous to TGA onset criterion.
@everywhere const PYROLYSIS_ONSET_THRESHOLD = 0.01
# Solid component indices (DF, int1-4, Char) — excludes moisture (7) and gases (8-13)
@everywhere const SOLID_COMPONENT_INDICES = 1:6

@everywhere """
	time_to_threshold(t, curve, fraction)

Find the time at which `curve` first reaches `fraction` of its total change,
defined as T₀ + fraction*(T_max - T₀) where T₀ = curve[1].
Returns Inf if the threshold is never reached.
"""
@everywhere function time_to_threshold(t::Vector{Float64}, curve::Vector{Float64}, fraction::Float64)
	T0 = curve[1]
	T_max = maximum(curve)
	threshold = T0 + fraction * (T_max - T0)
	idx = findfirst(x -> x >= threshold, curve)
	return idx !== nothing ? t[idx] : Inf
end

@everywhere """
	compute_dry_mass_series(sol) -> Vector{Float64}

Compute dry solid mass [kg/m²] at each output time from the concentration profiles.
Dry mass = sum of solid components (indices 1:6), excluding moisture and gases.
"""
@everywhere function compute_dry_mass_series(sol)
	nt = length(sol.times)
	A = sol.problem.mesh.cross_section_area
	dry_mass = Vector{Float64}(undef, nt)
	for i in 1:nt
		ξ = sol.concentrations[i]  # [NC × n_cells]
		n_cells = size(ξ, 2)
		m = 0.0
		for j in 1:n_cells
			v = sol.problem.mesh.cell_volumes[j]
			for k in SOLID_COMPONENT_INDICES
				m += ξ[k, j] * v
			end
		end
		dry_mass[i] = m / A
	end
	return dry_mass
end

@everywhere """
	find_first_pyrolysis_peak(t, mlr, global_peak_idx)

Find the first local MLR maximum after `FIRST_PEAK_T_START` and before the
global peak. Returns (peak_value, peak_time).
"""
@everywhere function find_first_pyrolysis_peak(t::Vector{Float64}, mlr::Vector{Float64},
	global_peak_idx::Int)
	i_start = findfirst(x -> x >= FIRST_PEAK_T_START, t)
	if i_start === nothing || i_start >= global_peak_idx
		return (NaN, NaN)
	end
	# Search for first local maximum in [i_start, global_peak_idx-1]
	for i in (i_start+1):(global_peak_idx-1)
		if mlr[i] >= mlr[i-1] && mlr[i] >= mlr[i+1]
			return (mlr[i], t[i])
		end
	end
	# Fallback: maximum in [i_start, midpoint to global peak]
	i_mid = div(i_start + global_peak_idx, 2)
	if i_mid > i_start
		local_range = i_start:i_mid
		idx = local_range[argmax(mlr[local_range])]
		return (mlr[idx], t[idx])
	end
	return (NaN, NaN)
end

@everywhere """
	find_pyrolysis_onset(t, dry_mass) -> Float64

Time at which cumulative dry mass loss first reaches `PYROLYSIS_ONSET_THRESHOLD`
(fraction of initial dry mass). Uses linear interpolation between time steps.
"""
@everywhere function find_pyrolysis_onset(t::Vector{Float64}, dry_mass::Vector{Float64})
	m0 = dry_mass[1]
	threshold_mass = m0 * (1.0 - PYROLYSIS_ONSET_THRESHOLD)
	idx = findfirst(m -> m <= threshold_mass, dry_mass)
	if idx === nothing
		return Inf
	elseif idx == 1
		return t[1]
	else
		# Linear interpolation for precise crossing
		f = (dry_mass[idx-1] - threshold_mass) / (dry_mass[idx-1] - dry_mass[idx])
		return t[idx-1] + f * (t[idx] - t[idx-1])
	end
end

@everywhere """
	find_back_temp_at_mass_loss(t, total_mass, back_temp, fraction) -> Float64

Back-surface temperature at the moment when total mass loss reaches `fraction`
of the initial total mass. Uses linear interpolation.
"""
@everywhere function find_back_temp_at_mass_loss(t::Vector{Float64}, total_mass::Vector{Float64},
	back_temp::Vector{Float64}, fraction::Float64)
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

@everywhere struct SensitivityResult
	# Scalars
	peak_mlr::Float64
	time_to_peak_mlr::Float64
	time_to_ignition::Float64
	total_mass_loss_pct::Float64
	t_surf_70pct::Float64
	t_surf_95pct::Float64
	t_back_5pct::Float64
	t_back_95pct::Float64
	final_thickness::Float64
	peak_mlr_first::Float64
	time_to_peak_mlr_first::Float64
	pyrolysis_onset_time::Float64
	T_back_50pct_mass_loss::Float64
	retcode::Symbol
	solve_time::Float64
	# Curves
	times::Vector{Float64}
	mlr_curve::Vector{Float64}
	surface_temp_curve::Vector{Float64}
	back_temp_curve::Vector{Float64}
	dry_mass_curve::Vector{Float64}
	total_mass_curve::Vector{Float64}
end

const SCALAR_QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr", "time_to_ignition",
	"total_mass_loss_pct", "t_surf_70pct", "t_surf_95pct",
	"t_back_5pct", "t_back_95pct",
	"final_thickness",
	"peak_mlr_first", "time_to_peak_mlr_first",
	"pyrolysis_onset_time", "T_back_50pct_mass_loss",
	"solve_time",
]

@everywhere function extract_qois(sol)
	t = sol.times
	mlr = sol.surface_mass_flux
	T_surf = sol.surface_temperature
	T_back = sol.back_temperature
	total_mass = sol.total_mass

	peak_mlr_val = maximum(mlr)
	peak_idx = argmax(mlr)
	time_to_peak = t[peak_idx]

	tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
	time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf

	total_mass_loss = (1.0 - total_mass[end] / total_mass[1]) * 100.0

	# Threshold-crossing times for temperature curves
	t_surf_70 = time_to_threshold(t, T_surf, 0.70)
	t_surf_95 = time_to_threshold(t, T_surf, 0.95)
	t_back_5  = time_to_threshold(t, T_back, 0.05)
	t_back_95 = time_to_threshold(t, T_back, 0.95)

	# --- New QoIs ---
	# 1. First pyrolysis peak MLR (skip moisture transient)
	peak_first, time_peak_first = find_first_pyrolysis_peak(t, mlr, peak_idx)

	# 2. Pyrolysis onset (1% dry mass loss)
	dry_mass = compute_dry_mass_series(sol)
	onset_time = find_pyrolysis_onset(t, dry_mass)

	# 3. Back temperature at 50% total mass loss
	T_back_50 = find_back_temp_at_mass_loss(t, total_mass, T_back, 0.50)

	return SensitivityResult(
		peak_mlr_val, time_to_peak, time_to_ignition,
		total_mass_loss,
		t_surf_70, t_surf_95,
		t_back_5, t_back_95,
		sol.thickness_history[end],
		peak_first, time_peak_first,
		onset_time, T_back_50,
		sol.retcode, sol.solve_time,
		t, mlr, T_surf, T_back, dry_mass, total_mass,
	)
end

function scalar_qoi_vector(r::SensitivityResult)
	return [r.peak_mlr, r.time_to_peak_mlr, r.time_to_ignition,
		r.total_mass_loss_pct, r.t_surf_70pct, r.t_surf_95pct,
		r.t_back_5pct, r.t_back_95pct,
		r.final_thickness,
		r.peak_mlr_first, r.time_to_peak_mlr_first,
		r.pyrolysis_onset_time, r.T_back_50pct_mass_loss,
		r.solve_time]
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

function save_envelope_csv(filepath, times, p5, p25, p50, p75, p95)
	open(filepath, "w") do io
		println(io, "time,p5,p25,p50,p75,p95")
		for i in eachindex(times)
			println(io, "$(times[i]),$(p5[i]),$(p25[i]),$(p50[i]),$(p75[i]),$(p95[i])")
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
	mkpath(joinpath(RESULTS_DIR, "dry_mass_curves"))
	mkpath(joinpath(RESULTS_DIR, "total_mass_curves"))

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
		dry_file = joinpath(RESULTS_DIR, "dry_mass_curves", @sprintf("run_%04d.csv", i))
		tm_file = joinpath(RESULTS_DIR, "total_mass_curves", @sprintf("run_%04d.csv", i))
		if isfile(mlr_file) && isfile(surf_file) && isfile(back_file) && isfile(dry_file) && isfile(tm_file)
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
					dry_data = readdlm(joinpath(RESULTS_DIR, "dry_mass_curves", @sprintf("run_%04d.csv", id)), ','; header = true)[1]
					tm_data = readdlm(joinpath(RESULTS_DIR, "total_mass_curves", @sprintf("run_%04d.csv", id)), ','; header = true)[1]

					t = Float64.(mlr_data[:, 1])
					mlr = Float64.(mlr_data[:, 2])
					T_surf = Float64.(surf_data[:, 2])
					T_back = Float64.(back_data[:, 2])
					dry_mass = Float64.(dry_data[:, 2])
					total_mass = Float64.(tm_data[:, 2])

					peak_mlr_val = maximum(mlr)
					peak_idx = argmax(mlr)
					time_to_peak = t[peak_idx]
					tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
					time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf
					total_mass_loss = (1.0 - total_mass[end] / total_mass[1]) * 100.0
					t_surf_70 = time_to_threshold(t, T_surf, 0.70)
					t_surf_95 = time_to_threshold(t, T_surf, 0.95)
					t_back_5 = time_to_threshold(t, T_back, 0.05)
					t_back_95 = time_to_threshold(t, T_back, 0.95)

					# New QoIs
					peak_first, time_peak_first = find_first_pyrolysis_peak(t, mlr, peak_idx)
					onset_time = find_pyrolysis_onset(t, dry_mass)
					T_back_50 = find_back_temp_at_mass_loss(t, total_mass, T_back, 0.50)

					final_thickness = NaN  # not recoverable from curves
					solve_time = NaN

					qois = [peak_mlr_val, time_to_peak, time_to_ignition,
						total_mass_loss, t_surf_70, t_surf_95,
						t_back_5, t_back_95, final_thickness,
						peak_first, time_peak_first,
						onset_time, T_back_50,
						solve_time]
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
				save_curve_csv(
					joinpath(RESULTS_DIR, "dry_mass_curves", @sprintf("run_%04d.csv", i)),
					result.times, result.dry_mass_curve)
				save_curve_csv(
					joinpath(RESULTS_DIR, "total_mass_curves", @sprintf("run_%04d.csv", i)),
					result.times, result.total_mass_curve)

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

		save_envelope_csv(joinpath(RESULTS_DIR, "mlr_envelope.csv"),
			common_times, p5 .* 1000, p25 .* 1000, p50 .* 1000, p75 .* 1000, p95 .* 1000)
		println("  Saved mlr_envelope.csv")

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

			save_envelope_csv(joinpath(RESULTS_DIR, "surface_temp_envelope.csv"),
				common_times, p5_ts, p25_ts, p50_ts, p75_ts, p95_ts)
			println("  Saved surface_temp_envelope.csv")
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

			save_envelope_csv(joinpath(RESULTS_DIR, "back_temp_envelope.csv"),
				common_times, p5_tb, p25_tb, p50_tb, p75_tb, p95_tb)
			println("  Saved back_temp_envelope.csv")
		end
	end

	# --- 4. Scalar QoI Histograms ---
	println("  Generating scalar QoI histograms...")
	hist_qois = ["peak_mlr", "time_to_ignition", "total_mass_loss_pct", "t_surf_70pct"]
	hist_labels = ["Peak MLR [kg/(m²·s)]", "Time to Ignition [s]",
		"Total Mass Loss [%]", "Surface Temp 70% Rise Time [s]"]
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
