"""
Group Sobol Sensitivity Indices via Polynomial Chaos Expansion

Computes variance-based sensitivity indices for the Douglas Fir cone calorimeter
simulation, decomposing output variance into contributions from each measurement
source (instrument group).

Method:
  1. Transform correlated parameter samples to independent standard normals.
	 For Normal groups: Cholesky inverse z = L⁻¹(θ - μ). For Student-t groups
	 (finite dof): Cholesky inverse then Rosenblatt sequential conditioning
	 (Rosenblatt, 1952) using the conditional-t property of the multivariate-t
	 distribution to produce exactly independent standard normals.
  2. Fit an adaptive-degree (2–4) Polynomial Chaos Expansion (PCE) surrogate to
	 each scalar QoI. Degree is selected by maximizing leave-one-out Q². For
	 degree ≥ 3 where the oversampling ratio is low, a sparse LASSO fit
	 (coordinate descent with 5-fold CV) is used instead of OLS.
  3. Compute first-order and total group Sobol indices analytically from the
	 PCE coefficients.
  4. Quantify index uncertainty via BCa (bias-corrected and accelerated)
	 bootstrap (1000 replicates, Efron 1987). BCa uses jackknife influence
	 values (efficient rank-1 OLS downdate) to correct for bias and skewness
	 in the bootstrap distribution of bounded [0,1] Sobol indices.
  5. Assess convergence by computing indices on nested subsamples, with a
	 three-criterion diagnostic (terminal stability, no oscillation, no drift).

The 16 uncertain parameters fall into 7 independent groups (correlations are
within-group only), so standard Sobol formulas apply without correlated-input
extensions.

Groups (by measurement source):
  1. HFM calibration bias (ε_HFM)           — 1 parameter
  2. DSC calibration bias (ε_DSC)           — 1 parameter
  3. Virgin Cp, low-T (HFM regression)      — 2 parameters (a₁, b₁)
  4. Virgin Cp, high-T (STA per-curve)      — 2 parameters (a₃, b₃)
  5. Char Cp (independent char STA)         — 2 parameters (a_char, b_char)
  6. Heat of sorption (total cov.)          — 3 parameters (a, b, c)
  7. Heats of reaction (conditional cov.)   — 5 parameters (ΔH₁–ΔH₅)

Reads:  results/parameters.csv, results/scalar_qois.csv
Writes: results/sobol_group_indices.csv, results/sobol_bootstrap_ci.csv,
		results/pce_validation.csv, results/pce_coefficients.csv,
		results/sobol_variance_decomposition.png,
		results/sobol_bar_<qoi>.png, results/pce_validation_<qoi>.png,
		results/sobol_convergence_<qoi>.png

Usage:  julia --threads=auto --project=. examples/sensitivity_analysis/sensitivity_sobol.jl
"""

using PolyChaos
using LinearAlgebra
using Statistics
using Distributions    # TDist, cdf, quantile, Normal — for Rosenblatt transform
using DelimitedFiles
using Printf
using Plots
using Random

const RESULTS_DIR = joinpath(@__DIR__, "results")

# ============================================================================
# Group Definitions (must match sensitivity_douglas_fir.jl exactly)
# ============================================================================

struct ParamGroup
    name::String
    indices::Vector{Int}     # indices into the 16-element θ vector
    μ::Vector{Float64}
    Σ::Matrix{Float64}
    ν::Float64               # degrees of freedom: Inf → Normal, finite → Student-t
end

const GROUPS = [
    # 1. HFM calibration bias (data-driven, SRM 1453/1453e, DF thickness)
    ParamGroup("HFM calibration", [1], [0.0], [0.0046154243860687^2;;], Inf),

    # 2. DSC calibration bias (certified-only RMS, 1.13%, Type B, ν=∞)
    ParamGroup("DSC calibration", [2], [0.0], [0.011322211895623097^2;;], Inf),

    # 3. Virgin Cp, low-T (HFM specimen-level regression, slope/intercept in J/(kg·K))
    # Source: gsa_parameter_config_douglas_fir.json (12 runs → 6 specimens DF_S07…DF_S12)
    # dof = n_specimens - 1 = 5
    ParamGroup(
        "Virgin Cp, low-T",
        [3, 4],
        [0.004160854771848128 * 1e3, -0.02289659782621535 * 1e3],
        [
            8.075333870288374e-05^2 * 1e6 -1.1398271736701971e-06 * 1e6;
            -1.1398271736701971e-06 * 1e6 0.01705528336860232^2 * 1e6
        ],
        5.0,
    ),

    # 4. Virgin Cp, high-T (STA per-curve, slope/intercept in J/(kg·K))
    # Source: gsa_parameter_config_douglas_fir.json
    ParamGroup(
        "Virgin Cp, high-T",
        [5, 6],
        [0.009136799235573225 * 1e3, -2.1393044552435705 * 1e3],
        [
            2.7703324106603792e-08 * 1e6 -9.7846678915291e-06 * 1e6;
            -9.7846678915291e-06 * 1e6 0.0038340012469450084 * 1e6
        ],
        9.0,
    ),

    # 5. Char Cp (independent STA, slope/intercept in J/(kg·K))
    # Source: douglas_fir_independent_char_sta/regression_results.json
    ParamGroup(
        "Char Cp",
        [7, 8],
        [0.0018372509265710314 * 1e3, 0.3124273785167445 * 1e3],
        [
            0.00040051156523772673^2 * 1e6 -9.472629865704193e-05 * 1e6;
            -9.472629865704193e-05 * 1e6 0.2369805093355783^2 * 1e6
        ],
        4.0,
    ),

    # 6. Heat of sorption (total covariance (CR2 + calibration))
    # Source: papers/energetics/figures/hos_fit_results.json -> covariance_matrix
    ParamGroup(
        "Heat of sorption",
        [9, 10, 11],
        [2450.420909987441, 875.6441221735336, 0.1121684273316275],
        [
            515.0312055400392 919.1341869013884 0.24363365219954056;
            919.1341869013885 4572.546938637243 1.054987399933727;
            0.24363365219954064 1.0549873999337271 0.00028091975868012163
        ],
        19.0,
    ),

    # 7. Heats of reaction (conditional covariance, kJ²/kg² -> J²/kg²)
    # Source: gsa_parameter_config_douglas_fir.json -> conditional_covariance
    ParamGroup(
        "Heats of reaction",
        [12, 13, 14, 15, 16],
        [
            -19.20066288032251 * 1e3,
            -5.263201466602739 * 1e3,
            264.42479037373835 * 1e3,
            -217.15169043762444 * 1e3,
            -0.32640885469229025 * 1e3,
        ],
        [
            12.188429198060906 24.748705356848376 29.002904740480126 88.4320951009737 -6.604044259616362;
            24.748705356848376 50.700350890015784 59.67848721837838 184.93304070273382 -13.71270702182141;
            29.002904740480123 59.67848721837838 80.78981196904051 246.88080651001275 -7.889710589815444;
            88.4320951009737 184.93304070273382 246.88080651001275 837.274115396477 0.9278934423897951;
            -6.604044259616367 -13.712707021821425 -7.889710589815465 0.9278934423897454 41.07609754085615
        ] .* 1e6,
        9.0,
    ),
]

const N_PARAMS = 16

# Scalar QoIs to analyze (skip solve_time — not physical,
# time_to_peak_mlr_first — discontinuous peak switching gives Q²≈0, use SRRC instead)
const QOI_NAMES = [
    "peak_mlr",
    "time_to_peak_mlr",
    "time_to_ignition",
    "total_mass_loss_pct",
    "t_surf_70pct",
    "t_surf_95pct",
    "t_back_5pct",
    "t_back_95pct",
    "final_thickness",
    "peak_mlr_first",
    "pyrolysis_onset_time",
    "T_back_50pct_mass_loss",
    "time_to_comb_ignition",
    "time_to_comb_extinction",
    "avg_mlr_comb_burning",
    "peak_surface_temp",
]

# Column indices in scalar_qois.csv (after sample_id and retcode columns)
const ALL_QOI_NAMES = [
    "peak_mlr",
    "time_to_peak_mlr",
    "time_to_ignition",
    "total_mass_loss_pct",
    "t_surf_70pct",
    "t_surf_95pct",
    "t_back_5pct",
    "t_back_95pct",
    "final_thickness",
    "peak_mlr_first",
    "time_to_peak_mlr_first",
    "pyrolysis_onset_time",
    "T_back_50pct_mass_loss",
    "time_to_comb_ignition",
    "time_to_comb_extinction",
    "avg_mlr_comb_burning",
    "peak_surface_temp",
    "solve_time",
]

# PCE configuration
#
# PCE_DEGREE_RANGE is PINNED to 3:3 based on the convergence analysis in
# `pce_degree_selection.md` (dated 2026-04-10, dataset n = 8192). All 16
# scalar QoIs have group Sobol indices that have converged at degree 3 to
# within 0.01 per group (the pre-registered tolerance), and raising to
# degree 4 changes no Sobol index by more than 0.0055 in any group while
# dropping oversampling from 8.45 to 1.69 and forcing sparse LASSO fits.
#
# DO NOT change this to a wider range without first re-running
# `compare_degrees.jl` on the current data. If compare_degrees still
# reports "USE DEGREE 3" then the pin is correct; if it reports a
# different degree, update both this constant AND the markdown decision
# record, and document the reason in the paper.
#
# To temporarily explore other degrees (e.g., for a sensitivity check in
# a supplementary analysis), change this to the desired UnitRange. The
# rest of the code handles arbitrary ranges unchanged.
const PCE_DEGREE_RANGE = 3:3      # see pce_degree_selection.md for rationale
const SPARSE_THRESHOLD = 3.0      # Use sparse (LASSO) fit when oversampling ratio < this
const N_BOOTSTRAP = 1000          # Bootstrap replicates for Sobol CI
const BOOTSTRAP_CI_LEVEL = 0.95   # Confidence level for bootstrap intervals
# Subsample sizes for convergence diagnostics. Each fraction must be a power
# of 1/2 (so the corresponding prefix 2^k of the Sobol sequence is a balanced
# (t, m, s)-net). At n = 8192 = 2^13 these correspond to [512, 1024, 2048,
# 4096, 8192] = [2^9, 2^10, 2^11, 2^12, 2^13]. Non-power-of-2 prefixes lose
# the low-discrepancy structure and give misleading convergence plots.
const CONVERGENCE_FRACTIONS = [0.0625, 0.125, 0.25, 0.5, 1.0]

# Fraction-of-variance threshold for the pre-bootstrap basis-contamination
# diagnostic. If verify_basis_normalization flags some basis functions as
# failing the 3σ orthonormality check AND more than this fraction of any
# QoI's PCE variance D = Σ_{α≠0} c_α² comes from those failing terms, the
# main loop aborts BEFORE the (hours-long) bootstrap phase with a detailed
# per-QoI report. The threshold is deliberately conservative: 5% means a
# Sobol index in that QoI could be biased by up to ~5% relative, which is
# typically within published CI widths but worth flagging.
#
# Override by setting `JULIA_SOBOL_FORCE_BIASED_RUN=1` in the environment
# to skip the abort (e.g., when you know the affected basis functions are
# not in the Sobol-index numerator for the groups you care about, or when
# running deliberately to measure the bias).
const BIASED_BASIS_ABORT_THRESHOLD = 0.05

# ============================================================================
# Step 1: Transform to independent standard normals
# ============================================================================

"""
	transform_to_standard_normal(param_data, groups)

Transform correlated parameter samples to independent standard normal space.

For Normal groups (ν=Inf): Cholesky inverse z = L⁻¹(θ - μ) gives standard
normal variates directly (components are independent by construction).

For Student-t groups (finite ν): the Cholesky inverse gives
  z_hat = √((ν-2)/V)·z  where z ~ N(0,I) and V ~ χ²(ν)/ν
which has a multivariate-t(ν, 0, I) distribution after rescaling by √(ν/(ν-2)).
Components share the latent mixing variable V, so they are marginally t(ν) but
NOT independent. A marginal (Nataf) transform would fix each marginal to N(0,1)
but leave the t-copula dependence intact, breaking Hermite polynomial
orthogonality.

Instead, we apply the Rosenblatt sequential conditioning transform
(Rosenblatt, 1952), which exploits the conditional-t property of the
multivariate-t distribution. For y ~ t_p(ν, 0, I):
  - y₁ ~ t(ν)
  - y_k | y_{1:k-1} ~ t(ν+k-1, 0, (ν + d²_{k-1})/(ν + k - 1))
where d²_{k-1} = Σⱼ₌₁ᵏ⁻¹ yⱼ² is the accumulated squared Mahalanobis distance.
Mapping each conditional CDF through Φ⁻¹ produces exactly independent standard
normals, because the sequential conditioning absorbs the shared mixing variable
through d².

For p=1 groups, the Rosenblatt reduces to the marginal (Nataf) transform.

References:
  Rosenblatt (1952), Ann. Math. Stat. 23(3), 470–472.
  Kotz & Nadarajah (2004), Multivariate t Distributions, Ch. 2.

Input:  param_data [n_samples × 16] in physical units
Output: Z [n_samples × 16] in standard normal space
"""
function transform_to_standard_normal(
    param_data::Matrix{Float64},
    groups::Vector{ParamGroup},
)
    n = size(param_data, 1)
    Z = similar(param_data)
    std_normal = Normal()

    for g in groups
        idx = g.indices
        L = cholesky(Symmetric(g.Σ)).L

        for i = 1:n
            z_hat = L \ (param_data[i, idx] - g.μ)

            if isfinite(g.ν)
                # Rosenblatt transform for multivariate-t(ν, 0, I)
                # After Cholesky whitening, z_hat ~ scaled t_p(ν, 0, (ν-2)/ν · I).
                # Rescale to standard multivariate-t: y = z_hat · √(ν/(ν-2)).
                ν = g.ν
                scale = sqrt(ν / (ν - 2))
                d_sq = 0.0  # accumulated squared Mahalanobis distance

                for j in eachindex(z_hat)
                    y_j = z_hat[j] * scale

                    # Conditional: y_k | y_{1:k-1} ~ t(ν+k-1) with scale √((ν+d²)/(ν+k-1))
                    ν_cond = ν + j - 1
                    s_cond = j == 1 ? 1.0 : sqrt((ν + d_sq) / ν_cond)

                    u = cdf(TDist(ν_cond), y_j / s_cond)
                    u = clamp(u, 1e-15, 1.0 - 1e-15)
                    z_hat[j] = quantile(std_normal, u)

                    d_sq += y_j^2
                end
            end

            Z[i, idx] = z_hat
        end
    end

    return Z
end

# ============================================================================
# Step 2: Build PCE design matrix and fit coefficients
# ============================================================================

"""
	build_design_matrix(Z, mop)

Evaluate all multivariate Hermite basis polynomials at all sample points.
Returns the orthonormalized design matrix Φ [n_samples × n_basis].

PolyChaos uses monic Hermite polynomials where E[He_n²] = n!.
We normalize to E[Ψ²] = 1 by dividing by √(∏ αᵢ!).

Uses cached univariate evaluations: precompute He_k(z_d) for each dimension d
and degree k once, then assemble multivariate basis via products. This avoids
calling PolyChaos.evaluate per-basis-function (which is O(n_basis) calls) and
gives O(n_samples × n_dims × degree + n_samples × n_basis) total work.
"""
function build_design_matrix(Z::Matrix{Float64}, mop::MultiOrthoPoly)
    n_samples = size(Z, 1)
    n_dims = size(Z, 2)
    n_basis = mop.dim
    deg = mop.deg  # maximum total degree

    # --- Precompute univariate monic Hermite values ---
    # uni_vals[d][k+1][i] = He_k(Z[i,d]) for dimension d, degree k, sample i
    # PolyChaos stores univariate polys in mop.uni; evaluate each degree 0..deg
    uni_vals = Vector{Vector{Vector{Float64}}}(undef, n_dims)
    for d = 1:n_dims
        uni_vals[d] = Vector{Vector{Float64}}(undef, deg + 1)
        for k = 0:deg
            # Evaluate degree-k univariate orthogonal poly at all Z[:,d]
            uni_vals[d][k+1] = PolyChaos.evaluate(k, Z[:, d], mop.uni[d])
        end
    end

    # --- Assemble multivariate basis via tensor products ---
    Φ = Matrix{Float64}(undef, n_samples, n_basis)
    for j = 1:n_basis
        α = mop.ind[j, :]  # multi-index: degree per dimension

        # Product of univariate evaluations
        col = ones(n_samples)
        for d = 1:n_dims
            k = α[d]
            if k > 0
                col .*= uni_vals[d][k+1]
            end
        end

        # Normalize to orthonormal: divide by √(∏ αᵢ!)
        norm_factor = sqrt(prod(factorial.(α)))
        if norm_factor > 0
            col ./= norm_factor
        end

        Φ[:, j] = col
    end

    return Φ
end

"""
	verify_basis_normalization(Φ, mop; z_threshold=3.0)

Empirically verify that the orthonormalized PCE basis satisfies E[Ψ_j²] ≈ 1
(using sample mean as proxy). This catches systematic normalization errors
(e.g., wrong norm factor, broken isoprobabilistic transform) AND finite-sample
QMC effects that bias the sample second moment for high-variation basis
functions at the current sample size.

Uses a z-score test rather than a fixed tolerance: each basis function's
deviation |E[Ψ²]-1| is compared against its empirical standard error
SE = std(Ψ²)/√n. This accounts for the fact that higher-order basis functions
have larger sampling variance (e.g., for a product of k degree-1 Hermites,
Var(Ψ²) = 3ᵏ - 1), so a fixed tolerance would produce false positives at
high PCE degrees or modest sample sizes.

Returns `(max_z, n_fail, pass, failing_idx, empirical_norms)`:
  - `max_z`           = max |E[Ψ²] − 1| / SE across non-constant basis functions
  - `n_fail`          = number of non-constant basis functions with z > z_threshold
  - `pass`            = true iff n_fail == 0
  - `failing_idx`     = list of basis-function column indices (into Φ/c) that
                        failed the test; used by the pre-bootstrap diagnostic to
                        determine whether failing terms are actually contributing
                        meaningfully to the PCE variance decomposition
  - `empirical_norms` = full vector of E[Ψ²_j] for all j (for reporting)
"""
function verify_basis_normalization(
    Φ::Matrix{Float64},
    mop::MultiOrthoPoly;
    z_threshold::Float64 = 3.0,
)
    n_samples, n_basis = size(Φ)
    Φ_sq = Φ .^ 2
    empirical_norms = vec(mean(Φ_sq, dims = 1))
    empirical_se = vec(std(Φ_sq, dims = 1)) ./ sqrt(n_samples)

    max_z = 0.0
    n_fail = 0
    worst_idx = 1
    failing_idx = Int[]
    for j = 2:n_basis  # skip constant term (always = 1)
        se_j = empirical_se[j]
        z_j = se_j > 0 ? abs(empirical_norms[j] - 1.0) / se_j : 0.0
        if z_j > max_z
            max_z = z_j
            worst_idx = j
        end
        if z_j > z_threshold
            n_fail += 1
            push!(failing_idx, j)
        end
    end

    pass = n_fail == 0
    if !pass
        α_worst = mop.ind[worst_idx, :]
        dev_worst = empirical_norms[worst_idx] - 1.0
        se_worst = empirical_se[worst_idx]
        @warn "Basis normalization check FAILED: $n_fail / $(n_basis-1) basis functions " *
              "exceed $(z_threshold)σ threshold. " *
              "Worst: basis $worst_idx (α=$α_worst), E[Ψ²]=$(round(empirical_norms[worst_idx], digits=4)) " *
              "(z=$(round(max_z, digits=2)), SE=$(round(se_worst, digits=4))). " *
              "Sobol indices may be biased — run the pre-bootstrap diagnostic " *
              "(emitted automatically before Step 4) to quantify the impact per QoI."
    end

    return max_z, n_fail, pass, failing_idx, empirical_norms
end

"""
	compute_ols_hat_diag(Φ)

Compute the diagonal of the OLS hat matrix H = Φ(ΦᵀΦ)⁻¹Φᵀ without
materializing the full n × n matrix.

For n = 8192, p = 969 (degree 3), the full H is ~512 MB per QoI. Since we
only need the diagonal and the p × n factor G = (ΦᵀΦ)⁻¹Φᵀ (which is already
used for BCa jackknife), we can compute h_ii = ⟨Φ[i, :], G[:, i]⟩ directly
from G at ~61 MB total memory.

Returns `(h_diag::Vector{Float64}, G::Matrix{Float64})`. The caller is
expected to pass G back into BCa jackknife to avoid a second (ΦᵀΦ)⁻¹Φᵀ
factorization.
"""
function compute_ols_hat_diag(Φ::Matrix{Float64})
    n, p = size(Φ)
    # G = (Φ'Φ)^{-1} Φ'  — p × n factor, also used by BCa jackknife.
    # Using `\` with a Symmetric wrapper lets Julia pick a Cholesky-based
    # solver when Φ'Φ is positive definite (the usual case for a
    # full-column-rank design matrix).
    G = Symmetric(Φ' * Φ) \ Φ'
    h_diag = Vector{Float64}(undef, n)
    @inbounds for i = 1:n
        s = 0.0
        for j = 1:p
            s += Φ[i, j] * G[j, i]
        end
        h_diag[i] = s
    end
    return h_diag, G
end

"""
	fit_pce(Φ, y)

Fit PCE coefficients via ordinary least squares.
Returns (c, R², Q², y_pred, h_diag, residuals). The hat matrix diagonal
and residuals are needed downstream for BCa bootstrap acceleration constants
via the efficient rank-1 jackknife downdate.
"""
function fit_pce(Φ::Matrix{Float64}, y::Vector{Float64})
    n, p = size(Φ)

    # Input validation
    if !all(isfinite, Φ)
        error("fit_pce: Φ contains $(count(!isfinite, Φ)) non-finite entries")
    end
    if !all(isfinite, y)
        error("fit_pce: y contains non-finite entries")
    end

    # OLS fit
    c = Φ \ y
    y_pred = Φ * c
    residuals = y - y_pred

    # R²
    ss_res = sum(residuals .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)

    # Leave-one-out Q² via analytical formula (no refitting needed):
    #   Q² = 1 - (1/n) Σ [eᵢ/(1-hᵢᵢ)]² / Var(y)
    # h_diag is computed via compute_ols_hat_diag, which uses the p × n
    # factor G = (ΦᵀΦ)⁻¹Φᵀ instead of materializing the n × n hat matrix
    # (saving ~512 MB at n = 8192).
    h_diag, _ = compute_ols_hat_diag(Φ)
    # Guard against hᵢᵢ ≈ 1 (leverage-1 points, degenerate when p = n)
    safe_denom = 1.0 .- h_diag
    any_bad_leverage = any(abs.(safe_denom) .< eps(Float64))
    loo_residuals = residuals ./ safe_denom
    ss_loo = sum(loo_residuals .^ 2)

    R2 = isfinite(ss_res) && ss_tot > 0 ? 1.0 - ss_res / ss_tot : -Inf
    Q2 = isfinite(ss_loo) && !any_bad_leverage && ss_tot > 0 ? 1.0 - ss_loo / ss_tot : -Inf

    return c, R2, Q2, y_pred, h_diag, residuals
end

# ============================================================================
# Step 3: Compute group Sobol indices from PCE coefficients
# ============================================================================

"""
	compute_group_sobol(c, mop, groups)

Compute first-order (S_g) and total (S_Tg) group Sobol indices from PCE
coefficients. Uses the orthonormality of the basis: D = Σ_{α≠0} c_α².
"""
function compute_group_sobol(
    c::Vector{Float64},
    mop::MultiOrthoPoly,
    groups::Vector{ParamGroup},
)

    n_basis = mop.dim
    n_groups = length(groups)

    # Build group membership lookup: param_index → group_index
    param_to_group = zeros(Int, N_PARAMS)
    for (gi, g) in enumerate(groups)
        for idx in g.indices
            param_to_group[idx] = gi
        end
    end

    # Total variance (exclude constant term α=0)
    D = sum(c[j]^2 for j = 2:n_basis)

    if D ≈ 0
        return zeros(n_groups), zeros(n_groups), 0.0
    end

    S_first = zeros(n_groups)   # first-order group indices
    S_total = zeros(n_groups)   # total group indices

    for j = 2:n_basis  # skip constant term
        α = mop.ind[j, :]
        c2 = c[j]^2

        # Which groups does this basis function "touch"?
        active_groups = Set{Int}()
        for (pi, deg) in enumerate(α)
            if deg > 0
                push!(active_groups, param_to_group[pi])
            end
        end

        # Total index: any term involving group g
        for gi in active_groups
            S_total[gi] += c2
        end

        # First-order index: terms involving ONLY one group
        if length(active_groups) == 1
            gi = first(active_groups)
            S_first[gi] += c2
        end
    end

    S_first ./= D
    S_total ./= D

    return S_first, S_total, D
end

# ============================================================================
# Step 3b: Sparse PCE via LASSO (coordinate descent with K-fold CV)
# ============================================================================

"""
	soft_threshold(z, λ)

Soft-thresholding operator for LASSO: sign(z) * max(|z| - λ, 0).
"""
soft_threshold(z::Float64, λ::Float64) = sign(z) * max(abs(z) - λ, 0.0)

# ----------------------------------------------------------------------------
# Global LASSO progress registry
#
# Each active fit_pce_lasso call registers its atomic progress counters here.
# An external monitor (e.g., in the degree-selection loop) can poll the
# registry to see what every concurrent LASSO call is doing, without needing
# a dedicated scheduler thread per call. This avoids the thread-starvation
# problem that would occur if each LASSO call spawned its own sleep() task.
# ----------------------------------------------------------------------------

struct LassoProgress
    id::Int
    lambda_done::Vector{Threads.Atomic{Int}}   # per fold
    outer_iters::Vector{Threads.Atomic{Int}}   # per fold
    n_lambda::Int
    p::Int
    t_start::Float64
end

const LASSO_PROGRESS_REGISTRY = LassoProgress[]
const LASSO_PROGRESS_LOCK = ReentrantLock()
const LASSO_PROGRESS_COUNTER = Ref(0)

function lasso_register!(lambda_done, outer_iters, n_lambda, p, t_start)
    lock(LASSO_PROGRESS_LOCK) do
        LASSO_PROGRESS_COUNTER[] += 1
        id = LASSO_PROGRESS_COUNTER[]
        push!(
            LASSO_PROGRESS_REGISTRY,
            LassoProgress(id, lambda_done, outer_iters, n_lambda, p, t_start),
        )
        id
    end
end

function lasso_unregister!(id::Int)
    lock(LASSO_PROGRESS_LOCK) do
        idx = findfirst(p -> p.id == id, LASSO_PROGRESS_REGISTRY)
        if idx !== nothing
            deleteat!(LASSO_PROGRESS_REGISTRY, idx)
        end
    end
end

function lasso_progress_snapshot()
    lock(LASSO_PROGRESS_LOCK) do
        # Return a copy so the caller doesn't need to hold the lock while printing
        [
            (
                p.id,
                p.p,
                p.n_lambda,
                p.t_start,
                [x[] for x in p.lambda_done],
                [x[] for x in p.outer_iters],
            ) for p in LASSO_PROGRESS_REGISTRY
        ]
    end
end

function format_elapsed(s::Float64)
    # floor(Int, ...) handles non-integer Float64 safely, unlike Int(...) which
    # throws InexactError on fractional values. s % 60 returns a Float64
    # remainder (e.g., 5.3), not an integer, so Int(s % 60) explodes — which
    # is exactly the bug that killed the monitor in previous runs.
    if s < 60
        return @sprintf("%.0fs", s)
    elseif s < 3600
        return @sprintf("%dm%02ds", floor(Int, s / 60), floor(Int, s % 60))
    else
        return @sprintf("%dh%02dm", floor(Int, s / 3600), floor(Int, (s % 3600) / 60))
    end
end

"""
	print_monitor_status(deg, t0, now, qois_done_atomic, valid_qois, prev_state)

Print a single monitor status block. Extracted from the monitor task body so
it can be wrapped in a targeted `try/catch` (a print failure shouldn't kill
the whole monitor loop).

`prev_state::Dict{Int, NamedTuple{...}}` tracks per-LASSO-call history
between cycles for stuck-call detection. Pass an empty dict on the first call;
the function mutates it in place. A call is considered STUCK when its slowest
fold's λ index hasn't advanced in `stuck_threshold_secs` seconds while its
cumulative iteration count is still climbing.
"""
function print_monitor_status(
    deg::Int,
    t0::Float64,
    now::Float64,
    qois_done_atomic::Threads.Atomic{Int},
    valid_qois::Vector{String},
    prev_state::Dict{
        Int,
        NamedTuple{(:min_lam, :sum_iter, :last_change_time),Tuple{Int,Int,Float64}},
    };
    stuck_threshold_secs::Float64 = 30.0,
)
    elapsed = now - t0
    n_done = qois_done_atomic[]
    snaps = lasso_progress_snapshot()
    n_active = length(snaps)
    println(
        @sprintf(
            "  [deg %d monitor %s] QoIs: %d/%d done, %d LASSO calls active",
            deg,
            format_elapsed(elapsed),
            n_done,
            length(valid_qois),
            n_active
        )
    )
    for (id, p_lasso, nlam, t_lasso, lam_done, iters) in snaps
        # Defensive: skip malformed entries (shouldn't happen but be safe)
        if isempty(lam_done) || isempty(iters)
            println(@sprintf("    ↳ lasso#%d (p=%d): <no fold data>", id, p_lasso))
            continue
        end
        min_lam, max_lam = extrema(lam_done)
        sum_iter = sum(iters)
        call_elapsed = now - t_lasso

        # Stuck-call detection: compare against previous observation.
        #
        # Two failure modes:
        #   ⚠ STUCK — slowest fold's λ-index hasn't advanced but iter count IS
        #             climbing (inner loop thrashing on a hard λ)
        #   ⚠ HANG  — neither λ-index nor iter count advancing (truly frozen,
        #             investigate: scheduler starvation, deadlock, or
        #             all threads mid-very-long inner iteration)
        stuck_marker = ""
        prev = get(prev_state, id, nothing)
        if prev === nothing
            # First observation for this id
            prev_state[id] =
                (min_lam = min_lam, sum_iter = sum_iter, last_change_time = now)
        else
            if min_lam != prev.min_lam
                # Slowest fold advanced — reset
                prev_state[id] =
                    (min_lam = min_lam, sum_iter = sum_iter, last_change_time = now)
            else
                # Slowest fold has NOT advanced since last observation
                stuck_secs = now - prev.last_change_time
                iters_added = sum_iter - prev.sum_iter
                if stuck_secs >= stuck_threshold_secs
                    if iters_added > 0
                        stuck_marker =
                            @sprintf("  ⚠ STUCK %.0fs (iters +%d)", stuck_secs, iters_added)
                    else
                        stuck_marker = @sprintf("  ⚠ HANG %.0fs (0 iters)", stuck_secs)
                    end
                end
                # Update sum_iter but NOT last_change_time (still stuck)
                prev_state[id] = (
                    min_lam = prev.min_lam,
                    sum_iter = sum_iter,
                    last_change_time = prev.last_change_time,
                )
            end
        end

        println(
            @sprintf(
                "    ↳ lasso#%d (p=%d): λ-range=[%d..%d]/%d, Σ iters=%d, t=%.0fs%s",
                id,
                p_lasso,
                min_lam,
                max_lam,
                nlam,
                sum_iter,
                call_elapsed,
                stuck_marker
            )
        )
    end
    # Prune stale entries from prev_state (LASSO calls that have finished)
    active_ids = Set{Int}(s[1] for s in snaps)
    for k in collect(keys(prev_state))
        if k ∉ active_ids
            delete!(prev_state, k)
        end
    end
    flush(stdout)
    return nothing
end

"""
	standardize_design_matrix(Φ; min_col_norm = sqrt(eps()))

Standardize the columns of a PCE design matrix for LASSO coordinate descent.

The first column (constant basis) is left unscaled — it represents the intercept
and is not penalized by LASSO. Columns 2..p are divided by their L²-norm / √n so
that each standardized column has unit sample variance (col_norm = 1.0).

**Dead column handling:** columns with L²-norm / √n below `min_col_norm` are
zeroed out rather than divided through. This prevents the numerical blow-up where
dividing Φ[:,j] by ε produces entries ~1/ε that propagate to Inf/NaN in later dot
products. A genuinely degenerate basis column (one that is essentially zero across
all samples) contributes nothing to the fit, so LASSO correctly assigns it zero
weight — zeroing it out here preserves that and avoids the singularity.

Returns a NamedTuple with:
  - `Φ_scaled`::Matrix{Float64}  — the standardized design matrix (full copy)
  - `col_norms`::Vector{Float64} — the scaling factors; col_norms[j] = 0 for
    dead columns (and for the unpenalized j=1 intercept). Used to rescale
    the LASSO coefficients back to the original Φ scale on exit.
  - `col_sq_norms`::Vector{Float64} — per-column sum of squares of `Φ_scaled`,
    precomputed for the coordinate-descent denominator.
  - `n_dead`::Int — number of columns zeroed out (for diagnostics).

The computation is deterministic and identical to the inline version previously
in `fit_pce_lasso`, so bit-for-bit output is preserved when this helper is used.
"""
function standardize_design_matrix(
    Φ::Matrix{Float64};
    min_col_norm::Float64 = sqrt(eps(Float64)),
)

    n, p = size(Φ)
    col_norms = zeros(p)
    n_dead = 0
    Φ_scaled = copy(Φ)
    for j = 2:p
        col_norms[j] = norm(Φ_scaled[:, j]) / sqrt(n)
        if col_norms[j] > min_col_norm
            Φ_scaled[:, j] ./= col_norms[j]
        else
            # Dead column: zero it out and mark for rescale skip
            Φ_scaled[:, j] .= 0.0
            col_norms[j] = 0.0
            n_dead += 1
        end
    end

    # Precompute column squared norms (denominator of coordinate-descent updates)
    col_sq_norms = vec(sum(Φ_scaled .^ 2, dims = 1))  # length p

    return (
        Φ_scaled = Φ_scaled,
        col_norms = col_norms,
        col_sq_norms = col_sq_norms,
        n_dead = n_dead,
    )
end

"""
	run_cv_fold!(fold, Φ_scaled, y_centered, col_sq_norms, fold_ids, lambdas,
	             max_iter, max_inner_iter, tol, outer_iters_used, lambda_done)

Run a single CV fold's coordinate-descent LASSO fit over the full λ path,
returning the per-λ test MSE and the fold's final coefficient vector (at
the smallest λ in the path, with warm-starting from previous λ values).

Extracted from `fit_pce_lasso` so that Julia compiles it as a standalone
function — keeping all locals properly isolated per task (the equivalent
code inlined inside a `map() do` + `@spawn` closure triggered what appeared
to be a compiler optimization bug that occasionally corrupted residuals on
non-first folds).

This function operates directly on the shared, read-only `Φ_scaled` (no
per-fold matrix copies). The residual `r_full` is maintained full-length with
zeros on test rows, so contiguous dot products `dot(Φ_scaled[:, j], r_full)`
equal the training-only dot products used in the previous implementation.
Residual updates iterate `train_idx` directly.

Returns a tuple `(mse_row::Vector{Float64}, c_fold::Vector{Float64})`:
  - `mse_row[li]` = per-sample test MSE at the li-th λ in `lambdas`
  - `c_fold` = the fold's coefficient vector at the final (smallest) λ,
               useful for KKT correctness testing against the LASSO loss
               on (Φ_scaled[train_idx, :], y_centered[train_idx]).
"""
function run_cv_fold!(
    fold::Int,
    Φ_scaled::Matrix{Float64},
    y_centered::Vector{Float64},
    col_sq_norms::Vector{Float64},
    fold_ids::Vector{Int},
    lambdas::Vector{Float64},
    max_iter::Int,
    max_inner_iter::Int,
    tol::Float64,
    outer_iters_used::Vector{Threads.Atomic{Int}},
    lambda_done::Vector{Threads.Atomic{Int}},
)
    n, p = size(Φ_scaled)
    n_lambda = length(lambdas)

    # Fold index vectors
    train_idx = Int[]
    test_idx = Int[]
    sizehint!(train_idx, n)
    sizehint!(test_idx, n)
    @inbounds for i = 1:n
        if fold_ids[i] == fold
            push!(test_idx, i)
        else
            push!(train_idx, i)
        end
    end
    n_train = length(train_idx)
    n_test = length(test_idx)

    # col_sq_train[j] = sum over training rows of Φ_scaled[i, j]².
    # Compute via col_sq_norms - (sum over test rows of Φ_scaled[i,j]²).
    col_sq_train = Vector{Float64}(undef, p)
    @inbounds for j = 1:p
        Φj = view(Φ_scaled, :, j)
        s_test = 0.0
        for ii = 1:n_test
            i = test_idx[ii]
            s_test += Φj[i]^2
        end
        col_sq_train[j] = col_sq_norms[j] - s_test
    end

    # Initialize: r_full = y_centered with test rows zeroed.
    c_fold = zeros(p)
    r_full = copy(y_centered)
    @inbounds for ii = 1:n_test
        r_full[test_idx[ii]] = 0.0
    end
    y_pred_full = Vector{Float64}(undef, n)  # reusable test-prediction buffer
    mse_row = Vector{Float64}(undef, n_lambda)

    # Active set for sparse inner sweeps
    active = Int[]
    sizehint!(active, p)

    # Count non-converged (λ) steps so the caller can warn when coordinate
    # descent hits max_iter without converging. Silent non-convergence would
    # produce biased cv_mse values at the affected λ.
    n_nonconverged = 0

    for (li, λ) in enumerate(lambdas)
        converged_at_li = false
        for outer_iter = 1:max_iter
            # Cooperative yield at the top of every outer iteration. This is
            # a safety net: even with a dedicated :interactive thread for the
            # monitor, additional yield points ensure the scheduler has
            # frequent opportunities to run the monitor task in case of thread
            # affinity or scheduling quirks. Cost is ~10ns when no task is
            # waiting, which is negligible compared to the full-pass dot
            # products (milliseconds each).
            yield()
            Threads.atomic_add!(outer_iters_used[fold], 1)
            max_change = 0.0

            # Full pass over all coordinates (skip intercept at j=1).
            #
            # LASSO coordinate-descent update for the per-fold loss
            #   (1/(2 n_train))·‖y_train − Φ_train c‖² + λ ‖c‖_1
            # is:
            #   c_j_new = soft_threshold(Φ_jᵀ r̃, n_train·λ) / col_sq_train[j]
            # where r̃ = r + Φ_j c_j_old (the "partial residual").
            #
            # `Φ_scaled` is standardized so ‖Φ_scaled[:, j]‖² = n over the *full*
            # sample, but ‖Φ_scaled,train[:, j]‖² = col_sq_train[j] varies with
            # the CV fold partition. Dividing by n_train instead of
            # col_sq_train[j] would produce a fixed point that is NOT the LASSO
            # optimum; the KKT condition Φ_jᵀ r* = n_train·λ would be violated
            # by a factor of (col_sq_train[j]/n_train − 1)·c_j. The violation is
            # small (~1-2%) for low-order Hermite basis functions and larger
            # (~10-15%) for high-order ones because col_sq_train[j] has higher
            # variance for the latter.
            @inbounds for j = 2:p
                Φj = view(Φ_scaled, :, j)
                # dot(Φj, r_full) is contiguous; test rows of r_full are zero
                # so they contribute 0, equivalent to dot over training rows.
                numerator = dot(Φj, r_full) + col_sq_train[j] * c_fold[j]
                c_new = soft_threshold(numerator, n_train * λ) / col_sq_train[j]
                if c_new != c_fold[j]
                    delta = c_fold[j] - c_new
                    @inbounds for ii = 1:n_train
                        i = train_idx[ii]
                        r_full[i] += Φj[i] * delta
                    end
                    max_change = max(max_change, abs(delta))
                    c_fold[j] = c_new
                end
            end

            # Unpenalized intercept (column 1)
            @inbounds begin
                Φ1 = view(Φ_scaled, :, 1)
                z_1 = (dot(Φ1, r_full) + col_sq_train[1] * c_fold[1]) / col_sq_train[1]
                if z_1 != c_fold[1]
                    delta = c_fold[1] - z_1
                    @inbounds for ii = 1:n_train
                        i = train_idx[ii]
                        r_full[i] += Φ1[i] * delta
                    end
                    max_change = max(max_change, abs(delta))
                    c_fold[1] = z_1
                end
            end

            # Build active set
            empty!(active)
            @inbounds for j = 2:p
                if c_fold[j] != 0.0
                    push!(active, j)
                end
            end

            # Active set sweeps — capped at max_inner_iter (typically 200) to
            # prevent pathological λ values from burning thousands of sweeps
            # without converging. If the inner loop hits the cap, the outer
            # loop's subsequent full pass will either fully converge (max_change
            # drops below tol) or continue iterating. Either way, progress on
            # the outer counter is visible every ~200 sweeps instead of every
            # ~10000 sweeps. Correctness is preserved because convergence is
            # governed by the outer loop's tol check, not the inner cap.
            for _ = 1:max_inner_iter
                as_change = 0.0
                # Same LASSO update as the full pass above: divide the
                # soft-threshold by col_sq_train[j], not n_train.
                @inbounds for j in active
                    Φj = view(Φ_scaled, :, j)
                    numerator = dot(Φj, r_full) + col_sq_train[j] * c_fold[j]
                    c_new = soft_threshold(numerator, n_train * λ) / col_sq_train[j]
                    if c_new != c_fold[j]
                        delta = c_fold[j] - c_new
                        @inbounds for ii = 1:n_train
                            i = train_idx[ii]
                            r_full[i] += Φj[i] * delta
                        end
                        as_change = max(as_change, abs(delta))
                        c_fold[j] = c_new
                    end
                end
                @inbounds begin
                    Φ1 = view(Φ_scaled, :, 1)
                    z_1 = (dot(Φ1, r_full) + col_sq_train[1] * c_fold[1]) / col_sq_train[1]
                    if z_1 != c_fold[1]
                        delta = c_fold[1] - z_1
                        @inbounds for ii = 1:n_train
                            i = train_idx[ii]
                            r_full[i] += Φ1[i] * delta
                        end
                        as_change = max(as_change, abs(delta))
                        c_fold[1] = z_1
                    end
                end
                if as_change < tol
                    break
                end
            end

            if max_change < tol
                converged_at_li = true
                break
            end
        end
        if !converged_at_li
            n_nonconverged += 1
        end

        # Test-set prediction via active coefficients
        @inbounds begin
            fill!(y_pred_full, c_fold[1])
            for j = 2:p
                cj = c_fold[j]
                if cj != 0.0
                    Φj = view(Φ_scaled, :, j)
                    @simd for i = 1:n
                        y_pred_full[i] += Φj[i] * cj
                    end
                end
            end
        end
        sum_sq = 0.0
        @inbounds for ii = 1:n_test
            i = test_idx[ii]
            diff = y_centered[i] - y_pred_full[i]
            sum_sq += diff * diff
        end
        mse_row[li] = sum_sq / n_test

        Threads.atomic_add!(lambda_done[fold], 1)
    end

    if n_nonconverged > 0
        @warn "run_cv_fold! fold $fold: coordinate descent hit max_iter ($(max_iter)) without converging at $n_nonconverged / $n_lambda λ values. cv_mse for those λ may be biased."
    end

    return mse_row, c_fold
end

"""
	fit_pce_lasso(Φ, y; n_lambda=50, n_folds=5, lambda_ratio=1e-4)

Fit sparse PCE coefficients via LASSO (L₁ regularization) using coordinate
descent, with K-fold cross-validation for regularization parameter selection.

Used when the oversampling ratio (n_samples / n_basis) is small (< SPARSE_THRESHOLD),
which occurs at higher PCE degrees. The LASSO automatically selects the most
important basis terms and zeros out the rest.

`lambda_ratio` controls the λ path: the smallest λ on the path is
`lambda_max · lambda_ratio`, where lambda_max is the smallest λ at which all
coefficients are zero. The default `1e-4` follows glmnet's convention for the
`n > p` regime (comfortable oversampling), giving the CV curve enough room to
find a minimum interior to the path. For `n ≤ p` (under-determined) use a
larger value like `1e-2` to avoid numerical problems near λ_min.

Performance: Uses the standard residual-tracking trick so each coordinate update
is O(n) rather than O(n×p), giving O(n×p) per full pass. Solutions are warm-started
along the λ path (previous λ's solution initializes the next), which typically
converges in a few iterations per λ step.

Reference: Blatman & Sudret (2011), "Adaptive sparse polynomial chaos expansion
based on least angle regression", J. Comput. Phys.

Returns coefficients c, R², Q² (from CV), and predictions y_pred.
"""
function fit_pce_lasso(
    Φ::Matrix{Float64},
    y::Vector{Float64};
    precomputed::Union{Nothing,NamedTuple} = nothing,
    n_lambda::Int = 50,
    n_folds::Int = 5,
    lambda_ratio::Float64 = 1e-4,
    max_iter::Int = 10000,
    max_inner_iter::Int = 1000,
    tol::Float64 = 1e-7,
    verbose::Bool = false,
)

    t_start = time()
    n, p = size(Φ)
    verbose && println(
        "      LASSO: n=$n samples, p=$p basis functions, $n_folds-fold CV over $n_lambda λ values",
    )

    # --- Input validation ---
    if !all(isfinite, Φ)
        n_bad = count(!isfinite, Φ)
        error(
            "fit_pce_lasso: Φ contains $n_bad non-finite entries (out of $(length(Φ))). " *
            "Likely a basis evaluation issue — check build_design_matrix / Hermite evaluation at extreme Z.",
        )
    end
    if !all(isfinite, y)
        error("fit_pce_lasso: y contains non-finite entries")
    end

    # Center y (don't regularize the intercept — but our Φ[:,1] is the constant
    # term, so we exclude it from penalization)
    y_mean = mean(y)
    y_centered = y .- y_mean

    # Standardize columns of Φ — either reuse a precomputed standardization
    # (shared across QoIs at the same degree, eliminates redundant copies) or
    # compute it internally. Either path gives mathematically identical results
    # since the standardization depends only on Φ, not y.
    precomp = if precomputed === nothing
        standardize_design_matrix(Φ)
    else
        # Sanity check: caller supplied precomputed must match Φ's dimensions.
        # This catches the bug where a caller passes Φ from one degree and
        # precomputed from another.
        size(precomputed.Φ_scaled) == (n, p) || error(
            "fit_pce_lasso: precomputed.Φ_scaled dimensions $(size(precomputed.Φ_scaled)) " *
            "do not match Φ dimensions $((n, p))",
        )
        length(precomputed.col_norms) == p || error(
            "fit_pce_lasso: precomputed.col_norms length $(length(precomputed.col_norms)) " *
            "does not match p=$p",
        )
        length(precomputed.col_sq_norms) == p || error(
            "fit_pce_lasso: precomputed.col_sq_norms length $(length(precomputed.col_sq_norms)) " *
            "does not match p=$p",
        )
        precomputed
    end
    Φ_scaled = precomp.Φ_scaled
    col_norms = precomp.col_norms
    col_sq_norms = precomp.col_sq_norms
    n_dead = precomp.n_dead
    verbose &&
        n_dead > 0 &&
        println(
            "      LASSO: zeroed $n_dead dead columns (norm < $(@sprintf("%.1e", sqrt(eps(Float64)))))",
        )

    # Lambda path: from lambda_max (smallest λ giving all zeros) down
    # lambda_max = max_j |Φ_j' y| / n
    lambda_max = maximum(abs.(Φ_scaled[:, 2:end]' * y_centered)) / n
    if !isfinite(lambda_max) || lambda_max <= 0
        error(
            "fit_pce_lasso: lambda_max = $lambda_max (non-finite or non-positive). " *
            "Something is wrong with the standardized design matrix or y_centered.",
        )
    end
    lambda_min = lambda_max * lambda_ratio
    lambdas = exp.(range(log(lambda_max), log(lambda_min), length = n_lambda))

    # K-fold CV split
    rng = MersenneTwister(42)
    fold_ids = shuffle(rng, repeat(1:n_folds, ceil(Int, n / n_folds))[1:n])

    cv_mse = zeros(n_lambda)

    # Each CV fold walks the full lambda path with warm starts (each lambda
    # reuses the previous lambda's coefficients). Folds run in parallel.
    t_cv_start = time()
    fold_mse = zeros(n_folds, n_lambda)
    fold_wall = Vector{Float64}(undef, n_folds)  # per-fold wall-clock time
    fold_tid = Vector{Int}(undef, n_folds)      # thread each fold ran on
    # Shared progress tracker: lambda_done[fold] = highest λ index finished.
    # Written to by each fold task; read by external monitor (see LASSO_PROGRESS).
    lambda_done = [Threads.Atomic{Int}(0) for _ = 1:n_folds]
    outer_iters_used = [Threads.Atomic{Int}(0) for _ = 1:n_folds]

    # Register this LASSO call with the global progress registry so an
    # external monitor (e.g., the degree-selection-level monitor) can peek
    # at fold progress without requiring a dedicated scheduler thread per call.
    lasso_id = lasso_register!(lambda_done, outer_iters_used, n_lambda, p, time())

    # Masked-residual coordinate descent. Instead of copying
    # Φ_scaled[train_idx, :] per fold (hundreds of MB at degree 4), each fold
    # operates directly on the shared, read-only Φ_scaled and maintains a
    # full-length residual with zeros on test rows. See `run_cv_fold!` for
    # the full inner loop.
    tasks = map(1:n_folds) do fold
        Threads.@spawn begin
            t_fold_start = time()
            fold_tid[fold] = Threads.threadid()
            # run_cv_fold! also returns the fold's final coefficient vector
            # (at the smallest λ on the path). The CV MSE curve uses only
            # mse_row; c_fold is discarded here but is exposed to tests
            # that need to verify per-fold KKT conditions.
            mse_row, _c_fold = run_cv_fold!(
                fold,
                Φ_scaled,
                y_centered,
                col_sq_norms,
                fold_ids,
                lambdas,
                max_iter,
                max_inner_iter,
                tol,
                outer_iters_used,
                lambda_done,
            )
            fold_wall[fold] = time() - t_fold_start
            mse_row
        end
    end
    for (fold, task) in enumerate(tasks)
        fold_mse[fold, :] .= fetch(task)
    end
    # Unregister from the progress registry (LASSO is finished)
    lasso_unregister!(lasso_id)
    cv_wall = time() - t_cv_start
    if verbose
        # Parallel efficiency check:
        #  - Serial time estimate = sum of per-fold times
        #  - Parallel wall time = cv_wall
        #  - Ideal speedup = n_folds; actual = sum(fold_wall) / cv_wall
        #  - If threads are working, each fold runs on a distinct thread ID
        sum_fold = sum(fold_wall)
        speedup = sum_fold / cv_wall
        tid_str = join(fold_tid, ",")
        println(
            @sprintf(
                "      CV: %d folds × %d λ, wall=%.1fs, Σfolds=%.1fs, speedup=%.1fx (tids=[%s], nthreads=%d)",
                n_folds,
                n_lambda,
                cv_wall,
                sum_fold,
                speedup,
                tid_str,
                Threads.nthreads()
            )
        )
    end

    # --- NaN/Inf diagnostic on CV results ---
    # argmin in Julia doesn't handle NaN safely (NaN fails isless comparisons),
    # so we need to find the smallest FINITE cv_mse. First, detect and report
    # any bad fold/λ pairs so we can see where coordinate descent broke down.
    n_bad_fold = count(!isfinite, fold_mse)
    if n_bad_fold > 0
        bad_locs = findall(!isfinite, fold_mse)
        # Report up to 3 examples
        examples =
            [(Tuple(loc), fold_mse[loc]) for loc in bad_locs[1:min(3, length(bad_locs))]]
        verbose && println(
            "      LASSO WARNING: $n_bad_fold non-finite fold MSE entries " *
            "(out of $(n_folds * n_lambda)). First examples: $examples",
        )
    end

    # Compute mean across folds; treat non-finite fold entries as Inf so the
    # average is well-defined. If ALL folds at a λ are bad, that λ's cv_mse = Inf.
    for li = 1:n_lambda
        finite_mask = isfinite.(fold_mse[:, li])
        if any(finite_mask)
            cv_mse[li] = mean(fold_mse[finite_mask, li])
        else
            cv_mse[li] = Inf
        end
    end

    # NaN-safe argmin: pick smallest finite cv_mse
    if !any(isfinite, cv_mse)
        error(
            "fit_pce_lasso: all CV MSE values are non-finite. " *
            "The LASSO fit broke down at every λ — likely an ill-conditioned " *
            "Φ (check dead columns, degree, or oversampling ratio).",
        )
    end
    finite_idx = findall(isfinite, cv_mse)
    best_li_local = argmin(cv_mse[finite_idx])
    best_li = finite_idx[best_li_local]
    λ_best = lambdas[best_li]

    if verbose && best_li == finite_idx[end]
        println(
            "      LASSO NOTE: selected λ is the smallest finite-CV λ — " *
            "CV curve may still be decreasing; consider raising lambda_ratio.",
        )
    end

    # Refit on full data with best lambda, warm-started along λ path.
    # Uses the same views + active set strategy as the CV loop.
    #
    # Note: col_sq_norms[j] == n for all j (by construction of
    # standardize_design_matrix), so dividing z_j by n IS the correct LASSO
    # update here. This is why the full-data refit is unaffected by the bug
    # that used to plague run_cv_fold! — see the LONG comment there.
    c = zeros(p)
    r_full = copy(y_centered)
    active_full = Int[]
    sizehint!(active_full, p)
    n_refit_nonconverged = 0
    for li = 1:best_li
        λ_curr = lambdas[li]
        converged_at_li = false
        for _ = 1:max_iter
            max_change = 0.0
            # Full pass
            @inbounds for j = 2:p
                Φj = view(Φ_scaled, :, j)
                z_j = (dot(Φj, r_full) + col_sq_norms[j] * c[j]) / n
                c_new = soft_threshold(z_j, λ_curr)
                if c_new != c[j]
                    delta = c[j] - c_new
                    @simd for i in eachindex(r_full)
                        r_full[i] += Φj[i] * delta
                    end
                    max_change = max(max_change, abs(delta))
                    c[j] = c_new
                end
            end
            @inbounds begin
                Φ1 = view(Φ_scaled, :, 1)
                z_1 = (dot(Φ1, r_full) + col_sq_norms[1] * c[1]) / col_sq_norms[1]
                if z_1 != c[1]
                    delta = c[1] - z_1
                    @simd for i in eachindex(r_full)
                        r_full[i] += Φ1[i] * delta
                    end
                    c[1] = z_1
                end
            end

            # Active set
            empty!(active_full)
            @inbounds for j = 2:p
                if c[j] != 0.0
                    push!(active_full, j)
                end
            end
            # Inner active-set sweep, capped at max_inner_iter (same rationale
            # as in run_cv_fold!): prevents pathological λ from burning thousands
            # of sweeps without convergence. Outer loop's tol check still
            # governs correctness.
            for _ = 1:max_inner_iter
                as_change = 0.0
                @inbounds for j in active_full
                    Φj = view(Φ_scaled, :, j)
                    z_j = (dot(Φj, r_full) + col_sq_norms[j] * c[j]) / n
                    c_new = soft_threshold(z_j, λ_curr)
                    if c_new != c[j]
                        delta = c[j] - c_new
                        @simd for i in eachindex(r_full)
                            r_full[i] += Φj[i] * delta
                        end
                        as_change = max(as_change, abs(delta))
                        c[j] = c_new
                    end
                end
                @inbounds begin
                    Φ1 = view(Φ_scaled, :, 1)
                    z_1 = (dot(Φ1, r_full) + col_sq_norms[1] * c[1]) / col_sq_norms[1]
                    if z_1 != c[1]
                        delta = c[1] - z_1
                        @simd for i in eachindex(r_full)
                            r_full[i] += Φ1[i] * delta
                        end
                        as_change = max(as_change, abs(delta))
                        c[1] = z_1
                    end
                end
                if as_change < tol
                    break
                end
            end

            if max_change < tol
                converged_at_li = true
                break
            end
        end
        if !converged_at_li
            n_refit_nonconverged += 1
        end
    end

    if n_refit_nonconverged > 0
        @warn "fit_pce_lasso full-data refit: hit max_iter ($(max_iter)) without converging at $n_refit_nonconverged / $best_li λ values. Final coefficients may be biased; consider raising max_iter or loosening tol."
    end

    # Rescale coefficients back to original Φ scale
    for j = 2:p
        if col_norms[j] > 0
            c[j] /= col_norms[j]
        end
    end
    # Adjust constant term for y centering
    c[1] += y_mean

    y_pred = Φ * c
    residuals = y .- y_pred

    ss_res = sum(residuals .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)

    # R² and Q² — return -Inf (not NaN) on failure so downstream degree
    # selection can rank degrees correctly (-Inf < any finite Q²).
    R2 = isfinite(ss_res) && ss_tot > 0 ? 1.0 - ss_res / ss_tot : -Inf
    Q2 = isfinite(cv_mse[best_li]) && ss_tot > 0 ? 1.0 - cv_mse[best_li] * n / ss_tot : -Inf

    n_nonzero = count(abs.(c[2:end]) .> 1e-12)
    t_total = round(time() - t_start, digits = 1)
    verbose && print("\r" * " "^60 * "\r")  # clear progress line
    println(
        @sprintf(
            "      LASSO: λ=%.4e, %d/%d nonzero, R²=%.4f, Q²=%.4f  [%.1fs]",
            λ_best,
            n_nonzero,
            p - 1,
            R2,
            Q2,
            t_total
        )
    )

    return c, R2, Q2, y_pred
end

# ============================================================================
# Step 4: Bootstrap Confidence Intervals for Sobol Indices
# ============================================================================

"""
	bootstrap_sobol_ci(Z, y, mop, groups, n_boot; ci_level=0.95, use_sparse=false,
		Φ=nothing, c_hat=nothing,
		jack_Φ=nothing, jack_c=nothing, jack_h_diag=nothing,
		jack_residuals=nothing, jack_support=nothing)

Compute BCa (bias-corrected and accelerated) bootstrap confidence intervals
for first-order and total group Sobol indices (Efron, JASA 1987).

BCa adjusts the simple percentile endpoints for:
  - Bias: z₀ = Φ⁻¹(proportion of bootstrap replicates < point estimate)
  - Acceleration: â = Σ Lᵢ³ / (6·(Σ Lᵢ²)^{3/2}), where Lᵢ are jackknife
	influence values computed via the efficient OLS rank-1 downdate on the
	`jack_*` inputs:
	c_{(-i)} = jack_c - (jack_Φᵀ jack_Φ)⁻¹ jack_Φᵢᵀ eᵢ / (1 - hᵢᵢ)

The adjusted percentile positions are:
  α_BCa = Φ( z₀ + (z₀ + z_α) / (1 - â·(z₀ + z_α)) )

Jackknife input separation
--------------------------
For the OLS path, the bootstrap replicate design matrix Φ and the jackknife
design matrix jack_Φ are the same. For the LASSO path they differ:
  - `Φ`       = full basis matrix (used for bootstrap replicate fits)
  - `c_hat`   = LASSO coefficient vector in full-basis space
  - `jack_Φ`  = submatrix Φ[:, S] on the LASSO-selected support S
  - `jack_c`  = OLS refit of `jack_Φ \\ y` (well-defined smooth statistic)
  - `jack_h_diag`  = diag(jack_Φ · (jack_Φᵀ jack_Φ)⁻¹ · jack_Φᵀ)
  - `jack_residuals` = y − jack_Φ · jack_c
  - `jack_support` = indices into the full-basis c vector identifying jack_Φ's columns

This "active-set OLS" jackknife is a principled approximation for the LASSO
leave-one-out estimator (post-selection OLS on the selected support). The
LASSO coefficient map is non-smooth at sparsity boundaries, so the raw
LASSO LOO is ill-defined for BCa; the active-set OLS version IS smooth as
long as the selected support is stable under LOO, which holds empirically
for moderate signal-to-noise.

For the OLS path, pass `jack_Φ = Φ`, `jack_c = c_hat`, `jack_support = nothing`
(equivalent to `1:p`), and the jackknife is the standard OLS rank-1 downdate.

When jack_* are not provided, falls back to simple percentile intervals.

Returns (S1_median, S1_lo, S1_hi, ST_median, ST_lo, ST_hi) each of length n_groups.
"""
function bootstrap_sobol_ci(
    Z::Matrix{Float64},
    y::Vector{Float64},
    mop::MultiOrthoPoly,
    groups::Vector{ParamGroup},
    n_boot::Int;
    ci_level::Float64 = 0.95,
    use_sparse::Bool = false,
    Φ::Union{Nothing,Matrix{Float64}} = nothing,
    c_hat::Union{Nothing,Vector{Float64}} = nothing,
    jack_Φ::Union{Nothing,Matrix{Float64}} = nothing,
    jack_c::Union{Nothing,Vector{Float64}} = nothing,
    jack_h_diag::Union{Nothing,Vector{Float64}} = nothing,
    jack_residuals::Union{Nothing,Vector{Float64}} = nothing,
    jack_support::Union{Nothing,Vector{Int}} = nothing,
)

    n = length(y)
    n_groups = length(groups)
    α_lo = (1.0 - ci_level) / 2.0
    α_hi = 1.0 - α_lo

    boot_S1 = Matrix{Float64}(undef, n_boot, n_groups)
    boot_ST = Matrix{Float64}(undef, n_boot, n_groups)

    # Pre-generate bootstrap indices with a single RNG for reproducibility,
    # then run replicates in parallel across threads.
    rng = MersenneTwister(123)
    boot_indices = [rand(rng, 1:n, n) for _ = 1:n_boot]

    # If the caller supplied Φ, use Φ[idx, :] instead of rebuilding the
    # design matrix via build_design_matrix(Z[idx, :], mop) per replicate.
    #
    # Rationale (verified empirically, not assumed):
    #   1. PolyChaos.build_design_matrix is thread-safe on shared mop
    #      (verified in test_polychaos_threadsafety.jl: 96.9M entries
    #      bit-for-bit identical between sequential and 12-thread parallel
    #      execution with shared MultiOrthoPoly).
    #   2. However, the original bootstrap loop used to hang at production
    #      scale (n=8192, degree 3, 1000 replicates). The likely cause is
    #      GC pressure from per-replicate allocations: each call to
    #      build_design_matrix allocates a ~64 MB Φ_b matrix plus ~1 MB
    #      Z_b = Z[idx, :] plus ~10-20 MB LAPACK workspace. Across 1000
    #      replicates that's ~80 GB of allocation churn hitting Julia's
    #      stop-the-world GC.
    #   3. Φ[idx, :] reuses the already-computed Φ (no rebuild) and is
    #      mathematically identical to build_design_matrix(Z[idx, :], mop)
    #      by the row-wise definition of the PCE basis. This equivalence
    #      is verified empirically in test_bootstrap_phi_equivalence.jl on
    #      the actual publication data: all 40M+ entries bit-for-bit equal
    #      at degrees 2, 3, and 4.
    #
    # The chunked per-task preallocation below (Fix A) further reduces GC
    # pressure by reusing Φ_b buffers across all replicates in a chunk.
    use_direct_phi = !isnothing(Φ)
    p_basis = mop.dim

    # Chunked task parallelism with per-chunk preallocated buffers.
    #
    # Rationale: the naive @spawn-per-replicate approach allocates a fresh
    # Φ_b (n × p_basis Float64 matrix, ~64 MB at degree 3) per replicate,
    # producing ~80 GB of GC churn per bootstrap call. Julia's stop-the-world
    # GC pauses eat most of the wall time. By partitioning 1:n_boot into
    # n_chunks contiguous ranges and giving each chunk task its own
    # preallocated Φ_b and y_b, we allocate n_chunks times (not n_boot times),
    # reducing GC pressure by ~100x.
    #
    # Each chunk is still processed by a separate task, preserving
    # parallelism across threads. The computation inside (Φ_b \ y_b via QR)
    # is unchanged — this fix is purely about where memory is allocated,
    # so results are bit-for-bit identical to the per-replicate approach.
    n_chunks = min(Threads.nthreads(), n_boot)
    chunk_len = cld(n_boot, n_chunks)

    t_boot = time()
    chunk_tasks = map(1:n_chunks) do chunk_i
        Threads.@spawn begin
            start_b = (chunk_i - 1) * chunk_len + 1
            end_b = min(chunk_i * chunk_len, n_boot)

            # Preallocate per-chunk buffers (reused across all replicates
            # in this chunk — allocated once instead of per replicate).
            #
            # The Cholesky normal-equations path (OLS only, not sparse)
            # additionally needs G (p×p), b_vec (p), and c_sol (p). These
            # are small compared to Φ_b and preallocating them avoids ~15 KB
            # of allocation per replicate.
            Φ_b_buf = use_direct_phi ? Matrix{Float64}(undef, n, p_basis) : nothing
            y_b_buf = Vector{Float64}(undef, n)
            G_buf = use_sparse ? nothing : Matrix{Float64}(undef, p_basis, p_basis)
            b_vec_buf = use_sparse ? nothing : Vector{Float64}(undef, p_basis)
            c_sol_buf = use_sparse ? nothing : Vector{Float64}(undef, p_basis)

            for b = start_b:end_b
                idx = boot_indices[b]

                # Copy y rows into preallocated buffer
                @inbounds for k = 1:n
                    y_b_buf[k] = y[idx[k]]
                end

                # Build Φ_b: either copy into preallocated buffer (fast path)
                # or rebuild from Z (legacy fallback for callers not passing Φ)
                Φ_b_local = if use_direct_phi
                    # Column-major friendly copy: iterate j (outer), k (inner).
                    # Inside column j of Φ_b_buf we write sequentially.
                    # Reads from Φ[idx[k], j] are strided (random within column)
                    # but column j of Φ is ≤64 KB and fits in L2/L3 cache.
                    @inbounds for j = 1:p_basis
                        for k = 1:n
                            Φ_b_buf[k, j] = Φ[idx[k], j]
                        end
                    end
                    Φ_b_buf
                else
                    # Legacy path: rebuild per replicate. NOT thread-safe
                    # with shared mop, so this branch only works correctly
                    # if the caller arranges for sequential execution.
                    build_design_matrix(Z[idx, :], mop)
                end

                if use_sparse
                    c_b, _, _, _ = fit_pce_lasso(Φ_b_local, y_b_buf)
                else
                    # Normal equations + Cholesky (Fix B).
                    #
                    # Replaces `Φ_b \ y_b` (LAPACK QR via dgels) with the
                    # symmetric-positive-definite normal-equations approach.
                    # Verified on real data at degree 3 with n=8192, p=969:
                    #
                    #   - κ(Φ'Φ) = 11.21 (extremely well-conditioned)
                    #   - Cholesky relative error bound κ × ε_mach ≈ 2.5e-15
                    #   - Empirical max |Δ| vs QR: 1.33e-15 on boot_S1,
                    #     1.55e-15 on boot_ST (i.e., 1.5× ε_mach — equivalent
                    #     to QR at the last bit of Float64 precision)
                    #   - Empirical speedup: 19.8× per bootstrap call at
                    #     n_boot=100 (124s QR → 6.3s Cholesky on 12 threads)
                    #
                    # See test_bootstrap_fix_b_equivalence.jl for the full
                    # verification run. At this condition number and with this
                    # level of empirical agreement, Cholesky is indistinguishable
                    # from QR for publication purposes while being ~20× faster.
                    #
                    # Implementation:
                    #   1. BLAS.syrk! computes G = Φ_b'Φ_b into the upper
                    #      triangle of G_buf (~½ the ops of full gemm by
                    #      exploiting symmetry).
                    #   2. mul! for b = Φ_b'y_b (cheap gemv).
                    #   3. copyto! into c_sol_buf (ldiv! overwrites its RHS).
                    #   4. cholesky! decomposes G in-place (reads :U, overwrites
                    #      it with the Cholesky factor). Throws PosDefException
                    #      loudly if the matrix is not SPD — no silent failure.
                    #   5. ldiv! performs both triangular solves in place,
                    #      writing the solution into c_sol_buf.
                    BLAS.syrk!('U', 'T', 1.0, Φ_b_local, 0.0, G_buf)
                    mul!(b_vec_buf, Φ_b_local', y_b_buf)
                    copyto!(c_sol_buf, b_vec_buf)
                    F = cholesky!(Symmetric(G_buf, :U))
                    ldiv!(F, c_sol_buf)
                    c_b = c_sol_buf
                end

                S1_b, ST_b, _ = compute_group_sobol(c_b, mop, groups)
                # Distinct b across chunks → no write conflict on these rows
                @inbounds boot_S1[b, :] .= S1_b
                @inbounds boot_ST[b, :] .= ST_b
            end
        end
    end
    # Fetch (not wait) to propagate any exceptions from the chunk tasks
    for task in chunk_tasks
        fetch(task)
    end
    println(
        @sprintf(
            "    bootstrap: %d replicates in %d chunks [%.1fs]",
            n_boot,
            n_chunks,
            time() - t_boot
        )
    )

    S1_median = vec(median(boot_S1, dims = 1))
    ST_median = vec(median(boot_ST, dims = 1))
    S1_lo = zeros(n_groups)
    S1_hi = zeros(n_groups)
    ST_lo = zeros(n_groups)
    ST_hi = zeros(n_groups)

    use_bca =
        !isnothing(c_hat) &&
        !isnothing(jack_Φ) &&
        !isnothing(jack_c) &&
        !isnothing(jack_h_diag) &&
        !isnothing(jack_residuals)

    if use_bca
        # --- BCa correction (Efron 1987) ---
        # The point estimate for the statistic is computed from c_hat in full
        # basis space. This is what the bootstrap replicates are targeting,
        # so the bias correction z₀ is referenced to this value.
        S1_hat, ST_hat, _ = compute_group_sobol(c_hat, mop, groups)

        # Jackknife Sobol indices via OLS rank-1 downdate in the jackknife
        # subspace (either full basis for OLS, or the LASSO-selected support
        # for LASSO). For LASSO, jack_c is the *active-set OLS refit* which
        # is a smooth statistic suitable for the BCa acceleration constant.
        # After computing the per-i downdated coefficient vector in the
        # jackknife subspace, we embed it back into full basis space (zero-
        # padding at non-support positions) before computing Sobol indices.
        n_jack = size(jack_Φ, 1)
        n_jack == n || error("bootstrap_sobol_ci: jack_Φ has $n_jack rows, expected $n")
        p_full = length(c_hat)
        p_jack = size(jack_Φ, 2)
        length(jack_c) == p_jack || error(
            "bootstrap_sobol_ci: jack_c length $(length(jack_c)) != jack_Φ columns $p_jack",
        )
        length(jack_h_diag) == n_jack || error(
            "bootstrap_sobol_ci: jack_h_diag length $(length(jack_h_diag)) != n $n_jack",
        )
        length(jack_residuals) == n_jack || error(
            "bootstrap_sobol_ci: jack_residuals length $(length(jack_residuals)) != n $n_jack",
        )
        if jack_support !== nothing
            length(jack_support) == p_jack || error(
                "bootstrap_sobol_ci: jack_support length $(length(jack_support)) != jack_Φ columns $p_jack",
            )
        end

        # p_jack × n factor used by the rank-1 OLS LOO downdate
        G = Symmetric(jack_Φ' * jack_Φ) \ jack_Φ'

        jack_S1 = Matrix{Float64}(undef, n_jack, n_groups)
        jack_ST = Matrix{Float64}(undef, n_jack, n_groups)
        # Embarrassingly parallel: each i is independent. Each iteration
        # writes to its own row of jack_S1 / jack_ST (different rows = no
        # write conflicts). Threads.@threads statically partitions the loop
        # across default-pool threads for ~12x speedup on this chip.
        t_jack = time()
        Threads.@threads for i = 1:n_jack
            # c_loo in jackknife subspace; use view to avoid a per-iteration
            # column copy of G (saves ~O(p_jack) allocations per i).
            Gi = view(G, :, i)
            scale = jack_residuals[i] / (1.0 - jack_h_diag[i])
            if jack_support === nothing
                # Jackknife subspace == full basis space
                c_loo_full = jack_c .- scale .* Gi
            else
                # Embed active-set LOO vector into full basis space
                c_loo_full = zeros(p_full)
                @inbounds for (k, j) in enumerate(jack_support)
                    c_loo_full[j] = jack_c[k] - scale * Gi[k]
                end
            end
            S1_i, ST_i, _ = compute_group_sobol(c_loo_full, mop, groups)
            jack_S1[i, :] = S1_i
            jack_ST[i, :] = ST_i
        end
        println(
            @sprintf(
                "    jackknife: %d iterations [%.1fs] (jackknife space: p=%d)",
                n_jack,
                time() - t_jack,
                p_jack
            )
        )

        std_normal = Normal()
        z_α_lo = quantile(std_normal, α_lo)
        z_α_hi = quantile(std_normal, α_hi)

        for g = 1:n_groups
            S1_lo[g], S1_hi[g] = _bca_interval(
                boot_S1[:, g],
                S1_hat[g],
                jack_S1[:, g],
                z_α_lo,
                z_α_hi,
                n_boot,
            )
            ST_lo[g], ST_hi[g] = _bca_interval(
                boot_ST[:, g],
                ST_hat[g],
                jack_ST[:, g],
                z_α_lo,
                z_α_hi,
                n_boot,
            )
        end
    else
        for g = 1:n_groups
            S1_lo[g] = quantile(boot_S1[:, g], α_lo)
            S1_hi[g] = quantile(boot_S1[:, g], α_hi)
            ST_lo[g] = quantile(boot_ST[:, g], α_lo)
            ST_hi[g] = quantile(boot_ST[:, g], α_hi)
        end
    end

    return S1_median, S1_lo, S1_hi, ST_median, ST_lo, ST_hi
end

"""
	_bca_interval(boot_samples, θ_hat, jack_samples, z_α_lo, z_α_hi, n_boot)

Compute BCa-adjusted confidence interval endpoints for a single scalar statistic.

References:
  Efron, B. (1987). "Better Bootstrap Confidence Intervals", JASA 82(397), 171–185.
  Efron, B. & Tibshirani, R. (1993). An Introduction to the Bootstrap, Ch. 14.
"""
function _bca_interval(
    boot_samples::Vector{Float64},
    θ_hat::Float64,
    jack_samples::Vector{Float64},
    z_α_lo::Float64,
    z_α_hi::Float64,
    n_boot::Int,
)

    std_normal = Normal()

    # Bias correction: z₀ = Φ⁻¹(#{boot < θ̂} / B)
    prop_below = count(boot_samples .< θ_hat) / n_boot
    prop_below = clamp(prop_below, 0.5 / n_boot, 1.0 - 0.5 / n_boot)
    z0 = quantile(std_normal, prop_below)

    # Acceleration: â from jackknife influence values Lᵢ = θ̄_{(·)} - θ_{(-i)}
    L = mean(jack_samples) .- jack_samples
    denom = sum(L .^ 2)
    if denom < eps()
        a_hat = 0.0  # degenerate: constant jackknife → BCa reduces to percentile
    else
        a_hat = sum(L .^ 3) / (6.0 * denom^1.5)
    end

    # Adjusted percentiles: α_BCa = Φ( z₀ + (z₀ + z_α) / (1 - â(z₀ + z_α)) )
    num_lo = z0 + z_α_lo
    α_bca_lo = cdf(std_normal, z0 + num_lo / (1.0 - a_hat * num_lo))
    num_hi = z0 + z_α_hi
    α_bca_hi = cdf(std_normal, z0 + num_hi / (1.0 - a_hat * num_hi))

    α_bca_lo = clamp(α_bca_lo, 0.5 / n_boot, 1.0 - 0.5 / n_boot)
    α_bca_hi = clamp(α_bca_hi, 0.5 / n_boot, 1.0 - 0.5 / n_boot)

    return quantile(boot_samples, α_bca_lo), quantile(boot_samples, α_bca_hi)
end

# ============================================================================
# Step 5: Convergence Diagnostics
# ============================================================================

"""
	convergence_analysis(Z, y, groups, degree, fractions)

Compute Sobol indices on nested Sobol-prefix subsamples of the data to assess
convergence. `fractions` specifies the fraction of the full dataset to use at
each level; these should be powers of two so each prefix is a well-distributed
(t, m, s)-net of the underlying Sobol sequence.

The fit method at each subsample is chosen based on the LOCAL oversampling
ratio `ns / mop_sub.dim`, not on the full-data method. This prevents OLS from
being used at rank-deficient subsample sizes (e.g., at degree 3 with
ns = 1024, n_basis = 969, oversampling ≈ 1.06 — ill-posed for OLS). Previously,
if the full-data fit used OLS (oversampling > SPARSE_THRESHOLD on the full
data), the subsample fits inherited OLS regardless of their own oversampling,
producing noisy/garbage indices at small ns and triggering false "not
converged" diagnostics.

Returns (sample_sizes, S1_by_size) where S1_by_size[i] is the n_groups-length
first-order index vector at subsample size sample_sizes[i].
"""
function convergence_analysis(
    Z::Matrix{Float64},
    y::Vector{Float64},
    groups::Vector{ParamGroup},
    degree::Int,
    fractions::Vector{Float64},
)

    n = length(y)
    n_groups = length(groups)
    sample_sizes = [round(Int, f * n) for f in fractions]

    S1_by_size = Vector{Vector{Float64}}(undef, length(sample_sizes))

    conv_tasks = map(enumerate(sample_sizes)) do (i, ns)
        Threads.@spawn begin
            Z_sub = Z[1:ns, :]
            y_sub = y[1:ns]

            ops = [GaussOrthoPoly(degree) for _ = 1:N_PARAMS]
            mop_sub = MultiOrthoPoly(ops, degree)
            Φ_sub = build_design_matrix(Z_sub, mop_sub)

            # Choose OLS or LASSO based on the LOCAL oversampling ratio so
            # that subsamples with few effective observations per basis
            # function get sparse fits regardless of how the full-data fit
            # was configured.
            oversampling_local = ns / mop_sub.dim
            if oversampling_local < SPARSE_THRESHOLD
                c, _, _, _ = fit_pce_lasso(Φ_sub, y_sub)
            else
                c = Φ_sub \ y_sub
            end

            S1, _, _ = compute_group_sobol(c, mop_sub, groups)
            S1
        end
    end
    for (i, task) in enumerate(conv_tasks)
        S1_by_size[i] = fetch(task)
    end

    return sample_sizes, S1_by_size
end

"""
	assess_convergence(S1_by_size, sample_sizes; tol=0.02)

Robust convergence diagnostic for Sobol index sequences computed at increasing
subsample sizes. Checks three criteria across ALL subsample points:

1. Terminal stability: max |ΔS₁| between the last two sizes < tol
2. No oscillation: ≤ 1 sign change in the first-difference sequence over the
   last 4 subsample sizes (for groups with S₁ > tol)
3. No monotonic drift: |ΔS₁[end]| ≤ |ΔS₁[end-1]| + tol/2 — the step size
   must be shrinking (or at least not growing beyond noise tolerance)

Returns (converged, max_terminal_change, n_oscillating, n_drifting, details).
"""
function assess_convergence(
    S1_by_size::Vector{Vector{Float64}},
    sample_sizes::Vector{Int};
    tol::Float64 = 0.02,
)

    n_sizes = length(S1_by_size)
    n_groups = length(S1_by_size[1])

    max_terminal_change = 0.0
    n_oscillating = 0
    n_drifting = 0
    all_converged = true

    for g = 1:n_groups
        trajectory = [S1_by_size[k][g] for k = 1:n_sizes]

        # Skip groups with negligible S1 (noise-dominated)
        if maximum(abs.(trajectory)) < tol
            continue
        end

        diffs = diff(trajectory)

        # Criterion 1: Terminal stability
        terminal_change = abs(diffs[end])
        max_terminal_change = max(max_terminal_change, terminal_change)
        terminal_ok = terminal_change < tol

        # Criterion 2: Oscillation (sign changes in recent first-differences)
        n_check = min(4, n_sizes)
        recent_diffs = diffs[max(1, end-n_check+2):end]
        oscillation_ok = true
        if length(recent_diffs) >= 2
            nonzero_signs = filter(!iszero, sign.(recent_diffs))
            n_sign_changes =
                length(nonzero_signs) >= 2 ?
                count(
                    nonzero_signs[i] != nonzero_signs[i-1] for i = 2:length(nonzero_signs)
                ) : 0
            oscillation_ok = n_sign_changes <= 1
        end
        if !oscillation_ok
            n_oscillating += 1
        end

        # Criterion 3: Drift (step size must not be growing)
        drift_ok = true
        if length(recent_diffs) >= 2
            abs_diffs = abs.(recent_diffs)
            drift_ok = abs_diffs[end] <= abs_diffs[end-1] + tol / 2
        end
        if !drift_ok
            n_drifting += 1
        end

        if !(terminal_ok && oscillation_ok && drift_ok)
            all_converged = false
        end
    end

    return (
        converged = all_converged,
        max_terminal_change = max_terminal_change,
        n_oscillating = n_oscillating,
        n_drifting = n_drifting,
    )
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("\n" * "="^70)
    println("  Group Sobol Indices via Polynomial Chaos Expansion")
    println("="^70 * "\n")

    # --- Threading configuration ---
    # We parallelize with Threads.@spawn at the CV fold, QoI, bootstrap, and
    # convergence-subsample levels. BLAS must be set to 1 thread to avoid
    # oversubscription: with N Julia threads each calling multi-threaded BLAS,
    # you get N × BLAS_threads OS threads fighting for cores, which typically
    # makes the run slower than sequential.
    #
    # The BLAS thread setting is a process-global side effect — if this file
    # is `include`d from another context rather than run as a script, we need
    # to restore the original setting on exit. Wrap the entire body in
    # try/finally so the restore runs even if main() throws.
    n_jl_threads = Threads.nthreads()
    blas_threads_before = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    println("Threading: Julia=$n_jl_threads, BLAS set to 1 (was $blas_threads_before)")
    if n_jl_threads == 1
        println(
            "  ⚠ Julia is running with 1 thread. Launch with --threads=auto " *
            "or JULIA_NUM_THREADS=N for parallelism.",
        )
    end

    try  # see matching `finally` at the end of main() that restores BLAS threads

        # --- Load data ---
        param_csv = readdlm(joinpath(RESULTS_DIR, "parameters.csv"), ','; header = true)
        param_matrix = Float64.(param_csv[1][:, 2:end])  # [n_samples × 16]
        n_total = size(param_matrix, 1)
        println("Loaded $n_total parameter samples")

        scalar_csv = readdlm(joinpath(RESULTS_DIR, "scalar_qois.csv"), ','; header = true)
        scalar_data = scalar_csv[1]

        success_mask = scalar_data[:, 2] .!= "FAILED"
        n_success = sum(success_mask)
        println("Successful runs: $n_success / $n_total")

        sample_ids = Int.(scalar_data[success_mask, 1])
        qoi_data = Float64.(scalar_data[success_mask, 3:end])
        param_data = param_matrix[sample_ids, :]

        # --- Step 1: Transform to standard normal space ---
        println("\nTransforming to independent standard normal space...")
        Z = transform_to_standard_normal(param_data, GROUPS)

        # Verify: columns should be approximately N(0,1)
        println("  Verification (should be ~0 mean, ~1 std):")
        for (gi, g) in enumerate(GROUPS)
            for idx in g.indices
                μ_z = mean(Z[:, idx])
                σ_z = std(Z[:, idx])
                println(
                    @sprintf("    Z_%d (%-22s): μ=%+.4f, σ=%.4f", idx, g.name, μ_z, σ_z)
                )
            end
        end

        # Verify independence between groups
        println("\n  Max cross-group |correlation|:")
        corr_Z = cor(Z)
        max_cross = 0.0
        for i = 1:N_PARAMS, j = (i+1):N_PARAMS
            gi = findfirst(g -> i ∈ g.indices, GROUPS)
            gj = findfirst(g -> j ∈ g.indices, GROUPS)
            if gi != gj
                max_cross = max(max_cross, abs(corr_Z[i, j]))
            end
        end
        println(
            @sprintf(
                "    %.4f (should be small, ~0.02–0.05 for n=%d)",
                max_cross,
                n_success
            )
        )

        # --- Step 2: Adaptive PCE degree selection ---
        # Try each degree in PCE_DEGREE_RANGE and pick the one that maximizes
        # the average leave-one-out Q² across all QoIs.
        println("\n" * "-"^70)
        println("  Adaptive PCE Degree Selection")
        println("-"^70)

        group_names = [g.name for g in GROUPS]
        n_groups = length(GROUPS)

        # Collect valid QoI indices
        valid_qois = String[]
        for qname in QOI_NAMES
            qi = findfirst(==(qname), ALL_QOI_NAMES)
            y = qoi_data[:, qi]
            if all(isfinite, y) && std(y) > 0
                push!(valid_qois, qname)
            end
        end

        best_degree = first(PCE_DEGREE_RANGE)
        best_avg_Q2 = -Inf
        degree_results = Dict{Int,NamedTuple}()

        for deg in PCE_DEGREE_RANGE
            ops_d = [GaussOrthoPoly(deg) for _ = 1:N_PARAMS]
            mop_d = MultiOrthoPoly(ops_d, deg)
            n_basis_d = mop_d.dim
            oversampling_d = n_success / n_basis_d

            if oversampling_d < 1.5
                println(
                    @sprintf(
                        "  Degree %d: %d basis functions, oversampling %.1f× — SKIPPED (insufficient samples)",
                        deg,
                        n_basis_d,
                        oversampling_d
                    )
                )
                continue
            end

            Φ_d = build_design_matrix(Z, mop_d)
            use_sparse_d = oversampling_d < SPARSE_THRESHOLD
            fit_method_name = use_sparse_d ? "LASSO" : "OLS"
            println(
                "\n  >> Degree $deg: $n_basis_d basis, $(round(oversampling_d, digits=1))× oversampling — " *
                "fitting $(length(valid_qois)) QoIs in parallel via $fit_method_name...",
            )
            t_deg = time()

            # Standardize the design matrix ONCE per degree (shared across all
            # parallel fit_pce_lasso calls via the `precomputed` kwarg). This
            # eliminates the redundant 317 MB Φ_scaled copy per QoI at degree 4,
            # reducing peak memory from ~20 GB to ~340 MB during degree selection.
            precomp_d = use_sparse_d ? standardize_design_matrix(Φ_d) : nothing
            if use_sparse_d
                println(
                    "    (shared Φ_scaled: $(round(sizeof(precomp_d.Φ_scaled) / 1024^2, digits=1)) MB, " *
                    "$(precomp_d.n_dead) dead columns)",
                )
            end

            # Per-QoI completion counter (read by monitor)
            qois_done_atomic = Threads.Atomic{Int}(0)
            monitor_stop = Threads.Atomic{Bool}(false)

            # Monitor task on the :interactive threadpool. Design notes:
            #
            # 1. Wrap everything in try/catch — if any inner call throws, we want
            #    the error to appear in stderr immediately via @error, not be
            #    swallowed until wait(monitor_task) at the very end.
            #
            # 2. Inner print is also try/wrapped so a single print failure doesn't
            #    kill the whole monitor loop.
            #
            # 3. Poll every 5s instead of sleeping for 15s. We track last_print_time
            #    separately so the visible cadence is still ~15s. Shorter sleeps
            #    give faster response to monitor_stop and are more resilient to
            #    scheduler hiccups or GC pauses that delay a long sleep.
            #
            # 4. prev_state is passed between iterations to support stuck-call
            #    detection inside print_monitor_status.
            monitor_task = Threads.@spawn :interactive begin
                try
                    t0 = time()
                    last_print_time = t0
                    # State for stuck-call detection across cycles
                    prev_state = Dict{
                        Int,
                        NamedTuple{
                            (:min_lam, :sum_iter, :last_change_time),
                            Tuple{Int,Int,Float64},
                        },
                    }()
                    poll_secs = 5.0
                    print_interval = 15.0
                    while !monitor_stop[]
                        sleep(poll_secs)
                        if monitor_stop[]
                            break
                        end
                        now = time()
                        if now - last_print_time >= print_interval
                            last_print_time = now
                            try
                                print_monitor_status(
                                    deg,
                                    t0,
                                    now,
                                    qois_done_atomic,
                                    valid_qois,
                                    prev_state,
                                )
                            catch print_err
                                # A single print failure shouldn't kill the loop.
                                @error "Monitor print failed (continuing)" exception=(
                                    print_err,
                                    catch_backtrace(),
                                )
                            end
                        end
                    end
                catch loop_err
                    # Any exception that escapes the inner try/catch means the loop
                    # itself is broken. Report loudly so the user sees it immediately
                    # instead of waiting for wait(monitor_task) at the end of the degree.
                    @error "Monitor loop crashed — degree-level progress will be silent from here on" exception=(
                        loop_err,
                        catch_backtrace(),
                    )
                end
            end

            # Process QoIs in batches. After the memory refactor, each fold task
            # allocates only O(n) working buffers (r_full, y_pred_full, col_sq_train,
            # c_fold — totalling <100 KB per fold at degree 4), and Φ_scaled is
            # shared across all QoIs via precomp_d. Memory pressure is no longer a
            # concern, but batching still helps bound GC churn from temporary
            # allocations. Override with JULIA_LASSO_BATCH=N if needed.
            batch_size_qoi = parse(
                Int,
                get(ENV, "JULIA_LASSO_BATCH", string(max(1, Threads.nthreads()))),
            )
            println("    (batch_size=$batch_size_qoi QoIs in flight at a time)")

            q2_values = Vector{Float64}(undef, length(valid_qois))
            for batch_start = 1:batch_size_qoi:length(valid_qois)
                batch_end = min(batch_start + batch_size_qoi - 1, length(valid_qois))
                batch_idxs = batch_start:batch_end

                batch_tasks = map(batch_idxs) do qi_idx
                    qname = valid_qois[qi_idx]
                    Threads.@spawn begin
                        t_qoi = time()
                        qi = findfirst(==(qname), ALL_QOI_NAMES)
                        y = qoi_data[:, qi]
                        println(
                            @sprintf(
                                "    [deg %d] QoI %d/%d (%s) start on tid=%d",
                                deg,
                                qi_idx,
                                length(valid_qois),
                                qname,
                                Threads.threadid()
                            )
                        )
                        if use_sparse_d
                            # Share precomputed standardization across all QoIs
                            _, _, Q2_d, _ = fit_pce_lasso(
                                Φ_d,
                                y;
                                precomputed = precomp_d,
                                verbose = true,
                            )
                        else
                            _, _, Q2_d, _, _, _ = fit_pce(Φ_d, y)
                        end
                        println(
                            @sprintf(
                                "    [deg %d] QoI %d/%d (%s) done  Q²=%.4f [%.1fs]",
                                deg,
                                qi_idx,
                                length(valid_qois),
                                qname,
                                Q2_d,
                                time() - t_qoi
                            )
                        )
                        Threads.atomic_add!(qois_done_atomic, 1)
                        Q2_d
                    end
                end

                # Wait for this batch to complete before starting the next.
                # This ensures previous batch's Φ_scaled copies are eligible for GC.
                for (i, qi_idx) in enumerate(batch_idxs)
                    q2_values[qi_idx] = fetch(batch_tasks[i])
                end

                # Force GC between batches so memory is reclaimed before the next
                # batch's allocations begin. Without this, lazy GC can let the
                # resident set grow indefinitely.
                GC.gc(false)  # incremental collection (faster than full)
            end
            # Stop monitor and wait for it
            monitor_stop[] = true
            try
                wait(monitor_task)
            catch wait_err
                @warn "wait(monitor_task) threw during degree $deg cleanup" exception=(
                    wait_err,
                    catch_backtrace(),
                )
            end
            # Belt-and-suspenders: surface any task failure that escaped the inner
            # try/catch inside the monitor body (shouldn't happen, but just in case).
            if istaskfailed(monitor_task)
                try
                    fetch(monitor_task)  # re-raises the stored exception
                catch task_err
                    @error "Monitor task failed during degree $deg (stored exception)" exception=(
                        task_err,
                        catch_backtrace(),
                    )
                end
            end
            println(
                @sprintf(
                    "  >> Degree %d: all %d QoIs finished in %.1fs wall",
                    deg,
                    length(valid_qois),
                    time() - t_deg
                )
            )

            # Count and report any non-finite Q² values (these signal LASSO
            # breakdown or OLS rank deficiency for this QoI+degree combination).
            n_bad_q2 = count(!isfinite, q2_values)
            if n_bad_q2 > 0
                bad_qois =
                    [valid_qois[i] for i in eachindex(q2_values) if !isfinite(q2_values[i])]
                println(
                    "  Degree $deg: $n_bad_q2/$(length(q2_values)) QoIs produced non-finite Q² " *
                    "($(join(bad_qois, ", ")))",
                )
            end

            # Robust aggregation: use the mean of FINITE Q² values. If every QoI
            # broke, avg_Q2 = -Inf and the degree is skipped by the comparison below.
            finite_q2 = filter(isfinite, q2_values)
            avg_Q2 = isempty(finite_q2) ? -Inf : mean(finite_q2)
            min_Q2 = isempty(finite_q2) ? -Inf : minimum(finite_q2)

            fit_method = use_sparse_d ? "LASSO" : "OLS"
            println(
                @sprintf(
                    "  Degree %d: %d basis, %.1f× oversampling, avg Q²=%.4f, min Q²=%.4f (%s, %d/%d usable)",
                    deg,
                    n_basis_d,
                    oversampling_d,
                    avg_Q2,
                    min_Q2,
                    fit_method,
                    length(finite_q2),
                    length(q2_values)
                )
            )

            degree_results[deg] = (
                mop = mop_d,
                Φ = Φ_d,
                n_basis = n_basis_d,
                oversampling = oversampling_d,
                avg_Q2 = avg_Q2,
                use_sparse = use_sparse_d,
                n_bad_q2 = n_bad_q2,
            )

            if avg_Q2 > best_avg_Q2
                best_avg_Q2 = avg_Q2
                best_degree = deg
            end
        end

        if best_avg_Q2 == -Inf
            error(
                "Degree selection failed: no degree in $PCE_DEGREE_RANGE produced a " *
                "finite average Q². Check the non-finite Q² warnings above for details.",
            )
        end

        println(
            @sprintf("\n  → Selected degree %d (avg Q² = %.4f)", best_degree, best_avg_Q2)
        )
        if degree_results[best_degree].n_bad_q2 > 0
            println(
                "  ⚠ Note: selected degree had $(degree_results[best_degree].n_bad_q2) " *
                "QoIs with non-finite Q² — main fit loop may fail for those.",
            )
        end

        dr = degree_results[best_degree]
        mop = dr.mop
        Φ = dr.Φ
        n_basis = dr.n_basis
        use_sparse = dr.use_sparse

        # --- Step 2b: Verify basis normalization ---
        println("\n  Verifying basis orthonormality (E[Ψ²] ≈ 1)...")
        max_z, n_norm_fail, norm_pass, failing_basis_idx, empirical_basis_norms =
            verify_basis_normalization(Φ, mop)
        if norm_pass
            println(@sprintf("    ✓ PASSED (max z-score = %.2f)", max_z))
        else
            println(
                @sprintf(
                    "    ⚠ FAILED (%d basis functions exceed 3σ, max z = %.2f) — will run pre-bootstrap diagnostic to quantify impact per QoI",
                    n_norm_fail,
                    max_z
                )
            )
        end

        # Check condition number
        κ = cond(Φ' * Φ)
        println(@sprintf("  Condition number (ΦᵀΦ): %.2e", κ))
        if κ > 1e10
            @warn "Design matrix is poorly conditioned (κ > 1e10). Results may be unreliable."
        end

        # --- Step 3: Fit PCE and compute Sobol indices for each QoI ---
        println("\n" * "-"^70)
        println("  PCE Fitting & Group Sobol Indices (degree $best_degree)")
        println("-"^70)

        # Storage for results
        all_S_first = Dict{String,Vector{Float64}}()
        all_S_total = Dict{String,Vector{Float64}}()
        all_R2 = Dict{String,Float64}()
        all_Q2 = Dict{String,Float64}()
        all_coeffs = Dict{String,Vector{Float64}}()
        all_predictions = Dict{String,Vector{Float64}}()
        # For BCa: per-QoI jackknife design matrix, OLS coefficients, h_diag,
        # residuals, and support indices. For OLS QoIs these point at the full Φ;
        # for LASSO QoIs they point at an active-set OLS refit (see below).
        all_jack_Φ = Dict{String,Matrix{Float64}}()
        all_jack_c = Dict{String,Vector{Float64}}()
        all_jack_h_diag = Dict{String,Vector{Float64}}()
        all_jack_residuals = Dict{String,Vector{Float64}}()
        all_jack_support = Dict{String,Union{Nothing,Vector{Int}}}()
        all_boot_ci = Dict{String,NamedTuple}()

        # Share a single standardized Φ_scaled across all QoI fits at the selected
        # degree — eliminates one redundant 317 MB copy per QoI at degree 4.
        precomp_main = use_sparse ? standardize_design_matrix(Φ) : nothing
        if use_sparse
            println(
                "  Shared Φ_scaled: $(round(sizeof(precomp_main.Φ_scaled) / 1024^2, digits=1)) MB, " *
                "$(precomp_main.n_dead) dead columns",
            )
        end

        for qname in QOI_NAMES
            qi = findfirst(==(qname), ALL_QOI_NAMES)
            y = qoi_data[:, qi]

            if !all(isfinite, y) || std(y) ≈ 0
                println("\n  $qname: SKIPPED (constant or non-finite)")
                continue
            end

            # Fit PCE (OLS or sparse LASSO depending on oversampling)
            if use_sparse
                c, R2, Q2, y_pred =
                    fit_pce_lasso(Φ, y; precomputed = precomp_main, verbose = true)

                # Active-set OLS refit for the BCa jackknife.
                #
                # The LASSO coefficient map is non-smooth at sparsity boundaries,
                # so the naive OLS-downdate-with-LASSO-residuals used in earlier
                # versions of this code produced an ill-defined jackknife. The
                # correct approach (post-selection OLS, or "active-set jackknife")
                # takes the LASSO-selected support S, refits OLS on Φ[:, S], and
                # uses the resulting smooth estimator's LOO rank-1 downdate as
                # the input to the BCa acceleration constant â. The embedded
                # full-length LOO coefficient vector (zero on non-support
                # coordinates) is then fed into compute_group_sobol.
                #
                # References: post-selection inference literature, e.g.,
                #   Dezeure et al. (2015), "High-dimensional inference...",
                #   Belloni & Chernozhukov (2013), "Least squares after LASSO".
                # Efron (1987) BCa assumes a smooth statistic; this construction
                # gives us exactly that for the LASSO-selected submodel.
                jack_support = findall(j -> abs(c[j]) > 1e-12, 1:length(c))
                if isempty(jack_support)
                    # Degenerate case: LASSO selected nothing. Use a single
                    # dummy constant column so the jackknife machinery still
                    # runs (the resulting S_loo will all be zero, which is the
                    # correct answer for an all-zero fit).
                    jack_support = [1]
                end
                Φ_S = Φ[:, jack_support]
                c_S_OLS = Φ_S \ y
                jack_residuals_q = y .- Φ_S * c_S_OLS
                jack_h_diag_q, _ = compute_ols_hat_diag(Φ_S)

                all_jack_Φ[qname] = Φ_S
                all_jack_c[qname] = c_S_OLS
                all_jack_h_diag[qname] = jack_h_diag_q
                all_jack_residuals[qname] = jack_residuals_q
                all_jack_support[qname] = jack_support
            else
                c, R2, Q2, y_pred, h_diag_fit, residuals_fit = fit_pce(Φ, y)
                # For OLS, the jackknife subspace IS the full basis, so no
                # active-set projection is needed.
                all_jack_Φ[qname] = Φ
                all_jack_c[qname] = c
                all_jack_h_diag[qname] = h_diag_fit
                all_jack_residuals[qname] = residuals_fit
                all_jack_support[qname] = nothing
            end
            all_R2[qname] = R2
            all_Q2[qname] = Q2
            all_coeffs[qname] = c
            all_predictions[qname] = y_pred

            # Compute group Sobol indices
            S_first, S_total, D = compute_group_sobol(c, mop, GROUPS)
            all_S_first[qname] = S_first
            all_S_total[qname] = S_total

            # Display
            println(@sprintf("\n  %s  (R²=%.4f, Q²_LOO=%.4f, D=%.4e)", qname, R2, Q2, D))

            sorted = sortperm(S_first, rev = true)
            for (rank, gi) in enumerate(sorted)
                bar_first = repeat("█", round(Int, S_first[gi] * 40))
                bar_interact =
                    repeat("░", round(Int, max(0, (S_total[gi] - S_first[gi])) * 40))
                println(
                    @sprintf(
                        "    %d. %-22s S=%.4f  S_T=%.4f  %s%s",
                        rank,
                        group_names[gi],
                        S_first[gi],
                        S_total[gi],
                        bar_first,
                        bar_interact
                    )
                )
            end
            println(
                @sprintf(
                    "    Σ S_first = %.4f  (remainder = %.4f = interactions)",
                    sum(S_first),
                    1 - sum(S_first)
                )
            )

            # Verify Sobol index closure properties for disjoint groups:
            #   Σ S_first ≤ 1 (with equality iff no cross-group interactions)
            #   Σ S_total ≥ 1 (with equality iff no cross-group interactions)
            # Both checks catch different failure modes — violating either
            # indicates a problem in the basis orthonormalization, the Sobol
            # computation, or the (ΦᵀΦ)-based variance estimate.
            sum_S1 = sum(S_first)
            sum_ST = sum(S_total)
            if sum_ST < 1.0 - 0.01
                @warn "$qname: Σ S_T = $(round(sum_ST, digits=4)) < 1 — unexpected for total indices"
            end
            if sum_S1 > 1.0 + 0.01
                @warn "$qname: Σ S_first = $(round(sum_S1, digits=4)) > 1 — closure violation (check basis orthonormality)"
            end

            if R2 < 0.95
                @warn "$qname: R² = $(round(R2, digits=4)) < 0.95 — PCE fit may be inadequate"
            end
            if Q2 < 0.90
                @warn "$qname: Q²_LOO = $(round(Q2, digits=4)) < 0.90 — PCE may be overfitting"
            end
        end

        # --- Step 3b: Pre-bootstrap basis-contamination diagnostic ---
        #
        # verify_basis_normalization in Step 2b may flag some basis functions as
        # failing the 3σ check on E[Ψ²_j] = 1. This can happen at high PCE degrees
        # (typically degree 4) because high-order interaction terms like
        # Z_i·Z_j·Z_k·Z_l have heavy-tailed Ψ², and 8192 Owen-scrambled Sobol
        # samples cannot faithfully represent the tails of that distribution
        # in the joint space of the relevant input dimensions. The basis itself
        # is mathematically correct; the issue is sample-size finite-n error on
        # the second-moment estimate for those specific terms.
        #
        # Whether this matters depends entirely on whether those failing basis
        # functions are actually carrying variance in the PCE fit. LASSO usually
        # zeros them out (they contribute little to in-sample fit and pay an L¹
        # penalty). If every QoI has negligible c_j for every failing j, the
        # Sobol indices are robust and we can proceed. If a failing basis
        # function has a large coefficient in some QoI, that QoI's Sobol indices
        # inherit the normalization bias and should not be trusted at the current
        # degree/sample-size combination.
        #
        # This diagnostic runs BEFORE the (hours-long) bootstrap phase so the
        # user can abort degenerate degree-4 runs early.
        if !norm_pass
            println("\n" * "-"^70)
            println("  Pre-bootstrap basis-contamination diagnostic")
            println("-"^70)
            println(
                "  $n_norm_fail of $(n_basis - 1) non-constant basis functions failed the 3σ",
            )
            println("  orthonormality check in Step 2b. This check quantifies how much of")
            println("  each QoI's PCE variance D = Σ_{α≠0} c_α² comes from those failing")
            println("  basis functions. If a large fraction does, the Sobol indices for")
            println("  that QoI are biased and the bootstrap CI will inherit the bias.")
            println()
            println(
                @sprintf(
                    "  Abort threshold: %.1f%% of variance from biased basis",
                    100 * BIASED_BASIS_ABORT_THRESHOLD
                )
            )
            println("  Override: set JULIA_SOBOL_FORCE_BIASED_RUN=1 to proceed anyway")
            println()
            println(
                @sprintf(
                    "  %-30s %12s %12s %8s  %s",
                    "QoI",
                    "D (var)",
                    "bad_var",
                    "% of D",
                    "status"
                )
            )
            println("  " * "-"^72)

            max_bad_frac = 0.0
            worst_qoi = ""
            worst_basis_per_qoi = Dict{String,Tuple{Int,Float64,Vector{Int}}}()
            sorted_qois = sort(collect(keys(all_coeffs)))
            for qname in sorted_qois
                c = all_coeffs[qname]
                D = sum(c[j]^2 for j = 2:length(c))
                bad_var = sum(c[j]^2 for j in failing_basis_idx; init = 0.0)
                bad_frac = D > 0 ? bad_var / D : 0.0

                if bad_frac > max_bad_frac
                    max_bad_frac = bad_frac
                    worst_qoi = qname
                end

                # Track the single worst basis function for this QoI's bad
                # contribution — useful for deciding whether the issue is a
                # single outlier or broadly distributed.
                worst_j_this_qoi = 0
                worst_c2_this_qoi = 0.0
                for j in failing_basis_idx
                    c2_j = c[j]^2
                    if c2_j > worst_c2_this_qoi
                        worst_c2_this_qoi = c2_j
                        worst_j_this_qoi = j
                    end
                end
                worst_basis_per_qoi[qname] = (
                    worst_j_this_qoi,
                    worst_c2_this_qoi,
                    worst_j_this_qoi == 0 ? Int[] : mop.ind[worst_j_this_qoi, :],
                )

                status = if bad_frac < 0.005
                    "✓ fine"
                elseif bad_frac < BIASED_BASIS_ABORT_THRESHOLD
                    "⚠ minor"
                else
                    "✗ BIAS RISK"
                end
                println(
                    @sprintf(
                        "  %-30s %12.4e %12.4e %7.3f%%  %s",
                        qname,
                        D,
                        bad_var,
                        100 * bad_frac,
                        status
                    )
                )
            end

            println()
            println(
                @sprintf(
                    "  Summary: worst QoI is %s with %.3f%% of variance from",
                    worst_qoi,
                    100 * max_bad_frac
                )
            )
            println("    basis functions that failed the 3σ orthonormality check.")

            # Report the worst individual bad-basis contribution for the worst QoI
            if haskey(worst_basis_per_qoi, worst_qoi)
                worst_j, worst_c2, worst_α = worst_basis_per_qoi[worst_qoi]
                if worst_j > 0
                    D_worst = sum(
                        all_coeffs[worst_qoi][j]^2 for j = 2:length(all_coeffs[worst_qoi])
                    )
                    frac_single = D_worst > 0 ? worst_c2 / D_worst : 0.0
                    emp_norm_worst = empirical_basis_norms[worst_j]
                    println(
                        @sprintf(
                            "    Single worst term in %s: basis %d (α=%s)",
                            worst_qoi,
                            worst_j,
                            worst_α
                        )
                    )
                    println(
                        @sprintf(
                            "      c²=%.4e, contributes %.3f%% of variance, E[Ψ²]=%.4f",
                            worst_c2,
                            100 * frac_single,
                            emp_norm_worst
                        )
                    )
                end
            end

            println()
            if max_bad_frac >= BIASED_BASIS_ABORT_THRESHOLD
                println(
                    "  ✗ BIAS RISK DETECTED — at least one QoI has >" *
                    @sprintf("%.1f%%", 100 * BIASED_BASIS_ABORT_THRESHOLD) *
                    " of its variance",
                )
                println(
                    "    coming from basis functions that failed the orthonormality check.",
                )
                println(
                    "    The bootstrap would spend hours producing CIs around biased point",
                )
                println("    estimates. Recommended actions (in order of preference):")
                println()
                println("      1. Force a lower PCE degree (e.g., PCE_DEGREE_RANGE = 3:3)")
                println(
                    "         and re-run. At degree 3, basis functions have lower variation",
                )
                println("         and the QMC sample covers their support much better.")
                println()
                println(
                    "      2. Increase N_SAMPLES in sensitivity_douglas_fir.jl from 8192 to",
                )
                println(
                    "         16384 or 32768 and regenerate parameters.csv / scalar_qois.csv.",
                )
                println(
                    "         This reduces the QMC finite-sample bias at the current degree.",
                )
                println()
                println(
                    "      3. Accept the bias and set JULIA_SOBOL_FORCE_BIASED_RUN=1. Do this",
                )
                println(
                    "         only if you intend to report degrees or QoIs where the bias is",
                )
                println("         documented as acceptable.")
                println()
                if get(ENV, "JULIA_SOBOL_FORCE_BIASED_RUN", "0") != "1"
                    error(
                        "Pre-bootstrap diagnostic: bias risk exceeds " *
                        "BIASED_BASIS_ABORT_THRESHOLD = $BIASED_BASIS_ABORT_THRESHOLD " *
                        "(worst QoI: $worst_qoi at $(round(100*max_bad_frac, digits=2))%). " *
                        "See recommendations above. To override, set " *
                        "JULIA_SOBOL_FORCE_BIASED_RUN=1 in the environment.",
                    )
                else
                    println(
                        "  ⚠ JULIA_SOBOL_FORCE_BIASED_RUN=1 is set — proceeding anyway.",
                    )
                    println(
                        "    Sobol indices for affected QoIs will be biased. Document this",
                    )
                    println("    limitation when reporting the results.")
                end
            elseif max_bad_frac >= 0.005
                println(
                    "  ⚠ MINOR BIAS: worst QoI has " *
                    @sprintf("%.3f%%", 100 * max_bad_frac) *
                    " of variance from biased basis,",
                )
                println(
                    "    below the $(round(100 * BIASED_BASIS_ABORT_THRESHOLD, digits=1))% abort threshold. " *
                    "Proceeding with the bootstrap.",
                )
                println("    Mention this limitation in the paper's methodology section.")
            else
                println(
                    "  ✓ All QoIs have <0.5% of their variance from biased basis functions.",
                )
                println(
                    "    The basis normalization warning is a cosmetic QMC finite-sample effect;",
                )
                println("    the Sobol indices and bootstrap CIs are robust. Proceeding.")
            end
        end

        # --- Step 4: Bootstrap confidence intervals ---
        println("\n" * "-"^70)
        println("  BCa Bootstrap Confidence Intervals ($N_BOOTSTRAP replicates)")
        println("-"^70)

        t_all_boot = time()
        n_valid_qois = length(valid_qois)
        for (qi_idx, qname) in enumerate(valid_qois)
            qi = findfirst(==(qname), ALL_QOI_NAMES)
            y = qoi_data[:, qi]

            t_this_qoi = time()
            println(
                @sprintf(
                    "  [%d/%d] %s: starting bootstrap + jackknife...",
                    qi_idx,
                    n_valid_qois,
                    qname
                )
            )
            flush(stdout)

            S1_med, S1_lo, S1_hi, ST_med, ST_lo, ST_hi = bootstrap_sobol_ci(
                Z,
                y,
                mop,
                GROUPS,
                N_BOOTSTRAP;
                ci_level = BOOTSTRAP_CI_LEVEL,
                use_sparse = use_sparse,
                Φ = Φ,
                c_hat = all_coeffs[qname],
                jack_Φ = all_jack_Φ[qname],
                jack_c = all_jack_c[qname],
                jack_h_diag = all_jack_h_diag[qname],
                jack_residuals = all_jack_residuals[qname],
                jack_support = all_jack_support[qname],
            )

            all_boot_ci[qname] = (
                S1_med = S1_med,
                S1_lo = S1_lo,
                S1_hi = S1_hi,
                ST_med = ST_med,
                ST_lo = ST_lo,
                ST_hi = ST_hi,
            )

            println(
                @sprintf(
                    "  [%d/%d] %s done [%.1fs wall, %.1fs total]",
                    qi_idx,
                    n_valid_qois,
                    qname,
                    time() - t_this_qoi,
                    time() - t_all_boot
                )
            )
            ci_pct = round(Int, BOOTSTRAP_CI_LEVEL * 100)
            for gi = 1:n_groups
                if all_S_first[qname][gi] > 0.01  # only show non-negligible
                    println(
                        @sprintf(
                            "    %-22s S₁=%.4f [%.4f, %.4f]  S_T=%.4f [%.4f, %.4f]  (%d%% CI)",
                            group_names[gi],
                            all_S_first[qname][gi],
                            S1_lo[gi],
                            S1_hi[gi],
                            all_S_total[qname][gi],
                            ST_lo[gi],
                            ST_hi[gi],
                            ci_pct
                        )
                    )
                end
            end
        end
        println(
            @sprintf(
                "  Bootstrap total: %d QoIs in %.1fs",
                n_valid_qois,
                time() - t_all_boot
            )
        )

        # --- Step 5: Convergence diagnostics ---
        println("\n" * "-"^70)
        println("  Convergence Diagnostics")
        println("-"^70)

        conv_results = Dict{String,Tuple{Vector{Int},Vector{Vector{Float64}}}}()

        for qname in valid_qois
            qi = findfirst(==(qname), ALL_QOI_NAMES)
            y = qoi_data[:, qi]

            sample_sizes, S1_by_size =
                convergence_analysis(Z, y, GROUPS, best_degree, CONVERGENCE_FRACTIONS)
            conv_results[qname] = (sample_sizes, S1_by_size)

            # Robust multi-criterion convergence check
            if length(S1_by_size) >= 2
                conv_diag = assess_convergence(S1_by_size, sample_sizes)
                if conv_diag.converged
                    println(
                        @sprintf(
                            "  %s: ✓ converged (max terminal |ΔS₁| = %.4f)",
                            qname,
                            conv_diag.max_terminal_change
                        )
                    )
                else
                    issues = String[]
                    conv_diag.max_terminal_change >= 0.02 && push!(
                        issues,
                        @sprintf("terminal |ΔS₁|=%.4f", conv_diag.max_terminal_change)
                    )
                    conv_diag.n_oscillating > 0 &&
                        push!(issues, "$(conv_diag.n_oscillating) groups oscillating")
                    conv_diag.n_drifting > 0 &&
                        push!(issues, "$(conv_diag.n_drifting) groups drifting")
                    println(
                        @sprintf("  %s: ⚠ not converged — %s", qname, join(issues, ", "))
                    )
                end
            end
        end

        # ====================================================================
        # Save results
        # ====================================================================

        # --- Group Sobol indices CSV ---
        sobol_file = joinpath(RESULTS_DIR, "sobol_group_indices.csv")
        open(sobol_file, "w") do io
            qoi_keys = sort(collect(keys(all_S_first)))
            # Header
            header_parts = ["group"]
            for q in qoi_keys
                push!(header_parts, "$(q)_S1", "$(q)_ST")
            end
            println(io, join(header_parts, ","))
            # Data
            for (gi, gname) in enumerate(group_names)
                parts = ["\"$gname\""]
                for q in qoi_keys
                    push!(parts, @sprintf("%.6f", all_S_first[q][gi]))
                    push!(parts, @sprintf("%.6f", all_S_total[q][gi]))
                end
                println(io, join(parts, ","))
            end
        end
        println("\n  Saved $sobol_file")

        # --- Bootstrap CI CSV ---
        boot_file = joinpath(RESULTS_DIR, "sobol_bootstrap_ci.csv")
        open(boot_file, "w") do io
            ci_pct = round(Int, BOOTSTRAP_CI_LEVEL * 100)
            println(
                io,
                "qoi,group,S1,S1_lo_$(ci_pct),S1_hi_$(ci_pct),ST,ST_lo_$(ci_pct),ST_hi_$(ci_pct)",
            )
            for qname in sort(collect(keys(all_boot_ci)))
                ci = all_boot_ci[qname]
                for (gi, gname) in enumerate(group_names)
                    println(
                        io,
                        @sprintf(
                            "%s,\"%s\",%.6f,%.6f,%.6f,%.6f,%.6f,%.6f",
                            qname,
                            gname,
                            all_S_first[qname][gi],
                            ci.S1_lo[gi],
                            ci.S1_hi[gi],
                            all_S_total[qname][gi],
                            ci.ST_lo[gi],
                            ci.ST_hi[gi]
                        )
                    )
                end
            end
        end
        println("  Saved $boot_file")

        # --- PCE validation CSV ---
        validation_file = joinpath(RESULTS_DIR, "pce_validation.csv")
        open(validation_file, "w") do io
            println(
                io,
                "qoi,degree,fit_method,R2,Q2_LOO,n_basis,n_samples,oversampling_ratio",
            )
            fit_method = use_sparse ? "LASSO" : "OLS"
            for qname in sort(collect(keys(all_R2)))
                println(
                    io,
                    @sprintf(
                        "%s,%d,%s,%.6f,%.6f,%d,%d,%.1f",
                        qname,
                        best_degree,
                        fit_method,
                        all_R2[qname],
                        all_Q2[qname],
                        n_basis,
                        n_success,
                        n_success / n_basis
                    )
                )
            end
        end
        println("  Saved $validation_file")

        # --- Per-sample PCE predictions CSV (for validation scatter plots) ---
        pred_file = joinpath(RESULTS_DIR, "pce_predictions.csv")
        open(pred_file, "w") do io
            pred_qoi_keys = sort(collect(keys(all_predictions)))
            header =
                join(["sample_id"; ["$(q)_true,$(q)_pred" for q in pred_qoi_keys]], ",")
            println(io, header)
            for i = 1:n_success
                parts = [string(sample_ids[i])]
                for q in pred_qoi_keys
                    qi = findfirst(==(q), ALL_QOI_NAMES)
                    push!(parts, @sprintf("%.10e", qoi_data[i, qi]))
                    push!(parts, @sprintf("%.10e", all_predictions[q][i]))
                end
                println(io, join(parts, ","))
            end
        end
        println("  Saved $pred_file")

        # --- PCE coefficients CSV ---
        coeff_file = joinpath(RESULTS_DIR, "pce_coefficients.csv")
        open(coeff_file, "w") do io
            qoi_keys = sort(collect(keys(all_coeffs)))
            # Header: multi-index columns + one column per QoI
            idx_header = join(["alpha_$i" for i = 1:N_PARAMS], ",")
            qoi_header = join(qoi_keys, ",")
            println(io, "basis_idx,$idx_header,$qoi_header")
            for j = 1:n_basis
                α_str = join(string.(mop.ind[j, :]), ",")
                c_str = join([@sprintf("%.10e", all_coeffs[q][j]) for q in qoi_keys], ",")
                println(io, "$j,$α_str,$c_str")
            end
        end
        println("  Saved $coeff_file")

        # --- Convergence data CSV (Sobol indices vs subsample size) ---
        conv_file = joinpath(RESULTS_DIR, "sobol_convergence_data.csv")
        open(conv_file, "w") do io
            println(io, "qoi,subsample_size,group,S1")
            for qname in sort(collect(keys(conv_results)))
                sample_sizes, S1_by_size = conv_results[qname]
                for (si, ns) in enumerate(sample_sizes)
                    for (gi, gname) in enumerate(group_names)
                        println(
                            io,
                            @sprintf(
                                "%s,%d,\"%s\",%.6f",
                                qname,
                                ns,
                                gname,
                                S1_by_size[si][gi]
                            )
                        )
                    end
                end
            end
        end
        println("  Saved $conv_file")

        # ====================================================================
        # Plots
        # ====================================================================

        println("\nGenerating plots...")

        qoi_keys = sort(collect(keys(all_S_first)))

        # --- 1. Stacked bar chart: variance decomposition across all QoIs ---
        qoi_labels = [replace(q, "_" => " ") for q in qoi_keys]
        n_qois = length(qoi_keys)

        # Build matrix [n_groups × n_qois] for first-order indices
        S_matrix = zeros(n_groups, n_qois)
        for (qi, q) in enumerate(qoi_keys)
            S_matrix[:, qi] = all_S_first[q]
        end

        # Add interaction row
        interaction_row = [1.0 - sum(all_S_first[q]) for q in qoi_keys]

        # Sort groups by average importance
        avg_importance = mean(S_matrix, dims = 2)[:]
        sort_order = sortperm(avg_importance, rev = true)

        # Color palette for groups
        colors = [
            :steelblue,
            :coral,
            :seagreen,
            :mediumpurple,
            :goldenrod,
            :lightblue,
            :pink,
            :sienna,
        ]

        # Build stacked bar chart using areaplot-style manual stacking
        plt_stacked = plot(
            ylabel = "Fraction of variance (S₁)",
            title = "Variance Decomposition by Measurement Source",
            ylims = (0, 1.08),
            legend = :outertopright,
            size = (1050, 500),
            bottom_margin = 12Plots.mm,
            left_margin = 5Plots.mm,
            xticks = (1:n_qois, qoi_labels),
            xrotation = 25,
        )

        # Stack rectangles from top to bottom so legend order matches visual order
        cumulative = ones(n_qois)  # start from the top
        for (rank, gi) in enumerate(sort_order)
            vals = S_matrix[gi, :]
            new_cumulative = cumulative .- vals
            bar!(
                1:n_qois,
                cumulative,
                fillrange = new_cumulative,
                label = group_names[gi],
                color = colors[rank],
                alpha = 0.85,
                lw = 0.5,
                linecolor = :white,
            )
            cumulative = new_cumulative
        end
        hline!([1.0], color = :black, lw = 0.5, ls = :dash, label = false)

        savefig(plt_stacked, joinpath(RESULTS_DIR, "sobol_variance_decomposition.png"))
        println("  Saved sobol_variance_decomposition.png")

        # --- 2. Individual bar charts per QoI (with bootstrap CI error bars) ---
        for qname in qoi_keys
            S1 = all_S_first[qname]
            ST = all_S_total[qname]
            sorted = sortperm(S1, rev = true)

            plt = bar(
                1:n_groups,
                S1[sorted],
                label = "First-order (S₁)",
                color = :steelblue,
                alpha = 0.85,
                yticks = (1:n_groups, group_names[sorted]),
                orientation = :horizontal,
                xlims = (0, 1.05),
                xlabel = "Sobol Index",
                title = "Variance Decomposition — $(replace(qname, "_" => " "))",
                size = (800, 450),
                left_margin = 18Plots.mm,
            )
            bar!(
                1:n_groups,
                ST[sorted],
                label = "Total (S_T)",
                color = :coral,
                alpha = 0.4,
                orientation = :horizontal,
            )

            # Add bootstrap CI whiskers if available
            if haskey(all_boot_ci, qname)
                ci = all_boot_ci[qname]
                for (rank, gi) in enumerate(sorted)
                    # S1 CI
                    plot!(
                        [ci.S1_lo[gi], ci.S1_hi[gi]],
                        [rank, rank],
                        color = :black,
                        lw = 1.5,
                        label = (rank == 1 ? "95% CI" : false),
                    )
                end
            end

            savefig(plt, joinpath(RESULTS_DIR, "sobol_bar_$(qname).png"))
            println("  Saved sobol_bar_$(qname).png")
        end

        # --- 3. Validation scatter plots ---
        for qname in qoi_keys
            qi = findfirst(==(qname), ALL_QOI_NAMES)
            y = qoi_data[:, qi]
            y_pred = all_predictions[qname]
            R2 = all_R2[qname]
            Q2 = all_Q2[qname]

            plt = scatter(
                y,
                y_pred,
                xlabel = "Simulation",
                ylabel = "PCE prediction",
                title = @sprintf(
                    "%s (R²=%.4f, Q²=%.4f)",
                    replace(qname, "_" => " "),
                    R2,
                    Q2
                ),
                legend = false,
                alpha = 0.3,
                ms = 2,
                color = :steelblue,
                size = (500, 500),
                aspect_ratio = :equal,
            )
            # Perfect prediction line
            yl = extrema(y)
            plot!([yl[1], yl[2]], [yl[1], yl[2]], color = :red, lw = 1.5, ls = :dash)

            savefig(plt, joinpath(RESULTS_DIR, "pce_validation_$(qname).png"))
            println("  Saved pce_validation_$(qname).png")
        end

        # --- 4. Convergence plots ---
        for qname in qoi_keys
            if !haskey(conv_results, qname)
                continue
            end
            sample_sizes, S1_by_size = conv_results[qname]

            # Build matrix [n_sizes × n_groups]
            n_sizes = length(sample_sizes)
            conv_matrix = zeros(n_sizes, n_groups)
            for i = 1:n_sizes
                conv_matrix[i, :] = S1_by_size[i]
            end

            # Only plot groups with S1 > 0.02 at the final size
            important_groups = findall(S1_by_size[end] .> 0.02)

            plt = plot(
                xlabel = "Number of samples",
                ylabel = "First-order Sobol index (S₁)",
                title = "Convergence — $(replace(qname, "_" => " "))",
                size = (700, 450),
                legend = :outerright,
            )

            for gi in important_groups
                plot!(
                    sample_sizes,
                    conv_matrix[:, gi],
                    label = group_names[gi],
                    lw = 2,
                    marker = :circle,
                    ms = 4,
                    color = colors[mod1(gi, length(colors))],
                )
            end

            savefig(plt, joinpath(RESULTS_DIR, "sobol_convergence_$(qname).png"))
            println("  Saved sobol_convergence_$(qname).png")
        end

        # ====================================================================
        # Summary
        # ====================================================================

        println("\n" * "="^70)
        println("  Summary")
        println("="^70)

        println(
            @sprintf("\n  PCE degree: %d (%s)", best_degree, use_sparse ? "LASSO" : "OLS")
        )

        println("\n  PCE fit quality:")
        for qname in qoi_keys
            status = all_R2[qname] > 0.95 && all_Q2[qname] > 0.90 ? "✓" : "⚠"
            println(
                @sprintf(
                    "    %s %-25s R²=%.4f  Q²_LOO=%.4f",
                    status,
                    qname,
                    all_R2[qname],
                    all_Q2[qname]
                )
            )
        end

        ci_pct = round(Int, BOOTSTRAP_CI_LEVEL * 100)
        println("\n  Dominant measurement sources (S₁ > 0.05, with $ci_pct% bootstrap CI):")
        for qname in qoi_keys
            S1 = all_S_first[qname]
            important = [(gi, group_names[gi], S1[gi]) for gi = 1:n_groups if S1[gi] > 0.05]
            sort!(important, by = x -> -x[3])
            if haskey(all_boot_ci, qname)
                ci = all_boot_ci[qname]
                parts = [
                    @sprintf(
                        "%s %.0f%% [%.0f, %.0f]",
                        name,
                        s * 100,
                        ci.S1_lo[gi] * 100,
                        ci.S1_hi[gi] * 100
                    ) for (gi, name, s) in important
                ]
            else
                parts =
                    [@sprintf("%s (%.0f%%)", name, s * 100) for (_, name, s) in important]
            end
            println("    $(qname): $(join(parts, ", "))")
        end

        println("\n  Results in: $RESULTS_DIR")
        println("="^70)
    finally
        # Restore BLAS thread count so that `include`-ing this file doesn't
        # leave a global side effect behind for the caller.
        BLAS.set_num_threads(blas_threads_before)
    end
end

main()
