"""
Group Sobol Sensitivity Indices via Polynomial Chaos Expansion

Computes variance-based sensitivity indices for the Douglas Fir cone calorimeter
simulation, decomposing output variance into contributions from each measurement
source (instrument group).

Method:
  1. Transform correlated parameter samples to independent standard normals
	 using the Cholesky factors of each group's covariance matrix.
  2. Fit an adaptive-degree (2–4) Polynomial Chaos Expansion (PCE) surrogate to
	 each scalar QoI. Degree is selected by maximizing leave-one-out Q². For
	 degree ≥ 3 where the oversampling ratio is low, a sparse LASSO fit
	 (coordinate descent with 5-fold CV) is used instead of OLS.
  3. Compute first-order and total group Sobol indices analytically from the
	 PCE coefficients.
  4. Quantify index uncertainty via nonparametric bootstrap (1000 replicates).
  5. Assess convergence by computing indices on nested subsamples.

The 17 uncertain parameters fall into 8 independent groups (correlations are
within-group only), so standard Sobol formulas apply without correlated-input
extensions.

Groups (by measurement source):
  1. HFM calibration bias (ε_HFM)           — 1 parameter
  2. DSC calibration bias (ε_DSC)           — 1 parameter
  3. HoS calibration bias (ε_HoS)          — 1 parameter
  4. Virgin Cp, low-T (HFM regression)      — 2 parameters (a₁, b₁)
  5. Virgin Cp, high-T (STA per-curve)      — 2 parameters (a₃, b₃)
  6. Char Cp (independent char STA)         — 2 parameters (a_char, b_char)
  7. Heat of sorption (CR2 cov, cond.)      — 3 parameters (a, b, c)
  8. Heats of reaction (conditional cov.)   — 5 parameters (ΔH₁–ΔH₅)

Reads:  results/parameters.csv, results/scalar_qois.csv
Writes: results/sobol_group_indices.csv, results/sobol_bootstrap_ci.csv,
		results/pce_validation.csv, results/pce_coefficients.csv,
		results/sobol_variance_decomposition.png,
		results/sobol_bar_<qoi>.png, results/pce_validation_<qoi>.png,
		results/sobol_convergence_<qoi>.png

Usage:  julia --project=. examples/sensitivity_analysis/sensitivity_sobol.jl
"""

using PolyChaos
using LinearAlgebra
using Statistics
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
	indices::Vector{Int}     # indices into the 17-element θ vector
	μ::Vector{Float64}
	Σ::Matrix{Float64}
end

const GROUPS = [
	# 1. HFM calibration bias (data-driven, SRM 1453/1453e, DF thickness)
	ParamGroup("HFM calibration", [1],
		[0.0],
		[0.0046154243860687^2;;]),

	# 2. DSC calibration bias (certified-only RMS, 1.13%, dof=20)
	ParamGroup("DSC calibration", [2],
		[0.0],
		[0.011322211895623097^2;;]),

	# 3. HoS asymptote calibration bias (u_rel = 0.60%)
	ParamGroup("HoS calibration", [3],
		[0.0],
		[0.005974653346763125^2;;]),

	# 4. Virgin Cp, low-T (HFM regression, slope/intercept in J/(kg·K))
	# Source: douglas_fir_hfm/regression_results.json
	ParamGroup("Virgin Cp, low-T", [4, 5],
		[0.004133836600273531 * 1e3, -0.01801630709338176 * 1e3],
		[7.183142841163425e-05^2 * 1e6    -1.1925585013168409e-06 * 1e6;
		 -1.1925585013168409e-06 * 1e6    0.01850672773953015^2 * 1e6]),

	# 5. Virgin Cp, high-T (STA per-curve, slope/intercept in J/(kg·K))
	# Source: gsa_parameter_config_douglas_fir.json
	ParamGroup("Virgin Cp, high-T", [6, 7],
		[0.009136799235573225 * 1e3, -2.1393044552435705 * 1e3],
		[2.7703324106603792e-08 * 1e6    -9.7846678915291e-06 * 1e6;
		 -9.7846678915291e-06 * 1e6       0.0038340012469450084 * 1e6]),

	# 6. Char Cp (independent STA, slope/intercept in J/(kg·K))
	# Source: douglas_fir_independent_char_sta/regression_results.json
	ParamGroup("Char Cp", [8, 9],
		[0.0018372509265710314 * 1e3, 0.3124273785167445 * 1e3],
		[0.00040051156523772673^2 * 1e6    -9.472629865704193e-05 * 1e6;
		 -9.472629865704193e-05 * 1e6       0.2369805093355783^2 * 1e6]),

	# 7. Heat of sorption (CR2 covariance, conditional on f_cal exact)
	# Source: papers/energetics/figures/hos_fit_results.json -> covariance_matrix_cr2
	ParamGroup("Heat of sorption", [10, 11, 12],
		[2450.420909987441, 875.6441221735336, 0.1121684273316275],
		[255.54778632642243     826.4091656069307     0.24363365219954056;
		 826.4091656069309     4539.412144899307      1.054987399933727;
		 0.24363365219954064    1.0549873999337271    0.00028091975868012163]),

	# 8. Heats of reaction (conditional covariance, kJ²/kg² -> J²/kg²)
	# Source: gsa_parameter_config_douglas_fir.json -> conditional_covariance
	ParamGroup("Heats of reaction", [13, 14, 15, 16, 17],
		[-19.20066288032251 * 1e3, -5.263201466602739 * 1e3, 264.42479037373835 * 1e3, -217.15169043762444 * 1e3, -0.32640885469229025 * 1e3],
		[12.188429198060906    24.748705356848376   29.002904740480126    88.4320951009737     -6.604044259616362;
		 24.748705356848376   50.700350890015784    59.67848721837838    184.93304070273382   -13.71270702182141;
		 29.002904740480123   59.67848721837838     80.78981196904051    246.88080651001275    -7.889710589815444;
		 88.4320951009737    184.93304070273382     246.88080651001275    837.274115396477       0.9278934423897951;
		 -6.604044259616367  -13.712707021821425    -7.889710589815465     0.9278934423897454  41.07609754085615] .* 1e6),
]

const N_PARAMS = 17

# Scalar QoIs to analyze (skip time_to_ignition — quantized, and solve_time — not physical)
const QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr",
	"total_mass_loss_pct", "t_surf_70pct", "t_surf_95pct",
	"t_back_5pct", "t_back_95pct",
	"final_thickness",
	"peak_mlr_first", "time_to_peak_mlr_first",
	"pyrolysis_onset_time", "T_back_50pct_mass_loss",
]

# Column indices in scalar_qois.csv (after sample_id and retcode columns)
const ALL_QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr", "time_to_ignition",
	"total_mass_loss_pct", "t_surf_70pct", "t_surf_95pct",
	"t_back_5pct", "t_back_95pct",
	"final_thickness",
	"peak_mlr_first", "time_to_peak_mlr_first",
	"pyrolysis_onset_time", "T_back_50pct_mass_loss",
	"solve_time",
]

# PCE configuration
const PCE_DEGREE_RANGE = 2:4      # Adaptive: try degrees 2, 3, 4; pick best Q²
const SPARSE_THRESHOLD = 3.0      # Use sparse (LASSO) fit when oversampling ratio < this
const N_BOOTSTRAP = 1000          # Bootstrap replicates for Sobol CI
const BOOTSTRAP_CI_LEVEL = 0.95   # Confidence level for bootstrap intervals
# Subsample sizes for convergence diagnostics (powers of 2 up to full sample)
const CONVERGENCE_FRACTIONS = [0.125, 0.25, 0.5, 0.75, 1.0]

# ============================================================================
# Step 1: Transform to independent standard normals
# ============================================================================

"""
	transform_to_standard_normal(param_data, groups)

Transform correlated parameter samples to independent standard normal space
using the Cholesky factors of each group's covariance matrix.

Input:  param_data [n_samples × 17] in physical units
Output: Z [n_samples × 17] in standard normal space
"""
function transform_to_standard_normal(param_data::Matrix{Float64}, groups::Vector{ParamGroup})
	n = size(param_data, 1)
	Z = similar(param_data)

	for g in groups
		idx = g.indices
		L = cholesky(Symmetric(g.Σ)).L
		for i in 1:n
			Z[i, idx] = L \ (param_data[i, idx] - g.μ)
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
	for d in 1:n_dims
		uni_vals[d] = Vector{Vector{Float64}}(undef, deg + 1)
		for k in 0:deg
			# Evaluate degree-k univariate orthogonal poly at all Z[:,d]
			uni_vals[d][k+1] = PolyChaos.evaluate(k, Z[:, d], mop.uni[d])
		end
	end

	# --- Assemble multivariate basis via tensor products ---
	Φ = Matrix{Float64}(undef, n_samples, n_basis)
	for j in 1:n_basis
		α = mop.ind[j, :]  # multi-index: degree per dimension

		# Product of univariate evaluations
		col = ones(n_samples)
		for d in 1:n_dims
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
	verify_basis_normalization(Φ, mop; tol=0.15)

Empirically verify that the orthonormalized PCE basis satisfies E[Ψ_j²] ≈ 1
(using sample mean as proxy). This catches normalization errors that would
bias the Sobol variance decomposition toward lower-order terms.

Returns (max_deviation, pass) where pass is true if all basis functions
have sample E[Ψ²] within `tol` of 1.0.
"""
function verify_basis_normalization(Φ::Matrix{Float64}, mop::MultiOrthoPoly; tol::Float64 = 0.15)
	n_basis = size(Φ, 2)
	empirical_norms = vec(mean(Φ .^ 2, dims = 1))

	max_dev = 0.0
	n_fail = 0
	worst_idx = 1
	for j in 2:n_basis  # skip constant term (always = 1)
		dev = abs(empirical_norms[j] - 1.0)
		if dev > max_dev
			max_dev = dev
			worst_idx = j
		end
		if dev > tol
			n_fail += 1
		end
	end

	pass = n_fail == 0
	if !pass
		α_worst = mop.ind[worst_idx, :]
		@warn "Basis normalization check FAILED: $n_fail / $(n_basis-1) basis functions " *
			  "deviate from E[Ψ²]=1 by more than $tol. " *
			  "Worst: basis $worst_idx (α=$α_worst), E[Ψ²]=$(round(empirical_norms[worst_idx], digits=4)). " *
			  "Sobol indices may be biased."
	end

	return max_dev, pass
end

"""
	fit_pce(Φ, y)

Fit PCE coefficients via ordinary least squares.
Returns coefficients c, R², and leave-one-out Q².
"""
function fit_pce(Φ::Matrix{Float64}, y::Vector{Float64})
	n, p = size(Φ)

	# OLS fit
	c = Φ \ y
	y_pred = Φ * c
	residuals = y - y_pred

	# R²
	ss_res = sum(residuals .^ 2)
	ss_tot = sum((y .- mean(y)) .^ 2)
	R2 = 1.0 - ss_res / ss_tot

	# Leave-one-out Q² via analytical formula (no refitting needed)
	# Q² = 1 - (1/n) Σ [eᵢ/(1-hᵢᵢ)]² / Var(y)
	H = Φ * ((Φ' * Φ) \ Φ')  # hat matrix
	h_diag = diag(H)
	loo_residuals = residuals ./ (1.0 .- h_diag)
	ss_loo = sum(loo_residuals .^ 2)
	Q2 = 1.0 - ss_loo / ss_tot

	return c, R2, Q2, y_pred
end

# ============================================================================
# Step 3: Compute group Sobol indices from PCE coefficients
# ============================================================================

"""
	compute_group_sobol(c, mop, groups)

Compute first-order (S_g) and total (S_Tg) group Sobol indices from PCE
coefficients. Uses the orthonormality of the basis: D = Σ_{α≠0} c_α².
"""
function compute_group_sobol(c::Vector{Float64}, mop::MultiOrthoPoly,
	groups::Vector{ParamGroup})

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
	D = sum(c[j]^2 for j in 2:n_basis)

	if D ≈ 0
		return zeros(n_groups), zeros(n_groups), 0.0
	end

	S_first = zeros(n_groups)   # first-order group indices
	S_total = zeros(n_groups)   # total group indices

	for j in 2:n_basis  # skip constant term
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

"""
	fit_pce_lasso(Φ, y; n_lambda=50, n_folds=5, lambda_ratio=1e-4)

Fit sparse PCE coefficients via LASSO (L₁ regularization) using coordinate
descent, with K-fold cross-validation for regularization parameter selection.

Used when the oversampling ratio (n_samples / n_basis) is small (< SPARSE_THRESHOLD),
which occurs at higher PCE degrees. The LASSO automatically selects the most
important basis terms and zeros out the rest.

Performance: Uses the standard residual-tracking trick so each coordinate update
is O(n) rather than O(n×p), giving O(n×p) per full pass. Solutions are warm-started
along the λ path (previous λ's solution initializes the next), which typically
converges in a few iterations per λ step.

Reference: Blatman & Sudret (2011), "Adaptive sparse polynomial chaos expansion
based on least angle regression", J. Comput. Phys.

Returns coefficients c, R², Q² (from CV), and predictions y_pred.
"""
function fit_pce_lasso(Φ::Matrix{Float64}, y::Vector{Float64};
	n_lambda::Int = 50, n_folds::Int = 5, lambda_ratio::Float64 = 1e-4,
	max_iter::Int = 10000, tol::Float64 = 1e-7, verbose::Bool = false)

	t_start = time()
	n, p = size(Φ)
	verbose && println("      LASSO: n=$n samples, p=$p basis functions, $n_folds-fold CV over $n_lambda λ values")

	# Center y (don't regularize the intercept — but our Φ[:,1] is the constant
	# term, so we exclude it from penalization)
	y_mean = mean(y)
	y_centered = y .- y_mean

	# Standardize columns of Φ (except column 1 = constant)
	col_norms = zeros(p)
	Φ_scaled = copy(Φ)
	for j in 2:p
		col_norms[j] = norm(Φ_scaled[:, j]) / sqrt(n)
		if col_norms[j] > 0
			Φ_scaled[:, j] ./= col_norms[j]
		end
	end

	# Precompute column squared norms (for denominator of coordinate updates)
	col_sq_norms = vec(sum(Φ_scaled .^ 2, dims = 1))  # length p

	# Lambda path: from lambda_max (smallest λ giving all zeros) down
	# lambda_max = max_j |Φ_j' y| / n
	lambda_max = maximum(abs.(Φ_scaled[:, 2:end]' * y_centered)) / n
	lambda_min = lambda_max * lambda_ratio
	lambdas = exp.(range(log(lambda_max), log(lambda_min), length = n_lambda))

	# K-fold CV split
	rng = MersenneTwister(42)
	fold_ids = shuffle(rng, repeat(1:n_folds, ceil(Int, n / n_folds))[1:n])

	cv_mse = zeros(n_lambda)

	t_cv_start = time()
	for (li, λ) in enumerate(lambdas)
		if verbose && (li == 1 || li % 10 == 0 || li == n_lambda)
			pct = round(100 * li / n_lambda, digits = 0)
			elapsed = round(time() - t_cv_start, digits = 1)
			print("\r      CV: λ $li/$n_lambda ($(pct)%) [$(elapsed)s]")
		end
		fold_errors = zeros(n_folds)
		for fold in 1:n_folds
			train_idx = findall(fold_ids .!= fold)
			test_idx = findall(fold_ids .== fold)
			Φ_train = Φ_scaled[train_idx, :]
			y_train = y_centered[train_idx]
			Φ_test = Φ_scaled[test_idx, :]
			y_test = y_centered[test_idx]
			n_train = length(train_idx)

			# Precompute train column squared norms
			col_sq_train = vec(sum(Φ_train .^ 2, dims = 1))

			# Coordinate descent with residual tracking
			# Warm start from previous lambda (li > 1) would require
			# per-fold state — not worth the complexity for CV. Start from zero.
			c_fold = zeros(p)
			r = copy(y_train)  # residual r = y - Φ*c, initially c=0 so r=y

			for _ in 1:max_iter
				max_change = 0.0
				for j in 2:p  # skip intercept (column 1)
					# Partial residual for coordinate j: r_j = r + Φ[:,j]*c[j]
					# z_j = Φ[:,j]'*r_j / n_train = (Φ[:,j]'*r + col_sq_train[j]*c[j]) / n_train
					z_j = (dot(Φ_train[:, j], r) + col_sq_train[j] * c_fold[j]) / n_train
					c_new = soft_threshold(z_j, λ)
					if c_new != c_fold[j]
						delta = c_fold[j] - c_new
						# Update residual: r += Φ[:,j] * (c_old - c_new)
						@inbounds for i in eachindex(r)
							r[i] += Φ_train[i, j] * delta
						end
						max_change = max(max_change, abs(delta))
						c_fold[j] = c_new
					end
				end
				# Unpenalized intercept-like update for column 1
				z_1 = (dot(Φ_train[:, 1], r) + col_sq_train[1] * c_fold[1]) / col_sq_train[1]
				if z_1 != c_fold[1]
					delta = c_fold[1] - z_1
					@inbounds for i in eachindex(r)
						r[i] += Φ_train[i, 1] * delta
					end
					max_change = max(max_change, abs(delta))
					c_fold[1] = z_1
				end

				if max_change < tol
					break
				end
			end

			y_pred_test = Φ_test * c_fold
			fold_errors[fold] = mean((y_test .- y_pred_test) .^ 2)
		end
		cv_mse[li] = mean(fold_errors)
	end

	# Select lambda with minimum CV MSE
	best_li = argmin(cv_mse)
	λ_best = lambdas[best_li]

	# Refit on full data with best lambda, warm-started along λ path
	c = zeros(p)
	r_full = copy(y_centered)  # residual for full-data refit
	for li in 1:best_li
		λ_curr = lambdas[li]
		for _ in 1:max_iter
			max_change = 0.0
			for j in 2:p
				z_j = (dot(Φ_scaled[:, j], r_full) + col_sq_norms[j] * c[j]) / n
				c_new = soft_threshold(z_j, λ_curr)
				if c_new != c[j]
					delta = c[j] - c_new
					@inbounds for i in eachindex(r_full)
						r_full[i] += Φ_scaled[i, j] * delta
					end
					max_change = max(max_change, abs(delta))
					c[j] = c_new
				end
			end
			# Unpenalized intercept
			z_1 = (dot(Φ_scaled[:, 1], r_full) + col_sq_norms[1] * c[1]) / col_sq_norms[1]
			if z_1 != c[1]
				delta = c[1] - z_1
				@inbounds for i in eachindex(r_full)
					r_full[i] += Φ_scaled[i, 1] * delta
				end
				c[1] = z_1
			end

			if max_change < tol
				break
			end
		end
	end

	# Rescale coefficients back to original Φ scale
	for j in 2:p
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
	R2 = 1.0 - ss_res / ss_tot

	# Q² approximation from CV
	Q2 = 1.0 - cv_mse[best_li] * n / ss_tot

	n_nonzero = count(abs.(c[2:end]) .> 1e-12)
	t_total = round(time() - t_start, digits = 1)
	verbose && print("\r" * " "^60 * "\r")  # clear progress line
	println(@sprintf("      LASSO: λ=%.4e, %d/%d nonzero, R²=%.4f, Q²=%.4f  [%.1fs]", λ_best, n_nonzero, p - 1, R2, Q2, t_total))

	return c, R2, Q2, y_pred
end

# ============================================================================
# Step 4: Bootstrap Confidence Intervals for Sobol Indices
# ============================================================================

"""
	bootstrap_sobol_ci(Z, y, mop, groups, n_boot; ci_level=0.95, use_sparse=false)

Compute bootstrap confidence intervals for first-order and total group Sobol
indices by resampling (with replacement) the (Z, y) pairs, refitting the PCE,
and recomputing the indices.

Returns (S1_median, S1_lo, S1_hi, ST_median, ST_lo, ST_hi) each of length n_groups.
"""
function bootstrap_sobol_ci(Z::Matrix{Float64}, y::Vector{Float64},
	mop::MultiOrthoPoly, groups::Vector{ParamGroup}, n_boot::Int;
	ci_level::Float64 = 0.95, use_sparse::Bool = false)

	n = length(y)
	n_groups = length(groups)
	α_lo = (1.0 - ci_level) / 2.0
	α_hi = 1.0 - α_lo

	boot_S1 = Matrix{Float64}(undef, n_boot, n_groups)
	boot_ST = Matrix{Float64}(undef, n_boot, n_groups)

	rng = MersenneTwister(123)

	t_boot = time()
	for b in 1:n_boot
		if b == 1 || b % 50 == 0 || b == n_boot
			elapsed = round(time() - t_boot, digits = 1)
			if b > 1
				rate = (b - 1) / (time() - t_boot)
				eta = round((n_boot - b) / rate, digits = 0)
				print(@sprintf("\r    bootstrap %d/%d [%.1fs elapsed, ~%.0fs remaining]", b, n_boot, elapsed, eta))
			else
				print(@sprintf("\r    bootstrap %d/%d ...", b, n_boot))
			end
		end

		# Resample with replacement
		idx = rand(rng, 1:n, n)
		Z_b = Z[idx, :]
		y_b = y[idx]

		# Rebuild design matrix for bootstrap sample
		Φ_b = build_design_matrix(Z_b, mop)

		# Fit PCE
		if use_sparse
			c_b, _, _, _ = fit_pce_lasso(Φ_b, y_b)
		else
			c_b = Φ_b \ y_b  # fast OLS without diagnostics
		end

		S1_b, ST_b, _ = compute_group_sobol(c_b, mop, groups)
		boot_S1[b, :] = S1_b
		boot_ST[b, :] = ST_b
	end
	print("\r" * " "^70 * "\r")  # clear progress line

	S1_median = vec(median(boot_S1, dims = 1))
	S1_lo = [quantile(boot_S1[:, g], α_lo) for g in 1:n_groups]
	S1_hi = [quantile(boot_S1[:, g], α_hi) for g in 1:n_groups]
	ST_median = vec(median(boot_ST, dims = 1))
	ST_lo = [quantile(boot_ST[:, g], α_lo) for g in 1:n_groups]
	ST_hi = [quantile(boot_ST[:, g], α_hi) for g in 1:n_groups]

	return S1_median, S1_lo, S1_hi, ST_median, ST_lo, ST_hi
end

# ============================================================================
# Step 5: Convergence Diagnostics
# ============================================================================

"""
	convergence_analysis(Z, y, groups, degree, fractions; use_sparse=false)

Compute Sobol indices on nested subsamples of the data to assess convergence.
`fractions` specifies the fraction of the full dataset to use at each level.

Returns (sample_sizes, S1_by_size) where S1_by_size[i] is the n_groups-length
first-order index vector at subsample size sample_sizes[i].
"""
function convergence_analysis(Z::Matrix{Float64}, y::Vector{Float64},
	groups::Vector{ParamGroup}, degree::Int, fractions::Vector{Float64};
	use_sparse::Bool = false)

	n = length(y)
	n_groups = length(groups)
	sample_sizes = [round(Int, f * n) for f in fractions]

	S1_by_size = Vector{Vector{Float64}}(undef, length(sample_sizes))

	for (i, ns) in enumerate(sample_sizes)
		print(@sprintf("\r    convergence: subsample %d/%d (n=%d)", i, length(sample_sizes), ns))

		# Use first ns samples (preserving Sobol ordering)
		Z_sub = Z[1:ns, :]
		y_sub = y[1:ns]

		# Build PCE at given degree
		ops = [GaussOrthoPoly(degree) for _ in 1:N_PARAMS]
		mop = MultiOrthoPoly(ops, degree)
		Φ = build_design_matrix(Z_sub, mop)

		oversampling = ns / mop.dim
		if use_sparse && oversampling < SPARSE_THRESHOLD
			c, _, _, _ = fit_pce_lasso(Φ, y_sub)
		else
			c = Φ \ y_sub
		end

		S1, _, _ = compute_group_sobol(c, mop, groups)
		S1_by_size[i] = S1
	end
	print("\r" * " "^70 * "\r")  # clear progress line

	return sample_sizes, S1_by_size
end

# ============================================================================
# Main
# ============================================================================

function main()
	println("\n" * "="^70)
	println("  Group Sobol Indices via Polynomial Chaos Expansion")
	println("="^70 * "\n")

	# --- Load data ---
	param_csv = readdlm(joinpath(RESULTS_DIR, "parameters.csv"), ','; header = true)
	param_matrix = Float64.(param_csv[1][:, 2:end])  # [n_samples × 17]
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
			println(@sprintf("    Z_%d (%-22s): μ=%+.4f, σ=%.4f", idx, g.name, μ_z, σ_z))
		end
	end

	# Verify independence between groups
	println("\n  Max cross-group |correlation|:")
	corr_Z = cor(Z)
	max_cross = 0.0
	for i in 1:N_PARAMS, j in (i+1):N_PARAMS
		gi = findfirst(g -> i ∈ g.indices, GROUPS)
		gj = findfirst(g -> j ∈ g.indices, GROUPS)
		if gi != gj
			max_cross = max(max_cross, abs(corr_Z[i, j]))
		end
	end
	println(@sprintf("    %.4f (should be small, ~0.02–0.05 for n=%d)", max_cross, n_success))

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
	degree_results = Dict{Int, NamedTuple}()

	for deg in PCE_DEGREE_RANGE
		ops_d = [GaussOrthoPoly(deg) for _ in 1:N_PARAMS]
		mop_d = MultiOrthoPoly(ops_d, deg)
		n_basis_d = mop_d.dim
		oversampling_d = n_success / n_basis_d

		if oversampling_d < 1.5
			println(@sprintf("  Degree %d: %d basis functions, oversampling %.1f× — SKIPPED (insufficient samples)",
				deg, n_basis_d, oversampling_d))
			continue
		end

		Φ_d = build_design_matrix(Z, mop_d)
		use_sparse_d = oversampling_d < SPARSE_THRESHOLD

		q2_values = Float64[]
		for qname in valid_qois
			qi = findfirst(==(qname), ALL_QOI_NAMES)
			y = qoi_data[:, qi]
			if use_sparse_d
				_, _, Q2_d, _ = fit_pce_lasso(Φ_d, y)
			else
				_, _, Q2_d, _ = fit_pce(Φ_d, y)
			end
			push!(q2_values, Q2_d)
		end
		avg_Q2 = mean(q2_values)
		min_Q2 = minimum(q2_values)

		fit_method = use_sparse_d ? "LASSO" : "OLS"
		println(@sprintf("  Degree %d: %d basis, %.1f× oversampling, avg Q²=%.4f, min Q²=%.4f (%s)",
			deg, n_basis_d, oversampling_d, avg_Q2, min_Q2, fit_method))

		degree_results[deg] = (mop = mop_d, Φ = Φ_d, n_basis = n_basis_d,
			oversampling = oversampling_d, avg_Q2 = avg_Q2, use_sparse = use_sparse_d)

		if avg_Q2 > best_avg_Q2
			best_avg_Q2 = avg_Q2
			best_degree = deg
		end
	end

	println(@sprintf("\n  → Selected degree %d (avg Q² = %.4f)", best_degree, best_avg_Q2))

	dr = degree_results[best_degree]
	mop = dr.mop
	Φ = dr.Φ
	n_basis = dr.n_basis
	use_sparse = dr.use_sparse

	# --- Step 2b: Verify basis normalization ---
	println("\n  Verifying basis orthonormality (E[Ψ²] ≈ 1)...")
	max_dev, norm_pass = verify_basis_normalization(Φ, mop)
	if norm_pass
		println(@sprintf("    ✓ PASSED (max |E[Ψ²]-1| = %.4f)", max_dev))
	else
		println(@sprintf("    ⚠ FAILED (max |E[Ψ²]-1| = %.4f) — Sobol indices may be biased", max_dev))
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
	all_S_first = Dict{String, Vector{Float64}}()
	all_S_total = Dict{String, Vector{Float64}}()
	all_R2 = Dict{String, Float64}()
	all_Q2 = Dict{String, Float64}()
	all_coeffs = Dict{String, Vector{Float64}}()
	all_predictions = Dict{String, Vector{Float64}}()
	all_boot_ci = Dict{String, NamedTuple}()

	for qname in QOI_NAMES
		qi = findfirst(==(qname), ALL_QOI_NAMES)
		y = qoi_data[:, qi]

		if !all(isfinite, y) || std(y) ≈ 0
			println("\n  $qname: SKIPPED (constant or non-finite)")
			continue
		end

		# Fit PCE (OLS or sparse LASSO depending on oversampling)
		if use_sparse
			c, R2, Q2, y_pred = fit_pce_lasso(Φ, y; verbose = true)
		else
			c, R2, Q2, y_pred = fit_pce(Φ, y)
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
		println(@sprintf("\n  %s  (R²=%.4f, Q²_LOO=%.4f, D=%.4e)",
			qname, R2, Q2, D))

		sorted = sortperm(S_first, rev = true)
		for (rank, gi) in enumerate(sorted)
			bar_first = repeat("█", round(Int, S_first[gi] * 40))
			bar_interact = repeat("░", round(Int, max(0, (S_total[gi] - S_first[gi])) * 40))
			println(@sprintf("    %d. %-22s S=%.4f  S_T=%.4f  %s%s",
				rank, group_names[gi], S_first[gi], S_total[gi], bar_first, bar_interact))
		end
		println(@sprintf("    Σ S_first = %.4f  (remainder = %.4f = interactions)",
			sum(S_first), 1 - sum(S_first)))

		# Verify Sobol index closure properties
		sum_ST = sum(S_total)
		if sum_ST < 1.0 - 0.01
			@warn "$qname: Σ S_T = $(round(sum_ST, digits=4)) < 1 — unexpected for total indices"
		end

		if R2 < 0.95
			@warn "$qname: R² = $(round(R2, digits=4)) < 0.95 — PCE fit may be inadequate"
		end
		if Q2 < 0.90
			@warn "$qname: Q²_LOO = $(round(Q2, digits=4)) < 0.90 — PCE may be overfitting"
		end
	end

	# --- Step 4: Bootstrap confidence intervals ---
	println("\n" * "-"^70)
	println("  Bootstrap Confidence Intervals ($N_BOOTSTRAP replicates)")
	println("-"^70)

	for qname in valid_qois
		qi = findfirst(==(qname), ALL_QOI_NAMES)
		y = qoi_data[:, qi]

		print("  $qname... ")
		S1_med, S1_lo, S1_hi, ST_med, ST_lo, ST_hi = bootstrap_sobol_ci(
			Z, y, mop, GROUPS, N_BOOTSTRAP;
			ci_level = BOOTSTRAP_CI_LEVEL, use_sparse = use_sparse)

		all_boot_ci[qname] = (S1_med = S1_med, S1_lo = S1_lo, S1_hi = S1_hi,
			ST_med = ST_med, ST_lo = ST_lo, ST_hi = ST_hi)

		println("done")
		ci_pct = round(Int, BOOTSTRAP_CI_LEVEL * 100)
		for gi in 1:n_groups
			if all_S_first[qname][gi] > 0.01  # only show non-negligible
				println(@sprintf("    %-22s S₁=%.4f [%.4f, %.4f]  S_T=%.4f [%.4f, %.4f]  (%d%% CI)",
					group_names[gi],
					all_S_first[qname][gi], S1_lo[gi], S1_hi[gi],
					all_S_total[qname][gi], ST_lo[gi], ST_hi[gi], ci_pct))
			end
		end
	end

	# --- Step 5: Convergence diagnostics ---
	println("\n" * "-"^70)
	println("  Convergence Diagnostics")
	println("-"^70)

	conv_results = Dict{String, Tuple{Vector{Int}, Vector{Vector{Float64}}}}()

	for qname in valid_qois
		qi = findfirst(==(qname), ALL_QOI_NAMES)
		y = qoi_data[:, qi]

		sample_sizes, S1_by_size = convergence_analysis(
			Z, y, GROUPS, best_degree, CONVERGENCE_FRACTIONS;
			use_sparse = use_sparse)
		conv_results[qname] = (sample_sizes, S1_by_size)

		# Check if indices have stabilized
		if length(S1_by_size) >= 2
			final = S1_by_size[end]
			prev = S1_by_size[end-1]
			max_change = maximum(abs.(final .- prev))
			status = max_change < 0.02 ? "✓ converged" : "⚠ not converged"
			println(@sprintf("  %s: max |ΔS₁| = %.4f (last two sizes) — %s",
				qname, max_change, status))
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
			parts = [gname]
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
		println(io, "qoi,group,S1,S1_lo_$(ci_pct),S1_hi_$(ci_pct),ST,ST_lo_$(ci_pct),ST_hi_$(ci_pct)")
		for qname in sort(collect(keys(all_boot_ci)))
			ci = all_boot_ci[qname]
			for (gi, gname) in enumerate(group_names)
				println(io, @sprintf("%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f",
					qname, gname,
					all_S_first[qname][gi], ci.S1_lo[gi], ci.S1_hi[gi],
					all_S_total[qname][gi], ci.ST_lo[gi], ci.ST_hi[gi]))
			end
		end
	end
	println("  Saved $boot_file")

	# --- PCE validation CSV ---
	validation_file = joinpath(RESULTS_DIR, "pce_validation.csv")
	open(validation_file, "w") do io
		println(io, "qoi,degree,fit_method,R2,Q2_LOO,n_basis,n_samples,oversampling_ratio")
		fit_method = use_sparse ? "LASSO" : "OLS"
		for qname in sort(collect(keys(all_R2)))
			println(io, @sprintf("%s,%d,%s,%.6f,%.6f,%d,%d,%.1f",
				qname, best_degree, fit_method, all_R2[qname], all_Q2[qname],
				n_basis, n_success, n_success / n_basis))
		end
	end
	println("  Saved $validation_file")

	# --- Per-sample PCE predictions CSV (for validation scatter plots) ---
	pred_file = joinpath(RESULTS_DIR, "pce_predictions.csv")
	open(pred_file, "w") do io
		pred_qoi_keys = sort(collect(keys(all_predictions)))
		header = join(["sample_id"; ["$(q)_true,$(q)_pred" for q in pred_qoi_keys]], ",")
		println(io, header)
		for i in 1:n_success
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
		idx_header = join(["alpha_$i" for i in 1:N_PARAMS], ",")
		qoi_header = join(qoi_keys, ",")
		println(io, "basis_idx,$idx_header,$qoi_header")
		for j in 1:n_basis
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
					println(io, @sprintf("%s,%d,%s,%.6f",
						qname, ns, gname, S1_by_size[si][gi]))
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
	colors = [:steelblue, :coral, :seagreen, :mediumpurple, :goldenrod, :lightblue, :pink, :sienna]

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
		bar!(1:n_qois, cumulative,
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

		plt = bar(1:n_groups, S1[sorted],
			label = "First-order (S₁)",
			color = :steelblue, alpha = 0.85,
			yticks = (1:n_groups, group_names[sorted]),
			orientation = :horizontal,
			xlims = (0, 1.05),
			xlabel = "Sobol Index",
			title = "Variance Decomposition — $(replace(qname, "_" => " "))",
			size = (800, 450),
			left_margin = 18Plots.mm,
		)
		bar!(1:n_groups, ST[sorted],
			label = "Total (S_T)",
			color = :coral, alpha = 0.4,
			orientation = :horizontal,
		)

		# Add bootstrap CI whiskers if available
		if haskey(all_boot_ci, qname)
			ci = all_boot_ci[qname]
			for (rank, gi) in enumerate(sorted)
				# S1 CI
				plot!([ci.S1_lo[gi], ci.S1_hi[gi]], [rank, rank],
					color = :black, lw = 1.5, label = (rank == 1 ? "95% CI" : false))
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

		plt = scatter(y, y_pred,
			xlabel = "Simulation", ylabel = "PCE prediction",
			title = @sprintf("%s (R²=%.4f, Q²=%.4f)", replace(qname, "_" => " "), R2, Q2),
			legend = false, alpha = 0.3, ms = 2, color = :steelblue,
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
		for i in 1:n_sizes
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
			plot!(sample_sizes, conv_matrix[:, gi],
				label = group_names[gi],
				lw = 2, marker = :circle, ms = 4,
				color = colors[mod1(gi, length(colors))])
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

	println(@sprintf("\n  PCE degree: %d (%s)", best_degree, use_sparse ? "LASSO" : "OLS"))

	println("\n  PCE fit quality:")
	for qname in qoi_keys
		status = all_R2[qname] > 0.95 && all_Q2[qname] > 0.90 ? "✓" : "⚠"
		println(@sprintf("    %s %-25s R²=%.4f  Q²_LOO=%.4f",
			status, qname, all_R2[qname], all_Q2[qname]))
	end

	ci_pct = round(Int, BOOTSTRAP_CI_LEVEL * 100)
	println("\n  Dominant measurement sources (S₁ > 0.05, with $ci_pct% bootstrap CI):")
	for qname in qoi_keys
		S1 = all_S_first[qname]
		important = [(gi, group_names[gi], S1[gi]) for gi in 1:n_groups if S1[gi] > 0.05]
		sort!(important, by = x -> -x[3])
		if haskey(all_boot_ci, qname)
			ci = all_boot_ci[qname]
			parts = [@sprintf("%s %.0f%% [%.0f, %.0f]", name, s * 100,
				ci.S1_lo[gi] * 100, ci.S1_hi[gi] * 100)
					 for (gi, name, s) in important]
		else
			parts = [@sprintf("%s (%.0f%%)", name, s * 100) for (_, name, s) in important]
		end
		println("    $(qname): $(join(parts, ", "))")
	end

	println("\n  Results in: $RESULTS_DIR")
	println("="^70)
end

main()
