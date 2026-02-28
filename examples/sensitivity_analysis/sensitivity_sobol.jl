"""
Group Sobol Sensitivity Indices via Polynomial Chaos Expansion

Computes variance-based sensitivity indices for the Douglas Fir cone calorimeter
simulation, decomposing output variance into contributions from each measurement
source (instrument group).

Method:
  1. Transform correlated parameter samples to independent standard normals
	 using the Cholesky factors of each group's covariance matrix.
  2. Fit a degree-2 Polynomial Chaos Expansion (PCE) surrogate to each scalar
	 QoI using least-squares regression on a Hermite polynomial basis.
  3. Compute first-order and total group Sobol indices analytically from the
	 PCE coefficients.

The 17 uncertain parameters fall into 8 independent groups (correlations are
within-group only), so standard Sobol formulas apply without correlated-input
extensions.

Groups (by measurement source):
  1. HFM calibration bias (ε_HFM)           — 1 parameter
  2. DSC calibration bias (ε_STA)           — 1 parameter
  3. Water calibration bias (ε_water)       — 1 parameter
  4. Virgin Cp, low-T (HFM regression)      — 2 parameters (a₁, b₁)
  5. Virgin Cp, high-T (DSC regression)     — 2 parameters (a₃, b₃)
  6. Char Cp (DSC regression)               — 2 parameters (a_char, b_char)
  7. Heat of sorption (WLS fit)             — 3 parameters (a, b, c)
  8. Heats of reaction (conditional cov.)   — 5 parameters (ΔH₁–ΔH₅)

Reads:  results/parameters.csv, results/scalar_qois.csv
Writes: results/sobol_group_indices.csv, results/pce_validation.csv,
		results/pce_coefficients.csv, results/sobol_variance_decomposition.png,
		results/sobol_bar_<qoi>.png, results/pce_validation_<qoi>.png

Usage:  julia --project=. examples/sensitivity_analysis/sensitivity_sobol.jl
"""

using PolyChaos
using LinearAlgebra
using Statistics
using DelimitedFiles
using Printf
using Plots

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
	# 1. HFM calibration bias
	ParamGroup("HFM calibration", [1],
		[0.0],
		[0.005774^2;;]),

	# 2. DSC calibration bias
	ParamGroup("DSC calibration", [2],
		[0.0],
		[0.02613^2;;]),

	# 3. Water calibration bias (isothermal HoS calibration)
	ParamGroup("Water calibration", [3],
		[0.0],
		[0.01732^2;;]),

	# 4. Virgin Cp, low-T (HFM regression)
	ParamGroup("Virgin Cp, low-T", [4, 5],
		[-53.05, 4.245],
		[0.21819^2 * 1e6   -1.5350e-4 * 1e6;
			-1.5350e-4 * 1e6   7.038e-4^2 * 1e6]),

	# 5. Virgin Cp, high-T (DSC regression)
	ParamGroup("Virgin Cp, high-T", [6, 7],
		[-2131.0, 9.119],
		[0.056399^2 * 1e6    -7.1999e-6 * 1e6;
			-7.1999e-6 * 1e6    1.2769e-4^2 * 1e6]),

	# 6. Char Cp (independent STA)
	ParamGroup("Char Cp", [8, 9],
		[347.2, 1.784],
		[3.486e-3^2 * 1e6    -2.0047e-8 * 1e6;
			-2.0047e-8 * 1e6    5.782e-6^2 * 1e6]),

	# 7. Heat of sorption (WLS exponential fit)
	ParamGroup("Heat of sorption", [10, 11, 12],
		[2280.0, 787.0, 0.114],
		[173.957     22.025     0.09442;
			 22.025    941.547     0.18431;
			  0.09442    0.18431   9.715e-5]),

	# 8. Heats of reaction (conditional covariance, kJ²/kg² → J²/kg²)
	ParamGroup("Heats of reaction", [13, 14, 15, 16, 17],
		[-19153.3, -5573.8, 266699.7, -219638.7, -12381.6],
		[3.42252    8.04561   11.3715    50.1601    12.3280;
			8.04561   20.7989    28.7027   127.515     31.2934;
			11.3715    28.7027    41.1541   178.306     45.5874;
			50.1601   127.515    178.306    805.909    176.265;
			12.3280    31.2934    45.5874   176.265    313.704] .* 1e6),
]

const N_PARAMS = 17

# Scalar QoIs to analyze (skip time_to_ignition — quantized, and solve_time — not physical)
const QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr",
	"total_mass_loss_pct", "peak_surface_temp", "peak_back_temp",
	"final_thickness",
]

# Column indices in scalar_qois.csv (after sample_id and retcode columns)
const ALL_QOI_NAMES = [
	"peak_mlr", "time_to_peak_mlr", "time_to_ignition",
	"total_mass_loss_pct", "peak_surface_temp", "peak_back_temp",
	"final_thickness", "solve_time",
]

# PCE degree
const PCE_DEGREE = 2

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
"""
function build_design_matrix(Z::Matrix{Float64}, mop::MultiOrthoPoly)
	n_samples = size(Z, 1)
	n_basis = mop.dim

	# Evaluate monic basis
	Φ = Matrix{Float64}(undef, n_samples, n_basis)
	for j in 1:n_basis
		Φ[:, j] = PolyChaos.evaluate(mop.ind[j, :], Z, mop)
	end

	# Normalize to orthonormal: divide by √(∏ αᵢ!)
	for j in 1:n_basis
		α = mop.ind[j, :]
		norm_factor = sqrt(prod(factorial.(α)))
		if norm_factor > 0
			Φ[:, j] ./= norm_factor
		end
	end

	return Φ
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

	# --- Step 2: Build PCE basis ---
	println("\nBuilding PCE basis (degree $PCE_DEGREE, $N_PARAMS dimensions)...")
	ops = [GaussOrthoPoly(PCE_DEGREE) for _ in 1:N_PARAMS]
	mop = MultiOrthoPoly(ops, PCE_DEGREE)
	n_basis = mop.dim
	println("  Basis functions: $n_basis")
	println("  Oversampling ratio: $(round(n_success / n_basis, digits=1))×")

	println("  Evaluating design matrix...")
	Φ = build_design_matrix(Z, mop)
	println("  Design matrix: $(size(Φ))")

	# Check condition number
	κ = cond(Φ' * Φ)
	println(@sprintf("  Condition number (ΦᵀΦ): %.2e", κ))
	if κ > 1e10
		@warn "Design matrix is poorly conditioned (κ > 1e10). Results may be unreliable."
	end

	# --- Step 3: Fit PCE and compute Sobol indices for each QoI ---
	println("\n" * "-"^70)
	println("  PCE Fitting & Group Sobol Indices")
	println("-"^70)

	group_names = [g.name for g in GROUPS]
	n_groups = length(GROUPS)

	# Storage for results
	all_S_first = Dict{String, Vector{Float64}}()
	all_S_total = Dict{String, Vector{Float64}}()
	all_R2 = Dict{String, Float64}()
	all_Q2 = Dict{String, Float64}()
	all_coeffs = Dict{String, Vector{Float64}}()
	all_predictions = Dict{String, Vector{Float64}}()

	for qname in QOI_NAMES
		qi = findfirst(==(qname), ALL_QOI_NAMES)
		y = qoi_data[:, qi]

		if !all(isfinite, y) || std(y) ≈ 0
			println("\n  $qname: SKIPPED (constant or non-finite)")
			continue
		end

		# Fit PCE
		c, R2, Q2, y_pred = fit_pce(Φ, y)
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
			bar_interact = repeat("░", round(Int, (S_total[gi] - S_first[gi]) * 40))
			println(@sprintf("    %d. %-22s S=%.4f  S_T=%.4f  %s%s",
				rank, group_names[gi], S_first[gi], S_total[gi], bar_first, bar_interact))
		end
		println(@sprintf("    Σ S_first = %.4f  (remainder = %.4f = interactions)",
			sum(S_first), 1 - sum(S_first)))

		if R2 < 0.95
			@warn "$qname: R² = $(round(R2, digits=4)) < 0.95 — PCE fit may be inadequate"
		end
		if Q2 < 0.90
			@warn "$qname: Q²_LOO = $(round(Q2, digits=4)) < 0.90 — PCE may be overfitting"
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

	# --- PCE validation CSV ---
	validation_file = joinpath(RESULTS_DIR, "pce_validation.csv")
	open(validation_file, "w") do io
		println(io, "qoi,R2,Q2_LOO,n_basis,n_samples,oversampling_ratio")
		for qname in sort(collect(keys(all_R2)))
			println(io, @sprintf("%s,%.6f,%.6f,%d,%d,%.1f",
				qname, all_R2[qname], all_Q2[qname], n_basis, n_success,
				n_success / n_basis))
		end
	end
	println("  Saved $validation_file")

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
	colors = [:steelblue, :coral, :seagreen, :mediumpurple, :goldenrod, :lightblue, :pink]

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

	# --- 2. Individual bar charts per QoI ---
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

	# ====================================================================
	# Summary
	# ====================================================================

	println("\n" * "="^70)
	println("  Summary")
	println("="^70)

	println("\n  PCE fit quality:")
	for qname in qoi_keys
		status = all_R2[qname] > 0.95 && all_Q2[qname] > 0.90 ? "✓" : "⚠"
		println(@sprintf("    %s %-25s R²=%.4f  Q²_LOO=%.4f",
			status, qname, all_R2[qname], all_Q2[qname]))
	end

	println("\n  Dominant measurement sources (S₁ > 0.05):")
	for qname in qoi_keys
		S1 = all_S_first[qname]
		important = [(group_names[gi], S1[gi]) for gi in 1:n_groups if S1[gi] > 0.05]
		sort!(important, by = x -> -x[2])
		parts = [Printf.@sprintf("%s (%.0f%%)", name, s * 100) for (name, s) in important]
		println("    $(qname): $(join(parts, ", "))")
	end

	println("\n  Results in: $RESULTS_DIR")
	println("="^70)
end

main()
