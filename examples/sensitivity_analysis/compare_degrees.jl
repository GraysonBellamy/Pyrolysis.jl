#
# Definitive PCE degree comparison: degree 2 vs 3 vs 4
#
# Purpose
# -------
# Provide a rigorous, reproducible, unambiguous comparison of PCE degrees for
# the Douglas Fir Sobol analysis. The output is a verdict based on objective
# convergence criteria, not heuristics or opinion.
#
# Method
# ------
# For each degree d ∈ {2, 3, 4}:
#   1. Build the PCE design matrix Φ_d on the full n = 8192 Sobol samples.
#   2. Fit each QoI (OLS if oversampling ≥ SPARSE_THRESHOLD, LASSO otherwise).
#   3. Record the leave-one-out Q² (analytical for OLS, CV-based for LASSO).
#   4. Compute the 7 group first-order and total Sobol indices.
#
# Then for each QoI, compute the inter-degree step sizes:
#   Δ_{2→3}(g) = |S_deg3(g) − S_deg2(g)|
#   Δ_{3→4}(g) = |S_deg4(g) − S_deg3(g)|
#
# Decision rule (applied per QoI, then aggregated):
#
#   1. CONVERGED (prefer degree 3 — simpler fully-resolved model)
#        max_g Δ_{3→4}(g) < 0.01 (Sobol indices agree within 1%-point between d=3 and d=4)
#        AND Q²(deg3) > 0.95 (surrogate is high-quality at degree 3)
#
#   2. DEGREE 4 JUSTIFIED (higher degree captures real structure)
#        max_g Δ_{3→4}(g) < max_g Δ_{2→3}(g) (the step sizes are SHRINKING)
#        AND max_g Δ_{3→4}(g) ≥ 0.01 (still enough disagreement that degree 4 is meaningful)
#        AND degree 4's Q²_LOO is not dramatically worse than degree 3's
#
#   3. OVERFIT / DIVERGING (degree 4 is unreliable; prefer degree 3)
#        max_g Δ_{3→4}(g) ≥ max_g Δ_{2→3}(g) (step sizes NOT shrinking — not converging)
#        OR degree 4's Q²_LOO < degree 3's Q²_LOO (degree 4 fits WORSE)
#
# The verdict is determined by counting QoIs in each category and applying a
# clear aggregation rule. The script prints the evidence (per-QoI deltas),
# the per-QoI verdicts, and the overall verdict with justification.
#
# Usage:
#   julia --threads=auto --project=. examples/sensitivity_analysis/compare_degrees.jl
#

using DelimitedFiles
using Statistics
using Printf
using LinearAlgebra

function load_functions()
	src = read(joinpath(@__DIR__, "sensitivity_sobol.jl"), String)
	src = replace(src, r"^using Plots\s*$"m => "# Plots stripped for comparison script")
	src = replace(src, r"^main\(\)\s*$"m => "# main() disabled for comparison script")
	Base.include_string(@__MODULE__, src, "sensitivity_sobol_loaded_for_compare.jl")
	return nothing
end

load_functions()
using PolyChaos

const DEGREES_TO_TEST = [2, 3, 4]
const RESULTS_DIR = joinpath(@__DIR__, "results")

# Convergence thresholds (all in absolute Sobol index units, which live in [0,1])
const CONVERGENCE_TOL     = 0.01   # "converged" between two degrees
const Q2_ACCEPTABLE       = 0.95   # minimum Q² for "degree d is a usable surrogate"

function load_data()
	println("Loading data from $RESULTS_DIR ...")
	param_csv = readdlm(joinpath(RESULTS_DIR, "parameters.csv"), ','; header = true)
	scalar_csv = readdlm(joinpath(RESULTS_DIR, "scalar_qois.csv"), ','; header = true)

	param_matrix = Float64.(param_csv[1][:, 2:end])
	scalar_data = scalar_csv[1]
	success_mask = scalar_data[:, 2] .!= "FAILED"
	sample_ids = Int.(scalar_data[success_mask, 1])
	qoi_data = Float64.(scalar_data[success_mask, 3:end])
	param_data = param_matrix[sample_ids, :]

	n_success = size(param_data, 1)
	println("  Loaded $n_success successful samples")

	print("  Transforming to standard-normal space... ")
	flush(stdout)
	Z = Main.transform_to_standard_normal(param_data, Main.GROUPS)
	println("done")

	return Z, qoi_data, n_success
end

function collect_valid_qois(qoi_data)
	valid = String[]
	for qname in Main.QOI_NAMES
		qi = findfirst(==(qname), Main.ALL_QOI_NAMES)
		y = qoi_data[:, qi]
		if all(isfinite, y) && std(y) > 0
			push!(valid, qname)
		end
	end
	return valid
end

"""
	fit_all_qois_at_degree(Z, qoi_data, valid_qois, degree)

Fit a PCE at the specified degree for every valid QoI, returning a
NamedTuple with Q², first-order Sobol indices, total Sobol indices, and
the selected fit method ("OLS" or "LASSO") per QoI.
"""
function fit_all_qois_at_degree(Z, qoi_data, valid_qois, degree)
	println()
	println("─"^70)
	println("Degree $degree")
	println("─"^70)

	t0 = time()
	ops = [GaussOrthoPoly(degree) for _ in 1:Main.N_PARAMS]
	mop = Main.MultiOrthoPoly(ops, degree)
	Φ = Main.build_design_matrix(Z, mop)
	p_basis = mop.dim
	n_samples = size(Z, 1)
	oversampling = n_samples / p_basis
	use_sparse = oversampling < Main.SPARSE_THRESHOLD
	method = use_sparse ? "LASSO" : "OLS"

	@printf("  n_basis = %d, oversampling = %.2f, method = %s\n", p_basis, oversampling, method)

	# Share standardized Φ across QoIs if LASSO is used
	precomp = use_sparse ? Main.standardize_design_matrix(Φ) : nothing

	Q2_map = Dict{String, Float64}()
	S1_map = Dict{String, Vector{Float64}}()
	ST_map = Dict{String, Vector{Float64}}()
	nonzero_map = Dict{String, Int}()

	for qname in valid_qois
		qi = findfirst(==(qname), Main.ALL_QOI_NAMES)
		y = qoi_data[:, qi]

		if use_sparse
			c, R2, Q2, y_pred = Main.fit_pce_lasso(Φ, y;
				precomputed = precomp, verbose = false)
		else
			c, R2, Q2, y_pred, _, _ = Main.fit_pce(Φ, y)
		end

		Q2_map[qname] = Q2
		nonzero_map[qname] = count(abs.(c[2:end]) .> 1e-12)

		S1, ST, _ = Main.compute_group_sobol(c, mop, Main.GROUPS)
		S1_map[qname] = S1
		ST_map[qname] = ST
	end

	elapsed = round(time() - t0, digits = 1)
	mean_Q2 = mean(values(Q2_map))
	min_Q2 = minimum(values(Q2_map))
	@printf("  Fit time: %.1fs; mean Q² = %.5f; min Q² = %.5f\n",
		elapsed, mean_Q2, min_Q2)

	return (mop = mop, Φ = Φ, p_basis = p_basis, oversampling = oversampling,
		method = method, Q2 = Q2_map, S1 = S1_map, ST = ST_map,
		nonzero = nonzero_map)
end

"""
	max_delta(Sa, Sb)

Maximum absolute difference between two Sobol index vectors (of equal length).
"""
function max_delta(Sa::Vector{Float64}, Sb::Vector{Float64})
	return maximum(abs.(Sa .- Sb))
end

function classify_qoi(Δ23_S1, Δ23_ST, Δ34_S1, Δ34_ST, Q2_3, Q2_4)
	# Maximum step size across first-order AND total indices, across groups
	step_23 = max(Δ23_S1, Δ23_ST)
	step_34 = max(Δ34_S1, Δ34_ST)

	converged_34 = step_34 < CONVERGENCE_TOL
	shrinking = step_34 < step_23

	# Degree 4 "fits worse" means Q²(deg4) < Q²(deg3). This should rarely
	# happen for LOO Q² (higher-degree models usually fit better), but it's
	# a red flag if it does.
	q2_regression = Q2_4 < Q2_3 - 0.001

	verdict = if converged_34 && Q2_3 >= Q2_ACCEPTABLE
		"converged"
	elseif q2_regression
		"deg4_worse"
	elseif shrinking && step_34 >= CONVERGENCE_TOL
		"deg4_needed"
	elseif !shrinking
		"diverging"
	else
		"borderline"
	end

	return (step_23 = step_23, step_34 = step_34,
		verdict = verdict, q2_regression = q2_regression)
end

function main()
	println("="^70)
	println("  PCE Degree Comparison: 2 vs 3 vs 4")
	println("="^70)
	println()
	println("Criterion: Sobol index convergence as degree increases.")
	println("Decision:  Report the lowest degree at which indices stop moving,")
	println("           provided Q² is acceptable and convergence is monotone.")
	println()

	Z, qoi_data, n_samples = load_data()
	valid_qois = collect_valid_qois(qoi_data)
	println("Comparing $(length(valid_qois)) QoIs at $(length(DEGREES_TO_TEST)) degrees")

	# Fit at each degree
	results = Dict{Int, NamedTuple}()
	for d in DEGREES_TO_TEST
		results[d] = fit_all_qois_at_degree(Z, qoi_data, valid_qois, d)
	end

	# ---------- Per-QoI convergence table ----------
	println()
	println("="^70)
	println("  Per-QoI Sobol index convergence")
	println("="^70)
	println()
	println("For each QoI, the table below shows:")
	println("  • Q²: leave-one-out R² at each degree (higher = better fit)")
	println("  • step_2→3 = max |S(d=3) - S(d=2)|  across all groups and S₁, S_T")
	println("  • step_3→4 = max |S(d=4) - S(d=3)|  across all groups and S₁, S_T")
	println("  • classification (see legend below)")
	println()
	@printf("  %-28s %6s %6s %6s %9s %9s  %s\n",
		"QoI", "Q²_d2", "Q²_d3", "Q²_d4", "step_2→3", "step_3→4", "verdict")
	println("  " * "─"^92)

	verdict_counts = Dict("converged" => 0, "deg4_needed" => 0,
		"deg4_worse" => 0, "diverging" => 0, "borderline" => 0)
	per_qoi_verdicts = Dict{String, String}()
	max_overall_step_34 = 0.0
	worst_qoi_step_34 = ""

	for qname in sort(valid_qois)
		Q2_2 = results[2].Q2[qname]
		Q2_3 = results[3].Q2[qname]
		Q2_4 = results[4].Q2[qname]

		S1_2 = results[2].S1[qname]; ST_2 = results[2].ST[qname]
		S1_3 = results[3].S1[qname]; ST_3 = results[3].ST[qname]
		S1_4 = results[4].S1[qname]; ST_4 = results[4].ST[qname]

		Δ23_S1 = max_delta(S1_3, S1_2)
		Δ23_ST = max_delta(ST_3, ST_2)
		Δ34_S1 = max_delta(S1_4, S1_3)
		Δ34_ST = max_delta(ST_4, ST_3)

		c = classify_qoi(Δ23_S1, Δ23_ST, Δ34_S1, Δ34_ST, Q2_3, Q2_4)
		per_qoi_verdicts[qname] = c.verdict
		verdict_counts[c.verdict] += 1
		if c.step_34 > max_overall_step_34
			max_overall_step_34 = c.step_34
			worst_qoi_step_34 = qname
		end

		marker = if c.verdict == "converged"
			"✓ converged"
		elseif c.verdict == "deg4_needed"
			"→ deg4 needed"
		elseif c.verdict == "deg4_worse"
			"✗ deg4 WORSE Q²"
		elseif c.verdict == "diverging"
			"✗ diverging"
		else
			"~ borderline"
		end

		@printf("  %-28s %6.4f %6.4f %6.4f %9.5f %9.5f  %s\n",
			qname, Q2_2, Q2_3, Q2_4, c.step_23, c.step_34, marker)
	end

	println()
	println("  Legend:")
	println("    ✓ converged     step_3→4 < $CONVERGENCE_TOL AND Q²(d=3) ≥ $Q2_ACCEPTABLE")
	println("                    (degree 3 is fully resolved — use degree 3)")
	println("    → deg4 needed   step_3→4 ≥ $CONVERGENCE_TOL AND step_3→4 < step_2→3")
	println("                    (still converging, need degree 4)")
	println("    ✗ deg4 WORSE Q² degree 4 Q² is worse than degree 3 — overfit/bad fit")
	println("    ✗ diverging     step_3→4 ≥ step_2→3 — not converging (degree 4 unreliable)")
	println("    ~ borderline    none of the above fit cleanly")

	# ---------- Aggregate verdict ----------
	println()
	println("="^70)
	println("  Aggregate verdict")
	println("="^70)
	println()
	println("QoI verdict counts:")
	for k in ["converged", "deg4_needed", "deg4_worse", "diverging", "borderline"]
		n = verdict_counts[k]
		if n > 0
			@printf("  %-15s %d / %d\n", k * ":", n, length(valid_qois))
		end
	end
	@printf("\n  Worst 3→4 step: %.5f at QoI '%s'\n",
		max_overall_step_34, worst_qoi_step_34)

	# Overall Q² summary
	mean_Q2 = Dict(d => mean(values(results[d].Q2)) for d in DEGREES_TO_TEST)
	min_Q2  = Dict(d => minimum(values(results[d].Q2)) for d in DEGREES_TO_TEST)
	println("\n  Aggregate Q²:")
	for d in DEGREES_TO_TEST
		@printf("    degree %d (%-5s, p=%d, over=%.2f): mean Q² = %.5f, min Q² = %.5f\n",
			d, results[d].method, results[d].p_basis, results[d].oversampling,
			mean_Q2[d], min_Q2[d])
	end

	# ---------- Final decision ----------
	n_conv = verdict_counts["converged"]
	n_need4 = verdict_counts["deg4_needed"]
	n_bad4 = verdict_counts["deg4_worse"] + verdict_counts["diverging"]
	n_border = verdict_counts["borderline"]
	n_qois = length(valid_qois)

	println()
	println("─"^70)
	println("  Final decision")
	println("─"^70)

	chosen_degree = 0
	reason = ""

	if n_bad4 > 0 && n_conv + n_bad4 >= n_qois - n_border
		chosen_degree = 3
		reason = "degree 4 is unreliable ($(n_bad4) QoIs show Q² regression or " *
		         "diverging Sobol indices), while degree 3 is converged on the rest"
	elseif n_need4 == 0
		chosen_degree = 3
		reason = "all $(n_conv) non-borderline QoIs show converged Sobol indices " *
		         "between degree 3 and degree 4 (max 3→4 step = " *
		         @sprintf("%.5f", max_overall_step_34) * " < $CONVERGENCE_TOL)"
	elseif n_need4 > n_conv
		chosen_degree = 4
		reason = "$(n_need4) QoIs need degree 4 (Sobol indices still converging " *
		         "from degree 3 to degree 4); only $(n_conv) are converged at degree 3"
	elseif n_need4 > 0 && n_need4 <= n_conv
		chosen_degree = 4
		reason = "$(n_need4) QoIs need degree 4 and the remaining $(n_conv) QoIs " *
		         "are also converged at degree 4 (degree 4 is consistent with both)"
	else
		chosen_degree = 3
		reason = "no clear evidence that degree 4 is needed; prefer simpler model"
	end

	println()
	println("  → USE DEGREE $chosen_degree")
	println()
	println("  Reason: $reason")

	if chosen_degree == 3
		println()
		println("  To commit this decision, set in sensitivity_sobol.jl line ~160:")
		println("    const PCE_DEGREE_RANGE = 3:3")
		println()
		println("  Then re-run the pipeline. Degree 3 uses OLS (oversampling 8.5),")
		println("  the bootstrap will be fast (Fix B Cholesky), and no post-hoc")
		println("  justification for degree 4 is required.")
	else
		println()
		println("  To commit this decision, leave PCE_DEGREE_RANGE as the default")
		println("  and continue the pipeline. Degree 4's higher complexity is")
		println("  justified by the Sobol convergence evidence above.")
	end

	# Diagnostic outputs even if main verdict is clear
	if n_border > 0 || verdict_counts["deg4_worse"] > 0 || verdict_counts["diverging"] > 0
		println()
		println("  Notable exceptions:")
		for (qname, v) in per_qoi_verdicts
			if v in ["borderline", "deg4_worse", "diverging"]
				Q2_3 = results[3].Q2[qname]
				Q2_4 = results[4].Q2[qname]
				S1_3 = results[3].S1[qname]
				S1_4 = results[4].S1[qname]
				step = max_delta(S1_4, S1_3)
				@printf("    %-28s (%s): Q²_d3=%.4f, Q²_d4=%.4f, step_3→4=%.5f\n",
					qname, v, Q2_3, Q2_4, step)
			end
		end
	end

	# Save the comparison table as CSV for the paper supplement
	out_file = joinpath(RESULTS_DIR, "degree_comparison.csv")
	open(out_file, "w") do io
		header_parts = ["qoi", "Q2_d2", "Q2_d3", "Q2_d4", "step_2to3", "step_3to4", "verdict"]
		# Add per-group first-order and total index columns for each degree
		for d in DEGREES_TO_TEST
			for g in Main.GROUPS
				push!(header_parts, @sprintf("S1_%s_d%d", replace(g.name, " " => "_"), d))
			end
			for g in Main.GROUPS
				push!(header_parts, @sprintf("ST_%s_d%d", replace(g.name, " " => "_"), d))
			end
		end
		println(io, join(header_parts, ","))
		for qname in sort(valid_qois)
			Q2_2 = results[2].Q2[qname]; Q2_3 = results[3].Q2[qname]; Q2_4 = results[4].Q2[qname]
			S1_2 = results[2].S1[qname]; S1_3 = results[3].S1[qname]; S1_4 = results[4].S1[qname]
			ST_2 = results[2].ST[qname]; ST_3 = results[3].ST[qname]; ST_4 = results[4].ST[qname]
			Δ23 = max(max_delta(S1_3, S1_2), max_delta(ST_3, ST_2))
			Δ34 = max(max_delta(S1_4, S1_3), max_delta(ST_4, ST_3))
			row = [qname,
				@sprintf("%.6f", Q2_2), @sprintf("%.6f", Q2_3), @sprintf("%.6f", Q2_4),
				@sprintf("%.6f", Δ23), @sprintf("%.6f", Δ34),
				per_qoi_verdicts[qname]]
			for d in DEGREES_TO_TEST
				S1 = results[d].S1[qname]
				ST = results[d].ST[qname]
				for v in S1
					push!(row, @sprintf("%.6f", v))
				end
				for v in ST
					push!(row, @sprintf("%.6f", v))
				end
			end
			println(io, join(row, ","))
		end
	end
	println()
	println("  Full per-QoI data written to: $out_file")
	println("  Include this in your paper's supplementary material.")
	println("="^70)
end

main()
