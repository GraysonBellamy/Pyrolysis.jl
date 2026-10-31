# The conversion onset's cost and effect: V-RX-06's deck at each order through
# burnout on the moving mesh, at conversion thresholds 1e-4, 1e-3 and 1e-2, with
# steps and evaluations, the conversion reached, and the largest relative
# error against the bare law at conversions at or above 100 x_th.
#
# Run from the repository root:
#   taskset -c 0-19 julia --project=test --startup-file=no --threads=1 \
#       test/gates/kinetics_threshold_ladder.jl

include(joinpath(@__DIR__, "kinetics_identity.jl"))

function ladder_main()
    m, _ = IdentityDecks.load_cases()
    build =
        (args...; kw...) ->
            Base.invokelatest(IdentityDecks.latest(m, :_vrx06_setup), args...; kw...)
    exact =
        (args...; kw...) ->
            Base.invokelatest(IdentityDecks.latest(m, :nth_order_xi), args...; kw...)
    solve_ =
        (args...; kw...) ->
            Base.invokelatest(IdentityDecks.latest(m, :suite_solve), args...; kw...)
    println(
        "| n | x_th | accepted | rejected | rhs | Jacobians | t_end [s] | last x | error, x ≥ 100 x_th |",
    )
    println("|---|---|---|---|---|---|---|---|---|")
    for (order, t_end) in ((0.0, 75.0), (0.5, 150.0), (2.0, 20000.0))
        for x_th in (1e-4, 1e-3, 1e-2)
            S = build(20, order; t_end, x_th)
            sol = solve_(S.problem; use_ale = true, abstol = 1e-10, reltol = 1e-8)
            x = [
                sum(view(sol.m, 1, :, j)) / sum(view(sol.m, 1, :, 1)) for
                j in eachindex(sol.t)
            ]
            error = maximum(
                abs(x[j] / exact(sol.t[j]; ξ0 = 1.0, k = S.k, n = order) - 1) for
                j in eachindex(sol.t) if
                exact(sol.t[j]; ξ0 = 1.0, k = S.k, n = order) >= 100x_th;
                init = 0.0,
            )
            s = sol.statistics
            println(
                "| $order | $x_th | $(s.accepted) | $(s.rejected) | $(s.rhs_evaluations) | " *
                "$(s.jacobians) | $(round(last(sol.t); sigdigits = 4)) | " *
                "$(round(last(x); sigdigits = 3)) | $(round(error; sigdigits = 3)) |",
            )
            flush(stdout)
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == @__FILE__() && exit(ladder_main())
