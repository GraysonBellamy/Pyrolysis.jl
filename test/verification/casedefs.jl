"""
Shared verification-case definitions (Test-free).

Single source consumed by both harnesses: assert mode (`runcases.jl`, run by
`Pkg.test()`) and report mode (`scripts/verification/run_reports.jl`); the
Technical Reference Guide §18.1.3 describes both modes. This file deliberately depends on nothing from `Test` so the report runner can
include it from its isolated plotting environment.

Coordinate convention (matches `test_analytical.jl`): z = 0 bottom/substrate,
z = L top/exposed surface; exact solutions use depth x = L − z.
"""

using Pyrolysis
using Pyrolysis.Solver: get_cell_positions
using Pyrolysis.Materials: DepletionLimiter, STEFAN_BOLTZMANN
using Pyrolysis.Discretization: BEER_LAMBERT
using Pyrolysis.Problem: step_profile
using Pyrolysis.Physics: ideal_gas_pressure, background_gas_concentration
using Pyrolysis.Properties: porosity, effective_gas_viscosity
using Statistics: mean
using SpecialFunctions: erfc, erfcx

if !isdefined(@__MODULE__, :VerificationCase)
    include(joinpath(@__DIR__, "exact", "exact_solutions.jl"))

    """
    One quantity of interest evaluated from a solved case: `error < tol` is the
    assertion. `comment` records tolerance provenance (the observed error when
    the tolerance was pinned).
    """
    Base.@kwdef struct QoIResult
        name::String
        value::Float64
        exact::Float64 = NaN
        error::Float64
        tol::Float64
        comment::String = ""
    end

    """
    Spatial convergence spot check: solve `build(n)` for each `n ∈ n_cells`,
    reduce each solution to a scalar error norm via `error_fn`, and assert the
    observed order between consecutive resolutions. `build === nothing` reuses
    the case's own builder. Δz ∝ 1/n is used for the order fit, so `n_cells`
    must refine the same mesh family.
    """
    Base.@kwdef struct ConvergenceSpec
        n_cells::Vector{Int}
        error_fn::Function
        expected_order::Float64
        min_order::Float64
        min_mean_order::Float64
        solve_kwargs::NamedTuple = (;)
        build::Union{Nothing,Function} = nothing
    end

    Base.@kwdef struct CIConfig
        n_cells::Int
        budget_s::Float64            # warm-case target (compilation excluded); asserted at 5× for CI variance
        convergence::Union{Nothing,ConvergenceSpec} = nothing
    end

    """
    One plottable data series for report mode: numeric samples plus an optional
    `exact(x)` callable on the same abscissa. When present, the renderer draws
    a densely-sampled reference line under the numeric markers and derives the
    pointwise-error panel `y .- exact.(x)`. Plain data — no plotting deps.

    `crosscode = true` marks a third-party solver's result (FDS, Gpyro, … from
    `reference/`): the renderer then emits each figure twice — a clean variant
    with our series only, and a `_vs_codes` variant where cross-code series
    are drawn muted (grays, per-code marker shapes) behind ours.
    """
    Base.@kwdef struct CurveSeries
        label::String
        x::Vector{Float64}
        y::Vector{Float64}
        exact::Union{Nothing,Function} = nothing
        crosscode::Bool = false
    end

    """
    One report-mode figure: a solver-vs-exact overlay built from `series`.
    `name` is the stable slug used in artifact filenames and matches the
    case's `qoi_curves` tag. Consumed only by `scripts/verification/`; inert
    in assert mode.
    """
    Base.@kwdef struct ReportCurve
        name::Symbol
        title::String
        xlabel::String
        ylabel::String
        series::Vector{CurveSeries}
        xscale::Symbol = :identity
        yscale::Symbol = :identity
    end

    """
    Single-source verification-case definition. Assert mode consumes `ci`;
    report mode consumes `report`. `build(n_cells) -> PyrolysisProblem`,
    `qois(sol) -> Vector{QoIResult}`.

    `solve_fn` overrides how a built problem is solved: `nothing` (default)
    uses the standard `solve(problem; progress=false, kwargs...)`; a function
    is called as `solve_fn(problem; kwargs...)` and must return an object with
    a `retcode` field (`:Success` on success). MMS cases use this to route
    through the wrapped-RHS `solve_mms` path (exact/mms.jl). The override
    applies to the convergence sweep as well.

    `expected_retcode` is the return code assert mode requires: `:Success`
    by default, `:Terminated` for a case whose answer is that a callback
    stops the solve.

    `report` keys understood by the report runner (all optional; `(;)` means
    the case is skipped in report mode):

    - `n_cells_ladder::Vector{Int}` — Richardson resolution ladder
      (fallback: single run at `ci.n_cells`).
    - `curves::Function` — `sol -> Vector{ReportCurve}`; the case-authored
      overlay/error figure data (fallback: tables only, no figures).
    - `error_norms::Function` — `sol -> NamedTuple` of named error norms,
      e.g. `(L2 = …, Linf = …)` (fallback: `ci.convergence.error_fn` as a
      single `(norm = …,)` entry).
    - `expected_order::Float64` — reference slope for the convergence plot
      (fallback: `ci.convergence.expected_order`).
    - `solve_kwargs::NamedTuple` — report-tightened solver settings
      (fallback: `case.solve_kwargs` merged with convergence `solve_kwargs`).
    - `temporal_sweep::NamedTuple` — `(n_cells, dts)` fixed-dt time-order
      sweep (V-MMS-01).
    - `param_sweep::NamedTuple` — `(name, values, build, measure)` with
      `build(n_cells, value) -> problem` and `measure(sol) -> NamedTuple`
      (V-CPL-02 activation-energy ladder). `measure(sol, value)` is also
      accepted, for when the reference depends on the swept value
      (V-CPL-05 heating-rate ladder).
    - informational: `profile_times`, `qoi_curves`, `heating`, `note`.
    """
    Base.@kwdef struct VerificationCase
        id::String
        title::String
        tags::Vector{Symbol}
        references::Vector{String}
        build::Function
        solve_kwargs::NamedTuple = (;)
        solve_fn::Union{Nothing,Function} = nothing
        expected_retcode::Symbol = :Success
        qois::Function
        ci::CIConfig
        report::NamedTuple = (;)
    end

    """
    The core the suite solves on: `PYROLYSIS_CORE=new` routes the cases of
    `NEW_CORE_CASES` (by id prefix) through the finite-volume core
    (`Pyrolysis.Integration.solve`), in the energy state named by
    `PYROLYSIS_ENERGY_STATE` (`temperature`, the default, or `enthalpy`);
    every other case is skipped there.
    """
    const VERIFICATION_CORE = get(ENV, "PYROLYSIS_CORE", "old")
    const VERIFICATION_ENERGY_STATE =
        Symbol(get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature"))
    const NEW_CORE_CASES = (
        "V-RX-01",
        "V-RX-02",
        "V-RX-03",
        "V-HC-01",
        "V-HC-02",
        "V-HC-03",
        "V-HC-04",
        "V-HC-05",
    )

    runs_on_core(case::VerificationCase) =
        VERIFICATION_CORE != "new" || any(p -> startswith(case.id, p), NEW_CORE_CASES)

    @inline _case_solve(case::VerificationCase, problem; kwargs...) =
        case.solve_fn !== nothing ? case.solve_fn(problem; kwargs...) :
        VERIFICATION_CORE == "new" ?
        Pyrolysis.Integration.solve(
            problem;
            energy_state = VERIFICATION_ENERGY_STATE,
            kwargs...,
        ) : solve(problem; progress = false, kwargs...)

    """
    Cell-center positions of a solved case [m], on either core.
    """
    cell_positions(sol::PyrolysisSolution) = get_cell_positions(sol)
    cell_positions(sol::Pyrolysis.Outputs.CoreSolution) =
        Pyrolysis.Outputs.cell_positions(sol)

    observed_orders(errors, hs) =
        [log(errors[i-1] / errors[i]) / log(hs[i-1] / hs[i]) for i = 2:length(errors)]

    """
    Third-party solver outputs committed for report-mode overlays (FDS, Gpyro,
    MaCFP submissions — provenance in `reference/*/README.md`). Loaded only
    inside report `curves` closures; assert mode never touches these files.
    """
    const REFERENCE_DIR = joinpath(@__DIR__, "reference")

    """
        read_reference_csv(relpath; header_rows = 2) -> (names, cols)

    Parse `reference/<relpath>`: `header_rows` header lines (name row + units
    row, in either order — some MaCFP institutions swap them, so the row with
    the most letters is taken as the names), then numeric data. Returns the
    trimmed column names and one `Vector{Float64}` per column. Short rows pad
    columns they lack (files with a trailing comma drop the phantom column).
    """
    function read_reference_csv(relpath::AbstractString; header_rows::Int = 2)
        lines = readlines(joinpath(REFERENCE_DIR, relpath))
        head, body = lines[1:header_rows], lines[(header_rows+1):end]
        # The names row is the one whose first token is time-like ("Time",
        # "t_s", …) — letter counting is ambiguous (a units row like
        # "s,kg/m3,kg/m3" can tie a names row like "Time,Ya,Yb").
        is_time(tok) = (
            t = lowercase(strip(c -> isspace(c) || c == '"', tok));
            startswith(t, "time") || t in ("t", "t_s")
        )
        scores = [
            any(is_time, split(h, ',')[1:1]) ? count(isletter, h) + 1000 :
            count(isletter, h) for h in head
        ]
        name_row = head[argmax(scores)]
        names = String[
            strip(c -> isspace(c) || c == '"', s) for
            s in split(name_row, ',') if !isempty(strip(c -> isspace(c) || c == '"', s))
        ]
        cols = [Float64[] for _ in names]
        for ln in body
            isempty(strip(ln)) && continue
            vals = split(ln, ',')
            for i in eachindex(names)
                i ≤ length(vals) &&
                    !isempty(strip(vals[i])) &&
                    push!(cols[i], parse(Float64, strip(vals[i])))
            end
        end
        return names, cols
    end

    """
        reference_columns(relpath, colnames...) -> NTuple{N,Vector{Float64}}

    Fetch named columns from a reference CSV (exact match on trimmed names).
    """
    function reference_columns(
        relpath::AbstractString,
        colnames::AbstractString...;
        header_rows::Int = 2,
    )
        names, cols = read_reference_csv(relpath; header_rows)
        return map(colnames) do cn
            i = findfirst(==(cn), names)
            i === nothing && error(
                "column \"$cn\" not in reference/$relpath (has: $(join(names, ", ")))",
            )
            cols[i]
        end
    end

    const VERIFICATION_CASE_BUILDERS = Function[]

    """
    Discover and build every case in `cases/`: empties the registry, includes
    each file in sorted order (each pushes builders onto
    `VERIFICATION_CASE_BUILDERS`), and returns the built `VerificationCase`s.
    """
    function load_verification_cases()
        empty!(VERIFICATION_CASE_BUILDERS)
        for f in sort(filter(endswith(".jl"), readdir(joinpath(@__DIR__, "cases"))))
            include(joinpath(@__DIR__, "cases", f))
        end
        # invokelatest: the builders were defined by `include` DURING this
        # function call, so they are newer than the function's world age.
        return [Base.invokelatest(builder) for builder in VERIFICATION_CASE_BUILDERS]
    end
end
