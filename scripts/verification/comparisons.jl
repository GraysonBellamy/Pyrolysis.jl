# Matched, native-input comparisons. A manifest identifies the exact case,
# curve, and analytical reference; it cannot silently borrow another variant.
using TOML

const MATCHED_COMPARISONS = vcat(
    [
        TOML.parsefile(joinpath(REFERENCE_DIR, "matched", "$(solver).toml"))["series"]
        for solver in ("fds", "gpyro")
    ]...,
)

function add_matched_comparisons!(case_id, curves, accuracy_path)
    entries = filter(s -> s["case"] == case_id, MATCHED_COMPARISONS)
    isempty(entries) && return curves
    audit = Any[]
    for entry in entries
        ci = findfirst(c -> string(c.name) == entry["curve"], curves)
        ci === nothing && error("$(case_id): missing comparison curve $(entry["curve"])")
        curve = curves[ci]
        si = findfirst(s -> !s.crosscode && s.label == entry["reference"], curve.series)
        si === nothing &&
            error("$(case_id): missing analytical reference $(entry["reference"])")
        exact = curve.series[si].exact
        exact === nothing &&
            error("$(case_id): comparison requires an analytical reference")
        for level in ("base", "refined", "time_half", "time_quarter", "space_half")
            prefix = level == "base" ? "matched/" : "matched/$(level)/"
            indices, x, y = reference_columns(
                prefix * entry["path"],
                "series",
                "x",
                "y";
                header_rows = 1,
            )
            keep = findall(==(entry["index"]), indices)
            isempty(keep) && error("$(case_id): missing retained comparison samples")
            xx, yy = x[keep], y[keep]
            all(isfinite, xx) && all(isfinite, yy) || error("nonfinite comparison data")
            errors = abs.(yy .- exact.(xx))
            all(isfinite, errors) || error("$(case_id): nonfinite analytical reference")
            push!(
                audit,
                (
                    entry["label"],
                    curve.name,
                    entry["reference"],
                    level,
                    length(keep),
                    maximum(errors),
                    sqrt(sum(abs2, errors) / length(errors)),
                ),
            )
            if level == "base"
                push!(
                    curve.series,
                    CurveSeries(
                        label = entry["label"],
                        x = xx,
                        y = yy,
                        exact = exact,
                        crosscode = true,
                    ),
                )
            end
        end
    end
    write_csv(
        accuracy_path,
        ["code", "curve", "reference", "resolution", "samples", "Linf", "RMS"],
        audit,
    )
    return curves
end

has_matched_comparison(case_id) = any(s -> s["case"] == case_id, MATCHED_COMPARISONS)
