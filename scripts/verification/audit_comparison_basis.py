"""Audit every retained cross-code pair; never infer a grid from probe count.

Run export_comparison_meshes.jl first. Native mesh metadata is produced by
matched/run_comparisons.py. Outputs are the review matrix and independent
time-refinement changes, not a performance ranking.
"""
from pathlib import Path
import bisect
import csv
import hashlib
import json
import math
import tomllib

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "test/verification/reference"
OUT = ROOT / "docs/guide/verification/assets"


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def series(path):
    out = {}
    for r in csv.DictReader(path.open()):
        out.setdefault(int(r["series"]), []).append((float(r["x"]), float(r["y"])))
    for xy in out.values():
        if not all(math.isfinite(v) for row in xy for v in row):
            raise ValueError(f"nonfinite samples: {path}")
        if any(b[0] <= a[0] for a, b in zip(xy, xy[1:])):
            raise ValueError(f"nonmonotone samples: {path}")
    return out


def change(coarse, fine):
    """Compare on common support; disclose truncation, never extrapolate.

    Linear interpolation is an observation operator, not a timestep estimate.
    A later qualification must bound its error using denser output.
    """
    xs = [x for x, _ in fine]
    differences = []
    for x, y in coarse:
        if x < xs[0] or x > xs[-1]:
            continue
        j = bisect.bisect_left(xs, x)
        if j == 0 or xs[j] == x:
            yf = fine[j][1]
        else:
            a, b = fine[j-1], fine[j]
            yf = a[1] + (b[1]-a[1]) * (x-a[0])/(b[0]-a[0])
        differences.append(abs(y-yf))
    if not differences:
        raise ValueError("no common observation support")
    return max(differences), len(differences), len(coarse)-len(differences)


def load_meshes():
    document = tomllib.loads((REF / "comparison_meshes.toml").read_text())
    for path, digest in document["source_hashes"].items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest:
            raise ValueError(f"stale Julia mesh export: {path}; rerun export_comparison_meshes.jl")
    data = document["meshes"]
    return {r["case"]: r for r in data if r["plotted"]}, data


# Time-only configurations at the base mesh, by their step relative to dt.
TIME_LEVELS = (("time_quadruple", "4dt"), ("time_double", "2dt"),
               ("time_half", "dt/2"), ("time_quarter", "dt/4"))


def make_row(case, code, source, native, julia, dt, time_status, note=""):
    j = julia[case]
    n, h0, h1 = native
    # Compare spacing and extent, not nodes vs cells. Scalar limits suffice
    # only for uniform grids; nonuniform meshes need full face coordinates.
    uniform = math.isclose(h0, h1, rel_tol=1e-5, abs_tol=1e-12)
    match = (uniform and j["cells"] == n and
             math.isclose(j["spacing_min_m"], h0, rel_tol=1e-5) and
             math.isclose(j["spacing_max_m"], h1, rel_tol=1e-5))
    homogeneous = case.startswith(("V-RX-01", "V-RX-02", "V-RX-05", "V-CPL-03"))
    status = "homogeneous ODE; count is not work equivalence" if homogeneous else (
        "spacing matched" if match else "mesh mismatch")
    if case in ("V-HC-04c", "V-HC-04d"):
        status = "material interface represented differently by each code"
    return dict(case=case, code=code, source=source, julia_plotted_cells=j["cells"],
                julia_h_min_m=j["spacing_min_m"], julia_h_max_m=j["spacing_max_m"],
                native_intervals=n, native_h_min_m=h0, native_h_max_m=h1,
                requested_dt_s=dt, spatial_status=status, temporal_status=time_status,
                notes=note)


def main():
    julia, _ = load_meshes()
    rows, changes = [], []
    for solver in ("fds", "gpyro", "thermakin"):
        root = REF / "matched"
        manifest = root / f"{solver}.toml"
        bindings = tomllib.loads(manifest.read_text())["series"]
        provenance = json.loads((root / f"{solver}_provenance.json").read_text())
        if hashlib.sha256(manifest.read_bytes()).hexdigest() != provenance["manifest_sha256"]:
            raise ValueError(f"stale manifest: {manifest}")
        runs = {(r["case"], r["level"]): r for r in provenance["runs"]}
        for run in runs.values():
            for filename, digest in run["files"].items():
                if hashlib.sha256((root / filename).read_bytes()).hexdigest() != digest:
                    raise ValueError(f"stale native evidence: {filename}")
        seen = set()
        for binding in bindings:
            case, path = binding["case"], binding["path"]
            native_case = Path(path).stem
            run = runs[native_case, "base"]
            grid = run["grid"]
            key = case, path
            if key not in seen:
                seen.add(key)
                rows.append(make_row(case, solver, "matched/" + path,
                    (grid["intervals"], grid["spacing_min_m"], grid["spacing_max_m"]),
                    julia, run["dt"], ", ".join(step for level, step in TIME_LEVELS
                                                if (native_case, level) in runs) +
                    " runs retained at fixed mesh",
                    grid["coordinate_convention"] + "; " + run["timestep_interpretation"] +
                    ("; " + grid["coordinate_exception"] if "coordinate_exception" in grid else "")))
            levels = [series(root / prefix / path)[binding["index"]]
                      for prefix in ("", "time_half", "time_quarter")]
            for a, b in (("base", "time_half"), ("time_half", "time_quarter")):
                if runs[native_case, a]["grid"] != runs[native_case, b]["grid"]:
                    raise ValueError(f"time refinement changed mesh: {native_case}")
            if runs[native_case, "time_quarter"]["dt"] != runs[native_case, "space_half"]["dt"]:
                raise ValueError(f"space refinement changed dt: {native_case}")
            d1, samples1, dropped1 = change(levels[0], levels[1])
            d2, samples2, dropped2 = change(levels[1], levels[2])
            changes.append(dict(case=case, code=solver, curve=binding["curve"],
                reference=binding["reference"], delta_dt_to_half=d1, delta_half_to_quarter=d2,
                samples_base=samples1, samples_half=samples2,
                omitted_base_endpoints=dropped1, omitted_half_endpoints=dropped2,
                observation="linear interpolation on common support; native units; no order inferred"))

    rows.sort(key=lambda r: (r["case"], r["code"], r["source"]))
    step_path = OUT / "comparison_julia_steps.csv"
    completed = {}
    if step_path.is_file():
        for r in csv.DictReader(step_path.open()):
            completed.setdefault((r["case"], int(r["cells"]), r["variant"]), set()).add(int(r["level"]))
    for r in rows:
        stem = Path(r["source"]).stem
        variant = "P" + stem[-1] if stem.startswith("vcpl03p") else ""
        ready = completed.get((r["case"], r["native_intervals"], variant)) == {0, 1, 2}
        r["julia_temporal_evidence"] = ("three tolerance/step-cap runs retained at native interval count"
                                        if ready else "pending")
    write_csv(OUT / "comparison_basis.csv", rows)
    write_csv(OUT / "comparison_time_changes.csv", changes)
    # A small TOML handoff lets Julia run the same interval-count workloads
    # without a CSV dependency. Geometry/observation warnings remain binding.
    targets = {}
    for r in rows:
        stem = Path(r["source"]).stem
        variant = "P" + stem[-1] if stem.startswith("vcpl03p") else ""
        targets[r["case"], r["native_intervals"], variant] = r
    with (REF / "comparison_targets.toml").open("w") as stream:
        for (case, n, variant), r in sorted(targets.items()):
            stream.write("[[targets]]\n")
            for key, value in dict(case=case, n_cells=n, variant=variant,
                                   spatial_status=r["spatial_status"]).items():
                stream.write(f"{key} = {json.dumps(value, ensure_ascii=False)}\n")
            stream.write("\n")
    print(f"Audited {len(rows)} case/code/subcase pairs across {len({r['case'] for r in rows})} case IDs")
    print(f"{sum(r['spatial_status'] == 'spacing matched' for r in rows)} have matched plotted spacing; "
          f"{sum(r['spatial_status'] == 'mesh mismatch' for r in rows)} have different plotted meshes")


if __name__ == "__main__":
    main()
