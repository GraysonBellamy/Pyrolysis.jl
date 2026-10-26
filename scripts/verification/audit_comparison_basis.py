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
import re
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


def make_row(case, code, source, native, julia, all_julia, dt, time_status, note=""):
    j = julia[case]
    n, h0, h1 = native
    # Compare spacing and extent, not nodes vs cells. Scalar limits suffice
    # only for uniform grids; nonuniform meshes need full face coordinates.
    uniform = math.isclose(h0, h1, rel_tol=1e-5, abs_tol=1e-12)
    match = (uniform and j["cells"] == n and
             math.isclose(j["spacing_min_m"], h0, rel_tol=1e-5) and
             math.isclose(j["spacing_max_m"], h1, rel_tol=1e-5))
    homogeneous = case.startswith(("V-RX-01", "V-RX-02", "V-RX-05", "V-CPL-03"))
    rung = any(r["case"] == case and r["cells"] == n and
               math.isclose(r["spacing_min_m"], h0, rel_tol=1e-5) and
               math.isclose(r["spacing_max_m"], h1, rel_tol=1e-5) for r in all_julia)
    status = "homogeneous ODE; count is not work equivalence" if homogeneous else (
        "spacing matched" if match else "mesh mismatch")
    action = "retain base mesh; qualify time and observations" if match or homogeneous else (
        f"rerender Julia at existing {n}-cell rung; qualify time/observations" if rung else
        f"rerun Julia with {n} cells; verify the complete native mesh")
    if case in ("V-HC-04c", "V-HC-04d"):
        status = "interface representation must be audited separately"
        action = "use native face coordinates and audit assigned-material interface; keep off-face test separate"
    return dict(case=case, code=code, source=source, julia_plotted_cells=j["cells"],
                julia_h_min_m=j["spacing_min_m"], julia_h_max_m=j["spacing_max_m"],
                native_intervals=n, native_h_min_m=h0, native_h_max_m=h1,
                requested_dt_s=dt, spatial_status=status, temporal_status=time_status,
                runtime_qualified=False, next_action=action, notes=note)


def migration_observables(path):
    """Use native sample coordinates; group profiles by time and probes by depth."""
    table = list(csv.DictReader(path.open()))
    keys = list(table[0])
    curves = {}
    for row in table:
        if "x_m" in keys:
            history = path.name.startswith("vcpl01")
            group, coordinate = ("x_m", "t_s") if history else ("t_s", "x_m")
            curves.setdefault(f"{group}={row[group]}:T_K", []).append(
                (float(row[coordinate]), float(row["T_K"])))
        else:
            for key in keys[1:]:
                curves.setdefault(key, []).append((float(row["t_s"]), float(row[key])))
    return curves


def audit_migration(julia, all_julia, rows, changes):
    root = REF / "gpyro/v2_migration"
    provenance = json.loads((root / "provenance.json").read_text())
    if provenance["revision"] != "bff0126c68f2a35c084ef04490ae327185f281a3":
        raise ValueError("unexpected Gpyro V2 migration revision")
    for filename, expected in provenance["archived_inputs_sha256"].items():
        if hashlib.sha256((root.parent / filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f"stale migration input: {filename}")
    for path, key in ((root/"run_migration.py", "generator_sha256"),
                      (root.parent/"extract_results.py", "extractor_sha256")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != provenance[key]:
            raise ValueError(f"stale migration script: {path}")
    for other in (REF/"matched/gpyro_provenance.json", REF/"gpyro/v2_accuracy/provenance.json"):
        if json.loads(other.read_text())["executable_sha256"] != provenance["executable_sha256"]:
            raise ValueError(f"V2 executable differs between comparison suites: {other}")
    hashes = dict(provenance["extracted_sha256"])
    runs = {(r["case"], r["level"]): r for r in provenance["runs"]}
    expected_runs = {(case, level) for case in ("vhc01", "vhc03", "vrad01", "vcpl01", "vrx01")
                     for level in ("base", "refined", "time_half", "time_quarter", "space_half")}
    if set(runs) != expected_runs or len(provenance["runs"]) != len(expected_runs):
        raise ValueError("incomplete or duplicate V2 migration study")
    for run in runs.values():
        if not run["successful_termination"]:
            raise ValueError(f"incomplete V2 run: {run['case']}")
        hashes.update(run["files"])
    for filename, expected in hashes.items():
        if hashlib.sha256((root / filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f"stale V2 migration evidence: {filename}")
    (OUT/"comparison_v2_migration_accuracy.csv").write_text((root/"accuracy.csv").read_text())
    for case in sorted({r["case"] for r in runs.values()}):
        base = runs[case, "base"]
        grid = base["grid"]
        for level, factor in (("time_half", 2), ("time_quarter", 4)):
            other = runs[case, level]
            if other["grid"] != grid or other["dt"] != base["dt"]/factor:
                raise ValueError(f"invalid migrated time refinement: {case}/{level}")
        finer = runs[case, "space_half"]
        if finer["dt"] != base["dt"]/4 or finer["grid"]["intervals"] != grid["intervals"]*(1 if case == "vrx01" else 2):
            raise ValueError(f"invalid migrated space refinement: {case}")
        rows.append(make_row(base["case_id"], "gpyro_v2",
            f"gpyro/v2_migration/raw/base/{case}/{case}.data",
            (grid["intervals"], grid["spacing_min_m"], grid["spacing_max_m"]),
            julia, all_julia, base["dt"], "independent dt, dt/2, dt/4 runs retained",
            grid["coordinate_convention"] + "; full native profiles; final-observation guard is one requested step; " +
            base["timestep_interpretation"]))
        for path in sorted(root.glob(case + "_*.csv")):
            curves = [migration_observables(root / prefix / path.name)
                      for prefix in ("", "time_half", "time_quarter")]
            if not curves[0].keys() == curves[1].keys() == curves[2].keys():
                raise ValueError(f"native profile observation times changed: {path}")
            for curve in curves[0]:
                d1, n1, omitted1 = change(curves[0][curve], curves[1][curve])
                d2, n2, omitted2 = change(curves[1][curve], curves[2][curve])
                changes.append(dict(case=base["case_id"], code="gpyro_v2",
                    curve=path.stem + ":" + curve, reference="native observable; time sensitivity only",
                    delta_dt_to_half=d1, delta_half_to_quarter=d2, samples_base=n1, samples_half=n2,
                    omitted_base_endpoints=omitted1, omitted_half_endpoints=omitted2,
                    observation="native coordinates/times; common support; no extrapolation"))


def main():
    julia, all_julia = load_meshes()
    rows, changes = [], []
    for solver in ("fds", "gpyro"):
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
                    julia, all_julia, run["dt"], "independent dt, dt/2, dt/4 runs retained",
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

    v2 = REF / "gpyro/v2_accuracy"
    v2_runs = json.loads((v2 / "provenance.json").read_text())["runs"]
    by_case = {}
    for run in v2_runs:
        raw = v2 / "raw" / run["level"] / run["case"]
        deck = raw / (run["case"] + ".data")
        if hashlib.sha256(deck.read_bytes()).hexdigest() != run["input_sha256"]:
            raise ValueError(f"stale V2 input: {deck}")
        for filename, digest in run["outputs_sha256"].items():
            if hashlib.sha256((raw / filename).read_bytes()).hexdigest() != digest:
                raise ValueError(f"stale V2 output: {filename}")
        text = deck.read_text()
        by_case.setdefault(run["case"], {})[run["level"]] = (
            int(re.search(r"NCELLZ\(1\)=(\d+)", text)[1]),
            float(re.search(r"DT0=([^,]+)", text)[1]),
        )
    for case, levels in by_case.items():
        n, dt = levels["runs"]
        if levels["time_half"] != (n, dt/2) or levels["time_quarter"] != (n, dt/4):
            raise ValueError(f"invalid V2 temporal refinement: {case}")
        if levels["space_half"][1] != dt/4:
            raise ValueError(f"V2 spatial refinement changed timestep: {case}")
        for base in sorted(v2.glob(case + "_*.csv")):
            tables = [list(csv.DictReader((v2 / prefix / base.name).open()))
                      for prefix in ("", "time_half", "time_quarter")]
            keys = list(tables[0][0])
            ncoord = 2 if base.stem.endswith("_profiles") else 1
            coordinates = [[tuple(float(r[k]) for k in keys[:ncoord]) for r in t] for t in tables]
            if not coordinates[0] == coordinates[1] == coordinates[2]:
                raise ValueError(f"V2 time-study observation grids differ: {base}")
            for column in keys[ncoord:]:
                values = [[float(r[column]) for r in t] for t in tables]
                if not all(math.isfinite(x) for vs in values for x in vs):
                    raise ValueError(f"nonfinite V2 time-study value: {base}")
                d1 = max(abs(a-b) for a, b in zip(values[0], values[1]))
                d2 = max(abs(a-b) for a, b in zip(values[1], values[2]))
                changes.append(dict(case=f"V-{case[1:3].upper()}-{case[3:]}", code="gpyro_v2",
                    curve=base.stem + ":" + column, reference="native observable; time sensitivity only",
                    delta_dt_to_half=d1, delta_half_to_quarter=d2,
                    samples_base=len(values[0]), samples_half=len(values[1]),
                    omitted_base_endpoints=0, omitted_half_endpoints=0,
                    observation="identical native coordinates/times; native units; no order inferred"))

    audit_migration(julia, all_julia, rows, changes)

    # Remaining older FDS overlays and the earlier V2 accuracy family.
    # FDS counts below are the archived timing-deck
    # counts; their sparse DEVC samples cannot independently verify a mesh.
    legacy = [
        ("V-HC-01", "fds", "fds/custom/vhc01_semi_infinite.fds", 133, .05, .1,
         "reported 133 cells; rerun to retain full grid and instantaneous probes"),
        ("V-RX-01", "fds", "fds/custom/two_step_solid_reaction_matched.fds", 1, .001, "native controller",
         "uniform kinetics; audit native outer step, solid substeps and averaging"),
    ]
    for suffix, native in zip("abcd", ("d_matched", "c_matched", "b_matched", "a_ext")):
        legacy.append((f"V-HC-02{suffix}", "fds", f"fds/custom/heat_conduction_{native}.fds",
                       101, .1, .5, "dimensionless equivalence: t_Julia=2*t_FDS, dt_Julia-equivalent=1 s; probes time averaged"))
    for suffix in "ab":
        legacy.extend([
            (f"V-GT-07{suffix}", "gpyro_v2", f"gpyro/v2_accuracy/decks/vgt07{suffix}.data", 160, 1., .05,
             "161 boundary nodes; spatial/time refinement coupled; point vs cell-average inventory references differ"),
            (f"V-GT-08{suffix}", "gpyro_v2", f"gpyro/v2_accuracy/decks/vgt08{suffix}.data", 80, .1, .005,
             "81 boundary nodes; spatial/time refinement coupled; steady reference does not verify startup transient"),
        ])
    for suffix in "abc":
        legacy.append((f"V-RX-05{suffix}", "gpyro_v2", f"gpyro/v2_accuracy/decks/vrx05{suffix}.data",
                       1, .001, .001, "two boundary nodes; existing refinement is time-only; add third dt; retained-gas constant mismatch in c"))
    for case, solver, source, n, length, dt, note in legacy:
        if not (REF / source).is_file():
            raise ValueError(f"missing legacy deck: {source}")
        time_status = ("independent dt, dt/2, dt/4 runs retained; common-observation audit pending"
                       if solver == "gpyro_v2" else "time-only three-level qualification required")
        rows.append(make_row(case, solver, source, (n, length/n, length/n), julia, all_julia,
                            dt, time_status, note))
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
        if ready and r["next_action"].startswith("rerun Julia with"):
            r["next_action"] = r["next_action"].replace("rerun Julia with", "render comparison curves with")
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
          f"{sum(r['spatial_status'] == 'mesh mismatch' for r in rows)} have different plotted meshes; "
          "all runtime comparisons require qualification")


if __name__ == "__main__":
    main()
