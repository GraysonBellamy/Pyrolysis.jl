#!/usr/bin/env python3
"""Reproduce the five former 0.8200 overlays with pinned, unmodified Gpyro V2.

All five configurations are accuracy runs; no runtime claims are made.
The archived extractor supplies the same independent analytical checks and
SI output schema as before. Native coordinates and actual output times are
retained, rather than shifted to the requested observation times.
"""
from pathlib import Path
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
REVISION = "bff0126c68f2a35c084ef04490ae327185f281a3"
SUCCESS = "End of simulation reached successfully. Shutting down."
LEVELS = {"base": (1, 1), "refined": (2, 2), "time_half": (1, 2),
          "time_quarter": (1, 4), "space_half": (2, 4)}
# Base intervals match the published Julia meshes. RX-01 is a 0D ODE.
CASES = {
    "vhc01": dict(id="V-HC-01", cells=500, length=.05, dt=.1, end=180., point=1., profile=30.),
    "vhc03": dict(id="V-HC-03", cells=320, length=.02, dt=.1, end=1200., point=2., profile=80.),
    "vrad01": dict(id="V-RAD-01", cells=400, length=.01, dt=.1, end=180., point=1., profile=30.),
    "vcpl01": dict(id="V-CPL-01", cells=800, length=.04, dt=.02, end=300., point=1., profile=1.),
    "vrx01": dict(id="V-RX-01", cells=1, length=.001, dt=.001, end=20., point=.1, profile=1000.),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_value(text, key, value):
    text, count = re.subn(r"(?m)^" + re.escape(key) + r"\s*=\s*[^,\n]+",
                          f"{key} = {value}", text)
    if count != 1:
        raise ValueError(f"expected one {key} assignment, found {count}")
    return text


def deck_text(case, level):
    p = CASES[case]
    space, time = LEVELS[level]
    dt = p["dt"] / time
    text = (ARCHIVE / f"{case}.data").read_text()
    text = text[text.index("&GPYRO_GENERAL"):]
    # IO2=2 in the archived single-gas decks is invalid in V2. There is no
    # oxidation in these cases, so the sole inert gas receives index 1.
    changes = {"DT0": f"{dt:.17g}", "IO2": 1,
               "TSTOP(1)": f"{p['end'] + dt:.17g}",
               "NCELLZ(1)": 1 if case == "vrx01" else space*p["cells"] + 1,
               "DTDUMP_POINT": f"{p['point'] - 1e-8:.17g}",
               "DTDUMP_PROFILE": f"{p['profile'] - 1e-8:.17g}",
               "TMPTOL": "1D-8", "HTOL": "1D-12", "YITOL": "1D-10",
               "DTMIN_KILL": "1D-12"}
    if case == "vrx01":
        # The chain has unit condensed yields and zero porosity. Its old
        # unused gas-species flag is rejected by V2's porosity validation.
        changes["SOLVE_GAS_YJ"] = ".FALSE."
    if case != "vrx01":
        # V2's dynamically sized profile writer permits all thermal nodes.
        changes.update({"PROFILE_ISKIP(1)": 1, "POINT_Z(1)": 0,
                        "POINT_Z(2)": p["length"]})
    for key, value in changes.items():
        text = set_value(text, key, value)
    text = text.replace("&GPYRO_GENERAL\n", "&GPYRO_GENERAL\nSOLVER_DEFORMATION_MODE = 1,\n")
    text = text.replace("&GPYRO_GEOM\n", "&GPYRO_GEOM\nHALF_CELLS_AT_BC(1) = .TRUE.,\n")
    return (f"! Gpyro V2 {REVISION}; {p['id']}; {level}.\n"
            "! Same physical parameters as the archived deck.\n"
            "! TSTOP includes one requested step solely to retain the final observation.\n" + text)


def numeric_table(path):
    with path.open() as stream:
        rows = list(csv.reader(stream))
    data = [[float(v) for v in row if v.strip()] for row in rows[1:] if row]
    if not data or not all(math.isfinite(x) for row in data for x in row):
        raise ValueError(f"empty/nonfinite native output: {path}")
    if any(b[0] <= a[0] for a, b in zip(data, data[1:])):
        raise ValueError(f"nonmonotone native time: {path}")
    return rows[0], data


def run_one(job, exe, work, extract_only):
    case, level = job
    p = CASES[case]
    space, time = LEVELS[level]
    path = work / level / "gpyro_run" / case
    path.mkdir(parents=True, exist_ok=True)
    deck = path / f"{case}.data"
    expected = deck_text(case, level)
    if extract_only:
        if deck.read_text() != expected:
            raise ValueError(f"deck differs from generator: {deck}")
    else:
        deck.write_text(expected)
        with (path / f"{case}_run.log").open("w") as log:
            subprocess.run([str(exe), deck.name], cwd=path,
                           env={**os.environ, "OMP_NUM_THREADS": "1"},
                           stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
    log = (path / f"{case}_run.log").read_text()
    if SUCCESS not in log or REVISION not in log or re.search(r"NaN|Infinity|Error:", log):
        raise ValueError(f"incomplete/wrong-version native run: {path}")
    summary = path / f"{case}_summary_01_0001.csv"
    _, history = numeric_table(summary)
    if abs(history[-1][0] - p["end"]) > p["dt"]/time + 5e-5:
        raise ValueError(f"missing final observation: {path}")
    if case == "vrx01":
        mass_error = max(abs(row[4]-1) for row in history)
        temperature_error = max(abs(row[5]+273.15-400) for row in history)
        if mass_error > 1e-6 or temperature_error > 1e-6:
            raise ValueError(f"0D mass/temperature invariant failed: {path}")
        grid = dict(intervals=1, spatial_unknowns=1, coordinates_m=[],
                    spacing_min_m=p["length"], spacing_max_m=p["length"],
                    coordinate_convention="ZEROD homogeneous kinetics; one lumped unknown")
    else:
        profile = path / f"{case}_profile_01_01_TEMPERATURE.csv"
        header, data = numeric_table(profile)
        z = [float(x) for x in header[1:] if x.strip()]
        intervals = space*p["cells"]
        if len(z) != intervals+1 or any(len(row) != len(z)+1 for row in data):
            raise ValueError(f"truncated profile or wrong cell count: {path}")
        h = p["length"]/intervals
        if any(not math.isclose(x, i*h, abs_tol=1e-10) for i, x in enumerate(z)):
            raise ValueError(f"native coordinates do not match requested mesh: {path}")
        dz = [b-a for a, b in zip(z, z[1:])]
        grid = dict(intervals=intervals, spatial_unknowns=len(z), coordinates_m=z,
                    spacing_min_m=min(dz), spacing_max_m=max(dz),
                    coordinate_convention="boundary nodes; half control volumes at both boundaries")
    raw = HERE / "raw" / level / case
    raw.mkdir(parents=True, exist_ok=True)
    files = {}
    inputs = [deck, path / f"{case}_run.log", summary]
    if case != "vrx01":
        inputs.append(profile)
    for src in inputs:
        dst = raw / src.name
        shutil.copyfile(src, dst)
        files[str(dst.relative_to(HERE))] = digest(dst)
    print(f"{case} {level}: successful V2 run; {grid['spatial_unknowns']} unknowns", flush=True)
    return dict(case=case, case_id=p["id"], level=level, dt=p["dt"]/time,
                observation_end_s=p["end"], native_stop_s=p["end"]+p["dt"]/time,
                grid=grid, successful_termination=True,
                rejected_step_messages=len(re.findall("Failed timestep", log)),
                timestep_interpretation="requested initial step; native reductions permitted; console is a sparse step log",
                files=files)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("--extract-only", action="store_true")
    args = parser.parse_args()
    exe, work = args.executable.resolve(), args.workspace.resolve()
    jobs = [(case, level) for level in LEVELS for case in CASES]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        runs = list(pool.map(lambda job: run_one(job, exe, work, args.extract_only), jobs))
    extracted = {}
    metrics = []
    for level in LEVELS:
        dst = HERE if level == "base" else HERE/level
        dst.mkdir(parents=True, exist_ok=True)
        metric_path = work/level/"accuracy.csv"
        result = subprocess.run([sys.executable, str(ARCHIVE / "extract_results.py"), str(work/level),
                                 "--metrics", str(metric_path)],
                                text=True, capture_output=True, check=True)
        for row in csv.DictReader(metric_path.open()):
            if not math.isfinite(float(row["error"])):
                raise ValueError(f"nonfinite analytical error: {row}")
            metrics.append(dict(row, resolution=level))
        (dst/"analytical_checks.txt").write_text(result.stdout)
        print(level + " analytical checks:\n" + result.stdout, flush=True)
        for src in sorted((work/level/"gpyro_out").glob("*.csv")):
            # Normalize legacy CSV line endings, preserving its SI schema.
            target = dst/src.name
            with src.open() as stream, target.open("w", newline="") as out:
                csv.writer(out, lineterminator="\n").writerows(csv.reader(stream))
            extracted[str(target.relative_to(HERE))] = digest(target)
    with (HERE/"accuracy.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(metrics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)
    extracted["accuracy.csv"] = digest(HERE/"accuracy.csv")
    provenance = dict(source_repository="https://github.com/Gedeben/GPYRO-V2",
                      revision=REVISION, executable_sha256=digest(exe),
                      compiler="GNU Fortran 14.4.0; upstream gnu_linux flags",
                      source_modifications=False, threads=1, timing_measurements_retained=False,
                      generator_sha256=digest(Path(__file__)), extractor_sha256=digest(ARCHIVE/"extract_results.py"),
                      archived_inputs_sha256={c+".data": digest(ARCHIVE/(c+".data")) for c in CASES},
                      extracted_sha256=extracted, runs=runs)
    (HERE/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")


if __name__ == "__main__":
    main()
