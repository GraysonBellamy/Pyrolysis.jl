"""Run and extract matched accuracy comparisons; see README.md for mappings.

Python 3.12+, standard library only. No solver source changes or timing claims.
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

DEST = Path(__file__).resolve().parent
CASES = {
    "gpyro": [
        "vhc02a",
        "vhc02b",
        "vhc02c",
        "vhc02d",
        "vhc04a",
        "vhc04b",
        "vhc04c",
        "vhc05",
        "vrx02a",
        "vrx02b",
        "vrx02c",
        "vrx04a",
        "vrx04b",
        "vcpl03p1",
        "vcpl03p2",
        "vcpl03p3",
        "vcpl04",
    ],
    "fds": [
        "vhc04a",
        "vhc04b",
        "vhc04c",
        "vhc05",
        "vrx02a",
        "vrx02b",
        "vrx02c",
        "vrx04a",
        "vrx04b",
        "vrx05a",
        "vrx05b",
        "vcpl03p1",
        "vcpl03p2",
        "vcpl03p3",
        "vcpl04",
    ],
}
VERSIONS = {
    "gpyro": "bff0126c68f2a35c084ef04490ae327185f281a3",
    "fds": "FDS-6.11.1-0-gff928db-release",
}

# Keep the original coupled refinement for reproducibility. The other three
# runs isolate time error at the base mesh, then space error at a fixed dt.
LEVELS = {
    "base": (0, 0),
    "refined": (1, 1),
    "time_half": (0, 1),
    "time_quarter": (0, 2),
    "space_half": (1, 2),
}


def refinement(level):
    return LEVELS[level] if isinstance(level, str) else (level, level)


def parameters(case):
    p = dict(
        case=case,
        L=0.02,
        rho=1000.0,
        cp=1000.0,
        k=0.5,
        ti=300.0,
        end=10000.0,
        dt=2.0,
        nodes=81,
        dump=100.0,
        eps=0.0,
        front="fixed",
        tf=400.0,
        back="fixed",
        tb=300.0,
    )
    if case.startswith("vhc02"):
        bi = dict(a=0.1, b=1.0, c=10.0, d=100.0)[case[-1]]
        p.update(
            L=0.1,
            front="convective",
            tf=500.0,
            h=bi * 0.5 / 0.1,
            back="insulated",
            dump=1000.0,
        )
    elif case.startswith("vhc04"):
        p.update(
            k=0.1 if case[-1] != "b" else 0.01,
            interface=0.0102 if case[-1] == "c" else 0.01,
        )
        p.update(end=10 * 0.02**2 * 1000000.0 / p["k"], dt=1 / p["k"])
        p["dump"] = p["end"] / 2
    elif case == "vhc05":
        p.update(tf=800.0, end=5000.0, dt=1.0, dump=2500.0)
    elif case == "vrx04a":
        p.update(
            L=0.0001,
            k=100.0,
            ti=900.0,
            front="radiative",
            tf=300.0,
            back="insulated",
            eps=0.9,
            h=0.0,
            q=0.0,
            end=8.0,
            dt=0.002,
            dump=0.1,
            nodes=5,
        )
    elif case == "vrx04b":
        p.update(
            L=0.005,
            k=1.0,
            front="radiative",
            tf=300.0,
            back="insulated",
            eps=0.85,
            h=15.0,
            q=20000.0,
            end=2000.0,
            dt=0.2,
            dump=20.0,
            nodes=21,
        )
    elif case.startswith("vrx02"):
        rho, n, a, end, temp = (
            (100.0, 0.0, 0.05, 18.0, 400.0)
            if case[-1] == "a"
            else (100.0, 2.0, 0.1, 100.0, 400.0)
            if case[-1] == "b"
            else (
                1000.0,
                2 / 3,
                340000000000.0 * math.exp(-103000 / (8.314462618 * 450)),
                None,
                450.0,
            )
        )
        if end is None:
            end = 0.9 * 3 / a
        p.update(
            L=0.001,
            rho=rho,
            n=n,
            a=a,
            end=end,
            ti=temp,
            tf=temp,
            tb=temp,
            dt=end / 10000,
            dump=end / 100,
            nodes=2,
        )
    elif case.startswith("vrx05"):
        p.update(
            L=0.001,
            k=0.2,
            cp=1200.0,
            ti=1000.0,
            front="insulated",
            back="insulated",
            end=20.0,
            dt=0.001,
            dump=0.5,
            nodes=2,
            n=1.0,
            a=0.2,
            heat=1000000.0 if case[-1] == "a" else -100000.0,
        )
    elif case.startswith("vcpl03"):
        p.update(
            L=0.002,
            k=0.5,
            ti=400.0,
            front="insulated",
            back="insulated",
            end=8.0,
            dt=0.001,
            dump=0.1,
            nodes=9,
            n=1.0,
            a=2.0,
            product_rho=500.0
            if case[-1] == "1"
            else 2000.0
            if case[-1] == "2"
            else 1000.0,
        )
    elif case == "vcpl04":
        p.update(
            L=0.02,
            k=0.2,
            ti=500.0,
            tf=500.0,
            tb=500.0,
            end=2500.0,
            dt=0.25,
            dump=125.0,
            nodes=81,
            n=0.0,
            a=0.44 * 0.2 * 12.5 / 0.01**2 / 100000.0 * math.exp(40) / 1000,
            heat=-100000.0,
            activation=40 * 8.314462618 * 500,
        )
    else:
        raise ValueError(case)
    return p


def section(name, lines):
    return "&GPYRO_" + name + "\n" + "\n".join(lines) + "\n/\n\n"


def gpyro(p, level):
    case = p["case"]
    space, time = refinement(level)
    refine = 2**space
    L, rho, cp, k, ti, end, dt, nodes, dump = (
        p[x] for x in ("L", "rho", "cp", "k", "ti", "end", "dt", "nodes", "dump")
    )
    dt /= 2**time
    nodes = (nodes - 1) * refine + 1
    reaction = case.startswith(("vrx02", "vrx05", "vcpl03", "vcpl04"))
    layered = case.startswith("vhc04")
    affine = case == "vhc05"
    ns = 1 if case == "vcpl04" else 2 if reaction or layered or affine else 1
    general = [
        f"DT0={dt:.17g}, TAMB={ti}, TREF=300, P0=101325,",
        "GX=0, GY=0, GZ=0, NTDMA_ITERATIONS=2000, TMPTOL=1D-8, HTOL=1D-12, YITOL=1D-10,",
        "SOLVE_PRESSURE=.FALSE., SOLVE_GAS_YJ=.FALSE., SOLVE_GAS_ENERGY=.FALSE., THERMAL_EQUILIBRIUM=.TRUE.,",
        "CONVENTIONAL_RXN_ORDER=.TRUE., CONSTANT_DHVOL=.TRUE., SHYI_CORRECTION=.FALSE.,",
        f"SOLVER_DEFORMATION_MODE={(3 if case.startswith('vcpl03') else 1)}, FULL_QSG=.FALSE., BLOWING=.FALSE.,",
    ]
    points = [("TEMPERATURE", 0, 0), ("TEMPERATURE", 0, L), ("THICKNESS", 0, 0)]
    if reaction:
        points += [("CI", 1, 0)] + ([] if case == "vcpl04" else [("CI", 2, 0)])
    output = [
        f"CASENAME='{case}', N_POINT_QUANTITIES={len(points)}, N_PROFILE_QUANTITIES=1,",
        f"DTDUMP_POINT={dump - 1e-08:.17g}, DTDUMP_PROFILE={dump - 1e-08:.17g}, DTMIN_KILL=1D-12,",
        "PROFILE_QUANTITY(1)='TEMPERATURE', PROFILE_QUANTITY_INDEX(1)=0, PROFILE_DIRECTION(1)='z', PROFILE_ISKIP(1)=1,",
    ]
    for i, (q, ix, z) in enumerate(points, 1):
        output.append(
            f"POINT_QUANTITY({i})='{q}', POINT_QUANTITY_INDEX({i})={ix}, POINT_Z({i})={z},"
        )
    props = [f"NSSPEC={ns},"]
    for i in range(1, ns + 1):
        kval = (
            1.0
            if layered and i == 2
            else 0.4
            if affine and i == 1
            else 1.2
            if affine
            else k
        )
        r = p.get("product_rho", rho) if i == 2 else rho
        c = 1900.0 if case.startswith("vrx05") and i == 2 else cp
        props.append(
            f"NAME({i})='{chr(64 + i)}', R0({i})={r}, RS0({i})={r}, C0({i})={c}, K0Z({i})={kval}, NKZ({i})={(1 if affine and i == 2 else 0)}, EMIS({i})={p['eps']},"
        )
    if layered:
        general.append("GEOMETRY_IS_UPSIDE_DOWN=.FALSE.,")
    rx = ["NRXNS=0,"]
    if reaction:
        rx = [
            "NRXNS=1,",
            f"CFROM(1)='A', CTO(1)='{('A' if case == 'vcpl04' else 'B')}', Z(1)={p['a']:.17g}, E(1)={p.get('activation', 0) * (0.00831400015205145 / 8.314462618)}, ORDER(1)={p['n']}, CHI(1)={(0 if case.startswith('vcpl03') else 1)}, DHS(1)={p.get('heat', 0)}, DHV(1)={p.get('heat', 0)},",
        ]
    ic = [
        f"NIC={(2 if layered else 1)}, TMP_INITIAL(1)={ti}, TMPG_INITIAL(1)={ti}, P_INITIAL(1)=101325, YI0(1,1)={(0.5 if affine else 1)}, YJ0(1,1)=1,"
    ]
    if affine:
        ic.append("YI0(1,2)=0.5,")
    if layered:
        ic.append(
            f"TMP_INITIAL(2)={ti}, TMPG_INITIAL(2)={ti}, P_INITIAL(2)=101325, YI0(2,2)=1, YJ0(2,1)=1,"
        )
    bc = ["NSURF_IDX=2,"]
    for i, face in enumerate(("front", "back"), 1):
        mode = p[face]
        temp = p["tf" if i == 1 else "tb"]
        h = p.get("h", 0) if i == 1 else 0
        q = p.get("q", 0) / p["eps"] if mode == "radiative" else 0
        bc.append(
            f"SURF_IDX({i})={i}, QE({i})={q}, HC({i})={h}, RERADIATION({i})={('.TRUE.' if mode == 'radiative' else '.FALSE.')}, TINF({i})={temp}, TINFG({i})={temp}, TFIXED({i})={(temp if mode == 'fixed' else -1)}, PRES({i})=-1, MDOTPP({i})=0, HM({i})=0, YJINF({i},1)=1,"
        )
    geom = [
        f"NMESH=1, ZDIM(1)={L}, NCELLZ(1)={nodes}, NCELLX(1)=1, NCELLY(1)=1, HALF_CELLS_AT_BC(1)=.TRUE., DEFAULT_IC(1)=1, DEFAULT_SURF_IDX(1,5)=1, DEFAULT_SURF_IDX(1,6)=2,"
    ]
    if layered:
        geom.append(
            f"NOBST=2, IMESH(1)=1, Z1(1)=0, Z2(1)={L - p['interface']}, X1(1)=0, X2(1)=1, Y1(1)=0, Y2(1)=1, ICNUM(1)=1, SURF_IDX2D(1,5)=1,"
        )
        geom.append(
            f"IMESH(2)=1, Z1(2)={L - p['interface']}, Z2(2)={L}, X1(2)=0, X2(2)=1, Y1(2)=0, Y2(2)=1, ICNUM(2)=2, SURF_IDX2D(2,6)=2,"
        )
    sections = [
        ("GENERAL", general),
        ("OUTPUT", output),
        ("SPROPS", props),
        ("GPROPS", ["NGSPEC=1, IBG=1, IO2=1, CPG=1000, NAME(1)='gas', M(1)=28,"]),
        ("RXNS", rx),
        ("HGRXNS", ["NHGRXNS=0,"]),
        ("GYIELDS", ["GYIELDS(1,1)=1,"] if reaction else []),
        ("HGYIELDS", []),
        (
            "CASES",
            [f"NCASES=1, IMESH(1)=1, TSTOP(1)={end + dt:.17g}, ZEROD(1)=.FALSE.,"],
        ),
        ("IC", ic),
        ("ALLBC", bc),
        ("GEOM", geom),
    ]
    return "".join((section(*s) for s in sections))


def fds(p, level):
    case = p["case"]
    space, time = refinement(level)
    refine = 2**space
    L, rho, cp, k, ti, end, dt, nodes, dump = (
        p[x] for x in ("L", "rho", "cp", "k", "ti", "end", "dt", "nodes", "dump")
    )
    dt /= 2**time
    lines = [
        f"&HEAD CHID='{case}', TITLE='Matched Pyrolysis.jl verification problem' /",
        "&MESH IJK=3,3,3, XB=-0.15,0.15,-0.15,0.15,0,0.3 /",
        f"&TIME T_END={end:.17g}, DT={dt:.17g}, WALL_INCREMENT=1 /",
        f"&MISC SOLID_PHASE_ONLY=T, TMPA={(p['tf'] if p['front'] == 'radiative' else ti) - 273.15} /",
        f"&RADI RADIATION={('T' if p['front'] == 'radiative' else 'F')} /",
        f"&DUMP DT_DEVC={dump:.17g}, DT_PROF={end:.17g} /",
    ]
    reaction = case.startswith(("vrx02", "vrx05", "vcpl03", "vcpl04"))
    layered = case.startswith("vhc04")
    ns = 1 if case == "vcpl04" else 2 if reaction or layered else 1
    for i in range(1, ns + 1):
        kval = 1.0 if layered and i == 2 else k
        c = 1900.0 if case.startswith("vrx05") and i == 2 else cp
        density = p.get("product_rho", rho) if i == 2 else rho
        props = f"&MATL ID='{chr(64 + i)}', DENSITY={density}, SPECIFIC_HEAT={c / 1000}, CONDUCTIVITY={kval}, EMISSIVITY={p['eps']}"
        if case.startswith("vrx05") or case == "vcpl04":
            props += ", ADJUST_H=F"
        if case == "vhc05":
            props += ", CONDUCTIVITY_RAMP='k'"
        if reaction and i == 1:
            props += f", A={p['a'] * rho ** (1 - p['n']):.17g}, E={p.get('activation', 0) * 8.314472 / 8.314462618}, N_S={p['n']}, NU_MATL=1, MATL_ID='{('A' if case == 'vcpl04' else 'B')}', HEAT_OF_REACTION={p.get('heat', 0) / 1000}"
        lines.append(props + " /")
    if case == "vhc05":
        for temp in (300, 800):
            lines.append(
                f"&RAMP ID='k', T={temp - 273.15}, F={0.2 * (1 + 0.01 * temp)} /"
            )
    surf = "&SURF ID='wall', MATL_ID=" + ("'A','B'" if layered else "'A'")
    surf += ", THICKNESS=" + (
        f"{L - p['interface']},{p['interface']}" if layered else str(L)
    )
    # These are per-layer arrays. A scalar sets only layer 1; leaving layer 2
    # at its defaults silently produced a much coarser, stretched back layer.
    for layer in range(1, 3 if layered else 2):
        surf += f", STRETCH_FACTOR({layer})=1, CELL_SIZE({layer})={L / ((nodes - 1) * refine):.17g}"
    surf += f", HEAT_TRANSFER_COEFFICIENT={p.get('h', 0)}"
    if p["front"] == "fixed":
        surf += f", TMP_FRONT={p['tf'] - 273.15}"
    if p["front"] == "convective":
        surf += f", TMP_GAS_FRONT={p['tf'] - 273.15}"
    if p["front"] == "radiative":
        surf += f", TMP_INNER={ti - 273.15}, SKIP_INRAD=T, EXTERNAL_FLUX={(p['q'] / p['eps'] + 5.670373e-08 * p['tf'] ** 4) / 1000}"
    surf += (
        f", TMP_BACK={p['tb'] - 273.15}"
        if p["back"] == "fixed"
        else ", BACKING='INSULATED'"
    )
    lines += [
        surf + " /",
        "&VENT XB=-0.05,0.05,-0.05,0.05,0,0, SURF_ID='wall' /",
        "&PROF XYZ=0,0,0, IOR=3, QUANTITY='TEMPERATURE', CELL_CENTERED=T, FORMAT_INDEX=2, ID='temperature' /",
    ]
    for ident, quantity, extras in [
        ("front", "WALL TEMPERATURE", ""),
        ("back", "BACK WALL TEMPERATURE", ""),
        ("mean", "INSIDE WALL TEMPERATURE", f", DEPTH={L / 2}"),
    ]:
        lines.append(
            f"&DEVC ID='{ident}', XYZ=0,0,0, IOR=3, QUANTITY='{quantity}', TIME_AVERAGED=F{extras} /"
        )
    if case == "vrx04a":
        # The same volume-mean observable as Julia and Gpyro: native cell
        # temperatures weighted by cell volumes. Surface/depth trapezoids
        # would add a different observation/quadrature error to the comparison.
        count = (nodes - 1) * refine
        for j in range(1, count + 1):
            lines.append(
                f"&DEVC ID='cell_{j}', XYZ=0,0,0, IOR=3, QUANTITY='INSIDE WALL TEMPERATURE', DEPTH={L * (j - .5) / count}, TIME_AVERAGED=F /"
            )
    if reaction:
        for name in ("A",) if case == "vcpl04" else ("A", "B"):
            lines.append(
                f"&DEVC ID='{name}', XYZ=0,0,0, IOR=3, QUANTITY='SOLID DENSITY', MATL_ID='{name}', DEPTH={L / 2}, TIME_AVERAGED=F /"
            )
    if case.startswith("vcpl03"):
        lines.append(
            "&DEVC ID='thickness', XYZ=0,0,0, IOR=3, QUANTITY='WALL THICKNESS', TIME_AVERAGED=F /"
        )
    lines.append("&TAIL /")
    return "\n".join(lines) + "\n"


def run(job, executable, workspace, extract_only):
    solver, case, level = job
    p = parameters(case)
    path = workspace / solver / level / case
    path.mkdir(parents=True, exist_ok=True)
    deck = path / (case + (".fds" if solver == "fds" else ".data"))
    expected = (fds if solver == "fds" else gpyro)(p, level)
    if extract_only:
        if deck.read_text() != expected:
            raise RuntimeError(f"{deck}: input differs from generated deck; rerun")
    else:
        deck.write_text(expected)
        env = {**os.environ, "OMP_NUM_THREADS": "1"}
        if solver == "fds":
            env["I_MPI_FABRICS"] = "shm"
            mpi = executable.parent / "intelmpi"
            if mpi.is_dir():
                env["LD_LIBRARY_PATH"] = (
                    str(mpi / "lib") + ":" + env.get("LD_LIBRARY_PATH", "")
                )
                env["FI_PROVIDER_PATH"] = str(mpi / "prov")
        with (path / "run.log").open("w") as log:
            subprocess.run(
                [str(executable), deck.name],
                cwd=path,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=600,
                check=True,
            )
    log = (path / "run.log").read_text(errors="replace")
    success = (
        "FDS completed successfully"
        if solver == "fds"
        else "End of simulation reached successfully."
    ) in log
    if not success:
        raise RuntimeError(f"{job}: unsuccessful solver termination: {log[-500:]}")
    banner = (path / f"{case}.out").read_text() if solver == "fds" else log
    if VERSIONS[solver] not in banner:
        raise RuntimeError(f"{job}: expected pinned solver version {VERSIONS[solver]}")
    return extract(solver, case, level, path, deck)


def numeric_csv(path, skip):
    with path.open() as stream:
        rows = list(csv.reader(stream))[skip:]
    data = [
        [float(x) for x in row if x.strip()]
        for row in rows
        if any((x.strip() for x in row))
    ]
    if not data or not all((math.isfinite(x) for row in data for x in row)):
        raise RuntimeError(f"{path}: missing or nonfinite data")
    if len({len(row) for row in data}) != 1:
        raise RuntimeError(f"{path}: inconsistent columns")
    return data


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_grid(solver, case, path):
    """Actual initial grid, not device/probe count or requested CELL_SIZE."""
    if solver == "fds":
        text = (path / f"{case}.out").read_text()
        block = text.split("Solid Phase Node, Layer, Coordinates(m):")[1]
        coordinates = []
        layers = []
        for line in block.splitlines()[1:]:
            match = re.fullmatch(r"\s*(\d+)\s+(\d+)\s+([\d.Ee+-]+)\s*", line)
            if not match:
                break
            coordinates.append(float(match[3]))
            layers.append(int(match[2]))
        convention = "cell faces in mesh table; native thermal unknowns at centers; matched PROF uses CELL_CENTERED=T"
        unknowns = len(coordinates) - 1
    else:
        profile = path / f"{case}_profile_01_01_TEMPERATURE.csv"
        coordinates = numeric_csv(profile, 0)[0][1:]
        layers = []
        convention = "boundary-centered nodes; half control volumes at boundaries"
        unknowns = len(coordinates)
    widths = [b-a for a, b in zip(coordinates, coordinates[1:])]
    if not widths or min(widths) <= 0:
        raise RuntimeError(f"{case}: invalid native grid")
    return dict(coordinates_m=coordinates, layer_indices=layers,
                coordinate_convention=convention, spatial_unknowns=unknowns,
                intervals=len(widths), spacing_min_m=min(widths), spacing_max_m=max(widths))


def case_id(case):
    if case.startswith("vcpl"):
        return "V-CPL-" + case[4:6]
    return f"V-{case[1:3].upper()}-{case[3:]}"


def extract(solver, case, level, path, deck):
    p = parameters(case)
    space, time = refinement(level)
    refined = level
    label = "FDS 6.11.1" if solver == "fds" else "Gpyro V2"
    native = path / (
        f"{case}_devc.csv" if solver == "fds" else f"{case}_summary_01_0001.csv"
    )
    data = numeric_csv(native, 2 if solver == "fds" else 1)
    tolerance = max(p["dt"], 0.0001) * 1.1
    if abs(data[-1][0] - p["end"]) > tolerance:
        raise RuntimeError(f"{native}: final time {data[-1][0]} != {p['end']}")
    if any((b[0] <= a[0] for a, b in zip(data, data[1:]))):
        raise RuntimeError(f"{native}: nonmonotone output times")
    files = [deck, native]
    samples = []
    metadata = []

    def add(curve, reference, xy, suffix=""):
        index = len(metadata) + 1
        metadata.append(
            dict(
                case=case_id(case),
                curve=curve,
                reference=reference,
                label=label + suffix,
                path=f"{solver}/{case}.csv",
                index=index,
            )
        )
        samples.extend(((index, float(x), float(y)) for x, y in xy))

    profile = None
    if case.startswith(("vhc", "vcpl04")) or case == "vrx04a":
        profile = path / (
            f"{case}_prof_1.csv"
            if solver == "fds"
            else f"{case}_profile_01_01_TEMPERATURE.csv"
        )
        files.append(profile)
        prof = numeric_csv(profile, 2 if solver == "fds" else 0)
    if case.startswith("vhc02"):
        depth = prof[0][1:]
        for fo in (0.05, 0.1, 0.2, 0.5):
            t = fo * 0.1**2 / (0.5 / 1000000.0)
            row = min(prof[1:], key=lambda r: abs(r[0] - t))
            if abs(row[0] - t) > tolerance:
                raise RuntimeError(f"{profile}: missing t={t}")
            add(
                "T_profiles",
                f"Fo = {fo}",
                sorted(
                    (
                        ((p["L"] - d) / p["L"], v + 273.15)
                        for d, v in zip(depth, row[1:])
                    )
                ),
                f", Fo = {fo}",
            )
        add(
            "surface_T_history",
            "front face",
            ((r[0], r[1] + 273.15) for r in data if r[0] > 0),
        )
        add(
            "back_T_history",
            "adiabatic back face",
            ((r[0], r[2] + 273.15) for r in data if r[0] > 0),
        )
    elif case.startswith(("vhc04", "vhc05", "vcpl04")):
        depth, temp = zip(*prof) if solver == "fds" else (prof[0][1:], prof[-1][1:])
        if solver == "gpyro" and abs(prof[-1][0] - p["end"]) > tolerance:
            raise RuntimeError(f"{profile}: final profile missing")
        add(
            "steady_profiles" if case == "vcpl04" else "T_profile_steady",
            "δ = 0.44 ≈ 0.50 δ_c" if case == "vcpl04" else "steady solver profile",
            sorted(((p["L"] - d, t + 273.15) for d, t in zip(depth, temp))),
            ", δ = 0.44" if case == "vcpl04" else "",
        )
    elif case.startswith("vrx02"):
        n = "0.0" if case[-1] == "a" else "2.0" if case[-1] == "b" else str(2 / 3)
        add("xi_histories", f"ξ_A (n = {n})", ((r[0], r[4]) for r in data))
    elif case == "vrx04a":
        if solver == "gpyro":
            depth = prof[0][1:]
            means = [
                (
                    r[0],
                    sum(
                        (
                            (b - a) * (u + v) / 2
                            for a, b, u, v in zip(depth, depth[1:], r[1:], r[2:])
                        )
                    )
                    / p["L"]
                    + 273.15,
                )
                for r in prof[1:]
            ]
        else:
            m = (p["nodes"] - 1) * 2**space
            means = [
                (r[0], sum(r[4:]) / m + 273.15) for r in data
            ]
            if len(prof) != m or abs(means[-1][1] - (sum(r[1] for r in prof)/m + 273.15)) > 5e-5:
                raise RuntimeError(f"{case}: mean-temperature probes disagree with native cell-center profile")
        add("T_history_vs_exact", "slab mean temperature", means)
    elif case == "vrx04b":
        add(
            "surface_T_history",
            "surface temperature",
            ((r[0], r[1] + 273.15) for r in data),
        )
    elif case.startswith("vrx05"):
        add("temperature", "temperature", ((r[0], r[3] + 273.15) for r in data))
        for name, col in [("A", 4), ("B", 5)]:
            add("species", name, ((r[0], r[col]) for r in data), " " + name)
        add(
            "enthalpy",
            "enthalpy",
            (
                (r[0], p["L"] * (1200 * r[4] + 1900 * r[5]) * (r[3] + 273.15 - 298.15))
                for r in data
            ),
        )
    elif case.startswith("vcpl03"):
        perm = "P" + case[-1]
        thickness_column = 6 if solver == "fds" else 3
        add(
            "thickness_histories",
            perm,
            ((r[0], r[thickness_column] / p["L"]) for r in data),
            " " + perm,
        )
    else:
        raise ValueError(case)
    if not samples or not all((math.isfinite(v) for row in samples for v in row)):
        raise RuntimeError(f"{case}: invalid extracted samples")
    if case == "vhc04c":
        # Same physical 10.2 mm interface; only the Julia mesh changes in d.
        metadata.append(dict(
            metadata[0], case="V-HC-04d",
            reuse_note="Same physical problem as V-HC-04c; Julia aligns its mesh interface",
        ))
    out = DEST / ("" if level == "base" else level) / solver / f"{case}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["series", "x", "y"])
        writer.writerows(samples)
    raw = DEST / "raw" / solver / refined / case
    raw.mkdir(parents=True, exist_ok=True)
    grid = native_grid(solver, case, path)
    expected_intervals = (p["nodes"] - 1) * 2**space
    if grid["intervals"] != expected_intervals:
        raise RuntimeError(f"{solver} {case}: requested {expected_intervals} intervals, "
                           f"native grid has {grid['intervals']}")
    grid["declared_domain_m"] = p["L"]
    grid["coordinates_span_domain"] = math.isclose(
        grid["coordinates_m"][-1], p["L"], abs_tol=5e-8, rel_tol=1e-8)
    if not grid["coordinates_span_domain"]:
        # Upstream COMPUTE_Z_CELL_POSITION's first/last-node formulas overlap
        # when NCELLZ=2: it prints z=[0,3L/4], although both DZN=L/2.
        # These spatially homogeneous ODE cases still have total volume L.
        # Preserve the actual coordinates and expose the exception explicitly.
        if (solver == "gpyro" and grid["spatial_unknowns"] == 2 and
                case.startswith("vrx02") and
                math.isclose(grid["coordinates_m"][-1], .75*p["L"], rel_tol=1e-8)):
            grid["coordinate_exception"] = "native two-node position formula: printed span 3L/4; control-volume sum L; homogeneous ODE only"
        else:
            raise RuntimeError(f"{solver} {case}: native grid does not span the physical domain")
    grid_file = raw / "grid.json"
    grid_file.write_text(json.dumps(grid, indent=2) + "\n")
    hashes = {}
    hashes[str(grid_file.relative_to(DEST))] = sha256(grid_file)
    for f in files:
        dest = raw / f.name
        shutil.copyfile(f, dest)
        hashes[str(dest.relative_to(DEST))] = sha256(dest)
    hashes[str(out.relative_to(DEST))] = sha256(out)
    print(f"{solver} {case} {refined}: verified {len(samples)} samples", flush=True)
    return (
        dict(
            case=case,
            level=refined,
            dt=p["dt"] / 2**time,
            requested_intervals=(p["nodes"] - 1) * 2**space,
            grid=grid,
            observed_step_history="not retained; output intervals are not integration steps",
            timestep_interpretation=("outer step; native solid substepping permitted"
                                     if solver == "fds" else
                                     "requested initial step; native step reductions permitted"),
            successful_termination=True,
            files=hashes,
        ),
        metadata,
    )


def main():
    global DEST
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", choices=CASES)
    parser.add_argument("executable", type=Path)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("--extract-only", action="store_true")
    parser.add_argument("--levels", default=",".join(LEVELS),
                        help="comma-separated levels: " + ",".join(LEVELS))
    parser.add_argument("--cases", default="all")
    parser.add_argument("--destination", type=Path, default=DEST,
                        help="separate destination for a qualification study")
    args = parser.parse_args()
    executable = args.executable.resolve()
    DEST = args.destination.resolve()
    DEST.mkdir(parents=True, exist_ok=True)
    levels = args.levels.split(",")
    cases = CASES[args.solver] if args.cases == "all" else args.cases.split(",")
    if not levels or any(level not in LEVELS for level in levels):
        parser.error("unknown refinement level")
    if not cases or any(case not in CASES[args.solver] for case in cases):
        parser.error("unknown case")
    if DEST == Path(__file__).resolve().parent and (
        set(cases) != set(CASES[args.solver]) or set(levels) != set(LEVELS)
    ):
        parser.error("partial studies require --destination to preserve the complete reference manifest")
    jobs = [(args.solver, c, level) for c in cases for level in levels]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(
            pool.map(
                lambda job: run(
                    job, executable, args.workspace.resolve(), args.extract_only
                ),
                jobs,
            )
        )
    manifest = DEST / f"{args.solver}.toml"
    with manifest.open("w") as stream:
        for audit, series in results:
            if audit["level"] != "base":
                continue
            for s in series:
                stream.write("[[series]]\n")
                for key, value in s.items():
                    stream.write(f"{key} = {json.dumps(value, ensure_ascii=False)}\n")
                stream.write("\n")
    provenance = dict(
        solver=args.solver,
        revision=VERSIONS[args.solver],
        executable_sha256=sha256(executable),
        threads=1,
        timing_measurements_retained=False,
        source_modifications=False,
        generator_sha256=sha256(Path(__file__)),
        manifest_sha256=sha256(manifest),
        runs=[r[0] for r in results],
    )
    (DEST / f"{args.solver}_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
