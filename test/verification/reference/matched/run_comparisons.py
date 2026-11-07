"""Run and extract matched accuracy comparisons; see README.md for mappings.

Python 3.12+, standard library only. No solver source changes; scripts/timing/
times these decks.
"""

from pathlib import Path
import argparse
import concurrent.futures
import csv
import gzip
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
        "vhc01",
        "vhc02a",
        "vhc02b",
        "vhc02c",
        "vhc02d",
        "vhc03",
        "vhc04a",
        "vhc04b",
        "vhc04c",
        "vhc05",
        "vrx01",
        "vrx02a",
        "vrx02b",
        "vrx02c",
        "vrx04a",
        "vrx04b",
        "vrx05a",
        "vrx05b",
        "vrx05c",
        "vcpl01",
        "vcpl03p1",
        "vcpl03p2",
        "vcpl03p3",
        "vcpl04",
        "vgt07a",
        "vgt07b",
        "vgt08a",
        "vgt08b",
    ],
    "fds": [
        "vhc01",
        "vhc02a",
        "vhc02b",
        "vhc02c",
        "vhc02d",
        "vhc04a",
        "vhc04b",
        "vhc04c",
        "vhc05",
        "vrx01",
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
    "thermakin": [
        "vhc01",
        "vhc02a",
        "vhc02b",
        "vhc02c",
        "vhc02d",
        "vhc03",
        "vhc04a",
        "vhc04b",
        "vhc04c",
        "vhc05",
        "vrx01",
        "vrx02b",
        "vrx04a",
        "vrx04b",
        "vrx05a",
        "vrx05b",
        "vrx05c",
        "vcpl03p1",
        "vcpl03p2",
        "vcpl03p3",
        "vcpl04",
        "vgt01",
        "vgt07a",
        "vgt07b",
    ],
}
VERSIONS = {
    "gpyro": "bff0126c68f2a35c084ef04490ae327185f281a3",
    "fds": "FDS-6.11.1-0-gff928db-release",
    "thermakin": "ThermaKin Program Version 7",
}
LABELS = {"fds": "FDS 6.11.1", "gpyro": "Gpyro V2", "thermakin": "ThermaKin2Ds 7"}

# ThermaKin's own constants (Constants.h). Arrhenius energies are rescaled so
# E/R matches; σ stays native, as for FDS and Gpyro.
TK_R = 8.31447
TK_SIGMA = 5.6704e-08
# A fixed-temperature face is a convective coefficient on the boundary
# element, whose initial temperature is the fixed value: Crank-Nicolson does
# not damp a large-h mode, so an initial jump would oscillate indefinitely.
TK_FIXED_H = 1e9
# Mass absorption coefficient (m²/kg) that saturates the exposed element's
# absorbed fraction, so external radiation is deposited at that element.
TK_ABSORPTION = 1e6
# Absorbed flux without re-radiation: the incident flux is q/ε at this
# emissivity, leaving an emission loss εσT⁴ below 1e-3 W/m² for T < 1200 K.
TK_FLUX_EMISSIVITY = 1e-8

# Keep the original coupled refinement for reproducibility. The other three
# runs isolate time error at the base mesh, then space error at a fixed dt.
LEVELS = {
    "base": (0, 0),
    "refined": (1, 1),
    "time_half": (0, 1),
    "time_quarter": (0, 2),
    "space_half": (1, 2),
    "time_double": (0, -1),
    "time_quadruple": (0, -2),
}
# The coarser steps look for a cheaper configuration of each code at its
# accuracy. A coarse configuration the code cannot express or complete is
# recorded as unavailable, with the reason, instead of failing the study.
COARSE = {"time_double", "time_quadruple"}


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
    if case == "vhc01":
        p.update(
            L=0.05,
            cp=1400.0,
            k=0.2,
            front="flux",
            q=25000.0,
            h=20.0,
            tf=300.0,
            back="insulated",
            end=180.0,
            dt=0.1,
            dump=5.0,
            nodes=501,
        )
    elif case == "vhc03":
        p.update(
            front="flux",
            q=5000.0,
            h=0.0,
            tf=300.0,
            back="insulated",
            end=1200.0,
            dt=0.1,
            dump=20.0,
            nodes=321,
        )
    elif case == "vrx01":
        p.update(
            L=0.001,
            ti=400.0,
            front="insulated",
            back="insulated",
            end=20.0,
            dt=0.001,
            dump=0.1,
            nodes=2,
            rates=(0.389, 0.262),
        )
    elif case == "vgt01":
        p.update(
            L=0.04,
            ti=400.0,
            front="insulated",
            back="insulated",
            end=8.0,
            dt=0.01,
            dump=0.5,
            nodes=161,
            diffusivity=4e-07 / 0.4,
            z0=0.02,
            sigma0=0.001,
            peak=0.1,
        )
    elif case.startswith("vgt07"):
        porosity, gas_density = 0.6, 100000.0 * 0.032 / (8.314462618 * 300.0)
        p.update(
            L=1.0,
            rho=400.0,
            cp=1200.0,
            k=0.2,
            front="insulated",
            back="insulated",
            end=400.0,
            dt=0.05,
            dump=10.0,
            nodes=161,
            diffusivity=2e-05,
            hm=0.0001 if case[-1] == "a" else 0.01,
            porosity=porosity,
            pore_gas_density=porosity * gas_density,
            ambient_fraction=0.22,
        )
    elif case.startswith("vhc02"):
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
        # ThermaKin's linearized step diverges for this k(T) at dt >= 0.1 s.
        p.update(tf=800.0, end=5000.0, dt=1.0, dump=2500.0, thermakin_dt=0.05)
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
            heat=dict(a=1000000.0, b=-100000.0, c=0.0)[case[-1]],
            solid_yield=0.4 if case[-1] == "c" else 1.0,
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
    elif case == "vcpl01":
        p.update(L=0.04, end=300.0, dt=0.02, nodes=801, melt_T=433.0)
    elif case.startswith("vgt08"):
        p.update(L=0.1, end=1.0, dt=0.005, nodes=81)
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


# Native V2 inputs for physics the parametric deck below does not express:
# a surface-absorbed flux, in-depth absorption, apparent-heat-capacity latent
# heat, a two-step chain, retained reaction gas, mass-transfer films and Darcy
# inflow. Each holds its base level; a level changes only DT0, NCELLZ(1) and
# TSTOP(1).
GPYRO_DECKS = DEST / "gpyro_decks"
GPYRO_TEMPLATED = {
    "vhc01", "vhc03", "vrx01", "vrx05a", "vrx05b", "vrx05c",
    "vcpl01", "vgt07a", "vgt07b", "vgt08a", "vgt08b",
}
# A homogeneous ZEROD run has one lumped unknown and no mesh to refine.
GPYRO_LUMPED = {"vrx01"}


def gpyro_templated(p, level):
    case = p["case"]
    space, time = refinement(level)
    dt = p["dt"] / 2**time
    nodes = 1 if case in GPYRO_LUMPED else (p["nodes"] - 1) * 2**space + 1
    text = (GPYRO_DECKS / f"{case}.data").read_text()
    base = 1 if case in GPYRO_LUMPED else p["nodes"]
    for key, value, expected in (
        ("DT0", f"{dt:.17g}", p["dt"]),
        ("NCELLZ(1)", str(nodes), base),
        ("TSTOP(1)", f"{p['end'] + dt:.17g}", p["end"] + p["dt"]),
    ):
        pattern = r"(?<![\w(])" + re.escape(key) + r"\s*=\s*([^,\s]+)"
        found = re.findall(pattern, text)
        if len(found) != 1 or not math.isclose(float(found[0]), expected, rel_tol=1e-12):
            raise ValueError(f"{case}: template {key} = {found}, expected one value {expected}")
        text = re.sub(pattern, f"{key}={value}", text)
    return text


def gpyro(p, level):
    case = p["case"]
    if case in GPYRO_TEMPLATED:
        return gpyro_templated(p, level)
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


def fds_transient_profiles(case):
    """FDS writes a profile at every output where the comparison reads
    transient profiles; otherwise only the final (steady) profile."""
    return case == "vhc01" or case.startswith("vhc02")


def fds(p, level):
    case = p["case"]
    space, time = refinement(level)
    refine = 2**space
    L, rho, cp, k, ti, end, dt, nodes, dump = (
        p[x] for x in ("L", "rho", "cp", "k", "ti", "end", "dt", "nodes", "dump")
    )
    dt /= 2**time
    transient = fds_transient_profiles(case)
    lines = [
        f"&HEAD CHID='{case}', TITLE='Matched Pyrolysis.jl verification problem' /",
        "&MESH IJK=3,3,3, XB=-0.15,0.15,-0.15,0.15,0,0.3 /",
        f"&TIME T_END={end:.17g}, DT={dt:.17g}, WALL_INCREMENT=1 /",
        f"&MISC SOLID_PHASE_ONLY=T, TMPA={(p['tf'] if p['front'] == 'radiative' else ti) - 273.15} /",
        f"&RADI RADIATION={('T' if p['front'] == 'radiative' else 'F')} /",
        f"&DUMP DT_DEVC={dump:.17g}, DT_PROF={(dump if transient else end):.17g} /",
    ]
    reaction = case.startswith(("vrx02", "vrx05", "vcpl03", "vcpl04"))
    chain = case == "vrx01"
    layered = case.startswith("vhc04")
    ns = 3 if chain else 1 if case == "vcpl04" else 2 if reaction or layered else 1
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
        if chain and i < ns:
            # A → B → C, first order in each, E = 0, unit yield.
            props += f", A={p['rates'][i - 1]:.17g}, E=0, N_S=1, NU_MATL=1, MATL_ID='{chr(65 + i)}', HEAT_OF_REACTION=0"
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
    if p["front"] == "flux":
        # An absorbed flux q with a convective loss h(T − T∞) is convection
        # to T∞ + q/h. A pure net flux (h = 0) has no FDS input on a
        # conducting material: NET_HEAT_FLUX is incompatible with MATL_ID.
        if not p.get("h"):
            raise ValueError(f"{case}: FDS cannot apply a pure net flux to a conducting material")
        surf += f", TMP_GAS_FRONT={p['tf'] + p['q'] / p['h'] - 273.15}"
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
        f"&PROF XYZ=0,0,0, IOR=3, QUANTITY='TEMPERATURE', CELL_CENTERED=T, FORMAT_INDEX={1 if transient else 2}, ID='temperature' /",
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
    if reaction or chain:
        for name in ("A",) if case == "vcpl04" else ("A", "B", "C") if chain else ("A", "B"):
            lines.append(
                f"&DEVC ID='{name}', XYZ=0,0,0, IOR=3, QUANTITY='SOLID DENSITY', MATL_ID='{name}', DEPTH={L / 2}, TIME_AVERAGED=F /"
            )
    if case.startswith("vcpl03"):
        lines.append(
            "&DEVC ID='thickness', XYZ=0,0,0, IOR=3, QUANTITY='WALL THICKNESS', TIME_AVERAGED=F /"
        )
    lines.append("&TAIL /")
    return "\n".join(lines) + "\n"


def thermakin_dt(p, level):
    return p.get("thermakin_dt", p["dt"]) / 2 ** refinement(level)[1]


def thermakin(p, level):
    """Components (.cmp) and conditions (.cnd) inputs, keyed by file suffix.

    The first listed layer is the top (front) face. Elements are cell
    centered, so N elements match N Julia cells of the same width.
    """
    case = p["case"]
    space, _ = refinement(level)
    n = (p["nodes"] - 1) * 2**space
    dt = thermakin_dt(p, level)
    L, ti = p["L"], p["ti"]
    dx = L / n
    num = lambda x: f"{x:.17g}"
    exposed = p["front"] in ("radiative", "flux")
    absorption = TK_ABSORPTION if exposed else 0.0
    trn = p.get("diffusivity", 1e-06)

    def component(name, state="S", rho=p["rho"], cp=p["cp"], k=p["k"], eps=0.0):
        law = lambda v: " ".join(num(x) for x in (tuple(v) + (0, 0, 0))[:4])
        return (
            f"COMPONENT: {name}\nSTATE: {state}\nDENSITY: {law((rho,))}\n"
            f"HEAT CAPACITY: {law((cp,))}\nCONDUCTIVITY: {law(k if isinstance(k, tuple) else (k,))}\n"
            f"TRANSPORT: {law((trn,))}\n"
            f"EMISSIVITY & ABSORPTION: {num(eps)} {num(absorption if state == 'S' else 0.0)}\n\n"
        )

    def reaction(r1, r2, p1, p2, stoich, a, e=0.0, heat=0.0, limit="L 0"):
        return (
            f"REACTION: {r1} + {r2} -> {p1} + {p2}\n"
            f"STOICHIOMETRY: {' '.join(num(s) for s in stoich)}\n"
            f"ARRHENIUS: {num(a)} {num(e)}\nHEAT: {num(heat)} 0 0 0\nTEMP LIMIT: {limit}\n\n"
        )

    comps, rxns = [], []
    if case == "vhc05":
        comps.append(component("A", k=(0.2, 0.002)))
    elif case == "vgt01":
        comps += [component("S"), component("G", "G", rho=1.0, k=0.03)]
    elif case.startswith("vgt07"):
        # With gas density φρ_g, ξ/ρ is the tracer mass fraction Y and the
        # linear boundary law a·ρ·(ξ/ρ - b) is the film flux h_m (Y - Y∞).
        # The inert carrier only fixes the pressure, which ThermaKin lacks.
        comps += [component("MATRIX"),
                  component("TRACER", "G", rho=p["pore_gas_density"], cp=1000.0, k=0.03)]
    else:
        comps.append(component("A"))
    if exposed:
        # Same material; only the exposed element emits and absorbs.
        eps = p["eps"] if p["front"] == "radiative" else TK_FLUX_EMISSIVITY
        comps.append(component("A_surface", eps=eps))
    if case.startswith("vhc04"):
        comps.append(component("B", k=1.0))
    elif case == "vrx01":
        comps += [component("B"), component("C")]
        rxns += [
            reaction("A", "NOCOMP", "B", "NOCOMP", (1, 0, 1, 0), p["rates"][0]),
            reaction("B", "NOCOMP", "C", "NOCOMP", (1, 0, 1, 0), p["rates"][1]),
        ]
    elif case.startswith("vrx02"):
        # A + A at half stoichiometry each: dξ_A/dt = -A ξ_A².
        if p["n"] != 2.0:
            raise ValueError(f"{case}: ThermaKin rates are first order in each reactant")
        comps.append(component("B"))
        rxns.append(reaction("A", "A", "B", "NOCOMP", (0.5, 0.5, 1, 0), p["a"] * p["rho"] ** (1 - p["n"])))
    elif case.startswith("vrx05"):
        y = p["solid_yield"]
        comps.append(component("B", cp=1900.0))
        if y < 1:
            # With solid swelling 1, gas volume counts at the gas density;
            # the solid density keeps the element volume fixed, as for gas
            # retained in the pores.
            comps.append(component("G", "G", rho=p["rho"], cp=1000.0, k=0.03))
        # ThermaKin's positive heat is exothermic.
        rxns.append(reaction("A", "NOCOMP", "B", "G" if y < 1 else "NOCOMP",
                             (1, 0, y, 1 - y), p["a"], heat=-p["heat"]))
    elif case.startswith("vcpl03"):
        comps.append(component("B", rho=p["product_rho"]))
        rxns.append(reaction("A", "NOCOMP", "B", "NOCOMP", (1, 0, 1, 0), p["a"]))
    elif case == "vcpl04":
        # Self-regenerating A keeps ξ_A = ρ, so A·ξ_A is the zero-order rate.
        rxns.append(reaction("A", "NOCOMP", "A", "NOCOMP", (1, 0, 1, 0), p["a"],
                             p["activation"] * TK_R / 8.314462618, -p["heat"], f"U {num(p['tf'] + 100)}"))
    swelling = 1 if case.startswith("vrx05") and p["solid_yield"] < 1 else 0
    mixtures = (f"MIXTURES\nS SWELLING: {swelling}\nL SWELLING: 0\nG SWELLING LIMIT: 1e-30\n"
                # Series conduction is exact for layers normal to the heat flow.
                "PARALL CONDUCTIVITY: 0\nPARALL TRANSPORT: 0\n\n")
    cmp_text = "".join(comps) + mixtures + "".join(rxns)

    # Layers, top first: (thickness, temperature, {component: mass fraction}).
    if case.startswith("vhc04"):
        layers = [[L - p["interface"], ti, {"A": 1}], [p["interface"], ti, {"B": 1}]]
    elif case == "vgt01":
        layers = []
        for i in range(n):
            z = L - (i + 0.5) * dx
            xi = p["peak"] * math.exp(-((z - p["z0"]) ** 2) / (2 * p["sigma0"] ** 2))
            layers.append([dx, ti, {"S": 1000, "G": float(f"{xi:.9g}")}])
    elif case.startswith("vgt07"):
        layers = [[L, ti, {"MATRIX": 1}]]
    else:
        layers = [[L, ti, {"A": 1}]]
    top_fixed = p["front"] == "fixed" and p["tf"] != ti
    if exposed or top_fixed:
        first = layers[0]
        fractions = {"A_surface": 1} if exposed else first[2]
        layers[0:1] = [[dx, p["tf"] if top_fixed else ti, fractions], [first[0] - dx, ti, first[2]]]
    if p["back"] == "fixed" and p["tb"] != ti:
        last = layers[-1]
        layers[-1:] = [[last[0] - dx, ti, last[2]], [dx, p["tb"], last[2]]]
    structure = "".join(
        f"THICKNESS: {num(t)}\nTEMPERATURE: {num(temp)}\nMASS FRACTIONS:\n"
        + "".join(f"{c} {num(v)}\n" for c, v in fractions.items()) + "\n"
        for t, temp, fractions in layers if t > 0
    )

    def boundary(name, face):
        mode = p[face]
        temp = p["tf" if face == "front" else "tb"]
        h = TK_FIXED_H if mode == "fixed" else p.get("h", 0.0) if face == "front" else 0.0
        # ThermaKin assembles gas transport only when a boundary lists a mass
        # flux, so the sealed V-GT-01 slab carries a zero-coefficient entry.
        mass = "NO\n"
        if face == "front" and case == "vgt01":
            mass = "YES\nG LIN 0 0\n"
        elif face == "front" and case.startswith("vgt07"):
            mass = (f"YES\nTRACER LIN {num(p['hm'] / p['pore_gas_density'])} "
                    f"{num(p['ambient_fraction'])}\n")
        text = (f"{name} BOUNDARY\n\nMASS TRANSPORT: {mass}\n"
                f"OUTSIDE INIT TEMP: {num(temp if mode != 'insulated' else ti)}\n"
                f"OUTSIDE HEAT RATE: 0 0 0 0\nCONVECTION COEFF: {num(h)}\n\n")
        if face == "front" and exposed:
            flux = (p["q"] / p["eps"] + TK_SIGMA * temp**4 if mode == "radiative"
                    else p["q"] / TK_FLUX_EMISSIVITY)
            text += (f"EXTERNAL RADIATION: YES\nTIME PROG1: {num(flux)} 0 1e30\n"
                     f"TIME PROG2: {num(flux)} 0 0\nREPEAT: NO\nABSORPTION MODE: MAX\n\n")
        else:
            text += "EXTERNAL RADIATION: NO\n\n"
        return text + "FLAME: NO\n\n"

    steps = p["dump"] / dt
    if abs(steps - round(steps)) > 1e-9 * steps:
        raise ValueError(f"{case}: output interval is not a whole number of steps")
    cnd_text = (
        "OBJECT TYPE: 1D\n\nOBJECT STRUCTURE\n****************\n\n" + structure
        + "OBJECT BOUNDARIES\n*****************\n\n"
        + boundary("TOP", "front") + boundary("BOTTOM", "back")
        + "INTEGRATION PARAMETERS\n**********************\n\n"
        f"ELEMENT SIZE: {num(dx)}\nTIME STEP: {num(dt)}\n"
        # Half a step short: the loop runs while t < DURATION with t = k·dt.
        f"DURATION: {num(p['end'] - dt / 2)}\n\n"
        f"OUTPUT FREQUENCY:\nELEMENTS: 1\nTIME STEPS: {round(steps)}\n"
    )
    return {".cmp": cmp_text, ".cnd": cnd_text}


def thermakin_snapshots(path):
    """Parse every reported state: time, element centers (depth), T, ξ and
    boundary mass outflows (kg/s per unit area)."""
    text = path.read_text()
    if "Calculations are complete." not in text or "Error:" in text:
        raise RuntimeError(f"{path}: unsuccessful ThermaKin termination")
    snapshots = []
    for block in text.split("Time [s] =")[1:]:
        lines = block.splitlines()
        names, rows, thickness, outflow = None, [], None, {}
        for line in lines[1:]:
            fields = [f for f in line.split("\t") if f.strip()]
            if line.startswith("BOUNDARY"):
                flow_names = fields[4:]
            elif line.startswith(("TOP", "BOTTOM")):
                values = [float(f) for f in fields[1:]]
                outflow[fields[0].strip()] = dict(zip(flow_names, values[2:]))
            elif line.startswith("FROM TOP"):
                names = fields[3:]
            elif line.startswith("Total thickness"):
                thickness = float(fields[-1])
            elif names is not None and fields and line[0] in "0123456789":
                rows.append([float(f) for f in fields])
        if names is None or not rows or thickness is None:
            raise RuntimeError(f"{path}: incomplete output block")
        snapshots.append(dict(
            t=float(lines[0]),
            depth=[r[0] for r in rows],
            T=[r[1] for r in rows],
            xi={c: [r[2 + i] for r in rows] for i, c in enumerate(names)},
            thickness=thickness,
            outflow=outflow,
        ))
    return snapshots


def run(job, executable, workspace, extract_only):
    solver, case, level = job
    if level not in COARSE:
        return _run(job, executable, workspace, extract_only)
    try:
        return _run(job, executable, workspace, extract_only)
    except (ValueError, RuntimeError, subprocess.CalledProcessError,
            subprocess.TimeoutExpired) as err:
        print(f"{solver} {case} {level}: unavailable: {err}", flush=True)
        return dict(case=case, level=level, unavailable=str(err)[:300]), None


def _run(job, executable, workspace, extract_only):
    solver, case, level = job
    p = parameters(case)
    path = workspace / solver / level / case
    path.mkdir(parents=True, exist_ok=True)
    if solver == "thermakin":
        decks = {path / (case + suffix): text for suffix, text in thermakin(p, level).items()}
    else:
        deck = path / (case + (".fds" if solver == "fds" else ".data"))
        decks = {deck: (fds if solver == "fds" else gpyro)(p, level)}
    if extract_only:
        for deck, expected in decks.items():
            if deck.read_text() != expected:
                raise RuntimeError(f"{deck}: input differs from generated deck; rerun")
    elif solver == "thermakin":
        for deck, expected in decks.items():
            deck.write_text(expected)
        out = path / f"{case}.out"
        out.unlink(missing_ok=True)  # ThermaKin refuses to overwrite its output
        with (path / "run.log").open("w") as log:
            subprocess.run(
                [str(executable)],
                input=f"{case}.cmp\n{case}.cnd\n{case}.out\n",
                cwd=path,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=3600,
                check=True,
            )
    else:
        deck.write_text(decks[deck])
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
    if solver == "thermakin":
        # Completion is checked in the output itself by thermakin_snapshots.
        success = "Error:" not in log and "Done" in log
    else:
        success = (
            "FDS completed successfully"
            if solver == "fds"
            else "End of simulation reached successfully."
        ) in log
    if not success:
        raise RuntimeError(f"{job}: unsuccessful solver termination: {log[-500:]}")
    banner = (path / f"{case}.out").read_text() if solver != "gpyro" else log
    if VERSIONS[solver] not in banner:
        raise RuntimeError(f"{job}: expected pinned solver version {VERSIONS[solver]}")
    return extract(solver, case, level, path, list(decks))


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
    elif solver == "thermakin":
        # Printed centers (7 significant digits) must equal the requested
        # uniform layout; faces are depths from the top (front) face.
        cnd = (path / f"{case}.cnd").read_text()
        h = float(re.search(r"ELEMENT SIZE: (\S+)", cnd)[1])
        first = thermakin_snapshots(path / f"{case}.out")[0]
        n = len(first["depth"])
        if any(abs(c - (i + 0.5) * h) > 1e-6 * n * h for i, c in enumerate(first["depth"])) or (
                not math.isclose(first["thickness"], n * h, rel_tol=1e-6)):
            raise RuntimeError(f"{case}: printed element centers differ from the requested layout")
        coordinates = [i * h for i in range(n + 1)]
        layers = []
        convention = ("uniform element-centered unknowns; faces from ELEMENT SIZE, verified "
                      "against printed centers; boundary conditions act on boundary-element centers")
        unknowns = n
    elif case in GPYRO_LUMPED:
        # ZEROD prints no profile: one lumped unknown spans the declared depth.
        coordinates = [0.0, parameters(case)["L"]]
        layers = []
        convention = "ZEROD homogeneous kinetics; one lumped unknown"
        unknowns = 1
    else:
        profile = sorted(path.glob(f"{case}_profile_01_01_*.csv"))[0]
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
    family = re.match(r"v([a-z]+)", case)[1]
    rest = case[1 + len(family):]
    # V-CPL-03's P1/P2/P3 decks are subcases of one case.
    return f"V-{family.upper()}-{rest[:2] if family == 'cpl' else rest}"


def extract(solver, case, level, path, decks):
    p = parameters(case)
    label = LABELS[solver]
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

    if solver == "thermakin":
        native = path / f"{case}.out"
        thermakin_samples(p, thermakin_snapshots(native), add)
        return retain(solver, case, level, path, decks + [native], samples, metadata)
    space, _ = refinement(level)
    native = path / (
        f"{case}_devc.csv" if solver == "fds" else f"{case}_summary_01_0001.csv"
    )
    data = numeric_csv(native, 2 if solver == "fds" else 1)
    tolerance = max(p["dt"], 0.0001) * 1.1
    if abs(data[-1][0] - p["end"]) > tolerance:
        raise RuntimeError(f"{native}: final time {data[-1][0]} != {p['end']}")
    if any((b[0] <= a[0] for a, b in zip(data, data[1:]))):
        raise RuntimeError(f"{native}: nonmonotone output times")
    files = decks + [native]
    profile = None
    transient = solver == "fds" and fds_transient_profiles(case)
    if case.startswith(("vhc", "vcpl01", "vcpl04", "vgt")) or case == "vrx04a":
        quantity = (
            "YJ(01)" if case.startswith("vgt07")
            else "PRESSURE" if case.startswith("vgt08")
            else "TEMPERATURE"
        )
        profile = path / (
            f"{case}_prof_1.csv"
            if solver == "fds"
            else f"{case}_profile_01_01_{quantity}.csv"
        )
        files.append(profile)
        prof = numeric_csv(profile, (3 if transient else 2) if solver == "fds" else 0)

    def profile_at(t):
        """(depth from the front face, native value) at output time t."""
        rows = prof if solver == "fds" else prof[1:]
        row = min(rows, key=lambda r: abs(r[0] - t))
        if abs(row[0] - t) > tolerance:
            raise RuntimeError(f"{profile}: missing t={t}")
        if solver == "fds":
            n = int(row[1])
            return list(zip(row[2:2 + n], row[2 + n:2 + 2 * n]))
        return list(zip(prof[0][1:], row[1:]))

    def node_mean(depth, values):
        """Volume mean of boundary-node values: trapezoids over half end cells."""
        return sum(
            (b - a) * (u + v) / 2
            for a, b, u, v in zip(depth, depth[1:], values, values[1:])
        ) / p["L"]

    def isothermal(column, temperature):
        if max(abs(r[column] + 273.15 - temperature) for r in data) > 1e-6:
            raise RuntimeError(f"{native}: temperature left {temperature} K")

    if case == "vhc01":
        for t in (30.0, 60.0, 90.0, 180.0):
            add(
                "T_profiles",
                f"t = {t} s",
                sorted((d, v + 273.15) for d, v in profile_at(t)),
                f", t = {t} s",
            )
        # The exact surface form is 0/0 at t = 0.
        add(
            "surface_T_history",
            "surface Newton solve",
            ((r[0], r[1] + 273.15) for r in data if r[0] > 0),
        )
    elif case == "vhc03":
        alpha = p["k"] / (p["rho"] * p["cp"])
        for fo in (0.1, 0.5, 1.5):
            add(
                "T_profiles",
                f"Fo = {fo}",
                sorted(
                    (1 - d / p["L"], v + 273.15)
                    for d, v in profile_at(fo * p["L"] ** 2 / alpha)
                ),
                f", Fo = {fo}",
            )
        depth = prof[0][1:]
        add(
            "mean_T_history",
            "volume mean",
            ((r[0], node_mean(depth, r[1:]) + 273.15) for r in prof[1:]),
        )
    elif case == "vcpl01":
        depth = prof[0][1:]
        for d in (0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3):
            i = min(range(len(depth)), key=lambda i: abs(depth[i] - d))
            if abs(depth[i] - d) > 1e-9:
                raise RuntimeError(f"{profile}: no node at depth {d}")
            add(
                "T_histories_at_depths",
                f"x = {d * 1e3} mm",
                ((r[0], r[1 + i] + 273.15) for r in prof[1:] if r[0] > 0),
                f", x = {d * 1e3} mm",
            )

        def front(values):
            # First crossing of the melting isotherm, linear between nodes.
            T = [v + 273.15 - p["melt_T"] for v in values]
            for a, b, u, v in zip(depth, depth[1:], T, T[1:]):
                if u * v <= 0 and u != v:
                    return a + (b - a) * u / (u - v)
            return None

        fronts = [(r[0], front(r[1:])) for r in prof[1:] if r[0] > 0]
        add(
            "front_trajectory",
            "T_m isotherm",
            ((t, s) for t, s in fronts if s is not None),
        )
    elif case == "vrx01":
        # Gpyro reports mass fractions of the initial density.
        columns, scale = ((4, 5, 6), 1.0) if solver == "fds" else ((1, 2, 3), p["rho"])
        for name, column in zip("ABC", columns):
            add(
                "xi_histories",
                f"ξ_{name}",
                ((r[0], scale * r[column]) for r in data),
                f" ξ_{name}",
            )
        if solver == "gpyro":
            isothermal(5, p["ti"])
            if max(abs(r[4] - 1) for r in data) > 1e-6:
                raise RuntimeError(f"{native}: condensed mass not conserved")
    elif case.startswith("vgt07"):
        for t in (50.0, 150.0, 400.0):
            add("profiles", f"t = {t} s", profile_at(t), f", t = {t} s")
        depth = prof[0][1:]
        # Gpyro's native gas density (RHOGOFT constant R = 8314/M), not the
        # benchmark gas constant.
        density = p["porosity"] * 1e5 * 32 / (8314 * 300)
        add(
            "uptake",
            "surface influx",
            ((r[0], p["hm"] * (p["ambient_fraction"] - r[1])) for r in prof[1:] if r[0] > 0),
        )
        add(
            "inventory",
            "stored tracer",
            ((r[0], density * node_mean(depth, r[1:]) * p["L"]) for r in prof[1:]),
        )
        diffusivity = path / f"{case}_profile_01_02_D12.csv"
        files.append(diffusivity)
        if any(abs(v / p["diffusivity"] - 1) > 1e-8
               for row in numeric_csv(diffusivity, 0)[2:] for v in row[1:]):
            raise RuntimeError(f"{diffusivity}: binary diffusivity is not {p['diffusivity']}")
        isothermal(1, 300.0)
    elif case.startswith("vgt08"):
        add(
            "pressure",
            "steady pressure",
            sorted((p["L"] - d, v + 1e5) for d, v in profile_at(p["end"])),
        )
        # Outward-positive kg/(m² s): the top-face flux is printed in g, and
        # MLR integrates both outward boundary fluxes.
        top = [(r[0], -r[3] / 1000) for r in data if r[0] > 0]
        bottom = [(t, r[2] / 1000 - f) for (t, f), r in zip(top, (r for r in data if r[0] > 0))]
        inlet, outlet = (("bottom", bottom), ("top", top)) if case[-1] == "a" else (
            ("top", top), ("bottom", bottom))
        add("boundary_flux", f"inlet ({inlet[0]})", inlet[1], " inlet")
        add("boundary_flux", f"outlet ({outlet[0]})", outlet[1], " outlet")
        add(
            "boundary_flux",
            "net (inlet + outlet)",
            ((t, a + b) for (t, a), (_, b) in zip(inlet[1], outlet[1])),
            " net",
        )
        isothermal(1, 300.0)
    elif case.startswith("vhc02"):
        for fo in (0.05, 0.1, 0.2, 0.5):
            add(
                "T_profiles",
                f"Fo = {fo}",
                sorted(
                    ((p["L"] - d) / p["L"], v + 273.15)
                    for d, v in profile_at(fo * 0.1**2 / (0.5 / 1000000.0))
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
            means = [(r[0], node_mean(depth, r[1:]) + 273.15) for r in prof[1:]]
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
        # FDS releases its gas; Gpyro retains it (CJ), with heat capacity 1000.
        T, A, B, G = (3, 4, 5, None) if solver == "fds" else (1, 2, 3, 4)
        species = [("A", A), ("B", B)] + ([("G", G)] if p["solid_yield"] < 1 else [])
        add("temperature", "temperature", ((r[0], r[T] + 273.15) for r in data))
        for name, col in species:
            add("species", name, ((r[0], r[col]) for r in data), " " + name)
        add(
            "enthalpy",
            "enthalpy",
            (
                (
                    r[0],
                    p["L"]
                    * (1200 * r[A] + 1900 * r[B] + (0.0 if G is None else 1000 * r[G]))
                    * (r[T] + 273.15 - 298.15),
                )
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
    return retain(solver, case, level, path, files, samples, metadata)


def thermakin_samples(p, snapshots, add):
    """Observables from element-center values (ThermaKin has no face values)."""
    case, L = p["case"], p["L"]
    times = [s["t"] for s in snapshots]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise RuntimeError(f"{case}: nonmonotone output times")
    if not math.isclose(times[-1], p["end"], rel_tol=1e-6):
        raise RuntimeError(f"{case}: final time {times[-1]} != {p['end']}")

    def at(t):
        s = min(snapshots, key=lambda s: abs(s["t"] - t))
        if not math.isclose(s["t"], t, rel_tol=1e-6):
            raise RuntimeError(f"{case}: missing output at t={t}")
        return s

    def widths(s):
        # Centers print with 7 significant digits, too coarse to difference
        # into widths; verify the uniform layout and use thickness/N.
        n = len(s["depth"])
        h = s["thickness"] / n
        if any(abs(c - (i + 0.5) * h) > 1e-6 * s["thickness"] for i, c in enumerate(s["depth"])):
            raise RuntimeError(f"{case}: nonuniform elements at t={s['t']}")
        return [h] * n

    def history(curve, reference, value, suffix="", initial=True):
        add(curve, reference,
            ((s["t"], value(s)) for s in snapshots if initial or s["t"] > 0), suffix)

    surface = lambda s: s["T"][0]
    back = lambda s: s["T"][-1]
    if case == "vhc01":
        for t in (30.0, 60.0, 90.0, 180.0):
            s = at(t)
            add("T_profiles", f"t = {t} s", zip(s["depth"], s["T"]), f", t = {t} s")
        history("surface_T_history", "surface Newton solve", surface)
    elif case.startswith("vhc02"):
        for fo in (0.05, 0.1, 0.2, 0.5):
            s = at(fo * L**2 / (p["k"] / (p["rho"] * p["cp"])))
            add("T_profiles", f"Fo = {fo}",
                sorted(((L - d) / L, T) for d, T in zip(s["depth"], s["T"])), f", Fo = {fo}")
        # The eigenfunction series is singular at the initial corner.
        history("surface_T_history", "front face", surface, initial=False)
        history("back_T_history", "adiabatic back face", back, initial=False)
    elif case == "vhc03":
        alpha = p["k"] / (p["rho"] * p["cp"])
        for fo in (0.1, 0.5, 1.5):
            s = at(fo * L**2 / alpha)
            add("T_profiles", f"Fo = {fo}",
                sorted((1 - d / L, T) for d, T in zip(s["depth"], s["T"])), f", Fo = {fo}")
        history("mean_T_history", "volume mean",
                lambda s: sum(w * T for w, T in zip(widths(s), s["T"])) / L)
    elif case.startswith(("vhc04", "vhc05", "vcpl04")):
        s = snapshots[-1]
        add(
            "steady_profiles" if case == "vcpl04" else "T_profile_steady",
            "δ = 0.44 ≈ 0.50 δ_c" if case == "vcpl04" else "steady solver profile",
            sorted((L - d, T) for d, T in zip(s["depth"], s["T"])),
            ", δ = 0.44" if case == "vcpl04" else "",
        )
    elif case == "vrx01":
        for name in ("A", "B", "C"):
            history("xi_histories", f"ξ_{name}", lambda s: s["xi"][name][0], f" ξ_{name}")
    elif case.startswith("vrx02"):
        history("xi_histories", f"ξ_A (n = {p['n']})", lambda s: s["xi"]["A"][0])
    elif case == "vrx04a":
        history("T_history_vs_exact", "slab mean temperature",
                lambda s: sum(w * T for w, T in zip(widths(s), s["T"])) / L)
    elif case == "vrx04b":
        history("surface_T_history", "surface temperature", surface)
    elif case.startswith("vrx05"):
        history("temperature", "temperature", lambda s: s["T"][0])
        names = ("A", "B", "G") if p["solid_yield"] < 1 else ("A", "B")
        for name in names:
            history("species", name, lambda s: s["xi"][name][0], " " + name)
        heat_capacity = dict(A=1200.0, B=1900.0, G=1000.0)
        history("enthalpy", "enthalpy", lambda s: L * sum(
            heat_capacity[c] * s["xi"][c][0] for c in names) * (s["T"][0] - 298.15))
    elif case.startswith("vcpl03"):
        perm = "P" + case[-1]
        history("thickness_histories", perm, lambda s: s["thickness"] / L, " " + perm)
    elif case == "vgt01":
        for t in (2.0, 4.0, 8.0):
            s = at(t)
            add("xi_profiles", f"t = {t} s",
                sorted((L - d, xi) for d, xi in zip(s["depth"], s["xi"]["G"])), f", t = {t} s")

        def variance(s):
            w, z, xi = widths(s), [L - d for d in s["depth"]], s["xi"]["G"]
            m = sum(a * b for a, b in zip(xi, w))
            mean = sum(a * b * c for a, b, c in zip(xi, w, z)) / m
            return sum(a * b * (c - mean) ** 2 for a, b, c in zip(xi, w, z)) / m

        history("variance_history", "discrete second moment", variance)
    elif case.startswith("vgt07"):
        for t in (50.0, 150.0, 400.0):
            s = at(t)
            add("profiles", f"t = {t} s",
                zip(s["depth"], (xi / p["pore_gas_density"] for xi in s["xi"]["TRACER"])),
                f", t = {t} s")
        # Unit area: the printed boundary flow is the mass flux.
        history("uptake", "surface influx",
                lambda s: -s["outflow"]["TOP"]["TRACER"], initial=False)
        history("inventory", "stored tracer",
                lambda s: sum(w * xi for w, xi in zip(widths(s), s["xi"]["TRACER"])))
    else:
        raise ValueError(case)


def retain(solver, case, level, path, files, samples, metadata):
    """Write extracted samples; copy and hash native inputs/outputs and grid."""
    p = parameters(case)
    space, time = refinement(level)
    refined = level
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
    lumped = solver == "gpyro" and case in GPYRO_LUMPED
    expected_intervals = 1 if lumped else (p["nodes"] - 1) * 2**space
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
                case.startswith(("vrx02", "vrx05")) and
                math.isclose(grid["coordinates_m"][-1], .75*p["L"], rel_tol=1e-8)):
            grid["coordinate_exception"] = "native two-node position formula: printed span 3L/4; control-volume sum L; homogeneous ODE only"
        else:
            raise RuntimeError(f"{solver} {case}: native grid does not span the physical domain")
    grid_file = raw / "grid.json"
    grid_file.write_text(json.dumps(grid, indent=2) + "\n")
    hashes = {}
    hashes[str(grid_file.relative_to(DEST))] = sha256(grid_file)
    for f in files:
        if solver == "thermakin" and f.suffix == ".out":
            # Full element profiles at every output. The start date and total
            # runtime lines are timing records, not published; without them
            # and with an mtime-free gzip the retained file is reproducible.
            dest = raw / (f.name + ".gz")
            text = "".join(line for line in f.read_text().splitlines(keepends=True)
                           if not line.startswith(("Started on:", "Total runtime")))
            with dest.open("wb") as stream, gzip.GzipFile(
                filename="", mode="wb", fileobj=stream, mtime=0
            ) as archive:
                archive.write(text.encode())
        else:
            dest = raw / f.name
            shutil.copyfile(f, dest)
        hashes[str(dest.relative_to(DEST))] = sha256(dest)
    hashes[str(out.relative_to(DEST))] = sha256(out)
    print(f"{solver} {case} {refined}: verified {len(samples)} samples", flush=True)
    return (
        dict(
            case=case,
            level=refined,
            dt=thermakin_dt(p, level) if solver == "thermakin" else p["dt"] / 2**time,
            requested_intervals=expected_intervals,
            grid=grid,
            observed_step_history="not retained; output intervals are not integration steps",
            timestep_interpretation=("outer step; native solid substepping permitted"
                                     if solver == "fds" else
                                     "fixed step; Crank-Nicolson with one linearization per step"
                                     if solver == "thermakin" else
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
    unavailable = [audit for audit, series in results if series is None]
    results = [r for r in results if r[1] is not None]
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
        unavailable=unavailable,
    )
    if args.solver == "thermakin":
        # build_thermakin.sh records the compiler and the source hash.
        build = dict(line.split(": ", 1) for line in
                     (executable.parent / "build_info.txt").read_text().splitlines())
        if build["executable_sha256"] != provenance["executable_sha256"]:
            raise RuntimeError("build_info.txt does not describe this executable")
        provenance.update(
            source_sha256=build["source_sha256"],
            compiler=build["compiler"],
            compiler_flags=build["flags"],
            source_modifications=("start-time banner only: POSIX localtime_r/asctime_r "
                                  "replace MSVC localtime_s/asctime_s"),
        )
    (DEST / f"{args.solver}_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
