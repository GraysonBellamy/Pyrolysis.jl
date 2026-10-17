#!/usr/bin/env python3
"""Newton/integrator baseline from the performance-audit CSVs (no new runs).

Reads docs/performance_audit/{solves,accuracy,burnout}.csv and writes
scripts/mesh_phase0/results/newton_baseline.md. Standard library only.

    python3 scripts/mesh_phase0/baselines/newton_baseline.py
"""
import csv
import os
import sys
from collections import OrderedDict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
AUDIT = os.path.join(ROOT, "docs", "performance_audit")
OUT = os.path.join(ROOT, "scripts", "mesh_phase0", "results", "newton_baseline.md")

CASE_ORDER = ["eulerian", "ale", "ale_darcy"]
INTEGRATOR_ORDER = ["KenCarp4", "TRBDF2", "Rodas5P", "FBDF"]
STRATEGY_LABEL = {"direct": "exact", "approx": "reduced"}


def read(name):
    with open(os.path.join(AUDIT, name), newline="") as f:
        return list(csv.DictReader(f))


def f(x):
    return float(x)


def fmt(x, nd=3):
    if isinstance(x, str):
        return x
    if x != x:  # NaN
        return "nan"
    if abs(x) >= 1e4 or (abs(x) < 1e-2 and x != 0):
        return "%.2e" % x
    return ("%%.%df" % nd) % x


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def deck_key(r):
    return (CASE_ORDER.index(r["case"]), int(r["cells"]))


def main():
    solves = read("solves.csv")
    accuracy = read("accuracy.csv")
    burnout = read("burnout.csv")
    acc = {
        (r["case"], r["cells"], r["strategy"], r["linear_solver"], r["integrator"]): r
        for r in accuracy
    }

    # ---- (a) exact vs reduced Jacobian, KenCarp4 ----------------------------
    kc = [r for r in solves if r["integrator"] == "KenCarp4"]
    groups = OrderedDict()
    for r in sorted(kc, key=lambda r: (deck_key(r), r["linear_solver"])):
        groups.setdefault((r["case"], r["cells"], r["linear_solver"]), {})[r["strategy"]] = r
    rows_a = []
    stats_a = []  # (case, cells, ls, exact row, reduced row)
    for (case, cells, ls), d in groups.items():
        ex, rd = d.get("direct"), d.get("approx")
        if ex is None or rd is None:
            continue
        stats_a.append((case, cells, ls, ex, rd))
        rows_a.append(
            [
                case,
                cells,
                ls,
                fmt(f(ex["nonliniter"]) / f(ex["steps"]), 2),
                ex["nonlinfail"],
                ex["steps"],
                fmt(f(ex["median_s"]), 3),
                fmt(f(rd["nonliniter"]) / f(rd["steps"]), 2),
                rd["nonlinfail"],
                rd["steps"],
                fmt(f(rd["median_s"]), 3),
            ]
        )
    header_a = [
        "case",
        "cells",
        "linsolve",
        "exact its/step",
        "exact nlfail",
        "exact steps",
        "exact median s",
        "reduced its/step",
        "reduced nlfail",
        "reduced steps",
        "reduced median s",
    ]

    # ---- (b) integrators on the same decks (reduced Jacobian, sparspak) -----
    ints = [r for r in solves if r["strategy"] == "approx" and r["linear_solver"] == "sparspak"]
    ints.sort(key=lambda r: (deck_key(r), INTEGRATOR_ORDER.index(r["integrator"])))
    rows_b = []
    for r in ints:
        a = acc.get((r["case"], r["cells"], r["strategy"], r["linear_solver"], r["integrator"]))
        rows_b.append(
            [
                r["case"],
                r["cells"],
                r["integrator"],
                r["steps"],
                r["rejects"],
                r["nf"],
                r["njacs"],
                fmt(f(r["nonliniter"]) / f(r["steps"]), 1),
                r["nonlinfail"],
                fmt(f(r["median_s"]), 3),
                fmt(f(a["max_temperature_error_K"])) if a else "-",
                fmt(f(a["final_state_relative_inf"])) if a else "-",
                fmt(f(a["max_final_node_error_m"])) if a else "-",
            ]
        )
    header_b = [
        "case",
        "cells",
        "integrator",
        "steps",
        "rejects",
        "nf",
        "njacs",
        "its/step",
        "nlfail",
        "median s",
        "max T err K",
        "final state rel inf",
        "max node err m",
    ]

    # ---- (c) burnout deck ---------------------------------------------------
    rows_c = []
    for r in burnout:
        rows_c.append(
            [
                r["total_cells"],
                STRATEGY_LABEL[r["strategy"]],
                r["linear_solver"],
                r["integrator"],
                r["outputs"],
                r["retcode"],
                r["steps"],
                r["rejects"],
                r["nf"],
                r["njacs"],
                fmt(f(r["nonliniter"]) / f(r["steps"]), 1),
                r["nonlinfail"],
                fmt(f(r["median_s"]), 2),
                fmt(f(r["peak_mlr_g_m2_s"]), 2),
                fmt(f(r["peak_time_s"]), 1),
                fmt(f(r["burnout_time_s"]), 1),
                fmt(f(r["mass_loss_percent"]), 3),
                fmt(f(r["mass_ledger_rel"])),
            ]
        )
    header_c = [
        "cells",
        "Jacobian",
        "linsolve",
        "integrator",
        "outputs",
        "retcode",
        "steps",
        "rejects",
        "nf",
        "njacs",
        "its/step",
        "nlfail",
        "median s",
        "peak MLR g/m2/s",
        "t_peak s",
        "t_burnout s",
        "mass loss %",
        "ledger rel",
    ]

    # ---- (d) factual sentences computed from the rows -----------------------
    # exact vs reduced: identical, fewer or more iterations/failures
    same = [s for s in stats_a if s[3]["nonliniter"] == s[4]["nonliniter"] and s[3]["nonlinfail"] == s[4]["nonlinfail"]]
    ale_pairs = [s for s in stats_a if s[0] != "eulerian"]
    its_lo, its_hi = (
        min(min(f(s[3]["nonliniter"]) / f(s[3]["steps"]), f(s[4]["nonliniter"]) / f(s[4]["steps"])) for s in stats_a),
        max(max(f(s[3]["nonliniter"]) / f(s[3]["steps"]), f(s[4]["nonliniter"]) / f(s[4]["steps"])) for s in stats_a),
    )
    fail_cmp = []
    for case, cells, ls, ex, rd in ale_pairs:
        fail_cmp.append("%s/%s/%s %s vs %s" % (case, cells, ls, ex["nonlinfail"], rd["nonlinfail"]))
    fewer = sum(1 for s in ale_pairs if int(s[3]["nonlinfail"]) < int(s[4]["nonlinfail"]))
    equal = sum(1 for s in ale_pairs if int(s[3]["nonlinfail"]) == int(s[4]["nonlinfail"]))
    more = len(ale_pairs) - fewer - equal

    # integrator ranks per deck
    decks = OrderedDict()
    for r in ints:
        decks.setdefault((r["case"], r["cells"]), []).append(r)
    time_rank = {k: [] for k in INTEGRATOR_ORDER}
    err_rank = {k: [] for k in INTEGRATOR_ORDER}
    step_rank = {k: [] for k in INTEGRATOR_ORDER}
    for key, rs in decks.items():
        by_t = sorted(rs, key=lambda r: f(r["median_s"]))
        for i, r in enumerate(by_t):
            time_rank[r["integrator"]].append(i + 1)
        by_e = sorted(
            rs,
            key=lambda r: f(acc[(r["case"], r["cells"], r["strategy"], r["linear_solver"], r["integrator"])]["max_temperature_error_K"]),
        )
        for i, r in enumerate(by_e):
            err_rank[r["integrator"]].append(i + 1)
        by_s = sorted(rs, key=lambda r: int(r["steps"]))
        for i, r in enumerate(by_s):
            step_rank[r["integrator"]].append(i + 1)
    n_decks = len(decks)

    def mean(xs):
        return sum(xs) / len(xs)

    def rank_sentence(name, ranks):
        order = sorted(INTEGRATOR_ORDER, key=lambda k: mean(ranks[k]))
        parts = []
        for k in order:
            wins = sum(1 for x in ranks[k] if x == 1)
            parts.append("%s (mean rank %.2f, first on %d/%d)" % (k, mean(ranks[k]), wins, n_decks))
        return "By %s the integrators rank %s." % (name, ", ".join(parts))

    sentences = [
        "On the Eulerian decks the exact and reduced Jacobians are the same matrix, so "
        "their KenCarp4 rows are identical (%d of %d exact/reduced pairs identical in Newton "
        "iterations and failures); on the ALE decks the exact Jacobian does not reduce the Newton "
        "iterations per accepted step, which stay within %.1f-%.1f for both strategies."
        % (len(same), len(stats_a), its_lo, its_hi),
        "On the %d ALE pairs the exact Jacobian has fewer nonlinear failures than the reduced one in %d, "
        "the same number in %d and more in %d (exact vs reduced: %s), at %.1f-%.1fx the reduced median time at 100 cells."
        % (
            len(ale_pairs),
            fewer,
            equal,
            more,
            "; ".join(fail_cmp),
            min(f(s[3]["median_s"]) / f(s[4]["median_s"]) for s in ale_pairs if s[1] == "100"),
            max(f(s[3]["median_s"]) / f(s[4]["median_s"]) for s in ale_pairs if s[1] == "100"),
        ),
        rank_sentence("median wall time on the six decks (reduced Jacobian, sparspak)", time_rank)
        + " "
        + rank_sentence("accepted steps", step_rank),
        rank_sentence(
            "maximum temperature error against the tight exact-Jacobian reference (abstol 1e-10, reltol 1e-8)",
            err_rank,
        )
        + " Rodas5P takes one Jacobian per step (njacs = steps) and no Newton iterations; on the burnout deck "
        "TRBDF2 needs %s steps and %s rejects against KenCarp4's %s and %s, at %.2f vs %.2f s."
        % (
            [r for r in burnout if r["integrator"] == "TRBDF2"][0]["steps"],
            [r for r in burnout if r["integrator"] == "TRBDF2"][0]["rejects"],
            [r for r in burnout if r["integrator"] == "KenCarp4" and r["linear_solver"] == "sparspak" and r["outputs"] == "1001"][0]["steps"],
            [r for r in burnout if r["integrator"] == "KenCarp4" and r["linear_solver"] == "sparspak" and r["outputs"] == "1001"][0]["rejects"],
            f([r for r in burnout if r["integrator"] == "TRBDF2"][0]["median_s"]),
            f([r for r in burnout if r["integrator"] == "KenCarp4" and r["linear_solver"] == "sparspak" and r["outputs"] == "1001"][0]["median_s"]),
        ),
    ]

    md = []
    md.append("# Newton and integrator baseline (performance-audit CSVs)\n")
    md.append(
        "Source: `docs/performance_audit/solves.csv`, `accuracy.csv`, `burnout.csv` "
        "(gitignored audit data; no new runs). Generated by "
        "`scripts/mesh_phase0/baselines/newton_baseline.py`.\n"
    )
    md.append(
        "Decks: `eulerian`, `ale`, `ale_darcy` at 25 and 100 cells, 80 s horizon, "
        "abstol 1e-8 / reltol 1e-6; `burnout` = V-BO-01 (60 PMMA + 122 Kaowool cells, 500 s). "
        "`exact` = `DirectSolve` (full Jacobian), `reduced` = `ApproxSolve`. "
        "`its/step` = `nonliniter / steps` (KenCarp4 has four implicit stages per step; "
        "`nonliniter` counts every Newton iteration of every stage). "
        "`median s` is the median warm solve time of the audit machine.\n"
    )
    md.append("## (a) Exact versus reduced Jacobian, KenCarp4\n")
    md.append(table(header_a, rows_a) + "\n")
    md.append("## (b) Integrators on the same decks (reduced Jacobian, sparspak)\n")
    md.append(
        "Accuracy columns are from `accuracy.csv`: errors against a tight exact-Jacobian "
        "reference (`DirectSolve`, LU, abstol 1e-10, reltol 1e-8). `max T err K` is the maximum "
        "temperature difference over cells and output times, `final state rel inf` the relative "
        "infinity-norm difference of the final state vector, `max node err m` the largest final "
        "node-position difference (ALE only). The audit README notes that `final_state_relative_inf` "
        "is dominated by one shared component and is not an accuracy metric on its own.\n"
    )
    md.append(table(header_b, rows_b) + "\n")
    md.append("## (c) Burnout deck V-BO-01 (182 cells, reduced Jacobian)\n")
    md.append(
        "The `outputs = 3` row saves only three times, so its peak MLR, peak time and "
        "burnout time are read off a 3-point output grid and are not comparable with the "
        "1001-output rows; its step statistics are identical to the sparspak 1001-output row.\n"
    )
    md.append(table(header_c, rows_c) + "\n")
    md.append("## (d) Findings\n")
    for i, s in enumerate(sentences, 1):
        md.append("%d. %s" % (i, s))
    md.append("")
    with open(OUT, "w") as fh:
        fh.write("\n".join(md))
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
