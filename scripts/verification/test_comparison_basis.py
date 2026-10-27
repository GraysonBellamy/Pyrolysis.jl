"""Regression checks for spatial indexing and independent refinement controls."""
import importlib.util
from pathlib import Path
import re
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


runner = module("matched_runner", ROOT / "test/verification/reference/matched/run_comparisons.py")
audit = module("basis_audit", Path(__file__).with_name("audit_comparison_basis.py"))
migration = module("v2_migration", ROOT / "test/verification/reference/gpyro/v2_migration/run_migration.py")


class ComparisonBasisTest(unittest.TestCase):
    def test_all_active_gpyro_overlays_use_v2(self):
        for path in (ROOT / "test/verification/cases").glob("*.jl"):
            text = path.read_text()
            self.assertNotIn("Gpyro 0.8200", text, str(path))
            self.assertIsNone(re.search(r'"gpyro/(?:vhc|vrad|vcpl|vrx|vbo)', text), str(path))

    def test_v2_migration_preserves_ode_dimension_and_independent_controls(self):
        for case in migration.CASES:
            decks = {level: migration.deck_text(case, level) for level in migration.LEVELS}
            nodes = {level: int(re.search(r"NCELLZ\(1\)\s*=\s*(\d+)", text)[1])
                     for level, text in decks.items()}
            dt = {level: float(re.search(r"DT0\s*=\s*([^,]+)", text)[1])
                  for level, text in decks.items()}
            self.assertEqual([nodes[k] for k in ("base", "time_half", "time_quarter")], [nodes["base"]]*3)
            self.assertEqual(dt["time_quarter"], dt["base"]/4)
            self.assertEqual(dt["space_half"], dt["time_quarter"])
            self.assertEqual(nodes["space_half"]-1, 2*(nodes["base"]-1))
            if case == "vrx01":
                self.assertEqual(set(nodes.values()), {1})
                self.assertRegex(decks["base"], r"SOLVE_GAS_YJ\s*=\s*\.FALSE\.")
                self.assertRegex(decks["base"], r"IO2\s*=\s*1,")

    def test_both_fds_layers_receive_spacing(self):
        text = runner.fds(runner.parameters("vhc04a"), "base")
        sizes = {int(layer): float(h) for layer, h in
                 re.findall(r"CELL_SIZE\((\d+)\)=([\d.e+-]+)", text)}
        self.assertEqual(sizes, {1: .02/80, 2: .02/80})

    def test_independent_controls(self):
        p = runner.parameters("vhc02a")
        decks = {level: runner.gpyro(p, level) for level in runner.LEVELS}
        values = lambda pattern: {k: float(re.search(pattern, v)[1]) for k, v in decks.items()}
        dt = values(r"DT0=([^,]+)")
        nodes = values(r"NCELLZ\(1\)=(\d+)")
        self.assertEqual([nodes[k] for k in ("base", "time_half", "time_quarter")], [81]*3)
        self.assertEqual([dt[k] for k in ("base", "time_half", "time_quarter")], [2, 1, .5])
        self.assertEqual(dt["space_half"], dt["time_quarter"])
        self.assertEqual(nodes["space_half"]-1, 2*(nodes["base"]-1))

    def test_thermakin_independent_controls(self):
        p = runner.parameters("vhc05")
        decks = {level: runner.thermakin(p, level)[".cnd"] for level in runner.LEVELS}
        values = lambda pattern: {k: float(re.search(pattern, v)[1]) for k, v in decks.items()}
        dt, h = values(r"TIME STEP: (\S+)"), values(r"ELEMENT SIZE: (\S+)")
        duration = values(r"DURATION: (\S+)")
        self.assertEqual([dt[k] for k in ("base", "time_half", "time_quarter")], [.05, .025, .0125])
        self.assertEqual({h[k] for k in ("base", "time_half", "time_quarter")}, {.02/80})
        self.assertEqual((h["space_half"], dt["space_half"]), (.02/160, dt["time_quarter"]))
        for level in runner.LEVELS:
            self.assertEqual(duration[level], p["end"] - dt[level]/2)

    def test_thermakin_boundary_elements(self):
        # A fixed face starts at its fixed value in a one-element top layer.
        cnd = runner.thermakin(runner.parameters("vhc04a"), "base")[".cnd"]
        layers = re.findall(r"THICKNESS: (\S+)\nTEMPERATURE: (\S+)", cnd)
        self.assertEqual([(float(a), float(b)) for a, b in layers[:1]], [(.02/80, 400.0)])
        self.assertRegex(cnd, r"CONVECTION COEFF: 1000000000\n")
        # Only the exposed element emits; absorption saturates every element.
        cmp = runner.thermakin(runner.parameters("vrx04a"), "base")[".cmp"]
        emissivity = dict(re.findall(r"COMPONENT: (\S+)(?:.|\n)*?ABSORPTION: (\S+) 1000000", cmp))
        self.assertEqual(emissivity, {"A": "0", "A_surface": "0.90000000000000002"})
        # Sealed diffusion still needs a listed boundary mass flux.
        cnd = runner.thermakin(runner.parameters("vgt01"), "base")[".cnd"]
        self.assertIn("MASS TRANSPORT: YES\nG LIN 0 0\n", cnd)
        self.assertEqual(cnd.count("THICKNESS:"), 160)

    def test_fds_observation_location(self):
        text = runner.fds(runner.parameters("vrx04a"), "base")
        self.assertIn("CELL_CENTERED=T, FORMAT_INDEX=2", text)
        centers = [float(x) for x in re.findall(r"ID='cell_\d+'.*?DEPTH=([^,]+)", text)]
        self.assertEqual(len(centers), 4)
        for i, x in enumerate(centers):
            self.assertAlmostEqual(x, (i+.5)*1e-4/4)

    def test_fds_faces_are_not_unknown_count_or_probe_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "test.out").write_text(
                "Solid Phase Node, Layer, Coordinates(m):\n"
                " 0 1 0.0000000\n 1 1 0.0005000\n 2 2 0.0010000\n Exposed Backing\n")
            grid = runner.native_grid("fds", "test", path)
            self.assertEqual(grid["intervals"], 2)
            self.assertEqual(grid["spatial_unknowns"], 2)
            self.assertEqual(len(grid["coordinates_m"]), 3)

    def test_common_support_does_not_extrapolate(self):
        d, n, omitted = audit.change([(0, 0), (.5, 1), (1, 2), (2, 4)], [(0, 0), (1, 2)])
        self.assertEqual((d, n, omitted), (0, 3, 1))
        with self.assertRaises(ValueError):
            audit.change([(2, 4)], [(0, 0), (1, 2)])


if __name__ == "__main__":
    unittest.main()
