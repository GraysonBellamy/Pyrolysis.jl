# Development

This page covers working on Pyrolysis.jl itself. For using the package, start
with the [User Guide](guide/user_guide/01_introduction-installation.md).

## Repository layout

| Path | Contents |
|:---|:---|
| `src/` | The package, organised into submodules in load order: `Materials`, `Properties`, `BoundaryConditions`, `Physics` (the shared physics layers), `Meshes`, `States`, `Kernels`, `Closures`, `Assembly`, `Deactivation`, `ColoredJacobian` (the discretization), `Problem`, `Outputs`, `Integration`, plus the `Internal` umbrella. |
| `ext/` | Package extension (Plots recipes). |
| `test/` | Test suite mirroring the `src/` module layout (a directory per module, plus `jacobian/`, `multid/`, `type_analysis/` and `code_scan/`), the verification cases (`verification/`), and measurement scripts (`gates/`); `runtests.jl` walks the directories. |
| `benchmark/` | PkgBenchmark.jl regression suite over the canonical problems. |
| `examples/` | End-to-end studies with their own project environment. |
| `docs/` | The two canonical guide monoliths, the split chapters (`docs/guide/`), and the Documenter site (`docs/make.jl`). |
| `scripts/` | Maintenance scripts (formatting, guide splitting, profiling). |

## Running the tests

```julia
using Pkg
Pkg.test("Pyrolysis")
```

`test/runtests.jl` runs every `test_*.jl` file under the module directories,
each in a module of its own, the type analysis in a fresh process, and the
verification cases (`test/verification/runcases.jl`). The module suites test
both energy forms themselves; `PYROLYSIS_ENERGY_STATE=enthalpy` runs the
verification cases in the enthalpy form. `PYROLYSIS_TEST_GROUP` selects one
part (`jacobian`, `burnout`, `other`, `domains`, or `all`, the default) and
`PYROLYSIS_TEST_FAST=1` skips the verification cases' convergence spot
checks; the `tests.yml` workflow runs the parts side by side on every pull
request to `main`, and the whole suite weekly and on dispatch.

## Quality tooling

Formatting and static-quality checks run in isolated environments (the CI
`quality.yml` workflow runs the same commands):

```bash
julia --project=scripts scripts/format.jl --check   # JuliaFormatter (gating)
julia --project=test/quality test/quality/aqua.jl   # Aqua package hygiene
julia --project=test/quality test/quality/jet.jl    # JET analysis (advisory)
```

Run `scripts/format.jl` without `--check` to apply formatting.

## Benchmarks

`benchmark/benchmarks.jl` defines a PkgBenchmark.jl `SUITE` covering residual
evaluation, the colored Jacobian, and end-to-end solves on the canonical
configurations (PMMA cone, white pine, charring polymer). See
[`benchmark/README.md`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/benchmark/README.md)
for `benchmarkpkg`/`judge` usage.

## Working on the guides

The two monoliths — `docs/User_Guide.md` and
`docs/Technical_Reference_Guide.md` — are the canonical, hand-edited
documents. Two derived views are generated from them:

1. `docs/guide/**` (GitHub-browsable chapters) — regenerate after editing a
   monolith with:
   ```bash
   julia scripts/split_guides.jl
   ```
2. `docs/src/guide/**` (Documenter site sources) — generated automatically by
   the docs build; never edit these by hand. The generator
   (`docs/generate_guide.jl`) also converts `$`-delimited math to Documenter
   syntax and fails the build if anything is left unconverted.

### Explanatory figures

`docs/assets/figures/` is the canonical asset directory. The SVGs retain editable
text and vector shapes; their reproducible drawing sources are in
`scripts/figures/`. Edit those Julia sources before regenerating the SVG exports.
No plotting package or raster image generator is needed. Use a technical-report
style: white backgrounds, square outlines, fine rules, conventional label
capitalization, and explicit mathematical subscripts and superscripts. Preserve
the case of physical symbols and API identifiers; render API identifiers in
monospace. Mesh connectivity is black, with restrained blue for thermal fields
and numerical entries, teal for inventory transport and correction, brown for
heat, muted purple for global coupling, and grey for reference/inactive material.
Schematics say so in their titles and captions; plots identify computed
quantities and units.

Each figure should make one visual comparison or mechanism clear beside the
text. Prefer direct labels, aligned geometry and transfer arrows. Keep derivations,
API details, procedural qualifications and explanatory paragraphs in the guide
or caption; do not repeat them as text panels inside an image.

| Figure | Subject | User Guide | Technical Reference |
|:---|:---|:---|:---|
| 1 | Finite-volume mesh anatomy | §6.1 | §13.1 |
| 2 | Mesh framework and residual evaluation | reference in §6.1 | §13.4 |
| 3 | Fixed, Lagrangian and corrected ALE | §6.7 | §11.1 |
| 4 | Equidistribution and thickness activation | §10.5 | §11.7 |
| 5 | Burnout, transfers and retained state slots | §10.2 | §14.3 |
| 6 | Surface energy balance | §5.4 | §12.2 |
| 7 | Surface versus in-depth absorption | §7.4 | §8.1 |
| 8 | Actual quick-start results | §2.9; compact copy on Home | — |
| 9 | Constructor coordinates | §6.2 | — |
| 10 | One-dimensional coordinate geometry | §6.4 | — |
| 11 | Axial and lateral shrinkage/swelling | §7.6 | §10.6 |
| 12 | Porous gas transport | §7.5 | Chapter 9 |
| 13 | Three-step wood reaction network | §4.7 | — |
| 14 | Mass and energy accounting | §9.8 | §16.5 |
| 15 | State layout and actual Jacobian patterns | — | §15.1; reference in §15.9 |

After setting up the documentation environment below, regenerate from the
repository root:

```bash
julia --project=docs scripts/figures/generate.jl
julia scripts/test_guide_assets.jl
julia scripts/split_guides.jl
julia --project=docs docs/make.jl
```

The figure generator supports `--first` (1–8) and `--second` (9–15). Figure 8
executes the complete runnable block in canonical User Guide §2.10, omitting
only the optional Plots import and progress display, then saves the actual
boundary MLR, gas-generation proxy, temperatures and current positions in
`docs/assets/figures/data/`. Its TOML provenance records the Julia version,
source revision, executed-deck hash and solve outcome. For appearance-only
changes, `--from-data` redraws that figure from the committed CSVs without
another solve. Figure 9 exports the actual constructor coordinates. Figure 15
exports the stored sparsity entries and initial dense-AD derivatives for three
six-cell setups, and asserts that every nonzero derivative lies in the stored
pattern. Figure 4 evaluates the implemented thickness-activation switch;
Figure 7's exponential inset illustrates the labelled uniform optical depth.

Canonical Markdown embeds `assets/figures/name.svg`. The shared
`rebase_figure_paths` helper rewrites those destinations to
`../../assets/figures/name.svg` for both chapter views, leaving code examples
and external URLs untouched. The Documenter generator copies only SVG exports
to the ignored `docs/src/assets/figures/`; numerical data and drawing sources
remain in their canonical locations. The Home page uses that copied asset
directory too. Do not hand-edit either set of derived chapters or copied SVGs.

After regeneration, inspect the SVGs and representative rendered pages at
their normal display width. `rsvg-convert` or a browser can render SVGs for
inspection. The focused asset checks cover rebasing (including fenced and
inline code) and the existence of all fifteen canonical figure destinations;
the documentation build also exercises the Markdown conversion and asset copy.

## Building the documentation site

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The site is written to `docs/build/`. On pushes to `main`, the
`docs.yml` workflow builds the site and deploys it to Cloudflare Pages at
[pyrolysis.graysonbellamy.dev](https://pyrolysis.graysonbellamy.dev); pull
requests build the site for validation without deploying.
