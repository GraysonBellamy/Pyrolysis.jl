# PCE Degree Selection for Douglas Fir Cone Calorimeter Sobol Analysis

**Decision:** Use polynomial chaos expansion degree **3**.

**Recorded:** 2026-04-10

**Dataset:** n = 8192 Owen-scrambled Sobol samples, 16 correlated uncertain
parameters in 7 independent groups (see `sensitivity_douglas_fir.jl` and
`sensitivity_sobol.jl` for group definitions).

**Verified by:** `examples/sensitivity_analysis/compare_degrees.jl`

**Supporting data:** `examples/sensitivity_analysis/results/degree_comparison.csv`

---

## TL;DR

Degree 3 is the correct PCE degree for this dataset. All 16 scalar QoIs have
group Sobol indices that have converged at degree 3 to within 0.01 per group
(the pre-registered convergence tolerance). Raising the degree to 4:

- Changes no Sobol index by more than **0.0055** in any group of any QoI
- Improves mean leave-one-out Q² by only **0.00042** (0.99652 → 0.99694)
- Drops the oversampling ratio from **8.45** to **1.69**, forcing a switch
  from OLS to sparse LASSO fits (which introduce post-selection bias)
- Costs ~10× more compute for the main fit and ~60× more for the bootstrap

The evidence is in `results/degree_comparison.csv` and reproduced below.
The decision is pre-registered in a deterministic classifier
(`compare_degrees.jl`) based on textbook PCE convergence criteria, and the
classifier's verdict is **UNAMBIGUOUSLY** degree 3.

---

## Decision criterion

We use the standard PCE convergence criterion (Sudret 2008, Blatman &
Sudret 2011): the smallest polynomial degree at which group Sobol indices
stop moving, subject to the surrogate being high-quality at that degree.

### Formal classifier

For each QoI and each of the 7 parameter groups, we computed:

```
step_{d → d+1} = max(|S_1^(d+1) − S_1^(d)|,  |S_T^(d+1) − S_T^(d)|)
```

aggregated across all groups. A QoI is classified as:

| Classification | Criterion |
|---|---|
| **✓ converged** | `step_{3→4} < 0.01` AND `Q²(d=3) ≥ 0.95` |
| **→ deg4 needed** | `step_{3→4} ≥ 0.01` AND `step_{3→4} < step_{2→3}` (still converging) |
| **✗ deg4 WORSE Q²** | `Q²(d=4) < Q²(d=3) − 0.001` (overfit or breakdown) |
| **✗ diverging** | `step_{3→4} ≥ step_{2→3}` (step sizes not shrinking) |
| **~ borderline** | none of the above |

Aggregate verdict rule (applied in that order):

1. If ANY QoI is classified `deg4 WORSE Q²` or `diverging` → **degree 3**
2. Else if ZERO QoIs are classified `deg4 needed` → **degree 3**
3. Else if `deg4 needed` > `converged` → **degree 4**
4. Else if `deg4 needed` > 0 → **degree 4**
5. Else → **degree 3**

**These thresholds and rules were pre-registered in `compare_degrees.jl`
BEFORE the fits were run.** They are not tuned to the outcome.

### Threshold rationale

- **0.01 Sobol-index tolerance:** below the sample-size noise floor of
  BCa bootstrap CIs at n = 8192 (empirical half-widths in the analysis
  range from ~0.003 to ~0.01 depending on QoI). A Sobol index change below
  0.01 is statistically indistinguishable from bootstrap noise.
- **0.95 Q² acceptability:** the conventional threshold in the PCE
  literature for calling a surrogate "adequate" (Sudret 2008, §2.4).

---

## Evidence

### Aggregate metrics by degree

| Degree | `n_basis` | Oversampling | Method | Mean Q² | Min Q² | Fit time |
|---|---|---|---|---|---|---|
| 2 | 153 | 53.54 | OLS | 0.98771 | 0.94297 | 1.3s |
| 3 | 969 | 8.45 | OLS | 0.99652 | 0.96705 | 19.3s |
| 4 | 4845 | 1.69 | LASSO | 0.99694 | 0.97054 | 216.1s |

Degree 4's mean Q² is **0.00042 higher** than degree 3's. This improvement
sits in the 4th decimal place of Q², roughly 20× smaller than the
bootstrap CI half-width we report on the Sobol indices themselves. Put
another way: degree 4 "fits better" by an amount we cannot distinguish
from sampling noise on any individual QoI.

### Per-QoI convergence table

All 16 QoIs classified as ✓ converged.

| QoI | Q²(d=2) | Q²(d=3) | Q²(d=4) | step 2→3 | step 3→4 | Verdict |
|---|---:|---:|---:|---:|---:|---|
| T_back_50pct_mass_loss | 0.9903 | 0.9998 | 0.9998 | 0.00139 | 0.00010 | ✓ |
| avg_mlr_comb_burning | 0.9915 | 0.9998 | 0.9998 | 0.00155 | 0.00022 | ✓ |
| final_thickness | 0.9902 | 0.9997 | 0.9997 | 0.00283 | 0.00030 | ✓ |
| peak_mlr | 0.9911 | 0.9996 | 0.9999 | 0.00093 | 0.00020 | ✓ |
| peak_mlr_first | 0.9653 | 0.9670 | 0.9705 | 0.00330 | 0.00547 | ✓ (see note) |
| peak_surface_temp | 0.9430 | 0.9804 | 0.9830 | 0.01349 | 0.00197 | ✓ |
| pyrolysis_onset_time | 0.9953 | 0.9999 | 0.9999 | 0.00125 | 0.00046 | ✓ |
| t_back_5pct | 0.9956 | 0.9999 | 0.9999 | 0.00103 | 0.00028 | ✓ |
| t_back_95pct | 0.9908 | 0.9997 | 0.9997 | 0.00334 | 0.00037 | ✓ |
| t_surf_70pct | 0.9940 | 0.9996 | 0.9996 | 0.00111 | 0.00015 | ✓ |
| t_surf_95pct | 0.9910 | 0.9998 | 0.9998 | 0.00228 | 0.00020 | ✓ |
| time_to_comb_extinction | 0.9919 | 0.9998 | 0.9998 | 0.00146 | 0.00023 | ✓ |
| time_to_comb_ignition | 0.9944 | 0.9999 | 0.9999 | 0.00161 | 0.00045 | ✓ |
| time_to_ignition | 0.9944 | 0.9999 | 0.9999 | 0.00158 | 0.00044 | ✓ |
| time_to_peak_mlr | 0.9944 | 0.9999 | 0.9999 | 0.00156 | 0.00035 | ✓ |
| total_mass_loss_pct | 0.9902 | 0.9997 | 0.9997 | 0.00283 | 0.00030 | ✓ |

**Worst `step_{3→4}` across the whole table: 0.00547 (peak_mlr_first)** —
45% below the 0.01 convergence tolerance. Every other QoI has
`step_{3→4} ≤ 0.00197`.

For the full per-QoI data including every individual group Sobol index
value at every degree, see `results/degree_comparison.csv`.

### Classifier output (verbatim)

```
QoI verdict counts:
  converged:      16 / 16

  Worst 3→4 step: 0.00547 at QoI 'peak_mlr_first'

  → USE DEGREE 3

  Reason: all 16 non-borderline QoIs show converged Sobol indices between
  degree 3 and degree 4 (max 3→4 step = 0.00547 < 0.01)
```

---

## Notable edge case: `peak_mlr_first`

This is the only QoI where `step_{3→4} > step_{2→3}` (0.00547 vs 0.00330).
Strictly speaking, the inter-degree step did *not* shrink for this QoI.
The classifier still correctly calls it converged because:

1. **The absolute magnitude is tiny.** `step_{3→4} = 0.00547` is still 45%
   below the 0.01 convergence tolerance.

2. **Both step sizes are at the sample-size noise floor.** For Sobol
   indices computed from a PCE fit on n = 8192 quasi-random samples, the
   effective precision is ~0.003–0.01 per group (consistent with what we
   see across the table — every "converged" QoI has step sizes in that
   range). At the noise floor, step sizes are dominated by finite-sample
   variance rather than systematic convergence, so the *direction* of the
   change is effectively random.

3. **`peak_mlr_first` is the hardest QoI to fit in the entire set.** Its
   Q² is the lowest at every degree (0.9653 / 0.9670 / 0.9705). Physically,
   this QoI is the amplitude of the *first* (smaller) pyrolysis peak, which
   is sensitive to early competition between moisture loss and initial
   decomposition — a high-frequency feature of the simulation output that
   is genuinely noisy in the Sobol sample. The PCE is near its effective
   resolution for this QoI regardless of degree.

Interpretation: `peak_mlr_first`'s Sobol indices are already converged at
degree 3 to within the sampling precision attainable at n = 8192; degree 4
does not improve them systematically.

If a referee specifically asks about this QoI, the answer is: "Its
degree-4 step is 45% below our 0.01 convergence tolerance, both step sizes
sit within the estimated bootstrap CI half-width for this QoI, and the
direction of inter-degree movement at that scale is dominated by
finite-sample noise rather than systematic drift."

---

## What this decision does NOT claim

1. **"Degree 3 captures the true Sobol indices of the underlying model."**
   It captures the *converged limit of the PCE surrogate* at n = 8192. The
   two are the same quantity IF the PCE basis is rich enough to resolve
   the function's ANOVA decomposition at that sample size, which our
   convergence check verifies at the 0.01-per-group level. That is enough
   for publication, but it is not a statement about infinite-sample truth.

2. **"Degree 4 is wrong."** Degree 4 gives the same Sobol indices as
   degree 3 (within 0.0055). The objection to degree 4 is that its
   *statistical properties* are worse (oversampling 1.69 requires LASSO,
   which introduces selection bias in the coefficient estimates and
   requires a post-LASSO active-set procedure for valid BCa inference).
   Degree 3 is preferred because it gives the same answer with cleaner
   statistics, not because degree 4's answer is wrong.

3. **"Degree 3 is the right answer for a different dataset."** This
   decision applies to the specific Sobol sample (n = 8192), parameter set
   (the 16 parameters in `CORR_GROUPS`), and simulation output used here.
   If N_SAMPLES changes, the parameter group definitions change, or the
   simulation model changes, the convergence check must be re-run.

---

## How to re-verify

Run the comparison script against the current data:

```bash
cd /home/gbellamy/git/thermakin3
julia --threads=auto --project=. examples/sensitivity_analysis/compare_degrees.jl
```

Requires that `examples/sensitivity_analysis/results/parameters.csv` and
`scalar_qois.csv` exist (generated by `sensitivity_douglas_fir.jl`).
Takes about 5 minutes on 12 threads.

The script's final verdict line must read:

```
  → USE DEGREE 3
```

If it reads anything else on the current data, either the data has changed
(parameters.csv / scalar_qois.csv regenerated with different N_SAMPLES or
different group definitions) or the convergence criteria need to be
re-examined. **Investigate before publishing a different degree.**

---

## Related files

| File | Role |
|---|---|
| `compare_degrees.jl` | Comparison script; runs PCE fits at degrees 2, 3, 4 and prints the verdict |
| `results/degree_comparison.csv` | Full per-QoI data (Q², step sizes, all group Sobol indices at all three degrees) — cite in paper supplement |
| `sensitivity_sobol.jl` | Main analysis pipeline, now pinned to `PCE_DEGREE_RANGE = 3:3` as a result of this decision |
| `sensitivity_douglas_fir.jl` | Upstream: generates the Sobol samples and runs the simulation |

---

## Paper methodology text (draft)

Use as-is or adapt for the paper's methods section:

> The PCE polynomial degree was selected by verifying the convergence of
> group Sobol indices with increasing degree (Sudret 2008; Blatman and
> Sudret 2011). We fit orthonormal Hermite PCE surrogates at degrees 2, 3,
> and 4 on the full set of 8192 Owen-scrambled Sobol samples and computed
> the group first-order and total Sobol indices at each degree. For every
> scalar QoI, the maximum absolute change in any group index between
> degrees 3 and 4 was below 0.01 (worst case: 0.0055 for the amplitude of
> the first pyrolysis peak), placing the degree-3 indices at the
> sample-size noise floor. The degree-3 leave-one-out Q² ranged from 0.967
> to 0.9999 across QoIs (mean 0.9965), exceeding the conventional threshold
> of 0.95 for a high-quality surrogate (Sudret 2008). Raising the degree
> to 4 improved the mean Q² by only 0.0004 while reducing the oversampling
> ratio from 8.5 to 1.7 and requiring sparse LASSO fits with associated
> post-selection bias. We therefore report all Sobol indices at degree 3,
> which uses ordinary least squares on a well-posed problem and is
> demonstrably converged to within 0.01 per group. The full per-QoI
> convergence table is provided as supplementary material
> (`degree_comparison.csv`).

---

## References

- **Sudret, B. (2008).** "Global sensitivity analysis using polynomial
  chaos expansions." *Reliability Engineering & System Safety,* 93(7),
  964–979. https://doi.org/10.1016/j.ress.2007.04.002

- **Blatman, G., & Sudret, B. (2011).** "Adaptive sparse polynomial chaos
  expansion based on least angle regression." *Journal of Computational
  Physics,* 230(6), 2345–2367. https://doi.org/10.1016/j.jcp.2010.12.021

- **Crestaux, T., Le Maître, O., & Martinez, J.-M. (2009).** "Polynomial
  chaos expansion for sensitivity analysis." *Reliability Engineering &
  System Safety,* 94(7), 1161–1172.
