# Literature Property Comparison Study — Implementation Plan

## Objective

Demonstrate the practical value of the coordinated thermochemical measurements in this
work by running the 1D Douglas fir cone calorimeter simulation (50 kW/m², 12 mm slab,
6% MC) with literature-sourced property values substituted **one category at a time**.
Each comparison holds all other properties at our measured values, isolating the effect
of the literature substitution.

Results are integrated into the existing results sections of the manuscript rather than
forming a standalone section — each property category's comparison lives alongside the
corresponding measurements discussion.

---

## Property Categories and Literature Selections

### 1. Virgin Wood Heat Capacity

**Our baseline (Douglas fir):**

| Range            | Correlation (kJ/kg/K)              |
|------------------|-------------------------------------|
| T ≤ 323 K (HFM)    | c_p = −0.0531 + 4.245 × 10⁻³ T   |
| 323 < T ≤ 425 K    | c_p = −0.0310 + 4.177 × 10⁻³ T   |
| T > 425 K (DSC)    | c_p = −2.131 + 9.119 × 10⁻³ T    |

**Literature substitutions (extrapolated over full temperature range):**

| Label | Source | Correlation (kJ/kg/K) | Original Range | Behaviour vs Ours |
|-------|--------|-----------------------|----------------|-------------------|
| **Low** | Harada et al. (1998), various softwoods | c_p = 1.200 + 2.45 × 10⁻³ (T − 273.15) | 298–523 K | −17% at 450 K |
| **High** | Dupont et al. (2014), avg. softwood | c_p = 5.2 × 10⁻³ T − 0.2847 | 313–353 K | +5–9% in overlap; extrapolation grows |

Both are single-slope linear correlations that depart from our piecewise model primarily
at elevated temperatures (>400 K) where our DSC slope steepens to 9.1 × 10⁻³ kJ/kg/K².

**Note on extrapolation:** These literature correlations were measured over narrower
temperature ranges than the simulation requires. Extrapolating them is part of the point —
a modeller forced to use them would have no choice but to extrapolate, and we are
quantifying the consequence of doing so.

**Where in the paper:** Add comparison plots and 1–2 paragraphs at the end of
Section 3.1 (Heat Capacity of Virgin Wood).

---

### 2. Char Heat Capacity

**Our baseline (Douglas fir char):**

| Correlation (kJ/kg/K)                           | Range       |
|--------------------------------------------------|-------------|
| c_p = 0.3472 + 1.784 × 10⁻³ T                  | 473–773 K   |

**Literature substitutions:**

| Label | Source | Correlation (kJ/kg/K) | Original Range | Behaviour vs Ours |
|-------|--------|-----------------------|----------------|-------------------|
| **Low** | Hostikka & Matala (2017), birch char | c_p ≈ −0.370 + 2.122 × 10⁻³ T (linear fit to 5 tabulated points) | 471–723 K | 25–44% below |
| **High** | Edrich (1980), pine char | c_p = −0.29227 + 4.1387 × 10⁻³ T − 1.58261 × 10⁻⁶ T² | 473–773 K | 11–15% above |

Hostikka & Matala is widely used in FDS fire modeling studies and represents the extreme
low end. Edrich pine char represents the high end for softwood chars.

**Where in the paper:** Add at the end of Section 3.2 (Heat Capacity of Char).

---

### 3. Heats of Reaction (Scaled to Literature Total Heat of Pyrolysis)

**Our baseline (Douglas fir):**

| Reaction | ΔH_j (kJ/kg) | Reactant | Y_i (cumulative yield) | ΔH_j × Y_i (kJ/kg initial) |
|----------|---------------|----------|------------------------|------------------------------|
| 1 (A→B)  | −19.15       | A        | 1.0000                 | −19.15                       |
| 2 (B→C)  | −5.57        | B        | 0.9500                 | −5.29                        |
| 3 (C→D)  | +266.70      | C        | 0.7600                 | +202.69                      |
| 4 (D→E)  | −219.64      | D        | 0.2812                 | −61.76                       |
| 5 (E→F)  | −12.38       | E        | 0.2221                 | −2.75                        |
| **Total** | | | | **ΔH_P = +113.74 kJ/kg** |

**Literature targets (softwood only):**

| Label | Source | Species | ΔH_P target (kJ/kg) | Scale factor f |
|-------|--------|---------|----------------------|----------------|
| **High endothermic** | Rath et al. (2003) | Spruce | +428 | 3.763 |
| **Moderate endothermic** | Jones et al. (2020) | Pine | +382 | 3.359 |

#### How the scaling works

The overall heat of pyrolysis ΔH_P is defined (Equation 12 of the paper) as:

$$
\Delta H_P = \sum_{j=1}^{5} \Delta H_j \cdot Y_{i(j)}
$$

where Y_{i(j)} is the cumulative solid yield of the **reactant** species of reaction j,
computed from the stoichiometric coefficients (ν_j) via:

$$
Y_{i(j)} = \prod_{k=1}^{j-1} \nu_k \qquad (\text{with } Y_A = 1)
$$

The key property is that ΔH_P is **linear** in the individual ΔH_j. If we multiply
every individual heat of reaction by a common scale factor f:

$$
\Delta H_j^{\text{scaled}} = f \cdot \Delta H_j \qquad \text{for all } j = 1, \ldots, 5
$$

then the total heat of pyrolysis scales identically:

$$
\Delta H_P^{\text{scaled}} = \sum_{j=1}^{5} (f \cdot \Delta H_j) \cdot Y_{i(j)}
= f \cdot \sum_{j=1}^{5} \Delta H_j \cdot Y_{i(j)}
= f \cdot \Delta H_P^{\text{ours}}
$$

The cumulative yields Y_{i(j)} depend only on the stoichiometric coefficients (which we
are **not** changing — the kinetics stay fixed), so they factor out cleanly. We choose:

$$
f = \frac{\Delta H_P^{\text{target}}}{\Delta H_P^{\text{ours}}}
$$

This gives us the desired total heat of pyrolysis **exactly**, while preserving the
relative distribution of endothermic and exothermic character among the five reactions.

#### Numerical verification

**Rath et al. target: ΔH_P = +428 kJ/kg**

f = 428 / 113.74 = **3.763**

| Reaction | ΔH_j (ours) | ΔH_j (scaled) | ΔH_j × Y_i (scaled) |
|----------|-------------|----------------|---------------------|
| 1        | −19.15      | −72.06         | −72.06              |
| 2        | −5.57       | −20.96         | −19.91              |
| 3        | +266.70     | +1003.61       | +762.74             |
| 4        | −219.64     | −826.52        | −232.42             |
| 5        | −12.38      | −46.59         | −10.34              |
| **Total** | | | **+428.00 kJ/kg** ✓ |

**Jones et al. target: ΔH_P = +382 kJ/kg**

f = 382 / 113.74 = **3.359**

| Reaction | ΔH_j (ours) | ΔH_j (scaled) | ΔH_j × Y_i (scaled) |
|----------|-------------|----------------|---------------------|
| 1        | −19.15      | −64.32         | −64.32              |
| 2        | −5.57       | −18.71         | −17.77              |
| 3        | +266.70     | +895.74        | +680.76             |
| 4        | −219.64     | −737.69        | −207.44             |
| 5        | −12.38      | −41.58         | −9.23               |
| **Total** | | | **+382.00 kJ/kg** ✓ |

#### Why uniform scaling is appropriate

1. **Preserves the reaction profile shape.** Exothermic reactions stay exothermic;
   endothermic reactions stay endothermic. The temporal distribution of heat
   absorption/release during pyrolysis maintains the same qualitative character as our
   measured curve — only the amplitude changes.

2. **Literature reports only total ΔH_P.** Rath, Jones, and Basile all report a single
   net heat obtained by integrating over the full pyrolysis temperature range. They do not
   decompose it into per-reaction values. Uniform scaling is the minimal assumption that
   maps their total onto our reaction mechanism.

3. **Separates the "how much" from the "when."** Since we keep the kinetic parameters
   (A, E, ν) unchanged, the mass loss profile and species evolution are identical across
   all runs. Only the energy absorbed/released at each moment changes — and it changes
   in proportion to our measured distribution.

**Where in the paper:** Add at the end of Section 3.4 (Heats of Pyrolysis Reactions).

---

### 4. Heat of Sorption

**Our baseline:**

HoS(MC) = 2280 + 787 · exp(−0.114 · MC) kJ/kg
→ At 6% MC, HoS_avg ≈ 2730 kJ/kg

**Literature substitution:**

| Label | Value | Rationale |
|-------|-------|-----------|
| **Free water only** | 2442 kJ/kg (constant, no MC dependence) | Standard fire model assumption; what most FDS/GPYRO users implement |

This is the latent heat of vaporization of water at 25°C. It ignores the bound-water
excess entirely, which at 6% MC represents a ~12% underestimate of the energy required.

**Where in the paper:** Add at the end of Section 3.3 (Water Desorption Energetics),
tying into the existing averaging discussion in Section 3.4.

---

## Simulation Matrix

All simulations use the 1D Douglas fir cone calorimeter configuration from Appendix B:
12 mm slab, 50 kW/m², 6% MC, 600 cells, 1200 s duration.

| Run ID | Property Varied | Substitution | Config Change |
|--------|----------------|--------------|---------------|
| 0 | — (baseline) | Our measured values | None |
| 1 | Virgin c_p | Harada (low) | Replace piecewise virgin c_p with single linear correlation for species A and all intermediates |
| 2 | Virgin c_p | Dupont (high) | Same as above |
| 3 | Char c_p | Hostikka (low) | Replace char c_p for species F and propagate to intermediates |
| 4 | Char c_p | Edrich (high) | Same as above |
| 5 | Heats of reaction | Rath +428 kJ/kg (f = 3.763) | Scale all 5 heat_of_reaction values by f |
| 6 | Heats of reaction | Jones +382 kJ/kg (f = 3.359) | Scale all 5 heat_of_reaction values by f |
| 7 | Heat of sorption | Free water (2442 kJ/kg) | Replace HoS model with constant latent heat |

**Total: 8 simulation runs** (including baseline).

---

## Intermediate Species Heat Capacity — How Substitution Propagates

When swapping virgin c_p or char c_p, the intermediate species B–E must also be updated
because their heat capacities are interpolated between virgin and char via the degree of
conversion α_i (Equation 10 of the paper):

c_{p,i} = (1 − α_i) · c_{p,virgin} + α_i · c_{p,char}

The conversion fractions for Douglas fir are:

| Species | α_i | Weight on virgin | Weight on char |
|---------|------|-----------------|----------------|
| A       | 0.000 | 1.000          | 0.000          |
| B       | 0.060 | 0.940          | 0.060          |
| C       | 0.287 | 0.713          | 0.287          |
| D       | 0.860 | 0.140          | 0.860          |
| E       | 0.930 | 0.070          | 0.930          |
| F       | 1.000 | 0.000          | 1.000          |

So swapping virgin c_p affects A heavily, B strongly, C moderately, and has almost no
effect on D, E, F. Conversely, swapping char c_p affects F heavily, E and D strongly,
and has negligible effect on A and B.

---

## Outputs and Plotting

### Quantities of Interest (per the existing GSA)

1. Peak mass loss rate (kg/m²/s)
2. Time to peak mass loss rate (s)
3. Total mass loss (%)
4. Peak surface temperature (K)
5. Peak back face temperature (K)
6. Final thickness (m)

### Figures

For each property category, a **single 2-panel figure**:
- **Top:** Mass loss rate vs time — baseline (solid black) + 2 literature curves (colored dashed)
- **Bottom:** Back face temperature vs time — same styling

Optionally overlay the 5th–95th percentile prediction envelope from the existing PCE/GSA
as a shaded band behind the baseline, showing where the literature curves fall relative to
our measurement uncertainty.

### Summary Table

A single summary table (possibly at the end of Section 3.5 or in the discussion) showing
the % deviation from baseline for each QoI and each substitution:

| Substitution | Δ(MLR_peak) | Δ(t_peak) | Δ(m_total) | Δ(T_s,max) | Δ(T_b,max) | Δ(ℓ_f) |
|-------------|-------------|-----------|------------|------------|------------|---------|
| Virgin c_p — Harada | ? | ? | ? | ? | ? | ? |
| Virgin c_p — Dupont | ? | ? | ? | ? | ? | ? |
| Char c_p — Hostikka | ? | ? | ? | ? | ? | ? |
| Char c_p — Edrich | ? | ? | ? | ? | ? | ? |
| ΔH — Rath (+428) | ? | ? | ? | ? | ? | ? |
| ΔH — Jones (+382) | ? | ? | ? | ? | ? | ? |
| HoS — Free water | ? | ? | ? | ? | ? | ? |

---

## Implementation Steps

1. **Create variant TOML configs** (or a Julia driver script that loads the baseline
   config and programmatically overrides the relevant properties for each run).

2. **Run all 8 simulations** and extract the 6 QoIs from each.

3. **Generate the 4 comparison figures** (one per property category, 2 panels each).

4. **Generate the summary table.**

5. **Write manuscript text** for each section (1–2 paragraphs per property category
   interpreting the results).

---

## Expected Narrative Arc

- **Heats of reaction:** Dominant effect (consistent with GSA). Literature ΔH_P values of
  +382–428 kJ/kg (3.4–3.8× our value) should dramatically slow the pyrolysis front, reduce
  peak MLR, and lower peak temperatures. This is the "if you get this wrong, nothing else
  matters" result.

- **Virgin c_p:** Moderate effect despite large literature spread. Harada's weaker
  temperature dependence means less energy required to heat the virgin material to
  pyrolysis onset → faster pyrolysis front → higher MLR. Dupont's steeper slope has the
  opposite effect. The key point: a single-slope correlation misses the steepening above
  425 K.

- **Char c_p:** Small to moderate effect. Hostikka's much lower values mean the char layer
  behind the pyrolysis front heats up faster (less thermal inertia in the insulating layer),
  potentially affecting back face temperature timing. But char is already a small mass
  fraction so the overall energy budget impact is modest.

- **Heat of sorption:** Negligible effect (consistent with GSA and Section 3.4's averaging
  result). Using 2442 kJ/kg instead of ~2730 kJ/kg represents a 12% underestimate of
  drying energy, but drying completes before pyrolysis onset, so the downstream effect is
  small. This is a reassuring negative result — modellers can use the free-water assumption
  without significant penalty.
