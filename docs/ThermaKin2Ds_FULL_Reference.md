
# Technical Reference for Numerical Pyrolysis Model ThermaKin2Ds

**May 2019 Edition**  
**Stanislav I. Stoliarov**  
University of Maryland, Department of Fire Protection Engineering  
College Park, MD 20742, USA

---

## Summary

This document provides a description of the governing equations, solution methodology and applications of a numerical pyrolysis solver, ThermaKin2Ds. ThermaKin2Ds is a result of further development of ThermaKin2D solver [1], which, in turn, is an extension of ThermaKin [2].

The main purpose of ThermaKin2Ds is to enable a detailed analysis and quantitative prediction of the processes that take place inside and at the surface of a burning material in a wide range of fire scenarios. The solver accounts for chemical reactions (chemical mechanisms of essentially arbitrary complexity can be modeled), transient conduction (governed by temperature dependent properties), in-depth thermal radiation absorption and emission, mass transport of gaseous species through the condensed phase and expansion or contraction of the pyrolyzing object.

ThermaKin2Ds can be used to model both anaerobic and oxidative pyrolysis decoupled from the gas-phase (flaming) combustion. ThermaKin2Ds can also be used to model interaction between pyrolyzing material and gas-phase flame, including various modes of flame spread. The gas-phase flame is represented in the ThermaKin2Ds through an empirical submodel that relates the rate of production of combustible gases by a material object to the spatially resolved heat feedback that this object receives from the flame fueled by these gases.

In addition to the solver description, this document contains examples of ThermaKin2Ds input and output used to model thermal analysis and standard material flammability tests. The section of the document containing these examples is continuously updated to include new results.

---

## 1. Governing Equations

In ThermaKin2Ds, the material is represented by a mixture of components. By default, the maximum number of components is restricted to 50. Every component is characterized by density, heat capacity, thermal conductivity, gas transfer coefficient, emissivity, and radiation absorption (or extinction) coefficient.

The first four properties are defined by a flexible function of temperature:

property(T) = p0 + p1·T + p2·T^n

Emissivity and absorption coefficient are constant. All components are divided into three categories: solids, liquids, and gases.

The specific heat of material is calculated as:

c = Σ mf_i · c_i

The density of material is defined by:

ρ = 1 / ( Σ_s mf_s/ρ_s + Σ_l mf_l/ρ_l + γ Σ_g mf_g/ρ_g )

The swelling factor γ describes reaction of a volume of material to the presence of gases. γ is computed by volume-weighted averaging of solid and liquid swelling parameters. A small τ ≪ 1 ensures convergence to gas density at high gas content.

---

### 1.1 Chemical Reactions

Components may undergo chemical reactions, which can also represent physical transitions such as melting. Up to 50 reactions may be specified.

General reaction form:

θ1 COMP1 + θ2 COMP2 → θ3 COMP3 + θ4 COMP4 + h

Reaction rate:

r = A · exp(-E / (R·T)) · ξ_COMP1 · ξ_COMP2

where ξ = mf · ρ. If the second reactant is absent, ξ_COMP2 = 1.

Reaction heat h is temperature-dependent and defined using the same polynomial form as material properties. Positive h denotes an exothermic reaction. Reaction parameters are not automatically constrained to satisfy mass or energy conservation.

Lower and upper temperature limits may be specified for each reaction.

---

### 1.2 Heat and Mass Transfer

Heat conduction is described by Fourier’s law:

q_x = -k · ∂T/∂x

Thermal conductivity depends on component distribution. For layers normal to heat flow:

k_n = ( Σ v_i / k_i )⁻¹

For layers parallel to heat flow:

k_p = Σ v_i · k_i

Effective conductivity:

k = β k_p + (1 − β) k_n

---

Thermal radiation from an external source is absorbed according to a generalized Beer–Lambert law:

∂I_ex/∂x = -I_ex Σ α_i ξ_i

To comply with the second law of thermodynamics, the material re-radiates energy:

∂I_rr/∂x = σ T⁴ Σ ε_i ξ_i

Radiative transfer inside the material is modeled using a radiation diffusion approximation.

---

Gas mass transfer is driven by concentration gradients:

J_g,x = -λ ρ ∂(ξ_g/ρ)/∂x

Using Boyle’s law:

J_g,x = -λ ρ φ ∂P_g/∂x

If γ > 0, transport is diffusion-like; if γ = 0, transport follows Darcy’s law.

---

### 1.3 Conservation Equations

Mass conservation (Cartesian):

∂ξ_j/∂t =
- Σ θ_ij r_i
- ∂J_j,x/∂x
- ∂J_j,y/∂y
- deformation term

Mass conservation (axisymmetric):

∂ξ_j/∂t =
- Σ θ_ij r_i
- (1/z) ∂(z J_j,z)/∂z
- ∂J_j,x/∂x
- deformation term

---

Energy conservation includes conduction, radiation absorption and emission, gas convection, reaction heat, and expansion/contraction work. Local gas–solid thermal equilibrium is assumed.

---

## 2. Boundary Conditions

ThermaKin2Ds supports three object geometries:

- 1D
- 2D
- 2Dax (axisymmetric)

Boundary mass flux of component i:

J_i^B = a_i (ξ_i^B − b_i ρ^B)
or
J_i^B = a_i ξ_i^B exp(−b_i / RT)

---

### 2.1 1D Object

The 1D object is defined along x with top and bottom boundaries. x = 0 is the reference plane for deformation.

Convective heat flux:

q^B = h_c (T^B − T_e)

Environmental temperature is defined by a multi-parameter heating rate function capable of reproducing thermal analysis temperature histories.

External radiative heat flux is defined using a piecewise-linear function of time.

Surface flame ignition criterion:

Σ (J_i^B / J_i,CI) ≥ 1

Upon ignition, flame radiation and convection parameters replace external boundary conditions.

---

### 2.2 2D Object

The 2D object is defined in Cartesian coordinates x and y. Heat exchange occurs at front and back boundaries; top and bottom boundaries are adiabatic.

External heat flux is defined using a piecewise-linear spatial distribution in y, optionally time-dependent.

Surface flame heat feedback is modeled using empirical correlations based on flame height derived from integrated ignition criterion.

A background temperature Tb defines equilibrium conditions and background radiation.

---

### 2.3 2Dax Object

The 2Dax object is defined in cylindrical coordinates x and z. Radial variation of surface area and exposure is explicitly accounted for.

External heating is defined using spatial functions of z and time-dependent modules.

CAPA II–specific convection and radiation models are supported, including angular dependence of incident radiation.

Background temperature Tb is used to define equilibrium radiation and convection.

---

## 3. Solution Methodology

Material objects are discretized into finite volumes.

- 1D: Δx
- 2D: Δx, Δy
- 2Dax: Δx, Δz

Primary state variables are component masses and temperature.

Time integration:
- Crank–Nicolson in x
- Explicit in y or z

Typical values:
- Δx ≈ 5×10⁻⁵ m
- Δt ≈ 0.01 s

Stability condition for explicit dimensions:

Δt < Δy² / (2A)

Dynamic remeshing is used to control element size during expansion or contraction.

Radiation absorption is handled using either:
- Maximum absorption approximation
- Monte Carlo random absorption

---

## 4. Examples of Simulations: Input and Output

ThermaKin2Ds is distributed as a console executable. It reads:

- Components file (.cmp)
- Conditions file (.cnd)
- Output file (.txt)

All files are ASCII text and case-sensitive.

---

### 4.1 Thermal Analysis (TGA, DSC, MCC)

Thermal analysis experiments are simulated using a 1D object with a single element.

Key assumptions:
- Very high convection coefficient
- Instantaneous gas removal
- Zero emissivity

**Figure 1 (textual description):**  
Example PMMA.cmp file defining component properties, reaction parameters, emissivity, absorption coefficients, swelling parameters, and Arrhenius kinetics.

**Figure 2 (textual description):**  
Example TA_PMMA.cnd file defining object geometry, heating rate, boundary conditions, and integration parameters.

**Figure 3 (textual description):**  
Excerpt of output file showing time, temperature, mass fractions, concentrations, mass loss, and heat flow.

**Figure 4 (textual description):**  
MATLAB visualization tool showing mass loss rate, HRR, integrated heat flow, and temperature profiles.

---

### 4.2 Cone Calorimetry

Cone calorimetry simulations include radiant heater flux and flame heat feedback.

Insulation is modeled explicitly using Kaowool components with negligible gas transport.

**Figure 5 (textual description):**  
Component definitions for Kaowool insulation layers.

**Figure 6 (textual description):**  
Example CC_PMMA.cnd file for a 49.7 kW·m⁻² cone calorimetry simulation.

**Figure 7 (textual description):**  
MATLAB post-processing comparing simulated HRR and mass loss rate with experimental data.

---

## References

1. Stoliarov S. I., Leventon I. T., Lyon R. E., *Two-dimensional Model of Burning for Pyrolyzable Solids*, Fire and Materials, 38, 391–408 (2014).
2. Stoliarov S. I., Lyon R. E., *Thermo-Kinetic Model of Burning*, FAA Technical Note DOT/FAA/AR-TN08/17 (2008).
3. Blinder S. M., *Advanced Physical Chemistry*, Macmillan, 1969.
4. Holman J. P., *Heat Transfer*, 9th Ed., McGraw-Hill, 2002.
5. Siegel R., Howell J., *Thermal Radiation Heat Transfer*, 4th Ed., Taylor & Francis, 2002.
6. Scheidegger A. E., *The Physics of Flow through Porous Media*, University of Toronto Press, 1974.
7. Swann J. D. et al., *CAPA II*, Fire Safety Journal, 91, 130–139 (2017).
8. McKinnon M. B., Stoliarov S. I., Materials, 8, 6117–6153 (2015).
9. Leventon I. T. et al., Combustion and Flame, 179, 338–353 (2017).
10. Lannon C. M. et al., Fire and Materials, 42, 710–724 (2018).
11. Press W. H. et al., *Numerical Recipes in C++*, Cambridge, 2002.
12. Swann J. D. et al., IJHMT, 132, 347–361 (2019).
13. Lyon R. E. et al., FAA Report DOT/FAA/TC-12/53 (2013).
14. ASTM E1354-16a.
15. McCoy C. G. et al., Fire Safety Journal, 103, 38–48 (2019).
