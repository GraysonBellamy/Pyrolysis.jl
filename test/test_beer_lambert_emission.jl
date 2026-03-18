"""
Quantitative analysis: Should Beer-Lambert in-depth absorption include re-radiation (emission)?

This script derives the answer from radiation physics first principles (the 1D RTE)
and computes the magnitude of the omitted emission term under pyrolysis conditions.

Key equations
=============
The 1D radiative transfer equation for an absorbing-emitting gray medium is:

    dI/dz = -alpha * I + alpha * Ib(T)                                   (RTE)

where Ib(T) = sigma*T^4/pi is the Planck function (isotropic blackbody intensity).
Kirchhoff's law gives emission coefficient = absorption coefficient for a gray body.

Beer-Lambert retains only the first (absorption) term:

    dI/dz = -alpha * I                                                   (BL)

The question is: when is the emission term alpha*Ib(T) negligible compared to alpha*I?

    emission/absorption = Ib(T) / I(z) = sigma*T^4 / (pi * I(z))

For a collimated external beam I(z) = I0*exp(-alpha*z), the ratio becomes:

    ratio(z) = sigma*T(z)^4 / (pi * I0 * exp(-alpha*z))

This ratio grows with depth (as I decays) and with temperature.
"""

using Printf

# ============================================================================
# Physical constants
# ============================================================================
const sigma = 5.670374419e-8  # W/(m^2 K^4)

# ============================================================================
# 1. Emission-to-absorption ratio for collimated Beer-Lambert beam
# ============================================================================

"""
    emission_absorption_ratio(T, I0, alpha, z)

Ratio of the RTE emission term to the absorption term for a collimated beam:
    ratio = sigma*T^4 / (pi * I0 * exp(-alpha*z))

When ratio << 1, the emission term is negligible and Beer-Lambert is valid.
When ratio ~ 1, emission is comparable and Beer-Lambert is deficient.
"""
emission_absorption_ratio(T, I0, alpha, z) = sigma * T^4 / (pi * I0 * exp(-alpha * z))

println("=" ^ 80)
println("ANALYSIS: Beer-Lambert emission term for pyrolysis conditions")
println("=" ^ 80)

# ============================================================================
# 2. Blackbody emission flux vs. external flux
# ============================================================================
println("\n--- Section 1: Blackbody emission flux sigma*T^4 vs external flux ---\n")
println("This sets the UPPER BOUND on the importance of emission.")
println("If sigma*T^4 << I0, emission is negligible everywhere.\n")

I0 = 50e3  # 50 kW/m^2 external flux
@printf("External flux I0 = %.1f kW/m^2\n\n", I0 / 1e3)
@printf("%-10s  %15s  %15s  %15s\n", "T [K]", "sigma*T^4 [kW/m^2]", "sigma*T^4/I0", "sigma*T^4/(pi*I0)")

for T in [300, 400, 500, 600, 700, 800, 900, 1000]
    E = sigma * T^4
    @printf("%-10d  %15.3f  %15.4f  %15.5f\n", T, E / 1e3, E / I0, E / (pi * I0))
end

println("\nKey insight: sigma*T^4/(pi*I0) is the emission/absorption ratio at the SURFACE (z=0).")
println("At 500 K, the surface ratio is 0.023 -- emission is 2.3% of absorption.")
println("At 700 K, the surface ratio is 0.087 -- emission is 8.7% of absorption.")
println("At 800 K, the surface ratio is 0.148 -- emission is 15% of absorption.")
println("But this is the WORST case (uniform T). With a realistic temperature gradient")
println("(hot only near the surface), the integrated error is much smaller (see Section 4).")

# ============================================================================
# 3. Depth-dependent ratio for typical PMMA parameters
# ============================================================================
println("\n" * "=" ^ 80)
println("--- Section 2: Ratio vs depth for different absorption coefficients ---\n")

# PMMA parameters: alpha_mass ~ 1.47 m^2/kg, rho ~ 1180 kg/m^3
# alpha_vol = alpha_mass * rho ~ 1735 m^-1  (but varies; use range 100-2000)

# Penetration depths
println("Optical penetration depth 1/alpha:")
for alpha in [100, 500, 1000, 1500, 2000]
    @printf("  alpha = %5d m^-1  =>  1/alpha = %.2f mm\n", alpha, 1.0 / alpha * 1e3)
end

println("\n--- Emission/absorption ratio at specific depths ---")
println("    (T assumed uniform for worst-case analysis)\n")

for T in [500, 700, 800]
    @printf("\nT = %d K (sigma*T^4 = %.2f kW/m^2):\n", T, sigma * T^4 / 1e3)
    @printf("%-12s", "alpha [m^-1]")
    for tau in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
        @printf("  tau=%-4.1f", tau)
    end
    println()
    for alpha in [500, 1000, 1500]
        @printf("%-12d", alpha)
        for tau in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
            z = tau / alpha
            r = emission_absorption_ratio(T, I0, alpha, z)
            @printf("  %8.4f", r)
        end
        println()
    end
end

println("\nPhysical interpretation:")
println("  tau = alpha*z is the optical depth from the surface.")
println("  At tau=3, the beam has decayed by exp(-3) = 0.05 (95% absorbed).")
println("  The emission/absorption ratio grows as exp(tau) with depth.")

# ============================================================================
# 4. Depth at which emission = absorption (ratio = 1)
# ============================================================================
println("\n" * "=" ^ 80)
println("--- Section 3: Optical depth where emission equals absorption ---\n")

println("The ratio = 1 when: sigma*T^4/(pi*I0) * exp(tau) = 1")
println("  =>  tau_crit = ln(pi*I0 / (sigma*T^4))\n")

@printf("%-10s  %12s  %18s  %20s\n",
    "T [K]", "tau_crit", "z_crit (alpha=1500)",
    "Fraction I absorbed")

for T in [300, 500, 600, 700, 800, 900, 1000]
    ratio_surface = sigma * T^4 / (pi * I0)
    if ratio_surface < 1.0
        tau_crit = log(1.0 / ratio_surface)
        z_crit = tau_crit / 1500.0
        frac_absorbed = 1.0 - exp(-tau_crit)
        @printf("%-10d  %12.2f  %15.3f mm  %18.1f%%\n",
            T, tau_crit, z_crit * 1e3, frac_absorbed * 100)
    else
        @printf("%-10d  %12s  %18s  %20s\n", T, "0 (always)", "surface", "0%")
    end
end

println("\nKey result: At 700K, emission exceeds absorption at tau=2.4 (91% absorbed).")
println("At 800K, emission exceeds absorption at tau=1.9 (85% absorbed).")
println("The emission term is small where Beer-Lambert deposits most energy (near surface).")
println("It only becomes dominant deep in the material where the beam is nearly exhausted.")

# ============================================================================
# 5. Integrated error: total energy deposited with vs without emission
# ============================================================================
println("\n" * "=" ^ 80)
println("--- Section 4: Integrated volumetric source with and without emission ---\n")

println("We solve the full 1D RTE numerically:")
println("  dI/dz = -alpha*I + alpha*sigma*T^4/pi")
println("and compare q'''(z) = alpha*I(z) to Beer-Lambert q'''_BL(z) = alpha*I_BL(z)\n")

"""
Solve the 1D RTE for collimated beam in absorbing-emitting medium with given T(z).
Returns (z_array, I_RTE, I_BL, q_RTE, q_BL) on the discrete grid.
"""
function solve_1d_rte(; I0, alpha, L, T_func, N=1000)
    dz = L / N
    z = range(0, L, length=N+1) |> collect

    # Beer-Lambert (absorption only)
    I_BL = zeros(N + 1)
    I_BL[end] = I0  # radiation enters at z=L (surface)

    # Full RTE
    I_RTE = zeros(N + 1)
    I_RTE[end] = I0

    # March from surface (z=L) into depth (z=0)
    for i in N:-1:1
        z_mid = 0.5 * (z[i] + z[i+1])
        T_local = T_func(z_mid)
        Ib = sigma * T_local^4 / pi

        # Beer-Lambert: dI/dz = -alpha*I  (note: we march in -z direction)
        I_BL[i] = I_BL[i+1] * exp(-alpha * dz)

        # Full RTE: dI/dz = -alpha*I + alpha*Ib
        # Analytical solution for constant alpha, Ib over dz:
        # I(z-dz) = I(z)*exp(-alpha*dz) + Ib*(1 - exp(-alpha*dz))
        tau = alpha * dz
        I_RTE[i] = I_RTE[i+1] * exp(-tau) + Ib * (1.0 - exp(-tau))
    end

    # Volumetric absorption rate = alpha * I(z)
    q_BL = alpha .* I_BL
    q_RTE = alpha .* I_RTE

    # The NET volumetric source for RTE is alpha*I - alpha*Ib = alpha*(I - Ib)
    # For energy equation: q_net = alpha*(I - sigma*T^4/pi)
    q_net_RTE = zeros(N + 1)
    for i in 1:N+1
        T_local = T_func(z[i])
        Ib = sigma * T_local^4 / pi
        q_net_RTE[i] = alpha * (I_RTE[i] - Ib)
    end

    return (z=z, I_RTE=I_RTE, I_BL=I_BL, q_BL=q_BL, q_RTE=q_RTE, q_net_RTE=q_net_RTE)
end

# Trapezoidal integration
function trapz(z, f)
    s = 0.0
    for i in 1:length(z)-1
        s += 0.5 * (f[i] + f[i+1]) * (z[i+1] - z[i])
    end
    return s
end

# Case 1: Uniform temperature
println("Case 1: Uniform temperature profile")
println("=" ^ 50)
alpha_ref = 1500.0  # m^-1 (typical PMMA)
L = 5.0 / alpha_ref  # 5 optical depths (covers 99.3% absorption)

@printf("  alpha = %d m^-1, L = %.2f mm (5 optical depths)\n\n", alpha_ref, L * 1e3)
@printf("%-8s  %14s  %14s  %14s  %10s\n",
    "T [K]", "Q_BL [kW/m^2]", "Q_RTE_net [kW/m^2]",
    "Difference", "Rel. Error")

for T_uniform in [300, 500, 700, 800, 1000]
    T_func = z -> Float64(T_uniform)
    sol = solve_1d_rte(I0=I0, alpha=alpha_ref, L=L, T_func=T_func, N=2000)

    Q_BL = trapz(sol.z, sol.q_BL) / 1e3
    Q_RTE_net = trapz(sol.z, sol.q_net_RTE) / 1e3
    diff = Q_BL - Q_RTE_net
    rel_err = diff / Q_BL * 100

    @printf("%-8d  %14.3f  %14.3f  %14.3f  %9.2f%%\n",
        T_uniform, Q_BL, Q_RTE_net, diff, rel_err)
end

# Case 2: Linear temperature gradient (more realistic)
println("\nCase 2: Linear temperature profile T(z) = T_surf at surface, 300K at back")
println("=" ^ 50)
@printf("  alpha = %d m^-1, L = %.2f mm\n\n", alpha_ref, L * 1e3)
@printf("%-12s  %14s  %14s  %14s  %10s\n",
    "T_surf [K]", "Q_BL [kW/m^2]", "Q_RTE_net [kW/m^2]",
    "Difference", "Rel. Error")

for T_surf in [500, 600, 700, 800, 1000]
    T_func = z -> 300.0 + (T_surf - 300.0) * (z / L)  # T increases toward surface (z=L)
    sol = solve_1d_rte(I0=I0, alpha=alpha_ref, L=L, T_func=T_func, N=2000)

    Q_BL = trapz(sol.z, sol.q_BL) / 1e3
    Q_RTE_net = trapz(sol.z, sol.q_net_RTE) / 1e3
    diff = Q_BL - Q_RTE_net
    rel_err = diff / Q_BL * 100

    @printf("%-12d  %14.3f  %14.3f  %14.3f  %9.2f%%\n",
        T_surf, Q_BL, Q_RTE_net, diff, rel_err)
end

# ============================================================================
# 6. Analysis of ThermaKin's Eq. 12 formulation
# ============================================================================
println("\n" * "=" ^ 80)
println("--- Section 5: Analysis of ThermaKin's re-radiation (Eq. 12) ---\n")

println("""ThermaKin's Eq. 12:
  dI_rr/dx = (sigma*T^4 * sum(eps_i * v_i) / I_ex^0) * dI_ex/dx

This can be rewritten using Eq. 11 (dI_ex/dx = -I_ex * sum(alpha_i * xi_i)):
  dI_rr/dx = -(sigma*T^4 * sum(eps_i * v_i) / I_ex^0) * I_ex * sum(alpha_i * xi_i)

For a single component with eps = 1 (opaque material), this simplifies to:
  dI_rr/dx = -(sigma*T^4 / I_ex^0) * alpha * I_ex

Compare to the standard RTE emission term:
  (dI/dz)_emission = +alpha * sigma * T^4 / pi

KEY DIFFERENCES:
1. TK scales by I_ex/I_ex^0 (local-to-surface intensity ratio), so emission
   decays exponentially into depth just like absorption. The RTE emission term
   does NOT decay -- it depends only on local T.

2. TK divides by I_ex^0 (surface flux), not by pi. This is dimensionally
   consistent but physically different from the RTE.

3. In the RTE, emission is ISOTROPIC (hence 1/pi for hemispherical intensity).
   TK's formulation makes emission DIRECTIONAL, proportional to the local
   Beer-Lambert absorption rate. This is not standard radiation physics.

Physical interpretation of TK's approach:
  TK computes the re-radiation flux as a CORRECTION to Beer-Lambert that
  reduces the net absorbed energy. The correction is proportional to:
  - sigma*T^4 (blackbody emission, correct)
  - the LOCAL absorption rate dI_ex/dx (non-standard)

  This makes the re-radiation vanish deep in the material where the beam
  has been fully absorbed, which is physically wrong -- a hot medium emits
  regardless of whether there is an external beam present.

  However, the deep emission in the RTE is handled separately by TK through
  the "radiation diffusion approximation" (k ~ T^3 conductivity term), so
  TK's Eq. 12 is specifically correcting the Beer-Lambert absorption near
  the surface where the beam is strong.
""")

# Compute TK's re-radiation vs RTE emission for comparison
println("Numerical comparison: TK Eq. 12 vs RTE emission term")
println("=" ^ 50)
println("(For single component, eps=1, uniform T=700K, alpha=1500 m^-1)\n")

T_comp = 700.0
alpha_comp = 1500.0

@printf("%-12s  %15s  %15s  %12s\n",
    "tau=alpha*z", "RTE emission", "TK re-rad rate", "TK/RTE ratio")
@printf("%-12s  %15s  %15s  %12s\n",
    "", "[MW/m^3]", "[MW/m^3]", "")

for tau in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
    # RTE emission term: alpha * sigma * T^4 / pi  (constant, independent of depth)
    rte_emission = alpha_comp * sigma * T_comp^4 / pi

    # TK re-radiation rate (volumetric):
    # dI_rr/dx = -(sigma*T^4 / I0) * alpha * I0 * exp(-tau)
    # |dI_rr/dx| = sigma*T^4 * alpha * exp(-tau)
    tk_rerad = sigma * T_comp^4 * alpha_comp * exp(-tau)

    @printf("%-12.1f  %15.4f  %15.4f  %12.3f\n",
        tau, rte_emission / 1e6, tk_rerad / 1e6, tk_rerad / rte_emission)
end

println("\nNote: TK's re-radiation is pi times the RTE emission at the surface (tau=0),")
println("then decays exponentially. The RTE emission is constant with depth.")
println("TK OVERESTIMATES emission near the surface and UNDERESTIMATES it at depth.")

# ============================================================================
# 7. Comparison: Beer-Lambert vs P1 vs Beer-Lambert+emission
# ============================================================================
println("\n" * "=" ^ 80)
println("--- Section 6: Beer-Lambert vs P1 diffusion approximation ---\n")

println("""The P1 model solves:
  div(D_rad * grad(G)) = alpha*(G - 4*sigma*T^4)

where G is the incident radiation (zeroth moment of intensity) and D_rad = 1/(3*alpha).

This is the COMPLETE diffusion approximation to the RTE. It includes:
  - Absorption: alpha*G  (energy removed from radiation field)
  - Emission:   4*alpha*sigma*T^4  (energy added to radiation field)
  - Diffusion:  radiation "conducts" through the medium

The net volumetric source for the energy equation is:
  q_rad = alpha*(G - 4*sigma*T^4)

When is P1 vs Beer-Lambert appropriate?
  - P1 is valid for optically THICK media (tau >> 1) where radiation
    undergoes many absorption/re-emission events.
  - Beer-Lambert is valid for COLLIMATED beams in media where the beam
    dominates over thermal emission (I_beam >> sigma*T^4/pi).

For PYROLYSIS of semi-transparent polymers (PMMA):
  - External flux: 50 kW/m^2 (collimated or near-collimated)
  - Material is semi-transparent: alpha ~ 1000-2000 m^-1
  - Thickness: ~few mm -> optical depth tau ~ 3-10
  - Temperature: 300-700K in the absorption zone

The Beer-Lambert model is the natural choice for the EXTERNAL beam because:
  1. The beam is collimated (not diffuse like P1 assumes)
  2. In the absorption zone near the surface, the beam dominates over emission
  3. The P1 model would smear the beam (diffusion cannot represent collimated transport)

The P1 model is the natural choice for INTERNAL radiation (re-emission) because:
  1. Emitted radiation IS isotropic and diffuse
  2. In optically thick media, the diffusion approximation is excellent
  3. This gives the effective "radiation conductivity" k_rad ~ 16*sigma*T^3/(3*alpha)
""")

# Compute the P1 "radiation conductivity" for comparison
println("Effective radiation conductivity k_rad = 16*sigma*T^3 / (3*alpha):")
println("=" ^ 50)
@printf("%-10s  %-15s", "T [K]", "alpha [m^-1]")
@printf("  %15s  %15s\n", "k_rad [W/m/K]", "k_cond (PMMA)")

k_pmma = 0.21  # W/(m K), typical PMMA conductivity
for T in [400, 600, 800, 1000]
    for alpha in [500, 1500]
        k_rad = 16.0 * sigma * T^3 / (3.0 * alpha)
        @printf("%-10d  %-15d  %15.4f  %15.3f\n", T, alpha, k_rad, k_pmma)
    end
end

println("\nAt T=600K, k_rad = 0.04-0.13 W/m/K: 21-62% of k_cond for PMMA.")
println("At T=800K, k_rad = 0.10-0.31 W/m/K: 49-147% of k_cond for PMMA.")
println("Radiation conductivity is significant at pyrolysis temperatures!")
println("This is the dominant internal radiation transport mechanism,")
println("handled properly by the P1 model (or by k ~ T^3 conductivity).")

# ============================================================================
# 8. Final quantitative summary
# ============================================================================
println("\n" * "=" ^ 80)
println("--- CONCLUSIONS ---\n")

println("""
1. MAGNITUDE OF THE MISSING EMISSION TERM IN BEER-LAMBERT

   For typical pyrolysis conditions (I0=50 kW/m^2, alpha=1500 m^-1):

   Uniform T (worst case -- entire material at one temperature):
   - At T=700K: emission/absorption ratio = 8.7% at surface, 64% at tau=2.
     Integrated error: Beer-Lambert overestimates net absorption by 8.7%.
   - At T=800K: ratio = 14.8% at surface, 109% at tau=2.
     Integrated error: 14.8%.

   Realistic temperature gradient (T_surf at surface, 300K at back):
   - At T_surf=700K: integrated overestimate is only 0.9%.
   - At T_surf=800K: integrated overestimate is only 1.3%.
   - At T_surf=1000K: integrated overestimate is 2.2%.

   The gradient case is far smaller because only a thin layer near the
   surface is hot, and the beam is strongest there (emission/absorption
   ratio is smallest near the surface).

2. IS THE MISSING EMISSION A REAL DEFICIENCY?

   For REALISTIC temperature gradients, it is a SMALL deficiency (1-2%)
   that is ACCEPTABLE for the following reasons:

   a) Beer-Lambert is inherently a simplified model for collimated beam transport.
      It neglects scattering, diffuse reflection, and multi-directional effects.
      The emission error is smaller than these other approximations.

   b) The emission that IS physically important -- thermal radiation transport
      between hot and cold regions of the material -- is properly handled by
      either:
      - The P1 model (which solves the full diffusion approximation including
        emission), or
      - The radiation conductivity k_rad ~ T^3 term in thermal conductivity
        (as ThermaKin does with its "radiation diffusion approximation")

   c) Beer-Lambert and the emission correction operate in different regimes:
      Beer-Lambert handles the collimated external beam (forward-peaked).
      Emission is isotropic and better handled by a diffusion model.

3. IS THERMAKIN'S Eq. 12 PHYSICALLY CORRECT?

   TK's re-radiation formulation (Eq. 12) is NON-STANDARD and has issues:

   a) It scales emission by dI_ex/dx, making it proportional to the local
      Beer-Lambert absorption rate. In the RTE, emission depends only on
      local temperature, NOT on the beam intensity.

   b) At the surface (tau=0), TK's emission is pi times the RTE value.
      At depth, TK's emission decays exponentially while the RTE emission
      stays constant (for uniform T).

   c) However, TK's approach can be interpreted as a PRAGMATIC correction:
      it reduces the net absorbed energy in the Beer-Lambert zone by an
      amount that is roughly correct in an integrated sense, even if the
      spatial distribution differs from the exact RTE.

   d) TK handles deep-material radiation transport separately via the
      radiation diffusion approximation (k ~ T^3), so Eq. 12 is specifically
      addressing surface-zone emission only.

4. RECOMMENDATION FOR THIS CODEBASE

   The current implementation (Beer-Lambert without emission) is ACCEPTABLE
   because:
   - The error is 1-2% for realistic temperature gradients
   - For uniform high T (worst case), errors reach 9-15%, but this scenario
     does not arise in practice during pyrolysis (only the surface is hot)
   - The P1 model is available for cases requiring full emission physics
   - Adding emission to Beer-Lambert would create a hybrid model that is
     neither the clean Beer-Lambert approximation nor the complete RTE

   If emission correction is ever desired, the CORRECT approach is:
   q_net(z) = alpha * I(z) - alpha * sigma * T(z)^4 / pi    (from the RTE)
   NOT ThermaKin's Eq. 12 formulation.
""")
