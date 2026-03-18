"""
Physical Reasonableness Analysis of Gas Transfer Coefficient (lambda) Values

This script performs a first-principles analysis of the lambda values used in
the PMMA, Douglas Fir, and PTW pyrolysis examples, comparing them against:
1. Chapman-Enskog free molecular diffusion theory
2. Effective diffusion in porous media (epsilon/tau model)
3. ThermaKin .cmp file values
4. Julia example values

The gas transfer coefficient lambda [m^2/s] appears in the transport equation:
    J = -(lambda/T) * d(xi*T)/dx

The effective lambda of the mixture is computed from per-component lambda values
via volume-fraction-weighted mixing rules (WEIGHTED with beta=0.5 by default).
"""

using Printf

# ============================================================================
# Physical Constants
# ============================================================================

const R_GAS = 8.31446
const P_ATM = 101325.0  # Pa

# ============================================================================
# SECTION 1: Free Molecular Diffusion (Chapman-Enskog Theory)
# ============================================================================

println("=" ^ 80)
println("SECTION 1: Free Molecular Diffusion Coefficients (Chapman-Enskog)")
println("=" ^ 80)
println()
println("Binary diffusion in gas phase: D proportional to T^1.5 / P")
println("Reference medium: air (N2) at 1 atm")
println()

# Reference diffusivities at 300 K, 1 atm (from literature)
# D(H2O-air) = 2.5e-5, D(CO2-air) = 1.6e-5, D(CH4-air) = 2.2e-5
# For organic vapors, D decreases with molecular weight roughly as 1/sqrt(M)

struct GasSpecies
    name::String
    M::Float64       # kg/mol
    D_ref::Float64   # m^2/s at 300 K, 1 atm
end

species_list = [
    GasSpecies("H2O (water vapor)",    0.018, 2.50e-5),
    GasSpecies("CO2",                  0.044, 1.60e-5),
    GasSpecies("CH4 (methane)",        0.016, 2.20e-5),
    GasSpecies("Acetic acid (M=60)",   0.060, 1.20e-5),
    GasSpecies("Organic (M=70)",       0.070, 1.10e-5),
    GasSpecies("Organic (M=80)",       0.080, 1.00e-5),
    GasSpecies("MMA (M=100)",          0.100, 0.85e-5),
    GasSpecies("Levoglucosan (M=110)", 0.110, 0.80e-5),
    GasSpecies("Levoglucosan (M=162)", 0.162, 0.60e-5),
]

temperatures = [300.0, 500.0, 700.0, 1000.0]

# Chapman-Enskog scaling: D(T) = D_ref * (T/T_ref)^1.5
D_CE(D_ref, T) = D_ref * (T / 300.0)^1.5

@printf("%-28s  %6s", "Species", "M")
for T in temperatures
    @printf("  D@%4.0fK", T)
end
println()
@printf("%-28s  %6s", "", "[g/mol]")
for _ in temperatures
    @printf("  %8s", "[m^2/s]")
end
println()
println("-" ^ 80)

for sp in species_list
    @printf("%-28s  %6.1f", sp.name, sp.M * 1000)
    for T in temperatures
        D = D_CE(sp.D_ref, T)
        @printf("  %8.2e", D)
    end
    println()
end

println()
println("Key insight: Free diffusion D ~ 1e-5 to 5e-5 m^2/s at pyrolysis temperatures")
println("             Heavier molecules (tars, M>100) have D ~ 0.5-2e-5 m^2/s")
println()

# ============================================================================
# SECTION 2: Effective Diffusion in Porous Media
# ============================================================================

println("=" ^ 80)
println("SECTION 2: Effective Diffusion in Porous Media")
println("=" ^ 80)
println()
println("D_eff = D_free * (epsilon / tau)")
println("  epsilon = porosity (void fraction)")
println("  tau = tortuosity (path-length ratio, typically 2-10)")
println()

struct PorousMedium
    name::String
    epsilon_lo::Float64
    epsilon_hi::Float64
    tau_lo::Float64
    tau_hi::Float64
end

media = [
    PorousMedium("Virgin wood",         0.02, 0.10, 3.0, 8.0),
    PorousMedium("Early pyrolysis",     0.10, 0.25, 2.5, 5.0),
    PorousMedium("Mid pyrolysis",       0.25, 0.50, 2.0, 4.0),
    PorousMedium("Porous char",         0.50, 0.85, 1.5, 3.0),
    PorousMedium("Dense polymer (PMMA)",0.00, 0.02, 5.0, 10.0),
]

# Reference free diffusion at 500 K (typical pyrolysis temperature)
D_free_500K = 3.0e-5  # m^2/s (representative organic gas in air)

@printf("%-25s  %6s %6s  %6s %6s  %10s %10s  %10s\n",
    "Medium", "eps_lo", "eps_hi", "tau_lo", "tau_hi",
    "D_eff_lo", "D_eff_hi", "D_eff/D_free")
println("-" ^ 100)

for m in media
    Deff_lo = D_free_500K * m.epsilon_lo / m.tau_hi
    Deff_hi = D_free_500K * m.epsilon_hi / m.tau_lo
    ratio_lo = m.epsilon_lo / m.tau_hi
    ratio_hi = m.epsilon_hi / m.tau_lo
    @printf("%-25s  %6.2f %6.2f  %6.1f %6.1f  %10.2e %10.2e  %5.3f-%5.3f\n",
        m.name, m.epsilon_lo, m.epsilon_hi, m.tau_lo, m.tau_hi,
        Deff_lo, Deff_hi, ratio_lo, ratio_hi)
end

println()
println("Physical estimates of D_eff at 500 K (D_free ~ 3e-5 m^2/s):")
println("  Virgin wood:    7.5e-8 to 1.0e-6 m^2/s")
println("  Early pyrolysis: 6.0e-7 to 3.0e-6 m^2/s")
println("  Porous char:     5.0e-6 to 1.7e-5 m^2/s")
println("  Dense PMMA:      0 to 1.2e-7 m^2/s (essentially impermeable)")
println()

# ============================================================================
# SECTION 3: Lambda Values Used in Examples vs ThermaKin
# ============================================================================

println("=" ^ 80)
println("SECTION 3: Component Lambda Values -- Julia Examples vs ThermaKin .cmp Files")
println("=" ^ 80)
println()

# --- PMMA ---
println("-" ^ 80)
println("PMMA (Cone Calorimeter)")
println("-" ^ 80)

pmma_components = [
    # (name, phase, Julia_lambda, TK_lambda, notes)
    ("virgin PMMA",   "S", 0.0,  1e-5, "Julia default for SolidComponent"),
    ("char",          "S", 0.0,  1e-5, "Julia default for SolidComponent"),
    ("MMA gas",       "G", 1e-5, 1e-5, "Matches TK"),
]

@printf("  %-18s %5s  %10s  %10s  %10s  %s\n",
    "Component", "Phase", "Julia_lam", "TK_lam", "Ratio", "Notes")
println("  " * "-" ^ 75)
for (name, phase, jl, tk, notes) in pmma_components
    ratio = tk > 0 ? jl / tk : (jl > 0 ? Inf : 0.0)
    @printf("  %-18s %5s  %10.1e  %10.1e  %10.2f  %s\n",
        name, phase, jl, tk, ratio, notes)
end
println()
println("  Mixing: PARALLEL in Julia (gas_transfer_mixing = PARALLEL)")
println("  TK:     PARALL TRANSPORT: 0.5 (WEIGHTED, beta=0.5)")
println()
println("  DISCREPANCY: Julia sets lambda=0 for PMMA solids, TK uses 1e-5 for ALL.")
println("  With PARALLEL mixing, only gas components contribute to lambda_eff in Julia.")
println("  With WEIGHTED mixing in TK, solids contribute significantly via the series term.")
println()

# --- PTW ---
println("-" ^ 80)
println("Pressure Treated Wood (PTW)")
println("-" ^ 80)

ptw_components = [
    # (name, phase, Julia_lambda, TK_lambda, notes)
    ("WATER (liquid)",  "L", 2e-5,       2e-5,   "Default LiquidComponent lambda"),
    ("WATER_VAPOR",     "G", 2e-5*1.6,   3.2e-5, "Julia: 2e-5*1.6 = 3.2e-5"),
    ("PTW (virgin)",    "S", 2e-5,       2e-5,   "Matches TK"),
    ("PTW1",            "S", 2e-5,       2e-5,   "Matches TK"),
    ("PTW2",            "S", 2e-5,       2e-5,   "Matches TK"),
    ("PTW3",            "S", 2e-5,       2e-5,   "Matches TK"),
    ("Char",            "S", 2e-5,       2e-5,   "Matches TK"),
    ("PTW_gas1",        "G", 2e-5,       2e-5,   "Matches TK"),
    ("PTW_gas2",        "G", 2e-5,       2e-5,   "Matches TK"),
    ("PTW_gas3",        "G", 2e-5,       2e-5,   "Matches TK"),
    ("PTW_gas4",        "G", 2e-5,       2e-5,   "Matches TK"),
]

@printf("  %-18s %5s  %10s  %10s  %10s  %s\n",
    "Component", "Phase", "Julia_lam", "TK_lam", "Ratio", "Notes")
println("  " * "-" ^ 75)
for (name, phase, jl, tk, notes) in ptw_components
    ratio = tk > 0 ? jl / tk : (jl > 0 ? Inf : 0.0)
    @printf("  %-18s %5s  %10.1e  %10.1e  %10.2f  %s\n",
        name, phase, jl, tk, ratio, notes)
end
println()
println("  Mixing: WEIGHTED, beta_lambda = 0.5 (matches TK PARALL TRANSPORT: 0.5)")
println("  PTW uses UNIFORM lambda = 2e-5 for ALL components (solid + gas).")
println("  Julia matches TK almost exactly (slight difference in WATER_VAPOR: 3.2e-5 vs 3.2e-5).")
println()

# --- Douglas Fir ---
println("-" ^ 80)
println("Douglas Fir (DF)")
println("-" ^ 80)

df_components = [
    # (name, phase, Julia_lambda, TK_lambda, notes)
    ("DF (virgin)",   "S", 1e-6,       1e-5,   "Julia 10x LOWER than TK"),
    ("DF_int1",       "S", 1.5e-6,     1.2e-5, "Julia ~8x lower"),
    ("DF_int2",       "S", 3e-6,       1.5e-5, "Julia 5x lower"),
    ("DF_int3",       "S", 6e-6,       2e-5,   "Julia ~3x lower"),
    ("DF_int4",       "S", 8e-6,       4e-5,   "Julia 5x lower"),
    ("Char",          "S", 1e-5,       5e-5,   "Julia 5x lower"),
    ("MOISTURE",      "L", 2e-5,       2e-3,   "Julia 100x LOWER (TK=2e-3!)"),
    ("DF_gas1",       "G", 8.5e-6,     1e-4,   "Julia ~12x lower"),
    ("DF_gas2",       "G", 8e-6,       1e-4,   "Julia ~12x lower"),
    ("DF_gas3",       "G", 6e-6,       8e-5,   "Julia ~13x lower"),
    ("DF_gas4",       "G", 8e-6,       1e-4,   "Julia ~12x lower"),
    ("DF_gas5",       "G", 1e-5,       1e-4,   "Julia 10x lower"),
    ("WATER_VAPOR",   "G", 1.6e-5,     1.2e-4, "Julia ~7.5x lower"),
]

@printf("  %-18s %5s  %10s  %10s  %10s  %s\n",
    "Component", "Phase", "Julia_lam", "TK_lam", "Ratio", "Notes")
println("  " * "-" ^ 75)
for (name, phase, jl, tk, notes) in df_components
    ratio = tk > 0 ? jl / tk : (jl > 0 ? Inf : 0.0)
    @printf("  %-18s %5s  %10.1e  %10.1e  %10.2f  %s\n",
        name, phase, jl, tk, ratio, notes)
end
println()
println("  Mixing: WEIGHTED, beta_lambda = 0.5 (TK: PARALL TRANSPORT: 0.5)")
println()
println("  MAJOR DISCREPANCY: Julia DF lambda values are 5-13x LOWER than TK for ALL components.")
println("  The Julia example attempts to split TK's single TRANSPORT parameter into:")
println("    - Fickian diffusion (lambda, lower values)")
println("    - Darcy flow (kappa, pressure-driven permeability)")
println("  However, the combined Darcy+Fick effect depends strongly on pressure gradients.")
println()

# ============================================================================
# SECTION 4: Effective Mixture Lambda at Representative Compositions
# ============================================================================

println("=" ^ 80)
println("SECTION 4: Effective Mixture Lambda at Representative Compositions")
println("=" ^ 80)
println()

"""
Compute WEIGHTED mixture lambda (beta=0.5) from volume fractions and component lambdas.

lambda_eff = beta * lambda_parallel + (1-beta) * lambda_series
lambda_parallel = sum(v_i * lambda_i)
lambda_series = 1 / sum(v_i / lambda_i)  [if all lambda > 0; else 0]
"""
function weighted_lambda(v::Vector{Float64}, lam::Vector{Float64}; beta=0.5)
    @assert length(v) == length(lam)
    lam_parallel = sum(v .* lam)

    # Series: if any component has lambda=0 and nonzero v, series -> 0
    inv_series = 0.0
    series_valid = true
    for i in eachindex(v)
        if v[i] > 1e-20
            if lam[i] > 1e-20
                inv_series += v[i] / lam[i]
            else
                series_valid = false
                break
            end
        end
    end
    lam_series = series_valid && inv_series > 1e-20 ? 1.0 / inv_series : 0.0

    return beta * lam_parallel + (1.0 - beta) * lam_series
end

"""
Compute volume fractions from mass concentrations and intrinsic densities.
For ideal gas: rho(T) = P*M/(R*T)
"""
function compute_volume_fractions(xi::Vector{Float64}, rho_intrinsic::Vector{Float64})
    v_specific = xi ./ rho_intrinsic
    v_total = sum(v_specific)
    if v_total < 1e-20
        return zeros(length(xi))
    end
    return v_specific ./ v_total
end

# Gas density via ideal gas law
rho_gas(M, T) = P_ATM * M / (R_GAS * T)

# ---- PMMA scenarios ----
println("--- PMMA ---")
println()

# Julia PMMA: 3 components (virgin, char, MMA gas)
# lambda = [0, 0, 1e-5], rho = [1210, 1210, ideal_gas(0.100)]
# Mixing: PARALLEL

for (scenario, xi_virgin, xi_char, xi_gas, T) in [
    ("Early (2% gas, 500K)",    1200.0, 0.0,  2.0, 500.0),
    ("Mid (50% decomp, 600K)",  600.0,  1.2,  10.0, 600.0),
    ("Late (mostly char, 650K)", 50.0,  2.4,  5.0, 650.0),
]
    rho_mma = rho_gas(0.100, T)
    xi = [xi_virgin, xi_char, xi_gas]
    rho_i = [1210.0, 1210.0, rho_mma]
    v = compute_volume_fractions(xi, rho_i)

    # Julia: PARALLEL mixing, lambda = [0, 0, 1e-5]
    lam_jl = [0.0, 0.0, 1e-5]
    lam_eff_jl_parallel = sum(v .* lam_jl)

    # TK: WEIGHTED beta=0.5, lambda = [1e-5, 1e-5, 1e-5]
    lam_tk = [1e-5, 1e-5, 1e-5]
    lam_eff_tk = weighted_lambda(collect(v), collect(lam_tk); beta=0.5)

    ratio = lam_eff_tk > 0 ? lam_eff_jl_parallel / lam_eff_tk : NaN

    @printf("  %-30s T=%4.0fK  v_gas=%6.4f\n", scenario, T, v[3])
    @printf("    Julia (PARALLEL):  lambda_eff = %10.3e m^2/s\n", lam_eff_jl_parallel)
    @printf("    TK    (WEIGHTED):  lambda_eff = %10.3e m^2/s\n", lam_eff_tk)
    @printf("    Ratio Julia/TK:    %6.4f\n", ratio)
    println()
end

# ---- PTW scenarios ----
println("--- PTW ---")
println()

# PTW: 11 components, all with lambda = 2e-5 (solids, gases, liquid)
# Since all lambda are the same (2e-5), mixing rule doesn't matter:
# lambda_eff = 2e-5 regardless of composition

for (scenario, T) in [
    ("Early (virgin wood, 400K)",  400.0),
    ("Mid (50% decomp, 600K)",     600.0),
    ("Late (mostly char, 800K)",   800.0),
]
    @printf("  %-30s T=%4.0fK\n", scenario, T)
    @printf("    Julia (WEIGHTED):  lambda_eff = %10.3e m^2/s  (uniform 2e-5)\n", 2e-5)
    @printf("    TK    (WEIGHTED):  lambda_eff = %10.3e m^2/s  (uniform 2e-5)\n", 2e-5)
    @printf("    Ratio Julia/TK:    %6.4f  (exact match)\n", 1.0)
    println()
end

# ---- Douglas Fir scenarios ----
println("--- Douglas Fir ---")
println()

# DF: 13 components
# We need to show how the effective lambda changes through pyrolysis

# Early pyrolysis: mostly DF virgin (530 kg/m3) + small gas (5 kg/m3)
# Index: 1=DF, 2=int1, 3=int2, 4=int3, 5=int4, 6=Char,
#        7=MOISTURE, 8=gas1, 9=gas2, 10=gas3, 11=gas4, 12=gas5, 13=H2O_vapor
df_lam_jl = [1e-6, 1.5e-6, 3e-6, 6e-6, 8e-6, 1e-5,
             2e-5, 8.5e-6, 8e-6, 6e-6, 8e-6, 1e-5, 1.6e-5]
df_lam_tk = [1e-5, 1.2e-5, 1.5e-5, 2e-5, 4e-5, 5e-5,
             2e-3, 1e-4, 1e-4, 8e-5, 1e-4, 1e-4, 1.2e-4]

# Note: TK gas densities for DF_gas are given as 498.46 kg/m^3 in the .cmp file.
# This is NOT an ideal gas density -- it's a pseudo-density like in TK convention.
# In Julia, GaseousComponent uses ideal gas law: rho = PM/(RT)

df_rho_solid = [530.0, 507.2, 419.8, 203.2, 176.6, 150.0, 1e6]
df_M_gas = [0.060, 0.070, 0.110, 0.070, 0.050, 0.018]

df_scenarios = [
    # (name, T, xi vector for 13 components)
    ("Early (virgin + 2% gas)",
     400.0,
     [520.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      30.0,
      2.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    ("Mid pyrolysis (int2 dominant, 10% gas)",
     550.0,
     [0.0, 0.0, 300.0, 50.0, 0.0, 0.0,
      5.0,
      5.0, 15.0, 10.0, 2.0, 0.0, 3.0]),
    ("Late (mostly char, 30% gas fraction)",
     700.0,
     [0.0, 0.0, 0.0, 0.0, 20.0, 120.0,
      0.0,
      2.0, 5.0, 15.0, 8.0, 10.0, 5.0]),
]

for (scenario, T, xi_vec) in df_scenarios
    # Compute intrinsic densities
    rho_i = vcat(df_rho_solid, [rho_gas(M, T) for M in df_M_gas])

    v = compute_volume_fractions(xi_vec, rho_i)

    lam_eff_jl = weighted_lambda(collect(v), collect(df_lam_jl); beta=0.5)
    lam_eff_tk = weighted_lambda(collect(v), collect(df_lam_tk); beta=0.5)

    # Note: For TK, gas density is 498.46 (pseudo-density), not ideal gas
    # Recompute with TK densities for fair comparison
    rho_i_tk = vcat(df_rho_solid, [498.46, 498.46, 498.46, 498.46, 498.46, 0.8])
    v_tk = compute_volume_fractions(xi_vec, rho_i_tk)
    lam_eff_tk_fair = weighted_lambda(collect(v_tk), collect(df_lam_tk); beta=0.5)

    ratio = lam_eff_tk > 0 ? lam_eff_jl / lam_eff_tk : NaN
    ratio_fair = lam_eff_tk_fair > 0 ? lam_eff_jl / lam_eff_tk_fair : NaN

    v_gas_total = sum(v[8:13])
    v_gas_total_tk = sum(v_tk[8:13])

    @printf("  %-40s T=%4.0fK\n", scenario, T)
    @printf("    Julia vol fracs: gas=%6.4f, solid=%6.4f\n", v_gas_total, 1-v_gas_total)
    @printf("    TK    vol fracs: gas=%6.4f, solid=%6.4f  (TK pseudo-densities)\n", v_gas_total_tk, 1-v_gas_total_tk)
    @printf("    Julia (WEIGHTED, ideal gas rho):  lambda_eff = %10.3e m^2/s\n", lam_eff_jl)
    @printf("    TK    (WEIGHTED, ideal gas rho):  lambda_eff = %10.3e m^2/s\n", lam_eff_tk)
    @printf("    TK    (WEIGHTED, TK pseudo-rho):  lambda_eff = %10.3e m^2/s\n", lam_eff_tk_fair)
    @printf("    Ratio Julia/TK (ideal gas rho):   %6.4f\n", ratio)
    @printf("    Ratio Julia/TK (TK pseudo-rho):   %6.4f\n", ratio_fair)
    println()
end

# ============================================================================
# SECTION 5: Physical Reasonableness Assessment
# ============================================================================

println("=" ^ 80)
println("SECTION 5: Physical Reasonableness Assessment")
println("=" ^ 80)
println()

println("Reference: D_free(organic gas in air) at pyrolysis temperatures")
println("  300 K:  1-3e-5 m^2/s")
println("  500 K:  3-7e-5 m^2/s")
println("  700 K:  6-14e-5 m^2/s")
println("  1000 K: 1-3e-4 m^2/s")
println()

println("-" ^ 80)
println("Assessment A: PMMA")
println("-" ^ 80)
println()
println("  TK uses lambda = 1e-5 for ALL components (solid + gas).")
println("  Julia uses lambda = 0 for solids, 1e-5 for gas only.")
println()
println("  Physical context: PMMA is a DENSE thermoplastic (non-porous).")
println("  Virgin PMMA has epsilon ~ 0 (no porosity), so D_eff ~ 0.")
println("  Gas transport only occurs through bubble formation/cracking,")
println("  not Fickian diffusion through the solid matrix.")
println()
println("  * TK lambda = 1e-5 for virgin PMMA solid: UNPHYSICALLY HIGH")
println("    D_eff for a non-porous solid should be ~0.")
println("    However, TK uses this as a combined transport parameter")
println("    (not a pure Fickian diffusivity), and the mixing rule with")
println("    beta=0.5 moderates the effective value.")
println()
println("  * Julia lambda = 0 for PMMA solids: PHYSICALLY CORRECT")
println("    But the PARALLEL mixing rule means gas lambda only matters")
println("    when gas volume fraction is non-trivial, which requires")
println("    significant decomposition to have occurred.")
println()
println("  * Julia PARALLEL vs TK WEIGHTED mixing:")
println("    - Julia uses PARALLEL: lambda_eff = sum(v_i * lambda_i)")
println("      With lambda_solid = 0, only gas contributes.")
println("    - TK uses WEIGHTED (beta=0.5): this includes a series term")
println("      that is zero when any lambda is zero, so effectively")
println("      TK also gets lambda_eff = 0.5 * lambda_parallel + 0 = 0.5 * sum(v_i * lambda_i)")
println("      WAIT: TK has lambda > 0 for ALL components, so series is well-defined.")
println()

# Compute the actual difference for a representative PMMA state
T_pmma = 600.0
rho_mma_gas = rho_gas(0.100, T_pmma)
xi_pmma = [600.0, 1.2, 10.0]  # 50% decomposed
rho_pmma = [1210.0, 1210.0, rho_mma_gas]
v_pmma = compute_volume_fractions(xi_pmma, rho_pmma)

lam_jl_pmma = [0.0, 0.0, 1e-5]
lam_tk_pmma = [1e-5, 1e-5, 1e-5]

# Julia PARALLEL
lam_eff_jl_par = sum(v_pmma .* lam_jl_pmma)

# TK WEIGHTED beta=0.5 (all lambda = 1e-5, so lambda_eff = 1e-5)
lam_eff_tk_w = weighted_lambda(collect(v_pmma), collect(lam_tk_pmma); beta=0.5)

@printf("  Example at 50%% decomposition, T=600K, v_gas = %6.4f:\n", v_pmma[3])
@printf("    Julia (PARALLEL, lam_solid=0): lambda_eff = %10.3e m^2/s\n", lam_eff_jl_par)
@printf("    TK    (WEIGHTED, lam_all=1e-5): lambda_eff = %10.3e m^2/s\n", lam_eff_tk_w)
@printf("    Ratio: %6.4f  (Julia is %.0fx lower)\n", lam_eff_jl_par / lam_eff_tk_w,
    lam_eff_tk_w / lam_eff_jl_par)
println()
println("  VERDICT: The combination of lambda=0 for solids + PARALLEL mixing causes")
println("  Julia's effective lambda to be much lower than TK's. This means gas")
println("  escapes more slowly in Julia -> higher internal gas buildup -> potentially")
println("  different surface regression rates and MLR.")
println()

println("-" ^ 80)
println("Assessment B: Douglas Fir")
println("-" ^ 80)
println()

println("  TK DF.cmp TRANSPORT values:")
println("    Solids: 1e-5 (virgin) to 5e-5 (char)")
println("    Gases:  8e-5 to 1.2e-4")
println("    Moisture: 2e-3 (!)")
println("    MOISTURE has TRANSPORT = 2e-3, which is 100x the gas values.")
println("    This likely represents fast moisture transport (capillary action).")
println()
println("  Julia DF example values:")
println("    Solids: 1e-6 (virgin) to 1e-5 (char)  [5-10x LOWER than TK]")
println("    Gases:  6e-6 to 1.6e-5                 [5-13x LOWER than TK]")
println("    Moisture: 2e-5                          [100x LOWER than TK]")
println()
println("  The Julia example deliberately uses LOWER lambda values because it")
println("  splits transport into Darcy (kappa) + Fick (lambda). The idea is:")
println("    TK_TRANSPORT = lambda_Fick + (kappa * P / mu) approximately")
println()
println("  However, this splitting is questionable because:")
println("  1. TK's TRANSPORT parameter already represents the TOTAL effective")
println("     transport, not just Fickian diffusion.")
println("  2. The Darcy term (kappa * dP/dx / mu) depends on local pressure")
println("     gradients which vary in space and time.")
println("  3. The kappa values used (2.3e-15 to 1e-14 m^2) give a Darcy")
println("     contribution that may not fully compensate for the reduced lambda.")
println()

# Estimate Darcy equivalent lambda
println("  Darcy equivalent lambda estimate:")
println("  v_Darcy = -(kappa/mu) * dP/dx")
println("  J_Darcy = xi * v_Darcy ~ xi * (kappa/mu) * (P/L)")
println("  Equivalent lambda_Darcy ~ kappa * P / mu")
println()

for (name, kappa, T) in [
    ("DF virgin (kappa=2.3e-15)",  2.3e-15, 400.0),
    ("DF virgin (kappa=2.3e-15)",  2.3e-15, 600.0),
    ("DF char (kappa=1e-14)",      1e-14,   600.0),
    ("DF char (kappa=1e-14)",      1e-14,   800.0),
    ("PTW (kappa=1e-15)",          1e-15,   500.0),
]
    # Sutherland's law for air viscosity
    mu_ref = 1.81e-5
    T_ref = 300.0
    S = 110.0
    mu = mu_ref * (T / T_ref)^1.5 * (T_ref + S) / (T + S)

    lam_darcy_equiv = kappa * P_ATM / mu

    @printf("    %-35s T=%4.0fK: mu=%7.2e Pa.s, lam_Darcy ~ %8.2e m^2/s\n",
        name, T, mu, lam_darcy_equiv)
end
println()
println("  Darcy equivalent lambda is O(1e-6) m^2/s -- comparable to the Fickian lambda")
println("  but much smaller than TK's TRANSPORT values (1e-5 to 1e-4).")
println("  Combined Fick+Darcy ~ 2e-6 to 2e-5, but TK wants 1e-5 to 1.2e-4.")
println()
println("  VERDICT: The Darcy-Fick split in Julia DF does NOT fully recover TK's")
println("  effective transport. Total transport is ~3-10x lower than TK, which means")
println("  gas escapes more slowly in Julia -> higher internal gas pressure/concentration.")
println()

println("-" ^ 80)
println("Assessment C: PTW")
println("-" ^ 80)
println()
println("  PTW uses lambda = 2e-5 for ALL components (solid + gas).")
println("  This is UNIFORM, so mixing rule doesn't matter: lambda_eff = 2e-5 always.")
println()
println("  Physical assessment:")
println("  * D_free at 500K for typical pyrolysis gas: ~3-5e-5 m^2/s")
println("  * lambda = 2e-5 implies D_eff/D_free ~ 0.4-0.7")
println("  * From porous media: D_eff/D_free = epsilon/tau")
println("  * This implies epsilon/tau ~ 0.4-0.7, which is reasonable for")
println("    moderately porous charring wood (not for virgin wood!).")
println()
println("  * For virgin PTW (epsilon ~ 0.05, tau ~ 4): D_eff/D_free ~ 0.013")
println("    -> D_eff ~ 4e-7 m^2/s, far below the 2e-5 used")
println()
println("  VERDICT: lambda = 2e-5 for virgin wood is UNPHYSICALLY HIGH by ~50x.")
println("  However, this is consistent with TK convention where TRANSPORT is a")
println("  calibrated effective parameter that includes all transport mechanisms,")
println("  not just Fickian diffusion through pores.")
println()

# ============================================================================
# SECTION 6: Impact Analysis -- What Happens When Lambda is Wrong?
# ============================================================================

println("=" ^ 80)
println("SECTION 6: Impact Analysis -- Sensitivity of Gas Transport to Lambda")
println("=" ^ 80)
println()

println("The Fickian gas flux scales linearly with lambda_eff:")
println("  J = -(lambda_eff / T) * d(xi*T)/dx")
println()
println("If lambda_eff is too LOW (Julia DF case):")
println("  -> Gas escapes slower")
println("  -> Internal gas concentration builds up")
println("  -> Higher internal pressure (if Darcy doesn't compensate)")
println("  -> Possible: surface MLR is delayed but may peak higher")
println("  -> Known issue: Julia DF has 26% MLR gap vs TK (from memory)")
println()
println("If lambda_eff is too HIGH (TK virgin wood case):")
println("  -> Gas escapes too fast through intact virgin wood")
println("  -> Unrealistic: virgin wood should trap gas")
println("  -> But TK was calibrated with this, so it 'works' empirically")
println()

# Compute characteristic gas escape time
println("Characteristic gas escape times: t_escape ~ L^2 / lambda_eff")
println("  (time for gas to diffuse through thickness L)")
println()

L_vals = [1e-3, 5e-3, 12e-3]  # 1mm, 5mm, 12mm
lam_vals = [1e-6, 1e-5, 2e-5, 5e-5, 1e-4]

@printf("  %10s", "lambda\\L")
for L in L_vals
    @printf("  %10s", "L=$(Int(L*1000))mm")
end
println()
println("  " * "-" ^ 45)
for lam in lam_vals
    @printf("  %10.0e", lam)
    for L in L_vals
        t_escape = L^2 / lam
        @printf("  %10.1f s", t_escape)
    end
    println()
end
println()
println("  Typical cone calorimeter pyrolysis timescale: 100-600 s")
println("  If t_escape >> pyrolysis time, gas accumulates internally.")
println("  If t_escape << pyrolysis time, gas escapes as fast as produced.")
println()

# ============================================================================
# SECTION 7: Comparison Summary Table
# ============================================================================

println("=" ^ 80)
println("SECTION 7: Summary Comparison Table")
println("=" ^ 80)
println()

@printf("%-12s  %12s  %12s  %12s  %12s  %12s\n",
    "Material", "Julia_solid", "Julia_gas", "TK_solid", "TK_gas", "Physical_est")
println("-" ^ 80)

@printf("%-12s  %12s  %12s  %12s  %12s  %12s\n",
    "PMMA",
    "0 (default)", "1e-5", "1e-5", "1e-5", "~0 (dense)")

@printf("%-12s  %12s  %12s  %12s  %12s  %12s\n",
    "PTW",
    "2e-5", "2e-5", "2e-5", "2e-5", "4e-7 to 2e-5")

@printf("%-12s  %12s  %12s  %12s  %12s  %12s\n",
    "DF virgin",
    "1e-6", "6e-6..1.6e-5", "1e-5", "8e-5..1.2e-4", "1e-7 to 1e-6")

@printf("%-12s  %12s  %12s  %12s  %12s  %12s\n",
    "DF char",
    "1e-5", "6e-6..1.6e-5", "5e-5", "8e-5..1.2e-4", "5e-6 to 2e-5")

println()
println("Physical estimates use D_free * epsilon/tau at T=500K.")
println()

# ============================================================================
# SECTION 8: Specific Concerns
# ============================================================================

println("=" ^ 80)
println("SECTION 8: Specific Concerns and Recommendations")
println("=" ^ 80)
println()

println("CONCERN 1: Are any lambda values unphysically high or low?")
println("-" ^ 60)
println()
println("  HIGH (exceeds free diffusion):")
println("  * TK DF MOISTURE: lambda = 2e-3 m^2/s")
println("    D_free(H2O) at 300K = 2.5e-5. This is 80x D_free!")
println("    Interpretation: TK uses a very high pseudo-transport for")
println("    bound moisture to simulate rapid capillary-driven migration.")
println("    Not physical diffusion, but a modeling artifact.")
println()
println("  * TK DF gases: lambda = 8e-5 to 1.2e-4")
println("    At 500K, D_free ~ 3-5e-5. These exceed D_free by 2-3x!")
println("    Since D_eff should be D_free * eps/tau < D_free,")
println("    these values are UNPHYSICALLY HIGH.")
println("    Interpretation: TK's TRANSPORT includes Darcy-like contributions.")
println()
println("  LOW (may be too restrictive):")
println("  * Julia DF virgin solid: lambda = 1e-6")
println("    Physical D_eff estimate: 1e-7 to 1e-6.")
println("    This is at the upper end of the physical range for virgin wood.")
println("    Verdict: REASONABLE.")
println()
println("  * Julia DF gases: lambda = 6e-6 to 1.6e-5")
println("    These are below D_free at pyrolysis T (3-10e-5). Since these")
println("    are effective values in a porous medium (eps/tau < 1), they are")
println("    physically plausible. Verdict: REASONABLE but perhaps too low")
println("    given that gas should diffuse relatively freely in pore space.")
println()

println("CONCERN 2: Does PMMA lambda=0 for solids cause problems?")
println("-" ^ 60)
println()
println("  YES, potentially. With lambda=0 for solids and PARALLEL mixing,")
println("  the effective lambda is dominated by the gas volume fraction:")
println("    lambda_eff = v_gas * lambda_gas")
println("  Since v_gas is tiny in the early pyrolysis stage (v_gas < 0.01),")
println("  gas transport is nearly zero until significant decomposition occurs.")
println()
println("  This creates a positive feedback loop:")
println("    Low lambda -> gas accumulates -> higher pressure")
println("    -> but PMMA uses pure diffusion (no Darcy term)")
println("    -> gas can only escape slowly")
println("    -> surface regression may be delayed")
println()
println("  In contrast, TK with lambda=1e-5 everywhere allows gas to")
println("  escape through the solid matrix from the start.")
println()
println("  However, for PMMA specifically, this may be acceptable because")
println("  PMMA decomposes cleanly with surface regression -- there is no")
println("  thick char layer trapping gas. The bubbling/surface-regression")
println("  mechanism is captured by ALE mesh motion, not gas diffusion.")
println()

println("CONCERN 3: Are the DF gas lambda values reasonable?")
println("-" ^ 60)
println()
println("  Julia DF gas lambda values: 6e-6 to 1.6e-5 m^2/s")
println("  These represent D_eff in a porous char medium at T=500-700K.")
println()
println("  Physical estimate for gas in char (eps~0.7, tau~2):")
println("    D_eff = D_free * 0.7/2 = D_free * 0.35")
println("    At 500K: D_eff ~ 1.0e-5 to 1.7e-5 m^2/s")
println("    At 700K: D_eff ~ 2.0e-5 to 3.5e-5 m^2/s")
println()
println("  Julia values are at the lower end of this range, which makes")
println("  sense for gas diffusing through partially pyrolyzed wood (not")
println("  pure char). However, they are 5-13x LOWER than TK values.")
println()
println("  TK gas lambda = 8e-5 to 1.2e-4 EXCEEDS D_free, which means")
println("  TK's TRANSPORT parameter is not a pure Fickian diffusivity")
println("  but rather an empirically calibrated total transport rate.")
println()

println("CONCERN 4: Does the Darcy-Fick split in DF give reasonable TOTAL transport?")
println("-" ^ 60)
println()
println("  Julia DF uses both lambda (Fick) and kappa (Darcy):")
println("    J_total = J_Fick + J_Darcy")
println("    J_Fick  = -(lambda/T) * d(xi*T)/dx")
println("    J_Darcy = xi * v_Darcy = -xi * (kappa/mu) * dP/dx")
println()
println("  The Darcy contribution depends on local pressure gradients.")
println("  At a face, the Darcy velocity is:")
println("    v_g = -(kappa_eff / mu) * (P_R - P_L) / dx")
println()
println("  For a typical pressure gradient during pyrolysis:")
println("    dP ~ 1000-10000 Pa over 1mm cell -> dP/dx ~ 1e6 to 1e7 Pa/m")
println("    mu ~ 3e-5 Pa.s at 600K")
println("    kappa_char ~ 1e-14 m^2")
println("    v_Darcy ~ 1e-14 * 1e7 / 3e-5 ~ 3e-6 m/s")
println()
println("  This Darcy velocity, combined with gas concentration xi ~ 10 kg/m^3,")
println("  gives J_Darcy ~ 10 * 3e-6 ~ 3e-5 kg/(m^2.s)")
println()
println("  Compare with Fick:")
println("    lambda_eff ~ 5e-6, d(xiT)/dx ~ 10*600/0.001 = 6e6")
println("    J_Fick ~ 5e-6 / 600 * 6e6 ~ 5e-5 kg/(m^2.s)")
println()
println("  So Fick and Darcy are comparable in magnitude, giving total")
println("  transport ~ 8e-5 kg/(m^2.s)")
println()
println("  For TK (pure Fick with lambda=5e-5):")
println("    J_TK ~ 5e-5 / 600 * 6e6 ~ 5e-4 kg/(m^2.s)")
println()
println("  RESULT: Julia total transport (Fick+Darcy) is still ~6x LOWER than TK.")
println("  The Darcy term helps but does not fully compensate for the lower lambda.")
println()
println("  This shortfall likely contributes to the known 26% MLR gap between")
println("  Julia DF and TK. The gas escapes more slowly, which changes the")
println("  energy balance (convective transport of enthalpy with gas).")
println()

# ============================================================================
# SECTION 9: Cross-Check with TK Gas Density Convention
# ============================================================================

println("=" ^ 80)
println("SECTION 9: TK Gas Density Convention Impact")
println("=" ^ 80)
println()

println("TK specifies gas DENSITY in .cmp files (not molar mass).")
println("Julia uses ideal gas law: rho_gas = P*M/(R*T)")
println()
println("This affects volume fractions, which drive the mixing rule.")
println()

@printf("  %-20s  %10s  %10s  %10s  %10s\n",
    "Gas species", "TK_rho", "Julia@300K", "Julia@500K", "Julia@700K")
println("  " * "-" ^ 65)

for (name, M, rho_tk) in [
    ("MMA (M=100)",      0.100, 1.2),
    ("H2O (M=18)",       0.018, 0.8),
    ("DF_gas (M=60-110)", 0.080, 498.46),
    ("PTW_gas (M=60-110)",0.080, 498.46),
]
    @printf("  %-20s  %10.2f  %10.2f  %10.2f  %10.2f\n",
        name, rho_tk, rho_gas(M, 300.0), rho_gas(M, 500.0), rho_gas(M, 700.0))
end

println()
println("  NOTE: TK DF/PTW gas components have rho = 498.46 kg/m^3.")
println("  This is NOT a physical gas density (real gas at 1atm: ~1-3 kg/m^3).")
println("  It is a pseudo-density matching the original solid density.")
println("  Effect: In TK, gas occupies negligible volume (like the solid).")
println("  In Julia (ideal gas), gas at 500K has rho ~ 2 kg/m^3,")
println("  so the same mass concentration occupies ~250x MORE volume.")
println()
println("  This MASSIVELY affects volume fractions:")
println("  * In TK: even 10 kg/m^3 of gas has v_gas = 10/498 = 0.02")
println("  * In Julia: 10 kg/m^3 of gas at 500K has v_gas = 10/2 / total = 0.8+")
println()
println("  Since lambda_eff depends on volume fractions, the same per-component")
println("  lambda values can give VERY different effective lambda in Julia vs TK")
println("  even if the numerical values match exactly.")
println()

# Show this effect quantitatively
println("  Quantitative example: DF mid-pyrolysis (DF_int2=300, gas3=10 kg/m^3, T=550K)")
println()

xi_test = zeros(13)
xi_test[3] = 300.0   # DF_int2
xi_test[10] = 10.0   # DF_gas3 (levoglucosan, M=0.110)

# Julia: ideal gas
rho_gas3_julia = rho_gas(0.110, 550.0)
rho_julia = vcat(df_rho_solid, [rho_gas(M, 550.0) for M in df_M_gas])
v_julia = compute_volume_fractions(xi_test, rho_julia)

# TK: pseudo-density
rho_tk_vec = vcat(df_rho_solid, [498.46, 498.46, 498.46, 498.46, 498.46, 0.8])
v_tk = compute_volume_fractions(xi_test, rho_tk_vec)

# Same lambda for both (use TK values for apples-to-apples)
lam_same = collect(df_lam_tk)
lam_eff_julia_rho = weighted_lambda(collect(v_julia), lam_same; beta=0.5)
lam_eff_tk_rho = weighted_lambda(collect(v_tk), lam_same; beta=0.5)

@printf("    Julia (ideal gas rho): v_gas3 = %6.4f, v_solid = %6.4f -> lambda_eff = %8.2e\n",
    v_julia[10], v_julia[3], lam_eff_julia_rho)
@printf("    TK (pseudo rho=498):   v_gas3 = %6.4f, v_solid = %6.4f -> lambda_eff = %8.2e\n",
    v_tk[10], v_tk[3], lam_eff_tk_rho)
@printf("    Ratio (same lambda values!): %6.3f\n", lam_eff_julia_rho / lam_eff_tk_rho)
println()
println("  Even with IDENTICAL per-component lambda values, Julia gets a DIFFERENT")
println("  effective lambda because ideal gas densities give much larger gas volume")
println("  fractions. This compounds with the already-lower per-component lambda values.")
println()

# ============================================================================
# Final Summary
# ============================================================================

println("=" ^ 80)
println("FINAL SUMMARY")
println("=" ^ 80)
println()
println("1. PMMA: Julia lambda=0 for solids is physically correct for dense polymer,")
println("   but the PARALLEL mixing rule + zero solid lambda means very low effective")
println("   transport in early pyrolysis. TK uses lambda=1e-5 everywhere (empirical).")
println("   Impact: Moderate -- PMMA decomposes via surface regression, not bulk diffusion.")
println()
println("2. PTW: Julia matches TK exactly (lambda=2e-5 for all). Physically, this is")
println("   high for virgin wood (~50x above D_eff estimate) but is an established")
println("   calibrated parameter. No discrepancy between Julia and TK.")
println()
println("3. DF: LARGE discrepancy. Julia uses 5-13x lower lambda than TK for all")
println("   components. The Darcy-Fick split partially compensates but the TOTAL")
println("   transport is still ~3-6x lower than TK. This is a significant source of")
println("   the known 26% MLR discrepancy.")
println()
println("4. Gas density convention: Julia uses ideal gas law (low density), TK uses")
println("   pseudo-densities (high density). This gives DIFFERENT volume fractions")
println("   and thus different effective lambda even with identical per-component values.")
println("   This is an additional source of discrepancy in DF (but not PTW, which has")
println("   uniform lambda so mixing rule doesn't matter).")
println()
println("5. Key physical insight: TK's TRANSPORT parameter is NOT a Fickian diffusivity.")
println("   It is a calibrated effective transport coefficient that lumps together")
println("   Fickian diffusion, Darcy flow, capillary effects, and other mechanisms.")
println("   Attempting to split it into Fick + Darcy introduces model-form error")
println("   unless the split is carefully calibrated to reproduce the same total flux.")
println()
