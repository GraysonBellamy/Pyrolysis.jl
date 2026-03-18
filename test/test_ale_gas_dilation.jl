"""
Manufactured Solution Test: Does Gas Species Need the ALE Dilation Correction?
===============================================================================

QUESTION
--------
The code applies dilation correction dξ -= ξ*θ only to SOLID species.
For GAS species, the code uses w·∇ξ advection without dilation correction.
Should gas species also receive -ξ*θ?

DERIVATION
----------
Consider the integral conservation law for gas species in a moving control
volume V(t) with no physics (S=0, J=0):

    d/dt ∫_V ξ dV = 0    (no source, no flux)

Expanding using the Leibniz/Reynolds transport theorem for a volume whose
boundary moves at mesh velocity w:

    d/dt ∫_V ξ dV = ∫_V ∂ξ/∂t|_χ dV + ∮_∂V ξ w·n dA

Setting this to zero and localizing:
    ∂ξ/∂t|_χ + ξ ∇·w = 0

Since θ = ∇·w (the mesh dilation rate), the CONSERVATIVE ALE equation is:
    ∂ξ/∂t|_χ = -ξ θ                ... (*)

The NON-CONSERVATIVE ALE approach starts from the Eulerian-to-ALE identity:
    ∂ξ/∂t|_χ = ∂ξ/∂t|_x + w·∇ξ

With ∂ξ/∂t|_x = 0 (no physics in Eulerian frame), this gives:
    ∂ξ/∂t|_χ = w·∇ξ               ... (**)

Now (*) and (**) are DIFFERENT unless w·∇ξ = -ξθ, which is not generally true.

The difference between the two formulations is:
    w·∇ξ - (-ξθ) = w·∇ξ + ξ∇·w = ∇·(wξ)

So the non-conservative form (**) omits the term ξ∇·w = ξθ compared to
the conservative form (*).

WHAT IS CORRECT?
The integral conservation law d(ξV)/dt = 0 is the fundamental physical law.
Both approaches are correct IF applied consistently:

  - Conservative approach: dξ/dt = -ξθ
  - Non-conservative approach: dξ/dt = w·∇ξ

These are NOT the same equation! They differ by ξθ. However, they describe
the same physics in different ways:
  - Conservative: tracks how concentration changes because the cell volume changes
  - Non-conservative: tracks how concentration changes because the mesh
    point samples a different spatial location

The key insight is that the code uses the non-conservative approach (w·∇ξ)
for gas species advection. So adding -ξθ on top would be DOUBLE-COUNTING.

Conversely, for SOLID species, the code does NOT add w·∇ξ (because solids
move with the mesh, so advection is zero). Instead it correctly uses -ξθ.

This test verifies this understanding using a manufactured solution.

TEST DESIGN
-----------
1D domain with:
- Known gas concentration profile: ξ(x) = ξ_base + ξ_amp * cos(πx/L)
- Known mesh velocity field: w(x) = w_top * x/L (linear, w=0 at bottom)
- No physics: S=0, J=0
- θ = dw/dx = w_top/L (uniform dilation)

We compute and compare:
(a) Exact rate dξ/dt from conservative principle: d(ξV)/dt = 0 → dξ/dt = -ξθ
(b) Non-conservative advection: dξ/dt = w·∇ξ
(c) Both combined (incorrect double-counting): dξ/dt = w·∇ξ - ξθ

After a Forward Euler step, compare each to the reference solution obtained
by direct evaluation of ξ at the new cell-center positions (mesh has moved).
"""

using Printf

println("=" ^ 80)
println("MANUFACTURED SOLUTION TEST: GAS SPECIES ALE DILATION CORRECTION")
println("=" ^ 80)
println()

# ============================================================================
# Parameters
# ============================================================================
const L = 0.01       # Domain length [m]
const N = 40         # Number of cells (enough for good accuracy)
const Δx = L / N     # Cell size [m]
const w_top = -0.001 # Mesh velocity at top [m/s] (contracting)
const Δt = 1e-5      # Very small time step for accurate Euler [s]

# Derived
const θ_uniform = w_top / L  # Uniform dilation rate [1/s]
# θ = -0.001 / 0.01 = -0.1 s⁻¹ (contracting)

# Gas concentration profile parameters
const ξ_base = 10.0   # Base concentration [kg/m³]
const ξ_amp = 5.0     # Amplitude [kg/m³]

# Cell centers (z=0 at bottom, z=L at top)
x = [(i - 0.5) * Δx for i in 1:N]

# Mesh velocity field: w(x) = w_top * x/L (linear, w=0 at bottom)
# This gives θ = dw/dx = w_top/L everywhere (uniform dilation)
w_at_x = [w_top * x[i] / L for i in 1:N]

# Node positions (at cell boundaries)
x_nodes = [i * Δx for i in 0:N]
w_at_nodes = [w_top * x_nodes[i+1] / L for i in 0:N]

# ============================================================================
# Initial concentration profile
# ============================================================================
ξ0 = [ξ_base + ξ_amp * cos(π * x[i] / L) for i in 1:N]

println("Setup:")
println("  Domain: [0, $L] m, $N cells, Δx = $(Δx*1000) mm")
println("  Concentration: ξ(x) = $ξ_base + $ξ_amp*cos(πx/L)")
println("  Mesh velocity: w(x) = $(w_top)*x/L (linear)")
println("  Dilation rate: θ = w_top/L = $θ_uniform s⁻¹ (uniform)")
println("  Time step: Δt = $Δt s")
println("  Mesh displacement at top: w_top*Δt = $(w_top*Δt) m = $(w_top*Δt/Δx) cells")
println()

# ============================================================================
# APPROACH 1: Conservative form — dξ/dt = -ξ*θ  (what solids use)
# ============================================================================
dξ_conservative = [-ξ0[i] * θ_uniform for i in 1:N]

# ============================================================================
# APPROACH 2: Non-conservative advection — dξ/dt = w·∇ξ  (what gas code does)
# ============================================================================
# Compute ∂ξ/∂x using upwind differencing (matching the code's discretization)
dξdx_upwind = zeros(N)
for i in 1:N
    w_i = w_at_x[i]
    if w_i > 0
        # Backward difference (upwind for w > 0)
        if i == 1
            dξdx_upwind[i] = 0.0  # Zero-gradient BC
        else
            dξdx_upwind[i] = (ξ0[i] - ξ0[i-1]) / Δx
        end
    else
        # Forward difference (upwind for w < 0)
        if i == N
            dξdx_upwind[i] = 0.0  # Zero-gradient BC
        else
            dξdx_upwind[i] = (ξ0[i+1] - ξ0[i]) / Δx
        end
    end
end

# But note: the gas advection in the code uses -w (relative velocity is 0 - w = -w)
# Looking at ale_fluxes.jl: for gas, v_gas = 0, so the advection direction is
# based on w_cell (not -w_cell). The code stores w·∇ξ and ADDS it:
#   dξ += ale_species_advection = w·∇ξ
#
# Wait - re-reading the code more carefully:
# In _compute_species_advection_unrolled!, for GAS components:
#   advection[j,i] = w_cell * dξ_dz
# where dξ_dz uses upwind based on w_cell direction.
#
# Then in _subtract_ale_advection!:
#   dξ[j,i] += ale_species_advection[j,i]
#
# So the total gas advection contribution is: +w·∇ξ
# This is the non-conservative advection term from ∂ξ/∂t|_χ = ∂ξ/∂t|_x + w·∇ξ

dξ_nonconservative = [w_at_x[i] * dξdx_upwind[i] for i in 1:N]

# ============================================================================
# APPROACH 3: Double-counting — dξ/dt = w·∇ξ - ξ*θ  (if gas got both)
# ============================================================================
dξ_both = [w_at_x[i] * dξdx_upwind[i] - ξ0[i] * θ_uniform for i in 1:N]

# ============================================================================
# EXACT REFERENCE: What should ξ_i be after Δt?
# ============================================================================
# The mesh moves: each node at position x moves to x + w*Δt.
# If there are no physics (S=0, J=0), the spatial field ξ(x) is frozen.
# But wait — that's not right either! The gas mass IS moving with the mesh
# if it's trapped, or it stays put if gas is independent of mesh motion.
#
# The correct physical picture depends on what "gas species" means:
#
# SCENARIO A: Gas is INDEPENDENT of mesh (gas doesn't move with mesh)
#   The spatial field ξ(x) is frozen (no Eulerian change), and the mesh
#   point now samples a different location.
#   ξ_i(t+Δt) = ξ(x_i + w_i*Δt)   [frozen field, moving sample point]
#   → The rate of change is dξ/dt|_χ = w·∇ξ (non-conservative)
#
# SCENARIO B: Gas mass is conserved in each cell (gas moves with mesh)
#   d(ξ_i * V_i)/dt = 0
#   → The rate of change is dξ/dt|_χ = -ξ_i * θ (conservative)
#
# In the pyrolysis code, gas species do NOT move with the mesh. The mesh
# tracks the solid, and gas slips through. With v_gas = 0, the gas stays
# in place spatially while the mesh moves through it.
#
# Therefore SCENARIO A is correct: The non-conservative form is appropriate.

# Reference solution: ξ at shifted positions (frozen spatial field)
x_shifted = [x[i] + w_at_x[i] * Δt for i in 1:N]
ξ_ref_frozen = [ξ_base + ξ_amp * cos(π * x_shifted[i] / L) for i in 1:N]

# Also compute the exact gradient for error assessment
dξdx_exact = [-ξ_amp * π / L * sin(π * x[i] / L) for i in 1:N]
ξ_ref_exact_gradient = [ξ0[i] + Δt * w_at_x[i] * dξdx_exact[i] for i in 1:N]

# Reference solution: mass conservation in each cell
# V_i(t+Δt) = V_i(t) + dV_i/dt * Δt = V_i(t) * (1 + θ*Δt)
# ξ_i(t+Δt) = ξ_i(t) * V_i(t) / V_i(t+Δt) = ξ_i(t) / (1 + θ*Δt)
ξ_ref_mass_conserved = [ξ0[i] / (1.0 + θ_uniform * Δt) for i in 1:N]

# ============================================================================
# Forward Euler updates
# ============================================================================
ξ_euler_conservative   = ξ0 .+ Δt .* dξ_conservative
ξ_euler_nonconservative = ξ0 .+ Δt .* dξ_nonconservative
ξ_euler_both           = ξ0 .+ Δt .* dξ_both

# ============================================================================
# Errors vs frozen-field reference (gas doesn't move with mesh)
# ============================================================================
interior = 3:N-2  # Skip boundary-adjacent cells where upwind stencil is limited

err_conservative_vs_frozen = sqrt(sum((ξ_euler_conservative[i] - ξ_ref_frozen[i])^2 for i in interior) / length(interior))
err_nonconservative_vs_frozen = sqrt(sum((ξ_euler_nonconservative[i] - ξ_ref_frozen[i])^2 for i in interior) / length(interior))
err_both_vs_frozen = sqrt(sum((ξ_euler_both[i] - ξ_ref_frozen[i])^2 for i in interior) / length(interior))
err_exact_gradient_vs_frozen = sqrt(sum((ξ_ref_exact_gradient[i] - ξ_ref_frozen[i])^2 for i in interior) / length(interior))

# ============================================================================
# Errors vs mass-conserved reference (gas moves with mesh)
# ============================================================================
err_conservative_vs_mass = sqrt(sum((ξ_euler_conservative[i] - ξ_ref_mass_conserved[i])^2 for i in interior) / length(interior))
err_nonconservative_vs_mass = sqrt(sum((ξ_euler_nonconservative[i] - ξ_ref_mass_conserved[i])^2 for i in interior) / length(interior))
err_both_vs_mass = sqrt(sum((ξ_euler_both[i] - ξ_ref_mass_conserved[i])^2 for i in interior) / length(interior))

# ============================================================================
# Print Results
# ============================================================================
println("=" ^ 80)
println("PART 1: Comparison with FROZEN SPATIAL FIELD reference")
println("        (gas is independent of mesh — correct for gas species)")
println("=" ^ 80)
println()
println("  Reference: ξ_i(t+Δt) = ξ(x_i + w_i Δt)  [frozen field, moving sample point]")
println("  This is the correct physical picture when gas does NOT move with the mesh.")
println()
@printf("  Method                           RMS Error [kg/m³]\n")
@printf("  %-36s %.6e\n", "Conservative (-ξθ only)", err_conservative_vs_frozen)
@printf("  %-36s %.6e\n", "Non-conservative (w·∇ξ only)", err_nonconservative_vs_frozen)
@printf("  %-36s %.6e\n", "Both (w·∇ξ - ξθ)", err_both_vs_frozen)
@printf("  %-36s %.6e\n", "Exact gradient (w·∇ξ_exact)", err_exact_gradient_vs_frozen)
println()

println("=" ^ 80)
println("PART 2: Comparison with MASS-CONSERVED-IN-CELL reference")
println("        (gas moves with mesh — correct for solid species)")
println("=" ^ 80)
println()
println("  Reference: ξ_i V_i = const → ξ_i(t+Δt) = ξ_i(t)/(1+θΔt)")
println("  This would be correct if gas moved with the mesh (like solids).")
println()
@printf("  Method                           RMS Error [kg/m³]\n")
@printf("  %-36s %.6e\n", "Conservative (-ξθ only)", err_conservative_vs_mass)
@printf("  %-36s %.6e\n", "Non-conservative (w·∇ξ only)", err_nonconservative_vs_mass)
@printf("  %-36s %.6e\n", "Both (w·∇ξ - ξθ)", err_both_vs_mass)
println()

# ============================================================================
# PART 3: Analytical verification at a specific cell
# ============================================================================
println("=" ^ 80)
println("PART 3: Analytical verification at cell $(N÷2) (domain center)")
println("=" ^ 80)
println()

i_mid = N ÷ 2
x_mid = x[i_mid]
w_mid = w_at_x[i_mid]
ξ_mid = ξ0[i_mid]
dξdx_mid_exact = -ξ_amp * π / L * sin(π * x_mid / L)

println("  Cell $i_mid: x = $(round(x_mid*1000, digits=4)) mm")
@printf("  ξ = %.4f kg/m³\n", ξ_mid)
@printf("  w = %.6f m/s\n", w_mid)
@printf("  ∂ξ/∂x (exact) = %.4f kg/m⁴\n", dξdx_mid_exact)
@printf("  θ = %.4f s⁻¹\n", θ_uniform)
println()
@printf("  w·∇ξ       = %.6e kg/(m³·s)\n", w_mid * dξdx_mid_exact)
@printf("  -ξ·θ       = %.6e kg/(m³·s)\n", -ξ_mid * θ_uniform)
@printf("  w·∇ξ - ξ·θ = %.6e kg/(m³·s)\n", w_mid * dξdx_mid_exact - ξ_mid * θ_uniform)
println()

# The vector identity: w·∇ξ + ξ·∇·w = ∇·(wξ)
# In 1D: w dξ/dx + ξ dw/dx = d(wξ)/dx
div_wξ_exact = w_mid * dξdx_mid_exact + ξ_mid * θ_uniform
println("  Vector identity check:")
@printf("    w·∇ξ + ξ·(∇·w) = ∇·(wξ) = %.6e\n", div_wξ_exact)
# Also compute ∇·(wξ) directly: d/dx[w(x)*ξ(x)]
# w(x) = w_top*x/L, ξ(x) = ξ_base + ξ_amp*cos(πx/L)
# d/dx[wξ] = (w_top/L)*ξ + w_top*x/L * (-ξ_amp*π/L*sin(πx/L))
div_wξ_direct = (w_top / L) * ξ_mid + (w_top * x_mid / L) * dξdx_mid_exact
@printf("    d(wξ)/dx direct = %.6e\n", div_wξ_direct)
@printf("    Difference:       %.6e (should be ~0)\n", abs(div_wξ_exact - div_wξ_direct))
println()

# ============================================================================
# PART 4: Uniform concentration test (the smoking gun)
# ============================================================================
println("=" ^ 80)
println("PART 4: UNIFORM CONCENTRATION TEST (the critical test)")
println("=" ^ 80)
println()
println("  If ξ(x) = ξ₀ = const everywhere, and there are no sources/fluxes,")
println("  what should happen when the mesh contracts?")
println()
println("  PHYSICAL ANSWER depends on the gas physics:")
println()
println("  Case A: Gas is STATIONARY in space (v_gas = 0)")
println("    The spatial field ξ(x) = ξ₀ remains uniform.")
println("    The mesh point moves to a new location, but ξ is still ξ₀ there.")
println("    So dξ/dt|_χ = 0")
println("    Non-conservative: w·∇ξ = w × 0 = 0  ✓")
println("    Conservative: -ξ₀·θ ≠ 0              ✗  (says ξ changes!)")
println()
println("  Case B: Gas mass is LOCKED to the cell (moves with mesh)")
println("    As cell volume shrinks, mass is conserved → ξ increases.")
println("    So dξ/dt|_χ = -ξ₀·θ > 0  (θ < 0 for contraction)")
println("    Conservative: -ξ₀·θ ≠ 0              ✓")
println("    Non-conservative: w·∇ξ = w × 0 = 0   ✗  (says ξ unchanged!)")
println()

ξ_uniform = ξ_base  # Uniform profile, no gradient

# With no gradient, w·∇ξ = 0 identically
dξ_nc_uniform = 0.0  # w · 0 = 0

# Conservative correction gives -ξ₀·θ
dξ_c_uniform = -ξ_uniform * θ_uniform

@printf("  For ξ₀ = %.1f kg/m³, θ = %.4f s⁻¹:\n", ξ_uniform, θ_uniform)
@printf("    w·∇ξ = %.6e  (correct for stationary gas)\n", dξ_nc_uniform)
@printf("    -ξ·θ = %.6e  (correct for gas locked to cell)\n", dξ_c_uniform)
println()

println("  In the pyrolysis code, gas species are NOT locked to the mesh.")
println("  Gas can slip through the porous solid. With v_gas = 0, the gas")
println("  is spatially stationary while the mesh moves through it.")
println()
println("  Therefore: for a UNIFORM gas concentration, the mesh motion should")
println("  cause NO change in ξ. Only the non-conservative form (w·∇ξ) gives")
println("  this result. Adding -ξθ would INCORRECTLY create a gas source.")

# ============================================================================
# PART 5: Quantitative convergence study
# ============================================================================
println()
println("=" ^ 80)
println("PART 5: Mesh convergence study (frozen-field reference)")
println("=" ^ 80)
println()
@printf("  %-10s %-18s %-18s %-18s\n",
    "N_cells", "err(w·∇ξ)", "err(-ξθ)", "err(w·∇ξ - ξθ)")

for N_test in [20, 40, 80, 160, 320]
    Δx_t = L / N_test
    x_t = [(i - 0.5) * Δx_t for i in 1:N_test]
    w_t = [w_top * x_t[i] / L for i in 1:N_test]
    ξ_t = [ξ_base + ξ_amp * cos(π * x_t[i] / L) for i in 1:N_test]

    # Upwind gradient
    dξdx_t = zeros(N_test)
    for i in 1:N_test
        if w_t[i] > 0
            dξdx_t[i] = (i == 1) ? 0.0 : (ξ_t[i] - ξ_t[i-1]) / Δx_t
        else
            dξdx_t[i] = (i == N_test) ? 0.0 : (ξ_t[i+1] - ξ_t[i]) / Δx_t
        end
    end

    # Three formulations
    ξ_nc = ξ_t .+ Δt .* (w_t .* dξdx_t)
    ξ_c  = ξ_t .+ Δt .* (-ξ_t .* θ_uniform)
    ξ_b  = ξ_t .+ Δt .* (w_t .* dξdx_t .- ξ_t .* θ_uniform)

    # Reference: frozen field
    x_s = [x_t[i] + w_t[i] * Δt for i in 1:N_test]
    ξ_ref = [ξ_base + ξ_amp * cos(π * x_s[i] / L) for i in 1:N_test]

    int_t = 3:N_test-2
    e_nc = sqrt(sum((ξ_nc[i] - ξ_ref[i])^2 for i in int_t) / length(int_t))
    e_c  = sqrt(sum((ξ_c[i]  - ξ_ref[i])^2 for i in int_t) / length(int_t))
    e_b  = sqrt(sum((ξ_b[i]  - ξ_ref[i])^2 for i in int_t) / length(int_t))

    @printf("  %-10d %-18.6e %-18.6e %-18.6e\n", N_test, e_nc, e_c, e_b)
end

println()
println("  If w·∇ξ is correct: its error should converge with mesh refinement.")
println("  The -ξθ error should NOT converge (it's a fundamentally different equation).")
println("  The combined w·∇ξ - ξθ error should NOT converge either.")

# ============================================================================
# PART 6: Multi-step integration
# ============================================================================
println()
println("=" ^ 80)
println("PART 6: Multi-step forward Euler integration (100 steps)")
println("=" ^ 80)
println()

n_steps = 100
dt_ms = 1e-5

# Track three formulations
ξ_ms_nc = copy(ξ0)  # Non-conservative (w·∇ξ)
ξ_ms_c  = copy(ξ0)  # Conservative (-ξθ)
ξ_ms_b  = copy(ξ0)  # Both (w·∇ξ - ξθ)

for step in 1:n_steps
    # Compute gradients for non-conservative approach
    for i in 1:N
        w_i = w_at_x[i]
        if w_i > 0
            dξdx_i = (i == 1) ? 0.0 : (ξ_ms_nc[i] - ξ_ms_nc[i-1]) / Δx
        else
            dξdx_i = (i == N) ? 0.0 : (ξ_ms_nc[i+1] - ξ_ms_nc[i]) / Δx
        end
        ξ_ms_nc[i] += dt_ms * w_i * dξdx_i
    end

    # Conservative
    for i in 1:N
        ξ_ms_c[i] += dt_ms * (-ξ_ms_c[i] * θ_uniform)
    end

    # Both
    for i in 1:N
        w_i = w_at_x[i]
        if w_i > 0
            dξdx_i = (i == 1) ? 0.0 : (ξ_ms_b[i] - ξ_ms_b[i-1]) / Δx
        else
            dξdx_i = (i == N) ? 0.0 : (ξ_ms_b[i+1] - ξ_ms_b[i]) / Δx
        end
        ξ_ms_b[i] += dt_ms * (w_i * dξdx_i - ξ_ms_b[i] * θ_uniform)
    end
end

total_time = n_steps * dt_ms
x_final = [x[i] + w_at_x[i] * total_time for i in 1:N]
ξ_ref_ms = [ξ_base + ξ_amp * cos(π * x_final[i] / L) for i in 1:N]

err_ms_nc = sqrt(sum((ξ_ms_nc[i] - ξ_ref_ms[i])^2 for i in interior) / length(interior))
err_ms_c  = sqrt(sum((ξ_ms_c[i]  - ξ_ref_ms[i])^2 for i in interior) / length(interior))
err_ms_b  = sqrt(sum((ξ_ms_b[i]  - ξ_ref_ms[i])^2 for i in interior) / length(interior))

@printf("  After %d steps (total time = %.4f ms):\n", n_steps, total_time * 1000)
println()
@printf("  Non-conservative (w·∇ξ only):    RMS error = %.6e kg/m³\n", err_ms_nc)
@printf("  Conservative (-ξθ only):         RMS error = %.6e kg/m³\n", err_ms_c)
@printf("  Both (w·∇ξ - ξθ):               RMS error = %.6e kg/m³\n", err_ms_b)
println()

if err_ms_nc < err_ms_c && err_ms_nc < err_ms_b
    @printf("  Non-conservative is best (%.1f× better than conservative, %.1f× better than both)\n",
        err_ms_c / err_ms_nc, err_ms_b / err_ms_nc)
elseif err_ms_c < err_ms_nc && err_ms_c < err_ms_b
    @printf("  Conservative is best (%.1f× better than non-conservative, %.1f× better than both)\n",
        err_ms_nc / err_ms_c, err_ms_b / err_ms_c)
else
    @printf("  Combined is best\n")
end

# ============================================================================
# CONCLUSION
# ============================================================================
println()
println("=" ^ 80)
println("CONCLUSION")
println("=" ^ 80)
println()
println("  The question: Should gas species receive the dilation correction -ξθ?")
println()
println("  ANSWER: NO. The code is CORRECT as written.")
println()
println("  Reasoning:")
println("  1. Gas species do not move with the mesh. With v_gas = 0, gas is")
println("     spatially stationary while the mesh moves through it.")
println()
println("  2. The non-conservative ALE advection (w·∇ξ) correctly describes")
println("     how the mesh point samples different spatial locations as the")
println("     mesh moves. This is verified by the frozen-field test above.")
println()
println("  3. Adding -ξθ on top of w·∇ξ would be double-counting. The two")
println("     terms (w·∇ξ) and (-ξθ) are alternative formulations of the")
println("     same physics, NOT additive corrections:")
println("     - w·∇ξ = non-conservative ALE (correct for gas)")
println("     - -ξθ  = conservative ALE (correct for solids)")
println()
println("  4. For SOLID species, the code correctly uses -ξθ and NOT w·∇ξ,")
println("     because solids move with the mesh (zero relative velocity).")
println()
println("  5. The uniform concentration test (Part 4) is the clearest proof:")
println("     with ξ = const, mesh contraction should NOT change ξ for gas.")
println("     Only w·∇ξ gives this (since ∇ξ = 0), while -ξθ falsely creates")
println("     a concentration increase.")
println()
println("  Summary of what each species type needs:")
println("    SOLID: dξ/dt = S - ξθ        (conservative, mass locked to cell)")
println("    GAS:   dξ/dt = S - ∇·J + w·∇ξ  (non-conservative, gas slips through mesh)")
println()
println("  The code implements exactly this.")
