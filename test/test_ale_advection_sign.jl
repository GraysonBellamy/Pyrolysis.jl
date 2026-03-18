"""
Manufactured Solution Test for ALE Advection Sign Convention

This test definitively determines whether the ALE advection term w·∇T should be
ADDED to or SUBTRACTED from the Eulerian rate of change.

==============================
MATHEMATICAL DERIVATION
==============================

The ALE (Arbitrary Lagrangian-Eulerian) framework relates the rate of change
following a mesh point to the Eulerian (spatial) rate of change.

Let χ denote the mesh coordinate and x the spatial coordinate.
A mesh point at position x₀ moves with velocity w, so:
    x(t) = x₀ + w·t

The rate of change seen by the mesh point is the material derivative following
the mesh motion:
    df/dt|_χ = ∂f/∂t|_x + w · ∂f/∂x

This is the standard ALE relation. With no physics (∂f/∂t|_x = 0):
    df/dt|_χ = w · ∂f/∂x

So the ODE for the mesh-point value should be:
    dT_i/dt = (Eulerian RHS) + w · ∂T/∂x

The ALE correction is ADDED (positive w·∇T).

==============================
TEST SETUP
==============================

- Domain: [0, L] with L = 0.01 m
- N = 20 uniform cells
- Initial profile: T(x) = 300 + 500·sin(π·x/L)
- Mesh velocity: w = -0.001 m/s (uniform, mesh contracting)
- No physics: only ALE advection
- Time step: Δt = 0.001 s (small enough for Euler to be accurate)

The mesh point at position x_i moves to x_i + w·Δt.
Since there's no physics, the temperature field T(x) is frozen in space.
The mesh point now samples a different location, so:

    T_i(t + Δt) = T(x_i + w·Δt)
                 ≈ T(x_i) + w·Δt · (∂T/∂x)|_{x_i}

Therefore: dT_i/dt = w · ∂T/∂x    (w·∇T is ADDED)

If the code SUBTRACTS w·∇T, we get: dT_i/dt = -w · ∂T/∂x (WRONG)
"""

using Printf

# ============================================================================
# Parameters
# ============================================================================
const L = 0.01      # Domain length [m]
const N = 20        # Number of cells
const Δx = L / N    # Cell size [m]
const w = -0.001    # Mesh velocity [m/s] (contracting)
const Δt = 0.001    # Time step [s]

# Cell centers
x = [(i - 0.5) * Δx for i in 1:N]

# ============================================================================
# Initial temperature profile
# ============================================================================
T0 = [300.0 + 500.0 * sin(π * x[i] / L) for i in 1:N]

# ============================================================================
# Compute ∂T/∂x using upwind finite differences (matching the code)
# ============================================================================
# With w < 0, information flows from right to left → forward difference
dTdx = zeros(N)
for i in 1:N
    if w > 0
        # Backward difference (upwind for w > 0)
        if i == 1
            dTdx[i] = 0.0  # Zero-gradient BC (matches code)
        else
            dTdx[i] = (T0[i] - T0[i-1]) / Δx
        end
    else
        # Forward difference (upwind for w < 0)
        if i == N
            dTdx[i] = 0.0  # Zero-gradient BC (matches code)
        else
            dTdx[i] = (T0[i+1] - T0[i]) / Δx
        end
    end
end

# ============================================================================
# Compute the ALE advection term (this is what the code stores)
# ============================================================================
ale_advection = w .* dTdx  # = w · ∂T/∂x

# ============================================================================
# Two candidate formulations
# ============================================================================
# Option A: SUBTRACT w·∇T  (what the code does for energy: dT -= ale_advection)
dT_subtract = -ale_advection  # dT/dt = -w·∂T/∂x

# Option B: ADD w·∇T (what the code does for species: dξ += ale_advection)
dT_add = ale_advection  # dT/dt = +w·∂T/∂x

# ============================================================================
# Forward Euler update
# ============================================================================
T_subtract = T0 .+ Δt .* dT_subtract
T_add = T0 .+ Δt .* dT_add

# ============================================================================
# Analytical solution: T_exact(x_i + w*Δt)
# Each mesh point moves to x_i + w*Δt and samples the frozen field
# ============================================================================
T_exact = [300.0 + 500.0 * sin(π * (x[i] + w * Δt) / L) for i in 1:N]

# Also compute with centered differences for a more accurate gradient
dTdx_centered = zeros(N)
for i in 1:N
    if i == 1
        dTdx_centered[i] = (T0[i+1] - T0[i]) / Δx
    elseif i == N
        dTdx_centered[i] = (T0[i] - T0[i-1]) / Δx
    else
        dTdx_centered[i] = (T0[i+1] - T0[i-1]) / (2 * Δx)
    end
end

T_add_centered = T0 .+ Δt .* (w .* dTdx_centered)
T_subtract_centered = T0 .+ Δt .* (-w .* dTdx_centered)

# ============================================================================
# Compute errors (interior cells only, excluding boundary cells with BCs)
# ============================================================================
interior = 2:N-1

err_subtract = sqrt(sum((T_subtract[i] - T_exact[i])^2 for i in interior) / length(interior))
err_add = sqrt(sum((T_add[i] - T_exact[i])^2 for i in interior) / length(interior))
err_add_centered = sqrt(sum((T_add_centered[i] - T_exact[i])^2 for i in interior) / length(interior))
err_subtract_centered = sqrt(sum((T_subtract_centered[i] - T_exact[i])^2 for i in interior) / length(interior))

# ============================================================================
# Print results
# ============================================================================
println("=" ^ 80)
println("MANUFACTURED SOLUTION TEST FOR ALE ADVECTION SIGN")
println("=" ^ 80)
println()
println("Setup:")
println("  Domain: [0, $L] m, $N cells, Δx = $Δx m")
println("  Profile: T(x) = 300 + 500·sin(πx/L)")
println("  Mesh velocity: w = $w m/s (contracting)")
println("  Time step: Δt = $Δt s")
println("  Mesh displacement: w·Δt = $(w*Δt) m = $(w*Δt/Δx) cells")
println()

println("-" ^ 80)
println("Per-cell comparison (upwind gradient, matching code discretization):")
println("-" ^ 80)
@printf("%-6s %-10s %-12s %-12s %-12s %-14s %-14s\n",
    "Cell", "x [mm]", "T_exact", "T_subtract", "T_add", "err_subtract", "err_add")
for i in 1:N
    marker = (i == 1 || i == N) ? " (BC)" : ""
    @printf("%-6d %-10.4f %-12.4f %-12.4f %-12.4f %-14.2e %-14.2e%s\n",
        i, x[i]*1000, T_exact[i], T_subtract[i], T_add[i],
        abs(T_subtract[i] - T_exact[i]),
        abs(T_add[i] - T_exact[i]),
        marker)
end

println()
println("=" ^ 80)
println("RMS ERRORS (interior cells only, excluding boundary effects):")
println("=" ^ 80)
println()
@printf("  SUBTRACT w·∇T (code's energy: dT -= ale)  : RMS = %.6e K\n", err_subtract)
@printf("  ADD      w·∇T (code's species: dξ += ale)  : RMS = %.6e K\n", err_add)
println()
@printf("  SUBTRACT w·∇T (centered gradient)          : RMS = %.6e K\n", err_subtract_centered)
@printf("  ADD      w·∇T (centered gradient)           : RMS = %.6e K\n", err_add_centered)

println()
println("=" ^ 80)
println("INTERPRETATION:")
println("=" ^ 80)
println()

if err_add < err_subtract
    ratio = err_subtract / err_add
    println("  ✓ ADD w·∇T is CORRECT ($(round(ratio, digits=1))× lower error)")
    println()
    println("  The code's SPECIES formulation (dξ += ale_species_advection) is CORRECT.")
    println("  The code's ENERGY formulation (dT -= ale_advection) is WRONG.")
    println()
    println("  FIX: Change line 895 from:")
    println("    dT[i] -= ale_advection[i]")
    println("  to:")
    println("    dT[i] += ale_advection[i]")
elseif err_subtract < err_add
    ratio = err_add / err_subtract
    println("  ✓ SUBTRACT w·∇T is CORRECT ($(round(ratio, digits=1))× lower error)")
    println()
    println("  The code's ENERGY formulation (dT -= ale_advection) is CORRECT.")
    println("  The code's SPECIES formulation (dξ += ale_species_advection) is WRONG.")
    println()
    println("  FIX: Change lines 902-903 from:")
    println("    dξ[j, i] += ale_species_advection[j, i]")
    println("  to:")
    println("    dξ[j, i] -= ale_species_advection[j, i]")
else
    println("  Errors are equal — inconclusive (this shouldn't happen).")
end

# ============================================================================
# ADDITIONAL VERIFICATION: Direct Taylor expansion check
# ============================================================================
println()
println("=" ^ 80)
println("VERIFICATION: Direct Taylor expansion at cell 10")
println("=" ^ 80)
println()

i_test = 10
x_test = x[i_test]
T_at_x = 300.0 + 500.0 * sin(π * x_test / L)
T_at_x_shifted = 300.0 + 500.0 * sin(π * (x_test + w * Δt) / L)
dTdx_exact = 500.0 * π / L * cos(π * x_test / L)

T_taylor_add = T_at_x + w * Δt * dTdx_exact
T_taylor_sub = T_at_x - w * Δt * dTdx_exact

println("  Mesh point at x = $(x_test*1000) mm")
println("  After Δt, mesh point is at x + w·Δt = $((x_test + w*Δt)*1000) mm")
println()
println("  T(x)            = $(T_at_x)")
println("  T(x + w·Δt)     = $(T_at_x_shifted)     [EXACT: what mesh point should see]")
println("  ∂T/∂x (exact)   = $(dTdx_exact)")
println()
println("  T + (+w·Δt)·∂T/∂x = $(T_taylor_add)  [ADD formulation]")
println("  T + (-w·Δt)·∂T/∂x = $(T_taylor_sub)  [SUBTRACT formulation]")
println()

add_err = abs(T_taylor_add - T_at_x_shifted)
sub_err = abs(T_taylor_sub - T_at_x_shifted)

@printf("  Error with ADD:      %.2e K\n", add_err)
@printf("  Error with SUBTRACT: %.2e K\n", sub_err)
println()

if add_err < sub_err * 0.01
    println("  CONCLUSION: ADD formulation matches Taylor expansion exactly (to machine precision).")
    println("              The ALE correction is: dT/dt = (Eulerian RHS) + w·∂T/∂x")
elseif sub_err < add_err * 0.01
    println("  CONCLUSION: SUBTRACT formulation matches Taylor expansion exactly.")
    println("              The ALE correction is: dT/dt = (Eulerian RHS) - w·∂T/∂x")
end

# ============================================================================
# ADDITIONAL: Multi-step integration to rule out sign-flip artifacts
# ============================================================================
println()
println("=" ^ 80)
println("MULTI-STEP VERIFICATION (100 Euler steps)")
println("=" ^ 80)
println()

n_steps = 100
dt_small = 0.0001  # Very small step for accuracy

# Integrate with ADD formulation
T_ms_add = copy(T0)
# Integrate with SUBTRACT formulation
T_ms_sub = copy(T0)

for step in 1:n_steps
    # Current positions (mesh has moved)
    t_now = step * dt_small

    # Compute gradients for current T fields (upwind for w < 0 → forward diff)
    for i in 1:N
        if w > 0
            dTdx_i = (i == 1) ? 0.0 : (T_ms_add[i] - T_ms_add[i-1]) / Δx
        else
            dTdx_i = (i == N) ? 0.0 : (T_ms_add[i+1] - T_ms_add[i]) / Δx
        end
        T_ms_add[i] += dt_small * (w * dTdx_i)  # ADD
    end

    for i in 1:N
        if w > 0
            dTdx_i = (i == 1) ? 0.0 : (T_ms_sub[i] - T_ms_sub[i-1]) / Δx
        else
            dTdx_i = (i == N) ? 0.0 : (T_ms_sub[i+1] - T_ms_sub[i]) / Δx
        end
        T_ms_sub[i] -= dt_small * (w * dTdx_i)  # SUBTRACT
    end
end

total_time = n_steps * dt_small
T_exact_ms = [300.0 + 500.0 * sin(π * (x[i] + w * total_time) / L) for i in 1:N]

err_ms_add = sqrt(sum((T_ms_add[i] - T_exact_ms[i])^2 for i in interior) / length(interior))
err_ms_sub = sqrt(sum((T_ms_sub[i] - T_exact_ms[i])^2 for i in interior) / length(interior))

@printf("  After %d steps (total time = %.4f s, mesh displacement = %.4f mm):\n",
    n_steps, total_time, w * total_time * 1000)
println()
@printf("  ADD formulation      RMS error: %.6e K\n", err_ms_add)
@printf("  SUBTRACT formulation RMS error: %.6e K\n", err_ms_sub)
println()

if err_ms_add < err_ms_sub
    println("  Multi-step confirms: ADD is correct.")
    # Check if SUBTRACT is diverging
    if err_ms_sub > 10 * err_ms_add
        println("  SUBTRACT formulation error is $(round(err_ms_sub/err_ms_add, digits=1))× larger — ")
        println("  it advects temperature in the WRONG direction.")
    end
else
    println("  Multi-step confirms: SUBTRACT is correct.")
end

# ============================================================================
# PHYSICAL INTUITION CHECK
# ============================================================================
println()
println("=" ^ 80)
println("PHYSICAL INTUITION CHECK")
println("=" ^ 80)
println()
println("  Consider a cell at the peak of the sine profile (x = L/2).")
println("  T = 800 K (maximum), ∂T/∂x = 0.")
println("  Neighboring cells have T < 800 K.")
println()
println("  With w = -0.001 m/s, the mesh point moves LEFT (toward lower x).")
println("  To the left, T < 800 K (descending slope of sine).")
println("  So the mesh point should see a DECREASING temperature.")
println()
println("  At cell to the LEFT of peak (e.g., x = 0.45·L):")
println("  T < 800, ∂T/∂x > 0 (ascending slope)")
println("  w < 0, so w·∂T/∂x < 0 → temperature DECREASES")

i_left_of_peak = 9  # x ≈ 0.425*L
dTdx_at_9 = (T0[10] - T0[9]) / Δx
println()
@printf("  Cell %d: x = %.4f mm, T = %.2f K\n", i_left_of_peak, x[i_left_of_peak]*1000, T0[i_left_of_peak])
@printf("  ∂T/∂x ≈ %.1f K/m (positive, ascending slope)\n", dTdx_at_9)
@printf("  w·∂T/∂x = %.4f K/s (negative → T decreases)\n", w * dTdx_at_9)
println()
println("  Mesh point moves LEFT into cooler region → T should DECREASE ✓")
println("  ADD formulation gives dT/dt = w·∂T/∂x < 0 → T decreases ✓")
println("  SUBTRACT formulation gives dT/dt = -w·∂T/∂x > 0 → T increases ✗")

println()
println("=" ^ 80)
println("FINAL VERDICT")
println("=" ^ 80)
println()
println("  The ALE time derivative relation is:")
println()
println("    df/dt|_mesh = ∂f/∂t|_spatial + w · ∂f/∂x")
println()
println("  So the correction w·∇f should be ADDED to the Eulerian RHS:")
println("    dT/dt = (Eulerian RHS) + w·∂T/∂x")
println()
println("  SIGN CONSISTENCY CHECK IN THE CODE:")
println("  • Energy  (line 895): dT[i] -= ale_advection[i]  →  dT = RHS - w·∇T  [WRONG]")
println("  • Species (line 903): dξ[j,i] += ale_advection[j,i]  →  dξ = RHS + w·∇ξ  [CORRECT*]")
println()
println("  *For gas species with v_gas=0: the species equation is")
println("    dξ/dt|_mesh = ∂ξ/∂t|_spatial + w·∂ξ/∂x")
println("  BUT the transport is: ∂ξ/∂t|_spatial = S - ∇·(v_gas · ξ) = S  (with v_gas=0)")
println("  And the ALE species equation from the code comment says:")
println("    ∂ξ/∂t|_χ + (v - w)·∇ξ = S")
println("  With v=0: ∂ξ/∂t|_χ = S + w·∇ξ  ← This is CORRECT")
println()
println("  HOWEVER, the energy equation comment (line 151 of ale_fluxes.jl) says:")
println("    ρcₚ(∂T/∂t + w·∇T) = -∇·q + Q")
println("  Rearranging: ∂T/∂t|_χ = (-∇·q + Q)/(ρcₚ) - w·∇T  ← This claims SUBTRACT!")
println()
println("  THE COMMENT IS WRONG. The correct ALE relation is:")
println("    ∂T/∂t|_χ = ∂T/∂t|_x + w·∇T")
println("  NOT:")
println("    ∂T/∂t|_χ + w·∇T = ∂T/∂t|_x")
println()
println("  Wait — let's be careful about which side has which derivative.")
println()
println("  Standard ALE identity (e.g., Donea & Huerta, 'Finite Elements for Flow Problems'):")
println("    ∂f/∂t|_χ = ∂f/∂t|_x + w · ∇f")
println()
println("  If the Eulerian PDE is:  ∂T/∂t|_x = RHS")
println("  Then in the ALE frame:   ∂T/∂t|_χ = RHS + w · ∇T")
println()
println("  So the code should do: dT[i] = eulerian_RHS + w·∇T")
println("  i.e., dT[i] += ale_advection[i]  (not -=)")
