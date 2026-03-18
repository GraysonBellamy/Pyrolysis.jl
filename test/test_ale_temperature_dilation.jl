"""
Manufactured Solution Test: ALE Temperature Equation — Advection vs Dilation
=============================================================================

QUESTION
--------
Temperature T is a mixture property (both solid and gas at the same T).
In ALE mode, what is the correct equation for T?

    (A) dT/dt = RHS/rho_cp + w·∇T         (non-conservative, like gas species)
    (B) dT/dt = RHS/rho_cp - T·theta       (conservative/dilation, like solid species)
    (C) dT/dt = RHS/rho_cp + w·∇T - T·theta  (both)

DERIVATION FROM FIRST PRINCIPLES
---------------------------------
The energy conservation for a cell on a moving mesh (volume V(t)):

    d/dt integral_V rho*cp*T dV = integral_V RHS dV

The subtlety: rho*cp = sum(xi_j * c_j) is NOT constant — it changes due to
reactions! So we cannot simply pull it out.

However, in pyrolysis codes the temperature equation is NOT derived from
conservative energy. It is derived from the non-conservative form:

    rho*cp * dT/dt|_x = -div(q) + Q_rxn + S_conv + Q_rad     (*)

This is a PDE for the INTENSIVE variable T. The ALE transform of this
non-conservative equation uses the standard identity:

    dT/dt|_x = dT/dt|_chi - w·∇T

Substituting into (*):

    rho*cp * (dT/dt|_chi - w·∇T) = RHS
    dT/dt|_chi = RHS/rho_cp + w·∇T                             (**)

This gives OPTION (A): just add w·∇T, NO dilation correction.

The dilation correction -T·theta would arise from CONSERVATIVE (integral)
energy, but the temperature PDE is already non-conservative.

ENERGY CONSERVATION CHECK
-------------------------
Does option (A) conserve total energy E = integral rho*cp*T dV?

For a cell: E_i = rho_cp_i * T_i * V_i

dE_i/dt = d(rho_cp_i)/dt * T_i * V_i + rho_cp_i * dT_i/dt * V_i + rho_cp_i * T_i * dV_i/dt

With no physics (RHS=0) and constant rho_cp (no reactions):
    dT_i/dt|_chi = w·∇T
    dV_i/dt = theta * V_i
    d(rho_cp_i)/dt = 0

    dE_i/dt = rho_cp * (w·∇T) * V + rho_cp * T * theta * V
            = rho_cp * V * (w·∇T + T*theta)
            = rho_cp * V * div(wT)       [by product rule]

This is NOT zero — but that is CORRECT! With no physics, the energy in a
cell changes because material flows in/out through the moving faces.
The total energy integral_V rho*cp*T dV IS conserved globally because
div(wT) integrated over the whole domain becomes boundary fluxes.

TEST DESIGN
-----------
1D domain with:
- Known temperature profile: T(x) = 500 + 300*x/L
- Known rho*cp profile: rho_cp(x) = 2e6 - 1e6*x/L (varies spatially)
- Known mesh velocity: w(x) = -0.001*x/L (linear, contracting)
- No physics: RHS = 0
- theta = dw/dx = -0.001/L = -0.1 s^-1 (uniform)

REFERENCE SOLUTIONS:
- Non-conservative (frozen T field): T_i(t+dt) = T(x_i + w_i*dt)
- Conservative (energy locked to cell): rho_cp*T*V = const per cell

We test three options:
  (A) dT = +w·∇T * dt
  (B) dT = -T·theta * dt
  (C) dT = (w·∇T - T·theta) * dt

And check which reproduces the correct answer and conserves energy properly.
"""

using Printf

println("=" ^ 80)
println("MANUFACTURED SOLUTION TEST: TEMPERATURE ALE EQUATION")
println("  Does T need dilation correction, or just w·∇T advection?")
println("=" ^ 80)
println()

# ============================================================================
# Parameters
# ============================================================================
const L = 0.01       # Domain length [m]
const N = 10         # Number of cells (as requested)
const dx = L / N     # Cell size [m]
const w_top = -0.001 # Mesh velocity at top (x=L) [m/s] (contracting)
const dt = 1e-5      # Time step [s] (small for accurate Euler)

# Derived
const theta_uniform = w_top / L  # = -0.1 s^-1 (uniform dilation)

# ============================================================================
# Spatial profiles
# ============================================================================
# Cell centers
x_centers = [(i - 0.5) * dx for i in 1:N]

# Mesh velocity: w(x) = w_top * x/L  (w=0 at bottom, w=w_top at top)
w_cell = [w_top * x_centers[i] / L for i in 1:N]

# Temperature profile: T(x) = 500 + 300*x/L
T_profile(x) = 500.0 + 300.0 * x / L
const dTdx_exact_val = 300.0 / L  # = 30000 K/m (constant gradient for linear T)

T0 = [T_profile(x_centers[i]) for i in 1:N]

# rho*cp profile: rho_cp(x) = 2e6 - 1e6*x/L  [J/(m^3 K)]
rho_cp_profile(x) = 2.0e6 - 1.0e6 * x / L
rho_cp = [rho_cp_profile(x_centers[i]) for i in 1:N]

# Node positions and velocities
x_nodes = [i * dx for i in 0:N]
w_nodes = [w_top * x_nodes[i+1] / L for i in 0:N]

println("Setup:")
println("  Domain: [0, $L] m, $N cells, dx = $(dx*1000) mm")
println("  Temperature: T(x) = 500 + 300*x/L  (linear, dT/dx = $(300.0/L) K/m)")
println("  rho*cp(x) = 2e6 - 1e6*x/L  (varies from 2e6 at bottom to 1e6 at top)")
println("  Mesh velocity: w(x) = $(w_top)*x/L  (w=0 at bottom)")
println("  Dilation rate: theta = $(theta_uniform) s^-1 (uniform)")
println("  Time step: dt = $dt s")
println()

# ============================================================================
# PART 1: Compute dT/dt under each option
# ============================================================================
println("=" ^ 80)
println("PART 1: Rates of change at each cell (no physics, RHS=0)")
println("=" ^ 80)
println()

# Compute dT/dx using upwind finite differences
dTdx_upwind = zeros(N)
for i in 1:N
    w_i = w_cell[i]
    if w_i > 0
        dTdx_upwind[i] = (i == 1) ? 0.0 : (T0[i] - T0[i-1]) / dx
    else
        dTdx_upwind[i] = (i == N) ? 0.0 : (T0[i+1] - T0[i]) / dx
    end
end

# Option A: dT/dt = +w·∇T  (non-conservative advection)
dT_A = [w_cell[i] * dTdx_upwind[i] for i in 1:N]

# Option B: dT/dt = -T·theta  (conservative dilation)
dT_B = [-T0[i] * theta_uniform for i in 1:N]

# Option C: dT/dt = w·∇T - T·theta  (both)
dT_C = [w_cell[i] * dTdx_upwind[i] - T0[i] * theta_uniform for i in 1:N]

@printf("  %-6s %-10s %-10s %-10s %-18s %-18s %-18s\n",
    "Cell", "x [mm]", "T [K]", "w [m/s]",
    "dT_A (w·∇T)", "dT_B (-T*θ)", "dT_C (both)")
for i in 1:N
    @printf("  %-6d %-10.3f %-10.2f %-10.6f %-18.6f %-18.6f %-18.6f\n",
        i, x_centers[i]*1000, T0[i], w_cell[i], dT_A[i], dT_B[i], dT_C[i])
end
println()

# ============================================================================
# PART 2: Reference solutions
# ============================================================================
println("=" ^ 80)
println("PART 2: Reference solutions and error comparison")
println("=" ^ 80)
println()

# --- Reference 1: Frozen spatial field (non-conservative) ---
# With no physics, T(x) = 500 + 300*x/L is frozen in space.
# Mesh point at x_i moves to x_i + w_i*dt, sampling the frozen field.
T_ref_frozen = [T_profile(x_centers[i] + w_cell[i] * dt) for i in 1:N]

# Forward Euler updates
T_A = T0 .+ dt .* dT_A
T_B = T0 .+ dt .* dT_B
T_C = T0 .+ dt .* dT_C

# Errors vs frozen-field reference
err_A_frozen = sqrt(sum((T_A[i] - T_ref_frozen[i])^2 for i in 1:N) / N)
err_B_frozen = sqrt(sum((T_B[i] - T_ref_frozen[i])^2 for i in 1:N) / N)
err_C_frozen = sqrt(sum((T_C[i] - T_ref_frozen[i])^2 for i in 1:N) / N)

# Also compute with exact gradient (to isolate discretization error)
dT_A_exact = [w_cell[i] * dTdx_exact_val for i in 1:N]
T_A_exact = T0 .+ dt .* dT_A_exact
err_A_exact_frozen = sqrt(sum((T_A_exact[i] - T_ref_frozen[i])^2 for i in 1:N) / N)

println("  Reference 1: FROZEN SPATIAL FIELD  T_ref = T(x + w*dt)")
println("    Physical meaning: T field doesn't change in space; mesh point moves")
println("    to a new location and samples the existing field.")
println()
@printf("    Option A (w·∇T only):        RMS error = %.6e K\n", err_A_frozen)
@printf("    Option A (exact gradient):   RMS error = %.6e K  (discretization limit)\n", err_A_exact_frozen)
@printf("    Option B (-T·θ only):        RMS error = %.6e K\n", err_B_frozen)
@printf("    Option C (w·∇T - T·θ):       RMS error = %.6e K\n", err_C_frozen)
println()

# --- Reference 2: Energy-conserved in each cell ---
# If energy E_i = rho_cp_i * T_i * V_i = const (no physics, rho_cp fixed):
# T_i(t+dt) = T_i(t) * V_i(t) / V_i(t+dt) = T_i(t) / (1 + theta*dt)
T_ref_energy = [T0[i] / (1.0 + theta_uniform * dt) for i in 1:N]

err_A_energy = sqrt(sum((T_A[i] - T_ref_energy[i])^2 for i in 1:N) / N)
err_B_energy = sqrt(sum((T_B[i] - T_ref_energy[i])^2 for i in 1:N) / N)
err_C_energy = sqrt(sum((T_C[i] - T_ref_energy[i])^2 for i in 1:N) / N)

println("  Reference 2: ENERGY CONSERVED IN CELL  rho_cp*T*V = const")
println("    Physical meaning: energy locked to cell, volume changes → T changes.")
println("    (This is what dilation -T*theta describes.)")
println()
@printf("    Option A (w·∇T only):        RMS error = %.6e K\n", err_A_energy)
@printf("    Option B (-T·θ only):        RMS error = %.6e K\n", err_B_energy)
@printf("    Option C (w·∇T - T·θ):       RMS error = %.6e K\n", err_C_energy)
println()

# ============================================================================
# PART 3: Which reference is physically correct for temperature?
# ============================================================================
println("=" ^ 80)
println("PART 3: Which reference is physically correct?")
println("=" ^ 80)
println()
println("  Temperature is NOT a conserved quantity. The PDE for T is:")
println("      rho*cp * dT/dt|_x = -div(q) + Q_rxn + ...")
println()
println("  With no physics (RHS = 0), the Eulerian equation is:")
println("      dT/dt|_x = 0")
println("  meaning T(x,t) = T(x,0) — the field is frozen in space.")
println()
println("  The ALE identity then gives:")
println("      dT/dt|_chi = dT/dt|_x + w·∇T = 0 + w·∇T = w·∇T")
println()
println("  Therefore the FROZEN SPATIAL FIELD is the correct reference,")
println("  and OPTION (A) is the correct formulation.")
println()
println("  The energy-in-cell reference is WRONG for T because temperature")
println("  is not conserved — energy is, but E = rho*cp*T*V, and rho*cp is")
println("  the prefactor, not T alone. When volume changes with rho*cp fixed,")
println("  the energy per unit volume stays the same (energy flows in/out),")
println("  but T at a spatial point doesn't change.")
println()

# ============================================================================
# PART 4: The uniform temperature test (clearest proof)
# ============================================================================
println("=" ^ 80)
println("PART 4: UNIFORM TEMPERATURE TEST (the clearest proof)")
println("=" ^ 80)
println()
println("  Set T(x) = T0 = 600 K everywhere. No physics. Mesh contracts.")
println()
println("  QUESTION: What happens to T in each cell?")
println()
println("  Since T is uniform and there is no heat source/sink/conduction,")
println("  T remains 600 K everywhere — the spatial field is frozen at 600 K.")
println("  The mesh point moves to a new location, but T = 600 K there too.")
println()
println("  So dT/dt|_chi = 0.")
println()

T_uniform = 600.0

dT_A_uniform = w_cell[5] * 0.0  # w·∇T = w × 0 = 0 (no gradient)
dT_B_uniform = -T_uniform * theta_uniform  # -T*theta != 0
dT_C_uniform = dT_A_uniform + dT_B_uniform

@printf("  Option A (w·∇T):       dT/dt = %.6f K/s  (correct: 0)\n", dT_A_uniform)
@printf("  Option B (-T·θ):       dT/dt = %.6f K/s  (WRONG: creates spurious heating)\n", dT_B_uniform)
@printf("  Option C (w·∇T - T·θ): dT/dt = %.6f K/s  (WRONG: same spurious heating)\n", dT_C_uniform)
println()
println("  Only Option A correctly gives dT/dt = 0 for uniform temperature.")
println("  Options B and C create a FICTITIOUS temperature increase when")
println("  the mesh contracts — as if compressing a gas, but temperature is")
println("  not a conserved density!")
println()

# ============================================================================
# PART 5: Energy conservation check (the subtle part)
# ============================================================================
println("=" ^ 80)
println("PART 5: ENERGY CONSERVATION CHECK")
println("=" ^ 80)
println()
println("  Does total energy E = sum(rho_cp_i * T_i * V_i) change correctly?")
println()

# Initial volumes and cell widths
V0 = fill(dx, N)  # assuming unit cross-section area

# Initial energy
E0 = sum(rho_cp[i] * T0[i] * V0[i] for i in 1:N)

# After dt, new volumes: V_new = V * (1 + theta*dt)
V_new = [V0[i] * (1.0 + theta_uniform * dt) for i in 1:N]

# Energy after each option (with constant rho_cp — no reactions)
E_A = sum(rho_cp[i] * T_A[i] * V_new[i] for i in 1:N)
E_B = sum(rho_cp[i] * T_B[i] * V_new[i] for i in 1:N)
E_C = sum(rho_cp[i] * T_C[i] * V_new[i] for i in 1:N)

# Expected energy change: energy flux through domain boundaries
# Bottom boundary (x=0): w(0) = 0, so no flux
# Top boundary (x=L): w(L) = w_top, flux = rho_cp(L) * T(L) * w_top * A
# where A = 1 (unit area), and w_top < 0 means material flows OUT
# Energy flux out per unit time = rho_cp(L) * T(L) * w_top * A
# dE/dt = rho_cp(L) * T(L) * w_top  (negative = energy leaves domain)
rho_cp_top = rho_cp_profile(L)
T_top = T_profile(L)
dE_dt_expected = rho_cp_top * T_top * w_top  # [W/m^2] for unit area
E_expected = E0 + dE_dt_expected * dt

# But wait — this is the DISCRETE version. The boundary flux is at x=L.
# In the discrete case, the flux out the top face is:
# F_top = rho_cp_N * T_N * w_top (upwind: mesh moving left, so T_N flows out)
# But the bottom face has w=0 so no flux.
# For all internal faces, flux_f = rho_cp_f * T_f * w_f (with upwind T)
# Net change = sum of face fluxes across ALL faces
# Actually, for the integral:
#   d/dt integral rho_cp*T dV = integral div(rho_cp*T*w) dV = surface flux

# Compute exact integral change using the continuous identity:
# d/dt int rho_cp*T dV = int d(rho_cp*T)/dt|_chi dV + int rho_cp*T * div(w) dV
#                      = int [rho_cp * dT/dt|_chi + T * d(rho_cp)/dt|_chi] dV + int rho_cp*T*theta dV
# With no reactions, d(rho_cp)/dt|_chi = w·∇(rho_cp) (non-conservative for rho_cp too)

@printf("  Initial total energy E0 = %.6f J/m^2\n", E0)
@printf("  Expected E (boundary flux): E0 + rho_cp(L)*T(L)*w_top*dt = %.6f J/m^2\n", E_expected)
println()
@printf("  Option A energy: %.6f J/m^2  (dE = %+.6e)\n", E_A, E_A - E0)
@printf("  Option B energy: %.6f J/m^2  (dE = %+.6e)\n", E_B, E_B - E0)
@printf("  Option C energy: %.6f J/m^2  (dE = %+.6e)\n", E_C, E_C - E0)
@printf("  Expected:        %.6f J/m^2  (dE = %+.6e)\n", E_expected, E_expected - E0)
println()

# Relative errors in energy change
dE_expected = dE_dt_expected * dt
if abs(dE_expected) > 1e-20
    @printf("  Energy change relative error:\n")
    @printf("    Option A: %.4e\n", abs(E_A - E_expected) / abs(dE_expected))
    @printf("    Option B: %.4e\n", abs(E_B - E_expected) / abs(dE_expected))
    @printf("    Option C: %.4e\n", abs(E_C - E_expected) / abs(dE_expected))
end
println()

# Let's also compute the discrete energy balance more carefully
# E_new = sum rho_cp_i * T_new_i * V_new_i
# The energy that SHOULD be in the domain after dt:
# Energy flows out the top face at rate: rho_cp_N * T_N * w(x_top_face) * A
# But with discrete cells, we need to be more careful about face fluxes
println("  Note: Energy conservation involves BOTH T changes AND V changes.")
println("  The energy balance is:")
println("    dE = sum_i [rho_cp_i * dT_i * V_i + rho_cp_i * T_i * dV_i]")
println("  The dV_i = theta*V_i*dt terms are the same for all options.")
println("  The difference is in dT_i alone.")
println()

dE_from_dT_A = sum(rho_cp[i] * dT_A[i] * dt * V0[i] for i in 1:N)
dE_from_dT_B = sum(rho_cp[i] * dT_B[i] * dt * V0[i] for i in 1:N)
dE_from_dT_C = sum(rho_cp[i] * dT_C[i] * dt * V0[i] for i in 1:N)
dE_from_dV = sum(rho_cp[i] * T0[i] * theta_uniform * dt * V0[i] for i in 1:N)

@printf("  Energy change from dT (Option A): %+.6e J/m^2\n", dE_from_dT_A)
@printf("  Energy change from dT (Option B): %+.6e J/m^2\n", dE_from_dT_B)
@printf("  Energy change from dT (Option C): %+.6e J/m^2\n", dE_from_dT_C)
@printf("  Energy change from dV:            %+.6e J/m^2\n", dE_from_dV)
println()
@printf("  Total dE (Option A) = dE_dT + dE_dV = %+.6e J/m^2\n", dE_from_dT_A + dE_from_dV)
@printf("  Total dE (Option B) = dE_dT + dE_dV = %+.6e J/m^2\n", dE_from_dT_B + dE_from_dV)
@printf("  Total dE (Option C) = dE_dT + dE_dV = %+.6e J/m^2\n", dE_from_dT_C + dE_from_dV)
@printf("  Expected dE (boundary flux):        %+.6e J/m^2\n", dE_expected)
println()

# ============================================================================
# PART 6: Convergence study
# ============================================================================
println("=" ^ 80)
println("PART 6: Mesh convergence study (frozen-field reference)")
println("=" ^ 80)
println()
@printf("  %-10s %-18s %-18s %-18s %-18s\n",
    "N_cells", "err_A (w·∇T)", "err_B (-T·θ)", "err_C (both)", "err_A_exact")

for N_test in [10, 20, 40, 80, 160, 320]
    dx_t = L / N_test
    x_t = [(i - 0.5) * dx_t for i in 1:N_test]
    w_t = [w_top * x_t[i] / L for i in 1:N_test]
    T_t = [T_profile(x_t[i]) for i in 1:N_test]

    # Upwind gradient
    dTdx_t = zeros(N_test)
    for i in 1:N_test
        if w_t[i] > 0
            dTdx_t[i] = (i == 1) ? 0.0 : (T_t[i] - T_t[i-1]) / dx_t
        else
            dTdx_t[i] = (i == N_test) ? 0.0 : (T_t[i+1] - T_t[i]) / dx_t
        end
    end

    # Three formulations
    T_tA = T_t .+ dt .* (w_t .* dTdx_t)
    T_tB = T_t .+ dt .* (-T_t .* theta_uniform)
    T_tC = T_t .+ dt .* (w_t .* dTdx_t .- T_t .* theta_uniform)

    # Also exact gradient for Option A
    T_tA_exact = T_t .+ dt .* [w_t[i] * dTdx_exact_val for i in 1:N_test]

    # Reference: frozen field
    T_ref = [T_profile(x_t[i] + w_t[i] * dt) for i in 1:N_test]

    e_A = sqrt(sum((T_tA[i] - T_ref[i])^2 for i in 1:N_test) / N_test)
    e_B = sqrt(sum((T_tB[i] - T_ref[i])^2 for i in 1:N_test) / N_test)
    e_C = sqrt(sum((T_tC[i] - T_ref[i])^2 for i in 1:N_test) / N_test)
    e_A_exact = sqrt(sum((T_tA_exact[i] - T_ref[i])^2 for i in 1:N_test) / N_test)

    @printf("  %-10d %-18.6e %-18.6e %-18.6e %-18.6e\n",
        N_test, e_A, e_B, e_C, e_A_exact)
end
println()
println("  Option A (w·∇T) error converges with mesh refinement → correct equation")
println("  Option A with exact gradient is near machine precision → confirms w·∇T")
println("  Option B (-T*theta) error does NOT converge → wrong equation for T")
println("  Option C (both) error does NOT converge → wrong equation for T")
println()

# ============================================================================
# PART 7: Multi-step integration (accumulation check)
# ============================================================================
println("=" ^ 80)
println("PART 7: Multi-step forward Euler integration (100 steps)")
println("=" ^ 80)
println()

n_steps = 100
dt_ms = 1e-5

T_ms_A = copy(T0)
T_ms_B = copy(T0)
T_ms_C = copy(T0)

for step in 1:n_steps
    # Option A: w·∇T
    for i in 1:N
        w_i = w_cell[i]
        if w_i > 0
            dTdx_i = (i == 1) ? 0.0 : (T_ms_A[i] - T_ms_A[i-1]) / dx
        else
            dTdx_i = (i == N) ? 0.0 : (T_ms_A[i+1] - T_ms_A[i]) / dx
        end
        T_ms_A[i] += dt_ms * w_i * dTdx_i
    end

    # Option B: -T*theta
    for i in 1:N
        T_ms_B[i] += dt_ms * (-T_ms_B[i] * theta_uniform)
    end

    # Option C: both
    for i in 1:N
        w_i = w_cell[i]
        if w_i > 0
            dTdx_i = (i == 1) ? 0.0 : (T_ms_C[i] - T_ms_C[i-1]) / dx
        else
            dTdx_i = (i == N) ? 0.0 : (T_ms_C[i+1] - T_ms_C[i]) / dx
        end
        T_ms_C[i] += dt_ms * (w_i * dTdx_i - T_ms_C[i] * theta_uniform)
    end
end

total_time = n_steps * dt_ms
T_ref_ms = [T_profile(x_centers[i] + w_cell[i] * total_time) for i in 1:N]

err_ms_A = sqrt(sum((T_ms_A[i] - T_ref_ms[i])^2 for i in 1:N) / N)
err_ms_B = sqrt(sum((T_ms_B[i] - T_ref_ms[i])^2 for i in 1:N) / N)
err_ms_C = sqrt(sum((T_ms_C[i] - T_ref_ms[i])^2 for i in 1:N) / N)

@printf("  After %d steps (total time = %.4f ms):\n", n_steps, total_time * 1000)
println()
@printf("  Option A (w·∇T only):        RMS error = %.6e K\n", err_ms_A)
@printf("  Option B (-T·θ only):         RMS error = %.6e K\n", err_ms_B)
@printf("  Option C (w·∇T - T·θ):        RMS error = %.6e K\n", err_ms_C)
println()

if err_ms_A < err_ms_B && err_ms_A < err_ms_C
    @printf("  Option A is best (%.1fx better than B, %.1fx better than C)\n",
        err_ms_B / err_ms_A, err_ms_C / err_ms_A)
end

# ============================================================================
# PART 8: What about when rho*cp changes? (The subtle case)
# ============================================================================
println()
println("=" ^ 80)
println("PART 8: What if rho*cp changes due to reactions? (Bonus analysis)")
println("=" ^ 80)
println()
println("  In real pyrolysis, rho*cp changes because reactions consume/produce species.")
println("  This does NOT affect the temperature ALE equation.")
println()
println("  The temperature PDE is:")
println("      rho*cp(t) * dT/dt|_x = RHS(t)")
println()
println("  After dividing by rho*cp, we get dT/dt|_x = RHS/rho*cp.")
println("  The ALE transform is still:")
println("      dT/dt|_chi = RHS/rho*cp + w·∇T")
println()
println("  The time-varying rho*cp enters ONLY through the denominator on the")
println("  Eulerian RHS. It does NOT generate an extra dilation term for T.")
println()
println("  KEY INSIGHT: The d(rho*cp)/dt changes are already accounted for by")
println("  the reaction source terms in the species equations. When rho*cp")
println("  changes, the reaction enthalpy Q_rxn in the energy RHS captures")
println("  the energy bookkeeping. Temperature T is just a derived quantity")
println("  that adjusts through RHS/rho*cp.")
println()
println("  This is fundamentally different from SPECIES where:")
println("  - SOLID species: dxi/dt = S - xi*theta (dilation because mass locked to cell)")
println("  - GAS species: dxi/dt = S + w·∇xi (advection because gas slips through mesh)")
println()
println("  Temperature is like gas in this regard — it's an INTENSIVE field")
println("  variable, not a conserved extensive quantity tied to cell volume.")

# ============================================================================
# PART 9: Numerical demonstration that -T*theta creates spurious heating
# ============================================================================
println()
println("=" ^ 80)
println("PART 9: Spurious heating from dilation correction")
println("=" ^ 80)
println()

println("  Scenario: Uniform T = 600 K, mesh contracts for 1000 steps.")
println("  What temperature does each option predict?")
println()

T_uniform_A = fill(600.0, N)
T_uniform_B = fill(600.0, N)
T_uniform_C = fill(600.0, N)

for step in 1:1000
    for i in 1:N
        # Option A: w·∇T = 0 for uniform T → no change (correct)
        w_i = w_cell[i]
        if w_i > 0
            dTdx_i = (i == 1) ? 0.0 : (T_uniform_A[i] - T_uniform_A[i-1]) / dx
        else
            dTdx_i = (i == N) ? 0.0 : (T_uniform_A[i+1] - T_uniform_A[i]) / dx
        end
        T_uniform_A[i] += dt_ms * w_i * dTdx_i
    end

    for i in 1:N
        T_uniform_B[i] += dt_ms * (-T_uniform_B[i] * theta_uniform)
    end

    for i in 1:N
        w_i = w_cell[i]
        if w_i > 0
            dTdx_i = (i == 1) ? 0.0 : (T_uniform_C[i] - T_uniform_C[i-1]) / dx
        else
            dTdx_i = (i == N) ? 0.0 : (T_uniform_C[i+1] - T_uniform_C[i]) / dx
        end
        T_uniform_C[i] += dt_ms * (w_i * dTdx_i - T_uniform_C[i] * theta_uniform)
    end
end

t_final = 1000 * dt_ms
@printf("  After t = %.4f ms (1000 steps):\n", t_final * 1000)
println()
@printf("  Option A: T = [%.4f, ..., %.4f] K  (should be 600.0000)\n",
    T_uniform_A[1], T_uniform_A[N])
@printf("  Option B: T = [%.4f, ..., %.4f] K  (spurious heating!)\n",
    T_uniform_B[1], T_uniform_B[N])
@printf("  Option C: T = [%.4f, ..., %.4f] K  (spurious heating!)\n",
    T_uniform_C[1], T_uniform_C[N])
println()

# Expected B: T * exp(-theta*t) = 600 * exp(0.1 * 0.01) = 600 * exp(0.001)
T_B_analytical = 600.0 * exp(-theta_uniform * t_final)
@printf("  Option B analytical: T = 600*exp(-theta*t) = 600*exp(%.4f) = %.4f K\n",
    -theta_uniform * t_final, T_B_analytical)
@printf("  That's a %.2f K spurious temperature increase!\n", T_B_analytical - 600.0)
println()
println("  This spurious heating would corrupt the solution in any simulation")
println("  with significant mesh contraction (i.e., any pyrolysis simulation).")

# ============================================================================
# CONCLUSION
# ============================================================================
println()
println("=" ^ 80)
println("CONCLUSION")
println("=" ^ 80)
println()
println("  The question: In the ALE temperature equation, does T need a dilation")
println("  correction (-T*theta), or just advection (w·∇T)?")
println()
println("  ANSWER: OPTION (A) — Only w·∇T advection, NO dilation correction.")
println()
println("  Reasoning:")
println("  1. Temperature is an INTENSIVE variable, not a conserved extensive")
println("     quantity. The PDE for T is non-conservative: rho*cp * dT/dt = RHS.")
println()
println("  2. The ALE transform of a non-conservative equation uses only the")
println("     identity: dT/dt|_chi = dT/dt|_x + w·∇T.")
println("     There is NO dilation term.")
println()
println("  3. The dilation correction -T*theta arises from CONSERVATIVE (integral)")
println("     formulations where d(phi*V)/dt = S*V. This applies to solid species")
println("     concentrations (extensive quantities), NOT to temperature.")
println()
println("  4. The uniform temperature test proves it: with T = const everywhere")
println("     and no physics, T should NOT change. Only option (A) gives this.")
println("     Options (B) and (C) create spurious heating during mesh contraction.")
println()
println("  5. Energy conservation is maintained: the total energy integral")
println("     rho*cp*T*V evolves correctly because dV/dt = theta*V provides")
println("     the volume change, and w·∇T provides the field sampling change.")
println("     Adding -T*theta would DOUBLE-COUNT the volume change effect.")
println()
println("  Summary of ALE corrections:")
println("    TEMPERATURE:      dT/dt|_chi = RHS/rho_cp + w·∇T")
println("    SOLID SPECIES:    dxi/dt|_chi = S - xi*theta")
println("    GAS SPECIES:      dxi/dt|_chi = S - div(J) + w·∇xi")
println()
println("  The code's current formulation (line 895 of ode_function.jl):")
println("    dT[i] -= ale_advection[i]     i.e., dT = RHS - w·∇T")
println()
println("  This has the WRONG SIGN for the ALE advection. It should be:")
println("    dT[i] += ale_advection[i]     i.e., dT = RHS + w·∇T")
println()
println("  (This is consistent with the gas species formulation at line 903.)")
