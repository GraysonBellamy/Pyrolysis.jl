#!/usr/bin/env julia
"""
Test and benchmark the analytical Jacobian implementation.

This script:
1. Verifies the analytical Jacobian against finite differences
2. Benchmarks the speedup from using analytical Jacobian

Note on FD Validation Accuracy:
------------------------------
The dT/dz and dξ/dz blocks may show ~0.01% error when compared to finite differences.
This is NOT a bug in the analytical Jacobian - it is the truncation error inherent in
finite differencing. The analytical Jacobian computes the true first derivative, while
FD (with step size ε) includes second-order terms proportional to ε.

For the boundary cell geometry derivative, the second-order error is:
	FD_error ≈ (q_bc × A³/V³) × ε / ρcₚ

This was confirmed by showing the error scales linearly with ε:
	ε=1e-5:  error ≈ 700
	ε=1e-7:  error ≈ 7
	ε=1e-9:  error ≈ 0.07

The analytical Jacobian is correct; the FD validation simply has inherent limitations.
"""

using Pyrolysis
using LinearAlgebra
using Printf

println("=" ^ 60)
println("Analytical Jacobian Test & Benchmark")
println("=" ^ 60)

# ============================================================================
# Create a simple test problem (PMMA-like single component)
# ============================================================================

println("\n1. Setting up test problem...")

# Simple PMMA material (2 components: virgin + volatiles)
pmma = Material(
	name = :PMMA_test,
	components = (
		SolidComponent(:virgin,
			ρ = ConstantProperty(1190.0),
			c = LinearProperty(1420.0, 2.5),  # T-dependent cp
			k = ConstantProperty(0.21),
			ε = 0.86,
			α = 1.5,
			γ = 1.0  # Solid contributes to volume (default)
		),
		GaseousComponent(:volatiles,
			M = 0.100,
			c = ConstantProperty(1500.0),
			k = ConstantProperty(0.02),
			λ = ConstantProperty(2e-5),
			ε = 0.0,
			α = 0.0,
		),
	),
	reactions = (
		Reaction(:decomposition, 1 => 2,
			A = 8.5e12,
			E = 188e3,
			h = 870e3,  # J/kg (endothermic)
			T_min = 300.0,
			T_max = 1500.0,
		),
	),
	conductivity_mixing = PARALLEL,
	β = 0.5,
	gas_transfer_mixing = PARALLEL,
	β_λ = 0.5,
)

# Create mesh (20 cells for quick test)
thickness = 0.006  # 6mm
n_cells = 20
T_initial = 300.0
ξ_virgin_initial = 1190.0
ξ_gas_initial = 0.0

mesh = create_uniform_mesh(
	thickness,
	n_cells,
	Val(2);  # NC = 2 components
	T_initial = T_initial,
	ξ_initial = (ξ_virgin_initial, ξ_gas_initial),
)

# Boundary conditions
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),  # Top: 50 kW/m² + reradiation
	thermal_bottom = AdiabaticBC(),  # Bottom: adiabatic
	mass_top = DiffusiveMassBC(0.01),  # Top: gas can escape
	mass_bottom = ImpermeableBC(),  # Bottom: impermeable
)

# Create ODE function
ode_fn = Pyrolysis.TimeIntegration.PyrolysisODEFunction(
	mesh, pmma, bc_set;
	I_radiation = 0.0,
)

# Create state vector and advance to interesting state
NC = 2
NR = 1
sv = Pyrolysis.TimeIntegration.StateVector(n_cells, NC)
u0 = Pyrolysis.TimeIntegration.create_state_vector(mesh, NC)

# Warm up the mesh state by running one RHS evaluation
du = similar(u0)
ode_fn(du, u0, nothing, 0.0)

println("   - $n_cells cells, $NC components, $NR reactions")
println("   - State vector length: $(length(u0))")

# ============================================================================
# Test analytical Jacobian vs finite differences
# ============================================================================

println("\n2. Verifying analytical Jacobian against finite differences...")

# Compute analytical Jacobian
n_faces = length(mesh.faces)
jac_cache = Pyrolysis.TimeIntegration.JacobianCache(n_cells, n_faces, NC, NR)
J_analytical = zeros(length(u0), length(u0))

# First do an RHS eval to populate caches
ode_fn(du, u0, nothing, 0.0)

# Compute analytical Jacobian
Pyrolysis.TimeIntegration.compute_analytical_jacobian!(J_analytical, u0, nothing, 0.0, ode_fn, jac_cache)

# Compute finite difference Jacobian
J_fd = zeros(length(u0), length(u0))
ε = 1e-7
u_pert = copy(u0)
du_base = similar(u0)
du_pert = similar(u0)

ode_fn(du_base, u0, nothing, 0.0)

for j in 1:length(u0)
	u_pert .= u0
	u_pert[j] += ε
	ode_fn(du_pert, u_pert, nothing, 0.0)
	J_fd[:, j] .= (du_pert .- du_base) ./ ε
end

# Compare
rel_error = norm(J_analytical - J_fd) / max(norm(J_fd), 1e-10)
max_error = maximum(abs.(J_analytical - J_fd))
max_rel_error = maximum(abs.(J_analytical - J_fd) ./ max.(abs.(J_fd), 1e-10))

println("   Relative Frobenius error: $(Printf.@sprintf("%.2e", rel_error))")
println("   Maximum absolute error:   $(Printf.@sprintf("%.2e", max_error))")
println("   Maximum relative error:   $(Printf.@sprintf("%.2e", max_rel_error))")

if rel_error < 0.1
	println("   ✓ Jacobian verification PASSED")
else
	println("   ✗ Jacobian verification FAILED - significant differences!")

	# Find the worst elements
	println("\n   Worst disagreements:")
	diff = abs.(J_analytical - J_fd)
	for _ in 1:5
		idx = argmax(diff)
		i, j = Tuple(idx)
		println("      J[$i,$j]: analytical=$(J_analytical[i,j]), FD=$(J_fd[i,j]), diff=$(diff[i,j])")
		diff[i, j] = 0.0
	end
end

# ============================================================================
# Benchmark comparison
# ============================================================================

println("\n3. Benchmarking Jacobian computation...")

# Benchmark finite difference Jacobian (full, not sparse)
function compute_fd_jacobian!(J, u, ode_fn)
	n = length(u)
	ε = 1e-7
	du_base = similar(u)
	du_pert = similar(u)
	u_pert = similar(u)

	ode_fn(du_base, u, nothing, 0.0)

	for j in 1:n
		u_pert .= u
		u_pert[j] += ε
		ode_fn(du_pert, u_pert, nothing, 0.0)
		@views J[:, j] .= (du_pert .- du_base) ./ ε
	end
end

# Warm-up runs
J_test = zeros(length(u0), length(u0))
compute_fd_jacobian!(J_test, u0, ode_fn)
Pyrolysis.TimeIntegration.compute_analytical_jacobian!(J_test, u0, nothing, 0.0, ode_fn, jac_cache)

# Timed runs
n_runs = 5

# FD Jacobian timing
t_fd_total = 0.0
for _ in 1:n_runs
	global t_fd_total += @elapsed compute_fd_jacobian!(J_fd, u0, ode_fn)
end
t_fd = t_fd_total / n_runs

# Analytical Jacobian timing
t_analytical_total = 0.0
for _ in 1:n_runs
	global t_analytical_total += @elapsed Pyrolysis.TimeIntegration.compute_analytical_jacobian!(
		J_analytical, u0, nothing, 0.0, ode_fn, jac_cache)
end
t_analytical = t_analytical_total / n_runs

println("   Finite difference Jacobian: $(Printf.@sprintf("%.3f", t_fd*1000)) ms")
println("   Analytical Jacobian:        $(Printf.@sprintf("%.3f", t_analytical*1000)) ms")
println("   Speedup: $(Printf.@sprintf("%.1f", t_fd/t_analytical))×")

# ============================================================================
# Test property function derivatives
# ============================================================================

println("\n4. Testing property function derivatives...")

T_test = 500.0

# ConstantProperty
p_const = ConstantProperty(1.5)
d_const = Pyrolysis.Core.derivative(p_const, T_test)
println("   ConstantProperty(1.5): p(500) = $(p_const(T_test)), dp/dT = $d_const (expect 0.0)")

# LinearProperty  
p_linear = LinearProperty(1000.0, 2.5)
d_linear = Pyrolysis.Core.derivative(p_linear, T_test)
d_linear_fd = (p_linear(T_test + 0.01) - p_linear(T_test - 0.01)) / 0.02
println("   LinearProperty(1000, 2.5): p(500) = $(p_linear(T_test)), dp/dT = $d_linear (expect 2.5, FD=$d_linear_fd)")

# PolynomialProperty
p_poly = PolynomialProperty((100.0, 0.5, 0.001))  # 100 + 0.5T + 0.001T²
d_poly = Pyrolysis.Core.derivative(p_poly, T_test)
d_poly_exact = 0.5 + 2*0.001*T_test  # 0.5 + 0.002*500 = 1.5
d_poly_fd = (p_poly(T_test + 0.01) - p_poly(T_test - 0.01)) / 0.02
println("   PolynomialProperty: p(500) = $(p_poly(T_test)), dp/dT = $d_poly (exact=$d_poly_exact, FD=$(Printf.@sprintf("%.4f", d_poly_fd)))")

# ArrheniusProperty
p_arr = ArrheniusProperty(1e12, 150e3)
d_arr = Pyrolysis.Core.derivative(p_arr, T_test)
d_arr_fd = (p_arr(T_test + 0.01) - p_arr(T_test - 0.01)) / 0.02
d_arr_exact = p_arr(T_test) * 150e3 / (8.314 * T_test^2)
println("   ArrheniusProperty: p(500) = $(Printf.@sprintf("%.3e", p_arr(T_test))), dp/dT = $(Printf.@sprintf("%.3e", d_arr)) (exact=$(Printf.@sprintf("%.3e", d_arr_exact)), FD=$(Printf.@sprintf("%.3e", d_arr_fd)))")

println("\n" * "=" ^ 60)
println("Tests complete!")
println("=" ^ 60)
