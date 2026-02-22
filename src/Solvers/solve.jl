"""
Main solve interface for PyrolysisProblem.

Reference: Design Document Section 7.1, 12.1
"""

using ..Problem: PyrolysisProblem, PyrolysisSolution, validate_problem,
	apply_initial_conditions!, set_initial_conditions!
using ..TimeIntegration: PyrolysisODEFunction, PyrolysisODEFunctionWithRadiation,
	StateVector, create_state_vector, unpack_state!,
	# ALE types
	PyrolysisODEFunctionALE, PyrolysisODEFunctionALEWithRadiation,
	StateVectorALE, create_state_vector_ale,
	get_node_positions, get_temperature_ale, get_concentration_ale,
	compute_volumes_from_state, get_thickness_from_state,
	# P1 radiation
	PyrolysisODEFunctionALEWithP1Radiation,
	# Sparse Jacobian support
	jacobian_sparsity, jacobian_sparsity_ale, matrix_colors,
	# Analytical Jacobian support
	JacobianCache, JacobianCacheALE, compute_analytical_jacobian!, compute_analytical_jacobian_ale!
using ..Geometry: cell_center
using ..Physics: surface_mass_flux, evaluate_reactions!
using ..Properties: RadiationCache, effective_emissivity, effective_conductivity,
	effective_gas_transfer
using ..Discretization: compute_boundary_heat_fluxes, compute_surface_mass_flux,
	solve_boundary_temperature
using ..Core: get_radiative_flux_bc, extract_radiative_flux_from_bcs, has_radiative_flux_bc,
	RadiativeBC, CombinedThermalBC, AbstractThermalBC, AbstractRadiationBC,
	MarshakRadiationBC, ReflectiveRadiationBC, RadiationBCSet, GAS,
	get_thermal_bc, get_mass_bc
using LinearSolve: LUFactorization, KLUFactorization  # LU for MKL workaround, KLU for sparse Jacobian
using SparseArrays: nnz  # For logging sparse Jacobian info

import CommonSolve: solve  # Extend CommonSolve.solve for our problem type

# ============================================================================
# Main Solve Function
# ============================================================================

"""
	solve(problem::PyrolysisProblem; kwargs...) -> PyrolysisSolution

Solve a pyrolysis problem using DifferentialEquations.jl.

# Arguments
- `problem`: The PyrolysisProblem to solve

# Keyword Arguments
- `solver=Rodas5P()`: ODE solver to use (see DifferentialEquations.jl)
- `abstol=1e-8`: Absolute tolerance
- `reltol=1e-6`: Relative tolerance
- `maxiters=1000000`: Maximum iterations
- `dtmin=1e-12`: Minimum timestep
- `dtmax=Inf`: Maximum timestep
- `saveat=problem.output_times`: Times to save output
- `callback=nothing`: Additional callbacks
- `progress=true`: Show progress bar during solve
- `validate=true`: Validate problem before solving
- `use_ale=false`: Enable ALE mesh motion for shrinking/swelling materials
- `min_thickness=1e-6`: Minimum thickness before termination (ALE only) [m]
- `handle_depletion=false`: Enable surface cell depletion handling (ALE only)
- `depletion_threshold=0.05`: Thickness fraction threshold for cell merging (0.05 = 5% of initial)
- `min_cells=2`: Minimum number of cells before terminating instead of merging

# Returns
A `PyrolysisSolution` containing:
- Temperature and concentration profiles at saved times
- Mass loss rate history
- Surface temperature history
- Thickness history (when use_ale=true)
- Solver return code and timing

# ALE Mode
When `use_ale=true`, the solver uses `PyrolysisODEFunctionALE` which integrates
mesh motion directly into the ODE function. This ensures that node positions
evolve simultaneously with temperature and concentration fields, maintaining
proper conservation (GCL compliance). The thickness history is extracted from
the state vector and stored in `solution.thickness_history`.

# Surface Cell Depletion Handling
When `handle_depletion=true` (requires `use_ale=true`), the solver monitors
the surface cell thickness. When the cell shrinks below `depletion_threshold`
fraction of its initial thickness, the surface cell is merged with its 
interior neighbor. This prevents numerical issues from extremely thin cells:
This prevents the ρcₚ→0 singularity that causes temperature blow-up in materials
that fully decompose. The mesh is dynamically resized during integration.

# Recommended Solvers
- `Rodas5P()`: Robust implicit solver, good default (stiff problems) - **fastest**
- `KenCarp4()`: IMEX solver, similar performance to Rodas5P
- `QNDF()`: BDF-based, good for very stiff problems (slower than Rodas5P)
- `Tsit5()`: Explicit, for non-stiff/mildly stiff problems
- `TRBDF2()`: Implicit, efficient for moderately stiff

# Sparse Jacobian (Performance Note)
The `use_sparse_jacobian=true` option uses matrix coloring to reduce finite 
difference evaluations (~40× fewer calls). However, the sparse assembly overhead
often exceeds the savings for typical problem sizes. Benchmarks show:

- Our Douglas Fir problem (1801 DOFs): Dense is **faster** than sparse
- Sparse FD assembly overhead: ~100ms vs 66ms saved from fewer RHS calls
- Sparse LU is 2× faster, but doesn't compensate for assembly overhead

Sparse Jacobian may help for:
- Very large problems (>5000 DOFs)
- Expensive RHS (>1ms per call)
- Very sparse systems (<0.5% fill)

For most pyrolysis problems, the default dense Jacobian is optimal.

# Example
```julia
problem = PyrolysisProblem(mesh, material, bcs, (0.0, 100.0))
solution = solve(problem)
solution = solve(problem, solver=QNDF(), abstol=1e-10)
solution = solve(problem, progress=false)  # Disable progress bar

# With ALE mesh motion for shrinking materials:
solution = solve(problem, use_ale=true, min_thickness=0.1e-3)
thickness_history = solution.thickness_history

# With surface cell depletion handling (prevents blow-up):
solution = solve(problem, use_ale=true, handle_depletion=true, depletion_threshold=0.01)

# Sparse Jacobian (only for very large problems):
solution = solve(problem, use_sparse_jacobian=true)
```
"""

# ============================================================================
# Radiation Model Validation
# ============================================================================

"""
	_validate_radiation_model(model::RadiationModel, material, radiation_bcs)

Validate radiation model selection and warn about ignored parameters.

Issues warnings when:
- SURFACE_ABSORPTION is used with components that have α values (will be ignored)
- BEER_LAMBERT is used without any α values (no absorption will occur)
- P1_QUASI_STEADY is used without explicit radiation_bcs (defaults will be used)
"""
function _validate_radiation_model(model::RadiationModel, material, radiation_bcs)
	# Warn if SURFACE_ABSORPTION is used with non-zero absorption coefficients
	if model == SURFACE_ABSORPTION
		if has_absorption_coefficients(material.components)
			@warn """
			SURFACE_ABSORPTION radiation model selected, but material components have
			non-zero absorption coefficients (α). These will be IGNORED.

			With SURFACE_ABSORPTION, all incident radiation is deposited in the surface
			cell regardless of material optical properties.

			If you want in-depth absorption using the specified α values, use:
				radiation_model = BEER_LAMBERT
			""" maxlog = 1
		end
	end

	# Warn if BEER_LAMBERT is used without absorption coefficients
	if model == BEER_LAMBERT
		if !has_absorption_coefficients(material.components)
			@warn """
			BEER_LAMBERT radiation model selected, but no components have absorption
			coefficients (α) set. All radiation will pass through without absorption.

			Set absorption coefficients in your material definition:
				SolidComponent(:name, ρ=..., c=..., k=..., α=1.47)  # m²/kg

			Or use SURFACE_ABSORPTION for opaque materials:
				radiation_model = SURFACE_ABSORPTION
			""" maxlog = 1
		end
	end

	# Warn if P1 has no explicit radiation boundary conditions
	if model == P1_QUASI_STEADY
		if radiation_bcs === nothing
			@warn """
			P1_QUASI_STEADY radiation model selected without explicit radiation_bcs.
			Default radiation boundary conditions will be created from thermal BCs.

			For full control, specify radiation boundary conditions:
				radiation_bcs = RadiationBCSet(
					top = MarshakRadiationBC(I_external = 50e3, T_ambient = 300.0),
					bottom = ReflectiveRadiationBC()
				)
			""" maxlog = 1
		end
	end

	return nothing
end

# ============================================================================
# Main Solve Function
# ============================================================================

function solve(problem::PyrolysisProblem;
	solver = Rodas5P(autodiff = AutoFiniteDiff()),  # FiniteDiff for Jacobian (ForwardDiff needs code changes)
	abstol::Real = 1e-8,
	reltol::Real = 1e-6,
	maxiters::Int = 1000000,
	dtmin::Real = 1e-12,
	dtmax::Real = Inf,
	saveat = problem.output_times,
	callback = nothing,
	progress::Bool = true,
	validate::Bool = true,
	use_ale::Bool = false,
	min_thickness::Real = 1e-6,
	# Surface cell depletion handling
	handle_depletion::Bool = false,
	depletion_threshold::Real = 0.05,
	min_cells::Int = 2,
	# Sparse Jacobian support (significant speedup for large problems)
	use_sparse_jacobian::Bool = false,
	# Analytical Jacobian (experimental, significant speedup)
	use_analytical_jacobian::Bool = false,
	# Radiation model selection
	radiation_model::RadiationModel = NO_RADIATION,
	radiation_bcs::Union{Nothing, RadiationBCSet} = nothing)

	# Validate problem
	if validate
		validate_problem(problem; strict = false)
	end

	# Get problem parameters (needed for validation)
	mesh = problem.mesh
	material = problem.material
	tspan = problem.tspan
	NC = n_components(material)

	# Validate radiation model selection
	_validate_radiation_model(radiation_model, material, radiation_bcs)

	# Apply initial conditions to mesh
	set_initial_conditions!(problem)

	# Create solution container
	solution = PyrolysisSolution(problem)

	# Start timing
	start_time = time()

	# Create ODE function and state vector based on ALE mode and radiation model
	local ode_fn, u0, sv_ale

	# Get radiation flux function from boundary conditions
	radiation_flux = _get_radiation_flux_function(problem)

	# Determine I_radiation based on radiation model
	# For P1, radiation is handled via Marshak BCs (I_radiation = 0)
	# For SURFACE_ABSORPTION and BEER_LAMBERT, use external radiation flux
	# For NO_RADIATION, no external flux (I_radiation = 0)
	I_radiation = if radiation_model == P1_QUASI_STEADY || radiation_model == NO_RADIATION
		0.0
	else
		radiation_flux(0.0)  # Initial radiation intensity
	end

	if use_ale
		# Create ALE ODE function with radiation model
		ode_fn = PyrolysisODEFunctionALE(
			mesh, material, problem.bc_set;
			I_radiation = I_radiation,
			radiation_model = radiation_model,
			min_thickness = Float64(min_thickness),
		)
		u0 = create_state_vector_ale(mesh)
		sv_ale = StateVectorALE(mesh)
	else
		# Fixed mesh mode
		if radiation_model == P1_QUASI_STEADY
			error("P1_QUASI_STEADY radiation currently requires use_ale=true")
		end

		# Create fixed-mesh ODE function with radiation model
		ode_fn = PyrolysisODEFunction(
			mesh, material, problem.bc_set;
			I_radiation = I_radiation,
			radiation_model = radiation_model,
		)
		u0 = create_state_vector(mesh, NC)
		sv_ale = nothing
	end

	# Wrap with P1 radiation if requested (requires ALE mode)
	if radiation_model == P1_QUASI_STEADY
		if !use_ale
			error("P1_QUASI_STEADY radiation currently requires use_ale=true")
		end
		# Create radiation cache and boundary conditions
		rad_cache = Properties.RadiationCache(mesh.n_active_cells)

		# Get or create radiation boundary conditions
		if radiation_bcs === nothing
			# Create default Marshak BCs from thermal BCs
			radiation_flux = _get_radiation_flux_function(problem)
			ε_top, T_amb_top = _extract_radiative_params(problem.bc_set.thermal_top)
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(
					I_external = radiation_flux,
					T_ambient = T_amb_top,
					emissivity = ε_top,
				),
				bottom = ReflectiveRadiationBC(),
			)
		end

		# Wrap ODE function with P1 radiation
		ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn, rad_cache, radiation_bcs)
		@info "Using P1_QUASI_STEADY volumetric radiation model"
	elseif radiation_model == SURFACE_ABSORPTION
		@info "Using SURFACE_ABSORPTION radiation model"
	elseif radiation_model == BEER_LAMBERT
		@info "Using BEER_LAMBERT radiation model"
	end

	# Build ODE problem with optional Jacobian specification
	local ode_problem
	local jac_cache = nothing  # For analytical Jacobian

	if use_analytical_jacobian && use_sparse_jacobian
		error("Cannot use both use_analytical_jacobian and use_sparse_jacobian simultaneously")
	end

	if use_analytical_jacobian
		n_cells_val = length(mesh.cells)
		n_faces_val = length(mesh.faces)
		NR = n_reactions(material)

		if use_ale
			# ALE analytical Jacobian with sparse storage
			n_nodes_val = length(mesh.nodes)
			jac_cache = JacobianCacheALE(n_cells_val, n_faces_val, n_nodes_val, NC, NR)

			# Create sparsity pattern for efficient LU factorization
			jac_prototype = jacobian_sparsity_ale(n_cells_val, NC)

			# Create wrapper function for ALE analytical Jacobian
			# Note: J will be a sparse matrix matching jac_prototype structure
			function analytical_jac_ale!(J, u, p, t)
				compute_analytical_jacobian_ale!(J, u, p, t, ode_fn, jac_cache)
			end

			# Create ODEFunction with analytical Jacobian AND sparsity pattern
			# The jac_prototype enables sparse LU factorization (much faster for large systems)
			ode_fn_with_jac = OrdinaryDiffEq.ODEFunction(ode_fn;
				jac = analytical_jac_ale!,
				jac_prototype = jac_prototype)
			ode_problem = ODEProblem(ode_fn_with_jac, u0, tspan)

			@info "Using ALE analytical Jacobian with sparse storage ($(nnz(jac_prototype)) non-zeros)"
		else
			# Fixed-mesh analytical Jacobian with sparse storage
			jac_cache = JacobianCache(n_cells_val, n_faces_val, NC, NR)

			# Create sparsity pattern for efficient LU factorization
			jac_prototype = jacobian_sparsity(n_cells_val, NC)

			# Create wrapper function for analytical Jacobian
			function analytical_jac!(J, u, p, t)
				compute_analytical_jacobian!(J, u, p, t, ode_fn, jac_cache)
			end

			# Create ODEFunction with analytical Jacobian AND sparsity pattern
			ode_fn_with_jac = OrdinaryDiffEq.ODEFunction(ode_fn;
				jac = analytical_jac!,
				jac_prototype = jac_prototype)
			ode_problem = ODEProblem(ode_fn_with_jac, u0, tspan)

			@info "Using analytical Jacobian with sparse storage ($(nnz(jac_prototype)) non-zeros)"
		end

		# Use solver without autodiff since we provide the Jacobian
		# Also use KLU for sparse LU factorization - 10× less allocations than UMFPACK
		if solver isa Rodas5P && solver.autodiff isa AutoFiniteDiff
			solver = Rodas5P(linsolve = KLUFactorization())  # No autodiff, KLU for sparse
		elseif solver isa Rodas5P && (solver.linsolve === nothing || !(solver.linsolve isa KLUFactorization))
			# User provided Rodas5P but without KLU - switch to KLU for performance
			solver = Rodas5P(linsolve = KLUFactorization())
		end
	elseif use_sparse_jacobian
		# Compute sparsity pattern and coloring for efficient finite differences
		n_cells_val = length(mesh.cells)
		if use_ale
			jac_prototype = jacobian_sparsity_ale(n_cells_val, NC)
		else
			jac_prototype = jacobian_sparsity(n_cells_val, NC)
		end
		colorvec = matrix_colors(jac_prototype)

		# Create ODEFunction with sparsity information
		# The jac_prototype tells the solver the structure, colorvec reduces FD evaluations
		# The default AutoFiniteDiff will automatically use colorvec for compressed FD
		ode_fn_sparse = OrdinaryDiffEq.ODEFunction(ode_fn;
			jac_prototype = jac_prototype,
			colorvec = colorvec)
		ode_problem = ODEProblem(ode_fn_sparse, u0, tspan)
	else
		ode_problem = ODEProblem(ode_fn, u0, tspan)
	end

	# ============================================================================
	# MKL LU Factorization Workaround for Depletion Handling
	# ============================================================================
	# When using handle_depletion=true, the state vector is dynamically resized
	# during integration via resize!(integrator, new_length). For large problems
	# (state vector length > ~500), LinearSolve.jl's DefaultLinearSolver switches
	# to MKLLUFactorization which has a bug where the pivot array (ipiv) is not
	# properly resized, causing "DimensionMismatch: ipiv has length X, but needs
	# to be Y" errors.
	#
	# Workaround: When depletion handling is enabled and the state vector is large,
	# we switch to Julia's native LUFactorization which handles resize! correctly.
	# This has a small performance cost but ensures correctness.
	#
	# Reference: https://github.com/SciML/LinearSolve.jl (potential upstream bug)
	# ============================================================================
	MKL_RESIZE_THRESHOLD = 500  # State vector lengths above this trigger MKL
	if use_ale && handle_depletion && length(u0) > MKL_RESIZE_THRESHOLD
		# Check if user provided a custom solver or we're using the default
		if solver isa Rodas5P && !hasproperty(solver, :linsolve) ||
		   (hasproperty(solver, :linsolve) && solver.linsolve === nothing)
			@warn """Large state vector ($(length(u0)) elements) with depletion handling enabled.
			Switching from MKL to Julia's native LUFactorization to avoid resize! bug.
			This may be slightly slower but ensures correctness during cell merging.
			To suppress this warning, explicitly pass: solver=Rodas5P(linsolve=LUFactorization())"""
			solver = Rodas5P(autodiff = AutoFiniteDiff(), linsolve = LUFactorization())
		end
	end

	# Build solver options
	solve_kwargs = Dict{Symbol, Any}(
		:abstol => abstol,
		:reltol => reltol,
		:maxiters => maxiters,
		:dtmin => dtmin,
		:dtmax => dtmax,
		:dt => problem.dt_initial,
		:progress => false,  # Disable OrdinaryDiffEq's progress, we use our own
	)

	# When handling depletion with resize!, we cannot use saveat with interpolation
	# because the interpolation doesn't handle state vector size changes.
	# Instead, we save every step and filter to desired times in post-processing.
	if use_ale && handle_depletion
		# Save every step (dense output without interpolation)
		solve_kwargs[:save_everystep] = true
		# Store the desired saveat times for post-processing
		requested_saveat = saveat
	else
		solve_kwargs[:saveat] = saveat
		requested_saveat = nothing
	end

	# Build callback list
	callbacks_list = []

	# Create shared callback state if using depletion handling (so both callbacks see updates)
	depletion_cb_state = nothing
	if use_ale && handle_depletion
		# Pass mesh so DepletionCallbackState can compute initial cell thickness
		depletion_cb_state = DepletionCallbackState(sv_ale, mesh)
	end

	# Add thickness termination callback for ALE mode
	# Use shared state if available so it sees updated sv after cell merges
	if use_ale
		sv_ref = depletion_cb_state !== nothing ? depletion_cb_state : sv_ale
		thickness_cb = create_thickness_termination_callback(sv_ref, Float64(min_thickness))
		push!(callbacks_list, thickness_cb)
	end

	# Add surface cell depletion callback for ALE mode
	if use_ale && handle_depletion
		depletion_cb = create_depletion_callback(ode_fn, mesh, material, depletion_cb_state;
			threshold = Float64(depletion_threshold),
			min_cells = min_cells)
		push!(callbacks_list, depletion_cb)
	end

	# Add progress callback if requested
	prog = nothing
	if progress
		t_start, t_end = tspan
		t_span_length = t_end - t_start
		last_progress = Ref(0)

		# Create progress bar
		prog = Progress(100; desc = "Solving: ", showspeed = true, enabled = true)

		# Use PeriodicCallback to update progress bar at regular time intervals
		check_interval = t_span_length / 100  # Check every 1% of simulation time

		progress_cb = PeriodicCallback(
			integrator -> begin
				t = integrator.t
				current_progress = floor(Int, 100 * (t - t_start) / t_span_length)
				if current_progress > last_progress[]
					update!(prog, current_progress)
					last_progress[] = current_progress
				end
			end,
			check_interval;
			initial_affect = false,
			save_positions = (false, false),  # Don't save solution at progress check times
		)
		push!(callbacks_list, progress_cb)
	end

	# Add user callback if provided
	if callback !== nothing
		push!(callbacks_list, callback)
	end

	# Combine callbacks
	if !isempty(callbacks_list)
		if length(callbacks_list) == 1
			solve_kwargs[:callback] = callbacks_list[1]
		else
			solve_kwargs[:callback] = CallbackSet(callbacks_list...)
		end
	end

	# Use sparse solver if configured, otherwise user-provided or default

	# Solve
	raw_solution = OrdinaryDiffEq.solve(ode_problem, solver; solve_kwargs...)

	# Finish progress bar if it was created
	if prog !== nothing
		finish!(prog)
	end

	# Record solve time
	solution.solve_time = time() - start_time

	# Extract and store results based on ALE mode
	if use_ale
		extract_solution_ale!(solution, raw_solution, problem, sv_ale;
			requested_times = requested_saveat,
			radiation_model = radiation_model)
	else
		extract_solution!(solution, raw_solution, problem;
			radiation_model = radiation_model)
	end

	return solution
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	_get_radiation_flux_function(problem::PyrolysisProblem) -> Function

Get the radiation flux function, preferring RadiativeFluxBC from thermal_bcs over external_flux.

Returns a function t -> I_absorbed that gives the absorbed radiation intensity [W/m²].

Priority:
1. RadiativeFluxBC on top boundary (boundary_id=1) with absorptivity factor
2. external_flux from problem (alternative interface for simple cases)

Note: The radiation_model passed to solve() determines how the absorbed flux is used:
- SURFACE_ABSORPTION: flux added to surface cell via thermal_bc_flux
- BEER_LAMBERT: flux used for volumetric Beer-Lambert absorption (I_radiation)
- P1_QUASI_STEADY: handled via Marshak BC
"""
function _get_radiation_flux_function(problem::PyrolysisProblem)
	# Check for RadiativeFluxBC on top surface (boundary 1)
	top_boundary_id = 1
	if has_radiative_flux_bc(problem.bc_set, top_boundary_id)
		# Extract RadiativeFluxBC using type-stable accessor
		bc = get_thermal_bc(problem.bc_set, Val(top_boundary_id))
		rad_bc = get_radiative_flux_bc(bc)

		# Return a function that computes absorbed flux
		# This is used by BEER_LAMBERT and P1 models for volumetric absorption
		return function (t)
			I, _ = extract_radiative_flux_from_bcs(problem.bc_set, top_boundary_id, t)
			return I
		end
	else
		# Fallback to external_flux
		return problem.external_flux
	end
end

"""
	_create_radiation_bcs_from_thermal(bc_set::BoundaryConditionSet, radiation_flux::Function) -> RadiationBCSet

Create radiation boundary conditions (for P1 solver) from thermal boundary conditions.

Converts thermal BCs to MarshakRadiationBC on top surface (boundary 1) and
ReflectiveRadiationBC on bottom surface (boundary 2), returning a type-stable
RadiationBCSet (no Dict allocations in hot paths).

# Arguments
- `bc_set`: BoundaryConditionSet with thermal boundary conditions
- `radiation_flux`: Function t -> I_external for external radiation source

# Returns
`RadiationBCSet` (top, bottom)

# Notes
- Extracts emissivity from RadiativeBC if present, otherwise uses default 0.9
- Extracts T_ambient from RadiativeBC if present, otherwise uses 300.0 K
- Top surface (id=1) gets MarshakRadiationBC with external radiation
- Bottom surface (id=2) gets ReflectiveRadiationBC (adiabatic for radiation)
"""
function _create_radiation_bcs_from_thermal(bc_set::BoundaryConditionSet, radiation_flux::Function)
	# Process top surface (boundary 1 - exposed)
	top_bc = get_thermal_bc(bc_set, Val(1))
	emissivity, T_ambient = _extract_radiative_params(top_bc)

	return RadiationBCSet(
		top = MarshakRadiationBC(
			I_external = radiation_flux,
			T_ambient = T_ambient,
			emissivity = emissivity,
		),
		bottom = ReflectiveRadiationBC(),
	)
end

"""
	_extract_radiative_params(bc::AbstractThermalBC) -> Tuple{Any, Any}

Extract emissivity and T_ambient from a (possibly combined) thermal BC.
Returns (emissivity, T_ambient) with defaults (0.9, 300.0) if not found.
"""
_extract_radiative_params(bc::AbstractThermalBC) = (0.9, 300.0)

function _extract_radiative_params(bc::RadiativeBC)
	return (bc.ε, bc.T_inf)
end

function _extract_radiative_params(bc::CombinedThermalBC)
	for inner_bc in bc.bcs
		if inner_bc isa RadiativeBC
			return (inner_bc.ε, inner_bc.T_inf)
		end
	end
	# Recursively check for nested combined BCs
	for inner_bc in bc.bcs
		if inner_bc isa CombinedThermalBC
			result = _extract_radiative_params(inner_bc)
			if result[1] != 0.9 || result[2] != 300.0
				return result
			end
		end
	end
	return (0.9, 300.0)
end

"""
	create_ode_function(problem::PyrolysisProblem)

Create the appropriate ODE function for a problem (fixed mesh mode).
"""
function create_ode_function(problem::PyrolysisProblem)
	mesh = problem.mesh
	material = problem.material
	bc_set = problem.bc_set

	# Get radiation flux function (prefers RadiativeFluxBC over external_flux)
	radiation_flux = _get_radiation_flux_function(problem)

	if radiation_flux(0.0) == 0.0 && radiation_flux(1.0) == 0.0
		# No external radiation, use basic ODE function
		return PyrolysisODEFunction(mesh, material, bc_set;
			I_radiation = 0.0)
	else
		# Time-varying radiation, use wrapper
		base_fn = PyrolysisODEFunction(mesh, material, bc_set;
			I_radiation = 0.0)
		return PyrolysisODEFunctionWithRadiation(base_fn, radiation_flux)
	end
end

"""
	create_ode_function_ale(problem::PyrolysisProblem; kwargs...)

Create the ALE ODE function for a problem (moving mesh mode).

This function integrates mesh motion directly into the ODE, ensuring
node positions evolve simultaneously with temperature and concentration.
"""
function create_ode_function_ale(problem::PyrolysisProblem;
	min_thickness::Float64 = 1e-6)

	mesh = problem.mesh
	material = problem.material

	# Use bc_set directly from problem (already type-stable)
	bc_set = problem.bc_set

	# Get radiation flux function (prefers RadiativeFluxBC over external_flux)
	radiation_flux = _get_radiation_flux_function(problem)

	# Check if external flux is time-varying
	has_varying_radiation = !(radiation_flux(0.0) == 0.0 && radiation_flux(1.0) == 0.0)

	# Create ALE ODE function with type-stable BoundaryConditionSet
	base_fn = PyrolysisODEFunctionALE(
		mesh, material, bc_set;
		I_radiation = radiation_flux(0.0),
		min_thickness = min_thickness,
	)

	if has_varying_radiation
		# Wrap with time-varying radiation
		return PyrolysisODEFunctionALEWithRadiation(base_fn, radiation_flux)
	else
		return base_fn
	end
end

"""
	extract_solution!(solution::PyrolysisSolution, raw_solution, problem; radiation_model)

Extract results from raw ODE solution into PyrolysisSolution structure.
"""
function extract_solution!(solution::PyrolysisSolution, raw_solution, problem;
	radiation_model::RadiationModel = SURFACE_ABSORPTION)
	mesh = problem.mesh
	material = problem.material
	NC = n_components(material)
	NR = n_reactions(material)
	n_cells_val = length(mesh.cells)

	# Store raw solution for interpolation
	solution.raw_solution = raw_solution

	# Get return code
	solution.retcode = Symbol(raw_solution.retcode)

	# Extract times
	times = raw_solution.t
	solution.times = copy(times)

	# Pre-allocate output arrays
	n_times = length(times)
	solution.temperatures = Vector{Vector{Float64}}(undef, n_times)
	solution.concentrations = Vector{Matrix{Float64}}(undef, n_times)
	solution.gas_generation_rate = zeros(n_times)
	solution.surface_temperature = zeros(n_times)
	solution.back_temperature = zeros(n_times)
	solution.total_mass = zeros(n_times)

	# Pre-allocate boundary flux arrays
	solution.surface_mass_flux = zeros(n_times)
	solution.surface_heat_flux_net = zeros(n_times)
	solution.surface_heat_flux_absorbed = zeros(n_times)
	solution.surface_heat_flux_reradiation = zeros(n_times)
	solution.surface_heat_flux_convection = zeros(n_times)
	solution.surface_heat_flux_conduction = zeros(n_times)
	solution.surface_emissivity = zeros(n_times)
	solution.surface_conductivity = zeros(n_times)
	solution.back_heat_flux = zeros(n_times)

	# Create state vector descriptor
	sv = StateVector(n_cells_val, NC)

	# Pre-allocate work arrays for reaction rate computation
	rates = zeros(NR)
	heats = zeros(NR)

	# Get boundary conditions for flux computation
	bc_set = problem.bc_set
	thermal_bc_top = get_thermal_bc(bc_set, Val(1))
	thermal_bc_bottom = get_thermal_bc(bc_set, Val(2))
	mass_bc_top = get_mass_bc(bc_set, Val(1))

	# Extract data at each saved time
	for (i, t) in enumerate(times)
		u = raw_solution.u[i]

		# Extract temperatures
		T_profile = zeros(n_cells_val)
		for j in 1:n_cells_val
			T_profile[j] = u[j]
		end
		solution.temperatures[i] = T_profile

		# Extract concentrations
		ξ_profile = zeros(NC, n_cells_val)
		for j in 1:n_cells_val
			for k in 1:NC
				ξ_profile[k, j] = u[n_cells_val+(k-1)*n_cells_val+j]
			end
		end
		solution.concentrations[i] = ξ_profile

		# ================================================================
		# Solve for true boundary temperatures (not cell center values)
		# ================================================================
		# Surface (top, boundary_id=1)
		T_cell_top = T_profile[n_cells_val]
		ξ_top = ntuple(k -> ξ_profile[k, n_cells_val], NC)
		k_top = effective_conductivity(ξ_top, T_cell_top, material)
		ε_top = effective_emissivity(ξ_top, T_cell_top, material)
		Δz_half_top = mesh.cell_volumes[n_cells_val] / mesh.face_areas[mesh.cells[n_cells_val].face_ids[2]] / 2
		state_top = solve_boundary_temperature(T_cell_top, k_top, Δz_half_top, thermal_bc_top, t, ε_top;
			radiation_model = radiation_model)
		solution.surface_temperature[i] = state_top.T_surface

		# Back (bottom, boundary_id=2)
		T_cell_bottom = T_profile[1]
		ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
		k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
		ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
		Δz_half_bottom = mesh.cell_volumes[1] / mesh.face_areas[mesh.cells[1].face_ids[1]] / 2
		state_bottom = solve_boundary_temperature(T_cell_bottom, k_bottom, Δz_half_bottom, thermal_bc_bottom, t, ε_bottom;
			radiation_model = radiation_model)
		solution.back_temperature[i] = state_bottom.T_surface

		# Total mass (integrate ξ over domain)
		total_m = 0.0
		for j in 1:n_cells_val
			cell_vol = mesh.cell_volumes[j]
			for k in 1:NC
				total_m += ξ_profile[k, j] * cell_vol
			end
		end
		# Normalize by cross-sectional area for kg/m²
		solution.total_mass[i] = total_m / mesh.cross_section_area

		# ================================================================
		# Boundary flux computation (using solved surface temperature)
		# ================================================================

		# Surface cell properties (top boundary, z=L, boundary_id=1)
		# Note: ξ_top, k_top, ε_top already computed above for boundary temperature solve
		T_surface_solved = state_top.T_surface
		ξ_surface = ξ_top

		# Store effective properties
		solution.surface_emissivity[i] = ε_top
		solution.surface_conductivity[i] = k_top

		# Compute heat flux components using SOLVED surface temperature
		heat_fluxes = compute_boundary_heat_fluxes(thermal_bc_top, T_surface_solved, t, ε_top, radiation_model)
		solution.surface_heat_flux_absorbed[i] = heat_fluxes.absorbed
		solution.surface_heat_flux_reradiation[i] = heat_fluxes.reradiation
		solution.surface_heat_flux_convection[i] = heat_fluxes.convection
		solution.surface_heat_flux_net[i] = heat_fluxes.net

		# Conductive heat flux at surface: use solved boundary state
		# This is consistent with the boundary temperature solve and ensures energy balance
		solution.surface_heat_flux_conduction[i] = state_top.q_conducted

		# Compute total mass flux leaving at surface (actual MLR)
		total_mass_flux = 0.0
		d_surface = mesh.cell_volumes[n_cells_val] / mesh.cross_section_area / 2  # Half cell thickness

		for k in 1:NC
			comp = material.components[k]
			if comp.phase == GAS
				ξ_k = ξ_surface[k]
				if ξ_k > 0.0
					λρ = effective_gas_transfer(ξ_surface, T_surface_solved, material) * sum(ξ_surface)
					ρ_surface = sum(ξ_surface)
					J_k = compute_surface_mass_flux(mass_bc_top, ξ_k, λρ, ρ_surface, d_surface, t)
					total_mass_flux += J_k
				end
			end
		end
		solution.surface_mass_flux[i] = total_mass_flux

		# Back face heat flux (boundary_id=2, z=0) - use solved boundary temperature
		T_back_solved = state_bottom.T_surface
		ξ_back = ξ_bottom
		back_fluxes = compute_boundary_heat_fluxes(thermal_bc_bottom, T_back_solved, t, ε_bottom, radiation_model)
		solution.back_heat_flux[i] = back_fluxes.net

		# Gas generation rate: compute from reaction rates (physics-based)
		# This is NOT the actual mass loss rate (surface flux), but rather the
		# volumetric production of gases throughout the sample.
		# Total gas generation = Σ_cells Σ_reactions (gas_produced * r * V)
		# This is much smoother than numerical differentiation of total mass
		gen_rate = 0.0
		for j in 1:n_cells_val
			T = T_profile[j]
			ξ = ntuple(k -> ξ_profile[k, j], NC)
			cell_vol = mesh.cell_volumes[j]

			# Evaluate reaction rates at this cell
			Physics.evaluate_reactions!(rates, heats, ξ, T, material)

			# Sum gas production from each reaction
			for r_idx in 1:NR
				rxn = material.reactions[r_idx]
				r = rates[r_idx]

				# Gas production: positive stoichiometry for gas products
				if rxn.product1 > 0 && material.components[rxn.product1].phase == GAS
					gen_rate += rxn.stoich[3] * r * cell_vol
				end
				if rxn.product2 > 0 && material.components[rxn.product2].phase == GAS
					gen_rate += rxn.stoich[4] * r * cell_vol
				end
			end
		end
		# Normalize by cross-sectional area for kg/(m²·s)
		solution.gas_generation_rate[i] = gen_rate / mesh.cross_section_area
	end

	# For non-ALE mode, thickness is constant (initial mesh thickness)
	initial_thickness = Geometry.domain_length(mesh)
	solution.thickness_history = fill(initial_thickness, n_times)

	return solution
end

"""
	extract_solution_ale!(solution::PyrolysisSolution, raw_solution, problem, sv::StateVectorALE)

Extract results from ALE ODE solution into PyrolysisSolution structure.

For ALE solutions, node positions are part of the state vector, so we can
reconstruct the geometry at each saved time and compute volumes accordingly.
This ensures proper mass accounting with the evolving mesh.

Note: When surface cell depletion handling is enabled, the state vector length
may vary between saved times. This function detects the current number of cells
from the state vector length and adjusts accordingly.

# Keyword Arguments
- `requested_times::Union{Nothing,Vector{Float64}}`: If provided, filter the raw solution
  to only include times closest to these requested times. This is used when
  depletion handling saves every step and we want specific output times.
- `radiation_model::RadiationModel`: Radiation model for boundary flux diagnostics.
"""
function extract_solution_ale!(solution::PyrolysisSolution, raw_solution,
	problem, sv_initial::StateVectorALE{NC};
	requested_times::Union{Nothing, AbstractVector{<:Real}} = nothing,
	radiation_model::RadiationModel = SURFACE_ABSORPTION) where NC
	mesh = problem.mesh
	material = problem.material
	NR = n_reactions(material)

	# Store raw solution for interpolation
	solution.raw_solution = raw_solution

	# Get return code
	solution.retcode = Symbol(raw_solution.retcode)

	# Determine which time indices to extract
	all_times = raw_solution.t
	if requested_times !== nothing && !isempty(requested_times)
		# Find closest actual time for each requested time
		indices_to_extract = Int[]
		for t_req in requested_times
			if t_req > all_times[end]
				break  # Skip requested times beyond simulation end
			end
			idx = argmin(abs.(all_times .- t_req))
			if isempty(indices_to_extract) || idx != indices_to_extract[end]
				push!(indices_to_extract, idx)
			end
		end
		# Always include the final time
		if isempty(indices_to_extract) || indices_to_extract[end] != length(all_times)
			push!(indices_to_extract, length(all_times))
		end
	else
		indices_to_extract = collect(1:length(all_times))
	end

	# Extract times at selected indices
	times = all_times[indices_to_extract]
	solution.times = copy(times)
	n_times = length(times)

	# Pre-allocate output arrays
	solution.temperatures = Vector{Vector{Float64}}(undef, n_times)
	solution.concentrations = Vector{Matrix{Float64}}(undef, n_times)
	solution.gas_generation_rate = zeros(n_times)
	solution.surface_temperature = zeros(n_times)
	solution.back_temperature = zeros(n_times)
	solution.total_mass = zeros(n_times)
	solution.thickness_history = zeros(n_times)

	# Pre-allocate boundary flux arrays
	solution.surface_mass_flux = zeros(n_times)
	solution.surface_heat_flux_net = zeros(n_times)
	solution.surface_heat_flux_absorbed = zeros(n_times)
	solution.surface_heat_flux_reradiation = zeros(n_times)
	solution.surface_heat_flux_convection = zeros(n_times)
	solution.surface_heat_flux_conduction = zeros(n_times)
	solution.surface_emissivity = zeros(n_times)
	solution.surface_conductivity = zeros(n_times)
	solution.back_heat_flux = zeros(n_times)

	# Pre-allocate work arrays for reaction rate computation
	rates = zeros(NR)
	heats = zeros(NR)

	# Get boundary conditions for flux computation
	bc_set = problem.bc_set
	thermal_bc_top = get_thermal_bc(bc_set, Val(1))
	thermal_bc_bottom = get_thermal_bc(bc_set, Val(2))
	mass_bc_top = get_mass_bc(bc_set, Val(1))

	# Extract data at each saved time
	for (out_idx, raw_idx) in enumerate(indices_to_extract)
		t = all_times[raw_idx]
		u = raw_solution.u[raw_idx]

		# Determine current number of cells from state vector length
		# State vector: [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
		# Length = n + NC*n + (n+1) = n*(1+NC) + n + 1
		# Solving: len = n*(1 + NC + 1) + 1 = n*(NC + 2) + 1
		# So: n = (len - 1) / (NC + 2)
		state_len = length(u)
		n_cells = (state_len - 1) ÷ (NC + 2)
		n_nodes = n_cells + 1

		# Create state vector descriptor for current size
		sv = StateVectorALE(n_cells, Val(NC))

		# Extract temperatures
		T_profile = [get_temperature_ale(u, sv, j) for j in 1:n_cells]
		solution.temperatures[out_idx] = T_profile

		# Extract concentrations
		ξ_profile = zeros(NC, n_cells)
		for j in 1:n_cells
			for k in 1:NC
				ξ_profile[k, j] = get_concentration_ale(u, sv, j, k)
			end
		end
		solution.concentrations[out_idx] = ξ_profile

		# Extract node positions and compute geometry
		z = get_node_positions(u, sv)
		thickness = z[end] - z[1]
		solution.thickness_history[out_idx] = thickness

		# Compute cell volumes from current node positions
		volumes = compute_volumes_from_state(u, sv, mesh.cross_section_area)

		# ================================================================
		# Solve for true boundary temperatures (not cell center values)
		# ================================================================
		# Surface (top, boundary_id=1)
		T_cell_top = T_profile[n_cells]
		ξ_top = ntuple(k -> ξ_profile[k, n_cells], NC)
		k_top = effective_conductivity(ξ_top, T_cell_top, material)
		ε_top = effective_emissivity(ξ_top, T_cell_top, material)
		Δz_half_top = (z[n_cells+1] - z[n_cells]) / 2  # Half of surface cell
		state_top = solve_boundary_temperature(T_cell_top, k_top, Δz_half_top, thermal_bc_top, t, ε_top;
			radiation_model = radiation_model)
		solution.surface_temperature[out_idx] = state_top.T_surface

		# Back (bottom, boundary_id=2)
		T_cell_bottom = T_profile[1]
		ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
		k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
		ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
		Δz_half_bottom = (z[2] - z[1]) / 2  # Half of bottom cell
		state_bottom = solve_boundary_temperature(T_cell_bottom, k_bottom, Δz_half_bottom, thermal_bc_bottom, t, ε_bottom;
			radiation_model = radiation_model)
		solution.back_temperature[out_idx] = state_bottom.T_surface

		# Total mass (integrate ξ × V using current geometry)
		total_m = 0.0
		for j in 1:n_cells
			for k in 1:NC
				total_m += ξ_profile[k, j] * volumes[j]
			end
		end
		# Normalize by cross-sectional area for kg/m²
		solution.total_mass[out_idx] = total_m / mesh.cross_section_area

		# ================================================================
		# Boundary flux computation (using solved surface temperature)
		# ================================================================

		# Surface cell properties - reuse variables from boundary temperature solve
		T_surface_solved = state_top.T_surface
		ξ_surface = ξ_top

		# Store effective properties
		solution.surface_emissivity[out_idx] = ε_top
		solution.surface_conductivity[out_idx] = k_top

		# Compute heat flux components using SOLVED surface temperature
		heat_fluxes = compute_boundary_heat_fluxes(thermal_bc_top, T_surface_solved, t, ε_top, radiation_model)
		solution.surface_heat_flux_absorbed[out_idx] = heat_fluxes.absorbed
		solution.surface_heat_flux_reradiation[out_idx] = heat_fluxes.reradiation
		solution.surface_heat_flux_convection[out_idx] = heat_fluxes.convection
		solution.surface_heat_flux_net[out_idx] = heat_fluxes.net

		# Conductive heat flux at surface: use solved boundary state
		# This is consistent with the boundary temperature solve and ensures energy balance
		solution.surface_heat_flux_conduction[out_idx] = state_top.q_conducted

		# Compute total mass flux leaving at surface (actual MLR)
		# Sum over all gas components
		total_mass_flux = 0.0
		d_surface = (z[n_cells+1] - z[n_cells]) / 2  # Distance from cell center to boundary

		for k in 1:NC
			comp = material.components[k]
			if comp.phase == GAS
				ξ_k = ξ_surface[k]
				if ξ_k > 0.0
					# Compute effective transport coefficient at surface
					λρ = effective_gas_transfer(ξ_surface, T_surface_solved, material) * sum(ξ_surface)
					ρ_surface = sum(ξ_surface)

					# Use the mass BC flux computation
					J_k = compute_surface_mass_flux(mass_bc_top, ξ_k, λρ, ρ_surface, d_surface, t)
					total_mass_flux += J_k
				end
			end
		end
		solution.surface_mass_flux[out_idx] = total_mass_flux

		# Back face heat flux (boundary_id=2, z=0) - use solved boundary temperature
		T_back_solved = state_bottom.T_surface
		ξ_back = ξ_bottom
		back_fluxes = compute_boundary_heat_fluxes(thermal_bc_bottom, T_back_solved, t, ε_bottom, radiation_model)
		solution.back_heat_flux[out_idx] = back_fluxes.net

		# Gas generation rate: compute from reaction rates (physics-based)
		# This is NOT the actual mass loss rate (surface flux), but rather the
		# volumetric production of gases throughout the sample.
		gen_rate = 0.0
		for j in 1:n_cells
			T = T_profile[j]
			ξ = ntuple(k -> ξ_profile[k, j], NC)
			cell_vol = volumes[j]

			# Evaluate reaction rates at this cell
			Physics.evaluate_reactions!(rates, heats, ξ, T, material)

			# Sum gas production from each reaction
			for r_idx in 1:NR
				rxn = material.reactions[r_idx]
				r = rates[r_idx]

				# Gas production: positive stoichiometry for gas products
				if rxn.product1 > 0 && material.components[rxn.product1].phase == GAS
					gen_rate += rxn.stoich[3] * r * cell_vol
				end
				if rxn.product2 > 0 && material.components[rxn.product2].phase == GAS
					gen_rate += rxn.stoich[4] * r * cell_vol
				end
			end
		end
		# Normalize by cross-sectional area for kg/(m²·s)
		solution.gas_generation_rate[out_idx] = gen_rate / mesh.cross_section_area
	end

	return solution
end

# ============================================================================
# Alternative Solve Methods
# ============================================================================

"""
	solve_steady_state(problem::PyrolysisProblem; kwargs...) -> PyrolysisSolution

Solve for steady-state solution (if one exists).

This attempts to find the equilibrium state where dT/dt ≈ 0 and dξ/dt ≈ 0.
Note: Many pyrolysis problems don't have a steady state!

# Keyword Arguments
- `solver=Rodas5P()`: ODE solver
- `ss_abstol=1e-8`: Steady-state absolute tolerance
- `ss_reltol=1e-6`: Steady-state relative tolerance
- `tmax=1e6`: Maximum time to simulate
"""
function solve_steady_state(problem::PyrolysisProblem;
	solver = Rodas5P(),
	ss_abstol::Real = 1e-8,
	ss_reltol::Real = 1e-6,
	tmax::Real = 1e6,
	kwargs...)
	# Create modified problem with extended time span
	tspan_ss = (problem.tspan[1], tmax)

	# Use steady-state termination callback
	condition = (u, t, integrator) -> begin
		# Check if derivatives are small
		du = similar(u)
		integrator.f(du, u, nothing, t)
		max_deriv = maximum(abs, du)
		return max_deriv < ss_abstol
	end

	affect! = (integrator) -> terminate!(integrator)
	ss_callback = DiscreteCallback(condition, affect!)

	# Solve with steady-state detection
	return solve(problem;
		solver = solver,
		saveat = [tmax],
		callback = ss_callback,
		kwargs...)
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
	get_profile(solution::PyrolysisSolution, field::Symbol, time_idx::Int)

Get a field profile at a specific saved time index.

# Arguments
- `solution`: The solution
- `field`: :T for temperature, :ξ for all concentrations, or :ξ1, :ξ2, etc.
- `time_idx`: Index into solution.times

# Returns
Vector of values at cell centers.
"""
function get_profile(solution::PyrolysisSolution, field::Symbol, time_idx::Int)
	if field == :T
		return solution.temperatures[time_idx]
	elseif field == :ξ
		return solution.concentrations[time_idx]
	else
		# Parse :ξ1, :ξ2, etc.
		field_str = string(field)
		if startswith(field_str, "ξ") || startswith(field_str, "xi")
			idx_str = replace(field_str, r"^(ξ|xi)" => "")
			comp_idx = parse(Int, idx_str)
			return solution.concentrations[time_idx][comp_idx, :]
		end
	end
	error("Unknown field: $field. Use :T, :ξ, :ξ1, :ξ2, etc.")
end

"""
	get_cell_positions(solution::PyrolysisSolution)

Get cell center positions from the solution's mesh.
"""
function get_cell_positions(solution::PyrolysisSolution)
	return solution.problem.mesh.cell_centers
end

"""
	interpolate_solution(solution::PyrolysisSolution, t::Real)

Interpolate the full state vector at time t.
"""
function interpolate_solution(solution::PyrolysisSolution, t::Real)
	if solution.raw_solution === nothing
		error("Raw solution not available for interpolation")
	end
	return solution.raw_solution(t)
end
