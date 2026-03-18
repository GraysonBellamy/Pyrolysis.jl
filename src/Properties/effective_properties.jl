"""
Effective property computations.

Reference: Design Document Section 2.4

PERFORMANCE NOTE: All functions that iterate over material.components use
dispatch-based unrolling (explicit overloads for small tuple sizes) to avoid
allocations from indexing heterogeneous tuples with runtime indices.
"""

# ============================================================================
# Gas Viscosity (for Darcy flow)
# ============================================================================

"""
	effective_gas_viscosity(T::Float64; μ_ref=1.81e-5, T_ref=300.0, S=110.0) -> Float64

Compute effective gas viscosity using Sutherland's law.

μ(T) = μ_ref × (T/T_ref)^(3/2) × (T_ref + S)/(T + S)

This correlation captures the temperature dependence of gas viscosity
and is valid for air-like pyrolysis gases.

# Arguments
- `T`: Temperature [K]
- `μ_ref`: Reference viscosity at T_ref [Pa·s] (default: 1.81e-5 for air at 300K)
- `T_ref`: Reference temperature [K] (default: 300.0)
- `S`: Sutherland constant [K] (default: 110.0 for air)

# Returns
- `μ`: Dynamic viscosity [Pa·s]

# Example
```julia
μ_300 = effective_gas_viscosity(300.0)   # ≈ 1.81e-5 Pa·s
μ_600 = effective_gas_viscosity(600.0)   # ≈ 3.06e-5 Pa·s
```

# Notes
- Viscosity increases with temperature (unlike liquids)
- Default parameters are for air; adjust S for other gases
- Valid range approximately 200-2000 K
"""
@inline function effective_gas_viscosity(T::Float64;
	μ_ref::Float64 = 1.81e-5,  # Pa·s at 300K (air)
	T_ref::Float64 = 300.0,    # K
	S::Float64 = 110.0)        # Sutherland constant for air
	T_ratio = T / T_ref
	return μ_ref * T_ratio^1.5 * (T_ref + S) / (T + S)
end

"""
	mixture_gas_viscosity(ξ::NTuple{NC,Float64}, T::Float64, material) -> Float64

Compute mixture gas viscosity.

Currently uses a single effective viscosity (Sutherland's law for air-like gases).
Could be extended to use Wilke's mixing rule for multi-component mixtures.

# Arguments
- `ξ`: Mass concentrations [kg/m³] (not currently used)
- `T`: Temperature [K]
- `material`: Material definition (not currently used)

# Returns
- `μ`: Dynamic viscosity [Pa·s]
"""
@inline function mixture_gas_viscosity(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	# For simplicity, use single effective viscosity
	# Could extend to Wilke mixture rule if needed for accuracy with
	# multi-component gas mixtures of very different molecular weights
	return effective_gas_viscosity(T)
end

# ============================================================================
# Chapman-Enskog Diffusivity
# ============================================================================

# Physical constant for Chapman-Enskog formula (Bird, Stewart, Lightfoot)
# D_AB = 1.8583×10⁻³ × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)
# where T[K], M[g/mol], P[atm], σ[Å] → D[cm²/s]
const CHAPMAN_ENSKOG_PREFACTOR = 1.8583e-3  # cm²/s × atm × Å² / (K^1.5 × (g/mol)^0.5)

# Default Lennard-Jones parameters for air (background gas)
const SIGMA_AIR = 3.617   # Å - collision diameter for air
const EPSILON_K_AIR = 97.0  # K - well depth / k_B for air
const M_AIR = 0.029       # kg/mol - molar mass of air

"""
	collision_integral_diffusion(T_star::Float64) -> Float64

Compute the diffusion collision integral Ω_D using the Neufeld correlation.

This is an empirical fit to tabulated collision integrals for the 
Lennard-Jones (12-6) potential.

Ω_D(T*) = A/(T*)^B + C/exp(D·T*) + E/exp(F·T*) + G/exp(H·T*)

Valid for 0.3 ≤ T* ≤ 400 with < 0.1% error.

# Reference
Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972)

# Arguments
- `T_star`: Reduced temperature T* = k_B T / ε_AB

# Returns
- `Ω_D`: Dimensionless diffusion collision integral
"""
@inline function collision_integral_diffusion(T_star::Float64)
	# Neufeld et al. correlation coefficients
	A = 1.06036
	B = 0.15610
	C = 0.19300
	D = 0.47635
	E = 1.03587
	F = 1.52996
	G = 1.76474
	H = 3.89411

	return A / T_star^B + C / exp(D * T_star) + E / exp(F * T_star) + G / exp(H * T_star)
end

"""
	chapman_enskog_diffusivity(T, M; P=P_REF, σ=nothing, ε_k=nothing) -> Float64

Compute binary diffusion coefficient using Chapman-Enskog kinetic theory.

The Chapman-Enskog expression for binary diffusion of species A in background gas B:

D_AB = (1.8583×10⁻⁷) × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)

where:
- T: Temperature [K]
- M_A, M_B: Molar masses [kg/mol]
- P: Pressure [Pa]
- σ_AB: Average collision diameter [Å]
- Ω_D: Diffusion collision integral (function of T*)

# Arguments
- `T`: Temperature [K]
- `M`: Molar mass of diffusing species [kg/mol]
- `P`: Total pressure [Pa] (default: 101325 Pa)
- `σ`: Collision diameter [Å] (default: estimated from Tc or molar mass)
- `ε_k`: Well depth / k_B [K] (default: estimated)

# Returns
- `D`: Binary diffusion coefficient [m²/s]

# Default Parameters
For pyrolysis gases, if σ and ε_k are not provided, uses:
- σ estimated from Tc-based correlation or default to 4 Å
- ε_k estimated as 0.77 × Tc or default to 300 K
- Background gas is air (σ_air = 3.617 Å, ε_k_air = 97 K, M_air = 0.029 kg/mol)

# Examples
```julia
# Water vapor at 500 K
D_H2O = chapman_enskog_diffusivity(500.0, 0.018)  # ≈ 4.0e-5 m²/s

# CO₂ at 600 K  
D_CO2 = chapman_enskog_diffusivity(600.0, 0.044)  # ≈ 3.2e-5 m²/s

# With explicit Lennard-Jones parameters
D_CH4 = chapman_enskog_diffusivity(500.0, 0.016, σ=3.758, ε_k=148.6)
```

# Physical Behavior
- D ∝ T^1.5 to T^1.75 (depending on collision integral)
- D ∝ 1/P (inversely proportional to pressure)
- D ∝ √(1/M) (lighter species diffuse faster)

# Notes
- Valid for dilute gas limit (low density)
- Assumes ideal gas behavior
- Uses combining rules: σ_AB = (σ_A + σ_B)/2, ε_AB = √(ε_A × ε_B)
"""
@inline function chapman_enskog_diffusivity(T::Float64, M::Float64;
	P::Float64 = P_REF,
	σ::Union{Float64, Nothing} = nothing,
	ε_k::Union{Float64, Nothing} = nothing)

	# Estimate Lennard-Jones parameters if not provided
	# Using correlations for organic vapors
	σ_species = if σ !== nothing
		σ
	else
		# Rough estimate: σ ≈ 1.18 × (Vc)^(1/3) where Vc ~ M/ρ_c
		# For organic gases, σ typically 3-6 Å, estimate from M
		3.5 + 0.5 * log10(max(M / 0.016, 1.0))  # Å, anchored to CH4 ≈ 3.8 Å
	end

	ε_k_species = if ε_k !== nothing
		ε_k
	else
		# Rough estimate: ε/k ≈ 0.77 × Tc
		# For organics, ε/k typically 100-600 K
		150.0 + 100.0 * (M / 0.044)  # K, scaled roughly to molecular weight
	end

	# Combining rules for binary interaction parameters
	σ_AB = 0.5 * (σ_species + SIGMA_AIR)  # Å
	ε_k_AB = sqrt(ε_k_species * EPSILON_K_AIR)  # K

	# Reduced temperature
	T_star = T / ε_k_AB

	# Collision integral
	Ω_D = collision_integral_diffusion(T_star)

	# Chapman-Enskog diffusivity
	# Standard formula (Bird, Stewart, Lightfoot):
	# D_AB = 1.8583×10⁻³ × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)
	# where T in K, M in g/mol, P in atm, σ in Å → D in cm²/s
	#
	# Convert units:
	# - M: kg/mol → g/mol (multiply by 1000)
	# - P: Pa → atm (divide by 101325)
	# - D: cm²/s → m²/s (multiply by 1e-4)

	M_gmol = M * 1000.0        # kg/mol → g/mol
	M_air_gmol = M_AIR * 1000.0  # kg/mol → g/mol
	P_atm = P / 101325.0      # Pa → atm

	# Reduced mass factor with M in g/mol
	reduced_mass_factor = sqrt(1.0 / M_gmol + 1.0 / M_air_gmol)

	# D in cm²/s
	D_cm2s = CHAPMAN_ENSKOG_PREFACTOR * T^1.5 * reduced_mass_factor / (P_atm * σ_AB^2 * Ω_D)

	# Convert to m²/s
	D = D_cm2s * 1e-4

	return D
end

"""
	chapman_enskog_diffusivity(T, comp::Component; P=P_REF) -> Float64

Compute Chapman-Enskog diffusivity for a gas component.

Uses component's molar mass and any stored Lennard-Jones parameters.
"""
@inline function chapman_enskog_diffusivity(T::Float64, comp; P::Float64 = P_REF)
	if comp.phase != GAS || comp.molar_mass < 1e-6
		return 0.0  # No diffusion for non-gas components
	end
	return chapman_enskog_diffusivity(T, comp.molar_mass; P = P)
end

"""
	chapman_enskog_λ(T, M; P=P_REF, σ=nothing, ε_k=nothing) -> Float64

Convert Chapman-Enskog diffusivity to equivalent gas transfer coefficient λ.

The relationship between diffusivity D and gas transfer coefficient λ is:
- For Fickian flux: J = -λ ρ ∂Y/∂x = -D ∂ξ/∂x
- Therefore: λ = D / ρ_gas ≈ D × (R T) / (P M)

For consistency with existing λ-based transport, this returns the equivalent λ.

# Arguments
- `T`: Temperature [K]  
- `M`: Molar mass [kg/mol]
- `P`: Pressure [Pa]
- `σ, ε_k`: Optional Lennard-Jones parameters [Å, K]

# Returns  
- `λ`: Equivalent gas transfer coefficient [m²/s]

# Note
This λ should NOT be multiplied by the gas density in flux calculations
when using the Fickian formulation J = -λ ρ ∂Y/∂x, since the density
is already accounted for in the conversion.
"""
@inline function chapman_enskog_λ(T::Float64, M::Float64;
	P::Float64 = P_REF,
	σ::Union{Float64, Nothing} = nothing,
	ε_k::Union{Float64, Nothing} = nothing)

	D = chapman_enskog_diffusivity(T, M; P = P, σ = σ, ε_k = ε_k)

	# Convert D [m²/s] to λ [m²/s]
	# D = λ × (ρ_gas / something) but our formulation is:
	# J = -λ × ρ_face × ∂Y/∂x = -λ × ρ_face × ∂(ξ/ρ)/∂x
	# Standard Fick: J = -D × ∂ξ/∂x
	# So the relationship is complex. For gas in gas diffusion, λ ≈ D is reasonable.
	return D
end

# ============================================================================
# Heat capacity helpers
# ============================================================================

@inline function _sum_heat_capacity(ξ, T, comps::Tuple{C}, ::Val{1}) where C
	return ξ[1] * comps[1].heat_capacity(T)
end

@inline function _sum_heat_capacity(ξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return ξ[1] * comps[1].heat_capacity(T) + ξ[2] * comps[2].heat_capacity(T)
end

@inline function _sum_heat_capacity(ξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return ξ[1] * comps[1].heat_capacity(T) + ξ[2] * comps[2].heat_capacity(T) + ξ[3] * comps[3].heat_capacity(T)
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _sum_heat_capacity(ξ, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		result = 0.0
	end
	for i in 1:NC
		push!(expr.args, :(result += ξ[$i] * comps[$i].heat_capacity(T)))
	end
	push!(expr.args, :(return result))
	return expr
end

"""
	effective_heat_capacity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective volumetric heat capacity [J/(m³·K)].

ρcₚ = Σⱼ ξⱼ cₚ,ⱼ(T)
"""
@inline function effective_heat_capacity(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	return _sum_heat_capacity(ξ, T, material.components, Val(NC))
end

# ============================================================================
# Specific heat (mass-averaged) helpers
# ============================================================================

@inline function _sum_specific_heat(Y, T, comps::Tuple{C}, ::Val{1}) where C
	return Y[1] * comps[1].heat_capacity(T)
end

@inline function _sum_specific_heat(Y, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return Y[1] * comps[1].heat_capacity(T) + Y[2] * comps[2].heat_capacity(T)
end

@inline function _sum_specific_heat(Y, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return Y[1] * comps[1].heat_capacity(T) + Y[2] * comps[2].heat_capacity(T) + Y[3] * comps[3].heat_capacity(T)
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _sum_specific_heat(Y, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		result = 0.0
	end
	for i in 1:NC
		push!(expr.args, :(result += Y[$i] * comps[$i].heat_capacity(T)))
	end
	push!(expr.args, :(return result))
	return expr
end

"""
	effective_specific_heat(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute mass-averaged specific heat capacity [J/(kg·K)].

cₚ = Σⱼ Yⱼ cₚ,ⱼ(T)
"""
@inline function effective_specific_heat(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	Y = mass_fractions(ξ)
	return _sum_specific_heat(Y, T, material.components, Val(NC))
end

# ============================================================================
# Conductivity helpers
# ============================================================================

@inline function _conductivity_sums(v, T, comps::Tuple{C}, ::Val{1}) where C
	k = comps[1].thermal_conductivity(T)
	k_parallel = v[1] * k
	inv_k_series = v[1] > 1e-20 && k > 1e-20 ? v[1] / k : 0.0
	return k_parallel, inv_k_series
end

@inline function _conductivity_sums(v, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	k1 = comps[1].thermal_conductivity(T)
	k2 = comps[2].thermal_conductivity(T)

	k_parallel = v[1] * k1 + v[2] * k2

	inv_k_series = 0.0
	if v[1] > 1e-20 && k1 > 1e-20
		;
		inv_k_series += v[1] / k1;
	end
	if v[2] > 1e-20 && k2 > 1e-20
		;
		inv_k_series += v[2] / k2;
	end

	return k_parallel, inv_k_series
end

@inline function _conductivity_sums(v, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	k1 = comps[1].thermal_conductivity(T)
	k2 = comps[2].thermal_conductivity(T)
	k3 = comps[3].thermal_conductivity(T)

	k_parallel = v[1] * k1 + v[2] * k2 + v[3] * k3

	inv_k_series = 0.0
	if v[1] > 1e-20 && k1 > 1e-20
		;
		inv_k_series += v[1] / k1;
	end
	if v[2] > 1e-20 && k2 > 1e-20
		;
		inv_k_series += v[2] / k2;
	end
	if v[3] > 1e-20 && k3 > 1e-20
		;
		inv_k_series += v[3] / k3;
	end

	return k_parallel, inv_k_series
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _conductivity_sums(v, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		k_parallel = 0.0
		inv_k_series = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			k = comps[$i].thermal_conductivity(T)
			k_parallel += v[$i] * k
			if v[$i] > 1e-20 && k > 1e-20
				inv_k_series += v[$i] / k
			end
		end)
	end
	push!(expr.args, :(return k_parallel, inv_k_series))
	return expr
end

# ============================================================================
# Bruggeman Effective Medium Theory
# ============================================================================

"""
	bruggeman_effective_property(v::NTuple{NC,Float64}, k::NTuple{NC,Float64};
								 tol=1e-10, maxiter=100) -> Float64

Solve for effective property using Bruggeman effective medium theory.

The Bruggeman (symmetric) equation for spherical inclusions is:
	Σⱼ vⱼ (kⱼ - k_eff) / (kⱼ + 2k_eff) = 0

This is solved iteratively using Newton-Raphson.

# Arguments
- `v`: Volume fractions of components (should sum to 1)
- `k`: Property values of components
- `tol`: Convergence tolerance (default: 1e-10)
- `maxiter`: Maximum iterations (default: 100)

# Returns
- `k_eff`: Effective property value

# Physical Interpretation
Bruggeman EMT assumes a self-consistent medium where each phase is surrounded
by the effective medium. It's particularly appropriate for:
- Random heterogeneous mixtures
- Self-similar microstructures (char formation)
- Media near percolation threshold

Unlike parallel/series bounds, Bruggeman can capture:
- Percolation transitions (property drop when continuous phase breaks)
- Non-linear property changes with composition

# Notes
- For two phases, there's an analytical solution (not used here for generality)
- Initial guess is volume-weighted average (parallel mixing)
- Converges quadratically for well-posed problems

# Example
```julia
v = (0.7, 0.3)  # 70% phase 1, 30% phase 2
k = (1.0, 10.0)  # k₁ = 1, k₂ = 10
k_eff = bruggeman_effective_property(v, k)  # ≈ 2.1 (between series=1.43 and parallel=3.7)
```
"""
@inline function bruggeman_effective_property(v::NTuple{NC, Float64}, k::NTuple{NC, Float64};
	tol::Float64 = 1e-10, maxiter::Int = 100) where NC

	# Handle degenerate cases
	total_v = sum(v)
	if total_v < 1e-20
		return 0.0
	end

	# Count non-zero volume components and compute parallel bound for initial guess
	k_eff = 0.0
	n_active = 0
	for i in 1:NC
		if v[i] > 1e-20
			k_eff += v[i] * k[i]  # Parallel mixing as initial guess
			n_active += 1
		end
	end
	k_eff /= total_v  # Normalize

	# Single component - return its property
	if n_active == 1
		return k_eff
	end

	# Newton-Raphson iteration to solve Σ vⱼ(kⱼ - k_eff)/(kⱼ + 2k_eff) = 0
	for _ in 1:maxiter
		f = 0.0      # Function value
		df = 0.0     # Derivative w.r.t. k_eff

		@inbounds for i in 1:NC
			if v[i] > 1e-20 && k[i] > 1e-20
				denom = k[i] + 2.0 * k_eff
				if denom > 1e-20
					f += v[i] * (k[i] - k_eff) / denom
					# d/dk_eff of [v(k - k_eff)/(k + 2k_eff)]
					# = v * [-(k + 2k_eff) - (k - k_eff)*2] / (k + 2k_eff)²
					# = v * [-k - 2k_eff - 2k + 2k_eff] / (k + 2k_eff)²
					# = v * [-3k] / (k + 2k_eff)²
					df += v[i] * (-3.0 * k[i]) / (denom * denom)
				end
			end
		end

		# Newton update with safeguard
		if abs(df) < 1e-20
			break  # Derivative too small, can't continue
		end

		dk = -f / df

		# Limit step size to prevent negative k_eff
		if k_eff + dk < 0.1 * k_eff
			dk = -0.9 * k_eff  # Don't reduce by more than 90%
		end

		k_eff += dk

		# Convergence check
		if abs(f) < tol
			break
		end
	end

	return max(k_eff, 0.0)
end

# Helper to extract conductivities as tuple for Bruggeman
@inline function _conductivity_tuple(T, comps::Tuple{C}, ::Val{1}) where C
	return (comps[1].thermal_conductivity(T),)
end

@inline function _conductivity_tuple(T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return (comps[1].thermal_conductivity(T), comps[2].thermal_conductivity(T))
end

@inline function _conductivity_tuple(T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return (comps[1].thermal_conductivity(T),
		comps[2].thermal_conductivity(T),
		comps[3].thermal_conductivity(T))
end

@generated function _conductivity_tuple(T, comps::Tuple, ::Val{NC}) where NC
	exprs = [:(comps[$i].thermal_conductivity(T)) for i in 1:NC]
	return :(($(exprs...),))
end

"""
	effective_conductivity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective thermal conductivity [W/(m·K)].

Uses the mixing rule specified in the material:
- PARALLEL: k = Σⱼ vⱼ kⱼ
- SERIES: k = 1 / Σⱼ (vⱼ / kⱼ)
- WEIGHTED: k = β·k_parallel + (1-β)·k_series
- BRUGGEMAN: Σⱼ vⱼ(kⱼ - k_eff)/(kⱼ + 2k_eff) = 0

Note: Uses solid_volume_fractions which excludes components with swelling_factor ≈ 0.
This matches ThermaKin behavior: "when γ = 0, gaseous components do not contribute
to the thermal conductivity" (ThermaKin2Ds Reference, Equations 8-10).
"""
@inline function effective_conductivity(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	# Use solid volume fractions - excludes gas with swelling_factor = 0
	v = solid_volume_fractions(ξ, T, material)

	mixing = material.conductivity_mixing
	if mixing == BRUGGEMAN
		k_tuple = _conductivity_tuple(T, material.components, Val(NC))
		return bruggeman_effective_property(v, k_tuple)
	end

	k_parallel, inv_k_series = _conductivity_sums(v, T, material.components, Val(NC))
	k_series = inv_k_series > 1e-20 ? 1.0 / inv_k_series : 0.0

	# Apply mixing rule
	if mixing == PARALLEL
		return k_parallel
	elseif mixing == SERIES
		return k_series
	else  # WEIGHTED
		β = material.β
		return β * k_parallel + (1.0 - β) * k_series
	end
end

# ============================================================================
# Emissivity helpers
# ============================================================================

@inline function _sum_emissivity(v, comps::Tuple{C}, ::Val{1}) where C
	return v[1] * comps[1].emissivity
end

@inline function _sum_emissivity(v, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return v[1] * comps[1].emissivity + v[2] * comps[2].emissivity
end

@inline function _sum_emissivity(v, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return v[1] * comps[1].emissivity + v[2] * comps[2].emissivity + v[3] * comps[3].emissivity
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _sum_emissivity(v, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		result = 0.0
	end
	for i in 1:NC
		push!(expr.args, :(result += v[$i] * comps[$i].emissivity))
	end
	push!(expr.args, :(return result))
	return expr
end

"""
	effective_emissivity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective surface emissivity [-].

ε_eff = Σⱼ vⱼ εⱼ

Uses solid_volume_fractions to exclude gaseous components with swelling_factor ≈ 0.
This matches ThermaKin's treatment of conductivity and is physically correct:
surface emissivity is determined by the solid/liquid material at the surface,
not by gases trapped within the porous structure. Gases with γ = 0 fill existing
pores and don't contribute to surface radiative properties.
"""
@inline function effective_emissivity(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	v = solid_volume_fractions(ξ, T, material)
	ε_eff = _sum_emissivity(v, material.components, Val(NC))
	return clamp(ε_eff, 0.0, 1.0)
end

# ============================================================================
# Absorption helpers
# ============================================================================

@inline function _sum_absorption(ξ, comps::Tuple{C}, ::Val{1}) where C
	return ξ[1] * comps[1].absorption_coeff
end

@inline function _sum_absorption(ξ, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return ξ[1] * comps[1].absorption_coeff + ξ[2] * comps[2].absorption_coeff
end

@inline function _sum_absorption(ξ, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return ξ[1] * comps[1].absorption_coeff + ξ[2] * comps[2].absorption_coeff + ξ[3] * comps[3].absorption_coeff
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _sum_absorption(ξ, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		result = 0.0
	end
	for i in 1:NC
		push!(expr.args, :(result += ξ[$i] * comps[$i].absorption_coeff))
	end
	push!(expr.args, :(return result))
	return expr
end

"""
	effective_absorption(ξ::NTuple{NC,Float64}, material::Material) -> Float64

Compute effective absorption coefficient [m⁻¹].

α_eff = Σⱼ ξⱼ αⱼ

where αⱼ is the mass absorption coefficient [m²/kg].
"""
@inline function effective_absorption(ξ::NTuple{NC, Float64}, material) where NC
	return _sum_absorption(ξ, material.components, Val(NC))
end

# ============================================================================
# Gas transfer helpers
# ============================================================================

@inline function _gas_transfer_sums(v, T, comps::Tuple{C}, ::Val{1}) where C
	λ = comps[1].gas_transfer_coeff(T)
	λ_parallel = v[1] * λ
	# For SERIES: if component has volume but zero transport, it blocks flow (infinite resistance)
	inv_λ_series = if v[1] > 1e-20
		λ > 1e-20 ? v[1] / λ : Inf
	else
		0.0
	end
	return λ_parallel, inv_λ_series
end

@inline function _gas_transfer_sums(v, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	λ1 = comps[1].gas_transfer_coeff(T)
	λ2 = comps[2].gas_transfer_coeff(T)

	λ_parallel = v[1] * λ1 + v[2] * λ2

	# For SERIES: if component has volume but zero transport, it blocks flow (infinite resistance)
	inv_λ_series = 0.0
	if v[1] > 1e-20
		if λ1 > 1e-20
			inv_λ_series += v[1] / λ1
		else
			inv_λ_series = Inf  # Zero λ with nonzero volume blocks transport
		end
	end
	if v[2] > 1e-20 && !isinf(inv_λ_series)
		if λ2 > 1e-20
			inv_λ_series += v[2] / λ2
		else
			inv_λ_series = Inf  # Zero λ with nonzero volume blocks transport
		end
	end

	return λ_parallel, inv_λ_series
end

@inline function _gas_transfer_sums(v, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	λ1 = comps[1].gas_transfer_coeff(T)
	λ2 = comps[2].gas_transfer_coeff(T)
	λ3 = comps[3].gas_transfer_coeff(T)

	λ_parallel = v[1] * λ1 + v[2] * λ2 + v[3] * λ3

	# For SERIES: if component has volume but zero transport, it blocks flow (infinite resistance)
	inv_λ_series = 0.0
	if v[1] > 1e-20
		if λ1 > 1e-20
			inv_λ_series += v[1] / λ1
		else
			inv_λ_series = Inf
		end
	end
	if v[2] > 1e-20 && !isinf(inv_λ_series)
		if λ2 > 1e-20
			inv_λ_series += v[2] / λ2
		else
			inv_λ_series = Inf
		end
	end
	if v[3] > 1e-20 && !isinf(inv_λ_series)
		if λ3 > 1e-20
			inv_λ_series += v[3] / λ3
		else
			inv_λ_series = Inf
		end
	end

	return λ_parallel, inv_λ_series
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _gas_transfer_sums(v, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		λ_parallel = 0.0
		inv_λ_series = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			λ = comps[$i].gas_transfer_coeff(T)
			λ_parallel += v[$i] * λ
			if v[$i] > 1e-20 && !isinf(inv_λ_series)
				if λ > 1e-20
					inv_λ_series += v[$i] / λ
				else
					inv_λ_series = Inf
				end
			end
		end)
	end
	push!(expr.args, :(return λ_parallel, inv_λ_series))
	return expr
end

# Helper to extract gas transfer coefficients as tuple for Bruggeman
@inline function _gas_transfer_tuple(T, comps::Tuple{C}, ::Val{1}) where C
	return (comps[1].gas_transfer_coeff(T),)
end

@inline function _gas_transfer_tuple(T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return (comps[1].gas_transfer_coeff(T), comps[2].gas_transfer_coeff(T))
end

@inline function _gas_transfer_tuple(T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return (comps[1].gas_transfer_coeff(T),
		comps[2].gas_transfer_coeff(T),
		comps[3].gas_transfer_coeff(T))
end

@generated function _gas_transfer_tuple(T, comps::Tuple, ::Val{NC}) where NC
	exprs = [:(comps[$i].gas_transfer_coeff(T)) for i in 1:NC]
	return :(($(exprs...),))
end

"""
	effective_gas_transfer(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective gas transfer coefficient [m²/s].

Uses the same mixing rules as conductivity:
- PARALLEL: λ = Σⱼ vⱼ λⱼ
- SERIES: λ = 1 / Σⱼ (vⱼ / λⱼ)
- WEIGHTED: λ = β·λ_parallel + (1-β)·λ_series
- BRUGGEMAN: Σⱼ vⱼ(λⱼ - λ_eff)/(λⱼ + 2λ_eff) = 0
"""
@inline function effective_gas_transfer(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	v = volume_fractions(ξ, T, material)

	mixing = material.gas_transfer_mixing
	if mixing == BRUGGEMAN
		λ_tuple = _gas_transfer_tuple(T, material.components, Val(NC))
		return bruggeman_effective_property(v, λ_tuple)
	end

	λ_parallel, inv_λ_series = _gas_transfer_sums(v, T, material.components, Val(NC))
	λ_series = inv_λ_series > 1e-20 ? 1.0 / inv_λ_series : 0.0

	# Apply mixing rule
	if mixing == PARALLEL
		return λ_parallel
	elseif mixing == SERIES
		return λ_series
	else  # WEIGHTED
		β = material.β_λ
		return β * λ_parallel + (1.0 - β) * λ_series
	end
end

# ============================================================================
# Effective Permeability (Component-Based)
# ============================================================================

# Helper: compute parallel and series sums for permeability
@inline function _permeability_sums(v, T, comps::Tuple{C}, ::Val{1}) where C
	κ = get_permeability(comps[1], T)
	κ_parallel = v[1] * κ
	inv_κ_series = v[1] > 1e-20 && κ > 1e-30 ? v[1] / κ : 0.0
	return κ_parallel, inv_κ_series
end

@inline function _permeability_sums(v, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	κ1 = get_permeability(comps[1], T)
	κ2 = get_permeability(comps[2], T)

	κ_parallel = v[1] * κ1 + v[2] * κ2

	inv_κ_series = 0.0
	if v[1] > 1e-20 && κ1 > 1e-30
		inv_κ_series += v[1] / κ1
	end
	if v[2] > 1e-20 && κ2 > 1e-30
		inv_κ_series += v[2] / κ2
	end

	return κ_parallel, inv_κ_series
end

@inline function _permeability_sums(v, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	κ1 = get_permeability(comps[1], T)
	κ2 = get_permeability(comps[2], T)
	κ3 = get_permeability(comps[3], T)

	κ_parallel = v[1] * κ1 + v[2] * κ2 + v[3] * κ3

	inv_κ_series = 0.0
	if v[1] > 1e-20 && κ1 > 1e-30
		inv_κ_series += v[1] / κ1
	end
	if v[2] > 1e-20 && κ2 > 1e-30
		inv_κ_series += v[2] / κ2
	end
	if v[3] > 1e-20 && κ3 > 1e-30
		inv_κ_series += v[3] / κ3
	end

	return κ_parallel, inv_κ_series
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _permeability_sums(v, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		κ_parallel = 0.0
		inv_κ_series = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			κ = get_permeability(comps[$i], T)
			κ_parallel += v[$i] * κ
			if v[$i] > 1e-20 && κ > 1e-30
				inv_κ_series += v[$i] / κ
			end
		end)
	end
	push!(expr.args, :(return κ_parallel, inv_κ_series))
	return expr
end

# Helper to extract permeabilities as tuple for Bruggeman
@inline function _permeability_tuple(T, comps::Tuple{C}, ::Val{1}) where C
	return (get_permeability(comps[1], T),)
end

@inline function _permeability_tuple(T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	return (get_permeability(comps[1], T), get_permeability(comps[2], T))
end

@inline function _permeability_tuple(T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	return (get_permeability(comps[1], T),
		get_permeability(comps[2], T),
		get_permeability(comps[3], T))
end

@generated function _permeability_tuple(T, comps::Tuple, ::Val{NC}) where NC
	exprs = [:(get_permeability(comps[$i], T)) for i in 1:NC]
	return :(($(exprs...),))
end

"""
	effective_permeability(ξ::NTuple{NC,Float64}, T::Float64, material) -> Float64

Compute effective permeability [m²] using component mixing rules.

Only solid and liquid components contribute to permeability (they form the matrix).
Gases flow THROUGH the matrix and do not contribute (has_permeability returns false
for gases, and get_permeability returns 0.0).

Uses the same mixing rules as conductivity:
- PARALLEL: κ = Σⱼ vⱼ κⱼ (upper bound)
- SERIES: κ = 1 / Σⱼ (vⱼ / κⱼ) (lower bound)
- WEIGHTED: κ = β·κ_parallel + (1-β)·κ_series
- BRUGGEMAN: Effective medium theory

# Arguments
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material with components and mixing rules

# Returns
Effective permeability [m²]. Returns 0.0 if no components have permeability.

# Physical Interpretation
Permeability determines resistance to pressure-driven (Darcy) gas flow:
	v_gas = -(κ/μ) ∇P

Typical values:
- Virgin wood: ~10⁻¹⁵ m² (tight cellular structure)
- Char: ~10⁻¹² m² (porous with open channels)
- Impermeable: 0.0 m²

As pyrolysis progresses and virgin material converts to char, the effective
permeability increases, allowing easier gas escape.
"""
@inline function effective_permeability(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	v = volume_fractions(ξ, T, material)

	mixing = material.permeability_mixing
	if mixing == BRUGGEMAN
		κ_tuple = _permeability_tuple(T, material.components, Val(NC))
		return bruggeman_effective_property(v, κ_tuple)
	end

	κ_parallel, inv_κ_series = _permeability_sums(v, T, material.components, Val(NC))
	κ_series = inv_κ_series > 1e-30 ? 1.0 / inv_κ_series : 0.0

	# Apply mixing rule
	if mixing == PARALLEL
		return κ_parallel
	elseif mixing == SERIES
		return κ_series
	else  # WEIGHTED
		β = material.β_κ
		return β * κ_parallel + (1.0 - β) * κ_series
	end
end

"""
	permeability_derivatives(ξ::NTuple{NC,Float64}, T::Float64, material) -> (∂κ/∂T, ∂κ/∂ξ)

Compute derivatives of effective permeability for Jacobian assembly.

Uses numerical differentiation of `effective_permeability()` to compute:
- `∂κ/∂T`: Temperature derivative of effective permeability
- `∂κ/∂ξⱼ`: Concentration derivatives (one per component)

This replaces the old chain-rule approach (∂κ/∂φ × ∂φ/∂ξⱼ) which required
a permeability model. The component-based approach directly differentiates
the mixture rule.

# Arguments
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material with components and permeability mixing rules

# Returns
- `∂κ_∂T::Float64`: Temperature derivative [m²/K]
- `∂κ_∂ξ::NTuple{NC,Float64}`: Concentration derivatives [m²/(kg/m³) = m⁵/kg]

# Notes
For gases (which don't contribute to permeability), ∂κ/∂ξⱼ will be small
but may not be exactly zero due to indirect effects on volume fractions.
"""
@inline function permeability_derivatives(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	# Compute volume fractions and their derivatives
	v = volume_fractions(ξ, T, material)
	∂v_∂ξ = volume_fraction_derivatives(ξ, T, material)

	# Get component permeabilities (κⱼ values are independent of concentration)
	κ_tuple = _permeability_tuple(T, material.components, Val(NC))

	mixing = material.permeability_mixing

	if mixing == BRUGGEMAN
		# For Bruggeman, fall back to numerical differentiation
		# TODO: Implement analytical Bruggeman derivatives
		ε_T = 1.0
		κ_plus_T = effective_permeability(ξ, T + ε_T, material)
		κ_minus_T = effective_permeability(ξ, T - ε_T, material)
		∂κ_∂T = (κ_plus_T - κ_minus_T) / (2.0 * ε_T)

		∂κ_∂ξ = ntuple(Val(NC)) do j
			ε_ξ = max(1e-6, abs(ξ[j]) * 1e-6)
			ξ_plus = ntuple(i -> i == j ? ξ[i] + ε_ξ : ξ[i], Val(NC))
			ξ_minus = ntuple(i -> i == j ? max(0.0, ξ[i] - ε_ξ) : ξ[i], Val(NC))
			κ_plus_ξ = effective_permeability(ξ_plus, T, material)
			κ_minus_ξ = effective_permeability(ξ_minus, T, material)
			(κ_plus_ξ - κ_minus_ξ) / (2.0 * ε_ξ)
		end

		return ∂κ_∂T, ∂κ_∂ξ
	end

	# Analytical derivatives for PARALLEL, SERIES, and WEIGHTED

	# PARALLEL: κ = Σⱼ vⱼκⱼ
	# ∂κ/∂ξᵢ = Σⱼ κⱼ × ∂vⱼ/∂ξᵢ
	∂κ_parallel_∂ξ = ntuple(Val(NC)) do i
		deriv = zero(Float64)
		for j in 1:NC
			deriv += κ_tuple[j] * ∂v_∂ξ[j][i]  # ∂vⱼ/∂ξᵢ
		end
		deriv
	end

	# SERIES: κ = 1/(Σⱼ vⱼ/κⱼ)
	# ∂κ/∂ξᵢ = -κ² × Σⱼ (1/κⱼ) × ∂vⱼ/∂ξᵢ
	inv_κ_series = zero(Float64)
	for j in 1:NC
		if κ_tuple[j] > 1e-30
			inv_κ_series += v[j] / κ_tuple[j]
		end
	end
	κ_series = inv_κ_series > 1e-30 ? 1.0 / inv_κ_series : 0.0

	∂κ_series_∂ξ = ntuple(Val(NC)) do i
		if κ_series < 1e-30
			return 0.0
		end
		deriv_sum = zero(Float64)
		for j in 1:NC
			if κ_tuple[j] > 1e-30
				deriv_sum += (1.0 / κ_tuple[j]) * ∂v_∂ξ[j][i]  # ∂vⱼ/∂ξᵢ
			end
		end
		-κ_series^2 * deriv_sum
	end

	# Compute final derivatives based on mixing rule
	if mixing == PARALLEL
		∂κ_∂ξ = ∂κ_parallel_∂ξ
	elseif mixing == SERIES
		∂κ_∂ξ = ∂κ_series_∂ξ
	else  # WEIGHTED
		β = material.β_κ
		∂κ_∂ξ = ntuple(i -> β * ∂κ_parallel_∂ξ[i] + (1.0 - β) * ∂κ_series_∂ξ[i], Val(NC))
	end

	# Temperature derivative: κⱼ(T) may depend on T for temperature-dependent components
	# For now, assume κⱼ is constant (most common case)
	# TODO: Add analytical temperature derivatives if needed
	∂κ_∂T = 0.0

	return ∂κ_∂T, ∂κ_∂ξ
end

"""
	permeability_derivatives!(∂κ_∂ξ::AbstractVector, ξ::NTuple{NC,Float64}, T::Float64, 
							  material, v_work::AbstractVector, κ_work::AbstractVector,
							  ρ_work::AbstractVector) where NC

Compute derivatives of effective permeability in-place (allocation-free).

Same as `permeability_derivatives` but writes to a pre-allocated vector instead of
returning a tuple. This eliminates heap allocations when called in a loop.

# Arguments
- `∂κ_∂ξ`: Output vector of length NC for concentration derivatives [m⁵/kg]
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material with components and permeability mixing rules
- `v_work`: Work vector of length NC for volume fractions
- `κ_work`: Work vector of length NC for component permeabilities
- `ρ_work`: Work vector of length NC for component densities

# Returns
- `∂κ_∂T::Float64`: Temperature derivative [m²/K] (typically 0)

# Notes
The output vector `∂κ_∂ξ` is modified in-place. Work vectors are used for
intermediate calculations and their contents are undefined on exit.
"""
@inline function permeability_derivatives!(∂κ_∂ξ::AbstractVector, ξ::NTuple{NC, Float64},
	T::Float64, material, v_work::AbstractVector, κ_work::AbstractVector,
	ρ_work::AbstractVector) where NC
	# Step 1: Compute densities, volume fractions, and component permeabilities
	v_total = 0.0
	@inbounds for j in 1:NC
		ρ_j = _get_comp_density(T, material.components, j, Val(NC))
		ρ_work[j] = ρ_j
		if ξ[j] > 0.0 && ρ_j > 0.0
			v_total += ξ[j] / ρ_j
		end
		κ_work[j] = _get_comp_permeability(T, material.components, j, Val(NC))
	end

	inv_v_total = v_total > 1e-20 ? 1.0 / v_total : 0.0

	# Compute volume fractions
	@inbounds for j in 1:NC
		ρ_j = ρ_work[j]
		v_work[j] = ξ[j] > 0.0 && ρ_j > 0.0 ? ξ[j] / ρ_j * inv_v_total : 0.0
	end

	# Step 2: Compute κ_series for chain rule
	inv_κ_series = 0.0
	@inbounds for j in 1:NC
		if κ_work[j] > 1e-30
			inv_κ_series += v_work[j] / κ_work[j]
		end
	end
	κ_series = inv_κ_series > 1e-30 ? 1.0 / inv_κ_series : 0.0
	κ_series_sq = κ_series * κ_series

	# Step 3: Compute ∂κ/∂ξᵢ for each component i
	mixing = material.permeability_mixing
	β = material.β_κ

	@inbounds for i in 1:NC
		inv_ρi = ρ_work[i] > 0.0 ? 1.0 / ρ_work[i] : 0.0

		# Accumulate derivatives over all j
		∂κ_parallel_∂ξi = 0.0
		∂κ_series_∂ξi = 0.0

		for j in 1:NC
			inv_ρj = ρ_work[j] > 0.0 ? 1.0 / ρ_work[j] : 0.0
			# ∂vⱼ/∂ξᵢ = (δᵢⱼ/ρⱼ - vⱼ/ρᵢ) × inv_v_total
			if i == j
				dvj_dxi = (inv_ρj - v_work[j] * inv_ρi) * inv_v_total
			else
				dvj_dxi = -v_work[j] * inv_ρi * inv_v_total
			end
			# Parallel contribution: κⱼ × ∂vⱼ/∂ξᵢ
			∂κ_parallel_∂ξi += κ_work[j] * dvj_dxi
			# Series contribution: -κ² / κⱼ × ∂vⱼ/∂ξᵢ
			if κ_work[j] > 1e-30
				∂κ_series_∂ξi -= κ_series_sq / κ_work[j] * dvj_dxi
			end
		end

		# Combine based on mixing rule
		if mixing == PARALLEL
			∂κ_∂ξ[i] = ∂κ_parallel_∂ξi
		elseif mixing == SERIES
			∂κ_∂ξ[i] = ∂κ_series_∂ξi
		else  # WEIGHTED
			∂κ_∂ξ[i] = β * ∂κ_parallel_∂ξi + (1.0 - β) * ∂κ_series_∂ξi
		end
	end

	return 0.0  # ∂κ/∂T (assuming constant component permeabilities)
end

# Generated helper to get component density (avoids heterogeneous tuple access)
@generated function _get_comp_density(T::Float64, comps::Tuple, j::Int, ::Val{NC}) where NC
	# Generate a branch for each possible j value
	branches = Expr[]
	for i in 1:NC
		push!(branches, :(j == $i && return comps[$i].density(T)))
	end
	return Expr(:block, branches..., :(return 0.0))
end

# Generated helper to get component permeability
@generated function _get_comp_permeability(T::Float64, comps::Tuple, j::Int, ::Val{NC}) where NC
	branches = Expr[]
	for i in 1:NC
		push!(branches, :(j == $i && return comps[$i].permeability === nothing ? 0.0 : comps[$i].permeability(T)))
	end
	return Expr(:block, branches..., :(return 0.0))
end

"""
	harmonic_mean_weighted(a::Float64, b::Float64, d_a::Float64, d_b::Float64) -> Float64

Distance-weighted harmonic mean for non-uniform meshes.

For a property φ across a face separating cells L and R with distances d_L and d_R
from the cell centers to the face, the face value using resistance-in-series is:

	φ_face = (d_L + d_R) / (d_L/φ_L + d_R/φ_R)

This reduces to the standard harmonic mean 2ab/(a+b) when d_L = d_R.

# Arguments
- `a`: Property value in left cell
- `b`: Property value in right cell  
- `d_a`: Distance from left cell center to face
- `d_b`: Distance from right cell center to face
"""
@inline function harmonic_mean_weighted(a::Float64, b::Float64, d_a::Float64, d_b::Float64)
	if a < 1e-20 || b < 1e-20
		return 0.0
	end
	return (d_a + d_b) / (d_a / a + d_b / b)
end

"""
	harmonic_mean_product_weighted(a1::Float64, b1::Float64, a2::Float64, b2::Float64, 
								   d_1::Float64, d_2::Float64) -> Float64

Distance-weighted harmonic mean of products for non-uniform meshes.

Computes the harmonic mean of (a1·b1) and (a2·b2) with distance weighting.
Used for face-centered (λ·ρ) in gas transport on non-uniform meshes.

For gas flux J = -λρ ∂Y/∂x, the product (λ·ρ) is the effective diffusion coefficient
for mass fraction transport. At an interface between two cells with different properties,
flux continuity requires: J = -(Y_R - Y_L) / (R_L + R_R) where R_i = d_i / (λ_i·ρ_i).

# Arguments
- `a1`, `b1`: Values to multiply for left cell (e.g., λ_L, ρ_L)
- `a2`, `b2`: Values to multiply for right cell (e.g., λ_R, ρ_R)
- `d_1`: Distance from left cell center to face
- `d_2`: Distance from right cell center to face
"""
@inline function harmonic_mean_product_weighted(a1::Float64, b1::Float64,
	a2::Float64, b2::Float64,
	d_1::Float64, d_2::Float64)
	p1 = a1 * b1
	p2 = a2 * b2
	if p1 < 1e-20 || p2 < 1e-20
		return 0.0
	end
	return (d_1 + d_2) / (d_1 / p1 + d_2 / p2)
end
# ============================================================================
# Species-Specific Transport Scaling
# ============================================================================

# Default reference molecular weight (air ≈ 29 g/mol)
const M_REF_TRANSPORT = 0.029  # kg/mol

"""
	species_transport_scaling(M_j::Float64; M_ref::Float64=M_REF_TRANSPORT) -> Float64

Compute molecular weight scaling factor for species-specific gas transport.

Gas diffusivity scales with molecular weight according to kinetic theory:
- Knudsen diffusion (small pores): D ∝ √(T/M)
- Molecular diffusion: D_AB ∝ √(1/M_A + 1/M_B)

For simplified species-specific transport, we use:
	λⱼ = λ_eff × √(M_ref / Mⱼ)

This means:
- Light species (M < M_ref) diffuse FASTER (scaling > 1)
- Heavy species (M > M_ref) diffuse SLOWER (scaling < 1)

# Arguments
- `M_j`: Molar mass of species j [kg/mol]
- `M_ref`: Reference molar mass [kg/mol] (default: 0.029, air)

# Returns
Dimensionless scaling factor for transport coefficient.

# Examples
```julia
# Water vapor (M = 18 g/mol) diffuses faster
species_transport_scaling(0.018)  # ≈ 1.27

# MMA (M = 100 g/mol) diffuses slower  
species_transport_scaling(0.100)  # ≈ 0.54

# CO₂ (M = 44 g/mol) diffuses slightly slower
species_transport_scaling(0.044)  # ≈ 0.81
```
"""
@inline function species_transport_scaling(M_j::Float64; M_ref::Float64 = M_REF_TRANSPORT)
	if M_j < 1e-6
		# Non-physical molar mass (probably a solid), return 1.0
		return 1.0
	end
	return sqrt(M_ref / M_j)
end

"""
	species_transport_scaling(comp::Component; M_ref::Float64=M_REF_TRANSPORT) -> Float64

Compute scaling factor for a component's transport.

Only applies to GAS phase components with non-zero molar mass.
Returns 1.0 for solids and liquids (their λ contributes to effective λ but doesn't 
represent actual diffusion).
"""
@inline function species_transport_scaling(comp; M_ref::Float64 = M_REF_TRANSPORT)
	if comp.phase != GAS || comp.molar_mass < 1e-6
		return 1.0
	end
	return species_transport_scaling(comp.molar_mass; M_ref = M_ref)
end

# ============================================================================
# Permeability Utility Functions
# ============================================================================

"""
	harmonic_mean_permeability(κ_L::Float64, κ_R::Float64) -> Float64

Compute harmonic mean of permeabilities for face-centered value.

Harmonic mean is appropriate for resistances in series (flux through
heterogeneous media). This is the correct averaging for Darcy flow
across a material interface.

# Arguments
- `κ_L`: Permeability in left cell [m²]
- `κ_R`: Permeability in right cell [m²]

# Returns
Harmonic mean permeability [m²]. Returns 0.0 if either input is near zero.

# Notes
The harmonic mean gives: 2κ_Lκ_R / (κ_L + κ_R)

This reduces to the geometric mean for equal values and approaches
the smaller value when there's a large mismatch (conservative behavior).
"""
@inline function harmonic_mean_permeability(κ_L::Float64, κ_R::Float64)
	if κ_L < 1e-30 || κ_R < 1e-30
		return 0.0
	end
	return 2.0 * κ_L * κ_R / (κ_L + κ_R)
end

"""
	harmonic_mean_permeability_weighted(κ_L::Float64, κ_R::Float64,
										d_L::Float64, d_R::Float64) -> Float64

Distance-weighted harmonic mean for permeability (resistance in series).

This is the correct formulation for non-uniform meshes where cell centers
are not equidistant from the face. The derivation comes from treating
permeability as resistances in series:

	Total resistance = R_L + R_R = d_L/κ_L + d_R/κ_R

Solving for the effective permeability at the face:

	κ_face = (d_L + d_R) / (d_L/κ_L + d_R/κ_R)

# Arguments
- `κ_L, κ_R`: Permeability in left and right cells [m²]
- `d_L, d_R`: Distance from face to left and right cell centers [m]

# Returns
- `κ_face`: Face-averaged permeability [m²]

# Notes
- For uniform meshes (d_L = d_R), this reduces to the standard harmonic mean
- Most significant improvement occurs on stretched meshes (e.g., boundary refinement)
- Typical accuracy improvement: 1-5% for gas velocities near boundaries

# Reference
Physics Audit Recommendations (2026), Section 3
"""
@inline function harmonic_mean_permeability_weighted(κ_L::Float64, κ_R::Float64,
													 d_L::Float64, d_R::Float64)
	if κ_L < 1e-30 || κ_R < 1e-30
		return 0.0
	end
	# Distance-weighted harmonic mean
	return (d_L + d_R) / (d_L / κ_L + d_R / κ_R)
end

# ============================================================================
# Fused Property Computation
# ============================================================================
#
# Evaluates each component's density, heat_capacity, conductivity, gas_transfer,
# emissivity, and absorption coefficient ONCE, then derives all effective properties
# from those cached values. This eliminates redundant density evaluations that occur
# when each effective_* function independently computes volume fractions.
#
# Before: density(T) called ~5× per component per cell (mixture_density,
#         conductivity, gas_transfer, emissivity, each via volume_fractions)
# After:  density(T) called 1× per component per cell
# ============================================================================

"""
	_compute_all_properties_fused(ξ, T, material) -> (ρ, ρcp, cp, k, λ, ε, α)

Compute all effective material properties in a single fused pass.
Uses @generated to unroll component access at compile time.

Returns: (mixture_density, volumetric_heat_capacity, specific_heat,
          conductivity, gas_transfer, emissivity, absorption)
"""
@generated function _compute_all_properties_fused(ξ::NTuple{NC, Float64}, T::Float64,
	material::Material{NC}) where NC

	exprs = Expr[]

	# Step 1: Evaluate all per-component properties ONCE
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		cp_sym = Symbol("cp_$j")
		k_sym = Symbol("k_$j")
		λ_sym = Symbol("λ_$j")
		push!(exprs, quote
			comp = material.components[$j]
			$ρ_sym = comp.density(T)
			$cp_sym = comp.heat_capacity(T)
			$k_sym = comp.thermal_conductivity(T)
			$λ_sym = comp.gas_transfer_coeff(T)
		end)
	end

	# Step 2: Compute mass fractions and specific volumes
	push!(exprs, quote
		ρ_total = 0.0
	end)
	for j in 1:NC
		push!(exprs, :(ρ_total += ξ[$j]))
	end
	push!(exprs, quote
		inv_ρ_total = ρ_total > 1e-20 ? 1.0 / ρ_total : 0.0
	end)

	# Step 3: Compute ρcp = Σ ξ_j * cp_j and cp = Σ Y_j * cp_j
	push!(exprs, quote
		ρcp_eff = 0.0
	end)
	for j in 1:NC
		cp_sym = Symbol("cp_$j")
		push!(exprs, :(ρcp_eff += ξ[$j] * $cp_sym))
	end
	push!(exprs, quote
		cp_eff = ρ_total > 1e-20 ? ρcp_eff * inv_ρ_total : 0.0
	end)

	# Step 4: Compute α = Σ ξ_j * α_j (absorption, no volume fractions needed)
	push!(exprs, quote
		α_eff = 0.0
	end)
	for j in 1:NC
		push!(exprs, :(α_eff += ξ[$j] * material.components[$j].absorption_coeff))
	end

	# Step 5: Compute mixture density via specific volume with swelling
	push!(exprs, quote
		v_specific = 0.0
	end)
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		Y_sym = Symbol("Y_$j")
		push!(exprs, quote
			$Y_sym = ξ[$j] * inv_ρ_total
			if material.components[$j].phase == GAS
				v_specific += material.components[$j].swelling_factor * $Y_sym / $ρ_sym
			else
				v_specific += $Y_sym / $ρ_sym
			end
		end)
	end
	push!(exprs, quote
		ρ_mix = v_specific > 1e-20 ? 1.0 / v_specific : 0.0
	end)

	# Step 6: Compute volume fractions (full) for gas_transfer
	push!(exprs, quote
		V_total = 0.0
	end)
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		push!(exprs, quote
			if ξ[$j] > 0.0 && $ρ_sym > 0.0
				V_total += ξ[$j] / $ρ_sym
			end
		end)
	end

	# Volume fractions
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		v_sym = Symbol("v_$j")
		push!(exprs, quote
			$v_sym = (V_total > 1e-20 && ξ[$j] > 0.0 && $ρ_sym > 0.0) ? (ξ[$j] / $ρ_sym) / V_total : 0.0
		end)
	end

	# Step 7: Compute solid volume fractions (excludes gas with γ≈0) for conductivity & emissivity
	push!(exprs, quote
		V_solid_total = 0.0
	end)
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		push!(exprs, quote
			if !(material.components[$j].phase == GAS && material.components[$j].swelling_factor < 1e-10)
				if ξ[$j] > 0.0 && $ρ_sym > 0.0
					V_solid_total += ξ[$j] / $ρ_sym
				end
			end
		end)
	end

	# Solid volume fractions
	for j in 1:NC
		ρ_sym = Symbol("ρ_$j")
		sv_sym = Symbol("sv_$j")
		push!(exprs, quote
			if !(material.components[$j].phase == GAS && material.components[$j].swelling_factor < 1e-10)
				$sv_sym = (V_solid_total > 1e-20 && ξ[$j] > 0.0 && $ρ_sym > 0.0) ? (ξ[$j] / $ρ_sym) / V_solid_total : 0.0
			else
				$sv_sym = 0.0
			end
		end)
	end

	# Step 8: Compute effective conductivity from solid volume fractions
	push!(exprs, quote
		k_parallel = 0.0
		inv_k_series = 0.0
	end)
	for j in 1:NC
		sv_sym = Symbol("sv_$j")
		k_sym = Symbol("k_$j")
		push!(exprs, quote
			k_parallel += $sv_sym * $k_sym
			if $sv_sym > 1e-20 && $k_sym > 1e-20
				inv_k_series += $sv_sym / $k_sym
			end
		end)
	end
	push!(exprs, quote
		k_series = inv_k_series > 1e-20 ? 1.0 / inv_k_series : 0.0
		mixing_k = material.conductivity_mixing
		if mixing_k == BRUGGEMAN
			# Fall back to standalone function for Bruggeman (iterative solve)
			k_eff = effective_conductivity(ξ, T, material)
		elseif mixing_k == PARALLEL
			k_eff = k_parallel
		elseif mixing_k == SERIES
			k_eff = k_series
		else  # WEIGHTED
			β_k = material.β
			k_eff = β_k * k_parallel + (1.0 - β_k) * k_series
		end
	end)

	# Step 9: Compute effective gas transfer from full volume fractions
	push!(exprs, quote
		λ_parallel = 0.0
		inv_λ_series = 0.0
		λ_series_inf = false
	end)
	for j in 1:NC
		v_sym = Symbol("v_$j")
		λ_sym = Symbol("λ_$j")
		push!(exprs, quote
			λ_parallel += $v_sym * $λ_sym
			if $v_sym > 1e-20 && !λ_series_inf
				if $λ_sym > 1e-20
					inv_λ_series += $v_sym / $λ_sym
				else
					λ_series_inf = true
				end
			end
		end)
	end
	push!(exprs, quote
		λ_series = (!λ_series_inf && inv_λ_series > 1e-20) ? 1.0 / inv_λ_series : 0.0
		mixing_λ = material.gas_transfer_mixing
		if mixing_λ == BRUGGEMAN
			λ_eff = effective_gas_transfer(ξ, T, material)
		elseif mixing_λ == PARALLEL
			λ_eff = λ_parallel
		elseif mixing_λ == SERIES
			λ_eff = λ_series
		else  # WEIGHTED
			β_λ = material.β_λ
			λ_eff = β_λ * λ_parallel + (1.0 - β_λ) * λ_series
		end
	end)

	# Step 10: Compute effective emissivity from solid volume fractions
	push!(exprs, quote
		ε_eff = 0.0
	end)
	for j in 1:NC
		sv_sym = Symbol("sv_$j")
		push!(exprs, :(ε_eff += $sv_sym * material.components[$j].emissivity))
	end
	push!(exprs, :(ε_eff = clamp(ε_eff, 0.0, 1.0)))

	# Return all properties
	push!(exprs, :(return (ρ_mix, ρcp_eff, cp_eff, k_eff, λ_eff, ε_eff, α_eff)))

	return Expr(:block, exprs...)
end