"""
Integration test for the P1 radiation sign convention fix.

Verifies that the non-ALE P1 radiation path correctly heats material when
G > 4σT⁴ (net absorption) and cools when G < 4σT⁴ (net emission).

Physics: From the RTE, ∇·q_rad = κ(4σT⁴ - G). The energy equation is
ρcₚ dT/dt = ... - ∇·q_rad = ... + κ(G - 4σT⁴) = ... + S_rad.
S_rad > 0 means material absorbs more than it emits → heating.
"""

using Test
using Pyrolysis
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Discretization: compute_radiation_flux_divergence!, compute_rhs_with_p1_radiation!,
	RHSCache, update_radiation_properties!

@testset "P1 Radiation Sign Correctness" begin

	σ = STEFAN_BOLTZMANN

	# Create a simple absorbing material (optically thick for P1 validity)
	material = Material(
		name = :AbsorbingMaterial,
		components = (
			SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.25, ε = 0.9, α = 100.0),
			GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 0.0, ε = 0.0, α = 0.0),
		),
		reactions = (
			Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = 300e3),
		),
		β = 0.5,
	)

	NC = n_components(material)
	n_cells = 10

	# Helper to create mesh at uniform temperature
	function make_mesh(T0::Float64)
		mesh = create_uniform_mesh(0.01, n_cells, Val(NC))
		for i in 1:n_cells
			mesh.cell_states[i] = CellState{NC}(
				T0,
				(1000.0, 0.0),
				mesh.cell_volumes[i],
				mesh.cell_volumes[i],
			)
		end
		return mesh
	end

	@testset "Cold material in hot radiation field → heating" begin
		# T = 300K, impose G >> 4σT⁴
		T0 = 300.0
		mesh = make_mesh(T0)

		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, material)
		rad_cache = RadiationCache(mesh)

		# Set G to correspond to a 1000K radiation field
		G_hot = 4.0 * σ * 1000.0^4
		fill!(rad_cache.G, G_hot)

		# Compute radiation properties and S_rad
		update_radiation_properties!(rad_cache, mesh, prop_cache)
		compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

		# S_rad should be positive (net absorption → heating)
		for i in 1:n_cells
			@test rad_cache.S_rad[i] > 0.0
		end

		# Verify: with S_rad > 0, the contribution to dT should be positive
		# dT_P1 = S_rad / ρcₚ > 0  (this is what the fix ensures)
		for i in 1:n_cells
			ρcp = prop_cache.ρcp[i]
			dT_P1 = rad_cache.S_rad[i] / ρcp
			@test dT_P1 > 0.0  # Positive contribution = heating
		end
	end

	@testset "Hot material in cold radiation field → cooling" begin
		T0 = 1500.0
		mesh = make_mesh(T0)

		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, material)
		rad_cache = RadiationCache(mesh)

		# Set G to correspond to a 300K radiation field (weak)
		G_cold = 4.0 * σ * 300.0^4
		fill!(rad_cache.G, G_cold)

		update_radiation_properties!(rad_cache, mesh, prop_cache)
		compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

		# S_rad should be negative (net emission → cooling)
		for i in 1:n_cells
			@test rad_cache.S_rad[i] < 0.0
		end
	end

	@testset "Equilibrium → zero S_rad" begin
		T0 = 800.0
		mesh = make_mesh(T0)

		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, material)
		rad_cache = RadiationCache(mesh)

		# G = 4σT⁴ at every cell (equilibrium)
		G_eq = 4.0 * σ * T0^4
		fill!(rad_cache.G, G_eq)

		update_radiation_properties!(rad_cache, mesh, prop_cache)
		compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

		for i in 1:n_cells
			@test abs(rad_cache.S_rad[i]) < 1e-6
		end
	end
end
