"""
AD-vs-FD Jacobian correctness invariant.

Locks in the contract that `:ad` produces the same Jacobian as a centered
finite-difference reference, on every canonical `SimulationMode`. This is the
gate that the pluggable Jacobian backends rely on — once it passes,
`:ad` is the trusted reference for `verify_jacobian` mode.

Threshold rationale:
- AD is exact up to floating-point roundoff (~1e-14 relative).
- Centered FD has truncation error O(h²); with adaptive h ~ √eps ≈ 1.5e-8
  on Float64, the leading error is ~h² ≈ 2e-16 plus condition-number
  amplification. We target 1e-6 relative on the columns (ξ blocks have
  small magnitudes — relative error there is dominated by FD noise).
"""

using Test
using LinearAlgebra
using Pyrolysis
using Pyrolysis: JacobianSpec, DenseAD_ForwardDiff, compute_jacobian!

include(joinpath(@__DIR__, "perf_helpers.jl"))

# Centered finite-difference reference. Uses adaptive step size keyed to the
# magnitude of each component to handle the wide variable-magnitude span of
# the ALE state vector (T ~ 300 K, ξ ~ 1000 kg/m³, z ~ 0.001 m, χ ~ 1).
function _fd_jacobian!(J::AbstractMatrix, def, ws, u::Vector{Float64}, t::Float64)
    n = length(u)
    du_p = similar(u)
    du_m = similar(u)
    u_pert = copy(u)
    @inbounds for j = 1:n
        h = max(1e-7 * abs(u[j]), 1e-9)
        u_pert[j] = u[j] + h
        residual!(du_p, u_pert, def, t, ws)
        u_pert[j] = u[j] - h
        residual!(du_m, u_pert, def, t, ws)
        u_pert[j] = u[j]
        inv_2h = 1.0 / (2.0 * h)
        for i = 1:n
            J[i, j] = (du_p[i] - du_m[i]) * inv_2h
        end
    end
    return J
end

@testset "AD vs FD Jacobian correctness" begin
    configs = [
        ("Eulerian / NoRadiation / Fickian / NoChi", build_baseline_eulerian, false),
        ("ALE / NoRadiation / Fickian / NoChi", build_ale_no_radiation, true),
    ]

    for (label, builder, is_ale) in configs
        @testset "$label" begin
            def, ws, u, _ = builder()
            n = length(u)
            t = 0.0

            J_fd = zeros(n, n)
            _fd_jacobian!(J_fd, def, ws, u, t)

            J_ad = zeros(n, n)
            spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
            compute_jacobian!(J_ad, spec, u, def, t, ws)

            scale = max(maximum(abs, J_fd), 1e-12)
            abs_err = maximum(abs, J_ad .- J_fd)
            rel_err = abs_err / scale
            @test rel_err < 1e-6
        end
    end
end
