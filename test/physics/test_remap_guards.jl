"""
Remap regression tests.

- A remap is an instantaneous re-grid at fixed time. Mixing A_old (overlap
  volumes) with A_new (division volumes) silently scaled every concentration
  by A_old/A_new — now both meshes must share one cross-section area or the
  remap throws.
- The dead remap layer (`remap_cell_states!`, which ignored its `method`
  argument, `apply_mass_correction!`, `remap_solution_conservative!`) is
  deleted outright.
"""

using Test
using Pyrolysis
using Pyrolysis.Physics: remap_solution!, FirstOrderRemap, SecondOrderRemap

function _remap_meshes(; A_new = 1.0)
    mesh_old = create_uniform_mesh(0.01, 4, Val(1))
    mesh_new = create_uniform_mesh(0.01, 5, Val(1); cross_section_area = A_new)
    return mesh_old, mesh_new
end

@testset "remap validates a shared cross-section" begin
    NC = 1
    # Same-instant re-grid: mass-exact for the first-order method.
    mesh_old, mesh_new = _remap_meshes()
    u_old = zeros(4 * (1 + NC))
    for i = 1:4
        u_old[(i-1)*2+1] = 300.0 + 50.0 * i   # T
        u_old[(i-1)*2+2] = 100.0 * i          # ξ
    end
    u_new = zeros(5 * (1 + NC))
    remap_solution!(u_new, mesh_new, u_old, mesh_old, NC, FirstOrderRemap())
    m_old = sum(u_old[(i-1)*2+2] * mesh_old.cell_volumes[i] for i = 1:4)
    m_new = sum(u_new[(i-1)*2+2] * mesh_new.cell_volumes[i] for i = 1:5)
    @test m_new ≈ m_old rtol = 1e-12

    # Mismatched areas: pre-fix this silently scaled ξ by A_old/A_new;
    # now both remap orders throw.
    mesh_old2, mesh_new2 = _remap_meshes(A_new = 2.0)
    @test_throws ArgumentError remap_solution!(
        zeros(10),
        mesh_new2,
        u_old,
        mesh_old2,
        NC,
        FirstOrderRemap(),
    )
    @test_throws ArgumentError remap_solution!(
        zeros(10),
        mesh_new2,
        u_old,
        mesh_old2,
        NC,
        SecondOrderRemap(),
    )
end

@testset "dead remap layer deleted" begin
    for sym in (:remap_cell_states!, :apply_mass_correction!, :remap_solution_conservative!)
        @test !isdefined(Pyrolysis.Physics, sym)
        @test !isdefined(Pyrolysis, sym)
    end
end
