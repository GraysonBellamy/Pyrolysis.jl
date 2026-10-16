using Test
using Pyrolysis

# `stretch_from` names the boundary that gets the FINEST cells (default :top
# puts the finest cells at the heated surface).

_node_positions(mesh) = [node.z for node in mesh.nodes]
_cell_widths(mesh) = diff(_node_positions(mesh))

@testset "Mesh Generation" begin
    @testset "Stretched mesh orientation" begin
        L = 0.01
        n = 50
        s = 1.1

        @testset "stretch_from = :top refines the exposed surface" begin
            mesh = create_stretched_mesh(L, n, Val(2), s; stretch_from = :top)
            Δz = _cell_widths(mesh)

            @test length(Δz) == n
            @test all(>(0.0), Δz)
            @test sum(Δz) ≈ L rtol = 1e-12
            @test _node_positions(mesh)[end] == L

            # Finest cell at z = L (top, cell n), coarsest at z = 0 (bottom)
            @test Δz[end] ≈ minimum(Δz)
            @test Δz[1] ≈ maximum(Δz)
            @test Δz[1] / Δz[end] ≈ s^(n - 1) rtol = 1e-10

            # Adjacent-cell ratio is the stretch factor, growing away from
            # the refined (top) end: Δz_i / Δz_{i+1} = s
            @test all(Δz[i]/Δz[i+1] ≈ s for i = 1:(n-1))
        end

        @testset "stretch_from = :bottom refines the substrate" begin
            mesh = create_stretched_mesh(L, n, Val(2), s; stretch_from = :bottom)
            Δz = _cell_widths(mesh)

            @test sum(Δz) ≈ L rtol = 1e-12
            @test Δz[1] ≈ minimum(Δz)
            @test Δz[end] ≈ maximum(Δz)
            @test all(Δz[i+1]/Δz[i] ≈ s for i = 1:(n-1))
        end

        @testset "the two orientations are mirror images" begin
            Δz_top =
                _cell_widths(create_stretched_mesh(L, n, Val(2), s; stretch_from = :top))
            Δz_bot =
                _cell_widths(create_stretched_mesh(L, n, Val(2), s; stretch_from = :bottom))
            @test Δz_top ≈ reverse(Δz_bot)
        end

        @testset "boundary tags are position-correct after reflection" begin
            mesh = create_stretched_mesh(L, n, Val(2), s; stretch_from = :top)
            z = _node_positions(mesh)
            @test issorted(z)
            # boundary_id convention: 1 = top (z = L), 2 = bottom (z = 0)
            @test mesh.nodes[1].boundary_id == 2
            @test mesh.nodes[1].z == 0.0
            @test mesh.nodes[end].boundary_id == 1
            @test mesh.nodes[end].z == L
        end

        @testset "stretch ≈ 1 falls back to uniform" begin
            mesh = create_stretched_mesh(L, n, Val(2), 1.0)
            Δz = _cell_widths(mesh)
            @test all(w ≈ L / n for w in Δz)
        end

        @testset "argument validation" begin
            @test_throws ArgumentError create_stretched_mesh(
                L,
                n,
                Val(2),
                s;
                stretch_from = :left,
            )
            # Orientation is selected by stretch_from only, never by a
            # sub-unity ratio
            @test_throws ArgumentError create_stretched_mesh(L, n, Val(2), 0.9)
        end
    end
end
