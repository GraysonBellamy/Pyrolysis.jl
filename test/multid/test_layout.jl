using Test
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length

# The line layout's indices in closed form: node n's block at (n − 1) s, cell
# i's at (i − 1) s + w_node, the globals after the last node, with the stride
# s = w_node + w_cell; a node block holds x, w and I, a cell block the masses,
# the energy state and the history pair.
function _line_indices(NC, n, moving, intensities, chi)
    wn = 2moving + intensities
    wc = NC + 1 + 2chi
    s = wn + wc
    node(nd) = (nd - 1) * s
    cell(i) = (i - 1) * s + wn
    rows = Dict{Any,Int}()
    for i = 1:n
        for j = 1:NC
            rows[(:mass, i, j)] = cell(i) + j
        end
        rows[(:energy, i)] = cell(i) + NC + 1
        chi && (rows[(:G, i)] = cell(i) + NC + 2; rows[(:M, i)] = cell(i) + NC + 3)
    end
    for nd = 1:(n+1)
        moving && (rows[(:x, nd)] = node(nd) + 1; rows[(:w, nd)] = node(nd) + 2)
        intensities && (rows[(:I, nd)] = node(nd) + 2moving + 1)
    end
    chi && (rows[(:A,)] = n * s + wn + 1; rows[(:θ,)] = n * s + wn + 2)
    return rows, n * s + wn + 2chi
end

function _indices(l::StateLayout, NC, n, moving, intensities, chi)
    rows = Dict{Any,Int}()
    for i = 1:n
        for j = 1:NC
            rows[(:mass, i, j)] = mass_index(l, i, j)
        end
        rows[(:energy, i)] = energy_index(l, i)
        chi &&
            (rows[(:G, i)] = history_index(l, i, 1); rows[(:M, i)] = history_index(l, i, 2))
    end
    for nd = 1:(n+1)
        moving &&
            (rows[(:x, nd)] = position_index(l, nd); rows[(:w, nd)] = velocity_index(l, nd))
        intensities && (rows[(:I, nd)] = intensity_index(l, nd))
    end
    chi && (rows[(:A,)] = area_index(l); rows[(:θ,)] = lateral_rate_index(l))
    return rows
end

@testset "line layouts keep every index of the interleaved closed form" begin
    for NC in (1, 3),
        n in (1, 4),
        energy in (TemperatureState(), EnthalpyState()),
        (moving, intensities, chi) in (
            (false, false, false),
            (true, false, false),
            (false, true, false),
            (true, true, false),
            (true, true, true),
            (true, false, true),
        )

        expected, len = _line_indices(NC, n, moving, intensities, chi)
        topo = line_mesh(collect(range(0.0, 1.0; length = n + 1))).topology
        for l in (
            StateLayout(NC, n, energy; moving, intensities, chi),
            StateLayout(NC, topo, energy; moving, intensities, chi),
        )
            @test _indices(l, NC, n, moving, intensities, chi) == expected
            @test state_length(l) == len
            @test dimension(l) == 1
            @test closure_unknowns(l) === :velocities
            for nd = 1:(n+1)
                moving && @test position_index(l, nd, 1) == position_index(l, nd)
                moving && @test velocity_index(l, nd, 1) == velocity_index(l, nd)
            end
        end
    end
end

@testset "block order: interleaved on a line, cell by cell in more dimensions" begin
    topo = line_mesh([0.0, 1.0, 2.0]).topology
    @test block_order(topo) == [(:node, 1), (:cell, 1), (:node, 2), (:cell, 2), (:node, 3)]
    topo2 = box_mesh(0:1.0:2.0, 0:1.0).topology
    @test block_order(topo2) == [
        (:node, 1),
        (:node, 2),
        (:node, 4),
        (:node, 5),
        (:cell, 1),
        (:node, 3),
        (:node, 6),
        (:cell, 2),
    ]
end

@testset "$(D)D layouts: the rows partition the state, $(unknowns)" for D in (2, 3),
    unknowns in (:velocities, :multipliers)

    mesh = D == 2 ? box_mesh(0:0.5:1.0, 0:0.25:1.0) : box_mesh(0:0.5:1.0, 0:1.0, 0:0.5:1.0)
    topo = mesh.topology
    NC = 2
    for moving in (false, true)
        l = StateLayout(NC, topo, EnthalpyState(); moving, unknowns)
        @test dimension(l) == D
        @test n_nodes(l) == n_nodes(topo)
        rows = Int[]
        for i = 1:n_cells(topo)
            append!(rows, (mass_index(l, i, j) for j = 1:NC))
            push!(rows, energy_index(l, i))
            moving && unknowns === :multipliers && push!(rows, multiplier_index(l, i))
        end
        for nd = 1:n_nodes(topo), k = 1:D
            moving && push!(rows, position_index(l, nd, k))
            moving && unknowns === :velocities && push!(rows, velocity_index(l, nd, k))
        end
        @test sort(rows) == collect(1:state_length(l))
        d = is_differential(l)
        algebraic =
            moving ? (unknowns === :velocities ? D * n_nodes(topo) : n_cells(topo)) : 0
        @test count(!, d) == algebraic
        # A cell's block follows the blocks of its nodes.
        for i = 1:n_cells(topo), nd in topo.cells[i].nodes
            moving && @test position_index(l, nd, D) < mass_index(l, i, 1)
        end
        if moving && unknowns === :multipliers
            @test_throws ArgumentError velocity_index(l, 1, 1)
            @test !d[multiplier_index(l, 1)]
        elseif moving
            @test_throws ArgumentError multiplier_index(l, 1)
            @test !d[velocity_index(l, 1, D)] && d[position_index(l, 1, D)]
        end
    end
    @test_throws ArgumentError StateLayout(NC, topo, TemperatureState(); intensities = true)
    @test_throws ArgumentError StateLayout(
        NC,
        topo,
        TemperatureState();
        moving = true,
        unknowns = :potential,
    )
end

@testset "carried fields: the history pair and the deformation gradient, $(D)D" for D in
                                                                                    (2, 3)
    mesh = D == 2 ? box_mesh(0:0.5:1.0, 0:0.5:1.0) : box_mesh(0:0.5:1.0, 0:1.0, 0:0.5:1.0)
    topo = mesh.topology
    NC = 2
    for (history, deformation) in ((true, false), (false, true), (true, true)),
        unknowns in (:velocities, :multipliers)

        l = StateLayout(
            NC,
            topo,
            TemperatureState();
            moving = true,
            history,
            deformation,
            unknowns,
        )
        @test has_history(l) == history
        @test has_deformation(l) == deformation
        @test !has_chi(l)
        @test carried_fields(l) ==
              ((history ? (:history,) : ())..., (deformation ? (:deformation,) : ())...)
        rows = Int[]
        for i = 1:n_cells(topo)
            append!(rows, (mass_index(l, i, j) for j = 1:NC))
            push!(rows, energy_index(l, i))
            history && append!(rows, (history_index(l, i, 1), history_index(l, i, 2)))
            deformation &&
                append!(rows, (deformation_index(l, i, a, b) for b = 1:D for a = 1:D))
            unknowns === :multipliers && push!(rows, multiplier_index(l, i))
            # The cell block in order: masses, energy, history, F by columns, φ.
            block = cell_rows(l, i)
            @test rows[(end-length(block)+1):end] == collect(block)
        end
        for nd = 1:n_nodes(topo), k = 1:D
            push!(rows, position_index(l, nd, k))
            unknowns === :velocities && push!(rows, velocity_index(l, nd, k))
        end
        @test sort(rows) == collect(1:state_length(l))
        d = is_differential(l)
        deformation && @test d[deformation_index(l, 1, D, D)]
        history && @test d[history_index(l, 1, 2)]
    end
    @test_throws ArgumentError StateLayout(NC, topo, TemperatureState(); deformation = true)
    @test_throws ArgumentError history_index(
        StateLayout(NC, topo, TemperatureState()),
        1,
        1,
    )
end
