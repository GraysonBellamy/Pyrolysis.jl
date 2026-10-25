# Figures 9–15: actual constructor coordinates and structural Jacobians,
# with editable vector schematics for geometry, transport and accounting.

function constructor_figure()
    L = 7e-3
    meshes = [
        create_uniform_mesh(L, 12),
        create_stretched_mesh(L, 12, 1.2),
        create_stretched_mesh(L, 12, 1.2; stretch_from = :bottom),
        create_mesh_from_nodes(
            [
                collect(range(0, 6e-3; length = 9));
                collect(range(6e-3, L; length = 5))[2:end]
            ],
        ),
    ]
    names = ["Uniform", "Surface-Refined", "Back-Refined", "Layered"]
    mkpath(joinpath(OUT, "data"))
    open(joinpath(OUT, "data", "constructor_nodes.csv"), "w") do io
        println(io, "mesh,node_id,z_m")
        for (k, mesh) in enumerate(meshes), (i, x) in enumerate(mesh.positions)
            println(io, join((k, i, x[1]), ','))
        end
    end
    @assert meshes[4].positions[9][1] == 6e-3
    svg(
        "09_mesh_constructors",
        "Figure 9. Mesh Constructor Comparison",
        "Actual node coordinates for uniform, surface-refined, back-refined and layered meshes over 7 millimetres. The layered mesh has a node exactly on the 6 millimetre material interface.";
        height = 620,
    ) do io
        for (k, mesh) in enumerate(meshes)
            y = 118 + 108(k - 1)
            label(io, 35, y + 30, names[k]; size = 20)
            xx = [360 + 760x[1] / L for x in mesh.positions]
            fills = k == 4 ? [fill(GREY, 8); fill(GREEN, 4)] : fill(PALE, 12)
            meshstrip(io, xx, y; h = 52, ids = false, fills)
        end
        interface = 360 + 760 * 6 / 7
        line(io, interface, 428, interface, 511; color = TEAL, width = 3)
        label(
            io,
            interface,
            421,
            "Material Interface: 6 mm";
            anchor = "middle",
            color = TEAL,
            size = 18,
        )
        line(io, 360, 530, 1120, 530; color = MUTED)
        for z = 0:7
            x = 360 + 760z / 7
            line(io, x, 525, x, 537; color = MUTED)
            label(io, x, 559, string(z); anchor = "middle", size = 18)
        end
        label(io, 740, 591, "Coordinate z [mm]"; anchor = "middle", size = 19)
    end
end

function geometry_figure()
    svg(
        "10_coordinate_geometry",
        "Figure 10. One-Dimensional Coordinate Geometry (Schematic)",
        "A planar slice, cylindrical radial shell and spherical radial shell illustrate cross section as area, axial length or solid angle. Each model resolves one coordinate only, with face area determined by radius.";
        height = 525,
    ) do io
        for (x, title) in (
            (30, "(a) Planar Slice"),
            (425, "(b) Cylindrical Shell"),
            (820, "(c) Spherical Shell"),
        )
            panel(io, x, 85, 350, 505, title)
        end
        rect(io, 95, 195, 155, 160; fill = PALE)
        path(io, "M 95 195 L 130 165 L 285 165 L 250 195 Z"; fill = GREEN, color = TEAL)
        path(io, "M 250 195 L 285 165 L 285 325 L 250 355 Z"; fill = GREEN, color = TEAL)
        line(io, 85, 395, 290, 395; arrow = "ink")
        label(io, 290, 422, "+z"; anchor = "end")
        label(io, 268, 270, "A"; color = TEAL, size = 23, anchor = "middle")
        lines(
            io,
            48,
            466,
            ["`cross_section` = A [m²]", "Face Area: A"];
            size = 19,
            step = 31,
        )
        # Cross section through a circular annulus; extrusion length is A.
        circle(io, 595, 260, 96; fill = PALE, stroke = BLUE)
        circle(io, 595, 260, 55; fill = "white", stroke = BLUE)
        line(io, 595, 260, 701, 260; arrow = "ink")
        label(io, 666, 245, "r₊"; size = 18)
        label(io, 624, 284, "r₋"; size = 18)
        label(io, 595, 395, "Axial Length A (Out of Plane)"; anchor = "middle", size = 18)
        lines(
            io,
            443,
            466,
            ["`cross_section` = A [m]", "Face Area: 2πr A"];
            size = 19,
            step = 31,
        )
        circle(io, 995, 260, 96; fill = PALE, stroke = BLUE)
        circle(io, 995, 260, 55; fill = "white", stroke = BLUE)
        path(
            io,
            "M 995 260 L 1063 192 A 96 96 0 0 1 1091 260 Z";
            fill = GREEN,
            color = TEAL,
        )
        line(io, 995, 260, 1103, 260; arrow = "ink")
        label(io, 1043, 229, "A"; color = TEAL, size = 23)
        label(io, 1100, 286, "r₊"; size = 18)
        label(io, 1025, 286, "r₋"; size = 18)
        label(
            io,
            1000,
            395,
            "A: Solid Angle of Radial Sector";
            anchor = "middle",
            size = 18,
        )
        lines(
            io,
            838,
            466,
            ["`cross_section` = A [sr]", "Face Area: A r²"];
            size = 19,
            step = 31,
        )
    end
end

function shrinkage_figure()
    svg(
        "11_shrinkage",
        "Figure 11. Axial and Lateral Volume Change (Schematic)",
        "Initial volume compared with axial shrinkage, lateral shrinkage and swelling. Thickness and cross-sectional area are independent geometric factors, and mass-loss rate may be normalized by current or initial area.";
        height = 495,
    ) do io
        for (x, title, w, h, note) in (
            (35, "Initial Geometry", 190, 155, "V₀ = A₀ Δz₀"),
            (330, "Axial Shrinkage", 130, 155, "A = A₀; Δz < Δz₀"),
            (625, "Lateral Shrinkage", 190, 100, "A < A₀"),
            (920, "Swelling", 220, 185, "A and/or Δz increase"),
        )
            label(io, x, 105, title; weight = "bold", size = 21)
            rect(io, x + 10, 150, 190, 155; fill = "none", stroke = MUTED, dash = "6 4")
            rect(io, x + 10, 150 + (155-h)/2, w, h; fill = x == 625 ? GREEN : PALE)
            line(io, x + 10, 340, x + 10 + w, 340; arrow = "ink")
            label(io, x + 10 + w/2, 371, "Δz"; anchor = "middle")
            label(io, x + 10 + w/2, 236, "A"; color = TEAL, anchor = "middle", size = 26)
            label(io, x, 413, note; size = 18)
        end
        label(io, 600, 472, "V = A Δz;  θ = θ_{L} + θ_{A}"; anchor = "middle", size = 23)
    end
end

function gas_figure()
    svg(
        "12_gas_transport",
        "Figure 12. Gas Storage and Transport (Schematic)",
        "A porous cell distinguishes gas mass per bulk volume from pore density. The relative flux follows a lumped or dusty-gas closure; skeleton motion carries pore gas. A boundary diagram matches the porous flux to a reservoir through an optional mass film and accounts for carried enthalpy once.";
        height = 510,
    ) do io
        rect(io, 65, 115, 430, 190; fill = GREY, stroke = MUTED)
        for (x, y, r) in
            ((120, 165, 27), (205, 240, 35), (295, 175, 31), (395, 240, 34), (450, 155, 22))
            circle(io, x, y, r; fill = "white", stroke = TEAL)
        end
        label(io, 205, 247, "Gas"; anchor = "middle", size = 18)
        label(io, 80, 285, "Solid Matrix"; size = 18)
        label(io, 65, 345, "Bulk: ξ = m_{gas}/V"; size = 20)
        label(io, 295, 345, "Pores: ρ_{gas} = ξ/φ"; size = 20)
        for (y, title, color, marker) in (
            (150, "Wall / Molecular Friction", BLUE, "blue"),
            (215, "Darcy Flow", TEAL, "teal"),
            (280, "Skeleton Motion vₛ", PURPLE, "purple"),
        )
            line(io, 540, y, 805, y; color, arrow = marker, width = 2.2)
            label(io, 830, y+7, title; color)
        end
        label(io, 600, 345, "Relative to Mesh: vₛ − w = −c"; size = 20)
        line(io, 35, 385, 1165, 385; color = MUTED, width = 0.8)
        label(io, 65, 430, "One Transmitted Species Flux"; weight = "bold")
        label(io, 65, 478, "Porous Medium"; color = TEAL)
        label(io, 290, 478, "↔"; size = 26)
        label(io, 340, 478, "Optional Film"; color = BLUE)
        line(io, 520, 472, 675, 472; arrow = "ink")
        label(io, 700, 478, "Reservoir: P, T, Composition")
    end
end

function reactions_figure()
    svg(
        "13_reaction_network",
        "Figure 13. Three-Step Wood Decomposition (Schematic)",
        "Three single-reactant steps convert virgin wood through two intermediates to char. Each arrow labels a mass yield; gas branches and signed reaction heats identify the products and thermal effects.";
        height = 455,
    ) do io
        label(io, 35, 83, "Mass Yields [kg/kg Reactant]"; size = 19)
        for (x, name) in (
            (35, "Virgin Wood"),
            (340, "Intermediate 1"),
            (645, "Intermediate 2"),
            (950, "Char"),
        )
            rect(io, x, 120, 200, 60; fill = "white")
            label(io, x+100, 157, name; anchor = "middle", weight = "bold")
        end
        for (k, solid, gas, heat, sign) in (
            (1, "0.95", "0.05", "−19.2", "Exothermic"),
            (2, "0.80", "0.20", "+264", "Endothermic"),
            (3, "0.74", "0.26", "−0.33", "Exothermic"),
        )
            x = 235 + 305(k-1)
            line(io, x+5, 150, x+95, 150; arrow = "ink")
            label(io, x+50, 130, solid; anchor = "middle", size = 18)
            line(io, x+50, 155, x+50, 275; arrow = "teal", color = TEAL)
            label(io, x+65, 225, gas; color = TEAL)
            rect(io, x-50, 285, 200, 60; fill = GREEN, stroke = TEAL)
            label(io, x+50, 322, "Gas $k"; anchor = "middle", weight = "bold")
            label(io, x+50, 390, "h_{$k} = $heat kJ/kg"; anchor = "middle")
            label(
                io,
                x+50,
                425,
                sign;
                anchor = "middle",
                color = k == 2 ? ORANGE : PURPLE,
                size = 19,
            )
        end
    end
end

function conservation_figure()
    svg(
        "14_conservation",
        "Figure 14. Conservation of the Live Inventory (Schematic)",
        "Mass and sensible enthalpy are counted over live cells. Internal transfers cancel; boundary outflow and event losses leave the domain. The energy sources include reaction heat, radiation and sensible-enthalpy formation.";
        height = 525,
    ) do io
        panel(io, 35, 80, 530, 410, "(a) Mass [kg]")
        panel(io, 615, 80, 550, 410, "(b) Sensible Enthalpy [J]")
        for (x, symbol, color, marker) in
            ((110, "M", TEAL, "teal"), (700, "H", ORANGE, "orange"))
            rect(io, x, 200, 400, 145; fill = PALE, stroke = MUTED)
            label(
                io,
                x+200,
                233,
                "Live Inventory $symbol";
                anchor = "middle",
                weight = "bold",
            )
            line(io, x+200, 250, x+200, 335; color = MUTED, dash = "4 4")
            line(io, x+100, 275, x+300, 275; color, arrow = marker)
            flux = symbol == "M" ? "F_{m} A" : "F_{E} A"
            label(io, x+100, 319, "−$flux"; anchor = "middle", color, size = 21)
            label(io, x+300, 319, "+$flux"; anchor = "middle", color, size = 21)
            line(io, x+405, 265, x+465, 265; color, arrow = marker, width = 2)
            label(io, x+380, 179, "Boundary Outflow"; anchor = "middle", size = 18)
            line(io, x+320, 350, x+320, 423; color = PURPLE, arrow = "purple", width = 2)
            label(io, x+320, 454, "Event Loss"; anchor = "middle", color = PURPLE)
        end
        label(io, 235, 395, "Σⱼ S_{mass,ij} = 0"; anchor = "middle", size = 21)
        label(io, 235, 427, "Mass-Balanced Reactions"; anchor = "middle", size = 18)
        line(io, 800, 420, 800, 350; color = ORANGE, arrow = "orange", width = 2)
        label(io, 790, 454, "Reaction + Radiation"; anchor = "middle", size = 18)
        label(io, 790, 489, "+ Sensible-Enthalpy Formation"; anchor = "middle", size = 18)
    end
end

function jacobian_figure()
    cases = [
        ("Lagrangian", NoCorrection(), nothing),
        ("Corrected ALE", EquidistributionCorrection(), nothing),
        ("Variable-Area ALE", EquidistributionCorrection(), x -> 1 - 0.3x),
    ]
    cps = [
        setup(
            small_problem(; law);
            use_ale = true,
            radiation_model = BEER_LAMBERT,
            correction = corr,
        ) for (_, corr, law) in cases
    ]
    jacobians = [dense_jacobian(cp.def, cp.ws, cp.u0, 0.0) for cp in cps]
    for (cp, J) in zip(cps, jacobians)
        @assert isempty(missing_entries(J, cp.jacobian.pattern))
    end
    open(joinpath(OUT, "data", "jacobian_entries.csv"), "w") do io
        println(io, "case,row,column,initial_derivative")
        for (k, (cp, J)) in enumerate(zip(cps, jacobians))
            rows, cols, _ = findnz(cp.jacobian.pattern)
            for (i, j) in zip(rows, cols)
                println(io, join((k, i, j, J[i, j]), ','))
            end
        end
    end
    svg(
        "15_state_jacobian",
        "Figure 15. State Layout and Jacobian Sparsity",
        "Interleaved node and cell blocks mark differential and algebraic unknowns. Three six-cell problems show the stored sparsity and initial-state automatic-differentiation nonzeros for Lagrangian, corrected and variable-area models.";
        height = 755,
    ) do io
        for (x, w, title, body, color, fill) in (
            (35, 205, "Node 1", ["z₁ [D], w₁ [A], I₁ [A]"], BLUE, PALE),
            (260, 290, "Cell 1", ["m₁, m₂, T, G, M [D]"], TEAL, GREEN),
            (570, 205, "Node 2", ["z₂ [D], w₂ [A], I₂ [A]"], BLUE, PALE),
            (900, 265, "Global", ["A, θ_{A} [A]"], PURPLE, "#f2f2f2"),
        )
            box(io, x, 80, w, 85, title, body; color, fill)
        end
        label(io, 830, 135, "⋯"; size = 36, anchor = "middle")
        label(io, 35, 205, "D: Differential     A: Algebraic"; size = 19)
        for (k, (cp, J)) in enumerate(zip(cps, jacobians))
            x = 80 + (k-1) * 385
            y = 295
            side = 290
            S = cp.jacobian.pattern
            n = size(S, 1)
            d = side / n
            label(
                io,
                x + side/2,
                252,
                cases[k][1];
                anchor = "middle",
                size = 21,
                weight = "bold",
            )
            label(io, x + side/2, 278, "$n × $n"; anchor = "middle", size = 17)
            rect(io, x, y, side, side; fill = "white", stroke = MUTED)
            rows, cols, _ = findnz(S)
            for (i, j) in zip(rows, cols)
                color = J[i, j] == 0 ? "#d8d8d8" : BLUE
                rect(
                    io,
                    x + (j-1)*d + 0.2,
                    y + (i-1)*d + 0.2,
                    d-0.4,
                    d-0.4;
                    fill = color,
                    stroke = "none",
                )
            end
            if has_chi(cp.def.layout)
                for row in (area_index(cp.def.layout), lateral_rate_index(cp.def.layout))
                    rect(io, x, y+(row-1)*d, side, d; fill = "none", stroke = PURPLE)
                end
            end
            for i in (1, div(n, 2), n)
                label(io, x-10, y+(i-0.5)*d+5, string(i); anchor = "end", size = 15)
                label(io, x+(i-0.5)*d, y+side+24, string(i); anchor = "middle", size = 15)
            end
            label(io, x+side/2, y+side+53, "State Column"; anchor = "middle", size = 18)
            label(
                io,
                x-45,
                y+side/2,
                "Residual Row";
                anchor = "middle",
                size = 18,
                rotate = -90,
            )
        end
        rect(io, 70, 712, 17, 17; fill = "#d8d8d8", stroke = "none")
        label(io, 100, 728, "Structural Entry"; size = 18)
        rect(io, 410, 712, 17, 17; fill = BLUE, stroke = "none")
        label(io, 440, 728, "Nonzero Initial-State AD Derivative"; size = 18)
        for (x, note) in (
            (225, "Local Coupling"),
            (610, "Wider Correction Stencil"),
            (995, "Global Area Rows"),
        )
            label(io, x, 677, note; anchor = "middle", weight = "bold", size = 19)
        end

    end
end

function second_batch()
    constructor_figure()
    geometry_figure()
    shrinkage_figure()
    gas_figure()
    reactions_figure()
    conservation_figure()
    jacobian_figure()
end
