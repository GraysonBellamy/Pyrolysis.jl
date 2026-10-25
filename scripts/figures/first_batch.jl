# Editable vector compositions for figures 1–7. Details belong in the captions;
# diagrams retain only labels needed to read the geometry and physical transfers.
# Figure 4's switch is evaluated by the implementation, supplied by generate.jl.

function first_batch(activation_curve)
    svg(
        "01_mesh_anatomy",
        "Figure 1. Finite-Volume Mesh (Schematic)",
        "Nonuniform cells run from the pinned back to the exposed surface. The highlighted interior face sends equal and opposite flux contributions to its owner and neighbour; dots mark cell centres and nodes.";
        height = 420,
    ) do io
        nodes = [130, 390, 650, 865, 1050]
        label(io, 130, 93, "Pinned Back")
        label(io, 1050, 93, "Exposed Surface"; anchor = "end")
        meshstrip(io, nodes, 160; h = 95, ids = false, centers = true)
        for (i, x) in enumerate(nodes)
            label(io, x, 285, string(i); anchor = "middle", size = 18)
        end
        for i = 1:4
            label(
                io,
                (nodes[i] + nodes[i+1])/2,
                185,
                "Cell $i";
                anchor = "middle",
                size = 18,
            )
        end
        line(io, 130, 125, 65, 125; arrow = "ink")
        label(io, 65, 148, "n = −1"; size = 17)
        line(io, 1050, 125, 1120, 125; arrow = "ink")
        label(io, 1120, 148, "n = +1"; anchor = "end", size = 17)
        line(io, 1150, 200, 1060, 200; color = ORANGE, arrow = "orange", width = 2.5)
        label(io, 1108, 230, "Heat"; color = ORANGE, anchor = "middle")
        line(io, 610, 137, 705, 137; color = TEAL, arrow = "teal")
        label(io, 650, 115, "Face 3: F₃, n₃ = +1"; anchor = "middle", size = 18)
        label(io, 520, 243, "−F₃ A₃"; color = TEAL, anchor = "middle")
        label(io, 758, 243, "+F₃ A₃"; color = TEAL, anchor = "middle")
        dimension(io, 390, 650, 315, "Owner: Δz₂")
        label(io, 758, 341, "Neighbour"; anchor = "middle", size = 18)
        label(io, 130, 382, "z = 0"; anchor = "middle")
        line(io, 450, 377, 750, 377; arrow = "ink")
        label(io, 770, 383, "+z")
        label(io, 1050, 382, "z = L(t)"; anchor = "middle")
        circle(io, 130, 316, 3.5)
        label(io, 144, 322, "Node / Face"; size = 17)
        circle(io, 130, 343; fill = TEAL)
        label(io, 144, 349, "Cell Centre"; size = 17)
    end

    svg(
        "02_mesh_framework",
        "Figure 2. Mesh Framework (Schematic)",
        "Setup retains connectivity and the state layout. Each residual evaluation proceeds from evolving state and activity through geometry and derived fields, correction, kernels, closure and residual assembly.";
        height = 365,
    ) do io
        rect(io, 45, 85, 1110, 60; fill = PALE, stroke = MUTED)
        label(
            io,
            600,
            123,
            "Setup: Fixed Connectivity and State Layout";
            anchor = "middle",
            weight = "bold",
        )
        line(io, 127, 145, 127, 201; arrow = "ink")
        for (k, title, body) in (
            (1, "State", ["x, m, T/H, A", "Live-cell mask"]),
            (2, "Geometry", ["Volumes; faces", "Fields: ξ, T, P"]),
            (3, "Correction", ["c = w − vₛ"]),
            (4, "Kernels", ["Fluxes; sources"]),
            (5, "Closure", ["Volume; motion"]),
            (6, "Residual", ["M u̇ = f(u, t)"]),
        )
            x = 45 + 189(k-1)
            box(io, x, 215, 165, 108, title, body)
            k < 6 && line(io, x+169, 267, x+185, 267; arrow = "ink")
        end
    end

    svg(
        "03_mesh_motion",
        "Figure 3. Mesh Motion (Schematic)",
        "Aligned initial and evolved states distinguish mesh nodes from material markers. Fixed cells lose condensed density; Lagrangian nodes follow material; corrected ALE nodes move relative to material.";
        height = 655,
    ) do io
        label(io, 245, 87, "Initial"; weight = "bold")
        label(io, 780, 87, "Evolved"; weight = "bold")
        circle(io, 380, 81, 3.5)
        label(io, 394, 87, "Mesh Node"; size = 18)
        circle(io, 535, 81, 5; fill = "white", stroke = ORANGE)
        label(io, 549, 87, "Material"; size = 18)
        for (k, title, after, material, note) in (
            (
                1,
                "Fixed",
                [0, 0.25, 0.5, 0.75, 1],
                [0.12, 0.37, 0.62, 0.87],
                "w = 0; ρₛ decreases",
            ),
            (
                2,
                "Lagrangian",
                [0, 0.20, 0.37, 0.51, 0.62],
                [0.10, 0.285, 0.44, 0.565],
                "w = vₛ; c = 0",
            ),
            (
                3,
                "Corrected ALE",
                [0, 0.155, 0.31, 0.465, 0.62],
                [0.10, 0.285, 0.44, 0.565],
                "w = vₛ + c; vₛ − w = −c",
            ),
        )
            y = 130 + 157(k-1)
            label(io, 35, y+27, title; weight = "bold")
            meshstrip(
                io,
                245 .+ 320 .* [0, 0.25, 0.5, 0.75, 1],
                y;
                h = 48,
                ids = false,
                fills = fill(k == 1 ? GREY : PALE, 4),
            )
            meshstrip(
                io,
                780 .+ 320 .* after,
                y;
                h = 48,
                ids = false,
                fills = fill(PALE, 4),
            )
            for x in [0.12, 0.37, 0.62, 0.87]
                circle(io, 245+320x, y+75, 6; fill = "white", stroke = ORANGE)
            end
            for x in material
                circle(io, 780+320x, y+75, 6; fill = "white", stroke = ORANGE)
            end
            line(io, 610, y+24, 720, y+24; arrow = "ink")
            for x in (780, 1100)
                line(io, x, y-8, x, y+88; color = MUTED, dash = "4 4")
            end
            label(io, 245, y+117, note; color = k == 3 ? TEAL : INK)
        end
        line(io, 800, 607, 1080, 607; arrow = "ink")
        label(io, 800, 636, "Pinned Back"; size = 18)
        label(io, 1100, 636, "+z → Surface"; size = 18, anchor = "end")
    end

    svg(
        "04_equidistribution",
        "Figure 4. Layer Equidistribution (Schematic)",
        "Nonuniform reacting cells relax toward equal widths while correction is zero in the substrate and at layer boundaries. The adjacent computed activation curve rises as remaining thickness falls.";
        height = 465,
    ) do io
        label(io, 70, 100, "Inert Substrate"; size = 19)
        label(io, 325, 100, "Reacting Layer"; size = 19, color = TEAL)
        for (y, nodes) in (
            (145, [0, 95, 180, 205, 250, 325, 470]),
            (300, [0, 95, 180, 252.5, 325, 397.5, 470]),
        )
            xx = 70 .+ nodes
            meshstrip(
                io,
                xx,
                y;
                h = 55,
                ids = false,
                fills = [GREY, GREY, GREEN, GREEN, GREEN, GREEN],
            )
        end
        for (x1, x2) in ((275, 322.5), (320, 395), (395, 467.5))
            line(io, x1, 213, x2, 287; color = TEAL, arrow = "teal")
        end
        for x in (70, 250, 540)
            line(io, x, 130, x, 365; color = MUTED, dash = "4 4")
        end
        label(io, 160, 390, "c = 0"; anchor = "middle")
        label(io, 395, 390, "Equal Cell Widths"; anchor = "middle", color = TEAL)
        label(io, 305, 438, "Layer Boundaries: c = 0, w = vₛ"; anchor = "middle", size = 19)
        X, Y = chart_axes(
            io,
            770,
            145,
            360,
            210;
            xlim = (0, 0.5),
            ylim = (0, 1),
            xticks = [(0, "0"), (0.1, "0.1"), (0.3, "0.3"), (0.5, "0.5")],
            yticks = [(0, "0"), (0.5, "0.5"), (1, "1")],
            xlabel = "Remaining Thickness hᵢ / hᵢ(0)",
            ylabel = "Activation aᵢ",
        )
        r, a = activation_curve
        polyline(io, X.(r), Y.(a); color = TEAL)
    end

    svg(
        "05_deactivation",
        "Figure 5. Burnout and Retained State (Schematic)",
        "A surface-connected band transfers inventory to cell 2 on a moving mesh; cells 3 and 4 become inactive but their state slots remain. Small insets distinguish fixed-mesh gas venting and an interior void blocked by live overlying material.";
        height = 610,
    ) do io
        label(io, 65, 100, "Moving Mesh: Before"; weight = "bold")
        label(io, 650, 100, "After Transfer"; weight = "bold")
        meshstrip(
            io,
            [65, 180, 295, 410, 550],
            190;
            h = 65,
            fills = [PALE, GREEN, GREY, GREY],
        )
        line(io, 510, 165, 235, 165; color = TEAL, arrow = "teal", width = 2.5)
        label(
            io,
            375,
            147,
            "Mass, Enthalpy, Volume";
            anchor = "middle",
            color = TEAL,
            size = 19,
        )
        dimension(io, 295, 550, 276, "Selected Band: rᵢ < ε_{band}")
        label(io, 65, 344, "Trigger: r₄ ≤ εₘ"; size = 19)
        line(io, 570, 222, 630, 222; arrow = "ink")
        meshstrip(io, [650, 765, 1135], 190; h = 65, fills = [PALE, GREEN])
        line(io, 1135, 178, 1135, 267; color = ORANGE, width = 3)
        label(io, 1125, 155, "Active Surface"; anchor = "end", color = ORANGE, size = 19)
        label(io, 650, 302, "Retained Slots"; size = 18)
        for (i, active) in enumerate((1, 1, 0, 0))
            x = 825 + 77(i-1)
            label(io, x+30, 288, string(i); anchor = "middle", size = 17)
            rect(io, x, 300, 60, 40; fill = active == 1 ? GREEN : GREY, stroke = MUTED)
            label(io, x+30, 327, string(active); anchor = "middle")
        end
        label(io, 650, 365, "Mask: 1 = Live; 0 = Inactive"; size = 18)
        panel(io, 35, 405, 530, 180, "Fixed Mesh: Gas Vents")
        meshstrip(
            io,
            [65, 170, 275, 380, 485],
            475;
            h = 45,
            fills = [PALE, GREEN, GREY, GREY],
        )
        line(io, 440, 467, 440, 440; arrow = "teal", color = TEAL)
        label(io, 490, 458, "Gas"; color = TEAL, size = 18)
        line(io, 430, 549, 222, 549; arrow = "teal", color = TEAL)
        label(io, 325, 581, "Condensed Inventory"; anchor = "middle", size = 18)
        panel(io, 615, 405, 550, 180, "Interior Void: Stop")
        meshstrip(
            io,
            [650, 765, 880, 995, 1110],
            475;
            h = 45,
            fills = [PALE, GREY, PALE, PALE],
        )
        label(io, 820, 555, "r₂ ≤ εₘ"; anchor = "middle", size = 18)
        dimension(io, 880, 1110, 540, "Live Barrier")
    end

    svg(
        "06_surface_energy",
        "Figure 6. Surface Energy Transfers (Schematic)",
        "Absorbed radiation balances inward conduction and outward radiative and convective losses. Gas-carried enthalpy is shown separately as a contribution to the outward face energy flux, not to the surface-temperature balance.";
        height = 555,
    ) do io
        rect(io, 90, 140, 510, 200; fill = PALE)
        label(io, 115, 172, "Surface Cell")
        circle(io, 250, 265, 6; fill = TEAL)
        label(io, 250, 298, "T_{cell}"; anchor = "middle", color = TEAL)
        line(io, 600, 125, 600, 350; color = ORANGE, width = 3)
        label(io, 600, 108, "Tₛ"; color = ORANGE, anchor = "middle", size = 25)
        line(io, 850, 165, 612, 165; color = ORANGE, arrow = "orange", width = 2.5)
        label(io, 870, 172, "Absorbed Radiation"; color = ORANGE)
        line(io, 612, 235, 850, 235; color = ORANGE, arrow = "orange", width = 2.5)
        label(io, 870, 242, "Re-radiation"; color = ORANGE)
        line(io, 612, 305, 850, 305; color = BLUE, arrow = "blue", width = 2.5)
        label(io, 870, 312, "Convection"; color = BLUE)
        line(io, 575, 220, 280, 220; color = ORANGE, arrow = "orange", width = 2.5)
        label(io, 415, 200, "Conduction q_{cond}"; anchor = "middle", color = ORANGE)
        dimension(io, 250, 600, 370, "d")
        line(io, 50, 425, 1150, 425; color = MUTED, width = 0.8)
        label(io, 80, 470, "Outward Face Flux"; weight = "bold")
        label(io, 80, 510, "F_{E} = −q_{cond} + Σⱼ Jⱼ hⱼ(Tₛ)"; size = 23)
        line(io, 610, 475, 850, 475; arrow = "teal", color = TEAL, width = 2.5)
        label(io, 870, 482, "Gas Enthalpy"; color = TEAL)
    end

    svg(
        "07_radiation",
        "Figure 7. Radiation Absorption (Schematic)",
        "Matching slabs receive radiation from the right: surface absorption deposits it at the boundary, while Beer–Lambert attenuation distributes absorption through cells and transmits the remainder. The plotted intensity uses uniform optical depth alpha L equal to two.";
        height = 480,
    ) do io
        label(io, 65, 95, "Surface Absorption"; weight = "bold")
        label(io, 65, 260, "Beer–Lambert Absorption"; weight = "bold")
        for (y, shades) in
            ((140, fill(PALE, 4)), (310, ["#f5f5f5", "#ebebeb", "#dddddd", "#cccccc"]))
            meshstrip(io, [90, 190, 290, 390, 490], y; h = 60, ids = false, fills = shades)
            line(io, 600, y+30, 502, y+30; arrow = "orange", color = ORANGE, width = 2.5)
        end
        line(io, 490, 130, 490, 210; color = ORANGE, width = 4)
        label(io, 550, 135, "a q_{inc}"; anchor = "middle", color = ORANGE)
        for (j, x) in enumerate([490, 390, 290, 190, 90])
            label(io, x, 296, "I_{$(j-1)}"; anchor = "middle", color = ORANGE, size = 18)
            j < 5 && line(
                io,
                x-8,
                340,
                x-88,
                340;
                arrow = "orange",
                color = ORANGE,
                width = 3-0.5j,
            )
        end
        line(io, 80, 340, 40, 340; arrow = "orange", color = ORANGE)
        label(io, 65, 403, "Transmitted"; size = 18, color = ORANGE)
        label(
            io,
            315,
            403,
            "Absorption per Area: I_{in} − I_{out}";
            anchor = "middle",
            size = 18,
        )
        line(io, 200, 443, 470, 443; arrow = "ink")
        label(io, 185, 450, "z = 0"; anchor = "end", size = 18)
        label(io, 490, 450, "z = L"; size = 18)
        X, Y = chart_axes(
            io,
            810,
            145,
            320,
            225;
            xlim = (0, 1),
            ylim = (0, 1),
            xticks = [(0, "0"), (0.5, "0.5"), (1, "1")],
            yticks = [(0, "0"), (0.5, "0.5"), (1, "1")],
            xlabel = "Inward Depth (L − z)/L",
            ylabel = "I / I₀",
        )
        r = range(0, 1; length = 101)
        polyline(io, X.(r), Y.(exp.(-2 .* r)); color = ORANGE)
        for r = 0:0.25:1
            circle(io, X(r), Y(exp(-2r)), 4; fill = ORANGE)
        end
        label(io, 970, 110, "αL = 2"; anchor = "middle", size = 21)
    end
end
