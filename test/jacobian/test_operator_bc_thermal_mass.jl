"""
`BCThermalOperator` + `BCMassOperator` correctness tests.

Three gates:

1. **`apply!` equivalence.** The operators' boundary-cell `(dT, dξ)`
   contributions match the residual's BC slice, reconstructed from
   `apply_thermal_bc_flux!` + `apply_mass_bc_flux!` followed by
   `divergence_1d!` / `divergence_species_1d!` at the boundary cells.
   The boundary cell's `S_conv` energy row is NOT part of this slice —
   it is owned whole by `BCConvectionOperator` (tested in
   `test_operator_bc_convection.jl`); `BCMassOperator`'s energy row
   carries only the `ΔS_gen` boundary-face slice, which is zero under
   the `:standard` energy form used here.

2. **`contribute_to_jacobian!` equivalence.** Per-boundary AD agrees with
   `ForwardDiff.jacobian` over the same reference closure on the operators'
   declared support; off-pattern entries are zero
   (support-completeness on the BC slice).

3. **`∂T_s/∂ξ` regression.** ConvectiveBC + ConvectiveMassBC at
   `hot+gas` state — the regime where the legacy `bcs.jl:343-353` no-blowing
   chain rule misses `∂T_s/∂ξ` via `k_eff(ξ)`. The AD-through-Newton path
   recovers full agreement at rtol=1e-6. The baseline rel error
   (2.1e-5, the gap measured before the fix) drops to ≤ rtol.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: n_components, Material
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Discretization:
    apply_thermal_bc_flux!,
    apply_mass_bc_flux!,
    divergence_1d!,
    divergence_species_1d!,
    compute_conductive_fluxes!,
    compute_radiation_source!
using Pyrolysis.Physics: compute_gas_fluxes!
using Pyrolysis.Problem: populate!, build_ad_workspace
using Pyrolysis.Jacobian:
    BCThermalOperator,
    BCMassOperator,
    build_bc_thermal_workspace,
    build_bc_mass_workspace,
    StructuredJacobian,
    build_residual_context

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Build a Eulerian convective-BC + convective-mass-BC problem (the
# convective-BC regime without the ALE / radiation / Darcy axes, so the
# BC operators can be tested in isolation).
function _build_eulerian_convective_problem()
    mat = let
        virgin = SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, λ = 1e-5)
        gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        Material(name = :wood_conv_eul, components = (virgin, gas), reactions = (rxn,))
    end
    mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(2))
    bc = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(50e3) + ConvectiveBC(h = 25.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0), bottom = ImpermeableBC()),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    def, ws, u = initialize_state(problem; use_ale = false)
    return def, ws, u
end

function _set_uniform_state!(u, n_cells, T, ξ_vals)
    for i = 1:n_cells
        u[i] = T
    end
    for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Reference BC contribution: matches the BC slice of compute_rhs!.
# Returns dT_bc, dξ_bc (NC × n_cells), populated only at the boundary cells.
function _reference_bc_contribution(
    state,
    mesh,
    material::Material{NC},
    prop_cache,
    J_gas_buf,
    q_buf,
    div_q_buf,
    div_J_buf,
    bc_set,
    top_face_ids,
    bottom_face_ids,
    t,
    radiation_model = Pyrolysis.Materials.NO_RADIATION,
) where {NC}

    update_properties!(prop_cache, state, mesh, material)

    # Match `compute_rhs!` ordering: populate J_gas (interior + BC) BEFORE
    # the thermal BC solve so the surface transpiration term sees the actual
    # boundary-face mass flux.
    compute_gas_fluxes!(J_gas_buf, state, mesh, prop_cache, material)
    apply_mass_bc_flux!(
        J_gas_buf,
        state,
        mesh,
        prop_cache,
        bc_set,
        material,
        t,
        top_face_ids,
        bottom_face_ids,
    )

    # Conductive heat flux at interior faces (zero at boundary).
    compute_conductive_fluxes!(q_buf, state, mesh, prop_cache)
    # Apply thermal BCs — overwrites q_buf at boundary faces. Transpiration
    # is OFF, matching production: `to_problem_def` rejects
    # `use_transpiration_bc = true`, so the residual's surface Newton always
    # sees ṁ_g = 0 and the Jacobian must differentiate that same balance.
    apply_thermal_bc_flux!(
        q_buf,
        state,
        mesh,
        bc_set,
        prop_cache,
        material,
        J_gas_buf,
        t,
        top_face_ids,
        bottom_face_ids,
        radiation_model;
        use_transpiration_bc = false,
    )
    # Isolate BC's contribution to dT via divergence of (zero interior, BC at
    # boundary).
    q_bc_only = zeros(eltype(q_buf), length(q_buf))
    @inbounds for f in top_face_ids
        q_bc_only[f] = q_buf[f]
    end
    @inbounds for f in bottom_face_ids
        q_bc_only[f] = q_buf[f]
    end

    divergence_1d!(div_q_buf, q_bc_only, state, mesh)
    # BC-only J: zero at interior, BC value at boundary.
    J_bc_only = zeros(eltype(J_gas_buf), size(J_gas_buf)...)
    @inbounds for f in top_face_ids
        for j = 1:NC
            J_bc_only[j, f] = J_gas_buf[j, f]
        end
    end
    @inbounds for f in bottom_face_ids
        for j = 1:NC
            J_bc_only[j, f] = J_gas_buf[j, f]
        end
    end
    divergence_species_1d!(div_J_buf, J_bc_only, state, mesh)

    n_cells = length(state.T)
    Tu = eltype(state.T)

    # No S_conv piece here: the boundary cell's S_conv energy row (whose
    # face-average includes the BC flux) is owned whole by
    # BCConvectionOperator, and the ΔS_gen slice BCMassOperator does own is
    # zero under the :standard energy form this fixture uses.
    dT_bc = zeros(Tu, n_cells)
    dξ_bc = zeros(Tu, NC, n_cells)
    @inbounds for i = 1:n_cells
        ρcp = prop_cache.ρcp[i]
        dT_bc[i] = -div_q_buf[i] / ρcp
        for j = 1:NC
            if Pyrolysis.Materials.is_gas(material.components[j])
                dξ_bc[j, i] = -div_J_buf[j, i]
            end
        end
    end
    return dT_bc, dξ_bc
end

function _reference_bc_jacobian(
    u,
    ws,
    def,
    ::Material{NC},
    bc_set,
    top_face_ids,
    bottom_face_ids,
    radiation_model = Pyrolysis.Materials.NO_RADIATION,
) where {NC}
    layout = ws.layout
    mesh = ws.mesh
    material = def.material
    n_cells = layout.n_cells
    n_total = length(u)
    out_dim = n_cells * (1 + NC)
    x0 = copy(u)
    f_buf = zeros(out_dim)
    cfg = ForwardDiff.JacobianConfig(
        nothing,
        f_buf,
        x0,
        ForwardDiff.Chunk{min(12, n_total)}(),
    )
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)

    J_gas_buf = zeros(Tu, NC, length(mesh.faces))
    q_buf = zeros(Tu, length(mesh.faces))
    div_q_buf = zeros(Tu, n_cells)
    div_J_buf = zeros(Tu, NC, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = mesh,
            material = material,
            def = def,
            J_gas_buf = J_gas_buf,
            q_buf = q_buf,
            div_q_buf = div_q_buf,
            div_J_buf = div_J_buf,
            bc_set = bc_set,
            top_face_ids = top_face_ids,
            bottom_face_ids = bottom_face_ids,
            radiation_model = radiation_model

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.Problem.update_state_geometry!(
                    ws_ad.state_cache,
                    mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )
                dT, dξ = _reference_bc_contribution(
                    ws_ad.state_cache,
                    mesh,
                    material,
                    ws_ad.prop_cache,
                    J_gas_buf,
                    q_buf,
                    div_q_buf,
                    div_J_buf,
                    bc_set,
                    top_face_ids,
                    bottom_face_ids,
                    0.0,
                    radiation_model,
                )
                @inbounds for i = 1:n_cells
                    out[i] = dT[i]
                end
                @inbounds for j = 1:NC
                    off = n_cells + (j - 1) * n_cells
                    for i = 1:n_cells
                        out[off+i] = dξ[j, i]
                    end
                end
                return nothing
            end
        end
    J = ForwardDiff.jacobian(f!, f_buf, x0, cfg)
    return J
end

@testset "BCThermalOperator + BCMassOperator" begin

    @testset "interface declarations" begin
        thermal = BCThermalOperator()
        mass = BCMassOperator()
        @test Pyrolysis.Jacobian.reads(thermal) == ()
        @test Pyrolysis.Jacobian.writes(thermal) == (:T_s,)
        @test Pyrolysis.Jacobian.stencil(thermal) == (cells = 0:0,)
        @test Pyrolysis.Jacobian.reads(mass) == ()
        @test Pyrolysis.Jacobian.writes(mass) == ()
        @test Pyrolysis.Jacobian.stencil(mass) == (cells = 0:0,)
    end

    def, ws, u = _build_eulerian_convective_problem()
    n = ws.layout.n_cells
    NC = n_components(ws.layout)
    material = def.material
    bc_set = def.bcs
    top_face_ids = ws.top_face_ids
    bottom_face_ids = ws.bottom_face_ids

    # Hot + gas-heavy state — convective-BC regime.
    _set_uniform_state!(u, n, 1000.0, (400.0, 5.0))

    # Reference scratch
    prop_cache_ref = PropertyCache(n, length(ws.mesh.faces))
    J_gas_ref = zeros(NC, length(ws.mesh.faces))
    q_ref = zeros(length(ws.mesh.faces))
    div_q_ref = zeros(n)
    div_J_ref = zeros(NC, n)

    @testset "apply! matches BC residual slice — hot+gas state" begin
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        dT_ref, dξ_ref = _reference_bc_contribution(
            ws.state_cache,
            ws.mesh,
            material,
            prop_cache_ref,
            J_gas_ref,
            q_ref,
            div_q_ref,
            div_J_ref,
            bc_set,
            top_face_ids,
            bottom_face_ids,
            0.0,
        )

        dT_op = zeros(n)
        dξ_op = zeros(NC, n)
        Pyrolysis.Jacobian.apply!(
            dT_op,
            BCThermalOperator(),
            ws.state_cache,
            ws.mesh,
            material,
            bc_set,
            0.0,
            Pyrolysis.Materials.NO_RADIATION,
        )
        Pyrolysis.Jacobian.apply!(
            dT_op,
            dξ_op,
            BCMassOperator(),
            ws.state_cache,
            ws.mesh,
            material,
            bc_set,
            0.0,
            false,  # Fickian fixture — no Darcy venting term
        )

        # Boundary cells are 1 (bottom) and n (top) — only those have nonzero
        # BC contribution. Compare directly.
        @test isapprox(dT_op, dT_ref; rtol = 1e-9, atol = 1e-9)
        @test isapprox(dξ_op, dξ_ref; rtol = 1e-9, atol = 1e-9)
    end

    @testset "contribute_to_jacobian! matches DenseAD on declared support" begin
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        # Build sparse pattern: 1×(1+NC) thermal/energy block at T-row[bdy] ×
        # {T-col[bdy], ξ-cols[bdy]} (both operators are cell-local now —
        # the neighbor-reaching S_conv row belongs to BCConvectionOperator),
        # plus NC×(1+NC) mass block at ξ-rows[bdy] × {T-col[bdy], ξ-cols[bdy]}.
        # Boundary cells are 1 and n.
        n_total = length(u)
        I = Int[]
        J_idx = Int[]
        V = Float64[]
        for cell_id in (1, n)
            # Thermal + BC-mass energy row: T-row[cell_id]
            push!(I, cell_id)
            push!(J_idx, cell_id)
            push!(V, 1.0)
            for j = 1:NC
                ξ_off = n + (j - 1) * n
                push!(I, cell_id)
                push!(J_idx, ξ_off + cell_id)
                push!(V, 1.0)
            end
            # Mass: ξ-rows[j, cell_id]
            for r_j = 1:NC
                gr = n + (r_j - 1) * n + cell_id
                push!(I, gr)
                push!(J_idx, cell_id)
                push!(V, 1.0)
                for j = 1:NC
                    ξ_off = n + (j - 1) * n
                    push!(I, gr)
                    push!(J_idx, ξ_off + cell_id)
                    push!(V, 1.0)
                end
            end
        end
        A = sparse(I, J_idx, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(:eulerian, A)

        ctx = build_residual_context([BCThermalOperator(), BCMassOperator()])
        ws_thermal = build_bc_thermal_workspace(material)
        ws_mass = build_bc_mass_workspace(material)

        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            BCThermalOperator(),
            u,
            ctx,
            0.0,
            ws_thermal,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
            bc_set,
            Pyrolysis.Materials.NO_RADIATION,
        )
        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            BCMassOperator(),
            u,
            ctx,
            0.0,
            ws_mass,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
            bc_set,
            false,  # Fickian fixture — no Darcy venting term
        )

        J_ref = _reference_bc_jacobian(
            u,
            ws,
            def,
            material,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        J_op = Matrix(jac.A)
        max_dev = 0.0
        for (r, c) in zip(I, J_idx)
            diff = abs(J_op[r, c] - J_ref[r, c])
            max_dev = max(max_dev, diff)
        end
        @test max_dev < 1e-6

        # Support-completeness on the (T, ξ) physics block — every reference
        # entry outside the declared BC support must be zero (the BC operator
        # only writes at boundary cells).
        n_phys = n * (1 + NC)
        Pmask = falses(n_phys, n_phys)
        for (r, c) in zip(I, J_idx)
            if r <= n_phys && c <= n_phys
                Pmask[r, c] = true
            end
        end
        @test maximum(abs.(J_ref[1:n_phys, 1:n_phys][.!Pmask])) < 1e-6
    end

    @testset "∂T_s/∂ξ via k_eff(ξ) recovered" begin
        # At hot+gas state, the legacy `_compute_boundary_flux_jacobian` uses
        # `zero_J` and misses the chain rule through k_eff(ξ). The baseline
        # records `max_rel_dev = 2.116e-5` on the analytical path.
        # The AD-through-Newton operator must agree with DenseAD at rtol=1e-6
        # (i.e., the baseline number drops by orders of magnitude).
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        # Full BC contribution to T-row[top]: the thermal surface Newton
        # (cell-local; the S_conv row with its neighbor reach belongs to
        # BCConvectionOperator).
        n_total = length(u)
        I = Int[]
        J_idx = Int[]
        V = Float64[]
        for cell_id in (1, n)
            push!(I, cell_id)
            push!(J_idx, cell_id)
            push!(V, 1.0)
            for j = 1:NC
                ξ_off = n + (j - 1) * n
                push!(I, cell_id)
                push!(J_idx, ξ_off + cell_id)
                push!(V, 1.0)
                # mass rows (needed because BCMassOperator scatters them too)
                gr = n + (j - 1) * n + cell_id
                push!(I, gr)
                push!(J_idx, cell_id)
                push!(V, 1.0)
                for k = 1:NC
                    ξk_off = n + (k - 1) * n
                    push!(I, gr)
                    push!(J_idx, ξk_off + cell_id)
                    push!(V, 1.0)
                end
            end
        end
        A = sparse(I, J_idx, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(:eulerian, A)
        ctx = build_residual_context([BCThermalOperator(), BCMassOperator()])
        ws_thermal = build_bc_thermal_workspace(material)
        ws_mass = build_bc_mass_workspace(material)

        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            BCThermalOperator(),
            u,
            ctx,
            0.0,
            ws_thermal,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
            bc_set,
            Pyrolysis.Materials.NO_RADIATION,
        )
        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            BCMassOperator(),
            u,
            ctx,
            0.0,
            ws_mass,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
            bc_set,
            false,
        )

        # Reference: ForwardDiff over the full BC contribution to dT.
        J_ref_full = _reference_bc_jacobian(
            u,
            ws,
            def,
            material,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        # Entries of interest at the top cell's T-row: the T column, and the
        # ξ-solid column — the ∂T_s/∂ξ chain through k_eff(ξ) that the legacy
        # analytical path missed. (The ξ-gas column is exactly zero here:
        # conductivity mixes on the renormalized solid basis, and with
        # transpiration hard-zeroed — matching the residual — no other
        # thermal-BC path reaches ξ_gas at a uniform-T state.)
        top = n
        ξ_solid_off = n + (1 - 1) * n   # solid is component 1
        ξ_gas_off = n + (2 - 1) * n     # gas is component 2

        thermal_T_dev = abs(jac.A[top, top] - J_ref_full[top, top])
        thermal_ξT_dev = abs(jac.A[top, ξ_solid_off+top] - J_ref_full[top, ξ_solid_off+top])
        thermal_ξgT_dev = abs(jac.A[top, ξ_gas_off+top] - J_ref_full[top, ξ_gas_off+top])

        # rtol=1e-6 against the reference scale.
        ref_scale_T = max(abs(J_ref_full[top, top]), 1.0)
        ref_scale_ξT = max(abs(J_ref_full[top, ξ_solid_off+top]), 1.0)
        @test thermal_T_dev / ref_scale_T < 1e-6
        @test thermal_ξT_dev / ref_scale_ξT < 1e-6
        @test thermal_ξgT_dev / ref_scale_T < 1e-6

        # Sanity: the ∂T_s/∂ξ_solid entry is genuinely non-zero (otherwise the
        # regression test is meaningless — it would pass even on the legacy
        # zero-J path).
        @test abs(J_ref_full[top, ξ_solid_off+top]) > 1e-6
    end
end
