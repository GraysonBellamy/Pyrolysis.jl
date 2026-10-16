"""
The `RadiationOperator` — covers the no-radiation / surface-absorption /
P1 modes plus the Beer-Lambert volumetric path.

Represents the volumetric radiation source's contribution to `dT`:

```
dT[i] += Q_rad(state, mesh, model) / ρcp_eff(T_i, ξ_i)
```

# Mode coverage

| Model               | apply!                                    | Jacobian contribution        |
|---------------------|-------------------------------------------|------------------------------|
| `NO_RADIATION`      | no-op                                     | no-op                        |
| `SURFACE_ABSORPTION`| no-op (handled by BC)                     | no-op                        |
| `P1_QUASI_STEADY`   | no-op (handled in radiation cache)        | no-op                        |
| `BEER_LAMBERT`      | wraps `compute_beer_lambert_source!`      | `EmbeddedScaledSuffixSumCoupling` per channel |

# Beer-Lambert math

`q_rad[i] = I_in[i+1] · (1 − exp(−α_i Δz_i)) / Δz_i` where
`I_in[i+1] = α_surf · I_surface · exp(−Σ_{j > i} α_j Δz_j)`. On the
Kirchhoff default path (`apply_kirchhoff_absorptivity = true`) the surface
transmissivity is the live surface cell's effective surface absorptivity,
`α_surf = a_eff(T_n, ξ_·,n)` (`effective_absorptivity`; identical to
`ε_eff(T_n, ξ_·,n)` unless a component sets `surface_absorptivity`); with an
explicit BC absorptivity it is already folded into `I_surface` and
`α_surf = 1`. Chain rule for `j > i`:

```
∂q_rad[i] / ∂x_j = −q_rad[i] · Δz_j · ∂α_j/∂x_j      (i < j < n)
∂q_rad[i] / ∂x_n = −q_rad[i] · (Δz_n · ∂α_n/∂x_n − (∂a_eff/∂x_n)/a_eff)
∂q_rad[i] / ∂x_i = (I_in[i+1] − q_rad[i] · Δz_i) · ∂α_i/∂x_i  (j == i < n)
∂q_rad[i] / ∂x_j = 0                                  (j < i)
```

(the `∂a_eff` terms appear only on the Kirchhoff path; the surface cell's
own diagonal additionally carries `∂a_eff/∂x_n · I_surface(1−e^{−τ_n})/Δz_n`,
captured by the per-cell AD closure below).

The `j > i` block factors as `Diag(−q_rad) · U · diag(v)` where
`v[j] = Δz_j · ∂α_j/∂x_j − δ_{j,n}·(∂ε_eff/∂x_n)/ε_eff` and `U` is the
strict upper-triangular all-ones matrix. It is pushed as an
[`EmbeddedScaledSuffixSumCoupling`](@ref) with row scale `−q_rad[i]`; the
diagonal `j == i` term is added separately.

# Known approximation under ALE — missing z-channel

This operator fills **T and ξ columns only**. Under ALE the residual's
Beer–Lambert deposit also depends on the live geometry (`Δz_i = V_i/A` from
the z block of the state): `∂q_rad[i]/∂z_node ≠ 0` through both the per-cell
`(1 − e^{−α_iΔz_i})/Δz_i` factor and the transmission suffix
`exp(−Σ_{j>i} α_jΔz_j)`. Those z-columns are **deliberately omitted** — the
Jacobian is exact in (T, ξ) and inexact in z whenever
`radiation_model == BEER_LAMBERT` and the mesh moves. This degrades only the
implicit integrator's Newton convergence (sensitivities do not consume this
Jacobian); the residual itself is exact. The omission cannot be expressed
through the backend-level `is_approximate` trait (it is mode- and
layout-dependent, not backend-dependent) — verifier runs that perturb the z
block under BEER_LAMBERT + ALE must expect mismatches in these slots. A full
z-channel needs geometry-suffix couplings analogous to the α suffix above.
"""

using ForwardDiff
using ..Materials:
    Material,
    RadiationModel,
    NO_RADIATION,
    SURFACE_ABSORPTION,
    BEER_LAMBERT,
    P1_QUASI_STEADY,
    n_components
using ..Properties:
    PropertyCache, effective_absorption, effective_heat_capacity, effective_absorptivity
using ..Physics: compute_beer_lambert_source!

"""
	RadiationOperator <: Operator

Stateless tag. No-op for `NO_RADIATION` / `SURFACE_ABSORPTION` /
`P1_QUASI_STEADY` modes; for `BEER_LAMBERT` it produces a non-local
volumetric source whose Jacobian factors as
`Diag(−q_rad) · U · diag(v)` plus a diagonal contribution.
"""
struct RadiationOperator <: Operator end

"""
	BeerLambertWorkspace

Holds the Beer-Lambert workspace: the per-cell `q_rad`, the cell-edge
incident intensity `I_in[i+1]`, and the per-cell `Δz · ∂α/∂x` channel
buffers consumed when building the α suffix-sum coupling.

# Fields
- `q_rad::Vector{Float64}`        — per-cell volumetric source [W/m³].
- `I_in::Vector{Float64}`         — incident intensity entering cell
                                    `i` from above. Length `n_cells`,
                                    with `I_in[N] = α_surf · I_surface`
                                    (`α_surf = ε_eff` on the Kirchhoff
                                    path, 1 otherwise).
- `dα_dT::Vector{Float64}`        — `∂α_i/∂T_i`, per cell.
- `dα_dξ::Matrix{Float64}`        — `∂α_i/∂ξ_{j,i}`, `(NC, n_cells)`.
- `x_buf`, `f_buf`, `J_local`, `cfg` — preallocated ForwardDiff scratch
  for the diagonal `(T_i, ξ_·,i)` AD pass in `contribute_to_jacobian!`.
- `row_scale::Vector{Float64}`    — `-q_rad[i]/ρcp_i` per cell; aliased
  by every suffix coupling pushed in one assembly.
- `v_channels::Vector{Vector{Float64}}` — per-channel `Δz·∂α/∂x` payloads
  (entry 1 = T, entry 1+k = ξ_k), refilled in place each assembly and
  aliased by the pushed couplings.
"""
mutable struct BeerLambertWorkspace{CFG}
    q_rad::Vector{Float64}
    I_in::Vector{Float64}
    dα_dT::Vector{Float64}
    dα_dξ::Matrix{Float64}

    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    row_scale::Vector{Float64}
    v_channels::Vector{Vector{Float64}}
end

"""
	build_beer_lambert_workspace(material, n_cells) -> BeerLambertWorkspace

Allocate the Beer-Lambert workspace once at problem-build
time. `α` derivatives default to zero — populate them via
[`compute_alpha_derivatives!`](@ref) when actually needed.
"""
function build_beer_lambert_workspace(
    material::Material{NC,NR},
    n_cells::Integer,
) where {NC,NR}
    dim = 1 + NC
    x_buf = zeros(Float64, dim)
    f_buf = zeros(Float64, 1)
    J_local = zeros(Float64, 1, dim)
    cfg = ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{dim}())
    return BeerLambertWorkspace{typeof(cfg)}(
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, NC, n_cells),
        x_buf,
        f_buf,
        J_local,
        cfg,
        zeros(Float64, n_cells),
        [zeros(Float64, n_cells) for _ = 1:dim],
    )
end

"""
	apply!(dT, op::RadiationOperator, state, mesh, prop_cache, material,
	       I_radiation, model, apply_kirchhoff_absorptivity) -> nothing

For `NO_RADIATION` / `SURFACE_ABSORPTION` / `P1_QUASI_STEADY`, the
volumetric source is zero — `dT` is unchanged. For `BEER_LAMBERT`,
calls the existing [`compute_beer_lambert_source!`](@ref) into `dT`
(treating `dT` as the `q_rad` output buffer; the caller is responsible
for the `/ ρcp_eff` divide and final sum into the residual).
"""
function apply!(
    dT::AbstractVector,
    ::RadiationOperator,
    state,
    mesh,
    prop_cache,
    material::Material{NC},
    I_radiation::Real,
    model::RadiationModel,
    apply_kirchhoff_absorptivity::Bool,
) where {NC}
    if model == NO_RADIATION || model == SURFACE_ABSORPTION || model == P1_QUASI_STEADY
        return nothing
    end
    if model == BEER_LAMBERT
        # Wrap the existing physics function — no math duplication.
        compute_beer_lambert_source!(
            dT,
            state,
            mesh,
            prop_cache,
            I_radiation,
            apply_kirchhoff_absorptivity,
        )
        return nothing
    end
    error("RadiationOperator: unknown radiation model $(model).")
end

"""
	apply!(ws::BeerLambertWorkspace, op::RadiationOperator,
	       state, mesh, prop_cache, material, I_radiation,
	       model::RadiationModel, apply_kirchhoff_absorptivity) -> nothing

BEER_LAMBERT-specific apply! that also populates the cell-edge
`I_in[i+1]` array on the workspace, in addition to filling `ws.q_rad`.
Both are needed to compose the structured-coupling chain rule.
`I_in[n]` is the flux that actually crosses the surface —
`ε_eff·I_radiation` on the Kirchhoff path, `I_radiation` otherwise.

For non-Beer-Lambert modes this is a no-op (q_rad, I_in are left
zero). Errors on unknown models.
"""
function apply!(
    ws::BeerLambertWorkspace,
    op::RadiationOperator,
    state,
    mesh,
    prop_cache,
    material::Material{NC},
    I_radiation::Real,
    model::RadiationModel,
    apply_kirchhoff_absorptivity::Bool,
) where {NC}
    if model == NO_RADIATION || model == SURFACE_ABSORPTION || model == P1_QUASI_STEADY
        fill!(ws.q_rad, 0.0)
        fill!(ws.I_in, 0.0)
        return nothing
    end
    if model == BEER_LAMBERT
        # 1) Fill q_rad via the shared physics function.
        compute_beer_lambert_source!(
            ws.q_rad,
            state,
            mesh,
            prop_cache,
            I_radiation,
            apply_kirchhoff_absorptivity,
        )

        # 2) Re-derive I_in[i+1] from the same monotone integration. Cheap
        #    (same loop structure as compute_beer_lambert_source!) and
        #    cleanly avoids leaking workspace into the physics function.
        n = mesh.n_active_cells
        fill!(ws.I_in, 0.0)
        A_section = state.A_section[]
        I_current = Float64(I_radiation)
        if apply_kirchhoff_absorptivity && n > 0
            I_current *= Float64(prop_cache.a[n])
        end
        @inbounds for i = n:-1:1
            ws.I_in[i] = I_current
            α = prop_cache.α[i]
            Δz = state.volume[i] / A_section
            if α > 1e-20
                τ = α * Δz
                I_current *= exp(-τ)
            end
        end
        return nothing
    end
    error("RadiationOperator: unknown radiation model $(model).")
end

"""
	compute_alpha_derivatives!(ws::BeerLambertWorkspace, state, mesh, prop_cache,
	                           material) -> nothing

Populate `ws.dα_dT` and `ws.dα_dξ` via cell-local AD over the
property's effective absorption coefficient `α(T_i, ξ_i)` for the live
cells; ghost entries are zero. Called once per Jacobian assembly, before
the α suffix-sum coupling is built.

Implementation: per-cell stack-Dual sweep over the property function
(no `ForwardDiff.gradient` machinery — the per-call closure + config it
builds dominated the Beer-Lambert Jacobian allocation profile). For
materials whose `α` is state-independent (constant per component,
weighted by ξ), the derivatives reduce to constants and this is cheap.
"""

# Tag for the Beer-Lambert property stack-Dual sweeps (α and the
# Kirchhoff ε_eff chain). Keeps these duals from colliding with any
# outer AD in flight.
struct _BLPropTag end

function compute_alpha_derivatives!(
    ws::BeerLambertWorkspace,
    state,
    mesh,
    prop_cache,
    material::Material{NC,NR},
) where {NC,NR}
    # The standard linear-mixing form `α = Σⱼ ξⱼ · αⱼ` makes ∂α/∂T = 0 and
    # ∂α/∂ξⱼ = αⱼ; the Dual sweep keeps this generic for any future
    # ξ-dependent absorption law.
    n_cells = mesh.n_active_cells
    @inbounds for i = (n_cells+1):length(ws.dα_dT)
        ws.dα_dT[i] = 0.0
        for j = 1:NC
            ws.dα_dξ[j, i] = 0.0
        end
    end
    @inbounds for i = 1:n_cells
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        _, ξ_d = _stack_dual_seed(_BLPropTag, state.T[i], ξ_i)
        α_d = effective_absorption(ξ_d, material)
        # `effective_absorption` takes no T argument, so ∂α/∂T is
        # structurally zero (the seed's slot 1 is reserved for T).
        ws.dα_dT[i] = 0.0
        if α_d isa ForwardDiff.Dual
            p = ForwardDiff.partials(α_d)
            for j = 1:NC
                ws.dα_dξ[j, i] = p[1+j]
            end
        else
            for j = 1:NC
                ws.dα_dξ[j, i] = 0.0
            end
        end
    end
    return nothing
end

"""
	contribute_to_jacobian!(jac, op::RadiationOperator,
	                        u, t, ws, model) -> nothing

No-op for `NO_RADIATION` / `SURFACE_ABSORPTION` / `P1_QUASI_STEADY`. For
`BEER_LAMBERT`, callers must dispatch to the workspace-carrying overload
so the diagonal absorption and suffix-sum terms are assembled.
"""
function contribute_to_jacobian!(
    ::StructuredJacobian,
    ::RadiationOperator,
    ::AbstractVector,
    ::Real,
    ::Nothing,
    model::RadiationModel,
)
    if model == NO_RADIATION || model == SURFACE_ABSORPTION || model == P1_QUASI_STEADY
        return nothing
    end
    if model == BEER_LAMBERT
        error(
            "RadiationOperator(BEER_LAMBERT): contribute_to_jacobian! requires " *
            "the Beer-Lambert workspace overload. Dispatcher should pass " *
            "`BeerLambertWorkspace` rather than `nothing`.",
        )
    end
    error("RadiationOperator: unknown radiation model $(model).")
end

"""
	contribute_to_jacobian!(jac, op::RadiationOperator, u, t, ws,
	                        state, mesh, prop_cache, material, layout, I_radiation,
	                        model::RadiationModel, apply_kirchhoff_absorptivity) -> nothing

BEER_LAMBERT-specific overload. Writes the volumetric absorption Jacobian
contribution to the energy row.

Decomposition:
- Diagonal `j == i` (cell-local in `(T_i, ξ_·,i)`): captured via
  ForwardDiff over the closure
  `(T, ξ_·) ↦ I_in[i+1] · (1 − exp(−α(T, ξ) · Δz)) / (Δz · ρcp_eff(T, ξ))`,
  where `I_in[i+1]` is held fixed (Float64 closure capture). This gives
  `∂(q_rad/ρcp)/∂x_i` in one pass, including both the α-diagonal term
  `(I_in − q_rad·Δz) · ∂α/∂x_i / ρcp` and the ρcp chain
  `−q_rad/ρcp² · ∂ρcp/∂x_i`. On the Kirchhoff path the surface cell's
  closure instead recomputes the live transmissivity,
  `a_eff(T, ξ)·I_radiation·(1 − e^{−τ})/(Δz·ρcp)`, so the local
  `∂a_eff/∂(T, ξ_·)` enters its diagonal exactly.
- Suffix `j > i`:
  `∂(q_rad[i]/ρcp[i])/∂x_j = (−q_rad[i]/ρcp[i]) · Δz_j · ∂α_j/∂x_j`.
  Pushed onto `jac.structured` as one
  [`EmbeddedScaledSuffixSumCoupling`](@ref) per channel (T + each ξⱼ),
  with `row_scale[i] = −q_rad[i]/ρcp_eff[i]` and `v[k] = Δz_k · ∂α_k/∂x_k`.
  On the Kirchhoff path every row `i < n` also depends on the surface
  cell through the injected `a_eff(T_n, ξ_·,n)·I_radiation`; since
  `q_rad[i] ∝ a_eff`, the chain `∂(q_rad[i]/ρcp_i)/∂x_n +=
  (q_rad[i]/ρcp_i)·(∂a_eff/∂x_n)/a_eff` folds into the same couplings as
  `v[n] -= (∂a_eff/∂x_n)/a_eff`.
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::RadiationOperator,
    u::AbstractVector,
    t::Real,
    ws::BeerLambertWorkspace,
    state,
    mesh,
    prop_cache,
    material::Material{NC,NR},
    layout,
    I_radiation::Real,
    model::RadiationModel,
    apply_kirchhoff_absorptivity::Bool;
    assembly_profile::Symbol = :full,
) where {NC,NR}

    if model == NO_RADIATION || model == SURFACE_ABSORPTION || model == P1_QUASI_STEADY
        return nothing
    end
    model == BEER_LAMBERT || error("RadiationOperator: unknown radiation model $(model).")

    n_cells = layout.n_cells
    n_live = mesh.n_active_cells
    A_section = Float64(state.A_section[])
    A_jac = jac.A
    smap = jac.smap
    is_full = assembly_profile === :full

    # ---------------------------------------------------------------
    # Diagonal contribution: cell-local AD over (T_i, ξ_·,i) with I_in
    # held fixed. Captures the α-diagonal `(I_in - q_rad·Δz) · ∂α/∂x_i`
    # and the ρcp chain `-q_rad/ρcp² · ∂ρcp/∂x_i` in one pass. The
    # Kirchhoff surface cell recomputes a_eff(T, ξ) live (with the bare
    # I_radiation as base) so its diagonal carries ∂a_eff/∂x_n too.
    # ---------------------------------------------------------------
    x_buf = ws.x_buf
    f_buf = ws.f_buf
    J_local = ws.J_local
    cfg = ws.cfg

    @inbounds for i = 1:n_live
        Δz_i = state.volume[i] / A_section
        surface_kirchhoff = apply_kirchhoff_absorptivity && i == n_live
        # ws.I_in[i] = intensity entering cell i from above (includes a_eff·I
        # at the surface on the Kirchhoff path). The Kirchhoff surface cell
        # uses the bare I_radiation as base and re-attaches a_eff(T, ξ) inside
        # the closure so AD sees the live dependence.
        I_base_i = surface_kirchhoff ? Float64(I_radiation) : ws.I_in[i]

        x_buf[1] = state.T[i]
        for j = 1:NC
            x_buf[1+j] = state.ξ[j, i]
        end

        f! =
            let I_base = I_base_i,
                Δz = Δz_i,
                material = material,
                surface_kirchhoff = surface_kirchhoff

                (out, x_inner) -> begin
                    T = x_inner[1]
                    ξ = ntuple(j -> x_inner[1+j], Val(NC))
                    α = effective_absorption(ξ, material)
                    ρcp = effective_heat_capacity(ξ, T, material)
                    τ = α * Δz
                    scale =
                        surface_kirchhoff ? effective_absorptivity(ξ, T, material) :
                        one(ρcp)
                    # Numerically stable q_rad = scale·I_base·(1 - exp(-τ)) / Δz; for
                    # very small τ this collapses to scale·I_base·α (Taylor expansion).
                    # Pure expression is fine for ForwardDiff.
                    q_rad = scale * I_base * (1 - exp(-τ)) / Δz
                    out[1] = q_rad / ρcp
                end
            end
        ForwardDiff.jacobian!(J_local, f!, f_buf, x_buf, cfg)

        T_row = i
        fast_sparse_add!(A_jac, smap, T_row, i, J_local[1, 1])
        for k = 1:NC
            ξk_col = n_cells + (k - 1) * n_cells + i
            fast_sparse_add!(A_jac, smap, T_row, ξk_col, J_local[1, 1+k])
        end
    end

    # ---------------------------------------------------------------
    # Suffix contribution: j > i via EmbeddedScaledSuffixSumCoupling.
    # T-channel + per-ξⱼ-channel. Skipped under :ale_reduced — the
    # `EmbeddedScaledSuffixSumCoupling`s would be dropped by
    # `_emit_jacobian!`, and the suffix sparsity slots that would receive
    # their contribution are not emitted under that profile.
    # ---------------------------------------------------------------
    if is_full
        n_global = size(A_jac, 1)
        # row_scale[i] = -q_rad[i] / ρcp_eff[i]. One workspace buffer aliased
        # by every suffix coupling pushed below (they only read it).
        row_scale = ws.row_scale
        @inbounds for i = 1:n_live
            ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
            ρcp_i = effective_heat_capacity(ξ_i, state.T[i], material)
            row_scale[i] = -ws.q_rad[i] / ρcp_i
        end
        @inbounds for i = (n_live+1):n_cells
            row_scale[i] = 0.0
        end

        # Kirchhoff a_eff chain at the surface column: q_rad[i] ∝ a_eff(x_n),
        # so v[n] -= (∂a_eff/∂x_n)/a_eff per channel (see docstring). Skipped
        # when a_eff = 0 — then q_rad ≡ 0 and the couplings vanish anyway.
        # ∂a_eff/∂T_n is non-zero for T-dependent component ε (Kirchhoff
        # components) or densities.
        da_dT_n = 0.0
        da_dξ_n = ntuple(_ -> 0.0, Val(NC))
        a_n = 1.0
        if apply_kirchhoff_absorptivity && n_live > 0
            a_n = Float64(prop_cache.a[n_live])
            if a_n > 0.0
                ξ_n = ntuple(j -> state.ξ[j, n_live], Val(NC))
                T_d, ξ_d = _stack_dual_seed(_BLPropTag, state.T[n_live], ξ_n)
                a_d = effective_absorptivity(ξ_d, T_d, material)
                if a_d isa ForwardDiff.Dual
                    p = ForwardDiff.partials(a_d)
                    da_dT_n = p[1]
                    da_dξ_n = ntuple(j -> p[1+j], Val(NC))
                end
            end
        end
        has_a_chain = apply_kirchhoff_absorptivity && a_n > 0.0

        # T-channel
        if any(!iszero, ws.dα_dT) || (has_a_chain && da_dT_n != 0.0)
            v_T = ws.v_channels[1]
            @inbounds for k = 1:n_live
                Δz_k = state.volume[k] / A_section
                v_T[k] = Δz_k * ws.dα_dT[k]
            end
            @inbounds for k = (n_live+1):n_cells
                v_T[k] = 0.0
            end
            has_a_chain && (v_T[n_live] -= da_dT_n / a_n)
            push!(
                jac.structured,
                EmbeddedScaledSuffixSumCoupling(v_T, row_scale, 0, 0, n_cells, n_global),
            )
        end
        # ξⱼ-channels
        for k_species = 1:NC
            v_ξ = ws.v_channels[1+k_species]
            any_nonzero = false
            @inbounds for k = 1:n_live
                Δz_k = state.volume[k] / A_section
                val = Δz_k * ws.dα_dξ[k_species, k]
                val != 0.0 && (any_nonzero = true)
                v_ξ[k] = val
            end
            @inbounds for k = (n_live+1):n_cells
                v_ξ[k] = 0.0
            end
            if has_a_chain && da_dξ_n[k_species] != 0.0
                v_ξ[n_live] -= da_dξ_n[k_species] / a_n
                any_nonzero = true
            end
            any_nonzero || continue
            ξk_col_off = n_cells + (k_species - 1) * n_cells
            push!(
                jac.structured,
                EmbeddedScaledSuffixSumCoupling(
                    v_ξ,
                    row_scale,
                    0,
                    ξk_col_off,
                    n_cells,
                    n_global,
                ),
            )
        end
    end

    return nothing
end
