"""
Sparsity-pattern construction.

Builds **two patterns** keyed on a `StateLayout` + `SimulationMode`:

- `sparsity_phys::SparseMatrixCSC{Bool,Int}` — every coupling that does not
  involve a `z` column. Block-tridiagonal physics + Beer–Lambert non-local
  + the χ block (rows and dense χ-column couplings on T/ξ/χ rows).
- `sparsity_geom::SparseMatrixCSC{Bool,Int}` — every coupling that touches
  a `z` column or row. Empty (no nnz, full-size) under Eulerian.

Both patterns share the full state-vector size so backends can scatter into
either with absolute indices, and `S_phys ∪ S_geom` recovers the full
Jacobian sparsity. The two patterns share the same `SparseMatrixCSC{Bool,Int}`
type — the JacobianSpec contract requires it (see `interface.jl` field
docs). Eulerian's geometry pattern is encoded as `spzeros(Bool, n, n)` rather
than `nothing`.

The four `_emit_*_pattern!` helpers each push (row, col) pairs into a single
`(I, J)` pair. The orchestrator `build_sparsity` calls them in the right
order to assemble the two output matrices. Splitting by emitter (rather
than by output target) keeps each helper's logic local to one physics
block; `build_sparsity` is the only file that knows which block goes where.
"""

using ..Materials: n_components
using ..Problem:
    SimulationMode, StateLayout, mesh_mode, radiation_mode, chi_mode, state_length
using ..Problem:
    Eulerian,
    ALE,
    NoRadiation,
    SurfaceAbsorption,
    BeerLambert,
    Fickian,
    DarcyFick,
    NoChi,
    WithChi

# ============================================================================
# Index helpers — single source of truth for the column/row arithmetic
# ============================================================================

@inline _temp_idx(::StateLayout, i::Int) = i

@inline _xi_idx(layout::StateLayout, j::Int, i::Int) =
    layout.n_cells + (j - 1) * layout.n_cells + i

@inline function _z_offset(layout::StateLayout)
    NC = n_components(layout)
    return layout.n_cells * (1 + NC)
end

@inline _z_idx(layout::StateLayout, n::Int) = _z_offset(layout) + n

@inline function _chi_offset(layout::StateLayout, mode::SimulationMode)
    off = _z_offset(layout)
    mesh_mode(mode) isa ALE && (off += layout.n_nodes)
    return off
end

@inline _chi_idx(layout::StateLayout, mode::SimulationMode, i::Int) =
    _chi_offset(layout, mode) + i

# ============================================================================
# Emitter 1: physics block-tridiagonal couplings
# ============================================================================

"""
	_emit_physics_pattern!(I, J, layout, mode)

Emit (row, col) indices for the block-tridiagonal physics coupling
(T↔T, T↔ξ, ξ↔T, ξ↔ξ within each cell's tridiagonal stencil). Always emits;
the pattern shape does not depend on `mode`.
"""
function _emit_physics_pattern!(
    I::Vector{Int},
    J::Vector{Int},
    layout::StateLayout,
    ::SimulationMode,
)
    n = layout.n_cells
    NC = n_components(layout)
    for i = 1:n
        T_row = i

        # T row ← T_neighbor and T row ← ξ_{j,neighbor}
        for neighbor = max(1, i-1):min(n, i+1)
            push!(I, T_row);
            push!(J, neighbor)
            for j = 1:NC
                push!(I, T_row);
                push!(J, _xi_idx(layout, j, neighbor))
            end
        end

        # ξ rows for this cell
        for j = 1:NC
            ξ_row = _xi_idx(layout, j, i)
            for neighbor = max(1, i-1):min(n, i+1)
                push!(I, ξ_row);
                push!(J, neighbor)
                for k = 1:NC
                    push!(I, ξ_row);
                    push!(J, _xi_idx(layout, k, neighbor))
                end
            end
        end
    end
    return nothing
end

# ============================================================================
# Emitter 2: ALE geometry couplings (physics ↔ z, z ↔ z)
# ============================================================================

"""
	_emit_geometry_pattern!(I, J, layout, mode)

Emit indices for ALE-only geometry couplings:
- physics → z: cell `i` reads nodes `1..min(n_nodes, i+3)`. The lower bound
  is `1` because ALE advection's `w_cell[i] = ½(w_node[n_left]+w_node[n_right])`
  with `w_node[k] = Σⱼ<k θⱼ Δzⱼ` gives an exact dependence on every upstream
  `z[1..i+1]`. The upper bound `i+3` covers the TVD-limited gradient's read
  of `cell_centers[i+2] = ½(z[i+2]+z[i+3])` when the donor-donor cell is
  `i+2` (the w_cell > 0 branch — the donor side follows c = -w); the
  symmetric `cell_centers[i-2]` read on the w_cell < 0 branch is subsumed
  by the lower bound `1`.
- z → physics: lower triangular. Node `n` velocity sums strain rates over
  cells `1..n-1`.
- z → z: lower triangular (including diagonal). `wₖ = Σᵢ₌₁ᵏ⁻¹ θᵢ Δzᵢ` and
  `Δzᵢ = z_{i+1} - zᵢ` give `∂wₖ/∂zₘ ≠ 0` for `m ≤ k`.

**WithChi widening.** When `chi_mode(mode) isa WithChi`, the column-global
`θ_A = law'(χ̄)·(dχ̄/dt)/law(χ̄)` enters every `w_node[k] = Σᵢ<k(θ_pyr,i − θ_A)·Δzᵢ`.
Since `dχ̄/dt = (1/M)·Σⱼ S_dry,j·V_j` with `V_j = A_j·Δz_j` depends on every
`z_m`, every `w_node[k]` (and hence every `w_cell[i]`) acquires a non-trivial
derivative w.r.t. every `z_m`. The same holds for the physics → χ → χ̄ → θ_A
column-global path. Net effect: physics → z, z → physics, and z → z all
become **full bands** rather than the bounded/triangular shapes above.

No-op under Eulerian.
"""
function _emit_geometry_pattern!(
    I::Vector{Int},
    J::Vector{Int},
    layout::StateLayout,
    mode::SimulationMode,
)
    mesh_mode(mode) isa ALE || return nothing

    n = layout.n_cells
    NC = n_components(layout)
    n_nodes = layout.n_nodes
    is_with_chi = chi_mode(mode) isa WithChi

    # physics → z: full cumulative reach 1..min(n_nodes, i+3); WithChi widens
    # this to 1..n_nodes via the θ_A path (column-global χ̄ derivative).
    for i = 1:n
        T_row = i
        n_end = is_with_chi ? n_nodes : min(n_nodes, i + 3)
        for nn = 1:n_end
            z_col = _z_idx(layout, nn)
            push!(I, T_row);
            push!(J, z_col)
        end
        for j = 1:NC
            ξ_row = _xi_idx(layout, j, i)
            for nn = 1:n_end
                z_col = _z_idx(layout, nn)
                push!(I, ξ_row);
                push!(J, z_col)
            end
        end
    end

    # z → physics (lower triangular; node 1 is fixed → no deps).
    # WithChi widens to full bands: every z_row k depends on every cell `i`'s
    # (T_i, ξ_·,i) via the θ_A column-global path, plus the standard
    # prefix-sum below k.
    for nn = 2:n_nodes
        z_row = _z_idx(layout, nn)
        i_hi = is_with_chi ? n : (nn - 1)
        for i = 1:i_hi
            push!(I, z_row);
            push!(J, i)
            for j = 1:NC
                push!(I, z_row);
                push!(J, _xi_idx(layout, j, i))
            end
        end
    end

    # z → z (lower triangular including diagonal). WithChi widens to full
    # bands: every `w_node[k]` depends on every `z_m` via θ_A's dependence on
    # `dχ̄/dt = (1/M)·Σⱼ S_dry,j·A_j·Δz_j`.
    for nn = 1:n_nodes
        z_row = _z_idx(layout, nn)
        m_hi = is_with_chi ? n_nodes : nn
        for m = 1:m_hi
            z_col = _z_idx(layout, m)
            push!(I, z_row);
            push!(J, z_col)
        end
    end

    return nothing
end

# ============================================================================
# Emitter 3: χ block (variable cross-section)
# ============================================================================

"""
	_emit_chi_pattern!(I, J, layout, mode, profile = :full)

Emit χ-row dependencies (T_i, ξ_·,i, z_i, z_{i+1}) and the dense χ-column
coupling (every T/ξ/z/χ row depends on every χⱼ when A depends on χ̄).
No-op when WithChi is not in `mode`.

Under `profile = :full` the full χ block is emitted — row deps plus dense
column couplings. Under `profile = :ale_reduced`:

- The cell-local χ-row × {T_i, ξ_·,i} entries are still emitted (written by
  `ChiOperator` regardless of profile).
- The χ-row × {z_i, z_{i+1}} entries are dropped: under `:ale_reduced`,
  `GeometryOperator` is skipped, and `ChiOperator` defers z-column writes
  to the geometry AD pass — so under this profile no operator writes those
  slots, and emitting them strands structurally-allocated zeros that the
  masked-compare verifier flags as missing entries.
- The dense χ-column block is dropped: it is populated only by the rank-1
  `EmbeddedRank1Coupling` that `ChiOperator` pushes onto `jac.structured`,
  which `_emit_jacobian!` itself drops under `:ale_reduced`.
"""
function _emit_chi_pattern!(
    I::Vector{Int},
    J::Vector{Int},
    layout::StateLayout,
    mode::SimulationMode,
    profile::Symbol = :full,
)
    chi_mode(mode) isa WithChi || return nothing

    n = layout.n_cells
    NC = n_components(layout)
    chi_off = _chi_offset(layout, mode)
    is_ale = mesh_mode(mode) isa ALE
    is_full = profile === :full

    # χ-row dependencies for each cell:
    #   χ_i ← T_i, ξ_·,i, z_i, z_{i+1}   (z entries only when ALE)
    for i = 1:n
        χ_row = chi_off + i
        push!(I, χ_row);
        push!(J, i)                       # T_i
        for j = 1:NC
            push!(I, χ_row);
            push!(J, _xi_idx(layout, j, i))  # ξ_{j,i}
        end
        if is_ale && is_full
            push!(I, χ_row);
            push!(J, _z_idx(layout, i))      # z_i
            push!(I, χ_row);
            push!(J, _z_idx(layout, i + 1))  # z_{i+1}
        end
    end

    # Dense χ-column coupling: T, ξ, z, χ rows all depend on every χⱼ.
    # (z node 1 is fixed → omit.) Populated by EmbeddedRank1Coupling on
    # jac.structured, which is dropped under :ale_reduced.
    if is_full
        for j = 1:n
            χ_col = chi_off + j

            # T rows
            for i = 1:n
                push!(I, i);
                push!(J, χ_col)
            end

            # ξ rows (each component × each cell)
            for k = 1:NC
                for i = 1:n
                    push!(I, _xi_idx(layout, k, i));
                    push!(J, χ_col)
                end
            end

            # z rows (k = 2..n_nodes; node 1 fixed)
            if is_ale
                for k = 2:layout.n_nodes
                    push!(I, _z_idx(layout, k));
                    push!(J, χ_col)
                end
            end

            # χ rows
            for i = 1:n
                push!(I, chi_off + i);
                push!(J, χ_col)
            end
        end
    end

    return nothing
end

# ============================================================================
# Emitter 4: Beer-Lambert non-local coupling
# ============================================================================

"""
	_emit_radiation_pattern!(I, J, layout, mode, profile = :full)

Emit non-local Beer–Lambert coupling (T row at cell `i` depends on ξ_{·,k}
for every cell `k ≥ i` — radiation enters from cell `n` (surface) and
attenuates downward, so cell `i` sees the optical depth contributed by
every upstream cell). No-op for `NoRadiation`, `SurfaceAbsorption`, or any
other `RadiationMode` that returns `false` from non-local coupling.

Under `profile = :full` the full T_row × ξ_{j, k ≥ i} block is emitted,
plus the T_row × T_col(n) surface column: on the Kirchhoff default path
the injected intensity is `ε_eff(T_n, ξ_·,n)·I`, so every energy row
additionally depends on the surface cell's temperature through
`∂ε_eff/∂T_n` (the ξ_·,n dependence is already covered by the `k = n`
ξ-columns).
Under `profile = :ale_reduced` only the diagonal `k = i` is emitted: the
suffix `k > i` (and the ε_eff surface column) is populated by
`EmbeddedScaledSuffixSumCoupling` pushes that `RadiationOperator`
performs only under `:full`.
"""
function _emit_radiation_pattern!(
    I::Vector{Int},
    J::Vector{Int},
    layout::StateLayout,
    mode::SimulationMode,
    profile::Symbol = :full,
)
    radiation_mode(mode) isa BeerLambert || return nothing

    n = layout.n_cells
    NC = n_components(layout)
    is_full = profile === :full
    for i = 1:n
        T_row = i
        k_hi = is_full ? n : i
        for k = i:k_hi
            for j = 1:NC
                push!(I, T_row);
                push!(J, _xi_idx(layout, j, k))
            end
        end
    end
    # Kirchhoff surface-transmissivity chain: T_row(i) × T_col(n) for i < n
    # (the i = n diagonal already exists in the block-tridiagonal pattern).
    if is_full
        for i = 1:(n-1)
            push!(I, i);
            push!(J, n)
        end
    end
    return nothing
end

# ============================================================================
# Emitter 5: ALE advection (5-cell limiter window + cumulative coupling)
# ============================================================================

"""
	_emit_ale_advection_pattern!(I, J, layout, mode, profile)

Emit the additional non-local couplings introduced by
[`Pyrolysis.Jacobian.ALEAdvectionOperator`](@ref):

1. **5-cell limiter window** — the direct chain rule pulls in `(i, upwind,
   upup)` × both upwind directions. The reach for row `i` is cells
   `max(1, i-2):min(n, i+2)` (the operator declares `stencil = (cells = -2:2,)`).
2. **Cumulative mesh-velocity coupling** — the indirect chain rule through
   `w_cell = 0.5·(w_node[n_left] + w_node[n_right])` connects row `i` to
   every upstream cell `j ≤ max(n_left[i], n_right[i]) - 1 ≈ i`. For a
   1-D sequentially-numbered mesh this is `1:i`.

Under `profile = :full` the union for row `i` is `1:min(n, i+2)` —
strictly lower-triangular plus two superdiagonals. **WithChi widening:**
the column-global `θ_A(χ̄, T, ξ)` enters every `w_cell[i]` through
`dχ̄/dt = Σⱼ S_dry,j(T_j, ξ_j)·V_j/M`, so every advection row acquires a
rank-1 dependence on every cell's `(T, ξ)` — the reach becomes the full
band `1:n`. Under `profile = :ale_reduced` only the bounded 5-cell
limiter window `max(1, i-2):min(n, i+2)` is emitted; the cumulative and
θ_A widenings are dropped because the consumer (`ALEAdvectionOperator`)
skips its `compose_into!` calls under the same profile, so those entries
would never receive a write anyway.

Applies to T rows and ξ rows of every component (over-emitting
solid-component ξ rows is safe; the operator writes zero there and the
storage cost is negligible).

ALE-only; no-op under Eulerian.
"""
function _emit_ale_advection_pattern!(
    I::Vector{Int},
    J::Vector{Int},
    layout::StateLayout,
    mode::SimulationMode,
    profile::Symbol = :full,
)
    mesh_mode(mode) isa ALE || return nothing

    n = layout.n_cells
    NC = n_components(layout)
    is_reduced = profile === :ale_reduced
    is_with_chi = chi_mode(mode) isa WithChi

    for i = 1:n
        # Column reach for row `i`:
        #   :full               → 1..min(n, i+2)  (cumulative + limiter)
        #   :full + WithChi     → 1..n            (θ_A rank-1 reaches every cell)
        #   :ale_reduced        → max(1, i-2)..min(n, i+2)  (limiter only)
        j_lo = is_reduced ? max(1, i - 2) : 1
        j_hi = (!is_reduced && is_with_chi) ? n : min(n, i + 2)

        # T row
        T_row = i
        for j = j_lo:j_hi
            push!(I, T_row);
            push!(J, j)                       # T col
            for k = 1:NC
                push!(I, T_row);
                push!(J, _xi_idx(layout, k, j))  # ξₖ col
            end
        end

        # ξ rows for this cell — over-emit for all components (gas + solid)
        # because the operator writes only to gas rows and `sparse(I, J, …, |)`
        # deduplicates the resulting union with `_emit_physics_pattern!`. The
        # extra solid-row entries land as structural zeros that the assembly
        # never touches.
        for kk = 1:NC
            ξ_row = _xi_idx(layout, kk, i)
            for j = j_lo:j_hi
                push!(I, ξ_row);
                push!(J, j)
                for k = 1:NC
                    push!(I, ξ_row);
                    push!(J, _xi_idx(layout, k, j))
                end
            end
        end
    end
    return nothing
end

# ============================================================================
# Orchestrator
# ============================================================================

"""
	build_sparsity(layout::StateLayout, mode::SimulationMode)
	    -> (S_phys, S_geom)

Compose the emitters into two patterns:

- `S_phys` carries every coupling whose **column** is a T, ξ, or χ index:
  block-tridiagonal physics (`_emit_physics_pattern!`), Beer-Lambert
  non-local (`_emit_radiation_pattern!`), the χ block (`_emit_chi_pattern!`),
  and (ALE only) the ALE-advection limiter+cumulative widening
  (`_emit_ale_advection_pattern!`).
- `S_geom` carries every coupling whose **column** is a z index: the ALE
  geometry block (`_emit_geometry_pattern!`).

For Eulerian, `S_geom` is an empty `SparseMatrixCSC{Bool,Int}` of full
size (no nnz). The shared `SparseMatrixCSC{Bool,Int}` type is required by
[`JacobianSpec`](@ref)'s field-type-sharing contract.

`S_phys ∪ S_geom` recovers the full Jacobian sparsity for any
`(layout, mode)`.
"""
# Lazy sparsity: the (S_phys, S_geom) pattern depends only on
# `layout` (n_cells, n_nodes, NC), `mode` (a singleton type), and the
# assembly profile. We memoize on a structural fingerprint so repeat
# `solve()` invocations on the same problem skip the rebuild entirely;
# AMR resizes the workspace (changing the fingerprint), so cache
# invalidation is automatic.
const _SPARSITY_CACHE = Dict{
    Tuple{DataType,Int,Int,Int,Symbol},
    Tuple{SparseMatrixCSC{Bool,Int},SparseMatrixCSC{Bool,Int}},
}()

@inline function _sparsity_key(
    layout::StateLayout{NC},
    mode::SimulationMode,
    profile::Symbol,
) where {NC}
    return (typeof(mode), NC, layout.n_cells, layout.n_nodes, profile)
end

function build_sparsity(layout::StateLayout, mode::SimulationMode; profile::Symbol = :full)
    profile === :full ||
        profile === :ale_reduced ||
        throw(
            ArgumentError(
                "build_sparsity: unknown profile :$profile (expected :full or :ale_reduced)",
            ),
        )
    key = _sparsity_key(layout, mode, profile)
    cached = get(_SPARSITY_CACHE, key, nothing)
    cached !== nothing && return cached

    result = _build_sparsity_impl(layout, mode, profile)
    _SPARSITY_CACHE[key] = result
    return result
end

function _build_sparsity_impl(layout::StateLayout, mode::SimulationMode, profile::Symbol)
    n_total = state_length(layout)

    # S_phys: physics + radiation + chi + ALE advection widening
    I_phys = Int[]
    J_phys = Int[]
    _emit_physics_pattern!(I_phys, J_phys, layout, mode)
    _emit_radiation_pattern!(I_phys, J_phys, layout, mode, profile)
    _emit_chi_pattern!(I_phys, J_phys, layout, mode, profile)
    _emit_ale_advection_pattern!(I_phys, J_phys, layout, mode, profile)
    S_phys = sparse(I_phys, J_phys, trues(length(I_phys)), n_total, n_total, |)

    # S_geom: pure z-column couplings (ALE only). Empty under :ale_reduced
    # regardless of mode — the geometry block is the entire reason the
    # reduced profile exists.
    if mesh_mode(mode) isa ALE && profile === :full
        I_geom = Int[]
        J_geom = Int[]
        _emit_geometry_pattern!(I_geom, J_geom, layout, mode)
        S_geom = sparse(I_geom, J_geom, trues(length(I_geom)), n_total, n_total, |)
    else
        S_geom = spzeros(Bool, Int, n_total, n_total)
    end

    return (S_phys, S_geom)
end

# ============================================================================
# Greedy column coloring
# ============================================================================

"""
	matrix_colors(sparsity::SparseMatrixCSC) -> Vector{Int}

Greedy column coloring for finite-difference / sparse-AD Jacobian assembly.
Columns sharing a color have no common non-zero rows, so they can be
perturbed simultaneously by a single residual evaluation.

Uses an O(nnz) algorithm: builds a CSR representation from the CSC sparsity
to efficiently find the columns that share a row with each new column,
then assigns the minimum-available color.

Lives in the Jacobian module so every backend shares one coloring
implementation.
"""
function matrix_colors(sparsity::SparseMatrixCSC)
    n = size(sparsity, 2)
    m = size(sparsity, 1)

    row_ptr = zeros(Int, m + 1)
    @inbounds for idx = 1:length(sparsity.rowval)
        row_ptr[sparsity.rowval[idx]+1] += 1
    end
    row_ptr[1] = 1
    @inbounds for i = 1:m
        row_ptr[i+1] += row_ptr[i]
    end

    col_indices = Vector{Int}(undef, length(sparsity.rowval))
    row_count = zeros(Int, m)
    @inbounds for j = 1:n
        for idx = sparsity.colptr[j]:(sparsity.colptr[j+1]-1)
            row = sparsity.rowval[idx]
            pos = row_ptr[row] + row_count[row]
            col_indices[pos] = j
            row_count[row] += 1
        end
    end

    colors = zeros(Int, n)
    max_colors = 64
    forbidden = falses(max_colors)
    for j = 1:n
        fill!(forbidden, false)
        @inbounds for idx = sparsity.colptr[j]:(sparsity.colptr[j+1]-1)
            row = sparsity.rowval[idx]
            for k = row_ptr[row]:(row_ptr[row+1]-1)
                other_col = col_indices[k]
                if other_col < j
                    c = colors[other_col]
                    if c > 0
                        if c > length(forbidden)
                            old_len = length(forbidden)
                            resize!(forbidden, 2 * c)
                            forbidden[(old_len+1):end] .= false
                        end
                        forbidden[c] = true
                    end
                end
            end
        end
        color = 1
        while color <= length(forbidden) && forbidden[color]
            color += 1
        end
        colors[j] = color
    end
    return colors
end
