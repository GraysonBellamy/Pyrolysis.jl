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

Every row of `S_phys ∪ S_geom` stores its diagonal, in every profile, so the
shifted matrix `W = J − M/(γΔt)` the integrator factors has exactly the
pattern of `J`. Where no operator writes the diagonal (the ALE node rows and
the χ rows under `:ale_reduced`) it holds a structural zero.

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

@inline _xi_idx(layout::StateLayout, j::Int, i::Int) =
    layout.n_cells + (j - 1) * layout.n_cells + i

# Node rows/columns exist only under ALE, where `z_offset` is the end of the
# ξ block.
@inline _z_idx(layout::StateLayout, n::Int) = z_offset(layout) + n

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
            push!(I, T_row)
            push!(J, neighbor)
            for j = 1:NC
                push!(I, T_row)
                push!(J, _xi_idx(layout, j, neighbor))
            end
        end

        # ξ rows for this cell
        for j = 1:NC
            ξ_row = _xi_idx(layout, j, i)
            for neighbor = max(1, i-1):min(n, i+1)
                push!(I, ξ_row)
                push!(J, neighbor)
                for k = 1:NC
                    push!(I, ξ_row)
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
- physics → z: cell `i` is given nodes `1..min(n_nodes, i+3)`. The residual
  rows read `z[i-1..i+2]` (the cell volume and the face distances to the
  neighbouring cell centres); the TVD-limited temperature advection of the
  verification `w_override` hook reads `cell_centers[i±2]`, i.e.
  `z[i-2..i+3]`. The emitted band is a superset of both.
- z → physics: lower triangular. Node `n` velocity sums strain rates over
  cells `1..n-1`.
- z → z: lower triangular (including diagonal). `wₖ = Σᵢ₌₁ᵏ⁻¹ θᵢ Δzᵢ` and
  `Δzᵢ = z_{i+1} - zᵢ` give `∂wₖ/∂zₘ ≠ 0` for `m ≤ k`.

**WithChi widening.** When `chi_mode(mode) isa WithChi`, the column-global
`θ_A = law'(χ̄)·(dχ̄/dt)/law(χ̄)` enters every `w_node[k] = Σᵢ<k(θ_pyr,i − θ_A)·Δzᵢ`.
Since `dχ̄/dt = (1/M)·Σⱼ S_dry,j·V_j` with `V_j = A_j·Δz_j` depends on every
`z_m`, every `w_node[k]` acquires a non-trivial
derivative w.r.t. every `z_m`. The same holds for the physics → χ → χ̄ → θ_A
column-global path. Net effect: physics → z, z → physics, and z → z all
become **full bands** rather than the bounded/triangular shapes above.

**Beer–Lambert widening.** Under `BeerLambert` the intensity reaching cell `i`
is attenuated by the optical depth `τᵢ = Σₖ>ᵢ αₖ·Δzₖ` of every cell above it,
so the energy row of cell `i` depends on every node `z[i+1..n_nodes]`. The T
rows therefore reach `1..n_nodes`. The ξ rows keep the bounded reach: volumetric
absorption enters only the energy equation.

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
    is_beer_lambert = radiation_mode(mode) isa BeerLambert

    # physics → z: full cumulative reach 1..min(n_nodes, i+3); WithChi widens
    # this to 1..n_nodes via the θ_A path (column-global χ̄ derivative), and
    # Beer–Lambert widens the T rows to 1..n_nodes via the optical depth above
    # the cell.
    for i = 1:n
        T_row = i
        n_end = is_with_chi ? n_nodes : min(n_nodes, i + 3)
        n_end_T = is_beer_lambert ? n_nodes : n_end
        for nn = 1:n_end_T
            z_col = _z_idx(layout, nn)
            push!(I, T_row)
            push!(J, z_col)
        end
        for j = 1:NC
            ξ_row = _xi_idx(layout, j, i)
            for nn = 1:n_end
                z_col = _z_idx(layout, nn)
                push!(I, ξ_row)
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
            push!(I, z_row)
            push!(J, i)
            for j = 1:NC
                push!(I, z_row)
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
            push!(I, z_row)
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
  `ChiOperator` regardless of profile), and so is the structural diagonal
  χ-row × χ_i, which no operator writes under this profile.
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
    chi_off = chi_offset(layout)
    is_full = profile === :full

    # χ-row dependencies for each cell (WithChi is ALE-only):
    #   χ_i ← T_i, ξ_·,i, z_i, z_{i+1}
    # plus the structural diagonal χ_i ← χ_i.
    for i = 1:n
        χ_row = chi_off + i
        push!(I, χ_row)
        push!(J, χ_row)                   # χ_i
        push!(I, χ_row)
        push!(J, i)                       # T_i
        for j = 1:NC
            push!(I, χ_row)
            push!(J, _xi_idx(layout, j, i))  # ξ_{j,i}
        end
        if is_full
            push!(I, χ_row)
            push!(J, _z_idx(layout, i))      # z_i
            push!(I, χ_row)
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
                push!(I, i)
                push!(J, χ_col)
            end

            # ξ rows (each component × each cell)
            for k = 1:NC
                for i = 1:n
                    push!(I, _xi_idx(layout, k, i))
                    push!(J, χ_col)
                end
            end

            # z rows (k = 2..n_nodes; node 1 fixed)
            for k = 2:layout.n_nodes
                push!(I, _z_idx(layout, k))
                push!(J, χ_col)
            end

            # χ rows
            for i = 1:n
                push!(I, chi_off + i)
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
the injected intensity is `a_eff(T_n, ξ_·,n)·I` (the surface absorptivity
mixture, equal to `ε_eff` unless a component sets its own), so every energy
row additionally depends on the surface cell's temperature through
`∂a_eff/∂T_n` (the ξ_·,n dependence is already covered by the `k = n`
ξ-columns).
Under `profile = :ale_reduced` only the diagonal `k = i` is emitted: the
suffix `k > i` (and the a_eff surface column) is populated by
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
                push!(I, T_row)
                push!(J, _xi_idx(layout, j, k))
            end
        end
    end
    # Kirchhoff surface-transmissivity chain: T_row(i) × T_col(n) for i < n
    # (the i = n diagonal already exists in the block-tridiagonal pattern).
    if is_full
        for i = 1:(n-1)
            push!(I, i)
            push!(J, n)
        end
    end
    return nothing
end

# ============================================================================
# Orchestrator
# ============================================================================

# Lazy sparsity: the (S_phys, S_geom) pattern depends only on
# `layout` (n_cells, n_nodes, NC), `mode` (a singleton type), and the
# assembly profile. We memoize on a structural fingerprint so repeat
# `solve()` invocations on the same problem skip the rebuild entirely.
# Concurrent solves share the cache, so every
# access holds `_SPARSITY_CACHE_LOCK`; callers only read the cached patterns.
const _SPARSITY_CACHE = Dict{
    Tuple{DataType,Int,Int,Int,Symbol},
    Tuple{SparseMatrixCSC{Bool,Int},SparseMatrixCSC{Bool,Int}},
}()
const _SPARSITY_CACHE_LOCK = ReentrantLock()

@inline function _sparsity_key(
    layout::StateLayout{NC},
    mode::SimulationMode,
    profile::Symbol,
) where {NC}
    return (typeof(mode), NC, layout.n_cells, layout.n_nodes, profile)
end

"""
	build_sparsity(layout::StateLayout, mode::SimulationMode; profile = :full)
	    -> (S_phys, S_geom)

Compose the emitters into two patterns:

- `S_phys` carries every coupling whose **column** is a T, ξ, or χ index:
  block-tridiagonal physics (`_emit_physics_pattern!`), Beer-Lambert
  non-local (`_emit_radiation_pattern!`), and the χ block
  (`_emit_chi_pattern!`). The ALE dilation `−ξⱼθᵢ` is cell-local and lies
  inside the physics block; the cumulative mesh velocity couples only the
  z rows, which live in `S_geom`.
- `S_geom` carries every coupling whose **column** is a z index: the ALE
  geometry block (`_emit_geometry_pattern!`), or under `:ale_reduced` only
  the structural diagonal of the node rows.

For Eulerian, `S_geom` is an empty `SparseMatrixCSC{Bool,Int}` of full
size (no nnz). The shared `SparseMatrixCSC{Bool,Int}` type is required by
[`JacobianSpec`](@ref)'s field-type-sharing contract.

`S_phys ∪ S_geom` recovers the full Jacobian sparsity for any
`(layout, mode)`.
"""
function build_sparsity(layout::StateLayout, mode::SimulationMode; profile::Symbol = :full)
    profile === :full ||
        profile === :ale_reduced ||
        throw(
            ArgumentError(
                "build_sparsity: unknown profile :$profile (expected :full or :ale_reduced)",
            ),
        )
    key = _sparsity_key(layout, mode, profile)
    return @lock _SPARSITY_CACHE_LOCK get!(_SPARSITY_CACHE, key) do
        _build_sparsity_impl(layout, mode, profile)
    end
end

function _build_sparsity_impl(layout::StateLayout, mode::SimulationMode, profile::Symbol)
    n_total = state_length(layout)

    # S_phys: physics + radiation + chi
    I_phys = Int[]
    J_phys = Int[]
    _emit_physics_pattern!(I_phys, J_phys, layout, mode)
    _emit_radiation_pattern!(I_phys, J_phys, layout, mode, profile)
    _emit_chi_pattern!(I_phys, J_phys, layout, mode, profile)
    S_phys = sparse(I_phys, J_phys, trues(length(I_phys)), n_total, n_total, |)

    # S_geom: pure z-column couplings (ALE only). Under :ale_reduced it
    # holds only the structural diagonal of the node rows — dropping the
    # geometry block is the entire reason the reduced profile exists.
    if mesh_mode(mode) isa ALE
        I_geom = Int[]
        J_geom = Int[]
        if profile === :full
            _emit_geometry_pattern!(I_geom, J_geom, layout, mode)
        else
            for nn = 1:layout.n_nodes
                push!(I_geom, _z_idx(layout, nn))
                push!(J_geom, _z_idx(layout, nn))
            end
        end
        S_geom = sparse(I_geom, J_geom, trues(length(I_geom)), n_total, n_total, |)
    else
        S_geom = spzeros(Bool, Int, n_total, n_total)
    end

    return (S_phys, S_geom)
end

# ============================================================================
# Pattern locality
# ============================================================================

"""
	materializes_geometry_block(mode, profile) -> Bool

Whether the pattern `build_sparsity` emits for `mode` under `profile` holds
the ALE geometry block: under `:full` every cell's rows reach the node
columns below it, and the node rows reach every cell and node below them,
O(N²) entries in all. `solve` factors this pattern with UMFPACK, whose
supernodal factorization handles the dense block, and every other pattern
with KLU.
"""
materializes_geometry_block(mode::SimulationMode, profile::Symbol) =
    profile === :full && mesh_mode(mode) isa ALE

"""
	structural_zero_diagonal(layout, mode, profile) -> Vector{Int}

Rows whose stored diagonal no operator writes under `profile`: the node rows
and the χ rows under `:ale_reduced`, which drops their couplings. The slot
exists only so that every row stores its diagonal, and its value stays zero.
[`verify_jacobian_at`](@ref) leaves these slots out of an approximate
backend's compared support, like every other dropped coupling.
"""
function structural_zero_diagonal(
    layout::StateLayout,
    mode::SimulationMode,
    profile::Symbol,
)
    rows = Int[]
    profile === :ale_reduced || return rows
    if mesh_mode(mode) isa ALE
        append!(rows, _z_idx(layout, 1):_z_idx(layout, layout.n_nodes))
    end
    if chi_mode(mode) isa WithChi
        chi_off = chi_offset(layout)
        append!(rows, (chi_off+1):(chi_off+layout.n_cells))
    end
    return rows
end
