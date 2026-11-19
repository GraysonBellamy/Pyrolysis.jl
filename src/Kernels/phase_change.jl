# ============================================================================
# Phase changes
# ============================================================================

"""
	PhaseChangeKernel <: Kernel

The two-way phase changes of a material's liquids into their vapours
([`PhaseChange`](@ref)): per live cell `i` and phase change, the rate

	ω_v = K (ξ_v^eq − ξ_v),   ξ_v^eq = φ M a(X, T) p_sat(T)/(R T),

with the moisture content `X = ξ_w/(Σₛ wₛ ξ⁺ₛ + 10⁻³ kg/m³)`, the signed
concentrations `ξ_w`, `ξ_v` of the liquid and the vapour, the positive-part
concentrations `ξ⁺ₛ` of the sorbing solids, the trial temperature `T` and the
guarded pore fraction `φ` of the derived state; it adds `−ω_v Vᵢ` to the
liquid's mass source, `+ω_v Vᵢ` to the vapour's and `−h_pc ω_v Vᵢ` to the
energy source, with the heat `h_pc = L_ref + h_v − h_w + Q_s(X, T)` built
from the cell's own sensible enthalpies `h_v`, `h_w`, so that the formation
term of the inventory form cancels them
([`phase_change_rate`](@ref), [`phase_change_heat`](@ref)).

The signed liquid and vapour make both ends restoring: a negative trial liquid
mass gives a negative activity and draws vapour back, a negative trial vapour
mass raises the rate. Dead cells contribute nothing. The kernel holds no
per-cell table; `setup` places it after the reaction kernel, for a material
with phase changes only.
"""
struct PhaseChangeKernel <: Kernel end

function contribute!(
    c::Contributions{NC},
    ::PhaseChangeKernel,
    ds::DerivedState{NC},
    topo::MeshTopology,
    material::Material{NC},
    t,
) where {NC}
    @inbounds for i = 1:n_cells(topo)
        ds.active[i] || continue
        ξ = ntuple(j -> ds.ξ[j, i], Val(NC))
        ξ⁺ = ntuple(j -> ds.ξ⁺[j, i], Val(NC))
        _phase_changes!(c, material.phase_changes, ds, i, ξ, ξ⁺)
    end
    return c
end

@inline _phase_changes!(c, ::Tuple{}, ds, i, ξ, ξ⁺) = nothing
@inline function _phase_changes!(c, changes::Tuple, ds, i, ξ, ξ⁺)
    _phase_change!(c, first(changes), ds, i, ξ, ξ⁺)
    _phase_changes!(c, Base.tail(changes), ds, i, ξ, ξ⁺)
    return nothing
end

@inline function _phase_change!(c, pc, ds, i, ξ, ξ⁺)
    w, v = pc.liquid, pc.vapour
    T = ds.T_trial[i]
    X = moisture_content(pc, ξ, ξ⁺)
    @inbounds begin
        ω = phase_change_rate(pc, X, T, ds.φ[i], ξ[v]) * ds.V_eval[i]
        h = phase_change_heat(pc, X, T, ds.h[v, i], ds.h[w, i])
        c.S_mass[w, i] -= ω
        c.S_mass[v, i] += ω
        c.S_energy[i] -= h * ω
    end
    return nothing
end
