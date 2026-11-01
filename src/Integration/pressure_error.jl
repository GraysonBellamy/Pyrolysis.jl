# Error control uses the proposed state, never changes tolerances inside the
# residual, and includes thermal, condensed-mass and geometric sensitivities.
# The three supported steppers pass their normalized error vector `cache.atmp`
# to internalnorm. Other norm calls (Newton residuals, initial step estimates)
# retain the ordinary RMS norm. FBDF's order-comparison estimates are controlled
# by the same pressure metric as its accepted-order estimate.
struct PressureErrorTag end
const PressureErrorDual=ForwardDiff.Dual{PressureErrorTag,Float64,1}
const _PressureNormWrapper=FunctionWrapper{Float64,Tuple{Vector{Float64},Float64}}
@inline _rms_error(x::AbstractArray, t) = sqrt(sum(abs2, x)/length(x))
@inline _rms_error(x::Number, t) = abs(x)

# Keep the solver-facing type independent of the material and boundary types,
# as for ErasedRHS and ErasedJacobian. The concrete workspace stays behind the
# fixed-signature callback, so a new material does not recompile the integrator.
struct ErasedPressureNorm <: Function
    callback::Base.RefValue{_PressureNormWrapper}
end
@inline (norm::ErasedPressureNorm)(x::Vector{Float64}, t) = norm.callback[](x, Float64(t))
@inline (norm::ErasedPressureNorm)(x, t) = _rms_error(x, t)

struct PressureErrorNorm{C,G}
    cp::C
    state::Vector{PressureErrorDual}
    positions::Vector{PressureErrorDual}
    geometry::G
    gas::Vector{Int}
    masses::Vector{Float64}
    callback::Base.RefValue{_PressureNormWrapper}
    evaluations::Base.RefValue{Int}
end
function PressureErrorNorm(cp)
    d=cp.def
    ref=d.setup.reference
    u=PressureErrorDual.(cp.u0)
    x=PressureErrorDual.(ref.positions)
    g=MeshGeometry{PressureErrorDual}(d.setup.topology)
    # The components differ in type, so the per-cell loop reads these plain
    # vectors: indexing the tuple at run time would dispatch at every cell.
    comps=d.setup.material.components
    gas=collect(Int, gas_component_indices(comps))
    masses=Float64[molar_mass(comps[j]) for j in gas]
    return PressureErrorNorm(
        cp,
        u,
        x,
        g,
        gas,
        masses,
        Ref(_PressureNormWrapper(_rms_error)),
        Ref(0),
    )
end

# The norm reads only the cell volumes. On a line each is one expression of its
# two nodes, so it skips the faces and quadrature of the full geometry; in more
# dimensions the volumes rest on the face sweeps.
@inline _proposed_geometry!(g::MeshGeometry{1}, x, A, topo) = g
@inline _proposed_geometry!(g, x, A, topo) = geometry!(g, x, A, topo)
@inline function _proposed_volume(::MeshGeometry{1}, x, A, topo, i)
    n0, n1=topo.cells[i].nodes
    return shell_volume(topo.coordinates, x[n0], x[n1], A)
end
@inline _proposed_volume(g, x, A, topo, i) = g.cell_volume[i]

# The volume the condensed components hold, by recursion over the tuple, whose
# entries differ in type.
@inline _condensed_volume(::Tuple{}, m, T, j, V) = V
@inline function _condensed_volume(cs::Tuple, m, T, j, V)
    comp=first(cs)
    S=is_gas(comp) ? V : V+m[j]/get_intrinsic_density(comp, T)
    return _condensed_volume(Base.tail(cs), m, T, j+1, S)
end

function bind_pressure_norm!(norm::PressureErrorNorm, integrator)
    has_gas(norm.cp.def.setup.material) || return norm
    hasproperty(integrator.cache, :atmp) || throw(
        ArgumentError(
            "the selected integrator does not expose its normalized local error vector for pressure-aware error control",
        ),
    )
    norm.callback[]=_PressureNormWrapper(
        (x, t)->begin
            base=_rms_error(x, t)
            x===integrator.cache.atmp || return base
            norm.evaluations[]+=1
            return max(
                base,
                pressure_error(
                    norm,
                    integrator.u,
                    integrator.uprev,
                    x,
                    integrator.opts.abstol,
                    integrator.opts.reltol,
                ),
            )
        end,
    )
    return norm
end

# The pore volume is a difference of node positions, resolved no better than
# this many units in the last place of the positions; the pressure of a thin,
# nearly closed cell is known no better than that, whatever the budget asks.
const PORE_VOLUME_ROUNDOFF = 16

# The share of the pore volume `W` of cell `i` that round-off leaves unresolved:
# zero where no node moves, on a line a share of what the cell's faces sweep.
@inline _pore_roundoff(moving::Bool, x, A, topo, i, W) = 0.0
@inline function _pore_roundoff(moving::Bool, x, A, topo::MeshTopology{1}, i, W)
    moving || return 0.0
    n0, n1=topo.cells[i].nodes
    r=max(abs(ForwardDiff.value(x[n0])), abs(ForwardDiff.value(x[n1])))
    area=face_area(topo.coordinates, r, ForwardDiff.value(A))
    return PORE_VOLUME_ROUNDOFF*eps(r)*area/W
end

"""
    pressure_error(norm,u,uprev,scaled_error,abstol,reltol)

Largest local pressure-error ratio using the actual proposed pore volume and
thermal state. The error direction reconstructs the stepper's unscaled error;
forward sensitivities propagate through geometry and the selected caloric map.
The absolute budget is divided among gas species; their absolute pressure
errors also share one total budget, and neither budget lies below the pressure
the round-off of the cell's pore volume leaves unresolved. This controls a local
estimate, not accumulated global solution error.
"""
function pressure_error(norm::PressureErrorNorm, u, uprev, scaled_error, abstol, reltol)
    cp=norm.cp
    def=cp.def
    setup=def.setup
    layout=def.layout
    mat=setup.material
    ref=setup.reference
    NC=n_components(mat)
    ng=length(norm.gas)
    ng==0 && return 0.0
    ud=norm.state
    @inbounds for j in eachindex(u)
        e=scaled_error[j]*(abstol[j]+reltol*max(abs(u[j]), abs(uprev[j])))
        ud[j]=ForwardDiff.Dual{PressureErrorTag}(u[j], e)
    end
    D=dimension(layout)
    @inbounds for nd = 1:n_nodes(layout), k = 1:D
        j=D*(nd-1)+k
        norm.positions[j]=is_moving(layout) ? ud[position_index(layout, nd, k)] :
                          PressureErrorDual(ref.positions[j])
    end
    A=has_chi(layout) ? ud[area_index(layout)] : PressureErrorDual(ref.cross_section)
    _proposed_geometry!(norm.geometry, norm.positions, A, setup.topology)
    worst=0.0
    @inbounds for i = 1:n_cells(layout)
        layout.mask[energy_index(layout, i)] || continue
        m=ntuple(j->ud[mass_index(layout, i, j)], Val(NC))
        T=if energy_state(layout) isa TemperatureState
            ud[energy_index(layout, i)]
        else
            first(
                invert_caloric(
                    setup.caloric,
                    setup.laws,
                    ud[energy_index(layout, i)],
                    m,
                    i,
                    cp.ws.derived.T_warm[i],
                ),
            )
        end
        V=_proposed_volume(norm.geometry, norm.positions, A, setup.topology, i)
        W=V-_condensed_volume(mat.components, m, T, 1, zero(T))
        ForwardDiff.value(W)>0 && ForwardDiff.value(T)>0 || return Inf
        sumerror=0.0
        P=0.0
        unresolved=_pore_roundoff(
            is_moving(layout),
            norm.positions,
            A,
            setup.topology,
            i,
            ForwardDiff.value(W),
        )
        for k = 1:ng
            p=m[norm.gas[k]]*R_GAS*T/(norm.masses[k]*W)
            ep=abs(ForwardDiff.partials(p)[1])
            pj=max(ForwardDiff.value(p), 0.0)
            isfinite(ep) || return Inf
            worst=max(worst, ep/max(ref.pressure_tolerance/ng+reltol*pj, unresolved*pj))
            sumerror+=ep
            P+=pj
        end
        worst=max(worst, sumerror/max(ref.pressure_tolerance+reltol*P, unresolved*P))
    end
    return worst
end
