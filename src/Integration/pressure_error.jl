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
    callback::Base.RefValue{_PressureNormWrapper}
    evaluations::Base.RefValue{Int}
end
function PressureErrorNorm(cp)
    d=cp.def
    ref=d.setup.reference
    u=PressureErrorDual.(cp.u0)
    x=PressureErrorDual.(ref.positions)
    g=MeshGeometry{PressureErrorDual}(d.setup.topology)
    return PressureErrorNorm(cp, u, x, g, Ref(_PressureNormWrapper(_rms_error)), Ref(0))
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

"""
    pressure_error(norm,u,uprev,scaled_error,abstol,reltol)

Largest local pressure-error ratio using the actual proposed pore volume and
thermal state. The error direction reconstructs the stepper's unscaled error;
forward sensitivities propagate through geometry and the selected caloric map.
The absolute budget is divided among gas species; their absolute pressure
errors also share one total budget. This controls a local estimate, not
accumulated global solution error.
"""
function pressure_error(norm::PressureErrorNorm, u, uprev, scaled_error, abstol, reltol)
    cp=norm.cp
    def=cp.def
    setup=def.setup
    layout=def.layout
    mat=setup.material
    ref=setup.reference
    NC=n_components(mat)
    ng=count(is_gas, mat.components)
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
    geometry!(norm.geometry, norm.positions, A, setup.topology)
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
        V=norm.geometry.cell_volume[i]
        W=V-sum(
            ntuple(
                j->is_gas(mat.components[j]) ? zero(T) :
                   m[j]/get_intrinsic_density(mat.components[j], T),
                Val(NC),
            ),
        )
        ForwardDiff.value(W)>0 && ForwardDiff.value(T)>0 || return Inf
        sumerror=0.0
        P=0.0
        for j = 1:NC
            comp=mat.components[j]
            is_gas(comp) || continue
            p=m[j]*R_GAS*T/(molar_mass(comp)*W)
            ep=abs(ForwardDiff.partials(p)[1])
            pj=max(ForwardDiff.value(p), 0.0)
            isfinite(ep) || return Inf
            worst=max(worst, ep/(ref.pressure_tolerance/ng+reltol*pj))
            sumerror+=ep
            P+=pj
        end
        worst=max(worst, sumerror/(ref.pressure_tolerance+reltol*P))
    end
    return worst
end
