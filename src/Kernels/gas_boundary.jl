# Coupled surface species/thermal evaluation. q is positive into the material;
# all species and assembled energy fluxes are positive outward.
struct GasBoundaryBinding{B,I,M,N,F}
    bc::B
    indices::I
    masses::M
    closed::Bool
    names::N
    film::F
end
function GasBoundaryBinding(bc, material)
    idx=gas_component_indices(material.components)
    M=SVector{length(idx)}(map(j->molar_mass(material.components[j]), idx))
    return GasBoundaryBinding(
        bc,
        idx,
        M,
        closed_transport_channel(material),
        map(c->c.name, material.components),
        _binding_film(bc),
    )
end

# The film a face declares, resolved once at binding: an isobaric film has no
# meaning against a vacuum, so an `AnalogyFilm` on a constant-vacuum reservoir
# leaves the reservoir in direct contact. Resolving here keeps the surface
# solve concretely typed.
@inline _binding_film(::AbstractMassBC) = nothing
@inline function _binding_film(bc::Union{AmbientGasBC,ConvectivePressureBC})
    P=bc.reservoir.P
    return bc.film isa AnalogyFilm && P isa Real && P==0 ? nothing : bc.film
end

@inline _has_thermal_film(::AbstractThermalBC) = false
@inline _has_thermal_film(::ConvectiveBC) = true
@inline _has_thermal_film(bc::CombinedThermalBC) = any(_has_thermal_film, bc.bcs)
@inline _stream_temperature_check(::AbstractThermalBC, T, t) = nothing
@inline function _stream_temperature_check(bc::ConvectiveBC, T, t)
    isapprox(evaluate_param(bc.T_inf, t), T; rtol = 1e-12, atol = 1e-10) || throw(
        ArgumentError(
            "thermal film and mass reservoir describe the same gas stream and require the same incoming temperature",
        ),
    )
    return nothing
end
@inline function _stream_temperature_check(bc::CombinedThermalBC, T, t)
    foreach(b->_stream_temperature_check(b, T, t), bc.bcs)
    return nothing
end
@inline _primal_boundary(x) = primal(x)
@inline _primal_boundary(x::Tuple) = map(_primal_boundary, x)
@inline _primal_boundary(x::NamedTuple) = map(_primal_boundary, x)
@inline _primal_boundary(x::StaticVector) = map(primal, x)

# The gas film a `PoolingFilm` states for its face: the molar film of one
# conductance ψ g s_v^{n_D}, with ψ from the solved surface composition and the
# incoming temperature of the convective condition that carries the model.
struct PoolingGasFilm{C<:ConvectiveBC}
    convective::C
end

# The convective member carrying a `PoolingFilm`, `nothing` without one, found
# through the concrete member types so the resolved film stays inferred.
@inline _pooling_member(::AbstractThermalBC) = nothing
@inline _pooling_member(bc::ConvectiveBC{<:Any,<:Any,<:PoolingFilm}) = bc
@inline _pooling_member(bc::CombinedThermalBC) = _first_pooling_member(bc.bcs...)
@inline _first_pooling_member() = nothing
@inline function _first_pooling_member(b, rest...)
    p=_pooling_member(b)
    return p===nothing ? _first_pooling_member(rest...) : p
end
@inline _n_convective(::AbstractThermalBC) = 0
@inline _n_convective(::ConvectiveBC) = 1
@inline _n_convective(bc::CombinedThermalBC) = sum(_n_convective, bc.bcs)

# The film a face carries. An `AnalogyFilm` is the twin of the convective heat
# film: without one the reservoir is in direct contact with the pore gas, and
# where the convective film is a `PoolingFilm` the twin is that model's film.
# An explicit conductance on such a face would leave its gas exchange without
# the pooling factor its heat exchange applies, so it is rejected.
@inline _resolve_film(film, thermal) = film
@inline function _resolve_film(film::MassTransferFilm, thermal)
    _pooling_member(thermal)===nothing || throw(
        ArgumentError(
            "a face whose convective condition carries a PoolingFilm takes its gas exchange from that model; give film = AnalogyFilm() (the default), not an explicit MassTransferFilm",
        ),
    )
    return film
end
@inline function _resolve_film(film::AnalogyFilm, thermal)
    _has_thermal_film(thermal) || return nothing
    member=_pooling_member(thermal)
    member===nothing && return film
    film.Le==1 || throw(
        ArgumentError(
            "PoolingFilm sets the gas conductance of its face through g and n_D; AnalogyFilm(Le = $(film.Le)) has no role there",
        ),
    )
    _n_convective(thermal)==1 || throw(
        ArgumentError(
            "a face whose convective condition carries a PoolingFilm must have no other convective condition",
        ),
    )
    return PoolingGasFilm(member)
end
# Bare convective coefficient of a face: the sum over its convective members.
@inline _bare_convective_h(::AbstractThermalBC, t) = 0.0
@inline _bare_convective_h(bc::ConvectiveBC, t) = evaluate_param(bc.h, t)
@inline _bare_convective_h(bc::CombinedThermalBC, t) =
    sum(b->_bare_convective_h(b, t), bc.bcs)
@inline _film_conductance(
    film::MassTransferFilm,
    Tt,
    Ys,
    data,
    bound,
    material,
    thermal,
    t,
) = evaluate_param_T(film.conductance, Tt)
# Heat–mass transfer analogy at the film temperature: g = h/(c_p,∞ Le^{2/3})
# with c_p,∞ the mass-weighted heat capacity of the reservoir gas.
@inline function _film_conductance(
    film::AnalogyFilm,
    Tt,
    Ys,
    data,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    material,
    thermal,
    t,
) where {B,I,NG}
    Tf=(Tt+data.Tin)/2
    props=BoundaryConditions._gas_props(material, Tf)
    cp=zero(Tf)
    for k = 1:NG
        cp+=data.Y∞[k]*props[bound.indices[k]].cp
    end
    primal(cp)>0 ||
        throw(ArgumentError("the analogy film needs a reservoir gas with a heat capacity"))
    return _bare_convective_h(thermal, t)/(cp*film.Le^(2/3))
end
# The pooling model's conductance at the surface composition of this iterate,
# so the gas exchange carries the ψ the heat exchange applies.
@inline _film_conductance(film::PoolingGasFilm, Tt, Ys, data, bound, material, thermal, t) =
    BoundaryConditions.pooling_film_conductance(
        film.convective.blowing,
        _component_fractions(Ys, bound, material),
        material,
        Tt,
        evaluate_param(film.convective.T_inf, t),
    )

# Gas-component values spread over every material component (zero for the
# condensed ones).
@inline _component_fractions(Y, bound, material::Material{NC}) where {NC} =
    ntuple(Val(NC)) do j
        k=findfirst(==(j), bound.indices)
        k===nothing ? zero(eltype(Y)) : Y[k]
    end

# Integrated film flux of every gas component: the mass-fraction film of one
# conductance, or for a pooling film the molar film of one conductance, with
# moles counted as the model counts them.
@inline _film_flux(film, s, data, bound, material, m) =
    mass_film_flux(s.Ys, data.Y∞, m, s.g)
@inline function _film_flux(
    film::PoolingGasFilm,
    s,
    data,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    material,
    m,
) where {B,I,NG}
    Mm=BoundaryConditions.film_molar_masses(film.convective.blowing, material)
    M=SVector{NG}(ntuple(k->Mm[bound.indices[k]], Val(NG)))
    ns=s.Ys ./ M
    n∞=data.Y∞ ./ M
    # G = g/M_∞ with 1/M_∞ = Σ Y_∞/M.
    return molar_film_flux(ns/sum(ns), n∞/sum(n∞), sum(s.N ./ M), s.g*sum(n∞), M)
end

# The empirical blowing closures take the incoming gas to be their free
# stream. On a solved film that is a condition on the reservoir: one of the
# free-stream gas alone admits no product ingress at the root.
@inline _check_blowing_reservoir(::AbstractThermalBC, bound, x∞) = nothing
@inline _check_blowing_reservoir(bc::CombinedThermalBC, bound, x∞) =
    foreach(b->_check_blowing_reservoir(b, bound, x∞), bc.bcs)
@inline function _check_blowing_reservoir(
    bc::ConvectiveBC{<:Any,<:Any,<:Union{PropertyRatioFilm,PoolingFilm}},
    bound,
    x∞,
)
    fs=bc.blowing.free_stream.name
    for k in eachindex(x∞)
        name=bound.names[bound.indices[k]]
        primal(x∞[k])>0 &&
            name!==fs &&
            throw(
                ArgumentError(
                    "$(nameof(typeof(bc.blowing))) on a solved mass film needs a reservoir of its free-stream gas :$fs alone; it holds :$name",
                ),
            )
    end
    return nothing
end
@inline _hydraulic(::AbstractMassBC) = false
@inline _hydraulic(::ConvectivePressureBC) = true
@inline _boundary_reservoir(bc::Union{AmbientGasBC,ConvectivePressureBC}, material, t) =
    BoundaryConditions._reservoir_values(bc.reservoir, material, t)
@inline _boundary_reservoir(bc::NeumannPressureBC, material, t) =
    BoundaryConditions._reservoir_values(
        (P = 1.0, T = bc.T, mole_fractions = bc.mole_fractions),
        material,
        t,
    )
@inline _surface_pressure(bc::Union{AmbientGasBC,ConvectivePressureBC}, data, t) = data.Pr
@inline _surface_pressure(bc::NeumannPressureBC, data, t) =
    sum(data.pl)+data.d*evaluate_param(bc.gradient, t)

@inline function _surface_data(
    film,
    z,
    data,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    thermal,
    material::Material,
    laws,
    t,
    rad,
    pivot,
    pscale,
    Tscale,
) where {B,I,NG}
    bc=bound.bc
    nf=film===nothing ? 0 : NG-1
    hydraulic=_hydraulic(bc)
    Ts=z[end]*Tscale
    Tt=trial_temperature(Ts)
    rawPs=hydraulic ? z[nf+1]*pscale : _surface_pressure(bc, data, t)
    # A prescribed reservoir is already a validated thermodynamic state.
    # Only the unknown hydraulic pressure needs a trial-state extension.
    Ps=hydraulic ? positive_part(rawPs, data.pbudget) : rawPs
    x, Ys = if film===nothing
        data.x∞, data.Y∞
    else
        fractions=SVector{NG}(
            ntuple(
                j->j==pivot ? 1-sum((z[k] for k = 1:nf); init = zero(z[end])) :
                   z[j<pivot ? j : j-1],
                Val(NG),
            ),
        )
        mole=fractions ./ bound.masses
        mole/sum(mole), fractions
    end
    ps=Ps .* x
    Tpath=(data.Ttrial+Tt)/2
    coeff=transport_coefficients(data.ξ, Tpath, material, data.φ)
    wall=coeff.K
    K=SVector{NG}(map(j->wall[j], bound.indices))
    mob=coeff.κ/material.gas_viscosity(Tpath)
    # A wall heat film places the gas in LTE at the physical surface.
    # Direct exchange without that film retains the incoming reservoir state.
    Tdonor=_has_thermal_film(thermal) ? Tt : data.Tin
    half=pore_half_flux(data.pl, ps, data.Ttrial, Tdonor, Tpath, data.d, K, coeff.t, mob)
    N=half.n .* bound.masses
    NC=n_components(material)
    J=ntuple(j->begin
        k=findfirst(==(j), bound.indices)
        k===nothing ? zero(Ts+sum(N)) : N[k]
    end, Val(NC))
    Yface=_component_fractions(Ys .+ zero(Ts+sum(N)), bound, material)
    face=BlowingFace(
        J,
        Tt,
        material;
        surface_mass_fractions = film===nothing ? nothing : Yface,
    )
    qsurf=_surface_heat(thermal, Ts, data, t, rad, face)
    h=component_enthalpies(laws, Ts)
    hin=component_enthalpies(laws, data.Tin)
    Hs=sum(j->J[j]*h[j], 1:NC)
    Hr=sum(j->J[j]*(primal(J[j])>=0 ? h[j] : hin[j]), 1:NC)
    # Convert the thermal wall-film balance to the direct-reservoir energy
    # control surface. A prescribed flux or adiabatic boundary has no such film.
    qres=qsurf+(_has_thermal_film(thermal) ? Hr-Hs : zero(Hr))
    g=film===nothing ? zero(Ts) :
      _film_conductance(film, Tt, Ys, data, bound, material, thermal, t)
    g>=0 && !isnan(g) ||
        throw(ArgumentError("mass-film conductance must remain nonnegative"))
    return (; Ts, N, J, Ys, ps, rawPs, v = half.v, qsurf, qres, Hs, Hr, g, face)
end

@inline _surface_heat(::AdiabaticBC, Ts, data, t, rad, face) = zero(Ts)
@inline _surface_heat(bc::HeatFluxBC, Ts, data, t, rad, face) =
    evaluate_param(bc.flux, t)+zero(Ts)
@inline _surface_heat(bc::DirichletThermalBC, Ts, data, t, rad, face) =
    data.k*(Ts-data.Tcell)/data.d
@inline function _surface_heat(bc::AbstractThermalBC, Ts, data, t, rad, face)
    q, _=evaluate_bc_flux_and_derivative(
        bc,
        trial_temperature(Ts),
        t,
        data.ε,
        rad,
        face;
        a_eff = data.a,
    )
    return q
end

@inline function _surface_residual(
    film,
    z,
    data,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    thermal,
    material,
    laws,
    t,
    rad,
    pivot,
    pscale,
    Tscale,
    massscale,
    heatscale,
) where {B,I,NG}
    bc=bound.bc
    nf=film===nothing ? 0 : NG-1
    hydraulic=_hydraulic(bc)
    s=_surface_data(
        film,
        z,
        data,
        bound,
        thermal,
        material,
        laws,
        t,
        rad,
        pivot,
        pscale,
        Tscale,
    )
    m=sum(s.N)
    F=film===nothing || isinf(primal(s.g)) ? zero(s.N) :
      _film_flux(film, s, data, bound, material, m)
    S=typeof(s.Ts+m+s.rawPs)
    return SVector{length(z),S}(
        ntuple(Val(length(z))) do k
            if k<=nf
                j=k<pivot ? k : k+1
                # With an identically closed pore channel and zero film conductance,
                # composition is unobservable. Pin it to the reservoir state.
                # A coefficient that merely vanishes at this trial is not closed.
                return (isinf(primal(s.g)) || (bound.closed && primal(s.g)==0)) ?
                       s.Ys[j]-data.Y∞[j]+zero(s.Ts+m+s.rawPs) : (s.N[j]-F[j])/massscale
            elseif hydraulic && k==nf+1
                R=evaluate_param_T(bc.resistance, trial_temperature(s.Ts))
                isfinite(primal(R)) && primal(R)>=0 || throw(
                    ArgumentError(
                        "hydraulic resistance must remain finite and nonnegative",
                    ),
                )
                flow=bc.flow===:darcy_velocity ? s.v : m
                return (s.rawPs-data.Pr-R*flow)/pscale
            elseif thermal isa DirichletThermalBC
                return (s.Ts-evaluate_param(thermal.T, t))/Tscale
            else
                return (data.k*(s.Ts-data.Tcell)/data.d-s.qsurf)/heatscale
            end
        end,
    )
end

coupled_gas_boundary(bound::GasBoundaryBinding, thermal, data, material, laws, t, rad) =
    _coupled_gas_boundary(
        _resolve_film(bound.film, thermal),
        bound,
        thermal,
        data,
        material,
        laws,
        t,
        rad,
    )

function _coupled_gas_boundary(
    film,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    thermal,
    data,
    material,
    laws,
    t,
    rad,
) where {B,I,NG}
    bc=bound.bc
    nf=film===nothing ? 0 : NG-1
    hydraulic=_hydraulic(bc)
    nunknown=nf+Int(hydraulic)+1
    data.Pr>0 && _stream_temperature_check(thermal, data.Tin, t)
    film===nothing || _check_blowing_reservoir(thermal, bound, data.x∞)
    film!==nothing &&
        data.Pr==0 &&
        throw(
            ArgumentError(
                "an isobaric mass film requires a positive reservoir pressure; give film = nothing for a vacuum",
            ),
        )
    pivot=argmax(data.Y∞)
    pscale=max(primal(sum(data.pl)), primal(data.Pr), data.pbudget)
    Tscale=max(primal(data.Ttrial), data.Tin, 1.0)
    z0=SVector{nunknown}(ntuple(Val(nunknown)) do k
        k<=nf && return Float64(data.Y∞[k<pivot ? k : k+1])
        hydraulic && k==nf+1 && return Float64(data.Pr/pscale)
        return Float64(primal(data.Tcell)/Tscale)
    end)
    plain=_primal_boundary(data)
    s0=_surface_data(
        film,
        z0,
        plain,
        bound,
        thermal,
        material,
        laws,
        t,
        rad,
        pivot,
        pscale,
        Tscale,
    )
    massscale=max(maximum(abs, s0.N), isfinite(s0.g) ? s0.g : 0.0, 1e-14)
    heatscale=max(
        plain.k*Tscale/plain.d,
        abs(s0.qsurf),
        1.0,
        64eps(Float64)*abs(plain.Tcell)*plain.k/plain.d/2e-12,
    )
    fp=z->_surface_residual(
        film,
        z,
        plain,
        bound,
        thermal,
        material,
        laws,
        t,
        rad,
        pivot,
        pscale,
        Tscale,
        massscale,
        heatscale,
    )
    lower=SVector{nunknown}(ntuple(k->k==nunknown ? -Inf : 0.0, Val(nunknown)))
    feasible=z->nf==0 || sum((z[k] for k = 1:nf); init = zero(z[end]))<=1
    root=bounded_transport_root(
        fp,
        z0;
        lower = lower,
        admissible = feasible,
        tolerance = 2e-12,
    )
    f=z->_surface_residual(
        film,
        z,
        data,
        bound,
        thermal,
        material,
        laws,
        t,
        rad,
        pivot,
        pscale,
        Tscale,
        massscale,
        heatscale,
    )
    z=implicit_transport_root(f, root)
    s=_surface_data(
        film,
        z,
        data,
        bound,
        thermal,
        material,
        laws,
        t,
        rad,
        pivot,
        pscale,
        Tscale,
    )
    return merge(s, (iterations = root.iterations, residual = root.residual))
end

function _exchange_patch!(
    c,
    thermal,
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    f,
    kernel,
    ds,
    topo,
    material,
    t,
) where {B<:Union{AmbientGasBC,NeumannPressureBC,ConvectivePressureBC},I,NG}
    i, d=_boundary_cell(kernel, ds, topo, f)
    pl, _=_cell_pore(ds, i, bound.indices)
    res=_boundary_reservoir(bound.bc, material, t)
    x∞=SVector{NG}(ntuple(k->res.x[bound.indices[k]], Val(NG)))
    weighted=x∞ .* bound.masses
    Y∞=sum(weighted)==0 ? zero(weighted) : weighted/sum(weighted)
    data=(
        pl = pl,
        ξ = ntuple(j->ds.ξ⁺[j, i], Val(n_components(material))),
        φ = ds.φ[i],
        Ttrial = ds.T_trial[i],
        Tcell = ds.T[i],
        d = d,
        k = ds.k[i],
        ε = ds.ε[i],
        a = ds.a[i],
        Pr = res.P,
        Tin = res.T,
        x∞ = x∞,
        Y∞ = Y∞,
        pbudget = kernel.pressure_tolerance,
    )
    s=coupled_gas_boundary(
        bound,
        thermal,
        data,
        material,
        kernel.laws,
        t,
        kernel.radiation_model,
    )
    σ=_outward(topo, f, i)
    for k = 1:NG
        j=bound.indices[k]
        c.F_mass[j, f]+=σ*s.N[k]
        c.p_surface[j, f]=s.ps[k]
    end
    c.F_energy[f]+=σ*(s.Hr-s.qres)
    c.T_face[f]=s.Ts
    c.darcy_velocity[f]=σ*s.v
    c.surface_residual[f]=s.residual
    c.surface_iterations[f]=s.iterations
    c.boundary_enthalpy[f]=s.Hr
    c.boundary_heat[f]=s.qres
    c.surface_pressure_raw[f]=s.rawPs
    if eltype(ds.T)===Float64
        st=kernel.states[f]
        st.T_surface=s.Ts
        st.q_applied=s.qsurf
        st.q_conducted=ds.k[i]*(s.Ts-ds.T[i])/d
        st.T_cell=ds.T[i]
        st.k_eff=ds.k[i]
        st.Δz_half=d
        st.converged=true
        st.iterations=s.iterations
        BoundaryConditions._record_blowing!(st, thermal, s.face, t)
    end
    return nothing
end

function _exchange_patch!(
    c,
    thermal,
    bound::GasBoundaryBinding{<:ImpermeableBC},
    f,
    kernel,
    ds,
    topo,
    material,
    t,
)
    q, Ts=boundary_energy_flux(thermal, f, kernel, c, ds, topo, material, t)
    c.F_energy[f]+=q
    c.T_face[f]=Ts
    c.boundary_heat[f]=-q
    return nothing
end

function _exchange_patch!(
    c::Contributions{NC},
    thermal,
    bound::GasBoundaryBinding{<:MassFluxBC},
    f,
    kernel,
    ds,
    topo,
    material,
    t,
) where {NC}
    bc=bound.bc
    i, _=_boundary_cell(kernel, ds, topo, f)
    σ=_outward(topo, f, i)
    J=ntuple(Val(NC)) do j
        v=evaluate_param(BoundaryConditions._named_value(bc.fluxes, bound.names[j], 0.0), t)
        isfinite(v) || throw(ArgumentError("prescribed species flux must remain finite"))
        v
    end
    for j = 1:NC
        c.F_mass[j, f]+=σ*J[j]
    end
    q, Ts=boundary_energy_flux(thermal, f, kernel, c, ds, topo, material, t)
    h=component_enthalpies(kernel.laws, Ts)
    Hs=sum(j->J[j]*h[j], 1:NC)
    Hr=zero(Hs)
    for j = 1:NC
        if primal(J[j])<0
            Tin=evaluate_param(
                BoundaryConditions._named_value(bc.inflow_temperature, bound.names[j], NaN),
                t,
            )
            isfinite(Tin) && Tin>0 ||
                throw(ArgumentError("inflow temperature must remain finite and positive"))
            _stream_temperature_check(thermal, Tin, t)
            Hr+=J[j]*component_enthalpies(kernel.laws, Tin)[j]
        else
            Hr+=J[j]*h[j]
        end
    end
    qres=-q+(_has_thermal_film(thermal) ? Hr-Hs : zero(Hr))
    c.F_energy[f]+=σ*(Hr-qres)
    c.T_face[f]=Ts
    c.boundary_heat[f]=qres
    c.boundary_enthalpy[f]=Hr
    return nothing
end

"""
	solved_surface_fractions(kernel::ThermalBoundaryKernel, tag, c, f, material)

The surface mass fractions the gas exchange solved at face `f` of patch `tag`,
one per material component (zero for condensed ones), as the face's blowing
model saw them; `nothing` where the face exchanges gas without a mass film, in
which case the model saw none either. Recovered from the surface partial
pressures `c.p_surface` with the material's molar masses, the basis they
were written in.
"""
function solved_surface_fractions(
    kernel::ThermalBoundaryKernel,
    tag::Symbol,
    c,
    f,
    material,
)
    kernel.exchange===nothing && return nothing
    i=findfirst(p->p.tag===tag, kernel.patches)
    i===nothing && return nothing
    bound=kernel.exchange[i]
    _resolve_film(bound.film, kernel.patches[i].bc)===nothing && return nothing
    w=[c.p_surface[j, f] for j in bound.indices] .* bound.masses
    return _component_fractions(w/sum(w), bound, material)
end
