# Coupled surface species/thermal evaluation. q is positive into the material;
# all species and assembled energy fluxes are positive outward.
struct GasBoundaryBinding{B,I,M,N}
    bc::B
    indices::I
    masses::M
    closed::Bool
    names::N
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
    )
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

@inline _mass_film(::NeumannPressureBC) = nothing
@inline _mass_film(bc::Union{AmbientGasBC,ConvectivePressureBC}) = bc.film
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
    film=_mass_film(bc)
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
    xs=ntuple(j->begin
        k=findfirst(==(j), bound.indices)
        k===nothing ? zero(Ts+sum(N)) : x[k]+zero(Ts+sum(N))
    end, Val(NC))
    face=BlowingFace(
        J,
        Tt,
        material;
        surface_mole_fractions = film===nothing ? nothing : xs,
    )
    qsurf=_surface_heat(thermal, Ts, data, t, rad, face)
    h=component_enthalpies(laws, Ts)
    hin=component_enthalpies(laws, data.Tin)
    Hs=sum(j->J[j]*h[j], 1:NC)
    Hr=sum(j->J[j]*(primal(J[j])>=0 ? h[j] : hin[j]), 1:NC)
    # Convert the thermal wall-film balance to the direct-reservoir energy
    # control surface. A prescribed flux or adiabatic boundary has no such film.
    qres=qsurf+(_has_thermal_film(thermal) ? Hr-Hs : zero(Hr))
    g=film===nothing ? zero(Ts) : evaluate_param_T(film.conductance, Tt)
    g>=0 && !isnan(g) ||
        throw(ArgumentError("mass-film conductance must remain nonnegative"))
    return (; Ts, N, J, xs, Ys, ps, rawPs, v = half.v, qsurf, qres, Hs, Hr, g, face)
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
    film=_mass_film(bc)
    nf=film===nothing ? 0 : NG-1
    hydraulic=_hydraulic(bc)
    s=_surface_data(z, data, bound, thermal, material, laws, t, rad, pivot, pscale, Tscale)
    m=sum(s.N)
    F=film===nothing || isinf(primal(s.g)) ? zero(s.N) :
      mass_film_flux(s.Ys, data.Y∞, m, s.g)
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
                R=evaluate_param(bc.resistance, t)
                isfinite(R) && R>=0 || throw(
                    ArgumentError(
                        "hydraulic resistance must remain finite and nonnegative",
                    ),
                )
                flow=bc.flow===:darcy_velocity ? s.v :
                     bc.flow===:mass ? m : sum(s.N ./ bound.masses)
                return (s.rawPs-data.Pr-R*flow)/pscale
            elseif thermal isa DirichletThermalBC
                return (s.Ts-evaluate_param(thermal.T, t))/Tscale
            else
                return (data.k*(s.Ts-data.Tcell)/data.d-s.qsurf)/heatscale
            end
        end,
    )
end

function coupled_gas_boundary(
    bound::GasBoundaryBinding{B,I,<:StaticVector{NG}},
    thermal,
    data,
    material,
    laws,
    t,
    rad,
) where {B,I,NG}
    bc=bound.bc
    film=_mass_film(bc)
    nf=film===nothing ? 0 : NG-1
    hydraulic=_hydraulic(bc)
    nunknown=nf+Int(hydraulic)+1
    data.Pr>0 && _stream_temperature_check(thermal, data.Tin, t)
    film!==nothing &&
        data.Pr==0 &&
        throw(ArgumentError("an isobaric mass film requires a positive reservoir pressure"))
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
    s=_surface_data(z, data, bound, thermal, material, laws, t, rad, pivot, pscale, Tscale)
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
