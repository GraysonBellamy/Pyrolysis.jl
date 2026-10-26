using Test, Pyrolysis, StaticArrays, LinearAlgebra, ForwardDiff, Random
using Pyrolysis.Physics: interface_flux, TransportFailure, bounded_transport_root
using Pyrolysis.Materials: R_GAS

# Independent full friction matrices. Every dependent quantity, including the
# face pressure, forcing and donor composition, is formed in high precision.
function reference_half(left, right, Tl, Tr, Tpath, d, K, t, mob)
    Pl, Pr=sum(left), sum(right)
    P=(Pl+Pr)/2
    g=(left-right)/(BigFloat(R_GAS)*Tpath*d)
    a=K ./ (1 .+ t*P .* K)
    positive=sum(a .* g)>=0
    donor=positive ? left : right
    Pup=positive ? Pl : Pr
    p=P .* donor ./ Pup
    A=[(i==j ? inv(K[i])+t*P : zero(P))-t*p[i] for i in eachindex(K), j in eachindex(K)]
    v=mob*(Pl-Pr)/d
    return A\g+v .* (v>=0 ? left/(BigFloat(R_GAS)*Tl) : right/(BigFloat(R_GAS)*Tr))
end

function reference_interface(left, right, Tl, Tr, dl, dr, cl, cr)
    setprecision(256) do
        l, r=BigFloat.(left), BigFloat.(right)
        tl, tr, dL, dR=BigFloat(Tl), BigFloat(Tr), BigFloat(dl), BigFloat(dr)
        kl, kr=BigFloat.(cl.K), BigFloat.(cr.K)
        al, ar, ml, mr=BigFloat(cl.t), BigFloat(cr.t), BigFloat(cl.mob), BigFloat(cr.mob)
        ts=(dR*tl+dL*tr)/(dL+dR)
        leftflux=p->reference_half(l, p, tl, ts, tl, dL, kl, al, ml)
        rightflux=p->reference_half(p, r, ts, tr, tr, dR, kr, ar, mr)
        f=p->leftflux(p)-rightflux(p)
        p=collect((dR*l+dL*r)/(dL+dR))
        for iteration = 1:80
            residual=f(p)
            if norm(residual, Inf)<big"1e-68"
                return (p = copy(p), n = (leftflux(p)+rightflux(p))/2)
            end
            step=ForwardDiff.jacobian(f, p)\residual
            α=one(BigFloat)
            while minimum(p-α*step)<0 || norm(f(p-α*step), Inf)>norm(residual, Inf)
                α/=2
                α>big"1e-40" || error("independent root line search failed")
            end
            p-=α*step
        end
        error("independent root failed to converge")
    end
end

@testset "qualified interface roots at strong counterdiffusion" begin
    rng=MersenneTwister(1951)
    for power in (-12, 0, 4, 8, 12, 16), trial = 1:4
        # Exact binary endpoint sums; no pre-rounded dependent reference data.
        l=SVector(0.75, 0.25)
        r=SVector(0.25, 0.75)
        cl=(K = SVector(0.25+rand(rng), 0.25+rand(rng)), t = 10.0^power, mob = 0.0)
        cr=(K = SVector(0.25+rand(rng), 0.25+rand(rng)), t = 10.0^power, mob = 0.0)
        reference=reference_interface(l, r, 500.0, 500.0, 0.3, 0.7, cl, cr)
        @test minimum(reference.p)>0
        result=try
            interface_flux(l, r, 500.0, 500.0, 0.3, 0.7, cl, cr)
        catch failure
            failure isa TransportFailure || rethrow()
            failure
        end
        if result isa TransportFailure
            # Once pressure corrections are below Float64 resolution, the
            # physical reference still exists. The production solver must
            # diagnose its numerical limitation instead of accepting a flux.
            @test power>=8
            @test result.reason in (
                :interface_line_search,
                :interface_iterations,
                :singular_interface,
                :ill_conditioned_local_root,
                :unresolved_local_root,
                :unresolved_flux_continuity,
            )
        else
            scale=maximum(abs, reference.n)
            @test maximum(abs, BigFloat.(result.n)-reference.n)<=2e-8scale
            @test maximum(abs, BigFloat.(result.p)-reference.p)<=2e-8
        end
    end
end

@testset "local root diagnostics and exact initial-root derivatives" begin
    @test_throws TransportFailure bounded_transport_root(x->zero(x), SVector(1.0, 1.0))
    @test_throws TransportFailure bounded_transport_root(x->SVector(x[1]^2+1), SVector(1.0))
    f=x->SVector(2x[1]-4)
    root=bounded_transport_root(f, SVector(2.0))
    @test root.x==SVector(2.0)
    @test root.iterations==0
    derivative=ForwardDiff.derivative(3.0) do parameter
        only(Pyrolysis.Physics.implicit_transport_root(x->SVector(2x[1]-parameter-1), root))
    end
    @test derivative==0.5
end

@testset "equilibrium and resolved ulp-sized interface gradients" begin
    for model in (:lumped, :dusty),
        P in (8e4, 1e5, 1.2e5),
        distances in ((0.0025, 0.0025), (0.0013, 0.0047))

        l=SVector(0.7P, 0.3P)
        cl=model===:dusty ? (K = SVector(5.1e-5, 4.1e-5), t = 0.8, mob = 0.0) :
           (K = SVector(1e-5, 1e-5), t = 0.0, mob = 0.0)
        cr=model===:dusty ? (K = SVector(2.4e-5, 7.8e-5), t = 1.2, mob = 0.0) : cl
        dl, dr=distances
        equilibrium=interface_flux(l, l, 500.0, 500.0, dl, dr, cl, cr)
        @test iszero(equilibrium.n)
        @test equilibrium.p==l
        # A pressure difference close to one ulp must not be hidden by a
        # relative pressure or composition floor. The reference constructs
        # every dependent quantity from these inputs in 256-bit arithmetic.
        r=SVector(prevfloat(l[1]), nextfloat(l[2]))
        ref=reference_interface(l, r, 500.0, 500.0, dl, dr, cl, cr)
        result=interface_flux(l, r, 500.0, 500.0, dl, dr, cl, cr)
        @test maximum(abs, BigFloat.(result.n)-ref.n)<2e-8maximum(abs, ref.n)
        J=ForwardDiff.jacobian(p->interface_flux(p, l, 500.0, 500.0, dl, dr, cl, cr).n, l)
        # Linearize the independent high-precision interface solve with a
        # symmetric perturbation whose quadratic remainder is negligible.
        Jref=setprecision(256) do
            h=big"1e-15"
            columns=map(1:2) do j
                delta=SVector(j==1 ? h : zero(h), j==2 ? h : zero(h))
                up=reference_interface(
                    BigFloat.(l)+delta,
                    BigFloat.(l),
                    500.0,
                    500.0,
                    dl,
                    dr,
                    cl,
                    cr,
                ).n
                down=reference_interface(
                    BigFloat.(l)-delta,
                    BigFloat.(l),
                    500.0,
                    500.0,
                    dl,
                    dr,
                    cl,
                    cr,
                ).n
                (up-down)/(2h)
            end
            hcat(columns...)
        end
        @test J≈Float64.(Jref) rtol=2e-8
    end
end

@testset "trace-species interface scaling at an ulp pressure jump" begin
    # Original face inputs from an atmospheric-air heating problem. A tiny
    # product concentration is equal on both sides; its nonzero driven flux
    # must survive without a composition floor or a false singularity.
    l=SVector(5.416273477600149e-23, prevfloat(101325.0))
    r=SVector(l[1], 101325.0)
    c=(
        K = SVector(4.221617965530974e-5, 5.2027267845421225e-5),
        t = 1.480384900074019,
        mob = 5.417363866595937e-11,
    )
    d=0.0004999999999999996
    ref=reference_interface(l, r, 300.0, 300.0, d, d, c, c)
    actual=interface_flux(l, r, 300.0, 300.0, d, d, c, c)
    @test actual.n[1]!=0
    @test all(abs.(BigFloat.(actual.n) .- ref.n) .< 2e-8 .* abs.(ref.n))
    J=ForwardDiff.jacobian(p->interface_flux(p, r, 300.0, 300.0, d, d, c, c).n, l)
    reference=setprecision(256) do
        columns=map(1:2) do j
            h=min(BigFloat(l[j])*big"1e-8", big"1e-20")
            delta=BigFloat[j==k ? h : 0 for k = 1:2]
            a=reference_interface(
                BigFloat.(l)+delta,
                BigFloat.(r),
                300.0,
                300.0,
                d,
                d,
                c,
                c,
            ).n
            b=reference_interface(
                BigFloat.(l)-delta,
                BigFloat.(r),
                300.0,
                300.0,
                d,
                d,
                c,
                c,
            ).n
            (a-b)/(2h)
        end
        hcat(columns...)
    end
    @test all(abs.(BigFloat.(J) .- reference) .< 2e-8 .* abs.(reference))
end
