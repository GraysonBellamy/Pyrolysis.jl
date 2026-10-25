using Test, StaticArrays, LinearAlgebra, Random, ForwardDiff
using Pyrolysis.Physics: interface_flux, pore_half_flux
using Pyrolysis.Materials: R_GAS

@testset "total-flux interface and mechanism switches" begin
    rng=MersenneTwister(932)
    for trial = 1:500
        pl=SVector{3}(1 .+ 10rand(rng, 3))
        pr=SVector{3}(1 .+ 10rand(rng, 3))
        Tl, Tr=300+600rand(rng), 300+600rand(rng)
        dl, dr=0.01+rand(rng), 0.01+rand(rng)
        λl, λr=rand(rng), rand(rng)
        ml, mr=rand(rng), rand(rng)
        trial%3==0 && (λl = 0.0; mr = 0.0)
        trial%3==1 && (λr = 0.0; ml = 0.0)
        cl=(K = SVector(λl, λl, λl), t = 0.0, mob = ml)
        cr=(K = SVector(λr, λr, λr), t = 0.0, mob = mr)
        f=interface_flux(pl, pr, Tl, Tr, dl, dr, cl, cr)
        Ts=(dr*Tl+dl*Tr)/(dl+dr)
        nl=pore_half_flux(pl, f.p, Tl, Ts, Tl, dl, cl.K, 0.0, ml).n
        nr=pore_half_flux(f.p, pr, Ts, Tr, Tr, dr, cr.K, 0.0, mr).n
        @test f.n ≈ nl rtol=2e-12 atol=2e-15
        @test f.n ≈ nr rtol=2e-12 atol=2e-15
        @test minimum(f.p)>=0
    end
    # Piecewise-constant nonisothermal wall resistance, independent quadrature.
    l, r=SVector(2.0, 3.0), SVector(1.0, 5.0)
    cl=(K = SVector(2.0, 2.0), t = 0.0, mob = 0.0)
    cr=(K = SVector(7.0, 7.0), t = 0.0, mob = 0.0)
    @test interface_flux(l, r, 300.0, 900.0, 0.3, 0.7, cl, cr).n ≈
          (l-r)/(R_GAS*(0.3*300/2+0.7*900/7)) rtol=2e-15
    # Exact mechanism switch is conducting. Pure gas donor-half discretization:
    # 2(2-P*) = 2P*(P*-1), P*=sqrt(2); continuum root differs and is tested on refinement.
    cl=(K = SVector(1.0), t = 0.0, mob = 0.0)
    cr=(K = SVector(0.0), t = 0.0, mob = 1.0)
    f=interface_flux(SVector(2.0), SVector(1.0), 1/R_GAS, 1/R_GAS, 0.5, 0.5, cl, cr)
    @test f.p[1]≈sqrt(2.0)
    @test f.n[1]≈2*(2-sqrt(2.0))
    @test f.n[1]>0
end

@testset "dusty interface root and implicit derivative" begin
    rng=MersenneTwister(933)
    for trial = 1:200
        pl=SVector{3}(1 .+ 10rand(rng, 3))
        pr=SVector{3}(1 .+ 10rand(rng, 3))
        cl=(K = SVector{3}(0.1 .+ rand(rng, 3)), t = 10rand(rng), mob = 0.1rand(rng))
        cr=(K = SVector{3}(0.1 .+ rand(rng, 3)), t = 10rand(rng), mob = 0.1rand(rng))
        f=interface_flux(pl, pr, 500.0, 500.0, 0.3, 0.7, cl, cr)
        nl=pore_half_flux(pl, f.p, 500.0, 500.0, 500.0, 0.3, cl.K, cl.t, cl.mob).n
        nr=pore_half_flux(f.p, pr, 500.0, 500.0, 500.0, 0.7, cr.K, cr.t, cr.mob).n
        @test nl≈nr rtol=2e-9 atol=2e-13
        @test sum(f.n .* log.(pr ./ pl))≤1e-13*sum(abs, f.n .* log.(pr ./ pl))
        @test minimum(f.p)>=0
        @test f.residual≤2e-10
        if trial≤10
            fun=z->interface_flux(z, pr, 500.0, 500.0, 0.3, 0.7, cl, cr).n
            J=ForwardDiff.jacobian(fun, pl)
            ε=1e-4
            fd=hcat(
                [
                    (fun(setindex(pl, pl[j]+ε, j))-fun(setindex(pl, pl[j]-ε, j)))/(2ε) for
                    j = 1:3
                ]...,
            )
            @test J≈fd rtol=2e-6 atol=1e-10
        end
    end
end
