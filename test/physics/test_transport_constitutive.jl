using Test, LinearAlgebra, Random, StaticArrays, ForwardDiff
using Pyrolysis.Physics: molar_wall_flux, pore_half_flux
using Pyrolysis.Materials: R_GAS

# Every dependent quantity is constructed after conversion of the ORIGINAL
# inputs. In particular, Float64 mole fractions and r=tP are never the oracle.
function dense_reference(K, t, p, g)
    setprecision(256) do
        kb, tb, pb, gb = BigFloat.(K), BigFloat(t), BigFloat.(p), BigFloat.(g)
        P = sum(pb)
        A = [
            ((j==i ? inv(kb[j])+tb*P : zero(tb))-tb*pb[j]) for
            j in eachindex(kb), i in eachindex(kb)
        ]
        A\gb
    end
end

@testset "conditioned rank-one friction solve" begin
    rng=MersenneTwister(714)
    for ng = 2:6, exponent in (-24, -16, -8, 0, 4, 8, 12, 16), trial = 1:20
        K=SVector{ng}(10.0 .^ (4rand(rng, ng) .- 2))
        p=SVector{ng}(10.0 .^ (4rand(rng, ng) .- 2))
        g=SVector{ng}(randn(rng, ng))
        t=10.0^exponent/(sum(p)*maximum(K))
        n=molar_wall_flux(K, t, p, g)
        ref=dense_reference(K, t, p, g)
        # Componentwise absolute budget also covers cancellation to zero.
        scale=maximum(abs, ref)
        @test all(j->abs(BigFloat(n[j])-ref[j]) ≤ 3e-13*(abs(ref[j])+scale), 1:ng)
        # Check sum and original friction equations with high-precision arithmetic.
        kb, pb, gb, nb=BigFloat.(K), BigFloat.(p), BigFloat.(g), BigFloat.(n)
        tb=BigFloat(t)
        P=sum(pb)
        residual=nb ./ kb+tb .* (P .* nb-pb .* sum(nb))-gb
        @test maximum(abs, residual) ≤
              2e-13*max(
            maximum(abs, gb),
            maximum(abs, nb ./ kb)+maximum(abs, tb .* P .* nb),
        )
        perm=SVector{ng}(randperm(rng, ng))
        @test molar_wall_flux(K[perm], t, p[perm], g[perm]) ≈ n[perm] rtol=1e-12
    end
    for t in (0.0, 1e-24, 1e8, 1e16, 1e100), p in (0.0, 1e-30, 1.0, 1e10)
        @test molar_wall_flux(SVector(3.0), t, SVector(p), SVector(2.0)) == SVector(6.0)
    end
    K=SVector(1.0, 2.0, 3.0)
    @test molar_wall_flux(K, 1e5, zero(K), K)==K .* K
    @test molar_wall_flux(zero(K), 1e5, K, K)==zero(K)
end

@testset "donor half path: vacuum, dissipation and derivatives" begin
    rng=MersenneTwister(815)
    for trial = 1:2000
        K=SVector{3}(10.0 .^ (6rand(rng, 3) .- 3))
        l=SVector{3}(10.0 .^ (8rand(rng, 3) .- 4))
        r=SVector{3}(10.0 .^ (8rand(rng, 3) .- 4))
        t=10.0^(12rand(rng)-6)
        mob=10.0^(8rand(rng)-4)
        n=pore_half_flux(l, r, 500.0, 500.0, 500.0, 1.0, K, t, mob).n
        terms=n .* log.(r ./ l)
        @test sum(terms) ≤ 1e-13*sum(abs, terms)
        # An absent species cannot leave its empty endpoint.
        j=rand(rng, 1:3)
        le=setindex(l, 0.0, j)
        @test pore_half_flux(le, r, 500.0, 500.0, 500.0, 1.0, K, t, mob).n[j] ≤ 1e-14
    end
    K=SVector(1.0, 2.0, 3.0)
    f =
        z->pore_half_flux(
            z[SVector(1, 2, 3)],
            z[SVector(4, 5, 6)],
            1/R_GAS,
            1/R_GAS,
            1/R_GAS,
            1.0,
            K,
            10.0,
            1.0,
        ).n
    @test f(zeros(SVector{6,Float64}))==zero(K)
    J=ForwardDiff.jacobian(f, zeros(SVector{6,Float64}))
    @test all(isfinite, J)
    @test J ≈ hcat(Diagonal(K), -Diagonal(K))
    for z in
        (SVector(1.0, 2.0, 3.0, 2.0, 1.0, 3.0), SVector(1e-10, 0.0, 0.0, 2e-10, 0.0, 0.0))
        @test all(isfinite, ForwardDiff.jacobian(f, z))
    end
end

@testset "donor face uses one pressure state: original-input reference" begin
    rng=MersenneTwister(717)
    for exponent in (-24, -8, 0, 8, 16), trial = 1:30
        kl=SVector{3}(10 .^ (-2 .+ 4rand(rng, 3)))
        pl=SVector{3}(10 .^ (-6 .+ 8rand(rng, 3)))
        pr=SVector{3}(10 .^ (-6 .+ 8rand(rng, 3)))
        Tl, Tr, d=300+600rand(rng), 300+600rand(rng), 0.01+rand(rng)
        t=10.0^exponent
        mob=1e-4
        numerical=pore_half_flux(pl, pr, Tl, Tr, (Tl+Tr)/2, d, kl, t, mob).n
        reference=setprecision(256) do
            # Every dependent pressure, donor fraction and matrix entry is
            # formed here from the original endpoint inputs.
            l=BigFloat.(pl)
            r=BigFloat.(pr)
            K=BigFloat.(kl)
            TL=BigFloat(Tl)
            TR=BigFloat(Tr)
            D=BigFloat(d)
            tt=BigFloat(t)
            R=BigFloat(R_GAS)
            TB=(TL+TR)/2
            P=(sum(l)+sum(r))/2
            grad=(l-r)/(R*TB*D)
            diag=K ./ (1 .+ tt*P .* K)
            donor=sum(diag .* grad)>=0 ? l : r
            p=P*donor/sum(donor)
            A=Matrix(Diagonal(inv.(K) .+ tt*P))-tt*collect(p)*ones(BigFloat, 1, 3)
            n=A\collect(grad)
            v=BigFloat(mob)*(sum(l)-sum(r))/D
            n+=v*collect(v>=0 ? l/(R*TL) : r/(R*TR))
            Float64.(n)
        end
        for j = 1:3
            @test abs(numerical[j]-reference[j])<=5e-13*max(
                abs(reference[j]),
                maximum(abs, reference),
            )
        end
    end
end
