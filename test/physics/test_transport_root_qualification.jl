using Test, Pyrolysis, StaticArrays, LinearAlgebra, ForwardDiff, Random
using Pyrolysis.Physics:
    interface_flux, TransportFailure, bounded_transport_root, pore_half_flux
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

function reference_interface(left, right, Tl, Tr, dl, dr, cl, cr; bits = 256)
    setprecision(bits) do
        l, r=BigFloat.(left), BigFloat.(right)
        tl, tr, dL, dR=BigFloat(Tl), BigFloat(Tr), BigFloat(dl), BigFloat(dr)
        kl, kr=BigFloat.(cl.K), BigFloat.(cr.K)
        al, ar, ml, mr=BigFloat(cl.t), BigFloat(cr.t), BigFloat(cl.mob), BigFloat(cr.mob)
        ts=(dR*tl+dL*tr)/(dL+dR)
        leftflux=p->reference_half(l, p, tl, ts, tl, dL, kl, al, ml)
        rightflux=p->reference_half(p, r, ts, tr, tr, dR, kr, ar, mr)
        f=p->leftflux(p)-rightflux(p)
        # Each species' continuity is measured against its own flux, so a
        # trace species is resolved however small its pressure.
        mismatch=p->maximum(
            map(
                (a, b)->a==b ? zero(a) : abs(a-b)/max(abs(a), abs(b)),
                leftflux(p),
                rightflux(p),
            ),
        )
        # An absent species carries no flux and stays at zero pressure. The
        # step of every other species is solved in its own pressure, with
        # each equation normalized, so elimination resolves a trace species
        # relative to itself rather than to the largest pressure.
        absent=(l .== 0) .& (r .== 0)
        σ=collect(map((a, b, z)->z ? one(a) : max(a, b), l, r, absent))
        p=collect((dR*l+dL*r)/(dL+dR))
        for iteration = 1:80
            e=mismatch(p)
            if e<big"1e-50"
                return (p = copy(p), n = (leftflux(p)+rightflux(p))/2)
            end
            A=ForwardDiff.jacobian(f, p) .* transpose(σ)
            rows=[maximum(abs, A[i, :]) for i in axes(A, 1)]
            step=σ .* ((A ./ rows)\(f(p) ./ rows))
            step[absent] .= 0
            α=one(BigFloat)
            while minimum(p-α*step)<0 || mismatch(p-α*step)>e
                α/=2
                α>big"1e-40" || error("independent root line search failed")
            end
            p-=α*step
        end
        error("independent root failed to converge")
    end
end

# The interface's acceptance contract holds each species' continuity to 2e-10
# of the gross size of its flux, the sum of its terms' magnitudes. The error of
# each species flux against the reference, relative to the larger of the
# reference flux and that gross size (and floatmin), checked with a hundredfold
# margin by the callers.
function contract_error(actual, ref, l, r, Tl, Tr, dl, dr, cl, cr)
    Ts=(dr*Tl+dl*Tr)/(dl+dr)
    gross=max.(
        pore_half_flux(l, actual.p, Tl, Ts, Tl, dl, cl.K, cl.t, cl.mob).gross,
        pore_half_flux(actual.p, r, Ts, Tr, Tr, dr, cr.K, cr.t, cr.mob).gross,
    )
    return maximum(
        abs.(BigFloat.(actual.n) .- ref.n) ./ max.(abs.(ref.n), gross, floatmin(Float64)),
    )
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

@testset "a subnormal trace partial pressure" begin
    # Face inputs from a wood deck with five product gases, water vapour and
    # air at the onset of heating. The third product's partial pressure is
    # subnormal, so is its flux; two products are absent. Dividing that
    # species' Jacobian row by its flux overflows it.
    l=SVector(
        1.388129665992509e-10,
        1.1179743361572907e-105,
        1.430105215e-315,
        0.0,
        0.0,
        9.021194907543943,
        101336.13005567773,
    )
    r=SVector(
        1.4649917862771915e-10,
        1.5096174164360192e-105,
        2.86342976e-315,
        0.0,
        0.0,
        8.693840087487288,
        101403.56224167546,
    )
    tl, tr=300.0001629922102, 300.34914587667856
    dl, dr=0.0002499999897919644, 0.00024999998967156246
    cl=(
        K = SVector(
            1.0860591807166755e-5,
            1.0054954190500078e-5,
            8.02107863077249e-6,
            1.0054954190500078e-5,
            1.1897182241282057e-5,
            1.982863706880343e-5,
            1.562985210998852e-5,
        ),
        t = 1.574700644525456,
        mob = 5.4173616031903004e-11,
    )
    cr=(
        K = SVector(
            1.0866906904765689e-5,
            1.0060800834789198e-5,
            8.025742639447286e-6,
            1.0060800834789198e-5,
            1.1904100084097624e-5,
            1.9840166806829373e-5,
            1.563894038466857e-5,
        ),
        t = 1.5715000947014357,
        mob = 5.412521117726331e-11,
    )
    ref=reference_interface(l, r, tl, tr, dl, dr, cl, cr)
    actual=interface_flux(l, r, tl, tr, dl, dr, cl, cr)
    @test minimum(actual.p)>=0
    @test contract_error(actual, ref, l, r, tl, tr, dl, dr, cl, cr)<=2e-8
    @test actual.n[4]==0 && actual.n[5]==0
    J=ForwardDiff.jacobian(p->interface_flux(p, r, tl, tr, dl, dr, cl, cr).n, l)
    @test all(isfinite, J)
end

@testset "a species flux that is a near cancellation" begin
    # Face inputs from the same wood deck. In the first, the second product is
    # uniform across the face; the diffusive total flux carries it one way and
    # Darcy flow the other, and its net flux is 3.5e5 times smaller than the
    # sum of its terms' magnitudes. In the second, the fifth product diffuses
    # against the outflow, its net flux 5.8e6 times smaller than that sum.
    # Forming such a net flux rounds it by about eps times its terms, so the
    # solve resolves it relative to that sum, not to itself.
    faces=(
        (
            l = SVector(
                5.41640721367963,
                1.7171907968380176e-65,
                1.3655027189980404e-224,
                0.0,
                0.0,
                28367.084311956154,
                94874.68347191087,
            ),
            r = SVector(
                5.486750687823508,
                1.7171908034306018e-65,
                1.3655027333213306e-224,
                0.0,
                0.0,
                28411.47258403331,
                94828.71410796158,
            ),
            tl = 299.64066818017795,
            tr = 299.640668189273,
            dl = 2.499765328020675e-5,
            dr = 2.4997653280206316e-5,
            cl = (
                K = SVector(
                    1.0854592870886702e-5,
                    1.0049400254714753e-5,
                    8.01664812270648e-6,
                    1.0049400254714753e-5,
                    1.1890610735838816e-5,
                    1.981768455973136e-5,
                    1.562121883396275e-5,
                ),
                t = 1.5778602528046815,
                mob = 5.4223607743092193e-11,
            ),
            cr = (
                K = SVector(
                    1.085459287105149e-5,
                    1.0049400254867315e-5,
                    8.016648122828183e-6,
                    1.0049400254867315e-5,
                    1.1890610736019332e-5,
                    1.981768456003222e-5,
                    1.5621218834199898e-5,
                ),
                t = 1.5778602527208534,
                mob = 5.4223607741826135e-11,
            ),
        ),
        (
            l = SVector(
                4804.866779798065,
                15188.3893486698,
                24383.753034092126,
                767.096539706449,
                11.502333392619017,
                23177.814634001104,
                49579.92291597831,
            ),
            r = SVector(
                4698.963018274633,
                14856.462794421062,
                23995.012248518586,
                777.1524038843719,
                11.696598095101796,
                22666.16553111233,
                50388.05995090284,
            ),
            tl = 703.2405191428177,
            tr = 705.5494775083631,
            dl = 1.892330159989665e-5,
            dr = 1.8572252936341696e-5,
            cl = (
                K = SVector(
                    7.014292878296442e-5,
                    6.493973332418308e-5,
                    5.180398611331346e-5,
                    6.493973332418308e-5,
                    7.6837728687816e-5,
                    0.00012806288114636,
                    0.00010094510712719557,
                ),
                t = 0.2630523006320036,
                mob = 3.02852729978399e-10,
            ),
            cr = (
                K = SVector(
                    7.120156555663104e-5,
                    6.5919840527602e-5,
                    5.2585841186570505e-5,
                    6.5919840527602e-5,
                    7.79974071700994e-5,
                    0.00012999567861683234,
                    0.0001024686278067672,
                ),
                t = 0.25808155772191116,
                mob = 3.1477038676972066e-10,
            ),
        ),
    )
    for f in faces
        ref=reference_interface(f.l, f.r, f.tl, f.tr, f.dl, f.dr, f.cl, f.cr)
        actual=interface_flux(f.l, f.r, f.tl, f.tr, f.dl, f.dr, f.cl, f.cr)
        @test minimum(actual.p)>=0
        @test contract_error(actual, ref, f.l, f.r, f.tl, f.tr, f.dl, f.dr, f.cl, f.cr)<=2e-8
        J=ForwardDiff.jacobian(
            p->interface_flux(p, f.r, f.tl, f.tr, f.dl, f.dr, f.cl, f.cr).n,
            f.l,
        )
        @test all(isfinite, J)
    end
end

@testset "species order does not change a face" begin
    # Wide-range faces whose total pressure difference is a cancellation of
    # large partial pressures: an exact equilibrium, and a difference of
    # 1.5e-37 Pa against partial pressures of 2e4 Pa. A sum rounded at each
    # step gives a nonzero flux at the equilibrium and a wrong total in other
    # species orders.
    equilibrium=SVector(
        4707.016104067481,
        1705.2674128762196,
        0.0,
        3.9666785116579826e-206,
        23894.282859413146,
        8410.931357168873,
        0.0,
    )
    ce=(
        K = SVector(
            5.441671518578719e-5,
            6.605093776828729e-5,
            6.516173384129243e-5,
            7.179633747327119e-5,
            7.132362834072358e-5,
            5.0543950741019136e-5,
            4.587824647121041e-5,
        ),
        t = 0.49918631964585064,
        mob = 2.6642307763299436e-10,
    )
    for order in (SVector(1, 2, 3, 4, 5, 6, 7), SVector(6, 1, 2, 5, 7, 4, 3))
        p=equilibrium[order]
        c=(K = ce.K[order], t = ce.t, mob = ce.mob)
        flux=interface_flux(
            p,
            p,
            674.2162139453269,
            674.2162139453269,
            1.4e-4,
            1.4e-4,
            c,
            c,
        )
        @test all(iszero, flux.n)
    end
    l=SVector(
        1195.1467747237953,
        23783.363092044074,
        0.0,
        8.135587761057449e-247,
        2.3236718242273295e-90,
        21073.15966036879,
        2.495580014507679e-36,
    )
    r=SVector(
        1195.1467747237953,
        23783.363092044074,
        0.0,
        9.765219107480036e-247,
        2.323671824430157e-90,
        21073.15966036879,
        2.6469636410949682e-36,
    )
    tl, tr=797.5546770094605, 793.5422361873897
    dl, dr=0.00011257745703750401, 0.00012119038084899428
    cl=(
        K = SVector(
            2.2603078547129284e-5,
            2.2410678311053414e-5,
            2.0817508033888874e-5,
            2.5191021759809883e-5,
            2.0964953374170812e-5,
            2.1549100278759034e-5,
            3.0735404888406795e-5,
        ),
        t = 2.904484734458638,
        mob = 2.0095183256870005e-10,
    )
    cr=(
        K = SVector(
            2.2624409367577443e-5,
            2.2431827560885812e-5,
            2.0837153788122e-5,
            2.5214794856069482e-5,
            2.0984738274497377e-5,
            2.1569436446149115e-5,
            3.076441029144113e-5,
        ),
        t = 2.9042802603416984,
        mob = 2.011291151141226e-10,
    )
    # The reference needs more than 256 bits to resolve this face.
    ref=reference_interface(l, r, tl, tr, dl, dr, cl, cr; bits = 512)
    for order in (SVector(1, 2, 3, 4, 5, 6, 7), SVector(7, 6, 1, 2, 3, 4, 5))
        permuted(c) = (K = c.K[order], t = c.t, mob = c.mob)
        flux=interface_flux(l[order], r[order], tl, tr, dl, dr, permuted(cl), permuted(cr))
        back=SVector{7}(invperm(order))
        unpermuted=(n = flux.n[back], p = flux.p[back])
        @test contract_error(unpermuted, ref, l, r, tl, tr, dl, dr, cl, cr)<=2e-8
    end
end

# A face of the lumped law reads signed pressures, so a trial state reaches it
# below vacuum.
@testset "the lumped law at and below vacuum" begin
    T, d=300.0, 1e-4
    one(K, mob) = (K = SVector(K), t = 0.0, mob = mob)
    two(K, mob) = (K = SVector(K, K), t = 0.0, mob = mob)
    flux(pl, pr, cl, cr = cl) = interface_flux(pl, pr, T, T, d, d, cl, cr)
    wall(K, Δ) = K*Δ/(R_GAS*T*2d)

    @testset "no viscous flow below vacuum" begin
        # Negative conductances would turn a tiny state into a flux of thousands.
        closed=flux(SVector(-1e3), SVector(-1e4), one(0.0, 1e-3))
        @test closed.n == SVector(0.0)
        @test flux(SVector(-1.0), SVector(-2.0), one(1e-6, 1e-8)).n[1] ≈ wall(1e-6, 1.0)
        # Flow runs down the gradient however strong the mobility.
        @test flux(SVector(-1.0), SVector(-2.0), one(1e-8, 5e-6)).n[1] ≈ wall(1e-8, 1.0)
        @test flux(SVector(-2.0), SVector(-1.0), one(1e-8, 5e-6)).n[1] ≈ -wall(1e-8, 1.0)
    end

    @testset "a total of exactly zero" begin
        for (pl, pr) in (
            (SVector(1.0, -1.0), SVector(-1.0, 0.5)),
            (SVector(1.0, -1.0), SVector(2.0, -2.0)),
        )
            r=flux(pl, pr, two(0.0, 1e-8))
            @test all(isfinite, r.n)
            @test all(isfinite, r.p)
        end
    end

    @testset "a closed half beside vacuum carries nothing" begin
        @test flux(SVector(1.0), SVector(0.0), one(0.0, 0.0), one(0.0, 1e-8)).n ==
              SVector(0.0)
    end

    @testset "the flux is continuous across each vacuum branch" begin
        λ, m=1e-6, 1e-7
        c=one(λ, m)
        continuous(Pl, Pr) = isapprox(
            flux(SVector(Pl*(1-1e-9)), SVector(Pr), c).n[1],
            flux(SVector(Pl*(1+1e-9)), SVector(Pr), c).n[1];
            rtol = 1e-7,
        )
        @test flux(SVector(-1e-9), SVector(-2.0), c).n[1] ≈
              flux(SVector(1e-9), SVector(-2.0), c).n[1] rtol=1e-6
        # The interface reaches vacuum where the left flux at P*=0, (λ+mPl)Pl
        # over Td, equals the right wall flux, 2λ over Td.
        threshold=(-λ+sqrt(λ^2+8m*λ))/(2m)
        @test continuous(threshold, -2.0)
        # The pure-Darcy interface closes to vacuum with the upstream total.
        @test abs(flux(SVector(1e-9), SVector(-2.0), one(0.0, m)).n[1]) <
              1e-8*flux(SVector(1.0), SVector(-2.0), one(0.0, m)).n[1]
    end

    @testset "derivatives stay finite on every branch" begin
        for x0 in ([-1.0, -2.0], [1.0, -2.0], [-1.0, 2.0], [2.0, 1.0])
            g=ForwardDiff.gradient(
                x->interface_flux(
                    SVector(x[1]),
                    SVector(x[2]),
                    T,
                    T,
                    d,
                    d,
                    one(1e-6, 1e-8),
                    one(1e-6, 1e-8),
                ).n[1],
                x0,
            )
            @test all(isfinite, g)
            @test g[1]>0 && g[2]<0
        end
    end
end
