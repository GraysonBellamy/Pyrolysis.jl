# Actual initialized fields, not a search for component keywords. This audit
# covers every verification builder at its actual CI mesh and composition.
# Run with --project=test; its TSV is an engineering qualification artifact.
using Pyrolysis, Test
using Pyrolysis.Problem: initial_fields
using Pyrolysis.Properties: skeletal_volume_fraction, transport_coefficients
using Pyrolysis.Integration: gas_can_move
include(joinpath(@__DIR__, "../verification/casedefs.jl"))

function audit_initial(io, label, p)
    T, ξ=initial_fields(p, p.mesh)
    m=p.material
    gas=any(is_gas, m.components)
    φ=[1-skeletal_volume_fraction(Tuple(ξ[:, i]), T[i], m) for i in eachindex(T)]
    κ=Float64[]
    K=Float64[]
    if gas
        @test minimum(φ)>0
        for i in eachindex(T)
            c=transport_coefficients(Tuple(ξ[:, i]), T[i], m, φ[i])
            push!(κ, c.κ)
            push!(K, maximum(c.K))
        end
    end
    println(
        io,
        join(
            (
                label,
                m.name,
                gas,
                gas_can_move(p),
                length(T),
                minimum(φ),
                maximum(φ),
                maximum(κ; init = 0.0),
                maximum(K; init = 0.0),
            ),
            '\t',
        ),
    )
    flush(io)
end

output=isempty(ARGS) ? "/tmp/pyrolysis-current-initial-pore-states.tsv" : ARGS[1]
open(output, "w") do io
    println(
        io,
        "case\tmaterial\thas_gas\tgas_can_move\tcells\tminimum_porosity\tmaximum_porosity\tmaximum_permeability\tmaximum_wall_coefficient",
    )
    @testset "initialized verification pore states" begin
        for c in load_verification_cases()
            p=Base.invokelatest(c.build, c.ci.n_cells)
            Base.invokelatest(audit_initial, io, c.id, p)
        end
    end
end
