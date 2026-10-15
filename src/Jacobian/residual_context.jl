"""
The `ResidualContext` — Layer-1's dependency-graph carrier.

The residual has ordered derived fields:
`properties → reactions → gas_fluxes → BC implicit T_s → mesh velocity`.
An [`Operator`](@ref) without declared reads / writes on
derived fields recreates the same hidden-coupling problem the rebuild is
designed to dissolve.
The `ResidualContext` is the single explicit declaration point: each
operator declares what it `reads` and `writes`, and a topological sort
over the graph yields the assembly order.

A `ResidualContext` is constructed once per problem (at
`_build_jacobian_spec` time) and reused across Newton steps. The
declaration-graph fields (`providers`, `consumers`, `order`) carry no
per-step state. The `derived_fields` map IS per-step state: providers
push bare [`StructuredCoupling`](@ref) primitives into it during
assembly, consumers read them back. Cleared at the start of every
`_compute_jacobian!` call by the assembly driver.

# Fields

- `providers::Dict{Symbol, Operator}` — derived field → its writing
  operator. Built by walking each operator's `writes(op)`.
- `consumers::Dict{Symbol, Vector{Operator}}` — derived field → list of
  reading operators. Built by walking each operator's `reads(op)`.
- `order::Vector{Operator}` — topological order. Operator `j` appears
  after every operator that writes a field `j` reads.
- `derived_fields::Dict{Symbol, Vector{StructuredCoupling}}` — per-step
  cache of bare prefix-sum / rank-1 / suffix-sum primitives that
  upstream operators push during their `contribute_to_structured!`
  pass. Consumers look the field up by name and use the framework's
  [`compose_into!`](@ref) to wrap each primitive in a globally-embedded
  coupling. Decouples provider-side construction from consumer-side row
  scaling.

A field `s` with no provider is fine — it means `s` is read directly
from `state`. A field with multiple providers is an error (at most one
operator can produce a given field).
"""

"""
	ResidualContext

Holds the dependency graph for an operator list. Constructed via
[`build_residual_context`](@ref) at problem-build time.
"""
struct ResidualContext
    providers::Dict{Symbol,Operator}
    consumers::Dict{Symbol,Vector{Operator}}
    order::Vector{Operator}
    derived_fields::Dict{Symbol,Vector{StructuredCoupling}}
end

"""
	build_residual_context(ops::Vector{<:Operator}) -> ResidualContext

Construct a `ResidualContext` from `ops`. Validates that each derived
field has at most one provider, then performs a topological sort on
the read/write dependencies. Errors on duplicate providers or
dependency cycles.
"""
function build_residual_context(ops::AbstractVector{<:Operator})
    providers = Dict{Symbol,Operator}()
    consumers = Dict{Symbol,Vector{Operator}}()

    for op in ops
        for f in writes(op)
            if haskey(providers, f)
                error(
                    "ResidualContext: derived field :$f has two providers " *
                    "($(typeof(providers[f])) and $(typeof(op))). At most " *
                    "one operator may write a given field.",
                )
            end
            providers[f] = op
        end
        for f in reads(op)
            push!(get!(consumers, f, Operator[]), op)
        end
    end

    order = _topological_sort(ops, providers)
    return ResidualContext(
        providers,
        consumers,
        order,
        Dict{Symbol,Vector{StructuredCoupling}}(),
    )
end

# Kahn's algorithm — operator j depends on op i iff some f ∈ reads(j) has
# providers[f] === i. Errors on cycles.
function _topological_sort(
    ops::AbstractVector{<:Operator},
    providers::Dict{Symbol,Operator},
)
    # Build adjacency: deps[op] = list of operators it must follow.
    indeg = Dict{Operator,Int}()
    out = Dict{Operator,Vector{Operator}}()
    for op in ops
        indeg[op] = 0
        out[op] = Operator[]
    end
    for op in ops
        for f in reads(op)
            provider = get(providers, f, nothing)
            provider === nothing && continue   # field comes from state
            provider === op && continue         # self-edge — ignore
            push!(out[provider], op)
            indeg[op] += 1
        end
    end

    queue = Operator[op for op in ops if indeg[op] == 0]
    order = Operator[]
    sizehint!(order, length(ops))
    while !isempty(queue)
        op = popfirst!(queue)
        push!(order, op)
        for downstream in out[op]
            indeg[downstream] -= 1
            if indeg[downstream] == 0
                push!(queue, downstream)
            end
        end
    end

    if length(order) != length(ops)
        # Some node still has indeg > 0 → cycle. Report which fields are involved.
        stuck = [op for op in ops if indeg[op] > 0]
        error(
            "ResidualContext: dependency cycle among operators " *
            "$(typeof.(stuck)). Check `reads`/`writes` declarations.",
        )
    end
    return order
end

"""
	provider_for(ctx::ResidualContext, field::Symbol) -> Union{Operator, Nothing}

Look up which operator writes `field`, or `nothing` if the field comes
directly from `state`. Lightweight accessor used by the AD-default
`contribute_to_jacobian!` to decide whether to compose with a
[`StructuredCoupling`](@ref) or fall back to local AD.
"""
@inline provider_for(ctx::ResidualContext, field::Symbol) =
    get(ctx.providers, field, nothing)
