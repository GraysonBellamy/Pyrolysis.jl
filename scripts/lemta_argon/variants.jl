# The comparison's variants (included by run_cases.jl and compare.jl after
# LemtaArgon.jl): the decks as they stand, their convergence checks, one lever
# each within its documented uncertainty, and one diagnostic. Each variant has
# a label, a `kind`, the tests it applies to (`runs(r)`) and the keyword
# arguments of `simulate` for a test (`kwargs(r)`).

"""
The conductivity the thesis fitted to these cone tests' mass loss, one k(T) for
wood, char and intermediates per species (Flity 2023 thesis, eq. 5.19 and Fig.
5.14, read from the figure's vector paths): W/(m K) at 27 and 200 °C, and above
200 °C k(200) + λ1 (T − 200) + λ2 (T − 200)², T in °C, with λ1, λ2 fitted to
the curve's 200–500 °C points (rms ≤ 1e-4 W/(m K)); held on the same law above
500 °C, where the figure stops. **Fitted to the data the comparison scores**:
a diagnostic, never a deck value.
"""
const THESIS_K_FITTED = Dict(
    "balsa" => (k27 = 0.0529, k200 = 0.0565, λ1 = 3e-7, λ2 = 1.86e-6),
    "poplar" => (k27 = 0.0842, k200 = 0.0958, λ1 = 3.06e-5, λ2 = 1.62e-6),
    "fir" => (k27 = 0.0927, k200 = 0.1066, λ1 = 2.51e-4, λ2 = 0.0),
    "larch" => (k27 = 0.1168, k200 = 0.1368, λ1 = 1.95e-5, λ2 = 1.23e-6),
    "oak" => (k27 = 0.1268, k200 = 0.1492, λ1 = 1.2e-6, λ2 = 9.49e-7),
    "beech" => (k27 = 0.1282, k200 = 0.1515, λ1 = 6.61e-4, λ2 = 0.0),
    "ipe" => (k27 = 0.1652, k200 = 0.1973, λ1 = 1.8e-6, λ2 = 4.03e-7),
)

function thesis_k(species)
    p = THESIS_K_FITTED[species]
    return T -> begin
        θ = T - 273.15
        θ <= 200 ? interp([27.0, 200.0], [p.k27, p.k200], θ) :
        p.k200 + p.λ1 * (θ - 200) + p.λ2 * (θ - 200)^2
    end
end

"""
The char-yield shift of the heats sensitivity [percentage points]: spruce's
DSC char yield goes from 19.06 % at 20 mL/min to 22.55 % at 5 mL/min (Flity
2023 thesis, Table 4.5, means of 3 and 2 tests); a 5 cm block holds its
volatiles longer than either.
"""
const HEATS_ΔYC = 3.5

# The convergence checks run one light, one mid and one dense species.
const CONVERGENCE_SPECIES = ("balsa", "larch", "ipe")

all_runs(r) = true
has_cone_shrinkage(r) = haskey(deck(r.species), "shrinkage_cone")

const VARIANTS = Dict{String,NamedTuple}(
    "deck" => (label = "decks", kind = :deck, runs = all_runs, kwargs = r -> (;)),
    "cells200" => (
        label = "200 cells",
        kind = :convergence,
        runs = r -> r.species in CONVERGENCE_SPECIES,
        kwargs = r -> (cells = 200,),
    ),
    "tight" => (
        label = "reltol = abstol = 1e-8, dtmax 1 s",
        kind = :convergence,
        runs = r -> r.species in CONVERGENCE_SPECIES,
        kwargs = r -> (reltol = 1e-8, abstol = 1e-8, dtmax = 1.0),
    ),
    "shrink_none" => (
        label = "no shrinkage (fixed mesh)",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (shrinkage = nothing,),
    ),
    "shrink_cone" => (
        label = "shrinkage of the cone tests (I9CL2R)",
        kind = :lever,
        runs = has_cone_shrinkage,
        kwargs = r -> (shrinkage = :cone,),
    ),
    "h_natural" => (
        label = "front h at natural convection (6 W/(m² K))",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (h_front = chamber()["window"]["h_window"],),
    ),
    "heats_yield" => (
        label = "heats at +3.5 points of char yield (six-species fits)",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (heats = shifted_heats(deck(r.species), HEATS_ΔYC),),
    ),
    "char_pooled" => (
        label = "char k: pooled law (poplar, larch)",
        kind = :lever,
        runs = r -> deck(r.species)["char_conductivity"] != "pooled",
        kwargs = r -> (char_conductivity = "pooled",),
    ),
    "char_low" => (
        label = "char k: least of the three measured chars",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (char_conductivity = "low",),
    ),
    "char_high" => (
        label = "char k: greatest of the three measured chars",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (char_conductivity = "high",),
    ),
    "alpha_published" => (
        label = "wood absorptivity as published (≈ 230 °C source)",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (wood_absorptivity = deck(r.species)["absorptivity"],),
    ),
    "back_tc40" => (
        label = "back held at the 40 mm record (temperature tests)",
        kind = :lever,
        runs = r -> r isa TemperatureRun && live(r, length(DEPTHS)),
        kwargs = r -> (back = :tc40,),
    ),
    "volatiles_c0" => (
        label = "volatiles without heat capacity",
        kind = :lever,
        runs = all_runs,
        kwargs = r -> (volatiles_c = 0.0,),
    ),
    "thesis_k" => (
        label = "DIAGNOSTIC: the thesis's k(T) fitted to these tests (Fig. 5.14)",
        kind = :diagnostic,
        runs = all_runs,
        kwargs = r -> (k_solid = thesis_k(r.species),),
    ),
)
