module MeshPhase0

# Executable candidates for the dated mesh design, not production replacements.
# Keep reference answers in the test files independent of these implementations.
using LinearAlgebra

export depletion_factor,
    reaction_rate,
    positive_mass,
    trial_temperature,
    temperature_from_inventory,
    merge_state,
    material_axes,
    structured_mesh,
    cell_volume,
    sweep_operator,
    closure_velocity

# The plan leaves the join scale open. This prototype uses x0=1. The cubic
# Hermite onset matches tanh and its first derivative at x0 and is flat at zero.
function depletion_factor(x; x0 = 1.0)
    x <= 0 && return zero(x)
    x >= x0 && return tanh(x)
    t = x / x0
    y, slope = tanh(x0), x0 / cosh(x0)^2
    return (slope - 2y) * t^3 + (3y - slope) * t^2
end

function softplus(x)
    return x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
end

# This C1 guard is constant outside the property bracket, identity away from
# one-kelvin endpoint transitions. The transition width is a prototype choice.
function trial_temperature(T; lo = 200.0, hi = 3000.0, width = 1.0)
    0 < 2width < hi - lo || throw(ArgumentError("invalid temperature transition"))
    T <= lo && return zero(T) + lo
    T >= hi && return zero(T) + hi
    if T < lo + width
        t = (T - lo) / width
        return lo + width * (2t^2 - t^3)
    elseif T > hi - width
        t = (hi - T) / width
        return hi - width * (2t^2 - t^3)
    end
    return T
end

function reaction_rate(x, order, T; threshold = 1.0, delta = 1e-8, x0 = 1.0)
    # Short circuit before evaluating fractional powers at an underflowed base.
    x <= 0 && return zero(x + T)
    base = delta * softplus(x / delta)
    return exp(-1000 / trial_temperature(T)) *
           base^order *
           depletion_factor(x / threshold; x0)
end

# One-sided C1 positive part: zero at and below zero mass, exactly m at and
# above delta, smoothstep onset between; slopes match at both joints.
function positive_mass(m, delta)
    delta > 0 || throw(ArgumentError("delta must be positive"))
    m <= 0 && return zero(m)
    m >= delta && return m
    t = m / delta
    return m * (3t^2 - 2t^3)
end

function capacity_floor(C, floor)
    C <= floor && return zero(C) + floor
    C >= 2floor && return C
    t = (C - floor) / floor
    return floor * (1 + 2t^2 - t^3)
end

# Constant-cp instance of the proposed extended caloric map. This is deliberately
# not presented as validation of nonlinear property-law inversion. For this
# instance the independent inverse and its derivatives are available analytically.
# The prototype chooses kappa from the regularized capacity, with an explicit
# C1 transition between floor and twice floor; this choice is not in production.
function temperature_from_inventory(
    H,
    masses,
    cp;
    delta = 1e-10,
    floor = 1e-12,
    Tref = 300.0,
)
    length(masses) == length(cp) || throw(DimensionMismatch("masses and cp"))
    all(>(0), cp) || throw(ArgumentError("heat capacities must be positive"))
    C = sum(positive_mass(m, delta) * c for (m, c) in zip(masses, cp))
    kappa = capacity_floor(C, floor) - C
    return Tref + H / (C + kappa)
end

# Minimal event contract with constant cp and temperature-dependent skeletal
# density. Gas is the last component and contributes to energy but not Vsolid.
# A merge never destroys pore space: the merged volume is the larger of the
# summed geometric volumes and the merged skeletal volume plus the summed
# pre-merge pore volumes, so created volume is the skeletal expansion only.
function merge_state(a, b, cp, densities; molar_mass = 0.028, R = 8.314462618)
    Tref = 300.0
    masses = a.m + b.m
    H = dot(a.m, cp) * (a.T - Tref) + dot(b.m, cp) * (b.T - Tref)
    T = Tref + H / dot(masses, cp)
    skeletal(m, temp) = sum(m[j] / densities[j](temp) for j in eachindex(densities))
    pore_pre = (a.V - skeletal(a.m, a.T)) + (b.V - skeletal(b.m, b.T))
    Vs = skeletal(masses, T)
    V = max(a.V + b.V, Vs + pore_pre)
    pore_volume = V - Vs
    pressure = masses[end] == 0 ? 0.0 : masses[end] / molar_mass * R * T / pore_volume
    return (;
        m = masses,
        T,
        H,
        V,
        Vs,
        phi = pore_volume / V,
        pressure,
        created_volume = V - a.V - b.V,
    )
end

# Ordered reference frame: grain, radial, tangential. The last axis enforces
# handedness. A singular F is outside this material-frame contract.
function material_axes(F, reference)
    size(F) == (3, 3) || throw(ArgumentError("three-dimensional frame required"))
    grain = F * reference[:, 1]
    grain /= norm(grain)
    radial = F * reference[:, 2]
    radial -= dot(radial, grain) * grain
    radial /= norm(radial)
    return hcat(grain, radial, cross(grain, radial))
end

function structured_mesh(counts::NTuple{D,Int}; lengths = ones(D)) where {D}
    D in (2, 3) || throw(ArgumentError("quadrilaterals and hexahedra only"))
    nodes = CartesianIndices(ntuple(k -> counts[k] + 1, D))
    ids = LinearIndices(nodes)
    x = hcat([[(I[k] - 1) * lengths[k] / counts[k] for k = 1:D] for I in nodes]...)
    corners = collect(Iterators.product(ntuple(_ -> (0, 1), D)...))[:]
    cells = [
        [ids[CartesianIndex(ntuple(k -> I[k] + c[k], D))] for c in corners] for
        I in CartesianIndices(counts)
    ][:]
    return (; x, cells, corners)
end

function shape_data(point, corners)
    d, n = length(point), length(corners)
    N = [prod(c[k] == 1 ? point[k] : 1 - point[k] for k = 1:d) for c in corners]
    grad = [
        (corners[a][k] == 1 ? 1 : -1) *
        prod(corners[a][j] == 1 ? point[j] : 1 - point[j] for j = 1:d if j != k) for
        a = 1:n, k = 1:d
    ]
    return N, grad
end

gauss_points(d) =
    Iterators.product(ntuple(_ -> ((1 - 1 / sqrt(3)) / 2, (1 + 1 / sqrt(3)) / 2), d)...)

# Independent volume path: volume quadrature, not the face operator below.
function cell_volume(x, cell, corners)
    d = size(x, 1)
    return sum(det(x[:, cell] * last(shape_data(q, corners))) for q in gauss_points(d)) /
           2^d
end

# Face quadrature of v.n, assembled directly into D. W integrates the nodal
# shape functions over incident cells (lumped nodal volume). No volume AD is
# used here; the tests compare this operator with the volume derivative.
function sweep_operator(mesh)
    x, cells, corners = mesh.x, mesh.cells, mesh.corners
    d, nn = size(x)
    D = zeros(eltype(x), length(cells), d * nn)
    weights = zeros(eltype(x), nn)
    for (i, cell) in enumerate(cells)
        xc = x[:, cell]
        for q in gauss_points(d)
            N, grad = shape_data(q, corners)
            weights[cell] += N * det(xc * grad) / 2^d
        end
        for axis = 1:d, side in (0, 1), q in gauss_points(d - 1)
            point = [k == axis ? Float64(side) : q[k < axis ? k : k - 1] for k = 1:d]
            N, grad = shape_data(point, corners)
            J = xc * grad
            normal_area = (2side - 1) * det(J) * inv(J)'[:, axis] / 2^(d - 1)
            for (a, node) in enumerate(cell), k = 1:d
                D[i, d*(node-1)+k] += N[a] * normal_area[k]
            end
        end
    end
    return D, repeat(weights; inner = d)
end

function closure_velocity(D, weights, source, free; weighted = true)
    A = D[:, free]
    rank(A) == size(A, 1) || throw(ArgumentError("closure has deficient row rank"))
    W = weighted ? weights[free] : ones(length(free))
    all(>(0), W) || throw(ArgumentError("nonpositive nodal volume"))
    gradient = A' ./ W
    phi = (A * gradient) \ source
    w = zeros(size(D, 2))
    w[free] = gradient * phi
    return w
end

end
