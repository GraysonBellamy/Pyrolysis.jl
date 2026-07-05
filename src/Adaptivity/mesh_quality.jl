"""
Mesh quality metrics.

Mesh quality metrics ensure the mesh remains valid and well-conditioned
for numerical computations. Poor quality cells can cause accuracy loss
and solver instability.

## 1D Quality Metrics
- Aspect ratio: Δz_i / Δz_{i±1} should be bounded (typically < 2)
- Smoothness: |Δz_i - Δz_{i±1}| / Δz_i should be small (typically < 0.5)
"""

using ..Geometry: Mesh1D, cell_size, n_cells

# ============================================================================
# Quality Metric Types
# ============================================================================

"""
	AspectRatioMetric <: AbstractQualityMetric

Limits the ratio of adjacent cell sizes.

# Fields
- `max_ratio::Float64`: Maximum allowed ratio Δz_i / Δz_{i±1}
"""
struct AspectRatioMetric <: AbstractQualityMetric
    max_ratio::Float64
end
AspectRatioMetric() = AspectRatioMetric(2.0)

"""
	SmoothnessMetric <: AbstractQualityMetric

Limits the size jump between adjacent cells.

# Fields
- `max_size_jump::Float64`: Maximum allowed |Δz_i - Δz_{i±1}| / Δz_i
"""
struct SmoothnessMetric <: AbstractQualityMetric
    max_size_jump::Float64
end
SmoothnessMetric() = SmoothnessMetric(0.5)

"""
	JacobianMetric <: AbstractQualityMetric

Ensures positive cell volumes (Jacobian > 0).

# Fields
- `min_jacobian::Float64`: Minimum allowed volume (typically 0 or small positive)
"""
struct JacobianMetric <: AbstractQualityMetric
    min_jacobian::Float64
end
JacobianMetric() = JacobianMetric(0.0)

"""
	MinCellSizeMetric <: AbstractQualityMetric

Ensures cells don't become too small.

# Fields
- `h_min::Float64`: Minimum allowed cell size
"""
struct MinCellSizeMetric <: AbstractQualityMetric
    h_min::Float64
end
MinCellSizeMetric() = MinCellSizeMetric(1e-8)

"""
	CompositeQualityMetric <: AbstractQualityMetric

Combines multiple quality metrics.

# Fields
- `metrics::Vector{AbstractQualityMetric}`
"""
struct CompositeQualityMetric <: AbstractQualityMetric
    metrics::Vector{AbstractQualityMetric}
end

# ============================================================================
# Quality Checking
# ============================================================================

"""
	check_quality(mesh::Mesh1D, metric::AbstractQualityMetric) -> Bool

Check if mesh satisfies quality metric.

# Returns
- `true` if mesh passes quality check, `false` otherwise
"""
function check_quality end

"""
	quality_violations(mesh::Mesh1D, metric::AbstractQualityMetric) -> Vector{Int}

Find cells that violate quality metric.

# Returns
- Vector of cell indices that fail the quality check
"""
function quality_violations end

# ============================================================================
# Aspect Ratio
# ============================================================================

function check_quality(mesh::Mesh1D, metric::AspectRatioMetric)
    return isempty(quality_violations(mesh, metric))
end

function quality_violations(mesh::Mesh1D, metric::AspectRatioMetric)
    violations = Int[]
    n = n_cells(mesh)

    # Build list of active cells sorted by position
    active_cells = [(i, mesh.cell_centers[i]) for i = 1:n if mesh.active_cell_mask[i]]
    sort!(active_cells, by = x -> x[2])

    for idx = 1:(length(active_cells)-1)
        i = active_cells[idx][1]
        j = active_cells[idx+1][1]

        h_i = cell_size(mesh, i)
        h_j = cell_size(mesh, j)

        ratio = max(h_i / h_j, h_j / h_i)

        if ratio > metric.max_ratio
            push!(violations, i)
            push!(violations, j)
        end
    end

    return unique(violations)
end

"""
	compute_aspect_ratios(mesh::Mesh1D) -> Vector{Float64}

Compute aspect ratio for each cell (ratio to smaller neighbor).
"""
function compute_aspect_ratios(mesh::Mesh1D)
    n = n_cells(mesh)
    ratios = ones(n)

    active_cells = [(i, mesh.cell_centers[i]) for i = 1:n if mesh.active_cell_mask[i]]
    sort!(active_cells, by = x -> x[2])

    for idx = 1:(length(active_cells)-1)
        i = active_cells[idx][1]
        j = active_cells[idx+1][1]

        h_i = cell_size(mesh, i)
        h_j = cell_size(mesh, j)

        ratio = max(h_i / h_j, h_j / h_i)
        ratios[i] = max(ratios[i], ratio)
        ratios[j] = max(ratios[j], ratio)
    end

    return ratios
end

# ============================================================================
# Smoothness
# ============================================================================

function check_quality(mesh::Mesh1D, metric::SmoothnessMetric)
    return isempty(quality_violations(mesh, metric))
end

function quality_violations(mesh::Mesh1D, metric::SmoothnessMetric)
    violations = Int[]
    n = n_cells(mesh)

    active_cells = [(i, mesh.cell_centers[i]) for i = 1:n if mesh.active_cell_mask[i]]
    sort!(active_cells, by = x -> x[2])

    for idx = 1:(length(active_cells)-1)
        i = active_cells[idx][1]
        j = active_cells[idx+1][1]

        h_i = cell_size(mesh, i)
        h_j = cell_size(mesh, j)

        jump_i = abs(h_j - h_i) / h_i
        jump_j = abs(h_j - h_i) / h_j

        if jump_i > metric.max_size_jump
            push!(violations, i)
        end
        if jump_j > metric.max_size_jump
            push!(violations, j)
        end
    end

    return unique(violations)
end

"""
	compute_smoothness(mesh::Mesh1D) -> Vector{Float64}

Compute smoothness metric for each cell.
"""
function compute_smoothness(mesh::Mesh1D)
    n = n_cells(mesh)
    smoothness = zeros(n)

    active_cells = [(i, mesh.cell_centers[i]) for i = 1:n if mesh.active_cell_mask[i]]
    sort!(active_cells, by = x -> x[2])

    for idx = 1:(length(active_cells)-1)
        i = active_cells[idx][1]
        j = active_cells[idx+1][1]

        h_i = cell_size(mesh, i)
        h_j = cell_size(mesh, j)

        jump = abs(h_j - h_i)
        smoothness[i] = max(smoothness[i], jump / h_i)
        smoothness[j] = max(smoothness[j], jump / h_j)
    end

    return smoothness
end

# ============================================================================
# Jacobian (Positive Volume)
# ============================================================================

function check_quality(mesh::Mesh1D, metric::JacobianMetric)
    return isempty(quality_violations(mesh, metric))
end

function quality_violations(mesh::Mesh1D, metric::JacobianMetric)
    violations = Int[]
    A = mesh.cross_section_area

    @inbounds for i = 1:n_cells(mesh)
        if !mesh.active_cell_mask[i]
            continue
        end

        # Signed volume from node order: the cached cell_volumes pass through
        # abs(), so an inverted cell would otherwise satisfy every V > 0 check.
        cell = mesh.cells[i]
        width = mesh.nodes[cell.node_ids[2]].z - mesh.nodes[cell.node_ids[1]].z
        if width * A ≤ metric.min_jacobian
            push!(violations, i)
        end
    end

    return violations
end

# ============================================================================
# Minimum Cell Size
# ============================================================================

function check_quality(mesh::Mesh1D, metric::MinCellSizeMetric)
    return isempty(quality_violations(mesh, metric))
end

function quality_violations(mesh::Mesh1D, metric::MinCellSizeMetric)
    violations = Int[]

    @inbounds for i = 1:n_cells(mesh)
        if !mesh.active_cell_mask[i]
            continue
        end

        if cell_size(mesh, i) < metric.h_min
            push!(violations, i)
        end
    end

    return violations
end

# ============================================================================
# Composite Metric
# ============================================================================

function check_quality(mesh::Mesh1D, metric::CompositeQualityMetric)
    for sub_metric in metric.metrics
        if !check_quality(mesh, sub_metric)
            return false
        end
    end
    return true
end

function quality_violations(mesh::Mesh1D, metric::CompositeQualityMetric)
    all_violations = Int[]
    for sub_metric in metric.metrics
        append!(all_violations, quality_violations(mesh, sub_metric))
    end
    return unique(all_violations)
end

# ============================================================================
# Quality Enforcement
# ============================================================================

"""
	enforce_quality!(mesh::Mesh1D, metric::AbstractQualityMetric; 
					 max_iter::Int=10) -> Int

Attempt to fix quality violations through node smoothing.

# Arguments
- `mesh`: Mesh to modify
- `metric`: Quality metric to satisfy
- `max_iter`: Maximum smoothing iterations

# Returns
- Number of violations remaining after enforcement
"""
function enforce_quality!(mesh::Mesh1D, metric::AbstractQualityMetric; max_iter::Int = 10)

    for iter = 1:max_iter
        violations = quality_violations(mesh, metric)

        if isempty(violations)
            return 0
        end

        # Apply Laplacian smoothing to violating cells' nodes
        smoothed = false
        for cell_id in violations
            cell = mesh.cells[cell_id]

            # Smooth interior nodes only
            for node_id in cell.node_ids
                node = mesh.nodes[node_id]
                if node.is_boundary
                    continue
                end

                # Find neighbors
                z_neighbors = Float64[]
                for other_cell_id = 1:n_cells(mesh)
                    if !mesh.active_cell_mask[other_cell_id]
                        continue
                    end
                    other_cell = mesh.cells[other_cell_id]
                    if node_id in other_cell.node_ids
                        for other_node_id in other_cell.node_ids
                            if other_node_id != node_id
                                push!(z_neighbors, mesh.nodes[other_node_id].z)
                            end
                        end
                    end
                end

                if !isempty(z_neighbors)
                    z_avg = sum(z_neighbors) / length(z_neighbors)
                    z_new = 0.5 * (node.z + z_avg)  # Relaxed smoothing

                    mesh.nodes[node_id] = Node1D(
                        node.id,
                        z_new,
                        node.z_ref,
                        node.is_boundary,
                        node.boundary_id,
                    )
                    smoothed = true
                end
            end
        end

        if !smoothed
            break
        end

        update_geometry!(mesh)
    end

    return length(quality_violations(mesh, metric))
end

# ============================================================================
# Quality Summary
# ============================================================================

"""
	mesh_quality_summary(mesh::Mesh1D) -> NamedTuple

Compute summary statistics of mesh quality.

# Returns
NamedTuple with:
- `n_active_cells`: Number of active cells
- `min_cell_size`: Minimum cell size
- `max_cell_size`: Maximum cell size
- `max_aspect_ratio`: Maximum aspect ratio
- `max_smoothness`: Maximum smoothness deviation
- `all_positive_volume`: Whether all volumes are positive
"""
function mesh_quality_summary(mesh::Mesh1D)
    n = n_cells(mesh)

    sizes = [cell_size(mesh, i) for i = 1:n if mesh.active_cell_mask[i]]
    volumes = [mesh.cell_volumes[i] for i = 1:n if mesh.active_cell_mask[i]]

    aspect_ratios = compute_aspect_ratios(mesh)
    smoothness = compute_smoothness(mesh)

    active_ratios = [aspect_ratios[i] for i = 1:n if mesh.active_cell_mask[i]]
    active_smooth = [smoothness[i] for i = 1:n if mesh.active_cell_mask[i]]

    return (
        n_active_cells = length(sizes),
        min_cell_size = minimum(sizes; init = Inf),
        max_cell_size = maximum(sizes; init = 0.0),
        max_aspect_ratio = maximum(active_ratios; init = 1.0),
        max_smoothness = maximum(active_smooth; init = 0.0),
        all_positive_volume = all(v -> v > 0, volumes),
    )
end
