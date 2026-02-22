"""
ThermaKin Output Parser (Julia)

Extracts data from ThermaKin output files (.txt) into structured DataFrames.

Structure of ThermaKin output:
- Header with version, components file, conditions file, simulation info
- Time blocks separated by asterisks (******************)
  - Boundary data (FRONT/BACK): Area, Heat Flow, Mass Flow for each component
  - Heat flux info (FROM AXIS, FRONT/BACK HEAT FLUX)
  - Spatial profiles: position, temperature, concentrations for each component
  - Summary: Thickness, Mass, Total radius, Total mass
"""

module ThermaKinParser

using DataFrames
using Printf
using CSV

export ThermaKinHeader, ThermaKinData
export parse_thermakin_output, get_front_surface_temperature, get_back_surface_temperature
export get_mass_loss_rate, get_mass_loss_rate_per_area
export get_heat_flow_in, get_heat_flux_in, get_total_mass, get_thickness
export get_temperature_profile, get_component_mass_flow
export save_thermakin_data, analyze_thermakin_output, print_header_info, print_data_summary

"""
Header information from ThermaKin output file.
"""
Base.@kwdef mutable struct ThermaKinHeader
	version::String = ""
	components_file::String = ""
	conditions_file::String = ""
	n_components::Int = 0
	n_reactions::Int = 0
	object_type::String = ""
	n_layers::Int = 0
	start_time::String = ""
end

"""
Container for all extracted ThermaKin data.
"""
struct ThermaKinData
	header::ThermaKinHeader
	components::Vector{String}
	boundary_data::DataFrame
	heat_flux_data::DataFrame
	spatial_data::DataFrame
	summary_data::DataFrame
end

# Empty constructor
ThermaKinData() = ThermaKinData(
	ThermaKinHeader(),
	String[],
	DataFrame(),
	DataFrame(),
	DataFrame(),
	DataFrame(),
)

"""
	parse_thermakin_output(filepath::String) -> ThermaKinData

Parse a ThermaKin output file and extract all data.

# Arguments
- `filepath::String`: Path to the ThermaKin output file

# Returns
- `ThermaKinData`: Container with header, boundary data, heat flux, spatial profiles, and summary
"""
function parse_thermakin_output(filepath::String)
	content = read(filepath, String)

	# Parse header
	header = parse_header(content)

	# Split into time blocks by asterisks
	blocks = split(content, "*" ^ 74)

	# Filter to time blocks (skip header block)
	time_blocks = filter(b -> occursin("Time [s]", b), blocks[2:end])

	if isempty(time_blocks)
		@warn "No time blocks found in file"
		return ThermaKinData()
	end

	# Extract component names from first time block
	components = extract_components(first(time_blocks))

	# Parse all time blocks
	boundary_records = Vector{Dict{String, Any}}()
	heat_flux_records = Vector{Dict{String, Any}}()
	spatial_records = Vector{Dict{String, Any}}()
	summary_records = Vector{Dict{String, Any}}()
	for block in time_blocks
		time = extract_time(block)
		isnothing(time) && continue

		# Parse boundary data
		boundary = parse_boundary_data(block, time, components)
		append!(boundary_records, boundary)

		# Parse heat flux data
		heat_flux = parse_heat_flux(block, time)
		!isnothing(heat_flux) && push!(heat_flux_records, heat_flux)

		# Parse spatial profiles
		spatial = parse_spatial_data(block, time, components)
		append!(spatial_records, spatial)

		# Parse summary
		summary = parse_summary(block, time)
		!isnothing(summary) && push!(summary_records, summary)
	end

	# Convert to DataFrames
	boundary_df = isempty(boundary_records) ? DataFrame() : DataFrame(boundary_records)
	heat_flux_df = isempty(heat_flux_records) ? DataFrame() : DataFrame(heat_flux_records)
	spatial_df = isempty(spatial_records) ? DataFrame() : DataFrame(spatial_records)
	summary_df = isempty(summary_records) ? DataFrame() : DataFrame(summary_records)

	# Convert from_top coordinates to from_back for consistency
	# This ensures all downstream analysis uses the same convention:
	# - from_back = 0 is the back (substrate, z=0)
	# - from_back = thickness is the front (heated surface, z=L)
	if !isempty(spatial_df) && "from_top" in names(spatial_df) && any(spatial_df.from_top)
		spatial_df = transform(groupby(spatial_df, :time_s)) do df
			if any(df.from_top)
				# Convert: from_back = max(from_top) - from_top
				# max(from_top) is the thickness at this time
				thickness = maximum(df.from_back_m)
				df_new = copy(df)
				df_new.from_back_m = thickness .- df.from_back_m
				df_new.from_top .= false  # Now it's from_back convention
				df_new
			else
				df
			end
		end
	end

	return ThermaKinData(header, components, boundary_df, heat_flux_df, spatial_df, summary_df)
end

"""
Parse header information from content.
"""
function parse_header(content::String)
	header = ThermaKinHeader()

	# Version
	m = match(r"ThermaKin Program Version (\d+)", content)
	!isnothing(m) && (header.version = m.captures[1])

	# Components file
	m = match(r"Components file:\s+(\S+)", content)
	!isnothing(m) && (header.components_file = m.captures[1])

	# Conditions file
	m = match(r"Conditions file:\s+(\S+)", content)
	!isnothing(m) && (header.conditions_file = m.captures[1])

	# Number of components
	m = match(r"Number of components:\s+(\d+)", content)
	!isnothing(m) && (header.n_components = parse(Int, m.captures[1]))

	# Number of reactions
	m = match(r"Number of reactions:\s+(\d+)", content)
	!isnothing(m) && (header.n_reactions = parse(Int, m.captures[1]))

	# Object type
	m = match(r"Object type:\s+(\S+)", content)
	!isnothing(m) && (header.object_type = m.captures[1])

	# Number of layers
	m = match(r"Number of layers:\s+(\d+)", content)
	!isnothing(m) && (header.n_layers = parse(Int, m.captures[1]))

	# Start time
	m = match(r"Started on:\s+(.+)", content)
	!isnothing(m) && (header.start_time = strip(m.captures[1]))

	return header
end

"""
Extract component names from the boundary header line.
"""
function extract_components(block::AbstractString)
	# Clean up Windows line endings
	block = replace(block, "\r\n" => "\n")
	block = replace(block, "\r" => "\n")

	# Find the BOUNDARY header line with component names
	# Format: BOUNDARY  AREA [m^2]  HEAT FLOW IN [W]  MASS FLOW OUT [kg/s]:  comp1  comp2...
	m = match(r"BOUNDARY\s+AREA \[m\^2\]\s+HEAT FLOW IN \[W\]\s+MASS FLOW OUT \[kg/s\]:\s+(.+)", block)
	if !isnothing(m)
		components_str = strip(m.captures[1])
		components = filter(!isempty, strip.(split(components_str, '\t')))
		return collect(String, components)
	end

	return String[]
end

"""
Extract time value from block.
"""
function extract_time(block::AbstractString)
	# Clean up Windows line endings
	block = replace(block, "\r" => "")
	m = match(r"Time \[s\] =\s*([\d.eE+-]+)", block)
	isnothing(m) && return nothing
	return parse(Float64, m.captures[1])
end

"""
Parse FRONT/BACK or TOP/BOTTOM boundary data.
"""
function parse_boundary_data(block::AbstractString, time::Float64, components::AbstractVector{<:AbstractString})
	records = Dict{String, Any}[]

	# Clean up Windows line endings
	block = replace(block, "\r" => "")

	# Try both naming conventions: FRONT/BACK and TOP/BOTTOM
	# Normalize to FRONT/BACK for output consistency
	for (boundary_type, output_name) in [("FRONT", "FRONT"), ("BACK", "BACK"), ("TOP", "FRONT"), ("BOTTOM", "BACK")]
		# Match boundary line: TYPE  value  value  value  value...
		# Make sure we don't match header lines like "Top Boundary"
		pattern = Regex("^$(boundary_type)\\s+([\\d.eE+-]+)\\s+([\\d.eE+-]+)\\s+(.+)", "m")
		m = match(pattern, block)
		if !isnothing(m)
			# Skip if we already have this boundary (avoid duplicates when both TOP and FRONT match)
			if any(r -> r["boundary"] == output_name, records)
				continue
			end

			record = Dict{String, Any}(
				"time_s" => time,
				"boundary" => output_name,
				"area_m2" => parse(Float64, m.captures[1]),
				"heat_flow_in_W" => parse(Float64, m.captures[2]),
			)

			# Parse mass flow values for each component
			mass_flows = filter(!isempty, strip.(split(strip(m.captures[3]), '\t')))
			for (i, comp) in enumerate(components)
				key = "mass_flow_$(comp)_kg_s"
				if i <= length(mass_flows)
					try
						record[key] = parse(Float64, mass_flows[i])
					catch
						record[key] = NaN
					end
				else
					record[key] = NaN
				end
			end

			push!(records, record)
		end
	end

	return records
end

"""
Parse heat flux data.
"""
function parse_heat_flux(block::AbstractString, time::Float64)
	# Clean up Windows line endings
	block = replace(block, "\r" => "")

	record = Dict{String, Any}("time_s" => time)
	found_any = false

	# FROM AXIS
	m = match(r"FROM AXIS \[m\] =\s*([\d.eE+-]+)", block)
	if !isnothing(m)
		record["from_axis_m"] = parse(Float64, m.captures[1])
		found_any = true
	end

	# FRONT HEAT FLUX
	m = match(r"FRONT HEAT FLUX IN \[W/m\^2\] =\s*([\d.eE+-]+)", block)
	if !isnothing(m)
		record["front_heat_flux_W_m2"] = parse(Float64, m.captures[1])
		found_any = true
	end

	# BACK HEAT FLUX
	m = match(r"BACK HEAT FLUX IN \[W/m\^2\] =\s*([\d.eE+-]+)", block)
	if !isnothing(m)
		record["back_heat_flux_W_m2"] = parse(Float64, m.captures[1])
		found_any = true
	end

	return found_any ? record : nothing
end

"""
Parse spatial profile data (position, temperature, concentrations).
"""
function parse_spatial_data(block::AbstractString, time::Float64, components::AbstractVector{<:AbstractString})
	records = Dict{String, Any}[]

	# Clean up Windows line endings
	block = replace(block, "\r" => "")

	lines = split(block, '\n')
	in_spatial = false
	is_from_top = false  # Track coordinate convention

	for line in lines
		line_stripped = strip(line)

		# Start of spatial data section - check both conventions
		if startswith(line_stripped, "FROM BACK [m]")
			in_spatial = true
			is_from_top = false
			continue
		elseif startswith(line_stripped, "FROM TOP [m]")
			in_spatial = true
			is_from_top = true
			continue
		end

		# End markers
		if startswith(line_stripped, "----") || startswith(line_stripped, "Thickness") ||
		   startswith(line_stripped, "Total thickness") || startswith(line_stripped, "Total mass")
			in_spatial = false
			continue
		end

		# Skip header lines
		if startswith(line_stripped, "FROM AXIS")
			continue
		end

		if in_spatial && !isempty(line_stripped)
			# Parse spatial data row - split by tabs
			parts = split(line, '\t')
			parts = filter(!isempty, map(strip, parts))

			if length(parts) >= 2
				try
					position = parse(Float64, parts[1])

					record = Dict{String, Any}(
						"time_s" => time,
						"from_back_m" => position,  # Store position (will be from_top if is_from_top)
						"from_top" => is_from_top,   # Track which convention was used
						"temperature_K" => parse(Float64, parts[2]),
					)

					# Parse concentrations - they start at index 3
					for (i, comp) in enumerate(components)
						key = "conc_$(comp)_kg_m3"
						if i + 2 <= length(parts)
							try
								record[key] = parse(Float64, parts[i+2])
							catch
								record[key] = NaN
							end
						else
							record[key] = NaN
						end
					end

					push!(records, record)
				catch e
					# Skip malformed lines
					continue
				end
			end
		end
	end

	return records
end

"""
Parse summary data (thickness, mass, etc.).
"""
function parse_summary(block::AbstractString, time::Float64)
	# Clean up Windows line endings
	block = replace(block, "\r" => "")

	record = Dict{String, Any}("time_s" => time)
	found_any = false

	# Thickness - try both formats
	m = match(r"Thickness \[m\] =\s*([\d.eE+-]+)", block)
	if isnothing(m)
		m = match(r"Total thickness \[m\] =\s*([\d.eE+-]+)", block)
	end
	if !isnothing(m)
		record["thickness_m"] = parse(Float64, m.captures[1])
		found_any = true
	end

	# Mass per area - try both formats
	m = match(r"Mass \[kg/m\^2\] =\s*([\d.eE+-]+)", block)
	if isnothing(m)
		m = match(r"Total mass \[kg/m\^2\] =\s*([\d.eE+-]+)", block)
	end
	if !isnothing(m)
		record["mass_kg_m2"] = parse(Float64, m.captures[1])
		found_any = true
	end

	# Total radius (for 2D)
	m = match(r"Total radius \[m\]\s*=\s*([\d.eE+-]+)", block)
	if !isnothing(m)
		record["total_radius_m"] = parse(Float64, m.captures[1])
		found_any = true
	end

	# Total mass (for 2D)
	m = match(r"Total mass \[kg\] =\s*([\d.eE+-]+)", block)
	if !isnothing(m)
		record["total_mass_kg"] = parse(Float64, m.captures[1])
		found_any = true
	end

	return found_any ? record : nothing
end

# ============================================================================
# Convenience functions for analysis
# ============================================================================

"""
	get_front_surface_temperature(data::ThermaKinData) -> DataFrame

Extract front surface temperature vs time.
"""
function get_front_surface_temperature(data::ThermaKinData)
	if isempty(data.spatial_data)
		return DataFrame()
	end

	# Get maximum position (front surface, largest from_back) for each time
	result = combine(groupby(data.spatial_data, :time_s)) do df
		idx = argmax(df.from_back_m)
		(from_back_m = df.from_back_m[idx],
			front_temperature_K = df.temperature_K[idx])
	end

	return result
end

"""
	get_back_surface_temperature(data::ThermaKinData) -> DataFrame

Extract back surface temperature vs time.
"""
function get_back_surface_temperature(data::ThermaKinData)
	if isempty(data.spatial_data)
		return DataFrame()
	end

	# Get minimum position (back surface, smallest from_back) for each time
	result = combine(groupby(data.spatial_data, :time_s)) do df
		idx = argmin(df.from_back_m)
		(from_back_m = df.from_back_m[idx],
			back_temperature_K = df.temperature_K[idx])
	end

	return result
end

"""
	get_mass_loss_rate(data::ThermaKinData) -> DataFrame

Calculate total mass loss rate from front boundary.
"""
function get_mass_loss_rate(data::ThermaKinData)
	if isempty(data.boundary_data)
		return DataFrame()
	end

	front = filter(row -> row.boundary == "FRONT", data.boundary_data)

	# Sum all mass flow columns
	mass_cols = filter(col -> startswith(String(col), "mass_flow_"), names(front))

	if isempty(mass_cols)
		return DataFrame()
	end

	result = DataFrame(
		time_s = front.time_s,
		total_mass_loss_rate_kg_s = [sum(row[col] for col in mass_cols; init = 0.0) for row in eachrow(front)],
	)

	return result
end

"""
	get_mass_loss_rate_per_area(data::ThermaKinData) -> DataFrame

Calculate mass loss rate per unit area from front boundary [kg/(m²·s)].
"""
function get_mass_loss_rate_per_area(data::ThermaKinData)
	if isempty(data.boundary_data)
		return DataFrame()
	end

	front = filter(row -> row.boundary == "FRONT", data.boundary_data)

	# Sum all mass flow columns
	mass_cols = filter(col -> startswith(String(col), "mass_flow_"), names(front))

	if isempty(mass_cols)
		return DataFrame()
	end

	result = DataFrame(
		time_s = front.time_s,
		mlr_kg_m2_s = [sum(row[col] for col in mass_cols; init = 0.0) / row.area_m2
					   for row in eachrow(front)],
	)

	return result
end

"""
	get_heat_flow_in(data::ThermaKinData; boundary::String="FRONT") -> DataFrame

Extract heat flow into the specified boundary vs time.
"""
function get_heat_flow_in(data::ThermaKinData; boundary::String = "FRONT")
	if isempty(data.boundary_data)
		return DataFrame()
	end

	filtered = filter(row -> row.boundary == boundary, data.boundary_data)

	return DataFrame(
		time_s = filtered.time_s,
		heat_flow_in_W = filtered.heat_flow_in_W,
	)
end

"""
	get_heat_flux_in(data::ThermaKinData; boundary::String="FRONT") -> DataFrame

Extract heat flux (W/m²) into the specified boundary vs time.
"""
function get_heat_flux_in(data::ThermaKinData; boundary::String = "FRONT")
	if isempty(data.boundary_data)
		return DataFrame()
	end

	filtered = filter(row -> row.boundary == boundary, data.boundary_data)

	return DataFrame(
		time_s = filtered.time_s,
		heat_flux_in_W_m2 = filtered.heat_flow_in_W ./ filtered.area_m2,
	)
end

"""
	get_total_mass(data::ThermaKinData) -> DataFrame

Extract total mass per unit area vs time from summary data.
"""
function get_total_mass(data::ThermaKinData)
	if isempty(data.summary_data) || !("mass_kg_m2" in names(data.summary_data))
		return DataFrame()
	end

	return DataFrame(
		time_s = data.summary_data.time_s,
		mass_kg_m2 = data.summary_data.mass_kg_m2,
	)
end

"""
	get_thickness(data::ThermaKinData) -> DataFrame

Extract sample thickness vs time from summary data.
"""
function get_thickness(data::ThermaKinData)
	if isempty(data.summary_data) || !("thickness_m" in names(data.summary_data))
		return DataFrame()
	end

	return DataFrame(
		time_s = data.summary_data.time_s,
		thickness_m = data.summary_data.thickness_m,
	)
end

"""
	get_temperature_profile(data::ThermaKinData, time::Float64; tol::Float64=0.01) -> DataFrame

Extract temperature profile at a specific time (with tolerance).
"""
function get_temperature_profile(data::ThermaKinData, time::Float64; tol::Float64 = 0.01)
	if isempty(data.spatial_data)
		return DataFrame()
	end

	# Find closest time
	times = unique(data.spatial_data.time_s)
	idx = argmin(abs.(times .- time))

	if abs(times[idx] - time) > tol
		@warn "Requested time $time not found. Closest available: $(times[idx])"
	end

	actual_time = times[idx]
	profile = filter(row -> row.time_s ≈ actual_time, data.spatial_data)

	return sort(profile, :from_back_m)
end

"""
	get_component_mass_flow(data::ThermaKinData, component::String; boundary::String="FRONT") -> DataFrame

Extract mass flow rate for a specific component.
"""
function get_component_mass_flow(data::ThermaKinData, component::String; boundary::String = "FRONT")
	if isempty(data.boundary_data)
		return DataFrame()
	end

	col_name = Symbol("mass_flow_$(component)_kg_s")

	if !(col_name in propertynames(data.boundary_data))
		available = filter(col -> startswith(String(col), "mass_flow_"), names(data.boundary_data))
		error("Component '$component' not found. Available: $available")
	end

	filtered = filter(row -> row.boundary == boundary, data.boundary_data)

	return DataFrame(
		time_s = filtered.time_s,
		mass_flow_kg_s = filtered[!, col_name],
	)
end

# ============================================================================
# I/O Functions
# ============================================================================

"""
	save_thermakin_data(data::ThermaKinData, output_dir::String; prefix::String="thermakin")

Save all DataFrames to CSV files.
"""
function save_thermakin_data(data::ThermaKinData, output_dir::String; prefix::String = "thermakin")
	mkpath(output_dir)

	files_saved = String[]

	if !isempty(data.boundary_data)
		path = joinpath(output_dir, "$(prefix)_boundary.csv")
		CSV.write(path, data.boundary_data)
		push!(files_saved, path)
		println("Boundary data saved to: $path")
		println("  Shape: $(size(data.boundary_data))")
	end

	if !isempty(data.heat_flux_data)
		path = joinpath(output_dir, "$(prefix)_heat_flux.csv")
		CSV.write(path, data.heat_flux_data)
		push!(files_saved, path)
		println("\nHeat flux data saved to: $path")
		println("  Shape: $(size(data.heat_flux_data))")
	end

	if !isempty(data.spatial_data)
		path = joinpath(output_dir, "$(prefix)_spatial.csv")
		CSV.write(path, data.spatial_data)
		push!(files_saved, path)
		println("\nSpatial data saved to: $path")
		println("  Shape: $(size(data.spatial_data))")
	end

	if !isempty(data.summary_data)
		path = joinpath(output_dir, "$(prefix)_summary.csv")
		CSV.write(path, data.summary_data)
		push!(files_saved, path)
		println("\nSummary data saved to: $path")
		println("  Shape: $(size(data.summary_data))")
	end

	# Derived quantities
	front_temp = get_front_surface_temperature(data)
	if !isempty(front_temp)
		path = joinpath(output_dir, "$(prefix)_front_temp.csv")
		CSV.write(path, front_temp)
		push!(files_saved, path)
		println("\nFront surface temperature saved to: $path")
	end

	mlr = get_mass_loss_rate(data)
	if !isempty(mlr)
		path = joinpath(output_dir, "$(prefix)_mlr.csv")
		CSV.write(path, mlr)
		push!(files_saved, path)
		println("Mass loss rate saved to: $path")
	end

	return files_saved
end

"""
	print_header_info(data::ThermaKinData)

Print header information in a formatted way.
"""
function print_header_info(data::ThermaKinData)
	h = data.header
	println("ThermaKin Version: $(h.version)")
	println("Components file: $(h.components_file)")
	println("Conditions file: $(h.conditions_file)")
	println("Number of components: $(h.n_components)")
	println("Number of reactions: $(h.n_reactions)")
	println("Components: $(join(data.components, ", "))")
end

"""
	print_data_summary(data::ThermaKinData)

Print summary statistics about the parsed data.
"""
function print_data_summary(data::ThermaKinData)
	if !isempty(data.boundary_data)
		println("\nBoundary data shape: $(size(data.boundary_data))")
	end

	if !isempty(data.heat_flux_data)
		println("Heat flux data shape: $(size(data.heat_flux_data))")
	end

	if !isempty(data.spatial_data)
		println("Spatial data shape: $(size(data.spatial_data))")
	end

	if !isempty(data.summary_data)
		println("Summary data shape: $(size(data.summary_data))")
	end
end

"""
	analyze_thermakin_output(filepath::String; verbose::Bool=true) -> ThermaKinData

Parse and analyze a ThermaKin output file with optional verbose output.
"""
function analyze_thermakin_output(filepath::String; verbose::Bool = true)
	println("Parsing $filepath...")
	data = parse_thermakin_output(filepath)

	if verbose
		print_header_info(data)
	end

	return data
end

end # module ThermaKinParser


# ============================================================================
# Command-line interface when run as script
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
	using .ThermaKinParser
	using CSV

	if length(ARGS) < 1
		println("Usage: julia parse_thermakin_output.jl <input_file> [--output-dir DIR] [--prefix PREFIX]")
		println()
		println("Arguments:")
		println("  input_file     Path to ThermaKin output file (.txt)")
		println("  --output-dir   Output directory for CSV files (default: same as input)")
		println("  --prefix       Prefix for output files (default: 'thermakin')")
		exit(1)
	end

	# Parse arguments
	local input_file = ARGS[1]
	local output_dir = dirname(input_file)
	local prefix = "thermakin"

	local i = 2
	while i <= length(ARGS)
		if ARGS[i] == "--output-dir" && i + 1 <= length(ARGS)
			output_dir = ARGS[i+1]
			i += 2
		elseif ARGS[i] == "--prefix" && i + 1 <= length(ARGS)
			prefix = ARGS[i+1]
			i += 2
		else
			i += 1
		end
	end

	# Use current directory if output_dir is empty
	if isempty(output_dir)
		output_dir = "."
	end

	# Parse and analyze
	local data = analyze_thermakin_output(input_file; verbose = true)

	# Save to CSV
	println()
	local files = save_thermakin_data(data, output_dir; prefix = prefix)

	println("\nDone!")
end
