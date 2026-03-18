"""
Douglas Fir Cone Calorimeter: ThermaKin vs Julia Comparison

Parses ThermaKin output (DF_sensitivity_Cone_output.txt) and Julia results
(douglas_fir_results_50kWm2_fit.csv) and generates comparison plots.

Usage:
  julia --project=. examples/df_sensitivity_comparison.jl
"""

using Plots
using DelimitedFiles

# ============================================================================
# ThermaKin Output Parser (reused from ptw_comparison.jl)
# ============================================================================

"""
    parse_thermakin_output(filepath) -> Dict

Parse ThermaKin text output file and extract time-series data.
"""
function parse_thermakin_output(filepath::String)
	lines = readlines(filepath)

	times = Float64[]
	total_mass = Float64[]
	thickness = Float64[]
	heat_flow_in = Float64[]
	surface_temp = Float64[]
	back_temp = Float64[]

	component_names = String[]
	mlr_components = Dict{String, Vector{Float64}}()
	total_mlr = Float64[]

	names_parsed = false

	i = 1
	while i <= length(lines)
		line = lines[i]

		if startswith(line, "Time [s] =")
			t = parse(Float64, split(line, "=")[2])
			push!(times, t)

			# Find BOUNDARY block
			i += 1
			while i <= length(lines) && !startswith(lines[i], "BOUNDARY")
				i += 1
			end

			# Parse component names from header (once)
			if !names_parsed && i <= length(lines)
				header = lines[i]
				parts = split(header, '\t')
				found_mass = false
				for p in parts
					p_stripped = strip(p)
					if occursin("MASS FLOW OUT", p_stripped)
						found_mass = true
						continue
					end
					if found_mass && !isempty(p_stripped)
						push!(component_names, p_stripped)
						mlr_components[p_stripped] = Float64[]
					end
				end
				names_parsed = true
			end

			# Parse TOP boundary line
			i += 1
			if i <= length(lines) && startswith(lines[i], "TOP")
				top_parts = split(lines[i], '\t')
				vals = Float64[]
				for p in top_parts
					p_stripped = strip(p)
					if p_stripped == "TOP" || isempty(p_stripped)
						continue
					end
					try
						push!(vals, parse(Float64, p_stripped))
					catch
					end
				end
				if length(vals) >= 2
					push!(heat_flow_in, vals[2])
					comp_flows = vals[3:end]
					for (j, name) in enumerate(component_names)
						if j <= length(comp_flows)
							push!(mlr_components[name], comp_flows[j])
						else
							push!(mlr_components[name], 0.0)
						end
					end
					push!(total_mlr, sum(comp_flows))
				end
			end

			# Skip to spatial profile header
			while i <= length(lines) && !startswith(lines[i], "FROM TOP")
				i += 1
			end

			# Parse spatial cells for surface and back temperatures
			i += 1
			cell_temps = Float64[]
			while i <= length(lines) && !isempty(strip(lines[i])) && !startswith(lines[i], "Total")
				cell_parts = split(lines[i], '\t')
				cell_vals = Float64[]
				for p in cell_parts
					p_stripped = strip(p)
					if isempty(p_stripped)
						continue
					end
					try
						push!(cell_vals, parse(Float64, p_stripped))
					catch
					end
				end
				if length(cell_vals) >= 2
					push!(cell_temps, cell_vals[2])
				end
				i += 1
			end

			if !isempty(cell_temps)
				push!(surface_temp, cell_temps[1])
				if length(cell_temps) >= 2
					push!(back_temp, cell_temps[end-1])
				else
					push!(back_temp, cell_temps[1])
				end
			end

			# Parse Total thickness and Total mass
			while i <= length(lines) && !startswith(lines[i], "Total thickness")
				i += 1
			end
			if i <= length(lines)
				push!(thickness, parse(Float64, split(lines[i], "=")[2]))
			end
			i += 1
			if i <= length(lines) && startswith(lines[i], "Total mass")
				push!(total_mass, parse(Float64, split(lines[i], "=")[2]))
			end
		end

		i += 1
	end

	return Dict(
		:time => times,
		:total_mass => total_mass,
		:thickness => thickness,
		:mlr => total_mlr,
		:heat_flow_in => heat_flow_in,
		:surface_temp => surface_temp,
		:back_temp => back_temp,
		:component_names => component_names,
		:mlr_components => mlr_components,
	)
end

# ============================================================================
# Load Julia CSV Results
# ============================================================================

function load_julia_results(filepath::String)
	header_line = readline(filepath)
	headers = split(header_line, ',')
	data = readdlm(filepath, ','; skipstart=1)
	result = Dict{String, Vector{Float64}}()
	for (j, h) in enumerate(headers)
		result[strip(h)] = data[:, j]
	end
	return result
end

# ============================================================================
# Main
# ============================================================================

tk_file = joinpath(@__DIR__, "..", "docs", "DF_sensitivity_Cone_output.txt")
julia_file = joinpath(@__DIR__, "..", "douglas_fir_results_50kWm2_fit.csv")

println("Parsing ThermaKin output...")
tk = parse_thermakin_output(tk_file)
println("  - $(length(tk[:time])) time steps, t = [$(tk[:time][1]), $(tk[:time][end])] s")
println("  - Components: $(join(tk[:component_names], ", "))")

println("Loading Julia results...")
jl = load_julia_results(julia_file)
println("  - $(length(jl["time"])) time steps, t = [$(jl["time"][1]), $(jl["time"][end])] s")

# ============================================================================
# Comparison Plots
# ============================================================================

println("\nGenerating comparison plots...")

default(
	fontfamily = "Computer Modern",
	linewidth = 2,
	framestyle = :box,
	grid = true,
	gridalpha = 0.3,
	legend = :best,
	size = (600, 400),
)

tk_style = (color = :blue, label = "ThermaKin")
jl_style = (color = :red, linestyle = :dash, label = "Julia")

# --- Plot 1: Mass Loss Rate ---
p1 = plot(
	tk[:time], tk[:mlr] .* 1000,
	xlabel = "Time [s]", ylabel = "MLR [g/(m²·s)]",
	title = "Mass Loss Rate",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["mlr"] .* 1000,
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Plot 2: Surface Temperature ---
p2 = plot(
	tk[:time], tk[:surface_temp],
	xlabel = "Time [s]", ylabel = "Temperature [K]",
	title = "Surface Temperature",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["front_temp"],
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Plot 3: Back Temperature ---
p3 = plot(
	tk[:time], tk[:back_temp],
	xlabel = "Time [s]", ylabel = "Temperature [K]",
	title = "Back Temperature",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["back_temp"],
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Plot 4: Total Mass ---
p4 = plot(
	tk[:time], tk[:total_mass] .* 1000,
	xlabel = "Time [s]", ylabel = "Mass [g/m²]",
	title = "Total Mass",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["mass"] .* 1000,
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Plot 5: Thickness ---
p5 = plot(
	tk[:time], tk[:thickness] .* 1000,
	xlabel = "Time [s]", ylabel = "Thickness [mm]",
	title = "Sample Thickness",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["thickness"] .* 1000,
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Plot 6: Heat Flow In ---
p6 = plot(
	tk[:time], tk[:heat_flow_in] ./ 1000,
	xlabel = "Time [s]", ylabel = "Heat Flow [kW/m²]",
	title = "Net Heat Flow In (Top Surface)",
	color = tk_style.color, label = tk_style.label,
)
plot!(jl["time"], jl["q_net_surface"] ./ 1000,
	color = jl_style.color, linestyle = jl_style.linestyle, label = jl_style.label,
)

# --- Combined figure ---
p_combined = plot(p1, p2, p3, p4, p5, p6,
	layout = (2, 3),
	size = (1500, 800),
	plot_title = "Douglas Fir Cone Calorimeter: ThermaKin vs Julia",
	margin = 5Plots.mm,
)

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

savefig(p_combined, joinpath(outdir, "df_sensitivity_tk_vs_julia.png"))
println("Saved: output/df_sensitivity_tk_vs_julia.png")

# ============================================================================
# Component-Wise MLR Comparison
# ============================================================================

println("Generating component MLR comparison...")

# Map TK component names to Julia CSV column names
gas_map = [
	("WATER_VAPOR", "mlr_gen_water", "Water vapor"),
	("DF_GAS1", "mlr_gen_gas1", "Gas 1"),
	("DF_GAS2", "mlr_gen_gas2", "Gas 2"),
	("DF_GAS3", "mlr_gen_gas3", "Gas 3"),
	("DF_GAS4", "mlr_gen_gas4", "Gas 4"),
	("DF_GAS5", "mlr_gen_gas5", "Gas 5"),
]
colors = [:cyan, :blue, :green, :orange, :purple, :brown]

p_comp = plot(
	xlabel = "Time [s]",
	ylabel = "MLR [g/(m²·s)]",
	title = "Component MLR: ThermaKin (solid) vs Julia (dashed)",
	legend = :topright,
	size = (800, 500),
)

for (k, (tk_name, jl_col, display_name)) in enumerate(gas_map)
	if haskey(tk[:mlr_components], tk_name) && haskey(jl, jl_col)
		plot!(tk[:time], tk[:mlr_components][tk_name] .* 1000,
			color = colors[k], label = display_name, linewidth = 2)
		plot!(jl["time"], jl[jl_col] .* 1000,
			color = colors[k], linestyle = :dash, label = "", linewidth = 1.5)
	elseif haskey(tk[:mlr_components], tk_name)
		plot!(tk[:time], tk[:mlr_components][tk_name] .* 1000,
			color = colors[k], label = "$(display_name) (TK only)", linewidth = 2)
	end
end

savefig(p_comp, joinpath(outdir, "df_sensitivity_component_mlr.png"))
println("Saved: output/df_sensitivity_component_mlr.png")

# ============================================================================
# Quantitative Summary
# ============================================================================

println("\n" * "="^70)
println("  Quantitative Comparison Summary — Douglas Fir")
println("="^70)

# Peak MLR
tk_peak_mlr = maximum(tk[:mlr])
jl_peak_mlr = maximum(jl["mlr"])
println("\nPeak MLR:")
println("  ThermaKin: $(round(tk_peak_mlr * 1000, digits=3)) g/(m²·s)")
println("  Julia:     $(round(jl_peak_mlr * 1000, digits=3)) g/(m²·s)")
println("  Diff:      $(round((jl_peak_mlr - tk_peak_mlr) / tk_peak_mlr * 100, digits=2))%")

# Time to peak MLR
tk_t_peak = tk[:time][argmax(tk[:mlr])]
jl_t_peak = jl["time"][argmax(jl["mlr"])]
println("\nTime to Peak MLR:")
println("  ThermaKin: $(round(tk_t_peak, digits=1)) s")
println("  Julia:     $(round(jl_t_peak, digits=1)) s")
println("  Diff:      $(round(jl_t_peak - tk_t_peak, digits=1)) s")

# Final mass
tk_final_mass = tk[:total_mass][end]
jl_final_mass = jl["mass"][end]
println("\nFinal Mass:")
println("  ThermaKin: $(round(tk_final_mass * 1000, digits=2)) g/m²")
println("  Julia:     $(round(jl_final_mass * 1000, digits=2)) g/m²")
println("  Diff:      $(round((jl_final_mass - tk_final_mass) / tk_final_mass * 100, digits=2))%")

# Mass loss %
tk_mass_loss = (1 - tk[:total_mass][end] / tk[:total_mass][1]) * 100
jl_mass_loss = (1 - jl["mass"][end] / jl["mass"][1]) * 100
println("\nTotal Mass Loss:")
println("  ThermaKin: $(round(tk_mass_loss, digits=2))%")
println("  Julia:     $(round(jl_mass_loss, digits=2))%")

# Peak surface temperature
tk_peak_T = maximum(tk[:surface_temp])
jl_peak_T = maximum(jl["front_temp"])
println("\nPeak Surface Temperature:")
println("  ThermaKin: $(round(tk_peak_T, digits=1)) K")
println("  Julia:     $(round(jl_peak_T, digits=1)) K")
println("  Diff:      $(round(jl_peak_T - tk_peak_T, digits=1)) K")

# Final thickness
tk_final_thick = tk[:thickness][end]
jl_final_thick = jl["thickness"][end]
println("\nFinal Thickness:")
println("  ThermaKin: $(round(tk_final_thick * 1000, digits=3)) mm")
println("  Julia:     $(round(jl_final_thick * 1000, digits=3)) mm")
println("  Diff:      $(round((jl_final_thick - tk_final_thick) / tk_final_thick * 100, digits=2))%")

println("\n" * "="^70)
println("  Comparison complete!")
println("="^70)
