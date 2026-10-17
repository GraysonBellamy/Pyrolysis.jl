"""
	Assembly

The residual of the finite-volume core: the state mapping, the kernels in a
fixed order, the external sources and the row assembly, as a pure function
of `(du, u, def, t, ws)`.
"""
module Assembly

using ..Materials
using ..Meshes
using ..States
using ..Kernels

include("definition.jl")
include("residual.jl")

export NoSources, VolumetricSources
export CoreDefinition, CoreWorkspace, buffers
export residual!, assemble_rows!, add_sources!

end # module Assembly
