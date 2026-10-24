```@meta
CurrentModule = Pyrolysis
```

# Meshes

A mesh is face-based: faces carry their nodes, their two cells and a
boundary tag, cells their faces and nodes, and the node positions are the
geometric state. A line's back face is tagged `:bottom` (z = 0) and its
exposed face `:top`. The coordinate system and the cross section are
properties of the mesh. See the User Guide's
[meshes chapter](../guide/user_guide/06_meshes-geometry.md).

## Constructors

```@docs
create_uniform_mesh
create_stretched_mesh
create_mesh_from_nodes
```

## Coordinate systems

```@docs
Planar
Cylindrical
Spherical
```

## Topology and geometry

```@docs
MeshTopology
MeshGeometry
Meshes.Mesh
Meshes.geometry!
Meshes.box_mesh
Meshes.cell_mesh
```
