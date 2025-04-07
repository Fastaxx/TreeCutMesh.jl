module TreeCutMesh

using CairoMakie
using LinearAlgebra
using BenchmarkTools

# Export main functionality
export ThreadedQuadTreeCell, 
       create_root_cell,
       refine_cell_whitney_with_constraints!,
       balance_quadtree!,
       build_constrained_quadtree,
       get_leaf_cells,
       level_set_circle,
       level_set_flower,
       is_mixed_cell,
       CellGeometry,
       compute_geometric_fractions,
       visualize_quadtree,
       visualize_geometric_fractions,
       visualize_interface_zoom

# Include submodules
include("quadtree_types.jl")
include("level_set_functions.jl")
include("mesh_operations.jl")
include("geometric_fractions.jl")
include("visualization.jl")

end