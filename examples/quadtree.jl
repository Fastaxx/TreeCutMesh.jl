using TreeCutMesh
using CairoMakie

"""
This example creates a mesh for a more complex "flower" shape
and demonstrates the advanced visualization capabilities.
"""

# Define domain size and mesh parameters
domain_size = (2.0, 2.0)
max_level = 7
min_cell_size = 0.001

# Define flower interface with 5 petals
flower_level_set = (x, y) -> level_set_flower(x, y, 1.0, 1.0, 0.5, 5, 0.3)

# Build quadtree with balanced refinement
println("Building quadtree mesh for flower shape...")
quadtree_flower = build_constrained_quadtree(0.0, 0.0, domain_size[1], domain_size[2], 
                                           flower_level_set, max_level=max_level, 
                                           min_cell_size=min_cell_size, lip_const=0.5)

# Get statistics about the mesh
leaves = get_leaf_cells(quadtree_flower)
println("Mesh contains $(length(leaves)) leaf cells")
println("Maximum refinement level: $(maximum(c.level for c in leaves))")

# Identify cells crossed by the interface
mixed_cells = filter(cell -> is_mixed_cell(cell, flower_level_set), leaves)
println("Number of mixed cells: $(length(mixed_cells))")

# Visualize the mesh
fig_mesh = visualize_quadtree(quadtree_flower, flower_level_set, domain_size, color_by_level=true)

# Compute geometric fractions
cell_geoms = compute_geometric_fractions(leaves, flower_level_set)

# Visualize geometric fractions
fig_fractions = visualize_geometric_fractions(quadtree_flower, flower_level_set, domain_size)

# Zoom view on the interface
fig_zoom = visualize_interface_zoom(quadtree_flower, flower_level_set, domain_size)

# Save figures
save("flower_mesh.png", fig_mesh)
save("flower_fractions.png", fig_fractions)
save("flower_interface_zoom.png", fig_zoom)

println("Done! Figures saved as PNG files.")

# Define circle level set function
circle_level_set = (x, y) -> level_set_circle(x, y, 1.0, 1.0, 0.5)

# Create a circle mesh
println("Building quadtree mesh for circle shape...")
quadtree_circle = build_constrained_quadtree(0.0, 0.0, domain_size[1], domain_size[2], 
                                           circle_level_set, max_level=max_level, 
                                           min_cell_size=min_cell_size, lip_const=0.5)

# Get statistics about the mesh
leaves_circle = get_leaf_cells(quadtree_circle)
println("Mesh contains $(length(leaves_circle)) leaf cells")
println("Maximum refinement level: $(maximum(c.level for c in leaves_circle))")

# Identify cells crossed by the interface
mixed_cells_circle = filter(cell -> is_mixed_cell(cell, circle_level_set), leaves_circle)
println("Number of mixed cells: $(length(mixed_cells_circle))")

# Visualize the mesh           
fig_mesh_circle = visualize_quadtree(quadtree_circle, circle_level_set, domain_size, color_by_level=true)

# Compute geometric fractions
cell_geoms_circle = compute_geometric_fractions(leaves_circle, circle_level_set)

# Visualize geometric fractions
fig_fractions_circle = visualize_geometric_fractions(quadtree_circle, circle_level_set, domain_size)

# Zoom view on the interface
fig_zoom_circle = visualize_interface_zoom(quadtree_circle, circle_level_set, domain_size)

# Save figures
save("circle_mesh.png", fig_mesh_circle)
save("circle_fractions.png", fig_fractions_circle)
save("circle_interface_zoom.png", fig_zoom_circle)
println("Done! Figures saved as PNG files.")