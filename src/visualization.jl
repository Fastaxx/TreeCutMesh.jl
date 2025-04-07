"""
    visualize_quadtree(root, level_set_func, domain_size; color_by_level=false)

Visualize the quadtree mesh with the level set interface.
- `color_by_level`: If true, cells are colored by their refinement level
"""
function visualize_quadtree(root, level_set_func, domain_size; color_by_level=false)
    leaves = get_leaf_cells(root)
    
    fig = Figure(size=(1000, 1000))
    ax = Axis(fig[1, 1], aspect=1, title="Quadtree mesh with Whitney decomposition (fully-threaded structure)")
    
    # Interface
    n_points = 2000
    x_range = range(0, domain_size[1], length=n_points)
    y_range = range(0, domain_size[2], length=n_points)
    z = [level_set_func(x, y) for y in y_range, x in x_range]
    contour!(ax, x_range, y_range, z, levels=[0], color=:red, linewidth=2)
    
    # Cells
    for cell in leaves
        rect = Rect(cell.x_min, cell.y_min, cell.width, cell.height)
        
        if color_by_level
            max_level = maximum([c.level for c in leaves])
            color_val = cell.level / max_level
            
            if level_set_func(cell.x_min + cell.width/2, cell.y_min + cell.height/2) < 0
                cell_color = RGBf(0.8*color_val, 0.2, 1-0.8*color_val)
                alpha = 0.4
            else
                cell_color = RGBf(0.2*color_val, 0.2+0.6*color_val, 0.9-0.5*color_val)
                alpha = 0.2
            end
            
            poly!(ax, rect, color=(cell_color, alpha), strokewidth=1, strokecolor=:black)
        else
            poly!(ax, rect, color=(:white, 0), strokewidth=1, strokecolor=:black)
        end
    end
    
    return fig
end

"""
    visualize_geometric_fractions(root, level_set_func, domain_size)

Visualize volume and face fractions for all cells in the mesh.
"""
function visualize_geometric_fractions(root, level_set_func, domain_size)
    leaves = get_leaf_cells(root)
    cell_geometries = compute_geometric_fractions(leaves, level_set_func, num_points=8)
    
    fig = Figure(size=(1800, 900))
    
    # Volume fraction
    ax1 = Axis(fig[1, 1], aspect=1, title="Volume fraction")
    
    # Interface
    n_points = 1000
    x_range = range(0, domain_size[1], length=n_points)
    y_range = range(0, domain_size[2], length=n_points)
    z = [level_set_func(x, y) for y in y_range, x in x_range]
    contour!(ax1, x_range, y_range, z, levels=[0], color=:black, linewidth=2)
    
    # Colormap for volume fractions
    for cell in leaves
        geom = cell_geometries[cell]
        rect = Rect(cell.x_min, cell.y_min, cell.width, cell.height)
        
        # Color gradient based on volume fraction
        color = RGBf(0.0, 0.0, geom.volume_fraction)
        alpha = 0.7
        
        poly!(ax1, rect, color=(color, alpha), strokewidth=1, strokecolor=:black)
        
        # For larger cells, display numerical value
        if cell.width > domain_size[1]/30 && cell.height > domain_size[2]/30
            text!(ax1, cell.x_min + cell.width/2, cell.y_min + cell.height/2, 
                  text="$(round(geom.volume_fraction, digits=2))", 
                  align=(:center, :center), color=:white, fontsize=12)
        end
    end
    
    # Surface fractions
    ax2 = Axis(fig[1, 2], aspect=1, title="Surface fractions")
    contour!(ax2, x_range, y_range, z, levels=[0], color=:black, linewidth=2)
    
    # Draw cell edges with color indicating surface fraction
    for cell in leaves
        geom = cell_geometries[cell]
        
        # Cell coordinates
        x_min, y_min = cell.x_min, cell.y_min
        x_max, y_max = x_min + cell.width, y_min + cell.height
        
        # North face (red)
        lines!(ax2, [x_min, x_max], [y_max, y_max], 
               color=RGBf(geom.face_fraction_north, 0.0, 0.0), 
               linewidth=3*geom.face_fraction_north + 0.5)
        
        # South face (green)
        lines!(ax2, [x_min, x_max], [y_min, y_min], 
               color=RGBf(0.0, geom.face_fraction_south, 0.0), 
               linewidth=3*geom.face_fraction_south + 0.5)
        
        # East face (blue)
        lines!(ax2, [x_max, x_max], [y_min, y_max], 
               color=RGBf(0.0, 0.0, geom.face_fraction_east), 
               linewidth=3*geom.face_fraction_east + 0.5)
        
        # West face (yellow)
        lines!(ax2, [x_min, x_min], [y_min, y_max], 
               color=RGBf(geom.face_fraction_west, geom.face_fraction_west, 0.0), 
               linewidth=3*geom.face_fraction_west + 0.5)
    end
    
    # Legend
    ax_legend = Axis(fig[1, 3], width=100, height=400)
    hidedecorations!(ax_legend)
    hidespines!(ax_legend)
    
    colors = [
        (RGBf(1.0, 0.0, 0.0), "North"),
        (RGBf(0.0, 1.0, 0.0), "South"),
        (RGBf(0.0, 0.0, 1.0), "East"),
        (RGBf(1.0, 1.0, 0.0), "West")
    ]
    
    for (i, (color, label)) in enumerate(colors)
        lines!(ax_legend, [0.1, 0.9], [i, i], color=color, linewidth=5)
        text!(ax_legend, 0.5, i + 0.2, text=label, align=(:center, :center))
    end
    
    return fig
end

"""
    visualize_interface_zoom(root, level_set_func, domain_size; zoom_factor=0.3, num_points=10)

Visualize the quadtree with zoomed view on the interface region.
- `zoom_factor`: Controls how much to zoom out from the interface cells
- `num_points`: Number of quadrature points for geometric fraction calculation
"""
function visualize_interface_zoom(root, level_set_func, domain_size; zoom_factor=0.3, num_points=10)
    leaves = get_leaf_cells(root)
    cell_geometries = compute_geometric_fractions(leaves, level_set_func, num_points=num_points)
    
    # Identify mixed cells (those crossed by the interface)
    mixed_cells = filter(cell -> is_mixed_cell(cell, level_set_func), leaves)
    
    if isempty(mixed_cells)
        error("No mixed cells found to zoom on interface.")
    end
    
    # Determine bounds of the region of interest around mixed cells
    x_min = minimum(cell.x_min for cell in mixed_cells)
    x_max = maximum(cell.x_min + cell.width for cell in mixed_cells)
    y_min = minimum(cell.y_min for cell in mixed_cells)
    y_max = maximum(cell.y_min + cell.height for cell in mixed_cells)
    
    # Add margin around the region
    width = x_max - x_min
    height = y_max - y_min
    x_min -= width * zoom_factor
    x_max += width * zoom_factor
    y_min -= height * zoom_factor
    y_max += height * zoom_factor
    
    # Ensure region is within the domain
    x_min = max(0.0, x_min)
    y_min = max(0.0, y_min)
    x_max = min(domain_size[1], x_max)
    y_max = min(domain_size[2], y_max)
    
    # Create figure with two visualizations: overview and zoom
    fig = Figure(size=(1800, 900))
    
    # Overview (small)
    ax1 = Axis(fig[1, 1], aspect=1, title="Overview")
    
    # Interface
    n_points = 500
    x_range = range(0, domain_size[1], length=n_points)
    y_range = range(0, domain_size[2], length=n_points)
    z = [level_set_func(x, y) for y in y_range, x in x_range]
    contour!(ax1, x_range, y_range, z, levels=[0], color=:black, linewidth=2)
    
    # Rectangle indicating zoom area
    poly!(ax1, Rect(x_min, y_min, x_max-x_min, y_max-y_min), 
          color=(RGBf(1.0, 0.0, 0.0), 0.2), 
          strokewidth=2, strokecolor=:red)
    
    # Colormap for volume fractions (overview)
    for cell in leaves
        geom = cell_geometries[cell]
        rect = Rect(cell.x_min, cell.y_min, cell.width, cell.height)
        
        # Color gradient based on volume fraction
        color = RGBf(0.0, 0.0, geom.volume_fraction)
        alpha = 0.7
        
        poly!(ax1, rect, color=(color, alpha), strokewidth=0.5, strokecolor=:black)
    end
    
    # Zoomed view (large)
    ax2 = Axis(fig[1, 2:3], aspect=1, 
               title="Zoom on interface - Volume fraction",
               xlabel="x", ylabel="y")
    
    # Limit view to region of interest
    limits!(ax2, x_min, x_max, y_min, y_max)
    
    # Interface on zoomed view (with more points for higher precision)
    n_points_zoom = 2000
    x_range_zoom = range(x_min, x_max, length=n_points_zoom)
    y_range_zoom = range(y_min, y_max, length=n_points_zoom)
    z_zoom = [level_set_func(x, y) for y in y_range_zoom, x in x_range_zoom]
    contour!(ax2, x_range_zoom, y_range_zoom, z_zoom, levels=[0], color=:black, linewidth=3)
    
    # Colormap for volume fractions (zoomed with more details)
    for cell in leaves
        # Only draw cells in the zoomed region
        if cell.x_min <= x_max && cell.x_min + cell.width >= x_min && 
           cell.y_min <= y_max && cell.y_min + cell.height >= y_min
           
            geom = cell_geometries[cell]
            rect = Rect(cell.x_min, cell.y_min, cell.width, cell.height)
            
            # Use more visible color scale
            # Blue for fluid (values close to 1)
            # Red for solid (values close to 0)
            color_val = geom.volume_fraction
            
            if color_val > 0.5
                # From white (0.5) to blue (1.0)
                intensity = 2 * (color_val - 0.5)
                color = RGBf(1-intensity, 1-intensity, 1.0)
            else
                # From red (0.0) to white (0.5)
                intensity = 2 * color_val
                color = RGBf(1.0, intensity, intensity)
            end
            
            alpha = 0.7
            
            poly!(ax2, rect, color=(color, alpha), strokewidth=1.5, strokecolor=:black)
            
            # Display numerical values for large enough cells
            # and those near the interface
            if (x_max - x_min) / 20 < cell.width && is_mixed_cell(cell, level_set_func)
                text!(ax2, cell.x_min + cell.width/2, cell.y_min + cell.height/2, 
                      text="$(round(color_val, digits=2))", 
                      align=(:center, :center), 
                      color=color_val > 0.3 ? :black : :white, 
                      fontsize=min(14, 6 + 40*cell.width/(x_max-x_min)))
            end
        end
    end
    
    # Color bar
    colorbar_axis = Axis(fig[1, 4], width=60)
    hidedecorations!(colorbar_axis)
    hidespines!(colorbar_axis)
    
    n_bars = 100
    for i in 1:n_bars
        frac = (i-1)/(n_bars-1)
        y = frac * domain_size[2]
        
        if frac > 0.5
            # From white (0.5) to blue (1.0)
            intensity = 2 * (frac - 0.5)
            color = RGBf(1-intensity, 1-intensity, 1.0)
        else
            # From red (0.0) to white (0.5)
            intensity = 2 * frac
            color = RGBf(1.0, intensity, intensity)
        end
        
        poly!(colorbar_axis, Rect(0, y, 1, domain_size[2]/n_bars), color=color)
    end
    
    # Color bar labels
    text!(colorbar_axis, 0.5, 0, text="0.0 (solid)", align=(:center, :bottom))
    text!(colorbar_axis, 0.5, domain_size[2]/2, text="0.5", align=(:center, :center))
    text!(colorbar_axis, 0.5, domain_size[2], text="1.0 (fluid)", align=(:center, :top))
    
    return fig
end