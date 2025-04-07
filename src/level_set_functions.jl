"""
    level_set_circle(x, y, center_x, center_y, radius)

Level set function for a circle centered at (center_x, center_y) with given radius.
Returns negative values inside the circle, positive values outside, and zero on the circle.
"""
function level_set_circle(x, y, center_x, center_y, radius)
    return sqrt((x - center_x)^2 + (y - center_y)^2) - radius
end

"""
    level_set_flower(x, y, center_x, center_y, radius, petals, amplitude)

Level set function for a flower/star shape centered at (center_x, center_y).
- `radius`: Base radius of the flower
- `petals`: Number of petals/points
- `amplitude`: How pronounced the petals are (0=circle, higher values give sharper points)
Returns negative values inside the shape, positive values outside, and zero on the boundary.
"""
function level_set_flower(x, y, center_x, center_y, radius, petals, amplitude)
    dx = x - center_x
    dy = y - center_y
    r = sqrt(dx^2 + dy^2)
    θ = atan(dy, dx) - deg2rad(9)
    
    # Radius that varies with angle
    r_interface = radius * (1 + amplitude * cos(petals * θ))
    
    return r - r_interface
end

"""
    is_mixed_cell(cell, level_set_func)

Determine if a cell is cut by the interface (mixed cell).
Returns true if the level set function changes sign within the cell.
"""
function is_mixed_cell(cell, level_set_func)
    # Corner points
    corners = [
        (cell.x_min, cell.y_min),                      
        (cell.x_min + cell.width, cell.y_min),         
        (cell.x_min, cell.y_min + cell.height),        
        (cell.x_min + cell.width, cell.y_min + cell.height)   
    ]
    
    # Evaluate level set function at corners
    values = [level_set_func(x, y) for (x, y) in corners]
    
    # Check if interface crosses the cell (change of sign)
    has_positive = any(values .> 0)
    has_negative = any(values .< 0)
    
    return has_positive && has_negative
end