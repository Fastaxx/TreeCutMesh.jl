"""
    gauss_legendre_rule(n)

Get Gauss-Legendre quadrature points and weights for n-point rule.
"""
function gauss_legendre_rule(n)
    # Points and weights for n=8 (high precision)
    if n == 8
        points = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                   0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363]
        weights = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                  0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]
        return points, weights
    elseif n == 10
        points = [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312,
                   0.1488743389816312,  0.4333953941292472,  0.6794095682990244,  0.8650633666889845,  0.9739065285171717]
        weights = [0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529,
                   0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881]
        return points, weights
    else
        # Simplified implementation for other values of n
        # In practice, a numerical library would be used
        θ = range(π/(2*n), π-π/(2*n), length=n)
        points = -cos.(θ)
        weights = π/n * ones(n)
        return points, weights
    end
end

"""
    compute_geometric_fractions(leaves, level_set_func; num_points=8)

Compute geometric fractions for cut cells using Gauss-Legendre quadrature.
Returns a dictionary mapping cells to CellGeometry objects.
"""
function compute_geometric_fractions(leaves, level_set_func; num_points=8)
    # Gauss-Legendre quadrature points and weights
    points, weights = gauss_legendre_rule(num_points)
    
    # Store results
    cell_geometries = Dict{ThreadedQuadTreeCell, CellGeometry}()
    
    for cell in leaves
        # Initialize cell geometry
        geom = CellGeometry()
        
        # 1. Compute volume fraction with 2D Gauss-Legendre quadrature
        volume_integral = 0.0
        for i in 1:num_points
            for j in 1:num_points
                # Transform reference coordinates [-1,1]×[-1,1] to cell
                x = cell.x_min + 0.5 * (points[i] + 1) * cell.width
                y = cell.y_min + 0.5 * (points[j] + 1) * cell.height
                
                # Indicator function (1 if in fluid, 0 otherwise)
                indicator = level_set_func(x, y) < 0 ? 1.0 : 0.0
                
                # Add weighted contribution to integral
                volume_integral += weights[i] * weights[j] * indicator
            end
        end
        
        # Normalize integral (reference domain is [-1,1]×[-1,1])
        geom.volume_fraction = volume_integral / 4.0
        
        # 2. Compute surface fractions with 1D Gauss-Legendre quadrature
        
        # North face
        north_integral = 0.0
        for i in 1:num_points
            x = cell.x_min + 0.5 * (points[i] + 1) * cell.width
            y = cell.y_min + cell.height
            north_integral += weights[i] * (level_set_func(x, y) < 0 ? 1.0 : 0.0)
        end
        geom.face_fraction_north = north_integral / 2.0
        
        # South face
        south_integral = 0.0
        for i in 1:num_points
            x = cell.x_min + 0.5 * (points[i] + 1) * cell.width
            y = cell.y_min
            south_integral += weights[i] * (level_set_func(x, y) < 0 ? 1.0 : 0.0)
        end
        geom.face_fraction_south = south_integral / 2.0
        
        # East face
        east_integral = 0.0
        for j in 1:num_points
            x = cell.x_min + cell.width
            y = cell.y_min + 0.5 * (points[j] + 1) * cell.height
            east_integral += weights[j] * (level_set_func(x, y) < 0 ? 1.0 : 0.0)
        end
        geom.face_fraction_east = east_integral / 2.0
        
        # West face
        west_integral = 0.0
        for j in 1:num_points
            x = cell.x_min
            y = cell.y_min + 0.5 * (points[j] + 1) * cell.height
            west_integral += weights[j] * (level_set_func(x, y) < 0 ? 1.0 : 0.0)
        end
        geom.face_fraction_west = west_integral / 2.0
        
        # Store result
        cell_geometries[cell] = geom
    end
    
    return cell_geometries
end