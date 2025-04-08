"""
    update_neighbors_after_refinement!(parent_cell, children)

Update neighbor references after refinement of a cell.
Handles connections between sibling cells and external neighbors.
"""
function update_neighbors_after_refinement!(parent_cell, children)
    sw, se, nw, ne = children  # South-West, South-East, North-West, North-East
    
    # 1. Configure links between sibling cells
    sw.east = se
    sw.north = nw
    
    se.west = sw
    se.north = ne
    
    nw.south = sw
    nw.east = ne
    
    ne.south = se
    ne.west = nw
    
    # 2. Connect with parent's external neighbors
    
    # South neighbor of parent
    if !isnothing(parent_cell.south)
        s_neighbor = parent_cell.south
        if !s_neighbor.is_leaf
            # Find the north children of neighbor that match our south children
            for child in s_neighbor.children
                if child.y_min + child.height == sw.y_min && 
                   child.x_min + child.width > sw.x_min && 
                   child.x_min < sw.x_min + sw.width
                    sw.south = child
                    child.north = sw
                end
                if child.y_min + child.height == se.y_min && 
                   child.x_min + child.width > se.x_min && 
                   child.x_min < se.x_min + se.width
                    se.south = child
                    child.north = se
                end
            end
        else
            # Unrefined neighbor, connect both south children to this neighbor
            sw.south = s_neighbor
            se.south = s_neighbor
        end
    end
    
    # North neighbor of parent
    if !isnothing(parent_cell.north)
        n_neighbor = parent_cell.north
        if !n_neighbor.is_leaf
            for child in n_neighbor.children
                if child.y_min == nw.y_min + nw.height && 
                   child.x_min + child.width > nw.x_min && 
                   child.x_min < nw.x_min + nw.width
                    nw.north = child
                    child.south = nw
                end
                if child.y_min == ne.y_min + ne.height && 
                   child.x_min + child.width > ne.x_min && 
                   child.x_min < ne.x_min + ne.width
                    ne.north = child
                    child.south = ne
                end
            end
        else
            nw.north = n_neighbor
            ne.north = n_neighbor
        end
    end
    
    # West neighbor of parent
    if !isnothing(parent_cell.west)
        w_neighbor = parent_cell.west
        if !w_neighbor.is_leaf
            for child in w_neighbor.children
                if child.x_min + child.width == sw.x_min && 
                   child.y_min + child.height > sw.y_min && 
                   child.y_min < sw.y_min + sw.height
                    sw.west = child
                    child.east = sw
                end
                if child.x_min + child.width == nw.x_min && 
                   child.y_min + child.height > nw.y_min && 
                   child.y_min < nw.y_min + nw.height
                    nw.west = child
                    child.east = nw
                end
            end
        else
            sw.west = w_neighbor
            nw.west = w_neighbor
        end
    end
    
    # East neighbor of parent
    if !isnothing(parent_cell.east)
        e_neighbor = parent_cell.east
        if !e_neighbor.is_leaf
            for child in e_neighbor.children
                if child.x_min == se.x_min + se.width && 
                   child.y_min + child.height > se.y_min && 
                   child.y_min < se.y_min + se.height
                    se.east = child
                    child.west = se
                end
                if child.x_min == ne.x_min + ne.width && 
                   child.y_min + child.height > ne.y_min && 
                   child.y_min < ne.y_min + ne.height
                    ne.east = child
                    child.west = ne
                end
            end
        else
            se.east = e_neighbor
            ne.east = e_neighbor
        end
    end
end

"""
    refine_cell_whitney_with_constraints!(cell, level_set_func, max_level, min_cell_size, lip_const=1.0)

Refine a cell according to Whitney's criterion and other constraints.
- `level_set_func`: Function defining the interface
- `max_level`: Maximum refinement level allowed
- `min_cell_size`: Minimum cell size allowed
- `lip_const`: Lipschitz constant for Whitney criterion
"""
function refine_cell_whitney_with_constraints!(cell, level_set_func, max_level, min_cell_size, lip_const=1.0)
    # Check if cell is already refined enough or too small
    if cell.level >= max_level || cell.width < min_cell_size || cell.height < min_cell_size
        return
    end
    
    # Calculate center coordinates
    x_center = cell.x_min + cell.width/2
    y_center = cell.y_min + cell.height/2
    
    # Dense sampling to better detect interface
    vertices = [
        (cell.x_min, cell.y_min),
        (cell.x_min + cell.width, cell.y_min),
        (cell.x_min, cell.y_min + cell.height),
        (cell.x_min + cell.width, cell.y_min + cell.height),
        (x_center, cell.y_min),
        (x_center, cell.y_min + cell.height),
        (cell.x_min, y_center),
        (cell.x_min + cell.width, y_center),
        (x_center, y_center)
    ]
    
    # Add intermediate points for complex shapes
    if cell.level <= 4
        for i in [0.25, 0.5, 0.75]
            for j in [0.25, 0.5, 0.75]
                push!(vertices, (cell.x_min + i*cell.width, cell.y_min + j*cell.height))
            end
        end
    end
    
    # Calculate level set values and absolute values
    values = [level_set_func(x, y) for (x, y) in vertices]
    abs_values = [abs(val) for val in values]
    
    has_positive = any(values .> 0)
    has_negative = any(values .< 0)
    is_interface_cell = has_positive && has_negative
    
    # Calculate cell diagonal
    diag_size = sqrt(cell.width^2 + cell.height^2)
    
    min_abs_value = minimum(abs_values)
    needs_refinement = min_abs_value <= lip_const * diag_size
    
    # Refine if necessary
    if needs_refinement || is_interface_cell
        cell.is_leaf = false
        
        half_width = cell.width / 2
        half_height = cell.height / 2
        
        # Create children with parent information
        sw = ThreadedQuadTreeCell(cell.x_min, cell.y_min, half_width, half_height, cell.level + 1)
        se = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min, half_width, half_height, cell.level + 1)
        nw = ThreadedQuadTreeCell(cell.x_min, cell.y_min + half_height, half_width, half_height, cell.level + 1)
        ne = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min + half_height, half_width, half_height, cell.level + 1)
        
        # Set parentage
        sw.parent = cell
        se.parent = cell
        nw.parent = cell
        ne.parent = cell
        
        # Add children to cell
        push!(cell.children, sw)
        push!(cell.children, se)
        push!(cell.children, nw)
        push!(cell.children, ne)
        
        # Update neighbor relationships
        update_neighbors_after_refinement!(cell, [sw, se, nw, ne])
        
        # Recursively refine children
        for child in cell.children
            refine_cell_whitney_with_constraints!(child, level_set_func, max_level, min_cell_size, lip_const)
        end
    end
end

"""
    get_leaf_cells(root)

Get all leaf cells in the quadtree in O(N) time without search.
"""
function get_leaf_cells(root)
    leaves = ThreadedQuadTreeCell[]
    stack = [root]
    
    # Iterative version more efficient than recursive
    while !isempty(stack)
        cell = pop!(stack)
        
        if cell.is_leaf
            push!(leaves, cell)
        else
            # Add all children to the stack
            append!(stack, cell.children)
        end
    end
    
    return leaves
end

"""
    are_neighbors(cell1, cell2)

Check if two cells are neighbors.
With the fully-threaded structure, this is an O(1) operation.
"""
function are_neighbors(cell1, cell2)
    return cell1.north == cell2 || cell1.south == cell2 || 
           cell1.east == cell2 || cell1.west == cell2
end
"""
    balance_quadtree!(root, max_level, min_cell_size)

Balance the quadtree so that neighboring cells differ by at most one level.
Uses efficient O(N) traversal with direct neighbor references.
"""
function balance_quadtree!(root, max_level, min_cell_size)
    changes_made = true
    
    while changes_made
        changes_made = false
        leaves = get_leaf_cells(root)
        cells_to_refine = ThreadedQuadTreeCell[]
        
        # Check all leaf cells against their neighbors
        for leaf in leaves
            # Check all 4 directions
            neighbors = [leaf.north, leaf.south, leaf.east, leaf.west]
            for neighbor in neighbors
                if !isnothing(neighbor)
                    # Case 1: Neighbor is a leaf node
                    if neighbor.is_leaf
                        # Check for difference > 1 in either direction
                        if abs(leaf.level - neighbor.level) > 1
                            # Refine the one with lower level
                            if leaf.level < neighbor.level
                                # Current leaf needs refinement
                                if leaf.level < max_level && leaf.width > min_cell_size && leaf.height > min_cell_size
                                    push!(cells_to_refine, leaf)
                                    changes_made = true
                                    break  # No need to check other neighbors of this leaf
                                end
                            else
                                # Neighbor needs refinement
                                if neighbor.level < max_level && neighbor.width > min_cell_size && 
                                   neighbor.height > min_cell_size
                                    push!(cells_to_refine, neighbor)
                                    changes_made = true
                                end
                            end
                        end
                    else
                        # Case 2: Neighbor is not a leaf - check its children that border this cell
                        for child in neighbor.children
                            # Only check children that are adjacent to the current leaf
                            if are_adjacent(leaf, child)
                                if child.is_leaf && abs(leaf.level - child.level) > 1
                                    if leaf.level < child.level
                                        # Current leaf needs refinement
                                        if leaf.level < max_level && leaf.width > min_cell_size && 
                                           leaf.height > min_cell_size
                                            push!(cells_to_refine, leaf)
                                            changes_made = true
                                            break  # No need to check other children
                                        end
                                    end
                                    # We don't need to handle the case where the child needs refinement,
                                    # as we'll catch that in a later iteration
                                end
                            end
                        end
                        if changes_made
                            break  # No need to check other neighbors of this leaf
                        end
                    end
                end
            end
        end
        
        # Refine identified cells 
        for cell in cells_to_refine
            if cell.is_leaf && cell.level < max_level && 
               cell.width > min_cell_size && cell.height > min_cell_size
                
                half_width = cell.width / 2
                half_height = cell.height / 2
                
                sw = ThreadedQuadTreeCell(cell.x_min, cell.y_min, half_width, half_height, cell.level + 1)
                se = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min, half_width, half_height, cell.level + 1)
                nw = ThreadedQuadTreeCell(cell.x_min, cell.y_min + half_height, half_width, half_height, cell.level + 1)
                ne = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min + half_height, half_width, half_height, cell.level + 1)
                
                # Set parentage
                sw.parent = cell
                se.parent = cell
                nw.parent = cell
                ne.parent = cell
                
                cell.is_leaf = false
                cell.children = [sw, se, nw, ne]
                
                # Update neighbor relationships
                update_neighbors_after_refinement!(cell, [sw, se, nw, ne])
            end
        end
    end
    
    return root
end

"""
    are_adjacent(cell1, cell2)

Helper function to determine if two cells share a boundary.
"""
function are_adjacent(cell1, cell2)
    # Check if cells touch on north-south boundary
    ns_adjacent = (abs(cell1.y_min + cell1.height - cell2.y_min) < 1e-10) || 
                 (abs(cell2.y_min + cell2.height - cell1.y_min) < 1e-10)
                 
    # Check if cells overlap in x direction when adjacent in y
    ns_overlap = (cell1.x_min < cell2.x_min + cell2.width) && 
                (cell2.x_min < cell1.x_min + cell1.width)
    
    # Check if cells touch on east-west boundary
    ew_adjacent = (abs(cell1.x_min + cell1.width - cell2.x_min) < 1e-10) || 
                 (abs(cell2.x_min + cell2.width - cell1.x_min) < 1e-10)
                 
    # Check if cells overlap in y direction when adjacent in x
    ew_overlap = (cell1.y_min < cell2.y_min + cell2.height) && 
                (cell2.y_min < cell1.y_min + cell1.height)
    
    return (ns_adjacent && ns_overlap) || (ew_adjacent && ew_overlap)
end

"""
    enforce_interface_level_constraint!(root, level_set_func, max_level, min_cell_size)

Ensure that all interface cells are at the same refinement level.
Returns true if any changes were made to the tree.
"""
function enforce_interface_level_constraint!(root, level_set_func, max_level, min_cell_size)
    leaves = get_leaf_cells(root)
    interface_cells = filter(cell -> is_mixed_cell(cell, level_set_func), leaves)
    
    if isempty(interface_cells)
        return false
    end
    
    # Determine maximum level of interface cells
    max_interface_level = maximum([cell.level for cell in interface_cells])
    
    # Find interface cells not at the maximum level
    cells_to_refine = filter(cell -> cell.level < max_interface_level, interface_cells)
    
    if isempty(cells_to_refine)
        return false
    end
    
    # Refine these cells to the maximum level
    changes_made = false
    for cell in cells_to_refine
        if cell.level < max_interface_level && cell.width > min_cell_size && cell.height > min_cell_size
            cell.is_leaf = false
            changes_made = true
            
            # Create and configure children with their neighbor relationships
            half_width = cell.width / 2
            half_height = cell.height / 2
            
            sw = ThreadedQuadTreeCell(cell.x_min, cell.y_min, half_width, half_height, cell.level + 1)
            se = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min, half_width, half_height, cell.level + 1)
            nw = ThreadedQuadTreeCell(cell.x_min, cell.y_min + half_height, half_width, half_height, cell.level + 1)
            ne = ThreadedQuadTreeCell(cell.x_min + half_width, cell.y_min + half_height, half_width, half_height, cell.level + 1)
            
            # Set parentage
            sw.parent = cell
            se.parent = cell
            nw.parent = cell
            ne.parent = cell
            
            cell.children = [sw, se, nw, ne]
            
            # Update neighbor relationships
            update_neighbors_after_refinement!(cell, [sw, se, nw, ne])
            
            # Recursively refine children if necessary
            for child in cell.children
                if child.level < max_interface_level
                    refine_cell_whitney_with_constraints!(child, level_set_func, max_level, min_cell_size, 1.0)
                end
            end
        end
    end
    
    return changes_made
end

"""
    build_constrained_quadtree(x_min, y_min, width, height, level_set_func; 
        max_level=8, min_cell_size=0.001, lip_const=1.0)

Build a complete constrained quadtree with Whitney refinement, balancing, and interface constraints.
"""
function build_constrained_quadtree(x_min, y_min, width, height, level_set_func; 
                                   max_level=8, min_cell_size=0.001, lip_const=1.0)
    # Create root cell
    root = create_root_cell(x_min, y_min, width, height)
    
    # First pass: refine with Whitney's criterion
    refine_cell_whitney_with_constraints!(root, level_set_func, max_level, min_cell_size, lip_const)
    
    # Second pass: balance the quadtree (O(N) performance)
    balance_quadtree!(root, max_level, min_cell_size)
    
    # Third pass: apply interface level constraint
    changes_made = true
    while changes_made
        changes_made = enforce_interface_level_constraint!(root, level_set_func, max_level, min_cell_size)
        if changes_made
            balance_quadtree!(root, max_level, min_cell_size)
        end
    end
    
    return root
end