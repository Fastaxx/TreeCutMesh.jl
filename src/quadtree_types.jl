"""
    ThreadedQuadTreeCell

Structure representing a "fully-threaded" quadtree cell with direct references to neighbors.

# Fields
- `x_min::Float64`: The minimum x-coordinate of the cell
- `y_min::Float64`: The minimum y-coordinate of the cell
- `width::Float64`: The width of the cell
- `height::Float64`: The height of the cell
- `level::Int`: The refinement level of the cell
- `children::Vector{ThreadedQuadTreeCell}`: Child cells (empty if leaf)
- `is_leaf::Bool`: True if the cell is a leaf (no children)
- `parent::Union{Nothing, ThreadedQuadTreeCell}`: Reference to parent cell (nothing for root)
- `north::Union{Nothing, ThreadedQuadTreeCell}`: North neighbor reference
- `south::Union{Nothing, ThreadedQuadTreeCell}`: South neighbor reference
- `east::Union{Nothing, ThreadedQuadTreeCell}`: East neighbor reference
- `west::Union{Nothing, ThreadedQuadTreeCell}`: West neighbor reference
"""
mutable struct ThreadedQuadTreeCell
    x_min::Float64
    y_min::Float64
    width::Float64
    height::Float64
    level::Int
    children::Vector{ThreadedQuadTreeCell}
    is_leaf::Bool
    
    # Reference to parent cell
    parent::Union{Nothing, ThreadedQuadTreeCell}
    
    # Direct references to neighbors (O(1) access)
    north::Union{Nothing, ThreadedQuadTreeCell}  # north neighbor
    south::Union{Nothing, ThreadedQuadTreeCell}  # south neighbor
    east::Union{Nothing, ThreadedQuadTreeCell}   # east neighbor
    west::Union{Nothing, ThreadedQuadTreeCell}   # west neighbor
    
    # Constructor for a cell without children (leaf)
    function ThreadedQuadTreeCell(x_min, y_min, width, height, level)
        return new(x_min, y_min, width, height, level, ThreadedQuadTreeCell[], true, nothing, nothing, nothing, nothing, nothing)
    end
end

"""
    CellGeometry

Structure for storing geometric fractions used in cut cell methods.

# Fields
- `volume_fraction::Float64`: Fraction of cell volume occupied by fluid
- `face_fraction_north::Float64`: Fraction of north face occupied by fluid
- `face_fraction_south::Float64`: Fraction of south face occupied by fluid
- `face_fraction_east::Float64`: Fraction of east face occupied by fluid
- `face_fraction_west::Float64`: Fraction of west face occupied by fluid
"""
mutable struct CellGeometry
    # Volume fraction (% of cell occupied by fluid)
    volume_fraction::Float64
    
    # Surface fractions for each face
    face_fraction_north::Float64
    face_fraction_south::Float64
    face_fraction_east::Float64
    face_fraction_west::Float64
    
    # Default constructor
    function CellGeometry()
        return new(0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

"""
    create_root_cell(x_min, y_min, width, height)

Create a root cell for the quadtree with the specified dimensions.
"""
function create_root_cell(x_min, y_min, width, height)
    root = ThreadedQuadTreeCell(x_min, y_min, width, height, 0)
    # Root has no neighbors at initialization
    return root
end