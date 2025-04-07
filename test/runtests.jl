using TreeCutMesh
using Test
using LinearAlgebra

@testset "TreeCutMesh.jl" begin
    # Test basic tree structure
    @testset "Quadtree Structure" begin
        # Create a root cell
        root = create_root_cell(0.0, 0.0, 1.0, 1.0)
        @test root.level == 0
        @test root.is_leaf == true
        @test isempty(root.children)
        
        # Simple circle level set
        circle_level_set = (x, y) -> level_set_circle(x, y, 0.5, 0.5, 0.25)
        
        # Refine once
        refine_cell_whitney_with_constraints!(root, circle_level_set, 1, 0.01)
        @test root.is_leaf == false
        @test length(root.children) == 4
        
        # Test child properties
        sw, se, nw, ne = root.children
        @test sw.level == 1
        @test sw.x_min ≈ 0.0 && sw.y_min ≈ 0.0
        @test sw.width ≈ 0.5 && sw.height ≈ 0.5
        
        # Test neighbor connections
        @test sw.east == se && sw.north == nw
        @test se.west == sw && se.north == ne
        @test nw.south == sw && nw.east == ne
        @test ne.south == se && ne.west == nw
    end
    
    @testset "Level Set Functions" begin
        # Circle
        @test level_set_circle(0.5, 0.5, 0.5, 0.5, 0.25) ≈ -0.25
        @test level_set_circle(0.8, 0.5, 0.5, 0.5, 0.25) ≈ 0.05
        
        # Flower
        flower = (x, y) -> level_set_flower(x, y, 0.5, 0.5, 0.25, 5, 0.3)
        # Test a point we know is inside
        @test flower(0.5, 0.5) < 0
        
        # Test mixed cell detection
        cell = ThreadedQuadTreeCell(0.4, 0.4, 0.2, 0.2, 1)
        circle_level_set = (x, y) -> level_set_circle(x, y, 0.5, 0.5, 0.25)
        @test is_mixed_cell(cell, circle_level_set) == false  # This cell should be fully inside
        
        cell_mixed = ThreadedQuadTreeCell(0.3, 0.3, 0.4, 0.4, 1)
        @test is_mixed_cell(cell_mixed, circle_level_set) == true  # This cell should be mixed
    end
    
    @testset "Tree Building" begin
        # Test full tree building
        circle_level_set = (x, y) -> level_set_circle(x, y, 0.5, 0.5, 0.25)
        tree = build_constrained_quadtree(0.0, 0.0, 1.0, 1.0, circle_level_set, max_level=4)
        leaves = get_leaf_cells(tree)
        
        # Test some basic properties
        @test length(leaves) > 0
        @test all(leaf -> leaf.is_leaf, leaves)
        @test all(leaf -> 0 <= leaf.level <= 4, leaves)
        
        # The tree should be balanced
        for leaf in leaves
            neighbors = [leaf.north, leaf.south, leaf.east, leaf.west]
            for n in neighbors
                if !isnothing(n)
                    @test abs(n.level - leaf.level) <= 1
                end
            end
        end
    end
    
    @testset "Geometric Fractions" begin
        # Test computation of geometric fractions
        circle_level_set = (x, y) -> level_set_circle(x, y, 0.5, 0.5, 0.25)
        cell_inside = ThreadedQuadTreeCell(0.4, 0.4, 0.2, 0.2, 1)
        cell_outside = ThreadedQuadTreeCell(0.0, 0.0, 0.2, 0.2, 1)
        cell_mixed = ThreadedQuadTreeCell(0.3, 0.3, 0.4, 0.4, 1)
        
        geoms = compute_geometric_fractions([cell_inside, cell_outside, cell_mixed], circle_level_set)
        
        # Cell entirely inside should have volume fraction ≈ 1
        @test geoms[cell_inside].volume_fraction ≈ 1.0 atol=0.01
        
        # Cell entirely outside should have volume fraction ≈ 0
        @test geoms[cell_outside].volume_fraction ≈ 0.0 atol=0.01
        
        # Mixed cell should have volume fraction between 0 and 1
        @test 0 < geoms[cell_mixed].volume_fraction < 1
        
        # Test Gauss-Legendre rule
        points, weights = TreeCutMesh.gauss_legendre_rule(8)
        @test length(points) == 8
        @test length(weights) == 8
        # Sum of weights should be 2 (for interval [-1,1])
        @test sum(weights) ≈ 2.0 atol=1e-14
    end
end