using PathIntersections

@testset "Cut meshes" begin
    
    cells_per_dimension = 4
    cells_per_dimension_x, cells_per_dimension_y = cells_per_dimension, cells_per_dimension
    circle = PresetGeometries.Circle(R=0.33, x0=0, y0=0)

    rd = RefElemData(Quad(), N=3; quad_rule_face=gauss_quad(0, 0, 10))
    md = MeshData(rd, (circle, ), cells_per_dimension_x, cells_per_dimension_y)

    A = 4 - pi * .33^2
    @test sum(md.wJq) ≈ A

end