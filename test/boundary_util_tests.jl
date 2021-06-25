@testset "Boundary utilities" begin
    rd = RefElemData(Tri(), N=3)
    md = MeshData(uniform_mesh(Tri(), 2)..., rd)
    on_bottom_boundary(x, y, tol = 1e-13) = abs(y+1) < tol
    on_top_boundary(x, y, tol = 1e-13) = abs(y-1) < tol

    boundary_dict = tag_boundary_faces(md, Dict(:bottom => on_bottom_boundary, :top => on_top_boundary))
    @test boundary_dict == Dict(:bottom=>[1],:top=>[4])

    boundary_dict = tag_boundary_faces(md, nothing)
    @test boundary_dict == Dict(:entire_boundary => [1,2,4,5])

    # test named tuple version
    boundary_dict = tag_boundary_faces(md,(; :bottom => on_bottom_boundary, :top => on_top_boundary))
    @test boundary_dict == (; :bottom=>[1],:top=>[4])

end