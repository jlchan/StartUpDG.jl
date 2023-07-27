@testset "HOHQMesh" begin 

    # Quad meshes
    filename = "testset_HOHQMesh_meshes/easy_example.mesh"
    hmd = read_HOHQMesh(filename)
    rd = RefElemData(Quad(), 4)
    md = MeshData(hmd, rd)
    area = 16 - 0.25^2 * pi
    @test all(md.J .> 0)
    @test abs(sum(md.wJq) - area) < 1e-4
    @test md.mesh_type isa StartUpDG.HOHQMeshType
    @test keys(md.mesh_type.boundary_faces) == (:Bottom, :Right, :Top, :circle, :Left)

    filename = "testset_HOHQMesh_meshes/trixi_example.mesh"
    hmd = read_HOHQMesh(filename)
    md = MeshData(hmd, rd)
    @test all(md.J .> 0)

    filename = "testset_HOHQMesh_meshes/GingerbreadMan.mesh"
    hmd = read_HOHQMesh(filename)
    rd = RefElemData(Quad(), 4)
    md = MeshData(hmd, rd)   
    @test all(md.J .> 0) 
    @test sum(md.wJq) ≈ 1278.9397481563008

    # Hex meshes
    filename = "testset_HOHQMesh_meshes/MSMappedHex4P4.mesh"
    hmd = read_HOHQMesh(filename)
    rd = RefElemData(Hex(), 4)
    md = MeshData(hmd, rd)    
    @test all(md.J .> 0)
    @test sum(md.wJq) ≈ 1.0

    # Tri meshes
    filename = "testset_HOHQMesh_meshes/MSMappedTri4P4.mesh"
    hmd = read_HOHQMesh(filename, Tri())
    rd = RefElemData(Tri(), 4)
    md = MeshData(hmd, rd)
    @test all(md.J .> 0)
    @test sum(md.wJq) ≈ 1.0

    # Tet meshes
    filename = "testset_HOHQMesh_meshes/TetMesh44.mesh"
    @test_nowarn hmd = read_HOHQMesh(filename, Tet())
    rd = RefElemData(Tet(), 4)
    md = MeshData(hmd, rd)
    @test all(md.J .> 0)
end
