@testset "HOHQMesh" begin 

    filename = "testset_HOHQMesh_meshes/easy_example.mesh"
    # filename = "test/testset_HOHQMesh_meshes/easy_example.mesh"
    hmd = read_HOHQMesh(filename)
    rd = RefElemData(Quad(), 4)
    md = MeshData(hmd, rd)
    area = 16 - 0.25^2 * pi
    @test all(md.J .> 0)
    @test abs(sum(md.wJq) - area) < 1e-4

    filename = "testset_HOHQMesh_meshes/GingerbreadMan.mesh"
    # filename = "test/testset_HOHQMesh_meshes/GingerbreadMan.mesh"
    hmd = read_HOHQMesh(filename)
    rd = RefElemData(Quad(), 4)
    md = MeshData(hmd, rd)    
    @test sum(md.wJq) ≈ 1060.558197978162

    # hex meshes
    filename = "testset_HOHQMesh_meshes/MSMappedHex4P4.mesh"
    hmd = read_HOHQMesh(filename)
    # md = MeshData(hmd, rd)    

    # Tri and Tet meshes
    filename = "testset_HOHQMesh_meshes/MSMappedTri4P4.mesh"
    # filename = "test/testset_HOHQMesh_meshes/MSMappedTri4P4.mesh"
    hmd = read_HOHQMesh(filename, Tri())
    rd = RefElemData(Tri(), 4)
    md = MeshData(hmd, rd)
    @test all(md.J .> 0)

    filename = "testset_HOHQMesh_meshes/MSMappedTet4P4.mesh"
    hmd = read_HOHQMesh(filename, Tet())

end
