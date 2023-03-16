@testset "HOHQMesh" begin 
    rd = RefElemData(Quad(), 3)

    filename = "testset_HOHQMesh_meshes/easy_example.mesh"
    hmd = read_HOHQMesh(filename)
    md = MeshData(hmd, rd)
    area = 16 - 0.25^2 * pi
    @test all(md.J .> 0)
    @test abs(sum(md.wJq) - area) < 1e-4

    filename = "testset_HOHQMesh_meshes/GingerbreadMan.mesh"
    hmd = read_HOHQMesh(filename)
    md = MeshData(hmd, rd)    
    @test sum(md.wJq) ≈ 1075.1016547806466
end
