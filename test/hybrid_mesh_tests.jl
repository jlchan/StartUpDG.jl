@testset "Hybrid meshes" begin
    @testset "Hybrid mesh utilities" begin
        u = (collect(1:3), collect(3:-1:1))
        v = (collect(3:-1:1), collect(1:3))
        @test StartUpDG.match_coordinate_vectors(u, v) == [3, 2, 1]
    end

    @testset "Hybrid mesh RefElemData and MeshData" begin

        rds = RefElemData((Tri(), Quad()), N = 3)
        @test rds.Tri.element_type == Tri()
        @test rds.Quad.element_type == Quad()
        out = (@capture_out Base.show(stdout, MIME"text/plain"(), rds))
        @test out[1:19] == "MultipleRefElemData"

        # Simple hybrid mesh for testing
        #   1  7______8______9
        #      |      | 3  / |
        #      |   4  |  / 5 |
        #   0  4 ---- 5 ---- 6 
        #      |      |      |
        #      |   1  |   2  |
        #   -1 1 ---- 2 ---- 3
        #     -1      0      1
        VX = [-1; 0; 1; -1; 0; 1; -1; 0; 1]
        VY = [-1; -1; -1; 0; 0; 0; 1; 1; 1]
        EToV = [[1 2 4 5], [2 3 5 6], [5 8 9], [4 5 7 8], [9 6 5]]

        md = MeshData(VX, VY, EToV, rds)
        @test typeof(md.mesh_type) <: typeof(StartUpDG.HybridMesh((Quad(), Tri()), (VX, VY), EToV))
        @test md.VX == md.mesh_type.VXYZ[1]
        @test md.VY == md.mesh_type.VXYZ[2]

        # test if all nodes on boundary are ±1
        @test all(@. abs(max(abs(md.xf[md.mapB]), abs(md.yf[md.mapB])) - 1) < 100 * eps() )

        ## test that the DG derivative of a polynomial recovers the exact derivative
        (; x, y  ) = md
        u = @. x^3 - x^2 * y + 2 * y^3
        dudx = @. 3 * x^2 - 2 * x * y

        (; rxJ, sxJ, J  ) = md
        dudr, duds = similar(md.x), similar(md.x)
        dudr.Quad .= rds.Quad.Dr * u.Quad
        duds.Quad .= rds.Quad.Ds * u.Quad
        dudr.Tri .= rds.Tri.Dr * u.Tri
        duds.Tri .= rds.Tri.Ds * u.Tri

        @test norm(@. dudx - (rxJ * dudr + sxJ * duds) / J) < 1e3 * eps()

        # compute jumps
        (; mapP  ) = md
        uf = NamedArrayPartition(Tri=rds.Tri.Vf * u.Tri, Quad=rds.Quad.Vf * u.Quad)
        uP = uf[mapP]    
        u_jump = similar(uf)
        u_jump .= uP - uf
        @test mapP !== md.mapM # make sure connectivity maps aren't just the same
        @test norm(u_jump) < 10 * length(uf) * eps() # jumps should be zero for a continuous function
    end

    @testset "Curved hybrid meshes" begin
        rds = RefElemData((Tri(), Quad()), N = 2)
        md = MeshData(HybridMeshExample(), rds)

        (; x, y  ) = md
        @. x = x + 0.1 * sin(pi * x) * sin(pi * y)
        @. y = y + 0.1 * sin(pi * x) * sin(pi * y)
        md = MeshData(rds, md, x, y)

        # test differentiation of a linear function 
        u = @. 2 * x + y
        dudr, duds = similar(u), similar(u)
        dudr.Tri .= rds.Tri.Dr * u.Tri
        duds.Tri .= rds.Tri.Ds * u.Tri
        dudr.Quad .= rds.Quad.Dr * u.Quad
        duds.Quad .= rds.Quad.Ds * u.Quad
        dudx = @. (md.rxJ * dudr + md.sxJ * duds) / md.J
        dudy = @. (md.ryJ * dudr + md.syJ * duds) / md.J

        @test all(@. dudx ≈ 2.0)
        @test all(@. dudy ≈ 1.0)
    end
end