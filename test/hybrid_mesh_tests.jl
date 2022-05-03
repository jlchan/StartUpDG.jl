@testset "Hybrid mesh utilities" begin
    A, B = randn(5, 5), randn(5, 5)
    x, y = randn(5), randn(5)
    AB = ArrayPartition(A, B)
    xy = ArrayPartition(x, y)
    @test norm(AB * xy - ArrayPartition(A * x, B * y)) < 10 * eps()

    u = (1:3, 3:-1:1)
    v = (3:-1:1, 1:3)
    @test StartUpDG.match_coordinate_vectors(u, v) == [3, 2, 1]
end

@testset "Hybrid mesh RefElemData and MeshData" begin

    rds = RefElemData((Tri(), Quad()), N = 3)
    @test rds[Tri()].element_type == Tri()
    @test rds[Quad()].element_type == Quad()

    #   1  3_______4______5
    #      |   4   |  6  / 
    #      |1     2|7  5
    #      |   3   |  /  
    #  -1  1 ----- 2        1

    VX = [-1; 0; -1; 0; 1]
    VY = [-1; -1; 1; 1; 1]
    EToV = [[1 2 3 4], [2 4 5]]
    md = MeshData(VX, VY, EToV, rds)
    @test md.FToF == vec([1  7  3  4  5  6  2])

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
    @test md.FToF == vec([1, 5, 3, 11, 2, 6, 7, 17, 9, 15, 4, 12, 16, 14, 10, 13, 8, 18])

    # test if all nodes on boundary are ±1
    @test all(@. abs(max(abs(md.xf[md.mapB]), abs(md.yf[md.mapB])) - 1) < 100 * eps() )

    ## test that the DG derivative of a polynomial recovers the exact derivative
    @unpack x, y = md
    u = @. x^3 - x^2 * y + 2 * y^3
    dudx = @. 3 * x^2 - 2 * x * y

    # compute derivatives
    @unpack rxJ, sxJ, J = md
    dudr = ArrayPartition(getproperty.(values(rds), :Dr)...) * u
    duds = ArrayPartition(getproperty.(values(rds), :Ds)...) * u
    @test norm(@. dudx - (rxJ * dudr + sxJ * duds) / J) < 10 * length(u) * eps()

    # compute jumps
    @unpack mapP = md
    uf = ArrayPartition(getproperty.(values(rds), :Vf)...) * u 
    uP = uf[mapP]    
    u_jump = similar(uf)
    u_jump .= uP - uf
    @test norm(u_jump) < 10 * length(uf) * eps()
    if norm(u_jump) > 10 * length(uf) * eps()
        @show size.(u.x)
        @show size.(uf.x)
    end
end