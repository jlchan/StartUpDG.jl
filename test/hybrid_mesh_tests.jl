@testset "Hybrid mesh utilities" begin

    #rds = Dict((elem => RefElemData(N=2, elem) for elem in (Tri(), Quad())))
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
    #     
    #     -1     0     1
    VX = [-1; 0; 1; -1; 0; 1; -1; 0; 1]
    VY = [-1; -1; -1; 0; 0; 0; 1; 1; 1]
    EToV = [[1 2 4 5], [2 3 5 6], [5 8 9], [4 5 7 8], [9 6 5]]

    md = MeshData(VX, VY, EToV, rds)
    @test md.FToF == vec([1, 5, 3, 11, 2, 6, 7, 17, 9, 15, 4, 12, 16, 14, 10, 13, 8, 18])

    # test if all nodes on boundary are ±1
    @test all(@. abs(max(abs(md.xf[md.mapB]), abs(md.yf[md.mapB])) - 1) < 100 * eps() )

end