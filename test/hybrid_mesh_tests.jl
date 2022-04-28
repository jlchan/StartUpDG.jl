@testset "Hybrid mesh utilities" begin
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

    rds = Dict((elem => RefElemData(N=2, elem) for elem in (Tri(), Quad())))

    fvs = Dict(Pair.(keys(rds), getproperty.(values(rds), :fv)))
    FToF = StartUpDG.connect_mesh(EToV, fvs)

    @test FToF == vec([1  5  3  14  2  6  7  17  16  10  13  12  11  4  15  9  8  18])

end