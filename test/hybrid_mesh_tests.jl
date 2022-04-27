@testset "Hybrid mesh utilities" begin
    # Simple hybrid mesh for testing
    #   1  7_____8_____9
    #      |     |   / |
    #      |     | /   |
    #   0  4 --- 5 --- 6 
    #      |     |     |
    #      |     |     |
    #   -1 1 --- 2 --- 3
    #     
    #     -1     0     1
    VX = [-1; 0; 1; -1; 0; 1; -1; 0; 1]
    VY = [-1; -1; -1; 0; 0; 0; 1; 1; 1]
    EToV = [[1 2 4 5], [2 3 5 6], [4 5 7 8], [2 3 6], [6 5 2]]

    rds = Dict((elem => RefElemData(N=2, elem) for elem in (Tri(), Quad())))
    fvs = Dict(Pair.(keys(rds), getproperty.(values(rds), :fv)))

    FToF = StartUpDG.connect_mesh(EToV, fvs)

    @test FToF == vec([1  5  3  11  17  14  15  18  9  10  4  12  16  6  7  13  2  8])

end