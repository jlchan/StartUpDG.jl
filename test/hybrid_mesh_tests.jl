@testset "Hybrid mesh utilities" begin
    # Simple hybrid mesh for testing
    #   1  4_____5_____6
    #      |     |   / |
    #      |     | /   |
    #   0  1 --- 2 --- 3 
    #     -1     0     1
    VX = [-1; 0; 1; -1; 0; 1]
    VY = [0; 0; 0; 1; 1; 1]
    EToV = [[1 2 5 4], [2 3 6], [6 5 2]]

    # rd = (RefElemData(N=2, Tri()), RefElemData(N=2, Quad()))
    # fvs = Dict(Pair.(getproperty.(rd, :elementType), getproperty.(rd, :fv)))
    rd = Dict((elem => RefElemData(N=2, elem) for elem in (Tri(), Quad())))
    fvs = Dict(Pair.(keys(rd), getproperty.(values(rd), :fv)))

    @test_nowarn_debug FToF = StartUpDG.connect_mesh(EToV, fvs)
end