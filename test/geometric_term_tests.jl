@testset "Geometric terms" begin
    @testset "Element $elem" for elem in [Tri() Quad() Hex() Tet()]
        tol = 5e2 * eps()
        N = 3
        rd = RefElemData(elem, N)
        geofacs = geometric_factors(rd.rst..., rd.Drst...)
        if ndims(elem) == 2
            rxJ, sxJ, ryJ, syJ, J = geofacs    
            @test all(rxJ .≈ 1)
            @test norm(sxJ) < tol
            @test norm(ryJ) < tol
            @test all(syJ .≈ 1)
            @test all(J .≈ 1)
        elseif ndims(elem) == 3
            rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = geofacs
            @test all(rxJ .≈ 1)
            @test norm(sxJ) < tol
            @test norm(txJ) < tol
            @test norm(ryJ) < tol
            @test all(syJ .≈ 1)
            @test norm(tyJ) < tol
            @test norm(rzJ) < tol
            @test norm(szJ) < tol
            @test all(tzJ .≈ 1)
            @test all(J .≈ 1)
        end
    end
end