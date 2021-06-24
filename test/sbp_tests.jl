@testset "Non-default SBP elements" begin

    tol = 200*eps()

    N = 3
    f(N,r,s) = r^N + s^N # function for testing 
    
    @testset "Hicken's 2N triangular SBP operators" begin
        rd = RefElemData(Tri(),SBP(), N, quadrature_strength=2*N)
        @unpack rq,sq,wq = rd
        rq2,sq2,wq2 = quad_nodes(Tri(),2*N)        
        @test sum(wq.*f.(2*N,rq,sq)) ≈ sum(wq2.*f.(2*N,rq2,sq2))
        @test rd.Dr*rd.r.^N ≈ N*rd.r.^(N-1)
        @test rd.Ds*rd.s.^N ≈ N*rd.s.^(N-1)
        @test norm(rd.Dr*rd.s + rd.Ds*rd.r) < tol
    end

    @testset "Kubatko's Legendre face node triangular SBP operators" begin
        rd = RefElemData(Tri(), SBP(), N, quadrature_strength=2*N-1, quad_rule_face=:Legendre)
        @unpack rq,sq,wq = rd
        rq2,sq2,wq2 = quad_nodes(Tri(),2*N-1)
        @test sum(wq.*f.(2*N-1,rq,sq)) ≈ sum(wq2.*f.(2*N-1,rq2,sq2))
        @test rd.Dr*rd.r.^N ≈ N*rd.r.^(N-1)
        @test rd.Ds*rd.s.^N ≈ N*rd.s.^(N-1)
        @test norm(rd.Dr*rd.s + rd.Ds*rd.r) < tol
    end
    
    @testset "Warning for N=6 Kubatko Lobatto SBP nodes" begin
        @test_logs (:warn,"N=6 SBP operators with quadrature strength 2N-1 and Lobatto face nodes may require very small timesteps.") RefElemData(Tri(),SBP(),6,quadrature_strength=11)
        @test_throws SystemError RefElemData(Tri(), SBP(), 7,quadrature_strength=14)
    end 

end