@testset "Multidimensional SBP tests" begin

    tol = 200*eps()

    N = 3
    f(N,r,s) = r^N + s^N # function for testing differentiation
    
    @testset "Hicken's 2N triangular SBP operators" begin
        rd = RefElemData(Tri(), SBP{Hicken}(), N)
        (; rq,sq,wq  ) = rd
        rq2,sq2,wq2 = quad_nodes(Tri(),2*N)        
        @test sum(wq.*f.(2*N,rq,sq)) ≈ sum(wq2.*f.(2*N,rq2,sq2))
        @test rd.Dr*rd.r.^N ≈ N*rd.r.^(N-1)
        @test rd.Ds*rd.s.^N ≈ N*rd.s.^(N-1)
        @test norm(rd.Dr*rd.s + rd.Ds*rd.r) < tol
        @test StartUpDG.eigenvalue_inverse_trace_constant(rd) ≈ inverse_trace_constant(rd)    
    end

    @testset "Kubatko's Legendre face node triangular SBP operators" begin
        rd = RefElemData(Tri(), SBP{Kubatko{LegendreFaceNodes}}(), N)
        (; rq,sq,wq  ) = rd
        rq2,sq2,wq2 = quad_nodes(Tri(),2*N-1)
        @test sum(wq.*f.(2*N-1,rq,sq)) ≈ sum(wq2.*f.(2*N-1,rq2,sq2))
        @test rd.Dr*rd.r.^N ≈ N*rd.r.^(N-1)
        @test rd.Ds*rd.s.^N ≈ N*rd.s.^(N-1)
        @test norm(rd.Dr*rd.s + rd.Ds*rd.r) < tol
        @test StartUpDG.eigenvalue_inverse_trace_constant(rd) ≈ inverse_trace_constant(rd)    
    end
    
    @testset "Warning for N=6 Kubatko Lobatto SBP nodes" begin
        @test_logs (:warn,"N=6 SBP operators with quadrature strength 2N-1 and Lobatto face nodes may require very small timesteps.") RefElemData(Tri(),SBP{Kubatko{LobattoFaceNodes}}(),6)
    end 

    @testset "TensorProductLobatto quad" begin
        rd = RefElemData(Quad(),SBP(),N)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 4
        @test abs(sum(rd.wf)) ≈ 8
        @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
        @test rd.Pq*rd.Vq ≈ I
        @test all(vec.(rd.rstf) .≈ (x->getindex(x,rd.Fmask)).(rd.rst))
        @test StartUpDG.eigenvalue_inverse_trace_constant(rd) ≈ inverse_trace_constant(rd)    
    end
    
    @testset "TensorProductLobatto Hex" begin
        rd = RefElemData(Hex(),SBP(),N)
        @test propertynames(rd)[1] == :element_type
        @test rd.rf == rd.r[rd.Fmask]
        @test rd.sf == rd.s[rd.Fmask]
        @test rd.tf == rd.t[rd.Fmask]
        @test rd.t == rd.rst[3]
        @test rd.tf == rd.rstf[3]    
        @test rd.tq == rd.rstq[3]
        @test rd.rp == rd.rstp[1]    
        @test rd.sp == rd.rstp[2]
        @test rd.tp == rd.rstp[3]    
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 8
        @test abs(sum(rd.wf)) ≈ 6*4
        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test abs(sum(rd.wf .* rd.nsJ)) < tol
        @test abs(sum(rd.wf .* rd.ntJ)) < tol
        @test rd.Pq*rd.Vq ≈ I
        @test StartUpDG.eigenvalue_inverse_trace_constant(rd) ≈ inverse_trace_constant(rd)
    end

end