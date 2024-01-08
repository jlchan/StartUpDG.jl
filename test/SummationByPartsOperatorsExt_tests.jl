
@testset "SummationByPartsOperators ext tests for nonperiodic operators" begin
    tol = 200*eps()
    N = 3
    f(N,r,s) = r^N + s^N # function for testing differentiation

    D = derivative_operator(MattssonNordström2004(), 
                            derivative_order=1, accuracy_order=N+1,
                            xmin=-1.0, xmax=1.0, N=20)
    @testset "Line" begin
        rd = RefElemData(Line(), D)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 2
        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test rd.Pq*rd.Vq ≈ I
        @test all(vec.(rd.rstf) .≈ (x->getindex(x,rd.Fmask)).(rd.rst))
        @test StartUpDG.eigenvalue_inverse_trace_constant(rd) ≈ inverse_trace_constant(rd)    

        md = MeshData(uniform_mesh(Line(), 2), rd)
        @test estimate_h(1, rd, md) ≈ 0.5
    end

    @testset "Quad" begin
        rd = RefElemData(Quad(), D)
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
    
    @testset "Hex" begin
        rd = RefElemData(Hex(), D)
        @test propertynames(rd)[1] == :element_type
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

@testset "SummationByPartsOperators ext tests for periodic operators" begin
    tol = 200*eps()
    N = 3
    D = periodic_derivative_operator(derivative_order=1, accuracy_order=N+1,
                                     xmin=0.0, xmax=2.0, N=20)

    @testset "Line" begin
        rd = RefElemData(Line(), D)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 2

        # diff matrices should be skew symmetric
        Qrst = (A -> rd.M * A).(rd.Drst)
        @test all((A -> norm(A+A')).(Qrst) .< tol)

        @test rd.Pq*rd.Vq ≈ I

    end

    @testset "Quad" begin
        rd = RefElemData(Quad(), D)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 4

        # diff matrices should be skew symmetric
        Qrst = (A -> rd.M * A).(rd.Drst)
        @test all((A -> norm(A+A')).(Qrst) .< tol)

        @test rd.Pq*rd.Vq ≈ I
    end
    
    @testset "Hex" begin
        rd = RefElemData(Hex(), D)
        @test propertynames(rd)[1] == :element_type
        @test rd.t == rd.rst[3]
        @test rd.tq == rd.rstq[3]
        @test rd.rp == rd.rstp[1]    
        @test rd.sp == rd.rstp[2]
        @test rd.tp == rd.rstp[3]    
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 8
        
        # diff matrices should be skew symmetric
        Qrst = (A -> rd.M * A).(rd.Drst)
        @test all((A -> norm(A+A')).(Qrst) .< tol)

        @test rd.Pq*rd.Vq ≈ I
    end
end