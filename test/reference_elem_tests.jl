# some code not tested to avoid redundancy from tests in NodesAndModes.
@testset "Reference elements" begin
    tol = 5e2*eps()

    N = 2
    @testset "Interval" begin
        rd = RefElemData(Line(),N)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)
        @test abs(sum(rd.rq .* rd.wq)) < tol 
        @test rd.nrJ ≈ [-1,1]
        @test rd.Pq * rd.Vq ≈ I
        @test rd.r[rd.Fmask[:]] ≈ rd.rf
        @test invoke(inverse_trace_constant, Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)
        @test propertynames(rd)[1] == :elementType
        @test rd.elementType == rd.element_type
        @test rd.approximationType == rd.approximation_type
    end

    @testset "Triangle" begin
        rd = RefElemData(Tri(),N)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 2
        @test abs(sum(rd.wf)) ≈ 6
        @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
        @test rd.Pq * rd.Vq ≈ I
        Vfp = vandermonde(Line(), N, quad_nodes(Line(), N)[1]) / vandermonde(Line(), N, nodes(Line(), N))
        rstf = (x->Vfp * x[reshape(rd.Fmask, rd.Nfq ÷ rd.Nfaces, rd.Nfaces)]).(rd.rst)
        @test all(vec.(rstf) .≈ rd.rstf)
        @test invoke(inverse_trace_constant,Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)
        @test propertynames(rd)[1] == :elementType

        @test StartUpDG.num_vertices(Tri()) == 3
        @test StartUpDG.num_faces(Tri()) == 3

    end

    @testset "Quad" begin
        rd = RefElemData(Quad(),N)
        @test rd.r == rd.rst[1]
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 4
        @test abs(sum(rd.wf)) ≈ 8
        @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
        @test rd.Pq * rd.Vq ≈ I
        Vfp = vandermonde(Line(), N, quad_nodes(Line(), N)[1]) / vandermonde(Line(), N, nodes(Line(), N))
        rstf = (x->Vfp * x[reshape(rd.Fmask,rd.Nfq÷rd.Nfaces,rd.Nfaces)]).(rd.rst)
        @test all(vec.(rstf) .≈ rd.rstf)
        @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)    

        @test StartUpDG.num_vertices(Quad()) == 4
        @test StartUpDG.num_faces(Quad()) == 4
    end

    @testset "Hex" begin
        rd = RefElemData(Hex(),N)
        @test propertynames(rd)[1] == :elementType
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
        @test rd.Pq * rd.Vq ≈ I
        @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)
        # TODO: test interpolation of Fmask matches rd.rstf.

        @test StartUpDG.num_vertices(Hex()) == 8
        @test StartUpDG.num_faces(Hex()) == 6
    end

    @testset "Tet" begin
        rd = RefElemData(Tet(),N)
        @test propertynames(rd)[1] == :elementType
        @test rd.t == rd.rst[3]
        @test rd.tf == rd.rstf[3]    
        @test rd.tq == rd.rstq[3]
        @test rd.rp == rd.rstp[1]    
        @test rd.sp == rd.rstp[2]
        @test rd.tp == rd.rstp[3]    
        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)    
        @test abs(sum(rd.wq)) ≈ 4/3
        @test abs(sum(rd.wf)) ≈ 8
        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test abs(sum(rd.wf .* rd.nsJ)) < tol
        @test abs(sum(rd.wf .* rd.ntJ)) < tol
        @test rd.Pq * rd.Vq ≈ I
        @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)

        @test StartUpDG.num_vertices(Tet()) == 4
        @test StartUpDG.num_faces(Tet()) == 4
    end

    @testset "Wedge" begin
        rd = RefElemData(Wedge(), N)
        @test rd.Pq * rd.Vq ≈ I
        @test propertynames(rd)[1] == :elementType
        @test rd.r == rd.rst[1]
        @test rd.s == rd.rst[2]
        @test rd.t == rd.rst[3]

        @test rd.rf == rd.rstf[1]    
        @test rd.sf == rd.rstf[2]
        @test rd.tf == rd.rstf[3] 

        @test rd.rq == rd.rstq[1]    
        @test rd.sq == rd.rstq[2]
        @test rd.tq == rd.rstq[3] 

        @test rd.rp == rd.rstp[1]    
        @test rd.sp == rd.rstp[2]
        @test rd.tp == rd.rstp[3] 

        @test isapprox(rd.rf, rd.Vf * rd.r)
        @test isapprox(rd.rq, rd.Vq * rd.r)

        @test isapprox(rd.sf, rd.Vf * rd.s)
        @test isapprox(rd.sq, rd.Vq * rd.s)

        @test isapprox(rd.tf, rd.Vf * rd.t)
        @test isapprox(rd.tq, rd.Vq * rd.t)

        @test rd.Np == length(rd.r)  
        @test rd.Nq == length(rd.rq)
        
        # 2 * tri_face + 3 * quad_face
        @test abs(sum(rd.wf)) - (2*2 + 3*4) < tol

        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test abs(sum(rd.wf .* rd.nsJ)) < tol
        @test abs(sum(rd.wf .* rd.ntJ)) < tol
        
        @test StartUpDG.num_faces(Wedge()) == 5
        @test StartUpDG.num_vertices(Wedge()) == 6
    end

    @testset "Misc Pyr" begin
        @test StartUpDG.num_faces(Pyr()) == 5
        @test StartUpDG.num_vertices(Pyr()) == 5
    end
end