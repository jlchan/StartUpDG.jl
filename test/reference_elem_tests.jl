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
        @suppress begin # suppress warnings
            @test invoke(inverse_trace_constant, Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)
        end
        @test propertynames(rd)[1] == :element_type

        # test for deprecated CamlCase approximationType usage
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
        @suppress begin # suppress warnings
            @test invoke(inverse_trace_constant,Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)
        end
        @test propertynames(rd)[1] == :element_type

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
        @suppress begin # suppress warnings
            @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)    
        end

        @test StartUpDG.num_vertices(Quad()) == 4
        @test StartUpDG.num_faces(Quad()) == 4
    end

    @testset "Hex" begin
        rd = RefElemData(Hex(),N)
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
        @test rd.Pq * rd.Vq ≈ I
        @suppress begin # suppress warnings
            @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)
        end
        # TODO: test interpolation of Fmask matches rd.rstf.

        @test StartUpDG.num_vertices(Hex()) == 8
        @test StartUpDG.num_faces(Hex()) == 6
    end

    @testset "Tet" begin
        rd = RefElemData(Tet(),N)
        @test propertynames(rd)[1] == :element_type
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
        @suppress begin # suppress warnings
            @test invoke(inverse_trace_constant, Tuple{RefElemData}, rd) ≈ inverse_trace_constant(rd)
        end

        @test StartUpDG.num_vertices(Tet()) == 4
        @test StartUpDG.num_faces(Tet()) == 4

        @test StartUpDG.face_type(Tet(), 3) == StartUpDG.face_type(Tet())
    end

    @testset "Wedge" begin
        rd = RefElemData(Wedge(), N)
        
        @test propertynames(rd)[1] == :element_type
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
        
        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test abs(sum(rd.wf .* rd.nsJ)) < tol
        @test abs(sum(rd.wf .* rd.ntJ)) < tol

        (; node_ids_by_face  ) = rd.element_type
        @test sum(rd.wf[node_ids_by_face[1]]) ≈ 4
        # Note: this is not the true area of face 2. Because we map 
        # all faces back to the reference face, there is a factor of 
        # sqrt(2) difference from the true area. 
        @test sum(rd.wf[node_ids_by_face[2]]) ≈ 4 
        @test sum(rd.wf[node_ids_by_face[3]]) ≈ 4
        @test sum(rd.wf[node_ids_by_face[4]]) ≈ 2
        @test sum(rd.wf[node_ids_by_face[5]]) ≈ 2

        @test rd.Pq * rd.Vq ≈ I
        
        # 1/2 of a hex
        @test sum(rd.wq) ≈ 4
        
        @test StartUpDG.num_faces(Wedge()) == 5
        @test StartUpDG.num_vertices(Wedge()) == 6

        @test face_type(Wedge(), 3) == Quad()
        @test inverse_trace_constant(rd) ≈ 18.56357670538197
    end

    @testset "Misc Pyr" begin
        @test StartUpDG.num_faces(Pyr()) == 5
        @test StartUpDG.num_vertices(Pyr()) == 5
    end
end

inverse_trace_constant_compare(rd::RefElemData{3, <:Wedge, <:TensorProductWedge}) = ((9.0, 13.898979485566365, 19.292060161853993, 26.999999999999808), (12.0, 16.898979485566365, 22.292060161853993, 29.999999999999808), (16.0, 20.898979485566365, 26.292060161853993, 33.99999999999981), (21.0, 25.898979485566365, 31.292060161853993, 38.99999999999981))


@testset "Tensor reference elements" begin
    tol = 5e2*eps()
    @testset "Degree $tri_grad triangle" for tri_grad = [1, 2, 3, 4]
        @testset "Degree $line_grad line" for line_grad = [1, 2, 3, 4]
        line = RefElemData(Line(), line_grad)
        tri  = RefElemData(Tri(), tri_grad)
        tensor = TensorProductWedge(tri, line)
        rd = RefElemData(Wedge(), tensor, 1)

        @test rd.approximation_type.line.N == line_grad
        @test rd.approximation_type.tri.N == tri_grad

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
        
        @test abs(sum(rd.wf .* rd.nrJ)) < tol
        @test abs(sum(rd.wf .* rd.nsJ)) < tol
        @test abs(sum(rd.wf .* rd.ntJ)) < tol

        @unpack node_ids_by_face = rd.element_type
        @test sum(rd.wf[node_ids_by_face[1]]) ≈ 4
        # Note: this is not the true area of face 2. Because we map 
        # all faces back to the reference face, there is a factor of 
        # sqrt(2) difference from the true area. 
        @test sum(rd.wf[node_ids_by_face[2]]) ≈ 4 
        @test sum(rd.wf[node_ids_by_face[3]]) ≈ 4
        @test sum(rd.wf[node_ids_by_face[4]]) ≈ 2
        @test sum(rd.wf[node_ids_by_face[5]]) ≈ 2

        @test rd.Pq * rd.Vq ≈ I
        
        # 1/2 of a hex
        @test sum(rd.wq) ≈ 4
        
        @test StartUpDG.num_faces(Wedge()) == 5
        @test StartUpDG.num_vertices(Wedge()) == 6

        @test face_type(Wedge(), 1) == Quad()
        @test face_type(Wedge(), 2) == Quad()
        @test face_type(Wedge(), 3) == Quad()
        @test face_type(Wedge(), 4) == Tri()
        @test face_type(Wedge(), 5) == Tri()

        @test inverse_trace_constant(rd) ≈ inverse_trace_constant_compare(rd)[line_grad][tri_grad]
        end
    end
end