@testset "SummationByParts.jl extension tests" begin
    tol = 200 * eps()

    @testset "SBP{T} refactor" begin
        @test SBP() == SBP(StartUpDG.DefaultSBPType())
        @test SBP{Hicken}() == SBP(Hicken())
        @test SBP{Kubatko{LobattoFaceNodes}}().data isa Kubatko{LobattoFaceNodes}
        @test SBP{TensorProductLobatto}() isa SBP{TensorProductLobatto}
        approx_type = SBP(SummationByPartsDiagE{LobattoFaceNodes}())
        @test approx_type.data.quadrature_degree === nothing
        @test approx_type.data.tol == 100 * eps()
    end

    # checks the SBP property Q + Q' = E' * B * E for each coordinate direction
    function check_SBP_property(rd)
        (; wq, wf, Vf, Drst, nrstJ) = rd
        return all(map((D, nJ) -> begin
                           Q = Diagonal(wq) * D
                           norm(Q + Q' - Vf' * Diagonal(wf .* nJ) * Vf) < 100 * tol
                       end, Drst, nrstJ))
    end

    @testset "Tri, $F, q = $q, N = $N" for F in (LobattoFaceNodes, LegendreFaceNodes),
                                           N in (3, 5), q in (2 * N - 1, 2 * N)

        approx_type = SBP(SummationByPartsDiagE{F}(; quadrature_degree = q))
        rd = RefElemData(Tri(), approx_type, N)
        (; r, s, rq, sq, wq, Dr, Ds) = rd

        # resolved quadrature degree should be stored
        @test rd.approximation_type.data.quadrature_degree == q
        @test rd.approximation_type isa SBP{<:SummationByPartsDiagE{F}}

        # quadrature exactness
        f(k, r, s) = r^k + s^k + r^(k ÷ 2) * s^(k - k ÷ 2)
        rq2, sq2, wq2 = quad_nodes(Tri(), q)
        @test sum(wq) ≈ 2
        @test sum(wq .* f.(q, rq, sq)) ≈ sum(wq2 .* f.(q, rq2, sq2))

        # nodes = quadrature nodes, face nodes = subset of volume nodes
        @test rd.rst == rd.rstq
        @test rd.Nq == rd.Np
        @test rd.M == Diagonal(wq)
        @test all(vec.(rd.rstf) .≈ (x -> getindex(x, rd.Fmask)).(rd.rst))
        @test all(sum(rd.Vf, dims = 2) .== 1)

        # face node counts: N+2 Lobatto or N+1 Legendre nodes per edge
        @test rd.Nfq ÷ rd.num_faces == (F == LobattoFaceNodes ? N + 2 : N + 1)

        # differentiation accuracy 
        @test Dr * r .^ N ≈ N * r .^ (N - 1)
        @test Ds * s .^ N ≈ N * s .^ (N - 1)
        @test norm(Dr * s + Ds * r) < tol * rd.Np # roundoff scales with the number of nodes

        @test check_SBP_property(rd)
        @test inverse_trace_constant(rd) ≈ StartUpDG.eigenvalue_inverse_trace_constant(rd)

        # check that MeshData can be constructed
        md = MeshData(uniform_mesh(Tri(), 2)..., rd)
        @test md.x ≈ rd.V1 * md.VX[transpose(md.EToV)]
        @test all(md.J .> 0)
    end

    @testset "Tri default quadrature degree" begin
        N = 3
        rd = RefElemData(Tri(), SBP(SummationByPartsDiagE{LobattoFaceNodes}()), N)
        @test rd.approximation_type.data.quadrature_degree == 2 * N - 1
        rd = RefElemData(Tri(), SBP{SummationByPartsDiagE{LegendreFaceNodes}}(), N)
        @test rd.approximation_type.data.quadrature_degree == 2 * N - 1
    end

    @testset "Tet, N = $N, q = $q" for N in (1, 2, 3), q in (2 * N - 1, 2 * N)
        approx_type = SBP(SummationByPartsDiagE{LobattoFaceNodes}(; quadrature_degree = q))
        rd = RefElemData(Tet(), approx_type, N)
        (; r, s, t, rq, sq, tq, wq, Dr, Ds, Dt) = rd

        @test rd.approximation_type.data.quadrature_degree == q

        # quadrature exactness
        f(k, r, s, t) = r^k + s^k + t^k + r^(k ÷ 2) * t^(k - k ÷ 2)
        rq2, sq2, tq2, wq2 = quad_nodes(Tet(), q)
        @test sum(wq) ≈ 4 / 3
        @test sum(wq .* f.(q, rq, sq, tq)) ≈ sum(wq2 .* f.(q, rq2, sq2, tq2))

        # face quadrature weights are unscaled (scaling is carried by nrstJ): 4 faces of area 2
        @test sum(rd.wf) ≈ 4 * 2

        @test rd.rst == rd.rstq
        @test rd.M == Diagonal(wq)
        @test all(vec.(rd.rstf) .≈ (x -> getindex(x, rd.Fmask)).(rd.rst))
        @test all(sum(rd.Vf, dims = 2) .== 1)

        @test Dr * r .^ N ≈ N * r .^ (N - 1)
        @test Ds * s .^ N ≈ N * s .^ (N - 1)
        @test Dt * t .^ N ≈ N * t .^ (N - 1)
        @test norm(Dr * s + Ds * r) < tol * rd.Np # roundoff scales with the number of nodes
        @test norm(Dr * t + Dt * r) < tol * rd.Np

        @test check_SBP_property(rd)
        @test inverse_trace_constant(rd) ≈ StartUpDG.eigenvalue_inverse_trace_constant(rd)

        md = MeshData(uniform_mesh(Tet(), 2)..., rd)
        @test md.x ≈ rd.V1 * md.VX[transpose(md.EToV)]
        @test all(md.J .> 0)
    end

    @testset "Tet default SBP type" begin
        rd = RefElemData(Tet(), SBP(), 2)
        @test rd.approximation_type isa SBP{<:SummationByPartsDiagE{LobattoFaceNodes}}
        @test rd.approximation_type.data.quadrature_degree == 3
    end

    @testset "Errors" begin
        N = 3
        approx_type = SBP(SummationByPartsDiagE{LobattoFaceNodes}(quadrature_degree = 2 *
                                                                                      N - 2))
        @test_throws ArgumentError RefElemData(Tri(), approx_type, N)
        approx_type = SBP(SummationByPartsDiagE{LegendreFaceNodes}())
        @test_throws ArgumentError RefElemData(Tet(), approx_type, N)
        # q = 21 is not available in SummationByParts.jl
        @test_throws Exception RefElemData(Tri(),
                                           SBP(SummationByPartsDiagE{LobattoFaceNodes}()),
                                           11)
    end
end
