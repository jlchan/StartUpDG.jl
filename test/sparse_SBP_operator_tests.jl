@testset "Sparse SBP operators for approx_type <: Polynomial" begin
    tol = 5e2*eps()
    N = 3    
    @testset "Line" begin 
        rd = RefElemData(Line(), Polynomial{Gauss}(), N)
        (Q,), E = sparse_low_order_SBP_operators(rd)
        @test norm(sum(Q, dims=2)) < tol
        @test Q + Q' ≈ E' * Diagonal([-1,1]) * E
    end

    @testset "Tri" begin
        rd = RefElemData(Tri(), N)
        (Qr, Qs), E = sparse_low_order_SBP_operators(rd)
        @test norm(sum(Qr, dims=2)) < tol
        @test norm(sum(Qs, dims=2)) < tol
        @test Qr + Qr' ≈ E' * Diagonal(rd.wf .* rd.nrJ) * E
        @test Qs + Qs' ≈ E' * Diagonal(rd.wf .* rd.nsJ) * E
    end

    @testset "Quad (SBP)" begin
        rd = RefElemData(Quad(), SBP(), N)
        (Qr, Qs), E = sparse_low_order_SBP_operators(rd)
        @test norm(sum(Qr, dims=2)) < tol
        @test norm(sum(Qs, dims=2)) < tol
        @test Qr + Qr' ≈ E' * Diagonal(rd.wf .* rd.nrJ) * E
        @test Qs + Qs' ≈ E' * Diagonal(rd.wf .* rd.nsJ) * E
    end

    @testset "Quad (Polynomial{Gauss})" begin
        rd = RefElemData(Quad(), Gauss(), N)
        (Qr, Qs), E = sparse_low_order_SBP_operators(rd)
        @test norm(sum(Qr, dims=2)) < tol
        @test norm(sum(Qs, dims=2)) < tol
        @test Qr + Qr' ≈ E' * Diagonal(rd.wf .* rd.nrJ) * E
        @test Qs + Qs' ≈ E' * Diagonal(rd.wf .* rd.nsJ) * E
    end

    @testset "Hex (Polynomial{Gauss})" begin
        rd = RefElemData(Hex(), Gauss(), N)
        (Qr, Qs, Qt), E = sparse_low_order_SBP_operators(rd)
        @test norm(sum(Qr, dims=2)) < tol
        @test norm(sum(Qs, dims=2)) < tol
        @test norm(sum(Qt, dims=2)) < tol
        @test Qr + Qr' ≈ E' * Diagonal(rd.wf .* rd.nrJ) * E
        @test Qs + Qs' ≈ E' * Diagonal(rd.wf .* rd.nsJ) * E
        @test Qt + Qt' ≈ E' * Diagonal(rd.wf .* rd.ntJ) * E
    end
end

@testset "Subcell limiting operators for approx_type <: SBP" begin    
    @testset "$elem_type" for elem_type in [Tri()] 
        N = 3
        rd = RefElemData(elem_type, SBP(), N)
        Qrst, E = sparse_low_order_SBP_operators(rd)
        Δrst, Rrst = subcell_limiting_operators(rd)
        for dim in eachindex(Rrst)
            @test Δrst[dim] * Rrst[dim] ≈ I

            # check conservation for a random set of limiting factors
            r = randn(size(Rrst[dim], 2))
            r = r .- sum(r) / length(r) # so that sum(r) = 0
            θ = rand(size(Rrst[dim], 1))
            @test abs(sum(Δrst[dim] * Diagonal(θ) * Rrst[dim] * r)) < 10 * length(r) * eps()
        end
    end

    @testset "Tensor product elements: $elem_type" for elem_type in [Quad(), Hex()] 
        N = 2
        rd = RefElemData(elem_type, SBP(), N)
        Qrst, E = sparse_low_order_SBP_operators(rd)
        Δrst, Rrst = subcell_limiting_operators(rd)

        for dim in eachindex(Rrst)
            conservation_vectors = nullspace(Matrix(Qrst[dim]))
            Δ = Δrst[dim]
            R = Rrst[dim]

            # make random conservative vector
            r = randn(length(rd.r))
            r = r - conservation_vectors * (conservation_vectors' * r)

            # create limiting factors
            θ = rand(size(R, 1))

            @test norm(conservation_vectors' * (Δ * (θ .* (R * r)))) < 100 * eps()
        end
    end
end