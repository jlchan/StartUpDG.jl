"""
    function sparse_low_order_SBP_operators(rd; factor=1.01)

Constructs sparse low order SBP operators given a `RefElemData`. 
Returns operators `Qrst..., E ≈ Vf * Pq` that satisfy a generalized 
summation-by-parts (GSBP) property:

        `Q_i + Q_i^T = E' * B_i * E`

`factor` is a scaling which determines how close a node must be to 
another node to be considered a neighbor.
"""
function sparse_low_order_SBP_operators(rd::RefElemData{NDIMS}; factor=1.01) where {NDIMS}
    (; Pq, Vf, rstq, wf, nrstJ) = rd

    # if element is a simplex, convert nodes to an equilateral triangle for symmetry
    rstq = map_nodes_to_symmetric_element(rd.element_type, rstq...)

    # build distance and adjacency matrices
    D = [norm(SVector{NDIMS}(getindex.(rstq, i)) - SVector{NDIMS}(getindex.(rstq, j))) 
            for i in eachindex(first(rstq)), j in eachindex(first(rstq))]
    A = zeros(Int, length(first(rstq)), length(first(rstq)))
    for i in axes(D, 1)
        # heuristic cutoff - we take NDIMS + 1 neighbors, but the smallest distance = 0, 
        # so we need to access the first NDIMS + 2 sorted distance entries.
        dist = sort(view(D, i, :))[NDIMS + 2] 
        for j in findall(view(D, i, :) .< factor * dist)
            A[i, j] = 1
        end
    end
    A = (A + A' .> 0) # symmetrize
    L = Diagonal(vec(sum(A, dims=2))) - A
    sorted_eigvals = sort(eigvals(L))

    # check that there's only one zero null vector
    if sorted_eigvals[2] < 10 * size(A, 1) * eps()
        error("Warning: the connectivity between nodes yields a graph with " * 
            "more than one connected component. Try increasing the value of `factor`.")
    end
    
    # note: this part doesn't do anything for a nodal SBP operator. In that case, E_dense = E = Vf 
    # and E should just reduce to Vf, which is a face extraction operator. 
    E_dense = Vf * Pq
    E = zeros(size(E_dense))
    for i in axes(E, 1)
        # find all j such that E[i,j] ≥ 0.5, e.g., points which positively contribute to at least half of the 
        # interpolation. These seem to be associated with volume points "j" that are close to face point "i".
        ids = findall(E_dense[i, :] .>= 0.5)
        E[i, ids] .= E_dense[i, ids] ./ sum(E_dense[i, ids]) # normalize so sum(E, dims=2) = [1, 1, ...] still.
    end
    Brst = (nJ -> diagm(wf .* nJ)).(nrstJ)

    # compute potential 
    e = ones(size(L, 2))
    right_hand_sides = map(B -> 0.5 * sum(E' * B, dims=2), Brst)
    psi_augmented = map(b -> [L e; e' 0] \ [b; 0], right_hand_sides)
    psi = map(x -> x[1:end-1], psi_augmented)

    # construct sparse skew part
    function construct_skew_matrix_from_potential(ψ)
        S = zeros(length(ψ), length(ψ))
        for i in axes(S, 1), j in axes(S, 2)
            if A[i,j] > 0
                S[i,j] = ψ[j] - ψ[i]
            end
        end
        return S
    end

    Srst = construct_skew_matrix_from_potential.(psi)
    Qrst = map((S, B) -> S + 0.5 * E' * B * E, Srst, Brst)
    return sparse.(Qrst), sparse(E)
end

function sparse_low_order_SBP_1D_operators(rd::RefElemData) 
    E = zeros(2, rd.N+1)
    E[1, 1] = 1
    E[2, end] = 1

    # create volume operators
    Q =  diagm(1 => ones(rd.N), -1 => -ones(rd.N))
    Q[1,1] = -1
    Q[end,end] = 1
    Q = 0.5 * Q

    return (sparse(Q),), sparse(E)
end

sparse_low_order_SBP_operators(rd::RefElemData{1, Line, <:Union{<:SBP, <:Polynomial{Gauss}}}) = 
    sparse_low_order_SBP_1D_operators(rd)

function diagonal_1D_mass_matrix(N, ::SBP)
    _, w1D = gauss_lobatto_quad(0, 0, N)
    return Diagonal(w1D)
end

function diagonal_1D_mass_matrix(N, ::Polynomial{Gauss})
    _, w1D = gauss_quad(0, 0, N)
    return Diagonal(w1D)
end

function sparse_low_order_SBP_operators(rd::RefElemData{2, Quad, <:Union{<:SBP, <:Polynomial{Gauss}}}) 
    (Q1D,), E1D = sparse_low_order_SBP_1D_operators(rd)

    # permute face node ordering for the first 2 faces
    ids = reshape(1:(rd.N+1) * 2, :, 2)
    Er = zeros((rd.N+1) * 2, rd.Np)
    Er[vec(ids'), :] .= kron(I(rd.N+1), E1D)
    Es = kron(E1D, I(rd.N+1))
    E = vcat(Er, Es)

    M1D = diagonal_1D_mass_matrix(rd.N, rd.approximation_type)
    Qr = kron(M1D, Q1D)
    Qs = kron(Q1D, M1D)

    return sparse.((Qr, Qs)), sparse(E)
end

function sparse_low_order_SBP_operators(rd::RefElemData{3, Hex, <:Union{<:SBP, <:Polynomial{Gauss}}}) 
    (Q1D,), E1D = sparse_low_order_SBP_1D_operators(rd)

    # permute face nodes for the faces in the ±r and ±s directions
    ids = reshape(1:(rd.N+1)^2 * 2, rd.N+1, :, 2) 
    Er, Es, Et = ntuple(_ -> zeros((rd.N+1)^2 * 2, rd.Np), 3)
    Er[vec(permutedims(ids, [2, 3, 1])), :] .= kron(I(rd.N+1), E1D, I(rd.N+1))
    Es[vec(permutedims(ids, [3, 2, 1])), :] .= kron(I(rd.N+1), I(rd.N+1), E1D)
    Et .= kron(E1D, I(rd.N+1), I(rd.N+1))

    # create boundary extraction matrix
    E = vcat(Er, Es, Et)

    M1D = diagonal_1D_mass_matrix(rd.N, rd.approximation_type)
    Qr = kron(M1D, Q1D, M1D)
    Qs = kron(M1D, M1D, Q1D)
    Qt = kron(Q1D, M1D, M1D)

    return sparse.((Qr, Qs, Qt)), sparse(E)
end

# builds subcell operators D * R such that for r such that sum(r) = 0, 
# D * Diagonal(θ) * R * r is also diagonal for any choice of θ. This is
# useful in subcell limiting (see, for example  
# https://doi.org/10.1016/j.compfluid.2022.105627 for a 1D example)
function subcell_limiting_operators(Qr::AbstractMatrix; tol = 100 * eps())
    Qr = sparse(Qr)
    Sr = droptol!(Qr - Qr', tol)
    Br = droptol!(Qr + Qr', tol)

    num_interior_fluxes = nnz(Sr) ÷ 2
    num_boundary_indices = nnz(Br)
    num_fluxes = num_interior_fluxes + num_boundary_indices

    interior_indices = findall(Sr .> tol)
    boundary_indices = findall(abs.(Br) .> tol)

    matrix_indices = zeros(Int, size(Sr))
    for (i, ij) in enumerate(boundary_indices)
        matrix_indices[ij] = i
        matrix_indices[ij.I[2], ij.I[1]] = i
    end

    for (i, ij) in enumerate(interior_indices)
        matrix_indices[ij] = i + length(boundary_indices)
        matrix_indices[ij.I[2], ij.I[1]] = i + length(boundary_indices)
    end

    D = zeros(size(Qr, 2), num_fluxes)
    for i in axes(matrix_indices, 1), j in axes(matrix_indices, 2)
        if matrix_indices[i, j] !== 0
            D[i, matrix_indices[i, j]] = Qr[i, j]    
        end
    end

    d = vec(ones(size(D, 1))' * D)
    ids = findall(@. abs(d) > 1e2 * eps())

    Z = nullspace(D)[ids, :]
    A = pinv(D)[ids,:]

    # compute R via linear algebra
    m, n = size(A)
    z = size(Z, 2)
    C = hcat(kron(I(n), Z), kron(ones(n), I(m)))
    sol = pinv(C) * -vec(A)
    X = reshape(sol[1:n*z], (z, n))
    R = pinv(D) + nullspace(D) * X

    # check if R satisfies our pseudoinverse condition
    @assert norm(D * R - I(size(D,1))) < tol * length(D)

    return D, R
end

"""
    Δrst, Rrst = subcell_limiting_operators(rd::RefElemData)

Returns tuples of subcell limiting operators Drst = (Δr, Δs, ...) and R = (Rr, Rs, ...) 
such that for r where sum(r) = 0, sum(D * Diagonal(θ) * R * r) = 0 for any choice of θ. 
These operators are useful for conservative subcell limiting techniques (see 
https://doi.org/10.1016/j.compfluid.2022.105627 for an example of such an approach on 
tensor product elements). 

Sparse SBP operators used in an intermediate step when buidling these subcell 
limiting operators; by default, these operators are constructed using
`sparse_low_order_SBP_operators`. To construct subcell limiting operators for a 
general SBP operator, one can use the following:

    Δ, R = subcell_limiting_operators(Q::AbstractMatrix; tol = 100 * eps())
"""
function subcell_limiting_operators(rd::RefElemData)
    Qrst, _ = sparse_low_order_SBP_operators(rd)
    Δrst, Rrst = subcell_limiting_operators.(Qrst)
    return zip(Δrst, Rrst)
end    

# specialize for single dimension to return tuples of operators
function subcell_limiting_operators(rd::RefElemData{1})
    Qrst, _ = sparse_low_order_SBP_operators(rd)
    Δrst, Rrst = first(subcell_limiting_operators.(Qrst))
    return (Δrst,), (Rrst,)
end

# TODO: specialize for quadrilaterals and hex? 

