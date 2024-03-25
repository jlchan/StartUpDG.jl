"""
    function RefElemData(elementType::Line, approxType::SBP, N)
    function RefElemData(elementType::Quad, approxType::SBP, N)
    function RefElemData(elementType::Hex,  approxType::SBP, N)
    function RefElemData(elementType::Tri,  approxType::SBP, N)
    
SBP reference element data for `Quad()`, `Hex()`, and `Tri()` elements. 

For `Line()`, `Quad()`, and `Hex()`, `approxType` is `SBP{TensorProductLobatto}`.

For `Tri()`, `approxType` can be `SBP{Kubatko{LobattoFaceNodes}}`, `SBP{Kubatko{LegendreFaceNodes}}`, or `SBP{Hicken}`. 
"""
function RefElemData(elementType::Line, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps(), kwargs...)

    rd = RefElemData(elementType, N; quad_rule_vol = gauss_lobatto_quad(0,0,N), kwargs...)        
    
    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Quad, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps(), kwargs...)

    # make 2D SBP nodes/weights
    r1D, w1D = gauss_lobatto_quad(0, 0, N)
    sq, rq = vec.(NodesAndModes.meshgrid(r1D)) # this is to match ordering of nrstJ
    wr, ws = vec.(NodesAndModes.meshgrid(w1D)) 
    wq = wr .* ws
    quad_rule_vol = (rq, sq, wq)
    quad_rule_face = (r1D, w1D)

    rd = RefElemData(elementType, N; quad_rule_vol = quad_rule_vol, quad_rule_face = quad_rule_face, kwargs...)

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Hex, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps(), kwargs...)

    # make 2D SBP nodes/weights
    r1D, w1D = gauss_lobatto_quad(0, 0, N)
    rf, sf = vec.(NodesAndModes.meshgrid(r1D, r1D))
    wr, ws = vec.(NodesAndModes.meshgrid(w1D, w1D))
    wf = wr .* ws
    sq, rq, tq = vec.(NodesAndModes.meshgrid(r1D, r1D, r1D)) # this is to match ordering of nrstJ
    wr, ws, wt = vec.(NodesAndModes.meshgrid(w1D, w1D, w1D)) 
    wq = wr .* ws .* wt
    quad_rule_vol = (rq, sq, tq, wq)
    quad_rule_face = (rf, sf, wf)

    rd = RefElemData(elementType, N; quad_rule_vol = quad_rule_vol, quad_rule_face = quad_rule_face, kwargs...)

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Tri, approxType::SBP, N; tol = 100*eps(), kwargs...)
    
    quad_rule_vol, quad_rule_face = diagE_sbp_nodes(elementType, approxType, N)

    # build polynomial reference element using quad rules; will be modified to create SBP RefElemData
    rd = RefElemData(elementType, Polynomial(), N; quad_rule_vol=quad_rule_vol, quad_rule_face=quad_rule_face, kwargs...)

    # determine Fmask = indices of face nodes among volume nodes
    Ef, Fmask = build_Ef_Fmask(rd)

    # Build traditional SBP operators from hybridized operators. See Section 3.2 of 
    # [High-order entropy stable dG methods for the SWE](https://arxiv.org/pdf/2005.02516.pdf)
    # by Wu and Chan 2021. [DOI](https://doi.org/10.1016/j.camwa.2020.11.006)
    (Qrh, Qsh), _ = hybridized_SBP_operators(rd)
    Nq = length(rd.wq)
    Vh_sbp = [I(Nq); Ef]
    Qr = Vh_sbp' * Qrh * Vh_sbp
    Qs = Vh_sbp' * Qsh * Vh_sbp
    Dr,Ds = (x -> diagm(1 ./ rd.wq) * x).((Qr, Qs))
    
    rd = @set rd.rst = quad_rule_vol[1:2]   # set nodes = SBP nodes
    rd = @set rd.rstq = quad_rule_vol[1:2]  # set quad nodes = SBP nodes
    rd = @set rd.Drst = (Dr, Ds)
    rd = @set rd.Fmask = vec(Fmask)

    # TODO: make these more efficient with custom operators?
    rd = @set rd.Vf = droptol!(sparse(Ef), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) 

    # make V1 the interpolation matrix from triangle vertices to SBP nodal points
    rd = @set rd.V1 = vandermonde(elementType, N, rd.rst...) / rd.VDM * rd.V1

    # Vp operator = projects SBP nodal vector onto degree N polynomial, then interpolates to plotting points
    rd = @set rd.Vp = vandermonde(elementType, N, rd.rstp...) / rd.VDM * rd.Pq

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end


#####
##### Utilities for SBP 
#####

# - HDF5 file created using MAT.jl and the following code:
# vars = matread("src/data/sbp_nodes/KubatkoQuadratureRules.mat")
# h5open("src/data/sbp_nodes/KubatkoQuadratureRules.h5", "w") do file
#     for qtype in ("Q_GaussLobatto", "Q_GaussLegendre")
#         group = create_group(file, qtype) # create a group
#         for fieldname in ("Points", "Domain", "Weights")
#             subgroup = create_group(group, fieldname)
#             for N in 1:length(vars[qtype])
#                 subgroup[string(N)] = vars[qtype][N][fieldname]
#             end
#         end
#     end
# end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LobattoFaceNodes}}, N)    
    
    if N==6
        @warn "N=6 SBP operators with quadrature strength 2N-1 and Lobatto face nodes may require very small timesteps."
    end
    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LobattoFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLobatto"]["Points"][string(N)][]
    r, s = (rs[:, i] for i = 1:size(rs, 2))
    w = vec(vars["Q_GaussLobatto"]["Weights"][string(N)][])
    quad_rule_face = gauss_lobatto_quad(0, 0, N+1)     

    return (r, s, w), quad_rule_face 
end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LegendreFaceNodes}}, N)    

    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LegendreFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLegendre"]["Points"][string(N)][]
    r, s = (rs[:, i] for i = 1:size(rs, 2))
    w = vec(vars["Q_GaussLegendre"]["Weights"][string(N)][])
    quad_rule_face = gauss_quad(0, 0, N)

    return (r, s, w), quad_rule_face 
end

parsevec(type, str) = str |>
  (x -> split(x, ", ")) |>
  (x -> map(y -> parse(type, y), x))
  
function diagE_sbp_nodes(elem::Tri, approxType::SBP{Hicken}, N)    
    
    if N > 4
        @error "N > 4 triangular `SBP{Hicken}` operators are not available."
    end

    # from Jason Hicken https://github.com/OptimalDesignLab/SummationByParts.jl/tree/work
    lines = readlines((@__DIR__)*"/data/sbp_nodes/tri_diage_p$N.dat") 
    r = parsevec(Float64,lines[11])
    s = parsevec(Float64,lines[12])
    w = parsevec(Float64,lines[13])

    # convert Hicken format to biunit right triangle
    r = @. 2*r-1 
    s = @. 2*s-1
    w = 2.0 * w/sum(w)

    quad_rule_face = gauss_lobatto_quad(0,0,N+1) 

    return (r,s,w), quad_rule_face 
end

function build_Ef_Fmask(rd_sbp::RefElemData; tol = 100*eps())
    
    (; rq, sq, rf, sf, Nfaces ) = rd_sbp   
    rf, sf = (x->reshape(x, length(rf) ÷ Nfaces, Nfaces)).((rf, sf))
    Fmask = zeros(Int, length(rf) ÷ Nfaces, Nfaces) # 
    Ef = zeros(length(rf), length(rq)) # extraction matrix
    for i in eachindex(rq)
        for f = 1:rd_sbp.Nfaces
            id = findall(@. abs(rq[i]-rf[:,f]) + abs(sq[i]-sf[:,f]) .< tol)
            Fmask[id,f] .= i
            Ef[id .+ (f-1) * size(rf,1), i] .= 1
        end
    end
    return Ef, Fmask
end

get_face_nodes(x::AbstractVector, Fmask) = view(x, Fmask)
get_face_nodes(x::AbstractMatrix, Fmask) = view(x, Fmask, :)

function _convert_RefElemData_fields_to_SBP(rd, approx_type::SBP)
    rd = @set rd.M = Diagonal(rd.wq)
    rd = @set rd.Pq = I
    rd = @set rd.Vq = I
    rd = @set rd.approximation_type = approx_type
    return rd
end

"""
    function hybridized_SBP_operators(rd::RefElemData{DIMS}) 

Constructs hybridized SBP operators given a `RefElemData`. Returns operators `Qrsth..., VhP, Ph`.
"""
function hybridized_SBP_operators(rd)
    (; M, Vq, Pq, Vf, wf, Drst, nrstJ ) = rd
    Qrst = (D->Pq' * M * D * Pq).(Drst)
    Ef = Vf * Pq
    Brst = (nJ->diagm(wf .* nJ)).(nrstJ)
    Qrsth = ((Q, B)->.5*[Q-Q' Ef'*B; -B*Ef B]).(Qrst, Brst)
    Vh = [Vq; Vf]
    Ph = M \ transpose(Vh)
    VhP = Vh * Pq
    return Qrsth, VhP, Ph, Vh
end

# default to doing nothing
map_nodes_to_symmetric_element(element_type, rst...) = rst

# for triangles, map to an equilateral triangle
function map_nodes_to_symmetric_element(::Tri, r, s)
    # biunit right triangular vertices
    v1, v2, v3 = SVector{2}.(zip(nodes(Tri(), 1)...))

    denom = (v2[2] - v3[2]) * (v1[1] - v3[1]) + (v3[1] - v2[1]) * (v1[2] - v3[2])
    L1 = @. ((v2[2] - v3[2]) * (r - v3[1]) + (v3[1] - v2[1]) * (s - v3[2])) / denom
    L2 = @. ((v3[2] - v1[2]) * (r - v3[1]) + (v1[1] - v3[1]) * (s - v3[2])) / denom
    L3 = @. 1 - L1 - L2

    # equilateral vertices
    v1 = SVector{2}(2 * [-.5, -sqrt(3) / 6])
    v2 = SVector{2}(2 * [.5, -sqrt(3)/6])
    v3 = SVector{2}(2 * [0, sqrt(3)/3])

    x = @. v1[1] * L1 + v2[1] * L2 + v3[1] * L3
    y = @. v1[2] * L1 + v2[2] * L2 + v3[2] * L3

    return x, y
end


"""
    function sparse_low_order_SBP_operators(rd; factor=1.01)

Constructs sparse low order SBP operators given a `RefElemData`. 
Returns operators `Qrst..., E ≈ Vf * Pq` that satisfy a generalized 
summation-by-parts (GSBP) property:

        `Q_i + Q_i^T = E' * B_i * E`
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
    @assert sorted_eigvals[2] > 100 * eps() 
    
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

    d = (ones(size(D, 1))' * D)'
    ids = findall(@. abs(d) > 1e2 * eps())

    Z = nullspace(D)[ids, :]
    A = pinv(D)[ids,:]
    X = randn(size(Z, 2), size(D, 1))
    y = randn(length(ids))
    e = ones(size(D, 1))

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

# specialize for single dimension
function subcell_limiting_operators(rd::RefElemData{1})
    Qrst, _ = sparse_low_order_SBP_operators(rd)
    Δrst, Rrst = subcell_limiting_operators(Qrst[1])
    return Δrst, Rrst
end    
