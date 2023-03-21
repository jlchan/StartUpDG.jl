"""
    function RefElemData(elementType::Line, approxType::SBP, N)
    function RefElemData(elementType::Quad, approxType::SBP, N)
    function RefElemData(elementType::Hex,  approxType::SBP, N)
    function RefElemData(elementType::Tri,  approxType::SBP, N)
    
SBP reference element data for `Quad()`, `Hex()`, and `Tri()` elements. 

For `Line()`, `Quad()`, and `Hex()`, `approxType` is `SBP{TensorProductLobatto}`.

For `Tri()`, `approxType` can be `SBP{Kubatko{LobattoFaceNodes}}`, `SBP{Kubatko{LegendreFaceNodes}}`, or `SBP{Hicken}`. 
"""
function RefElemData(elementType::Line, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps())

    rd = RefElemData(elementType, N, quad_rule_vol = gauss_lobatto_quad(0,0,N))        
    
    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Quad, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps())

    # make 2D SBP nodes/weights
    r1D, w1D = gauss_lobatto_quad(0, 0, N)
    sq, rq = vec.(NodesAndModes.meshgrid(r1D)) # this is to match ordering of nrstJ
    wr, ws = vec.(NodesAndModes.meshgrid(w1D)) 
    wq = wr .* ws
    quad_rule_vol = (rq, sq, wq)
    quad_rule_face = (r1D, w1D)

    rd = RefElemData(elementType, N, quad_rule_vol = quad_rule_vol, quad_rule_face = quad_rule_face)

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Hex, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps())

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

    rd = RefElemData(elementType, N, quad_rule_vol = quad_rule_vol, quad_rule_face = quad_rule_face)

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Tri, approxType::SBP, N; tol = 100*eps())
    
    quad_rule_vol, quad_rule_face = diagE_sbp_nodes(elementType, approxType, N)

    # build polynomial reference element using quad rules; will be modified to create SBP RefElemData
    rd = RefElemData(elementType, Polynomial(), N; quad_rule_vol=quad_rule_vol, quad_rule_face=quad_rule_face)

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

# for triangles and tets, map to an equilateral triangle/tet
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
    function sparse_low_order_SBP_operators(rd)

Constructs sparse low order SBP operators given a `RefElemData`. 
Returns operators `Qrst..., E ≈ Vf*Pq` that satisfy the GSBP property
        `Q_i + Q_i^T = E' * B_i * E`
"""
function sparse_low_order_SBP_operators(rd::RefElemData{NDIMS}) where {NDIMS}
    (; Pq, Vf, rstq, wf, nrstJ) = rd

    # if element is a simplex, convert nodes to an equilateral triangle for symmetry
    rstq = map_nodes_to_symmetric_element(rd.element_type, rstq...)

    # build distance and adjacency matrices
    D = [norm(SVector{NDIMS}(getindex.(rstq, i)) - SVector{NDIMS}(getindex.(rstq, j))) 
            for i in eachindex(first(rstq)), j in eachindex(first(rstq))]
    A = zeros(Int, length(first(rstq)), length(first(rstq)))
    for i in axes(D, 1)
        dist = sort(view(D, i, :))[NDIMS + 2] # heuristic cutoff - smallest distance = 0, and we take NDIMS more
        for j in findall(view(D, i, :) .< 1.01 * dist)
            A[i, j] = 1
        end
    end
    A = (A + A' .> 0) # symmetrize
    L = Diagonal(vec(sum(A, dims=2))) - A
    sorted_eigvals = sort(eigvals(L))
    @assert sorted_eigvals[2] > 100 * eps() # check that there's only one zero null vector

    # TODO: sparsify the E matrix
    E = Vf * Pq    
    Brst = (nJ -> diagm(wf .* nJ)).(nrstJ)

    # compute potential 
    e = ones(size(L, 2))
    right_hand_sides = map(B -> 0.5 * sum(E' * B, dims=2), Brst)
    psi_augmented = map(b -> [L e; e' 0] \ [b; 0], right_hand_sides)
    psi = map(x -> x[1:end-1], psi_augmented)

    # construct sparse skew part
    function construct_skew_matrix(ψ)
        S = zeros(length(ψ), length(ψ))
        for i in axes(S, 1), j in axes(S, 2)
            if A[i,j] > 0
                S[i,j] = ψ[j] - ψ[i]
            end
        end
        return S
    end
    Srst = construct_skew_matrix.(psi)
    Qrst = map((S,B) -> S + 0.5 * E' * B * E, Srst, Brst)
    return Qrst, E
end

