# dimension of a cut cell polynomial space
@inline Np_cut(N) = (N + 1) * (N + 2) ÷ 2 

"""
    `PhysicalFrame{NDIMS} <: AbstractElemShape{NDIMS}`
    
`PhysicalFrame` element type. Uses a total degree N approximation space, but is 
computed with a tensor product Legendre basis as opposed to a triangular PKDO basis.
Stores fields `shifting` and `scaling` to shift/scale physical coordinates so that 
they are on the reference element. 

    PhysicalFrame()
    PhysicalFrame(x, y)
    PhysicalFrame(x, y, vx, vy): stores coordinates `vx, vy` of background Cartesian cell 

Constructors for a PhysicalFrame object (optionally uses arrays of points `x`, `y` on a cut element).
"""
struct PhysicalFrame{NDIMS, Shifting <: SVector{NDIMS}, Scaling <: SVector{NDIMS}, VXYZ} <: AbstractElemShape{NDIMS}
    shifting::Shifting
    scaling::Scaling
    vxyz::VXYZ # coordinates of background Cartesian cell
end

# defaults to 2D for now
PhysicalFrame() = PhysicalFrame(Val{2}())

# This is an alternative method of computing shifting/scaling factors, 
# but results in nodes on [-1,1]^2 being mapped outside the background element. 
# This is currently unused. 
function get_shifting_and_scaling_centroid(x,y)
    shifting = SVector(mean(x), mean(y))
    scaling = SVector( 1 / maximum(abs.(x .- shifting[1])), 1 / maximum(abs.(y .- shifting[2])) )

    return shifting, scaling
end

function get_shifting_and_scaling_maxfill(x,y)
    scaling = SVector(map(x -> 2 / (x[2] - x[1]), (extrema(x), extrema(y))))
    shifting = SVector( mean(extrema(x)), mean(extrema(y)) )

    return shifting, scaling
end

# default shifting and scaling
function PhysicalFrame(ndims::Val{2}) 
    shifting = SVector(0.0, 0.0)
    scaling = SVector(1.0, 1.0)
    vxyz = (SVector(-1., 1.), SVector(-1., 1.)) # set background coordinates to reference cell
    return PhysicalFrame(shifting, scaling, vxyz) 
end

function PhysicalFrame(x, y)
    shifting, scaling = get_shifting_and_scaling_maxfill(x,y)
    return PhysicalFrame(shifting, scaling, nothing)
end

function PhysicalFrame(x, y, vx, vy)
    shifting, scaling = get_shifting_and_scaling_maxfill(x,y)
    vxyz = (vx, vy)
    return PhysicalFrame(shifting, scaling, vxyz)
end

function shift_and_scale(elem::PhysicalFrame{2}, x, y)
    (; shifting, scaling ) = elem
    r = @. (x - shifting[1]) * scaling[1]
    s = @. (y - shifting[2]) * scaling[2]
    return r, s
end

function map_nodes_to_background_cell(elem::PhysicalFrame{2}, r, s)
    (; vxyz ) = elem
    vx, vy = vxyz
    dx, dy = diff(vx), diff(vy)
    x = @. 0.5 * (1 + r) * dx + vx[1]
    y = @. 0.5 * (1 + s) * dy + vy[1]
    return x, y
end

function map_nodes_to_cutcell_boundingbox(elem::PhysicalFrame{2}, r, s)
    (; shifting, scaling ) = elem 
    x = @. r / scaling[1] + shifting[1]
    y = @. s / scaling[2] + shifting[2]
    return x, y
end


function NodesAndModes.basis(elem::PhysicalFrame{2}, N, x, y)
    Np = (N + 1) * (N + 2) ÷ 2

    r, s = shift_and_scale(elem, x, y)

    (; scaling ) = elem

    sk = 1
    V, Vr, Vs = ntuple(x->zeros(length(r), Np), 3)
    for j = 0:N
        P_j = jacobiP(s, 0, 0, j)
        dP_j = grad_jacobiP(s, 0, 0, j)
        for i = 0:N-j
            P_i = jacobiP(r, 0, 0, i)
            dP_i = grad_jacobiP(r, 0, 0, i)
            @. V[:, sk]  = P_i * P_j
            @. Vr[:, sk] = dP_i * P_j * scaling[1]
            @. Vs[:, sk] = P_i * dP_j * scaling[2]
            sk += 1
        end
    end    
    return V, Vr, Vs
end

# integrated Legendre polynomials. basis = degree N 1D basis
# see also https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620280813
function antidiff_basis_1D(basis_1D, N, r)
    rq, wq = gauss_quad(0, 0, N+1) 
    Vout = zeros(eltype(r), length(r), N+1)
    for i in eachindex(r)
        h = r[i] + 1
        Vout[i, :] .= h/2 .* vec(wq' * basis_1D(N, map_to_interval.(rq, -1, r[i])))
    end
    return Vout
end

function antidiff_basis_2D(N, r, s)
    Ix, Iy = ntuple(_ -> zeros(eltype(r), length(r), Np_cut(N)), 2)
    Vr, _ = basis(Line(), N, r)
    Vs, _ = basis(Line(), N, s)
    basis_1D = (N, r) -> basis(Line(), N, r)[1]
    Ir = antidiff_basis_1D(basis_1D, N, r)
    Is = antidiff_basis_1D(basis_1D, N, s)

    sk = 1
    for j in 0:N, i in 0:N-j
        Ix[:, sk] .= Ir[:, i+1] .* Vs[:, j+1]
        Iy[:, sk] .= Vr[:, i+1] .* Is[:, j+1]
        sk += 1
    end
    return Ix, Iy
end

# computes antidifferentiation operators for total degree N bases 
# on the reference quadrilateral element. These operators map from 
# the space of total degree N to total degree N+1 polynomials. 
function antidiff_operators(N)
    # antidiff basis: P^N -> nodes
    # inv(VNp1): nodes -> P^{N+1}
    #   => inv(VNp1) * antidiff basis: P^N -> P^{N+1}
    r, s = nodes(Quad(), N+1)
    VNp1, VNp1r, VNp1s = basis(PhysicalFrame(), N+1, r, s)
    Ir, Is = map(A -> VNp1 \ A, antidiff_basis_2D(N, r, s))

    # check that the antidiff operator is the pseudo-inverse of the P^{N+1} -> P^N differentiation matrix
    V, _ = basis(PhysicalFrame(), N, r, s)
    Dr, Ds = V \ VNp1r, V \ VNp1s 
    @assert norm(Dr * Ir - I) + norm(Ds * Is - I) < 100 * eps() * size(Ir, 1)

    return Ir, Is
end

import NodesAndModes: equi_nodes
"""
    function NodesAndModes.equi_nodes(elem::PhysicalFrame, curve, N)

Returns back `Np(N)` equally spaced nodes on the background quadrilateral corresponding 
to `elem`, with points inside of `curve` removed.
"""
function NodesAndModes.equi_nodes(elem::PhysicalFrame{2}, curve, N)
    r, s = equi_nodes(Quad(), N)
    x, y = map_nodes_to_cutcell_boundingbox(elem, r, s)
    ids = .!PathIntersections.is_contained.(curve, zip(x, y))
    return x[ids], y[ids]
end

function NodesAndModes.equi_nodes(elem::PhysicalFrame{2}, 
                                  curves::Union{<:Tuple, <:AbstractArray}, N)
    r, s = equi_nodes(Quad(), N)
    x, y = map_nodes_to_cutcell_boundingbox(elem, r, s)
    ids = .!PathIntersections.is_contained.(first(curves), zip(x, y))
    for curve in Base.tail(curves)
        ids = ids .&& .!PathIntersections.is_contained.(curve, zip(x, y))
    end
    return x[ids], y[ids]
end

function triangulate_points(coordinates::AbstractMatrix)
    triin=Triangulate.TriangulateIO()
    triin.pointlist = coordinates
    triout, _ = triangulate("Q", triin)
    VX, VY = (triout.pointlist[i,:] for i = 1:size(triout.pointlist,1))
    EToV = permutedims(triout.trianglelist)
    return (VX, VY), EToV
end

"""
    caratheodory_pruning_qr(V, w_in)

This performs Caratheodory pruning using a naive QR-based algorithm. 
Returns (w, inds), where `inds` denotes sub-selected indices for a 
reduced quadrature rule, and `w` is a vector of reduced positive weights.

The original Matlab code this was based on was authored by Akil Narayan.
""" 
function caratheodory_pruning_qr(V, w_in)

    if length(w_in) <= size(V, 2)
        return w_in, eachindex(w_in)
    end
    w = copy(w_in)
    M, N = size(V)
    inds = collect(1:M)
    m = M-N
    Q, _ = qr(V)
    for _ in 1:m
        kvec = view(Q, :, size(Q, 2))

        # for subtracting the kernel vector
        idp = findall(@. kvec > 0)
        alphap, k0p = findmin(view(w, inds[idp]) ./ view(kvec, idp))
        k0p = idp[k0p]
    
        # for adding the kernel vector
        idn = findall(@. kvec < 0)
        alphan, k0n = findmax(view(w, inds[idn]) ./ view(kvec, idn))
        k0n = idn[k0n]
    
        alpha, k0 = abs(alphan) < abs(alphap) ? (alphan, k0n) : (alphap, k0p)
        @. w[inds] = w[inds] - alpha * kvec
        deleteat!(inds, k0)
        Q, _ = qr(V[inds, :])
    end
    return w, inds
end