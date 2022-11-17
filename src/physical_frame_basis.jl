# dimension of a cut cell polynomial space
@inline Np_cut(N) = (N + 1) * (N + 2) ÷ 2 

@inline mean(x) = sum(x) / length(x)

"""
    `PhysicalFrame <: AbstractElemShape`
    
`PhysicalFrame` element type. Uses a total degree N approximation space, but is 
computed with a tensor product Legendre basis as opposed to a triangular PKDO basis.
Stores fields `shifting` and `scaling` to shift/scale physical coordinates so that 
they are on the reference element. 

    PhysicalFrame()
    PhysicalFrame(x, y)

Constructors for a PhysicalFrame object (optionally uses arrays of points `x`, `y` on a cut element).
"""
struct PhysicalFrame{Shifting <: Union{<:SVector, <:NTuple}, Scaling <: Union{<:SVector, <:NTuple}, VXYZ} <: AbstractElemShape
    shifting::Shifting
    scaling::Scaling
    vxyz::VXYZ # coordinates of background Cartesian cell
end

# default shifting and scaling, restricted to 2D for now
function PhysicalFrame()     
    return PhysicalFrame(SVector(0.0, 0.0), SVector(1.0, 1.0), 
                         (SVector(-1., 1.), SVector(-1., 1.))) # set background coordinates to reference cell
end

function PhysicalFrame(x, y)
    shifting = (mean(x), mean(y))
    scaling = map(x -> 2 / (x[2] - x[1]), (extrema(x), extrema(y)))
    return PhysicalFrame(shifting, scaling, nothing)
end

function PhysicalFrame(x, y, vx, vy)
    shifting = (mean(x), mean(y))
    scaling = map(x -> 2 / (x[2] - x[1]), (extrema(x), extrema(y)))
    return PhysicalFrame(shifting, scaling, (vx, vy))
end

function shift_and_scale(elem::PhysicalFrame, x, y)
    @unpack shifting, scaling = elem
    r = @. (x - shifting[1]) * scaling[1]
    s = @. (y - shifting[2]) * scaling[2]
    return r, s
end

function NodesAndModes.basis(elem::PhysicalFrame, N, x, y)
    Np = (N + 1) * (N + 2) ÷ 2

    r, s = shift_and_scale(elem, x, y)

    @unpack scaling = elem

    sk = 1
    V, Vr, Vs = ntuple(x->zeros(length(r), Np), 3)
    for j = 0:N
        P_j = jacobiP(s, 0, 0, j)
        for i = 0:N-j
            P_i = jacobiP(r, 0, 0, i)
            V[:, sk]  = P_i .* P_j
            Vr[:, sk] = grad_jacobiP(r, 0, 0, i) .* P_j * scaling[1]
            Vs[:, sk] = P_i .* grad_jacobiP(s, 0, 0, j) * scaling[2]
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
function NodesAndModes.equi_nodes(elem::PhysicalFrame, curve, N)
    @unpack vxyz = elem
    r, s = equi_nodes(Quad(), N)
    x, y = map_nodes_to_background_cell(elem, r, s)
    ids = .!is_contained.(curve, zip(x, y))
    return x[ids], y[ids]
end

function map_nodes_to_background_cell(elem::PhysicalFrame, r, s)
    @unpack vxyz = elem
    vx, vy = vxyz
    dx, dy = diff(vx), diff(vy)
    x = @. 0.5 * (1 + r) * dx + vx[1]
    y = @. 0.5 * (1 + s) * dy + vy[1]
    return x, y
end