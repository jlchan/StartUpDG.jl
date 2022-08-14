@inline mean(x) = sum(x) / length(x)

"""
    `PhysicalFrame <: AbstractElemShape`

`PhysicalFrame` element type. Uses a total degree N approximation space, but is 
computed with a tensor product Legendre basis as opposed to a triangular PKDO basis.
"""
struct PhysicalFrame{shifting <: Union{<:SVector, <:NTuple}, scaling <: Union{<:SVector, <:NTuple}} <: AbstractElemShape
    shifting::shifting
    scaling::scaling
end

"""
    PhysicalFrame()
    PhysicalFrame(CutCell(), x, y)

Constructor for a PhysicalFrame object given arrays of points `x`, `y` on a cut element. 
"""
# default shiftinging and scaling
PhysicalFrame() = PhysicalFrame(SVector(0.0, 0.0), SVector(1.0, 1.0))

function PhysicalFrame(x, y)
    shifting = (mean(x), mean(y))
    scaling = map(x -> 2 / (x[2] - x[1]), (extrema(x), extrema(y)))
    return PhysicalFrame(shifting, scaling)
end

function NodesAndModes.vandermonde(elem::PhysicalFrame, N, x, y)
    Np = (N + 1) * (N + 2) ÷ 2

    @unpack shifting, scaling = elem
    r = @. (x - shifting[1]) * scaling[1]
    s = @. (y - shifting[2]) * scaling[2]

    sk = 1
    # V, Vr, Vs = ntuple(x->zeros(length(r), Np), 3)
    V = zeros(length(r), Np)
    for j = 0:N
        P_j = jacobiP(s, 0, 0, j)
        for i = 0:N-j
            P_i = jacobiP(r, 0, 0, i)
            V[:, sk]  = P_i .* P_j
            # Vr[:, sk] = grad_jacobiP(r, 0, 0, i) .* P_j * scaling[1]
            # Vs[:, sk] = P_i .* grad_jacobiP(s, 0, 0, j) * scaling[2]
            sk += 1
        end
    end
    return V
    # return V, Vr, Vs
end