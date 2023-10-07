module StartUpDGSummationByPartsOperatorsExt

using LinearAlgebra: LinearAlgebra, Diagonal, diag, norm, UniformScaling
using SparseArrays: sparse, droptol!, spzeros

using StartUpDG

# Required for visualization code
if isdefined(Base, :get_extension)
    using SummationByPartsOperators:
        SummationByPartsOperators,
        DerivativeOperator,
        grid,
        AbstractDerivativeOperator,
        AbstractNonperiodicDerivativeOperator,
        PeriodicDerivativeOperator,
        AbstractPeriodicDerivativeOperator
else
    # Until Julia v1.9 is the minimum required version for Trixi.jl, we still support Requires.jl
    using ..SummationByPartsOperators
    using ..SummationByPartsOperators:
        AbstractDerivativeOperator,
        AbstractPeriodicDerivativeOperator,
        AbstractNonperiodicDerivativeOperator,
        DerivativeOperator,
        PeriodicDerivativeOperator,
        grid
end

function construct_1d_operators(D::AbstractDerivativeOperator, tol)
    nodes_1d = collect(grid(D))
    M = SummationByPartsOperators.mass_matrix(D)
    if M isa UniformScaling
        weights_1d = M * ones(Bool, length(nodes_1d))
    else
        weights_1d = diag(M)
    end

    # StartUpDG assumes nodes from -1 to +1. Thus, we need to re-scale everything.
    # We can adjust the grid spacing as follows.
    xmin = SummationByPartsOperators.xmin(D)
    xmax = SummationByPartsOperators.xmax(D)
    factor = 2 / (xmax - xmin)
    @. nodes_1d = factor * (nodes_1d - xmin) - 1
    @. weights_1d = factor * weights_1d

    D_1d = droptol!(inv(factor) * sparse(D), tol)
    I_1d = Diagonal(ones(Bool, length(nodes_1d)))

    return nodes_1d, weights_1d, D_1d, I_1d
end

function StartUpDG.RefElemData(element_type::Line,
    D::AbstractDerivativeOperator;
    tol = 100 * eps(),)
    approximation_type = D
    N = SummationByPartsOperators.accuracy_order(D) # kind of polynomial degree

    # 1D operators
    nodes_1d, weights_1d, D_1d = construct_1d_operators(D, tol)

    # volume
    rq = r = nodes_1d
    wq = weights_1d
    Dr = D_1d
    M = Diagonal(wq)
    Pq = LinearAlgebra.I
    Vq = LinearAlgebra.I

    VDM = nothing # unused generalized Vandermonde matrix

    rst = (r,)
    rstq = (rq,)
    Drst = (Dr,)

    # face
    face_vertices = StartUpDG.face_vertices(element_type)
    face_mask = [1, length(nodes_1d)]

    rf = [-1.0; 1.0]
    nrJ = [-1.0; 1.0]
    wf = [1.0; 1.0]
    if D isa AbstractPeriodicDerivativeOperator
        # we do not need any face stuff for periodic operators
        Vf = spzeros(length(wf), length(wq))
    else
        Vf = sparse([1, 2], [1, length(nodes_1d)], [1.0, 1.0])
    end
    LIFT = Diagonal(wq) \ (Vf' * Diagonal(wf))

    rstf = (rf,)
    nrstJ = (nrJ,)

    # low order interpolation nodes
    r1 = StartUpDG.nodes(element_type, 1)
    V1 = StartUpDG.vandermonde(element_type, 1, r) /
         StartUpDG.vandermonde(element_type, 1, r1)

    return RefElemData(element_type,
        approximation_type,
        N,
        face_vertices,
        V1,
        rst,
        VDM,
        face_mask,
        rst,
        LinearAlgebra.I, # plotting
        rstq,
        wq,
        Vq, # quadrature
        rstf,
        wf,
        Vf,
        nrstJ, # faces
        M,
        Pq,
        Drst,
        LIFT)
end

function StartUpDG.RefElemData(element_type::Quad,
    D::AbstractDerivativeOperator;
    tol = 100 * eps(),)
    approximation_type = D
    N = SummationByPartsOperators.accuracy_order(D) # kind of polynomial degree

    # 1D operators
    nodes_1d, weights_1d, D_1d, I_1d = construct_1d_operators(D, tol)

    # volume
    s, r = vec.(StartUpDG.NodesAndModes.meshgrid(nodes_1d)) # this is to match
    # ordering of nrstJ
    rq = r
    sq = s
    wr, ws = vec.(StartUpDG.NodesAndModes.meshgrid(weights_1d))
    wq = wr .* ws
    Dr = kron(I_1d, D_1d)
    Ds = kron(D_1d, I_1d)
    M = Diagonal(wq)
    Pq = LinearAlgebra.I
    Vq = LinearAlgebra.I

    VDM = nothing # unused generalized Vandermonde matrix

    rst = (r, s)
    rstq = (rq, sq)
    Drst = (Dr, Ds)

    # face
    face_vertices = StartUpDG.face_vertices(element_type)
    face_mask = vcat(StartUpDG.find_face_nodes(element_type, r, s)...)

    rf, sf, wf, nrJ, nsJ = StartUpDG.init_face_data(element_type,
        quad_rule_face = (nodes_1d, weights_1d))
    if D isa AbstractPeriodicDerivativeOperator
        # we do not need any face stuff for periodic operators
        Vf = spzeros(length(wf), length(wq))
    else
        Vf = sparse(eachindex(face_mask), face_mask, ones(Bool, length(face_mask)))
    end
    LIFT = Diagonal(wq) \ (Vf' * Diagonal(wf))

    rstf = (rf, sf)
    nrstJ = (nrJ, nsJ)

    # low order interpolation nodes
    r1, s1 = StartUpDG.nodes(element_type, 1)
    V1 = StartUpDG.vandermonde(element_type, 1, r, s) /
         StartUpDG.vandermonde(element_type, 1, r1, s1)

    return RefElemData(element_type,
        approximation_type,
        N,
        face_vertices,
        V1,
        rst,
        VDM,
        face_mask,
        rst,
        LinearAlgebra.I, # plotting
        rstq,
        wq,
        Vq, # quadrature
        rstf,
        wf,
        Vf,
        nrstJ, # faces
        M,
        Pq,
        Drst,
        LIFT)
end

function StartUpDG.RefElemData(element_type::Hex,
    D::AbstractDerivativeOperator;
    tol = 100 * eps(),)
    approximation_type = D
    N = SummationByPartsOperators.accuracy_order(D) # kind of polynomial degree

    # 1D operators
    nodes_1d, weights_1d, D_1d, I_1d = construct_1d_operators(D, tol)

    # volume
    # to match ordering of nrstJ
    s, r, t = vec.(StartUpDG.NodesAndModes.meshgrid(nodes_1d, nodes_1d, nodes_1d))
    rq = r
    sq = s
    tq = t
    wr, ws, wt = vec.(StartUpDG.NodesAndModes.meshgrid(weights_1d, weights_1d, weights_1d))
    wq = wr .* ws .* wt
    Dr = kron(I_1d, I_1d, D_1d)
    Ds = kron(I_1d, D_1d, I_1d)
    Dt = kron(D_1d, I_1d, I_1d)
    M = Diagonal(wq)
    Pq = LinearAlgebra.I
    Vq = LinearAlgebra.I

    VDM = nothing # unused generalized Vandermonde matrix

    rst = (r, s, t)
    rstq = (rq, sq, tq)
    Drst = (Dr, Ds, Dt)

    # face
    face_vertices = StartUpDG.face_vertices(element_type)
    face_mask = vcat(StartUpDG.find_face_nodes(element_type, r, s, t)...)

    rf, sf, tf, wf, nrJ, nsJ, ntJ = let
        rf, sf = vec.(StartUpDG.NodesAndModes.meshgrid(nodes_1d, nodes_1d))
        wr, ws = vec.(StartUpDG.NodesAndModes.meshgrid(weights_1d, weights_1d))
        wf = wr .* ws
        StartUpDG.init_face_data(element_type, quad_rule_face = (rf, sf, wf))
    end
    Vf = sparse(eachindex(face_mask), face_mask, ones(Bool, length(face_mask)))
    LIFT = Diagonal(wq) \ (Vf' * Diagonal(wf))

    rstf = (rf, sf, tf)
    nrstJ = (nrJ, nsJ, ntJ)

    # low order interpolation nodes
    r1, s1, t1 = StartUpDG.nodes(element_type, 1)
    V1 = StartUpDG.vandermonde(element_type, 1, r, s, t) /
         StartUpDG.vandermonde(element_type, 1, r1, s1, t1)

    return RefElemData(element_type,
        approximation_type,
        N,
        face_vertices,
        V1,
        rst,
        VDM,
        face_mask,
        rst,
        LinearAlgebra.I, # plotting
        rstq,
        wq,
        Vq, # quadrature
        rstf,
        wf,
        Vf,
        nrstJ, # faces
        M,
        Pq,
        Drst,
        LIFT)
end

# specialized Hex constructor in 3D to reduce memory usage.
function StartUpDG.RefElemData(element_type::Hex,
    D::AbstractPeriodicDerivativeOperator;
    tol = 100 * eps(),)
    approximation_type = D
    N = SummationByPartsOperators.accuracy_order(D) # kind of polynomial degree

    # 1D operators
    nodes_1d, weights_1d, D_1d, I_1d = construct_1d_operators(D, tol)

    # volume
    # to match ordering of nrstJ
    s, r, t = vec.(StartUpDG.NodesAndModes.meshgrid(nodes_1d, nodes_1d, nodes_1d))
    rq = r
    sq = s
    tq = t
    wr, ws, wt = vec.(StartUpDG.NodesAndModes.meshgrid(weights_1d, weights_1d, weights_1d))
    wq = wr .* ws .* wt
    Dr = kron(I_1d, I_1d, D_1d)
    Ds = kron(I_1d, D_1d, I_1d)
    Dt = kron(D_1d, I_1d, I_1d)
    M = Diagonal(wq)
    Pq = LinearAlgebra.I
    Vq = LinearAlgebra.I

    VDM = nothing # unused generalized Vandermonde matrix

    rst = (r, s, t)
    rstq = (rq, sq, tq)
    Drst = (Dr, Ds, Dt)

    # face
    # We do not need any face data for periodic operators. Thus, we just
    # pass `nothing` to save memory.
    face_vertices = ntuple(_ -> nothing, 3)
    face_mask = nothing
    wf = nothing
    rstf = ntuple(_ -> nothing, 3)
    nrstJ = ntuple(_ -> nothing, 3)
    Vf = nothing
    LIFT = nothing

    # low order interpolation nodes
    V1 = nothing # do not need to store V1, since we specialize StartUpDG.MeshData to avoid using it.

    return RefElemData(element_type,
        approximation_type,
        N,
        face_vertices,
        V1,
        rst,
        VDM,
        face_mask,
        rst,
        LinearAlgebra.I, # plotting
        rstq,
        wq,
        Vq, # quadrature
        rstf,
        wf,
        Vf,
        nrstJ, # faces
        M,
        Pq,
        Drst,
        LIFT)
end

function Base.show(io::IO,
    mime::MIME"text/plain",
    rd::RefElemData{NDIMS, ElementType, ApproximationType}) where {
    NDIMS,
    ElementType <: StartUpDG.AbstractElemShape,
    ApproximationType <: AbstractDerivativeOperator,
}
    @nospecialize rd
    print(io, "RefElemData for an approximation using an ")
    show(IOContext(io, :compact => true), rd.approximation_type)
    print(io, " on $(rd.element_type) element")
end

function Base.show(io::IO,
    rd::RefElemData{NDIMS, ElementType, ApproximationType}) where {
    NDIMS,
    ElementType <: StartUpDG.AbstractElemShape,
    ApproximationType <: AbstractDerivativeOperator,
}
    @nospecialize rd
    print(io, "RefElemData{", summary(rd.approximation_type), ", ", rd.element_type, "}")
end

function StartUpDG.inverse_trace_constant(rd::RefElemData{
    NDIMS,
    ElementType,
    ApproximationType,
}) where {
    NDIMS,
    ElementType <: Union{Line, Quad, Hex},
    ApproximationType <: AbstractDerivativeOperator,
}
    D = rd.approximation_type

    # the inverse trace constant is the maximum eigenvalue corresponding to
    #       M_f * v = λ * M * v
    # where M_f is the face mass matrix and M is the volume mass matrix.
    # Since M is diagonal and since M_f is just the boundary "mask" matrix
    # (which extracts the first and last entries of a vector), the maximum
    # eigenvalue is the inverse of the first or last mass matrix diagonal.
    left_weight = SummationByPartsOperators.left_boundary_weight(D)
    right_weight = SummationByPartsOperators.right_boundary_weight(D)
    max_eigenvalue = max(inv(left_weight), inv(right_weight))

    # For tensor product elements, the trace constant for higher dimensional
    # elements is the one-dimensional trace constant multiplied by `NDIMS`. See
    #     "GPU-accelerated discontinuous Galerkin methods on hybrid meshes."
    #     Chan, Jesse, et al (2016), https://doi.org/10.1016/j.jcp.2016.04.003
    # for more details (specifically, Appendix A.1, Theorem A.4).
    return NDIMS * max_eigenvalue
end

# This is used in `estimate_dt`. `estimate_h` uses that `Jf / J = O(h^{NDIMS-1}) / O(h^{NDIMS}) = O(1/h)`.
# However, since we do not initialize `Jf` for periodic FDSBP operators, we specialize `estimate_h`
# based on the reference grid provided by SummationByPartsOperators.jl and information about the domain size
# provided by `md::MeshData``.
function StartUpDG.estimate_h(e,
    rd::RefElemData{NDIMS, ElementType, ApproximationType},
    md::MeshData) where {
    NDIMS,
    ElementType <: StartUpDG.AbstractElemShape,
    ApproximationType <: SummationByPartsOperators.AbstractPeriodicDerivativeOperator,
}
    D = rd.approximation_type
    x = grid(D)

    # we assume all SummationByPartsOperators.jl reference grids are rescaled to [-1, 1]
    xmin = SummationByPartsOperators.xmin(D)
    xmax = SummationByPartsOperators.xmax(D)
    factor = 2 / (xmax - xmin)

    # If the domain has size L^NDIMS, then `minimum(md.J)^(1 / NDIMS) = L`.
    # WARNING: this is not a good estimate on anisotropic grids.
    return minimum(diff(x)) * factor * minimum(md.J)^(1 / NDIMS)
end

end
