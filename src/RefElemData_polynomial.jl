# The following functions determine what default quadrature type to use for each element.
# For tensor product elements, we default to TensorProductQuadrature. 
# For simplices, wedges, and pyramids, we default to MultidimensionalQuadrature

# simplices and pyramids default to multidimensional quadrature
RefElemData(elem::Union{Line, Tri, Tet, Wedge, Pyr}, 
            approx_type::Polynomial{DefaultPolynomialType}, 
            N; kwargs...) = 
    RefElemData(elem, Polynomial{MultidimensionalQuadrature}(), N; kwargs...)

# on quad and hex elements, default to a tensor product quadrature 
RefElemData(elem::Union{Quad, Hex}, approximation_type::Polynomial{DefaultPolynomialType}, N; kwargs...) = 
    RefElemData(elem, Polynomial(TensorProductQuadrature(gauss_quad(0, 0, N))), N; kwargs...)

# special case: for lines, tensor product and multidimensional quadrature are the same
RefElemData(elem::Line, approx_type::Polynomial{<:TensorProductQuadrature}, N; kwargs...) = 
    RefElemData(elem, Polynomial{MultidimensionalQuadrature}(), N; 
                quad_rule_vol=approx_type.data.quad_rule_1D, kwargs...)

"""
    RefElemData(elem::Line, approximation_type, N;
                quad_rule_vol = quad_nodes(elem, N+1))
    RefElemData(elem, approximation_type, N;
                quad_rule_vol = quad_nodes(elem, N),
                quad_rule_face = quad_nodes(face_type(elem), N))

Constructor for `RefElemData` for different element types.
"""
function RefElemData(elem::Line, 
                     approx_type::Polynomial{MultidimensionalQuadrature},
                     N; quad_rule_vol=quad_nodes(elem, N+1), Nplot=10)

    fv = face_vertices(elem)

    # reference element nodes
    r = nodes(elem, N)
    Fmask = [1 N+1]

    # compute operators
    VDM = vandermonde(elem, N, r)
    Dr = grad_vandermonde(elem, N, r)/VDM
    V1 = vandermonde(elem, 1, r) / vandermonde(elem, 1, nodes(elem, 1))

    # quadrature operators
    rq, wq = quad_rule_vol
    rf, wf  = [-1.0; 1.0], [1.0; 1.0]
    nrJ = [-1.0; 1.0]
    Vq = vandermonde(elem, N, rq) / VDM
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))
    Vf = vandermonde(elem, N, rf) / VDM
    LIFT = M \ (Vf') # lift matrix

    # plotting nodes
    rp = equi_nodes(elem, Nplot) 
    Vp = vandermonde(elem, N, rp) / VDM

    return RefElemData(elem, approx_type, N, fv, V1,
                       tuple(r), VDM, vec(Fmask),
                       tuple(rp), Vp,
                       tuple(rq), wq, Vq,
                       tuple(rf), wf, Vf, tuple(nrJ),
                       M, Pq, tuple(Dr), LIFT)
end


function RefElemData(elem::Union{Tri, Quad},
                     approx_type::Polynomial{MultidimensionalQuadrature}, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face=quad_nodes(face_type(elem), N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle

    # Construct matrices on reference elements
    r, s = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s)...)

    # low order interpolation nodes
    r1, s1 = nodes(elem, 1) 
    V1 = vandermonde(elem, 1, r, s) / vandermonde(elem, 1, r1, s1)

    # differentiation operators
    VDM, Vr, Vs = basis(elem, N, r, s)
    Dr = Vr / VDM
    Ds = Vs / VDM

    # quadrature nodes
    rq, sq, wq = quad_rule_vol
    rf, sf, wf, nrJ, nsJ = init_face_data(elem; quad_rule_face)

    # quadrature-based operators
    Vq = vandermonde(elem, N, rq, sq) / VDM
    Vf = vandermonde(elem, N, rf, sf) / VDM # interpolates from nodes to face nodes
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))
    LIFT = M \ (Vf' * diagm(wf)) # lift matrix used in rhs evaluation

    rp, sp = equi_nodes(elem, Nplot) # plotting nodes
    Vp = vandermonde(elem, N, rp, sp) / VDM

    return RefElemData(elem, approx_type, N, fv, V1,
                       tuple(r, s), VDM, vec(Fmask),
                       tuple(rp, sp), Vp,
                       tuple(rq, sq), wq, Vq,
                       tuple(rf, sf), wf, Vf, tuple(nrJ, nsJ),
                       M, Pq, (Dr, Ds), LIFT)
end

function RefElemData(elem::Union{Tet, Hex},
                     approx_type::Polynomial{MultidimensionalQuadrature}, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face=quad_nodes(face_type(elem), N),
                     Nplot=10)

    if elem isa Hex && N > 4
        @warn "Since N > 4, we suggest using `RefElemData(Hex(), Polynomial(TensorProductQuadrature(gauss_quad(0, 0, $N+1))), $N)`, " * 
              "which is more efficient."
    end      

    fv = face_vertices(elem) 

    # Construct matrices on reference elements
    r, s, t = nodes(elem, N)
    tol = 1e2 * eps() * length(r) # loosen the tolerance if N >> 1
    Fmask = hcat(find_face_nodes(elem, r, s, t, tol)...)
    VDM, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> A / VDM).((Vr, Vs, Vt))

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)

    # Nodes on faces, and face node coordinate
    rq, sq, tq, wq = quad_rule_vol
    rf, sf, tf, wf, nrJ, nsJ, ntJ = init_face_data(elem; quad_rule_face)

    # quadrature operators
    Vq = vandermonde(elem, N, rq, sq, tq) / VDM
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))
    Vf = vandermonde(elem, N, rf, sf, tf) / VDM
    LIFT = M \ (Vf' * diagm(wf))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp, tp) / VDM

    return RefElemData(elem, approx_type, N, fv, V1,
                       tuple(r, s, t), VDM, vec(Fmask),
                       tuple(rp, sp, tp), Vp,
                       tuple(rq, sq, tq), wq, Vq,
                       tuple(rf, sf, tf), wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, (Dr, Ds, Dt), LIFT)
end

"""
    RefElemData(elem::Wedge, approximation_type::Polynomial, N;
                quad_rule_vol=quad_nodes(elem, N),
                quad_rule_face_quad=quad_nodes(Quad(), N), 
                quad_rule_face_tri=quad_nodes(Tri(), N), 
                quad_rule_face=(quad_rule_face_quad, quad_rule_face_tri),
                Nplot=10)

Builds operators for prisms/wedges
"""
function RefElemData(elem::Wedge, 
                     approximation_type::Polynomial{MultidimensionalQuadrature}, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face_quad=quad_nodes(Quad(), N), 
                     quad_rule_face_tri=quad_nodes(Tri(), N), 
                     quad_rule_face=(quad_rule_face_quad, quad_rule_face_tri),
                     Nplot=10)

    #Find the vertices of the faces
    fv = face_vertices(elem)

    #Get interpolation nodes of degree N 
    r, s, t = nodes(elem, N)
    
    VDM, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> A / VDM).((Vr, Vs, Vt))
    Drst = (Dr, Ds, Dt)
    
    Fmask = find_face_nodes(elem, r, s, t)

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)
    
    # build face quadrature nodes
    rquad, squad, wquad = quad_rule_face[1]
    rtri, stri, wtri = quad_rule_face[2]
    rf = vcat( rquad,         -rquad,   -one.(wquad),  rtri      , rtri)
    sf = vcat(-one.(wquad),   rquad,    rquad      ,  stri      , stri)
    tf = vcat( squad,         squad,    squad      , -one.(wtri), one.(wtri))
    wf = vcat(wquad, wquad, wquad, wtri, wtri)

    # index into the face nodes     
    quad_face_ids(f) = (1:length(wquad)) .+ (f-1) * length(wquad)
    tri_face_ids(f) = (1:length(wtri)) .+ (f-1) * length(wtri) .+ 3 * length(wquad)
    node_ids_by_face = (quad_face_ids(1), quad_face_ids(2), quad_face_ids(3), 
                        tri_face_ids(1), tri_face_ids(2))

    rstf = tuple(rf, sf, tf)
    Vf = vandermonde(elem, N, rf, sf, tf) / VDM
        
    # for nrJ and nsJ normal on face 1-3 coincide with the triangular normals
    zt, zq = zeros(length(wtri)), zeros(length(wquad))
    et, eq = ones(length(wtri)), ones(length(wquad))

    nrJ = [zq; eq; -eq; zt; zt]
    nsJ = [-eq; eq; zq; zt; zt]
    ntJ = [zq; zq; zq; -et; et]

    rq, sq, tq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq, sq, tq) / VDM
    M  = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp, tp) / VDM

    LIFT = M \ (Vf' * diagm(wf))

    return RefElemData(Wedge(node_ids_by_face), approximation_type, N, fv, V1,
                       tuple(r, s, t), VDM, Fmask,
                       tuple(rp, sp, tp), Vp,
                       tuple(rq, sq, tq), wq, Vq,
                       rstf, wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, Drst, LIFT)
end

"""
    RefElemData(elem::Pyr, 
                approximation_type::Polynomial, N;
                quad_rule_vol=quad_nodes(elem, N),
                quad_rule_face_quad=quad_nodes(Quad(), N), 
                quad_rule_face_tri=quad_nodes(Tri(), N), 
                quad_rule_face=(quad_rule_face_quad, quad_rule_face_tri),
                Nplot=10)

Builds operators for pyramids.
"""
function RefElemData(elem::Pyr, 
                     approximation_type::Polynomial{MultidimensionalQuadrature}, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face_quad=quad_nodes(Quad(), N), 
                     quad_rule_face_tri=quad_nodes(Tri(), N), 
                     quad_rule_face=(quad_rule_face_quad, quad_rule_face_tri),
                     Nplot=10)

    #Find the vertices of the faces
    fv = face_vertices(elem)

    #Get interpolation nodes of degree N 
    r, s, t = nodes(elem, N)        
    VDM = vandermonde(elem, N, r, s, t)    
    Fmask = find_face_nodes(elem, r, s, t)

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)
    
    # build face quadrature nodes
    rquad, squad, wquad = quad_rule_face[1]
    rtri, stri, wtri = quad_rule_face[2]
    sf = vcat(rtri,         rtri,   -one.(wtri), -stri,     squad)
    rf = vcat(-one.(wtri), -stri,   rtri,         rtri,     rquad)
    tf = vcat(stri,         stri,   stri,         stri,     -one.(wquad))
    wf = vcat(wtri, wtri, wtri, wtri, wquad)

    # Index into the face nodes. 
    # Faces are ordered tri faces (+/-r, then +/-s), then the quad face
    tri_face_ids(f) = (1:length(wtri)) .+ (f-1) * length(wtri) 
    quad_face_ids = (1:length(wquad)) .+ 4 * length(wtri)
    node_ids_by_face = (tri_face_ids(1), tri_face_ids(2), 
                        tri_face_ids(3), tri_face_ids(4), 
                        quad_face_ids)

    rstf = tuple(rf, sf, tf)
    Vf = vandermonde(elem, N, rf, sf, tf) / VDM
        
    # for nrJ and nsJ normal on face 1-3 coincide with the triangular normals
    zt, zq = zeros(length(wtri)), zeros(length(wquad))
    et, eq = ones(length(wtri)),  ones(length(wquad))

    nrJ = [-et;  et;  zt;  zt;  zq]
    nsJ = [ zt;  zt; -et;  et;  zq]
    ntJ = [ zt;  et;  zt;  et; -eq]

    rq, sq, tq, wq = quad_rule_vol
    Vq, Vrq, Vsq, Vtq = map(A -> A / VDM, basis(elem, N, rq, sq, tq))
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    # We define nodal differentiation matrices using quadrature instead of using 
    # interpolation nodes and linear algebra. This is because the basis is not 
    # polynomial, so derivative(pyramid_space) ∉ pyramid_space. Instead, we define 
    # a weak gradient `D` s.t. for any u ∈ pyramid_space, we have:
    #       (Du, v) = (du/dx, v) ∀v ∈ pyramid_space
    Drst = map(Vderiv -> M \ (Vq' * diagm(wq) * Vderiv), (Vrq, Vsq, Vtq))

    LIFT = M \ (Vf' * diagm(wf))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp, tp) / VDM

    return RefElemData(Pyr(node_ids_by_face), approximation_type, N, fv, V1,
                       tuple(r, s, t), VDM, Fmask,
                       tuple(rp, sp, tp), Vp,
                       tuple(rq, sq, tq), wq, Vq,
                       rstf, wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, Drst, LIFT)
end

"""
    RefElemData(elem::Union{Quad, Hex}, approximation_type::TensorProductQuadrature, N)
    RefElemData(elem::Union{Quad, Hex}, approximation_type::Polynomial{<:TensorProductQuadrature}, N)

Constructor for quadrilateral and hexahedral `RefElemData` where the quadrature is assumed to have a 
tensor product structure. 
"""

function RefElemData(elem::Quad, 
                     approximation_type::Polynomial{<:TensorProductQuadrature}, N;
                     quad_rule_face = approximation_type.data.quad_rule_1D,
                     Nplot = 10)

    fv = face_vertices(elem) 

    # Construct matrices on reference elements
    r, s = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s)...)

    # construct 1D operator for faster Kronecker solves
    r1D = nodes(Line(), N)
    rq1D, wq1D = approximation_type.data.quad_rule_1D
    VDM_1D = vandermonde(Line(), N, r1D)
    Vq1D = vandermonde(Line(), N, rq1D) / VDM_1D
    invVDM_1D = inv(VDM_1D)
    invM_1D = VDM_1D * VDM_1D'
    M1D = Vq1D' * diagm(wq1D) * Vq1D

    # form kronecker products of multidimensional matrices to invert/multiply
    VDM = kron(VDM_1D, VDM_1D)
    invVDM = kron(invVDM_1D, invVDM_1D)
    invM = kron(invM_1D, invM_1D)

    M = kron(M1D, M1D)

    _, Vr, Vs = basis(elem, N, r, s)
    Dr, Ds = (A -> A * invVDM).((Vr, Vs))

    # low order interpolation nodes
    r1, s1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s) / vandermonde(elem, 1, r1, s1)

    # Nodes on faces, and face node coordinate
    rf, sf, wf, nrJ, nsJ = init_face_data(elem, quad_rule_face=quad_rule_face)

    # quadrature nodes - build from 1D nodes.    
    rq, sq, wq = tensor_product_quadrature(elem, approximation_type.data.quad_rule_1D...)
    Vq = kron(Vq1D, Vq1D) # vandermonde(elem, N, rq, sq, tq) * invVDM
    Pq = invM * (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf) * invVDM
    LIFT = invM * (Vf' * diagm(wf))

    # plotting nodes
    rp1D = LinRange(-1, 1, Nplot + 1)
    Vp1D = vandermonde(Line(), N, rp1D) / VDM_1D
    Vp = kron(Vp1D, Vp1D)
    rp, sp = vec.(StartUpDG.NodesAndModes.meshgrid(rp1D, rp1D))
    
    return RefElemData(elem, approximation_type, N, fv, V1,
                       tuple(r, s), VDM, vec(Fmask),
                       tuple(rp, sp), Vp, 
                       tuple(rq, sq), wq, Vq,
                       tuple(rf, sf), wf, Vf, tuple(nrJ, nsJ),
                       M, Pq, (Dr, Ds), LIFT)
end

function RefElemData(elem::Hex, 
                     approximation_type::Polynomial{<:TensorProductQuadrature}, N;
                     quad_rule_face =
                        tensor_product_quadrature(face_type(elem), 
                                                  approximation_type.data.quad_rule_1D...),
                     Nplot = 10)

    fv = face_vertices(elem) 

    # Construct matrices on reference elements
    r, s, t = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s, t)...)

    # construct 1D operator for faster Kronecker solves
    r1D = nodes(Line(), N)
    rq1D, wq1D = approximation_type.data.quad_rule_1D
    VDM_1D = vandermonde(Line(), N, r1D)
    Vq1D = vandermonde(Line(), N, rq1D) / VDM_1D
    invVDM_1D = inv(VDM_1D)
    invM_1D = VDM_1D * VDM_1D'
    M1D = Vq1D' * diagm(wq1D) * Vq1D

    # form kronecker products of multidimensional matrices to invert/multiply
    # use dense matrix "kron" if N is small.
    # use memory-saving "kronecker" if N is large.
    build_kronecker_product = (N > 10) ? kronecker : kron

    VDM = build_kronecker_product(VDM_1D, VDM_1D, VDM_1D)
    invM = build_kronecker_product(invM_1D, invM_1D, invM_1D)

    # always use the more efficient Kronecker product to compute invVDM 
    # since we multiply by it to compute Dr, Ds, Dt
    invVDM = kronecker(invVDM_1D, invVDM_1D, invVDM_1D) 
    
    M = build_kronecker_product(M1D, M1D, M1D)
    
    _, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> Matrix(A * invVDM)).((Vr, Vs, Vt))

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)

    # Nodes on faces, and face node coordinate
    rf, sf, tf, wf, nrJ, nsJ, ntJ = init_face_data(elem, quad_rule_face=quad_rule_face)

    # quadrature nodes - build from 1D nodes.    
    rq, sq, tq, wq = tensor_product_quadrature(elem, approximation_type.data.quad_rule_1D...)
    Vq = build_kronecker_product(Vq1D, Vq1D, Vq1D) # vandermonde(elem, N, rq, sq, tq) * invVDM
    Pq = invM * (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf, tf) * invVDM
    LIFT = invM * (Vf' * diagm(wf))

    # plotting nodes
    rp1D = LinRange(-1, 1, Nplot + 1)
    Vp1D = vandermonde(Line(), N, rp1D) / VDM_1D
    Vp = build_kronecker_product(Vp1D, Vp1D, Vp1D)
    rp, sp, tp = vec.(StartUpDG.NodesAndModes.meshgrid(rp1D, rp1D, rp1D))
    
    return RefElemData(elem, approximation_type, N, fv, V1,
                       tuple(r, s, t), VDM, vec(Fmask),
                       tuple(rp, sp, tp), Vp, 
                       tuple(rq, sq, tq), wq, Vq,
                       tuple(rf, sf, tf), wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, (Dr, Ds, Dt), LIFT)
end

RefElemData(elem::Hex, approximation_parameter::TensorProductQuadrature, N; Nplot = 10) = 
    RefElemData(elem, Polynomial(approximation_parameter), N; Nplot)

    
function tensor_product_quadrature(::Line, r1D, w1D)
    return r1D, w1D
end
  
function tensor_product_quadrature(::Union{Tri, Quad}, r1D, w1D)
    sq, rq = vec.(StartUpDG.NodesAndModes.meshgrid(r1D))
    ws, wr = vec.(StartUpDG.NodesAndModes.meshgrid(w1D))
    wq = wr .* ws
    return rq, sq, wq
end
  
function tensor_product_quadrature(::Union{Tet, Hex}, r1D, w1D)
    rq, sq, tq = vec.(StartUpDG.NodesAndModes.meshgrid(r1D, r1D, r1D))
    wr, ws, wt = vec.(StartUpDG.NodesAndModes.meshgrid(w1D, w1D, w1D))
    wq = wr .* ws .* wt
    return rq, sq, tq, wq
end

"""
    RefElemData(elem::Union{Tri, Tet, Pyr}, approx_type::Polynomial{<:TensorProductQuadrature}, N; kwargs...)
    RefElemData(elem::Union{Wedge}, 
                     approx_type::Polynomial{<:TensorProductQuadrature}, N; 
                     quad_rule_tri = stroud_quad_nodes(Tri(), 2 * N),
                     quad_rule_line = gauss_quad(0, 0, N),
                     kwargs...)

Uses collapsed coordinate volume quadrature. Should be called via
```julia
RefElemData(Tri(), Polynomial(TensorProductQuadrature()), N)
```
"""
function RefElemData(elem::Union{Tri, Tet, Pyr}, approx_type::Polynomial{<:TensorProductQuadrature}, N; kwargs...)
    rd = RefElemData(elem, Polynomial{MultidimensionalQuadrature}(), N; 
                     quad_rule_vol=stroud_quad_nodes(elem, 2 * N), kwargs...)
    @set rd.approximation_type = approx_type
    return rd
end

function RefElemData(elem::Union{Wedge}, 
                     approx_type::Polynomial{<:TensorProductQuadrature}, N; 
                     quad_rule_tri = stroud_quad_nodes(Tri(), 2 * N),
                     quad_rule_line = gauss_quad(0, 0, N),
                     kwargs...)

    rq_tri, sq_tri, wq_tri = quad_rule_tri
    rq_1D, wq_1D = quad_rule_line
    rq = repeat(rq_tri, length(rq_1D))
    sq = repeat(sq_tri, length(rq_1D))
    tq = repeat(rq_1D, length(rq_tri))
    wq = repeat(wq_tri, length(wq_1D)) .* repeat(wq_1D, length(wq_tri))
    quad_rule_vol = (rq, sq, tq, wq)

    rd = RefElemData(elem, Polynomial{MultidimensionalQuadrature}(), N; 
                     quad_rule_vol, kwargs...)
    @set rd.approximation_type = approx_type
    return rd          
end


"""
    RefElemData(elem::Union{Line, Quad, Hex}, approximation_type::Polynomial{Gauss}, N)

Builds a `rd::RefElemData` with (N+1)-point Gauss quadrature in each dimension. 
"""
function RefElemData(element_type::Union{Line, Quad, Hex}, 
                     approximation_type::Polynomial{<:TensorProductGaussCollocation}, 
                     N; kwargs...) 

    quadrature_type = TensorProductQuadrature(gauss_quad(0, 0, N))
    rd = RefElemData(element_type, Polynomial(quadrature_type), N; kwargs...)
    return @set rd.approximation_type = approximation_type 
end
  
RefElemData(element_type::Union{Line, Quad, Hex}, ::TensorProductGaussCollocation, N; kwargs...) = 
    RefElemData(element_type, Polynomial{TensorProductGaussCollocation}(), N; kwargs...)
