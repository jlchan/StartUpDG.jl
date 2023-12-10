function init_face_data(elem::Tri; quad_rule_face = gauss_quad(0,0,N))
    r1D, w1D = quad_rule_face
    e = ones(size(r1D)) 
    z = zeros(size(r1D)) 
    rf, sf = map_face_nodes(elem, r1D)
    wf = vec(repeat(w1D, 3, 1));
    nrJ = [z; e; -e]
    nsJ = [-e; e; z]
    return rf,sf,wf,nrJ,nsJ
end

function init_face_data(elem::Quad; quad_rule_face=gauss_quad(0, 0, N))
    Nfaces = 4
    r1D, w1D = quad_rule_face
    e = ones(size(r1D))
    z = zeros(size(r1D))
    rf, sf = map_face_nodes(elem, r1D)
    wf = vec(repeat(w1D, Nfaces, 1)) 
    nrJ = [-e; e; z; z]
    nsJ = [z; z; -e; e]

    return rf, sf, wf, nrJ, nsJ
end

function init_face_data(elem::Hex; quad_rule_face=quad_nodes(Quad(), N))
    rquad, squad, wquad = quad_rule_face
    e = ones(size(rquad))
    zz = zeros(size(rquad))
    rf, sf, tf = map_face_nodes(elem, rquad, squad)
    Nfaces = 6
    wf = vec(repeat(wquad, Nfaces, 1));
    nrJ = [-e;  e; zz; zz; zz; zz]
    nsJ = [zz; zz; -e;  e; zz; zz]
    ntJ = [zz; zz; zz; zz; -e;  e]
    return rf, sf, tf, wf, nrJ, nsJ, ntJ
end

function init_face_data(elem::Tet; quad_rule_face=quad_nodes(Tri(), N))
    rquad, squad, wquad = quad_rule_face
    e = ones(size(rquad))
    zz = zeros(size(rquad))
    rf, sf, tf = map_face_nodes(elem, rquad, squad)
    Nfaces = 4
    wf = vec(repeat(wquad, Nfaces, 1));
    nrJ = [zz; e; -e; zz]
    nsJ = [-e; e; zz; zz]
    ntJ = [zz; e; zz; -e]
    return rf, sf, tf, wf, nrJ, nsJ, ntJ
end


"""
    RefElemData(elem::Line, approximation_type, N;
                quad_rule_vol = quad_nodes(elem, N+1))
    RefElemData(elem, approximation_type, N;
                quad_rule_vol = quad_nodes(elem, N),
                quad_rule_face = quad_nodes(face_type(elem), N))

Constructor for `RefElemData` for different element types.
"""
function RefElemData(elem::Line, 
                     approx_type::Polynomial{DefaultPolynomialType}, N; 
                     quad_rule_vol=quad_nodes(elem, N+1), 
                     Nplot=10)

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

_default_approx_type(elem, approximation_type, N) = approximation_type
_default_approx_type(::Union{Quad, Hex}, ::Polynomial, N) = 
    Polynomial(TensorProductQuadrature(gauss_quad(0, 0, N+1)))

function RefElemData(elem::Union{Tri, Quad}, 
                     approx_type::Polynomial{DefaultPolynomialType}, N;
                     kwargs...)    
    return RefElemData2D(elem, _default_approx_type(elem, approx_type, N), 
                         N; kwargs...)
end

# internal constructor used for dispatch
function RefElemData2D(elem::Union{Tri, Quad}, 
                       approx_type, N;
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


function RefElemData(elem::Union{Hex, Tet}, 
                     approx_type::Polynomial{DefaultPolynomialType}, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face=quad_nodes(face_type(elem), N),
                     Nplot=10)

    if elem isa Hex && N > 4
        @warn "Since N > 4, we suggest using `RefElemData(Hex(), TensorProductQuadrature(gauss_quad(0, 0, $N+1)), $N)`, " * 
              "which is more efficient."
    end      

    fv = face_vertices(elem) 

    # Construct matrices on reference elements
    r, s, t = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s, t)...)
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
function RefElemData(elem::Wedge, approximation_type::Polynomial, N;
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
    RefElemData(elem::Pyr, approximation_type::Polynomial, N;
                quad_rule_vol=quad_nodes(elem, N),
                quad_rule_face_quad=quad_nodes(Quad(), N), 
                quad_rule_face_tri=quad_nodes(Tri(), N), 
                quad_rule_face=(quad_rule_face_quad, quad_rule_face_tri),
                Nplot=10)

Builds operators for pyramids.
"""
function RefElemData(elem::Pyr, approximation_type::Polynomial, N;
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
    RefElemData(elem::Hex, ::TensorProductQuadrature, N; Nplot = 10)
    RefElemData(elem::Hex, approximation_type::Polynomial{<:TensorProductQuadrature}, N; Nplot = 10)

Constructor for hexahedral `RefElemData` where the quadrature is assumed to have tensor product structure.
"""
RefElemData(elem::Hex, approximation_parameter::TensorProductQuadrature, N; Nplot = 10) = 
    RefElemData(elem, Polynomial(approximation_parameter), N; Nplot)

function RefElemData(elem::Hex, 
                     approximation_type::Polynomial{<:TensorProductQuadrature}, N;
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
    VDM = kronecker(VDM_1D, VDM_1D, VDM_1D)
    invVDM = kronecker(invVDM_1D, invVDM_1D, invVDM_1D)
    invM = kronecker(invM_1D, invM_1D, invM_1D)

    M = kronecker(M1D, M1D, M1D)

    _, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> A * invVDM).((Vr, Vs, Vt))

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)

    # Nodes on faces, and face node coordinate
    quad_rule_face = tensor_product_quadrature(face_type(elem), approximation_type.data.quad_rule_1D...)
    rf, sf, tf, wf, nrJ, nsJ, ntJ = init_face_data(elem, quad_rule_face=quad_rule_face)

    # quadrature nodes - build from 1D nodes.    
    rq, sq, tq, wq = tensor_product_quadrature(elem, approximation_type.data.quad_rule_1D...)
    Vq = kronecker(Vq1D, Vq1D, Vq1D) # vandermonde(elem, N, rq, sq, tq) * invVDM
    Pq = invM * (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf, tf) * invVDM
    LIFT = invM * (Vf' * diagm(wf))

    # plotting nodes
    rp1D = LinRange(-1, 1, Nplot + 1)
    Vp1D = vandermonde(Line(), N, rp1D) / VDM_1D
    Vp = kronecker(Vp1D, Vp1D, Vp1D)
    rp, sp, tp = vec.(StartUpDG.NodesAndModes.meshgrid(rp1D, rp1D, rp1D))
    
    return RefElemData(elem, approximation_type, N, fv, V1,
                       tuple(r, s, t), VDM, vec(Fmask),
                       tuple(rp, sp, tp), Vp, 
                       tuple(rq, sq, tq), wq, Vq,
                       tuple(rf, sf, tf), wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, (Dr, Ds, Dt), LIFT)
end


function tensor_product_quadrature(::Line, r1D, w1D)
    return r1D, w1D
end
  
function tensor_product_quadrature(::Quad, r1D, w1D)
    sq, rq = vec.(StartUpDG.NodesAndModes.meshgrid(r1D))
    ws, wr = vec.(StartUpDG.NodesAndModes.meshgrid(w1D))
    wq = wr .* ws
    return rq, sq, wq
end
  
function tensor_product_quadrature(::Hex, r1D, w1D)
    rq, sq, tq = vec.(StartUpDG.NodesAndModes.meshgrid(r1D, r1D, r1D))
    wr, ws, wt = vec.(StartUpDG.NodesAndModes.meshgrid(w1D, w1D, w1D))
    wq = wr .* ws .* wt
    return rq, sq, tq, wq
end

"""
    RefElemData(elem::Union{Line, Quad, Hex}, approximation_type::Polynomial{Gauss}, N)

Builds `rd::RefElemData` with (N+1)-point Gauss quadrature in each dimension. 
"""
function RefElemData(element_type::Union{Line, Quad, Hex}, 
                     approximation_type::Polynomial{<:Gauss}, N; kwargs...) 
    # explicitly specify Gauss quadrature rule with (N+1)^d points 
    quad_rule_vol = tensor_product_quadrature(element_type, StartUpDG.gauss_quad(0, 0, N)...)

    # hacky fix to call the default constructor
    rd = RefElemData(element_type, Polynomial(), N; quad_rule_vol)    
    return @set rd.approximation_type = approximation_type 
end
  
RefElemData(element_type::Union{Line, Quad, Hex}, ::Gauss, N; kwargs...) = 
    RefElemData(element_type, Polynomial{Gauss}(), N; kwargs...)
