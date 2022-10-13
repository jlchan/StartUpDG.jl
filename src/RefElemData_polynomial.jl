function init_face_data(elem::Tri; quad_rule_face = gauss_quad(0,0,N))
    r1D, w1D = quad_rule_face
    e = ones(size(r1D)) 
    z = zeros(size(r1D)) 
    rf,sf = map_face_nodes(elem,r1D)
    wf = vec(repeat(w1D,3,1));
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
    RefElemData(elem::Line, N;
                quad_rule_vol = quad_nodes(elem, N+1))
    RefElemData(elem::Union{Tri, Quad}, N;
                 quad_rule_vol = quad_nodes(elem, N),
                 quad_rule_face = gauss_quad(0, 0, N))
    RefElemData(elem::Union{Hex, Tet}, N;
                 quad_rule_vol = quad_nodes(elem, N),
                 quad_rule_face = quad_nodes(Quad(), N))
    RefElemData(elem; N, kwargs...) # version with keyword args

Constructor for `RefElemData` for different element types.
"""
function RefElemData(elem::Line, approxType::Polynomial, N; 
                     quad_rule_vol=quad_nodes(elem, N+1), 
                     Nplot=10)

    fv = face_vertices(elem)

    # Construct matrices on reference elements
    r = nodes(elem, N)
    Fmask = [1 N+1]
    VDM = vandermonde(elem, N, r)
    Dr = grad_vandermonde(elem, N, r)/VDM

    V1 = vandermonde(elem, 1, r)/vandermonde(elem, 1, [-1; 1])

    rq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq) / VDM
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    rf  = [-1.0; 1.0]
    nrJ = [-1.0; 1.0]
    wf  = [1.0; 1.0]
    Vf = vandermonde(elem, N, rf) / VDM
    LIFT = M \ (Vf') # lift matrix

    # plotting nodes
    rp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp) / VDM

    return RefElemData(elem, approxType, N, fv, V1,
                       tuple(r), VDM, vec(Fmask),
                       Nplot, tuple(rp), Vp,
                       tuple(rq), wq, Vq,
                       tuple(rf), wf, Vf, tuple(nrJ),
                       M, Pq, tuple(Dr), LIFT)
end

function RefElemData(elem::Union{Tri, Quad},  approxType::Polynomial, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face=quad_nodes(face_type(elem), N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle

    # Construct matrices on reference elements
    r,s = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s)...)

    VDM, Vr, Vs = basis(elem, N, r, s)
    Dr = Vr / VDM
    Ds = Vs / VDM

    # low order interpolation nodes
    r1, s1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s) / vandermonde(elem, 1, r1, s1)

    rf, sf, wf, nrJ, nsJ = init_face_data(elem, quad_rule_face = quad_rule_face)
    rq, sq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq, sq) / VDM
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf) / VDM # interpolates from nodes to face nodes
    LIFT = M \ (Vf' * diagm(wf)) # lift matrix used in rhs evaluation

    # plotting nodes
    rp, sp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp) / VDM

    return RefElemData(elem, approxType, N, fv, V1,
                       tuple(r, s), VDM, vec(Fmask),
                       Nplot, tuple(rp, sp), Vp,
                       tuple(rq, sq), wq, Vq,
                       tuple(rf, sf), wf, Vf, tuple(nrJ, nsJ),
                       M, Pq, (Dr, Ds), LIFT)
end

function RefElemData(elem::Tet, approxType::Polynomial, N;
                     quad_rule_vol=quad_nodes(elem, N),
                     quad_rule_face=quad_nodes(face_type(elem), N),
                     Nplot=10)

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
    rf, sf, tf, wf, nrJ, nsJ, ntJ = init_face_data(elem, quad_rule_face = quad_rule_face)
    # quadrature nodes - build from 1D nodes.
    rq, sq, tq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq, sq, tq) / VDM
    M = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf, tf) / VDM
    LIFT = M \ (Vf' * diagm(wf))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp, tp) / VDM

    return RefElemData(elem, approxType, N, fv, V1,
                       tuple(r, s, t), VDM, vec(Fmask),
                       Nplot, tuple(rp, sp, tp), Vp,
                       tuple(rq, sq, tq), wq, Vq,
                       tuple(rf, sf, tf), wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, (Dr, Ds, Dt), LIFT)
end

# specialize constructor for `Hex` to allow for higher polynomial degrees `N`
function RefElemData(elem::Hex, approxType::Polynomial, N;
                     quad_rule_vol = quad_nodes(elem, N),
                     quad_rule_face = quad_nodes(face_type(elem), N),
                     Nplot = 10)

    fv = face_vertices(elem) 

    # Construct matrices on reference elements
    r, s, t = nodes(elem, N)
    Fmask = hcat(find_face_nodes(elem, r, s, t)...)

    # construct 1D operator for faster Kronecker solves
    r1D = nodes(Line(), N)
    rq1D, wq1D = quad_nodes(Line(), N)
    VDM_1D = vandermonde(Line(), N, r1D)
    Vq1D = vandermonde(Line(), N, rq1D) / VDM_1D
    invVDM_1D = inv(VDM_1D)
    invM_1D = VDM_1D * VDM_1D'
    M1D = Vq1D' * diagm(wq1D) * Vq1D

    # form kronecker products of multidimensional matrices to invert/multiply
    M = kronecker(M1D, M1D, M1D)
    VDM = kronecker(VDM_1D, VDM_1D, VDM_1D)
    invVDM = kronecker(invVDM_1D, invVDM_1D, invVDM_1D)
    invM = kronecker(invM_1D, invM_1D, invM_1D)

    _, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> A * invVDM).((Vr, Vs, Vt))

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)

    # Nodes on faces, and face node coordinate
    rf, sf, tf, wf, nrJ, nsJ, ntJ = init_face_data(elem, quad_rule_face=quad_rule_face)

    # quadrature nodes - build from 1D nodes.
    rq, sq, tq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq, sq, tq) * invVDM
    # M = Vq' * diagm(wq) * Vq
    Pq = invM * (Vq' * diagm(wq))

    Vf = vandermonde(elem, N, rf, sf, tf) * invVDM
    LIFT = invM * (Vf' * diagm(wf))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    # Vp = vandermonde(elem, N, rp, sp, tp) * invVDM 
    rp1D = LinRange(-1, 1, Nplot + 1)
    Vp1D = vandermonde(Line(), N, rp1D) / VDM_1D
    Vp = kronecker(Vp1D, Vp1D, Vp1D)

    Drst = (Dr, Ds, Dt)

    return RefElemData(elem, approxType, N, fv, V1,
                       tuple(r, s, t), VDM, vec(Fmask),
                       Nplot, tuple(rp, sp, tp), Vp, 
                       tuple(rq, sq, tq), wq, Vq,
                       tuple(rf, sf, tf), wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, Drst, LIFT)
end


"""
    RefElemData(elem::Wedge, approxType::Polynomial, N;
    quad_rule_vol=quad_nodes(elem, N),
    quad_rule_face=quad_nodes(face_type(elem), N),
    Nplot=10)

Build operators for prisms/wedges
"""
function RefElemData(elem::Wedge, approxType::Polynomial, N;
    quad_rule_vol=quad_nodes(elem, N),
    quad_rule_face=quad_nodes(face_type(elem, 1), N),
    Nplot=10)

    #Find the vertices of the faces
    fv = face_vertices(elem)
    #Get interpolation nodes of degree N 
    r, s, t = nodes(elem, N)
    
    VDM, Vr, Vs, Vt = basis(elem, N, r, s, t)
    Dr, Ds, Dt = (A -> A/VDM).((Vr, Vs, Vt))
    Drst = (Dr, Ds, Dt)

    fn1, fn2, fn3, fn4, fn5  = find_face_nodes(elem, r, s, t)
    Fmask = tuple(fn1, fn2, fn3, fn4, fn5)

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)
    
    rqf, sqf, wqf = quad_nodes(Quad(), N)
    rtf, stf, wtf = quad_nodes(Tri(), N)
    tri = RefElemData(Tri(), N)


    otf = ones(size(rtf))
    #sqf refers to the height of quad-face nodes. 
    rf = vcat(repeat(tri.rf[1:(N+1)],N+1),repeat(tri.rf[N+2:2*(N+1)],N+1), repeat(tri.rf[(2*(N+1)+1):3*(N+1)],N+1), rtf, rtf)
    sf = vcat(repeat(tri.sf[1:(N+1)],N+1),repeat(tri.sf[N+2:2*(N+1)],N+1), repeat(tri.sf[(2*(N+1)+1):3*(N+1)],N+1), stf, stf)
    tf = vcat(sqf, sqf, sqf, -otf, otf)
    wf = vec(append!(repeat(wqf, 3), repeat(wtf, 2)))

    rstf = tuple(rf, sf, tf)
    Vf = vandermonde(elem, N, rf, sf, tf) / VDM
    
    #quad face quadrature
    zt = zeros(size(fn4)) 
    et = ones(size(fn4))
    zq = zeros(size(fn1))
    eq = ones(size(fn1))

    #for nrJ and nsJ normal on face 1-3 coincide with the triangular normals
    nrJ = [zq; eq; -eq; zt; zt]
    nsJ = [-eq; eq; zq; zt; zt]
    ntJ = [zq; zq; zq; -et; et]     #normal on the face??


    rq, sq, tq, wq = quad_rule_vol
    Vq = vandermonde(elem, N, rq, sq, tq) / VDM
    M  = Vq' * diagm(wq) * Vq
    Pq = M \ (Vq' * diagm(wq))

    # plotting nodes
    rp, sp, tp = equi_nodes(elem, Nplot)
    Vp = vandermonde(elem, N, rp, sp, tp) / VDM

    LIFT = M \ (Vf' *diagm(wf))

    return RefElemData(elem, approxType, N, fv, V1,
                        tuple(r, s, t), VDM, Fmask,
                        Nplot, tuple(rp, sp, tp), Vp,
                        tuple(rq, sq, tq), wq, Vq,
                        rstf, wf, Vf, tuple(nrJ, nsJ, ntJ),
                        M, Pq, Drst, LIFT)
end

