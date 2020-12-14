#####
##### ordering of faces in terms of vertices
#####
function face_vertices(elem::Line)
    return 1,2
end

function face_vertices(elem::Tri)
        #return [1,2],[2,3],[3,1]
        # return [1,3],[3,2],[2,1]
        tol = 5e2*eps()
        r,s = nodes(Tri(),1)
        e1 = findall(@. abs(s+1)<tol)
        e2 = findall(@. abs(r+s)<tol)
        e3 = findall(@. abs(r+1)<tol)
        return e1,e2,e3
end

function face_vertices(elem::Quad)
        #return [1,2],[2,4],[3,4],[1,3] # ordering matters
        tol = 5e2*eps()
        r,s = nodes(Quad(),1)
        e1 = findall(@. abs(s+1)<tol)
        e2 = findall(@. abs(r-1)<tol)
        e3 = findall(@. abs(s-1)<tol)
        e4 = findall(@. abs(r+1)<tol)
        return e1,e2,e3,e4
end

function face_vertices(elem::Hex)
        x1D = LinRange(-1,1,2)
        r, s, t = vec.(meshgrid(x1D,x1D,x1D))
        fv1 = map(x->x[1], findall(@. abs(r+1) < 1e-10))
        fv2 = map(x->x[1], findall(@. abs(r-1) < 1e-10))
        fv3 = map(x->x[1], findall(@. abs(s+1) < 1e-10))
        fv4 = map(x->x[1], findall(@. abs(s-1) < 1e-10))
        fv5 = map(x->x[1], findall(@. abs(t+1) < 1e-10))
        fv6 = map(x->x[1], findall(@. abs(t-1) < 1e-10))
        return fv1,fv2,fv3,fv4,fv5,fv6
        # return fv3,fv4,fv1,fv2,fv5,fv6
end

#####
##### face data for diff elements
#####

function init_face_data(elem::Tri, N; quad_nodes_face=gauss_quad(0,0,N))
    #Nodes on faces, and face node coordinate
    r1D, w1D = quad_nodes_face
    e = ones(size(r1D)) # vector of all ones
    z = zeros(size(r1D)) # vector of all zeros
    rf = [r1D; -r1D; -e];
    sf = [-e; r1D; -r1D];
    wf = vec(repeat(w1D,3,1));
    nrJ = [z; e; -e]
    nsJ = [-e; e; z]
    return rf,sf,wf,nrJ,nsJ
end

function init_face_data(elem::Quad,N; quad_nodes_face=gauss_quad(0,0,N))
    r1D,w1D = quad_nodes_face
    e = ones(size(r1D))
    z = zeros(size(r1D))
    rf = [r1D; e; -r1D; -e]
    sf = [-e; r1D; e; -r1D]
    wf = vec(repeat(w1D,4,1)); # 4 faces
    nrJ = [z; e; z; -e]
    nsJ = [-e; z; e; z]
    return rf,sf,wf,nrJ,nsJ
end

function init_face_data(elem::Hex, N)
    rquad,squad,wquad = quad_nodes(Quad(),N)
    e = ones(size(rquad))
    zz = zeros(size(rquad))
    rf = [-e; e; rquad; rquad; rquad; rquad]
    sf = [rquad; rquad; -e; e; squad; squad]
    tf = [squad; squad; squad; squad; -e; e]
    Nfaces = 6
    wf = vec(repeat(wquad,Nfaces,1));
    nrJ = [-e; e; zz;zz; zz;zz]
    nsJ = [zz;zz; -e; e; zz;zz]
    ntJ = [zz;zz; zz;zz; -e; e]
    return rf,sf,tf,wf,nrJ,nsJ,ntJ
end

#####
##### initialization of RefElemData
#####

function init_reference_elem(elem::Line,N;Nq=N+1)
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = face_vertices(elem)
    Nfaces = length(fv)
    @pack! rd = Nfaces,fv

    # Construct matrices on reference elements
    r = nodes(elem,N)
    VDM = vandermonde(elem, N, r)
    Dr = grad_vandermonde(elem, N, r)/VDM
    @pack! rd = r,VDM,Dr

    V1 = vandermonde(elem,1,r)/vandermonde(elem,1,[-1;1])
    @pack! rd = V1

    rq,wq = quad_nodes(elem,Nq)
    Vq = vandermonde(elem,N,rq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,wq,Vq,M,Pq

    rf = [-1.0;1.0]
    nrJ = [-1.0;1.0]
    wf = [1.0;1.0]
    Vf = vandermonde(elem,N,rf)/VDM
    LIFT = M\(Vf') # lift matrix
    @pack! rd = rf,wf,nrJ,Vf,LIFT

    # plotting nodes
    rp = equi_nodes(elem,10)
    Vp = vandermonde(elem,N,rp)/VDM
    @pack! rd = rp,Vp

    return rd
end

function init_reference_elem(elem::Union{Tri,Quad}, N; Nq=N,
                             quad_nodes_face=gauss_quad(0,0,N))
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = face_vertices(elem) # set faces for triangle
    Nfaces = length(fv)
    @pack! rd = fv, Nfaces

    # Construct matrices on reference elements
    r, s = nodes(elem,N)
    VDM,Vr,Vs = basis(elem,N, r, s)
    Dr = Vr/VDM
    Ds = Vs/VDM
    @pack! rd = r,s,VDM,Dr,Ds

    # low order interpolation nodes
    r1,s1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s)/vandermonde(elem,1,r1,s1)
    @pack! rd = V1

    rf,sf,wf,nrJ,nsJ = init_face_data(elem,N,quad_nodes_face=quad_nodes_face)
    @pack! rd = rf,sf,wf,nrJ,nsJ

    rq,sq,wq = quad_nodes(elem,Nq)
    Vq = vandermonde(elem,N,rq,sq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,sq,wq,Vq,M,Pq

    Vf = vandermonde(elem,N,rf,sf)/VDM # interpolates from nodes to face nodes
    LIFT = M\(Vf'*diagm(wf)) # lift matrix used in rhs evaluation
    @pack! rd = Vf,LIFT

    # plotting nodes
    rp, sp = equi_nodes(elem,10)
    Vp = vandermonde(elem,N,rp,sp)/VDM
    @pack! rd = rp,sp,Vp

    return rd
end

function init_reference_elem(elem::Hex, N; Nq = N)
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = face_vertices(elem) # set faces for triangle
    Nfaces = length(fv)
    @pack! rd = fv, Nfaces

    # Construct matrices on reference elements
    r,s,t = nodes(elem,N)
    VDM,Vr,Vs,Vt = basis(elem,N,r,s,t)
    Dr,Ds,Dt = (A->A/VDM).((Vr,Vs,Vt))
    @pack! rd = r,s,t,VDM

    # low order interpolation nodes
    r1,s1,t1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s,t)/vandermonde(elem,1,r1,s1,t1)
    @pack! rd = V1

    #Nodes on faces, and face node coordinate
    rf,sf,tf,wf,nrJ,nsJ,ntJ = init_face_data(elem,N)
    @pack! rd = rf,sf,tf,wf,nrJ,nsJ,ntJ

    # quadrature nodes - build from 1D nodes.
    rq,sq,tq,wq = quad_nodes(elem,Nq)
    Vq = vandermonde(elem,N,rq,sq,tq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,sq,tq,wq,Vq,M,Pq

    Vf = vandermonde(elem,N,rf,sf,tf)/VDM
    LIFT = M\(Vf'*diagm(wf))
    @pack! rd = Dr,Ds,Dt,Vf,LIFT

    # plotting nodes
    rp,sp,tp = equi_nodes(elem,15)
    Vp = vandermonde(elem,N,rp,sp,tp)/VDM
    @pack! rd = rp,sp,tp,Vp

    return rd
end
