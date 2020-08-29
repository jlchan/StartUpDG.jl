#####
##### ordering of faces in terms of vertices
#####

function tri_face_vertices()
        return [1,2],[2,3],[3,1]
end

function quad_face_vertices()
        return [1,2],[2,4],[3,4],[1,3] # ordering matters
end

function hex_face_vertices()
        x1D = LinRange(-1,1,2)
        r, s, t = meshgrid(x1D,x1D,x1D)
        fv1 = map(x->x[1], findall(@. abs(r+1) < 1e-10))
        fv2 = map(x->x[1], findall(@. abs(r-1) < 1e-10))
        fv3 = map(x->x[1], findall(@. abs(s+1) < 1e-10))
        fv4 = map(x->x[1], findall(@. abs(s-1) < 1e-10))
        fv5 = map(x->x[1], findall(@. abs(t+1) < 1e-10))
        fv6 = map(x->x[1], findall(@. abs(t-1) < 1e-10))
        return fv1,fv2,fv3,fv4,fv5,fv6
end

#####
##### initialization of RefElemData
#####

function init_reference_interval(N;Nq=N+1)
    # initialize a new reference element data struct
    rd = RefElemData()

    # Construct matrices on reference elements
    r,_ = gauss_lobatto_quad(0,0,N)
    VDM = vandermonde_1D(N, r)
    Dr = grad_vandermonde_1D(N, r)/VDM
    @pack! rd = r,VDM,Dr

    V1 = vandermonde_1D(1,r)/vandermonde_1D(1,[-1;1])
    @pack! rd = V1

    rq,wq = gauss_quad(0,0,Nq)
    Vq = vandermonde_1D(N, rq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,wq,Vq,M,Pq

    rf = [-1.0;1.0]
    nrJ = [-1.0;1.0]
    Vf = vandermonde_1D(N,rf)/VDM
    LIFT = M\(Vf') # lift matrix
    @pack! rd = rf,nrJ,Vf,LIFT

    # plotting nodes
    rp = LinRange(-1,1,50)
    Vp = vandermonde_1D(N,rp)/VDM
    @pack! rd = rp,Vp

    return rd
end

function init_reference_tri(N;Nq=2*N)
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = tri_face_vertices() # set faces for triangle
    Nfaces = length(fv)
    @pack! rd = fv

    # Construct matrices on reference elements
    r, s = Tri.nodes_2D(N)
    VDM = Tri.vandermonde_2D(N, r, s)
    Vr, Vs = Tri.grad_vandermonde_2D(N, r, s)
    Dr = Vr/VDM
    Ds = Vs/VDM
    @pack! rd = r,s,VDM,Dr,Ds

    # low order interpolation nodes
    r1,s1 = Tri.nodes_2D(1)
    V1 = Tri.vandermonde_2D(1,r,s)/Tri.vandermonde_2D(1,r1,s1)
    @pack! rd = V1

    #Nodes on faces, and face node coordinate
    r1D, w1D = gauss_quad(0,0,N)
    Nfp = length(r1D) # number of points per face
    e = ones(Nfp) # vector of all ones
    z = zeros(Nfp) # vector of all zeros
    rf = [r1D; -r1D; -e];
    sf = [-e; r1D; -r1D];
    wf = vec(repeat(w1D,3,1));
    nrJ = [z; e; -e]
    nsJ = [-e; e; z]
    @pack! rd = rf,sf,wf,nrJ,nsJ

    rq,sq,wq = Tri.quad_nodes_2D(Nq)
    Vq = Tri.vandermonde_2D(N,rq,sq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,sq,wq,Vq,M,Pq

    Vf = Tri.vandermonde_2D(N,rf,sf)/VDM # interpolates from nodes to face nodes
    LIFT = M\(Vf'*diagm(wf)) # lift matrix used in rhs evaluation
    @pack! rd = Vf,LIFT

    # plotting nodes
    rp, sp = Tri.equi_nodes_2D(10)
    Vp = Tri.vandermonde_2D(N,rp,sp)/VDM
    @pack! rd = rp,sp,Vp

    return rd
end

# default to full quadrature nodes
# if quad_nodes_1D=tuple of (r1D,w1D) is supplied, use those nodes
function init_reference_quad(N,quad_nodes_1D = gauss_quad(0,0,N))
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = quad_face_vertices() # set faces for triangle
    Nfaces = length(fv)
    @pack! rd = fv

    # Construct matrices on reference elements
    r, s = Quad.nodes_2D(N)
    VDM = Quad.vandermonde_2D(N, r, s)
    Vr, Vs = Quad.grad_vandermonde_2D(N, r, s)
    Dr = Vr/VDM
    Ds = Vs/VDM
    @pack! rd = r,s,VDM

    # low order interpolation nodes
    r1,s1 = Quad.nodes_2D(1)
    V1 = Quad.vandermonde_2D(1,r,s)/Quad.vandermonde_2D(1,r1,s1)
    @pack! rd = V1

    #Nodes on faces, and face node coordinate
    r1D,w1D = quad_nodes_1D
    Nfp = length(r1D)
    e = ones(size(r1D))
    z = zeros(size(r1D))
    rf = [r1D; e; -r1D; -e]
    sf = [-e; r1D; e; -r1D]
    wf = vec(repeat(w1D,Nfaces,1));
    nrJ = [z; e; z; -e]
    nsJ = [-e; z; e; z]
    @pack! rd = rf,sf,wf,nrJ,nsJ

    # quadrature nodes - build from 1D nodes.
    # can also use "rq,sq,wq = Quad.quad_nodes_2D(2*N)"
    rq,sq = vec.(meshgrid(r1D))
    wr,ws = vec.(meshgrid(w1D))
    wq = wr .* ws
    Vq = Quad.vandermonde_2D(N,rq,sq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,sq,wq,Vq,M,Pq

    Vf = Quad.vandermonde_2D(N,rf,sf)/VDM # interpolates from nodes to face nodes
    LIFT = M\(Vf'*diagm(wf)) # lift matrix used in rhs evaluation
    @pack! rd = Dr,Ds,Vf,LIFT

    # plotting nodes
    rp, sp = Quad.equi_nodes_2D(15)
    Vp = Quad.vandermonde_2D(N,rp,sp)/VDM
    @pack! rd = rp,sp,Vp

    return rd
end

function init_reference_hex(N,quad_nodes_1D=gauss_quad(0,0,N))
    # initialize a new reference element data struct
    rd = RefElemData()

    fv = hex_face_vertices() # set faces for triangle
    Nfaces = length(fv)
    @pack! rd = fv

    # Construct matrices on reference elements
    r,s,t = Hex.nodes_3D(N)
    VDM = Hex.vandermonde_3D(N,r,s,t)
    Vr,Vs,Vt = Hex.grad_vandermonde_3D(N,r,s,t)
    Dr,Ds,Dt = (A->A/VDM).(Hex.grad_vandermonde_3D(N,r,s,t))
    @pack! rd = r,s,t,VDM

    # low order interpolation nodes
    r1,s1,t1 = Hex.nodes_3D(1)
    V1 = Hex.vandermonde_3D(1,r,s,t)/Hex.vandermonde_3D(1,r1,s1,t1)
    @pack! rd = V1

    #Nodes on faces, and face node coordinate
    r1D,w1D = quad_nodes_1D
    rquad,squad = vec.(meshgrid(r1D,r1D))
    wr,ws = vec.(meshgrid(w1D,w1D))
    wquad = wr.*ws
    e = ones(size(rquad))
    zz = zeros(size(rquad))
    rf = [-e; e; rquad; rquad; rquad; rquad]
    sf = [rquad; rquad; -e; e; squad; squad]
    tf = [squad; squad; squad; squad; -e; e]
    wf = vec(repeat(wquad,Nfaces,1));
    nrJ = [-e; e; zz;zz; zz;zz]
    nsJ = [zz;zz; -e; e; zz;zz]
    ntJ = [zz;zz; zz;zz; -e; e]

    @pack! rd = rf,sf,tf,wf,nrJ,nsJ,ntJ

    # quadrature nodes - build from 1D nodes.
    rq,sq,tq = vec.(meshgrid(r1D,r1D,r1D))
    wr,ws,wt = vec.(meshgrid(w1D,w1D,w1D))
    wq = wr.*ws.*wt
    Vq = Hex.vandermonde_3D(N,rq,sq,tq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))
    @pack! rd = rq,sq,tq,wq,Vq,M,Pq

    Vf = Hex.vandermonde_3D(N,rf,sf,tf)/VDM
    LIFT = M\(Vf'*diagm(wf))
    @pack! rd = Dr,Ds,Dt,Vf,LIFT

    # plotting nodes
    rp,sp,tp = Hex.equi_nodes_3D(15)
    Vp = Hex.vandermonde_3D(N,rp,sp,tp)/VDM
    @pack! rd = rp,sp,tp,Vp

    return rd
end
