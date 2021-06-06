# define ApproximationType
struct Polynomial end 

# Type alias where all operators are just Matrix{Tv}
const PolynomialRefElemData{Dim,ElemShape,Nfaces,Tv} = RefElemData{Dim,ElemShape,Polynomial,Nfaces,Tv,Matrix{Tv},Matrix{Tv},Matrix{Tv},Matrix{Tv},Matrix{Tv},Matrix{Tv}}

function init_face_data(elem::Tri, N; quad_nodes_face=gauss_quad(0,0,N))
    #Nodes on faces, and face node coordinate
    r1D, w1D = quad_nodes_face
    e = ones(size(r1D)) # vector of all ones
    z = zeros(size(r1D)) # vector of all zeros
    rf,sf = map_face_nodes(elem,r1D)
    wf = vec(repeat(w1D,3,1));
    nrJ = [z; e; -e]
    nsJ = [-e; e; z]
    return rf,sf,wf,nrJ,nsJ
end

function init_face_data(elem::Quad,N; quad_nodes_face=gauss_quad(0,0,N))
    r1D,w1D = quad_nodes_face
    e = ones(size(r1D))
    z = zeros(size(r1D))
    rf,sf = map_face_nodes(elem,r1D)
    wf = vec(repeat(w1D,4,1)); # 4 faces
    nrJ = [z; e; z; -e]
    nsJ = [-e; z; e; z]
    return rf,sf,wf,nrJ,nsJ
end

function init_face_data(elem::Hex, N)
    rquad,squad,wquad = quad_nodes(Quad(),N)
    e = ones(size(rquad))
    zz = zeros(size(rquad))
    rf,sf,tf = map_face_nodes(elem,rquad,squad)
    Nfaces = 6
    wf = vec(repeat(wquad,Nfaces,1));
    nrJ = [-e; e; zz;zz; zz;zz]
    nsJ = [zz;zz; -e; e; zz;zz]
    ntJ = [zz;zz; zz;zz; -e; e]
    return rf,sf,tf,wf,nrJ,nsJ,ntJ
end


"""
    RefElemData(elem::Line, N;
                quad_rule_vol = quad_nodes(elem,N+1))
    RefElemData(elem::Union{Tri,Quad}, N;
                 quad_rule_vol = quad_nodes(elem,N),
                 quad_rule_face = gauss_quad(0,0,N))
    RefElemData(elem::Hex,N;
                 quad_rule_vol = quad_nodes(elem,N),
                 quad_rule_face = quad_nodes(Quad(),N))
    RefElemData(elem; N, kwargs...) # version with keyword arg

Constructor for RefElemData for different element types.
"""
function RefElemData(elem::Line, approxType::Polynomial, N; quad_rule_vol = quad_nodes(elem,N+1), Nplot=10)

    fv = face_vertices(elem)

    # Construct matrices on reference elements
    r = nodes(elem,N)
    Fmask = [1 N+1]
    VDM = vandermonde(elem, N, r)
    Dr = grad_vandermonde(elem, N, r)/VDM

    V1 = vandermonde(elem,1,r)/vandermonde(elem,1,[-1;1])

    rq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    rf = [-1.0;1.0]
    nrJ = [-1.0;1.0]
    wf = [1.0;1.0]
    Vf = vandermonde(elem,N,rf)/VDM
    LIFT = M\(Vf') # lift matrix

    # plotting nodes
    rp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp)/VDM

    return RefElemData(elem,approxType,N,fv,V1,
                       tuple(r),VDM,vec(Fmask),
                       Nplot,tuple(rp),Vp,
                       tuple(rq),wq,Vq,
                       tuple(rf),wf,Vf,tuple(nrJ),
                       M,Pq,tuple(Dr),LIFT)
end

function RefElemData(elem::Union{Tri,Quad},  approxType::Polynomial, N;
                     quad_rule_vol = quad_nodes(elem,N),
                     quad_rule_face = gauss_quad(0,0,N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle

    # Construct matrices on reference elements
    r,s = nodes(elem,N)
    Fmask = hcat(find_face_nodes(elem,r,s)...)

    VDM,Vr,Vs = basis(elem,N,r,s)
    Dr = Vr/VDM
    Ds = Vs/VDM

    # low order interpolation nodes
    r1,s1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s)/vandermonde(elem,1,r1,s1)

    rf,sf,wf,nrJ,nsJ = init_face_data(elem,N,quad_nodes_face=quad_rule_face)

    rq,sq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq,sq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    Vf = vandermonde(elem,N,rf,sf)/VDM # interpolates from nodes to face nodes
    LIFT = M\(Vf'*diagm(wf)) # lift matrix used in rhs evaluation

    # plotting nodes
    rp, sp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp,sp)/VDM

    Drs = (Dr,Ds)

    # sparsify for Quad
    # tol = 1e-13
    # Drs = typeof(elem)==Quad ? droptol!.(sparse.((Dr,Ds)),tol) : (Dr,Ds)
    # Vf = typeof(elem)==Quad ? droptol!(sparse(Vf),tol) : Vf
    # LIFT = typeof(elem)==Quad ? droptol!(sparse(LIFT),tol) : LIFT

    return RefElemData(elem,approxType,N,fv,V1,
                       tuple(r,s),VDM,vec(Fmask),
                       Nplot,tuple(rp,sp),Vp,
                       tuple(rq,sq),wq,Vq,
                       tuple(rf,sf),wf,Vf,tuple(nrJ,nsJ),
                       M,Pq,Drs,LIFT)
end

function RefElemData(elem::Hex, approxType::Polynomial, N;
                     quad_rule_vol = quad_nodes(elem,N),
                     quad_rule_face = quad_nodes(Quad(),N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle

    # Construct matrices on reference elements
    r,s,t = nodes(elem,N)
    Fmask = hcat(find_face_nodes(elem,r,s,t)...)
    VDM,Vr,Vs,Vt = basis(elem,N,r,s,t)
    Dr,Ds,Dt = (A->A/VDM).((Vr,Vs,Vt))

    # low order interpolation nodes
    r1,s1,t1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s,t)/vandermonde(elem,1,r1,s1,t1)

    #Nodes on faces, and face node coordinate
    rf,sf,tf,wf,nrJ,nsJ,ntJ = init_face_data(elem,N)

    # quadrature nodes - build from 1D nodes.
    rq,sq,tq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq,sq,tq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    Vf = vandermonde(elem,N,rf,sf,tf)/VDM
    LIFT = M\(Vf'*diagm(wf))

    # plotting nodes
    rp,sp,tp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp,sp,tp)/VDM

    # Drst = sparse.((Dr,Ds,Dt))
    Drst = (Dr,Ds,Dt)
    # Vf = sparse(Vf)

    return RefElemData(elem,approxType,N,fv,V1,
                       tuple(r,s,t),VDM,vec(Fmask),
                       Nplot,tuple(rp,sp,tp),Vp,
                       tuple(rq,sq,tq),wq,Vq,
                       tuple(rf,sf,tf),wf,Vf,tuple(nrJ,nsJ,ntJ),
                       M,Pq,Drst,LIFT)
end
