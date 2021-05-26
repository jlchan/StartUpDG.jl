#####
##### ordering of faces in terms of vertices
#####
function face_vertices(elem::Line)
    return 1,2
end

function face_vertices(elem::Tri)
        tol = 5e2*eps()
        r,s = nodes(Tri(),1)
        e1 = findall(@. abs(s+1)<tol)
        e2 = findall(@. abs(r+s)<tol)
        e3 = findall(@. abs(r+1)<tol)
        return e1,e2,e3
end

function face_vertices(elem::Quad)
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
