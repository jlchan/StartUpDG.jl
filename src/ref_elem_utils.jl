#####
##### ordering of faces in terms of vertices
#####
function face_vertices(elem::Line)
    return 1,2
end

function find_face_nodes(elem::Tri,r,s,tol=50*eps())
    e1 = findall(@. abs(s+1)<tol)
    e2 = findall(@. abs(r+s)<tol)
    e3 = findall(@. abs(r+1)<tol)
    return e1,reverse(e2),reverse(e3)
end

function find_face_nodes(elem::Quad,r,s,tol=50*eps())
    e1 = findall(@. abs(s+1)<tol)
    e2 = findall(@. abs(r-1)<tol)
    e3 = findall(@. abs(s-1)<tol)
    e4 = findall(@. abs(r+1)<tol)
    return e1,e2,reverse(e3),reverse(e4)
end

function find_face_nodes(elem::Hex,r,s,t,tol=50*eps())
    fv1 = findall(@. abs(r+1) < tol)
    fv2 = findall(@. abs(r-1) < tol)
    fv3 = findall(@. abs(s+1) < tol)
    fv4 = findall(@. abs(s-1) < tol)
    fv5 = findall(@. abs(t+1) < tol)
    fv6 = findall(@. abs(t-1) < tol)
    return fv1,fv2,fv3,fv4,fv5,fv6
end

# face vertices = face nodes of degree 1
face_vertices(elem) = find_face_nodes(elem,nodes(elem,1)...)

# function face_vertices(elem::Hex)
#     x1D = LinRange(-1,1,2)
#     r, s, t = vec.(meshgrid(x1D,x1D,x1D))
#     return find_face_nodes(elem,r,s,t)
#     # return fv3,fv4,fv1,fv2,fv5,fv6
# end

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
