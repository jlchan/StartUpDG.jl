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

#####
##### face data for diff elements
#####

function map_face_nodes(elem::Tri, face_nodes)
    r1D = face_nodes
    e = ones(size(r1D)) # vector of all ones
    rf = [r1D; -r1D; -e];
    sf = [-e; r1D; -r1D];
    return rf,sf
end
function map_face_nodes(elem::Quad, face_nodes)
    r1D = face_nodes
    e = ones(size(r1D))
    rf = [r1D; e; -r1D; -e]
    sf = [-e; r1D; e; -r1D]
    return rf,sf
end
function map_face_nodes(elem::Hex, face_nodes...)
    rquad,squad = face_nodes
    e = ones(size(rquad))
    rf = [-e; e; rquad; rquad; rquad; rquad]
    sf = [rquad; rquad; -e; e; squad; squad]
    tf = [squad; squad; squad; squad; -e; e]
    return rf,sf,tf
end

