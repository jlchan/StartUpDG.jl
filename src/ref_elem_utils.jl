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



"""
    function inverse_trace_constant(rd::RefElemData)

Returns inverse trace constant as reported in ["GPU-accelerated dG methods on hybrid meshes"](https://doi.org/10.1016/j.jcp.2016.04.003)
by Chan, Wang, Modave, Remacle, Warburton 2016. 
"""
inverse_trace_constant(rd::RefElemData{1}) = (rd.N+1)*(rd.N+2)/2
inverse_trace_constant(rd::RefElemData{1,Line,SBP}) = rd.N*(rd.N+1)/2 # assumes SBP <=> DGSEM
inverse_trace_constant(rd::RefElemData{2,Quad}) = (rd.N+1)*(rd.N+2)
inverse_trace_constant(rd::RefElemData{2,Quad,SBP}) = rd.N*(rd.N+1) # assumes SBP <=> DGSEM
inverse_trace_constant(rd::RefElemData{3,Hex}) = 3*(rd.N+1)*(rd.N+2)/2
inverse_trace_constant(rd::RefElemData{3,Hex,SBP}) = 3*rd.N*(rd.N+1)/2 # assumes SBP <=> DGSEM

# precomputed
_inverse_trace_constants(rd::RefElemData{2,Tri,Polynomial}) = (5.999999999999999, 10.898979485566365, 16.292060161853993, 23.999999999999808, 31.884512140579055, 42.42373503225737, 52.88579066878113, 66.25284319164409, 79.3535377715693, 95.53911875636945)
inverse_trace_constant(rd::RefElemData{2,Tri,Polynomial}) = _inverse_trace_constants(rd)[rd.N]

# generic fallback
function inverse_trace_constant(rd::RefElemData)
    return maximum(eigvals(Matrix(rd.Vf'*diagm(rd.wf)*rd.Vf),Matrix(rd.M)))
end
