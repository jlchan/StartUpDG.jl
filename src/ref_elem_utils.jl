#####
##### face data for diff elements
#####

function map_face_nodes(::Tri, face_nodes)
    r1D = face_nodes
    e = ones(size(r1D)) # vector of all ones
    rf = [r1D; -r1D; -e];
    sf = [-e; r1D; -r1D];
    return rf, sf
end

function map_face_nodes(::Quad, face_nodes)
    r1D = face_nodes
    e = ones(size(r1D))
    rf = [-e; e; r1D; r1D]
    sf = [r1D; r1D; -e; e]
    return rf, sf
end

function map_face_nodes(::Hex, face_nodes...)
    r, s = face_nodes
    e = ones(size(r))
    rf = [-e; e; r; r; r; r]
    sf = [r; r; -e; e; s; s]
    tf = [s; s; s; s; -e; e]
    return rf, sf, tf
end

function map_face_nodes(::Tet, face_nodes...)
    r, s = face_nodes
    e = ones(size(r))
    rf = [r; -(e + r + s); -e; r]
    sf = [-e; r; r; s]
    tf = [s; s; s; -e]
    return rf, sf, tf
end

function init_face_data(elem::Tri; quad_rule_face = gauss_quad(0, 0, N))
    r1D, w1D = quad_rule_face
    e = ones(size(r1D)) 
    z = zeros(size(r1D)) 
    rf, sf = map_face_nodes(elem, r1D)
    wf = vec(repeat(w1D, 3, 1));
    nrJ = [z; e; -e]
    nsJ = [-e; e; z]
    return rf, sf, wf, nrJ, nsJ
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

# default to doing nothing
map_nodes_to_symmetric_element(element_type, rst...) = rst

# for triangles, map to an equilateral triangle
function map_nodes_to_symmetric_element(::Tri, r, s)
    # biunit right triangular vertices
    v1, v2, v3 = SVector{2}.(zip(nodes(Tri(), 1)...))

    denom = (v2[2] - v3[2]) * (v1[1] - v3[1]) + (v3[1] - v2[1]) * (v1[2] - v3[2])
    L1 = @. ((v2[2] - v3[2]) * (r - v3[1]) + (v3[1] - v2[1]) * (s - v3[2])) / denom
    L2 = @. ((v3[2] - v1[2]) * (r - v3[1]) + (v1[1] - v3[1]) * (s - v3[2])) / denom
    L3 = @. 1 - L1 - L2

    # equilateral vertices
    v1 = SVector{2}(2 * [-.5, -sqrt(3) / 6])
    v2 = SVector{2}(2 * [.5, -sqrt(3)/6])
    v3 = SVector{2}(2 * [0, sqrt(3)/3])

    x = @. v1[1] * L1 + v2[1] * L2 + v3[1] * L3
    y = @. v1[2] * L1 + v2[2] * L2 + v3[2] * L3

    return x, y
end


"""
    function inverse_trace_constant(rd::RefElemData)

Returns the degree-dependent constant in the inverse trace equality over the reference element (as
reported in ["GPU-accelerated dG methods on hybrid meshes"](https://doi.org/10.1016/j.jcp.2016.04.003)
by Chan, Wang, Modave, Remacle, Warburton 2016).

Can be used to estimate dependence of maximum stable timestep on degree of approximation.
"""
inverse_trace_constant(rd::RefElemData{1, Line, <:Polynomial}) = (rd.N+1) * (rd.N+2) / 2
inverse_trace_constant(rd::RefElemData{2, Quad, <:Polynomial}) = (rd.N+1) * (rd.N + 2)
inverse_trace_constant(rd::RefElemData{3, Hex, <:Polynomial}) = 3 * (rd.N + 1) * (rd.N + 2) / 2
inverse_trace_constant(rd::RefElemData{1, Line, SBP{TensorProductLobatto}}) = rd.N * (rd.N + 1) / 2
inverse_trace_constant(rd::RefElemData{2, Quad, SBP{TensorProductLobatto}}) = rd.N * (rd.N + 1)
inverse_trace_constant(rd::RefElemData{3, Hex, SBP{TensorProductLobatto}}) = 3 * rd.N * (rd.N + 1) / 2

# precomputed
_inverse_trace_constants(rd::RefElemData{2, Tri, <:Polynomial}) = (6.0, 10.898979485566365, 16.292060161853993, 23.999999999999808, 31.884512140579055, 42.42373503225737, 52.88579066878113, 66.25284319164409, 79.3535377715693, 95.53911875636945)
_inverse_trace_constants(rd::RefElemData{3, Tet, <:Polynomial}) = (10., 16.892024376045097, 23.58210016200093, 33.828424659883034, 43.40423356477473, 56.98869932201791, 69.68035962892684)
_inverse_trace_constants(rd::RefElemData{3, <:Wedge, <:Polynomial}) = (9.92613593327531, 18.56357670538197, 29.030325215439625, 42.98834597283998, 58.802145509223536, 78.00615833786019, 99.27149051377008, 123.76230676427465, 150.48304574455943)
_inverse_trace_constants(rd::RefElemData{3, <:Pyr, <:Polynomial}) = (10.000000000000007, 17.524350232967805, 27.133709492299054, 38.83465837771858, 52.6289353625632, 68.51745536932947, 86.50095449182132, 106.58014715058823, 128.75578358915527)
_inverse_trace_constants(rd::RefElemData{2, Tri, SBP{Kubatko{LegendreFaceNodes}}}) = (3.0, 6.6666666666667105, 11.519123865805936, 17.945430732284063, 25.73017849975611, 47.67304311231776)
_inverse_trace_constants(rd::RefElemData{2, Tri, SBP{Kubatko{LobattoFaceNodes}}}) = (4.0, 10.0, 13.606721028333366, 17.948880228130577, 41.14350198287034, 340.3588047925354)
_inverse_trace_constants(rd::RefElemData{2, Tri, SBP{Hicken}}) = (6.666666666666666, 13.309638971217003, 21.905170262973158, 30.569992349262947)

inverse_trace_constant(rd::RefElemData{2, Tri, <:Polynomial})     = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{2, Tri, <:SBP})            = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{3, Tet, <:Polynomial})     = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{3, <:Wedge, <:Polynomial}) = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{3, <:Pyr, <:Polynomial})   = _inverse_trace_constants(rd)[rd.N]

# generic fallback
function inverse_trace_constant(rd::RefElemData) 
    @warn "Computing the inverse trace constant using an eigenvalue problem; this may be expensive."
    return eigenvalue_inverse_trace_constant(rd)
end

function eigenvalue_inverse_trace_constant(rd::RefElemData)
    return maximum(real(eigvals(Matrix(rd.Vf' * diagm(rd.wf) * rd.Vf), Matrix(rd.Vq' * diagm(rd.wq) * rd.Vq))))
end
