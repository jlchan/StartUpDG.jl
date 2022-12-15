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
    rf = [r; r; -e; r]
    sf = [-e; s; r; s]
    tf = [s; -(e + r + s); s; -e]
    return rf, sf, tf
end

"""
    function inverse_trace_constant(rd::RefElemData)

Returns the degree-dependent constant in the inverse trace equality over the reference element (as
reported in ["GPU-accelerated dG methods on hybrid meshes"](https://doi.org/10.1016/j.jcp.2016.04.003)
by Chan, Wang, Modave, Remacle, Warburton 2016).

Can be used to estimate dependence of maximum stable timestep on degree of approximation.
"""
inverse_trace_constant(rd::RefElemData{1, Line, Polynomial}) = (rd.N+1) * (rd.N+2) / 2
inverse_trace_constant(rd::RefElemData{2, Quad, Polynomial}) = (rd.N+1) * (rd.N + 2)
inverse_trace_constant(rd::RefElemData{3, Hex, Polynomial}) = 3 * (rd.N + 1) * (rd.N + 2) / 2
inverse_trace_constant(rd::RefElemData{1, Line, SBP{TensorProductLobatto}}) = rd.N * (rd.N + 1) / 2
inverse_trace_constant(rd::RefElemData{2, Quad, SBP{TensorProductLobatto}}) = rd.N * (rd.N + 1)
inverse_trace_constant(rd::RefElemData{3, Hex, SBP{TensorProductLobatto}}) = 3 * rd.N * (rd.N + 1) / 2

# precomputed
_inverse_trace_constants(rd::RefElemData{2, Tri, Polynomial})   = (6.0, 10.898979485566365, 16.292060161853993, 23.999999999999808, 31.884512140579055, 42.42373503225737, 52.88579066878113, 66.25284319164409, 79.3535377715693, 95.53911875636945)
_inverse_trace_constants(rd::RefElemData{3, Tet, Polynomial})   = (10., 16.892024376045097, 23.58210016200093, 33.828424659883034, 43.40423356477473, 56.98869932201791, 69.68035962892684)
_inverse_trace_constants(rd::RefElemData{3, <:Wedge, Polynomial}) = (9.92613593327531, 18.56357670538197, 29.030325215439625, 42.98834597283998, 58.802145509223536, 78.00615833786019, 99.27149051377008, 123.76230676427465, 150.48304574455943)

inverse_trace_constant(rd::RefElemData{2, Tri, Polynomial})     = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{3, Tet, Polynomial})     = _inverse_trace_constants(rd)[rd.N]
inverse_trace_constant(rd::RefElemData{3, <:Wedge, Polynomial}) = _inverse_trace_constants(rd)[rd.N]

# generic fallback
function inverse_trace_constant(rd::RefElemData)
    @warn "Computing the inverse trace constant using an eigenvalue problem; this may be expensive."
    return maximum(eigvals(Matrix(rd.Vf' * diagm(rd.wf) * rd.Vf), Matrix(rd.Vq' * diagm(rd.wq) * rd.Vq)))
end
