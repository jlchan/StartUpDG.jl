module StartUpDGSummationByPartsExt

using StartUpDG
using StartUpDG: SBP, SummationByPartsDiagE, LobattoFaceNodes, LegendreFaceNodes, Tri, Tet,
                 gauss_quad, gauss_lobatto_quad

using SummationByParts: SummationByParts

# submodules of SummationByParts.jl which we use to generate cubature rules
const Cubature = SummationByParts.Cubature
const SymCubatures = SummationByParts.SymCubatures

_include_vertices(::Type{LobattoFaceNodes}) = true
_include_vertices(::Type{LegendreFaceNodes}) = false
function _include_vertices(::Type{F}) where {F}
    throw(ArgumentError("Unknown face node type $F for `SummationByPartsDiagE`. " *
                        "Expected `LobattoFaceNodes` or `LegendreFaceNodes`."))
end

# the face quadrature must be exact for degree 2N polynomials for the resulting 
# operators to satisfy the SBP property and be exact for degree N polynomials. 
function _check_face_quadrature_degree(face_quadrature_degree, N)
    if face_quadrature_degree < 2 * N
        throw(ArgumentError("The face quadrature rule (degree $face_quadrature_degree) " *
                            "is not accurate enough for a degree N = $N SBP operator."))
    end
end

# Triangles: `getTriCubatureDiagE` returns a cubature whose boundary nodes are 
# Gauss-Lobatto (`vertices = true`) or Gauss-Legendre (`vertices = false`) nodes on each edge. 
function StartUpDG._summation_by_parts_diagE_nodes(::Tri,
                                                   approxType::SBP{<:SummationByPartsDiagE{F}},
                                                   N) where {F}
    q, approxType = StartUpDG._resolve_quadrature_degree(approxType, N)
    vertices = _include_vertices(F)
    (; tol) = approxType.data

    cub, vtx = Cubature.getTriCubatureDiagE(q, Float64; vertices, tol)
    nodes = SymCubatures.calcnodes(cub, vtx) # 2 x num_nodes matrix
    r, s = nodes[1, :], nodes[2, :]
    w = SymCubatures.calcweights(cub)

    # determine the 1D face quadrature rule from the number of nodes on each edge
    num_face_nodes = 2 * cub.vertices + cub.midedges + 2 * cub.numedge
    if vertices
        quad_rule_face = gauss_lobatto_quad(0, 0, num_face_nodes - 1)
        face_quadrature_degree = 2 * num_face_nodes - 3
    else
        quad_rule_face = gauss_quad(0, 0, num_face_nodes - 1)
        face_quadrature_degree = 2 * num_face_nodes - 1
    end
    _check_face_quadrature_degree(face_quadrature_degree, N)

    return (r, s, w), quad_rule_face, approxType
end

# Tetrahedra: `getTetCubatureDiagE` with `faceopertype = :DiagE` returns a cubature whose 
# face nodes coincide with a triangular cubature from `getTriCubatureForTetFaceDiagE`. 
function StartUpDG._summation_by_parts_diagE_nodes(::Tet,
                                                   approxType::SBP{<:SummationByPartsDiagE{LobattoFaceNodes}},
                                                   N)
    q, approxType = StartUpDG._resolve_quadrature_degree(approxType, N)
    (; tol) = approxType.data

    cub, vtx = Cubature.getTetCubatureDiagE(q, Float64; faceopertype = :DiagE, tol)
    nodes = SymCubatures.calcnodes(cub, vtx) # 3 x num_nodes matrix
    r, s, t = nodes[1, :], nodes[2, :], nodes[3, :]
    w = SymCubatures.calcweights(cub)

    # determine the triangular face cubature rule by matching the number of face nodes
    num_face_nodes = SymCubatures.getnumfacenodes(cub)
    face_quadrature_degree = 0
    face_cub, face_vtx = nothing, nothing
    for qf in 2:2:10
        face_cub, face_vtx = Cubature.getTriCubatureForTetFaceDiagE(qf, Float64;
                                                                    faceopertype = :DiagE,
                                                                    tol)
        if face_cub.numnodes == num_face_nodes
            face_quadrature_degree = qf
            break
        end
    end
    if face_quadrature_degree == 0
        error("Could not find a triangular face cubature rule with $num_face_nodes nodes " *
              "matching the tetrahedral cubature rule of degree $q.")
    end
    _check_face_quadrature_degree(face_quadrature_degree, N)

    face_nodes = SymCubatures.calcnodes(face_cub, face_vtx) # 2 x num_face_nodes matrix
    rf, sf = face_nodes[1, :], face_nodes[2, :]
    wf = SymCubatures.calcweights(face_cub)

    return (r, s, t, w), (rf, sf, wf), approxType
end

function StartUpDG._summation_by_parts_diagE_nodes(::Tet,
                                                   ::SBP{<:SummationByPartsDiagE{LegendreFaceNodes}},
                                                   N)
    throw(ArgumentError("`SummationByPartsDiagE{LegendreFaceNodes}` is not supported on `Tet()` " *
                        "elements; use `SummationByPartsDiagE{LobattoFaceNodes}` instead."))
end

end # module
