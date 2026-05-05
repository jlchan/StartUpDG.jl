module StartUpDGT8codeExt

import T8code
import MPI
using LinearAlgebra: Diagonal, norm
using NodesAndModes: Tri
using StaticArrays: SMatrix, SVector
using StartUpDG
using StartUpDG: MeshData, RefElemData,
                 connect_mesh, build_node_maps, geometric_factors, compute_normals
using T8code: ForestWrapper
using T8code.Libt8:
                    SC_LP_ESSENTIAL, SC_LP_PRODUCTION,
                    T8_ECLASS_TRIANGLE,
                    t8_init, sc_init, sc_is_initialized,
                    t8_cmesh_new_hypercube, t8_cmesh_commit,
                    t8_scheme_new_default,
                    t8_forest_new_uniform, t8_forest_is_committed,
                    t8_forest_get_local_num_leaf_elements, t8_forest_get_num_local_trees,
                    t8_forest_get_tree_num_leaf_elements,
                    t8_forest_get_leaf_element_in_tree,
                    t8_forest_element_coordinate, t8_forest_get_scheme,
                    t8_forest_get_tree_class,
                    t8_element_get_num_faces, t8_element_get_face_corner,
                    t8_forest_leaf_face_neighbors,
                    sc_free, t8_get_package_id,
                    t8_forest_t

# Exported functions and types
export T8codeForestMeshType

# --- one-time SC/t8 init (single-process v1) ---
const _sc_t8_ready = Ref(false)

function __init__()
    println("Initializing StartUpDGT8codeExt...")
    if !MPI.Initialized()
        MPI.Init()
    end
    if !_sc_t8_ready[]
        if sc_is_initialized() == 0
            sc_init(MPI.COMM_WORLD.val, 0, 1, C_NULL, SC_LP_ESSENTIAL)
        end
        t8_init(SC_LP_PRODUCTION)
        _sc_t8_ready[] = true
    end

    # Clean up t8code before MPI shuts down.
    MPI.add_finalize_hook!() do
        T8code.clean_up()
        status = T8code.Libt8.sc_finalize_noabort()
        if status != 0
            @warn("Inconsistent state detected after finalizing t8code.")
        end
    end
end

mutable struct T8codeForestMeshType{F <: ForestWrapper}
    forest_wrapper::F
    function T8codeForestMeshType(forest_wrapper::F) where {F <: ForestWrapper}
        mesh_type = new{F}(forest_wrapper)

        finalizer(mesh_type) do mesh_type
            finalize(mesh_type.forest_wrapper)
        end

        return mesh_type
    end
end

"""
    t8code_uniform_mesh(element_type, level::Integer; comm = MPI.COMM_WORLD)
"""
function t8code_uniform_mesh(element_type, level::Integer;
                             comm::MPI.Comm = MPI.COMM_WORLD)
    mpicom = comm.val

    # Build the uniform forest, it is automatically partitioned among the processes.
    cmesh = t8_cmesh_new_hypercube(T8code_Eclass(element_type), mpicom, 0, 0, 0)
    scheme = t8_scheme_new_default()
    forest = t8_forest_new_uniform(cmesh, scheme, level, 0, mpicom)

    return T8codeForestMeshType(ForestWrapper(forest))
end

T8code_Eclass(::Tri) = T8_ECLASS_TRIANGLE
T8code_Eclass(::Quad) = T8_ECLASS_QUAD
T8code_Eclass(::Hex) = T8_ECLASS_HEX
T8code_Eclass(::Tet) = T8_ECLASS_TET
T8code_Eclass(::Pyr) = T8_ECLASS_PYR
T8code_Eclass(::Wedge) = T8_ECLASS_WEDGE
T8code_Eclass(::Line) = T8_ECLASS_LINE

function _merge_vertices(points::Vector{NTuple{2, Float64}};
                         tol::Float64 = 100 * eps(Float64))
    verts = NTuple{2, Float64}[]
    id = zeros(Int, length(points))
    for (i, p) in enumerate(points)
        found = 0
        for (j, v) in enumerate(verts)
            if norm(SVector(v[1] - p[1], v[2] - p[2])) <
               tol * max(1.0, max(norm(SVector(v[1], v[2])), norm(SVector(p[1], p[2]))))
                found = j
                break
            end
        end
        if found == 0
            push!(verts, p)
            id[i] = length(verts)
        else
            id[i] = found
        end
    end
    return verts, id
end

function _collect_triangle_mesh_arrays(forest::t8_forest_t)
    @assert t8_forest_is_committed(forest) == 1
    K = Int(t8_forest_get_local_num_leaf_elements(forest))
    ntrees = Int(t8_forest_get_num_local_trees(forest))
    buf = zeros(3)
    pts = NTuple{2, Float64}[]
    elem_tree_elem = Vector{Tuple{Int, Int}}(undef, K)
    k = 0
    for itree in 0:(ntrees - 1)
        ne = Int(t8_forest_get_tree_num_leaf_elements(forest, itree))
        for ielement in 0:(ne - 1)
            k += 1
            elem_tree_elem[k] = (itree, ielement)
            elem = t8_forest_get_leaf_element_in_tree(forest, itree, ielement)
            for v in 0:2
                t8_forest_element_coordinate(forest, itree, elem, v, pointer(buf))
                push!(pts, (buf[1], buf[2]))
            end
        end
    end
    verts, pid = _merge_vertices(pts)
    EToV = Matrix{Int}(undef, K, 3)
    for e in 1:K
        o = 3 * (e - 1)
        EToV[e, 1] = pid[o + 1]
        EToV[e, 2] = pid[o + 2]
        EToV[e, 3] = pid[o + 3]
    end
    VX = [v[1] for v in verts]
    VY = [v[2] for v in verts]
    return VX, VY, EToV, elem_tree_elem
end

function _startupdg_local_face_for_t8_face(scheme, tree_class,
                                           elem::Ptr{T8code.Libt8.t8_element_t},
                                           iface, EToV, e, rd)
    c0 = T8code.Libt8.t8_element_get_face_corner(scheme, tree_class, elem, iface, Cint(0))
    c1 = T8code.Libt8.t8_element_get_face_corner(scheme, tree_class, elem, iface, Cint(1))
    g1 = EToV[e, c0 + 1]
    g2 = EToV[e, c1 + 1]
    sg = sort((g1, g2))
    for lf in 1:length(rd.fv)
        u = rd.fv[lf]
        if sort((EToV[e, u[1]], EToV[e, u[2]])) == sg
            return lf
        end
    end
    error("StartUpDGT8codeExt: could not map t8 face $iface to StartUpDG local face")
end

function _linear_face_node_indices(Nfp::Int, Nfaces::Int, e::Int, lf::Int)
    base = (e - 1) * Nfp * Nfaces + (lf - 1) * Nfp
    return (base + 1):(base + Nfp)
end

function _mapB_from_t8_boundaries(forest::t8_forest_t, rd::RefElemData{2},
                                  EToV::Matrix{Int},
                                  elem_tree_elem::Vector{Tuple{Int, Int}}, K::Int,
                                  Nfp::Int, Nfaces::Int)
    scheme = t8_forest_get_scheme(forest)
    forest_balanced = Cint(1)
    boundary_mask = falses(Nfp * Nfaces * K)
    for e in 1:K
        itree, ielement = elem_tree_elem[e]
        tree_class = t8_forest_get_tree_class(forest, itree)
        elem = t8_forest_get_leaf_element_in_tree(forest, itree, ielement)
        nface = Int(t8_element_get_num_faces(scheme, tree_class, elem))
        for iface in 0:(nface - 1)
            neighids_ref = Ref{Ptr{T8code.Libt8.t8_locidx_t}}()
            neighbors_ref = Ref{Ptr{Ptr{T8code.Libt8.t8_element_t}}}()
            neigh_scheme_ref = Ref{T8code.Libt8.t8_eclass_t}()
            dual_faces_ref = Ref{Ptr{Cint}}()
            num_neighbors_ref = Ref{Cint}()
            t8_forest_leaf_face_neighbors(forest, itree, elem, neighbors_ref, iface,
                                          dual_faces_ref, num_neighbors_ref, neighids_ref,
                                          neigh_scheme_ref, forest_balanced)
            num_neighbors = Int(num_neighbors_ref[])
            if neighbors_ref[] != C_NULL
                sc_free(t8_get_package_id(), neighbors_ref[])
            end
            if dual_faces_ref[] != C_NULL
                sc_free(t8_get_package_id(), dual_faces_ref[])
            end
            if neighids_ref[] != C_NULL
                sc_free(t8_get_package_id(), neighids_ref[])
            end
            if num_neighbors > 0
                continue
            end
            lf = _startupdg_local_face_for_t8_face(scheme, tree_class, elem, iface, EToV, e,
                                                   rd)
            boundary_mask[_linear_face_node_indices(Nfp, Nfaces, e, lf)] .= true
        end
    end
    return findall(boundary_mask)
end

function StartUpDG.MeshData(fw::T8codeForestMeshType, rd::RefElemData{2, <:Tri};
                            kwargs...)
    isempty(kwargs) ||
        throw(ArgumentError("MeshData(::T8codeForestMeshType, rd): unexpected keywords"))

    forest = pointer(fw.forest_wrapper)::t8_forest_t
    VX, VY, EToV, elem_tree_elem = _collect_triangle_mesh_arrays(forest)
    K = size(EToV, 1)
    (; fv) = rd
    FToF = connect_mesh(EToV, fv)
    Nfaces, K2 = size(FToF)
    @assert K == K2

    (; V1) = rd
    x = V1 * VX[transpose(EToV)]
    y = V1 * VY[transpose(EToV)]
    (; Vf) = rd
    xf = Vf * x
    yf = Vf * y
    mapM, mapP, mapB_geom = build_node_maps(FToF, (xf, yf))
    Nfp = size(Vf, 1) ÷ Nfaces
    mapM = reshape(mapM, Nfp * Nfaces, K)
    mapP = reshape(mapP, Nfp * Nfaces, K)

    mapB = _mapB_from_t8_boundaries(forest, rd, EToV, elem_tree_elem, K, Nfp, Nfaces)
    sort!(mapB)
    mapB_geom_sorted = sort(collect(mapB_geom))
    if mapB != mapB_geom_sorted
        error("StartUpDGT8codeExt: T8 boundary node set disagrees with connectivity mapP == mapM")
    end

    (; Dr, Ds) = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)
    (; Vq, wq) = rd
    xq, yq = (z -> Vq * z).((x, y))
    wJq = Diagonal(wq) * (Vq * J)
    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)
    periodicity = (false, false)

    md = MeshData(fw, FToF, tuple(x, y), tuple(xf, yf), tuple(xq, yq), wJq,
                  mapM, mapP, mapB,
                  rstxyzJ, J, tuple(nxJ, nyJ), Jf, periodicity)
    return md
end

end # module
