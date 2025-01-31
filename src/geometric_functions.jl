"""
    geometric_factors(x, y, Dr, Ds)
    geometric_factors(x, y, z, Dr, Ds, Dt, Filters=(I, I, I))

Compute metrics of mappings between "real" elements and reference elements,
outward pointing normals on faces of every elements, and Jacobian.

x,y,z are arrays of coordinates, and Dr, Ds, Dt are nodal differentiation matrices
Filters = tuple of filtering matrices (e.g., to reduce degree in r,s, and t for GCL)

Geometric terms in 3D are constructed to ensure satisfaction of free-stream preservation
using the curl-based construction from 'Metric identities and the DG-SEM on curvilinear
meshes' (Kopriva 2006).
"""

function geometric_factors(x, y, Dr, Ds)
    xr, xs = Dr * x, Ds * x
    yr, ys = Dr * y, Ds * y

    rxJ, sxJ =  ys, -yr
    ryJ, syJ = -xs,  xr

    J = @. -xs * yr + xr * ys

    return rxJ, sxJ, ryJ, syJ, J
end

function geometric_factors(x, y, z, Dr, Ds, Dt, Filters=(I, I, I))

    xr, xs, xt = Dr * x, Ds * x, Dt * x
    yr, ys, yt = Dr * y, Ds * y, Dt * y
    zr, zs, zt = Dr * z, Ds * z, Dt * z

    Fr, Fs, Ft = (Dr * y) .* z, (Ds * y) .* z, (Dt * y) .* z
    Fr, Fs, Ft = ((A, x)-> A * x).(Filters, (Fr, Fs, Ft))
    rxJ = Dt * Fs - Ds * Ft
    sxJ = Dr * Ft - Dt * Fr
    txJ = Ds * Fr - Dr * Fs

    Fr, Fs, Ft = (Dr * x) .* z, (Ds * x) .* z, (Dt * x) .* z
    Fr, Fs, Ft = ((A, x) -> A * x).(Filters, (Fr, Fs, Ft))
    ryJ = -(Dt * Fs - Ds * Ft)
    syJ = -(Dr * Ft - Dt * Fr)
    tyJ = -(Ds * Fr - Dr * Fs)

    Fr, Fs, Ft = (Dr * y) .* x, (Ds * y) .* x, (Dt * y) .* x
    Fr, Fs, Ft = ((A, x) -> A * x).(Filters, (Fr, Fs, Ft))
    rzJ = -(Dt * Fs - Ds * Ft)
    szJ = -(Dr * Ft - Dt * Fr)
    tzJ = -(Ds * Fr - Dr * Fs)

    J = @. xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) + zr * (xs * yt - ys * xt)

    return rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J
end

# physical outward normals are computed via Nanson's formula: G * nhatJ, where 
# G = matrix of J-scaled geometric terms. Here, Vf is a face interpolation matrix 
# which maps interpolation nodes to face nodes. 
function compute_normals(geo::SMatrix{Dim, Dim}, Vf, nrstJ...) where {Dim}
    nxyzJ = ntuple(x -> zeros(size(Vf, 1), size(first(geo), 2)), Dim)
    for i = 1:Dim, j = 1:Dim
        nxyzJ[i] .+= (Vf * geo[i,j]) .* nrstJ[j]
    end
    Jf = sqrt.(sum(map(x -> x.^2, nxyzJ)))
    return nxyzJ..., Jf
end

"""
    estimate_h(rd::RefElemData, md::MeshData)
    estimate_h(e, rd::RefElemData, md::MeshData) # e = element index

Estimates the mesh size via min size_of_domain * |J|/|Jf|, since |J| = O(hᵈ) and |Jf| = O(hᵈ⁻¹). 
"""
function estimate_h(rd::RefElemData{DIM}, md::MeshData{DIM}) where {DIM}
    hmin = Inf
    for e in 1:md.num_elements
        h_e = estimate_h(e, rd, md)
        hmin = min(hmin, h_e)
    end
    domain_size = sum(rd.M * md.J)
    return hmin * domain_size^(1/DIM)
end

function estimate_h(e, rd::RefElemData, md::MeshData)
    Jf_e = reshape(view(md.Jf, :, e), :, rd.Nfaces)
    Jf_face = zero(eltype(md.Jf))
    for f in 1:rd.Nfaces
        Jf_face = max(Jf_face, minimum(view(Jf_e, :, f)) / face_scaling(rd, f))
    end
    h_e = minimum(view(md.J, :, e)) / Jf_face
    return h_e
end

# specialization for elements with different face types
function estimate_h(e, rd::RefElemData{3, <:Union{Wedge, Pyr}}, md::MeshData)
    (; node_ids_by_face ) = rd.element_type
    Jf_e = view(md.Jf, :, e)
    Jf_face = zero(eltype(md.Jf))
    for f in 1:rd.num_faces
        Jf_face = max(Jf_face, minimum(view(Jf_e, node_ids_by_face[f])) / face_scaling(rd, f))
    end
    h_e = minimum(view(md.J, :, e)) / Jf_face
    return h_e
end

face_scaling(rd, f) = 1.0
face_scaling(rd::RefElemData{2, Tri}, f) = f == 3 ? sqrt(2) : 1.0 # Jf incorporates length of long triangle edge
face_scaling(rd::RefElemData{3, Tet}, f) = f == 2 ? sqrt(3) : 1.0 # Jf incorporates area of larger triangle face
face_scaling(rd::RefElemData{3, Wedge}, f) = f == 2 ? sqrt(2) : 1.0 # Jf incorporates area of larger triangle face

