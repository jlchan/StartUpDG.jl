"""
    geometric_factors(x, y, Dr, Ds)
    geometric_factors(x, y, z, Dr, Ds, Dt, Filters=(I,I,I))

Compute metrics of mappings between "real" elements and reference elements,
outward pointing normals on faces of every elements, and Jacobian.

x,y,z are arrays of coordinates, and Dr, Ds, Dt are nodal differentiation matrices
Filters = tuple of filtering matrices (e.g., to reduce degree in r,s, and t for GCL)

Geometric terms in 3D are constructed to ensure satisfaction of free-stream preservation
using the curl-based construction from 'Metric identities and the DG-SEM on curvilinear
meshes' (Kopriva 2006).
"""

function geometric_factors(x, y, Dr, Ds)
    xr, xs = Dr*x, Ds*x
    yr, ys = Dr*y, Ds*y

    J = @. -xs * yr + xr * ys
    rxJ, sxJ =  ys, -yr
    ryJ, syJ = -xs,  xr

    return rxJ, sxJ, ryJ, syJ, J
end

function geometric_factors(x, y, z, Dr, Ds, Dt, Filters=(I,I,I))

    xr, xs, xt = Dr*x, Ds*x, Dt*x
    yr, ys, yt = Dr*y, Ds*y, Dt*y
    zr, zs, zt = Dr*z, Ds*z, Dt*z

    Fr, Fs, Ft = (Dr*y).*z, (Ds*y).*z, (Dt*y).*z
    Fr,Fs,Ft = ((A,x)->A*x).(Filters,(Fr,Fs,Ft))
    rxJ = Dt*(Fs) - Ds*(Ft)
    sxJ = Dr*(Ft) - Dt*(Fr)
    txJ = Ds*(Fr) - Dr*(Fs)

    Fr, Fs, Ft = (Dr*x).*z, (Ds*x).*z, (Dt*x).*z
    Fr,Fs,Ft = ((A,x)->A*x).(Filters,(Fr,Fs,Ft))
    ryJ = -(Dt*(Fs) - Ds*(Ft))
    syJ = -(Dr*(Ft) - Dt*(Fr))
    tyJ = -(Ds*(Fr) - Dr*(Fs))

    Fr, Fs, Ft = (Dr*y).*x, (Ds*y).*x, (Dt*y).*x
    Fr,Fs,Ft = ((A,x)->A*x).(Filters,(Fr,Fs,Ft))
    rzJ = -(Dt*(Fs) - Ds*(Ft))
    szJ = -(Dr*(Ft) - Dt*(Fr))
    tzJ = -(Ds*(Fr) - Dr*(Fs))

    J = @. xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) + zr * (xs * yt - ys * xt)

    return rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J
end

"""
    estimate_h(rd::RefElemData,md::MeshData)

Estimates the mesh size via min size_of_domain * |J|/|sJ|, since |J| = O(hᵈ) and |sJ| = O(hᵈ⁻¹). 
"""
function estimate_h(rd::RefElemData{DIM}, md::MeshData{DIM}) where {DIM}
    hmin = Inf
    for e in 1:md.num_elements
        sJ_e = reshape(view(md.sJ, :, e), rd.Nfq ÷ rd.Nfaces, rd.Nfaces)
        sJ_face = zero(eltype(md.sJ))
        for f in 1:rd.Nfaces
            sJ_face = max(sJ_face, minimum(view(sJ_e, :, f)) / face_scaling(rd, f))
        end
        h_e = minimum(view(md.J, :, e)) / sJ_face
        hmin = min(hmin, h_e)
    end
    return hmin * compute_domain_size(rd, md)^(1/DIM)
end
face_scaling(rd, f) = 1.0
face_scaling(rd::RefElemData{2, Tri}, f) = f == 3 ? sqrt(2) : 1.0 # sJ incorporates length of long triangle edge
face_scaling(rd::RefElemData{3, Tet}, f) = f == 2 ? sqrt(3) : 1.0 # sJ incorporates area of larger triangle face
compute_domain_size(rd::RefElemData, md::MeshData) = sum(rd.M * md.J)

