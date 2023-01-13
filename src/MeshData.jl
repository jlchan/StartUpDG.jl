"""
    struct MeshData{Dim, Tv, Ti}

MeshData: contains info for a high order piecewise polynomial discretization on an
unstructured mesh. 

Example:
```julia
N, K1D = 3, 2
rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(Tri(), K1D)
md = MeshElemData(VXY, EToV, rd)
@unpack x, y = md
```
"""
Base.@kwdef struct MeshData{Dim, MeshType, 
                            VolumeType, FaceType, VolumeQType, FToFType, 
                            VolumeWeightType, VolumeGeofacsType, VolumeJType, 
                            ConnectivityTypeM, ConnectivityTypeP, BoundaryMapType}

    # this field defaults to the element shape, but can be used to specify 
    # different mesh types (e.g., vertex-mapped, cut-cell, hybrid, non-conforming meshes. 
    mesh_type::MeshType 

    FToF::FToFType                  # face connectivity

    xyz::NTuple{Dim, VolumeType}    # physical points
    xyzf::NTuple{Dim, FaceType}     # face nodes
    xyzq::NTuple{Dim, VolumeQType}  # phys quad points, Jacobian-scaled weights
    wJq::VolumeWeightType

    # arrays of connectivity indices between face nodes
    mapM::ConnectivityTypeM
    mapP::ConnectivityTypeP
    mapB::BoundaryMapType

    # volume geofacs Gij = dx_i/dxhat_j
    rstxyzJ::VolumeGeofacsType
    J::VolumeJType

    # normalized surface geofacs
    nxyz::NTuple{Dim, FaceType}

    # surface geofacs
    nxyzJ::NTuple{Dim, FaceType}
    Jf::FaceType

    is_periodic::NTuple{Dim, Bool}
end

# MeshData constructor where we do not specify `nxyz` and instead compute it from `nrstJ` and `Jf`
function MeshData(mesh_type, FToF, xyz, xyzf, xyzq, wJq, mapM, mapP, mapB, rstxyzJ, J, nxyzJ, Jf, is_periodic) 

    nxyz = map(nJ -> nJ ./ Jf, nxyzJ)        
                         
    return MeshData(mesh_type, FToF, xyz, xyzf, xyzq, wJq, mapM, mapP, mapB, rstxyzJ, J, nxyz, nxyzJ, Jf, is_periodic)
end

# TODO: remove in next breaking release after v0.15. `MeshData` constructor where `VXYZ`, `EToV` are specified 
@deprecate MeshData(element_type, VXYZ, EToV, FToF, xyz, xyzf, xyzq, wJq, mapM, mapP, mapB, rstxyzJ, J, nxyzJ, Jf, periodicity) MeshData(VertexMappedMesh(element_type, VXYZ, EToV), FToF, xyz, xyzf, xyzq, wJq, mapM, mapP, mapB, rstxyzJ, J, nxyzJ, Jf, periodicity)

function ConstructionBase.setproperties(md::MeshData, patch::NamedTuple)
    fields = (haskey(patch, symbol) ? getproperty(patch, symbol) : getproperty(md, symbol) for symbol in fieldnames(typeof(md)))
    return MeshData(fields...)
end

ConstructionBase.getproperties(md::MeshData) = 
    (; mesh_type=md.mesh_type, FToF=md.FToF, xyz=md.xyz, xyzf=md.xyzf, xyzq=md.xyzq, wJq=md.wJq,
       mapM=md.mapM, mapP=md.mapP, mapB=md.mapB, rstxyzJ=md.rstxyzJ, J=md.J, nxyz = md.nxyz, nxyzJ=md.nxyzJ, Jf=md.Jf,
       is_periodic=md.is_periodic)

function Base.show(io::IO, md::MeshData{DIM}) where {DIM}
    @nospecialize md
    print(io,"MeshData{$DIM}")
end

function Base.show(io::IO, ::MIME"text/plain", md::MeshData{DIM, MeshType}) where {DIM, MeshType}
    @nospecialize md
    print(io,"MeshData of dimension $DIM with $(md.num_elements) elements")
end

function Base.propertynames(x::MeshData{1}, private::Bool = false)
    return (fieldnames(MeshData)...,
            :num_elements, :VX, :x, :xq, :xf, :nx, :rxJ)
end
function Base.propertynames(x::MeshData{2}, private::Bool = false) 
    return (fieldnames(MeshData)...,
            :num_elements, :VX, :VY, :x, :y, :xq, :yq, :xf, :yf, 
            :nx, :ny, :rxJ, :sxJ, :ryJ, :syJ)
end
function Base.propertynames(x::MeshData{3}, private::Bool = false) 
    return (fieldnames(MeshData)...,
            :num_elements, :VX, :VY, :VZ, :x, :y, :z, :xq, :yq, :zq, :xf, :yf, :zf, 
            :nx, :ny, :nz, :rxJ, :sxJ, :txJ, :ryJ, :syJ, :tyJ, :rzJ, :szJ, :tzJ)
end

# seems like we need to use @inline to ensure type stability
@inline function meshdata_getproperty(x::MeshData, s::Symbol)
    if s==:x
        return getfield(x, :xyz)[1]
    elseif s==:y
        return getfield(x, :xyz)[2]
    elseif s==:z
        return getfield(x, :xyz)[3]

    elseif s==:xq
        return getfield(x, :xyzq)[1]
    elseif s==:yq
        return getfield(x, :xyzq)[2]
    elseif s==:zq
        return getfield(x, :xyzq)[3]

    elseif s==:xf
        return getfield(x, :xyzf)[1]
    elseif s==:yf
        return getfield(x, :xyzf)[2]
    elseif s==:zf
        return getfield(x, :xyzf)[3]

    elseif s==:nxJ
        return getfield(x, :nxyzJ)[1]
    elseif s==:nyJ
        return getfield(x, :nxyzJ)[2]
    elseif s==:nzJ
        return getfield(x, :nxyzJ)[3]

    elseif s==:nx
        return getfield(x, :nxyz)[1]
    elseif s==:ny
        return getfield(x, :nxyz)[2]
    elseif s==:nz
        return getfield(x, :nxyz)[3]        

    elseif s==:rxJ
        return getfield(x, :rstxyzJ)[1,1]
    elseif s==:sxJ
        return getfield(x, :rstxyzJ)[1,2]
    elseif s==:txJ
        return getfield(x, :rstxyzJ)[1,3]
    elseif s==:ryJ
        return getfield(x, :rstxyzJ)[2,1]
    elseif s==:syJ
        return getfield(x, :rstxyzJ)[2,2]
    elseif s==:tyJ
        return getfield(x, :rstxyzJ)[2,3]
    elseif s==:rzJ
        return getfield(x, :rstxyzJ)[3,1]
    elseif s==:szJ
        return getfield(x, :rstxyzJ)[3,2]
    elseif s==:tzJ
        return getfield(x, :rstxyzJ)[3,3]
        
    # old behavior where K = num_elements        
    elseif s==:K || s==:num_elements 
        return num_elements(x)

    # old notation in the NDG book where sJ (surface Jacobian) is 
    # used instead of Jf (Jacobian for the face)
    elseif s==:sJ 
        return getfield(x, :Jf)

    else
        return getfield(x, s)
    end
end

# generic fallback
Base.getproperty(x::MeshData, s::Symbol) = meshdata_getproperty(x, s)

"""
    struct VertexMappedMesh

The default `MeshData` mesh type, represents a mesh which is defined purely by 
vertex locations and element-to-vertex connectivities. For example, these include 
affine triangular meshes or bilinear quadrilateral or trilinear hexahedral meshes.

# Fields
element_type :: TE <: AbstractElemShape \\
VXYZ :: TV \\
EToV :: TEV
"""
struct VertexMappedMesh{TE <: AbstractElemShape, TV, TEV}
    element_type::TE
    VXYZ::TV
    EToV::TEV
end

# convenience accessor routines for `VertexMappedMesh` types (lets you do `md.VX` instead of `md.mesh_type.VX`)
function Base.getproperty(x::MeshData{Dim, <:VertexMappedMesh}, s::Symbol) where {Dim}

    if s===:VX
        return getfield(getfield(x, :mesh_type), :VXYZ)[1]
    elseif s===:VY
        return getfield(getfield(x, :mesh_type), :VXYZ)[2]
    elseif s===:VZ
        return getfield(getfield(x, :mesh_type), :VXYZ)[3]
    elseif s===:EToV
        return getfield(getfield(x, :mesh_type), s)
    else
        meshdata_getproperty(x, s)
    end
end

num_elements(md::MeshData) = size(md.x, 2) # number of columns in the "x" coordinate array
num_elements(md::MeshData{Dim, <:VertexMappedMesh}) where {Dim} = size(md.mesh_type.EToV, 1)

"""
    MeshData(VXYZ, EToV, rd::RefElemData)

Returns a MeshData struct with high order DG mesh information from the unstructured
mesh information (VXYZ..., EToV).

    MeshData(rd::RefElemData, md::MeshData, xyz...)

Given new nodal positions `xyz...` (e.g., from mesh curving), recomputes geometric terms
and outputs a new MeshData struct. Only fields modified are the coordinate-dependent terms
    `xyz`, `xyzf`, `xyzq`, `rstxyzJ`, `J`, `nxyzJ`, `Jf`.
"""

# splats VXYZ 
MeshData(VXYZ::T, EToV, other_args...) where {NDIMS, T <: NTuple{NDIMS}} = 
    MeshData(VXYZ..., EToV, other_args...)

function MeshData(VX::AbstractVector{Tv}, EToV, rd::RefElemData{1}) where {Tv}

    # Construct global coordinates
    @unpack V1 = rd
    x = V1 * VX[transpose(EToV)]
    K = size(EToV, 1)
    Nfaces = 2

    FToF = zeros(Int, Nfaces, K)
    sk = 1
    for e = 1:K
        l = 2 * e-1
        r = 2 * e
        FToF[1:2, e] .= [l-1; r+1]
        sk += 1
    end
    FToF[1, 1] = 1
    FToF[Nfaces, K] = Nfaces * K

    # Connectivity maps
    @unpack Vf = rd
    xf = Vf * x
    mapM = reshape(1:2*K, 2, K)
    mapP = copy(mapM)
    mapP[1, 2:end] .= mapM[2, 1:end-1]
    mapP[2, 1:end-1] .= mapM[1, 2:end]
    mapB = findall(@. mapM[:] == mapP[:])

    # Geometric factors and surface normals
    J = repeat(transpose(diff(VX) / 2), length(rd.r), 1)
    rxJ = one.(J)
    nxJ = repeat([-1.0; 1.0], 1, K)
    Jf = abs.(nxJ)

    @unpack Vq, wq = rd
    xq = Vq * x
    wJq = diagm(wq) * (Vq * J)

    is_periodic = (false,)

    return MeshData(VertexMappedMesh(rd.element_type, tuple(VX), EToV), FToF,
                    tuple(x), tuple(xf), tuple(xq), wJq,
                    collect(mapM), mapP, mapB,
                    SMatrix{1,1}(tuple(rxJ)), J,
                    tuple(nxJ), Jf,
                    is_periodic)
end

function MeshData(VX, VY, EToV, rd::RefElemData{2})

    @unpack fv = rd
    FToF = connect_mesh(EToV, fv)
    Nfaces, K = size(FToF)

    #Construct global coordinates
    @unpack V1 = rd
    x = V1 * VX[transpose(EToV)]
    y = V1 * VY[transpose(EToV)]

    # Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    @unpack Vf = rd
    xf = Vf * x
    yf = Vf * y
    mapM, mapP, mapB = build_node_maps(FToF, (xf, yf))
    Nfp = size(Vf, 1) ÷ Nfaces
    mapM = reshape(mapM, Nfp * Nfaces, K)
    mapP = reshape(mapP, Nfp * Nfaces, K)

    #Compute geometric factors and surface normals
    @unpack Dr, Ds = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)

    @unpack Vq, wq = rd
    xq, yq = (x -> Vq * x).((x, y))
    wJq = diagm(wq) * (Vq * J)

    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)

    is_periodic = (false, false)
    return MeshData(VertexMappedMesh(rd.element_type, tuple(VX, VY), EToV), FToF,
                    tuple(x, y), tuple(xf, yf), tuple(xq, yq), wJq,
                    mapM, mapP, mapB,
                    SMatrix{2, 2}(tuple(rxJ, ryJ, sxJ, syJ)), J,
                    tuple(nxJ, nyJ), Jf,
                    is_periodic)

end

function MeshData(VX, VY, VZ, EToV, rd::RefElemData{3})

    @unpack fv = rd
    FToF = connect_mesh(EToV, fv)
    Nfaces, K = size(FToF)

    #Construct global coordinates
    @unpack V1 = rd
    x, y, z = (x -> V1 * x[transpose(EToV)]).((VX, VY, VZ))

    #Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    @unpack r, s, t, Vf = rd
    xf, yf, zf = (x -> Vf * x).((x, y, z))
    mapM, mapP, mapB = build_node_maps(rd, FToF, (xf, yf, zf))
    mapM = reshape(mapM, :, K)
    mapP = reshape(mapP, :, K)
    
    #Compute geometric factors and surface normals
    @unpack Dr, Ds, Dt = rd
    rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = geometric_factors(x, y, z, Dr, Ds, Dt)
    rstxyzJ = SMatrix{3, 3}(rxJ, ryJ, rzJ, sxJ, syJ, szJ, txJ, tyJ, tzJ)

    @unpack Vq, wq = rd
    xq, yq, zq = (x -> Vq * x).((x, y, z))
    wJq = diagm(wq) * (Vq * J)

    nxJ, nyJ, nzJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)

    is_periodic = (false, false, false)
    return MeshData(VertexMappedMesh(rd.element_type, tuple(VX, VY, VZ), EToV), FToF,
                    tuple(x, y, z), tuple(xf, yf, zf), tuple(xq, yq, zq), wJq,
                    mapM, mapP, mapB,
                    rstxyzJ, J, tuple(nxJ, nyJ, nzJ), Jf,
                    is_periodic)
end

function recompute_geometry(rd::RefElemData{Dim}, xyz) where {Dim}
    # compute new quad and plotting points
    xyzf = map(x -> rd.Vf * x, xyz)
    xyzq = map(x -> rd.Vq * x, xyz)

    #Compute geometric factors and surface normals
    geo = geometric_factors(xyz..., rd.Drst...)
    J = last(geo)
    wJq = diagm(rd.wq) * (rd.Vq * J)

    if Dim==1
        rstxyzJ = SMatrix{Dim, Dim}(geo[1])
    elseif Dim==2
        rstxyzJ = SMatrix{Dim, Dim}(geo[1], geo[3],
                                    geo[2], geo[4])
    elseif Dim==3
        rstxyzJ = SMatrix{Dim, Dim}(geo[1], geo[4], geo[7],
                                    geo[2], geo[5], geo[8],
                                    geo[3], geo[6], geo[9])
    end
    geof = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)
    nxyzJ = geof[1:Dim]
    Jf = last(geof)
    return xyzf, xyzq, rstxyzJ, J, wJq, nxyzJ, Jf
end

"""
    struct CurvedMesh

Mesh type indicating that the mesh has been curved. Stores the original mesh type as a field.

# Fields
original_mesh_type :: T\\
"""
struct CurvedMesh{T}
    original_mesh_type::T
end 

function MeshData(rd::RefElemData, md::MeshData{Dim}, xyz...) where {Dim}

    xyzf, xyzq, rstxyzJ, J, wJq, nxyzJ, Jf = recompute_geometry(rd, xyz)

    nxyz = map(n -> n ./ Jf, nxyzJ)

    mesh_type = CurvedMesh(md.mesh_type)

    # TODO: should we warp VXYZ as well? Or create a new mesh type?
    return setproperties(md, (; mesh_type, xyz, xyzq, xyzf, rstxyzJ, J, wJq, nxyz, nxyzJ, Jf))
end
