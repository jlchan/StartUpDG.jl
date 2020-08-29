"""
connect_mesh(EToV,fv)

Initialize element connectivity matrices, element to element and
element to face connectivity. Works on general element types, so long as

        fv = array of arrays containing (unordered) indices of face vertices

is provided.

# Examples
- connect_mesh(EToV,UniformTriMesh.fv())
- connect_mesh(EToV,[[1,2],[2,3],[3,1]])
```jldoctest
```
"""
function connect_mesh(EToV,fv)
        Nfaces = length(fv)
        K = size(EToV,1)

        # sort and find matches
        fnodes = [[sort(EToV[e,ids]) for ids = fv, e = 1:K]...]
        p = sortperm(fnodes) # sorts by lexicographic ordering by default
        fnodes = fnodes[p,:]

        FToF = reshape(collect(1:Nfaces*K),Nfaces,K)
        for f = 1:size(fnodes,1)-1
                if fnodes[f,:]==fnodes[f+1,:]
                        f1 = FToF[p[f]]
                        f2 = FToF[p[f+1]]
                        FToF[p[f]] = f2
                        FToF[p[f+1]] = f1
                end
        end
        return FToF
end


"""
    function readGmsh2D(filename)

reads GMSH 2D file format 2.2 0 8
returns EToV,VX,VY

# Examples

EToV,VX,VY = readGmsh2D("eulerSquareCylinder2D.msh")

```jldoctest
"""
function readGmsh2D(filename)
    f = open(filename)
    lines = readlines(f)

    function findline(name,lines)
        for (i,line) in enumerate(lines)
            if line==name
                return i
            end
        end
    end

    node_start = findline("\$Nodes",lines)+1
    Nv = parse(Int64,lines[node_start])
    VX,VY,VZ = ntuple(x->zeros(Float64,Nv),3)
    for i = 1:Nv
        vals = [parse(Float64,c) for c in split(lines[i+node_start])]
        # first entry =
        VX[i] = vals[2]
        VY[i] = vals[3]
    end

    elem_start = findline("\$Elements",lines)+1
    K_all      = parse(Int64,lines[elem_start])
    K = 0
    for e = 1:K_all
        if length(split(lines[e+elem_start]))==8
            K = K + 1
        end
    end
    EToV = zeros(Int64,K,3)
    sk = 1
    for e = 1:K_all
        if length(split(lines[e+elem_start]))==8
            vals = [parse(Int64,c) for c in split(lines[e+elem_start])]
            EToV[sk,:] .= vals[6:8]
            sk = sk + 1
        end
    end

    EToV = EToV[:,vec([1 3 2])] # permute for Gmsh ordering

    return EToV,VX,VY
end


# annotate types for geofacs + connectivity arrays for speed in RHS evals
mutable struct MeshData
    VX; VY; VZ              # vertex coordinates
    K::Int                  # num elems
    EToV                    # mesh vertex array
    FToF::Array{Int64,2}    # face connectivity

    x; y; z                 # physical points
    xf; yf; zf
    xq; yq; zq;             # phys quad points, Jacobian-scaled weights
    wJq::Array{Float64,2}

    # arrays of connectivity indices between face nodes
    mapM
    mapP::Array{Int64,2}
    mapB::Array{Int64,1}

    # volume geofacs
    rxJ::Array{Float64,2}
    sxJ::Array{Float64,2}
    txJ::Array{Float64,2}
    ryJ::Array{Float64,2}
    syJ::Array{Float64,2}
    tyJ::Array{Float64,2}
    rzJ::Array{Float64,2}
    szJ::Array{Float64,2}
    tzJ::Array{Float64,2}
    J::Array{Float64,2}

    # surface geofacs
    nxJ::Array{Float64,2}
    nyJ::Array{Float64,2}
    nzJ::Array{Float64,2}
    sJ::Array{Float64,2}

    MeshData() = new() # empty initializer
end


# dispatch to 2D or 3D version if tuple called
function init_mesh(VXYZ,EToV,rd::RefElemData)
    return init_mesh(VXYZ...,EToV,rd)
end

function init_mesh(VX,VY,EToV,rd::RefElemData)
    # initialize a new mesh data struct
    md = MeshData()

    @unpack fv = rd
    FToF = connect_mesh(EToV,fv)
    Nfaces,K = size(FToF)
    @pack! md = FToF,K,VX,VY,EToV

    #Construct global coordinates
    @unpack V1 = rd
    x = V1*VX[transpose(EToV)]
    y = V1*VY[transpose(EToV)]
    @pack! md = x,y

    #Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    @unpack r,s,Vf = rd
    xf,yf = (x->Vf*x).((x,y))
    mapM,mapP,mapB = build_node_maps((xf,yf),FToF)
    Nfp = convert(Int,size(Vf,1)/Nfaces)
    mapM = reshape(mapM,Nfp*Nfaces,K)
    mapP = reshape(mapP,Nfp*Nfaces,K)
    @pack! md = xf,yf,mapM,mapP,mapB

    #Compute geometric factors and surface normals
    @unpack Dr,Ds = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    @pack! md = rxJ, sxJ, ryJ, syJ, J

    @unpack Vq,wq = rd
    xq,yq = (x->Vq*x).((x,y))
    wJq = diagm(wq)*(Vq*J)
    @pack! md = xq,yq,wJq

    #physical normals are computed via G*nhatJ, G = matrix of geometric terms
    @unpack nrJ,nsJ = rd
    nxJ = (Vf*rxJ).*nrJ + (Vf*sxJ).*nsJ
    nyJ = (Vf*ryJ).*nrJ + (Vf*syJ).*nsJ
    sJ = @. sqrt(nxJ^2 + nyJ^2)
    @pack! md = nxJ,nyJ,sJ

    return md
end

"init mesh in 3D"
function init_mesh(VX,VY,VZ,EToV,rd::RefElemData)
    # initialize a new mesh data struct
    md = MeshData()

    @unpack fv = rd
    FToF = connect_mesh(EToV,fv)
    Nfaces,K = size(FToF)
    @pack! md = FToF,K,VX,VY,VZ,EToV

    #Construct global coordinates
    @unpack V1 = rd
    x = V1*VX[transpose(EToV)]
    y = V1*VY[transpose(EToV)]
    z = V1*VZ[transpose(EToV)]
    @pack! md = x,y,z

    #Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    @unpack r,s,t,Vf = rd
    xf,yf,zf = (x->Vf*x).((x,y,z))
    mapM,mapP,mapB = build_node_maps((xf,yf,zf),FToF)
    Nfp = convert(Int,size(Vf,1)/Nfaces)
    mapM = reshape(mapM,Nfp*Nfaces,K)
    mapP = reshape(mapP,Nfp*Nfaces,K)
    @pack! md = xf,yf,zf,mapM,mapP,mapB

    #Compute geometric factors and surface normals
    @unpack Dr,Ds,Dt = rd
    rxJ,sxJ,txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ,J = geometric_factors(x,y,z,Dr,Ds,Dt)
    @pack! md = rxJ,sxJ,txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ,J

    @unpack Vq,wq = rd
    xq,yq,zq = (x->Vq*x).((x,y,z))
    wJq = diagm(wq)*(Vq*J)
    @pack! md = xq,yq,zq,wJq

    #physical normals are computed via G*nhatJ, G = matrix of geometric terms
    @unpack nrJ,nsJ,ntJ = rd
    nxJ = nrJ.*(Vf*rxJ) + nsJ.*(Vf*sxJ) + ntJ.*(Vf*txJ)
    nyJ = nrJ.*(Vf*ryJ) + nsJ.*(Vf*syJ) + ntJ.*(Vf*tyJ)
    nzJ = nrJ.*(Vf*rzJ) + nsJ.*(Vf*szJ) + ntJ.*(Vf*tzJ)
    sJ = @. sqrt(nxJ.^2 + nyJ.^2 + nzJ.^2)
    @pack! md = nxJ,nyJ,nzJ,sJ

    return md
end
