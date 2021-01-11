"""
    connect_mesh(EToV,fv)

Initialize element connectivity matrices, element to element and element to face
connectivity. `EToV` is a `K` by `Nv` matrix whose rows identify the `Nv` vertices
which make up an element.

`fv` (an array of arrays containing unordered indices of face vertices).
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
    init_DG_mesh(VX::AbstractArray,EToV,rd::RefElemData)
    init_DG_mesh(VX,VY,EToV,rd::RefElemData)
    init_DG_mesh(VX,VY,VZ,EToV,rd::RefElemData)

Initializes mesh data (node connectivity, geometric terms, etc).
"""
# distinguish 1D mesh construction from 2D/3D
# assumes VX is ordered left-to-right
function init_DG_mesh(VX::AbstractArray,EToV,rd::RefElemData)
    md = MeshData()

    # Construct global coordinates
    @unpack V1 = rd
    x = V1*VX[transpose(EToV)]
    K = size(EToV,1)
    @pack! md = K,x

    # Connectivity maps
    @unpack Vf = rd
    xf = Vf*x
    mapM = reshape(1:2*K,2,K)
    mapP = copy(mapM)
    mapP[1,2:end] .= mapM[2,1:end-1]
    mapP[2,1:end-1] .= mapM[1,2:end]
    mapB = findall(@. mapM[:]==mapP[:])
    @pack! md = xf,mapM,mapP,mapB

    # # Make periodic
    # mapP[1] = mapM[end]
    # mapP[end] = mapM[1]

    # Geometric factors and surface normals
    J = repeat(transpose(diff(VX)/2),length(rd.r),1)
    rxJ = one.(J)
    nxJ = repeat([-1;1],1,K)
    sJ = abs.(nxJ)
    @pack! md = J,rxJ,nxJ,sJ

    @unpack Vq,wq = rd
    xq = Vq*x
    wJq = diagm(wq)*(Vq*J)
    @pack! md = xq,wJq

    return md
end

function init_DG_mesh(VX,VY,EToV,rd::RefElemData)
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


function init_DG_mesh(VX,VY,VZ,EToV,rd::RefElemData)
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
