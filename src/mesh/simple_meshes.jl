"""
    function readGmsh2D(filename)

reads triangular GMSH 2D file format 2.2 0 8. returns VX,VY,EToV

# Examples
```julia
VXY,EToV = readGmsh2D("eulerSquareCylinder2D.msh")
```
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

    return VX,VY,EToV
end

"""
        uniform_mesh(elem::Line,Kx)
        uniform_mesh(elem::Tri,Kx,Ky)
        uniform_mesh(elem::Quad,Kx,Ky)
        uniform_mesh(elem::Hex,Kx,Ky,Kz)        

Uniform Kx (by Ky by Kz) mesh on ``[-1,1]^d``, where `d` is the spatial dimension.
Returns VX,VY,VZ,EToV. Can also use kwargs via `uniform_mesh(elem; K1D=16) `
"""
function uniform_mesh(elem::Line,K1D)
    VX = collect(LinRange(-1,1,K1D+1))
    EToV = transpose(reshape(sort([1:K1D; 2:K1D+1]),2,K1D))
    return VX,Matrix(EToV)
end

function uniform_mesh(elem::Tri,Kx,Ky)

        (VY, VX) = meshgrid(LinRange(-1,1,Ky+1),LinRange(-1,1,Kx+1))
        sk = 1
        EToV = zeros(Int,2*Kx*Ky,3)
        for ey = 1:Ky
                for ex = 1:Kx
                        id(ex,ey) = ex + (ey-1)*(Kx+1) # index function
                        id1 = id(ex,ey)
                        id2 = id(ex+1,ey)
                        id3 = id(ex+1,ey+1)
                        id4 = id(ex,ey+1)
                        EToV[2*sk-1,:] = [id1 id3 id2]
                        EToV[2*sk,:]   = [id3 id1 id4]
                        sk += 1
                end
        end
        return VX[:],VY[:],EToV
end

uniform_mesh(elem::Union{Tri,Quad},Kx) = uniform_mesh(elem,Kx,Kx)
uniform_mesh(elem::Hex,Kx) = uniform_mesh(elem,Kx,Kx,Kx)

function uniform_mesh(elem::Quad,Nx,Ny)

        Nxp = Nx+1;
        Nyp = Ny+1;
        Nv = convert(Int,Nxp*Nyp);
        K = convert(Int,Nx*Ny);

        x1D = LinRange(-1,1,Nxp);
        y1D = LinRange(-1,1,Nyp);
        x, y = meshgrid(x1D,y1D);
        I, J = meshgrid(collect(1:Nxp),collect(1:Nyp));
        inds = @. (I-1)*Ny + (J+I-1);
        EToV = zeros(Int,K,4);
        k = 1;
        for i = 1:Ny
                for j = 1:Nx
                        EToV[k,:] = [inds[i,j] inds[i,j+1] inds[i+1,j] inds[i+1,j+1]];
                        k += 1;
                end
        end

        VX = x[:];
        VY = y[:];

        return VX[:],VY[:],EToV
end

function uniform_mesh(elem::Hex,Nx,Ny,Nz)
        Nxp = Nx+1
        Nyp = Ny+1
        Nzp = Nz+1
        Nv = convert(Int,Nxp*Nyp*Nzp)
        K  = convert(Int,Nx*Ny*Nz)

        x1D = LinRange(-1,1,Nxp)
        y1D = LinRange(-1,1,Nyp)
        z1D = LinRange(-1,1,Nzp)

        # hack to work with old vertex ordering
        function mymeshgrid(x1D,y1D,z1D)
                Np = length(x1D)*length(y1D)*length(z1D)
                x = zeros(Np)
                y = zeros(Np)
                z = zeros(Np)
                sk = 1
                for k = 1:length(z1D)
                        for j = 1:length(y1D)
                                for i = 1:length(x1D)
                                        x[sk] = x1D[i]
                                        y[sk] = y1D[j]
                                        z[sk] = z1D[k]
                                        sk += 1
                                end
                        end
                end
                return x,y,z
        end
        x, y, z = mymeshgrid(x1D,y1D,z1D)

        EToV = zeros(Int,K,8)
        k = 1
        for e = 1:K
                em = e-1
                k = div(em,(Nx*Ny))
                j = div(em - k*Nx*Ny,Nx)
                i = em % Nx

                EToV[e,1] = i     + Nxp*j     + Nxp*Nyp*k
                EToV[e,2] = i     + Nxp*(j+1) + Nxp*Nyp*k
                EToV[e,3] = (i+1) + Nxp*j     + Nxp*Nyp*k
                EToV[e,4] = (i+1) + Nxp*(j+1) + Nxp*Nyp*k
                EToV[e,5] = i     + Nxp*j     + Nxp*Nyp*(k+1)
                EToV[e,6] = i     + Nxp*(j+1) + Nxp*Nyp*(k+1)
                EToV[e,7] = (i+1) + Nxp*j     + Nxp*Nyp*(k+1)
                EToV[e,8] = (i+1) + Nxp*(j+1) + Nxp*Nyp*(k+1)
        end
        EToV = @. EToV + 1 # re-index to 1 index

        VX = x[:];
        VY = y[:];
        VZ = z[:];
        return VX[:],VY[:],VZ[:],EToV
end

# keyword argument version
uniform_mesh(elem::T; K1D) where {T <: AbstractElemShape} = uniform_mesh(elem,K1D)

