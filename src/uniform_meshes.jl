###########################
### Triangular meshes #####
###########################

"""
uniform_tri_mesh(Kx::Int,Ky::Int)

Matlab uniform triangular mesh.

# Examples
```jldoctest
```
"""

function uniform_tri_mesh(Kx,Ky)
        (VY, VX) = meshgrid(LinRange(-1,1,Ky+1),LinRange(-1,1,Kx+1))
        sk = 1
        EToV = zeros(Int,2*Kx*Ky,3)
        for ey = 1:Ky
                for ex = 1:Kx
                        id(ex,ey) = ex + (ey-1)*(Kx+1) # index function
                        id1 = id(ex,ey);
                        id2 = id(ex+1,ey);
                        id3 = id(ex+1,ey+1);
                        id4 = id(ex,ey+1);
                        EToV[2*sk-1,:] = [id1 id2 id3];
                        EToV[2*sk,:] = [id3 id4 id1];
                        sk += 1
                end
        end
        return (VX[:],VY[:],EToV)
end

function uniform_tri_mesh(Kx)
        return uniform_tri_mesh(Kx,Kx)
end

function tri_face_vertices()
        return [1,2],[2,3],[3,1]
end


##############################
### Quadrilateral meshes #####
##############################

"""
uniform_quad_mesh(Kx,Ky)

Matlab uniform triangular mesh.

# Examples
```jldoctest
```
"""

function uniform_quad_mesh(Nx,Ny)

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

function uniform_quad_mesh(Kx)
        return uniform_quad_mesh(Kx,Ky)
end

function quad_face_vertices()
        return [1,2],[2,4],[3,4],[1,3] # ordering matters
end

#############################
##### Hexahedral meshes #####
#############################

export uniform_hex_mesh
export hex_face_vertices

"""
uniform_hex_mesh(Kx::Int,Ky::Int)

Matlab uniform hexahedral mesh.

# Examples
```jldoctest
```
"""

function uniform_hex_mesh(Nx,Ny,Nz)
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
                EToV[e,2] = (i+1) + Nxp*j     + Nxp*Nyp*k
                EToV[e,3] = i     + Nxp*(j+1) + Nxp*Nyp*k
                EToV[e,4] = (i+1) + Nxp*(j+1) + Nxp*Nyp*k
                EToV[e,5] = i     + Nxp*j     + Nxp*Nyp*(k+1)
                EToV[e,6] = (i+1) + Nxp*j     + Nxp*Nyp*(k+1)
                EToV[e,7] = i     + Nxp*(j+1) + Nxp*Nyp*(k+1)
                EToV[e,8] = (i+1) + Nxp*(j+1) + Nxp*Nyp*(k+1)
        end
        EToV = @. EToV + 1 # re-index to 1 index

        VX = x[:];
        VY = y[:];
        VZ = z[:];
        return VX[:],VY[:],VZ[:],EToV
end

function uniform_hex_mesh(Kx)
        return uniform_hex_mesh(Kx,Kx,Kx)
end


function hex_face_vertices()
        x1D = LinRange(-1,1,2)
        r, s, t = meshgrid(x1D,x1D,x1D)
        fv1 = map(x->x[1], findall(@. abs(r+1) < 1e-10))
        fv2 = map(x->x[1], findall(@. abs(r-1) < 1e-10))
        fv3 = map(x->x[1], findall(@. abs(s+1) < 1e-10))
        fv4 = map(x->x[1], findall(@. abs(s-1) < 1e-10))
        fv5 = map(x->x[1], findall(@. abs(t+1) < 1e-10))
        fv6 = map(x->x[1], findall(@. abs(t-1) < 1e-10))
        return fv1,fv2,fv3,fv4,fv5,fv6
end
