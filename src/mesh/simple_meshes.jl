"""
    uniform_mesh(elem::Line,Kx)
    uniform_mesh(elem::Tri,Kx,Ky)
    uniform_mesh(elem::Quad,Kx,Ky)
    uniform_mesh(elem::Hex,Kx,Ky,Kz)
    uniform_mesh(elem, K)

Uniform `Kx` (by `Ky` by `Kz`) mesh on ``[-1,1]^d``, where `d` is the spatial dimension.
Returns `(VX,VY,VZ)`, `EToV`. When only one `K` is specified, it assumes a uniform mesh with
`K` elements in each coordinate direction.

`K` can also be specified using a keyword argument `K1D`, e.g., `uniform_mesh(elem; K1D = 16)`.
"""
function uniform_mesh(elem::Line, K1D)
    VX = collect(LinRange(-1, 1, K1D + 1))
    EToV = transpose(reshape(sort([1:K1D; 2:K1D+1]), 2, K1D))
    return (VX,), Matrix(EToV)
end

function uniform_mesh(elem::Tri, Kx, Ky)

    (VY, VX) = meshgrid(LinRange(-1, 1, Ky + 1), LinRange(-1, 1, Kx + 1))
    sk = 1
    EToV = zeros(Int, 2 * Kx * Ky, 3)
    for ey = 1:Ky
        for ex = 1:Kx
            id(ex, ey) = ex + (ey - 1) * (Kx + 1) # index function
            id1 = id(ex, ey)
            id2 = id(ex + 1, ey)
            id3 = id(ex + 1, ey + 1)
            id4 = id(ex, ey + 1)
            EToV[2*sk-1, :] = [id1 id2 id3]
            EToV[2*sk, :] = [id3 id4 id1]
            sk += 1
        end
    end
    return (VX[:], VY[:]), EToV
end

uniform_mesh(elem::Union{Tri,Quad}, Kx) = uniform_mesh(elem, Kx, Kx)

function uniform_mesh(elem::Quad, Nx, Ny)

    Nxp = Nx + 1
    Nyp = Ny + 1
    Nv = convert(Int, Nxp * Nyp)
    K = convert(Int, Nx * Ny)

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    x, y = meshgrid(x1D, y1D)
    I, J = meshgrid(collect(1:Nxp), collect(1:Nyp))
    inds = @. (I - 1) * Ny + (J + I - 1)
    EToV = zeros(Int, K, 4)
    k = 1
    for i = 1:Ny
        for j = 1:Nx
            EToV[k, :] = [inds[i, j] inds[i, j+1] inds[i+1, j] inds[i+1, j+1]]
            k += 1
        end
    end

    VX = x[:]
    VY = y[:]

    return (VX[:], VY[:]), EToV
end


# hack to work with old vertex ordering
function mymeshgrid(x1D, y1D, z1D)
    Np = length(x1D) * length(y1D) * length(z1D)
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
    return x, y, z
end

function uniform_mesh(elem::Hex, Nx, Ny, Nz)
    Nxp = Nx + 1
    Nyp = Ny + 1
    Nzp = Nz + 1
    Nv = convert(Int, Nxp * Nyp * Nzp)
    K = convert(Int, Nx * Ny * Nz)

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    z1D = LinRange(-1, 1, Nzp)

    x, y, z = mymeshgrid(x1D, y1D, z1D)

    EToV = zeros(Int, K, 8)
    k = 1
    for e = 1:K
        em = e - 1
        k = div(em, (Nx * Ny))
        j = div(em - k * Nx * Ny, Nx)
        i = em % Nx

        EToV[e, 1] = i + Nxp * j + Nxp * Nyp * k
        EToV[e, 2] = i + Nxp * (j + 1) + Nxp * Nyp * k
        EToV[e, 3] = (i + 1) + Nxp * j + Nxp * Nyp * k
        EToV[e, 4] = (i + 1) + Nxp * (j + 1) + Nxp * Nyp * k
        EToV[e, 5] = i + Nxp * j + Nxp * Nyp * (k + 1)
        EToV[e, 6] = i + Nxp * (j + 1) + Nxp * Nyp * (k + 1)
        EToV[e, 7] = (i + 1) + Nxp * j + Nxp * Nyp * (k + 1)
        EToV[e, 8] = (i + 1) + Nxp * (j + 1) + Nxp * Nyp * (k + 1)
    end
    EToV = @. EToV + 1 # re-index to 1 index

    return (x[:], y[:], z[:]), EToV
end


function uniform_mesh(::Wedge, Kx, Ky, Kz)
    (VY, VX, VZ) = meshgrid(LinRange(-1, 1, Ky + 1), LinRange(-1, 1, Kx + 1), LinRange(-1, 1, Kz + 1))
    sk = 1
    shift = (Kx + 1) * (Ky + 1)
    id(ex, ey, ez) = ex + (ey - 1) * (Kx + 1)  + (ez - 1) * (shift) # index function
    EToV = zeros(Int, 2 * Kx * Ky * Kz, 6)
    for ey = 1:Ky
        for ex = 1:Kx
            for ez = 1:Kz
                id1 = id(ex, ey, ez)
                id2 = id(ex + 1, ey, ez)
                id3 = id(ex, ey + 1, ez)
                id4 = id(ex + 1, ey + 1, ez)
                id5 = id(ex, ey, ez + 1)
                id6 = id(ex + 1, ey, ez + 1)
                id7 = id(ex, ey + 1, ez + 1)
                id8 = id(ex + 1, ey + 1, ez + 1)
                EToV[2 * sk - 1, :] = [id1 id2 id4 id5 id6 id8] 
                EToV[2 * sk, :]     = [id1 id4 id3 id5 id8 id7]   
                sk += 1
            end
        end
    end
    return vec.((VX, VY, VZ)), EToV
end

# split each cube into 6 pyramids. Pyramid 1: 
#      _________________
#     /.               /
#    / .              / |  
#   /  .             /  |
#  /___.____________/   |   
# |    .     5      |   |
# |    .     ∘      |   |
# |    ._ _ _ _ _ _ | _ |
# |   /2            |  /4
# |  .              | /
# |/________________|/
# 1                 3

function uniform_mesh(elem::Pyr, Nx, Ny, Nz)

    Nxp = Nx + 1
    Nyp = Ny + 1
    Nzp = Nz + 1
    Nv = Nxp * Nyp * Nzp
    K = Nx * Ny * Nz

    id(i, j, k) = i + Nxp * j + Nxp * Nyp * k

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    z1D = LinRange(-1, 1, Nzp)

    x, y, z = vec.(mymeshgrid(x1D, y1D, z1D))

    # 6 pyramids per cube    
    EToV = zeros(Int, 6 * K, 5)
    ee = 1
    for e = 1:K
        em = e - 1
        k = div(em, (Nx * Ny))
        j = div(em - k * Nx * Ny, Nx)
        i = em % Nx

        # add a new center vertex for each element
        ids = [id(i + ii, j + jj, k + kk) for ii in 0:1 for jj in 0:1 for kk in 0:1] .+ 1
        vx, vy, vz = map(x -> mean(view(x, ids)), (x, y, z))
        x = vcat(x, vx)
        y = vcat(y, vy)
        z = vcat(z, vz)
        centroid_id = length(x) - 1 # zero indexing 

        # pyramid 1: bottom face
        EToV[ee, 1] = id(i, j, k)
        EToV[ee, 2] = id(i, j + 1, k)
        EToV[ee, 3] = id(i + 1, j, k)
        EToV[ee, 4] = id(i + 1, j + 1, k)
        EToV[ee, 5] = centroid_id # last entry is the cell centroid vertex

        # pyramid 2: +r face
        EToV[ee+1, 1] = id(i + 1, j, k + 1)
        EToV[ee+1, 2] = id(i + 1, j, k)
        EToV[ee+1, 3] = id(i + 1, j + 1, k + 1)
        EToV[ee+1, 4] = id(i + 1, j + 1, k)
        EToV[ee+1, 5] = centroid_id # last entry is the cell centroid vertex

        # pyramid 3: +s face
        EToV[ee+2, 1] = id(i, j + 1, k + 1)
        EToV[ee+2, 2] = id(i + 1, j + 1, k + 1)
        EToV[ee+2, 3] = id(i, j + 1, k)
        EToV[ee+2, 4] = id(i + 1, j + 1, k)
        EToV[ee+2, 5] = centroid_id # last entry is the cell centroid vertex

        # pyramid 4: -r face
        EToV[ee+3, 1] = id(i, j, k + 1)
        EToV[ee+3, 2] = id(i, j + 1, k + 1)
        EToV[ee+3, 3] = id(i, j, k)
        EToV[ee+3, 4] = id(i, j + 1, k)
        EToV[ee+3, 5] = centroid_id # last entry is the cell centroid vertex

        # pyramid 5: -s face
        EToV[ee+4, 1] = id(i, j, k)
        EToV[ee+4, 2] = id(i + 1, j, k)
        EToV[ee+4, 3] = id(i, j, k + 1)
        EToV[ee+4, 4] = id(i + 1, j, k + 1)
        EToV[ee+4, 5] = centroid_id # last entry is the cell centroid vertex

        # face 6: top face
        EToV[ee+5, 1] = id(i, j, k + 1)
        EToV[ee+5, 2] = id(i + 1, j, k + 1)
        EToV[ee+5, 3] = id(i, j + 1, k + 1)
        EToV[ee+5, 4] = id(i + 1, j + 1, k + 1)
        EToV[ee+5, 5] = centroid_id # last entry is the cell centroid vertex

        ee += 6
    end

    EToV = @. EToV + 1 # re-index to 1 index
    return (x[:], y[:], z[:]), EToV
end

function uniform_mesh(elem::Tet, Nx, Ny, Nz)

    Nxp = Nx + 1
    Nyp = Ny + 1
    Nzp = Nz + 1
    Nv = Nxp * Nyp * Nzp
    K = Nx * Ny * Nz

    id(i, j, k) = i + Nxp * j + Nxp * Nyp * k

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    z1D = LinRange(-1, 1, Nzp)

    x, y, z = mymeshgrid(x1D, y1D, z1D)

    # 6 tets per cube
    EToV = zeros(Int, 6 * K, 4)
    ee = 1
    for e = 1:K
        em = e - 1
        k = div(em, (Nx * Ny))
        j = div(em - k * Nx * Ny, Nx)
        i = em % Nx

        EToV[ee, 1] = id(i, j, k)
        EToV[ee, 2] = id(i + 1, j, k)
        EToV[ee, 3] = id(i + 1, j + 1, k)
        EToV[ee, 4] = id(i, j, k + 1)

        EToV[ee+1, 1] = id(i + 1, j, k)
        EToV[ee+1, 2] = id(i + 1, j, k + 1)
        EToV[ee+1, 3] = id(i + 1, j + 1, k)
        EToV[ee+1, 4] = id(i, j, k + 1)

        EToV[ee+2, 1] = id(i, j, k + 1)
        EToV[ee+2, 2] = id(i + 1, j + 1, k + 1)
        EToV[ee+2, 3] = id(i + 1, j, k + 1)
        EToV[ee+2, 4] = id(i + 1, j + 1, k)

        EToV[ee+3, 1] = id(i, j + 1, k + 1)
        EToV[ee+3, 2] = id(i + 1, j + 1, k + 1)
        EToV[ee+3, 3] = id(i, j, k + 1)
        EToV[ee+3, 4] = id(i + 1, j + 1, k)

        EToV[ee+4, 1] = id(i, j + 1, k)
        EToV[ee+4, 2] = id(i, j + 1, k + 1)
        EToV[ee+4, 3] = id(i, j, k + 1)
        EToV[ee+4, 4] = id(i + 1, j + 1, k)

        EToV[ee+5, 1] = id(i, j, k)
        EToV[ee+5, 2] = id(i + 1, j + 1, k)
        EToV[ee+5, 3] = id(i, j + 1, k)
        EToV[ee+5, 4] = id(i, j, k + 1)

        ee += 6
    end

    EToV = @. EToV + 1 # re-index to 1 index
    return (x[:], y[:], z[:]), EToV
end

uniform_mesh(elem::Union{Hex, Wedge, Pyr, Tet}, Kx) = uniform_mesh(elem, Kx, Kx, Kx)

# keyword argument version
uniform_mesh(elem::AbstractElemShape; K1D) = uniform_mesh(elem, K1D)
