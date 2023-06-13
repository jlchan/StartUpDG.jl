"""
!!! warning "Experimental implementation"
    This is an experimental feature and may change without warning in future releases.

This is a proof of concept implementation of a non-conforming mesh in StartUpDG.jl. 
The intended usage is as follows:

```julia
rd = RefElemData(Quad(), N=7)
md = MeshData(NonConformingQuadMeshExample(), rd)

(; x, y) = md
u = @. sin(pi * x) * sin(pi * y)

# interpolate to faces
num_total_faces = num_faces(rd.element_type) * md.num_elements
uf = reshape(rd.Vf * u, :, num_total_faces)

# interpolate faces to mortars (`uf` denotes mortar faces for `NonConformingMesh` types)
(; nonconforming_faces, mortar_interpolation_matrix) = md.mesh_type

u_mortar = reshape(mortar_interpolation_matrix * uf[:, nonconforming_faces], :, 2 * length(nonconforming_faces))

# construct interior (uM = u⁻ "minus") values and exterior (uP = u⁺ "plus") values
uM = hcat(uf, u_mortar)
uP = uM[md.mapP]
```
The `mortar_projection_matrix` similarly maps values from 2 mortar faces back to values on the 
original non-conforming face. These can be used to create DG solvers on non-conforming meshes.

"""
struct NonConformingMesh{TV, TE, NCF, I, P}
    VXYZ::TV
    EToV::TE
    nonconforming_faces::NCF
    mortar_interpolation_matrix::I
    mortar_projection_matrix::P
end

#  Example non-conforming quad mesh
# 
# Vertex ordering
#       2_____________4 ________8
#    1 |               |        |
#      |               |    3   |
#      |       1      5|________7
#      |               |        |
#      |               |    2   |
#   -1 |1_____________3|________6
#      -1              1       1.5
# 
# Face ordering
#       _______________ ________
#      |       4       |    12  |
#      |               |9  3  10|
#      |1      1      2|____11__|
#      |               |    8   |
#      |               |5   2  6|
#      |_______3_______|____7___|
#                           
# Mortar ordering: first the unsplit faces, then the split faces
#       _______________ ________
#      |       3       |    11  |
#      |             13|8  3   9|
#      |1      1       |____10__|
#      |               |    7   |
#      |             12|4   2  5|
#      |_______2_______|____6___|
#                           
# custom type to trigger this demo
struct NonConformingQuadMeshExample{T1, T2, T3, T4}
    VXY::T1
    EToV::T2
    FToF::T3
    nonconforming_faces::T4
end 

function NonConformingQuadMeshExample()
    VX = [-1, -1,  1,  1, 1,  2, 2, 2]
    VY = [-1,  1, -1,  1, 0, -1, 0, 1]
    EToV = [1 3 2 4; 3 6 5 7; 5 7 4 8]
    #FToF = [1, 2, 3, 12, 5, 6, 7, 13, 9, 7, 11, 4, 8] # FToF[mortar_face] = exterior_mortar_face
    FToF = [1,  2,  3,  4,  13,  6,  7,  11,  14,  10,  8,  12,  5,  9]
    nonconforming_faces = [2]
    return NonConformingQuadMeshExample((VX, VY), EToV, FToF, nonconforming_faces)
end

# one non-conforming quad face is split into 2 mortar faces
num_mortars_per_face(rd::RefElemData) = num_mortars_per_face(rd.element_type)
num_mortars_per_face(::Union{Quad, Tri}) = 2

num_elements(md::MeshData{Dim, <:NonConformingMesh}) where {Dim} = size(getfield(getfield(md, :mesh_type), :EToV), 1)

MeshData(mesh::NonConformingQuadMeshExample, rd::RefElemData{2, Quad}) = 
    MeshData(mesh.VXY, mesh.EToV, mesh.FToF, mesh.nonconforming_faces, rd)

function MeshData(VXY, EToV, FToF, nonconforming_faces, rd::RefElemData{2, Quad})    

    VX, VY = VXY
    num_elements = size(EToV, 1)

    # assume each mortar face is a uniform subdivision
    num_face_points = length(rd.rf) ÷ num_faces(rd.element_type)
    r1D = rd.sf[1:num_face_points]
    w1D = rd.wf[1:num_face_points]
    r_split = vcat(0.5 * (1 .+ r1D) .- 1, 0.5 * (1 .+ r1D)) 
    w_split = 0.5 * vcat(w1D, w1D)

    # maps interpolation to face node onto mortars
    mortar_interpolation_matrix = vandermonde(Line(), rd.N, r_split) / vandermonde(Line(), rd.N, r1D)
    mortar_mass_matrix = mortar_interpolation_matrix' * diagm(w_split) * mortar_interpolation_matrix

    # L2 projection of mortars back to face nodes
    mortar_projection_matrix = mortar_mass_matrix \ (mortar_interpolation_matrix' * diagm(w_split))
    
    # construct element nodal coordinates
    (; V1 ) = rd
    x = V1 * VX[transpose(EToV)]
    y = V1 * VY[transpose(EToV)]

    # construct face coordinates
    (; Vf ) = rd
    xf = Vf * x
    yf = Vf * y

    # reshape xf, yf into arrays of size (num_nodes_per_face x num_faces)
    num_element_faces = num_faces(rd.element_type) * num_elements
    xf, yf = reshape.((xf, yf), num_face_points, num_element_faces)

    x_mortar = reshape(mortar_interpolation_matrix * xf[:, nonconforming_faces], :, 
                       num_mortars_per_face(rd) * length(nonconforming_faces))
    y_mortar = reshape(mortar_interpolation_matrix * yf[:, nonconforming_faces], :, 
                       num_mortars_per_face(rd) * length(nonconforming_faces))                       

    # each nonconforming face is split into 2 mortar faces. we append the mortar faces to the existing element faces. 
    # TODO: should we use LazyArrays.Hcat for this in the future?
    xM = hcat(xf, x_mortar)
    yM = hcat(yf, y_mortar)

    # compute node maps between mortars
    mapM, mapP, mapB = build_node_maps(FToF, (xM, yM))

    #Compute geometric factors and surface normals
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, rd.Drst...)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)

    (; Vq, wq) = rd
    xq, yq = (x -> Vq * x).((x, y))
    wJq = diagm(wq) * (Vq * J)

    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)

    is_periodic = (false, false)

    mesh_type = NonConformingMesh(tuple(VX, VY), EToV, 
                                  nonconforming_faces, 
                                  mortar_interpolation_matrix, 
                                  mortar_projection_matrix)

    return MeshData(mesh_type, FToF,
                    tuple(x, y), tuple(x_mortar, y_mortar), tuple(xq, yq), wJq,
                    mapM, mapP, mapB,
                    SMatrix{2, 2}(tuple(rxJ, ryJ, sxJ, syJ)), J,
                    tuple(nxJ, nyJ), Jf,
                    is_periodic)
end
