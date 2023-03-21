"""
!!! warning "Experimental implementation"
    This is an experimental feature and may change in future releases.

This is a proof of concept implementation of a non-conforming mesh in StartUpDG.jl. 
The intended usage is as follows:

```julia
rd = RefElemData(Quad(), N=7)
md = MeshData(NonConformingQuadMeshExample(), rd)

(; x, y ) = md
u = @. sin(pi * x) * sin(pi * y)

# interpolate to faces
num_total_faces = num_faces(rd.element_type) * md.num_elements
u_face = reshape(rd.Vf * u, :, num_total_faces)

# interpolate faces to mortars (`uf` denotes mortar faces for `NonConformingMesh` types)
(; conforming_faces, non_conforming_faces, mortar_interpolation_matrix ) = md.mesh_type
u_mortar = similar(md.xf)
view(u_mortar, :, 1:length(conforming_faces)) .= view(u_face, :, conforming_faces)

# interpolate to non-conforming faces, which are stored after the conforming faces
for (i, f) in enumerate(non_conforming_faces)
    mortar_face_ids = (1:num_mortars_per_face(rd)) .+ (i-1) * num_mortars_per_face(rd) .+ length(conforming_faces)
    u_mortar[:, mortar_face_ids] .= reshape(mortar_interpolation_matrix * u_face[:, f], :, num_mortars_per_face(rd))
end

# get exterior values
uP = u_mortar[md.mapP]
```
The `mortar_projection_matrix` similarly maps values from 2 mortar faces back to values on the 
original non-conforming face. These can be used to create DG solvers on non-conforming meshes.

"""
struct NonConformingMesh{TV, TE, CF, NCF, I, P}
    VXYZ::TV
    EToV::TE
    conforming_faces::CF
    non_conforming_faces::NCF
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
    non_conforming_faces::T4
end 

function NonConformingQuadMeshExample()
    VX = [-1, -1,  1,  1, 1,  2, 2, 2]
    VY = [-1,  1, -1,  1, 0, -1, 0, 1]
    EToV = [1 3 2 4; 3 6 5 7; 5 7 4 8]
    FToF = [1, 2, 3, 12, 5, 6, 7, 13, 9, 7, 11, 4, 8] # FToF[mortar_face] = exterior_mortar_face
    non_conforming_faces = [2]
    return NonConformingQuadMeshExample((VX, VY), EToV, FToF, non_conforming_faces)
end

# one non-conforming quad face is split into 2 mortar faces
num_mortars_per_face(rd::RefElemData{2, Quad}) = 2

num_elements(md::MeshData{Dim, <:NonConformingMesh}) where {Dim} = size(getproperty(md.mesh_type, :EToV), 1)

function MeshData(mesh::NonConformingQuadMeshExample, rd::RefElemData{2, Quad})

    (VX, VY) = mesh.VXY
    EToV = mesh.EToV
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
    
    #Construct global coordinates
    (; V1 ) = rd
    x = V1 * VX[transpose(EToV)]
    y = V1 * VY[transpose(EToV)]

    #Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    (; Vf ) = rd
    xf = Vf * x
    yf = Vf * y

    # TODO: fix hard-coding of connectivity
    FToF = mesh.FToF
    non_conforming_faces = mesh.non_conforming_faces
    num_total_faces = num_faces(rd.element_type) * num_elements
    conforming_faces = setdiff(1:num_total_faces, non_conforming_faces)

    # one non-conforming face is divided into 2 in 2D
    num_mortar_faces = num_mortars_per_face(rd) * length(non_conforming_faces) + length(conforming_faces)
   
    # copy over conforming face data
    x_mortar = similar(xf, (num_face_points, num_mortar_faces))
    y_mortar = similar(yf, (num_face_points, num_mortar_faces))
    xf, yf = reshape.((xf, yf), num_face_points, num_total_faces)
    view(x_mortar, :, 1:length(conforming_faces)) .= view(xf, :, conforming_faces)
    view(y_mortar, :, 1:length(conforming_faces)) .= view(yf, :, conforming_faces)

    # interpolate to non-conforming faces
    for (i, f) in enumerate(non_conforming_faces)
        x_interpolated = reshape(mortar_interpolation_matrix * xf[:, f], size(xf, 1), 2)
        y_interpolated = reshape(mortar_interpolation_matrix * yf[:, f], size(xf, 1), 2)        
        
        mortar_face_ids = (1:num_mortars_per_face(rd)) .+ (i-1) * num_mortars_per_face(rd) .+ length(conforming_faces)
        view(x_mortar, :, mortar_face_ids) .= x_interpolated
        view(y_mortar, :, mortar_face_ids) .= y_interpolated
    end

    # compute node maps between mortars
    mapM, mapP, mapB = build_node_maps(FToF, (x_mortar, y_mortar))

    #Compute geometric factors and surface normals
    (; Dr, Ds ) = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)

    (; Vq, wq ) = rd
    xq, yq = (x -> Vq * x).((x, y))
    wJq = diagm(wq) * (Vq * J)

    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)

    is_periodic = (false, false)

    mesh_type = NonConformingMesh(tuple(VX, VY), EToV, conforming_faces, non_conforming_faces, 
                                  mortar_interpolation_matrix, mortar_projection_matrix)

    return MeshData(mesh_type, FToF,
                    tuple(x, y), tuple(x_mortar, y_mortar), tuple(xq, yq), wJq,
                    mapM, mapP, mapB,
                    SMatrix{2, 2}(tuple(rxJ, ryJ, sxJ, syJ)), J,
                    tuple(nxJ, nyJ), Jf,
                    is_periodic)
end
