
@testset "Nonconforming MeshData" begin
    rd = RefElemData(Quad(), N=3)
    md = MeshData(NonConformingQuadMeshExample(), rd)

    (; x, y) = md
    u = @. sin(pi * x) * sin(pi * y)

    # interpolate to faces
    num_total_faces = num_faces(rd.element_type) * md.num_elements
    u_face = reshape(rd.Vf * u, :, num_total_faces)

    # interpolate faces to mortars (`uf` denotes mortar faces for `NonConformingMesh` types)
    (; conforming_faces, non_conforming_faces, mortar_interpolation_matrix) = md.mesh_type

    u_mortar = similar(md.xf)
    view(u_mortar, :, 1:length(conforming_faces)) .= view(u_face, :, conforming_faces)    
    # interpolate to non-conforming faces
    for (i, f) in enumerate(non_conforming_faces)
        mortar_face_ids = (1:num_mortars_per_face(rd)) .+ (i-1) * num_mortars_per_face(rd) .+ length(conforming_faces)
        u_mortar[:, mortar_face_ids] .= reshape(mortar_interpolation_matrix * u_face[:, f], :, num_mortars_per_face(rd))
    end

    # get exterior values
    uP = u_mortar[md.mapP]

    # check that the jumps are zero for continuous u with zero boundary values
    @test maximum(abs.(uP - u_mortar)) < 100 * eps()
end