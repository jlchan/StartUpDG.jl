
@testset "Nonconforming MeshData" begin
    rd = RefElemData(Quad(), N=3)
    md = MeshData(NonConformingQuadMeshExample(), rd)

    (; x, y) = md
    u = @. sin(pi * x) * sin(pi * y)

    # interpolate to faces
    num_total_faces = num_faces(rd.element_type) * md.num_elements
    uf = reshape(rd.Vf * u, :, num_total_faces)

    # interpolate faces to mortars (`uf` denotes mortar faces for `NonConformingMesh` types)
    (; nonconforming_faces, mortar_interpolation_matrix) = md.mesh_type
    u_mortar = reshape(mortar_interpolation_matrix * uf[:, nonconforming_faces], :, 
                       num_mortars_per_face(rd) * length(nonconforming_faces))
    
    # construct interior (uM = u⁻ "minus") values and exterior (uP = u⁺ "plus") values
    uM = hcat(uf, u_mortar)
    uP = uM[md.mapP]

    # check that the jumps are zero for continuous u with zero boundary values
    @test maximum(abs.(uP - uM)) < 100 * eps()
end