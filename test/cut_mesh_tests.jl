@testset "Cut meshes" begin
    
    cells_per_dimension = 4
    cells_per_dimension_x, cells_per_dimension_y = cells_per_dimension, cells_per_dimension
    circle = PresetGeometries.Circle(R=0.33, x0=0, y0=0)

    rd = RefElemData(Quad(), N=3)
    md = MeshData(rd, (circle, ), cells_per_dimension_x, cells_per_dimension_y)

    # check the volume of the domain
    A = 4 - pi * .33^2
    @test sum(md.wJq) ≈ A

    # check the length of the boundary of the domain
    face_weights = reshape(rd.wf, :, num_faces(rd.element_type))[:, 1]
    wJf = vec(Diagonal(face_weights) * reshape(md.Jf, length(face_weights), :))
    @test sum(wJf[md.mapB]) ≈ (8 + 2 * pi * .33)

    # check continuity of a function that's in the global polynomial space
    @unpack physical_frame_elements = md.mesh_type
    @unpack x, y = md
    u = @. x^rd.N - x * y^(rd.N-1) - x^(rd.N-1) * y + y^rd.N
    uf = similar(md.xf)
    uf.cartesian .= rd.Vf * u.cartesian
    for e in 1:size(md.x.cut, 2)
        ids = md.mesh_type.cut_face_nodes[e]
        VDM = vandermonde(physical_frame_elements[e], rd.N, md.x.cut[:, e], md.y.cut[:, e])
        Vf = vandermonde(physical_frame_elements[e], rd.N, md.xf.cut[ids], md.yf.cut[ids]) / VDM
        uf.cut[ids] .= Vf * u.cut[:, e]
    end
    @test all(uf .≈ vec(uf[md.mapP]))

    # TODO: add tests on taking derivatives 
end