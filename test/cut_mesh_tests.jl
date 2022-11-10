using PathIntersections

@testset "Cut meshes" begin
    
    cells_per_dimension = 4
    cells_per_dimension_x, cells_per_dimension_y = cells_per_dimension, cells_per_dimension
    circle = PresetGeometries.Circle(R=0.33, x0=0, y0=0)

    rd = RefElemData(Quad(), N=3)
    md = MeshData(rd, (circle, ), cells_per_dimension_x, cells_per_dimension_y)

    @test_throws ErrorException("Face index f = 5 > 4; too large.") StartUpDG.neighbor_across_face(5, nothing, nothing)

    @test StartUpDG.num_cartesian_elements(md) + StartUpDG.num_cut_elements(md) == md.num_elements

    @test (@capture_out Base.show(stdout, MIME"text/plain"(), md)) == "Cut-cell MeshData of dimension 2 with 16 elements (12 Cartesian, 4 cut)"

    # test constructor with only one "cells_per_dimension" argument
    @test_nowarn MeshData(rd, (circle, ), cells_per_dimension_x)

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

    # test creation of equispaced nodes on cut cells
    x, y = equi_nodes(physical_frame_elements[1], circle, 10)
    # shouldn't have more points than equispaced points on a quad
    @test 0 < length(x) <= length(first(equi_nodes(Quad(), 10))) 
    # no points should be contained in the circle
    @test !all(is_contained.(circle, zip(x, y))) 

    # TODO: add tests on taking derivatives 

end
    
@testset "State redistribution" begin     
    # test state redistribution 
    cells_per_dimension = 4
    circle = PresetGeometries.Circle(R=0.6, x0=0, y0=0)
    rd = RefElemData(Quad(), N=4)
    md = MeshData(rd, (circle, ), cells_per_dimension, cells_per_dimension)

    srd = StateRedistribution(rd, md)
    e = @. 0 * md.x + 1 # constant
    u = @. md.x + md.x^3 .* md.y # degree 4 polynomial
    ecopy, ucopy = copy.((e, u))

    # two ways of applying SRD
    apply!(u, srd)
    srd(e) # in-place application of SRD functor

    # test exactness
    @test norm(e .- ecopy) < 1e3 * eps()
    @test norm(u .- ucopy) < 1e3 * eps()
end