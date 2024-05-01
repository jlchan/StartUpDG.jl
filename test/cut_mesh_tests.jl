using StartUpDG: PathIntersections

@testset "Cut meshes ($(typeof(quadrature_type)))" for quadrature_type = [Subtriangulation(), MomentFitting()]
    @testset "Correctness" begin
        
        cells_per_dimension = 4
        cells_per_dimension_x, cells_per_dimension_y = cells_per_dimension, cells_per_dimension
        circle = PresetGeometries.Circle(R=0.33, x0=0, y0=0)

        rd = RefElemData(Quad(), N=3)
        md = MeshData(rd, (circle, ), 
                    cells_per_dimension_x, cells_per_dimension_y,
                    quadrature_type; 
                    precompute_operators=true)

        @test_throws ErrorException("Face index f = 5 > 4; too large.") StartUpDG.neighbor_across_face(5, nothing, nothing)

        @test StartUpDG.num_cartesian_elements(md) + StartUpDG.num_cut_elements(md) == md.num_elements

        @test (@capture_out Base.show(stdout, MIME"text/plain"(), md)) == "Cut-cell MeshData of dimension 2 with 16 elements (12 Cartesian, 4 cut)"

        # test constructor with only one "cells_per_dimension" argument
        @test_nowarn MeshData(rd, (circle, ), cells_per_dimension_x, quadrature_type)

        # check the volume of the domain
        A = 4 - pi * .33^2
        @test abs(sum(md.wJq) - A) < 7e-5 # should be 6.417e-5

        # check the length of the boundary of the domain
        wJf = md.mesh_type.cut_cell_data.wJf
        @test abs(sum(wJf[md.mapB]) - (8 + 2 * pi * .33)) < 8e-5 # should be 7.66198e-5

        # check continuity of a function that's in the global polynomial space
        (; physical_frame_elements) = md.mesh_type
        (; x, y) = md
        u = @. x^rd.N - x * y^(rd.N-1) - x^(rd.N-1) * y + y^rd.N
        uf = similar(md.xf)
        uf.cartesian .= rd.Vf * u.cartesian
        for e in 1:size(md.x.cut, 2)
            ids = md.mesh_type.cut_face_nodes[e]
            Vf = md.mesh_type.cut_cell_operators.face_interpolation_matrices[e]
            uf.cut[ids] .= Vf * u.cut[:, e]
        end
        @test all(uf .≈ vec(uf[md.mapP]))

        dudx_exact = @. rd.N * x^(rd.N-1) - y^(rd.N-1) - (rd.N-1) * x^(rd.N-2) * y
        dudy_exact = @. -(rd.N-1) * x * y^(rd.N-2) - x^(rd.N-1) + rd.N * y^(rd.N-1)
        (; physical_frame_elements, cut_face_nodes) = md.mesh_type
        dudx, dudy = similar(md.x), similar(md.x)     
        dudx.cartesian .= (md.rxJ.cartesian .* (rd.Dr * u.cartesian)) ./ md.J.cartesian
        dudy.cartesian .= (md.syJ.cartesian .* (rd.Ds * u.cartesian)) ./ md.J.cartesian
        for (e, elem) in enumerate(physical_frame_elements)
            VDM = vandermonde(elem, rd.N, x.cut[:, e], y.cut[:, e])
            Vq, Vxq, Vyq = map(A -> A / VDM, basis(elem, rd.N, md.xq.cut[:,e], md.yq.cut[:, e]))
        
            M  = Vq' * diagm(md.wJq.cut[:, e]) * Vq
            Qx = Vq' * diagm(md.wJq.cut[:, e]) * Vxq
            Qy = Vq' * diagm(md.wJq.cut[:, e]) * Vyq   
            Dx, Dy = M \ Qx, M \ Qy
            # LIFT = M \ (Vf' * diagm(wJf))

            # TODO: add interface flux terms into test
            dudx.cut[:, e] .= Dx * u.cut[:,e] # (md.rxJ.cut[:,e] .* (Dr * u.cut[:,e])) 
            dudy.cut[:, e] .= Dy * u.cut[:,e] # (md.rxJ.cut[:,e] .* (Dr * u.cut[:,e])) 
        end

        @test dudx ≈ dudx_exact
        @test dudy ≈ dudy_exact

        # test normals are unit and non-zero
        @test all(@. md.nx^2 + md.ny^2 ≈ 1)

        # test creation of equispaced nodes on cut cells
        x, y = equi_nodes(physical_frame_elements[1], circle, 10)
        # shouldn't have more points than equispaced points on a quad
        @test 0 < length(x) <= length(first(equi_nodes(Quad(), 10))) 
        # no points should be contained in the circle
        @test !all(PathIntersections.is_contained.(circle, zip(x, y))) 

    end
        
    @testset "State redistribution" begin     
        # test state redistribution 
        cells_per_dimension = 4
        circle = PresetGeometries.Circle(R=0.6, x0=0, y0=0)
        rd = RefElemData(Quad(), N=4)
        md = MeshData(rd, (circle, ), cells_per_dimension, quadrature_type)

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
end

#######################################
#        test weak SBP property       #
#######################################
@testset "Cut mesh weak SBP property" begin

    cells_per_dimension = 2
    circle = PresetGeometries.Circle(R=0.66, x0=.1, y0=0)

    rd = RefElemData(Quad(), N=3)

    md = MeshData(rd, (circle, ), cells_per_dimension; 
                  precompute_operators=true)

    mt = md.mesh_type
    (; wJf) = md.mesh_type.cut_cell_data
    wf = wJf ./ md.Jf

    for (e, elem) in enumerate(mt.physical_frame_elements)
        xq, yq, wq = md.xq.cut[:,e], md.yq.cut[:,e], md.wJq.cut[:,e]

        Vq, Vxq, Vyq = basis(elem, rd.N, xq, yq)
        face_ids = mt.cut_face_nodes[e]
        Vf = vandermonde(elem, rd.N, md.xf.cut[face_ids], md.yf.cut[face_ids]) 
        Qx, Qy = Vq' * diagm(wq) * Vxq, Vq' * diagm(wq) * Vyq

        Bx = Diagonal(wf.cut[face_ids] .* md.nxJ.cut[face_ids])
        By = Diagonal(wf.cut[face_ids] .* md.nyJ.cut[face_ids])

        # represent constant vector in basis
        ee = Vq \ ones(size(Vq, 1))
        @test norm(ee' * Qx - ee' * Vf' * Bx * Vf) < 100 * eps()
        @test norm(ee' * Qy - ee' * Vf' * By * Vf) < 100 * eps()
    end
end