@testset "Mesh, timestep, and Triangulate utils" begin
    tol = 5e2*eps()
    @testset "Timestep utils" begin
        rk4a,rk4b,rk4c = ck45()
        @test rk4c[1] ≈ 0.0 && rk4c[end] ≈ 1.0
        @test rk4a[1] ≈ 0.0
        @test all(rk4b .> 0)
    end

    @testset "Gmsh reading" begin
        VXY, EToV = readGmsh2D("squareCylinder2D.msh")
        @test size(EToV)==(3031,3)

        # malpasset data taken from 
        # https://github.com/li12242/NDG-FEM/blob/master/Application/SWE/SWE2d/Benchmark/%40Malpasset2d/mesh/triMesh/malpasset.msh
        VXY, EToV = readGmsh2D("malpasset.msh")
        rd = RefElemData(Tri(), 1)
        md = MeshData(VXY, EToV, rd)
        @test all(md.J .> 0) 
    end

    @testset "Triangulate utils/example meshes" begin
        # test triangulate
        meshIO = triangulate_domain(SquareDomain())
        (VX,VY),EToV = triangulateIO_to_VXYEToV(meshIO)
        rd = RefElemData(Tri(),2)
        md = MeshData((VX,VY),EToV,rd)
        @test size(EToV,1)==md.num_elements==620
        @test length(VX)==length(VY)==338
        @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2,3,4]

        h = .1
        meshIO = triangulate_domain(RectangularDomainWithHole()) # square_hole_domain(h)
        (VX,VY),EToV = triangulateIO_to_VXYEToV(meshIO)
        rd = RefElemData(Tri(),2)
        md = MeshData((VX,VY),EToV,rd)
        @test size(EToV,1)==md.num_elements==598
        @test length(VX)==length(VY)==327
        @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2]
        meshIO2 = refine(meshIO,h)
        @test sort(unique(get_node_boundary_tags(meshIO2,rd,MeshData(triangulateIO_to_VXYEToV(meshIO2)...,rd))))==[0,1,2]

        meshIO = triangulate_domain(Scramjet())
        (VX,VY),EToV = triangulateIO_to_VXYEToV(meshIO)
        rd = RefElemData(Tri(),2)
        md = MeshData((VX,VY),EToV,rd)
        @test size(EToV,1)==md.num_elements==1550
        @test length(VX)==length(VY)==871
        @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2,3]

        boundary_faces = tag_boundary_faces(meshIO, rd, md, Dict(:wall=>1, :inflow=>2, :outflow=>3))
        @test length(boundary_faces[:wall])==173
        @test length(boundary_faces[:inflow])==14
        @test length(boundary_faces[:outflow])==5
    end
end
