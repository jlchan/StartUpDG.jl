# gmsh test files are located in StartUpDG/test/testset_Gmsh_meshes/
# malpasset.msh was previously a 2.2 file. Exported version 
# of 4.1 has been added for testing

@testset "Gmsh" begin 
    @testset "Gmsh reading" begin
        VXY, EToV = read_Gmsh_2D("testset_Gmsh_meshes/squareCylinder2D.msh")
        @test size(EToV) == (3031, 3)

        # malpasset data taken from 
        # https://github.com/li12242/NDG-FEM/blob/master/Application/SWE/SWE2d/Benchmark/%40Malpasset2d/mesh/triMesh/malpasset.msh
        VXY, EToV = read_Gmsh_2D("testset_Gmsh_meshes/malpasset.msh")
        rd = RefElemData(Tri(), 1)
        md = MeshData(VXY, EToV, rd)
        @test all(md.J .> 0)
    end
    
    @testset "test gmsh element id remapping" begin
        testvec = [16, 15, 16, 17, 18]
        @test StartUpDG.remap_element_grouping(testvec) == [1, 2, 1, 3, 4]
    end

    @testset "$approxType MeshData initialization with gmsh import" for approxType = [Polynomial(), SBP()]
        @testset "2D tri gmsh import verison $version" for version = [2.2, 4.1] 
            file = "pert_mesh"
            tol = 5e7*eps() # higher tolerance due to floating point issues?
            N = 3

            rd = RefElemData(Tri(), approxType, N)
            if version == 2.2
                VXY, EToV = read_Gmsh_2D("testset_Gmsh_meshes/$file"*".msh");
            elseif version == 4.1
                VXY, EToV = read_Gmsh_2D("testset_Gmsh_meshes/$file"*"_v4.msh");
            end
            md = MeshData(VXY, EToV, rd)

            (; wq, Dr, Ds, Vq, Vf, wf  ) = rd
            Nfaces = length(rd.fv)
            (; x, y, xq, yq, xf, yf, K  ) = md
            (; rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ, wJq  ) = md
            (; FToF, mapM, mapP, mapB  ) = md

            @test md.x == md.xyz[1]

            @testset  "check positivity of Jacobian" begin
                @test all(J .> 0)
            end

            @testset  "check differentiation" begin
                u = @. x^2 + 2 * x * y - y^2
                dudx_exact = @. 2*x + 2*y
                dudy_exact = @. 2*x - 2*y
                dudr,duds = (D->D*u).((Dr, Ds))
                dudx = @. (rxJ * dudr + sxJ * duds) / J
                dudy = @. (ryJ * dudr + syJ * duds) / J
                @test dudx ≈ dudx_exact
                @test dudy ≈ dudy_exact
            end

            @testset "check volume integration" begin
                @test Vq*x ≈ xq
                @test Vq*y ≈ yq
                @test diagm(wq) * (Vq * J) ≈ wJq
                #@test abs(sum(xq .* wJq)) < tol skip=true 
                #@test abs(sum(yq .* wJq)) < tol skip=true
            end

            @testset "check surface integration begin" begin
                @test Vf * x ≈ xf
                @test Vf * y ≈ yf
                @test abs(sum(wf .* nxJ)) < tol 
                @test abs(sum(wf .* nyJ)) < tol 
                #@test sum(@. wf * nxJ * (1 + xf) / 2) ≈ 2.0 skip=true # check sign of normals
            end

            @testset "check connectivity and boundary maps" begin
                u = @. (1-x) * (1+x) * (1-y) * (1+y)
                uf = Vf * u
                @test uf ≈ uf[mapP]        
                #@test norm(uf[mapB]) < tol skip=true
            end

            @testset "check periodic node connectivity maps" begin
                md = make_periodic(md, (true, true))
                (; mapP  ) = md
                u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y))
                uf = Vf * u
                #@test uf ≈ uf[mapP]
            end

            @testset "check MeshData struct copying" begin
                xyz = (x->x .+ 1).(md.xyz) # affine shift
                md2 = MeshData(rd, md, xyz...)
                @test sum(norm.(md2.rstxyzJ .- md.rstxyzJ)) < tol
                @test sum(norm.(md2.nxyzJ .- md.nxyzJ)) < tol
                @test all(md2.xyzf .≈ (x->x .+ 1).(md.xyzf))
            end
        end
    end
    
    @testset "gmsh version 4.1 file with one data grouping" begin
        file = "testset_Gmsh_meshes/one_group_v4.msh"
        if isfile(file)
            VXY_1, EToV_1 = read_Gmsh_2D(file)
            VXY_2, EToV_2, group_2 = read_Gmsh_2D(file, true) 
            @test VXY_1 == VXY_2 
            @test EToV_1 == EToV_2 
            f = open(file)
            lines = readlines(f)
            num_elements = StartUpDG.get_num_elements(lines)
            @test length(unique(group_2))==1
	    close(f)
        else
            @warn "file for this test is missing"
        end
    end;

    @testset "gmsh version 4.1 file with no grouping data" begin
        # test promps for grouping data from the file. 
        # suppose such grouping data however is not in the file
        file = "testset_Gmsh_meshes/no_group_v4.msh"
        if isfile(file)
            VXY_1, EToV_1 = read_Gmsh_2D(file)
            VXY_2, EToV_2, group_2 = read_Gmsh_2D(file, true) 
            @test VXY_1 == VXY_2 
            @test EToV_1 == EToV_2 
            f = open(file)
            lines = readlines(f)
            num_elements = StartUpDG.get_num_elements(lines)
            @test group_2 == zeros(Int, num_elements)
	    close(f)
        else
            @warn "file for this test is missing"
        end
    end;

    @testset "Compare output between v2.2 and v4.1" begin
        @testset "file:$file" for file in ["mesh_no_pert","pert_mesh","malpasset"]
            if isfile("testset_Gmsh_meshes/$file.msh") && isfile("testset_Gmsh_meshes/$file"*"_v4.msh")
                VXY_v2, EToV_v2 = read_Gmsh_2D("testset_Gmsh_meshes/$file.msh");
                VXY_v4, EToV_v4 = read_Gmsh_2D("testset_Gmsh_meshes/$file" * "_v4.msh");
                @test VXY_v2 == VXY_v4
                @test EToV_v2 == EToV_v4 
            end
        end
    end
end
