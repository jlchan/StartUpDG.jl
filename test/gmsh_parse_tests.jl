# gmsh test files are located in StartUpDG/test/testset_mesh/
# malpasset.msh was previously a 2.2 file. Exported version 
# of 4.1 has been added for testing

@testset "file:$file" for file in readdir("testset_mesh",join=true) 
    io = open("stderr.txt","w")
    redirect_stderr(io)
    f = open(file)
    lines = readlines(f)
    format_line = StartUpDG.findline("\$MeshFormat",lines)+1
    version,_,dataSize = split(lines[format_line])
    gmsh_version = parse(Float64,version)

    if gmsh_version == 2.2
        VXY, EToV = readGmsh2D(file);
    elseif gmsh_version == 4.1
        VXY, EToV = readGmsh2D_v4(file);
        VXY, EToV, grouping = readGmsh2D_v4(file,true);
    end
    @test 1==1
end

@testset "gmsh version 4.1 file with one data grouping" begin
    io = open("stderr.txt","w")
    redirect_stderr(io)
    file = "testset_mesh/one_group_v4.msh"
    if isfile(file)
        VXY_1,EToV_1 = readGmsh2D_v4(file)
        VXY_2,EToV_2,group_2 = readGmsh2D_v4(file,true) 
        @test VXY_1 == VXY_2 
        @test EToV_1 == EToV_2 
        f = open(file)
        lines = readlines(f)
        num_elements = StartUpDG.getNumElements(lines)
        @test length(unique(group_2))==1
    else
        @info "file for this test is missing"
    end
    close(io)
    rm("stderr.txt")
end;

@testset "gmsh version 4.1 file with no grouping data" begin
    io = open("stderr.txt","w")
    redirect_stderr(io)
    file = "testset_mesh/no_group_v4.msh"
    if isfile(file)
        VXY_1,EToV_1 = readGmsh2D_v4(file)
        VXY_2,EToV_2,group_2 = readGmsh2D_v4(file,true) 
        @test VXY_1 == VXY_2 
        @test EToV_1 == EToV_2 
        f = open(file)
        lines = readlines(f)
        num_elements = StartUpDG.getNumElements(lines)
        @test group_2 == zeros(Int,num_elements)
    else
        @info "file for this test is missing"
    end
    close(io)
    rm("stderr.txt")
end;

@testset "Compare output between v2.2 and v4.1" begin
@testset "file:$file" for file in ["mesh_no_pert","pert_mesh","malpasset"]
    if isfile("testset_mesh/$file.msh") && isfile("testset_mesh/$file"*"_v4.msh")
        VXY_v2, EToV_v2 = readGmsh2D("testset_mesh/$file.msh");
        VXY_v4, EToV_v4 = readGmsh2D_v4("testset_mesh/$file" * "_v4.msh");
        @test VXY_v2 == VXY_v4
        @test EToV_v2 == EToV_v4 
    end
end
end