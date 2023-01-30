using WriteVTK

function quad(x, y)
    return (x + y)^2
end

@testset "VTKWriter test" for N = [1, 2, 3, 4]
    @testset "Write Triangles" begin
        rd = RefElemData(Tri(), N)
        md = MeshData(uniform_mesh(Tri(), 4)..., rd)
        interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) /rd.VDM
        pdata = [quad.(vec(interpolate*md.x), vec(interpolate*md.y))]
        filename = "tri_md_" * string(N)
        check = filename * ".vtu"
        #Todo: Can we implement a better check?
        vtu_name = Meshdata_to_vtk(md, rd, pdata, ["(x+y)^2"], filename, true)
        @test vtu_name[1] == check
    end
end

@testset "Type to vtk" begin
    @test type_to_vtk(Tri()) == VTKCellTypes.VTK_LAGRANGE_TRIANGLE
end

@testset "VTK-Node" for order=[-1, 0, 1]
    @testset "VTK order" begin
        tri_sud_vertices = [-1.0 1.0 -1.0; -1.0 -1.0 1.0]
        if order == -1
            @test triangle_vtk_order(tri_sud_vertices, order, 2) == nothing
        elseif order == 0
            @test triangle_vtk_order(tri_sud_vertices, order, 2) == [-1.0; -1.0]
        else
            @test triangle_vtk_order(tri_sud_vertices, order, 2) == tri_sud_vertices
        end
    end
end