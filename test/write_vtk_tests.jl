@testset "VTK tests" begin

    function quad(x, y)
        return (x + y)^2
    end

    # test expected output of `vtk_order`
    deg_one_order(::Tri) = permutedims(hcat(nodes(Tri(), 1)...))
    deg_zero_order(::Tri) = [-1.0; -1.0]
    deg_one_order(::Quad) = [-1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0]
    deg_one_order(::Wedge) = [-1.0 1.0 -1.0 -1.0 1.0 -1.0; -1.0 -1.0 1.0 -1.0 -1.0 1.0; -1.0 -1.0 -1.0 1.0 1.0 1.0]
    deg_zero_order(elem::Union{Quad, Wedge}) = deg_one_order(elem)


    @testset "VTKWriter test for $elem" for elem in [Tri(), Quad(), Wedge()]
        N = 3 # test only N=3 for CI time
        @testset "Write Mesh" begin
            rd = RefElemData(elem, N)
            md = MeshData(uniform_mesh(elem, 2)..., rd)
            interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
            pdata = [quad.(interpolate * md.x, interpolate * md.y)]
            filename = replace(string(elem), "()" => "") * "_" * string(N)
            check = filename * ".vtu"
            # Todo: Can we implement a better check?
            vtu_name = MeshData_to_vtk(md, rd, pdata, ["(x+y)^2"], filename, true)
            @test vtu_name[1] == check
            rm(check) # remove created file after test is done
        end

        @testset "VTK-Node degree $order" for order=[0, 1]
            @testset "VTK order" begin
                @test_throws AssertionError StartUpDG.vtk_order(elem, -1)
                if order == 0
                    @test StartUpDG.vtk_order(elem, order) ≈ deg_zero_order(elem)
                else
                    @test StartUpDG.vtk_order(elem, order) ≈ deg_one_order(elem)
                end
            end
        end
    end

end # VTK tests