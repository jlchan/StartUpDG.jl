@testset "VTK tests" begin

    # test expected output of `vtk_order`
    deg_one_order(::Tri) = permutedims(hcat(nodes(Tri(), 1)...))
    deg_zero_order(::Tri) = [-1.0; -1.0]
    deg_one_order(::Quad) = [-1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0]
    deg_one_order(::Hex) = [-1.0 -1.0  1.0  1.0 -1.0 -1.0  1.0  1.0;
                            -1.0  1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0;
                            -1.0 -1.0 -1.0 -1.0  1.0  1.0  1.0  1.0]
    deg_one_order(::Wedge) = [-1.0  -1.0   1.0  -1.0  -1.0   1.0
                              -1.0   1.0  -1.0  -1.0   1.0  -1.0
                              -1.0  -1.0  -1.0   1.0   1.0   1.0]
    deg_one_order(::Tet) = [-1.0  1.0 -1.0 -1.0;
                            -1.0 -1.0  1.0 -1.0;
                            -1.0 -1.0 -1.0  1.0]
    deg_zero_order(elem::Union{Quad, Hex, Wedge, Tet}) = deg_one_order(elem)


    @testset "VTKWriter test for $elem" for elem in [Tri(), Quad(), Hex(), Wedge(), Tet()]
        @testset "Write Mesh" begin
            function quad(x, y)
                return (x + y)^2
            end
            rd = RefElemData(elem, 2) # test only N=2 for CI time
            md = MeshData(uniform_mesh(elem, 2)..., rd)
            interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
            pdata = [quad.(interpolate * md.x, interpolate * md.y)]
            filename = replace(string(elem), "()" => "") * "_" * string(rd.N)
            check = filename * ".vtu"
            # Todo: Can we implement a better check?
            vtu_name = MeshData_to_vtk(md, rd, pdata, ["(x+y)^2"], filename, true)
            @test vtu_name[1] == check
            rm(check) # remove created file after test is done

            # ======= check `export_to_vtk` interfaces
            
            # copy of MeshData_to_vtk but with rd, md reversed for consistency
            vtu_name = export_to_vtk(rd, md, pdata, ["(x+y)^2"], filename; 
                                     equi_dist_nodes=false)
            @test vtu_name[1] == check
            rm(check) # remove created file after test is done

            # without filename specified
            vtu_name = export_to_vtk(rd, md, pdata, filename)            
            @test vtu_name[1] == check
            rm(check) # remove created file after test is done

            # with data/dataname specified via Dict
            data = quad.(interpolate * md.x, interpolate * md.y)
            vtu_name = export_to_vtk(rd, md, Dict("(x+y)^2" => data), filename)
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

    @testset "TensorProduct VTKWriter" begin
        @testset "Degree ($tri_N, $line_N)" for tri_N in 2, line_N in 3
            function quad(x, y)
                return (x + y)^2
            end
            line = RefElemData(Line(), line_N)
            tri  = RefElemData(Tri(), tri_N)
            tensor = TensorProductWedge(tri, line)
            rd = RefElemData(Wedge(), tensor)
            md = MeshData(uniform_mesh(Wedge(), 2)..., rd)
            pdata = [quad.(md.x, md.y)]
            filename = "TensorProductWedge" * "_" * string(tri_N) * "_" * string(line_N)
            check = filename * ".vtu"
            # Todo: Can we implement a better check?
            vtu_name = MeshData_to_vtk(md, rd, pdata, ["(x+y)^2"], filename, true, true)
            @test vtu_name[1] == check
            rm(check) # remove created file after test is done
        end
    end

end # VTK tests
