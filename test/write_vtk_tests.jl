function quad(x, y)
    return (x + y)^2
end

# test expected output of `vtk_order`
function deg_one_order(elem)
    if elem == Tri()
        return permutedims(hcat(nodes(Tri(), 1)...))
    elseif elem == Quad()
        return [-1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0]
    end
end

function deg_zero_order(elem)
    #If the order of the tri-vertices changes in StartUpDG update this line.
    if elem == Tri()
        return [-1.0; -1.0]
    elseif elem == Quad()
        return [-1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0]
    end
end

@testset "VTKWriter test for $elem" for elem in [Tri(), Quad()]
    @testset "Polygrad $N" for N = [1, 2, 3, 4]
        @testset "Write Mesh" begin
            rd = RefElemData(elem, N)
            md = MeshData(uniform_mesh(elem, 4)..., rd)
            interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
            pdata = [quad.(vec(interpolate*md.x), vec(interpolate*md.y))]
            filename = replace(string(elem), "()"=>"") * "_" * string(2)
            check = filename * ".vtu"
            #Todo: Can we implement a better check?
            vtu_name = MeshData_to_vtk(md, rd, pdata, ["(x+y)^2"], filename, true)
            @test vtu_name[1] == check
        end
    end

    @testset "VTK-Node degree $order" for order=[-1, 0, 1]
        @testset "VTK order" begin
            if order == -1
                @test vtk_order(elem, order) == nothing
            elseif order == 0
                @test vtk_order(elem, order) ≈ deg_zero_order(elem)
            else
                @test vtk_order(elem, order) ≈ deg_one_order(elem)
            end
        end
    end
end