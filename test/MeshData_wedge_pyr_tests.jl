approx_elem_types_to_test = [(Polynomial(), Wedge()),
                             (Polynomial(), Pyr())]

@testset "3D MeshData tests for wedges and pyramids" begin 

    # tests initialization (mostly checking if tolerances)
    @testset "Initialization test for N = 1 for $element_type" for element_type in [Pyr(), Wedge()]
        K1D = 2
        rd = RefElemData(element_type, 1)        
        md = MeshData(uniform_mesh(element_type, K1D)..., rd; is_periodic=true)
    end

    @testset "$approximation_type $element_type MeshData initialization" for (approximation_type, element_type) in approx_elem_types_to_test
        tol = 5e2*eps()

        N = 3
        K1D = 2
        rd = RefElemData(element_type, approximation_type, N)        
        md = MeshData(uniform_mesh(element_type, K1D)..., rd)
        (; wq, Dr, Ds, Dt, Vq, Vf, wf  ) = rd
        Nfaces = length(rd.fv)
        (; x, y, z, xq, yq, zq, wJq, xf, yf, zf, K  ) = md
        (; rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J  ) = md
        (; nxJ, nyJ, nzJ, sJ  ) = md
        (; FToF, mapM, mapP, mapB  ) = md

        @test typeof(md.mesh_type) <: StartUpDG.VertexMappedMesh{<:typeof(rd.element_type)}
        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        @test all(J .> 0)
        h = estimate_h(rd, md)        
        @test h <= 2 / K1D + tol

        # check differentiation
        u = @. x^2 + 2 * x * y - y^2 + x * y * z
        dudx_exact = @. 2*x + 2*y + y*z
        dudy_exact = @. 2*x - 2*y + x*z
        dudz_exact = @. x*y
        dudr,duds,dudt = (D->D*u).((Dr, Ds, Dt))
        dudx = @. (rxJ * dudr + sxJ * duds + txJ * dudt) / J
        dudy = @. (ryJ * dudr + syJ * duds + tyJ * dudt) / J
        dudz = @. (rzJ * dudr + szJ * duds + tzJ * dudt) / J
        @test dudx ≈ dudx_exact
        @test dudy ≈ dudy_exact
        @test dudz ≈ dudz_exact

        # check volume integration
        @test Vq * x ≈ xq
        @test Vq * y ≈ yq
        @test Vq * z ≈ zq
        @test diagm(wq) * (Vq * J) ≈ wJq
        @test abs(sum(xq .* wJq)) < tol
        @test abs(sum(yq .* wJq)) < tol
        @test abs(sum(zq .* wJq)) < tol

        # check surface integration
        @test Vf * x ≈ xf
        @test Vf * y ≈ yf
        @test Vf * z ≈ zf
        @test abs(sum(diagm(wf) * nxJ)) < tol
        @test abs(sum(diagm(wf) * nyJ)) < tol
        @test abs(sum(diagm(wf) * nzJ)) < tol
        @test md.nx .* md.Jf ≈ md.nxJ
        @test md.ny .* md.Jf ≈ md.nyJ
        @test md.nz .* md.Jf ≈ md.nzJ

        # check connectivity and boundary maps
        u = @. (1-x) * (1+x) * (1-y) * (1+y) * (1-z) * (1+z)
        uf = Vf * u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        md_periodic = make_periodic(md, (true, true, true))
        @test md_periodic.mapP != md.mapP # check that the node mapping actually changed

        u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y)) * sin(pi * (.5 + z))
        (; mapP  ) = md_periodic
        uf = Vf * u
        @test uf ≈ uf[mapP] 

        md = MeshData(uniform_mesh(rd.element_type, K1D)..., rd; is_periodic=true)
        @test isempty(md.mapB)
        
    end
    
    @testset "TensorProductWedge MeshData" begin
        element_type = Wedge()
        tol = 5e2*eps()
        @testset "Degree $tri_grad triangle" for tri_grad = [2, 3]
            @testset "Degree $line_grad line" for line_grad = [2, 3]
                line = RefElemData(Line(), line_grad)
                tri  = RefElemData(Tri(), tri_grad)
                tensor = TensorProductWedge(tri, line)
                
                rd = RefElemData(element_type, tensor) 
                K1D = 2       
                md = MeshData(uniform_mesh(element_type, K1D)..., rd)
                (; wq, Dr, Ds, Dt, Vq, Vf, wf  ) = rd
                Nfaces = length(rd.fv)
                (; x, y, z, xq, yq, zq, wJq, xf, yf, zf, K  ) = md
                (; rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J  ) = md
                (; nxJ, nyJ, nzJ, sJ  ) = md
                (; FToF, mapM, mapP, mapB  ) = md

                @test StartUpDG._short_typeof(rd.approximation_type) == "TensorProductWedge{Polynomial, Polynomial}"
                @test typeof(md.mesh_type) <: StartUpDG.VertexMappedMesh{<:typeof(rd.element_type)}
                @test md.x == md.xyz[1]
                @test md.y == md.xyz[2]
                @test md.z == md.xyz[3]

                # check positivity of Jacobian
                @test all(J .> 0)
                h = estimate_h(rd, md)        
                @test h <= 2 / K1D + tol

                # check differentiation
                u = @. x^2 + 2 * x * y - y^2 + x * y * z
                dudx_exact = @. (2*x + 2*y + y*z)
                dudy_exact = @. 2*x - 2*y + x*z
                dudz_exact = @. x*y
                dudr,duds,dudt = (D->D*u).((Dr, Ds, Dt))
                dudx = @. (rxJ * dudr + sxJ * duds + txJ * dudt) / J
                dudy = @. (ryJ * dudr + syJ * duds + tyJ * dudt) / J
                dudz = @. (rzJ * dudr + szJ * duds + tzJ * dudt) / J
                @test dudx ≈ dudx_exact
                @test dudy ≈ dudy_exact
                @test dudz ≈ dudz_exact

                # check volume integration
                @test Vq * x ≈ xq
                @test Vq * y ≈ yq
                @test Vq * z ≈ zq
                @test diagm(wq) * (Vq * J) ≈ wJq
                @test abs(sum(xq .* wJq)) < tol
                @test abs(sum(yq .* wJq)) < tol
                @test abs(sum(zq .* wJq)) < tol

                # check surface integration
                @test Vf * x ≈ xf
                @test Vf * y ≈ yf
                @test Vf * z ≈ zf
                @test abs(sum(diagm(wf) * nxJ)) < tol
                @test abs(sum(diagm(wf) * nyJ)) < tol
                @test abs(sum(diagm(wf) * nzJ)) < tol
                @test md.nx .* md.Jf ≈ md.nxJ
                @test md.ny .* md.Jf ≈ md.nyJ
                @test md.nz .* md.Jf ≈ md.nzJ

                # check connectivity and boundary maps
                u = @. (1-x) * (1+x) * (1-y) * (1+y) * (1-z) * (1+z)
                uf = Vf * u
                @test uf ≈ uf[mapP]
                @test norm(uf[mapB]) < tol

                # check periodic node connectivity maps
                md_periodic = make_periodic(md, (true, true, true))
                @test md_periodic.mapP != md.mapP # check that the node mapping actually changed

                u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y)) * sin(pi * (.5 + z))
                (; mapP  ) = md_periodic
                uf = Vf * u
                @test uf ≈ uf[mapP] 

                md = MeshData(uniform_mesh(rd.element_type, K1D)..., rd; is_periodic=true)
                @test isempty(md.mapB)
    
            end
        end
    end
end
