@testset "1D and 2D MeshData tests" begin 
    @testset "MeshData constructors" begin
        rd = RefElemData(Line(), 3)
        @test_nowarn MeshData(2, rd; coordinates_min=(-1, ), coordinates_max=(1,), is_periodic=true)
        
        for elem in [Tri(), Quad()]
            rd = RefElemData(elem, 3)
            @test_nowarn MeshData(2, 4, rd; coordinates_min=(-1, -2), coordinates_max=(1, 2), is_periodic=true)
            @test_nowarn MeshData(2, rd; coordinates_min=(-1, -2), coordinates_max=(1, 2), is_periodic=true)
        end

        for elem in [Hex(), Tet(), Pyr(), Wedge()]
            rd = RefElemData(elem, 3)
            @test_nowarn MeshData(2, 4, 8, rd; coordinates_min=(-1, -2, -4), coordinates_max=(1, 2, 4), is_periodic=true)
            @test_nowarn MeshData(2, rd; coordinates_min=(-1, -2, -4), coordinates_max=(1, 2, 4), is_periodic=true)
        end
    end
    
    @testset "$approx_type MeshData initialization" for approx_type = [Polynomial(), SBP()]
        @testset "1D mesh initialization" begin
            tol = 5e2*eps()

            N = 3
            K1D = 2
            rd = RefElemData(Line(), approx_type, N)        
            md = MeshData(uniform_mesh(Line(), K1D)..., rd)
            (; wq, Dr, Vq, Vf, wf  ) = rd
            (; Nfaces  ) = rd
            (; x, xq, xf, K  ) = md
            (; rxJ, J, nxJ, wJq  ) = md
            (; mapM, mapP, mapB  ) = md

            @test typeof(md.mesh_type) <: StartUpDG.VertexMappedMesh{<:typeof(rd.element_type)}

            @test md.x == md.xyz[1]

            # check positivity of Jacobian
            @test all(J .> 0)
            h = estimate_h(rd, md)
            @test h ≈ 2 / K1D 

            # check differentiation
            u = @. x^2 + 2*x
            dudx_exact = @. 2*x + 2
            dudr = Dr * u
            dudx = @. (rxJ * dudr) / J
            @test dudx ≈ dudx_exact

            # check volume integration
            @test Vq * x ≈ xq
            @test diagm(wq) * (Vq * J) ≈ wJq
            @test abs(sum(xq .* wJq)) < tol

            # check surface integration
            @test Vf * x ≈ xf
            @test abs(sum(nxJ)) < tol

            # check connectivity and boundary maps
            u = @. (1-x) * (1+x)
            uf = Vf * u
            @test uf ≈ uf[mapP]
            @test norm(uf[mapB]) < tol

            # check periodic node connectivity maps
            md = make_periodic(md)
            (; mapP  ) = md
            u = @. sin(pi * (.5 + x))
            uf = Vf * u
            @test uf ≈ uf[mapP]
            
            md = MeshData(uniform_mesh(Line(), K1D)..., rd; is_periodic=true)
            @test isempty(md.mapB)
        end

        @testset "2D tri mesh initialization" begin
            tol = 5e3*eps() # higher tolerance due to floating point issues?

            N = 3
            K1D = 2
            rd = RefElemData(Tri(), approx_type, N)
            md = MeshData(uniform_mesh(Tri(), K1D)..., rd)
            (; wq, Dr, Ds, Vq, Vf, wf  ) = rd
            Nfaces = length(rd.fv)
            (; x, y, xq, yq, xf, yf, K  ) = md
            (; rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ, wJq  ) = md
            (; FToF, mapM, mapP, mapB  ) = md

            @test typeof(md.mesh_type) <: StartUpDG.VertexMappedMesh{<:typeof(rd.element_type)}
            @test md.x == md.xyz[1]

            # check positivity of Jacobian
            @test all(J .> 0)
            h = estimate_h(rd, md)
            @test h ≈ 2 / K1D

            # check differentiation
            u = @. x^2 + 2 * x * y - y^2
            dudx_exact = @. 2*x + 2*y
            dudy_exact = @. 2*x - 2*y
            dudr,duds = (D->D*u).((Dr, Ds))
            dudx = @. (rxJ * dudr + sxJ * duds) / J
            dudy = @. (ryJ * dudr + syJ * duds) / J
            @test dudx ≈ dudx_exact
            @test dudy ≈ dudy_exact

            # check volume integration
            @test Vq*x ≈ xq
            @test Vq*y ≈ yq
            @test diagm(wq) * (Vq * J) ≈ wJq
            @test abs(sum(xq .* wJq)) < tol
            @test abs(sum(yq .* wJq)) < tol

            # check surface integration
            @test Vf * x ≈ xf
            @test Vf * y ≈ yf
            @test abs(sum(wf .* nxJ)) < tol
            @test abs(sum(wf .* nyJ)) < tol
            @test md.nx .* md.Jf ≈ md.nxJ
            @test md.ny .* md.Jf ≈ md.nyJ
            @test sum(@. wf * nxJ * (1 + xf) / 2) ≈ 2.0 # check sign of normals

            # check connectivity and boundary maps
            u = @. (1-x) * (1+x) * (1-y) * (1+y)
            uf = Vf * u
            @test uf ≈ uf[mapP]
            @test norm(uf[mapB]) < tol

            # check periodic node connectivity maps
            md = make_periodic(md, (true, true))
            (; mapP  ) = md
            u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y))
            uf = Vf * u
            @test uf ≈ uf[mapP]

            # check MeshData struct copying
            xyz = (x->x .+ 1).(md.xyz) # affine shift
            md2 = MeshData(rd, md, xyz...)
            @test sum(norm.(md2.rstxyzJ .- md.rstxyzJ)) < tol
            @test sum(norm.(md2.nxyzJ .- md.nxyzJ)) < tol
            @test all(md2.xyzf .≈ (x->x .+ 1).(md.xyzf))

            md = MeshData(uniform_mesh(rd.element_type, K1D)..., rd; is_periodic=true)
            @test isempty(md.mapB)

        end

        @testset "2D quad mesh initialization" begin
            tol = 5e2*eps()

            N, K1D = 3, 2
            rd = RefElemData(Quad(), approx_type, N)        
            md = MeshData(uniform_mesh(Quad(), K1D)..., rd)
            (; wq, Dr, Ds, Vq, Vf, wf  ) = rd
            Nfaces = length(rd.fv)
            (; x, y, xq, yq, xf, yf, K  ) = md
            (; rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ, wJq  ) = md
            (; FToF, mapM, mapP, mapB  ) = md

            @test typeof(md.mesh_type) <: StartUpDG.VertexMappedMesh{<:typeof(rd.element_type)}
            @test md.x == md.xyz[1]

            # check positivity of Jacobian
            @test all(J .> 0)
            h = estimate_h(rd, md)
            @test h ≈ 2 / K1D 

            # check differentiation
            u = @. x^2 + 2 * x * y - y^2
            dudx_exact = @. 2*x + 2*y
            dudy_exact = @. 2*x - 2*y
            dudr,duds = (D->D * u).((Dr, Ds))
            dudx = @. (rxJ * dudr + sxJ * duds) / J
            dudy = @. (ryJ * dudr + syJ * duds) / J
            @test dudx ≈ dudx_exact
            @test dudy ≈ dudy_exact

            # check volume integration
            @test Vq * x ≈ xq
            @test Vq * y ≈ yq
            @test diagm(wq) * (Vq * J) ≈ wJq
            @test abs(sum(xq .* wJq)) < tol
            @test abs(sum(yq .* wJq)) < tol

            # check surface integration
            @test Vf * x ≈ xf
            @test Vf * y ≈ yf
            @test abs(sum(diagm(wf) * nxJ)) < tol
            @test abs(sum(diagm(wf) * nyJ)) < tol
            @test md.nx .* md.Jf ≈ md.nxJ
            @test md.ny .* md.Jf ≈ md.nyJ
            @test sum(@. wf * nxJ * (1 + xf) / 2) ≈ 2.0 # check sign of normals

            # check connectivity and boundary maps
            u = @. (1-x) * (1+x) * (1-y) * (1+y)
            uf = Vf * u
            @test uf ≈ uf[mapP]
            @test norm(uf[mapB]) < tol

            # check periodic node connectivity maps
            md = make_periodic(md, (true, true))
            (; mapP  ) = md
            u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y))
            uf = Vf * u
            @test uf ≈ uf[mapP]

            md = MeshData(uniform_mesh(rd.element_type, K1D)..., rd; is_periodic=true)
            @test isempty(md.mapB)

        end

        @testset "2D curved tests" begin
            rd = RefElemData(Quad(), N=4)
            md = MeshData(uniform_mesh(Quad(), 4, 4)..., rd)
            (; x, y  ) = md
            x_curved = @. x + 0.1 * sin(pi * x) * sin(pi * y)
            y_curved = @. y + 0.1 * sin(pi * x) * sin(pi * y)
            md2 = MeshData(rd, md, x_curved, y_curved)
            @test sum(@. md2.wJq) ≈ 4
            @test sum(@. md2.wJq * md2.xq^2) ≈ 4/3
            @test md2.nx ≈ md2.nxJ ./ md2.Jf
            @test md2.ny ≈ md2.nyJ ./ md2.Jf
            @test typeof(md2) <: MeshData{2, <:StartUpDG.CurvedMesh}
            @test num_elements(md) == size(md.x, 2)
        end
    end
end

# Note: we test wedge and pyr elements in a separate file
approx_elem_types_to_test = [(Polynomial(), Hex()), 
                             (SBP(), Hex()), 
                             (Polynomial(), Tet()),
                             ]
@testset "3D MeshData tests" begin 
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
end

@testset "Very high orders" begin
    @testset "Tet" begin 
        N = 12
        rd = RefElemData(Tet(), N)
        md = MeshData(uniform_mesh(Tet(), 2), rd; is_periodic=true)

        (; x, y, z) = md
        u = @. sin(pi * x) * sin(pi * y) * sin(pi * z)
        dudx_exact = @. pi * cos(pi * x) * sin(pi * y) * sin(pi * z)
        dudy_exact = @. pi * sin(pi * x) * cos(pi * y) * sin(pi * z)
        dudz_exact = @. pi * sin(pi * x) * sin(pi * y) * cos(pi * z)

        dudr, duds, dudt = (D -> D * u).(rd.Drst)
        dudx = @. (md.rxJ * dudr + md.sxJ * duds + md.txJ * dudt) / md.J
        dudy = @. (md.ryJ * dudr + md.syJ * duds + md.tyJ * dudt) / md.J
        dudz = @. (md.rzJ * dudr + md.szJ * duds + md.tzJ * dudt) / md.J

        @test norm(dudx - dudx_exact, Inf) < 1e-3
        @test norm(dudy - dudy_exact, Inf) < 1e-3
        @test norm(dudz - dudz_exact, Inf) < 1e-3
    end
end
