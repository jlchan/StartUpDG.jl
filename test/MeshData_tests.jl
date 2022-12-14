@testset "1D and 2D MeshData tests" begin 
    @testset "$approxType MeshData initialization" for approxType = [Polynomial(), SBP()]
        @testset "1D mesh initialization" begin
            tol = 5e2*eps()

            N = 3
            K1D = 2
            rd = RefElemData(Line(), approxType, N)        
            md = MeshData(uniform_mesh(Line(), K1D)..., rd)
            @unpack wq, Dr, Vq, Vf, wf = rd
            @unpack Nfaces = rd
            @unpack x, xq, xf, K = md
            @unpack rxJ, J, nxJ, wJq = md
            @unpack mapM, mapP, mapB = md

            @test md.mesh_type==rd.element_type

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
            @unpack mapP = md
            u = @. sin(pi * (.5 + x))
            uf = Vf * u
            @test uf ≈ uf[mapP]
        end

        @testset "2D tri mesh initialization" begin
            tol = 5e3*eps() # higher tolerance due to floating point issues?

            N = 3
            K1D = 2
            rd = RefElemData(Tri(), approxType, N)
            md = MeshData(uniform_mesh(Tri(), K1D)..., rd)
            @unpack wq, Dr, Ds, Vq, Vf, wf = rd
            Nfaces = length(rd.fv)
            @unpack x, y, xq, yq, xf, yf, K = md
            @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ, wJq = md
            @unpack FToF, mapM, mapP, mapB = md

            @test md.mesh_type==rd.element_type
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
            @unpack mapP = md
            u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y))
            uf = Vf * u
            @test uf ≈ uf[mapP]

            # check MeshData struct copying
            xyz = (x->x .+ 1).(md.xyz) # affine shift
            md2 = MeshData(rd, md, xyz...)
            @test sum(norm.(md2.rstxyzJ .- md.rstxyzJ)) < tol
            @test sum(norm.(md2.nxyzJ .- md.nxyzJ)) < tol
            @test all(md2.xyzf .≈ (x->x .+ 1).(md.xyzf))
        end

        @testset "2D quad mesh initialization" begin
            tol = 5e2*eps()

            N, K1D = 3, 2
            rd = RefElemData(Quad(), approxType, N)        
            md = MeshData(uniform_mesh(Quad(), K1D)..., rd)
            @unpack wq, Dr, Ds, Vq, Vf, wf = rd
            Nfaces = length(rd.fv)
            @unpack x, y, xq, yq, xf, yf, K = md
            @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ, wJq = md
            @unpack FToF, mapM, mapP, mapB = md

            @test md.mesh_type==rd.element_type
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
            @unpack mapP = md
            u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y))
            uf = Vf * u
            @test uf ≈ uf[mapP]
        end

        @testset "2D curved tests" begin
            rd = RefElemData(Quad(), N=4)
            md = MeshData(uniform_mesh(Quad(), 4, 4)..., rd)
            @unpack x, y = md
            x_curved = @. x + 0.1 * sin(pi * x) * sin(pi * y)
            y_curved = @. y + 0.1 * sin(pi * x) * sin(pi * y)
            md = MeshData(rd, md, x_curved, y_curved)
            @test sum(@. md.wJq) ≈ 4
            @test sum(@. md.wJq * md.xq^2) ≈ 4/3
        end
    end
end

approx_elem_types_to_test = [(Polynomial(), Hex()), 
                             (SBP(), Hex()), 
                             (Polynomial(), Tet()),
                             (Polynomial(), Wedge())]
@testset "3D MeshData tests" begin 
    @testset "$approximation_type $element_type MeshData initialization" for (approximation_type, element_type) in approx_elem_types_to_test
        tol = 5e2*eps()

        N = 3
        K1D = 2
        rd = RefElemData(element_type, approximation_type, N)        
        md = MeshData(uniform_mesh(element_type, K1D)..., rd)
        @unpack wq, Dr, Ds, Dt, Vq, Vf, wf = rd
        Nfaces = length(rd.fv)
        @unpack x, y, z, xq, yq, zq, wJq, xf, yf, zf, K = md
        @unpack rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = md
        @unpack nxJ, nyJ, nzJ, sJ = md
        @unpack FToF, mapM, mapP, mapB = md

        @test md.mesh_type==rd.element_type
        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        # @show J[1,:]
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
        @unpack mapP = md_periodic
        uf = Vf * u
        @test uf ≈ uf[mapP] 
        
    end
end
