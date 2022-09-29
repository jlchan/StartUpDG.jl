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

    @testset "3D hex mesh initialization" begin
        tol = 5e2*eps()

        N = 2
        K1D = 2
        rd = RefElemData(Hex(), approxType, N)        
        md = MeshData(uniform_mesh(Hex(), K1D)..., rd)
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
        @test h ≈ 2 / K1D 

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

        # check connectivity and boundary maps
        u = @. (1-x) * (1+x) * (1-y) * (1+y) * (1-z) * (1+z)
        uf = Vf * u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        md = make_periodic(md, (true, true, true))
        @unpack mapP = md
        u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y)) * sin(pi * (.5 + z))
        uf = Vf * u
        @test uf ≈ uf[mapP]
    end
end

@testset "3D polynomial tet mesh initialization" begin
    tol = 5e2 * eps()

    N = 3
    K1D = 2
    rd = RefElemData(Tet(), Polynomial(), N)
    md = MeshData(uniform_mesh(Tet(), K1D)..., rd)
    @unpack wq, Dr, Ds, Dt, Vq, Vf, wf = rd
    Nfaces = length(rd.fv)
    @unpack x, y, z, xq, yq, zq, wJq, xf, yf, zf, K = md
    @unpack rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = md
    @unpack nxJ, nyJ, nzJ, sJ = md
    @unpack FToF, mapM, mapP, mapB = md

    @test md.mesh_type==rd.element_type
    @test md.x == md.xyz[1]

    # check positivity of Jacobian

    @test all(J .> 0)
    # h = estimate_h(rd,md)
    # @test h ≈ 2/K1D 

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
    @test Vq*x ≈ xq
    @test Vq*y ≈ yq
    @test Vq*z ≈ zq
    @test diagm(wq)*(Vq*J) ≈ wJq
    @test abs(sum(xq.*wJq)) < tol
    @test abs(sum(yq.*wJq)) < tol
    @test abs(sum(zq.*wJq)) < tol

    # check surface integration
    @test Vf * x ≈ xf
    @test Vf * y ≈ yf
    @test Vf * z ≈ zf
    @test abs(sum(diagm(wf) * nxJ)) < tol
    @test abs(sum(diagm(wf) * nyJ)) < tol
    @test abs(sum(diagm(wf) * nzJ)) < tol

    # check connectivity and boundary maps
    u = @. (1 - x) * (1 + x) * (1 - y) * (1 + y) * (1 - z) * (1 + z)
    uf = Vf * u
    @test uf ≈ uf[mapP]
    @test norm(uf[mapB]) < tol

    # check periodic node connectivity maps
    md = make_periodic(md, (true, true, true))
    @unpack mapP = md
    u = @. sin(pi * (.5 + x)) * sin(pi * (.5 + y)) * sin(pi * (.5 + z))
    uf = Vf * u
    @test uf ≈ uf[mapP]
end
