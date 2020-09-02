using StartUpDG
using Test
using UnPack
using LinearAlgebra

@testset "Other utils" begin
    @test eye(4)≈I

    EToV,VX,VY = readGmsh2D("squareCylinder2D.msh")
    @test size(EToV)==(3031,3)
end

# some code not tested to avoid redundancy from tests in NodesAndModes.
@testset "Reference elements" begin
    tol = 5e2*eps()

    N = 2

    #####
    ##### interval
    #####
    rd = init_reference_interval(N)
    @test abs(sum(rd.rq.*rd.wq)) < tol
    @test rd.nrJ ≈ [-1,1]
    @test rd.Pq*rd.Vq ≈ I

    #####
    ##### triangles
    #####
    rd = init_reference_tri(N)
    @test abs(sum(rd.wq)) ≈ 2
    @test abs(sum(rd.wf)) ≈ 6
    @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
    @test rd.Pq*rd.Vq ≈ I

    #####
    ##### quads
    #####
    rd = init_reference_quad(N)
    @test abs(sum(rd.wq)) ≈ 4
    @test abs(sum(rd.wf)) ≈ 8
    @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
    @test rd.Pq*rd.Vq ≈ I

    #####
    ##### hexes
    #####
    rd = init_reference_hex(N)
    @test abs(sum(rd.wq)) ≈ 8
    @test abs(sum(rd.wf)) ≈ 6*4
    @test abs(sum(rd.wf .* rd.nrJ)) < tol
    @test abs(sum(rd.wf .* rd.nsJ)) < tol
    @test abs(sum(rd.wf .* rd.ntJ)) < tol
    @test rd.Pq*rd.Vq ≈ I
end

@testset "2D tri/quad mesh initialization" begin
    tol = 5e2*eps()

    N = 2
    K1D = 2
    init_ref_elem = [init_reference_tri, init_reference_quad]
    unif_mesh = [uniform_tri_mesh, uniform_quad_mesh]
    for (init_ref_elem,unif_mesh) in zip(init_ref_elem,unif_mesh)
        rd = init_ref_elem(N)
        VX,VY,EToV = unif_mesh(K1D)
        md = init_DG_mesh(VX,VY,EToV,rd)
        @unpack wq,Dr,Ds,Vq,Vf,wf = rd
        Nfaces = length(rd.fv)
        @unpack x,y,xq,yq,xf,yf,K = md
        @unpack rxJ,sxJ,ryJ,syJ,J,nxJ,nyJ,sJ,wJq = md
        @unpack FToF,mapM,mapP,mapB = md

        # check differentiation
        u = @. x^2 + 2*x*y - y^2
        dudx_exact = @. 2*x + 2*y
        dudy_exact = @. 2*x - 2*y
        dudr,duds = (D->D*u).((Dr,Ds))
        dudx = (rxJ.*dudr + sxJ.*duds)./J
        dudy = (ryJ.*dudr + syJ.*duds)./J
        @test dudx ≈ dudx_exact
        @test dudy ≈ dudy_exact

        # check volume integration
        @test Vq*x ≈ xq
        @test Vq*y ≈ yq
        @test diagm(wq)*(Vq*J) ≈ wJq
        @test abs(sum(xq.*wJq)) < tol
        @test abs(sum(yq.*wJq)) < tol

        # check surface integration
        @test Vf*x ≈ xf
        @test Vf*y ≈ yf
        @test abs(sum(diagm(wf)*nxJ)) < tol
        @test abs(sum(diagm(wf)*nyJ)) < tol

        # check connectivity and boundary maps
        u = @. (1-x)*(1+x)*(1-y)*(1+y)
        uf = Vf*u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        LX,LY = 2,2
        build_periodic_boundary_maps!(md,rd,LX,LY)
        @unpack mapP = md
        #mapPB = build_periodic_boundary_maps(xf,yf,LX,LY,Nfaces*K,mapM,mapP,mapB)
        #mapP[mapB] = mapPB
        u = @. sin(pi*(.5+x))*sin(pi*(.5+y))
        uf = Vf*u
        @test uf ≈ uf[mapP]
    end
end

@testset "3D hex mesh initialization" begin
    tol = 5e2*eps()

    N = 2
    K1D = 2
    init_ref_elem = [init_reference_hex]
    unif_mesh = [uniform_hex_mesh]
    for (init_ref_elem,unif_mesh) in zip(init_ref_elem,unif_mesh)
        rd = init_ref_elem(N)
        VX,VY,VZ,EToV = unif_mesh(K1D)
        md = init_DG_mesh(VX,VY,VZ,EToV,rd)
        @unpack wq,Dr,Ds,Dt,Vq,Vf,wf = rd
        Nfaces = length(rd.fv)
        @unpack x,y,z,xq,yq,zq,wJq,xf,yf,zf,K = md
        @unpack rxJ,sxJ,txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ,J = md
        @unpack nxJ,nyJ,nzJ,sJ = md
        @unpack FToF,mapM,mapP,mapB = md

        # check differentiation
        u = @. x^2 + 2*x*y - y^2 + x*y*z
        dudx_exact = @. 2*x + 2*y + y*z
        dudy_exact = @. 2*x - 2*y + x*z
        dudz_exact = @. x*y
        dudr,duds,dudt = (D->D*u).((Dr,Ds,Dt))
        dudx = (rxJ.*dudr + sxJ.*duds + txJ.*dudt)./J
        dudy = (ryJ.*dudr + syJ.*duds + tyJ.*dudt)./J
        dudz = (rzJ.*dudr + szJ.*duds + tzJ.*dudt)./J
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
        @test Vf*x ≈ xf
        @test Vf*y ≈ yf
        @test Vf*z ≈ zf
        @test abs(sum(diagm(wf)*nxJ)) < tol
        @test abs(sum(diagm(wf)*nyJ)) < tol
        @test abs(sum(diagm(wf)*nzJ)) < tol

        # check connectivity and boundary maps
        u = @. (1-x)*(1+x)*(1-y)*(1+y)*(1-z)*(1+z)
        uf = Vf*u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        LX,LY,LZ = 2,2,2
        build_periodic_boundary_maps!(md,rd,LX,LY,LZ)
        @unpack mapP = md
        # mapPB = build_periodic_boundary_maps(
        #     xf,yf,zf,LX,LY,LZ,Nfaces*K,mapM,mapP,mapB)
        # mapP[mapB] = mapPB
        u = @. sin(pi*(.5+x))*sin(pi*(.5+y))*sin(pi*(.5+z))
        uf = Vf*u
        @test uf ≈ uf[mapP]
    end
end
