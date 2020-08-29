using StartUpDG
using Test
using UnPack
using LinearAlgebra

# avoids redundancy from tests in NodesAndModes.
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
end

@testset "2D tri/quad mesh initialization" begin
    tol = 5e2*eps()

    N = 2
    K1D = 2
    init_ref_fxns = [init_reference_tri,init_reference_quad]
    unif_mesh_fxns = [uniform_tri_mesh,uniform_quad_mesh]
    for (init_ref_elem,unif_mesh) in zip(init_ref_fxns,unif_mesh_fxns)
        rd = init_ref_elem(N)
        VX,VY,EToV = unif_mesh(K1D)
        md = init_DG_mesh(VX,VY,EToV,rd)
        @unpack wq,Dr,Ds,Vq,Vf,wf= rd
        @unpack x,y,xq,yq,xf,yf = md
        @unpack rxJ,sxJ,ryJ,syJ,J,nxJ,nyJ,sJ,wJq = md

        u = @. x^2 + 2*x*y - y^2
        dudx_exact = @. 2*x + 2*y
        dudy_exact = @. 2*x - 2*y
        dudr,duds = (A->A*u).((Dr,Ds))
        dudx = (rxJ.*dudr + sxJ.*duds)./J
        dudy = (ryJ.*dudr + syJ.*duds)./J
        @test dudx ≈ dudx_exact
        @test dudy ≈ dudy_exact

        @test Vq*x ≈ xq
        @test Vq*y ≈ yq
        @test diagm(wq)*(Vq*J) ≈ wJq
        @test abs(sum(xq.*wJq)) < tol
        @test abs(sum(yq.*wJq)) < tol

        @test Vf*x ≈ xf
        @test Vf*y ≈ yf
        @test abs(sum(diagm(wf)*nxJ)) < tol
        @test abs(sum(diagm(wf)*nyJ)) < tol
    end
end
