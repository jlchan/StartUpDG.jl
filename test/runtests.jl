using Test: testset_beginend
using Suppressor
using StartUpDG
using Test
using LinearAlgebra
using Triangulate

@testset "Mesh, timestep, and Triangulate utils" begin
    tol = 5e2*eps()

    rk4a,rk4b,rk4c = ck45()
    @test rk4c[1] ≈ 0.0 && rk4c[end] ≈ 1.0
    @test rk4a[1] ≈ 0.0
    @test all(rk4b .> 0)

    VX,VY,EToV = readGmsh2D("squareCylinder2D.msh")
    @test size(EToV)==(3031,3)

    # test triangulate
    meshIO = square_domain()
    VX,VY,EToV = triangulateIO_to_VXYEToV(meshIO)
    rd = RefElemData(Tri(),2)
    md = MeshData(VX,VY,EToV,rd)
    @test size(EToV,1)==md.num_elements==620
    @test length(VX)==length(VY)==338
    @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2,3,4]

    h = .1
    meshIO = square_hole_domain(h)
    VX,VY,EToV = triangulateIO_to_VXYEToV(meshIO)
    rd = RefElemData(Tri(),2)
    md = MeshData(VX,VY,EToV,rd)
    @test size(EToV,1)==md.num_elements==598
    @test length(VX)==length(VY)==327
    @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2]
    meshIO2 = refine(meshIO,h)
    @test sort(unique(get_node_boundary_tags(meshIO2,rd,MeshData(triangulateIO_to_VXYEToV(meshIO2)...,rd))))==[0,1,2]

    meshIO = scramjet()
    VX,VY,EToV = triangulateIO_to_VXYEToV(meshIO)
    rd = RefElemData(Tri(),2)
    md = MeshData(VX,VY,EToV,rd)
    @test size(EToV,1)==md.num_elements==1550
    @test length(VX)==length(VY)==871
    @test sort(unique(get_node_boundary_tags(meshIO,rd,md)))==[0,1,2,3]
end

@testset "Geometric terms for $elem elements" for elem in [Tri() Quad() Hex()]
    tol = 5e2*eps()
    N = 3
    rd = RefElemData(elem,N)
    geofacs = geometric_factors(rd.rst...,rd.Drst...)
    if elem != Hex()
        rxJ,sxJ,ryJ,syJ,J = geofacs    
        @test all(rxJ .≈ 1)
        @test norm(sxJ) < tol
        @test norm(ryJ) < tol
        @test all(syJ .≈ 1)
        @test all(J .≈ 1)
    else
        rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = geofacs
        @test all(rxJ .≈ 1)
        @test norm(sxJ) < tol
        @test norm(txJ) < tol
        @test norm(ryJ) < tol
        @test all(syJ .≈ 1)
        @test norm(tyJ) < tol
        @test norm(rzJ) < tol
        @test norm(szJ) < tol
        @test all(tzJ .≈ 1)
        @test all(J .≈ 1)
    end

end

# some code not tested to avoid redundancy from tests in NodesAndModes.
@testset "Reference elements" begin
    tol = 5e2*eps()

    N = 2

    #####
    ##### interval
    #####
    rd = RefElemData(Line(),N)
    @test rd.r == rd.rst[1]
    @test abs(sum(rd.rq.*rd.wq)) < tol
    @test rd.nrJ ≈ [-1,1]
    @test rd.Pq*rd.Vq ≈ I
    @test rd.r[rd.Fmask[:]] ≈ rd.rf
    @test invoke(inverse_trace_constant,Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)

    #####
    ##### triangles
    #####
    rd = RefElemData(Tri(),N)
    @test rd.r == rd.rst[1]
    @test abs(sum(rd.wq)) ≈ 2
    @test abs(sum(rd.wf)) ≈ 6
    @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
    @test rd.Pq*rd.Vq ≈ I
    Vfp = vandermonde(Line(),N,quad_nodes(Line(),N)[1])/vandermonde(Line(),N,nodes(Line(),N))
    rstf = (x->Vfp*x[reshape(rd.Fmask,rd.Nfq÷rd.Nfaces,rd.Nfaces)]).(rd.rst)
    @test all(vec.(rstf) .≈ rd.rstf)
    @test invoke(inverse_trace_constant,Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)

    #####
    ##### quads
    #####
    rd = RefElemData(Quad(),N)
    @test rd.r == rd.rst[1]
    @test abs(sum(rd.wq)) ≈ 4
    @test abs(sum(rd.wf)) ≈ 8
    @test abs(sum(rd.wf .* rd.nrJ)) + abs(sum(rd.wf .* rd.nsJ)) < tol
    @test rd.Pq*rd.Vq ≈ I
    Vfp = vandermonde(Line(),N,quad_nodes(Line(),N)[1])/vandermonde(Line(),N,nodes(Line(),N))
    rstf = (x->Vfp*x[reshape(rd.Fmask,rd.Nfq÷rd.Nfaces,rd.Nfaces)]).(rd.rst)
    @test all(vec.(rstf) .≈ rd.rstf)
    @test invoke(inverse_trace_constant,Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)

    #####
    ##### hexes
    #####
    rd = RefElemData(Hex(),N)
    @test rd.r == rd.rst[1]
    @test abs(sum(rd.wq)) ≈ 8
    @test abs(sum(rd.wf)) ≈ 6*4
    @test abs(sum(rd.wf .* rd.nrJ)) < tol
    @test abs(sum(rd.wf .* rd.nsJ)) < tol
    @test abs(sum(rd.wf .* rd.ntJ)) < tol
    @test rd.Pq*rd.Vq ≈ I
    @test invoke(inverse_trace_constant,Tuple{RefElemData},rd) ≈ inverse_trace_constant(rd)
    # TODO: test interpolation of Fmask matches rd.rstf.
end

@testset "$approxType MeshData initialization" for approxType = [Polynomial(),SBP()]
    @testset "1D mesh initialization" begin
        tol = 5e2*eps()

        N = 3
        K1D = 2
        rd = RefElemData(Line(),approxType,N)
        VX,EToV = uniform_mesh(Line(),K1D)
        md = MeshData(VX,EToV,rd)
        @unpack wq,Dr,Vq,Vf,wf = rd
        @unpack Nfaces = rd
        @unpack x,xq,xf,K = md
        @unpack rxJ,J,nxJ,wJq = md
        @unpack mapM,mapP,mapB = md

        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        @test all(J .> 0)
        h = estimate_h(rd,md)
        @test h ≈ 2/K1D 

        # check differentiation
        u = @. x^2 + 2*x
        dudx_exact = @. 2*x + 2
        dudr = Dr*u
        dudx = (rxJ.*dudr)./J
        @test dudx ≈ dudx_exact

        # check volume integration
        @test Vq*x ≈ xq
        @test diagm(wq)*(Vq*J) ≈ wJq
        @test abs(sum(xq.*wJq)) < tol

        # check surface integration
        @test Vf*x ≈ xf
        @test abs(sum(nxJ)) < tol

        # check connectivity and boundary maps
        u = @. (1-x)*(1+x)
        uf = Vf*u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        md = make_periodic(md)
        @unpack mapP = md
        u = @. sin(pi*(.5+x))
        uf = Vf*u
        @test uf ≈ uf[mapP]
    end

    @testset "2D tri mesh initialization" begin
        tol = 5e3*eps() # higher tolerance due to floating point issues?

        N = 3
        K1D = 2
        rd = RefElemData(Tri(),approxType,N)
        VX,VY,EToV = uniform_mesh(Tri(),K1D)
        md = MeshData(VX,VY,EToV,rd)
        @unpack wq,Dr,Ds,Vq,Vf,wf = rd
        Nfaces = length(rd.fv)
        @unpack x,y,xq,yq,xf,yf,K = md
        @unpack rxJ,sxJ,ryJ,syJ,J,nxJ,nyJ,sJ,wJq = md
        @unpack FToF,mapM,mapP,mapB = md

        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        # @show J[1,:]
        @test all(J .> 0)
        h = estimate_h(rd,md)
        @test h ≈ 2/K1D

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
        @test abs(sum(wf.*nxJ)) < tol
        @test abs(sum(wf.*nyJ)) < tol
        @test sum(@. wf*nxJ*(1+xf)/2) ≈ 2.0 # check sign of normals

        # check connectivity and boundary maps
        u = @. (1-x)*(1+x)*(1-y)*(1+y)
        uf = Vf*u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        md = make_periodic(md,(true,true))
        @unpack mapP = md
        u = @. sin(pi*(.5+x))*sin(pi*(.5+y))
        uf = Vf*u
        @test uf ≈ uf[mapP]

        # check MeshData struct copying
        xyz = (x->x .+ 1).(md.xyz) # affine shift
        md2 = MeshData(md,rd,xyz...)
        @test sum(norm.(md2.rstxyzJ .- md.rstxyzJ)) < tol
        @test sum(norm.(md2.nxyzJ .- md.nxyzJ)) < tol
        @test all(md2.xyzf .≈ (x->x .+ 1).(md.xyzf))
    end

    @testset "2D quad mesh initialization" begin
        tol = 5e2*eps()

        N = 3
        K1D = 2
        rd = RefElemData(Quad(),approxType,N)
        VX,VY,EToV = uniform_mesh(Quad(),K1D)
        md = MeshData(VX,VY,EToV,rd)
        @unpack wq,Dr,Ds,Vq,Vf,wf = rd
        Nfaces = length(rd.fv)
        @unpack x,y,xq,yq,xf,yf,K = md
        @unpack rxJ,sxJ,ryJ,syJ,J,nxJ,nyJ,sJ,wJq = md
        @unpack FToF,mapM,mapP,mapB = md

        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        @test all(J .> 0)
        h = estimate_h(rd,md)
        @test h ≈ 2/K1D 

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
        @test sum(@. wf*nxJ*(1+xf)/2) ≈ 2.0 # check sign of normals

        # check connectivity and boundary maps
        u = @. (1-x)*(1+x)*(1-y)*(1+y)
        uf = Vf*u
        @test uf ≈ uf[mapP]
        @test norm(uf[mapB]) < tol

        # check periodic node connectivity maps
        md = make_periodic(md,(true,true))
        @unpack mapP = md
        u = @. sin(pi*(.5+x))*sin(pi*(.5+y))
        uf = Vf*u
        @test uf ≈ uf[mapP]
    end

    @testset "3D hex mesh initialization" begin
        tol = 5e2*eps()

        N = 2
        K1D = 2
        rd = RefElemData(Hex(),approxType,N)
        VX,VY,VZ,EToV = uniform_mesh(Hex(),K1D)
        md = MeshData(VX,VY,VZ,EToV,rd)
        @unpack wq,Dr,Ds,Dt,Vq,Vf,wf = rd
        Nfaces = length(rd.fv)
        @unpack x,y,z,xq,yq,zq,wJq,xf,yf,zf,K = md
        @unpack rxJ,sxJ,txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ,J = md
        @unpack nxJ,nyJ,nzJ,sJ = md
        @unpack FToF,mapM,mapP,mapB = md

        @test md.x == md.xyz[1]

        # check positivity of Jacobian
        # @show J[1,:]
        @test all(J .> 0)
        h = estimate_h(rd,md)
        @test h ≈ 2/K1D 


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
        md = make_periodic(md,(true,true,true))
        @unpack mapP = md
        u = @. sin(pi*(.5+x))*sin(pi*(.5+y))*sin(pi*(.5+z))
        uf = Vf*u
        @test uf ≈ uf[mapP]
    end
end

@testset "Boundary condition utils" begin
    rd = RefElemData(Tri(),N=3)
    md = MeshData(uniform_mesh(Tri(),1)...,rd)
    on_bottom_boundary(x,y,tol=1e-13) = abs(y+1) < tol
    on_top_boundary(x,y,tol=1e-13) = abs(y-1) < tol
    boundary_dict = tag_boundary_faces(md,Dict(:bottom=>on_bottom_boundary,:top=>on_top_boundary))
    @test boundary_dict == Dict(:bottom=>[1],:top=>[4])
    boundary_dict = tag_boundary_faces(md,nothing)
    @test boundary_dict == Dict(:entire_boundary=>[1,2,4,5])
end

@testset "Base.show tests" begin
    rd = RefElemData(Tri(),N=3)
    md = MeshData(uniform_mesh(Tri(),1)...,rd)
    @test (@capture_out Base.show(stdout,MIME"text/plain"(),rd)) == "RefElemData for a degree 3 Polynomial() approximation on Tri() element."
    @test (@capture_out Base.show(stdout,rd)) == "RefElemData{N=3,Polynomial(),Tri()}."
    @test (@capture_out Base.show(stdout,MIME"text/plain"(),md)) == "MeshData of dimension 2 with 2 elements"
    @test (@capture_out Base.show(stdout,md)) == "MeshData{2}"
end