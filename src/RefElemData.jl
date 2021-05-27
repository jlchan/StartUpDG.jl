"""
    struct RefElemData{Dim, ElemShape <: AbstractElemShape, Tv, Nfaces} 

RefElemData: contains info (interpolation points, volume/face quadrature, operators)
for a high order nodal polynomial basis on a given reference element.

Use `@unpack` to extract fields. Example:
```julia
N = 3
rd = RefElemData(Tri(),N)
@unpack r,s = rd
```
"""
struct RefElemData{Dim, ElemShape <: AbstractElemShape, Nfaces, Tv} 

    elemShape::ElemShape

    N::Int         # degree
    # Nfaces::Int    # num faces - redundant, remove later (breaking change)
                   # replace with "nfaces(rd::RefElemData{DIM,ElemShape,Tv,Nfaces}) = Nfaces"
    fv::Union{NTuple{Nfaces,Int},NTuple{Nfaces,Vector{Int}}} # list of vertices defining faces, e.g., ([1,2],[2,3],[3,1]) for a triangle
    V1::Matrix{Tv} # low order interpolation matrix

    rst::NTuple{Dim,Vector{Tv}}
    VDM::Matrix{Tv}     # generalized Vandermonde matrix

    # plotting nodes
    Nplot::Int
    rstp::NTuple{Dim,Vector{Tv}}
    Vp::Matrix{Tv}      # interpolation matrix to plotting nodes

    rstq::NTuple{Dim,Vector{Tv}}
    wq::Vector{Tv}
    Vq::Matrix{Tv}           # quad interp mat

    rstf::NTuple{Dim}
    wf::Vector{Tv}      # quad weights
    Vf::Matrix{Tv}           # face quad interp mat

    # reference normals, quad weights
    nrstJ::NTuple{Dim,Vector{Tv}}

    M::Matrix{Tv}          # mass matrix
    Pq::Matrix{Tv}               # L2 projection matrix

    # specialize diff and lift (dense, sparse, Bern, etc)
    Drst::NTuple{Dim,Matrix{Tv}} # differentiation operators
    LIFT::Matrix{Tv}             # lift matrix
end



function Base.show(io::IO, rd::RefElemData)
    @nospecialize rd
    println("Degree $(rd.N) RefElemData on $(rd.elemShape) element.")
end

# convenience unpacking routines
function Base.getproperty(x::RefElemData{Dim,ElemShape,Nfaces}, s::Symbol) where {Dim,ElemShape,Nfaces}
    if s==:r
        return getfield(x,:rst)[1]
    elseif s==:s
        return getfield(x,:rst)[2]
    elseif s==:t
        return getfield(x,:rst)[3]

    elseif s==:rq
        return getfield(x,:rstq)[1]
    elseif s==:sq
        return getfield(x,:rstq)[2]
    elseif s==:tq
        return getfield(x,:rstq)[3]

    elseif s==:rf
        return getfield(x,:rstf)[1]
    elseif s==:sf
        return getfield(x,:rstf)[2]
    elseif s==:tf
        return getfield(x,:rstf)[3]

    elseif s==:rp
        return getfield(x,:rstp)[1]
    elseif s==:sp
        return getfield(x,:rstp)[2]
    elseif s==:tp
        return getfield(x,:rstp)[3]

    elseif s==:nrJ
        return getfield(x,:nrstJ)[1]
    elseif s==:nsJ
        return getfield(x,:nrstJ)[2]
    elseif s==:ntJ
        return getfield(x,:nrstJ)[3]

    elseif s==:Dr
        return getfield(x,:Drst)[1]
    elseif s==:Ds
        return getfield(x,:Drst)[2]
    elseif s==:Dt
        return getfield(x,:Drst)[3]

    elseif s==:Nfaces
        return Nfaces

    else
        return getfield(x,s)
    end
end

"""
    RefElemData(elem::Line, N;
                quad_rule_vol = quad_nodes(elem,N+1))
    RefElemData(elem::Union{Tri,Quad}, N;
                 quad_rule_vol = quad_nodes(elem,N),
                 quad_rule_face = gauss_quad(0,0,N))
    RefElemData(elem::Hex,N;
                 quad_rule_vol = quad_nodes(elem,N),
                 quad_rule_face = quad_nodes(Quad(),N))
    RefElemData(elem; N, kwargs...) # version with keyword arg

Constructor for RefElemData for different element types.
"""
RefElemData(elem; N, kwargs...) = RefElemData(elem, N; kwargs...)

function RefElemData(elem::Line, N; quad_rule_vol = quad_nodes(elem,N+1), Nplot=10)

    fv = face_vertices(elem)
    Nfaces = length(fv)

    # Construct matrices on reference elements
    r = nodes(elem,N)
    VDM = vandermonde(elem, N, r)
    Dr = grad_vandermonde(elem, N, r)/VDM

    V1 = vandermonde(elem,1,r)/vandermonde(elem,1,[-1;1])

    rq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    rf = [-1.0;1.0]
    nrJ = [-1.0;1.0]
    wf = [1.0;1.0]
    Vf = vandermonde(elem,N,rf)/VDM
    LIFT = M\(Vf') # lift matrix

    # plotting nodes
    rp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp)/VDM

    return RefElemData(elem,N,fv,V1,
                       tuple(r),VDM,
                       Nplot,tuple(rp),Vp,
                       tuple(rq),wq,Vq,
                       tuple(rf),wf,Vf,tuple(nrJ),
                       M,Pq,tuple(Dr),LIFT)
end

function RefElemData(elem::Union{Tri,Quad}, N;
                     quad_rule_vol = quad_nodes(elem,N),
                     quad_rule_face = gauss_quad(0,0,N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle
    Nfaces = length(fv)

    # Construct matrices on reference elements
    r,s = nodes(elem,N)
    VDM,Vr,Vs = basis(elem,N,r,s)
    Dr = Vr/VDM
    Ds = Vs/VDM

    # low order interpolation nodes
    r1,s1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s)/vandermonde(elem,1,r1,s1)

    rf,sf,wf,nrJ,nsJ = init_face_data(elem,N,quad_nodes_face=quad_rule_face)

    rq,sq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq,sq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    Vf = vandermonde(elem,N,rf,sf)/VDM # interpolates from nodes to face nodes
    LIFT = M\(Vf'*diagm(wf)) # lift matrix used in rhs evaluation

    # plotting nodes
    rp, sp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp,sp)/VDM

    Drs = (Dr,Ds)

    # sparsify for Quad
    # tol = 1e-13
    # Drs = typeof(elem)==Quad ? droptol!.(sparse.((Dr,Ds)),tol) : (Dr,Ds)
    # Vf = typeof(elem)==Quad ? droptol!(sparse(Vf),tol) : Vf
    # LIFT = typeof(elem)==Quad ? droptol!(sparse(LIFT),tol) : LIFT

    return RefElemData(elem,N,fv,V1,
                       tuple(r,s),VDM,
                       Nplot,tuple(rp,sp),Vp,
                       tuple(rq,sq),wq,Vq,
                       tuple(rf,sf),wf,Vf,tuple(nrJ,nsJ),
                       M,Pq,Drs,LIFT)
end

function RefElemData(elem::Hex,N;
                     quad_rule_vol = quad_nodes(elem,N),
                     quad_rule_face = quad_nodes(Quad(),N),
                     Nplot=10)

    fv = face_vertices(elem) # set faces for triangle
    Nfaces = length(fv)

    # Construct matrices on reference elements
    r,s,t = nodes(elem,N)
    VDM,Vr,Vs,Vt = basis(elem,N,r,s,t)
    Dr,Ds,Dt = (A->A/VDM).((Vr,Vs,Vt))

    # low order interpolation nodes
    r1,s1,t1 = nodes(elem,1)
    V1 = vandermonde(elem,1,r,s,t)/vandermonde(elem,1,r1,s1,t1)

    #Nodes on faces, and face node coordinate
    rf,sf,tf,wf,nrJ,nsJ,ntJ = init_face_data(elem,N)

    # quadrature nodes - build from 1D nodes.
    rq,sq,tq,wq = quad_rule_vol
    Vq = vandermonde(elem,N,rq,sq,tq)/VDM
    M = Vq'*diagm(wq)*Vq
    Pq = M\(Vq'*diagm(wq))

    Vf = vandermonde(elem,N,rf,sf,tf)/VDM
    LIFT = M\(Vf'*diagm(wf))

    # plotting nodes
    rp,sp,tp = equi_nodes(elem,Nplot)
    Vp = vandermonde(elem,N,rp,sp,tp)/VDM

    # Drst = sparse.((Dr,Ds,Dt))
    Drst = (Dr,Ds,Dt)
    # Vf = sparse(Vf)

    return RefElemData(elem,N,fv,V1,
                       tuple(r,s,t),VDM,
                       Nplot,tuple(rp,sp,tp),Vp,
                       tuple(rq,sq,tq),wq,Vq,
                       tuple(rf,sf,tf),wf,Vf,tuple(nrJ,nsJ,ntJ),
                       M,Pq,Drst,LIFT)
end
