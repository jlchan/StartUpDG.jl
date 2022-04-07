"""
    struct RefElemData{Dim, ElemShape <: AbstractElemShape, ApproximationType, Nfaces, Tv, VQ, VF, MM, P, D, L} 

RefElemData: contains info (interpolation points, volume/face quadrature, operators)
for a high order nodal basis on a given reference element. 

Example:
```julia
N = 3
rd = RefElemData(Tri(), N)
@unpack r, s = rd
```
""" 
struct RefElemData{Dim, ElemShape <: AbstractElemShape, ApproximationType, Nfaces, 
                   Tv, TVDM, VQ, VF, MM, P, D, L, VP,
                   RST, RSTP, RSTQ, WQ, RSTF, WF, FMASK, NRSTJ} 

    elementType::ElemShape
    approximationType::ApproximationType # Polynomial / SBP{...}

    N::Int               # polynomial degree of accuracy
    fv::Union{NTuple{Nfaces, Int}, NTuple{Nfaces, Vector{Int}}} # list of vertices defining faces, e.g., ([1,2],[2,3],[3,1]) for a triangle
    V1::Matrix{Tv}       # low order interpolation matrix

    rst::RST
    VDM::TVDM     # generalized Vandermonde matrix
    Fmask::FMASK   # indices of face nodes

    # plotting nodes: TODO - remove? Probably doesn't need to be in RefElemData
    Nplot::Int
    rstp::RSTP
    Vp::VP      # interpolation matrix to plotting nodes

    # quadrature 
    rstq::RSTQ
    wq::WQ
    Vq::VQ              # quad interp mat

    # face quadrature 
    rstf::RSTF
    wf::WF      # quad weights
    Vf::VF              # face quad interp mat
    nrstJ::NRSTJ    # reference normals, quad weights

    M::MM                # mass matrix
    Pq::P               # L2 projection matrix

    # specialize diff and lift (dense, sparse, Bern, etc)
    Drst::NTuple{Dim, D} # differentiation operators
    LIFT::L             # lift matrix
end

function Base.show(io::IO, ::MIME"text/plain", rd::RefElemData)
    @nospecialize rd
    print(io,"RefElemData for a degree $(rd.N) $(rd.approximationType) approximation on $(rd.elementType) element.")
end

function Base.show(io::IO, rd::RefElemData)
    @nospecialize basis # reduce precompilation time
    print(io,"RefElemData{N=$(rd.N),$(rd.approximationType),$(rd.elementType)}.")
end

_propertynames(::Type{RefElemData}, private::Bool = false) = (:Nfaces, :Np, :Nq, :Nfq)
function Base.propertynames(x::RefElemData{1}, private::Bool=false) 
    return (fieldnames(RefElemData)..., _propertynames(RefElemData)...,
            :r, :rq, :rf, :rp, :nrJ, :Dr)
end
function Base.propertynames(x::RefElemData{2}, private::Bool = false)
    return (fieldnames(RefElemData)..., _propertynames(RefElemData)...,
            :r, :s, :rq, :sq, :rf, :sf, :rp, :sp, :nrJ, :nsJ, :Dr, :Ds)
end
function Base.propertynames(x::RefElemData{3}, private::Bool = false)
    return (fieldnames(RefElemData)..., _propertynames(RefElemData)...,
            :r, :s, :t, :rq, :sq, :tq, :rf, :sf, :tf, 
            :rp, :sp, :tp, :nrJ, :nsJ, :ntJ, :Dr, :Ds, :Dt)
end

# convenience unpacking routines
function Base.getproperty(x::RefElemData{Dim, ElementType, ApproxType, Nfaces}, s::Symbol) where {Dim, ElementType, ApproxType, Nfaces}
    if s==:r
        return getfield(x, :rst)[1]
    elseif s==:s
        return getfield(x, :rst)[2]
    elseif s==:t
        return getfield(x, :rst)[3]

    elseif s==:rq
        return getfield(x, :rstq)[1]
    elseif s==:sq
        return getfield(x, :rstq)[2]
    elseif s==:tq
        return getfield(x, :rstq)[3]

    elseif s==:rf
        return getfield(x, :rstf)[1]
    elseif s==:sf
        return getfield(x, :rstf)[2]
    elseif s==:tf
        return getfield(x, :rstf)[3]

    elseif s==:rp
        return getfield(x, :rstp)[1]
    elseif s==:sp
        return getfield(x, :rstp)[2]
    elseif s==:tp
        return getfield(x, :rstp)[3]

    elseif s==:nrJ
        return getfield(x, :nrstJ)[1]
    elseif s==:nsJ
        return getfield(x, :nrstJ)[2]
    elseif s==:ntJ
        return getfield(x, :nrstJ)[3]

    elseif s==:Dr
        return getfield(x, :Drst)[1]
    elseif s==:Ds
        return getfield(x, :Drst)[2]
    elseif s==:Dt
        return getfield(x, :Drst)[3]
        
    elseif s==:Nfaces || s==:num_faces
        return Nfaces
    elseif s==:Np
        return length(getfield(x, :rst)[1])
    elseif s==:Nq
        return length(getfield(x, :rstq)[1])
    elseif s==:Nfq
        return length(getfield(x, :rstf)[1])
    elseif s==:elemShape || s==:element_type # for compatibility with Trixi formatting
        return getfield(x, :elementType)
    elseif s==:approximation_type
        return getfield(x, :approximationType)
    else
        return getfield(x, s)
    end
end

"""
    function RefElemData(elem; N, kwargs...)
    function RefElemData(elem, approxType; N, kwargs...)

Keyword argument constructor for RefElemData (to "label" `N` via `rd = RefElemData(Line(), N=3)`)
"""
RefElemData(elem; N, kwargs...) = RefElemData(elem, N; kwargs...)
RefElemData(elem, approxType; N, kwargs...) = RefElemData(elem, approxType, N; kwargs...)

# default to Polynomial-type RefElemData
RefElemData(elem, N::Int; kwargs...) = RefElemData(elem, Polynomial(), N; kwargs...)


@inline Base.ndims(::Line) = 1
@inline Base.ndims(::Union{Tri,Quad}) = 2
@inline Base.ndims(::Union{Tet,Hex}) = 3

@inline face_type(::Union{Tri,Quad}) = Line()
@inline face_type(::Hex) = Quad()
@inline face_type(::Tet) = Tri()

# ====================================================
#          RefElemData approximation types
# ====================================================

struct Polynomial end 

# ========= SBP approximation types ============

struct DefaultSBPType end

# line/quad/hex nodes
struct TensorProductLobatto end

# triangle node types
struct Hicken end 
struct Kubatko{FaceNodeType} end

# face node types for Kubatko
struct LegendreFaceNodes end
struct LobattoFaceNodes end

# SBP ApproximationType: the more common diagonal E diagonal-norm SBP operators on tri/quads.
struct SBP{Type}
    SBP() = new{DefaultSBPType}() # no-parameter default
    SBP{T}() where {T} = new{T}()  # default constructor
end 

# sets default to TensorProductLobatto on Quads 
RefElemData(elem::Union{Line,Quad,Hex}, approxT::SBP{DefaultSBPType}, N) = RefElemData(elem, SBP{TensorProductLobatto}(), N)

# sets default to Kubatko{LobattoFaceNodes} on Tris
RefElemData(elem::Tri, approxT::SBP{DefaultSBPType}, N) = RefElemData(elem, SBP{Kubatko{LobattoFaceNodes}}(), N)

