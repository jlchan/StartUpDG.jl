"RefElemData: contains info on a general 3D reference element.
For 1D and 2D reference elements, some fields may be unused.

Fields:
    fv # list of vertices defining faces, e.g., ([1,2],[2,3],[3,1]) for a triangle

    V1      # low order interpolation matrix
    VDM     # Vandermonde matrix
    Vp      # interpolation matrix to plotting nodes

    r; s; t         # interpolation nodes
    rq; sq; tq      # volume quadrature
    rf; sf; tf;     # surface quadrature
    rp; sp; tp      # plotting nodes

    # quadrature weights
    wq::Array{Float64,1}
    wf::Array{Float64,1}

    nrJ; nsJ; ntJ # reference scaled normals

    # differentiation matrices
    Dr::Array{Float64,2}
    Ds::Array{Float64,2}
    Dt::Array{Float64,2}
    Vq::Array{Float64,2}        # quadrature interpolation matrices
    Vf::Array{Float64,2}
    M::Array{Float64,2}         # mass matrix
    Pq::Array{Float64,2}        # L2 projection matrix
    LIFT::Array{Float64,2}      # quadrature-based lift matrix

Use @unpack to extract fields. Example:

N = 4
rd = init_reference_tri(N)
@unpack r,s = rd

"
mutable struct RefElemData
    # annotate types for all arrays involved in RHS evaluation

    Nfaces # num faces
    fv # list of vertices defining faces, e.g., ([1,2],[2,3],[3,1]) for a triangle

    V1      # low order interpolation matrix
    VDM     # Vandermonde matrix
    Vp      # interpolation matrix to plotting nodes

    r; s; t         # interpolation nodes
    rq; sq; tq      # volume quadrature
    rf; sf; tf;     # surface quadrature
    rp; sp; tp      # plotting nodes

    # quadrature weights
    wq::Array{Float64,1}
    wf::Array{Float64,1}

    nrJ; nsJ; ntJ # reference scaled normals

    # differentiation matrices
    Dr::Array{Float64,2}
    Ds::Array{Float64,2}
    Dt::Array{Float64,2}
    Vq::Array{Float64,2}        # quadrature interpolation matrices
    Vf::Array{Float64,2}
    M::Array{Float64,2}         # mass matrix
    Pq::Array{Float64,2}        # L2 projection matrix
    LIFT::Array{Float64,2}      # quadrature-based lift matrix

    RefElemData() = new() # empty initializer
end

"
mutable struct MeshData: container for info related to a physical mesh
Fields:
    VX; VY; VZ              # vertex coordinates
    K::Int                  # num elems
    EToV                    # mesh vertex array
    FToF::Array{Int64,2}    # face connectivity

    x; y; z                 # physical points
    xf; yf; zf
    xq; yq; zq;             # phys quad points, Jacobian-scaled weights
    wJq::Array{Float64,2}

    # arrays of connectivity indices between face nodes
    mapM
    mapP::Array{Int64,2}
    mapB::Array{Int64,1}

    # volume geofacs
    rxJ::Array{Float64,2}
    sxJ::Array{Float64,2}
    txJ::Array{Float64,2}
    ryJ::Array{Float64,2}
    syJ::Array{Float64,2}
    tyJ::Array{Float64,2}
    rzJ::Array{Float64,2}
    szJ::Array{Float64,2}
    tzJ::Array{Float64,2}
    J::Array{Float64,2}

    # surface geofacs
    nxJ::Array{Float64,2}
    nyJ::Array{Float64,2}
    nzJ::Array{Float64,2}
    sJ::Array{Float64,2}
"
# annotate types for geofacs + connectivity arrays for speed in RHS evals
mutable struct MeshData
    VX; VY; VZ              # vertex coordinates
    K::Int                  # num elems
    EToV                    # mesh vertex array
    FToF::Array{Int64,2}    # face connectivity

    x; y; z                 # physical points
    xf; yf; zf
    xq; yq; zq;             # phys quad points, Jacobian-scaled weights
    wJq::Array{Float64,2}

    # arrays of connectivity indices between face nodes
    mapM
    mapP::Array{Int64,2}
    mapB::Array{Int64,1}

    # volume geofacs
    rxJ::Array{Float64,2}
    sxJ::Array{Float64,2}
    txJ::Array{Float64,2}
    ryJ::Array{Float64,2}
    syJ::Array{Float64,2}
    tyJ::Array{Float64,2}
    rzJ::Array{Float64,2}
    szJ::Array{Float64,2}
    tzJ::Array{Float64,2}
    J::Array{Float64,2}

    # surface geofacs
    nxJ::Array{Float64,2}
    nyJ::Array{Float64,2}
    nzJ::Array{Float64,2}
    sJ::Array{Float64,2}

    MeshData() = new() # empty initializer
end
