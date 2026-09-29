"""
    RefElemData(elementType::Line, approxType::SBP, N)
    RefElemData(elementType::Quad, approxType::SBP, N)
    RefElemData(elementType::Hex,  approxType::SBP, N)
    RefElemData(elementType::Tri,  approxType::SBP, N)
    RefElemData(elementType::Tet,  approxType::SBP, N)

SBP reference element data for `Line()`, `Quad()`, `Hex()`, `Tri()`, and `Tet()` elements.

For `Line()`, `Quad()`, and `Hex()`, `approxType` is `SBP{TensorProductLobatto}`.

For `Tri()`, `approxType` can be `SBP{Kubatko{LobattoFaceNodes}}`, `SBP{Kubatko{LegendreFaceNodes}}`,
`SBP{Hicken}`, or `SBP{<:SummationByPartsDiagE}`.

For `Tet()`, `approxType` is `SBP{<:SummationByPartsDiagE{LobattoFaceNodes}}`.

`SBP{<:SummationByPartsDiagE}` types require loading SummationByParts.jl (`using SummationByParts`).
"""
function RefElemData(elementType::Line, approxType::SBP{TensorProductLobatto}, N;
                     tol = 100 * eps(), kwargs...)
    rd = RefElemData(elementType, N; quad_rule_vol = gauss_lobatto_quad(0, 0, N), kwargs...)

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Union{Quad, Hex}, approxType::SBP{TensorProductLobatto},
                     N; tol = 100 * eps(), kwargs...)
    rd = RefElemData(elementType,
                     Polynomial(TensorProductQuadrature(gauss_lobatto_quad(0, 0, N))),
                     N; kwargs...)

    # brute-force determine Fmask so that rd.rf = rd.r[rd.Fmask], etc.
    new_Fmask = copy(rd.Fmask)
    for fid in eachindex(rd.rf)
        rstf = SVector(getindex.(rd.rstf, fid))
        for i in eachindex(rd.r)
            rst = SVector(getindex.(rd.rst, i))
            if rstf ≈ rst
                new_Fmask[fid] = i
                break
            end
        end
    end
    rd = @set rd.Fmask = new_Fmask

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Union{Tri, Tet}, approxType::SBP, N; tol = 100 * eps(),
                     kwargs...)
    # `diagE_sbp_nodes` may return a modified `approxType` (e.g., with default parameters resolved)
    quad_rule_vol, quad_rule_face, approxType = diagE_sbp_nodes(elementType, approxType, N)

    # build polynomial reference element using quad rules; will be modified to create SBP RefElemData
    rd = RefElemData(elementType, Polynomial(), N; quad_rule_vol = quad_rule_vol,
                     quad_rule_face = quad_rule_face, kwargs...)

    # determine Fmask = indices of face nodes among volume nodes
    Ef, Fmask = build_Ef_Fmask(rd; tol)

    # Build traditional SBP operators from hybridized operators. See Section 3.2 of
    # [High-order entropy stable dG methods for the SWE](https://arxiv.org/pdf/2005.02516.pdf)
    # by Wu and Chan 2021. [DOI](https://doi.org/10.1016/j.camwa.2020.11.006)
    Qrsth, _ = hybridized_SBP_operators(rd)
    Nq = length(rd.wq)
    Vh_sbp = [I(Nq); Ef]
    Drst = map(Qh -> diagm(1 ./ rd.wq) * (Vh_sbp' * Qh * Vh_sbp), Qrsth)

    rst_sbp = quad_rule_vol[1:(end - 1)] # last entry of `quad_rule_vol` = weights
    rd = @set rd.rst = rst_sbp   # set nodes = SBP nodes
    rd = @set rd.rstq = rst_sbp  # set quad nodes = SBP nodes
    rd = @set rd.Drst = Drst
    rd = @set rd.Fmask = vec(Fmask)

    # TODO: make these more efficient with custom operators?
    rd = @set rd.Vf = droptol!(sparse(Ef), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf))

    # make V1 the interpolation matrix from element vertices to SBP nodal points
    rd = @set rd.V1 = vandermonde(elementType, N, rd.rst...) / rd.VDM * rd.V1

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

#####
##### Utilities for SBP 
#####

# - HDF5 file created using MAT.jl and the following code:
# vars = matread("src/data/sbp_nodes/KubatkoQuadratureRules.mat")
# h5open("src/data/sbp_nodes/KubatkoQuadratureRules.h5", "w") do file
#     for qtype in ("Q_GaussLobatto", "Q_GaussLegendre")
#         group = create_group(file, qtype) # create a group
#         for fieldname in ("Points", "Domain", "Weights")
#             subgroup = create_group(group, fieldname)
#             for N in 1:length(vars[qtype])
#                 subgroup[string(N)] = vars[qtype][N][fieldname]
#             end
#         end
#     end
# end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LobattoFaceNodes}}, N)
    if N == 6
        @warn "N=6 SBP operators with quadrature strength 2N-1 and Lobatto face nodes may require very small timesteps."
    end
    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LobattoFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLobatto"]["Points"][string(N)][]
    r, s = (rs[:, i] for i in 1:size(rs, 2))
    w = vec(vars["Q_GaussLobatto"]["Weights"][string(N)][])
    quad_rule_face = gauss_lobatto_quad(0, 0, N + 1)

    return (r, s, w), quad_rule_face, approxType
end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LegendreFaceNodes}}, N)
    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LegendreFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLegendre"]["Points"][string(N)][]
    r, s = (rs[:, i] for i in 1:size(rs, 2))
    w = vec(vars["Q_GaussLegendre"]["Weights"][string(N)][])
    quad_rule_face = gauss_quad(0, 0, N)

    return (r, s, w), quad_rule_face, approxType
end

parsevec(type, str) = str |>
                      (x -> split(x, ", ")) |>
                      (x -> map(y -> parse(type, y), x))

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Hicken}, N)
    if N > 4
        @error "N > 4 triangular `SBP{Hicken}` operators are not available."
    end

    # from Jason Hicken https://github.com/OptimalDesignLab/SummationByParts.jl/tree/work
    lines = readlines((@__DIR__) * "/data/sbp_nodes/tri_diage_p$N.dat")
    r = parsevec(Float64, lines[11])
    s = parsevec(Float64, lines[12])
    w = parsevec(Float64, lines[13])

    # convert Hicken format to biunit right triangle
    r = @. 2 * r - 1
    s = @. 2 * s - 1
    w = 2.0 * w / sum(w)

    quad_rule_face = gauss_lobatto_quad(0, 0, N + 1)

    return (r, s, w), quad_rule_face, approxType
end

# SBP nodes generated by SummationByParts.jl. More specific methods are defined in the
# `StartUpDGSummationByPartsExt` package extension, which is loaded via `using SummationByParts`.
function diagE_sbp_nodes(elem::Union{Tri, Tet}, approxType::SBP{<:SummationByPartsDiagE}, N)
    error("`SBP{<:SummationByPartsDiagE}` approximation types require SummationByParts.jl. " *
          "Please load it via `using SummationByParts` before constructing the `RefElemData`.")
end

# resolves `quadrature_degree = nothing` to the default `2N-1` and checks that the quadrature
# is accurate enough for a degree N SBP operator.
function _resolve_quadrature_degree(approxType::SBP{<:SummationByPartsDiagE{F}},
                                    N) where {F}
    (; quadrature_degree, tol) = approxType.data
    q = something(quadrature_degree, 2 * N - 1)
    if q < 2 * N - 1
        throw(ArgumentError("`quadrature_degree = $q` is too small for a degree N = $N SBP operator; " *
                            "`quadrature_degree` must be at least 2N-1 = $(2 * N - 1)."))
    end
    return q, SBP(SummationByPartsDiagE{F}(; quadrature_degree = q, tol))
end

# Determines the indices of face nodes among the volume nodes by matching coordinates.
# Returns the extraction matrix `Ef` (with `Ef * u = u[Fmask]`) and `Fmask`.
function build_Ef_Fmask(rd_sbp::RefElemData; tol = 100 * eps())
    (; rstq, rstf, Nfaces) = rd_sbp
    Nfq = length(rstf[1])
    Nfp = Nfq ÷ Nfaces
    rstf = map(x -> reshape(x, Nfp, Nfaces), rstf)
    Fmask = zeros(Int, Nfp, Nfaces)
    Ef = zeros(Nfq, length(rstq[1])) # extraction matrix
    for i in eachindex(rstq[1])
        for f in 1:Nfaces
            distance = sum(((xf, xq),) -> abs.(xf[:, f] .- xq[i]), zip(rstf, rstq))
            id = findall(distance .< tol)
            Fmask[id, f] .= i
            Ef[id .+ (f - 1) * Nfp, i] .= 1
        end
    end
    if any(iszero, Fmask)
        num_unmatched = count(iszero, Fmask)
        error("$num_unmatched face node(s) could not be matched to a volume node within " *
              "tolerance $tol. The face quadrature nodes must coincide with a subset of the " *
              "volume quadrature nodes.")
    end
    return Ef, Fmask
end

get_face_nodes(x::AbstractVector, Fmask) = view(x, Fmask)
get_face_nodes(x::AbstractMatrix, Fmask) = view(x, Fmask, :)

function _convert_RefElemData_fields_to_SBP(rd, approx_type::SBP)
    rd = @set rd.M = Diagonal(rd.wq)
    rd = @set rd.Pq = I
    rd = @set rd.Vq = I
    rd = @set rd.approximation_type = approx_type
    # Vp operator = projects SBP nodal vector onto degree N polynomial, then interpolates to plotting points
    rd = @set rd.Vp = vandermonde(rd.element_type, rd.N, rd.rstp...) /
                      vandermonde(rd.element_type, rd.N, rd.rst...)
    return rd
end

"""
    function hybridized_SBP_operators(rd::RefElemData{DIMS}) 

Constructs hybridized SBP operators given a `RefElemData`. Returns operators `Qrsth..., VhP, Ph`.
"""
function hybridized_SBP_operators(rd)
    (; M, Vq, Pq, Vf, wf, Drst, nrstJ) = rd
    Qrst = (D -> Pq' * M * D * Pq).(Drst)
    Ef = Vf * Pq
    Brst = (nJ -> diagm(wf .* nJ)).(nrstJ)
    Qrsth = ((Q, B) -> 0.5 * [Q-Q' Ef'*B; -B*Ef B]).(Qrst, Brst)
    Vh = [Vq; Vf]
    Ph = M \ transpose(Vh)
    VhP = Vh * Pq
    return Qrsth, VhP, Ph, Vh
end
