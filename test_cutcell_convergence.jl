using Plots
using PathIntersections
using StartUpDG

function compute_L2_error(N, cells_per_dimension; 
                          u_exact = (x,y) -> sin(pi * x) * sin(pi * y),
                          use_srd=false, use_interp=false)
    rd = RefElemData(Quad(), N, quad_rule_vol=quad_nodes(Quad(), N+1))    
    objects = (PresetGeometries.Circle(R=0.3),)
    md = MeshData(rd, objects, cells_per_dimension, cells_per_dimension; precompute_operators=true)

    srd = StateRedistribution(rd, md)

    (; physical_frame_elements, cut_face_nodes) = md.mesh_type
    Vq_cut, Pq_cut, M_cut = ntuple(_ -> Matrix{Float64}[], 3)
    for (e, elem) in enumerate(physical_frame_elements)

        VDM = vandermonde(elem, rd.N, md.x.cut[:, e], md.y.cut[:, e]) # TODO: should these be md.x, md.y?
        Vq, _ = map(A -> A / VDM, basis(elem, rd.N, md.xq.cut[:,e], md.yq.cut[:, e]))

        M = Vq' * diagm(md.wJq.cut[:, e]) * Vq        
        push!(M_cut, M)
        push!(Vq_cut, Vq)
        push!(Pq_cut, M \ (Vq' * diagm(md.wJq.cut[:, e])))
    end

    if use_interp==true 
        u = u_exact.(md.xyz...)
    else # projection
        uq = u_exact.(md.xyzq...)
        u = similar(md.x)
        u.cartesian .= rd.Pq * uq.cartesian
        for e in eachindex(physical_frame_elements)
            u.cut[:, e] .= Pq_cut[e] * uq.cut[:, e]
        end
    end

    if use_srd == true
        srd(u)
    end

    # eval solution at quad points
    uq = similar(md.xq)
    uq.cartesian .= rd.Vq * u.cartesian
    for e in eachindex(physical_frame_elements)
        uq.cut[:, e] .= Vq_cut[e] * u.cut[:, e]
    end

    L2err = sqrt(abs(sum(md.wJq .* (uq - u_exact.(md.xyzq...)).^2)))
    return L2err
end

N = 3
num_cells = [4, 8, 16, 32, 64]

L2_error, L2_error_srd, L2_error_interp, L2_error_interp_srd  = ntuple(_ -> Float64[], 4)
for cells_per_dimension in num_cells
    @show cells_per_dimension
    use_interp = true
    push!(L2_error_interp, compute_L2_error(N, cells_per_dimension; use_interp))
    push!(L2_error_interp_srd, compute_L2_error(N, cells_per_dimension; use_interp, use_srd=true))
    use_interp = false
    push!(L2_error, compute_L2_error(N, cells_per_dimension; use_interp))
    push!(L2_error_srd, compute_L2_error(N, cells_per_dimension; use_interp, use_srd=true))
end

h = 2 ./ num_cells
plot()
plot!(h, L2_error_interp, marker=:dot, label="Interp")
plot!(h, L2_error_interp_srd, marker=:dot, label="Interp (SRD)")
plot!(h, L2_error, marker=:dot, label="Projection")
plot!(h, L2_error_srd, marker=:dot, label="Projection (SRD)")
plot!(h, 1e-1*h.^(N+1), linestyle=:dash, label="h^{N+1}")
plot!(xaxis=:log, yaxis=:log, legend = :topleft)