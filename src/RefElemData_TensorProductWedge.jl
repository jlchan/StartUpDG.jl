struct TensorProductWedge{TTri <: RefElemData{2, <:Tri}, TLine <: RefElemData{1, <:Line}}
    tri::TTri
    line::TLine
end

# for pretty printing
function _short_typeof(approx_type::TensorProductWedge) 
    T1 = _short_typeof(approx_type.line.approximation_type)
    T2 = _short_typeof(approx_type.tri.approximation_type)
    return "TensorProductWedge{$T1, $T2}"
end

# for clarity that we're taking a tensor product of nodes
_wedge_tensor_product(line, tri) = vec.(meshgrid(line, tri))

function RefElemData(elem::Wedge, approximation_type::TensorProductWedge; kwargs...)

    (; tri, line) = approximation_type

    # Find the vertices of the faces
    fv = face_vertices(elem)

    # create tensor product from interpolation nodes
    t, r = _wedge_tensor_product(line.r, tri.r)
    _, s = _wedge_tensor_product(line.r, tri.s)

    # low order interpolation nodes
    r1, s1, t1 = nodes(elem, 1)
    V1 = vandermonde(elem, 1, r, s, t) / vandermonde(elem, 1, r1, s1, t1)
    
    # isnothing(line.VDM) accounts for `RefElemData` related to Trixi.jl's implementation for SummationByPartsOperators.jl. 
    VDM = isnothing(line.VDM) ? nothing : kron(line.VDM, tri.VDM)    
    Dr  = kron(I(line.Np), tri.Dr)
    Ds  = kron(I(line.Np), tri.Ds)
    Dt  = kron(line.Dr, I(tri.Np))
    Drst = (Dr, Ds, Dt)
    
    # assumes interpolation nodes contain face nodes
    Fmask = find_face_nodes(elem, r, s, t)
    
    # build face quadrature nodes
    rft, sft = map(x->reshape(x, :, 3), tri.rstf)
    tf1, rf1 = vec.(meshgrid(line.rq, view(rft, :, 1)))
    _, sf1   = vec.(meshgrid(line.rq, view(sft, :, 1)))
    tf2, rf2 = vec.(meshgrid(line.rq, view(rft, :, 2)))
    _, sf2   = vec.(meshgrid(line.rq, view(sft, :, 2)))
    tf3, rf3 = vec.(meshgrid(line.rq, view(rft, :, 3)))
    _, sf3   = vec.(meshgrid(line.rq, view(sft, :, 3)))
    rf = vcat(rf1, rf2, rf3, tri.rq, tri.rq)
    sf = vcat(sf1, sf2, sf3, tri.sq, tri.sq)
    tf = vcat(tf1, tf2, tf3, -ones(length(tri.wq)), ones(length(tri.wq)))
    rstf = (rf, sf, tf)

    wft = reshape(tri.wf, :, 3)
    wf1 = (x->x[1] .* x[2])(vec.(meshgrid(line.wq, view(wft, :, 1))))
    wf2 = (x->x[1] .* x[2])(vec.(meshgrid(line.wq, view(wft, :, 2))))
    wf3 = (x->x[1] .* x[2])(vec.(meshgrid(line.wq, view(wft, :, 3))))    
    wf = vcat(wf1, wf2, wf3, tri.wq, tri.wq)

    # index into the face nodes     
    num_line_nodes = length(line.wq)
    num_tri_single_face_nodes = size(wft, 1)
    num_quad_face_nodes = num_line_nodes * num_tri_single_face_nodes
    num_tri_nodes = length(tri.wq)
    quad_face_ids(f) = (1:num_quad_face_nodes) .+ (f-1) * num_quad_face_nodes
    tri_face_ids(f) = (1:num_tri_nodes) .+ (f-1) * num_tri_nodes .+ 3 * num_quad_face_nodes
    node_ids_by_face = (quad_face_ids(1), quad_face_ids(2), quad_face_ids(3), 
                        tri_face_ids(1), tri_face_ids(2))
                        
    # for nrJ and nsJ normal on face 1-3 coincide with the triangular normals
    zt, zq = zeros(num_tri_nodes), zeros(num_quad_face_nodes)
    et, eq = ones(num_tri_nodes), ones(num_quad_face_nodes)
    
    nrJ = [zq; eq; -eq; zt; zt]
    nsJ = [-eq; eq; zq; zt; zt]
    ntJ = [zq; zq; zq; -et; et] 
    
    # Create face interpolation matrix
    if tri.approximation_type isa Polynomial
        Vf_tri, _, _ = basis(tri.element_type, tri.N, rf, sf)
        Vf_tri = Vf_tri / tri.VDM
    elseif tri.approximation_type isa SBP
        tri_Vq = tri.Vq isa UniformScaling ? I(num_tri_nodes) : tri.Vq
        
        # the rows of tri.Vf are ordered by nodes, then by faces. we construct the 
        # interpolation operator to quadrilateral faces by extracting rows corresponding 
        # to each face of the triangle and repeating it line.Nq times 
        num_points_per_tri_face = size(tri.Vf, 1) ÷ 3
        interp_to_quad_faces = 
            vcat(kron(ones(line.Nq), tri.Vf[(1:num_points_per_tri_face), :]), 
                 kron(ones(line.Nq), tri.Vf[num_points_per_tri_face .+ (1:num_points_per_tri_face), :]), 
                 kron(ones(line.Nq), tri.Vf[2 * num_points_per_tri_face .+ (1:num_points_per_tri_face), :]))

        # the interpolation to triangular faces is simply the interpolation to 
        # triangular quadrature nodes, once for the top triangular face and once 
        # for the bottom triangular face. 
        interp_to_tri_faces = vcat(tri_Vq, tri_Vq)
        Vf_tri = vcat(interp_to_quad_faces, interp_to_tri_faces)
    else 
        error("approximation type $(tri.approximation_type) not yet supported")
    end

    if line.approximation_type isa Polynomial
        Vf_line, _ = basis(line.element_type, line.N, tf)
        Vf_line = Vf_line / line.VDM
    elseif line.approximation_type isa SBP
        # 3 = number of quad faces
        line_Vq = line.Vq isa UniformScaling ? I(num_line_nodes) : line.Vq
        interp_to_quad_faces = repeat(kron(line_Vq, ones(tri.Nfq ÷ 3)), 3, 1) 

        interp_to_tri_faces = kron(line.Vf, ones(tri.Nq)) 
        Vf_line = vcat(interp_to_quad_faces, interp_to_tri_faces)
    else 
        error("approximation type $(line.approximation_type) not yet supported")
    end

    Vf = spzeros(length(rf), length(r))
    id = 1
    for j in axes(Vf_line, 2)
        for i in axes(Vf_tri, 2)
            @. Vf[:, id] = Vf_tri[:, i] * Vf_line[:, j]
            id += 1
        end
    end

    if tri.approximation_type isa SBP && line.approximation_type isa SBP
        Fmask = []
        for i in eachindex(rf)
            for j in eachindex(r)
                if Vf[i, j] != 0
                    push!(Fmask, j)
                end
            end
        end
    end

    # create tensor product quadrature rule
    tq, rq  = _wedge_tensor_product(line.rq, tri.rq)
    _,  sq  = _wedge_tensor_product(line.rq, tri.sq)
    wt, wrs = _wedge_tensor_product(line.wq, tri.wq)
    wq = wt .* wrs

    if line.Vq isa UniformScaling && tri.Vq isa UniformScaling
        Vq = I
    else
        Vq = kron(line.Vq isa UniformScaling ? I(num_line_nodes) : line.Vq,
                  tri.Vq isa UniformScaling ? I(num_tri_nodes) : tri.Vq)
    end
    M  = Vq' * diagm(wq) * Vq
    Pq = tri.Pq isa UniformScaling && line.Pq isa UniformScaling ? 
            I : M \ (Vq' * diagm(wq))
    LIFT = M \ (Vf' * diagm(wf))

    if line.Nplot != tri.Nplot
        error("Unequal Nplot values for line and triangle are not currently supported.")
    end

    # tensor product plotting nodes
    tp, rp  = _wedge_tensor_product(line.rp, tri.rp)
    _,  sp  = _wedge_tensor_product(line.rp, tri.sp)
    # `line.Vp` is a `UniformScaling` type for `RefElemData` from SummationByPartsOperators.jl
    Vp = line.Vp isa UniformScaling ? kron(I(num_line_nodes), tri.Vp) : kron(line.Vp, tri.Vp)

    # set the polynomial degree as the tuple of the line and triangle degree for now
    return RefElemData(Wedge(node_ids_by_face), approximation_type, (line.N, tri.N), fv, V1,
                       tuple(r, s, t), VDM, Fmask,
                       tuple(rp, sp, tp), Vp,
                       tuple(rq, sq, tq), wq, Vq,
                       rstf, wf, Vf, tuple(nrJ, nsJ, ntJ),
                       M, Pq, Drst, LIFT, (line.Nplot, tri.Nplot))
end

# TODO: add link to proof when we write it up
inverse_trace_constant(rd::RefElemData{3, <:Wedge, <:TensorProductWedge}) = 
    inverse_trace_constant(rd.approximation_type.line) + inverse_trace_constant(rd.approximation_type.tri)
