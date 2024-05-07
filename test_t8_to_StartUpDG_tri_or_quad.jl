using StartUpDG

# computes: VX, VY, levels, neighbor_ids, 
#           dual_faces, orientations
include("generate_t8_arrays.jl")

old_orientations = deepcopy(orientations)

# flip orientation of the second t8 face
# to ensure consistent CCW orientation
for e in eachindex(neighbor_ids)
    if length(neighbor_ids[e]) == 3 # if element "e" is a triangle
        # reverse the orientation of face 2
        @. orientations[e][2] = !Bool(old_orientations[e][2])

        # if the neighbor face isn't face 2 (e.g., also flipped), 
        # reverse orientations of face 2 neighbors.
        for nbr in eachindex(neighbor_ids[e][2])
            enbr = neighbor_ids[e][2][nbr]
            fnbr = dual_faces[e][2][nbr]
            if fnbr != 2
                @. orientations[enbr][fnbr] = !Bool(old_orientations[enbr][fnbr])
            end
        end
        
        # reverse the order in which the neighbors and 
        # neighboring faces are specified. This should
        # only impact split (non-conforming) faces. 
        reverse!(neighbor_ids[e][2])
        reverse!(dual_faces[e][2])
    end
end

# next, convert t8 face ordering to StartUpDG
for e in eachindex(neighbor_ids)  
    if length(neighbor_ids[e]) == 3 # if element "e" is a triangle      
        # Here, p[StUpDG_face] = t8_face 
        p = SVector(3, 1, 2)
        # Here, p_inverse[t8_face] = StUpDG_face
        p_inverse = SVector(2, 3, 1) 
    else # if it's a quad
        p = SVector(1, 2, 3, 4)
        p_inverse = SVector(1, 2, 3, 4)
    end

    permute!(neighbor_ids[e], p)

    permute!(dual_faces[e], p)
    permute!(orientations[e], p)
    for f in eachindex(dual_faces[e])
        if length(dual_faces[e][f]) > 0
            dual_faces[e][f] = p_inverse[dual_faces[e][f]]
        end
    end
end      

num_elements = length(levels)
faces_per_element = length.(neighbor_ids)
face_offsets = cumsum(faces_per_element) .- faces_per_element[1] 
num_faces = sum(faces_per_element)

# compute indices of nonconforming faces
# nonconforming_to_mortar_face = Dict{Int, Vector{Int}}()
nonconforming_to_mortar_face = [[f] for f in 1:num_faces]
mortar_offset = 0
for e in eachindex(neighbor_ids)
    for f in eachindex(neighbor_ids[e])
        if length(neighbor_ids[e][f]) > 1

            # store which faces are split
            face_index = f + face_offsets[e]
            @show face_index
            mortar_index = eachindex(neighbor_ids[e][f]) .+ mortar_offset .+ num_faces
            nonconforming_to_mortar_face[face_index] = mortar_index

            # increment mortar count
            mortar_offset += length(neighbor_ids[e][f]) 
        end
    end
end

# count new mortar faces produced by split faces
nonconforming_faces = findall(length.(nonconforming_to_mortar_face) .> 1)
num_mortar_faces = sum(length.(nonconforming_to_mortar_face[nonconforming_faces]))
mortar_faces = vcat(nonconforming_to_mortar_face[nonconforming_faces]...)

# create face-by-face orientation vectors
face_orientations = ones(Int, num_faces + num_mortar_faces) 
for e in eachindex(orientations) 
    for f in eachindex(orientations[e])
        face_global_index = f + face_offsets[e]
        if length(orientations[e][f]) > 1
            @. face_orientations[nonconforming_to_mortar_face[face_global_index]] = orientations[e][f]
        else
            face_orientations[face_global_index] = isempty(orientations[e][f]) ? 1 : first(orientations[e][f])
        end
    end
end

# determine face indices of mortar neighbors
FToF = collect(1:(num_faces + num_mortar_faces))
for e in eachindex(neighbor_ids)
    for f in eachindex(neighbor_ids[e])
        # if it's not a boundary node
        if length(neighbor_ids[e][f]) > 0
            enbr = neighbor_ids[e][f]
            fnbr = dual_faces[e][f]

            # if face is on a non-conforming interface
            if any(levels[e] .!= levels[enbr])
                if any(levels[e] .< levels[enbr]) # big face
                    big_face = f + face_offsets[e]
                    dual_face_global_indices = @. fnbr + face_offsets[enbr]
                    face_global_indices = nonconforming_to_mortar_face[big_face]                    
                    @. FToF[face_global_indices] = dual_face_global_indices
                else # small face
                    # since it's a small face, (enbr, fnbr) should both be 
                    # one-element vectors with only a `first` entry.
                    dual_big_face = first(dual_faces[e][f] + face_offsets[enbr])
                    nbrs = neighbor_ids[first(enbr)][first(fnbr)]

                    # determine whether e is the first or second neighbor
                    # !!! todo: 3D
                    inbr = nbrs[1] == e ? 1 : 2

                    dual_face_global_indices = nonconforming_to_mortar_face[dual_big_face]
                    face_global_indices = f .+ face_offsets[e]
                    FToF[face_global_indices] = dual_face_global_indices[inbr]
                end
            else # if face is conforming
                dual_face_global_indices = fnbr + face_offsets[enbr]
                face_global_indices = f + face_offsets[e]
                FToF[first(face_global_indices)] = first(dual_face_global_indices)
            end            
        end
    end
end


# function annotate_mesh(etype, VX, VY, neighbor_ids, dual_faces, orientations)

#     # fv = ([2, 3], [1, 3], [1, 2]) # t8 Tri ordering
#     if etype isa Tri
#         fv = ([1, 2], [2, 3], [3, 1]) # StartUpDG Tri ordering
#         p_vertex = SVector(1, 2, 3, 1)
#     elseif etype isa Quad
#         fv = ([1, 3], [2, 4], [1, 2], [3, 4]) # StartUpDG Quad ordering
#         p_vertex = SVector(1, 2, 4, 3, 1)
#     end
    
#     avg(x) = sum(x) / length(x)

#     # t8_to_StartUpDG_face_ordering 
#     plot(size = 1000 .* (1, 1))
#     for e in eachindex(VX, VY)
#         xc, yc = avg(VX[e]), avg(VY[e])
#         annotate!([xc], [yc], [string(e)], markersize=8)
#         plot!(VX[e][p_vertex], VY[e][p_vertex])
#         for v in eachindex(VX[e])
#             xv, yv = VX[e][v], VY[e][v]            
#             annotate!([xc + 0.9 * (xv - xc)], [yc + 0.9 * (yv - yc)], string(v))
#         end
#         for f in eachindex(neighbor_ids[e])
#             global_f = f + face_offsets[e]
#             fids = fv[f]
#             plot!(VX[e][fids], VY[e][fids])
#             xfc = avg(VX[e][fids])
#             yfc = avg(VY[e][fids])
#             annotate!([xc + 0.85 * (xfc - xc)], [yc + 0.85 * (yfc - yc)], string(f) * "(" * string(global_f) * ")")
#         end
#     end
#     display(plot!(leg=false))
# end

# annotate_mesh(etype, VX, VY, neighbor_ids, dual_faces, orientations)

# Create a StartUpDG mesh. 
# etype is defined in "generate_t8_arrays.jl"
rd = RefElemData(etype, Polynomial(), 3)

# construct element nodal coordinates
VX_local = VX
VY_local = VY

(; V1) = rd
x = zeros(size(V1, 1), length(VX_local))
y = zeros(size(V1, 1), length(VX_local))
for e in eachindex(VX_local)
    view(x, :, e) .= V1 * VX_local[e]
    view(y, :, e) .= V1 * VY_local[e]
end

xf, yf = (x -> reshape(rd.Vf * x, rd.Nfq ÷ rd.num_faces, :)).((x, y))

mortar_interpolation_matrix, mortar_projection_matrix = StartUpDG.compute_mortar_operators(rd)

if length(nonconforming_faces) > 0
    xm, ym = (x -> reshape(mortar_interpolation_matrix * x, :, 2 * length(nonconforming_faces))).((xf[:, nonconforming_faces], yf[:, nonconforming_faces]))
    xM, yM = [xf xm], [yf ym]
else
    xM, yM = xf, yf
end

mapM = reshape(1:length(xM), size(xM))
mapP = copy(mapM)
for (f, fnbr) in enumerate(FToF)
    if f != fnbr # if it's not a boundary face        
        @assert face_orientations[f] == face_orientations[fnbr]
        if face_orientations[f] == 1
            @. mapP[end:-1:1, f] = mapM[:, fnbr]
        else
            @. mapP[:, f] = mapM[:, fnbr]
        end
    end
end

# check that all node coordinates match
xy = [[xM[i, j], yM[i, j]] for i in axes(xM, 1), j in axes(xM, 2)]
xydiff = norm.(xy .- xy[mapP])
@show norm(xydiff)

anim = @animate for i in eachindex(xy) 
    scatter(xM, yM, color=:black, markersize=2)
    scatter!(SVector.(xy[i])..., leg=false, marker=:star, markersize=8)
    scatter!(SVector.(xy[mapP[i]])..., leg=false, marker=:star5, markersize=6)
end when i !== mapP[i] 
gif(anim, "check_mapP.gif", fps=4)

rxJ, sxJ, ryJ, syJ, J = StartUpDG.geometric_factors(x, y, rd.Drst...)
rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)
nxJ, nyJ, Jf = StartUpDG.compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)

nxJf = reshape(nxJ, :, num_faces)
nyJf = reshape(nyJ, :, num_faces)
nxJm, nyJm = map(x -> reshape(mortar_interpolation_matrix * x, :, 2 * length(nonconforming_faces)), 
                 (nxJf[:, nonconforming_faces], nyJf[:, nonconforming_faces]))
nxJ = [nxJf nxJm]
nyJ = [nyJf nyJm]
Jf = @. sqrt(nxJ^2 + nyJ^2)

nx, ny = nxJ ./ Jf, nyJ ./ Jf

p = (; mapP, rxJ, sxJ, ryJ, syJ, J, nx, ny, Jf, rd, 
       num_elements, num_faces, mortar_interpolation_matrix, mortar_projection_matrix,
       nonconforming_faces, nonconforming_to_mortar_face, mortar_faces)

function rhs!(du, u, p, t)
    (; mapP, rxJ, sxJ, ryJ, syJ, J, nx, ny, Jf, rd, 
       num_elements, num_faces, mortar_interpolation_matrix, mortar_projection_matrix,
       nonconforming_faces, nonconforming_to_mortar_face, mortar_faces) = p

    # interpolate to local faces
    uf = reshape(rd.Vf * u, :, num_faces)

    # interpolate to mortars 
    um = reshape(mortar_interpolation_matrix * uf[:, nonconforming_faces], :, 2 * length(nonconforming_faces))
    uM = [uf um] # mapM 
    uP = uM[mapP]
    
    # compute fluxes on mortars (and inactive faces) 
    # TODO: keep track of active faces to reduce work
    mortar_fluxes = @. 0.5 * (uP - uM) * nxJ - 0.5 * abs(nx) * (uP - uM) * Jf
    
    projected_mortar_fluxes = 
        mortar_projection_matrix * reshape(mortar_fluxes[:, mortar_faces], :, length(nonconforming_faces))

    fluxes = mortar_fluxes[:, 1:num_faces]
    fluxes[:, nonconforming_faces] .= projected_mortar_fluxes
    
    fluxes = reshape(fluxes, :, num_elements)    
    du .= -(rxJ .* (rd.Dr * u) + sxJ .* (rd.Ds * u) + rd.LIFT * fluxes) ./ J
end

using OrdinaryDiffEq

u = @. exp(-50 * ((x-0.5)^2 + (y-0.5)^2))
tspan = (0, 2.0)
ode = ODEProblem(rhs!, u, tspan, p)

sol = solve(ode, Tsit5(), saveat=LinRange(tspan..., 25))
# scatter(vec(rd.Vp * x), vec(rd.Vp * y), zcolor=vec(rd.Vp * sol.u[end]))
anim = @animate for u in sol.u
    scatter(vec(rd.Vp * x), vec(rd.Vp * y), zcolor=vec(rd.Vp * u), msw=0)
end
gif(anim, "advec.gif")