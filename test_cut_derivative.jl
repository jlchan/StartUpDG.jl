using StartUpDG

cells_per_dimension = 32
circle = PresetGeometries.Circle(R=0.66, x0=0, y0=0)

rd = RefElemData(Quad(), N=3)
md = MeshData(rd, (circle, ), cells_per_dimension; precompute_operators=true)

(; differentiation_matrices, lift_matrices, face_interpolation_matrices) = 
    md.mesh_type.cut_cell_operators

(; x, y) = md

u_exact(x, y) = sin(pi * x) * sin(pi * y)
dudx_exact(x, y) = pi * cos(pi * x) * sin(pi * y)
dudy_exact(x, y) = pi * sin(pi * x) * cos(pi * y)

(; N) = rd
u_exact(x, y) = x^N + y^N
dudx_exact(x, y) = N * x^(N-1)
dudy_exact(x, y) = N * y^(N-1)

u = u_exact.(x, y)
(; physical_frame_elements, cut_face_nodes) = md.mesh_type

uf = similar(md.xf)
uf.cartesian = rd.Vf * u.cartesian
for e in eachindex(face_interpolation_matrices)
    ids = cut_face_nodes[e]
    Vf = face_interpolation_matrices[e]
    uf.cut[ids] = Vf * u.cut[:, e]
end

uP = vec(uf[md.mapP])
flux = @. 0.5 * (uP - uf) 

dudx, dudy = similar(md.x), similar(md.x)     
dudx.cartesian .= (md.rxJ.cartesian .* (rd.Dr * u.cartesian)) ./ md.J
dudy.cartesian .= (md.syJ.cartesian .* (rd.Ds * u.cartesian)) ./ md.J
for (e, elem) in enumerate(physical_frame_elements)
    Dx, Dy = differentiation_matrices[e]
    LIFT = lift_matrices[e]
    ids = cut_face_nodes[e]
    dudx.cut[:, e] .= Dx * u.cut[:,e] + LIFT * (flux[ids] .* md.nxJ.cut[ids])
    dudy.cut[:, e] .= Dy * u.cut[:,e] + LIFT * (flux[ids] .* md.nyJ.cut[ids])
end

@show norm(dudx - dudx_exact.(x,y), Inf)
@show norm(dudy - dudy_exact.(x,y), Inf)

# scatter(md.xyz..., dudx - dudx_exact.(x,y), zcolor=dudx - dudx_exact.(x,y), leg=false)