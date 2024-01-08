using NodesAndModes: face_basis
using StartUpDG

N=7
rd = RefElemData(Tri(), N)
md_ref = MeshData(uniform_mesh(Tri(), 1)..., rd)
(; x, y) = md_ref
# x = @. x + 0.5 * cos(pi/2 * x) * cos(pi/2 * y)
# y = @. y + 0.5 * cos(pi/2 * x) * cos(pi/2 * y)

x = @. x + 0.25 * cos(pi/2 * x) * sin(pi/2 * y + 1)
y = @. y + 0.25 * cos(pi/2 * x + 2) * cos(pi/2 * y)

# x = x + 0.1 * randn(size(x))
# y = y + 0.1 * randn(size(x))
(; Dr, Ds) = rd
xr, xs = Dr * x, Ds * x
yr, ys = Dr * y, Ds * y

# interpolate to quadrature directly 
xrq, xsq, yrq, ysq = (x -> rd.Vq * x).((xr, xs, yr, ys))
Jq = @. -xsq * yrq + xrq * ysq
wJq = Diagonal(rd.wq) * Jq

md = MeshData(rd, md_ref, x, y)

@show sum(wJq), sum(md.wJq)

Vf_face = face_basis(Tri(), rd.N, rd.rstf...)
norm(Vf_face * (Vf_face \ md.xf) - md.xf)

@show sum(md_ref.wJq), sum(md.wJq)
@show sum(md_ref.wJq .* md_ref.xq), sum(md.wJq .* md.xq)
@show sum(md_ref.wJq .* md_ref.xq.^2), sum(md.wJq .* md.xq.^2)


