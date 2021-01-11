# Example: computing DG derivatives

`MeshData` can be used to compute DG derivatives. Suppose ``f`` is a differentiable function and the domain ``\Omega`` can be decomposed into non-overlapping elements ``D^k``. The approximation of ``\frac{\partial f}{\partial x}`` can be approximated using the following formulation: find piecewise polynomial ``u`` such that for all piecewise polynomials ``v``
```math
\int_{\Omega} u v = \sum_k \left( \int_{D^k} \frac{\partial u}{\partial x}v + \int_{\partial D^k} \frac{1}{2} \left[u\right]n_x v \right)
```
Here, ``\left[u\right] = u^+ - u`` denotes the jump across an element interface, and ``n_x`` is the ``x``-component of the outward unit normal on ``D^k``.

Discretizing the left-hand side of this formulation yields a mass matrix. Inverting this mass matrix to the right hand side yields the DG derivative. We show how to compute it for a uniform triangular mesh using `MeshData` and `StartUpDG.jl`.

We first construct the triangular mesh and initialize `md::MeshData`.
```julia
using StartUpDG
using Plots

N = 3
K1D = 8
rd = RefElemData(Tri(),N)
VX,VY,EToV = uniform_mesh(Tri(),K1D)
md = MeshData(VX,VY,EToV,rd)
```
We can approximate a function using interpolation
```julia
f(x,y) = exp(-5*(x^2+y^2))*sin(1+pi*x)*sin(2+pi*y)
@unpack x,y = md
u = @. f(x,y)
```
or using quadrature-based projection
```julia
@unpack Pq = rd
@unpack x,y,xq,yq = md
u = Pq*f.(xq,yq)
```
To visualize the resulting approximation, we can use using `scatter` in Plots.jl. There are also more efficient approaches, such as outputting the result to VTK or using `Triplot.jl`. 
```julia
@unpack Vp = rd
xp,yp,up = Vp*x,Vp*y,Vp*u # interp to plotting points
plot_opt = (msw=0,leg=false,ratio=1,cam=(0,90))
scatter(xp,yp,up,zcolor=up,plot_opt...)
```
Given the nodal values of the polynomial approximation `u`, we can compute its DG derivative as follows
```julia
function dg_deriv_x(u,md::MeshData,rd::RefElemData)
  @unpack Vf,Dr,Ds,LIFT = rd
  @unpack rxJ,sxJ,J,nxJ,mapP = md
  uf = Vf*u
  ujump = uf[mapP]-uf

  # local derivatives using chain rule + lifted flux contributions
  ux = rxJ.*(Dr*u) + rxJ.*(Dr*u)  
  dudxJ = ux + .5*LIFT*(ujump.*nxJ)

  return dudxJ./J
end
```
We can visualize the result as follows:
```julia
dudx = dg_deriv_x(u,md,rd)
uxp = Vp*dudx
scatter(xp,yp,uxp,zcolor=uxp,plot_opt...)
```

![alt text](src/figs/u.png "u")
![alt text](src/figs/dudx.png "dudx")
