# Reference elements

The struct `rd::RefElemData` contains data for a given element type. Currently, four types of reference elements are supported: `Line`, `Tri`, `Quad`, and `Hex`.

To initalize a `RefElemData`, just specify the element type and polynomial degree.
```julia
N = 3
rd = RefElemData(Line(),N)
rd = RefElemData(Tri(),N)
rd = RefElemData(Quad(),N)
rd = RefElemData(Hex(),N)
```

## Specifying different quadrature rules.

By default, `RefElemData` initializes volume and surface quadrature rules to be the minimum rules which exactly integrate the unweighted volume and surface mass matrices. If
```julia
N = 3

# create degree N tensor product Gauss-Lobatto rule
r1D,w1D = gauss_lobatto_quad(0,0,N)
rq,sq = vec.(StartUpDG.meshgrid(r1D))
wr,ws = vec.(StartUpDG.meshgrid(w1D))
wq = @. wr*ws

rd = RefElemData(Quad(),N; quad_rule_vol =(rq,sq,wq),  
                           quad_rule_face=(r1D,w1D))
```
This results in a DG spectral element method (DG-SEM) discretization, with a diagonal lumped mass matrix and differentiation matrices which satisfy a summation-by-parts property. 
