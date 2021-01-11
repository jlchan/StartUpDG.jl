# Time-stepping

For convenience, we include utilities for two simple Runge-Kutta schemes. For more advanced functionality, we recommend the interface to [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).

## Carpenter and Kennedy's (4,5) method

`ck45()` returns coefficients for a low-storage Runge-Kutta method. Example usage:
```julia
function rk_run(u,rhs,Nsteps)
  rk4a,rk4b,rk4c = ck45()
  for i = 1:Nsteps  
    for INTRK = 1:5    
      res = rk4a[INTRK]*res + dt*rhs(u) # i = RK stage
      @. u += rk4b[INTRK]*res
    end
  end
  return u
end
```

## PI error control with Dormand-Prince (5,6)

Todo: finish, add example.
