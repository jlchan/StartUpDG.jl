# Time-stepping

For convenience, we include utilities for two explicit Runge-Kutta schemes. For more advanced time-stepping functionality, we recommend [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl). 

## Carpenter and Kennedy's (4,5) method

`ck45()` returns coefficients for a low-storage Runge-Kutta method. Example usage:
```julia
function rk_run(u,rhs,dt,Nsteps)
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

## Dormand-Prince (5,6) with PI error control

`dp56()` returns coefficients for an embedded Runge-Kutta method.

Example usage (TODO: FIX):
```julia

rka,rkE,rkc = dp56()
PI = PIparams(order=5)
Qtmp = similar.(Q)
rhsQrk = ntuple(x->zero.(Q),length(rkE))
prevErrEst = nothing
rhsQ = rhs(Q,rd,md)
bcopy!.(rhsQrk[1],rhsQ) # initialize DOPRI rhs (FSAL property)

while t < T
    for INTRK = 2:7
        k = zero.(Qtmp)
        for s = 1:INTRK-1
            bcopy!.(k, @. k + rka[INTRK,s]*rhsQrk[s])
        end
        bcopy!.(Qtmp, @. Q + dt*k)
        rhsQ = rhs(Qtmp,rd,md)
        bcopy!.(rhsQrk[INTRK],rhsQ)
    end

    global t,dt,i,prevErrEst
    accept_step, dtnew, prevErrEst = compute_adaptive_dt(Q,rhsQrk,dt,rkE,PI,prevErrEst)
    if accept_step
        t += dt
        bcopy!.(Q, Qtmp)
        bcopy!.(rhsQrk[1], rhsQrk[7]) # use FSAL property            
    end
    dt = min(T-t,dtnew)
    i = i + 1  # number of total steps attempted
end
```
