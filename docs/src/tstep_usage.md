# Time-stepping

For convenience, we include utilities for two explicit Runge-Kutta schemes. For more advanced time-stepping functionality, we recommend [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).

## Carpenter and Kennedy's (4,5) method

`ck45()` returns coefficients for a low-storage Runge-Kutta method.

### Example usage:
```julia
using Plots
using StartUpDG

# Brusselator
B = 3
f(y1,y2) = [1+y1^2*y2-(B+1)*y1, B*y1-y1^2*y2]
f(Q) = f(Q[1],Q[2])
p,q = 1.5, 3.0
Q = [p;q]

dt = .1
FinalTime = 20

res = zero.(Q) # init RK residual
rk4a,rk4b,rk4c = ck45()
Nsteps = ceil(FinalTime/dt)
dt = FinalTime/Nsteps

plot() # init plot
for i = 1:Nsteps
    global res # yes, I know...this is just for simplicty
    for INTRK=1:5
        rhsQ = f(Q)
        @. res = rk4a[INTRK]*res + dt*rhsQ # i = RK stage
        @. Q =  Q + rk4b[INTRK]*res
    end
    scatter!([i*dt;i*dt],Q)
end
display(plot!(leg=false))
```

## Dormand-Prince (5,6) with PI error control

`dp56()` returns coefficients for an embedded Runge-Kutta method.

### Example usage:
```julia
using Plots
using StartUpDG

# Brusselator
B = 3
f(y1,y2) = [1+y1^2*y2-(B+1)*y1, B*y1-y1^2*y2]
f(Q) = f(Q[1],Q[2])
p,q = 1.5, 3.0
Q = [p;q]

dt = .01
FinalTime = 20

rka,rkE,rkc = dp56()
PI = PIparams(order=5,errTol=1e-6)

Qtmp = copy(Q) # temp storage for RK stage vecs y_i
rhsQrk = (f(Q),ntuple(x->zero(Q),length(rkE)-1)...) # RK f(y_i), set first elem (FSAL)
prevErrEst = nothing
t = 0.0

plot()
while t < FinalTime
    global Q,rhsQrk,t,dt,prevErrEst
    for INTRK = 2:7
        fill!(Qtmp,zero(eltype(first(Q))))
        for s = 1:INTRK-1
            @. Qtmp = Qtmp + rka[INTRK,s]*rhsQrk[s]
        end
        @. Qtmp = Q + dt*Qtmp
        rhsQrk[INTRK] .= f(Qtmp)
    end
    # Hairer seminorm error estimator
    errEstVec = sum(rkE.*rhsQrk) ./ (PI.errTol * (@. 1 + abs(Q)))
    errEst    = mapreduce(x->x*x,+,errEstVec) / sum(length.(Q))

    accept_step,dtnew,prevErrEst = compute_adaptive_dt(errEst,dt,PI,prevErrEst)
    if accept_step
        t += dt
        @. Q = Qtmp
        @. rhsQrk[1] = rhsQrk[7] # use FSAL property

        scatter!([t;t],Q)
    end
    dt = min(FinalTime-t,dtnew)
end
display(plot!(leg=false))
```
