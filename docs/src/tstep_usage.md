# Time-stepping

For convenience, we include a commonly used explicit Runge-Kutta scheme. For more advanced time-stepping functionality, we recommend [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl). 

## Carpenter and Kennedy's (4,5) method

`ck45()` returns coefficients for a low-storage Runge-Kutta method.

### Example usage:
```julia
using Plots
using StartUpDG

# Brusselator
B = 3
f(y1, y2) = [1 + y1^2 * y2 - (B+1) * y1, B * y1 - y1^2 * y2]
f(Q) = f(Q[1], Q[2])
p,q = 1.5, 3.0
Q = [p; q]

dt = .1
FinalTime = 20

res = zero.(Q) # init RK residual
rk4a, rk4b, rk4c = ck45()
Nsteps = ceil(FinalTime / dt)
dt = FinalTime / Nsteps

plot() # init plot
for i = 1:Nsteps
    global res # yes, I know...this is just for simplicty
    for INTRK = 1:5
        rhsQ = f(Q)
        @. res = rk4a[INTRK] * res + dt * rhsQ # i = RK stage
        @. Q =  Q + rk4b[INTRK] * res
    end
    scatter!([i*dt;i*dt],Q)
end
display(plot!(leg=false))
```