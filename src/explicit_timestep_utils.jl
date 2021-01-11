bcopy!(x,y) = x .= y

"""
    ck45()

Returns coefficients rka,rkb,rkc for the 4th order 5-stage low storage Carpenter/Kennedy
Runge Kutta method. Coefficients evolve the residual, solution, and local time, e.g.,
```julia
res = rk4a[INTRK]*res + dt*rhs
@. u += rk4b[INTRK]*res
```
"""
function ck45()
    rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];

    rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0]

    rk4c = [ 0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0 ...
    1.0];
    return rk4a,rk4b,rk4c
end

"""
    dp56()

Dormand-prince 5th order 7 stage Runge-Kutta method (6 function evals via the FSAL property)
Returns Butcher table arrays `A`,`c` and error evolution coefficients `rkE`.

Note there is no `b` array needed due to the FSAL property and because the last stage
is used to compute the error estimator.
"""
function dp56()
    rka = [0.0             0.0             0.0             0.0             0.0             0.0         0.0
            0.2             0.0             0.0             0.0             0.0             0.0         0.0
            3.0/40.0        9.0/40.0        0.0             0.0             0.0             0.0         0.0
            44.0/45.0      -56.0/15.0       32.0/9.0        0.0             0.0             0.0         0.0
            19372.0/6561.0 -25360.0/2187.0  64448.0/6561.0  -212.0/729.0    0.0             0.0         0.0
            9017.0/3168.0  -355.0/33.0      46732.0/5247.0  49.0/176.0      -5103.0/18656.0 0.0         0.0
            35.0/384.0      0.0             500.0/1113.0    125.0/192.0     -2187.0/6784.0  11.0/84.0   0.0 ]

    rkc = vec([0.0 0.2 0.3 0.8 8.0/9.0 1.0 1.0 ])

    # coefficients to evolve error estimator = b1-b2
    rkE = vec([71.0/57600.0  0.0 -71.0/16695.0 71.0/1920.0 -17253.0/339200.0 22.0/525.0 -1.0/40.0 ])

    return rka,rkE,rkc
end

"""
    struct PIparams

Struct containing PI controller parameters.
"""
Base.@kwdef struct PIparams{T}
    order
    errTol::T = 5e-4
    dtmax::T = 1e6
    dtmin::T = 1e-12
end

"""
    compute_adaptive_dt(Q,rhsQrk,dt,rkE,PI::PIparams,prevErrEst=nothing)

returns accept_step (true/false), dt_new, errEst.
uses PI error control method copied from Paranumal library (Warburton et al).

Inputs:
    Q: container of arrays, Q[i] = ith solution field
    rhsQrk: container whose entries are type(Q) for RK rhs evaluations
"""
function compute_adaptive_dt(Q,rhsQrk,dt,rkE,PI::PIparams,prevErrEst=nothing)

    @unpack errTol,order,dtmax,dtmin = PI

    # assemble error estimate using rkE = error est coefficients
    errEstVec = zero.(Q)
    for s = 1:7
        map(bcopy!,errEstVec, @. errEstVec + rkE[s]*rhsQrk[s])
        #bcopy!.(errEstVec, @. errEstVec + rkE[s]*rhsQrk[s])
    end

    # compute scalar error estimate
    errEst = 0.0
    for field = 1:length(Q)
        errEstScale = @. abs(errEstVec[field]) / (errTol*(1+abs(Q[field])))
        errEst += sum(errEstScale.^2) # hairer seminorm
    end
    errEst = sqrt(errEst/sum(length.(Q)))
    accept_step = errEst < 1.0

    # default values taken from Paranumal library by Tim Warburton
    dtnew = .8*dt*(.9/errEst)^(.4/(order+1)) # P controller
    if isnothing(prevErrEst)==false # use PI controller if prevErrEst available
        dtnew *= (prevErrEst/max(1e-14,errEst))^(.3/(order+1))
    end
    return accept_step, max(min(dtmax,dtnew),dtmin), errEst # max/min dt
end
