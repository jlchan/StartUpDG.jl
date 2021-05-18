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
PI = PIparams(order=5,errTol=1e-5)

# store RK f(y_i), init first entry via FSAL
rhsQrk = (f(Q),ntuple(x->zero(Q),length(rkE)-1)...)
Qtmp = copy(Q) # scratch storage for RK stage vecs y_i
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
