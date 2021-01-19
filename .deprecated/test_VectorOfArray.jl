using StaticArrays
using BenchmarkTools

using RecursiveArrayTools

N = 5
K = 64*64

a,b = 1.0,1.0
Np = (N+1)*(N+2)÷2
Q = SVector{4}([randn(Np,K) for i = 1:4])
rhsQ = SVector{4}([randn(Np,K) for i = 1:4])
function inplace1!(Q,rhsQ,a,b)
    for i = 1:length(Q)
        @. Q[i] = a*Q[i] + b*rhsQ[i]
    end
end
inplace2!(Q,rhsQ,a,b) = Q = @. a*Q + b*rhsQ
inplace3!(Q,rhsQ,a,b) = @. Q = a*Q + b*rhsQ

# bcopy!(x,y) = x.=y

Q2 = VectorOfArray(Q)
rhsQ2 = VectorOfArray(rhsQ)

@btime inplace1!($Q,$rhsQ,$a,$b)   # 208.795 μs (0 allocations: 0 bytes)
@btime inplace2!.($Q,$rhsQ,$a,$b)  # 263.483 μs (8 allocations: 2.63 MiB)
@btime inplace3!.($Q,$rhsQ,$a,$b)  # 216.073 μs (0 allocations: 0 bytes)
@btime inplace3!($Q2,$rhsQ2,$a,$b) # 217.012 μs (0 allocations: 0 bytes)
@btime axpby!.($a, $Q, $b, $rhsQ)  # 113.659 μs (0 allocations: 0 bytes)
