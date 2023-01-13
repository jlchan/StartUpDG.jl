@testset "NamedArrayPartition tests" begin
    x = NamedArrayPartition(a = ones(10), b = rand(20))
    @test typeof(@. sin(x * x^2 / x - 1)) <: NamedArrayPartition
    @test typeof(x.^2) <: NamedArrayPartition
    @test x.a ≈ ones(10)
    @test typeof(x .+ x[1:end]) <: Vector # test broadcast precedence 
    @test all(x .== x[1:end]) 
    y = copy(x)
    @test zero(x, (10, 20)) == zero(x) # test that ignoring dims works
    @test typeof(zero(x)) <: NamedArrayPartition
    @test (y .*= 2).a[1] ≈ 2 # test in-place bcast

    @test length(Array(x))==30
    @test typeof(Array(x)) <: Array
    @test propertynames(x) == (:a, :b)

    x = NamedArrayPartition(a = ones(1), b = 2*ones(1))
    @test Base.summary(x) == "NamedArrayPartition{Float64, RecursiveArrayTools.ArrayPartition{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, NamedTuple{(:a, :b), Tuple{Int64, Int64}}} with arrays:"
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), x)) == "(a = [1.0], b = [2.0])"

    # TODO: get FillArrays broadcast working
    # using FillArrays
    # c = Fill(2., (2, 5))
    # x = NamedArrayPartition(a = ones(10), b = rand(20), c=c)

    using StructArrays, StaticArrays
    x = NamedArrayPartition(a = StructArray{SVector{2, Float64}}((ones(5), 2*ones(5))),
                            b = StructArray{SVector{2, Float64}}((3 * ones(2,2), 4*ones(2,2))))
    @test typeof(x.a) <: StructVector{<:SVector{2}}
    @test typeof(x.b) <: StructArray{<:SVector{2}, 2}
    @test typeof((x->x[1]).(x)) <: NamedArrayPartition
    @test typeof(map(x->x[1], x)) <: NamedArrayPartition
end

# x = NamedArrayPartition(a = ones(10), b = rand(20)) 
# x_ap = ArrayPartition(x)
# @btime @. x_ap * x_ap; #   498.836 ns (5 allocations: 2.77 KiB) 
# @btime @. x * x;  # 2.032 μs (5 allocations: 2.84 KiB) - 5x slower than ArrayPartition
