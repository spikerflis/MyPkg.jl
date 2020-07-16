
using Test
using MyPkg

x = rand(100)
y = rand(100)
q = 0.1

@test lls(x, y) isa Tuple{Float64, Float64}
@test bootstrap_lls(x, y) isa Tuple{Vector{Float64,}, Vector{Float64}}

@test mymean(x) isa Float64
@test mystd(x) isa Float64
@test mymedian(x) isa Float64
@test myhist(x) isa Tuple

@test pearson(x, y) isa Float64
@test bootstrap_pearson(x, y, 100) isa Vector{Float64}

@test myIQR(x) isa Float64

@test P2() isa QuantileMethod
@test quantile(x, 0.2, P2()) isa Float64