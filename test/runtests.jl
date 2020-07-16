
using Test
using MyPkg
using Distributions

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


D = MixtureModel(Normal[Normal(1, 0.4), Normal(-3, 0.3)])
pts = rand(D, 100)
gr = -6.0:0.1:6.0

@test mykde(gr, pts, h = 0.5) isa Vector{Float64}