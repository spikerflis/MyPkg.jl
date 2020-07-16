using Statistics, StatsBase
export pearson, bootstrap_pearson

"""
    pearson(x, y)

Compute the Pearson correlation coefficient between `x` and `y`.
"""
function pearson(x, y)
    n, m = length(x), length(y)
    @assert n == m
    mX, mY = mean(x), mean(y)
    num = n*sum(x .* y) - sum(x)*sum(y)
    den = sqrt(n*sum(x.^2)-sum(x)^2)*sqrt(n*sum(y.^2)-sum(y)^2)
    r = num/den
end

"""
    bootstrap_pearson(x, y, n_boot::Int = 10000)

Bootstrap `n_boot` estimates of the Pearson correlation coefficient for 
`x` and `y`.
"""
function bootstrap_pearson(x, y, n_boot::Int = 10000)
    N, M = length(x), length(y)
    rs = zeros(Float64, n_boot)
    for i = 1:n_boot
        inds = StatsBase.sample(1:N, N, replace = true)
        xsample = x[inds]
        ysample = y[inds]
        
        rs[i] = pearson(xsample, ysample)
    end
    
    return rs
end