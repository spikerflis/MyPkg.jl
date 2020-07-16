export lls, bootstrap_lls

"""
    lls(x, y)

Linear least squares regression between `x` and `y`. Returns 
a tuple (intercept, slope).
"""
function lls(x, y)
    n, m = length(x), length(y)
    @assert n == m
    
    mX = Statistics.mean(x)
    mY = Statistics.mean(y)
    
    X = sum(x)
    Y = sum(y)
    X2 = sum(x.^2)
    Y2 = sum(x.^2)
    XY = sum(x .* y)
    
    intcpt = (mY*X2 - mX*XY) / (X2 - n*mX^2)
    slope = (XY - n*mX*mY) / (X2 - n*mX^2)
    
    intcpt, slope
end

"""
    bootstrap_lls(x, y, n_boot::Int = 1000)

Compute `n_boot` different bootstrap estimates 
for the linear regression line between `x` and `y`.

Returns a tuple  `(intercepts, slopes)` where 
`intercepts` and `slopes` both are column vectors.
"""
function bootstrap_lls(x, y, n_boot::Int = 1000)
    N, M = length(x), length(y)
    @assert N == M
    
    intcpts = zeros(Float64, n_boot)
    slopes = zeros(Float64, n_boot)
    
    for i = 1:n_boot
        idxs = sample(1:N, N, replace = true)
        xsample = x[idxs]
        ysample = y[idxs]
        intcpts[i], slopes[i] = lls(xsample, ysample)
    end

    return intcpts, slopes
end