import StatsBase.quantile 
export quantile

import StatsBase.quantile
export QuantileMethod, P2, quantile, myIQR

""" Abstract type for quantile methods. """
abstract type QuantileMethod end

""" 
    P2()

The P-squared algorithm for estimating the quantile of a dataset [^JainChlamtac]. The 
algorithm does this in one sweep, avoiding excessive memory allocation for 
large data sets.

[^JainChlamtac]: Jain, Raj, and Imrich Chlamtac. "The P2 algorithm for dynamic 
    calculation of quantiles and histograms without storing observations." 
    Communications of the ACM 28.10 (1985): 1076-1085.
"""
struct P2 <: QuantileMethod end 

"""
    quantile(x, p::Real, method::QuantileMethod)

Estimate the `p`-quantile of the iterable collection `x` using the provided 
`method`.

## Example 

```julia
x = rand(10000)
quantile(x, 0.2, P2())
```
"""
quantile(x, p::Real, method::QuantileMethod)

""" Update quantile using the quadratic P2 formula. """
function p2_quadratic(d, nprev, ncurr, nnext, qprev, qcurr, qnext)
    # split for clarity
    left = d/(nnext - nprev)
    right = (ncurr - nprev + d)*(qnext - qcurr)/(nnext - ncurr) + 
            (nnext - ncurr - d)*(qcurr - qprev)/(ncurr - nprev)
    return qcurr + left * right 
end

""" Update quantile using the linear P2 formula. """
function p2_linear(d, ncurr, nnext, qcurr, qnext)
    return qcurr + d*(qnext-qcurr)/(nnext - ncurr)
end

function quantile(x, p::Real, method::P2)
    @assert 0.0 <= p <= 1.0
    @assert length(x) > 6
    s = sort(x)
    # qs is an array where
    # q[1] is the minimum of observations so far
    # q[2] current estimate of the p/2 quantile
    # q[3] current estimate of the p quantile
    # q[4] current estimate of the (p+1)/2 quantile
    # q[5] maximum of observations so far
    qs = s[1:5]
    n = [i for i in 1:5]
    ndash = [1, 1+2p, 1+4p, 3+2p, 5]
    dn = [0, p/2, p, (1+p)/2, 1]
    
    for (j, xⱼ) in enumerate(s[6:end])
        if xⱼ < qs[1]
            qs[1] = xⱼ; k = 1
        elseif qs[1] <= xⱼ < qs[2]
            k = 1
        elseif qs[2] <= xⱼ < qs[3]
            k = 2
        elseif qs[3] <= xⱼ < qs[4]
            k = 3
        elseif qs[4] <= xⱼ < qs[5]
            k = 4
        elseif qs[5] < xⱼ 
            qs[5] = xⱼ; k = 4
        else
           error("condition not met!") 
        end
        
        for i = (k+1):5
            n[i] += 1
        end
        
        for i = 1:5
            ndash[i] = ndash[i] + dn[i]
        end
        
        # Adjust quantile estimates.
        for i = 2:4
            dᵢ = ndash[i] - n[i]
            
            if ((dᵢ >= 1) && (n[i+1] - n[i] > 1)) || ((dᵢ <= -1) && (n[i-1] - n[i] < -1))
                dᵢ = sign(dᵢ)
                
                # Use P2 parabolic formula if suited.
                qadj = p2_quadratic(dᵢ, n[i-1], n[i], n[i+1], qs[i-1], qs[i], qs[i+1])
                if qs[i-1] < qadj < qs[i+1]
                    qs[i] = qadj
                else # Otherwise, use linear formula
                    qs[i] = p2_linear(dᵢ, n[i], n[i+1], qs[i], qs[i+1])
                end
                n[i] = n[i] + dᵢ
            end
        end
    end
    
    # Current estimate of the p-quantile
    return qs[3]
end


""" 
    myIQR(x)

Compute the interquartile range (IQR) for `x`
"""
function myIQR(x)
    Q1 = quantile(x, 0.25)
    Q3 = quantile(x, 0.75)
    
    Q3 - Q1
end