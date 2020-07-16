export mymean, mymedian, mystd, myhist

""" 
    mymean(x)

Compute the mean of `x`.
"""
mymean(x) = sum(x)/length(x)

""" 
    mymedian(x)

Compute the median of `x`.
"""
function mymedian(x)
    n = length(x)
    s = sort(x)
    
    if n % 2 == 1
        s[ceil(Int, n/2)]
    else
        idxs = [floor(Int, (n+1)/2), ceil(Int, (n+1)/2)]
        sum(s[idxs]) / 2
    end
end

"""
    mystd(x)

Compute the corrected sample standard deviation of `x`.
"""
function mystd(x)
    n = length(x)
    sqrt((1/(n - 1) * sum((x .- mymean(x)) .^ 2)))
end

"""
    myhist(x::AbstractVector{T}; 
        bmin = minimum(x), 
        # slightly increase max value to allow inclusion of rightmost point 
        bmax = maximum(x)+(maximum(x)-minimum(x))*0.00001, 
        nbins::Int = floor(Int, length(x)/8))

Compute an equispaced histogram of the values in `x`. 

Keyword arguments `bmin` and `bmax` specify the lower and upper bounds 
for the histogram,while `nbins` specify the number of bins.

Returns the left bin edges, and the counts in each bin (left edge inclusive, 
right edge exclusive).
"""
function myhist end

function myhist(x::AbstractVector; 
        bmin = minimum(x), 
        # slightly increase max to allow inclusion of the rightmost point 
        bmax = maximum(x)+(maximum(x)-minimum(x))*0.00001, 
        nbins::Int = floor(Int, length(x)/8))
    
    s = sort(x)
    bin_counts = zeros(Int, nbins)
    i = 1
    j = 1
    
    # bin size
    bs = (bmax - bmin)/nbins

    while i <= length(s)
        #intvl = Interval(bmin + (j-1)*bs, bmin + (j)*bs, true, false)
        if bmin + (j-1)*bs <= s[i] < bmin + (j)*bs
            bin_counts[j] += 1
            i += 1
        else
            j += 1
        end
        
        if j == nbins + 1; break; end
    end
    left_bin_edges = [bmin + i*bs for i = 0:nbins-1]
    return left_bin_edges, bin_counts
end
