using Distances, NearestNeighbors
export mykde

"""
    mykde(gr, pts::AbstractVector, h = 0.1, metric::Metric = Chebyshev(), normalize = true)

Compute a kernel density estimate to the distribution of `pts` over the grid `gr` using 
a box kernel with bandwith `h`.
"""
function mykde(gr, pts::AbstractVector; h = 0.1, metric::Metric = Chebyshev(), normalize = false)
    tree = KDTree(transpose(pts), metric)
    n = length(pts)
    fx = zeros(Float64, length(gr))
    
    # Box scaling
    s = 1/(2*n*h)
    
    for i = 1:length(gr)
        fX = (inrange(tree, [gr[i]], h, false) |> length) - 1
        fx[i] = s .* fX
    end
    
    if normalize
        fx .= fx ./ sum(fx)
    end
    
    fx
end