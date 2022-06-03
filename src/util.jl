export movie2nodes, nodes2movie

function _movie2nodes(X::AbstractArray{T, 3} where T)
    idxs = Base.compute_itspace(X, Val((1, 2)))
    _X = [X[i...] for i in idxs]
    filter!(x->!all(isnan.(x)), _X)
    return hcat(_X...)
end

function movie2nodes(X::AbstractArray, x::AbstractVector, y::AbstractVector)
    nodes = collect(Iterators.product(x, y))[:]
    idxs = Base.compute_itspace(X, Val((1, 2)))
    _X = [X[i...] for i in idxs]
    nans = .![all(isnan.(x)) for x in _X]
    nodes = nodes[nans] .|> collect
    _X = _X[nans]
    return hcat(_X...), hcat(nodes...)'
end

function nodes2movie(X::AbstractArray, _X, nodes)
    Y = deepcopy(X)
    Y[:] .= NaN
    for (n, i) in nodes |> eachrow |> enumerate
        Y[i..., :] .= _X[:, n]
    end
    return Y
end
