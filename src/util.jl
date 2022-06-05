using ScatteredInterpolation
using ProgressLogging

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
    # Assume the NaNs are consistent over the third dimension, and interpolate over the good channels
    idxs = CartesianIndices(Y[:, :, 1])[.!isnan.(Y[:, :, 1])]
    Y[:] .= NaN

    threadlog = 0
    @withprogress for (i, t) in enumerate(1:lastindex(Y, 3))
        is = .!isnan.(_X[i, :])
        itp = ScatteredInterpolation.interpolate(Multiquadratic(), nodes[is, :]', _X[i, is])
        out = ScatteredInterpolation.evaluate(itp, hcat(collect.(Tuple.(idxs))...))
        Y[idxs, t] .= out
        (threadlog += 1)%1000 == 0 && @logprogress threadlog/size(Y, 3)
    end
    return Y
end
