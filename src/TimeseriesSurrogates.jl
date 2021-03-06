import .TimeseriesSurrogates: surrogenerator, Surrogate, SurrogateGenerator
using Statistics
using Random

export Murray

struct Murray <: Surrogate
# `x` should be a tuple containing an TÃN array of measurements and an NÃd matrix of node locations.
    discretization::Int
    sample_rate::Real
    highpass_freq::Real
end

function surrogenerator(x::Tuple{AbstractMatrix, AbstractMatrix}, rf::Murray, rng=Random.default_rng())
    X, V = x
    s = similar(X)
    if any(isnan, X)
        error("Can't do nans yet")
    else
        ð· = distance_matrix_euclidean(V)
        ð² = cor(X, dims=1)
        SA_Î», SA_â = spatial_autocorrelation(ð², ð·, rf.discretization)
        TA_Î = temporal_autocorrelation(X)[:]
    end
    N = size(X, 1)

    init = (ð· = ð·, ð² = ð², SA_Î» = SA_Î», SA_â = SA_â, TA_Î = TA_Î, N = N)
    return SurrogateGenerator(rf, x, s, init, rng)
end

function (sg::SurrogateGenerator{<:Murray})()
    ð·, SA_Î», SA_â, TA_Î, N = getfield.(Ref(sg.init),
        (:ð·, :SA_Î», :SA_â, :TA_Î, :N))
    s, rng, sample_rate, highpass_freq = sg.s, sg.rng, sg.method.sample_rate, sg.method.highpass_freq
    # Find a way to set the python rng
    s .= spatiotemporal_model_timeseries(ð·; SA_Î», SA_â, TA_Î, N, sample_rate, highpass_freq)
    return s
end
