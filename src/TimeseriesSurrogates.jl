import TimeseriesSurrogates: surrogenerator, Surrogate
using MurraySurrogates
using Statistics

struct Murray <: Surrogate
# `x` should be a tuple containing an N-vector of measurements and an N×d matrix of node locations.
    discretization::Int
    sample_rate::Real
    highpass_freq::Real
end

function surrogenerator(X::Tuple{AbstractVector, AbstractMatrix}, rf::Murray, rng=Random.default_rng())
    x, V = X
    s = similar(x)
    if any(isnan, x)
        error("Can't do nans yet")
    else
        𝐷 = distance_matrix_euclidean(V)
        𝛲 = cor(X, dims=1)
        SA_λ, SA_∞ = spatial_autocorrelation(𝛲, 𝐷, rf.discretization)
        TA_Δ = temporal_autocorrelation(X)
    end
    N = size(length(x))

    init = (𝐷 = 𝐷, 𝛲 = 𝛲, SA_λ = SA_λ, SA_∞ = SA_∞, TA_Δ = TA_Δ, N = N)
    return SurrogateGenerator(rf, x, s, init, rng)
end

function (sg::SurrogateGenerator{<:RPNFT})()
    𝐷, SA_λ, SA_∞, TA_Δ, N = getfield.(Ref(sg.init),
        (:𝐷, :SA_λ, :SA_∞, :TA_Δ, :N))
    s, rng, sample_rate, highpass_freq = sg.s, sg.rng, sg.sample_rate, sg.highpass_freq
    # Find a way to set the python rng
    s .= spatiotemporal_model_timeseries(𝐷; SA_λ, SA_∞, TA_Δ, N, sample_rate, highpass_freq)
    return s
end
