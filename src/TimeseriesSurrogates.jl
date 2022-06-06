import .TimeseriesSurrogates: surrogenerator, Surrogate, SurrogateGenerator
using Statistics
using Random

export Murray

struct Murray <: Surrogate
# `x` should be a tuple containing an T×N array of measurements and an N×d matrix of node locations.
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
        𝐷 = distance_matrix_euclidean(V)
        𝛲 = cor(X, dims=1)
        SA_λ, SA_∞ = spatial_autocorrelation(𝛲, 𝐷, rf.discretization)
        TA_Δ = temporal_autocorrelation(X)[:]
    end
    N = size(X, 1)

    init = (𝐷 = 𝐷, 𝛲 = 𝛲, SA_λ = SA_λ, SA_∞ = SA_∞, TA_Δ = TA_Δ, N = N)
    return SurrogateGenerator(rf, x, s, init, rng)
end

function (sg::SurrogateGenerator{<:Murray})()
    𝐷, SA_λ, SA_∞, TA_Δ, N = getfield.(Ref(sg.init),
        (:𝐷, :SA_λ, :SA_∞, :TA_Δ, :N))
    s, rng, sample_rate, highpass_freq = sg.s, sg.rng, sg.method.sample_rate, sg.method.highpass_freq
    # Find a way to set the python rng
    s .= spatiotemporal_model_timeseries(𝐷; SA_λ, SA_∞, TA_Δ, N, sample_rate, highpass_freq)
    return s
end
