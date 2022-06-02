export get_spectrum_ta,
       how_much_noise,
       make_spectrum,
       spatial_temporal_timeseries,
       spatiotemporal_model_timeseries,
       spatiotemporal_noiseless_model_timeseries,
       ta_to_alpha,
       ta_to_alpha_fast,
       temporal_autocorrelation

get_spectrum_ta(spectrum::AbstractVector)::Float64 = spatiotemporal.get_spectrum_ta(spectrum)

how_much_noise(spectrum::AbstractVector, target_ta::Number)::Float64 = spatiotemporal.how_much_noise(spectrum, target_ta)

make_spectrum(tslen::Int, sample_rate::Number, alpha::Number, highpass_freq::Number)::AbstractVector = spatiotemporal.make_spectrum(tslen, sample_rate, alpha, highpass_freq)

spatial_temporal_timeseries(cm::AbstractMatrix, spectra; dims=1)::AbstractMatrix = (LinearAlgebra.checksquare(cm); st = spatiotemporal.spatial_temporal_timeseries(cm, spectra); dims == 1 ? st' : st)

function spatiotemporal_model_timeseries(distance_matrix::AbstractMatrix;
            SA_λ ::Number,
            SA_∞::Number,
            TA_Δ ::AbstractVector,
            N::Int,
            sample_rate::Number,
            highpass_freq::Number=0.0,
            dims=1)
    LinearAlgebra.checksquare(distance_matrix)
    st = spatiotemporal.spatiotemporal_model_timeseries(distance_matrix::AbstractMatrix, SA_λ, SA_∞, TA_Δ, N, sample_rate, highpass_freq)
    dims == 1 ? st' : st
end

function spatiotemporal_noiseless_model_timeseries(distance_matrix::AbstractMatrix;
            SA_λ::Number,
            SA_∞::Number,
            TA_Δ::AbstractVector,
            N::Int,
            sample_rate::Number,
            highpass_freq::Number=0.0,
            dims=1)
    LinearAlgebra.checksquare(distance_matrix)
    st = spatiotemporal.spatiotemporal_noiseless_model_timeseries(distance_matrix::AbstractMatrix, SA_λ, SA_∞, TA_Δ, N, sample_rate, highpass_freq)
    dims == 1 ? st' : st
end

function ta_to_alpha(tslen::Int, sample_rate::Number, target_ta::Number, highpass_freq=0.0)::Float64
    ta_to_alpha(tslen, sample_rate, highpass_freq, target_ta)
end

function ta_to_alpha_fast(tslen::Int, sample_rate::Number, target_ta::Number, highpass_freq=0.0)::Float64
    ta_to_alpha_fast(tslen, sample_rate, highpass_freq, target_ta)
end
