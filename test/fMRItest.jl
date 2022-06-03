using MurraySurrogates
using Downloads
using MAT
using Statistics

# ? --------------------------- Download the data -------------------------- ? #
file = Downloads.download("https://unisydneyedu-my.sharepoint.com/:u:/g/personal/bhar9988_uni_sydney_edu_au/EYeMzQvHhuxLiwVly0qH8p4BVXQpW2SsJi9n2vn6yuJLnQ?e=KEXJat&download=1")
_X = matread(file)["sigBPass"]

# ? ------------- Determine parameters for the surrogate model ------------- ? #
"""
The general Murray procedure has four steps:
    1. Use SA_λ, SA_∞ and a spatial distance matrix to generate a spatial correlation matrix
    2. Choose a power spectrum of 1/𝑓²
    3. Generate a time series for each node with spatial correlations from (1) and a power spectrum from (2)
    4. Add white noise to the time series to reduce the temporal autocorrelation (i.e. power spectrum)
However, determining SA_λ and SA_∞ is apparently difficult. We will try it anyway.

See https://spatiotemporal.readthedocs.io/ for usage.
"""
X, nodes = movie2nodes(_X, 1:lastindex(_X, 1), 1:lastindex(_X, 2))
𝐷 = distance_matrix_euclidean(nodes)
𝛲 = cor(X, dims=1) # Correlation matrix
discret = 2 # This is a parameter (the bin size) that you have to manually set.
sub = 1:5:lastindex(𝐷, 1) # Computing SA is really slow, so we only use some of 𝛲 and 𝐷. More points would be more accurate.
SA_λ, SA_∞ = spatial_autocorrelation(𝛲[sub, sub], 𝐷[sub, sub], discret)
TA_Δ = temporal_autocorrelation(X)[sub] # I think this is a dubious quantity. It depends on the sampling period.

# ? ------------------ Generate the surrogate time series ------------------ ? #
sample_rate = 1 # Set this to the fMRI TR
highpass_freq = 0 # Don't apply a high-pass filter

S = spatiotemporal_model_timeseries(𝐷[sub, sub]; SA_λ, SA_∞, TA_Δ, N=size(_X, 3), sample_rate, highpass_freq) # This is also really slow

# ? ---------- Compare the surrogate and the original time series ---------- ? #
new_TA_Δ = temporal_autocorrelation(S)[:]
TA_error = filter(!isnan, (new_TA_Δ - TA_Δ)) |> x->x.^2 |> mean |> sqrt

notnanS = hcat(filter(x->!all(isnan.(x)), eachcol(S)|>collect)...)
SA_error = mean(cor(notnanS)) - SA_∞

S = nodes2movie(_X, S, nodes[sub, :])
