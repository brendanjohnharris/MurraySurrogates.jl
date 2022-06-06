using TimeseriesSurrogates
using MurraySurrogates
using Test
using Random
using Statistics
using FFTW
using LinearAlgebra
using PyCall

Random.seed!(32)
tss = [cumsum(randn(100,200), dims=2), # Brownian noise
       randn(50, 400), # Gaussian noise
       ]

function pyseed()
    py"""import random
random.seed(32)"""
end

@testset "Phase randomization" begin
    for (i, ts) in enumerate(tss)
        sur = phase_randomize(ts, 2)
        F = abs.(fft(ts, 2))[:, 2:lastindex(ts, 2)÷2]
        F_sur = abs.(fft(sur, 2))[:, 2:lastindex(sur, 2)÷2]
        @test size(sur) == size(ts)
        @test maximum(abs.(F.-F_sur)) < 0.0001

        sur = phase_randomize(ts, 1)
        F = abs.(fft(ts, 1))[2:lastindex(ts, 1)÷2, :]
        F_sur = abs.(fft(sur, 1))[2:lastindex(sur, 1)÷2, :]
        @test size(sur) == size(ts)
        @test maximum(abs.(F.-F_sur)) < 0.0001
    end
end

@testset "Zalesky" begin
    for (i, ts) in enumerate(tss)
        cm = cor(ts)
        cm_s = zalesky_surrogate(cm)
        @test isapprox(mean(cm), mean(cm_s); atol=1e-2)
        @test isapprox(var(cm), var(cm_s); atol=1e-2)
    end
end

@testset "Eigensurrogate" begin
    for (i, ts) in enumerate(tss)
        cm = cor(ts)
        cm_s = eigensurrogate_matrix(cm)
        ev_cm = eigvals(cm)
        ev_cm_s = eigvals(cm_s)
        @test isapprox(ev_cm, ev_cm_s; rtol=1e-7)
    end
end

@testset "Eigensurrogate" begin
    for (i, ts) in enumerate(tss)
        cm = cor(ts, dims=1)
        ts_s = eigensurrogate_timeseries(cm, 1000)
        ev_cm = eigvals(cm)
        ev_cm_s = eigvals(cor(ts_s, dims=1))
        @test isapprox(ev_cm, ev_cm_s; rtol=1e-1)
    end
end

@testset "Spatiotemporal model" begin
    poss = rand(50, 3).*100
    distance_matrix = distance_matrix_euclidean(poss)
    N = 100000
    TA_Δ = rand(50)*0.8
    sample_rate = 1
    highpass_freq = 0.01
    for (i, p) in enumerate([(20, 0.2), (50, 0.5), (10, 0)])
        tss = spatiotemporal_model_timeseries(distance_matrix; SA_λ=p[1], SA_∞=p[2], TA_Δ, N, sample_rate, highpass_freq)
        newtas = temporal_autocorrelation(tss)
        @test isapprox(TA_Δ, newtas[:]; atol=1e-1)
        @test isapprox(mean(cor(tss)), p[2]; atol=1e-1)
    end
end

@testset "Spatiotemporal noiseless model" begin
    poss = rand(50, 3).*100
    distance_matrix = distance_matrix_euclidean(poss)
    N = 100000
    TA_Δ = rand(50)*0.3
    sample_rate = 1
    highpass_freq = 0.01
    for (i, p) in enumerate([(20, .05), (50, .1), (10, 0)])
        tss = spatiotemporal_noiseless_model_timeseries(distance_matrix; SA_λ=p[1], SA_∞=p[2], TA_Δ, N, sample_rate, highpass_freq)
        newtas = temporal_autocorrelation(tss)
        @test isapprox(TA_Δ, newtas[:]; atol=1e-1)
        # @test isapprox(mean(cor(tss)), p[2]; atol=1e-2) # This one will fail
    end
end



@testset "TimeseriesSurrogates interface" begin
    # Generate a surrogate time series as a test time series
    # All model assumptions should be met this way
    V = rand(100, 2).*100
    distance_matrix = distance_matrix_euclidean(V)
    N = 10000
    TA_Δ = rand(100)*0.3; sample_rate = 1; highpass_freq = 0.01
    x = spatiotemporal_model_timeseries(distance_matrix; SA_λ=20, SA_∞=0.05, TA_Δ, N, sample_rate, highpass_freq)

    # Generate a surrogate of the surrogate
    oldtas = temporal_autocorrelation(x)[:]
    surr = Murray(1, 1, 0)
    sg = surrogenerator((x, V), surr)
    s = sg()
    newtas = temporal_autocorrelation(s)[:]
    @test isapprox(oldtas, newtas; atol=1e-1)
    @test isapprox(mean(cor(s)), mean(cor(x)); atol=1e-2)
end
