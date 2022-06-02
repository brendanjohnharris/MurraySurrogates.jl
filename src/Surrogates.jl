using LinearAlgebra
using Statistics

export eigensurrogate_matrix,
       eigensurrogate_timeseries,
       phase_randomize,
       zalesky_surrogate

eigensurrogate_matrix(cm::AbstractMatrix) = (LinearAlgebra.checksquare(cm); spatiotemporal.eigensurrogate_matrix(cm))

eigensurrogate_timeseries(cm::AbstractMatrix, N_timepoints; dims=1) = (LinearAlgebra.checksquare(cm); st = spatiotemporal.eigensurrogate_timeseries(cm, N_timepoints); dims == 1 ? st' : st)
eigensurrogate_timeseries(ts; dims=1) = (cm = cor(ts; dims); eigensurrogate_timeseries(cm, size(ts, dims); dims))

# Dim correpsonds to the dmension containing time, i.e. the dimension to perform the surrogate transformation on.
function phase_randomize(tss::AbstractArray, dim::Int=1)
    f = spatiotemporal.phase_randomize
    if dim == 1
        f(tss')'
    elseif dim == 2
        f(tss)
    else
        mapslices(f, tss; dims=1)
    end
end

phase_randomize(tss::AbstractVector, dim::Int=1) = phase_randomize(tss[:, :])

zalesky_surrogate(cm::AbstractMatrix) = (LinearAlgebra.checksquare(cm); spatiotemporal.zalesky_surrogate(cm))
