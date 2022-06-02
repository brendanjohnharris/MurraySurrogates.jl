using LinearAlgebra

export temporal_autocorrelation,
       spatial_autocorrelation,
        distance_matrix_euclidean

temporal_autocorrelation(x::AbstractVector) = spatiotemporal.temporal_autocorrelation(x)

temporal_autocorrelation(X::AbstractMatrix; dims=1) = mapslices(temporal_autocorrelation, X; dims)

spatial_autocorrelation(cm::AbstractMatrix, dist::AbstractMatrix, discretization::Int) = (LinearAlgebra.checksquare.((cm, dist)); spatiotemporal.spatial_autocorrelation(cm, dist, discretization))

distance_matrix_euclidean(d; dims=1) = [norm(a-b) for a in eachslice(d; dims), b in eachslice(d; dims)]
