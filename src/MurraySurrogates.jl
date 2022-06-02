module MurraySurrogates

using PyCall
using Conda

Conda.pip_interop(true)
Conda.pip("install", "spatiotemporal")

const spatiotemporal = PyNULL()
function __init__()
    copy!(spatiotemporal, pyimport_conda("spatiotemporal", "spatiotemporal"))
end


include("Surrogates.jl")
include("Models.jl")
include("Statistics.jl")


end
