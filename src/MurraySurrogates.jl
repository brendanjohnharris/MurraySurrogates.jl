module MurraySurrogates

using PyCall
using Conda

const spatiotemporal = PyNULL()
function __init__()
    Conda.pip_interop(true)
    Conda.pip("install", "spatiotemporal")
    copy!(spatiotemporal, pyimport_conda("spatiotemporal", "spatiotemporal"))
end


include("Surrogates.jl")
include("Models.jl")
include("Statistics.jl")
include("util.jl")

end
