module MurraySurrogates

using PyCall
using Conda

using Requires
const spatiotemporal = PyNULL()
function __init__()
    Conda.pip_interop(true)
    Conda.pip("install", "spatiotemporal")
    copy!(spatiotemporal, pyimport_conda("spatiotemporal", "spatiotemporal"))

    @require TimeseriesSurrogates="c804724b-8c18-5caa-8579-6025a0767c70" include("TimeseriesSurrogates.jl")
end


include("Surrogates.jl")
include("Models.jl")
include("Statistics.jl")
include("util.jl")

end
