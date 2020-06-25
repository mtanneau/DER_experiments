module DemandResponse

const DR = DemandResponse

using LinearAlgebra

include("DER/DER.jl")  # Distributed Energy Resources models
using .DER

# include("aggregator.jl")  # Aggregator

include("oracle.jl")

end # module
