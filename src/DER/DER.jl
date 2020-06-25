module DER

import MathOptInterface
const MOI = MathOptInterface

"""
    Resource

Abstract Class for Energy Resources.
"""
abstract type Resource end

function add_resource_to_model! end

include("house.jl")

include("fixedLoad.jl")
include("curtailableLoad.jl")
include("deferrableLoad.jl")
include("battery.jl")
include("thermalLoad.jl")
include("shiftableLoad.jl")

end  # module