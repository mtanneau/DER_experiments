module DER

using LinearAlgebra
using Random
using SparseArrays

import MathProgBase
const MPB = MathProgBase

import Linda.Oracle: LindaOracleMIP

"""
    Resource

Abstract Class for Energy Resources.
"""
abstract type Resource end

"""
    LindaOracleMIP

Create a MIP oracle from a Resource
"""
function LindaOracleMIP end


function addmodel! end


include("fixedLoad.jl")
include("curtailableLoad.jl")
include("deferrableLoad.jl")
include("battery.jl")
include("thermalLoad.jl")
include("house.jl")

end  # module