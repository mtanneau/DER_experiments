import GLPK

@testset "DERs" begin
    include("fixedLoad.jl")
    include("deferrableLoad.jl")
    include("curtailableLoad.jl")
    include("shiftableLoad.jl")
    include("battery.jl")
    include("thermalLoad.jl")
end