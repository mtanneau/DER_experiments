import GLPKMathProgInterface: GLPKSolverMIP

@testset "DERs" begin
    include("fixedLoad.jl")
    include("deferrableLoad.jl")
    include("curtailableLoad.jl")
    include("uninterruptibleLoad.jl")
    include("battery.jl")
    include("thermalLoad.jl")
    include("house.jl")
end