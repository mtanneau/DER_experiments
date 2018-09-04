using DemandResponse
const DR = DemandResponse

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end



# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "DER",
    "aggregator"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end
