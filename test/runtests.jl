using Pkg
if !haskey(Pkg.installed(), "Linda")
    Pkg.clone("https://github.com/ds4dm/Linda.jl")
end

using LinearAlgebra
using Random
using Printf
using Test

using DemandResponse
const DR = DemandResponse

import MathProgBase
const MPB = MathProgBase

import GLPKMathProgInterface: GLPKSolverMIP
import Linda

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "DER/DER",
    "aggregator"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end