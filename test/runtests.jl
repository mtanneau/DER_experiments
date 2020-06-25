using LinearAlgebra
using Random
using Printf
using Test

include(joinpath(@__DIR__, "../src/DemandResponse.jl"))
const DR = DemandResponse

import MathOptInterface
const MOI = MathOptInterface

import GLPK
import Linda

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "DER/DER",
    "colgen"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end