using QPSReader
using Logging
using Printf

const RMPDIR = joinpath(@__DIR__, "rmp")
const RMPS = readdir(RMPDIR)

CG = Dict{Tuple{Int,Int}, Int}()
res = Dict{Tuple{Int,Int}, Dict}()

# Go through each RMP file to find number of CG iterations
for frmp in RMPS
    # File name has form DER_<T>_<R>_<iter>.mps
    fname, T, R, n = split(frmp, r"_|\.")[1:4]
    T = parse(Int, T)
    R = parse(Int, R)
    niter = parse(Int, n)

    # Check if such file already exists
    n_ = get(CG, (T, R), 0)
    CG[T, R] = max(niter, n_)
end

function parse_cg_log(flog)
    tmp = NaN
    tsp = NaN
    open(flog, "r") do f
        for ln in eachline(f)
            if occursin("Total time / MP", ln)
                tmp = parse(Float64, split(ln)[end][1:end-1])
                continue
            end
            if occursin("Total time / SP", ln)
                tsp = parse(Float64, split(ln)[end][1:end-1])
                continue
            end
        end
    end
    return tmp, tsp
end

# Go through the files again to access number of columns and linking constraints
for ((T, R), niter) in CG
    r_ = Dict()
    r_["cgiter"] = niter

    # Read 
    frmp = "DER_$(T)_$(R)_$(niter).mps"
    qps = with_logger(Logging.NullLogger()) do
        readqps(joinpath(RMPDIR, frmp))
    end
    
    # Problem stats:
    #   * linking constraints
    #   * linking variables
    #   * Number of sub-problems
    #   * Density of lower blocks
    M = qps.ncon
    N = qps.nvar
    m0 = M - R
    
    nz = 0
    blockidx = zeros(Int, N)
    for (i, j, v) in zip(qps.arows, qps.acols, qps.avals)
        if 1 <= i <= m0
            nz += 1 
        else
            @assert blockidx[j] == 0
            blockidx[j] = i - m0
        end
    end
    n0 = sum(blockidx .== 0)

    r_["m0"] = m0
    r_["n0"] = n0
    r_["R"] = R
    r_["N"] = N
    r_["density"] = nz / (m0 * N)

    # Parse log to get MP and SP total time
    flog = "cg_$(T)_$(R).log"
    tmp, tsp = parse_cg_log(joinpath(@__DIR__, "log", flog))

    r_["tmp"] = tmp
    r_["tsp"] = tsp

    @info frmp r_

    res[T, R] = r_
end

# Sort keys
KEYS = sort(collect(keys(res)))

# Print the table
for (T, R) in KEYS
    r = res[T, R]
    @printf "DER-%d & %5d & %3d & %5.1f & %5.1f & %3d & %6d & %.1f \$\\%%\$ \\\\\n" T R r["cgiter"] r["tmp"] r["tsp"] r["m0"]  r["N"] 100*r["density"]
end