using ArgParse

include("experiments.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--T", "-T"
            help = "Length of the time horizon"
            arg_type = Int
            default = 24
        "--R", "-R"
            help = "Number of resources"
            arg_type = Int
            default = 1024
        "--seed", "-s"
            help = "Random seed"
            arg_type = Int
            default = 42
    end

    return parse_args(s, as_symbols=true)
end

function main()

    args = parse_commandline()

    rmp_solvers = [:TLP, :MSK, :GRB, :CPX, :SGRB, :SCPX]

    # Compilation round
    run_experiment(;
        T0=408+7, T=6, R=16, ξ=0.33, seed=0,
        solvers=rmp_solvers,
        cg_verbose=false,
        export_res=false, res_folder="res/"
    )

    # Real business
    run_experiment(;
        T0=408+7,
        T=args[:T],
        R=args[:R],
        ξ=0.33,
        seed=args[:seed],
        solvers=rmp_solvers,
        cg_iter_max=200,
        time_limit=7200.0,
        cg_verbose=true,
        export_res=true,
        res_folder="res/"
    )

    println("Done.")
    flush(Base.stdout)

    return nothing
end

main()