# Generate list of jobs
println("Creating jobs files")
S = [42, 84, 21, 93, 5]                     # Random seeds
T = [24, 48, 96]                            # Length of time horizon
R = [1024, 2048, 4096, 8192, 16384, 32768]  # Number of resources

open("jobs.txt", "w") do jobs_file
    for s in S, t in T, r in R
        write(jobs_file, "julia --project=. scripts/benchmark.jl -T $t -R $r -s $s --export > scripts/out/output_$(t)_$(r)_$(s).out 2>&1\n")
    end
end