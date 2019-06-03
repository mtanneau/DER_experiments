using CSV
using DataFrames
using Statistics

CUR_DIR = "$(@__DIR__)"
LOG_DIR = "$(@__DIR__)/log/"
OUT_DIR = "$(@__DIR__)/out/"
RES_DIR = "$(@__DIR__)/res/"

"""
    gmean1(u)

Compute the geometric mean with a shift of 1.
"""
gmean1(u) = exp(mean(log.(u .+ oneunit(eltype(u))))) - oneunit(eltype(u))

"""

"""
function export_results_table(
    solvers, metrics, df_,
    res_file_name::String
)
    u = []

    for s in solvers
        df_tmp = df_[df_.solver .== s, metrics]
        rename!(sb -> (sb == :T || sb== :R) ? sb : Symbol(String(sb) * "_" * s), df_tmp)
        push!(u, df_tmp)
    end
    
    df_merged = u[1]

    for d in u[2:end]
        df_merged = join(df_merged, d, on=[:T, :R])
    end

    sort!(df_merged, (:T, :R))
    
    return df_merged 
end


# Parse log files to get solver-specific logs
out_files = readdir(OUT_DIR)
@info "Parsing $(length(out_files)) output files"
for f in out_files

    s = open(OUT_DIR * f) do file
        read(file, String)
    end

    s = split(s, "~~FLAG_EXPERIMENT~~")[2]
    # println(s)
    s_ = split(s, "#"^72)
    
    println(f)
    for log in s_[2:end]
        g_ = split(log, '`')[2][1:end-4]*".out"
        # println("\t", g_)
        
        open(LOG_DIR * g_, "w") do logfile
            write(logfile, log)
        end
        
    end
    
    # println()
end


# Read results files and concatenate into single dataframe
@info "Importing results"
res_files = [f for f in readdir(RES_DIR) if f[end-3:end] == ".csv"]
df = vcat(CSV.read.(RES_DIR .* res_files)...)

# Find list of all jobs that failed
files = readdir(RES_DIR);
res_files = [f for f in files if f[end-3:end] == ".csv"]
out_files = [f for f in files if f[end-3:end] == ".out"]

no_res = []

for f in out_files
    if ! (f[1:end-4] * ".csv" in res_files)
        push!(no_res, f)
    end
end

#======================
    Process results
======================#
@info "Computing aggregate metrics"
# Compute average time per inner iteration
df.t_inner_iter = df.time_mp_total ./ (df.num_iter_bar .+ df.num_iter_spx);

# Change status from categorical to numerical
df.status_ind = [s == "Optimal" ? 100.0 : 0.0 for s in df.status]
deletecols!(df, :status);

# `dual_bound` may be infinite, so we remove it
deletecols!(df, :dual_bound)

# Aggregate results
df_ = aggregate(df, [:T, :R, :solver], gmean1);

@info "Exporting results"
df_all = export_results_table(
    ["TLP", "MSK", "GRB", "CPX", "SCPX", "SGRB"],
    [
        :T, :R, :time_sp_total_gmean1, :time_mp_total_gmean1, :time_cg_total_gmean1,
        :t_inner_iter_gmean1,
        :n_cg_iter_gmean1, :status_ind_gmean1,
        :num_iter_bar_gmean1, :num_iter_spx_gmean1,
        :num_cols_tot_gmean1, :nsp_priced_gmean1
    ],
    df_,
    ""
);
CSV.write("$(CUR_DIR)/results.csv", df_all);

@info "Done"