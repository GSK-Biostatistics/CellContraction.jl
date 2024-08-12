using CSV
using DataFrames
using Distributed


n_procs = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(n_procs)

@everywhere using CellContraction

df_in = CSV.read(joinpath(pwd(), "data", "population.csv"), DataFrame; header=1, types=Float64)
df = CSV.read(joinpath(pwd(), "data", "tool_compounds.csv"), DataFrame; header=1)

cmpd_idx = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
cmpd = df.Compound[cmpd_idx]

mkpath(joinpath(pwd(), "output", cmpd))

concentrations = CSV.read(joinpath(pwd(), "data", "concentrations.csv"), DataFrame; header=1).Multipliers

for conc in concentrations
    df_in_ = copy(df_in)
    poreBlock!(df_in_, df, conc, cmpd_idx)

    inputs = [values(row) for row in eachrow(df_in_[!, 2:end])]
    n_beats = 1000
    X = pmap(p -> runSimulCompBio(p, n_beats), inputs)

    df_out = initOutData()
    for x in X
        push!(df_out, x)
    end

    CSV.write(joinpath(pwd(), "output", cmpd, "biomarkers_$(conc).csv"), df_out)
end

for w in workers()
    rmprocs(w)
end
