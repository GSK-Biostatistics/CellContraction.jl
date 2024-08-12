using CSV
using DataFrames
using Distributed


n_procs = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(n_procs)

@everywhere using CellContraction

df_in = CSV.read(joinpath(pwd(), "data", "population.csv"), DataFrame; header=1, types=Float64)
df = CSV.read(joinpath(pwd(), "data", "mechanisms.csv"), DataFrame; header=1, types=String)

mech_idx = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
mech = df.Mechanism[mech_idx]

d_factors = getMechScaling(df, mech_idx)
d_combinations = genMechCombs(d_factors)

mkpath(joinpath(pwd(), "output", mech))

combinations = mechAlter(df_in, d_combinations)

for comb in keys(combinations)
    df_in_ = combinations[comb]
    
    inputs = [values(row) for row in eachrow(df_in_[!, 2:end])]
    n_beats = 1000
    X = pmap(p -> runSimulCompBio(p, n_beats), inputs)

    df_out = initOutData()
    for x in X
        push!(df_out, x)
    end

    CSV.write(joinpath(pwd(), "output", mech, "biomarkers_$(comb).csv"), df_out)
end

for w in workers()
    rmprocs(w)
end
