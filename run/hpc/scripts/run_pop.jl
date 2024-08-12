using CSV
using DataFrames
using Distributed


n_procs = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(n_procs)

@everywhere using CellContraction

df_in = CSV.read(joinpath(pwd(), "data", "population.csv"), DataFrame; header=1, types=Float64)

mkpath(joinpath(pwd(), "output"))

inputs = [values(row) for row in eachrow(df_in[!, 2:end])]
n_beats = 1000
X = pmap(p -> runSimulCompBio(p, n_beats), inputs)

df_out = initOutData()
for x in X
    push!(df_out, x)
end

CSV.write(joinpath(pwd(), "output", "biomarkers.csv"), df_out)

for w in workers()
    rmprocs(w)
end
