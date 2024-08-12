using CellContraction

using CSV
using DataFrames
using Glob


function main()
    # Loading simulation results
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # assuming you have run all the simulations on the HPC
    mechanisms = setdiff(readdir(path_to_outsim), ["biomarkers.csv", "Beta_increase"])

    df0 = CSV.read(joinpath(path_to_outsim, "biomarkers.csv"), DataFrame; header=1, types=Float64)

    # Collecting biomarker values for the full population of models for all concentrations for each simulated compound
    df = Dict()
    for mech in mechanisms
        df[mech] = Dict()
        path_to_mech = joinpath(path_to_outsim, mech)
        files = glob("*.csv", path_to_mech)
        for f in files
            idx_param_value = split(split(f, "biomarkers_m")[2], ".csv")[1]
            df[mech][idx_param_value] = CSV.read(f, DataFrame; header=1, types=Float64)
        end
    end
    
    # Run perturbation-response analysis for all simulated mechanisms of action
    at_least_models_n = 50  # we want at least 50 converging models for all simulated perturbations of a given mechanism 
    bio = "T_peak"
    for (i, mech) in enumerate(mechanisms)
        println("Analysing mechanism ($(i)/$(length(mechanisms))): $(mech)")
        runPerturbResponseAnalysis(df, df0, mech, bio, at_least_models_n, path_to_outsim)
    end

end

main()
