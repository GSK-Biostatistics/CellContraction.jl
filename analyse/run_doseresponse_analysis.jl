using CellContraction

using CSV
using DataFrames
using Glob
using OrderedCollections
using Random
using Statistics


function assembleFittedParameters(method::String, path_to_outsim::String, compounds::Vector{String})
    # Assemble the fitted parameter values for all the compounds within a single csv file
    parameters = ["IC50", "Hill", "B"]

    suffix = " mean"
    if method == "Bayesian"
        suffix = suffix * " (std)"
    end
    col_names = cat("Compound", parameters .* suffix, dims=1)

    df_fitted = DataFrame(col_names .=> [T[] for T in [String]])

    for cmpd in compounds
        path_to_doseresponse = joinpath(path_to_outsim, cmpd, method)
        df = CSV.read(joinpath(path_to_doseresponse, "params_stats.csv"), DataFrame; header=1)

        data_i = (cmpd,)
        for param in parameters
            df_p = filter(row -> row.Column1 == param, df)
            mean = df_p.mean[1]
            valstring = string(round(mean, digits=4))
            if method == "Bayesian"        
                std = df_p.std[1]
                valstring = valstring * " ($(round(std, digits=4)))"
            end
            data_i = (data_i..., valstring)
        end
        push!(df_fitted, data_i)
    end

    path_to_out = joinpath(pwd(), "analyse", "output")
    mkpath(path_to_out)
    CSV.write(joinpath(path_to_out, "dr_fitted_parameters_$(method).csv"), df_fitted)
end

function main(method::String)
    Random.seed!(8)  # enforcing reproducibility

    # Loading simulation results
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # assuming you have run all the simulations on the HPC
    compounds = setdiff(readdir(path_to_outsim), ["biomarkers.csv"])

    df0 = CSV.read(joinpath(path_to_outsim, "biomarkers.csv"), DataFrame; header=1, types=Float64)

    path_to_data = joinpath(pwd(), "run", "hpc", "data")
    df_props = CSV.read(joinpath(path_to_data, "tool_compounds.csv"), DataFrame; header=1)

    # Collecting biomarker values for the full population of models for all concentrations for each simulated compound
    df = Dict()
    for cmpd in compounds
        df[cmpd] = Dict()
        path_to_cmpd = joinpath(path_to_outsim, cmpd)
        files = glob("*.csv", path_to_cmpd)
        for f in files
            conc = parse(Float64, split(split(f, "_")[end], ".csv")[1])
            df[cmpd][conc] = CSV.read(f, DataFrame; header=1, types=Float64)
        end
    end
    concentrations = sort(collect(keys(df[compounds[1]])))

    # Filtering out models that did not lead to a fully repolarised AP, CaT, and AT.
    # If less than 50 models are left in total at a given concentration, remove that concentration.
    # If for a concentration we are left with a subset of models, select the very same (and only those) models for the other concentrations too
    df_new = OrderedDict()
    conc_max = concentrations[1]
    idx_shared = collect(1:size(df0)[1])
    for cmpd in compounds
        df_new[cmpd] = OrderedDict()
        for conc in reverse(concentrations)
            idx = getConverging(df[cmpd][conc])
            if length(idx) > 50  # we want at least 50 'converging' models from the whole population
                conc_max = conc
                idx_shared = idx
                break
            end
        end
        for conc in concentrations
            if conc <= conc_max
                push!(df_new[cmpd], conc => df[cmpd][conc][idx_shared, :].T_peak ./ median(df0.T_peak))  # normalise by unperturbed population values
            end
        end
    end

    # Run dose-response analysis for all simulated compounds
    method_module = eval(Symbol(method))
    for (i, cmpd) in enumerate(compounds)
        println("Analysing compound ($(i)/$(length(compounds))): $(cmpd)")
        method_module.runDoseResponseAnalysis(df_new, df_props, cmpd, path_to_outsim)
    end

    # Join all results
    assembleFittedParameters(method, path_to_outsim, compounds)
end

method = "Median"  # use either "Median" or "Bayesian"
main(method)
