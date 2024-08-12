using CellContraction

using CairoMakie
using CSV
using DataFrames
using Glob
using OrderedCollections
using Statistics


new_label = Dict("Ca50"=>"Ca50","Cao"=>"Cao","GCaL"=>"GCaL","Gncx"=>"GNCX","GNaK"=>"GNaK","Myo"=>"Kuw","GJrel"=>"Jrel","GJup"=>"Jup")

function plotFigure(f, i, j, x, Y, param, idx, mech)
    delta = minimum([x[i+1]-x[i] for i=1:length(x)-1])    
    xlab = (i == 4) ? "Scaling Factor" : ""
    ylab = (j == 1) ? rich("AT", subscript("peak"), "\n(fraction of contr.)") : ""
    axis = Axis(f[i, j], xlabel=xlab, ylabel=ylab, limits=(x[1]-delta, x[end]+delta, nothing, nothing), title="$(mech) ($(new_label[param[1]]))", xticks=x)
    for i in idx
        boxplot!(axis, fill(x[i], length(Y[i])), Y[i], width=0.5*delta, color=Makie.wong_colors()[1], show_outliers=false)
    end
    f
end

function main()
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # note: in this folder you must have only mechanism perturbation simulations, NOT drug block simulations
    mechanisms = setdiff(readdir(path_to_outsim), ["biomarkers.csv", "Beta_increase"])

    df0 = CSV.read(joinpath(path_to_outsim, "biomarkers.csv"), DataFrame; header=1, types=Float64)

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

    at_least_models_n = 50
    bio = "T_peak"

    size_inches = (2 * 8.27, 2 * 0.65 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=18)
    
    idx_preferred_order = [2, 1, 3, 6, 4, 7, 5, 8]
    n_cols = 2
    for (k, mech) in enumerate(mechanisms[idx_preferred_order])
        i = (k - 1) ÷ n_cols + 1
        j = mod(k - 1, n_cols) + 1
       
        idx_list, sc_list = get_baseline(df, mech)
        idx_shared, idx_best = get_shared(idx_list, at_least_models_n)
        df_bio = build_mech(idx_shared, sc_list, idx_best, df0, df, mech, bio)

        labels = collect(keys(df_bio))
        M = reduce(hcat, map(x->split(x, "_"), labels))
        dims = div(size(M)[1], 2)
        param_idx = [2*i for i=1:dims]
        param = M[param_idx, 1]
        value_idx = [2*i+1 for i=1:dims]

        x = parse.(Float64, M[value_idx, :])
        Y = collect(values(df_bio))
        y_median = vcat(median.(values(df_bio), dims=1)...)
        x = vec(x)
        idx = deleteat!(collect(1:length(y_median)), isnan.(y_median))
        f = plotFigure(f, i, j, x, Y, param, idx, mech)
    end
    
    path_to_output = joinpath(pwd(), "analyse", "output")
    save(joinpath(path_to_output, "$(bio)_all_mechanisms.png"), f, px_per_unit=5)
end

main()
