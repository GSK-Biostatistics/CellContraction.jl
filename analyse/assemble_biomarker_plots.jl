using CellContraction

using CairoMakie
using CSV
using DataFrames
using Glob
using OrderedCollections
using Statistics


function plotFigure(f, i, j, x, Y, Y0, bio, idx, mech)
    delta = minimum([x[i+1]-x[i] for i=1:length(x)-1])
    if i == 6
        axis = CairoMakie.Axis(f[i, j], xlabel=mech, limits=(x[idx[1]]-delta, x[idx[end]]+delta, nothing, nothing), xticks=([0, 0.5, 1], ["0%", "50%", "100%"]), ylabel="$(bio)")
    else
        axis = CairoMakie.Axis(f[i, j], limits=(x[idx[1]]-delta, x[idx[end]]+delta, nothing, nothing), xticks=([0, 0.5, 1], ["0%", "50%", "100%"]), ylabel="$(bio)")
    end
    for i in idx
        CairoMakie.boxplot!(axis, fill(1-x[i], length(Y[i])), Y[i]*median(Y0), width=0.5*delta, color=Makie.wong_colors()[1], show_outliers=false)
    end
    [f, axis]
end

function main()
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # note: in this folder you must have only mechanism perturbation simulations, NOT drug block simulations
    mechanisms = setdiff(readdir(path_to_outsim), ["biomarkers.csv"])

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
    mech = "IKr_block"  # choose which mechanism to plot (a folder with same name needs to be available in the specified output path)
    idx_list, sc_list = get_baseline(df, mech)
    idx_shared, idx_best = get_shared(idx_list, at_least_models_n)

    biomarkers = OrderedDict(  # renaming to more commonly known acronyms
        "RMP"=>"RMP (mV)",
        "V_peak"=>"Vpeak (mV)",
        "APD40"=>"APD40 (ms)",
        "APD90"=>"APD90 (ms)",
        "dVdt_max"=>"dVdtMax (mV/ms)",
        "Ca_diast"=>"CaiD (nM)",
        "Ca_peak"=>"CTpeak (nM)",
        "CTD50"=>"CTD50 (ms)",
        "CTD90"=>"CTD90 (ms)",
        "Tttp"=>"ATttp (ms)",
        "T_peak"=>"ATpeak (kPa)",
        "Trt50"=>"ATrt50 (ms)",
        "Trt90"=>"ATrt90 (ms)",
        "dTdt_max"=>"dATdtMax (kPa/ms)",
        "dTdt_min"=>"dATdtMin (kPa/ms)",
        "EMw90"=>"EMw (ms)",
        "Tri9040"=>"Tri9040 (ms)",
        "qNet"=>"qNet (μC/μF)"
    )

    size_inches = (2 * 8.27, 2 * 0.75 * 11.69)
    size_pt = 72 .* size_inches
    f = CairoMakie.Figure(size=size_pt, fontsize=18)

    labels = nothing
    param = nothing
    idx = nothing
    axes = []

    n_cols = 3
    for (k, bio) in enumerate(keys(biomarkers))
        i = (k - 1) ÷ n_cols + 1
        j = mod(k - 1, n_cols) + 1

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

        Y0 = df0[idx_shared, :][!, bio]
        f, axis = plotFigure(f, i, j, x, Y, Y0, biomarkers[bio], idx, mech)

        if i != Int(length(keys(biomarkers)) / n_cols)
            hidexdecorations!(axis, grid=false, ticks=false)
        end
        push!(axes, axis)
    end
    linkxaxes!(axes...)

    path_to_output = joinpath(pwd(), "analyse", "output")
    save(joinpath(path_to_output, "$(mech).png"), f, px_per_unit=5)
end

main()
