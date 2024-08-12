module Median
export runDoseResponseAnalysis

using ForceImport
@force using ..Drug
@force using ..Utils

using CairoMakie
using CSV
using DataFrames
using LinearAlgebra
using LsqFit
using OrderedCollections
using Statistics


function runDoseResponseAnalysis(df_new, df_props, cmpd, path_to_outsim)
    df_cmpd = df_new[cmpd]
    EFTPCmax = getCmpdPropLess(df_props, cmpd)
    x = getDoses(EFTPCmax, collect(keys(df_cmpd)))
    Y = collect(values(df_cmpd))
    xr, xl, _ = prior_params(x, [[Y[i][j] for i in eachindex(Y)] for j in eachindex(Y[1])])
    y_median = dropdims(median(reduce(hcat, Y), dims=1), dims=1)

    idx_all = findall(y_median .< 1e-2)
    idx_end = length(x)
    if ~isnothing(idx_all) & (length(idx_all) > 3)
        idx_end = idx_all[3]
    end
    x = x[1:idx_end]
    y_median = y_median[1:idx_end]

    path_to_median = joinpath(path_to_outsim, cmpd, "Median")
    mkpath(path_to_median)

    @. sigmoid(x, p) = p[3] + (1 - p[3]) * (1 + 10 ^ (p[2]*(x - p[1]))) ^ (-1)
    p0 = [0.5*(xl + xr), 1, 0.1]
    fit = curve_fit(sigmoid, x, y_median, p0)
    best_fit = coef(fit)
    save_stats(best_fit, path_to_median)

    N = Int64(1e2 * length(x))
    x_test = range(x[1], x[end], length=N)
    y_median = sigmoid(x_test, best_fit[1:end-1])
    fig = with_theme(() -> plot_fit(df_cmpd, EFTPCmax, x_test, y_median, length(x)), theme_minimal())
    save(joinpath(path_to_median, "fit_to_median_values.pdf"), fig, pt_per_unit=1)
end


function plot_fit(df_cmpd, EFTPCmax, x_test, y_median, n_doses)
    size_inches = (2 * 2/3 * 8.27, 2 * 1/3 * 11.69)
    size_pt = 72 .* size_inches
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(keys(df_cmpd)))]
    f = CairoMakie.Figure(resolution=size_pt, fontsize=12)
    axis = CairoMakie.Axis(f[1, 1], xlabel=rich("log", subscript("10"), "(Dose [M])"), ylabel="Peak Tension\n(fraction of control)")
    CairoMakie.lines!(axis, x_test, y_median, color=Makie.wong_colors()[1], linewidth=3, label="Fit to median values")
    axislegend("Deterministic Model", position=:rt)
    boxes = []
    for (i, key) in enumerate(keys(df_cmpd))
        if i > n_doses
            break
        end
        b = CairoMakie.boxplot!(axis, fill(getDoses(EFTPCmax, key), length(df_cmpd[key])), df_cmpd[key], width=0.4, color=colors[i])
        push!(boxes, b)
    end
    labels = map(x->string(x)*"x", collect(keys(df_cmpd)))
    CairoMakie.Legend(f[1, 2], boxes, labels[1:n_doses], rich("Doses as multipliers\nof EFTPC", subscript("max")))
    f
end


function save_stats(best_fit, path_to_median)
    push!(best_fit, 1e6 * 10 ^ best_fit[1])
    df_params_stats = DataFrame(col => val for (col, val) in zip(["pIC50", "Hill", "B", "IC50"], best_fit))
    rename!(df_params_stats, :first => "", :second => "mean")
    CSV.write(joinpath(path_to_median, "params_stats.csv"), df_params_stats)
end


function prior_params(x, Y)
    massimo = maximum(Y, dims=1)[1]
    minimo = minimum(Y, dims=1)[1]
    ir = findfirst(massimo .< 0.5)
    if isnothing(ir)
        ir = length(x)
    end
    il = findfirst(minimo .< 0.5)
    if isnothing(il)
        il = 1
    end
    sdy = 5 * maximum(std(Y))
    [x[ir], x[il], sdy]
end


end  # module
