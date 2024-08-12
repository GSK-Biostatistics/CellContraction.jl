module Mech
export runPerturbResponseAnalysis, get_baseline, get_shared, build_mech

using ForceImport
@force using ..Utils

using CairoMakie
using Combinatorics
using CSV
using DataFrames
using OrderedCollections
using Polynomials
using Statistics


function runPerturbResponseAnalysis(df, df0, mech, bio, at_least_models_n, path_to_outsim)
    path_to_mech = joinpath(path_to_outsim, mech, "LR")
    mkpath(path_to_mech)

    idx_list, sc_list = get_baseline(df, mech)
    idx_shared, idx_best = get_shared(idx_list, at_least_models_n)
    df_mech = build_mech(idx_shared, sc_list, idx_best, df0, df, mech, bio)

    labels = collect(keys(df_mech))
    M = reduce(hcat, map(x->split(x, "_"), labels))

    dims = div(size(M)[1], 2)
    param_idx = [2*i for i=1:dims]
    param = M[param_idx, 1]

    value_idx = [2*i+1 for i=1:dims]
    x = parse.(Float64, M[value_idx, :])
    Y = collect(values(df_mech))
    y_median = vcat(median.(values(df_mech), dims=1)...)
    y_std = vcat(std.(values(df_mech), dims=1)...)
    n_models = length(Y[1])

    if size(param)[1] == 1
        x = vec(x)
        idx = deleteat!(collect(1:length(y_median)), isnan.(y_median))
        deg = 2
        poly = Polynomials.fit(x[idx], y_median[idx], deg)
        save_stats(coeffs(poly), x, param, path_to_mech)
    end

    if size(param)[1] == 1
        fig = with_theme(() -> plot_figure_1D(x, Y, param, bio, poly, idx, n_models), theme_minimal())
    else
        fig = with_theme(() -> plot_figure_2D(x, y_median, y_std, param, bio, n_models), theme_minimal())
    end
    save(joinpath(path_to_mech, "fit_to_median_values.pdf"), fig, pt_per_unit=1)
end


function get_baseline(df, mech)
    idx_list = []
    sc_list = []
    for sc in sort(names(df[mech]), by=x->parse(Int64, x[1:findfirst("_", x)[1]-1]))
        idx = getConverging(df[mech][sc])
        push!(idx_list, idx)
        push!(sc_list, sc)
    end
    return idx_list, sc_list
end


function get_shared(idx_list, at_least_models_n)
    n_scales = length(idx_list)
    idx_shared = sort(intersect(idx_list...))
    idx_best = collect(1:n_scales)
    if length(idx_shared) >= at_least_models_n
        return idx_shared, idx_best
    end
    
    for i=1:div(n_scales, 2)
        combs_i = collect(combinations(1:n_scales, i))
        idx_common_combs = []
        for comb in combs_i
            idx_left = setdiff(idx_best, idx_best[comb])
            idx_common = sort(intersect(idx_list[idx_left]...))
            push!(idx_common_combs, idx_common)
        end
        i_sort_common = sortperm(idx_common_combs, by=length)
        best_comb = combs_i[i_sort_common[end]]
        idx_shared = idx_common_combs[i_sort_common[end]]
        if length(idx_shared) >= at_least_models_n
            return idx_shared, setdiff(idx_best, best_comb)
        end
    end
end


function build_mech(idx_shared, sc_list, idx_best, df0, df, mech, bio)
    sc_list_new = [sc_list[i] for i in idx_best]
    df_mech = OrderedDict()
    for sc in sc_list
        if sc in sc_list_new
            push!(df_mech, sc => df[mech][sc][idx_shared, :][!, Symbol(bio)] ./ median(df0[idx_shared, :][!, Symbol(bio)]))  # normalise by unperturbed population values
        else
            push!(df_mech, sc => fill(NaN, length(idx_shared)))
        end
    end
    return df_mech
end


function save_stats(coefficients, x, param, path_to_mech)
    n = length(coefficients)
    headers = ["x^$(i)" for i in 0:n-1]
    df_params_stats = DataFrame(headers .=> [T[] for T in [Float64]])
    push!(df_params_stats, coefficients)
    insertcols!(df_params_stats, :domain => [x])
    insertcols!(df_params_stats, :param => [String(param[1])])
    CSV.write(joinpath(path_to_mech, "params_stats.csv"), df_params_stats)
end


function plot_figure_1D(x, Y, param, bio, poly, idx, n_models)
    if bio == "T_peak"
        bio = "Peak Tension"
    end
    size_inches = (2 * 2/3* 8.27, 2 * 1/3 * 11.69)
    size_pt = 72 .* size_inches
  
    delta = minimum([x[i+1]-x[i] for i=1:length(x)-1])
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(Y))]
    f = CairoMakie.Figure(resolution=size_pt, fontsize=12)
    axis = CairoMakie.Axis(f[1, 1], limits=(x[1]-delta, x[end]+delta, nothing, nothing), xlabel="$(param[1]) scaling factor", ylabel="$(bio)\n(fraction of control)")
    boxes = []
    for i in idx
        b = CairoMakie.boxplot!(axis, fill(x[i], length(Y[i])), Y[i], width=0.5*delta, color=colors[i], show_outliers=false)
        push!(boxes, b)
    end

    x_test = range(x[1], x[end], length=Int64(1e2 * length(x)))
    CairoMakie.lines!(axis, x_test, poly.(x_test), color=Makie.wong_colors()[1], linewidth=3, label="LR (deg=2) fit to sub-pop. medians")

    axislegend("#models = $(n_models)", position=:ct)
    labels = map(x->param[1]*"="*string(x), x)
    CairoMakie.Legend(f[1, 2], boxes, labels[idx], "Factors")
    f
end


function plot_figure_2D(x, y_median, y_std, param, bio, n_models)
    if bio == "T_peak"
        bio = "Peak Tension"
    end
    size_inches = (2 * 8.27, 2/5 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(resolution=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], title="#models = $(n_models)", xlabel="$(param[1]) scaling factor", ylabel="$(param[2]) scaling factor")
    dx = vec(x[1, :])
    dy = vec(x[2, :])
    hm = heatmap!(ax1, dx, dy, y_median, interpolate=false)
    Colorbar(f[1, 2], hm, label="$(bio) - Median\n(fraction of control)")
    ax2 = Axis(f[1, 3], xlabel="$(param[1]) scaling factor", ylabel="$(param[2]) scaling factor")
    hm = heatmap!(ax2, dx, dy, y_std, interpolate=false)
    Colorbar(f[1, 4], hm, label="$(bio) - SD\n(fraction of control)")
    f
end


end  # module
