module Bayesian
export runDoseResponseAnalysis

using ForceImport
@force using ..Drug
@force using ..Utils

using CairoMakie
using CSV
using DataFrames
using LinearAlgebra
using OrderedCollections
using Turing
using Statistics
using StatsPlots


function runDoseResponseAnalysis(df_new, df_props, cmpd, path_to_outsim)
    df_cmpd = df_new[cmpd]
    EFTPCmax = getCmpdPropLess(df_props, cmpd)
    x = getDoses(EFTPCmax, collect(keys(df_cmpd)))
    Y = collect(values(df_cmpd))

    y_median = dropdims(median(reduce(hcat, Y), dims=1), dims=1)
    idx_all = findall(y_median .< 1e-2)
    idx_end = length(x)
    if ~isnothing(idx_all) & (length(idx_all) > 3)
        idx_end = idx_all[3]
    end
    x = x[1:idx_end]

    Y = [[Y[i][j] for i in 1:idx_end] for j in eachindex(Y[1])]

    path_to_bayesian = joinpath(path_to_outsim, cmpd, "Bayesian")
    mkpath(path_to_bayesian)

    xr, xl, sdy = prior_params(x, Y)

    @. sigmoid(x, p) = p[3] + (1 - p[3]) * (1 + 10 ^ (p[2]*(x - p[1]))) ^ (-1)

    @model function bayesianModel(x, y)
        pIC50 ~ Normal(0.5*(xl+xr), 5)
        Hill ~ truncated(Normal(1, 5); lower=0)
        B ~ Normal(0.5, 0.5)
        σ² ~ truncated(Normal(0, sdy); lower=0)
        mu = sigmoid(x, [pIC50, Hill, B])
        return y ~ MvNormal(mu, σ² * I)
    end

    model_train = bayesianModel(x, Y)
    chain = sample(model_train, NUTS(0.95), MCMCSerial(), 1_000, 4)
    save_stats(chain, path_to_bayesian)
    chainplt = StatsPlots.plot(chain)
    savefig(chainplt, joinpath(path_to_bayesian, "params_distr_chains.png"))

    x_test, y_mean, y_lower, y_upper = predict_test(x, chain, bayesianModel)

    fig = with_theme(() -> plot_fit(df_cmpd, EFTPCmax, x_test, y_mean, y_lower, y_upper, length(x)), theme_minimal())
    save(joinpath(path_to_bayesian, "sample_posterior_predictive_mean.pdf"), fig, pt_per_unit=1)
end


function plot_fit(df_cmpd, EFTPCmax, x_test, y_mean, y_lower, y_upper, n_doses)
    size_inches = (2 * 2/3 * 8.27, 2 * 1/3 * 11.69)
    size_pt = 72 .* size_inches
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(keys(df_cmpd)))]
    f = CairoMakie.Figure(resolution=size_pt, fontsize=12)
    axis = CairoMakie.Axis(f[1, 1], xlabel=rich("log", subscript("10"), "(Dose [M])"), ylabel="Peak Tension\n(fraction of control)")
    CairoMakie.lines!(axis, x_test, y_mean, color=Makie.wong_colors()[1], label="Posterior samples\npredicted mean")
    CairoMakie.band!(axis, x_test, y_lower, y_upper, color=(Makie.wong_colors()[1], 0.3), label="95% CI")
    axislegend("Bayesian Model", position=:rt)
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


function save_stats(chain, path_to_bayesian)
    pIC50_samples = reshape(Array(chain[:pIC50]), :)
    IC50_samples = map(x -> 1e6*10^x, pIC50_samples)
    IC50_mean = mean(IC50_samples)
    IC50_sd = std(IC50_samples)
    df_summary = summarize(chain)
    data = collect(values(df_summary[:, 2:end]))
    data = [data; [IC50_mean, IC50_sd, NaN, NaN, NaN, NaN, NaN]']
    headers = map(string, names(df_summary))[2:end]
    rows = map(string, df_summary[:, :parameters])
    pop!(rows)
    push!(rows, "sigma2")
    push!(rows, "IC50")
    df1 = DataFrame("" => rows)
    df2 = DataFrame(data, headers)
    df_params_stats = hcat(df1, df2)
    CSV.write(joinpath(path_to_bayesian, "params_stats.csv"), df_params_stats)
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


function predict_test(x, chain, bayesianModel)
    N = Int64(1e2 * length(x))
    x_test = range(x[1], x[end], length=N)
    model_test = bayesianModel(x_test, missing)
    predictions = predict(model_test, chain)
    Y_pred = Array(predictions)
    y_mean = vec(mapslices(x -> mean(x), Y_pred; dims=1))
    y_lower = vec(mapslices(x -> quantile(x, 0.025), Y_pred; dims=1))
    y_upper = vec(mapslices(x -> quantile(x, 0.975), Y_pred; dims=1))
    [x_test, y_mean, y_lower, y_upper]
end


end  # module
