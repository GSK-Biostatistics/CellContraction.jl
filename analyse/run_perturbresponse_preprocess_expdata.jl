using CairoMakie
using CSV
using DataFrames
using Glob
using LsqFit


function plotFigure(X, Y, compounds)
    size_inches = (1.2 * 8.27, 1.2 * 0.36 * 11.69)
    size_pt = 72 .* size_inches
    colors = [cgrad(:seaborn_colorblind, [0, 1])[z] for z ∈ range(0, 1, length=length(compounds))]
    f = Figure(size=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], xlabel=rich("log", subscript("10"), "(Concentration [M])"), ylabel="Sarcomere shortening\n(% control)")
    l = []
    for (i, (x, y)) in enumerate(zip(X, Y))
        li = lines!(ax1, x, y, color=colors[i], linewidth=3)
        push!(l, li)
    end
    Legend(f[1, 2], l, compounds, "Compounds")
    f
end

function main()
    path_to_expdata = joinpath(pwd(), "analyse", "data", "digitized_figures")
    path_to_ec50 = joinpath(pwd(), "analyse", "data")

    df_EC50 = CSV.read(joinpath(path_to_ec50, "AbiGerges2020.csv"), DataFrame; header=1)
    EC50 = Dict()
    for (key, val) in zip(df_EC50[!, "Compound"], df_EC50[!, "EC50"])
        EC50[key] = val
    end

    files = glob("*.csv", path_to_expdata)
    data = Dict()
    for f in files
        cmpd = split(basename(f), ".csv")[1]
        data[cmpd] = CSV.read(f, DataFrame; header=["x", "y"], types=Float64)
    end

    X = []
    Y = []

    col_names = ["Compound", "Hill", "B", "T", "EC50", "EC50_exp"]
    col_types = [String[], Float64[], Float64[], Float64[], Float64[], Float64[]]
    df_fitted = DataFrame(col_names .=> col_types)

    for cmpd in keys(data)
        ec50 = EC50[cmpd]
        x = data[cmpd][!, "x"]
        idx = sortperm(x)
        x = x[idx]
        y = data[cmpd][!, "y"]
        y = y[idx]
        x = log10.(1e-6 .* x)
        pEC50 = log10(1e-6 * ec50)
    
        @. sigmoid(x, p) = p[2] + (p[3] - p[2]) * (1 + 10 ^ (p[1]*(-x + p[4]))) ^ (-1)
        p0 = [1, 100, y[end], pEC50]
        fit = curve_fit(sigmoid, x, y, p0)
        best_fit = coef(fit)

        N = Int64(1e2 * length(x))
        x_test = range(x[1], x[end], length=N)
        yy = sigmoid(x_test, best_fit)
        
        push!(X, x_test)
        push!(Y, yy)

        best_fit[end] = 1e6 * 10 ^ (best_fit[end])
        push!(df_fitted, (cmpd, best_fit..., ec50))
    end

    cols_to_transform = setdiff(names(df_fitted), ["Compound", "EC50_exp"])
    for col in cols_to_transform
        transform!(df_fitted, col => ByRow(x -> "$(round(x, digits=4))") => col)
    end
    path_to_out = joinpath(pwd(), "analyse", "output")
    CSV.write(joinpath(path_to_out, "pr_fitted_parameters.csv"), df_fitted)

    fig = with_theme(() -> plotFigure(X, Y, collect(keys(data))), theme_minimal())
    save(joinpath(path_to_out, "pr_hill_fit_to_input_exp_data.png"), fig, px_per_unit=5)
end

main()
