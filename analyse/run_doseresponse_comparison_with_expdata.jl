using CairoMakie
using CSV
using DataFrames
using Printf


function plotFigure(df, df_, method)
    size_inches = (2 * 8.27, 2 * 0.225 * 11.69)
    size_pt = 72 .* size_inches
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=nrow(df))]
    f = Figure(size=size_pt, fontsize=14)
    axis = Axis(f[1, 1], xticks=1:nrow(df), ylabel="Delta (Order of Magnitude)")
    ooms = [-1, 0, 1, 2]
    for oom in ooms
        hlines!(axis, [oom], color=:gray, linewidth=0.8)
    end
    s = []
    for (i, val) in enumerate(df_[!, :oom])
        si = scatter!(axis, i, val, color=colors[i], markersize=20)
        push!(s, si)
    end
    legends = ["$(lpad(i, 2, ' ')). " for i=1:length(df_[!,:Compound])] .* df_[!,:Compound]
    idx = [i for i in vcat([[i, i+11] for i in 1:11]...)]
    Legend(f[1, 2], s[idx], legends[idx], "Compounds", nbanks=2)
    title = rich("IC50", subscript("pred"), " (active tension peak) vs. IC50", subscript("exp"), " (sarcomere shortening)")
    Label(f[0, :], title, fontsize=20)
    f
end

function main(method::String)
    path_to_expdata = joinpath(pwd(), "analyse", "data")
    df_exp = CSV.read(joinpath(path_to_expdata, "Nguyen2017.csv"), DataFrame; header=1)

    path_to_fitted_params = joinpath(pwd(), "analyse", "output")
    df_pred = CSV.read(joinpath(path_to_fitted_params, "dr_fitted_parameters_$(method).csv"), DataFrame; header=1)

    col = filter(col -> startswith(col, "IC50"), names(df_pred))[1]

    df_ = df_pred[!, [Symbol(col)]]
    if method == "Bayesian"
        df_mean = select(df_pred, Symbol(col) => (x -> parse.(Float64, map(x -> x[1], split.(x, " (")))) => Symbol("IC50 mean"))
        df_std = select(df_pred, Symbol(col) => (x -> parse.(Float64, map(x -> x[2][1:end-1], split.(x, " (")))) => Symbol("IC50 std"))
        df_ = hcat(df_mean, df_std)
    end
    df = hcat(df_exp, df_)

    # Format summary table to export
    df_to_csv = copy(df)
    if method == "Bayesian"
        df_to_csv[!, :IC50pred] = string.(df_to_csv[!, Symbol("IC50 mean")]) .* "(" .* string.(df_to_csv[!, Symbol("IC50 std")]) .* ")"
        select!(df_to_csv, Not(Symbol("IC50 mean")))
        select!(df_to_csv, Not(Symbol("IC50 std")))
    else
        rename!(df_to_csv, Symbol("IC50 mean") => :IC50pred)
    end
    rename!(
        df_to_csv,
        :IC50 => Symbol("Observed IC50 (μM) (sarcomere shortening)"),
        :IC50pred => Symbol("Predicted IC50 (μM) (peak active tension)")
    )
    CSV.write(joinpath(path_to_fitted_params, "dr_comparison_summary_$(method).csv"), df_to_csv)

    # Compare with experiments only a subset of all the 28 tool compounds
    discard = ["Cisapride", "Dofetilide", "Moxifloxacin", "Sotalol"]
    df = filter(row -> !startswith(row.IC50, ">"),  df)  # filtering out compounds which did not show a change in sarcomere shortening above at least 25% from baseline in the experiments (if there is an effect, the authors suggest it must be above max tested conc.)
    df = filter(row -> row.Compound ∉ discard, df)  # filtering out compounds which, from simulations, resulted to have no effect on contractility (did not reach at least 25% inhibition from baseline - mimicking strategy used in the paper)
    df.IC50 = parse.(Float64, df.IC50)  # now all elements are numbers so can convert String to Float64

    # Compute mismatch in order of magnitude (base 10) of IC50_pred with respect to order of magnitude of IC50_exp
    mismatch(x1, x2) = parse.(Int64, [(@sprintf "%.2E" x2i)[end-2:end] for x2i in x2]) - parse.(Int64, [(@sprintf "%.2E" x1i)[end-2:end] for x1i in x1])
    df_ = select(df, :Compound, [:IC50, Symbol("IC50 mean")] => ((x1, x2) -> mismatch(x1, x2)) => :oom)
    sort!(df_, [:oom])

    fig = with_theme(() -> plotFigure(df, df_, method), theme_minimal())
    save(joinpath(path_to_fitted_params, "dr_simul_vs_expdata_$(method).png"), fig, px_per_unit=5)
end

method = "Median"  # use either "Median" or "Bayesian"
main(method)
