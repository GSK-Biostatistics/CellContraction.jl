using CairoMakie
using CSV
using DataFrames
using Glob
using OrderedCollections
using Polynomials
using Statistics


function plotFigure(x, param, y_query, poly)
    size_inches = (2 * 2/3* 8.27, 2 * 1/3 * 11.69)
    size_pt = 72 .* size_inches
    delta = minimum([x[i+1]-x[i] for i=1:length(x)-1])
    f = CairoMakie.Figure(size=size_pt, fontsize=12)
    axis = CairoMakie.Axis(f[1, 1], limits=(x[1]-delta, x[end]+delta, nothing, nothing), xlabel="$(param) scaling factor", ylabel="Peak Tension\n(fraction of control)")

    x_test = range(x[1], x[end], length=Int64(1e2 * length(x)))
    CairoMakie.lines!(axis, x_test, poly.(x_test), color=Makie.wong_colors()[1], linewidth=3)
    r = roots(poly - y_query)
    local x_query = 0
    try
        x_query = minimum(r)
        CairoMakie.vlines!(axis, x_query, color=:red, label="$(param) factor = $(round(x_query, digits=2))")
        CairoMakie.hlines!(axis, y_query, color=:purple, label="Sarcomere shortening\nat the EC50 = $(round(y_query, digits=2))")
        axislegend("Root found", position=:rc)
    catch
        CairoMakie.hlines!(axis, y_query, color=:purple, label="Sarcomere shortening\nat the EC50 = $(round(y_query, digits=2))")
        axislegend("Root NOT found", position=:rc)
    end
    [f, x_query]
end

function main()
    path_to_fitted_params = joinpath(pwd(), "analyse", "output")
    df1 = CSV.read(joinpath(path_to_fitted_params, "pr_fitted_parameters.csv"), DataFrame)
    compounds = setdiff(df1.Compound, ["Isoproterenol", "Epinephrine", "Dobutamine"])  # removing beta agonists given that this analysis is only for 1-parameter mechanism perturbations

    path_to_compounds = joinpath(pwd(), "analyse", "data")
    df2 = CSV.read(joinpath(path_to_compounds, "AbiGerges2020.csv"), DataFrame)

    df_comp = DataFrame(
        "Compound"=>String[],
        "Experimental sarcomere shortening (fraction of control) at the EC50"=>Float64[],
        "Mechanism"=>String[],
        "Predicted scaling factor"=>Float64[]
    )

    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # assuming you have run all the simulations on the HPC

    # For each compound (positive inotrope) we retrieve the associated simulated mechanism perturbation
    for cmpd in compounds
        equals_cmpd(name) = (name == cmpd)
        df1_cmpd = filter(:Compound => equals_cmpd, df1)
        upper_asymptote = df1_cmpd.T
        y_query = 0.5 * (1 + upper_asymptote[1] / 100.0)

        df2_cmpd = filter(:Compound => equals_cmpd, df2)
        mech = df2_cmpd[!, "Mechanism"][1]

        path_to_mech = joinpath(path_to_outsim, mech)
        df_moa = CSV.read(joinpath(path_to_mech, "LR", "params_stats.csv"), DataFrame; header=1)
        x = eval(Meta.parse(df_moa[1, :domain]))
        param = df_moa[1, :param]
        deg = 2
        coefficients = [df_moa[1, "x^$(i)"] for i in 0:deg]
        poly = Polynomial(coefficients)

        path_to_cmpd = joinpath(path_to_mech, cmpd)
        mkpath(path_to_cmpd)

        # We calculate the perturbation (scaling factor) needed for a mechanism to achieve the x% increase in tension observed experimentally
        x_query = 0
        try
            fig, x_query = with_theme(() -> plotFigure(x, param, y_query, poly), theme_minimal())
            save(joinpath(path_to_cmpd, "root.pdf"), fig, pt_per_unit=1)
        catch
            println("Root finding failed!")
        end
        
        push!(df_comp, (cmpd, round(y_query, digits=2), mech, round(x_query, digits=2)))
    end

    CSV.write(joinpath(path_to_fitted_params, "pr_comparison_summary.csv"), df_comp)
end

main()
