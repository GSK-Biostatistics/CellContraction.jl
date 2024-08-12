using CellContraction

using CairoMakie
using CSV
using DataFrames
using Glob
using OrderedCollections
using Statistics


function plotFigure(f, i, j, df_cmpd, EFTPCmax, n_doses, cmpd, pIC50)
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(keys(df_cmpd)))]
    axis = Axis(f[i, j], title="$(cmpd)", xticks=[-12, -10, -8, -6, -4, -2, 0])
    boxes = []
    elements = []
    for (i, key) in enumerate(keys(df_cmpd))
        if i > n_doses
            break
        end
        b = boxplot!(axis, fill(getDoses(EFTPCmax, key), length(df_cmpd[key])), df_cmpd[key], width=0.4, color=colors[i], show_outliers=false)
        push!(boxes, b)
    end
    if pIC50[1] == "line"
        vl = vlines!(pIC50[2], color=:red, linewidth=0.8)
        push!(elements, vl)
    else
        vs = vspan!(pIC50[2], 1.5, color=(:red, 0.15))
        push!(elements, vs)
    end
    xlims!(-12.5, 1.5)
    num_to_str(x) = string(x == floor(x) ? Int(x) : x)
    labels = map(x->num_to_str(x)*"x", collect(keys(df_cmpd)))
    [f, axis, labels[1:n_doses], boxes, elements[1]]
end

function main()
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # note: in this folder you must have only drug block simulations, NOT mechanism perturbation simulations
    compounds = setdiff(readdir(path_to_outsim), ["biomarkers.csv"])

    df0 = CSV.read(joinpath(path_to_outsim, "biomarkers.csv"), DataFrame; header=1, types=Float64)

    path_to_data = joinpath(pwd(), "run", "hpc", "data")
    df_props = CSV.read(joinpath(path_to_data, "tool_compounds.csv"), DataFrame; header=1)

    path_to_analysis = joinpath(pwd(), "analyse")
    
    path_to_exp = joinpath(path_to_analysis, "data")
    df_exp = CSV.read(joinpath(path_to_exp, "Nguyen2017.csv"), DataFrame; header=1)

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

    at_least_models_n = 50
    bio = "T_peak"

    df_new = OrderedDict()
    conc_max = concentrations[1]
    idx_shared = collect(1:size(df0)[1])
    for cmpd in compounds
        df_new[cmpd] = OrderedDict()
        for conc in reverse(concentrations)
            idx = getConverging(df[cmpd][conc])
            if length(idx) > at_least_models_n
                conc_max = conc
                idx_shared = idx
                break
            end
        end
        for conc in concentrations
            if conc <= conc_max
                push!(df_new[cmpd], conc => df[cmpd][conc][idx_shared, :][!, Symbol(bio)] ./ median(df0[!, Symbol(bio)]))
            end
        end
    end

    size_inches = (2 * 8.27, 2 * 0.65 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=14)
    
    axs = []
    labs = []
    boxs = []
    elems = []

    n_cols = 4
    for (k, cmpd) in enumerate(compounds)
        i = (k - 1) ÷ n_cols + 1
        j = mod(k - 1, n_cols) + 1

        equals_cmpd(name) = (name == cmpd)
        df_exp_cmpd = filter(:Compound => equals_cmpd, df_exp)

        ic50_split = collect(eachsplit(df_exp_cmpd[1, "IC50"], ">"))
        if length(ic50_split) > 1
            pIC50 = ("vspan", log10(1e-6 * parse(Float64, ic50_split[2])))
        else
            pIC50 = ("line", log10(1e-6 * parse(Float64, ic50_split[1])))
        end

        df_cmpd = df_new[cmpd]
        EFTPCmax = getCmpdPropLess(df_props, cmpd)
        x = getDoses(EFTPCmax, collect(keys(df_cmpd)))
        Y = collect(values(df_cmpd))
        y_median = dropdims(median(reduce(hcat, Y), dims=1), dims=1)
        idx_all = findall(y_median .< 1e-2)
        idx_end = length(x)
        if ~isnothing(idx_all) & (length(idx_all) > 2)
            idx_end = idx_all[2]
        end
        x = x[1:idx_end]
        f, axis, labels, boxes, elements = plotFigure(f, i, j, df_cmpd, EFTPCmax, length(x), cmpd, pIC50)
        if i != Int(length(compounds) / n_cols)
            hidexdecorations!(axis, grid=false, ticks=false)
        end
        if j != 1
            hideydecorations!(axis, grid=false, ticks=false)
        end
        push!(axs, axis)
        push!(labs, labels)
        push!(boxs, boxes)
        push!(elems, elements)
    end

    elems = elems[1:2]
    linkaxes!(axs...)
    Label(f[0, :], rich("AT", subscript("peak"), " | x-axis: log", subscript("10"), "(Drug Concentration [M]) | y-axis: Fraction of Control"), fontsize=20)
    idx = sortperm(labs, by=length)[end]

    lab1 = "Could not observe at least 25%\ndecrease in sarcomere shortening\nwhen dosing outside this region"
    lab2 = "Observed IC50\n(sarcomere shortening)"
    t1 = rich("Doses as multipliers\nof EFTPC", subscript("max"))
    t2 = "Experimental Data"
    Legend(
        f[end+1, :],
        [boxs[idx], elems],
        [labs[idx], [lab1; lab2]],
        [t1, t2],
        nbanks=2,
        titleposition=:left,
        orientation=:horizontal,
        tellheight=true,
        tellwidth=false,
        framevisible=false
    )
    rowsize!(f.layout, 8, Relative(0.1))

    path_to_output = joinpath(path_to_analysis, "output")
    save(joinpath(path_to_output, "$(bio)_all_negative_inotropes.png"), f, px_per_unit=5)
end

main()
