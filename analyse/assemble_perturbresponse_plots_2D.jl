using CellContraction

using CairoMakie
using CSV
using DataFrames
using Glob
using Statistics


function fancy(labels, n_params)
    for (i, lab) in enumerate(labels)
        k = [split(lab, "_")[2*i] for i=1:n_params]
        v = [split(lab, "_")[2*i+1] for i=1:n_params]
        lab_new = ""
        for (ki, vi) in zip(k, v)
            lab_new = lab_new * ki * "=$(vi), "
        end
        labels[i] = lab_new[1:end-2]
    end
    labels
end


function plotFigure2D(d, R_mean, R_std, labels)
    size_inches = (2 * 8.27, 2 * 0.225 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=18)

    m1 = split(labels[1], "_")[2]
    m2 = split(labels[1], "_")[4]
    dx = map(x -> x[1], d)
    dy = map(x -> x[2], d)

    ax1 = Axis(f[1, 1], xlabel="Scaling Factor ($(m1))", ylabel="Scaling Factor ($(m2))")
    hm1 = heatmap!(ax1, dx, dy, R_mean, interpolate=true, colormap=:Blues)
    # contour!(ax1, dx, dy, R_mean, levels = [1.71, 1.93, 2.58], color=:red)
    Colorbar(f[1, 2], hm1, label=rich("AT", subscript("peak"), "  Mean\n(fraction of contr.)"))
    ax2 = Axis(f[1, 4], xlabel="Scaling Factor ($(m1))", ylabel="Scaling Factor ($(m2))")
    hm2 = heatmap!(ax2, dx, dy, R_std, interpolate=true, colormap=:Blues)
    Colorbar(f[1, 5], hm2, label=rich("AT", subscript("peak"), "  STD\n(fraction of contr.)"))
    Label(f[0, :], "Beta_increase (GKs + GCaL)", fontsize=18, font=:bold)
    colsize!(f.layout, 3, Auto(0.05))
    f
end


function main()
    path_to_outsim = joinpath(pwd(), "run", "hpc", "output")  # folder containing all the machanism perturbation simulations
    df0 = CSV.read(joinpath(path_to_outsim, "biomarkers.csv"), DataFrame; header=1, types=Float64)
  
    mech = "Beta_increase"  # the only mechanism that can be really plotted as a heatmap given its 2D nature
    path_to_mech = joinpath(path_to_outsim, mech)
    
    df = Dict()
    files = glob("*.csv", path_to_mech)
    for f in files
        idx_param_value = split(split(f, "biomarkers_m")[2], ".csv")[1]
        df[idx_param_value] = CSV.read(f, DataFrame; header=1, types=Float64)
    end
    
    x = collect(keys(df))
    sort!(x, by=x->parse(Int64, split(x, "_")[1]))

    idx = [3, 5]
    d = [parse.(Float64, split(xi, "_")[idx]) for xi in x]

    idx_list = []
    for xi in x
        idx = getConverging(df[xi])
        push!(idx_list, idx)
    end
    idx_shared = intersect(idx_list...)

    R_mean = [mean(df[xi][idx_shared, :].T_peak ./ median(df0[idx_shared, :].T_peak)) for xi in x]
    R_std = [std(df[xi][idx_shared, :].T_peak ./ median(df0[idx_shared, :].T_peak)) for xi in x]
    
    labels = ["m$(xi)" for xi in x]
    f = plotFigure2D(d, R_mean, R_std, labels)
    
    path_to_output = joinpath(pwd(), "analyse", "output")
    save(joinpath(path_to_output, "T_peak_$(mech).png"), f, px_per_unit=5)
end

main()
