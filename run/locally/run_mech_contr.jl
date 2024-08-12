using CellContraction

using CairoMakie
using CSV
using DataFrames


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

function plotFigure1D(t, Y, d, R, labels, n_params)
    size_inches = (2 * 8.27, 2/5 * 11.69)
    size_pt = 72 .* size_inches
    m = split(labels[1], "_")[2]
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(Y))]
    f = Figure(size=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], xlabel="Time (ms)", ylabel="Tension (kPa)")
    for (i, y) in enumerate(Y)
        lines!(ax1, t, y, color=colors[i], linewidth=2)
    end
    ax2 = Axis(f[1, 2], xlabel="$(m) scaling factor", ylabel="Peak Tension\n(fraction of control)")
    s = []
    dx = map(x -> x[1], d)
    for i in eachindex(dx)
        si = scatter!(ax2, dx[i], R[i], color=colors[i])
        si.markersize = 20
        push!(s, si)
    end
    Legend(f[1, 3], s, fancy(labels, n_params), "Factors")
    f
end

function plotFigure2D(t, Y, d, R, labels, n_params)
    size_inches = (2 * 8.27, 2/5 * 11.69)
    size_pt = 72 .* size_inches
    m1 = split(labels[1], "_")[2]
    m2 = split(labels[1], "_")[4]
    f = Figure(size=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], xlabel="Time (ms)", ylabel="Tension (kPa)")
    l = []
    for y in Y
        li = lines!(ax1, t, y, linewidth=2)
        push!(l, li)
    end
    Legend(f[1, 2], l, fancy(labels, n_params), "Factors")
    ax3 = Axis(f[1, 3], xlabel="$(m1) scaling factor", ylabel="$(m2) scaling factor")
    dx = map(x -> x[1], d)
    dy = map(x -> x[2], d)
    hm = heatmap!(ax3, dx, dy, R, interpolate=true)  # if points in dx and points in dy are not equally-spaced, use interpolate=false
    Colorbar(f[1, 4], hm, label="Peak Tension\n(fraction of control)")
    f
end

function plotFigureND(t, Y, labels, n_params)
    size_inches = (8.27, 2/5 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], xlabel="Time (ms)", ylabel="Tension (kPa)")
    l = []
    for y in Y
        li = lines!(ax1, t, y, linewidth=2)
        push!(l, li)
    end
    Legend(f[1, 2], l, fancy(labels, n_params), "Factors")
    f
end

function main()
    # There are two possible ways of applying a specific perturbation to a/multiple mechanism/s
    # (1) Read mechanism from .csv file using row number (excluding the header, so for mechanism in row 2 use mech_idx=1)
    df = CSV.read(joinpath(pwd(), "run", "locally", "data", "mechanisms.csv"), DataFrame; header=1, types=String)
    mech_idx = 1
    d_factors = getMechScaling(df, mech_idx)

    # (2) Write by hand a dictionary containing the scaling factors for the specified mechanism(s)
    # Note: it is possible to specify as many mechanisms as needed.
    # However, keep in mind that this will result in a high number of factors' combinations (parameter sets) to be simulated!
    # d_factors = Dict("Ca50"=>[0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])  # 1D example
    # d_factors = Dict("GCaL"=>[1.5, 2.0, 2.5], "Cao"=>[0.8, 1.0, 1.2])  # 2D example
    # d_factors = Dict("Gncx"=>[1.5, 1.7], "GCaL"=>[0.8, 1.2], "GNa"=>[0.8, 0.9])  # ND (N=3) example

    df_in = CSV.read(joinpath(pwd(), "run", "locally", "data", "population.csv"), DataFrame; header=1, types=Float64)
    d_combinations = genMechCombs(d_factors)  # generate all combinations from given scaling factors
    df_in_dict = mechAlter(df_in, d_combinations)  # apply scaling to population coefficients' matrix, this will return as many new matrices as the number of combinations

    n_beats = 1000
    n_curves = 1

    p0 = values(df_in[1, 2:end])  # control model - mechanism(s) unperturbed
    df_out0 = initOutData()
    B0 = runSimulCompBio(p0, n_beats)
    push!(df_out0, B0)

    t = Vector{Float64}()
    Y = Vector{Vector{Float64}}()
    df_out = initOutData()

    for comb_id in keys(df_in_dict)
        df_in_ = df_in_dict[comb_id]

        p = values(df_in_[1, 2:end])  # control model - mechanism(s) perturbed
        p, sol = runSimulation(p, n_beats)
    
        stim_period = p[68]
        t = collect((n_beats-n_curves)*stim_period:0.01:n_beats*stim_period)
        C = getCurves(t, sol, p)
        push!(Y, C["T"])

        B = computeBiomarkers(t, sol, p)
        push!(df_out, B)

        t = t .- t[1]
    end

    println(df_out)

    d = map(x -> collect(values(x)), values(d_combinations))  # mechanism(s) scaling factors
    R = df_out.T_peak ./ df_out0.T_peak  # selected biomarker response to mechanism(s) perturbation (normalised by control value)
    labels = collect(keys(df_in_dict))
    n_params = length(keys(d_factors))

    # Plotting
    theme = theme_light()
    if n_params == 1
        fig = with_theme(() -> plotFigure1D(t, Y, d, R, labels, n_params), theme)
    elseif n_params == 2
        fig = with_theme(() -> plotFigure2D(t, Y, d, R, labels, n_params), theme)
    else
        fig = with_theme(() -> plotFigureND(t, Y, labels, n_params), theme)
    end
    fig
    # save("Figure3.pdf", fig, pt_per_unit=1)
end

main()
