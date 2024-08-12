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

function plotFigure(f, i, j, t, Y, c, labels)
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(Y))]
    axis = Axis(f[i, j], title=c, xlabel="Time (ms)", ylabel="Current (mV/ms)")
    for (i, y) in enumerate(Y)
        lines!(axis, t, y, color=colors[i], linewidth=2, label=labels[i])
    end
    axislegend("Factors", position=:rc, framevisible=false)
    [f, axis]
end

function main()
    # As a bonus script, we show how to use mechanism scaling to simulate pure ion channel blocks.
    # For this purpose, we directly scale ion channel conductances exposed at the user-level as parameters.
    df_in = CSV.read(joinpath(pwd(), "run", "locally", "data", "population.csv"), DataFrame; header=1, types=Float64)
    n_beats = 1000
    n_curves = 1

    size_inches = (2 * 8.27, 2 * 2/5 * 11.69)
    size_pt = 72 .* size_inches
    fig = Figure(size=size_pt, fontsize=12)

    currents = ["INaL", "ICaL", "Ito", "IKr", "IKs", "IK1"]
    n_cols = 3
    axs = []
    for (k, current) in enumerate(currents)  # simulate pure block of all ion channels, one at a time
        i = (k - 1) ÷ n_cols + 1
        j = mod(k - 1, n_cols) + 1

        d_factors = Dict(replace(current, "I"=>"G")=>[0.25, 0.5, 0.75, 1.0])  # parameter controlling current 'Ix' is normally called 'Gx'
        d_combinations = genMechCombs(d_factors)
        df_in_dict = mechAlter(df_in, d_combinations)

        t = Vector{Float64}()
        Y = Vector{Vector{Float64}}()

        for comb_id in keys(df_in_dict)
            df_in_ = df_in_dict[comb_id]
            p = values(df_in_[1, 2:end])
            p, sol = runSimulation(p, n_beats)
            stim_period = p[68]
            t = collect((n_beats-n_curves)*stim_period:0.01:n_beats*stim_period)
            C = getCurrents(t, sol, p)  # similarly to 'getCurves', this function returns all the ionic currents (discrete solution) at the chosen time points
            push!(Y, C[current])
            t = t .- t[1]
        end

        # Plotting
        labels = fancy(collect(keys(df_in_dict)), length(keys(d_factors)))
        fig, axis = with_theme(() -> plotFigure(fig, i, j, t, Y, current, labels), theme_light())
        push!(axs, axis)
    end
    linkxaxes!(axs...)
    fig
    # save("Figure3_bonus.pdf", fig, pt_per_unit=1)
end

main()
