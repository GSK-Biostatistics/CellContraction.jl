using CellContraction

using CairoMakie
using CSV
using DataFrames


function plotFigure(t, Y, labels)
    size_inches = (2 * 8.27, 2/5 * 11.69)  # a proportion of the A4 paper size
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=12)
    c = Makie.wong_colors()[1:length(Y)]
    for (i, y) in enumerate(Y)
        axis = Axis(f[1, i], xlabel="Time (ms)", ylabel=labels[i])
        lines!(axis, t, y, color=c[i], linewidth=3)
    end
    f
end

function main()
    # Each row of this dataset is a set of scaling coefficients for the ion channels.
    # Different scaling coefficients values encode ion channel composition variability among different cardiac cells,
    # so each row of this dataset constitutes a different model. The whole dataset thus represents a population of models.
    df_in = CSV.read(joinpath(pwd(), "run", "locally", "data", "population.csv"), DataFrame; header=1, types=Float64)
    inputs = [values(row) for row in eachrow(df_in[!, 2:end])]

    p = inputs[1]  # control model scaling coefficients are stored in the first row
    n_beats = 1000  # number of beats to simulate to reach a steady-state
    p, sol = runSimulation(p, n_beats)  # running control simulation (this also returns again the parameters used)

    n_curves = 1  # choosing only the last-beat curve to be visualised and to extract biomarkers from
    stim_period = p[68]  # needed to reconstruct the time vector
    t = collect((n_beats-n_curves)*stim_period:0.01:n_beats*stim_period)  # reconstructing time vector with resolution of 1e-2 milliseconds
    C = getCurves(t, sol, p)  # returning solution dictionary with curves (discrete) from a dense solution which can interpolate at any time point
    println(names(C))  # available solution curves

    df_out = initOutData()  # initialising empty dataframe for output data
    B = computeBiomarkers(t, sol, p)  # calculating biomarker values from dense solution
    push!(df_out, B)  # updating output dataframe with calculated values
    println(df_out)

    t = t .- t[1]  # limit cycle does not start from 0 milliseconds: now it does
    Y = [C[key] for key in ["V", "Ca", "T"]]
    labels = ["Transmembrane Voltage (mV)", "Intracellular Calcium (nM)", "Active Tension (kPa)"]

    # Plotting
    fig = with_theme(() -> plotFigure(t, Y, labels), theme_light())
    fig  # output solution plot
    # save("Figure1.pdf", fig, pt_per_unit=1)
end

main()
