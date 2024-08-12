using CellContraction

using CairoMakie
using CSV
using DataFrames


function plotFigure(t, Y, d, R, labels)
    size_inches = (2 * 8.27, 2/5 * 11.69)
    size_pt = 72 .* size_inches
    colors = [cgrad(:viridis, [0, 1])[z] for z ∈ range(0, 1, length=length(Y))]
    f = Figure(size=size_pt, fontsize=12)
    ax1 = Axis(f[1, 1], xlabel="Time (ms)", ylabel="Tension (kPa)")
    for (i, y) in enumerate(Y)
        lines!(ax1, t, y, color=colors[i], linewidth=2)
    end
    ax2 = Axis(f[1, 2], xlabel="log10[Dose (M)]", ylabel="Peak Tension\n(fraction of control)")
    s = []
    for i in eachindex(d)
        si = scatter!(ax2, d[i], R[i], color=colors[i])
        si.markersize = 20
        push!(s, si)
    end
    Legend(f[1, 3], s, labels, "Doses as multipliers\nof EFTPCmax")
    f
end

function main()
    df_in = CSV.read(joinpath(pwd(), "run", "locally", "data", "population.csv"), DataFrame; header=1, types=Float64)

    # Tool compounds are characterised by their half-maximal inhibitory concentration (IC50 (uM)) and Hill coefficient
    # for each ion channel they block and by their effective free therapeutic plasma concentration (EFTPCmax)
    df = CSV.read(joinpath(pwd(), "run", "locally", "data", "tool_compounds.csv"), DataFrame; header=1)

    cmpd_idx = 28  # let's simulate the compound stored in row n.28 as an example
    # println(df.Compound[cmpd_idx])  # compound name

    # Loading drug concentration values to be simulated. These are given as multipliers of the EFTPCmax
    concentrations = CSV.read(joinpath(pwd(), "run", "locally", "data", "concentrations.csv"), DataFrame; header=1).Multipliers

    n_beats = 1000
    n_curves = 1

    p0 = values(df_in[1, 2:end])  # control model - without drug
    df_out0 = initOutData()
    B0 = runSimulCompBio(p0, n_beats)
    push!(df_out0, B0)

    t = Vector{Float64}()
    Y = Vector{Any}()
    df_out = initOutData()

    for conc in concentrations
        df_in_ = copy(df_in)
        poreBlock!(df_in_, df, conc, cmpd_idx)  # scaling full population coefficients by the current drug concentration using the pore-block model 

        p = values(df_in_[1, 2:end])  # control model - scaled by drug
        p, sol = runSimulation(p, n_beats)
       
        stim_period = p[68]
        t = collect((n_beats-n_curves)*stim_period:0.01:n_beats*stim_period)
        C = getCurves(t, sol, p)
        push!(Y, C["T"])  # collecting active tension curves at each different simulated drug concentration
    
        B = computeBiomarkers(t, sol, p)
        push!(df_out, B)  # updating biomarker dataframe at each concentration

        t = t .- t[1]
    end

    println(df_out)

    EFTPCmax, _ = getCmpdProp(df, cmpd_idx)  # get effective free therapeutic plasma concentration from input dataset
    d = log10.(1e-6 * (EFTPCmax * concentrations))  # actually simulated doses
    R = df_out.T_peak / df_out0.T_peak  # selected biomarker dose-dependent response (normalised by control value)
    labels = ["$(conc)x" for conc in concentrations]

    # Plotting
    fig = with_theme(() -> plotFigure(t, Y, d, R, labels), theme_light())
    fig
    # save("Figure2.pdf", fig, pt_per_unit=1)
end

main()
