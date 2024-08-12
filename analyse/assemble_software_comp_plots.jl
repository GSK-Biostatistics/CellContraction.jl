using CairoMakie
using CSV
using DataFrames
using Glob
using OrderedCollections


function main()
    path_to_analysis = joinpath(pwd(), "analyse")

    path_to_data = joinpath(path_to_analysis, "data", "julia_matlab_va_datasets")
    files = glob("*.csv", path_to_data)
    df = OrderedDict()
    for f in files
        println(f)
        name = split(split(f, "\\")[end], ".csv")[1]
        df[name] = CSV.read(f, DataFrame; header=1, types=Float64)
    end
    
    biomarkers = names(collect(values(df))[1])

    add_units = OrderedDict(
       "RMP"=>"RMP (mV)",
       "V_peak"=>"Vpeak (mV)",
       "APD40"=>"APD40 (ms)",
       "APD90"=>"APD90 (ms)",
       "dVdt_max"=>"dVdtMax (mV/ms)",
       "Ca_diast"=>"CaiD (nM)",
       "Ca_peak"=>"CTpeak (nM)",
       "CTD50"=>"CTD50 (ms)",
       "CTD90"=>"CTD90 (ms)",
       "Tttp"=>"ATttp (ms)",
       "T_peak"=>"ATpeak (kPa)",
       "Trt50"=>"ATrt50 (ms)",
       "Trt90"=>"ATrt90 (ms)",
       "dTdt_max"=>"dATdtMax (kPa/ms)",
       "dTdt_min"=>"dATdtMin (kPa/ms)",
       "EMw90"=>"EMw (ms)",
       "Tri9040"=>"Tri9040 (ms)",
       "qNet"=>"qNet (μC/μF)"
    )

    size_inches = (2 * 8.27, 2 * 0.75 * 11.69)
    size_pt = 72 .* size_inches
    f = Figure(size=size_pt, fontsize=18)

    xpos = [0.95, 1.0, 1.05]
    axes = []
    dataset_names = keys(df)
    
    n_cols = 3
    for (k, b) in enumerate(keys(add_units))
        i = (k - 1) ÷ n_cols + 1
        j = mod(k - 1, n_cols) + 1
        axis = Axis(f[i, j], ylabel=add_units[b], xticks=(xpos, [rich(si, color=ci, font=:bold) for (si, ci) in zip(replace.(dataset_names, "-"=>" "), Makie.wong_colors()[1:length(xpos)])]))
        push!(axes, axis)
        for (c, name) in enumerate(dataset_names)
            y = df[name][!, b]
            boxplot!(axis, fill(xpos[c], length(y)), y, width=0.05, color=Makie.wong_colors()[c], show_outliers=true)
        end
        if i != Int(length(biomarkers) / n_cols)
            hidexdecorations!(axis, grid=false, ticks=false)
        end
    end
    linkxaxes!(axes...)

    path_to_output = joinpath(path_to_analysis, "output")
    save(joinpath(path_to_output, "software_comparison.png"), f, px_per_unit=5)
end

main()
