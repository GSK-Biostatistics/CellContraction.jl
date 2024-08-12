module Utils
export checkAP, checkSS, getConverging, reduceByMedianValue

using DataFrames
using Peaks
using Statistics


function checkAP(t::Vector{Float64}, V::Vector{Float64}, dVdt::Vector{Float64})::Bool
    # Constants
    peak_threshold = 10  # min voltage to be considered a peak
	# DEP
    V_DEP = 0
    i_DEP = findfirst(V .>= V_DEP)
    t_DEP = 100
    # REP
    V_REP = -60
    V_REP2 = -20
    V_REP3 = -30
    # EADs
    t_EADs = 150
    i_EADs = findfirst(t .> t_EADs)

    # Looking for peaks
    ipeak, Vpeak = findmaxima(V)

    # Looking for big peaks (V>0) after 150 ms
    big_peak = (Vpeak .> V_DEP) .& (t[ipeak] .>= t_EADs)
    peakBIG = count(big_peak) > 1 && length(Vpeak) > 1

    # Check if there are any peaks above a threshold
    big_peak2 = (Vpeak .> peak_threshold)
    big_ipeak = ipeak[big_peak2]
    i_tt = isempty(big_ipeak) ? length(t) : big_ipeak[1]
    pos_dVdt = findfirst((dVdt[i_EADs:end-1] .> 0.01) .& (t[i_EADs:end-1] .> max(t[i_tt], t_EADs)) .& (V[i_EADs:end-1] .> V_REP))
    EADs_flag = !isnothing(pos_dVdt)

    # Check conditions for different types of APs
    if maximum(V) ≤ V_DEP || t[i_DEP] > t_DEP && V[end] ≤ V_REP
        # 5) SD or 6) DF
        return false
    elseif peakBIG && (V[1] > V_REP2 && V[end] > V_REP2)
        # 4) RF
        return false
    elseif peakBIG || EADs_flag && (V[1] > V_REP2 && V[end] > V_REP2)
        # 2) EADs
        return false
    elseif EADs_flag && V[1] > V_REP3 && V[end] <= V_REP
        # 3) EADs_2
        return false
    else
        # 1) OK
        return true
    end
end


function checkSS(x::Vector{Float64})::Bool
    threshold = 5  # hardcoded!!
	return 100 * abs(x[end] - x[1]) / abs(x[1]) < threshold
end


function getConverging(df)
    df_sum = DataFrames.transform(df, [:OK_V, :OK_Ca, :OK_T] => (+) => :OK)
    idx = findall(df_sum[!, :OK] .== 3)
    return idx
end


function reduceByMedianValue(df_cmpd, threshold=1e-2, n_elems=2)
    x = collect(keys(df_cmpd))
    y_median = map(median, collect(values(df_cmpd)))

    idx_all = findall(y_median .< threshold)
    idx_end = length(x)
    if ~isnothing(idx_all) & (length(idx_all) > n_elems)
        idx_end = idx_all[n_elems]
    end
    x = x[1:idx_end]

    for conc in keys(df_cmpd)
        if conc ∉ x
            delete!(df_cmpd, conc)
        end
    end
end


end  # module
