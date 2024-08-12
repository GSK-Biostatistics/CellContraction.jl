module Drug
export getCmpdProp, getCmpdPropLess, getDoses, genMechCombs, getMechScaling, mechAlter, poreBlock!

using DataFrames
using OrderedCollections


function getCmpdProp(df, cmpd_idx)
    cmpd = df.Compound[cmpd_idx]
    row = findfirst(==(cmpd), df.Compound)
    channels = map(x -> collect(eachsplit(x, "_"))[1], names(df)[3:2:end])
    d = Dict()
    for cha in channels
        cha_columns = names(df, x -> startswith(x, cha*"_"))
        Hill, IC50 = values(df[row, cha_columns])
        if !isnan(IC50)
            push!(d, cha=>Dict("IC50"=>IC50, "Hill"=>Hill))
        end
    end
    EFTPCmax = df.EFTPCmax[row]
    return EFTPCmax, d
end


function getCmpdPropLess(df, cmpd)
    row = findfirst(==(cmpd), df.Compound)
    EFTPCmax = df.EFTPCmax[row]
    return EFTPCmax
end


function getDoses(EFTPCmax, x)
    return log10.(1e-6 * (EFTPCmax * x))
end


function sigmoid(C::Float64, IC50::Float64, Hill::Float64)::Float64
    factor = (1 + (C/IC50) ^ Hill) ^ (-1)
    return factor
end


function poreBlock!(df_in, df, conc, cmpd_idx)
    EFTPCmax, IC50_Hill_dct = getCmpdProp(df, cmpd_idx)
    for cha in keys(IC50_Hill_dct)
        df_tmp = mapcols(x -> x*sigmoid(EFTPCmax*conc, IC50_Hill_dct[cha]["IC50"], IC50_Hill_dct[cha]["Hill"]), df_in)
        df_in[!, cha] = df_tmp[!, cha]
    end
    nothing
end


function getMechScaling(df, mech_idx)
    cmpd = df.Mechanism[mech_idx]
    row = findfirst(==(cmpd), df.Mechanism)
    mechanisms = names(df)[2:end]
    d_factors = Dict{String, Vector{Float64}}()
    for mec in mechanisms
        val = values(df[row, mec])
        if !contains(val, "NaN")
            factors = parse.(Float64, split(val[2:end-1], ", "))
            push!(d_factors, mec=>factors)
        end
    end
    return d_factors
end


function genMechCombs(d_factors)
    function generate_combinations(lists)
        return collect(Base.Iterators.product(lists...))
    end
    combinations = generate_combinations(values(d_factors))
    d_combinations = OrderedDict()
    for (comb_id, comb) in enumerate(combinations)
        factors = Dict(mec=>val for (mec, val) in zip(keys(d_factors), comb))
        push!(d_combinations, comb_id=>factors)
    end
    return d_combinations
end


function mechAlter(df_in, d_combinations)
    d_data = OrderedDict()
    for comb_id in keys(d_combinations)
        str = "m" * string(comb_id)
        df_in_ = copy(df_in)
        factor = d_combinations[comb_id]
        for mec in keys(factor)
            str = str * "_" * mec * "_$(factor[mec])"
            df_tmp = mapcols(x -> x*factor[mec], df_in_)
            df_in_[!, mec] = df_tmp[!, mec]
        end
        push!(d_data, str=>df_in_)
    end
    return d_data
end


end  # module
