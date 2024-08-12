module Runner
export computeBiomarkers, getCurrents, getCurves, runSimulation, runSimulCompBio, initOutData

using ForceImport
@force using ..Currents
@force using ..Initializer
@force using ..Solver
@force using ..Utils

using DataFrames
using OrderedCollections


function active_tension(U::Any, dU::Any, p::Vector{Float64})::Vector{Vector{Float64}}
    lambda     = p[34]
    dr         = p[39]
    beta_0     = p[44]
    Tref       = p[50]
    lambda_min = p[55]
    lambda_max = p[56]

    XS    = U[3, :]
    XW    = U[4, :]
    ZETAS = U[5, :]
    ZETAW = U[6, :]

    dXS    = dU[3, :]
    dXW    = dU[4, :]
    dZETAS = dU[5, :]
    dZETAW = dU[6, :]

    lambda0 = min(lambda_max, lambda)
    Lfac = max(0, 1 + beta_0 * (lambda0 + min(lambda_min, lambda0) - (1 + lambda_min)))
    
    T = Lfac * (Tref/dr) * ((1 .+ ZETAS) .* XS + ZETAW .* XW)
    dTdt = Lfac * (Tref/dr) * ((1 .+ ZETAS) .* dXS + dZETAS .* XS + ZETAW .* dXW + dZETAW .* XW)
    [T, dTdt]
end


function currents(U::Any, p::Vector{Float64})::Vector{Vector{Float64}}
    v      = U[1, :]
    nai    = U[2, :]
    nass   = U[3, :]
    ki     = U[4, :]
    kss    = U[5, :]
    cai    = U[6, :]
    cass   = U[7, :]
    mL     = U[8, :]
    hL     = U[9, :]
    hLp    = U[10, :]
    a      = U[11, :]
    iF     = U[12, :]
    iS     = U[13, :]
    ap     = U[14, :]
    iFp    = U[15, :]
    iSp    = U[16, :]
    d      = U[17, :]
    ff     = U[18, :]
    fs     = U[19, :]
    fcaf   = U[20, :]
    fcas   = U[21, :]
    jca    = U[22, :]
    nca    = U[23, :]
    nca_i  = U[24, :]
    ffp    = U[25, :]
    fcafp  = U[26, :]
    xs1    = U[27, :]
    xs2    = U[28, :]
    CaMKt  = U[29, :]
    ikr_c0 = U[30, :]
    ikr_c1 = U[31, :]
    ikr_c2 = U[32, :]
    ikr_o  = U[33, :]
    ikr_i  = U[34, :]
    cli    = U[35, :]
    clss   = U[36, :]

    cao         = p[1]
    nao         = p[2]
    ko          = p[3]
    clo         = p[4]
    R           = p[5]
    T           = p[6]
    F           = p[7]
    ICaL_fracSS = p[17]
    PKNa        = p[28]
    KmCaMK      = p[29]
    CaMKo       = p[32]
    KmCaM       = p[33]
    GNaL        = p[70]
    Gto         = p[71]
    PCa         = p[72]
    GKr         = p[73]
    GKs         = p[74]
    GK1         = p[75]

    vffrt = v .* (F * F / (R * T))
    vfrt = v .* (F / (R * T))

    ENa = (R * T / F) .* log.(nao ./ nai)
    EK = (R * T / F) .* log.(ko ./ ki)
    EKs = (R * T / F) .* log.((ko + PKNa * nao) ./ (ki .+ PKNa .* nai))

    CaMKb = CaMKo .* (1.0 .- CaMKt) ./ (1.0 .+ KmCaM ./ cass)
    CaMKa = CaMKb .+ CaMKt
    fINaLp = (1.0 ./ (1.0 .+ KmCaMK ./ CaMKa))
    fItop = (1.0 ./ (1.0 .+ KmCaMK ./ CaMKa))
    fICaLp = (1.0 ./ (1.0 .+ KmCaMK ./ CaMKa))

    INaL_out = getINaL_ORd2011.(v, mL, hL, hLp, fINaLp, ENa, GNaL)
    INaL = [y[1] for y in INaL_out]
    ICaL_out = getICaL_ORd2011_jt.(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, clss, ICaL_fracSS, vffrt, vfrt, PCa)
    ICaL_ss = [y[1] for y in ICaL_out]
    ICaL_i = [y[4] for y in ICaL_out]
    ICaL = ICaL_ss .+ ICaL_i
    Ito_out = getITo_ORd2011.(v, a, iF, iS, ap, iFp, iSp, fItop, EK, Gto)
    Ito = [y[1] for y in Ito_out]
    IKr_out = getIKr_ORd2011_MM.(v, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i, ko, EK, vfrt, GKr)
    IKr = [y[1] for y in IKr_out]
    IKs_out = getIKs_ORd2011.(v, xs1, xs2, cai, EKs, GKs)
    IKs = [y[1] for y in IKs_out]
    IK1 = getIK1_CRLP.(v, ko, EK, GK1)
    [INaL, ICaL, Ito, IKr, IKs, IK1]
end


function compute_repolarisation_time_p_percent(t::Vector{Float64}, x::Vector{Float64}, p::Int64)::Float64
    m = x[1]
    M = maximum(x)
    threshold = m + (1 - p/100) * (M - m)
    idx = findlast(x .>= threshold)
    APD = t[end]
    if !isnothing(idx)
        APD = t[idx] - t[1]
    end
end


function computeBiomarkers(t::Vector{Float64}, sol::Any, p::Vector{Float64})::Vector{Float64}
    V, dVdt, Ca, T, dTdt = values(getCurves(t, sol, p))
    Vttp     = t[argmax(V)] - t[1]
    RMP      = minimum(V)
    V_peak   = maximum(V)
    APD40    = compute_repolarisation_time_p_percent(t, V, 40) - Vttp
    APD90    = compute_repolarisation_time_p_percent(t, V, 90) - Vttp
    Tri9040  = APD90 - APD40
    dVdt_max = maximum(dVdt)
    qNet     = sol(t, idxs=52).u[end] - sol(t, idxs=52).u[1]
    Cattp    = t[argmax(Ca)] - t[1]
    Ca_diast = minimum(Ca)
    Ca_peak  = maximum(Ca)
    Ca_ampl  = Ca_peak - Ca_diast
    CTD50    = compute_repolarisation_time_p_percent(t, Ca, 50) - Cattp
    CTD90    = compute_repolarisation_time_p_percent(t, Ca, 90) - Cattp
    EMw90    = CTD90 - APD90
    Tttp     = t[argmax(T)] - t[1]
    T_peak   = maximum(T)
    Trt50    = compute_repolarisation_time_p_percent(t, T, 50) - Tttp
    Trt90    = compute_repolarisation_time_p_percent(t, T, 90) - Tttp
    dTdt_max = maximum(dTdt)
    dTdt_min = minimum(dTdt)
    OK_V     = checkAP(t.-t[1], V, dVdt)
    OK_Ca    = checkSS(Ca)
    OK_T     = checkSS(T)
    [RMP, V_peak, APD40, APD90, Tri9040, dVdt_max, qNet, Ca_diast, Ca_peak, Ca_ampl, CTD50, CTD90, EMw90, Tttp, T_peak, Trt50, Trt90, dTdt_max, dTdt_min, OK_V, OK_Ca, OK_T]
end


function runSimulation(scaling_coefficients::NTuple{14, Float64}, n_beats::Int64)
    p = initConstants(scaling_coefficients)
    sol = solveEP(n_beats, p)
    return p, sol
end


function getCurves(t::Vector{Float64}, sol::Any, p::Vector{Float64})::OrderedDict{String, Vector{Float64}}
    U = sol(t, idxs=[1, 6, 44, 45, 48, 49])
    dU = sol(t, Val{1}, idxs=[1, 6, 44, 45, 48, 49])
    V, dVdt = U[1, :], dU[1, :]
    Ca = U[2, :] * 1e6
    T, dTdt = active_tension(U, dU, p)
    return OrderedDict("V"=>V, "dVdt"=>dVdt, "Ca"=>Ca, "T"=>T, "dTdt"=>dTdt)
end


function getCurrents(t::Vector{Float64}, sol::Any, p::Vector{Float64})::OrderedDict{String, Vector{Float64}}
    U = sol(t, idxs=[1, 2, 3, 4, 5, 6, 7, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 50, 51])
    INaL, ICaL, Ito, IKr, IKs, IK1 = currents(U, p)
    return OrderedDict("INaL"=>INaL, "ICaL"=>ICaL, "Ito"=>Ito, "IKr"=>IKr, "IKs"=>IKs, "IK1"=>IK1)
end


function runSimulCompBio(scaling_coefficients::NTuple{14, Float64}, n_beats::Int64)::Vector{Float64}
    bio = zeros(Float64, 22)
    try
        p, sol = runSimulation(scaling_coefficients, n_beats)
        stim_period = p[68]
        t = collect((n_beats-1)*stim_period:0.01:n_beats*stim_period)
        bio = computeBiomarkers(t, sol, p)
    catch e
        println("ODE solver failed due to the following error:\n$e")
    end
    return bio
end


function initOutData()
    df_out = DataFrame(
        RMP=Float64[],
        V_peak=Float64[],
        APD40=Float64[],
        APD90=Float64[],
        Tri9040=Float64[],
        dVdt_max=Float64[],
        qNet=Float64[],
        Ca_diast=Float64[],
        Ca_peak=Float64[],
        Ca_ampl=Float64[],
        CTD50=Float64[],
        CTD90=Float64[],
        EMw90=Float64[],
        Tttp=Float64[],
        T_peak=Float64[],
        Trt50=Float64[],
        Trt90=Float64[],
        dTdt_max=Float64[],
        dTdt_min=Float64[],
        OK_V=Float64[],
        OK_Ca=Float64[],
        OK_T=Float64[]
    )
    return df_out
end


end  # module
