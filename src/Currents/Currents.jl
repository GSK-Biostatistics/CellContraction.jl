module Currents
export getINa_Grandi, getINaL_ORd2011, getITo_ORd2011, getICaL_ORd2011_jt, getIKr_ORd2011_MM, getIKs_ORd2011, getIK1_CRLP, getINaK_ORd2011, getINaCa_ORd2011, getJrel_ORd2011


# INa
function getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, GNa)
    mss = 1 / ((1 + exp( -(56.86 + v) / 9.03 ))^2)
    taum = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2)
    dm = (mss - m) / taum
    ah = (v >= -40) * (0) + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 ))
    bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1e5 * exp(0.3485 * v)))
    tauh = 1 / (ah + bh)
    hss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2)
    dh = (hss - h) / tauh
    aj = (v >= -40) * (0) + (v < -40) * (((-2.5428e4*exp(0.2444*v) - 6.948e-6 * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )))
    bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )))
    tauj = 1 / (aj + bj)
    jss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2)
    dj = (jss - j) / tauj
    hssp = 1 / ((1 + exp( (v + 71.55 + 6)/7.43 ))^2)
    dhp = (hssp - hp) / tauh
    taujp = 1.46 * tauj
    djp = (jss - jp) / taujp
    INa = GNa*(v-ENa)*m^3.0*((1.0-fINap)*h*j+fINap*hp*jp)
    return INa, dm, dh, dhp, dj, djp
end

# INaL
function getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, GNaL)
    mLss = 1.0 / (1.0 + exp((- (v + 42.85)) / 5.264))
    tm = 0.1292 * exp(-((v + 45.79) / 15.54) ^ 2) + 0.06487 * exp(-((v - 4.823) / 51.12) ^ 2)
    tmL = tm
    dmL = (mLss - mL) / tmL
    hLss = 1.0 / (1.0 + exp((v + 87.61) / 7.488))
    thL = 200.0
    dhL = (hLss - hL) / thL
    hLssp = 1.0 / (1.0 + exp((v + 93.81) / 7.488))
    thLp = 3.0 * thL
    dhLp = (hLssp - hLp) / thLp
    INaL = GNaL * (v - ENa) * mL * ((1.0 - fINaLp) * hL + fINaLp * hLp)
    return INaL, dmL, dhL, dhLp
end

# Ito
function getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, Gto)
    ass = 1.0 / (1.0 + exp((-(v - 14.34)) / 14.82))
    ta = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(v - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((v + 100.0) / 29.3814)))
    da = (ass - a) / ta
    iss = 1.0 / (1.0 + exp((v + 43.94) / 5.711))
    delta_epi = 1.0
    tiF = 4.562 + 1 / (0.3933 * exp((-(v + 100.0)) / 100.0) + 0.08004 * exp((v + 50.0) / 16.59))
    tiS = 23.62 + 1 / (0.001416 * exp((-(v + 96.52)) / 59.05) + 1.780e-8 * exp((v + 114.1) / 8.079))
    tiF *= delta_epi
    tiS *= delta_epi
    AiF = 1.0 / (1.0 + exp((v - 213.6) / 151.2))
    AiS = 1.0 - AiF
    diF = (iss - iF) / tiF
    diS = (iss - iS) / tiS
    i = AiF * iF + AiS * iS
    assp = 1.0 / (1.0 + exp((-(v - 24.34)) / 14.82))
    dap = (assp - ap) / ta
    dti_develop = 1.354 + 1.0e-4 / (exp((v - 167.4) / 15.89) + exp(-(v - 12.23) / 0.2154))
    dti_recover = 1.0 - 0.5 / (1.0 + exp((v + 70.0) / 20.0))
    tiFp = dti_develop * dti_recover * tiF
    tiSp = dti_develop * dti_recover * tiS
    diFp = (iss - iFp) / tiFp
    diSp = (iss - iSp) / tiSp
    ip = AiF * iFp + AiS * iSp
    Ito = Gto * (v - EK) * ((1.0 - fItop) * a * i + fItop * ap * ip)
    return Ito, da, diF, diS, dap, diFp, diSp
end

# ICaL
function getICaL_ORd2011_jt(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, clss, ICaL_fracSS, vffrt, vfrt, PCa)
    dss = 1.0763 * exp(-1.0070 * exp(-0.0829 * v))
    if v > 31.4978
        dss = 1
    end
    td = 0.6 + 1.0 / (exp(-0.05 * (v + 6.0)) + exp(0.09 * (v + 14.0)))
    dd = (dss - d) / td
    fss = 1.0 / (1.0 + exp((v + 19.58) / 3.696))
    tff = 7.0 + 1.0 / (0.0045 * exp(-(v + 20.0) / 10.0) + 0.0045 * exp((v + 20.0) / 10.0))
    tfs = 1000.0 + 1.0 / (0.000035 * exp(-(v + 5.0) / 4.0) + 0.000035 * exp((v + 5.0) / 6.0))
    Aff = 0.6
    Afs = 1.0 - Aff
    dff = (fss - ff) / tff
    dfs = (fss - fs) / tfs
    f = Aff * ff + Afs * fs
    fcass = fss
    tfcaf = 7.0 + 1.0 / (0.04 * exp(-(v - 4.0) / 7.0) + 0.04 * exp((v - 4.0) / 7.0))
    tfcas = 100.0 + 1.0 / (0.00012 * exp(-v / 3.0) + 0.00012 * exp(v / 7.0))
    Afcaf = 0.3 + 0.6 / (1.0 + exp((v - 10.0) / 10.0))
    Afcas = 1.0 - Afcaf
    dfcaf = (fcass - fcaf) / tfcaf
    dfcas = (fcass - fcas) / tfcas
    fca = Afcaf * fcaf + Afcas * fcas
    tjca = 72.5
    jcass = 1.0 / (1.0 + exp((v + 18.08) / (2.7916)))
    djca = (jcass - jca) / tjca
    tffp = 2.5 * tff
    dffp = (fss - ffp) / tffp
    fp = Aff * ffp + Afs * fs
    tfcafp = 2.5 * tfcaf    
    dfcafp = (fcass - fcafp) / tfcafp
    fcap = Afcaf * fcafp + Afcas * fcas
    Kmn = 0.002
    k2n = 500.0
    km2n = jca * 1
    anca = 1.0 / (k2n / km2n + (1.0 + Kmn / cass)^4.0)
    dnca = anca * k2n - nca * km2n
    anca_i = 1.0 / (k2n / km2n + (1.0 + Kmn / cai)^4.0)
    dnca_i = anca_i * k2n - nca_i * km2n
    Io = 0.5 * (nao + ko + clo + 4 * cao) / 1000
    Ii = 0.5 * (nass + kss + clss + 4 * cass) / 1000
    dielConstant = 74
    temp = 310
    constA = 1.82 * 10^6 * (dielConstant * temp)^(-1.5)
    gamma_cai = 10^(-constA * 4 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_cao = 10^(-constA * 4 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    gamma_nai = 10^(-constA * 1 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_nao = 10^(-constA * 1 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    gamma_ki = 10^(-constA * 1 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_kao = 10^(-constA * 1 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    PhiCaL_ss = 4.0 * vffrt * (gamma_cai * cass * exp(2.0 * vfrt) - gamma_cao * cao) / (exp(2.0 * vfrt) - 1.0)
    PhiCaNa_ss = 1.0 * vffrt * (gamma_nai * nass * exp(1.0 * vfrt) - gamma_nao * nao) / (exp(1.0 * vfrt) - 1.0)
    PhiCaK_ss = 1.0 * vffrt * (gamma_ki * kss * exp(1.0 * vfrt) - gamma_kao * ko) / (exp(1.0 * vfrt) - 1.0)
    Ii = 0.5 * (nai + ki + cli + 4 * cai) / 1000
    gamma_cai = 10^(-constA * 4 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_cao = 10^(-constA * 4 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    gamma_nai = 10^(-constA * 1 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_nao = 10^(-constA * 1 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    gamma_ki = 10^(-constA * 1 * (sqrt(Ii) / (1 + sqrt(Ii)) - 0.3 * Ii))
    gamma_kao = 10^(-constA * 1 * (sqrt(Io) / (1 + sqrt(Io)) - 0.3 * Io))
    gammaCaoMyo = gamma_cao
    gammaCaiMyo = gamma_cai
    PhiCaL_i = 4.0 * vffrt * (gamma_cai * cai * exp(2.0 * vfrt) - gamma_cao * cao) / (exp(2.0 * vfrt) - 1.0)
    PhiCaNa_i = 1.0 * vffrt * (gamma_nai * nai * exp(1.0 * vfrt) - gamma_nao * nao) / (exp(1.0 * vfrt) - 1.0)
    PhiCaK_i = 1.0 * vffrt * (gamma_ki * ki * exp(1.0 * vfrt) - gamma_kao * ko) / (exp(1.0 * vfrt) - 1.0)
    PCap = 1.1 * PCa
    PCaNa = 0.00125 * PCa * 1.1737 / 1.8969
    PCaK = 3.574e-4 * PCa * 1.1737 / 1.8969
    PCaNap = 0.00125 * PCap
    PCaKp = 3.574e-4 * PCap
    ICaL_ss = (1.0 - fICaLp) * PCa * PhiCaL_ss * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCap * PhiCaL_ss * d * (fp * (1.0 - nca) + jca * fcap * nca)
    ICaNa_ss = (1.0 - fICaLp) * PCaNa * PhiCaNa_ss * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaNap * PhiCaNa_ss * d * (fp * (1.0 - nca) + jca * fcap * nca)
    ICaK_ss = (1.0 - fICaLp) * PCaK * PhiCaK_ss * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaKp * PhiCaK_ss * d * (fp * (1.0 - nca) + jca * fcap * nca)
    ICaL_i = (1.0 - fICaLp) * PCa * PhiCaL_i * d * (f * (1.0 - nca_i) + jca * fca * nca_i) + fICaLp * PCap * PhiCaL_i * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i)
    ICaNa_i = (1.0 - fICaLp) * PCaNa * PhiCaNa_i * d * (f * (1.0 - nca_i) + jca * fca * nca_i) + fICaLp * PCaNap * PhiCaNa_i * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i)
    ICaK_i = (1.0 - fICaLp) * PCaK * PhiCaK_i * d * (f * (1.0 - nca_i) + jca * fca * nca_i) + fICaLp * PCaKp * PhiCaK_i * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i)
    ICaL_i = ICaL_i * (1 - ICaL_fracSS)
    ICaNa_i = ICaNa_i * (1 - ICaL_fracSS)
    ICaK_i = ICaK_i * (1 - ICaL_fracSS)
    ICaL_ss = ICaL_ss * ICaL_fracSS
    ICaNa_ss = ICaNa_ss * ICaL_fracSS
    ICaK_ss = ICaK_ss * ICaL_fracSS
    return ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, gammaCaoMyo, gammaCaiMyo
end

# IKr
function getIKr_ORd2011_MM(V, c0, c1, c2, o, i, ko, EK, vfrt, GKr)
    alpha = 0.1161 * exp(0.2990 * vfrt)
    beta =  0.2442 * exp(-1.604 * vfrt)
    alpha1 = 1.25 * 0.1235
    beta1 =  0.1911
    alpha2 = 0.0578 * exp(0.9710 * vfrt)
    beta2 = 0.349e-3 * exp(-1.062 * vfrt)
    alphai = 0.2533 * exp(0.5953 * vfrt)
    betai = 1.25 * 0.0522 * exp(-0.8209 * vfrt)
    alphac2ToI = 0.52e-4 * exp(1.525 * vfrt)
    betaiToC2 = (beta2 * betai * alphac2ToI) / (alpha2 * alphai)
    dc0 = c1 * beta - c0 * alpha
    dc1 = c0 * alpha + c2 * beta1 - c1 * (beta + alpha1)
    dc2 = c1 * alpha1 + o * beta2 + i * betaiToC2 - c2 * (beta1 + alpha2 + alphac2ToI)
    doo = c2 * alpha2 + i * betai - o * (beta2 + alphai)
    di = c2 * alphac2ToI + o * alphai - i * (betaiToC2 + betai)
    IKr = GKr * o * (V - EK)
    return IKr, dc0, dc1, dc2, doo, di
end

# IKs
function getIKs_ORd2011(v, xs1, xs2, cai, EKs, GKs)
    xs1ss = 1.0 / (1.0 + exp((-(v + 11.60)) / 8.932))
    txs1 = 817.3 + 1.0 / (2.326e-4 * exp((v + 48.28) / 17.80) + 0.001292 * exp((-(v + 210.0)) / 230.0))
    dxs1 = (xs1ss - xs1) / txs1
    xs2ss = xs1ss
    txs2 = 1.0 / (0.01 * exp((v - 50.0) / 20.0) + 0.0193 * exp((-(v + 66.54)) / 31.0))
    dxs2 = (xs2ss - xs2) / txs2
    KsCa = 1.0 + 0.6 / (1.0 + (3.8e-5 / max(cai, 1e-6)) ^ 1.4)  # max here helps mitigate the numerical instability when cai becomes negative for extreme mechanisms' perturbations, however, it is important to address the underlying issues causing negative Ca2+ concentrations to ensure the accuracy and reliability of simulation results
    IKs = GKs * KsCa * xs1 * xs2 * (v - EKs)
    return IKs, dxs1, dxs2
end

# IK1
function getIK1_CRLP(v, ko, EK, GK1)
    aK1 = 4.094 / (1 + exp(0.1217 * (v - EK - 49.934)))
    bK1 = (15.72 * exp(0.0674 * (v - EK - 3.257)) + exp(0.0618 * (v - EK - 594.31))) / (1 + exp(-0.1629 * (v - EK + 14.207)))
    K1ss = aK1 / (aK1 + bK1)
    IK1 = GK1 * sqrt(ko / 5) * K1ss * (v - EK)
    return IK1
end

# INaCa
function getINaCa_ORd2011(v, F, R, T, nass, nai, nao, cass, cai, cao, INaCa_fracSS, Gncx)
    zca = 2.0
    kna1 = 15.0
    kna2 = 5.0
    kna3 = 88.12
    kasymm = 12.5
    wna = 6.0e4
    wca = 6.0e4
    wnaca = 5.0e3
    kcaon = 1.5e6
    kcaoff = 5.0e3
    qna = 0.5224
    qca = 0.1670
    hca = exp((qca * v * F)/(R * T))
    hna = exp((qna * v * F)/(R * T))
    h1 = 1 + nai/kna3 * (1 + hna)
    h2 = (nai * hna)/(kna3 * h1)
    h3 = 1.0/h1
    h4 = 1.0 + nai/kna1 * (1 + nai/kna2)
    h5 = nai * nai/(h4 * kna1 * kna2)
    h6 = 1.0/h4
    h7 = 1.0 + nao/kna3 * (1.0 + 1.0/hna)
    h8 = nao/(kna3 * hna * h7)
    h9 = 1.0/h7
    h10 = kasymm + 1.0 + nao/kna1 * (1.0 + nao/kna2)
    h11 = nao * nao/(h10 * kna1 * kna2)
    h12 = 1.0/h10
    k1 = h12 * cao * kcaon
    k2 = kcaoff
    k3p = h9 * wca
    k3pp = h8 * wnaca
    k3 = k3p + k3pp
    k4p = h3 * wca/hca
    k4pp = h2 * wnaca
    k4 = k4p + k4pp
    k5 = kcaoff
    k6 = h6 * cai * kcaon
    k7 = h5 * h2 * wna
    k8 = h8 * h11 * wna
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
    E1 = x1/(x1 + x2 + x3 + x4)
    E2 = x2/(x1 + x2 + x3 + x4)
    E3 = x3/(x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    KmCaAct = 150.0e-6
    allo = 1.0 / (1.0 + (KmCaAct / cai) ^ 2.0)
    zna = 1.0
    JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
    JncxCa = E2 * k2 - E1 * k1
    INaCa_i = (1 - INaCa_fracSS) * Gncx * allo * (zna * JncxNa + zca * JncxCa)
    h1 = 1 + nass / kna3 * (1 + hna)
    h2 = (nass * hna) / (kna3 * h1)
    h3 = 1.0 / h1
    h4 = 1.0 + nass / kna1 * (1 + nass / kna2)
    h5 = nass * nass / (h4 * kna1 * kna2)
    h6 = 1.0 / h4
    h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna)
    h8 = nao / (kna3 * hna * h7)
    h9 = 1.0 / h7
    h10 = kasymm + 1.0 + nao / kna1 * (1 + nao / kna2)
    h11 = nao * nao / (h10 * kna1 * kna2)
    h12 = 1.0 / h10
    k1 = h12 * cao * kcaon
    k2 = kcaoff
    k3p = h9 * wca
    k3pp = h8 * wnaca
    k3 = k3p + k3pp
    k4p = h3 * wca / hca
    k4pp = h2 * wnaca
    k4 = k4p + k4pp
    k5 = kcaoff
    k6 = h6 * cass * kcaon
    k7 = h5 * h2 * wna
    k8 = h8 * h11 * wna
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    KmCaAct = 150.0e-6
    allo = 1.0 / (1.0 + (KmCaAct / cass)^2.0)
    JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
    JncxCa = E2 * k2 - E1 * k1
    INaCa_ss = INaCa_fracSS * Gncx * allo * (zna * JncxNa + zca * JncxCa)
    return INaCa_i, INaCa_ss
end

# INaK
function getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, Pnak)
    zna = 1.0
    k1p = 949.5
    k1m = 182.4
    k2p = 687.2
    k2m = 39.4
    k3p = 1899.0
    k3m = 79300.0
    k4p = 639.0
    k4m = 40.0
    Knai0 = 9.073
    Knao0 = 27.78
    delta = -0.1550
    Knai = Knai0 * exp((delta * v * F)/(3.0 * R * T))
    Knao = Knao0 * exp(((1.0 - delta) * v * F)/(3.0 * R * T))
    Kki = 0.5
    Kko = 0.3582
    MgADP = 0.05
    MgATP = 9.8
    Kmgatp = 1.698e-7
    H = 1.0e-7
    eP = 4.2
    Khp = 1.698e-7
    Knap = 224.0
    Kxkur = 292.0
    P = eP / (1.0 + H / Khp + nai / Knap + ki / Kxkur)
    a1 = (k1p * (nai / Knai)^3.0) / ((1.0 + nai / Knai)^3.0 + (1.0 + ki / Kki)^2.0 - 1.0)
    b1 = k1m * MgADP
    a2 = k2p
    b2 = (k2m * (nao / Knao)^3.0) / ((1.0 + nao / Knao)^3.0 + (1.0 + ko / Kko)^2.0 - 1.0)
    a3 = (k3p * (ko / Kko)^2.0) / ((1.0 + nao / Knao)^3.0 + (1.0 + ko / Kko)^2.0 - 1.0)
    b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp)
    a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp)
    b4 = (k4m * (ki / Kki)^2.0) / ((1.0 + nai / Knai)^3.0 + (1.0 + ki / Kki)^2.0 - 1.0)
    x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2
    x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4
    x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1
    x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    zk = 1.0
    JnakNa = 3.0 * (E1 * a3 - E2 * b3)
    JnakK = 2.0 * (E4 * b1 - E3 * a1)
    INaK = Pnak * (zna * JnakNa + zk * JnakK)
    return INaK
end

# Jrel
function getJrel_ORd2011(Jrelnp, Jrelp, ICaL, cajsr, fJrelp, JrelG)
    jsrMidpoint = 1.7
    bt = 4.75
    a_rel = 0.5 * bt
    Jrel_inf = a_rel * (-ICaL) / (1.0 + (jsrMidpoint / cajsr)^8.0)
    tau_rel = bt / (1.0 + 0.0123 / cajsr)
    if tau_rel < 0.001
        tau_rel = 0.001
    end
    dJrelnp = (Jrel_inf - Jrelnp) / tau_rel
    btp = 1.25 * bt
    a_relp = 0.5 * btp
    Jrel_infp = a_relp * (-ICaL) / (1.0 + (jsrMidpoint / cajsr)^8.0)
    tau_relp = btp / (1.0 + 0.0123 / cajsr)
    if tau_relp < 0.001
        tau_relp = 0.001
    end
    dJrelp = (Jrel_infp - Jrelp) / tau_relp
    Jrel = JrelG * ((1.0 - fJrelp) * Jrelnp + fJrelp * Jrelp)
    return Jrel, dJrelnp, dJrelp
end


end  # module
