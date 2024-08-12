module Model
export computeRates!

using ForceImport
@force using ..Currents


function computeRates!(du, u, p, t)
    # get states
    v         = u[1]
    nai       = u[2]
    nass      = u[3]
    ki        = u[4]
    kss       = u[5]
    cai       = u[6]
    cass      = u[7]
    cansr     = u[8]
    cajsr     = u[9]
    m         = u[10]
    hp        = u[11]
    h         = u[12]
    j         = u[13]
    jp        = u[14]
    mL        = u[15]
    hL        = u[16]
    hLp       = u[17]
    a         = u[18]
    iF        = u[19]
    iS        = u[20]
    ap        = u[21]
    iFp       = u[22]
    iSp       = u[23]
    d         = u[24]
    ff        = u[25]
    fs        = u[26]
    fcaf      = u[27]
    fcas      = u[28]
    jca       = u[29]
    nca       = u[30]
    nca_i     = u[31]
    ffp       = u[32]
    fcafp     = u[33]
    xs1       = u[34]
    xs2       = u[35]
    Jrel_np   = u[36]
    CaMKt     = u[37]
    ikr_c0    = u[38]
    ikr_c1    = u[39]
    ikr_c2    = u[40]
    ikr_o     = u[41]
    ikr_i     = u[42]
    Jrel_p    = u[43]
    XS        = max(0, u[44])
    XW        = max(0, u[45])
    Ca_TRPN   = max(0, u[46])
    TmBlocked = u[47]
    ZETAS     = u[48]
    ZETAW     = u[49]
    cli       = u[50]
    clss      = u[51]
  # qNet      = u[52]

    # get parameters
    cao           = p[1]
    nao           = p[2]
    ko            = p[3]
    clo           = p[4]
    R             = p[5]
    T             = p[6]
    F             = p[7]
  # L             = p[8]
  # rad           = p[9]
  # vcell         = p[10]
  # Ageo          = p[11]
    Acap          = p[12]
    vmyo          = p[13]
    vnsr          = p[14]
    vjsr          = p[15]
    vss           = p[16]
    ICaL_fracSS   = p[17]
    INaCa_fracSS  = p[18]
    cmdnmax       = p[19]
    kmcmdn        = p[20]
    trpnmax       = p[21]
    BSRmax        = p[22]
    KmBSR         = p[23]
    BSLmax        = p[24]
    KmBSL         = p[25]
    csqnmax       = p[26]
    kmcsqn        = p[27]
    PKNa          = p[28]
    KmCaMK        = p[29]
    aCaMK         = p[30]
    bCaMK         = p[31]
    CaMKo         = p[32]
    KmCaM         = p[33]
    lambda        = p[34]
    lambda_rate   = p[35]
    perm50        = p[36]
    TRPN_n        = p[37]
    koff          = p[38]
    dr            = p[39]
    wfrac         = p[40]
    TOT_A         = p[41]
    ktm_unblock   = p[42]
    beta_1        = p[43]
  # beta_0        = p[44]
    gamma         = p[45]
    gamma_wu      = p[46]
    phi           = p[47]
    nperm         = p[48]
    ca50          = p[49]
  # Tref          = p[50]
    nu            = p[51]
    mu            = p[52]
    k_ws          = p[53]
    k_uw          = p[54]
  # lambda_min    = p[55]
  # lambda_max    = p[56]
    GKb           = p[57]
    PNab          = p[58]
    PCab          = p[59]
    GpCa          = p[60]
    Fjunc         = p[61]
    Fsl           = p[62]
    GClCa         = p[63]
    GClB          = p[64]
    KdClCa        = p[65]
  # stim_amp      = p[66]
  # stim_duration = p[67]
  # stim_period   = p[68]
    GNa           = p[69]
    GNaL          = p[70]
    Gto           = p[71]
    PCa           = p[72]
    GKr           = p[73]
    GKs           = p[74]
    GK1           = p[75]
    Gncx          = p[76]
    Pnak          = p[77]
    JrelG         = p[78]
    JupG          = p[79]
    MyoG          = p[80]
    Istim         = p[81]

    # compute algebraics
    vffrt = v * F * F / (R * T)
    vfrt = v * F / (R * T)

    ENa = (R * T / F) * log(nao / nai)
    EK = (R * T / F) * log(ko / ki)
    EKs = (R * T / F) * log((ko + PKNa * nao) / (ki + PKNa * nai))
    ECl = (R * T / F) * log(cli / clo)
    EClss = (R * T / F) * log(clss / clo)

    CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass)
    CaMKa = CaMKb + CaMKt
    dCaMKt = aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt
    fINap = (1.0 / (1.0 + KmCaMK / CaMKa))
    fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa))
    fItop = (1.0 / (1.0 + KmCaMK / CaMKa))
    fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa))
    fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa))
    fJupp = (1.0 / (1.0 + KmCaMK / CaMKa))

    k_ws = k_ws * mu
    k_uw = k_uw * nu

    cdw = phi * k_uw * (1 - dr) * (1 - wfrac) / ((1 - dr) * wfrac)
    cds = phi * k_ws * (1 - dr) * wfrac / dr
    k_wu = k_uw * (1 / wfrac - 1) - k_ws
    k_su = k_ws * (1 / dr - 1) * wfrac
    A = (0.25 * TOT_A) / ((1 - dr) * wfrac + dr) * (dr / 0.25)

    # lambda0 = min(lambda_max, lambda)
    # Lfac = max(0, 1 + beta_0 * (lambda0 + min(lambda_min, lambda0) - (1 + lambda_min)))

    XU = (1 - TmBlocked) - XW - XS
    xb_ws = k_ws * XW
    xb_uw = k_uw * XU * MyoG
    xb_wu = k_wu * XW
    xb_su = k_su * XS

    gamma_rate = gamma * max((ZETAS .> 0) .* ZETAS, (ZETAS .< -1) .* (-ZETAS - 1))
    xb_su_gamma = gamma_rate .* XS
    gamma_rate_w = gamma_wu * abs.(ZETAW)
    xb_wu_gamma = gamma_rate_w .* XW

    dXS = xb_ws - xb_su - xb_su_gamma
    dXW = xb_uw - xb_wu - xb_ws - xb_wu_gamma

    ca50 = ca50 + beta_1 * min(0.2, lambda - 1)
    dCa_TRPN = koff * (((cai * 1000) / ca50) ^ TRPN_n * (1 - Ca_TRPN) - Ca_TRPN)

    XSSS = dr * 0.5
    XWSS = (1-dr) * wfrac * 0.5
    ktm_block = ktm_unblock * (perm50^nperm) * 0.5 / (0.5 - XSSS - XWSS)

    dTmBlocked = ktm_block * min(100, (Ca_TRPN^-(nperm/2))) * XU - ktm_unblock * (Ca_TRPN^(nperm/2)) * TmBlocked
    dZETAS = A * lambda_rate - cds * ZETAS
    dZETAW = A * lambda_rate - cdw * ZETAW

    # Ta = Lfac * (Tref/dr) * ((ZETAS+1)*XS + (ZETAW)*XW)

    INa, dm, dh, dhp, dj, djp = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, GNa)

    INaL, dmL, dhL, dhLp = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, GNaL)

    Ito, da, diF, diS, dap, diFp, diSp = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, Gto)

    (ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, gammaCaoMyo, gammaCaiMyo) = getICaL_ORd2011_jt(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, clss, ICaL_fracSS, vffrt, vfrt, PCa)
    ICaL = ICaL_ss + ICaL_i
    ICaNa = ICaNa_ss + ICaNa_i
    ICaK = ICaK_ss + ICaK_i
    # ICaL_tot = ICaL + ICaNa + ICaK

    IKr, dikr_c0, dikr_c1, dikr_c2, dikr_o, dikr_i = getIKr_ORd2011_MM(v, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i, ko, EK, vfrt, GKr)

    IKs, dxs1, dxs2 = getIKs_ORd2011(v, xs1, xs2, cai, EKs, GKs)

    IK1 = getIK1_CRLP(v, ko, EK, GK1)

    INaCa_i, INaCa_ss = getINaCa_ORd2011(v, F, R, T, nass, nai, nao, cass, cai, cao, INaCa_fracSS, Gncx)
    # INaCa_TOT = INaCa_i + INaCa_ss

    INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, Pnak)

    xkb = 1.0 / (1.0 + exp(-(v - 10.8968) / (23.9871)))
    IKb = GKb * xkb * (v - EK)
    INab = PNab * vffrt * (nai * exp(vfrt) - nao) / (exp(vfrt) - 1.0)
    ICab = PCab * 4.0 * vffrt * (gammaCaiMyo * cai * exp(2.0 * vfrt) - gammaCaoMyo * cao) / (exp(2.0 * vfrt) - 1.0)
    IpCa = GpCa * cai / (0.0005 + cai)

    I_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / cass) * (v - EClss)
    I_ClCa_sl = Fsl * GClCa / (1 + KdClCa / cai) * (v - ECl)
    I_ClCa = I_ClCa_junc + I_ClCa_sl
    I_Clbk = GClB * (v - ECl)

    Jrel, dJrel_np, dJrel_p = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss, cajsr, fJrelp, JrelG)

    Jleak = 0.0048825 * cansr / 15.0

    Jupnp = 0.005425 * cai / (cai + 0.00092)
    Jupp = 2.75 * 0.005425 * cai / (cai + 0.00092 - 0.00017)

    Jup = JupG * ((1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak)
    Jtr = (cansr - cajsr) / 60

    INet = INaL + ICaL + Ito + IKr + IKs + IK1

    JdiffNa = (nass - nai) / 2.0
    JdiffK = (kss - ki) / 2.0
    JdiffCl = (clss - cli) / 2.0
    Jdiff = (cass - cai) / 0.2

    Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / (kmcmdn + cai)^2.0)
    dcai = Bcai * (-(ICaL_i + IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo - dCa_TRPN * trpnmax)

    Bcass = 1.0 / (1.0 + BSRmax * KmBSR / (KmBSR + cass)^2.0 + BSLmax * KmBSL / (KmBSL + cass)^2.0)
    dcass = Bcass * (-(ICaL_ss - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff)

    dcansr = Jup - Jtr * vjsr / vnsr

    Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / (kmcsqn + cajsr)^2.0)
    dcajsr = Bcajsr * (Jtr - Jrel)

    # compute rates
    du[1]  = -(INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i + INaCa_ss + INaK + INab + IKb + IpCa + ICab + I_ClCa + I_Clbk + Istim)
    du[2]  = -(ICaNa_i + INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo
    du[3]  = -(ICaNa_ss + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa
    du[4]  = -(ICaK_i + Ito + IKr + IKs + IK1 + IKb + Istim - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo
    du[5]  = -(ICaK_ss) * Acap / (F * vss) - JdiffK
    du[6]  = dcai
    du[7]  = dcass
    du[8]  = dcansr
    du[9]  = dcajsr
    du[10] = dm
    du[11] = dhp
    du[12] = dh
    du[13] = dj
    du[14] = djp
    du[15] = dmL
    du[16] = dhL
    du[17] = dhLp
    du[18] = da
    du[19] = diF
    du[20] = diS
    du[21] = dap
    du[22] = diFp
    du[23] = diSp
    du[24] = dd
    du[25] = dff
    du[26] = dfs
    du[27] = dfcaf
    du[28] = dfcas
    du[29] = djca
    du[30] = dnca
    du[31] = dnca_i
    du[32] = dffp
    du[33] = dfcafp
    du[34] = dxs1
    du[35] = dxs2
    du[36] = dJrel_np
    du[37] = dCaMKt
    du[38] = dikr_c0
    du[39] = dikr_c1
    du[40] = dikr_c2
    du[41] = dikr_o
    du[42] = dikr_i
    du[43] = dJrel_p
    du[44] = dXS
    du[45] = dXW
    du[46] = dCa_TRPN
    du[47] = dTmBlocked
    du[48] = dZETAS
    du[49] = dZETAW
    du[50] = -(I_Clbk + I_ClCa_sl) * Acap / (-1 * F * vmyo) + JdiffCl * vss / vmyo
    du[51] = -I_ClCa_junc * Acap / (-1 * F * vss) - JdiffCl
    du[52] = INet * 1.0e-3
    nothing
end


end  # module
