module Initializer
export initConstants, initStates


function initConstants(scaling)
    cao           = 1.8 * scaling[14]
    nao           = 140.0
    ko            = 5.0
    clo           = 150.0
    R             = 8314.0
    T             = 310.0
    F             = 96485.0
    L             = 0.01
    rad           = 0.0011
    vcell         = 1000.0 * 3.14 * rad^2 * L
    Ageo          = 2.0 * 3.14 * rad^2 + 2.0 * 3.14 * rad * L
    Acap          = 2.0 * Ageo
    vmyo          = 0.68 * vcell
    vnsr          = 0.0552 * vcell
    vjsr          = 0.0048 * vcell
    vss           = 0.02 * vcell
    ICaL_fracSS   = 0.8
    INaCa_fracSS  = 0.35
    cmdnmax       = 0.05
    kmcmdn        = 0.00238
    trpnmax       = 0.07
    BSRmax        = 0.047
    KmBSR         = 0.00087
    BSLmax        = 1.124
    KmBSL         = 0.0087
    csqnmax       = 10.0
    kmcsqn        = 0.8
    PKNa          = 0.01833
    KmCaMK        = 0.15
    aCaMK         = 0.05
    bCaMK         = 0.00068
    CaMKo         = 0.05
    KmCaM         = 0.0015
    lambda        = 1.0
    lambda_rate   = 0.0
    perm50        = 0.35
    TRPN_n        = 2.0
    koff          = 0.1
    dr            = 0.25
    wfrac         = 0.5
    TOT_A         = 25.0
    ktm_unblock   = 0.021 
    beta_1        = -2.4
    beta_0        = 2.3
    gamma         = 0.0085
    gamma_wu      = 0.615
    phi           = 2.23  
    nperm         = 2.036 
    ca50          = 0.805 * scaling[12]
    Tref          = 120.0
    nu            = 7.0
    mu            = 3.0
    k_ws          = 0.004
    k_uw          = 0.026
    lambda_min    = 0.87
    lambda_max    = 1.2
    GKb           = 0.0189
    PNab          = 1.9239e-09
    PCab          = 5.9194e-08 * 1.8969
    GpCa          = 5e-04
    Fjunc         = 1.0
    Fsl           = 1.0 - Fjunc
    GClCa         = 0.2843
    GClB          = 1.98e-3
    KdClCa        = 0.1
    stim_amp      = -53.0
    stim_duration = 1.0
    stim_period   = 1000.0
    GNa           = 11.7802 * scaling[1]
    GNaL          = 0.0279 * scaling[2]
    Gto           = 0.16 * scaling[3]
    PCa           = 8.3757e-05 * 1.8969 * scaling[9]
    GKr           = 0.0321 * sqrt(ko/5) * scaling[4]
    GKs           = 0.0011 * scaling[5]
    GK1           = 0.6992 * scaling[6]
    Gncx          = 0.0034 * scaling[7]
    Pnak          = 15.4509 * scaling[8]
    JrelG         = 1.5378 * scaling[10]
    JupG          = 1.0 * scaling[11]
    MyoG          = 1.0 * scaling[13]
    Istim         = 0.0
    [cao, nao, ko, clo, R, T, F, L, rad, vcell, Ageo, Acap, vmyo, vnsr, vjsr, vss, ICaL_fracSS, INaCa_fracSS, cmdnmax, kmcmdn, trpnmax, BSRmax, KmBSR, BSLmax, KmBSL, csqnmax, kmcsqn, PKNa, KmCaMK, aCaMK, bCaMK, CaMKo, KmCaM, lambda, lambda_rate, perm50, TRPN_n, koff, dr, wfrac, TOT_A, ktm_unblock, beta_1, beta_0, gamma, gamma_wu, phi, nperm, ca50, Tref, nu, mu, k_ws, k_uw, lambda_min, lambda_max, GKb, PNab, PCab, GpCa, Fjunc, Fsl, GClCa, GClB, KdClCa, stim_amp, stim_duration, stim_period, GNa, GNaL, Gto, PCa, GKr, GKs, GK1, Gncx, Pnak, JrelG, JupG, MyoG, Istim   ]
end

function initStates()
    [-87, 7, 7, 145, 145, 0.0001, 0.0001, 1.2, 1.2, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 29, 29, 0]
end


end  # module
