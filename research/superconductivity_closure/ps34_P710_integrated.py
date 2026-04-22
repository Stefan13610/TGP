#!/usr/bin/env python3
# =============================================================================
#  ps34_P710_integrated.py
# -----------------------------------------------------------------------------
#  Pelna integracja P7.10 (N(E_F)-factor) z kalibracja ps29c (N=45, P7.5a+b).
#
#  Roznica vs ps33:
#    - ps33: uproszczone parametry, baseline RMS ~1.5 (niewiarygodne)
#    - ps34: import ps29c MATERIALS + P7.5a+b corrections, baseline RMS ~0.21
#
#  Test: dodaj g(N_F/N_F_ref) factor do istniejacego T_pred_corr i fituj p
#  minimalizujac:
#     (A) RMS_SC na phonon-mediated SC materialach (N=36 po off-model filter)
#     (B) liczba FAIL w negatywnej walidacji (23 non-SC z ps32)
#
#  Oczekiwania:
#    - p > 0 (monotonicznie rosnace z N_F)
#    - RMS_SC w akceptowalnym zakresie (<0.30)
#    - FAIL w non-SC redukowane z 20 do < 10
#    - Kandydat: power-law p_fit ~ 1 (BCS-like) lub sigmoid
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import pearsonr

# === Import dokladnie to samo jak ps29c ===
K_B = 8.617333e-5
A_BOHR = 0.52917721067
A_STAR = 7.725 * A_BOHR
C_0 = 48.8222
SIGMA_A = 2.5856
ALPHA_P6B = 1.04
LAMBDA_0_P6B = 0.0962
OMEGA_0 = 15.0
BETA_P6D = 2.527
ALPHA_PB = 0.2887
K_DW = 3.498
LAMBDA_E_CUP = 0.0513
A_ZR_SQ = 0.181
A_MAP = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

ALPHA2 = 1.0 / 137.036 ** 2
ETA_F = -1.308       # P7.5a
ETA_6P = +0.331      # P7.5b
N_F_REF = 0.90       # Nb (referencja P7.10)


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)


def M_gauss(a_A):
    n = max(1, round(a_A / A_STAR))
    d = a_A - n * A_STAR
    return np.exp(-d ** 2 / SIGMA_A ** 2)


def B_mag(lam_sf):
    return 1.0 / (1.0 + BETA_P6D * lam_sf)


def B_PB(mu_eff):
    return np.exp(-ALPHA_PB * mu_eff ** 2)


def Tc_phonon(a_A, orb, z, omega, lam_sf, mu_eff=0.0):
    A = A_MAP[orb] if isinstance(orb, str) else float(orb) * A_MAP["d"]
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * A ** 2 * M_gauss(a_A)
            * (Lambda_eff * 1e-3) / K_B
            * B_mag(lam_sf) * B_PB(mu_eff))


def Tc_cuprate(a_A, n_layers, z_planar=8, mu_eff=0.0):
    return (K_DW * k_d(z_planar) * C_0 * A_ZR_SQ * M_gauss(a_A)
            * np.sqrt(n_layers)
            * (LAMBDA_E_CUP * 1e-3) / K_B
            * B_PB(mu_eff))


def apply_P75_correction(orb_cls, geom, name, T_raw):
    Zh, Ah, zh, dHH = geom
    factor = 1.0
    if orb_cls == "f":
        factor *= (1.0 + ETA_F * Zh ** 2 * ALPHA2)
    if (orb_cls == "d_cu" and Zh >= 80) or name in {"Pb", "Hg_elem"}:
        factor *= (1.0 + ETA_6P * Zh ** 2 * ALPHA2)
    return T_raw * max(factor, 1e-6) ** 2


def g_powerlaw(x, p):
    if x <= 0: return 0.0
    return min(1.0, x ** p)


def g_sigmoid(x, x0, n):
    if x <= 0: return 0.0
    return 1.0 / (1.0 + (x0 / x) ** n)


def apply_P710(T, N_F, p):
    """P7.10 czynnik na T_pred przez g(N_F/N_F_ref)."""
    x = N_F / N_F_REF
    return T * g_powerlaw(x, p)


# =============================================================================
#   MATERIALS z ps29c + nowa kolumna N_F (states/eV/atom, DFT literature)
# =============================================================================
# format: (name, T_obs, kind, params, orb_cls, geom, N_F)
MATERIALS = [
    # ---- cuprates (N=13) ----
    ("La2CuO4",    38.00, "cuprate", (3.780, 1), "d_cu", (57, 139, 4, 3.78),   0.50),
    ("YBCO",       92.00, "cuprate", (3.820, 2), "d_cu", (56, 137, 4, 3.82),   0.70),
    ("BiSCCO2212", 85.00, "cuprate", (3.820, 2), "d_cu", (83, 209, 4, 3.82),   0.60),
    ("Tl2212",    108.00, "cuprate", (3.850, 2), "d_cu", (81, 205, 4, 3.85),   0.65),
    ("Hg1223",    138.00, "cuprate", (3.855, 3), "d_cu", (80, 201, 4, 3.855),  0.75),
    ("Tl2223",    125.00, "cuprate", (3.820, 3), "d_cu", (81, 205, 4, 3.82),   0.70),
    ("Nd2CuO4",    24.00, "cuprate", (3.945, 1), "d_cu", (60, 144, 4, 3.945),  0.45),
    ("Bi2201",     34.00, "cuprate", (3.810, 1), "d_cu", (83, 209, 4, 3.81),   0.55),
    ("Hg1201",     94.00, "cuprate", (3.880, 1), "d_cu", (80, 201, 1, 3.88),   0.60),
    ("Hg1212",    127.00, "cuprate", (3.860, 2), "d_cu", (80, 201, 2, 3.86),   0.70),
    ("Hg1234",    125.00, "cuprate", (3.850, 4), "d_cu", (80, 201, 4, 3.85),   0.80),
    ("Tl2201",     90.00, "cuprate", (3.860, 1), "d_cu", (81, 205, 2, 3.86),   0.55),
    ("Bi2223",    110.00, "cuprate", (3.820, 3), "d_cu", (83, 209, 4, 3.82),   0.70),
    # ---- phonon SC (N=32) ----
    ("Al",          1.18, "phonon",  (4.046, "s",  12, 15.0, 0.0), "s",       (0,0,0,0), 0.42),
    ("Pb",          7.20, "phonon",  (4.950, "sp", 12,  8.3, 0.0), "s",       (82,207,12,3.50), 0.36),
    ("Nb",          9.26, "phonon",  (3.301, "d",   8, 22.0, 0.20),"d_other", (0,0,0,0), 0.97),
    ("V",           5.30, "phonon",  (3.024, "d",   8, 31.0, 0.60),"d_other", (0,0,0,0), 1.32),
    ("Hg_elem",     4.15, "phonon",  (2.989, "sp",  6,  9.0, 0.0), "s",       (80,201,6,2.99), 0.27),
    ("MgB2",       39.00, "phonon",  (3.086, "sp",  5, 75.0, 0.0), "sp",      (0,0,0,0), 0.72),
    ("FeSe_bulk",   8.00, "phonon",  (3.770, "d",   8, 25.0, 0.90),"d_fe",    (0,0,0,0), 1.10),
    ("FeSe/STO",   65.00, "phonon",  (3.770, "d",   8, 80.0, 0.20),"d_fe",    (0,0,0,0), 1.20),
    ("Ba122-Co",   22.00, "phonon",  (3.960, "d",   8, 30.0, 0.30),"d_fe",    (56,137,4,3.96), 1.00),
    ("LaFeAsO",    26.00, "phonon",  (4.040, "d",   8, 35.0, 0.50),"d_fe",    (57,139,4,4.04), 1.05),
    ("NdFeAsO-F",  55.00, "phonon",  (3.970, "d",   8, 40.0, 0.30),"d_fe",    (60,144,4,3.97), 1.15),
    ("Nb3Sn",      18.30, "phonon",  (5.290, "d",  12, 22.0, 0.10),"d_other", (0,0,0,0), 0.85),
    ("NbTi",       10.00, "phonon",  (3.300, "d",   8, 25.0, 0.15),"d_other", (0,0,0,0), 0.88),
    ("H3S",       203.00, "phonon",  (3.100, "sp",  8, 175.0, 0.0),"sp",      (0,0,0,0), 0.42),
    ("LaH10",     250.00, "phonon",  (5.100, 1.0,  12, 250.0, 0.0),"f",       (57,139,12,3.61), 0.62),
    ("CeH9",      100.00, "phonon",  (3.500, 1.0,   8, 135.0, 0.0),"f",       (58,140,8,3.50), 0.50),
    ("CeH10",     115.00, "phonon",  (3.500, 1.0,   8, 175.0, 0.0),"f",       (58,140,12,3.50), 0.55),
    ("La_amb",      6.00, "phonon",  (3.770, 1.0,  12,  12.0, 0.10),"f",      (57,139,12,3.77), 0.45),
    ("Y_amb",       1.30, "phonon",  (3.648, 1.0,  12,  17.0, 0.15),"d_other",(0,0,0,0), 0.78),
    ("Th_amb",      1.38, "phonon",  (5.080, 1.0,  12,  10.0, 0.10),"f",      (90,232,12,3.59), 0.40),
    ("Ce_5GPa",     1.70, "phonon",  (4.900, 0.8,  12,  12.0, 0.30),"f",      (58,140,12,3.46), 0.35),
    ("Ta_elem",     4.47, "phonon",  (3.300, "d",   8, 21.0, 0.0), "d_other", (73,181,8,3.30), 0.88),
    ("W_elem",      0.012,"phonon",  (3.160, "d",   8, 26.0, 0.0), "d_other", (74,184,8,3.16), 0.32),
    ("Re_elem",     1.70, "phonon",  (2.760, "d",  12, 30.0, 0.0), "d_other", (75,186,12,2.76), 0.38),
    ("Os_elem",     0.65, "phonon",  (2.730, "d",  12, 33.0, 0.0), "d_other", (76,190,12,2.73), 0.30),
    ("Ir_elem",     0.11, "phonon",  (2.710, "d",  12, 30.0, 0.0), "d_other", (77,192,12,2.71), 0.45),
    ("YH6",       224.00, "phonon",  (3.680, 1.0,  12, 180.0, 0.0),"d_other", (0,0,0,0), 0.57),
    ("ThH10",     161.00, "phonon",  (3.810, 1.0,  12, 170.0, 0.0),"f",       (90,232,12,3.81), 0.60),
    ("LiFeAs",     18.00, "phonon",  (3.780, "d",   8, 25.0, 0.30),"d_fe",    (0,0,0,0), 1.00),
    ("FeSeTe",     14.00, "phonon",  (3.800, "d",   8, 22.0, 0.40),"d_fe",    (52,128,4,3.80), 1.10),
    ("UTe2",        1.60, "phonon",  (4.160, "f",   8,  10.0, 0.30),"f",      (92,238,4,4.16), 0.50),
    ("CeCoIn5",     2.30, "phonon",  (4.630, "f",   8,  10.0, 0.40),"f",      (58,140,4,4.63), 1.50),
]


# Off-model (znane ze sa poza domena: 5d + heavy-f)
OFF_MODEL = {"Ta_elem","W_elem","Re_elem","Os_elem","Ir_elem","UTe2","CeCoIn5"}


# =============================================================================
#   NON_SC (z ps32 + N_F)
# =============================================================================
NON_SC = [
    ("Cu",       0.00,    3.615, "d",  12, 28.0, 0.0,  0.14),
    ("Ag",       0.00035, 4.086, "d",  12, 21.5, 0.0,  0.10),
    ("Au",       0.00005, 4.078, "d",  12, 17.0, 0.0,  0.07),
    ("Pt",       0.00,    3.924, "d",  12, 23.9, 0.0,  2.00),
    ("Pd",       0.00,    3.890, "d",  12, 25.2, 0.0,  2.30),
    ("Li",       0.0004,  3.510, "s",   8, 30.0, 0.0,  0.44),
    ("Na",       0.00,    4.290, "s",   8, 13.0, 0.0,  0.28),
    ("K",        0.00,    5.330, "s",   8,  8.6, 0.0,  0.42),
    ("Fe",       0.00,    2.870, "d",   8, 35.0, 3.0,  1.50),
    ("Co",       0.00,    2.510, "d",  12, 40.0, 3.0,  1.10),
    ("Ni",       0.00,    3.520, "d",  12, 38.0, 2.0,  2.90),
    ("Ca_amb",   0.00,    5.588, "s",  12, 19.0, 0.0,  0.14),
    ("Sr_amb",   0.00,    6.080, "s",  12, 12.0, 0.0,  0.15),
    ("Ba_amb",   0.00,    5.030, "s",  12,  7.0, 0.0,  0.20),
    ("Mg",       0.00,    3.210, "s",  12, 40.0, 0.0,  0.11),
    ("Zn",       0.85,    2.660, "sp", 12, 26.0, 0.0,  0.13),
    ("Cd",       0.52,    2.980, "sp", 12, 16.0, 0.0,  0.10),
    ("Diamond",  0.00,    3.570, "sp",  4, 165.0,0.0,  0.00),
    ("Si_amb",   0.00,    5.431, "sp",  4,  55.0,0.0,  0.00),
    ("Ge_amb",   0.00,    5.658, "sp",  4,  35.0,0.0,  0.00),
    ("NaCl",     0.00,    5.640, "sp",  6, 26.0, 0.0,  0.00),
    ("LiF",      0.00,    4.030, "sp",  6, 80.0, 0.0,  0.00),
    ("Graphite", 0.00,    2.460, "sp",  3, 170.0,0.0,  0.03),
]


# =============================================================================
def compute_SC_tpred(mat, p_P710=None):
    """T_pred po P7.5a+b (+ P7.10 jesli p_P710 podane)."""
    name, Tobs, kind, params, orb_cls, geom, N_F = mat
    if kind == "cuprate":
        a_A, n_layers = params
        T = Tc_cuprate(a_A, n_layers)
    else:
        a_A, orb, z, om, lam, *rest = list(params) + [0.0]
        mu = rest[0] if rest else 0.0
        T = Tc_phonon(a_A, orb, z, om, lam, mu_eff=mu)
    T = apply_P75_correction(orb_cls, geom, name, T)
    if p_P710 is not None:
        T = apply_P710(T, N_F, p_P710)
    return T


def compute_nonsc_tpred(row, p_P710=None):
    """T_pred non-SC po P7.5 (brak 6p/f heavy) + opt. P7.10."""
    name, Tub, a, orb, z, om, lam, N_F = row
    T = Tc_phonon(a, orb, z, om, lam)
    if p_P710 is not None:
        T = apply_P710(T, N_F, p_P710)
    return T


def sc_rms(p_P710, exclude=None):
    """RMS(dlog) na SC CORE (bez off-model i bez exclude)."""
    resid = []
    for m in MATERIALS:
        name = m[0]
        if name in OFF_MODEL: continue
        if exclude and name in exclude: continue
        Tp = compute_SC_tpred(m, p_P710)
        if Tp <= 0: return 1e9
        resid.append(np.log10(Tp) - np.log10(m[1]))
    return float(np.sqrt(np.mean(np.asarray(resid) ** 2))), len(resid)


def nonsc_fails(p_P710):
    """Liczba FAIL na nieSC."""
    fails, borders = 0, 0
    Tpreds = []
    for row in NON_SC:
        name, Tub = row[0], row[1]
        Tp = compute_nonsc_tpred(row, p_P710)
        Tpreds.append(Tp)
        if Tub < 1e-4:
            if Tp > 1.0: fails += 1
            elif Tp > 0.1: borders += 1
        else:
            r = Tp / Tub
            if r > 20: fails += 1
            elif r > 5: borders += 1
    return fails, borders, float(np.median(Tpreds)), float(np.percentile(Tpreds, 95))


def fit_p(lambda_fails=0.05, lambda_median=0.3):
    """Fit p maksymalizujacy: keep SC + kill non-SC."""
    def obj(p):
        if p <= 0.05: return 1e9
        rms_sc, _ = sc_rms(p)
        fails, _, med, _ = nonsc_fails(p)
        return rms_sc + lambda_median * med + lambda_fails * fails
    res = minimize_scalar(obj, bounds=(0.1, 5.0), method="bounded")
    return float(res.x), float(res.fun)


def main():
    print("=" * 78)
    print("  ps34_P710_integrated.py — pełna integracja P7.10 z ps29c")
    print("=" * 78)

    # ---- (1) Baseline (bez P7.10)
    print("\n  (1) Baseline: P7.5a+b, bez P7.10 (g=1 all materials):")
    rms0, N_core = sc_rms(None)
    fails0, borders0, med0, p95_0 = nonsc_fails(None)
    print(f"       SC CORE (N={N_core}, bez off-model): RMS(dlog) = {rms0:.4f}")
    print(f"       non-SC (N={len(NON_SC)}): median T = {med0:.3f}K, "
          f"P95 = {p95_0:.1f}K, FAIL={fails0}, BORDER={borders0}")

    # ---- (2) Fit p minimalizujac lacznie
    print("\n  (2) Fit pojedynczego p w g(x) = min(1, x^p):")
    p_fit, _ = fit_p(lambda_fails=0.05, lambda_median=0.3)
    rms_fit, _ = sc_rms(p_fit)
    fails_fit, borders_fit, med_fit, p95_fit = nonsc_fails(p_fit)
    print(f"       p = {p_fit:.3f}")
    print(f"       SC: RMS = {rms_fit:.4f}  ({rms_fit-rms0:+.4f} vs baseline)")
    print(f"       non-SC: median = {med_fit:.3f}K, P95 = {p95_fit:.1f}K, "
          f"FAIL={fails_fit} (z {fails0}), BORDER={borders_fit}")

    # ---- (3) Sweep po p i zobacz trade-off
    print("\n  (3) Sweep p (trade-off):")
    print(f"  {'p':<6s} {'RMS_SC':>8s} {'medNonSC':>10s} {'FAIL':>4s} {'BORDER':>6s}")
    for p in [0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        rms_p, _ = sc_rms(p)
        fails_p, borders_p, med_p, _ = nonsc_fails(p)
        flag = " <-- fit" if abs(p - p_fit) < 0.1 else ""
        print(f"  {p:<6.2f} {rms_p:8.4f} {med_p:10.3f} {fails_p:4d}"
              f" {borders_p:6d}{flag}")

    # ---- (4) Per-material dla najlepszego p
    print(f"\n  (4) SC materialy po P7.10 (p = {p_fit:.3f}):")
    print(f"  {'name':<12s} {'T_obs':>7s} {'T_0':>8s} {'T_P710':>9s} "
          f"{'d0':>7s} {'d710':>7s} {'N_F':>5s} {'g(x)':>6s}")
    worst = []
    for m in MATERIALS:
        name, Tobs, kind, params, orb_cls, geom, N_F = m
        if name in OFF_MODEL:
            T0 = compute_SC_tpred(m, None)
            T710 = compute_SC_tpred(m, p_fit)
            d0 = np.log10(T0) - np.log10(Tobs)
            d710 = np.log10(T710) - np.log10(Tobs)
            g = g_powerlaw(N_F / N_F_REF, p_fit)
            print(f"  {name:<12s} {Tobs:7.2f} {T0:8.2f} {T710:9.2f} "
                  f"{d0:+7.3f} {d710:+7.3f} {N_F:5.2f} {g:6.3f} [off]")
            continue
        T0 = compute_SC_tpred(m, None)
        T710 = compute_SC_tpred(m, p_fit)
        d0 = np.log10(T0) - np.log10(Tobs)
        d710 = np.log10(T710) - np.log10(Tobs)
        g = g_powerlaw(N_F / N_F_REF, p_fit)
        worst.append((abs(d710), name, d710))
        mark = " <-- out" if abs(d710) > 0.4 else ""
        print(f"  {name:<12s} {Tobs:7.2f} {T0:8.2f} {T710:9.2f} "
              f"{d0:+7.3f} {d710:+7.3f} {N_F:5.2f} {g:6.3f}{mark}")

    # ---- (5) non-SC po P7.10
    print(f"\n  (5) non-SC materialy po P7.10 (p = {p_fit:.3f}):")
    print(f"  {'name':<10s} {'T_ub':>8s} {'T_base':>8s} {'T_P710':>9s} "
          f"{'N_F':>5s} {'g':>6s} {'verdict':>8s}")
    for row in NON_SC:
        name, Tub = row[0], row[1]
        N_F = row[7]
        T0 = compute_nonsc_tpred(row, None)
        T710 = compute_nonsc_tpred(row, p_fit)
        g = g_powerlaw(N_F / N_F_REF, p_fit)
        if Tub < 1e-4:
            verdict = "PASS" if T710 < 0.1 else ("BORDER" if T710 < 1.0 else "FAIL")
        else:
            r = T710 / Tub
            verdict = ("PASS" if r < 5 else
                       ("BORDER" if r < 20 else "FAIL"))
        Tub_s = f"{Tub:.4f}" if Tub > 0 else "~0"
        T0_s = f"{T0:.2f}" if T0 < 100 else f"{T0:.1e}"
        T710_s = f"{T710:.3f}" if T710 < 100 else f"{T710:.1e}"
        print(f"  {name:<10s} {Tub_s:>8s} {T0_s:>8s} {T710_s:>9s} "
              f"{N_F:5.2f} {g:6.3f} {verdict:>8s}")

    # ---- (6) Per-klasa po P7.10
    print(f"\n  (6) Per-klasa RMS po P7.10 (p = {p_fit:.3f}):")
    classes = {}
    for m in MATERIALS:
        if m[0] in OFF_MODEL: continue
        orb_cls = m[4]
        T710 = compute_SC_tpred(m, p_fit)
        d = np.log10(T710) - np.log10(m[1])
        classes.setdefault(orb_cls, []).append((m[0], d))
    print(f"  {'klasa':<10s} {'N':>3s} {'mean':>8s} {'RMS':>8s}")
    for cls in sorted(classes):
        ds = [d for _, d in classes[cls]]
        mean_d = np.mean(ds)
        rms_d = np.sqrt(np.mean(np.asarray(ds) ** 2))
        print(f"  {cls:<10s} {len(ds):3d} {mean_d:+8.3f} {rms_d:8.3f}")

    # ---- (7) Correlation check
    Tobs_v = [m[1] for m in MATERIALS if m[0] not in OFF_MODEL]
    Tp_v = [compute_SC_tpred(m, p_fit) for m in MATERIALS if m[0] not in OFF_MODEL]
    r_P710, p_P710 = pearsonr(np.log10(Tobs_v), np.log10(Tp_v))

    Tp_v_0 = [compute_SC_tpred(m, None) for m in MATERIALS if m[0] not in OFF_MODEL]
    r_base, p_base = pearsonr(np.log10(Tobs_v), np.log10(Tp_v_0))

    print(f"\n  (7) Pearsonr (log-log) CORE N={len(Tobs_v)}:")
    print(f"       Baseline (P7.5 only):    r = {r_base:.4f}  p = {p_base:.3e}")
    print(f"       With P7.10 (p={p_fit:.2f}): r = {r_P710:.4f}  p = {p_P710:.3e}")

    # ---- (8) Verdict
    print("\n" + "=" * 78)
    print("  VERDICT P7.10 integrated")
    print("=" * 78)
    print(f"""
    Baseline (P7.5a+b):
      SC: RMS = {rms0:.4f}, r = {r_base:.4f}
      non-SC: FAIL = {fails0}/{len(NON_SC)} (87%), median = {med0:.2f} K

    + P7.10 (p = {p_fit:.3f}):
      SC: RMS = {rms_fit:.4f} ({'UTRZYMANE' if rms_fit < rms0+0.05 else 'POGORSZONE'})
      SC: r = {r_P710:.4f} ({'UTRZYMANE' if r_P710 > r_base-0.02 else 'POGORSZONE'})
      non-SC: FAIL = {fails_fit}/{len(NON_SC)} (-{fails0-fails_fit}), median = {med_fit:.3f} K

    Status: {'AKCEPTOWANE jako uniwersalna poprawka Eq. 5' if (rms_fit < rms0+0.05 and fails_fit < fails0 - 5) else 'WYMAGA re-fit globalny'}

    p = {p_fit:.3f} to nowa stala TGP:
      - Fizycznie: T_c ~ N_F^p ~ N_F (BCS-like dla p ~ 1)
      - Numerycznie: zamyka {fails0-fails_fit} spurious non-SC predictions
      - Koszt: {rms_fit-rms0:+.4f} dRMS na SC CORE
""")


if __name__ == "__main__":
    main()
