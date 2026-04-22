#!/usr/bin/env python3
# =============================================================================
#  ps29c_expanded_dataset.py
# -----------------------------------------------------------------------------
#  Rozszerzenie datasetu SC z N=29 do N=45 przez dodanie 16 starannie
#  wybranych materialow pokrywajacych luki:
#
#    (A) multi-layer Hg/Tl cuprates (z_h, d_HH niezaleznie zmieniaja sie)
#        -> Hg1201, Hg1212, Hg1234, Tl2201, Bi2223
#    (B) heavy d-metals pierwiastkowe (6p+ SOC test)
#        -> Ta, W, Re, Os, Ir
#    (C) nowe hydrydy (2019-2025): YH6, ThH10
#    (D) Fe-SC: LiFeAs, FeSeTe
#    (E) heavy fermion (f): UTe2, CeCoIn5
#
#  Cel: rozstrzygnac P7.5c (H-C vs H-D), poprawic statystyke η_f i η_6p,
#  zamknac N < 10 bariere F-testu.
#
#  Kanoniczna formula T_pred wg ps17 (tgp_sc paper, Eq. 5). Po T_pred
#  stosujemy P7.5a + P7.5b, a potem testujemy H-* na residuach.
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import pearsonr, f as f_dist

# ==============================================================================
#  TGP Tc CORE (konsystentne z ps17)
# ==============================================================================
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


def Tc_phonon(a_A, orb, z, omega, lam_sf, mu_eff=0.0, eta=1.0):
    """P6.B/C/D + P7 pathway (non-cuprates)."""
    if isinstance(orb, str):
        A = A_MAP[orb]
    else:
        A = float(orb) * A_MAP["d"]
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * A ** 2 * M_gauss(a_A)
            * (Lambda_eff * 1e-3) / K_B
            * B_mag(lam_sf) * B_PB(mu_eff))


def Tc_cuprate(a_A, n_layers, z_planar=8, mu_eff=0.0):
    """P6.A cuprate pathway."""
    return (K_DW * k_d(z_planar) * C_0 * A_ZR_SQ * M_gauss(a_A)
            * np.sqrt(n_layers)
            * (LAMBDA_E_CUP * 1e-3) / K_B
            * B_PB(mu_eff))


# ==============================================================================
#  DATASET: (name, T_obs, kind, params_for_Tpred, heavy_geom, orb_class)
#    kind = "phonon" or "cuprate"
#    params = (a_A, orb_or_eta, z, omega, lam_sf, mu_eff)  for phonon
#           = (a_A, n_layers)                              for cuprate
#    heavy_geom = (Z_h, A_h, z_h_coord, d_HH_angstrom)
#    orb_class = string used for P7.5a/b routing
# ==============================================================================

MATERIALS = [
    # -----------------  ISTNIEJACE N=29 (z ps17)  -----------------
    # cuprates
    ("La2CuO4",    38.00, "cuprate", (3.780, 1),                "d_cu", (57, 139,  4, 3.78)),
    ("YBCO",       92.00, "cuprate", (3.820, 2),                "d_cu", (56, 137,  4, 3.82)),
    ("BiSCCO2212", 85.00, "cuprate", (3.820, 2),                "d_cu", (83, 209,  4, 3.82)),
    ("Tl2212",    108.00, "cuprate", (3.850, 2),                "d_cu", (81, 205,  4, 3.85)),
    ("Hg1223",    138.00, "cuprate", (3.855, 3),                "d_cu", (80, 201,  4, 3.855)),
    ("Tl2223",    125.00, "cuprate", (3.820, 3),                "d_cu", (81, 205,  4, 3.82)),
    ("Nd2CuO4",    24.00, "cuprate", (3.945, 1),                "d_cu", (60, 144,  4, 3.945)),
    ("Bi2201",     34.00, "cuprate", (3.810, 1),                "d_cu", (83, 209,  4, 3.81)),
    # elementy
    ("Al",          1.18, "phonon",  (4.046, "s",  12, 15.0, 0.0), "s",       ( 0,   0,  0, 0.0)),
    ("Pb",          7.20, "phonon",  (4.950, "sp", 12,  8.3, 0.0), "s",       (82, 207, 12, 3.50)),
    ("Nb",          9.26, "phonon",  (3.301, "d",   8, 22.0, 0.20),"d_other", ( 0,   0,  0, 0.0)),
    ("V",           5.30, "phonon",  (3.024, "d",   8, 31.0, 0.60),"d_other", ( 0,   0,  0, 0.0)),
    ("Hg_elem",     4.15, "phonon",  (2.989, "sp",  6,  9.0, 0.0), "s",       (80, 201,  6, 2.99)),
    ("MgB2",       39.00, "phonon",  (3.086, "sp",  5, 75.0, 0.0), "sp",      ( 0,   0,  0, 0.0)),
    # Fe-SC
    ("FeSe_bulk",   8.00, "phonon",  (3.770, "d",   8, 25.0, 0.90),"d_fe",    ( 0,   0,  0, 0.0)),
    ("FeSe/STO",   65.00, "phonon",  (3.770, "d",   8, 80.0, 0.20),"d_fe",    ( 0,   0,  0, 0.0)),
    ("Ba122-Co",   22.00, "phonon",  (3.960, "d",   8, 30.0, 0.30),"d_fe",    (56, 137,  4, 3.96)),
    ("LaFeAsO",    26.00, "phonon",  (4.040, "d",   8, 35.0, 0.50),"d_fe",    (57, 139,  4, 4.04)),
    ("NdFeAsO-F",  55.00, "phonon",  (3.970, "d",   8, 40.0, 0.30),"d_fe",    (60, 144,  4, 3.97)),
    ("Nb3Sn",      18.30, "phonon",  (5.290, "d",  12, 22.0, 0.10),"d_other", ( 0,   0,  0, 0.0)),
    ("NbTi",       10.00, "phonon",  (3.300, "d",   8, 25.0, 0.15),"d_other", ( 0,   0,  0, 0.0)),
    ("H3S",       203.00, "phonon",  (3.100, "sp",  8, 175.0, 0.0),"sp",      ( 0,   0,  0, 0.0)),
    # 4f hydridy + elementy
    ("LaH10",     250.00, "phonon",  (5.100, 1.0,  12, 250.0, 0.0),"f",       (57, 139, 12, 3.61)),
    ("CeH9",      100.00, "phonon",  (3.500, 1.0,   8, 135.0, 0.0),"f",       (58, 140,  8, 3.50)),
    ("CeH10",     115.00, "phonon",  (3.500, 1.0,   8, 175.0, 0.0),"f",       (58, 140, 12, 3.50)),
    ("La_amb",      6.00, "phonon",  (3.770, 1.0,  12,  12.0, 0.10),"f",      (57, 139, 12, 3.77)),
    ("Y_amb",       1.30, "phonon",  (3.648, 1.0,  12,  17.0, 0.15),"d_other",( 0,   0,  0, 0.0)),
    ("Th_amb",      1.38, "phonon",  (5.080, 1.0,  12,  10.0, 0.10),"f",      (90, 232, 12, 3.59)),
    ("Ce_5GPa",     1.70, "phonon",  (4.900, 0.8,  12,  12.0, 0.30),"f",      (58, 140, 12, 3.46)),
    # -----------------  NOWE N=16 (ps29c)  -----------------
    # (A) multi-layer Hg/Tl/Bi cuprates (heavy 6p, rozne d_HH & n_layers)
    ("Hg1201",     94.00, "cuprate", (3.880, 1),                "d_cu", (80, 201,  1, 3.88)),
    ("Hg1212",    127.00, "cuprate", (3.860, 2),                "d_cu", (80, 201,  2, 3.86)),
    ("Hg1234",    125.00, "cuprate", (3.850, 4),                "d_cu", (80, 201,  4, 3.85)),
    ("Tl2201",     90.00, "cuprate", (3.860, 1),                "d_cu", (81, 205,  2, 3.86)),
    ("Bi2223",    110.00, "cuprate", (3.820, 3),                "d_cu", (83, 209,  4, 3.82)),
    # (B) heavy d-metals (6p-ish SOC test: Ta, W, Re, Os, Ir maja Z>=73)
    ("Ta_elem",     4.47, "phonon",  (3.300, "d",   8, 21.0, 0.0), "d_other", (73, 181,  8, 3.30)),
    ("W_elem",      0.012,"phonon",  (3.160, "d",   8, 26.0, 0.0), "d_other", (74, 184,  8, 3.16)),
    ("Re_elem",     1.70, "phonon",  (2.760, "d",  12, 30.0, 0.0), "d_other", (75, 186, 12, 2.76)),
    ("Os_elem",     0.65, "phonon",  (2.730, "d",  12, 33.0, 0.0), "d_other", (76, 190, 12, 2.73)),
    ("Ir_elem",     0.11, "phonon",  (2.710, "d",  12, 30.0, 0.0), "d_other", (77, 192, 12, 2.71)),
    # (C) nowe hydrydy 2019-2025
    ("YH6",       224.00, "phonon",  (3.680, 1.0,  12, 180.0, 0.0),"d_other", ( 0,   0,  0, 0.0)),   # Y=39 <50
    ("ThH10",     161.00, "phonon",  (3.810, 1.0,  12, 170.0, 0.0),"f",       (90, 232, 12, 3.81)),
    # (D) Fe-SC
    ("LiFeAs",     18.00, "phonon",  (3.780, "d",   8, 25.0, 0.30),"d_fe",    ( 0,   0,  0, 0.0)),
    ("FeSeTe",     14.00, "phonon",  (3.800, "d",   8, 22.0, 0.40),"d_fe",    (52, 128,  4, 3.80)),  # Te=52 threshold
    # (E) heavy fermion f-class
    ("UTe2",        1.60, "phonon",  (4.160, "f",   8,  10.0, 0.30),"f",      (92, 238,  4, 4.16)),
    ("CeCoIn5",     2.30, "phonon",  (4.630, "f",   8,  10.0, 0.40),"f",      (58, 140,  4, 4.63)),
]


OUTLIER_NAMES = {"H3S", "Y_amb"}  # nadal nie w pelni wyjasnione


# ==============================================================================
#  Oblicz T_pred i zastosuj P7.5a + P7.5b
# ==============================================================================

def compute_tpred(mat):
    name, Tobs, kind, params, orb_cls, geom = mat
    if kind == "cuprate":
        a_A, n_layers = params
        return Tc_cuprate(a_A, n_layers)
    else:
        a_A, orb, z, omega, lam, *rest = list(params) + [0.0]
        mu = rest[0] if rest else 0.0
        return Tc_phonon(a_A, orb, z, omega, lam, mu_eff=mu)


def apply_P75_correction(mat, T_pred_raw):
    """P7.5a dla f-klasy, P7.5b dla 6p-heavy cuprate + Pb/Hg_elem."""
    name, Tobs, kind, params, orb_cls, geom = mat
    Zh, Ah, zh, dHH = geom
    factor = 1.0
    if orb_cls == "f":
        factor *= (1.0 + ETA_F * Zh ** 2 * ALPHA2)
    if (orb_cls == "d_cu" and Zh >= 80) or name in {"Pb", "Hg_elem"}:
        factor *= (1.0 + ETA_6P * Zh ** 2 * ALPHA2)
    factor_sq = max(factor, 1e-6) ** 2   # T_c ~ A^2
    return T_pred_raw * factor_sq


# ==============================================================================
#  H-* scalars (jak ps27)
# ==============================================================================

def X_H_Z(Z, A, z, d):    return float(Z) ** 2 if Z > 0 else 0.0
def X_H_A(Z, A, z, d):    return float(A) ** 2 if A > 0 else 0.0
def X_H_B(Z, A, z, d):    return float(A) ** (4.0/3.0) if A > 0 else 0.0
def X_H_C(Z, A, z, d):    return float(z) * (float(Z) ** 2) if Z > 0 else 0.0

def X_H_D(Z, A, z, d, xi):
    if Z == 0: return 0.0
    return (float(Z) ** 2) * (1.0 + float(z) * np.exp(-d / max(xi, 1e-3)))

def X_H_D2(Z, A, z, d, xi):
    if A == 0: return 0.0
    return (float(A) ** 2) * (1.0 + float(z) * np.exp(-d / max(xi, 1e-3)))


XI_LOWER, XI_UPPER = 0.3, 50.0


def fit_eta_static(rows, scalar_fn):
    dlog = np.array([r["dlog"] for r in rows])
    Xs = np.array([scalar_fn(r["Zh"], r["Ah"], r["zh"], r["dHH"]) for r in rows])
    X_ref = max(Xs.max(), 1.0)
    Xn = Xs / X_ref

    def loss(eta):
        factor = 1.0 + eta * Xn
        factor = np.where(factor > 1e-6, factor, 1e-6)
        return float(np.sum((dlog + 2.0 * np.log10(factor)) ** 2))

    res = minimize_scalar(loss, bounds=(-10, 10), method="bounded")
    return float(res.x), float(res.fun), float(X_ref)


def fit_eta_xi(rows, scalar_fn):
    dlog = np.array([r["dlog"] for r in rows])

    def loss(p):
        eta, xi = p
        if xi < XI_LOWER or xi > XI_UPPER:
            return 1e6 + 100.0 * (max(XI_LOWER - xi, 0) ** 2
                                  + max(xi - XI_UPPER, 0) ** 2)
        Xs = np.array([scalar_fn(r["Zh"], r["Ah"], r["zh"], r["dHH"], xi)
                       for r in rows])
        X_ref = max(Xs.max(), 1.0)
        Xn = Xs / X_ref
        factor = 1.0 + eta * Xn
        factor = np.where(factor > 1e-6, factor, 1e-6)
        return float(np.sum((dlog + 2.0 * np.log10(factor)) ** 2))

    best = (1e6, 0.0, 5.0)
    for xi_init in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 35.0]:
        for eta_init in [-2.0, -1.0, -0.3, 0.0, 0.3, 1.0, 2.0]:
            res = minimize(loss, [eta_init, xi_init], method="Nelder-Mead",
                           options={"xatol": 1e-5, "fatol": 1e-7,
                                    "maxiter": 400})
            if res.fun < best[0] and XI_LOWER <= res.x[1] <= XI_UPPER:
                best = (float(res.fun), float(res.x[0]), float(res.x[1]))

    at_bound = best[2] > 0.95 * XI_UPPER
    return best[1], best[0], best[2], at_bound


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))


def f_test(rss0, rss1, n, k0, k1):
    df_n, df_d = k1 - k0, n - k1
    if rss1 <= 0 or df_d <= 0 or df_n <= 0:
        return float("nan"), float("nan")
    F = ((rss0 - rss1) / df_n) / (rss1 / df_d)
    return float(F), float(1.0 - f_dist.cdf(F, df_n, df_d))


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("  ps29c_expanded_dataset.py")
    print("  Rozszerzenie SC datasetu do N=45 (16 nowych materialow)")
    print("=" * 78)

    # Oblicz T_pred i zastosuj P7.5a+b
    enriched = []
    for m in MATERIALS:
        name, Tobs, kind, params, orb_cls, geom = m
        Zh, Ah, zh, dHH = geom
        Tp_raw = compute_tpred(m)
        Tp_corr = apply_P75_correction(m, Tp_raw)
        dlog_raw = np.log10(Tp_raw) - np.log10(Tobs)
        dlog_corr = np.log10(Tp_corr) - np.log10(Tobs)
        enriched.append({
            "name": name, "Tobs": Tobs, "Tp_raw": Tp_raw, "Tp_corr": Tp_corr,
            "dlog_raw": dlog_raw, "dlog_corr": dlog_corr,
            "Zh": Zh, "Ah": Ah, "zh": zh, "dHH": dHH,
            "orb_cls": orb_cls, "kind": kind,
            "new": name not in {
                "La2CuO4","YBCO","BiSCCO2212","Tl2212","Hg1223","Tl2223",
                "Nd2CuO4","Bi2201","Al","Pb","Nb","V","Hg_elem","MgB2",
                "FeSe_bulk","FeSe/STO","Ba122-Co","LaFeAsO","NdFeAsO-F",
                "Nb3Sn","NbTi","H3S","LaH10","CeH9","CeH10","La_amb",
                "Y_amb","Th_amb","Ce_5GPa"
            },
        })

    # -------- print full table
    print(f"\n  Kompletna tabela N={len(enriched)}:")
    print(f"  {'name':<12s} {'T_obs':>7s} {'T_pred0':>8s} {'T_pred(7.5)':>11s}"
          f" {'dlog_0':>7s} {'dlog_7.5':>8s} {'class':<8s} {'Zh':>3s}"
          f" {'new':>4s}")
    for r in enriched:
        marker = "*" if r["new"] else " "
        print(f"  {r['name']:<12s} {r['Tobs']:7.2f} {r['Tp_raw']:8.2f} "
              f"{r['Tp_corr']:11.2f} {r['dlog_raw']:+7.3f} "
              f"{r['dlog_corr']:+8.3f} {r['orb_cls']:<8s} "
              f"{r['Zh']:3d} {marker:>4s}")

    # -------- baseline stats
    N = len(enriched)
    Tobs_all = np.array([r["Tobs"] for r in enriched])
    Tpr_all = np.array([r["Tp_raw"] for r in enriched])
    Tpc_all = np.array([r["Tp_corr"] for r in enriched])
    r_raw, p_raw = pearsonr(np.log10(Tobs_all), np.log10(Tpr_all))
    r_corr, p_corr = pearsonr(np.log10(Tobs_all), np.log10(Tpc_all))
    rms_raw = rms([r["dlog_raw"] for r in enriched])
    rms_corr = rms([r["dlog_corr"] for r in enriched])

    print(f"\n  Baseline (pre-P7.5) na N={N}:")
    print(f"    r(log) = {r_raw:.4f}  (p = {p_raw:.3e})")
    print(f"    RMS_log = {rms_raw:.4f}")
    print(f"  Po P7.5a + P7.5b na N={N}:")
    print(f"    r(log) = {r_corr:.4f}  (p = {p_corr:.3e})")
    print(f"    RMS_log = {rms_corr:.4f}  (delta = {rms_corr-rms_raw:+.4f})")

    # -------- outliery
    outliers_raw = [r["name"] for r in enriched if abs(r["dlog_raw"]) > 0.4]
    outliers_corr = [r["name"] for r in enriched if abs(r["dlog_corr"]) > 0.4]
    print(f"\n  Outliers (|dlog|>0.4):")
    print(f"    pre-P7.5:  N={len(outliers_raw)}  {outliers_raw}")
    print(f"    post-P7.5: N={len(outliers_corr)}  {outliers_corr}")

    # -------- klasy
    classes = {}
    for r in enriched:
        classes.setdefault(r["orb_cls"], []).append(r["dlog_corr"])
    print(f"\n  Klasy po P7.5 (dlog):")
    print(f"  {'klasa':<10s} {'N':>3s} {'mean':>8s} {'RMS':>8s} {'max|d|':>8s}")
    for cls, dls in sorted(classes.items()):
        m = float(np.mean(dls))
        rm = rms(dls)
        mx = float(max(abs(d) for d in dls))
        print(f"  {cls:<10s} {len(dls):3d} {m:+8.3f} {rm:8.3f} {mx:8.3f}")

    # =====================================================================
    #  Rozszerzone P7.5c (H-C / H-D) — teraz z wiekszym N
    # =====================================================================
    print("\n" + "=" * 78)
    print("  P7.5c rerun: test H-Z/A/B/C/D na rozszerzonym datasecie")
    print("=" * 78)

    # === CHECKPOINT: 5d pierwiastkowe SC (Ta/W/Re/Os/Ir) sa OFF-MODEL
    # o 2-3 rzedow wielkosci. To ogromny sygnal, ale pochodzi z innej fizyki
    # (McMillan lambda ~ 1/(M omega^2) w phonon coupling dla ciezkich jader).
    # Heavy fermions (UTe2, CeCoIn5) tez sa off o ~100x (non-phonon pairing).
    # Dla sprawiedliwego testu P7.5c odfiltrowujemy je.
    OFF_MODEL_5D = {"Ta_elem", "W_elem", "Re_elem", "Os_elem", "Ir_elem"}
    OFF_MODEL_HF = {"UTe2", "CeCoIn5"}
    OFF_MODEL = OFF_MODEL_5D | OFF_MODEL_HF

    # Subsety
    # 6p-heavy CLEAN = cuprates Z>=80 + Pb + Hg_elem (oryginalna P7.5b domain)
    SUBSET_6P_CLEAN = [r for r in enriched
                       if (r["orb_cls"] == "d_cu" and r["Zh"] >= 80)
                       or r["name"] in {"Pb", "Hg_elem"}]
    # 6p/5d-heavy FULL = + 5d elementy (pokazuje ze rosnie szum)
    SUBSET_6P_5D = [r for r in enriched if r["orb_cls"] in {"d_cu"}
                    and r["Zh"] >= 73
                    or r["name"] in {"Pb", "Hg_elem",
                                     "Ta_elem","W_elem","Re_elem",
                                     "Os_elem","Ir_elem"}]
    # 4f/5f extended (class f, Zh>=50) — BEZ heavy fermions
    SUBSET_4F_EXT = [r for r in enriched if r["orb_cls"] == "f"
                     and r["Zh"] >= 50 and r["name"] not in OFF_MODEL_HF]
    # heavy-union CLEAN = 6p_clean + f_ext
    SUBSET_HU_CLEAN = list({id(x): x for x
                            in SUBSET_6P_CLEAN + SUBSET_4F_EXT}.values())
    # CORE: wszystko z Zh>0 BEZ off-model
    SUBSET_CORE = [r for r in enriched if r["Zh"] > 0
                   and r["name"] not in OFF_MODEL]

    # residuals uzywamy POST-P7.5 (zeby testowac geometrie jako ORTOGONALNE
    # wzgledem P7.5a+b)
    for r in enriched:
        r["dlog"] = r["dlog_corr"]

    def analyze(rows_raw, label):
        print(f"\n  SUBSET: {label}  (N={len(rows_raw)})")
        if len(rows_raw) < 3:
            print("    za maly subset — skipping")
            return
        dl = np.array([r["dlog"] for r in rows_raw])
        bias0 = float(dl.mean())
        rss_int = float(np.sum((dl - bias0) ** 2))
        rms_int = rms(dl - bias0)
        print(f"  bias0={bias0:+.4f}  RMS(intercept)={rms_int:.4f}")

        results = {}
        for lbl, fn in [("H-Z  Z^2", X_H_Z), ("H-A  A^2", X_H_A),
                        ("H-B  A^(4/3)", X_H_B), ("H-C  z*Z^2", X_H_C)]:
            eta, rss, _ = fit_eta_static(rows_raw, fn)
            F, p = f_test(rss_int, rss, len(rows_raw), 1, 2)
            rms_a = float(np.sqrt(rss / len(rows_raw)))
            results[lbl] = (rms_a, eta, F, p, None)
            mark = "*" if p < 0.05 else " "
            print(f"    {lbl:<15s} eta={eta:+.4f} RMS={rms_a:.4f} "
                  f"red={(1-rms_a/rms_int)*100:+.1f}% F={F:.2f} p={p:.3f}{mark}")

        for lbl, fn in [("H-D  pair(Z)", X_H_D), ("H-D2 pair(A)", X_H_D2)]:
            eta, rss, xi, at_bd = fit_eta_xi(rows_raw, fn)
            F, p = f_test(rss_int, rss, len(rows_raw), 1, 3)
            rms_a = float(np.sqrt(rss / len(rows_raw)))
            results[lbl] = (rms_a, eta, F, p, xi)
            mark = "*" if p < 0.05 else " "
            bnd = " [xi@bound]" if at_bd else ""
            print(f"    {lbl:<15s} eta={eta:+.4f} xi={xi:.2f}A RMS={rms_a:.4f} "
                  f"red={(1-rms_a/rms_int)*100:+.1f}% F={F:.2f} p={p:.3f}{mark}{bnd}")

        # ranking
        print(f"  Ranking:")
        items = sorted(results.items(), key=lambda kv: kv[1][0])
        for k, v in items:
            rms_a, eta, F, p, xi = v
            xi_s = f"{xi:.2f}" if xi is not None else "-"
            mark = "*" if p < 0.05 else " "
            print(f"    {k:<15s} RMS={rms_a:.4f}  eta={eta:+.4f}  "
                  f"xi={xi_s:>7s}  F={F:5.2f}  p={p:.3f}{mark}")

    analyze(SUBSET_6P_CLEAN, f"6p-heavy CLEAN (cuprates Z>=80 + Pb + Hg_elem)")
    analyze(SUBSET_6P_5D,    f"6p/5d-heavy FULL (incl. Ta..Ir off-model)")
    analyze(SUBSET_4F_EXT,   f"f-class extended (Z>=50, bez UTe2/CeCoIn5)")
    analyze(SUBSET_HU_CLEAN, f"heavy-union CLEAN")
    analyze(SUBSET_CORE,     f"CORE (wszystkie Zh>0 bez off-model)")

    # =====================================================================
    #  Interesujacy przypadek: ThH10
    # =====================================================================
    print("\n" + "=" * 78)
    print("  ThH10 analiza szczegolowa (hold-out dla P7.5a)")
    print("=" * 78)
    for r in enriched:
        if r["name"] == "ThH10":
            print(f"  T_obs = {r['Tobs']:.1f}")
            print(f"  T_pred(raw) = {r['Tp_raw']:.2f}   dlog_raw = {r['dlog_raw']:+.3f}")
            print(f"  T_pred(P7.5a) = {r['Tp_corr']:.2f}  dlog_corr = {r['dlog_corr']:+.3f}")
            print(f"  P7.5a z eta_f={ETA_F} i Z=90 to OVER-CORRECTION (-0.475)")
            print(f"  Sugestia: eta_f(Z) moze wymaga nasycenia przy bardzo ciezkich")
            print(f"  aktynowcach — Z^2/137^2 = 0.431 dla Z=90 jest w obszarze gdzie")
            print(f"  linearna aproksymacja QED sie zalamuje.")

    # =====================================================================
    #  Per-klasa residua po P7.5 — kierunkowskazy dla P7.6/P7.7
    # =====================================================================
    print("\n" + "=" * 78)
    print("  Residualne kanaly po P7.5 (identyfikacja dalszych subproblemow)")
    print("=" * 78)
    print("  Klasy z mean(dlog) istotnie !=0 wskazuja gdzie jest jeszcze fizyka")
    print("  niewlaczona w model.")

    for cls, dls in sorted(classes.items()):
        if len(dls) < 2: continue
        m = float(np.mean(dls))
        rm = rms(dls)
        se = rm / np.sqrt(len(dls))
        t_stat = m / se if se > 0 else 0.0
        sig = "" if abs(t_stat) < 1.0 else (" <-- SIGNAL" if abs(t_stat) < 2.0
                                             else " <-- STRONG SIGNAL")
        print(f"  {cls:<10s} N={len(dls):2d}  mean={m:+.3f}  RMS={rm:.3f}  "
              f"t={t_stat:+.2f}{sig}")


if __name__ == "__main__":
    main()
