#!/usr/bin/env python3
# =============================================================================
#  ps27_nuclear_geometry_scaling.py
# -----------------------------------------------------------------------------
#  Pelny test fizycznej intuicji uzytkownika:
#    "Gestsze upakowanie ladunku/masy w jadrze silniej moduluje pole substratu
#     Phi. Jadra ciezkich atomow moga rezonowac ze soba tworzac nieoczywiste
#     studnie — sprzezenie nielokalne."
#
#  Forma (1 + eta Z^2/137^2) z ps25 jest PRZYBLIZENIEM pierwszego-rzedu
#  (lokalny Darwin term dla pojedynczego jadra). Pelna fizyka moze wymagac:
#
#  H-Z   (baseline P7.5b):  X = Z_h^2                                     [1 par.]
#  H-A   (nukleon masa):     X = A_h^2                                    [1 par.]
#  H-B   (objetosc jadra):   X = A_h^(4/3)                                [1 par.]
#  H-C   (liczba-skalowane): X = z_h * Z_h^2                              [1 par.]
#  H-D   (pairwise geom.):   X(xi) = Z_h^2 * (1 + z_h * exp(-d_HH/xi))    [2 par.]
#  H-D2  (pairwise z A):     X(xi) = A_h^2 * (1 + z_h * exp(-d_HH/xi))    [2 par.]
#
#  W kazdym wariancie fitujemy:
#      dlog_new = dlog + 2 log10(1 + eta * X / X_ref)
#  na subsetach: 6p-heavy (N=7), 4f (N=4), heavy-union (N=11), CORE N=25.
#
#  Cel: zobaczyc ktory X najlepiej opisuje sygnal. Jesli H-D z rozsadnym
#  xi_Phi (3-20 A) wygrywa, to potwierdza teze uzytkownika o nielokalnym
#  sprzezeniu substratu przez ciezkie jadra.
#
#  Test unifikacji: czy pojedyncze (eta, xi) dziala na 4f i 6p lacznie?
#  Gdyby tak, to byloby to unified substrate coherence effect — jedna stala
#  TGP zamiast dwoch. Niemal napewno NIE (rozne znaki w ps25), ale test i tak
#  warto zrobic: pokaze jak daleko jest od unifikacji.
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import pearsonr, f as f_dist

ALPHA2 = 1.0 / 137.036 ** 2

# ==============================================================================
#  DANE: pelny N=29 z polami krystalograficznymi
#  (name, T_obs, T_pred, orb_class, Z_h, A_h, z_h_coord, d_HH_angstrom)
#  Z_h = 0 => brak ciezkiego (< 50) => X = 0 dla wszystkich H-*
#  z_h_coord = liczba nearest-neighbor heavy-heavy (approximate)
#  d_HH = odleglosc najblizszych ciezkich-ciezkich (A)
# ==============================================================================
MATERIALS = [
    # cuprates — Z_h = cation (La, Ba, Bi, Tl, Hg, Nd)
    ("La2CuO4",     38.00,  53.27, "d_cu",   57, 139,  4, 3.78),
    ("YBCO",        92.00,  75.59, "d_cu",   56, 137,  4, 3.82),
    ("BiSCCO2212",  85.00,  75.59, "d_cu",   83, 209,  4, 3.82),
    ("Tl2212",     108.00,  75.76, "d_cu",   81, 205,  4, 3.85),
    ("Hg1223",     138.00,  92.82, "d_cu",   80, 201,  4, 3.855),
    ("Tl2223",     125.00,  92.58, "d_cu",   81, 205,  4, 3.82),
    ("Nd2CuO4",     24.00,  53.86, "d_cu",   60, 144,  4, 3.945),
    ("Bi2201",      34.00,  53.41, "d_cu",   83, 209,  4, 3.81),
    # elementy lekkie / posrednie
    ("Al",           1.18,   2.96, "s",       0,   0,  0, 0.0),
    ("Pb",           7.20,   4.96, "s",      82, 207, 12, 3.50),
    ("Nb",           9.26,  13.83, "d_other", 0,   0,  0, 0.0),
    ("V",            5.30,  10.95, "d_other", 0,   0,  0, 0.0),
    ("Hg_elem",      4.15,   2.52, "s",      80, 201,  6, 2.99),
    ("MgB2",        39.00,  31.37, "sp",      0,   0,  0, 0.0),
    # Fe-SC (Z_h = La/Nd/Ba jesli obecne; Fe Z=26 < 50 wiec 0)
    ("FeSe_bulk",    8.00,   7.85, "d_fe",    0,   0,  0, 0.0),
    ("FeSe/STO",    65.00,  57.23, "d_fe",    0,   0,  0, 0.0),
    ("Ba122-Co",    22.00,  17.90, "d_fe",   56, 137,  4, 3.96),
    ("LaFeAsO",     26.00,  16.35, "d_fe",   57, 139,  4, 4.04),
    ("NdFeAsO-F",   55.00,  24.15, "d_fe",   60, 144,  4, 3.97),
    ("Nb3Sn",       18.30,  22.03, "d_other", 0,   0,  0, 0.0),
    ("NbTi",        10.00,  17.24, "d_other", 0,   0,  0, 0.0),
    ("H3S",        203.00,  76.05, "sp",      0,   0,  0, 0.0),
    # 4f/5f hydridy i elementy
    ("LaH10",      250.00, 368.08, "f",      57, 139, 12, 3.61),
    ("CeH9",       100.00, 143.13, "f",      58, 140,  8, 3.50),
    ("CeH10",      115.00, 187.48, "f",      58, 140, 12, 3.50),
    ("La_amb",       6.00,  14.34, "f",      57, 139, 12, 3.77),
    ("Y_amb",        1.30,  18.46, "d_other", 0,   0,  0, 0.0),  # Y=39 <50
    ("Th_amb",       1.38,  10.40, "f",      90, 232, 12, 3.59),
    ("Ce_5GPa",      1.70,   6.02, "f",      58, 140, 12, 3.46),
]

OUTLIER_NAMES = {"H3S", "Y_amb", "Th_amb", "Ce_5GPa"}

# subsety
def subset(names):
    return [m for m in MATERIALS if m[0] in names]

SUBSET_6P = {"BiSCCO2212", "Tl2212", "Hg1223", "Tl2223", "Bi2201", "Pb", "Hg_elem"}
SUBSET_4F = {"LaH10", "CeH9", "CeH10", "La_amb"}
SUBSET_4F_EXT = {"LaH10", "CeH9", "CeH10", "La_amb", "Th_amb", "Ce_5GPa"}
SUBSET_CORE = {m[0] for m in MATERIALS if m[0] not in OUTLIER_NAMES}

# ==============================================================================
#  Scalar formulas dla kazdej hipotezy
# ==============================================================================

def X_H_Z(Z, A, z, d):
    return float(Z) ** 2 if Z > 0 else 0.0

def X_H_A(Z, A, z, d):
    return float(A) ** 2 if A > 0 else 0.0

def X_H_B(Z, A, z, d):
    return float(A) ** (4.0/3.0) if A > 0 else 0.0

def X_H_C(Z, A, z, d):
    return float(z) * (float(Z) ** 2) if Z > 0 else 0.0

def X_H_D(Z, A, z, d, xi):
    """Pairwise geometric: local (1) + coherent coupling (z neighbors at d)."""
    if Z == 0:
        return 0.0
    return (float(Z) ** 2) * (1.0 + float(z) * np.exp(-d / max(xi, 1e-3)))

def X_H_D2(Z, A, z, d, xi):
    """Pairwise geometric z nukleon mass."""
    if A == 0:
        return 0.0
    return (float(A) ** 2) * (1.0 + float(z) * np.exp(-d / max(xi, 1e-3)))


# ==============================================================================
#  Fit routines
# ==============================================================================

def dlog_baseline(rows):
    return np.array([np.log10(r[2]) - np.log10(r[1]) for r in rows])


def fit_eta_static(rows, scalar_fn, xi_val=None):
    """Fit pojedynczy eta przy ustalonym X(Z,A,z,d) [lub X(Z,A,z,d,xi)]."""
    dlog = dlog_baseline(rows)
    if xi_val is None:
        Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7]) for r in rows])
    else:
        Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7], xi_val) for r in rows])
    # normalizacja do X_ref = max aby eta bylo ~O(1)
    X_ref = max(Xs.max(), 1.0)
    Xn = Xs / X_ref

    def loss(eta):
        factor = 1.0 + eta * Xn
        factor = np.where(factor > 1e-6, factor, 1e-6)
        return float(np.sum((dlog + 2.0 * np.log10(factor)) ** 2))

    res = minimize_scalar(loss, bounds=(-10.0, 10.0), method="bounded")
    return float(res.x), float(res.fun), float(X_ref)


#  Fizyczne granice xi:
#    - lower: 0.3 A (mniej niz twardy r_nucleon ~ fm, ale < 1 wiazaniowej)
#    - upper: 50 A (kilka komorek elementarnych; powyzej zdegenerowane z
#      count-scaling bo exp(-d/xi) -> 1 i X = Z^2 * (1 + z) ~ H-C)
XI_LOWER = 0.3
XI_UPPER = 50.0


def fit_eta_xi(rows, scalar_fn):
    """Fit eta AND xi jednoczesnie (dla H-D, H-D2). Bounded xi w [0.3, 50] A."""
    dlog = dlog_baseline(rows)

    def loss(params):
        eta, xi = params
        # bariera penalna dla xi poza [XI_LOWER, XI_UPPER]
        if xi < XI_LOWER or xi > XI_UPPER:
            return 1e6 + 100.0 * (
                max(XI_LOWER - xi, 0.0) ** 2 +
                max(xi - XI_UPPER, 0.0) ** 2
            )
        Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7], xi) for r in rows])
        X_ref = max(Xs.max(), 1.0)
        Xn = Xs / X_ref
        factor = 1.0 + eta * Xn
        factor = np.where(factor > 1e-6, factor, 1e-6)
        return float(np.sum((dlog + 2.0 * np.log10(factor)) ** 2))

    # grid search nad xi, potem dofit
    best = (1e6, 0.0, 5.0)
    for xi_init in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 35.0]:
        for eta_init in [-2.0, -1.0, -0.3, 0.0, 0.3, 1.0, 2.0]:
            res = minimize(loss, [eta_init, xi_init],
                           method="Nelder-Mead",
                           options={"xatol": 1e-5, "fatol": 1e-7,
                                    "maxiter": 400})
            # wymuszaj, ze koncowe xi jest w bounds (Nelder-Mead nie ma bounds,
            # ale nasza loss ma penalty; zrzucamy te rozwiazania ktore wyplynely)
            if res.fun < best[0] and XI_LOWER <= res.x[1] <= XI_UPPER:
                best = (float(res.fun), float(res.x[0]), float(res.x[1]))

    # X_ref dla raportu z najlepszym xi
    best_xi = best[2]
    Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7], best_xi) for r in rows])
    X_ref = max(Xs.max(), 1.0)
    # flaga: czy xi uciekalo do gornej granicy?
    xi_at_bound = (best_xi > 0.95 * XI_UPPER)
    return best[1], best[0], X_ref, best_xi, xi_at_bound


def compute_dlog_corrected(rows, scalar_fn, eta, X_ref, xi=None):
    dlog = dlog_baseline(rows)
    if xi is None:
        Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7]) for r in rows])
    else:
        Xs = np.array([scalar_fn(r[4], r[5], r[6], r[7], xi) for r in rows])
    Xn = Xs / X_ref
    factor = 1.0 + eta * Xn
    factor = np.where(factor > 1e-6, factor, 1e-6)
    return dlog + 2.0 * np.log10(factor)


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))


def f_test(rss_null, rss_alt, n, k_null, k_alt):
    df_num = k_alt - k_null
    df_den = n - k_alt
    if rss_alt <= 0 or df_den <= 0 or df_num <= 0:
        return float("nan"), float("nan")
    F = ((rss_null - rss_alt) / df_num) / (rss_alt / df_den)
    p = 1.0 - f_dist.cdf(F, df_num, df_den)
    return float(F), float(p)


# ==============================================================================
#  Analiza jednego subsetu — wszystkie hipotezy
# ==============================================================================

def analyze_subset(subset_names, label, full_print=True):
    rows = subset(subset_names)
    n = len(rows)
    dlog = dlog_baseline(rows)
    bias0 = float(dlog.mean())
    rss_intercept = float(np.sum((dlog - bias0) ** 2))
    rms_intercept = rms(dlog - bias0)
    rms_bare = rms(dlog)

    if full_print:
        print(f"\n{'='*78}")
        print(f"  SUBSET: {label}  (N={n})")
        print(f"{'='*78}")
        print(f"  baseline:  bias0={bias0:+.4f}  RMS(bez korekt)={rms_bare:.4f}  "
              f"RMS(z intercept)={rms_intercept:.4f}")

    results = {}

    # H-Z, H-A, H-B, H-C
    for label_h, fn, kpar in [
        ("H-Z  Z^2",        X_H_Z, 1),
        ("H-A  A^2",        X_H_A, 1),
        ("H-B  A^(4/3)",    X_H_B, 1),
        ("H-C  z*Z^2",      X_H_C, 1),
    ]:
        eta, rss, X_ref = fit_eta_static(rows, fn)
        rms_alt = float(np.sqrt(rss / n))
        F, p = f_test(rss_intercept, rss, n, 1, 2)
        results[label_h] = (eta, rms_alt, F, p, kpar)
        if full_print:
            print(f"  {label_h:<15s}  eta={eta:+.4f}  RMS={rms_alt:.4f}  "
                  f"redukcja={(1-rms_alt/rms_intercept)*100:+.1f}%  "
                  f"F={F:.2f} p={p:.3f}{'*' if p<0.05 else ''}")

    # H-D, H-D2 (2 parametry: eta + xi)
    for label_h, fn, kpar in [
        ("H-D  pair(Z)",  X_H_D,  2),
        ("H-D2 pair(A)",  X_H_D2, 2),
    ]:
        eta, rss, X_ref, xi, xi_bound = fit_eta_xi(rows, fn)
        rms_alt = float(np.sqrt(rss / n))
        F, p = f_test(rss_intercept, rss, n, 1, 3)
        results[label_h] = (eta, rms_alt, F, p, kpar, xi, xi_bound)
        bound_flag = "  [xi@upper-bound:degenerate->H-C]" if xi_bound else ""
        if full_print:
            print(f"  {label_h:<15s}  eta={eta:+.4f}  xi={xi:.2f} A  "
                  f"RMS={rms_alt:.4f}  "
                  f"redukcja={(1-rms_alt/rms_intercept)*100:+.1f}%  "
                  f"F={F:.2f} p={p:.3f}{'*' if p<0.05 else ''}{bound_flag}")

    return results, rms_intercept, rms_bare


# ==============================================================================
#  Test jednoparametrowy: ranking wariantow
# ==============================================================================

def rank_results(results, label):
    print(f"\n  Ranking ({label}):")
    items = sorted(results.items(), key=lambda kv: kv[1][1])  # by RMS
    print(f"  {'hipoteza':<15s}  {'RMS':>8s}  {'eta':>8s}  "
          f"{'xi[A]':>7s}  {'F':>6s}  {'p':>7s}")
    for k, v in items:
        if len(v) == 7:  # with xi + bound flag
            eta, rms_alt, F, p, _, xi, xi_bound = v
            xi_str = f"{xi:.2f}" + ("^" if xi_bound else "")
        else:
            eta, rms_alt, F, p, _ = v
            xi_str = "-"
        p_str = f"{p:.3f}" + ("*" if p < 0.05 else "")
        print(f"  {k:<15s}  {rms_alt:8.4f}  {eta:+8.4f}  {xi_str:>7s}  "
              f"{F:6.2f}  {p_str:>7s}")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("  ps27_nuclear_geometry_scaling.py")
    print("  Test intuicji: nielokalne sprzezenie jader ciezkich atomow")
    print("  przez pole substratu Phi (forma pairwise z dlugoscia koherencji).")
    print("=" * 78)

    # -------- SUBSET 6p-heavy (N=7, glowny test P7.5b)
    res_6p, rms_6p_null, _ = analyze_subset(SUBSET_6P, "6p-heavy (N=7)")
    rank_results(res_6p, "6p-heavy")

    # -------- SUBSET 4f core (N=4, P7.5a baseline)
    res_4f, rms_4f_null, _ = analyze_subset(SUBSET_4F, "4f core (N=4)")
    rank_results(res_4f, "4f core")

    # -------- SUBSET 4f-extended (N=6, z Th_amb, Ce_5GPa)
    res_4fe, rms_4fe_null, _ = analyze_subset(SUBSET_4F_EXT,
                                              "4f extended (N=6, +Th_amb, Ce_5GPa)")
    rank_results(res_4fe, "4f extended")

    # -------- SUBSET heavy-union (N=11)
    heavy_union = SUBSET_6P | SUBSET_4F_EXT
    res_union, rms_union_null, _ = analyze_subset(heavy_union,
                                          "heavy-union (N=11, 6p + 4f-ext)")
    rank_results(res_union, "heavy-union")

    # -------- SUBSET CORE N=25 (full model test)
    res_core, rms_core_null, rms_core_bare = analyze_subset(
        SUBSET_CORE, "CORE N=25 (pelny)")
    rank_results(res_core, "CORE N=25")

    # ------------------------------------------------------------------
    #  Test unifikacji: czy pojedyncze (eta, xi) dziala na 4f + 6p?
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  TEST UNIFIKACJI: czy jedna stala eta + jedna xi_Phi opisuje")
    print("  lacznie 4f + 6p heavy? Jesli TAK -> unified substrate coherence.")
    print("=" * 78)

    rows_union = subset(heavy_union)
    dlog_u = dlog_baseline(rows_union)
    rss_u_intercept = float(np.sum((dlog_u - dlog_u.mean()) ** 2))
    rms_u_intercept = rms(dlog_u - dlog_u.mean())

    # H-D z jedna para (eta, xi) na cala union
    eta_u, rss_u, _, xi_u, xi_u_bound = fit_eta_xi(rows_union, X_H_D)
    rms_u = float(np.sqrt(rss_u / len(rows_union)))
    F, p = f_test(rss_u_intercept, rss_u, len(rows_union), 1, 3)
    bound_note = " [xi@upper-bound => degeneracja z H-C]" if xi_u_bound else ""
    print(f"\n  Unified H-D (eta + xi) na N={len(rows_union)}:")
    print(f"    eta = {eta_u:+.4f}   xi = {xi_u:.2f} A{bound_note}")
    print(f"    RMS = {rms_u:.4f}  (intercept baseline = {rms_u_intercept:.4f})")
    print(f"    redukcja = {(1-rms_u/rms_u_intercept)*100:+.1f}%")
    print(f"    F(2,{len(rows_union)-3}) = {F:.3f}  p = {p:.4f}"
          f"{'*' if p<0.05 else ''}")

    # Porownanie per-material co sie dzieje
    dlog_u_corr = compute_dlog_corrected(rows_union, X_H_D, eta_u,
                                          X_ref=max(np.array([X_H_D(r[4],r[5],r[6],r[7],xi_u) for r in rows_union]).max(),1.0),
                                          xi=xi_u)
    print(f"\n  Per-material po unified H-D:")
    print(f"  {'material':<14s}  {'dlog_base':>10s}  {'dlog_corr':>10s}  {'class':<8s}")
    for row, d0, d1 in zip(rows_union, dlog_baseline(rows_union), dlog_u_corr):
        print(f"  {row[0]:<14s}  {d0:+10.4f}  {d1:+10.4f}  {row[3]:<8s}")

    # ------------------------------------------------------------------
    #  Porownanie: H-D (unified, 2 par) vs P7.5a+b (split, 2 par)
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Porownanie: H-D unified vs P7.5a + P7.5b (dwa subsety, 2 par.)")
    print("=" * 78)

    # zastosuj P7.5a (eta_f=-1.308) + P7.5b (eta_6p=+0.331) do union subsetu
    ETA_F, ETA_6P = -1.308, +0.331

    rows_union_list = subset(heavy_union)
    dlogs_split = []
    for r in rows_union_list:
        name, Tobs, Tpred, orb, Zh, Ah, zh, dHH = r
        d0 = np.log10(Tpred) - np.log10(Tobs)
        factor = 1.0
        if orb == "f":
            factor *= (1.0 + ETA_F * Zh**2 * ALPHA2) ** 2
        if (orb == "d_cu" and Zh >= 80) or name in {"Pb", "Hg_elem"}:
            factor *= (1.0 + ETA_6P * Zh**2 * ALPHA2) ** 2
        d_corr = d0 + np.log10(max(factor, 1e-6))
        dlogs_split.append(d_corr)

    rms_split = rms(np.array(dlogs_split))
    print(f"  P7.5a+b split (fixed eta_f=-1.308, eta_6p=+0.331):")
    print(f"    RMS = {rms_split:.4f}  (redukcja vs intercept "
          f"{(1-rms_split/rms_u_intercept)*100:+.1f}%)")
    print(f"  H-D unified (fitted eta, xi):")
    print(f"    RMS = {rms_u:.4f}  (redukcja vs intercept "
          f"{(1-rms_u/rms_u_intercept)*100:+.1f}%)")

    if rms_u < rms_split:
        print(f"  => H-D unified WYGRYWA ({rms_u:.4f} < {rms_split:.4f}) "
              f"pomimo 1 para. mniej niz split.")
    else:
        print(f"  => P7.5a+b split jest lepsze ({rms_split:.4f} < {rms_u:.4f}).")
        print(f"     Unifikacja przez pojedyncze (eta, xi) nie zamyka")
        print(f"     roznicy znakow miedzy 4f (screening) i 6p (enhancement).")

    # ------------------------------------------------------------------
    #  Podsumowanie: ktory wariant jest fizyczniew najbardziej uprawniony?
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Podsumowanie — fizyczna interpretacja wynikow")
    print("=" * 78)

    print("""
  Gluczowe pytania:
    1. Czy H-A (nukleon mass) daje lepszy fit niz H-Z (ladunek jadra)?
       -> Jesli TAK: sygnal pochodzi z gestosci baryonowej, nie Coulomb.
    2. Czy H-D z pairwise geometrycznym poprawia nad H-Z?
       -> Jesli TAK: nielokalne sprzezenie substratu przez ciezkie jadra
          jest realne, a xi_Phi to nowa stala TGP (dlugosc koherencji).
    3. Czy fitted xi_Phi jest fizycznie rozsadne (~ 3-20 A)?
       -> Jesli xi << a (lattice): lokalny efekt, nie pairwise.
       -> Jesli xi ~ c (cell): koherentne sprzezenie miedzywarstwowe.
       -> Jesli xi >> cell: efekt jest globalny, inny mechanizm.
""")

    # finalna tabela rank
    print("\n  Finalna tabela — najlepsze RMS per subset:")
    print(f"  {'subset':<20s}  {'N':>3s}  {'best H':<12s}  {'RMS':>7s}  "
          f"{'redukcja':>9s}  {'eta':>8s}  {'xi':>6s}")
    for sublabel, subn, rsn, intercept in [
        ("6p-heavy", 7, res_6p, rms_6p_null),
        ("4f core", 4, res_4f, rms_4f_null),
        ("4f ext.", 6, res_4fe, rms_4fe_null),
        ("heavy-union", 11, res_union, rms_union_null),
        ("CORE N=25", 25, res_core, rms_core_null),
    ]:
        best_h = min(rsn.items(), key=lambda kv: kv[1][1])
        name_h, v = best_h
        rms_b = v[1]
        eta_b = v[0]
        if len(v) == 7:
            xi_b = v[5]
            xi_bound_b = v[6]
            xi_str = f"{xi_b:.2f}" + ("^" if xi_bound_b else "")
        else:
            xi_str = "-"
        red = (1 - rms_b / intercept) * 100.0
        print(f"  {sublabel:<20s}  {subn:3d}  {name_h:<12s}  {rms_b:7.4f}  "
              f"{red:+8.1f}%  {eta_b:+8.3f}  {xi_str:>7s}")


if __name__ == "__main__":
    main()
