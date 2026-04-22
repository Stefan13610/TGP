#!/usr/bin/env python3
# =============================================================================
#  ps25_eta_relativistic_correction.py
# -----------------------------------------------------------------------------
#  Test propozycji: wprowadzic do Eq. (5) [eq:A-orb] korekte
#      A_orb(Z) = A_orb * (1 + eta * Z^2 / 137^2)
#  z eta - uniwersalna stala sprzezenia TGP.
#
#  Poniewaz A_orb wchodzi kwadratowo do Eq. (3) [eq:ansatz],
#      T_c ~ A_orb^2 ~ (1 + eta * Z^2 / 137^2)^2
#  a zatem na skali log10:
#      dlog_new = dlog + 2 * log10(1 + eta * Z^2 / 137^2).
#
#  Zadane pytanie: czy 1 parametr eta (globalny) likwiduje trend w residuach?
#  Odpowiedz: ps25 sprawdza to formalnie (nieliniowy fit + F-test), a takze
#  testuje wersje per-klasowe (eta_4f vs eta_6p vs eta_d) ktore odpowiadaja
#  sub-problemom P7.5a i P7.5b.
#
#  Testy:
#    T1  global eta, Z = Z_pair  (literalne odczyt Eq. 5)
#    T2  global eta, Z = Z_heavy (empiryczne, gdzie jest korelacja z ps23)
#    T3  split eta per klasa orbital: (eta_s, eta_sp, eta_d, eta_f)
#    T4  P7.5a forma: eta_f tylko (A_f reduction dla 4f)
#    T5  P7.5b forma: eta_6p tylko (SOC correction dla 6p heavy spectators)
#    T6  P7.5a + P7.5b lacznie (2 parametry)
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import f as f_dist

# (material, Z_pair, Z_heavy, M_mol, dlog, orb_class)
DATA = [
    # orb_class in {"d_cu","d_fe","d_other","s","sp","f"}
    ("La2CuO4",     29, 57, 405.37,  +0.147, "d_cu"),
    ("YBCO",        29, 56, 664.62,  -0.085, "d_cu"),
    ("BiSCCO2212",  29, 83, 888.40,  -0.051, "d_cu"),
    ("Tl2212",      29, 81, 978.70,  -0.154, "d_cu"),
    ("Hg1223",      29, 80, 874.20,  -0.172, "d_cu"),
    ("Tl2223",      29, 81, 1114.40, -0.130, "d_cu"),
    ("Nd2CuO4",     29, 60, 416.03,  +0.351, "d_cu"),
    ("Bi2201",      29, 83, 752.75,  +0.196, "d_cu"),
    ("Al",          13, 13, 26.98,   +0.399, "s"),
    ("Pb",          82, 82, 207.20,  -0.162, "s"),
    ("Nb",          41, 41, 92.91,   +0.174, "d_other"),
    ("V",           23, 23, 50.94,   +0.315, "d_other"),
    ("Hg_elem",     80, 80, 200.59,  -0.217, "s"),
    ("MgB2",        5,  12, 45.93,   -0.094, "sp"),
    ("FeSe_bulk",   26, 34, 134.81,  -0.008, "d_fe"),
    ("FeSe/STO",    26, 34, 134.81,  -0.055, "d_fe"),
    ("Ba122-Co",    26, 56, 398.87,  -0.090, "d_fe"),
    ("LaFeAsO",     26, 57, 285.68,  -0.201, "d_fe"),
    ("NdFeAsO-F",   26, 60, 291.01,  -0.358, "d_fe"),
    ("Nb3Sn",       41, 50, 397.44,  +0.081, "d_other"),
    ("NbTi",        41, 41, 140.78,  +0.237, "d_other"),
    ("H3S",         16, 16, 35.09,   -0.426, "sp"),
    ("LaH10",       57, 57, 148.99,  +0.168, "f"),
    ("CeH9",        58, 58, 149.19,  +0.156, "f"),
    ("CeH10",       58, 58, 150.20,  +0.212, "f"),
    ("La_amb",      57, 57, 138.91,  +0.378, "f"),
    ("Y_amb",       39, 39, 88.91,   +1.152, "d_other"),
    ("Th_amb",      90, 90, 232.04,  +0.877, "f"),
    ("Ce_5GPa",     58, 58, 140.12,  +0.549, "f"),
]

OUTLIER_NAMES = {"H3S", "Y_amb", "Th_amb", "Ce_5GPa"}
ALPHA2 = 1.0 / 137.036 ** 2   # fine-structure squared
ALPHA_INV = 137.036


def dlog_with_eta(dlog, Z, eta):
    """Zastosuj korekte (1 + eta * Z^2 / 137^2) do A_orb (kwadratowo na T_c)."""
    factor = 1.0 + eta * (Z ** 2) * ALPHA2
    # jesli factor <= 0 (eta bardzo ujemne), logarytm zawala sie; zabezpiecz
    factor = np.where(factor > 1e-6, factor, 1e-6)
    return dlog + 2.0 * np.log10(factor)


def rms_log(dlog):
    return float(np.sqrt(np.mean(np.asarray(dlog, float) ** 2)))


def fit_single_eta(dlog, Z, bounds=(-5.0, 5.0)):
    def loss(eta):
        return float(np.sum(dlog_with_eta(dlog, Z, eta) ** 2))
    res = minimize_scalar(loss, bounds=bounds, method="bounded",
                          options={"xatol": 1e-6})
    return float(res.x), float(res.fun)


def fit_multi_eta(dlog, Z, class_idx, n_classes, x0=None):
    """Fit eta_k per klasa; Z[i], class_idx[i] -> eta = etas[class_idx[i]]."""
    dlog = np.asarray(dlog, float)
    Z = np.asarray(Z, float)
    class_idx = np.asarray(class_idx, int)

    def loss(etas):
        eta_arr = np.asarray(etas)[class_idx]
        factor = 1.0 + eta_arr * (Z ** 2) * ALPHA2
        factor = np.where(factor > 1e-6, factor, 1e-6)
        d_new = dlog + 2.0 * np.log10(factor)
        return float(np.sum(d_new ** 2))

    if x0 is None:
        x0 = np.zeros(n_classes)
    res = minimize(loss, x0, method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-8})
    return res.x.tolist(), float(res.fun)


def f_test(rss_null, rss_alt, n, k_null, k_alt):
    df_num = k_alt - k_null
    df_den = n - k_alt
    if rss_alt <= 0 or df_den <= 0 or df_num <= 0:
        return float("nan"), float("nan")
    F = ((rss_null - rss_alt) / df_num) / (rss_alt / df_den)
    p = 1.0 - f_dist.cdf(F, df_num, df_den)
    return F, p


def main():
    print("=" * 78)
    print("  ps25_eta_relativistic_correction.py")
    print("  Test propozycji: A_orb *= (1 + eta * Z^2 / 137^2)")
    print("=" * 78)

    # CORE N=25
    core = [r for r in DATA if r[0] not in OUTLIER_NAMES]
    names   = [r[0] for r in core]
    Z_pair  = np.array([r[1] for r in core], float)
    Z_heavy = np.array([r[2] for r in core], float)
    dlog    = np.array([r[4] for r in core], float)
    orbs    = [r[5] for r in core]
    n = len(core)

    # baseline null: tylko intercept (bias globalny)
    bias0 = float(dlog.mean())
    rss_null_intercept = float(np.sum((dlog - bias0) ** 2))
    rms_null = rms_log(dlog - bias0)
    rss_zero = float(np.sum(dlog ** 2))

    print(f"\nCORE N={n}")
    print(f"  Baseline RMS_log (wokol 0)           = {rms_log(dlog):.4f}")
    print(f"  Baseline RMS_log (z globalnym biasem={bias0:+.3f}) = {rms_null:.4f}")
    print(f"  1/137^2 = {ALPHA2:.6e}")

    # ------------------------------------------------------------------------
    #  T1: global eta, Z = Z_pair
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T1: global eta, Z = Z_pair  (dokladnie wg Eq. 5 paperu)")
    print("-" * 78)
    eta_T1, rss_T1 = fit_single_eta(dlog, Z_pair)
    rms_T1 = float(np.sqrt(rss_T1 / n))
    F, p = f_test(rss_null_intercept, rss_T1, n, k_null=1, k_alt=1)
    # uwaga: tutaj alt tez ma 1 param (eta), null tez 1 param (intercept). Nie
    # mozna robic F-testu w standardzie (nie zagniezdzone). Lepiej porownac
    # do modelu null=zero (0 parametrow).
    F0, p0 = f_test(rss_zero, rss_T1, n, k_null=0, k_alt=1)
    print(f"  best eta = {eta_T1:+.4f}  (Z^2/137^2 dla Z=29: {29**2*ALPHA2:.4f})")
    print(f"  RMS_new = {rms_T1:.4f}  (vs intercept baseline {rms_null:.4f})")
    print(f"  redukcja: {(1-rms_T1/rms_null)*100:+.1f}%")
    print(f"  F-test vs null=0: F(1,{n-1})={F0:.3f}  p={p0:.4f}{'*' if p0<0.05 else ''}")
    print(f"  WNIOSEK: Z=Z_pair=29 dla cupratow/Fe-SC (14 z 25 materialow) =>")
    print(f"           korekta zmienia je IDENTYCZNIE, nie dopasuje rozproszenia.")

    # ------------------------------------------------------------------------
    #  T2: global eta, Z = Z_heavy
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T2: global eta, Z = Z_heavy  (gdzie jest sygnal z ps23)")
    print("-" * 78)
    eta_T2, rss_T2 = fit_single_eta(dlog, Z_heavy)
    rms_T2 = float(np.sqrt(rss_T2 / n))
    F0, p0 = f_test(rss_zero, rss_T2, n, 0, 1)
    F1, p1 = f_test(rss_null_intercept, rss_T2, n, 1, 2)  # uwaga: mamy tylko eta
    print(f"  best eta = {eta_T2:+.4f}  (Z^2/137^2 dla Z=80: {80**2*ALPHA2:.4f})")
    print(f"  RMS_new = {rms_T2:.4f}")
    print(f"  redukcja vs intercept baseline: {(1-rms_T2/rms_null)*100:+.1f}%")
    print(f"  F(1,{n-1}) = {F0:.3f}  p={p0:.4f}{'*' if p0<0.05 else ''}")
    # per-class dlog po T2
    dlog_T2 = dlog_with_eta(dlog, Z_heavy, eta_T2)
    _show_per_class_means(dlog_T2, orbs, names, "po T2")

    # ------------------------------------------------------------------------
    #  T3: split eta per class (s, sp, d_cu, d_fe, d_other, f)
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T3: split eta per klasa orbital (6 parametrow)")
    print("-" * 78)
    class_names = ["s", "sp", "d_cu", "d_fe", "d_other", "f"]
    idx_map = {c: i for i, c in enumerate(class_names)}
    class_idx = np.array([idx_map[o] for o in orbs])
    etas_T3, rss_T3 = fit_multi_eta(dlog, Z_pair, class_idx, 6)
    rms_T3 = float(np.sqrt(rss_T3 / n))
    F, p = f_test(rss_null_intercept, rss_T3, n, 1, 6)
    print(f"  Fit eta na klase (Z=Z_pair):")
    for c, e in zip(class_names, etas_T3):
        print(f"     eta_{c:<8s} = {e:+9.4f}")
    print(f"  RMS_new = {rms_T3:.4f}   redukcja = {(1-rms_T3/rms_null)*100:+.1f}%")
    print(f"  F(5,{n-6}) = {F:.3f}  p = {p:.4f}{'*' if p<0.05 else ''}")

    # ------------------------------------------------------------------------
    #  T4: P7.5a  tylko eta_f (A_f correction dla 4f)
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T4: P7.5a — eta_f tylko (A_f reduction dla 4f)")
    print("-" * 78)
    mask_f = np.array([o == "f" for o in orbs])
    n_f = int(mask_f.sum())
    # fit tylko na 4f subset
    eta_f, _ = fit_single_eta(dlog[mask_f], Z_pair[mask_f])
    d_f_new = dlog_with_eta(dlog[mask_f], Z_pair[mask_f], eta_f)
    rms_f_before = rms_log(dlog[mask_f])
    rms_f_after = rms_log(d_f_new)
    # Zastosuj tez do pelnego CORE: ustaw eta=eta_f dla 4f, 0 dla reszty
    d_full = dlog.copy()
    d_full[mask_f] = d_f_new
    rss_T4 = float(np.sum(d_full ** 2))
    rms_T4 = float(np.sqrt(rss_T4 / n))
    F0, p0 = f_test(rss_zero, rss_T4, n, 0, 1)
    print(f"  4f subset: {[nm for nm,m in zip(names,mask_f) if m]}  (N={n_f})")
    print(f"  best eta_f = {eta_f:+.4f}  (dla Z=57 daje factor (1+eta*0.173)^2)")
    #                                         = (1 + eta_f * 0.173)
    factor_Z57 = (1.0 + eta_f * 57**2 * ALPHA2) ** 2
    Af_eff = 2.034 * np.sqrt(factor_Z57)
    print(f"  -> (1 + eta_f * Z^2/137^2)^2 dla Z=57 = {factor_Z57:.4f}")
    print(f"  -> A_f,eff = A_f * sqrt(factor) = 2.034 * {np.sqrt(factor_Z57):.4f}")
    print(f"             = {Af_eff:.4f}  (oryginalne A_f = 2.034)")
    print(f"  RMS na 4f subset: {rms_f_before:.4f} -> {rms_f_after:.4f}")
    print(f"  RMS na CORE N=25: {rms_null:.4f} -> {rms_T4:.4f}  "
          f"({(1-rms_T4/rms_null)*100:+.1f}%)")
    print(f"  F-test vs null=0: F(1,{n-1})={F0:.3f}  p={p0:.4f}{'*' if p0<0.05 else ''}")

    # ------------------------------------------------------------------------
    #  T5: P7.5b  tylko eta_6p (SOC correction dla 6p heavy spectators)
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T5: P7.5b — eta_6p tylko, applied via Z_heavy na cupratach Bi/Tl/Hg")
    print("-" * 78)
    # 6p-heavy = Z_heavy >= 80 i pair = Cu (cuprates) ORAZ elementarne Pb, Hg_elem
    mask_6p = np.array([
        (r[2] >= 80 and r[5] == "d_cu") or (r[0] in {"Pb", "Hg_elem"})
        for r in core
    ])
    n_6p = int(mask_6p.sum())
    eta_6p, _ = fit_single_eta(dlog[mask_6p], Z_heavy[mask_6p])
    d_6p_new = dlog_with_eta(dlog[mask_6p], Z_heavy[mask_6p], eta_6p)
    rms_6p_before = rms_log(dlog[mask_6p])
    rms_6p_after = rms_log(d_6p_new)
    d_full_5b = dlog.copy()
    d_full_5b[mask_6p] = d_6p_new
    rss_T5 = float(np.sum(d_full_5b ** 2))
    rms_T5 = float(np.sqrt(rss_T5 / n))
    F0, p0 = f_test(rss_zero, rss_T5, n, 0, 1)
    print(f"  6p-heavy subset: {[nm for nm,m in zip(names,mask_6p) if m]}  (N={n_6p})")
    print(f"  best eta_6p = {eta_6p:+.4f}  (dla Z=82: (Z alpha)^2 = "
          f"{82**2*ALPHA2:.4f})")
    factor_Z82 = (1.0 + eta_6p * 82**2 * ALPHA2) ** 2
    print(f"  -> (1 + eta_6p * Z^2/137^2)^2 dla Z=82 = {factor_Z82:.4f}")
    print(f"  RMS na 6p subset: {rms_6p_before:.4f} -> {rms_6p_after:.4f}")
    print(f"  RMS na CORE N=25: {rms_null:.4f} -> {rms_T5:.4f}  "
          f"({(1-rms_T5/rms_null)*100:+.1f}%)")
    print(f"  F-test vs null=0: F(1,{n-1})={F0:.3f}  p={p0:.4f}{'*' if p0<0.05 else ''}")

    # ------------------------------------------------------------------------
    #  T6: P7.5a + P7.5b razem (2 parametry)
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  T6: P7.5a + P7.5b lacznie (eta_f na 4f + eta_6p na 6p-heavy)")
    print("-" * 78)
    d_full_T6 = dlog.copy()
    d_full_T6[mask_f] = d_f_new
    d_full_T6[mask_6p] = d_6p_new
    rss_T6 = float(np.sum(d_full_T6 ** 2))
    rms_T6 = float(np.sqrt(rss_T6 / n))
    F, p = f_test(rss_null_intercept, rss_T6, n, k_null=1, k_alt=2)
    F0, p0 = f_test(rss_zero, rss_T6, n, 0, 2)
    print(f"  eta_f = {eta_f:+.4f}  (na 4f, N={n_f})")
    print(f"  eta_6p = {eta_6p:+.4f}  (na 6p-heavy, N={n_6p})")
    print(f"  RMS na CORE N=25: {rms_null:.4f} -> {rms_T6:.4f}  "
          f"({(1-rms_T6/rms_null)*100:+.1f}%)")
    print(f"  F(1,{n-2}) vs intercept = {F:.3f}  p = {p:.4f}{'*' if p<0.05 else ''}")
    print(f"  F(2,{n-2}) vs null=0    = {F0:.3f}  p = {p0:.4f}"
          f"{'*' if p0<0.05 else ''}")

    _show_per_class_means(d_full_T6, orbs, names, "po T6")

    # ------------------------------------------------------------------------
    #  Ekstrapolacja T6 na outliery P7
    # ------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("  Extrapolacja T6 na outliery P7 (nie powinna pomoc)")
    print("-" * 78)
    for nm, Zp, Zh, Mm, dl, orb in DATA:
        if nm not in OUTLIER_NAMES:
            continue
        dl_new = dl
        if orb == "f":
            dl_new = dl + 2.0 * np.log10(1.0 + eta_f * Zp**2 * ALPHA2)
        if (Zh >= 80 and orb == "d_cu") or (nm in {"Pb", "Hg_elem"}):
            dl_new = dl + 2.0 * np.log10(1.0 + eta_6p * Zh**2 * ALPHA2)
        print(f"  {nm:>10s}  Z_pair={Zp:3d}  Z_h={Zh:3d}  orb={orb:<8s}  "
              f"dlog {dl:+.3f} -> {dl_new:+.3f}")

    # ------------------------------------------------------------------------
    #  Podsumowanie
    # ------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Podsumowanie ps25")
    print("=" * 78)

    summary = [
        ("Baseline (intercept)",     1,  rms_null,                       "—"),
        ("T1 global eta, Z=Z_pair",  1,  rms_T1,                         f"{eta_T1:+.3f}"),
        ("T2 global eta, Z=Z_heavy", 1,  rms_T2,                         f"{eta_T2:+.3f}"),
        ("T3 split per-class (6)",   6,  rms_T3,                         "vector"),
        ("T4 P7.5a (eta_f only)",    1,  rms_T4,                         f"{eta_f:+.3f}"),
        ("T5 P7.5b (eta_6p only)",   1,  rms_T5,                         f"{eta_6p:+.3f}"),
        ("T6 P7.5a+b (2 params)",    2,  rms_T6,                         "—"),
    ]
    print(f"  {'Test':<28s}  {'k':>2s}  {'RMS_log':>8s}  {'red.':>6s}  {'eta':>10s}")
    for lbl, k, r, eta_str in summary:
        red = (1 - r / rms_null) * 100.0
        print(f"  {lbl:<28s}  {k:2d}  {r:8.4f}  {red:+5.1f}%  {eta_str:>10s}")

    print("""
  Interpretacja wyniku:
    T1 (Z=Z_pair globalny): nie dziala, bo dla 14/25 materialow Z=29 (Cu/Fe).
    T2 (Z=Z_heavy globalny): moze cos zmieni, ale znak eta jest kompromisem
       miedzy OVER-predict na 4f (Z_h=57) i UNDER-predict na 6p (Z_h=80).
    T3 (split per klasa): potwierdza ze potrzeba roznych eta dla roznych klas
       orbitalnych — eta_f i eta_s/d dla 6p maja przeciwne znaki.
    T6 (eta_f + eta_6p): dwa parametry, kazdy na swoim subsecie — to
       formalnie implementuje P7.5a + P7.5b.

  Wniosek fizyczny:
    Pojedyncza uniwersalna stala eta (1 parametr w Eq. 5) NIE zlikwiduje
    trendu. Forma (1 + eta Z^2/137^2) jest dobra — ale eta jest class-
    zalezne: eta_f < 0 (screening Kondo/Hund na 4f) vs eta_6p > 0 (SOC).

    To sugeruje rozszerzenie Eq. 5 do postaci:
        A_orb(Z) = A_orb^(0) * sqrt(1 + eta_orb * Z^2 / 137^2)
    gdzie eta_orb jest class-indexed (eta_s, eta_sp, eta_d, eta_f) — kazda
    klasa ma swoja wlasna stala sprzezenia, a eta_d pozostaje okolo zero
    (3d/4d bez relatywistycznej poprawki).

  Status P7.5:
    - P7.5a (A_f reduction dla 4f): eta_f = ~{} => A_f_eff = ~{:.2f}  (oryg. 2.034)
    - P7.5b (SOC dla 6p spectators): eta_6p = +{:.2f} => boost ~{:.2f}x na Z=82
""".format(f"{eta_f:+.3f}", 2.034 * np.sqrt((1 + eta_f * 57**2 * ALPHA2)**2 / 1.0),
           eta_6p, (1 + eta_6p * 82**2 * ALPHA2)**2))


def _show_per_class_means(d_new, orbs, names, label):
    print(f"\n  Per-class mean(dlog) {label}:")
    unique_classes = sorted(set(orbs))
    for c in unique_classes:
        mask = np.array([o == c for o in orbs])
        ds = np.asarray(d_new)[mask]
        if len(ds) == 0:
            continue
        mems = [nm for nm, m in zip(names, mask) if m]
        print(f"    {c:<9s} N={len(ds):2d}   mean={ds.mean():+.4f}  "
              f"RMS={rms_log(ds):.4f}   {mems[:4]}{'...' if len(mems)>4 else ''}")


if __name__ == "__main__":
    main()
