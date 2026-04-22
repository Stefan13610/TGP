#!/usr/bin/env python3
# =============================================================================
#  ps26_validation_P75ab_N29.py
# -----------------------------------------------------------------------------
#  Walidacja korekt P7.5a (eta_f=-1.308) + P7.5b (eta_6p=+0.331) na pelnym
#  N=29 datasecie. Baseline z ps17: r=0.8752, RMS_log=0.3571.
#
#  Oczekiwanie:
#    - 4f-class (LaH10, CeH9, CeH10, La_amb + Th_amb + Ce_5GPa) przesuwa sie
#      do zera dzieki P7.5a
#    - 6p-heavy (BSCCO2212, Tl2212, Hg1223, Tl2223, Bi2201, Pb, Hg_elem)
#      przesuwa sie blizej zera dzieki P7.5b
#    - Y_amb i H3S pozostaja bez zmian (nie sa w ani jednej klasie)
#    - Global r powinien wzrosnac, RMS_log zmalec
#
#  Dodatkowo:
#    - Leave-one-out cross-validation eta_f na 4f core subset (N=4)
#      -> czy eta_f jest stabilne?
#    - Per-class residualna analiza (d_fe, d_other) -> co jeszcze zostalo?
# =============================================================================
import numpy as np
from scipy.stats import pearsonr
from scipy.optimize import minimize_scalar

ALPHA2 = 1.0 / 137.036 ** 2
ETA_F  = -1.308   # P7.5a fit z ps25
ETA_6P = +0.331   # P7.5b fit z ps25

# Pelny dataset N=29 z ps17_results.txt
# (name, Z_pair, Z_heavy, orb_class, T_obs, T_pred)
DATA = [
    # cuprates
    ("La2CuO4",     29, 57, "d_cu",    38.00,  53.27),
    ("YBCO",        29, 56, "d_cu",    92.00,  75.59),
    ("BiSCCO2212",  29, 83, "d_cu",    85.00,  75.59),
    ("Tl2212",      29, 81, "d_cu",   108.00,  75.76),
    ("Hg1223",      29, 80, "d_cu",   138.00,  92.82),
    ("Tl2223",      29, 81, "d_cu",   125.00,  92.58),
    ("Nd2CuO4",     29, 60, "d_cu",    24.00,  53.86),
    ("Bi2201",      29, 83, "d_cu",    34.00,  53.41),
    # phonon + elementary + Fe-SC + hydridy + P7 outliery
    ("Al",          13, 13, "s",        1.18,   2.96),
    ("Pb",          82, 82, "s",        7.20,   4.96),
    ("Nb",          41, 41, "d_other",  9.26,  13.83),
    ("V",           23, 23, "d_other",  5.30,  10.95),
    ("Hg_elem",     80, 80, "s",        4.15,   2.52),
    ("MgB2",         5, 12, "sp",      39.00,  31.37),
    ("FeSe_bulk",   26, 34, "d_fe",     8.00,   7.85),
    ("FeSe/STO",    26, 34, "d_fe",    65.00,  57.23),
    ("Ba122-Co",    26, 56, "d_fe",    22.00,  17.90),
    ("LaFeAsO",     26, 57, "d_fe",    26.00,  16.35),
    ("NdFeAsO-F",   26, 60, "d_fe",    55.00,  24.15),
    ("Nb3Sn",       41, 50, "d_other", 18.30,  22.03),
    ("NbTi",        41, 41, "d_other", 10.00,  17.24),
    ("H3S",         16, 16, "sp",     203.00,  76.05),
    ("LaH10",       57, 57, "f",      250.00, 368.08),
    ("CeH9",        58, 58, "f",      100.00, 143.13),
    ("CeH10",       58, 58, "f",      115.00, 187.48),
    ("La_amb",      57, 57, "f",        6.00,  14.34),
    ("Y_amb",       39, 39, "d_other",  1.30,  18.46),
    ("Th_amb",      90, 90, "f",        1.38,  10.40),
    ("Ce_5GPa",     58, 58, "f",        1.70,   6.02),
]


def correction_P75a(Z, eta_f=ETA_F):
    """Kwadrat czynnika mnozacego T_pred (bo A wchodzi kwadratowo)."""
    return (1.0 + eta_f * Z**2 * ALPHA2) ** 2


def correction_P75b(Z, eta_6p=ETA_6P):
    return (1.0 + eta_6p * Z**2 * ALPHA2) ** 2


def is_6p_heavy(name, Z_heavy, orb):
    """6p-heavy: cuprat z ciezkim kationem (Bi/Tl/Hg spectator) lub
    elementarne Pb/Hg_elem."""
    if name in {"Pb", "Hg_elem"}:
        return True
    if orb == "d_cu" and Z_heavy >= 80:
        return True
    return False


def apply_corrections(data, eta_f=ETA_F, eta_6p=ETA_6P):
    """Zwraca liste (name, Z_p, Z_h, orb, T_obs, T_pred_new, dlog_new,
    applied_flag)."""
    out = []
    for name, Zp, Zh, orb, Tobs, Tpred in data:
        factor = 1.0
        flags = []
        # P7.5a na f-orbital (dotyczy Z_pair 4f/5f)
        if orb == "f":
            factor *= correction_P75a(Zp, eta_f)
            flags.append("7.5a")
        # P7.5b na 6p-heavy
        if is_6p_heavy(name, Zh, orb):
            factor *= correction_P75b(Zh, eta_6p)
            flags.append("7.5b")
        Tpred_new = Tpred * factor
        dlog_new = np.log10(Tpred_new) - np.log10(Tobs)
        out.append((name, Zp, Zh, orb, Tobs, Tpred_new,
                    float(dlog_new), ",".join(flags) if flags else "-"))
    return out


def global_metrics(corrected):
    """Wylicz Pearson r i RMS_log na log10(T_obs) vs log10(T_pred)."""
    log_obs  = np.array([np.log10(r[4]) for r in corrected])
    log_pred = np.array([np.log10(r[5]) for r in corrected])
    dlog     = log_pred - log_obs
    r, p_r = pearsonr(log_obs, log_pred)
    rms = float(np.sqrt(np.mean(dlog ** 2)))
    return float(r), float(p_r), rms


def per_class_stats(corrected):
    by_class = {}
    for row in corrected:
        orb = row[3]
        by_class.setdefault(orb, []).append(row[6])
    print(f"\n  {'klasa':<10s}  {'N':>3s}  {'mean':>8s}  {'RMS':>8s}  "
          f"{'max|dlog|':>10s}")
    order = ["d_cu", "d_fe", "d_other", "s", "sp", "f"]
    for orb in order:
        if orb not in by_class:
            continue
        ds = np.array(by_class[orb])
        print(f"  {orb:<10s}  {len(ds):3d}  {ds.mean():+.4f}  "
              f"{np.sqrt((ds**2).mean()):.4f}  {np.abs(ds).max():10.4f}")


def loo_eta_f(data):
    """Leave-one-out cross-validation dla eta_f na 4f core subset."""
    core_4f = [r for r in data if r[3] == "f" and r[0] not in
               {"Th_amb", "Ce_5GPa"}]
    # LaH10, CeH9, CeH10, La_amb -- ten sam co w ps25
    names = [r[0] for r in core_4f]
    Zps   = np.array([r[1] for r in core_4f], float)
    Tobs  = np.array([r[4] for r in core_4f], float)
    Tpred = np.array([r[5] for r in core_4f], float)
    dlog0 = np.log10(Tpred) - np.log10(Tobs)

    print(f"\n  LOO on 4f core (N={len(core_4f)}):")
    print(f"  {'held-out':<10s}  {'eta_f fit':>10s}  {'A_f_eff@Z':>12s}  "
          f"{'dlog held (pred)':>18s}")
    fits = []
    for i in range(len(core_4f)):
        mask = np.ones(len(core_4f), dtype=bool)
        mask[i] = False
        # fit eta na reszcie
        def loss(eta):
            factor = 1.0 + eta * Zps[mask]**2 * ALPHA2
            factor = np.where(factor > 1e-6, factor, 1e-6)
            d = dlog0[mask] + 2.0 * np.log10(factor)
            return float(np.sum(d ** 2))
        res = minimize_scalar(loss, bounds=(-5, 5), method="bounded")
        eta_loo = float(res.x)
        fits.append(eta_loo)
        # predykcja na held-out
        Zh = Zps[i]
        factor_held = (1.0 + eta_loo * Zh**2 * ALPHA2) ** 2
        Af_eff = 2.034 * np.sqrt(factor_held)
        dlog_held_corr = dlog0[i] + np.log10(factor_held)
        print(f"  {names[i]:<10s}  {eta_loo:+10.4f}  {Af_eff:12.4f}  "
              f"{dlog_held_corr:+18.4f}")
    mean_fit = float(np.mean(fits))
    std_fit  = float(np.std(fits, ddof=1))
    print(f"\n  LOO mean eta_f = {mean_fit:+.4f}  +- {std_fit:.4f}")
    print(f"  full fit eta_f = {ETA_F:+.4f}  (ps25)")
    print(f"  spread/ful = {std_fit/abs(mean_fit)*100:.1f}%  "
          f"(nizej = stabilniejszy)")


def print_big_table(corrected, label):
    print(f"\n{'-'*78}\n  {label}\n{'-'*78}")
    print(f"  {'material':<12s}  {'Z_p':>4s}  {'Z_h':>4s}  {'orb':<8s}  "
          f"{'T_obs':>7s}  {'T_pred':>8s}  {'dlog':>7s}  {'flag':<10s}")
    for name, Zp, Zh, orb, Tobs, Tpred, dlog, flag in corrected:
        mark = "*" if abs(dlog) > 0.4 else " "
        print(f"  {name:<12s}  {Zp:4d}  {Zh:4d}  {orb:<8s}  "
              f"{Tobs:7.2f}  {Tpred:8.2f}  {dlog:+.3f}{mark} {flag:<10s}")


def outliers_list(corrected, thresh=0.4):
    return [r for r in corrected if abs(r[6]) > thresh]


def main():
    print("=" * 78)
    print("  ps26_validation_P75ab_N29.py")
    print("  Walidacja P7.5a+b na pelnym N=29 datasecie (baseline ps17)")
    print("=" * 78)
    print(f"\n  eta_f  = {ETA_F:+.4f}  (P7.5a, 4f Kondo/Hund screening)")
    print(f"  eta_6p = {ETA_6P:+.4f}  (P7.5b, 6p SOC boost)")

    # ------------------------------------------------------------------
    #  Baseline N=29 (bez korekty)
    # ------------------------------------------------------------------
    baseline = [(name, Zp, Zh, orb, Tobs, Tpred,
                 float(np.log10(Tpred) - np.log10(Tobs)), "-")
                for name, Zp, Zh, orb, Tobs, Tpred in DATA]
    r0, p0, rms0 = global_metrics(baseline)
    print(f"\n  Baseline N=29:")
    print(f"    r(log)   = {r0:.4f}   (p = {p0:.3e})")
    print(f"    RMS_log  = {rms0:.4f}")
    out0 = outliers_list(baseline)
    print(f"    outliers (|dlog|>0.4): {len(out0)}")
    for r in out0:
        print(f"       {r[0]:<12s}  dlog = {r[6]:+.3f}")

    per_class_stats(baseline)

    # ------------------------------------------------------------------
    #  P7.5a+b correction applied
    # ------------------------------------------------------------------
    corrected = apply_corrections(DATA)
    r1, p1, rms1 = global_metrics(corrected)
    print(f"\n" + "=" * 78)
    print(f"  Po zastosowaniu P7.5a + P7.5b")
    print("=" * 78)
    print(f"\n  N=29 (wszystkie materialy):")
    print(f"    r(log)   = {r1:.4f}   (baseline {r0:.4f}, zmiana {r1-r0:+.4f})")
    print(f"    RMS_log  = {rms1:.4f}   (baseline {rms0:.4f}, "
          f"redukcja {(1-rms1/rms0)*100:+.1f}%)")
    out1 = outliers_list(corrected)
    print(f"    outliers (|dlog|>0.4): {len(out1)}  "
          f"(bylo {len(out0)})")
    for r in out1:
        print(f"       {r[0]:<12s}  dlog = {r[6]:+.3f}  flags: {r[7]}")

    per_class_stats(corrected)

    # ------------------------------------------------------------------
    #  Szczegolowa tabela
    # ------------------------------------------------------------------
    print_big_table(corrected, "Pelna tabela po P7.5a+b")

    # ------------------------------------------------------------------
    #  LOO cross-validation eta_f
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Leave-one-out cross-validation eta_f na 4f core (N=4)")
    print("=" * 78)
    loo_eta_f(DATA)

    # ------------------------------------------------------------------
    #  Analiza residualna — co jeszcze zostalo?
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Residualne kanaly po P7.5a+b")
    print("=" * 78)
    print("""
  Klasy, gdzie mean(dlog) wciaz daleko od 0:
""")
    # wyswietl residualne sygnaly
    by_class = {}
    for row in corrected:
        by_class.setdefault(row[3], []).append(row[6])
    signals = []
    for orb, ds in by_class.items():
        ds = np.array(ds)
        if len(ds) < 2:
            continue
        if abs(ds.mean()) > 0.08:
            signals.append((orb, len(ds), float(ds.mean()),
                            float(np.sqrt((ds**2).mean()))))

    if signals:
        print(f"  {'klasa':<10s}  {'N':>3s}  {'mean':>8s}  {'RMS':>8s}  "
              f"sugestia")
        suggestions = {
            "d_cu":    "cuprat bias po P7.5b, sprawdzic dodatkowy Fuchs-Kliewer",
            "d_fe":    "Fe-SC pod-szacowane  nesting/Hund/multi-orbital -> P7.6?",
            "d_other": "phonon-TM nad-szacowane  Eliashberg mu* ? -> P7.7?",
            "s":       "s-orbital 2 skale (Al vs Pb/Hg) - rozproszenie",
        }
        for orb, N, m, rms_c in signals:
            sug = suggestions.get(orb, "(niezdefiniowane)")
            print(f"  {orb:<10s}  {N:3d}  {m:+.4f}  {rms_c:.4f}  {sug}")
    else:
        print("  Brak residualnych sygnalow > 0.08 mean.")

    # ------------------------------------------------------------------
    #  Podsumowanie
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  Podsumowanie ps26")
    print("=" * 78)
    print(f"""
  Baseline (ps17):
    r = {r0:.4f}   RMS_log = {rms0:.4f}   outliers = {len(out0)}

  Po P7.5a + P7.5b:
    r = {r1:.4f}   RMS_log = {rms1:.4f}   outliers = {len(out1)}

  Zmiana:
    delta r    = {r1-r0:+.4f}
    delta RMS  = {rms1-rms0:+.4f}   ({(1-rms1/rms0)*100:+.1f}%)
    outliers:  {len(out0)} -> {len(out1)}  ({len(out0)-len(out1)} zamknietych)

  Werdykt:
    - Jesli r wzroslo o >= +0.03 i outliers spadly o >= 1 przy JEDNEJ nowej
      stalej (eta_f alone, bo eta_6p ma marginalne efekty) -> P7.5a jest
      zaakceptowane do paperu SC v2.
    - Jesli Th_amb przestal byc outlierem po P7.5a -> silna walidacja.
    - Rezidualne sygnaly w d_fe i d_other wskazuja kolejne potencjalne
      sub-problemy P7.6/P7.7 (spoza P7.5 Z-relatywistycznego).

  Uwaga: ta walidacja jest self-consistent (eta_f fitowane na 4 z 6
  materialow f-klasy, potem zastosowane do wszystkich 6 vs baseline).
  Th_amb i Ce_5GPa NIE byly uzyte do fitu eta_f w ps25 (byly outlierami)
  -> ich poprawa jest WALIDACJA HOLD-OUT, nie fit.
""")


if __name__ == "__main__":
    main()
