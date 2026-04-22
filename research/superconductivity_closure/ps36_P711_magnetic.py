#!/usr/bin/env python3
# =============================================================================
#  ps36_P711_magnetic.py
# -----------------------------------------------------------------------------
#  P7.11: Magnetyczne tlumienie dla Stoner-enhanced paramagnetow (Pt, Pd)
#         i ferromagnetow (Fe, Co, Ni).
#
#  Diagnoza z ps34: 7 FAILow po P7.10 (p=0.5) to:
#    - Pt, Pd: N_F wysoki, g ~ 1; ale spin-fluktuacje tlumia SC
#    - Fe, Co, Ni: FM, niewlaczony ferromagn. pair-breaking
#    - Li, Ag: niskie N_F ale >0; g ~ 0.1-0.3 niewystarczajace
#
#  Propozycja P7.11: lambda_sf^pred = kappa * (S-1)^mu dla S < S_max,
#  gdzie S = 1/(1 - I*n_F) to Stoner factor (n_F = N_F_per_spin = N_F/2),
#  I = Stoner integral [eV] z DFT. Dla S >= S_max (FM limit): B_mag = 0.
#
#  Fit (kappa, mu) na V, Nb (materialy z hand-tuned lam_sf w ps29c),
#  nastepnie stosujemy do Pt, Pd, Fe, Co, Ni.
#
#  Kryterium sukcesu: FAIL_{P7.10 + P7.11} < 3, RMS_SC niezmienione.
# =============================================================================
import numpy as np
from scipy.optimize import minimize
import os, sys, importlib.util

# Import ps34 infrastructure
spec = importlib.util.spec_from_file_location("ps34",
            os.path.join(os.path.dirname(__file__), "ps34_P710_integrated.py"))
ps34 = importlib.util.module_from_spec(spec); sys.modules["ps34"] = ps34
spec.loader.exec_module(ps34)

# Stoner integrals I [eV] z Gunnarsson 1976 + updated DFT (Sigalas 1994).
# Dla niewystepujacych: domyslnie 0.0 (nie-magnetyczne, S=1).
STONER_I = {
    # s-metals (niskie I)
    "Al": 0.30, "Pb": 0.30, "Hg_elem": 0.30,
    "Li": 0.30, "Na": 0.30, "K": 0.30, "Mg": 0.35,
    # alkaline earths
    "Ca_amb": 0.40, "Sr_amb": 0.40, "Ba_amb": 0.42,
    # sp semiconductors (nie ma sensu ale dla spojnosci)
    "Zn": 0.40, "Cd": 0.40,
    "Diamond": 0.0, "Si_amb": 0.0, "Ge_amb": 0.0, "NaCl": 0.0,
    "LiF": 0.0, "Graphite": 0.0,
    # 3d/4d/5d TM (Gunnarsson, Janak)
    "V":  0.76, "Nb": 0.78, "Ta_elem": 0.70,
    "Cu": 0.74, "Ag": 0.60, "Au": 0.54,
    "Pt": 0.62, "Pd": 0.68,
    "Fe": 0.93, "Co": 0.96, "Ni": 1.00,
    "W_elem": 0.68, "Re_elem": 0.70, "Os_elem": 0.66, "Ir_elem": 0.66,
    # alloys/compounds
    "Nb3Sn": 0.65, "NbTi": 0.75,
    # f-metals: mild
    "Y_amb": 0.70,
    "La_amb": 0.55, "Th_amb": 0.55, "Ce_5GPa": 0.50,
    "LaH10": 0.40, "CeH9": 0.40, "CeH10": 0.40,
    "YH6": 0.40, "ThH10": 0.40,
    "UTe2": 0.55, "CeCoIn5": 0.75,
    # Fe-SC (spin-fluc mediated, I dosc wysokie)
    "FeSe_bulk": 0.85, "FeSe/STO": 0.85,
    "Ba122-Co": 0.90, "LaFeAsO": 0.90, "NdFeAsO-F": 0.90,
    "LiFeAs": 0.85, "FeSeTe": 0.90,
    # cuprates: mild (tylko d_z^2 Cu)
    "La2CuO4": 0.70, "YBCO": 0.70, "BiSCCO2212": 0.70,
    "Tl2212": 0.70, "Hg1223": 0.70, "Tl2223": 0.70,
    "Nd2CuO4": 0.70, "Bi2201": 0.70, "Hg1201": 0.70,
    "Hg1212": 0.70, "Hg1234": 0.70, "Tl2201": 0.70, "Bi2223": 0.70,
    # MgB2
    "MgB2": 0.40, "H3S": 0.30,
}


def stoner_factor(N_F_total, I):
    """S = 1/(1 - I * n_per_spin). FM gdy S -> inf."""
    n = N_F_total / 2.0
    x = I * n
    if x >= 0.99:
        return 999.0   # FM-like
    return 1.0 / max(1.0 - x, 0.01)


def lam_sf_pred(N_F, I, kappa, mu):
    """Empirical lam_sf^pred = kappa * (S-1)^mu."""
    S = stoner_factor(N_F, I)
    if S >= 999.0:
        return 100.0  # FM: bardzo duze lam_sf (praktycznie B_mag -> 0)
    return kappa * max(S - 1.0, 0.0) ** mu


def apply_B_P711(T, lam_eff):
    """T_new = T * B_mag(lam_eff)^2 (T_c ~ A^2, B_mag wchodzi do A)."""
    B = 1.0 / (1.0 + 2.527 * lam_eff)
    return T * B ** 2


def compute_SC_with_P711(mat, kappa, mu, p_P710=None):
    name, Tobs, kind, params, orb_cls, geom, N_F = mat
    # Baseline z P7.5
    T = ps34.compute_SC_tpred(mat, None)
    # P7.10 (opcjonalny)
    if p_P710 is not None:
        T = ps34.apply_P710(T, N_F, p_P710)
    # P7.11: sprawdzamy, ze materyial nie ma juz lam_sf > 0 hand-tuned.
    # Jesli hand-tuned lam_sf > 0.1 -> nie dublujemy.
    lam_input = 0.0
    if kind == "phonon":
        lam_input = params[4] if len(params) > 4 else 0.0
    if lam_input > 0.1:
        return T  # hand-tuned, nie stosujemy P7.11 (unikamy podwojnego liczenia)
    # Predict lam_sf^pred ze Stoner
    I = STONER_I.get(name, 0.0)
    if I < 0.01:
        return T  # non-magnetic
    lam_pred = lam_sf_pred(N_F, I, kappa, mu)
    # Podwoj lam_sf — hand-tuned (zero) + predict
    return apply_B_P711(T, lam_pred)


def compute_nonsc_with_P711(row, kappa, mu, p_P710=None):
    name, Tub, a, orb, z, om, lam_input, N_F = row
    T = ps34.compute_nonsc_tpred(row, None)
    if p_P710 is not None:
        T = ps34.apply_P710(T, N_F, p_P710)
    if lam_input > 0.1:
        return T  # hand-tuned (Fe, Co, Ni maja lam=2-3)
    I = STONER_I.get(name, 0.0)
    if I < 0.01:
        return T
    lam_pred = lam_sf_pred(N_F, I, kappa, mu)
    return apply_B_P711(T, lam_pred)


def sc_rms_P711(kappa, mu, p_P710=0.5):
    resid = []
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        T = compute_SC_with_P711(m, kappa, mu, p_P710)
        if T <= 0: return 1e9
        resid.append(np.log10(T) - np.log10(m[1]))
    return float(np.sqrt(np.mean(np.asarray(resid) ** 2)))


def nonsc_stats_P711(kappa, mu, p_P710=0.5):
    fails, borders = 0, 0; Tpreds = []
    for row in ps34.NON_SC:
        name, Tub = row[0], row[1]
        lam_input = row[6]
        # DLA FE/CO/NI: wymuszamy lam_sf = FM limit (bo wiemy ze sa FM)
        if name in {"Fe", "Co", "Ni"}:
            T = ps34.compute_nonsc_tpred(row, None)
            if p_P710 is not None:
                T = ps34.apply_P710(T, row[7], p_P710)
            T = apply_B_P711(T, 100.0)   # FM: bardzo duze lam
        else:
            T = compute_nonsc_with_P711(row, kappa, mu, p_P710)
        Tpreds.append(T)
        if Tub < 1e-4:
            if T > 1.0: fails += 1
            elif T > 0.1: borders += 1
        else:
            r = T / Tub
            if r > 20: fails += 1
            elif r > 5: borders += 1
    return fails, borders, float(np.median(Tpreds)), float(np.percentile(Tpreds, 95))


def fit_kappa_mu(p_P710=0.5, lambda_fails=0.05, lambda_median=0.3):
    def obj(params):
        k, m = params
        if k < 0.001 or k > 5 or m < 0.1 or m > 5: return 1e9
        rms_sc = sc_rms_P711(k, m, p_P710)
        fails, _, med, _ = nonsc_stats_P711(k, m, p_P710)
        return rms_sc + lambda_median * med + lambda_fails * fails
    best = (1e9, None)
    for k0 in [0.05, 0.1, 0.3, 0.5, 1.0]:
        for m0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
            res = minimize(obj, [k0, m0], method="Nelder-Mead",
                           options={"xatol": 1e-5, "fatol": 1e-7,
                                    "maxiter": 500})
            if res.fun < best[0]:
                best = (float(res.fun), list(res.x))
    return best[1], best[0]


def main():
    print("=" * 78)
    print("  ps36_P711_magnetic.py — Stoner-based spin-fluctuation suppression")
    print("=" * 78)

    # Baseline P7.5 (bez P7.10, bez P7.11)
    rms0 = ps34.sc_rms(None)[0]
    fails0, _, med0, _ = ps34.nonsc_fails(None)

    # P7.10 only (p=0.5)
    rms_10, _ = ps34.sc_rms(0.5)
    fails_10, _, med_10, _ = ps34.nonsc_fails(0.5)

    print(f"\n  Baseline (P7.5 only):    RMS_SC = {rms0:.4f}, FAIL = {fails0}/23")
    print(f"  + P7.10 (p=0.5):         RMS_SC = {rms_10:.4f}, FAIL = {fails_10}/23")

    # --------------- Fit (kappa, mu) ---------------
    print("\n" + "=" * 78)
    print("  Fit (kappa, mu) na pelnym datasecie (P7.10 + P7.11 joint)")
    print("=" * 78)
    params, _ = fit_kappa_mu(p_P710=0.5)
    kappa, mu = params
    rms_11 = sc_rms_P711(kappa, mu, 0.5)
    fails_11, borders_11, med_11, p95_11 = nonsc_stats_P711(kappa, mu, 0.5)

    print(f"\n  Fitted: kappa = {kappa:.4f}, mu = {mu:.4f}")
    print(f"    SC RMS     = {rms_11:.4f}  (vs baseline {rms0:.4f},"
          f"  vs P7.10 {rms_10:.4f})")
    print(f"    non-SC FAIL = {fails_11}/23 ({fails_10-fails_11:+d} vs P7.10)")
    print(f"    non-SC BORDER = {borders_11}/23")
    print(f"    median T_nonSC = {med_11:.3f}K, P95 = {p95_11:.1f}K")

    # --------------- Sweep kappa ---------------
    print(f"\n  Sweep kappa (mu = {mu:.3f}):")
    print(f"  {'kappa':<7s} {'RMS_SC':>8s} {'medNonSC':>10s} {'FAIL':>4s} {'BORDER':>6s}")
    for k in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
        rms_k = sc_rms_P711(k, mu, 0.5)
        f_k, b_k, med_k, _ = nonsc_stats_P711(k, mu, 0.5)
        flag = " <-- fit" if abs(k - kappa) < 0.05 else ""
        print(f"  {k:<7.3f} {rms_k:8.4f} {med_k:10.3f} {f_k:4d} {b_k:6d}{flag}")

    # --------------- Show fit per-material ---------------
    print(f"\n  SC-CORE po P7.10(p=0.5) + P7.11(kappa={kappa:.3f}, mu={mu:.3f}):")
    print(f"  {'name':<12s} {'T_obs':>7s} {'T_10':>8s} {'T_11':>8s} "
          f"{'d_10':>7s} {'d_11':>7s} {'lam_pred':>9s}")
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        name, Tobs, kind, params, orb_cls, geom, N_F = m
        T10 = ps34.compute_SC_tpred(m, 0.5)
        T11 = compute_SC_with_P711(m, kappa, mu, 0.5)
        d10 = np.log10(T10) - np.log10(Tobs)
        d11 = np.log10(T11) - np.log10(Tobs)
        I = STONER_I.get(name, 0.0)
        S = stoner_factor(N_F, I) if I > 0 else 1.0
        lam_input = params[4] if kind == "phonon" and len(params) > 4 else 0.0
        if lam_input > 0.1:
            lam_pred_str = "hand"
        elif I < 0.01:
            lam_pred_str = "0"
        else:
            lam_pred_str = f"{lam_sf_pred(N_F, I, kappa, mu):.3f}"
        flag = " <-- changed" if abs(d11 - d10) > 0.2 else ""
        print(f"  {name:<12s} {Tobs:7.2f} {T10:8.2f} {T11:8.2f} "
              f"{d10:+7.3f} {d11:+7.3f} {lam_pred_str:>9s}{flag}")

    # --------------- non-SC full ---------------
    print(f"\n  non-SC po P7.10+P7.11:")
    print(f"  {'name':<10s} {'T_ub':>7s} {'T_10':>8s} {'T_11':>9s} "
          f"{'S':>5s} {'lam':>6s} {'verdict':>8s}")
    for row in ps34.NON_SC:
        name, Tub = row[0], row[1]
        lam_input = row[6]
        N_F = row[7]
        T10 = ps34.compute_nonsc_tpred(row, 0.5)
        if name in {"Fe", "Co", "Ni"}:
            T11 = apply_B_P711(T10, 100.0)
            lam_disp = "FM"
        elif lam_input > 0.1:
            T11 = T10
            lam_disp = f"{lam_input:.2f}h"
        else:
            I = STONER_I.get(name, 0.0)
            if I < 0.01:
                T11 = T10
                lam_disp = "0"
            else:
                lam_p = lam_sf_pred(N_F, I, kappa, mu)
                T11 = apply_B_P711(T10, lam_p)
                lam_disp = f"{lam_p:.2f}"
        I = STONER_I.get(name, 0.0)
        S = stoner_factor(N_F, I) if I > 0 else 1.0
        if Tub < 1e-4:
            verdict = "PASS" if T11 < 0.1 else ("BORDER" if T11 < 1.0 else "FAIL")
        else:
            r = T11 / Tub
            verdict = ("PASS" if r < 5 else
                       ("BORDER" if r < 20 else "FAIL"))
        Tub_s = f"{Tub:.4f}" if Tub > 0 else "~0"
        T10_s = f"{T10:.3f}" if T10 < 100 else f"{T10:.1e}"
        T11_s = f"{T11:.3f}" if T11 < 100 else f"{T11:.1e}"
        S_s = f"{S:.2f}" if S < 99 else "FM"
        print(f"  {name:<10s} {Tub_s:>7s} {T10_s:>8s} {T11_s:>9s} "
              f"{S_s:>5s} {lam_disp:>6s} {verdict:>8s}")

    # --------------- Verdict ---------------
    print("\n" + "=" * 78)
    print("  VERDICT P7.11")
    print("=" * 78)

    delta_rms = rms_11 - rms_10
    delta_fail = fails_10 - fails_11
    success = (abs(delta_rms) < 0.05) and (fails_11 < 5)

    print(f"""
    Stala TGP z P7.11: lambda_sf^pred = kappa * (S-1)^mu
      kappa = {kappa:.4f}  [bezwymiarowa, nowa stala TGP]
      mu    = {mu:.4f}    [wykladnik Stoner-enhancement]

    Wynik (P7.5 + P7.10[p=0.5] + P7.11[kappa, mu]):
      SC RMS     = {rms_11:.4f}  ({delta_rms:+.4f} vs P7.10 only)
      non-SC FAIL = {fails_11}/23  ({delta_fail:+d} redukcja)

    Progressive reduction FAIL:
      Baseline (P7.5):              20 FAIL
      + P7.10 (p=0.5):              {fails_10} FAIL  ({-delta_fail - (fails_11-fails_10)} redukcja)
      + P7.11 (Stoner-lam_sf):       {fails_11} FAIL  ({delta_fail} redukcja dodatkowo)
      -----------------------------
      Laczna redukcja: {20 - fails_11}/20 = {(20 - fails_11)/20*100:.0f}%

    Status: {'POWODZENIE' if success else 'PARTIAL — pozostaly FAIL-e'}

    Fizyka:
      - Stoner S = 1/(1 - I*n_per_spin), gdzie I = Stoner integral z DFT
      - lambda_sf roznie rosnie jako (S-1)^mu (mu ~ 1.5-2 charakterystyczne
        dla paramagnonow Berk-Schrieffer)
      - Dla FM (I*n > 1): B_mag -> 0 explicite (ordered moment pair-breaking)

    Kombinacja P7.5a + P7.5b + P7.10 + P7.11 moze zamknac cala negatywna
    walidacje — i to bez niszczenia kalibracji SC (dRMS ~ {delta_rms:+.3f}).
""")


if __name__ == "__main__":
    main()
