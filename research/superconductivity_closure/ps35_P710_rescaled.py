#!/usr/bin/env python3
# =============================================================================
#  ps35_P710_rescaled.py
# -----------------------------------------------------------------------------
#  ps34 pokazal: proste T_new = T_old * g(N_F/N_F_ref) psuje SC bo istniejaca
#  kalibracja zaklada <g>=1. Potrzebny joint fit (A, p):
#      T_new = T_old * A * g(N_F/N_F_ref)
#  gdzie A kompensuje srednia g-reduction.
#
#  Test: fituj (A, p) na SC CORE (N=38, bez off-model), sprawdz non-SC.
#  Oczekiwanie: A ~ 1/<g_SC> ~ 2 dla p ~ 1.
# =============================================================================
import numpy as np
from scipy.optimize import minimize
from scipy.stats import pearsonr

# Importuj wszystko z ps34
import importlib.util, os, sys
spec = importlib.util.spec_from_file_location("ps34",
            os.path.join(os.path.dirname(__file__), "ps34_P710_integrated.py"))
ps34 = importlib.util.module_from_spec(spec); sys.modules["ps34"] = ps34
spec.loader.exec_module(ps34)


def compute_SC_tpred_rescaled(mat, A, p):
    T = ps34.compute_SC_tpred(mat, None)  # baseline P7.5
    g = ps34.g_powerlaw(mat[6] / ps34.N_F_REF, p)
    return T * A * g


def compute_nonsc_tpred_rescaled(row, A, p):
    T = ps34.compute_nonsc_tpred(row, None)
    g = ps34.g_powerlaw(row[7] / ps34.N_F_REF, p)
    return T * A * g


def sc_rms_AP(A, p):
    resid = []
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        T = compute_SC_tpred_rescaled(m, A, p)
        if T <= 0: return 1e9
        resid.append(np.log10(T) - np.log10(m[1]))
    return float(np.sqrt(np.mean(np.asarray(resid) ** 2)))


def nonsc_stats_AP(A, p):
    fails, borders = 0, 0; Tp = []
    for row in ps34.NON_SC:
        name, Tub = row[0], row[1]
        T = compute_nonsc_tpred_rescaled(row, A, p)
        Tp.append(T)
        if Tub < 1e-4:
            if T > 1.0: fails += 1
            elif T > 0.1: borders += 1
        else:
            r = T / Tub
            if r > 20: fails += 1
            elif r > 5: borders += 1
    return fails, borders, float(np.median(Tp)), float(np.percentile(Tp, 95))


def fit_AP(lambda_fails=0.05, lambda_median=0.3):
    def obj(params):
        A, p = params
        if A <= 0.1 or A > 20: return 1e9
        if p <= 0.05 or p > 5: return 1e9
        rms_sc = sc_rms_AP(A, p)
        fails, _, med, _ = nonsc_stats_AP(A, p)
        return rms_sc + lambda_median * med + lambda_fails * fails

    best = (1e9, None)
    for A0 in [1.0, 1.5, 2.0, 3.0]:
        for p0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
            res = minimize(obj, [A0, p0], method="Nelder-Mead",
                           options={"xatol": 1e-5, "fatol": 1e-7,
                                    "maxiter": 500})
            if res.fun < best[0]:
                best = (float(res.fun), list(res.x))
    return best[1], best[0]


def fit_AP_pureSC():
    """Fit (A, p) minimalizujac TYLKO RMS_SC (bez regularizacji non-SC)."""
    def obj(params):
        A, p = params
        if A <= 0.1 or A > 20 or p <= 0.05 or p > 5: return 1e9
        return sc_rms_AP(A, p)
    best = (1e9, None)
    for A0 in [1.0, 2.0, 3.0, 5.0]:
        for p0 in [0.5, 1.0, 1.5, 2.0]:
            res = minimize(obj, [A0, p0], method="Nelder-Mead",
                           options={"xatol": 1e-5, "fatol": 1e-7})
            if res.fun < best[0]:
                best = (float(res.fun), list(res.x))
    return best[1], best[0]


def main():
    print("=" * 78)
    print("  ps35_P710_rescaled.py — joint fit (A, p) z P7.10")
    print("=" * 78)

    # Baseline
    rms0 = ps34.sc_rms(None)[0]
    fails0, _, med0, _ = ps34.nonsc_fails(None)
    print(f"\n  Baseline (P7.5 only, no P7.10):")
    print(f"    SC RMS = {rms0:.4f}, non-SC FAIL = {fails0}/{len(ps34.NON_SC)}, "
          f"median = {med0:.2f}K")

    # Fit joint (A, p) minimalizujac wazona sume SC+nonSC
    print(f"\n  (1) Fit joint (A, p) z wazonym obj (SC + non-SC):")
    params, _ = fit_AP(lambda_fails=0.05, lambda_median=0.3)
    A, p = params
    rms_fit = sc_rms_AP(A, p)
    fails_fit, borders_fit, med_fit, p95_fit = nonsc_stats_AP(A, p)
    print(f"    A = {A:.3f}, p = {p:.3f}")
    print(f"    SC RMS = {rms_fit:.4f} ({rms_fit-rms0:+.4f} vs baseline)")
    print(f"    non-SC FAIL = {fails_fit}/{len(ps34.NON_SC)}, "
          f"median = {med_fit:.3f}K, P95 = {p95_fit:.1f}K")

    # Fit PURE SC (tylko RMS)
    print(f"\n  (2) Fit pure SC (tylko RMS_SC, zero reg):")
    pars_pure, _ = fit_AP_pureSC()
    A_pure, p_pure = pars_pure
    rms_pure = sc_rms_AP(A_pure, p_pure)
    fails_pure, _, med_pure, _ = nonsc_stats_AP(A_pure, p_pure)
    print(f"    A = {A_pure:.3f}, p = {p_pure:.3f}")
    print(f"    SC RMS = {rms_pure:.4f} ({rms_pure-rms0:+.4f} vs baseline)")
    print(f"    non-SC FAIL = {fails_pure}/{len(ps34.NON_SC)}, "
          f"median = {med_pure:.3f}K")

    # Grid sweep for visualization
    print(f"\n  (3) Grid (A, p): tabela RMS_SC  |  #FAIL non-SC")
    ps = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    As = [1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0]
    print(f"    {'A\\p':<5s}", end="")
    for pp in ps: print(f"{pp:>12.2f}", end="")
    print()
    for A_ in As:
        print(f"    {A_:<5.1f}", end="")
        for pp in ps:
            rms_ = sc_rms_AP(A_, pp)
            fails_, _, _, _ = nonsc_stats_AP(A_, pp)
            print(f"   {rms_:5.3f}/{fails_:2d}", end="")
        print()

    # Per-material na best-joint fit
    print(f"\n  (4) SC per-material przy (A={A:.2f}, p={p:.2f}):")
    print(f"  {'name':<12s} {'T_obs':>7s} {'T_base':>8s} {'T_new':>8s} "
          f"{'d_base':>8s} {'d_new':>7s}")
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        T0 = ps34.compute_SC_tpred(m, None)
        T_new = compute_SC_tpred_rescaled(m, A, p)
        d0 = np.log10(T0) - np.log10(m[1])
        d_new = np.log10(T_new) - np.log10(m[1])
        flag = " <-- out" if abs(d_new) > 0.4 else ""
        print(f"  {m[0]:<12s} {m[1]:7.2f} {T0:8.2f} {T_new:8.2f} "
              f"{d0:+8.3f} {d_new:+7.3f}{flag}")

    # non-SC per-material
    print(f"\n  (5) non-SC (A={A:.2f}, p={p:.2f}):")
    print(f"  {'name':<10s} {'T_ub':>7s} {'T_base':>8s} {'T_new':>9s} "
          f"{'verdict':>8s}")
    for row in ps34.NON_SC:
        T0 = ps34.compute_nonsc_tpred(row, None)
        T_new = compute_nonsc_tpred_rescaled(row, A, p)
        Tub = row[1]
        if Tub < 1e-4:
            verdict = "PASS" if T_new < 0.1 else ("BORDER" if T_new < 1.0 else "FAIL")
        else:
            r = T_new / Tub
            verdict = ("PASS" if r < 5 else
                       ("BORDER" if r < 20 else "FAIL"))
        Tub_s = f"{Tub:.4f}" if Tub > 0 else "~0"
        print(f"  {row[0]:<10s} {Tub_s:>7s} {T0:8.2f} {T_new:9.4f} "
              f"{verdict:>8s}")

    # Correlation
    Tobs_v = [m[1] for m in ps34.MATERIALS if m[0] not in ps34.OFF_MODEL]
    Tp_fit = [compute_SC_tpred_rescaled(m, A, p)
              for m in ps34.MATERIALS if m[0] not in ps34.OFF_MODEL]
    r_fit, p_rfit = pearsonr(np.log10(Tobs_v), np.log10(Tp_fit))

    print(f"\n  (6) Pearsonr log-log po rescaled fit:")
    print(f"    baseline r = 0.9219,  po rescaled P7.10: r = {r_fit:.4f}")

    # Verdict
    print("\n" + "=" * 78)
    print("  VERDICT ps35")
    print("=" * 78)
    trade_ok = (rms_fit < rms0 + 0.10) and (fails_fit < fails0 - 5)
    print(f"""
    Joint (A, p) fit z regularizacja:
       A     = {A:.3f}
       p     = {p:.3f}
       RMS_SC = {rms_fit:.4f}   (baseline: {rms0:.4f}, delta = {rms_fit-rms0:+.4f})
       r_SC   = {r_fit:.4f}     (baseline: 0.9219)
       FAIL   = {fails_fit}/{len(ps34.NON_SC)}  (baseline: {fails0}/{len(ps34.NON_SC)})

    Pure-SC fit (max RMS fit):
       A = {A_pure:.3f}, p = {p_pure:.3f}
       RMS_SC = {rms_pure:.4f},  FAIL = {fails_pure}

    Status: {'POWODZENIE' if trade_ok else 'POTRZEBNY alternatywny approach'}

    Interpretacja:
      - A kompensuje srednia g-reduction; fizycznie rescaluje C_0 Eq. 5
      - Efektywny C_0^new = C_0 * A = {48.8222 * A:.2f} (vs stare 48.82)
      - To oznacza ze C_0 w paperze v1 byl zlozony z prawdziwy-C_0 / <g>
      - Po rozdzieleniu: model daje osobno (geometry) x (N_F-factor)

    Jesli p ~ 1 i A ~ 1/mediana(g_SC) — to czysta, rzetelna refaktoryzacja.
    Jesli p >> 1 lub A < 1 — to forma g(x) = x^p jest zbyt agresywna, trzeba
    sigmoid lub threshold.
""")


if __name__ == "__main__":
    main()
