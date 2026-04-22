#!/usr/bin/env python3
# =============================================================================
#  ps37_P712_sigmoid.py
# -----------------------------------------------------------------------------
#  P7.12: sigmoid cutoff dla N_F factor. Cel: zamknac 7 pozostalych FAILow
#  (Cu, Ag, Au, Li, Ca, Mg, Graphite) bez psucia SC CORE (zwl. Hg_elem,
#  Os, W przy N_F ~ 0.3).
#
#  Testowane formy g(x), x = N_F / N_F_ref:
#    S-1: Sigmoid  g(x) = s(k, x0)(x) / s(k, x0)(1)
#    S-2: Threshold piecewise g(x) = max(0, (x-x_t)/(1-x_t))^p
#    S-3: Stretched-exp  g(x) = exp(-alpha * (1/x - 1)^beta) dla x <= 1
#    S-4: Hill  g(x) = x^n / (x^n + x_h^n)  (Michaelis-Menten pod lewej)
#
#  Fit (params) na joint obj (RMS_SC + 0.3*median_NonSC + 0.05*#FAIL).
#  Cel: FAIL <= 4 bez dRMS > +0.02.
# =============================================================================
import numpy as np
from scipy.optimize import minimize
from scipy.stats import pearsonr

import importlib.util, os, sys
spec = importlib.util.spec_from_file_location(
    "ps34", os.path.join(os.path.dirname(__file__), "ps34_P710_integrated.py"))
ps34 = importlib.util.module_from_spec(spec); sys.modules["ps34"] = ps34
spec.loader.exec_module(ps34)

# Importuj tez ps36 dla P7.11 Stoner override
spec36 = importlib.util.spec_from_file_location(
    "ps36", os.path.join(os.path.dirname(__file__), "ps36_P711_magnetic.py"))
ps36 = importlib.util.module_from_spec(spec36); sys.modules["ps36"] = ps36
spec36.loader.exec_module(ps36)


# -----------------------------------------------------------------------------
# Formy g(x)
# -----------------------------------------------------------------------------
def g_sigmoid(x, k, x0):
    """S-1: S-shaped g(x) = sigmoid((x-x0)*k), normalized g(1)=1"""
    x = max(x, 1e-6)
    raw = 1.0 / (1.0 + np.exp(-k * (x - x0)))
    norm = 1.0 / (1.0 + np.exp(-k * (1.0 - x0)))
    return min(1.0, raw / norm)


def g_threshold(x, x_t, p):
    """S-2: Zero poniżej x_t, power-law powyzej"""
    if x < x_t: return 0.0
    y = (x - x_t) / (1.0 - x_t) if (1 - x_t) > 1e-6 else 0.0
    return min(1.0, max(0.0, y) ** p)


def g_stretched(x, alpha, beta):
    """S-3: exp(-alpha * (1/x - 1)^beta) dla x <= 1"""
    if x <= 1e-6: return 0.0
    if x >= 1.0: return 1.0
    return float(np.exp(-alpha * (1.0/x - 1.0) ** beta))


def g_hill(x, x_h, n):
    """S-4: Hill function x^n / (x^n + x_h^n)"""
    if x <= 0: return 0.0
    return x**n / (x**n + x_h**n)


# -----------------------------------------------------------------------------
# Wrappery: uzywaja compute_*_with_P711 z ps36 (P7.5 + P7.11), potem nakladamy
# g(x) dla P7.12. Fe/Co/Ni: jawny FM-override (jak w ps36.nonsc_stats_P711).
# -----------------------------------------------------------------------------
KAPPA = 0.247
MU = 5.0
FM_SET = {"Fe", "Co", "Ni"}


def compute_SC_with_P712(mat, g_func, g_params):
    """P7.5 + P7.11 (Stoner) + nasze g_func(x) dla P7.12."""
    T_p711 = ps36.compute_SC_with_P711(mat, KAPPA, MU, p_P710=None)
    if T_p711 <= 0:
        return 0.0
    x = mat[6] / ps34.N_F_REF
    g = g_func(x, *g_params)
    return T_p711 * g


def compute_nonSC_with_P712(row, g_func, g_params):
    """P7.5 + P7.11 + P7.12 z jawnym FM-override dla Fe/Co/Ni."""
    name = row[0]
    if name in FM_SET:
        T_base = ps34.compute_nonsc_tpred(row, None)
        T_p711 = ps36.apply_B_P711(T_base, 100.0)  # FM -> B_mag ~ 0
    else:
        T_p711 = ps36.compute_nonsc_with_P711(row, KAPPA, MU, p_P710=None)
    if T_p711 < 0:
        return 0.0
    x = row[7] / ps34.N_F_REF
    g = g_func(x, *g_params)
    return T_p711 * g


# -----------------------------------------------------------------------------
# Stats
# -----------------------------------------------------------------------------
def sc_rms(g_func, params):
    resid = []
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        T = compute_SC_with_P712(m, g_func, params)
        if T <= 0: return 1e9
        resid.append(np.log10(T) - np.log10(m[1]))
    return float(np.sqrt(np.mean(np.asarray(resid)**2)))


def nonsc_stats(g_func, params):
    fails, borders = 0, 0; Tp = []
    for row in ps34.NON_SC:
        Tub = row[1]
        T = compute_nonSC_with_P712(row, g_func, params)
        Tp.append(T)
        if Tub < 1e-4:
            if T > 1.0: fails += 1
            elif T > 0.1: borders += 1
        else:
            r = T / Tub if Tub > 1e-6 else 1e9
            if r > 20: fails += 1
            elif r > 5: borders += 1
    return fails, borders, float(np.median(Tp)), float(np.percentile(Tp, 95))


def joint_obj(g_func, params, lam_fail=0.05, lam_med=0.3):
    rms = sc_rms(g_func, params)
    fails, _, med, _ = nonsc_stats(g_func, params)
    return rms + lam_med * med + lam_fail * fails


# -----------------------------------------------------------------------------
# Fit-ery
# -----------------------------------------------------------------------------
def fit_form(g_func, name, param_bounds, param_starts):
    """Fit (params) dla danego g_func."""
    best = (1e9, None)
    for start in param_starts:
        try:
            res = minimize(
                lambda p: (
                    1e9 if any(
                        (p[i] < param_bounds[i][0] or p[i] > param_bounds[i][1])
                        for i in range(len(p))
                    )
                    else joint_obj(g_func, p)
                ),
                start, method="Nelder-Mead",
                options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 300}
            )
            if res.fun < best[0]:
                best = (float(res.fun), list(res.x))
        except Exception as e:
            continue
    return best[1], best[0]


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("  ps37_P712_sigmoid.py - sigmoid/threshold/hill g(N_F) formy")
    print("=" * 78)

    # Baseline: P7.5 only
    rms_0 = ps34.sc_rms(None)[0]
    fails_0, _, med_0, _ = ps34.nonsc_fails(None)
    print("\n  Baseline (P7.5 only):")
    print("    SC RMS =", round(rms_0, 4),
          "  non-SC FAIL =", fails_0, "/", len(ps34.NON_SC))

    # Reference point: power-law p=0.5 + P7.11 (from ps36)
    def g_ref(x, p=0.5): return min(1.0, x ** p) if x > 0 else 0.0
    rms_ref = sc_rms(g_ref, (0.5,))
    fails_ref, _, med_ref, _ = nonsc_stats(g_ref, (0.5,))
    print("\n  Ref (P7.10 p=0.5 + P7.11 kappa=0.247):")
    print("    SC RMS =", round(rms_ref, 4),
          "  non-SC FAIL =", fails_ref, "/", len(ps34.NON_SC),
          "  median =", round(med_ref, 3))

    # -- S-1: Sigmoid
    print("\n" + "-" * 78)
    print("  S-1: Sigmoid  g(x) = [1+exp(-k(x-x0))]^-1 / norm")
    print("-" * 78)
    params, fval = fit_form(
        g_sigmoid, "S-1",
        param_bounds=[(1.0, 50.0), (0.0, 0.9)],
        param_starts=[[5.0, 0.3], [10.0, 0.3], [15.0, 0.2], [20.0, 0.25]],
    )
    k_opt, x0_opt = params
    rms_s1 = sc_rms(g_sigmoid, params)
    fails_s1, bord_s1, med_s1, p95_s1 = nonsc_stats(g_sigmoid, params)
    print("    Fitted: k =", round(k_opt, 3), ", x0 =", round(x0_opt, 3))
    print("    SC RMS =", round(rms_s1, 4),
          "  FAIL =", fails_s1, "/", len(ps34.NON_SC),
          "  BORDER =", bord_s1)
    print("    median =", round(med_s1, 3), "K  P95 =", round(p95_s1, 2), "K")

    # -- S-2: Threshold
    print("\n" + "-" * 78)
    print("  S-2: Threshold  g(x) = max(0, (x-x_t)/(1-x_t))^p")
    print("-" * 78)
    params, fval = fit_form(
        g_threshold, "S-2",
        param_bounds=[(0.01, 0.5), (0.1, 5.0)],
        param_starts=[[0.1, 1.0], [0.15, 0.5], [0.2, 1.0], [0.25, 1.5]],
    )
    xt_opt, p_opt = params
    rms_s2 = sc_rms(g_threshold, params)
    fails_s2, bord_s2, med_s2, p95_s2 = nonsc_stats(g_threshold, params)
    print("    Fitted: x_t =", round(xt_opt, 3), ", p =", round(p_opt, 3))
    print("    SC RMS =", round(rms_s2, 4),
          "  FAIL =", fails_s2, "/", len(ps34.NON_SC),
          "  BORDER =", bord_s2)
    print("    median =", round(med_s2, 3), "K  P95 =", round(p95_s2, 2), "K")

    # -- S-3: Stretched-exp
    print("\n" + "-" * 78)
    print("  S-3: Stretched-exp  g(x) = exp(-alpha (1/x - 1)^beta)")
    print("-" * 78)
    params, fval = fit_form(
        g_stretched, "S-3",
        param_bounds=[(0.1, 20.0), (0.3, 4.0)],
        param_starts=[[1.0, 1.0], [2.0, 1.5], [5.0, 0.7], [3.0, 2.0]],
    )
    a_opt, b_opt = params
    rms_s3 = sc_rms(g_stretched, params)
    fails_s3, bord_s3, med_s3, p95_s3 = nonsc_stats(g_stretched, params)
    print("    Fitted: alpha =", round(a_opt, 3), ", beta =", round(b_opt, 3))
    print("    SC RMS =", round(rms_s3, 4),
          "  FAIL =", fails_s3, "/", len(ps34.NON_SC),
          "  BORDER =", bord_s3)
    print("    median =", round(med_s3, 3), "K  P95 =", round(p95_s3, 2), "K")

    # -- S-4: Hill
    print("\n" + "-" * 78)
    print("  S-4: Hill  g(x) = x^n / (x^n + x_h^n)")
    print("-" * 78)
    params, fval = fit_form(
        g_hill, "S-4",
        param_bounds=[(0.05, 0.9), (1.0, 20.0)],
        param_starts=[[0.2, 3.0], [0.25, 5.0], [0.15, 4.0], [0.3, 8.0]],
    )
    xh_opt, n_opt = params
    rms_s4 = sc_rms(g_hill, params)
    fails_s4, bord_s4, med_s4, p95_s4 = nonsc_stats(g_hill, params)
    print("    Fitted: x_h =", round(xh_opt, 3), ", n =", round(n_opt, 3))
    print("    SC RMS =", round(rms_s4, 4),
          "  FAIL =", fails_s4, "/", len(ps34.NON_SC),
          "  BORDER =", bord_s4)
    print("    median =", round(med_s4, 3), "K  P95 =", round(p95_s4, 2), "K")

    # -- Wybor najlepszego
    print("\n" + "=" * 78)
    print("  Porownanie form (wszystkie obejmuja P7.11 Stoner)")
    print("=" * 78)
    results = [
        ("REF power p=0.5", rms_ref, fails_ref, med_ref),
        ("S-1 sigmoid", rms_s1, fails_s1, med_s1),
        ("S-2 threshold", rms_s2, fails_s2, med_s2),
        ("S-3 stretched", rms_s3, fails_s3, med_s3),
        ("S-4 hill", rms_s4, fails_s4, med_s4),
    ]
    print("  {:<18s} {:>8s} {:>8s} {:>10s}".format(
        "Form", "SC RMS", "FAIL", "medNonSC"))
    for name, rms, fails, med in results:
        print("  {:<18s} {:>8.4f} {:>8d} {:>10.3f}".format(
            name, rms, fails, med))

    # Find best by FAIL then RMS
    best_form = min(results[1:], key=lambda r: (r[2], r[1]))
    print("\n  Best: {}  (FAIL={}, RMS={:.4f})".format(
        best_form[0], best_form[2], best_form[1]))

    # Per-material dla najlepszego wariantu
    name2fn = {
        "S-1 sigmoid":   (g_sigmoid, [k_opt, x0_opt]),
        "S-2 threshold": (g_threshold, [xt_opt, p_opt]),
        "S-3 stretched": (g_stretched, [a_opt, b_opt]),
        "S-4 hill":      (g_hill, [xh_opt, n_opt]),
    }
    g_best, pars_best = name2fn[best_form[0]]

    print("\n  Per-material, best form on SC CORE:")
    print("  {:<12s} {:>7s} {:>8s} {:>8s} {:>7s}".format(
        "name", "T_obs", "T_ref", "T_p712", "d_p712"))
    for m in ps34.MATERIALS:
        if m[0] in ps34.OFF_MODEL: continue
        T_ref = compute_SC_with_P712(m, g_ref, (0.5,))
        T_new = compute_SC_with_P712(m, g_best, pars_best)
        d_new = np.log10(max(T_new, 1e-6)) - np.log10(m[1])
        print("  {:<12s} {:>7.2f} {:>8.2f} {:>8.2f} {:>+7.3f}".format(
            m[0], m[1], T_ref, T_new, d_new))

    print("\n  Per-material, best form on non-SC:")
    print("  {:<12s} {:>7s} {:>8s} {:>9s} {:>8s}".format(
        "name", "T_ub", "T_ref", "T_p712", "verdict"))
    for row in ps34.NON_SC:
        T_ref = compute_nonSC_with_P712(row, g_ref, (0.5,))
        T_new = compute_nonSC_with_P712(row, g_best, pars_best)
        Tub = row[1]
        if Tub < 1e-4:
            if T_new < 0.1: verdict = "PASS"
            elif T_new < 1.0: verdict = "BORDER"
            else: verdict = "FAIL"
        else:
            r = T_new / Tub if Tub > 0 else 1e9
            if r < 5: verdict = "PASS"
            elif r < 20: verdict = "BORDER"
            else: verdict = "FAIL"
        Tub_s = "{:.4f}".format(Tub) if Tub > 0 else "~0"
        print("  {:<12s} {:>7s} {:>8.3f} {:>9.4f} {:>8s}".format(
            row[0], Tub_s, T_ref, T_new, verdict))

    # Verdict
    print("\n" + "=" * 78)
    print("  VERDICT ps37 / P7.12")
    print("=" * 78)
    d_rms = best_form[1] - rms_0
    red = fails_0 - best_form[2]
    print("    Najlepsza forma: ", best_form[0])
    print("    RMS_SC = {:.4f} (baseline {:.4f}, dRMS = {:+.4f})".format(
        best_form[1], rms_0, d_rms))
    print("    non-SC FAIL = {}/{}  (baseline {})".format(
        best_form[2], len(ps34.NON_SC), fails_0))
    print("    Redukcja FAILow od P7.5 only: {}".format(red))
    status = "AKCEPTOWANE" if (d_rms < 0.02 and red >= 15) else (
        "CZESCIOWE" if red >= 10 else "NIEWYSTARCZAJACE")
    print("    Status: {}".format(status))
    print()
    print("    Progressive reduction FAIL:")
    print("      P7.5:                 20")
    print("      + P7.10 (p=0.5):      12  (-8)")
    print("      + P7.11 (Stoner):      7  (-5)")
    print("      + P7.12 ({:<12s}): {:>3d}  ({:+d})".format(
        best_form[0][:12], best_form[2], best_form[2] - fails_ref))


if __name__ == "__main__":
    main()
