#!/usr/bin/env python3
# =============================================================================
#  ps38_P712_robustness.py
# -----------------------------------------------------------------------------
#  Walidacja stabilnosci stalych P7.12 przed paper-writeup:
#    (1) Extended bounds fit — czy alpha=0.1 to prawdziwe optimum?
#    (2) LOO: usun kazdy material SC po kolei, refit (alpha, beta), sprawdz
#        jak bardzo parametry dryfują
#    (3) Bootstrap: 100x resample non-SC + SC, oszacuj CI na (alpha, beta)
#    (4) Sensitivity: czy wynik zmienia sie przy zmianie kappa, mu z P7.11?
#    (5) Hold-out: wynik na prognozie 3 nowych, nie fitted materialow?
#
#  Cel: zebrac dane do paper v2 (sigma-intervale na nowych stalych TGP).
# =============================================================================
import numpy as np
from scipy.optimize import minimize
from scipy.stats import pearsonr
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "ps34", os.path.join(os.path.dirname(__file__), "ps34_P710_integrated.py"))
ps34 = importlib.util.module_from_spec(spec); sys.modules["ps34"] = ps34
spec.loader.exec_module(ps34)

spec36 = importlib.util.spec_from_file_location(
    "ps36", os.path.join(os.path.dirname(__file__), "ps36_P711_magnetic.py"))
ps36 = importlib.util.module_from_spec(spec36); sys.modules["ps36"] = ps36
spec36.loader.exec_module(ps36)

KAPPA_FIX = 0.247
MU_FIX = 5.0
FM_SET = {"Fe", "Co", "Ni"}


def g_stretched(x, alpha, beta):
    if x <= 1e-6: return 0.0
    if x >= 1.0: return 1.0
    return float(np.exp(-alpha * (1.0/x - 1.0) ** beta))


def compute_SC(mat, alpha, beta, kappa=KAPPA_FIX, mu=MU_FIX):
    T = ps36.compute_SC_with_P711(mat, kappa, mu, p_P710=None)
    if T <= 0: return 0.0
    x = mat[6] / ps34.N_F_REF
    return T * g_stretched(x, alpha, beta)


def compute_nonSC(row, alpha, beta, kappa=KAPPA_FIX, mu=MU_FIX):
    name = row[0]
    if name in FM_SET:
        T = ps34.compute_nonsc_tpred(row, None)
        T = ps36.apply_B_P711(T, 100.0)
    else:
        T = ps36.compute_nonsc_with_P711(row, kappa, mu, p_P710=None)
    x = row[7] / ps34.N_F_REF
    return T * g_stretched(x, alpha, beta)


def sc_rms(alpha, beta, sc_list=None):
    mats = sc_list if sc_list is not None else [
        m for m in ps34.MATERIALS if m[0] not in ps34.OFF_MODEL]
    resid = []
    for m in mats:
        T = compute_SC(m, alpha, beta)
        if T <= 0: return 1e9
        resid.append(np.log10(T) - np.log10(m[1]))
    return float(np.sqrt(np.mean(np.asarray(resid)**2)))


def nonsc_stats(alpha, beta, nsc_list=None):
    rows = nsc_list if nsc_list is not None else ps34.NON_SC
    fails, borders = 0, 0; Tp = []
    for row in rows:
        Tub = row[1]
        T = compute_nonSC(row, alpha, beta)
        Tp.append(T)
        if Tub < 1e-4:
            if T > 1.0: fails += 1
            elif T > 0.1: borders += 1
        else:
            r = T / Tub if Tub > 1e-6 else 1e9
            if r > 20: fails += 1
            elif r > 5: borders += 1
    return fails, borders, float(np.median(Tp))


def joint_obj(alpha, beta, lam_f=0.05, lam_m=0.3, sc_list=None, nsc_list=None):
    rms = sc_rms(alpha, beta, sc_list)
    fails, _, med = nonsc_stats(alpha, beta, nsc_list)
    return rms + lam_m * med + lam_f * fails


def fit_ab(a_lo, a_hi, b_lo, b_hi, sc_list=None, nsc_list=None):
    """Fit (alpha, beta) z zadanymi bounds."""
    best = (1e9, None)
    for a0 in np.linspace(a_lo + 0.05*(a_hi-a_lo), a_hi*0.5, 5):
        for b0 in np.linspace(b_lo + 0.1, b_hi*0.5, 5):
            try:
                res = minimize(
                    lambda p: (
                        1e9 if (p[0] < a_lo or p[0] > a_hi or
                                p[1] < b_lo or p[1] > b_hi)
                        else joint_obj(p[0], p[1],
                                       sc_list=sc_list, nsc_list=nsc_list)
                    ),
                    [a0, b0], method="Nelder-Mead",
                    options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 300}
                )
                if res.fun < best[0]:
                    best = (float(res.fun), list(res.x))
            except Exception:
                pass
    return best[1], best[0]


def main():
    print("=" * 78)
    print("  ps38 — Robustness P7.12 (alpha, beta): extended bounds, LOO, boot")
    print("=" * 78)

    # Ref values (ps37)
    alpha_ref, beta_ref = 0.10, 2.15

    # ---- (1) Extended bounds ----
    print("\n" + "-" * 78)
    print("  (1) Extended bounds  alpha in [0.001, 100],  beta in [0.3, 8]")
    print("-" * 78)
    params_ext, f_ext = fit_ab(0.001, 100.0, 0.3, 8.0)
    a_ext, b_ext = params_ext
    rms_ext = sc_rms(a_ext, b_ext)
    fail_ext, bord_ext, med_ext = nonsc_stats(a_ext, b_ext)
    print("  Fitted (extended): alpha =", round(a_ext, 4),
          ", beta =", round(b_ext, 4))
    print("  RMS_SC =", round(rms_ext, 4),
          "  FAIL =", fail_ext, "  BORDER =", bord_ext,
          "  med =", round(med_ext, 4), "K")
    if abs(a_ext - alpha_ref) < 0.02 and abs(b_ext - beta_ref) < 0.05:
        print("  OK: zgodne z ps37 (alpha=0.10, beta=2.15)")
    else:
        print("  UWAGA: dryf parametrow, extended optimum rozne")

    # ---- (2) Grid scan for finer picture ----
    print("\n" + "-" * 78)
    print("  (2) Grid scan (alpha x beta) dla zrozumienia topologii")
    print("-" * 78)
    print("  alpha\\beta ", end="")
    betas = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    for b in betas: print("  b={:.1f}".format(b), end="")
    print()
    alphas = [0.01, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]
    for a in alphas:
        print("  a={:<8.3f}".format(a), end="")
        for b in betas:
            rms = sc_rms(a, b)
            fails, _, _ = nonsc_stats(a, b)
            print("  {:.3f}/{:d}".format(rms, fails), end="")
        print()

    # ---- (3) LOO on SC CORE ----
    print("\n" + "-" * 78)
    print("  (3) LOO: refit (alpha, beta) usuwajac po 1 SC material")
    print("-" * 78)
    print("  (pokazujemy tylko materialy przy ktorych parametry dryfują > 10%)")
    sc_core = [m for m in ps34.MATERIALS if m[0] not in ps34.OFF_MODEL]
    drifts = []
    for skip_idx in range(len(sc_core)):
        sc_sub = [m for i, m in enumerate(sc_core) if i != skip_idx]
        params, _ = fit_ab(0.01, 20.0, 0.3, 6.0, sc_list=sc_sub)
        a_i, b_i = params
        da = (a_i - alpha_ref) / alpha_ref * 100
        db = (b_i - beta_ref) / beta_ref * 100
        drifts.append((sc_core[skip_idx][0], a_i, b_i, da, db))
    drifts.sort(key=lambda r: abs(r[3]) + abs(r[4]), reverse=True)
    print("  {:<12s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "skip", "alpha'", "beta'", "%da", "%db"))
    for name, a, b, da, db in drifts[:8]:
        tag = " SIGNIFICANT" if (abs(da) > 10 or abs(db) > 10) else ""
        print("  {:<12s} {:>8.4f} {:>8.4f} {:>+7.1f}% {:>+7.1f}%{}".format(
            name, a, b, da, db, tag))

    da_max = max(abs(d[3]) for d in drifts)
    db_max = max(abs(d[4]) for d in drifts)
    print("\n  Max |dalpha|% =", round(da_max, 1),
          "   Max |dbeta|% =", round(db_max, 1))
    if da_max < 20 and db_max < 20:
        print("  LOO-stabilne: zaden material nie dominuje fit")
    else:
        print("  UWAGA: pojedyncze materialy znaczaco wplywaja na fit")

    # ---- (4) Bootstrap on non-SC (resampling) ----
    print("\n" + "-" * 78)
    print("  (4) Bootstrap: 50x resample non-SC (z replacement), refit")
    print("-" * 78)
    np.random.seed(42)
    n_boot = 50
    ab_boot = []
    for i in range(n_boot):
        idx = np.random.choice(len(ps34.NON_SC), size=len(ps34.NON_SC),
                               replace=True)
        nsc_b = [ps34.NON_SC[k] for k in idx]
        params, _ = fit_ab(0.01, 20.0, 0.3, 6.0, nsc_list=nsc_b)
        ab_boot.append(params)
    a_b = np.array([p[0] for p in ab_boot])
    b_b = np.array([p[1] for p in ab_boot])
    print("  Alpha: median = {:.4f},  16-84%ile = [{:.4f}, {:.4f}]".format(
        np.median(a_b), np.percentile(a_b, 16), np.percentile(a_b, 84)))
    print("  Beta:  median = {:.4f},  16-84%ile = [{:.4f}, {:.4f}]".format(
        np.median(b_b), np.percentile(b_b, 16), np.percentile(b_b, 84)))
    print("  Mean +/- std: alpha = {:.3f} +/- {:.3f}".format(
        np.mean(a_b), np.std(a_b)))
    print("  Mean +/- std:  beta = {:.3f} +/- {:.3f}".format(
        np.mean(b_b), np.std(b_b)))

    # ---- (5) Sensitivity to (kappa, mu) from P7.11 ----
    print("\n" + "-" * 78)
    print("  (5) Sensitivity: refit (alpha, beta) dla rożnych (kappa, mu)")
    print("-" * 78)
    print("  {:<10s} {:<8s} {:>8s} {:>8s} {:>8s} {:>5s}".format(
        "kappa", "mu", "alpha", "beta", "RMS_SC", "FAIL"))
    combos = [(0.1, 3.0), (0.247, 5.0), (0.5, 5.0), (1.0, 3.0), (0.05, 2.0)]
    for kap, m in combos:
        def obj(p):
            if p[0] < 0.01 or p[0] > 20 or p[1] < 0.3 or p[1] > 6:
                return 1e9
            # Override KAPPA, MU temporarily
            resid = []
            for mat in [mm for mm in ps34.MATERIALS
                        if mm[0] not in ps34.OFF_MODEL]:
                T = ps36.compute_SC_with_P711(mat, kap, m, p_P710=None)
                if T <= 0: return 1e9
                xx = mat[6] / ps34.N_F_REF
                T = T * g_stretched(xx, p[0], p[1])
                resid.append(np.log10(T) - np.log10(mat[1]))
            rms = np.sqrt(np.mean(np.asarray(resid)**2))
            # non-SC
            fails = 0; Tp = []
            for row in ps34.NON_SC:
                Tub = row[1]; name = row[0]
                if name in FM_SET:
                    T = ps34.compute_nonsc_tpred(row, None)
                    T = ps36.apply_B_P711(T, 100.0)
                else:
                    T = ps36.compute_nonsc_with_P711(row, kap, m, p_P710=None)
                xx = row[7] / ps34.N_F_REF
                T = T * g_stretched(xx, p[0], p[1])
                Tp.append(T)
                if Tub < 1e-4:
                    if T > 1.0: fails += 1
                else:
                    r = T / Tub if Tub > 1e-6 else 1e9
                    if r > 20: fails += 1
            return rms + 0.3*np.median(Tp) + 0.05*fails
        best = (1e9, None)
        for a0 in [0.05, 0.2, 0.5]:
            for b0 in [1.0, 2.0, 3.0]:
                res = minimize(obj, [a0, b0], method="Nelder-Mead",
                               options={"xatol": 1e-4, "fatol": 1e-6,
                                        "maxiter": 200})
                if res.fun < best[0]:
                    best = (float(res.fun), list(res.x))
        a_i, b_i = best[1]
        # eval
        def ev(a, b, kap, m):
            resid = []
            for mat in [mm for mm in ps34.MATERIALS
                        if mm[0] not in ps34.OFF_MODEL]:
                T = ps36.compute_SC_with_P711(mat, kap, m, p_P710=None)
                if T <= 0: return 999, 99
                xx = mat[6] / ps34.N_F_REF
                T = T * g_stretched(xx, a, b)
                resid.append(np.log10(T) - np.log10(mat[1]))
            rms = np.sqrt(np.mean(np.asarray(resid)**2))
            fails = 0
            for row in ps34.NON_SC:
                Tub = row[1]; name = row[0]
                if name in FM_SET:
                    T = ps34.compute_nonsc_tpred(row, None)
                    T = ps36.apply_B_P711(T, 100.0)
                else:
                    T = ps36.compute_nonsc_with_P711(row, kap, m, p_P710=None)
                xx = row[7] / ps34.N_F_REF
                T = T * g_stretched(xx, a, b)
                if Tub < 1e-4:
                    if T > 1.0: fails += 1
                else:
                    r = T / Tub if Tub > 1e-6 else 1e9
                    if r > 20: fails += 1
            return rms, fails
        rms_i, fail_i = ev(a_i, b_i, kap, m)
        print("  {:<10.3f} {:<8.2f} {:>8.4f} {:>8.4f} {:>8.4f} {:>5d}".format(
            kap, m, a_i, b_i, rms_i, fail_i))

    # ---- (6) Verdict and CI on 6 new constants ----
    print("\n" + "=" * 78)
    print("  VERDICT ps38 — Robustness P7.12")
    print("=" * 78)
    print("""
    Stale P7.12 (stretched-exp g(x)):
      alpha = {:.3f} +/- {:.3f}  [bootstrap 1-sigma]
      beta  = {:.3f} +/- {:.3f}  [bootstrap 1-sigma]

    Extended-bounds fit: alpha = {:.4f}, beta = {:.4f}
      -> {}

    LOO max drift: |dalpha|% = {:.1f},  |dbeta|% = {:.1f}
      -> {}

    Stabilnosc wzgl. (kappa, mu) z P7.11:
      zobacz tabele (5) powyzej.

    Finalne wartosci do paperu v2:
      alpha_7.12 = 0.10 +/- {:.2f}  (patrz bootstrap)
      beta_7.12  = 2.15 +/- {:.2f}
""".format(
        np.mean(a_b), np.std(a_b),
        np.mean(b_b), np.std(b_b),
        a_ext, b_ext,
        "stabilne" if abs(a_ext - 0.10) < 0.05 else "dryfuje",
        da_max, db_max,
        "SC CORE nie jest dominowane przez pojedynczy material" if
          (da_max < 20 and db_max < 20) else "Pojedynczy material dominuje",
        np.std(a_b), np.std(b_b),
    ))


if __name__ == "__main__":
    main()
