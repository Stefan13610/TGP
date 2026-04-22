#!/usr/bin/env python3
# =============================================================================
#  r13b_loo_sweep.py
# -----------------------------------------------------------------------------
#  r13 pokazalo: U_full wygrywa AIC ale LOO overfit 2.62x. Potrzebujemy
#  znalezc najbardziej ROBUSTNY model (najlepszy LOO), nie najlepszy fit.
#
#  Wykonujemy LOO CV na wszystkich 6 modelach U0..U_full + porownujemy
#  LOO profile RMS.
# =============================================================================
import numpy as np
import importlib.util, os, sys

for mod in ["r00", "r01", "r02", "r06", "r10", "r12", "r13_unified_all_classes"]:
    fn = os.path.join(os.path.dirname(__file__), mod + ".py")
    if not os.path.exists(fn): continue
    spec = importlib.util.spec_from_file_location(mod, fn)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod] = m
    spec.loader.exec_module(m)

import r00, r01, r02, r06, r10, r12
r13 = sys.modules["r13_unified_all_classes"]


def fit_ols(X, y):
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coefs
    res = y - pred
    n, k = X.shape
    ss = np.sum(res**2)
    rms = np.sqrt(ss / n)
    aic = n * np.log(ss/n + 1e-30) + 2 * k
    bic = n * np.log(ss/n + 1e-30) + k * np.log(n)
    return coefs, pred, rms, aic, bic


def build_models(logTh, logNF, v_arr, delta_s, delta_sp):
    """Return dict of X matrices."""
    X_U0 = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                            logTh, logNF])
    X_U1 = np.column_stack([X_U0, v_arr])
    X_U2 = np.column_stack([X_U1, v_arr*delta_s, v_arr*delta_sp])
    X_U3 = np.column_stack([X_U1, delta_s*logTh, delta_sp*logTh])
    X_U4 = np.column_stack([X_U1, delta_s*logNF, delta_sp*logNF])
    X_F = np.column_stack([X_U1,
                            delta_s*logTh, delta_sp*logTh,
                            delta_s*logNF, delta_sp*logNF,
                            v_arr*delta_s, v_arr*delta_sp])
    # Plus intermediate:
    # U24 = U1 + v*cls + logNF*cls (U2 union U4)
    X_U24 = np.column_stack([X_U1, v_arr*delta_s, v_arr*delta_sp,
                              delta_s*logNF, delta_sp*logNF])
    return {"U0": X_U0, "U1": X_U1, "U2": X_U2, "U3": X_U3, "U4": X_U4,
            "U24": X_U24, "U_full": X_F}


def main():
    print("=" * 92)
    print("  r13b_loo_sweep.py - LOO robustness across U0..U_full")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    data = r12.build_dataset_N37()

    def get_class(name):
        return r02.ORB_CLASS.get(name, "?")

    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = get_class(d["name"])
        if d["name"] in NOBLE or d["name"] in ALKALINE: cls = "s"
        v = r13.V_COUNT.get(d["name"])
        if v is None: continue
        rows.append({
            "name": d["name"], "cls": cls, "R": R, "v": v,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
        })

    N = len(rows)
    print("  N = {}".format(N))

    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    cls_arr = np.array([r["cls"] for r in rows])
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)

    Xs = build_models(logTh, logNF, v_arr, delta_s, delta_sp)

    # -------- (A) AIC + LOO RMS + LOO profile RMS --------
    print("\n" + "-" * 92)
    print("  (A) Model quality: AIC, full RMS, LOO RMS, LOO profile RMS")
    print("-" * 92)
    print("  {:<8s}  {:>3s}  {:>8s}  {:>8s}  {:>8s}  {:>9s}  {:>8s}  {:>8s}".format(
        "model", "k", "RMS", "AIC", "LOO_RMS", "overfit", "prof_m", "prof_max"))

    results = {}
    for key in ["U0", "U1", "U2", "U3", "U4", "U24", "U_full"]:
        X = Xs[key]
        k = X.shape[1]
        # Full fit
        coefs, pred_full, rms_full, aic, bic = fit_ols(X, logR)
        # LOO CV
        loo_logR = np.zeros(N)
        for i in range(N):
            mask = np.ones(N, bool); mask[i] = False
            try:
                coefs_i, _, _, _ = np.linalg.lstsq(X[mask], logR[mask], rcond=None)
                loo_logR[i] = X[i] @ coefs_i
            except Exception:
                loo_logR[i] = np.nan
        loo_dlog = loo_logR - logR
        loo_rms = np.sqrt(np.nanmean(loo_dlog**2))
        overfit = loo_rms / max(rms_full, 1e-3)

        # LOO profile RMS (use LOO predictions of R in BG formula)
        prof_rms = np.zeros(N)
        for i, r in enumerate(rows):
            R_pred = 10**loo_logR[i]
            Th = r["Theta_D"]; rho0 = r["rho_0"]
            pts = []
            for T in r00.T_POINTS:
                rho_p = r01.rho_bg(T, rho0, R_pred, Th)
                rho_o = r["rho_obs"][T]
                pts.append(np.log10(rho_p) - np.log10(rho_o))
            prof_rms[i] = np.sqrt(np.mean(np.array(pts)**2))

        results[key] = {
            "k": k, "coefs": coefs, "rms_full": rms_full,
            "aic": aic, "loo_rms": loo_rms, "overfit": overfit,
            "prof_rms": prof_rms, "loo_logR": loo_logR,
        }
        print("  {:<8s}  {:>3d}  {:>8.4f}  {:>8.3f}  {:>8.4f}  {:>8.2f}x  {:>8.4f}  {:>8.4f}".format(
            key, k, rms_full, aic, loo_rms, overfit,
            np.mean(prof_rms), np.max(prof_rms)))

    # -------- (B) Pick robust model --------
    # Criterium: min mean LOO-profile-RMS (s nepoval overfit ratio <3)
    print("\n" + "-" * 92)
    print("  (B) Robust model selection (kryteria: min LOO-profile-RMS)")
    print("-" * 92)
    keys = list(results.keys())
    sorted_keys = sorted(keys, key=lambda k: np.mean(results[k]["prof_rms"]))
    best_robust = sorted_keys[0]
    print("\n  Ranking po mean LOO profile RMS:")
    for k in sorted_keys:
        r = results[k]
        mark = " <-- ROBUST" if k == best_robust else ""
        print("    {:<8s}  mean prof RMS = {:.4f}  (overfit {:.2f}x){}".format(
            k, np.mean(r["prof_rms"]), r["overfit"], mark))

    # Detail best robust
    print("\n  Best robust model: {}".format(best_robust))
    br = results[best_robust]
    print("    k = {}, full RMS = {:.4f}, LOO RMS = {:.4f}".format(
        br["k"], br["rms_full"], br["loo_rms"]))
    print("    mean LOO profile RMS = {:.4f}".format(np.mean(br["prof_rms"])))

    # Coefficients best robust
    names_map = {
        "U0": ["a(d)", "as", "asp", "b", "c"],
        "U1": ["a(d)", "as", "asp", "b", "c", "d(v)"],
        "U2": ["a(d)", "as", "asp", "b", "c", "d(v)", "ds", "dsp"],
        "U3": ["a(d)", "as", "asp", "b", "c", "d(v)", "bs", "bsp"],
        "U4": ["a(d)", "as", "asp", "b", "c", "d(v)", "cs", "csp"],
        "U24": ["a(d)", "as", "asp", "b", "c", "d(v)", "ds", "dsp", "cs", "csp"],
        "U_full": ["a(d)", "as", "asp", "b", "c", "d(v)", "bs", "bsp", "cs", "csp", "ds", "dsp"],
    }
    print("\n    Coefficients:")
    for nm, v in zip(names_map[best_robust], br["coefs"]):
        print("      {:<8s} = {:+.4f}".format(nm, v))

    # -------- (C) Per-klasa breakdown dla best robust --------
    print("\n" + "-" * 92)
    print("  (C) Per-klasa LOO profile RMS z best robust")
    print("-" * 92)
    br_prof = br["prof_rms"]
    print("  Per-klasa:")
    for cls in ["s", "sp", "d"]:
        mask_cls = cls_arr == cls
        vs = br_prof[mask_cls]
        if len(vs) > 0:
            print("    {:<3s} (N={:2d}): mean = {:.4f}, median = {:.4f}, max = {:.4f}".format(
                cls, len(vs), np.mean(vs), np.median(vs), np.max(vs)))

    # Top outliers
    print("\n  Top 7 LOO profile RMS outliers (best robust):")
    idx = np.argsort(-br_prof)
    for i in idx[:7]:
        print("    {:<5s} ({:<3s}): prof RMS = {:.4f}, dlog_R = {:+.3f}".format(
            rows[i]["name"], rows[i]["cls"], br_prof[i],
            br["loo_logR"][i] - logR[i]))

    # -------- (D) Compare robust vs U_full predictions per material --------
    print("\n" + "-" * 92)
    print("  (D) Porownanie {}-robust vs U_full (LOO profile RMS per mat.)".format(best_robust))
    print("-" * 92)
    full_prof = results["U_full"]["prof_rms"]
    wins_r = 0; wins_f = 0
    print("  {:<5s}  {:<3s}  {:>10s}  {:>10s}  {:>8s}".format(
        "name", "cls", "robust", "U_full", "diff"))
    for i, r in enumerate(rows):
        diff = full_prof[i] - br_prof[i]
        if br_prof[i] < full_prof[i] - 0.01: wins_r += 1
        elif full_prof[i] < br_prof[i] - 0.01: wins_f += 1
        # Only print noteworthy
        if abs(diff) > 0.05:
            tag = "R>F" if diff > 0 else "F>R"
            print("  {:<5s}  {:<3s}  {:>10.4f}  {:>10.4f}  {:>+8.4f}  {}".format(
                r["name"], r["cls"], br_prof[i], full_prof[i], diff, tag))

    print("\n  Wins: robust = {}, U_full = {}, tie = {}".format(
        wins_r, wins_f, N - wins_r - wins_f))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r13b")
    print("=" * 92)

    print("\n  Robust model: {}  (k={})".format(best_robust, br["k"]))
    print("  Mean LOO profile RMS = {:.4f}  (full fit: {:.4f})".format(
        np.mean(br_prof), np.mean(results["U_full"]["prof_rms"])))

    if best_robust == "U_full":
        print("\n  U_full jest najlepsza rowniez w LOO -> class-specific b,c są REAL.")
    elif best_robust in ["U2", "U4", "U24"]:
        print("\n  Class-specific intercept + ONE slope-modifier jest optimum:")
        if best_robust == "U2":
            print("    Formula: log R = a_cls + b*logTh + c*logNF + d_cls*v")
            print("    Slopes b, c uniwersalne; v-slope class-specific.")
        elif best_robust == "U4":
            print("    Formula: log R = a_cls + b*logTh + c_cls*logNF + d*v")
            print("    c (DOS-slope) class-specific; b (phonon), d uniwersalne.")
        elif best_robust == "U24":
            print("    Formula: log R = a_cls + b*logTh + c_cls*logNF + d_cls*v")
            print("    b uniwersalne; c, d class-specific.")
    else:
        print("\n  Minimum LOO model: {} (prosciejszy niz U_full) -> wybrac dla paperu.".format(
            best_robust))

    overfit_r = br["overfit"]
    if overfit_r < 1.5:
        print("\n  STABILNY: overfit ratio {:.2f}x < 1.5x -> model genaralizuje.".format(overfit_r))
    else:
        print("\n  Overfit ratio {:.2f}x -> ostrozna interpretacja.".format(overfit_r))


if __name__ == "__main__":
    main()
