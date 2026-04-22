#!/usr/bin/env python3
# =============================================================================
#  r15_minimal_model.py
# -----------------------------------------------------------------------------
#  Etap 10: Minimalny model oparty TYLKO na istotnych parametrach z bootstrap
#  r14 (3/10 parametrow U24 ma 95% CI poza 0: b(logTh), as(s-d), d(v)).
#
#  Testowane modele:
#    M_min3:  log R = a + b * log Theta + as * delta_s + d * v
#             (4 parametry: uniwersalna temp, klasa s vs d/sp, valence)
#
#    M_min4:  log R = a + b * log Theta + c * log NF + as * delta_s + d * v
#             (5 parametrow: dodaje DOS)
#
#    M_min5:  log R = a + as*delta_s + asp*delta_sp + b*logTh + c*logNF + d*v
#             (6 parametrow = U1 z r13b - baseline klasy+v bez interakcji)
#
#  Hipoteza: skoro 7/10 U24 to szum, minimalna forma powinna dac LOO RMS
#  wystarczajaco blisko U24 (rozn < 25%) przy duzo mniejszym k.
#
#  Po fitach: porownanie AIC/BIC/LOO z U24 + per-klasa profile RMS.
# =============================================================================
import numpy as np
import importlib.util, os, sys

for mod_name in ["r00", "r01", "r02", "r06", "r10", "r12",
                 "r13_unified_all_classes", "r13b_loo_sweep"]:
    fn = os.path.join(os.path.dirname(__file__), mod_name + ".py")
    if not os.path.exists(fn): continue
    spec = importlib.util.spec_from_file_location(mod_name, fn)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = m
    spec.loader.exec_module(m)

r00 = sys.modules["r00"]; r01 = sys.modules["r01"]
r02 = sys.modules["r02"]; r06 = sys.modules["r06"]
r10 = sys.modules["r10"]; r12 = sys.modules["r12"]
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


def loo_rms_and_profile(X, logR, rows):
    """Return LOO-predicted logR array + LOO profile RMS per material."""
    N = len(rows)
    loo_logR = np.zeros(N)
    for i in range(N):
        mask = np.ones(N, bool); mask[i] = False
        try:
            coefs_i, _, _, _ = np.linalg.lstsq(X[mask], logR[mask], rcond=None)
            loo_logR[i] = X[i] @ coefs_i
        except Exception:
            loo_logR[i] = np.nan
    loo_rms = np.sqrt(np.nanmean((loo_logR - logR)**2))

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

    return loo_logR, loo_rms, prof_rms


def main():
    print("=" * 92)
    print("  r15_minimal_model.py - minimalny model na istotnych parametrach")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    data = r12.build_dataset_N37()
    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = r02.ORB_CLASS.get(d["name"], "?")
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

    # ------------------------------------------------------------
    # Build models
    # ------------------------------------------------------------
    # M_min3: only b(logTh), as, d(v) - najmniejszy mozliwy
    X_min3 = np.column_stack([np.ones_like(logTh), delta_s, logTh, v_arr])
    names_min3 = ["a", "as(s-d)", "b(logTh)", "d(v)"]

    # M_min4: + NF (c uniwersalne)
    X_min4 = np.column_stack([np.ones_like(logTh), delta_s, logTh, logNF, v_arr])
    names_min4 = ["a", "as(s-d)", "b(logTh)", "c(logNF)", "d(v)"]

    # M_min4sp: s i sp rozdzielone
    X_min4sp = np.column_stack([np.ones_like(logTh), delta_s, delta_sp, logTh, logNF, v_arr])
    names_min4sp = ["a", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)", "d(v)"]

    # U24 reference - 10 params
    X_U24 = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                              logTh, logNF, v_arr,
                              v_arr*delta_s, v_arr*delta_sp,
                              delta_s*logNF, delta_sp*logNF])
    names_U24 = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)",
                  "d(v)", "ds(v*s)", "dsp(v*sp)", "cs(s*logNF)", "csp(sp*logNF)"]

    # U1 (z r13b) - klasy + v, brak interakcji - 6 params
    X_U1 = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                             logTh, logNF, v_arr])
    names_U1 = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)", "d(v)"]

    models = {
        "M_min3":  (X_min3, names_min3),
        "M_min4":  (X_min4, names_min4),
        "M_min4sp": (X_min4sp, names_min4sp),
        "U1":      (X_U1, names_U1),
        "U24":     (X_U24, names_U24),
    }

    # ------------------------------------------------------------
    # Fit & LOO
    # ------------------------------------------------------------
    print("\n" + "-" * 92)
    print("  (A) AIC/BIC/LOO comparison")
    print("-" * 92)
    print("  {:<10s}  {:>3s}  {:>8s}  {:>9s}  {:>9s}  {:>8s}  {:>10s}  {:>8s}".format(
        "model", "k", "RMS", "AIC", "BIC", "LOO_RMS", "LOO_profRMS", "overfit"))

    results = {}
    for key, (X, names) in models.items():
        coefs, pred, rms, aic, bic = fit_ols(X, logR)
        loo_logR, loo_rms, prof_rms = loo_rms_and_profile(X, logR, rows)
        overfit = loo_rms / max(rms, 1e-3)
        mean_prof = float(np.mean(prof_rms))
        results[key] = {
            "X": X, "names": names, "coefs": coefs, "rms": rms,
            "aic": aic, "bic": bic, "loo_rms": loo_rms,
            "prof_rms": prof_rms, "mean_prof": mean_prof, "overfit": overfit,
            "k": X.shape[1],
        }
        print("  {:<10s}  {:>3d}  {:>8.4f}  {:>+9.3f}  {:>+9.3f}  {:>8.4f}  {:>10.4f}  {:>7.2f}x".format(
            key, X.shape[1], rms, aic, bic, loo_rms, mean_prof, overfit))

    # ------------------------------------------------------------
    # Coefficients detail
    # ------------------------------------------------------------
    print("\n" + "-" * 92)
    print("  (B) Coefficients for each model")
    print("-" * 92)
    for key in ["M_min3", "M_min4", "M_min4sp", "U1", "U24"]:
        r = results[key]
        print("\n  {} (k={}):".format(key, r["k"]))
        for nm, v in zip(r["names"], r["coefs"]):
            print("    {:<18s} = {:+.4f}".format(nm, v))

    # ------------------------------------------------------------
    # Per-class breakdown for M_min4sp (najciekawsze)
    # ------------------------------------------------------------
    print("\n" + "-" * 92)
    print("  (C) Per-class LOO profile RMS for all models")
    print("-" * 92)
    print("  {:<10s}  {:<4s}  {:>4s}  {:>10s}  {:>10s}  {:>10s}".format(
        "model", "cls", "N", "mean", "median", "max"))
    for key in ["M_min3", "M_min4", "M_min4sp", "U1", "U24"]:
        r = results[key]
        for cls in ["s", "sp", "d"]:
            mask_cls = cls_arr == cls
            vs = r["prof_rms"][mask_cls]
            if len(vs) > 0:
                print("  {:<10s}  {:<4s}  {:>4d}  {:>10.4f}  {:>10.4f}  {:>10.4f}".format(
                    key, cls, int(len(vs)),
                    float(np.mean(vs)), float(np.median(vs)), float(np.max(vs))))

    # ------------------------------------------------------------
    # Decision
    # ------------------------------------------------------------
    print("\n" + "=" * 92)
    print("  VERDICT r15")
    print("=" * 92)

    u24_prof = results["U24"]["mean_prof"]
    print("\n  U24 baseline (10 param): LOO profile RMS = {:.4f}".format(u24_prof))
    print("\n  Minimal models vs U24:")
    for key in ["M_min3", "M_min4", "M_min4sp", "U1"]:
        r = results[key]
        ratio = r["mean_prof"] / u24_prof
        verdict = "OK" if ratio < 1.25 else ("borderline" if ratio < 1.50 else "FAILS")
        print("    {:<10s} (k={}): prof RMS = {:.4f} ({:.2f}x U24) -- {}".format(
            key, r["k"], r["mean_prof"], ratio, verdict))

    # Rekomendacja
    candidates = sorted(
        [(k, results[k]) for k in ["M_min3", "M_min4", "M_min4sp", "U1"]],
        key=lambda kv: kv[1]["bic"]
    )
    best_min = candidates[0]
    print("\n  Best BIC-minimal model: {}".format(best_min[0]))
    print("    Parameters: k = {}".format(best_min[1]["k"]))
    print("    BIC = {:+.3f}  (vs U24 BIC = {:+.3f})".format(
        best_min[1]["bic"], results["U24"]["bic"]))
    print("    LOO profile RMS = {:.4f}  (vs U24 = {:.4f})".format(
        best_min[1]["mean_prof"], u24_prof))

    print("\n  Nastepny krok:")
    print("    - jesli best-min ma LOO < 1.25x U24: to JEST realna forma rownania")
    print("    - r16: ostrzal stopy (Cu-Ni, Nichrome, Brass, Constantan)")
    print("    - r17: Matthiessen rule validation (disorder component)")


if __name__ == "__main__":
    main()
