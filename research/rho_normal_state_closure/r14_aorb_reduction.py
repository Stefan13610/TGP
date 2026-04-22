#!/usr/bin/env python3
# =============================================================================
#  r14_aorb_reduction.py
# -----------------------------------------------------------------------------
#  Etap 9 projektu rho(T): ostrzal rownania U24 - test czy klasowe parametry
#  (a_cls, c_cls, d_cls) MOZNA WYPROWADZIC Z TGP A_orb (rdzenia teorii).
#
#  U24 z r13b (N=37):
#     log R = a_cls + 0.79 * log Theta + c_cls * log N_F + d_cls * v
#
#     cls   A_orb    a_cls    c_cls   d_cls
#     s    -0.111   -1.51    +0.03   +0.50
#     sp   +0.207   -1.21    -0.39   +0.19
#     d    +0.310   +0.58    +0.41   -0.10
#
#  Hipoteza: (a, c, d)_cls sa liniowymi funkcjami A_orb_cls.
#
#  Testy:
#    (A) Bootstrap (N=2000) CI dla U24 parametrow
#    (B) Linear regression a_cls vs A_orb, c_cls vs A_orb, d_cls vs A_orb
#    (C) Fit SIMPLIFIED formula:
#          log R = alpha_0 + alpha_1*A_orb + b*log Theta
#                  + (gamma_0 + gamma_1*A_orb)*log N_F
#                  + (delta_0 + delta_1*A_orb)*v
#        -- 7 parametrow zamiast 10
#    (D) Fit HYPER-SIMPLIFIED:
#          log R = beta_0 + beta_1*A_orb + b*log Theta
#                  + gamma*log N_F + delta*v
#        -- 5 parametrow (intercept liniowy w A_orb, reszta uniwersalna)
#    (E) AIC + LOO porownanie U0, U1, ..., U24, U_red7, U_red5
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
r13b = sys.modules["r13b_loo_sweep"]


# TGP A_orb values (z paper v1)
A_ORB = {
    "s":  -0.111,
    "sp": +0.207,
    "d":  +0.310,
    "f":  +2.034,  # nie uzywane tu
}


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


def build_U24(logTh, logNF, v_arr, delta_s, delta_sp):
    X = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                          logTh, logNF, v_arr,
                          v_arr*delta_s, v_arr*delta_sp,
                          delta_s*logNF, delta_sp*logNF])
    names = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)",
              "d(v)", "ds(v*s)", "dsp(v*sp)", "cs(s*logNF)", "csp(sp*logNF)"]
    return X, names


def build_Ured7(logTh, logNF, v_arr, A_arr):
    """U_red7: 7 params - a(linear in A_orb), b univ, c(linear in A_orb),
       d(linear in A_orb).
       log R = alpha0 + alpha1*A + b*logTh
               + (gamma0 + gamma1*A)*logNF
               + (delta0 + delta1*A)*v
    """
    X = np.column_stack([np.ones_like(logTh), A_arr,
                          logTh,
                          logNF, A_arr*logNF,
                          v_arr, A_arr*v_arr])
    names = ["alpha0", "alpha1(A)", "b(logTh)",
              "gamma0(logNF)", "gamma1(A*logNF)",
              "delta0(v)", "delta1(A*v)"]
    return X, names


def build_Ured5(logTh, logNF, v_arr, A_arr):
    """U_red5: 5 params - only intercept depends on A_orb.
       log R = beta0 + beta1*A + b*logTh + gamma*logNF + delta*v
    """
    X = np.column_stack([np.ones_like(logTh), A_arr,
                          logTh, logNF, v_arr])
    names = ["beta0", "beta1(A)", "b(logTh)", "gamma(logNF)", "delta(v)"]
    return X, names


def build_Ured6_c(logTh, logNF, v_arr, A_arr):
    """U_red6: like U_red5 but c also linear in A.
       log R = beta0 + beta1*A + b*logTh
               + (gamma0 + gamma1*A)*logNF + delta*v
    """
    X = np.column_stack([np.ones_like(logTh), A_arr,
                          logTh,
                          logNF, A_arr*logNF,
                          v_arr])
    names = ["beta0", "beta1(A)", "b(logTh)",
              "gamma0(logNF)", "gamma1(A*logNF)", "delta(v)"]
    return X, names


def build_Ured6_d(logTh, logNF, v_arr, A_arr):
    """U_red6: like U_red5 but d also linear in A.
       log R = beta0 + beta1*A + b*logTh
               + gamma*logNF + (delta0 + delta1*A)*v
    """
    X = np.column_stack([np.ones_like(logTh), A_arr,
                          logTh, logNF,
                          v_arr, A_arr*v_arr])
    names = ["beta0", "beta1(A)", "b(logTh)",
              "gamma(logNF)", "delta0(v)", "delta1(A*v)"]
    return X, names


def loo_profile_rms(X, logR, rows, r00_mod, r01_mod):
    N = len(rows)
    loo_logR = np.zeros(N)
    for i in range(N):
        mask = np.ones(N, bool); mask[i] = False
        coefs_i, _, _, _ = np.linalg.lstsq(X[mask], logR[mask], rcond=None)
        loo_logR[i] = X[i] @ coefs_i
    loo_dlog = loo_logR - logR
    loo_rms = np.sqrt(np.mean(loo_dlog**2))
    prof_rms = np.zeros(N)
    for i, r in enumerate(rows):
        R_pred = 10**loo_logR[i]
        Th = r["Theta_D"]; rho0 = r["rho_0"]
        pts = []
        for T in r00_mod.T_POINTS:
            rho_p = r01_mod.rho_bg(T, rho0, R_pred, Th)
            rho_o = r["rho_obs"][T]
            pts.append(np.log10(rho_p) - np.log10(rho_o))
        prof_rms[i] = np.sqrt(np.mean(np.array(pts)**2))
    return loo_logR, loo_rms, prof_rms


def main():
    print("=" * 92)
    print("  r14_aorb_reduction.py - ostrzal U24 + prob uproszczenia przez A_orb")
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
            "A_orb": A_ORB[cls],
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
        })

    N = len(rows)
    print("  N = {}".format(N))

    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    A_arr = np.array([r["A_orb"] for r in rows])
    logR = np.array([np.log10(r["R"]) for r in rows])
    cls_arr = np.array([r["cls"] for r in rows])
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)

    # -------- (A) Bootstrap U24 --------
    print("\n" + "-" * 92)
    print("  (A) Bootstrap U24 (N_boot=2000) - confidence intervals na parametrach")
    print("-" * 92)

    X_U24, names_U24 = build_U24(logTh, logNF, v_arr, delta_s, delta_sp)
    coefs_U24, pred_U24, rms_U24, aic_U24, _ = fit_ols(X_U24, logR)

    np.random.seed(42)
    Nboot = 2000
    boots = np.zeros((Nboot, len(names_U24)))
    for b in range(Nboot):
        idx = np.random.randint(0, N, size=N)
        try:
            cb, _, _, _ = np.linalg.lstsq(X_U24[idx], logR[idx], rcond=None)
            boots[b] = cb
        except Exception:
            boots[b] = coefs_U24

    print("  {:<16s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}  {:>8s}".format(
        "param", "OLS", "boot_mean", "boot_std", "CI2.5", "CI97.5"))
    for i, nm in enumerate(names_U24):
        bm = np.mean(boots[:, i]); bs = np.std(boots[:, i])
        c25 = np.percentile(boots[:, i], 2.5)
        c975 = np.percentile(boots[:, i], 97.5)
        sig = "*" if (c25 * c975 > 0) else "  "  # CI nie przechodzi przez 0
        print("  {:<16s}  {:>+10.4f}  {:>+10.4f}  {:>10.4f}  {:>+10.4f}  {:>+8.4f} {}".format(
            nm, coefs_U24[i], bm, bs, c25, c975, sig))

    print("\n  (*) = 95% CI nie zawiera zera - parametr ISTOTNY")

    # Efektywne b, c, d per klasa z bootstrap
    print("\n  Efektywne exponenty per klasa (+ 95% CI):")
    def eff_from_boot(cls_target, boot_row):
        idx_map = {n: i for i, n in enumerate(names_U24)}
        a = boot_row[idx_map["a(d)"]]
        if cls_target == "s":  a += boot_row[idx_map["as(s-d)"]]
        if cls_target == "sp": a += boot_row[idx_map["asp(sp-d)"]]
        b = boot_row[idx_map["b(logTh)"]]
        c = boot_row[idx_map["c(logNF)"]]
        if cls_target == "s":  c += boot_row[idx_map["cs(s*logNF)"]]
        if cls_target == "sp": c += boot_row[idx_map["csp(sp*logNF)"]]
        d = boot_row[idx_map["d(v)"]]
        if cls_target == "s":  d += boot_row[idx_map["ds(v*s)"]]
        if cls_target == "sp": d += boot_row[idx_map["dsp(v*sp)"]]
        return a, b, c, d

    print("  {:<4s}  {:>14s}  {:>14s}  {:>14s}  {:>14s}".format(
        "cls", "a_eff", "b_eff", "c_eff", "d_eff"))
    cls_params = {}
    for cls in ["s", "sp", "d"]:
        vals = np.array([eff_from_boot(cls, boots[b]) for b in range(Nboot)])
        means = vals.mean(axis=0)
        los = np.percentile(vals, 2.5, axis=0)
        his = np.percentile(vals, 97.5, axis=0)
        print("  {:<4s}  [{:>+.3f}..{:>+.3f}]  [{:>+.3f}..{:>+.3f}]  [{:>+.3f}..{:>+.3f}]  [{:>+.3f}..{:>+.3f}]".format(
            cls, los[0], his[0], los[1], his[1], los[2], his[2], los[3], his[3]))
        cls_params[cls] = {
            "a": (means[0], los[0], his[0]),
            "b": (means[1], los[1], his[1]),
            "c": (means[2], los[2], his[2]),
            "d": (means[3], los[3], his[3]),
        }

    # -------- (B) a_cls, c_cls, d_cls vs A_orb --------
    print("\n" + "-" * 92)
    print("  (B) Korelacje (a_cls, c_cls, d_cls) vs A_orb")
    print("-" * 92)

    A_list = np.array([A_ORB["s"], A_ORB["sp"], A_ORB["d"]])
    a_list = np.array([cls_params["s"]["a"][0], cls_params["sp"]["a"][0], cls_params["d"]["a"][0]])
    c_list = np.array([cls_params["s"]["c"][0], cls_params["sp"]["c"][0], cls_params["d"]["c"][0]])
    d_list = np.array([cls_params["s"]["d"][0], cls_params["sp"]["d"][0], cls_params["d"]["d"][0]])

    print("  {:<4s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}".format(
        "cls", "A_orb", "a_cls", "c_cls", "d_cls"))
    for cls, A, a, c, d in zip(["s", "sp", "d"], A_list, a_list, c_list, d_list):
        print("  {:<4s}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}".format(
            cls, A, a, c, d))

    # Linear regression
    def lin_fit(x, y):
        X = np.column_stack([np.ones_like(x), x])
        coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coefs
        res = y - pred
        ss = np.sum(res**2); tss = np.sum((y - y.mean())**2)
        r2 = 1 - ss/max(tss, 1e-30)
        return coefs[0], coefs[1], r2

    a_int, a_slp, a_r2 = lin_fit(A_list, a_list)
    c_int, c_slp, c_r2 = lin_fit(A_list, c_list)
    d_int, d_slp, d_r2 = lin_fit(A_list, d_list)

    print("\n  Linear fit (y = intercept + slope * A_orb):")
    print("    a_cls = {:+.3f} + {:+.3f} * A_orb    R^2 = {:.4f}".format(a_int, a_slp, a_r2))
    print("    c_cls = {:+.3f} + {:+.3f} * A_orb    R^2 = {:.4f}".format(c_int, c_slp, c_r2))
    print("    d_cls = {:+.3f} + {:+.3f} * A_orb    R^2 = {:.4f}".format(d_int, d_slp, d_r2))

    # -------- (C) Fit U_red modeli --------
    print("\n" + "-" * 92)
    print("  (C) Redukowane modele: U_red5, U_red6c, U_red6d, U_red7 vs U24")
    print("-" * 92)

    models = {}

    X_U5, n_U5 = build_Ured5(logTh, logNF, v_arr, A_arr)
    coefs, pred, rms, aic, bic = fit_ols(X_U5, logR)
    models["U_red5"] = {"X": X_U5, "names": n_U5, "coefs": coefs,
                         "rms": rms, "aic": aic, "bic": bic, "k": X_U5.shape[1]}

    X_U6c, n_U6c = build_Ured6_c(logTh, logNF, v_arr, A_arr)
    coefs, pred, rms, aic, bic = fit_ols(X_U6c, logR)
    models["U_red6c"] = {"X": X_U6c, "names": n_U6c, "coefs": coefs,
                          "rms": rms, "aic": aic, "bic": bic, "k": X_U6c.shape[1]}

    X_U6d, n_U6d = build_Ured6_d(logTh, logNF, v_arr, A_arr)
    coefs, pred, rms, aic, bic = fit_ols(X_U6d, logR)
    models["U_red6d"] = {"X": X_U6d, "names": n_U6d, "coefs": coefs,
                          "rms": rms, "aic": aic, "bic": bic, "k": X_U6d.shape[1]}

    X_U7, n_U7 = build_Ured7(logTh, logNF, v_arr, A_arr)
    coefs, pred, rms, aic, bic = fit_ols(X_U7, logR)
    models["U_red7"] = {"X": X_U7, "names": n_U7, "coefs": coefs,
                         "rms": rms, "aic": aic, "bic": bic, "k": X_U7.shape[1]}

    models["U24"] = {"X": X_U24, "names": names_U24, "coefs": coefs_U24,
                      "rms": rms_U24, "aic": aic_U24, "bic": 0.0, "k": X_U24.shape[1]}

    # LOO profile RMS dla kazdego
    print("\n  LOO sweep (profile RMS):")
    print("  {:<10s}  {:>3s}  {:>8s}  {:>9s}  {:>10s}  {:>10s}  {:>7s}".format(
        "model", "k", "RMS_full", "AIC", "LOO_RMS", "LOO_profRMS", "overfit"))

    loo_results = {}
    for key in ["U_red5", "U_red6c", "U_red6d", "U_red7", "U24"]:
        m = models[key]
        loo_logR, loo_rms, prof_rms = loo_profile_rms(m["X"], logR, rows, r00, r01)
        loo_results[key] = {"loo_logR": loo_logR, "loo_rms": loo_rms,
                             "prof_rms": prof_rms}
        overfit = loo_rms / max(m["rms"], 1e-3)
        print("  {:<10s}  {:>3d}  {:>8.4f}  {:>9.3f}  {:>10.4f}  {:>10.4f}  {:>6.2f}x".format(
            key, m["k"], m["rms"], m["aic"], loo_rms,
            np.mean(prof_rms), overfit))

    # -------- (D) Best reduced model details --------
    print("\n" + "-" * 92)
    print("  (D) Best reduced model: coefficients + efektywne exponenty")
    print("-" * 92)

    # Pick best LOO prof RMS (among reduced)
    reduced_keys = ["U_red5", "U_red6c", "U_red6d", "U_red7"]
    best_red = min(reduced_keys, key=lambda k: np.mean(loo_results[k]["prof_rms"]))
    print("\n  Best reduced: {}".format(best_red))
    m = models[best_red]
    for nm, v in zip(m["names"], m["coefs"]):
        print("    {:<16s} = {:+.4f}".format(nm, v))

    # U24 benchmark
    print("\n  U24 benchmark (full param) coefs:")
    for nm, v in zip(names_U24, coefs_U24):
        print("    {:<16s} = {:+.4f}".format(nm, v))

    # -------- (E) Per-klasa efective exponents from best reduced --------
    print("\n" + "-" * 92)
    print("  (E) Efektywne (a, c, d)_cls z best reduced model ({})".format(best_red))
    print("-" * 92)

    def eff_from_reduced(key, cls_target):
        m = models[key]
        cs = m["coefs"]; nm = m["names"]
        A_cls = A_ORB[cls_target]
        a = cs[0]
        if "alpha1(A)" in nm: a += cs[nm.index("alpha1(A)")] * A_cls
        if "beta1(A)" in nm: a += cs[nm.index("beta1(A)")] * A_cls
        b = cs[nm.index("b(logTh)")]
        c = 0.0
        if "gamma0(logNF)" in nm: c = cs[nm.index("gamma0(logNF)")]
        elif "gamma(logNF)" in nm: c = cs[nm.index("gamma(logNF)")]
        if "gamma1(A*logNF)" in nm: c += cs[nm.index("gamma1(A*logNF)")] * A_cls
        d = 0.0
        if "delta0(v)" in nm: d = cs[nm.index("delta0(v)")]
        elif "delta(v)" in nm: d = cs[nm.index("delta(v)")]
        if "delta1(A*v)" in nm: d += cs[nm.index("delta1(A*v)")] * A_cls
        return a, b, c, d

    print("  {:<4s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}  (vs U24 boot mean)".format(
        "cls", "a_eff", "b_eff", "c_eff", "d_eff"))
    for cls in ["s", "sp", "d"]:
        a_r, b_r, c_r, d_r = eff_from_reduced(best_red, cls)
        a_u = cls_params[cls]["a"][0]
        c_u = cls_params[cls]["c"][0]
        d_u = cls_params[cls]["d"][0]
        print("  {:<4s}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}  ({:+.3f}, {:+.3f}, {:+.3f})".format(
            cls, a_r, b_r, c_r, d_r, a_u, c_u, d_u))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r14")
    print("=" * 92)

    U24_prof = np.mean(loo_results["U24"]["prof_rms"])
    best_prof = np.mean(loo_results[best_red]["prof_rms"])
    ratio = best_prof / U24_prof

    print("\n  LOO profile RMS:")
    print("    U24 (10 param):       {:.4f}".format(U24_prof))
    print("    {} ({} param):      {:.4f}".format(best_red, models[best_red]["k"], best_prof))
    print("    Ratio {}/U24:       {:.2f}x".format(best_red, ratio))

    if ratio < 1.10:
        print("\n  *** SUCCES: Reduced model rozni sie o <10% LOO RMS od U24 !!")
        print("      Klasowe parametry SA derivable z A_orb. TGP rdzen wystarczy.")
    elif ratio < 1.25:
        print("\n  Reduced model jest akceptowalny (do 25% pogorszenia).")
        print("  Warto go uwzglednic dla parsymonia - paper v3 scope zredukowany.")
    else:
        print("\n  Reduced model TRACI zbyt duzo - class-specific parametry")
        print("  zawieraja cos wiecej niz A_orb.")

    print("\n  R^2 klasowych parametrow vs A_orb:")
    print("    a_cls vs A_orb: R^2 = {:.4f}  {}".format(
        a_r2, "silna zaleznosc" if a_r2 > 0.8 else "slaba"))
    print("    c_cls vs A_orb: R^2 = {:.4f}  {}".format(
        c_r2, "silna zaleznosc" if c_r2 > 0.8 else "slaba"))
    print("    d_cls vs A_orb: R^2 = {:.4f}  {}".format(
        d_r2, "silna zaleznosc" if d_r2 > 0.8 else "slaba"))

    # Bootstrap significance
    print("\n  Bootstrap stability U24:")
    sig_count = 0
    for i, nm in enumerate(names_U24):
        c25 = np.percentile(boots[:, i], 2.5)
        c975 = np.percentile(boots[:, i], 97.5)
        if c25 * c975 > 0: sig_count += 1
    print("  {}/{} parametrow istotnych (95% CI poza 0)".format(sig_count, len(names_U24)))

    if sig_count < len(names_U24) - 2:
        print("\n  {} parametrow NIE sa istotne - probably nadmiarowe.".format(
            len(names_U24) - sig_count))
    else:
        print("\n  Wiekszosc parametrow istotnych - U24 nie przeparametryzowane.")

    print("""
  Nastepne kroki:
    - Jesli reduced model dzialal: paper v3 z 5-7 parametrami
    - Jesli nie: zatrzymac 10 param, ale z bootstrap CI
    - r15: ostrzal na stopy Cu-Ni, nichrome, brass (test uniwersalnosci)
""")


if __name__ == "__main__":
    main()
