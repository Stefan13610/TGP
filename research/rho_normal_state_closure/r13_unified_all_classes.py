#!/usr/bin/env python3
# =============================================================================
#  r13_unified_all_classes.py
# -----------------------------------------------------------------------------
#  Etap 8 projektu rho(T): unifikacja M1 (z r12 dla d-class) na wszystkie
#  3 klasy TGP (s, sp, d).
#
#  r12 dla d-class (N=24):
#     log R = -0.293 + 1.17 * log(Theta) + 0.53 * log(N_F) - 0.114 * g
#
#  Pytanie: czy b = +1.17 i c = +0.53 SA UNIWERSALNE dla wszystkich klas,
#  czy s-class i sp-class maja wlasne wykladniki?
#
#  Testowane modele na N=37 (5 s + 8 sp + 24 d):
#
#    U0 -- log R = a_cls + b*log Theta + c*log NF
#          (kazda klasa ma wlasny intercept, ale jedno b, jedno c)
#
#    U1 -- + d*v (universal valence-electron intercept shift)
#
#    U2 -- + d_cls*v (class-specific v-slopes)
#
#    U3 -- + b_cls*log Theta (class-specific Theta-slopes, ADDITIVE)
#
#    U4 -- + c_cls*log NF (class-specific NF-slopes)
#
#    U_full -- wszystkie oddzielnie (9 parametrow class x 3 slopes)
#
#  Porownanie AIC/BIC + LOO CV na best model.
#
#  g (valence electrons):
#    s-class (noble Cu/Ag/Au): g=1
#    s-class (alkaline-earth Ca/Sr): g=2
#    sp-class: varies (Be/Mg/Zn/Cd=2, Al=3, Sn/Pb=4, Bi=5)
#    d-class: grupa PT = 3..10 (Sc=3, Ti=4, ..., Ni=10)
# =============================================================================
import numpy as np
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "r00", os.path.join(os.path.dirname(__file__), "r00_dataset.py"))
r00 = importlib.util.module_from_spec(spec); sys.modules["r00"] = r00
spec.loader.exec_module(r00)

spec1 = importlib.util.spec_from_file_location(
    "r01", os.path.join(os.path.dirname(__file__), "r01_bg_baseline.py"))
r01 = importlib.util.module_from_spec(spec1); sys.modules["r01"] = r01
spec1.loader.exec_module(r01)

spec2 = importlib.util.spec_from_file_location(
    "r02", os.path.join(os.path.dirname(__file__), "r02_tgp_formula.py"))
r02 = importlib.util.module_from_spec(spec2); sys.modules["r02"] = r02
spec2.loader.exec_module(r02)

spec6 = importlib.util.spec_from_file_location(
    "r06", os.path.join(os.path.dirname(__file__), "r06_extension.py"))
r06 = importlib.util.module_from_spec(spec6); sys.modules["r06"] = r06
spec6.loader.exec_module(r06)

spec10 = importlib.util.spec_from_file_location(
    "r10", os.path.join(os.path.dirname(__file__), "r10_extension_v2.py"))
r10 = importlib.util.module_from_spec(spec10); sys.modules["r10"] = r10
spec10.loader.exec_module(r10)

spec12 = importlib.util.spec_from_file_location(
    "r12", os.path.join(os.path.dirname(__file__), "r12_gparam_model.py"))
r12 = importlib.util.module_from_spec(spec12); sys.modules["r12"] = r12
spec12.loader.exec_module(r12)


# Valence electron count (main-shell conducting) per material:
V_COUNT = {
    # s-class (noble 1-valent)
    "Cu": 1, "Ag": 1, "Au": 1,
    # s-class (alkaline-earth 2-valent)
    "Ca": 2, "Sr": 2,
    # sp-class
    "Be": 2, "Mg": 2, "Zn": 2, "Cd": 2,  # s² sub-sp
    "Al": 3,
    "Sn": 4, "Pb": 4,
    "Bi": 5,
    # d-class (grupa w PT = valence)
    "Sc": 3, "Y": 3, "La": 3, "Lu": 3,
    "Ti": 4, "Zr": 4, "Hf": 4,
    "V": 5, "Nb": 5, "Ta": 5,
    "Cr": 6, "Mo": 6, "W": 6,
    "Mn": 7, "Tc": 7, "Re": 7,
    "Fe": 8, "Ru": 8, "Os": 8,
    "Co": 9, "Rh": 9, "Ir": 9,
    "Ni": 10, "Pd": 10, "Pt": 10,
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


def main():
    print("=" * 92)
    print("  r13_unified_all_classes.py - M1 formalism na s+sp+d (N=37)")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    data = r12.build_dataset_N37()
    print("\n  N_total = {}".format(len(data)))

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
        v = V_COUNT.get(d["name"])
        if v is None: continue
        rows.append({
            "name": d["name"], "cls": cls, "R": R, "v": v,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
        })

    # Breakdown per class
    class_counts = {}
    for r in rows: class_counts.setdefault(r["cls"], []).append(r)
    print("\n  Podzial klasowy:")
    for cls in ["s", "sp", "d"]:
        lst = class_counts.get(cls, [])
        print("    {:<3s} (N={:2d}): {}".format(
            cls, len(lst),
            ", ".join("{}(v={})".format(r["name"], r["v"]) for r in lst)))

    # Arrays
    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    cls_arr = np.array([r["cls"] for r in rows])

    # Dummies: d = reference class (implicit intercept a_d)
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)
    # delta_d implicit

    # -------- (A) Model comparison --------
    print("\n" + "-" * 92)
    print("  (A) Porownanie modeli unifikacji klas")
    print("-" * 92)
    print("  Nazwy param: a=intercept(d), as=da_s, asp=da_sp, b=logTh, c=logNF,")
    print("                d=v (univ), ds=dd_s, dsp=dd_sp, bs=db_s, bsp=db_sp,")
    print("                cs=dc_s, csp=dc_sp")
    print()

    models = {}

    # U0: log R = a + a_s*delta_s + a_sp*delta_sp + b*logTh + c*logNF
    X_U0 = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                            logTh, logNF])
    names_U0 = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)"]
    coefs, pred, rms, aic, bic = fit_ols(X_U0, logR)
    models["U0"] = {"X": X_U0, "coefs": coefs, "pred": pred, "rms": rms,
                     "aic": aic, "bic": bic, "names": names_U0, "k": X_U0.shape[1]}

    # U1: + d*v
    X_U1 = np.column_stack([X_U0, v_arr])
    names_U1 = names_U0 + ["d(v)"]
    coefs, pred, rms, aic, bic = fit_ols(X_U1, logR)
    models["U1"] = {"X": X_U1, "coefs": coefs, "pred": pred, "rms": rms,
                     "aic": aic, "bic": bic, "names": names_U1, "k": X_U1.shape[1]}

    # U2: + d_s*v*delta_s + d_sp*v*delta_sp (class-specific v-slopes)
    X_U2 = np.column_stack([X_U1, v_arr*delta_s, v_arr*delta_sp])
    names_U2 = names_U1 + ["ds(v*s)", "dsp(v*sp)"]
    coefs, pred, rms, aic, bic = fit_ols(X_U2, logR)
    models["U2"] = {"X": X_U2, "coefs": coefs, "pred": pred, "rms": rms,
                     "aic": aic, "bic": bic, "names": names_U2, "k": X_U2.shape[1]}

    # U3: U1 + b_s*delta_s*logTh + b_sp*delta_sp*logTh (class-specific Theta-slopes)
    X_U3 = np.column_stack([X_U1, delta_s*logTh, delta_sp*logTh])
    names_U3 = names_U1 + ["bs(s*logTh)", "bsp(sp*logTh)"]
    coefs, pred, rms, aic, bic = fit_ols(X_U3, logR)
    models["U3"] = {"X": X_U3, "coefs": coefs, "pred": pred, "rms": rms,
                     "aic": aic, "bic": bic, "names": names_U3, "k": X_U3.shape[1]}

    # U4: U1 + c_s*delta_s*logNF + c_sp*delta_sp*logNF (class-specific NF-slopes)
    X_U4 = np.column_stack([X_U1, delta_s*logNF, delta_sp*logNF])
    names_U4 = names_U1 + ["cs(s*logNF)", "csp(sp*logNF)"]
    coefs, pred, rms, aic, bic = fit_ols(X_U4, logR)
    models["U4"] = {"X": X_U4, "coefs": coefs, "pred": pred, "rms": rms,
                     "aic": aic, "bic": bic, "names": names_U4, "k": X_U4.shape[1]}

    # U_full: wszystko interactive (U1 + class*log Theta + class*log NF + class*v)
    X_F = np.column_stack([X_U1,
                            delta_s*logTh, delta_sp*logTh,
                            delta_s*logNF, delta_sp*logNF,
                            v_arr*delta_s, v_arr*delta_sp])
    names_F = names_U1 + ["bs(s*logTh)", "bsp(sp*logTh)",
                          "cs(s*logNF)", "csp(sp*logNF)",
                          "ds(v*s)", "dsp(v*sp)"]
    coefs, pred, rms, aic, bic = fit_ols(X_F, logR)
    models["U_full"] = {"X": X_F, "coefs": coefs, "pred": pred, "rms": rms,
                         "aic": aic, "bic": bic, "names": names_F, "k": X_F.shape[1]}

    # Print table
    aic_min = min(m["aic"] for m in models.values())
    bic_min = min(m["bic"] for m in models.values())
    print("  {:<8s}  {:>3s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}".format(
        "model", "k", "RMS", "AIC", "BIC", "dAIC"))
    for key in ["U0", "U1", "U2", "U3", "U4", "U_full"]:
        m = models[key]
        tags = []
        if m["aic"] == aic_min: tags.append("best AIC")
        if m["bic"] == bic_min: tags.append("best BIC")
        tag = " (" + ", ".join(tags) + ")" if tags else ""
        print("  {:<8s}  {:>3d}  {:>9.4f}  {:>9.3f}  {:>9.3f}  {:>+9.3f}{}".format(
            key, m["k"], m["rms"], m["aic"], m["bic"],
            m["aic"] - aic_min, tag))

    # Best model coefficients
    best_key = min(models.keys(), key=lambda k: models[k]["aic"])
    best = models[best_key]
    print("\n  Best model ({}) coefs:".format(best_key))
    for nm, v in zip(best["names"], best["coefs"]):
        print("    {:<16s} = {:+.4f}".format(nm, v))

    # -------- (B) Efektywne b, c, d per klasa --------
    print("\n" + "-" * 92)
    print("  (B) Efektywne wykladniki per klasa (z best modelu)")
    print("-" * 92)

    def get_effective(model_key, cls_target):
        """Zwroc (a_eff, b_eff, c_eff, d_eff) dla klasy cls_target."""
        m = models[model_key]
        cs = m["coefs"]; nm = m["names"]
        idx = {n: i for i, n in enumerate(nm)}

        a_eff = cs[idx["a(d)"]]
        if cls_target == "s" and "as(s-d)" in idx:
            a_eff += cs[idx["as(s-d)"]]
        if cls_target == "sp" and "asp(sp-d)" in idx:
            a_eff += cs[idx["asp(sp-d)"]]

        b_eff = cs[idx["b(logTh)"]]
        if cls_target == "s" and "bs(s*logTh)" in idx:
            b_eff += cs[idx["bs(s*logTh)"]]
        if cls_target == "sp" and "bsp(sp*logTh)" in idx:
            b_eff += cs[idx["bsp(sp*logTh)"]]

        c_eff = cs[idx["c(logNF)"]]
        if cls_target == "s" and "cs(s*logNF)" in idx:
            c_eff += cs[idx["cs(s*logNF)"]]
        if cls_target == "sp" and "csp(sp*logNF)" in idx:
            c_eff += cs[idx["csp(sp*logNF)"]]

        d_eff = cs[idx["d(v)"]] if "d(v)" in idx else 0.0
        if cls_target == "s" and "ds(v*s)" in idx:
            d_eff += cs[idx["ds(v*s)"]]
        if cls_target == "sp" and "dsp(v*sp)" in idx:
            d_eff += cs[idx["dsp(v*sp)"]]

        return a_eff, b_eff, c_eff, d_eff

    print("  Z best modelu ({}):".format(best_key))
    print("  {:<4s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}".format(
        "cls", "a_eff", "b_eff", "c_eff", "d_eff"))
    for cls in ["d", "s", "sp"]:
        a, b, c, d = get_effective(best_key, cls)
        print("  {:<4s}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}".format(cls, a, b, c, d))

    # For U_full (upper bound), show what each class would want
    print("\n  Z U_full (max flexibility):")
    print("  {:<4s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}".format(
        "cls", "a_eff", "b_eff", "c_eff", "d_eff"))
    for cls in ["d", "s", "sp"]:
        a, b, c, d = get_effective("U_full", cls)
        print("  {:<4s}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}  {:>+9.3f}".format(cls, a, b, c, d))

    # -------- (C) LOO CV on best model --------
    print("\n" + "-" * 92)
    print("  (C) LOO CV na best modelu ({})".format(best_key))
    print("-" * 92)
    X_b = best["X"]
    loo_dlog = []
    for i in range(len(rows)):
        mask = np.ones(len(rows), bool); mask[i] = False
        coefs_i, _, _, _ = np.linalg.lstsq(X_b[mask], logR[mask], rcond=None)
        pred_i = X_b[i] @ coefs_i
        dlog = pred_i - logR[i]
        loo_dlog.append(dlog)

    loo_rms = np.sqrt(np.mean(np.array(loo_dlog)**2))
    print("  LOO RMS = {:.4f}".format(loo_rms))
    print("  Full RMS = {:.4f}".format(best["rms"]))
    print("  Overfit ratio = {:.2f}x".format(loo_rms / max(best["rms"], 1e-3)))

    # Top 5 outliers
    loo_abs = [(rows[i]["name"], rows[i]["cls"], loo_dlog[i])
               for i in range(len(rows))]
    loo_abs.sort(key=lambda x: -abs(x[2]))
    print("\n  Top 5 LOO outliers:")
    for nm, cls, dl in loo_abs[:5]:
        print("    {:<5s} ({}): dlog = {:+.3f}".format(nm, cls, dl))

    # -------- (D) Full profile rho(T_i) na best modelu --------
    print("\n" + "-" * 92)
    print("  (D) Pelny profil rho(T_i) z best modelu")
    print("-" * 92)
    print("  {:<5s}  {:<3s}  {:>3s}  {:>9s}  {:>9s}  {:>9s}  {:>8s}".format(
        "name", "cls", "v", "R_obs", "R_pred", "dlogR", "profRMS"))
    profs = []
    for i, r in enumerate(rows):
        R_pred = 10**best["pred"][i]
        Th = r["Theta_D"]; rho0 = r["rho_0"]
        rms_pts = []
        for T in r00.T_POINTS:
            rho_p = r01.rho_bg(T, rho0, R_pred, Th)
            rho_o = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_p) - np.log10(rho_o))
        prms = np.sqrt(np.mean(np.array(rms_pts)**2))
        profs.append((r["name"], r["cls"], prms))
        dlog = best["pred"][i] - logR[i]
        print("  {:<5s}  {:<3s}  {:>3d}  {:>9.2f}  {:>9.2f}  {:>+9.3f}  {:>8.4f}".format(
            r["name"], r["cls"], r["v"], r["R"], R_pred, dlog, prms))

    vals = [x[2] for x in profs]
    print("\n  Profile RMS ({}):".format(best_key))
    print("    mean   = {:.4f}".format(np.mean(vals)))
    print("    median = {:.4f}".format(np.median(vals)))
    print("    max    = {:.4f}  ({})".format(max(vals),
        profs[np.argmax(vals)][0]))

    # Per-class breakdown
    print("\n  Per-class profile RMS:")
    for cls in ["s", "sp", "d"]:
        vs = [x[2] for x in profs if x[1] == cls]
        if vs:
            print("    {:<3s} (N={:2d}): mean = {:.4f}, median = {:.4f}, max = {:.4f}".format(
                cls, len(vs), np.mean(vs), np.median(vs), max(vs)))

    # Histogram
    print("\n  Histogram profile RMS:")
    bins = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 5.0)]
    for lo, hi in bins:
        members = [(nm, cls, v) for nm, cls, v in profs if lo <= v < hi]
        cnt = len(members)
        mem_str = ", ".join("{}({})".format(nm, cls) for nm, cls, _ in members)
        if len(mem_str) > 80: mem_str = mem_str[:77] + "..."
        print("    [{:.1f}, {:.1f}): {:2d}  {}".format(lo, hi, cnt, mem_str))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r13")
    print("=" * 92)

    dAIC_U0 = models["U0"]["aic"] - aic_min
    dAIC_U1 = models["U1"]["aic"] - aic_min
    dAIC_U_full = models["U_full"]["aic"] - aic_min

    print("\n  Model selection (AIC):")
    print("    U0 (klasy-intercept, univ b,c):           dAIC = {:+.2f}".format(dAIC_U0))
    print("    U1 (+ v intercept, univ b,c,d):           dAIC = {:+.2f}".format(dAIC_U1))
    print("    U_full (class-specific b,c,d):            dAIC = {:+.2f}".format(dAIC_U_full))

    if best_key == "U1" and dAIC_U_full > 5:
        print("\n  -> WYKŁADNIKI b, c SĄ UNIWERSALNE (inter-class)")
        print("     Klasa wplywa jedynie przez intercept (a_cls) i d_cls * v.")
        print("     Jedna formula TGP dla calego normal-state transport!")
        cs = best["coefs"]; nm = best["names"]
        idx = {n: i for i, n in enumerate(nm)}
        print("\n     FORMULA (baseline d-class):")
        print("       log R = {:+.3f} + {:+.3f}*log(Theta) + {:+.3f}*log(N_F) + {:+.3f}*v".format(
            cs[idx["a(d)"]], cs[idx["b(logTh)"]], cs[idx["c(logNF)"]], cs[idx["d(v)"]]))
        print("     z shift'ami:")
        print("       s-class:  a += {:+.3f}".format(cs[idx["as(s-d)"]]))
        print("       sp-class: a += {:+.3f}".format(cs[idx["asp(sp-d)"]]))

    elif best_key == "U_full":
        print("\n  -> UNIWERSALNE b, c ODRZUCONE: klasy maja rozne slopes.")
        print("     Potrzeba osobnych 3 fitow per klasa.")
    else:
        print("\n  -> Mieszany wynik: {} wygralo.".format(best_key))

    overfit = loo_rms / max(best["rms"], 1e-3)
    if overfit < 2.0:
        print("\n  LOO OK (overfit ratio {:.2f}x).".format(overfit))
    else:
        print("\n  WARNING: LOO overfit ratio {:.2f}x >= 2x.".format(overfit))

    mean_prof = np.mean(vals)
    if mean_prof < 0.25:
        status = "READY do paperu v3 (unified TGP formula)"
    elif mean_prof < 0.35:
        status = "BLISKO - wymaga drobnej weryfikacji"
    else:
        status = "WYMAGANA dalsza praca - profile RMS zbyt wysokie"
    print("\n  Mean profile RMS = {:.4f}  ->  {}".format(mean_prof, status))


if __name__ == "__main__":
    main()
