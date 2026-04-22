#!/usr/bin/env python3
# =============================================================================
#  r18_nordheim.py
# -----------------------------------------------------------------------------
#  Etap 10: proba unifikacji pure+alloys przez rozszerzenie U24 o term
#  Nordheim'owski / disorder-dependent.
#
#  Z r16 wiemy:
#    - d-class alloys: R_obs << R_pred (U24 overshoots 8-30x)
#    - s-class alloys: R_obs > R_pred (U24 undershoots 2-3x)
#    - sp-class alloys: OK
#    - dlaczego? W stopach disorder saturuje phonon scattering (Ioffe-Regel):
#      kF * l_e ~ 1 -> ρ_i(T) nie rosnie liniowo z T, bo jest limit
#
#  Hipotezy do testu:
#
#  H1 Nordheim-basic:
#      R_alloy = R_U24_pure_weighted  +  alpha * x*(1-x)*(dZ)^2
#      gdzie x = fraction, dZ = Z_A - Z_B
#
#  H2 Saturation (Ioffe-Regel):
#      R_alloy = R_U24_pure / (1 + lambda * rho_0 / rho_sat)
#      rho_sat ~ h/(e^2 * k_F) ~ 150 uOhm*cm (rough universal limit)
#
#  H3 Disorder-logarithmic:
#      log R_alloy = log R_U24_pure - kappa * log(1 + rho_0/rho_0_ref)
#
#  H4 Class-specific saturation:
#      Osobno dla s/sp/d class.
#
#  Test: fit H1-H4 na UNIFIED dataset (37 pure + 13 alloys = 50) i sprawdzic:
#    - czy pure still RMS<0.2
#    - czy alloy logR RMS spada ponizej 0.3 (vs 0.70 dla czystego U24)
#    - AIC/BIC preferencje
# =============================================================================
import numpy as np
import importlib.util, os, sys

for mod_name in ["r00", "r01", "r02", "r06", "r10", "r12",
                 "r13_unified_all_classes", "r13b_loo_sweep", "r16_alloys_test"]:
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
r16 = sys.modules["r16_alloys_test"]


# Z (atomic number) dla stopowych constituentow
Z_TABLE = {
    "Cu": 29, "Ni": 28, "Zn": 30, "Sn": 50, "Pb": 82, "Ag": 47, "Au": 79,
    "Cr": 24, "Mn": 25, "Fe": 26,
}


def compute_nordheim_term(alloy):
    """x*(1-x)*(dZ)^2 averaged over all pair contributions."""
    comp = alloy.get("composition", None)
    if comp is None:
        # Try to parse from key
        return 0.0
    # For binary: straightforward
    items = list(comp.items())
    if len(items) < 2:
        return 0.0
    # Most common: dominant pair
    items.sort(key=lambda kv: -kv[1])
    A, xA = items[0]
    B, xB = items[1]
    ZA = Z_TABLE.get(A, 0)
    ZB = Z_TABLE.get(B, 0)
    if ZA == 0 or ZB == 0:
        return 0.0
    return xA * xB * (ZA - ZB)**2


# For alloys we stored composition implicitly in the name; let's reconstruct:
ALLOY_COMPOSITIONS = {
    "CuNi_70_30": {"Cu": 0.70, "Ni": 0.30},
    "CuNi_55_45": {"Cu": 0.55, "Ni": 0.45},
    "NiCr_80_20": {"Ni": 0.80, "Cr": 0.20},
    "Brass_Cu70Zn30": {"Cu": 0.70, "Zn": 0.30},
    "Bronze_Cu90Sn10": {"Cu": 0.90, "Sn": 0.10},
    "Manganin":       {"Cu": 0.86, "Mn": 0.12, "Ni": 0.02},
    "Monel_400":      {"Ni": 0.65, "Cu": 0.33, "Fe": 0.02},
    "AgAu_50_50":     {"Ag": 0.50, "Au": 0.50},
    "AgAu_90_10":     {"Ag": 0.90, "Au": 0.10},
    "Cu3Au_disord":   {"Cu": 0.75, "Au": 0.25},
    "Cu3Au_order":    {"Cu": 0.75, "Au": 0.25},
    "PbSn_63_37":     {"Pb": 0.37, "Sn": 0.63},
    "Brass_Cu63Zn37": {"Cu": 0.63, "Zn": 0.37},
}


def nordheim_from_comp(comp):
    """Full Nordheim: sum_{i<j} x_i x_j (Z_i - Z_j)^2"""
    items = list(comp.items())
    total = 0.0
    for i in range(len(items)):
        for j in range(i+1, len(items)):
            A, xA = items[i]
            B, xB = items[j]
            ZA = Z_TABLE.get(A, 0); ZB = Z_TABLE.get(B, 0)
            if ZA > 0 and ZB > 0:
                total += xA * xB * (ZA - ZB)**2
    return total


def main():
    print("=" * 100)
    print("  r18_nordheim.py - proba unifikacji pure+alloys przez Nordheim term")
    print("=" * 100)

    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    # ------- Build unified dataset (pure + alloys) -------
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
            "nordheim": 0.0, "is_alloy": False,
        })
    N_pure = len(rows)
    print("  N_pure  = {}".format(N_pure))

    # Add alloys
    for key, alloy in r16.ALLOYS.items():
        R_obs, _, Th = r16.fit_R_alloy(alloy)
        if R_obs is None or R_obs <= 0: continue
        comp = ALLOY_COMPOSITIONS.get(key, {})
        N_nord = nordheim_from_comp(comp)
        rows.append({
            "name": key, "cls": alloy["cls"], "R": R_obs,
            "v": alloy["v_eff"],
            "Theta_D": alloy["Theta_D"], "N_F": alloy["N_F"],
            "rho_obs": alloy["rho_T"], "rho_0": alloy["rho_0"],
            "nordheim": N_nord, "is_alloy": True,
        })
    N = len(rows)
    N_alloy = N - N_pure
    print("  N_alloy = {}".format(N_alloy))
    print("  N_total = {}".format(N))

    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    rho0_arr = np.array([r["rho_0"] for r in rows])
    nord_arr = np.array([r["nordheim"] for r in rows])
    is_alloy = np.array([r["is_alloy"] for r in rows], dtype=float)
    cls_arr = np.array([r["cls"] for r in rows])
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)

    print("\n  Nordheim values:")
    for r in rows:
        if r["is_alloy"]:
            print("    {:<20s}  Nord = {:>8.1f}  rho_0 = {:>6.1f}".format(
                r["name"], r["nordheim"], r["rho_0"]))

    # --------- Baseline: U24 on PURE only ----------
    print("\n" + "-" * 100)
    print("  (Baseline) U24 fit on PURE only (r13b), predict alloys")
    print("-" * 100)
    mask_p = ~is_alloy.astype(bool)
    X_U24_all = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                                   logTh, logNF, v_arr,
                                   v_arr*delta_s, v_arr*delta_sp,
                                   delta_s*logNF, delta_sp*logNF])
    c_pure, _, _, _ = np.linalg.lstsq(X_U24_all[mask_p], logR[mask_p], rcond=None)
    pred_pure_only = X_U24_all @ c_pure
    res_p = pred_pure_only[mask_p] - logR[mask_p]
    res_a = pred_pure_only[~mask_p] - logR[~mask_p]
    print("  Pure RMS  = {:.4f}".format(float(np.sqrt(np.mean(res_p**2)))))
    print("  Alloy RMS = {:.4f}  (catastrophic)".format(float(np.sqrt(np.mean(res_a**2)))))

    # ------- Helper --------
    def fit_and_eval(X, logR_t, mask_p, label, names):
        coefs, _, _, _ = np.linalg.lstsq(X, logR_t, rcond=None)
        pred = X @ coefs
        res = pred - logR_t
        n, k = X.shape
        ss = np.sum(res**2)
        rms = np.sqrt(ss / n)
        aic = n * np.log(ss/n + 1e-30) + 2 * k
        bic = n * np.log(ss/n + 1e-30) + k * np.log(n)
        rms_p = float(np.sqrt(np.mean(res[mask_p]**2)))
        rms_a = float(np.sqrt(np.mean(res[~mask_p]**2)))
        print("  {:<14s}  k={:>2d}  pure RMS={:.4f}  alloy RMS={:.4f}  total RMS={:.4f}  AIC={:+.2f}".format(
            label, k, rms_p, rms_a, rms, aic))
        return coefs, rms, rms_p, rms_a, aic, names

    # --------- H1 Nordheim-basic ---------
    print("\n" + "-" * 100)
    print("  (H1) Unified U24 + alpha * Nordheim (N_norm = x(1-x)(dZ)^2)")
    print("-" * 100)
    # Normalize Nordheim to avoid scaling issues: log(1 + N/100)
    logNord = np.log10(1 + nord_arr / 100.0)
    X_H1 = np.column_stack([X_U24_all, logNord])
    names_H1 = ["a", "as", "asp", "b", "c", "d", "ds", "dsp", "cs", "csp", "alpha_N"]
    c_H1, _, _, _, _, _ = fit_and_eval(X_H1, logR, mask_p, "H1", names_H1)
    for nm, v in zip(names_H1, c_H1):
        print("    {:<10s} = {:+.4f}".format(nm, v))

    # --------- H2 Saturation (rho_0-dependent) ---------
    print("\n" + "-" * 100)
    print("  (H2) Unified U24 + lambda * log10(rho_0)  (disorder saturation)")
    print("-" * 100)
    logrho0 = np.log10(rho0_arr + 1e-3)
    X_H2 = np.column_stack([X_U24_all, logrho0])
    names_H2 = ["a", "as", "asp", "b", "c", "d", "ds", "dsp", "cs", "csp", "lam_rho0"]
    c_H2, _, _, _, _, _ = fit_and_eval(X_H2, logR, mask_p, "H2", names_H2)
    for nm, v in zip(names_H2, c_H2):
        print("    {:<10s} = {:+.4f}".format(nm, v))

    # --------- H3 Disorder-log + class-specific ---------
    print("\n" + "-" * 100)
    print("  (H3) Unified + class-specific rho_0 response")
    print("-" * 100)
    X_H3 = np.column_stack([X_U24_all, logrho0, logrho0*delta_s, logrho0*delta_sp])
    names_H3 = ["a", "as", "asp", "b", "c", "d", "ds", "dsp", "cs", "csp",
                 "lam_rho0_d", "lam_rho0_s", "lam_rho0_sp"]
    c_H3, _, _, _, _, _ = fit_and_eval(X_H3, logR, mask_p, "H3", names_H3)
    for nm, v in zip(names_H3, c_H3):
        print("    {:<10s} = {:+.4f}".format(nm, v))

    # --------- H4 is_alloy indicator + class (simplest) ---------
    print("\n" + "-" * 100)
    print("  (H4) Unified + is_alloy binary (per-class offset)")
    print("-" * 100)
    X_H4 = np.column_stack([X_U24_all, is_alloy, is_alloy*delta_s, is_alloy*delta_sp])
    names_H4 = ["a", "as", "asp", "b", "c", "d", "ds", "dsp", "cs", "csp",
                 "alloy_d", "alloy_s", "alloy_sp"]
    c_H4, _, _, _, _, _ = fit_and_eval(X_H4, logR, mask_p, "H4", names_H4)
    for nm, v in zip(names_H4, c_H4):
        print("    {:<10s} = {:+.4f}".format(nm, v))

    # --------- Combined H2 + H1 ---------
    print("\n" + "-" * 100)
    print("  (H5) H2 + H1 combined (rho_0 + Nordheim)")
    print("-" * 100)
    X_H5 = np.column_stack([X_U24_all, logrho0, logNord])
    names_H5 = names_H2 + ["alpha_N"]
    c_H5, _, _, _, _, _ = fit_and_eval(X_H5, logR, mask_p, "H5", names_H5)
    for nm, v in zip(names_H5, c_H5):
        print("    {:<10s} = {:+.4f}".format(nm, v))

    # --------- Big picture summary ---------
    print("\n" + "-" * 100)
    print("  Summary (wszystkie hipotezy)")
    print("-" * 100)
    print("  {:<18s}  {:>3s}  {:>10s}  {:>10s}  {:>10s}  {:>8s}".format(
        "model", "k", "pure RMS", "alloy RMS", "total RMS", "AIC"))
    models = {
        "U24 (pure only)": (X_U24_all, c_pure),
        "H1 Nordheim":   (X_H1, c_H1),
        "H2 rho_0-log":  (X_H2, c_H2),
        "H3 rho_0 x cls": (X_H3, c_H3),
        "H4 is_alloy x cls": (X_H4, c_H4),
        "H5 H1+H2 comb": (X_H5, c_H5),
    }
    for label, (X, c) in models.items():
        if label == "U24 (pure only)":
            # pure-only calib
            pred = X @ c_pure
        else:
            pred = X @ c
        res = pred - logR
        n, k = X.shape
        rms_p = float(np.sqrt(np.mean(res[mask_p]**2)))
        rms_a = float(np.sqrt(np.mean(res[~mask_p]**2)))
        rms = float(np.sqrt(np.mean(res**2)))
        ss = np.sum(res**2)
        aic = n * np.log(ss/n + 1e-30) + 2 * k
        print("  {:<18s}  {:>3d}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>+8.2f}".format(
            label, k, rms_p, rms_a, rms, aic))

    # --------- Verdict ---------
    print("\n" + "=" * 100)
    print("  VERDICT r18")
    print("=" * 100)
    print()
    # Find best alloy RMS
    best_alloy = min(["H1", "H2", "H3", "H4", "H5"],
                     key=lambda lbl: {
                         "H1": np.sqrt(np.mean((X_H1 @ c_H1 - logR)[~mask_p]**2)),
                         "H2": np.sqrt(np.mean((X_H2 @ c_H2 - logR)[~mask_p]**2)),
                         "H3": np.sqrt(np.mean((X_H3 @ c_H3 - logR)[~mask_p]**2)),
                         "H4": np.sqrt(np.mean((X_H4 @ c_H4 - logR)[~mask_p]**2)),
                         "H5": np.sqrt(np.mean((X_H5 @ c_H5 - logR)[~mask_p]**2)),
                     }[lbl])
    print("  Najlepszy model alloy-fit: {}".format(best_alloy))
    print()
    print("  Kluczowa interpretacja:")
    print("    - Jesli lambda_rho0 < 0: zwiekszenie disorder (wieksze rho_0) zmniejsza R")
    print("      = phonon scattering saturuje w obecnosci disorder.")
    print("    - Jesli lambda_rho0 > 0: R rosnie z disorder (Matthiessen additive, nie niezalezne)")
    print("    - Nordheim term alpha_N > 0: concentration-fluctuation scattering dodaje do R")
    print()
    print("  Praktycznie:")
    print("    - jesli H2 daje alloy RMS < 0.25: SUKCES, mamy unified pure+alloy formula")
    print("    - jesli zadna hipoteza < 0.3: klasa stopow wymaga osobnej teorii")


if __name__ == "__main__":
    main()
