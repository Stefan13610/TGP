#!/usr/bin/env python3
# =============================================================================
#  r17_sensitivity.py
# -----------------------------------------------------------------------------
#  Etap 9 - test wrazliwosci U24 na perturbacje Theta_D i N_F
#
#  Tezownienie: dane literaturowe Theta_D i N_F maja niepewnosci ~5-15%
#  (rozne zrodla daja rozne wartosci). Czy U24 coefs zmieniaja sie istotnie
#  gdy perturbujemy input?
#
#  Testy:
#    (A) Perturbate Theta_D o +- 10% (symetrycznie) i sprawdzic zmiane coefs
#    (B) Perturbate N_F o +- 15% (wieksza niepewnosc w DFT DOS)
#    (C) Monte Carlo: randomize obie o gaussian noise N(0, sigma) i zbierac
#        rozklad coefs
#    (D) Specific outlier test: sprawdzic jak zachowania "outliery" (Bi, Be, Ca)
#        zmienia coefs (leave-outlier-out vs full)
#
#  Oczekiwany wynik: coefs powinny byc stabilne bo LOO overfit ~1.68x
#  oznacza ze model nie jest uber-dopasowany do pojedynczych punktow.
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


def build_U24(logTh, logNF, v_arr, delta_s, delta_sp):
    X = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                          logTh, logNF, v_arr,
                          v_arr*delta_s, v_arr*delta_sp,
                          delta_s*logNF, delta_sp*logNF])
    return X


def fit_U24(logTh, logNF, v_arr, delta_s, delta_sp, logR):
    X = build_U24(logTh, logNF, v_arr, delta_s, delta_sp)
    coefs, _, _, _ = np.linalg.lstsq(X, logR, rcond=None)
    pred = X @ coefs
    rms = float(np.sqrt(np.mean((pred - logR)**2)))
    return coefs, rms


def rho_bg_predict(R, rho0, Th, T):
    return rho0 + R * (T/Th)**5 * r01.J5(Th/T)


def loo_profile_rms(logTh, logNF, v_arr, delta_s, delta_sp, logR, rows):
    N = len(rows)
    X = build_U24(logTh, logNF, v_arr, delta_s, delta_sp)
    prof = np.zeros(N)
    for i in range(N):
        mask = np.ones(N, bool); mask[i] = False
        coefs_i, _, _, _ = np.linalg.lstsq(X[mask], logR[mask], rcond=None)
        R_pred = 10**(X[i] @ coefs_i)
        Th = rows[i]["Theta_D"]; rho0 = rows[i]["rho_0"]
        pts = []
        for T in r00.T_POINTS:
            rho_p = rho_bg_predict(R_pred, rho0, Th, T)
            rho_o = rows[i]["rho_obs"][T]
            pts.append(np.log10(rho_p) - np.log10(rho_o))
        prof[i] = np.sqrt(np.mean(np.array(pts)**2))
    return float(np.mean(prof))


def main():
    print("=" * 100)
    print("  r17_sensitivity.py - perturbacje Theta_D / N_F i stabilnosc U24")
    print("=" * 100)

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

    Theta_arr = np.array([r["Theta_D"] for r in rows])
    NF_arr = np.array([r["N_F"] for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    cls_arr = np.array([r["cls"] for r in rows])
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)

    names = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)",
             "d(v)", "ds(v*s)", "dsp(v*sp)", "cs(s*logNF)", "csp(sp*logNF)"]

    # ------------------------------------------------------------
    # Baseline
    # ------------------------------------------------------------
    logTh_base = np.log10(Theta_arr)
    logNF_base = np.log10(NF_arr + 1e-30)
    c_base, rms_base = fit_U24(logTh_base, logNF_base, v_arr, delta_s, delta_sp, logR)
    prof_base = loo_profile_rms(logTh_base, logNF_base, v_arr, delta_s, delta_sp, logR, rows)
    print("\n  Baseline (no perturbation): RMS = {:.4f}, LOO prof RMS = {:.4f}".format(rms_base, prof_base))

    # ------------------------------------------------------------
    # (A) Theta_D perturbation: +10%, -10%
    # ------------------------------------------------------------
    print("\n" + "-" * 100)
    print("  (A) Theta_D perturbation (+- 10%)")
    print("-" * 100)
    print("  {:<18s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}".format(
        "param", "base", "Th+10%", "Th-10%", "mean_dev", "rel_dev%"))
    perts = []
    for scale in [1.10, 0.90]:
        logTh_p = np.log10(Theta_arr * scale)
        c_p, rms_p = fit_U24(logTh_p, logNF_base, v_arr, delta_s, delta_sp, logR)
        perts.append(c_p)
    c_plus, c_minus = perts
    for i, nm in enumerate(names):
        base = c_base[i]
        plus = c_plus[i]
        minus = c_minus[i]
        dev = max(abs(plus - base), abs(minus - base))
        rel = dev / max(abs(base), 0.01) * 100
        print("  {:<18s}  {:>+10.4f}  {:>+10.4f}  {:>+10.4f}  {:>10.4f}  {:>9.1f}%".format(
            nm, base, plus, minus, dev, rel))

    # ------------------------------------------------------------
    # (B) N_F perturbation: +15%, -15%
    # ------------------------------------------------------------
    print("\n" + "-" * 100)
    print("  (B) N_F perturbation (+- 15%)")
    print("-" * 100)
    print("  {:<18s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}".format(
        "param", "base", "NF+15%", "NF-15%", "mean_dev", "rel_dev%"))
    perts = []
    for scale in [1.15, 0.85]:
        logNF_p = np.log10(NF_arr * scale + 1e-30)
        c_p, rms_p = fit_U24(logTh_base, logNF_p, v_arr, delta_s, delta_sp, logR)
        perts.append(c_p)
    c_plus, c_minus = perts
    for i, nm in enumerate(names):
        base = c_base[i]
        plus = c_plus[i]
        minus = c_minus[i]
        dev = max(abs(plus - base), abs(minus - base))
        rel = dev / max(abs(base), 0.01) * 100
        print("  {:<18s}  {:>+10.4f}  {:>+10.4f}  {:>+10.4f}  {:>10.4f}  {:>9.1f}%".format(
            nm, base, plus, minus, dev, rel))

    # ------------------------------------------------------------
    # (C) Monte Carlo: gaussian noise sigma=0.10 on log Theta, 0.15 on log N_F
    # ------------------------------------------------------------
    print("\n" + "-" * 100)
    print("  (C) Monte Carlo (N=500) - gaussian noise sigma_Th=0.04, sigma_NF=0.06")
    print("      (odchylenia w log10: 0.04 = 10% w rzeczywistym scale)")
    print("-" * 100)
    rng = np.random.default_rng(42)
    N_mc = 500
    sigma_logTh = 0.04
    sigma_logNF = 0.06
    all_coefs = np.zeros((N_mc, 10))
    for k in range(N_mc):
        noise_Th = rng.normal(0, sigma_logTh, N)
        noise_NF = rng.normal(0, sigma_logNF, N)
        logTh_p = logTh_base + noise_Th
        logNF_p = logNF_base + noise_NF
        c_p, _ = fit_U24(logTh_p, logNF_p, v_arr, delta_s, delta_sp, logR)
        all_coefs[k] = c_p
    mc_mean = np.mean(all_coefs, axis=0)
    mc_std = np.std(all_coefs, axis=0)
    print("  {:<18s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}".format(
        "param", "base", "MC_mean", "MC_std", "|dev|", "|dev|/|base|"))
    for i, nm in enumerate(names):
        base = c_base[i]
        mn = mc_mean[i]
        st = mc_std[i]
        dev = abs(mn - base)
        rel = st / max(abs(base), 0.01) * 100
        print("  {:<18s}  {:>+10.4f}  {:>+10.4f}  {:>10.4f}  {:>10.4f}  {:>9.1f}%".format(
            nm, base, mn, st, dev, rel))

    # ------------------------------------------------------------
    # (D) Outlier-removal test
    # ------------------------------------------------------------
    print("\n" + "-" * 100)
    print("  (D) Leave-outlier-out test: Bi, Be, Ca, Sr, Tc")
    print("-" * 100)
    outliers = ["Bi", "Be", "Ca", "Sr", "Tc"]
    mask_ok = np.array([r["name"] not in outliers for r in rows])
    c_clean, rms_clean = fit_U24(
        logTh_base[mask_ok], logNF_base[mask_ok], v_arr[mask_ok],
        delta_s[mask_ok], delta_sp[mask_ok], logR[mask_ok])
    print("  N (without outliers) = {}".format(mask_ok.sum()))
    print("  {:<18s}  {:>10s}  {:>10s}  {:>10s}".format(
        "param", "base (N=37)", "clean (N=32)", "shift"))
    for i, nm in enumerate(names):
        base = c_base[i]
        clean = c_clean[i]
        shift = clean - base
        print("  {:<18s}  {:>+10.4f}  {:>+10.4f}  {:>+10.4f}".format(
            nm, base, clean, shift))

    # ------------------------------------------------------------
    # VERDICT
    # ------------------------------------------------------------
    print("\n" + "=" * 100)
    print("  VERDICT r17")
    print("=" * 100)
    print()

    # Kluczowe coefs: b(logTh) i d(v) - najbardziej robust w bootstrap
    key_coef_names = ["b(logTh)", "d(v)", "as(s-d)"]
    key_ok = True
    print("  Stabilnosc kluczowych parametrow (rel_dev < 30%):")
    for nm in key_coef_names:
        i = names.index(nm)
        # MC std
        rel_mc = mc_std[i] / max(abs(c_base[i]), 0.01) * 100
        # Max perturbation from (A)
        mark = "OK" if rel_mc < 30 else "UWAGA"
        print("    {:<18s}: MC std/|base| = {:.1f}%  [{}]".format(nm, rel_mc, mark))
        if rel_mc > 30: key_ok = False

    print()
    # Check whether b(logTh) still positive (i.e., Bloch-Mott direction)
    i_b = names.index("b(logTh)")
    n_negative_b = int((all_coefs[:, i_b] < 0).sum())
    print("  b(logTh) > 0 w {} / {} prob MC ({}%)".format(
        N_mc - n_negative_b, N_mc, (N_mc - n_negative_b)/N_mc*100))

    i_d = names.index("d(v)")
    n_pos_d = int((all_coefs[:, i_d] > 0).sum())
    print("  d(v) < 0 w {} / {} prob MC ({}%)".format(
        N_mc - n_pos_d, N_mc, (N_mc - n_pos_d)/N_mc*100))

    print()
    if key_ok:
        print("  [OK] Kluczowe parametry U24 stabilne na perturbacje Theta/N_F w granicach 10-15%.")
        print("       Formula JEST robustna, niepewnosci literaturowe nie psuja wynikow.")
    else:
        print("  [UWAGA] Ktoras z kluczowych parametrow niestabilna. Sprawdzic szczegoly.")

    print()
    print("  Wniosek: ostrzal sensitivity ZAKONCZONY. Gotowe do r18 (Nordheim).")


if __name__ == "__main__":
    main()
