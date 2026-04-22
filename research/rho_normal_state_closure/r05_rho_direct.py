#!/usr/bin/env python3
# =============================================================================
#  r05_rho_direct.py
# -----------------------------------------------------------------------------
#  r04 pokazal: F_TGP NIE przewiduje lambda_ep (RMS>0.3, r=-0.06). Unifikacja
#  przez lambda_TGP nie dziala.
#
#  Pivot: testujemy BEZPOSREDNIO czy kombinacja (F_TGP, Theta_D, N_F) przewiduje
#  AMPLITUDE BG (R) oraz PELNY profil rho(T) na wszystkich 4 T-points.
#
#  Hipotezy:
#     D0:  R = C * F_TGP                          (baseline z r02, RMS~0.34)
#     D1:  R = C * F_TGP * Theta_D                (dodaje energy scale)
#     D2:  R = C * F_TGP * Theta_D^p              (fit p)
#     D3:  R = C * F_TGP^a * Theta_D^b            (full power law)
#     D4:  R = C * F_TGP^a * Theta_D^b * N_F^c    (N_F jawnie)
#     D5:  R = C * Theta_D * (lam_lit)            (classic BG, "cheat" z lit)
#          -- porownanie jako upper bound dla RMS moze (dostosuj tylko C)
#
#  Dla najlepszej formy:
#     - RMS_log(R_TGP, R_fit) - celujemy <0.25
#     - Test pelnego profilu: predykcja rho(T_i) dla i=1..4 z wspolnym rho_0
#       i TGP R, oceń RMS_log(rho_pred, rho_obs) per-material.
#     - LOO dla stabilnosci.
#
#  Jesli D5 >> D4: brakujace lambda jest dominujace - TGP nie wystarcza
#  Jesli D4 ~ D5: TGP wystarcza (z F_TGP, Theta_D, N_F odzyskujemy lambda·Theta)
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import pearsonr
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

spec4 = importlib.util.spec_from_file_location(
    "r04", os.path.join(os.path.dirname(__file__), "r04_lambda_tgp.py"))
r04 = importlib.util.module_from_spec(spec4); sys.modules["r04"] = r04
spec4.loader.exec_module(r04)


def main():
    print("=" * 92)
    print("  r05_rho_direct.py - R_BG = f(F_TGP, Theta_D, N_F) direct, bypass lambda")
    print("=" * 92)

    data = r00.load_dataset()

    # Fit R_BG per material (z r01)
    rows = []
    for d in data:
        R, rms, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        z = r02.COORD[d["name"]]
        cls = r02.ORB_CLASS[d["name"]]
        F_TGP = r02.tgp_factor_rho(d, z, use_g=False)
        lam_lit = r04.LAMBDA_EP_LIT.get(d["name"], None)
        rows.append({
            "name":    d["name"],
            "cls":     cls,
            "R":       R,
            "Theta_D": d["Theta_D"],
            "N_F":     d["N_F"],
            "F_TGP":   F_TGP,
            "a_latt":  d["a_latt"],
            "rho_0":   d["rho_0"],
            "rho_obs": d["rho"],
            "lam_lit": lam_lit,
            "bg_rms":  rms,
        })

    # Array shortcuts
    F_arr = np.array([r["F_TGP"] for r in rows])
    Th_arr = np.array([r["Theta_D"] for r in rows])
    NF_arr = np.array([r["N_F"] for r in rows])
    R_arr = np.array([r["R"] for r in rows])
    lam_arr = np.array([r["lam_lit"] for r in rows])
    logR = np.log10(R_arr)
    logF = np.log10(F_arr + 1e-30)
    logTh = np.log10(Th_arr)
    logNF = np.log10(NF_arr + 1e-30)

    # -------- D0: R = C * F_TGP (baseline) --------
    print("\n" + "-" * 92)
    print("  D0: R = C * F_TGP  (baseline z r02)")
    print("-" * 92)
    def obj_D0(logC):
        pred = 10**logC * F_arr
        return np.sqrt(np.mean((np.log10(pred) - logR)**2))
    res_D0 = minimize_scalar(obj_D0, bounds=(-5, 10), method="bounded")
    rms_D0 = res_D0.fun
    C_D0 = 10**res_D0.x
    r_D0, _ = pearsonr(logF, logR)
    print("    C = {:.3e}, RMS_log = {:.4f}, r = {:.3f}".format(
        C_D0, rms_D0, r_D0))

    # -------- D1: R = C * F_TGP * Theta_D --------
    print("\n" + "-" * 92)
    print("  D1: R = C * F_TGP * Theta_D  (Ziman high-T: rho_inf ~ lambda*Theta_D)")
    print("-" * 92)
    FTh = F_arr * Th_arr
    def obj_D1(logC):
        pred = 10**logC * FTh
        return np.sqrt(np.mean((np.log10(pred) - logR)**2))
    res_D1 = minimize_scalar(obj_D1, bounds=(-10, 10), method="bounded")
    rms_D1 = res_D1.fun
    C_D1 = 10**res_D1.x
    r_D1, p_D1 = pearsonr(np.log10(FTh), logR)
    print("    C = {:.3e}, RMS_log = {:.4f}, r = {:.3f}, p = {:.4f}".format(
        C_D1, rms_D1, r_D1, p_D1))

    # -------- D2: R = C * F_TGP * Theta_D^p (fit p) --------
    print("\n" + "-" * 92)
    print("  D2: R = C * F_TGP * Theta_D^p  (fit p)")
    print("-" * 92)
    def obj_D2(params):
        logC, p = params
        pred = 10**logC * F_arr * Th_arr**p
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - logR)**2))
    best_D2 = (1e9, None)
    for lc0 in [-5, 0, 5]:
        for p0 in [0, 1, 2]:
            r_ = minimize(obj_D2, [lc0, p0], method="Nelder-Mead",
                          options={"xatol": 1e-5, "fatol": 1e-7})
            if r_.fun < best_D2[0]:
                best_D2 = (r_.fun, r_.x)
    logC_D2, p_D2 = best_D2[1]
    rms_D2 = best_D2[0]
    print("    C = {:.3e}, p = {:.3f}, RMS_log = {:.4f}".format(
        10**logC_D2, p_D2, rms_D2))

    # -------- D3: R = C * F_TGP^a * Theta_D^b (full power-law) --------
    print("\n" + "-" * 92)
    print("  D3: R = C * F_TGP^a * Theta_D^b  (full power-law)")
    print("-" * 92)
    # multilinear regression: logR = logC + a*logF + b*logTh
    X = np.column_stack([np.ones_like(logF), logF, logTh])
    coefs, _, _, _ = np.linalg.lstsq(X, logR, rcond=None)
    logC_D3, a_D3, b_D3 = coefs
    pred_D3 = 10**logC_D3 * F_arr**a_D3 * Th_arr**b_D3
    rms_D3 = np.sqrt(np.mean((np.log10(pred_D3) - logR)**2))
    print("    C = {:.3e}, a = {:.3f}, b = {:.3f}, RMS_log = {:.4f}".format(
        10**logC_D3, a_D3, b_D3, rms_D3))

    # -------- D4: R = C * F_TGP^a * Theta_D^b * N_F^c --------
    print("\n" + "-" * 92)
    print("  D4: R = C * F_TGP^a * Theta_D^b * N_F^c  (pelna 3-param)")
    print("-" * 92)
    X4 = np.column_stack([np.ones_like(logF), logF, logTh, logNF])
    coefs4, _, _, _ = np.linalg.lstsq(X4, logR, rcond=None)
    logC_D4, a_D4, b_D4, c_D4 = coefs4
    pred_D4 = 10**logC_D4 * F_arr**a_D4 * Th_arr**b_D4 * NF_arr**c_D4
    rms_D4 = np.sqrt(np.mean((np.log10(pred_D4) - logR)**2))
    print("    C = {:.3e}, a = {:.3f}, b = {:.3f}, c = {:.3f}".format(
        10**logC_D4, a_D4, b_D4, c_D4))
    print("    RMS_log = {:.4f}".format(rms_D4))

    # -------- D5: R = C * lambda_lit * Theta_D (with lit lambda) ------------
    print("\n" + "-" * 92)
    print("  D5: R = C * lambda_lit * Theta_D  (UPPER BOUND - uses lit lambda)")
    print("-" * 92)
    lamTh = lam_arr * Th_arr
    def obj_D5(logC):
        pred = 10**logC * lamTh
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - logR)**2))
    res_D5 = minimize_scalar(obj_D5, bounds=(-10, 10), method="bounded")
    rms_D5 = res_D5.fun
    C_D5 = 10**res_D5.x
    r_D5, p_D5 = pearsonr(np.log10(lamTh), logR)
    print("    C = {:.3e}, RMS_log = {:.4f}, r = {:.3f}, p = {:.4f}".format(
        C_D5, rms_D5, r_D5, p_D5))

    # -------- D6: R = C * lambda_lit^a * Theta_D^b (fit ich wag) ------------
    print("\n" + "-" * 92)
    print("  D6: R = C * lambda_lit^a * Theta_D^b  (UPPER BOUND full)")
    print("-" * 92)
    loglam = np.log10(lam_arr + 1e-30)
    X6 = np.column_stack([np.ones_like(logF), loglam, logTh])
    coefs6, _, _, _ = np.linalg.lstsq(X6, logR, rcond=None)
    logC_D6, a_D6, b_D6 = coefs6
    pred_D6 = 10**logC_D6 * lam_arr**a_D6 * Th_arr**b_D6
    rms_D6 = np.sqrt(np.mean((np.log10(pred_D6) - logR)**2))
    print("    C = {:.3e}, a = {:.3f}, b = {:.3f}, RMS_log = {:.4f}".format(
        10**logC_D6, a_D6, b_D6, rms_D6))

    # -------- D7: R = C * Theta_D^p + N_F^q (no F_TGP) ----------------------
    print("\n" + "-" * 92)
    print("  D7: R = C * Theta_D^p * N_F^q  (bez F_TGP - kontrola)")
    print("-" * 92)
    X7 = np.column_stack([np.ones_like(logTh), logTh, logNF])
    coefs7, _, _, _ = np.linalg.lstsq(X7, logR, rcond=None)
    logC_D7, p_D7, q_D7 = coefs7
    pred_D7 = 10**logC_D7 * Th_arr**p_D7 * NF_arr**q_D7
    rms_D7 = np.sqrt(np.mean((np.log10(pred_D7) - logR)**2))
    print("    C = {:.3e}, p = {:.3f}, q = {:.3f}, RMS_log = {:.4f}".format(
        10**logC_D7, p_D7, q_D7, rms_D7))

    # -------- Summary --------
    print("\n" + "=" * 92)
    print("  PODSUMOWANIE FITOW R_BG")
    print("=" * 92)
    fits = [
        ("D0: C * F_TGP",              rms_D0, 2),
        ("D1: C * F_TGP * Theta_D",    rms_D1, 2),
        ("D2: C * F_TGP * Theta_D^p",  rms_D2, 2),
        ("D3: C * F_TGP^a * Th^b",     rms_D3, 3),
        ("D4: C * F^a * Th^b * NF^c",  rms_D4, 4),
        ("D5: C * lam_lit * Th",       rms_D5, 2),
        ("D6: C * lam^a * Th^b",       rms_D6, 3),
        ("D7: C * Th^p * NF^q",        rms_D7, 3),
    ]
    print("\n  {:<30s} {:>10s} {:>8s}".format("Form", "RMS_log", "params"))
    for name, rms, p in fits:
        print("  {:<30s} {:>10.4f} {:>8d}".format(name, rms, p))

    best = min(fits, key=lambda x: x[1])
    print("\n  Best: {} (RMS = {:.4f})".format(best[0], best[1]))

    # -------- Per-material comparison dla D4 (full TGP model) vs D5 (lit upper bound)
    print("\n" + "-" * 92)
    print("  D4 vs D5: ile tracimy bo nie znamy lambda?")
    print("-" * 92)
    print("  {:<5s} {:>4s} {:>10s} {:>10s} {:>10s} {:>8s} {:>8s}".format(
        "name", "cls", "R_obs", "R_D4", "R_D5", "dlog D4", "dlog D5"))
    resid_D4 = []
    resid_D5 = []
    for i, r in enumerate(rows):
        Rp4 = pred_D4[i]
        Rp5 = C_D5 * lam_arr[i] * Th_arr[i]
        d4 = np.log10(Rp4) - np.log10(r["R"])
        d5 = np.log10(Rp5) - np.log10(r["R"])
        resid_D4.append(d4); resid_D5.append(d5)
        print("  {:<5s} {:>4s} {:>10.3f} {:>10.3f} {:>10.3f} {:>+8.3f} {:>+8.3f}".format(
            r["name"], r["cls"], r["R"], Rp4, Rp5, d4, d5))

    # -------- Full rho(T_i) profile test: uses D4 predicted R --------
    print("\n" + "-" * 92)
    print("  PROFIL rho(T): fit BG z R_TGP (D4), porownaj rho_obs(T_i) na 4 T")
    print("-" * 92)
    print("  {:<5s} {:>4s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>8s}".format(
        "name", "cls", "r(77)obs", "r(77)TGP", "r(295)obs", "r(295)TGP",
        "r(1000)", "RMS_log"))
    profile_rms = []
    for i, r in enumerate(rows):
        R_pred = pred_D4[i]
        Th = r["Theta_D"]
        rho0 = r["rho_0"]
        rms_pts = []
        rp_print = []
        for T in r00.T_POINTS:
            rho_pred = r01.rho_bg(T, rho0, R_pred, Th)
            rho_obs = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_pred) - np.log10(rho_obs))
            rp_print.append((T, rho_pred))
        rms_mat = np.sqrt(np.mean(np.array(rms_pts)**2))
        profile_rms.append(rms_mat)
        # print tylko 77, 295, 1000
        o77 = r["rho_obs"][77.0]; p77 = rp_print[0][1]
        o295 = r["rho_obs"][295.0]; p295 = rp_print[1][1]
        o1000 = r["rho_obs"][1000.0]; p1000 = rp_print[3][1]
        print("  {:<5s} {:>4s} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>8.4f}".format(
            r["name"], r["cls"], o77, p77, o295, p295, o1000, rms_mat))
    print("\n  Per-material profile RMS_log statistics:")
    print("    mean = {:.4f}".format(np.mean(profile_rms)))
    print("    median = {:.4f}".format(np.median(profile_rms)))
    print("    max = {:.4f} ({})".format(
        max(profile_rms), rows[np.argmax(profile_rms)]["name"]))
    print("    CoV = {:.1%}".format(np.std(profile_rms)/np.mean(profile_rms)))

    # -------- LOO dla D4 --------
    print("\n" + "-" * 92)
    print("  LOO dla D4: stabilnosc wykladnikow (a, b, c)?")
    print("-" * 92)
    abc_loo = []
    for i in range(len(rows)):
        mask = np.ones(len(rows), dtype=bool)
        mask[i] = False
        X_sub = X4[mask]
        logR_sub = logR[mask]
        coef_sub, _, _, _ = np.linalg.lstsq(X_sub, logR_sub, rcond=None)
        abc_loo.append(coef_sub[1:])  # (a, b, c)
    abc_arr = np.array(abc_loo)
    print("    a LOO: mean {:.3f} +/- {:.3f}  (range {:.3f}..{:.3f})".format(
        abc_arr[:, 0].mean(), abc_arr[:, 0].std(),
        abc_arr[:, 0].min(), abc_arr[:, 0].max()))
    print("    b LOO: mean {:.3f} +/- {:.3f}  (range {:.3f}..{:.3f})".format(
        abc_arr[:, 1].mean(), abc_arr[:, 1].std(),
        abc_arr[:, 1].min(), abc_arr[:, 1].max()))
    print("    c LOO: mean {:.3f} +/- {:.3f}  (range {:.3f}..{:.3f})".format(
        abc_arr[:, 2].mean(), abc_arr[:, 2].std(),
        abc_arr[:, 2].min(), abc_arr[:, 2].max()))
    # sensitivity
    max_d = [(abs(abc_loo[i][j] - coefs4[j+1])) for i in range(len(rows))
             for j in range(3)]

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r05")
    print("=" * 92)
    gap_lit = rms_D4 - rms_D6  # D4 = TGP-only, D6 = lambda-lit-based
    gap_D4_D5 = rms_D4 - rms_D5  # D5 = lam*Theta
    print("\n  Gap TGP-only (D4) vs lit-based (D6): {:+.4f}".format(gap_lit))
    print("  Gap TGP-only (D4) vs D5 (lam*Th): {:+.4f}".format(gap_D4_D5))
    if gap_lit < 0.05:
        verdict = "TGP z Theta_D (D4) osiaga upper bound lit (D6). MOCNY!"
    elif gap_lit < 0.15:
        verdict = "TGP ma moderat niedomiar vs lit - brakuje drobnych efektow."
    else:
        verdict = "TGP daje znaczny niedomiar - wymaga bezposredniego lambda."
    print("  " + verdict)

    print("""
  Nastepne kroki:
    - Analiza: ktore materialy dominuja residuum D4? (Cu/Ag/Au? Fe/Ni?)
    - Test na rozszerzonym dataset (Mo, Ta, W, Be, Ti, Bi, Sb) - r06
    - Rozpoznanie czy D4 wspolczynniki (a, b, c) sa FIZYCZNE czy tylko empiryczne:
         * a (F_TGP) ~ ?
         * b (Theta_D) ~ ? (oczekujemy ~1 dla Ziman)
         * c (N_F) ~ ? (oczekujemy ~0 bo F_TGP juz tego nie wymaga)
""")


if __name__ == "__main__":
    main()
