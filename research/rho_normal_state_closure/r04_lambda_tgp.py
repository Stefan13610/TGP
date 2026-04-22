#!/usr/bin/env python3
# =============================================================================
#  r04_lambda_tgp.py
# -----------------------------------------------------------------------------
#  Etap 2 (unifikacja): czy TGP Eq.5 struktura przewiduje lambda_ep?
#
#  Kluczowa hipoteza po r03:
#     R_BG ~ lambda_ep * Theta_D  (Pearson r = 0.79, p = 0.0004)
#     gdzie R_BG to amplituda BG a lambda_ep z literatury (Allen-Mitrovic).
#
#  r04 pyta bezposrednio: CZY LAMBDA_EP = C_lam * F_TGP?
#  Jesli tak, to TGP jest uniwersalnym generatorem lambda_ep, a ten z kolei
#  napedza zarowno T_c (McMillan) jak i rho(T) (Bloch-Gruneisen).
#
#  F_TGP = k_d(z) * A_orb^2 * M_gauss(a) * Lambda_eff(omega)
#
#  Uwaga: NIE uzywamy g(N_F) z P7.12 - r02 pokazal ze jest mode-specific
#  (dziala dla SC gap, nie dla transportu). Wiec:
#
#     lambda_ep^TGP = C_lambda * F_TGP   (hipoteza A0)
#     lambda_ep^TGP = C_lambda * F_TGP * N_F^q   (hipoteza A1, fit q)
#     lambda_ep^TGP = C_lambda * F_TGP^alpha   (hipoteza A2, fit alpha)
#
#  Testy:
#     - RMS_log(lambda_TGP, lambda_lit) - celujemy w <0.15
#     - Pearson r - celujemy w >0.90
#     - Residuum per klasa (noble/sp/d/FM/Stoner)
#     - LOO dla robust check
#
#  Dalej (r05): sprawdz czy lambda_TGP napedza rho(T) przez pelny BG na wszystkie
#  T_i simultaneously (nie tylko amplitude R).
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import pearsonr
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "r00", os.path.join(os.path.dirname(__file__), "r00_dataset.py"))
r00 = importlib.util.module_from_spec(spec); sys.modules["r00"] = r00
spec.loader.exec_module(r00)

spec2 = importlib.util.spec_from_file_location(
    "r02", os.path.join(os.path.dirname(__file__), "r02_tgp_formula.py"))
r02 = importlib.util.module_from_spec(spec2); sys.modules["r02"] = r02
spec2.loader.exec_module(r02)


# -----------------------------------------------------------------------------
# Literatura lambda_ep (Allen-Mitrovic 1982, Grimvall 1981, DFT reviews)
# -----------------------------------------------------------------------------
# Wartosci z aktualnej literatury teorii - najczesciej McMillan-Allen inversion
# albo direct DFT.
LAMBDA_EP_LIT = {
    # noble (slabe s-p)
    "Cu": 0.13,   # Grimvall; Allen 1987
    "Ag": 0.12,   # Grimvall
    "Au": 0.15,   # Grimvall
    # sp metals (SC)
    "Al": 0.43,   # Allen-Mitrovic 1982
    "Pb": 1.55,   # Allen-Mitrovic (strong coupling)
    "Sn": 0.72,   # Allen-Mitrovic
    # bcc TM SC
    "Nb": 1.12,   # Savrasov DFT, Allen
    "V":  0.82,   # Allen-Mitrovic
    # FM (d elektrony maja duze sprzezenie ale tez spin-fluct)
    # UWAGA: w FM lambda_ep ma sens tylko w parafazie (T>T_C) albo para-orbital
    # Te wartosci z Savrasov dla paramagnonowego limita (bez spin-fluctuations)
    "Fe": 0.25,   # DFT bare (bez spin-fluct)
    "Ni": 0.18,   # DFT bare
    # Stoner paramagnetyki
    "Pt": 0.66,   # Allen, Savrasov
    "Pd": 0.35,   # Allen, Savrasov
    # hcp sp
    "Cd": 0.40,   # Allen-Mitrovic
    "Zn": 0.38,   # Allen-Mitrovic
    "Mg": 0.35,   # Allen-Mitrovic
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 92)
    print("  r04_lambda_tgp.py - TGP Eq.5 struktura vs lambda_ep(literatura)")
    print("=" * 92)
    print("""
  Cel: Czy lambda_ep = C_lambda * F_TGP(k_d, A_orb, M_gauss, Lambda_eff)?
  Jesli tak, unifikacja: lambda_TGP -> T_c (McMillan) i lambda_TGP -> rho(T) (BG).
""")

    data = r00.load_dataset()

    # Zbuduj tabele: dla kazdego materialu F_TGP (no g) i lambda_ep z lit
    rows = []
    for d in data:
        name = d["name"]
        if name not in LAMBDA_EP_LIT:
            continue
        z = r02.COORD[name]
        cls = r02.ORB_CLASS[name]
        F_TGP = r02.tgp_factor_rho(d, z, use_g=False)
        lam_lit = LAMBDA_EP_LIT[name]
        rows.append({
            "name":    name,
            "cls":     cls,
            "z":       z,
            "Theta_D": d["Theta_D"],
            "a_latt":  d["a_latt"],
            "N_F":     d["N_F"],
            "Z":       d["Z"],
            "F_TGP":   F_TGP,
            "lam_lit": lam_lit,
        })

    # -------- Tabela diagnostyczna --------
    print("-" * 92)
    print("  Input data (F_TGP = k_d * A_orb^2 * M_gauss * Lambda_eff)")
    print("-" * 92)
    print("  {:<5s} {:>4s} {:>3s} {:>7s} {:>7s} {:>6s} {:>3s} {:>12s} {:>9s}".format(
        "name", "cls", "z", "Th_D", "a[A]", "N_F", "Z", "F_TGP", "lam_lit"))
    for r in rows:
        print("  {:<5s} {:>4s} {:>3d} {:>7.0f} {:>7.3f} {:>6.2f} {:>3d} {:>12.5e} {:>9.3f}".format(
            r["name"], r["cls"], r["z"], r["Theta_D"], r["a_latt"],
            r["N_F"], r["Z"], r["F_TGP"], r["lam_lit"]))

    # -------- (A0) lambda = C * F_TGP  (simplest) --------
    print("\n" + "-" * 92)
    print("  (A0) lambda_ep = C_lambda * F_TGP  (proporcjonalnosc)")
    print("-" * 92)
    F_arr = np.array([r["F_TGP"] for r in rows])
    lam_arr = np.array([r["lam_lit"] for r in rows])
    logF = np.log10(F_arr + 1e-30)
    loglam = np.log10(lam_arr)

    # Linear in log: loglam = logC + 1 * logF
    # Fit C fixing slope=1
    def obj_A0(logC):
        pred = 10**logC * F_arr
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - loglam)**2))
    res_A0 = minimize_scalar(obj_A0, bounds=(-10, 10), method="bounded")
    C_A0 = 10**res_A0.x
    rms_A0 = res_A0.fun
    r_A0, p_A0 = pearsonr(logF, loglam)
    print("    C_lambda  = {:.3e}".format(C_A0))
    print("    RMS_log   = {:.4f}".format(rms_A0))
    print("    Pearson r(log F, log lam) = {:.3f},  p = {:.5f}".format(
        r_A0, p_A0))

    # -------- (A1) lambda = C * F_TGP * N_F^q --------
    print("\n" + "-" * 92)
    print("  (A1) lambda_ep = C * F_TGP * N_F^q  (DOS-modulowana)")
    print("-" * 92)
    NF_arr = np.array([r["N_F"] for r in rows])

    def obj_A1(params):
        logC, q = params
        pred = 10**logC * F_arr * NF_arr**q
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - loglam)**2))

    best_A1 = (1e9, None)
    for lc0 in [-3, -1, 1, 3]:
        for q0 in [-1, 0, 0.5, 1, 2]:
            r_ = minimize(obj_A1, [lc0, q0], method="Nelder-Mead",
                          options={"xatol": 1e-5, "fatol": 1e-7})
            if r_.fun < best_A1[0]:
                best_A1 = (r_.fun, r_.x)
    logC_A1, q_A1 = best_A1[1]
    rms_A1 = best_A1[0]
    print("    C = {:.3e},  q = {:.3f}".format(10**logC_A1, q_A1))
    print("    RMS_log = {:.4f}".format(rms_A1))

    # -------- (A2) lambda = C * F_TGP^alpha (power law) --------
    print("\n" + "-" * 92)
    print("  (A2) lambda_ep = C * F_TGP^alpha  (power law, fit alpha)")
    print("-" * 92)
    slope, intercept = np.polyfit(logF, loglam, 1)
    C_A2 = 10**intercept
    pred_A2 = C_A2 * F_arr**slope
    rms_A2 = np.sqrt(np.mean((np.log10(pred_A2) - loglam)**2))
    print("    alpha = {:.3f},  C = {:.3e}".format(slope, C_A2))
    print("    RMS_log = {:.4f}".format(rms_A2))

    # -------- (A3) lambda = C * F_TGP * N_F (fix q=1, Ziman-like) --------
    print("\n" + "-" * 92)
    print("  (A3) lambda_ep = C * F_TGP * N_F  (q = 1 fixed, fizyczny limit)")
    print("-" * 92)
    pred_FN = F_arr * NF_arr
    def obj_A3(logC):
        pred = 10**logC * pred_FN
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - loglam)**2))
    res_A3 = minimize_scalar(obj_A3, bounds=(-10, 10), method="bounded")
    C_A3 = 10**res_A3.x
    rms_A3 = res_A3.fun
    r_A3, _ = pearsonr(np.log10(pred_FN + 1e-30), loglam)
    print("    C = {:.3e}".format(C_A3))
    print("    RMS_log = {:.4f}".format(rms_A3))
    print("    Pearson r = {:.3f}".format(r_A3))

    # -------- (A4) Reassign Cu/Ag/Au to "s" class (filled d below E_F) ------
    print("\n" + "-" * 92)
    print("  (A4) Cu/Ag/Au reklasyfikacja 'sp' -> 's' (A_orb = -0.111)")
    print("-" * 92)
    # A_orb = -0.111 dla noble -> A_orb^2 = 0.01232 (vs 0.0428 dla sp)
    # factor reduction = 0.01232 / 0.0428 = 0.288
    F_noble_adj = F_arr.copy()
    for i, r in enumerate(rows):
        if r["name"] in ("Cu", "Ag", "Au"):
            F_noble_adj[i] *= (r02.A_ORB["s"]**2 / r02.A_ORB["sp"]**2)
    # Fit prop
    def obj_A4(logC):
        pred = 10**logC * F_noble_adj
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - loglam)**2))
    res_A4 = minimize_scalar(obj_A4, bounds=(-10, 10), method="bounded")
    C_A4 = 10**res_A4.x
    rms_A4 = res_A4.fun
    r_A4, p_A4 = pearsonr(np.log10(F_noble_adj + 1e-30), loglam)
    print("    C = {:.3e},  RMS_log = {:.4f}".format(C_A4, rms_A4))
    print("    Pearson r = {:.3f}, p = {:.4f}".format(r_A4, p_A4))

    # -------- (A5) Hopfield: lambda ~ F_TGP * N_F / Theta_D^2 ----------------
    print("\n" + "-" * 92)
    print("  (A5) Hopfield: lambda ~ N_F * <I^2> / (M * <omega^2>)")
    print("       Test: lambda = C * F_TGP * N_F / Theta_D^2  (dim-less)")
    print("-" * 92)
    Th_arr = np.array([r["Theta_D"] for r in rows])
    hop = F_arr * NF_arr / Th_arr**2
    def obj_A5(logC):
        pred = 10**logC * hop
        return np.sqrt(np.mean((np.log10(pred + 1e-30) - loglam)**2))
    res_A5 = minimize_scalar(obj_A5, bounds=(-20, 20), method="bounded")
    C_A5 = 10**res_A5.x
    rms_A5 = res_A5.fun
    r_A5, p_A5 = pearsonr(np.log10(hop + 1e-30), loglam)
    print("    C = {:.3e},  RMS_log = {:.4f}".format(C_A5, rms_A5))
    print("    Pearson r = {:.3f}, p = {:.4f}".format(r_A5, p_A5))

    # -------- (A6) Odwrot: lambda ~ 1/Theta_D ------------------------------
    print("\n" + "-" * 92)
    print("  (A6) Sanity: lambda_ep vs 1/Theta_D (Pb ma Theta_D=105, Cu=343)")
    print("-" * 92)
    inv_th = 1.0 / Th_arr
    r_inv, p_inv = pearsonr(np.log10(inv_th), loglam)
    print("    Pearson r(log(1/Th), log lam) = {:.3f}, p = {:.4f}".format(
        r_inv, p_inv))
    # Fit proportionality
    slope_inv, intercept_inv = np.polyfit(np.log10(inv_th), loglam, 1)
    C_inv = 10**intercept_inv
    pred_inv = C_inv * inv_th**slope_inv
    rms_inv = np.sqrt(np.mean((np.log10(pred_inv) - loglam)**2))
    print("    alpha = {:.3f}, C = {:.3e}, RMS_log = {:.4f}".format(
        slope_inv, C_inv, rms_inv))

    # -------- (A7) lambda ~ N_F alone --------
    print("\n" + "-" * 92)
    print("  (A7) Sanity: lambda_ep vs N_F alone (DOS dominuje?)")
    print("-" * 92)
    r_NF, p_NF = pearsonr(np.log10(NF_arr), loglam)
    slope_NF, intercept_NF = np.polyfit(np.log10(NF_arr), loglam, 1)
    print("    Pearson r(log NF, log lam) = {:.3f}, p = {:.4f}".format(
        r_NF, p_NF))
    print("    alpha_NF = {:.3f}".format(slope_NF))
    pred_NF_only = 10**intercept_NF * NF_arr**slope_NF
    rms_NF_only = np.sqrt(np.mean((np.log10(pred_NF_only) - loglam)**2))
    print("    RMS_log = {:.4f}".format(rms_NF_only))

    # -------- Per-material predykcja (best model) --------
    fits = [
        ("A0: C * F_TGP",           rms_A0, C_A0, None, None),
        ("A1: C * F_TGP * N_F^q",   rms_A1, 10**logC_A1, q_A1, None),
        ("A2: C * F_TGP^alpha",     rms_A2, C_A2, slope, None),
        ("A3: C * F_TGP * N_F",     rms_A3, C_A3, 1.0, None),
        ("A4: A0 with Cu/Ag/Au=s",  rms_A4, C_A4, None, None),
        ("A5: C * F * NF / Th^2",   rms_A5, C_A5, None, None),
        ("A6: C * (1/Th)^alpha",    rms_inv, C_inv, slope_inv, None),
        ("A7: C * N_F^alpha",       rms_NF_only, 10**intercept_NF, slope_NF, None),
    ]
    best_fit = min(fits, key=lambda x: x[1])
    print("\n" + "=" * 92)
    print("  Best model: {} (RMS = {:.4f})".format(best_fit[0], best_fit[1]))
    print("=" * 92)

    # Per-material comparison tabelka dla best modelu
    print("\n  Predykcja {}:".format(best_fit[0]))
    print("  {:<5s} {:>4s} {:>9s} {:>9s} {:>8s}".format(
        "name", "cls", "lam_lit", "lam_TGP", "dlog"))
    if best_fit[0].startswith("A0"):
        pred_best = C_A0 * F_arr
    elif best_fit[0].startswith("A1"):
        pred_best = 10**logC_A1 * F_arr * NF_arr**q_A1
    elif best_fit[0].startswith("A2"):
        pred_best = C_A2 * F_arr**slope
    elif best_fit[0].startswith("A3"):
        pred_best = C_A3 * F_arr * NF_arr
    elif best_fit[0].startswith("A4"):
        pred_best = C_A4 * F_noble_adj
    elif best_fit[0].startswith("A5"):
        pred_best = C_A5 * hop
    elif best_fit[0].startswith("A6"):
        pred_best = C_inv * inv_th**slope_inv
    elif best_fit[0].startswith("A7"):
        pred_best = 10**intercept_NF * NF_arr**slope_NF
    else:
        pred_best = C_A0 * F_arr

    dlogs = []
    for r, pred in zip(rows, pred_best):
        dlog = np.log10(pred + 1e-30) - np.log10(r["lam_lit"])
        dlogs.append(dlog)
        flag = " <- out" if abs(dlog) > 0.5 else ""
        print("  {:<5s} {:>4s} {:>9.3f} {:>9.3f} {:>+8.3f}{}".format(
            r["name"], r["cls"], r["lam_lit"], pred, dlog, flag))

    # Per-class residuum
    print("\n  Residuum per klasa (best model):")
    classes = {}
    for r, dlog in zip(rows, dlogs):
        classes.setdefault(r["cls"], []).append((r["name"], dlog))
    for cls, lst in classes.items():
        vals = [v for _, v in lst]
        print("    {:<5s} (N={}):  mean dlog = {:+.3f},  RMS = {:.3f}".format(
            cls, len(lst), np.mean(vals), np.sqrt(np.mean(np.asarray(vals)**2))))
        print("           vals: {}".format(
            ", ".join("{}{:+.2f}".format(n, v) for n, v in lst)))

    # -------- LOO dla A0 (najprostszy test stabilnosci) --------
    print("\n" + "-" * 92)
    print("  LOO (Leave-One-Out) dla A0: C_lambda stabilne?")
    print("-" * 92)
    C_loo = []
    for i in range(len(rows)):
        mask = np.ones(len(rows), dtype=bool)
        mask[i] = False
        F_sub = F_arr[mask]
        lam_sub = lam_arr[mask]
        def _obj(logC):
            pred = 10**logC * F_sub
            return np.sqrt(np.mean((np.log10(pred) - np.log10(lam_sub))**2))
        res = minimize_scalar(_obj, bounds=(-10, 10), method="bounded")
        C_loo.append(10**res.x)
    print("    C_lambda LOO range: {:.3e} .. {:.3e}".format(
        min(C_loo), max(C_loo)))
    print("    mean +/- sd: {:.3e} +/- {:.3e}  (CoV = {:.1%})".format(
        np.mean(C_loo), np.std(C_loo), np.std(C_loo)/np.mean(C_loo)))
    # Ktory usuniety material zmienia C najbardziej?
    print("\n    Sensitivity per material (|C_LOO - C_full|/C_full):")
    for i, r in enumerate(rows):
        delta = abs(C_loo[i] - C_A0) / C_A0
        print("      {:<5s} {:>6.1%}".format(r["name"], delta))

    # -------- Physical sanity: lambda_TGP ranges? --------
    print("\n" + "-" * 92)
    print("  Sanity: zakres lambda_TGP vs literatura")
    print("-" * 92)
    lam_TGP_best = pred_best
    print("    lambda_lit: [{:.2f}, {:.2f}]".format(lam_arr.min(), lam_arr.max()))
    print("    lambda_TGP: [{:.2f}, {:.2f}]".format(
        lam_TGP_best.min(), lam_TGP_best.max()))
    # ratio
    ratio = lam_TGP_best / lam_arr
    print("    ratio lam_TGP/lam_lit:")
    for r, ra in zip(rows, ratio):
        print("      {:<5s}: {:.2f}".format(r["name"], ra))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r04")
    print("=" * 92)
    print("\n  Fit                             RMS_log  params")
    for name, rms, *_ in fits:
        print("  {:<28s} {:>8.4f}  (best={})".format(
            name, rms, "YES" if (name == best_fit[0]) else ""))
    print("\n  Best RMS_log = {:.4f}".format(best_fit[1]))

    # Target ocena
    target = 0.15
    if best_fit[1] < target:
        conclusion = "SILNY SIGNAL: lambda_ep lezy w TGP Eq.5 strukturze"
    elif best_fit[1] < 0.30:
        conclusion = "UMIARKOWANY: TGP dostarcza czesc wariancji lambda_ep"
    else:
        conclusion = "SLABY: TGP nie zamyka lambda_ep - brakuje czynnikow"
    print("\n  " + conclusion)

    print("""
  Interpretacja:
    - Jesli A0 zamyka (RMS<0.15): TGP struktura DIREKTNIE daje lambda_ep.
      Wtedy unifikacja SC/rho(T) przez McMillan + BG.

    - Jesli A3 (N_F=q=1) zamyka: lambda_ep ~ F_TGP * N_F (fizyczne bo
      Eliashberg lambda = 2 * N_F * <I^2> / (M * omega^2); N_F jest tam jawny).

    - Jesli A1 daje q ~ 0.5-1: TGP pokrywa dynamike (fonon + orbita), ale
      DOS dodaje dodatkowe skalowanie.

    - Jesli A2 alpha != 1: nieliniowosc - sprzezenie non-perturbacyjne.

  Nastepne kroki:
    - r05: pelny fit rho(T_i) na wszystkich 4 T simultaneously uzywajac
           lambda_TGP jako input
    - r06: rozszerzyc dataset (Mo, Ta, W, Be, Ti, Bi, Sb) - zrobi LOO
           bardziej wiarygodnym
    - r07: porownanie lambda_TGP * Theta_D (Eq.5 + Debye) z rho_slope
""")


if __name__ == "__main__":
    main()
