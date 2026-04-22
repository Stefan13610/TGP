#!/usr/bin/env python3
# =============================================================================
#  r03_class_split.py
# -----------------------------------------------------------------------------
#  r02 pokazal: bez g: r=0.57, RMS=0.34. Residuum ma strukture klasowa.
#
#  Tu testujemy:
#    (A) Class-specific C_R: R = C_R^(cls) * F_TGP. Czy CoV per class < 30%?
#    (B) Non-linear fit: R = C * F_TGP^alpha. Czy alpha =/= 1? (1 = Eq.5-like,
#        0.5 = mean-field, 2 = strong-coupling).
#    (C) N_F korekta dla transport: R = C * F_TGP * N_F^q. Fit q.
#        (W przeciwienstwie do P7.12 SC, moze q dla transportu jest male +/-.)
#    (D) Add lambda_tr = lambda_ep * (1 - <cos theta>) - jak to zmienia wynik?
#    (E) Lambda_ep czyste z SC: R = C * lambda_ep^p. Fit na metalach gdzie znamy lambda.
#
#  Cel: znalezc najprostszy universal factor.
# =============================================================================
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar, minimize
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


# Literature lambda_ep (McMillan-Allen) dla metali gdzie znamy
# Zrodlo: Allen & Mitrovic 1982, Grimvall 1981, DFT reviews
LAMBDA_EP_LIT = {
    "Cu": 0.13, "Ag": 0.12, "Au": 0.15,  # slabe s-p sprzezenie
    "Al": 0.43, "Pb": 1.55, "Sn": 0.72,
    "Nb": 1.12, "V":  0.82,
    "Fe": 0.25, "Ni": 0.18,  # niskie lambda_ep, bo FM ma inne scattery
    "Pt": 0.66, "Pd": 0.35,
    "Cd": 0.40, "Zn": 0.38, "Mg": 0.35,
}


def main():
    print("=" * 92)
    print("  r03_class_split.py - klasowa struktura residuow + alternatywne")
    print("  parametryzacje")
    print("=" * 92)

    data = r00.load_dataset()
    R_fit = {}
    for d in data:
        R, _, _, _ = r01.fit_R_only(d)
        if R is not None: R_fit[d["name"]] = R

    # Zbuduj tabele
    rows = []
    for d in data:
        if d["name"] not in R_fit: continue
        z = r02.COORD[d["name"]]
        cls = r02.ORB_CLASS[d["name"]]
        F_TGP = r02.tgp_factor_rho(d, z, use_g=False)
        lam = LAMBDA_EP_LIT.get(d["name"], 0.0)
        rows.append({
            "name": d["name"], "R": R_fit[d["name"]],
            "F_TGP": F_TGP, "cls": cls,
            "N_F": d["N_F"], "Theta_D": d["Theta_D"],
            "lambda_ep": lam,
        })

    # -------- (A) Class-specific C_R --------
    print("\n" + "-" * 92)
    print("  (A) C_R per klasa  (R = C_R^(cls) * F_TGP)")
    print("-" * 92)
    classes = {}
    for r in rows:
        classes.setdefault(r["cls"], []).append(r)

    class_CR = {}
    print("  {:<6s} {:>3s} {:>13s} {:>13s} {:>8s}".format(
        "cls", "N", "median C_R", "mean+/-sd", "CoV"))
    for cls, lst in classes.items():
        CRs = [r["R"] / r["F_TGP"] for r in lst]
        med = np.median(CRs); mu = np.mean(CRs); sd = np.std(CRs)
        class_CR[cls] = med
        print("  {:<6s} {:>3d} {:>13.2e} {:>13.2e} {:>7.1%}".format(
            cls, len(lst), med, mu, sd/mu if mu > 0 else 0))

    # RMS z class-specific
    resid = []
    for r in rows:
        pred = class_CR[r["cls"]] * r["F_TGP"]
        resid.append(np.log10(pred) - np.log10(r["R"]))
    rms_cs = np.sqrt(np.mean(np.asarray(resid)**2))
    print("\n  Per-class C_R RMS_log = {:.4f}".format(rms_cs))
    print("  (vs H0 global C_R = {:.4f})".format(0.3391))

    # -------- (B) Non-linear R = C * F_TGP^alpha --------
    print("\n" + "-" * 92)
    print("  (B) Non-linear: R = C * F_TGP^alpha")
    print("-" * 92)
    F_arr = np.array([r["F_TGP"] for r in rows])
    R_arr = np.array([r["R"] for r in rows])
    logF = np.log10(F_arr + 1e-20)
    logR = np.log10(R_arr)
    # linear fit in log: logR = logC + alpha * logF
    slope, intercept = np.polyfit(logF, logR, 1)
    C_fit = 10**intercept
    print("    Fit: alpha = {:.3f}, logC = {:.3f} -> C = {:.3e}".format(
        slope, intercept, C_fit))
    pred = C_fit * F_arr**slope
    rms_nl = np.sqrt(np.mean((np.log10(pred) - logR)**2))
    print("    RMS_log = {:.4f}".format(rms_nl))

    # -------- (C) R = C * F_TGP * N_F^q --------
    print("\n" + "-" * 92)
    print("  (C) R = C * F_TGP * N_F^q  (dodajemy jawny N_F exponent)")
    print("-" * 92)
    NF_arr = np.array([r["N_F"] for r in rows])

    def obj_Cq(params):
        logC, q = params
        C = 10**logC
        pred = C * F_arr * NF_arr**q
        return np.sqrt(np.mean((np.log10(pred + 1e-20) - logR)**2))

    best = (1e9, None)
    for lc0 in [2, 3, 4, 5]:
        for q0 in [-1, 0, 1, 2]:
            res = minimize(obj_Cq, [lc0, q0], method="Nelder-Mead",
                           options={"xatol": 1e-4, "fatol": 1e-6})
            if res.fun < best[0]:
                best = (res.fun, res.x)
    logC_opt, q_opt = best[1]
    print("    Fit: q = {:.3f}, logC = {:.3f} -> C = {:.3e}".format(
        q_opt, logC_opt, 10**logC_opt))
    print("    RMS_log = {:.4f}".format(best[0]))

    # -------- (D) lambda_ep fit --------
    print("\n" + "-" * 92)
    print("  (D) R = C * lambda_ep^p  (literatura lambda_ep)")
    print("-" * 92)
    lam_arr = np.array([r["lambda_ep"] for r in rows])
    ok = lam_arr > 0
    R_ok = R_arr[ok]
    lam_ok = lam_arr[ok]
    slope2, intercept2 = np.polyfit(np.log10(lam_ok), np.log10(R_ok), 1)
    C_lam = 10**intercept2
    print("    Fit: p = {:.3f}, C = {:.3e}".format(slope2, C_lam))
    pred_lam = C_lam * lam_ok**slope2
    rms_lam = np.sqrt(np.mean((np.log10(pred_lam) - np.log10(R_ok))**2))
    print("    RMS_log = {:.4f}  (N = {})".format(rms_lam, np.sum(ok)))
    r_lam, p_lam = pearsonr(np.log10(lam_ok), np.log10(R_ok))
    print("    Pearson r(log R, log lam_ep) = {:.3f}, p = {:.4f}".format(
        r_lam, p_lam))

    # -------- (E) Lambda_ep + Theta_D skala (fizyczne rho_ph limit) --------
    print("\n" + "-" * 92)
    print("  (E) R = C * lambda_ep * Theta_D  (high-T BG scaling)")
    print("-" * 92)
    Th_arr = np.array([r["Theta_D"] for r in rows])
    Th_ok = Th_arr[ok]
    prod_ok = lam_ok * Th_ok
    slope3, intercept3 = np.polyfit(np.log10(prod_ok), np.log10(R_ok), 1)
    C_lt = 10**intercept3
    print("    Fit: p = {:.3f}, C = {:.3e}".format(slope3, C_lt))
    pred_lt = C_lt * prod_ok**slope3
    rms_lt = np.sqrt(np.mean((np.log10(pred_lt) - np.log10(R_ok))**2))
    print("    RMS_log = {:.4f}".format(rms_lt))
    r_lt, p_lt = pearsonr(np.log10(prod_ok), np.log10(R_ok))
    print("    Pearson r(log R, log(lam*Theta)) = {:.3f}, p = {:.4f}".format(
        r_lt, p_lt))

    # -------- (F) The physical picture: R = C * lambda_ep * Theta_D / N_F? --------
    print("\n" + "-" * 92)
    print("  (F) Fizyczna hipoteza: R = C * lambda_ep / N_F  (Ziman Eq.9.59)")
    print("-" * 92)
    NF_ok = NF_arr[ok]
    ratio_ok = lam_ok / NF_ok
    slope4, intercept4 = np.polyfit(np.log10(ratio_ok), np.log10(R_ok), 1)
    C_r = 10**intercept4
    pred_r = C_r * ratio_ok**slope4
    rms_r = np.sqrt(np.mean((np.log10(pred_r) - np.log10(R_ok))**2))
    r_r, p_r = pearsonr(np.log10(ratio_ok), np.log10(R_ok))
    print("    Fit: p = {:.3f}, C = {:.3e}".format(slope4, C_r))
    print("    RMS_log = {:.4f}".format(rms_r))
    print("    Pearson r(log R, log(lam/NF)) = {:.3f}, p = {:.4f}".format(
        r_r, p_r))

    # Verdict
    print("\n" + "=" * 92)
    print("  VERDICT r03")
    print("=" * 92)
    fits = [
        ("H0 global C_R", 0.3391, 2),
        ("(A) class C_R", rms_cs, 2 + len(class_CR) - 1),
        ("(B) C * F^alpha", rms_nl, 2),
        ("(C) C * F * N_F^q", best[0], 2),
        ("(D) C * lambda_ep^p", rms_lam, 2),
        ("(E) C * (lambda*Theta)^p", rms_lt, 2),
        ("(F) C * (lambda/NF)^p", rms_r, 2),
    ]
    print("\n  {:<28s} {:>10s} {:>8s}".format("Fit", "RMS_log", "params"))
    for name, rms, p in fits:
        print("  {:<28s} {:>10.4f} {:>8d}".format(name, rms, p))

    best_fit = min(fits, key=lambda x: x[1])
    print("\n  Best: {} (RMS = {:.4f})".format(best_fit[0], best_fit[1]))

    print("""
    Interpretacja:
      - H0 = TGP Eq.5 structure (bez g). 0.34 log-RMS jest moderate.
      - (D) lambda_ep alone: test czy R jest przede wszystkim funkcja lambda.
      - (F) Ziman: R ~ lambda/N_F (fizyczne BG).

    Jesli (F) wygrywa:
      TGP_transport = (powerful Ziman formula) + N_F normalization
      Dla paper v3: "TGP daje lambda_ep, a ta daje R(T) przez uniwersalny BG".

    Jesli (A) wygrywa o klasa-specyficznym C_R:
      Orbital-class separacja jest fizyczna; potrzebujemy wiecej materialow.

    Nastepny krok:
      - Rozszerzyc dataset (Mo, Ta, W, Be, Ti) + niskie T (20K, 50K)
      - Sprawdzic czy model r03 najlepszy predykuje rho(T) na kazdym T, nie
        tylko amplitude R
""")


if __name__ == "__main__":
    main()
