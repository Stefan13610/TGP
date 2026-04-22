#!/usr/bin/env python3
# =============================================================================
#  r07_outlier_class.py
# -----------------------------------------------------------------------------
#  r06 pokazal: D4 niestabilny (dryf 60-120%), RMS 0.20 -> 0.38 na N=23.
#  Dominujace problemy:
#    (i)   Bi - semimetal z N_F=0.006 (100x nizsze niz reszta) - ekstremalny
#          leverage point (LOO |d|=2.07)
#    (ii)  Cu/Ag/Au - pure s-transport (d-shell inertny), a my dajemy im
#          klasa "sp" z A_orb=0.207 (tez sp)
#    (iii) Ti/Be - ekstremalne Theta_D (420/1440 K) + niesp. TM
#
#  r07 pyta:
#    (A) Czy reklasyfikacja Cu/Ag/Au "sp" -> "s" (A_orb=-0.111) stabilizuje fit?
#    (B) Czy wykluczenie Bi (semimetal) daje uniwersalny BG?
#    (C) Czy class-specific D4 (noble / sp / d) ma stabilne wykladniki?
#    (D) Czy jest jakis uniwersalny sub-form ktory generalizuje dobrze na N=22?
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

spec6 = importlib.util.spec_from_file_location(
    "r06", os.path.join(os.path.dirname(__file__), "r06_extension.py"))
r06 = importlib.util.module_from_spec(spec6); sys.modules["r06"] = r06
spec6.loader.exec_module(r06)


def fit_D4(F, Th, NF, R):
    """Fit R = C * F^a * Th^b * NF^c w log-log multilinear."""
    X = np.column_stack([np.ones_like(F), np.log10(F + 1e-30),
                         np.log10(Th), np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, a, b, c = coefs
    pred = 10**logC * F**a * Th**b * NF**c
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, a, b, c, rms, pred


def fit_D3(F, Th, R):
    """Fit R = C * F^a * Th^b."""
    X = np.column_stack([np.ones_like(F), np.log10(F + 1e-30), np.log10(Th)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, a, b = coefs
    pred = 10**logC * F**a * Th**b
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, a, b, rms, pred


def fit_D_ThNF(Th, NF, R):
    """Fit R = C * Th^b * NF^c."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th), np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, b, c = coefs
    pred = 10**logC * Th**b * NF**c
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, b, c, rms, pred


def fit_D_Theta_only(Th, R):
    """Fit R = C * Th^b."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, b = coefs
    pred = 10**logC * Th**b
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, b, rms, pred


def main():
    print("=" * 92)
    print("  r07_outlier_class.py - analiza outlierow + class-specific D4")
    print("=" * 92)

    # Hot-patch r02
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)

    data = r06.build_extended_dataset()

    # Wyodrebnic wszystkie rows z R_fit
    rows = []
    for d in data:
        R, rms, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        z = r06.get_coord(d["name"])
        cls = r06.get_class(d["name"])
        lam = r06.get_lambda(d["name"])
        F = r02.tgp_factor_rho(d, z, use_g=False)
        rows.append({
            "name": d["name"], "cls": cls, "R": R,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"], "F_TGP": F,
            "lam_lit": lam, "rho_obs": d["rho"], "rho_0": d["rho_0"],
            "Z": d["Z"], "a_latt": d["a_latt"],
        })

    print("\n  N_total = {}".format(len(rows)))

    # -------- (A) Reklasyfikacja Cu/Ag/Au "sp" -> "s" --------
    print("\n" + "-" * 92)
    print("  (A) Reklasyfikacja Cu/Ag/Au: 'sp' -> 's' (A_orb -> -0.111)")
    print("-" * 92)
    NOBLE = {"Cu", "Ag", "Au"}
    rows_A = []
    for r in rows:
        new_r = dict(r)
        if r["name"] in NOBLE:
            # Skala F_TGP przez (A_s / A_sp)^2
            factor = (r02.A_ORB["s"] / r02.A_ORB["sp"])**2
            new_r["F_TGP"] = r["F_TGP"] * factor
            new_r["cls"] = "s"
        rows_A.append(new_r)

    F_A = np.array([r["F_TGP"] for r in rows_A])
    Th_arr = np.array([r["Theta_D"] for r in rows])
    NF_arr = np.array([r["N_F"] for r in rows])
    R_arr = np.array([r["R"] for r in rows])

    C_A, a_A, b_A, c_A, rms_A, pred_A = fit_D4(F_A, Th_arr, NF_arr, R_arr)
    print("    D4 z reklasyfikacja noble='s':")
    print("    C = {:.3e}, a = {:+.3f}, b = {:+.3f}, c = {:+.3f}".format(
        C_A, a_A, b_A, c_A))
    print("    RMS_log = {:.4f}  (baseline N=23 D4: 0.3841)".format(rms_A))

    # -------- (B) Wykluczenie Bi (semimetal) --------
    print("\n" + "-" * 92)
    print("  (B) Wykluczenie Bi (semimetal, N_F=0.006 = ekstremalny leverage)")
    print("-" * 92)
    mask_B = np.array([r["name"] != "Bi" for r in rows])
    F_B = np.array([r["F_TGP"] for r in rows])[mask_B]
    Th_B = Th_arr[mask_B]; NF_B = NF_arr[mask_B]; R_B = R_arr[mask_B]
    C_B, a_B, b_B, c_B, rms_B, _ = fit_D4(F_B, Th_B, NF_B, R_B)
    print("    D4 bez Bi (N={}):".format(mask_B.sum()))
    print("    C = {:.3e}, a = {:+.3f}, b = {:+.3f}, c = {:+.3f}".format(
        C_B, a_B, b_B, c_B))
    print("    RMS_log = {:.4f}".format(rms_B))

    # -------- (C) Wykluczenie Bi + reklasyfikacja noble --------
    print("\n" + "-" * 92)
    print("  (C) Bi excluded + noble reklasyfikacja")
    print("-" * 92)
    F_C = F_A[mask_B]
    C_C, a_C, b_C, c_C, rms_C, _ = fit_D4(F_C, Th_B, NF_B, R_B)
    print("    D4:  C = {:.3e}, a = {:+.3f}, b = {:+.3f}, c = {:+.3f}".format(
        C_C, a_C, b_C, c_C))
    print("    RMS_log = {:.4f}".format(rms_C))

    # -------- (D) Wykluczenie tez Be (ekstremalny Theta_D=1440) --------
    print("\n" + "-" * 92)
    print("  (D) Bi + Be excluded + noble reklasyfikacja")
    print("-" * 92)
    mask_D = np.array([r["name"] not in ("Bi", "Be") for r in rows])
    F_D = F_A[mask_D]; Th_D = Th_arr[mask_D]
    NF_D = NF_arr[mask_D]; R_D = R_arr[mask_D]
    C_D, a_D, b_D, c_D, rms_D, _ = fit_D4(F_D, Th_D, NF_D, R_D)
    print("    D4 (N={}):  a = {:+.3f}, b = {:+.3f}, c = {:+.3f}".format(
        mask_D.sum(), a_D, b_D, c_D))
    print("    RMS_log = {:.4f}".format(rms_D))

    # -------- (E) Prostsze formy: R = C * Theta^b * N_F^c --------
    print("\n" + "-" * 92)
    print("  (E) Prostsze formy (bez F_TGP): R = C * Theta^b * N_F^c")
    print("-" * 92)
    # Na N=23 (full)
    C_e, b_e, c_e, rms_e, _ = fit_D_ThNF(Th_arr, NF_arr, R_arr)
    print("    N=23 full:  b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
        b_e, c_e, rms_e))
    # Bez Bi
    C_e2, b_e2, c_e2, rms_e2, _ = fit_D_ThNF(Th_B, NF_B, R_B)
    print("    N=22 no-Bi: b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
        b_e2, c_e2, rms_e2))
    # Bez Bi+Be
    C_e3, b_e3, c_e3, rms_e3, _ = fit_D_ThNF(Th_D, NF_D, R_D)
    print("    N=21 no-Bi,Be: b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
        b_e3, c_e3, rms_e3))

    # -------- (F) Class-specific fit: osobny fit dla noble / sp / d --------
    print("\n" + "-" * 92)
    print("  (F) Class-specific: osobny D_ThNF dla each class (N=23)")
    print("-" * 92)
    classes = {}
    for r in rows_A:  # z reklasyfikacja noble
        classes.setdefault(r["cls"], []).append(r)
    for cls, lst in classes.items():
        Th_c = np.array([x["Theta_D"] for x in lst])
        NF_c = np.array([x["N_F"] for x in lst])
        R_c = np.array([x["R"] for x in lst])
        if len(lst) >= 3:
            Cc, bc, cc, rmsc, _ = fit_D_ThNF(Th_c, NF_c, R_c)
            print("    {:<5s} (N={}):  C = {:.2e}, b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
                cls, len(lst), Cc, bc, cc, rmsc))
            print("         members: {}".format([x["name"] for x in lst]))
        else:
            print("    {:<5s} (N={}):  za maly sample".format(cls, len(lst)))

    # -------- (G) R = C * Theta_D alone (najprostszy) --------
    print("\n" + "-" * 92)
    print("  (G) Najprostsza: R = C * Theta_D^b na rozszerzonej")
    print("-" * 92)
    for mask_name, mask_arr in [("N=23 full", np.ones(len(rows), bool)),
                                 ("N=22 no Bi", mask_B),
                                 ("N=21 no Bi,Be", mask_D)]:
        Cg, bg, rmsg, _ = fit_D_Theta_only(Th_arr[mask_arr], R_arr[mask_arr])
        r_, p_ = pearsonr(np.log10(Th_arr[mask_arr]),
                          np.log10(R_arr[mask_arr]))
        print("    {:<14s}:  b = {:+.3f},  C = {:.2e},  RMS = {:.4f},  r = {:.3f}, p = {:.4f}".format(
            mask_name, bg, Cg, rmsg, r_, p_))

    # -------- (H) R = C * lambda_lit * Theta_D^b (poza-out) --------
    print("\n" + "-" * 92)
    print("  (H) Z literaturowym lambda_ep: R = C * lam^a * Theta^b")
    print("-" * 92)
    lam_arr = np.array([r["lam_lit"] if r["lam_lit"] else np.nan
                        for r in rows])
    for mask_name, mask_arr in [("N=23 full", ~np.isnan(lam_arr)),
                                 ("N=22 no Bi", mask_B & ~np.isnan(lam_arr)),
                                 ("N=21 no Bi,Be",
                                  mask_D & ~np.isnan(lam_arr))]:
        la = lam_arr[mask_arr]; Tha = Th_arr[mask_arr]; Ra = R_arr[mask_arr]
        X = np.column_stack([np.ones_like(la), np.log10(la), np.log10(Tha)])
        y = np.log10(Ra)
        coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        predh = 10**coefs[0] * la**coefs[1] * Tha**coefs[2]
        rmsh = np.sqrt(np.mean((np.log10(predh) - y)**2))
        print("    {:<14s}:  a = {:+.3f}, b = {:+.3f}, C = {:.3e}, RMS = {:.4f}".format(
            mask_name, coefs[1], coefs[2], 10**coefs[0], rmsh))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r07")
    print("=" * 92)
    print("""
  Testowe scenariusze (D4 forma):
    Baseline N=23 full:            a=+0.42, b=+0.46, c=-0.00, RMS=0.38
    (A) noble reclass (s):         a={:+.3f}, b={:+.3f}, c={:+.3f}, RMS={:.4f}
    (B) bez Bi:                    a={:+.3f}, b={:+.3f}, c={:+.3f}, RMS={:.4f}
    (C) bez Bi + noble s:          a={:+.3f}, b={:+.3f}, c={:+.3f}, RMS={:.4f}
    (D) bez Bi,Be + noble s:       a={:+.3f}, b={:+.3f}, c={:+.3f}, RMS={:.4f}

  Interpretacja:
    - Jesli (D) reprodukuje N=15 D4 (a=-1.89, b=+1.21, c=+1.05): wyniki r05
      byly faktyczne, ale Bi/Be byly poza BG domeny.
    - Jesli dalej daleko: forma D4 nie uniwersalna - potrzebujemy
      class-specific.

  Najprostsze formy bez F_TGP (E, G):
    - Jesli R = C*Th^b*NF^c zamyka na ~0.30 RMS: prostszy model wygrywa.
    - Jesli niet: F_TGP jest potrzebny (po oczyszczeniu outlierow).
""".format(a_A, b_A, c_A, rms_A,
           a_B, b_B, c_B, rms_B,
           a_C, b_C, c_C, rms_C,
           a_D, b_D, c_D, rms_D))


if __name__ == "__main__":
    main()
