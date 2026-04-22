#!/usr/bin/env python3
# =============================================================================
#  r08_class_profile.py
# -----------------------------------------------------------------------------
#  r07 pokazal: class-specific D_ThNF fit (R = C_cls * Theta_D^b * N_F^c):
#    sp (N=8): b=+0.25, c=-0.48, RMS=0.161
#    d  (N=12): b=+0.84, c=+0.29, RMS=0.219
#
#  Hipoteza: wykladnik sp c=-0.48 jest NAPRAWDE fizyczny albo driven by
#  Be+Bi (low-N_F outliers).
#
#  r08 testuje:
#    (A) Class-specific profil pelny rho(T_i) - czy R_pred przez BG daje
#        dobre rho(77, 295, 500, 1000)?
#    (B) Sub-class sp: "normal" sp (N_F > 0.1) vs "low-carrier" (N_F<0.05)
#    (C) LOO per class: stabilnosc wykladnikow?
#    (D) Cross-class prediction: fitujemy na sp, predykujemy d? Nie, to bez sensu
#        ale LOO powinien sprawdzic czy fit solidny.
#    (E) Porownanie z literature lambda per class
#
#  Cel: ocenic czy class-specific form jest GOTOWA do paperu v3.
# =============================================================================
import numpy as np
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

spec6 = importlib.util.spec_from_file_location(
    "r06", os.path.join(os.path.dirname(__file__), "r06_extension.py"))
r06 = importlib.util.module_from_spec(spec6); sys.modules["r06"] = r06
spec6.loader.exec_module(r06)


def fit_class(Th, NF, R):
    """Fit R = C * Th^b * NF^c (multilinear log-log)."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th),
                         np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, b, c = coefs
    pred = 10**logC * Th**b * NF**c
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, b, c, rms, pred


def fit_class_theta_only(Th, R):
    """Fit R = C * Th^b (prostszy, 2 parametry)."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = 10**coefs[0] * Th**coefs[1]
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**coefs[0], coefs[1], rms, pred


def main():
    print("=" * 92)
    print("  r08_class_profile.py - class-specific pelny profil rho(T_i)")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)

    data = r06.build_extended_dataset()

    # Zbuduj rows z reklasyfikacja noble "sp" -> "s"
    NOBLE = {"Cu", "Ag", "Au"}
    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        z = r06.get_coord(d["name"]); cls = r06.get_class(d["name"])
        lam = r06.get_lambda(d["name"])
        # Reklasyfikacja noble do 's'
        if d["name"] in NOBLE:
            cls = "s"
        F = r02.tgp_factor_rho(d, z, use_g=False)
        rows.append({
            "name": d["name"], "cls": cls, "R": R,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"], "F_TGP": F,
            "lam_lit": lam, "rho_obs": d["rho"], "rho_0": d["rho_0"],
            "Z": d["Z"],
        })

    print("\n  N_total = {}".format(len(rows)))
    class_counts = {}
    for r in rows:
        class_counts[r["cls"]] = class_counts.get(r["cls"], 0) + 1
    print("  Klasy:", class_counts)

    # -------- (A) Class-specific fit + profile rho(T_i) --------
    print("\n" + "-" * 92)
    print("  (A) Class-specific fit R = C_cls * Theta^b * N_F^c + profil rho(T_i)")
    print("-" * 92)

    # Fit per class
    classes = {}
    for r in rows:
        classes.setdefault(r["cls"], []).append(r)

    fits = {}
    for cls, lst in classes.items():
        Th = np.array([r["Theta_D"] for r in lst])
        NF = np.array([r["N_F"] for r in lst])
        R = np.array([r["R"] for r in lst])
        if len(lst) >= 4:
            C, b, c, rms, pred = fit_class(Th, NF, R)
            fits[cls] = {"C": C, "b": b, "c": c, "rms": rms}
            print("  {:<5s} (N={}): C = {:.3e}, b = {:+.3f}, c = {:+.3f}, RMS_R = {:.4f}".format(
                cls, len(lst), C, b, c, rms))
        elif len(lst) >= 3:
            # tylko Th fit
            C, b, rms, pred = fit_class_theta_only(Th, R)
            fits[cls] = {"C": C, "b": b, "c": 0.0, "rms": rms}
            print("  {:<5s} (N={}): C = {:.3e}, b = {:+.3f} (tylko Th), RMS_R = {:.4f}".format(
                cls, len(lst), C, b, rms))
        else:
            # Fallback: C = median(R/(Th * NF))
            ratios = [r["R"] / (r["Theta_D"] * r["N_F"]) for r in lst]
            C = np.median(ratios)
            fits[cls] = {"C": C, "b": 1.0, "c": 1.0, "rms": 0.0}
            print("  {:<5s} (N={}): za maly sample, C = median ratio = {:.3e}".format(
                cls, len(lst), C))

    # -------- (A.cd) Full profile rho(T_i) test --------
    print("\n  Full profil rho(T_i) z R_class_pred:")
    print("  {:<5s} {:>4s} {:>8s} {:>10s} {:>10s} {:>8s} {:>10s} {:>8s}".format(
        "name", "cls", "R_obs", "R_pred", "rho(77)o", "rho(77)p",
        "rho(1000)", "profRMS"))
    profile_rms_all = []
    for r in rows:
        cls = r["cls"]
        f = fits[cls]
        R_pred = f["C"] * r["Theta_D"]**f["b"] * r["N_F"]**f["c"]
        rho0 = r["rho_0"]; Th = r["Theta_D"]
        rms_pts = []
        for T in r00.T_POINTS:
            rho_pred = r01.rho_bg(T, rho0, R_pred, Th)
            rho_obs = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_pred) - np.log10(rho_obs))
        rms_mat = np.sqrt(np.mean(np.array(rms_pts)**2))
        profile_rms_all.append((r["name"], cls, rms_mat))
        r_o = r["rho_obs"]
        p_77 = r01.rho_bg(77.0, rho0, R_pred, Th)
        p_1000 = r01.rho_bg(1000.0, rho0, R_pred, Th)
        print("  {:<5s} {:>4s} {:>8.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.2f} {:>8.4f}".format(
            r["name"], cls, r["R"], R_pred, r_o[77.0], p_77,
            r_o[1000.0], rms_mat))

    # Stats
    prof_vals = [x[2] for x in profile_rms_all]
    print("\n  Profile RMS per-material statystyki:")
    print("    mean   = {:.4f}".format(np.mean(prof_vals)))
    print("    median = {:.4f}".format(np.median(prof_vals)))
    print("    max    = {:.4f}  ({})".format(
        max(prof_vals), profile_rms_all[np.argmax(prof_vals)][0]))

    print("\n  Histogram per-material profile RMS:")
    bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 5.0)]
    for lo, hi in bins:
        members = [(n, c, v) for n, c, v in profile_rms_all if lo <= v < hi]
        cnt = len(members)
        mem_str = ", ".join("{}({})".format(n, c) for n, c, _ in members)
        print("    [{:.1f}, {:.1f}): {:2d}  {}".format(lo, hi, cnt, mem_str))

    print("\n  Per-class profile RMS:")
    for cls in fits:
        vs = [v for _, c, v in profile_rms_all if c == cls]
        print("    {:<5s} (N={}):  mean = {:.3f}, max = {:.3f}".format(
            cls, len(vs), np.mean(vs), max(vs) if vs else 0))

    # -------- (B) Sub-class sp: normalne vs low-carrier --------
    print("\n" + "-" * 92)
    print("  (B) Sub-class sp: normal (N_F > 0.08) vs low-carrier (N_F < 0.05)")
    print("-" * 92)
    sp_rows = [r for r in rows if r["cls"] == "sp"]
    sp_normal = [r for r in sp_rows if r["N_F"] > 0.08]
    sp_low = [r for r in sp_rows if r["N_F"] < 0.05]
    print("    sp_normal (N={}): {}".format(
        len(sp_normal), [r["name"] for r in sp_normal]))
    print("    sp_low    (N={}): {}".format(
        len(sp_low), [r["name"] for r in sp_low]))

    if len(sp_normal) >= 4:
        Th_n = np.array([r["Theta_D"] for r in sp_normal])
        NF_n = np.array([r["N_F"] for r in sp_normal])
        R_n = np.array([r["R"] for r in sp_normal])
        C_n, b_n, c_n, rms_n, _ = fit_class(Th_n, NF_n, R_n)
        print("    sp_normal: C = {:.2e}, b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
            C_n, b_n, c_n, rms_n))

    if len(sp_low) >= 2:
        Th_l = np.array([r["Theta_D"] for r in sp_low])
        NF_l = np.array([r["N_F"] for r in sp_low])
        R_l = np.array([r["R"] for r in sp_low])
        # Only 2 points - fit constant
        ratio = R_l / (Th_l * NF_l)
        print("    sp_low ratios R/(Th*NF): {}".format(
            ", ".join("{}:{:.2e}".format(r["name"], r["R"]/(r["Theta_D"]*r["N_F"]))
                      for r in sp_low)))

    # -------- (C) LOO per class --------
    print("\n" + "-" * 92)
    print("  (C) LOO per class: stabilnosc b, c")
    print("-" * 92)
    for cls, lst in classes.items():
        if len(lst) < 5: continue
        Th_all = np.array([r["Theta_D"] for r in lst])
        NF_all = np.array([r["N_F"] for r in lst])
        R_all = np.array([r["R"] for r in lst])
        b_loo = []; c_loo = []
        for i in range(len(lst)):
            mask = np.ones(len(lst), bool); mask[i] = False
            _, b_i, c_i, _, _ = fit_class(Th_all[mask], NF_all[mask], R_all[mask])
            b_loo.append(b_i); c_loo.append(c_i)
        print("    {:<5s}: b LOO = {:+.3f} +/- {:.3f}  (range {:+.3f}..{:+.3f})".format(
            cls, np.mean(b_loo), np.std(b_loo), min(b_loo), max(b_loo)))
        print("    {:<5s}: c LOO = {:+.3f} +/- {:.3f}  (range {:+.3f}..{:+.3f})".format(
            cls, np.mean(c_loo), np.std(c_loo), min(c_loo), max(c_loo)))
        # Max leverage
        for i, r in enumerate(lst):
            db = abs(b_loo[i] - fits[cls]["b"])
            dc = abs(c_loo[i] - fits[cls]["c"])
            if db > 0.15 or dc > 0.15:
                print("       {:<5s} leverage: db={:+.3f}, dc={:+.3f}".format(
                    r["name"], b_loo[i] - fits[cls]["b"],
                    c_loo[i] - fits[cls]["c"]))

    # -------- (D) LOO cross-validation: predykcja held-out --------
    print("\n" + "-" * 92)
    print("  (D) LOO CV: fit bez materialu X, przewidz X, profil RMS")
    print("-" * 92)
    print("  {:<5s} {:>4s} {:>10s} {:>10s} {:>10s} {:>10s}".format(
        "name", "cls", "R_obs", "R_LOO", "dlog_R", "prof_RMS"))
    loo_prof = []
    for i, r_test in enumerate(rows):
        cls = r_test["cls"]
        # Fit na innych w tej klasie
        same_cls = [x for x in rows if x["cls"] == cls and x["name"] != r_test["name"]]
        if len(same_cls) < 3:
            print("  {:<5s} {:>4s}  skip (za maly sample)".format(
                r_test["name"], cls))
            continue
        Th_s = np.array([x["Theta_D"] for x in same_cls])
        NF_s = np.array([x["N_F"] for x in same_cls])
        R_s = np.array([x["R"] for x in same_cls])
        if len(same_cls) >= 4:
            C_s, b_s, c_s, _, _ = fit_class(Th_s, NF_s, R_s)
        else:
            C_s, b_s, _, _ = fit_class_theta_only(Th_s, R_s)
            c_s = 0.0
        R_pred = C_s * r_test["Theta_D"]**b_s * r_test["N_F"]**c_s
        dlog = np.log10(R_pred) - np.log10(r_test["R"])
        # Profile RMS z predykcja
        rms_pts = []
        for T in r00.T_POINTS:
            rho_pred = r01.rho_bg(T, r_test["rho_0"], R_pred, r_test["Theta_D"])
            rho_obs = r_test["rho_obs"][T]
            rms_pts.append(np.log10(rho_pred) - np.log10(rho_obs))
        prof_rms = np.sqrt(np.mean(np.array(rms_pts)**2))
        loo_prof.append((r_test["name"], cls, dlog, prof_rms))
        print("  {:<5s} {:>4s} {:>10.3f} {:>10.3f} {:>+10.3f} {:>10.4f}".format(
            r_test["name"], cls, r_test["R"], R_pred, dlog, prof_rms))

    # LOO CV stats
    if loo_prof:
        dlogs = [x[2] for x in loo_prof]
        profs = [x[3] for x in loo_prof]
        print("\n  LOO CV statystyki:")
        print("    mean |dlog R| = {:.3f}".format(np.mean(np.abs(dlogs))))
        print("    RMS  profile = {:.3f}".format(np.sqrt(np.mean(np.array(profs)**2))))
        print("    mean profile = {:.3f}".format(np.mean(profs)))
        print("    max  profile = {:.3f}  ({})".format(
            max(profs), loo_prof[np.argmax(profs)][0]))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r08")
    print("=" * 92)
    print("\n  Class-specific R = C * Theta^b * N_F^c  (na fit full data):")
    for cls, f in sorted(fits.items(), key=lambda x: -classes[x[0]].__len__()):
        print("    {:<5s} (N={}):  C={:.2e}  b={:+.3f}  c={:+.3f}  RMS_R={:.4f}".format(
            cls, class_counts[cls], f["C"], f["b"], f["c"], f["rms"]))

    print("\n  Per-class profile rho(T_i) performance:")
    for cls in fits:
        vs = [v for _, c, v in profile_rms_all if c == cls]
        good = sum(1 for v in vs if v < 0.25)
        bad = sum(1 for v in vs if v > 0.50)
        print("    {:<5s} (N={}): mean={:.3f}, {}/{} <0.25 RMS, {} outliers (>0.5)".format(
            cls, len(vs), np.mean(vs), good, len(vs), bad))

    # Verdict wygląd
    mean_prof = np.mean(prof_vals)
    if mean_prof < 0.25:
        status = "GOTOWE do paperu v3 (mean profile RMS < 0.25)"
    elif mean_prof < 0.35:
        status = "BLISKO gotowego - potrzebne oczyszczenie 2-3 outlierow"
    else:
        status = "WYMAGANA dalsza praca - za duze RMS"
    print("\n  " + status)

    print("""
  Uwagi:
    - Jesli sp_normal vs sp_low rozni sie znacznie (c flipuje znak po
      wykluczeniu Be/Bi): c=-0.48 to artefakt low-carrier outlierow
      a nie prawdziwy dla 'zwyklych' sp metali.

    - Jesli LOO CV profile RMS ~ full fit RMS: forma class-specific jest
      robust, nie overfit.

    - Następnym krokiem (jesli r08 zamyka): r09 weryfikacja physical
      meaning wykladnika d c=+0.29, oraz paper draft.
""")


if __name__ == "__main__":
    main()
