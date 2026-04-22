#!/usr/bin/env python3
# =============================================================================
#  r11_dshell_subclass.py
# -----------------------------------------------------------------------------
#  Etap 6 projektu rho(T): diagnoza collapse b_d z r10.
#
#  r08 (N=12 d): b=+0.84
#  r10 (N=21 d): b=+0.02
#
#  Hipoteza: d-class nie jest jednorodna -- wypelnienie podpasma d (group III..X)
#  wplywa na R(Theta_D, N_F). Rozbijmy d-class na podgrupy wedlug grupy w ukladzie
#  okresowym, fitujmy osobno i porownajmy wykladniki.
#
#  Podzial d-materialow (N=21) wedlug grupy:
#    III  (d1+f mixing): Sc, Y, La                 N=3
#    IV   (d2):          Ti, Zr, Hf                N=3
#    V    (d3):          V, Nb, Ta                 N=3
#    VI   (d4/d5):       Mo, W                     N=2
#    VII  (d5/d6):       Re                        N=1
#    VIII (d6/d7):       Fe, Ru, Os                N=3
#    IX   (d7/d8):       Co, Rh, Ir                N=3
#    X    (d8/d9/d10):   Ni, Pd, Pt                N=3
#
#  Dla N>=3 grup testujemy: R = C * Theta^b * N_F^c  (lub Theta-only jezeli N=3).
#
#  Cel: Czy w podgrupach b jest STABILNY i rozni od 0?
#      Jesli TAK: "class-specific" rozdziela na d3..d10; Paper v3 wymaga 8 podklas.
#      Jesli NIE: b naprawde ~= 0 dla d-class, trzeba inszej formy.
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


# Mapa materialow na grupe ukladu okresowego
D_GROUP = {
    # Grupa III (d1 + s2 or f-mixing)
    "Sc":  "III",
    "Y":   "III",
    "La":  "III",
    # Grupa IV (d2 s2)
    "Ti":  "IV",
    "Zr":  "IV",
    "Hf":  "IV",
    # Grupa V (d3 s2)
    "V":   "V",
    "Nb":  "V",
    "Ta":  "V",
    # Grupa VI (d4-d5 s2-s1)
    "Mo":  "VI",
    "W":   "VI",
    # Grupa VII (d5-d6 s2)
    "Re":  "VII",
    # Grupa VIII (d6-d7)
    "Fe":  "VIII",
    "Ru":  "VIII",
    "Os":  "VIII",
    # Grupa IX (d7-d8)
    "Co":  "IX",
    "Rh":  "IX",
    "Ir":  "IX",
    # Grupa X (d8-d9-d10)
    "Ni":  "X",
    "Pd":  "X",
    "Pt":  "X",
}

# Z-ladder: okres w tabeli (3d=4, 4d=5, 5d=6)
D_ROW = {
    # 3d (okres 4)
    "Sc": "3d", "Ti": "3d", "V": "3d", "Fe": "3d", "Co": "3d", "Ni": "3d",
    # 4d (okres 5)
    "Y": "4d", "Zr": "4d", "Nb": "4d", "Mo": "4d", "Ru": "4d", "Rh": "4d", "Pd": "4d",
    # 5d (okres 6)
    "La": "5d", "Hf": "5d", "Ta": "5d", "W": "5d", "Re": "5d", "Os": "5d",
    "Ir": "5d", "Pt": "5d",
}


def fit_class(Th, NF, R):
    X = np.column_stack([np.ones_like(Th), np.log10(Th), np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, b, c = coefs
    pred = 10**logC * Th**b * NF**c
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, b, c, rms


def fit_class_theta_only(Th, R):
    X = np.column_stack([np.ones_like(Th), np.log10(Th)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = 10**coefs[0] * Th**coefs[1]
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**coefs[0], coefs[1], rms


def fit_class_nf_only(NF, R):
    X = np.column_stack([np.ones_like(NF), np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = 10**coefs[0] * NF**coefs[1]
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**coefs[0], coefs[1], rms


def fit_class_const(R):
    """R = C_const (constant)."""
    C = 10**np.mean(np.log10(R))
    rms = np.sqrt(np.mean((np.log10(R) - np.log10(C))**2))
    return C, rms


def main():
    print("=" * 92)
    print("  r11_dshell_subclass.py - diagnoza collapse b_d przez grupowanie d-class")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)

    data = r10.build_full_dataset()
    print("\n  N_total = {}, test na d-class only".format(len(data)))

    def get_class(name):
        if name in r02.ORB_CLASS: return r02.ORB_CLASS[name]
        return "?"

    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    # Build rows
    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = get_class(d["name"])
        if d["name"] in NOBLE or d["name"] in ALKALINE: cls = "s"
        rows.append({
            "name": d["name"], "cls": cls, "R": R,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
        })

    d_rows = [r for r in rows if r["cls"] == "d"]
    print("\n  d-class N = {}".format(len(d_rows)))

    # -------- (A) Grupy --------
    print("\n" + "-" * 92)
    print("  (A) Podgrupy d-class wg grupy ukladu okresowego")
    print("-" * 92)
    groups = {}
    for r in d_rows:
        g = D_GROUP.get(r["name"], "?")
        r["group"] = g
        groups.setdefault(g, []).append(r)

    print("  {:<6s} {:>3s}   {:<30s}".format("grupa", "N", "czlonkowie"))
    for g in sorted(groups.keys()):
        lst = groups[g]
        print("  {:<6s} {:>3d}   {}".format(
            g, len(lst), ", ".join(r["name"] for r in lst)))

    # -------- (B) Fit per grupa --------
    print("\n" + "-" * 92)
    print("  (B) Fit R = C * Theta^b * N_F^c per grupa (jezeli N>=4)")
    print("  Fallback: Theta-only (N=3) lub const (N<3)")
    print("-" * 92)
    print("  {:<6s} {:>3s}   {:>12s}  {:>9s}  {:>9s}  {:>9s}  {:>6s}".format(
        "grupa", "N", "C", "b", "c", "RMS_R", "forma"))
    group_fits = {}
    for g in sorted(groups.keys()):
        lst = groups[g]
        Th = np.array([r["Theta_D"] for r in lst])
        NF = np.array([r["N_F"] for r in lst])
        R = np.array([r["R"] for r in lst])
        if len(lst) >= 4:
            C, b, c, rms = fit_class(Th, NF, R)
            form = "TN"
        elif len(lst) == 3:
            # Try both Theta-only and NF-only, pick better
            C_t, b_t, rms_t = fit_class_theta_only(Th, R)
            C_n, b_n, rms_n = fit_class_nf_only(NF, R)
            if rms_t < rms_n:
                C, b, c, rms = C_t, b_t, 0.0, rms_t; form = "T"
            else:
                C, b, c, rms = C_n, 0.0, b_n, rms_n; form = "N"
        elif len(lst) == 2:
            # Median ratio
            C, rms = fit_class_const(R)
            b = c = 0.0; form = "const"
        else:
            C = R[0]; b = c = 0.0; rms = 0.0; form = "single"
        group_fits[g] = {"C": C, "b": b, "c": c, "rms": rms, "form": form, "N": len(lst)}
        print("  {:<6s} {:>3d}   {:>12.3e}  {:>+9.3f}  {:>+9.3f}  {:>9.4f}  {:>6s}".format(
            g, len(lst), C, b, c, rms, form))

    # -------- (C) Fit per okres (3d vs 4d vs 5d) --------
    print("\n" + "-" * 92)
    print("  (C) Fit per okres (3d vs 4d vs 5d)")
    print("-" * 92)
    rows_3d = [r for r in d_rows if D_ROW.get(r["name"]) == "3d"]
    rows_4d = [r for r in d_rows if D_ROW.get(r["name"]) == "4d"]
    rows_5d = [r for r in d_rows if D_ROW.get(r["name"]) == "5d"]
    print("    3d (N={}): {}".format(len(rows_3d), [r["name"] for r in rows_3d]))
    print("    4d (N={}): {}".format(len(rows_4d), [r["name"] for r in rows_4d]))
    print("    5d (N={}): {}".format(len(rows_5d), [r["name"] for r in rows_5d]))
    print()
    row_fits = {}
    for lab, lst in [("3d", rows_3d), ("4d", rows_4d), ("5d", rows_5d)]:
        Th = np.array([r["Theta_D"] for r in lst])
        NF = np.array([r["N_F"] for r in lst])
        R = np.array([r["R"] for r in lst])
        if len(lst) >= 4:
            C, b, c, rms = fit_class(Th, NF, R)
            row_fits[lab] = (C, b, c, rms, len(lst))
            print("    {:<3s} (N={}): C = {:.3e}, b = {:+.3f}, c = {:+.3f}, RMS = {:.4f}".format(
                lab, len(lst), C, b, c, rms))

    # -------- (D) LOO per grupa --------
    print("\n" + "-" * 92)
    print("  (D) LOO per grupa (N>=4): stabilnosc wykladnikow")
    print("-" * 92)
    for g, lst in sorted(groups.items()):
        if len(lst) < 4: continue
        Th_all = np.array([r["Theta_D"] for r in lst])
        NF_all = np.array([r["N_F"] for r in lst])
        R_all = np.array([r["R"] for r in lst])
        b_loo = []; c_loo = []
        for i in range(len(lst)):
            mask = np.ones(len(lst), bool); mask[i] = False
            _, b_i, c_i, _ = fit_class(Th_all[mask], NF_all[mask], R_all[mask])
            b_loo.append(b_i); c_loo.append(c_i)
        print("    {:<4s} (N={:2d}):  b = {:+.3f} +/- {:.3f}    c = {:+.3f} +/- {:.3f}".format(
            g, len(lst), np.mean(b_loo), np.std(b_loo),
            np.mean(c_loo), np.std(c_loo)))

    # -------- (E) R_subclass predict vs R_d_full predict --------
    print("\n" + "-" * 92)
    print("  (E) Porownanie: R_subclass_pred vs R_d_full_pred (N=21 d)")
    print("-" * 92)
    # r10 d-class (N=21): C=7.122e+01, b=+0.021, c=+0.247
    C_full, b_full, c_full = 7.122e+01, 0.021, 0.247
    print("\n  {:<5s} {:<4s}  {:>9s}  {:>9s}  {:>7s}  {:>9s}  {:>7s}  {:>8s}".format(
        "name", "grp", "R_obs", "R_full", "dlog_f", "R_sub", "dlog_s", "better"))
    sub_wins = 0; full_wins = 0; ties = 0
    sub_res = []; full_res = []
    for r in d_rows:
        g = r["group"]
        f = group_fits.get(g)
        R_obs = r["R"]
        R_full_pred = C_full * r["Theta_D"]**b_full * r["N_F"]**c_full
        if f and f["form"] in ("TN", "T", "N", "const"):
            if f["form"] == "TN":
                R_sub_pred = f["C"] * r["Theta_D"]**f["b"] * r["N_F"]**f["c"]
            elif f["form"] == "T":
                R_sub_pred = f["C"] * r["Theta_D"]**f["b"]
            elif f["form"] == "N":
                R_sub_pred = f["C"] * r["N_F"]**f["c"]
            else:
                R_sub_pred = f["C"]
        else:
            R_sub_pred = R_obs  # single point
        dlog_f = np.log10(R_full_pred) - np.log10(R_obs)
        dlog_s = np.log10(R_sub_pred) - np.log10(R_obs)
        sub_res.append(dlog_s); full_res.append(dlog_f)
        if abs(dlog_s) < abs(dlog_f) - 0.01:
            better = "sub"; sub_wins += 1
        elif abs(dlog_f) < abs(dlog_s) - 0.01:
            better = "full"; full_wins += 1
        else:
            better = "tie"; ties += 1
        print("  {:<5s} {:<4s}  {:>9.2f}  {:>9.2f}  {:>+7.3f}  {:>9.2f}  {:>+7.3f}  {:>8s}".format(
            r["name"], g, R_obs, R_full_pred, dlog_f, R_sub_pred, dlog_s, better))

    rms_full = np.sqrt(np.mean(np.array(full_res)**2))
    rms_sub = np.sqrt(np.mean(np.array(sub_res)**2))
    print("\n  RMS(d-full N=21 fit) = {:.4f}".format(rms_full))
    print("  RMS(sub-grupy fit)  = {:.4f}".format(rms_sub))
    print("  Win counts:  sub={}  full={}  tie={}".format(sub_wins, full_wins, ties))

    # -------- (F) Full BG profile rho(T_i) z sub-class pred --------
    print("\n" + "-" * 92)
    print("  (F) BG profile rho(T_i) na N=21 d-class z sub-class R_pred")
    print("-" * 92)
    print("  {:<5s} {:<4s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}".format(
        "name", "grp", "R_obs", "R_sub", "profRMS", "dR_log"))
    prof_rms_all = []
    for r in d_rows:
        g = r["group"]; f = group_fits.get(g)
        if not f: continue
        if f["form"] == "TN":
            R_pred = f["C"] * r["Theta_D"]**f["b"] * r["N_F"]**f["c"]
        elif f["form"] == "T":
            R_pred = f["C"] * r["Theta_D"]**f["b"]
        elif f["form"] == "N":
            R_pred = f["C"] * r["N_F"]**f["c"]
        elif f["form"] == "const":
            R_pred = f["C"]
        else:
            R_pred = r["R"]
        Th = r["Theta_D"]; rho0 = r["rho_0"]
        rms_pts = []
        for T in r00.T_POINTS:
            rho_p = r01.rho_bg(T, rho0, R_pred, Th)
            rho_o = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_p) - np.log10(rho_o))
        prms = np.sqrt(np.mean(np.array(rms_pts)**2))
        dr = np.log10(R_pred) - np.log10(r["R"])
        prof_rms_all.append((r["name"], g, prms))
        print("  {:<5s} {:<4s}  {:>8.2f}  {:>8.2f}  {:>8.4f}  {:>+8.3f}".format(
            r["name"], g, r["R"], R_pred, prms, dr))

    vals = [x[2] for x in prof_rms_all]
    print("\n  Sub-class d-profile RMS: mean={:.4f}, median={:.4f}, max={:.4f}  ({})".format(
        np.mean(vals), np.median(vals), max(vals),
        prof_rms_all[np.argmax(vals)][0]))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r11")
    print("=" * 92)
    print("\n  Stabilnosc b per grupa (z LOO jesli N>=4):")
    consistent = True
    for g, lst in sorted(groups.items()):
        f = group_fits.get(g)
        if f is None or f["form"] not in ("TN",):
            continue
        if len(lst) >= 4:
            # Rebuild LOO
            Th_all = np.array([r["Theta_D"] for r in lst])
            NF_all = np.array([r["N_F"] for r in lst])
            R_all = np.array([r["R"] for r in lst])
            b_loo = []
            for i in range(len(lst)):
                mask = np.ones(len(lst), bool); mask[i] = False
                _, b_i, _, _ = fit_class(Th_all[mask], NF_all[mask], R_all[mask])
                b_loo.append(b_i)
            sigma_b = np.std(b_loo)
            if sigma_b > abs(f["b"]): consistent = False

    if rms_sub < 0.8 * rms_full:
        print("\n  SUB-GRUPY dzialaja (RMS_sub << RMS_full). Paper v3 = 'per-group' theory.")
        print("    RMS_full = {:.4f}, RMS_sub = {:.4f} ({:.0f}% redukcja)".format(
            rms_full, rms_sub, 100*(rms_full-rms_sub)/rms_full))
    elif rms_sub < rms_full:
        print("\n  Sub-grupy tylko marginalnie lepsze. Nie warto rozbijac na 8 podklas.")
    else:
        print("\n  Sub-grupy NIE poprawiaja: single d-class form (nawet z b~0) jest OK.")

    print("""
  Wnioski metodologiczne:
    1. Jesli sub-grupy daja b mocno rozne (III:+0.8, VI:-0.3, X:+0.1), to znaczy ze
       w d-class jest schowana grupa-zalezna fizyka (d-shell filling).
    2. Jesli wszystkie grupy maja b ~ 0 (zgodne z r10), to d-class naprawde nie ma
       zaleznosci R od Theta_D -- to moze byc prawdziwe odkrycie.
    3. Okres (3d vs 4d vs 5d) moze byc lepszym podzialem niz grupa -- z fizyki SOC
       i spin-orbit coupling wzmacniajacym sie od 3d do 5d.

  Nastepny krok:
    - r12: dodac 4-ci parametr g (d-shell filling, 0..10) do pelnego fitu
      R = C * Theta^b * N_F^c * g^d  -- jesli d rozne od 0, to g istotne
    - r13: sprawdzic czy Z (liczba atomowa) albo okres moga zastapic grupe
""")


if __name__ == "__main__":
    main()
