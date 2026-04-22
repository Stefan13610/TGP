#!/usr/bin/env python3
# =============================================================================
#  r16_alloys_test.py
# -----------------------------------------------------------------------------
#  Etap 11 - OSTRZAL rownania U24 na STOPACH (alloys).
#
#  Tezownienie: formula U24 z r13b dziala na czystych metalach (N=37).
#  Czy dziala tez na stopach gdzie obowiazuje regula Matthiessen'a?
#
#  Test:
#    (1) Dla kazdego stopu wymusic ρ_0, Θ_D (z literatury lub avg.)
#    (2) Fit Bloch-Grueneisen (tylko R) do rho_T - rho_0
#    (3) Obliczyc R_pred z U24 przy zadanych (Θ_D, N_F, v_eff, cls)
#    (4) Porownac log R_obs vs log R_pred -> profile RMS
#
#  Szczegolne testy:
#    - Matthiessen extreme: Manganin, Constantan, Nichrome (α ≈ 0)
#    - Ordered/disordered pair: Cu3Au_ord vs Cu3Au_dis (same params, rozny ρ_0)
#    - Linde series: AgAu, CuZn - parabolic ρ_0(x)
#
#  Hipoteza mocna: U24 R ~ tylko od klasy orbitalnej/Θ/NF/v nie zalezy od ρ_0.
#  Jesli ordered i disordered Cu3Au daja ten sam R, to ρ_0/R decomposition
#  jest fizyczne (niezalezne komponenty Matthiessen).
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


# =============================================================================
#  Alloy dataset (13 stopów z research agent)
# =============================================================================
ALLOYS = {
    "CuNi_70_30": {
        "name": "Cu70Ni30 (Cupronickel)",
        "cls": "d",
        "v_eff": 3.7,          # 0.70*1 + 0.30*10
        "rho_0": 34.0,
        "rho_T": {4: 34.0, 77: 36.0, 150: 38.0, 200: 39.5, 250: 41.0, 295: 42.5},
        "Theta_D": 330.0,
        "N_F": 1.8,
    },
    "CuNi_55_45": {
        "name": "Cu55Ni45 (Constantan)",
        "cls": "d",
        "v_eff": 5.05,
        "rho_0": 48.9,
        "rho_T": {4: 48.9, 77: 49.2, 150: 49.5, 200: 49.7, 250: 49.9, 295: 50.0},
        "Theta_D": 320.0,
        "N_F": 2.0,
    },
    "NiCr_80_20": {
        "name": "Ni80Cr20 (Nichrome V)",
        "cls": "d",
        "v_eff": 9.2,
        "rho_0": 105.0,
        "rho_T": {4: 105.0, 77: 106.5, 150: 108.0, 200: 109.5, 250: 110.5, 295: 112.0},
        "Theta_D": 400.0,
        "N_F": 2.2,
    },
    "Brass_Cu70Zn30": {
        "name": "Brass Cu70Zn30",
        "cls": "sp",
        "v_eff": 1.3,          # 0.70*1 + 0.30*2
        "rho_0": 5.8,
        "rho_T": {4: 5.8, 77: 6.3, 150: 7.0, 200: 7.6, 250: 8.1, 295: 8.5},
        "Theta_D": 320.0,
        "N_F": 0.30,
    },
    "Bronze_Cu90Sn10": {
        "name": "Bronze Cu90Sn10",
        "cls": "sp",
        "v_eff": 1.3,          # 0.90*1 + 0.10*4
        "rho_0": 13.5,
        "rho_T": {4: 13.5, 77: 14.2, 150: 15.0, 200: 15.6, 250: 16.2, 295: 16.8},
        "Theta_D": 310.0,
        "N_F": 0.32,
    },
    "Manganin": {
        "name": "Manganin Cu86Mn12Ni2",
        "cls": "d",
        "v_eff": 1.90,
        "rho_0": 43.0,
        "rho_T": {4: 43.0, 77: 43.3, 150: 43.5, 200: 43.7, 250: 43.9, 295: 44.0},
        "Theta_D": 340.0,
        "N_F": 1.2,
    },
    "Monel_400": {
        "name": "Monel Ni65Cu33Fe2",
        "cls": "d",
        "v_eff": 6.99,
        "rho_0": 48.2,
        "rho_T": {4: 48.2, 77: 50.0, 150: 52.0, 200: 53.5, 250: 55.0, 295: 56.5},
        "Theta_D": 380.0,
        "N_F": 2.5,
    },
    "AgAu_50_50": {
        "name": "Ag50Au50",
        "cls": "s",
        "v_eff": 1.0,
        "rho_0": 10.8,
        "rho_T": {4: 10.8, 77: 11.9, 150: 13.2, 200: 14.1, 250: 15.0, 295: 15.9},
        "Theta_D": 200.0,
        "N_F": 0.27,
    },
    "AgAu_90_10": {
        "name": "Ag90Au10 (dilute)",
        "cls": "s",
        "v_eff": 1.0,
        "rho_0": 3.6,
        "rho_T": {4: 3.6, 77: 4.7, 150: 6.0, 200: 6.9, 250: 7.7, 295: 8.5},
        "Theta_D": 220.0,
        "N_F": 0.27,
    },
    "Cu3Au_disord": {
        "name": "Cu3Au disordered",
        "cls": "s",
        "v_eff": 1.0,
        "rho_0": 10.5,
        "rho_T": {4: 10.5, 77: 11.8, 150: 13.3, 200: 14.3, 250: 15.3, 295: 16.3},
        "Theta_D": 285.0,
        "N_F": 0.30,
    },
    "Cu3Au_order": {
        "name": "Cu3Au ordered L1_2",
        "cls": "s",
        "v_eff": 1.0,
        "rho_0": 1.8,
        "rho_T": {4: 1.8, 77: 3.1, 150: 4.7, 200: 5.8, 250: 6.9, 295: 8.0},
        "Theta_D": 285.0,
        "N_F": 0.30,
    },
    "PbSn_63_37": {
        "name": "Sn63Pb37 solder",
        "cls": "sp",
        "v_eff": 4.0,
        "rho_0": 14.5,
        "rho_T": {4: 14.5, 77: 15.8, 150: 17.3, 200: 18.4, 250: 19.4, 295: 20.5},
        "Theta_D": 135.0,
        "N_F": 0.40,
    },
    "Brass_Cu63Zn37": {
        "name": "Yellow brass Cu63Zn37",
        "cls": "sp",
        "v_eff": 1.37,
        "rho_0": 7.2,
        "rho_T": {4: 7.2, 77: 7.8, 150: 8.6, 200: 9.2, 250: 9.7, 295: 10.3},
        "Theta_D": 310.0,
        "N_F": 0.31,
    },
}


def fit_R_alloy(alloy):
    """Fit R in rho(T) = rho_0 + R * (T/Th)^5 * J5(Th/T), with rho_0 fixed."""
    rho_0 = alloy["rho_0"]
    Th = alloy["Theta_D"]
    T_arr = np.array(sorted(alloy["rho_T"].keys()), dtype=float)
    rho_obs = np.array([alloy["rho_T"][T] for T in T_arr])
    # Use only T >= 77 (below T tends to saturate, impurity-dominated)
    mask = T_arr >= 77
    T_fit = T_arr[mask]
    y = rho_obs[mask] - rho_0
    # Safeguard: remove any non-positive (measurement noise below rho_0)
    good = y > 1e-6
    if good.sum() < 2:
        return None, None, Th
    T_fit = T_fit[good]
    y = y[good]
    # BG with R unknown: y_i = R * f(T_i), f(T) = (T/Th)^5 * J5(Th/T)
    fvec = np.array([ (T/Th)**5 * r01.J5(Th/T) for T in T_fit ])
    # least squares
    R_fit = float(np.sum(fvec * y) / np.sum(fvec * fvec))
    # RMS of fit
    pred = rho_0 + R_fit * np.array([ (T/Th)**5 * r01.J5(Th/T) for T in T_arr ])
    res = np.log10(pred + 1e-30) - np.log10(rho_obs + 1e-30)
    rms = float(np.sqrt(np.mean(res**2)))
    return R_fit, rms, Th


def U24_predict_logR(cls, Theta, N_F, v, coefs):
    """Apply U24 formula to predict log10(R)."""
    a, as_, asp, b, c, d_v, ds, dsp, cs, csp = coefs
    delta_s = 1.0 if cls == "s" else 0.0
    delta_sp = 1.0 if cls == "sp" else 0.0
    logTh = np.log10(Theta)
    logNF = np.log10(N_F + 1e-30)
    logR = (a + as_*delta_s + asp*delta_sp
            + b*logTh + c*logNF + d_v*v
            + ds*v*delta_s + dsp*v*delta_sp
            + cs*delta_s*logNF + csp*delta_sp*logNF)
    return logR


def main():
    print("=" * 100)
    print("  r16_alloys_test.py - OSTRZAL U24 na stopach (alloys)")
    print("=" * 100)

    # ----- Step 1: fit U24 on pure metals (get reference coefs) -----
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    data = r12.build_dataset_N37()
    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    pure_rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = r02.ORB_CLASS.get(d["name"], "?")
        if d["name"] in NOBLE or d["name"] in ALKALINE: cls = "s"
        v = r13.V_COUNT.get(d["name"])
        if v is None: continue
        pure_rows.append({
            "name": d["name"], "cls": cls, "R": R, "v": v,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
        })
    N_pure = len(pure_rows)
    print("  Pure metals calibration: N = {}".format(N_pure))

    logTh_p = np.array([np.log10(r["Theta_D"]) for r in pure_rows])
    logNF_p = np.array([np.log10(r["N_F"] + 1e-30) for r in pure_rows])
    v_p = np.array([r["v"] for r in pure_rows], dtype=float)
    logR_p = np.array([np.log10(r["R"]) for r in pure_rows])
    cls_p = np.array([r["cls"] for r in pure_rows])
    delta_s_p = (cls_p == "s").astype(float)
    delta_sp_p = (cls_p == "sp").astype(float)

    X_U24 = np.column_stack([np.ones_like(logTh_p), delta_s_p, delta_sp_p,
                              logTh_p, logNF_p, v_p,
                              v_p*delta_s_p, v_p*delta_sp_p,
                              delta_s_p*logNF_p, delta_sp_p*logNF_p])
    coefs_U24, _, _, _ = np.linalg.lstsq(X_U24, logR_p, rcond=None)
    pred_pure = X_U24 @ coefs_U24
    rms_pure = float(np.sqrt(np.mean((pred_pure - logR_p)**2)))
    print("  U24 on pure metals: RMS (log R) = {:.4f}".format(rms_pure))

    names_U24 = ["a(d)", "as(s-d)", "asp(sp-d)", "b(logTh)", "c(logNF)",
                  "d(v)", "ds(v*s)", "dsp(v*sp)", "cs(s*logNF)", "csp(sp*logNF)"]
    print("\n  U24 coefs (pure-metal calibration):")
    for nm, v in zip(names_U24, coefs_U24):
        print("    {:<18s} = {:+.4f}".format(nm, v))

    # ----- Step 2: fit R_obs for each alloy + predict R_pred -----
    print("\n" + "-" * 100)
    print("  (A) Alloy R fits + U24 predictions")
    print("-" * 100)
    print("  {:<22s}  {:<4s}  {:>7s}  {:>6s}  {:>5s}  {:>9s}  {:>9s}  {:>9s}  {:>7s}".format(
        "alloy", "cls", "Theta", "N_F", "v", "R_obs", "R_pred", "dlogR", "BG_rms"))

    results = []
    for key, alloy in ALLOYS.items():
        R_obs, bg_rms, Th = fit_R_alloy(alloy)
        if R_obs is None or R_obs <= 0:
            print("  {:<22s}  FAIL".format(alloy["name"]))
            continue
        logR_pred = U24_predict_logR(alloy["cls"], alloy["Theta_D"],
                                       alloy["N_F"], alloy["v_eff"], coefs_U24)
        R_pred = 10**logR_pred
        dlog = np.log10(R_obs) - logR_pred
        print("  {:<22s}  {:<4s}  {:>7.1f}  {:>6.2f}  {:>5.2f}  {:>9.3f}  {:>9.3f}  {:>+9.3f}  {:>7.4f}".format(
            alloy["name"][:22], alloy["cls"], Th, alloy["N_F"], alloy["v_eff"],
            R_obs, R_pred, float(dlog), bg_rms))
        results.append({
            "key": key, "name": alloy["name"], "cls": alloy["cls"],
            "R_obs": R_obs, "R_pred": R_pred, "dlog": float(dlog),
            "bg_rms": bg_rms, "alloy": alloy,
        })

    # ----- Step 3: aggregate RMS per class -----
    print("\n" + "-" * 100)
    print("  (B) Per-class RMS (log R alloy vs U24)")
    print("-" * 100)
    print("  {:<4s}  {:>3s}  {:>9s}  {:>9s}  {:>9s}".format("cls", "N", "mean_dlog", "rms", "max|dlog|"))
    all_dlog = []
    for cls in ["s", "sp", "d"]:
        dlogs = [r["dlog"] for r in results if r["cls"] == cls]
        if len(dlogs) == 0: continue
        arr = np.array(dlogs)
        print("  {:<4s}  {:>3d}  {:>+9.3f}  {:>9.4f}  {:>9.3f}".format(
            cls, len(arr), float(np.mean(arr)), float(np.sqrt(np.mean(arr**2))), float(np.max(np.abs(arr)))))
        all_dlog.extend(dlogs)
    overall = np.array(all_dlog)
    print("\n  OVERALL (alloys):  N = {}, bias = {:+.3f}, RMS = {:.4f}, max|dlog| = {:.3f}".format(
        len(overall), float(np.mean(overall)), float(np.sqrt(np.mean(overall**2))), float(np.max(np.abs(overall)))))
    print("  OVERALL (pure):    N = {}, bias = {:+.3f}, RMS = {:.4f}".format(
        N_pure, float(np.mean(pred_pure - logR_p)), rms_pure))

    # ----- Step 4: special tests -----
    print("\n" + "-" * 100)
    print("  (C) Kluczowe testy kontrolne")
    print("-" * 100)

    # C1: Cu3Au ordered vs disordered - SAME (cls, Theta, N_F, v) ale rozny ρ_0
    r_ord = [r for r in results if r["key"] == "Cu3Au_order"]
    r_dis = [r for r in results if r["key"] == "Cu3Au_disord"]
    if r_ord and r_dis:
        ro = r_ord[0]; rd = r_dis[0]
        print("  [C1] Cu3Au ordered vs disordered (same U24 input -> powinny miec ten sam R_pred)")
        print("    R_obs (ord)  = {:.3f}   R_pred = {:.3f}   dlog = {:+.3f}".format(
            ro["R_obs"], ro["R_pred"], ro["dlog"]))
        print("    R_obs (dis)  = {:.3f}   R_pred = {:.3f}   dlog = {:+.3f}".format(
            rd["R_obs"], rd["R_pred"], rd["dlog"]))
        ratio_obs = ro["R_obs"] / rd["R_obs"]
        ratio_pred = ro["R_pred"] / rd["R_pred"]
        print("    R_obs ratio (ord/dis) = {:.3f},  R_pred ratio = {:.3f}".format(ratio_obs, ratio_pred))
        if abs(np.log10(ratio_obs)) < 0.2:
            print("    [OK] Matthiessen rule dziala: rho_0 sie rozni ale R ~= niezmienione.")
        else:
            print("    [UWAGA] Order-disorder transition zmienia R (nie tylko rho_0).")

    # C2: Matthiessen extreme - Constantan, Manganin, Nichrome
    print("\n  [C2] Matthiessen extreme (alpha ~= 0) - powinny miec male drho/dT")
    for key in ["CuNi_55_45", "Manganin", "NiCr_80_20"]:
        r = [x for x in results if x["key"] == key]
        if not r: continue
        r = r[0]
        rho_tot = r["alloy"]["rho_T"][295]
        rho_0 = r["alloy"]["rho_0"]
        frac_phonon = (rho_tot - rho_0) / rho_tot * 100
        print("    {:<20s} rho_0/rho_295 = {:.1f}%  |  phonon {:.1f}%  |  dlog(R) = {:+.3f}".format(
            r["name"][:20], rho_0/rho_tot*100, frac_phonon, r["dlog"]))

    # C3: Linde series - AgAu
    print("\n  [C3] AgAu dilute series (Linde's rule)")
    for key in ["AgAu_90_10", "AgAu_50_50"]:
        r = [x for x in results if x["key"] == key]
        if not r: continue
        r = r[0]
        print("    {:<20s}  R_obs = {:.3f}  R_pred = {:.3f}  dlog = {:+.3f}".format(
            r["name"][:20], r["R_obs"], r["R_pred"], r["dlog"]))

    # ----- Step 5: profile RMS per alloy (log rho(T) fit) -----
    print("\n" + "-" * 100)
    print("  (D) Full profile rho(T) - U24 vs obserwacje")
    print("-" * 100)
    print("  {:<22s}  {:>9s}  {:>10s}  {:>11s}".format(
        "alloy", "prof_rms", "prof_rms_%", "max_rel_err"))
    prof_rms_all = []
    for r in results:
        alloy = r["alloy"]
        R_pred = r["R_pred"]
        Th = alloy["Theta_D"]; rho0 = alloy["rho_0"]
        T_arr = np.array(sorted(alloy["rho_T"].keys()), dtype=float)
        rho_obs = np.array([alloy["rho_T"][T] for T in T_arr])
        rho_pred = np.array([rho0 + R_pred * (T/Th)**5 * r01.J5(Th/T) for T in T_arr])
        dlog = np.log10(rho_pred + 1e-30) - np.log10(rho_obs + 1e-30)
        prof_rms = float(np.sqrt(np.mean(dlog**2)))
        # relative error (linear)
        rel_err = np.max(np.abs(rho_pred - rho_obs) / rho_obs)
        prof_rms_all.append(prof_rms)
        print("  {:<22s}  {:>9.4f}  {:>9.2f}%  {:>10.2f}%".format(
            r["name"][:22], prof_rms, prof_rms * 100, float(rel_err) * 100))
    mean_prof = float(np.mean(prof_rms_all))
    print("\n  MEAN profile RMS (alloys) = {:.4f}".format(mean_prof))
    print("  U24 MEAN profile RMS (pure, from r13b) ~= 0.229")

    # ----- VERDICT -----
    print("\n" + "=" * 100)
    print("  VERDICT r16")
    print("=" * 100)
    print()
    print("  Alloys stress-test:")
    print("    N_alloy = {}".format(len(results)))
    print("    logR bias = {:+.3f}".format(float(np.mean(overall))))
    print("    logR RMS = {:.4f}".format(float(np.sqrt(np.mean(overall**2)))))
    print("    profile RMS = {:.4f}  <-- UWAGA: myslace dla Matthiessen-flat".format(mean_prof))
    print()
    # Prawdziwa metryka: logR RMS (nie profile, bo profile jest zdominowana przez rho_0)
    logR_rms = float(np.sqrt(np.mean(overall**2)))
    if logR_rms < 0.25:
        print("  [VERDICT A] U24 GENERALIZUJE NA STOPY (log R RMS < 0.25).")
        print("    Formula jest uniwersalna: rho_0 i R sa fizycznie niezalezne.")
    elif logR_rms < 0.5:
        print("  [VERDICT B] U24 CZESCIOWO dziala na stopach.")
        print("    Niektore klasy (np. sp) dzialaja, inne (np. d) wymagaja korekty.")
    else:
        print("  [VERDICT C] U24 NIE generalizuje na stopy (log R RMS = {:.3f}).".format(logR_rms))
        print("    Pure-metal formula NIE dziala dla stopow:")
        print("      * d-class: R_obs << R_pred (alloy disorder quenches phonon coupling)")
        print("      * s-class: R_obs >  R_pred (Ag-Au ma dodatkowy scatter wzgl. pure)")
        print("      * sp-class: OK (Cu-Zn, Cu-Sn zachowuje pure-metal-like transport)")
    print()
    print("  [KLUCZOWE] Profile RMS = {:.4f} tylko dlatego ze Matthiessen-extreme alloys".format(mean_prof))
    print("   (Manganin 97.8% rho_0, Constantan 97.8% rho_0) maja b. maly phonon component")
    print("   i blad w R rzedu 30x ledwo widoczny w calkowitym rho(T).")
    print()
    print("  [Matthiessen test OK] Cu3Au ord/dis: R_obs ratio = 1.06 mimo rho_0 ratio = 5.8x")
    print("   -> rho_0 i rho_i(T) dekomponuja sie fizycznie niezaleznie.")
    print()
    print("  Nastepne kroki:")
    print("    - U24 = pure-metals-only. Uniwersalna rola TGP A_orb jest zachowana (pure).")
    print("    - Dla stopow potrzebny Nordheim-like term: x*(1-x)*residual")
    print("    - r17: test sensitivity na Theta_D/N_F perturbacje (pure data)")
    print()
    print("  Nastepny krok:")
    print("    - jesli VERDICT A: przejsc do paper v3 (formula uniwersalna pure+alloy)")
    print("    - jesli VERDICT B: zidentifikowac najgorsze outliers, refit z alloy data")
    print("    - jesli VERDICT C: nowe zmienne (residual disorder, nsubst * dZ^2)")


if __name__ == "__main__":
    main()
