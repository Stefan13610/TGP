#!/usr/bin/env python3
# =============================================================================
#  r12_gparam_model.py
# -----------------------------------------------------------------------------
#  Etap 7 projektu rho(T): weryfikacja per-grupa wzorca z r11 poprzez:
#
#  (1) Rozszerzenie datasetu o brakujace d-metale: Cr (VI), Tc (VII), Lu (III)
#      -> d-class z N=21 do N=24
#      -> grupa VI z N=2 do N=3
#      -> grupa VII z N=1 do N=2
#      -> grupa III z N=3 do N=4 (Lu jak Y i La)
#
#  (2) Continuous g-parameter: Testujemy czy jeden model z grupa jako
#      ciagly parametr (g = numer grupy w PT, 3..10 dla Sc..Ni) odtwarza
#      per-grupa pattern z r11.
#
#  Modele (ANOVA-style):
#     M0:  log R = a + b * log(Theta) + c * log(N_F)                    [3 param]
#     M1:  log R = a + b * log(Theta) + c * log(N_F) + d * g           [4 param]
#     M2:  + e * g * log(Theta)                                         [5 param]
#     M3:  + f * g * log(N_F)                                           [6 param]
#     M4:  + h * g^2                                                    [7 param]
#     M5:  + i * g^2 * log(Theta) + j * g^2 * log(N_F)                  [9 param]
#
#  AIC/BIC porownanie. Jesli M1-M3 istotnie lepsze niz M0 -> continuous g
#  zastepuje per-grupa. Jesli M3-M5 nie poprawiaja istotnie -> liniowe
#  g efekty wystarcza.
#
#  LOO CV na best model: sprawdzenie czy nie overfit.
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

spec11 = importlib.util.spec_from_file_location(
    "r11", os.path.join(os.path.dirname(__file__), "r11_dshell_subclass.py"))
r11 = importlib.util.module_from_spec(spec11); sys.modules["r11"] = r11
spec11.loader.exec_module(r11)


# -----------------------------------------------------------------------------
# Nowe materialy: Cr (VI), Tc (VII), Lu (III)
# Schema: (name, Theta_D, rho_0, r77, r295, r500, r1000, a_latt, N_F, Z)
# -----------------------------------------------------------------------------
RHO_EXT3 = [
    # Cr: bcc 3d5 4s1, group VI. CRC. SDW ponizej 311K powoduje anomalia
    # ale ponad T_N pure BG. rho(295)=12.7.
    # Theta_D zmienne wg zrodla: CRC=630, Desai=460. Uzywamy 630 (from spec heat).
    # N_F ~ 0.55 (DFT).
    ("Cr",  630, 0.050,  1.600, 12.700, 28.00,  80.00,  2.884, 0.55, 24),
    # Tc: hcp 4d5 5s2, group VII. Radioaktywny, ale ISS dane istnieja:
    # Theta_D=454, rho(295)=22.6 (literature), a=2.740
    # N_F ~ 0.57 (DFT McMillan data).
    ("Tc",  454, 0.050,  3.200, 22.600, 44.00,  95.00,  2.740, 0.57, 43),
    # Lu: hcp 5d1 4f14 6s2, group III. CRC Handbook.
    # Theta_D=210, rho(295)=58.2. Lu pelnym 4f^14 ale walencja 5d1 6s2.
    # N_F ~ 0.95 (DFT-Materials Project).
    ("Lu",  210, 1.000, 13.000, 58.200,105.00, 185.00,  3.505, 0.95, 71),
]

COORD_EXT3 = {"Cr": 8, "Tc": 12, "Lu": 12}

ORB_CLASS_EXT3 = {"Cr": "d", "Tc": "d", "Lu": "d"}

LAMBDA_EP_EXT3 = {
    "Cr": 0.44,   # literature (McMillan dla SDW state)
    "Tc": 0.90,   # MRI/Savrasov (Tc=7.8K!)
    "Lu": 0.42,
}

D_GROUP_EXT = {
    "Cr": "VI",
    "Tc": "VII",
    "Lu": "III",
}

D_ROW_EXT = {
    "Cr": "3d",
    "Tc": "4d",
    "Lu": "5d",
}

# Mapowanie grupa -> numer kolumny w PT
GROUP_NUM = {
    "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7,
    "VIII": 8, "IX": 9, "X": 10,
}


def build_dataset_N37():
    full = (list(r00.RHO_DATA) + list(r06.RHO_EXT) +
            list(r10.RHO_EXT2) + list(RHO_EXT3))
    names = [x[0] for x in full]
    assert len(names) == len(set(names)), "duplicate names"
    out = []
    for row in full:
        d = {
            "name":    row[r00.COL_NAME],
            "Theta_D": row[r00.COL_THETA],
            "rho_0":   row[r00.COL_RHO_0],
            "rho":     {T: row[c]
                        for T, c in zip(r00.T_POINTS, r00.RHO_COLS)},
            "a_latt":  row[r00.COL_A_LATT],
            "N_F":     row[r00.COL_NF],
            "Z":       row[r00.COL_Z],
        }
        out.append(d)
    return out


def fit_model(X, y):
    """OLS fit. Return (coefs, pred, rms, aic, bic) where AIC/BIC are
    Gaussian-likelihood estimators."""
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coefs
    res = y - pred
    n = len(y)
    k = X.shape[1]
    rms = np.sqrt(np.mean(res**2))
    # Gaussian AIC: 2k - 2 * ln(L)  ~= n * ln(sum(res^2)/n) + 2k + const
    ss = np.sum(res**2)
    aic = n * np.log(ss / n + 1e-30) + 2 * k
    bic = n * np.log(ss / n + 1e-30) + k * np.log(n)
    return coefs, pred, rms, aic, bic, res


def model_matrix(logTh, logNF, g, order=0):
    """Zbuduj macierz X dla modelu:
        order=0: a + b*logTh + c*logNF                            (3 param)
        order=1: + d*g                                            (4 param)
        order=2: + e*g*logTh                                      (5 param)
        order=3: + f*g*logNF                                      (6 param)
        order=4: + h*g^2                                          (7 param)
        order=5: + i*g^2*logTh + j*g^2*logNF                      (9 param)
    """
    cols = [np.ones_like(logTh), logTh, logNF]
    if order >= 1: cols += [g]
    if order >= 2: cols += [g * logTh]
    if order >= 3: cols += [g * logNF]
    if order >= 4: cols += [g**2]
    if order >= 5: cols += [g**2 * logTh, g**2 * logNF]
    return np.column_stack(cols)


def main():
    print("=" * 92)
    print("  r12_gparam_model.py - d-class z continuous g-parameter (N=24)")
    print("=" * 92)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(ORB_CLASS_EXT3)
    r02.COORD.update(COORD_EXT3)

    data = build_dataset_N37()
    print("\n  N_total = {}".format(len(data)))

    # Rozszerz D_GROUP z r11
    D_GROUP = dict(r11.D_GROUP)
    D_GROUP.update(D_GROUP_EXT)
    D_ROW = dict(r11.D_ROW)
    D_ROW.update(D_ROW_EXT)

    def get_class(name):
        return r02.ORB_CLASS.get(name, "?")

    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    # Build d-class rows
    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = get_class(d["name"])
        if d["name"] in NOBLE or d["name"] in ALKALINE: cls = "s"
        if cls != "d": continue
        g_lab = D_GROUP.get(d["name"])
        if g_lab is None: continue
        rows.append({
            "name": d["name"], "cls": cls, "R": R,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "group": g_lab, "g_num": GROUP_NUM[g_lab],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
            "row": D_ROW.get(d["name"], "?"),
        })

    print("  d-class N = {}".format(len(rows)))
    # Breakdown per group
    groups = {}
    for r in rows:
        groups.setdefault(r["group"], []).append(r)
    print("\n  Per-grupa breakdown:")
    for g, lst in sorted(groups.items(), key=lambda x: GROUP_NUM[x[0]]):
        print("    {:<5s} (g={:d}, N={:d}): {}".format(
            g, GROUP_NUM[g], len(lst), ", ".join(r["name"] for r in lst)))

    # Arrays
    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    g_arr = np.array([r["g_num"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    names = [r["name"] for r in rows]

    # -------- (A) Model comparison (AIC/BIC) --------
    print("\n" + "-" * 92)
    print("  (A) Porownanie modeli continuous-g (N={})".format(len(rows)))
    print("-" * 92)
    print("  {:>5s}  {:>3s}  {:>9s}  {:>9s}  {:>9s}  {:>9s}".format(
        "model", "k", "RMS", "AIC", "BIC", "dAIC"))

    models = []
    for order in range(6):
        X = model_matrix(logTh, logNF, g_arr, order=order)
        coefs, pred, rms, aic, bic, res = fit_model(X, logR)
        models.append({
            "order": order, "k": X.shape[1], "coefs": coefs,
            "pred": pred, "rms": rms, "aic": aic, "bic": bic, "res": res,
        })

    aic_min = min(m["aic"] for m in models)
    bic_min = min(m["bic"] for m in models)
    for m in models:
        tag = " (best AIC)" if m["aic"] == aic_min else ""
        if m["bic"] == bic_min: tag += " (best BIC)"
        print("  M{:<4d}  {:>3d}  {:>9.4f}  {:>9.3f}  {:>9.3f}  {:>+9.3f}{}".format(
            m["order"], m["k"], m["rms"], m["aic"], m["bic"],
            m["aic"] - aic_min, tag))

    # Show coefficients
    print("\n  Coefficients (best AIC model):")
    best = [m for m in models if m["aic"] == aic_min][0]
    coef_names = ["a(const)", "b(logTh)", "c(logNF)", "d(g)",
                  "e(g*logTh)", "f(g*logNF)", "h(g^2)",
                  "i(g^2*logTh)", "j(g^2*logNF)"][:best["k"]]
    for nm, v in zip(coef_names, best["coefs"]):
        print("    {:<14s} = {:+.4f}".format(nm, v))

    # -------- (B) Per-group effective b(g), c(g) from best model --------
    print("\n" + "-" * 92)
    print("  (B) Efektywne b(g), c(g) z best modelu (M{})".format(best["order"]))
    print("-" * 92)
    print("  Wyciagamy przeciwnik fitu zaleznie od grupy:")
    # M3 has: a + b*logTh + c*logNF + d*g + e*g*logTh + f*g*logNF
    # For a specific g, effective b_eff = b + e*g, c_eff = c + f*g
    coefs = best["coefs"]
    k = best["k"]
    print("    {:<5s}  {:>3s}  {:>9s}  {:>9s}".format("grupa", "g", "b_eff", "c_eff"))
    for g_lab in sorted(groups.keys(), key=lambda g: GROUP_NUM[g]):
        g_val = GROUP_NUM[g_lab]
        # Base
        b_eff = coefs[1] if k >= 2 else 0.0
        c_eff = coefs[2] if k >= 3 else 0.0
        # Linear g
        if best["order"] >= 2:
            b_eff += coefs[4] * g_val
        if best["order"] >= 3:
            c_eff += coefs[5] * g_val
        if best["order"] >= 5:
            b_eff += coefs[7] * g_val**2
            c_eff += coefs[8] * g_val**2
        print("    {:<5s}  {:>3d}  {:>+9.3f}  {:>+9.3f}".format(g_lab, g_val, b_eff, c_eff))

    # Porownanie z r11 per-group
    print("\n  Porownanie z r11 per-grupa fittem (tylko grupy N>=3):")
    r11_fits = {
        "III": ("T", 1.00), "IV": ("T", 1.25), "V": ("T", 2.08),
        "X":   ("T", 1.13),
        "VIII": ("N", 0.83), "IX": ("N", 0.57),
    }
    print("    {:<5s}  {:<10s}  {:<14s}  {:<14s}".format(
        "grupa", "r11_form", "r12_b_eff", "r12_c_eff"))
    for g_lab in ["III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]:
        g_val = GROUP_NUM[g_lab]
        b_eff = coefs[1]
        c_eff = coefs[2]
        if best["order"] >= 2: b_eff += coefs[4] * g_val
        if best["order"] >= 3: c_eff += coefs[5] * g_val
        if best["order"] >= 5:
            b_eff += coefs[7] * g_val**2
            c_eff += coefs[8] * g_val**2
        r11_info = r11_fits.get(g_lab, ("--", 0.0))
        print("    {:<5s}  {:<10s}  {:>+9.3f}      {:>+9.3f}".format(
            g_lab, "{}={:+.2f}".format(r11_info[0], r11_info[1]),
            b_eff, c_eff))

    # -------- (C) LOO CV on best model --------
    print("\n" + "-" * 92)
    print("  (C) LOO CV na best modelu (M{})".format(best["order"]))
    print("-" * 92)
    print("  {:<5s}  {:<4s}  {:>9s}  {:>9s}  {:>9s}".format(
        "name", "grp", "R_obs", "R_LOO", "dlog_R"))
    loo_dlog = []
    order_best = best["order"]
    X_full = model_matrix(logTh, logNF, g_arr, order=order_best)
    for i in range(len(rows)):
        mask = np.ones(len(rows), bool); mask[i] = False
        coefs_i, _, _, _ = np.linalg.lstsq(X_full[mask], logR[mask], rcond=None)
        pred_i = X_full[i] @ coefs_i
        dlog = pred_i - logR[i]
        loo_dlog.append(dlog)
        print("  {:<5s}  {:<4s}  {:>9.2f}  {:>9.2f}  {:>+9.3f}".format(
            rows[i]["name"], rows[i]["group"],
            rows[i]["R"], 10**pred_i, dlog))

    loo_rms = np.sqrt(np.mean(np.array(loo_dlog)**2))
    full_rms = best["rms"]
    print("\n  LOO dlog_R RMS = {:.4f}".format(loo_rms))
    print("  Full fit  RMS = {:.4f}".format(full_rms))
    print("  Overfit ratio  = {:.2f}x".format(loo_rms / max(full_rms, 1e-3)))

    # -------- (D) Full BG profile rho(T_i) z best modelu --------
    print("\n" + "-" * 92)
    print("  (D) Pelny profil rho(T_i) z M{} predictions".format(order_best))
    print("-" * 92)
    print("  {:<5s}  {:<4s}  {:>9s}  {:>9s}  {:>9s}".format(
        "name", "grp", "R_obs", "R_pred", "profRMS"))
    prof_rms_all = []
    for i, r in enumerate(rows):
        R_pred = 10**best["pred"][i]
        Th = r["Theta_D"]; rho0 = r["rho_0"]
        rms_pts = []
        for T in r00.T_POINTS:
            rho_p = r01.rho_bg(T, rho0, R_pred, Th)
            rho_o = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_p) - np.log10(rho_o))
        prms = np.sqrt(np.mean(np.array(rms_pts)**2))
        prof_rms_all.append((r["name"], r["group"], prms))
        print("  {:<5s}  {:<4s}  {:>9.2f}  {:>9.2f}  {:>9.4f}".format(
            r["name"], r["group"], r["R"], R_pred, prms))
    vals = [x[2] for x in prof_rms_all]
    print("\n  Profile RMS (M{}):  mean = {:.4f}, median = {:.4f}, max = {:.4f}".format(
        order_best, np.mean(vals), np.median(vals), max(vals)))

    # -------- (E) Compare full profile RMS: M0 vs M_best vs per-grupa r11 --------
    print("\n" + "-" * 92)
    print("  (E) Profile RMS: M0 (universal) vs M{} (g-model) vs r11 per-grupa".format(
        order_best))
    print("-" * 92)
    # M0 predictions
    M0 = models[0]
    prof_m0 = []
    for i, r in enumerate(rows):
        R_pred = 10**M0["pred"][i]
        Th = r["Theta_D"]; rho0 = r["rho_0"]
        rms_pts = []
        for T in r00.T_POINTS:
            rho_p = r01.rho_bg(T, rho0, R_pred, Th)
            rho_o = r["rho_obs"][T]
            rms_pts.append(np.log10(rho_p) - np.log10(rho_o))
        prms = np.sqrt(np.mean(np.array(rms_pts)**2))
        prof_m0.append(prms)

    print("    M0 (universal):    mean profile RMS = {:.4f}".format(np.mean(prof_m0)))
    print("    M{} (best):         mean profile RMS = {:.4f}".format(
        order_best, np.mean(vals)))
    print("    r11 (per-grupa):   mean profile RMS = 0.1053 (N=21)")

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r12")
    print("=" * 92)

    m_best_aic = best["order"]
    dAIC_M0 = models[0]["aic"] - aic_min
    dAIC_M1 = models[1]["aic"] - aic_min
    print("\n  Best model (AIC): M{}  (k={}, RMS={:.4f})".format(
        m_best_aic, best["k"], best["rms"]))
    print("  dAIC M0 (universal): {:+.2f}  {}".format(
        dAIC_M0, "(silne dowody przeciw)" if dAIC_M0 > 10 else "(slabe)"))
    print("  dAIC M1 (+g linear): {:+.2f}".format(dAIC_M1))

    print("\n  Interpretacja:")
    if m_best_aic >= 2 and dAIC_M0 > 10:
        print("    g-zalezne wykladniki sa ISTOTNE -> per-grupa pattern z r11 jest REALNY")
        print("    Continuous g-model opisuje go jednym fitem z {} parametrow.".format(best["k"]))
    elif m_best_aic <= 1 and dAIC_M0 < 10:
        print("    g-zaleznosci wykladnikow NIE istotne statystycznie.")
        print("    M0 (universal) jest OK - r11 per-grupa pattern byl artefaktem N=3.")
    else:
        print("    Mixed: g-intercept istotne, ale g-slope pod presja AIC.")

    # Check if LOO not too much bigger than full fit
    overfit = loo_rms / max(full_rms, 1e-3)
    if overfit < 2.0:
        print("\n  LOO CV OK (overfit ratio {:.2f}x < 2x) -> model GENERALIZUJE.".format(overfit))
    else:
        print("\n  LOO CV WARNING (overfit ratio {:.2f}x >= 2x) -> model overfits.".format(overfit))

    print("\n  Nastepne kroki:")
    print("    - Jesli M0 wystarczy: reintegracja d-class, szukamy lepszego modelu")
    print("    - Jesli M2/M3 best: Paper v3 = continuous g-model z 5-6 parametrami")
    print("    - Jesli M5 best: potrzeba kwadratowych efektow g^2")
    print("    - r13: zbadac s-class (5 mat) i sp-class (8) tak samo solidnie")


if __name__ == "__main__":
    main()
