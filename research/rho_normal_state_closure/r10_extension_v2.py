#!/usr/bin/env python3
# =============================================================================
#  r10_extension_v2.py
# -----------------------------------------------------------------------------
#  Etap 5 projektu rho(T): rozszerzenie datasetu z N=23 do N=34 celem
#  weryfikacji class-specific wykladnikow z r08:
#    d-class (N=12):  R = 4.71e-01 * Theta^+0.84 * N_F^+0.29   (RMS=0.22)
#    sp-class (N=8): R = 3.88e+00 * Theta^+0.25 * N_F^-0.48   (RMS=0.16)
#    s-class (N=3):   R = 6.93e-02 * Theta^+0.84              (RMS=0.03)
#
#  Nowe materialy (N=11):
#  ---- d-class (9 sztuk) ----
#    Ir  (fcc  5d)  - group IX, like Rh
#    Re  (hcp  5d)  - group VII
#    Ru  (hcp  4d)  - group VIII
#    Os  (hcp  5d)  - group VIII
#    Zr  (hcp  4d)  - group IV
#    Hf  (hcp  5d)  - group IV
#    Sc  (hcp  3d)  - group III
#    Y   (hcp  4d)  - group III
#    La  (dhcp 5d)  - group III (low Theta_D=142K)
#  ---- alkaline-earth (2 sztuki - kandydaci do s-class) ----
#    Ca  (fcc  s2)
#    Sr  (fcc  s2)
#
#  Alkalie (K, Rb, Cs, Na, Li) WYKLUCZONE - topia sie nizej 500K,
#  nie mamy rho(500)/rho(1000) w fazie stalej.
#  Ba WYKLUCZONE - m.p.=1000K, rho(1000K) nieodgadnelne.
#
#  Cel: zweryfikowac czy class-specific wykladniki stabilizuja sie na N=34.
#  Jesli d-class b=+0.84, c=+0.29 trzyma sie przy doubled sample, to forma
#  jest poprawna i mozemy ruszac z paperem v3.
#
#  Zrodla:
#    CRC Handbook 97th ed. (2016), Table 12-41 "Resistivity of Pure Metals"
#    Desai et al. 1984 JPCRD 13 (refractory TMs)
#    White & Minges 1997 NIST recommended values
#    Vafek/Mazin DFT (Materials Project) for N(E_F)
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


# -----------------------------------------------------------------------------
# Nowe materialy (N=11)  schema: (name, Th, rho0, r77, r295, r500, r1000, a, N_F, Z)
# -----------------------------------------------------------------------------
RHO_EXT2 = [
    # ---- d-class 4d/5d refractory & early TM ----
    # Ir: fcc 5d group IX. CRC/Desai. rho(295)=5.07. Like Rh ale wyzszy Z.
    ("Ir",  420, 0.010,  0.260,  5.070, 10.50,  25.00,  3.839, 1.09, 77),
    # Re: hcp 5d group VII. Matula/Desai. rho(295)=17.22 Matula.
    ("Re",  416, 0.020,  2.520, 17.220, 33.00,  74.00,  2.761, 0.45, 75),
    # Ru: hcp 4d group VIII. CRC.
    ("Ru",  600, 0.020,  0.420,  7.100, 15.50,  35.00,  2.706, 0.92, 44),
    # Os: hcp 5d group VIII. CRC.
    ("Os",  500, 0.040,  0.810,  8.100, 15.50,  36.00,  2.734, 0.85, 76),
    # Zr: hcp 4d group IV (alpha phase). CRC. m.p.=2128K.
    # rho(1000K) still in alpha (alpha->beta at 1135K).
    ("Zr",  250, 0.500, 10.000, 42.100, 77.00, 140.00,  3.232, 1.30, 40),
    # Hf: hcp 5d group IV (alpha). CRC.
    ("Hf",  213, 0.700,  7.700, 33.000, 63.00, 120.00,  3.195, 0.90, 72),
    # Sc: hcp 3d group III. CRC.
    ("Sc",  360, 1.000, 14.000, 56.200, 98.00, 180.00,  3.309, 1.35, 21),
    # Y: hcp 4d group III. CRC. m.p.=1799K.
    ("Y",   248, 1.000, 14.200, 59.600,100.00, 180.00,  3.647, 1.45, 39),
    # La: dhcp 5d group III. CRC. m.p.=1193K; przechodzi alpha->beta 613K,
    # beta->gamma 1141K; rho(1000) w fazie beta, juz blisko melt.
    ("La",  142, 2.000, 18.000, 61.500,105.00, 180.00,  3.770, 2.50, 57),

    # ---- alkaline-earth (s2) - test s-class extension ----
    # Ca: fcc s2. CRC. m.p.=1115K.
    ("Ca",  229, 0.300,  0.910,  3.360,  6.20,  13.50,  5.588, 0.30, 20),
    # Sr: fcc s2. CRC. m.p.=1050K. rho(1000) blisko melt.
    ("Sr",  147, 1.000,  3.500, 13.200, 23.00,  50.00,  6.084, 0.30, 38),
]

# Coord for new materials
COORD_EXT2 = {
    "Ir": 12, "Re": 12, "Ru": 12, "Os": 12,          # fcc(Ir) + hcp(Re,Ru,Os) = 12
    "Zr": 12, "Hf": 12, "Sc": 12, "Y": 12, "La": 12,  # hcp
    "Ca": 12, "Sr": 12,                               # fcc alkaline-earth
}

# ORB_CLASS for new (default) - nobles/alkaline-earth reclassifiable later
ORB_CLASS_EXT2 = {
    "Ir": "d", "Re": "d", "Ru": "d", "Os": "d",
    "Zr": "d", "Hf": "d", "Sc": "d", "Y": "d", "La": "d",
    "Ca": "sp",  # pierwotnie sp (s2 valence ale mozliwe p mixing)
    "Sr": "sp",  # ^^^
}

# Lambda_ep dla nowych - z Allen/Mitrovic 1982 + Savrasov 1996 + McMillan
LAMBDA_EP_EXT2 = {
    "Ir": 0.34,   # McMillan (niski Tc=0.11K -> lambda ~ 0.34)
    "Re": 0.46,   # Allen-Mitrovic 1982
    "Ru": 0.38,
    "Os": 0.44,
    "Zr": 0.41,
    "Hf": 0.34,
    "Sc": 0.40,
    "Y":  0.38,
    "La": 0.80,   # f-level hybridization boost, McMillan Tc=6K
    "Ca": 0.32,
    "Sr": 0.40,
}


def build_full_dataset():
    """r00 + r06.RHO_EXT + RHO_EXT2 -> lista dict N=34."""
    full = list(r00.RHO_DATA) + list(r06.RHO_EXT) + list(RHO_EXT2)
    names = [x[0] for x in full]
    assert len(names) == len(set(names)), "Duplicate nazwy"
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


def fit_class(Th, NF, R):
    """Fit R = C * Th^b * NF^c (multilinear log-log). Return (C, b, c, rms)."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th),
                         np.log10(NF + 1e-30)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    logC, b, c = coefs
    pred = 10**logC * Th**b * NF**c
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**logC, b, c, rms


def fit_class_theta_only(Th, R):
    """Fit R = C * Th^b (2-param)."""
    X = np.column_stack([np.ones_like(Th), np.log10(Th)])
    y = np.log10(R)
    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = 10**coefs[0] * Th**coefs[1]
    rms = np.sqrt(np.mean((np.log10(pred) - y)**2))
    return 10**coefs[0], coefs[1], rms


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 92)
    print("  r10_extension_v2.py - extension do N=34, weryfikacja class fits")
    print("=" * 92)

    # Hot-patch z obu rozszerzen
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(ORB_CLASS_EXT2)
    r02.COORD.update(COORD_EXT2)

    data = build_full_dataset()
    print("\n  N_total = {} (r00:{} + r06:{} + r10:{})".format(
        len(data), len(r00.RHO_DATA), len(r06.RHO_EXT), len(RHO_EXT2)))

    # -------- Print nowych materialow --------
    print("\n  Nowe w r10 (N={}):".format(len(RHO_EXT2)))
    print("  {:<5s} {:>4s} {:>7s} {:>6s} {:>7s} {:>8s} {:>7s} {:>4s} {:>4s}".format(
        "name", "cls", "Th_D", "a[A]", "N_F", "rho295", "lam_ep", "z", "Z"))
    for row in RHO_EXT2:
        nm = row[0]
        cls = ORB_CLASS_EXT2[nm]
        z = COORD_EXT2[nm]
        lam = LAMBDA_EP_EXT2[nm]
        print("  {:<5s} {:>4s} {:>7.0f} {:>6.3f} {:>7.3f} {:>8.3f} {:>7.3f} {:>4d} {:>4d}".format(
            nm, cls, row[1], row[7], row[8], row[4], lam, z, row[9]))

    # -------- Fit R per material --------
    def get_class(name):
        if name in r02.ORB_CLASS: return r02.ORB_CLASS[name]
        return "?"

    NOBLE = {"Cu", "Ag", "Au"}
    ALKALINE = {"Ca", "Sr"}

    rows = []
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        cls = get_class(d["name"])
        # Reclass: noble and alkaline-earth do 's'
        if d["name"] in NOBLE or d["name"] in ALKALINE:
            cls = "s"
        rows.append({
            "name": d["name"], "cls": cls, "R": R,
            "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_obs": d["rho"], "rho_0": d["rho_0"],
            "Z": d["Z"],
        })

    class_counts = {}
    for r in rows:
        class_counts[r["cls"]] = class_counts.get(r["cls"], 0) + 1
    print("\n  Klasy po reklasyfikacji (noble+alkaline -> s):")
    for cls, n in sorted(class_counts.items(), key=lambda x: -x[1]):
        mem = [r["name"] for r in rows if r["cls"] == cls]
        print("    {:<3s} (N={:2d}): {}".format(cls, n, ", ".join(mem)))

    # -------- (A) Class fit N=34 vs N=23 r08 --------
    print("\n" + "-" * 92)
    print("  (A) Class-specific fit na N=34 vs porownanie z r08 (N=23)")
    print("-" * 92)
    # r08 baseline values
    R08 = {
        "d":  {"C": 4.711e-01, "b": +0.837, "c": +0.293, "rms": 0.2191, "N": 12},
        "sp": {"C": 3.878e+00, "b": +0.254, "c": -0.479, "rms": 0.1612, "N": 8},
        "s":  {"C": 6.929e-02, "b": +0.839, "c":  0.000, "rms": 0.0324, "N": 3},
    }

    classes = {}
    for r in rows:
        classes.setdefault(r["cls"], []).append(r)

    fits = {}
    print("  {:<5s}  {:>4s}  {:>12s}  {:>9s}  {:>9s}  {:>9s}   {:>6s} {:>6s} {:>6s}".format(
        "cls", "N", "C", "b", "c", "RMS_R", "db", "dc", "dRMS"))
    for cls in ["d", "sp", "s"]:
        lst = classes.get(cls, [])
        if len(lst) < 3: continue
        Th = np.array([r["Theta_D"] for r in lst])
        NF = np.array([r["N_F"] for r in lst])
        R = np.array([r["R"] for r in lst])
        if len(lst) >= 4:
            C, b, c, rms = fit_class(Th, NF, R)
        else:
            C, b, rms = fit_class_theta_only(Th, R); c = 0.0
        fits[cls] = {"C": C, "b": b, "c": c, "rms": rms, "N": len(lst)}
        # Porownaj z r08
        o = R08.get(cls, None)
        if o:
            db = b - o["b"]; dc = c - o["c"]; drms = rms - o["rms"]
            print("  {:<5s}  {:>4d}  {:>12.3e}  {:>+9.3f}  {:>+9.3f}  {:>9.4f}   {:>+6.3f} {:>+6.3f} {:>+6.3f}".format(
                cls, len(lst), C, b, c, rms, db, dc, drms))
        else:
            print("  {:<5s}  {:>4d}  {:>12.3e}  {:>+9.3f}  {:>+9.3f}  {:>9.4f}".format(
                cls, len(lst), C, b, c, rms))
    print("\n  (baseline r08:")
    for cls, o in R08.items():
        print("    {:<3s} (N={:2d}):  C={:.3e}  b={:+.3f}  c={:+.3f}  RMS={:.4f}".format(
            cls, o["N"], o["C"], o["b"], o["c"], o["rms"]))

    # -------- (B) Full profile rho(T_i) test na N=34 --------
    print("\n" + "-" * 92)
    print("  (B) Pelny profil rho(T_i) z class-specific R_pred")
    print("-" * 92)
    print("  {:<5s} {:>4s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>8s}".format(
        "name", "cls", "R_obs", "R_pred", "rho(77)", "pred(77)",
        "rho(1000)", "profRMS"))
    profile_rms_all = []
    for r in rows:
        cls = r["cls"]
        f = fits.get(cls, None)
        if f is None: continue
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
        mark = "  " if r["name"] in [x[0] for x in RHO_EXT2] else "  "
        if r["name"] in [x[0] for x in RHO_EXT2]: mark = "* "
        print("{}{:<5s} {:>4s} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.2f} {:>8.4f}".format(
            mark, r["name"], cls, r["R"], R_pred, r_o[77.0], p_77,
            r_o[1000.0], rms_mat))

    # Stats
    print("\n  Profile RMS per-material statystyki (N={}):".format(len(profile_rms_all)))
    prof_vals = np.array([x[2] for x in profile_rms_all])
    print("    mean   = {:.4f}".format(np.mean(prof_vals)))
    print("    median = {:.4f}".format(np.median(prof_vals)))
    print("    max    = {:.4f}  ({})".format(
        np.max(prof_vals), profile_rms_all[np.argmax(prof_vals)][0]))

    # Histogram
    print("\n  Histogram profile RMS:")
    bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 5.0)]
    for lo, hi in bins:
        members = [(n, c, v) for n, c, v in profile_rms_all if lo <= v < hi]
        cnt = len(members)
        mem_str = ", ".join("{}({})".format(n, c) for n, c, _ in members)
        # skrot jesli duzo
        if len(mem_str) > 80: mem_str = mem_str[:77] + "..."
        print("    [{:.1f}, {:.1f}): {:2d}  {}".format(lo, hi, cnt, mem_str))

    # Per-class profile RMS
    print("\n  Per-class profile RMS:")
    for cls in ["d", "sp", "s"]:
        vs = [v for _, c, v in profile_rms_all if c == cls]
        if not vs: continue
        print("    {:<3s} (N={:2d}):  mean = {:.3f}, median = {:.3f}, max = {:.3f}".format(
            cls, len(vs), np.mean(vs), np.median(vs), max(vs)))

    # -------- (C) LOO per class --------
    print("\n" + "-" * 92)
    print("  (C) LOO per class: stabilnosc b, c na N=34")
    print("-" * 92)
    for cls in ["d", "sp", "s"]:
        lst = classes.get(cls, [])
        if len(lst) < 5: continue
        Th_all = np.array([r["Theta_D"] for r in lst])
        NF_all = np.array([r["N_F"] for r in lst])
        R_all = np.array([r["R"] for r in lst])
        b_loo = []; c_loo = []
        max_leverage = []
        for i in range(len(lst)):
            mask = np.ones(len(lst), bool); mask[i] = False
            if len(lst) - 1 >= 4:
                _, b_i, c_i, _ = fit_class(Th_all[mask], NF_all[mask], R_all[mask])
            else:
                _, b_i, _ = fit_class_theta_only(Th_all[mask], R_all[mask])
                c_i = 0.0
            b_loo.append(b_i); c_loo.append(c_i)
            db = abs(b_i - fits[cls]["b"])
            dc = abs(c_i - fits[cls]["c"])
            if db > 0.1 or dc > 0.1:
                max_leverage.append((lst[i]["name"], b_i - fits[cls]["b"], c_i - fits[cls]["c"]))
        print("    {:<3s} (N={:2d}):  b = {:+.3f} +/- {:.3f}  (range {:+.3f}..{:+.3f})".format(
            cls, len(lst), np.mean(b_loo), np.std(b_loo), min(b_loo), max(b_loo)))
        print("    {:<3s} (N={:2d}):  c = {:+.3f} +/- {:.3f}  (range {:+.3f}..{:+.3f})".format(
            cls, len(lst), np.mean(c_loo), np.std(c_loo), min(c_loo), max(c_loo)))
        # Top 3 leverage
        max_leverage.sort(key=lambda x: -abs(x[1])-abs(x[2]))
        for nm, db, dc in max_leverage[:3]:
            print("       {:<5s} db={:+.3f}  dc={:+.3f}".format(nm, db, dc))

    # -------- (D) Sub-class sp: nowe Ca/Sr jako s?  --------
    # Zobacz czy Ca/Sr pasuja do s-class (re-klasyfikacja
    # noble-like) lub do sp (z mieszaniem p)
    print("\n" + "-" * 92)
    print("  (D) Ca, Sr: kandydat s vs sp-normal")
    print("-" * 92)
    # Predict z s-class fit (C, b); predict z sp-class fit (C, b, c)
    for nm in ["Ca", "Sr"]:
        r = next(x for x in rows if x["name"] == nm)
        Th = r["Theta_D"]; NF = r["N_F"]; R_obs = r["R"]
        # Fit s = 3 base noble (Cu Ag Au) + extra sec - ale tu juz zreclass'owalismy
        # do s wiec fits[s] zawiera je. Zeby wyciagnac "czysty" s-fit bez Ca/Sr:
        base_s = [x for x in classes.get("s", []) if x["name"] in NOBLE]
        Th_s = np.array([x["Theta_D"] for x in base_s])
        R_s = np.array([x["R"] for x in base_s])
        if len(base_s) >= 3:
            C_s, b_s, _ = fit_class_theta_only(Th_s, R_s)
            R_s_pred = C_s * Th**b_s
        else:
            R_s_pred = np.nan
        # sp fit bez Be/Bi (sp_normal) - uzywamy z r08
        # r08 sp_normal (N=6): C=5.46e+01, b=-0.132, c=-0.107
        C_sp, b_sp, c_sp = 5.46e+01, -0.132, -0.107
        R_sp_pred = C_sp * Th**b_sp * NF**c_sp
        res_s = np.log10(R_s_pred) - np.log10(R_obs)
        res_sp = np.log10(R_sp_pred) - np.log10(R_obs)
        print("  {:<4s}: R_obs={:.2f}, R_s_pred={:.2f} (dlog={:+.3f}), "
              "R_sp_pred={:.2f} (dlog={:+.3f})  -> {}".format(
            nm, R_obs, R_s_pred, res_s, R_sp_pred, res_sp,
            "s lepszy" if abs(res_s) < abs(res_sp) else "sp lepszy"))

    # -------- Verdict --------
    print("\n" + "=" * 92)
    print("  VERDICT r10")
    print("=" * 92)

    # Class stability
    print("\n  Class-specific R = C * Theta^b * N_F^c (N=34 vs N=23):")
    print("  {:<3s}  {:>4s}  {:>12s}  {:>9s}  {:>9s}  {:>9s}".format(
        "cls", "N", "C", "b", "c", "RMS_R"))
    for cls in ["d", "sp", "s"]:
        f = fits.get(cls)
        if not f: continue
        o = R08.get(cls)
        print("  {:<3s}  {:>4d}  {:>12.3e}  {:>+9.3f}  {:>+9.3f}  {:>9.4f}   (r08: b={:+.3f}, c={:+.3f}, rms={:.4f})".format(
            cls, f["N"], f["C"], f["b"], f["c"], f["rms"],
            o["b"], o["c"], o["rms"]))

    # Stability verdict
    print("\n  Stabilnosc wykladnikow (|db|,|dc|):")
    stable = True
    for cls in ["d", "sp", "s"]:
        f = fits.get(cls); o = R08.get(cls)
        if not (f and o): continue
        db = abs(f["b"] - o["b"]); dc = abs(f["c"] - o["c"])
        # threshold 0.2 dla b, 0.15 dla c (s-class bez c wiec pomin)
        if cls == "s":
            ok = db < 0.3
        else:
            ok = db < 0.25 and dc < 0.20
        if not ok: stable = False
        tag = "STABLE" if ok else "DRIFT"
        print("    {:<3s}: db={:+.3f}, dc={:+.3f}  [{}]".format(cls, f["b"]-o["b"], f["c"]-o["c"], tag))

    mean_prof = np.mean(prof_vals)
    print("\n  Mean profile RMS na N={}: {:.4f}".format(len(prof_vals), mean_prof))
    if mean_prof < 0.25 and stable:
        status = "WERYFIKACJA POZYTYWNA: class-specific D trzyma na N=34, paper v3 READY"
    elif mean_prof < 0.35 and stable:
        status = "WERYFIKACJA WARUNKOWA: wykladniki stabilne ale profile RMS rosnie"
    elif stable:
        status = "WERYFIKACJA MIESZANA: exponents stabilne, ale RMS za duze na paper"
    else:
        status = "NEGATIVE: exponents drifted - forma class-specific nie generalizuje"
    print("\n  " + status)

    print("""
  Nastepne kroki:
    - Jesli POZYTYWNE: przyklejac tabele do paperu v3, draftowa abstrakt
      i Methods section.
    - Jesli WARUNKOWE: analizowac outliery z nowych materialow (Zr, Sc, Y, La
      czesto maja problemy z czystoscia); moze extra ochrona.
    - Jesli NEGATYWNE: wrocic do class-by-class z jeszcze mniejszymi
      podklasami (grupa III, IV, VII, VIII, IX, itp.).
""")


if __name__ == "__main__":
    main()
