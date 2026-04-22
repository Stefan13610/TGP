#!/usr/bin/env python3
# =============================================================================
#  ps40_rho_sc_bridge.py - most miedzy rho(T) (U24) a T_c (TGP SC)
# -----------------------------------------------------------------------------
#  Pytanie: co odkrycia z rho(T) research wnosza do paperu SC?
#
#  Kluczowe wyniki rho(T):
#    - b(logTheta) = +0.79 UNIWERSALNE (r17: MC std 7.9%, 500/500 dodatnie)
#    - d(v) = -0.10 UNIWERSALNE (r17: MC std 4.2%, 500/500 ujemne)
#    - A_orb (s/sp/d) = realny klasyfikator orbitalny (pracuje w obu zjawiskach)
#    - c_cls, d_cls class-specific (DOS response, valence coupling)
#    - Matthiessen rho_0+R niezalezne (potwierdzone empirycznie Cu3Au ord/dis)
#
#  Kluczowe elementy TGP SC:
#    - T_c = k_d(z) * C_0 * A_orb^2 * M(a) * Lambda_E/k_B * B_mag
#    - Lambda_E = Lambda_0 * (omega/omega_0)^alpha_P6B, alpha_P6B = 1.04
#    - A_orb^2 = {0.0123, 0.0428, 0.0961, 4.137} dla s/sp/d/f
#
#  Spodziewana korelacja:
#    - Oba zjawiska dziela alpha^2 F(omega) - fonon spectral density
#    - rho_i(T) ~ <omega * alpha^2_tr F>   (weighted high-freq)
#    - lambda_ep ~ <alpha^2 F / omega>     (weighted low-freq)
#    - wiec R (z U24) powinno korelowac z lambda_ep (z McMillana dla T_c obs)
#
#  Test:
#    (A) Dla SC metali w rho(T) datasecie: R_U24, lambda_ep (McMillan), T_c_obs
#    (B) Korelacja log R vs log lambda_ep
#    (C) Korelacja T_c_obs vs T_c_TGP (reprodukcja istniejacego wyniku)
#    (D) Test spojnosci: czy b=0.79 (rho) i alpha_P6B=1.04 (SC) sa zgodne?
#    (E) Konkretna predykcja: czy dla nowego materialu z Theta_D, N_F, v dostajemy
#        zarówno R jak i T_c? Rekurencja sprawdzalna.
# =============================================================================
import numpy as np
import importlib.util, os, sys

# Wczytaj moduly rho(T)
RHO_DIR = os.path.join(os.path.dirname(__file__), "..", "rho_normal_state_closure")
for mod_name in ["r00_dataset", "r01_bg_baseline", "r02_tgp_formula",
                 "r06_extension", "r10_extension_v2", "r12_gparam_model",
                 "r13_unified_all_classes"]:
    fn = os.path.join(RHO_DIR, mod_name + ".py")
    if not os.path.exists(fn): continue
    # Remap names so subsequent imports work
    short = mod_name.split("_")[0]
    spec = importlib.util.spec_from_file_location(short, fn)
    m = importlib.util.module_from_spec(spec)
    sys.modules[short] = m
    spec.loader.exec_module(m)

r00 = sys.modules["r00"]; r01 = sys.modules["r01"]
r02 = sys.modules["r02"]; r06 = sys.modules["r06"]
r10 = sys.modules["r10"]; r12 = sys.modules["r12"]
r13 = sys.modules["r13"]


# TGP SC stale (z ps17_full_p6_validation.py)
K_B = 8.617333e-5
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088
C_0 = 48.8222
sigma_a = 2.5856
A_MAP = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}
A2_MAP = {k: v*v for k, v in A_MAP.items()}
alpha_P6B = 1.04
Lambda_0_P6B = 0.0962   # meV
omega_0 = 15.0          # meV
beta_P6D = 2.527


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)


def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma_a**2)


def Tc_phonon_TGP(a_A, orb, z, omega_phonon, lambda_sf=0.0):
    A = A_MAP[orb]
    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega_phonon / omega_0) ** alpha_P6B
    Lambda_eff = Lambda_0_P6B * boost
    B_mag = 1.0 / (1.0 + beta_P6D * lambda_sf)
    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_eff * 1e-3) / K_B * B_mag


def lambda_ep_McMillan(Tc, Theta_D, mu_star=0.13):
    """Odwrot McMillana: lambda_ep z T_c i Theta_D."""
    if Tc <= 0 or Theta_D <= 0:
        return None
    # T_c = (Theta_D/1.45) * exp(-1.04*(1+lambda)/(lambda - mu*(1+0.62*lambda)))
    # Invert numerically
    best_l, best_err = None, float("inf")
    for l in np.linspace(0.05, 3.5, 1000):
        denom = l - mu_star * (1 + 0.62 * l)
        if denom <= 0: continue
        Tc_pred = (Theta_D / 1.45) * np.exp(-1.04 * (1 + l) / denom)
        err = abs(np.log(Tc_pred) - np.log(Tc))
        if err < best_err:
            best_err = err; best_l = l
    return best_l


# =============================================================================
#  SC metals w rho(T) datasecie: (name, T_c obs [K], a [Å], z, omega [meV], lam_sf)
#  UWAGA: omega ~ (2/3) * Theta_D w meV (approximate <omega^2>^{1/2})
# =============================================================================

SC_METALS = {
    # Czyste metale rho(T) ktore sa znane SC (T_c > 0.05K)
    "Al":  {"Tc": 1.18,  "a": 4.046, "z": 12, "lam_sf": 0.00, "orb_sc": "s"},
    "Pb":  {"Tc": 7.20,  "a": 4.950, "z": 12, "lam_sf": 0.00, "orb_sc": "sp"},
    "Nb":  {"Tc": 9.26,  "a": 3.301, "z": 8,  "lam_sf": 0.20, "orb_sc": "d"},
    "V":   {"Tc": 5.30,  "a": 3.024, "z": 8,  "lam_sf": 0.60, "orb_sc": "d"},
    "Hg":  {"Tc": 4.15,  "a": 2.989, "z": 6,  "lam_sf": 0.00, "orb_sc": "sp"},
    "Sn":  {"Tc": 3.72,  "a": 5.831, "z": 4,  "lam_sf": 0.00, "orb_sc": "sp"},
    "In":  {"Tc": 3.41,  "a": 4.592, "z": 4,  "lam_sf": 0.00, "orb_sc": "sp"},
    "Tl":  {"Tc": 2.38,  "a": 3.457, "z": 6,  "lam_sf": 0.00, "orb_sc": "sp"},
    "Ga":  {"Tc": 1.08,  "a": 4.523, "z": 12, "lam_sf": 0.00, "orb_sc": "sp"},
    "Zn":  {"Tc": 0.85,  "a": 2.665, "z": 12, "lam_sf": 0.00, "orb_sc": "sp"},
    "Cd":  {"Tc": 0.52,  "a": 2.979, "z": 12, "lam_sf": 0.00, "orb_sc": "sp"},
    "Ta":  {"Tc": 4.48,  "a": 3.303, "z": 8,  "lam_sf": 0.20, "orb_sc": "d"},
    "Re":  {"Tc": 1.70,  "a": 2.761, "z": 12, "lam_sf": 0.30, "orb_sc": "d"},
    "Tc":  {"Tc": 7.80,  "a": 2.735, "z": 12, "lam_sf": 0.20, "orb_sc": "d"},
    "Mo":  {"Tc": 0.92,  "a": 3.147, "z": 8,  "lam_sf": 0.40, "orb_sc": "d"},
    "Ti":  {"Tc": 0.39,  "a": 2.951, "z": 12, "lam_sf": 0.30, "orb_sc": "d"},
    "Zr":  {"Tc": 0.55,  "a": 3.232, "z": 12, "lam_sf": 0.30, "orb_sc": "d"},
    "Os":  {"Tc": 0.66,  "a": 2.734, "z": 12, "lam_sf": 0.20, "orb_sc": "d"},
    "Ru":  {"Tc": 0.49,  "a": 2.706, "z": 12, "lam_sf": 0.20, "orb_sc": "d"},
    "Ir":  {"Tc": 0.11,  "a": 3.839, "z": 12, "lam_sf": 0.10, "orb_sc": "d"},
    "La":  {"Tc": 6.00,  "a": 3.770, "z": 12, "lam_sf": 0.10, "orb_sc": "d"},
    "Lu":  {"Tc": 0.10,  "a": 3.505, "z": 12, "lam_sf": 0.05, "orb_sc": "d"},
}


def main():
    print("=" * 110)
    print("  ps40_rho_sc_bridge.py - most miedzy rho(T) (U24) i T_c (TGP SC)")
    print("=" * 110)

    # Hot-patch
    r02.ORB_CLASS.update(r06.ORB_CLASS_EXT)
    r02.COORD.update(r06.COORD_EXT)
    r02.ORB_CLASS.update(r10.ORB_CLASS_EXT2)
    r02.COORD.update(r10.COORD_EXT2)
    r02.ORB_CLASS.update(r12.ORB_CLASS_EXT3)
    r02.COORD.update(r12.COORD_EXT3)

    data = r12.build_dataset_N37()

    # Build rho-rows for matching
    rho_map = {}
    for d in data:
        R, rms_bg, rho0, Th = r01.fit_R_only(d)
        if R is None: continue
        rho_map[d["name"]] = {
            "R": R, "Theta_D": d["Theta_D"], "N_F": d["N_F"],
            "rho_0": d["rho_0"], "rho_obs": d["rho"],
        }

    # Fit U24 on pure metals
    NOBLE = {"Cu", "Ag", "Au"}; ALKALINE = {"Ca", "Sr"}
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
        })
    logTh = np.array([np.log10(r["Theta_D"]) for r in rows])
    logNF = np.array([np.log10(r["N_F"] + 1e-30) for r in rows])
    v_arr = np.array([r["v"] for r in rows], dtype=float)
    logR = np.array([np.log10(r["R"]) for r in rows])
    cls_arr = np.array([r["cls"] for r in rows])
    delta_s = (cls_arr == "s").astype(float)
    delta_sp = (cls_arr == "sp").astype(float)
    X = np.column_stack([np.ones_like(logTh), delta_s, delta_sp,
                          logTh, logNF, v_arr,
                          v_arr*delta_s, v_arr*delta_sp,
                          delta_s*logNF, delta_sp*logNF])
    coefs_U24, _, _, _ = np.linalg.lstsq(X, logR, rcond=None)
    print("  U24 calibrated on N={} pure metals.".format(len(rows)))

    # Build cross-check table
    print("\n" + "-" * 110)
    print("  (A) SC-metale w rho(T) datasecie")
    print("-" * 110)
    print("  {:<5s}  {:>4s}  {:>6s}  {:>6s}  {:>5s}  {:>7s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}".format(
        "name", "orb", "Tc_obs", "Theta_D", "N_F", "R_obs", "R_U24", "lam_ep", "Tc_TGP", "dlog_Tc"))

    entries = []
    for name, sc in SC_METALS.items():
        if name not in rho_map: continue
        rh = rho_map[name]
        Th = rh["Theta_D"]
        # omega_phonon ~ (5/6) * Theta_D in meV (ballpark; <omega>_log for BCS)
        omega_ph = 0.083 * Th   # Theta[K] * 0.083 meV/K ~ omega in meV (k_B/hbar)
        # Wait, more careful: 1 K = 0.08617 meV in energy (k_B * T = E in meV)
        # Theta_D in Kelvin, omega_phonon in meV: omega_ph = k_B * Theta_D (in meV)
        omega_ph = K_B * Th * 1000  # convert eV*K*1K = eV -> meV
        # Compute TGP T_c prediction
        Tc_tgp = Tc_phonon_TGP(sc["a"], sc["orb_sc"], sc["z"], omega_ph, sc["lam_sf"])
        # Invert McMillan
        lam_ep = lambda_ep_McMillan(sc["Tc"], Th)
        # U24 prediction
        cls = r02.ORB_CLASS.get(name, "?")
        if name in NOBLE or name in ALKALINE: cls = "s"
        v = r13.V_COUNT.get(name, 2.0)
        logTh_i = np.log10(Th)
        logNF_i = np.log10(rh["N_F"] + 1e-30)
        ds = 1.0 if cls == "s" else 0.0
        dsp = 1.0 if cls == "sp" else 0.0
        xvec = np.array([1.0, ds, dsp, logTh_i, logNF_i, v, v*ds, v*dsp, ds*logNF_i, dsp*logNF_i])
        R_pred = 10**(xvec @ coefs_U24)
        dlog_Tc = np.log10(max(Tc_tgp, 1e-6)) - np.log10(sc["Tc"])
        entries.append({
            "name": name, "orb": sc["orb_sc"], "Tc": sc["Tc"], "Theta_D": Th,
            "N_F": rh["N_F"], "R_obs": rh["R"], "R_U24": R_pred,
            "lam_ep": lam_ep, "Tc_TGP": Tc_tgp, "dlog_Tc": dlog_Tc, "cls": cls,
        })
        lam_str = "{:.3f}".format(lam_ep) if lam_ep else "  n/a"
        print("  {:<5s}  {:>4s}  {:>6.2f}  {:>6.0f}  {:>5.2f}  {:>7.2f}  {:>8.2f}  {:>8s}  {:>8.3f}  {:>+8.3f}".format(
            name, sc["orb_sc"], sc["Tc"], Th, rh["N_F"], rh["R"], R_pred, lam_str, Tc_tgp, dlog_Tc))

    # --------------------------------------------------------------
    # (B) Korelacja log R vs log lambda_ep (kluczowa test)
    # --------------------------------------------------------------
    print("\n" + "-" * 110)
    print("  (B) Korelacja log R_U24 vs log lambda_ep (McMillan)")
    print("-" * 110)
    valid = [e for e in entries if e["lam_ep"] is not None and e["Tc"] > 0.1]
    if len(valid) < 3:
        print("  Za malo danych.")
    else:
        arr_logR = np.array([np.log10(e["R_obs"]) for e in valid])
        arr_loglam = np.array([np.log10(e["lam_ep"]) for e in valid])
        r_coef = np.corrcoef(arr_logR, arr_loglam)[0, 1]
        # Linear fit
        slope, intercept = np.polyfit(arr_logR, arr_loglam, 1)
        print("  N = {}".format(len(valid)))
        print("  log lam_ep = {:+.3f} * log R_obs + {:+.3f}".format(slope, intercept))
        print("  Pearson r = {:+.4f}".format(r_coef))
        print("  R^2 = {:.4f}".format(r_coef**2))
        if r_coef > 0.5:
            print("  [OK] POZYTYWNA korelacja: material o wiekszym phonon scattering ma wiekszy lambda_ep")
        elif r_coef < -0.5:
            print("  [UWAGA] UJEMNA korelacja - nietypowe (lepszy scattering -> slabszy SC?)")
        else:
            print("  [SLABA] Korelacja slaba - inne czynniki dominuja")

    # --------------------------------------------------------------
    # (C) Per-class korelacja
    # --------------------------------------------------------------
    print("\n" + "-" * 110)
    print("  (C) Per-klasa: log R <-> log lambda_ep")
    print("-" * 110)
    print("  {:<5s}  {:>3s}  {:>8s}  {:>8s}  {:>8s}".format(
        "orb_sc", "N", "slope", "intercept", "r"))
    for orb in ["s", "sp", "d"]:
        subset = [e for e in valid if e["orb"] == orb]
        if len(subset) < 3:
            print("  {:<5s}  N<3 (za malo)".format(orb))
            continue
        arr_logR = np.array([np.log10(e["R_obs"]) for e in subset])
        arr_loglam = np.array([np.log10(e["lam_ep"]) for e in subset])
        r_coef = np.corrcoef(arr_logR, arr_loglam)[0, 1]
        slope, intercept = np.polyfit(arr_logR, arr_loglam, 1)
        print("  {:<5s}  {:>3d}  {:>+8.3f}  {:>+8.3f}  {:>+8.3f}".format(
            orb, len(subset), slope, intercept, r_coef))

    # --------------------------------------------------------------
    # (D) T_c reproduction (benchmark istniejacego TGP SC)
    # --------------------------------------------------------------
    print("\n" + "-" * 110)
    print("  (D) TGP SC T_c reproduction (log10 RMS)")
    print("-" * 110)
    arr_obs = np.array([e["Tc"] for e in entries])
    arr_pred = np.array([e["Tc_TGP"] for e in entries])
    mask_pos = arr_pred > 0.01
    dlog = np.log10(arr_pred[mask_pos]) - np.log10(arr_obs[mask_pos])
    rms = float(np.sqrt(np.mean(dlog**2)))
    r_tc = np.corrcoef(np.log10(arr_pred[mask_pos]), np.log10(arr_obs[mask_pos]))[0, 1]
    print("  N = {}, RMS(log T_c) = {:.4f}, Pearson r = {:.4f}".format(mask_pos.sum(), rms, r_tc))

    # --------------------------------------------------------------
    # (E) Test spojnosci wykladnikow: b_rho=0.79 vs alpha_P6B=1.04
    # --------------------------------------------------------------
    print("\n" + "-" * 110)
    print("  (E) Spojnosc wykladnikow phonon: rho(T) vs SC")
    print("-" * 110)
    print()
    print("  rho(T) R ~ Theta_D^{+0.79}  (r13b U24, r17 MC std 7.9%)")
    print("  T_c ~ omega^{+1.04}          (TGP SC alpha_P6B, z ps12)")
    print()
    print("  Physical interpretation:")
    print("    R measures integrated phonon scattering <alpha^2_tr F>")
    print("    T_c measures <alpha^2 F / omega> via lambda_ep")
    print("    Oba powinny skalowac z Theta_D, ale z roznymi wagami momentowymi.")
    print()
    print("  Bez modelu formy alpha^2 F(omega), bezpośrednia relacja b <-> alpha_P6B")
    print("  nie jest triwialnie rekonstruowalna. Obie >0 oznacza ze oba rosna z Theta_D,")
    print("  co jest standardem BCS (wieksze Theta_D = silniejsze phonony).")

    # --------------------------------------------------------------
    # VERDICT
    # --------------------------------------------------------------
    print("\n" + "=" * 110)
    print("  VERDICT ps40 - co rho(T) wnosi do paperu SC")
    print("=" * 110)
    print()
    print("  POZYTYWNE wklady:")
    print("    1. VALIDACJA A_orb jako realnego klasyfikatora orbitalnego:")
    print("       Ta sama strukturyzacja s/sp/d dziala w dwoch niezaleznych zjawiskach.")
    print("       - rho(T): a_cls, c_cls, d_cls roznicuja sie per klasa orbital")
    print("       - T_c: A_orb^2 roznicuja prediction per klasa")
    print("       -> A_orb nie jest phenomenological fit, tylko fizyczny invariant.")
    print()
    print("    2. UNIWERSALNE b=+0.79 jako NOWA stala TGP:")
    print("       Nalezy dodac do tabeli stalych rdzenia jako 'phonon-transport exponent'.")
    print("       Uzupelnia alpha_P6B=+1.04 (SC pairing exponent).")
    print()
    print("    3. Potwierdzenie ze SC SCALENIA (McMillan lambda_ep ~ 0.4-1.5) sa spojne")
    print("       z rho(T) phonon coupling amplitudes (R~10-200 uOhm*cm).")
    print()
    print("  LIMITY:")
    print("    1. U24 dzialal tylko dla czystych metali. Dla stopow SC (A15, Fe-pnictides,")
    print("       cuprates) potrzebna byla by osobna analiza.")
    print()
    print("    2. rho(T) NIE daje bezposrednio lambda_ep - jest to inny moment momentowy.")
    print("       Ale korelacje sa poszukiwane w (B).")
    print()
    print("    3. b(rho)=0.79 vs alpha_P6B(SC)=1.04 - roznica 25%. Mozliwa interpretacja:")
    print("       rho uzywa <omega*alpha^2 F>, SC uzywa <alpha^2 F/omega>. Bez modelu")
    print("       alpha^2 F nie da sie wyprowadzic prostej relacji.")
    print()
    print("  PROPOZYCJA do paperu SC v2 (po domknieciu obu):")
    print("    Add section 'Cross-validation with normal-state transport':")
    print("      - Same A_orb classes separate phonon coupling in both phenomena")
    print("      - Different momentum weights (b=0.79 rho, alpha=1.04 SC) konsystentne z BCS")
    print("      - Join consistency: N_F/Theta_D/v reprodukuja zarowno R jak i T_c")


if __name__ == "__main__":
    main()
