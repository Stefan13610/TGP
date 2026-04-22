#!/usr/bin/env python3
# =============================================================================
#  r02_tgp_formula.py
# -----------------------------------------------------------------------------
#  Kluczowy test: czy R (BG amplitude) == f(TGP constants z paperu SC)?
#
#  Strategia:
#    (1) Importuj TGP constants z ps29c (C_0, A_orb per class, ...).
#    (2) Zaklad: R = C_R * k_d(z) * A_orb^2 * M_gauss(a) * Lambda_eff * g(N_F)
#       gdzie C_R jest jedna nowa stala TGP, a reszta jest identyczna jak w
#       predykcji T_c (paper v1 Eq. 5).
#    (3) Dla kazdego materialu policz TGP_factor, fit liniowy R = C_R * factor.
#    (4) Zrewiduj residuum po klasach.
#
#  Jesli R = C_R * (Eq. 5 structure) zbiega w RMS < 0.3 -> TGP ma predykcje
#  transportu, model jest spojny z normalnym stanem.
#
#  Jesli nie -> musimy wprowadzic alfa_tr (transport-weighted) niezaleznie od
#  lambda_ep z SC.
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import pearsonr
import importlib.util, os, sys

# r00 dataset + r01 BG
spec = importlib.util.spec_from_file_location(
    "r00", os.path.join(os.path.dirname(__file__), "r00_dataset.py"))
r00 = importlib.util.module_from_spec(spec); sys.modules["r00"] = r00
spec.loader.exec_module(r00)

spec1 = importlib.util.spec_from_file_location(
    "r01", os.path.join(os.path.dirname(__file__), "r01_bg_baseline.py"))
r01 = importlib.util.module_from_spec(spec1); sys.modules["r01"] = r01
spec1.loader.exec_module(r01)


# -----------------------------------------------------------------------------
# TGP constants (z paperu SC v1 + P7.5)
# -----------------------------------------------------------------------------
# A_orb klasy:
A_ORB = {
    "s":  -0.111,  # czysto s (alkalies)
    "sp":  0.207,  # sp-metals (Al, Pb, Sn, ...)
    "d":   0.310,  # d-band (Nb, V, Fe, Ni, Pt, Pd)
    "f":   2.034,  # f-band (lantanowce, aktynowce)
}

# P7.5 relativistic corrections
ETA_F = -1.308   # 4f
ETA_6P = 0.331   # 6p (Pb, Tl, Bi)
ALPHA = 137.036  # fine-structure

# Paper v1 constants
C_0 = 48.8222            # TGP prefactor
A_STAR = 4.68            # Gaussian center (A)
SIGMA = 2.15             # Gaussian width (A)
LAMBDA_0 = 0.175         # substrate coupling (meV at ref omega)
OMEGA_0 = 30.0           # ref omega (meV)
ALPHA_P6B = 0.5          # freq-dep exponent

# P7.10/P7.12 N_F factor
N_F_REF = 0.90
ALPHA_712 = 0.10         # ps37 value (soft deg)
BETA_712 = 2.15


# -----------------------------------------------------------------------------
# TGP factor functions (z ps29c)
# -----------------------------------------------------------------------------
def k_d(z):
    """Coordination factor k_d(z) ~ ln(z+1) / ln(13)."""
    return np.log(z + 1.0) / np.log(13.0)


def M_gauss(a, a_star=A_STAR, sigma=SIGMA):
    """Lattice-matching Gaussian."""
    return float(np.exp(-(a - a_star)**2 / (2 * sigma**2)))


def Lambda_eff(omega_meV, Lambda_0=LAMBDA_0, omega_0=OMEGA_0, alpha=ALPHA_P6B):
    """Frequency-dependent substrate coupling."""
    return Lambda_0 * (omega_meV / omega_0)**alpha


def g_stretched(x, alpha=ALPHA_712, beta=BETA_712):
    """P7.12 stretched-exp g(x), x = N_F/N_F_ref."""
    if x <= 1e-6: return 0.0
    if x >= 1.0: return 1.0
    return float(np.exp(-alpha * (1.0/x - 1.0) ** beta))


def theta_D_to_omega_meV(Theta_D_K):
    """Theta_D [K] -> omega [meV]. k_B * Theta_D / hbar."""
    # k_B = 0.08617 meV/K
    return 0.08617 * Theta_D_K


# -----------------------------------------------------------------------------
# Orbital class assignment dla 15 metali
# -----------------------------------------------------------------------------
ORB_CLASS = {
    # Noble: d11 s1 -- d decyduje (band blisko E_F) ale full-d daje s/sp
    # Cu, Ag, Au: d-band ponizej E_F, nosniki z s-p -> klasa sp
    "Cu":  "sp", "Ag":  "sp", "Au":  "sp",
    # sp metals - czyste
    "Al":  "sp", "Pb":  "sp", "Sn":  "sp",
    # bcc TM - d-band dominuje
    "Nb":  "d",  "V":   "d",
    # FM: d-band
    "Fe":  "d",  "Ni":  "d",
    # Stoner paramagnetyki: d-band
    "Pt":  "d",  "Pd":  "d",
    # hcp sp
    "Cd":  "sp", "Zn":  "sp", "Mg":  "sp",
}

# Nie ma f-metali w pilocie, nie ma 6p. P7.5a/b nie aktywne.


def A_orb_eff(name):
    """Aktywna klasa A_orb z ewentualna P7.5a/b korekta."""
    cls = ORB_CLASS[name]
    A_base = A_ORB[cls]
    # Pb (Z=82) moze miec lekka 6p SOC korekte, ale ~0.3% -> pomijamy
    return A_base


# -----------------------------------------------------------------------------
# TGP factor dla ρ(T): f_TGP(material)
# -----------------------------------------------------------------------------
def tgp_factor_rho(d, z_coord, use_g=True):
    """
    F_TGP = k_d(z) * A_orb^2 * M_gauss(a) * Lambda_eff(omega) * g(N_F/N_F_ref)

    Uzywamy omega = k_B * Theta_D, traktujac Theta_D jako input.
    Zwracamy liczbe - F_TGP(material).
    """
    kd = k_d(z_coord)
    Aorb = A_orb_eff(d["name"])
    Aorb2 = Aorb ** 2
    Mg = M_gauss(d["a_latt"])
    omega_meV = theta_D_to_omega_meV(d["Theta_D"])
    Le = Lambda_eff(omega_meV)
    F = kd * Aorb2 * Mg * Le
    if use_g:
        x = d["N_F"] / N_F_REF
        F *= g_stretched(x)
    return F


# Koordynacja z struktury krystalograficznej
COORD = {
    "Cu":  12, "Ag":  12, "Au":  12,  # fcc
    "Al":  12, "Pb":  12,              # fcc
    "Sn":   6,                          # tetr beta
    "Nb":   8, "V":    8,              # bcc
    "Fe":   8, "Ni":  12,              # Fe bcc, Ni fcc
    "Pt":  12, "Pd":  12,              # fcc
    "Cd":  12, "Zn":  12, "Mg":  12,   # hcp ~ 12
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 92)
    print("  r02_tgp_formula.py  -  czy R(BG) == C_R * f_TGP(material)?")
    print("=" * 92)

    data = r00.load_dataset()
    R_fit = {}
    for d in data:
        R, rms, _, _ = r01.fit_R_only(d)
        if R is not None:
            R_fit[d["name"]] = R

    # Compute TGP factor for each
    rows = []
    for d in data:
        if d["name"] not in R_fit: continue
        z = COORD[d["name"]]
        F_with_g = tgp_factor_rho(d, z, use_g=True)
        F_no_g = tgp_factor_rho(d, z, use_g=False)
        R = R_fit[d["name"]]
        rows.append((d["name"], d["Theta_D"], d["a_latt"], d["N_F"],
                     z, ORB_CLASS[d["name"]], R, F_with_g, F_no_g))

    # Tabela
    print("\n  {:<5s} {:>7s} {:>7s} {:>5s} {:>4s} {:>4s} {:>10s} {:>14s} {:>14s}".format(
        "name", "Theta_D", "a[A]", "N_F", "z", "cls", "R_BG", "F_TGP (g)", "F_TGP (no g)"))
    for name, Th, a, NF, z, cls, R, Fg, F0 in rows:
        print("  {:<5s} {:>7.0f} {:>7.3f} {:>5.2f} {:>4d} {:>4s} {:>10.3f} {:>14.5f} {:>14.5f}".format(
            name, Th, a, NF, z, cls, R, Fg, F0))

    # Compute ratio R / F_TGP = effective C_R
    print("\n  R / F_TGP(with g):")
    print("  {:<5s} {:>10s} {:>14s} {:>12s} {:>5s}".format(
        "name", "R", "F_TGP", "R/F_TGP", "cls"))
    ratios_g = []
    for name, Th, a, NF, z, cls, R, Fg, F0 in rows:
        if Fg < 1e-10: continue
        rr = R / Fg
        ratios_g.append((name, rr, cls))
        print("  {:<5s} {:>10.3f} {:>14.5f} {:>12.2e} {:>5s}".format(
            name, R, Fg, rr, cls))

    vals = [r[1] for r in ratios_g]
    print("\n  Statystyki R/F_TGP (with g):")
    print("    median =", "{:.2e}".format(np.median(vals)))
    print("    mean   =", "{:.2e}".format(np.mean(vals)),
          "+/-", "{:.2e}".format(np.std(vals)))
    print("    min/max=", "{:.2e}".format(min(vals)),
          "/", "{:.2e}".format(max(vals)))
    print("    CoV    = {:.1%}".format(np.std(vals)/np.mean(vals)))

    # Per klasa
    print("\n  R/F_TGP per klasa:")
    classes = set(r[2] for r in ratios_g)
    for cls in classes:
        vs = [r[1] for r in ratios_g if r[2] == cls]
        ns = [r[0] for r in ratios_g if r[2] == cls]
        print("    {:<7s} (N={}): median = {:.2e}, vals = [{}]".format(
            cls, len(vs),
            np.median(vs),
            ", ".join("{} {:.1e}".format(n, v) for n, v in zip(ns, vs))))

    # Correlation test
    print("\n  Correlation R vs F_TGP:")
    R_arr = np.array([r[6] for r in rows])
    F_arr = np.array([r[7] for r in rows])  # F_TGP with g
    F_nog = np.array([r[8] for r in rows])  # F_TGP no g
    r_log_g, p_val_g = pearsonr(np.log10(R_arr), np.log10(F_arr + 1e-20))
    r_log_0, p_val_0 = pearsonr(np.log10(R_arr), np.log10(F_nog + 1e-20))
    print("    log(R) vs log(F_TGP with g):   r = {:.3f}, p = {:.4f}".format(
        r_log_g, p_val_g))
    print("    log(R) vs log(F_TGP without g): r = {:.3f}, p = {:.4f}".format(
        r_log_0, p_val_0))

    # Alternative 1-param fit
    print("\n  Fit C_R globalny: R = C_R * F_TGP (log-log)")
    def obj(logC):
        C = 10**logC
        resid = [np.log10(C * F) - np.log10(R)
                 for R, F in zip(R_arr, F_arr + 1e-20)]
        return np.sqrt(np.mean(np.array(resid)**2))
    res = minimize_scalar(obj, bounds=(-5, 10), method="bounded")
    C_R_opt = 10**res.x
    print("    C_R (F with g)  = {:.3e}".format(C_R_opt))
    print("    RMS_log         = {:.4f}".format(res.fun))

    def obj0(logC):
        C = 10**logC
        resid = [np.log10(C * F) - np.log10(R)
                 for R, F in zip(R_arr, F_nog + 1e-20)]
        return np.sqrt(np.mean(np.array(resid)**2))
    res0 = minimize_scalar(obj0, bounds=(-5, 10), method="bounded")
    C_R_0 = 10**res0.x
    print("    C_R (F no g)    = {:.3e}".format(C_R_0))
    print("    RMS_log         = {:.4f}".format(res0.fun))

    # Per-material predykcja
    print("\n  Predykcja R_TGP (with g, C_R = {:.2e}):".format(C_R_opt))
    print("  {:<5s} {:>10s} {:>10s} {:>8s} {:>5s}".format(
        "name", "R_obs", "R_TGP", "dlog", "cls"))
    for r in rows:
        name, _, _, _, _, cls, R, Fg, _ = r
        R_pred = C_R_opt * Fg
        dlog = np.log10(R_pred + 1e-20) - np.log10(R)
        flag = " <- out" if abs(dlog) > 0.5 else ""
        print("  {:<5s} {:>10.3f} {:>10.3f} {:>+8.3f} {:>5s}{}".format(
            name, R, R_pred, dlog, cls, flag))

    # Verdict
    print("\n" + "=" * 92)
    print("  VERDICT r02")
    print("=" * 92)
    print("""
    Hipotezy testowane:
      H1:  R = C_R * (k_d * A_orb^2 * M_gauss * Lambda_eff * g(N_F))
           r_log = {:.3f}  (p = {:.4f})
           C_R = {:.2e}, RMS_log = {:.4f}

      H0:  R = C_R * (k_d * A_orb^2 * M_gauss * Lambda_eff)  [bez g]
           r_log = {:.3f}  (p = {:.4f})
           C_R = {:.2e}, RMS_log = {:.4f}

    {}

    Nastepny krok:
      - r03: szersza analiza per-T residuum (czy BG tail/T^5 jest reprodukowany)
      - r04: dodac inne formy TGP_factor (np. lambda_tr = lambda_ep * (1 -<cos theta>))
      - r05: rozszerzyc dataset o Mo, Ta, W, Be, Ti
""".format(
        r_log_g, p_val_g, C_R_opt, res.fun,
        r_log_0, p_val_0, C_R_0, res0.fun,
        "H1 lepsze" if res.fun < res0.fun else "H0 lepsze"
    ))


if __name__ == "__main__":
    main()
