#!/usr/bin/env python3
# =============================================================================
#  r01_bg_baseline.py
# -----------------------------------------------------------------------------
#  Baseline: klasyczny Bloch-Gruneisen fit dla 15 metali.
#
#  Wzor BG-Wilson:
#     rho_ph(T) = R * (T/Theta_D)^5 * J_5(Theta_D/T)
#     J_5(y)    = int_0^y x^5 / [(exp(x)-1)(1 - exp(-x))] dx
#
#  Total: rho(T) = rho_0 + rho_ph(T)
#
#  Testy:
#    (A) Fit (rho_0, R, Theta_D) na 4 pkt T — czy fit Theta_D zgadza sie z lit?
#    (B) Fix Theta_D = lit, fit (rho_0_fit, R) — czy R jest uniwersalne?
#    (C) Residuum per klasa materialow (noble / SC-elem / FM / Stoner / hcp)
#
#  Expected: R ~ specific to material (zalezy od lambda_ep), Theta_D zbiega
#  w 10%. Outliery: Fe (magnon scatter), Pb (bardzo nisko Theta_D).
# =============================================================================
import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "r00", os.path.join(os.path.dirname(__file__), "r00_dataset.py"))
r00 = importlib.util.module_from_spec(spec); sys.modules["r00"] = r00
spec.loader.exec_module(r00)


# -----------------------------------------------------------------------------
# Bloch-Gruneisen function
# -----------------------------------------------------------------------------
def _integrand(x):
    """x^5 / [(exp(x)-1)(1-exp(-x))] — stabilny dla x>0."""
    if x < 1e-6:
        return x**3  # limit for small x: x^5 / x^2 = x^3
    e = np.exp(x)
    return x**5 / ((e - 1.0) * (1.0 - 1.0/e))


def J5(y):
    """Integral from 0 to y of _integrand. Cache-able."""
    if y <= 0: return 0.0
    if y > 50: y = 50  # zatrzymanie dla bardzo malego T/Theta_D (saturation ~124.4)
    val, _ = quad(_integrand, 0.0, y, limit=200)
    return float(val)


def rho_bg(T, rho_0, R, Theta_D):
    """Bloch-Gruneisen: rho(T) = rho_0 + R * (T/Th)^5 * J_5(Th/T)."""
    if T <= 0 or Theta_D <= 0: return rho_0
    y = Theta_D / T
    return rho_0 + R * (T / Theta_D)**5 * J5(y)


def rho_bg_vec(Ts, rho_0, R, Theta_D):
    return np.array([rho_bg(T, rho_0, R, Theta_D) for T in Ts])


# -----------------------------------------------------------------------------
# Fitting routines
# -----------------------------------------------------------------------------
def fit_free(d):
    """Fit (rho_0, R, Theta_D) na 4 T-points (77, 295, 500, 1000)."""
    Ts = np.array(r00.T_POINTS)
    rhos = np.array([d["rho"][T] for T in r00.T_POINTS])
    rho_0_guess = d["rho_0"]
    Th_guess = d["Theta_D"]
    R_guess = d["rho"][295] * 4  # rough
    try:
        popt, pcov = curve_fit(
            rho_bg_vec, Ts, rhos,
            p0=[rho_0_guess, R_guess, Th_guess],
            bounds=([0.0, 0.1, 50.0], [1.0, 1000.0, 1000.0]),
            maxfev=5000,
        )
        # predictions at all T + rho_0
        pred = rho_bg_vec(Ts, *popt)
        rms = np.sqrt(np.mean((np.log10(pred) - np.log10(rhos))**2))
        return popt, rms
    except Exception as e:
        return None, 99.0


def fit_R_only(d, Theta_D_fixed=None):
    """Fix (rho_0 = d['rho_0'], Theta_D = lit), fit R tylko."""
    Ts = np.array(r00.T_POINTS)
    rhos = np.array([d["rho"][T] for T in r00.T_POINTS])
    Th = Theta_D_fixed if Theta_D_fixed else d["Theta_D"]
    rho_0 = d["rho_0"]

    def model(Ts, R):
        return np.array([rho_bg(T, rho_0, R, Th) for T in Ts])

    R_guess = d["rho"][295] * 4
    try:
        popt, pcov = curve_fit(model, Ts, rhos, p0=[R_guess],
                               bounds=([0.1], [1000.0]), maxfev=5000)
        pred = model(Ts, *popt)
        rms = np.sqrt(np.mean((np.log10(pred) - np.log10(rhos))**2))
        return float(popt[0]), rms, rho_0, Th
    except Exception:
        return None, 99.0, rho_0, Th


# -----------------------------------------------------------------------------
# Klasy materialow (z r00)
# -----------------------------------------------------------------------------
CLASSES = {
    "noble":    {"Cu", "Ag", "Au"},
    "sp":       {"Al", "Pb", "Sn"},
    "sc_bcc":   {"Nb", "V"},
    "FM":       {"Fe", "Ni"},
    "stoner":   {"Pt", "Pd"},
    "hcp_sp":   {"Cd", "Zn", "Mg"},
}


def get_class(name):
    for c, s in CLASSES.items():
        if name in s: return c
    return "?"


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 92)
    print("  r01_bg_baseline.py  -  Bloch-Gruneisen fit dla 15 metali")
    print("=" * 92)

    data = r00.load_dataset()

    # Test A: fit (rho_0, R, Theta_D) free
    print("\n" + "-" * 92)
    print("  (A) Free fit:  (rho_0, R, Theta_D)")
    print("-" * 92)
    print("  {:<5s} {:>7s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>7s} {:>8s}".format(
        "name", "Th_lit", "Th_fit", "%dTh", "R_fit",
        "rho_0_l", "rho_0_f", "cls", "RMS_log"))
    th_diffs = []
    R_free = {}
    th_fit_all = {}
    for d in data:
        popt, rms = fit_free(d)
        cls = get_class(d["name"])
        if popt is None:
            print("  {:<5s}  FIT FAILED".format(d["name"]))
            continue
        rho0_f, R_f, Th_f = popt
        R_free[d["name"]] = R_f
        th_fit_all[d["name"]] = Th_f
        d_pct = (Th_f - d["Theta_D"]) / d["Theta_D"] * 100
        th_diffs.append(d_pct)
        print("  {:<5s} {:>7.0f} {:>8.0f} {:>+7.1f}% {:>8.2f} {:>8.4f} {:>8.4f} {:>7s} {:>8.4f}".format(
            d["name"], d["Theta_D"], Th_f, d_pct, R_f,
            d["rho_0"], rho0_f, cls, rms))

    print("\n  Th_D % drift statistics:")
    print("    mean |dTh|% = {:.1f}".format(np.mean(np.abs(th_diffs))))
    print("    max  |dTh|% = {:.1f}".format(np.max(np.abs(th_diffs))))

    # Test B: fix Theta_D = lit, fit R only
    print("\n" + "-" * 92)
    print("  (B) Fix rho_0, Theta_D = literatura, fit R only")
    print("-" * 92)
    print("  {:<5s} {:>7s} {:>10s} {:>8s} {:>8s} {:>8s}".format(
        "name", "Th_D", "R_fit", "cls", "RMS_log", "R/rho(295)"))
    R_fixed = {}
    rms_list = []
    for d in data:
        R, rms, rho0, Th = fit_R_only(d)
        if R is None: continue
        R_fixed[d["name"]] = R
        rms_list.append(rms)
        cls = get_class(d["name"])
        ratio = R / d["rho"][295]
        print("  {:<5s} {:>7.0f} {:>10.2f} {:>8s} {:>8.4f} {:>8.3f}".format(
            d["name"], d["Theta_D"], R, cls, rms, ratio))

    print("\n  Statystyki R_fixed:")
    R_vals = list(R_fixed.values())
    print("    median R =", round(np.median(R_vals), 2))
    print("    mean R   =", round(np.mean(R_vals), 2),
          "+/-", round(np.std(R_vals), 2))
    print("    min/max  =", round(min(R_vals), 2), "/", round(max(R_vals), 2))
    print("    CoV      = {:.1%}".format(np.std(R_vals)/np.mean(R_vals)))

    # Test C: residual analysis per class
    print("\n" + "-" * 92)
    print("  (C) R_fixed po klasach materialow")
    print("-" * 92)
    for cls, names in CLASSES.items():
        Rs = [R_fixed[n] for n in names if n in R_fixed]
        if not Rs: continue
        rho295 = [next(d for d in data if d["name"] == n)["rho"][295]
                  for n in names if n in R_fixed]
        ratio = [R / r for R, r in zip(Rs, rho295)]
        print("  {:<8s}: N={}  R=[{}]  R/rho(295)=[{}]".format(
            cls, len(Rs),
            ", ".join("{:.1f}".format(r) for r in Rs),
            ", ".join("{:.2f}".format(r) for r in ratio)))

    # Test D: check relacja R vs rho(295). Klasyczne: R ~ 4.225 * rho(Theta_D)
    # Dla T >> Theta_D: rho(T) ~ R * T / (4 * Theta_D)
    # Czyli R ~ 4 * Theta_D * drho/dT (high-T slope)
    print("\n" + "-" * 92)
    print("  (D) Klasyczny test: R predykcja z slope rho(T) w T = 1000K")
    print("-" * 92)
    print("  High-T limit BG: drho/dT ~ R / (4*Theta_D) dla T >> Theta_D")
    print("  Wiec R_slope = 4 * Theta_D * drho/dT (oczekujemy ~R_fit)")
    print()
    print("  {:<5s} {:>8s} {:>8s} {:>10s} {:>10s} {:>6s}".format(
        "name", "Th_D", "R_fit", "R_slope", "R_s/R_f", "cls"))
    for d in data:
        if d["name"] not in R_fixed: continue
        slope = (d["rho"][1000] - d["rho"][500]) / 500
        R_slope = 4 * d["Theta_D"] * slope
        R_f = R_fixed[d["name"]]
        print("  {:<5s} {:>8.0f} {:>8.2f} {:>10.2f} {:>10.3f} {:>6s}".format(
            d["name"], d["Theta_D"], R_f, R_slope, R_slope / R_f,
            get_class(d["name"])))

    # Test E: 1-parametrowy globalny fit na wszystkich danych z fixed rho_0, Theta_D
    print("\n" + "-" * 92)
    print("  (E) Globalny fit: jedno R dla wszystkich metali?")
    print("-" * 92)
    # Skalar: rho_predicted(T_i, M_i) = rho_0(M_i) + R_global * f(T, Theta_D_M)
    def global_rms(R_g):
        resid = []
        for d in data:
            for T in r00.T_POINTS:
                pred = rho_bg(T, d["rho_0"], R_g, d["Theta_D"])
                obs = d["rho"][T]
                resid.append(np.log10(pred) - np.log10(obs))
        return np.sqrt(np.mean(np.asarray(resid)**2))
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(global_rms, bounds=(1, 500), method="bounded")
    print("  Globalny R opt =", round(res.x, 2))
    print("  Globalny RMS_log =", round(res.fun, 4))
    print("  (per-material avg RMS =", round(np.mean(rms_list), 4), ")")
    print("  Wniosek: jesli globalny >> per-material -> R nie jest uniwersalne")

    # Verdict
    print("\n" + "=" * 92)
    print("  VERDICT r01")
    print("=" * 92)
    print("""
    (A) Fit free Theta_D:
        mean |dTh|% = {:.1f}, max = {:.1f}
        {}

    (B) Fit R (fixed Theta_D):
        median R = {:.2f}, CoV = {:.1%}
        {}

    (D) R_slope vs R_fit:
        jesli stosunek ~ 1 dla wiekszosci to BG jest dobrym modelem high-T.

    (E) Globalny R:
        {:.2f} (RMS {:.3f}) vs per-material RMS {:.3f}
        {}
""".format(
        np.mean(np.abs(th_diffs)), np.max(np.abs(th_diffs)),
        "OK" if np.mean(np.abs(th_diffs)) < 15 else "DRYFT",
        np.median(R_vals), np.std(R_vals)/np.mean(R_vals),
        "R nie jest uniwersalne (CoV > 30%)" if np.std(R_vals)/np.mean(R_vals) > 0.3
            else "R blisko uniwersalne",
        res.x, res.fun, np.mean(rms_list),
        "jedno R nie wystarcza" if res.fun > 2*np.mean(rms_list)
            else "jedno R dziala do pewnego stopnia"
    ))


if __name__ == "__main__":
    main()
