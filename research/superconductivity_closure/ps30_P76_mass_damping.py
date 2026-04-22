#!/usr/bin/env python3
# =============================================================================
#  ps30_P76_mass_damping.py
# -----------------------------------------------------------------------------
#  P7.6: Mass-damping korekta A_d dla ciezkich pierwiastkow d-klasy.
#
#  Problem z ps29c: Ta, W, Re, Os, Ir (pierwiastkowe 5d) sa
#  over-predicted 4x - 2000x. Przyczyna fizyczna:
#     lambda_ep ~ N(E_F) * <I^2> / (M * omega^2)
#  Ciezkie jadra z porownywalnym omega daja male lambda_ep, a nasza
#  Eq. 5 ma A_d = 0.310 uniwersalne bez mass-dependence.
#
#  Propozycja P7.6: A_d^eff(M) = A_d * (M_0/M)^(gamma_M / 2)
#  co daje T_c^eff = T_c_raw * (M_0/M)^gamma_M
#  (bo T_c ~ A^2)
#
#  Fit gamma_M na subset {Ta, Re, Os} (srodkowe 5d bez skrajnych anomalii
#  W, Ir — te moga wymagac dodatkowej korekty bandowej).
#
#  Test: czy jednym gamma_M mozna zamknac ~3 materialy? Jesli TAK, to
#  mamy nowa stala TGP z domenami (mass-damping mnoznik dla ciezkich d).
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import f as f_dist

K_B = 8.617333e-5
A_BOHR = 0.52917721067
A_STAR = 7.725 * A_BOHR
C_0 = 48.8222
SIGMA_A = 2.5856
ALPHA_P6B = 1.04
LAMBDA_0_P6B = 0.0962
OMEGA_0 = 15.0
BETA_P6D = 2.527
A_MAP = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a):
    n = max(1, round(a / A_STAR))
    d = a - n * A_STAR
    return np.exp(-d ** 2 / SIGMA_A ** 2)

def B_mag(lam): return 1.0 / (1.0 + BETA_P6D * lam)


def Tc_phonon(a, orb, z, omega, lam, mass_factor=1.0):
    """Tc_phonon z opcjonalnym mass_factor (P7.6) na A^2."""
    A = A_MAP[orb] if isinstance(orb, str) else float(orb) * A_MAP["d"]
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * (A ** 2) * mass_factor * M_gauss(a)
            * (Lambda_eff * 1e-3) / K_B * B_mag(lam))


# ==============================================================================
#  Dane: 3d/4d (baseline d) vs 5d (cele P7.6)
# ==============================================================================
# (name, T_obs, a, orb, z, omega, lam, M_atomic)
D_METALS = [
    # 3d/4d — fitowalne z A_d=0.310 uniwersalnie
    ("V",       5.30, 3.024, "d",  8, 31.0, 0.60,  51),
    ("Nb",      9.26, 3.301, "d",  8, 22.0, 0.20,  93),
    # 5d — cele P7.6
    ("Ta",      4.47, 3.300, "d",  8, 21.0, 0.00, 181),
    ("W",       0.012,3.160, "d",  8, 26.0, 0.00, 184),
    ("Re",      1.70, 2.760, "d", 12, 30.0, 0.00, 186),
    ("Os",      0.65, 2.730, "d", 12, 33.0, 0.00, 190),
    ("Ir",      0.11, 2.710, "d", 12, 30.0, 0.00, 192),
]

M_REF = 93.0  # Nb (referencyjny 4d)


def compute_Tpred(row, gamma_M=0.0):
    name, Tobs, a, orb, z, om, lam, M = row
    # mass_factor = (M_REF / M)^gamma_M (tylko gdy M > M_REF; inaczej 1)
    if M > M_REF and gamma_M > 0:
        mf = (M_REF / M) ** gamma_M
    else:
        mf = 1.0
    return Tc_phonon(a, orb, z, om, lam, mass_factor=mf)


def fit_gamma(rows_fit):
    """Fit gamma_M minimalizujac RMS(dlog) na podanym subsecie."""
    def loss(gamma):
        resid = []
        for r in rows_fit:
            Tp = compute_Tpred(r, gamma)
            resid.append(np.log10(Tp) - np.log10(r[1]))
        return float(np.sum(np.asarray(resid) ** 2))

    res = minimize_scalar(loss, bounds=(0.0, 10.0), method="bounded")
    return float(res.x), float(res.fun)


def main():
    print("=" * 78)
    print("  ps30_P76_mass_damping.py  (nowy subproblem P7.6)")
    print("  Mass-damping korekta A_d dla ciezkich d-klasowych pierwiastkow")
    print("=" * 78)

    # -------- baseline (gamma=0, brak korekty P7.6)
    print(f"\n  Baseline (bez P7.6, gamma_M = 0):")
    print(f"  {'name':<8s} {'M':>4s} {'T_obs':>7s} {'T_pred':>8s} "
          f"{'dlog':>7s}")
    for r in D_METALS:
        Tp = compute_Tpred(r, 0.0)
        dlog = np.log10(Tp) - np.log10(r[1])
        flag = ""
        if r[7] > M_REF and abs(dlog) > 0.5:
            flag = " <-- 5d OUT"
        print(f"  {r[0]:<8s} {r[7]:4d} {r[1]:7.2f} {Tp:8.2f} "
              f"{dlog:+7.3f}{flag}")

    # -------- Fit na Ta, Re, Os (srodkowe 5d)
    TARGETS_CORE = ["Ta", "Re", "Os"]
    TARGETS_FULL = ["Ta", "W", "Re", "Os", "Ir"]

    rows_core = [r for r in D_METALS if r[0] in TARGETS_CORE]
    rows_full = [r for r in D_METALS if r[0] in TARGETS_FULL]

    print("\n" + "-" * 78)
    print(f"  FIT gamma_M na subset = {TARGETS_CORE} (srodkowe 5d)")
    print("-" * 78)
    gamma_core, rss_core = fit_gamma(rows_core)
    rms_core = np.sqrt(rss_core / len(rows_core))
    print(f"  gamma_M = {gamma_core:.3f}   RMS(dlog) = {rms_core:.4f}")

    print("\n  Po P7.6 (gamma_M={:.3f}):".format(gamma_core))
    print(f"  {'name':<8s} {'M':>4s} {'T_obs':>7s} {'T_pred':>8s} "
          f"{'dlog':>7s}  {'w fit?':>7s}")
    for r in D_METALS:
        Tp = compute_Tpred(r, gamma_core)
        dlog = np.log10(Tp) - np.log10(r[1])
        inf = "FIT" if r[0] in TARGETS_CORE else ("HOLD" if r[7] > M_REF
                                                   else "ref")
        print(f"  {r[0]:<8s} {r[7]:4d} {r[1]:7.2f} {Tp:8.2f} "
              f"{dlog:+7.3f}  {inf:>7s}")

    # -------- Fit na pelny {Ta, W, Re, Os, Ir}
    print("\n" + "-" * 78)
    print(f"  FIT gamma_M na pelny subset = {TARGETS_FULL}")
    print("-" * 78)
    gamma_full, rss_full = fit_gamma(rows_full)
    rms_full = np.sqrt(rss_full / len(rows_full))
    print(f"  gamma_M = {gamma_full:.3f}   RMS(dlog) = {rms_full:.4f}")

    print("\n  Po P7.6 (gamma_M={:.3f}, fit na 5):".format(gamma_full))
    print(f"  {'name':<8s} {'M':>4s} {'T_obs':>7s} {'T_pred':>8s} "
          f"{'dlog':>7s}")
    for r in D_METALS:
        Tp = compute_Tpred(r, gamma_full)
        dlog = np.log10(Tp) - np.log10(r[1])
        print(f"  {r[0]:<8s} {r[7]:4d} {r[1]:7.2f} {Tp:8.2f} "
              f"{dlog:+7.3f}")

    # -------- LOO
    print("\n" + "-" * 78)
    print("  LOO na {Ta, Re, Os}:")
    print("-" * 78)
    for held in TARGETS_CORE:
        rows_loo = [r for r in rows_core if r[0] != held]
        g, _ = fit_gamma(rows_loo)
        hel = next(r for r in D_METALS if r[0] == held)
        Tp = compute_Tpred(hel, g)
        dlog = np.log10(Tp) - np.log10(hel[1])
        print(f"    held-out {held:<6s}  gamma_M_fit = {g:.3f}  "
              f"pred dlog = {dlog:+.3f}")

    # -------- Verdict
    print("\n" + "=" * 78)
    print("  P7.6 verdict")
    print("=" * 78)
    # max |dlog| na core po korekcie
    max_core = max(abs(np.log10(compute_Tpred(r, gamma_core))
                         - np.log10(r[1])) for r in rows_core)
    max_full = max(abs(np.log10(compute_Tpred(r, gamma_full))
                         - np.log10(r[1])) for r in rows_full)
    print(f"\n  fit na core {{Ta,Re,Os}}:  gamma_M = {gamma_core:.3f}  "
          f"max|dlog_core| = {max_core:.3f}")
    print(f"  fit na full {{Ta..Ir}}:    gamma_M = {gamma_full:.3f}  "
          f"max|dlog_full| = {max_full:.3f}")
    print(f"\n  Interpretacja:")
    if max_core < 0.3:
        print(f"    core-fit zamyka {{Ta,Re,Os}} w |dlog|<0.3 — OK, P7.6 dziala.")
    else:
        print(f"    core-fit NIE zamyka nawet srednich 5d — jeden gamma_M"
              f" nie wystarcza.")
    if max_full > 1.0:
        print(f"    full-fit NIE zamyka W i Ir — sa anomalnie ponad trend mass.")
        print(f"    To sugeruje dodatkowy czynnik band-structure/filling dla"
              f" W (d^4) i Ir (d^7).")
    else:
        print(f"    full-fit zamyka wszystkie 5 w |dlog|<1.0.")
    print(f"\n  Nowa stala TGP: gamma_M = {gamma_core:.3f} (z rekomendacja")
    print(f"  'domeny zastosowania: ciezkie d-pierwiastki M > {M_REF:.0f} amu')")


if __name__ == "__main__":
    main()
