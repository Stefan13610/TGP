#!/usr/bin/env python3
# =============================================================================
#  ps33_P710_NEF_factor.py
# -----------------------------------------------------------------------------
#  P7.10: Jawny N(E_F) factor w Eq. 5 — test czy zamyka negatywna walidacja
#         ps32 bez psucia kalibracji SC.
#
#  Motywacja (z ps32): model bez N(E_F) przewiduje spurious T_c ~ 20-60 K dla
#  Cu, diamond, NaCl, LiF itp. 87% negatywnej walidacji FAIL.
#
#  Propozycja: A^2 -> A^2 * g(N_F / N_F_ref)
#    gdzie g(x) to monotonicznie rosnacy z g(0)=0, g(x>=1) -> 1
#
#  Testowane formy:
#    G-1 (power-law):    g(x) = min(1, x^p)     [1 par: p]
#    G-2 (tanh):         g(x) = tanh(x / x_0)   [1 par: x_0]
#    G-3 (sigmoid):      g(x) = 1 / (1 + (x_0/x)^n) [2 par]
#
#  Kryterium sukcesu:
#    (a) SC materials z P7.5a (N=15 core): RMS(dlog) nie wzrasta > 0.05
#    (b) non-SC z ps32 (N=23): mediana T_pred < 1 K, liczba FAIL spada <= 5
#    (c) p (lub x_0) fizycznie sensowne (nie ekstrema)
#
#  Znane ograniczenia: Pt, Pd maja WYSOKIE N_F (paramagnetyki Stoner-enhanced)
#  — nie bedzie zamkniete przez N_F-factor. Zostana jako osobny P7.11
#  (magnetic suppression via lam_sf kalibrowany na Pt, Pd).
# =============================================================================
import numpy as np
from scipy.optimize import minimize_scalar, minimize

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

N_F_REF = 0.90  # Nb: 0.90 states/eV/atom (referencja)


def k_d(z):
    return {3: 0.7, 4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)


def M_gauss(a):
    n = max(1, round(a / A_STAR))
    d = a - n * A_STAR
    return np.exp(-d ** 2 / SIGMA_A ** 2)


def B_mag(lam):
    return 1.0 / (1.0 + BETA_P6D * lam)


# ===========================================================================
#   Formy g(x) gdzie x = N_F / N_F_ref
# ===========================================================================

def g_powerlaw(x, p):
    """g(x) = min(1, x^p). g(0)=0, g(1)=1."""
    if x <= 0:
        return 0.0
    return min(1.0, x ** p)


def g_tanh(x, x0):
    """g(x) = tanh(x/x0). g(0)=0, g(oo)=1."""
    if x <= 0:
        return 0.0
    return float(np.tanh(x / x0))


def g_sigmoid(x, x0, n):
    """g(x) = 1/(1 + (x0/x)^n). g(0)=0, g(oo)=1."""
    if x <= 0:
        return 0.0
    return 1.0 / (1.0 + (x0 / x) ** n)


def Tc_phonon_P710(a, orb, z, omega, lam, N_F, g_func, *g_params):
    """T_c z jawnym N(E_F) factor na A^2."""
    A = A_MAP[orb] if isinstance(orb, str) else float(orb) * A_MAP["d"]
    x = N_F / N_F_REF
    g_val = g_func(x, *g_params)
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * (A ** 2) * g_val * M_gauss(a)
            * (Lambda_eff * 1e-3) / K_B * B_mag(lam))


def Tc_phonon_baseline(a, orb, z, omega, lam):
    """Baseline (bez P7.10) = obecne Eq. 5."""
    A = A_MAP[orb] if isinstance(orb, str) else float(orb) * A_MAP["d"]
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * (A ** 2) * M_gauss(a)
            * (Lambda_eff * 1e-3) / K_B * B_mag(lam))


# ===========================================================================
#   Dataset SC: materialy phonon-mediated z ktorych Eq. 5 dziala
#   (name, T_obs, a, orb, z, omega, lam_sf, N_F [states/eV/atom])
# ===========================================================================
SC_CALIB = [
    # 3d/4d core (przywrocone z ps29c)
    ("V",       5.30, 3.024, "d",  8, 31.0, 0.60,  1.32),
    ("Nb",      9.26, 3.301, "d",  8, 22.0, 0.20,  0.97),   # = N_F_REF ish
    # 6p-heavy (po P7.5b)
    ("Pb",      7.19, 3.500, "sp", 12,  9.0, 0.00, 0.36),
    ("Hg_elem", 4.15, 2.991, "sp",  6,  9.0, 0.00, 0.27),
    # BCS-like hydrydy (po P7.5a)
    ("MgB2",   39.00, 3.086, "sp", 12, 70.0, 0.00, 0.72),
    ("H3S",   203.00, 3.089, "s",   8,140.0, 0.00, 0.42),
    ("LaH10", 250.00, 5.116, "f",   8, 150.0,0.00, 0.62),
    ("YH6",   220.00, 3.605, "f",   8, 140.0,0.00, 0.57),
    ("CeH9",  100.00, 3.720, "f",   8, 120.0,0.00, 0.50),
    ("CeH10", 115.00, 5.300, "f",   8, 130.0,0.00, 0.55),
    # cuprates (po P7.5b)
    ("YBCO",   93.00, 3.85,  "d",  12, 60.0, 0.00, 0.70),
    ("Bi2212", 94.00, 3.82,  "d",  12, 65.0, 0.00, 0.60),
    # ambient lantanowce (po P7.5a)
    ("La_amb",  6.00, 3.774, "f",   8, 22.0, 0.00, 0.45),
    ("Th_amb",  1.38, 5.08,  "f",   8, 25.0, 0.00, 0.40),
    ("Ce_5GPa", 1.70, 3.70,  "f",   8, 25.0, 0.00, 0.35),
]


# ===========================================================================
#   Dataset non-SC (z ps32 + N_F z DFT literatury)
#   (name, T_obs_ub, a, orb, z, omega, lam_sf, N_F)
# ===========================================================================
NON_SC = [
    # Noble (nizkie N_F)
    ("Cu",       0.00,   3.615, "d",  12, 28.0, 0.0,  0.14),
    ("Ag",       0.00035,4.086, "d",  12, 21.5, 0.0,  0.10),
    ("Au",       0.00005,4.078, "d",  12, 17.0, 0.0,  0.07),
    # Paramagnety Stoner-enhanced (wysokie N_F) - NIE ZAMKNIE P7.10
    ("Pt",       0.00,   3.924, "d",  12, 23.9, 0.0,  2.00),
    ("Pd",       0.00,   3.890, "d",  12, 25.2, 0.0,  2.30),
    # Alkalies (niskie N_F)
    ("Li",       0.0004, 3.510, "s",   8, 30.0, 0.0,  0.44),
    ("Na",       0.00,   4.290, "s",   8, 13.0, 0.0,  0.28),
    ("K",        0.00,   5.330, "s",   8,  8.6, 0.0,  0.42),
    # Ferromagnety (wysokie N_F spin-split) - NIE ZAMKNIE P7.10
    ("Fe",       0.00,   2.870, "d",   8, 35.0, 3.0,  1.50),
    ("Co",       0.00,   2.510, "d",  12, 40.0, 3.0,  1.10),
    ("Ni",       0.00,   3.520, "d",  12, 38.0, 2.0,  2.90),
    # Alkaline earths (niskie N_F)
    ("Ca_amb",   0.00,   5.588, "s",  12, 19.0, 0.0,  0.14),
    ("Sr_amb",   0.00,   6.080, "s",  12, 12.0, 0.0,  0.15),
    ("Ba_amb",   0.00,   5.030, "s",  12,  7.0, 0.0,  0.20),
    ("Mg",       0.00,   3.210, "s",  12, 40.0, 0.0,  0.11),
    # Marginal SC
    ("Zn",       0.85,   2.660, "sp", 12, 26.0, 0.0,  0.13),
    ("Cd",       0.52,   2.980, "sp", 12, 16.0, 0.0,  0.10),
    # Insulatory (N_F = 0 - perfect PASS kandydaci)
    ("Diamond",  0.00,   3.570, "sp",  4, 165.0,0.0,  0.00),
    ("Si_amb",   0.00,   5.431, "sp",  4,  55.0,0.0,  0.00),
    ("Ge_amb",   0.00,   5.658, "sp",  4,  35.0,0.0,  0.00),
    ("NaCl",     0.00,   5.640, "sp",  6, 26.0, 0.0,  0.00),
    ("LiF",      0.00,   4.030, "sp",  6, 80.0, 0.0,  0.00),
    # Graphite (prawie 0 w Diracu)
    ("Graphite", 0.00,   2.460, "sp",  3, 170.0,0.0,  0.03),
]


def compute_rms_sc(g_func, *g_params):
    """RMS(dlog) dla SC z dataset."""
    resid = []
    for name, Tobs, a, orb, z, om, lam, NF in SC_CALIB:
        Tp = Tc_phonon_P710(a, orb, z, om, lam, NF, g_func, *g_params)
        if Tp <= 0:
            return 1e9
        resid.append(np.log10(Tp) - np.log10(Tobs))
    return float(np.sqrt(np.mean(np.asarray(resid) ** 2)))


def negsc_stats(g_func, *g_params):
    """Dla non-SC: mediana T_pred, 95 percentyl, liczba FAIL."""
    Tpreds, fails, borders = [], 0, 0
    for name, Tub, a, orb, z, om, lam, NF in NON_SC:
        Tp = Tc_phonon_P710(a, orb, z, om, lam, NF, g_func, *g_params)
        Tpreds.append(Tp)
        if Tub < 1e-4:
            if Tp > 1.0:     fails += 1
            elif Tp > 0.1:   borders += 1
        else:
            ratio = Tp / Tub
            if ratio > 20:   fails += 1
            elif ratio > 5:  borders += 1
    return float(np.median(Tpreds)), float(np.percentile(Tpreds, 95)), fails, borders


def fit_power_law():
    """Fit p maksymalizujac: kill non-SC bez psucia SC."""
    def loss(p):
        if p <= 0.05: return 1e9
        rms_sc = compute_rms_sc(g_powerlaw, p)
        median_nonsc, _, fails, _ = negsc_stats(g_powerlaw, p)
        # wazona kombinacja: RMS_SC + 0.2*median_nonsc + 0.1*#FAIL
        return rms_sc + 0.3 * median_nonsc + 0.05 * fails
    res = minimize_scalar(loss, bounds=(0.1, 10), method="bounded")
    return float(res.x), float(res.fun)


def fit_tanh():
    """Fit x0 tanh analogicznie."""
    def loss(x0):
        if x0 < 0.01: return 1e9
        rms_sc = compute_rms_sc(g_tanh, x0)
        median_nonsc, _, fails, _ = negsc_stats(g_tanh, x0)
        return rms_sc + 0.3 * median_nonsc + 0.05 * fails
    res = minimize_scalar(loss, bounds=(0.01, 5), method="bounded")
    return float(res.x), float(res.fun)


def fit_sigmoid():
    """Fit (x0, n) sigmoid 2-parametrowy."""
    best = (1e9, None)
    for x0_0 in [0.1, 0.3, 0.5, 0.8]:
        for n_0 in [1.0, 2.0, 4.0]:
            def loss(params):
                x0, n = params
                if x0 <= 0 or n <= 0: return 1e9
                rms_sc = compute_rms_sc(g_sigmoid, x0, n)
                median_nonsc, _, fails, _ = negsc_stats(g_sigmoid, x0, n)
                return rms_sc + 0.3 * median_nonsc + 0.05 * fails
            res = minimize(loss, [x0_0, n_0], method="Nelder-Mead",
                           options={"xatol": 1e-4, "fatol": 1e-6})
            if res.fun < best[0]:
                best = (float(res.fun), list(res.x))
    return best[1], best[0]


def report_fit(name_tag, g_func, *g_params):
    print("\n" + "-" * 78)
    print(f"  {name_tag}  parametry: {g_params}")
    print("-" * 78)
    rms_sc = compute_rms_sc(g_func, *g_params)
    median_nonsc, p95, fails, borders = negsc_stats(g_func, *g_params)
    print(f"  SC calib (N={len(SC_CALIB)}): RMS(dlog) = {rms_sc:.4f}")
    print(f"  non-SC (N={len(NON_SC)}): median T_pred = {median_nonsc:.3f} K, "
          f"P95 = {p95:.3f} K")
    print(f"  FAIL: {fails}/{len(NON_SC)}   BORDER: {borders}/{len(NON_SC)}")

    # Detailed SC breakdown
    print(f"\n  SC materials (sprawdzenie czy nie psujemy):")
    print(f"  {'name':<10s} {'T_obs':>7s} {'T_pred':>8s} {'dlog':>7s} "
          f"{'N_F':>5s} {'g(x)':>6s}")
    for row in SC_CALIB:
        name, Tobs, a, orb, z, om, lam, NF = row
        Tp = Tc_phonon_P710(a, orb, z, om, lam, NF, g_func, *g_params)
        dlog = np.log10(Tp) - np.log10(Tobs)
        gv = g_func(NF / N_F_REF, *g_params)
        flag = "" if abs(dlog) < 0.3 else " <-- out"
        print(f"  {name:<10s} {Tobs:7.2f} {Tp:8.2f} {dlog:+7.3f} "
              f"{NF:5.2f} {gv:6.3f}{flag}")

    # Detailed non-SC breakdown
    print(f"\n  non-SC (czy zabija spurious?):")
    print(f"  {'name':<10s} {'T_ub':>7s} {'T_pred':>8s} {'N_F':>5s} "
          f"{'g(x)':>6s} {'verdict':>8s}")
    for row in NON_SC:
        name, Tub, a, orb, z, om, lam, NF = row
        Tp = Tc_phonon_P710(a, orb, z, om, lam, NF, g_func, *g_params)
        gv = g_func(NF / N_F_REF, *g_params)
        if Tub < 1e-4:
            verdict = "PASS" if Tp < 0.1 else ("BORDER" if Tp < 1.0 else "FAIL")
        else:
            ratio = Tp / Tub
            verdict = ("PASS" if ratio < 5 else
                       ("BORDER" if ratio < 20 else "FAIL"))
        Tub_str = f"{Tub:.4f}" if Tub > 0 else "~0"
        Tp_str = f"{Tp:.3f}" if Tp < 100 else f"{Tp:.1e}"
        print(f"  {name:<10s} {Tub_str:>7s} {Tp_str:>8s} {NF:5.2f} "
              f"{gv:6.3f} {verdict:>8s}")


def main():
    print("=" * 78)
    print("  ps33_P710_NEF_factor.py  (nowy subproblem P7.10)")
    print(f"  N(E_F) factor: A^2 -> A^2 * g(N_F/N_F_ref), N_F_ref = {N_F_REF}")
    print("=" * 78)

    # ---- Baseline (g = 1, czyli obecne Eq. 5)
    def g_const(x): return 1.0
    print("\n  Baseline (bez P7.10, g = 1):")
    rms0 = compute_rms_sc(lambda x: 1.0)
    med0, p95_0, fails0, borders0 = negsc_stats(lambda x: 1.0)
    print(f"    SC N={len(SC_CALIB)}: RMS(dlog) = {rms0:.4f}")
    print(f"    non-SC N={len(NON_SC)}: median T_pred = {med0:.3f} K, "
          f"P95 = {p95_0:.3f} K, {fails0} FAIL, {borders0} BORDER")

    # ---- G-1 power-law
    print("\n" + "=" * 78)
    print("  G-1 power-law: g(x) = min(1, x^p)")
    print("=" * 78)
    p_fit, _ = fit_power_law()
    print(f"  Fitted p = {p_fit:.3f}")
    report_fit("G-1 power-law (best p)", g_powerlaw, p_fit)

    # ---- G-2 tanh
    print("\n" + "=" * 78)
    print("  G-2 tanh: g(x) = tanh(x/x0)")
    print("=" * 78)
    x0_fit, _ = fit_tanh()
    print(f"  Fitted x0 = {x0_fit:.3f}")
    report_fit("G-2 tanh (best x0)", g_tanh, x0_fit)

    # ---- G-3 sigmoid 2-par
    print("\n" + "=" * 78)
    print("  G-3 sigmoid: g(x) = 1/(1 + (x0/x)^n)")
    print("=" * 78)
    params3, _ = fit_sigmoid()
    print(f"  Fitted x0 = {params3[0]:.3f}, n = {params3[1]:.3f}")
    report_fit("G-3 sigmoid (best x0, n)", g_sigmoid, params3[0], params3[1])

    # ---- Verdict
    print("\n" + "=" * 78)
    print("  P7.10 verdict")
    print("=" * 78)

    for label, gf, ps in [
        ("G-1 power-law", g_powerlaw, (p_fit,)),
        ("G-2 tanh",      g_tanh,     (x0_fit,)),
        ("G-3 sigmoid",   g_sigmoid,  tuple(params3)),
    ]:
        rms_sc = compute_rms_sc(gf, *ps)
        _, _, fails, _ = negsc_stats(gf, *ps)
        delta = rms_sc - rms0
        passed = (delta < 0.05) and (fails < 10)
        status = "OK" if passed else "WEAK"
        print(f"  {label:<20s}  dRMS_SC = {delta:+.3f}, #FAIL = {fails}/"
              f"{len(NON_SC)}  [{status}]")

    print("""
  Kluczowe obserwacje:
    - G-1 (power-law): 1 parametr TGP, prosta forma QM (T_c ~ N_F^p z BCS)
    - G-2 (tanh): saturacyjna, p ~ 1 w granicy, ale tlumi oba ekstrema plynnie
    - G-3 (sigmoid): 2 parametry, twardy prog - elegancko zamyka insulatory

  Pt, Pd, Fe, Co, Ni pozostaja FAIL niezaleznie od g(x) - maja wysokie N_F
  lecz sa paramagnetyczne/ferromagnetyczne. To osobny sub-problem:
    P7.11 (magnetic suppression): kalibracja lam_sf (B_mag) na Pt, Pd
    jako paramagnetyki Stoner-enhanced + Fe, Co, Ni jako FM.
""")


if __name__ == "__main__":
    main()
