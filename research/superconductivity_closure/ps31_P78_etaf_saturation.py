#!/usr/bin/env python3
# =============================================================================
#  ps31_P78_etaf_saturation.py
# -----------------------------------------------------------------------------
#  P7.8: Saturacja eta_f(Z) dla bardzo ciezkich aktynowcow.
#
#  Problem: P7.5a ma eta_f = -1.308 (liniowa poprawka 1+eta·Z^2/137^2).
#    Przy Z=90 (Th): Z^2/137^2 = 0.431, factor = (1-1.308·0.431)^2 = 0.191
#    ThH10: T_raw=284, T_corr=54, T_obs=161 -> dlog = -0.475 (OVER-CORRECTION)
#
#  Liniowe QED przyblizenie (Darwin ~ (Z alpha)^2) zalamuje sie gdy
#  (Z alpha)^2 > 0.3. Potrzeba saturacji.
#
#  Proponowane formy:
#     F-0 (linear, baseline):     (1 + eta_f · x)
#     F-1 (linear z wspol.):      (1 + eta_f^0 · (1 - beta · x) · x)
#     F-2 (tanh):                 (1 + eta_f^0 · x / (1 + beta·x^2))
#     F-3 (exp):                  exp(eta_f · x)                       [1 par.]
#  gdzie x = Z^2/137^2.
#
#  Fit na f-klasie z Z>=50: LaH10, CeH9, CeH10, La_amb, Th_amb, Ce_5GPa,
#  ThH10. ZeroWa korekta dla Z<50 (nieistotne).
#
#  Cel: forma ktora daje ThH10 dlog ~ 0 bez psucia La/Ce.
# =============================================================================
import numpy as np
from scipy.optimize import minimize, minimize_scalar

ALPHA2 = 1.0 / 137.036 ** 2

# f-klasa: (name, T_obs, T_raw, Z)
F_MATERIALS = [
    ("LaH10",   250.00, 368.08, 57),
    ("CeH9",    100.00, 143.13, 58),
    ("CeH10",   115.00, 187.48, 58),
    ("La_amb",    6.00,  14.34, 57),
    ("Th_amb",    1.38,  10.40, 90),
    ("Ce_5GPa",   1.70,   6.02, 58),
    ("ThH10",   161.00, 283.98, 90),   # nowy hold-out
]


# ==============================================================================
#  Formy eta_f(Z)
# ==============================================================================

def factor_F0(eta_f, Z):
    """Linear (P7.5a)."""
    x = Z ** 2 * ALPHA2
    return 1.0 + eta_f * x


def factor_F1(eta0, beta, Z):
    """Linear z saturacja: eta_eff(Z) = eta0 * (1 - beta * x)."""
    x = Z ** 2 * ALPHA2
    eta_eff = eta0 * (1.0 - beta * x)
    return 1.0 + eta_eff * x


def factor_F2(eta0, beta, Z):
    """tanh-like: eta_eff(Z) = eta0 / (1 + beta * x)."""
    x = Z ** 2 * ALPHA2
    eta_eff = eta0 / (1.0 + beta * x)
    return 1.0 + eta_eff * x


def factor_F3(eta, Z):
    """Exponential: factor = exp(eta * x)."""
    x = Z ** 2 * ALPHA2
    return np.exp(eta * x)


def apply_factor_to_Tc(T_raw, factor):
    """T_c ~ A^2, wiec T_new = T_raw * factor^2."""
    return T_raw * max(factor, 1e-6) ** 2


def rss_for_form(form_name, params, rows):
    """Zwroc sum of squares residuow log10."""
    total = 0.0
    for name, Tobs, Traw, Z in rows:
        if form_name == "F0":
            f = factor_F0(params[0], Z)
        elif form_name == "F1":
            f = factor_F1(params[0], params[1], Z)
        elif form_name == "F2":
            f = factor_F2(params[0], params[1], Z)
        elif form_name == "F3":
            f = factor_F3(params[0], Z)
        else:
            raise ValueError(form_name)
        Tnew = apply_factor_to_Tc(Traw, f)
        dlog = np.log10(Tnew) - np.log10(Tobs)
        total += dlog ** 2
    return total


def fit_form(form_name, rows, n_params):
    if n_params == 1:
        res = minimize_scalar(
            lambda p: rss_for_form(form_name, [p], rows),
            bounds=(-5.0, 5.0), method="bounded")
        return [float(res.x)], float(res.fun)
    else:
        best = (1e9, None)
        for e0 in [-3, -2, -1, -0.5, 0, 0.5, 1]:
            for b0 in [-3, -1, 0, 1, 3]:
                res = minimize(lambda p: rss_for_form(form_name, list(p), rows),
                               [e0, b0], method="Nelder-Mead",
                               options={"xatol": 1e-5, "fatol": 1e-8})
                if res.fun < best[0]:
                    best = (float(res.fun), list(res.x))
        return best[1], best[0]


def main():
    print("=" * 78)
    print("  ps31_P78_etaf_saturation.py  (nowy subproblem P7.8)")
    print("  Saturacja eta_f(Z) dla aktynowcow — fix ThH10")
    print("=" * 78)

    # ---------- baseline (F-0 linear, eta_f = -1.308)
    print("\n  Baseline F-0 (P7.5a, liniowa, eta_f = -1.308):")
    print(f"  {'name':<10s} {'Z':>3s} {'T_obs':>7s} {'T_raw':>8s} "
          f"{'T_corr':>8s} {'dlog':>7s}")
    for name, Tobs, Traw, Z in F_MATERIALS:
        f = factor_F0(-1.308, Z)
        Tnew = apply_factor_to_Tc(Traw, f)
        dlog = np.log10(Tnew) - np.log10(Tobs)
        flag = " <-- OVER" if dlog < -0.4 else (" <-- UNDER" if dlog > 0.4 else "")
        print(f"  {name:<10s} {Z:3d} {Tobs:7.2f} {Traw:8.2f} {Tnew:8.2f} "
              f"{dlog:+7.3f}{flag}")
    rss_F0 = rss_for_form("F0", [-1.308], F_MATERIALS)
    print(f"  RSS(F-0, eta=-1.308) = {rss_F0:.4f}  "
          f"RMS = {np.sqrt(rss_F0/len(F_MATERIALS)):.4f}")

    # ---------- refit F0 na wszystkich (z ThH10)
    print("\n" + "-" * 78)
    print("  F-0 refit (1 par.) na N=7 (z ThH10):")
    params0, rss0 = fit_form("F0", F_MATERIALS, 1)
    rms0 = np.sqrt(rss0 / len(F_MATERIALS))
    print(f"    eta_f = {params0[0]:+.4f}   RSS = {rss0:.4f}  RMS = {rms0:.4f}")
    show_fit("F0", params0, F_MATERIALS)

    # ---------- F1 linear correction
    print("\n" + "-" * 78)
    print("  F-1 (linear saturation, 2 par.):  eta(Z) = eta0 · (1 - beta·x)")
    params1, rss1 = fit_form("F1", F_MATERIALS, 2)
    rms1 = np.sqrt(rss1 / len(F_MATERIALS))
    print(f"    eta0 = {params1[0]:+.4f}  beta = {params1[1]:+.4f}"
          f"   RSS = {rss1:.4f}  RMS = {rms1:.4f}")
    show_fit("F1", params1, F_MATERIALS)

    # ---------- F2 tanh-like (actually Lorentzian saturation)
    print("\n" + "-" * 78)
    print("  F-2 (Lorentzian saturation, 2 par.):  eta(Z) = eta0 / (1 + beta·x)")
    params2, rss2 = fit_form("F2", F_MATERIALS, 2)
    rms2 = np.sqrt(rss2 / len(F_MATERIALS))
    print(f"    eta0 = {params2[0]:+.4f}  beta = {params2[1]:+.4f}"
          f"   RSS = {rss2:.4f}  RMS = {rms2:.4f}")
    show_fit("F2", params2, F_MATERIALS)

    # ---------- F3 exponential
    print("\n" + "-" * 78)
    print("  F-3 (exponential, 1 par.):  factor = exp(eta · x)")
    params3, rss3 = fit_form("F3", F_MATERIALS, 1)
    rms3 = np.sqrt(rss3 / len(F_MATERIALS))
    print(f"    eta = {params3[0]:+.4f}   RSS = {rss3:.4f}  RMS = {rms3:.4f}")
    show_fit("F3", params3, F_MATERIALS)

    # ---------- comparison
    print("\n" + "=" * 78)
    print("  Porownanie form na N=7 (f-klasa, Z>=50):")
    print("=" * 78)
    print(f"  {'forma':<30s} {'RSS':>6s} {'RMS':>6s} {'parametry':<s}")
    for lab, rss, rms, p, k in [
        ("F-0 linear (original -1.308)", rss_F0, np.sqrt(rss_F0/7), "-1.308", 0),
        ("F-0 linear (refit)",   rss0, rms0, f"{params0[0]:+.3f}", 1),
        ("F-1 linear saturation", rss1, rms1, f"eta0={params1[0]:+.3f} beta={params1[1]:+.3f}", 2),
        ("F-2 Lorentzian sat.",   rss2, rms2, f"eta0={params2[0]:+.3f} beta={params2[1]:+.3f}", 2),
        ("F-3 exponential",       rss3, rms3, f"{params3[0]:+.3f}", 1),
    ]:
        print(f"  {lab:<30s} {rss:6.3f} {rms:6.3f}  {p}")

    # ---------- verdict
    print("\n  LOO test dla F-1 (2-par) i F-0 (1-par) — hold-out Th_amb i ThH10")
    for held in ["Th_amb", "ThH10"]:
        rows_fit = [r for r in F_MATERIALS if r[0] != held]
        p0, _ = fit_form("F0", rows_fit, 1)
        p1, _ = fit_form("F1", rows_fit, 2)
        held_row = next(r for r in F_MATERIALS if r[0] == held)
        f0 = factor_F0(p0[0], held_row[3])
        f1 = factor_F1(p1[0], p1[1], held_row[3])
        T0 = apply_factor_to_Tc(held_row[2], f0)
        T1 = apply_factor_to_Tc(held_row[2], f1)
        d0 = np.log10(T0) - np.log10(held_row[1])
        d1 = np.log10(T1) - np.log10(held_row[1])
        print(f"    held {held:<10s}  F-0: eta={p0[0]:+.3f}, "
              f"held dlog = {d0:+.3f}")
        print(f"    held {held:<10s}  F-1: eta0={p1[0]:+.3f} beta={p1[1]:+.3f}, "
              f"held dlog = {d1:+.3f}")


def show_fit(form_name, params, rows):
    print(f"  {'name':<10s} {'Z':>3s} {'T_obs':>7s} {'T_corr':>8s} "
          f"{'dlog':>7s}")
    for name, Tobs, Traw, Z in rows:
        if form_name == "F0":
            f = factor_F0(params[0], Z)
        elif form_name == "F1":
            f = factor_F1(params[0], params[1], Z)
        elif form_name == "F2":
            f = factor_F2(params[0], params[1], Z)
        elif form_name == "F3":
            f = factor_F3(params[0], Z)
        Tnew = apply_factor_to_Tc(Traw, f)
        dlog = np.log10(Tnew) - np.log10(Tobs)
        print(f"    {name:<10s} {Z:3d} {Tobs:7.2f} {Tnew:8.2f} {dlog:+7.3f}")


if __name__ == "__main__":
    main()
