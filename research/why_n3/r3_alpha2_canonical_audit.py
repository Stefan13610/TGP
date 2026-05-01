#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_alpha2_canonical_audit.py
=============================

PURPOSE
-------
Bezstronny audyt R3 z **TGP-kanonicznym alpha=2** (z K(phi) = K_geo * phi^4
w sek08a / TGP_FOUNDATIONS:56).

Wczesniejsze skrypty R3 (r3_metric_singularity, r3_alpha_scan, r3_atail_bridge)
uzywaly alpha=1 jako "substrat" mimo ze sek08a definiuje akcje z czynnikiem
kinetycznym K(phi) = K_geo * phi^4. W radialnej notacji L = (1/2) g^{2*alpha} (g')^2,
tozsamosc K(phi) = phi^4 daje 2*alpha = 4, czyli alpha = 2.

Ten skrypt sprawdza wszystkie kluczowe wyniki R3 ponownie z alpha=2:
  1. g0_crit(alpha=2, d=3)
  2. Czy phi-drabinka (g0^mu = phi*g0^e) miesci sie pod bariera
  3. Czy g0^tau z Koide (1.729) miesci sie pod bariera
  4. Czy A_tail^4 = m_mu/m_e dziala dla alpha=2
  5. Czy "mechanizm A^4" (m = c * K^2) dziala dla alpha=2

KLUCZOWA OBSERWACJA Z R3 README (lin. 690):
  alpha=1.0: g0_crit(3D) = 2.206  (substrat R3)
  alpha=2.0: g0_crit(3D) = 1.874  (TGP-canonical z phi^4)

Z r3_atail_bridge:
  alpha=1.0: (A_mu/A_e)^4 = 206.6  (matche PDG 206.77, diff 0.10%)
  alpha=2.0: (A_mu/A_e)^4 = 1221   (FAIL, factor 6 off)

Pytanie audytu: czy alpha=2 (TGP-canonical) daje spojny obraz N=3 + mass ratio?

WYNIK SPODZIEWANY:
  alpha=2 daje N=3 z bariery, ALE LAMIE mass formula A^4.
  alpha=1 daje correct mass formula A^4, ALE LAMIE TGP-canonical kinetic.
  R3 jest fundamentalnie sprzeczny z TGP-canonical alpha=2.

Autor: Audyt 2026-05-01
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
import math

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS {name}  {detail}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")
    return condition

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941   # R3 calibration of electron
R21_PDG = 206.7682   # m_mu/m_e PDG 2024
R31_PDG = 3477.23    # m_tau/m_e PDG 2024


# ================================================================
# UNIVERSAL ODE SOLVER
# Lagrangian: L = (1/2) g^{2 alpha} (g')^2 + g^3/3 - g^4/4
# (Note: R3 uses U(g) = g^3/3 - g^4/4, NOT g^3/3 - g^4/2 as stated
#  in some places. Universal conservation law uses U_cons = 2g^3/3 - g^4/2,
#  which is 2*L_potential_part. Both consistent.)
# ODE: g'' + (alpha/g)(g')^2 + ((d-1)/r)g' = (1-g) * g^{2-2*alpha}
# ================================================================

def solve_ode(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    singular = [False]
    g_min_val = [g0]

    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        if g < g_min_val[0]:
            g_min_val[0] = g
        rhs_val = (1.0 - g) * g**(2 - 2*alpha)
        if r < 1e-12:
            gpp = rhs_val / max(d, 1.0)
        else:
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]

    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol, singular[0], g_min_val[0]


def find_g0_crit(alpha, d=3, g0_lo=1.01, g0_hi=5.0, tol=1e-7):
    g_threshold = 0.005
    sol, sing, g_min = solve_ode(g0_hi, alpha, d)
    is_bad = sing or g_min < g_threshold or not sol.success
    if not is_bad:
        g0_hi = 10.0
        sol, sing, g_min = solve_ode(g0_hi, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success
        if not is_bad:
            return None

    for _ in range(70):
        g0_mid = (g0_lo + g0_hi) / 2
        sol, sing, g_min = solve_ode(g0_mid, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success
        if is_bad:
            g0_hi = g0_mid
        else:
            g0_lo = g0_mid
        if g0_hi - g0_lo < tol:
            break
    return (g0_lo + g0_hi) / 2


def extract_atail(r, g, r_min=80.0, r_max=250.0):
    """A_tail z far-field: (g-1)*r ~= B*cos(r) + C*sin(r) -> A = sqrt(B^2 + C^2)"""
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 10:
        return None
    u_f = (g[mask] - 1.0) * r_f

    def model(r, B, C):
        return B * np.cos(r) + C * np.sin(r)

    try:
        popt, _ = curve_fit(model, r_f, u_f, p0=[0.01, 0.01])
        return math.sqrt(popt[0]**2 + popt[1]**2)
    except:
        return None


def get_atail(g0, alpha, d=3):
    sol, sing, _ = solve_ode(g0, alpha, d)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])


def compute_K_action(r, g, alpha, d=3, r_max_int=200.0):
    """
    Calka kinetyczna K = int_0^R_max (1/2) g^{2 alpha} (g')^2 r^{d-1} dr
    """
    mask = r <= r_max_int
    r_f = r[mask]
    g_f = g[mask]
    gp = np.gradient(g_f, r_f)
    integrand = 0.5 * g_f**(2*alpha) * gp**2 * r_f**(d-1)
    return np.trapezoid(integrand, r_f)


# ================================================================
# RUN
# ================================================================

print("=" * 76)
print("  R3 AUDIT z TGP-canonical alpha=2 (z K(phi) = K_geo * phi^4)")
print("=" * 76)
print()
print("Backstory: sek08a definiuje akcje z K(phi) = K_geo * phi^4. W radialnej")
print("notacji to oznacza alpha = 2 w L = (1/2) g^{2*alpha} (g')^2. Wczesniejsze")
print("skrypty R3 uzywaly alpha = 1 jako 'substrat' bez wyjasnienia.")
print()

# ----------------------------------------------------------------
# SECTION 1: g0_crit dla alpha = 2 (TGP-canonical) vs alpha = 1 (R3 'substrat')
# ----------------------------------------------------------------
print("=" * 76)
print("  SEKCJA 1: g0_crit(alpha) — bariera topologiczna")
print("=" * 76)

g0_crit_alpha1 = find_g0_crit(alpha=1.0, d=3)
g0_crit_alpha2 = find_g0_crit(alpha=2.0, d=3)
g0_crit_alpha075 = find_g0_crit(alpha=0.75, d=3)

print(f"\n  g0_crit(alpha=0.75, d=3) = {g0_crit_alpha075:.4f}  (R3 'geometryczne')")
print(f"  g0_crit(alpha=1.00, d=3) = {g0_crit_alpha1:.4f}  (R3 'substrat')")
print(f"  g0_crit(alpha=2.00, d=3) = {g0_crit_alpha2:.4f}  (TGP-canonical phi^4)")

check("S1.1 alpha=1 g0_crit ~ 2.206 (R3 reproduces)",
      abs(g0_crit_alpha1 - 2.206) < 0.05,
      f"= {g0_crit_alpha1:.4f}")
check("S1.2 alpha=2 g0_crit < 2.206 (lower barrier)",
      g0_crit_alpha2 < g0_crit_alpha1,
      f"= {g0_crit_alpha2:.4f} < {g0_crit_alpha1:.4f}")

# ----------------------------------------------------------------
# SECTION 2: Czy g0^tau z Koide (1.729) miesci sie pod bariera dla alpha=2?
# ----------------------------------------------------------------
print()
print("=" * 76)
print("  SEKCJA 2: Czy phi-drabinka + g0^tau(Koide) miesci sie pod bariera?")
print("=" * 76)

g0_e = G0_E
g0_mu = g0_e * PHI            # phi-drabinka (R3 INPUT!)
g0_tau_koide = 1.72931        # Z Koide K=2/3 + g0_e + phi-drabinka (R3 INPUT!)

print(f"\n  g0^e   = {g0_e:.4f}  (R3 calibration, INPUT)")
print(f"  g0^mu  = {g0_mu:.4f}  (phi-drabinka g0^e * phi, INPUT)")
print(f"  g0^tau = {g0_tau_koide:.4f}  (Koide K=2/3 derived, INPUT)")
print()

# alpha = 1
print(f"  alpha=1.0 (R3 substrat): g0_crit = {g0_crit_alpha1:.4f}")
print(f"     g0^e   = {g0_e:.4f} < {g0_crit_alpha1:.4f}?  {g0_e < g0_crit_alpha1}")
print(f"     g0^mu  = {g0_mu:.4f} < {g0_crit_alpha1:.4f}?  {g0_mu < g0_crit_alpha1}")
print(f"     g0^tau = {g0_tau_koide:.4f} < {g0_crit_alpha1:.4f}?  {g0_tau_koide < g0_crit_alpha1}")
print(f"     g0^4 (=phi*g0^tau) = {g0_tau_koide*PHI:.4f} > {g0_crit_alpha1:.4f}?  {g0_tau_koide*PHI > g0_crit_alpha1}")
print()

# alpha = 2
print(f"  alpha=2.0 (TGP-canonical phi^4): g0_crit = {g0_crit_alpha2:.4f}")
print(f"     g0^e   = {g0_e:.4f} < {g0_crit_alpha2:.4f}?  {g0_e < g0_crit_alpha2}")
print(f"     g0^mu  = {g0_mu:.4f} < {g0_crit_alpha2:.4f}?  {g0_mu < g0_crit_alpha2}")
print(f"     g0^tau = {g0_tau_koide:.4f} < {g0_crit_alpha2:.4f}?  {g0_tau_koide < g0_crit_alpha2}")
print(f"     g0^4 (=phi*g0^tau) = {g0_tau_koide*PHI:.4f} > {g0_crit_alpha2:.4f}?  {g0_tau_koide*PHI > g0_crit_alpha2}")

margin_alpha1 = g0_crit_alpha1 - g0_tau_koide
margin_alpha2 = g0_crit_alpha2 - g0_tau_koide

print()
print(f"  Margin do bariery (g0_crit - g0^tau):")
print(f"    alpha=1.0: {margin_alpha1:+.4f}  (komfortowy)")
print(f"    alpha=2.0: {margin_alpha2:+.4f}  ({'wciaz dziala' if margin_alpha2 > 0 else 'TAU NAD BARIERA!'})")

check("S2.1 alpha=2: g0^tau pod bariera",
      g0_tau_koide < g0_crit_alpha2,
      f"margin = {margin_alpha2:+.4f}")
check("S2.2 alpha=2: g0^4 nad bariera (4. zakazana)",
      g0_tau_koide * PHI > g0_crit_alpha2,
      f"g0^4={g0_tau_koide*PHI:.4f} > {g0_crit_alpha2:.4f}")
check("S2.3 alpha=2: N=3 wciaz spojne z bariera",
      g0_tau_koide < g0_crit_alpha2 < g0_tau_koide * PHI,
      "")

# ----------------------------------------------------------------
# SECTION 3: A_tail extraction dla alpha = 2
# ----------------------------------------------------------------
print()
print("=" * 76)
print("  SEKCJA 3: A_tail i mass formula A^4 dla alpha = 2")
print("=" * 76)

print()
print("  alpha = 1.0 (R3 substrat) — REPRODUKCJA:")
A_e_a1 = get_atail(g0_e, alpha=1.0)
A_mu_a1 = get_atail(g0_mu, alpha=1.0)
A_tau_a1 = get_atail(g0_tau_koide, alpha=1.0)
if all(x is not None for x in [A_e_a1, A_mu_a1, A_tau_a1]):
    r21_a1 = (A_mu_a1 / A_e_a1)**4
    r31_a1 = (A_tau_a1 / A_e_a1)**4
    print(f"    A_e   = {A_e_a1:.5f}")
    print(f"    A_mu  = {A_mu_a1:.5f}")
    print(f"    A_tau = {A_tau_a1:.5f}")
    print(f"    (A_mu/A_e)^4  = {r21_a1:.2f}  vs PDG {R21_PDG} (diff {(r21_a1/R21_PDG-1)*100:+.2f}%)")
    print(f"    (A_tau/A_e)^4 = {r31_a1:.1f}  vs PDG {R31_PDG} (diff {(r31_a1/R31_PDG-1)*100:+.2f}%)")
    check("S3.1 alpha=1: m_mu/m_e ~ 207 z A^4",
          abs(r21_a1 / R21_PDG - 1) < 0.02,
          f"diff {(r21_a1/R21_PDG-1)*100:+.2f}%")

print()
print("  alpha = 2.0 (TGP-canonical phi^4):")
A_e_a2 = get_atail(g0_e, alpha=2.0)
A_mu_a2 = get_atail(g0_mu, alpha=2.0)
A_tau_a2 = get_atail(g0_tau_koide, alpha=2.0)
if all(x is not None for x in [A_e_a2, A_mu_a2, A_tau_a2]):
    r21_a2 = (A_mu_a2 / A_e_a2)**4
    r31_a2 = (A_tau_a2 / A_e_a2)**4
    print(f"    A_e   = {A_e_a2:.5f}")
    print(f"    A_mu  = {A_mu_a2:.5f}")
    print(f"    A_tau = {A_tau_a2:.5f}")
    print(f"    (A_mu/A_e)^4  = {r21_a2:.2f}  vs PDG {R21_PDG} (diff {(r21_a2/R21_PDG-1)*100:+.2f}%)")
    print(f"    (A_tau/A_e)^4 = {r31_a2:.1f}  vs PDG {R31_PDG} (diff {(r31_a2/R31_PDG-1)*100:+.2f}%)")
    check("S3.2 alpha=2: m_mu/m_e ~ 207 z A^4",
          abs(r21_a2 / R21_PDG - 1) < 0.10,
          f"diff {(r21_a2/R21_PDG-1)*100:+.2f}%")

# ----------------------------------------------------------------
# SECTION 4: Mechanism A^4 via K^2 (m = c * K^2)
# ----------------------------------------------------------------
print()
print("=" * 76)
print("  SEKCJA 4: Mass mechanism m = c * K^2 dla alpha = 1 vs 2")
print("=" * 76)

def compute_mass_K2(g0, alpha, d=3):
    sol, sing, _ = solve_ode(g0, alpha, d)
    if sing or not sol.success:
        return None
    K = compute_K_action(sol.t, sol.y[0], alpha, d)
    return K

print()
print("  alpha = 1.0 (R3 substrat):")
K_e_a1 = compute_mass_K2(g0_e, alpha=1.0)
K_mu_a1 = compute_mass_K2(g0_mu, alpha=1.0)
K_tau_a1 = compute_mass_K2(g0_tau_koide, alpha=1.0)
if all(x is not None for x in [K_e_a1, K_mu_a1, K_tau_a1]):
    rK21_a1 = (K_mu_a1 / K_e_a1)**2
    rK31_a1 = (K_tau_a1 / K_e_a1)**2
    print(f"    K_e   = {K_e_a1:.5f}")
    print(f"    K_mu  = {K_mu_a1:.5f}")
    print(f"    K_tau = {K_tau_a1:.5f}")
    print(f"    (K_mu/K_e)^2  = {rK21_a1:.2f}  vs PDG {R21_PDG} (diff {(rK21_a1/R21_PDG-1)*100:+.2f}%)")
    print(f"    (K_tau/K_e)^2 = {rK31_a1:.1f}  vs PDG {R31_PDG} (diff {(rK31_a1/R31_PDG-1)*100:+.2f}%)")

print()
print("  alpha = 2.0 (TGP-canonical):")
K_e_a2 = compute_mass_K2(g0_e, alpha=2.0)
K_mu_a2 = compute_mass_K2(g0_mu, alpha=2.0)
K_tau_a2 = compute_mass_K2(g0_tau_koide, alpha=2.0)
if all(x is not None for x in [K_e_a2, K_mu_a2, K_tau_a2]):
    rK21_a2 = (K_mu_a2 / K_e_a2)**2
    rK31_a2 = (K_tau_a2 / K_e_a2)**2
    print(f"    K_e   = {K_e_a2:.5f}")
    print(f"    K_mu  = {K_mu_a2:.5f}")
    print(f"    K_tau = {K_tau_a2:.5f}")
    print(f"    (K_mu/K_e)^2  = {rK21_a2:.2f}  vs PDG {R21_PDG} (diff {(rK21_a2/R21_PDG-1)*100:+.2f}%)")
    print(f"    (K_tau/K_e)^2 = {rK31_a2:.1f}  vs PDG {R31_PDG} (diff {(rK31_a2/R31_PDG-1)*100:+.2f}%)")
    check("S4.1 alpha=2: m_mu/m_e ~ 207 z K^2",
          abs(rK21_a2 / R21_PDG - 1) < 0.10,
          f"diff {(rK21_a2/R21_PDG-1)*100:+.2f}%")
    check("S4.2 alpha=2: m_tau/m_e ~ 3477 z K^2",
          abs(rK31_a2 / R31_PDG - 1) < 0.10,
          f"diff {(rK31_a2/R31_PDG-1)*100:+.2f}%")

# ----------------------------------------------------------------
# DIAGNOSE
# ----------------------------------------------------------------
print()
print("=" * 76)
print("  DIAGNOZA")
print("=" * 76)
print()
print("OBSERWACJE z numerycznego audytu:")
print()
print("  1. Bariera topologiczna g0_crit JEST UNIWERSALNA strukturalnie")
print("     (zalezy od alpha, ale MECHANIZM 'soliton ODE konczy sie")
print("     singularnoscia g_min=0' istnieje dla wszystkich alpha).")
print()
print("  2. Dla alpha=1 (R3 'substrat'): g0_crit ~ 2.21, mass formula A^4 OK.")
print("     Dla alpha=2 (TGP-canonical phi^4): g0_crit ~ 1.87, mass formula?")
print()
print("  3. Czy alpha=2 daje N=3? TAK strukturalnie:")
print("     g0^e=0.87 < g0^mu=1.41 < g0^tau=1.73 < g0_crit=1.87 < g0^4=2.80")
print("     ALE margin g0_crit - g0^tau = 0.14 (wezszy niz dla alpha=1: 0.48)")
print()
print("  4. Czy mass formula A^4 dziala dla alpha=2? NUMERYCZNIE SPRAWDZONE")
print("     wyzej w SEKCJI 3.")
print()
print("WNIOSEK: alpha = 2 (TGP-canonical) JEST kompatybilne z N=3 z bariery.")
print("Wczesniejsze R3 wnioski 'substrat alpha=1 daje N=2 deficit 3.1%' bylo")
print("BLEDNE z perspektywy TGP-canonical, bo R3 'substrat' to nie TGP-substrat.")
print()
print("Mass formula: jesli A^4 / K^2 daje OK ratio dla alpha=2, R3 jest")
print("fundamentalnie zgodne z TGP. Jesli LAMIE - R3 jest niespojne z TGP")
print("i wymaga albo zmiany TGP (sek08a K(phi) = phi^2 zamiast phi^4),")
print("albo zmiany R3 (porzucenie A^4 mass formula).")

print()
print("=" * 76)
print(f"  WYNIK: PASS={PASS}  FAIL={FAIL}")
print("=" * 76)
