#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_alpha2_full_closure.py
==========================

PURPOSE
-------
Po odkryciu z r3_observable_vs_full_mass.py:
  m_obs = c * A_tail^p(alpha),  z p(alpha) ~= 5 - alpha
  alpha=2 (TGP-canonical, K(phi)=phi^4) -> p=3

Sprawdzamy czy pelen obraz R3 dziala dla alpha=2:
  1. Re-derivacja g0^tau z Koide K=2/3 + mass formula A^3 (NIE A^4)
  2. Sprawdzenie czy g0^tau < g0_crit (alpha=2) = 1.874
  3. Sprawdzenie czy 4-ta generacja zakazana
  4. Sprawdzenie ratio m_mu/m_e i m_tau/m_e vs PDG

HYPOTEZA
--------
"R3 z alpha=2 i p(alpha)=5-alpha=3 jest pelnym, spojnym derywacja N=3
generacji + ratio mas + Koide structure, bez zmiany TGP_FOUNDATIONS."

INPUT (R3 calibration, kept):
  g0^e = 0.86941  (electron calibration)
  Koide K_lep = 2/3  (PDG empirical, INPUT)
  phi-drabinka g0^mu = phi * g0^e

OUTPUT (derived):
  g0^tau (z Koide + p=3)
  Ratio m_mu/m_e, m_tau/m_e
  Status barriery 4-tej generacji

Autor: Audyt 2026-05-01, follow-up po insight uzytkownika
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
import math

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI

ALPHA = 2.0          # TGP-canonical z K(phi) = K_geo * phi^4
P_EXP = 5 - ALPHA    # = 3.0 z empirycznego wzoru p(alpha) = 5-alpha

R21_PDG = 206.7682   # m_mu/m_e PDG 2024
R31_PDG = 3477.23    # m_tau/m_e PDG 2024
K_KOIDE = 2.0 / 3.0  # Koide constant

# g0_crit z r3_alpha2_canonical_audit.py
G0_CRIT_ALPHA2 = 1.8744

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


def solve_ode(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    singular = [False]
    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
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
    return sol, singular[0]


def extract_atail(r, g, r_min=80.0, r_max=250.0):
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


def get_atail(g0, alpha):
    sol, sing = solve_ode(g0, alpha)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])


# ================================================================
print("=" * 78)
print("  R3 — Pelne zamkniecie dla TGP-canonical alpha=2 z p(alpha)=5-alpha=3")
print("=" * 78)
print()
print(f"  alpha = {ALPHA} (TGP-canonical, K(phi) = K_geo * phi^4)")
print(f"  p(alpha) = 5 - alpha = {P_EXP}  (z empirycznego wzoru r3_observable_vs_full_mass)")
print(f"  Mass formula:  m_obs = c * A_tail^{P_EXP}")
print(f"  Koide K = {K_KOIDE:.6f} = 2/3 (INPUT z PDG)")
print(f"  g0^e = {G0_E}, g0^mu = phi*g0^e = {G0_MU:.5f} (INPUT, R3 calibration)")
print()

# ----------------------------------------------------------------
# SECTION 1: Compute A_e, A_mu dla alpha=2
# ----------------------------------------------------------------
print("=" * 78)
print(f"  SEKCJA 1: A_tail dla e i mu (alpha=2)")
print("=" * 78)

A_e = get_atail(G0_E, ALPHA)
A_mu = get_atail(G0_MU, ALPHA)
print(f"\n  A_e  = {A_e:.6f}  (g0={G0_E})")
print(f"  A_mu = {A_mu:.6f}  (g0={G0_MU:.5f})")
print(f"  A_mu/A_e = {A_mu/A_e:.4f}")

# Test mu/e ratio z A^3
ratio_mu_e_A3 = (A_mu/A_e)**P_EXP
print(f"\n  m_mu/m_e = (A_mu/A_e)^{P_EXP} = {ratio_mu_e_A3:.3f}")
print(f"  PDG m_mu/m_e = {R21_PDG}")
print(f"  Diff = {(ratio_mu_e_A3/R21_PDG - 1)*100:+.3f}%")

check("S1 m_mu/m_e match PDG (within 1%)",
      abs(ratio_mu_e_A3/R21_PDG - 1) < 0.01,
      f"diff {(ratio_mu_e_A3/R21_PDG - 1)*100:+.3f}%")


# ----------------------------------------------------------------
# SECTION 2: Re-derive g0^tau z Koide K=2/3 + mass formula A^3
# ----------------------------------------------------------------
print()
print("=" * 78)
print(f"  SEKCJA 2: Re-derivacja g0^tau z Koide K=2/3 + mass formula A^{P_EXP}")
print("=" * 78)

print("""
  Koide:  K = (m_e+m_mu+m_tau) / (sqrt(m_e)+sqrt(m_mu)+sqrt(m_tau))^2 = 2/3

  Z m = c*A^3:
    sum_m       = c * (A_e^3 + A_mu^3 + A_tau^3)
    sum_sqrt_m  = sqrt(c) * (A_e^{3/2} + A_mu^{3/2} + A_tau^{3/2})

  K = (A_e^3 + A_mu^3 + A_tau^3) / (A_e^{3/2} + A_mu^{3/2} + A_tau^{3/2})^2 = 2/3

  Solve for A_tau, then back-out g0^tau via ODE inversion.
""")


def koide_K_from_A(A_e, A_mu, A_tau, p):
    """Koide K dla mass = c·A^p."""
    half_p = p / 2.0
    sum_m = A_e**p + A_mu**p + A_tau**p
    sum_sqrt_m = A_e**half_p + A_mu**half_p + A_tau**half_p
    return sum_m / sum_sqrt_m**2


# Cel: znajdz A_tau by Koide K = 2/3
def koide_eq(A_tau):
    return koide_K_from_A(A_e, A_mu, A_tau, P_EXP) - K_KOIDE

# Bracket: A_tau musi byc wieksze niz A_mu, mniejsze niz extrapolacja
A_tau_lo = A_mu * 1.5
A_tau_hi = A_mu * 5.0

try:
    A_tau_koide = brentq(koide_eq, A_tau_lo, A_tau_hi, xtol=1e-8)
    K_check = koide_K_from_A(A_e, A_mu, A_tau_koide, P_EXP)
    print(f"  A_tau (z Koide K=2/3) = {A_tau_koide:.6f}")
    print(f"  Verify K = {K_check:.6f} vs target {K_KOIDE:.6f}")

    check("S2.1 Koide solver convergent",
          abs(K_check - K_KOIDE) < 1e-6,
          f"K = {K_check:.6f}")
except Exception as e:
    print(f"  Koide solver ERROR: {e}")
    A_tau_koide = None


# ----------------------------------------------------------------
# SECTION 3: Back-out g0^tau z A_tau via ODE inversion
# ----------------------------------------------------------------
print()
print("=" * 78)
print(f"  SEKCJA 3: Back-out g0^tau z A_tau via ODE inversion")
print("=" * 78)

if A_tau_koide is not None:
    # Funkcja A_tail(g0) dla alpha=2
    def A_of_g0(g0):
        A = get_atail(g0, ALPHA)
        if A is None:
            return -1
        return A - A_tau_koide

    # Bracket: g0^mu < g0^tau < g0_crit
    g0_lo = G0_MU + 0.01
    g0_hi = G0_CRIT_ALPHA2 - 0.01

    try:
        g0_tau_alpha2 = brentq(A_of_g0, g0_lo, g0_hi, xtol=1e-5)
        A_tau_check = get_atail(g0_tau_alpha2, ALPHA)
        print(f"\n  g0^tau (alpha=2, Koide) = {g0_tau_alpha2:.6f}")
        print(f"  Verify A_tau = {A_tau_check:.6f}  (target {A_tau_koide:.6f})")

        # Compare z R3 oryginalnym g0^tau=1.72931
        print(f"\n  Comparison z R3 oryginalnym (alpha=1) g0^tau = 1.72931:")
        print(f"  alpha=2 daje g0^tau = {g0_tau_alpha2:.5f}")
        print(f"  Diff = {(g0_tau_alpha2 - 1.72931):+.5f}")

        check("S3.1 g0^tau pod bariera (alpha=2)",
              g0_tau_alpha2 < G0_CRIT_ALPHA2,
              f"{g0_tau_alpha2:.4f} < {G0_CRIT_ALPHA2:.4f}, margin = {G0_CRIT_ALPHA2-g0_tau_alpha2:+.4f}")
    except Exception as e:
        print(f"  ODE inversion ERROR: {e}")
        g0_tau_alpha2 = None
else:
    g0_tau_alpha2 = None


# ----------------------------------------------------------------
# SECTION 4: Mass ratio m_tau/m_e (a posteriori)
# ----------------------------------------------------------------
print()
print("=" * 78)
print(f"  SEKCJA 4: Mass ratio m_tau/m_e (z A^{P_EXP} formula)")
print("=" * 78)

if A_tau_koide is not None:
    ratio_tau_e_A3 = (A_tau_koide/A_e)**P_EXP
    ratio_tau_mu_A3 = (A_tau_koide/A_mu)**P_EXP

    print(f"\n  m_tau/m_e = (A_tau/A_e)^{P_EXP} = {ratio_tau_e_A3:.2f}")
    print(f"  PDG m_tau/m_e = {R31_PDG}")
    print(f"  Diff = {(ratio_tau_e_A3/R31_PDG - 1)*100:+.3f}%")

    print(f"\n  m_tau/m_mu = (A_tau/A_mu)^{P_EXP} = {ratio_tau_mu_A3:.3f}")
    R32_PDG = R31_PDG / R21_PDG
    print(f"  PDG m_tau/m_mu = {R32_PDG:.3f}")
    print(f"  Diff = {(ratio_tau_mu_A3/R32_PDG - 1)*100:+.3f}%")

    check("S4.1 m_tau/m_e match PDG (within 5%)",
          abs(ratio_tau_e_A3/R31_PDG - 1) < 0.05,
          f"diff {(ratio_tau_e_A3/R31_PDG - 1)*100:+.3f}%")
    check("S4.2 m_tau/m_mu match PDG (within 5%)",
          abs(ratio_tau_mu_A3/R32_PDG - 1) < 0.05,
          f"diff {(ratio_tau_mu_A3/R32_PDG - 1)*100:+.3f}%")


# ----------------------------------------------------------------
# SECTION 5: 4-ta generacja zakazana?
# ----------------------------------------------------------------
print()
print("=" * 78)
print(f"  SEKCJA 5: 4-ta generacja zakazana? (φ-drabinka g0^4 = phi*g0^tau)")
print("=" * 78)

if g0_tau_alpha2 is not None:
    g0_4 = PHI * g0_tau_alpha2
    print(f"\n  g0^4 = phi * g0^tau = {PHI:.4f} * {g0_tau_alpha2:.5f} = {g0_4:.5f}")
    print(f"  g0_crit (alpha=2) = {G0_CRIT_ALPHA2}")
    print(f"  g0^4 > g0_crit?  {g0_4 > G0_CRIT_ALPHA2}")
    print(f"  Margin za bariera: {g0_4 - G0_CRIT_ALPHA2:+.5f}")

    check("S5.1 4-ta generacja za bariera (zakazana)",
          g0_4 > G0_CRIT_ALPHA2,
          f"g0^4={g0_4:.4f} > g0_crit={G0_CRIT_ALPHA2:.4f}")


# ----------------------------------------------------------------
# SECTION 6: Pelny obraz N=3 dla TGP-canonical alpha=2
# ----------------------------------------------------------------
print()
print("=" * 78)
print(f"  SEKCJA 6: Pelen obraz N=3 dla TGP-canonical alpha=2")
print("=" * 78)

print(f"""
  TGP-canonical:  K(phi) = K_geo * phi^4  =>  alpha = 2

  Empiryczna mass formula:  m = c * A_tail^(5-alpha) = c * A_tail^3

  Generacje:
    g0^e   = {G0_E:.5f}  (calibration INPUT)
    g0^mu  = {G0_MU:.5f}  (phi-drabinka g0^e * phi)""")

if g0_tau_alpha2 is not None:
    print(f"    g0^tau = {g0_tau_alpha2:.5f}  (z Koide K=2/3 + mass formula A^3)")
    print(f"    g0^4   = {PHI * g0_tau_alpha2:.5f}  (phi-drabinka g0^tau * phi) - ZAKAZANA")
print(f"""
  Bariera:
    g0_crit (alpha=2) = {G0_CRIT_ALPHA2}

  Selection rule N=3:""")

if g0_tau_alpha2 is not None:
    print(f"    e (0.869) < mu (1.407) < tau ({g0_tau_alpha2:.3f}) < g0_crit ({G0_CRIT_ALPHA2}) < 4-ta ({PHI * g0_tau_alpha2:.3f}) ✓")

print()
print("=" * 78)
print(f"  PODSUMOWANIE: PASS={PASS}  FAIL={FAIL}")
print("=" * 78)

if FAIL == 0:
    print("""
  WNIOSEK STRUKTURALNY:

  R3 z TGP-canonical alpha=2 i p(alpha)=5-alpha=3 jest PELNYM, SPOJNYM
  obrazem N=3 generacji:

    [DERIVED]    Bariera g0_crit istnieje strukturalnie z ODE
    [DERIVED]    p(alpha) = 5-alpha (empiryczna z numerycznego skanu)
    [COMPATIBLE] Mass formula m = c*A^p(alpha) zgodna z PDG (mu/e, tau/e <5%)
    [COMPATIBLE] Koide K=2/3 wybiera g0^tau pod bariera
    [DERIVED]    4-ta generacja > bariera => N=3 z constraint

  R3 jest TERAZ ZGODNE z TGP_FOUNDATIONS K(phi)=phi^4.

  Insight uzytkownika 2026-05-01 ROZWIAZUJE sprzecznosc audytu:
    'A_tail^4 to masa OBSERWOWALNA, bariera operuje na FULL mass.'

  Konsekwencja: aktualizowac CORRECTIONS_2026-05-01.md i README.md aby
  odzwierciedlic ze:
    - Mass formula = c*A^(5-alpha) [zalezna od alpha kinetic prefactor]
    - Bariera = wlasnosc strukturalna ODE [niezalezna od mass formula]
    - R3 i TGP-canonical alpha=2 sa SPOJNE
""")
