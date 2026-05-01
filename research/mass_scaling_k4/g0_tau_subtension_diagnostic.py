#!/usr/bin/env python3
"""
g0_tau_subtension_diagnostic.py
================================

PURPOSE
-------
Diagnostyka sub-tension residualnej -0.085% w m_tau/m_e dla Phase 2 (alpha=2).

Hipotezy do przetestowania:
  H1: Numeryczna precyzja ODE solver (extract_atail curve_fit dla cos/sin tail)
  H2: K_PDG != 2/3 dokladnie (KILLED: K_PDG = 2/3 do 0.001%, w PDG error bars)
  H3: Empiryczne A^p(alpha) = A^3 != pelna formula c*A^2*g0^(e^2/2)
       => Section 4 r3_alpha2_full_closure.py uzywa SKROTU, nie pelnej formuly
  H4: Propagacja bledu z g0_mu = phi*g0_e (sama phi-drabinka jest empiryczna)

Test glowny: replikuje Section 2-4 z r3_alpha2_full_closure.py ALE z pelna formula:
  m_obs(g0, alpha) = c_M * A_tail^2 * g0^[e^2*(1-alpha/4)]
  dla alpha=2: m_obs = c_M * A_tail^2 * g0^(e^2/2)

Jesli residue spada <<0.085% -> sub-tension CLOSES jako artefakt A^3 approximation.
Jesli residue nadal ~0.085% -> H1 lub H4.

Autor: Audyt 2026-05-02 (g0_tau sub-tension closure attempt)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI

ALPHA = 2.0
P_EXP = 5 - ALPHA  # = 3.0

# Phase 2 universal exponent dla g0
E_FUND = math.sqrt(4 * math.pi / 137.036)  # e in Gaussian-like normalization?
# UWAGA: w PHASE2_n_alpha_derivation.md "e^2" to oznaczenie symboliczne kombinatorycznej
# stalej rzedu O(1), NIE alpha_em. Sprawdze bezposrednio z r3_alpha2_full_closure.txt.
# W skrypcie why_n3 e^2 ~ 0.0073... nie. Czytamy z dokumentu: n(alpha) = e^2*(1-alpha/4)
# z fitu numerycznego, gdzie "e" to constant ~ 1 fit-derived.
# Backsolve from Phase 2 fit: m_obs = c*A^p(alpha), p(alpha) = 5-alpha
# implicates n = e^2*(1-alpha/4) reproducuje p(alpha) iff... let me work out
#
# Z PHASE2: m_obs = c_M * A_tail^2 * g0^[e^2*(1-alpha/4)]
# Symbolicznie: e^2 = 2  =>  n(alpha) = 2*(1-alpha/4) = 2 - alpha/2
# Dla alpha=2: n = 1
# Wtedy m_obs = c*A^2*g0^1
# Empirical scaling test: jesli g0 ~ A^q (z phi-drabinki?), to m ~ A^(2+q)
# By bylo m ~ A^3, potrzeba q=1, czyli g0 ~ A
# To jest hipoteza testowalna numerycznie.

# Alternatywnie: e^2 to fit-parameter w Phase 2; sprawdzimy 2 i alpha_em^-1
# Use whatever Phase 2 paper says
E_SQUARED_HYP1 = 2.0      # symboliczne "e^2" ~ rank-2 invariant
E_SQUARED_HYP2 = 1.0      # naive
E_SQUARED_HYP3 = math.e**2  # natural log e^2 = 7.389

R21_PDG = 206.7682
R31_PDG = 3477.23
K_KOIDE = 2.0 / 3.0
G0_CRIT_ALPHA2 = 1.8744


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


print("=" * 78)
print("  g0_tau sub-tension diagnostic: A^3 approx vs pelna m_obs formula")
print("=" * 78)
print()
print(f"  alpha = {ALPHA}")
print(f"  g0_e = {G0_E}, g0_mu = phi*g0_e = {G0_MU:.6f}")
print()

# ----------------------------------------------------------------
# STEP 1: Compute A_e, A_mu via ODE
# ----------------------------------------------------------------
print("=" * 78)
print("  STEP 1: A_tail dla e i mu (alpha=2)")
print("=" * 78)
A_e = get_atail(G0_E, ALPHA)
A_mu = get_atail(G0_MU, ALPHA)
print(f"  A_e  = {A_e:.6f}")
print(f"  A_mu = {A_mu:.6f}")
print(f"  A_mu/A_e = {A_mu/A_e:.6f}")
print()

# ----------------------------------------------------------------
# STEP 2: Test rozne formuly mass dla mu/e ratio (kalibracja c_M)
# ----------------------------------------------------------------
print("=" * 78)
print("  STEP 2: Kalibracja exponent n dla g0 (z mu/e ratio)")
print("=" * 78)
print()
print("  Test: m_obs = c_M * A^2 * g0^n  =>  m_mu/m_e = (A_mu/A_e)^2 * (g0_mu/g0_e)^n")
print(f"  Cel: PDG m_mu/m_e = {R21_PDG}")
print()

# Solve for n: (A_mu/A_e)^2 * (g0_mu/g0_e)^n = R21_PDG
# (g0_mu/g0_e)^n = R21_PDG / (A_mu/A_e)^2
# n = log(R21_PDG / (A_mu/A_e)^2) / log(g0_mu/g0_e)
ratio_A_squared = (A_mu/A_e)**2
log_arg = R21_PDG / ratio_A_squared
n_fit = math.log(log_arg) / math.log(G0_MU/G0_E)
print(f"  (A_mu/A_e)^2 = {ratio_A_squared:.6f}")
print(f"  R21_PDG / (A_mu/A_e)^2 = {log_arg:.6f}")
print(f"  log(g0_mu/g0_e) = log({G0_MU/G0_E:.6f}) = {math.log(G0_MU/G0_E):.6f}")
print(f"  n (fit dla mu/e) = {n_fit:.6f}")
print()
print(f"  Phase 2 prescription: n = e^2*(1-alpha/4) = e^2 * 0.5 dla alpha=2")
print(f"  Therefore e^2 (z fitu) = {2*n_fit:.6f}")
print()

# Test rozne kandydaty
print("  Test kandydatow dla e^2:")
for label, e2_val in [("e^2 = 2 (rank-2 invariant)", 2.0),
                      ("e^2 = pi^2/g (Casimir-like)", math.pi**2/3),
                      ("alpha_em^-1 = 137.036", 137.036),
                      ("e (Euler) = 2.718", math.e),
                      ("e^2 (Euler^2) = 7.389", math.e**2),
                      (f"FIT (numerical) = {2*n_fit:.4f}", 2*n_fit)]:
    n_test = e2_val * (1 - ALPHA/4)
    pred = ratio_A_squared * (G0_MU/G0_E)**n_test
    print(f"    {label}:  n={n_test:.4f}  =>  m_mu/m_e = {pred:.4f}  (PDG {R21_PDG}, diff {(pred/R21_PDG-1)*100:+.4f}%)")
print()

# Use fit n
e2_FIT = 2 * n_fit
n_CANON = n_fit
print(f"  USE: n_canonical = {n_CANON:.6f} (from mu/e exact fit)")
print()


# ----------------------------------------------------------------
# STEP 3: Bisekt g0_tau wedlug DWOCH formul
# ----------------------------------------------------------------
print("=" * 78)
print("  STEP 3: Bisekt g0_tau dla DWOCH formul mass")
print("=" * 78)
print()

# Formula A: m ~ A^p(alpha) = A^3 (skrot empiryczny, Section 4 r3_alpha2_full_closure.py)
def mass_A3(g0, A_value):
    return A_value**P_EXP

# Formula B: m = c*A^2*g0^n_canon (pelna Phase 2)
def mass_full(g0, A_value):
    return A_value**2 * g0**n_CANON

# Calibracja c_M tak by m_e = 1 dla obu
c_A3 = 1.0 / A_e**P_EXP
c_full = 1.0 / (A_e**2 * G0_E**n_CANON)

print("  Calibracja c_M (by m_e = 1 unit dla obu formul):")
print(f"    Formula A (A^3):     c_M = {c_A3:.6e}")
print(f"    Formula B (A^2*g0^n): c_M = {c_full:.6e}")
print()

# Bisekt g0_tau wedlug Koide K=2/3 dla obu formul
def koide_K_from_masses(m1, m2, m3):
    sum_m = m1 + m2 + m3
    sum_sqrt = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return sum_m / sum_sqrt**2

def koide_eq_A3(g0_tau):
    A_tau = get_atail(g0_tau, ALPHA)
    if A_tau is None:
        return 1e10
    m_e = mass_A3(G0_E, A_e) * c_A3
    m_mu = mass_A3(G0_MU, A_mu) * c_A3
    m_tau = mass_A3(g0_tau, A_tau) * c_A3
    return koide_K_from_masses(m_e, m_mu, m_tau) - K_KOIDE

def koide_eq_full(g0_tau):
    A_tau = get_atail(g0_tau, ALPHA)
    if A_tau is None:
        return 1e10
    m_e = mass_full(G0_E, A_e) * c_full
    m_mu = mass_full(G0_MU, A_mu) * c_full
    m_tau = mass_full(g0_tau, A_tau) * c_full
    return koide_K_from_masses(m_e, m_mu, m_tau) - K_KOIDE

print("  Bisekt g0_tau (Koide K=2/3) dla Formula A (A^3 approximation):")
g0_tau_A3 = brentq(koide_eq_A3, 1.45, G0_CRIT_ALPHA2 - 0.01, xtol=1e-7)
A_tau_A3 = get_atail(g0_tau_A3, ALPHA)
m_e_A3 = mass_A3(G0_E, A_e) * c_A3
m_mu_A3 = mass_A3(G0_MU, A_mu) * c_A3
m_tau_A3 = mass_A3(g0_tau_A3, A_tau_A3) * c_A3
ratio_tau_e_A3 = m_tau_A3 / m_e_A3
ratio_mu_e_A3 = m_mu_A3 / m_e_A3
print(f"    g0_tau (A^3) = {g0_tau_A3:.6f}")
print(f"    A_tau        = {A_tau_A3:.6f}")
print(f"    m_mu/m_e     = {ratio_mu_e_A3:.4f}  (PDG {R21_PDG}, diff {(ratio_mu_e_A3/R21_PDG-1)*100:+.4f}%)")
print(f"    m_tau/m_e    = {ratio_tau_e_A3:.4f}  (PDG {R31_PDG}, diff {(ratio_tau_e_A3/R31_PDG-1)*100:+.4f}%)")
print()

print("  Bisekt g0_tau (Koide K=2/3) dla Formula B (pelna A^2*g0^n):")
g0_tau_full = brentq(koide_eq_full, 1.45, G0_CRIT_ALPHA2 - 0.01, xtol=1e-7)
A_tau_full = get_atail(g0_tau_full, ALPHA)
m_e_full = mass_full(G0_E, A_e) * c_full
m_mu_full = mass_full(G0_MU, A_mu) * c_full
m_tau_full = mass_full(g0_tau_full, A_tau_full) * c_full
ratio_tau_e_full = m_tau_full / m_e_full
ratio_mu_e_full = m_mu_full / m_e_full
print(f"    g0_tau (full) = {g0_tau_full:.6f}")
print(f"    A_tau         = {A_tau_full:.6f}")
print(f"    m_mu/m_e      = {ratio_mu_e_full:.4f}  (PDG {R21_PDG}, diff {(ratio_mu_e_full/R21_PDG-1)*100:+.4f}%)")
print(f"    m_tau/m_e     = {ratio_tau_e_full:.4f}  (PDG {R31_PDG}, diff {(ratio_tau_e_full/R31_PDG-1)*100:+.4f}%)")
print()


# ----------------------------------------------------------------
# STEP 4: Konkluzja
# ----------------------------------------------------------------
print("=" * 78)
print("  STEP 4: Konkluzja sub-tension")
print("=" * 78)
print()
print("  Phase 2 a posteriori test (Section 4 r3_alpha2_full_closure.py):")
print(f"    Formula A (A^3 skrot)     -> m_tau/m_e diff = {(ratio_tau_e_A3/R31_PDG-1)*100:+.4f}%")
print(f"    Formula B (full A^2*g0^n) -> m_tau/m_e diff = {(ratio_tau_e_full/R31_PDG-1)*100:+.4f}%")
print()
diff_A3 = abs(ratio_tau_e_A3/R31_PDG - 1)
diff_full = abs(ratio_tau_e_full/R31_PDG - 1)
print(f"  Improvement: {(1 - diff_full/diff_A3)*100:.1f}% redukcja residue")
print()
if diff_full < 0.0001:
    print("  WNIOSEK: Sub-tension CLOSES — pelna formula daje EXACT PDG.")
    print("           Residue -0.085% w r3_alpha2_full_closure.py to ARTEFAKT")
    print("           empirycznego A^3 approximation w Section 4.")
elif diff_full < diff_A3:
    print("  WNIOSEK: Sub-tension PARTIALLY closes — pelna formula daje lepsze")
    print("           dopasowanie ale nie EXACT. Pozostaje residual rzedu")
    print(f"           {diff_full*100:.4f}% — wymaga dalszych analiz (NLO corrections?).")
else:
    print("  WNIOSEK: Sub-tension NIE closes — pelna formula NIE poprawia residue.")
    print("           Inne hipotezy: H1 (numeryczna precyzja ODE), H4 (g0_mu kalibracja).")
print()
print("=" * 78)
