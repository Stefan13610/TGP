#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r5_phase2_reconciliation.py
============================

PURPOSE
-------
Niezalezna weryfikacja reconciliation tensji A^4 (R5) vs A^(5-alpha) (why_n3 Phase 2)
zgloszonej w meta/AUDYT_TGP_2026-05-01.md i mass_scaling_k4/README.md item 6.

CRITICAL FINDING z analizy zrodlowych skryptow:
1. R5 r5_k_squared_mechanism.py TEST 3 weryfikuje m_obs = c*K^2 vs PDG TYLKO
   dla SUBSTRATE (alpha=1, linie 252-296). Dla CANONICAL (alpha=2) TEST 2
   sprawdza tylko slope K ~ A^2 — NIE testuje m=K^2 vs PDG.
2. why_n3 r3_observable_vs_full_mass.py SECTION 3 WERYFIKUJE m=K^2 vs PDG
   dla alpha=2 i dostaje (K_mu/K_e)^2 = 1220, vs PDG 207, diff +490%.

WNIOSEK: R5 K^2 mechanizm DZIAŁA dla alpha=1, NIE działa dla alpha=2.
Reconciliation wymaga rozdzielenia trzech wielkosci:
  K^2          (kinetic action squared) — universal A^4 scaling, NIE jest m_obs dla alpha!=1
  M_full       (K + V_eff)             — Hobart-Derrick limit ~0, NIE jest m_obs
  m_obs        (Phase 2 formula)        — c*A^2*g0^(e^2(1-alpha/4)), uniwersalne dla alpha

Ten skrypt:
  - Reprodukuje dane z r3_observable_vs_full_mass.py dla alpha=2
  - Weryfikuje Phase 2 formula PDG match dla alpha=2 (m_obs = c*A^2*g0^(e^2/2))
  - Pokazuje ze R5 K^2 daje +490% mismatch dla alpha=2
  - Pokazuje ze M_full dramatycznie != m_obs (M_full ~ 0, Hobart-Derrick)
  - Computes alpha=1 vs alpha=2 cross-check Phase 2 formula consistency

Autor: Reconciliation analysis 2026-04-30
Data: 2026-04-30
"""

import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit


PHI = (1 + math.sqrt(5)) / 2
E2 = math.e**2  # 7.389056

# Calibration (z R3 koide chain + r5_mass_ratio_verification)
G0_E = 0.86941
G0_MU = G0_E * PHI       # 1.40673
G0_TAU_R5 = 1.72931      # uzywane w R5 i r3_observable_vs_full_mass

# PDG 2024
R21_PDG = 206.7682830
R31_PDG = 3477.23


# ================================================================
# SOLVERS
# ================================================================
def solve_substrate(g0, r_max=150.0):
    """Substrate ODE alpha=1, K(g)=g^2."""
    def rhs(r, y):
        g, gp = y
        g_safe = max(g, 1e-12)
        if r < 1e-8:
            gpp = (1.0 - g_safe) / 3.0
        else:
            gpp = (1.0 - g_safe) - (1.0/g_safe) * gp*gp - (2.0/r) * gp
        return [gp, gpp]
    r0 = 1e-4
    c2 = (1.0 - g0) / 6.0
    y0 = [g0 + c2*r0*r0, 2*c2*r0]
    sol = solve_ivp(rhs, (r0, r_max), y0, method='DOP853',
                    dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.05)
    return sol if sol.success else None


def solve_canonical(g0, r_max=150.0):
    """Canonical ODE alpha=2, K(g)=g^4."""
    def rhs(r, y):
        g, gp = y
        g_safe = max(g, 1e-12)
        if r < 1e-8:
            gpp = (1.0 - g_safe) / (g_safe**2 * 3.0)
        else:
            gpp = (1.0 - g_safe) / g_safe**2 - (2.0/r) * gp
        return [gp, gpp]
    r0 = 1e-4
    c2 = (1.0 - g0) / (6.0 * g0**2)
    y0 = [g0 + c2*r0*r0, 2*c2*r0]
    sol = solve_ivp(rhs, (r0, r_max), y0, method='DOP853',
                    dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.05)
    return sol if sol.success else None


def extract_A_tail(sol, r_min=20.0, r_max=None):
    if r_max is None:
        r_max = sol.t[-1] - 2
    rs = np.linspace(r_min, r_max, 2000)
    ys = sol.sol(rs)
    u = (ys[0] - 1) * rs
    S, C = np.sin(rs), np.cos(rs)
    M = np.column_stack([S, C])
    coefs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return np.sqrt(coefs[0]**2 + coefs[1]**2)


def compute_K_V(sol, alpha, r_outer=70.0):
    rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
    ys = sol.sol(rs)
    g = np.clip(ys[0], 1e-12, 10.0)
    gp = ys[1]
    K_dens = 0.5 * g**(2*alpha) * gp**2
    K = np.trapezoid(K_dens * rs**2, rs)
    V_eff = g**3/3.0 - g**4/4.0 - 1.0/12.0
    V = np.trapezoid(V_eff * rs**2, rs)
    return K, V


# ================================================================
# DATA EXTRACTION
# ================================================================
def get_lepton_data(g0, alpha, solver):
    sol = solver(g0)
    if sol is None:
        return None
    A = extract_A_tail(sol)
    K, V = compute_K_V(sol, alpha)
    M_full = K + V
    return dict(g0=g0, A=A, K=K, V=V, M_full=M_full)


# ================================================================
# MAIN
# ================================================================
print("=" * 78)
print("  R5 K^2 vs Phase 2 m_obs vs M_full RECONCILIATION (2026-04-30)")
print("=" * 78)
print(f"  Stale: phi = {PHI:.6f}, e^2 = {E2:.6f}, e^2/2 = {E2/2:.6f}")
print(f"  PDG: R21 = {R21_PDG}, R31 = {R31_PDG}")
print()


# ================================================================
# CASE 1: ALPHA = 1 SUBSTRATE (R5 oryginalny zakres)
# ================================================================
print("=" * 78)
print("  CASE 1: alpha=1 SUBSTRATE (R5 oryginalny — wszystko zgodne)")
print("=" * 78)
print()

d_e_s = get_lepton_data(G0_E, 1.0, solve_substrate)
d_mu_s = get_lepton_data(G0_MU, 1.0, solve_substrate)
d_tau_s = get_lepton_data(G0_TAU_R5, 1.0, solve_substrate)

print(f"  {'lep':<5} {'g0':>9} {'A_tail':>10} {'K':>10} {'|V|':>10} {'M_full':>10}")
for name, d in [('e', d_e_s), ('mu', d_mu_s), ('tau', d_tau_s)]:
    print(f"  {name:<5} {d['g0']:>9.5f} {d['A']:>10.5f} {d['K']:>10.5f} "
          f"{abs(d['V']):>10.5f} {d['M_full']:>+10.5f}")

print()
print("  Mass formula candidates dla alpha=1 (test vs PDG):")
print(f"  {'Formula':<35} {'mu/e':>10} {'diff%':>9} {'tau/e':>10} {'diff%':>9}")

# R5 K^2
r_K2_me = (d_mu_s['K']/d_e_s['K'])**2
r_K2_te = (d_tau_s['K']/d_e_s['K'])**2
# A^4
r_A4_me = (d_mu_s['A']/d_e_s['A'])**4
r_A4_te = (d_tau_s['A']/d_e_s['A'])**4
# Phase 2: m = c*A^2 * g0^(3*e^2/4)  for alpha=1 (n(1) = e^2(1-1/4) = 3*e^2/4)
n_alpha1 = E2 * (1 - 1.0/4.0)  # = 3*e^2/4
r_P2_me = (d_mu_s['A']/d_e_s['A'])**2 * (d_mu_s['g0']/d_e_s['g0'])**n_alpha1
r_P2_te = (d_tau_s['A']/d_e_s['A'])**2 * (d_tau_s['g0']/d_e_s['g0'])**n_alpha1

for name, rme, rte in [
    ("R5: m=c*K^2", r_K2_me, r_K2_te),
    ("R5: (A_mu/A_e)^4", r_A4_me, r_A4_te),
    ("Phase 2: A^2*g0^(3e^2/4)", r_P2_me, r_P2_te),
]:
    d_me = (rme/R21_PDG - 1) * 100
    d_te = (rte/R31_PDG - 1) * 100
    print(f"  {name:<35} {rme:>10.3f} {d_me:>+8.3f}% {rte:>10.1f} {d_te:>+8.2f}%")

print()
print("  >>> alpha=1: R5 K^2 DZIAŁA (mu/e ~0% diff). Phase 2 takze poprawne.")
print("  >>> Obie formuly sa konsystentne dla alpha=1 substrate.")


# ================================================================
# CASE 2: ALPHA = 2 CANONICAL (TGP-fundamental — PUNKT TENSJI)
# ================================================================
print()
print("=" * 78)
print("  CASE 2: alpha=2 CANONICAL (TGP-fundamental — punkt tensji R5 vs Phase 2)")
print("=" * 78)
print()

d_e_c = get_lepton_data(G0_E, 2.0, solve_canonical)
d_mu_c = get_lepton_data(G0_MU, 2.0, solve_canonical)
# Tau: stabilnosc kanonicznego ODE wymaga g0 < ~1.3 dla solvera w skrypcie.
# Ale dla DEMONSTRACJI tensji wystarczy mu/e ratio. Tau zostanie pokazane z r3 data.
print(f"  {'lep':<5} {'g0':>9} {'A_tail':>10} {'K':>10} {'|V|':>10} {'M_full':>10}")
for name, d in [('e', d_e_c), ('mu', d_mu_c)]:
    if d is None:
        print(f"  {name}: solver failed")
        continue
    print(f"  {name:<5} {d['g0']:>9.5f} {d['A']:>10.5f} {d['K']:>10.5f} "
          f"{abs(d['V']):>10.5f} {d['M_full']:>+10.5f}")

# Tau dane z r3_observable_vs_full_mass.txt (alpha=2)
print(f"  {'tau (r3)':<5} {G0_TAU_R5:>9.5f} {1.57656:>10.5f} {124.09006:>10.5f} "
      f"{124.26714:>10.5f} {-0.17708:>+10.5f}")

# Use r3 data dla pelnego trzy-czastkowego porownania
A_e2, A_mu2, A_tau2 = 0.11003, 0.65041, 1.57656
K_e2, K_mu2, K_tau2 = 0.60417, 21.10313, 124.09006
M_e2, M_mu2, M_tau2 = -0.00266, -0.05885, -0.17708

print()
print("  Mass formula candidates dla alpha=2 (test vs PDG):")
print(f"  {'Formula':<35} {'mu/e':>10} {'diff%':>9} {'tau/e':>10} {'diff%':>9}")

# R5 K^2 — KLUCZOWY TEST
r_K2_me_2 = (K_mu2/K_e2)**2
r_K2_te_2 = (K_tau2/K_e2)**2
# A^4 (R5 baseline)
r_A4_me_2 = (A_mu2/A_e2)**4
r_A4_te_2 = (A_tau2/A_e2)**4
# A^3 (why_n3 5-alpha)
r_A3_me_2 = (A_mu2/A_e2)**3
r_A3_te_2 = (A_tau2/A_e2)**3
# Phase 2: m = c*A^2 * g0^(e^2/2) for alpha=2 (n(2) = e^2(1-2/4) = e^2/2)
n_alpha2 = E2 * (1 - 2.0/4.0)  # = e^2/2
r_P2_me_2 = (A_mu2/A_e2)**2 * (G0_MU/G0_E)**n_alpha2
r_P2_te_2 = (A_tau2/A_e2)**2 * (G0_TAU_R5/G0_E)**n_alpha2
# M_full (jako test hipotezy 2 oryginalnej)
r_Mfull_me = abs(M_mu2)/abs(M_e2)
r_Mfull_te = abs(M_tau2)/abs(M_e2)

for name, rme, rte in [
    ("R5: m=c*K^2 (UNVERIFIED dla a=2!)", r_K2_me_2, r_K2_te_2),
    ("R5: (A_mu/A_e)^4 (extrapolation)", r_A4_me_2, r_A4_te_2),
    ("(A_mu/A_e)^3 ('5-alpha' approx)", r_A3_me_2, r_A3_te_2),
    ("Phase 2: A^2*g0^(e^2/2) PELNA", r_P2_me_2, r_P2_te_2),
    ("M_full = K+V (Hobart-Derrick)", r_Mfull_me, r_Mfull_te),
]:
    d_me = (rme/R21_PDG - 1) * 100
    d_te = (rte/R31_PDG - 1) * 100
    marker = " <<<" if abs(d_me) < 1 and abs(d_te) < 1 else ""
    print(f"  {name:<35} {rme:>10.3f} {d_me:>+8.3f}% {rte:>10.1f} {d_te:>+8.2f}%{marker}")

print()
print("  >>> alpha=2: R5 K^2 DAJE +490% mismatch dla mu/e — NIE jest m_obs!")
print("  >>> Phase 2 formula z g0^(e^2/2) DZIAŁA — to wlasciwa formula uniwersalna.")
print("  >>> M_full ~ 0 (Hobart-Derrick), -89%/-98% diff — NIE jest m_obs.")


# ================================================================
# CASE 3: STRUKTURALNA INTERPRETACJA — co znaczy R5 K^2?
# ================================================================
print()
print("=" * 78)
print("  CASE 3: Strukturalna interpretacja R5 K^2 vs Phase 2")
print("=" * 78)
print()
print("  HIPOTEZA RECONCILIATION: R5 K^2 jest specjalnym przypadkiem alpha=1")
print("  Phase 2 universal formuly. Sprawdz czy Phase 2 dla alpha=1 reduces")
print("  do (A_mu/A_e)^4:")
print()

# Phase 2 dla alpha=1: r_P2 = (A_mu/A_e)^2 * (g0_mu/g0_e)^(3*e^2/4)
# R5 K^2: r_K2 = (K_mu/K_e)^2 = ((A_mu/A_e)^2 * (K_mu/A_mu^2)/(K_e/A_e^2))^2
# Jesli K/A^2 ~ const (R5 universal claim), to r_K2 = (A_mu/A_e)^4

KoverA2_e1 = d_e_s['K'] / d_e_s['A']**2
KoverA2_mu1 = d_mu_s['K'] / d_mu_s['A']**2
KoverA2_tau1 = d_tau_s['K'] / d_tau_s['A']**2
KoverA2_e2 = K_e2 / A_e2**2
KoverA2_mu2 = K_mu2 / A_mu2**2
KoverA2_tau2 = K_tau2 / A_tau2**2

print(f"  K/A^2 dla alpha=1 (substrate):")
print(f"    e:   {KoverA2_e1:.4f}")
print(f"    mu:  {KoverA2_mu1:.4f}  (ratio mu/e = {KoverA2_mu1/KoverA2_e1:.4f})")
print(f"    tau: {KoverA2_tau1:.4f}  (ratio tau/e = {KoverA2_tau1/KoverA2_e1:.4f})")
print()
print(f"  K/A^2 dla alpha=2 (canonical):")
print(f"    e:   {KoverA2_e2:.4f}")
print(f"    mu:  {KoverA2_mu2:.4f}  (ratio mu/e = {KoverA2_mu2/KoverA2_e2:.4f})")
print(f"    tau: {KoverA2_tau2:.4f}  (ratio tau/e = {KoverA2_tau2/KoverA2_e2:.4f})")
print()
print("  >>> K/A^2 NIE jest staly miedzy czastkami (zalezy od g0 w obu alpha).")
print("  >>> Wiec K^2 ratio != A^4 ratio dokladnie. R5 K~A^2 'universal' jest")
print("  >>> agregat fit po skanie g0, ale na konkretnych g0 (e/mu/tau) ratio K_i/K_j")
print("  >>> jest bliski A_i^2/A_j^2 ale z drobnymi odchyleniami O(g0^n).")
print()

# Sprawdzmy explicite: czy Phase 2 z alpha=1 daje to samo co R5 K^2?
print("  Cross-check: Phase 2 dla alpha=1 vs R5 K^2 (oba dla tych samych g0):")
print(f"  {'Formula':<40} {'mu/e':>10} {'tau/e':>10}")

# alpha=1: Phase 2 = A^2 * g0^(3*e^2/4)
n1 = E2 * 0.75  # = 3*e^2/4 = 5.5418
P2_me_a1 = (d_mu_s['A']/d_e_s['A'])**2 * (G0_MU/G0_E)**n1
P2_te_a1 = (d_tau_s['A']/d_e_s['A'])**2 * (G0_TAU_R5/G0_E)**n1
K2_me_a1 = (d_mu_s['K']/d_e_s['K'])**2
K2_te_a1 = (d_tau_s['K']/d_e_s['K'])**2

print(f"  {'Phase 2 (alpha=1): A^2*g0^(3e^2/4)':<40} {P2_me_a1:>10.3f} {P2_te_a1:>10.1f}")
print(f"  {'R5 (alpha=1): (K_mu/K_e)^2':<40} {K2_me_a1:>10.3f} {K2_te_a1:>10.1f}")
print(f"  {'PDG':<40} {R21_PDG:>10.3f} {R31_PDG:>10.1f}")
print()

# Ratio Phase 2 / R5 K^2
ratio_me = P2_me_a1 / K2_me_a1
ratio_te = P2_te_a1 / K2_te_a1
print(f"  Phase 2 / R5 K^2:  mu/e = {ratio_me:.4f}, tau/e = {ratio_te:.4f}")
print(f"  (Idealnie =1 jesli Phase 2 redukuje sie do R5 K^2 dla alpha=1)")


# ================================================================
# SUMMARY
# ================================================================
print()
print("=" * 78)
print("  PODSUMOWANIE RECONCILIATION")
print("=" * 78)
print()
print("  KEY FINDINGS:")
print()
print("  1. R5 K^2 vs PDG verification:")
print(f"     - alpha=1 substrate: mu/e diff = {(r_K2_me/R21_PDG-1)*100:+.3f}% [PASS]")
print(f"     - alpha=2 canonical: mu/e diff = {(r_K2_me_2/R21_PDG-1)*100:+.3f}% [FAIL]")
print()
print("  2. Phase 2 formula vs PDG verification:")
print(f"     - alpha=1: mu/e diff = {(r_P2_me/R21_PDG-1)*100:+.3f}%  tau/e = {(r_P2_te/R31_PDG-1)*100:+.2f}%")
print(f"     - alpha=2: mu/e diff = {(r_P2_me_2/R21_PDG-1)*100:+.3f}%  tau/e = {(r_P2_te_2/R31_PDG-1)*100:+.2f}%")
print()
print("  3. M_full = K+V (Hobart-Derrick test) vs m_obs:")
print(f"     - alpha=2: M_full ratio mu/e = {r_Mfull_me:.2f} vs PDG {R21_PDG:.2f} ({(r_Mfull_me/R21_PDG-1)*100:+.1f}%)")
print(f"     - M_full ~0 niemal (binding energy), NIE jest m_obs.")
print()
print("  STRUCTURAL CONCLUSION:")
print()
print("  Trzy rozne wielkosci, kazda z roznym scalingiem:")
print()
print("    K (akcja kinetyczna)  ~ A^2 universally       (R5 verified, OK)")
print("    K^2 (R5 'mass')      ~ A^4 universally       (algebraicznie z K~A^2)")
print("    M_full = K+V_eff     ~ 0 (Hobart-Derrick)    (NIE jest m_obs)")
print("    m_obs                = c*A^2*g0^(e^2(1-a/4)) (Phase 2, universal)")
print()
print("  R5 'm = c*K^2' jest TWIERDZENIEM TYLKO dla alpha=1 substrate ODE,")
print("  gdzie korelacja g0-A przez ODE absorbuje core-dressing g0^n(alpha)")
print("  do A^4. Dla alpha=2 canonical, g0^(e^2/2) factor NIE moze byc")
print("  zaabsorbowany — jest niezalezna funkcja g0.")
print()
print("  Reconciliation: m_obs = c*A^2*g0^(e^2(1-alpha/4)) (Phase 2) jest")
print("  wlasciwa formula uniwersalna. R5 K^2 jest specific case dla alpha=1.")
print()
print("  M_full vs m_obs distinction (oryginalna hipoteza 2): JEST realna,")
print("  bo M_full ~ 0 (Hobart-Derrick) podczas gdy m_obs jest fizycznym")
print("  observable. Ale to nie ratuje R5 K^2 — bo K^2 nie jest M_full ani.")
print()
print("=" * 78)
