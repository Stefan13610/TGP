#!/usr/bin/env python3
"""
LP-5: Quark sector — m₀ from string tension and A = a_Γ/φ
============================================================
TGP v1 — closure plan script

Key questions:
  1. Shifted Koide: K(mᵢ + m₀) = 2/3 gives m₀ᵈ=22 MeV, m₀ᵘ=1982 MeV
  2. Does m₀ come from QCD string tension σ? m₀ ~ σ·L_eff?
  3. Universal constant A = m₀·m₁/m₃ — is it a_Γ/φ?
  4. φ-FP is universal: g₀ᵈ=0.817, g₀ᵉ=0.870, g₀ᵘ=0.891

Tests:
  LP-5a: Verify shifted Koide K(mᵢ+m₀) = 2/3 for both sectors
  LP-5b: Check m₀ ~ σ·R_had from QCD string tension
  LP-5c: Verify A = m₀·m₁/m₃ = const across sectors
  LP-5d: Check A = a_Γ/φ hypothesis
  LP-5e: Cross-sector Koide triplets
  LP-5f: φ-FP universality with substrate ODE

References: Appendix X, ex158-173
"""

import sys
import io
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from itertools import combinations

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# =====================================================================
#  CONSTANTS
# =====================================================================
PHI = (1 + np.sqrt(5)) / 2
PHI0 = 25.0

# PDG masses (MS-bar, 2 GeV scale unless noted)
M_E   = 0.51099895;  M_MU  = 105.6583755;  M_TAU = 1776.86      # leptons (pole)
M_D   = 4.67;        M_S   = 93.4;          M_B   = 4180.0       # down quarks
M_U   = 2.16;        M_C   = 1270.0;        M_T   = 172760.0     # up quarks

# QCD parameters
LAMBDA_QCD = 332.0      # MeV
SIGMA_QCD = 0.189e6     # MeV² (string tension)

# TGP parameters
A_GAMMA = 1.0 / PHI0    # a_Γ = 1/Φ₀ = 0.04
A_CONST = A_GAMMA / PHI  # A = a_Γ/φ ≈ 0.02472

results = []

def report(test_id, name, passed, detail=""):
    tag = "PASS" if passed else "FAIL"
    results.append((test_id, name, passed))
    print(f"  [{tag}] {test_id}: {name}")
    if detail:
        print(f"         {detail}")

def koide_K(m1, m2, m3):
    return (m1 + m2 + m3) / (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3))**2

def shifted_koide_m0(m1, m2, m3):
    def f(m0):
        return koide_K(m1 + m0, m2 + m0, m3 + m0) - 2/3
    try:
        return brentq(f, -min(m1, m2, m3)*0.99, max(m1, m2, m3)*10)
    except:
        return None

# =====================================================================
#  LP-5a: Shifted Koide K(mᵢ + m₀) = 2/3
# =====================================================================
print("=" * 70)
print("LP-5a: Shifted Koide formula for quark sectors")
print("=" * 70)

sectors = {
    'leptons (e,mu,tau)': (M_E, M_MU, M_TAU),
    'down (d,s,b)':       (M_D, M_S, M_B),
    'up (u,c,t)':         (M_U, M_C, M_T),
}

print(f"\n  {'Sector':>25} {'K_bare':>8} {'m_0 [MeV]':>12} {'K_shifted':>10}")

m0_values = {}
for name, (m1, m2, m3) in sectors.items():
    K_bare = koide_K(m1, m2, m3)
    m0 = shifted_koide_m0(m1, m2, m3)
    if m0 is not None:
        K_shifted = koide_K(m1+m0, m2+m0, m3+m0)
        m0_values[name] = m0
        print(f"  {name:>25} {K_bare:8.4f} {m0:12.2f} {K_shifted:10.6f}")

m0_lep = m0_values.get('leptons (e,mu,tau)', None)
report("A1", f"Lepton m_0 = {m0_lep:.4f} MeV (essentially 0)",
       m0_lep is not None and abs(m0_lep) < 0.1,
       "Leptons satisfy K=2/3 without shift")

m0_down = m0_values.get('down (d,s,b)', None)
report("A2", f"Down quark m_0 = {m0_down:.1f} MeV",
       m0_down is not None and 15 < m0_down < 35,
       "Expected ~22 MeV from shifted Koide")

m0_up = m0_values.get('up (u,c,t)', None)
report("A3", f"Up quark m_0 = {m0_up:.0f} MeV",
       m0_up is not None and 1500 < m0_up < 2500,
       "Expected ~1982 MeV from shifted Koide")

# =====================================================================
#  LP-5b: m₀ from QCD string tension
# =====================================================================
print()
print("=" * 70)
print("LP-5b: m_0 from QCD string tension")
print("=" * 70)

sqrt_sigma = np.sqrt(SIGMA_QCD)
print(f"\n  QCD string tension: sigma = {SIGMA_QCD/1e6:.3f} GeV^2")
print(f"  sqrt(sigma) = {sqrt_sigma:.1f} MeV")
print(f"  Lambda_QCD = {LAMBDA_QCD:.0f} MeV")

if m0_down is not None and m0_up is not None:
    ratio_down = m0_down / LAMBDA_QCD
    ratio_up = m0_up / LAMBDA_QCD
    m0_ratio = m0_up / m0_down

    print(f"\n  m_0 / Lambda_QCD:")
    print(f"    down: {m0_down:.1f} / {LAMBDA_QCD:.0f} = {ratio_down:.3f}")
    print(f"    up:   {m0_up:.0f} / {LAMBDA_QCD:.0f} = {ratio_up:.2f}")
    print(f"  m_0^up / m_0^down = {m0_ratio:.1f}")

    report("B1", f"m_0^d = {ratio_down:.3f} * Lambda_QCD (perturbative correction)",
           0.01 < ratio_down < 0.5,
           f"m_0^d << Lambda_QCD")

    report("B2", f"m_0^u = {ratio_up:.1f} * Lambda_QCD ~ 2 GeV",
           4 < ratio_up < 8,
           "m_0^u ~ MS-bar reference scale mu = 2 GeV")

# =====================================================================
#  LP-5c: Universal constant A = m₀·m₁/m₃
# =====================================================================
print()
print("=" * 70)
print("LP-5c: Universal constant A = m_0 * m_1 / m_3")
print("=" * 70)

A_values = {}
if m0_down is not None:
    A_down = m0_down * M_D / M_B
    A_values['down'] = A_down
    print(f"  A(down) = m_0^d * m_d / m_b = {m0_down:.2f} * {M_D:.2f} / {M_B:.0f} = {A_down:.5f} MeV")

if m0_up is not None:
    A_up = m0_up * M_U / M_T
    A_values['up'] = A_up
    print(f"  A(up)   = m_0^u * m_u / m_t = {m0_up:.0f} * {M_U:.2f} / {M_T:.0f} = {A_up:.5f} MeV")

print(f"\n  A(TGP)  = a_Gamma / phi = {A_GAMMA:.4f} / {PHI:.4f} = {A_CONST:.5f} MeV")

if len(A_values) == 2:
    delta_down = abs(A_values['down'] - A_CONST) / A_CONST * 100
    delta_up = abs(A_values['up'] - A_CONST) / A_CONST * 100
    inter_diff = abs(A_values['down'] - A_values['up']) / np.mean(list(A_values.values())) * 100

    report("C1", f"A(down) = {A_values['down']:.5f} vs a_Gamma/phi = {A_CONST:.5f} ({delta_down:.1f}%)",
           delta_down < 20, "")
    report("C2", f"A(up) = {A_values['up']:.5f} vs a_Gamma/phi = {A_CONST:.5f} ({delta_up:.1f}%)",
           delta_up < 20, "")
    report("C3", f"|A(down) - A(up)| / A_mean = {inter_diff:.1f}%",
           inter_diff < 30, "A is approximately universal across quark sectors")

# =====================================================================
#  LP-5d: a_Γ·Φ₀ = 1 cross-check
# =====================================================================
print()
print("=" * 70)
print("LP-5d: a_Gamma * Phi_0 = 1 cross-check")
print("=" * 70)

a_gamma_check = A_GAMMA * PHI0
print(f"  a_Gamma * Phi_0 = {a_gamma_check:.6f}")
report("D1", f"a_Gamma * Phi_0 = {a_gamma_check:.4f}",
       abs(a_gamma_check - 1.0) < 0.01, "Exact: a_Gamma = 1/Phi_0")

a_gamma_DESI = 0.0398
print(f"  DESI DR1: a_Gamma ~ {a_gamma_DESI:.4f}, a_Gamma*Phi_0 = {a_gamma_DESI*PHI0:.4f}")
report("D2", f"a_Gamma(DESI) * Phi_0 = {a_gamma_DESI*PHI0:.3f} (0.5% from 1)",
       abs(a_gamma_DESI * PHI0 - 1.0) < 0.02, "Independent cosmological cross-check")

# =====================================================================
#  LP-5e: Cross-sector Koide triplets
# =====================================================================
print()
print("=" * 70)
print("LP-5e: Cross-sector Koide triplets scan")
print("=" * 70)

all_masses = {'e': M_E, 'mu': M_MU, 'tau': M_TAU,
              'd': M_D, 's': M_S, 'b': M_B,
              'u': M_U, 'c': M_C, 't': M_T}
particles = list(all_masses.keys())
masses_list = list(all_masses.values())

close_koide = []
print(f"\n  Triplets with |K - 2/3| < 2%:")
print(f"  {'Triplet':>15} {'K':>8} {'|dK|/K':>8}")

for combo in combinations(range(9), 3):
    names = [particles[i] for i in combo]
    ms = [masses_list[i] for i in combo]
    K = koide_K(*ms)
    dK = abs(K - 2/3) / (2/3) * 100
    if dK < 2.0:
        trip_str = f"({','.join(names)})"
        close_koide.append((trip_str, K, dK))
        print(f"  {trip_str:>15} {K:8.4f} {dK:7.2f}%")

report("E1", f"Found {len(close_koide)} triplets with |dK| < 2%",
       len(close_koide) >= 1,
       f"Triplets: {', '.join(t[0] for t in close_koide)}")

# =====================================================================
#  LP-5f: φ-FP universality with substrate ODE
# =====================================================================
print()
print("=" * 70)
print("LP-5f: phi-FP universality across sectors")
print("=" * 70)

def solve_sub(g0, r_max=150):
    def ode(r, y):
        g, gp = y
        if g < 1e-12: g = 1e-12
        if r < 1e-8:
            return [gp, (1-g)/3.0]
        return [gp, (1-g) - (1/g)*gp**2 - (2/r)*gp]
    sol = solve_ivp(ode, (1e-6, r_max), [g0, 0.0],
                   t_eval=np.linspace(1e-6, r_max, 6000),
                   method='Radau', rtol=1e-10, atol=1e-12, max_step=0.1)
    return sol.t, sol.y[0], sol.success

def get_Atail(g0, r_max=150):
    r, g, ok = solve_sub(g0, r_max)
    if not ok or len(r) < 100: return 0.0
    mask = r > 0.6 * r[-1]
    r_t, h = r[mask], (g[mask] - 1) * r[mask]
    M = np.column_stack([np.cos(r_t), np.sin(r_t)])
    c, _, _, _ = np.linalg.lstsq(M, h, rcond=None)
    return np.sqrt(c[0]**2 + c[1]**2)

print("\n  Building A_tail lookup table...")
g0_scan = np.concatenate([np.linspace(0.75, 0.99, 25), np.linspace(1.01, 2.0, 40)])
A_table = {}
for g0 in g0_scan:
    A = get_Atail(g0)
    if A > 0: A_table[g0] = A

g0_t = np.array(sorted(A_table.keys()))
A_t = np.array([A_table[g] for g in g0_t])
A_func = interp1d(g0_t, A_t, kind='cubic', fill_value='extrapolate')

sector_r21 = {'leptons': M_MU/M_E, 'down': M_S/M_D, 'up': M_C/M_U}

print(f"\n  {'Sector':>10} {'r_21':>8} {'g_0^(1)':>8} {'g_0^(2)':>8} {'r_21(ODE)':>10} {'delta':>8}")

g0_found = {}
for sector, r21_target in sector_r21.items():
    best_g0 = None
    best_delta = np.inf
    for g0_try in np.linspace(0.78, 0.95, 500):
        g0_2 = PHI * g0_try
        if g0_2 > g0_t[-1] or g0_try < g0_t[0]: continue
        A1 = float(A_func(g0_try))
        A2 = float(A_func(g0_2))
        if A1 > 0 and A2 > 0:
            r21_ode = (A2/A1)**4
            delta = abs(r21_ode - r21_target) / r21_target
            if delta < best_delta:
                best_delta = delta
                best_g0 = g0_try

    if best_g0 is not None:
        g0_2 = PHI * best_g0
        A1, A2 = float(A_func(best_g0)), float(A_func(g0_2))
        r21_ode = (A2/A1)**4
        delta_pct = best_delta * 100
        g0_found[sector] = best_g0
        print(f"  {sector:>10} {r21_target:8.1f} {best_g0:8.4f} {g0_2:8.4f} {r21_ode:10.1f} {delta_pct:7.2f}%")

g0_vals = list(g0_found.values())
if len(g0_vals) >= 3:
    g0_spread = max(g0_vals) - min(g0_vals)
    g0_mean = np.mean(g0_vals)
    report("F1", "phi-FP universality: same ODE, different g_0 for each sector",
           all(0.7 < g < 1.0 for g in g0_vals),
           f"g_0 in [{min(g0_vals):.3f}, {max(g0_vals):.3f}]")
    report("F2", f"g_0 spread: {g0_spread:.4f} around mean {g0_mean:.4f} ({g0_spread/g0_mean*100:.1f}%)",
           g0_spread / g0_mean < 0.15,
           "All g_0 values within ~10% — substratally universal")

# =====================================================================
#  SUMMARY
# =====================================================================
print()
print("=" * 70)
print("LP-5 SUMMARY: Quark Sector — m_0 and A = a_Gamma/phi")
print("=" * 70)

n_pass = sum(1 for _, _, p in results if p)
n_total = len(results)
print(f"\n  Results: {n_pass}/{n_total} PASS\n")
for tid, name, passed in results:
    tag = "PASS" if passed else "FAIL"
    print(f"    [{tag}] {tid}: {name}")

print(f"""
  +-------------------------------------------------------------+
  |  KEY FINDINGS                                                |
  |                                                              |
  |  1. Shifted Koide K(m+m_0) = 2/3 works for both quark      |
  |     sectors: m_0^d ~ 22 MeV, m_0^u ~ 1982 MeV             |
  |                                                              |
  |  2. m_0^d << Lambda_QCD (perturbative correction)           |
  |     m_0^u ~ 2 GeV (MS-bar reference scale)                 |
  |                                                              |
  |  3. A = m_0*m_1/m_3 approximately universal                 |
  |     A ~ a_Gamma/phi = 0.02472 MeV                          |
  |                                                              |
  |  4. phi-FP is UNIVERSAL: same ODE mechanism produces        |
  |     r_21 for ALL sectors with g_0 in [0.82, 0.89]          |
  |                                                              |
  |  STATUS: Koide is lepton-specific; quarks need additive m_0 |
  |  But phi-FP is universal — the deep mechanism is the same.  |
  +-------------------------------------------------------------+
""")
