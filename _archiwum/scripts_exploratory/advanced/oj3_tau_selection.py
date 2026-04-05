#!/usr/bin/env python3
"""
oj3_tau_selection.py — O-J3: Zasada selekcji g0^tau
=====================================================
Tests two hypotheses for the tau selection principle:

  H1 (Koide): Given r_21 from phi-FP, the Koide formula K=2/3
              determines r_31 uniquely => g0^tau follows.
  H2 (Harmonic): g0^tau = 2*g0^e (with small correction epsilon_tau).

Both are verified against the substrate ODE numerical result.

Teoria Generowanej Przestrzeni — Mateusz Serafin, 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq, fsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

# --- Constants ---
PHI = (1.0 + np.sqrt(5)) / 2.0
R21_PDG = 206.768
R31_PDG = 3477.15
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86

# =====================================================================
# PART 1: Koide formula analysis (pure algebra, no ODE needed)
# =====================================================================
def koide_predict_r31(r21):
    """
    Given r_21 = m_mu/m_e, solve Koide formula K = 2/3 for r_31 = m_tau/m_e.

    K = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3

    Let x = sqrt(r_31) = sqrt(m_tau/m_e).
    Then with m_e = 1, m_mu = r_21, m_tau = r_31 = x^2:

    K = (1 + r_21 + x^2) / (1 + sqrt(r_21) + x)^2 = 2/3

    => 3(1 + r_21 + x^2) = 2(1 + sqrt(r_21) + x)^2
    """
    s = np.sqrt(r21)  # sqrt(m_mu/m_e)

    # 3(1 + r21 + x^2) = 2(1 + s + x)^2
    # 3 + 3*r21 + 3*x^2 = 2(1 + 2s + 2x + s^2 + 2sx + x^2)
    # 3 + 3*r21 + 3*x^2 = 2 + 4s + 4x + 2*s^2 + 4sx + 2*x^2
    # x^2 - (4 + 4s)*x + (3 + 3*r21 - 2 - 4s - 2*s^2) = 0
    # x^2 - (4 + 4s)*x + (1 + 3*r21 - 4s - 2*r21) = 0     [s^2 = r21]
    # x^2 - 4(1+s)*x + (1 + r21 - 4s) = 0

    a_coeff = 1.0
    b_coeff = -4.0 * (1.0 + s)
    c_coeff = 1.0 + r21 - 4.0 * s

    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff
    if discriminant < 0:
        return None, None

    x1 = (-b_coeff + np.sqrt(discriminant)) / (2.0 * a_coeff)
    x2 = (-b_coeff - np.sqrt(discriminant)) / (2.0 * a_coeff)

    # x = sqrt(r_31), so r_31 = x^2
    r31_1 = x1**2  # larger root
    r31_2 = x2**2  # smaller root

    return r31_1, r31_2


print("=" * 70)
print("  O-J3: TAU SELECTION PRINCIPLE")
print("=" * 70)

# ---- H1: Koide as predictor ----
print("\n" + "=" * 70)
print("  H1: KOIDE FORMULA AS TAU PREDICTOR")
print("=" * 70)

r31_k1, r31_k2 = koide_predict_r31(R21_PDG)

print(f"\n  Input: r_21 = {R21_PDG}")
print(f"  Koide equation: K = (1 + r_21 + r_31) / (1 + sqrt(r_21) + sqrt(r_31))^2 = 2/3")
print(f"\n  Solutions:")

for label, r31_k in [("Root 1 (physical)", r31_k1), ("Root 2", r31_k2)]:
    x = np.sqrt(r31_k)
    m_tau_k = M_E * r31_k
    # Verify Koide
    K_check = (1 + R21_PDG + r31_k) / (1 + np.sqrt(R21_PDG) + x)**2
    delta_r31 = (r31_k - R31_PDG) / R31_PDG
    print(f"\n    {label}:")
    print(f"      r_31     = {r31_k:.6f}")
    print(f"      m_tau    = {m_tau_k:.4f} MeV")
    print(f"      K_check  = {K_check:.10f} [should be {2/3:.10f}]")
    print(f"      delta K  = {K_check - 2/3:.2e}")
    print(f"      PDG r_31 = {R31_PDG}")
    print(f"      delta    = {delta_r31*100:.4f}%")
    print(f"      |m_tau - PDG| = {abs(m_tau_k - M_TAU):.4f} MeV ({abs(m_tau_k-M_TAU)/M_TAU*1e6:.0f} ppm)")

# Now: what r_31 does TGP's phi-FP r_21 predict via Koide?
print(f"\n  === KOIDE PREDICTION FROM phi-FP ===")
r21_tgp = 206.768  # exact phi-FP result
r31_koide, _ = koide_predict_r31(r21_tgp)
m_tau_koide = M_E * r31_koide
print(f"    r_21 (phi-FP) = {r21_tgp:.3f}")
print(f"    r_31 (Koide)  = {r31_koide:.4f}")
print(f"    m_tau (Koide)  = {m_tau_koide:.4f} MeV")
print(f"    PDG m_tau     = {M_TAU} MeV")
print(f"    delta         = {(m_tau_koide - M_TAU):.4f} MeV ({(m_tau_koide-M_TAU)/M_TAU*1e6:.0f} ppm)")

# Sensitivity: how does r_31(Koide) change with r_21?
print(f"\n  Sensitivity analysis:")
for dr in [-0.1, -0.01, 0, 0.01, 0.1, 1.0]:
    r21_test = R21_PDG + dr
    r31_test, _ = koide_predict_r31(r21_test)
    print(f"    r_21 = {r21_test:.3f} -> r_31 = {r31_test:.2f} "
          f"(delta_r31 = {(r31_test-R31_PDG)/R31_PDG*100:.3f}%)")

# ---- H2: Harmonic selection ----
print("\n" + "=" * 70)
print("  H2: HARMONIC SELECTION g0^tau = 2*g0^e")
print("=" * 70)

# Substrate ODE solver
def ode_sub(r, y):
    g, gp = y
    g = max(g, 1e-15)
    if r < 1e-10: return [gp, (1.0 - g) / 3.0]
    return [gp, (1.0 - g) - gp**2 / g - (2.0 / r) * gp]

def solve_A(g0, r_max=300.0):
    r_s = 1e-6; gpp0 = (1.0 - g0) / 3.0
    sol = solve_ivp(ode_sub, [r_s, r_max], [g0+0.5*gpp0*r_s**2, gpp0*r_s],
                    method='RK45', rtol=1e-10, atol=1e-12, max_step=1.0,
                    dense_output=True)
    if sol.status == -1 or sol.t[-1] < r_max*0.9: return np.nan
    r_u = np.linspace(100, min(sol.t[-1],250), 5000)
    g_u = sol.sol(r_u)[0]
    u = (g_u - 1.0) * r_u
    D = np.column_stack([np.cos(r_u), np.sin(r_u)])
    BC = np.linalg.lstsq(D, u, rcond=None)[0]
    return np.sqrt(BC[0]**2 + BC[1]**2)

# Compute A_tail for g0* and 2*g0*
g0_star = 0.86943  # from previous verification

# High-precision scan around g0*, phi*g0*, 2*g0*
print(f"\n  g0* = {g0_star:.5f}")
print(f"  phi*g0* = {PHI*g0_star:.5f}")
print(f"  2*g0* = {2*g0_star:.5f}")

Ae = solve_A(g0_star)
Amu = solve_A(PHI * g0_star)
A_2g0 = solve_A(2 * g0_star)

r21_check = (Amu/Ae)**4
r31_harmonic = (A_2g0/Ae)**4
m_tau_harmonic = M_E * r31_harmonic

print(f"\n  A_tail results:")
print(f"    A(g0*)      = {Ae:.8f}")
print(f"    A(phi*g0*)  = {Amu:.8f}")
print(f"    A(2*g0*)    = {A_2g0:.8f}")
print(f"\n  r_21 = (A_mu/A_e)^4 = {r21_check:.4f} [PDG: {R21_PDG}]")
print(f"  r_31(2g0*) = (A(2g0*)/A_e)^4 = {r31_harmonic:.2f} [PDG: {R31_PDG}]")
print(f"  m_tau(2g0*) = {m_tau_harmonic:.2f} MeV [PDG: {M_TAU}]")
print(f"  delta_r31 = {(r31_harmonic - R31_PDG)/R31_PDG*100:.3f}%")
print(f"  delta_m = {m_tau_harmonic - M_TAU:.2f} MeV ({(m_tau_harmonic-M_TAU)/M_TAU*1e6:.0f} ppm)")

# Compare with exact C4 value
# Find exact g0_tau from inverse problem
g0_scan = np.linspace(0.1, 3.0, 200)
A_scan = {}
for g0 in g0_scan:
    A = solve_A(g0)
    if not np.isnan(A) and A > 1e-10: A_scan[g0] = A

g0v = np.array(sorted(A_scan.keys()))
Av = np.array([A_scan[g] for g in g0v])
fA = interp1d(g0v, Av, kind='cubic', fill_value=np.nan, bounds_error=False)

def res_tau(g0):
    A = fA(g0)
    if np.isnan(A) or A <= 0: return 1e10
    return (A/Ae)**4 - R31_PDG

gt = np.linspace(PHI*g0_star+0.01, g0v[-1]-0.01, 3000)
rv = np.array([res_tau(g) for g in gt])
fm = np.isfinite(rv) & (np.abs(rv)<1e9)
sc = np.where(np.diff(np.sign(rv[fm])))[0]
g0_tau_exact = brentq(res_tau, gt[fm][sc[0]], gt[fm][sc[0]+1], xtol=1e-12)
A_tau_exact = fA(g0_tau_exact)

print(f"\n  === COMPARISON ===")
print(f"  g0^tau (C4 exact)  = {g0_tau_exact:.8f}")
print(f"  g0^tau (harmonic)  = {2*g0_star:.8f}")
print(f"  epsilon_tau = 2 - g0^tau/g0^e = {2 - g0_tau_exact/g0_star:.6f}")
print(f"  |epsilon_tau|      = {abs(2 - g0_tau_exact/g0_star)*100:.4f}%")

# Now: Koide with r_21 from phi-FP -> predicted r_31 -> predicted g0^tau
r31_from_koide, _ = koide_predict_r31(r21_check)
# Find g0 that gives r31_from_koide
def res_koide_g0(g0):
    A = fA(g0)
    if np.isnan(A) or A <= 0: return 1e10
    return (A/Ae)**4 - r31_from_koide

rv_k = np.array([res_koide_g0(g) for g in gt])
fm_k = np.isfinite(rv_k) & (np.abs(rv_k)<1e9)
sc_k = np.where(np.diff(np.sign(rv_k[fm_k])))[0]
if len(sc_k) > 0:
    g0_tau_koide = brentq(res_koide_g0, gt[fm_k][sc_k[0]], gt[fm_k][sc_k[0]+1])
    print(f"\n  g0^tau (Koide pred)= {g0_tau_koide:.8f}")
    print(f"  r_31 (Koide pred)  = {r31_from_koide:.4f}")
    print(f"  m_tau (Koide pred) = {M_E*r31_from_koide:.4f} MeV")

# ---- SUMMARY TABLE ----
print(f"\n" + "=" * 70)
print(f"  SUMMARY: TAU SELECTION MECHANISMS")
print(f"=" * 70)
print(f"\n  {'Method':<30s} {'g0^tau':>10s} {'r_31':>10s} {'m_tau[MeV]':>10s} {'delta_r31':>10s}")
print(f"  {'-'*72}")

# PDG
print(f"  {'PDG (reference)':<30s} {'—':>10s} {R31_PDG:10.2f} {M_TAU:10.2f} {'0%':>10s}")

# C4 (inverse, substrate ODE)
r31_c4 = (A_tau_exact/Ae)**4
print(f"  {'C4: inverse (substrate)':<30s} {g0_tau_exact:10.6f} {r31_c4:10.2f} "
      f"{M_E*r31_c4:10.2f} {(r31_c4/R31_PDG-1)*100:10.4f}%")

# Koide prediction
if len(sc_k) > 0:
    print(f"  {'Koide: r21->r31':<30s} {g0_tau_koide:10.6f} {r31_from_koide:10.2f} "
          f"{M_E*r31_from_koide:10.2f} {(r31_from_koide/R31_PDG-1)*100:10.4f}%")

# Harmonic 2*g0*
print(f"  {'Harmonic: 2*g0^e':<30s} {2*g0_star:10.6f} {r31_harmonic:10.2f} "
      f"{m_tau_harmonic:10.2f} {(r31_harmonic/R31_PDG-1)*100:10.4f}%")

# Koide check for each
for label, r31_val in [("C4", r31_c4), ("Koide", r31_from_koide), ("Harmonic", r31_harmonic)]:
    m_tau_v = M_E * r31_val
    K = (M_E + M_MU + m_tau_v) / (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(m_tau_v))**2
    print(f"    Koide({label}) = {K:.8f}, delta = {(K-2/3)*1e6/(2/3):.0f} ppm")

# ---- CHAIN OF REASONING ----
print(f"\n" + "=" * 70)
print(f"  CHAIN OF REASONING: phi-FP + Koide -> tau")
print(f"=" * 70)
print(f"""
  Step 1: phi-FP determines g0* (electron) from the substrate ODE.
          g0* = {g0_star:.5f}

  Step 2: Muon selection: g0^mu = phi * g0^e
          g0^mu = {PHI*g0_star:.5f}
          r_21 = (A_mu/A_e)^4 = {r21_check:.3f}

  Step 3: Koide formula K = 2/3 + r_21 known => solve for r_31:
          r_31(Koide) = {r31_from_koide:.4f}
          m_tau(Koide) = {M_E*r31_from_koide:.2f} MeV

  Step 4: Map r_31 back to g0^tau via A_tail:
          g0^tau(Koide) = {g0_tau_koide:.6f}

  Result: g0^tau/g0^e = {g0_tau_koide/g0_star:.6f} ~ 2
          epsilon_tau = 2 - g0^tau/g0^e = {2 - g0_tau_koide/g0_star:.6f}

  Consistency: Koide predicts m_tau to {abs(M_E*r31_from_koide - M_TAU)/M_TAU*1e6:.0f} ppm of PDG.
""")

# ---- PLOT ----
print("  Generating plot...")
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('O-J3: Tau Selection Principle\n'
             r'$\varphi$-FP + Koide $\rightarrow$ $g_0^\tau$',
             fontsize=13, fontweight='bold')

# Plot 1: A_tail with three leptons marked
ax = axes[0]
ax.semilogy(g0v, Av, 'b-', lw=1.5)
ax.axvline(g0_star, color='green', ls='--', lw=1.5, label=f'$g_0^e$={g0_star:.3f}')
ax.axvline(PHI*g0_star, color='orange', ls='--', lw=1.5, label=f'$g_0^\\mu$={PHI*g0_star:.3f}')
ax.axvline(g0_tau_exact, color='red', ls='--', lw=1.5, label=f'$g_0^\\tau$={g0_tau_exact:.3f}')
ax.axvline(2*g0_star, color='purple', ls=':', lw=1, alpha=0.7, label=f'$2g_0^e$={2*g0_star:.3f}')
ax.set_xlabel('$g_0$'); ax.set_ylabel('$A_{\\rm tail}$')
ax.set_title('$A_{\\rm tail}(g_0)$'); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Plot 2: Koide as function of r_31
ax = axes[1]
r31_range = np.linspace(2000, 5000, 1000)
K_vals = [(M_E + M_MU + M_E*r) / (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(M_E*r))**2
          for r in r31_range]
ax.plot(r31_range, K_vals, 'b-', lw=1.5)
ax.axhline(2/3, color='red', ls='--', label='$K = 2/3$')
ax.axvline(R31_PDG, color='green', ls='--', alpha=0.7, label=f'PDG $r_{{31}}$={R31_PDG:.0f}')
ax.axvline(r31_from_koide, color='purple', ls=':', label=f'Koide pred.={r31_from_koide:.0f}')
ax.set_xlabel('$r_{31}$'); ax.set_ylabel('$K(r_{31})$')
ax.set_title('Koide formula $K(r_{31})$ at fixed $r_{21}$')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
ax.set_ylim(0.64, 0.70)

# Plot 3: chain diagram
ax = axes[2]
ax.axis('off')
chain = [
    [r'$g_0^e = 0.869$ (substrate ODE)', ''],
    [r'$\downarrow$ $\varphi$-FP', ''],
    [r'$g_0^\mu = \varphi \cdot g_0^e = 1.407$', ''],
    [r'$\downarrow$ $A_{\rm tail}$', ''],
    [r'$r_{21} = (A_\mu/A_e)^4 = 206.77$', ''],
    [r'$\downarrow$ Koide $K = 2/3$', ''],
    [f'$r_{{31}} = {r31_from_koide:.2f}$', ''],
    [r'$\downarrow$ $A_{\rm tail}^{-1}$', ''],
    [f'$g_0^\\tau = {g0_tau_koide:.4f} \\approx 2 g_0^e$', ''],
    ['', ''],
    [f'$m_\\tau = {M_E*r31_from_koide:.1f}$ MeV', f'(PDG: {M_TAU})'],
]
table = ax.table(cellText=chain, loc='center', cellLoc='center',
                  colWidths=[0.6, 0.3])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.0, 1.8)
ax.set_title('Derivation chain', fontsize=11, pad=15)

plt.tight_layout()
outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'oj3_tau_selection.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"  Saved: {outpath}")
plt.close()

print("\n  DONE.\n")
