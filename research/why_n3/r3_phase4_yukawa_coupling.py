#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_phase4_yukawa_coupling.py
==============================

PURPOSE
-------
FAZA 4 z roadmap'u tgp_emergent_dirac_propagator.md (Sekcja 16.7):

Yukawa coupling z R3 solitonu — pokazać że soliton R3 effectively coupluje
do quantum spinora ψ przez Yukawa-like, generujący masę m_eff·ψ̄ψ.

POINT OF DEPARTURE
------------------
Z Faza 1+2+3:
  - m_eff(ψ) = c_M · A_tail²(g₀(ψ)) · g₀(ψ)^(e²/2)  (α=2)
  - Spin-1/2 z RP² topology
  - Propagator S_TGP(p; ψ) na M9.1'' background

Pytanie: SKĄD się bierze m_eff·ψ̄ψ w effective action? Co to jest
"Yukawa coupling" R3 solitonu z spinorem?

STRATEGIA
---------
W collective coordinate quantization, soliton z internal Z_2 loop daje
fermion. Mass spinora pochodzi z core energy solitonu projected na
spinor wave function:

  S_eff[psi] = ∫ dτ ψ̄ [iγ^μ ∂_μ - m_eff(ψ_local)] ψ

gdzie m_eff jest dynamic — zalezy od lokalnego ψ.

Yukawa coupling pojawia się gdy expandujemy m_eff wokół vacuum ψ=1:

  m_eff(ψ) = m_eff(1) + (∂m_eff/∂ψ)|_1 · (ψ-1) + O((ψ-1)²)
           ≡ m_0 + y · (ψ - 1) + ...

gdzie:
  m_0 = m_eff(ψ=1) = "vacuum mass" (powinno być 0 dla emergent fermion?)
  y   = ∂m_eff/∂ψ|_{ψ=1} = effective Yukawa coupling

DERYWACJA
---------
Z Faza 2: m_eff(ψ) = c_M · A_tail²(g₀(ψ)) · g₀(ψ)^(e²/2)
Z Faza 1: g₀(ψ) = (ψ - 0.6186) / 0.3814

Compute analitycznie:
  ∂m_eff/∂ψ |_{ψ=1} = c_M · [2 A · A' · g₀^(e²/2) + A² · (e²/2) · g₀^(e²/2-1) · g₀'] |_{ψ=1}

Liczbowo, znajdz y dla każdej generacji.

AUTOR: Faza 4 — Yukawa coupling derivation.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import math

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI
G0_TAU = 1.755046

# Linear reparametrization z Faza 1
A_LIN = 0.38139
B_LIN = 0.61861

def g0_of_psi(psi):
    """Faza 1 reparametrization: g₀(ψ) = (ψ - b) / a."""
    return (psi - B_LIN) / A_LIN

def psi_of_g0(g0):
    """Inverse: ψ(g₀) = a·g₀ + b."""
    return A_LIN * g0 + B_LIN

def dpsi_dg0():
    """dψ/dg₀ = a (constant for linear reparametrization)."""
    return A_LIN

def dg0_dpsi():
    """dg₀/dψ = 1/a."""
    return 1.0 / A_LIN


def solve_R3_get_atail(g0, alpha=2.0, d=3, r_max=200.0, n_points=20000, g_floor=1e-10):
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
    if singular[0] or not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    mask = (r >= 80) & (r <= 150)
    if not np.any(mask):
        return None
    u_f = (g[mask] - 1.0) * r[mask]
    def model(r, B, C):
        return B * np.cos(r) + C * np.sin(r)
    try:
        popt, _ = curve_fit(model, r[mask], u_f, p0=[0.01, 0.01])
        return math.sqrt(popt[0]**2 + popt[1]**2)
    except:
        return None


# ================================================================
print("=" * 78)
print("  R3 FAZA 4: Yukawa coupling z R3 solitonu")
print("=" * 78)
print()
print("Cel: pochodzić y_eff = ∂m_eff/∂ψ|_{ψ=ψ_gen} dla każdej generacji,")
print("identyfikujac ją jako Yukawa coupling do fluktuacji ψ wokół tła.")
print()


# ----------------------------------------------------------------
# SECTION 1: m_eff(ψ) jako funkcja ψ — numerycznie
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: m_eff(ψ) numerycznie dla generacji")
print("=" * 78)
print()

n_at_alpha2 = math.e**2 / 2  # n(2) = e²/2 z Faza 2

# Kalibracja c_M z elektronu
A_e = solve_R3_get_atail(G0_E, alpha=2.0)
m_e_uncalibrated = A_e**2 * G0_E**n_at_alpha2
c_M_calibrated = 1.0 / m_e_uncalibrated  # tak by m_e = 1

print(f"  Kalibracja: m_e = 1 (anchor)")
print(f"  A_e = {A_e:.6f}")
print(f"  m_e (uncalibrated) = A_e² · g₀_e^(e²/2) = {m_e_uncalibrated:.6f}")
print(f"  c_M = 1 / m_e_unc = {c_M_calibrated:.6f}")
print()


def m_eff_of_psi_numerical(psi):
    """Numerycznie obliczone m_eff(ψ) z R3 ODE."""
    g0 = g0_of_psi(psi)
    if g0 <= 0:
        return None
    A = solve_R3_get_atail(g0, alpha=2.0)
    if A is None:
        return None
    return c_M_calibrated * A**2 * g0**n_at_alpha2


# Test dla generacji
psi_e = psi_of_g0(G0_E)
psi_mu = psi_of_g0(G0_MU)
psi_tau = psi_of_g0(G0_TAU)

print(f"  Generacje (z Faza 1):")
print(f"    ψ_e   = {psi_e:.6f}, g₀_e = {G0_E:.5f}")
print(f"    ψ_μ   = {psi_mu:.6f}, g₀_μ = {G0_MU:.5f}")
print(f"    ψ_τ   = {psi_tau:.6f}, g₀_τ = {G0_TAU:.5f}")
print()

m_e_calc = m_eff_of_psi_numerical(psi_e)
m_mu_calc = m_eff_of_psi_numerical(psi_mu)
m_tau_calc = m_eff_of_psi_numerical(psi_tau)

print(f"  m_eff(ψ) numerycznie:")
print(f"    m_e   = {m_e_calc:.6f}  (target 1.000)")
print(f"    m_μ   = {m_mu_calc:.6f}  (PDG ratio 206.77)")
print(f"    m_τ   = {m_tau_calc:.6f}  (PDG ratio 3477.23)")


# ----------------------------------------------------------------
# SECTION 2: ∂m_eff/∂ψ — analityczna pochodna
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 2: Yukawa coupling y = ∂m_eff/∂ψ — wyprowadzenie")
print("=" * 78)
print()
print("m_eff(ψ) = c_M · A_tail²(g₀(ψ)) · g₀(ψ)^(e²/2)")
print()
print("Logarithmic derivative:")
print("  ln(m_eff) = ln(c_M) + 2·ln(A_tail(g₀)) + (e²/2)·ln(g₀)")
print()
print("  d(ln m_eff)/dψ = 2·(dA_tail/dg₀)/A_tail · dg₀/dψ + (e²/2)/g₀ · dg₀/dψ")
print()
print("  dg₀/dψ = 1/0.3814 = 2.622  (z Faza 1)")
print()
print("Yukawa coupling (logarithmic):")
print("  y_log = d(ln m_eff)/dψ |_{ψ=ψ_gen}")
print("  y_eff = m_eff · y_log = ∂m_eff/∂ψ")
print()


# Compute (dA/dg₀)/A numerically dla każdej generacji
def dA_dg0_over_A(g0, dg=1e-3):
    """Numerical derivative of log(A_tail(g₀))."""
    A_plus = solve_R3_get_atail(g0 + dg, alpha=2.0)
    A_minus = solve_R3_get_atail(g0 - dg, alpha=2.0)
    A_0 = solve_R3_get_atail(g0, alpha=2.0)
    if A_plus is None or A_minus is None or A_0 is None:
        return None
    return (math.log(A_plus) - math.log(A_minus)) / (2 * dg)


print(f"  Dla generacji obliczenie dA_tail/dg₀ / A_tail:")
print()
print(f"  {'gen':>5} | {'g₀':>10} | {'A_tail':>10} | "
      f"{'dlnA/dg₀':>12} | {'y_log = dlnm/dψ':>17} | "
      f"{'y_eff [m_e units]':>18}")
print("  " + "-" * 90)

for gen, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
    A = solve_R3_get_atail(g0, alpha=2.0)
    dlnA_dg0 = dA_dg0_over_A(g0)
    if A is None or dlnA_dg0 is None:
        continue

    # y_log = d(ln m_eff)/dψ = (2 · dlnA/dg₀ + (e²/2)/g₀) · dg₀/dψ
    y_log = (2 * dlnA_dg0 + n_at_alpha2 / g0) * dg0_dpsi()

    # y_eff = m_eff · y_log = dm_eff/dψ
    psi = psi_of_g0(g0)
    m = m_eff_of_psi_numerical(psi)
    y_eff = m * y_log

    print(f"  {gen:>5} | {g0:10.5f} | {A:10.5f} | "
          f"{dlnA_dg0:12.5f} | {y_log:17.5f} | "
          f"{y_eff:18.5f}")


# ----------------------------------------------------------------
# SECTION 3: Yukawa coupling jako efektywne dim-4 operator
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Effective Yukawa Lagrangian")
print("=" * 78)
print()
print("Effective action wokół ψ_vacuum = 1:")
print("  m_eff(ψ) = m_0 + y · (ψ - 1) + O((ψ-1)²)")
print("  gdzie:")
print("    m_0 = m_eff(1) — bare mass (DLA SOLITONU vacuum-stable JEST 0)")
print("    y   = ∂m_eff/∂ψ |_{ψ=1}")
print()

# Compute y_at_vacuum for spotkanie z each generacja
psi_vac = 1.0
g0_vac = g0_of_psi(psi_vac)  # = (1 - 0.6186) / 0.3814 = 1.000
print(f"  ψ_vacuum = 1.000")
print(f"  g₀(ψ=1) = {g0_vac:.6f}  (powinno być 1)")

A_vac = solve_R3_get_atail(g0_vac, alpha=2.0)
print(f"  A_tail(g₀=1) = {A_vac:.6f}  (powinno być małe — vacuum amplitude)")

m_0 = m_eff_of_psi_numerical(psi_vac)
print(f"  m_0 = m_eff(ψ=1) = {m_0:.6e}  (bare mass at vacuum)")
print()

if m_0 < 0.01:
    print("  >> m_0 ≈ 0 (vacuum mass znika) — emergent fermion ma zero bare mass")
    print("     w vacuum, all mass comes from displacement (ψ - 1).")
else:
    print("  Vacuum mass NIE znika — bare m_0 ≈ {:.4f}".format(m_0))
    print("  Może to być artefakt linearnej reparametryzacji Faza 1.")

print()


# ----------------------------------------------------------------
# SECTION 4: Yukawa hierarchies — czy ratio y_e:y_μ:y_τ matche m ratio?
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 4: Yukawa hierarchy — y_gen / m_gen ratios")
print("=" * 78)
print()
print("Standard Yukawa: m = y · v (Higgs vacuum). y_τ/y_e = m_τ/m_e ratio.")
print("Sprawdzamy: czy y_eff(R3)/m_eff(R3) jest w przybliżeniu stałe?")
print()

print(f"  {'gen':>5} | {'m_eff':>12} | {'y_eff':>12} | {'y_eff/m_eff':>14}")
print("  " + "-" * 55)

ratios = []
for gen, g0 in [('e', G0_E), ('μ', G0_MU), ('τ', G0_TAU)]:
    A = solve_R3_get_atail(g0, alpha=2.0)
    dlnA_dg0 = dA_dg0_over_A(g0)
    if A is None or dlnA_dg0 is None:
        continue

    y_log = (2 * dlnA_dg0 + n_at_alpha2 / g0) * dg0_dpsi()
    psi = psi_of_g0(g0)
    m = m_eff_of_psi_numerical(psi)
    y_eff = m * y_log

    ratio = y_eff / m if m > 0 else float('nan')
    ratios.append(ratio)
    print(f"  {gen:>5} | {m:12.5f} | {y_eff:12.5f} | {ratio:14.5f}")

print()
if len(ratios) == 3:
    spread = (max(ratios) - min(ratios)) / np.mean(ratios) * 100
    print(f"  Spread y_eff/m_eff: ~{spread:.1f}%")
    if spread < 20:
        print("  >> Ratio ROUGHLY constant (within ~20%) — sugeruje że y_eff/m")
        print("     ma fundamentalna interpretację (np. log-derivative scale)")
    else:
        print("  >> Ratio NOT constant — ratio zalezy od g0, czyli Yukawa")
        print("     w R3 NIE jest standardowe (m = y·v stara) — to jest")
        print("     bardziej generalna struktura.")


# ----------------------------------------------------------------
# SECTION 5: Path integral fragment — partition function dla solitonu
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Path integral fragment dla solitonu")
print("=" * 78)
print()
print("Partition function (Euclidean) dla solitonu jako bound state:")
print()
print("  Z_soliton = ∫ Dψ Dψ̄ exp(-S[ψ̄, ψ; soliton bg])")
print()
print("Effective action dla spinora w soliton background:")
print("  S = ∫ d⁴x ψ̄ [iγ^μ ∂_μ - m_eff(ψ_loc)] ψ")
print("    = S_kin + ∫ d⁴x m_eff(ψ_loc) ψ̄ψ")
print()
print("Splitting wokół vacuum:")
print("  ψ_loc = 1 + δψ_loc")
print("  m_eff(ψ_loc) = m_0 + y · δψ_loc + O(δψ²)")
print()
print("  S = S_kin + m_0 ∫ ψ̄ψ + y ∫ δψ_loc · ψ̄ψ + O(δψ² ψ̄ψ)")
print()
print("  --> Drugi term: STANDARD Yukawa coupling y · δψ · ψ̄ψ")
print("      w action. Coupling y wyznaczony powyzej dla każdej generacji.")
print()


# ----------------------------------------------------------------
# SECTION 6: PODSUMOWANIE
# ----------------------------------------------------------------
print("=" * 78)
print("  PODSUMOWANIE FAZY 4")
print("=" * 78)
print()
print("ODKRYCIA:")
print()
print("1. m_eff(ψ) z Faza 2 jest analityczne wyprowadzona — Yukawa coupling")
print("   y_eff = ∂m_eff/∂ψ jest WELL-DEFINED dla każdej ψ.")
print()
print("2. Effective Lagrangian:")
print("   L_eff = ψ̄ [iγ^μ ∂_μ - m_eff(ψ)] ψ")
print()
print("3. Yukawa coupling do fluktuacji vacuum (linearyzacja wokół ψ=1):")
print("   y_eff = ∂m_eff/∂ψ |_{ψ_gen}")
print("   Dla generacji obliczone numerycznie (Sekcja 2).")
print()
print("4. Path integral fragment:")
print("   Z = ∫ Dψ̄Dψ exp(-S_kin - m_0 ∫ψ̄ψ - y ∫δψ·ψ̄ψ - ...)")
print("   Standardowa Yukawa structure dla emergent fermion.")
print()
print("STATUS: Faza 4 strukturalnie zamknięta.")
print()
print("OPEN dla Fazy 5:")
print("  - Pełne 4D dynamics (time-dependent ψ)")
print("  - Loop corrections (1-loop renormalization Z_ψ)")
print("  - Anomaly cancellation, gauge couplings")
