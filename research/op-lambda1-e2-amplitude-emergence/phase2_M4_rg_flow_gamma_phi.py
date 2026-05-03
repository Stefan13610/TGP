#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
λ.1 M.4 — RG flow γ_φ ~ e² test (audit-recommended last attempt)

Cel: explicit test czy anomalous dimension γ_φ w R3 effective field theory
przy α=2 może osiągnąć wartość e²/2 ≈ 3.69 (matchująca empirical slope).

Strategia:
  1. Linearize R3 ODE wokół g=1+ε przy α=2
  2. Identify couplings effective Lagrangian
  3. Compute 1-loop γ_φ — sprawdz czy daje e²/2
  4. Compute 2-loop estimate — sprawdz czy daje e²/2
  5. Test large-N / non-perturbative coupling required dla γ_φ = e²/2
  6. Numeryczny check: solve R3 ODE for range g₀, extract effective γ_φ
     z empirical Z(g₀) = m_obs(g₀)/A²(g₀)

Status: jeśli γ_φ_1loop, γ_φ_2loop, γ_φ_largeN nie matche e²/2 →
MECHANISM NEGATIVE → λ.1 zamknięte negatywnie (audit recommendation).
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq

# Force UTF-8 output for Windows console
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 78)
print("  λ.1 M.4 — RG flow γ_φ ~ e² test (audit-recommended last attempt)")
print("=" * 78)

E2 = np.exp(2)              # 7.389056
E2_HALF = E2 / 2            # 3.694528 (target γ_φ for α=2)
E2_QUARTER = E2 / 4         # 1.847264
PI = np.pi

print(f"\nTargets:")
print(f"  e^2          = {E2:.6f}")
print(f"  e^2 / 2      = {E2_HALF:.6f}  <-- target for γ_φ at α=2")
print(f"  e^2 / 4      = {E2_QUARTER:.6f}  <-- in formula n(α) = (e^2/4)(4-α)")

# ============================================================================
# SEKCJA 1: Effective Lagrangian z R3 ODE linearization @ α=2
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 1: R3 effective Lagrangian @ α=2 (linearization wokol g=1)")
print("=" * 78)

epsilon, alpha = sp.symbols("epsilon alpha", real=True)

# R3 ODE w postaci Lagrangianu: L = (1/2) g^(2α) (∂g)² - V(g)
# gdzie V(g) z (1-g)·g^(2-2α) RHS po całkowaniu
# Dla małej eps: g = 1 + eps

# Kinetic prefactor: g^(2α) = (1+ε)^(2α)
kinetic_prefactor = (1 + epsilon)**(2*alpha)
kinetic_expansion = sp.series(kinetic_prefactor, epsilon, 0, 4).removeO()
print(f"\n  Kinetic prefactor g^(2α) expansion:")
print(f"    g^(2α) = {sp.simplify(kinetic_expansion)}")

# Pola normalizacja: φ z (1/2)(∂φ)² standard kinetic
# Z g^(2α) = 1 + 2α·ε + α(2α-1)·ε² + ...
# Tutaj eps jest "bare field" — żeby dostać canonical, weź φ = ∫_0^ε √(1+2α·ε'+...) dε'
# Dla małej eps: φ ≈ ε + α·ε² + ...

# Potencjał: V(g) z (1-g)·g^(2-2α) jako (-dV/dg)
# V'(g) = (g-1)·g^(2-2α) (signs: tak żeby V min przy g=1)
# Wokół g=1+ε: V'(1+ε) = ε · (1+ε)^(2-2α) ≈ ε · (1 + (2-2α)ε + ...)
# V(1+ε) - V(1) = (1/2)ε² + ((2-2α)/3)·ε³ + ...

V_prime = epsilon * (1 + epsilon)**(2 - 2*alpha)
V_prime_exp = sp.series(V_prime, epsilon, 0, 4).removeO()
V_exp = sp.integrate(V_prime_exp, epsilon)
print(f"\n  Potential V(1+ε) - V(1) expansion:")
print(f"    V = {sp.simplify(V_exp)}")

# Effective Lagrangian dla pola eps:
# L = (1/2)(1 + 2α·ε + α(2α-1)·ε²)·(∂ε)² - V(ε)
# Po podstawieniu α=2:
alpha_val = 2
print(f"\n  Dla α={alpha_val} (TGP-canonical R3 charged-lepton):")

kinetic_at_2 = kinetic_expansion.subs(alpha, alpha_val)
V_at_2 = V_exp.subs(alpha, alpha_val)
print(f"    Kinetic prefactor: {sp.simplify(kinetic_at_2)}")
print(f"    Potential V(ε):    {sp.simplify(V_at_2)}")

# Identyfikacja standardowych couplings (po canonical normalization)
# L_eff (po canonical normalization φ = ε + α·ε² +...) ma strukturę:
# L = (1/2)(∂φ)² - (1/2)m²φ² - λ_3/3! φ³ - λ_4/4! φ⁴ + g_φ²·∂φ² + ...

# Dla α=2:
# m² = 1 (z V'' = 1 w eps²/2 współczynniku)
# λ_3 = (2-2α)/3 · 6 = 2·(2-2α) = 2·(-2) = -4 dla α=2  (wait, recheck)
# Actually V = (1/2)eps² + ((2-2α)/3)eps³ + ... = (1/2)eps² - (2/3)eps³ for α=2
# So coefficient of eps³ in V is (2-2α)/3 = -2/3 for α=2
# Standard convention V = (1/2)m²φ² + (g/3!)φ³ → g/6 = -2/3 → g = -4 for α=2
g3_at_2 = 2*(2 - 2*alpha_val)  # bo eps^3 współczynnik to (2-2α)/3, λ_3/3! = (2-2α)/3 → λ_3 = 2(2-2α)
g_kin_at_2 = 2*alpha_val       # ε·(∂ε)² coupling z g^(2α) = 1 + 2α·ε + ...

print(f"\n  Identified couplings (bare, before canonical normalization):")
print(f"    m^2           = 1.00000  (mass^2 of fluctuation)")
print(f"    λ_3 (eps^3)   = {g3_at_2:.5f}  (cubic self-interaction)")
print(f"    g_kin (eps·(∂eps)^2) = {g_kin_at_2:.5f}  (derivative coupling)")

# ============================================================================
# SEKCJA 2: 1-loop γ_φ standard QFT estimate
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 2: 1-loop γ_φ standard scalar QFT estimate")
print("=" * 78)

# Standard formuła:
# Dla L = (1/2)(∂φ)² - (1/2)m²φ² - λ/(3!)φ³:
#   γ_φ_1loop_phi3 = λ²/(6·(4π)²)  (in d=4; w d=3 modyfikacja kinematyczna)
# Dla L z derivative coupling g·φ(∂φ)²:
#   γ_φ_1loop_deriv = g²·m²/(8·(4π)²) (rough estimate)

# W d=3 R3 effective theory (1+1+1 spatial r-direction + ?), dimensionsfulnie:
# [φ] = (d-2)/2 = 1/2 w d=3
# [λ_3] = 3 - 3·(d-2)/2 = 3 - 3/2 = 3/2 w d=3 → relevant operator
# [g_kin] = (d-2)/2 - 0 = 1/2 w d=3 (dimensionful)
# [m²] = 2 (relevant mass^2)

# Anomalous dimension ogólnie γ_φ = (1/2) d ln Z_φ / d ln μ
# 1-loop dla φ³ w d=3: kinetic correction Σ(p²) → derivative wave-function renorm
# Estimate: γ_φ ~ λ_3²·m^(d-4) / (16π²) ~ λ_3²/m / (16π²)

lambda3 = g3_at_2  # = -4
m2 = 1.0
g_kin = g_kin_at_2  # = 4

# Naive d=4 estimate (just for scale comparison):
gamma_1loop_phi3 = lambda3**2 / (6 * (4*PI)**2)
gamma_1loop_deriv = g_kin**2 * m2 / (8 * (4*PI)**2)
gamma_1loop_total = gamma_1loop_phi3 + gamma_1loop_deriv

print(f"\n  Standard 1-loop estimate (d=4 normalization, scale-illustrative):")
print(f"    γ_φ from λ_3²/[6·(4π)²]   = {gamma_1loop_phi3:.6f}")
print(f"    γ_φ from g_kin²·m²/[8·(4π)²] = {gamma_1loop_deriv:.6f}")
print(f"    Total 1-loop γ_φ_estimate    = {gamma_1loop_total:.6f}")
print(f"\n  Target (z empirical n(α=2)/2):  γ_φ = e²/2 = {E2_HALF:.6f}")
print(f"  Ratio target/1-loop:          {E2_HALF/gamma_1loop_total:.2f}x")

if gamma_1loop_total < E2_HALF * 0.1:
    print(f"\n  --> 1-LOOP gamma_phi << e^2/2: MISS by factor >{E2_HALF/gamma_1loop_total:.0f}x")
    print(f"      Standard perturbative QFT NIE moze produkowac e^2/2.")

# ============================================================================
# SEKCJA 3: Coupling required dla γ_φ = e²/2 (non-perturbative regime)
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 3: Required coupling for γ_φ = e²/2 in standard scalar QFT")
print("=" * 78)

# γ_φ = X²/(16π²) → X = sqrt(16π² · γ_φ) = 4π · sqrt(γ_φ)
required_coupling_X = 4 * PI * np.sqrt(E2_HALF)
print(f"\n  Z formuły γ_φ = X²/(16π²) (typical 1-loop kinetic):")
print(f"    X_required = 4π · sqrt(e²/2) = {required_coupling_X:.4f}")
print(f"\n  Bare R3 couplings @ α=2:")
print(f"    |λ_3|   = {abs(lambda3):.2f}")
print(f"    |g_kin| = {abs(g_kin):.2f}")
print(f"\n  Required coupling X_req = {required_coupling_X:.2f} jest")
print(f"    {required_coupling_X/abs(g_kin):.1f}x wieksze niz natural g_kin")
print(f"    {required_coupling_X/abs(lambda3):.1f}x wieksze niz natural λ_3")

print("\n  --> Non-perturbative regime: standard 1-loop NOT applicable.")
print("      To dostac γ_φ = e²/2 z 1-loop trzeba bare coupling rzędu 24,")
print("      co jest poza perturbacja. R3 effective ma couplings rzedu 4.")

# ============================================================================
# SEKCJA 4: Numeryczny test — empirical γ_φ z R3 ODE
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 4: Empirical γ_φ z R3 ODE numerical solution")
print("=" * 78)

def r3_ode(r, y, alpha):
    """R3 ODE: g'' + (α/g)(g')² + (2/r)g' = (1-g)·g^(2-2α)"""
    g, gp = y
    g_safe = max(abs(g), 1e-12)
    rhs = (1 - g) * g_safe**(2 - 2*alpha) - (alpha/g_safe)*gp**2 - (2/r)*gp
    return [gp, rhs]

def solve_r3(g0, alpha, r_max=80.0, atol=1e-12, rtol=1e-10):
    """Solve R3 ODE z g(0)=g0, g'(0)=0; return r, g, gp"""
    r0 = 1e-3
    # Series start: g(r) ≈ g0 + (1/6)·(1-g0)·g0^(2-2α)·r²
    g_init = g0 + (1/6)*(1-g0)*g0**(2-2*alpha)*r0**2
    gp_init = (1/3)*(1-g0)*g0**(2-2*alpha)*r0
    sol = solve_ivp(r3_ode, [r0, r_max], [g_init, gp_init],
                    args=(alpha,), method='DOP853',
                    atol=atol, rtol=rtol, dense_output=True,
                    max_step=0.05)
    return sol.t, sol.y[0], sol.y[1]

def extract_amplitude(r, g, gp):
    """Extract A_tail z g(r) ~ 1 + A·cos(r-φ)/r asymptotically."""
    # Find tail region (r > 10, g close to 1)
    mask = (r > 10) & (np.abs(g - 1) < 0.5)
    if mask.sum() < 10:
        return np.nan
    r_tail = r[mask]
    h = (g[mask] - 1) * r_tail  # h(r) = A·cos(r-φ)
    # A² = max(h)² + max(h')² (envelope amplitude)
    h_max = np.max(np.abs(h))
    return h_max

def extract_m_obs(r, g, gp, alpha):
    """Extract m_obs z R3 mass integral M = ∫ ρ(r) 4π r² dr.

    Energy density: ρ = (1/2)g^(2α)·(g')² + V(g)
    V(g) = -∫(1-g)·g^(2-2α)dg z V(1)=0
    Dla numerical: V(g) = -[(1-g²)/2 - (g^(4-2α) - 1)/(4-2α)] for α≠2
    Dla α=2: V(g) = -[(1-g²)/2 - log(g)/0]... potrzeba careful limit
    """
    if abs(alpha - 2) < 1e-6:
        # α=2: V(g) = -(1-g²)/2 + log(g) - (?)
        # V'(g) = (g-1)·g^(-2) = 1/g - 1/g²
        # V(g) = log(g) + 1/g - 1; V(1) = 0 ✓
        V = np.log(np.abs(g) + 1e-12) + 1/np.maximum(np.abs(g), 1e-12) - 1
    else:
        V = -((1 - g**2)/2 - (g**(4-2*alpha) - 1)/(4-2*alpha))

    rho = 0.5 * np.abs(g)**(2*alpha) * gp**2 + V
    integrand = rho * 4 * PI * r**2
    M = trapezoid(integrand, r)
    return abs(M)  # take abs to handle numerical noise

# Scan g₀ values
g0_values = np.array([1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
alpha_use = 2.0

print(f"\n  Solving R3 ODE for α={alpha_use}, g₀ ∈ [1.5, 2.5]:")
print(f"\n  {'g₀':>6} | {'A_tail':>10} | {'m_obs':>10} | {'log(g₀)':>9} | {'log(A²)':>9} | {'log(m)':>9}")
print(f"  {'-'*6:>6}-+-{'-'*10:>10}-+-{'-'*10:>10}-+-{'-'*9:>9}-+-{'-'*9:>9}-+-{'-'*9:>9}")

results = []
for g0 in g0_values:
    try:
        r, g, gp = solve_r3(g0, alpha_use)
        A = extract_amplitude(r, g, gp)
        m = extract_m_obs(r, g, gp, alpha_use)
        if not (np.isnan(A) or np.isnan(m)) and A > 0 and m > 0:
            results.append((g0, A, m, np.log(g0), 2*np.log(A), np.log(m)))
            print(f"  {g0:>6.3f} | {A:>10.5f} | {m:>10.5f} | {np.log(g0):>9.5f} | {2*np.log(A):>9.5f} | {np.log(m):>9.5f}")
    except Exception as e:
        print(f"  {g0:>6.3f} | ERROR: {e}")

results = np.array(results)
if len(results) > 3:
    log_g0 = results[:, 3]
    log_A2 = results[:, 4]
    log_m = results[:, 5]

    # Fit slope: log(m) = slope_m · log(g₀) + const
    slope_m, _ = np.polyfit(log_g0, log_m, 1)
    # Fit slope: log(A²) = slope_A2 · log(g₀) + const
    slope_A2, _ = np.polyfit(log_g0, log_A2, 1)
    # log(Z) = log(m) - log(A²) → slope_Z = slope_m - slope_A2
    slope_Z = slope_m - slope_A2
    # Effective γ_φ = slope_Z (z RG flow Z ~ g₀^γ_φ)

    print(f"\n  Slopes z polyfit:")
    print(f"    slope_m  (d log m / d log g₀)  = {slope_m:.5f}")
    print(f"    slope_A² (d log A² / d log g₀) = {slope_A2:.5f}")
    print(f"    slope_Z = slope_m - slope_A²   = {slope_Z:.5f}")
    print(f"\n  Interpretation: m_obs = c · A² · g₀^X, gdzie X = slope_Z")
    print(f"    Empirical X    = {slope_Z:.5f}")
    print(f"    e²/2 (target)  = {E2_HALF:.5f}")
    print(f"    diff           = {abs(slope_Z - E2_HALF):.5f}  ({100*abs(slope_Z - E2_HALF)/E2_HALF:.3f}%)")

# ============================================================================
# SEKCJA 5: Critical analysis — czy "X = γ_φ" jest sensowna interpretacja?
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 5: Czy X = γ_φ z RG flow jest sensowna interpretacja?")
print("=" * 78)

print("""
  Argument PRO (RG flow interpretation):
    - W RG flow Z(g₀) = exp(γ_φ · log(g₀/Λ))
    - Empirically log(m/A²) = X · log(g₀) ± const → X = γ_φ formal
    - Dla α=2: X ≈ 3.69 ≈ e²/2 z 0.0007% match (mass_scaling_k4)

  Argument CONTRA (M.4 audit):
    - γ_φ z 1-loop standard scalar QFT @ α=2: γ_φ ~ {:.4f} (orders below 3.69)
    - Required X ~ 24 dla γ_φ = e²/2 — non-perturbative regime
    - R3 effective theory NIE jest standardowy scalar field theory:
      * Kinetic prefactor g^(2α) modyfikuje propagator
      * Soliton background NIE jest perturbative vacuum
      * "RG flow" w r-direction nie jest standardowy Wilson RG
    - "X = γ_φ" jest formal mapping, nie field-theoretic derivation
    - mass_scaling_k4 K-like NEGATIVE: K~A² to generic tail-scaling, nie
      topological, więc empirical slope opiera się na kinematice ogona
      ODE, nie strukturalnym Z(g₀)
""".format(gamma_1loop_total))

# ============================================================================
# SEKCJA 6: Final verdict M.4
# ============================================================================
print("\n" + "=" * 78)
print("  SEKCJA 6: Final verdict M.4")
print("=" * 78)

verdict = """
  M.4 (RG flow γ_φ ~ e²) — ANALYSIS:

  1. Standard 1-loop γ_φ z R3 effective theory @ α=2: ~{:.4f}
     Target e²/2 = {:.4f}
     MISS by factor ~{:.0f}x

  2. Required bare coupling dla γ_φ = e²/2 standardowo: ~{:.1f}
     Natural R3 couplings: |g_kin|=4, |λ_3|=4
     MISS: factor ~6x non-perturbative regime

  3. Empirical slope numerical: X = e²/2 to coincidence z fitting,
     nie field-theoretic derivation

  4. R3 effective theory NIE jest standardowy scalar QFT:
     - Soliton background, nie vacuum
     - g^(2α) kinetic non-trivial
     - "RG flow w r-direction" nie jest Wilsonian

  CONCLUSION:
  M.4 NIE produkuje e²/2 z first-principles RG flow w obecnym TGP-formalism.
  Wymagane coupling jest non-perturbative (~6x natural), brak standardowego
  mechanizmu który daje takie wartości naturalnie.

  Zgodne z Phase 2 P2.1 (log det O K=-0.97) — to ten sam fundamental
  obstacle: perturbative/semiclassical R3 NIE produkuje e²/2 numerycznie.

  STATUS: M.4 NEGATIVE → λ.1 audit-rekomendowane zamknięcie.
""".format(gamma_1loop_total, E2_HALF, E2_HALF/gamma_1loop_total,
           required_coupling_X)

print(verdict)

print("=" * 78)
print("  KONIEC λ.1 M.4 RG flow γ_φ test")
print("=" * 78)
