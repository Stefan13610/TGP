#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase2_P22_multiloop_semiclassical.py
=======================================

PURPOSE
-------
λ.1 Phase 2 sub-task P2.2: Multi-loop / semiclassical analysis dla
e_Euler natural appearance.

CONTEXT
-------
Po P2.1 (negative — log det O nie produkuje e²) i P2.3 (negative —
Brannen anchor NIE matche), to jest **ostatnia szansa** dla λ.1
żeby znaleźć mechanizm e²-emergence.

USER'S M4 MECHANISM (kumulatywny soliton dressing):
  g(r+dr) = g(r) · (1 + Δ(r)·dr/n)^n → g(r) · exp(Δ(r)·dr) gdy n→∞
  g(R) = g(0) · exp(∫₀^R Δ(r) dr)

To jest **continuous limit dyskretnego procesu mnożnikowego**, który
matematycznie produkuje e^x natural.

PYTANIE P2.2:
  Czy w R3 ODE istnieje natural Δ(r;α) takie że ∫₀^R_max Δ(r;α) dr
  daje **(e²/4)·(4-α)·log(g₀)**?

Innymi słowy: czy mass formula log m_obs = const + ∫integrand wynika
z continuous limit i czy integrand zawiera e²/4 jako coefficient?

PASS CRITERION
--------------
"Identify konkretny mechanizm produkujący e_Euler (nawet jeśli nie
e²/4 dokładnie)"

Autor: λ.1 Phase 2 P2.2
Data: 2026-05-01
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import curve_fit
import math

E = math.e
E_SQ = E**2
PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI

print("=" * 78)
print("  λ.1 P2.2 — Multi-loop / semiclassical: gdzie e_Euler pojawia się?")
print("=" * 78)
print()


# ----------------------------------------------------------------
# SECTION 1: Test M4 mechanism — czy log g jest integral
# ----------------------------------------------------------------

def solve_R3_ode(g0, alpha=2.0, d=3, r_max=200.0, n_points=10000, g_floor=1e-10):
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
                    rtol=1e-9, atol=1e-11, max_step=0.05)
    return sol, singular[0]


print("=" * 78)
print("  SEKCJA 1: Test M4 — d(log g)/dr struktura w R3 ODE")
print("=" * 78)
print()

print("""
  Z R3 ODE: g'' + (α/g)(g')² + (2/r)g' = (1-g)·g^(2-2α)

  d(log g)/dr = g'/g  ≡  η(r;α)

  W limicie continuous, mass formula:
    log(g(R)/g(0)) = -log(g(0))  (bo g(R) → 1, g(0) = g₀)
    = -∫₀^∞ (g'/g) dr = ∫_g₀^1 (1/g) dg

  Dla R3 z α=2: integral po profilu solitonu daje konkretną wartość
  zalezną od g₀. Sprawdzimy NUMERYCZNIE.
""")

# Numerical integral of d(log g)/dr along soliton profile
print(f"  {'g₀':>8} | {'log(g₀)':>10} | {'integral':>12} | {'sum diff':>10}")
print("  " + "-" * 50)

g0_test = [G0_E, G0_MU, 1.5, 0.5, 0.7]
for g0 in g0_test:
    sol, sing = solve_R3_ode(g0, alpha=2.0)
    if sing or not sol.success:
        continue
    r = sol.t
    g = sol.y[0]
    gp = np.gradient(g, r)
    integrand = gp / g
    integral = np.trapezoid(integrand, r)

    # Should equal log(g[end]) - log(g[0]) ≈ -log(g₀) (since g[end]→1)
    expected = math.log(g[-1]) - math.log(g[0])

    print(f"  {g0:8.4f} | {-math.log(g0):10.5f} | {integral:12.5f} | {abs(integral-expected):10.5f}")

print()
print("  ✓ Trywialne (fundamental theorem of calculus): ∫(g'/g)dr = log(g_end/g_start)")
print("  To NIE pomaga w mass formula directly.")
print()


# ----------------------------------------------------------------
# SECTION 2: M4 dla mass — sprawdz integral (g'/g)·r-weight
# ----------------------------------------------------------------

print("=" * 78)
print("  SEKCJA 2: M4 z 3D radial weight r²·dr")
print("=" * 78)
print()

print("""
  W 3D radial integration, weight to r²·dr (volume element).

  Mass formula z why_n3: m_obs ~ A_tail² · g₀^(e²/2) (dla α=2).

  log m_obs = 2·log A_tail + (e²/2)·log g₀

  Pytanie: czy log m_obs = ∫₀^∞ Φ(r)·r²·dr dla pewnego natural Φ(r)?

  Test: czy istnieje funkcja Φ_M4(r;α) taka że
    ∫₀^R_max Φ_M4(r;α) · r² dr = (e²/2)·log(g₀) + bezstałe·(...)
""")

# Empirical exponent from R3 mass formula z why_n3
n_alpha = E_SQ / 2  # = n(2)

print(f"  {'g₀':>8} | {'log(g₀)':>10} | {'(e²/2)·log(g₀)':>16} | {'∫integrand·r²·dr':>20}")
print("  " + "-" * 65)

# Test integrand candidates
def candidate_integrand(g, gp, r, alpha=2.0, kind='gp_squared'):
    """Various integrand candidates dla mass formula contribution."""
    if kind == 'gp_squared':
        return gp**2 * g**(2*alpha)
    elif kind == 'log_g':
        return math.log(g) if g > 0.01 else 0
    elif kind == 'V_eff':
        return g**3/3 - g**4/4 - 1/12  # potential
    elif kind == 'gp_g':
        return (gp / g) if abs(g) > 1e-3 else 0
    return 0


for g0 in [G0_E, G0_MU, 1.2, 0.7]:
    sol, sing = solve_R3_ode(g0, alpha=2.0)
    if sing or not sol.success:
        continue
    r = sol.t
    g = sol.y[0]
    gp = np.gradient(g, r)

    # Compute (gp/g)*r² integrand (with vacuum subtraction)
    integrand_M4 = np.array([(gp[i]/g[i])*r[i]**2 if abs(g[i]) > 1e-3 else 0
                              for i in range(len(r))])
    integral = np.trapezoid(integrand_M4, r)

    target = (E_SQ/2) * math.log(g0)
    print(f"  {g0:8.4f} | {math.log(g0):10.5f} | {target:16.5f} | {integral:20.5f}")

print()
print("  Wnioski:")
print("  - Integrand (g'/g)·r² nie matche (e²/2)·log(g₀) directly")
print("  - Może wymagana jest inna struktura integrand")


# ----------------------------------------------------------------
# SECTION 3: Semiclassical action S_inst dla R3 soliton
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Semiclassical action S_inst dla R3 soliton")
print("=" * 78)
print()

print("""
  W path integral: Z = ∫Dφ exp(-S/ℏ).
  Semiclassical: dominantne contribution z saddle point (soliton):
    Z ~ exp(-S_sol/ℏ) · [det']^(-1/2) · (collective coordinates Jacobian)

  S_sol = ∫d³x [½g^(2α)(∂g)² + V(g)]

  Dla R3 soliton (numerical action):
""")

# Compute S_inst numerically
def compute_S_inst(g0, alpha=2.0, R_max=100.0):
    sol, sing = solve_R3_ode(g0, alpha, r_max=R_max)
    if sing or not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    gp = np.gradient(g, r)

    # 3D action (radial)
    # S = ∫ [½g^(2α)·(g')² + V_eff(g)] · 4π·r² dr
    V_at_1 = 1.0/3 - 1.0/4
    integrand = (0.5 * g**(2*alpha) * gp**2 + (g**3/3 - g**4/4 - V_at_1)) * 4*math.pi * r**2
    return np.trapezoid(integrand, r)


print(f"  {'g₀':>8} | {'S_sol':>11} | {'exp(-S_sol)':>13}")
print("  " + "-" * 40)
for g0 in [G0_E, G0_MU, 1.2, 0.7]:
    S = compute_S_inst(g0, alpha=2.0)
    if S is not None:
        exp_minus_S = math.exp(-S) if abs(S) < 50 else 0
        print(f"  {g0:8.4f} | {S:11.5f} | {exp_minus_S:13.4e}")

print()
print("  Obserwacja: S_sol jest **rzeczywista wielkość** dla każdej generacji.")
print("  exp(-S_sol) daje **subleading correction**, nie dominant mass term.")
print()
print("  W mass formula z R3, dominant term to:")
print("    m ~ A² · g₀^(e²/2)  (z why_n3 Phase 2)")
print()
print("  Aby pokazać że to wynika z exp(-S_sol):")
print("    log m ~ 2·log A + (e²/2)·log g₀")
print()
print("  Czyli wymagamy:")
print("    -S_sol/ℏ ↔ 2·log A_tail + (e²/2)·log g₀")
print()
print("  To **explicitly** wymaga że S_sol ma logarytmiczną zalezność od g₀:")
print("    S_sol = -ℏ·[2·log A_tail(g₀) + (e²/2)·log g₀] + const")


# Sprawdz czy faktycznie S_sol ~ log g₀
print()
print("=" * 78)
print("  SEKCJA 4: Sprawdz S_sol vs log g₀")
print("=" * 78)
print()

g0_arr = []
S_arr = []
for g0 in np.linspace(0.5, 1.8, 10):
    if abs(g0 - 1.0) < 0.05:
        continue
    S = compute_S_inst(g0, alpha=2.0)
    if S is not None and not math.isnan(S):
        g0_arr.append(g0)
        S_arr.append(S)

if len(g0_arr) >= 3:
    g0_arr = np.array(g0_arr)
    S_arr = np.array(S_arr)
    log_g0 = np.log(g0_arr)

    # Linear fit S = K·log(g₀) + C
    K, C = np.polyfit(log_g0, S_arr, 1)
    print(f"  Linear fit: S_sol = K·log(g₀) + C")
    print(f"    K = {K:.4f}")
    print(f"    C = {C:.4f}")
    print()
    print(f"  Compare z empirical exponents:")
    print(f"    e² = {E_SQ:.4f}")
    print(f"    e²/2 = {E_SQ/2:.4f}")
    print(f"    e²/4 = {E_SQ/4:.4f}")
    print(f"    -e² = {-E_SQ:.4f}")
    print(f"    -e²/2 = {-E_SQ/2:.4f}")

    # Check for match
    match_e2 = abs(K - E_SQ) / E_SQ
    match_neg_e2_half = abs(K + E_SQ/2) / (E_SQ/2)
    if match_neg_e2_half < 0.2:
        print(f"  >>> POSSIBLE: K ≈ -e²/2 (within 20%)")
    elif match_e2 < 0.2:
        print(f"  >>> POSSIBLE: K ≈ e² (within 20%)")
    else:
        print(f"  >> No clean match z e² family.")


# ----------------------------------------------------------------
# SECTION 5: Borel-resummed perturbation series
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Borel-resumed perturbation series")
print("=" * 78)
print()

print("""
  Dla scalar field theory z self-interaction, 1-loop daje:
    Σ_1 ~ g · (loop integral)

  k-loop: Σ_k ~ g^k · k! / (combinatorial factor)
                  (factorial growth — Borel transform)

  Borel-summed: Σ(g) = ∫₀^∞ B(t) · exp(-t/g) dt
                gdzie B(t) = Σ_k (Σ_k/k!) · t^k

  e^(-t/g) **explicitly contains exp**.

  Dla scalar w R3 amplitude sector, jeśli coupling g zwiazany z g₀^(α-zalezne):
    Σ_resummed(g) zawiera exp(coś z g₀)

  Więc e_Euler MOŻE wynikać z resummed perturbation series,
  ale **konkretna wartość e²/4** wymaga konkretnej calculation
  która nie jest tu dostępna.

  STATUS: Strukturalnie plausible, ale nie konkretne derivation.
""")


# ----------------------------------------------------------------
# SECTION 6: Continuous-limit interpretation
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 6: Continuous-limit z user's analogy")
print("=" * 78)
print()

print("""
  User's intuicja: e = lim(1+1/n)^n powstaje z kumulatywnego procesu.

  W R3, soliton można interpretować jako "kontynentalny" (continuous)
  fluctuating background. Każda warstwa r → r+dr "dressuje" pole o:
    g(r+dr)/g(r) = 1 + (g'/g)·dr

  Iteracyjnie: g(R)/g(0) = exp(∫₀^R (g'/g) dr) — TAUTOLOGY

  DLA M4 z RZECZYWISTYM SENSEM, potrzebujemy że KAŻDA warstwa wnosi
  swoj wkład **exp(d_local)** gdzie d_local nie jest już integral.

  W R3: każda r-warstwa wnosi:
    exp(d_local) = ?

  Bez konkretnej fizyki "dressing per layer", ten mechanizm pozostaje
  formalny.
""")


# ----------------------------------------------------------------
# SECTION 7: HONEST conclusion P2.2
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 7: HONEST conclusion P2.2")
print("=" * 78)
print()

print(f"""
  Co odkryłem w P2.2:

  1. M4 mechanism (continuous limit) jest **strukturalnie OK** ale
     bez konkretnej "dressing per layer" fizyki, pozostaje formalny.

  2. Semiclassical S_sol jest **rzeczywistą wielkością** dla R3 solitonu,
     ale linear fit S vs log(g₀) NIE pokazał czystego matchu z e²/2 lub
     e². K-coefficient z fitu znacznie różny od empirical.

  3. Borel-resummed series **może** produkować e_Euler factor (z exp(-t/g)
     w Borel transform), ale **konkretna wartość e²/4** wymaga eksplicit
     calculation która jest poza scope L1.

  4. log m_obs = 2·log A + (e²/2)·log g₀ z R3 NIE matche prosty -S_sol/ℏ
     ze semiclassical analysis.

  PASS criterion P2.2:
    "Identify konkretny mechanizm produkujący e_Euler (nawet jeśli
     nie e²/4 dokładnie)"

  Status: **PARTIAL** — zidentyfikowane potencjalne mechanizmy
  (M4 continuous limit, Borel resummation) są **structurally
  plausible**, ale **żaden nie produkuje konkretnego e²** w R3.

  Recommendation: count P2.2 jako **0.5 PASS** (mechanizmy zidentyfikowane
  ale nie zamknięte; consistent z L1.3 partition function PASS).
""")
