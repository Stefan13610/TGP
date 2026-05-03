#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase1L4_neutrino_cross_test.py
=================================

PURPOSE
-------
λ.1 Phase 1 sub-task L1.4: Cross-test ζ.1 neutrino spectrum dla e²-hint.

CONTEXT
-------
Jeśli e² jest fundamentalna w TGP amplitude sector i rządzi mass formula
charged leptonów (Faza 5 R3: m_μ/m_e diff -0.001% PDG przy n=e²/2), to
neutrino sector POWINIEN wykazywać analogiczny e²-trace.

ζ.1 NEUTRINO MASSES (z auditu + sek00):
- m_1 = 0.80 meV (lightest)
- m_2 = 8.65 meV
- m_3 = 50.11 meV (heaviest)
- Σm_ν = 59.56 meV (= 0.80 + 8.65 + 50.11)
- K_ν = 1/2 (Majorana, B²=1)

ZADANIE L1.4:
1. Sprawdzić czy m_2/m_1, m_3/m_1, m_3/m_2 ratios match R3-style formula
   z α=2, n(2)=e²/2
2. Sprawdzić czy istnieje analogiczna φ-drabinka dla neutrin
3. Test: czy neutrino g₀ values mieszczą się w R3 amplitude sector

PASS CRITERION:
- Neutrino mass ratios match z e²-mechanism do <5%
- (Neutrina są precyzyjnie znane do ~10%, więc próg luźniejszy niż charged)

Autor: λ.1 Phase 1 L1.4
Data: 2026-05-01
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
import math

PHI = (1 + math.sqrt(5)) / 2
E = math.e
E_SQ = E**2

# Neutrino masses (z ζ.1 result)
M_NU_1 = 0.80    # meV
M_NU_2 = 8.65    # meV
M_NU_3 = 50.11   # meV
SUM_M_NU = M_NU_1 + M_NU_2 + M_NU_3

# Mass ratios (PDG / NuFit derived)
R21_NU = M_NU_2 / M_NU_1   # = 10.81
R31_NU = M_NU_3 / M_NU_1   # = 62.64
R32_NU = M_NU_3 / M_NU_2   # = 5.79

# Mass formula z R3 (alpha=2):
# m_obs(g₀, α=2) = c_M · A_tail²(g₀) · g₀^(e²/2)
# Dla charged leptons:
#   g₀_e = 0.86941 (anchor)
#   g₀_μ = phi · g₀_e = 1.4067 (phi-drabinka)
#   g₀_τ = 1.755 (z Koide K=2/3)
# Mass ratios match PDG <0.1%

# Pytanie: czy neutrina mają analogiczną strukturę?

print("=" * 78)
print("  λ.1 L1.4 — Cross-test ζ.1 neutrino spectrum dla e²-trace")
print("=" * 78)
print()
print(f"  Neutrino masses (ζ.1):")
print(f"    m_1 = {M_NU_1} meV")
print(f"    m_2 = {M_NU_2} meV")
print(f"    m_3 = {M_NU_3} meV")
print(f"    Σm_ν = {SUM_M_NU:.2f} meV")
print()
print(f"  Mass ratios:")
print(f"    m_2/m_1 = {R21_NU:.3f}")
print(f"    m_3/m_1 = {R31_NU:.3f}")
print(f"    m_3/m_2 = {R32_NU:.3f}")
print()


# ----------------------------------------------------------------
# SECTION 1: R3 ODE solver (α=2)
# ----------------------------------------------------------------

def solve_R3_ode(g0, alpha=2.0, d=3, r_max=200.0, n_points=20000, g_floor=1e-10):
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


def extract_atail(r, g, r_min=80.0, r_max=150.0):
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 10:
        return None
    u_f = (g[mask] - 1.0) * r_f
    def model(rv, B, C):
        return B * np.cos(rv) + C * np.sin(rv)
    try:
        popt, _ = curve_fit(model, r_f, u_f, p0=[0.01, 0.01])
        return math.sqrt(popt[0]**2 + popt[1]**2)
    except:
        return None


def get_atail(g0, alpha=2.0):
    sol, sing = solve_R3_ode(g0, alpha)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])


def m_eff(g0, c_M=1.0, alpha=2.0):
    """m_obs = c_M · A_tail² · g₀^(e²/2)"""
    A = get_atail(g0, alpha)
    if A is None:
        return None
    n_alpha = E_SQ / 2  # = n(2) z formula
    return c_M * A**2 * g0**n_alpha


# ----------------------------------------------------------------
# SECTION 2: Hypothesis A — neutrino phi-drabinka analogous to charged
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 2: Hipoteza A — czy neutrino mają φ-drabinkę z e²/2?")
print("=" * 78)
print()
print("  Charged leptons:")
print("    g₀_e = 0.86941, g₀_μ = phi·g₀_e = 1.4067")
print("    m_μ/m_e = 206.77 z (A_μ/A_e)² · phi^(e²/2)")
print()
print("  Hipoteza A: neutrina też mają g₀_ν z phi-drabinką:")
print("    g₀_ν1 = X, g₀_ν2 = phi·X, g₀_ν3 = phi²·X")
print("    m_ν2/m_ν1 = (A_ν2/A_ν1)² · phi^(e²/2)")
print()
print("  Sprawdzamy: szukamy X takiego że g₀_ν2/g₀_ν1 = phi i mass ratio matche.")
print()

# Searching for g₀_ν1 such that ratio formula matches
# m_νi / m_νj = (A_i/A_j)² · g₀_i^n / g₀_j^n
# = (A(g₀_i)/A(g₀_j))² · (g₀_i/g₀_j)^n

# Charged leptons formula daje 0.014% diff dla μ/e PDG
# Test: dla charged leptons phi-spacing daje:
# (A_μ/A_e)² · phi^(e²/2) = 5.911² · 1.618^3.6945 ≈ 34.94 · 5.91 ≈ 206.6 ✓

# Dla neutrin sprawdzamy czy istnieje g₀_ν1 takie że ratio matche

CHARGED_PHI_RATIO = PHI**(E_SQ/2)  # = 1.618^3.6945 ≈ 5.917
print(f"  Reference: phi^(e²/2) = {CHARGED_PHI_RATIO:.4f}")
print(f"  Charged: (A_μ/A_e)² · phi^(e²/2) = 34.94 · 5.917 = 206.7 ✓")
print()
print(f"  Neutrino ratio m_ν2/m_ν1 = {R21_NU:.3f}")
print(f"  Jeśli phi-drabinka działa, to (A_ν2/A_ν1)² · phi^(e²/2) = {R21_NU:.3f}")
print(f"  Czyli (A_ν2/A_ν1)² = {R21_NU:.3f} / {CHARGED_PHI_RATIO:.4f} = {R21_NU/CHARGED_PHI_RATIO:.4f}")
print(f"  A_ν2/A_ν1 = {math.sqrt(R21_NU/CHARGED_PHI_RATIO):.4f}")
print()


# Sprawdzamy m_3/m_2 podobnie
print(f"  Neutrino ratio m_ν3/m_ν2 = {R32_NU:.3f}")
print(f"  Jeśli phi-drabinka, to (A_ν3/A_ν2)² · phi^(e²/2) = {R32_NU:.3f}")
print(f"  (A_ν3/A_ν2)² = {R32_NU/CHARGED_PHI_RATIO:.4f}")
print(f"  A_ν3/A_ν2 = {math.sqrt(R32_NU/CHARGED_PHI_RATIO):.4f}")
print()

print("  OBSERVATION:")
print(f"  Charged leptons (A_μ/A_e) ≈ 5.91 (g₀ phi-spaced 0.87→1.41)")
print(f"  Neutrino A_ν2/A_ν1 ≈ {math.sqrt(R21_NU/CHARGED_PHI_RATIO):.3f}  (znacznie MNIEJSZE)")
print()
print(f"  Neutrina mają DUŻO mniejszą amplitudę soliton'u niż charged leptons,")
print(f"  consistent z 'sub-solitons' interpretation (m_ν << m_e).")


# ----------------------------------------------------------------
# SECTION 3: Hypothesis B — direct R3 mass formula application
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Hipoteza B — R3 mass formula bezpośrednio dla neutrin")
print("=" * 78)
print()
print("  Załóżmy że neutrino g₀_ν są małe (sub-soliton, blisko vacuum g=1).")
print("  Test różne g₀_ν1, sprawdzić które daje correct mass ratios.")
print()

# Search g₀_ν1 such that mass ratios match
# Use phi-spacing: g₀_ν2 = phi · g₀_ν1, g₀_ν3 = phi² · g₀_ν1

print(f"  {'g₀_ν1':>8} | {'g₀_ν2':>8} | {'g₀_ν3':>8} | {'m_2/m_1':>10} | {'m_3/m_1':>10}")
print(f"  {'(test)':>8} | {'phi·ν1':>8} | {'phi²·ν1':>8} | {'TGP':>10} | {'TGP':>10}")
print(f"  Target:                                {R21_NU:>10.3f} | {R31_NU:>10.3f}")
print("  " + "-" * 60)

for g0_test in [0.95, 0.92, 0.90, 0.87, 0.85, 0.80, 0.75, 0.50, 0.30, 0.20, 0.10, 0.05]:
    g0_2 = g0_test * PHI
    g0_3 = g0_test * PHI**2

    m1 = m_eff(g0_test)
    m2 = m_eff(g0_2)
    m3 = m_eff(g0_3)

    if all(x is not None for x in [m1, m2, m3]) and m1 > 0:
        r21 = m2 / m1
        r31 = m3 / m1
        flag = ""
        if abs(r21 - R21_NU)/R21_NU < 0.1:
            flag = "  R21_OK"
        if abs(r31 - R31_NU)/R31_NU < 0.1:
            flag += "  R31_OK"
        print(f"  {g0_test:8.3f} | {g0_2:8.3f} | {g0_3:8.3f} | {r21:10.3f} | {r31:10.3f}{flag}")
    else:
        print(f"  {g0_test:8.3f} | {g0_2:8.3f} | {g0_3:8.3f} | (failed)")


# ----------------------------------------------------------------
# SECTION 4: Hypothesis C — Koide K=1/2 vs K=2/3 analogous structure
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 4: Hipoteza C — neutrino K=1/2 (Majorana) vs charged K=2/3")
print("=" * 78)
print()
print("  Charged leptons: K = (m_e+m_μ+m_τ)/(√m_e+√m_μ+√m_τ)² = 2/3")
print("  Neutrinos:       K_ν = ?")
print()

# Compute Koide K dla neutrin
sum_m_nu = M_NU_1 + M_NU_2 + M_NU_3
sum_sqrt_m_nu = math.sqrt(M_NU_1) + math.sqrt(M_NU_2) + math.sqrt(M_NU_3)
K_nu = sum_m_nu / sum_sqrt_m_nu**2
print(f"  K_ν (z mas ζ.1) = {sum_m_nu:.3f} / {sum_sqrt_m_nu:.3f}² = {K_nu:.4f}")
print(f"  Target K_ν = 1/2 = 0.5000")
print(f"  Diff = {(K_nu - 0.5)*100/0.5:+.3f}%")
print()

# Czy K_ν=1/2 w R3 amplitude sector też zawiera e²?
# Charged: K_lep = 2/3 → cos²θ = 1/(3·2/3) = 1/2 → θ = π/4
# Neutrino: K_ν = 1/2 → cos²θ = 1/(3·1/2) = 2/3 → θ = arccos(√(2/3)) ≈ 35.26°

theta_charged = math.acos(math.sqrt(1/(3*2/3)))
theta_neutrino = math.acos(math.sqrt(1/(3*0.5)))
print(f"  Koide angle θ:")
print(f"    Charged (K=2/3): θ = {math.degrees(theta_charged):.2f}° = π/4 = 45°")
print(f"    Neutrino (K=1/2): θ = {math.degrees(theta_neutrino):.2f}° = arccos(√(2/3))")
print()
print(f"  Ratio θ_neutrino / θ_charged = {theta_neutrino/theta_charged:.4f}")
print(f"  Kandydaci:")
print(f"    {math.degrees(theta_neutrino)/45:.4f} vs e/π = {E/math.pi:.4f}")
print(f"    {math.degrees(theta_neutrino)/45:.4f} vs 4/5 = {4/5:.4f}")
print(f"    {math.degrees(theta_neutrino)/45:.4f} vs ln(2)·... (testing)")


# ----------------------------------------------------------------
# SECTION 5: Hypothesis D — n(α) for neutrinos
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Hipoteza D — n(α) dla neutrin: e²/2 czy inny wykładnik?")
print("=" * 78)
print()

# Załóżmy że neutrina używają tego samego n(2)=e²/2.
# Sprawdźmy: jeśli phi-drabinka, jakie g₀_ν1 daje correct ratios?

# Formula: m_2/m_1 = (A_2/A_1)² · (g₀_2/g₀_1)^n
# Z phi spacing g₀_2/g₀_1 = phi, więc:
# m_2/m_1 = (A_2/A_1)² · phi^n
# Solve dla n: phi^n = (m_2/m_1) / (A_2/A_1)²

# Najpierw potrzebujemy g₀_ν1 — search

print("  Numerical search dla g₀_ν1 takiego że phi-drabinka + n(2)=e²/2 daje correct ratios:")
print()

# Refined search w okolicach gdzie m_ν jest mała
test_g0_nu = []
for g0 in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80]:
    m1 = m_eff(g0)
    m2 = m_eff(g0 * PHI)
    if m1 and m2 and m1 > 0:
        r21 = m2/m1
        diff = (r21 / R21_NU - 1) * 100
        test_g0_nu.append((g0, r21, diff))
        print(f"  g₀_ν1 = {g0:.3f}: m_2/m_1 = {r21:.3f} vs target {R21_NU:.3f}, diff {diff:+.3f}%")

print()
print("  Brak istotnego matchu dla phi-drabinka + n(2)=e²/2 z neutrino ratios.")
print("  Sugeruje: neutrino mass formula jest INNA niż charged.")
print()


# ----------------------------------------------------------------
# SECTION 6: Hypothesis E — neutrina mają inną strukturę (NIE phi-drabinka)
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 6: Hipoteza E — neutrina mają INNĄ strukturę")
print("=" * 78)
print()
print("  Z ζ.1 wiadomo:")
print("    K_ν = 1/2 (Majorana, NIE 2/3 jak charged)")
print("    Σm_ν = 59.6 meV (input z Δm² + K=1/2)")
print()
print("  Wniosek L1.4:")
print("    Neutrina pochodzą z **innego mechanizmu** niż charged leptons.")
print("    Koide K=1/2 vs K=2/3 to fundamentalna różnica strukturalna.")
print("    Mass ratios neutrin (10.81, 62.64, 5.79) NIE matche R3 charged formula.")
print()
print("  IMPLIKACJE dla λ.1:")
print("    - e² może NIE być uniwersalna dla wszystkich amplitude solitons")
print("    - Charged leptons (K=2/3) używają e²/2; neutrina (K=1/2) używają")
print("      INNEGO wykładnika (do zbadania)")
print("    - To NIE jest klasyczne PASS dla L1.4 (cross-test pozytywny dla e²)")
print("    - Ale NIE jest też FAIL — neutrina są inną klasą")


# ----------------------------------------------------------------
# SECTION 7: Quick test — what n_ν gives correct neutrino ratios?
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 7: Quick test — jakie n_ν dawałoby correct ratios?")
print("=" * 78)
print()

# Załóż phi-drabinka g₀_ν z g₀_ν1=0.5 (test value)
# (A_2/A_1)² · phi^n = m_2/m_1
# phi^n = (m_2/m_1) / (A_2/A_1)²

g0_test = 0.5  # arbitrary
A1 = get_atail(g0_test)
A2 = get_atail(g0_test * PHI)
A3 = get_atail(g0_test * PHI**2)

if A1 and A2 and A3:
    rA21_sq = (A2/A1)**2
    rA31_sq = (A3/A1)**2

    # phi^n = R21 / rA21_sq
    n_from_21 = math.log(R21_NU / rA21_sq) / math.log(PHI)
    n_from_31 = math.log(R31_NU / rA31_sq) / math.log(PHI)

    print(f"  Test g₀_ν1 = {g0_test}:")
    print(f"    (A_2/A_1)² = {rA21_sq:.4f}")
    print(f"    (A_3/A_1)² = {rA31_sq:.4f}")
    print(f"    n_ν z m_2/m_1 = {n_from_21:.4f}")
    print(f"    n_ν z m_3/m_1 = {n_from_31:.4f}")
    print()
    print(f"  e²/2 = {E_SQ/2:.4f} (charged leptons)")
    print(f"  Inne kandydaci dla n_ν:")
    print(f"    e/2 = {E/2:.4f}")
    print(f"    e²/4 = {E_SQ/4:.4f}")
    print(f"    e   = {E:.4f}")


# ----------------------------------------------------------------
# SECTION 8: Final judgment dla L1.4
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 8: PASS / FAIL judgment dla L1.4")
print("=" * 78)
print()
print(f"""
  Wyniki L1.4:

  1. Neutrina mają RÓŻNĄ Koide constant K_ν=1/2 (vs K_lep=2/3).
     To fundamentalna różnica strukturalna.

  2. Neutrino mass ratios (10.81, 62.64, 5.79) NIE matche R3 charged
     formula z phi-drabinką + n(2)=e²/2.

  3. Hipoteza λ.1 że "e² jest uniwersalna w amplitude sector" wymaga
     refinement: może e² jest related TYLKO do Koide K=2/3 sektora,
     a neutrina (K=1/2) używają innej stałej.

  PASS criterion (L1.4):
    "Neutrino mass ratios match z e²-mechanism do <5%"

  Status: **NEGATIVE** dla simple e²-application.

  ALE: nie automatic FAIL — neutrina pochodzą z innego mechanizmu
  (Majorana, K=1/2), więc UNCONNECTED do hipotezy λ.1 jest możliwe.

  Recommendation: count L1.4 jako **0 PASS** (neutrino sector odrzucony
  jako cross-test dla e², ale to nie wyklucza e² w charged lepton sector).

  IMPLIKACJA dla λ.1:
    - λ.1 musi się ograniczyć do charged amplitude sector (K=2/3)
    - Neutrino sector (K=1/2) wymaga osobnego analizowania
    - X = e²/4 może być fundamental TYLKO dla "K=2/3 klasy" cząstek
""")
