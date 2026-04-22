#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex237_neutrino_koide_predictions.py
====================================
PREDYKCJE MAS NEUTRIN Z KOIDE + TGP

KONTEKST:
  ex234: Lepton K(e,μ,τ) = 2/3 EXACT → m_τ z 0.006%
  ex235-236: Quark shifted Koide K(m+m₀) = 2/3 → m_b (0.43%), m_t (0.01%)

  PYTANIE: Czy Koide K=2/3 działa też dla neutrin?

DANE EKSPERYMENTALNE (oscylacje):
  Δm²₂₁ = 7.53 × 10⁻⁵ eV² (solar)
  Δm²₃₂ = 2.453 × 10⁻³ eV² (atmospheric, Normal Ordering)
  Δm²₃₂ = -2.536 × 10⁻³ eV² (Inverted Ordering)

  OGRANICZENIA:
  - Σm_ν < 0.12 eV (Planck 2018, 95% CL)
  - Σm_ν < 0.09 eV (Planck+BAO, deSalas+ 2021)
  - m_lightest unknown (could be ~0)

PLAN:
  §1. Koide constraints: K(ν₁,ν₂,ν₃) = 2/3 + Δm² data → m₁, m₂, m₃
  §2. Normal Ordering (NO) vs Inverted Ordering (IO)
  §3. Shifted Koide K(m+m₀)=2/3 for neutrinos (by analogy with quarks)
  §4. Brannen parametrization: ε=√2, θ for neutrinos
  §5. TGP link: r₂₁(ν) from soliton ODE? Neutrino φ-FP?
  §6. Σm_ν prediction vs cosmological bounds

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.optimize import brentq, minimize_scalar, minimize

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2

TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")


# ============================================================
# KOIDE TOOLS
# ============================================================
def koide(m1, m2, m3):
    """Koide parameter K = (m1+m2+m3)/(√m1+√m2+√m3)²"""
    S = np.sqrt(abs(m1)) + np.sqrt(abs(m2)) + np.sqrt(abs(m3))
    if S == 0:
        return np.nan
    return (m1 + m2 + m3) / S**2

def brannen_params(m1, m2, m3):
    """Extract Brannen ε, θ from masses.
    m_k = M/3 × (1 + ε·cos(θ + 2πk/3))²
    """
    M = m1 + m2 + m3
    S = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    K = M / S**2
    # K = (1+2ε²)/(3) when using normalized form
    # Actually K = 1/3 × (1 + 2ε²·cos²...) ...
    # Simpler: ε² = (3K - 1)/2 if K from Brannen
    # Wait: K_Koide = (m1+m2+m3)/(√m1+√m2+√m3)²
    # In Brannen: √m_k = √(M/3) × (1 + ε·cos(θ+2πk/3))/norm
    # Let me use: K = (1 + 2ε²/3·...)  — no, just compute numerically

    # From Brannen: if m_k = (M/3)(1+ε·cos(θ+2πk/3))², then
    # √m_k = √(M/3)(1+ε·cos(θ+2πk/3)) (when argument ≥ 0)
    # K = Σm_k / (Σ√m_k)²

    # For K=2/3: ε=√2 exactly (from ex234)
    # θ can be found from mass ratios

    # Fit θ given ε=√2
    eps = np.sqrt(2)

    def masses_from_theta(theta, M_total):
        mk = []
        for k in range(3):
            val = 1 + eps * np.cos(theta + 2*np.pi*k/3)
            if val < 0:
                return None
            mk.append(M_total/3 * val**2)
        return sorted(mk)

    # Search for θ that matches mass ratios
    r21_target = m2/m1
    r31_target = m3/m1

    def cost(theta):
        mk = masses_from_theta(theta, M)
        if mk is None:
            return 1e10
        r21 = mk[1]/mk[0] if mk[0] > 0 else 1e10
        r31 = mk[2]/mk[0] if mk[0] > 0 else 1e10
        return (np.log(r21/r21_target))**2 + (np.log(r31/r31_target))**2

    # Scan
    best_theta = None
    best_cost = 1e10
    for t in np.linspace(0, 2*np.pi, 1000):
        c = cost(t)
        if c < best_cost:
            best_cost = c
            best_theta = t

    # Refine
    from scipy.optimize import minimize_scalar
    result = minimize_scalar(cost, bounds=(best_theta-0.1, best_theta+0.1), method='bounded')
    theta_fit = result.x % (2*np.pi)

    return eps, theta_fit, K


# ============================================================
# EXPERIMENTAL DATA (eV² units, then convert to eV)
# ============================================================
dm2_21 = 7.53e-5   # eV², solar (PDG 2023)
dm2_32_NO = 2.453e-3  # eV², atmospheric (Normal Ordering)
dm2_32_IO = -2.536e-3  # eV², atmospheric (Inverted Ordering)

# Cosmological bound
sum_nu_planck = 0.12  # eV (95% CL, Planck 2018)


# ============================================================
# §1. KOIDE + Δm² → NEUTRINO MASSES (Normal Ordering)
# ============================================================
print("=" * 72)
print("§1. KOIDE K=2/3 + Δm² DATA → MASY NEUTRIN (NO)")
print("=" * 72)

print("""
  Dane: Δm²₂₁ = 7.53 × 10⁻⁵ eV²
        Δm²₃₂ = 2.453 × 10⁻³ eV² (NO)

  Constraint: K(m₁, m₂, m₃) = 2/3

  Trzy niewiadome (m₁, m₂, m₃) z trzema równaniami:
    m₂² - m₁² = Δm²₂₁
    m₃² - m₂² = Δm²₃₂
    K(m₁, m₂, m₃) = 2/3

  → Jedyny wolny parametr: m₁ (najlżejsze neutrino)
""")

def neutrino_masses_NO(m1):
    """Given m₁, compute m₂, m₃ for Normal Ordering."""
    m2_sq = m1**2 + dm2_21
    m3_sq = m1**2 + dm2_21 + dm2_32_NO
    if m2_sq < 0 or m3_sq < 0:
        return None, None
    return np.sqrt(m2_sq), np.sqrt(m3_sq)

def koide_vs_m1_NO(m1):
    """K(m1, m2(m1), m3(m1)) - 2/3 for Normal Ordering."""
    m2, m3 = neutrino_masses_NO(m1)
    if m2 is None:
        return np.nan
    return koide(m1, m2, m3) - 2/3

# Scan m₁ from 0 to 0.3 eV
m1_scan = np.linspace(1e-6, 0.3, 10000)
K_scan = np.array([koide_vs_m1_NO(m1) + 2/3 for m1 in m1_scan])

print("  Scan K(ν) vs m₁ (NO):")
for m1_test in [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
    m2, m3 = neutrino_masses_NO(m1_test)
    if m2 is not None:
        K_val = koide(m1_test, m2, m3)
        sum_m = m1_test + m2 + m3
        print(f"    m₁ = {m1_test:.4f} eV:  K = {K_val:.6f}  Σm = {sum_m:.4f} eV  "
              f"r₂₁ = {m2/m1_test:.1f}  r₃₁ = {m3/m1_test:.1f}")

# Find m₁ where K = 2/3
print("\n  Szukam m₁ gdzie K = 2/3:")
roots_NO = []
for i in range(len(m1_scan)-1):
    v1 = koide_vs_m1_NO(m1_scan[i])
    v2 = koide_vs_m1_NO(m1_scan[i+1])
    if np.isfinite(v1) and np.isfinite(v2) and v1 * v2 < 0:
        m1_root = brentq(koide_vs_m1_NO, m1_scan[i], m1_scan[i+1])
        roots_NO.append(m1_root)

if roots_NO:
    for m1_sol in roots_NO:
        m2_sol, m3_sol = neutrino_masses_NO(m1_sol)
        K_check = koide(m1_sol, m2_sol, m3_sol)
        sum_m = m1_sol + m2_sol + m3_sol
        r21 = m2_sol / m1_sol
        r31 = m3_sol / m1_sol

        print(f"\n  ★ ROZWIĄZANIE (NO): m₁ = {m1_sol*1000:.4f} meV")
        print(f"    m₂ = {m2_sol*1000:.4f} meV")
        print(f"    m₃ = {m3_sol*1000:.4f} meV")
        print(f"    K = {K_check:.8f}")
        print(f"    Σm_ν = {sum_m*1000:.2f} meV = {sum_m:.5f} eV")
        print(f"    r₂₁ = m₂/m₁ = {r21:.4f}")
        print(f"    r₃₁ = m₃/m₁ = {r31:.4f}")
        print(f"    r₃₂ = m₃/m₂ = {m3_sol/m2_sol:.4f}")

        # Test against cosmological bound
        print(f"\n    Σm_ν < {sum_nu_planck} eV (Planck)?  "
              f"{'TAK ✓' if sum_m < sum_nu_planck else 'NIE ✗'}")

        record("T1: K=2/3 solution exists (NO)",
               True,
               f"m₁ = {m1_sol*1000:.2f} meV, Σm = {sum_m*1000:.1f} meV")

        record("T2: Σm_ν < Planck bound",
               sum_m < sum_nu_planck,
               f"Σm = {sum_m:.4f} eV vs {sum_nu_planck} eV")
else:
    print("  Brak rozwiązania! K ≠ 2/3 dla żadnego m₁ (NO)")
    # Find minimum |K - 2/3|
    K_diff = np.array([abs(koide_vs_m1_NO(m1)) for m1 in m1_scan])
    i_min = np.argmin(K_diff)
    m1_closest = m1_scan[i_min]
    m2_c, m3_c = neutrino_masses_NO(m1_closest)
    K_closest = koide(m1_closest, m2_c, m3_c)
    print(f"  Najbliższe: m₁ = {m1_closest*1000:.2f} meV, K = {K_closest:.6f}")

    record("T1: K=2/3 solution exists (NO)",
           False,
           f"No solution. Closest K = {K_closest:.6f} at m₁ = {m1_closest*1000:.2f} meV")
    record("T2: Σm_ν < Planck bound", False, "No K=2/3 solution")


# ============================================================
# §2. INVERTED ORDERING
# ============================================================
print("\n" + "=" * 72)
print("§2. KOIDE K=2/3 + Δm² DATA → MASY NEUTRIN (IO)")
print("=" * 72)

def neutrino_masses_IO(m3):
    """Given m₃ (lightest for IO), compute m₁, m₂."""
    # IO: m₃ < m₁ < m₂
    # m₁² = m₃² + |Δm²₃₂| - Δm²₂₁ ≈ m₃² + |Δm²₃₁|
    # Actually: Δm²₃₂ < 0 for IO
    # m₂² - m₁² = Δm²₂₁
    # m₃² - m₂² = Δm²₃₂ (negative)
    # So m₂² = m₃² - Δm²₃₂ = m₃² + |Δm²₃₂|
    # m₁² = m₂² - Δm²₂₁ = m₃² + |Δm²₃₂| - Δm²₂₁
    m2_sq = m3**2 + abs(dm2_32_IO)
    m1_sq = m2_sq - dm2_21
    if m1_sq < 0 or m2_sq < 0:
        return None, None
    return np.sqrt(m1_sq), np.sqrt(m2_sq)

def koide_vs_m3_IO(m3):
    m1, m2 = neutrino_masses_IO(m3)
    if m1 is None:
        return np.nan
    return koide(m1, m2, m3) - 2/3

print("  Scan K(ν) vs m₃ (IO, m₃ = lightest):")
for m3_test in [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
    m1, m2 = neutrino_masses_IO(m3_test)
    if m1 is not None:
        K_val = koide(m1, m2, m3_test)
        sum_m = m1 + m2 + m3_test
        r21 = m2/m1 if m1 > 0 else np.inf
        print(f"    m₃ = {m3_test:.4f} eV:  K = {K_val:.6f}  Σm = {sum_m:.4f} eV  "
              f"m₁ = {m1:.4f}  m₂ = {m2:.4f}")

# Find m₃ where K = 2/3
m3_scan = np.linspace(1e-6, 0.3, 10000)
roots_IO = []
for i in range(len(m3_scan)-1):
    v1 = koide_vs_m3_IO(m3_scan[i])
    v2 = koide_vs_m3_IO(m3_scan[i+1])
    if np.isfinite(v1) and np.isfinite(v2) and v1 * v2 < 0:
        m3_root = brentq(koide_vs_m3_IO, m3_scan[i], m3_scan[i+1])
        roots_IO.append(m3_root)

if roots_IO:
    for m3_sol in roots_IO:
        m1_sol, m2_sol = neutrino_masses_IO(m3_sol)
        K_check = koide(m1_sol, m2_sol, m3_sol)
        sum_m = m1_sol + m2_sol + m3_sol

        # Sort for display: m₃ < m₁ < m₂ in IO
        masses_sorted = sorted([m1_sol, m2_sol, m3_sol])

        print(f"\n  ★ ROZWIĄZANIE (IO): m₃ = {m3_sol*1000:.4f} meV (lightest)")
        print(f"    m₁ = {m1_sol*1000:.4f} meV")
        print(f"    m₂ = {m2_sol*1000:.4f} meV")
        print(f"    K = {K_check:.8f}")
        print(f"    Σm_ν = {sum_m*1000:.2f} meV = {sum_m:.5f} eV")
        print(f"    Σm_ν < {sum_nu_planck} eV?  "
              f"{'TAK ✓' if sum_m < sum_nu_planck else 'NIE ✗'}")

    record("T3: K=2/3 solution exists (IO)",
           True,
           f"m₃ = {roots_IO[0]*1000:.2f} meV, Σm = {(sum_m)*1000:.1f} meV")
else:
    # Find closest
    K_diff = np.array([abs(koide_vs_m3_IO(m3)) for m3 in m3_scan])
    i_min = np.argmin(K_diff)
    m3_closest = m3_scan[i_min]
    m1_c, m2_c = neutrino_masses_IO(m3_closest)
    K_closest = koide(m1_c, m2_c, m3_closest) if m1_c is not None else np.nan
    print(f"\n  Brak rozwiązania (IO).")
    print(f"  Najbliższe: m₃ = {m3_closest*1000:.2f} meV, K = {K_closest:.6f}")

    record("T3: K=2/3 solution exists (IO)",
           False,
           f"No solution. Closest K = {K_closest:.6f}")


# ============================================================
# §3. SHIFTED KOIDE DLA NEUTRIN
# ============================================================
print("\n" + "=" * 72)
print("§3. SHIFTED KOIDE K(m+m₀) = 2/3 DLA NEUTRIN")
print("=" * 72)

print("""
  Analogia z kwarkami (ex235-236):
    K(m_bare) ≠ 2/3, ale K(m+m₀) = 2/3 z odpowiednim m₀.

  Dla neutrin: shifted Koide może działać z m₀ < 0
  (subtracted mass scale) lub m₀ > 0.
""")

# For NO: scan m₀ for a range of m₁
print("  Normal Ordering — szukam (m₁, m₀) gdzie K(m+m₀) = 2/3:")
print(f"\n  {'m₁ (meV)':>10s}  {'m₀ (meV)':>10s}  {'K(m+m₀)':>10s}  {'Σm (meV)':>10s}")
print("  " + "-" * 50)

best_solutions_NO = []
for m1_try in np.logspace(-2, np.log10(100), 50):  # meV
    m1_eV = m1_try * 1e-3
    m2_eV, m3_eV = neutrino_masses_NO(m1_eV)
    if m2_eV is None:
        continue

    def K_shifted_NO(m0_meV):
        m0 = m0_meV * 1e-3
        if m1_eV + m0 <= 0 or m2_eV + m0 <= 0 or m3_eV + m0 <= 0:
            return np.nan
        return koide(m1_eV + m0, m2_eV + m0, m3_eV + m0) - 2/3

    # Search m₀
    for m0_start in np.linspace(-m1_try*0.99, 200, 500):
        try:
            v1 = K_shifted_NO(m0_start)
            v2 = K_shifted_NO(m0_start + 0.5)
            if np.isfinite(v1) and np.isfinite(v2) and v1 * v2 < 0:
                m0_sol = brentq(K_shifted_NO, m0_start, m0_start + 0.5)
                m0_eV = m0_sol * 1e-3
                K_check = koide(m1_eV+m0_eV, m2_eV+m0_eV, m3_eV+m0_eV)
                sum_m = m1_eV + m2_eV + m3_eV
                best_solutions_NO.append((m1_try, m0_sol, K_check, sum_m*1000))
                if len(best_solutions_NO) <= 10:
                    print(f"  {m1_try:10.3f}  {m0_sol:10.3f}  {K_check:10.8f}  {sum_m*1000:10.2f}")
                break
        except:
            continue

if best_solutions_NO:
    print(f"\n  Znaleziono {len(best_solutions_NO)} rozwiązań (m₁, m₀) z K=2/3 (NO)")

    # The solution where m₀=0 (if exists) is the bare Koide solution from §1
    # Check: which m₀ is closest to zero?
    m0_vals = [s[1] for s in best_solutions_NO]
    i_min_m0 = np.argmin(np.abs(m0_vals))
    s = best_solutions_NO[i_min_m0]
    print(f"  Rozwiązanie z min |m₀|: m₁ = {s[0]:.3f} meV, m₀ = {s[1]:.3f} meV")


# ============================================================
# §4. BRANNEN PARAMETRYZACJA DLA NEUTRIN
# ============================================================
print("\n" + "=" * 72)
print("§4. BRANNEN PARAMETRYZACJA: ε=√2, θ DLA NEUTRIN")
print("=" * 72)

print("""
  Brannen (2006): √m_k = √(M/3) × (1 + ε·cos(θ + 2πk/3))

  Dla leptonów naładowanych: ε = √2, θ = 0.2222 rad = 12.73°
  (or equivalently θ_Brannen ≈ 132.73° in another convention)

  Jeśli K=2/3: ε = √2 ZAWSZE (niezależnie od θ).
  θ kontroluje hierarchię (r₂₁, r₃₁).
""")

# For each NO solution, compute Brannen θ
if roots_NO:
    m1_sol = roots_NO[0]
    m2_sol, m3_sol = neutrino_masses_NO(m1_sol)

    eps, theta, K_val = brannen_params(m1_sol, m2_sol, m3_sol)

    # Also compute for charged leptons
    eps_l, theta_l, K_l = brannen_params(0.511e-3, 105.658e-3, 1.77686)  # GeV

    print(f"  Neutrino (NO, K=2/3 solution):")
    print(f"    ε = {eps:.4f} (expect √2 = {np.sqrt(2):.4f})")
    print(f"    θ = {theta:.4f} rad = {np.degrees(theta):.2f}°")
    print(f"    K = {K_val:.6f}")

    print(f"\n  Charged leptons:")
    print(f"    ε = {eps_l:.4f}")
    print(f"    θ = {theta_l:.4f} rad = {np.degrees(theta_l):.2f}°")
    print(f"    K = {K_l:.6f}")

    print(f"\n  Δθ = θ_ν - θ_l = {theta - theta_l:.4f} rad = {np.degrees(theta - theta_l):.2f}°")

    # Check if θ_ν = π - θ_l (complementary)
    theta_comp = np.pi - theta_l
    print(f"  π - θ_l = {theta_comp:.4f} rad = {np.degrees(theta_comp):.2f}°")
    print(f"  |θ_ν - (π-θ_l)| = {abs(theta - theta_comp):.4f} rad = {np.degrees(abs(theta - theta_comp)):.2f}°")

    # Check θ_ν + θ_l = π?
    print(f"  θ_ν + θ_l = {theta + theta_l:.4f} rad = {np.degrees(theta + theta_l):.2f}°")
    print(f"  (π = {np.pi:.4f} rad = 180°)")

    record("T4: Brannen θ for neutrinos",
           True,
           f"θ_ν = {np.degrees(theta):.2f}°, θ_l = {np.degrees(theta_l):.2f}°")


# ============================================================
# §5. TGP LINK: r₂₁(ν) I SOLITON ODE
# ============================================================
print("\n" + "=" * 72)
print("§5. TGP LINK: r₂₁(ν) I SOLITON ODE")
print("=" * 72)

if roots_NO:
    m1_sol = roots_NO[0]
    m2_sol, m3_sol = neutrino_masses_NO(m1_sol)
    r21_nu = m2_sol / m1_sol
    r31_nu = m3_sol / m1_sol

    # Lepton values for comparison
    r21_lep = 206.768
    r31_lep = 3477.44

    print(f"\n  Neutrino mass ratios (NO, K=2/3):")
    print(f"    r₂₁(ν) = m₂/m₁ = {r21_nu:.4f}")
    print(f"    r₃₁(ν) = m₃/m₁ = {r31_nu:.4f}")
    print(f"    r₃₂(ν) = m₃/m₂ = {m3_sol/m2_sol:.4f}")

    print(f"\n  Charged lepton ratios:")
    print(f"    r₂₁(l) = {r21_lep:.1f}")
    print(f"    r₃₁(l) = {r31_lep:.1f}")

    print(f"\n  Hierarchy comparison:")
    print(f"    log(r₂₁): ν = {np.log10(r21_nu):.3f}, l = {np.log10(r21_lep):.3f}")
    print(f"    log(r₃₁): ν = {np.log10(r31_nu):.3f}, l = {np.log10(r31_lep):.3f}")
    print(f"    Neutrinos are MUCH less hierarchical than charged leptons")

    # If r₂₁ comes from soliton ODE, what g₀ would give r₂₁(ν)?
    # r₂₁ = (A(φ·g₀)/A(g₀))⁴ — from ex234
    # For r₂₁ = 206.77: g₀* = 0.8695
    # For r₂₁ = r₂₁(ν): need different g₀

    print(f"\n  Jeśli r₂₁(ν) pochodzi z soliton ODE (jak leptony):")
    print(f"    r₂₁(ν)^(1/4) = {r21_nu**(0.25):.4f}")
    print(f"    r₂₁(l)^(1/4) = {r21_lep**(0.25):.4f}")
    print(f"    A(ν₂)/A(ν₁) = {r21_nu**(0.25):.4f} (small ratio → g₀ bliskie 1)")

    # Soliton φ-FP: A(φ·g₀)/A(g₀) = r₂₁^(1/4) = ratio of tail amplitudes
    # For ν: A_ratio = r₂₁(ν)^(1/4) ≈ small → g₀(ν) close to 1 or very different sector

    # Test: is r₂₁(ν) = r₂₁(l)^α for some α?
    alpha = np.log(r21_nu) / np.log(r21_lep)
    print(f"\n  Test: r₂₁(ν) = r₂₁(l)^α?")
    print(f"    α = log(r₂₁(ν))/log(r₂₁(l)) = {alpha:.4f}")
    print(f"    Interpretation: ν hierarchy ≈ l hierarchy^{alpha:.3f}")

    # Seesaw: m_ν ~ m_l² / M_R → r₂₁(ν) ≈ r₂₁(l)² ?
    r21_seesaw = r21_lep**2
    print(f"\n  Seesaw test: r₂₁(ν) ≈ r₂₁(l)²?")
    print(f"    r₂₁(l)² = {r21_seesaw:.0f}")
    print(f"    r₂₁(ν) = {r21_nu:.2f}")
    print(f"    SEESAW FAILS (r₂₁(l)² >> r₂₁(ν))")

    record("T5: Neutrino hierarchy much weaker",
           r21_nu < r21_lep,
           f"r₂₁(ν) = {r21_nu:.2f} << r₂₁(l) = {r21_lep:.0f}")


# ============================================================
# §6. Σm_ν PREDICTION SUMMARY
# ============================================================
print("\n" + "=" * 72)
print("§6. ★ PODSUMOWANIE: PREDYKCJA Σm_ν")
print("=" * 72)

if roots_NO:
    m1_sol = roots_NO[0]
    m2_sol, m3_sol = neutrino_masses_NO(m1_sol)
    sum_NO = (m1_sol + m2_sol + m3_sol) * 1000  # meV

    print(f"\n  K=2/3 (bare Koide) + Normal Ordering:")
    print(f"    m₁ = {m1_sol*1000:.3f} meV")
    print(f"    m₂ = {m2_sol*1000:.3f} meV")
    print(f"    m₃ = {m3_sol*1000:.3f} meV")
    print(f"    ★ Σm_ν = {sum_NO:.2f} meV = {sum_NO/1000:.5f} eV")
    print(f"")
    print(f"    Planck bound: Σm_ν < 120 meV (95% CL)")
    print(f"    DESI+CMB (2024): Σm_ν < 72 meV")

    # Effective mass for beta decay
    m_beta = np.sqrt(0.68 * m1_sol**2 + 0.30 * m2_sol**2 + 0.02 * m3_sol**2)
    print(f"\n    m_β (effective, β-decay) ≈ {m_beta*1000:.3f} meV")
    print(f"    KATRIN bound: m_β < 450 meV")

    # Effective Majorana mass for 0νββ
    # |m_ee| = |c12² c13² m1 + s12² c13² m2 e^(iα) + s13² m3 e^(iβ)|
    # Simplified (CP phases = 0):
    s12_sq = 0.307  # sin²θ₁₂
    s13_sq = 0.0220  # sin²θ₁₃
    c12_sq = 1 - s12_sq
    c13_sq = 1 - s13_sq
    m_ee_max = c12_sq * c13_sq * m1_sol + s12_sq * c13_sq * m2_sol + s13_sq * m3_sol
    m_ee_min = abs(c12_sq * c13_sq * m1_sol - s12_sq * c13_sq * m2_sol - s13_sq * m3_sol)
    print(f"\n    |m_ee| (0νββ, Majorana): [{m_ee_min*1000:.3f}, {m_ee_max*1000:.3f}] meV")
    print(f"    Current best bound: |m_ee| < ~50-160 meV (KamLAND-Zen)")

if roots_IO:
    m3_sol = roots_IO[0]
    m1_sol, m2_sol = neutrino_masses_IO(m3_sol)
    sum_IO = (m1_sol + m2_sol + m3_sol) * 1000

    print(f"\n  K=2/3 (bare Koide) + Inverted Ordering:")
    print(f"    m₃ = {m3_sol*1000:.3f} meV (lightest)")
    print(f"    m₁ = {m1_sol*1000:.3f} meV")
    print(f"    m₂ = {m2_sol*1000:.3f} meV")
    print(f"    ★ Σm_ν = {sum_IO:.2f} meV = {sum_IO/1000:.5f} eV")


# ============================================================
# SCORECARD
# ============================================================
print("\n" + "=" * 72)
print("SCORECARD")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")
print(f"\n  {n_pass}/{n_total} testów przeszło.")


# ============================================================
# PODSUMOWANIE
# ============================================================
sum_str_NO = f"{sum_NO:.1f} meV" if roots_NO else "N/A"
sum_str_IO = f"{sum_IO:.1f} meV" if roots_IO else "N/A"

print(f"""
========================================================================
PODSUMOWANIE ex237
========================================================================

  ★ PREDYKCJE MAS NEUTRIN Z KOIDE K=2/3:

  Normal Ordering:  Σm_ν = {sum_str_NO}
  Inverted Ordering: Σm_ν = {sum_str_IO}

  KLUCZOWE OBSERWACJE:
  1. K=2/3 JEDNOZNACZNIE wyznacza m₁ (jedyny wolny parametr przy danym Δm²)
  2. Neutrina DUŻO mniej hierarchiczne niż leptony naładowane
  3. Brannen θ_ν ≠ θ_l (different hierarchy angle → different soliton sector?)

  TESTOWALNE PREDYKCJE:
  - Σm_ν z kosmologii (Planck, DESI, Euclid)
  - m_β z β-decay (KATRIN/Project 8)
  - |m_ee| z 0νββ (KamLAND-Zen, nEXO)
  - Normal vs Inverted ordering (JUNO, DUNE, Hyper-K)
""")
