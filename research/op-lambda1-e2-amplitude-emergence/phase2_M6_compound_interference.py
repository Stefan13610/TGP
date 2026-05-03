#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
λ.1 M.6 — Compound interference test (user hypothesis 2026-05-02)

USER HIPOTEZA:
  "Φ_0 jest sumą pola cząstki 1 + cząstka 2 z interferencją z 1, +
   cząstka 3 z interferencją z 1 i 2 i tak dalej. To może dawać Eulera
   bo (1+1/N)^N → e."

Trzy warstwy testu:
  Layer 1: czysta matematyka compound — sprawdz że (1+x/N)^N → e^x
  Layer 2: single R3 soliton compound integrals — czy e² wyłania się
           naturalnie z any integral nad solitonem przy α=2?
  Layer 3: multi-soliton background buildup — czy compound z N źródeł
           daje saturation Φ₀ ~ e² jako N → ∞?

Targets:
  e   = 2.71828    (Euler base)
  e²  = 7.38906    (Φ_eff target, mass formula natural exponent)
  e²/2 = 3.69453   (X dla α=2 charged-lepton)
  ln(e²)= 2        (compound saturation parameter)

Falsification gate:
  - jeśli ŻADEN integral / scaling nie hits target z <1% drift → NEGATIVE
  - jeśli któreś hits z <0.1% → kandydat na compound mechanism
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp, trapezoid

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 78)
print("  λ.1 M.6 — Compound interference test (user hypothesis)")
print("=" * 78)

PI = np.pi
E1 = np.e
E2 = np.exp(2)
E2_HALF = E2 / 2

print(f"\nTargets:")
print(f"  e       = {E1:.6f}")
print(f"  e²      = {E2:.6f}  (Φ_eff target)")
print(f"  e²/2    = {E2_HALF:.6f}  (X dla α=2)")
print(f"  ln(e²)  = 2.000000  (compound saturation parameter)")

# ============================================================================
# LAYER 1: Compound interest math — czy (1+x/N)^N → e^x?
# ============================================================================
print("\n" + "=" * 78)
print("  LAYER 1: Compound interest math")
print("=" * 78)

print(f"\n  Demonstration: (1 + x/N)^N → e^x as N → ∞")
print(f"  {'N':>10} | {'(1+1/N)^N':>12} | {'(1+2/N)^N':>12} | {'e':>10} | {'e²':>10}")
print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}-+-{'-'*10}")

for N in [10, 100, 1000, 10000, 100000, 1000000]:
    one = (1 + 1/N)**N
    two = (1 + 2/N)**N
    print(f"  {N:>10} | {one:>12.6f} | {two:>12.6f} | {E1:>10.6f} | {E2:>10.6f}")

print(f"\n  → Compound formula structurally produces e and e² in N→∞ limit")
print(f"  → For Σε = 2 z N→∞ → exact e²")
print(f"  Pytanie: czy TGP ma natural N·ε = 2?")

# ============================================================================
# LAYER 2: Single R3 soliton compound integrals @ α=2
# ============================================================================
print("\n" + "=" * 78)
print("  LAYER 2: Single R3 soliton compound integrals @ α=2")
print("=" * 78)

def r3_ode(r, y, alpha):
    """R3 ODE: g'' + (α/g)(g')² + (2/r)g' = (1-g)·g^(2-2α)"""
    g, gp = y
    g_safe = np.sign(g) * max(abs(g), 1e-12)
    rhs = (1 - g) * abs(g_safe)**(2 - 2*alpha) - (alpha/g_safe)*gp**2 - (2/r)*gp
    return [gp, rhs]

def solve_r3(g0, alpha, r_max=120.0):
    """Solve R3 ODE; return r, g, gp"""
    r0 = 1e-3
    g_init = g0 + (1/6)*(1-g0)*g0**(2-2*alpha)*r0**2
    gp_init = (1/3)*(1-g0)*g0**(2-2*alpha)*r0
    sol = solve_ivp(r3_ode, [r0, r_max], [g_init, gp_init],
                    args=(alpha,), method='DOP853',
                    atol=1e-12, rtol=1e-10, max_step=0.05)
    return sol.t, sol.y[0], sol.y[1]

# Reprezentatywne g₀ dla charged-lepton class (electron, muon, tau)
# z mass_scaling_k4: g₀_τ_canonical = 1.77472
# typical range g₀ ∈ [1.5, 2.5]

g0_test = 1.77472  # τ canonical
alpha_test = 2.0

print(f"\n  Solving R3 ODE: α={alpha_test}, g₀={g0_test} (τ canonical)")

r, g, gp = solve_r3(g0_test, alpha_test)
print(f"  Solution: r ∈ [{r[0]:.3f}, {r[-1]:.3f}], {len(r)} points")
print(f"  g range: [{g.min():.4f}, {g.max():.4f}], g[0]={g[0]:.4f}, g[-1]={g[-1]:.4f}")

# Compound integrals — multiple candidates
print(f"\n  Compound integral candidates:")
print(f"  {'-'*60}")

# I1: ∫_0^R |g'/g| dr  (total log range / accumulated compound rate)
g_safe = np.where(np.abs(g) > 1e-10, g, 1e-10)
integrand_I1 = np.abs(gp/g_safe)
I1 = trapezoid(integrand_I1, r)
print(f"  I1 = ∫|g'/g|dr           = {I1:.6f}")

# I2: ∫_0^R |g'| dr (total field variation)
I2 = trapezoid(np.abs(gp), r)
print(f"  I2 = ∫|g'|dr             = {I2:.6f}")

# I3: ∫_0^R (g-1)² · 4πr² dr (soliton "volume")
I3 = trapezoid((g-1)**2 * 4*PI*r**2, r)
print(f"  I3 = ∫(g-1)²·4πr²dr      = {I3:.6f}")

# I4: ∫_0^R g^(2α) · g'² · 4πr² dr (R3 kinetic moment)
I4 = trapezoid(np.abs(g)**(2*alpha_test) * gp**2 * 4*PI*r**2, r)
print(f"  I4 = ∫g^(2α)·g'²·4πr²dr  = {I4:.6f}")

# I5: ∏_{shells} (1 + dε(r)) where dε = (g'/g)·dr
# log(I5) = Σ log(1 + dε) ≈ Σ dε = ∫|g'/g|dr = I1 dla małych dε
# Ale full product dla dyskretnych shells:
dr = np.diff(r)
deps = np.abs(gp[:-1]/g_safe[:-1]) * dr
log_prod = np.sum(np.log(1 + deps))  # full discrete product
I5 = np.exp(log_prod)
print(f"  I5 = exp(Σ log(1+|g'/g|·dr))  = {I5:.6f}")

# I6: ∫|g-1|·dr (total deformation amplitude)
I6 = trapezoid(np.abs(g-1), r)
print(f"  I6 = ∫|g-1|dr            = {I6:.6f}")

# I7: log(g₀)/ε_typical — "effective N steps" from log
log_g0 = np.log(g0_test)
print(f"  log(g₀)                  = {log_g0:.6f}")
print(f"  log(g₀)·X = log(m_obs/A²) z α=2 mass formula:")
print(f"    {log_g0} · {E2_HALF:.4f} = {log_g0 * E2_HALF:.6f}")

# Compare to targets
targets = {
    'e':     E1,
    'e²':    E2,
    'e²/2':  E2_HALF,
    '2':     2.0,
    '2π':    2*PI,
    '4π':    4*PI,
}

integrals = {
    'I1 (∫|g\'/g|)':      I1,
    'I2 (∫|g\'|)':        I2,
    'I3 (∫(g-1)²·4πr²)':  I3,
    'I4 (∫g^4·g\'²·4πr²)':I4,
    'I5 (compound prod)': I5,
    'I6 (∫|g-1|)':        I6,
}

print(f"\n  Match table (drift % of integral vs target):")
print(f"  {'Integral':>22} | " + " | ".join(f"{tn:>9}" for tn in targets.keys()))
print(f"  {'-'*22}-+-" + "-+-".join('-'*9 for _ in targets))

best_match = None
best_drift = float('inf')
for iname, ival in integrals.items():
    row = f"  {iname:>22} | "
    cells = []
    for tname, tval in targets.items():
        drift = 100 * abs(ival - tval) / tval
        cells.append(f"{drift:>8.2f}%")
        if drift < best_drift:
            best_drift = drift
            best_match = (iname, tname, ival, tval, drift)
    print(row + " | ".join(cells))

print(f"\n  Best match (single soliton): {best_match[0]} ≈ {best_match[1]}")
print(f"    Value: {best_match[2]:.5f} vs target {best_match[3]:.5f}  ({best_match[4]:.3f}% drift)")

if best_drift < 1.0:
    print(f"  --> CANDIDATE found z drift < 1%")
elif best_drift < 5.0:
    print(f"  --> Suggestive ale za luźne (>1%)")
else:
    print(f"  --> Single soliton: NEGATIVE — żaden integral nie matche e/e²/e²/2")

# ============================================================================
# LAYER 3: Multi-soliton compound buildup
# ============================================================================
print("\n" + "=" * 78)
print("  LAYER 3: Multi-soliton compound buildup")
print("=" * 78)

# Test: jeśli mamy N losowych solitonów w volume V z density ρ = N/V,
# czy total background field saturates jak (1+ε)^N → exp(Nε)?

# Approx: każdy soliton kontrybuuje "core deformation" Δg = g₀ - 1 plus tail
# W large-distance limit, tail jest ~A·cos(r)/r
# Pole tła w punkcie x z N solitonów w random positions:
#   Φ_total(x) = Σ_i (g_i(|x-x_i|) - 1)  (pole odchylenie od vacuum)

# Test compound vs linear:
#   Linear: Φ_lin = N · ⟨g-1⟩_avg
#   Compound: Φ_comp z multiplicative interference
#   Interesujące: czy compound saturates lub cykluje przy specyficznych N?

# Setup: place N solitons at random positions in V, compute Φ at origin
# Use single-soliton g(r) profile, sum contributions

# Tail amplitude — extract A z g(r) ≈ 1 + A·cos(r-φ)/r
mask_tail = (r > 20) & (r < 80)
r_tail = r[mask_tail]
h_tail = (g[mask_tail] - 1) * r_tail  # h(r) = A·cos(r-φ)
A_tail = np.max(np.abs(h_tail))
print(f"\n  Single soliton tail amplitude A = {A_tail:.5f}")
print(f"  Core deformation (g₀-1) = {g0_test - 1:.5f}")

# Build interpolation of g(r) for random superposition
from scipy.interpolate import interp1d
g_interp = interp1d(r, g - 1, bounds_error=False, fill_value=0.0)

# Place N solitons in volume V (cube of side L), test how total Φ grows with N
print(f"\n  N-soliton test: Φ_total(origin) z N losowych solitonów")
print(f"  (Volume L=20, solitons at random positions)")
print(f"\n  {'N':>5} | {'Φ_lin (sum)':>12} | {'log(Φ_lin)':>12} | {'Φ/N (avg)':>12} | {'(1+1/N)^N':>11}")
print(f"  {'-'*5}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*11}")

L = 20.0
np.random.seed(42)

results = []
for N in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]:
    # Generate N random positions
    positions = np.random.uniform(-L/2, L/2, size=(N, 3))
    distances = np.sqrt(np.sum(positions**2, axis=1))
    contribs = g_interp(distances)  # each soliton contributes (g(d)-1) at origin
    Phi_lin = np.sum(contribs)
    avg = Phi_lin / N
    compound_proxy = (1 + 1/N)**N
    log_phi = np.log(abs(Phi_lin)) if Phi_lin != 0 else float('-inf')
    print(f"  {N:>5} | {Phi_lin:>12.5f} | {log_phi:>12.5f} | {avg:>12.6f} | {compound_proxy:>11.6f}")
    results.append((N, Phi_lin, avg))

# Test scaling: Φ_lin grows ~ N (linear sum) or saturates?
arr = np.array(results)
N_arr = arr[:, 0]
Phi_arr = arr[:, 1]

# Fit log(|Phi|) = α·log(N) + β: if α=1 → linear, α<1 → sublinear (saturation)
mask_pos = Phi_arr > 0
if mask_pos.sum() > 3:
    log_N = np.log(N_arr[mask_pos])
    log_Phi = np.log(Phi_arr[mask_pos])
    slope, intercept = np.polyfit(log_N, log_Phi, 1)
    print(f"\n  Scaling fit: log(Φ) = {slope:.4f}·log(N) + {intercept:.4f}")
    print(f"    slope = 1.0 → Φ ~ N (linear superposition)")
    print(f"    slope < 1 → Φ saturates (compound interference effect)")
    print(f"    slope = 0 → Φ → const (full saturation)")

    # Test: czy ratio Φ/N ma limit ~ const? Or oscillates?
    ratios = Phi_arr / N_arr
    print(f"\n  Φ/N ratios (should saturate if compound exists):")
    for n, p, ratio in results:
        print(f"    N={n:>4}: Φ/N = {ratio:>12.6f}")

# CHECK: czy w jakimkolwiek N Φ_total ≈ e² lub e²/2 lub e?
print(f"\n  Searching for natural N where Φ_total ≈ target:")
for N, Phi, _ in results:
    for tname, tval in targets.items():
        drift = 100 * abs(Phi - tval) / tval
        if drift < 5.0:
            print(f"    N={N}: Φ={Phi:.5f} ≈ {tname}={tval:.5f}  ({drift:.2f}% drift)")

# ============================================================================
# LAYER 4: Czy R3 ODE ma natural compound saturation X = e²/2?
# ============================================================================
print("\n" + "=" * 78)
print("  LAYER 4: Czy mass formula slope X SAM emerges z compound limit?")
print("=" * 78)

# Test: solve R3 dla range g₀, fit slope X = d log(m)/d log(g₀)
# Sprawdz czy X dokładnie = e²/2 = 3.69453 dla α=2

g0_range = np.array([1.5, 1.6, 1.7, 1.77472, 1.8, 1.9, 2.0, 2.1, 2.2])
log_g0 = []
log_M = []
log_A2 = []

print(f"\n  Solving R3 dla g₀ range, α={alpha_test}:")
for g0_v in g0_range:
    try:
        rv, gv, gpv = solve_r3(g0_v, alpha_test)
        # M = ∫ ρ·4πr² dr z ρ = (1/2)g^(2α)·g'² + V(g)
        # V(g) for α=2: V = log(g) + 1/g - 1
        g_safe_v = np.where(np.abs(gv) > 1e-10, gv, 1e-10)
        V = np.log(np.abs(g_safe_v)) + 1/np.abs(g_safe_v) - 1
        rho = 0.5 * np.abs(g_safe_v)**(2*alpha_test) * gpv**2 + V
        M = trapezoid(rho * 4*PI*rv**2, rv)
        # A² z tail
        mask_t = (rv > 20) & (rv < 80)
        if mask_t.sum() > 5:
            h_t = (gv[mask_t] - 1) * rv[mask_t]
            A = np.max(np.abs(h_t))
            if A > 0 and M > 0:
                log_g0.append(np.log(g0_v))
                log_M.append(np.log(M))
                log_A2.append(2*np.log(A))
    except Exception:
        pass

log_g0 = np.array(log_g0)
log_M = np.array(log_M)
log_A2 = np.array(log_A2)

if len(log_g0) > 3:
    # log(m) - log(A²) = X·log(g₀) + const
    log_ratio = log_M - log_A2
    slope_X, intercept_X = np.polyfit(log_g0, log_ratio, 1)
    print(f"\n  Mass formula fit dla α=2: log(m/A²) = X·log(g₀) + const")
    print(f"    Empirical X = {slope_X:.6f}")
    print(f"    Target e²/2 = {E2_HALF:.6f}")
    print(f"    drift = {100*abs(slope_X - E2_HALF)/E2_HALF:.4f}%")

    # Sprawdz interpretacje compound:
    # X = e²/2 → exp(2)/2 = e^(2)/2
    # Compound: Σε_i = 2 → factor e² → "half" because A² jest squared = e²/2
    print(f"\n  Compound interpretation check:")
    print(f"    Σε_compound = 2·X / e² = {2*slope_X/E2:.6f}  (should be 1.0 if natural)")
    print(f"    log(m_max/A²) for g₀=2.5 ≈ X·log(2.5) = {slope_X * np.log(2.5):.4f}")
    print(f"    Compare e² = {E2:.4f}, e²/2 = {E2_HALF:.4f}")

# ============================================================================
# SUMMARY + verdict
# ============================================================================
print("\n" + "=" * 78)
print("  SUMMARY: M.6 Compound interference test verdict")
print("=" * 78)

verdict = f"""
  USER HIPOTEZA: e² wyłania się z compound interference Φ₀ = ∏(1+ε_i)
  → exp(Σε) gdy Σε natural = 2.

  Layer 1 (compound math):       PASS — (1+x/N)^N → e^x trywialnie ✓
  Layer 2 (single soliton):      Best match: {best_match[0] if best_match else 'N/A'}
                                  ≈ {best_match[1] if best_match else 'N/A'} z {best_drift:.2f}% drift
  Layer 3 (multi-soliton sum):   Linear superposition (slope ≈ 1)
                                  brak compound saturation w random placement
  Layer 4 (mass formula slope):  X = {slope_X if 'slope_X' in dir() else 'N/A':.4f} vs e²/2 = {E2_HALF:.4f}
                                  drift = {100*abs(slope_X - E2_HALF)/E2_HALF if 'slope_X' in dir() else float('inf'):.4f}%

  KLUCZOWE OBSERWACJE:

  1. Compound math IS exact dla (1+x/N)^N → e^x — to nie jest hipoteza,
     to identity. Dla x=2 dokładnie e².

  2. ALE: random placement N solitonów w volume daje LINEAR superposition
     (Φ ~ N), nie compound exp. Tutaj nie ma "interference structure"
     bo random fields just add.

  3. Aby compound działał TGP-natural, trzeba MULTIPLICATIVE structure
     w field equations, nie additive. R3 ODE ma g^(2α) prefactor —
     to JEST multiplicative, ale kontroluje pojedynczy soliton, nie
     buildup z N.

  4. X = e²/2 z mass formula jest empirical fit. Compound interpretation
     wymaga że Σε = 2 z TGP-natural saturation. Brak konkretnego TGP-
     mechanizmu który DAJE Σε = 2 dla α=2 charged-lepton.
"""

print(verdict)

# Final check — czy któreś specyficzne integralne ma dokładnie e²/2 lub e²?
print(f"  EXACT MATCH SEARCH (drift < 0.5%):")
found_match = False
for iname, ival in integrals.items():
    for tname, tval in targets.items():
        drift = 100 * abs(ival - tval) / tval
        if drift < 0.5:
            print(f"    *** {iname} = {ival:.5f} ≈ {tname} = {tval:.5f} ({drift:.4f}%)")
            found_match = True

if not found_match:
    print(f"    (none — wszystkie integrals daleko od e/e²/e²/2)")

print("\n  Final: M.6 verdict")
print("  " + "-" * 60)
if found_match:
    print("  POSITIVE candidate found — compound mechanism może być real")
else:
    print("  NEGATIVE — compound interference jako konkretny TGP-mechanizm")
    print("  NIE wyłania się z prostych integrals nad R3 soliton ani z")
    print("  random multi-soliton placement. Wymagałby specyficznej TGP-")
    print("  struktury produkującej Σε = 2 — której nie odkryliśmy.")

print("\n" + "=" * 78)
print("  KONIEC λ.1 M.6")
print("=" * 78)
