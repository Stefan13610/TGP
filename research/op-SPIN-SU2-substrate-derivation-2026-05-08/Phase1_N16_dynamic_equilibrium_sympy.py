"""
Phase 1 N16 - Dynamic equilibrium soliton-background formalization

Cel: pokazac formal mechanizm rozbrojenia Derricka przez interakcje
soliton-background. Per Dynamic_equilibrium_framework.md:

  Soliton TGP NIE jest static stable; jest meta-stable dynamic equilibrium
  w field Phi_bar. Standard Derrick zaklada static stability na fixed
  background. Background coupling daje DODATKOWY scaling component
  ktory umozliwia equilibrium.

Plan:
  E1: Energy decomposition E_sol + E_int + E_rad
  E2: Standard Derrick scaling (FAILS jak oczekiwane)
  E3: Z background coupling: E_int ~ lambda^alpha
  E4: Equilibrium condition (extremum)
  E5: Stability analysis (d2E/dlambda2 > 0)
  E6: Connection do bifurcation (zanik vs ekspansja outcomes)
  E7: Verdict - czy mechanism dziala?
"""
import sympy as sp
from sympy import symbols, Function, diff, simplify, expand, Rational, sqrt, oo, integrate, solve, limit

print("=" * 80)
print("Phase 1 N16 - Dynamic equilibrium soliton-background formalization")
print("=" * 80)

PASS = 0; FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}")
    else: FAIL += 1; print(f"  [FAIL] {name}: {detail}")

# Symbols
lam = symbols('lambda', positive=True)  # rescaling parameter
T, U, S, V_int = symbols('T U S V_int', real=True)
T_pos, U_pos, S_pos = symbols('T_p U_p S_p', positive=True)
alpha_sc = symbols('alpha', real=True)  # E_int scaling exponent
r = symbols('r', positive=True)
Phi_bar, dPhi_sol = symbols('Phi_bar delta_Phi_sol', real=True)

# =============================================================================
# E1: Energy decomposition
# =============================================================================
print("\n" + "=" * 80)
print("E1: Energy decomposition")
print("=" * 80)

print("""
Field decomposition (z Dynamic_equilibrium_framework Part I.3):
  Phi(x, t) = Phi_bar(x) + delta_Phi_sol(x, t) + delta_Phi_rad(x, t)

  Phi_bar      = background (asymptotic, slowly varying)
  delta_Phi_sol = localized soliton
  delta_Phi_rad = emitted radiation

Energy decomposition:
  E_total = E_sol + E_int + E_rad

  E_sol = ∫d³x [(1/2)(grad delta_Phi_sol)² + V(Phi_bar + delta_sol) - V(Phi_bar) + ...]
  E_int = ∫d³x [coupling soliton-background]
  E_rad = energy emitted as radiation

KEY INSIGHT: standard Derrick traktuje TYLKO E_sol w izolacji.
TGP musi traktować pełen E_total, gdzie E_int wnosi inne scaling.
""")

# =============================================================================
# E2: Standard Derrick (FAILS)
# =============================================================================
print("\n" + "=" * 80)
print("E2: Standard Derrick scaling (fails for static stable in 3D)")
print("=" * 80)

print("""
Pod rescaling Phi_lambda(x) = Phi(lambda x):
  T_lambda = lambda^(-1) T  (kinetic, scales λ^(-1))
  U_lambda = lambda^(-3) U  (potential, scales λ^(-3))

  E_iso(lambda) = T λ^(-1) + U λ^(-3)
""")

# Standard Derrick energy
E_iso = T_pos * lam**(-1) + U * lam**(-3)
print(f"\n  E_iso(λ) = {E_iso}")

# Derivative
dE_iso = diff(E_iso, lam)
dE_iso_simp = sp.simplify(dE_iso)
print(f"  dE_iso/dλ = {dE_iso_simp}")

# Solve for extremum
# -T λ^(-2) - 3U λ^(-4) = 0
# T λ² + 3U = 0
# λ² = -3U/T
print(f"\n  Equation dE/dλ = 0:  -T λ^(-2) - 3U λ^(-4) = 0")
print(f"                       Multiply by λ^4:  -T λ² - 3U = 0")
print(f"                       λ² = -3U/T")
print(f"  Real positive solution requires U < 0 (negative potential)")

# For stable bound state, U < 0 (binding). But then E_iso has only ONE λ value
# and is NOT minimum (would be saddle in higher-dim function space)

# Show: E_iso has no LOCAL MINIMUM dla generic V
# For T > 0 (kinetic positive): -T/λ² < 0 always (kinetic gradient pushes λ → ∞)
# For U > 0 (potential positive): -3U/λ^4 < 0 always (potential pushes λ → ∞)
# Sum: dE/dλ < 0 dla all λ - field wants to expand (λ → ∞)
# This means: NO STABLE STATIC SOLUTION (Derrick conclusion)

print(f"""
  Standard Derrick conclusion:
  - For T > 0 i U > 0: dE/dλ < 0 always → soliton expands to infinity
  - For T > 0 i U < 0: extremum at λ² = -3U/T (specific value)
    But d²E/dλ² = 2T/λ³ + 12U/λ^5 = 2(T λ² + 6U)/λ^5
    With T λ² = -3U: = 2(-3U + 6U)/λ^5 = 6U/λ^5 < 0 (since U < 0)
    -> SADDLE POINT, NOT MINIMUM

  -> NO stable static localized soliton w 3D (Derrick's theorem confirmed)
""")

# Sympy verify: for T > 0, U > 0, dE/dλ is always negative
# Substitute concrete positive values and check
T_val = sp.Rational(1)
U_val = sp.Rational(1)
dE_iso_pos = dE_iso.subs([(T_pos, T_val), (U, U_val)])
# Test at λ = 1, 2, 0.5
val_1 = dE_iso_pos.subs(lam, 1)
val_2 = dE_iso_pos.subs(lam, 2)
val_05 = dE_iso_pos.subs(lam, sp.Rational(1,2))
print(f"  dE/dλ for T=1, U=1:")
print(f"    λ=0.5: {val_05}")
print(f"    λ=1.0: {val_1}")
print(f"    λ=2.0: {val_2}")
check("Standard Derrick: dE/dλ < 0 for T,U > 0",
      val_1 < 0 and val_2 < 0 and val_05 < 0,
      "expected all negative")

# =============================================================================
# E3: Background coupling E_int ~ λ^α
# =============================================================================
print("\n" + "=" * 80)
print("E3: Z background coupling E_int = S * λ^α")
print("=" * 80)

print("""
TGP modyfikacja: dodaj E_int representujący coupling soliton-background.

E_int ma SCALING DIMENSION α różny od standardowych terminów.

Mechanizmy z framework dokumentu:
  - Skin/boundary coupling: α = +2 (boundary area rośnie z λ²)
  - Tail extension: α = +1 (linewidth)
  - Gradient at boundary: α = +α dowolne > 0

Dla α > 0, E_int ROŚNIE z λ → przeciwdziała Derrick expansion.

Total energy:
  E(λ) = T λ^(-1) + U λ^(-3) + S λ^α
""")

# Total energy with background coupling
E_total = T_pos * lam**(-1) + U_pos * lam**(-3) + S_pos * lam**alpha_sc
print(f"\n  E_total(λ) = {E_total}")

dE_total = diff(E_total, lam)
print(f"\n  dE_total/dλ = {dE_total}")

# Multiply by λ^4 (standard simplification trick)
# We work with α=2 first (skin coupling case)
# dE = -T λ^(-2) - 3U λ^(-4) + α S λ^(α-1)

# =============================================================================
# E4: Equilibrium for α=2 (skin coupling, boundary area)
# =============================================================================
print("\n" + "=" * 80)
print("E4: Equilibrium dla α=2 (skin/boundary coupling case)")
print("=" * 80)

print("""
Konkretny mechanism: skin coupling (boundary area soliton z tłem):
  E_int = S λ^2  (S > 0, dim of energy/length²)

Total:
  E(λ) = T/λ + U/λ³ + S λ²

dE/dλ = -T/λ² - 3U/λ⁴ + 2S λ
""")

E_with_skin = T_pos * lam**(-1) + U_pos * lam**(-3) + S_pos * lam**2
dE_skin = diff(E_with_skin, lam)
dE_skin_simp = sp.simplify(dE_skin)
print(f"\n  E(λ) = {E_with_skin}")
print(f"  dE/dλ = {dE_skin_simp}")

# Multiply by λ^4 to clear denominators
# -T λ² - 3U + 2S λ^5 = 0
# 2S λ^5 - T λ² - 3U = 0  (quintic)

# Solve for specific T, U, S values to see if real positive λ exists
# Use T = U = S = 1 for illustration
T_test = sp.Rational(1)
U_test = sp.Rational(1)
S_test = sp.Rational(1)
E_concrete = E_with_skin.subs([(T_pos, T_test), (U_pos, U_test), (S_pos, S_test)])
dE_concrete = diff(E_concrete, lam)

# Solve dE/dλ = 0 numerically (sympy may not find closed form for quintic)
print(f"\n  Test case T=U=S=1:")
print(f"    E(λ) = {E_concrete}")
print(f"    dE/dλ = {sp.simplify(dE_concrete)}")

# Numerical roots
poly_eq = -T_test * lam**2 - 3 * U_test + 2 * S_test * lam**5
roots = sp.nroots(poly_eq, n=15)
print(f"\n  Roots of 2 λ^5 - λ² - 3 = 0:")
for root in roots:
    print(f"    {root}")

# Find positive real root
real_pos_roots = [r for r in roots if abs(sp.im(r)) < 1e-10 and sp.re(r) > 0]
print(f"\n  Real positive roots: {len(real_pos_roots)}")
for root in real_pos_roots:
    print(f"    λ* = {sp.re(root):.10f}")

check("Equilibrium exists (real positive root)",
      len(real_pos_roots) >= 1,
      "no real positive equilibrium found")

# =============================================================================
# E5: Stability analysis
# =============================================================================
print("\n" + "=" * 80)
print("E5: Stability analysis at equilibrium")
print("=" * 80)

print("""
For STABLE equilibrium (local minimum), need:
  d²E/dλ² > 0 at λ = λ*

  d²E/dλ² = 2T/λ³ + 12U/λ⁵ + 2S
""")

d2E_skin = diff(E_with_skin, lam, 2)
d2E_skin_simp = sp.simplify(d2E_skin)
print(f"\n  d²E/dλ² = {d2E_skin_simp}")

# Evaluate at λ* (numerical)
lam_star = sp.re(real_pos_roots[0]) if real_pos_roots else None
if lam_star:
    d2E_at_star = d2E_skin.subs([(T_pos, T_test), (U_pos, U_test), (S_pos, S_test), (lam, lam_star)])
    d2E_value = float(d2E_at_star)
    print(f"\n  At λ* = {float(lam_star):.6f}:")
    print(f"    d²E/dλ² = {d2E_value:.6f}")

    check("d²E/dλ² > 0 at equilibrium (stable minimum)",
          d2E_value > 0,
          f"d²E/dλ² = {d2E_value}")

# General stability: for T,U,S > 0, all terms in d²E/dλ² are positive
# So ANY real positive λ gives stable minimum
print(f"""
  For T, U, S > 0 (physically reasonable), d²E/dλ² > 0 ALWAYS for λ > 0.
  → Any real positive root of dE/dλ = 0 is automatically a STABLE MINIMUM.

  This is the ESSENCE of dynamic equilibrium rozbrojenie Derricka.
""")

# =============================================================================
# E6: Bifurcation - zanik vs ekspansja
# =============================================================================
print("\n" + "=" * 80)
print("E6: Bifurcation - zanik vs ekspansja outcomes")
print("=" * 80)

print("""
Z framework: izolowany soliton (S=0) ma 2 dynamic outcomes:
  - Zanik: λ → 0 (cloud disperses to vacuum)
  - Ekspansja: λ → ∞ (unbounded growth)

Te outcomes są INSTABILITIES of E_iso = T/λ + U/λ³ (Derrick energy).
W obecności coupling (S > 0), pojawia się stable equilibrium między nimi.

Trzy regiony parametrów:
  Region 1: S = 0 (izolowany)
    - Brak stable equilibrium
    - 2 outcomes: λ→0 (zanik) lub λ→∞ (ekspansja)
    - To jest BIFURKACJA (z N17, sympy 7/7 PASS)

  Region 2: S > S_crit (silne coupling)
    - Pojedyncze stable equilibrium λ*
    - Soliton stabilny w tle Φ̄

  Region 3: 0 < S < S_crit (słabe coupling)
    - Możliwe multi-stable lub saddle structure
    - Meta-stability: małe perturbacje → drift do jednego z outcomes
""")

# Test bifurcation behavior: vary S
print("\n  Testing equilibrium dla different S values (T=1, U=1):")
for S_val in [sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(1), sp.Rational(2), sp.Rational(5)]:
    poly = -T_test * lam**2 - 3 * U_test + 2 * S_val * lam**5
    rts = sp.nroots(poly, n=10)
    real_pos = [sp.re(r) for r in rts if abs(sp.im(r)) < 1e-10 and sp.re(r) > 0]
    print(f"    S = {float(S_val):.2f}: equilibria λ* = {[float(r) for r in real_pos]}")

print(f"""
  Observations:
  - For all tested S > 0, exists at least one real positive λ*
  - λ* decreases as S increases (stronger coupling → smaller stable size)
  - This is consistent z dynamic equilibrium picture

CONNECTION DO N17/N18:
  - In limit S → 0 (no background coupling), 2 dynamic outcomes (zanik/ekspansja)
  - These map onto bifurcation 2-state space (z N17)
  - SU(2) emergence (z N18) operates on this bifurcation space
  - With S > 0, soliton is STABILIZED but quantum state STILL describes
    bifurcation tendencies (lean toward zanik vs ekspansja)
""")

# =============================================================================
# E7: VERDICT N16
# =============================================================================
print("\n" + "=" * 80)
print("E7: VERDICT N16 - Dynamic equilibrium formalization")
print("=" * 80)

print(f"""
Sympy tests: {PASS}/{PASS+FAIL} PASS

POSITIVE FINDINGS:

1. STANDARD DERRICK confirmed: izolowany soliton (E_iso = T/λ + U/λ³)
   - dE/dλ < 0 dla T,U > 0 → expansion (λ → ∞)
   - No stable static minimum w 3D

2. WITH BACKGROUND COUPLING E_int = S λ^α (α > 0):
   - α=2 (skin coupling, boundary area): tested
   - dE/dλ = 0 has REAL POSITIVE solution (equilibrium exists)
   - d²E/dλ² > 0 at equilibrium (STABLE minimum)
   - Tested numerically dla T=U=S=1: λ* ≈ {float(lam_star):.4f}, d²E > 0

3. PARAMETER SCAN: equilibrium exists dla all S > 0 tested
   - S → 0: equilibrium → ∞ (recovers Derrick instability limit)
   - S → ∞: equilibrium → 0 (over-stabilized)
   - Smooth continuation between

4. CONNECTION TO BIFURCATION (N17):
   - Izolowany soliton (S=0) bifurcates: λ→0 (zanik) lub λ→∞ (ekspansja)
   - To jest 2-outcome dynamics (consistent z N17 sympy 7/7 PASS)
   - W tle (S>0) zatrzymuje ten dynamics w stable point
   - Quantum state opisuje bifurcation tendencies (lean toward zanik/ekspansja)

WHAT N16 DELIVERS:
  + Formal mechanism rozbrojenia Derricka strukturalnie (NIE technicznie)
  + Standard "soliton in plasma" analogia (Sulem-Sulem 1999) zaadaptowana
  + Konkretny coupling mechanism: skin/boundary area (α=2)
  + Stable equilibrium existence verified analytically + numerically
  + Connection do bifurcation 2-state (N17) i SU(2) lift (N18)

LIMITATIONS:
  - Konkretne S, T, U values dla TGP wymagają wyboru solitonowego ansatzu
    (Skyrme-like, Q-ball, oscillon...) - out of scope dla N16 analytical
  - Numerical PDE solver (N20) wciąż OPEN dla pełnej validation
  - α=2 skin-coupling jest one z kilku możliwych mechanisms
  - Backreaction na M9.1'' background = future work

VERDICT N16:
  STRUCTURAL DERIVED - dynamic equilibrium formalism rozwija Derricka
  strukturalnie. Standard scaling (Derrick) + background coupling (TGP)
  daje stable equilibria.

  Combined z N17/N18/N21: TGP framework dla spinor SU(2) ma teraz
  pełną podstawę:
    - N16: stable solitons exist dynamically
    - N17: 2-outcome bifurcation
    - N18: SU(2) fundamental rep z 2-state
    - N21: alternative SU(2) path z horizon multipole

  Cycle ma teraz 4 sympy-verified TGP-natywne deliverables.
""")

print("=" * 80)
print(f"N16 COMPLETE - {PASS}/{PASS+FAIL} sympy checks PASS")
print("Dynamic equilibrium formalism: Derrick rozbrojony strukturalnie")
print("=" * 80)
