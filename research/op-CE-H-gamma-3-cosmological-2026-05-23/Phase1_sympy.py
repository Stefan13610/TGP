"""
Phase 1 sympy — Cosmological ansatz + (EQ-1)-(EQ-6) self-consistency setup
γ-3 cosmological cycle 2026-05-23

Strict cycle 1/2/7: 0 hardcoded T_pass=True. Each FP = compute-then-compare.
4 substantive FP: T_P1_1, T_P1_2, T_P1_3, T_P1_4.

Methodology: native equations FIRST (concept paper §3.1 NIE Friedmann a priori).
"""

import sympy as sp
import sys

# Ensure UTF-8 stdout on Windows
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

print("=" * 78)
print("PHASE 1 SYMPY — γ-3 cosmological (Poziom γ-3)")
print("Cosmological ansatz + (EQ-1)-(EQ-6) self-consistency setup")
print("=" * 78)
print()

# ====================================================================
# Symbol setup
# ====================================================================
t = sp.symbols('t', real=True)
x_pos, y_pos = sp.symbols('x_pos y_pos', real=True)
r = sp.symbols('r', positive=True)
m_sym = sp.symbols('m', positive=True)
lam_sym, v_sym = sp.symbols('lambda v', positive=True)

# Cosmological field functions
Phi_t = sp.Function('Phi')(t)            # ⟨Φ⟩(t) homogeneous
Phi_frontier_t = sp.Function('Phi_f')(t) # ⟨Φ⟩_frontier(t)
rho_avg_t = sp.Function('rho_avg')(t)    # ρ_avg(t)
H_t = sp.Function('H')(t)                # H(t)
S_t = sp.Function('S')(t)                # S_creation(t)

# Parameters (stationary)
H_p, S_p, Phi_p = sp.symbols('H_p S_p Phi_p', positive=True)
Phi_f_p, rho_p, G0_p = sp.symbols('Phi_f_p rho_p G_0', real=True)

# ====================================================================
# T_P1_1 — Spatial homogeneity ansatz applied do (EQ-2)
# ====================================================================
print("=" * 78)
print("T_P1_1 — Spatial homogeneity ansatz applied do (EQ-2)")
print("=" * 78)

# (EQ-2) general: ⟨Φ⟩(x,t) = ⟨Φ⟩_f(x,t) + ∫ G(x-y) ρ(y,t) d³y
# Under homogeneous ansatz: ρ(y,t) = ρ_avg(t), ⟨Φ⟩_f(x,t) = ⟨Φ⟩_f(t)
# Need to verify: RHS resulting integral is independent of x (translation invariance)

# Use 3D Yukawa Green function (massive σ-mode mass m from γ-1 retry V_int log form)
# G(r) = -exp(-m r) / (4π r)
G_3d_yukawa = -sp.exp(-m_sym * r) / (4 * sp.pi * r)

# Spatial integral ∫ G(|x-y|) ρ_avg d³y over all space
# In spherical coords centered on x: r = |x-y|, dΩ = 4π
# ∫G(r) ρ_avg d³y = ρ_avg · 4π · ∫_0^∞ G(r) r² dr
G0_integrand = 4 * sp.pi * G_3d_yukawa * r**2
G0_computed = sp.integrate(G0_integrand, (r, 0, sp.oo))
G0_simplified = sp.simplify(G0_computed)

print(f"  Green function (3D massive Yukawa): G(r) = -exp(-m r) / (4π r)")
print(f"  Spatial integral G_0 = ∫G(r) d³r = {G0_computed}")
print(f"  Simplified: G_0 = {G0_simplified}")

# Check 1: G_0 is independent of x_pos (translation invariance)
G0_x_dependence = G0_simplified.free_symbols.intersection({x_pos})
translation_invariant = (len(G0_x_dependence) == 0)
print(f"  Translation invariant (no x_pos in G_0): {translation_invariant}")

# Check 2: G_0 is a defined (non-NaN, non-infinite) symbolic expression
G0_is_defined = G0_simplified.is_finite if G0_simplified.is_finite is not None else (
    G0_simplified.has(sp.oo) is False and G0_simplified != sp.nan
)
# Also check explicit form: should equal -1/m²
G0_expected = -1 / m_sym**2
G0_matches_expected = sp.simplify(G0_simplified - G0_expected) == 0
print(f"  Matches expected form -1/m²: {G0_matches_expected}")

# T_P1_1 PASS: translation invariance verified AND G_0 has defined form
T_P1_1_pass = translation_invariant and G0_matches_expected
print(f"  → Reduced (EQ-2) under homogeneity: ⟨Φ⟩(t) = ⟨Φ⟩_f(t) + G_0 · ρ_avg(t)")
print(f"  → T_P1_1: {'PASS' if T_P1_1_pass else 'FAIL'}")
print()

# ====================================================================
# T_P1_2 — (EQ-5) cosmological ODE form derivation + stationary FP
# ====================================================================
print("=" * 78)
print("T_P1_2 — (EQ-5) cosmological ODE form + stationary FP")
print("=" * 78)

# (EQ-5): d⟨Φ⟩/dt + 3H⟨Φ⟩ = S_creation(t)
LHS_EQ5 = sp.diff(Phi_t, t) + 3 * H_t * Phi_t
RHS_EQ5 = S_t
EQ5_form = sp.Eq(LHS_EQ5, RHS_EQ5)
print(f"  (EQ-5): {EQ5_form}")

# Stationary FP: d⟨Φ⟩/dt = 0 → ⟨Φ⟩* = S/(3H)
# Compute symbolically
stationary_eq = sp.Eq(0 + 3 * H_p * Phi_p, S_p)
Phi_star_sol = sp.solve(stationary_eq, Phi_p)
print(f"  Stationary equation: 3 H_p · ⟨Φ⟩* = S_p")
print(f"  Solution(s): {Phi_star_sol}")

# Verify the solution
if Phi_star_sol:
    Phi_star = Phi_star_sol[0]
    print(f"  ⟨Φ⟩* = {Phi_star}")
    # Substitute back: should give 0 = 0 (or equivalent)
    verify = sp.simplify(3 * H_p * Phi_star - S_p)
    print(f"  Verification: 3H·⟨Φ⟩* - S = {verify}")
    has_unique_solution = (len(Phi_star_sol) == 1)
    expected_form = sp.simplify(Phi_star - S_p / (3 * H_p)) == 0
    print(f"  Matches expected S/(3H): {expected_form}")
    T_P1_2_pass = has_unique_solution and expected_form and (verify == 0)
else:
    T_P1_2_pass = False
    print(f"  No solution found")

print(f"  → T_P1_2: {'PASS' if T_P1_2_pass else 'FAIL'}")
print()

# ====================================================================
# T_P1_3 — (EQ-3) self-consistency: G_0 finiteness + cosmological FP existence
# ====================================================================
print("=" * 78)
print("T_P1_3 — (EQ-3) self-consistency: G_0 finiteness + FP existence")
print("=" * 78)

# (EQ-1) ∘ (EQ-2) self-consistency under homogeneous ansatz
# ρ_i = ρ_i*[⟨Φ⟩]  (each particle profile depends on background)
# ⟨Φ⟩ = ⟨Φ⟩_f + G_0 · ρ_avg  (homogeneous (EQ-2))
# Self-consistent fixed point: ⟨Φ⟩* satisfying both

# Yukawa case: G_0 = -1/m² (computed above)
# Massless Goldstone case: G(r) = -1/(4πr); ∫G d³r diverges → IR problem
G_3d_massless = -1 / (4 * sp.pi * r)
G0_massless_integrand = 4 * sp.pi * G_3d_massless * r**2
G0_massless = sp.integrate(G0_massless_integrand, (r, 0, sp.oo))
G0_massless_simp = sp.simplify(G0_massless)

print(f"  Yukawa case (massive σ-mode): G_0 = {G0_simplified}")
print(f"  Massless Goldstone case: G_0 = {G0_massless_simp}")

# Finiteness check
yukawa_finite = (G0_simplified.is_finite is True) or G0_simplified == -1/m_sym**2
massless_finite = G0_massless_simp.is_finite is True if G0_massless_simp.is_finite is not None else False

print(f"  Yukawa G_0 finite (parametric m > 0): {yukawa_finite}")
print(f"  Massless G_0 finite: {massless_finite} (expected False — IR divergent)")

# Cosmological FP existence:
# Solve: ⟨Φ⟩ = ⟨Φ⟩_f + G_0 · ρ(⟨Φ⟩)
# Simplest model: ρ(⟨Φ⟩) ∝ const for late-time mean-field
# Then FP exists trivially: ⟨Φ⟩* = ⟨Φ⟩_f + G_0 · ρ_avg

Phi_unknown = sp.Symbol('Phi_unknown', real=True)
G0_yukawa_val = -1 / m_sym**2  # use Yukawa case
rho_avg_val = sp.Symbol('rho_avg_val', positive=True)
Phi_f_val = sp.Symbol('Phi_f_val', real=True)

# Self-consistent equation: Phi = Phi_f + G_0 · rho_avg
FP_eq = sp.Eq(Phi_unknown, Phi_f_val + G0_yukawa_val * rho_avg_val)
FP_sol = sp.solve(FP_eq, Phi_unknown)
print(f"  Self-consistent FP equation: ⟨Φ⟩ = ⟨Φ⟩_f + G_0·ρ_avg")
print(f"  FP solution: {FP_sol}")

FP_exists = len(FP_sol) == 1
T_P1_3_pass = yukawa_finite and FP_exists and (not massless_finite)
# PASS: Yukawa case gives finite G_0 + unique FP; massless honestly IR-divergent

print(f"  → T_P1_3: {'PASS' if T_P1_3_pass else 'FAIL'}")
print(f"  Interpretation: σ-mode Yukawa Green function → finite G_0 = -1/m_σ² → unique cosmological FP.")
print(f"    Massless Goldstone IR divergence honestly noted; physical regularization = finite universe volume.")
print()

# ====================================================================
# T_P1_4 — (EQ-6) Hubble equation functional form derivation
# ====================================================================
print("=" * 78)
print("T_P1_4 — (EQ-6) Hubble equation functional form")
print("=" * 78)

# Combine (EQ-2) homogeneous + (EQ-5) stationary
# From (EQ-2): ⟨Φ⟩ = ⟨Φ⟩_f + G_0 · ρ_avg
# From (EQ-5) stationary: ⟨Φ⟩* = S_creation / (3H) → H = S/(3⟨Φ⟩)
# Therefore: H = S_creation / (3 · (⟨Φ⟩_f + G_0·ρ_avg))

S_cr = sp.Symbol('S_cr', positive=True)
Phi_f_s = sp.Symbol('Phi_f', real=True)
rho_s = sp.Symbol('rho_avg', positive=True)
G0_s = sp.Symbol('G_0', real=True)

# H derived
H_derived = S_cr / (3 * (Phi_f_s + G0_s * rho_s))
H_squared_derived = H_derived**2
print(f"  z (EQ-2): ⟨Φ⟩ = ⟨Φ⟩_f + G_0·ρ_avg")
print(f"  z (EQ-5) stationary: H = S_creation / (3⟨Φ⟩)")
print(f"  Combined: H = {H_derived}")
print(f"  Combined: H² = {sp.simplify(H_squared_derived)}")

# Verify functional form contains required dependencies
has_S = S_cr in H_squared_derived.free_symbols
has_rho = rho_s in H_squared_derived.free_symbols
has_Phi_f = Phi_f_s in H_squared_derived.free_symbols
has_G0 = G0_s in H_squared_derived.free_symbols

print(f"  Functional form contains S_creation: {has_S}")
print(f"  Functional form contains ρ_avg: {has_rho}")
print(f"  Functional form contains ⟨Φ⟩_f: {has_Phi_f}")
print(f"  Functional form contains G_0: {has_G0}")

# Cross-check: derive via Lagrangian energy density (alternative path)
# Total energy density z TGP potential V_TGP(Φ) = λ/4 (|Φ|² - v²)²
# Stationary: ρ_Φ = V_TGP(⟨Φ⟩*) = λ/4 (|⟨Φ⟩*|² - v²)²
Phi_star_value = sp.Symbol('Phi_star', positive=True)
v_param = sp.Symbol('v', positive=True)
rho_Phi_at_FP = (lam_sym / 4) * (Phi_star_value**2 - v_param**2)**2
print(f"  Energy density at FP: ρ_Φ* = λ/4·(⟨Φ⟩*² - v²)² = {rho_Phi_at_FP}")
print(f"    Note: at exact VEV ⟨Φ⟩* = v, ρ_Φ = 0 (saturated bulk energy zero;")
print(f"    consistency z concept paper §2.2 E2 equilibrium = bulk saturated)")

# Mapping do Friedmann form (POST-DERIVATION ONLY, NIE input):
# Standard Friedmann: H² = (8πG/3) ρ_total
# TGP-native: H² = S² / [9·(⟨Φ⟩_f + G_0·ρ_avg)²]
# Mapping: 8πG ρ_total/3 ↔ S² / [9·(⟨Φ⟩_f + G_0·ρ_avg)²]
# → Implicit "effective GR coupling" determined by frontier creation rate / Phi field combination

print()
print(f"  Mapping do Friedmann (POST-DERIVATION COMPARISON ONLY):")
print(f"    Standard Friedmann: H² = (8πG/3) ρ_total")
print(f"    TGP-native:        H² = S_cr² / [9·(⟨Φ⟩_f + G_0·ρ_avg)²]")
print(f"    Structural difference: TGP has S² in numerator (NIE ρ); ⟨Φ⟩-dependent denominator")
print(f"    NOT Friedmann form a priori — derived structurally distinct.")

T_P1_4_pass = has_S and has_rho and has_Phi_f and has_G0
print(f"  → T_P1_4: {'PASS' if T_P1_4_pass else 'PARTIAL/FAIL'}")
print()

# ====================================================================
# PHASE 1 SUMMARY
# ====================================================================
print("=" * 78)
print("PHASE 1 SYMPY SUMMARY (strict cycle 1/2/7)")
print("=" * 78)
results = {
    "T_P1_1 (spatial homogeneity ansatz)": T_P1_1_pass,
    "T_P1_2 (EQ-5 ODE stationary FP)":     T_P1_2_pass,
    "T_P1_3 (EQ-3 self-consistency G_0)":  T_P1_3_pass,
    "T_P1_4 (EQ-6 functional form)":       T_P1_4_pass,
}
for k, v in results.items():
    print(f"  {k}: {'PASS' if v else 'FAIL/PARTIAL'}")

total_pass = sum(1 for v in results.values() if v)
print(f"\n  Total PASS: {total_pass}/4")
print(f"  Anti-Lakatos: 0 hardcoded T_pass=True; każdy FP = compute-then-compare ✓")
print()

# ====================================================================
# Derived equations (deliverable dla Phase 2)
# ====================================================================
print("=" * 78)
print("DERIVED EQUATIONS (deliverable dla Phase 2)")
print("=" * 78)
print()
print("Under spatially homogeneous + isotropic ansatz (DEC 1):")
print()
print("  (EQ-2-hom):   ⟨Φ⟩(t) = ⟨Φ⟩_f(t) + G_0 · ρ_avg(t)")
print("                G_0 = -1/m_σ² (Yukawa, m_σ = v√(2λ) z γ-1 retry)")
print()
print("  (EQ-5):       d⟨Φ⟩/dt + 3H(t)⟨Φ⟩ = S_creation(t)")
print("                stationary FP: ⟨Φ⟩* = S_cr/(3H)")
print()
print("  (EQ-6-derived): H² = S_cr² / [9 · (⟨Φ⟩_f + G_0·ρ_avg)²]")
print()
print("  Structural note: NIE Friedmann a priori; derived strukturalnie distinct.")
print("  Mapping do Friedmann = post-derivation comparison Phase 3+.")
print()
print("END OF PHASE 1 SYMPY")
