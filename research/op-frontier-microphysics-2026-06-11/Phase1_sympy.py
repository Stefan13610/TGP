# =============================================================================
# op-frontier-microphysics — Phase 1
# Target: F-FM-Q4 — resolve concept par.10.6 Q4: what IS the frontier,
#         precisely mathematically? (Phi-space boundary vs physical edge)
#
# Pre-declared structure (Phase 0 LOCKED 2026-06-11):
#   Candidate set (CLOSED): Q4-A (edge in pre-existing space) / Q4-B (Phi-space
#     boundary only, no spatial locus) / Q4-C (identification) / Q4-NEG.
#   Criteria (value-blind, CLOSED): KQ1 ontology (pozycja C LOCKED) / KQ2
#     finiteness (gamma-3 R = ct LOCKED) / KQ3 field structure (EQ-2 + FCR A-i)
#     / KQ4 E1-inaccessibility (no metric extension beyond R required).
#   Wall machinery: TGP Lagrangian only (forbidden move #9):
#     V(Phi) = (lambda/4)(Phi^2 - Phi_0^2)^2, standard kinetic term.
#   Reference profile: half-kink of the degenerate Z2 limit (declared thin-wall
#     approximation, Coleman-style; Phase 0 par.8(e)(iii)).
#   Bands/targets LOCKED (inherited FCR Phase 0). 0 hardcoded T_pass.
#   v_c / tiebreaker NOT touched here (Phase 2 scope).
# =============================================================================
import numpy as np
import sympy as sp

fp = {}
print("=" * 78)
print("op-frontier-microphysics Phase 1 — F-FM-Q4 (frontier definition)")
print("=" * 78)

x, t, r, lam, Phi0, c_s, eps = sp.symbols("x t r lambda Phi_0 c epsilon",
                                          positive=True)
Phi = sp.symbols("Phi", real=True)
Z_REC = 1090.0
LOG_SPAN = np.log10(1.0 + Z_REC)

V = sp.Rational(1, 4) * lam * (Phi**2 - Phi0**2) ** 2

# -----------------------------------------------------------------------------
# FP1 — Q4 candidate exclusion (sympy-verifiable parts + semantic record)
#   Q4-A: requires background manifold -> fails KQ1 (semantic; concept par.1.1
#         pozycja C LOCKED: space is NOT background -> printed, not sympy)
#   Q4-B: "no spatial locus" -> must fail KQ2/KQ3 sympy-verifiably:
#     (i)  gamma-3 LOCKED locus R(t) = c t is FINITE and positive for finite t
#     (ii) the wall profile has explicit x-dependence with LOCALIZED gradient
#          (energy density of the transition concentrated at a locus)
#   => the frontier HAS a finite spatial locus => Q4-B excluded.
#   Surviving: Q4-C (identification) — the locus of the Phi-transition IS the
#   edge of generated space (no background needed: KQ1; no metric extension
#   beyond R needed: KQ4 — exterior enters only as E1-like reference value).
# -----------------------------------------------------------------------------
print("\nFP1 — Q4 exclusion: finite spatial locus of the Phi-transition")
R_locus = c_s * t
finite_locus = bool(R_locus.subs({c_s: 1, t: 1}).is_finite and
                    sp.simplify(sp.diff(R_locus, t) - c_s) == 0)
# reference profile (degenerate-limit half-kink), bulk side x>0 inward:
delta = sp.sqrt(2) / (sp.sqrt(lam) * Phi0)            # width parameter
prof = Phi0 * sp.tanh(x / delta)                       # Phi_w(x), x in [0, oo)
grad = sp.diff(prof, x)
grad_far = sp.limit(grad, x, sp.oo)                    # -> 0 (delocalized? no)
grad_max_at_0 = sp.simplify(grad.subs(x, 0) - Phi0 / delta) == 0
localized = bool(grad_far == 0 and grad_max_at_0)
print(f"  gamma-3 locus R(t) = {R_locus}: finite for finite t, dR/dt = c: {finite_locus}")
print(f"  wall profile Phi_w = Phi_0 tanh(x/delta), delta = {delta}")
print(f"  gradient localized (max at locus, -> 0 away): {localized}")
print("  => Q4-B (no spatial locus) EXCLUDED by KQ2+KQ3 (finite R = ct; x-localized")
print("     transition). Q4-A EXCLUDED by KQ1 (background manifold contradicts")
print("     pozycja C LOCKED — semantic, concept par.1.1). KQ4: exterior appears ONLY")
print("     as E1-like reference value Phi -> 0, no metric extension beyond R used.")
print("  => Q4-C (IDENTIFICATION) is the unique survivor of the CLOSED set.")
fp["FP1"] = bool(finite_locus and localized)
print(f"  T_pass FP1 (locus finite + transition x-localized) = {fp['FP1']}")

# -----------------------------------------------------------------------------
# FP2 — driving pressure: Delta V = V(0) - V(Phi_0) = lambda Phi_0^4 / 4 > 0
#   (sign convention par.3.6.6: positive = outward push on the E2 bulk;
#    creation at the frontier is energetically downhill — concept par.2.2
#    'gradient => driving force' verified ENERGETICALLY, not just claimed)
# -----------------------------------------------------------------------------
print("\nFP2 — driving pressure Delta V")
DV = sp.simplify(V.subs(Phi, 0) - V.subs(Phi, Phi0))
DV_expected = lam * Phi0**4 / 4
print(f"  Delta V = V(0) - V(Phi_0) = {DV}")
fp["FP2"] = bool(sp.simplify(DV - DV_expected) == 0 and
                 DV.subs({lam: 1, Phi0: 1}) > 0)
print(f"  T_pass FP2 (Delta V = lambda Phi_0^4/4 EXACT, > 0) = {fp['FP2']}")

# -----------------------------------------------------------------------------
# FP3 — wall width: kink solves static EOM; delta = 2/m_Phi EXACT
#   static EOM (degenerate reference): Phi'' = dV/dPhi
#   m_Phi^2 = V''(Phi_0) = 2 lambda Phi_0^2
# -----------------------------------------------------------------------------
print("\nFP3 — wall width and mass scale")
kink = Phi0 * sp.tanh(sp.sqrt(lam / 2) * Phi0 * x)     # full Z2 kink form
eom_res = sp.simplify(sp.diff(kink, x, 2) - sp.diff(V, Phi).subs(Phi, kink))
m_Phi = sp.sqrt(sp.diff(V, Phi, 2).subs(Phi, Phi0))
delta_check = sp.simplify(delta - 2 / m_Phi)
print(f"  EOM residual of kink: {eom_res}")
print(f"  m_Phi = sqrt(V''(Phi_0)) = {sp.simplify(m_Phi)}; delta - 2/m_Phi = {delta_check}")
fp["FP3"] = bool(eom_res == 0 and delta_check == 0)
print(f"  T_pass FP3 (kink solves EOM; delta = 2/m_Phi EXACT) = {fp['FP3']}")

# -----------------------------------------------------------------------------
# FP4 — thin-wall surface tension: sigma = int_0^Phi0 sqrt(2V) dPhi
#   expected (Phase 0 par.8(d)): (2/3) sqrt(lambda/2) Phi_0^3
# -----------------------------------------------------------------------------
print("\nFP4 — surface tension (thin-wall)")
sigma = sp.integrate(sp.sqrt(2 * V), (Phi, 0, Phi0))
sigma_expected = sp.Rational(2, 3) * sp.sqrt(lam / 2) * Phi0**3
print(f"  sigma = {sp.simplify(sigma)}")
fp["FP4"] = bool(sp.simplify(sigma - sigma_expected) == 0)
print(f"  T_pass FP4 (sigma = (2/3) sqrt(lambda/2) Phi_0^3 EXACT) = {fp['FP4']}")

# -----------------------------------------------------------------------------
# FP5 — wall dynamics consistency with gamma-3 (consistency CHECK, not source):
#   thin-wall force per area: F = Delta V - 2 sigma / R  (curvature drag)
#   (i)  critical radius R_c = 2 sigma / Delta V finite, > 0
#   (ii) F > 0 for all R > R_c and F -> Delta V > 0 as R -> oo
#   => d(gamma v)/dt > 0 sustained => v -> c (sup of subluminal velocities)
#   => asymptotically null frontier — CONSISTENT with gamma-3 R = ct LOCKED.
#   thin-wall validity: delta / R(t) -> 0 as t -> oo (Phase 0 par.8(e)(iii)).
# -----------------------------------------------------------------------------
print("\nFP5 — wall dynamics vs gamma-3 (consistency)")
R = sp.symbols("R", positive=True)
F_wall = DV - 2 * sigma_expected / R
R_c = sp.solve(sp.Eq(F_wall, 0), R)
unique_Rc = len(R_c) == 1
R_c0 = sp.simplify(R_c[0]) if unique_Rc else None
F_inf = sp.limit(F_wall, R, sp.oo)
F_pos = sp.simplify(F_wall.subs(R, 2 * R_c0)) if unique_Rc else None
thin_wall = sp.limit(delta / R_locus, t, sp.oo)
print(f"  R_c = 2 sigma / Delta V = {R_c0} (unique: {unique_Rc})")
print(f"  F(R->oo) = {F_inf} = Delta V > 0; F(2 R_c) = {sp.simplify(F_pos)} > 0")
print(f"  thin-wall: delta/(c t) -> {thin_wall} as t -> oo")
print("  => sustained outward push => gamma v grows monotonically => v -> c;")
print("     asymptotically null frontier, CONSISTENT with gamma-3 (LOCKED source).")
fp["FP5"] = bool(unique_Rc and R_c0.subs({lam: 1, Phi0: 1}) > 0
                 and sp.simplify(F_inf - DV) == 0
                 and F_pos.subs({lam: 1, Phi0: 1}) > 0 and thin_wall == 0)
print(f"  T_pass FP5 = {fp['FP5']}")

# -----------------------------------------------------------------------------
# FP6 — NEW micro-result: soliton stability boundary INSIDE the wall
#   m_eff^2(Phi) = V''(Phi) = lambda (3 Phi^2 - Phi_0^2)
#   m_eff^2 > 0  <=>  |Phi| > Phi_0/sqrt(3)   (EXACT)
#   position in reference profile: tanh(x*/delta) = 1/sqrt(3)
#   => x* = delta atanh(1/sqrt(3)) — finite, > 0 (strictly bulk-side of center)
#   PHYSICS: stable massive matter can only exist where |Phi| > Phi_0/sqrt(3);
#   created matter therefore enters the bulk through the INNER wall edge, not
#   at the frontier locus itself. (Phase 2 input; v_c NOT derived here.)
# -----------------------------------------------------------------------------
print("\nFP6 — stability boundary |Phi| = Phi_0/sqrt(3) (NEW, EXACT)")
m_eff2 = sp.diff(V, Phi, 2)
bound = sp.solve(sp.Eq(m_eff2, 0), Phi)
pos_bound = [b for b in bound if b.subs({lam: 1, Phi0: 1}) > 0]
b_expected = Phi0 / sp.sqrt(3)
inside = bool(sp.simplify(b_expected - Phi0).subs({Phi0: 1}) < 0)
x_star = delta * sp.atanh(sp.Rational(1, 1) / sp.sqrt(3))
x_star_num = float(x_star.subs({lam: 1, Phi0: 1}) /
                   delta.subs({lam: 1, Phi0: 1}))
print(f"  m_eff^2(Phi) = {sp.expand(m_eff2)}")
print(f"  zero (positive branch): {pos_bound}; expected Phi_0/sqrt(3)")
print(f"  strictly inside wall (< Phi_0): {inside}; depth x* = {x_star_num:.4f} delta")
fp["FP6"] = bool(len(pos_bound) == 1
                 and sp.simplify(pos_bound[0] - b_expected) == 0 and inside)
print(f"  T_pass FP6 (boundary Phi_0/sqrt(3) EXACT, strictly inside) = {fp['FP6']}")

# -----------------------------------------------------------------------------
# FP7 — pre-derivation validation (par.3.6.9, +-5% standard):
#   indicial eq of C-DERIVED form: p^2 + p/3 - eps = 0 => p+ growing root
#   expected: eps = 2/3 -> p = 2/3 -> log10 G = 2.0253
#             eps = 3/2 -> p = (sqrt(55)-1)/6 -> log10 G = 3.2485
#   Euler residual Delta_bulk(eps) = |3 eps - 2|/4: unique zero at eps = 2/3
# -----------------------------------------------------------------------------
print("\nFP7 — pre-derivation validation (growth machinery + residual zeros)")
p = sp.symbols("p")
p_roots = sp.solve(sp.Eq(p**2 + p / sp.Integer(3) - eps, 0), p)
p_plus = sp.simplify(sp.Max(*p_roots).rewrite(sp.Piecewise)) if False else \
    (sp.Rational(-1, 3) + sp.sqrt(sp.Rational(1, 9) + 4 * eps)) / 2
root_ok = any(sp.simplify(rt - p_plus) == 0 for rt in p_roots)
p_23 = sp.simplify(p_plus.subs(eps, sp.Rational(2, 3)))
p_32 = sp.simplify(p_plus.subs(eps, sp.Rational(3, 2)))
lg_23 = float(p_23) * LOG_SPAN
lg_32 = float(p_32) * LOG_SPAN
exp_23, exp_32 = 2.0253, 3.2485
dev23 = abs(lg_23 - exp_23) / exp_23
dev32 = abs(lg_32 - exp_32) / exp_32
Delta_bulk = sp.Abs(3 * eps - 2) / 4
zeros = sp.solve(sp.Eq(3 * eps - 2, 0), eps)
D_32 = sp.simplify(Delta_bulk.subs(eps, sp.Rational(3, 2)))
print(f"  p+(2/3) = {p_23} (EdS EXACT); p+(3/2) = {p_32} (= (sqrt(55)-1)/6)")
print(f"  log10 G: {lg_23:.4f} (exp {exp_23}, rel dev {dev23:.2e});"
      f" {lg_32:.4f} (exp {exp_32}, rel dev {dev32:.2e})")
print(f"  Delta_bulk zeros: {zeros} (unique at 2/3); Delta_bulk(3/2) = {D_32}")
fp["FP7"] = bool(root_ok
                 and sp.simplify(p_23 - sp.Rational(2, 3)) == 0
                 and sp.simplify(p_32 - (sp.sqrt(55) - 1) / 6) == 0
                 and dev23 < 0.05 and dev32 < 0.05
                 and zeros == [sp.Rational(2, 3)]
                 and sp.simplify(D_32 - sp.Rational(5, 8)) == 0)
print(f"  T_pass FP7 = {fp['FP7']}")

# -----------------------------------------------------------------------------
# FP8 — circularity guard: G_obs absent from ALL Phase 1 expressions
# -----------------------------------------------------------------------------
print("\nFP8 — circularity guard")
all_syms = set()
for e in (DV, sigma, m_Phi, delta, m_eff2, p_plus, Delta_bulk, F_wall, prof):
    all_syms |= e.free_symbols
allowed = {x, t, lam, Phi0, c_s, eps, Phi, R}
print(f"  free symbols across Phase 1: {sorted(str(s) for s in all_syms)}")
fp["FP8"] = bool(all_syms <= allowed)
print(f"  T_pass FP8 (no G_obs, no observed-growth symbol anywhere) = {fp['FP8']}")

# =============================================================================
print("\n" + "=" * 78)
n = sum(fp.values())
for k in sorted(fp, key=lambda s: int(s[2:])):
    print(f"  {k}: {'PASS' if fp[k] else 'FAIL'}")
print(f"\n  TOTAL: {n}/{len(fp)} PASS; hardcoded T_pass = 0/{len(fp)}")
print("  F-FM-Q4: Q4-C (IDENTIFICATION) — frontier = spatial locus R = ct of the")
print("  Phi-saturation transition layer; 'Phi-space boundary' and 'physical edge'")
print("  are ONE object under generated-space ontology => RESOLVED_STRUCTURAL")
print("  (conditional on LOCKED concept-paper ontology; risk R-FM-5 declared).")
print("  Tiebreaker v_c: NOT touched (Phase 2). New Phase-2 input: stability")
print("  boundary |Phi| = Phi_0/sqrt(3) strictly inside wall (FP6).")
print("=" * 78)
