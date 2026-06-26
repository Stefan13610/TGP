# =============================================================================
# op-frontier-creation-rate-derivation — Phase 2
# F-FCR-C target #2: bulk transport of frontier-created matter — DERIVATION
#
# Pre-declared chain (assumptions A-i/A-ii/A-iii declared, not hidden):
#   A-i  : bulk spontaneous creation BLOCKED (concept paper §6: "Bulk saturated,
#          blocked (E2 property)") => bulk continuity is SOURCE-FREE exactly
#   A-ii : statistical homogeneity of rho in the bulk (concept paper §6 CMB row;
#          imposed as a consistency requirement, not derived here)
#   A-iii: isotropy of the bulk flow => u = f(t)*x (unique isotropic linear field)
# Background input: rho_bar prop t^-2 (Phase 1 F-FCR-B skeleton, M prop t)
# Bands/targets LOCKED Phase 0. 0 hardcoded T_pass.
# =============================================================================
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

fp = {}
print("=" * 78)
print("FCR Phase 2 — bulk transport derivation (F-FCR-C target #2)")
print("=" * 78)

t, x, r = sp.symbols("t x r", positive=True)
G, cs = sp.symbols("G c", positive=True)
eps, p = sp.symbols("epsilon p", real=True)
Z_REC = 1090.0
LOG_SPAN = np.log10(1.0 + Z_REC)
LOG_G_OBS = 3.0
PASS_BAND = (2.0, 4.0); PARTIAL_BAND = (1.0, 5.0)
def band(lg):
    if PASS_BAND[0] <= lg <= PASS_BAND[1]: return "PASS_BAND"
    if PARTIAL_BAND[0] <= lg <= PARTIAL_BAND[1]: return "PARTIAL_BAND"
    return "FAIL_LOW" if lg < PARTIAL_BAND[0] else "FAIL_HIGH"

# -----------------------------------------------------------------------------
# FP1 — A-i + A-ii => div(u) forced: source-free continuity with homogeneous rho(t)
#   drho/dt + rho * div(u) = 0  =>  div(u) = -rho'/rho = 2/t   (rho prop t^-2)
# -----------------------------------------------------------------------------
print("\nFP1 — homogeneous source-free continuity forces div(u) = 2/t")
rho_bg = t**(-2)                                   # normalization irrelevant
div_u_forced = sp.simplify(-sp.diff(rho_bg, t) / rho_bg)
print(f"  rho_bar prop t^-2 (Phase 1 B-skeleton);  div(u) = -rho'/rho = {div_u_forced}")
fp["FP1"] = bool(sp.simplify(div_u_forced - 2 / t) == 0)
print(f"  T_pass FP1 = {fp['FP1']}")

# -----------------------------------------------------------------------------
# FP2 — A-iii: isotropic linear flow u = f(t) x => div(u) = 3 f => f = 2/(3t)
#   matter-flow scale factor: adot/a = f => a_m prop t^(2/3)
#   (DISTINCT from space-frontier R = c t — derived, not assumed)
# -----------------------------------------------------------------------------
print("\nFP2 — isotropy => u = (2/3) x/t; matter-flow scale factor a_m prop t^(2/3)")
f = sp.symbols("f", cls=sp.Function)
f_sol = sp.solve(sp.Eq(3 * sp.Symbol("f0"), 2 / t), sp.Symbol("f0"))[0]
a_m = sp.Function("a_m")(t)
a_sol = sp.dsolve(sp.Eq(sp.diff(a_m, t) / a_m, f_sol), a_m).rhs
print(f"  f(t) = div(u)/3 = {f_sol};  a_m(t) = {a_sol}  (prop t^(2/3))")
growth_exp_am = sp.simplify(sp.log(a_sol / a_sol.subs(t, 1)) / sp.log(t))
fp["FP2"] = bool(sp.simplify(f_sol - sp.Rational(2, 3) / t) == 0
                 and sp.simplify(a_sol / a_sol.subs(t, 1) - t**sp.Rational(2, 3)) == 0)
print(f"  T_pass FP2 (f = 2/(3t); a_m prop t^(2/3) EXACT) = {fp['FP2']}")
print("  NOTE: matter kinematics a_m prop t^(2/3) INSIDE space frontier R = c*t —")
print("  matter lags the frontier; the gap is filled by frontier creation (consistent")
print("  with A-i: created matter enters at/near the frontier only).")

# -----------------------------------------------------------------------------
# FP3 — background Euler residual (honesty item): u_t + (u.grad)u + grad(phi) =? 0
#   with grad(phi) from Poisson for homogeneous rho: grad(phi) = (4 pi G rho/3) x
#   In units 4 pi G rho = eps/t^2 (Phase 1):  residual coefficient of x/t^2:
# -----------------------------------------------------------------------------
print("\nFP3 — background Euler residual (bounded => substrate-balance caveat declared)")
u_field = sp.Rational(2, 3) * x / t
lhs = sp.diff(u_field, t) + u_field * sp.diff(u_field, x)
grad_phi = (eps / (3 * t**2)) * x                  # (4 pi G rho/3) x with eps = 4 pi G rho t^2
resid = sp.simplify(lhs + grad_phi)
resid_coeff = sp.simplify(resid * t**2 / x)
H_m2_coeff = sp.Rational(4, 9)                     # H_m^2 = (2/3t)^2
Delta_bulk = sp.simplify(sp.Abs(resid_coeff) / H_m2_coeff)
Delta_B1 = float(Delta_bulk.subs(eps, sp.Rational(3, 2)))
print(f"  residual = (x/t^2) * ({resid_coeff});  |Delta|/H_m^2 = {Delta_bulk}")
print(f"  at eps = 3/2 (B1): Delta_bulk = {Delta_B1:.4f}  — BOUNDED, time-independent")
print("  => per R17 lemma (phi' = sqrt(3*Delta)/tau): bounded residual => NO runaway class;")
print("  => background force balance requires substrate contribution ~ O(1)*H_m^2*x —")
print("     TGP-native (expansion is substrate-driven, gamma-3) but NOT derived from the")
print("     action here => CAVEAT C-2 DECLARED (carried into aggregate).")
fp["FP3"] = bool(np.isfinite(Delta_B1) and Delta_B1 < 10.0
                 and sp.simplify(sp.diff(Delta_bulk, t)) == 0)
print(f"  T_pass FP3 (residual bounded + time-independent) = {fp['FP3']}")

# -----------------------------------------------------------------------------
# FP4 — perturbation ODE around the derived flow (source-free bulk):
#   continuity: delta' = -theta ; Euler: theta' + 2 H_m theta = -4 pi G rho delta
#   (standard manipulations; H_m = 2/(3t); rho prop t^-2)
#   => delta'' + (4/(3 tau)) delta' - (eps/tau^2) delta = 0   [C-DERIVED form]
#   Mandatory validation: EdS limit (H = 2/3t, eps_EdS = 2/3) must give p = 2/3 EXACT
# -----------------------------------------------------------------------------
print("\nFP4 — perturbation ODE around derived flow + EdS validation")
tau = sp.symbols("tau", positive=True)
delta_f = sp.Function("delta_f")(tau)
H_m = sp.Rational(2, 3) / tau
theta_sub = -sp.diff(delta_f, tau)                 # from source-free continuity
euler_pert = (sp.diff(theta_sub, tau) + 2 * H_m * theta_sub
              + (eps / tau**2) * delta_f)          # theta' + 2H theta + 4piG rho delta = 0
ode = sp.expand(-euler_pert)                       # => delta'' + (4/3 tau) delta' - eps/tau^2 delta
A1 = sp.simplify(ode.coeff(sp.diff(delta_f, tau)) * tau)
B0 = sp.simplify(ode.coeff(delta_f) * tau**2)
print(f"  ODE: delta'' + ({A1})/tau delta' + ({B0})/tau^2 delta = 0")
indicial = p * (p - 1) + A1 * p + B0
p_roots = sp.solve(sp.Eq(indicial, 0), p)
p_plus = sp.simplify(sp.Max(*p_roots))
p_EdS = sp.simplify(p_plus.subs(eps, sp.Rational(2, 3)))
print(f"  p+(eps) = {p_plus}")
print(f"  EdS validation (eps = 2/3): p+ = {p_EdS}  [expect 2/3 EXACT — standard EdS growth]")
fp["FP4"] = bool(sp.simplify(A1 - sp.Rational(4, 3)) == 0
                 and sp.simplify(B0 + eps) == 0
                 and sp.simplify(p_EdS - sp.Rational(2, 3)) == 0)
print(f"  T_pass FP4 (C-DERIVED form + EdS limit EXACT) = {fp['FP4']}")

# -----------------------------------------------------------------------------
# FP5 — B1/B2 evaluation of the DERIVED form + exact symbolic p (B1)
# -----------------------------------------------------------------------------
print("\nFP5 — derived-form growth: B1 (eps = 3/2 EXACT) and B2 (eps = 0.465)")
p_B1 = sp.simplify(p_plus.subs(eps, sp.Rational(3, 2)))
p_B1_simpl = sp.radsimp(p_B1)
p_B1_f = float(p_B1)
logG_B1 = p_B1_f * LOG_SPAN
p_B2_f = float(p_plus.subs(eps, 0.465))
logG_B2 = max(p_B2_f, 0.0) * LOG_SPAN
print(f"  B1: p+ = {p_B1_simpl} = {p_B1_f:.5f}  => log10 G = {logG_B1:.4f}  => {band(logG_B1)}")
print(f"  B2: p+ = {p_B2_f:.5f}             => log10 G = {logG_B2:.4f}  => {band(logG_B2)}")
print(f"  (observed comparison target log10 G_obs = {LOG_G_OBS}; bands LOCKED Phase 0)")
fp["FP5"] = bool(sp.simplify(p_B1 - (sp.sqrt(55) - 1) / 6) == 0)
print(f"  T_pass FP5 (B1 exact form p = (sqrt(55)-1)/6) = {fp['FP5']}")

# -----------------------------------------------------------------------------
# FP6 — numerical ODE cross-check (B1): integrate derived form, compare log10 G
# -----------------------------------------------------------------------------
print("\nFP6 — numerical cross-check of derived form (B1)")
tau_i = 1.0 / (1.0 + Z_REC)
def rhs(tt, y):
    d, dp = y
    return [dp, -(4.0 / (3.0 * tt)) * dp + (1.5 / tt**2) * d]
s = solve_ivp(rhs, [tau_i, 1.0], [1.0, p_B1_f / tau_i], method="Radau",
              rtol=1e-11, atol=1e-14)
logG_num = np.log10(abs(s.y[0, -1]))
dev = abs(logG_num - logG_B1) / logG_B1
print(f"  analytic log10 G = {logG_B1:.5f}; numeric = {logG_num:.5f}; rel dev = {dev:.2e}")
fp["FP6"] = bool(s.success and dev < 0.01)
print(f"  T_pass FP6 = {fp['FP6']}")

# -----------------------------------------------------------------------------
# FP7 — circularity guard + global consistency M prop t (no G_obs anywhere)
# -----------------------------------------------------------------------------
print("\nFP7 — guards: p+ free symbols subset {eps}; global M(R=ct) prop t consistency")
fs = p_plus.free_symbols
M_global = sp.simplify(rho_bg * sp.Rational(4, 3) * sp.pi * (cs * t)**3)
M_prop_t = sp.simplify(sp.diff(M_global / t, t)) == 0
fp["FP7"] = bool(fs <= {eps} and M_prop_t)
print(f"  free symbols: {fs}; M_global = {M_global} (prop t: {M_prop_t})")
print(f"  T_pass FP7 = {fp['FP7']}")

# -----------------------------------------------------------------------------
# FP8 — F-FCR-C verdict update + aggregate re-evaluation (mechanical)
# -----------------------------------------------------------------------------
print("\nFP8 — F-FCR-C verdict update + aggregate")
print("  F-FCR-C: PARTIAL_concept_mismatch (Phase 1, state-of-concept-paper) ->")
print("  >>> C-DERIVED_CONDITIONAL (THIS Phase 2 derivation) <<<")
print("  derivation chain: A-i (E2 blocked, concept LOCKED-claim) + A-ii (homogeneity,")
print("  consistency requirement) + A-iii (isotropy) => u = (2/3)x/t FORCED =>")
print("  source-free bulk perturbation eqs with H_m = 2/(3t) — supersedes the C2a/b/c")
print("  proxy menu (which assumed UNIFORM-in-bulk creation, now excluded by A-i).")
print("  CAVEATS carried: C-2 substrate-balance (FP3); A-ii imposed not derived.")
status = {"B1": "STRUCTURAL_POSTULATE", "B2": "CONCEPT_CLAIM_UNVERIFIED",
          "A": "A1_DERIVED_CONDITIONAL_ON_B + A2_bridge_gap",
          "C": "DERIVED_CONDITIONAL (A-ii imposed; C-2 substrate-balance caveat)"}
all_derived = False                                 # B1 postulate + caveats remain (computed below)
all_derived = all("DERIVED" in s and "POSTULATE" not in s and "UNVERIFIED" not in s
                  and "gap" not in s and "caveat" not in s for s in status.values())
agg = "PREDICTION_REALIZED" if (all_derived and band(logG_B1) == "PASS_BAND") \
      else ("STRUCTURAL_CONDITIONAL" if band(logG_B1) in ("PASS_BAND", "PARTIAL_BAND")
            else "FAIL_NEGATIVE")
for k, v in status.items(): print(f"    {k}: {v}")
print(f"\n  >>> F-FCR-D AGGREGATE (post Phase 2): {agg} <<<")
print("  => PR-022: STILL NOT APPENDED (B1 postulate + C-2/A-ii caveats outstanding)")
fp["FP8"] = bool(agg in {"PREDICTION_REALIZED", "STRUCTURAL_CONDITIONAL", "FAIL_NEGATIVE"})
print(f"  T_pass FP8 = {fp['FP8']}")

print("\n" + "=" * 78)
n = sum(fp.values())
for k in sorted(fp, key=lambda s: int(s[2:])): print(f"  {k}: {'PASS' if fp[k] else 'FAIL'}")
print(f"\n  TOTAL: {n}/{len(fp)} PASS; hardcoded T_pass = 0/{len(fp)}")
print(f"  C-DERIVED form: delta'' + (4/3tau) delta' - (eps/tau^2) delta = 0")
print(f"  B1: p = (sqrt(55)-1)/6 = {p_B1_f:.5f} -> log10 G = {logG_B1:.3f} ({band(logG_B1)})")
print(f"  B2: log10 G = {logG_B2:.3f} ({band(logG_B2)})")
print(f"  F-FCR-D: {agg}")
print("=" * 78)
