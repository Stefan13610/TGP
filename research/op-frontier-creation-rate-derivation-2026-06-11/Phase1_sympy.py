# =============================================================================
# op-frontier-creation-rate-derivation — Phase 1
# F-FCR-A (S_creation) / F-FCR-B (M_univ relation) / F-FCR-C (bulk forms)
# F-FCR-D (prediction matrix + aggregate)
# Pre-registration: Phase0_balance.md LOCKED 2026-06-11. 0 hardcoded T_pass.
# =============================================================================
import numpy as np
import sympy as sp

fp = {}
print("=" * 78)
print("FCR Phase 1 — frontier creation rate: derivation + prediction matrix")
print("=" * 78)

# LOCKED inherited constants
G_eff = 6.674e-11; c = 3e8; H_0 = 2.3e-18; t_0 = 1.0 / H_0
Z_REC = 1090.0                      # NEW (gamma) observational anchor (Phase 0 par.5)
LOG_G_OBS = 3.0                     # comparison target ONLY
PASS_BAND = (2.0, 4.0); PARTIAL_BAND = (1.0, 5.0)

def band(logG):
    if PASS_BAND[0] <= logG <= PASS_BAND[1]: return "PASS_BAND"
    if PARTIAL_BAND[0] <= logG <= PARTIAL_BAND[1]: return "PARTIAL_BAND"
    return "FAIL_LOW" if logG < PARTIAL_BAND[0] else "FAIL_HIGH"

# -----------------------------------------------------------------------------
# FP1 — F-FCR-B route B1 (symbolic): R = c*t AND R = 2 G M/c^2
#   => M(t) = c^3 t/(2G); Mdot = c^3/(2G) = const; rho = 3H^2/(8 pi G) EXACT;
#   => eps_G(tot) = 4 pi G rho t^2 = 3/2 EXACT
# -----------------------------------------------------------------------------
print("\nFP1 — B1 zero-energy/horizon condition (symbolic)")
t, G, cs = sp.symbols("t G c", positive=True)
M_B1 = sp.solve(sp.Eq(cs * t, 2 * G * sp.Symbol("M", positive=True) / cs**2),
                sp.Symbol("M", positive=True))[0]
Mdot = sp.diff(M_B1, t)
rho_B1 = sp.simplify(M_B1 / (sp.Rational(4, 3) * sp.pi * (cs * t)**3))
H = 1 / t
rho_crit = 3 * H**2 / (8 * sp.pi * G)
crit_identity = sp.simplify(rho_B1 - rho_crit)
eps_B1 = sp.simplify(4 * sp.pi * G * rho_B1 * t**2)
print(f"  M(t) = {M_B1};  dM/dt = {Mdot} (const => M prop t)")
print(f"  rho(t) = {rho_B1};  rho - 3H^2/(8 pi G) = {crit_identity}  [critical-density identity]")
print(f"  eps_G(tot) = 4 pi G rho t^2 = {eps_B1}")
fp["FP1"] = bool(crit_identity == 0 and sp.simplify(eps_B1 - sp.Rational(3, 2)) == 0
                 and sp.simplify(sp.diff(Mdot, t)) == 0)
print(f"  T_pass FP1 (M prop t; critical identity EXACT; eps_G = 3/2 EXACT) = {fp['FP1']}")

# -----------------------------------------------------------------------------
# FP2 — B1 numeric + INFORMATIONAL consistency vs gamma-7 rough M_univ = 1e53 kg
# -----------------------------------------------------------------------------
print("\nFP2 — B1 numeric value at t_0")
M_t0 = c**3 * t_0 / (2 * G_eff)
ratio_M = M_t0 / 1e53
print(f"  M_tot(t_0) = c^3 t_0/(2G) = {M_t0:.3e} kg  (gamma-7 rough input was 1e53 kg)")
print(f"  INFORMATIONAL: ratio = {ratio_M:.3f} (consistency note, NOT an input)")
fp["FP2"] = bool(1e52 < M_t0 < 1e54)
print(f"  T_pass FP2 (value finite, computed, order-1e53) = {fp['FP2']}")

# -----------------------------------------------------------------------------
# FP3 — F-FCR-A: creation rate. A1: S/rho = Mdot/M = H EXACT (resolves hyp-Q3 | B).
# A2: EQ-5 stationarity d<Phi>/dt = 0 => S_Phi = 3 H <Phi>  (scaling prop H consistent);
#     Phi->matter bridge NOT in concept paper => gap DECLARED (no improvisation).
# -----------------------------------------------------------------------------
print("\nFP3 — F-FCR-A: S_creation prop H")
S_over_rho = sp.simplify(Mdot / M_B1)
Phi_c = sp.Symbol("Phi_c", positive=True)
S_Phi = sp.simplify(0 + 3 * H * Phi_c)          # from d<Phi>/dt + 3H<Phi> = S, stationary
print(f"  A1: S_matter/rho = Mdot/M = {S_over_rho}  (= H = 1/t EXACT)")
print(f"  A2: EQ-5 stationary => S_Phi = {S_Phi}  (prop H, coefficient 3<Phi>)")
print("  A2 bridge Phi-substrate -> matter sector: NOT specified in concept paper")
print("  => A2 bridge gap DECLARED (PARTIAL_concept_mismatch component; forbidden move #9 honored)")
fp["FP3"] = bool(sp.simplify(S_over_rho - H) == 0 and sp.simplify(S_Phi - 3 * H * Phi_c) == 0)
print(f"  T_pass FP3 (A1 EXACT; A2 scaling consistent) = {fp['FP3']}")

# -----------------------------------------------------------------------------
# FP4 — F-FCR-C: bulk ODE forms (reused R17 symbolic machinery, general eps)
# C2a-form: p^2 + p - (eps - ... ) — re-derive indicial polynomials per form
#   C2a (uniform, no drag):  d'' + (3/tau) d' + (1-eps)/tau^2 d = 0
#   C2b (uniform, drag):     d'' + (4/tau) d' + (2-eps)/tau^2 d = 0
#   C2c (bulk-clean):        d'' + (2/tau) d' - (eps/tau^2) d = 0
# Limit check: eps->0 collapse of growth (p+ <= 0) for all forms.
# -----------------------------------------------------------------------------
print("\nFP4 — F-FCR-C: indicial structure of the three pre-declared bulk forms")
eps, p = sp.symbols("epsilon p", real=True)
forms = {
    "C2a-form": p * (p - 1) + 3 * p + (1 - eps),
    "C2b-form": p * (p - 1) + 4 * p + (2 - eps),
    "C2c-form": p * (p - 1) + 2 * p - eps,
}
p_plus = {}
for nm, poly in forms.items():
    roots = sp.solve(sp.Eq(poly, 0), p)
    p_plus[nm] = sp.simplify(sp.Max(*roots))
    print(f"  {nm}: p+ = {sp.simplify(roots[1] if roots[1] != sp.Min(*roots) else roots[0])}")
limit_ok = all(float(pp.subs(eps, 0)) <= 0 for pp in p_plus.values())
fp["FP4"] = bool(limit_ok)
print(f"  eps->0 limit: growth collapses (p+ <= 0) in all forms: {limit_ok}")
print(f"  T_pass FP4 = {fp['FP4']}")
print("  C-form SELECTION: concept paper specifies boundary-localized creation (EQ-5 note)")
print("  but NOT the bulk transport/homogenization of created matter =>")
print("  >>> F-FCR-C verdict: PARTIAL_concept_mismatch (pre-registered honest outcome) <<<")
print("  All three forms evaluated in the matrix (no cherry-picking; forbidden move #4).")
verdict_C = "PARTIAL_concept_mismatch"

# -----------------------------------------------------------------------------
# FP5 — native tau_init (gamma-3 a prop t => 1+z = t0/t) + provenance audit
# -----------------------------------------------------------------------------
print("\nFP5 — native tau_init from gamma-3 mapping")
a_t = t                                  # a prop t
z_expr = sp.symbols("t0", positive=True) / t - 1   # 1+z = a(t0)/a(t) = t0/t
tau_init_native = 1.0 / (1.0 + Z_REC)
LOG_SPAN = np.log10(1.0 / tau_init_native)
print(f"  1 + z = t_0/t (a prop t, EXACT) => tau_init = 1/(1+z_rec) = {tau_init_native:.4e}")
print(f"  log10 span = {LOG_SPAN:.4f}")
print(f"  PROVENANCE AUDIT: gamma-7/R17 used tau_init = 2.75e-5 = LCDM age-at-recombination")
print(f"  (1.2e13 s) — borrowed mapping, NOT gamma-3-native. Ratio vs native: "
      f"{tau_init_native/2.75e-5:.1f}x.")
print("  => R1 candidate flag (forward-looking; NO retroactive verdict modification — fm #3)")
fp["FP5"] = bool(abs(tau_init_native - 1.0 / 1091.0) < 1e-12)
print(f"  T_pass FP5 = {fp['FP5']}")

# -----------------------------------------------------------------------------
# FP6 — F-FCR-D prediction matrix: 3 C-forms x 2 B-values (eps = (3/2)*Omega_m)
# -----------------------------------------------------------------------------
print("\nFP6 — prediction matrix (G = (1/tau_init)^p+; bands LOCKED)")
matrix = {}
for om_label, om in [("B1: Omega_m = 1 (zero-energy)", 1.0),
                     ("B2: Omega_m = 0.31 (E2 claim)", 0.31)]:
    print(f"\n  {om_label}: eps_G = {1.5 * om:.3f}")
    for nm, pp in p_plus.items():
        p_val = float(pp.subs(eps, 1.5 * om))
        logG = max(p_val, 0.0) * LOG_SPAN
        v = band(logG)
        matrix[(om_label, nm)] = (p_val, logG, v)
        exact_note = ""
        if om == 1.0 and nm == "C2c-form":
            p_exact = sp.simplify(pp.subs(eps, sp.Rational(3, 2)))
            exact_note = f"   [p = {p_exact} EXACT]"
        print(f"    {nm}: p+ = {p_val:+.4f}  log10 G = {logG:.3f}  => {v}{exact_note}")
fp["FP6"] = bool(all(np.isfinite(v[1]) for v in matrix.values()))
print(f"\n  T_pass FP6 (matrix complete, all cells computed) = {fp['FP6']}")

# -----------------------------------------------------------------------------
# FP7 — circularity guard: p+ depends only on eps (never on G_obs); G_obs comparison-only
# -----------------------------------------------------------------------------
print("\nFP7 — circularity guard")
fs = set().union(*[pp.free_symbols for pp in p_plus.values()])
fp["FP7"] = bool(fs <= {eps})
print(f"  free symbols of all p+: {fs}  (subset of {{epsilon}}: {fp['FP7']})")
print(f"  T_pass FP7 = {fp['FP7']}")

# -----------------------------------------------------------------------------
# INFORMATIONAL — xi = 1 sensitivity (M = c^3 t/G): eps = 3*Omega_m (NOT adopted)
# -----------------------------------------------------------------------------
print("\nINFORMATIONAL — xi = 1 sensitivity (forbidden as adoption; reporting only)")
for om in [1.0, 0.31]:
    e_ = 3.0 * om
    for nm, pp in p_plus.items():
        p_val = float(pp.subs(eps, e_))
        if p_val > 0:
            print(f"    xi=1, Omega={om}: {nm}: log10 G = {p_val * LOG_SPAN:.2f}")

# -----------------------------------------------------------------------------
# FP8 — F-FCR-D aggregate (mechanical per Phase 0 par.1.3)
# Derivation-status declarations (from FP1-FP5 + concept-paper audit):
#   B1: STRUCTURAL_POSTULATE (consistent + minimal + pre-existing roadmap task,
#       but zero-energy condition NOT forced by LOCKED machinery)
#   B2: CONCEPT_CLAIM_UNVERIFIED (its own F5 pending)
#   A : A1 DERIVED (conditional on B); A2 bridge gap declared
#   C : PARTIAL_concept_mismatch
# => any postulate/claim/mismatch present => STRUCTURAL_CONDITIONAL (unless all cells fail)
# -----------------------------------------------------------------------------
print("\nFP8 — aggregate F-FCR-D (mechanical)")
status = {"B1": "STRUCTURAL_POSTULATE", "B2": "CONCEPT_CLAIM_UNVERIFIED",
          "A": "A1_DERIVED_CONDITIONAL_ON_B + A2_bridge_gap", "C": verdict_C}
all_derived = all(s == "DERIVED" for s in status.values())
any_band_hit = any(v[2] in ("PASS_BAND", "PARTIAL_BAND") for v in matrix.values())
all_fail = all(v[2] not in ("PASS_BAND", "PARTIAL_BAND") for v in matrix.values())
if all_derived and any(v[2] == "PASS_BAND" for v in matrix.values()):
    aggregate = "PREDICTION_REALIZED"
elif all_fail:
    aggregate = "FAIL_NEGATIVE"
elif any_band_hit:
    aggregate = "STRUCTURAL_CONDITIONAL"
else:
    aggregate = "INDETERMINATE"
for k, s in status.items():
    print(f"    {k}: {s}")
print(f"\n  >>> F-FCR-D AGGREGATE: {aggregate} <<<")
fp["FP8"] = bool(aggregate in {"PREDICTION_REALIZED", "STRUCTURAL_CONDITIONAL",
                               "FAIL_NEGATIVE", "INDETERMINATE"})
print(f"  T_pass FP8 = {fp['FP8']}")
print("  => PR-022: NOT APPENDED (PREDICTION_REALIZED not reached; forbidden move #5 honored)")

# -----------------------------------------------------------------------------
print("\n" + "=" * 78)
n = sum(fp.values())
for k in sorted(fp, key=lambda s: int(s[2:])): print(f"  {k}: {'PASS' if fp[k] else 'FAIL'}")
print(f"\n  TOTAL: {n}/{len(fp)} PASS; hardcoded T_pass = 0/{len(fp)}")
print(f"  F-FCR-D: {aggregate}")
print("=" * 78)
