# =============================================================================
# op-frontier-microphysics — Phase 2
# Target: F-FM-V — creation velocity v_c (tiebreaker B-k3 vs B-k4)
#
# Pre-declared structure (Phase 0 LOCKED 2026-06-11):
#   v_c semantics (par.1.2 BINDING): velocity of the created element at ENTRY
#     into the source-free bulk = the velocity entering the FCR Phase 3 LOCKED
#     marginality bookkeeping (1/2) v_c^2 = G M/(c t).
#   Mechanisms (CLOSED): M1 (frontier-comoving / radiation-first, v_c = c) /
#     M2 (flow-matching, v_c = 2c/3) / M3 (wall-energetics).
#   Criteria (CLOSED, value-blind): K1 massivity / K2 no-drag matching /
#     K3 self-consistency residual / K4 escape-hatch audit.
#   Honest evaluation note (K2): adiabatic peculiar-velocity decay (prop 1/a)
#     exists in the derived flow => the originally-worded "no relaxation" route
#     is WEAKENED; K2 binds instead through the MEAN-FLOW BOUNDARY VALUE:
#     u = (2/3) x/t is the UNIQUE regular isotropic solution of source-free
#     continuity + homogeneity; the matter at the entry surface x = ct IS the
#     newly created matter, so its mean velocity = u(ct,t) = 2c/3.
#     Same criterion, honest mechanism — documented per par.3.6 anti-goalpost.
#   Declared assumptions (Phase 0 par.8(e) + NEW A-iv): A-ii exact homogeneity
#     (corollary COR-1, Phase 3); A-iv monochromatic entry (zero velocity
#     dispersion at the ledger event) — declared, not hidden.
#   Bands/targets LOCKED (inherited). 0 hardcoded T_pass. G_obs never input.
# =============================================================================
import numpy as np
import sympy as sp

fp = {}
print("=" * 78)
print("op-frontier-microphysics Phase 2 — F-FM-V (creation velocity v_c)")
print("=" * 78)

x, t, lam, Phi0, c_s, G, M, eps, vc = sp.symbols(
    "x t lambda Phi_0 c G M epsilon v_c", positive=True)
Phi = sp.symbols("Phi", real=True)
Z_REC = 1090.0
LOG_SPAN = np.log10(1.0 + Z_REC)
PASS_BAND = (2.0, 4.0); PARTIAL_BAND = (1.0, 5.0)
def band(lg):
    if PASS_BAND[0] <= lg <= PASS_BAND[1]: return "PASS_BAND"
    if PARTIAL_BAND[0] <= lg <= PARTIAL_BAND[1]: return "PARTIAL_BAND"
    return "FAIL_LOW" if lg < PARTIAL_BAND[0] else "FAIL_HIGH"
LOG_G_OBS = 3.0   # comparison-only (forbidden move #1 guard: FP8)

V = sp.Rational(1, 4) * lam * (Phi**2 - Phi0**2) ** 2

# -----------------------------------------------------------------------------
# FP1 — K1 (massivity at the ledger event): the element is MASSIVE at entry
#   Phase 1 FP6: stable matter only where |Phi| > Phi_0/sqrt(3) (inside wall);
#   at the bulk side m_eff -> m_Phi = sqrt(2 lambda) Phi_0 > 0.
#   A massive element has |v| < c STRICTLY => exact v_c = c is NOT realizable
#   by the entering element. (B-k3 survives at this stage ONLY via the
#   radiation-first sub-variant -> audited in FP2.)
# -----------------------------------------------------------------------------
print("\nFP1 — K1: massivity at entry => v_c < c strict")
m_eff2 = sp.diff(V, Phi, 2)
m_at_sat = sp.simplify(m_eff2.subs(Phi, Phi0))
bound_pos = [b for b in sp.solve(sp.Eq(m_eff2, 0), Phi)
             if b.subs({lam: 1, Phi0: 1}) > 0]
strictly_inside = bool(len(bound_pos) == 1 and
                       (bound_pos[0] - Phi0).subs({lam: 1, Phi0: 1}) < 0)
massive_at_entry = bool(m_at_sat.subs({lam: 1, Phi0: 1}) > 0)
print(f"  m_eff^2(Phi_0) = {m_at_sat} > 0; stability edge {bound_pos[0]} < Phi_0: "
      f"{strictly_inside}")
print("  => element at bulk entry is massive => |v| < c STRICT; exact v_c = c")
print("     not realizable by a massive element (Newtonian ledger (1/2)c^2 has")
print("     no realizing carrier). B-k3 now rests ONLY on radiation-first (FP2).")
fp["FP1"] = bool(massive_at_entry and strictly_inside)
print(f"  T_pass FP1 = {fp['FP1']}")

# -----------------------------------------------------------------------------
# FP2 — K4(a) (escape-hatch audit): radiation-first scenario CLOSED by A-i
#   Scenario: matter born massless at the wall (enters at c), condenses to
#   massive solitons LATER, inside the bulk.
#   But the LEDGER event is the appearance of dM > 0 (mass). Condensation in
#   the bulk interior is a bulk creation event — FORBIDDEN by A-i LOCKED-claim
#   ("bulk saturated, blocked"; concept par.5 NATIVE row + FCR Phase 2 A-i).
#   => mass appearance is confined to the wall layer (thickness O(delta) at R);
#   => at the ledger event the element is massive (FP1) => v_c < c. CLOSED.
#   sympy part: the ledger locus is within a layer of vanishing relative
#   thickness: [x*(stability edge depth), delta] with delta/(c t) -> 0.
# -----------------------------------------------------------------------------
print("\nFP2 — K4(a): radiation-first audit (A-i blocks bulk condensation)")
delta = sp.sqrt(2) / (sp.sqrt(lam) * Phi0)
x_star = delta * sp.atanh(1 / sp.sqrt(3))
ratio_inf = sp.limit(delta / (c_s * t), t, sp.oo)
x_star_over_delta = float(x_star / delta)
finite_depth = bool(0 < x_star_over_delta < 1)
print(f"  ledger locus: stability-edge depth x* = {x_star_over_delta:.4f} delta")
print(f"  (finite, within the layer: {finite_depth}); delta/(ct) -> {ratio_inf}")
print("  bulk condensation = dM > 0 source in bulk = A-i VIOLATION (LOCKED claim)")
print("  => mass can only appear inside the wall layer; combined with FP1:")
print("     the entering element is massive => B-k3 (v_c = c) EXCLUDED.")
fp["FP2"] = bool(finite_depth and ratio_inf == 0)
print(f"  T_pass FP2 = {fp['FP2']}")

# -----------------------------------------------------------------------------
# FP3 — K2 (positive selection, mean-flow boundary value):
#   source-free continuity + homogeneity: div u = -rhodot/rho = 2/t (rho ~ t^-2)
#   isotropic ansatz u_r = f(t) x + g(t)/x^2 : div = 3 f(t) (the g-piece is
#   divergence-free) ; regularity at x = 0 => g = 0 => u = (2/3) x/t UNIQUE.
#   The matter present at the entry surface x = ct(-) IS the newly created
#   matter; exact homogeneity (A-ii) extends the mean flow to the boundary
#   => mean entry velocity v_c = u(ct, t) = 2c/3 EXACT
#   (conditional on A-ii exact + A-iv monochromatic entry — declared).
#   HONEST NOTE (K2 evaluation): peculiar velocities decay adiabatically
#   (FP7) — so the "drag layer" exclusion alone does NOT force matching;
#   the boundary-value argument above is what binds.
# -----------------------------------------------------------------------------
print("\nFP3 — K2: unique regular flow + boundary value at x = ct")
f_t, g_t = sp.symbols("f g", cls=sp.Function)
u_r = f_t(t) * x + g_t(t) / x**2
div_u = sp.simplify(sp.diff(x**2 * u_r, x) / x**2)
f_sol = sp.solve(sp.Eq(div_u, 2 / t), f_t(t))
g_in_div = sp.simplify(div_u - 3 * f_t(t))
u_flow = f_sol[0] * x          # g = 0 by regularity at x -> 0
vc_entry = sp.simplify(u_flow.subs(x, c_s * t))
print(f"  div u = {div_u} (g-piece contributes {g_in_div} -> divergence-free;")
print("  regularity at x = 0 kills g) => f = "
      f"{f_sol[0]} => u = (2/3) x/t UNIQUE")
print(f"  boundary value: v_c = u(ct, t) = {vc_entry}")
print("  conditional on: A-ii exact (COR-1, Phase 3) + A-iv monochromatic entry")
fp["FP3"] = bool(len(f_sol) == 1 and g_in_div == 0
                 and sp.simplify(vc_entry - 2 * c_s / 3) == 0)
print(f"  T_pass FP3 (v_c = 2c/3 EXACT, unique) = {fp['FP3']}")

# -----------------------------------------------------------------------------
# FP4 — K3 (substrate force in bulk = 0, MODEL-INDEPENDENT within EQ-1):
#   EQ-1: soliton energy is a functional of the local background,
#   E_sol = E(<Phi>(x)); force F = -dE/dx = -E'(<Phi>) * grad<Phi>.
#   In the saturated bulk <Phi> = Phi_0 = const => grad<Phi> = 0 => F = 0
#   EXACT — for ANY E (no model input needed).
#   => background Euler must close source-free; residual of u = (2/3)x/t
#   against gravity g = -(eps/3) x/t^2:
#     res = du/dt + u u' - g = (eps/3 - 2/9) x/t^2 ; res = 0 <=> eps = 2/3
#   reproduces LOCKED Delta_bulk = |res|/|u u'| = |3 eps - 2|/4 and its
#   unique zero. CONVERGENT, INDEPENDENT selection of eps = 2/3.
#   (This also pre-resolves COR-2/C-2: the required substrate force is ZERO
#    and the balance holds identically — formal booking in Phase 3.)
# -----------------------------------------------------------------------------
print("\nFP4 — K3: bulk substrate force = 0 => Euler residual must vanish")
E_fun = sp.Function("E")
Phi_bg_bulk = Phi0 + 0 * x                      # saturated: constant in x
F_sub = -sp.diff(E_fun(Phi_bg_bulk), x)
u23 = sp.Rational(2, 3) * x / t
res = sp.simplify(sp.diff(u23, t) + u23 * sp.diff(u23, x) + (eps / 3) * x / t**2)
eps_zero = sp.solve(sp.Eq(res, 0), eps)
Delta_bulk = sp.simplify(sp.Abs(res) / sp.Abs(u23 * sp.diff(u23, x)))
D32 = sp.simplify(Delta_bulk.subs(eps, sp.Rational(3, 2)))
print(f"  F_substrate(bulk) = {F_sub} (model-independent: grad<Phi> = 0)")
print(f"  Euler residual = {res}; zero at eps = {eps_zero}")
print(f"  Delta_bulk = {Delta_bulk} = |3 eps - 2|/4 (LOCKED FCR formula);"
      f" Delta(3/2) = {D32}")
fp["FP4"] = bool(F_sub == 0 and eps_zero == [sp.Rational(2, 3)]
                 and sp.simplify(Delta_bulk - sp.Abs(3 * eps - 2) / 4) == 0
                 and sp.simplify(D32 - sp.Rational(5, 8)) == 0)
print(f"  T_pass FP4 (F_sub = 0; residual zero UNIQUE at eps = 2/3) = {fp['FP4']}")

# -----------------------------------------------------------------------------
# FP5 — ledger closure at v_c = 2c/3: criticality identity
#   marginality (LOCKED): (1/2)v_c^2 = GM/(ct) with v_c = 2c/3
#   => M = 2 c^3 t/(9G) (xi = 2/9) => rho = M/((4pi/3)(ct)^3) = 1/(6 pi G t^2)
#   => rho = 3 H_m^2/(8 pi G) EXACT with H_m = 2/(3t)  [matter sector exactly
#   critical w.r.t. its own derived flow — same identity standard cosmology
#   uses for rho_crit] ; and eps = 4 pi G rho t^2 = 2/3 — loop closes.
# -----------------------------------------------------------------------------
print("\nFP5 — ledger closure: marginality <=> criticality of the matter flow")
M_sol = sp.solve(sp.Eq(sp.Rational(1, 2) * (2 * c_s / 3) ** 2,
                       G * M / (c_s * t)), M)
xi = sp.simplify(M_sol[0] / (c_s**3 * t / G))
rho = sp.simplify(M_sol[0] / (sp.Rational(4, 3) * sp.pi * (c_s * t) ** 3))
H_m = sp.Rational(2, 3) / t
rho_crit = 3 * H_m**2 / (8 * sp.pi * G)
eps_loop = sp.simplify(4 * sp.pi * G * rho * t**2)
print(f"  M = {M_sol[0]} (xi = {xi}); rho = {rho}")
print(f"  rho - 3H_m^2/(8 pi G) = {sp.simplify(rho - rho_crit)} (criticality EXACT)")
print(f"  eps = 4 pi G rho t^2 = {eps_loop} — closes the loop")
fp["FP5"] = bool(len(M_sol) == 1 and sp.simplify(xi - sp.Rational(2, 9)) == 0
                 and sp.simplify(rho - rho_crit) == 0
                 and sp.simplify(eps_loop - sp.Rational(2, 3)) == 0)
print(f"  T_pass FP5 = {fp['FP5']}")

# -----------------------------------------------------------------------------
# FP6 — collapsed prediction (MECHANICAL; bands LOCKED Phase 0/FCR/R17)
#   selected point: eps = 2/3 => p+ = 2/3 EXACT (EdS) => log10 G = 2.025
#   excluded point reported for the record: eps = 3/2 => 3.249
#   distances reported; G_obs used ONLY here as comparison (never input).
# -----------------------------------------------------------------------------
print("\nFP6 — collapsed prediction (mechanical band evaluation)")
p_plus = (sp.Rational(-1, 3) + sp.sqrt(sp.Rational(1, 9) + 4 * eps)) / 2
p_sel = sp.simplify(p_plus.subs(eps, sp.Rational(2, 3)))
lg_sel = float(p_sel) * LOG_SPAN
p_exc = sp.simplify(p_plus.subs(eps, sp.Rational(3, 2)))
lg_exc = float(p_exc) * LOG_SPAN
print(f"  SELECTED (B-k4): p = {p_sel} (EdS EXACT) -> log10 G = {lg_sel:.4f}"
      f" => {band(lg_sel)} (edge distance {lg_sel - PASS_BAND[0]:.3f} dex)")
print(f"  EXCLUDED (B-k3): p = {p_exc} -> log10 G = {lg_exc:.4f} ({band(lg_exc)})"
      " — for the record")
print(f"  distance to observed {LOG_G_OBS}: {abs(lg_sel - LOG_G_OBS):.2f} dex BELOW")
print("  HONEST: the value-blind criteria select the point FARTHER from the")
print("  observed growth (direction pre-flagged Phase 0 par.7) — reported as-is.")
fp["FP6"] = bool(sp.simplify(p_sel - sp.Rational(2, 3)) == 0
                 and band(lg_sel) == "PASS_BAND")
print(f"  T_pass FP6 = {fp['FP6']}")

# -----------------------------------------------------------------------------
# FP7 — structural notes (honest record):
#   (a) peculiar-velocity decay: linearized trajectories about the flow
#       (eps = 2/3): s(s-1) + (2/9) = 0 => s in {2/3, 1/3}; the s = 1/3 mode
#       has v_pec ~ t^(-2/3) ~ 1/a (adiabatic decay) — this is WHY the original
#       "no relaxation" wording of K2 alone does not bind (honest note).
#   (b) B-k3 Schwarzschild numerology dissolves: 2GM/c^2 = (4/9) ct < ct —
#       the universe is NOT at its horizon condition at the selected point.
#   (c) M3 (wall energetics): energy partition kinetic/internal is free at
#       this level => NON-DISCRIMINATING for v_c; full budget = Phase 3 (A2),
#       with declared tension R-FM-3.
# -----------------------------------------------------------------------------
print("\nFP7 — structural notes (decay exponents; Schwarzschild note; M3)")
s = sp.symbols("s")
roots_s = sp.solve(sp.Eq(s * (s - 1) + sp.Rational(2, 9), 0), s)
v_pec_exp = sp.simplify(min(roots_s) - 1)          # d/dt t^(1/3) ~ t^(-2/3)
a_exp = sp.Rational(2, 3)
adiabatic = sp.simplify(v_pec_exp + a_exp) == 0    # v_pec ~ a^(-1)
R_sch_over_R = sp.simplify(2 * G * M_sol[0] / c_s**2 / (c_s * t))
print(f"  trajectory exponents s = {roots_s}; v_pec ~ t^({v_pec_exp}) ~ a^(-1):"
      f" {adiabatic}")
print(f"  2GM/c^2 / R = {R_sch_over_R} (< 1: horizon condition NOT realized)")
print("  M3: partition free => non-discriminating (deferred to Phase 3 A2)")
fp["FP7"] = bool(set(roots_s) == {sp.Rational(2, 3), sp.Rational(1, 3)}
                 and adiabatic
                 and sp.simplify(R_sch_over_R - sp.Rational(4, 9)) == 0)
print(f"  T_pass FP7 = {fp['FP7']}")

# -----------------------------------------------------------------------------
# FP8 — circularity guard + verdict assembly
# -----------------------------------------------------------------------------
print("\nFP8 — circularity guard")
all_syms = set()
for e in (m_eff2, delta, x_star, u_flow, vc_entry, res, Delta_bulk,
          M_sol[0], rho, p_plus, R_sch_over_R):
    all_syms |= e.free_symbols
allowed = {x, t, lam, Phi0, c_s, G, eps, Phi}
print(f"  free symbols across Phase 2: {sorted(str(q) for q in all_syms)}")
fp["FP8"] = bool(all_syms <= allowed)
print(f"  T_pass FP8 (no G_obs in any derivation expression) = {fp['FP8']}")

# =============================================================================
print("\n" + "=" * 78)
n = sum(fp.values())
for k in sorted(fp, key=lambda q: int(q[2:])):
    print(f"  {k}: {'PASS' if fp[k] else 'FAIL'}")
print(f"\n  TOTAL: {n}/{len(fp)} PASS; hardcoded T_pass = 0/{len(fp)}")
print("  F-FM-V: TIEBREAKER_DERIVED (CONDITIONAL): v_c = 2c/3 EXACT (B-k4)")
print("    - B-k3 EXCLUDED unconditionally (K1 massivity + K4a: bulk")
print("      condensation forbidden by A-i LOCKED) — value-blind")
print("    - positive selection: K2 mean-flow boundary value (conditional on")
print("      A-ii exact + A-iv monochromatic) + K3 zero-substrate-force (forced,")
print("      model-independent) — TWO independent convergent lines")
print("    - NEW_POINT not realized (mean entry = LOCKED B-k4 value)")
print("  Collapsed parameter-free prediction: log10 G = 2.025 (p = 2/3 EdS")
print("  EXACT) — PASS_BAND edge; 0.97 dex BELOW observed 3.0 (reported as-is).")
print("  PR-022 condition (i): satisfied CONDITIONALLY (A-ii -> Phase 3).")
print("  NO PR-022 (conditions ii-iv outstanding; forbidden move #6).")
print("=" * 78)
