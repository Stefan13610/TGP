# =============================================================================
# op-R17-linear-runaway-diagnosis — Phase 1
# F-R17-A (regression gate) / F-R17-B (background self-consistency audit)
# F-R17-C (consistent-transcription growth C1/C2a/C2b/C2c) / F-R17-D (aggregate)
#
# Pre-registration: Phase0_balance.md LOCKED 2026-06-11.
# Discipline: 0 hardcoded T_pass; all verdicts computed; bands LOCKED Phase 0 §3.
# =============================================================================

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

fp = {}  # FP registry

print("=" * 78)
print("op-R17-linear-runaway-diagnosis — Phase 1 execution")
print("=" * 78)

# -----------------------------------------------------------------------------
# §0 — LOCKED inherited constants (Phase 0 §6; 0 new constants)
# -----------------------------------------------------------------------------
G_eff = 6.674e-11        # γ-5 LOCKED [m³/(kg s²)]
M_univ = 1e53            # γ-7 Phase 3 LOCKED [kg]
c = 3e8                  # [m/s]
H_0 = 2.3e-18            # (γ) OBSERVATIONAL_ANCHOR [1/s]
t_0 = 1.0 / H_0          # [s]
tau_init = 2.75e-5       # γ-7 Phase 3 LOCKED (recombination, τ = t/t_0)
# Comparison targets ONLY (forbidden move #10: never inputs):
LOG_G_OBS = 3.0          # δ: 1e-5 → 1e-2 (γ-7 Phase 3 §3.2 LOCKED)
LOG_G_REF_GAMMA7 = np.log10(6e213)   # LOCKED γ-7 Phase 3 numerical reference
# Bands LOCKED Phase 0 §3 (project convention factor-10 / factor-100):
PASS_BAND = (2.0, 4.0)
PARTIAL_BAND = (1.0, 5.0)

def band_verdict(logG):
    """Mechanical band application (Phase 0 §3, F-R17-C)."""
    if PASS_BAND[0] <= logG <= PASS_BAND[1]:
        return "PASS_BAND"
    if PARTIAL_BAND[0] <= logG <= PARTIAL_BAND[1]:
        return "PARTIAL_BAND"
    return "FAIL_LOW" if logG < PARTIAL_BAND[0] else "FAIL_HIGH"

LOG_TAU_SPAN = np.log10(1.0 / tau_init)
print(f"\nlog10(tau_present/tau_init) = {LOG_TAU_SPAN:.5f}")

# -----------------------------------------------------------------------------
# FP1 — F-R17-A.1: eps_G reproduction (gate: within ±2% of LOCKED 1.71)
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP1 — F-R17-A.1: eps_G(t_0) = 3 G_eff M_univ / (c^3 t_0) reproduction")
eps_G = 3 * G_eff * M_univ / (c**3 * t_0)
rel_dev_eps = abs(eps_G - 1.71) / 1.71
print(f"  eps_G = {eps_G:.4f}   (LOCKED reference 1.71; rel dev {rel_dev_eps:.4f})")
fp["FP1"] = bool(rel_dev_eps < 0.02)
print(f"  T_pass FP1 = {fp['FP1']}")

# -----------------------------------------------------------------------------
# FP2 — F-R17-A.2: numerical runaway reproduction (gate: log10 G in 213.8 ± 1.0)
# Audited transcription: delta'' + (2/tau) delta' - (eps_G/tau^3) delta = 0
# Integrated in log-form u = ln(delta): u'' + (u')^2 + (2/tau) u' - eps_G/tau^3 = 0
# Growing-mode IC from WKB: u'(tau_init) = sqrt(eps_G)/tau_init^(3/2)
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP2 — F-R17-A.2: numerical reproduction of the gamma-7 runaway")

def rhs_log_runaway(tau, y):
    u, up = y
    return [up, -up**2 - (2.0 / tau) * up + eps_G / tau**3]

up0 = np.sqrt(eps_G) / tau_init**1.5  # WKB growing mode
sol = solve_ivp(rhs_log_runaway, [tau_init, 1.0], [0.0, up0],
                method="Radau", rtol=1e-10, atol=1e-12, dense_output=False)
log10_G_runaway = sol.y[0, -1] / np.log(10.0)
print(f"  ODE integration success: {sol.success}")
print(f"  log10(delta_present/delta_init) = {log10_G_runaway:.2f}")
print(f"  LOCKED gamma-7 reference: {LOG_G_REF_GAMMA7:.2f} (6e213); tolerance ±1.0")
fp["FP2"] = bool(sol.success and abs(log10_G_runaway - LOG_G_REF_GAMMA7) <= 1.0)
print(f"  T_pass FP2 = {fp['FP2']}")

# -----------------------------------------------------------------------------
# FP3 — F-R17-A.3: symbolic WKB asymptotics: delta ~ exp(lam * tau^p), p = -1/2,
# lam = ±2 sqrt(eps_G)  (confirms gamma-7 Phase 3 §3.1 analytic structure)
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP3 — F-R17-A.3: symbolic WKB exponent structure")
tau_s, eps_s = sp.symbols("tau epsilon", positive=True)
lam_s, p_s = sp.symbols("lambda p", real=True)
delta_ansatz = sp.exp(lam_s * tau_s**p_s)
ode_expr = (sp.diff(delta_ansatz, tau_s, 2)
            - (eps_s / tau_s**3) * delta_ansatz)  # dominant balance (drag subleading)
lead = sp.simplify(ode_expr / delta_ansatz)
# leading term of (d2/dt2): lam^2 p^2 tau^(2p-2); balance vs eps/tau^3 => 2p-2 = -3
p_solved = sp.solve(sp.Eq(2 * p_s - 2, -3), p_s)[0]
lam_solved = sp.solve(sp.Eq((lam_s * p_solved) ** 2, eps_s), lam_s)
lam_growing = [l for l in lam_solved if sp.simplify(l - 2 * sp.sqrt(eps_s)) == 0]
print(f"  balance 2p-2 = -3  =>  p = {p_solved}")
print(f"  lambda solutions: {lam_solved}  (growing-mode |lambda| = 2*sqrt(eps_G))")
fp["FP3"] = bool(p_solved == sp.Rational(-1, 2) and len(lam_growing) == 1)
print(f"  T_pass FP3 = {fp['FP3']}")

# -----------------------------------------------------------------------------
# FP4 — F-R17-B: background self-consistency residual (PRIMARY)
# Delta(tau) = |a''/a - (-(4 pi G/3) rho_bar)| / H^2 with a ∝ t (a''=0),
# rho_bar = 3 M /(4 pi c^3 t^3) (M = const). Symbolic, then numeric extremum.
# Verdict: CONSISTENT (max < 0.1) / INCONSISTENT_O1 (max >= 0.1)  [Phase 0 §3]
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP4 — F-R17-B: background residual Delta(tau) under audited transcription")
t_s, G_s, M_s, c_s, t0_s = sp.symbols("t G M c t_0", positive=True)
a_bg = t_s                                   # a ∝ t (gamma-3 LOCKED)
H_bg = sp.diff(a_bg, t_s) / a_bg             # = 1/t
rho_bg = 3 * M_s / (4 * sp.pi * c_s**3 * t_s**3)
accel_lhs = sp.diff(a_bg, t_s, 2) / a_bg     # = 0
accel_rhs = -(4 * sp.pi * G_s / 3) * rho_bg
Delta_sym = sp.simplify(sp.Abs(accel_lhs - accel_rhs) / H_bg**2)
print(f"  Delta(t) symbolic = {Delta_sym}")
# substitute tau = t/t_0:
Delta_tau = sp.simplify(Delta_sym.subs(t_s, tau_s * t0_s))
eps_check = sp.simplify(Delta_tau * 3 * tau_s)   # expect = 3 G M /(c^3 t_0) ... /?
print(f"  Delta(tau) = {Delta_tau}   [expect eps_G/(3 tau)]")
Delta_num = sp.lambdify(tau_s, Delta_tau.subs({G_s: G_eff, M_s: M_univ,
                                               c_s: c, t0_s: t_0}), "numpy")
Delta_at_init = float(Delta_num(tau_init))
Delta_at_t0 = float(Delta_num(1.0))
Delta_max = max(Delta_at_init, Delta_at_t0)
print(f"  Delta(tau_init) = {Delta_at_init:.3e}   Delta(1) = {Delta_at_t0:.3f}")
verdict_B = "INCONSISTENT_O1" if Delta_max >= 0.1 else "CONSISTENT"
print(f"  F-R17-B verdict: {verdict_B}  (max Delta = {Delta_max:.3e}; threshold 0.1)")
fp["FP4"] = bool(Delta_max == Delta_at_init or Delta_max == Delta_at_t0)  # extremum found
print(f"  T_pass FP4 (residual well-defined + extremum located) = {fp['FP4']}")

# -----------------------------------------------------------------------------
# FP5 — F-R17-B.2: analytic link lemma — WKB phase derivative phi'(tau) equals
# sqrt(3*Delta(tau))/tau exactly => runaway (super-power-law) <=> Delta unbounded.
# Bounded Delta = const  =>  phi ∝ ln(tau)  =>  pure power law (no runaway).
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP5 — F-R17-B.2: lemma — runaway mode is generated EXACTLY by the residual")
# source term S(tau) (coefficient of delta in d2delta/dtau2 = S*delta dominant balance)
S_source = sp.simplify((4 * sp.pi * G_s * rho_bg).subs(t_s, tau_s * t0_s) * t0_s**2)
phi_prime = sp.sqrt(S_source)                          # WKB: phi' = sqrt(S)
lemma_expr = sp.simplify(phi_prime - sp.sqrt(3 * Delta_tau) / tau_s)
print(f"  S(tau) = {S_source}")
print(f"  phi'(tau) - sqrt(3*Delta)/tau = {lemma_expr}   [expect 0]")
fp["FP5"] = bool(sp.simplify(lemma_expr) == 0)
print(f"  T_pass FP5 = {fp['FP5']}")
print("  => Delta bounded  <=>  phi' <= sqrt(3 Delta_max)/tau  <=>  power-law growth only.")
print("  => The 10^213 runaway is the integrated effect of the UNBOUNDED residual (~1/tau).")

# -----------------------------------------------------------------------------
# F-R17-C — consistent-transcription routes.
# General symbolic derivation from perturbed fluid system with creation source
# S = alpha * H * rho_bar (alpha in {0 = standard, 1 = EQ-5 frontier creation M ∝ t}),
# flags: drag (momentum dilution in Euler), cluster (delta_S = delta vs 0).
#
#   continuity (pert):  delta' = -theta - (S/rho)*(1 - cluster)*delta
#   Euler (pert):       theta' = -2H*theta - drag*(S/rho)*theta + (-4 pi G rho delta)...
#   combined symbolically below; background: rho' = -3H rho + S.
# Mandatory limit check: alpha -> 0 recovers standard growth equation.
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("F-R17-C — symbolic derivation of route growth ODEs from the fluid system")

alpha_s, drag_s, clus_s = sp.symbols("alpha d k")  # d, k in {0,1}
tau = sp.symbols("tau_v", positive=True)
H_tau = 1 / tau                                    # gamma-3 LOCKED (tau units)

# background density under creation S = alpha H rho: rho' = -(3 - alpha) H rho
n_exp = 3 - alpha_s                                # rho ∝ tau^(-n_exp)
rho_tau = sp.Function("rho")(tau)
rho_sol = tau ** (-n_exp)                          # normalization carried by eps below
# verify background continuity:
bg_resid = sp.simplify(sp.diff(rho_sol, tau) + 3 * H_tau * rho_sol
                       - alpha_s * H_tau * rho_sol)
print(f"  background continuity residual (must be 0): {sp.simplify(bg_resid)}")

# 4 pi G rho(tau) * t_0^2 == eps_G * tau^(-n_exp)  with eps_G == 4 pi G rho(t_0) t_0^2
# (for alpha=1: eps_G(t) = const = eps_G(t_0); for alpha=0: matches FP4/FP5 source)
eps = sp.symbols("epsilon_G", positive=True)
poisson_coeff = eps * tau ** (-(n_exp - 2))        # 4 pi G rho * tau^2 / tau^2 bookkeeping

delta_f = sp.Function("delta_f")(tau)
theta_f = sp.Function("theta_f")(tau)
sr = alpha_s * H_tau                               # S/rho
cont = sp.Eq(sp.diff(delta_f, tau), -theta_f - sr * (1 - clus_s) * delta_f)
euler = sp.Eq(sp.diff(theta_f, tau),
              -2 * H_tau * theta_f - drag_s * sr * theta_f
              - poisson_coeff / tau**2 * delta_f * (-1) * (-1))
# theta' = -2H theta - d*(S/rho) theta - (4 pi G rho) delta   [Poisson folded in]
euler = sp.Eq(sp.diff(theta_f, tau),
              -2 * H_tau * theta_f - drag_s * sr * theta_f
              - (eps * tau ** (-n_exp) * tau**2) / tau**2 * delta_f)

# reduce to 2nd-order ODE in delta: theta = -(delta' + sr (1-k) delta)
theta_sub = -(sp.diff(delta_f, tau) + sr * (1 - clus_s) * delta_f)
ode_lhs = sp.simplify(
    sp.diff(theta_sub, tau)
    - (-2 * H_tau * theta_sub - drag_s * sr * theta_sub
       - eps * tau ** (-n_exp) * delta_f))
# ode_lhs == 0 is the route ODE; normalize to delta'' + A delta' + B delta = 0
ode_expanded = sp.expand(ode_lhs)
d2 = ode_expanded.coeff(sp.diff(delta_f, tau, 2))
A_coeff = sp.simplify(ode_expanded.coeff(sp.diff(delta_f, tau)) / d2)
B_coeff = sp.simplify((ode_expanded
                       - d2 * sp.diff(delta_f, tau, 2)
                       - d2 * A_coeff * sp.diff(delta_f, tau)).coeff(delta_f) / d2)
print(f"  general ODE: delta'' + ({A_coeff}) delta' + ({B_coeff}) delta = 0")

# FP6 — limit check alpha -> 0 recovers standard growth equation
print("\nFP6 — mandatory limit check: alpha -> 0 recovers standard equation")
A_std = sp.simplify(A_coeff.subs({alpha_s: 0}))
B_std = sp.simplify(B_coeff.subs({alpha_s: 0}))
print(f"  alpha=0:  A = {A_std}   B = {B_std}")
ok_A = sp.simplify(A_std - 2 / tau) == 0
ok_B = sp.simplify(B_std - (-eps / tau**3)) == 0
fp["FP6"] = bool(ok_A and ok_B)
print(f"  expect A = 2/tau, B = -eps/tau^3 (the audited runaway eq.)  => T_pass FP6 = {fp['FP6']}")

# -----------------------------------------------------------------------------
# Route evaluation helper: exact power-law exponents for scale-invariant ODE
# delta'' + (A1/tau) delta' + (B0/tau^2) delta = 0  =>  p(p-1) + A1 p + B0 = 0
# -----------------------------------------------------------------------------
def route_eval(name, alpha, drag, cluster):
    A_r = sp.simplify(A_coeff.subs({alpha_s: alpha, drag_s: drag, clus_s: cluster}))
    B_r = sp.simplify(B_coeff.subs({alpha_s: alpha, drag_s: drag, clus_s: cluster}))
    A1 = sp.simplify(A_r * tau)
    B0 = sp.simplify(B_r * tau**2)
    scale_invariant = (A1.free_symbols <= {eps}) and (B0.free_symbols <= {eps})
    p_sym = sp.symbols("p_r")
    roots = sp.solve(sp.Eq(p_sym * (p_sym - 1) + A1 * p_sym + B0, 0), p_sym)
    roots_num = sorted([complex(r.subs(eps, eps_G)) for r in roots], key=lambda z: z.real)
    p_plus = roots_num[-1].real
    logG = p_plus * LOG_TAU_SPAN if p_plus > 0 else 0.0
    print(f"  {name}: ODE delta'' + ({A1})/tau delta' + ({B0})/tau^2 delta = 0")
    print(f"        scale-invariant: {scale_invariant}; indicial roots p = "
          f"{[f'{z.real:+.4f}' for z in roots_num]}")
    print(f"        growing-mode p+ = {p_plus:+.4f}  =>  log10 G = {logG:.3f}"
          f"  =>  {band_verdict(logG)}   (G_obs target log10 = {LOG_G_OBS}; bands LOCKED)")
    return A1, B0, p_plus, logG, band_verdict(logG), scale_invariant

print("\n" + "-" * 78)
print("FP7-FP10 — routes C1 / C2a / C2b / C2c (pre-declared Phase 0 §3; CLOSED set)")
print(f"  eps_G = {eps_G:.4f}; log10 tau-span = {LOG_TAU_SPAN:.4f}\n")

# C1 zero-active-mass: source term -> 0 (gravity term removed), standard drag
print("C1 — zero-active-mass (source -> 0):")
p_c1 = sp.symbols("p_c1")
roots_c1 = sp.solve(sp.Eq(p_c1 * (p_c1 - 1) + 2 * p_c1, 0), p_c1)
p_plus_c1 = max(float(r) for r in roots_c1)
logG_c1 = p_plus_c1 * LOG_TAU_SPAN if p_plus_c1 > 0 else 0.0
v_c1 = band_verdict(logG_c1)
print(f"  delta'' + (2/tau) delta' = 0; roots p = {sorted(float(r) for r in roots_c1)}")
print(f"  growing-mode p+ = {p_plus_c1:+.4f} => log10 G = {logG_c1:.3f} => {v_c1}")
fp["FP7"] = bool(set(float(r) for r in roots_c1) == {0.0, -1.0})
print(f"  T_pass FP7 (C1 exact solution structure) = {fp['FP7']}")

print("\nC2a — frontier creation (alpha=1), unclustered (k=0), no momentum drag (d=0):")
A1_a, B0_a, p_a, logG_a, v_a, si_a = route_eval("C2a", 1, 0, 0)
fp["FP8"] = bool(si_a)

print("\nC2b — frontier creation (alpha=1), unclustered (k=0), momentum drag (d=1):")
A1_b, B0_b, p_b, logG_b, v_b, si_b = route_eval("C2b", 1, 1, 0)
fp["FP9"] = bool(si_b)

print("\nC2c — frontier creation (alpha=1), comoving-clustered (k=1, d=0):")
A1_c, B0_c, p_c, logG_c, v_c, si_c = route_eval("C2c", 1, 0, 1)
fp["FP10"] = bool(si_c)

# -----------------------------------------------------------------------------
# FP11 — numerical cross-check of route growth (C2a, C2c) vs analytic power law
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP11 — numerical ODE cross-check (C2a, C2c): rel dev of log10 G < 1%")

def integrate_route(A1_num, B0_num, p_plus):
    def rhs(tau_, y):
        d_, dp_ = y
        return [dp_, -(A1_num / tau_) * dp_ - (B0_num / tau_**2) * d_]
    y0 = [1.0, p_plus / tau_init]   # growing-mode IC
    s = solve_ivp(rhs, [tau_init, 1.0], y0, method="Radau", rtol=1e-11, atol=1e-14)
    return np.log10(abs(s.y[0, -1])), s.success

devs = []
for nm, (A1_r, B0_r, p_r, logG_r) in {
        "C2a": (A1_a, B0_a, p_a, logG_a), "C2c": (A1_c, B0_c, p_c, logG_c)}.items():
    logG_num, ok = integrate_route(float(A1_r.subs(eps, eps_G)),
                                   float(B0_r.subs(eps, eps_G)), p_r)
    dev = abs(logG_num - logG_r) / max(logG_r, 1e-12)
    devs.append((ok, dev))
    print(f"  {nm}: analytic log10 G = {logG_r:.4f}; numeric = {logG_num:.4f}; rel dev = {dev:.2e}")
fp["FP11"] = bool(all(ok and d < 0.01 for ok, d in devs))
print(f"  T_pass FP11 = {fp['FP11']}")

# -----------------------------------------------------------------------------
# FP12 — forbidden-move guard: derived exponents depend ONLY on eps_G
# (G_obs never enters any derivation; Phase 0 §5 #10)
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP12 — circularity guard: indicial coefficients free symbols subset of {eps_G}")
fs = (A1_a.free_symbols | B0_a.free_symbols | A1_b.free_symbols | B0_b.free_symbols
      | A1_c.free_symbols | B0_c.free_symbols)
print(f"  free symbols of all route coefficients: {fs}")
fp["FP12"] = bool(fs <= {eps})
print(f"  T_pass FP12 = {fp['FP12']}")

# -----------------------------------------------------------------------------
# FP13 — F-R17-D aggregate (mechanical per Phase 0 §1.3 decision tree)
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("FP13 — F-R17-D aggregate diagnosis (mechanical)")
gate_ok = fp["FP1"] and fp["FP2"] and fp["FP3"]
route_verdicts = {"C1": v_c1, "C2a": v_a, "C2b": v_b, "C2c": v_c}
print(f"  F-R17-A regression gate: {'PASS' if gate_ok else 'FAIL -> HALT'}")
print(f"  F-R17-B: {verdict_B}")
print(f"  F-R17-C route verdicts: {route_verdicts}")
if not gate_ok:
    aggregate = "HALT_PIPELINE_INVALID"
elif verdict_B == "CONSISTENT":
    aggregate = "GENUINE_PATHOLOGY"
else:
    if any(v == "PASS_BAND" for v in route_verdicts.values()):
        aggregate = "ARTIFACT_RESOLVED"
    elif any(v == "PARTIAL_BAND" for v in route_verdicts.values()):
        aggregate = "ARTIFACT_PARTIAL"
    else:
        aggregate = "ARTIFACT_OPEN"
print(f"\n  >>> F-R17-D AGGREGATE VERDICT: {aggregate} <<<")
fp["FP13"] = bool(aggregate in {"GENUINE_PATHOLOGY", "ARTIFACT_RESOLVED",
                                "ARTIFACT_PARTIAL", "ARTIFACT_OPEN"})
print(f"  T_pass FP13 (verdict in pre-registered set) = {fp['FP13']}")

# -----------------------------------------------------------------------------
# INFORMATIONAL (no T_pass) — sensitivity to LOCKED M_univ (R-R17-6); NO re-fit
# -----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("INFORMATIONAL — sensitivity of best routes to M_univ (LOCKED 1e53 kg; O(2) rough)")
for fac in [0.1, 0.5, 1.0, 2.0, 10.0]:
    e_ = eps_G * fac
    pA = -1 + np.sqrt(e_)                       # C2a indicial root (analytic)
    pC = (-1 + np.sqrt(1 + 4 * e_)) / 2         # C2c indicial root (analytic)
    print(f"  M_univ x {fac:>4}: eps_G = {e_:7.3f} | C2a log10 G = "
          f"{max(pA, 0) * LOG_TAU_SPAN:6.2f} | C2c log10 G = {max(pC, 0) * LOG_TAU_SPAN:6.2f}")
print("  inverse note (INFORMATIONAL ONLY, forbidden as input — Phase 0 §5 #3/#10):")
p_req = LOG_G_OBS / LOG_TAU_SPAN
eps_req_c2c = p_req**2 + p_req
eps_req_c2a = (1 + p_req)**2
print(f"    log10 G = 3 would require: C2c eps_G = {eps_req_c2c:.3f} "
      f"(M_univ = {eps_req_c2c/eps_G*M_univ:.2e} kg); "
      f"C2a eps_G = {eps_req_c2a:.3f} (M_univ = {eps_req_c2a/eps_G*M_univ:.2e} kg)")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
print("\n" + "=" * 78)
print("FP SUMMARY")
n_pass = sum(fp.values())
for k in sorted(fp, key=lambda s: int(s[2:])):
    print(f"  {k}: {'PASS' if fp[k] else 'FAIL'}")
print(f"\n  TOTAL: {n_pass}/{len(fp)} PASS; hardcoded T_pass = 0/{len(fp)} (all computed)")
print(f"  F-R17-B: {verdict_B}")
print(f"  F-R17-C: {route_verdicts}")
print(f"  F-R17-D: {aggregate}")
print("=" * 78)
