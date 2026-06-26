# -*- coding: ascii -*-
# Phase 2 sympy -- op-frontier-bridge-and-asymmetry-2026-06-12
# Modul A: F-BA-2 (GAP-2 bottom-up J_source + regulator) and F-BA-3 (GAP-3 A-iv monochromaticity)
#
# Pre-registration: Phase0_balance.md LOCKED 2026-06-12 (par. 1.3 F-BA-2/F-BA-3; par. 8(e),(i) pre-derivations)
# LOCKED inputs (read-only; forbidden move #4): Mdot = 2c^3/(9G); v_c = 2c/3; u = (2/3)x/t;
#   rho_bar = 1/(6 pi G t^2); marginality (1/2)v_c^2 = GM/(ct) (FCR P3 trichotomy);
#   wall ledger: DeltaV = lam Phi0^4/4, delta = 2/m_Phi, stability edge |Phi| = Phi0/sqrt(3) (FM P1);
#   supply = DeltaV * 4 pi R^2 c, t_* (FM P3); trajectory exponents {2/3, 1/3} (FM P2 FP7).
# Declared approximations: thin-wall (delta << ct); Newtonian frontier bookkeeping (inherited LOCKED).
# 0 hardcoded T_pass; circularity guard FP11.

import sympy as sp

OUT = []
def log(s=""):
    OUT.append(str(s))
    print(s)

RESULTS = []
def verdict(name, desc, computed_pass, detail):
    status = "PASS" if computed_pass else "FAIL"
    RESULTS.append((name, status))
    log("[%s] %s -- %s" % (status, name, desc))
    log("        %s" % detail)

log("=" * 78)
log("Phase 2 sympy -- op-frontier-bridge-and-asymmetry (modul A: GAP-2 + GAP-3)")
log("=" * 78)

t, c, G, lam, Phi0, x, rho_e, eps = sp.symbols('t c G lambda Phi0 x rho_e epsilon', positive=True)
mPhi = sp.sqrt(2 * lam) * Phi0
delta = 2 / mPhi
v_c = sp.Rational(2, 3) * c                       # LOCKED FM P2
R = c * t                                          # LOCKED gamma-3
rho_bar = 1 / (6 * sp.pi * G * t**2)               # LOCKED FM (COR-1)
Mdot_top = 2 * c**3 / (9 * G)                      # LOCKED FCR P3 + FM (read-only)
DeltaV = lam * Phi0**4 / 4                         # LOCKED FM P1 ledger

# ================================================================ F-BA-2 (GAP-2)
log("")
log("--- F-BA-2: bottom-up J_source + regulator (GAP-2) ---")

# FP1: bottom-up kinematic functional: J = rho_e (Rdot - v_c) per area; Mdot = 4 pi R^2 J
J_per_area = rho_bar * (sp.diff(R, t) - v_c)
Mdot_bu = sp.simplify(4 * sp.pi * R**2 * J_per_area)
fp1_ok = sp.simplify(Mdot_bu - Mdot_top) == 0
verdict("FP1", "bottom-up Mdot = 4 pi R^2 rho_e (Rdot - v_c) == 2c^3/(9G) top-down EXACT",
        fp1_ok, "Mdot_bottom_up = %s ; difference vs LOCKED top-down = %s"
        % (Mdot_bu, sp.simplify(Mdot_bu - Mdot_top)))

# FP2: regulator: marginality at the entry boundary FORCES rho_e
# (1/2)v_c^2 = G M / (c t), M = (4 pi/3) rho_e (c t)^3 -> unique rho_e
M_of_rho = sp.Rational(4, 3) * sp.pi * rho_e * (c * t)**3
marg_eq = sp.Eq(sp.Rational(1, 2) * v_c**2, G * M_of_rho / (c * t))
rho_sol = sp.solve(marg_eq, rho_e)
fp2_ok = (len(rho_sol) == 1) and (sp.simplify(rho_sol[0] - rho_bar) == 0)
verdict("FP2", "marginality boundary condition has UNIQUE solution rho_e = 1/(6 pi G t^2) = rho_bar",
        fp2_ok, "solutions: %s ; matches LOCKED rho_bar: %s"
        % (rho_sol, sp.simplify(rho_sol[0] - rho_bar) == 0))

# FP3: trichotomy directions (regulator mechanism): cost(rho_e) = (1/2)v_c^2 - G M(rho_e)/(ct)
cost = sp.Rational(1, 2) * v_c**2 - G * M_of_rho / (c * t)
dcost = sp.simplify(sp.diff(cost, rho_e))
cost_at_bar = sp.simplify(cost.subs(rho_e, rho_bar))
# over-deposition (rho_e > rho_bar) -> cost < 0 (runaway branch, excluded by LOCKED trichotomy)
# under-deposition -> cost > 0 (blocked branch, excluded by M ~ t)
cost_over = sp.simplify(cost.subs(rho_e, 2 * rho_bar))
cost_under = sp.simplify(cost.subs(rho_e, rho_bar / 2))
fp3_ok = (dcost.is_negative is True or sp.simplify(dcost + sp.Rational(4,3)*sp.pi*G*c**2*t**2) == 0) \
         and cost_at_bar == 0 and (cost_over.is_negative is True) and (cost_under.is_positive is True)
verdict("FP3", "trichotomy regulator: cost(rho_bar) = 0 EXACT; over-deposition -> runaway branch (cost<0); under -> blocked (cost>0)",
        fp3_ok, "d(cost)/d(rho_e) = %s < 0 ; cost(rho_bar) = %s ; cost(2 rho_bar) = %s ; cost(rho_bar/2) = %s"
        % (dcost, cost_at_bar, cost_over, cost_under))

# FP4: conversion efficiency eta(t) = demand/supply = (t_*/t)^2 ; t_* == FM P3 LOCKED form
supply = DeltaV * 4 * sp.pi * R**2 * c
eta = sp.simplify((Mdot_top * c**2) / supply)
tstar_FM = (sp.sqrt(2) / (3 * sp.sqrt(sp.pi))) * c / (sp.sqrt(G * lam) * Phi0**2)
eta_form_ok = sp.simplify(eta - (tstar_FM / t)**2) == 0
eta_at_tstar = sp.simplify(eta.subs(t, tstar_FM))
surplus_limit = sp.limit(1 - eta, t, sp.oo)
fp4_ok = eta_form_ok and (eta_at_tstar == 1) and (surplus_limit == 1)
verdict("FP4", "efficiency eta(t) = (t_*/t)^2 EXACT (t_* == FM P3 LOCKED); eta(t_*) = 1; surplus -> wall kinetics",
        fp4_ok, "eta = %s ; eta == (t_*/t)^2: %s ; eta(t_*) = %s ; lim(1-eta) = %s"
        % (eta, eta_form_ok, eta_at_tstar, surplus_limit))

# FP5: time-exponent audit: capacity ~ t^2 regulated to rate ~ t^0 (pre-registered DERIVED condition)
exp_rate = sp.simplify(sp.diff(sp.log(Mdot_bu), sp.log(t))) if Mdot_bu.has(t) else 0
exp_capacity = sp.simplify(t * sp.diff(supply, t) / supply)
fp5_ok = (exp_rate == 0) and (sp.simplify(exp_capacity - 2) == 0)
verdict("FP5", "exponent audit: d ln(Mdot_bu)/d ln t = 0 (regulated) vs capacity exponent = 2",
        fp5_ok, "rate exponent = %s ; capacity exponent = %s" % (exp_rate, exp_capacity))

# ================================================================ F-BA-3 (GAP-3)
log("")
log("--- F-BA-3: A-iv monochromaticity from wall dynamics (GAP-3) ---")

# FP6: entry surface unique: |Phi| = Phi0/sqrt3 crossed once (monotone profile); depth = delta*atanh(1/sqrt3)
prof = Phi0 * sp.tanh(mPhi * x / 2)
xs = sp.solve(sp.Eq(prof, Phi0 / sp.sqrt(3)), x)
x_star = delta * sp.atanh(1 / sp.sqrt(3))
prof_monotone = sp.simplify(sp.diff(prof, x))
# identity check via exp-map (atanh/log rewrite plumbing, same class as Phase 1 par. 5.6):
# exp(m_Phi * x) for both forms must equal 2 + sqrt(3) exactly
exp_solve = sp.radsimp(sp.simplify(sp.powsimp(sp.exp(xs[0] * mPhi), force=True)))
exp_star = sp.radsimp(sp.simplify(sp.exp(x_star * mPhi).rewrite(sp.log)))
ident_solve = sp.simplify(sp.radsimp(exp_solve - (2 + sp.sqrt(3)))) == 0
ident_star = sp.simplify(sp.radsimp(exp_star - (2 + sp.sqrt(3)))) == 0
fp6_ok = (len(xs) == 1) and ident_solve and ident_star and (prof_monotone.is_positive is True)
verdict("FP6", "entry surface UNIQUE: single crossing of |Phi| = Phi0/sqrt3 (monotone profile); depth x* = delta*atanh(1/sqrt3)",
        fp6_ok, "crossings: %d ; exp(m_Phi x): solved = %s, atanh-form = %s (both == 2+sqrt3: %s/%s) ; dPhi/dx > 0: %s"
        % (len(xs), exp_solve, exp_star, ident_solve, ident_star, prof_monotone.is_positive is True))

# FP7: geometric dispersion: deposition spread across layer of thickness <= delta
u_flow = sp.Rational(2, 3) * x / t                 # LOCKED FCR P2 derived flow
dv_geom = sp.simplify(u_flow.subs(x, c * t) - u_flow.subs(x, c * t - delta))
ratio_geom = sp.simplify(dv_geom / v_c)
# exponent in (delta/(c t)): d ln(ratio)/d ln(delta)
w = sp.symbols('w', positive=True)                 # w = delta/(c t) bookkeeping variable
ratio_in_w = sp.simplify(ratio_geom.subs(delta, w * c * t))
exponent_geom = sp.simplify(w * sp.diff(ratio_in_w, w) / ratio_in_w)
thinwall_limit = sp.limit(ratio_in_w, w, 0)
fp7_ok = (sp.simplify(ratio_geom - delta / (c * t)) == 0) and (exponent_geom == 1) and (thinwall_limit == 0)
verdict("FP7", "geometric entry dispersion: dv/v_c = delta/(ct) EXACT -- exponent 1 (pre-registered >= 1); -> 0 thin-wall",
        fp7_ok, "dv/v_c = %s ; exponent in (delta/ct) = %s ; thin-wall limit = %s"
        % (ratio_geom, exponent_geom, thinwall_limit))

# FP8: temporal dispersion: formation time delta/c; flow changes by |du/dt|*(delta/c) at x = ct
dv_temp = sp.simplify(sp.Abs(sp.diff(u_flow, t)).subs(x, c * t) * (delta / c))
ratio_temp = sp.simplify(dv_temp / v_c)
fp8_ok = sp.simplify(ratio_temp - delta / (c * t)) == 0
verdict("FP8", "temporal entry dispersion: |du/dt|*(delta/c)/v_c = delta/(ct) -- same order, exponent 1",
        fp8_ok, "dv_t/v_c = %s" % ratio_temp)

# FP9: residual (recoil) dispersion DECAYS: dv_pec/dt = -(du/dx) v_pec => v_pec ~ t^(-2/3) = 1/a_m
vp = sp.Function('v_pec')
ode = sp.Eq(sp.Derivative(vp(t), t), -sp.diff(u_flow, x) * vp(t))
sol = sp.dsolve(ode, vp(t))
vpec_t = sol.rhs
a_m = t**sp.Rational(2, 3)                         # LOCKED FCR P2: a_m ~ t^(2/3)
invariant = sp.simplify(vpec_t * a_m)
fp9_ok = (not invariant.has(t))
# consistency with LOCKED FM P2 FP7 exponents {2/3, 1/3}: v_pec ~ d/dt t^(1/3) ~ t^(-2/3) -- same law
verdict("FP9", "recoil/peculiar dispersion decays: v_pec ~ 1/a_m ~ t^(-2/3) (consistent with LOCKED FM exponents {2/3,1/3})",
        fp9_ok, "v_pec(t) = %s ; v_pec * a_m = %s (t-free: %s)" % (vpec_t, invariant, fp9_ok))

# FP10: dust attractor: effective pressure term decays faster than density
# P_eff ~ rho <v_pec^2> ~ t^(-2) * t^(-4/3); w_eff = P/(rho c^2) ~ t^(-4/3) -> 0
w_eff = sp.simplify((vpec_t / c)**2)
w_eff_norm = sp.simplify(w_eff / w_eff.subs(t, 1))
w_exponent = sp.simplify(t * sp.diff(w_eff_norm, t) / w_eff_norm)
w_limit = sp.limit(w_eff_norm, t, sp.oo)
fp10_ok = (sp.simplify(w_exponent + sp.Rational(4, 3)) == 0) and (w_limit == 0)
verdict("FP10", "dust attractor: w_eff = <v_pec^2>/c^2 ~ t^(-4/3) -> 0 (growth form protected asymptotically)",
        fp10_ok, "w_eff exponent = %s (expected -4/3) ; limit t->oo = %s" % (w_exponent, w_limit))

# ================================================================ FP11: circularity guard
blacklist = {'G_obs', 'eta_B', 'etaB', 'f_max', 'fbar_max', 'asym_obs', 'log10G_obs'}
free = set()
for expr in [Mdot_bu, rho_sol[0], cost, eta, ratio_geom, ratio_temp, vpec_t, w_eff, x_star, supply]:
    free |= {str(fs) for fs in expr.free_symbols}
hits = free & blacklist
fp11_ok = (len(hits) == 0)
verdict("FP11", "circularity guard: no G_obs / eta_B / observed-asymmetry symbol in any verdict expression",
        fp11_ok, "free symbols audited: %s ; blacklist hits: %s" % (sorted(free), sorted(hits)))

# ================================================================ summary
log("=" * 78)
npass = sum(1 for _, st in RESULTS if st == "PASS")
log("SUMMARY: %d/%d PASS ; 0 hardcoded T_pass (all verdicts computed)" % (npass, len(RESULTS)))
log("F-BA-2 (GAP-2): bottom-up functional J = rho_e (Rdot - v_c); regulator = LOCKED marginality")
log("  trichotomy applied at entry boundary (unique rho_e = rho_bar; over->runaway, under->blocked);")
log("  Mdot_bu == 2c^3/(9G) EXACT; exponent t^0 vs capacity t^2; eta = (t_*/t)^2 EXACT.")
log("F-BA-3 (GAP-3): geometric + temporal entry dispersion = delta/(ct) EXACT (exponent 1, -> 0);")
log("  recoil dispersion NOT derived from wall dynamics (creation kinematics unknown -- declared);")
log("  but decays as 1/a_m and w_eff ~ t^(-4/3) -> 0 (dust attractor protects growth machinery).")

with open(__file__.replace('.py', '.txt'), 'w') as f:
    f.write("\n".join(OUT) + "\n")
