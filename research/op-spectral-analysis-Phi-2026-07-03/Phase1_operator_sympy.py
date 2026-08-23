#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-spectral-analysis-Phi (CP-7 / L03) -- Phase 1: EXACT fluctuation operator.

Derives (sympy, no hand-waving) the second variation of
    E[u] = int [ (1/2) F(u) |grad u|^2 + W(u) ] d^3x
around a radial background u0(r), reduces it to Sturm-Liouville form,
and specializes to the three formulations locked in Phase0_balance.md:

  F-A  (canonical grav, alpha=2):  F = K_geo * psi^4,  U_A' = K_geo*gamma*(psi^7 - psi^6)
  F-S  (crown soliton log-form) :  F = 1 + 4*ln g,     W'   = g^2 (1-g)
  F-S' (substrate, alpha=1)     :  F = K_geo * g^2,    W'   = g^2 (1-g)

Tests (LOCKED in Phase 0):
  T1  EL equation of E_A == eq:full-field (thm:field-eq, alpha=2, beta=gamma)
  T2  Q-form: Q = W''(u0) - (1/2) F''(u0) u0'^2 - F'(u0) [u0'' + (2/r) u0']
      derived independently via delta^2 E + integration by parts == locked form
  T3  C1: vacuum psi=1 in F-A: m_sp^2 = gamma  (exact)
  T4  EL of E_S == a3d ODE:  f*(g'' + 2 g'/r) + (2/g) g'^2 = g^2 (1-g)
  T5  vacuum g=1 in F-S: continuum edge = W''(1)/f(1) = -1 (exact)  [dualism L04, measured]
  T6  vacuum g=1 in F-S': edge = W''(1)/(K_geo*1^2) = -1/K_geo      [same sign]
  T7  ghost-wall F-A endpoint psi->0: canonical variable chi = psi^3/3,
      effective potential -> classification data (limit point?)

All results printed; PASS/FAIL per test. ASCII only.
"""
import sympy as sp

PASS = 0
FAIL = 0
def test(name, cond, detail=""):
    global PASS, FAIL
    ok = bool(cond)
    PASS += ok
    FAIL += (not ok)
    print(("  [OK] " if ok else "  [FAIL] ") + name + ((" :: " + detail) if detail else ""))

r = sp.symbols('r', positive=True)
eps = sp.symbols('epsilon')
gamma, beta, Kgeo, q, rho, lam = sp.symbols('gamma beta K_geo q rho lambda', positive=True)

u0 = sp.Function('u0')(r)      # background
v  = sp.Function('v')(r)       # fluctuation
F  = sp.Function('F')
W  = sp.Function('W')

print("=" * 72)
print("Phase 1: exact second variation  E[u] = int (1/2)F|grad u|^2 + W(u)")
print("=" * 72)

# ----------------------------------------------------------------------
# Generic energy density (radial, weight r^2), expand to O(eps^2)
# ----------------------------------------------------------------------
u = u0 + eps * v
dens = (sp.Rational(1, 2) * F(u) * sp.diff(u, r)**2 + W(u)) * r**2
dens2 = sp.diff(dens, eps, 2).subs(eps, 0) / 2          # eps^2 coefficient
dens2 = sp.expand(dens2)

# Integrate the cross term F'(u0) u0' v v' by parts:
#   r^2 F' u0' v v' = d/dr[ (1/2) r^2 F' u0' v^2 ] - (1/2) d/dr(r^2 F' u0') v^2
cross_coeff = sp.Rational(1, 1) * r**2 * sp.diff(F(u0), u0.func(r)) if False else None
# do it explicitly on the expression:
Fp  = sp.Derivative(F(u0), u0)  # placeholder; use diff via subs function
# Simpler: rebuild expected pieces and compare coefficient-wise.
u0p = sp.diff(u0, r)
u0pp = sp.diff(u0, r, 2)
vp = sp.diff(v, r)

F0   = F(u0)
F1_  = sp.diff(F(u0), r) / u0p          # F'(u0) (chain rule, formal)
F1_  = sp.simplify(F1_)
# safer: define F', F'' as derivatives of F at u0 using sp.Derivative with .doit through dummy
x = sp.symbols('x')
Fx = F(x); Wx = W(x)
F1 = Fx.diff(x).subs(x, u0)
F2 = Fx.diff(x, 2).subs(x, u0)
W2 = Wx.diff(x, 2).subs(x, u0)

expected_dens2 = (sp.Rational(1, 2) * F0 * vp**2
                  + F1 * u0p * v * vp
                  + sp.Rational(1, 2) * (W2 + sp.Rational(1, 2) * F2 * u0p**2) * v**2) * r**2

test("T2a: eps^2 density == (1/2)F v'^2 + F' u0' v v' + (1/2)[W'' + (1/2)F'' u0'^2] v^2",
     sp.simplify(dens2 - expected_dens2) == 0)

# Integration by parts of the cross term:
#   int r^2 F' u0' v v' dr = -(1/2) int d/dr(r^2 F' u0') v^2 dr   (+ boundary)
ibp_extra = -sp.Rational(1, 2) * sp.diff(r**2 * F1 * u0p, r)      # coefficient of v^2
Qtimesr2 = sp.expand(sp.Rational(1, 2) * (W2 + sp.Rational(1, 2) * F2 * u0p**2) * r**2 * 2
                     + 2 * ibp_extra) / 2 * 2
# Effective quadratic form: int [ (1/2) r^2 F v'^2 + (1/2) Qeff r^2 v^2 ] dr with
Qeff = sp.simplify((W2 + sp.Rational(1, 2) * F2 * u0p**2) - sp.diff(r**2 * F1 * u0p, r) / r**2)
Q_locked = W2 - sp.Rational(1, 2) * F2 * u0p**2 - F1 * (u0pp + 2 * u0p / r)
test("T2b: Q == W''(u0) - (1/2)F''(u0) u0'^2 - F'(u0)[u0'' + (2/r)u0']  (locked form)",
     sp.simplify(Qeff - Q_locked) == 0,
     "operator: L[v] = -(1/r^2)(r^2 F v')' + Q v")

# ----------------------------------------------------------------------
# F-A: canonical gravitational formulation (alpha=2)
# ----------------------------------------------------------------------
print("-" * 72)
print("F-A: F = K_geo psi^4, U_A' = K_geo gamma (psi^7 - psi^6)   [beta=gamma]")
psi = sp.Function('psi')(r)
FA = Kgeo * x**4
UAp = Kgeo * (gamma * x**7 - beta * x**6)            # U_A'(x)
UA = sp.integrate(UAp, x)

# EL equation: -(1/r^2)(r^2 F psi')' + (1/2)F' psi'^2 + U_A'(psi) = 0  (source-free)
psip = sp.diff(psi, r)
psipp = sp.diff(psi, r, 2)
EL_A = (-(1 / r**2) * sp.diff(r**2 * FA.subs(x, psi) * psip, r)
        + sp.Rational(1, 2) * FA.diff(x).subs(x, psi) * psip**2
        + UAp.subs(x, psi))
# eq:full-field (alpha=2, source-free), multiplied by -K_geo psi^4:
eq_full_field = psipp + 2 * psip / r + 2 * psip**2 / psi + beta * psi**2 - gamma * psi**3
test("T1: EL(E_A) == -K_geo psi^4 * [eq:full-field(alpha=2)]",
     sp.simplify(EL_A + Kgeo * psi**4 * eq_full_field) == 0,
     "canonical action reproduces thm:field-eq exactly")

# C1: vacuum psi=1 (beta=gamma): mass gap
UA2 = UAp.diff(x)
m2_vac = sp.simplify((UA2.subs(x, 1) / FA.subs(x, 1)).subs(beta, gamma))
test("T3 (C1): F-A vacuum m_sp^2 = U_A''(1)/K(1) = gamma", sp.simplify(m2_vac - gamma) == 0,
     "m_sp^2 = %s" % m2_vac)

# Q on Yukawa background psi0 = 1 + A exp(-sqrt(gamma) r)/r  (printed for Phase 2 use)
A = sp.symbols('A', positive=True)
psi0 = 1 + A * sp.exp(-sp.sqrt(gamma) * r) / r
FA_1 = FA.diff(x)
FA_2 = FA.diff(x, 2)
Q_A = (UA2.subs(x, psi0)
       - sp.Rational(1, 2) * FA_2.subs(x, psi0) * sp.diff(psi0, r)**2
       - FA_1.subs(x, psi0) * (sp.diff(psi0, r, 2) + 2 * sp.diff(psi0, r) / r))
Q_A = Q_A.subs(beta, gamma)
Q_A_at_vac = sp.limit(Q_A.subs(A, 0), r, sp.oo) if False else sp.simplify(Q_A.subs(A, 0))
test("T3b: Q_A(A=0) == gamma*K_geo (vacuum limit consistency)",
     sp.simplify(Q_A_at_vac - gamma * Kgeo) == 0)

# ----------------------------------------------------------------------
# F-S: crown log-form (alpha=2): f = 1 + 4 ln g, W' = g^2(1-g)
# ----------------------------------------------------------------------
print("-" * 72)
print("F-S: f = 1 + 4 ln g, W' = g^2(1-g)   [crown scripts a3d/ls10, ALPHA=2]")
g = sp.Function('g')(r)
fS = 1 + 4 * sp.log(x)
WSp = x**2 * (1 - x)
gp = sp.diff(g, r)
gpp = sp.diff(g, r, 2)
# EL: -(1/r^2)(r^2 f g')' + (1/2) f' g'^2 + W'(g) = 0
EL_S = (-(1 / r**2) * sp.diff(r**2 * fS.subs(x, g) * gp, r)
        + sp.Rational(1, 2) * fS.diff(x).subs(x, g) * gp**2
        + WSp.subs(x, g))
# a3d ODE: f*(g'' + 2g'/r) + (2/g) g'^2 = g^2 (1-g)
a3d_ode = fS.subs(x, g) * (gpp + 2 * gp / r) + (2 / g) * gp**2 - g**2 * (1 - g)
test("T4: EL(E_S) == -[a3d ODE]", sp.simplify(EL_S + a3d_ode) == 0,
     "crown soliton ODE is EL of E_S with f=1+4ln g, W'=g^2(1-g)")

WS2 = WSp.diff(x)
edge_S = sp.simplify(WS2.subs(x, 1) / fS.subs(x, 1))
test("T5: F-S vacuum edge = W''(1)/f(1) = -1 (dualism L04, to be MEASURED in Phase 2)",
     sp.simplify(edge_S + 1) == 0, "edge = %s" % edge_S)

# ----------------------------------------------------------------------
# F-S': substrate (alpha=1): F = K_geo g^2, same W
# ----------------------------------------------------------------------
print("-" * 72)
print("F-S': F = K_geo g^2, W' = g^2(1-g)   [sek08b cor:alpha1-preferred]")
fSp_ = Kgeo * x**2
edge_Sp = sp.simplify(WS2.subs(x, 1) / fSp_.subs(x, 1))
test("T6: F-S' vacuum edge = -1/K_geo (same negative sign as F-S)",
     sp.simplify(edge_Sp + 1 / Kgeo) == 0, "edge = %s" % edge_Sp)
# EL check vs formulation-dictionary ODE: g'' = (1-g) - (1/g)g'^2 - (2/r)g'  (K_geo=1)
EL_Sp = (-(1 / r**2) * sp.diff(r**2 * (x**2).subs(x, g) * gp, r)
         + sp.Rational(1, 2) * (x**2).diff(x).subs(x, g) * gp**2
         + WSp.subs(x, g))
dict_ode = gpp - ((1 - g) - (1 / g) * gp**2 - 2 * gp / r)
test("T6b: EL(E_S') == -g^2 * [dictionary ODE alpha=1]",
     sp.simplify(EL_Sp + g**2 * dict_ode) == 0)

# ----------------------------------------------------------------------
# T7: ghost-wall / degenerate endpoint psi -> 0 in F-A (K = psi^4)
# ----------------------------------------------------------------------
print("-" * 72)
print("T7: F-A endpoint psi->0 (K=psi^4 degenerate): canonical variable analysis")
# Canonical variable: dchi/dpsi = sqrt(K) = psi^2  => chi = psi^3/3.
chi = sp.symbols('chi', positive=True)
psi_of_chi = (3 * chi) ** sp.Rational(1, 3)
test("T7a: chi(psi)=psi^3/3 invertible, psi(chi) smooth for chi>0",
     sp.simplify(psi_of_chi**3 / 3 - chi) == 0)
# Potential in canonical variable: U_A(psi(chi)); mass-type curvature near chi->0:
UA_g = UA.subs(beta, gamma)
U_chi = UA_g.subs(x, psi_of_chi)
d2U = sp.diff(U_chi, chi, 2)
lead = sp.limit(sp.simplify(d2U * chi**sp.Rational(2, 3)), chi, 0)
print("    U(chi) = %s" % sp.simplify(U_chi))
print("    d2U/dchi2 ~ chi^(-2/3) * (%s)  as chi->0" % lead)
# Sign of the leading curvature coefficient decides attract/repel nature at the
# degenerate endpoint; positive => repulsive barrier => endpoint inaccessible
# (limit-point, no BC needed, spectrum unaffected).
lead_sign_positive = sp.simplify(lead) != 0 and (lead / gamma / Kgeo).is_positive
test("T7b: leading curvature coefficient at chi->0 is POSITIVE (repulsive barrier)",
     bool(lead_sign_positive), "coeff = %s" % lead)

print("=" * 72)
print("Phase 1 result: %d PASS, %d FAIL" % (PASS, FAIL))
print("=" * 72)
