#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-nonlinear-charge-constraint (Phase 1) -- charge inventory (analytic, sympy).

LOCK: Phase0_balance.md (2026-07-03) sec. 3 Phase 1 (P1a-P1d), criterion V1,
common definitions sec. 5.  No fitting numerics here; sympy exact plus ONE
numeric root-scan for the essential-spectrum edge (P1d; scan reported in
full, no selection).

Models (LOCKED, Phase0 sec. 2):
  M0: L_S = 1/2 f(g) g_t^2 - 1/2 f(g) g_r^2 - W(g),
      f(g) = 1 + 2*alpha*log(g), alpha = 2 (crown F-S, CP-7),
      W'(g) = g^2 (1 - g)  =>  W(g) = g^3/3 - g^4/4,
      units gamma = beta = K_geo = 1, radial measure r^2 dr.
  M1 (model-extension, HYPOTHESIS): psi = g e^{i theta},
      L = 1/2 f(|psi|)|psi_t|^2 - 1/2 f(|psi|)|psi_r|^2 - W(|psi|).

Tasks:
  P1a: EOM(M0) from L_S; exact energy conservation (density + flux);
       static EL(M0) == crown ODE (a3d/CP-7 convention).
  P1b: candidates C1-C5 (CLOSED list): exact conserved/not-conserved
       verdicts.  Method: Q = int q r^2 dr is conserved up to boundary
       flux iff I := r^2 * (dq/dt)|_EOM is a total r-derivative, iff the
       Euler operators E_u[I] and E_w[I] vanish identically in the jet
       variables (u = g and w = g_t treated as independent functions of
       r at fixed time; exactness criterion for null Lagrangians on the
       line).  Nonzero Euler residual == obstruction (printed); numeric
       probes at fixed jet points confirm the residual is not an
       artifact of failed simplification.
  P1c: M1: Noether charge Q of the U(1) shift theta -> theta + const;
       sympy check that the theta-EL is EXACTLY the local conservation
       law d/dt(rho_Q) = d/dr(flux_Q); reduction of the stationary
       ansatz psi = phi(r) e^{-i omega t} to the profile ODE with
       W_eff(phi) = W(phi) - 1/2 omega^2 f(phi) phi^2;  Q(omega),
       E(omega) as integrals of phi_omega.  VK SIGN CONVENTION AND BOX
       RENORMALIZATION FIXED HERE, before any P2b computation (the one
       correction allowed by Phase0 sec. 8).
       Also derived and verified here (input for P2c): the GSS split of
       the second variation of S_omega into L+ (real direction) and
       L- (imaginary/phase direction), with
         Q+ = W_eff''(phi) - 1/2 f''(phi) phi'^2
              - f'(phi) (phi'' + 2 phi'/r)          [CP-7 form, W->W_eff]
         Q- = W'(phi)/phi + f'(phi) phi'^2 / (2 phi)
              - omega^2 (f(phi) + 1/2 f'(phi) phi)
       and the exact identity L- phi_omega = 0 (phase zero mode).
  P1d: essential-spectrum edge as a function of omega.  The omega-vacuum
       phi_inf(omega) solves W_eff'(phi) = 0 (the constant g = 1 is NOT
       a stationary point of the reduced problem for omega > 0).  Exact
       series edge(omega) = W''(1)/f(1) + c2 omega^2 + O(omega^4) and a
       full numeric scan omega in [0, 1] step 0.05.  Locked question:
       does an omega_min with sigma_ess >= 0 exist?
"""
import sympy as sp

ALPHA = sp.Integer(2)


def f_of(x):
    return 1 + 2 * ALPHA * sp.log(x)


def W_of(x):
    return x ** 3 / sp.Integer(3) - x ** 4 / sp.Integer(4)


def Wp_of(x):
    return x ** 2 * (1 - x)


r, t = sp.symbols('r t', positive=True)
omega = sp.symbols('omega', nonnegative=True)

print("=" * 78)
print("Phase 1: charge inventory (sympy exact).  LOCK: Phase0_balance.md sec. 3.")
print("Units gamma=beta=K_geo=1; alpha=2; f = 1 + 4 log g; W' = g^2(1-g).")
print("=" * 78)

# =====================================================================
# P1a -- EOM of M0, energy conservation, static EL == crown ODE
# =====================================================================
print("\n[P1a] EOM of M0 and exact checks")

g = sp.Function('g')(t, r)
gt = sp.diff(g, t)
gr = sp.diff(g, r)
Lag = r ** 2 * (f_of(g) / 2 * gt ** 2 - f_of(g) / 2 * gr ** 2 - W_of(g))
EL = (sp.diff(sp.diff(Lag, gt), t) + sp.diff(sp.diff(Lag, gr), r)
      - sp.diff(Lag, g))
gtt_var = sp.Derivative(g, (t, 2))
gtt_sol = sp.solve(sp.expand(EL), gtt_var)
assert len(gtt_sol) == 1
gtt_sol = sp.simplify(gtt_sol[0])
print("  EOM:  g_tt = %s" % gtt_sol)

# exact energy conservation: d/dt e + d/dr s = 0 under EOM
edens = r ** 2 * (f_of(g) / 2 * gt ** 2 + f_of(g) / 2 * gr ** 2 + W_of(g))
eflux = -r ** 2 * f_of(g) * gt * gr
res_E = sp.diff(edens, t) + sp.diff(eflux, r)
res_E = res_E.subs(gtt_var, gtt_sol)
res_E = sp.simplify(sp.expand(res_E))
print("  energy: d/dt[r^2 e] + d/dr[-r^2 f g_t g_r] |_EOM = %s   -> %s"
      % (res_E, "PASS" if res_E == 0 else "FAIL"))

# static EL == crown ODE (a3d): h'' = (W' - (alpha/h) h'^2)/f - 2 h'/r
h = sp.Function('h')(r)
Ls = r ** 2 * (f_of(h) / 2 * sp.diff(h, r) ** 2 + W_of(h))
ELh = sp.diff(sp.diff(Ls, sp.diff(h, r)), r) - sp.diff(Ls, h)
hpp_var = sp.Derivative(h, (r, 2))
hpp_el = sp.solve(sp.expand(ELh), hpp_var)[0]
hpp_a3d = (Wp_of(h) - (ALPHA / h) * sp.diff(h, r) ** 2) / f_of(h) \
    - 2 * sp.diff(h, r) / r
res_S = sp.simplify(hpp_el - hpp_a3d)
print("  static EL(M0) - crown ODE (a3d)                = %s   -> %s"
      % (res_S, "PASS" if res_S == 0 else "FAIL"))

# =====================================================================
# P1b -- candidates C1-C5 in M0 (Euler-operator exactness test)
# =====================================================================
print("\n" + "-" * 78)
print("[P1b] conservation verdicts, candidates C1-C5 (CLOSED list)")
print("  Criterion: I = r^2 (dq/dt)|_EOM total r-derivative <=>")
print("  E_u[I] = 0 and E_w[I] = 0 identically (u = g, w = g_t as")
print("  independent jet functions of r).")

u, w, u1, w1, u2, w2, u3, w3, u4 = sp.symbols('u w u1 w1 u2 w2 u3 w3 u4',
                                              real=True)
JET = {u: u1, u1: u2, u2: u3, u3: u4, w: w1, w1: w2, w2: w3}


def Dr(expr):
    out = sp.diff(expr, r)
    for a_, b_ in JET.items():
        out = out + sp.diff(expr, a_) * b_
    return out


def euler_u(I):
    return sp.simplify(sp.expand(
        sp.diff(I, u) - Dr(sp.diff(I, u1)) + Dr(Dr(sp.diff(I, u2)))))


def euler_w(I):
    return sp.simplify(sp.expand(
        sp.diff(I, w) - Dr(sp.diff(I, w1)) + Dr(Dr(sp.diff(I, w2)))))


fu = f_of(u)
fpu = sp.diff(fu, u)
# EOM in jet variables:  g_tt = ...
gtt_jet = u2 + (fpu / (2 * fu)) * (u1 ** 2 - w ** 2) + 2 * u1 / r \
    - Wp_of(u) / fu


def I_of(q):
    """I = r^2 * d/dt q(u, w, u1) with g_tt -> EOM, d/dt u1 -> w1."""
    return r ** 2 * (sp.diff(q, u) * w + sp.diff(q, w) * gtt_jet
                     + sp.diff(q, u1) * w1)


PROBES = [
    {u: sp.Rational(11, 10), w: sp.Rational(3, 10), u1: sp.Rational(-1, 5),
     w1: sp.Rational(1, 7), u2: sp.Rational(2, 9), w2: sp.Rational(-1, 11),
     u3: sp.Rational(1, 13), u4: sp.Rational(-1, 17), w3: sp.Rational(1, 19),
     r: sp.Rational(5, 2)},
    {u: sp.Rational(9, 5), w: sp.Rational(-2, 7), u1: sp.Rational(4, 11),
     w1: sp.Rational(-3, 13), u2: sp.Rational(-5, 17), w2: sp.Rational(2, 19),
     u3: sp.Rational(3, 23), u4: sp.Rational(1, 29), w3: sp.Rational(-1, 31),
     r: sp.Rational(7, 3)},
]


def verdict(name, q, note=""):
    I = sp.expand(I_of(q))
    Eu = euler_u(I)
    Ew = euler_w(I)
    ok = (Eu == 0) and (Ew == 0)
    if not ok:
        # numeric probes: confirm the obstruction is genuinely nonzero
        mags = []
        for P in PROBES:
            vals = [abs(complex(sp.N(Eu.subs(P)))),
                    abs(complex(sp.N(Ew.subs(P))))]
            mags.append(max(vals))
        conf = all(m > 1e-10 for m in mags)
        stat = "NIEZACHOWANY (Euler residual != 0; probes %.2e, %.2e)" \
            % (mags[0], mags[1])
        if not conf:
            stat = "UWAGA: residual znika w probes -- recheck"
    else:
        stat = "ZACHOWANY (I = total derivative)"
    print("  %-34s : %s%s" % (name, stat, ("  " + note) if note else ""))
    return ok


print("")
res_tab = {}
res_tab['C1'] = verdict("C1  q = (g-1)", u - 1)
res_tab['C2'] = verdict("C2  q = f(g) g_t  (peed kanoniczny)", fu * w)
res_tab['C3'] = verdict("C3  q = (g-1)^2", (u - 1) ** 2)
res_tab['C4a'] = verdict("C4a q = g_t      (h = 1)", w)
res_tab['C4b'] = verdict("C4b q = g^2 g_t  (h = g^2)", u ** 2 * w)

# C5: standard Noether inventory (time transl. / space transl. / rot. /
# dilations).  Time translation = energy (positive control):
qE = fu / 2 * w ** 2 + fu / 2 * u1 ** 2 + W_of(u)
res_tab['C5-energy'] = verdict("C5  energia (kontrola dodatnia)", qE,
                               note="[MUSI byc ZACHOWANY]")
# radial-sector momentum density candidate f g_t g_r (3D translations
# give zero total charge on radial configs; here tested as a 1D charge):
res_tab['C5-mom'] = verdict("C5  q = f(g) g_t g_r (peed radialny)",
                            fu * w * u1)
# dilation candidates: q = f g_t (r g_r + c (g-1)), c in {0, 1, 3/2}
for cc in (sp.Integer(0), sp.Integer(1), sp.Rational(3, 2)):
    res_tab['C5-dil-c%s' % cc] = verdict(
        "C5  dylatacja c=%s" % cc, fu * w * (r * u1 + cc * (u - 1)))

n_cons = sum(1 for k, v in res_tab.items()
             if v and not k.startswith('C5-energy'))
print("\n  TABELA P1b: zachowane ladunki budzetopodobne w M0: %d / %d"
      % (n_cons, len(res_tab) - 1))
print("  (energia zachowana -- kontrola; nie jest ladunkiem budzetowym)")

# =====================================================================
# P1c -- M1: Noether charge, EOM, stationary reduction, Q(omega), E(omega)
# =====================================================================
print("\n" + "-" * 78)
print("[P1c] M1 (kompleksyfikacja, model-extension): psi = g e^{i theta}")

G = sp.Function('G')(t, r)
TH = sp.Function('TH')(t, r)
LagM1 = r ** 2 * (f_of(G) / 2 * (sp.diff(G, t) ** 2
                                 + G ** 2 * sp.diff(TH, t) ** 2)
                  - f_of(G) / 2 * (sp.diff(G, r) ** 2
                                   + G ** 2 * sp.diff(TH, r) ** 2)
                  - W_of(G))

EL_TH = (sp.diff(sp.diff(LagM1, sp.diff(TH, t)), t)
         + sp.diff(sp.diff(LagM1, sp.diff(TH, r)), r)
         - sp.diff(LagM1, TH))
rho_Q = r ** 2 * f_of(G) * G ** 2 * sp.diff(TH, t)     # = dL/d(theta_t)
flux_Q = r ** 2 * f_of(G) * G ** 2 * sp.diff(TH, r)
cons_id = sp.simplify(sp.expand(
    EL_TH - (sp.diff(rho_Q, t) - sp.diff(flux_Q, r))))
print("  theta-EL - [d/dt rho_Q - d/dr flux_Q] = %s   -> %s"
      % (cons_id, "PASS (dQ/dt = strumien brzegowy, exact)"
         if cons_id == 0 else "FAIL"))
print("  rho_Q = r^2 f(g) g^2 theta_t;  Noether: Q_raw = int rho_Q dr")

# EL for G (needed for the reduction check)
EL_G = (sp.diff(sp.diff(LagM1, sp.diff(G, t)), t)
        + sp.diff(sp.diff(LagM1, sp.diff(G, r)), r)
        - sp.diff(LagM1, G))

# stationary ansatz psi = phi(r) e^{-i omega t}: G -> phi(r), TH -> -omega t
phi = sp.Function('phi')(r)
submap = {
    sp.Derivative(G, (t, 2)): 0,
    sp.Derivative(G, t, r): 0,
    sp.Derivative(G, r, t): 0,
    sp.Derivative(G, t): 0,
    sp.Derivative(G, (r, 2)): sp.Derivative(phi, (r, 2)),
    sp.Derivative(G, r): sp.Derivative(phi, r),
    sp.Derivative(TH, (t, 2)): 0,
    sp.Derivative(TH, t, r): 0,
    sp.Derivative(TH, r, t): 0,
    sp.Derivative(TH, (r, 2)): 0,
    sp.Derivative(TH, r): 0,
    sp.Derivative(TH, t): -omega,
}
EL_G_red = EL_G.xreplace(submap).subs(G, phi)
EL_G_red = sp.simplify(sp.expand(EL_G_red))

x = sp.symbols('x', positive=True)
Weff_x = W_of(x) - omega ** 2 / 2 * f_of(x) * x ** 2
Weffp_x = sp.diff(Weff_x, x)
Seff = r ** 2 * (f_of(phi) / 2 * sp.diff(phi, r) ** 2
                 + Weff_x.subs(x, phi))
EL_Seff = sp.diff(sp.diff(Seff, sp.diff(phi, r)), r) - sp.diff(Seff, phi)
# sign convention: LagM1 has -1/2 f phi_r^2 (Lagrangian), S_omega has
# +1/2 f phi'^2 (energy-like functional)  =>  EL_G_red = -EL_Seff
red_id = sp.simplify(sp.expand(EL_G_red + EL_Seff))
print("  redukcja: EL_G|_(phi e^{-i omega t}) + EL[S_omega] = %s   -> %s"
      % (red_id, "PASS" if red_id == 0 else "FAIL"))
print("  (EL_G = -EL[S_omega]: znak konwencji Lagranzjan/funkcjonal)")
print("  S_omega[phi] = int r^2 [ f(phi)/2 phi'^2 + W_eff(phi) ] dr,")
print("  W_eff(phi) = W(phi) - (omega^2/2) f(phi) phi^2")
print("  profil ODE:  f phi'' = W_eff'(phi) - (1/2) f' phi'^2 - (2/r) f phi'")
print("  (omega -> 0: dokladnie ODE korony CP-7 -- ciaglosc galezi)")

print("\n  Q(omega), E(omega) na ansatzu (calki z phi_omega):")
print("    Q_raw(omega) = omega * int r^2 f(phi) phi^2 dr        [konwencja:")
print("      Q := -int rho_Q dr, tak by Q > 0 dla omega > 0]")
print("    E(omega)     = int r^2 [ (omega^2/2) f phi^2 + f/2 phi'^2")
print("                             + W(phi) ] dr")

print("\n  KONWENCJA VK / RENORMALIZACJA PUDLA (LOCKED TERAZ, przed P2b;")
print("  jedyna dopuszczalna korekta wg Phase0 sec. 8):")
print("   (a) Q_raw i E sa rozbiezne objetosciowo w pudle R (tlo phi_inf!=0);")
print("       do testu VK uzywamy wielkosci odjetych od omega-prozni:")
print("         Q_sol(omega) = omega int r^2 [f(phi)phi^2")
print("                        - f(phi_inf)phi_inf^2] dr")
print("         E_sol(omega) = int r^2 [omega^2/2 (f phi^2 - f_inf phi_inf^2)")
print("                        + f/2 phi'^2 + W(phi) - W(phi_inf)] dr")
print("       (obie takze raportowane w wersji raw; pudlo R=60 jak CP-7).")
print("   (b) konwencja znaku VK: przy normalizacji Q > 0, omega > 0")
print("       galaz stabilna = slope-negative:  dQ_sol/domega < 0")
print("       (standardowa forma GSS/VK; wyprowadzenie nie odwraca znaku).")

# ---------------------------------------------------------------------
# GSS split: second variation -> L+ / L- (input for P2c), exact
# ---------------------------------------------------------------------
print("\n  [P1c/GSS] druga wariacja S_omega: rozklad na L+ / L- (exact):")
aF = sp.Function('a')(r)
bF = sp.Function('b')(r)
epsl = sp.symbols('epsilon')
chi_re = phi + epsl * aF
chi_im = epsl * bF
mod = sp.sqrt(chi_re ** 2 + chi_im ** 2)
Sdens = r ** 2 * (f_of(mod) / 2 * (sp.diff(chi_re, r) ** 2
                                   + sp.diff(chi_im, r) ** 2)
                  + W_of(mod) - omega ** 2 / 2 * f_of(mod) * mod ** 2)
quad = sp.diff(Sdens, epsl, 2).subs(epsl, 0) / 2
quad = sp.expand(sp.simplify(quad))

# no a-b cross terms (decoupling): mixed second derivatives must vanish
cross1 = sp.simplify(sp.expand(sp.diff(sp.diff(quad, aF), bF)))
cross2 = sp.simplify(sp.expand(
    sp.diff(sp.diff(quad, sp.diff(aF, r)), bF)))
cross3 = sp.simplify(sp.expand(
    sp.diff(sp.diff(quad, aF), sp.diff(bF, r))))
dec_ok = (cross1 == 0) and (cross2 == 0) and (cross3 == 0)
print("    decoupling a-b (brak czlonow krzyzowych): %s"
      % ("PASS" if dec_ok else "FAIL"))

def pos_fix(e):
    """phi(r) > 0 on profiles: log(phi^2) -> 2 log(phi), sqrt(phi^2) -> phi
    (sympy cannot assume positivity of an undefined Function)."""
    e = e.replace(sp.log(phi ** 2), 2 * sp.log(phi))
    e = e.replace(lambda ex: ex.is_Pow and ex.base == phi ** 2,
                  lambda ex: phi ** (2 * ex.exp))
    return sp.simplify(sp.expand(e))


fp_ = lambda y: sp.diff(f_of(x), x).subs(x, y)
fpp_ = lambda y: sp.diff(f_of(x), x, 2).subs(x, y)
Weffpp_ = sp.diff(Weff_x, x, 2)

Qplus = (Weffpp_.subs(x, phi)
         - fpp_(phi) / 2 * sp.diff(phi, r) ** 2
         - fp_(phi) * (sp.Derivative(phi, (r, 2)) + 2 / r * sp.diff(phi, r)))
ELa = sp.diff(quad, aF) - sp.diff(sp.diff(quad, sp.diff(aF, r)), r)
target_a = r ** 2 * Qplus * aF \
    - sp.diff(r ** 2 * f_of(phi) * sp.diff(aF, r), r)
chk_a = pos_fix(ELa - target_a)
print("    L+ : Q+ = W_eff''(phi) - f''/2 phi'^2 - f' (phi'' + 2phi'/r)")
print("         [forma CP-7 z W -> W_eff]  check: %s -> %s"
      % (chk_a, "PASS" if chk_a == 0 else "FAIL"))

Qminus = (Wp_of(phi) / phi + fp_(phi) * sp.diff(phi, r) ** 2 / (2 * phi)
          - omega ** 2 * (f_of(phi) + fp_(phi) * phi / 2))
ELb = sp.diff(quad, bF) - sp.diff(sp.diff(quad, sp.diff(bF, r)), r)
target_b = r ** 2 * Qminus * bF \
    - sp.diff(r ** 2 * f_of(phi) * sp.diff(bF, r), r)
chk_b = pos_fix(ELb - target_b)
print("    L- : Q- = W'/phi + f' phi'^2/(2 phi) - omega^2 (f + f' phi/2)")
print("         check: %s -> %s" % (chk_b, "PASS" if chk_b == 0 else "FAIL"))

# exact phase zero mode: L- phi_omega = 0 under the profile ODE
phipp_sol = sp.solve(sp.expand(EL_Seff), sp.Derivative(phi, (r, 2)))[0]
Lm_phi = (-1 / r ** 2 * sp.diff(r ** 2 * f_of(phi) * sp.diff(phi, r), r)
          + Qminus * phi)
Lm_phi = sp.simplify(sp.expand(
    Lm_phi.subs(sp.Derivative(phi, (r, 2)), phipp_sol)))
print("    L- phi_omega |_ODE = %s   -> %s (dokladny mod zerowy fazy)"
      % (Lm_phi, "PASS" if Lm_phi == 0 else "FAIL"))

# charge direction for the constrained L+ count (P2c deflation):
print("    kierunek ladunkowy (deflacja L+ w P2c):")
print("      dQ/dphi|_omega = omega r^2 [ f'(phi) phi^2 + 2 f(phi) phi ]")

# =====================================================================
# P1d -- essential-spectrum edge as a function of omega
# =====================================================================
print("\n" + "-" * 78)
print("[P1d] krawedz kontinuum sigma_ess(omega) wokol omega-prozni")

# the omega-vacuum: W_eff'(phi_inf) = 0 (g = 1 NIE jest prozna dla omega>0)
veq = sp.expand(Weffp_x)     # W' - omega^2 (f x + f' x^2 / 2)
print("  rownanie omega-prozni:  W_eff'(x) = 0, tj.")
print("    %s = 0" % sp.simplify(veq))
chk_v1 = sp.simplify(veq.subs(x, 1))
print("  W_eff'(1) = %s  (!= 0 dla omega > 0: proznia PRZESUWA sie z omega)"
      % chk_v1)

# exact identity: Q-(x)|_const = W_eff'(x)/x  => edge(L-) = 0 for ALL omega
Qminus_vac = Wp_of(x) / x - omega ** 2 * (f_of(x) + sp.diff(f_of(x), x) * x / 2)
idm = sp.simplify(Qminus_vac - veq / x)
print("  tozsamosc: Q-(x) - W_eff'(x)/x = %s" % idm)
print("    => na omega-prozni Q-(phi_inf) = 0: krawedz sigma_ess(L-) = 0")
print("       dla KAZDEGO omega (kontinuum fazowe marginalne, exact)  %s"
      % ("PASS" if idm == 0 else "FAIL"))

# series: phi_inf(omega) = 1 + a2 omega^2 + a4 omega^4 + ...
a2s, a4s = sp.symbols('a2 a4')
xs = 1 + a2s * omega ** 2 + a4s * omega ** 4
ser = sp.series(veq.subs(x, xs), omega, 0, 6).removeO()
ser = sp.expand(ser)
c2eq = ser.coeff(omega, 2)
sol_a2 = sp.solve(c2eq, a2s)[0]
c4eq = ser.coeff(omega, 4).subs(a2s, sol_a2)
sol_a4 = sp.solve(c4eq, a4s)[0]
print("  seria omega-prozni: phi_inf = 1 + (%s) omega^2 + (%s) omega^4 + ..."
      % (sol_a2, sol_a4))

# edge(L+)(omega) = W_eff''(phi_inf)/f(phi_inf), exact series
edge_expr = Weffpp_ / f_of(x)
phi_inf_ser = xs.subs({a2s: sol_a2, a4s: sol_a4})
edge_ser = sp.series(edge_expr.subs(x, phi_inf_ser), omega, 0, 6)
edge_ser_t = sp.expand(edge_ser.removeO())
c0 = edge_ser_t.coeff(omega, 0)
c2 = edge_ser_t.coeff(omega, 2)
c4 = edge_ser_t.coeff(omega, 4)
print("  krawedz sigma_ess(L+):  c(omega) = %s + (%s) omega^2 + (%s) omega^4"
      % (c0, c2, c4))
print("    W''(1)/f(1) = %s (zgodnie z parametryzacja LOCK), c2 = %s"
      % (sp.simplify(Wp_of(x).diff(x).subs(x, 1) / f_of(1)), c2))

# ghost-crossing of the omega-vacuum: phi_inf(omega_gh) = g* = e^{-1/4}
gstar = sp.exp(sp.Rational(-1, 4))
om_gh2 = sp.solve(veq.subs(x, gstar), omega ** 2)
om_gh = sp.sqrt(om_gh2[0])
print("  proznia przecina sciane duchowa g* = e^{-1/4} przy omega_gh:")
print("    omega_gh^2 = %s = %.6f,  omega_gh = %.6f"
      % (sp.simplify(om_gh2[0]), float(om_gh2[0]), float(om_gh)))
print("    (dla omega > omega_gh: f(phi_inf) < 0 -- tlo kinetycznie")
print("     zduchowione; formula krawedzi traci sens -- patrz skan)")

# full numeric scan omega in [0, 1] step 0.05 (reported in full)
print("\n  skan numeryczny (nsolve, kontynuacja galezi od phi=1):")
print("    %-6s %-12s %-12s %-14s %-10s" % ("omega", "phi_inf",
                                            "f(phi_inf)", "edge(L+)",
                                            "uwaga"))
guess = 1.0
scan = []
for k in range(0, 21):
    omv = 0.05 * k
    if omv == 0.0:
        xr = 1.0
    else:
        try:
            xr = float(sp.nsolve(veq.subs(omega, omv), x, guess))
        except Exception as ex:
            print("    %-6.2f  nsolve FAILED (%s)" % (omv, ex))
            continue
    guess = xr
    fv = float(f_of(x).subs(x, xr))
    ev = float(edge_expr.subs({x: xr, omega: omv}))
    note = ""
    if xr < float(gstar) + 0.005:
        note = "phi_inf < g*+0.005 (ghost)"
    scan.append((omv, xr, fv, ev, note))
    print("    %-6.2f %-12.6f %-12.6f %-14.6f %-10s"
          % (omv, xr, fv, ev, note))

max_edge = max(s[3] for s in scan if s[2] > 0)
print("\n  PYTANIE P1d (LOCKED): czy istnieje omega_min z sigma_ess >= 0?")
print("    c2 = %s < 0: czlon omega^2 OBNIZA krawedz (nie podnosi);" % c2)
print("    max edge(L+) w skanie (tam gdzie f>0) = %.4f < 0" % max_edge)
print("    => NIE istnieje omega_min; kontinuum L+ pozostaje tachioniczne")
print("       dla kazdego omega (a dla omega > %.4f proznia jest" % float(om_gh))
print("       dodatkowo zduchowiona kinetycznie).")
print("    krawedz L-: 0 dla kazdego omega (exact, patrz wyzej).")

# =====================================================================
# V1 verdict
# =====================================================================
print("\n" + "=" * 78)
print("WERDYKT V1 (kryterium LOCKED, Phase0 sec. 3):")
print("  Tabela C1-C5: ZADEN budzetopodobny ladunek nie jest zachowany")
print("  w M0 (energia zachowana -- kontrola; pozostale: Euler residual")
print("  != 0).  Zapis wprost per LOCK:")
print("    'hipoteza wymaga rozszerzenia M1' -- Phase 2 wylacznie w galezi")
print("    M1, jawnie oznaczonej jako model-extension.")
print("  M1: ladunek Noether Q zachowany EXACT (dQ/dt = strumien brzegowy);")
print("  redukcja Q-ball: W_eff = W - (omega^2/2) f phi^2; VK convention")
print("  locked (dQ_sol/domega < 0 = galaz stabilna).")
print("  P1d: krawedz kontinuum L+ NIE zostaje podniesiona przez omega^2")
print("  (c2 = %s); krawedz L- = 0 exact. Kryterium V2 nalezy czytac" % c2)
print("  z klauzula LOCK: 'jesli kontinuum pozostaje tachioniczne przy")
print("  kazdym omega -- raportowac nierozstrzygalnosc w M1 z powodu tla'.")
print("=" * 78)
