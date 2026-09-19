#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- Phase 1: analityka pre-rejestrowana
(LOCK sec. 2; MD sec. 4). ZERO ewolucji numerycznej.

P1-H1: wyprowadzenie U_mat z literalnego sqrt(-g)*(q/Phi0)*psi*rho
       (eq:L-mat-unified + eq:vol-element-M911, sek08a -- TYLKO ODCZYT);
       gate tozsamosci z lam*rho_hat*psi^2/(4-3psi) oraz pochodnej
       lam*rho_hat*psi*(8-3psi)/(4-3psi)^2; gate tozsamosci prozniowo
       odjetej (bez kancelacji) lam*rho_hat*(psi-1)*(psi+4)/(4-3psi);
       kontrola numeryczna 1e-12 w psi in {0.5, 1, 1.2} [INPUT-MD].
P1-H2: odpowiedz zlinearyzowana (-lap+1) dpsi = -5 lam rho_hat;
       dpsi = -5 lam (G_Yuk * rho_hat); redukcja sferyczna (kernel
       Yukawy) + kwadratura scipy (wzorzec gate'u P3-H1a);
       kontrola rezydualna ODE (deskryptywnie).
P1-H3: fakty brzegowe: lim U_mat przy psi->4/3- = +oo; psi->0+ = 0.

Formy CYTAT (MD sec. 1; dziedziczone op-r3-stationary-states
Phase1_output.txt): M=psi^6/(4-3psi)^2, K=psi^4, U=psi^4/4-psi^3/3,
K_geo=gamma=c0=1 [LOCK].
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

OUT = []


def w(s=""):
    OUT.append(str(s))


BAR = "=" * 78
SEP = "-" * 78

w(BAR)
w("PHASE 1 -- analityka pre-rejestrowana; K_geo=gamma=c0=1 [INPUT]")
w("Zrodlo czlonu materii (TYLKO ODCZYT, core sek08a):")
w("  eq:L-mat-unified: L_mat = -(q/Phi0)*psi*rho ;")
w("  rho = -T^mu_mu/c0^2 >= 0 (L01 formal definition 2026-05-04);")
w("  czynnik psi = konsekwencja elementu objetosci (audit A4")
w("  2026-05-01, Option 2 -- preserve axiom);")
w("  eq:vol-element-M911: sqrt(-g) = c0*psi/(4-3psi).")
w("  rho = rho0 * rho_hat(r), rho_hat = exp(-r^2/18) [LOCK, FROZEN];")
w("  lam := (q*c0/Phi0)*rho0 > 0 (bezwymiarowe).")
w(BAR)

# ---------------------------------------------------------------- P1-H1
w("P1-H1: wyprowadzenie U_mat z zapisu literalnego (sympy)")
psi, r, rp = sp.symbols("psi r rp", positive=True)
q, Phi0, c0, rho0, lam = sp.symbols("q Phi0 c0 rho0 lam", positive=True)
rho_hat = sp.exp(-r**2 / 18)

sqrtg = c0 * psi / (4 - 3 * psi)
L_mat = -(q / Phi0) * psi * (rho0 * rho_hat)
# gestosc potencjalna (Lagranzjan L = T - U): U_mat = -sqrt(-g)*L_mat
U_mat_lit = sp.expand(-sqrtg * L_mat)
w("  U_mat literalnie = -sqrt(-g)*L_mat = " + str(sp.simplify(U_mat_lit)))

lam_def = (q * c0 / Phi0) * rho0
U_mat_exp = lam * rho_hat * psi**2 / (4 - 3 * psi)
g1 = sp.simplify(U_mat_lit - U_mat_exp.subs(lam, lam_def))
w("  gate tozsamosci: simplify(U_mat_lit - lam*rho_hat*psi^2/(4-3psi))"
  " = " + str(g1) + ("  PASS" if g1 == 0 else "  FAIL"))

dU_target = lam * rho_hat * psi * (8 - 3 * psi) / (4 - 3 * psi) ** 2
g2 = sp.simplify(sp.diff(U_mat_exp, psi) - dU_target)
w("  gate pochodnej: simplify(dU_mat/dpsi - lam*rho_hat*psi*(8-3psi)/"
  "(4-3psi)^2) = " + str(g2) + ("  PASS" if g2 == 0 else "  FAIL"))

dUm_vac = lam * rho_hat * (psi - 1) * (psi + 4) / (4 - 3 * psi)
g3 = sp.simplify(U_mat_exp - U_mat_exp.subs(psi, 1) - dUm_vac)
w("  gate tozsamosci prozniowo odjetej (bez kancelacji):")
w("    simplify(U_mat(psi)-U_mat(1) - lam*rho_hat*(psi-1)*(psi+4)/"
  "(4-3psi)) = " + str(g3) + ("  PASS" if g3 == 0 else "  FAIL"))

# kontrola numeryczna 1e-12 (lam=1, rho_hat=1 tj. r=0) [INPUT-MD]
w("  kontrola numeryczna (lam=1, r=0): prog |d| <= 1e-12*max(1,|v|)")


def um_f(p):
    return p * p / (4.0 - 3.0 * p)


def ump_f(p):
    return p * (8.0 - 3.0 * p) / (4.0 - 3.0 * p) ** 2


def dum_f(p):
    return (p - 1.0) * (p + 4.0) / (4.0 - 3.0 * p)


numgate_ok = True
for pv in (0.5, 1.0, 1.2):
    pb = sp.Rational(*np.float64(pv).as_integer_ratio())
    subs = {psi: pb, lam: 1, r: 0}
    for nm, expr, ff in (("Um ", U_mat_exp, um_f),
                         ("Ump", dU_target, ump_f),
                         ("dUm", dUm_vac, dum_f)):
        vs = float(sp.N(expr.subs(subs), 30))
        vf = ff(np.float64(pv))
        d = abs(vs - vf)
        ok = d <= 1e-12 * max(1.0, abs(vs))
        numgate_ok = numgate_ok and ok
        w("    psi=%.1f %s sympy(bin)=%+.16e float=%+.16e |d|=%.2e %s"
          % (pv, nm, vs, vf, d, "PASS" if ok else "FAIL"))

p1h1 = (g1 == 0) and (g2 == 0) and (g3 == 0) and numgate_ok
w("  P1-H1: " + ("PASS" if p1h1 else "FAIL"))
w(SEP)

# ---------------------------------------------------------------- P1-H2
w("P1-H2: odpowiedz zlinearyzowana wokol psi=1 (statycznie)")
eps, lamh = sp.symbols("eps lamh", positive=True)
f = sp.Function("f")
psi_e = 1 + eps * f(r)
Kpsi = psi_e**4
Upr = psi_e**2 * (psi_e - 1)
dUm_dpsi = (lamh * eps) * rho_hat * psi_e * (8 - 3 * psi_e) \
    / (4 - 3 * psi_e) ** 2
RHS = (sp.diff(r**2 * Kpsi * sp.diff(psi_e, r), r) / r**2
       - sp.Rational(1, 2) * 4 * psi_e**3 * sp.diff(psi_e, r) ** 2
       - Upr - dUm_dpsi)
o1 = sp.expand(sp.series(RHS, eps, 0, 2).removeO().coeff(eps, 1))
target = sp.diff(f(r), r, 2) + 2 * sp.diff(f(r), r) / r - f(r) \
    - 5 * lamh * rho_hat
gl = sp.simplify(o1 - target)
w("  EOM statyczne, lam = eps*lamh, psi = 1 + eps*f(r); O(eps):")
w("    " + str(sp.simplify(o1)) + " = 0")
w("  gate: simplify(O(eps) - [f'' + 2f'/r - f - 5*lamh*rho_hat]) = "
  + str(gl) + ("  PASS" if gl == 0 else "  FAIL"))
w("  => (-lap + 1) dpsi = -5*lam*rho_hat  (dpsi=eps*f, lam=eps*lamh)")
w("  => dpsi = -5*lam*(G_Yuk * rho_hat), G_Yuk = exp(-r)/(4 pi r)")
w("  wsp. 5 = dU_mat/dpsi|_{psi=1} / (lam rho_hat) = 1*(8-3)/(4-3)^2")

w("  redukcja sferyczna (kernel Yukawy, wyprowadzenie kontrolowane):")
s = sp.symbols("s", positive=True)
mu_of_s = (r**2 + rp**2 - s**2) / (2 * r * rp)
jac = sp.simplify(sp.diff(mu_of_s, s))
w("    s^2 = r^2+rp^2-2 r rp mu ; dmu/ds = " + str(jac)
  + "  (= -s/(r rp): " + str(sp.simplify(jac + s / (r * rp))) + ")")
s_at_p1 = sp.sqrt(sp.expand((r - rp) ** 2))
s_at_m1 = sp.sqrt(sp.expand((r + rp) ** 2))
w("    s(mu=+1) = " + str(sp.simplify(s_at_p1)) + " = |r-rp| ;"
  + " s(mu=-1) = " + str(sp.simplify(s_at_m1)))
kern = sp.integrate(sp.exp(-s), (s, sp.Abs(r - rp), r + rp)) \
    / (2 * r * rp)
w("    (1/4pi) int dOmega' e^{-s}/s = (1/(2 r rp)) int_{|r-rp|}^{r+rp}"
  " e^{-s} ds")
w("      = " + str(sp.simplify(kern)))
w("    => dpsi(r) = -(5 lam / r) int_0^oo rp rho_hat(rp)"
  " [e^{-|r-rp|} - e^{-(r+rp)}]/2 drp")
w("    => dpsi(0) = -5 lam int_0^oo rp rho_hat(rp) e^{-rp} drp")


def rho_np(t):
    return np.exp(-t * t / 18.0)


RMAX_Q = 60.0  # gorna granica kwadratury [INPUT-MD]


def dpsi_lin(rv, lamv):
    """Wzorzec gate'u P3-H1a: kwadratura scipy (MD sec. 4)."""
    if rv == 0.0:
        val = quad(lambda t: t * rho_np(t) * np.exp(-t), 0.0, RMAX_Q,
                   limit=200)[0]
        return -5.0 * lamv * val

    def integ(t):
        return t * rho_np(t) * (np.exp(-abs(rv - t))
                                - np.exp(-(rv + t))) / 2.0

    val = quad(integ, 0.0, RMAX_Q, points=[rv], limit=400)[0]
    return -5.0 * lamv * val / rv


I0 = quad(lambda t: t * rho_np(t) * np.exp(-t), 0.0, RMAX_Q,
          limit=200)[0]
w("  I0 = int_0^60 rp exp(-rp^2/18 - rp) drp = %.12f" % I0)
w("  PREDYKCJA pre-rejestrowana [LOCK sec. 2]: dpsi < 0 wszedzie,")
w("  ogon ~ e^{-r}/r; glebokosci liniowe dpsi_lin(0) = -5*lam*I0:")
for lv in (0.01, 0.05, 0.2, 0.5):
    w("    lam=%-5g  dpsi_lin(0) = %+.8f" % (lv, -5.0 * lv * I0))

w("  profil dpsi_lin(r)/lam (kwadratura):")
for rv in (0.025, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 40.0):
    w("    r=%-6.3f dpsi_lin/lam = %+.8e" % (rv, dpsi_lin(rv, 1.0)))

# kontrola rezydualna ODE (deskryptywnie)
hg = 0.025
rg = (np.arange(int(50.0 / hg)) + 0.5) * hg
fg = np.array([dpsi_lin(x, 1.0) for x in rg])
lap = np.zeros_like(fg)
lap[1:-1] = ((fg[2:] - 2 * fg[1:-1] + fg[:-2]) / hg**2
             + (fg[2:] - fg[:-2]) / hg / rg[1:-1])
res = -lap[1:-1] + fg[1:-1] + 5.0 * rho_np(rg[1:-1])
scale = float(np.max(np.abs(5.0 * rho_np(rg))))
w("  kontrola rezydualna ODE (siatka h=0.025, r<50, FD 2. rzedu):")
w("    max|(-lap+1)dpsi_lin + 5 rho_hat| / max|5 rho_hat| = %.3e"
  % (float(np.max(np.abs(res))) / scale))
mono = bool(np.all(fg < 0.0))
w("  znak: dpsi_lin(r) < 0 na calej siatce: " + str(mono))
p1h2 = (gl == 0) and mono
w("  P1-H2: " + ("PASS (wyprowadzenie + kwadratura-wzorzec)"
                 if p1h2 else "FAIL"))
w(SEP)

# ---------------------------------------------------------------- P1-H3
w("P1-H3: fakty brzegowe U_mat (sympy limit; lam,rho_hat>0)")
rhos = sp.symbols("rhos", positive=True)
Um = lam * rhos * psi**2 / (4 - 3 * psi)
lim_hi = sp.limit(Um, psi, sp.Rational(4, 3), "-")
lim_lo = sp.limit(Um, psi, 0, "+")
w("  lim_{psi->4/3-} U_mat = " + str(lim_hi)
  + ("  (= +oo: zrodlo ODPYCHA od sufitu tam gdzie rho_hat>0)"
     if lim_hi == sp.oo else "  NIEOCZEKIWANE"))
w("  lim_{psi->0+}  U_mat = " + str(lim_lo)
  + ("  (= 0: podloga BEZ bezposredniej bariery od materii)"
     if lim_lo == 0 else "  NIEOCZEKIWANE"))
w("  pochodna na podlodze: lim_{psi->0+} dU_mat/dpsi = "
  + str(sp.limit(lam * rhos * psi * (8 - 3 * psi) / (4 - 3 * psi) ** 2,
                 psi, 0, "+")) + "  (sila materii znika na podlodze)")
w("  ASYMETRIA pre-rejestrowana [LOCK sec. 2]: stabilizacja")
w("  latwiejsza dla kolapsow gornych (sufit) niz dolnych (podloga)")
w("  -- konfrontacja w Q-H2.")
p1h3 = (lim_hi == sp.oo) and (lim_lo == 0)
w("  P1-H3: " + ("PASS" if p1h3 else "FAIL"))
w(BAR)
w("PHASE 1 PODSUMOWANIE: P1-H1 %s ; P1-H2 %s ; P1-H3 %s"
  % tuple("PASS" if x else "FAIL" for x in (p1h1, p1h2, p1h3)))
w("  U_mat(psi,r) = lam * rho_hat(r) * psi^2/(4-3psi)   [WYPROWADZONE]")
w("  dU_mat/dpsi  = lam * rho_hat(r) * psi(8-3psi)/(4-3psi)^2")
w("  U_mat - U_mat(1) = lam * rho_hat * (psi-1)(psi+4)/(4-3psi)")
w("  dpsi_lin(r) = -(5 lam/r) int rp rho_hat [e^{-|r-rp|}-e^{-(r+rp)}]"
  "/2 drp ; dpsi_lin(0) = -5 lam I0, I0=%.9f" % I0)
w(BAR)

with open("TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14/"
          "Phase1_output.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
print("written", len(OUT), "lines")
