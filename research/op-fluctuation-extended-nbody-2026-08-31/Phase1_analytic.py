# -*- coding: utf-8 -*-
"""
Phase 1 (analitycznie, sympy) — op-fluctuation-extended-nbody-2026-08-31
LOCK: Phase0_balance.md §3, §4 (QN-1), §8.

P1-1: det macierzy korelacji 3x3 = 1 - sum g_ij^2 + 2 g12 g13 g23 (exact).
P1-2: Delta F3 do rzedu 4: g12 g13 g23 - (1/2)(g12^2 g13^2 + g12^2 g23^2
      + g13^2 g23^2) + O(g^5) (sympy series, exact).
P1-3: znak czlonu wiodacego dla g_ij > 0: DODATNI (3-cialo OSLABIA
      przyciaganie parowe) — probka numeryczna + granica symboliczna.
P1-4: tozsamosc monopolowa (rank-one): dla K = c*J (J = macierz jedynek)
      det(I - SigA^{-1} K SigB^{-1} K^T) = 1 - c^2 * C_A * C_B,
      C_X = 1^T SigX^{-1} 1 ("pojemnosc" klastra) — sympy exact (n=2).
P1-5: forma dalekiego pola na krytycznosci: F ~ -1/2 C_A C_B (4 pi d)^{-2}
      => d ln|F| / d ln d = -2 NIEZALEZNIE od R (R tylko w amplitudzie).
P1-6 (deskryptywnie): wiodace Delta F3 dla konfiguracji T i C przy g=c/d.
"""
import sympy as sp

results = {}
def check(name, cond, detail=""):
    results[name] = bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")

print("=" * 72)
print("Phase 1 — analityka (sympy exact)")
print("=" * 72)

g12, g13, g23, eps = sp.symbols('g12 g13 g23 epsilon', positive=True)

# ------------------------------------------------------------------ P1-1
Sig3 = sp.Matrix([[1, g12, g13], [g12, 1, g23], [g13, g23, 1]])
det3 = sp.expand(Sig3.det())
det3_expected = 1 - g12**2 - g13**2 - g23**2 + 2 * g12 * g13 * g23
check("P1-1: det Sigma3/G0^3 = 1 - sum g^2 + 2 g12 g13 g23",
      sp.simplify(det3 - det3_expected) == 0, f"det = {det3}")

# ------------------------------------------------------------------ P1-2
# Delta F3 = 1/2 [ ln det3 - sum_par ln(1 - g_ij^2) ]; skala eps -> szereg
subs_eps = {g12: eps * g12, g13: eps * g13, g23: eps * g23}
dF3 = sp.Rational(1, 2) * (sp.log(det3.subs(subs_eps))
                           - sp.log(1 - (eps * g12)**2)
                           - sp.log(1 - (eps * g13)**2)
                           - sp.log(1 - (eps * g23)**2))
series = sp.series(dF3, eps, 0, 5).removeO()
series = sp.expand(series)
expected = (eps**3 * g12 * g13 * g23
            - sp.Rational(1, 2) * eps**4 * (g12**2 * g13**2
                                            + g12**2 * g23**2
                                            + g13**2 * g23**2))
check("P1-2: Delta F3 = g12 g13 g23 - 1/2 sum g^2 g^2 + O(g^5)",
      sp.simplify(series - expected) == 0,
      f"szereg(eps^3..4) = {series}")

# ------------------------------------------------------------------ P1-3
# znak: czlon wiodacy g12 g13 g23 > 0 dla g_ij > 0 (G_m > 0) => Delta F3 > 0
# przy malych g: 3-cialo DODATNIE = OSLABIA laczne przyciaganie parowe.
val = dF3.subs({eps: 1, g12: sp.Rational(1, 10), g13: sp.Rational(1, 10),
                g23: sp.Rational(1, 12)})
check("P1-3: Delta F3 > 0 (probka exact g=(0.1,0.1,1/12))",
      sp.simplify(val) > 0, f"Delta F3 = {sp.N(val, 8)}")

# ------------------------------------------------------------------ P1-4
# rank-one (monopol): K = c*J_{2x2}; SigA, SigB symetryczne dodatnie 2x2
a0, a1, ar, b0, b1, br, c = sp.symbols('a0 a1 a_r b0 b1 b_r c')
SigA = sp.Matrix([[a0, ar], [ar, a1]])
SigB = sp.Matrix([[b0, br], [br, b1]])
J = sp.ones(2, 2)
M = sp.eye(2) - c**2 * SigA.inv() * J * SigB.inv() * J
one = sp.ones(2, 1)
CA = (one.T * SigA.inv() * one)[0, 0]
CB = (one.T * SigB.inv() * one)[0, 0]
lhs = sp.simplify(M.det())
rhs = sp.simplify(1 - c**2 * CA * CB)
check("P1-4: det(I - c^2 SigA^-1 J SigB^-1 J) = 1 - c^2 C_A C_B (n=2 exact)",
      sp.simplify(lhs - rhs) == 0)

# ------------------------------------------------------------------ P1-5
# daleko-polowa forma krytyczna: F(d) = -1/2 C_A C_B / (4 pi d)^2
d_, CAs, CBs = sp.symbols('d C_A C_B', positive=True)
F_far = -sp.Rational(1, 2) * CAs * CBs / (4 * sp.pi * d_)**2
# d ln|F| / d ln d = d * d/dd ln|F|   (incydent P1-5a: pierwotna proba
# sp.diff(., sp.log(d_)) to blad skladniowy sympy — poprawka maszynerii,
# kryterium bez zmian; patrz Phase_FINAL_close)
slope2 = sp.simplify(d_ * sp.diff(sp.log(-F_far), d_))
check("P1-5: wykladnik dalekiego pola = -2 NIEZALEZNIE od C_A, C_B (od R)",
      sp.simplify(slope2 + 2) == 0, f"d ln|F|/d ln d = {slope2}")

# ------------------------------------------------------------------ P1-6
# deskryptywnie: krytycznie g(d) = kappa/d;
# T: (d,d,d*sqrt2) -> DF3 ~ kappa^3/(sqrt2 d^3);
# C: (d,d,2d)      -> DF3 ~ kappa^3/(2 d^3)
kap = sp.symbols('kappa', positive=True)
DF3_T = (kap / d_) * (kap / d_) * (kap / (sp.sqrt(2) * d_))
DF3_C = (kap / d_) * (kap / d_) * (kap / (2 * d_))
print(f"    P1-6 (deskryptywnie): DF3_T = {sp.simplify(DF3_T)}, "
      f"DF3_C = {sp.simplify(DF3_C)}, T/C = {sp.simplify(DF3_T/DF3_C)}")
print("    obie konfiguracje: DF3 ~ 1/d^3 (wolniej o d^-1 od czlonu"
      " parowego 1/d^2 * g^2 ~ 1/d^4... raport w Phase 3 numerycznie)")

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 1: {npass}/{ntot} PASS")
for k, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k}")
if npass == ntot:
    print("WERDYKT QN-1: TAK — czlon wiodacy Delta F3 = g12*g13*g23 (rzad 3),")
    print("  DODATNI dla g_ij>0: nieaddytywnosc 3-cialowa OSLABIA przyciaganie")
    print("  parowe (znak okreslony symbolicznie). Monopol: wykladnik -2,")
    print("  R tylko w amplitudzie (C_A C_B) — hipoteza do testu w Phase 2.")
