#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-action-audit-spectrum-insert (Phase 2) -- Q-D1: dyspersja prozni
z JEDNEJ akcji (sympy + kontrola float).

LOCK sec. 2 Phase 2; MD sec. 3, 5. Linearyzacja EOM (z Phase 1,
P1b PASS) wokol psi*=1: psi = 1 + eps e^{i(kx - om t)} -- WYPROWADZIC,
nie postulowac; omega^2(k) = [Keff(1) k^2 + Ueff''(1)] / M(1);
m^2 = Ueff''(1)/M(1); c_s^2 = Keff(1)/M(1); znaki M, Keff na (0,4/3).
Odczyt A (K_A = psi^5/(4-3psi)) rownolegle -- odnotowanie, NIE wplywa
na werdykt (PRIMARY = B). Deskryptywnie: wariant ZNAKOWANEGO g^tt
(MD sec. 3 kaweat) -- NIE werdyktotworczy.

Q-D1-PASS: M(1)>0 i Ueff''(1)>0 i omega^2(k)>0 dla k>=0 ORAZ
M(psi)>0, Keff(psi)>0 na (0,4/3). Q-D1-FAIL: przeciwnie.

REJESTR WEJSC [INPUT]: K_geo = gamma = c0 = 1 [LOCK sec. 1]. Formy
CYTAT MD sec. 1 (w [eq:vol-element-M911], V [eq:V-M911],
K [eq:K-coupling-unified], metryka [eq:metric-M911-canonical]).
"""
import numpy as np
import sympy as sp

OUT = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-action-audit-spectrum-insert-2026-09-13/Phase2_output.txt")
lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


psi = sp.Symbol('psi', positive=True)
M = psi ** 6 / (4 - 3 * psi) ** 2      # WYNIK Phase 1 (P1a)
Keff = psi ** 4                        # odczyt B (DZIEDZICZONY)
Ueff = psi ** 4 / 4 - psi ** 3 / 3     # = w*V (tozsamosc P1b)
KA = psi ** 5 / (4 - 3 * psi)          # odczyt A (rownolegle)

em("=" * 78)
em("PHASE 2 -- Q-D1: dyspersja prozni z jednej akcji; K_geo=gamma=c0=1")
em("  [INPUT]; M, Keff, Ueff z Phase 1 (P1b PASS); psi* = 1")
em("=" * 78)

# ---------------------------------------------------- linearyzacja EOM
t, x, k, om, eps = sp.symbols('t x k omega epsilon', real=True)
f = sp.Function('psi')(t, x)
ft, fx = sp.diff(f, t), sp.diff(f, x)
Mf, Kf, Uf = M.subs(psi, f), Keff.subs(psi, f), Ueff.subs(psi, f)
L = sp.Rational(1, 2) * Mf * ft ** 2 - sp.Rational(1, 2) * Kf * fx ** 2 - Uf
EL = (sp.diff(sp.diff(L, ft), t) + sp.diff(sp.diff(L, fx), x)
      - sp.diff(L, f))

mode = 1 + eps * sp.exp(sp.I * (k * x - om * t))
EL_lin = EL.subs(f, mode).doit()
EL_lin = sp.expand(EL_lin)
coef1 = sp.simplify(EL_lin.coeff(eps, 1).subs(eps, 0)
                    / sp.exp(sp.I * (k * x - om * t)))
em()
em("Linearyzacja EOM wokol psi*=1 (psi = 1 + eps e^{i(kx-om t)}),")
em("  wspolczynnik O(eps) EL: %s = 0" % coef1)
sol = sp.solve(sp.Eq(coef1, 0), om ** 2)
em("  => omega^2(k) = %s   [WYPROWADZONE]" % sol)
om2 = sp.simplify(sol[0])
M1 = M.subs(psi, 1)
K1 = Keff.subs(psi, 1)
U2 = sp.diff(Ueff, psi, 2)
U2_1 = U2.subs(psi, 1)
om2_formula = (K1 * k ** 2 + U2_1) / M1
em("  formula LOCKa: [Keff(1)k^2 + Ueff''(1)]/M(1) = %s ;"
   % sp.expand(om2_formula))
em("  simplify(omega^2_wyprowadzone - formula) = %s"
   % sp.simplify(om2 - om2_formula))
m2 = sp.simplify(U2_1 / M1)
cs2 = sp.simplify(K1 / M1)
em()
em("  M(1) = %s ; Keff(1) = %s ; Ueff''(psi) = %s ; Ueff''(1) = %s"
   % (M1, K1, sp.factor(U2), U2_1))
em("  m^2   = Ueff''(1)/M(1) = %s" % m2)
em("  c_s^2 = Keff(1)/M(1)   = %s" % cs2)
em("  omega^2(k) = %s  -- omega^2(0) = %s > 0; brak tachionu"
   % (sp.expand(om2), om2.subs(k, 0)))

# ------------------------------------------------- znaki M, Keff na (0,4/3)
em()
em("Znaki na dziedzinie (0, 4/3):")
dom = sp.Interval.open(0, sp.Rational(4, 3))
zer_M = sp.solveset(sp.Eq(M, 0), psi, dom)
zer_K = sp.solveset(sp.Eq(Keff, 0), psi, dom)
zer_KA = sp.solveset(sp.Eq(KA, 0), psi, dom)
em("  M = psi^6/(4-3psi)^2: licznik psi^6>0, mianownik (4-3psi)^2>0")
em("    na (0,4/3); zera w dziedzinie: %s => M > 0 na (0,4/3)" % zer_M)
em("  Keff = psi^4: zera w dziedzinie: %s => Keff > 0 na (0,4/3)" % zer_K)
em("  odczyt A: K_A = psi^5/(4-3psi): zera: %s; 4-3psi>0 na (0,4/3)"
   % zer_KA)
em("    => K_A > 0 na (0,4/3)")
# kontrola float na siatce dziedziny
grid = np.linspace(1e-6, 4.0 / 3.0 - 1e-6, 200001)
Mg = grid ** 6 / (4.0 - 3.0 * grid) ** 2
Kg = grid ** 4
KAg = grid ** 5 / (4.0 - 3.0 * grid)
em("  kontrola float (siatka 2e5 pkt na (0,4/3)): min M = %.3e > 0: %s;"
   % (Mg.min(), Mg.min() > 0))
em("    min Keff = %.3e > 0: %s; min K_A = %.3e > 0: %s"
   % (Kg.min(), Kg.min() > 0, KAg.min(), KAg.min() > 0))
em("  Ueff''(psi) = psi(3psi-2): Ueff''>0 dla psi>2/3; Ueff''(1)=%s>0"
   % U2_1)

# --------------------------------------------------- odczyt A rownolegle
em()
em("Odczyt A rownolegle (odnotowanie; NIE wplywa na werdykt, PRIMARY=B):")
KA1 = KA.subs(psi, 1)
om2_A = (KA1 * k ** 2 + U2_1) / M1
em("  K_A(1) = %s = Keff(1) -- w psi*=1 odczyty A i B NIE roznicuja"
   % KA1)
em("  omega^2_A(k) = [K_A(1)k^2 + Ueff''(1)]/M(1) = %s (identyczna)"
   % sp.expand(om2_A))
em("  c_s^2(psi) roznicuje poza prozni: B: Keff/M = (4-3psi)^2/psi^2;")
em("    A: K_A/M = (4-3psi)/psi -- obie dodatnie na (0,4/3)")

# ------------------------------- deskryptywnie: wariant znakowanego g^tt
em()
em("-" * 78)
em("DESKRYPTYWNIE (MD sec.3; NIE werdyktotworcze): wariant ZNAKOWANEGO")
em("g^tt (literalna kontrakcja w sygnaturze (-,+,+,+)):")
L_s = (-sp.Rational(1, 2) * Mf * ft ** 2
       + sp.Rational(1, 2) * Kf * fx ** 2 - Uf)
EL_s = (sp.diff(sp.diff(L_s, ft), t) + sp.diff(sp.diff(L_s, fx), x)
        - sp.diff(L_s, f))
coef1_s = sp.simplify(sp.expand(EL_s.subs(f, mode).doit())
                      .coeff(eps, 1).subs(eps, 0)
                      / sp.exp(sp.I * (k * x - om * t)))
sol_s = sp.solve(sp.Eq(coef1_s, 0), om ** 2)
em("  omega^2_sgn(k) = %s -- tachion dla k<1 (omega^2(0) = %s < 0);"
   % (sol_s, sol_s[0].subs(k, 0)))
em("  dokladnie omega^2 = k^2 - 1 audytu (rozdz. 2). LOCK zamraza czlon")
em("  czasowy przez |g^tt| => wariant znakowany NIE jest PRIMARY.")

# ------------------------------------------------------------- werdykt
em()
em("=" * 78)
cond1 = (M1 > 0)
cond2 = (U2_1 > 0)
cond3 = sp.simplify(om2.subs(k, 0)) > 0  # omega^2 rosnace w k^2, min w k=0
cond4 = (zer_M == sp.EmptySet and Mg.min() > 0)
cond5 = (zer_K == sp.EmptySet and Kg.min() > 0)
em("WERDYKT Q-D1 (litera LOCKa sec. 2):")
em("  M(1)=%s>0: %s; Ueff''(1)=%s>0: %s; omega^2(k)=k^2+%s>0 dla k>=0:"
   " %s" % (M1, bool(cond1), U2_1, bool(cond2), om2.subs(k, 0),
            bool(cond3)))
em("  M>0 na (0,4/3): %s; Keff>0 na (0,4/3): %s"
   % (bool(cond4), bool(cond5)))
QD1 = "Q-D1-PASS" if all([cond1, cond2, cond3, cond4, cond5]) \
    else "Q-D1-FAIL"
em("  => %s" % QD1)
em("  m^2 = %s, c_s^2 = %s (odczyt B = PRIMARY);" % (m2, cs2))
em("  odczyt A w psi*=1 identyczny (K_A(1)=1).")
em("=" * 78)

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("zapisano:", OUT)
