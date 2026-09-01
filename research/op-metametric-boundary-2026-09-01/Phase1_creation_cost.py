#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metametric-boundary (Phase 1) -- Q1: rownanie stanu kreacji (kontinuum).

LOCK: Phase0_balance.md sec. 3 Phase 1; decyzje FROZEN:
Phase_method_decisions.md sec. 7.
P1a (sympy, formalnie): U(0)=0 < U(1)=1/56, U''(1)=-1, brak minimum
    lokalnego U w (0,1) (U' = g^6(beta-gamma g) > 0 na (0,1)).
P1b (kwadratura radialna): Delta E_create(soliton mu | proznia) i
    (| stan pusty) na R=60, h z tabeli kotwicy {0.05,0.025,0.0125}.
P1c: eps(2pi) = (E_cell - E_cell[g=1])/V dla tla 2pi z npz
    (READ-ONLY, bez relaksowania), 4 tablice (A1.0/A0.7 x N32/N48).
Werdykt Q1 wg litery: granica metametryczna w kontinuum ISTNIEJE <=>
    w policzonym zbiorze jest konfiguracja z mu=0 oddzielajaca rezim
    oplacalnej kreacji (mu<0) od kosztownej (mu>0).

REJESTR WEJSC: g0_mu = phi*0.90548 = 1.4650974 [INPUT], beta=gamma=1
[INPUT], kotwica lam_min=-1.646589 [INPUT, tabela h], tlo 2pi z
../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz [READ-ONLY].
"""
import os
import time
import numpy as np
import sympy as sy
from scipy.integrate import solve_ivp

BETA = 1.0
GAMMA = 1.0
PHI_GOLD = (1 + 5 ** 0.5) / 2
G0_MU = PHI_GOLD * 0.90548           # INPUT (definicja: g0_mu := phi*g0_e)
R_RAD = 60.0
H_LIST = (0.05, 0.025, 0.0125)
NPZ = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz")

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


K_fun = lambda g: g ** 4
U_fun = lambda g: BETA * g ** 7 / 7 - GAMMA * g ** 8 / 8

print("=" * 78)
print("Phase 1 (op-metametric-boundary) -- Q1: rownanie stanu kreacji.")
print("LOCK sec. 3 Phase 1; FROZEN: Phase_method_decisions.md sec. 7.")
print("INPUT: g0_mu=phi*0.90548=%.7f, beta=gamma=1, kotwica lam_min="
      "-1.646589 [INPUT], progi QB-2 {0.197,0.298,0.331} [INPUT,"
      " nieuzywane w P1]." % G0_MU)
print("=" * 78)

# ------------------------------------------------------------------ P1a
stamp("[P1a] sympy: struktura U(g) = beta g^7/7 - gamma g^8/8 (beta=gamma=1)")
g = sy.symbols('g', real=True)
beta_s, gamma_s = sy.Integer(1), sy.Integer(1)
U = beta_s * g ** 7 / 7 - gamma_s * g ** 8 / 8
U0 = sy.simplify(U.subs(g, 0))
U1v = sy.simplify(U.subs(g, 1))
U2 = sy.diff(U, g, 2)
U2_1 = sy.simplify(U2.subs(g, 1))
Up = sy.factor(sy.diff(U, g))
crit = sy.solveset(sy.diff(U, g), g, domain=sy.Interval.open(0, 1))
ok_a1 = (U0 == 0)
ok_a2 = (U1v == sy.Rational(1, 56))
ok_a3 = (U2_1 == -1)
ok_a4 = (crit == sy.EmptySet)
print("  U(0) = %s (wymog 0): %s" % (U0, "PASS" if ok_a1 else "FAIL"))
print("  U(1) = %s (wymog 1/56): %s" % (U1v, "PASS" if ok_a2 else "FAIL"))
print("  U''(1) = %s (wymog -1): %s" % (U2_1, "PASS" if ok_a3 else "FAIL"))
print("  U'(g) = %s; solveset(U'=0, (0,1)) = %s (wymog: zbior pusty"
      " -> brak minimum lokalnego w (0,1)): %s"
      % (Up, crit, "PASS" if ok_a4 else "FAIL"))
sgn = sy.simplify(sy.Piecewise((1, Up > 0), (0, True)))
print("  formalnie: U' = g^6(1-g) > 0 na (0,1) => U scisle ROSNACE na"
      " (0,1);")
print("  wniosek o kierunku relaksacji: U(0)=0 < U(g<1) < U(1)=1/56 --")
print("  proznia g=1 lezy na maksimum lokalnym (U''(1)=-1<0), stan pusty")
print("  g=0 energetycznie nizej; relaksacja monotoniczna g->0 BEZ")
print("  posredniego minimum (zadnej granicy w kontinuum miedzy 1 a 0).")
p1a_pass = ok_a1 and ok_a2 and ok_a3 and ok_a4
print("  P1a: %s" % ("PASS" if p1a_pass else "FAIL"))

# ------------------------------------------------------------------ P1b
stamp("[P1b] profil radialny mu kanoniczny (g0=%.7f), R=60" % G0_MU)


def profile_canonical(g0, R, rtol, atol):
    def rhs(r_, y):
        gv, gp = y
        gv = max(gv, 1e-3)
        drv = gv ** 2 * (BETA - GAMMA * gv) - (2.0 / gv) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / 3.0]
        return [gp, drv - 2 * gp / r_]
    sol = solve_ivp(rhs, [1e-6, R + 0.05], [g0, 0.0], method='DOP853',
                    max_step=0.02, rtol=rtol, atol=atol, dense_output=True)
    assert float(sol.t[-1]) >= R, "profil niekompletny -- STOP"
    return sol


sol12 = profile_canonical(G0_MU, R_RAD, rtol=1e-12, atol=1e-14)
tab = {}
for h in H_LIST:
    N = int(round(R_RAD / h))
    r = (np.arange(N) + 0.5) * h
    vals = sol12.sol(r)
    gv, gp = vals[0], vals[1]
    dens = 0.5 * K_fun(gv) * gp ** 2 + U_fun(gv)
    E_vs_vac = float(np.sum(4 * np.pi * r ** 2 * (dens - U_fun(1.0))) * h)
    E_vs_empty = float(np.sum(4 * np.pi * r ** 2 * (dens - 0.0)) * h)
    tab[h] = (E_vs_vac, E_vs_empty)
    print("  h=%.4f: DeltaE_create(sol|proznia) = %+.6f;"
          " DeltaE_create(sol|pusty) = %+.6f"
          % (h, E_vs_vac, E_vs_empty))
E_vac = tab[0.0125][0]
E_emp = tab[0.0125][1]
conv_vac = abs(tab[0.025][0] - tab[0.0125][0])
conv_emp = abs(tab[0.025][1] - tab[0.0125][1])
print("  wartosci glowne (h=0.0125): DeltaE(sol|vac) = %+.6f"
      " (|delta h/2| = %.2e); DeltaE(sol|pusty) = %+.6f (|delta| = %.2e)"
      % (E_vac, conv_vac, E_emp, conv_emp))
print("  uwaga skalowa: DeltaE(sol|pusty) zawiera objetosciowy koszt"
      " prozni U(1)*V(R) ~ (4pi/3)R^3/56 = %.1f -- rosnie z R^3"
      % (4 * np.pi / 3 * R_RAD ** 3 / 56))
print("  znaki: mu(sol|vac) = %s, mu(sol|pusty) = %s"
      % ("+" if E_vac > 0 else "-", "+" if E_emp > 0 else "-"))

# ------------------------------------------------------------------ P1c
stamp("[P1c] eps(2pi) z npz (READ-ONLY, bez relaksowania)")
mtime_before = os.path.getmtime(NPZ)
data = np.load(NPZ)
d2pi = 2 * np.pi
eps_tab = {}
for tag in ("A1.0", "A0.7"):
    for N in (32, 48):
        g3 = data["2pi__%s__N%d" % (tag, N)]
        h = d2pi / N
        E = float(np.sum(U_fun(g3)) * h ** 3)
        for ax in range(3):
            gn = np.roll(g3, -1, axis=ax)
            dg = (gn - g3) / h
            E += 0.5 * float(np.sum(K_fun(0.5 * (g3 + gn)) * dg ** 2)
                             * h ** 3)
        E_vac_cell = U_fun(1.0) * d2pi ** 3
        eps = (E - E_vac_cell) / d2pi ** 3
        eps_tab[(tag, N)] = eps
        print("  tlo 2pi/%s N=%d: E_cell = %.6f; E_cell[g=1] = %.6f;"
              " eps(2pi) = %+.8f (znak %s)"
              % (tag, N, E, E_vac_cell, eps, "+" if eps > 0 else "-"))
data.close()
mtime_after = os.path.getmtime(NPZ)
ok_ro = (mtime_before == mtime_after)
print("  npz READ-ONLY: mtime niezmieniony: %s"
      % ("PASS" if ok_ro else "FAIL"))
eps48 = eps_tab[("A1.0", 48)]
dconv = abs(eps_tab[("A1.0", 32)] - eps48)
print("  glowna (A1.0, N=48): eps(2pi) = %+.8f (|delta 32->48| = %.2e)"
      % (eps48, dconv))

# ------------------------------------------------------------ werdykt Q1
print("\n" + "=" * 78)
print("WERDYKT Q1 (litera LOCKa sec. 3 Phase 1):")
mus = {"DeltaE(sol|vac)": E_vac, "DeltaE(sol|pusty)": E_emp,
       "eps(2pi) A1.0 N48": eps_tab[("A1.0", 48)],
       "eps(2pi) A1.0 N32": eps_tab[("A1.0", 32)],
       "eps(2pi) A0.7 N48": eps_tab[("A0.7", 48)],
       "eps(2pi) A0.7 N32": eps_tab[("A0.7", 32)]}
for k, v in mus.items():
    print("  %-22s = %+.8f" % (k, v))
signs = set(np.sign(v) for v in mus.values())
zero_cross = (len(signs) > 1) or (0.0 in signs)
print("\n  czy w policzonym zbiorze istnieje konfiguracja z mu=0"
      " ODDZIELAJACA rezim mu<0 od mu>0: %s"
      % ("TAK" if zero_cross else "NIE"))
if zero_cross:
    print("  Q1-POS: granica metametryczna w kontinuum ISTNIEJE w"
          " policzonym zbiorze -> charakteryzacja + Q2 nadal (drzewo"
          " LOCK sec. 6)")
else:
    print("  Q1-NEG: wszystkie policzone koszty kreacji maja ten sam znak"
          " (+) --")
    print("  zadna konfiguracja z mu=0 nie oddziela rezimow; formalizacja")
    print("  P1a: U scisle rosnace na (0,1), wiec relaksacja ku g->0 NIE")
    print("  napotyka granicy w kontinuum (spodziewany negatyw LOCKa,")
    print("  pelnoprawny wynik) -> przejscie do Q2 (podloga substratowa).")
print("  P1a: %s | P1b: policzone (3 siatki) | P1c: policzone (4 tablice,"
      " npz nietkniety: %s)" % ("PASS" if p1a_pass else "FAIL",
                                "TAK" if ok_ro else "NIE"))
print("=" * 78)
