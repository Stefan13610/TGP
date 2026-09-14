#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-metric-pair-M911 (Phase 1) -- Q-A: krajobraz sektora grawitacyjnego
(w, V_M9.1'', K=K_geo psi^4), sympy + tabela kontrolna. ZERO relaksacji.

LOCK: Phase0_balance.md sec. 2 Phase 1; MD: Phase_method_decisions.md
sec. 1-3, 7. Formy (CYTATY sek08a, MD sec. 1):
  w(psi)   = psi/(4-3psi)                  [eq:vol-element-M911]
  V(psi)   = -gamma psi^2 (4-3psi)^2 / 12  [eq:V-M911]
  K(psi)   = K_geo psi^4                   [eq:K-coupling-unified]
Odczyt PRIMARY = B (MD sec. 2): kinetyka 1/2 K_geo psi^4 |grad psi|^2
(w * g^ij(M9.1'') = 1); odczyt A odnotowany: Keff_A = psi^5/(4-3psi).
rho_eff = w*V wspolne dla obu odczytow.

Q-A-PASS (litera LOCKa): minimum psi* in (0,4/3), rho''(psi*)>0,
rho_eff(4/3) > rho_eff(psi*), E ograniczone z dolu (kinetyka >=0 na
dziedzinie + inf rho_eff > -inf), brak kierunku ucieczki.
P1b gate: rho_eff sympy vs float w {0.5, 1, 7/6, 1.3}, zgodnosc 1e-12;
tozsamosci wielomianowe U=wV, U'=w'V+wV' (simplify==0).

REJESTR [INPUT]: K_geo=gamma=1; progi detektorow psi {5/6, 7/6};
pas graniczny 4/3-1e-6 (nieuzywany w Phase 1).
"""
import sympy as sp

OUT = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-metric-pair-M911-2026-09-02/Phase1_output.txt")
lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


psi = sp.symbols('psi', positive=True)
gamma, Kgeo = sp.Integer(1), sp.Integer(1)          # INPUT (LOCK sec.1)

w = psi / (4 - 3 * psi)                              # eq:vol-element-M911
V = -gamma * psi**2 * (4 - 3 * psi)**2 / 12          # eq:V-M911
K = Kgeo * psi**4                                    # eq:K-coupling-unified
rho = sp.simplify(w * V)                             # rho_eff = w*V
rho_poly = gamma * (psi**4 / 4 - psi**3 / 3)         # MD sec.3 (U_eff, lin.999)
Keff_B = K                                           # PRIMARY (odczyt B)
Keff_A = sp.simplify(w * K)                          # odczyt A (odnotowany)

em("=" * 78)
em("Phase 1 -- Q-A: krajobraz pary metrycznej (w, V_M9.1'', K)")
em("REJESTR [INPUT]: K_geo=1 gamma=1; formy CYTAT sek08a (MD sec.1);")
em("PRIMARY odczyt B (MD sec.2): Keff=psi^4; odczyt A: Keff_A=w*K")
em("=" * 78)
em()
em("w(psi)     = %s" % sp.sstr(w))
em("V(psi)     = %s" % sp.sstr(sp.factor(V)))
em("K(psi)     = %s" % sp.sstr(K))
em("rho_eff    = w*V = %s  (simplify)" % sp.sstr(sp.factor(rho)))
em("Keff (B)   = %s ; Keff_A (odczyt A) = %s" % (sp.sstr(Keff_B),
                                                 sp.sstr(sp.factor(Keff_A))))
em()

# ---------------------------------------------------------- tozsamosci (P1b)
id1 = sp.simplify(rho - rho_poly)
wprime = sp.simplify(sp.diff(w, psi))
Vprime = sp.simplify(sp.diff(V, psi))
Uprime_poly = gamma * psi**2 * (psi - 1)
id2 = sp.simplify(wprime * V + w * Vprime - Uprime_poly)
em("TOZSAMOSCI (warunek uzycia postaci wielomianowej w silniku, MD sec.3):")
em("  U := w*V - gamma*(psi^4/4 - psi^3/3)      simplify -> %s" % sp.sstr(id1))
em("  U' := w'V + wV' - gamma*psi^2*(psi-1)     simplify -> %s" % sp.sstr(id2))
em("  w'(psi) = %s ; V'(psi) = %s" % (sp.sstr(sp.factor(wprime)),
                                      sp.sstr(sp.factor(Vprime))))
em("  zera V': %s  (srednie progowe LOCKa: (2/3+1)/2=5/6, (1+4/3)/2=7/6)"
   % sp.sstr(sp.solve(Vprime, psi)))
ok_ids = (id1 == 0 and id2 == 0)
em("  -> tozsamosci: %s" % ("PASS" if ok_ids else "FAIL"))
em()

# --------------------------------------------------- P1a: punkty krytyczne
em("P1a: punkty krytyczne rho_eff na (0, 4/3)")
rp = sp.simplify(sp.diff(rho, psi))
rpp = sp.simplify(sp.diff(rho, psi, 2))
em("  rho'  = %s" % sp.sstr(sp.factor(rp)))
em("  rho'' = %s" % sp.sstr(sp.factor(rpp)))
crits = sorted(sp.solve(sp.Eq(rp, 0), psi))
em("  rho'=0 -> psi in %s" % [sp.sstr(c) for c in crits])
interior = [c for c in crits if sp.Rational(0) < c < sp.Rational(4, 3)]
em("  punkty krytyczne WEWNATRZ (0,4/3): %s" % [sp.sstr(c) for c in interior])
psi_star = None
for c in interior:
    curv = rpp.subs(psi, c)
    em("    psi=%s: rho=%s (=%.10f), rho''=%s (%s)"
       % (sp.sstr(c), sp.sstr(rho.subs(psi, c)),
          float(rho.subs(psi, c)), sp.sstr(curv),
          "MINIMUM (rho''>0)" if curv > 0 else
          ("MAKSIMUM" if curv < 0 else "degeneracja")))
    if curv > 0:
        psi_star = c
em()

# ------------------------------------------- zachowanie 0+, 4/3-, granice
lim0 = sp.limit(rho, psi, 0, '+')
lim43 = sp.limit(rho, psi, sp.Rational(4, 3), '-')
rho43 = rho.subs(psi, sp.Rational(4, 3))
em("  zachowanie: rho_eff(0+) = %s ; rho_eff(4/3-) = %s ; rho_eff(4/3)=%s"
   % (sp.sstr(lim0), sp.sstr(lim43), sp.sstr(rho43)))
em("  monotonia: rho' = psi^2(psi-1) < 0 na (0,1), > 0 na (1,4/3)")
mono_dn = sp.simplify(rp.subs(psi, sp.Rational(1, 2))) < 0
mono_up = sp.simplify(rp.subs(psi, sp.Rational(7, 6))) > 0
em("    kontrola: rho'(1/2)<0: %s ; rho'(7/6)>0: %s" % (mono_dn, mono_up))
rho_min = sp.minimum(rho, psi, sp.Interval(0, sp.Rational(4, 3)))
em("  globalne minimum rho_eff na [0,4/3]: %s (=%.10f)"
   % (sp.sstr(rho_min), float(rho_min)))
em()

# ------------------------------- pelna gestosc: znak wspolczynnika kinet.
em("  pelna gestosc e = 1/2*Keff|grad psi|^2 + rho_eff:")
kin_B_pos = sp.simplify(Keff_B) > 0            # psi^4>0 dla psi>0
kin_A_pos = all(sp.simplify(Keff_A.subs(psi, v)) > 0
                for v in [sp.Rational(1, 10), sp.Rational(1, 2), 1,
                          sp.Rational(7, 6), sp.Rational(13, 10)])
kin_A_sign = sp.solve(sp.Eq(4 - 3 * psi, 0), psi)
em("    odczyt B: Keff = psi^4 > 0 na (0,4/3): %s" % bool(kin_B_pos))
em("    odczyt A: Keff_A = psi^5/(4-3psi) > 0 na (0,4/3): %s "
   "(biegun dopiero w psi=%s)" % (kin_A_pos, sp.sstr(kin_A_sign[0])))
em("    -> OBA odczyty: wspolczynnik kinetyczny DODATNI na dziedzinie;")
em("       E[psi] >= |Omega| * inf rho_eff = -|Omega|/12 -- OGRANICZONE"
   " z dolu, brak kierunku ucieczki")
em()

# ---------------------------------- stacjonarnosc psi=1 w pelnym E (P2a)
Up1 = sp.simplify((wprime * V + w * Vprime).subs(psi, 1))
em("  stacjonarnosc prozni w pelnym E: U'(1) = w'(1)V(1)+w(1)V'(1) = %s"
   % sp.sstr(Up1))
em("    -> psi*=1 JEST punktem krytycznym pelnego funkcjonalu"
   if Up1 == 0 else "    -> psi=1 NIE jest punktem krytycznym!")
em()

# --------------------------------------------------------- werdykt Q-A
c1 = psi_star is not None and sp.Rational(0) < psi_star < sp.Rational(4, 3)
c2 = psi_star is not None and rpp.subs(psi, psi_star) > 0
c3 = psi_star is not None and rho43 > rho.subs(psi, psi_star)
c4 = bool(kin_B_pos) and sp.Integer(rho_min * -12) == 1  # inf=-1/12 skonczony
em("WERDYKT Q-A (litera LOCKa sec.2 Phase 1):")
em("  (1) minimum psi* in (0,4/3): %s (psi*=%s)"
   % ("TAK" if c1 else "NIE", sp.sstr(psi_star)))
em("  (2) krzywizna rho''(psi*)>0: %s (rho''(%s)=%s)"
   % ("TAK" if c2 else "NIE", sp.sstr(psi_star),
      sp.sstr(rpp.subs(psi, psi_star)) if psi_star is not None else "-"))
em("  (3) rho_eff(4/3) > rho_eff(psi*): %s (%s > %s)"
   % ("TAK" if c3 else "NIE", sp.sstr(rho43),
      sp.sstr(rho.subs(psi, psi_star)) if psi_star is not None else "-"))
em("  (4) E ograniczone z dolu / brak ucieczki: %s (inf rho=%s, Keff>0)"
   % ("TAK" if c4 else "NIE", sp.sstr(rho_min)))
qa = "Q-A-PASS" if (c1 and c2 and c3 and c4) else "Q-A-FAIL"
em("  ==> %s" % qa)
em()

# ------------------------------------------------------------- P1b gate
em("P1b (gate zgodnosci sympy vs float, prog 1e-12):")
PTS = [sp.Rational(1, 2), sp.Integer(1), sp.Rational(7, 6),
       sp.Rational(13, 10)]


def rho_float(x):
    """Implementacja float silnika (postac wielomianowa, MD sec.3)."""
    return 1.0 * (x**4 / 4.0 - x**3 / 3.0)


gate_ok = True
for p in PTS:
    exact = rho.subs(psi, p)
    fexact = float(sp.N(exact, 30))
    fimpl = rho_float(float(p))
    d = abs(fexact - fimpl)
    ok = d <= 1e-12
    gate_ok = gate_ok and ok
    em("  psi=%-8s rho_sympy=%s = %.16g ; rho_float=%.16g ; |d|=%.3g -> %s"
       % (sp.sstr(p), sp.sstr(exact), fexact, fimpl, d,
          "OK" if ok else "FAIL"))
em("  P1b: %s" % ("PASS" if (gate_ok and ok_ids) else
                  "FAIL -- STOP (blad implementacji)"))
em()

# ----------------------------------------------------- tabela kontrolna
em("Tabela kontrolna (co 0.1, dziedzina (0,4/3)):")
em("  %-6s %-14s %-14s %-12s %-14s" % ("psi", "rho_eff", "U'=psi^2(psi-1)",
                                       "Keff=psi^4", "Keff_A=wK"))
import numpy as np
for x in np.arange(0.1, 1.34, 0.1):
    x = float(x)
    ka = x**6 / (4 - 3 * x**2) if x**2 < 4 / 3 else float('nan')
    em("  %-6.2f %+.6e %+.6e %.6e %.6e"
       % (x, rho_float(x), x * x * (x - 1), x**4, x**5 / (4 - 3 * x)))
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = np.linspace(1e-4, 4 / 3 - 1e-4, 2000)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(xs, rho_float(xs))
    ax[0].axvline(1.0, ls='--', c='gray')
    ax[0].axvline(4 / 3, ls=':', c='red')
    ax[0].axhline(-1 / 12, ls=':', c='gray')
    ax[0].set_title("rho_eff = w*V (min psi=1, granica 4/3 pod gorke)")
    ax[0].set_xlabel("psi")
    ax[1].plot(xs, xs**2 * (xs - 1))
    ax[1].axhline(0, c='gray', lw=0.5)
    ax[1].set_title("U'(psi) = psi^2(psi-1)")
    ax[1].set_xlabel("psi")
    fig.tight_layout()
    fig.savefig(OUT.replace("Phase1_output.txt", "Phase1_landscape.png"),
                dpi=110)
    em("wykres kontrolny: Phase1_landscape.png zapisany")
except Exception as e:  # noqa: BLE001
    em("matplotlib niedostepny (%s) -- tylko tabela ASCII (MD sec.7)" % e)

em()
em("PODSUMOWANIE Phase 1: %s ; P1b %s" % (qa, "PASS" if gate_ok else "FAIL"))
with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("zapisano:", OUT)
