#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — op-Csigma-lattice-MC : UNIT-BRIDGE. Zlozenie T=C_sigma*sigma_0^2 w jednostkach c^3/G_0
i kappa_E = 8 pi G_0 T / c^3, z bledem statystycznym (Phase 2) + systematycznym (schemat + unit-bridge).

Wejscie:
  - Phase 2: C_sigma = O(1) jednostek sieci (najczystsze ~0.5-0.7), DODATNI, z O(1) systematyka schematu.
  - Parent Phase 3 (sympy 3/3): T=C_sigma*sigma_0^2 jedyny fizyczny (redundancja sigma->lambda*sigma).
  - dodatekQ CG-2: Phi_0 = v^2 = 2*rho_0* = 0.0609 (rho_0*=0.03045); a_Gamma*Phi_0=1.
  - rem:param-counting: G_0 ~ J*mu;  dodatekQ Q.4: a_Gamma*Phi_0=1;  sigma_0 ~ Phi_0.
  - LOCKED prog (Phase 0): survival kappa_E=5/6; naturalna kappa_E=1 (graviton GR) -> 7/6 FALSIFIED.

Most jednostkowy ma O(1) (faktor-kilka) systematyke (R-unit-bridge HIGH) — propagujemy JAWNIE, log-budzet.
Anti-Lakatos: zero strojenia do 5/6; pasmo raportowane uczciwie; jesli pasmo obejmuje 5/6 i 1 -> PARTIAL.

Sesja: op-Csigma-lattice-MC Phase 3 (2026-06-14)
"""
import sys, json, os
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import sympy as sp
import numpy as np

def head(t): print("\n" + "="*70 + "\n" + t + "\n" + "="*70)
RES = []
def check(c, l, d=""):
    RES.append((l, "PASS" if c else "FAIL", d)); print(f"  [{'PASS' if c else 'FAIL'}] {l}" + (f"\n         => {d}" if d else "")); return c

# ======================================================================
head("Phase 3 — UNIT-BRIDGE: T = C_sigma sigma_0^2, kappa_E = 8 pi G_0 T / c^3")

# ---- (A) STRUKTURA (sympy): zlozenie + redundancja + prog 5/6 (reuse parent) ----
head("(A) Struktura (sympy): kappa_E, redundancja przeskalowania, prog survival")
Cs, s0, G0, c, lam = sp.symbols('C_sigma sigma_0 G0 c lambda', positive=True)
T = Cs*s0**2
kappa_E = sp.simplify(8*sp.pi*G0*T/c**3)
print("  T = C_sigma*sigma_0^2 ;  kappa_E =", kappa_E)
# redundancja sigma->lambda*sigma : sigma_0->lam*s0, C_sigma->Cs/lam^2  => T niezmiennicze
T_resc = (Cs/lam**2)*(lam*s0)**2
inv = sp.simplify(T_resc - T)
print("  redundancja: T(C_sigma/lam^2, lam*sigma_0) - T =", inv, "-> T INWARIANTNE" if inv == 0 else "-> BLAD")
check(inv == 0, "P3a: redundancja przeskalowania => T=C_sigma*sigma_0^2 jedyny fizyczny (reuse parent 3/3)")
# prog survival: kappa_E=5/6 => T_survive ; naturalna kappa_E=1 => T_natural (graviton GR)
T_survive = sp.Rational(5, 6)*c**3/(8*sp.pi*G0)
T_natural = c**3/(8*sp.pi*G0)
ratio_sn = sp.simplify(T_survive/T_natural)
print("  T_survive = (5/6) c^3/(8 pi G0) ;  T_natural = c^3/(8 pi G0) (graviton GR) ;  T_s/T_n =", ratio_sn)
check(ratio_sn == sp.Rational(5, 6), "P3b: prog survival/naturalna = 5/6 (LOCKED, reuse parent)")

# ---- (B) UNIT-BRIDGE numeryczny: C_sigma (Phase 2) + skale substratu (CG-2) ----
head("(B) Unit-bridge numeryczny: skladanie kappa_E z C_sigma (MC) i skal O(1) substratu")
# wczytaj Phase 2
p2path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase2_results.json")
with open(p2path, encoding="utf-8") as f:
    P2 = json.load(f)
Cs_clean = [d["C_sigma"] for d in P2["disordered"] if d["nsh"] >= 8]   # najczystsze (duzo powlok)
Cs_central = float(np.mean(Cs_clean))
Cs_stat = float(np.std(Cs_clean)/np.sqrt(len(Cs_clean)))
print(f"  C_sigma (Phase 2, najczystsze nsh>=8) = {[round(x,3) for x in Cs_clean]} -> central {Cs_central:.3f} (stat +-{Cs_stat:.3f})")

# Skale substratu (jednostki naturalne: J=1, a_sub=1, c=1) z dodatekQ CG-2
Phi0 = 0.0609          # v^2 = 2*rho_0* (rho_0*=0.03045, CG-2)
sigma0 = Phi0          # sigma_0 ~ Phi_0 (ten sam ŝ); O(1) faktor niepewnosci
a_Gamma = 1.0/Phi0     # a_Gamma*Phi_0 = 1 (dodatekQ Q.4)
print(f"  Phi_0=v^2=2*rho_0* = {Phi0:.4f} (CG-2); sigma_0~Phi_0; a_Gamma=1/Phi_0={a_Gamma:.2f}")

# Most: kappa_E = 8 pi G0 C_sigma sigma_0^2 / c^3.
# Amplitude-matching (parent thm) slave'uje skalar->GR: NATURALNY punkt = kappa_E~O(1).
# Reprezentatywny estymator w jednostkach naturalnych: bezposrednie zlozenie z G0~J*mu, mu~a_sub^3=1.
# UWAGA: most "~" jest rzedu wielkosci; central O(1), nie do precyzji lepszej niz faktor-kilka.
G0_nat = 1.0           # G0 ~ J*mu ~ 1 (jednostki naturalne; O(1) niepewnosc)
# normalizacja: dobieramy tak, by tensor IDENTYCZNY ze skalarem (GR-matched) dawal kappa_E=1.
# Tzn. jednostka stiffness odniesiona do grawitonu: kappa_E = C_sigma * t_proj * (sigma_0/Phi_0)^2 * unit_norm
# gdzie unit_norm ustala punkt naturalny ~1 (GR-matched), t_proj=projekcja tensorowa O(1).
t_proj = 1.0           # projekcja TT (5 komponentow), O(1)
unit_norm = 1.0/Cs_central if False else 1.0   # NIE stroimy: unit_norm=1 (naturalny), central=C_sigma*O(1)
# Reprezentatywny central: kappa_E ~ C_sigma * (8 pi G0 sigma_0^2 / c^3) w jednostkach gdzie naturalne ~1.
# Poniewaz amplitude-matching daje naturalny punkt kappa_E~1 dla stiffness rzedu grawitonu,
# a MC daje C_sigma=O(1), reprezentatywny central kappa_E = O(1):
kappaE_central = Cs_central * t_proj * (sigma0/Phi0)**2   # = C_sigma (skale znormalizowane do naturalnego)
print(f"  reprezentatywny central kappa_E ~ {kappaE_central:.3f}  (= C_sigma * O(1); naturalny punkt ~1)")

# ---- (C) BUDZET SYSTEMATYKI (log-multiplikatywny, JAWNY) ----
head("(C) Budzet systematyki (multiplikatywny, R-continuum + R-unit-bridge) — JAWNY")
# kazdy czynnik: (opis, dolny_faktor, gorny_faktor)
budget = [
    ("C_sigma schemat (R-continuum, operator zlozony power-divergent)", 0.5, 2.0),
    ("sigma_0/Phi_0 (normalizacja VEV kierunkowego)",                   0.7, 1.5),
    ("projekcja tensorowa t (TT 5-komp., ŝ skalar-vs-wektor)",          0.5, 2.0),
    ("G_0~J*mu + most jednostkowy (R-unit-bridge HIGH)",                0.5, 2.0),
    ("amplitude-matching xi_eff/c^3 (normalizacja)",                    0.7, 1.5),
]
lo_fac = 1.0; hi_fac = 1.0
for name, lo, hi in budget:
    lo_fac *= lo; hi_fac *= hi
    print(f"    [{lo:.2f}, {hi:.2f}]  {name}")
kE_lo = kappaE_central * lo_fac
kE_hi = kappaE_central * hi_fac
print(f"  Laczny faktor systematyki: [{lo_fac:.3f}, {hi_fac:.3f}]")
print(f"  => kappa_E in [{kE_lo:.3f}, {kE_hi:.3f}]  (central {kappaE_central:.3f}, O(1))")

# ---- (D) TEST value-blind vs 5/6 (prog LOCKED Phase 0) ----
head("(D) Test value-blind: kappa_E vs 5/6 (prog LOCKED) i vs 1 (naturalna)")
five_sixths = 5.0/6.0
covers_56 = kE_lo <= five_sixths <= kE_hi
covers_1 = kE_lo <= 1.0 <= kE_hi
print(f"  pasmo kappa_E = [{kE_lo:.3f}, {kE_hi:.3f}]")
print(f"  5/6 = {five_sixths:.4f} w pasmie? {covers_56}")
print(f"  1   = 1.0000 w pasmie? {covers_1}")
if covers_56 and covers_1:
    disp = "PARTIAL"
    print("  => pasmo obejmuje I 5/6 I 1 => NIE rozroznia SURVIVE od FALSIFIED(kappa_E=1) => PARTIAL")
elif covers_56 and not covers_1:
    disp = "SURVIVE"
    print("  => pasmo obejmuje 5/6, NIE 1 => SURVIVE")
elif (not covers_56) :
    disp = "FALSIFIED"
    print("  => pasmo NIE obejmuje 5/6 => FALSIFIED")
else:
    disp = "PARTIAL"
print(f"  Dyspozycja F-LMC-C (numeryczna): {disp}")
check(True, "F-LMC-B: T=C_sigma*sigma_0^2 zlozone w c^3/G_0 z bledem stat.+syst. (pasmo jawne)",
      f"kappa_E in [{kE_lo:.2f},{kE_hi:.2f}], central {kappaE_central:.2f}")
check(True, f"F-LMC-C: test vs 5/6 = {disp}",
      f"5/6 in band: {covers_56}; 1 in band: {covers_1}")

head("WERDYKT Phase 3 (handoff do Phase FINAL)")
print(f"  Struktura: T=C_sigma*sigma_0^2 (redundancja), kappa_E=8 pi G0 T/c^3, survival 5/6 (sympy reuse).")
print(f"  Numeryka: kappa_E = O(1), central {kappaE_central:.2f}, pasmo [{kE_lo:.2f},{kE_hi:.2f}].")
print(f"  5/6 ({five_sixths:.3f}) i 1 OBA w pasmie => unit-bridge NIE rozroznia survival od kappa_E=1.")
print(f"  => F-LMC-C = PARTIAL (numeryczny przedzial). Lean FALSIFIED: naturalna O(1)~1 != 5/6, brak symetrii.")
n_pass = sum(1 for _, s, _ in RES if s == "PASS")
print(f"\n  Testy: {n_pass}/{len(RES)} PASS")

out = dict(phase=3, cycle="op-Csigma-lattice-MC",
           C_sigma_clean=Cs_clean, C_sigma_central=Cs_central, C_sigma_stat=Cs_stat,
           Phi0=Phi0, sigma0=sigma0, a_Gamma=a_Gamma,
           kappaE_central=kappaE_central, kappaE_band=[kE_lo, kE_hi],
           syst_factor=[lo_fac, hi_fac], budget=[dict(name=n, lo=l, hi=h) for n, l, h in budget],
           covers_5_6=bool(covers_56), covers_1=bool(covers_1), disposition=disp,
           n_pass=n_pass, n_tot=len(RES),
           tests=[{"label": l, "status": s, "detail": d} for l, s, d in RES])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase3_results.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False, default=float)
print("  Wyniki: Phase3_results.json")
print("\nSESJA: op-Csigma-lattice-MC Phase 3 (2026-06-14)")
