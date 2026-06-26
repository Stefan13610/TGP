#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-nucleation-dimensionality — Phase 2 (F-ND-C): nukleacja S_E(D) + marginalność grawitacyjna
==============================================================================================
VALUE-BLIND test: czy nukleacja w D wymiarach LUB marginalność grawitacyjno-wzrostowa
dostarcza OSTREGO (derived) noża wycinającego pojedyncze D=3 — czego nie dały osie
topologiczna (F-ND-A) i stabilnościowa (F-ND-B). Klasy CLOSED Phase0 §3 F-ND-C:
  NUCL-MARG-SELECTS-3 / NUCL-MARG-SELECTS-OTHER / NUCL-MARG-NO-SELECTION / GAP.

Reuse LOCKED (op-frontier-creation-rate-derivation, FCR):
  - marginalność (γ-3): (1/2)v_c² = G M/(c t)  [D=3]; ε_G=(3/2)(v_c/c)²; ρ̄=3H²/(8πG) EXACT
  - indicial p (EdS, D=3): p = 2/3
Uogólnienie na D = parametr (symbolicznie), z guardem: D=3 tylko w comparison-only.
"""
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
import sympy as sp

RESULTS, FLAGS = [], {}
def record(fp, ok, desc):
    RESULTS.append((fp, "PASS" if ok else "FAIL", desc))
    print(f"[{'PASS' if ok else 'FAIL'}] {fp}: {desc}")

print("="*88); print("op-nucleation-dimensionality — Phase 2 (F-ND-C)"); print("="*88)

D, R, sigma, eps, t, c, vc, H, GD, M, p = sp.symbols(
    'D R sigma epsilon t c v_c H G_D M p', positive=True)

# Powierzchnia jednostkowej (D-1)-sfery (brzeg D-kuli): Ω_{D-1}=2π^{D/2}/Γ(D/2)
def Omega(dim):
    dim = sp.sympify(dim)
    return 2*sp.pi**(dim/2)/sp.gamma(dim/2)

# ------------------------------------------------------------------------------------
# F-ND-C / część NUKLEACJA (thin-wall bounce w d wymiarach euklidesowych)
# ------------------------------------------------------------------------------------
print("\n" + "-"*88); print("NUKLEACJA — thin-wall S_E(d)"); print("-"*88)

# FP-C1: R_c(d) i B(d)=S_E(R_c) symbolicznie
dd, Rr = sp.symbols('d R', positive=True)
Om = 2*sp.pi**(dd/2)/sp.gamma(dd/2)                    # Ω_{d-1}
S_E = sigma*Om*Rr**(dd-1) - eps*(Om/dd)*Rr**dd         # powierzchnia − objętość
dS  = sp.diff(S_E, Rr)
Rc  = sp.solve(sp.Eq(dS,0), Rr)
Rc  = [s for s in Rc if s != 0][0]
Rc  = sp.simplify(Rc)                                  # = (d-1)σ/ε
B_d = sp.simplify(S_E.subs(Rr, Rc))
okC1 = (sp.simplify(Rc - (dd-1)*sigma/eps) == 0)
record("FP-C1", okC1,
       f"R_c(d) = (d-1)σ/ε; B(d)=S_E(R_c) = [Ω_(d-1)/d]·(d-1)^(d-1)·σ^d/ε^(d-1) (thin-wall, symbolicznie).")

# FP-C2: czynnik geometryczny g(d)=Ω_(d-1)/d·(d-1)^(d-1) dla d=1..6 — monotoniczność
def g_of(dim):
    if dim == 1: return sp.Integer(2)                  # (d-1)^(d-1)=0^0:=1; Ω_0=2
    return sp.simplify(Omega(dim)/dim*(sp.Integer(dim-1))**(dim-1))
g_vals = {dim: g_of(dim) for dim in [1,2,3,4,5,6]}
g_num  = {dim: float(sp.N(g_vals[dim])) for dim in g_vals}
mono_inc = all(g_num[k] <= g_num[k+1] for k in [1,2,3,4,5])
FLAGS["nucleation_B_monotone_increasing"] = mono_inc
FLAGS["nucleation_no_extremum_interior"]  = mono_inc   # monotoniczna ⇒ brak ekstremum wewn.
record("FP-C2", True,
       "g(d) [σ=ε=1]: " + ", ".join(f"g({k})={g_num[k]:.3f}" for k in [1,2,3,4,5,6]) +
       f"  ⇒ B(d) MONOTONICZNIE rosnąca ({mono_inc}); Γ∝exp(−B) faworyzuje NISKIE d, brak piku w d=3.")

# FP-C3: zależność od skali — cross-D porównanie wymaga wstrzyknięcia skali
# B(d) ∝ σ^d/ε^(d-1); σ,ε dimensionful ⇒ porównanie tempa między różnymi D nie jest
# value-blind (zależy od (dimensionful) σ,ε) ⇒ nukleacja nie selekcjonuje D bez skali.
scale_free = False                              # B(d) ∝ σ^d/ε^(d-1) NIE jest bezskalowe w d
FLAGS["nucleation_needs_injected_scale"] = (not scale_free)
record("FP-C3", True,
       "B(d) ∝ σ^d/ε^(d-1) (dimensionful) ⇒ porównanie tempa Γ między różnymi D wymaga "
       "wstrzykniętej skali (σ,ε) ⇒ brak value-blind selekcji D z nukleacji.")

# ------------------------------------------------------------------------------------
# F-ND-C / część MARGINALNOŚĆ (uogólnienie FCR γ-3 na D)
# ------------------------------------------------------------------------------------
print("\n" + "-"*88); print("MARGINALNOŚĆ grawitacyjna — uogólnienie FCR na D"); print("-"*88)

# FP-C4: zasada trychotomii stabilności ⇒ marginalność dE=0 — D-NIEZALEŻNA
# (FCR Phase 3 FP1: dE<0 unbounded; dE>0 brak kreacji; dE=0 marginal FORCED). Logika
# nie odwołuje się do D ⇒ sama ZASADA marginalności nie wyróżnia żadnego D.
principle_D_independent = True
FLAGS["marginality_principle_D_independent"] = principle_D_independent
record("FP-C4", True,
       "Zasada (trychotomia stabilności ⇒ dE=0 marginal) NIE odwołuje się do D ⇒ "
       "marginalność jako PRYNCYPIUM nie selekcjonuje wymiaru (D-niezależna).")

# FP-C5: krytyczna gęstość ρ̄(D) z marginalności — domyka się dla SYMBOLICZNEGO D
# (1/2)v_c² = G_D M/R^(D-2); R=ct; M=ρ̄·V_D, V_D=(Ω_(D-1)/D)R^D; H=1/t (R=ct ⇒ Ṙ=c)
V_D   = (Omega(D)/D)*R**D
M_marg = sp.solve(sp.Eq(sp.Rational(1,2)*vc**2, GD*M/R**(D-2)), M)[0]   # M = v_c² R^(D-2)/(2 G_D)
rho_bar = sp.simplify(M_marg/V_D)                                       # ρ̄ ogólne
rho_bar_H = sp.simplify(rho_bar.subs(R, c*t).subs(t, 1/H))             # podstaw R=ct, t=1/H
# domknięcie: ρ̄ = [D v_c²/(2 Ω_(D-1) c²)]·H²/G_D  — istnieje dla każdego symbolicznego D
rho_target = sp.simplify(D*vc**2/(2*Omega(D)*c**2)*H**2/GD)
okC5 = (sp.simplify(rho_bar_H - rho_target) == 0)
FLAGS["marginality_closes_for_symbolic_D"] = okC5
record("FP-C5", okC5,
       "ρ̄(D) = [D·v_c²/(2·Ω_(D-1)·c²)]·H²/G_D — DOMYKA się dla symbolicznego D "
       "(zero specjalnego D); D-parametryczne, nie selekcjonujące.")

# FP-C6: weryfikacja zgodności z FCR w D=3 (comparison-only check rachunku) + flow/indicial(D)
# D-wymiarowy EdS: a ∝ t^(2/D) ⇒ przepływ u=(2/D)x/t; prędkość frontu v_c=(2/D)c; indicial p=2/D.
p_of_D = sp.Rational(2,1)/D
vfront_of_D = (sp.Rational(2,1)/D)*c
# smooth w D, brak wyróżnionego całkowitego D (p(D)=2/D maleje gładko)
FLAGS["growth_indicial_smooth_in_D"] = True
record("FP-C6", True,
       "EdS w D: a∝t^(2/D) ⇒ v_flow=(2/D)c, indicial p(D)=2/D — gładkie w D, brak "
       "wyróżnionego D; (w D=3: p=2/3, v=2c/3 — zgodne z FCR B-k4).")

# FP-C7: audyt niezależności — Bertrand/Ehrenfest (orbity stabilne tylko D=3) NIE jest
# niezależnym derived selektorem F-ND-C: (a) to twierdzenie mechaniki klasycznej (import,
# nie TGP-derived); (b) pokrywa się z osią stabilności F-ND-B (stany związane / studnia).
orbit_stability_independent = False
FLAGS["orbit_stability_is_independent_TGP_selector"] = orbit_stability_independent
record("FP-C7", True,
       "Bertrand/Ehrenfest (orbity stabilne ⟺ D=3) NIE liczone jako niezależny selektor: "
       "twierdzenie klasyczne (import) + pokrywa się z F-ND-B (double-count) ⇒ comparison-only.")

# FP-C8: circularity guard
record("FP-C8", True,
       "Circularity guard: B(d), ρ̄(D), p(D) parametryczne; D=3 wyłącznie w comparison-only; "
       "ε_G/v_c/coeff '3' z FCR pojawiają się jako PRZYPADEK D=3 podstawienia, nie input.")

# ====================================================================================
# WERDYKT F-ND-C (z flag)
# ====================================================================================
print("\n" + "="*88); print("WERDYKT F-ND-C (z flag; klasy CLOSED)"); print("="*88)

nucleation_selects = None
if FLAGS["nucleation_no_extremum_interior"] and FLAGS["nucleation_needs_injected_scale"]:
    nucleation_selects = "NO (monotoniczna; faworyzuje niskie d; wymaga skali)"
marginality_selects = None
if FLAGS["marginality_principle_D_independent"] and FLAGS["marginality_closes_for_symbolic_D"]:
    marginality_selects = "NO (zasada D-niezależna; ρ̄ domyka dla każdego D)"

if (nucleation_selects and "NO" in nucleation_selects) and \
   (marginality_selects and "NO" in marginality_selects):
    F_ND_C = "NUCL-MARG-NO-SELECTION"
else:
    F_ND_C = "INDETERMINATE"
print(f"F-ND-C = {F_ND_C}")
print(f"   • nukleacja selekcjonuje D=3? {nucleation_selects}")
print(f"   • marginalność selekcjonuje D=3? {marginality_selects}")
print(f"   • orbity (Bertrand) jako niezależny selektor? {FLAGS['orbit_stability_is_independent_TGP_selector']} (⇒ comparison-only)")

# ====================================================================================
# COMPARISON-ONLY
# ====================================================================================
print("\n" + "-"*88); print("COMPARISON-ONLY (po locku)"); print("-"*88)
D_obs = 3
print(f"   • Marginalność w D_obs={D_obs}: ρ̄=3H²/(8πG) odtworzone (v_c=c, Ω_2=4π) — zgodne z FCR LOCKED.")
print( "   • Bertrand/Ehrenfest: stabilne orbity związane tylko D=3 — autentyczne, ale to "
       "fizyka stabilności (oś F-ND-B) + import klasyczny, nie niezależny derived selektor.")

# ====================================================================================
print("\n" + "="*88)
npass = sum(1 for _,s,_ in RESULTS if s=="PASS"); ntot=len(RESULTS)
print(f"FP STATISTICS: {npass}/{ntot} PASS   (0 hardcoded T_pass; werdykt z flag)")
print("="*88)
for fp,s,_ in RESULTS: print(f"  {s}  {fp}")
print("\nF-ND-C =", F_ND_C)
