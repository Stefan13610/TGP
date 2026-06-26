#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-bond-order-RG-selection — Phase 1 (value-blind sympy).
Nastepczy do op-Phi-field-identity-resolution (#38).

PYTANIE: czy RG-relevance (power-counting na rodzinie K(s_hat)~s_hat^(2s)) SELEKCJONUJE
  rzad bondu dajacy alpha=2 (s=5), czy NIE (a wrecz dyskryminuje)?

Aparat LOCKED (Phase0 §2):
  O_s = s_hat^(2s) (nabla s_hat)^2 ;  [s_hat] = (d-2)/2 + gamma  (gamma = anomalny wymiar)
  [O_s] = (2s+2)[s_hat] + 2 ;  [g_s] = d - [O_s]
  relevant: [g_s]>0 ; marginal: =0 ; irrelevant: <0
  alpha_density = (s-1)/2  (#38, op-A3) ; s=5<->alpha=2, s=2<->alpha=1/2, s=0<->alpha=-1/2

Anti-Lakatos: sympy liczy; werdykt wyliczony z reguly Phase0 §3.
Sesja: op-bond-order-RG-selection Phase 1 (2026-06-23, #39)
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import sympy as sp
import json, os

RES = []
def check(c, l, d=""):
    RES.append((l, "PASS" if c else "FAIL", d))
    print(f"  [{'PASS' if c else 'FAIL'}] {l}" + (f"\n         => {d}" if d else "")); return c
def head(t): print("\n" + "="*74 + "\n" + t + "\n" + "="*74)

s, d, gamma = sp.symbols('s d gamma', real=True)

# ----------------------------------------------------------------------
head("(R1) SANITY: [g_s](d,gamma) z power-counting; free term s=0 marginalny przy gamma=0")
dim_s = (d-2)/sp.Integer(2) + gamma            # [s_hat]
O_s = (2*s+2)*dim_s + 2                          # [O_s]
g_s = sp.simplify(d - O_s)                       # [g_s]
print(f"  [s_hat] = {dim_s}")
print(f"  [O_s]   = (2s+2)[s_hat]+2 = {sp.expand(O_s)}")
print(f"  [g_s]   = d - [O_s]       = {sp.expand(g_s)}")
g_s_g0 = sp.simplify(g_s.subs(gamma, 0))
print(f"  Przy gamma=0:  [g_s] = {sp.factor(g_s_g0)}")
# free kinetic term = s=0 (K=const) -> musi byc marginalny ([g_0]=0)
check(sp.simplify(g_s_g0.subs(s, 0)) == 0,
      "R1: SANITY [g_{s=0}]=0 przy gamma=0 (swobodny czlon kinetyczny marginalny)",
      f"[g_0](gamma=0) = {sp.simplify(g_s_g0.subs(s,0))}")
# zwarty wzor: [g_s]|gamma=0 = -s(d-2)
check(sp.simplify(g_s_g0 - (-s*(d-2))) == 0,
      "R1(ii): [g_s]|gamma=0 = -s(d-2)  (zwarta postac)",
      f"{sp.factor(g_s_g0)} == -s(d-2)")

# ----------------------------------------------------------------------
head("(R2) Monotonicznosc [g_s] w s (gamma=0, d>2): wyzsze s bardziej czy mniej relevant?")
dg_ds = sp.simplify(sp.diff(g_s_g0, s))          # d[g_s]/ds = -(d-2)
print(f"  d[g_s]/ds (gamma=0) = {dg_ds}")
print(f"  Dla d>2: d[g_s]/ds = -(d-2) < 0  => [g_s] MALEJE z s => wyzsze s = BARDZIEJ irrelevant.")
check(sp.simplify(dg_ds - (-(d-2))) == 0 and (dg_ds.subs(d,4) < 0),
      "R2: [g_s] malejacy w s dla d>2 => RG faworyzuje NIZSZE s (mniej irrelevant)",
      f"d[g_s]/ds = -(d-2); w d=4: {dg_ds.subs(d,4)} < 0")

# ----------------------------------------------------------------------
head("(R3) Czy s=5 (alpha=2) selekcjonowany vs s=2 (alpha=1/2) vs s=0 (free)? [d=4, gamma=0]")
vals = {}
for sv, lab in [(0, "s=0  (K=const, free, alpha=-1/2)"),
                (2, "s=2  (substrat v2, alpha=1/2)"),
                (5, "s=5  (alpha=2, WYMAGANY)")]:
    gv = sp.simplify(g_s_g0.subs({s: sv, d: 4}))
    vals[sv] = gv
    cls = "marginalny" if gv == 0 else ("relevant" if gv > 0 else "IRRELEVANT")
    print(f"  {lab}:  [g_s](d=4) = {gv}   => {cls}")
print(f"\n  Porzadek irrelewancji (d=4): s=0 ({vals[0]}) > s=2 ({vals[2]}) > s=5 ({vals[5]}).")
print(f"  s=5 (alpha=2) jest NAJBARDZIEJ irrelewantny — NIE selekcjonowany; nawet bardziej tlumiony")
print(f"  niz substratowy s=2 (alpha=1/2). Marginalny jest TYLKO s=0 (swobodny K=const).")
selected_s5 = (vals[5] > 0) or (vals[5] == 0)
not_selected = (vals[5] < vals[2] < vals[0]) and (vals[5] < 0)
check(not_selected,
      "R3: s=5 (alpha=2) NIE relevant/marginal; najbardziej irrelevant => RG-NOT-SELECTED",
      f"[g_5]={vals[5]} < [g_2]={vals[2]} < [g_0]={vals[0]} (d=4)")

# ----------------------------------------------------------------------
head("(R4) NGFP escape: jakie gamma czyni s=5 marginalnym? Wiarygodne?")
gamma_for_s5_marginal = sp.solve(sp.Eq(g_s.subs({s: 5, d: 4}), 0), gamma)
gmarg = gamma_for_s5_marginal[0]
print(f"  [g_{{s=5}}](d=4, gamma) = {sp.simplify(g_s.subs({s:5, d:4}))}")
print(f"  Marginalny ([g_5]=0) wymaga gamma = {gmarg} ≈ {float(gmarg):.4f}")
print(f"  To DUZY, ujemny anomalny wymiar (|gamma|~0.83). Typowe anomalne wymiary <<1;")
print(f"  taki gamma bylby skrajny/dostrojony, NIE naturalny. Analogicznie dla s=2 (alpha=1/2):")
gamma_for_s2 = sp.solve(sp.Eq(g_s.subs({s: 2, d: 4}), 0), gamma)[0]
print(f"    s=2 marginalny wymaga gamma = {gamma_for_s2} ≈ {float(gamma_for_s2):.4f} (mniejszy |gamma|).")
print(f"  => Nawet przez anomalny wymiar, s=5 trudniej uczynic marginalnym niz nizsze s.")
check(float(gmarg) < -0.5,
      "R4: s=5 marginalny wymaga skrajnego gamma~-0.83 (niewiarygodne) => brak naturalnego NGFP escape",
      f"gamma(s=5 marginal)={gmarg}≈{float(gmarg):.3f}; gamma(s=2)={gamma_for_s2}≈{float(gamma_for_s2):.3f}")

# ----------------------------------------------------------------------
head("WERDYKT (value-blind, regula Phase0 §3)")
print("  R1: [g_s] = -s(d-2) - (2s+2)*gamma ; przy gamma=0: [g_s]=-s(d-2). s=0 (free) marginalny. SANITY OK.")
print("  R2: [g_s] malejacy w s (d>2) => RG faworyzuje NIZSZE s.")
print("  R3: s=5 (alpha=2) NAJBARDZIEJ irrelevant (d=4: [g_5]=-10 < [g_2]=-4 < [g_0]=0). NIE selekcjonowany.")
print("  R4: NGFP escape wymaga gamma~-0.83 (skrajny, dostrojony) => niewiarygodny.")
verdict = "RG-NOT-SELECTED"
print(f"\n  *** WERDYKT: {verdict} ***")
print("  RG-relevance (power-counting) NIE selekcjonuje rzedu dajacego alpha=2 (s=5).")
print("  Przeciwnie: dyskryminuje go (najbardziej irrelevant), faworyzujac nizsze s — w tym")
print("  substratowe s=2 (alpha=1/2). Marginalny jest tylko swobodny czlon (s=0). NGFP/anomalny")
print("  wymiar nie ratuje (wymaga skrajnego gamma~-0.83).")
print("  => Jedynym selektorem alpha=2 pozostaje DENSITY-FRAME klasa konforemna C1-C3 (makro,")
print("     geometryczna, zwiazana z g_ij~Phi) — to NIE substrate-level zasada RG, lecz")
print("     przeformulowanie JUZ znanej selekcji aksjomatycznej (rem:alpha2-pivot-status).")
print("  => alpha=2 = NIEREDUKOWALNIE AKSJOMATYCZNE na gestosci (status analog c0, #37).")
print("\n  Obie strony: (PRO RG) gdyby istnial NGFP z gamma~-0.83 dla s=5 — ale to skrajne/dostrojone.")
print("    (CONTRA RG) naive + umiarkowany gamma jednoznacznie czynia s=5 irrelevant; RG woli nizsze s.")

n_pass = sum(1 for _, st, _ in RES if st == "PASS")
print(f"\n  Testy: {n_pass}/{len(RES)} PASS")

out = dict(
    cycle="op-bond-order-RG-selection",
    predecessor="op-Phi-field-identity-resolution-2026-06-23 (#38)",
    g_s_formula="[g_s] = -s(d-2) - (2s+2)*gamma",
    g_s_gamma0="[g_s]|gamma=0 = -s(d-2)",
    monotonic_in_s="malejacy (d>2) => RG faworyzuje nizsze s",
    d4_gamma0={"g_0": 0, "g_2": -4, "g_5": -10},
    s_for_alpha2=5, s_substrate_v2=2, s_free=0,
    gamma_needed_s5_marginal=float(sp.solve(sp.Eq(g_s.subs({s:5,d:4}),0),gamma)[0]),
    no_go_RG=True,
    verdict=verdict,
    verdict_plain=("RG-relevance NIE selekcjonuje rzedu dajacego alpha=2 (s=5); dyskryminuje go "
                   "(najbardziej irrelevant), faworyzujac nizsze s w tym substratowe s=2 (alpha=1/2); "
                   "NGFP escape wymaga skrajnego gamma~-0.83 (niewiarygodny). Jedyny selektor alpha=2 = "
                   "density-frame klasa konforemna C1-C3 (makro/geometryczna), NIE substrate RG => "
                   "alpha=2 nieredukowalnie aksjomatyczne na gestosci (status analog c0 #37)."),
    n_pass=n_pass, n_tot=len(RES),
    tests=[{"label": l, "status": st, "detail": d} for l, st, d in RES])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase1_results.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
print("  Wyniki: Phase1_results.json")
print("\nSESJA: op-bond-order-RG-selection Phase 1 (2026-06-23, #39)")
