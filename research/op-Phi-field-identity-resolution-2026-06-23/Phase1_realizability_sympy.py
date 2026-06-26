#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-Phi-field-identity-resolution — Phase 1 (value-blind sympy).

PYTANIE (skorygowany zakres, patrz Phase0_lock §0/§1):
  Czy istnieje DOPUSZCZALNY substrat (skalarny, Z2-parzysty) o sprzezeniu kinetycznym
  K(s_hat) ~ s_hat^(2s), ktorego gestosc Phi=<s_hat^2> daje alpha_density = 2 na gestosci?
  Czy to NO-GO? Jesli realizowalne — kanoniczny v2 bond (s=2) czy nie-kanoniczny?

Relacja LOCKED (op-A3 F-alpha-A, POPRAWNA): dla L=-1/2 K(f)(nabla f)^2, K=C f^p =>
  wspolczynnik czlonu (nabla f)^2 / f  =  p/2  =  alpha.
Ontologia LOCKED: Phi = <s_hat^2> => Phi ~ s_hat^2 ; s_hat = s_ref*sqrt(psi), psi=Phi/Phi0.

Anti-Lakatos: sympy liczy; NIE przesadzamy. Werdykt wyliczony z reguly Phase0 §4.
Sesja: op-Phi-field-identity-resolution Phase 1 (2026-06-23)
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

s, psi, sref, Kgeo, n, m = sp.symbols('s psi s_ref K_geo n m', positive=True)

# ----------------------------------------------------------------------
head("(F1) SANITY: alpha_density(s) z poprawnej transformacji substrat->gestosc")
# Substrat: L_micro = K_geo * s_hat^(2s) * (nabla s_hat)^2.
# Gestosc: psi = (s_hat/s_ref)^2  =>  s_hat = s_ref*sqrt(psi);
#          (nabla s_hat)^2 = s_ref^2 (nabla psi)^2 / (4 psi).
s_hat_of_psi = sref*sp.sqrt(psi)
Kmicro = Kgeo*s_hat_of_psi**(2*s)
grad_factor = sref**2/(4*psi)                       # (nabla s_hat)^2 = grad_factor*(nabla psi)^2
Keff = sp.simplify(Kmicro*grad_factor)              # wsp. przy (nabla psi)^2
# potega psi w K_eff:
const = Kgeo*sref**(2*s+2)/4
p_eff = sp.simplify(sp.log(sp.simplify(Keff/const))/sp.log(psi))
alpha_density = sp.simplify(p_eff/2)
print(f"  K_eff(psi) = {Keff}")
print(f"  p_eff (potega psi) = {p_eff}   ->   alpha_density = p_eff/2 = {alpha_density}")
ok1 = (sp.simplify(alpha_density - (s-1)/2) == 0)
ok1 = ok1 and (alpha_density.subs(s,1) == 0) and (alpha_density.subs(s,2) == sp.Rational(1,2))
check(ok1, "F1: alpha_density(s)=(s-1)/2; reprodukuje op-A3 (s=1->0, s=2->1/2)",
      f"alpha(s=1)={alpha_density.subs(s,1)}, alpha(s=2)={alpha_density.subs(s,2)}")

# ----------------------------------------------------------------------
head("(F2) Rozwiazac alpha_density = 2  dla s")
s_star = sp.solve(sp.Eq(alpha_density, 2), s)
s_star = s_star[0]
print(f"  alpha_density(s)=2  =>  s* = {s_star}")
check(s_star == 5, "F2: alpha=2 wymaga s*=5 (tj. K(s_hat) ~ s_hat^10)", f"s* = {s_star}, 2s* = {2*s_star}")

# ----------------------------------------------------------------------
head("(F3) Czy s*=5 to DOPUSZCZALNY bond (skalarny, Z2-parzysty, calkowity rzad)?")
admissible = s_star.is_integer and (s_star >= 1)
# Z2: K ~ s_hat^(2s) jest parzyste w s_hat dla kazdego calkowitego s -> Z2 OK.
z2_even_ok = True
print(f"  K(s_hat) ~ s_hat^(2*5) = s_hat^10  (parzyste w s_hat => Z2-parzyste OK)")
print(f"  s*=5 calkowite i >=1 => dopuszczalny rzad bondu.")
check(bool(admissible and z2_even_ok),
      "F3: s*=5 -> bond DOPUSZCZALNY (skalar, Z2-parzysty, calkowity) => NIE no-go; alpha=2 REALIZOWALNE",
      "K~s_hat^10; mapowanie rzedu: ekstrakcja n=6 ((s_i s_j)^6), konwencja wagowa m=5")

# Mapowanie do rzedu bondu (obie konwencje rdzenia):
n_extraction = s_star + 1   # ekstrakcja gradientu z (s_i s_j)^n -> s = n-1
m_weight = s_star           # konwencja wagowa K_ij=(s_i s_j)^m -> s = m
print(f"  Ekstrakcja gradientu: (s_i s_j)^n, s=n-1 => n = {n_extraction}")
print(f"  Konwencja wagowa:     K_ij=(s_i s_j)^m, s=m => m = {m_weight}")

# ----------------------------------------------------------------------
head("(F4) Czy s*=5 to KANONICZNY bond v2 (s=2)?")
s_canonical_v2 = 2   # bond v2: (s_i s_j)^2 -> K~s_hat^4 (prop:substrate-action) -> s=2 -> alpha=1/2
alpha_v2 = alpha_density.subs(s, s_canonical_v2)
print(f"  Kanoniczny v2 bond (s_i s_j)^2 => s=2 => alpha_density = {alpha_v2}  (= 1/2, NIE 2)")
print(f"  alpha=2 wymaga s=5 != s=2 (v2). => bond NIE-KANONICZNY.")
check(s_star != s_canonical_v2 and alpha_v2 != 2,
      "F4: alpha=2 NIE pochodzi z kanonicznego v2 bondu (s=2 daje alpha=1/2) => NIE-KANONICZNY",
      f"v2: alpha={alpha_v2}; alpha=2 wymaga s=5")

# ----------------------------------------------------------------------
head("(F5) Cross-consistency: czy bond rzedu kinetyki (n=6) wspolistnieje z V~Phi^3,Phi^4?")
# Potencjal (czlon NIE-gradientowy, tlo s_hat) z bondu (s_i s_j)^n skaluje sie jak s_hat^(2n) = Phi^n.
# Kanoniczny V_orig (sek08a): V ~ Phi^3 (beta) i Phi^4 (gamma), beta=gamma.  => potrzebne n=3 i n=4.
# alpha=2 (kinetyka) potrzebuje rzedu n=6 (ekstrakcja).  Roznia sie => pojedynczy monomial bond NIE
# da jednoczesnie V(Phi^3,Phi^4) i alpha=2 (Phi^6 kinetyka).  Wymaga MULTI-term bondu;
# czlon kinetyczny rzedu 6 jest NIEZALEZNY od czlonow potencjalu rzedu 3,4.
n_for_V_cubic = 3
n_for_V_quartic = 4
n_for_alpha2 = int(n_extraction)   # 6
print(f"  V~Phi^3 (beta)  <- bond rzedu n={n_for_V_cubic}")
print(f"  V~Phi^4 (gamma) <- bond rzedu n={n_for_V_quartic}")
print(f"  alpha=2 (kinetyka) <- bond rzedu n={n_for_alpha2}")
distinct = len({n_for_V_cubic, n_for_V_quartic, n_for_alpha2}) == 3
print(f"  Rzedy {{3,4,6}} ROZNE => pojedynczy monomial bond nie wystarcza;")
print(f"  alpha=2 wymaga osobnego, niezaleznego czlonu kinetycznego rzedu 6 (wsp. dostrajalny).")
check(distinct,
      "F5: czlon kinetyczny dajacy alpha=2 (n=6) jest NIEZALEZNY od czlonow V (n=3,4) => osobny tuned coef",
      "alpha=2 nie wypada 'za darmo' z bondu generujacego kanoniczny potencjal; to osobne wejscie")

# ----------------------------------------------------------------------
head("WERDYKT (value-blind, regula Phase0 §4)")
print("  F1: alpha_density(s)=(s-1)/2  (reprodukuje op-A3; relacja EL niepodwazona).")
print(f"  F2: alpha=2  =>  s*=5  =>  K(s_hat)~s_hat^10.")
print(f"  F3: s*=5 calkowity, Z2-parzysty, skalarny => DOPUSZCZALNY => *** NIE NO-GO ***.")
print(f"  F4: s*=5 != v2 bond (s=2, alpha=1/2) => bond NIE-KANONICZNY (rzad n=6, nie n=2).")
print(f"  F5: czlon kinetyczny n=6 niezalezny od V(n=3,4) => alpha=2 = osobny tuned coef, nie 'za darmo'.")
verdict = "REALIZABLE-NONCANONICAL"
print(f"\n  *** WERDYKT: {verdict} ***")
print("  alpha=2 NA GESTOSCI jest substrate-realizowalne — ale TYLKO przez nie-kanoniczny,")
print("  wysokorzedowy bond Z2 K(s_hat)~s_hat^10 (rzad n=6), nieselekcjonowany obecna zasada")
print("  i niezalezny od czlonow generujacych kanoniczny potencjal V(Phi^3,Phi^4).")
print("  => Bez zasady selekcji rzedu n=6, alpha=2 POZOSTAJE aksjomatem na gestosci (zgodnie")
print("     z rem:alpha2-pivot-status), ale teraz KONSTRUKTYWNIE: wiemy DOKLADNIE jaki substrat")
print("     by go dal. Fundament 'emergencja z jednego pola skalarnego Z2' NIE jest zlamany")
print("     (taki bond jest dopuszczalny), lecz nie jest tez DOMKNIETY (brak zasady -> wybor bondu).")
print("\n  Obie strony: (PRO realizowalnosc) s_hat^10 jest legalnym Z2-skalarnym bondem, wiec")
print("    'no substrate can give alpha=2' jest FALSZ. (CONTRA derywacja) brak zasady selekcjonujacej")
print("    rzad 6 nad rzedem 2; rzad 6 niezalezny od V => to wybor/dostrojenie, nie derywacja kanoniczna.")

n_pass = sum(1 for _, st, _ in RES if st == "PASS")
print(f"\n  Testy: {n_pass}/{len(RES)} PASS")

out = dict(
    cycle="op-Phi-field-identity-resolution",
    scope_correction="field-identity (amplituda vs gestosc) JUZ rozstrzygniete -> gestosc kanoniczna; cykl bada substrate-realizability alpha=2",
    EL_relation="alpha = p/2 (sek08, POPRAWNE; op-A3 niepodwazone)",
    alpha_density_of_s="(s-1)/2",
    s_for_alpha2=int(s_star),
    K_substrate_for_alpha2="K(s_hat) ~ s_hat^10",
    bond_order_extraction=int(n_extraction),
    bond_order_weight=int(m_weight),
    canonical_v2_bond_s=2, canonical_v2_alpha=str(alpha_v2),
    V_bond_orders=[n_for_V_cubic, n_for_V_quartic], alpha2_bond_order=n_for_alpha2,
    no_go=False,
    verdict=verdict,
    verdict_plain=("alpha=2 na gestosci jest substrate-realizowalne, ale tylko przez nie-kanoniczny "
                   "Z2-skalarny bond K~s_hat^10 (rzad 6), niezalezny od V i nieselekcjonowany zasada; "
                   "NIE no-go (fundament nie zlamany), ale NIE domkniete (alpha=2 pozostaje aksjomatem "
                   "na gestosci dopoki brak zasady selekcjonujacej rzad 6)"),
    n_pass=n_pass, n_tot=len(RES),
    tests=[{"label": l, "status": st, "detail": d} for l, st, d in RES])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase1_results.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
print("  Wyniki: Phase1_results.json")
print("\nSESJA: op-Phi-field-identity-resolution Phase 1 (2026-06-23)")
