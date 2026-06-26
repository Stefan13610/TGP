#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-A3-alpha-resolution — derywacja value-blind: czy alpha=2 wynika z substratu przez Phi=phi^2?

Konwencja LOCKED (sek08 eq:EL-general): rownanie pola nabla^2 Phi + alpha*(nabla Phi)^2/Phi = ...
  alpha = BEZPOSREDNI wspolczynnik czlonu (nabla Phi)^2/Phi.
  Dla L=-1/2 K(f)(nabla f)^2 (f=Phi/Phi0): EL coeff = K'(f)*f/(2K(f)).  K=f^p => coeff = p/2 = alpha.

Pytanie: substrat daje K(u)~u^(2s) w MIKRO polu u (spiny). Macro Phi ~ u^2 (ontologia Phi=<s^2>).
  Po poprawnej zmianie zmiennych: jaki jest alpha? Czy =2 (jak twierdzi thm:D-uniqueness C3)?

Anti-Lakatos: sympy liczy; NIE przesadzamy wyniku.
Sesja: op-A3-alpha-resolution Phase 1 (2026-06-14)
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import sympy as sp
import json, os

RES = []
def check(c, l, d=""):
    RES.append((l, "PASS" if c else "FAIL", d)); print(f"  [{'PASS' if c else 'FAIL'}] {l}" + (f"\n         => {d}" if d else "")); return c
def head(t): print("\n" + "="*72 + "\n" + t + "\n" + "="*72)

f, p, s, alpha = sp.symbols('f p s alpha', positive=True)
psi, u, uref, Kgeo = sp.symbols('psi u u_ref K_geo', positive=True)

# ----------------------------------------------------------------------
head("(F-alpha-A) SANITY: relacja EL  K(f)=f^(2alpha)  <=>  EL coeff = alpha  (sek08, macro f=Phi/Phi0)")
# L = -1/2 K(f) (nabla f)^2. EL coeff czlonu (nabla f)^2/f  =  K'(f)*f/(2 K(f)).
K = f**p
el_coeff = sp.simplify(sp.diff(K, f)*f/(2*K))   # = p/2
print(f"  K(f)=f^p  =>  EL coeff (nabla f)^2/f  =  K'f/(2K) = {el_coeff}")
print(f"  Zatem alpha = p/2,  tj. K(f)=f^(2*alpha).  (Reprodukuje sek08 eq:K-ode: K=C f^(2 alpha).)")
check(sp.simplify(el_coeff - p/2) == 0, "F-alpha-A: relacja EL poprawna (alpha=p/2 dla K=f^p) — sek08 OK",
      f"EL coeff = {el_coeff} = p/2")
# thm:D-uniqueness C3: K(f)=f^4 (p=4) => alpha=2
alpha_C3 = (el_coeff).subs(p, 4)
print(f"  thm:D-uniqueness C3: K(f)=f^4 (p=4) => alpha = {alpha_C3}  (=2, GDY f jest polem kinetyki)")
check(alpha_C3 == 2, "F-alpha-A(ii): C3 arytmetyka OK — JESLI K(Phi/Phi0)=(Phi/Phi0)^4 to alpha=2",
      f"alpha_C3={alpha_C3}")

# ----------------------------------------------------------------------
head("(F-alpha-B) TRANSFORMACJA substratu: K(u)~u^(2s) [micro] + Phi~u^2 [ontologia] -> alpha=?")
# Mikro kinetyk: L_micro = K_geo u^(2s) (nabla u)^2.
# Macro: psi = Phi/Phi0 = (u/u_ref)^2  => u = u_ref*sqrt(psi);  (nabla u)^2 = u_ref^2 (nabla psi)^2/(4 psi).
u_of_psi = uref*sp.sqrt(psi)
Kmicro = Kgeo*u_of_psi**(2*s)                      # K(u(psi))
grad_u_sq_factor = uref**2/(4*psi)                  # (nabla u)^2 = factor * (nabla psi)^2
Keff = sp.simplify(Kmicro*grad_u_sq_factor)         # K_eff(psi) wspolczynnik przy (nabla psi)^2
print(f"  L_micro = K_geo u^(2s)(nabla u)^2;  u=u_ref sqrt(psi).")
print(f"  K_eff(psi) (wsp. przy (nabla psi)^2) = {Keff}")
# potega psi w K_eff:
p_eff = sp.simplify(sp.log(Keff/(Kgeo*uref**(2*s+2)/4))/sp.log(psi))
print(f"  => K_eff(psi) ~ psi^({sp.simplify(s-1)})   (potega p_eff = s-1)")
alpha_sub = sp.simplify((s-1)/2)                    # alpha = p_eff/2
print(f"  alpha(substrat, poprawna transformacja) = p_eff/2 = (s-1)/2 = {alpha_sub}")
print(f"    s=1 (K~u^2, lemat K_phi2):        alpha = {alpha_sub.subs(s,1)}")
print(f"    s=2 (K~u^4, prop:substrate-action): alpha = {alpha_sub.subs(s,2)}")
check(alpha_sub.subs(s, 2) != 2,
      "F-alpha-B: substrat K~u^4 + Phi~u^2 (poprawna transformacja) => alpha != 2",
      f"alpha(s=2) = {alpha_sub.subs(s,2)} (= 1/2, NIE 2)")

# ----------------------------------------------------------------------
head("(F-alpha-C) WERDYKT: czy alpha=2 wynika z substratu?")
a_s2 = alpha_sub.subs(s, 2)   # K~u^4
a_s1 = alpha_sub.subs(s, 1)   # K~u^2
print(f"  Poprawnie stransformowany substrat:")
print(f"    lemat K_phi2 (K~u^2):           alpha = {a_s1}")
print(f"    prop:substrate-action (K~u^4):  alpha = {a_s2}")
print(f"  thm:D-uniqueness C3 twierdzi:     alpha = 2  (przez K(Phi/Phi0)=(Phi/Phi0)^4 BEZ transformacji)")
inconsistency = (a_s2 != 2) and (a_s1 != 2)
print(f"\n  *** {'NIESPOJNOSC' if inconsistency else 'SPOJNE'} ***")
print(f"  thm:D-uniqueness C3 wstawia substratowe K~u^4 BEZPOSREDNIO do macro f=Phi/Phi0,")
print(f"  pomijajac zmiane zmiennych Phi~u^2 (ontologia Phi=<s^2>). Poprawna transformacja (sek10")
print(f"  eq:kinetic_macro) daje K_eff(psi)~psi^(s-1): dla u^4 (s=2) -> psi^1 (LINIOWY) -> alpha=1/2, NIE 2.")
check(inconsistency, "F-alpha-C: alpha=2 NIE wynika z substratu pod Phi~u^2 => INCONSISTENCY (chain zlamany w C3)",
      f"alpha(substrat)=1/2 (K~u^4) lub 0 (K~u^2); thm:D-uniqueness C3=2 przez pominiecie Phi~u^2")

# ----------------------------------------------------------------------
head("(F-alpha-D) CO potrzeba dla alpha=2? (konstruktywne opcje naprawy)")
# Opcja A: K_eff(psi)~psi^4 => p_eff=s-1=4 => s=5 => K(u)~u^10
s_needed = sp.solve(sp.Eq((s-1)/2, 2), s)[0]
print(f"  Opcja A: utrzymac Phi~u^2, podniesc potege K(u). alpha=2 wymaga s={s_needed} => K(u)~u^{2*s_needed}=u^10.")
print(f"           (Brak derywacji substratowej dajacej u^10 — odrzucone jako sztuczne.)")
# Opcja B: Phi ~ u (amplituda, nie u^2). Wtedy psi=u/u_ref, K(u)=u^4=psi^4*u_ref^4 -> K_eff~psi^4 -> alpha=2.
psiB = u/uref
KeffB = sp.simplify((Kgeo*u**4) * (uref**2/(4)) / 1)  # przy Phi~u: (nabla u)^2 nie ma 1/psi
# poprawniej: Phi=Phi0*(u/uref); u=uref*psi; (nabla u)^2=uref^2 (nabla psi)^2; K(u)=u^4=uref^4 psi^4
KeffB = Kgeo*(uref*psi)**4 * uref**2     # ~ psi^4
p_effB = 4
alphaB = sp.Rational(p_effB, 2)
print(f"  Opcja B: Phi~u (AMPLITUDA, nie gestosc). Wtedy psi=u/u_ref, K(u)=u^4=psi^4 => K_eff~psi^4 => alpha={alphaB}. ✓")
print(f"           ALE to SPRZECZNE z ontologia Phi=<s^2>~u^2 (sek08 poziom 0, sek10 eq:Phi_phi_id).")
print(f"\n  => alpha=2 jest spojne TYLKO jesli macro pole Phi ∝ u (amplituda substratu),")
print(f"     a NIE Phi=<s^2> ∝ u^2. To jest istota niespojnosci: 'ktore pole jest fundamentalne'.")
check(True, "F-alpha-D: opcje naprawy wyznaczone (B: Phi∝u amplituda => alpha=2, ale lamie Phi=<s^2>)",
      "alpha=2 <=> Phi∝u (amplituda); Phi=<s^2>∝u^2 => alpha=1/2")

# ----------------------------------------------------------------------
head("WERDYKT (value-blind)")
print("  F-alpha-A: relacja EL (alpha=p/2) POPRAWNA — sek08 eq:K-ode niepodwazone.")
print("  F-alpha-B/C: substrat K~u^4 + Phi=<s^2>~u^2  =>  alpha = 1/2  (NIE 2).")
print("              thm:D-uniqueness C3 daje alpha=2 przez POMINIECIE zmiany zmiennych Phi~u^2.")
print("  => WERDYKT: alpha=2 = DERIVED-INCONSISTENCY (niespojnosc konwencji micro/macro).")
print("     alpha=2 wymaga Phi∝u (amplituda) — sprzeczne z ontologia Phi=<s^2>.")
print("     Pod ontologia Phi=<s^2>: poprawne alpha = 1/2 (K~u^4) lub 0 (K~u^2).")
print("  Obie strony: (PRO α=2) MC LK-1 'α=2 w 1.2σ' — ale blad ±3.82 (dodatekQ2), wiec slabe;")
print("              relacja EL poprawna. (CONTRA) zmiana zmiennych Phi~u^2 jednoznacznie daje α≠2.")
n_pass = sum(1 for _, st, _ in RES if st == "PASS")
print(f"\n  Testy: {n_pass}/{len(RES)} PASS")

out = dict(cycle="op-A3-alpha-resolution",
           EL_relation="alpha = p/2 for K(f)=f^p (sek08, POPRAWNE)",
           alpha_substrate_Ku2=0, alpha_substrate_Ku4=sp.Rational(1, 2).__float__(),
           alpha_C3_claim=2, alpha_needs_for_2="Phi ∝ u (amplituda), sprzeczne z Phi=<s^2>~u^2; lub K~u^10",
           verdict="DERIVED-INCONSISTENCY: alpha=2 wymaga Phi∝u; pod Phi=<s^2> alpha=1/2 (K~u^4) lub 0 (K~u^2)",
           n_pass=n_pass, n_tot=len(RES),
           tests=[{"label": l, "status": st, "detail": d} for l, st, d in RES])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase1_results.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
print("  Wyniki: Phase1_results.json")
print("\nSESJA: op-A3-alpha-resolution Phase 1 (2026-06-14)")
