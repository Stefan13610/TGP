#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mini-test R1 v2 - KOREKTA artefaktow v1 + ostra postac pytania R1.
==================================================================
v1 (R1_test.py) = nienaruszony zapis value-blind. Tu poprawiam dwa artefakty:
  - T3 v1: 'masa' liczona w psi=1, ktore NIE jest ekstremum V_M911 (V'(1)=gamma/3 != 0,
           cecha TGP: proznia podtrzymywana zrodlem). => nierownosc byla artefaktem.
  - T4 v1: sigma zadeklarowane positive => sympy wykluczyl brzeg sigma=0. Lokus JEST w Phi=0.

OSTRA POSTAC R1 (klucz): pole-redefinicja Phi=F(sigma) zachowuje CALA fizyke. Wiec
'crossover' bedzie REALNY tylko jesli (K_sub, V_sub) i (K_TGP, V_TGP) NIE sa rownowazne
przez redefinicje pola. Test: czy obie ramy daja TO SAMO pole kanoniczne chi?
"""
import sympy as sp

psi, sigma, chi, gamma = sp.symbols('psi sigma chi gamma', positive=True)
sig = sp.symbols('sig', real=True)   # bez zalozenia dodatniosci (dla T4 brzeg)

print("="*70)
print("mini-test R1 v2 - korekta + ostra postac")
print("="*70)

# ---- T3' : czy redefinicja pola zachowuje fizyke? (pole kanoniczne chi) ----
# Rama gestosci: L = K(psi)(d psi)^2, K=1/(4 psi).  g := 2K = 1/(2 psi).
#   chi_density = int sqrt(g) d psi
# Rama amplitudy: L = (d sigma)^2 = (1/2)*2*(d sigma)^2, g=2, psi=sigma^2.
#   chi_amp = int sqrt(2) d sigma, z sigma=sqrt(psi)
print("\n[T3'] Czy obie ramy daja TO SAMO pole kanoniczne chi(psi)?")
K = 1/(4*psi)
g_density = 2*K
chi_density = sp.simplify(sp.integrate(sp.sqrt(g_density), psi))
# rama amplitudy:
chi_amp_sigma = sp.integrate(sp.sqrt(2), sigma)            # = sqrt(2)*sigma
chi_amp = sp.simplify(chi_amp_sigma.subs(sigma, sp.sqrt(psi)))
print(f"    chi_density(psi) = {chi_density}")
print(f"    chi_amp(psi)     = {chi_amp}")
same_canon = sp.simplify(chi_density - chi_amp) == 0
print(f"    T3' chi_density == chi_amp : {same_canon}")
print("      => redefinicja Phi=sigma^2 zachowuje fizyke (to SAMO pole kanoniczne).")
print("      => 'goly crossover' (samo K, ta sama V przepisana) = TAUTOLOGIA. Potwierdzone.")

# ---- T3'' : ostra postac - kiedy crossover JEST realny? ----
# Realny <=> (K_sub,V_sub) i (K_TGP,V_TGP) NIE sa redef-rownowazne.
# Probe: dwie RÓŻNE pary, ktorych nie da sie zmapowac. Pokazujemy mechanizm:
# jesli V_sub w ramie kanonicznej != V_TGP w ramie kanonicznej (po staloej/skali) => rozne teorie.
print("\n[T3''] Ostra postac: dyskryminator redef-rownowaznosci par (K,V)")
print("    Mechanizm (do Phase 1, wymaga V_sub substratu):")
print("      1. policz chi_sub z K_sub=Phi^-1 ; chi_TGP z K_TGP=Phi^+4")
print("      2. wyraz V_sub(chi_sub) i V_TGP(chi_TGP)")
print("      3. rowne (mod stala/skala) => ta sama teoria (crossover fake);")
print("         rozne => RÓŻNE teorie (substrat != TGP; #49 stoi jako 'nie ta sama teoria').")
# demonstracja na chi dla K_TGP=Phi^4:
KTGP = psi**4
chi_TGP = sp.simplify(sp.integrate(sp.sqrt(2*KTGP), psi))
print(f"    chi_TGP(Phi) dla K=Phi^4 : {chi_TGP}   (chi ~ Phi^3 - silnie nieliniowe)")
print(f"    chi_sub(Phi) dla K=Phi^-1: {chi_density}  (chi ~ Phi^(1/2))")
print("    => mapy kanoniczne RÓŻNE funkcyjnie (Phi^3 vs Phi^1/2); zgodnosc par zalezy")
print("       wylacznie od tego, czy potencjaly skladaja sie do tej samej V(chi). [Phase 1]")

# ---- T4' : lokus nieodwracalnosci (poprawnie, z brzegiem) ----
print("\n[T4'] Lokus nieodwracalnosci mapy F (z brzegiem sigma=0):")
p = sp.symbols('p', positive=True)
F = sig**(2*p)
Fp = sp.diff(F, sig)
sols = sp.solve(sp.Eq(Fp, 0), sig)             # sig real, brzeg dozwolony
print(f"    F(sig)=sig^(2p), F'(sig)={Fp}")
print(f"    F'=0 przy sig = {sols}   (brzeg Phi=0)")
print("    => jedyny punkt, gdzie redefinicja TRACI odwracalnosc = Phi=0 (nukleacja).")
print("       Tam i tylko tam 'crossover' moze byc czyms wiecej niz zmiana zmiennej.")

print("\n" + "="*70)
print("WNIOSEK v2 (poprawiony):")
print("  T3' : redef Phi=sigma^2 daje to samo pole kanoniczne -> goly crossover = TAUTOLOGIA.")
print("  T3'': realny crossover <=> pary (K,V) NIE-redef-rownowazne -> wymaga V_sub [PHASE 1].")
print("  T4' : jedyny lokus nieodwracalnosci = Phi=0 (nukleacja).")
print("  RAZEM z v1: R1 realne; teza zywa tylko jako (K,V)-niezgodnosc zakotwiczona w Phi=0.")
print("="*70)
