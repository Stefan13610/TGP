#!/usr/bin/env python3
"""
Phase 1 — op-sigma-ab-pole-residue. Czy istnieje izolowany biegun sigma_ab z residuum
ustalajacym C_sigma? Anti-Lakatos: wyliczone sympy, nie zakladane.

Bloki:
  P1 spektrum wolnego <sigma sigma>: kontinuum vs biegun (branch cut babla)
  P2 projekcja jadra interakcji na fale L=2 (spin-2): kontakt vs pochodna
  P3 warunek bieguna Bethego-Salpetera: 1 = K_{L=2} * G(M^2)
  P4 OPE coeff 2 m_s^2 vs pozycja bieguna
"""
import sympy as sp

print("="*70); print("P1 — SPEKTRUM wolnego <sigma_ab sigma_cd> (kanal spin-2)"); print("="*70)
p,m = sp.symbols('p m', positive=True)
Pi = sp.atan(p/(2*m))/(4*sp.pi*p)              # babel 3D (op-CG4 Phase 3)
print("Pi(p) =", Pi)
# analityczna struktura: atan(z) ma punkty rozgalezienia w z=+-i, tj. p=+-2 i m
print("atan(p/2m): punkty rozgalezienia przy p = +- 2 i m  => p^2 = -4 m^2 (euklides)")
print(" => CIECIE (kontinuum dwuczastkowe) od progu p^2 = -(2 m_s)^2; BRAK izolowanego bieguna w wolnej teorii.")
print(" (biegun bylby prostym 1/(p^2+M^2); atan to funkcja z cieciem, nie biegun.)")

print("\n"+"="*70); print("P2 — PROJEKCJA jadra interakcji na fale L=2 (spin-2)"); print("="*70)
x = sp.symbols('x', real=True)                  # x = cos(theta_pp')
P2 = sp.legendre(2, x)                           # fala spin-2 w 3D = L=2
print("P_2(x) =", P2)
# (a) kontakt phi^4: jadro = const (niezalezne od kata)
proj_contact = sp.integrate(1*P2, (x,-1,1))
print("kontakt phi^4 (M0): proj_{L=2} = int_{-1}^{1} 1 * P_2 dx =", proj_contact,
      "  => ZERO (kontakt = czysta L=0)")
# (b) interakcja pochodna ~ (k.k') = k k' x  -> P_1; rzut na L=2:
proj_deriv1 = sp.integrate(x*P2, (x,-1,1))
print("pochodna ~ (k.k') ~ x:   proj_{L=2} = int x*P_2 dx =", proj_deriv1, "  => ZERO (L=1, nie L=2)")
# (c) kwadrupolowa pochodna ~ (k.k')^2 ~ x^2 -> ma L=0 i L=2:
proj_deriv2 = sp.integrate(x**2*P2, (x,-1,1))
print("kwadrupol ~ (k.k')^2 ~ x^2: proj_{L=2} = int x^2*P_2 dx =", proj_deriv2,
      "  => NIEZEROWE (tylko interakcja >=2 pochodne ma L=2)")

print("\n"+"="*70); print("P3 — WARUNEK BIEGUNA Bethego-Salpetera w kanale spin-2"); print("="*70)
K2, G = sp.symbols('K2 G')                       # K2 = jadro L=2, G = petla
print("Warunek bound-state: 1 = K_{L=2} * G(M^2)")
print(" - substrat M0 (kontakt): K_{L=2} = 0  => 1 = 0 SPRZECZNE => BRAK bieguna spin-2.")
print(" - interakcja kwadrupolowa: K_{L=2} != 0, ale to czlon >=2-pochodne = IRRELEVANT (UV),")
print("   nie czesc zwalidowanego substratu M0; sila wiazania niegwarantowana.")

print("\n"+"="*70); print("P4 — OPE coeff 2 m_s^2 vs pozycja bieguna"); print("="*70)
print("Heredity EOM coeff M^2 = 2 m_s^2 (closure Path B) = wsp. operatorowy (krotki dystans/OPE),")
print("NIE pozycja bieguna funkcji spektralnej (kontinuum od 4 m_s^2).")
print("Aby 2 m_s^2 bylo biegunem, musialby istniec bound state ~1.41 m_s PONIZEJ progu 2 m_s -")
print("a P3 pokazuje: kontakt M0 NIE wiaze spin-2 => taki biegun NIE powstaje.")

print("\n"+"="*70)
print("PODSUMOWANIE Phase 1 (value-blind, wyliczone):")
print(" C-POLE   = FAIL : wolne <sigma sigma> = kontinuum (ciecie atan), brak izolowanego bieguna.")
print(" C-KERNEL = FAIL : kontakt M0 ma proj_{L=2}=0 (i ~x daje 0); tylko >=2-pochodne (irrelevant) daja L=2.")
print(" C-RESIDUE= FAIL : brak bieguna => brak residuum on-shell => C_sigma NIE ustalone.")
print(" C-MATCH  = USTALENIE : 2 m_s^2 = coeff OPE, NIE masa bieguna.")
print(" => Framework NIE dostarcza pole-residue dla sigma_ab. C_sigma pozostaje WOLNYM parametrem (op-CG4 Ph3).")
print(" => kappa_E = WOLNY PARAMETR sektora radiacyjnego (opcja b uczciwa).")
print("="*70)
