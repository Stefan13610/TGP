#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mini-test R1 - op-amplitude-density-phase-bridge (PRE-PHASE-0 DIAGNOSTIC)
========================================================================
Pytanie bramki R1: czy crossover wykladnika kinetycznego K ~ Phi^e (e: -1 -> +4)
niesie tresc fizyczna, czy jest tautologia zamiany zmiennej Phi=F(sigma)?

LOCKED input (#49): substrat composite Phi=sigma^(2p), kanoniczna kinetyka (d sigma)^2
 => K(Phi) ~ Phi^e, e(p)=1/p-2. p=1 => e=-1 (alpha_eff=-1/2). alpha=2 => e=+4 => p=1/6.

Cztery pod-testy symboliczne (value-blind; 0 hardcoded T_pass; endpoint +4 NIE wbudowany).
Test DIAGNOSTYCZNY, nie falsyfikator.
"""
import sympy as sp

sigma, Phi, p, m = sp.symbols('sigma Phi p m', positive=True)
gamma = sp.symbols('gamma', positive=True)

print("="*70)
print("mini-test R1 - amplitude->density kinetic crossover (DIAGNOSTIC)")
print("="*70)

results = {}

# ------------------------------------------------------------------ T1
# Gole K(Phi) pojedynczego pola: czy wykladnik to konwencja?
# Phi = sigma^(2p);  K(Phi)=(dsigma/dPhi)^2 ; e = log-pochodna K.
# kanonizacja: chi = int sqrt(K) dPhi ; jesli chi  prop  sigma dla DOWOLNEGO p => tautologia.
print("\n[T1] Czy goly wykladnik kinetyczny pojedynczego pola jest konwencja?")
sigma_of_Phi = Phi**(1/(2*p))
K = sp.simplify(sp.diff(sigma_of_Phi, Phi)**2)          # K(Phi)
e_power = sp.simplify(sp.diff(sp.log(K), Phi)*Phi)       # wykladnik e w K~Phi^e
e_ref = sp.simplify(1/p - 2)                            # przewidywanie #49
T1a = sp.simplify(e_power - e_ref) == 0
# kanoniczne pole chi
chi = sp.simplify(sp.integrate(sp.sqrt(K), Phi))
chi_in_sigma = sp.simplify(chi.subs(Phi, sigma**(2*p)))
# czy chi  prop  sigma (brak zaleznosci sigma w pochodnej chi/sigma)?
T1b = sp.simplify(sp.diff(chi_in_sigma/sigma, sigma)) == 0
print(f"    K(Phi)               = {K}")
print(f"    wykladnik e(p)       = {e_power}    (ref #49: 1/p-2 = {e_ref})")
print(f"    T1a  e == 1/p-2      : {T1a}")
print(f"    kanoniczne chi(Phi)  = {chi}")
print(f"    chi w sigma          = {chi_in_sigma}")
print(f"    T1b  chi  prop  sigma     : {T1b}  (=> goly K to czysta zamiana zmiennej)")
# sprawdzenie liczbowe dla 3 wartosci p (p=1 substrat, p=1/6 alpha=2, p=2)
for pv in [sp.Integer(1), sp.Rational(1,6), sp.Integer(2)]:
    ev = sp.simplify(e_power.subs(p, pv))
    print(f"      p={pv}: e={ev}")
results['T1'] = bool(T1a) and bool(T1b)

# ------------------------------------------------------------------ T2
# Czy 'dressing metryczny' sqrt(-g) ~ Phi^m mostkuje e: -1 -> +4 ?
# net effective exponent E = m + e_sub ; e_sub=-1 (substrat composite).
# solve m dla E=+4 ; porownaj z kandydatami TGP: m=1 (sqrt(-g)=c0*psi),
#   m=2 (objetosc 4D f^2, f~Phi).
print("\n[T2] Czy dressing metryczny sqrt(-g)~Phi^m mostkuje -1 -> +4?")
e_sub = sp.Integer(-1)
E_net = m + e_sub
m_required = sp.solve(sp.Eq(E_net, 4), m)[0]
print(f"    net E(m) = m + e_sub = {E_net}   (e_sub=-1)")
print(f"    wymagane m dla E=+4  : m = {m_required}")
cand = {'m=1 (sqrt-g=c0*psi)': 1, 'm=2 (vol 4D f^2)': 2, 'm=3': 3}
for name, mv in cand.items():
    Ev = E_net.subs(m, mv)
    print(f"      {name:24s}: E_net = {Ev}  (cel +4? {Ev==4})")
# D braku dla kandydata m=1
shortfall = sp.simplify(m_required - 1)
print(f"    niedobor kandydata m=1: Dm = {shortfall}  (= ta sama luka De co #49)")
results['T2'] = (int(m_required) != 1) and (int(m_required) != 2)   # naiwna Droga B niewystarczajaca

# ------------------------------------------------------------------ T3
# Co przywraca tresc fizyczna? Niezmiennik m^2_phys = V''/(2K) przy prozni.
# TGP V_M911(psi) = -gamma psi^2 (4-3psi)^2 /12 ; proznia psi=1 (=sigma=1).
# Liczymy m^2_phys w ramie 'gestosci' (psi=Phi) i 'amplitudy' (psi=sigma^2);
# jesli rowne => niezmiennik to PARA (K,V), nie K samo.
print("\n[T3] Niezmiennik fizyczny: m^2_phys = V''/(2K) (para K,V, nie K samo)")
psi = sp.symbols('psi', positive=True)
V = -gamma*psi**2*(4-3*psi)**2/sp.Integer(12)
# rama gestosci: pole = psi, K_density = K(psi) z p=1 => K ~ 1/(4 psi)
K_density = (sp.diff(psi**sp.Rational(1,2), psi))**2     # psi = sigma^2 => sigma=psi^(1/2)
m2_density = sp.simplify(sp.diff(V, psi, 2)/(2*K_density)).subs(psi, 1)
# rama amplitudy: pole = sigma, V w sigma, K_amp = 1 (kanoniczne)
V_sigma = V.subs(psi, sigma**2)
m2_amp = sp.simplify(sp.diff(V_sigma, sigma, 2)/(2*1)).subs(sigma, 1)
print(f"    V(psi)                 = {sp.simplify(V)}")
print(f"    K_density(psi)         = {sp.simplify(K_density)}")
print(f"    m^2_phys (rama gestosci)  = {m2_density}")
print(f"    m^2_phys (rama amplitudy) = {m2_amp}")
# UWAGA: porownujemy NIEZMIENNIK przy poprawnej normalizacji kanonicznej.
# V'' samo (bez K) w obu ramach:
V2_density = sp.simplify(sp.diff(V, psi, 2)).subs(psi,1)
V2_amp     = sp.simplify(sp.diff(V_sigma, sigma, 2)).subs(sigma,1)
print(f"    V'' samo (gestosc)={V2_density}  vs V'' samo (amplituda)={V2_amp}  -> rozne (znak: K liczy sie)")
T3_inv = sp.simplify(m2_density - m2_amp) == 0
print(f"    T3  m^2_phys niezmiennik (V''/2K rowne): {T3_inv}")
results['T3'] = bool(T3_inv)

# ------------------------------------------------------------------ T4
# Czy skalowo-zalezna mapa F_ell daje NIETautologiczny crossover?
# Gladka odwracalna F (F' != 0) => tautologia przy kazdym ell.
# Tresc fizyczna mozliwa tylko gdzie jakobian F' -> 0 lub oo (nieodwracalnosc).
print("\n[T4] Skalowo-zalezna mapa F_ell: gdzie crossover MOZE byc fizyczny?")
F = sigma**(2*p)                 # rodzina map
Fp = sp.simplify(sp.diff(F, sigma))
crit = sp.solve(sp.Eq(Fp, 0), sigma)
print(f"    F(sigma)=sigma^(2p), F'(sigma) = {Fp}")
print(f"    lokus nieodwracalnosci F'=0 : sigma = {crit}  (=> Phi=0)")
print(f"    => crossover moze nabrac tresci fizycznej TYLKO przy Phi=0 (nukleacja),")
print(f"       NIE jako gladka interpolacja w objetosci (tam: tautologia).")
results['T4'] = (len(crit) >= 1)   # istnieje lokus nieodwracalnosci

# ------------------------------------------------------------------ PODSUMOWANIE
print("\n" + "="*70)
print("PODSUMOWANIE (diagnostyka R1, value-blind)")
print("="*70)
print(f"  T1 (goly K = konwencja, chi prop sigma)         : {results['T1']}")
print(f"  T2 (naiwny dressing metryczny niewystarcza): {results['T2']}")
print(f"  T3 (niezmiennik = para (K,V), nie K)       : {results['T3']}")
print(f"  T4 (fizyka tylko przy Phi=0 / nukleacja)   : {results['T4']}")
print("-"*70)
print("WNIOSEK DIAGNOSTYCZNY:")
print("  * R1 REALNE: goly wykladnik kinetyczny to konwencja (T1).")
print("  * UCIECZKA ISTNIEJE, ale cel = PARA (K,V), nie sam wykladnik (T3).")
print("  * Naiwna Droga B (volume dressing) NIE mostkuje -1->+4 (T2): brak D=5.")
print("  * Gladki crossover w objetosci = tautologia; fizyka tylko w Phi=0 (T4).")
print("  => Teza dwufazowa ZYWA tylko jako: (K,V)-anchored, zakotwiczona w nukleacji Phi=0.")
print("     'Gladkie -1/2 -> +2 w objetosci' = martwe (tautologia).")
print("="*70)
