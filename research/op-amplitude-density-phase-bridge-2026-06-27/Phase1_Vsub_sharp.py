#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 (first-pass) - ostry test R1 z ODNALEZIONYM V_sub.
==========================================================
Wejscia z rdzenia (axioms/substrat/dodatekB_substrat.tex, sek08b_ghost_resolution.tex):
  - mikro Hamiltonian (eq:B-H):
        H_Gamma = sum_i[(m0^2/2) s^2 + (lam0/4) s^4] - J sum_<ij> s_i s_j   (s = s-hat, Z2)
      => na poziomie s: KANONICZNA kinetyka + kwartyk (zwykla teoria phi^4 w s).
  - slownik:  phi=(Phi/Phi0)^(1/2),  Phi=<s^2>  (gestosc = parzysta w s).
  - V_sub (GL w amplitudzie phi):  U(phi)=(beta/3)phi^3-(gamma/4)phi^4, beta=gamma w prozni.
  - sprzezenie geometryczne v2 (eq:K-geometric):  K_ij=J(phi_i phi_j)^2
      => F_kin=int (Kgeo/2) phi^4 (grad phi)^2  => rdzen TWIERDZI alpha=2 'z pierwszych zasad'.
  - #49 (LOCKED): naturalny coarse-graining s->Phi=<s^2> daje alpha_eff=-1/2 (K~Phi^-1).

PYTANIE OSTRE: czy kwartyczne sprzezenie K_ij=J(phi_i phi_j)^2 (zrodlo alpha=2)
JEST obecne w mikro H_Gamma (bilinearny bond -J s_i s_j), czy to NIEZALEZNY postulat
na poziomie parametru porzadku? Rozstrzyga, czy alpha=2 jest first-principles, czy selekcja.
"""
import sympy as sp

# pole i jego pochodne w rozwinieciu gradientowym 1D: s(x+a) = s + a*sx + a^2/2*sxx
s, sx, sxx, a, J = sp.symbols('s sx sxx a J', real=True)
phi, beta, gamma = sp.symbols('phi beta gamma', positive=True)

print("="*70)
print("Phase 1 first-pass - ostry test R1 (V_sub odnaleziony)")
print("="*70)

s_j = s + a*sx + a**2*sxx/2     # sasiad

# ---- T_A: bilinearny bond H_Gamma  -J s_i s_j  -> kinetyka kanoniczna w s ----
print("\n[T_A] Bond mikro H_Gamma: -J s_i s_j  -> jaka kinetyka w s?")
# tozsamosc: -J s_i s_j = (J/2)(s_i-s_j)^2 - (J/2)(s_i^2+s_j^2)
diff2 = sp.expand((s - s_j)**2)
coef_grad_A = sp.simplify(diff2.coeff(a, 2))     # wspolczynnik a^2 -> (grad s)^2
print(f"    (s_i - s_j)^2  |_a^2  = {coef_grad_A}     (=> (grad s)^2)")
K_A = sp.simplify(sp.Rational(1,2)*J*coef_grad_A/ (sx**2))   # prefaktor przy (sx)^2
print(f"    prefaktor kinetyczny K_s = {K_A}")
print(f"    zaleznosc od pola s? d/ds = {sp.diff(K_A, s)}  -> {'STALA (kanoniczna)' if sp.diff(K_A,s)==0 else 'pole-zalezna'}")
e_A = 0 if sp.diff(K_A, s)==0 else None
print(f"    => wykladnik e_s = {e_A} (kanoniczna)  => alpha=1/2 w s; w Phi=s^2: K~Phi^-1 (=#49!)")

# ---- T_B: sprzezenie geometryczne v2  K_ij=J(phi_i phi_j)^2 na (phi_i-phi_j)^2 ----
print("\n[T_B] Sprzezenie geometryczne v2: K_ij=J(phi_i phi_j)^2 * (phi_i-phi_j)^2")
phi_j = phi + a*sx + a**2*sxx/2     # tu sx := grad phi (reuse symbolu)
Kij = J*(phi*phi_j)**2
term_B = sp.expand(Kij*(phi - phi_j)**2)
coef_grad_B = sp.simplify(term_B.coeff(a, 2))
K_B = sp.simplify(coef_grad_B/(sx**2))
print(f"    prefaktor kinetyczny K(phi) = {K_B}")
e_B = sp.simplify(sp.diff(sp.log(K_B), phi)*phi)
print(f"    wykladnik e_phi w K~phi^e = {e_B}   => alpha = e/2 = {sp.simplify(e_B/2)}")

# ---- T_C: czy te dwa sprzezenia to TO SAMO? (struktura na poziomie s) ----
print("\n[T_C] Czy K_ij=J(phi_i phi_j)^2 wystepuje w mikro H_Gamma?")
print("    H_Gamma bond  : -J s_i s_j         (BILINEARNY w polu)")
print("    geom coupling : +J (phi_i phi_j)^2 (KWARTYCZNY w polu; phi~|s|)")
print("    -> phi_i phi_j ~ |s_i||s_j| ; (phi_i phi_j)^2 ~ s_i^2 s_j^2  (kwartyk bond)")
print("    -> s_i^2 s_j^2  !=  s_i s_j  : RÓŻNE operatory bond.")
struct_same = sp.simplify((s**2)*(s_j**2) - (s*s_j))  # != 0 generycznie
print(f"    s_i^2 s_j^2 - s_i s_j (test rownosci) = {sp.simplify(struct_same)}  (!=0 => RÓŻNE)")
print("    WNIOSEK T_C: kwartyczne K_ij (zrodlo alpha=2) NIE jest obecne w eq:B-H.")
print("                 alpha=2 = POSTULAT v2 na poziomie parametru porzadku, nie z H_Gamma.")

# ---- T_D: czy mapa s -> Phi=<s^2> to redefinicja pola? (tautologia czy CG?) ----
print("\n[T_D] Mapa s -> Phi=s^2 : redefinicja (tautologia) czy coarse-graining?")
Phi_of_s = s**2
dPhi = sp.diff(Phi_of_s, s)
print(f"    dPhi/ds = {dPhi}  -> zero przy s=0 (nieodwracalna); s i -s -> ta sama Phi (sign-blind)")
print("    => NIE jest redefinicja pola (gubi znak + degeneruje w s=0).")
print("    => crossover przez s->Phi=<s^2> NIE jest tautologia (R1 omijalne TU),")
print("       ALE to dokladnie CG z #49 -> daje alpha_eff=-1/2, nie +2.")

# ---- V_sub jako para (K,V) do nastepnego kroku ----
print("\n[V_sub] Potencjal substratu (z rdzenia), w amplitudzie phi:")
U = sp.Rational(1,3)*beta*phi**3 - sp.Rational(1,4)*gamma*phi**4
Uvac = sp.simplify(sp.diff(U, phi).subs(phi,1))
print(f"    U(phi) = {U}")
print(f"    U'(1) = {Uvac}  -> =0 wymusza beta=gamma (proznia)")
U_bg = U.subs(beta, gamma)
print(f"    U(phi)|_(beta=gamma) = {sp.factor(U_bg)}")

print("\n" + "="*70)
print("SYNTEZA Phase 1 first-pass:")
print(f"  T_A: mikro bond -J s_i s_j -> KANONICZNA kinetyka (e_s=0) -> w Phi: K~Phi^-1 (#49).")
print(f"  T_B: geom v2 K_ij=J(phi phi)^2 -> K(phi)~phi^{int(e_B)} -> alpha={int(sp.simplify(e_B/2))} (manuskrypt).")
print(f"  T_C: kwartyczny K_ij NIE jest w eq:B-H (s^2 s^2 != s s) -> alpha=2 = POSTULAT v2.")
print(f"  T_D: s->Phi=<s^2> = coarse-graining (nie redef) -> R1 omijalne, ale to droga #49 (-1/2).")
print("  ----")
print("  WNIOSEK: V_sub ZNALEZIONY. alpha=2 NIE wynika z mikro H_Gamma; pochodzi z osobnego")
print("  postulatu geometrycznego v2 (K_ij=J(phi_i phi_j)^2). To RATYFIKUJE #49 od strony")
print("  analitycznej i USTAWIA teze usera: -1/2 = mikro (s, H_Gamma), +2 = postulat v2 (phi).")
print("  Luka -1/2 -> +2 = status sprzezenia geometrycznego v2 (postulat vs derywacja).")
print("="*70)
