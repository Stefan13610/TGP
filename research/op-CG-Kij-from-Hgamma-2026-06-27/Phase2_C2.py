#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 - op-CG-Kij-from-Hgamma : F-CGK-C2 (czteropolowy/wyzszy bond) + agregat F-CGK-D.
=========================================================================================
F-CGK-C2: jaki bond mikroskopowy -J sum h(s_i)h(s_j) generuje field-zalezna sztywnosc
          K(s) ~ s^k ? Wyprowadzenie sztywnosci przez rozwiniecie gradientowe + IBP.
          Czy bond dajacy manuskrypt (K~phi^4, alpha=2) ma wspolczynnik FIKSOWANY czy WOLNY?
Reguly LOCKED w Phase0_LOCK.md.
"""
import sympy as sp

s, a, J, q = sp.symbols('s a J q', positive=True)

print("="*70)
print("Phase 2 - F-CGK-C2 + agregat F-CGK-D (value-blind)")
print("="*70)

# =====================================================================
# F-CGK-C2 : bond separowalny -J h(s_i)h(s_j), h(s)=s^q -> sztywnosc K(s)
# rozwiniecie gradientowe 1D + calkowanie przez czesci (IBP):
#   K(s) coeff przy s'^2  =  J a^2 * h'(s)^2   (po IBP czlonu s'')
# =====================================================================
print("\n[F-CGK-C2] bond -J (s_i s_j)^q -> sztywnosc K(s) (rozwiniecie + IBP)")
h = s**q
hp = sp.diff(h, s)
# sztywnosc (z derywacji w komentarzu naglowka): K ~ J a^2 (h')^2
K_stiff = sp.simplify(J*a**2*hp**2)
e_stiff = sp.simplify(sp.diff(sp.log(K_stiff), s)*s)   # wykladnik K~s^e
print(f"    h(s)=s^q,  h'={hp}")
print(f"    K_stiff(s) = J a^2 (h')^2 = {K_stiff}")
print(f"    wykladnik sztywnosci e_K = {e_stiff}   (= 2(q-1))")

# mapowanie: stopien bondu 2q -> wykladnik amplitudowy e_K -> alpha_ampl=e_K/2
# (KONWENCJA: alpha amplitudowe, K(phi)~phi^{2 alpha_ampl}; phi=amplituda)
# UWAGA: -1/2 z #49 to alpha_EFF w GESTOSCI (Phi=phi^2), alpha_eff=alpha_ampl-1/2.
print("\n    Mapowanie stopnia bondu -> e_K -> alpha (KONWENCJA AMPLITUDOWA; phi=s):")
print("    e_K = 2(q-1) = 2 alpha_ampl  =>  alpha_ampl = q-1 ;  stopien bondu = 2q")
for q_v, tag in [(1, "bilinear  (s_i s_j)^1 -- W eq:B-H!"),
                 (2, "czteropolowy (s_i s_j)^2 -- sek08b H_Gamma"),
                 (3, "szesciopolowy (s_i s_j)^3 -- manuskrypt alpha=2")]:
    deg = 2*q_v
    e_v = sp.simplify(e_stiff.subs(q, q_v))      # wykladnik amplitudowy
    a_ampl = sp.Rational(e_v, 2) if e_v != 0 else sp.Integer(0)
    a_eff = sp.nsimplify(a_ampl - sp.Rational(1,2))   # gestosc
    inB = (q_v == 1)
    print(f"      q={q_v} (stopien {deg}): e_K={e_v}, alpha_ampl={a_ampl}, alpha_eff(gest)={a_eff}  "
          + ("[OBECNY w eq:B-H]" if inB else "[NIEOBECNY; extra aksjomat]"))

# wniosek C2: alpha_ampl=2 -> q=3 -> bond stopnia 6 = (s_i s_j)^3
print("\n    alpha_ampl=2 (manuskrypt) wymaga q=3 => bond szesciopolowy -J (s_i s_j)^3:")
print("      (i)  NIEOBECNY w eq:B-H (tam bilinear, q=1).")
print("      (ii) wspolczynnik J' WOLNY (nowy aksjomat; nie fiksowany symetria Z2).")
print("      (iii) operator wysokowymiarowy => RG-irrelewantny (spojne z C1: wszystkie irrelewantne).")
C2 = "C-AXIOM"
print(f"    F-CGK-C2 => {C2} (nawet augmentacja H_Gamma = nowy wolny aksjomat, nie derywacja)")

# bonus R5: trzy wykladniki amplitudowe = trzy konstrukcje kinetyczne (H_Gamma != F_kin)
print("\n    [bonus R5] trzy alpha_ampl z trzech konstrukcji (nie jednej rodziny!):")
print("       * bilinear s_i s_j (eq:B-H, JEDYNY obecny): kanon w s -> alpha_eff=-1/2 (gestosc, #49)")
print("       * (phi_i phi_j)^2 jako ENERGIA bondu (sek08b H_Gamma): K~phi^2 -> alpha_ampl=1")
print("       * (phi_i phi_j)^2 jako WSPOLCZYNNIK sztywnosci na (phi_i-phi_j)^2 (manuskrypt F_kin): K~phi^4 -> alpha_ampl=2")
print("       => R5 = artefakt H_Gamma != F_kin + rama ampl/gest, NIE czysta rodzina 2/4/6.")

# =====================================================================
# AGREGAT F-CGK-D (z flag Phase 1 + Phase 2)
# =====================================================================
print("\n[F-CGK-D] AGREGAT (value-blind, z flag):")
B_flag = "B-REFUTED"          # Phase 1: eta_lit<<1
C1_flag = "all-irrelevant"    # Phase 1: Delta[O_n]>3 dla n=-1,0,1,2
C2_flag = C2                  # C-AXIOM
print(f"    B  = {B_flag}")
print(f"    C1 = {C1_flag}  (=> C-AXIOM-contribution)")
print(f"    C2 = {C2_flag}")
# regula LOCKED: NON-DERIVABLE jesli (B-REFUTED AND C-AXIOM)
C_axiom = (C1_flag == "all-irrelevant") and (C2_flag == "C-AXIOM")
if (B_flag != "B-REFUTED") or (not C_axiom):
    # check DERIVABLE branches
    verdict = "DERIVABLE or UNDETERMINED (sprawdz reguly)"
else:
    verdict = "NON-DERIVABLE"
print(f"    regula: NON-DERIVABLE <= (B-REFUTED AND C-AXIOM)")
print(f"    => F-CGK-D = {verdict}")

print("\n" + "="*70)
print("WERDYKT op-CG-Kij-from-Hgamma:")
print(f"  F-CGK-D = {verdict}")
print("  alpha=2 (K_ij=J(phi_i phi_j)^2) NIE jest wyprowadzalne z mikro H_Gamma:")
print("   * Gaussian CG bilinearnego bondu -> alpha_eff=-1/2 (#49, anchor potwierdzony);")
print("   * eta_lit<<1 -> brak ucieczki przez wymiar anomalny (B-REFUTED);")
print("   * caly sektor kinetyczny kompozytu Phi=eps RG-irrelewantny (C1; zgodne #39);")
print("   * alpha=2 wymaga szesciopolowego bondu (s_i s_j)^3: nieobecny w eq:B-H,")
print("     wspolczynnik wolny, irrelewantny => NOWY AKSJOMAT, nie derywacja (C2).")
print("  ----")
print("  BONUS: trzy wykladniki (R5) z trzech konstrukcji (H_Gamma != F_kin): bilinear->(-1/2 gest.),")
print("         (phi phi)^2-energia->(alpha=1), (phi phi)^2-sztywnosc->(alpha=2). Tylko bilinear w eq:B-H.")
print("  => alpha=2 RATYFIKOWANE jako nieredukowalny aksjomat (status_map l.72, HONEST_FRAMING).")
print("="*70)
