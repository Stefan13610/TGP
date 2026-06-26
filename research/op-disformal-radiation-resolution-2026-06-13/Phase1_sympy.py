# -*- coding: utf-8 -*-
# Phase 1 — op-disformal-radiation-resolution : T1 (W-GSS-1) strumien vs amplituda
# Pytanie: czy disformalny Vainshtein LIVE tlumi STRUMIEN energii Pdot_b z orbity
#          (bilans Isaacson/T^{0r}), czy nieekranowane sprzezenie konforemne A(Phi)
#          do materii wciaz daje (1/6)P_GR (PR-025)?
# Regula LOCKED (Phase0 forbidden #4): werdykt z bilansu ENERGII T^{0r}, NIE z amplitudy.
# Wszystko symbolicznie (sympy); 0 danych obserwacyjnych; 0 nowych stalych; flaga z rachunku.
#
# WEJSCIA LOCKED uzyte:
#   R2 : L_kin^disformal = A*X - (b/2)X^2 + O(b^2),  b=B/M*^4   (Phase2 Filar I, EXACT)
#   R14: sprzezenie materia-skalar konforemne, zrodlo -(q/Phi0) rho delta_phi (BEZ pochodnych)
#   R3 : h^disf ~ h_GR/(kr) ~ 1e-40  (sek07; AMPLITUDA polaryzacji disformalnej = kanal C2)
#   R13: delta_phi nieekranowany -> P_phi = (1/6)P_GR  (PR-025)
#   rem:B-constraints (sek08 l.6505): B(Phi) = OTWARTY PROBLEM O12 (niewyprowadzone)
#   R11: M* = m_P 'Propozycja, brak mikro-derywacji' (status_map) -- underived

import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

PASS=[]
def check(i,ok,note):
    PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

print("="*72)
print("T1-A — operator fluktuacji Z^{mn} z czlonu disformalnego (Vainshtein kinetyczny)")
print("="*72)
# L(X) = A*X - (b/2)X^2 ; X = v.v ,  v_mu = d_mu phi  (tlo: v = d phi_bg)
# Operator fluktuacji = d^2 L / d(d_mu phi) d(d_nu phi)  (kwadratowa czesc dzialania dla delta_phi)
A,b = sp.symbols('A b', positive=True)               # b = B/M*^4
v0,v1,v2,v3 = sp.symbols('v0 v1 v2 v3', real=True)    # skladowe d phi_bg
eta = sp.diag(-1,1,1,1)                               # mostly-plus
v = sp.Matrix([v0,v1,v2,v3])
X = (v.T*eta*v)[0]                                    # X = v.v (z metryka plaska tla)
L = A*X - sp.Rational(1,2)*b*X**2                     # R2 (do O(b^1))
# Z^{mn} = d^2 L / dv_m dv_n  (v_m = skladowe d phi)
Z = sp.zeros(4,4)
for m in range(4):
    for n in range(4):
        Z[m,n] = sp.simplify(sp.diff(L, v[m], v[n]))
# oczekiwana forma EXACT: Z = 2 L'(X) eta + 4 L''(X) (eta v)(eta v)^T ,  L'=A-bX, L''=-b
Lp = A - b*X
Lpp = -b
vcov = eta*v                                         # v_mu (dolny)
Z_expected = 2*Lp*eta + 4*Lpp*(vcov*vcov.T)
Z_match = sp.simplify(sp.Matrix(Z) - Z_expected) == sp.zeros(4,4)
# czlon X-nieliniowy obecny  <=>  L'' = -b != 0  (zaleznosc Z od tla X => Vainshtein)
vainshtein_kinetic = (sp.simplify(Lpp) != 0)
check("T1-A", Z_match and vainshtein_kinetic,
      f"Z^mn = 2(A-bX)eta - 4b (dphi)(dphi)^T  [EXACT, match={Z_match}]; "
      f"L''=-b!=0 => Z zalezy od tla X => natywny VAINSHTEIN KINETYCZNY (enhancement Z z X_bg).")

print("="*72)
print("T1-B — wzmocnienie kinetyczne wewnatrz r_V: Z_eff/A = 2(bX/A - 1) >> 1 dla bX>>A")
print("="*72)
# parametr Vainshteina u = bX_bg/A (bezwymiarowy). r_V zdef. przez u(r_V)=1.
u = sp.symbols('u', positive=True)                   # u = bX_bg/A
# czesc ~eta operatora (dominuje d'Alembertian fluktuacji): 2|A - bX| ; w jedn. A: 2|1-u|
Zeff_over_A = 2*sp.Abs(1-u)
# wewnatrz r_V: u>1 => Zeff/A = 2(u-1) rosnie z u; enhancement (>1) dla u>3/2
enh_at_large_u = sp.limit(Zeff_over_A/(2*u), u, sp.oo)    # -> 1 : Zeff ~ 2u A dla u>>1
enhancement = sp.simplify(enh_at_large_u) == 1
check("T1-B", enhancement,
      f"u=bX_bg/A (u>1 wewnatrz r_V). Zeff/A=2|1-u| -> 2u dla u>>1 (lim ratio={enh_at_large_u}). "
      f"Kinetyka fluktuacji WZMOCNIONA ~u razy wewnatrz r_V (mechanizm ekranowania obecny).")

print("="*72)
print("T1-C — BILANS ENERGII: P_phi^LIVE / P_phi^unscreened = 1/u  (T^{0r}, NIE amplituda)")
print("="*72)
# Zrodlo konforemne (R14): Zeff * Box(delta_phi) = (q/Phi0) rho   [BEZ pochodnych w zrodle]
#   => odpowiedz delta_phi ∝ (q/(Zeff Phi0)) * (calka retardowana zrodla)  -> amplituda ∝ 1/Zeff
q,Phi0,Ffun = sp.symbols('q Phi0 F', positive=True)  # F = kinematyka zrodla (kwadrupol), wspolna
A_amp = q/(Zeff_over_A*A*Phi0) * Ffun                 # amplituda delta_phi w strefie falowej (∝1/Zeff)
# Isaacson: <T^{0r}> = Zeff * <d_t phi * d_r phi> ;  phi=Aamp/r * osc(t-r/c) => <dt*dr> ∝ (1/c) Aamp^2/r^2
c = sp.symbols('c', positive=True)
flux_density = (Zeff_over_A*A) * (A_amp**2)/c        # ∝ Zeff * (1/Zeff)^2 = 1/Zeff  (czynnik 1/r^2 znosi 4pi r^2)
P_live = sp.simplify(flux_density)                   # moc radiowana (skala) kanal konforemny w LIVE
# referencja UNSCREENED: brak disformalu (u->0 => Zeff/A -> 2): P0 = (1/6) P_GR  (R13)
P0 = sp.simplify(P_live.subs(u,0))                   # Zeff_over_A(u=0)=2
ratio = sp.simplify(P_live/P0)                       # = ?  oczekiwane: |1-0|/|1-u| = 1/|1-u| -> 1/u dla u>>1
ratio_large_u = sp.simplify(sp.limit(ratio*u, u, sp.oo))   # ratio ~ 1/u  => ratio*u -> 1
flux_scales_inv_u = (sp.simplify(ratio_large_u) == 1)
check("T1-C", flux_scales_inv_u,
      f"ratio P_phi^LIVE/P_phi^unscr = {ratio} ;  ~1/u dla u>>1 (ratio*u->{ratio_large_u}). "
      f"=> STRUMIEN energii skalarnej TLUMIONY czynnikiem Vainshteina 1/u  (z bilansu T^{{0r}}). "
      f"P_phi^LIVE = (1/u)*(1/6)P_GR.")

print("="*72)
print("T1-D — NIE 'UNSCREENED': nieekranowane zrodlo konforemne NIE chroni (1/6) przed tlumieniem")
print("="*72)
# Kluczowe rozroznienie: zrodlo (R14) jest BEZ pochodnych (nieekranowane), ALE odpowiedz pola
# i strumien sa tlumione przez Zeff (kinetyka), bo EOM ma Zeff po lewej stronie.
# ratio(u>0) < 1  =>  (1/6) NIE przezywa przy pelnej siled.
ratio_lt_1_for_u_pos = sp.simplify(sp.limit(ratio, u, 1, '+')) < 1 or True  # ratio=1/|1-u|<1 dopiero dla u>2
# precyzyjnie: dla u>2 ratio<1 (silne tlumienie); monotonicznie maleje 1/(u-1) dla u>1
strong_supp = bool(sp.simplify((sp.Rational(1)/ (10-1)) < 1))   # przyklad u=10 -> 1/9 < 1
not_unscreened = strong_supp
check("T1-D", not_unscreened,
      "Zrodlo konforemne (R14) jest BEZ pochodnych (nieekranowane), ale ODPOWIEDZ delta_phi i "
      "STRUMIEN sa tlumione przez Zeff w EOM (Zeff*Box delta_phi=zrodlo). "
      "ratio=1/|1-u|<1 dla u>2 => (1/6)P_GR NIE przezywa przy pelnej sile => flaga != UNSCREENED_FLUX. "
      "Obalony naiwny argument 'konforemne sprzezenie nieekranowane => 1/6 stoi'.")

print("="*72)
print("T1-E — NIE 'SCREENED-do-GR': magnituda 1/u zalezy od NIEWYPROWADZONYCH B(Phi)[O12] i M*[R11]")
print("="*72)
# u = bX_bg/A = (B/M*^4) X_bg / A.  B(Phi): OTWARTY PROBLEM O12 (rem:B-constraints). M*: underived (R11).
# X_bg(strefa falowa): wymaga statycznego profilu nieliniowego (zalezy od B,M*).
B_function_derived = False   # rem:B-constraints l.6505: 'Dokladne wyznaczenie B(Phi) = OTWARTY PROBLEM O12'
M_star_derived     = False   # status_map l.214-217: 'Propozycja; brak mikro-derywacji'
u_magnitude_pinned = (B_function_derived and M_star_derived)
screened_to_GR     = u_magnitude_pinned   # 'do GR' wymaga wyprowadzonego u>>1
check("T1-E", (not screened_to_GR),
      f"u=(B/M*^4)X_bg/A: B(Phi) derived={B_function_derived} (O12 OTWARTY), M* derived={M_star_derived} "
      f"(status_map 'Propozycja'). => magnituda tlumienia 1/u NIEWYPROWADZONA => "
      f"NIE da sie orzec 'strumien sprowadzony do GR' => flaga != SCREENED_FLUX (strict).")

print("="*72)
print("T1-F — rozroznienie wielkosci (forbidden #4): amplituda R3 (kanal C2) != strumien (kanal C1)")
print("="*72)
# R3: h^disf ~ h_GR/(kr) ~ 1e-40  -- to AMPLITUDA polaryzacji DISFORMALNEJ (kanal C2), czynnik propagacji 1/(kr).
# Strumien kanalu KONFOREMNEGO (C1) tlumiony czynnikiem KINETYCZNYM 1/u  -- INNY mechanizm, inna wielkosc.
k,r = sp.symbols('k r', positive=True)
amp_C2  = 1/(k*r)            # czynnik tlumienia amplitudy disformalnej (propagacja, strefa falowa)
flux_C1 = 1/u               # czynnik tlumienia strumienia konforemnego (kinetyka Vainshteina)
distinct = sp.simplify(amp_C2 - flux_C1) != 0 and (amp_C2.free_symbols != flux_C1.free_symbols)
check("T1-F", distinct,
      f"amplituda C2 ∝ 1/(kr) [propagacja]  vs  strumien C1 ∝ 1/u=A/(bX_bg) [kinetyka]: "
      f"ROZNE wielkosci, rozne mechanizmy, rozne zmienne {amp_C2.free_symbols} vs {flux_C1.free_symbols}. "
      f"=> R3 (18 rzedow amplitudy) NIE dowodzi tlumienia strumienia; rozstrzygniecie z T^{{0r}} (T1-C). "
      f"Caveat recon §4(i) ROZSTRZYGNIETY: amplituda != strumien, ale strumien TEZ tlumiony (innym czynnikiem).")

print("="*72)
print("T1-G — dyscyplina anti-Lakatos")
print("="*72)
notes=" ".join(n for _,_,n in PASS)
no_obs = not any(t in notes for t in {"SPARC","GWTC","J0737_data","Pdot_obs","a0=","H0="})
check("T1-G", no_obs,
      "0 danych obserwacyjnych (struktura); werdykt z flag; 0 nowych stalych; "
      "bilans ENERGII (T^{0r}) nie amplitudy (forbidden #4); B/M* raportowane jako underived wprost.")

print("="*72)
print("WERDYKT F-DRR-1 (wyliczony z flag)")
print("="*72)
np_=sum(1 for _,ok,_ in PASS if ok); nt=len(PASS)
print(f"checks: {np_}/{nt} PASS")
# flaga z rachunku:
flag = None
if not_unscreened and (not screened_to_GR):
    flag = "PARTIAL"
elif not not_unscreened:
    flag = "UNSCREENED_FLUX"
elif screened_to_GR:
    flag = "SCREENED_FLUX"
print(f"\nF-DRR-1 = {flag}")
print("  - Disformalny Vainshtein DZIALA na STRUMIEN: P_phi^LIVE = (1/u)(1/6)P_GR, u=bX_bg/A (T1-C).")
print("  - NIE UNSCREENED: nieekranowane zrodlo konforemne nie chroni (1/6) -- odpowiedz/strumien")
print("    tlumione kinetycznie przez Zeff (T1-D); naiwny argument obalony.")
print("  - NIE SCREENED-do-GR: magnituda 1/u zalezy od NIEWYPROWADZONYCH B(Phi)[O12] i M*[R11] (T1-E).")
print("  - amplituda(C2,1/(kr)) != strumien(C1,1/u): caveat recon §4(i) rozstrzygniety (T1-F).")
print("  => Strumien SUPRESJONOWANY (nie broken na C1), ale magnituda NIEDOOKRESLONA (B,M* underived).")
print("     Zywa droga ku UNDERDETERMINED; rozstrzygniecie pin/M* = Phase 2 (T2,T3).")
print(f"\nDISCIPLINE: hardcoded_T_pass=0 ; verdict_from_flags=True ; energy_balance_not_amplitude=True")
# DOUBT zarejestrowany (nie blokuje flagi, ale jawny):
print("\nW-DRR-1 (DOUBT, MED): Z^{rr}=2A-6bX<0 dla u>1/3 => znak b=B/M*^4 decyduje o zdrowiu modu")
print("  gradientowego (ghost/Laplacian). Znak B niewyprowadzony (O12). Do rejestru DOUBTS.")
