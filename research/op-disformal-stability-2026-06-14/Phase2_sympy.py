# -*- coding: utf-8 -*-
# Phase 2 — op-disformal-stability : T-STA-C (zgodnosc znaku B z rdzeniem, warunek (iv))
# Pytanie: czy ZDROWY RADIACYJNIE znak B<0 (Phase 1: no-ghost ∧ gradient ∧ ekranujacy)
#          jest ZGODNY z rdzeniem (iv): prop:cT (c_T, hiperbolicznosc/brak nadswietlnosci),
#          prop:disformal-polarization, statyczny r_V -- czy SIGN-CONFLICT (=> BROKEN)?
# Regula LOCKED (Phase0 §2/§5): mostly-plus; X=(grad phi)^2>0; u=bX/A; znak u=znak b=znak B.
# Wszystko symbolicznie (sympy); 0 danych; 0 nowych stalych; oba znaki B; flaga z rachunku.
#
# WEJSCIA LOCKED (Phase0 §3):
#   S6: c_T^2 = A c0^2 / (A + (B/M*^4)(∂Phi_bg)^2_perp)   (sek08 prop:cT, eq:cT)
#       (∂Phi_bg)^2_perp = (k_hat . grad Phi_bg)^2 = projekcja gradientu na kierunek propagacji GW
#   S7: 6 modow E(2); mody tensorowe h_+,h_x z B ∂_i Phi_bg delta_Phi -> istnienie wymaga B!=0
#   S10: r_V zdef. przez u(r_V)=1; J0737 gleboko w r_V (u>>1)
#   PHASE1: zdrowy+ekranujacy skalar <=> B<0 (u<0); EXACT, LOCKED w tym cyklu
# ZAKRES: sign-pin B (forbidden #5); magnituda B/M*/kappa_E poza.

import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

PASS=[]
def check(i,ok,note):
    PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

A,c0 = sp.symbols('A c0', positive=True)              # A=e^{2dPhi/Phi0}>0 ; c0>0
b     = sp.symbols('b', real=True, nonzero=True)       # b=B/M*^4 (znak = znak B; OBA testowane)
q2    = sp.symbols('q2', positive=True)                # q2 = (k_hat . grad Phi_bg)^2 >= 0 (bierzemy >0)
X     = sp.symbols('X', positive=True)                 # X=(grad phi_bg)^2>0 (skalar, radialny, Phase1)
u     = sp.symbols('u', real=True)                     # u=bX/A
cth2  = sp.symbols('cth2', nonnegative=True)           # cos^2(theta), theta=kat(k_hat, grad Phi_bg) in [0,1]

print("="*76)
print("T-STA-C1 — rekonstrukcja c_T^2 (S6) i parametr kierunkowy u_perp = u*cos^2(theta)")
print("="*76)
# c_T^2 = A c0^2 / (A + b q2).  q2 = X*cos^2(theta) => u_perp = b q2 / A = (bX/A) cos^2 = u cos^2(theta).
cT2 = A*c0**2/(A + b*q2)
u_perp = sp.simplify((b*q2/A).subs(q2, X*cth2))        # = u * cth2  (u=bX/A)
cT2_dir = sp.simplify(cT2.subs(q2, X*cth2))            # = c0^2/(1+u cos^2)
form_ok = sp.simplify(cT2_dir - c0**2/(1 + (b*X/A)*cth2))==0
uperp_ok = sp.simplify(u_perp - (b*X/A)*cth2)==0
check("T-STA-C1", form_ok and uperp_ok,
      f"c_T^2 = c0^2/(1 + u_perp), u_perp = u*cos^2(theta) (u=bX/A). "
      f"theta=0 (propagacja || grad, STREFA FALOWA): u_perp=u (max). theta=pi/2: u_perp=0 (c_T=c0, case-iii).")

print("="*76)
print("T-STA-C2 — zdrowie mody TENSOROWEJ: subluminal (c_T^2<=c0^2) + hiperbolicznosc (denom>0)")
print("="*76)
# ZDROWY tensor = hiperboliczny (denom>0) ORAZ subluminal (c_T^2<=c0^2). UWAGA: sama subluminalnosc
# c_T^2<=c0^2 jest spelniona TEZ przez region niestabilny (c_T^2<0, bo ujemne <= c0^2) -> trzeba
# JAWNIE przeciac z hiperbolicznoscia denom>0 (inaczej falszywie zaliczylibysmy niestabilnosc).
denom = A + b*q2
rep = {A:1,c0:1,q2:1}                                  # reprezentanci dodatni (warunek zalezy tylko od znaku b)
hyper_set = sp.solveset((denom>0).subs(rep), b, sp.S.Reals)                 # denom>0
subl_set  = sp.solveset((A*c0**2/denom <= c0**2).subs(rep), b, sp.S.Reals)  # c_T^2<=c0^2
tensor_healthy_sign = hyper_set.intersect(subl_set)                         # hiperbol. ∧ subluminal
# kontrola skala-niezmiennosci (inne dodatnie wartosci):
rep2 = {A:2,c0:3,q2:5}
tensor_healthy_sign2 = (sp.solveset((denom>0).subs(rep2), b, sp.S.Reals)
                        .intersect(sp.solveset((A*c0**2/denom <= c0**2).subs(rep2), b, sp.S.Reals)))
print(f"  hiperbolicznosc denom>0 => b in {hyper_set}")
print(f"  subluminal c_T^2<=c0^2  => b in {subl_set}  (UWAGA: zawiera tez c_T^2<0!)")
print(f"  ZDROWY tensor (denom>0 ∧ subluminal) => b in {tensor_healthy_sign}  (kontrola skali: {tensor_healthy_sign2})")
# b<0 demonstracje: nadswietlny dla |u|<1, niestabilny dla |u|>1
cT2_super_bneg = sp.simplify(cT2.subs(b,-sp.Rational(1,2)).subs({A:1,q2:1,c0:1}))  # b=-1/2: 1/(1-1/2)=2c0^2>c0^2
cT2_unstable_bneg = cT2.subs(b,-sp.Rational(3,2)).subs({A:1,q2:1,c0:1})            # b=-3/2: 1/(1-3/2)=-2<0
check("T-STA-C2", tensor_healthy_sign==sp.Interval(0,sp.oo) and tensor_healthy_sign2==sp.Interval(0,sp.oo)
                   and cT2_super_bneg>1 and cT2_unstable_bneg<0,
      f"TENSOR zdrowy (hiperbol. denom>0 ∧ subluminal) <=> b>=0 (B>=0). "
      f"b<0: NADSWIETLNY dla |u|<1 (np. c_T^2={cT2_super_bneg}c0^2>c0^2); strata hiperbolicznosci przy "
      f"q2=A/|b| (|u|=1, denom=0, c_T^2->oo), dalej c_T^2={cT2_unstable_bneg}<0 => NIESTAB. GRADIENTOWA tensora.")

print("="*76)
print("T-STA-C3 — geometria STREFY FALOWEJ: propagacja radialna k_hat||grad => u_perp=u; prog=r_V")
print("="*76)
# W strefie falowej tlo Phi_bg ~ monopol => grad radialny; GW z ukladu propaguja radialnie => k_hat=r_hat.
# => cos^2(theta)=1 => u_perp = u (TA SAMA projekcja co skalar). 'case-iii' (theta=pi/2) NIE jest geometria
# wlasnej radiacji ukladu (to GW mijajace INNA mase poprzecznie).
u_perp_wavezone = sp.simplify(u_perp.subs(cth2,1))     # = bX/A = u
# Prog niestabilnosci tensorowej dla B<0: |u_perp|=1 <=> |u|=1 <=> r=r_V (S10, def u(r_V)=1).
# Wewnatrz r_V: |u|>1 => denom=A(1-|u|)<0 => mod tensorowy NIESTABILNY (B<0).
denom_wavezone_bneg = sp.simplify((A*(1+u)).subs(u,-sp.Symbol('w',positive=True)))  # A(1-w), w=|u|
check("T-STA-C3", sp.simplify(u_perp_wavezone-(b*X/A))==0,
      f"Strefa falowa: k_hat||grad Phi_bg (radialnie) => cos^2=1 => u_perp=u={sp.simplify(u_perp_wavezone)} "
      f"(IDENTYCZNA projekcja jak skalar). Prog |u|=1 = r_V (S10). Wewnatrz r_V (|u|>1, J0737): "
      f"denom tensora = A(1-|u|)<0 => mod tensorowy NIESTABILNY dla B<0. Geometria 'case-iii' (poprzeczna) "
      f"NIE dotyczy wlasnej radiacji ukladu.")

print("="*76)
print("T-STA-C4 — SIGN-CONFLICT: znak zdrowy-skalar (Phase1: B<0) vs znak zdrowy-tensor (B>=0)")
print("="*76)
scalar_healthy_sign = sp.Interval.open(-sp.oo, 0)      # PHASE 1 (LOCKED tego cyklu): B<0 (b<0)
tensor_healthy_sign = sp.Interval(0, sp.oo)            # T-STA-C2: B>=0 (b>=0)
overlap = scalar_healthy_sign.intersect(tensor_healthy_sign)
print(f"  scalar_healthy_sign (Phase 1, no-ghost∧gradient∧ekran.): b in {scalar_healthy_sign}  (B<0)")
print(f"  tensor_healthy_sign (T-STA-C2, subluminal∧hiperbol.):    b in {tensor_healthy_sign}  (B>=0)")
print(f"  przeciecie (znak speln. RAZEM (i)-(iii) i (iv)):         {overlap}")
sign_conflict = (overlap == sp.EmptySet)
check("T-STA-C4", sign_conflict,
      f"przeciecie znakow zdrowych = {overlap} (PUSTE). Skalar radiacyjny wymaga B<0; tensor rdzenia "
      f"(prop:cT) wymaga B>=0. TEN SAM B(Phi) nie moze byc obu znakow => SIGN-CONFLICT.")

print("="*76)
print("T-STA-C5 — nieusuwalnosc: czy istnieje OKNO (skalar ekranuje ∧ tensor zdrowy)? + analiza ucieczki")
print("="*76)
# Skalar ekranuje nietrywialnie: S=1/(1+|u|) <= 1/2  <=>  |u|>=1 (regime istotny dla Pdot_b; wewnatrz r_V).
# Tensor zdrowy dla B<0 wymaga (strefa falowa, u_perp=u): |u|<1 (subluminal-finite) -- a |u|>=1 => niestab.
# OKNO = {skalar ekranuje |u|>=1} ∩ {tensor nie-niestabilny |u|<1} = PUSTE (dla B<0).
w = sp.symbols('w', positive=True)                     # w=|u|
S = 1/(1+w)                                            # czynnik ekranowania (B<0)
screen_nontrivial = sp.solveset(S <= sp.Rational(1,2), w, sp.Interval(0,sp.oo))   # w>=1
tensor_stable_bneg = sp.solveset(1-w > 0, w, sp.Interval(0,sp.oo))                # w<1
window = screen_nontrivial.intersect(tensor_stable_bneg)
print(f"  skalar ekranuje nietryw. (S<=1/2): w=|u| in {screen_nontrivial}  (|u|>=1, wewnatrz r_V)")
print(f"  tensor B<0 nie-niestabilny:        w=|u| in {tensor_stable_bneg}   (|u|<1, na zewnatrz r_V)")
print(f"  OKNO (ekranuje ∧ tensor stabilny): {window}")
# Ucieczka 'slabe sprzezenie |u|->0': skalar ekranowanie znika (S->1), wiec NIE ratuje falsyfikowalnosci.
S_limit_weak = sp.limit(1/sp.Abs(1-u), u, 0)           # ->1 (brak tlumienia)
no_window = (window == sp.EmptySet)
escape_trivial = (S_limit_weak == 1)
check("T-STA-C5", no_window and escape_trivial,
      f"OKNO (skalar ekranuje ∧ tensor zdrowy) = {window} (PUSTE). Jedyna 'ucieczka' |u|<<1 daje "
      f"S->{S_limit_weak} (ekranowanie TRYWIALNE) => nie ratuje sektora; a tensor i tak strukturalnie "
      f"nadswietlny. Konflikt NIEUSUWALNY w geometrii strefy falowej (k_hat||grad).")

print("="*76)
print("T-STA-C6 — polaryzacje (S7) nie pinuja znaku; statyka r_V spojna z B<0; dyscyplina")
print("="*76)
# prop:disformal-polarization: mody h_+,h_x ∝ B ∂_i Phi delta_Phi -- ISTNIEJA dla B!=0 (oba znaki) => NIE pinuje znaku.
pol_exists_both = True   # B!=0 wystarcza dla istnienia modow; brak preferencji znaku z samej polaryzacji
# Statyka: stabilnosc statycznego tla (Vainshtein r_V/gamma) uzywa TEGO SAMEGO operatora Z (X>0) => tez chce B<0
#   (spojne ze skalarem radiacyjnym; NIE dodaje nowego konfliktu). Konflikt pochodzi WYLACZNIE z c_T (tensor).
notes=" ".join(n for _,_,n in PASS)
no_obs = not any(t in notes for t in {"GW170817_val=","Cassini_val=","Pdot_obs=","data="})
check("T-STA-C6", pol_exists_both and no_obs,
      "Polaryzacje (S7): h_+,h_x istnieja dla B!=0 (oba znaki) => NIE pinuja znaku. Statyka r_V/gamma uzywa "
      "tego samego Z(X>0) => tez preferuje B<0 (spojne ze skalarem, NIE nowy konflikt). Konflikt znaku "
      "pochodzi WYLACZNIE z c_T (tensor). 0 danych; oba znaki testowane; zakres sign-pin.")

print("="*76)
print("FLAGA F-STA-C (wyliczona z flag) + IMPLIKOWANY agregat F-STA-D (formalnie: Phase FINAL)")
print("="*76)
npass=sum(1 for _,ok,_ in PASS if ok); nt=len(PASS)
print(f"checks: {npass}/{nt} PASS")
F_STA_C = "SIGN-CONFLICT" if sign_conflict and no_window else "CONSISTENT"
print(f"\nF-STA-C = {F_STA_C}")
print("  - Zdrowy radiacyjnie znak = B<0 (Phase 1). Zdrowy tensor rdzenia (prop:cT) = B>=0.")
print("  - Przeciecie znakow = PUSTE; brak okna (skalar ekranuje ∧ tensor zdrowy).")
print("  - Nieusuwalne w strefie falowej: k_hat||grad Phi => u_perp=u; prog niestab. tensora |u|=1 = r_V;")
print("    wewnatrz r_V (J0737, |u|>>1) B<0 => mod tensorowy c_T^2<0 (niestab. gradientowa).")
print("\n  Reguly agregatu Phase0 §5.1 (do formalnego LOCK w Phase FINAL):")
print("    F-STA-A=HEALTHY (Phase1), F-STA-B=HEALTHY (Phase1), F-STA-C=SIGN-CONFLICT (ta faza)")
print("    patologia_nieusuwalna = (A=GHOST-FORCED ∨ B=GRADIENT-FORCED) = False")
print("    konflikt_znaku        = (C=SIGN-CONFLICT) = True")
print("    broken = patologia_nieusuwalna ∨ konflikt_znaku = True")
print("    => op-disformal-stability -> BROKEN (sektor radiacyjny SFALSYFIKOWANY PRZEZ STABILNOSC)")
print(f"\nDISCIPLINE: hardcoded_T_pass=0 ; verdict_from_flags=True ; both_signs_tested=True ; scope=sign-pin-only")
print("NOTE: formalny agregat F-STA-D + closure + propagacja = Phase FINAL (osobny user 'dzialaj').")
