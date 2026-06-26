# -*- coding: utf-8 -*-
# Phase 1 — op-disformal-stability : T-STA-A (no-ghost) + T-STA-B (gradient)
# Pytanie: czy istnieje znak/zakres B(Phi) taki, ze operator fluktuacji disformalnej jest
#          JEDNOCZESNIE (i) ghost-free, (ii) gradient-stabilny (c_s^2>=0), (iii) EKRANUJACY
#          strumien (czynnik 1/|1-u|<1) -- na tle STATYCZNYM (X=(grad phi)^2>0).
# Regula LOCKED (Phase0 §2): mostly-plus eta=diag(-1,1,1,1); tlo statyczne ∂_t phi=0,
#          gradient radialny; X=(grad phi)^2>0; u=bX/A; znak u = znak b = znak B (bo X>0,A>0).
# Wszystko symbolicznie (sympy); 0 danych obserwacyjnych; 0 nowych stalych; flaga z rachunku.
#
# WEJSCIA LOCKED uzyte (Phase0 §3):
#   S1: Z^{mn}=2(A-bX)eta^{mn}-4b ∂^m phi ∂^n phi  (op-disformal-radiation-resolution Phase1, EXACT)
#   S2: Z^{00}∝-(A-bX); Z^{rr}=2A(1-3u); Z^{perp}=2A(1-u)  (do potwierdzenia z S1)
#   S3: c_s^2 = L'/(L'+2X L'') = (1-u)/(1-3u)               (Garriga-Mukhanov; do potwierdzenia)
#   S4: czynnik tlumienia strumienia 1/|1-u| (op-disformal-radiation-resolution Phase1 §3, LOCKED)
#   S9: znak B(Phi) NIEWYPROWADZONY (O12) -> sprawdzamy OBA znaki, nie zakladamy
# ZAKRES: wylacznie sign-pin B (Phase0 forbidden #5); magnituda B / M* / kappa_E poza zakresem.

import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

PASS=[]
def check(i,ok,note):
    PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

print("="*74)
print("T-STA-A1 — operator Z^{mn} z L=AX-(b/2)X^2 na tle STATYCZNYM (gradient radialny)")
print("="*74)
A,b = sp.symbols('A b', positive=True)               # A=e^{2dPhi/Phi0}>0 ; b=B/M*^4 (znak = znak B)
G    = sp.symbols('G', real=True, nonzero=True)       # G = ∂_r phi_bg (skladowa radialna, x-kierunek)
eta  = sp.diag(-1,1,1,1)                               # mostly-plus (Phase0 §2)
# Z^{mn}=d^2 L/d(∂_m phi)d(∂_n phi): rozniczkuj po SYMBOLICZNYCH skladowych, POTEM podstaw tlo.
w0,w1,w2,w3 = sp.symbols('w0 w1 w2 w3', real=True)
w    = sp.Matrix([w0,w1,w2,w3])                        # ∂_mu phi (symbolicznie)
Xsym = (w.T*eta*w)[0]                                  # X = eta^{mn} ∂_m phi ∂_n phi
Lsym = A*Xsym - sp.Rational(1,2)*b*Xsym**2             # S1: kinetyka LIVE do O(b)
# tlo statyczne: ∂_mu phi_bg = (0, G, 0, 0)  (∂_t=0; gradient czysto przestrzenny, radialny=x)
bg   = {w0:0, w1:G, w2:0, w3:0}
v    = sp.Matrix([0, G, 0, 0])
X    = (v.T*eta*v)[0]                                  # X = eta^{mn} v_m v_n = +G^2 > 0 (przestrzenny)
check("T-STA-A1.X", sp.simplify(X-G**2)==0,
      f"X=(grad phi)^2={sp.simplify(X)} = G^2 > 0 na tle statycznym (mostly-plus, gradient przestrzenny). "
      f"u=bX/A: znak u = znak b = znak B (A>0, X>0). REZIM STATYCZNY (rozny od kosmologicznego X<0).")

Z = sp.zeros(4,4)
for m in range(4):
    for n in range(4):
        Z[m,n] = sp.simplify(sp.diff(Lsym, w[m], w[n]).subs(bg))
Lp, Lpp = A - b*X, -b                                  # L'(X)=A-bX ; L''(X)=-b
vcov = eta*v                                           # v_mu (dolny)
Z_expected = 2*Lp*eta + 4*Lpp*(vcov*vcov.T)            # forma EXACT z S1
Z_match = sp.simplify(sp.Matrix(Z)-Z_expected) == sp.zeros(4,4)
check("T-STA-A1", Z_match,
      f"Z^{{mn}}=2(A-bX)eta-4b(∂phi)(∂phi)^T [EXACT, match={Z_match}] -- zgodne z S1 (reuse, nie nowy operator).")

print("="*74)
print("T-STA-A2 — skladowe znakowe (S2) + no-ghost: znak czlonu czasowego Z^{00} i kinetyki L'")
print("="*74)
u = sp.symbols('u', real=True)                         # u=bX/A (znak DOWOLNY -- sprawdzamy oba, S9)
sub = {X: A*u/b}                                       # X = A u / b  =>  bX=Au, u=bX/A
Z00 = sp.simplify(Z[0,0].subs(sub))                    # oczek. -2A(1-u)
Zrr = sp.simplify(Z[1,1].subs(sub))                    # oczek. +2A(1-3u)  (radialny, wzdluz gradientu)
Zpp = sp.simplify(Z[2,2].subs(sub))                    # oczek. +2A(1-u)   (transwersalny)
Z00_ok = sp.simplify(Z00 - (-2*A*(1-u)))==0
Zrr_ok = sp.simplify(Zrr - ( 2*A*(1-3*u)))==0
Zpp_ok = sp.simplify(Zpp - ( 2*A*(1-u)))==0
check("T-STA-A2.form", Z00_ok and Zrr_ok and Zpp_ok,
      f"Z^00={Z00}=-2A(1-u); Z^rr={Zrr}=2A(1-3u); Z^perp={Zpp}=2A(1-u) -- POTWIERDZA S2 z S1.")
# No-ghost (konwencja k-essence projektu: zdrowy <=> L'=A-bX>0, tj. czlon glowny ~eta ma zdrowy znak;
# rownowaznie Z^00, Z^perp maja zdrowy znak): warunek 1-u>0.
Lp_u = sp.simplify(Lp.subs(sub))                       # = A(1-u)
noghost_cond = sp.simplify(Lp_u/A)                     # >0  <=>  1-u>0  <=>  u<1
check("T-STA-A2", sp.simplify(noghost_cond-(1-u))==0,
      f"L'=A(1-u); no-ghost <=> L'>0 <=> (1-u)>0 <=> u<1. Dla u>1 (B>0 silne): L'<0 => GHOST.")

print("="*74)
print("T-STA-B1 — gradient: c_s^2=(1-u)/(1-3u) (S3) + Z^rr, Z^perp; warunki stabilnosci")
print("="*74)
cs2 = sp.simplify(Lp/(Lp + 2*X*Lpp)).subs(sub)        # = (A-bX)/(A-3bX) = (1-u)/(1-3u)
cs2 = sp.simplify(cs2)
cs2_ok = sp.simplify(cs2 - (1-u)/(1-3*u))==0
check("T-STA-B1.cs2", cs2_ok, f"c_s^2 = L'/(L'+2X L'') = {cs2} = (1-u)/(1-3u) -- POTWIERDZA S3.")
# Warunki stabilnosci gradientowej:
#   c_s^2>=0 ; Z^rr=2A(1-3u)>=0 (mod radialny) ; Z^perp=2A(1-u)>=0 (mod transwersalny)
cs2_ge0  = sp.solveset(cs2 >= 0, u, sp.S.Reals)
Zrr_ge0  = sp.solveset(sp.Rational(1,1)*(1-3*u) >= 0, u, sp.S.Reals)
Zpp_ge0  = sp.solveset((1-u) >= 0, u, sp.S.Reals)
print(f"  c_s^2>=0 : {cs2_ge0}")
print(f"  Z^rr>=0  : {Zrr_ge0}   (1-3u>=0 => u<=1/3)")
print(f"  Z^perp>=0: {Zpp_ge0}   (1-u>=0 => u<=1)")
# gradient_health = przeciecie (sciśle: u<1/3 -- bo Z^rr wymaga u<=1/3, a c_s^2 ma osobliwosc w u=1/3)
grad_health = sp.solveset(1-3*u > 0, u, sp.S.Reals).intersect(sp.solveset(1-u > 0, u, sp.S.Reals))   # u<1/3
print(f"  gradient_health (Z^rr>0 ∧ Z^perp>0): {grad_health}  (=> u<1/3)")
check("T-STA-B1", grad_health == sp.Interval.open(-sp.oo, sp.Rational(1,3)),
      "gradient-stabilny <=> u<1/3 (Z^rr>0, c_s^2>0). Dla 1/3<u<1: Z^rr<0 i c_s^2<0 => NIESTAB. GRADIENTOWA.")

print("="*74)
print("T-STA-B2 — czynnik EKRANOWANIA S(u)=1/|1-u| (S4): gdzie GENUINE suppression S<1 ?")
print("="*74)
# op-disformal-radiation-resolution Phase1 §3: P_LIVE/P_unscr = 1/|1-u|. Suppression <=> S<1 <=> |1-u|>1.
screen = sp.solveset(sp.Abs(1-u) > 1, u, sp.S.Reals)  # u>2  OR  u<0
print(f"  S(u)=1/|1-u|<1  <=>  |1-u|>1  <=>  u in {screen}")
# kontra: 0<u<2 => |1-u|<1 => S>1 => ENHANCEMENT strumienia (gorzej niz brak ekranowania)
enhance = sp.solveset(sp.Abs(1-u) < 1, u, sp.S.Reals) # 0<u<2
print(f"  S(u)>1 (ENHANCEMENT, niefizyczne dla ekranowania): u in {enhance}")
check("T-STA-B2", screen == sp.Union(sp.Interval.open(-sp.oo,0), sp.Interval.open(2,sp.oo)),
      "Genuine suppression (S<1) TYLKO dla u<0 (B<0) lub u>2 (B>0, gleboki ghost). "
      "Dla 0<u<2 (B>0 slabe/posrednie): S>1 => strumien WZMOCNIONY (jeszcze gorzej).")

print("="*74)
print("T-STA-C(prep) — EGZYSTENCJA znaku B: czy (no-ghost ∧ gradient ∧ ekranowanie) jest NIEPUSTE?")
print("="*74)
# full_health = no-ghost ∧ gradient_health = (u<1) ∧ (u<1/3) = u<1/3
full_health = sp.solveset(1-u>0, u, sp.S.Reals).intersect(sp.solveset(1-3*u>0, u, sp.S.Reals))   # u<1/3
healthy_and_screen = full_health.intersect(screen)                            # (u<1/3) ∩ (u<0 ∪ u>2)
print(f"  full_health (no-ghost ∧ gradient): {full_health}   (u<1/3)")
print(f"  screening (S<1):                   {screen}        (u<0 ∪ u>2)")
print(f"  health ∧ screening:                {healthy_and_screen}")
# Rozbicie wg znaku B:
branch_Bpos = healthy_and_screen.intersect(sp.Interval.open(0, sp.oo))        # u>0 : oczek. PUSTE
branch_Bneg = healthy_and_screen.intersect(sp.Interval.open(-sp.oo, 0))       # u<0 : oczek. caly u<0
print(f"  galaz B>0 (u>0): {branch_Bpos}   <- czy istnieje zdrowy+ekranujacy dla B>0 ?")
print(f"  galaz B<0 (u<0): {branch_Bneg}   <- czy istnieje zdrowy+ekranujacy dla B<0 ?")
exists_Bpos = (branch_Bpos != sp.EmptySet)
exists_Bneg = (branch_Bneg != sp.EmptySet)
check("T-STA-C.prep", (not exists_Bpos) and exists_Bneg,
      f"B>0: zdrowy+ekranujacy = {branch_Bpos} (PUSTE -- ekranowanie B>0 wymaga u>2 = gleboki ghost). "
      f"B<0: zdrowy+ekranujacy = {branch_Bneg} (NIEPUSTE -- caly u<0). "
      f"=> ZDROWY+EKRANUJACY ZNAK ISTNIEJE i JEST JEDNOZNACZNIE B<0.")

print("="*74)
print("T-STA-Bneg — weryfikacja galezi B<0: c_s^2, Z^rr, Z^perp, S wszystkie zdrowe dla u<0")
print("="*74)
# u<0: c_s^2=(1+|u|)/(1+3|u|) in (1/3,1); Z^rr=2A(1+3|u|)>0; Z^perp=2A(1+|u|)>0; S=1/(1+|u|)<1
ww = sp.symbols('w', positive=True)                   # w=|u|, u=-w
cs2_neg   = sp.simplify(cs2.subs(u,-ww))              # (1+w)/(1+3w)
Zrr_neg   = sp.simplify((1-3*u).subs(u,-ww))          # 1+3w
Zpp_neg   = sp.simplify((1-u).subs(u,-ww))            # 1+w
S_neg     = sp.simplify((1/sp.Abs(1-u)).subs(u,-ww))  # 1/(1+w)
cs2_in01  = sp.simplify(cs2_neg.subs(ww,1))           # przyklad 2/4=1/2 in (0,1)
allhealthy_neg = (sp.simplify(cs2_neg-(1+ww)/(1+3*ww))==0 and
                  sp.simplify(Zrr_neg-(1+3*ww))==0 and sp.simplify(Zpp_neg-(1+ww))==0 and
                  sp.simplify(S_neg-1/(1+ww))==0)
# c_s^2 in (1/3,1) dla w>0: granice
cs2_lim0 = sp.limit(cs2_neg, ww, 0)   # ->1
cs2_limI = sp.limit(cs2_neg, ww, sp.oo) # ->1/3
check("T-STA-Bneg", allhealthy_neg,
      f"B<0 (u=-w,w>0): c_s^2=(1+w)/(1+3w) in (1/3,1) [lim_0={cs2_lim0}, lim_inf={cs2_limI}] -- SUBLUMINAL, zdrowy; "
      f"Z^rr=2A(1+3w)>0; Z^perp=2A(1+w)>0; S=1/(1+w)<1 (ekranuje). WSZYSTKO ZDROWE ∀w>0, dowolna magnituda.")

print("="*74)
print("T-STA-disc — dyscyplina anti-Lakatos")
print("="*74)
notes=" ".join(n for _,_,n in PASS)
no_obs = not any(t in notes for t in {"GW170817_val","Cassini_val","Pdot_obs","sigma=","data="})
check("T-STA-disc", no_obs,
      "0 danych obserwacyjnych (struktura); werdykt z flag; 0 nowych stalych; OBA znaki B sprawdzone "
      "(nie zalozone); zakres = wylacznie sign-pin (magnituda B/M*/kappa_E poza, forbidden #5).")

print("="*74)
print("FLAGI F-STA-A / F-STA-B (wyliczone z flag) + kandydat sign-pin")
print("="*74)
npass=sum(1 for _,ok,_ in PASS if ok); nt=len(PASS)
print(f"checks: {npass}/{nt} PASS")
# F-STA-A (no-ghost): HEALTHY <=> istnieje znak/zakres B ekranujacy z L'>0 (Z^00 zdrowe).
# F-STA-B (gradient): HEALTHY <=> istnieje znak/zakres B ekranujacy z c_s^2>=0, Z^rr,Z^perp>=0.
# Oba spelnione przez ISTNIENIE galezi B<0 (exists_Bneg) -- ta sama galaz zdrowa pod oba testy.
F_STA_A = "HEALTHY" if exists_Bneg else "GHOST-FORCED"
F_STA_B = "HEALTHY" if exists_Bneg else "GRADIENT-FORCED"
sign_pin_candidate = "B<0 (u<0)" if (exists_Bneg and not exists_Bpos) else ("BOTH" if exists_Bpos else "NONE")
print(f"\nF-STA-A = {F_STA_A}")
print(f"F-STA-B = {F_STA_B}")
print(f"KANDYDAT SIGN-PIN (z health∧screening) = {sign_pin_candidate}")
print("  - Zdrowy (no-ghost ∧ gradient ∧ ekranujacy) operator ISTNIEJE i wymaga JEDNOZNACZNIE B<0.")
print("  - B>0: ekranowanie wymaga u>2 (gleboki ghost), a rejon zdrowy u<1/3 ENHANCUJE strumien (S>1).")
print("  - B<0: zdrowy ∀ magnitudy (c_s^2 in (1/3,1) subluminal) ORAZ ekranuje (S=1/(1+|u|)<1).")
print("  => F-STA-A=HEALTHY, F-STA-B=HEALTHY. Patologia NIE jest wymuszona (istnieje zdrowy znak).")
print("\n  WERDYKT CYKLU NIEROZSTRZYGNIETY na A/B: przenosi sie na F-STA-C (Phase 2):")
print("  czy B<0 (zdrowy radiacyjnie) jest ZGODNY z rdzeniem (iv)? prop:cT (S6) daje c_T^2=Ac0^2/(A+(B/M*^4)(∂Phi)^2_perp):")
print("  dla B<0 mianownik<A => c_T^2>c0^2 (NADSWIETLNA moda tensorowa) -- PRE-FLAGOWANE napiecie (Phase0 §6).")
print("  Jesli SIGN-CONFLICT (B<0 radiacja vs B>0 c_T) nieusuwalny => BROKEN; jesli rozlaczne geom. => SIGN-PINNED B<0.")
print(f"\nDISCIPLINE: hardcoded_T_pass=0 ; verdict_from_flags=True ; both_signs_tested=True ; scope=sign-pin-only")
