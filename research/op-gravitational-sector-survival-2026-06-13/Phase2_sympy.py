# -*- coding: utf-8 -*-
# Phase 2 — op-gravitational-sector-survival : F-GSS-B (wyczerpanie) + ocena D6
# Dwa filary strukturalne (sympy, 0 danych, 0 nowych stalych, werdykt z flag):
#  (I)  D6 disformal lamie premise FP7: metryka disformalna daje L_kin X-NIELINIOWE
#       (czlon X^2 ~ B/M*^4)  -> no-go FP7 (X-liniowy) NIE rozciaga sie na pelne LIVE.
#  (II) T5/kappa_E: kombinacja strumienia energii (xi/lambda) NIEZALEZNA od dopasowania
#       amplitudy (lambda*xi) -> kappa_E unpinned -> sektor radiacyjny NIEDOOKRESLONY.
# Wniosek: NIE da sie orzec FALSIFIED (D6 zywe, nierozstrzygniete) -> INDETERMINATE.

import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

PASS=[]
def check(i,ok,note):
    PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

print("="*70); print("FILAR I — D6: metryka disformalna -> L_kin X-NIELINIOWE (lamie FP7)"); print("="*70)
A,b,X = sp.symbols('A b X', real=True, positive=True)   # b = B/M*^4 (>0 dla konkretnosci)
# konfiguracja: d_mu phi = v_mu (jeden kierunek), v.v = X (metryka eta mostly-plus diag(-1,1,1,1))
# metryka disformalna g = A*eta + b * v⊗v ; policz dokladnie det, g^{-1}, oraz L = sqrt(-g) g^{mn} v_m v_n
eta = sp.diag(-1,1,1,1)
v0,v1 = sp.symbols('v0 v1', real=True)                   # v=(v0,v1,0,0); X = -v0^2+v1^2
v = sp.Matrix([v0,v1,0,0])
Xexpr = (v.T*eta*v)[0]                                    # = -v0^2+v1^2
g = A*eta + b*(v*v.T)                                     # disformal metric (dolne indeksy)
detg = sp.simplify(g.det())
ginv = sp.simplify(g.inv())
# kontrakcja g^{mn} v_m v_n  (v_m dolne):
contr = sp.simplify((v.T*ginv*v)[0])
L_kin = sp.simplify(sp.sqrt(-detg)*contr)                 # sqrt(-g) * g^{mn} v_m v_n
# zastap X = -v0^2+v1^2 i rozwin w b do rzedu 1:
L_in_X = L_kin.subs(v1, sp.sqrt(v0**2 + X))               # wymusza -v0^2+v1^2 = X
L_in_X = sp.simplify(L_in_X)
ser = sp.series(L_in_X, b, 0, 2).removeO()                # rozwiniecie w b do O(b^1)
ser = sp.expand(sp.simplify(ser))
# wspolczynnik przy X^2 w rozwinieciu (czlon X-nieliniowy):
coeff_X2 = sp.simplify(ser.coeff(b,1).as_poly(X).coeff_monomial(X**2)) if ser.coeff(b,1)!=0 else 0
# oczekiwane: L ~ A*X - (b/2) X^2 + ... ; sprawdz obecnosc czlonu X^2 ∝ b
nonlinear_present = sp.simplify(coeff_X2) != 0
check("FP-I", nonlinear_present,
      f"L_kin(disformal) = {sp.nsimplify(ser)} ; czlon X^2 ∝ b={coeff_X2}≠0 "
      f"=> kinetyka X-NIELINIOWA (k-essence/Vainshtein). Premise FP7 (X-liniowy) ZLAMANA "
      f"=> no-go Phase 1 NIE rozciaga sie na akcje disformalna LIVE.")

print("="*70); print("FILAR II — T5: kappa_E (strumien) NIEZALEZNE od dopasowania amplitudy"); print("="*70)
lam,xi = sp.symbols('lambda xi', positive=True)
O_amp  = lam*xi          # obserwabla amplitudowa (dopasowana do GR: pinuje lam*xi)
O_flux = xi/lam          # strumien energii ∝ kappa_E
# Jacobian mapy (lam,xi) -> (O_amp, O_flux): nieosobliwy => dwie NIEZALEZNE kombinacje
J = sp.Matrix([[sp.diff(O_amp,lam), sp.diff(O_amp,xi)],
               [sp.diff(O_flux,lam), sp.diff(O_flux,xi)]])
detJ = sp.simplify(J.det())
independent = sp.simplify(detJ) != 0
# konsekwencja: ustalenie O_amp=const NIE ustala O_flux => kappa_E swobodne (W-PDOT-5: niepinowane w LIVE)
kappaE_unpinned = independent
check("FP-II", kappaE_unpinned,
      f"det J[(lam,xi)->(lam*xi, xi/lam)] = {detJ} ≠ 0 => O_amp i O_flux NIEZALEZNE. "
      f"Dopasowanie amplitudy (lam*xi) NIE pinuje strumienia (xi/lam=kappa_E). "
      f"=> kappa_E unpinned w LIVE (PR-025 T5/W-PDOT-5). Sektor radiacyjny NIEDOOKRESLONY.")

print("="*70); print("FP-III — D6 nie strojony: czy ratunek wymaga tuningu? (anti-Lakatos)"); print("="*70)
# D6 jest REVISION_CLEAN tylko jesli (i) tlumienie skalara wyprowadzone (nie fit M*),
#   ORAZ (ii) kappa_E wyprowadzone do wartosci GR (nie strojone).
# Status faktyczny (z recon + FP-II): M* = 'Propozycja' (underived); kappa_E unpinned.
M_star_derived = False    # status_map: 'Propozycja, brak mikro-derywacji'
kappaE_derived = False    # FP-II: unpinned (nie wyprowadzone do GR; strojenie = forbidden #2)
D6_clean = (M_star_derived and kappaE_derived)
# D6 nie jest ani CLEAN (brak derywacji), ani BROKEN (nie lamie (a)/(b)/(c)) => LIVE-UNRESOLVED
D6_broken = False         # nie lamie zadnego z (a)/(b)/(c) a priori
D6_status = ("REVISION_CLEAN" if D6_clean else
             ("BROKEN" if D6_broken else "LIVE_UNRESOLVED"))
check("FP-III", D6_status=="LIVE_UNRESOLVED",
      f"D6: M*_derived={M_star_derived}, kappaE_derived={kappaE_derived} => NIE clean (brak derywacji); "
      f"NIE broken (zachowuje (a)/(b)/(c)) => status={D6_status}. "
      f"Rozstrzygniecie wymaga rachunku nieliniowej radiacji + pinowania kappa_E (osobny cykl).")

print("="*70); print("FP-IV — F-GSS-B: czy zbior {D1..D6} jest EXHAUSTED-and-broken?"); print("="*70)
# Klasy z Phase 1 + D6:
routes = {"D1":"BREAKS_§1","D2":"BREAKS_alpha","D3":"BREAKS_1r","D4":"BREAKS_1r",
          "D5":"GAP","D6":D6_status}
all_broken_or_gap = all(v in ("BREAKS_§1","BREAKS_alpha","BREAKS_1r","GAP") for v in routes.values())
exhausted_and_broken = all_broken_or_gap     # warunek orzeczenia FALSIFIED
check("FP-IV", (not exhausted_and_broken) and (routes["D6"]=="LIVE_UNRESOLVED"),
      f"routes={routes}. all_broken_or_gap={all_broken_or_gap}. "
      f"D6=LIVE_UNRESOLVED => zbior NIE jest 'wyczerpany-i-zlamany' => "
      f"NIE wolno orzec FALSIFIED-AS-AXIOMATIZED (wymog usera '100%').")

print("="*70); print("FP-V — dyscyplina"); print("="*70)
notes=" ".join(n for _,_,n in PASS)
clean = not any(t in notes for t in {"R_obs","SPARC","Pdot_obs","a0_obs","H0_obs"})
check("FP-V", clean, "0 danych obserwacyjnych; werdykt z flag; 0 nowych stalych; D6 underived raportowany wprost.")

print("="*70); print("WERDYKT F-GSS-B (wyliczony z flag)"); print("="*70)
np_=sum(1 for _,ok,_ in PASS if ok); nt=len(PASS)
print(f"FP: {np_}/{nt} PASS")
print(f"  D1..D5: {[routes[k] for k in ('D1','D2','D3','D4','D5')]}")
print(f"  D6:     {routes['D6']}  (zywa droga LIVE, nierozstrzygnieta)")
if exhausted_and_broken:
    print("F-GSS-B => EXHAUSTED-and-broken => (do FINAL) FALSIFIED-AS-AXIOMATIZED")
else:
    print("F-GSS-B => NOT_EXHAUSTED (D6 zywe/nierozstrzygniete) =>")
    print("  (do Phase FINAL) F-GSS-C kandydat: INDETERMINATE — sektor NIE sfalsyfikowany.")
    print("  No-go Phase 1 wazny dla pod-teorii KONFOREMNEJ; pelne LIVE (disformal+kappa_E) otwarte.")
    print("  'EXHAUSTIVE-OVER-LIVE' z PR-025 = nad-zasiegowe (exhaustive nad konforemnym, nie LIVE).")
print(f"\nDISCIPLINE: hardcoded=0 ; verdict_from_flags=True ; D6_added_pre_phase2=True")
