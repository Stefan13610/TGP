# -*- coding: utf-8 -*-
# AUDYT WERDYKTU (reviewer, Claudian 2026-06-14) — op-disformal-stability
# Cel: sprawdzic argument Phase 2 (BROKEN via prop:cT induced-TT) i podac wlasciwy obiekt.
# Trzy rachunki:
#  (1) fizyczna dyspersja skalara z Z^{mn} (proper k-essence) -> zdrowy dla B<0
#  (2) induced-TT c_T z prop:cT (naive operator) -> c_T^2<0 dla B<0,|u|>1 (ale to NIE-DOF)
#  (3) g_eff signatura/wyznacznik -> DISFORMAL VIABILITY: flip przy |u|=1 dla B<0 (wlasciwy argument)
import sympy as sp
A,G,b,c = sp.symbols('A G b c', positive=True)
u,W = sp.symbols('u W', real=True)
eta = sp.diag(-1,1,1,1); v = sp.Matrix([0,G,0,0])

print("="*68); print("(1) FIZYCZNY skalar: Z^{mn} k_m k_n = 0 (radialnie, k||grad)"); print("="*68)
w,kr = sp.symbols('omega k_r', positive=True)
eta_kk = -(w/c)**2 + kr**2; dphi_k = G*kr
expr = 2*(A-b*G**2)*eta_kk - 4*b*dphi_k**2
sol = sp.solve(sp.Eq(expr,0), w**2)[0]
cs2 = sp.simplify((sol/kr**2)/c**2).subs(b, u*A/G**2)
print("  c_scalar^2/c0^2 =", sp.simplify(cs2), " ; B<0 (u=-W):", sp.simplify(cs2.subs(u,-W)))
print("  -> dla B<0 ZAWSZE >0 (zdrowy, choc nadswietlny). Argument 'scalar wants B<0' POPRAWNY.")

print("="*68); print("(2) induced-TT c_T (prop:cT, NAIVE operator A*box + b dPhi dPhi dd)"); print("="*68)
print("  c_T^2/c0^2 = 1+u (proof-form);  B<0 (u=-W):", sp.simplify(1-W), " (<0 dla W>1)")
print("  ALE rem:GW-scope-2026: dPhi^TT to NIE niezalezny spin-2 DOF -> c_T^2<0 NIE jest")
print("  niestabilnoscia propagujacego modu. To obiekt slaved. Argument Phase 2 NIEROBUSTNY.")

print("="*68); print("(3) g_eff VIABILITY (wlasciwy argument): det + sygnatura"); print("="*68)
g = sp.simplify(A*eta + b*(v*v.T))
det = sp.factor(g.det())
print("  g_eff = diag(-A, A+bG^2, A, A) ;  det =", det, "= -A^4(1+u)")
print("  eigenvalues:", g.eigenvals())
print("  radial eig = A(1+u): B<0 (u=-W) -> A(1-W); <0 dla W>1 => SYGNATURA FLIP (drugi czas).")
print("  => g_eff Lorentzowska <=> 1+u>0 <=> (B>0 dowolne) lub (B<0, |u|<1).  Prog |u|=1 = r_V.")

print("="*68); print("WNIOSEK AUDYTU"); print("="*68)
print("  screening S=1/|1-u| male <=> |u|>>1.")
print("  B<0,|u|>>1: skalar zdrowy ALE g_eff sygnatura flip (niefizyczne).")
print("  B>0,|u|>2 (screening): L'=A(1-u)<0 => ghost skalara.")
print("  B<0,|u|<1: g_eff OK, ale S~1 (brak screeningu) => PR-025 konforemny stoi.")
print("  => (g_eff Lorentz)  AND  (skalar zdrowy)  AND  (strong screening) = ZBIOR-PUSTY.")
print("  WERDYKT: BROKEN PRAWDOPODOBNIE POPRAWNY, ale via DISFORMAL VIABILITY, nie via induced-TT.")
print("  Prog |u|=1, ktory Phase 2 przypisala 'niestab. tensora', to NAPRAWDE degeneracja g_eff.")
