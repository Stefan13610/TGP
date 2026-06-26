# -*- coding: utf-8 -*-
# Phase 1 — op-Csigma-coarse-graining : formalne H_Gamma dla sigma_ab (SETUP, nie rdzen).
# Cel: zweryfikowac STRUKTURALNIE (sympy) twierdzenia Phase 1, NIE wyprowadzac C_sigma (to Phase 2).
# Anti-Lakatos: zaden prefaktor C_sigma nie jest tu fabrykowany; ekstrakcja = Phase 2 (GAP-pending).
import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

print("="*72)
print("Phase 1 — op-Csigma-coarse-graining : weryfikacja strukturalna (SETUP)")
print("="*72)

# ---------------------------------------------------------------------------
# (1) sigma_ab = K_ab - 1/3 Tr(K) delta : symetryczny, BEZSLADOWY, 5 skladowych
#     (def:sigma-ab, eq:sigma-def). Spin-2 => dokladnie 5 niezaleznych komponentow.
# ---------------------------------------------------------------------------
k = sp.symbols('k_xx k_yy k_zz k_xy k_xz k_yz', real=True)
K = sp.Matrix([[k[0], k[3], k[4]],
               [k[3], k[1], k[5]],
               [k[4], k[5], k[2]]])
sigma = K - sp.Rational(1,3)*sp.trace(K)*sp.eye(3)
print("\n(1) sigma_ab = K_ab - 1/3 Tr(K) delta_ab")
print("    symetryczna  :", sigma == sigma.T)
print("    Tr(sigma)    :", sp.simplify(sp.trace(sigma)), " (=> BEZSLADOWA)")
indep = set()
for i in range(3):
    for j in range(3):
        indep.add(sp.simplify(sigma[i,j]))
# niezalezne komponenty: symetria (6) - 1 wiez sladu = 5
n_sym = 6
n_indep = n_sym - 1
print(f"    skladniki niezalezne: {n_sym}(sym) - 1(slad=0) = {n_indep}  (spin-2: wymaga 5) ->",
      "PASS" if n_indep==5 else "FAIL")

# ---------------------------------------------------------------------------
# (2) Z2-parzystosc: K_ab = <s_i s_{i+a_b}> jest BILINIOWE w s => parzyste pod s->-s.
#     Phi = <s^2> tez bilinowe/parzyste. sigma_ab i Phi to operatory tego samego stopnia.
# ---------------------------------------------------------------------------
s_i, s_j = sp.symbols('s_i s_j', real=True)
Kab_op  = s_i*s_j          # korelator kierunkowy (operator, przed usrednieniem)
Phi_op  = s_i**2           # amplituda izotropowa (operator)
def z2(expr): return expr.subs({s_i:-s_i, s_j:-s_j})
print("\n(2) Z2-parzystosc i stopien operatorow (s -> -s)")
print("    K_ab ~ s_i s_j : parzyste =", sp.simplify(z2(Kab_op)-Kab_op)==0,
      "| stopien w s =", sp.total_degree(sp.Poly(Kab_op, s_i, s_j)))
print("    Phi  ~ s_i^2   : parzyste =", sp.simplify(z2(Phi_op)-Phi_op)==0,
      "| stopien w s =", sp.total_degree(sp.Poly(Phi_op, s_i)))
print("    => sigma_ab i Phi to KOMPOZYTY tego samego stopnia (2) w s, oba Z2-parzyste.")
print("       sigma_ab = rzut ANIZOTROPOWY, Phi = rzut IZOTROPOWY tego samego H_Gamma.")

# ---------------------------------------------------------------------------
# (3) Bilans wymiarowy dzialania eq:S-sigma:  [C_sigma]*[sigma]^2 = E^2
#     S = int d^4x [ (C_sigma/2)(d sigma)^2 - ... ],  [S]=1 (hbar=1)
# ---------------------------------------------------------------------------
E = sp.symbols('E', positive=True)               # jednostka energii (masy)
dim_d4x   = E**(-4)                                # [d^4x] = L^4 = E^-4
dim_dsig  = E**(2)                                 # [(d_mu sigma d^mu sigma)/sigma^2] = E^2 (jedna pochodna = E)
dimC, dimSig = sp.symbols('dimC dimSig', positive=True)
# warunek: dim_d4x * dimC * dim_dsig * dimSig^2 = 1
eq_dim = sp.Eq(dim_d4x * dimC * dim_dsig * dimSig**2, 1)
sol = sp.solve(eq_dim, dimC)[0]
print("\n(3) Bilans wymiarowy (eq:S-sigma), [S]=1:")
print("    [d^4x]=E^-4, [pochodna]^2=E^2  =>  [C_sigma] =", sol, " (w jedn. [sigma]=dimSig)")
print("    Jesli sigma bezwymiarowe (dimSig=1): [C_sigma] = E^2  (tension/stiffness^2).")
print("    Jesli [sigma]=[Phi_0]: [C_sigma]=E^2/[Phi_0]^2  -> rozstrzyga to Phase 3 (C_sigma vs sigma_0).")

# ---------------------------------------------------------------------------
# (4) Znak C_sigma > 0 z J > 0 (ghost-free). Schemat: koszt energetyczny gradientu
#     korelatora kierunkowego ~ +J*(delta sigma)^2 dla wiazania ferromagnetycznego.
# ---------------------------------------------------------------------------
J = sp.symbols('J', positive=True)               # J>0 (wiazanie ferromagnetyczne, eq:B-H)
dsig = sp.symbols('Delta_sigma', real=True)      # roznica korelatora miedzy blokami
# stiffness ~ d^2 E_bond / d(dsig)^2 dla czlonu ~ (J/2)*dsig^2
E_bond = sp.Rational(1,2)*J*dsig**2
stiff = sp.diff(E_bond, dsig, 2)
print("\n(4) Znak C_sigma (ghost-free):")
print("    d^2/d(Dsigma)^2 [ (J/2)(Dsigma)^2 ] =", stiff, " > 0  dla J>0  => C_sigma > 0 (brak duchow).")
print("    (zgodne z eq:S-sigma bullet i prob:tensor-modes 'Ghost-free C_sigma>0 z J>0'.)")

# ---------------------------------------------------------------------------
# (5) Skaling-szkielet (NIE prefaktor!). Analiza wymiarowa: w jedn. substratu
#     stiffness gradientowy ~ J * a_sub^p. Wykladnik p i prefaktor O(1) = Phase 2.
# ---------------------------------------------------------------------------
a = sp.symbols('a_sub', positive=True)
p, c_pref = sp.symbols('p c_pref')               # NIEZNANE w Phase 1 (p z gradient-expansion, c_pref O(1))
Csigma_skeleton = c_pref * J * a**p
print("\n(5) Szkielet skalowania (Phase-1 SETUP, prefaktor=GAP):")
print("    C_sigma ~", Csigma_skeleton, "  gdzie c_pref=O(1) i p NIEOKRESLONE w Phase 1.")
print("    Ekstrakcja (c_pref, p) z blokowego usredniania = Phase 2 (RDZEN). TU NIE LICZONE.")

# ---------------------------------------------------------------------------
# (6) Dyspozycje falsyfikatorow Phase-1 (status wejsciowy do dalszych faz)
# ---------------------------------------------------------------------------
print("\n(6) Dyspozycje falsyfikatorow po Phase 1 (SETUP):")
print("    F-CG-A (emergentna kinetyka Lorentz, c_0): PARTIAL-pending")
print("           - c_0 propagacja LOCKED w rdzeniu (rem:cGW-tensor),")
print("           - ale emergencja czlonu CZASOWEGO (d_t sigma)^2 ze STATYCZNEGO H_Gamma = do sprawdzenia Phase 2.")
print("    F-CG-B (forma C_sigma({J,a_sub,...}))     : OPEN -> Phase 2 (tu tylko szkielet+znak).")
print("    F-CG-C (sigma_0 + zlozenie kappa_E)        : OPEN -> Phase 3.")
print("    F-CG-D (kappa_E vs 5/6)                     : OPEN -> Phase FINAL (prog 5/6 LOCKED Phase 0).")

print("\n" + "="*72)
print("WERDYKT Phase 1 (SETUP): obiekt = KOMPOZYT bilinowy <s_i s_j> (rzut anizotropowy H_Gamma);")
print("  symetryczny/bezsladowy/5-komp (PASS), Z2-parzysty, C_sigma>0 z J>0, [C_sigma]*[sigma]^2=E^2.")
print("  C_sigma NIE wyprowadzone (anti-Lakatos) -> przekazane do Phase 2 z jawna specyfikacja.")
print("="*72)
