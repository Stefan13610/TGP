# -*- coding: utf-8 -*-
# Phase 2 — op-Csigma-coarse-graining : RDZEN. Ekstrakcja C_sigma przez propagator KOMPOZYTU
# (sigma_ab = <s_i s_j>_TF jest operatorem zlozonym -> jego kinetyka = babel/composite propagator).
# Anti-Lakatos: liczymy STRUKTURE + SKALING + znak w kontrolowanym przyblizeniu Gaussa;
#               prefaktor O(1) zalezny od schematu/normalizacji sigma_0 -> JAWNIE NIE pinowany do 5/6.
import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def head(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

head("Phase 2 — RDZEN: C_sigma z coarse-grainingu korelatora kierunkowego")

# ---------------------------------------------------------------------------
# (A) METODA: kinetyka kompozytu sigma_ab ~ <s s>_TF = inverse babla 2-punktowego.
#     W 3D substracie (dodatekQ: d=3, n=1, Z2) babel z propagatorem G(q)=1/(q^2+m^2):
#     Pi(p) = \int d^3q/(2pi)^3  1/[(q^2+m^2)((p+q)^2+m^2)]
#     Weryfikujemy rozwiniecie malego p PRZEZ calki beta-funkcyjne (bez zgadywania wzoru).
# ---------------------------------------------------------------------------
q, m = sp.symbols('q m', positive=True)
def radial(num_pow, den_pow):   # (1/2pi^2) \int_0^inf q^2 * q^num_pow / (q^2+m^2)^den_pow dq
    a = sp.Rational(2+num_pow+1, 2)        # 2a-1 = 2+num_pow
    integ = sp.integrate(q**(2+num_pow)/(q**2+m**2)**den_pow, (q, 0, sp.oo))
    return sp.simplify(integ/(2*sp.pi**2))

Pi0   = radial(0, 2)                         # Pi(0) = \int q^2/D^2
# p^2-coeff: angular avg <(p.q)^2>=p^2 q^2/3 => coeff = (1/2pi^2)\int q^2[ -1/D^3 + (4/3)q^2/D^4 ]
c_p2  = sp.simplify(-radial(0,3) + sp.Rational(4,3)*radial(2,4))
print("(A) Babel 3D (zweryfikowany calkami beta-funkcyjnymi):")
print("    Pi(0)          =", Pi0, "   (oczek. 1/(8*pi*m))  ->", "PASS" if sp.simplify(Pi0-1/(8*sp.pi*m))==0 else "FAIL")
print("    coeff p^2      =", c_p2, "   (oczek. -1/(96*pi*m^3)) ->", "PASS" if sp.simplify(c_p2+1/(96*sp.pi*m**3))==0 else "FAIL")
ratio = sp.simplify(-c_p2/Pi0)
print("    -coeff/Pi(0)   =", ratio, "  => sztywnosc/masa ~ 1/(12 m^2)  ->", "PASS" if sp.simplify(ratio-1/(12*m**2))==0 else "FAIL")
print("    WNIOSEK: babel jest SKONCZONY i ANALITYCZNY w p=0 => emergentny czlon p^2 (sztywnosc) ISTNIEJE,")
print("             dodatni, z C_sigma ~ (sprzezenie^2) * 1/(96*pi*m^3). To realizuje F-CG-A (czesc przestrzenna).")

# ---------------------------------------------------------------------------
# (B) ZNAK i SKALING C_sigma (z J>0). C_sigma > 0 (ghost-free, zgodne z Phase 1 i rdzeniem).
#     Skaling wymiarowy: C_sigma ~ c_pref * J * a_sub^p / (m a_sub)^k. Wykladniki/prefaktor = residual.
# ---------------------------------------------------------------------------
J, a, xi = sp.symbols('J a_sub xi_corr', positive=True)  # xi_corr = 1/m (dlugosc korelacji substratu)
c_pref, p, kk = sp.symbols('c_pref p k', real=True)
Csigma_scaling = c_pref * J * a**p * (xi/a)**kk
print("\n(B) Skaling C_sigma (kontrolowany, prefaktor=residual):")
print("    C_sigma ~", Csigma_scaling, "  ; znak: C_sigma>0 <= J>0 (babel dodatni). c_pref,p,k = GAP.")
print("    (k>0: bliskosc krytycznosci m->0 WZMACNIA sztywnosc ~ xi^3 z babla 1/m^3 -- struktura, nie liczba.)")

# ---------------------------------------------------------------------------
# (C) ZLOZENIE kappa_E (F-CG-C/D). h^TT=2 sigma/sigma_0 => podstawienie sigma=sigma_0 h/2.
#     Sztywnosc tensorowa TGP: (C_sigma sigma_0^2/8)(d h)^2 ; graviton GR: (c^3/(64 pi G0))(d h)^2.
#     kappa_E = (sztywnosc TGP)/(sztywnosc GR) = 8 pi G0 C_sigma sigma_0^2 / c^3.
# ---------------------------------------------------------------------------
Cs, s0, G0, c, h = sp.symbols('C_sigma sigma_0 G0 c h', positive=True)
stiff_TGP = sp.Rational(1,8)*Cs*s0**2          # wspolczynnik (d h)^2 z S_sigma po podstawieniu sigma=sigma0 h/2
stiff_GR  = c**3/(64*sp.pi*G0)                  # wspolczynnik (d h)^2 grawitonu GR
kappa_E   = sp.simplify(stiff_TGP/stiff_GR)
print("\n(C) Zlozenie kappa_E = (sztywnosc tensorowa TGP)/(sztywnosc grawitonu GR):")
print("    kappa_E =", kappa_E, " = 8*pi*G0*C_sigma*sigma_0^2 / c^3")

# ---------------------------------------------------------------------------
# (D) TEST value-blind vs 5/6 (prog LOCKED Phase 0) i vs 1 (naturalna). NIE stroimy.
# ---------------------------------------------------------------------------
print("\n(D) Test progu (value-blind):")
print("    survival  <=> kappa_E = 5/6 EXACT (miara zero, parent).")
print("    naturalna  =>  kappa_E = 1  => total 1+1/6 = 7/6 = galaz B PR-025 (FALSIFIED).")
print("    Pytanie rozstrzygajace: czy 8 pi G0 C_sigma sigma_0^2 / c^3 = 5/6 ?")
# G0 z sektora skalarnego (rem:param-counting: G0 ~ J*mu), C_sigma z sektora tensorowego (J*a_sub^p) -- NIEZALEZNE projekcje.
mu = sp.symbols('mu', positive=True)
print("    G0 ~ J*mu (skalar, <s^2>);  C_sigma ~ J*a_sub^p (tensor, <s s>).  RÓŻNE projekcje TEGO SAMEGO s.")
print("    => kappa_E ~ (J mu)(J a_sub^p) sigma_0^2 / c^3 : iloczyn NIEZALEZNYCH prefaktorow substratu.")

# ---------------------------------------------------------------------------
# (E) BRAK SYMETRII LOCKUJACEJ kappa_E. W GR jedno G lockuje amplitude<->strumien (1 parametr).
#     W TGP: amplituda zlockowana (xi_eff=4 pi G0 sigma0 Phi0, thm:amplitude-matching),
#     ale strumien ~ C_sigma sigma0^2 NIEZALEZNY (det J=-xi/C_sigma != 0, parent). 2 parametry.
# ---------------------------------------------------------------------------
xi_eff, Phi0 = sp.symbols('xi_eff Phi0', positive=True)
O_amp  = xi_eff/(Cs*s0)     # ~ amplituda (pinned do GR)
O_flux = Cs*s0**2           # ~ strumien (= kappa_E * GR)
Jac = sp.Matrix([[sp.diff(O_amp,Cs), sp.diff(O_amp,s0)],
                 [sp.diff(O_flux,Cs), sp.diff(O_flux,s0)]])
detJ = sp.simplify(Jac.det())
print("\n(E) Czy jakas symetria lockuje kappa_E do 5/6 (lub 1)?")
print("    det J(amp,flux | C_sigma,sigma_0) =", detJ, "!= 0  => amplituda _|_ strumien (2 param TGP vs 1 GR).")
print("    BRAK toznosci Warda lockujacej C_sigma sigma_0^2 do sztywnosci grawitonu => kappa_E NIE jest chronione.")
print("    => generycznie kappa_E = O(1), ale ANI 1 ANI 5/6 nie sa wymuszone. 5/6 pozostaje MIARA ZERO.")

head("WERDYKT Phase 2")
print("F-CG-A (Lorentz/kinetyka): PASS-CONDITIONAL")
print("   - babel daje DODATNI, skonczony czlon p^2 => sztywnosc PRZESTRZENNA emerguje z H_Gamma (C_sigma>0).")
print("   - czlon CZASOWY (d_t sigma)^2 i c_0 = dziedziczone z metryki emergentnej (rem:cGW-tensor LOCKED),")
print("     jak w sektorze skalarnym. Lorentz NIE wyprowadzony de novo z H_Gamma -> flaga zaleznosci.")
print("F-CG-B (forma C_sigma): PARTIAL")
print("   - metoda (composite babel) + skaling C_sigma~c_pref*J*a_sub^p/(m a_sub)^k + znak C_sigma>0: DERIVED.")
print("   - prefaktor O(1) i wykladniki: schemat-/normalizacja-zalezne -> NIE pinowane (anti-Lakatos).")
print("F-CG-C/D (kappa_E vs 5/6): PARTIAL-lean-FALSIFIED")
print("   - kappa_E = 8 pi G0 C_sigma sigma_0^2 / c^3 ZLOZONE (struktura jawna).")
print("   - BRAK symetrii lockujacej (det J != 0) => kappa_E generyczne O(1), 5/6 = miara zero (NIE chronione).")
print("   - skaling nie odroznia 5/6 od 1 bez prefaktora => survival ani potwierdzone ani wykluczone liczbowo,")
print("     ale 'naturalna' O(1) != 5/6 => lean FALSIFIED (wzmocniony vs parent: kappa_E teraz O(1)-bounded).")
print("AGREGAT F-CG-E: UNDERDETERMINED-fine-tuned (WEZSZY) ; lean FALSIFIED.")
print("   Residual GAP: prefaktor O(1) C_sigma -> wymaga lattice-MC kierunkowego babla + Phase 3 (sigma_0).")
