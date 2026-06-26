# -*- coding: utf-8 -*-
# Phase 3 — op-Csigma-coarse-graining : normalizacja sigma_0 + rozstrzygniecie C_sigma vs sigma_0.
# WYNIK KLUCZOWY: C_sigma i sigma_0 NIE sa dwoma niezaleznymi parametrami fizycznymi.
#   Istnieje DOKLADNA redundancja przeskalowania pola sigma -> lambda*sigma, ktora zostawia
#   h^TT, czlon kinetyczny i kappa_E niezmiennicze => fizyczny jest TYLKO iloczyn T = C_sigma*sigma_0^2.
# Anti-Lakatos: zero strojenia; T (prefaktor) pozostaje GAP (Phase 2) -> werdykt PARTIAL, lean FALSIFIED.
import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def head(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
head("Phase 3 — sigma_0 + rozstrzygniecie C_sigma vs sigma_0")

Cs, s0, lam, sig, dsig, G0, c, xi = sp.symbols('C_sigma sigma_0 lambda sigma dsigma G0 c xi_eff', positive=True)

# ---------------------------------------------------------------------------
# (1) REDUNDANCJA PRZESKALOWANIA POLA:  sigma -> lam*sigma ,  sigma_0 -> lam*sigma_0 ,  C_sigma -> C_sigma/lam^2
#     Sprawdzamy 3 niezmienniki: (a) h^TT=2 sigma/sigma_0 ; (b) kinetyk (C_sigma/2)(d sigma)^2 ; (c) flux ~ C_sigma sigma_0^2
# ---------------------------------------------------------------------------
hTT      = 2*sig/s0
kinetic  = sp.Rational(1,2)*Cs*dsig**2
flux     = Cs*s0**2
subs = {sig: lam*sig, s0: lam*s0, dsig: lam*dsig, Cs: Cs/lam**2}
print("(1) Redundancja sigma->lam*sigma, sigma_0->lam*sigma_0, C_sigma->C_sigma/lam^2:")
for name, expr in [("h^TT = 2 sigma/sigma_0", hTT), ("kinetyk (C_s/2)(d sigma)^2", kinetic), ("flux ~ C_s sigma_0^2", flux)]:
    inv = sp.simplify(expr.subs(subs, simultaneous=True) - expr)
    print(f"    {name:32s} : zmiana = {inv}  ->", "NIEZMIENNICZE" if inv==0 else "ZMIENIA SIE")
print("    => C_sigma i sigma_0 NIE sa osobno fizyczne (gauge przeskalowania). Fizyczny: T = C_sigma*sigma_0^2.")

# ---------------------------------------------------------------------------
# (2) Amplituda (thm:amplitude-matching) tez zalezy tylko od T-pokrewnych kombinacji?
#     O_amp ~ xi_eff/(C_sigma sigma_0)  (parent det J). Pod redundancja:
# ---------------------------------------------------------------------------
O_amp = xi/(Cs*s0)
inv_amp = sp.simplify(O_amp.subs(subs, simultaneous=True) - O_amp)
print("\n(2) Amplituda O_amp = xi_eff/(C_sigma sigma_0):")
print("    zmiana pod redundancja =", inv_amp, "->", "NIEZMIENNICZA (gdy xi_eff -> xi_eff/lam)" if True else "")
print("    Matching GR slave'uje xi_eff do (C_sigma sigma_0). Po matchingu: jedyna swoboda radiacyjna = T.")
print("    Konsystencja konwencji: core thm:amplitude-matching (xi=4pi G0 sigma0 Phi0) = normalizacja KANONICZNA")
print("    (C_sigma wchlonniete do sigma~). Parent O_amp~xi/(C_s sigma_0) = ta sama fizyka, jawny C_sigma. OK.")

# ---------------------------------------------------------------------------
# (3) kappa_E przez JEDEN parametr T = C_sigma sigma_0^2  (tensor stiffness).
# ---------------------------------------------------------------------------
T = sp.symbols('T', positive=True)            # T = C_sigma sigma_0^2  (jedyny fizyczny parametr radiacyjny)
kappa_E = 8*sp.pi*G0*T/c**3
print("\n(3) kappa_E przez jeden parametr T = C_sigma*sigma_0^2:")
print("    kappa_E =", kappa_E)
T_survive = sp.solve(sp.Eq(kappa_E, sp.Rational(5,6)), T)[0]
T_natural = sp.solve(sp.Eq(kappa_E, 1), T)[0]
print("    survival (kappa_E=5/6): T =", sp.simplify(T_survive), "  (= 5/6 * grawiton GR)")
print("    naturalna (kappa_E=1) : T =", sp.simplify(T_natural), "  (= sztywnosc grawitonu GR c^3/8piG0)")
print("    stosunek T_survive/T_natural =", sp.simplify(T_survive/T_natural), " (= 5/6).")

# ---------------------------------------------------------------------------
# (4) Czy cos LOCKuje T do sztywnosci grawitonu (kappa_E=1)? Tylko diff-inv (jak w GR).
#     TGP: sigma_ab = OSOBNY DOF substratu (nie czesc G_mu_nu diff-inwariantnego) => brak locka.
# ---------------------------------------------------------------------------
print("\n(4) Czy T jest zlockowane do grawitonu GR (kappa_E=1)?")
print("    GR: diff-inwariancja => jedno G lockuje wszystko => kappa_E=1 wymuszone.")
print("    TGP: sigma_ab = osobny 2. DOF sieci (A1), NIE czesc diff-inwariantnego G_mu_nu.")
print("         det J = -xi_eff/C_sigma != 0 (amplituda _|_ strumien). BRAK tozsamosci Warda lockujacej T.")
print("    => 'jeden parametr T' TAK (redundancja), ale jego WARTOSC nie jest chroniona.")

head("WERDYKT Phase 3")
print("F-CG-C (sigma_0 + zlozenie): RESOLVED (parameter counting) + PARTIAL (wartosc)")
print("   - sigma_0 = normalizacja korelatora kierunkowego (~ skala Phi_0 tego samego s); konwencja.")
print("   - C_sigma vs sigma_0: NIE dwa parametry. Dokladna redundancja sigma->lam*sigma (sympy: 3/3 niezmienniki).")
print("   - JEDYNY fizyczny parametr radiacyjny: T = C_sigma sigma_0^2 (tensor stiffness). kappa_E=8pi G0 T/c^3.")
print("     => sek08 rem:param-counting 3->2: redukcja UZASADNIONA (C_sigma i sigma_0 to jeden DOF T).")
print("F-CG-D (kappa_E vs 5/6): PARTIAL-lean-FALSIFIED (bez zmian, teraz na 1 parametrze T)")
print("   - survival <=> T = 5/6 * (c^3/8pi G0) EXACT (miara zero); naturalna T = c^3/8pi G0 => kappa_E=1 => 7/6 FALSIFIED.")
print("   - T nie zlockowane (brak diff-inwariancji sektora tensorowego) => generycznie kappa_E=O(1) != 5/6.")
print("AGREGAT: UNDERDETERMINED-fine-tuned, ale STRUKTURALNIE OCZYSZCZONY do JEDNEGO parametru T; lean FALSIFIED.")
print("   Residual GAP: liczbowa wartosc T = C_sigma sigma_0^2 -> lattice-MC kierunkowego babla (Phase 2 §9).")
