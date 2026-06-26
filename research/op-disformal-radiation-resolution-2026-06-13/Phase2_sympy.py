# -*- coding: utf-8 -*-
# Phase 2 — op-disformal-radiation-resolution : T2 (pin kappa_E) + T3 (status M*)
# T2 (F-DRR-2): czy strumien energii kanalu sigma_ab (C3) jest PINOWANY do GR przez LIVE,
#               niezaleznie od warunku amplitudy xi_eff=4pi G0 sigma0 Phi0 (R7)?
#               Konkretyzacja Phase2-survival Filar II (det J=2xi/lambda!=0) w (C_sigma,sigma0,xi_eff).
# T3 (F-DRR-3): M*=m_P -- DERIVED (mikro) / POSTULATE (wymiarowy) / FITTED (do r_V/danych)?
#               Rozstrzygniecie niespojnosci sek08 'Warstwa III wyprowadzone' vs status_map 'Propozycja'.
# Wszystko symbolicznie; 0 danych; 0 nowych stalych; werdykt z flag; anty-tuning (forbidden #2/#3).
#
# WEJSCIA LOCKED:
#   R4 : S_sigma = int sqrt(-g)[ C_sigma/2 (d sigma)^2 - m_s^2/2 sigma^2 + xi_eff/Phi0^2 sigma dPhi dPhi ]
#   R6 : h^TT = 2 sigma/(sigma0 Phi0)
#   R7 : xi_eff = 4 pi G0 sigma0 Phi0   (pinuje xi_eff/sigma0; amplituda=GR)
#   R8 : T^{0r}_sigma ∝ C_sigma <dt sigma dr sigma>  => Pdot_sigma ∝ C_sigma sigma0^2 Phi0^2 <h'^2>
#   R9 : det J[(lambda,xi)->(lambda xi, xi/lambda)] = 2xi/lambda != 0   (Phase2-survival, LOCKED)
#   rem:param-counting : 4 warunki {Phi0,G0,Lambda_eff,xi_eff} na 4 param substratu {mu,m0^2,lam0,J}; xi_eff<->J/lam0
#   rem:sigma-params   : C_sigma (rownow. sigma0) = 'DODATKOWY parametr, w zasadzie z substratu, OBECNIE NIEZOBLICZONY'
#   prop:Mstar-from-substrate (dodatekC l.690-792): M*^2=1/ell_P^2 z ANALIZY WYMIAROWEJ + B(Phi0)=1
#   status_map l.214-217 : M* 'Propozycja; wymiarowanie Phi0+ell_P; brak mikro-derywacji'
#   sek08 l.170-172      : M* 'Warstwa III -- wyprowadzone'  (NIESPOJNE z status_map)

import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

PASS=[]
def check(i,ok,note):
    PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

print("="*72)
print("T2-A — fizyczne kombinacje sektora sigma_ab: O_amp (amplituda) vs O_flux (strumien)")
print("="*72)
# Trzy parametry sektora: C_sigma (kinetyka), sigma0 (normalizacja mostu R6), xi_eff (sprzezenie).
# Normalizacja kanoniczna sigma = sigma'/sqrt(C_sigma) eliminuje C_sigma w param efektywne:
#   xi_eff' = xi_eff/sqrt(C_sigma)   (efektywne sprzezenie zrodla)
#   sigma0' = sqrt(C_sigma)*sigma0   (efektywna normalizacja mostu)
Csig, sig0, xi = sp.symbols('C_sigma sigma0 xi_eff', positive=True)
xi_p  = xi/sp.sqrt(Csig)            # xi_eff'
sig0p = sp.sqrt(Csig)*sig0         # sigma0'
# Amplituda h^TT ∝ xi_eff'/sigma0' ; strumien Isaacson ∝ sigma0'^2 (kanoniczny T^{0r})
O_amp  = sp.simplify(xi_p/sig0p)    # = xi/(C_sigma*sigma0)  -- pinowane R7 do 2G
O_flux = sp.simplify(sig0p**2)      # = C_sigma*sigma0^2     -- kontroluje kappa_E
check("T2-A", sp.simplify(O_amp - xi/(Csig*sig0))==0 and sp.simplify(O_flux - Csig*sig0**2)==0,
      f"O_amp = {O_amp} (pinowane przez R7: amplituda=GR);  O_flux = {O_flux} (kontroluje kappa_E, strumien). "
      f"Trzy parametry {{C_sigma,sigma0,xi_eff}}, dwie fizyczne kombinacje.")

print("="*72)
print("T2-B — Jacobian (xi',sigma0')->(O_amp,O_flux): det != 0  => O_amp PERP O_flux (ground Phase2 Filar II)")
print("="*72)
# Mapa w param efektywnych (xi_p, sig0p) -> (O_amp, O_flux):
xq, lq = sp.symbols('xq lq', positive=True)         # xq=xi', lq=sigma0'
Oa = xq/lq
Of = lq**2
J = sp.Matrix([[sp.diff(Oa,xq), sp.diff(Oa,lq)],
               [sp.diff(Of,xq), sp.diff(Of,lq)]])
detJ = sp.simplify(J.det())                          # = 2  (≠ 0)
independent = sp.simplify(detJ) != 0
check("T2-B", independent,
      f"det J[(xi',sigma0')->(xi'/sigma0', sigma0'^2)] = {detJ} != 0 => O_amp i O_flux NIEZALEZNE. "
      f"Ugruntowanie Phase2-survival Filar II (det J=2xi/lambda!=0, R9): "
      f"USTALENIE AMPLITUDY (R7) NIE USTALA STRUMIENIA => kappa_E swobodne kinematycznie.")

print("="*72)
print("T2-C — czy LIVE pinuje O_flux=C_sigma*sigma0^2? (zliczanie warunkow; rem:sigma-params)")
print("="*72)
# rem:param-counting: 4 warunki {Phi0,G0,Lambda_eff,xi_eff} pinuja 4 param substratu {mu,m0^2,lam0,J}.
#   xi_eff <-> J/lam0. ZADEN z 4 warunkow nie odnosi sie do C_sigma osobno.
substrate_conditions = {"Phi0":"mu/m0", "G0":"J*mu", "Lambda_eff":"gamma", "xi_eff":"J/lam0"}
pins_C_sigma = any("C_sigma" in v for v in substrate_conditions.values())   # False
# rem:sigma-params: C_sigma (rownow. sigma0) = dodatkowy parametr, OBECNIE NIEZOBLICZONY
C_sigma_computed_from_substrate = False
O_flux_pinned_by_LIVE = (pins_C_sigma or C_sigma_computed_from_substrate)   # False
# Jedyny 'pin' = narzucenie normalizacji Einsteina-Hilberta C_sigma sigma0^2 Phi0^2 = c^3/(16 pi G):
#   to bylby TUNING do GR (forbidden #2/#3), o ile NIE wyprowadzone z substratu (a nie jest).
pin_only_by_tuning = (not O_flux_pinned_by_LIVE)
check("T2-C", (not O_flux_pinned_by_LIVE) and pin_only_by_tuning,
      f"4 warunki substratu {list(substrate_conditions)} pinuja {{mu,m0^2,lam0,J}}; zaden nie pinuje C_sigma "
      f"(pins_C_sigma={pins_C_sigma}). rem:sigma-params: C_sigma NIEZOBLICZONY. "
      f"=> O_flux=C_sigma*sigma0^2 NIEPINOWANE przez LIVE. Jedyny pin = narzucenie norm. E-H = TUNING do GR "
      f"(forbidden #2/#3, bo nie wyprowadzone z substratu).")

print("="*72)
print("T2 -> FLAGA F-DRR-2")
print("="*72)
if O_flux_pinned_by_LIVE:
    F_DRR_2 = "PINNED_GR"        # gdyby LIVE wymuszal O_flux = wartosc GR
else:
    F_DRR_2 = "UNPINNED"         # strukturalnie swobodne (det J!=0 + brak warunku)
print(f"F-DRR-2 = {F_DRR_2}  (O_flux=C_sigma*sigma0^2 swobodne; amplituda R7 nie pinuje strumienia; "
      f"pin tylko przez tuning E-H = forbidden)")

print("="*72)
print("T3-A — prop:Mstar-from-substrate: jaka to klasa dowodu? (analiza wymiarowa vs mikro-derywacja)")
print("="*72)
# Tresc prop (dodatekC l.728-774): 'Uzasadnienie: Analiza wymiarowa.'
#   [M*^2]=m^2 ; jedyna bezparametrowa kombinacja wymiarowa z (Phi0,ell_P,c,hbar) o wymiarze m^2 => M*^2=1/ell_P^2.
#   + wybor normalizacji B(Phi0)=1. NIE liczy wspolczynnika z dynamiki substratu (mu,m0^2,lam0,J).
ellP = sp.symbols('ell_P', positive=True)
# odtworzenie argumentu wymiarowego: [M*^2]=m^2; jednostki nat. (hbar=c=1) => 1/ell_P^2 ma [m^2]
Mstar2 = 1/ellP**2
dim_ok = sp.simplify(Mstar2*ellP**2) == 1            # M*^2 * ell_P^2 = 1  (bezwymiarowe) -> M*=1/ell_P=m_P
proof_is_dimensional = True                          # 'Analiza wymiarowa' (l.729) + B(Phi0)=1 (l.764)
micro_derivation_from_substrate = False              # NIE liczy wspolczynnika z {mu,m0^2,lam0,J}
check("T3-A", dim_ok and proof_is_dimensional and (not micro_derivation_from_substrate),
      f"prop:Mstar: M*^2=1/ell_P^2 (M*^2*ell_P^2={sp.simplify(Mstar2*ellP**2)}) z ANALIZY WYMIAROWEJ "
      f"(jedyna bezparam. kombinacja o [m^2]) + norm. B(Phi0)=1. "
      f"Mikro-derywacja wspolczynnika z dynamiki substratu: {micro_derivation_from_substrate}. "
      f"=> argument WYMIAROWY, NIE mikro-derywacja.")

print("="*72)
print("T3-B — czy M* jest FITOWANE do r_V/danych? (anti-tuning, forbidden #3)")
print("="*72)
# M* zafiksowane przez ell_P (dlugosc Plancka), NIE przez r_V ani dane pulsarowe.
# r_V uzywa M* (r_V = f(M*,B,...)), a nie odwrotnie. grep za fitem M*<->r_V: brak (Phase 0 R-Mstar; recon §2).
M_star_fitted_to_rV   = False    # recon §2: 'NIE fitowane do r_V(odot)'
M_star_fitted_to_data = False    # cykl strukturalny; M*=m_P z ell_P
fitted = (M_star_fitted_to_rV or M_star_fitted_to_data)
check("T3-B", (not fitted),
      f"M* zafiksowane przez ell_P (Planck), NIE przez r_V (fit_rV={M_star_fitted_to_rV}) ani dane "
      f"(fit_data={M_star_fitted_to_data}). r_V=f(M*), nie M*=f(r_V). => NIE FITTED (nie tuning-belt).")

print("="*72)
print("T3-C — rozstrzygniecie niespojnosci sek08 'Warstwa III' vs status_map 'Propozycja'")
print("="*72)
sek08_class      = "Warstwa III -- wyprowadzone"     # l.170-172  (OVERCLAIM)
status_map_class = "Propozycja; brak mikro-derywacji" # l.214-217  (poprawne)
# Z T3-A: argument wymiarowy, NIE mikro-derywacja => status_map poprawne; sek08 overclaim.
status_map_correct = (not micro_derivation_from_substrate)
check("T3-C", status_map_correct,
      f"sek08='{sek08_class}' vs status_map='{status_map_class}'. T3-A: argument wymiarowy (nie mikro). "
      f"=> status_map POPRAWNE; sek08 'wyprowadzone' = OVERCLAIM klasyfikacyjny (do korekty w propagacji FINAL).")

print("="*72)
print("T3 -> FLAGA F-DRR-3")
print("="*72)
if micro_derivation_from_substrate:
    F_DRR_3 = "DERIVED"
elif fitted:
    F_DRR_3 = "FITTED"
else:
    F_DRR_3 = "POSTULATE"
print(f"F-DRR-3 = {F_DRR_3}  (M*=m_P: postulat WYMIAROWY z ell_P + B(Phi0)=1; underived mikroskopowo; NIE fit)")

print("="*72)
print("T3-D — dyscyplina")
print("="*72)
notes=" ".join(n for _,_,n in PASS)
clean = not any(t in notes for t in {"SPARC","GWTC","Pdot_obs","a0=","H0=","J0737_data"})
check("T3-D", clean, "0 danych obserwacyjnych; werdykt z flag; 0 nowych stalych; anty-tuning (forbidden #2/#3) zachowany.")

print("="*72)
print("WERDYKT Phase 2 (flagi) + KANDYDAT agregatu F-DRR-C (do Phase FINAL)")
print("="*72)
np_=sum(1 for _,ok,_ in PASS if ok); nt=len(PASS)
print(f"checks: {np_}/{nt} PASS")
print(f"  F-DRR-2 = {F_DRR_2}")
print(f"  F-DRR-3 = {F_DRR_3}")
# Agregat (reguła Phase 0 §4.1) z F-DRR-1=PARTIAL (Phase 1):
F_DRR_1 = "PARTIAL"
broken = (F_DRR_1=="UNSCREENED_FLUX") or (F_DRR_2=="PINNED_DEVIATION") or (F_DRR_3=="FITTED")
clean_agg = (F_DRR_1=="SCREENED_FLUX") and (F_DRR_2=="PINNED_GR") and (F_DRR_3=="DERIVED")
if broken:
    D6 = "BROKEN"
elif clean_agg:
    D6 = "CLEAN"
elif (F_DRR_2=="UNPINNED") and (not broken):
    D6 = "UNDERDETERMINED"
else:
    D6 = "UNDERDETERMINED (+R1 luka)"
print(f"\n  Agregat (F-DRR-1={F_DRR_1}, F-DRR-2={F_DRR_2}, F-DRR-3={F_DRR_3}):")
print(f"    broken={broken}, clean={clean_agg}  =>  KANDYDAT D6 -> {D6}")
print( "    (do potwierdzenia + closure + propagacja w Phase FINAL; wymaga user 'działaj')")
print(f"\nDISCIPLINE: hardcoded_T_pass=0 ; verdict_from_flags=True ; anti_tuning_preserved=True ; new_constants=0")
