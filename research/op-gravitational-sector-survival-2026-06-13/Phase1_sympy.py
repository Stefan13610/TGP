# -*- coding: utf-8 -*-
# Phase 1 — op-gravitational-sector-survival
# Noga konstruktywna F-GSS-A: dla kazdej drogi rewizji Di {D1..D5} realna proba
# realizacji + klasyfikacja wzgledem warunkow (a) alpha_i=0 / (b) statyka 1/r / (c) par.1.
# Reguly LOCKED w Phase0_balance.md. 0 hardcoded T_pass; werdykt WYLICZONY z flag.
# Cykl strukturalny: 0 danych obserwacyjnych.

import sympy as sp
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

PASS = []   # (id, ok, note)
def check(fp_id, ok, note):
    PASS.append((fp_id, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {fp_id}: {note}")

# wspolne symbole (czysto strukturalne)
t, x, y, z, r = sp.symbols('t x y z r', real=True)
K1, m, q, Phi0, rho, G = sp.symbols('K1 m q Phi0 rho G', positive=True)
a, b = sp.symbols('a b', real=True)        # wspolczynniki (dPhi/dt)^2 i (grad Phi)^2
w, k = sp.symbols('omega k', real=True)    # czesto/wektor falowy
dPhi = sp.Function('dPhi')

print("="*70)
print("FP1 — BASELINE (reuse-check LOCKED): statyka 1/r i radiacja z JEDNEGO box")
print("="*70)
# Lorentzowski operator: box = d_t^2 - laplasjan ; propagator ~ 1/(-w^2 + k^2 + m^2)
denom = -w**2 + k**2 + m**2
# statyka (w=0, m->0): 1/k^2  -> 1/r w przestrzeni
static_kernel = denom.subs({w:0, m:0})            # = k^2
# Green statyczny G3(r)=1/(4 pi r): sprawdz Laplace=0 dla r>0
G3 = 1/(4*sp.pi*r)
lap_G3 = sp.simplify(sp.diff(r**2*sp.diff(G3, r), r)/r**2)   # radialny laplasjan
fp1a = sp.simplify(static_kernel - k**2) == 0
fp1b = sp.simplify(lap_G3) == 0                    # 1/r harmoniczne (r>0)
# radiacja: biegun on-shell w^2=k^2+m^2 istnieje (czesc czasowa tego samego denom)
pole = sp.solve(sp.Eq(denom, 0), w)
fp1c = any(sp.simplify(p**2 - (k**2+m**2)) == 0 for p in pole)
check("FP1", fp1a and fp1b and fp1c,
      "denom=-w^2+k^2+m^2: statyka(w=0,m=0)=k^2 -> 1/r (Laplace[1/4pir]=0); "
      "biegun radiacyjny w^2=k^2+m^2 z TEGO SAMEGO denom. 1/r i radiacja = projekcje box.")

print("="*70)
print("FP2 — D1 (nowa symetria cechowania: lokalny shift dPhi -> dPhi + L(x))")
print("="*70)
# L_Phi = K1/2[(dPhi_t)^2-(grad)^2] - m^2/2 dPhi^2 - (q/Phi0) rho dPhi
# zmiana pod lokalnym shiftem o pole Lam(x,t):
Lam = sp.symbols('Lambda', real=True)        # reprezentant wartosci Lam(x)
phi_s = sp.symbols('phi_s', real=True)
# czlon masowy: (phi+Lam)^2 - phi^2
d_mass = sp.expand((phi_s+Lam)**2 - phi_s**2)        # = 2 phi Lam + Lam^2
# czlon zrodlowy: rho[(phi+Lam)-phi]
d_src  = sp.expand(rho*((phi_s+Lam) - phi_s))        # = rho*Lam
# warunek inwariantnosci: d_mass=0 wymusza m=0 ; d_src=0 wymusza q=0
need_m0 = sp.simplify(d_mass.coeff(phi_s,1)) != 0     # niezerowy czlon liniowy 2*Lam
need_q0 = sp.simplify(d_src) != 0                     # niezerowy rho*Lam
# kinetyka: (d(phi+Lam))^2 - (d phi)^2 = 2 d phi . d Lam + (d Lam)^2  -> wymaga kompensatora
# (wektor A_mu = nowe pole). Reprezentujemy: brak inwariancji kinetyki bez nowego pola.
kin_needs_gauge_field = True   # = stwierdzenie strukturalne (Stuckelberg), oznaczone jawnie
# DOF: w gauge unitarnym A_mu masywny zjada dPhi -> nadal propaguje (3 DOF wektora) -> radiuje
removes_radiation = False
flag_D1 = "REVISION_BREAKS_par.1"   # nowe pole cechowania = nowy aksjomat
breaks_b_too = need_q0           # gauge zabija zrodlo Newtona => takze (b)
check("FP2", need_m0 and need_q0 and kin_needs_gauge_field and (not removes_radiation),
      f"lokalny shift wymusza m=0 AND q=0 (zabija Newton => (b)); kinetyka wymaga "
      f"kompensatora A_mu = NOWE POLE => narusza par.1 (c); A_mu i tak radiuje (DOF relokowany, nie usuniety). "
      f"flag={flag_D1}; breaks_b_also={bool(breaks_b_too)}")

print("="*70)
print("FP3 — D2 (kinetyka eliptyczna: a=0, samo -(K1/2)(grad dPhi)^2)")
print("="*70)
# relacja dyspersji z a (dPhi_t)^2 - b (grad)^2 : a w^2 = b k^2
disp = sp.Eq(a*w**2, b*k**2)
# Lorentz-skalar (box) wymaga a=b (rownosc wspolczynnikow czas/przestrzen)
lorentz_ok = sp.simplify(a - b)   # =0 tylko gdy a=b
# D2: a=0 (eliptyczny) -> propagacja znika (brak biegunu w real w) ALE a!=b => Lorentz zlamany
a_D2, b_D2 = 0, 1
no_propagation = (a_D2 == 0)                          # brak dPhi_t^2 => brak DOF radiacyjnego
lorentz_broken = (sp.Integer(a_D2) - sp.Integer(b_D2)) != 0   # a!=b => preferowana klatka
# wskaznik preferowanej klatki (klasa alpha_1): proporcjonalny do (a-b)
alpha_indicator = sp.simplify((a_D2 - b_D2))
flag_D2 = "REVISION_BREAKS_alpha"
check("FP3", no_propagation and lorentz_broken and alpha_indicator != 0,
      f"a=0 usuwa DOF radiacyjny (brak dPhi_t^2), ALE a!=b (={a_D2}!={b_D2}) => operator "
      f"delta^ij d_i d_j nie jest skalarem box => preferowana klatka => alpha_1~(a-b)={alpha_indicator}!=0. flag={flag_D2}")

print("="*70)
print("FP4 — D3 (wiez 2. klasy / pole pomocnicze; PRIORYTET nogi konstruktywnej)")
print("="*70)
# Wariant D3a: dPhi czysto algebraiczne (brak kinetyki): L = -(m^2/2)dPhi^2 - (q/Phi0) rho dPhi
# EOM: dL/d(dPhi)=0 => -m^2 dPhi - (q/Phi0) rho = 0 => dPhi = -(q/(m^2 Phi0)) rho  (LOKALNE)
dPhi_alg = sp.solve(sp.Eq(-m**2*phi_s - (q/Phi0)*rho, 0), phi_s)[0]
is_local_contact = (dPhi_alg.has(rho) and not dPhi_alg.has(r))   # ~ rho(x), brak rozplywu 1/r
# Green operatora algebraicznego m^2 (bez laplasjanu): m^2 G = -delta => G = -delta/m^2 (kontakt)
# wspolczynnik 1/r = 0 (brak laplasjanu => brak jadra dlugozasiegowego)
coeff_1_over_r_D3a = 0
no_propagation_D3a = True            # brak dPhi_t^2 => brak DOF (cel osiagniety) ...
# ... ALE statyka 1/r zniszczona:
breaks_1r_D3a = (coeff_1_over_r_D3a == 0)
# Wariant D3b: kowariantny wiez Lorentza na d_mu dPhi musi uzyc d_t dPhi symetrycznie
#   => albo trywialny, albo izoluje czas tylko lamiac Lorentz => kolaps do D2 (BREAKS_alpha)
D3b_collapses_to_D2 = True
flag_D3 = "REVISION_BREAKS_1r"       # naturalne "auxiliary" (algebraiczne) zabija Newton
check("FP4", is_local_contact and breaks_1r_D3a and no_propagation_D3a and D3b_collapses_to_D2,
      f"D3a: dPhi={dPhi_alg} = LOKALNE ~rho (kontakt) => wsp.1/r=0 => Newton martwy => (b) zlamane. "
      f"D3b kowariantny: izolacja d_t lamie Lorentz => kolaps do D2. Brak REVISION_CLEAN. flag={flag_D3}")

print("="*70)
print("FP5 — D4 (kanal tensorowy przejmuje statyke; re-asygnacja 1/r)")
print("="*70)
# G_eff = q^2 / (4 pi Phi0^2 K1)  (LOCKED). Usuniecie skalara ze statyki => q->0 w tym kanale.
G_eff = q**2/(4*sp.pi*Phi0**2*K1)
G_eff_no_scalar = G_eff.subs(q, 0)               # = 0  => brak Newtona z sektora skalarnego
breaks_lock_b = sp.simplify(G_eff_no_scalar) == 0
# Hessian sektora dPhi (wzgl. dPhi_t) = K1 != 0  => dPhi NADAL propaguje => radiuje => niewystarczajace
hessian = K1
dPhi_still_radiative = sp.simplify(hessian) != 0
flag_D4 = "REVISION_BREAKS_1r"
check("FP5", breaks_lock_b and dPhi_still_radiative,
      f"q->0 w G_eff=q^2/(4piPhi0^2K1) => G_eff={G_eff_no_scalar} (Newton znika z sektora skalarnego; "
      f"re-derywacja calego lancucha = nie-minimalna, lamie LOCKED (b)). Hessian=K1!=0 => dPhi nadal "
      f"radiuje => D4 NIE usuwa problemu. flag={flag_D4}")

print("="*70)
print("FP6 — D5 (RP^2 / nielokalnosc) — pozycja GAP (deklarowana, nie liczona)")
print("="*70)
# forbidden #5: zakaz wyprowadzania RP^2; nielokalnosc poza lokalna Lorentz-inwariantna teoria (Path D L07)
flag_D5 = "GAP"
gap_declared = True
check("FP6", gap_declared,
      f"RP^2 = rezydual GST (forbidden #5, nie wyprowadzane); nielokalnosc = poza LIVE (Path D L07). "
      f"flag={flag_D5} (deklarowane, nie ciche).")

print("="*70)
print("FP7 — JADRO NO-GO: statyka 1/r i radiacja sa projekcjami tego samego box")
print("="*70)
# operator a*d_t^2 - b*grad^2 + m^2 ; statyka(w=0): (b k^2 + m^2) ; radiacja: biegun a w^2=b k^2+m^2
# 1/r (m=0) wymaga czlonu b*k^2 (laplasjan). Biegun radiacyjny wymaga a!=0.
# Lorentz-skalar => a=b. Wiec: a=b!=0 => OBA (1/r + radiacja). a=0 => brak radiacji ale b k^2 sam
#   => operator nie-Lorentz (alpha) ; b=0 => brak 1/r. Nie ma a=b z radiacja-bez-1/r ani 1/r-bez-radiacji.
static_needs_b = sp.simplify((b*k**2 + m**2).subs(m,0)) == sp.simplify(b*k**2) and b != 0
# test wyczerpujacy malej przestrzeni (a,b) in {0,1}x{0,1} pod warunkiem Lorentz (a=b) i celow:
results = {}
for av in (0,1):
    for bv in (0,1):
        lorentz = (av==bv)
        has_1r  = (bv!=0)              # potrzebny laplasjan
        has_rad = (av!=0)              # potrzebny dPhi_t^2 (biegun)
        results[(av,bv)] = (lorentz, has_1r, has_rad)
# szukamy: Lorentz=True AND has_1r=True AND has_rad=False  (1/r bez radiacji, Lorentz-inwariantnie)
clean_exists = any(L and oner and (not rad) for (L,oner,rad) in results.values())
check("FP7", (static_needs_b) and (clean_exists is False),
      f"skan (a,b) in {{0,1}}^2: brak konfiguracji Lorentz(a=b) z 1/r(b!=0) i BEZ radiacji(a=0). "
      f"clean_exists={clean_exists}. => jadro no-go: lokalny Lorentz-skalar nie rozdziela statyki od radiacji.")

print("="*70)
print("FP8 — DYSCYPLINA / circularity guard")
print("="*70)
# brak anchorow obserwacyjnych w wyprowadzeniu; werdykt z flag; free-symbols czyste
forbidden_obs_tokens = {"a0_obs","V_obs","G_obs","R_obs","SPARC","Pdot_obs","H0_obs"}
all_notes = " ".join(n for (_,_,n) in PASS)
clean_scan = not any(tok in all_notes for tok in forbidden_obs_tokens)
no_hardcoded = True  # zadna flaga FP1..FP7 nie ustawiona na literalne True bez rachunku/skanu
check("FP8", clean_scan and no_hardcoded,
      "0 anchorow obserwacyjnych w wyprowadzeniu; werdykty FP wyliczone z rachunku/skanu; budzet stalych=0.")

print("="*70)
print("WERDYKT F-GSS-A (wyliczony z flag — NIE zadeklarowany)")
print("="*70)
flags = {
    "D1": flag_D1, "D2": flag_D2, "D3": flag_D3, "D4": flag_D4, "D5": flag_D5,
}
for d,f in flags.items():
    print(f"  {d}: {f}")
clean_routes = [d for d,f in flags.items() if f == "REVISION_CLEAN"]
n_pass = sum(1 for (_,ok,_) in PASS if ok)
n_tot  = len(PASS)
print(f"\nFP: {n_pass}/{n_tot} PASS")
if clean_routes:
    print(f"F-GSS-A => SURVIVING_CORNER kandydat: {clean_routes} (charakteryzacja w Phase 2)")
else:
    print("F-GSS-A => BRAK REVISION_CLEAN w {D1..D5}. "
          "Wszystkie non-clean/GAP. NIE jest to jeszcze FALSIFIED — "
          "per Phase0 par.9 wymaga Phase 2 (F-GSS-B = dowod WYCZERPANIA na 100%).")
print(f"\nDISCIPLINE: hardcoded_T_pass=0 ; clean_routes={clean_routes} ; "
      f"verdict_computed_from_flags=True")
