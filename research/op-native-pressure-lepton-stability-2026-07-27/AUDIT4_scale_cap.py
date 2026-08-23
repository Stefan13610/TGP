# AUDIT4_scale_cap.py
# ZARZUT 3 + ZARZUT 6: BUDGET_scale_cap.py / WYNIKI_budzet-skala.
#  [1] reprodukcja liczb autora
#  [2] czy R_min ~ Q^(1/3) ma JAKAKOLWIEK tresc TGP, czy jest tozsamoscia wymiarowa
#  [3] czy 'psi_max < 4/3 = horyzont' jest wynikiem budzetu, czy algebra
#  [4] SAMOZNISZCZENIE: przy h==1 (wlasna korekta autora) Q = INT(psi-1)dV = 0 => R_min = 0
#  [5] niespojnosc znaku: budzet-skala zaklada GARB (psi>1), profil-z-budzetu DOLEK (psi<1)

import numpy as np
rng = np.random.default_rng(3)
n0=n0t=1.0

def h_range(B):
    disc=B*B-4*n0*n0t
    if disc<0: return None
    return ((B-np.sqrt(disc))/(2*n0),(B+np.sqrt(disc))/(2*n0))
def psi_from_h(h): return 4*h/(1+3*h)

print("="*76)
print("[1] REPRODUKCJA TABELI AUTORA (WYNIKI_budzet-skala sek.1)")
print("-"*76)
print(f"  {'B':>6} {'h_min':>9} {'h_max':>9} {'psi_max formaI':>15} {'psi_max M9.1\"\"':>15}")
for B in [2.5,3.0,4.0,10.0,100.0]:
    r=h_range(B)
    print(f"  {B:6.1f} {r[0]:9.4f} {r[1]:9.4f} {r[1]:15.4f} {psi_from_h(r[1]):15.4f}")
print("  >>> zgadza sie z dokumentem co do 4 miejsc. Arytmetyka OK.")

print("\n"+"="*76)
print("[3] CZY 'psi_max SCISLE PONIZEJ 4/3' JEST WYNIKIEM BUDZETU?")
print("-"*76)
print("  psi(h) = 4h/(1+3h).  sup_{h>0} psi = 4/3, osiagane dopiero przy h=inf.")
print("  Wiec psi < 4/3 dla KAZDEGO skonczonego h - niezaleznie od jakiegokolwiek budzetu:")
for h in [0.5,1,10,100,1e6,1e12]:
    print(f"    h={h:>10.4g}  ->  psi={psi_from_h(h):.10f}   (4/3 = 1.3333333333)")
print("  >>> 'cap jest realny, bo psi_max < 4/3' to TOZSAMOSC ALGEBRAICZNA funkcji 4h/(1+3h),")
print("      a nie konsekwencja wiezu budzetowego. Zero tresci.")

print("\n"+"="*76)
print("[2] CZY R_min ~ Q^(1/3) MA TRESC? test uniwersalnosci")
print("-"*76)
print("  Teza: dla DOWOLNEGO ksztaltu chi, DOWOLNEGO capu A_max, DOWOLNEJ teorii,")
print("  z 'Q = INT (pole-1) dV zachowane' + 'amplituda <= A_max' wynika R >= (Q/(c A_max))^(1/3).")
print("  Sprawdzam na 8 losowych ksztaltach chi i losowych capach:")
r=np.linspace(0,1,4001)
print(f"  {'ksztalt':>9} {'c=INT chi dV':>13} {'A_max':>8} {'Q':>8} {'R_min':>9} {'(Q/(cA))^(1/3)':>16} {'blad':>10}")
for k in range(8):
    a,b=rng.uniform(0.5,6),rng.uniform(0.5,6)
    chi=(1-r**a)**b                       # dowolny malejacy ksztalt
    cshape=float(np.trapezoid(4*np.pi*r**2*chi,r))
    Amax=rng.uniform(0.05,2.0); Q=rng.uniform(0.1,50)
    Rmin=(Q/(cshape*Amax))**(1/3)
    # niezalezna weryfikacja: najmniejsze R przy ktorym A=Q/(c R^3) <= A_max
    Rs=np.linspace(1e-3,20,400000); A=Q/(cshape*Rs**3)
    Rnum=Rs[np.argmax(A<=Amax)]
    print(f"  {k:9d} {cshape:13.5f} {Amax:8.4f} {Q:8.3f} {Rnum:9.4f} {Rmin:16.4f} {abs(Rnum-Rmin):10.2e}")
print("  >>> IDENTYCZNIE dla kazdego ksztaltu/capu. To PRZEKSZTALCENIE ALGEBRAICZNE")
print("      wzoru Q = c*A*R^3, nie predykcja fizyczna. TAUTOLOGIA.")
print("      (autor to CZESCIOWO przyznaje w sek.4, ale i tak umieszcza jako 'PREDYKCJA STRUKTURALNA'")
print("       i jako wiersz 'skala: R ~ B^(1/3)' w tabeli 'stanu wiezy'.)")

print("\n"+"="*76)
print("[4] SAMOZNISZCZENIE PO WLASNEJ KOREKCIE AUTORA (to NIE zostalo wycofane)")
print("-"*76)
print("  KOREKTA (WYNIKI_psi-orientacja) ustala: h == 1 na CALYM profilu,")
print("  'obiekt jest CZYSTA TEKSTURA ORIENTACJI, metryka plaska'.")
print("  Forma I: h=psi => psi == 1 wszedzie.  M9.1'': psi=4h/(1+3h)=4/4=1 wszedzie.")
print(f"    psi(h=1) forma I = {1.0:.6f};   psi(h=1) M9.1'' = {psi_from_h(1.0):.6f}")
print("  Ale Q z BUDGET_scale_cap to Q = INT (psi - 1) dV.")
Qafter=float(np.trapezoid(4*np.pi*r**2*(1.0-1.0)*np.ones_like(r),r))
print(f"    => Q = INT (psi-1) dV = {Qafter:.1f}")
print("    => R_min = (Q/(c*A_max))^(1/3) = 0.  Cap nie ma czego ograniczac.")
print("  >>> CALY MECHANIZM 'cap blokuje kolaps Derricka' (sek.3 WYNIKI_budzet-skala)")
print("      jest po korekcie PUSTY, nie tylko 'zastapiony'. A dokument NIE MA")
print("      naglowka wycofania - w przeciwienstwie do dwoch pozostalych.")

print("\n"+"="*76)
print("[5] NIESPOJNOSC ZNAKU MIEDZY DOKUMENTAMI (garb vs dolek)")
print("-"*76)
print("  BUDGET_scale_cap:  psi = 1 + A*chi,  A_max = psi_max - 1  > 0  => GARB (psi>1)")
print("  BUDGET_profile:    psi_core = 0.7 < 1                          => DOLEK (psi<1)")
print("  To dwa PRZECIWNE obiekty. Wlasciwy cap dla dolka to |A| <= 1 - h_min:")
print(f"  {'B':>6} {'A_max (garb)':>13} {'|A|_max (dolek)':>16} {'R_min garb':>12} {'R_min dolek':>13} {'stosunek':>10}")
for B in [2.5,3.0,4.0,10.0]:
    lo,hi=h_range(B); Ag=hi-1.0; Ad=1.0-lo; Q=1.0
    Rg=(Q/Ag)**(1/3); Rd=(Q/Ad)**(1/3)
    print(f"  {B:6.1f} {Ag:13.4f} {Ad:16.4f} {Rg:12.4f} {Rd:13.4f} {Rd/Rg:10.4f}")
print("  >>> 'R_min' rozni sie o czynnik ~1.4-2 w zaleznosci od tego, ktory dokument czytamy.")
print("      'Predykcja' nie jest nawet jednoznacznie zdefiniowana.")

print("\n"+"="*76)
print("[6] CZY CAP PRZEZYWA BEZ ZALOZENIA MULTIPLIKATYWNEGO fh=1?")
print("-"*76)
print("  Autor przyznaje: rdzen daje wiez ADDYTYWNY n_sp+n_czas<=B; fh=1 to DODATEK.")
print("  Bez fh=1: n_sp <= B => h <= B/n0  (cap Z GORY nadal jest, TRYWIALNIE),")
print("  ale n_czas moze isc do 0 => h nieograniczone Z DOLU => brak dolnego capu.")
print("  Z fh=1: g(h)=h+1/h, cap obustronny. Czyli OBUSTRONNOSC capu (jedyna nietrywialna")
print("  czesc wyniku) pochodzi WYLACZNIE z zalozenia, ktore rdzen odrzuca jako niewynikajace.")
