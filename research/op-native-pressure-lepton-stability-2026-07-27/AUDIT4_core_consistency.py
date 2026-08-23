# AUDIT4_core_consistency.py
# Konfrontacja calej warstwy budzetowej z RDZENIEM: core/sek08c, prop:antipodal-from-budget.
#
# CYTATY (dowod prop. antypodycznej, linie 319-347 sek08c):
#   (i)  "Budzet informacyjny narzuca ADDYTYWNY wiaz:  n_sp + n_czas = B."    <- ROWNOSC
#   (ii) "n_sp * n_czas = n_0 * n_0t = const"                                  <- ROWNOSC
#   (iii) teza prop.: "Jezeli budzet B jest dzielony WYLACZNIE miedzy przestrzen (h)
#         i czas (f), PRZY BRAKU UKRYTYCH STOPNI SWOBODY, ORAZ dekompozycja jest
#         multiplikatywna, to fh = 1."

import numpy as np
print("="*76)
print("[1] NADOKRESLENIE: rdzen naklada OBIE ROWNOSCI naraz. Co z tego wynika?")
print("-"*76)
print("  n_sp + n_czas = B   i   n_sp * n_czas = C   =>   n_sp, n_czas = pierwiastki")
print("  x^2 - Bx + C = 0  =>  OBA SA STALE  =>  h = n_sp/n0 = CONST.")
B_core, C_core = 2.0, 1.0    # w normalizacji n0=n0t=1: B=2, C=1
r = np.roots([1,-B_core,C_core])
print(f"  przy n0=n0t=1 (czyli B=2, C=1): pierwiastki = {np.round(r,6)}  => n_sp=n_czas=1")
print(f"  => h == 1 IDENTYCZNIE.  Pole Phi(x) NIE MOZE SIE ZMIENIAC.")
print("  >>> Rdzen, wziety doslownie, jest NADOKRESLONY: obie rownosci naraz zabijaja Phi.")
print("      Autor po cichu zamienia (i) na NIEROWNOSC (n_sp+n_czas <= B).")
print("      To NAPRAWIA rdzen, ale caly 'cap' jest wtedy KONSTRUKCJA AUTORA, nie rdzeniem.")
print("      (autor zaznacza 'moja interpretacja' - ale nie zaznacza, ze wersja rdzenia")
print("       jest wewnetrznie sprzeczna z istnieniem zmiennego Phi.)")

print("\n"+"="*76)
print("[2] SPRZECZNOSC W KOREKCIE: n_orient LAMIE PRZESLANKE, ktora daje fh=1")
print("-"*76)
print("  BUDGET_psi_orientation.py (linie 11-15) pisze:")
print("      'budzet dzieli sie na TRZY czesci:  n_sp + n_czas + n_orient <= B'")
print("      'z fh=1:  n_sp = h,  n_czas = 1/h  =>  g(h)=h+1/h'")
print("  Ale prop:antipodal-from-budget wymaga, by budzet byl dzielony")
print("  'WYLACZNIE miedzy przestrzen i czas, PRZY BRAKU UKRYTYCH STOPNI SWOBODY'.")
print("  n_orient JEST dodatkowym stopniem swobody, ktory zjada budzet.")
print("  >>> Wprowadzenie n_orient UNIEWAZNIA fh=1, a wiec i g(h)=h+1/h,")
print("      a wiec i 'g_min=2 przy h=1', a wiec i caly wzor |df|_max=(B-2)/c,")
print("      ktory jest reklamowany jako JEDYNY ZYSK NETTO korekty.")
print("      Wynik uzywa przeslanki, ktora sam falsyfikuje.")

print("\n"+"="*76)
print("[3] 'PROZNIA = MINIMUM BUDZETU' - czy to wynik, czy wybor n0 = n0t?")
print("-"*76)
print("  g(h) = n0*h + n0t/h  ma minimum w  h* = sqrt(n0t/n0), NIE w h=1.")
print("  h*=1 wymaga n0 = n0t, czyli: w prozni tyle samo budzetu na PRZESTRZEN co na CZAS.")
print("  Rdzen mowi tylko n0 = n_sp(Phi0), n0t = n_czas(Phi0) - nic o ich rownosci.")
print(f"\n  {'n0':>6} {'n0t':>6} {'h* = argmin g':>15} {'g_min':>9} {'psi*(M9.1\"\")':>14}  interpretacja")
def psi_of_h(h): return 4*h/(1+3*h)
for n0,n0t,tag in [(1,1,"wybor autora (arbitralny)"),
                   (3,1,"3 kierunki przestrz. vs 1 czasowy"),
                   (1,3,"odwrotnie"),
                   (3,3,"n0=n0t, inna skala"),
                   (4,1,"4 dof czasoprzestrzeni?")]:
    hs=np.linspace(1e-3,60,2000000); g=n0*hs+n0t/hs; i=np.argmin(g)
    print(f"  {n0:6.1f} {n0t:6.1f} {hs[i]:15.4f} {g[i]:9.4f} {psi_of_h(hs[i]):14.4f}  {tag}")
print("  >>> Minimum wypada w prozni (h=1) WYLACZNIE gdy n0 = n0t.")
print("      Przy n0=3, n0t=1 minimum jest w h=0.577 - PROZNIA NIE JEST minimum budzetu,")
print("      a stan najtanszy to SKURCZONA przestrzen. Wynik z WYNIKI_budzet-skala sek.2")
print("      ('efekt uboczny o duzej wadze: proznia = minimum budzetu') to WYBOR NORMALIZACJI.")

print("\n"+"="*76)
print("[4] KONSEKWENCJA DLA 'CAPU': jak bardzo psi_max zalezy od wolnych wyborow?")
print("-"*76)
print(f"  {'n0':>5} {'n0t':>5} {'B':>6} {'h_max':>9} {'psi_max (M9.1\"\")':>17}")
for n0,n0t in [(1,1),(3,1),(1,3)]:
    for Bv in [3.0,6.0]:
        disc=Bv*Bv-4*n0*n0t
        if disc<0: print(f"  {n0:5.1f} {n0t:5.1f} {Bv:6.1f} {'BRAK (B<g_min)':>9}"); continue
        hmax=(Bv+np.sqrt(disc))/(2*n0)
        print(f"  {n0:5.1f} {n0t:5.1f} {Bv:6.1f} {hmax:9.4f} {psi_of_h(hmax):17.4f}")
print("  >>> psi_max zalezy od (n0, n0t, B) - trzech nieustalonych liczb.")
print("      'Cap' istnieje jakosciowo (AM-GM), ale jego WARTOSC nie jest niczym wyznaczona.")
