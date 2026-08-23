# PRESSURE_wigner_seitz.py
# TEZA (uzytkownik): pojedynczy soliton NIE moze istniec; istnieje tylko dzieki
# cisnieniu innych. Sciana = dowod, ze konfiguracja wielu obiektow jest
# KORZYSTNIEJSZA ENERGETYCZNIE niz jej rozpad/kolaps.
#
# Realizacja: komorka Wignera-Seitza. Obiekt w kuli promienia L, na brzegu
# u'(L)=0 (symetria z sasiadami) ZAMIAST u(inf)=1 (proznia). To usuwa premise
# Derricka/testu izolowanego (audyt DODATEK A liczyl wlasnie u(inf)->1).
#
# Sektor: KANONICZNY F-A (alpha=2, phi^4), z audytu:
#   E = 4pi INT [ 1/2 u^4 u'^2 + W(u) ] r^2 dr,  W(u)=u^8/8 - u^7/7
#   EL:  u'' + 2u'/r = u^3 - u^2 - 2u'^2/u
#
# ANTY-FIT (zadeklarowane PRZED liczeniem):
#   - zero skanowania parametru do celu; pytamy TYLKO o istnienie minimum E(L)
#   - falsyfikator: brak minimum przy skonczonym L => teza pada w tym sektorze
#   - zero dopasowania do mas/PDG

import numpy as np
from scipy.integrate import solve_ivp, trapezoid

def W(u):  return u**8/8 - u**7/7
def Wp(u): return u**7 - u**6          # = u^6(u-1)
def Wpp(u):return 7*u**6 - 6*u**5      # = u^5(7u-6)

print("="*66)
print("[0] STRUKTURA POTENCJALU KANONICZNEGO W(u)=u^8/8-u^7/7")
print("-"*66)
us=np.linspace(1e-6,3,300000)
i=np.argmin(W(us))
print(f"  globalne minimum W na u>0:  u*={us[i]:.6f}   W(u*)={W(us[i]):.8f}")
print(f"  W''(1) = {Wpp(1.0):+.4f}  (>0 => proznia u=1 STABILNA)")
print(f"  W'(u)=u^6(u-1): jedyne zero na u>0 to u=1 => JEDNA studnia")
print(f"  czy W(u) >= W(1) dla wszystkich u>0? {np.all(W(us) >= W(1.0)-1e-12)}")

print("\n"+"="*66)
print("[1] TWIERDZENIE ENERGETYCZNE (analityczne, decydujace)")
print("-"*66)
print("  E[u] = 4pi INT [ 1/2 u^4 u'^2 + (W(u)-W(1)) ] r^2 dr")
print("  - czlon kinetyczny 1/2 u^4 u'^2 >= 0  (dodatnio okreslony)")
print("  - W(u)-W(1) >= 0  (u=1 jest GLOBALNYM minimum, sprawdzone wyzej)")
print("  => E[u] >= 0, z rownoscia TYLKO dla u == 1 (jednorodna proznia).")
print("  WNIOSEK: w sektorze F-A zadna niejednorodna konfiguracja NIE moze byc")
print("           korzystniejsza energetycznie od jednorodnej - przy ZADNEJ gestosci.")

print("\n"+"="*66)
print("[2] TEST NUMERYCZNY: czy w komorce W-S istnieje NIEJEDNORODNE rozwiazanie?")
print("-"*66)
def rhs(r,y):
    u,up=y
    u=max(u,1e-9)
    return [up, u**3-u**2-2*up**2/u - (2*up/r if r>1e-12 else 0.0)]

def shoot(g0,L):
    """start u(0)=g0, u'(0)=0; zwraca u'(L) (szukamy zera = Neumann na brzegu)"""
    s=solve_ivp(rhs,[1e-8,L],[g0,0.0],rtol=1e-10,atol=1e-12,dense_output=True,max_step=L/200)
    if s.status!=0 or not np.isfinite(s.y[1,-1]): return np.nan, s
    return s.y[1,-1], s

print("  {:>5} {:>7} {:>12} {:>12}  komentarz".format("L","g0","u(L)","up(L)"))
for L in [1.0, 2.0, 5.0, 10.0]:
    for g0 in [0.5, 0.8, 0.95, 1.05, 1.2, 2.0]:
        upL,s=shoot(g0,L)
        uL=s.y[0,-1] if np.isfinite(s.y[0,-1]) else np.nan
        note=""
        if not np.isfinite(upL): note="rozbieglo (runaway)"
        elif abs(upL)<1e-6: note="<-- Neumann spelniony"
        else: note="u' != 0 na brzegu"
        print(f"  {L:5.1f} {g0:7.3f} {uL:12.4g} {upL:12.4g}  {note}")

print("\n  Analiza mechaniczna (dlaczego): rownanie to czastka w potencjale -W")
print("  z tarciem 2u'/r. u=1 jest MINIMUM W => MAKSIMUM -W => rownowaga NIESTABILNA.")
print("  Czastka wypuszczona w spoczynku przy u!=1 stacza sie MONOTONICZNIE i nigdy")
print("  nie wraca do u'=0. => brak niejednorodnego rozwiazania w komorce, dla kazdego L.")

print("\n"+"="*66)
print("[3] E(L) - czy istnieje minimum przy skonczonym L?")
print("-"*66)
print("  Jedyne rozwiazanie w komorce to u==1 (jednorodne), dla ktorego E(L)=0")
print("  identycznie dla KAZDEGO L. => E(L) jest STALA, nie ma minimum przy")
print("  skonczonym L. FALSYFIKATOR ZADEKLAROWANY W NAGLOWKU: teza PADA w tym sektorze.")

print("\n"+"="*66)
print("[4] KONTROLA DRUGIEGO SEKTORA: M9.1'' (gravity) z ELEMENTEM OBJETOSCIOWYM")
print("-"*66)
# sek08a/sek08c: V_M911(psi) = -gamma*psi^2*(4-3psi)^2/12 ; sqrt(-g)=c0*psi/(4-3psi)
# efektywny potencjal w energii = sqrt(-g)*V:
#   V_eff = -c0*gamma*psi^3*(4-3psi)/12
ps=np.linspace(1e-6,2.5,400000)
Veff=-(ps**3)*(4-3*ps)/12.0          # c0=gamma=1
j=np.argmin(Veff)
print(f"  V_eff(psi) = -psi^3(4-3psi)/12   (c0=gamma=1)")
print(f"  globalne minimum: psi*={ps[j]:.6f}  V_eff={Veff[j]:.8f}")
print(f"  V_eff'(psi) = -psi^2(1-psi) => zera w psi=0 i psi=1 (psi=1 krytyczne) ")
print(f"  czy V_eff(psi) >= V_eff(1) dla wszystkich psi>0? {np.all(Veff >= Veff[np.argmin(abs(ps-1))]-1e-12)}")
print("  => TAKZE jedno globalne minimum przy psi=1. Twierdzenie z [1] stosuje sie")
print("     do OBU kanonicznych sektorow (F-A materii i M9.1'' grawitacji).")

print("\n"+"="*66)
print("WERDYKT:")
print("  Cisnienie sasiadow NIE ratuje lokalizacji w kanonicznym sektorze SKALARNYM F-A.")
print("  Przeszkoda NIE jest warunkiem brzegowym ani Derrickiem - jest STRUKTURALNA:")
print("  W ma JEDNO globalne minimum i czlon kinetyczny jest dodatni, wiec jednorodna")
print("  proznia jest globalnym minimum energii przy kazdej gestosci.")
print("  Aby cokolwiek bylo zlokalizowane, potrzebny jest WIEZ zabraniajacy u==1:")
print("  ladunek topologiczny (sektor wstazek 2T) albo zrodlo materii rho.")
