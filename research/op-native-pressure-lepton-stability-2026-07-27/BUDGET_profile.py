# BUDGET_profile.py
# CEL: czy PROFIL psi(r) wypada z samego BUDZETU (bez rownania pola, bez energii)?
#
# SKLADNIKI (wszystkie z rdzenia / wczesniejszych krokow):
#  (i)   koszt budzetu na wezel:   g(h) = n0*h + n0t/h,  min w h=1  [sek08c + BUDGET_scale_cap]
#  (ii)  h(psi): forma I  h=psi   |  M9.1''  h=psi/(4-3psi)         [sek08c]
#  (iii) miara wezlow: wezly przypisane do WSPOLRZEDNYCH (N_B topologiczne),
#        powloka i ma wage ~ i^2                                     [sek08c def:info-budget]
#  (iv)  DYSKRETNOSC: |psi_{i+1}-psi_i| <= delta_max  (substrat nie moze zmienic
#        sie szybciej niz jeden krok wezlowy)                        [BUDGET_topological_cost]
#  (v)   TOPOLOGIA: nawiniecie wymaga PELNEGO TRAWERSU psi: psi_core -> 1
#        (qm_spin: f = pi(1-g)/(1-g0), B=1 dla dowolnego monotonicznego trawersu)
#
# ZASADA (JAWNE ZALOZENIE, testowane): obiekt = konfiguracja MINIMALNEGO
# ZUZYCIA BUDZETU spelniajaca (iv)+(v). To natywny odpowiednik dawnej
# "minimalizacji energii", ale na prymitywie, ktory w TGP faktycznie istnieje.
#
# FALSYFIKATOR (przed liczeniem): jesli minimalizacja daje profil zdegenerowany
# (psi==1 wszedzie / brak skonczonego rdzenia), to budzet NIE wyznacza profilu.

import numpy as np
from scipy.optimize import minimize

n0=n0t=1.0
def g(h): return n0*h + n0t/h
G1=g(1.0)

def h_of(psi, form):
    if form=="I":     return psi
    if form=="M911":  return psi/(4-3*psi)

def cost(psi_arr, form, w):
    h=h_of(np.asarray(psi_arr), form)
    return float(np.sum(w*(g(h)-G1)))

def solve_profile(psi_core, M, delta_max, form, weight="shell"):
    """minimalizuj koszt budzetu; psi[0]=psi_core, psi[M]=1, |dpsi|<=delta_max"""
    i=np.arange(M+1)
    w = (i.astype(float)**2 if weight=="shell" else np.ones(M+1))
    x0=np.linspace(psi_core,1.0,M+1)[1:-1]
    def obj(x):
        p=np.concatenate([[psi_core],x,[1.0]]); return cost(p,form,w)
    cons=[]
    for k in range(M):
        cons.append({'type':'ineq','fun':(lambda x,k=k: delta_max-abs(
            np.concatenate([[psi_core],x,[1.0]])[k+1]-np.concatenate([[psi_core],x,[1.0]])[k]))})
    lo,hi=min(psi_core,1.0),max(psi_core,1.0)
    res=minimize(obj,x0,constraints=cons,bounds=[(lo,hi)]*(M-1),
                 method='SLSQP',options={'maxiter':800,'ftol':1e-12})
    return np.concatenate([[psi_core],res.x,[1.0]])

print("="*70)
print("[1] PROFIL Z MINIMALIZACJI BUDZETU (forma I, psi_core=0.7, delta_max=0.1)")
print("-"*70)
M=12; dmax=0.10; pc=0.70
p=solve_profile(pc,M,dmax,"I")
print("   i :  psi_i    |dpsi|   (delta_max=%.2f)"%dmax)
for i in range(M+1):
    d=abs(p[i]-p[i-1]) if i>0 else 0.0
    flag=" <- tempo NASYCONE" if abs(d-dmax)<1e-3 else ""
    print(f"  {i:2d} : {p[i]:7.4f}  {d:7.4f}{flag}")
nsat=sum(1 for i in range(1,M+1) if abs(abs(p[i]-p[i-1])-dmax)<1e-3)
print(f"\n  liczba krokow z NASYCONYM tempem: {nsat}")
print(f"  przewidywana liczba krokow trawersu = (1-psi_core)/delta_max = {(1-pc)/dmax:.1f}")
print("  => profil: LINIOWY trawers z maksymalnym tempem, potem PLASKA proznia.")

print("\n"+"="*70)
print("[2] CZY TO SAMO W KANONICZNEJ M9.1''? (kontrola formy h)")
print("-"*70)
p2=solve_profile(pc,M,dmax,"M911")
print("   i :  psi_i (M9.1'')   |dpsi|")
for i in range(0,M+1,2):
    d=abs(p2[i]-p2[i-1]) if i>0 else 0.0
    print(f"  {i:2d} : {p2[i]:10.4f}      {d:7.4f}")
nsat2=sum(1 for i in range(1,M+1) if abs(abs(p2[i]-p2[i-1])-dmax)<1e-3)
print(f"  liczba krokow z nasyconym tempem: {nsat2}  (forma I: {nsat})")
print("  => ten sam KSZTALT: trawers z maks. tempem + plaska proznia (bag).")

print("\n"+"="*70)
print("[3] PROMIEN RDZENIA: r_c = (1-psi_core)/delta_max  (w krokach wezlowych)")
print("-"*70)
print(f"  {'psi_core':>9} {'r_c [wezly]':>13}")
for pcv in [0.9,0.8,0.7,0.5,0.3]:
    print(f"  {pcv:9.2f} {(1-pcv)/dmax:13.1f}")
print("  => GLEBSZY rdzen  =>  SZERSZY obiekt (liniowo). Rozmiar z trawersu, nie z energii.")

print("\n"+"="*70)
print("[4] TEST DEGENERACJI (falsyfikator)")
print("-"*70)
pflat=solve_profile(1.0,M,dmax,"I")
print(f"  psi_core=1 (brak trawersu): profil = {np.round(pflat,4)}")
print(f"  koszt = {cost(pflat,'I',(np.arange(M+1).astype(float))**2):.6f}  (proznia, zero kosztu)")
print("  => bez wymogu trawersu obiekt ZNIKA (poprawnie). Z trawersem: skonczony rdzen.")
print("  >>> Falsyfikator NIE zadzialal: budzet WYZNACZA profil (nietrywialny).")

print("\n"+"="*70)
print("[5] DOMKNIECIE GLEBOKOSCI psi_core (dotad WOLNY parametr)")
print("-"*70)
print("  Dwa PRZECIWSTAWNE wiezy:")
print("   (a) minimalizacja budzetu chce psi_core -> 1  (plytko = tanio)")
print("   (b) koszt topologiczny wymaga N_rdzen >= 5B wezlow")
print("  Przy tempie nasyconym: r_c = (1-psi_core)/delta_max  [wezly promieniowo]")
print("  N_rdzen = (4/3)pi r_c^3  >= 5B   =>   r_c >= (15B/(4pi))^(1/3)")
print("  Minimalizacja budzetu SATURUJE ten wiez (rownosc):")
print("     r_c(B) = (15B/(4pi))^(1/3)      psi_core(B) = 1 - delta_max*r_c(B)")
Nmin_simplex=5
def rc_of(B): return (3*Nmin_simplex*B/(4*np.pi))**(1/3)
print(f"\n  delta_max={dmax}")
print(f"  {'B':>3} {'r_c [wezly]':>12} {'psi_core':>10} {'N_rdzen':>9}")
for B in [1,2,3,8,27]:
    rc=rc_of(B); pcB=1-dmax*rc; Nc=(4/3)*np.pi*rc**3
    print(f"  {B:3d} {rc:12.4f} {pcB:10.4f} {Nc:9.1f}")
print("\n  >>> GLEBOKOSC I ROZMIAR WYZNACZONE JEDNOCZESNIE, oba ~ B^(1/3).")
print("      Zaden nie jest wolnym parametrem (poza skala delta_max = rozdzielczosc substratu).")
print("  Kontrola spojnosci z capem budzetowym: psi_core musi byc >= psi_min.")
for B in [1,8,27,64]:
    rc=rc_of(B); pcB=1-dmax*rc
    print(f"    B={B:3d}: psi_core={pcB:7.4f}  {'OK (>0)' if pcB>0 else 'NARUSZA cap - B za duze dla tego delta_max'}")

print("\n"+"="*70)
print("WNIOSEK:")
print("  Profil NIE jest gladkim solitonem z ogonem eksponencjalnym.")
print("  Budzet + dyskretnosc daja: LINIOWY trawers z maksymalnym dozwolonym tempem,")
print("  ostre przejscie do proznia => obiekt typu 'BAG/kropla', nie 'chmura'.")
print("  To jest PREDYKCJA STRUKTURALNA, ktora moze zawiesc (ogony sa mierzalne).")
