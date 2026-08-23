# AUDIT4_localization.py
# ZARZUT 2: czy "trawers SKUPIONY" w BUDGET_psi_orientation [D] to wynik czy artefakt?
#
# Po ustaleniu h==1 (AUDIT4_h_decoupling) zadanie autora redukuje sie DOKLADNIE do:
#      min  SUM_{i=1..M} w_i * d_i        d_i = |df_i| >= 0
#      s.t. SUM d_i = pi,  d_i <= dmax = (B-2)/c
# To jest LINIOWY program z rosnaca waga w_i -> rozwiazanie ZACHLANNE:
# wypelnij do capu najmniejsze w_i.  "Lokalizacja" = wybor MONOTONICZNEJ wagi, nic wiecej.
#
# TESTY:
#  [1] zaleznosc od WAGI w (1, i, i^2, i^3) - czy skupienie przezywa w=1?
#  [2] zaleznosc od M (rozmiar domeny) - czy promien rdzenia dryfuje z pudlem?
#  [3] zaleznosc od WYKLADNIKA kosztu p (L1 vs L2) - czy skupienie przezywa gradient kwadratowy?
#  [4] czy WERDYKT SLOWNY autora ('SKUPIONY'/'ROZLOZONY') flipuje z M?
#  [5] granica kontinuum: promien FIZYCZNY przy ustalonej domenie a stale gestszej siatce

import numpy as np
np.set_printoptions(suppress=True)

TOT = np.pi

def solve_L1(w, dmax, tot=TOT):
    """min sum w_i d_i, sum d=tot, 0<=d<=dmax -> zachlannie na najmniejszych w"""
    M = len(w); d = np.zeros(M)
    order = np.argsort(w, kind='stable')
    rem = tot
    for j in order:
        take = min(dmax, rem); d[j] = take; rem -= take
        if rem <= 1e-15: break
    return (d, rem <= 1e-9)

def solve_L2(w, dmax, tot=TOT):
    """min sum w_i d_i^2, sum d=tot, 0<=d<=dmax -> water-filling d_i=min(dmax, lam/(2w_i))"""
    w = np.asarray(w, float)
    if np.sum(np.full(len(w), dmax)) < tot - 1e-12: return (None, False)
    lo, hi = 0.0, 1e12
    for _ in range(300):
        lam = 0.5*(lo+hi)
        d = np.minimum(dmax, lam/(2*w))
        if d.sum() < tot: lo = lam
        else: hi = lam
    return (np.minimum(dmax, 0.5*(lo+hi)/(2*w)), True)

def weights(kind, M):
    i = np.arange(1, M+1).astype(float)
    return {"w=1": np.ones(M), "w=i": i, "w=i^2": i**2, "w=i^3": i**3}[kind]

def r_q(d, q):
    """indeks promieniowy, przy ktorym skumulowane nawiniecie osiaga q*total"""
    c = np.cumsum(d); tot = c[-1]
    k = int(np.searchsorted(c, q*tot)) + 1
    return k

dmax = 1.0   # = (B-2)/c dla B=3, c=1 (wartosci autora)

print("="*78)
print("[1]+[2] PROMIEN RDZENIA vs WAGA i vs M   (koszt L1 - dokladnie jak u autora)")
print("-"*78)
print(f"  dmax=(B-2)/c={dmax},  trawers=pi.  r90 = ile krokow miesci 90% nawiniecia")
print(f"  {'waga':>7} | " + " | ".join(f"M={M:<4d} r90 (akt.)" for M in [10,20,40,80,160]))
for kind in ["w=1","w=i","w=i^2","w=i^3"]:
    row=[]
    for M in [10,20,40,80,160]:
        w = weights(kind, M); d,ok = solve_L1(w, dmax)
        row.append(f"{r_q(d,0.9):>6d} ({int((d>1e-9).sum()):>3d})" if ok else "   INFEAS   ")
    print(f"  {kind:>7} | " + " | ".join(row))
print("  UWAGA: przy w=1 WSZYSTKIE rozwiazania sa optymalne (cel = const = pi).")
print("         Powyzsze 'r90' dla w=1 to tylko to, co zwraca moj tie-break - LP jest ZDEGENEROWANY.")

print("\n"+"="*78)
print("[1b] DEGENERACJA przy w=1 (jednorodna miara): czy wynik jest w ogole okreslony?")
print("-"*78)
M=20; w=np.ones(M)
for name,d in [("zachlanne (skrajnie skupione)", np.array([1.,1.,1.,0.1416]+[0.]*16)),
               ("jednorodne (rozlane)",          np.full(M,np.pi/M)),
               ("na BRZEGU domeny",              np.array([0.]*16+[0.1416,1.,1.,1.]))]:
    print(f"    {name:32s} koszt = {float(np.sum(w*d)):.6f}   suma={d.sum():.4f}")
print("  >>> IDENTYCZNY koszt. Bez wagi rosnacej model NIE lokalizuje w ogole.")
print("      Cala 'lokalizacja' pochodzi WYLACZNIE z wyboru w = i^2.")

print("\n"+"="*78)
print("[3] KOSZT KWADRATOWY (jakikolwiek czlon gradientowy) - czy skupienie przezywa?")
print("-"*78)
print(f"  {'waga':>7} | " + " | ".join(f"M={M:<4d} r90 (akt.)" for M in [10,20,40,80,160]))
for kind in ["w=1","w=i","w=i^2","w=i^3"]:
    row=[]
    for M in [10,20,40,80,160]:
        w=weights(kind,M); d,ok=solve_L2(w,dmax)
        row.append(f"{r_q(d,0.9):>6d} ({int((d>1e-9).sum()):>3d})" if ok else "   INFEAS   ")
    print(f"  {kind:>7} | " + " | ".join(row))
print("  >>> L2: przy w=1 i w=i promien r90 ROSNIE Z M => ARTEFAKT DOMENY (pudlo).")
print("      Tylko w=i^2 / i^3 (suma 1/w zbiezna) daje r90 niezalezne od M.")
print("      Czyli 'lokalizacja' zalezy od NIETESTOWANEGO wyboru wagi i wykladnika.")

print("\n"+"="*78)
print("[4] CZY WERDYKT SLOWNY AUTORA FLIPUJE Z M?  (regula: aktywne >= M-1 => 'ROZLOZONY')")
print("-"*78)
print("  autor: if np.sum(df_s>1e-4) >= M-1:  'trawers ROZLOZONY ... rozmiar wyznaczaja SASIEDZI'")
print("         else:                          'trawers SKUPIONY - budzet lokalizuje sam z siebie'")
print(f"  {'M':>4} {'aktywnych':>10} {'M-1':>5}  WERDYKT WYPISANY PRZEZ SKRYPT AUTORA")
for M in [4,5,6,8,10,20,40]:
    w=weights("w=i^2",M); d,ok=solve_L1(w,dmax)
    if not ok:
        print(f"  {M:4d} {'-':>10} {M-1:5d}  NIEDOPUSZCZALNE (M*dmax < pi)"); continue
    act=int((d>1e-4).sum())
    verdict = "ROZLOZONY (rozmiar od SASIADOW)" if act>=M-1 else "SKUPIONY (budzet lokalizuje)"
    print(f"  {M:4d} {act:10d} {M-1:5d}  {verdict}")
print("  >>> Werdykt SLOWNY flipuje przy M=4..5, choc fizyka (4 aktywne kroki) sie NIE zmienia.")
print("      Autor uruchomil TYLKO M=20. Test na 'skupienie' nie byl testem, tylko M-zaleznym progiem.")

print("\n"+"="*78)
print("[5] GRANICA KONTINUUM: czy obiekt ma ROZMIAR FIZYCZNY, czy rozmiar SIATKI?")
print("-"*78)
print("  Ustalam domene fizyczna L=1. Odstep wezlowy a=L/M. Budzet NA WEZEL staly => dmax staly.")
print("  Promien fizyczny rdzenia R_phys = a * r90 = r90/M.")
print(f"  {'M':>5} {'r90 [wezly]':>13} {'R_phys/L':>12}")
for M in [10,20,40,80,160,320,640]:
    w=weights("w=i^2",M); d,ok=solve_L1(w,dmax)
    print(f"  {M:5d} {r_q(d,0.9):13d} {r_q(d,0.9)/M:12.5f}")
print("  >>> R_phys/L ~ 1/M -> 0.  Obiekt to ZAWSZE ~pi/dmax KOMOREK SIATKI.")
print("      To nie jest skala dynamiczna - to obcieciowy (UV) artefakt dyskretyzacji.")
print("      Nie ma granicy kontinuum; kazdy obiekt ma ten sam rozmiar ~3 wezly,")
print("      niezaleznie od B_topologicznego, ladunku, czy czegokolwiek innego.")

print("\n"+"="*78)
print("[6] KONTROLA: czy rozmiar zalezy od ladunku topologicznego? (trawers B*pi)")
print("-"*78)
print(f"  {'B_top':>6} {'trawers':>9} {'r90 [wezly]':>13} {'r90/r90(B=1)':>14} {'B^(1/3)':>9}")
M=400; w=weights("w=i^2",M)
base=None
for Bt in [1,2,3,8,27]:
    d,ok=solve_L1(w,dmax,tot=Bt*np.pi); r=r_q(d,0.9)
    if base is None: base=r
    print(f"  {Bt:6d} {Bt*np.pi:9.4f} {r:13d} {r/base:14.4f} {Bt**(1/3):9.4f}")
print("  >>> w tym modelu r ~ B LINIOWO (nie B^(1/3)) - zgodnie z wlasna korekta autora.")
