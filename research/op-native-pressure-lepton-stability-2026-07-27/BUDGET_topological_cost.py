# BUDGET_topological_cost.py
# LUKA DO DOMKNIECIA: czym jest zachowane Q w R_min ~ Q^(1/3)?
# HIPOTEZA: Q = LADUNEK TOPOLOGICZNY B (nawiniecie), a zwiazek z budzetem jest
# taki, ze REPREZENTACJA nawiniecia KOSZTUJE WEZLY (substrat jest DYSKRETNY).
#
# ZRODLO: sek08c def:info-budget - "N_B jest liczba TOPOLOGICZNA (ilosc wezlow
# w kostce, nie pole fizyczne)". Substrat = siec wezlow Gamma (dodatek app:substrat).
#
# ARGUMENT (rygorystyczny, kombinatoryczny - NIE heurystyka):
#   Odwzorowanie stopnia B miedzy triangulowanymi S^3 jest SYMPLICJALNE:
#   kazdy sympleks dziedziny pokrywa co najwyzej JEDEN sympleks celu (raz).
#   Aby pokryc cel B razy, trzeba >= B * (liczba sympleksow celu) sympleksow.
#   MINIMALNA triangulacja S^3 = brzeg 4-sympleksu (dS^4) = 5 czworoscianow.
#   => N_sympleksow(dziedzina) >= 5B.
#
#   Wezly = budzet (staly na objetosc) => objetosc >= V_wezla * 5B => R >= c*B^(1/3).
#
# FALSYFIKATOR (przed liczeniem): jesli minimalna triangulacja S^3 nie jest
#   skonczona / jesli stopien nie jest ograniczony liczba sympleksow => brak
#   kosztu topologicznego => Q NIE jest ladunkiem topologicznym, hipoteza pada.

import itertools, math

print("="*68)
print("[1] MINIMALNA TRIANGULACJA S^3 = brzeg 4-sympleksu (weryfikacja)")
print("-"*68)
V=list(range(5))                                   # 5 wierzcholkow 4-sympleksu
faces={k:list(itertools.combinations(V,k+1)) for k in range(4)}
nV,nE,nF,nT=len(faces[0]),len(faces[1]),len(faces[2]),len(faces[3])
print(f"  wierzcholki={nV}  krawedzie={nE}  trojkaty={nF}  czworosciany={nT}")
chi=nV-nE+nF-nT
print(f"  charakterystyka Eulera chi = {nV}-{nE}+{nF}-{nT} = {chi}   (S^3 ma chi=0)")
print(f"  weryfikacja chi(S^3)=0: {'OK' if chi==0 else 'BLAD'}")
print(f"  => minimalna triangulacja S^3 ma N_min = {nT} czworoscianow.")
print("  >>> SKONCZONA. Falsyfikator NIE zadzialal.")

print("\n"+"="*68)
print("[2] KOSZT TOPOLOGICZNY: ile wezlow/sympleksow potrzeba na nawiniecie B?")
print("-"*68)
Nmin=nT
print(f"  Odwzorowanie symplicjalne stopnia B: kazdy sympleks dziedziny pokrywa")
print(f"  co najwyzej 1 sympleks celu (ze znakiem). Cel ma {Nmin} czworoscianow.")
print(f"  Pokrycie B-krotne wymaga:   N(B) >= {Nmin} * B")
print(f"\n  {'B':>3} {'N_min(B)':>10}  interpretacja")
for B in [1,2,3,4,8,27]:
    tag = "  <- fundamentalny fermion (filar spin-1/2, B=1)" if B==1 else ""
    print(f"  {B:3d} {Nmin*B:10d}{tag}")
print("\n  => KOSZT ROSNIE LINIOWO z ladunkiem topologicznym: N ~ 5B.")
print("     Topologia KOSZTUJE BUDZET. To jest szukany zwiazek Q <-> nawiniecie.")

print("\n"+"="*68)
print("[3] STAD ROZMIAR: R_min(B) przy STALEJ gestosci wezlow (budzet)")
print("-"*68)
print("  gestosc wezlow n = const (budzet)  =>  V >= N(B)/n = 5B/n")
print("  V = (4/3)pi R^3  =>  R_min(B) = (3*5B/(4 pi n))^(1/3)  ~  B^(1/3)")
n=1.0
def Rmin(B): return (3*Nmin*B/(4*math.pi*n))**(1/3)
R1=Rmin(1)
print(f"\n  {'B':>3} {'R_min':>9} {'R_min/R_min(1)':>16} {'B^(1/3)':>10}")
for B in [1,2,3,8,27,64]:
    print(f"  {B:3d} {Rmin(B):9.4f} {Rmin(B)/R1:16.4f} {B**(1/3):10.4f}")
print("\n  >>> R_min(B)/R_min(1) == B^(1/3) DOKLADNIE (skalowanie nasyceniowe).")

print("\n"+"="*68)
print("[4] DOMKNIECIE LUKI")
print("-"*68)
print("  Poprzednio: R_min ~ Q^(1/3) z ZALOZONYM Q = INT(psi-1)dV (luka).")
print("  Teraz:      Q = B (ladunek topologiczny), a koszt N ~ 5B jest")
print("              KOMBINATORYCZNY (bound na stopien odwzorowania symplicjalnego),")
print("              nie zalozony. Budzet = wezly => R_min ~ B^(1/3).")
print("  => 'zachowana ilosc' NIE jest dowolna: to nawiniecie, ktore JUZ mamy")
print("     (pi_3=Z, filar spin-1/2 B=1) i ktore JUZ jest chronione topologicznie.")

print("\n"+"="*68)
print("OGRANICZENIA (uczciwie):")
print("  - bound N>=5B dotyczy odwzorowan SYMPLICJALNYCH; substrat niekoniecznie")
print("    jest symplicjalny => to IDEALIZACJA. Skalowanie (~B) jest odporne,")
print("    STALA (5) zalezy od dyskretyzacji.")
print("  - to daje DOLNE ograniczenie rozmiaru (anty-kolaps), nie pelny profil.")
print("  - gestosc wezlow n przyjeta stala; jej zwiazek z psi (n ~ 1/h) NIE")
print("    zostal tu uzyty - pelne sprzezenie budzet<->profil pozostaje OTWARTE.")
