# BUDGET_scale_cap.py
# HIPOTEZA: skala obiektu pochodzi z WIEZU BUDZETOWEGO (sek08c def:info-budget),
# a NIE z minimum energii (energia = relacyjna, nie prymityw substratu).
#
# ZRODLO (core/sek08c, cytat):
#   B = N_B * s_0,  s_0 = m_sp^2 = gamma (entropia na wezel)
#   "B jest stala - nie zalezy od Phi(x), poniewaz N_B jest liczba TOPOLOGICZNA"
#   budzet dzielony: n_sp + n_czas = B          (wiez ADDYTYWNY, prop:antipodal proof)
#   h = n_sp/n_0,  f = n_czas/n_0t              (definicje udzialow)
#   n_sp * n_czas = const  =>  f*h = 1          (zalozenie MULTIPLIKATYWNE, tamze)
#
# UWAGA UCZCIWOSCI: rdzen sam zaznacza, ze wiez addytywny NIE daje fh=1 i ze
# multiplikatywnosc jest DODATKOWYM zalozeniem. Ja uzywam OBU naraz:
#   fh = 1  (ksztalt metryki)  ORAZ  n_sp + n_czas <= B  (budzet nieprzekraczalny).
# To jest MOJA INTERPRETACJA, nie cytat z twierdzenia. Jawnie oznaczone.
#
# FALSYFIKATOR (deklarowany PRZED liczeniem):
#   Jesli n_sp + n_czas jest NIEOGRANICZONE od gory przy fh=1 (brak cap na psi),
#   to budzet NIE daje skali => hipoteza pada.
#   Jesli cap istnieje => sprawdzamy, czy blokuje kolaps Derricka.

import numpy as np

print("="*68)
print("[1] CZY BUDZET DAJE CAP? Funkcja zuzycia budzetu g(h)=n0*h + n0t/h")
print("-"*68)
print("  z fh=1: n_sp = n0*h,  n_czas = n0t*f = n0t/h")
print("  zuzycie calkowite: g(h) = n0*h + n0t/h")
n0, n0t = 1.0, 1.0     # normalizacja referencyjna (h=f=1 przy psi=1)
hs=np.linspace(0.01,50,500000)
g=n0*hs+n0t/hs
imin=np.argmin(g)
print(f"  g(h) ma MINIMUM w h*={hs[imin]:.4f}, g_min={g[imin]:.4f}")
print(f"  g(h) -> inf gdy h->0 ORAZ h->inf   => zuzycie ROSNIE w obie strony")
print(f"  => dla budzetu B > g_min istnieje SKONCZONY przedzial [h_min,h_max]")
print("  >>> CAP ISTNIEJE. Falsyfikator NIE zadzialal.")

def h_range(B):
    # n0*h + n0t/h <= B  =>  n0 h^2 - B h + n0t <= 0
    disc=B**2-4*n0*n0t
    if disc<0: return None
    return ((B-np.sqrt(disc))/(2*n0), (B+np.sqrt(disc))/(2*n0))

print("\n  Dozwolony przedzial h dla roznych budzetow B (jednostki g_min=2):")
print(f"  {'B':>6} {'h_min':>10} {'h_max':>10}")
for B in [2.0,2.5,3.0,4.0,6.0,10.0]:
    r=h_range(B)
    print(f"  {B:6.2f} {r[0]:10.4f} {r[1]:10.4f}" if r else f"  {B:6.2f}   brak (B<g_min)")

print("\n"+"="*68)
print("[2] PRZELOZENIE NA psi (cap na amplitude pola)")
print("-"*68)
print("  Forma I (g_ij=psi*delta):        h=psi        => psi_max = h_max")
print("  M9.1'' (kanoniczna):             h=psi/(4-3psi) => psi = 4h/(1+3h)")
def psi_from_h(h): return 4*h/(1+3*h)
print(f"  {'B':>6} {'psi_max(formaI)':>16} {'psi_max(M9.1\"\")':>16}")
for B in [2.5,3.0,4.0,6.0,10.0,100.0]:
    r=h_range(B)
    print(f"  {B:6.2f} {r[1]:16.4f} {psi_from_h(r[1]):16.4f}")
print("  UWAGA: w M9.1'' psi_max -> 4/3 = 1.3333 gdy B->inf (horyzont f=0).")
print("  Czyli budzet daje cap STRICTLY PONIZEJ horyzontu - cap jest realny,")
print("  a psi=4/3 to granica absolutna (caly budzet na przestrzen, zegar staje).")
print("  [NIE uzywam tu obalonego mostu psi<->g0 (fit 2-punktowy) - to inny argument.]")

print("\n"+"="*68)
print("[3] CZY CAP BLOKUJE KOLAPS DERRICKA? (rachunek skalowania)")
print("-"*68)
print("  Obiekt: psi(r) = 1 + A*chi(r/R), amplituda A, rozmiar R.")
print("  Zachowana 'ilosc substratu' Q ~ INT (psi-1) dV = c*A*R^3  (staly ksztalt chi)")
print("  Kolaps Derricka: R->0. Przy STALYM Q wymaga A = Q/(c R^3) -> INF.")
print("  ALE budzet narzuca A <= A_max = psi_max - 1.  Stad:")
print("     Q/(c R^3) <= A_max   =>   R >= R_min = (Q/(c*A_max))^(1/3)")
c=1.0
print(f"\n  {'psi_max':>9} {'A_max':>8} {'R_min (Q=1)':>13} {'R_min (Q=8)':>13}")
for pm in [1.05,1.10,1.20,1.30,4/3]:
    Am=pm-1.0
    print(f"  {pm:9.4f} {Am:8.4f} {(1.0/(c*Am))**(1/3):13.4f} {(8.0/(c*Am))**(1/3):13.4f}")
print("\n  >>> KOLAPS ZABLOKOWANY: istnieje SKONCZONY minimalny rozmiar R_min > 0.")
print("      Mechanizm: cap na amplitude + zachowana ilosc => nie da sie scisnac.")
print("      ZERO minimalizacji energii. ZERO czlonu Skyrme'a.")

print("\n"+"="*68)
print("[4] STRUKTURALNA PREDYKCJA (falsyfikowalna)")
print("-"*68)
print("  R_min ~ Q^(1/3)  <=>  GESTOSC NASYCENIA (stala max gestosc).")
print("  Q=1 -> R=1 ; Q=8 -> R=2 ; Q=27 -> R=3  (przy A=A_max):")
for Q in [1,8,27,64]:
    print(f"     Q={Q:3d}  R_min={(Q/(c*(4/3-1)))**(1/3)/ (1/(c*(4/3-1)))**(1/3):6.3f} (znormalizowane)")
print("  To jest podpis materii NASYCONEJ (jak R ~ A^(1/3) dla jader).")
print("  Uwaga: samo skalowanie R~Q^(1/3) jest generyczne dla kazdego capu;")
print("  NATYWNE jest tu ZRODLO capu (budzet), nie sam wykladnik.")

print("\n"+"="*68)
print("BILANS:")
print("  CAP na psi: WYPADA z budzetu (n_sp+n_czas<=B przy fh=1) - falsyfikator nie zadzialal.")
print("  KOLAPS: zablokowany capem + zachowana iloscia => R_min > 0 (dolne ograniczenie).")
print("  ROZPAD/rozplyniecie (gorne ogr.): NIE z budzetu - potrzebna PRZESZKODA")
print("     TOPOLOGICZNA (sektor wstazek 2T: wstazka nie domyka sie sama) - MAMY.")
print("  => rozmiar OBUSTRONNIE ograniczony: budzet od dolu, topologia od gory.")
