# AUDIT4_simplicial.py
# ZARZUT 4: BUDGET_topological_cost.py - bound N >= 5B.
#  (a) czy "weryfikacja" minimalnej triangulacji S^3 przez chi=0 cokolwiek weryfikuje?
#  (b) czy bound jest poprawnym twierdzeniem?
#  (c) czy bound ogranicza LICZBE WEZLOW (=budzet wg sek08c), czy liczbe CZWOROSCIANOW?

import itertools, math
from math import comb

print("="*76)
print("[a] CZY chi = 0 WERYFIKUJE COKOLWIEK? (autor: 'chi(S^3)=0: OK')")
print("-"*76)
print("  Fakt: dla KAZDEJ zamknietej rozmaitosci 3-wym. chi = 0 (dualnosc Poincare).")
print("  Wiec chi=0 NIE odroznia S^3 od T^3, RP^3, ani od niczego innego,")
print("  i tym bardziej NIE dowodzi MINIMALNOSCI triangulacji.")
print("\n  Dowod operacyjny: ruch stellarny 1->4 (podzial czworoscianu) zmienia")
print("  (V,E,F,T) o (+1,+4,+6,+3) => d(chi) = 1-4+6-3 = 0. chi jest slepe.")
V,E,F,T = 5,10,10,5
print(f"  {'krok':>5} {'V':>5} {'E':>5} {'F':>5} {'T':>5} {'chi':>5}   (wszystko to S^3, rozne triangulacje)")
for k in range(6):
    print(f"  {k:5d} {V:5d} {E:5d} {F:5d} {T:5d} {V-E+F-T:5d}")
    V,E,F,T = V+1,E+4,F+6,T+3
print("  >>> chi=0 dla WSZYSTKICH. 'Weryfikacja' w [1] skryptu autora jest PUSTA.")
print("      (Sama teza 'min. triangulacja S^3 = dDelta^4 = 5 czworoscianow' jest PRAWDZIWA,")
print("       ale to znany fakt z literatury, a NIE cos, co ten skrypt sprawdzil.)")

print("\n"+"="*76)
print("[b] CZY BOUND T_dom >= 5|B| JEST POPRAWNY? (audyt logiki, nie numeryki)")
print("-"*76)
print("  Odwzorowanie symplicjalne f: K -> L. Dla generycznego punktu p wewnatrz")
print("  czworoscianu sigma z L:  deg f = SUMA znakow po czworoscianach tau z K,")
print("  ktore odwzorowuja sie NIEZDEGENEROWANIE NA sigma.")
print("  => #{tau : f(tau)=sigma} >= |deg f| = |B|   dla KAZDEGO sigma z L.")
print("  Kazde tau odwzorowuje sie na CO NAJWYZEJ jeden sigma (albo degeneruje).")
print("  => T_dom >= |B| * T_cel. Biorac NAJMNIEJSZY mozliwy cel T_cel=5 (dDelta^4)")
print("     dostajemy NAJSLABSZY, czyli BEZPIECZNY bound: T_dom >= 5|B|.")
print("  >>> TWIERDZENIE POPRAWNE. Uzycie '5' (minimum po CELU) jest w dobra strone.")

print("\n"+"="*76)
print("[c] ALE: BOUND DOTYCZY CZWOROSCIANOW, A BUDZET LICZY WEZLY")
print("-"*76)
print("  sek08c (cytowane przez autora): 'N_B = ilosc WEZLOW w kostce'.")
print("  Skrypt autora robi podstawienie:  'Wezly = budzet ... V >= V_wezla * 5B'")
print("  czyli utozsamia LICZBE CZWOROSCIANOW z LICZBA WEZLOW. To NIE jest to samo.")
print("\n  Twarde ograniczenie na wierzcholki: T <= C(V,4)  (kazdy czworoscian to 4-elem. podzbior)")
print(f"  {'B':>4} {'T_min=5B':>9} {'V_min (T<=C(V,4))':>19} {'R~V^(1/3), znorm.':>19} {'B^(1/3)':>9} {'B^(1/12)':>10}")
def Vmin_for_T(T):
    v=4
    while comb(v,4) < T: v+=1
    return v
base=None
for B in [1,2,3,8,27,64,1000,10**6]:
    Tm=5*B; vm=Vmin_for_T(Tm); r=vm**(1/3)
    if base is None: base=r
    print(f"  {B:4d} {Tm:9d} {vm:19d} {r/base:19.4f} {B**(1/3):9.4f} {B**(1/12):10.4f}")
print("  >>> W NAJGORSZYM PRZYPADKU liczba WEZLOW rosnie jak ~B^(1/4), wiec R ~ B^(1/12),")
print("      a NIE B^(1/3). Wykladnik 1/3 wymaga DODATKOWEGO zalozenia T = O(V)")
print("      (triangulacja o ograniczonym stopniu wierzcholka), ktorego nigdzie nie ma.")

print("\n"+"="*76)
print("[c2] Ile realnie wynosi T/V w triangulacjach S^3? (relacja Eulera + 2F=4T)")
print("-"*76)
print("  Zamknieta rozmaitosc 3-wym.: kazdy trojkat lezy w dokladnie 2 czworoscianach => F=2T.")
print("  chi=0 => V-E+F-T=0 => E = V + T.  Ale E <= C(V,2), wiec T <= C(V,2)-V.")
V0,E0,F0,T0=5,10,10,5
print(f"  {'V':>5} {'T':>6} {'T/V':>7} {'E=V+T':>7} {'C(V,2)':>8}   (ciag 1->4 od dDelta^4)")
V,E,F,T=5,10,10,5
for k in range(8):
    print(f"  {V:5d} {T:6d} {T/V:7.3f} {V+T:7d} {comb(V,2):8d}")
    V,E,F,T=V+1,E+4,F+6,T+3
print("  >>> T/V -> 3 w tym ciagu, ale E=V+T <= C(V,2) dopuszcza T ~ V^2/2.")
print("      Zwiazek T ~ V jest ZALOZENIEM o regularnosci siatki, nie twierdzeniem.")

print("\n"+"="*76)
print("[d] CZY FALSYFIKATOR MOGL ZADZIALAC?")
print("-"*76)
print("  Autor: 'jesli minimalna triangulacja S^3 NIE JEST SKONCZONA ... hipoteza pada'.")
print("  Kazda zwarta rozmaitosc PL ma skonczona triangulacje (fakt podrecznikowy).")
print("  => falsyfikator mial prawdopodobienstwo zadzialania 0. To NIE jest falsyfikator.")
print("  Analogicznie w BUDGET_scale_cap: 'jesli h+1/h jest NIEOGRANICZONE z gory'.")
print("  h+1/h -> inf przy h->0 i h->inf to AM-GM/rachunek elementarny.")
print("  => tez prawdopodobienstwo zadzialania 0.")
print("  >>> OBA 'falsyfikatory zadeklarowane przed liczeniem' sa PSEUDO-falsyfikatorami:")
print("      testowaly fakty matematyczne znane a priori, nie hipotezy o TGP.")

print("\n"+"="*76)
print("[e] REPRODUKCJA POZOSTALYCH LICZB AUTORA (kontrola arytmetyki)")
print("-"*76)
Vv=list(range(5))
faces={k:list(itertools.combinations(Vv,k+1)) for k in range(4)}
print(f"  dDelta^4: V={len(faces[0])} E={len(faces[1])} F={len(faces[2])} T={len(faces[3])}  chi={len(faces[0])-len(faces[1])+len(faces[2])-len(faces[3])}")
n=1.0; Nmin=5
Rmin=lambda B:(3*Nmin*B/(4*math.pi*n))**(1/3)
print(f"  R_min(B)/R_min(1) dla B=1,2,3,8,27,64: "
      + ", ".join(f"{Rmin(B)/Rmin(1):.4f}" for B in [1,2,3,8,27,64]))
print("  >>> arytmetyka autora sie zgadza; R/R(1)=B^(1/3) jest TOZSAMOSCIA (V ~ B, R ~ V^(1/3)).")
print("      Zero tresci fizycznej: to definicja objetosci, nie predykcja.")
