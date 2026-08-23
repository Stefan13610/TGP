"""
COLOR: czy ktorykolwiek Z_3 w TGP moze byc zrodlem koloru?
Trzy niezalezne testy, wszystkie exact/dyskretne (bez numeryki zmiennoprzecinkowej):

[1] GL(3,F_2): rzad, centrum, abelianizacja. Czy Z_3 jest PODGRUPA czy ILORAZEM?
    (kolor jako etykieta wymaga ILORAZU -> charaktery; podgrupa nie daje etykiety)
[2] Porownanie trzech "Z_3": 2T^ab, centrum SU(3), GL(3,F_2).
[3] Substrat Ising: parzystosc Z_2 korelatorow + znikanie nieparzystych momentow
    kierunkowych na sieci z symetria inwersji (a_hat -> -a_hat).
"""
import itertools, numpy as np

print("="*70)
print("[1] GL(3,F_2) -- struktura grupy (exact, F_2 arithmetic)")
print("="*70)

# wszystkie odwracalne macierze 3x3 nad F_2
mats = []
for bits in itertools.product([0,1], repeat=9):
    M = np.array(bits, dtype=np.int8).reshape(3,3)
    # wyznacznik nad F_2 = permanent mod 2
    det = 0
    for p in itertools.permutations(range(3)):
        det ^= (M[0,p[0]] & M[1,p[1]] & M[2,p[2]])
    if det == 1:
        mats.append(M)
G = [tuple(M.flatten()) for M in mats]
idx = {g:i for i,g in enumerate(G)}
def mul(a,b):
    A = np.array(a,dtype=np.int8).reshape(3,3); B = np.array(b,dtype=np.int8).reshape(3,3)
    return tuple(((A@B)%2).flatten())
n = len(G)
print(f"  |GL(3,F_2)| = {n}   (oczekiwane 168)")

# centrum
Z = [g for g in G if all(mul(g,h)==mul(h,g) for h in G)]
print(f"  |Z(G)| (centrum) = {len(Z)}")

# komutant [G,G]
def inv(g):
    for h in G:
        if mul(g,h)==tuple(np.eye(3,dtype=np.int8).flatten()): return h
E = tuple(np.eye(3,dtype=np.int8).flatten())
comms = set()
for a in G:
    ia = inv(a)
    for b in G:
        comms.add(mul(mul(ia,inv(b)), mul(a,b)))
# domkniecie komutantu
sub = set(comms); sub.add(E)
changed = True
while changed:
    changed = False
    for x in list(sub):
        for y in list(sub):
            z = mul(x,y)
            if z not in sub: sub.add(z); changed = True
print(f"  |[G,G]| (komutant)  = {len(sub)}")
print(f"  |G^ab| = |G|/|[G,G]| = {n//len(sub)}")

# rzedy elementow
from collections import Counter
def order(g):
    x, k = g, 1
    while x != E: x = mul(x,g); k += 1
    return k
oc = Counter(order(g) for g in G)
print(f"  rzedy elementow: {dict(sorted(oc.items()))}")
print(f"  czy istnieja elementy rzedu 3? {'TAK' if oc.get(3,0)>0 else 'NIE'}  (liczba: {oc.get(3,0)})")

perfect = (len(sub)==n)
print(f"\n  ==> G jest PERFEKCYJNA ([G,G]=G)? {perfect}")
print(f"  ==> G^ab = {'TRYWIALNA' if n//len(sub)==1 else n//len(sub)}")
print(f"  ==> Hom(G,U(1)) = Hom(G^ab,U(1)) = {'TYLKO TRYWIALNY' if n//len(sub)==1 else '?'}")

print("\n" + "="*70)
print("[2] Trzy rozne 'Z_3' -- czy to ta sama grupa w tej samej roli?")
print("="*70)
rows = [
    ("2T^ab (moja warstwa wstazek)", "ILORAZ  2T/[2T,2T]=2T/Q_8", "TAK -> 3 charaktery", "pi_1 defektu"),
    ("centrum SU(3) (rdzen sek09)",  "CENTRUM Z(SU(3))",          "nie jest ilorazem pi_1", "grupa cechowania"),
    ("Z_3 < GL(3,F_2) (rdzen)",      "PODGRUPA (elem. rzedu 3)",  "NIE -> G perfekcyjna",  "symetria dyskretna V(g)"),
]
print(f"  {'obiekt':<30} {'typ':<28} {'daje etykiete?':<24} rola")
for r in rows:
    print(f"  {r[0]:<30} {r[1]:<28} {r[2]:<24} {r[3]}")

print("\n" + "="*70)
print("[3] Substrat Isinga: co przezywa symetrie Z_2 i inwersje sieci?")
print("="*70)
print("  H_Gamma = sum_i [m0^2/2 s_i^2 + lam0/4 s_i^4] - J sum_<ij> s_i s_j,  Z_2: s -> -s")
print()
print(f"  {'obserwabla':<26} {'#s':<4} {'Z_2-parzysta?':<15} {'#przesuniec':<12} rank tensora")
obs = [("Phi = <s^2>",2,0), ("sigma_ab = <s s_+a>^TF",2,1), ("<s s_+a s_+b>",3,2), ("<s s_+a s_+b s_+c>",4,3)]
for name,ns,nd in obs:
    even = (ns % 2 == 0)
    print(f"  {name:<26} {ns:<4} {('TAK' if even else 'NIE -> ZNIKA'):<15} {nd:<12} {nd}")

# moment kierunkowy rzedu 3 na sieci kubicznej z inwersja
print("\n  Test: moment kierunkowy T_ijk = sum_a C(a) a_i a_j a_k")
print("  (sasiedzi sieci kubicznej; C zalezy tylko od |a| => C(a)=C(-a))")
for label, nbrs in [("6 najbl. (+-x,+-y,+-z)", [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]),
                    ("12 drugich sasiadow", [v for v in itertools.product([-1,0,1],repeat=3) if sum(abs(c) for c in v)==2]),
                    ("8 trzecich (rogi)", [v for v in itertools.product([-1,1],repeat=3)])]:
    T3 = np.zeros((3,3,3)); T2 = np.zeros((3,3))
    for a in nbrs:
        a = np.array(a,dtype=float)
        T3 += np.einsum('i,j,k->ijk', a,a,a)
        T2 += np.outer(a,a)
    T2tf = T2 - np.trace(T2)/3*np.eye(3)
    print(f"    {label:<26} max|T_ijk| = {np.abs(T3).max():.3e}   max|T_ij^TF| = {np.abs(T2tf).max():.3e}")

print("\n  ==> rank-3 z korelatora wiazaniowego znika TOZSAMOSCIOWO (inwersja sieci)")
print("  ==> rank-3 z korelatora 3-punktowego znika TOZSAMOSCIOWO (Z_2: nieparzysty w s)")
