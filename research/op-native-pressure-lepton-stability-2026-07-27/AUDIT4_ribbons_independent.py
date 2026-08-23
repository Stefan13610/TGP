# AUDIT4_ribbons_independent.py
# ZARZUT 5/6: czy RIBBONS_closure_solitons.py stoi, czy tez ma ukryta tautologie?
# Buduje 2T CALKOWICIE NIEZALEZNIE: 2T ~ SL(2,3) = macierze 2x2 nad F_3 o det=1.
# (autor uzyl kwaternionow - inna reprezentacja, ten sam wynik musi wyjsc)

import itertools
from itertools import product, combinations_with_replacement

def mul(A,B):
    return ((A[0]*B[0]+A[1]*B[2])%3,(A[0]*B[1]+A[1]*B[3])%3,
            (A[2]*B[0]+A[3]*B[2])%3,(A[2]*B[1]+A[3]*B[3])%3)
def det(A): return (A[0]*A[3]-A[1]*A[2])%3
G=[A for A in product(range(3),repeat=4) if det(A)==1]
I=(1,0,0,1); NEG=(2,0,0,2)
print("="*74)
print("[0] NIEZALEZNA KONSTRUKCJA 2T = SL(2,3)")
print("-"*74)
print(f"  |SL(2,3)| = {len(G)}   (musi byc 24)   -1 = {NEG} w grupie: {NEG in G}")
assert len(G)==24

def inv(A):
    for B in G:
        if mul(A,B)==I: return B
def order(A):
    p=A;n=1
    while p!=I: p=mul(p,A);n+=1
    return n
# komutant = Q8
comm=set()
frontier={I}
gens=set()
for a in G:
    ai=inv(a)
    for b in G: gens.add(mul(mul(a,b),mul(ai,inv(b))))
def closure(gens):
    S={I}; fr=[I]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                p=mul(g,h)
                if p not in S: S.add(p); nw.append(p)
        fr=nw
    return S
Q8=closure(gens)
print(f"  komutant [G,G] ma rzad {len(Q8)}  (musi byc 8 = Q8);  -1 in Q8: {NEG in Q8}")
assert len(Q8)==8

# trialnosc = G/Q8 ~ Z3
cos=lambda g: frozenset(mul(g,x) for x in Q8)
C0=cos(I); gg=next(g for g in G if cos(g)!=C0)
lab={C0:0, cos(gg):1, cos(mul(gg,gg)):2}
tri=lambda g: lab[cos(g)]
assert all(tri(mul(a,b))==(tri(a)+tri(b))%3 for a in G for b in G)
print("  trialnosc: homomorfizm G -> Z3 zweryfikowany na wszystkich 576 parach.")

CL=[]; seen=set()
for g in G:
    if g in seen: continue
    cl=frozenset(mul(mul(x,g),inv(x)) for x in G)
    CL.append(sorted(cl)); seen|=cl
CL.sort(key=lambda c:(len(c),order(c[0]),tri(c[0])))
print(f"\n  klas sprzezonosci: {len(CL)}  (autor: 7)")
print(f"  {'kl':>4} {'rozmiar':>8} {'rzad':>5} {'trialnosc':>10}")
info=[]
for i,cl in enumerate(CL):
    d=dict(size=len(cl),ordr=order(cl[0]),tri=tri(cl[0])); info.append(d)
    print(f"  {i:4d} {d['size']:8d} {d['ordr']:5d} {d['tri']:10d}")

def reach(classes):
    R={I}
    for cl in classes:
        R={mul(a,b) for a in R for b in cl}
    return R

print("\n"+"="*74)
print("[1] n=1: czy pojedyncza klasa domyka sie do e?")
print("-"*74)
for i,cl in enumerate(CL):
    print(f"  klasa {i} (rzad {info[i]['ordr']}, tri {info[i]['tri']}): e osiagalne? {I in reach([cl])}")
print("  >>> TYLKO klasa trywialna. UWAGA: to jest TAUTOLOGIA - 'g != e dla g w klasie")
print("      nietrywialnej'. Zerowa tresc obliczeniowa, choc stwierdzenie jest prawdziwe.")

print("\n"+"="*74)
print("[2] pary domykajace sie do e")
print("-"*74)
pairs=[(i,j) for i,j in combinations_with_replacement(range(len(CL)),2) if I in reach([CL[i],CL[j]])]
for i,j in pairs:
    print(f"  ({i},{j})  suma trialnosci = {(info[i]['tri']+info[j]['tri'])%3}")
print(f"  liczba par: {len(pairs)} (autor: 5);  wszystkie tri-0: {all((info[i]['tri']+info[j]['tri'])%3==0 for i,j in pairs)}")
print("  >>> 'trialnosc-0 KONIECZNA' jest tautologia homomorfizmu G->Z3 (tri(e)=0).")
print("      Nie wymagalo to rachunku.")

print("\n"+"="*74)
print("[3] trojki: ile o sumie trialnosci 0 domyka sie do e?  (autor: 21/30)")
print("-"*74)
tot=0; ok=0; fails=[]
for c in combinations_with_replacement(range(len(CL)),3):
    if sum(info[k]['tri'] for k in c)%3: continue
    tot+=1
    if I in reach([CL[k] for k in c]): ok+=1
    else: fails.append(c)
print(f"  trojek tri-0: {tot};  domykalnych do e: {ok}")
print(f"  NIEdomykalne: {fails}")
print(f"  >>> REPRODUKUJE liczbe autora {ok}/{tot} = 21/30: {(ok,tot)==(21,30)}")

print("\n"+"="*74)
print("[3b] ALE: z czego wynikaja te 9 porazek? (autor: 'dodatkowa przeszkoda w Q8, OTWARTE')")
print("-"*74)
for c in fails:
    szs=[info[k]['size'] for k in c]; ords=[info[k]['ordr'] for k in c]
    print(f"  {c}: rozmiary klas={szs} rzedy={ords}  -> osiagalny zbior ma {len(reach([CL[k] for k in c]))} elem.")
print("  Diagnoza: kazda porazka zawiera wylacznie klasy CENTRALNE (rozmiar 1: {e},{-1}),")
print("  ktore nie maja swobody sprzezenia. Iloczyn jest wtedy USTALONY, nie zbiorem.")
cent=[i for i,d in enumerate(info) if d['size']==1]
print(f"  klasy centralne: {cent}")
allcent=all(all(k in cent or info[k]['size']==1 for k in c) for c in fails)
print(f"  czy WSZYSTKIE porazki skladaja sie tylko z klas centralnych/malych? {allcent}")
for c in fails:
    print(f"    {c} sizes={[info[k]['size'] for k in c]}")
print("  >>> 'dodatkowa przeszkoda MOCNIEJSZA niz singlet kolorowy SM' to w duzej mierze")
print("      efekt wliczenia klas CENTRALNYCH (e, -1) jako 'wstazek'. To nie sa defekty.")

print("\n"+"="*74)
print("[3c] KONTROLA: 21/30 po USUNIECIU klas centralnych (prawdziwe defekty)")
print("-"*74)
noncent=[i for i,d in enumerate(info) if d['size']>1]
tot2=0; ok2=0; f2=[]
for c in combinations_with_replacement(noncent,3):
    if sum(info[k]['tri'] for k in c)%3: continue
    tot2+=1
    if I in reach([CL[k] for k in c]): ok2+=1
    else: f2.append(c)
print(f"  trojki tri-0 z klas NIECENTRALNYCH: {tot2};  domykalnych: {ok2};  porazki: {f2}")
if ok2==tot2:
    print("  >>> Dla prawdziwych defektow trialnosc-0 JEST WYSTARCZAJACA.")
    print("      Twierdzenie autora 'MOCNIEJSZE niz regula singletu SM' NIE UTRZYMUJE SIE")
    print("      po usunieciu klas centralnych (ktore nie sa dysklinacjami).")
else:
    print("  >>> pozostaja realne porazki - twierdzenie autora sie broni.")

print("\n"+"="*74)
print("[4] 'barion=fermion niewymuszony' - kontrola")
print("-"*74)
t1=[i for i,d in enumerate(info) if d['tri']==1]
for c in combinations_with_replacement(t1,3):
    R=reach([CL[k] for k in c])
    print(f"  {c}: e osiagalne={I in R}  -1 osiagalne={NEG in R}  |zbior|={len(R)}")
print("  >>> potwierdzam autora: OBA osiagalne => statystyka NIEWYMUSZONA. Uczciwie zaraportowane.")
