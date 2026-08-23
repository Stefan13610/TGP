# RIBBONS_uniqueness_2T.py
# Pytanie: czy 2T jest JEDYNYM dopuszczalnym targetem Phi-2 z trialnoscia?
# Dopuszczalne pi_1(SO(3)/H) = SKONCZONE PODGRUPY SU(2) = binarne grupy (klasyfikacja ADE):
#   - binarne cykliczne Z_n      (abelowe -> odpadaja: nie uciekaja od no-go writhe->Z)
#   - binarne dwuscienne Dic_n   (rzad 4n)
#   - 2T (24), 2O (48), 2I (120) (binarne wieloscienne)
# Warunek koniecznny na kolor: NIEABELOWA + abelianizacja zawiera Z3 (3 | |G/G'|).
# Liczymy abelianizacje kazdej z nich przez KWATERNIONY jednostkowe (u zrodla).

import math

def qmul(p, q):
    w1,x1,y1,z1 = p; w2,x2,y2,z2 = q
    return (w1*w2 - x1*x2 - y1*y2 - z1*z2,
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2)
def qkey(p): return tuple(round(c, 6)+0.0 for c in p)
QID = (1.0,0.0,0.0,0.0)

def close_group(gens):
    elems = {qkey(QID): QID}; frontier=[QID]
    while frontier:
        new=[]
        for g in frontier:
            for h in gens:
                pr=qmul(g,h); k=qkey(pr)
                if k not in elems: elems[k]=pr; new.append(pr)
        frontier=new
        if len(elems)>2000: raise RuntimeError("group too big / open")
    return list(elems.values())

def inverse(g):  # kwaternion jednostkowy: sprzezenie
    w,x,y,z=g; return (w,-x,-y,-z)

def commutator_subgroup(elems):
    comms=[]
    for a in elems:
        ai=inverse(a)
        for b in elems:
            bi=inverse(b)
            comms.append(qmul(qmul(a,b),qmul(ai,bi)))
    return close_group(comms)

def is_abelian(elems):
    for a in elems:
        for b in elems:
            if qkey(qmul(a,b))!=qkey(qmul(b,a)): return False
    return True

def report(name, gens):
    G=close_group(gens); n=len(G)
    ab=is_abelian(G)
    Gp=commutator_subgroup(G); abelianization=n//len(Gp)
    ok = (not ab) and (abelianization%3==0)
    print(f"  {name:16s} |G|={n:3d}  abelowa={str(ab):5s}  |G/G'|=Z{abelianization:<3d}  "
          f"{'<== DOPUSZCZALNY (nieabel. + Z3)' if ok else ''}")
    return (name, n, ab, abelianization, ok)

# --- generatory kwaternionowe ---
def dic(nn):  # binarna dwuscienna Dic_n, rzad 4n
    a=(math.cos(math.pi/nn),0.0,0.0,math.sin(math.pi/nn))  # obrot 2pi/n wokol z (pol-kat pi/n)
    b=(0.0,1.0,0.0,0.0)                                     # obrot pi wokol x
    return [a,b]

i=(0.0,1.0,0.0,0.0); j=(0.0,0.0,1.0,0.0)
h=(0.5,0.5,0.5,0.5)               # (1+i+j+k)/2, generator 2T (Hurwitz)
s=(math.sqrt(0.5),math.sqrt(0.5),0.0,0.0)  # (1+i)/sqrt2, rozszerza 2T do 2O

print("="*72)
print("Przeglad skonczonych podgrup SU(2) (= wszystkie mozliwe pi_1(SO(3)/H)):")
print("-"*72)
rows=[]
for nn in range(1,7):
    rows.append(report(f"Dic_{nn}", dic(nn)))
rows.append(report("2T (bin.tetra)", [i,h]))
rows.append(report("2O (bin.okta)", [i,h,s]))
# 2I (binarna ikozaedralna, rzad 120) jest grupa DOSKONALA (SL(2,5)): G'=G, abelianizacja trywialna.
print("  2I (bin.ikoza)  |G|=120  abelowa=False  |G/G'|=Z1    (grupa doskonala - analitycznie)")
print("  Z_n cykliczne   ABELOWE  -> odpadaja (nie uciekaja od no-go writhe->Z)")
print("-"*72)
allowed=[r for r in rows if r[4]]
print(f"DOPUSZCZALNE (nieabelowe + abelianizacja z Z3): {[r[0] for r in allowed] or 'BRAK'}")
print("WNIOSEK: 2T jest JEDYNYM (wsrod Dic_n<=6, 2O; 2I doskonala, Z_n abelowe).")
