# RIBBONS_s54_frame_structure.py
# Domkniecie §5.4 sciezka (A): kolor 2T z reperu samo-generowanej metryki.
# M = S^3/2T = przestrzen konfiguracji REPERU (frame) metryki g_ij=psi*delta_ij
# z tetraedryczna symetria substratu T. Wtedy:
#   pi_3(M)=Z  -> nawiniecie B (SPIN, jak filar),   pi_1(M)=2T -> dysklinacje (KOLOR).
# Oba z JEDNEGO M -> os §0.
#
# Nowy fakt do weryfikacji (dotyka §0 = czy spin i kolor sa niezalezne):
#   jak centralne -1 (spin-framing) i klasy rzedu-3 (kolor) siedza w 2T?
#   Twierdzenie: 2T = Q8 |x| Z3 (Q8 normalny, Z3 iloraz=abelianizacja=kolor),
#   centrum Z2={+-1} lezy w Q8 => -1 (spin) jest "bezbarwny" (trywialny w Z3),
#   a elementy rzedu 3 (kolor) NIE sa centralne. => spin i kolor to ROZNE
#   skladniki jednej grupy, powiazane struktura poluproduktu (hook do §0 lock).

import math
def qmul(p,q):
    w1,x1,y1,z1=p; w2,x2,y2,z2=q
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
def qk(p): return tuple(round(c,6)+0.0 for c in p)
QID=(1.,0.,0.,0.); NEG=(-1.,0.,0.,0.)
def close(gens):
    E={qk(QID):QID}; fr=[QID]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                pr=qmul(g,h);k=qk(pr)
                if k not in E:E[k]=pr;nw.append(pr)
        fr=nw
    return list(E.values())
def inv(g):w,x,y,z=g;return(w,-x,-y,-z)
def order(g):
    p=g;n=1
    while qk(p)!=qk(QID):p=qmul(p,g);n+=1
    return n
def comm_sub(E):
    cs=[]
    for a in E:
        ai=inv(a)
        for b in E:
            cs.append(qmul(qmul(a,b),qmul(ai,inv(b))))
    return close(cs)

G=close([(0.,1.,0.,0.),(0.5,0.5,0.5,0.5)])  # 2T
assert len(G)==24
Q8=comm_sub(G)
Q8k=set(qk(x) for x in G if qk(x) in set(qk(y) for y in Q8))
Q8set=set(qk(x) for x in Q8)
center=[g for g in G if all(qk(qmul(g,h))==qk(qmul(h,g)) for h in G)]

print("2T struktura (weryfikacja u zrodla):")
print(f"  |2T|={len(G)},  komutator |Q8|={len(Q8)} (oczek. 8),  centrum |Z|={len(center)}")
print(f"  -1 w Q8 (komutatorze)?  {qk(NEG) in Q8set}   => spin-framing lezy w bezbarwnym rdzeniu")
# kolor = obraz w abelianizacji 2T/Q8 = Z3. Warstwa elementu g: g*Q8.
def coset(g): return frozenset(qk(qmul(g,x)) for x in Q8)
cos={}
for g in G: cos.setdefault(coset(g),[]).append(g)
print(f"  liczba warstw 2T/Q8 = {len(cos)} (oczek. 3 = Z3 = kolor)")
# kolor elementu: ktora z 3 warstw. -1 i elementy rzedu 3:
def color_class(g):
    ck=coset(g)
    return list(cos.keys()).index(ck)
print(f"  kolor(-1) = warstwa {color_class(NEG)} (trywialna=0 => -1 BEZBARWNY)")
ord3=[g for g in G if order(g)==3]
colors3=sorted(set(color_class(g) for g in ord3))
print(f"  elementy rzedu 3 (kolor): warstwy {colors3} (nietrywialne => nosza kolor)")
print(f"  centralnosc elementow rzedu 3: {[all(qk(qmul(g,h))==qk(qmul(h,g)) for h in G) for g in ord3[:2]]} (False => kolor NIE centralny)")
print()
print("WNIOSEK (§5.4 A + hook §0):")
print(" M=S^3/2T = przestrzen reperu: pi_3=Z (spin B), pi_1=2T (kolor). Jedno M.")
print(" 2T = Q8|x|Z3: spin(-1) w bezbarwnym Q8, kolor=Z3 iloraz. Rozne skladniki,")
print(" powiazane poluproduktem => KANDYDAT na lock §0 (do rozstrzygniecia nastepna strzalka).")
print(" Zero czlonu Skyrme, zero Derricka, zero plaskiego tla => background-free.")
