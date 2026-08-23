# RIBBONS_closure_solitons.py
# Strzalka kwarki -> solitony: czy DOMYKANIE wpada automatycznie?
# BEZ ladunku (decyzja: ladunek jest relacyjny, warstwa wyzej).
#
# Pytanie topologiczne (klasyczne dla defektow nieabelowych):
#   kolekcja wstazek (dysklinacji) o klasach C1..Cn moze sie DOMKNAC/anihilowac
#   <=> istnieja reprezentanci gi in Ci takie, ze g1*g2*...*gn = e.
#   (dla defektow nieabelowych klasa jest okreslona do sprzezenia -> pytamy o ISTNIENIE)
# Jesli produkt = -1 (a nie e), konfiguracja domyka sie do obiektu niosacego
# CENTRALNY element -1 = spin/fermion (FR (-1)^B) - to NIE jest prozniowe domkniecie.
#
# Liczymy:
#  [1] n=1: czy pojedyncza wstazka moze sie domknac? (test konfinementu topologicznego)
#  [2] n=2: pary (kandydat: mezon = kolor + antykolor)
#  [3] n=3: trojki (kandydat: barion = 3 kolory)
#  [4] co zostaje jako RESZTA (e vs -1) dla trojek trialnosci (1,1,1)

import itertools
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
        for b in E: cs.append(qmul(qmul(a,b),qmul(ai,inv(b))))
    return close(cs)

G=close([(0.,1.,0.,0.),(0.5,0.5,0.5,0.5)]); assert len(G)==24
Q8=comm_sub(G); Q8k=set(qk(x) for x in Q8)
def coset(g): return frozenset(qk(qmul(g,x)) for x in Q8)
C0=coset(QID); gen=next(g for g in G if coset(g)!=C0)
lab={C0:0, coset(gen):1, coset(qmul(gen,gen)):2}
def tri(g): return lab[coset(g)]
assert all(tri(qmul(a,b))==(tri(a)+tri(b))%3 for a in G for b in G)

# klasy sprzezonosci
def conj_classes(E):
    seen=set(); out=[]
    for e in E:
        if qk(e) in seen: continue
        cl={qk(qmul(qmul(g,e),inv(g))) for g in E}
        out.append(sorted(cl)); seen|=cl
    return out
CL=conj_classes(G)
byk={qk(g):g for g in G}
def cls_info(cl):
    rep=byk[cl[0]]
    return dict(size=len(cl), ordr=order(rep), tri=tri(rep))
info=[cls_info(c) for c in CL]
print("KLASY WSTAZEK w 2T (7 typow):")
for i,(c,d) in enumerate(zip(CL,info)):
    print(f"  C{i}: rozmiar={d['size']:2d}  rzad={d['ordr']}  trialnosc(kolor)={d['tri']}")

def can_reach(classes, target):
    """czy istnieja reprezentanci g_i in C_i z iloczynem = target"""
    reach={qk(QID)}
    for c in classes:
        nxt=set()
        for r in reach:
            for gk in c:
                nxt.add(qk(qmul(byk[r], byk[gk])))
        reach=nxt
    return qk(target) in reach

print("\n[1] n=1: czy POJEDYNCZA wstazka moze sie domknac (anihilowac do prozni)?")
for i,(c,d) in enumerate(zip(CL,info)):
    ok=can_reach([c], QID)
    print(f"  C{i} (rzad {d['ordr']}, kolor {d['tri']}): domyka sie? {ok}")
print("  => tylko klasa trywialna. KAZDA nietrywialna wstazka jest topologicznie ZWIAZANA.")

print("\n[2] n=2 pary (kandydat: mezon = kolor+antykolor). Domkniecie do e:")
found2=[]
for i,j in itertools.combinations_with_replacement(range(len(CL)),2):
    if can_reach([CL[i],CL[j]], QID):
        t=(info[i]['tri']+info[j]['tri'])%3
        found2.append((i,j,t))
for i,j,t in found2:
    print(f"  C{i}+C{j}: domyka sie do e   (suma trialnosci = {t})")
print(f"  => par domykajacych: {len(found2)}; WSZYSTKIE maja sume trialnosci 0?",
      all(t==0 for _,_,t in found2))

print("\n[3] n=3 trojki o trialnosciach (1,1,1) - kandydat BARION:")
tri1=[i for i,d in enumerate(info) if d['tri']==1]
print(f"  klasy trialnosci 1: {['C%d'%i for i in tri1]}")
for combo in itertools.combinations_with_replacement(tri1,3):
    to_e   = can_reach([CL[k] for k in combo], QID)
    to_neg = can_reach([CL[k] for k in combo], NEG)
    names="+".join("C%d"%k for k in combo)
    print(f"  {names}: ->e? {to_e}   ->(-1)? {to_neg}")

print("\n[4] Czy trialnosc-0 WYSTARCZA do domkniecia? (test koniecznosci vs dostatecznosci)")
# wszystkie trojki o sumie trialnosci 0
tot=0; closable=0
for combo in itertools.combinations_with_replacement(range(len(CL)),3):
    t=sum(info[k]['tri'] for k in combo)%3
    if t!=0: continue
    tot+=1
    if can_reach([CL[k] for k in combo], QID): closable+=1
print(f"  trojki o sumie trialnosci 0: {tot};  z nich domykalnych do e: {closable}")
print(f"  => trialnosc-0 {'WYSTARCZA' if closable==tot else 'jest KONIECZNA ale NIE WYSTARCZA'}")

print("\n[5] Co zostaje po zneutralizowaniu koloru? (reszta w Q8 = rdzen bezbarwny)")
print(f"  suma trialnosci 0  <=>  iloczyn calkowity lezy w Q8 (|Q8|={len(Q8)})")
print(f"  Q8 zawiera -1 (spin/fermion): {qk(NEG) in Q8k}")
print("  => domkniecie koloru NIE usuwa sektora Q8: pozostaje etykieta spinowa (e vs -1).")
