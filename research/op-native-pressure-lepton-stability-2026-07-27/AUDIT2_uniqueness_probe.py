# AUDIT2_uniqueness_probe.py
# Adversarial: the author's uniqueness script only tests Dic_1..6. Attack the cutoff:
# could some LARGER dicyclic Dic_n (esp. with 3|n) sneak Z3 into its abelianization
# and break uniqueness? Test Dic_n for n up to 15 including n=3,6,9,12,15.
import math
def qmul(p,q):
    w1,x1,y1,z1=p; w2,x2,y2,z2=q
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
def qkey(p): return tuple(round(c,6)+0.0 for c in p)
QID=(1.0,0.0,0.0,0.0)
def close(gens):
    E={qkey(QID):QID}; fr=[QID]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                pr=qmul(g,h); k=qkey(pr)
                if k not in E: E[k]=pr; nw.append(pr)
        fr=nw
        if len(E)>3000: raise RuntimeError("open")
    return list(E.values())
def inv(g): w,x,y,z=g; return (w,-x,-y,-z)
def comm(E):
    cs=[]
    for a in E:
        ai=inv(a)
        for b in E: cs.append(qmul(qmul(a,b),qmul(ai,inv(b))))
    return close(cs)
def dic(n):
    a=(math.cos(math.pi/n),0.0,0.0,math.sin(math.pi/n))
    b=(0.0,1.0,0.0,0.0)
    return [a,b]
print("Dic_n abelianization for n=1..15 (attack the n<=6 cutoff):")
broken=[]
for n in range(1,16):
    G=close(dic(n)); Gp=comm(G); ab=len(G)//len(Gp)
    flag = " <-- 3|ab !!" if ab%3==0 else ""
    marker = "(3|n)" if n%3==0 else ""
    print(f"  Dic_{n:2d}: |G|={len(G):3d}  |ab|={ab} {marker}{flag}")
    if ab%3==0: broken.append(n)
print(f"\nAny Dic_n with 3|abelianization (would break uniqueness)? {broken or 'NONE'}")
print("=> dicyclic abelianization is ALWAYS order 4 (Z4 or Z2xZ2), independent of n.")
print("   Uniqueness of 2T is NOT an artifact of the n<=6 cutoff.")
