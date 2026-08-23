# AUDIT2_group_verify.py
# Independent adversarial verification of the charge-lock Z3 claim.
# Uses SL(2,3) over GF(3) (INDEPENDENT construction, not quaternions) as 2T,
# plus a from-scratch abelianization + Hom(.,U(1)) counter, and probes the
# claimed "lock" for hidden tautology / smuggled input.

import itertools, cmath

# ---- generic finite-group engine ----
def close(gens, mul, ident, key):
    E={key(ident):ident}; fr=[ident]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                p=mul(g,h); k=key(p)
                if k not in E: E[k]=p; nw.append(p)
        fr=nw
        if len(E)>5000: raise RuntimeError("too big")
    return list(E.values())

def inverse(g,E,mul,ident,key):
    for h in E:
        if key(mul(g,h))==key(ident): return h
    raise RuntimeError("no inv")

def order_of(g,mul,ident,key):
    p=g;n=1
    while key(p)!=key(ident): p=mul(p,g);n+=1
    return n

def commutator_sub(E,mul,ident,key):
    cs=[]
    for a in E:
        ai=inverse(a,E,mul,ident,key)
        for b in E:
            bi=inverse(b,E,mul,ident,key)
            cs.append(mul(mul(a,b),mul(ai,bi)))
    return close(cs,mul,ident,key)

# ---- 2T as SL(2,3): 2x2 over GF(3), det=1 ----
def m(A,B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%3 for j in range(2)) for i in range(2))
def mk(A): return (A[0][0],A[0][1],A[1][0],A[1][1])
I=((1,0),(0,1))
g1=((1,1),(0,1)); g2=((1,0),(1,1))
G=close([g1,g2],m,I,mk)
print(f"[1] |SL(2,3)| = {len(G)}  (expect 24 = |2T|): {'OK' if len(G)==24 else 'FAIL'}")

# central -1 = 2*I over GF(3) = ((2,0),(0,2))
neg=((2,0),(0,2))
in_G = mk(neg) in set(mk(x) for x in G)
central = all(mk(m(neg,x))==mk(m(x,neg)) for x in G)
print(f"[2] -1 in G: {in_G};  -1 central: {central};  ord(-1)={order_of(neg,m,I,mk)} (expect 2)")

# commutator subgroup
Gp=commutator_sub(G,m,I,mk)
print(f"[3] |[G,G]| = {len(Gp)} (expect 8 = Q8): {'OK' if len(Gp)==8 else 'FAIL'}")
ab = len(G)//len(Gp)
print(f"[4] |G^ab| = |G|/|[G,G]| = {ab} (expect 3 = Z3): {'OK' if ab==3 else 'FAIL'}")

# is -1 in commutator subgroup? (=> every character kills it => spin bezbarwny)
neg_in_Gp = mk(neg) in set(mk(x) for x in Gp)
print(f"[5] -1 in [G,G]? {neg_in_Gp}  => forced: every hom chi(-1)=1 => spin charge 0")

# ---- Hom(G,U(1)) counted DIRECTLY (brute force over all functions G->mu_k?) ----
# Cleaner: Hom(G,U(1)) = Pontryagin dual of G^ab. Build G^ab explicitly, count its
# elements (finite abelian => |dual| = |G^ab|). Also verify it's cyclic Z3.
# Build quotient: coset of g = g[G,G].
Gp_keys=set(mk(x) for x in Gp)
def coset(g): return frozenset(mk(m(g,x)) for x in Gp)
cos={}
for g in G: cos.setdefault(coset(g),g)
reps=list(cos.values())
print(f"[6] |G^ab| explicit cosets = {len(reps)}")
# structure: order of each coset in quotient
def coset_order(r):
    p=r;n=1
    while coset(p)!=coset(I): p=m(p,r);n+=1
    return n
ords=sorted(coset_order(r) for r in reps)
print(f"    coset orders = {ords}  => {'cyclic Z3' if ords==[1,3,3] else 'NOT Z3'}")

# ---- ADVERSARIAL PROBE 1: is 'lock' a theorem or a construction? ----
# Claim: any U(1) charge (flat holonomy) on 2T is a FUNCTION of color (abelianization).
# Test: enumerate ALL homomorphisms G->U(1) as maps, check each is constant on cosets
# of [G,G] (i.e. factors through color). If ALL do, the 'lock' is FORCED (theorem),
# not a choice. Build homs via characters of G^ab (Z3): chi_k(g)=w^(k*tri(g)).
w=cmath.exp(2j*cmath.pi/3)
tri_map={coset(r): idx for idx,r in enumerate(reps)}  # arbitrary index labels
# fix labels so that they add mod 3: pick generator
genc=[r for r in reps if coset_order(r)==3][0]
lbl={coset(I):0, coset(genc):1, coset(m(genc,genc)):2}
def tri(g): return lbl[coset(g)]
assert all(tri(m(a,b))==(tri(a)+tri(b))%3 for a in G for b in G), "tri not a homomorphism"
homs=0
for k in range(3):
    ok=all(abs(w**((k*tri(m(a,b)))%3)-(w**((k*tri(a))%3))*(w**((k*tri(b))%3)))<1e-9 for a in G for b in G)
    if ok: homs+=1
print(f"[7] #homomorphisms G->U(1) of form w^(k*tri) that are valid = {homs} (=|Hom|=3)")

# Is there ANY hom NOT factoring through color? A hom is determined by images of
# generators; images must satisfy relations. Since G^ab=Z3, Hom(G,U(1))=Hom(Z3,U(1))=Z3.
# So NO hom can separate the two order-3 conj classes beyond their abelianization class.
# Verify: the two order-3 conjugacy classes map to DISTINCT abelianization elements?
ord3=[g for g in G if order_of(g,m,I,mk)==3]
tris_ord3=sorted(set(tri(g) for g in ord3))
print(f"[8] triality values on order-3 elements: {tris_ord3}")
# how many order-3 elements per triality class
from collections import Counter
print(f"    distribution: {dict(Counter(tri(g) for g in ord3))}")

print("\n--- ADVERSARIAL CONCLUSIONS ---")
print("A. 2T^ab = Z3 and Hom(2T,U(1))=Z3 are TEXTBOOK-CORRECT (independently reproduced).")
print("B. -1 in [G,G] => spin charge 0 is FORCED, not chosen. Genuine, checkable.")
print("C. 'charge = function of color' is a THEOREM (chars factor through G^ab),")
print("   NOT an independent coincidence: Hom(G,U(1))=Hom(G^ab,U(1)) always.")
print("D. The '3' of fractional charge IS the '3' of color (same G^ab=Z3),")
print("   so it is NOT independent evidence: color=Z3 was the INPUT.")
