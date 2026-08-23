#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TIER-B (honest, parameter-free): stable particle = PERIODIC braid (finite order
# in B3/center = PSL(2,Z)). NO fitting, NO reparametrization. Count them, see structure.
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

def matmul(M,N):
    a,b,c,d=M; e,f,g,h=N
    return (a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)
def inv(M):  # det=1
    a,b,c,d=M; return (d,-b,-c,a)
def neg(M):
    a,b,c,d=M; return (-a,-b,-c,-d)
def canon(M):            # PSL: identify M ~ -M
    return min(M, neg(M))
I=(1,0,0,1)
S1=(1,1,0,1); S1i=(1,-1,0,1); S2=(1,0,-1,1); S2i=(1,0,1,1)
GEN=[('s1',S1,(1,0)),('s1i',S1i,(1,0)),('s2',S2,(2,1)),('s2i',S2i,(2,1))]
def swap(p,ab):
    a,b=ab; p=list(p); p[a-1],p[b-1]=p[b-1],p[a-1]; return tuple(p)

# BFS enumerate distinct braids (matrix) with min complexity + a permutation
LMAX=10
best={I:(0,(1,2,3))}
frontier=[(I,(1,2,3))]
for L in range(1,LMAX+1):
    nxt=[]
    for (M,perm) in frontier:
        for (_,G,ab) in GEN:
            M2=matmul(M,G); p2=swap(perm,ab)
            if M2 not in best:
                best[M2]=(L,p2); nxt.append((M2,p2))
    frontier=nxt
print(f"distinct SL(2,Z) elements enumerated (|complexity|<= {LMAX}): {len(best)}")

# elliptic (finite order in PSL): trace in {-1,0,1}
def trace(M): return M[0]+M[3]
def psl_order(M):
    t=trace(M)
    if t==0: return 2      # SL order 4 -> PSL order 2
    if t in (1,-1): return 3  # SL order 6 or 3 -> PSL order 3
    return None
ell = {}   # canon-matrix -> (complexity, perm, order)
for M,(L,p) in best.items():
    if psl_order(M) is not None:
        cM=canon(M)
        if cM not in ell or L<ell[cM][0]:
            ell[cM]=(L,p,psl_order(M))
print(f"elliptic (periodic/stable) PSL elements found: {len(ell)}")

# conjugacy classes in PSL: g M g^-1 ~ +-M2, using enumerated matrices as conjugators
conjs=list(best.keys())
elems=list(ell.keys())
idx={m:i for i,m in enumerate(elems)}
parent=list(range(len(elems)))
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb
for g in conjs:
    gi=inv(g)
    for m in elems:
        conj=canon(matmul(matmul(g,m),gi))
        if conj in idx:
            union(idx[m], idx[conj])
# gather classes
from collections import defaultdict
classes=defaultdict(list)
for m in elems:
    classes[find(idx[m])].append(m)

print("\n"+"="*70)
print("RESULT: conjugacy classes of stable (periodic) braids in PSL(2,Z)")
print("="*70)
def perm_class(p):
    fx=sum(1 for i in range(3) if p[i]==i+1)
    return 'id' if fx==3 else ('transposition' if fx==1 else '3-cycle')
summary=[]
for k,ms in classes.items():
    # representative: smallest complexity
    rep=min(ms, key=lambda m: ell[m][0])
    L,p,order=ell[rep]
    summary.append((order, perm_class(p), L, len(ms)))
summary.sort()
print(f"{'PSL order':>10} {'permutation':>14} {'min complexity':>15} {'#elements':>10}")
for order,pc,L,n in summary:
    print(f"{order:>10} {pc:>14} {L:>15} {n:>10}")
print(f"\nTOTAL stable conjugacy classes (excl. identity): {len(summary)}")

print("\n"+"="*70)
print("HONEST CHECK: can these be '3 generations of ONE particle'?")
print("  (generations must share IDENTICAL external quantum numbers)")
print("="*70)
from collections import Counter
byqn=Counter((pc,) for order,pc,L,n in summary)
for qn,cnt in byqn.items():
    print(f"  permutation={qn[0]:>14}: {cnt} stable class(es)")
print("""
Reading:
  If '3 generations' were real here, we'd need 3 classes with the SAME
  permutation (same external QN) differing only by internal complexity (mass).
""")
