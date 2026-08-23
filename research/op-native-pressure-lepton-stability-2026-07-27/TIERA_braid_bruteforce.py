#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TIER-A brute force: relational combinatorics of 3-strand braids (background-free).
# B3 = <s1,s2 | s1 s2 s1 = s2 s1 s2>. Faithful-mod-center handle: B3 -> SL(2,Z).
#   s1 -> [[1,1],[0,1]],  s2 -> [[1,0],[-1,1]]   (braid relation verified: ABA=BAB).
# Metastability proxy = geodesic complexity (min crossings to build the braid),
#   i.e. irreducibility under the braid moves -- NO metric, NO energy, NO 'attraction'.
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from itertools import product
from fractions import Fraction

# generators as 2x2 int matrices (row-major tuples)
def matmul(M,N):
    a,b,c,d = M; e,f,g,h = N
    return (a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)
I = (1,0,0,1)
S1  = (1,1,0,1);  S1i = (1,-1,0,1)
S2  = (1,0,-1,1); S2i = (1,0,1,1)
GEN = [('s1',S1,+1,(1,0)), ('s1i',S1i,-1,(1,0)),
       ('s2',S2,+1,(2,1)), ('s2i',S2i,-1,(2,1))]  # (name,matrix,exp,transposition idx a<->b in 1..3)

def apply_perm_transposition(p, ab):
    a,b = ab; p=list(p)
    p[a-1],p[b-1]=p[b-1],p[a-1]; return tuple(p)

# sanity: braid relation
assert matmul(matmul(S1,S2),S1)==matmul(matmul(S2,S1),S2), "braid relation FAILS"
print("braid relation s1s2s1=s2s1s2 in SL(2,Z): OK")

# BFS over Cayley graph -> geodesic length (complexity) of each distinct braid.
# state key = (matrix, expsum) uniquely fixes braid (matrix mod center + winding).
LMAX = 8
start = (I, 0, (1,2,3))       # (matrix, expsum, permutation)
best = {(I,0):0}
frontier = [start]
by_len = {0:[start]}
for L in range(1, LMAX+1):
    nxt=[]; seen_this=[]
    for (M,e,perm) in frontier:
        for (nm,G,de,ab) in GEN:
            M2 = matmul(M,G); e2 = e+de
            perm2 = apply_perm_transposition(perm, ab)
            key=(M2,e2)
            if key not in best:
                best[key]=L; st=(M2,e2,perm2); nxt.append(st); seen_this.append(st)
    by_len[L]=seen_this; frontier=nxt

print("\n" + "="*72)
print("STEP 1: braid spectrum by complexity (geodesic # of crossings)")
print("  distinct braids that FIRST appear at each complexity level")
print("="*72)
# permutation classes in S3
def perm_class(p):
    p=tuple(p)
    if p==(1,2,3): return 'id'
    ncyc = p.count  # not needed
    # count fixed points
    fx = sum(1 for i in range(3) if p[i]==i+1)
    if fx==1: return 'transposition'
    if fx==0: return '3-cycle'
    return 'id'
print(f"{'complexity':>10} {'#braids':>8} {'id':>5} {'transp':>7} {'3cyc':>5}")
for L in range(0,LMAX+1):
    lst=by_len[L]
    cid=sum(1 for s in lst if perm_class(s[2])=='id')
    ctr=sum(1 for s in lst if perm_class(s[2])=='transposition')
    c3 =sum(1 for s in lst if perm_class(s[2])=='3-cycle')
    print(f"{L:>10} {len(lst):>8} {cid:>5} {ctr:>7} {c3:>5}")

print("\n" + "="*72)
print("STEP 2: the 2-crossing braids (Bilson-Thompson regime)")
print("="*72)
two = by_len[2]
print(f"# distinct 2-crossing braids: {len(two)}")
for (M,e,perm) in sorted(two, key=lambda s:(s[1],s[0])):
    print(f"   matrix={M} expsum={e:+d} perm={perm} class={perm_class(perm)}")

print("\n" + "="*72)
print("STEP 3: charge layer -- twists t_i in {-1,0,+1} per strand, Q = (sum t)/3")
print("  Does charge-in-thirds emerge? Do multiplicities match ONE SM generation?")
print("="*72)
from collections import Counter
mult = Counter()
for t in product((-1,0,1), repeat=3):
    mult[sum(t)] += 1
print(f"{'sum t':>6} {'Q':>7} {'# twist-vectors':>16}")
for s in range(-3,4):
    print(f"{s:>6} {str(Fraction(s,3)):>7} {mult[s]:>16}")
print("\n  SM one generation (charges, with color): "
      "nu:0(x1) e:-1(x1) u:+2/3(x3) d:-1/3(x3)  [+ antiparticles]")
print("  twist-only multiplicities: Q=+-1:1, +-2/3:3, +-1/3:6, 0:7")
print("  -> 'thirds' appear for free; '3' appears at +-2/3 (echo of 3 colors);")
print("     but +-1/3 gives 6 not 3 -> twist ALONE undercounts; needs braid x twist.")

print("\n" + "="*72)
print("STEP 4: is there a natural 'stable + ladder' and a natural 3?")
print("="*72)
print(f"  complexity-0: {len(by_len[0])} braid (identity) = unique simplest object")
print(f"  complexity-1: {len(by_len[1])} braids")
print(f"  complexity-2: {len(by_len[2])} braids  <- BT particle regime")
print(f"  growth ratio 2->3->4: {len(by_len[3])/max(1,len(by_len[2])):.2f}, "
      f"{len(by_len[4])/max(1,len(by_len[3])):.2f}")
print("="*72)
