#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TEST: does the SM color<->charge lock (3Q + t = 0 mod 3) fall out of braids?
#   color/triality t = A3 content of permutation (id->0, (123)->1, (132)->2)
#   charge candidate A: intrinsic writhe  Q = writhe/3  (structural, NOT free)
#   charge candidate B: added twists Sum(t_i)/3 (independent -> trivially no lock)
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from collections import defaultdict, Counter

def matmul(M,N):
    a,b,c,d=M; e,f,g,h=N
    return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
I=(1,0,0,1)
S1=(1,1,0,1); S1i=(1,-1,0,1); S2=(1,0,-1,1); S2i=(1,0,1,1)
# (name, matrix, writhe-delta, transposition)
GEN=[('s1',S1,+1,(1,0)),('s1i',S1i,-1,(1,0)),('s2',S2,+1,(2,1)),('s2i',S2i,-1,(2,1))]
def swap(p,ab):
    a,b=ab; p=list(p); p[a-1],p[b-1]=p[b-1],p[a-1]; return tuple(p)

# BFS: distinct braids (key = matrix,writhe) with permutation
LMAX=9
start=(I,0,(1,2,3))
best={(I,0):(1,2,3)}
frontier=[start]
allbraids=[start]
for L in range(1,LMAX+1):
    nxt=[]
    for (M,w,perm) in frontier:
        for (_,G,dw,ab) in GEN:
            M2=matmul(M,G); w2=w+dw; p2=swap(perm,ab)
            key=(M2,w2)
            if key not in best:
                best[key]=p2; st=(M2,w2,p2); nxt.append(st); allbraids.append(st)
    frontier=nxt
print(f"distinct braids enumerated: {len(allbraids)}")

def triality(p):
    return {(1,2,3):0,(2,3,1):1,(3,1,2):2}.get(p, None)  # None = odd (transposition)

# --- Sanity: Z2 lock  sign(perm) == writhe mod 2 ? ---
z2=Counter()
for (M,w,p) in allbraids:
    even = triality(p) is not None    # even permutation
    z2[(even, w%2)] += 1
print("\n[sanity] Z2 lock: even-permutation vs writhe mod 2")
print("  even-perm & writhe even :", z2[(True,0)])
print("  even-perm & writhe odd  :", z2[(True,1)])
print("  odd-perm  & writhe even :", z2[(False,0)])
print("  odd-perm  & writhe odd  :", z2[(False,1)])
print("  -> if only diagonal nonzero, sign is LOCKED to writhe parity")

# --- Main test: triality t vs writhe mod 3 (structural charge) ---
print("\n" + "="*60)
print("MAIN TEST: triality t  vs  writhe mod 3   (even-perm braids)")
print("  SM lock would require: writhe = -t (mod 3), i.e. a DIAGONAL table")
print("="*60)
joint=Counter()
for (M,w,p) in allbraids:
    t=triality(p)
    if t is not None:
        joint[(t, w%3)] += 1
print(f"\n{'t\\w%3':>6} {'0':>8} {'1':>8} {'2':>8}")
for t in (0,1,2):
    row=[joint[(t,r)] for r in (0,1,2)]
    print(f"{t:>6} {row[0]:>8} {row[1]:>8} {row[2]:>8}")

# assess: is it diagonal (locked) or spread (free)?
diag = sum(joint[(t,(-t)%3)] for t in (0,1,2))
total = sum(joint.values())
offdiag = total-diag
print(f"\n  cells consistent with lock (w = -t mod 3): {diag}")
print(f"  cells violating lock:                      {offdiag}")
print(f"  -> lock holds only if violating = 0")

# --- Candidate B note: free twists ---
print("\n" + "="*60)
print("Candidate B (approved mapping): charge = Sum(twist_i)/3, twists free")
print("="*60)
print("  twists are assigned INDEPENDENTLY of permutation ->")
print("  for ANY triality t, Sum(twist) can be any value mod 3.")
print("  => color and charge are independent DOF => NO lock (by construction).")
print("  (analytic, no computation needed)")
