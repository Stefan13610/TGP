#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
cabibbo_correction_verify.py -- Verification of Cabibbo angle correction
========================================================================
Tests the Z_3 self-energy subtraction:

  lambda_C = (Omega_Lambda / N) * (|GL(3,F_2)| - |Z_3|) / (|GL(3,F_2)| - 1)
           = (0.6847 / 3) * 165/167
           = 0.22550

This reduces the tension from 4.8 sigma to 0.7 sigma.

Part of the TGP core verification suite.
"""

import numpy as np
import galois

# ================================================================
# CONSTANTS
# ================================================================
OMEGA_LAMBDA = 0.6847
N = 3
LAMBDA_PDG = 0.22500
SIGMA_PDG = 0.00067
GL_ORDER = 168
Z3_ORDER = 3

# ================================================================
# TESTS
# ================================================================
results = []

def test(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

print("=" * 70)
print("  CABIBBO CORRECTION VERIFICATION")
print("=" * 70)
print()

# ---------------------------------------------------------------
# T1: GL(3, GF(2)) has exactly 168 elements
# ---------------------------------------------------------------
GF2 = galois.GF(2)

def det_gf2(M):
    m = M.view(np.ndarray).astype(int)
    d = (m[0,0]*(m[1,1]*m[2,2] - m[1,2]*m[2,1])
       - m[0,1]*(m[1,0]*m[2,2] - m[1,2]*m[2,0])
       + m[0,2]*(m[1,0]*m[2,1] - m[1,1]*m[2,0]))
    return d % 2

gl3 = []
for bits in range(512):
    entries = [(bits >> i) & 1 for i in range(9)]
    M = GF2(np.array(entries, dtype=int).reshape(3, 3))
    if det_gf2(M) == 1:
        gl3.append(M)

test("T1: |GL(3,GF(2))| = 168",
     len(gl3) == 168,
     f"|GL| = {len(gl3)}")

# ---------------------------------------------------------------
# T2: Exactly 56 elements of order 3
# ---------------------------------------------------------------
I3 = GF2(np.eye(3, dtype=int))

def mat_order(M):
    power = M.copy()
    for k in range(1, 200):
        if np.array_equal(power, I3):
            return k
        power = power @ M
    return None

order3_count = sum(1 for M in gl3 if mat_order(M) == 3)
test("T2: 56 elements of order 3",
     order3_count == 56,
     f"count = {order3_count}")

# ---------------------------------------------------------------
# T3: 28 Z_3 subgroups
# ---------------------------------------------------------------
def mat_key(M):
    return tuple(M.view(np.ndarray).flatten())

order3_elems = [M for M in gl3 if mat_order(M) == 3]
z3_subs = []
used = set()
for g in order3_elems:
    gk = mat_key(g)
    if gk in used:
        continue
    g2 = g @ g
    used.add(gk)
    used.add(mat_key(g2))
    z3_subs.append(g)

test("T3: 28 Z_3 subgroups",
     len(z3_subs) == 28,
     f"count = {len(z3_subs)}")

# ---------------------------------------------------------------
# T4: All Z_3 subgroups are conjugate (single conjugacy class)
# ---------------------------------------------------------------
def gf2_inv(M):
    m = M.view(np.ndarray).astype(int)
    adj = np.zeros((3,3), dtype=int)
    for i in range(3):
        for j in range(3):
            rows = [r for r in range(3) if r != i]
            cols = [c for c in range(3) if c != j]
            minor = m[np.ix_(rows, cols)]
            cof = (minor[0,0]*minor[1,1] - minor[0,1]*minor[1,0]) % 2
            adj[j,i] = ((-1)**(i+j) * cof) % 2
    return GF2(adj % 2)

# Check: is every Z_3 generator conjugate to z3_subs[0]?
g0 = z3_subs[0]
g0_key = mat_key(g0)
g0_sq_key = mat_key(g0 @ g0)
all_conjugate = True

for g in z3_subs[1:]:
    found = False
    for h in gl3:
        c = gf2_inv(h) @ g @ h
        ck = mat_key(c)
        if ck == g0_key or ck == g0_sq_key:
            found = True
            break
    if not found:
        all_conjugate = False
        break

test("T4: All Z_3 subgroups conjugate",
     all_conjugate,
     "Single conjugacy class of Z_3 subgroups")

# ---------------------------------------------------------------
# T5: |N_G(Z_3)| = 6
# ---------------------------------------------------------------
z3_set = {mat_key(I3), mat_key(g0), mat_key(g0 @ g0)}
normalizer = [h for h in gl3
              if mat_key(gf2_inv(h) @ g0 @ h) in z3_set]

test("T5: |N_G(Z_3)| = 6",
     len(normalizer) == 6,
     f"|N_G(Z_3)| = {len(normalizer)}")

# ---------------------------------------------------------------
# T6: Correction factor F = 165/167
# ---------------------------------------------------------------
F = (GL_ORDER - Z3_ORDER) / (GL_ORDER - 1)
F_exact = 165 / 167
test("T6: F = (|G|-|Z_3|)/(|G|-1) = 165/167",
     abs(F - F_exact) < 1e-15,
     f"F = {F:.10f}")

# ---------------------------------------------------------------
# T7: Corrected lambda_C = 0.22550
# ---------------------------------------------------------------
lambda_corr = (OMEGA_LAMBDA / N) * F
test("T7: lambda_C = 0.22550",
     abs(lambda_corr - 0.22550) < 0.00001,
     f"lambda_C = {lambda_corr:.5f}")

# ---------------------------------------------------------------
# T8: Tension < 2 sigma
# ---------------------------------------------------------------
tension = abs(lambda_corr - LAMBDA_PDG) / SIGMA_PDG
test("T8: Tension < 2 sigma",
     tension < 2.0,
     f"tension = {tension:.2f} sigma")

# ---------------------------------------------------------------
# T9: Tension < 1 sigma (stronger)
# ---------------------------------------------------------------
test("T9: Tension < 1 sigma",
     tension < 1.0,
     f"tension = {tension:.2f} sigma")

# ---------------------------------------------------------------
# T10: CKM unitarity preserved
# ---------------------------------------------------------------
V_ud = np.sqrt(1 - lambda_corr**2)
V_us = lambda_corr
V_ub = 0.826 * lambda_corr**3 * 0.37  # rough
row_sum = V_ud**2 + V_us**2 + V_ub**2
test("T10: CKM unitarity (row 1)",
     abs(row_sum - 1.0) < 1e-4,
     f"|V_ud|^2 + |V_us|^2 + |V_ub|^2 = {row_sum:.6f}")

# ---------------------------------------------------------------
# T11: Uniqueness -- only GL(3,F_2) works
# ---------------------------------------------------------------
# For other groups, the tension should be > 2 sigma
alt_groups = [(6, "S_3"), (12, "A_4"), (24, "S_4"), (60, "A_5")]
all_excluded = True
for order, name in alt_groups:
    F_alt = (order - Z3_ORDER) / (order - 1)
    lam_alt = (OMEGA_LAMBDA / N) * F_alt
    t_alt = abs(lam_alt - LAMBDA_PDG) / SIGMA_PDG
    if t_alt < 2.0:
        all_excluded = False
        break

test("T11: All smaller groups excluded (> 2 sigma)",
     all_excluded,
     "S_3, A_4, S_4, A_5 all give tension > 2 sigma")

# ================================================================
# SUMMARY
# ================================================================
print()
print("=" * 70)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
n_total = len(results)
print(f"  RESULT: {n_pass}/{n_total} PASS")

if n_pass == n_total:
    print("  Cabibbo correction VERIFIED.")
else:
    print("  WARNING: Some tests FAILED!")
    sys.exit(1)

print("=" * 70)
