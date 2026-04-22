"""
Attempt to derive m_0 from TGP first principles.

KEY INSIGHT from previous analysis:
  For up-type quarks: K=2/3 and r_21=phi give nearly the same m_0
  (1981 vs 2049, 3.4% difference)

This suggests a DEEP connection between Koide and golden ratio
in the up sector. Can we exploit this?

Physical picture: m_0 represents the CONFINEMENT ENERGY per quark.
In TGP, confinement arises from the soliton potential V(g).
The confining part is the cubic term: (beta/3) g^3.
For quarks (color charge), there's an additional confinement
from the gluon field, which in TGP maps to the substrate
N-body interaction.

Hypothesis: m_0 = (constituent mass) - (current mass) averaged
over the sector, i.e. the dynamical chiral symmetry breaking contribution.
"""
import numpy as np
from scipy.optimize import brentq

phi = (1 + np.sqrt(5)) / 2

# PDG masses
m_e, m_mu, m_tau = 0.511, 105.658, 1776.86
m_d, m_s, m_b = 4.67, 93.4, 4180.0
m_u, m_c, m_t = 2.16, 1270.0, 172760.0

# TGP parameters
g0e = 0.869
alpha_s = 0.1174  # TGP prediction
N_c = 3
N_f = 6

def koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2

def find_m0(m1, m2, m3):
    def obj(m0):
        return koide(m1+m0, m2+m0, m3+m0) - 2.0/3
    lo = -min(m1, m2, m3) + 1e-10
    return brentq(obj, lo, 1e8, xtol=1e-12)

m0_d = find_m0(m_d, m_s, m_b)
m0_u = find_m0(m_u, m_c, m_t)

# ============================================================
# 1. Constituent quark model comparison
# ============================================================
print("=" * 65)
print("1. CONSTITUENT QUARK MODEL COMPARISON")
print("=" * 65)

# Constituent masses (approximate, from hadron spectroscopy)
m_const_u = 336  # MeV
m_const_d = 340  # MeV
m_const_s = 486  # MeV
m_const_c = 1550  # MeV (varies by model)
m_const_b = 4730  # MeV (varies by model)
m_const_t = 174000  # MeV (close to pole mass)

# Dynamical mass = constituent - current
dm_u = m_const_u - m_u   # 334 MeV
dm_d = m_const_d - m_d   # 335 MeV
dm_s = m_const_s - m_s   # 393 MeV
dm_c = m_const_c - m_c   # 280 MeV
dm_b = m_const_b - m_b   # 550 MeV

print(f"  Dynamical mass (constituent - current):")
print(f"    u: {dm_u:.0f} MeV  d: {dm_d:.0f} MeV  s: {dm_s:.0f} MeV")
print(f"    c: {dm_c:.0f} MeV  b: {dm_b:.0f} MeV")
print(f"  Average (u,d): {(dm_u+dm_d)/2:.0f} MeV ~ Lambda_QCD")
print(f"  But m_0(down) = {m0_d:.0f} MeV << {dm_d:.0f} MeV")
print(f"  And m_0(up) = {m0_u:.0f} MeV >> {dm_u:.0f} MeV")
print(f"  -> Constituent mass model does NOT explain m_0")

# ============================================================
# 2. Geometric mean hypothesis
# ============================================================
print("\n" + "=" * 65)
print("2. GEOMETRIC MEAN HYPOTHESES")
print("=" * 65)

# Try: m_0 = sqrt(m_1 * m_2) or similar
gm_d = np.sqrt(m_d * m_s)
gm_u = np.sqrt(m_u * m_c)
print(f"  sqrt(m_1 * m_2):")
print(f"    down: {gm_d:.1f} MeV  (m_0 = {m0_d:.1f}, ratio = {m0_d/gm_d:.3f})")
print(f"    up:   {gm_u:.1f} MeV  (m_0 = {m0_u:.0f}, ratio = {m0_u/gm_u:.3f})")

gm_d3 = (m_d * m_s * m_b)**(1/3)
gm_u3 = (m_u * m_c * m_t)**(1/3)
print(f"\n  (m_1 * m_2 * m_3)^(1/3):")
print(f"    down: {gm_d3:.1f} MeV  (m_0 = {m0_d:.1f}, ratio = {m0_d/gm_d3:.3f})")
print(f"    up:   {gm_u3:.0f} MeV  (m_0 = {m0_u:.0f}, ratio = {m0_u/gm_u3:.3f})")

# ============================================================
# 3. TGP-specific: m_0 from alpha_s and sector parameters
# ============================================================
print("\n" + "=" * 65)
print("3. LOOKING FOR PATTERN: m_0 / geometric_mean")
print("=" * 65)

# m_0(down)/gm = 1.050, m_0(up)/gm = 37.86
# Very different -> geometric mean doesn't work simply

# What about: m_0 = f(r_21) * m_1?
# m_0/m_1: down = 4.698, up = 917.4
# r_21: down = 20.0, up = 588.0
# m_0/m_1 vs r_21: not linear

# INSIGHT: try m_0 = m_2 * h(K_PDG)
# where K_PDG is the un-shifted Koide value
K_d = koide(m_d, m_s, m_b)
K_u = koide(m_u, m_c, m_t)

print(f"  K(down) = {K_d:.6f},  m_0/m_2 = {m0_d/m_s:.4f}")
print(f"  K(up) = {K_u:.6f},  m_0/m_2 = {m0_u/m_c:.4f}")
print(f"  K(lep) = {koide(m_e, m_mu, m_tau):.6f},  m_0/m_2 = ~0")

# K_d = 0.731, m_0/m_2 = 0.235
# K_u = 0.849, m_0/m_2 = 1.56
# These don't have an obvious simple relation either

# ============================================================
# 4. THE SOLITON COUPLING APPROACH
# ============================================================
print("\n" + "=" * 65)
print("4. SOLITON COUPLING: N-BODY FORCE INTERPRETATION")
print("=" * 65)
print("""
In the N-body sector of TGP, quarks interact via:
  V_2 = C_eff * e^{-m_sp r} / r  (pairwise Yukawa)
  V_3 = irreducible 3-body from Phi^4

For a quark inside a hadron:
  - 3 quarks confined within R_had ~ 1/Lambda_QCD ~ 1 fm
  - Each quark is a soliton with A_tail coupling
  - The confinement energy ~ V_2(R_had) per pair
  - Total binding: ~ N_c*(N_c-1)/2 * V_2(R_had)

The ADDITIVE mass shift m_0 would then be:
  m_0 = [confinement energy per quark] = (N_c-1) * C_eff * e^{-m_sp R} / R

But this is the SAME for all quarks in a sector (same C_eff)...
unless C_eff differs between up and down types.
""")

# ============================================================
# 5. THE ISOSPIN-DEPENDENT SOLITON
# ============================================================
print("=" * 65)
print("5. ALTERNATIVE: m_0 = m_2 * (K_PDG - 2/3) / (2/3)")
print("=" * 65)

# This is a shot in the dark but let's check
# The "Koide excess" DK = K - 2/3 tells us how far from Koide
DK_d = K_d - 2/3
DK_u = K_u - 2/3

m0_test_d = m_s * DK_d / (2/3)
m0_test_u = m_c * DK_u / (2/3)
print(f"  m_2 * (K-2/3)/(2/3):")
print(f"    down: {m0_test_d:.1f} MeV  (actual m_0 = {m0_d:.1f})")
print(f"    up:   {m0_test_u:.0f} MeV  (actual m_0 = {m0_u:.0f})")

# Close-ish for down (9.1 vs 21.9) but not great
# Let me try other combinations

m0_test_d2 = m_s * DK_d
m0_test_u2 = m_c * DK_u
print(f"\n  m_2 * (K-2/3):")
print(f"    down: {m0_test_d2:.1f} MeV  (actual m_0 = {m0_d:.1f})")
print(f"    up:   {m0_test_u2:.0f} MeV  (actual m_0 = {m0_u:.0f})")

# ============================================================
# 6. SYSTEMATIC FORMULA SEARCH
# ============================================================
print("\n" + "=" * 65)
print("6. SYSTEMATIC SEARCH: m_0 = A * m_1^a * m_2^b * m_3^c")
print("=" * 65)

# m_0(down) = A * 4.67^a * 93.4^b * 4180^c = 21.94
# m_0(up) = A * 2.16^a * 1270^b * 172760^c = 1981.5
# Two equations, four unknowns (A, a, b, c) — underdetermined

# Fix A=1: 4.67^a * 93.4^b * 4180^c = 21.94
#           2.16^a * 1270^b * 172760^c = 1981.5
# Take log:
#   a*log(4.67) + b*log(93.4) + c*log(4180) = log(21.94)
#   a*log(2.16) + b*log(1270) + c*log(172760) = log(1981.5)

import itertools

# Try simple integer/half-integer exponents
best = None
best_err = 1e10

for a_num in range(-4, 5):
    for b_num in range(-4, 5):
        for c_num in range(-4, 5):
            a, b, c = a_num/2, b_num/2, c_num/2
            if a == 0 and b == 0 and c == 0:
                continue
            try:
                pred_d = m_d**a * m_s**b * m_b**c
                pred_u = m_u**a * m_c**b * m_t**c
                if pred_d <= 0 or pred_u <= 0:
                    continue
                err_d = abs(np.log(pred_d/m0_d))
                err_u = abs(np.log(pred_u/m0_u))
                total_err = err_d + err_u
                if total_err < best_err:
                    best_err = total_err
                    best = (a, b, c, pred_d, pred_u, err_d, err_u)
            except:
                continue

a, b, c, pd, pu, ed, eu = best
print(f"  Best fit: m_0 = m_1^{a} * m_2^{b} * m_3^{c}")
print(f"    down: {pd:.1f} MeV  (actual: {m0_d:.1f}, err: {ed*100:.1f}%)")
print(f"    up:   {pu:.0f} MeV  (actual: {m0_u:.0f}, err: {eu*100:.1f}%)")

# Also search with coefficients involving alpha_s, N_c, phi
print("\n  Top 5 best fits (A=1):")
results = []
for a_num in range(-4, 5):
    for b_num in range(-4, 5):
        for c_num in range(-4, 5):
            a, b, c = a_num/2, b_num/2, c_num/2
            if a == 0 and b == 0 and c == 0:
                continue
            try:
                pred_d = m_d**a * m_s**b * m_b**c
                pred_u = m_u**a * m_c**b * m_t**c
                if pred_d <= 0 or pred_u <= 0 or np.isnan(pred_d) or np.isnan(pred_u):
                    continue
                err = abs(np.log(pred_d/m0_d)) + abs(np.log(pred_u/m0_u))
                results.append((err, a, b, c, pred_d, pred_u))
            except:
                continue

results.sort()
for err, a, b, c, pd, pu in results[:10]:
    print(f"    m1^{a:+.1f} m2^{b:+.1f} m3^{c:+.1f}: "
          f"d={pd:.1f} ({pd/m0_d:.3f}x)  u={pu:.0f} ({pu/m0_u:.3f}x)  "
          f"err={err:.3f}")

# ============================================================
# 7. WITH PREFACTOR SEARCH
# ============================================================
print("\n" + "=" * 65)
print("7. WITH PREFACTOR: m_0 = A * m_1^a * m_2^b * m_3^c")
print("=" * 65)

# For each (a,b,c), find optimal A
results2 = []
for a_num in range(-4, 5):
    for b_num in range(-4, 5):
        for c_num in range(-4, 5):
            a, b, c = a_num/2, b_num/2, c_num/2
            if a == 0 and b == 0 and c == 0:
                continue
            try:
                xd = m_d**a * m_s**b * m_b**c
                xu = m_u**a * m_c**b * m_t**c
                if xd <= 0 or xu <= 0 or np.isnan(xd) or np.isnan(xu):
                    continue
                # Optimal A: minimize (log(A*xd/m0_d))^2 + (log(A*xu/m0_u))^2
                # => log A = (log(m0_d/xd) + log(m0_u/xu)) / 2
                logA = (np.log(m0_d/xd) + np.log(m0_u/xu)) / 2
                A = np.exp(logA)
                err = abs(np.log(A*xd/m0_d)) + abs(np.log(A*xu/m0_u))
                results2.append((err, A, a, b, c, A*xd, A*xu))
            except:
                continue

results2.sort()
print("  Top 10 fits:")
for err, A, a, b, c, pd, pu in results2[:10]:
    print(f"    {A:.4f} * m1^{a:+.1f} m2^{b:+.1f} m3^{c:+.1f}: "
          f"d={pd:.1f} ({pd/m0_d:.3f}x)  u={pu:.0f} ({pu/m0_u:.3f}x)  "
          f"err={err:.4f}")
