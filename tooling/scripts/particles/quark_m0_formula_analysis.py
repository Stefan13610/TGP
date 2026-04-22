"""
Deep analysis of the best m_0 formula: m_0 = A * m_3/m_1

This is circular since m_3 is what we want to predict.
But let's check: what does it MEAN physically?

m_0 = A * r_31 means the Koide shift is proportional to the
third-to-first generation mass ratio.

Combined with shifted Koide K(m_i + m_0) = 2/3, this gives
a self-consistent equation for r_31.

Also explore: the NON-circular formulae.
"""
import numpy as np
from scipy.optimize import brentq

# PDG masses
m_d, m_s, m_b = 4.67, 93.4, 4180.0
m_u, m_c, m_t = 2.16, 1270.0, 172760.0

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
# 1. Self-consistent equation: m_0 = A * m_3/m_1
# ============================================================
print("=" * 65)
print("1. SELF-CONSISTENT: m_0 = A * r_31, K(m+m_0) = 2/3")
print("=" * 65)

# Given m_1, m_2, and the formula m_0 = A * m_3/m_1:
# K(m_1 + A*m_3/m_1, m_2 + A*m_3/m_1, m_3 + A*m_3/m_1) = 2/3
# This is one equation in one unknown (m_3), with A determined from data.

A_fit = 0.0246
print(f"  A = {A_fit}")

# For down sector: solve for m_3 = m_b
def self_consistent_down(m3):
    m0 = A_fit * m3 / m_d
    return koide(m_d + m0, m_s + m0, m3 + m0) - 2/3

m_b_pred = brentq(self_consistent_down, 100, 50000)
print(f"\n  Down sector (m_d={m_d}, m_s={m_s}):")
print(f"    Predicted m_b = {m_b_pred:.1f} MeV  (PDG: {m_b:.0f})")
print(f"    Accuracy: {abs(m_b_pred - m_b)/m_b*100:.2f}%")

# For up sector: solve for m_3 = m_t
def self_consistent_up(m3):
    m0 = A_fit * m3 / m_u
    return koide(m_u + m0, m_c + m0, m3 + m0) - 2/3

m_t_pred = brentq(self_consistent_up, 1000, 500000)
print(f"\n  Up sector (m_u={m_u}, m_c={m_c}):")
print(f"    Predicted m_t = {m_t_pred:.0f} MeV  (PDG: {m_t:.0f})")
print(f"    Accuracy: {abs(m_t_pred - m_t)/m_t*100:.2f}%")

# ============================================================
# 2. What value of A gives EXACT predictions?
# ============================================================
print("\n" + "=" * 65)
print("2. EXACT A FOR EACH SECTOR")
print("=" * 65)

A_d = m0_d * m_d / m_b
A_u = m0_u * m_u / m_t
print(f"  A(down) = m_0 * m_1 / m_3 = {A_d:.6f}")
print(f"  A(up)   = m_0 * m_1 / m_3 = {A_u:.6f}")
print(f"  Ratio A_u/A_d = {A_u/A_d:.4f}")
print(f"  Average: {(A_d+A_u)/2:.6f}")

# These are very close! A_d = 0.02452, A_u = 0.02481
# A ~ 0.0247

# What is 0.0247? Let's check known constants:
print(f"\n  A ~ {(A_d+A_u)/2:.4f}")
print(f"  1/(4*pi) = {1/(4*np.pi):.4f}")
print(f"  alpha_s/(2*pi) = {0.1174/(2*np.pi):.4f}")
print(f"  alpha_s^2 = {0.1174**2:.4f}")
print(f"  1/(8*pi^2) = {1/(8*np.pi**2):.4f}")
print(f"  alpha_EM = {1/137.036:.6f}")

# A ~ 0.0247 is close to:
# - 1/40 = 0.025
# - alpha_s^2 / (4*pi) = 0.1174^2 / (4*pi) = 0.001096 -- no
# - alpha_s / (4*pi) = 0.1174/(4*pi) = 0.00934 -- no

# ============================================================
# 3. NON-CIRCULAR APPROACH: m_0 from (m_1, m_2) only
# ============================================================
print("\n" + "=" * 65)
print("3. NON-CIRCULAR: m_0 from (m_1, m_2) ONLY")
print("=" * 65)
print("  Formula: m_0 = B * m_2^2 / m_1  (closest non-circular)")

B = m0_d * m_d / m_s**2
B2 = m0_u * m_u / m_c**2
print(f"  B(down) = m_0 * m_1 / m_2^2 = {B:.6f}")
print(f"  B(up)   = m_0 * m_1 / m_2^2 = {B2:.6f}")
print(f"  Ratio B_u/B_d = {B2/B:.4f}")

# B_down = 0.01177, B_up = 0.002655 — factor 4.4 different, not universal

# Try: m_0 = C * m_2 * sqrt(r_21)
# m_0/m_2 = C * sqrt(m_2/m_1)
C_d = m0_d / (m_s * np.sqrt(m_s/m_d))
C_u = m0_u / (m_c * np.sqrt(m_c/m_u))
print(f"\n  m_0 = C * m_2 * sqrt(r_21):")
print(f"    C(down) = {C_d:.6f}")
print(f"    C(up)   = {C_u:.6f}")
print(f"    Ratio = {C_u/C_d:.4f}")

# Try: m_0 = D * sqrt(m_1 * m_2)
D_d = m0_d / np.sqrt(m_d * m_s)
D_u = m0_u / np.sqrt(m_u * m_c)
print(f"\n  m_0 = D * sqrt(m_1 * m_2):")
print(f"    D(down) = {D_d:.6f}")
print(f"    D(up)   = {D_u:.6f}")
print(f"    Ratio = {D_u/D_d:.4f}")

# ============================================================
# 4. THE ANSWER: m_0 / m_1 ~ r_21^p
# ============================================================
print("\n" + "=" * 65)
print("4. SCALING: m_0 / m_1 = f(r_21)")
print("=" * 65)

# m_0/m_1: down = 4.698, up = 917.4
# r_21: down = 20.0, up = 588.0
# Is m_0/m_1 a power of r_21?

ratio_d = m0_d / m_d  # 4.698
ratio_u = m0_u / m_u  # 917.4
r21_d = m_s / m_d  # 20.0
r21_u = m_c / m_u  # 588.0

# m_0/m_1 = K * r_21^p
# log(ratio) = log(K) + p * log(r_21)
# Two equations:
p = np.log(ratio_u / ratio_d) / np.log(r21_u / r21_d)
K = ratio_d / r21_d**p
print(f"  m_0/m_1 = {K:.6f} * r_21^{p:.4f}")
print(f"  Check: down = {K * r21_d**p:.2f} (actual: {ratio_d:.2f})")
print(f"  Check: up   = {K * r21_u**p:.1f} (actual: {ratio_u:.1f})")

# p ~ 1.6! Very close to phi = 1.618!
print(f"\n  Exponent p = {p:.4f}")
print(f"  phi = {(1+np.sqrt(5))/2:.4f}")
print(f"  Difference: {abs(p - (1+np.sqrt(5))/2):.4f} ({abs(p - (1+np.sqrt(5))/2)/((1+np.sqrt(5))/2)*100:.1f}%)")

phi = (1 + np.sqrt(5)) / 2

# If p = phi:
K_phi = ratio_d / r21_d**phi
print(f"\n  With p = phi exactly:")
print(f"  m_0/m_1 = {K_phi:.6f} * r_21^phi")
print(f"  Predictions:")
pred_d = K_phi * r21_d**phi * m_d
pred_u = K_phi * r21_u**phi * m_u
print(f"    m_0(down) = {pred_d:.1f} MeV  (actual: {m0_d:.1f})")
print(f"    m_0(up)   = {pred_u:.0f} MeV  (actual: {m0_u:.0f})")
err_d = abs(pred_d - m0_d) / m0_d * 100
err_u = abs(pred_u - m0_u) / m0_u * 100
print(f"    Errors: {err_d:.1f}%, {err_u:.1f}%")

# ============================================================
# 5. COMBINE: m_0 = K * m_1 * r_21^phi, then predict m_3
# ============================================================
print("\n" + "=" * 65)
print("5. FULL PREDICTION: m_0 = K * m_1 * r_21^phi -> m_3")
print("=" * 65)

# Using m_0 = K_phi * m_1 * r_21^phi with K_phi from down sector:
# Then K(m + m_0) = 2/3 -> predict m_3

for name, m1, m2, m3_pdg in [("down", m_d, m_s, m_b), ("up", m_u, m_c, m_t)]:
    r21 = m2 / m1
    m0_pred = K_phi * m1 * r21**phi

    # Predict m3 from shifted Koide
    def obj(m3):
        return koide(m1 + m0_pred, m2 + m0_pred, m3 + m0_pred) - 2/3

    try:
        m3_pred = brentq(obj, m2, 1e7)
    except:
        m3_pred = None

    print(f"\n  {name} sector:")
    print(f"    r_21 = {r21:.2f}")
    print(f"    m_0(predicted) = {m0_pred:.1f} MeV  (exact: {find_m0(m1, m2, m3_pdg):.1f})")
    if m3_pred:
        print(f"    m_3(predicted) = {m3_pred:.1f} MeV  (PDG: {m3_pdg:.1f})")
        print(f"    Accuracy: {abs(m3_pred - m3_pdg)/m3_pdg*100:.2f}%")

# ============================================================
# 6. LEPTON CHECK: does this formula give m_0 ~ 0?
# ============================================================
print("\n" + "=" * 65)
print("6. LEPTON CHECK: m_0 formula applied to leptons")
print("=" * 65)

m_e = 0.511
m_mu = 105.658
r21_lep = m_mu / m_e
m0_lep_pred = K_phi * m_e * r21_lep**phi
print(f"  r_21(lepton) = {r21_lep:.2f}")
print(f"  m_0(lepton, predicted) = {m0_lep_pred:.2f} MeV")
print(f"  m_0(lepton, actual) = ~0 MeV")
print(f"  PROBLEM: formula gives m_0 = {m0_lep_pred:.0f} MeV for leptons!")
print(f"  The formula is NOT universal — it only works for quarks.")
print(f"  NEED: a mechanism that gives m_0 = 0 for colorless particles")

# ============================================================
# 7. COLOR FACTOR: m_0 = K * m_1 * r_21^phi * delta_color
# ============================================================
print("\n" + "=" * 65)
print("7. WITH COLOR FACTOR: m_0 = K' * C_F * m_1 * r_21^phi")
print("=" * 65)
print("  Where delta_color = 0 for leptons, C_F = 4/3 for quarks")
# This doesn't help — it just rescales K by C_F

# The REAL question: why does m_0 = 0 for leptons but not quarks?
# Answer: leptons DON'T have QCD confinement energy
# The formula should be: m_0 = alpha_s * [something] for quarks
#                         m_0 = 0 for leptons (alpha_EM contribution negligible)

# Revised approach:
# m_0 = alpha_s * m_1 * r_21^p / (4*pi)
print("\n  Revised: m_0 = (alpha_s/4pi) * m_1 * r_21^p")
K_rev = m0_d / (0.1174/(4*np.pi) * m_d * r21_d**phi)
print(f"  K_rev(down) = {K_rev:.4f}")
K_rev2 = m0_u / (0.1174/(4*np.pi) * m_u * r21_u**phi)
print(f"  K_rev(up) = {K_rev2:.4f}")
print(f"  Ratio = {K_rev2/K_rev:.4f}")
# Still sector-dependent — alpha_s alone doesn't make it universal

print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"""
RESULT: m_0/m_1 = K * r_21^phi  (phi = golden ratio = 1.618...)

  K = {K_phi:.6f}  (fitted from down sector)
  p = {p:.4f} ~ phi = 1.618  (fitted from both sectors)

  Predictions (from m_1, m_2 only):
    m_b = [solved from shifted Koide]
    m_t = [solved from shifted Koide]

  The golden ratio appears AGAIN in TGP mass physics:
    1. g0^mu = phi * g0^e  (phi-FP mechanism)
    2. m_0/m_1 ~ r_21^phi  (Koide shift scaling)

  OPEN ISSUE: formula gives m_0 >> 0 for leptons (should be ~0).
  Need: color-dependent prefactor that vanishes for leptons.

  STATUS: Hypothesis — promising scaling law, not yet derived
  from first principles. The golden ratio exponent is striking
  but may be numerological (only 2 data points for fit).
""")
