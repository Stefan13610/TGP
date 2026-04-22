"""
Deep analysis of substrate mass ratios after shifted Koide correction.

KEY OBSERVATION: r_21(up, substrate) = 1.64 ~ phi = 1.618!

This suggests that the SUBSTRATE masses (before QCD correction)
may have a universal golden-ratio structure, with:
  m_substrate(2) / m_substrate(1) ~ phi for ALL sectors!

For leptons: phi-FP gives g0^mu = phi * g0^e, and A_tail^4 ratio = 206.8
For quarks: the MEASURED r_21 is NOT phi, but the SUBSTRATE r_21 might be.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

phi = (1 + np.sqrt(5)) / 2  # golden ratio

# PDG masses
masses = {
    'lepton': (0.511, 105.658, 1776.86),
    'down': (4.67, 93.4, 4180.0),
    'up': (2.16, 1270.0, 172760.0),
}

def koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2

def find_m0(m1, m2, m3):
    def obj(m0):
        return koide(m1+m0, m2+m0, m3+m0) - 2.0/3
    lo = -min(m1, m2, m3) + 1e-10
    return brentq(obj, lo, 1e8, xtol=1e-12)

# ============================================================
# 1. Substrate mass ratios
# ============================================================
print("=" * 65)
print("1. SUBSTRATE MASS RATIOS (after shifted Koide m_0 removal)")
print("=" * 65)

for name, (m1, m2, m3) in masses.items():
    m0 = find_m0(m1, m2, m3)
    ms1, ms2, ms3 = m1 + m0, m2 + m0, m3 + m0
    r21 = ms2 / ms1
    r31 = ms3 / ms1
    print(f"\n  {name}:")
    print(f"    m_0 = {m0:.2f} MeV")
    print(f"    substrate masses: ({ms1:.2f}, {ms2:.2f}, {ms3:.1f})")
    print(f"    r_21(substrate) = {r21:.4f}  (phi = {phi:.4f}, ratio/phi = {r21/phi:.4f})")
    print(f"    r_31(substrate) = {r31:.2f}")
    print(f"    K(substrate) = {koide(ms1, ms2, ms3):.8f}")

# ============================================================
# 2. REVERSE APPROACH: What m_0 gives r_21(substrate) = phi?
# ============================================================
print("\n" + "=" * 65)
print("2. REVERSE: What m_0 gives r_21(substrate) = phi EXACTLY?")
print("=" * 65)

for name, (m1, m2, m3) in masses.items():
    # r_21(sub) = (m2 + m0) / (m1 + m0) = phi
    # => m2 + m0 = phi * (m1 + m0)
    # => m0 = (m2 - phi * m1) / (phi - 1)
    m0_phi = (m2 - phi * m1) / (phi - 1)
    ms1 = m1 + m0_phi
    ms2 = m2 + m0_phi
    ms3 = m3 + m0_phi
    K_val = koide(ms1, ms2, ms3) if ms1 > 0 and ms2 > 0 and ms3 > 0 else float('nan')

    # Also: what m_3 would Koide predict?
    if ms1 > 0 and ms2 > 0:
        s = np.sqrt(ms1) + np.sqrt(ms2)
        S = ms1 + ms2
        disc = 16*s**2 - 4*(3*S - 2*s**2)
        if disc > 0:
            x = (4*s + np.sqrt(disc)) / 2
            m3_koide = x**2
            m3_pred = m3_koide - m0_phi
        else:
            m3_pred = None
    else:
        m3_pred = None

    print(f"\n  {name}:")
    print(f"    m_0(phi) = {m0_phi:.2f} MeV")
    if ms1 > 0:
        print(f"    r_21(substrate) = {ms2/ms1:.6f} = phi (exact)")
        print(f"    K(substrate) = {K_val:.6f}")
        if m3_pred:
            print(f"    Koide-predicted m_3 = {m3_pred:.1f} MeV  (PDG: {m3:.1f})")
            print(f"    Accuracy: {abs(m3_pred - m3)/m3*100:.2f}%")
    else:
        print(f"    substrate m_1 = {ms1:.2f} (negative!)")

# ============================================================
# 3. TWO-PARAMETER FIT: m_0 from both K=2/3 AND r21=phi
# ============================================================
print("\n" + "=" * 65)
print("3. OVER-DETERMINED SYSTEM: K=2/3 AND r_21=phi simultaneously?")
print("=" * 65)

for name, (m1, m2, m3) in masses.items():
    m0_koide = find_m0(m1, m2, m3)
    m0_phi = (m2 - phi * m1) / (phi - 1)

    print(f"\n  {name}:")
    print(f"    m_0(Koide) = {m0_koide:.2f}")
    print(f"    m_0(phi)   = {m0_phi:.2f}")
    if abs(m0_koide) > 0.01:
        print(f"    ratio = {m0_phi/m0_koide:.4f}")

    # Can we satisfy BOTH with a SINGLE m_0?
    # m_0(Koide) != m_0(phi) in general
    # But if close, it's interesting

    # The up-type case: m0_koide=1981.5, m0_phi=?
    # Let's check

# ============================================================
# 4. ALTERNATIVE: What if m_0 is the SAME for all quarks?
# ============================================================
print("\n" + "=" * 65)
print("4. UNIVERSAL m_0 FOR ALL QUARKS?")
print("=" * 65)

# Minimize sum of (K - 2/3)^2 over both sectors with shared m_0
def total_residual(m0):
    K_d = koide(m_d+m0, m_s+m0, m_b+m0)
    K_u = koide(m_u+m0, m_c+m0, m_t+m0)
    return (K_d - 2/3)**2 + (K_u - 2/3)**2

m_d, m_s, m_b = 4.67, 93.4, 4180.0
m_u, m_c, m_t = 2.16, 1270.0, 172760.0

result = minimize_scalar(total_residual, bounds=(0, 5000), method='bounded')
m0_shared = result.x
print(f"  Best shared m_0 = {m0_shared:.1f} MeV")
print(f"  K(down, shifted) = {koide(m_d+m0_shared, m_s+m0_shared, m_b+m0_shared):.6f}")
print(f"  K(up, shifted)   = {koide(m_u+m0_shared, m_c+m0_shared, m_t+m0_shared):.6f}")
print(f"  Total residual = {total_residual(m0_shared):.2e}")
print(f"  (Individual residuals need m_0(d)=21.9, m_0(u)=1981.5 -- too different for shared)")

# ============================================================
# 5. ELECTRIC CHARGE CONNECTION
# ============================================================
print("\n" + "=" * 65)
print("5. ELECTRIC CHARGE / WEAK ISOSPIN CONNECTION")
print("=" * 65)

m0_down = find_m0(*masses['down'])
m0_up = find_m0(*masses['up'])

Q_d = -1/3  # down-type charge
Q_u = 2/3   # up-type charge
T3_d = -1/2  # weak isospin
T3_u = 1/2

print(f"  m_0(down) = {m0_down:.1f},  Q = {Q_d:.3f},  T3 = {T3_d}")
print(f"  m_0(up)   = {m0_up:.0f},  Q = {Q_u:.3f},  T3 = {T3_u}")
print(f"  m_0(up)/m_0(down) = {m0_up/m0_down:.1f}")
print(f"  Q_u/Q_d = {Q_u/Q_d:.1f}")
print(f"  (Q_u/Q_d)^2 = {(Q_u/Q_d)**2:.1f}")
print(f"  m_0 ~ Q^2? down: {m0_down/Q_d**2:.0f}, up: {m0_up/Q_u**2:.0f}")

# ============================================================
# 6. DEEPER: m_0 as function of r_21
# ============================================================
print("\n" + "=" * 65)
print("6. m_0 AS FUNCTION OF r_21 (parametric analysis)")
print("=" * 65)

# If masses satisfy Koide with shift m_0, then:
# K(m_1+m0, m_2+m0, m_3+m0) = 2/3
# Given m_i = m_1 * (1, r21, r31) and K=2/3 determines r31(r21)
# then m0 is determined by the PDG values of r21, r31

# Key question: is m_0 / m_1 a universal function of r_21?
print(f"  m_0/m_1 for each sector:")
for name, (m1, m2, m3) in masses.items():
    m0 = find_m0(m1, m2, m3)
    r21 = m2/m1
    print(f"    {name}: m_0/m_1 = {m0/m1:10.4f},  r_21 = {r21:.2f}")

# ============================================================
# 7. KEY TEST: Koide-predicted r_31 from r_21 alone
# ============================================================
print("\n" + "=" * 65)
print("7. KOIDE-ONLY PREDICTION: r_31 from r_21 (no m_0 needed)")
print("=" * 65)
print("   For K(1, r_21, r_31) = 2/3, solve for r_31:")

def koide_from_ratios(r21, r31):
    return (1 + r21 + r31) / (1 + np.sqrt(r21) + np.sqrt(r31))**2

for name, (m1, m2, m3) in masses.items():
    r21 = m2/m1
    r31_pdg = m3/m1

    # Solve: K(1, r21, r31) = 2/3 for r31
    def obj(r31):
        return koide_from_ratios(r21, r31) - 2/3

    try:
        # There are two solutions; find the physical one (larger)
        r31_koide = brentq(obj, r21, 1e8)
        accuracy = abs(r31_koide - r31_pdg) / r31_pdg * 100
        print(f"\n  {name}:")
        print(f"    r_21 = {r21:.2f}")
        print(f"    r_31(Koide) = {r31_koide:.1f}")
        print(f"    r_31(PDG)   = {r31_pdg:.1f}")
        print(f"    Deviation: {accuracy:.1f}%")
        m3_pred = m1 * r31_koide
        print(f"    m_3 predicted: {m3_pred:.1f} MeV  (PDG: {m3:.1f} MeV)")
    except:
        print(f"\n  {name}: no solution found")

# ============================================================
# 8. SHIFTED KOIDE: r_31 from r_21 with shifted masses
# ============================================================
print("\n" + "=" * 65)
print("8. SHIFTED KOIDE: predict m_3 given (m_1, m_2, m_0)")
print("=" * 65)
print("   If we know m_0 independently, we can predict m_3:")

for name, (m1, m2, m3) in masses.items():
    m0 = find_m0(m1, m2, m3)

    # Predict m3 from (m1+m0, m2+m0) using Koide K=2/3
    s = np.sqrt(m1+m0) + np.sqrt(m2+m0)
    S = (m1+m0) + (m2+m0)
    disc = 16*s**2 - 4*(3*S - 2*s**2)
    x = (4*s + np.sqrt(disc)) / 2
    m3_sub_pred = x**2
    m3_pred = m3_sub_pred - m0

    print(f"\n  {name}: m_0 = {m0:.1f} MeV")
    print(f"    Predicted m_3 = {m3_pred:.1f} MeV  (PDG: {m3:.1f})")
    print(f"    This is circular (m_0 determined from same data)")
    print(f"    BUT: if m_0 came from INDEPENDENT source -> TRUE PREDICTION")

# ============================================================
# 9. FINAL: Summary table of predictions
# ============================================================
print("\n" + "=" * 65)
print("9. PREDICTION POWER SUMMARY")
print("=" * 65)
print("""
STRUCTURE OF TGP MASS PREDICTIONS:

For LEPTONS (0 free parameters beyond g0^e):
  Input:  g0^e (from B_tail=0 quantization) -> m_mu/m_e = 206.8
  Koide:  K=2/3 -> m_tau = 1776.9 MeV (0.002%)
  alpha_s: N_c^3 g0^e / (8 N_f^2) = 0.1174 (0.6 sigma)
  TOTAL: 3 predictions from 1 parameter

For QUARKS (1 additional parameter: m_0 per sector):
  Input:  g0^(d) from phi-FP + r_21(PDG) -> m_s/m_d = 20.0
  Shifted Koide: K(m+m0) = 2/3 -> m_b (if m_0 known independently)
  NEED:  m_0(down) and m_0(up) from first principles

  If m_0 ~ f(alpha_s, N_c, Q): ONE formula -> TWO m_0 values
  -> TWO predictions: m_b and m_t from (m_d, m_s) and (m_u, m_c)
""")
