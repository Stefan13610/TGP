#!/usr/bin/env python3
"""
r31_universality_proof_v47b.py -- WHY is r31 independent of alpha?

OBSERVATION (F_alpha_canonical_v47b.py):
  For ALL alpha in [0.5, 3.0], Koide Q_K = 3/2 gives r31 = 3477.4.

HYPOTHESIS: This is a PURELY ALGEBRAIC consequence of:
  1. Koide Q_K = (sum sqrt(m))^2 / (sum m) = 3/2
  2. r_21 = m_mu/m_e = 206.768  (calibration input)
  3. phi-FP: g0_mu = phi * g0_e  (selection rule)

If r31 is determined by these three conditions ALONE
(without reference to alpha or the ODE), then r31 universality
is automatic.

TEST: Solve the Koide relation algebraically given r_21.
"""
import numpy as np
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2

# Physical constants
m_e = 0.51099895   # MeV
m_mu = 105.6583755  # MeV
m_tau = 1776.86     # MeV
r_21 = m_mu / m_e   # 206.768
r_31_pdg = m_tau / m_e  # 3477.18

print("=" * 70)
print("WHY IS r31 INDEPENDENT OF alpha?")
print("=" * 70)

# ================================================================
# Section 1: Pure Koide algebra
# ================================================================
print("\n1. PURE KOIDE ALGEBRA")
print("-" * 50)

# Koide relation: (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2
#                / (m_e + m_mu + m_tau) = 3/2
#
# Let x = sqrt(m_e), y = sqrt(m_mu), z = sqrt(m_tau)
# Then: (x + y + z)^2 / (x^2 + y^2 + z^2) = 3/2
#
# Equivalently: 2(x+y+z)^2 = 3(x^2+y^2+z^2)
# Expand: 2(x^2+y^2+z^2+2xy+2xz+2yz) = 3(x^2+y^2+z^2)
# => 4(xy+xz+yz) = x^2+y^2+z^2
# => x^2+y^2+z^2 - 4xy - 4xz - 4yz = 0
#
# Divide by x^2: 1 + (y/x)^2 + (z/x)^2 - 4(y/x) - 4(z/x) - 4(y/x)(z/x) = 0
#
# Let u = y/x = sqrt(r_21), v = z/x = sqrt(r_31)
# Then: 1 + u^2 + v^2 - 4u - 4v - 4uv = 0
# => v^2 - 4v(1+u) + (1 + u^2 - 4u) = 0
# Quadratic in v:
# v = [4(1+u) +/- sqrt(16(1+u)^2 - 4(1+u^2-4u))] / 2
# v = 2(1+u) +/- sqrt(4(1+u)^2 - (1+u^2-4u))
# Discriminant: 4(1+2u+u^2) - 1 - u^2 + 4u = 3 + 12u + 3u^2 = 3(1+4u+u^2)
# v = 2(1+u) +/- sqrt(3(1+4u+u^2))
# = 2(1+u) +/- sqrt(3)*sqrt((1+u)^2 + 2u)
# Hmm, let me simplify: 1+4u+u^2 = (1+u)^2 + 2u = (u+2)^2 - 3

u = np.sqrt(r_21)
print(f"  u = sqrt(r_21) = {u:.8f}")

disc = 3 * (1 + 4*u + u**2)
print(f"  Discriminant: 3*(1+4u+u^2) = {disc:.6f}")

v_plus = 2*(1+u) + np.sqrt(disc)
v_minus = 2*(1+u) - np.sqrt(disc)

r31_plus = v_plus**2
r31_minus = v_minus**2

print(f"\n  v+ = {v_plus:.8f}, r31+ = v+^2 = {r31_plus:.4f}")
print(f"  v- = {v_minus:.8f}, r31- = v-^2 = {r31_minus:.4f}")
print(f"  PDG: r31 = {r_31_pdg:.2f}")
print(f"  Match: r31+ vs PDG: {abs(r31_plus - r_31_pdg)/r_31_pdg*100:.4f}%")

# ================================================================
# Section 2: Does it match the ODE result?
# ================================================================
print("\n\n2. COMPARISON WITH ODE RESULT")
print("-" * 50)

r31_ode = 3477.44  # from F_alpha_canonical_v47b.py
print(f"  r31 (Koide algebra): {r31_plus:.4f}")
print(f"  r31 (ODE, any alpha): {r31_ode:.4f}")
print(f"  r31 (PDG):            {r_31_pdg:.2f}")
print(f"  Algebra vs ODE diff:  {abs(r31_plus - r31_ode):.4f} ({abs(r31_plus-r31_ode)/r31_ode*100:.4f}%)")

# ================================================================
# Section 3: The formula
# ================================================================
print("\n\n3. EXACT FORMULA FOR r31")
print("-" * 50)

print(f"""
  Given r_21 = m_mu/m_e, the Koide relation Q_K = 3/2 gives:

    sqrt(r_31) = 2(1 + sqrt(r_21)) + sqrt(3(1 + 4*sqrt(r_21) + r_21))

  This is a CLOSED-FORM algebraic expression.
  It depends ONLY on r_21, not on alpha, not on the ODE.

  Numerical check:
    r_21 = {r_21:.6f}
    sqrt(r_21) = {u:.8f}
    sqrt(r_31) = {v_plus:.8f}
    r_31 = {r31_plus:.6f}
    PDG  = {r_31_pdg:.2f}
    ODE  = {r31_ode:.2f}
""")

# ================================================================
# Section 4: Sensitivity analysis
# ================================================================
print("4. SENSITIVITY dr31/dr21")
print("-" * 50)

# From the formula: v(u) = 2(1+u) + sqrt(3(1+4u+u^2))
# dv/du = 2 + (4+2u)/(2*sqrt(3(1+4u+u^2))) * sqrt(3)
#       = 2 + sqrt(3)*(2+u)/sqrt(1+4u+u^2)
# dr31/dr21 = 2v*dv/du / (2u*du/du)... no, simpler:
# r31 = v^2, r21 = u^2
# dr31/dr21 = d(v^2)/d(u^2) = (2v*dv/du) / (2u) = v/u * dv/du

dvdu = 2 + np.sqrt(3) * (2 + u) / np.sqrt(1 + 4*u + u**2)
dr31_dr21 = v_plus / u * dvdu

print(f"  dv/du = {dvdu:.6f}")
print(f"  dr31/dr21 = {dr31_dr21:.4f}")
print(f"  For 0.01% change in r21: dr31/r31 = {dr31_dr21*r_21/r31_plus*0.01:.4f}%")

# ================================================================
# Section 5: Verify numerically for different r21 values
# ================================================================
print("\n\n5. r31 vs r21 (algebraic)")
print("-" * 50)

print(f"  {'r21':>10s} {'sqrt(r21)':>10s} {'r31':>10s} {'r31/r21':>10s}")
for r21_test in [100, 150, 200, 206.768, 250, 300, 500]:
    u_t = np.sqrt(r21_test)
    disc_t = 3 * (1 + 4*u_t + u_t**2)
    v_t = 2*(1+u_t) + np.sqrt(disc_t)
    r31_t = v_t**2
    print(f"  {r21_test:10.3f} {u_t:10.4f} {r31_t:10.2f} {r31_t/r21_test:10.4f}")


# ================================================================
# Section 6: Connection to ODE
# ================================================================
print("\n\n6. WHY THE ODE RESPECTS THIS")
print("-" * 50)

print("""
  The ODE produces A_tail(g0) for each alpha.
  We calibrate: (A_mu/A_e)^4 = r_21, where A_i = A_tail(g0_i).
  Then find g0_tau from Koide: (sum A^2)^2 / (sum A^4) = 3/2.
  Then compute r31 = (A_tau/A_e)^4.

  But the Koide relation is ALGEBRAIC in (A_e, A_mu, A_tau).
  Given A_mu/A_e = r_21^(1/4), the Koide relation determines
  A_tau/A_e = r31^(1/4) UNIQUELY.

  The A-values come from the ODE, but their RATIO r31/r21
  is determined by ALGEBRA, not by ODE dynamics.

  CONCLUSION: r31 universality is TAUTOLOGICAL.
  The Koide relation + r_21 input ALGEBRAICALLY determines r31.
  The ODE dynamics determine g0_tau, but r31 is fixed by algebra.

  This explains the F_alpha_canonical_v47b.py result:
  r31 = 3477.4 for all alpha, because the Koide algebra
  doesn't know about alpha.
""")

# ================================================================
# Section 7: What DOES depend on alpha?
# ================================================================
print("7. WHAT DEPENDS ON alpha")
print("-" * 50)

print("""
  alpha-DEPENDENT quantities:
  - g0_e(alpha):  electron initial value (calibrated)
  - g0_mu(alpha): muon initial value (= phi * g0_e)
  - g0_tau(alpha): tau initial value (from Koide + ODE dynamics)
  - A_e, A_mu, A_tau: individual amplitudes
  - margin g0_crit - g0_tau: how close tau is to collapse

  alpha-INDEPENDENT quantities:
  - r_21 = 206.77: INPUT (calibration)
  - r_31 = 3477.4: OUTPUT (Koide algebra)
  - Q_K = 3/2: IMPOSED (Koide relation)
  - A_mu/A_e = r_21^(1/4): follows from calibration
  - A_tau/A_e = r_31^(1/4): follows from Koide algebra

  The alpha-dependence is hidden in the INDIVIDUAL amplitudes
  and in g0_tau, but the MASS RATIOS are universal.
""")


# ================================================================
# Section 8: The other solution v_minus
# ================================================================
print("8. THE OTHER KOIDE SOLUTION")
print("-" * 50)

print(f"  v- = sqrt(r31_minus) = {v_minus:.8f}")
print(f"  r31- = {r31_minus:.6f}")
print(f"  m_tau- = {r31_minus * m_e:.2f} MeV")
print(f"  This would give m_tau = {r31_minus * m_e:.2f} MeV (vs PDG {m_tau} MeV)")

# Is v_minus positive?
if v_minus > 0:
    print(f"  v- > 0: second solution is PHYSICAL (lighter tau)")
    print(f"  r31-/r21 = {r31_minus/r_21:.4f}")
else:
    print(f"  v- < 0: second solution is UNPHYSICAL")


# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  THEOREM: Given r_21, the Koide relation Q_K = 3/2 determines r_31
  algebraically:

    r_31 = [2(1 + sqrt(r_21)) + sqrt(3(1 + 4*sqrt(r_21) + r_21))]^2

  Numerical value: r_31 = {r31_plus:.4f} (PDG: {r_31_pdg:.2f}, diff {abs(r31_plus-r_31_pdg)/r_31_pdg*100:.3f}%)

  This formula is INDEPENDENT of:
  - alpha (ODE kinetic parameter)
  - the specific ODE dynamics
  - the form (A vs B) of the equation
  - the phi-FP spacing rule

  It requires ONLY r_21 as input.

  The r_31 universality observed in F_alpha_canonical_v47b.py
  is therefore a TAUTOLOGY: it's pure algebra, not ODE dynamics.

  What the ODE DOES determine: g0_tau(alpha), i.e., which initial
  value produces the Koide-compatible tau soliton. This depends
  on alpha. But the resulting mass ratio r_31 does not.
""")
