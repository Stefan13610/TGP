#!/usr/bin/env python3
"""
koide_joint_stability_v47b.py -- PATH 17: JOINT STABILITY OF FULL FERMION CONTENT

USER'S INSIGHT (corrected understanding):
  Don't analyze charges and masses separately!
  The FULL system (all particle types, all charges, all masses)
  must be stable SIMULTANEOUSLY. The stability condition
  constrains BOTH the charge assignments AND the mass ratios.

  The "unnatural" charges (-1, +2/3, -1/3) and Q_K = 3/2
  are BOTH outputs of the same stability principle.

APPROACH:
  1. Parametrize a general fermion content per generation:
     - Lepton: charge Q_l, mass m_l, N_c = 1
     - Up-type quark: charge Q_u, mass m_u, N_c = 3
     - Down-type quark: charge Q_d, mass m_d, N_c = 3
     - Neutrino: charge 0, mass ~ 0

  2. Physical constraints:
     (a) Anomaly cancellation: Q_l + 3*Q_u + 3*Q_d = 0
     (b) Charge quantization: Q_u - Q_d = 1 (W-boson coupling)
         => Q_u = Q_d + 1  => Q_l = -3*(2*Q_d + 1)
         For Q_d = -1/3: Q_l = -1, Q_u = 2/3 (SM!)
     (c) Hypercharge: Y = Q - T_3

  3. Stability functionals (ALL depending on masses AND charges jointly):
     S1: Vacuum polarization energy (QED)
     S2: Running coupling (asymptotic freedom threshold)
     S3: Hydrogen atom binding (requires e-p bound state)
     S4: Nuclear stability (requires residual strong force)
     S5: Combined TGP functional

  4. Scan over (masses, charges) parameter space simultaneously
     to find stability regions and check if Q_K = 3/2 is selected.
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.optimize import minimize, minimize_scalar

# Physical constants
ALPHA_EM = 1/137.036
ALPHA_S = 0.118
M_PLANCK = 1.22e19  # GeV
V_HIGGS = 246.0  # GeV

# Observed masses (GeV)
M_E_OBS = 0.000511
M_MU_OBS = 0.10566
M_TAU_OBS = 1.77686
M_U_OBS = 0.00216
M_C_OBS = 1.27
M_T_OBS = 172.76
M_D_OBS = 0.00467
M_S_OBS = 0.0934
M_B_OBS = 4.18


def Q_K(masses):
    """Koide Q for 3 masses."""
    sq = np.sqrt(np.abs(masses))
    S = np.sum(sq)
    P = np.sum(masses)
    if P == 0:
        return 3.0
    return S**2 / P


# ================================================================
print("=" * 70)
print("PATH 17: JOINT STABILITY OF FULL FERMION CONTENT")
print("=" * 70)


# ================================================================
# Section 1: CHARGE STRUCTURE FROM ANOMALY CANCELLATION
# ================================================================
print("\nSECTION 1: CHARGE STRUCTURE")
print("-" * 50)

print("""
  Per generation, with N_c = 3 colors:
  Anomaly cancellation (U(1)_EM):
    N_c * Q_u + N_c * Q_d + Q_l + Q_nu = 0
    3*Q_u + 3*Q_d + Q_l = 0  (Q_nu = 0)

  W-boson coupling (SU(2) doublets):
    Q_u - Q_d = 1  (up and down in same doublet)
    Q_nu - Q_l = 1  (same for leptons, but Q_nu = 0 => Q_l = -1... no)
    Actually: Q_l = Q_nu - 1 => Q_l = -1 (for Q_nu = 0)

  Wait: the charge quantization comes from:
    Q_u - Q_d = 1 (quark doublet)
    Q_nu - Q_l = 1 (lepton doublet)
    Q_nu = 0 => Q_l = -1

  Combined with anomaly: 3*Q_u + 3*Q_d - 1 = 0
    Q_u + Q_d = 1/3
    Q_u = Q_d + 1
    => 2*Q_d + 1 = 1/3 => Q_d = -1/3, Q_u = 2/3

  So charges are FULLY DETERMINED by:
    [1] SU(2) doublet structure (Q_u - Q_d = 1, Q_nu - Q_l = 1)
    [2] Anomaly cancellation (3*Q_u + 3*Q_d + Q_l = 0)
    [3] Q_nu = 0

  There is NO freedom in charges. They are (-1, 0, +2/3, -1/3).
  The "unnaturalness" is just SU(2) x U(1) structure + anomaly freedom.
""")

print("  Charges are fixed by gauge structure + anomaly cancellation.")
print("  The free parameters are ONLY the masses.")
print()

# So the question refines to:
# Given the FIXED charges (-1, 0, +2/3, -1/3),
# what mass configuration minimizes some stability functional?

# But the user's deeper point is:
# Maybe the ENTIRE structure (gauge group + charges + masses)
# is selected jointly. In TGP, the "gauge group" might emerge
# from the soliton dynamics, not be assumed.


# ================================================================
# Section 2: STABILITY FUNCTIONALS
# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: JOINT STABILITY FUNCTIONALS")
print("=" * 70)

print("""
  We define stability functionals that depend on ALL fermion
  parameters simultaneously. The stable universe is the one
  that optimizes these functionals.

  S1: VACUUM ENERGY (Coleman-Weinberg)
      E_vac = sum_i (-1)^(2s_i) * N_c,i * m_i^4 * ln(m_i^2/mu^2) / (64*pi^2)
      For fermions (s=1/2): contribution is NEGATIVE.
      Total vacuum energy should be minimized (most negative = most stable??)
      Or: total vacuum energy = 0 (cosmological constant problem).

  S2: RUNNING ALPHA_EM (Landau pole avoidance)
      alpha_EM(mu) diverges at Landau pole.
      Position depends on fermion content: sum N_c,i * Q_i^2.
      Farther Landau pole = more stable.

  S3: ATOMIC STABILITY
      Hydrogen atom exists if: alpha_EM * m_e / m_p < 1
      (Bohr radius > Compton wavelength)
      Chemistry requires: alpha_EM * (m_e/m_p)^(1/2) << 1

  S4: NUCLEAR STABILITY
      Nuclei exist if: m_d - m_u > 0 (neutron heavier than proton)
      But not too heavy: m_d - m_u < ~10 MeV (beta stability)

  S5: TGP SOLITON COEXISTENCE
      All solitons must have g0 < g0_crit = 8/5.
      The phi-FP rule spaces generations by phi.
      The "field balance" requires some joint condition.
""")


# ================================================================
# Section 3: THE JOINT FUNCTIONAL
# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: COMPOSITE STABILITY FUNCTIONAL")
print("=" * 70)

def stability_functional(masses_gen, charges, colors, verbose=False):
    """
    Compute joint stability score for a generation of fermions.

    masses_gen: [m_l, m_u, m_d] for one generation (GeV)
    charges: [Q_l, Q_u, Q_d]
    colors: [N_l, N_u, N_d] = [1, 3, 3]

    Returns: dict of stability measures.
    """
    ml, mu, md = masses_gen
    Ql, Qu, Qd = charges
    Nl, Nu, Nd = colors

    results = {}

    # S1: Vacuum energy contribution
    # E_vac ~ -sum N_c * m^4 * ln(m^2/v^2) / (64*pi^2)
    # (simplified; sign is for fermions)
    mu_scale = V_HIGGS  # renormalization scale
    E_vac = 0
    for m, N, Q in zip([ml, mu, md], [Nl, Nu, Nd], [Ql, Qu, Qd]):
        if m > 0:
            E_vac -= N * m**4 * np.log(m**2 / mu_scale**2) / (64 * np.pi**2)
    results['E_vac'] = E_vac

    # S2: QED beta function coefficient (per generation)
    # b_EM = -(4/3) * sum N_c,i * Q_i^2
    b_EM = -(4/3) * (Nl*Ql**2 + Nu*Qu**2 + Nd*Qd**2)
    results['b_EM'] = b_EM

    # Landau pole scale: mu_L = mu_0 * exp(3*pi / (|b_EM|*alpha_EM(mu_0)))
    # Higher = more stable
    if b_EM != 0:
        log_landau = 3*np.pi / (abs(b_EM) * ALPHA_EM)
        results['log_landau'] = log_landau
    else:
        results['log_landau'] = np.inf

    # S3: Atomic stability parameter
    # eta_atom = alpha_EM * m_l / (m_u + m_d)  (lepton/quark mass ratio)
    # Needs to be small for bound atoms
    m_nucleon_eff = 3*(mu + md)  # very rough
    if m_nucleon_eff > 0:
        eta_atom = ALPHA_EM * ml / m_nucleon_eff
    else:
        eta_atom = np.inf
    results['eta_atom'] = eta_atom

    # S4: Nuclear stability: need m_d > m_u (neutron heavier)
    # and not too much: (m_d - m_u) / (m_d + m_u) < threshold
    delta_q = (md - mu) / (md + mu) if (md + mu) > 0 else 0
    results['delta_q'] = delta_q

    # S5: Hierarchy parameters
    if ml > 0:
        h_lu = mu / ml  # up/lepton mass ratio
        h_ld = md / ml  # down/lepton mass ratio
    else:
        h_lu = h_ld = 0
    results['h_lu'] = h_lu
    results['h_ld'] = h_ld

    # S6: Koide Q for this generation's masses
    QK_gen = Q_K([ml, mu, md])
    results['QK_gen'] = QK_gen

    # S7: Charge-mass coupling parameter
    # C = sum |Q_i| * sqrt(m_i) (charge-weighted amplitude)
    C = Nl*abs(Ql)*np.sqrt(ml) + Nu*abs(Qu)*np.sqrt(mu) + Nd*abs(Qd)*np.sqrt(md)
    results['C_Qm'] = C

    # S8: "Anomaly" with mass: sum Q_i * m_i^(1/2)
    anom_mass = Nl*Ql*np.sqrt(ml) + Nu*Qu*np.sqrt(mu) + Nd*Qd*np.sqrt(md)
    results['anom_mass_sqrt'] = anom_mass

    # S9: sum Q_i^2 * m_i (related to vacuum polarization)
    vac_pol = Nl*Ql**2*ml + Nu*Qu**2*mu + Nd*Qd**2*md
    results['vac_pol'] = vac_pol

    if verbose:
        for k, v in results.items():
            print(f"    {k}: {v:.6g}")

    return results


# Standard Model charges
Q_SM = [-1, 2/3, -1/3]
N_SM = [1, 3, 3]

print("\n  Generation 1 (SM values):")
res1 = stability_functional([M_E_OBS, M_U_OBS, M_D_OBS], Q_SM, N_SM, verbose=True)

print("\n  Generation 2 (SM values):")
res2 = stability_functional([M_MU_OBS, M_C_OBS, M_S_OBS], Q_SM, N_SM, verbose=True)

print("\n  Generation 3 (SM values):")
res3 = stability_functional([M_TAU_OBS, M_T_OBS, M_B_OBS], Q_SM, N_SM, verbose=True)


# ================================================================
# Section 4: SCAN MASS SPACE WITH FIXED CHARGES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: SCAN MASS SPACE -- WHAT Q_K MINIMIZES COMBINED FUNCTIONAL?")
print("=" * 70)

# Fix m_e and m_mu (from phi-FP), vary m_tau.
# For each m_tau, compute:
# (a) Q_K for leptons
# (b) The "matching" quark masses assuming same generation structure
# (c) Joint stability functional

# Physical picture: for a given lepton mass hierarchy,
# the quark masses adjust to maintain stability.
# Does the stable point give Q_K = 3/2?

# Approach: assume the mass ratios within each sector are
# constrained by Q_K, but Q_K can vary.
# For a given Q_K, compute r31 from the Koide formula,
# then compute ALL fermion masses for all 3 generations.

u_val = np.sqrt(M_MU_OBS / M_E_OBS)  # sqrt of lepton mass ratio

def masses_from_QK(QK_val, m1, r21):
    """Given Q_K and m1, r21, compute m2 and m3."""
    u = np.sqrt(r21)
    # (1+u+v)^2 = QK*(1+u^2+v^2)
    a_coef = QK_val - 1
    b_coef = -2*(1+u)
    c_coef = (QK_val-1)*(1+u**2) - 2*u
    disc = b_coef**2 - 4*a_coef*c_coef
    if disc < 0 or a_coef == 0:
        return None
    v = (-b_coef + np.sqrt(disc)) / (2*a_coef)
    if v <= 0:
        return None
    m2 = r21 * m1
    m3 = v**2 * m1
    return [m1, m2, m3]


# For quarks, we need their r21 values too
r21_lep = M_MU_OBS / M_E_OBS    # 206.77
r21_up = M_C_OBS / M_U_OBS       # 587.96
r21_down = M_S_OBS / M_D_OBS     # 20.00

print(f"\n  Mass ratios (r21):")
print(f"    Leptons: {r21_lep:.2f}")
print(f"    Up quarks: {r21_up:.2f}")
print(f"    Down quarks: {r21_down:.2f}")

# Now scan Q_K for ALL sectors simultaneously
# Assume all three sectors share the SAME Q_K (universality hypothesis)

print(f"\n  Scan: ALL sectors with same Q_K (universality)")
print(f"  {'Q_K':>8s} {'m_tau':>8s} {'m_t':>10s} {'m_b':>8s} "
      f"{'E_vac_1':>10s} {'E_vac_3':>10s} {'anom_1':>10s} {'anom_3':>10s}")

for QK_scan in np.arange(1.05, 2.50, 0.05):
    m_lep = masses_from_QK(QK_scan, M_E_OBS, r21_lep)
    m_up = masses_from_QK(QK_scan, M_U_OBS, r21_up)
    m_down = masses_from_QK(QK_scan, M_D_OBS, r21_down)

    if m_lep is None or m_up is None or m_down is None:
        continue

    # Generation 3 masses
    m_tau_s = m_lep[2]
    m_t_s = m_up[2]
    m_b_s = m_down[2]

    # Generation 1 stability
    res1_s = stability_functional([M_E_OBS, M_U_OBS, M_D_OBS], Q_SM, N_SM)

    # Generation 3 stability (masses depend on Q_K)
    res3_s = stability_functional([m_tau_s, m_t_s, m_b_s], Q_SM, N_SM)

    marker = " <-- Koide" if abs(QK_scan - 1.5) < 0.03 else ""
    if abs(QK_scan - 1.5) < 0.03 or abs(res3_s['anom_mass_sqrt']) < 0.5:
        marker += " *"
    print(f"  {QK_scan:8.4f} {m_tau_s:8.4f} {m_t_s:10.2f} {m_b_s:8.2f} "
          f"{res1_s['E_vac']:10.4g} {res3_s['E_vac']:10.4g} "
          f"{res1_s['anom_mass_sqrt']:10.4f} {res3_s['anom_mass_sqrt']:10.4f}{marker}")


# ================================================================
# Section 5: THE KEY TEST -- CHARGE-MASS ANOMALY = 0
# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: CHARGE-MASS ANOMALY CANCELLATION")
print("=" * 70)

print("""
  The anomaly cancellation for CHARGES:
    Q_l + 3*Q_u + 3*Q_d = -1 + 2 - 1 = 0

  The MASS ANALOGUE would be:
    Q_l*f(m_l) + 3*Q_u*f(m_u) + 3*Q_d*f(m_d) = 0

  For f(m) = sqrt(m):
    -sqrt(m_l) + 2*sqrt(m_u) - sqrt(m_d) = 0

  This means: sqrt(m_u) = (sqrt(m_l) + sqrt(m_d)) / 2
  i.e., sqrt of up-quark mass = arithmetic mean of sqrt of lepton and down-quark masses!

  Does this hold PER GENERATION? And if imposed for all 3 generations,
  does it constrain Q_K?
""")

# For each generation, find what m_u would satisfy the balance
print("  Charge-mass balance: -sqrt(m_l) + 2*sqrt(m_u) - sqrt(m_d) = 0")
print("  => sqrt(m_u) = (sqrt(m_l) + sqrt(m_d)) / 2")
print()

for gname, ml, mu_obs, md in [("Gen 1", M_E_OBS, M_U_OBS, M_D_OBS),
                                ("Gen 2", M_MU_OBS, M_C_OBS, M_S_OBS),
                                ("Gen 3", M_TAU_OBS, M_T_OBS, M_B_OBS)]:
    sq_u_pred = (np.sqrt(ml) + np.sqrt(md)) / 2
    m_u_pred = sq_u_pred**2
    print(f"  {gname}: predicted m_u = {m_u_pred:.4f} GeV, actual = {mu_obs:.4f} GeV, "
          f"ratio = {m_u_pred/mu_obs:.4f}")

# Now: if this balance holds for ALL 3 generations,
# AND each sector has Koide Q_K, what value of Q_K is selected?

print(f"\n  If charge-mass balance holds for all generations:")
print(f"  AND each sector has same Q_K:")
print(f"  What Q_K is self-consistent?")

# The balance: sqrt(m_u^(i)) = (sqrt(m_l^(i)) + sqrt(m_d^(i))) / 2
# for i = 1, 2, 3 (generations)
#
# Combined with Koide for leptons: sqrt(m_l^(3)) = f(sqrt(m_l^(1)), sqrt(m_l^(2)), Q_K)
# Combined with Koide for down quarks: sqrt(m_d^(3)) = f(sqrt(m_d^(1)), sqrt(m_d^(2)), Q_K)
# Then: sqrt(m_u^(3)) = (sqrt(m_l^(3)) + sqrt(m_d^(3))) / 2
# And Koide for up quarks: sqrt(m_u^(3)) = f(sqrt(m_u^(1)), sqrt(m_u^(2)), Q_K)
#
# This gives an OVER-DETERMINED system. For a general Q_K,
# the balance + Koide conditions may be INCOMPATIBLE.
# There might be only a DISCRETE set of Q_K values where they're consistent!

print(f"\n  Self-consistency check:")
print(f"  For each Q_K, compute:")
print(f"    - m_tau from Koide (leptons)")
print(f"    - m_b from Koide (down quarks)")
print(f"    - m_t_predicted from balance: sqrt(m_t) = (sqrt(m_tau) + sqrt(m_b))/2")
print(f"    - m_t_Koide from Koide (up quarks)")
print(f"    - CONSISTENCY: m_t_predicted = m_t_Koide?")
print()

print(f"  {'Q_K':>8s} {'m_tau':>8s} {'m_b':>8s} {'m_t(bal)':>10s} {'m_t(Koide)':>12s} "
      f"{'ratio':>8s} {'diff%':>8s}")

best_QK = None
best_diff = 1e20

for QK_scan in np.arange(1.05, 2.80, 0.01):
    m_lep = masses_from_QK(QK_scan, M_E_OBS, r21_lep)
    m_up = masses_from_QK(QK_scan, M_U_OBS, r21_up)
    m_down = masses_from_QK(QK_scan, M_D_OBS, r21_down)

    if m_lep is None or m_up is None or m_down is None:
        continue

    m_tau_s = m_lep[2]
    m_b_s = m_down[2]
    m_t_koide = m_up[2]

    # Balance prediction for m_t:
    sq_t_bal = (np.sqrt(m_tau_s) + np.sqrt(m_b_s)) / 2
    m_t_bal = sq_t_bal**2

    ratio = m_t_bal / m_t_koide if m_t_koide > 0 else 0
    diff = abs(ratio - 1.0) * 100

    if diff < best_diff:
        best_diff = diff
        best_QK = QK_scan

    marker = ""
    if abs(QK_scan - 1.5) < 0.006:
        marker = " <-- Koide 3/2"
    if diff < 5:
        marker += " <-- CLOSE!"

    if abs(QK_scan - 1.5) < 0.02 or diff < 10 or abs(QK_scan - best_QK) < 0.02:
        print(f"  {QK_scan:8.4f} {m_tau_s:8.4f} {m_b_s:8.4f} {m_t_bal:10.2f} "
              f"{m_t_koide:12.2f} {ratio:8.4f} {diff:8.2f}%{marker}")

print(f"\n  BEST MATCH: Q_K = {best_QK:.4f} (ratio closest to 1.0, diff = {best_diff:.2f}%)")
print(f"  Is this 3/2? {'YES!' if abs(best_QK - 1.5) < 0.05 else 'NO -- Q_K = ' + f'{best_QK:.4f}'}")


# ================================================================
# Section 6: GENERALIZED BALANCE -- try different f(m)
# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: WHICH f(m) GIVES SELF-CONSISTENCY AT Q_K = 3/2?")
print("=" * 70)

# Instead of assuming f(m) = sqrt(m), search for the exponent p
# such that: -m_l^p + 2*m_u^p - m_d^p = 0 per generation
# AND Koide holds with Q_K = 3/2.

QK_target = 1.5
m_lep_K = masses_from_QK(QK_target, M_E_OBS, r21_lep)
m_up_K = masses_from_QK(QK_target, M_U_OBS, r21_up)
m_down_K = masses_from_QK(QK_target, M_D_OBS, r21_down)

if m_lep_K and m_up_K and m_down_K:
    print(f"\n  Masses at Q_K = 3/2:")
    print(f"    Leptons: {m_lep_K[0]:.6f}, {m_lep_K[1]:.6f}, {m_lep_K[2]:.4f} GeV")
    print(f"    Up:      {m_up_K[0]:.6f}, {m_up_K[1]:.4f}, {m_up_K[2]:.2f} GeV")
    print(f"    Down:    {m_down_K[0]:.6f}, {m_down_K[1]:.6f}, {m_down_K[2]:.4f} GeV")

    print(f"\n  Testing: -m_l^p + 2*m_u^p - m_d^p = 0 for each generation")
    print(f"  {'p':>8s} {'Gen1_res':>12s} {'Gen2_res':>12s} {'Gen3_res':>12s} {'max|res|':>12s}")

    for p in np.arange(0.05, 2.01, 0.05):
        residuals = []
        for i in range(3):
            ml_p = m_lep_K[i]**p
            mu_p = m_up_K[i]**p
            md_p = m_down_K[i]**p
            res = -ml_p + 2*mu_p - md_p
            # Normalize by scale
            scale = max(ml_p, mu_p, md_p)
            residuals.append(res / scale if scale > 0 else 0)

        max_res = max(abs(r) for r in residuals)
        marker = " <-- GOOD" if max_res < 0.3 else ""
        if max_res < 0.5 or abs(p - 0.5) < 0.03 or abs(p - 1.0) < 0.03:
            print(f"  {p:8.4f} {residuals[0]:12.6f} {residuals[1]:12.6f} "
                  f"{residuals[2]:12.6f} {max_res:12.6f}{marker}")


# ================================================================
# Section 7: THE DEEPEST TEST -- COMBINED ANOMALY + KOIDE
# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: COMBINED ANOMALY + KOIDE SELF-CONSISTENCY")
print("=" * 70)

print("""
  The ultimate question: is there a value of Q_K such that
  BOTH of the following hold simultaneously:

  (A) Koide Q_K is the SAME for all three charge sectors
  (B) The charge-mass balance is satisfied per generation:
      sum_flavors N_c * Q * m^p = 0

  This is an OVER-CONSTRAINED system:
  - 3 sectors x 1 parameter each = same Q_K (constraint A)
  - 3 generations x 1 equation each = balance (constraint B)
  - Free parameters: only Q_K and p

  Let's find the (Q_K, p) pair that minimizes the total residual.
""")

def total_residual(params):
    """Total residual for (Q_K, p) self-consistency."""
    QK_val, p = params
    if QK_val <= 1.01 or QK_val >= 2.99 or p <= 0.01 or p > 3.0:
        return 1e10

    m_l = masses_from_QK(QK_val, M_E_OBS, r21_lep)
    m_u = masses_from_QK(QK_val, M_U_OBS, r21_up)
    m_d = masses_from_QK(QK_val, M_D_OBS, r21_down)

    if m_l is None or m_u is None or m_d is None:
        return 1e10

    # Balance residual per generation
    total = 0
    for i in range(3):
        if m_l[i] <= 0 or m_u[i] <= 0 or m_d[i] <= 0:
            return 1e10
        ml_p = m_l[i]**p
        mu_p = m_u[i]**p
        md_p = m_d[i]**p
        scale = max(ml_p, mu_p, md_p)
        if scale <= 0:
            return 1e10
        res = (-ml_p + 2*mu_p - md_p) / scale
        total += res**2

    return total


# Grid search
print(f"\n  Grid search over (Q_K, p):")
print(f"  {'Q_K':>8s} {'p':>8s} {'residual':>12s}")

best_params = None
best_total = 1e20

for QK_g in np.arange(1.1, 2.5, 0.02):
    for p_g in np.arange(0.05, 2.0, 0.02):
        r = total_residual([QK_g, p_g])
        if r < best_total:
            best_total = r
            best_params = (QK_g, p_g)

print(f"\n  GRID BEST: Q_K = {best_params[0]:.4f}, p = {best_params[1]:.4f}, "
      f"residual = {best_total:.6f}")

# Refine with optimizer
from scipy.optimize import minimize
res_opt = minimize(total_residual, best_params, method='Nelder-Mead',
                   options={'xatol': 1e-6, 'fatol': 1e-10})

QK_opt, p_opt = res_opt.x
print(f"  OPTIMIZED: Q_K = {QK_opt:.6f}, p = {p_opt:.6f}, "
      f"residual = {res_opt.fun:.10f}")
print(f"  Q_K - 3/2 = {QK_opt - 1.5:+.6f}")
print(f"  Is Q_K ~ 3/2? {'YES!' if abs(QK_opt - 1.5) < 0.05 else 'NO'}")

# Verify at the optimum
if res_opt.fun < 1.0:
    m_l = masses_from_QK(QK_opt, M_E_OBS, r21_lep)
    m_u = masses_from_QK(QK_opt, M_U_OBS, r21_up)
    m_d = masses_from_QK(QK_opt, M_D_OBS, r21_down)
    if m_l and m_u and m_d:
        print(f"\n  At optimum (Q_K={QK_opt:.4f}, p={p_opt:.4f}):")
        for i, gname in enumerate(["Gen 1", "Gen 2", "Gen 3"]):
            ml_p = m_l[i]**p_opt
            mu_p = m_u[i]**p_opt
            md_p = m_d[i]**p_opt
            balance = -ml_p + 2*mu_p - md_p
            scale = max(ml_p, mu_p, md_p)
            print(f"    {gname}: balance/scale = {balance/scale:.6f}")
            print(f"      m_l={m_l[i]:.4g}, m_u={m_u[i]:.4g}, m_d={m_d[i]:.4g}")


# ================================================================
# Section 8: BROADER BALANCE -- not just -1 + 2 - 1
# ================================================================
print("\n" + "=" * 70)
print("SECTION 8: GENERALIZED BALANCE COEFFICIENTS")
print("=" * 70)

print("""
  Maybe the balance isn't exactly -1 + 2 - 1.
  Let the coefficients be determined by the actual charges and colors:
    c_l * m_l^p + c_u * m_u^p + c_d * m_d^p = 0

  where c_l, c_u, c_d are functions of (Q, N_c).

  Candidate coefficient sets:
    (a) Anomaly: c = N_c * Q = (-1, +2, -1)
    (b) Q^2 weighted: c = N_c * Q^2 = (1, 4/3, 1/3)
    (c) Hypercharge: c = N_c * Y
    (d) Sign of Q only: c = N_c * sign(Q) = (-1, +3, -3)
    (e) Just Q: c = Q = (-1, 2/3, -1/3) [no color]
""")

coefficient_sets = {
    'N_c*Q':       [-1, 2, -1],
    'N_c*Q^2':     [1, 4/3, 1/3],
    'Q (no color)': [-1, 2/3, -1/3],
    'N_c*|Q|':     [1, 2, 1],
    'N_c*sign(Q)': [-1, 3, -3],
}

for cname, coeffs in coefficient_sets.items():
    def total_res_c(params, c=coeffs):
        QK_val, p = params
        if QK_val <= 1.01 or QK_val >= 2.99 or p <= 0.01 or p > 3.0:
            return 1e10
        m_l = masses_from_QK(QK_val, M_E_OBS, r21_lep)
        m_u = masses_from_QK(QK_val, M_U_OBS, r21_up)
        m_d = masses_from_QK(QK_val, M_D_OBS, r21_down)
        if m_l is None or m_u is None or m_d is None:
            return 1e10
        total = 0
        for i in range(3):
            if m_l[i] <= 0 or m_u[i] <= 0 or m_d[i] <= 0:
                return 1e10
            val = c[0]*m_l[i]**p + c[1]*m_u[i]**p + c[2]*m_d[i]**p
            scale = max(abs(c[0])*m_l[i]**p, abs(c[1])*m_u[i]**p, abs(c[2])*m_d[i]**p)
            if scale > 0:
                total += (val/scale)**2
        return total

    # Grid search
    best_c = (1.5, 0.5)
    best_r = 1e20
    for QK_g in np.arange(1.1, 2.5, 0.05):
        for p_g in np.arange(0.05, 2.0, 0.05):
            r = total_res_c([QK_g, p_g])
            if r < best_r:
                best_r = r
                best_c = (QK_g, p_g)

    # Refine
    res_c = minimize(total_res_c, best_c, method='Nelder-Mead',
                     options={'xatol': 1e-6, 'fatol': 1e-10})
    QK_c, p_c = res_c.x

    marker = " <-- KOIDE!" if abs(QK_c - 1.5) < 0.05 else ""
    print(f"  {cname:>15s}: Q_K = {QK_c:.4f}, p = {p_c:.4f}, "
          f"residual = {res_c.fun:.6f}{marker}")


# ================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
  PATH 17 RESULTS:

  The user's hypothesis: stability of the FULL fermion system
  (charges + masses together) constrains Q_K.

  We tested: combined Koide + charge-mass balance:
    Koide: all sectors have same Q_K
    Balance: sum c_i * m_i^p = 0 per generation

  BEST FIT: Q_K = {QK_opt:.4f}, p = {p_opt:.4f} (with N_c*Q coefficients)

  The key finding depends on whether Q_K ~ 3/2 emerges
  from the self-consistency condition. If it does, this would
  be the DERIVATION of Koide from cross-sector balance --
  exactly the user's hypothesis.
""")
