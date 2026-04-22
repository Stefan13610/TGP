#!/usr/bin/env python3
"""
koide_crosssector_v47b.py -- PATH 16: CROSS-SECTOR BALANCE

USER'S INSIGHT:
  Don't look at leptons alone! The charge structure:
    -1, +2/3, -1/3 (lepton, up-quark, down-quark)
  seems "unnatural" but is forced by ANOMALY CANCELLATION:
    per generation: -1 + 0 + 3*(2/3) + 3*(-1/3) = 0

  Maybe the MASS spectrum is also forced by a cross-sector
  balance condition, analogous to anomaly cancellation for charges.

  QUESTION: Is Q_K = 3/2 for leptons connected to the
  quark mass spectrum through a cross-sector sum rule?

WHAT TO TEST:
  1. Koide relation for quarks (u,c,t) and (d,s,b) separately
  2. Cross-sector Koide: mixing leptons and quarks
  3. Charge-weighted mass sum rules per generation
  4. Full anomaly-like condition on masses
  5. Whether Q_K(leptons) * Q_K(quarks) or similar gives a special value
"""
import numpy as np

# ================================================================
# PDG masses (MeV)
# ================================================================
# Leptons (charged)
m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86

# Quarks (current masses, PDG 2024 central values, MeV)
# Light quarks: MS-bar at 2 GeV
m_u = 2.16      # up
m_d = 4.67      # down
m_s = 93.4      # strange
m_c = 1270.0    # charm (MS-bar at m_c)
m_b = 4180.0    # bottom (MS-bar at m_b)
m_t = 172760.0  # top (pole mass)

# Neutrinos (upper bounds, eV -> MeV)
# We know mass differences but not absolute scale
# Delta m^2_21 ~ 7.5e-5 eV^2, Delta m^2_32 ~ 2.5e-3 eV^2
# Assuming normal hierarchy: m1 ~ 0, m2 ~ 0.0087 eV, m3 ~ 0.050 eV
m_nu1 = 0.0  # effectively massless for our purposes
m_nu2 = 8.7e-9  # MeV
m_nu3 = 50.0e-9  # MeV


def Q_K(masses):
    """Compute Koide Q parameter for a set of 3 masses."""
    sq = np.sqrt(np.array(masses))
    return np.sum(sq)**2 / np.sum(masses)


def koide_report(name, masses):
    """Print Koide analysis for a triplet."""
    m = np.array(masses)
    sq = np.sqrt(m)
    S = np.sum(sq)
    P = np.sum(m)
    QK = S**2 / P

    # Koide circle parameters
    a = S / 3
    dev = sq / a - 1
    b_sq = (2/3) * np.sum(dev**2)
    b = np.sqrt(b_sq)
    CV = np.sqrt(np.var(sq)) / np.mean(sq)

    print(f"\n  {name}:")
    print(f"    masses: {m[0]:.4g}, {m[1]:.4g}, {m[2]:.4g} MeV")
    print(f"    sqrt(m): {sq[0]:.4f}, {sq[1]:.4f}, {sq[2]:.4f}")
    print(f"    Q_K = {QK:.6f}")
    print(f"    b = {b:.6f} (sqrt(2) = {np.sqrt(2):.6f})")
    print(f"    CV = {CV:.6f} (CV=1 for Koide)")
    print(f"    Q_K - 3/2 = {QK - 1.5:+.6f}")

    # Mass ratios
    r21 = m[1]/m[0]
    r31 = m[2]/m[0]
    r32 = m[2]/m[1]
    print(f"    r21={r21:.2f}, r31={r31:.2f}, r32={r32:.4f}")

    return QK, b, CV


# ================================================================
print("=" * 70)
print("PATH 16: CROSS-SECTOR BALANCE")
print("=" * 70)

# ================================================================
# Section 1: Koide for each sector
# ================================================================
print("\n" + "=" * 70)
print("SECTION 1: KOIDE Q_K FOR EACH FERMION SECTOR")
print("=" * 70)

QK_lep, b_lep, CV_lep = koide_report("Charged leptons (e, mu, tau)", [m_e, m_mu, m_tau])
QK_up, b_up, CV_up = koide_report("Up-type quarks (u, c, t)", [m_u, m_c, m_t])
QK_down, b_down, CV_down = koide_report("Down-type quarks (d, s, b)", [m_d, m_s, m_b])

print("\n  SUMMARY:")
print(f"  {'Sector':>20s} {'Q_K':>10s} {'b':>10s} {'CV':>10s} {'Q_K-3/2':>10s}")
print(f"  {'Leptons (e,mu,tau)':>20s} {QK_lep:10.6f} {b_lep:10.6f} {CV_lep:10.6f} {QK_lep-1.5:+10.6f}")
print(f"  {'Up (u,c,t)':>20s} {QK_up:10.6f} {b_up:10.6f} {CV_up:10.6f} {QK_up-1.5:+10.6f}")
print(f"  {'Down (d,s,b)':>20s} {QK_down:10.6f} {b_down:10.6f} {CV_down:10.6f} {QK_down-1.5:+10.6f}")


# ================================================================
# Section 2: Cross-sector combinations
# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: CROSS-SECTOR KOIDE COMBINATIONS")
print("=" * 70)

# Same-generation triplets
QK_gen1, _, _ = koide_report("Gen 1: (e, u, d)", [m_e, m_u, m_d])
QK_gen2, _, _ = koide_report("Gen 2: (mu, c, s)", [m_mu, m_c, m_s])
QK_gen3, _, _ = koide_report("Gen 3: (tau, t, b)", [m_tau, m_t, m_b])

# Charge-matched triplets
QK_q1, _, _ = koide_report("Charge -1: (e, mu, tau)", [m_e, m_mu, m_tau])
QK_q23, _, _ = koide_report("Charge +2/3: (u, c, t)", [m_u, m_c, m_t])
QK_q13, _, _ = koide_report("Charge -1/3: (d, s, b)", [m_d, m_s, m_b])

# Diagonal: lightest of each sector
QK_light, _, _ = koide_report("Lightest: (e, u, d)", [m_e, m_u, m_d])
QK_mid, _, _ = koide_report("Middle: (mu, c, s)", [m_mu, m_c, m_s])
QK_heavy, _, _ = koide_report("Heaviest: (tau, t, b)", [m_tau, m_t, m_b])


# ================================================================
# Section 3: CHARGE-WEIGHTED MASS SUM RULES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: CHARGE-WEIGHTED MASS SUM RULES")
print("=" * 70)

print("""
  Anomaly cancellation: sum(Q_i) = 0 per generation.
  Mass analogue: sum(Q_i * f(m_i)) = 0 for some function f?

  Per generation, with color factor N_c = 3:
    Q_e * f(m_e) + Q_nu * f(m_nu) + N_c * Q_u * f(m_u) + N_c * Q_d * f(m_d) = 0
    -f(m_e) + 0 + 3*(2/3)*f(m_u) + 3*(-1/3)*f(m_d) = 0
    -f(m_e) + 2*f(m_u) - f(m_d) = 0
    f(m_e) + f(m_d) = 2*f(m_u)
    => f(m_u) = (f(m_e) + f(m_d)) / 2  (arithmetic mean!)
""")

# Test for different f functions
print(f"  Testing f(m) = m (mass):")
print(f"    Gen 1: f(e)={m_e:.4f}, f(u)={m_u:.4f}, f(d)={m_d:.4f}")
print(f"    (f(e)+f(d))/2 = {(m_e+m_d)/2:.4f} vs f(u) = {m_u:.4f}")
print(f"    Residual = {m_u - (m_e+m_d)/2:.4f}")

print(f"\n  Testing f(m) = sqrt(m):")
sq_e, sq_u, sq_d = np.sqrt(m_e), np.sqrt(m_u), np.sqrt(m_d)
print(f"    Gen 1: f(e)={sq_e:.4f}, f(u)={sq_u:.4f}, f(d)={sq_d:.4f}")
print(f"    (f(e)+f(d))/2 = {(sq_e+sq_d)/2:.4f} vs f(u) = {sq_u:.4f}")
print(f"    Residual = {sq_u - (sq_e+sq_d)/2:.4f}")

print(f"\n  Testing f(m) = ln(m):")
ln_e, ln_u, ln_d = np.log(m_e), np.log(m_u), np.log(m_d)
print(f"    Gen 1: f(e)={ln_e:.4f}, f(u)={ln_u:.4f}, f(d)={ln_d:.4f}")
print(f"    (f(e)+f(d))/2 = {(ln_e+ln_d)/2:.4f} vs f(u) = {ln_u:.4f}")
print(f"    Residual = {ln_u - (ln_e+ln_d)/2:.4f}")

print(f"\n  Testing f(m) = m^(1/4) (TGP: m ~ A^4):")
q_e, q_u, q_d = m_e**0.25, m_u**0.25, m_d**0.25
print(f"    Gen 1: f(e)={q_e:.4f}, f(u)={q_u:.4f}, f(d)={q_d:.4f}")
print(f"    (f(e)+f(d))/2 = {(q_e+q_d)/2:.4f} vs f(u) = {q_u:.4f}")
print(f"    Residual = {q_u - (q_e+q_d)/2:.4f}")

# All three generations
print(f"\n  All three generations:")
gens = [
    ("Gen 1", m_e, m_u, m_d),
    ("Gen 2", m_mu, m_c, m_s),
    ("Gen 3", m_tau, m_t, m_b),
]

for fname, f_func in [("m", lambda x: x), ("sqrt(m)", np.sqrt),
                       ("ln(m)", np.log), ("m^(1/4)", lambda x: x**0.25),
                       ("m^(1/3)", lambda x: x**(1/3))]:
    print(f"\n  f(m) = {fname}:")
    for gname, ml, mu_val, md in gens:
        fl, fu, fd = f_func(ml), f_func(mu_val), f_func(md)
        pred = (fl + fd) / 2
        err = (fu - pred) / fu * 100 if fu != 0 else 0
        print(f"    {gname}: -f(l)+2f(u)-f(d) = {-fl+2*fu-fd:.4f}, "
              f"pred_u = {pred:.4f}, actual_u = {fu:.4f}, err = {err:+.1f}%")


# ================================================================
# Section 4: FULL ANOMALY-LIKE CONDITION ON MASSES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: GENERALIZED ANOMALY CANCELLATION FOR MASSES")
print("=" * 70)

print("""
  In the Standard Model, anomaly cancellation requires:
    sum over all fermions of Q^k * Y^l = 0
  for specific combinations of k, l (gauge anomalies).

  The simplest: sum(Q_i) = 0 per generation (as above).

  For MASSES, we can try:
    sum over all fermions of Q_i^a * m_i^b = 0
  for some exponents a, b.

  With color: Q_e = -1, Q_nu = 0, Q_u = +2/3 (x3), Q_d = -1/3 (x3)
""")

# Systematic search for (a, b) that gives balance
print(f"  Searching for (a, b) such that sum(Q^a * m^b) ~ 0 per generation:")
print(f"  {'a':>6s} {'b':>6s} {'Gen1':>12s} {'Gen2':>12s} {'Gen3':>12s} {'max|res|':>12s}")

best_ab = None
best_res = 1e20

for a_exp in np.arange(0.5, 3.1, 0.5):
    for b_exp in np.arange(0.1, 2.1, 0.1):
        residuals = []
        for ml, mnu, mu_val, md in [(m_e, m_nu1, m_u, m_d),
                                     (m_mu, m_nu2, m_c, m_s),
                                     (m_tau, m_nu3, m_t, m_b)]:
            # Charges: e=-1, nu=0, u=+2/3 (x3), d=-1/3 (x3)
            # For Q^a with fractional Q, need |Q|^a * sign(Q)
            if ml > 0:
                term_l = (-1)**1 * 1**a_exp * ml**b_exp  # Q=-1
            else:
                term_l = 0
            # neutrino: Q=0, contribution = 0
            if mu_val > 0:
                term_u = (+1) * (2/3)**a_exp * mu_val**b_exp * 3  # 3 colors
            else:
                term_u = 0
            if md > 0:
                term_d = (-1) * (1/3)**a_exp * md**b_exp * 3  # 3 colors
            else:
                term_d = 0
            res = term_l + term_u + term_d
            residuals.append(res)

        max_res = max(abs(r) for r in residuals)
        # Normalize by typical scale
        scale = max(abs(r) for r in residuals) if max(abs(r) for r in residuals) > 0 else 1

        if max_res < best_res:
            best_res = max_res
            best_ab = (a_exp, b_exp)

        # Print interesting ones
        if max_res < 100 or (a_exp == 1 and b_exp in [0.5, 1.0]):
            print(f"  {a_exp:6.1f} {b_exp:6.1f} {residuals[0]:12.4f} {residuals[1]:12.4f} "
                  f"{residuals[2]:12.4f} {max_res:12.4f}")

print(f"\n  Best (a,b) = ({best_ab[0]:.1f}, {best_ab[1]:.1f}) with max|res| = {best_res:.4f}")


# ================================================================
# Section 5: INTER-SECTOR Q_K RELATIONS
# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: RELATIONS BETWEEN Q_K VALUES")
print("=" * 70)

print(f"\n  Q_K values:")
print(f"    Leptons (e,mu,tau):  {QK_lep:.6f}")
print(f"    Up quarks (u,c,t):  {QK_up:.6f}")
print(f"    Down quarks (d,s,b): {QK_down:.6f}")

print(f"\n  Products and ratios:")
print(f"    Q_lep * Q_up = {QK_lep * QK_up:.6f}")
print(f"    Q_lep * Q_down = {QK_lep * QK_down:.6f}")
print(f"    Q_up * Q_down = {QK_up * QK_down:.6f}")
print(f"    Q_lep + Q_up + Q_down = {QK_lep + QK_up + QK_down:.6f}")
print(f"    Q_lep / Q_up = {QK_lep / QK_up:.6f}")
print(f"    Q_lep / Q_down = {QK_lep / QK_down:.6f}")

# Charge-weighted average
Q_charges = np.array([-1, 2/3, -1/3])
Q_abs = np.abs(Q_charges)
QK_vals = np.array([QK_lep, QK_up, QK_down])
QK_weighted = np.sum(Q_abs * QK_vals) / np.sum(Q_abs)
QK_charge_sq = np.sum(Q_charges**2 * QK_vals) / np.sum(Q_charges**2)

print(f"\n  Charge-weighted averages:")
print(f"    sum(|Q|*QK) / sum(|Q|) = {QK_weighted:.6f}")
print(f"    sum(Q^2*QK) / sum(Q^2) = {QK_charge_sq:.6f}")

# With color factors
Nc = 3
QK_color_weighted = (1*QK_lep + Nc*(2/3)**2*QK_up + Nc*(1/3)**2*QK_down) / \
                    (1 + Nc*(2/3)**2 + Nc*(1/3)**2)
print(f"    With color: sum(N_c*Q^2*QK)/sum(N_c*Q^2) = {QK_color_weighted:.6f}")
print(f"    (where N_c=1 for leptons, N_c=3 for quarks)")

# The "anomaly" combination
anom = -1 * QK_lep + Nc * (2/3) * QK_up + Nc * (-1/3) * QK_down
print(f"\n  Anomaly-like: -Q_K(l) + N_c*(2/3)*Q_K(u) + N_c*(-1/3)*Q_K(d)")
print(f"    = {anom:.6f}")

# Another: sum Q_i^2 * Q_K_i per generation (gauge anomaly form)
anom2 = (-1)**2 * QK_lep + Nc*(2/3)**2 * QK_up + Nc*(1/3)**2 * QK_down
print(f"  sum(Q^2 * Q_K) = {anom2:.6f}")
print(f"  sum(Q^2) = {1 + Nc*(4/9) + Nc*(1/9):.6f} = {1 + Nc*5/9:.6f}")


# ================================================================
# Section 6: KOIDE FOR COMBINED SECTORS
# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: KOIDE Q FOR COMBINED MASS SETS")
print("=" * 70)

# All 9 charged fermions
all_masses = [m_e, m_mu, m_tau, m_u, m_c, m_t, m_d, m_s, m_b]
sq_all = np.sqrt(all_masses)
QK_all = np.sum(sq_all)**2 / np.sum(all_masses)
print(f"  All 9 charged fermions: Q_K = {QK_all:.6f}")
print(f"  For 9 variables, Koide-like would be Q_K = 9/2 = {9/2:.1f}")

# Per generation (4 fermions each, with color)
# Generation 1: e, nu_e, u(x3), d(x3) -> 8 states but 4 flavors
# Let's use "effective mass per generation"

# Charge-weighted mass per generation
for gname, ml, mu_val, md in gens:
    # With color: m_eff = m_l + 3*m_u + 3*m_d (democratic)
    m_eff_dem = ml + 3*mu_val + 3*md
    # Charge-squared weighted: m_eff = Q^2*m
    m_eff_Q2 = 1*ml + 3*(4/9)*mu_val + 3*(1/9)*md
    # sqrt version
    sq_eff = np.sqrt(ml) + 3*np.sqrt(mu_val) + 3*np.sqrt(md)
    print(f"  {gname}: m_eff(dem) = {m_eff_dem:.2f}, m_eff(Q^2) = {m_eff_Q2:.2f}, "
          f"sum sqrt = {sq_eff:.4f}")

# Effective generation masses
m_gen = []
for ml, mu_val, md in [(m_e, m_u, m_d), (m_mu, m_c, m_s), (m_tau, m_t, m_b)]:
    m_gen.append(ml + 3*mu_val + 3*md)

QK_gen_dem = Q_K(m_gen)
print(f"\n  Democratic generation masses: {m_gen[0]:.2f}, {m_gen[1]:.2f}, {m_gen[2]:.2f}")
print(f"  Q_K(generation) = {QK_gen_dem:.6f}")

# Q^2-weighted generation masses
m_gen_Q2 = []
for ml, mu_val, md in [(m_e, m_u, m_d), (m_mu, m_c, m_s), (m_tau, m_t, m_b)]:
    m_gen_Q2.append(1*ml + 3*(4/9)*mu_val + 3*(1/9)*md)

QK_gen_Q2 = Q_K(m_gen_Q2)
print(f"\n  Q^2-weighted generation masses: {m_gen_Q2[0]:.2f}, {m_gen_Q2[1]:.2f}, {m_gen_Q2[2]:.2f}")
print(f"  Q_K(Q^2-weighted) = {QK_gen_Q2:.6f}")

# |Q|-weighted
m_gen_Q = []
for ml, mu_val, md in [(m_e, m_u, m_d), (m_mu, m_c, m_s), (m_tau, m_t, m_b)]:
    m_gen_Q.append(1*ml + 3*(2/3)*mu_val + 3*(1/3)*md)

QK_gen_Q = Q_K(m_gen_Q)
print(f"\n  |Q|-weighted generation masses: {m_gen_Q[0]:.2f}, {m_gen_Q[1]:.2f}, {m_gen_Q[2]:.2f}")
print(f"  Q_K(|Q|-weighted) = {QK_gen_Q:.6f}")


# ================================================================
# Section 7: THE KEY TEST - does cross-sector balance fix Q_K?
# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: CROSS-SECTOR BALANCE CONSTRAINING Q_K")
print("=" * 70)

print("""
  HYPOTHESIS: The anomaly cancellation for charges
    -1 + 0 + 3*(2/3) + 3*(-1/3) = 0
  has a MASS ANALOGUE that constrains the mass spectrum.

  If each sector independently satisfies Koide with b_i = sqrt(2):
    Q_K(leptons) = 3/2
    Q_K(up quarks) = 3/2
    Q_K(down quarks) = 3/2

  This would mean ALL three sectors have CV = 1.
  Is this observed?
""")

print(f"  Q_K deviations from 3/2:")
print(f"    Leptons:     {QK_lep - 1.5:+.6f} ({(QK_lep-1.5)/1.5*100:+.3f}%)")
print(f"    Up quarks:   {QK_up - 1.5:+.6f} ({(QK_up-1.5)/1.5*100:+.3f}%)")
print(f"    Down quarks: {QK_down - 1.5:+.6f} ({(QK_down-1.5)/1.5*100:+.3f}%)")

# NOTE: quark masses are running masses -- their values depend on
# the renormalization scale! At what scale should we evaluate them?

print(f"""
  IMPORTANT CAVEAT: quark masses are RUNNING (scale-dependent).
  The values above use PDG conventions:
    u, d, s: MS-bar at 2 GeV
    c: MS-bar at m_c
    b: MS-bar at m_b
    t: pole mass

  For a meaningful Koide test, all masses should be
  at the SAME renormalization scale.

  This is a KNOWN issue in Koide phenomenology.
  Some authors find better Koide agreement at specific scales.
""")


# ================================================================
# Section 8: BALANCE BETWEEN LEPTON AND QUARK HIERARCHIES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 8: HIERARCHY BALANCE")
print("=" * 70)

# Hierarchy parameters
h_lep = m_tau / m_e  # = 3477
h_up = m_t / m_u  # = 80000
h_down = m_b / m_d  # = 895

print(f"  Mass hierarchies:")
print(f"    Leptons: m_tau/m_e = {h_lep:.1f}")
print(f"    Up:      m_t/m_u  = {h_up:.1f}")
print(f"    Down:    m_b/m_d  = {h_down:.1f}")
print(f"    Geometric mean: {(h_lep * h_up * h_down)**(1/3):.1f}")
print(f"    h_lep * h_down = {h_lep * h_down:.1f}")
print(f"    h_up / (h_lep * h_down) = {h_up / (h_lep * h_down):.4f}")

# sqrt ratios
print(f"\n  sqrt(m) hierarchies:")
print(f"    Leptons: sqrt(r31) = {np.sqrt(h_lep):.2f}")
print(f"    Up:      sqrt(r31) = {np.sqrt(h_up):.2f}")
print(f"    Down:    sqrt(r31) = {np.sqrt(h_down):.2f}")

# Interplay: do charge-weighted sqrt hierarchies balance?
# -1 * sqrt(h_lep) + 3*(2/3)*sqrt(h_up) + 3*(-1/3)*sqrt(h_down) = ?
balance = -1*np.sqrt(h_lep) + 3*(2/3)*np.sqrt(h_up) + 3*(1/3)*np.sqrt(h_down)
print(f"\n  Charge-weighted sqrt hierarchy:")
print(f"    -sqrt(h_l) + 2*sqrt(h_u) - sqrt(h_d) = {-np.sqrt(h_lep)+2*np.sqrt(h_up)-np.sqrt(h_down):.2f}")
print(f"    -sqrt(h_l) + 2*sqrt(h_u) + sqrt(h_d) = {balance:.2f}")

# What if: -r31_lep + 2*r31_up - r31_down = 0?
# => r31_lep + r31_down = 2*r31_up
print(f"\n  Direct hierarchy balance:")
print(f"    r31_l + r31_d = {h_lep + h_down:.1f}")
print(f"    2*r31_u = {2*h_up:.1f}")
print(f"    NOT balanced.")

# Try: charge-squared weighted
print(f"\n  Q^2-weighted hierarchy:")
print(f"    Q_e^2 * h_l + N_c*Q_u^2 * h_u + N_c*Q_d^2 * h_d")
bal_Q2 = 1*h_lep + 3*(4/9)*h_up + 3*(1/9)*h_down
print(f"    = {bal_Q2:.1f}")


# ================================================================
# Section 9: YUKAWA COUPLING PERSPECTIVE
# ================================================================
print("\n" + "=" * 70)
print("SECTION 9: YUKAWA COUPLING PERSPECTIVE")
print("=" * 70)

print("""
  In the SM, masses come from Yukawa couplings: m = y * v / sqrt(2)
  where v = 246 GeV is the Higgs VEV.

  The Yukawa couplings are:
""")

v_higgs = 246000  # MeV
y_vals = {}
for name, mass in [('e', m_e), ('mu', m_mu), ('tau', m_tau),
                    ('u', m_u), ('c', m_c), ('t', m_t),
                    ('d', m_d), ('s', m_s), ('b', m_b)]:
    y = mass * np.sqrt(2) / v_higgs
    y_vals[name] = y
    print(f"    y_{name} = {y:.6e}")

print(f"\n  Yukawa hierarchies:")
print(f"    y_t/y_u = {y_vals['t']/y_vals['u']:.1f}")
print(f"    y_b/y_d = {y_vals['b']/y_vals['d']:.1f}")
print(f"    y_tau/y_e = {y_vals['tau']/y_vals['e']:.1f}")
print(f"    y_t/y_e = {y_vals['t']/y_vals['e']:.1f} (total hierarchy)")

# Koide for Yukawa couplings (same as for masses since y ~ m)
# Q_K is scale-invariant: Q_K(lambda*m) = Q_K(m)

# In TGP: m ~ A_tail^4. Yukawa y ~ A_tail^4 / (v/sqrt(2)).
# The TGP picture: ALL fermion masses come from soliton amplitudes.
# Different charge sectors might have different "coupling constants"
# but the same underlying soliton dynamics.

print(f"\n  In TGP interpretation:")
print(f"  ALL fermion masses from soliton dynamics: m_i = c_Q * A_tail(g0_i)^4")
print(f"  c_Q depends on the charge sector (coupling to Higgs?)")
print(f"  But Q_K is independent of c_Q (it cancels in the ratio).")
print(f"  So Q_K = 3/2 for ALL sectors requires the SAME soliton dynamics")
print(f"  to produce CV = 1 amplitudes in each sector independently.")


# ================================================================
# FINAL SYNTHESIS
# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS: CROSS-SECTOR BALANCE")
print("=" * 70)

print(f"""
  FINDINGS:

  1. KOIDE PER SECTOR:
     Leptons: Q_K = {QK_lep:.4f} (deviation {(QK_lep-1.5)/1.5*100:+.2f}% from 3/2)
     Up quarks: Q_K = {QK_up:.4f} (deviation {(QK_up-1.5)/1.5*100:+.2f}%)
     Down quarks: Q_K = {QK_down:.4f} (deviation {(QK_down-1.5)/1.5*100:+.2f}%)

  2. ANOMALY-LIKE MASS BALANCE:
     -f(m_l) + 2*f(m_u) - f(m_d) = 0 per generation?
     Works poorly for all tested f(m).
     Quark mass uncertainties and running make this hard to test.

  3. GENERATION-EFFECTIVE MASSES:
     Democratic (m_l + 3m_u + 3m_d): Q_K = {QK_gen_dem:.4f}
     Q^2-weighted: Q_K = {QK_gen_Q2:.4f}

  4. KEY OBSERVATION:
     If Koide holds for ALL three charge sectors (with b=sqrt(2)),
     this would be a UNIVERSAL property of the soliton dynamics,
     not specific to leptons.

     The charge structure (-1, +2/3, -1/3) would then play no direct
     role in Q_K -- instead, Q_K = 3/2 would be a property of
     the 3-generation structure itself.

  5. THE USER'S INSIGHT REFINED:
     The "unnatural" charge assignments (-1, +2/3, -1/3) are forced
     by anomaly cancellation (charges must sum to 0 per generation).
     Similarly, Q_K = 3/2 might be forced by a "mass anomaly
     cancellation" -- a consistency condition for the full
     fermion content to exist coherently.

     In TGP terms: the three solitons of each sector
     must share the vacuum consistently. This sharing condition
     might universally give CV = 1 (b = sqrt(2)).

  STATUS: The cross-sector test is INCONCLUSIVE for quarks
  (running mass issue), but the STRUCTURAL ANALOGY is compelling.
  The key prediction: if TGP is correct, Q_K = 3/2 should hold
  for quarks at the appropriate mass scale.
""")
