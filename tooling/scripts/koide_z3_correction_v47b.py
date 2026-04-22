#!/usr/bin/env python3
"""
koide_z3_correction_v47b.py -- PROPOSAL #4: Higher-order Z3 corrections

QUESTION:
  The zero-order Z3 model (rigid 120-degree phase spacing) gives Q_K -> 1
  (or doesn't match 3/2). Proposal #4 asks whether CORRECTIONS to the
  Z3 model (from amplitude ratios A_i/A_j) restore Q_K = 3/2.

APPROACH:
  1. Extract actual phases delta(g0) for the three TGP solitons
  2. Compute the overlap energy E_int for a triplet with these phases
  3. Find the equilibrium phases by minimizing E_int
  4. Check if the equilibrium gives Q_K = 3/2

  The key idea: if three soliton tails overlap as
    h_i(r) = A_i * sin(r - delta_i) / r
  then the interaction energy between i and j goes as
    E_ij ~ A_i * A_j * cos(delta_i - delta_j) / d_ij
  (for large separation d_ij).

  For a BOUND STATE (single composite), d_ij -> 0, and the
  relevant quantity is the total field energy from superposition.

  The total tail field: h_total(r) = sum_i A_i * sin(r - delta_i) / r
  The overlap energy ~ integral h_i * h_j * r^2 dr involves
  cos(delta_i - delta_j) terms.

  At zero order (equal amplitudes, Z3 phases 0, 2pi/3, 4pi/3):
    sum cos(delta_i - delta_j) = -3/2 for all pairs.
  This minimizes overlap -> orthogonality.

  But A_e << A_mu << A_tau. The unequal amplitudes mean the
  "effective" phases that minimize E_int may differ from pure Z3.

  QUESTION: does the amplitude-weighted phase optimization
  naturally produce the phase structure that gives Q_K = 3/2?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0


def solver_full(g0, r_max=400):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0]


def extract_A_delta(g0, r_min=50, r_max=300):
    """Extract amplitude A and phase delta from tail."""
    r, g = solver_full(g0, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0, 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    A = np.sqrt(bc[0]**2 + bc[1]**2)
    delta = np.arctan2(bc[1], bc[0])
    return A, delta


def QK_masses(A_vals):
    A = np.array(A_vals)
    A2 = A**2
    A4 = A**4
    if np.sum(A4) < 1e-30:
        return 0.0
    return np.sum(A2)**2 / np.sum(A4)


# ================================================================
print("=" * 70)
print("PROPOSAL #4: HIGHER-ORDER Z3 PHASE CORRECTIONS")
print("=" * 70)

# ================================================================
print("\n" + "=" * 70)
print("SECTION 1: ACTUAL PHASES OF TGP SOLITONS")
print("=" * 70)

# Koide-enforced triplet
g0_e = 0.86770494
g0_mu = PHI * g0_e
# Find g0_tau for Koide
from scipy.optimize import brentq

A_e_full, delta_e = extract_A_delta(g0_e)
A_mu_full, delta_mu = extract_A_delta(g0_mu)

def find_g0_tau():
    def resid(g0t):
        At, _ = extract_A_delta(g0t)
        xe, xm, xt = A_e_full**2, A_mu_full**2, At**2
        return (xe+xm+xt)**2/(xe**2+xm**2+xt**2) - 1.5
    return brentq(resid, 1.50, 1.59)

g0_tau = find_g0_tau()
A_tau_full, delta_tau = extract_A_delta(g0_tau)

print(f"\n  Soliton tails:")
print(f"  {'':>6s} {'g0':>10s} {'A':>12s} {'delta (rad)':>12s} {'delta (deg)':>12s}")
print(f"  {'e':>6s} {g0_e:10.6f} {A_e_full:12.8f} {delta_e:12.6f} {np.degrees(delta_e):12.2f}")
print(f"  {'mu':>6s} {g0_mu:10.6f} {A_mu_full:12.8f} {delta_mu:12.6f} {np.degrees(delta_mu):12.2f}")
print(f"  {'tau':>6s} {g0_tau:10.6f} {A_tau_full:12.8f} {delta_tau:12.6f} {np.degrees(delta_tau):12.2f}")

# Phase differences
d_emu = delta_mu - delta_e
d_etau = delta_tau - delta_e
d_mutau = delta_tau - delta_mu

# Normalize to [0, 2pi)
def norm_angle(a):
    return a % (2*np.pi)

print(f"\n  Phase differences:")
print(f"  delta_mu - delta_e   = {norm_angle(d_emu):.4f} rad = {np.degrees(norm_angle(d_emu)):.2f} deg")
print(f"  delta_tau - delta_e  = {norm_angle(d_etau):.4f} rad = {np.degrees(norm_angle(d_etau)):.2f} deg")
print(f"  delta_tau - delta_mu = {norm_angle(d_mutau):.4f} rad = {np.degrees(norm_angle(d_mutau)):.2f} deg")
print(f"  Z3 would give: 120.00 deg each")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: OVERLAP ENERGY MODEL")
print("=" * 70)

# Model: three soliton tails h_i(r) = A_i * sin(r - delta_i) / r
# Total field deviation: H(r) = sum h_i(r)
# Energy density: rho ~ H^2
# Overlap terms: E_ij = integral A_i*A_j * sin(r-di)*sin(r-dj) / r^2 * r^2 dr
#              = A_i*A_j * integral sin(r-di)*sin(r-dj) dr
#              = A_i*A_j * L/2 * cos(di - dj)  (for integration length L >> 1)
# (The sin*sin product averages to cos(di-dj)/2 over many oscillations)

# So overlap energy is proportional to:
# E_overlap ~ sum_{i<j} A_i * A_j * cos(delta_i - delta_j)

# Total energy model:
# E_total = sum_i A_i^2/2  (self-energy)
#         + sum_{i<j} A_i * A_j * cos(delta_i - delta_j)  (overlap)

# Given FIXED amplitudes A_e, A_mu, A_tau, what phases MINIMIZE E_total?

A_vals = [A_e_full, A_mu_full, A_tau_full]

def E_overlap(deltas, A_vals):
    """Overlap energy for given phases and amplitudes."""
    d = np.array(deltas)
    A = np.array(A_vals)
    E = 0.0
    for i in range(3):
        for j in range(i+1, 3):
            E += A[i] * A[j] * np.cos(d[i] - d[j])
    return E

# Fix delta_e = 0 (gauge freedom), optimize delta_mu and delta_tau
def E_to_minimize(params):
    d_mu, d_tau = params
    return E_overlap([0.0, d_mu, d_tau], A_vals)

# Grid search first
print(f"\n  Overlap energy landscape (delta_e = 0 fixed):")
print(f"  {'delta_mu':>10s} {'delta_tau':>10s} {'E_overlap':>12s}")

best_E = 1e10
best_dm, best_dt = 0, 0
for dm in np.arange(0, 2*np.pi, 0.2):
    for dt in np.arange(0, 2*np.pi, 0.2):
        E = E_to_minimize([dm, dt])
        if E < best_E:
            best_E = E
            best_dm, best_dt = dm, dt

# Refine with optimizer
result = minimize(E_to_minimize, [best_dm, best_dt], method='Nelder-Mead')
opt_dm, opt_dt = result.x % (2*np.pi)

print(f"\n  Optimal phases (minimizing overlap energy):")
print(f"  delta_e   = 0.000 rad = 0.00 deg")
print(f"  delta_mu  = {opt_dm:.4f} rad = {np.degrees(opt_dm):.2f} deg")
print(f"  delta_tau = {opt_dt:.4f} rad = {np.degrees(opt_dt):.2f} deg")
print(f"  E_overlap = {result.fun:.8f}")

print(f"\n  Phase differences:")
print(f"  delta_mu - delta_e   = {np.degrees(opt_dm):.2f} deg")
print(f"  delta_tau - delta_e  = {np.degrees(opt_dt):.2f} deg")
print(f"  delta_tau - delta_mu = {np.degrees((opt_dt - opt_dm) % (2*np.pi)):.2f} deg")

# Compare with Z3 and actual
print(f"\n  Comparison:")
print(f"  {'':>20s} {'dm-de':>8s} {'dt-de':>8s} {'dt-dm':>8s}")
print(f"  {'Z3 (equal A)':>20s} {'120.00':>8s} {'240.00':>8s} {'120.00':>8s}")
print(f"  {'Optimal (real A)':>20s} {np.degrees(opt_dm):8.2f} {np.degrees(opt_dt):8.2f} {np.degrees((opt_dt-opt_dm)%(2*np.pi)):8.2f}")
print(f"  {'Actual ODE':>20s} {np.degrees(norm_angle(d_emu)):8.2f} {np.degrees(norm_angle(d_etau)):8.2f} {np.degrees(norm_angle(d_mutau)):8.2f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: DOES OPTIMAL PHASE STRUCTURE IMPLY Q_K = 3/2?")
print("=" * 70)

# The overlap energy model tells us the OPTIMAL phases for given amplitudes.
# But do the optimal phases have any relation to Koide?

# Key insight: the overlap energy is:
# E = A_e*A_mu*cos(d1) + A_e*A_tau*cos(d2) + A_mu*A_tau*cos(d2-d1)
# where d1 = delta_mu, d2 = delta_tau (delta_e = 0)

# For EQUAL amplitudes A: E = A^2 * sum cos(di-dj)
# Minimum: phases at 120 deg (Z3), E = -3A^2/2

# For UNEQUAL amplitudes: heavier pair wants to be anti-aligned (180 deg)
# Lighter particles have less influence.

# The actual Koide relation doesn't depend on phases directly.
# Q_K is about AMPLITUDES only. The phases are a separate degree of freedom.

# But: if the ODE COUPLES phases to amplitudes through the potential,
# then knowing the optimal phases constrains the amplitudes!

# Test: for the optimal phases, compute the "effective" Koide parameter
# that would result from an overlap-based mass formula.

# If mass ~ self_energy + overlap:
# m_i ~ A_i^2 + 2 * sum_{j!=i} A_j * cos(delta_i - delta_j)
# This is the "dressed mass" from overlap.

print(f"\n  'Dressed' masses from self + overlap energy:")
for i, (name, Ai) in enumerate(zip(['e', 'mu', 'tau'], A_vals)):
    phases = [0, opt_dm, opt_dt]
    self_E = Ai**2
    overlap_E = 0
    for j in range(3):
        if i != j:
            overlap_E += A_vals[j] * np.cos(phases[i] - phases[j])
    # The overlap contribution to mass
    m_bare = Ai**4
    m_overlap = Ai * overlap_E  # first order correction
    print(f"  {name:>4s}: A = {Ai:.6f}, A^4 = {m_bare:.8f}, overlap = {m_overlap:+.8f}, ratio = {m_overlap/m_bare:+.4f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: REVERSE APPROACH -- WHAT PHASES GIVE Q_K = 3/2?")
print("=" * 70)

# Instead of: given phases -> compute Q_K
# Try: given Q_K = 3/2 -> what phase structure is needed?

# The standard Koide formula doesn't involve phases.
# But if masses are MODIFIED by overlap, then:
# m_i = A_i^4 * (1 + epsilon_i)
# where epsilon_i depends on phases and amplitude ratios.

# For Q_K(dressed masses) = 3/2 with small epsilon:
# Q_K(A^4(1+eps)) ~ Q_K(A^4) + corrections
# Since we already enforce Q_K(A^4) = 3/2 for bare masses,
# the overlap corrections would BREAK Koide unless they vanish
# or cancel.

# Test: what is the overlap correction to Q_K?
print(f"\n  Overlap corrections to masses (using optimal phases):")
print(f"  Parametrize: m_i = A_i^4 * (1 + eps_i)")
print(f"  eps_i = [overlap energy at site i] / A_i^4")

for alpha in [0.01, 0.05, 0.1, 0.5]:
    phases = [0, opt_dm, opt_dt]
    m_dressed = []
    for i in range(3):
        m_bare = A_vals[i]**4
        overlap = 0
        for j in range(3):
            if i != j:
                overlap += A_vals[i] * A_vals[j] * np.cos(phases[i] - phases[j])
        m_dressed.append(m_bare + alpha * overlap)

    # Ensure all positive
    if all(m > 0 for m in m_dressed):
        qk = QK_masses([m**0.25 for m in m_dressed])  # back to A
        # Actually QK_masses expects A, computes A^4
        # Let me use the dressed masses directly
        m = np.array(m_dressed)
        qk = np.sum(np.sqrt(m))**2 / np.sum(m) if np.all(m > 0) else 0
        print(f"  alpha = {alpha:.2f}: m_dressed = [{m_dressed[0]:.6e}, {m_dressed[1]:.6e}, {m_dressed[2]:.6e}], Q_K = {qk:.6f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: ORTHOGONALITY CONDITION")
print("=" * 70)

# From ex151_Z3_orthogonality.py:
# The condition sum_{i<j} cos(delta_i - delta_j) = -3/2
# is the Z3 orthogonality condition.
# For unequal amplitudes, the WEIGHTED version is:
# sum_{i<j} A_i * A_j * cos(delta_i - delta_j) = 0
# This is "weighted orthogonality" or "zero total overlap energy."

# Check actual ODE phases
print(f"\n  Unweighted: sum cos(di-dj) = ", end="")
s = (np.cos(d_emu) + np.cos(d_etau) + np.cos(d_mutau))
print(f"{s:.6f}  (Z3 value: -1.500000)")

print(f"  Weighted:   sum A_i*A_j*cos(di-dj) = ", end="")
w = (A_e_full*A_mu_full*np.cos(d_emu)
   + A_e_full*A_tau_full*np.cos(d_etau)
   + A_mu_full*A_tau_full*np.cos(d_mutau))
print(f"{w:.6f}")

print(f"  Weighted (optimal): sum A_i*A_j*cos(di-dj) = {result.fun:.6f}")

# The optimal phases minimize overlap energy.
# Is zero overlap energy a natural condition?
# It would mean: the solitons are ORTHOGONAL in the field sense.

# Check: if we FORCE zero overlap (not minimum), what phases?
# E_overlap = A_e*A_mu*cos(d1) + A_e*A_tau*cos(d2) + A_mu*A_tau*cos(d2-d1) = 0
# This is one equation in two unknowns (d1, d2).
# There's a curve of solutions.

print(f"\n  Finding phases with ZERO overlap energy:")
from scipy.optimize import fsolve

def zero_overlap(params):
    d1, d2 = params
    return [A_e_full*A_mu_full*np.cos(d1) + A_e_full*A_tau_full*np.cos(d2)
            + A_mu_full*A_tau_full*np.cos(d2-d1),
            # Second condition: symmetrize in some way
            # Let's minimize |d1 - 2pi/3|^2 + |d2 - 4pi/3|^2 as regularization
            # Actually, let's just scan d1 and solve for d2
            0]  # dummy

# Scan: for each d1, find d2 such that E_overlap = 0
print(f"  {'d1 (deg)':>10s} {'d2 (deg)':>10s} {'d2-d1 (deg)':>12s} {'E_overlap':>12s}")
for d1 in np.arange(0.5, 2*np.pi, 0.3):
    # E = c1*cos(d2) + c2*cos(d2-d1) + c3 = 0
    # where c1 = A_e*A_tau, c2 = A_mu*A_tau, c3 = A_e*A_mu*cos(d1)
    c1 = A_e_full * A_tau_full
    c2 = A_mu_full * A_tau_full
    c3 = A_e_full * A_mu_full * np.cos(d1)
    # E = c1*cos(d2) + c2*(cos(d2)*cos(d1) + sin(d2)*sin(d1)) + c3
    #   = (c1 + c2*cos(d1))*cos(d2) + c2*sin(d1)*sin(d2) + c3
    R = np.sqrt((c1 + c2*np.cos(d1))**2 + (c2*np.sin(d1))**2)
    if R > abs(c3):
        phi0 = np.arctan2(c2*np.sin(d1), c1 + c2*np.cos(d1))
        d2 = np.arccos(-c3/R) + phi0
        d2 = d2 % (2*np.pi)
        E_check = E_overlap([0, d1, d2], A_vals)
        if abs(E_check) < 0.01:
            print(f"  {np.degrees(d1):10.2f} {np.degrees(d2):10.2f} {np.degrees((d2-d1)%(2*np.pi)):12.2f} {E_check:12.6f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: PHASE-AMPLITUDE COUPLING IN ODE")
print("=" * 70)

# The key question: does the ODE COUPLE delta(g0) to A(g0) in a way
# that the phase differences track the amplitude ratios?

# Compute delta(g0) for many g0 values
print(f"\n  Phase delta(g0) across the ODE:")
print(f"  {'g0':>8s} {'A':>12s} {'delta (deg)':>12s} {'d(delta)/dg0':>14s}")

g0_scan = np.arange(0.50, 1.58, 0.05)
deltas_scan = []
A_scan = []
g0_valid = []

for g0 in g0_scan:
    if abs(g0 - 1.0) < 0.03:
        continue
    A, d = extract_A_delta(g0)
    if A > 0.01:
        deltas_scan.append(d)
        A_scan.append(A)
        g0_valid.append(g0)

for i, (g0, A, d) in enumerate(zip(g0_valid, A_scan, deltas_scan)):
    dd = 0
    if i > 0:
        dd = (deltas_scan[i] - deltas_scan[i-1]) / (g0_valid[i] - g0_valid[i-1])
    print(f"  {g0:8.4f} {A:12.8f} {np.degrees(d):12.2f} {dd:14.4f}")

# Compute the phase slope a = d(delta)/d(ln g0) at the lepton points
print(f"\n  Phase slope a = d(delta)/d(ln g0) near lepton solitons:")
# Numerical derivative
for name, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau)]:
    eps = 0.001
    _, d1 = extract_A_delta(g0 - eps)
    _, d2 = extract_A_delta(g0 + eps)
    a = (d2 - d1) / (np.log(g0+eps) - np.log(g0-eps))
    # Z3 condition: a * ln(phi) = 2pi/3
    z3_ratio = a * np.log(PHI) / (2*np.pi/3)
    print(f"  {name:>4s}: g0 = {g0:.6f}, a = {a:.4f}, a*ln(phi) = {a*np.log(PHI):.4f}, ratio to 2pi/3 = {z3_ratio:.4f}")


# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS: PROPOSAL #4")
print("=" * 70)
print("""
  FINDINGS:

  1. ACTUAL PHASES from the ODE are NOT Z3 (not 120-degree spacing).
     The phase differences depend on g0 through the ODE solution.

  2. OPTIMAL PHASES (minimizing overlap energy with actual amplitudes)
     differ from both Z3 and actual ODE phases.
     The overlap energy minimum depends strongly on A_tau >> A_mu >> A_e.

  3. ZERO OVERLAP ENERGY (orthogonality) defines a curve in (d1,d2) space,
     not a unique point. Additional condition needed to select phases.

  4. PHASE-AMPLITUDE COUPLING: the slope a = d(delta)/d(ln g0) varies
     significantly across the corridor. The Z3 condition a*ln(phi) = 2pi/3
     is not generically satisfied.

  5. OVERLAP CORRECTIONS TO MASSES: even small overlap (alpha ~ 0.01-0.1)
     shifts Q_K away from 3/2, not toward it. This means:
     - Koide is a property of BARE masses (A^4), not dressed masses
     - OR the overlap exactly vanishes (orthogonality)
     - OR the overlap somehow preserves Q_K (requires fine-tuning)
""")
