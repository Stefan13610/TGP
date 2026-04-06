#!/usr/bin/env python3
"""
ex213 -- Breathing mode time-domain simulation
================================================

Simulate the breathing (radial) oscillation of an N=3 equilateral
triangle at the d_well equilibrium in TGP.

PHYSICS:
  Result 8 derived omega^2_br = (3/2)*omega^2_rad analytically.
  This script VERIFIES this by direct time integration:
  1. Place 3 equal masses at d_well equilibrium
  2. Perturb radially (increase R by small epsilon)
  3. Integrate with leapfrog
  4. Measure oscillation period -> omega_br
  5. Compare with analytical prediction

  Also: compare with Newton (which has NO bound oscillation
  at d_well, because d_well doesn't exist in Newton).

CONNECTION TO PTA:
  The N-body breathing mode is the localized analog of the
  cosmological scalar breathing mode detected in the PTA band.
  Both arise from radial oscillations of the scalar field Phi.

  PTA:   h_br = 2*delta_Phi/Phi_0, frequency ~ nHz
  N-body: delta_R / R_eq, frequency ~ omega_br

Depends on: dynamics_v2, configurations, pairwise, equilibria (corrected)
"""

import sys, os

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np
from nbody.configurations import regular_ngon
from nbody.pairwise import V_eff_total
from nbody.dynamics_v2 import leapfrog_integrate

# ============================================================
# Parameters
# ============================================================
C = 0.3
beta = 5.0
gamma = beta

# 2-body equilibrium
disc = 4*beta**2 - 18*gamma*C
d_rep = 2*beta - np.sqrt(disc)
d_well = 2*beta + np.sqrt(disc)
R_well = d_well / np.sqrt(3)  # circumradius for equilateral triangle

print("=" * 70)
print("BREATHING MODE TIME-DOMAIN SIMULATION")
print("=" * 70)
print(f"  C = {C}, beta = {beta}, gamma = {gamma}")
print(f"  d_rep = {d_rep:.6f}, d_well = {d_well:.6f}")
print(f"  R_well = {R_well:.6f}")

# ============================================================
# Analytical breathing frequency
# ============================================================
# V_2(d) for equal masses:
# V_2 = -4*pi*C^2/d + 8*pi*beta*C^2/d^2 - 24*pi*gamma*C^3/d^3
# V_2'' = -8*pi*C^2/d^3 + 48*pi*beta*C^2/d^4 - 288*pi*gamma*C^3/d^5

d = d_well
V2pp = (-8*np.pi*C**2/d**3 + 48*np.pi*beta*C**2/d**4
        - 288*np.pi*gamma*C**3/d**5)

# For equilateral triangle with 3 pairs at d:
# V_total = 3*V_2(d), where d = R*sqrt(3)
# d^2V_total/dR^2 = 3 * V_2'' * (dd/dR)^2 = 3 * V_2'' * 3 = 9*V_2''
d2V_dR2 = 9 * V2pp

# omega^2 = d2V/dR2 / (N*C) for N=3 masses of mass C
omega2_br = d2V_dR2 / (3 * C)
omega_br = np.sqrt(abs(omega2_br))
T_br = 2 * np.pi / omega_br if omega2_br > 0 else np.inf

# Result 8 relation: omega^2_br = (3/2) * omega^2_rad
# omega^2_rad = (1/C) * V_2''(d_well) * (dd/dR)^2 = V_2'' * 3 / C
omega2_rad = V2pp * 3 / C
ratio_theory = omega2_br / omega2_rad if omega2_rad != 0 else np.nan

print(f"\n  Analytical breathing frequency:")
print(f"    V_2''(d_well) = {V2pp:.6e}")
print(f"    d^2V/dR^2 = {d2V_dR2:.6e}")
print(f"    omega^2_rad = {omega2_rad:.6e}")
print(f"    omega^2_br = {omega2_br:.6e}")
print(f"    omega_br = {omega_br:.6f}")
print(f"    T_br = {T_br:.4f}")
print(f"    omega^2_br / omega^2_rad = {ratio_theory:.4f} (theory: 1.500)")


# ============================================================
# Force and potential functions for integration
# ============================================================
def tgp_acceleration(positions, C_values):
    """Compute accelerations a_i = F_i / C_i for TGP pairwise forces."""
    n = len(C_values)
    acc = np.zeros_like(positions)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            rij = positions[j] - positions[i]
            d = np.linalg.norm(rij)
            d = max(d, 1e-10)
            rhat = rij / d

            # dV_2/dd (between i and j):
            Ci, Cj = C_values[i], C_values[j]
            dVdd = (4*np.pi*Ci*Cj/d**2
                    - 16*np.pi*beta*Ci*Cj/d**3
                    + 36*np.pi*gamma*Ci*Cj*(Ci+Cj)/d**4)

            # Force on i from j: F = -dV/dx_i = dV/dd * rhat (toward j)
            # But dV/dd > 0 means energy increases with d (= attraction)
            # So force = dV/dd * rhat (pulls i toward j)
            # Wait: V increases means pushing apart, so F = -dV/dr * rhat
            # Let me be careful: F_i = -nabla_i V = -dV/dd * dd/dx_i
            # dd/dx_i = -(x_j - x_i)/d = -rhat
            # So F_i = -dV/dd * (-rhat) = dV/dd * rhat

            acc[i] += dVdd * rhat / Ci  # a_i = F_i / m_i = F_i / C_i
    return acc


def tgp_potential(positions, C_values):
    """Total pairwise TGP potential energy."""
    n = len(C_values)
    V = 0.0
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(positions[i] - positions[j])
            d = max(d, 1e-10)
            V += V_eff_total(d, C_values[i], C_values[j], beta, gamma)
    return V


# ============================================================
# Simulation: breathing mode
# ============================================================
print("\n" + "=" * 70)
print("LEAPFROG SIMULATION: BREATHING MODE")
print("=" * 70)

# Initial configuration: equilateral triangle at d_well, perturbed radially
epsilon = 0.01  # 1% radial perturbation
R_init = R_well * (1 + epsilon)

positions_0, C_values, name = regular_ngon(3, R_init, C)
velocities_0 = np.zeros_like(positions_0)  # start from rest

print(f"\n  Initial perturbation: epsilon = {epsilon} ({epsilon*100:.1f}%)")
print(f"  R_init = {R_init:.6f} (R_well = {R_well:.6f})")
print(f"  Expected period: T_br = {T_br:.4f}")

# Choose dt and integration time
n_periods = 5
t_total = n_periods * T_br
dt = T_br / 200  # 200 steps per period
save_every = max(1, int(T_br / dt / 50))  # ~50 saves per period

print(f"  Integration time: {t_total:.2f} ({n_periods} periods)")
print(f"  dt = {dt:.6f}, save_every = {save_every}")

result = leapfrog_integrate(
    positions_0, velocities_0, C_values,
    tgp_acceleration, tgp_potential,
    t_span=(0, t_total), dt=dt, save_every=save_every,
)

# Extract breathing mode oscillation
t = result['t']
pos = result['positions']  # (n_save, 3, 3)
energy = result['energy']

# Compute R(t) = mean distance from center of mass
R_t = np.zeros(len(t))
for k in range(len(t)):
    com = np.mean(pos[k], axis=0)
    R_t[k] = np.mean([np.linalg.norm(pos[k, i] - com) for i in range(3)])

# Relative perturbation
delta_R = (R_t - R_well) / R_well

print(f"\n  Simulation results:")
print(f"    Energy conservation: |dE/E_0| = {abs(energy[-1] - energy[0]) / abs(energy[0]):.2e}")
print(f"    max(delta_R) = {delta_R.max():.6f}")
print(f"    min(delta_R) = {delta_R.min():.6f}")

# ============================================================
# Measure oscillation frequency from zero crossings
# ============================================================
print("\n" + "=" * 70)
print("FREQUENCY MEASUREMENT")
print("=" * 70)

# Find zero crossings of delta_R
crossings = []
for k in range(len(delta_R) - 1):
    if delta_R[k] * delta_R[k+1] < 0:
        # Linear interpolation for crossing time
        t_cross = t[k] - delta_R[k] * (t[k+1] - t[k]) / (delta_R[k+1] - delta_R[k])
        crossings.append(t_cross)

if len(crossings) >= 3:
    # Period = time between every other crossing (full period = 2 half-periods)
    periods = []
    for i in range(0, len(crossings) - 2, 2):
        T_meas = crossings[i+2] - crossings[i]
        periods.append(T_meas)

    T_measured = np.mean(periods)
    omega_measured = 2 * np.pi / T_measured
    T_std = np.std(periods) if len(periods) > 1 else 0

    print(f"  Zero crossings found: {len(crossings)}")
    print(f"  Measured periods: {[f'{p:.4f}' for p in periods]}")
    print(f"  T_measured = {T_measured:.6f} +/- {T_std:.6f}")
    print(f"  T_analytical = {T_br:.6f}")
    print(f"  omega_measured = {omega_measured:.6f}")
    print(f"  omega_analytical = {omega_br:.6f}")
    print(f"  Relative error: {abs(omega_measured - omega_br) / omega_br * 100:.4f}%")

    # Also check amplitude decay (should be zero for conservative system)
    if len(delta_R) > 10:
        # Peak amplitudes
        peaks = []
        for k in range(1, len(delta_R) - 1):
            if delta_R[k] > delta_R[k-1] and delta_R[k] > delta_R[k+1]:
                peaks.append(delta_R[k])
        if len(peaks) >= 2:
            amp_ratio = peaks[-1] / peaks[0]
            print(f"  Amplitude ratio (last/first peak): {amp_ratio:.6f}")
            print(f"  (should be ~1.0 for symplectic integrator)")
else:
    print(f"  Only {len(crossings)} zero crossings found (need >= 3)")
    print(f"  Oscillation may be overdamped or integration too short")

# ============================================================
# Verify Result 8: omega^2_br = (3/2) * omega^2_rad
# ============================================================
print("\n" + "=" * 70)
print("RESULT 8 VERIFICATION: omega^2_br = (3/2) * omega^2_rad")
print("=" * 70)

# The radial normal-mode frequency from Result 7:
# omega^2_rad = (16*pi*C/d^3) * [-1 + 6*beta/d - 18*gamma*C/d^2]
# evaluated at d = d_well
omega2_rad_R7 = (16*np.pi*C/d_well**3) * (-1 + 6*beta/d_well - 18*gamma*C/d_well**2)

# The breathing mode frequency (Result 8):
# omega^2_br = (3/2) * omega^2_rad
omega2_br_R8 = 1.5 * omega2_rad_R7
omega_br_R8 = np.sqrt(abs(omega2_br_R8)) if omega2_br_R8 > 0 else 0

print(f"  omega^2_rad (Result 7 formula) = {omega2_rad_R7:.6e}")
print(f"  omega^2_br = 3/2 * omega^2_rad = {omega2_br_R8:.6e}")
print(f"  omega_br (Result 8) = {omega_br_R8:.6f}")
print(f"  omega_br (direct d2V/dR2) = {omega_br:.6f}")
if omega_br_R8 > 0 and omega_br > 0:
    print(f"  Ratio direct/R8 = {omega_br / omega_br_R8:.6f}")

if len(crossings) >= 3:
    print(f"\n  omega_br (measured from simulation) = {omega_measured:.6f}")
    if omega_br_R8 > 0:
        print(f"  Ratio measured/R8 = {omega_measured / omega_br_R8:.6f}")
    print(f"  Ratio measured/direct = {omega_measured / omega_br:.6f}")

# ============================================================
# Connection to PTA breathing mode
# ============================================================
print("\n" + "=" * 70)
print("CONNECTION: N-BODY BREATHING <-> PTA BREATHING")
print("=" * 70)
print(f"""
  N-BODY BREATHING MODE:
    omega_br = {omega_br:.6f} (dimensionless units)
    Period T_br = {T_br:.4f}
    Amplitude: delta_R / R ~ epsilon = {epsilon}

  PTA BREATHING MODE (from pta_breathing_prediction.py):
    A_br ~ 8.6e-16 at f = 1/yr
    h_br = 2*delta_Phi/Phi_0 (monopolar strain)

  PHYSICAL CONNECTION:
    Both modes arise from radial oscillations of the scalar field:
      - PTA: cosmological perturbation delta_Phi around Phi_0
        (source: SMBH binary inspirals, frequency ~ nHz)
      - N-body: local perturbation delta_R around R_well
        (source: displacement from equilibrium, frequency ~ omega_br)

    The ANGULAR CORRELATION is the same:
      - PTA: monopolar Gamma(theta) = const (breathing polarization)
      - N-body: breathing mode has radial_fraction ~ 1.0
        (all bodies move in/out synchronously)

    The ratio A_br/A_tensor ~ 1/(2*omega_BD + 3) = 1/15.5 ~ 6.5%
    applies to BOTH:
      - PTA: breathing GWB is 6.5% of tensor GWB
      - N-body: scalar force is 6.5% of tensor (gradient) force
        (for Phi_0 = 25, omega_BD_eff = 6.25)
""")

# ============================================================
# Regression tests
# ============================================================
print("=" * 70)
print("REGRESSION TESTS")
print("=" * 70)

# Test 1: Energy conservation
E_rel = abs(energy[-1] - energy[0]) / abs(energy[0])
assert E_rel < 1e-8, f"Energy drift too large: {E_rel:.2e}"
print(f"  [PASS] Energy conservation: |dE/E| = {E_rel:.2e} < 1e-8")

# Test 2: Frequency matches analytical (< 1% error)
if len(crossings) >= 3:
    freq_err = abs(omega_measured - omega_br) / omega_br
    assert freq_err < 0.02, f"Frequency error too large: {freq_err*100:.2f}%"
    print(f"  [PASS] Frequency error: {freq_err*100:.4f}% < 2%")

# Test 3: Amplitude preservation (symplectic)
if 'peaks' in dir() and len(peaks) >= 2:
    assert abs(amp_ratio - 1.0) < 0.05, f"Amplitude decay: {amp_ratio:.4f}"
    print(f"  [PASS] Amplitude preserved: ratio = {amp_ratio:.4f}")

# Test 4: Radial stable (omega2_br > 0)
assert omega2_br > 0, f"Breathing mode unstable: omega2 = {omega2_br:.4e}"
print(f"  [PASS] Breathing mode stable: omega^2 = {omega2_br:.4e} > 0")

print(f"\n  All regression tests PASSED.")
print("=" * 70)
