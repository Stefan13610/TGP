"""
ex31_rotation_curve.py
======================
Galactic rotation curves: TGP (Path B) vs CDM (NFW) vs MOND vs Baryons only.

PHYSICAL SETUP:
  TGP Path B: 2-body Yukawa potential V_2B = -(C1*C2/r)*exp(-m_sp*r)
  For ordinary baryonic matter: C = m/(2*sqrt(pi)*m_Pl)
    => C^2 = m1*m2/(4*pi*m_Pl^2) = G*m1*m2  (Newtonian gravity recovered!)
  TGP modifies Newtonian gravity: V_TGP = -G*m1*m2/r * exp(-m_sp*r)
  Rotation curve: v_circ^2(r) = r * |dPhi/dr|

MODEL GALAXY (Milky Way-like):
  - Stellar disk: exponential, M_disk = 5e10 M_sun, R_d = 3.5 kpc
  - Bulge: Hernquist profile, M_bulge = 1e10 M_sun, r_b = 0.7 kpc
  - Gas disk: exponential, M_gas = 0.5e10 M_sun, R_gas = 7.0 kpc
  - NFW dark matter halo: M_200 = 1e12 M_sun, c_NFW = 15
  - MOND: Milgrom mu(x) = x/sqrt(1+x^2), a0 = 1.2e-10 m/s^2

MODELS COMPARED:
  1. Baryons only (disk + bulge + gas), Newtonian
  2. CDM: Baryons + NFW dark matter halo
  3. TGP: Yukawa-screened gravity, m_sp = 1/kpc, 1/(5kpc), 1/(15kpc)
  4. MOND: Milgrom interpolation

KEY FINDING:
  TGP Path B with Yukawa screening REDUCES gravity at r > 1/m_sp.
  => TGP CANNOT produce flat rotation curves (it falls below Newtonian).
  => MOND and CDM can match observed flat rotation curves.
  => TGP Path B (ordinary matter) reproduces Newton at r << 1/m_sp
     and is BELOW Newtonian at r >> 1/m_sp => cannot explain dark matter.
"""

import sys
import os
import numpy as np
from scipy.special import i0, i1, k0, k1

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
G      = 6.674e-11    # m^3 kg^-1 s^-2
M_sun  = 1.989e30     # kg
kpc    = 3.086e19     # m
km_s   = 1e3          # m/s
a0     = 1.2e-10      # m/s^2  (MOND Milgrom constant)

# ---------------------------------------------------------------------------
# Galaxy parameters (Milky Way-like)
# ---------------------------------------------------------------------------
M_disk  = 5.0e10 * M_sun     # stellar disk mass [kg]
R_d     = 3.5 * kpc          # stellar disk scale radius [m]

M_bulge = 1.0e10 * M_sun     # bulge mass [kg]
r_b     = 0.7 * kpc          # Hernquist bulge scale radius [m]

M_gas   = 0.5e10 * M_sun     # gas disk mass [kg]
R_gas   = 7.0 * kpc          # gas disk scale radius [m]

# NFW halo parameters
M_200   = 1.0e12 * M_sun     # virial mass [kg]
c_NFW   = 15.0               # NFW concentration

# Derive NFW parameters
# r_200: radius where mean density = 200 * rho_crit
rho_crit = 9.47e-27          # kg/m^3 (Planck 2018)
r_200 = (3.0 * M_200 / (4.0 * np.pi * 200.0 * rho_crit))**(1.0/3.0)
r_s = r_200 / c_NFW          # NFW scale radius [m]
# rho_s from M_200 = 4*pi*rho_s*r_s^3*[ln(1+c) - c/(1+c)]
rho_s = M_200 / (4.0 * np.pi * r_s**3 * (np.log(1.0 + c_NFW) - c_NFW/(1.0 + c_NFW)))

# TGP Yukawa screening masses [1/m]
m_sp_vals = {
    "1/kpc"    : 1.0 / (1.0  * kpc),
    "1/(5kpc)" : 1.0 / (5.0  * kpc),
    "1/(15kpc)": 1.0 / (15.0 * kpc),
}

# ---------------------------------------------------------------------------
# Radii grid
# ---------------------------------------------------------------------------
r_kpc_table = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0])
r_m_table   = r_kpc_table * kpc

# Fine grid for smooth curves (used internally for MOND, diagnostics)
r_kpc_fine = np.linspace(0.1, 35.0, 500)
r_m_fine   = r_kpc_fine * kpc

# ===========================================================================
# Component rotation curves
# ===========================================================================

def v_exp_disk(r_m, M_tot, R_scale):
    """
    Circular velocity for a razor-thin exponential disk (Freeman 1970).
      Sigma(R) = M_tot/(2*pi*R_scale^2) * exp(-R/R_scale)
      v^2(R) = 4*pi*G*Sigma0*R_scale * y^2 * [I0(y)*K0(y) - I1(y)*K1(y)]
      y = R / (2*R_scale)
    Uses modified Bessel functions I0, I1, K0, K1.
    """
    Sigma0 = M_tot / (2.0 * np.pi * R_scale**2)
    y = r_m / (2.0 * R_scale)
    y = np.where(np.asarray(y) < 1e-10, 1e-10, y)
    v2 = (4.0 * np.pi * G * Sigma0 * R_scale
          * y**2 * (i0(y) * k0(y) - i1(y) * k1(y)))
    return np.sqrt(np.maximum(v2, 0.0))


def v_hernquist_bulge(r_m, M_b, r_b_scale):
    """
    Circular velocity for a Hernquist bulge.
      rho(r) = M_b / (2*pi) * r_b / (r*(r+r_b)^3)
      M(<r) = M_b * r^2 / (r + r_b)^2
      v^2(r) = G * M(<r) / r
    """
    M_enc = M_b * r_m**2 / (r_m + r_b_scale)**2
    v2    = G * M_enc / r_m
    return np.sqrt(np.maximum(v2, 0.0))


def v_nfw(r_m):
    """
    Circular velocity for the NFW dark matter halo.
      v^2(r) = 4*pi*G*rho_s*r_s^3/r * [ln(1+r/r_s) - (r/r_s)/(1+r/r_s)]
    """
    x  = r_m / r_s
    v2 = 4.0 * np.pi * G * rho_s * r_s**3 / r_m * (np.log(1.0 + x) - x / (1.0 + x))
    return np.sqrt(np.maximum(v2, 0.0))


def v_baryons(r_m):
    """Total baryonic circular velocity (disk + bulge + gas), Newtonian."""
    v2_disk  = v_exp_disk(r_m, M_disk,  R_d)**2
    v2_bulge = v_hernquist_bulge(r_m, M_bulge, r_b)**2
    v2_gas   = v_exp_disk(r_m, M_gas,  R_gas)**2
    return np.sqrt(v2_disk + v2_bulge + v2_gas)


def v_cdm(r_m):
    """CDM: baryons + NFW halo."""
    return np.sqrt(v_baryons(r_m)**2 + v_nfw(r_m)**2)


def f_yukawa(r_m, m_sp):
    """
    Yukawa screening factor for a point mass (used as approximation for disk).
      f(r) = exp(-m_sp*r) * (1 + m_sp*r)
    This is the ratio of Yukawa radial acceleration to Newtonian:
      a_Yukawa / a_Newton = exp(-m_sp*r)*(1 + m_sp*r)
    = 1 at r -> 0 (Newton recovered)
    -> 0 at r -> inf (gravity screened)
    """
    t = m_sp * r_m
    return np.exp(-t) * (1.0 + t)


def v_tgp(r_m, m_sp):
    """
    TGP rotation curve: Yukawa-screened gravity applied to each baryonic component.
      v_TGP^2(R) = v_Newton^2(R) * f_yukawa(R, m_sp)
    The screening factor f_yukawa is the point-mass approximation:
      valid for the spherical bulge exactly, and as an approximation for disks.
    Each component is screened independently then combined in quadrature.
    """
    f = f_yukawa(r_m, m_sp)
    v2_disk  = v_exp_disk(r_m, M_disk,  R_d)**2 * f
    v2_bulge = v_hernquist_bulge(r_m, M_bulge, r_b)**2 * f
    v2_gas   = v_exp_disk(r_m, M_gas,  R_gas)**2 * f
    return np.sqrt(np.maximum(v2_disk + v2_bulge + v2_gas, 0.0))


def v_mond(r_m):
    """
    MOND circular velocity using Milgrom interpolation mu(x) = x/sqrt(1+x^2).
      mu(a/a0) * a = a_Newton
      => a^2 / sqrt(a^2 + a0^2) = a_Newton  (rearranged)
    Solve analytically:
      Let u = (a_Newton/a0)^2
      => x^2 = a^2/a0^2 satisfies x^2 / sqrt(1 + x^2) = u^{1/2}
      => x^4 = u*(1 + x^2)
      => x^4 - u*x^2 - u = 0   [quadratic in x^2]
      => x^2 = (u + sqrt(u^2 + 4u)) / 2
      => a_MOND = x * a0
      => v^2 = a_MOND * r
    """
    a_N = v_baryons(r_m)**2 / r_m    # Newtonian baryonic acceleration [m/s^2]
    u   = (a_N / a0)**2
    x2  = 0.5 * (u + np.sqrt(u**2 + 4.0 * u))
    a_M = np.sqrt(x2) * a0
    v2  = a_M * r_m
    return np.sqrt(np.maximum(v2, 0.0))


# ===========================================================================
# Compute all models at table radii
# ===========================================================================

print("=" * 72)
print("ex31_rotation_curve.py")
print("Galactic Rotation Curves: TGP (Path B) vs CDM vs MOND")
print("=" * 72)
print()
print("Galaxy parameters (Milky Way-like):")
print(f"  Stellar disk : M_disk  = {M_disk/M_sun:.2e} M_sun, R_d    = {R_d/kpc:.1f} kpc")
print(f"  Bulge        : M_bulge = {M_bulge/M_sun:.2e} M_sun, r_b    = {r_b/kpc:.1f} kpc")
print(f"  Gas disk     : M_gas   = {M_gas/M_sun:.2e} M_sun, R_gas  = {R_gas/kpc:.1f} kpc")
print(f"  NFW halo     : M_200   = {M_200/M_sun:.2e} M_sun, c_NFW  = {c_NFW:.0f}")
print(f"                 r_200   = {r_200/kpc:.1f} kpc,  r_s    = {r_s/kpc:.1f} kpc")
print(f"                 rho_s   = {rho_s*kpc**3/M_sun:.3e} M_sun/kpc^3")
print(f"  MOND a0      : {a0:.2e} m/s^2")
print()

# Evaluate all models
v_bar   = v_baryons(r_m_table) / km_s          # km/s
v_CDM   = v_cdm(r_m_table)     / km_s
v_MOND  = np.array([v_mond(r)   / km_s for r in r_m_table])

v_TGP = {}
for label, m_sp in m_sp_vals.items():
    v_TGP[label] = np.array([v_tgp(r, m_sp) / km_s for r in r_m_table])

# ===========================================================================
# Output table
# ===========================================================================
print("-" * 72)
header = (f"{'r [kpc]':>8s}  {'Bar-only':>9s}  "
          + "  ".join(f"{'TGP '+lbl:>11s}" for lbl in m_sp_vals)
          + f"  {'CDM+NFW':>9s}  {'MOND':>8s}")
print(header)
print("-" * 72)

for i, r in enumerate(r_kpc_table):
    tgp_cols = "  ".join(f"{v_TGP[lbl][i]:>11.2f}" for lbl in m_sp_vals)
    print(f"{r:>8.1f}  {v_bar[i]:>9.2f}  {tgp_cols}  {v_CDM[i]:>9.2f}  {v_MOND[i]:>8.2f}")

print("-" * 72)
print("All velocities in km/s.")
print()

# ===========================================================================
# Yukawa screening factors at key radii
# ===========================================================================
print("=" * 72)
print("Yukawa screening factor f(r) = exp(-m_sp*r)*(1+m_sp*r)")
print("(ratio of TGP to Newtonian acceleration; 1=Newton, 0=fully screened)")
print("-" * 72)
hdr2 = f"{'r [kpc]':>8s}  " + "  ".join(f"{'m_sp='+lbl:>14s}" for lbl in m_sp_vals)
print(hdr2)
print("-" * 72)
for i, r in enumerate(r_kpc_table):
    cols = "  ".join(
        f"{f_yukawa(r_m_table[i], m_sp):>14.4f}"
        for m_sp in m_sp_vals.values()
    )
    print(f"{r:>8.1f}  {cols}")
print("-" * 72)
print()

# ===========================================================================
# MOND interpolation check
# ===========================================================================
print("=" * 72)
print("MOND interpolation function mu(x) = x/sqrt(1+x^2) at each radius:")
print("  (x = a_Newton/a0; mu~x for x<<1 deep-MOND, mu~1 for x>>1 Newton)")
print("-" * 72)
print(f"{'r [kpc]':>8s}  {'a_N [m/s^2]':>14s}  {'x=a_N/a0':>12s}  {'mu(x)':>8s}  {'regime':>12s}")
print("-" * 72)
for i, r in enumerate(r_kpc_table):
    a_N = v_baryons(r_m_table[i])**2 / r_m_table[i]
    x   = a_N / a0
    mu  = x / np.sqrt(1.0 + x**2)
    regime = "deep-MOND" if x < 0.3 else ("Newton" if x > 3 else "transition")
    print(f"{r:>8.1f}  {a_N:>14.3e}  {x:>12.4f}  {mu:>8.4f}  {regime:>12s}")
print("-" * 72)
print()

# ===========================================================================
# Velocity ratios: TGP / Baryons-only
# ===========================================================================
print("=" * 72)
print("Velocity ratio v_TGP / v_Baryons (shows Yukawa suppression):")
print("  Ratio < 1: TGP predicts LESS gravity than Newton")
print("  Ratio = 1: No modification (small r limit)")
print("-" * 72)
hdr3 = f"{'r [kpc]':>8s}  " + "  ".join(f"{'TGP '+lbl:>14s}" for lbl in m_sp_vals)
print(hdr3)
print("-" * 72)
for i, r in enumerate(r_kpc_table):
    v_b = v_bar[i]
    ratios = "  ".join(
        f"{v_TGP[lbl][i]/v_b:>14.4f}" if v_b > 0 else f"{'---':>14s}"
        for lbl in m_sp_vals
    )
    print(f"{r:>8.1f}  {ratios}")
print("-" * 72)
print()

# ===========================================================================
# Physical analysis
# ===========================================================================
print("=" * 72)
print("PHYSICAL ANALYSIS: TGP Path B and Galactic Rotation Curves")
print("=" * 72)
print()
print("TGP Path B: The 2-body Yukawa potential for baryonic matter is")
print("  V_TGP(r) = -(C1*C2/r) * exp(-m_sp*r)")
print("  with C = m/(2*sqrt(pi)*m_Pl)  =>  C^2 = G*m1*m2  (Newton!)")
print()
print("So TGP RECOVERS Newtonian gravity with an exponential cutoff:")
print("  V_TGP = -G*m1*m2/r * exp(-m_sp*r)")
print()
print("Radial acceleration from Yukawa potential:")
print("  a_TGP = G*M/r^2 * exp(-m_sp*r) * (1 + m_sp*r)")
print("  => Screening factor f(r) = exp(-m_sp*r)*(1+m_sp*r)")
print()
print("KEY RESULT (from the table above):")
for lbl, m_sp in m_sp_vals.items():
    lam = 1.0 / m_sp / kpc
    f5  = f_yukawa(5.0  * kpc, m_sp)
    f20 = f_yukawa(20.0 * kpc, m_sp)
    f30 = f_yukawa(30.0 * kpc, m_sp)
    print(f"  m_sp = {lbl:>10s}  (lambda = {lam:5.1f} kpc):"
          f"  f(5kpc)={f5:.3f},  f(20kpc)={f20:.3f},  f(30kpc)={f30:.3f}")
print()
print("CONCLUSION:")
print()
print("  1. TGP with m_sp = 1/kpc:")
print("     Gravity is screened within ~1 kpc. At r=5 kpc, gravity is already")
print("     strongly suppressed. The rotation curve falls STEEPLY below Newtonian.")
print("     => Cannot produce a flat rotation curve.")
print()
print("  2. TGP with m_sp = 1/(5 kpc):")
print("     Screening occurs at ~5 kpc scale. The rotation curve peaks near")
print("     r~5 kpc and declines at larger radii.")
print("     => Cannot sustain a flat curve out to 30 kpc.")
print()
print("  3. TGP with m_sp = 1/(15 kpc):")
print("     Milder screening. The curve is close to Newtonian at small r,")
print("     but still declines significantly beyond 15 kpc.")
print("     => Still cannot produce a flat rotation curve at 20-30 kpc.")
print()
print("  4. CDM (NFW halo):")
print("     Baryons + NFW halo produces a nearly flat curve at 10-30 kpc,")
print("     consistent with observed Milky Way rotation (~220 km/s).")
print("     => CDM CAN explain flat rotation curves.")
print()
print("  5. MOND:")
print("     Deep-MOND limit (x << 1) gives v^4 = G*M_bar*a0 (Tully-Fisher!)")
print("     MOND boosts the baryonic acceleration in the low-acceleration regime.")
print("     => MOND CAN produce approximately flat rotation curves.")
print()
print("  FUNDAMENTAL ISSUE WITH TGP PATH B:")
print("  ------------------------------------")
print("  The Yukawa potential exp(-m_sp*r)/r is ALWAYS smaller than 1/r.")
print("  Therefore, TGP gravity is ALWAYS less than or equal to Newtonian.")
print("  A flat rotation curve REQUIRES extra gravity (or equivalent effect)")
print("  beyond the Newtonian baryonic prediction at large radii.")
print("  TGP Path B provides LESS gravity at large r, moving in the WRONG")
print("  direction to explain flat rotation curves.")
print()
print("  => TGP Path B CANNOT explain galactic rotation curves without")
print("     additional mass (dark matter) or modification (e.g., m_sp^2 < 0).")
print()

# ===========================================================================
# Comparison of asymptotic behaviors
# ===========================================================================
print("=" * 72)
print("ASYMPTOTIC BEHAVIOR COMPARISON (at r = 30 kpc):")
print("=" * 72)
r_ref = 30.0 * kpc
v_b30   = v_baryons(r_ref) / km_s
v_cdm30 = v_cdm(r_ref) / km_s
v_mon30 = v_mond(r_ref) / km_s
print(f"  Baryons only (Newton) : {v_b30:.1f} km/s")
for lbl, m_sp in m_sp_vals.items():
    v_t30 = v_tgp(r_ref, m_sp) / km_s
    diff  = v_t30 - v_b30
    print(f"  TGP m_sp={lbl:>10s}   : {v_t30:.1f} km/s  (delta = {diff:+.1f} km/s vs baryons)")
print(f"  CDM + NFW halo        : {v_cdm30:.1f} km/s  (delta = {v_cdm30-v_b30:+.1f} km/s vs baryons)")
print(f"  MOND                  : {v_mon30:.1f} km/s  (delta = {v_mon30-v_b30:+.1f} km/s vs baryons)")
print()
print("  TGP corrections are NEGATIVE (velocity reduced below Newtonian).")
print("  CDM/MOND corrections are POSITIVE (velocity boosted above Newtonian).")
print("  Observed flat rotation curves need POSITIVE corrections at large r.")
print()

# ===========================================================================
# Tully-Fisher relation check (MOND prediction)
# ===========================================================================
print("=" * 72)
print("TULLY-FISHER RELATION: MOND prediction v^4 = G*M_bar*a0")
print("=" * 72)
M_bar_total = M_disk + M_bulge + M_gas
v_TF_pred = (G * M_bar_total * a0)**0.25 / km_s
print(f"  Total baryonic mass : M_bar = {M_bar_total/M_sun:.2e} M_sun")
print(f"  MOND TF prediction  : v_TF  = (G*M_bar*a0)^(1/4) = {v_TF_pred:.1f} km/s")
print(f"  MOND v_circ(30kpc)  : {v_mon30:.1f} km/s  (should match v_TF in deep-MOND limit)")
print(f"  CDM  v_circ(30kpc)  : {v_cdm30:.1f} km/s")
print()
print("  Note: At 30 kpc, baryonic acceleration is not yet fully in deep-MOND")
print("  regime for MW-like parameters, so MOND v < v_TF (asymptotic).")
print()

# ===========================================================================
# Summary table for quick reference
# ===========================================================================
print("=" * 72)
print("SUMMARY: v_circ [km/s] at selected radii")
print("=" * 72)
r_sel = [5.0, 10.0, 20.0, 30.0]
r_sel_m = np.array(r_sel) * kpc

print(f"{'Model':>22s}", end="")
for r in r_sel:
    print(f"  {r:.0f} kpc", end="")
print()
print("-" * 72)

# Baryons only
print(f"{'Baryons only (Newton)':>22s}", end="")
for rm in r_sel_m:
    print(f"  {v_baryons(rm)/km_s:6.1f}", end="")
print()

# TGP models
for lbl, m_sp in m_sp_vals.items():
    print(f"{'TGP m_sp='+lbl:>22s}", end="")
    for rm in r_sel_m:
        print(f"  {v_tgp(rm, m_sp)/km_s:6.1f}", end="")
    print()

# CDM
print(f"{'CDM + NFW':>22s}", end="")
for rm in r_sel_m:
    print(f"  {v_cdm(rm)/km_s:6.1f}", end="")
print()

# MOND
print(f"{'MOND':>22s}", end="")
for rm in r_sel_m:
    print(f"  {v_mond(rm)/km_s:6.1f}", end="")
print()

print("=" * 72)
print()
print("FINAL CONCLUSION:")
print("  TGP Path B (baryonic matter with C=m/(2*sqrt(pi)*m_Pl)) REPRODUCES")
print("  Newtonian gravity at r << 1/m_sp, but with Yukawa screening it gives")
print("  LESS gravitational acceleration at r >> 1/m_sp than pure Newton.")
print()
print("  Observed galactic rotation curves are FLAT (constant v) at large r,")
print("  requiring MORE gravity than baryons alone provide in Newtonian theory.")
print()
print("  Therefore: TGP Path B CANNOT explain galactic rotation curves.")
print("             It cannot substitute for dark matter or MOND.")
print()
print("  This is a fundamental prediction: any Yukawa potential with m_sp > 0")
print("  is a SCREENED (weakened) version of Newtonian gravity.")
print("  Only m_sp -> 0 recovers Newton; m_sp > 0 always gives less force.")
