#!/usr/bin/env python3
"""
stage1_prime_G_running.py — Stage 1 prime: PROPER TGP Friedmann z G(psi)

Stage 1 NULL verdict assumed standard FRW: 3H^2 = 8 pi G_0 rho_total.
But sek04 ax:G says: G(Phi) = G_0 (Phi_0/Phi).
So PROPER TGP Friedmann is: 3H^2 = (8 pi G_0/3) rho_total / psi.

This Stage 1 prime computes WITH G(psi) running included.

Also includes c(psi) effect on distance integrals where applicable.
For light ray null geodesic on g_eff: c_eff(psi) ~ c_0/sqrt(psi).
But sek04 self-consistency principle says LOCAL observer measures c_0.
For COSMOLOGICAL distance ladder, integration uses c(psi(z)) along path.

Test: psi_today value AND H_today/H_LCDM_today ratio.

Pass criteria:
- H_today_TGP / H_today_LCDM in [1.05, 1.12] (5-12% boost, target 8.4%)
- BBN G drift (G(z=1e9)/G_today) in [0.95, 1.05] (within typical constraint)
- Ratio direction right (today H higher than LCDM today)
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.integrate import solve_ivp


# ============================================================================
# Cosmological parameters (Planck 2018)
# ============================================================================
OMEGA_M0 = 0.315
OMEGA_R0 = 9.1e-5
OMEGA_DE0_OBS = 1.0 - OMEGA_M0 - OMEGA_R0  # ~0.6849

PHI_0 = 24.6302
g_tilde = 5.0 * np.exp(2)/(12 * np.pi)
S_PRECISE = 18 * OMEGA_DE0_OBS / (PHI_0**2 * g_tilde)  # ~0.0207

DELTA_H_REQUIRED = 0.0837


# ============================================================================
# TGP physics with sek04 axioms
# ============================================================================
def G_TGP(psi):
    """G(psi) = G_0/psi (sek04 ax:G)."""
    return 1.0 / psi  # in units of G_0

def c_TGP(psi):
    """c(psi) = c_0/sqrt(psi) (sek04 ax:c)."""
    return 1.0 / np.sqrt(psi)  # in units of c_0


# ============================================================================
# Modified Friedmann z G(psi)
# ============================================================================
def V(psi):
    return 4.0*psi**3 - 3.0*psi**4

def dVdpsi(psi):
    return 12.0 * psi**2 * (1.0 - psi)


def H_TGP_modified(psi, u, a, V_0):
    """
    TGP Friedmann z G(psi):
    3H^2 = (8 pi G_0/3) rho_total / psi
    In dimensionless H_0=1, 8 pi G_0/3 = 1, rho_total in units rho_crit_today:
    H^2 = rho_total / psi (relative to today's H_0 normalization)

    BUT: today's normalization H_0 itself is in TGP units,
    so we need: H(today) = 1 by definition. This means
    rho_total(today) / psi_today = 1.
    """
    rho_m = OMEGA_M0 / a**3
    rho_r = OMEGA_R0 / a**4
    rho_psi = 0.5 * u**2 + V_0 * V(psi)
    rho_total = rho_m + rho_r + rho_psi
    H_squared = rho_total / psi  # G(psi) = G_0/psi modification
    return np.sqrt(max(H_squared, 1e-30))


def solve_FRW_modified(s_coupling, V_0, psi_init, a_init=1e-5, n_eval=8000):
    """Solve FRW Phi-EOM z G(psi) modified Friedmann."""
    def eom_N(N, y):
        psi, u = y
        a = np.exp(N)
        H = H_TGP_modified(psi, u, a, V_0)
        if H < 1e-30:
            return [0, 0]
        rho_m = OMEGA_M0 / a**3
        # EOM: psi_ddot + 3H psi_dot + V_0 V'(psi) + s rho_m = 0
        dpsi_dN = u / H
        du_dN = (-3.0*H*u - V_0*dVdpsi(psi) - s_coupling*rho_m) / H
        return [dpsi_dN, du_dN]

    N_init = np.log(a_init)
    N_eval = np.linspace(N_init, 0.0, n_eval)
    sol = solve_ivp(eom_N, [N_init, 0.0], [psi_init, 0.0],
                    t_eval=N_eval, method='RK45',
                    rtol=1e-10, atol=1e-13, max_step=0.01)
    if not sol.success:
        return None

    a_vals = np.exp(N_eval)
    z_vals = 1.0/a_vals - 1.0
    psi_vals = sol.y[0]
    u_vals = sol.y[1]

    rho_m_vals = OMEGA_M0 / a_vals**3
    rho_r_vals = OMEGA_R0 / a_vals**4
    rho_psi_vals = 0.5 * u_vals**2 + V_0 * V(psi_vals)
    H_vals = np.array([H_TGP_modified(psi_vals[i], u_vals[i], a_vals[i], V_0)
                       for i in range(len(a_vals))])
    G_vals = 1.0 / psi_vals  # G(psi) = G_0/psi
    c_vals = 1.0 / np.sqrt(psi_vals)  # c(psi) = c_0/sqrt(psi)

    return {
        'a': a_vals, 'z': z_vals, 'psi': psi_vals, 'u': u_vals,
        'rho_psi': rho_psi_vals, 'rho_m': rho_m_vals,
        'H': H_vals, 'G': G_vals, 'c': c_vals
    }


def shoot_V0_modified(s_coupling, psi_init, target_H_today=1.0,
                      V0_low=0.05, V0_high=3.0, tol=1e-5, max_iter=60):
    """
    Bisection: find V_0 such that H(today) = target.
    target_H_today = 1.0 means H_TGP_today = H_0 (Planck normalization).
    """
    for it in range(max_iter):
        V_mid = 0.5 * (V0_low + V0_high)
        result = solve_FRW_modified(s_coupling, V_mid, psi_init)
        if result is None:
            V0_high = V_mid
            continue
        H_today = result['H'][-1]
        diff = H_today - target_H_today
        if abs(diff) < tol:
            return V_mid, result, True, it
        if diff > 0:
            V0_high = V_mid
        else:
            V0_low = V_mid
    return V_mid, result, False, max_iter


# ============================================================================
# Run scenarios
# ============================================================================
print("=" * 78)
print("STAGE 1 PRIME: Proper TGP Friedmann z G(psi) running included")
print("=" * 78)
print(f"\n  Cosmology: Omega_m = {OMEGA_M0}, Omega_DE0_obs = {OMEGA_DE0_OBS:.4f}")
print(f"  s_precise = {S_PRECISE:.4e}")
print(f"  Required Hubble shift: |dH/H| = {DELTA_H_REQUIRED:.4f}")

print("\n" + "=" * 78)
print("Test 1: H_TGP / H_LCDM ratio with G(psi) included")
print("=" * 78)
print("""
  Pre-shooting test: just compute H_TGP using Stage 1 psi values,
  with G(psi) modification.
""")

# Stage 1 found psi_today ≈ 0.78 with V_0 = 0.83 shot.
# Let's compute H_TGP/H_LCDM with G(psi) = 1/psi included.

# At today (a=1):
psi_today_S1 = 0.78
H_LCDM_today = 1.0  # by normalization
H_TGP_today_with_Grun = np.sqrt(1.0 / psi_today_S1)  # rho_total normalized to 1
boost_today = H_TGP_today_with_Grun / H_LCDM_today
print(f"  At today (psi = {psi_today_S1}): G(psi)/G_0 = {1/psi_today_S1:.4f}")
print(f"    H_TGP_today / H_LCDM_today = sqrt(1/psi) = {boost_today:.4f}")
print(f"    Boost = {(boost_today - 1)*100:.2f}% (required: {DELTA_H_REQUIRED*100:.2f}%)")
print(f"    Ratio TGP/required: {(boost_today - 1)/DELTA_H_REQUIRED:.2f}x")

# At recomb:
psi_recomb_S1 = 0.984
H_LCDM_recomb = np.sqrt(OMEGA_M0/(1/1101)**3 + OMEGA_R0/(1/1101)**4 + OMEGA_DE0_OBS)
rho_total_recomb_TGP = OMEGA_M0*1101**3 + OMEGA_R0*1101**4 + OMEGA_DE0_OBS  # if Lambda const
H_TGP_recomb_with_Grun = np.sqrt(rho_total_recomb_TGP / psi_recomb_S1)
boost_recomb = H_TGP_recomb_with_Grun / H_LCDM_recomb
print(f"\n  At recomb (psi = {psi_recomb_S1}): G(psi)/G_0 = {1/psi_recomb_S1:.4f}")
print(f"    H_TGP_recomb / H_LCDM_recomb = {boost_recomb:.4f}")
print(f"    Boost at recomb = {(boost_recomb - 1)*100:.2f}%")


# ============================================================================
# Test 2: Re-shoot V_0 with modified Friedmann
# ============================================================================
print("\n" + "=" * 78)
print("Test 2: Re-shoot V_0 with modified Friedmann (G(psi) included)")
print("=" * 78)

initial_conditions = [
    ('psi_init = 1 - 1e-3', 1.0 - 1e-3),
    ('psi_init = 1 + 1e-3', 1.0 + 1e-3),
    ('psi_init = 1 + 1e-2', 1.0 + 1e-2),
]

results_summary = []

for label, psi_init in initial_conditions:
    print(f"\n  IC: {label}")
    V_0_shot, result, converged, iters = shoot_V0_modified(S_PRECISE, psi_init)

    if not converged:
        print(f"    SHOOT FAILED ({iters} iters; final V_0={V_0_shot:.4f})")
        continue

    z = result['z']
    psi = result['psi']
    H = result['H']
    G_vals = result['G']

    psi_today = psi[-1]
    H_today = H[-1]
    G_today = G_vals[-1]

    idx_recomb = np.argmin(np.abs(z - 1100.0))
    psi_recomb = psi[idx_recomb]
    H_recomb = H[idx_recomb]
    G_recomb = G_vals[idx_recomb]

    # BBN check (z = 1e9, deep matter era)
    a_BBN = 1.0/(1 + 1e9)
    if a_BBN > result['a'][0]:
        idx_BBN = np.argmin(np.abs(result['a'] - a_BBN))
        psi_BBN = psi[idx_BBN]
        G_BBN = G_vals[idx_BBN]
    else:
        psi_BBN = psi[0]  # use earliest available
        G_BBN = G_vals[0]

    # Reference: pure LCDM (no TGP modifications)
    a_arr = result['a']
    H_LCDM_arr = np.sqrt(OMEGA_M0/a_arr**3 + OMEGA_R0/a_arr**4 + OMEGA_DE0_OBS)

    H_LCDM_today_ref = np.sqrt(1.0)  # = 1 by normalization
    H_LCDM_recomb_ref = H_LCDM_arr[idx_recomb]

    # Hubble tension implication: ratio H_today_TGP / H_today_LCDM
    # Note: H_TGP_today is what SH0ES would measure (local Solar System)
    # H_LCDM_today is what Planck would infer assuming LCDM
    # Required: H_today_TGP / H_today_LCDM = 1.084 (8.4%)
    boost = H_today / H_LCDM_today_ref
    delta_H = boost - 1.0

    # BBN G drift
    G_BBN_drift = G_BBN / G_today
    BBN_OK = abs(G_BBN_drift - 1.0) < 0.05  # 5% tolerance

    print(f"    Shot V_0 = {V_0_shot:.6f} (converged in {iters} iter)")
    print(f"    psi_today  = {psi_today:.6f}, psi_recomb = {psi_recomb:.6f}")
    print(f"    H_today_TGP = {H_today:.6f}, H_today_LCDM = {H_LCDM_today_ref:.6f}")
    print(f"    Boost (TGP/LCDM): {boost:.6f} = +{delta_H*100:.2f}%")
    print(f"    Required: +{DELTA_H_REQUIRED*100:.2f}%")
    print(f"    Match: {delta_H/DELTA_H_REQUIRED*100:.1f}% of required")
    print(f"    G_BBN/G_today = {G_BBN_drift:.4f} (BBN OK if in [0.95, 1.05])")
    print(f"    BBN check: {'PASS' if BBN_OK else 'FAIL ({:+.0f}% drift)'.format((G_BBN_drift-1)*100)}")

    results_summary.append({
        'label': label, 'psi_init': psi_init, 'V_0_shot': V_0_shot,
        'psi_today': psi_today, 'psi_recomb': psi_recomb,
        'H_today': H_today, 'boost': boost, 'delta_H': delta_H,
        'G_BBN_drift': G_BBN_drift, 'BBN_OK': BBN_OK,
    })


# ============================================================================
# Verdict
# ============================================================================
print("\n" + "=" * 78)
print("STAGE 1' VERDICT")
print("=" * 78)

print(f"\n  {'IC':<25} {'V_0':>8} {'psi_today':>11} {'boost':>8} {'%':>7} {'G_BBN':>9} {'BBN':>5}")
print(f"  {'-'*25} {'-'*8} {'-'*11} {'-'*8} {'-'*7} {'-'*9} {'-'*5}")
for r in results_summary:
    bbn_str = "PASS" if r['BBN_OK'] else "FAIL"
    print(f"  {r['label']:<25} {r['V_0_shot']:>8.4f} {r['psi_today']:>11.6f} {r['boost']:>8.4f} {r['delta_H']*100:>+6.2f} {r['G_BBN_drift']:>9.4f} {bbn_str:>5}")

print(f"\n  Required Hubble boost: +{DELTA_H_REQUIRED*100:.2f}%")

# Find scenarios that match Hubble tension AND pass BBN
sweet_spots = [r for r in results_summary
               if 0.05 <= r['delta_H'] <= 0.12 and r['BBN_OK']]
overshoots = [r for r in results_summary if r['delta_H'] > 0.12]
undershoots = [r for r in results_summary if r['delta_H'] < 0.05 and r['delta_H'] > 0]
bbn_fail = [r for r in results_summary if not r['BBN_OK']]

print(f"\n  Sweet spot (5-12% AND BBN OK): {len(sweet_spots)} scenarios")
print(f"  Overshoot (>12%): {len(overshoots)} scenarios")
print(f"  Undershoot (<5%): {len(undershoots)} scenarios")
print(f"  BBN fail: {len(bbn_fail)} scenarios")

if len(sweet_spots) > 0:
    print("\n  >>> STAGE 1' PARTIAL PASS — TGP solves Hubble tension AND BBN safe")
    for r in sweet_spots:
        print(f"      {r['label']}: dH/H = {r['delta_H']*100:+.2f}%, V_0 = {r['V_0_shot']:.4f}")
elif len(overshoots) > 0 and len(bbn_fail) == 0:
    print("\n  >>> STAGE 1' OVERSHOOT — TGP gives more boost than needed; needs fine-tuning")
elif len(bbn_fail) > 0:
    print("\n  >>> STAGE 1' BBN_PROBLEM — G(psi) running gives correct H_0 but BBN tension")
    print("      Self-consistency principle (sek04 A4) may save this — needs further analysis")
else:
    print("\n  >>> STAGE 1' inconclusive")

print("\n" + "=" * 78)
print("Done.")
