#!/usr/bin/env python3
"""
stage1_reshoot_V0.py — Stage 1 verification of Z-test (Hubble tension solver candidate)

Stage 0 result (PRELIMINARY): TGP gives |dH/H| = 8.93%, ostensibly matching tension at 107%.
But Stage 0 had:
  - V_0 fixed at 0.685 (de2 baseline) — possibly inconsistent with source term ON
  - psi_init = 1 - 1e-3 — possibly wrong direction (push down vs up)
  - Used formula |dH/H| = 0.5 |dL/L| Omega_L — magnitude only, sign not tracked

Stage 1 tasks:
  1. Test MULTIPLE initial conditions: psi_init in {1-1e-3, 1, 1+1e-3, 1+1e-2}
  2. For each, re-shoot V_0 to enforce Omega_DE0(today) = 0.685
  3. Track Lambda(z) trajectory: is Lambda_recomb > Lambda_today (EDE-like, fixes tension)
     or Lambda_recomb < Lambda_today (LDE-like, makes tension WORSE)?
  4. Compute H(z) properly with shifted Lambda(z); compare with LCDM
  5. Final verdict: TGP solves Hubble tension (PASS) or makes worse / wrong direction (FAIL)

Pass criteria (from ROADMAP.md):
  - V_0 shoot succeeds: Omega_DE0(today) = 0.685 +/- 0.001
  - |dH/H| in [0.05, 0.12] (50-150% Hubble tension)
  - DIRECTION: Lambda_recomb > Lambda_today (EDE-like, makes H_0 inferred from CMB > 67.4)

If direction wrong (Lambda decay) -> TGP makes tension WORSE. Stage 0 conclusion REVERSED.
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp


# ============================================================================
# Cosmological parameters
# ============================================================================
OMEGA_M0 = 0.315
OMEGA_R0 = 9.1e-5
OMEGA_DE0_OBS = 1.0 - OMEGA_M0 - OMEGA_R0  # ~0.6849

PHI_0 = 24.6302
g_tilde = 5.0 * np.exp(2)/(12 * np.pi)  # delta.1+delta.2 closure
S_PRECISE = 18 * OMEGA_DE0_OBS / (PHI_0**2 * g_tilde)  # ~ 0.0207

DELTA_H_REQUIRED = 0.0837  # |H_SH0ES - H_Planck| / H_Planck

print("=" * 78)
print("STAGE 1: Re-shoot V_0 + initial conditions sensitivity test")
print("=" * 78)
print(f"\n  Cosmology: Omega_m = {OMEGA_M0}, Omega_DE0_obs = {OMEGA_DE0_OBS:.4f}")
print(f"  Phi_0 = {PHI_0}, g_tilde = {g_tilde:.4f}")
print(f"  s_precise = {S_PRECISE:.4e}")
print(f"  Required Hubble shift: |dH/H| = {DELTA_H_REQUIRED:.4f}")


# ============================================================================
# Potential
# ============================================================================
def V(psi):
    """V(psi) = 4 psi^3 - 3 psi^4 (de2 normalized)."""
    return 4.0*psi**3 - 3.0*psi**4

def dVdpsi(psi):
    """V'(psi) = 12 psi^2 (1-psi)."""
    return 12.0 * psi**2 * (1.0 - psi)


# ============================================================================
# FRW solver
# ============================================================================
def solve_FRW(s_coupling, V_0, psi_init, a_init=1e-5, n_eval=8000):
    """Integrate FRW Phi-EOM from a_init to a=1."""
    def eom_N(N, y):
        psi, u = y
        a = np.exp(N)
        rho_m = OMEGA_M0 / a**3
        rho_r = OMEGA_R0 / a**4
        rho_psi = 0.5 * u**2 + V_0 * V(psi)
        rho_total = rho_m + rho_r + rho_psi
        H = np.sqrt(max(rho_total, 1e-30))
        if H < 1e-30:
            return [0, 0]
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

    # Compute Lambda(z) and H(z)
    rho_m_vals = OMEGA_M0 / a_vals**3
    rho_r_vals = OMEGA_R0 / a_vals**4
    rho_psi_vals = 0.5 * u_vals**2 + V_0 * V(psi_vals)
    H_vals = np.sqrt(np.maximum(rho_m_vals + rho_r_vals + rho_psi_vals, 1e-30))

    return {
        'a': a_vals, 'z': z_vals, 'psi': psi_vals, 'u': u_vals,
        'rho_psi': rho_psi_vals, 'rho_m': rho_m_vals, 'H': H_vals
    }


def shoot_V0(s_coupling, psi_init, target=OMEGA_DE0_OBS, V0_low=0.05, V0_high=3.0,
             tol=1e-5, max_iter=60):
    """Bisection: find V_0 such that rho_psi(today) = target."""
    for it in range(max_iter):
        V_mid = 0.5 * (V0_low + V0_high)
        result = solve_FRW(s_coupling, V_mid, psi_init)
        if result is None:
            V0_high = V_mid
            continue
        rho_psi_today = result['rho_psi'][-1]
        diff = rho_psi_today - target
        if abs(diff) < tol:
            return V_mid, result, True, it
        if diff > 0:
            V0_high = V_mid
        else:
            V0_low = V_mid
    return V_mid, result, False, max_iter


# ============================================================================
# Test multiple initial conditions
# ============================================================================
print("\n" + "=" * 78)
print("Testing multiple initial conditions (with shooting V_0):")
print("=" * 78)

initial_conditions = [
    ('psi_init = 1 - 1e-3 (Stage 0)', 1.0 - 1e-3),
    ('psi_init = 1 (vacuum exact)',   1.0),
    ('psi_init = 1 + 1e-3',            1.0 + 1e-3),
    ('psi_init = 1 + 1e-2',            1.0 + 1e-2),
    ('psi_init = 1 + 1e-1',            1.0 + 1e-1),
]

results_summary = []

for label, psi_init in initial_conditions:
    print(f"\n  IC: {label}")
    V_0_shot, result, converged, iters = shoot_V0(S_PRECISE, psi_init)

    if not converged:
        print(f"    SHOOT FAILED ({iters} iters; final V_0={V_0_shot:.4f})")
        continue

    print(f"    Shot V_0 = {V_0_shot:.6f} (converged in {iters} iter)")

    # Extract values at z=0 and z=1100
    z = result['z']
    psi = result['psi']
    u = result['u']
    rho_psi = result['rho_psi']
    H = result['H']

    psi_today = psi[-1]
    rho_psi_today = rho_psi[-1]
    H_today = H[-1]

    idx_recomb = np.argmin(np.abs(z - 1100.0))
    psi_recomb = psi[idx_recomb]
    rho_psi_recomb = rho_psi[idx_recomb]
    H_recomb = H[idx_recomb]

    # LCDM comparison: same Omega_DE_const = OMEGA_DE0_OBS
    rho_psi_recomb_LCDM = OMEGA_DE0_OBS  # constant
    H_recomb_LCDM = np.sqrt(OMEGA_M0/result['a'][idx_recomb]**3
                            + OMEGA_R0/result['a'][idx_recomb]**4
                            + OMEGA_DE0_OBS)

    Lambda_recomb = rho_psi_recomb
    Lambda_today  = rho_psi_today

    print(f"    psi_today  = {psi_today:.6f}, psi_recomb = {psi_recomb:.6f}")
    print(f"    Lambda_today = {Lambda_today:.6f} (should = {OMEGA_DE0_OBS:.4f})")
    print(f"    Lambda_recomb = {Lambda_recomb:.6f}")
    print(f"    Lambda_recomb / Lambda_today = {Lambda_recomb/Lambda_today:.4f}")

    direction = "EDE-like (recomb > today: SOLVES tension)" if Lambda_recomb > Lambda_today else "LDE-like (recomb < today: WORSENS tension)"
    print(f"    Direction: {direction}")

    # Compare H at recomb
    H_ratio_recomb = H_recomb / H_recomb_LCDM
    print(f"    H_TGP(recomb) / H_LCDM(recomb) = {H_ratio_recomb:.6f}")

    # Hubble tension implication via sound horizon
    # r_s ~ integral c_s/H(z) dz from recomb to infinity
    # Approximation: r_s shift ~ -1/2 dlnH/dz integrated weight
    # For EDE: H_recomb higher -> r_s smaller -> H_0 inferred higher
    # Crude: dH_0_inferred/H_0 ~ +(H_TGP/H_LCDM - 1)|_recomb (positive if recomb H higher)
    # But also need to account for D_A effect (post-recomb integration)

    # Simpler estimate: integrate angular comoving distance
    # D_A(recomb) = a_recomb * integral_0^z_recomb c dz / H(z)
    a_arr = result['a']
    H_arr = result['H']
    H_LCDM_arr = np.sqrt(OMEGA_M0/a_arr**3 + OMEGA_R0/a_arr**4 + OMEGA_DE0_OBS)

    # Comoving distance integration: ∫ dz/H(z) from z=0 to z=1100
    # In terms of a (decreasing from 1 to a_recomb): dz/H = -da/(a^2 H)
    a_today_idx = -1
    a_recomb_idx = idx_recomb
    # Integrate from a_recomb to a_today (i.e., today value of integral)
    da = np.diff(a_arr[a_recomb_idx:])
    a_mid = 0.5 * (a_arr[a_recomb_idx:-1] + a_arr[a_recomb_idx+1:])
    H_mid = 0.5 * (H_arr[a_recomb_idx:-1] + H_arr[a_recomb_idx+1:])
    H_LCDM_mid = 0.5 * (H_LCDM_arr[a_recomb_idx:-1] + H_LCDM_arr[a_recomb_idx+1:])

    D_A_TGP = np.sum(da / (a_mid**2 * H_mid))
    D_A_LCDM = np.sum(da / (a_mid**2 * H_LCDM_mid))
    D_A_ratio = D_A_TGP / D_A_LCDM

    print(f"    D_A(z_recomb,today) ratio TGP/LCDM = {D_A_ratio:.6f}")

    # For fixed observed theta* = r_s/D_A:
    # If TGP r_s same (only Lambda_recomb shifts which is small at high z) and D_A_TGP > D_A_LCDM,
    # then theta* should shift, breaking CMB observation OR requiring different H_0.
    # Approximate H_0 inferred shift: H_0_inferred_TGP / H_0_inferred_LCDM ~ D_A_LCDM / D_A_TGP
    H0_inferred_ratio = 1.0 / D_A_ratio
    delta_H_implied = (H0_inferred_ratio - 1.0)

    print(f"    Implied H_0_inferred shift via D_A: {delta_H_implied:+.4f} ({delta_H_implied*100:+.2f}%)")

    results_summary.append({
        'label': label,
        'psi_init': psi_init,
        'V_0_shot': V_0_shot,
        'psi_today': psi_today,
        'psi_recomb': psi_recomb,
        'Lambda_today': Lambda_today,
        'Lambda_recomb': Lambda_recomb,
        'L_ratio': Lambda_recomb/Lambda_today,
        'H_ratio_recomb': H_ratio_recomb,
        'D_A_ratio': D_A_ratio,
        'delta_H': delta_H_implied,
        'direction': direction,
    })


# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 78)
print("SUMMARY OF ALL INITIAL CONDITIONS")
print("=" * 78)

print(f"\n  {'IC':<25} {'V_0_shot':>10} {'psi_today':>11} {'psi_recomb':>12} {'L_ratio':>9} {'dH/H':>9}")
print(f"  {'-'*25} {'-'*10} {'-'*11} {'-'*12} {'-'*9} {'-'*9}")
for r in results_summary:
    label_short = r['label'][:24]
    print(f"  {label_short:<25} {r['V_0_shot']:>10.4f} {r['psi_today']:>11.6f} {r['psi_recomb']:>12.6f} {r['L_ratio']:>9.4f} {r['delta_H']:>+9.4f}")

print(f"\n  Required Hubble tension: dH/H = +{DELTA_H_REQUIRED:.4f}")
print(f"  (POSITIVE delta_H means TGP increases inferred H_0, helping resolve tension)")

# Find the best scenario
best_match = min(results_summary, key=lambda r: abs(r['delta_H'] - DELTA_H_REQUIRED))
print(f"\n  Best match initial condition: {best_match['label']}")
print(f"  delta_H = {best_match['delta_H']:+.4f} vs required {DELTA_H_REQUIRED:+.4f}")
print(f"  Match: {best_match['delta_H']/DELTA_H_REQUIRED*100:.1f}% of required")


# ============================================================================
# Verdict
# ============================================================================
print("\n" + "=" * 78)
print("STAGE 1 VERDICT")
print("=" * 78)

# Count direction classifications
ede_like = sum(1 for r in results_summary if r['L_ratio'] > 1.0)
lde_like = sum(1 for r in results_summary if r['L_ratio'] < 1.0)
total = len(results_summary)

print(f"\n  Out of {total} initial conditions tested:")
print(f"    {ede_like} show EDE-like (Lambda_recomb > Lambda_today)")
print(f"    {lde_like} show LDE-like (Lambda_recomb < Lambda_today)")

# Check Stage 0 PASS criteria
sample_results = [r for r in results_summary if 0.05 <= abs(r['delta_H']) <= 0.12 and r['delta_H'] > 0]
print(f"\n  Initial conditions with dH/H in [0.05, 0.12] AND positive direction: {len(sample_results)} / {total}")

if len(sample_results) > 0:
    print("\n  >>> STAGE 1 PARTIAL PASS — at least one scenario solves H_0 tension correctly")
    for r in sample_results:
        print(f"      {r['label']}: dH/H = {r['delta_H']:+.4f}, V_0 = {r['V_0_shot']:.4f}")
elif ede_like == 0 and lde_like > 0:
    print("\n  >>> STAGE 1 FAIL — ALL scenarios show LDE-like behavior")
    print("      TGP precise s gives Lambda DECAY which makes Hubble tension WORSE not better")
    print("      Stage 0 conclusion REVERSED: TGP precise NOT a Hubble tension solver")
elif all(abs(r['delta_H']) < 0.01 for r in results_summary):
    print("\n  >>> STAGE 1 NULL — all scenarios give dH/H < 1%")
    print("      Stage 0 magnitude estimate was wrong; tension impact negligible")
else:
    print("\n  >>> STAGE 1 INCONCLUSIVE — needs deeper analysis")

print("\n" + "=" * 78)
print("Done.")
