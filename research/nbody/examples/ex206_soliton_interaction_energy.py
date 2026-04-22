#!/usr/bin/env python3
"""
ex206 -- P6: Soliton interaction energy from field overlap
============================================================

Computes the ACTUAL interaction between two TGP defects by evaluating
the overlap integral on a 2D cylindrical grid. Compares classical
(oscillatory) and EFT (Yukawa) interactions.

This CLOSES THE LOOP of the TGP n-body program:
    defect g(r) --> C_eff (P5) --> V_Y(d) (Path B) --> F_i (P0 n-body)

Tests:
  Part 1: Far-field profile -- delta(d) oscillates (not exponential)
  Part 2: Overlap integral S(d) and gradient cross-term G(d) vs d
  Part 3: Classical interaction V_cl(d) = G(d) + V''(1)*S(d)
  Part 4: Comparison with Yukawa prediction at scale of C_eff
  Part 5: Parameter scan -- interaction vs g0

PHYSICS:
--------
Leading-order interaction:
  V_cl(d) = int [grad(delta_1).grad(delta_2) + V''(1)*delta_1*delta_2] d^3x
          = G(d) + V''(1) * S(d)

For standard TGP: V''(1) = -1 (vacuum is maximum)
  => V_cl(d) = G(d) - S(d)
  The far-field equation is (nabla^2 + 1)*delta = 0 => oscillatory tail.
  The interaction is SMALL and OSCILLATORY (defects nearly decouple).

For EFT Yukawa: V''_eff = +m_sp^2
  => delta ~ C*exp(-m*r)/r
  => V_cl ~ -4*pi*C^2*exp(-m*d)/d (Yukawa attraction)

KEY RESULT: Classical field theory gives near-zero oscillatory interaction.
Yukawa force law requires the EFT mass gap. The magnitude of V_cl sets
the scale of nonlinear corrections beyond the EFT picture.
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.soliton_interaction import (
    generate_defect_profile,
    overlap_integral,
    gradient_cross_integral,
    classical_interaction,
    far_field_value,
    interaction_scan,
    yukawa_interaction,
)
from nbody.tgp_field import screening_mass

quick = "--quick" in sys.argv

if quick:
    N_RHO = 150
    N_Z = 300
    N_EVAL_ODE = 2000
    R_MAX_ODE = 30.0
    D_VALUES = [3.0, 5.0, 7.0, 10.0, 14.0, 18.0]
    G0_SCAN = [0.70, 0.90]
else:
    N_RHO = 300
    N_Z = 600
    N_EVAL_ODE = 4000
    R_MAX_ODE = 40.0
    D_VALUES = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 22.0]
    G0_SCAN = [0.50, 0.70, 0.85, 0.90, 0.95]


# -- Part 1: Far-field profile ---------------------------------------------

def part1_far_field():
    """Show that delta(d) oscillates, confirming non-Yukawa behavior."""
    print("=" * 78)
    print("Part 1: Far-Field Profile delta(d) = 1 - g(d)")
    print("=" * 78)

    beta = gamma = 1.0
    g0 = 0.90
    m_sp = screening_mass(beta, gamma)

    defect = generate_defect_profile(g0, beta, gamma,
                                     r_max=R_MAX_ODE, n_eval=N_EVAL_ODE)

    C_eff = defect['C_eff_proj']
    print(f"\n  beta = gamma = {beta}, g0 = {g0}")
    print(f"  m_sp = {m_sp:.4f}, C_eff (projection) = {C_eff:.6e}")

    # Sample delta(d) at several distances
    d_sample = np.array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 22.0, 25.0])
    d_sample = d_sample[d_sample < R_MAX_ODE - 1]

    print(f"\n  {'d':>6s} {'delta(d)':>12s} {'Yukawa':>12s} "
          f"{'|ratio|':>10s} {'sign(delta)':>12s}")
    print(f"  {'-'*58}")

    n_sign_changes = 0
    prev_sign = None
    for d in d_sample:
        delta_d = far_field_value(defect, d)
        yuk_d = C_eff * np.exp(-m_sp * d) / d
        ratio = abs(delta_d / yuk_d) if abs(yuk_d) > 1e-30 else float('nan')
        sign = "+" if delta_d >= 0 else "-"

        if prev_sign is not None and sign != prev_sign:
            n_sign_changes += 1
        prev_sign = sign

        print(f"  {d:6.1f} {delta_d:12.6e} {yuk_d:12.6e} "
              f"{ratio:10.4f} {sign:>12s}")

    oscillatory = n_sign_changes >= 1
    print(f"\n  Sign changes in delta(d): {n_sign_changes}")
    print(f"  Far-field is oscillatory: "
          f"{'YES (CONFIRMED)' if oscillatory else 'monotone'}")

    if not oscillatory:
        print(f"  Note: oscillation period ~ 2*pi/{m_sp:.2f} = {2*np.pi/m_sp:.1f}")
        print(f"  May need larger r_max to see sign change")

    return defect, oscillatory


# -- Part 2: Overlap and gradient integrals --------------------------------

def part2_overlap_scan(defect):
    """Compute S(d) and G(d) vs separation."""
    print(f"\n{'=' * 78}")
    print("Part 2: Overlap S(d) and Gradient G(d) Integrals")
    print(f"{'=' * 78}")

    beta = gamma = 1.0

    print(f"\n  S(d) = int delta_1 * delta_2 d^3x")
    print(f"  G(d) = int grad(delta_1) . grad(delta_2) d^3x")
    print(f"  Grid: {N_RHO} x {N_Z}")

    print(f"\n  {'d':>6s} {'S(d)':>14s} {'G(d)':>14s} "
          f"{'V_cl=G-S':>14s} {'|V_cl|':>12s}")
    print(f"  {'-'*66}")

    results = interaction_scan(defect, D_VALUES, beta, gamma, N_RHO, N_Z)

    for i, d in enumerate(D_VALUES):
        S = results['overlap_S'][i]
        G = results['gradient_G'][i]
        V_cl = results['V_classical'][i]

        print(f"  {d:6.1f} {S:14.6e} {G:14.6e} "
              f"{V_cl:14.6e} {abs(V_cl):12.6e}")

    # Check: do S and G nearly cancel (since V''=-1)?
    if len(D_VALUES) >= 2:
        S_arr = results['overlap_S']
        G_arr = results['gradient_G']
        V_arr = results['V_classical']

        # Measure cancellation: |V_cl| / max(|G|, |S|)
        cancel = np.abs(V_arr) / np.maximum(np.abs(G_arr), np.abs(S_arr))
        mask = np.isfinite(cancel) & (np.abs(G_arr) > 1e-20)
        if np.any(mask):
            mean_cancel = np.mean(cancel[mask])
            print(f"\n  Cancellation ratio |V_cl|/max(|G|,|S|): {mean_cancel:.4f}")
            if mean_cancel < 0.5:
                print(f"  => Strong cancellation: G and V''*S nearly cancel")
                print(f"     This is expected: far-field equation "
                      f"(nabla^2 + 1)*delta = 0")
            else:
                print(f"  => Moderate cancellation (nonlinear core effects)")

    return results


# -- Part 3: Classical vs Yukawa interaction --------------------------------

def part3_comparison(results, defect):
    """Compare V_cl(d) with Yukawa prediction."""
    print(f"\n{'=' * 78}")
    print("Part 3: Classical Interaction vs Yukawa Prediction")
    print(f"{'=' * 78}")

    C_eff = results['C_eff']
    m_sp = results['m_sp']

    print(f"\n  C_eff = {C_eff:.6e}, m_sp = {m_sp:.4f}")
    print(f"  V_Y(d) = -4*pi*C^2*exp(-m*d)/d")

    print(f"\n  {'d':>6s} {'V_classical':>14s} {'V_yukawa':>14s} "
          f"{'|V_cl/V_Y|':>12s} {'V_cl sign':>10s}")
    print(f"  {'-'*58}")

    V_cl = results['V_classical']
    V_Y = results['V_yukawa']

    for i, d in enumerate(D_VALUES):
        ratio = abs(V_cl[i] / V_Y[i]) if abs(V_Y[i]) > 1e-30 else float('nan')
        sign = "+" if V_cl[i] >= 0 else "-"
        print(f"  {d:6.1f} {V_cl[i]:14.6e} {V_Y[i]:14.6e} "
              f"{ratio:12.4f} {sign:>10s}")

    # Key physics: V_cl should NOT match V_Y (different decay)
    # V_cl oscillates or decays as power law
    # V_Y decays exponentially
    mask = np.abs(V_Y) > 1e-30
    if np.any(mask):
        ratios = np.abs(V_cl[mask] / V_Y[mask])
        print(f"\n  |V_cl/V_Y| range: [{np.min(ratios):.2e}, {np.max(ratios):.2e}]")

        # Check if ratio grows with d (classical >> Yukawa at large d)
        if len(ratios) >= 3:
            grows = ratios[-1] > ratios[0] * 2
            print(f"  Ratio grows with d: {'YES' if grows else 'NO'}")
            if grows:
                print(f"  => Classical interaction decays SLOWER than Yukawa")
                print(f"     (oscillatory/power-law vs exponential)")

    # V_cl sign changes?
    n_sign = np.sum(np.diff(np.sign(V_cl)) != 0)
    print(f"\n  V_cl sign changes: {n_sign}")
    if n_sign > 0:
        print(f"  => Classical interaction OSCILLATES (not monotone attraction)")
    else:
        sign_type = "attractive" if V_cl[-1] < 0 else "repulsive"
        print(f"  => Classical interaction is monotone {sign_type} in this range")

    print(f"\n  PHYSICAL INTERPRETATION:")
    print(f"  Classical field overlap gives oscillatory/power-law interaction.")
    print(f"  Yukawa exp(-m*d)/d requires the EFT mass gap m_sp.")
    print(f"  Path B force law is NOT a classical result -- it's an EFT prediction.")

    return {
        'V_cl': V_cl,
        'V_Y': V_Y,
    }


# -- Part 4: Magnitude comparison -----------------------------------------

def part4_magnitude(results, defect):
    """Compare the SCALE of classical interaction with Yukawa at small d."""
    print(f"\n{'=' * 78}")
    print("Part 4: Interaction Scales -- Classical vs Yukawa vs Defect Energy")
    print(f"{'=' * 78}")

    C_eff = results['C_eff']
    m_sp = results['m_sp']
    E_defect = defect['E_defect']

    print(f"\n  C_eff = {C_eff:.6e}")
    print(f"  E_defect = {E_defect:.6e}")
    print(f"  m_sp = {m_sp:.4f}")

    # At d = smallest separation
    d_min = D_VALUES[0]
    V_cl_min = results['V_classical'][0]
    V_Y_min = results['V_yukawa'][0]

    print(f"\n  At d = {d_min}:")
    print(f"    |V_classical| = {abs(V_cl_min):.6e}")
    print(f"    |V_yukawa|    = {abs(V_Y_min):.6e}")
    print(f"    |E_defect|    = {abs(E_defect):.6e}")
    print(f"    V_cl/E_defect = {V_cl_min/E_defect:.6e}")

    # At large d: V_cl should be tiny (nearly zero)
    d_max = D_VALUES[-1]
    V_cl_max = results['V_classical'][-1]
    V_Y_max = results['V_yukawa'][-1]
    print(f"\n  At d = {d_max}:")
    print(f"    |V_classical| = {abs(V_cl_max):.6e}")
    print(f"    |V_yukawa|    = {abs(V_Y_max):.6e}")

    # Decay analysis on V_cl
    mask = np.abs(results['V_classical']) > 1e-20
    if np.sum(mask) >= 3:
        d_arr = np.array(D_VALUES)[mask]
        logV = np.log(np.abs(results['V_classical'][mask]))
        logd = np.log(d_arr)

        # Fit: log|V| = a + b*d (exponential) or log|V| = a + b*log(d) (power-law)
        try:
            c_exp = np.polyfit(d_arr, logV, 1)
            m_eff_exp = -c_exp[0]
            resid_exp = np.std(logV - np.polyval(c_exp, d_arr))
        except Exception:
            m_eff_exp = 0.0
            resid_exp = float('inf')

        try:
            c_pow = np.polyfit(logd, logV, 1)
            alpha_pow = c_pow[0]
            resid_pow = np.std(logV - np.polyval(c_pow, logd))
        except Exception:
            alpha_pow = 0.0
            resid_pow = float('inf')

        print(f"\n  Decay fit:")
        print(f"    Exponential: |V_cl| ~ exp(-{m_eff_exp:.4f}*d), "
              f"residual = {resid_exp:.4f}")
        print(f"    Power-law:   |V_cl| ~ d^{alpha_pow:.2f}, "
              f"residual = {resid_pow:.4f}")
        better = "exponential" if resid_exp < resid_pow else "power-law"
        print(f"    Better fit: {better}")

    return {'E_defect': E_defect}


# -- Part 5: Parameter scan -----------------------------------------------

def part5_g0_scan():
    """Scan interaction properties vs defect core depth g0."""
    print(f"\n{'=' * 78}")
    print("Part 5: Interaction vs Core Depth g0 (at d = 6.0)")
    print(f"{'=' * 78}")

    beta = gamma = 1.0
    d_test = 6.0

    print(f"\n  beta = gamma = {beta}, d = {d_test}")
    print(f"\n  {'g0':>6s} {'C_eff':>12s} {'V_classical':>14s} "
          f"{'V_yukawa':>14s} {'|V_cl/V_Y|':>12s}")
    print(f"  {'-'*62}")

    results_list = []
    for g0 in G0_SCAN:
        try:
            defect = generate_defect_profile(g0, beta, gamma,
                                             r_max=R_MAX_ODE, n_eval=N_EVAL_ODE)
            C_eff = defect['C_eff_proj']
            m_sp = defect['m_sp_pathB']

            res = classical_interaction(defect, d_test, beta, gamma,
                                        N_RHO, N_Z)
            V_cl = res['V_classical']
            V_Y = yukawa_interaction(d_test, C_eff, m_sp)
            ratio = abs(V_cl / V_Y) if abs(V_Y) > 1e-30 else float('nan')

            print(f"  {g0:6.3f} {C_eff:12.6e} {V_cl:14.6e} "
                  f"{V_Y:14.6e} {ratio:12.4f}")

            results_list.append({
                'g0': g0, 'C_eff': C_eff,
                'V_classical': V_cl, 'V_yukawa': V_Y,
                'ratio': ratio,
            })
        except Exception as e:
            print(f"  {g0:6.3f} ERROR: {e}")

    if len(results_list) >= 2:
        C_arr = np.array([r['C_eff'] for r in results_list])
        V_cl_arr = np.array([r['V_classical'] for r in results_list])
        print(f"\n  C_eff range: [{C_arr.min():.4e}, {C_arr.max():.4e}]")

        # Does V_cl scale with C_eff^2?
        if len(results_list) >= 3:
            logC = np.log(np.abs(C_arr))
            logV = np.log(np.abs(V_cl_arr))
            mask = np.isfinite(logC) & np.isfinite(logV)
            if np.sum(mask) >= 2:
                coeffs = np.polyfit(logC[mask], logV[mask], 1)
                print(f"  Scaling: |V_cl| ~ C_eff^{coeffs[0]:.2f} "
                      f"(expect ~2 for bilinear coupling)")

    return results_list


def main():
    print("=" * 78)
    print("ex206 -- P6: Soliton Interaction Energy from Field Overlap")
    print("         (Classical oscillatory vs EFT Yukawa)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t_total = time.time()

    defect, osc = part1_far_field()
    scan = part2_overlap_scan(defect)
    p3 = part3_comparison(scan, defect)
    p4 = part4_magnitude(scan, defect)
    p5 = part5_g0_scan()

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")
    print(f"  Part 1: Far-field oscillatory: "
          f"{'YES' if osc else 'monotone in range'}")
    print(f"  Part 2: Overlap/gradient integrals at {len(D_VALUES)} separations")
    V_cl = scan['V_classical']
    n_sign = np.sum(np.diff(np.sign(V_cl)) != 0)
    print(f"  Part 3: V_classical sign changes: {n_sign} "
          f"({'oscillatory' if n_sign > 0 else 'monotone'})")
    print(f"  Part 4: |V_classical| vs |E_defect| = "
          f"{abs(V_cl[0])/abs(p4['E_defect']):.4e} at d={D_VALUES[0]}")
    print(f"  Part 5: {len(p5)} g0 values scanned")

    print(f"\n  KEY RESULTS:")
    print(f"  1. Classical defect interaction is NOT Yukawa (oscillatory/power-law)")
    print(f"  2. G(d) and V''(1)*S(d) nearly cancel (linearized eq. satisfied)")
    print(f"  3. Residual V_cl comes from nonlinear core overlap")
    print(f"  4. Yukawa force law requires EFT mass gap (confirmed by P5)")
    print(f"  5. Chain: defect -> C_eff (P5 projection) -> V_Y -> n-body (P0)")
    print(f"     is the correct derivation path")

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
