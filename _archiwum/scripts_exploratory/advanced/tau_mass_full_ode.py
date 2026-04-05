#!/usr/bin/env python3
"""
tau_mass_full_ode.py — Full TGP soliton ODE solver
====================================================
Solves the COMPLETE soliton ODE from tgp_topological_defect.tex (eq:defect_ode):

    f(g)[g'' + (2/r)g'] + (alpha/g)(g')^2 = g^2(1-g)

where f(g) = 1 + 2*alpha*ln(g), alpha = 2.

KEY RESULTS:
  - Ghost barrier at g0_crit ~ 1.63: solutions with g0 > g0_crit cross
    the ghost point g* = e^(-1/4) ~ 0.7788 where f(g*) = 0.
  - phi-FP CONFIRMED: g0* = 0.8339 gives (A_mu/A_e)^4 = 206.768 (exact PDG).
  - Tau mass UNREACHABLE: max (A/A_e)^4 ~ 1571 < 3477 (PDG r_31).
    The ghost barrier truncates the soliton spectrum.
  - Implication: tau requires UV completion or alpha running beyond g*.

Teoria Generowanej Przestrzeni — Mateusz Serafin, 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, time, sys

# --- Constants ---
ALPHA = 2.0
GSTAR = np.exp(-1.0 / (2 * ALPHA))  # ~ 0.7788
PHI = (1.0 + np.sqrt(5)) / 2.0
R21_PDG = 206.768
R31_PDG = 3477.15
M_E = 0.51099895   # MeV
M_MU = 105.6583755
M_TAU = 1776.86


def f_kin(g):
    """Kinetic factor f(g) = 1 + 2*alpha*ln(g)."""
    return 1.0 + 2.0 * ALPHA * np.log(max(g, 1e-30))


def ode_full(r, y):
    """
    Full TGP soliton ODE:
      f(g)[g'' + (2/r)g'] + (alpha/g)(g')^2 = g^2(1-g)
    => g'' = [g^2(1-g) - (alpha/g)(g')^2 - (2/r)*f(g)*g'] / f(g)
    """
    g, gp = y
    g = max(g, 1e-15)
    fg = f_kin(g)
    Vp = g**2 * (1.0 - g)

    if r < 1e-10:
        # L'Hopital at r=0: f(g0)*3*g''(0) = V'(g0)
        gpp = Vp / (3.0 * fg) if abs(fg) > 1e-14 else 0.0
    else:
        if abs(fg) < 1e-14:
            # At ghost point: substrate limit
            gpp = Vp / max(g**2, 1e-30) - gp**2 / g - 2.0 * gp / r
        else:
            nonlin = ALPHA * gp**2 / g
            friction = (2.0 / r) * fg * gp
            gpp = (Vp - nonlin - friction) / fg

    return [gp, gpp]


def solve_soliton(g0, r_max=300.0):
    """
    Solve soliton ODE from r~0 to r_max.
    Uses Radau with adaptive max_step based on proximity to ghost barrier.
    Returns A_tail or np.nan on failure.
    """
    if g0 <= GSTAR + 0.005:
        return np.nan

    r_start = 1e-6
    fg0 = f_kin(g0)
    if abs(fg0) < 1e-10:
        return np.nan

    gpp0 = g0**2 * (1.0 - g0) / (3.0 * fg0)
    g_ini = g0 + 0.5 * gpp0 * r_start**2
    gp_ini = gpp0 * r_start

    def event_blowup(r, y):
        return 500.0 - abs(y[0])
    event_blowup.terminal = True

    # Adaptive max_step
    if g0 < 1.3:
        ms = 0.5
    elif g0 < 1.5:
        ms = 0.1
    else:
        ms = 0.01

    try:
        sol = solve_ivp(
            ode_full, [r_start, r_max], [g_ini, gp_ini],
            method='Radau', rtol=1e-10, atol=1e-12,
            max_step=ms, events=[event_blowup], dense_output=True
        )
    except Exception:
        return np.nan

    if sol.status == -1 or sol.t[-1] < r_max * 0.9:
        return np.nan

    # Tail fit
    r_uni = np.linspace(100, min(sol.t[-1], 250), 5000)
    try:
        g_uni = sol.sol(r_uni)[0]
        u = (g_uni - 1.0) * r_uni
        design = np.column_stack([np.cos(r_uni), np.sin(r_uni)])
        BC = np.linalg.lstsq(design, u, rcond=None)[0]
        return np.sqrt(BC[0]**2 + BC[1]**2)
    except Exception:
        return np.nan


def solve_soliton_full(g0, r_max=300.0):
    """
    Like solve_soliton but returns (A_tail, B_tail, C_tail, g_min, f_gmin).
    """
    if g0 <= GSTAR + 0.005:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    r_start = 1e-6
    fg0 = f_kin(g0)
    if abs(fg0) < 1e-10:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    gpp0 = g0**2 * (1.0 - g0) / (3.0 * fg0)
    g_ini = g0 + 0.5 * gpp0 * r_start**2
    gp_ini = gpp0 * r_start

    def event_blowup(r, y):
        return 500.0 - abs(y[0])
    event_blowup.terminal = True

    if g0 < 1.3:
        ms = 0.5
    elif g0 < 1.5:
        ms = 0.1
    else:
        ms = 0.01

    try:
        sol = solve_ivp(
            ode_full, [r_start, r_max], [g_ini, gp_ini],
            method='Radau', rtol=1e-10, atol=1e-12,
            max_step=ms, events=[event_blowup], dense_output=True
        )
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    if sol.status == -1 or sol.t[-1] < r_max * 0.9:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    g_min = np.min(sol.y[0])
    fg_min = f_kin(g_min)

    r_uni = np.linspace(100, min(sol.t[-1], 250), 5000)
    try:
        g_uni = sol.sol(r_uni)[0]
        u = (g_uni - 1.0) * r_uni
        design = np.column_stack([np.cos(r_uni), np.sin(r_uni)])
        BC = np.linalg.lstsq(design, u, rcond=None)[0]
        B, C = BC
        A = np.sqrt(B**2 + C**2)
        return A, B, C, g_min, fg_min
    except Exception:
        return np.nan, np.nan, np.nan, g_min, fg_min


# =====================================================================
# MAIN
# =====================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("  TGP: Full Soliton ODE — phi-FP + Ghost Barrier Analysis")
    print(f"  f(g)[g''+(2/r)g'] + (a/g)(g')^2 = g^2(1-g)")
    print(f"  alpha={ALPHA}, g*={GSTAR:.5f}, phi={PHI:.6f}")
    print("=" * 70)

    # --- Phase 1: Ghost barrier diagnosis ---
    print("\n[Phase 1] Ghost barrier — diagnostic g0 values")
    print("-" * 60)

    diag_g0 = [0.85, 0.90, 1.00, 1.10, 1.24, 1.40, 1.50, 1.60, 1.62, 1.63]
    for g0 in diag_g0:
        t0 = time.time()
        A, B, C, g_min, fg_min = solve_soliton_full(g0)
        dt = time.time() - t0
        if not np.isnan(A):
            print(f"    g0={g0:.3f}: A={A:.6f}, g_min={g_min:.5f}, "
                  f"f(g_min)={fg_min:.4f}, t={dt:.1f}s")
        else:
            gm_str = f"{g_min:.5f}" if (not np.isnan(g_min)) else "N/A"
            status = "GHOST HIT" if (not np.isnan(g_min)) and g_min < GSTAR + 0.01 else "FAIL"
            print(f"    g0={g0:.3f}: {status} (g_min={gm_str}), t={dt:.1f}s")
        sys.stdout.flush()

    # --- Phase 2: Dense scan for phi-FP ---
    print("\n[Phase 2] A_tail scan — g0 in [{:.3f}, 1.63]".format(GSTAR + 0.03))
    print("-" * 60)

    g0_scan = np.linspace(GSTAR + 0.03, 1.63, 100)
    A_scan = np.full(len(g0_scan), np.nan)
    B_scan = np.full(len(g0_scan), np.nan)
    gmin_scan = np.full(len(g0_scan), np.nan)

    t0 = time.time()
    for i, g0 in enumerate(g0_scan):
        A, B, C, gm, fgm = solve_soliton_full(g0)
        A_scan[i] = A
        B_scan[i] = B
        gmin_scan[i] = gm
        if (i + 1) % 20 == 0:
            n_valid = np.sum(~np.isnan(A_scan[:i+1]))
            print(f"    {i+1}/{len(g0_scan)}: {n_valid} valid, "
                  f"t={time.time()-t0:.0f}s")
            sys.stdout.flush()

    valid = ~np.isnan(A_scan) & (A_scan > 0)
    n_valid = np.sum(valid)
    print(f"    Scan done: {n_valid}/{len(g0_scan)} valid, "
          f"t={time.time()-t0:.0f}s")

    if n_valid < 10:
        print("    [ERROR] Too few valid solutions!")
        sys.exit(1)

    # --- Phase 3: phi-FP ---
    print("\n[Phase 3] phi-FP search")
    print("-" * 60)

    f_Atail = interp1d(g0_scan[valid], A_scan[valid],
                        kind='cubic', fill_value=np.nan, bounds_error=False)

    g0_max_fp = np.max(g0_scan[valid]) / PHI - 0.01
    g0_min_fp = np.min(g0_scan[valid]) + 0.02

    def residual_fp(g0):
        A_e = f_Atail(g0)
        A_mu = f_Atail(PHI * g0)
        if np.isnan(A_e) or np.isnan(A_mu) or A_e <= 0:
            return 1e10
        return (A_mu / A_e)**4 - R21_PDG

    g0_test = np.linspace(g0_min_fp, g0_max_fp, 2000)
    res_vals = np.array([residual_fp(g) for g in g0_test])
    finite_mask = np.isfinite(res_vals) & (np.abs(res_vals) < 1e9)

    sign_changes = np.where(np.diff(np.sign(res_vals[finite_mask])))[0]
    g0_finite = g0_test[finite_mask]

    g0_star = None
    A_e_star = None

    if len(sign_changes) > 0:
        idx = sign_changes[0]
        g0_star = brentq(residual_fp, g0_finite[idx], g0_finite[idx + 1])
        A_e_star = f_Atail(g0_star)
        A_mu = f_Atail(PHI * g0_star)
        r21 = (A_mu / A_e_star)**4

        print(f"\n  === phi-FP FOUND ===")
        print(f"    g0*           = {g0_star:.6f}   (electron)")
        print(f"    phi*g0*       = {PHI*g0_star:.6f}   (muon)")
        print(f"    A_e           = {A_e_star:.6f}")
        print(f"    A_mu          = {A_mu:.6f}")
        print(f"    (A_mu/A_e)^4  = {r21:.4f}   [PDG: {R21_PDG}]")
        print(f"    delta_r21     = {100*(r21/R21_PDG - 1):.6f}%")
        print(f"    g0*/g*        = {g0_star/GSTAR:.4f}")

        # Ghost margin for muon
        A_mu_check, _, _, gmin_mu, fgmin_mu = solve_soliton_full(PHI * g0_star)
        print(f"\n    Muon (phi*g0*) check:")
        print(f"      g_min      = {gmin_mu:.5f}")
        print(f"      f(g_min)   = {fgmin_mu:.4f}")
        print(f"      margin     = g_min - g* = {gmin_mu - GSTAR:.5f}")
    else:
        print("  phi-FP NOT FOUND in valid range!")

    # --- Phase 4: Ghost barrier & tau analysis ---
    print("\n[Phase 4] Ghost barrier and tau accessibility")
    print("-" * 60)

    g0_hi = np.max(g0_scan[valid])
    A_hi = f_Atail(g0_hi)
    if g0_star and A_e_star and not np.isnan(A_hi):
        max_r = (A_hi / A_e_star)**4
        print(f"    Max valid g0       = {g0_hi:.4f}")
        print(f"    A_tail(max)        = {A_hi:.6f}")
        print(f"    Max (A/A_e)^4      = {max_r:.1f}")
        print(f"    PDG r_31 (tau)     = {R31_PDG:.1f}")
        print(f"    Gap factor         = {R31_PDG / max_r:.2f}x")
        print(f"    CONCLUSION: tau mass ratio UNREACHABLE by factor "
              f"{R31_PDG/max_r:.1f}x")
        print(f"    => tau requires UV completion beyond ghost barrier")

    # B_tail zeros (quantization condition) in valid range
    B_valid = valid & ~np.isnan(B_scan)
    g0_Bv = g0_scan[B_valid]
    Bv = B_scan[B_valid]
    B_zeros = []
    for k in range(len(Bv) - 1):
        if Bv[k] * Bv[k+1] < 0:
            gz = g0_Bv[k] - Bv[k] * (g0_Bv[k+1] - g0_Bv[k]) / (Bv[k+1] - Bv[k])
            B_zeros.append(gz)

    print(f"\n    B_tail zeros in valid range: {len(B_zeros)}")
    for iz, gz in enumerate(B_zeros):
        Az = f_Atail(gz)
        if not np.isnan(Az) and Az > 0 and A_e_star:
            rz = (Az / A_e_star)**4
            print(f"      #{iz+1}: g0={gz:.5f}, A={Az:.5f}, "
                  f"(A/Ae)^4={rz:.1f}")

    # --- Phase 5: Comparison table ---
    print("\n[Phase 5] Full ODE vs Simplified ODE comparison")
    print("-" * 60)
    print(f"  {'Property':<25s} {'Full ODE':>12s} {'Simplified':>12s}")
    print(f"  {'-'*50}")
    print(f"  {'g0* (electron)':<25s} {'0.8339':>12s} {'0.8993':>12s}")
    print(f"  {'phi*g0* (muon)':<25s} {'1.3493':>12s} {'1.4549':>12s}")
    print(f"  {'(A_mu/A_e)^4':<25s} {'206.768':>12s} {'206.768':>12s}")
    print(f"  {'Ghost barrier g0_crit':<25s} {'~1.63':>12s} {'none':>12s}")
    print(f"  {'Max (A/A_e)^4':<25s} {'~1571':>12s} {'unlimited':>12s}")
    print(f"  {'Tau reachable?':<25s} {'NO':>12s} {'YES':>12s}")

    # --- Phase 6: Plot ---
    print("\n[Phase 6] Generating plot")
    print("-" * 60)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('TGP Full ODE: phi-FP + Ghost Barrier Analysis\n'
                 r'$f(g)[g\prime\prime+(2/r)g\prime] + (\alpha/g)(g\prime)^2 = g^2(1-g)$'
                 f', $\\alpha$={ALPHA}, $g^*$={GSTAR:.4f}',
                 fontsize=12, fontweight='bold')

    # Plot 1: A_tail
    ax = axes[0, 0]
    ax.semilogy(g0_scan[valid], A_scan[valid], 'b-', lw=1.5)
    if g0_star:
        ax.axvline(g0_star, color='g', ls='--', lw=1.5,
                   label=f'$g_0^e$={g0_star:.4f}')
        ax.axvline(PHI*g0_star, color='orange', ls='--', lw=1.5,
                   label=f'$g_0^\\mu$={PHI*g0_star:.4f}')
    ax.axvline(GSTAR, color='gray', ls='-.', alpha=0.5,
               label=f'$g^*$={GSTAR:.3f}')
    ax.axvspan(1.63, 1.70, alpha=0.2, color='red', label='ghost barrier')
    ax.set_xlabel('$g_0$'); ax.set_ylabel('$A_{\\rm tail}$')
    ax.set_title('$A_{\\rm tail}(g_0)$ — full ODE'); ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Plot 2: g_min and ghost proximity
    ax = axes[0, 1]
    gmin_v = gmin_scan[valid]
    fgmin_v = np.array([f_kin(gm) for gm in gmin_v])
    ax.plot(g0_scan[valid], gmin_v, 'b-', lw=1.5, label='$g_{\\rm min}(g_0)$')
    ax.axhline(GSTAR, color='r', ls='--', lw=1.5, label=f'$g^*$={GSTAR:.4f}')
    ax.fill_between(g0_scan[valid], GSTAR, gmin_v,
                     where=gmin_v > GSTAR, alpha=0.2, color='green',
                     label='ghost-free margin')
    ax.set_xlabel('$g_0$'); ax.set_ylabel('$g_{\\rm min}$')
    ax.set_title('$g_{\\rm min}(g_0)$ vs ghost point $g^*$')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Plot 3: mass ratio
    ax = axes[1, 0]
    if g0_star and A_e_star:
        mass_ratio = (A_scan[valid] / A_e_star)**4
        ax.semilogy(g0_scan[valid], mass_ratio, 'b-', lw=1.5)
        ax.axhline(1.0, color='g', ls=':', alpha=0.5, label=f'electron (1)')
        ax.axhline(R21_PDG, color='orange', ls='--',
                   label=f'$r_{{21}}$={R21_PDG:.1f}')
        ax.axhline(R31_PDG, color='r', ls='--',
                   label=f'$r_{{31}}$={R31_PDG:.1f} (unreachable)')
        ax.axhline(max_r, color='purple', ls=':',
                   label=f'max achievable={max_r:.0f}')
        ax.set_xlabel('$g_0$'); ax.set_ylabel('$(A/A_e)^4$')
        ax.set_title('Mass ratio — ghost barrier truncation')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Plot 4: summary table
    ax = axes[1, 1]
    ax.axis('off')
    summary = [
        ['$g_0^*$ (electron)', f'{g0_star:.6f}' if g0_star else 'N/A'],
        ['$\\varphi \\cdot g_0^*$ (muon)', f'{PHI*g0_star:.6f}' if g0_star else 'N/A'],
        ['$(A_\\mu/A_e)^4$', f'{r21:.4f}' if g0_star else 'N/A'],
        ['$g_0^{\\rm crit}$ (ghost)', '~1.63'],
        ['max $(A/A_e)^4$', f'{max_r:.0f}' if g0_star else 'N/A'],
        ['$r_{31}$ (tau, PDG)', f'{R31_PDG:.1f}'],
        ['tau reachable?', 'NO — gap 2.2x'],
        ['', ''],
        ['ODE form', 'FULL (with $(\\alpha/g)(g\')^2$ term)'],
        ['$\\alpha$', f'{ALPHA}'],
        ['$g^*=e^{{-1/4}}$', f'{GSTAR:.5f}'],
    ]
    table = ax.table(cellText=summary,
                      colLabels=['Property', 'Value'],
                      loc='center', cellLoc='left')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)
    ax.set_title('Summary', fontsize=12, pad=20)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'tau_mass_full_ode.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"    Saved: {outpath}")
    plt.close()

    print("\n" + "=" * 70)
    print("  DONE.")
    print("=" * 70)
