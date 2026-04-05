#!/usr/bin/env python3
"""
tau_uv_completion.py — Three UV completion paths for the tau mass
==================================================================
The ghost barrier at g0_crit~1.63 prevents the full ODE from reaching r_31=3477.
This script tests three UV completion strategies:

  Path A: Substrate ODE (ghost-free) — phi-FP + tau scan
  Path B: Running alpha(g) — does lowering alpha near UV shift the barrier?
  Path C: Blended ODE — full ODE for g>g*, substrate for g<g*

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
GSTAR = np.exp(-1.0 / (2 * ALPHA))
PHI = (1.0 + np.sqrt(5)) / 2.0
R21_PDG = 206.768
R31_PDG = 3477.15
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86


def f_kin(g, alpha=ALPHA):
    return 1.0 + 2.0 * alpha * np.log(max(g, 1e-30))


# =====================================================================
# PATH A: Substrate ODE (ghost-free)
# g^2*g'' + g*(g')^2 + (2/r)*g^2*g' = V'(g) = g^2(1-g)
# => g'' = (1-g) - (g')^2/g - (2/r)*g'
# =====================================================================
def ode_substrate(r, y):
    g, gp = y
    g = max(g, 1e-15)
    Vp_over_g2 = 1.0 - g  # V'(g)/g^2 = g^2(1-g)/g^2 = 1-g
    if r < 1e-10:
        gpp = (1.0 - g) / 3.0  # L'Hopital: 3*g'' = V'/g^2 = 1-g
    else:
        gpp = Vp_over_g2 - gp**2 / g - (2.0 / r) * gp
    return [gp, gpp]


# =====================================================================
# PATH B: Running alpha(g)
# Full ODE but alpha depends on g:
#   alpha(g) = alpha_IR + (alpha_UV - alpha_IR) * sigma(g)
# where sigma smoothly transitions near g*
# =====================================================================
def alpha_running(g, alpha_ir=2.0, alpha_uv=1.0, g_trans=0.85, width=0.05):
    """Smooth transition: alpha_ir for g>g_trans, alpha_uv for g<g_trans."""
    return alpha_uv + (alpha_ir - alpha_uv) / (1.0 + np.exp(-(g - g_trans) / width))


def ode_running_alpha(r, y, alpha_ir=2.0, alpha_uv=1.0):
    g, gp = y
    g = max(g, 1e-15)
    al = alpha_running(g, alpha_ir, alpha_uv)
    fg = 1.0 + 2.0 * al * np.log(g)
    Vp = g**2 * (1.0 - g)
    if r < 1e-10:
        gpp = Vp / (3.0 * fg) if abs(fg) > 1e-14 else 0.0
    else:
        if abs(fg) < 1e-14:
            gpp = Vp / max(g**2, 1e-30) - gp**2 / g - 2.0 * gp / r
        else:
            nonlin = al * gp**2 / g
            friction = (2.0 / r) * fg * gp
            gpp = (Vp - nonlin - friction) / fg
    return [gp, gpp]


# =====================================================================
# PATH C: Blended ODE — smooth transition full->substrate near g*
# =====================================================================
def ode_blended(r, y, eps_blend=0.03):
    g, gp = y
    g = max(g, 1e-15)
    fg = f_kin(g)
    Vp = g**2 * (1.0 - g)

    # Weight: w=1 for full ODE (|f|>>eps), w=0 for substrate (|f|->0)
    w = min(abs(fg) / eps_blend, 1.0)

    if r < 1e-10:
        # Full: g'' = V'/(3f), Substrate: g'' = (1-g)/3
        gpp_full = Vp / (3.0 * fg) if abs(fg) > 1e-14 else (1.0 - g) / 3.0
        gpp_sub = (1.0 - g) / 3.0
        gpp = w * gpp_full + (1.0 - w) * gpp_sub
    else:
        gpp_sub = (1.0 - g) - gp**2 / g - (2.0 / r) * gp
        if abs(fg) < 1e-14:
            gpp = gpp_sub
        else:
            nonlin = ALPHA * gp**2 / g
            friction = (2.0 / r) * fg * gp
            gpp_full = (Vp - nonlin - friction) / fg
            gpp = w * gpp_full + (1.0 - w) * gpp_sub

    return [gp, gpp]


# =====================================================================
# GENERIC SOLVER
# =====================================================================
def solve_generic(g0, ode_func, r_max=300.0, method='Radau',
                  rtol=1e-10, atol=1e-12, max_step=0.5, use_substrate_ic=False):
    """Solve any ODE variant. Returns (A_tail, B_tail, g_min) or nans."""
    if g0 <= 0.01:
        return np.nan, np.nan, np.nan

    r_start = 1e-6

    if use_substrate_ic:
        # Substrate IC: g''(0) = (1-g0)/3
        gpp0 = (1.0 - g0) / 3.0
    else:
        fg0 = f_kin(g0)
        if abs(fg0) < 1e-10:
            return np.nan, np.nan, np.nan
        gpp0 = g0**2 * (1.0 - g0) / (3.0 * fg0)

    g_ini = g0 + 0.5 * gpp0 * r_start**2
    gp_ini = gpp0 * r_start

    def event_blowup(r, y):
        return 500.0 - abs(y[0])
    event_blowup.terminal = True

    try:
        sol = solve_ivp(
            ode_func, [r_start, r_max], [g_ini, gp_ini],
            method=method, rtol=rtol, atol=atol,
            max_step=max_step, events=[event_blowup], dense_output=True
        )
    except Exception:
        return np.nan, np.nan, np.nan

    if sol.status == -1 or sol.t[-1] < r_max * 0.9:
        return np.nan, np.nan, np.nan

    g_min = np.min(sol.y[0])

    r_uni = np.linspace(100, min(sol.t[-1], 250), 5000)
    try:
        g_uni = sol.sol(r_uni)[0]
        u = (g_uni - 1.0) * r_uni
        design = np.column_stack([np.cos(r_uni), np.sin(r_uni)])
        BC = np.linalg.lstsq(design, u, rcond=None)[0]
        A = np.sqrt(BC[0]**2 + BC[1]**2)
        return A, BC[0], g_min
    except Exception:
        return np.nan, np.nan, g_min


def scan_and_find_fp(ode_func, g0_range, label, use_substrate_ic=False, max_step=0.5):
    """Scan A_tail and find phi-FP."""
    g0_scan = g0_range
    A_scan = np.full(len(g0_scan), np.nan)
    B_scan = np.full(len(g0_scan), np.nan)
    gmin_scan = np.full(len(g0_scan), np.nan)

    t0 = time.time()
    for i, g0 in enumerate(g0_scan):
        A, B, gm = solve_generic(g0, ode_func, use_substrate_ic=use_substrate_ic,
                                  max_step=max_step)
        A_scan[i] = A
        B_scan[i] = B
        gmin_scan[i] = gm
        if (i + 1) % 25 == 0:
            n_v = np.sum(~np.isnan(A_scan[:i+1]))
            print(f"    [{label}] {i+1}/{len(g0_scan)}: {n_v} valid, "
                  f"t={time.time()-t0:.0f}s")
            sys.stdout.flush()

    valid = ~np.isnan(A_scan) & (A_scan > 1e-8)
    n_valid = np.sum(valid)
    print(f"    [{label}] Scan done: {n_valid}/{len(g0_scan)} valid, "
          f"t={time.time()-t0:.0f}s")

    if n_valid < 10:
        print(f"    [{label}] Too few valid solutions!")
        return None

    f_Atail = interp1d(g0_scan[valid], A_scan[valid],
                        kind='cubic', fill_value=np.nan, bounds_error=False)

    g0_max_fp = np.max(g0_scan[valid]) / PHI - 0.01
    g0_min_fp = np.min(g0_scan[valid]) + 0.02

    if g0_max_fp <= g0_min_fp:
        print(f"    [{label}] Range too narrow for phi-FP search!")
        return {'A_scan': A_scan, 'B_scan': B_scan, 'gmin_scan': gmin_scan,
                'valid': valid, 'g0_scan': g0_scan, 'f_Atail': f_Atail,
                'fp_found': False}

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

    result = {'A_scan': A_scan, 'B_scan': B_scan, 'gmin_scan': gmin_scan,
              'valid': valid, 'g0_scan': g0_scan, 'f_Atail': f_Atail,
              'fp_found': False}

    if len(sign_changes) > 0:
        idx = sign_changes[0]
        try:
            g0_star = brentq(residual_fp, g0_finite[idx], g0_finite[idx + 1])
        except:
            print(f"    [{label}] brentq failed")
            return result

        A_e = f_Atail(g0_star)
        A_mu = f_Atail(PHI * g0_star)
        r21 = (A_mu / A_e)**4

        print(f"\n    [{label}] phi-FP FOUND:")
        print(f"      g0*       = {g0_star:.6f}")
        print(f"      phi*g0*   = {PHI*g0_star:.6f}")
        print(f"      A_e       = {A_e:.6f}")
        print(f"      A_mu      = {A_mu:.6f}")
        print(f"      (A_mu/A_e)^4 = {r21:.4f}  [PDG: {R21_PDG}]")
        print(f"      delta     = {100*(r21/R21_PDG - 1):.5f}%")

        result['fp_found'] = True
        result['g0_star'] = g0_star
        result['A_e'] = A_e
        result['A_mu'] = A_mu
        result['r21'] = r21

        # Max achievable ratio
        g0_hi = np.max(g0_scan[valid])
        A_hi = f_Atail(g0_hi)
        if not np.isnan(A_hi) and A_e > 0:
            max_r = (A_hi / A_e)**4
            result['max_ratio'] = max_r
            print(f"      max (A/A_e)^4 = {max_r:.1f} at g0={g0_hi:.3f}")
            print(f"      tau r31={R31_PDG:.1f}: "
                  f"{'REACHABLE' if max_r > R31_PDG else f'NOT reachable (gap {R31_PDG/max_r:.1f}x)'}")

            # Try to find tau candidate
            if max_r > R31_PDG:
                def res_r31(g0t):
                    At = f_Atail(g0t)
                    if np.isnan(At) or At <= 0:
                        return 1e10
                    return (At / A_e)**4 - R31_PDG

                g0_tau_range = np.linspace(PHI * g0_star + 0.01,
                                            g0_hi - 0.01, 500)
                res_tau = np.array([res_r31(g) for g in g0_tau_range])
                tau_fm = np.isfinite(res_tau) & (np.abs(res_tau) < 1e9)
                tau_sc = np.where(np.diff(np.sign(res_tau[tau_fm])))[0]
                if len(tau_sc) > 0:
                    g0_tau_f = g0_tau_range[tau_fm]
                    try:
                        g0_tau = brentq(res_r31, g0_tau_f[tau_sc[0]],
                                        g0_tau_f[tau_sc[0]+1])
                        A_tau = f_Atail(g0_tau)
                        r31 = (A_tau / A_e)**4
                        m_tau = M_E * r31
                        # Koide
                        k_num = M_E + M_MU + m_tau
                        k_den = (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(m_tau))**2
                        koide = k_num / k_den
                        delta_tau = 100 * (r31 / R31_PDG - 1)
                        print(f"\n      *** TAU CANDIDATE ***")
                        print(f"      g0^tau    = {g0_tau:.6f}")
                        print(f"      A_tau     = {A_tau:.6f}")
                        print(f"      r_31      = {r31:.2f} [PDG: {R31_PDG}]")
                        print(f"      delta_r31 = {delta_tau:.3f}%")
                        print(f"      m_tau     = {m_tau:.2f} MeV [PDG: {M_TAU}]")
                        print(f"      Koide     = {koide:.6f} [exact: 0.666667]")
                        result['g0_tau'] = g0_tau
                        result['r31'] = r31
                        result['m_tau'] = m_tau
                        result['koide'] = koide
                    except:
                        pass
    else:
        print(f"    [{label}] NO phi-FP found!")
        if np.any(finite_mask):
            print(f"      Residual range: [{np.nanmin(res_vals[finite_mask]):.1f}, "
                  f"{np.nanmax(res_vals[finite_mask]):.1f}]")

    return result


# =====================================================================
# MAIN
# =====================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("  TGP: UV Completion Paths for Tau Mass")
    print(f"  Ghost barrier at g0_crit ~ 1.63 (full ODE)")
    print(f"  Need: r_31 = {R31_PDG} (tau/electron mass ratio)")
    print("=" * 70)

    results = {}

    # ---- PATH A: Substrate ODE ----
    print("\n" + "=" * 70)
    print("  PATH A: Substrate ODE (ghost-free)")
    print("  g^2*g'' + g*(g')^2 + (2/r)*g^2*g' = V'(g)")
    print("=" * 70)

    g0_range_A = np.linspace(0.05, 5.0, 150)
    results['A'] = scan_and_find_fp(ode_substrate, g0_range_A, "SubODE",
                                     use_substrate_ic=True, max_step=0.5)

    # ---- PATH B: Running alpha ----
    print("\n" + "=" * 70)
    print("  PATH B: Running alpha(g) — alpha_UV < alpha_IR")
    print("  alpha = 2.0 (IR, g~1) -> alpha_UV (UV, g < 0.85)")
    print("=" * 70)

    for alpha_uv in [1.5, 1.0, 0.5]:
        lbl = f"RunA(uv={alpha_uv})"
        print(f"\n  --- alpha_UV = {alpha_uv} ---")
        g_star_uv = np.exp(-1.0 / (2 * alpha_uv))
        print(f"  g*(alpha_UV) = {g_star_uv:.4f}")

        ode_B = lambda r, y, auv=alpha_uv: ode_running_alpha(r, y, alpha_ir=2.0,
                                                               alpha_uv=auv)
        g0_range_B = np.linspace(g_star_uv + 0.05, 3.5, 100)
        results[f'B_{alpha_uv}'] = scan_and_find_fp(ode_B, g0_range_B, lbl,
                                                     max_step=0.1)

    # ---- PATH C: Blended ODE ----
    print("\n" + "=" * 70)
    print("  PATH C: Blended ODE (full->substrate near ghost)")
    print("=" * 70)

    g0_range_C = np.linspace(GSTAR + 0.03, 3.5, 120)
    results['C'] = scan_and_find_fp(ode_blended, g0_range_C, "Blend",
                                     max_step=0.1)

    # ---- SUMMARY ----
    print("\n" + "=" * 70)
    print("  SUMMARY: UV Completion Paths")
    print("=" * 70)
    print(f"\n  {'Path':<30s} {'phi-FP?':>8s} {'g0*':>8s} {'r21':>10s} "
          f"{'max r':>10s} {'tau?':>6s}")
    print("  " + "-" * 75)

    for key, label in [('A', 'Substrate ODE'),
                        ('B_1.5', 'Running alpha(UV=1.5)'),
                        ('B_1.0', 'Running alpha(UV=1.0)'),
                        ('B_0.5', 'Running alpha(UV=0.5)'),
                        ('C', 'Blended full/substrate')]:
        r = results.get(key)
        if r is None:
            print(f"  {label:<30s}   FAILED")
            continue
        if not r.get('fp_found', False):
            print(f"  {label:<30s}      NO")
            continue
        fp = "YES"
        g0s = f"{r['g0_star']:.4f}"
        r21 = f"{r['r21']:.2f}"
        mr = f"{r.get('max_ratio', np.nan):.0f}" if 'max_ratio' in r else "?"
        tau = "YES" if r.get('g0_tau') else "NO"
        print(f"  {label:<30s} {fp:>8s} {g0s:>8s} {r21:>10s} {mr:>10s} {tau:>6s}")

    # ---- PLOT ----
    print("\n  Generating plot...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('TGP: UV Completion Paths for Tau Mass\n'
                 f'Ghost barrier at g0_crit~1.63, need r_31={R31_PDG:.0f}',
                 fontsize=12, fontweight='bold')

    colors = {'A': 'blue', 'B_1.5': 'green', 'B_1.0': 'orange',
              'B_0.5': 'red', 'C': 'purple'}
    labels_map = {'A': 'Substrate', 'B_1.5': r'$\alpha_{UV}$=1.5',
                   'B_1.0': r'$\alpha_{UV}$=1.0', 'B_0.5': r'$\alpha_{UV}$=0.5',
                   'C': 'Blended'}

    # Plot 1: A_tail comparison
    ax = axes[0, 0]
    for key in ['A', 'B_1.0', 'C']:
        r = results.get(key)
        if r and np.any(r['valid']):
            ax.semilogy(r['g0_scan'][r['valid']], r['A_scan'][r['valid']],
                        '-', color=colors[key], lw=1.5, label=labels_map[key])
    ax.set_xlabel('$g_0$'); ax.set_ylabel('$A_{\\rm tail}$')
    ax.set_title('$A_{\\rm tail}(g_0)$ — all paths')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Plot 2: g_min comparison
    ax = axes[0, 1]
    for key in ['A', 'B_1.0', 'C']:
        r = results.get(key)
        if r and np.any(r['valid']):
            gm = r['gmin_scan'][r['valid']]
            ax.plot(r['g0_scan'][r['valid']], gm,
                    '-', color=colors[key], lw=1.5, label=labels_map[key])
    ax.axhline(GSTAR, color='red', ls='--', label=f'$g^*$={GSTAR:.3f}')
    ax.set_xlabel('$g_0$'); ax.set_ylabel('$g_{\\rm min}$')
    ax.set_title('$g_{\\rm min}(g_0)$ — ghost proximity')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Plot 3: mass ratios
    ax = axes[1, 0]
    for key in ['A', 'B_1.0', 'C']:
        r = results.get(key)
        if r and r.get('fp_found') and np.any(r['valid']):
            A_e = r['A_e']
            ratio = (r['A_scan'][r['valid']] / A_e)**4
            ax.semilogy(r['g0_scan'][r['valid']], ratio,
                        '-', color=colors[key], lw=1.5, label=labels_map[key])
    ax.axhline(R21_PDG, color='orange', ls='--', label=f'$r_{{21}}$={R21_PDG:.0f}')
    ax.axhline(R31_PDG, color='red', ls='--', label=f'$r_{{31}}$={R31_PDG:.0f}')
    ax.set_xlabel('$g_0$'); ax.set_ylabel('$(A/A_e)^4$')
    ax.set_title('Mass ratios — tau accessibility')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Plot 4: summary table
    ax = axes[1, 1]
    ax.axis('off')
    rows = []
    for key, label in [('A', 'Substrate'), ('B_1.0', 'Run. alpha(1.0)'),
                         ('C', 'Blended')]:
        r = results.get(key)
        if r and r.get('fp_found'):
            tau_str = f"{r['m_tau']:.0f}" if r.get('m_tau') else "N/A"
            rows.append([label, f"{r['g0_star']:.4f}", f"{r['r21']:.1f}",
                         f"{r.get('max_ratio', 0):.0f}", tau_str])
    if rows:
        table = ax.table(cellText=rows,
                          colLabels=['Path', 'g0*', 'r21', 'max r', 'm_tau'],
                          loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.0, 1.8)
    ax.set_title('UV Completion Summary', fontsize=12, pad=20)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'tau_uv_completion.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"  Saved: {outpath}")
    plt.close()

    print("\n  DONE.\n")
