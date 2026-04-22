#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ringdown_qnm_shooting.py -- TGP QNM via Shooting Method
=========================================================
Finds quasinormal mode frequencies by solving the eigenvalue problem
directly in the complex omega plane.

Problem:
    d^2 Phi/dr*^2 + [omega^2 - V_eff(r*)] Phi = 0

Boundary conditions for QNM:
    r* -> r*_min:  Phi -> 0   (evanescent into rising potential barrier)
    r* -> r*_max:  Phi ~ exp(+i k_inf r*)  (outgoing wave)
    where k_inf = sqrt(omega^2 - V_inf)

Method:
    1. Solve background profile f(r) via BVP
    2. Build V_eff(r*) on tortoise grid, interpolate with cubic spline
    3. For trial complex omega:
       a. Integrate ODE from r*_min outward with Phi(r*_min)=0, Phi'(r*_min)=1
       b. At r*_max, decompose into incoming + outgoing waves:
          Phi = A_out * e^{+ik r*} + A_in * e^{-ik r*}
       c. QNM condition: A_in(omega) = 0
    4. Find zeros of A_in(omega) in complex plane via Newton-Raphson

Reference: dodatekC_ringdown.tex, thm:ringdown
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import os
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===================================================================
# PARAMETERS
# ===================================================================
ALPHA = 2.0
SOURCE_SIGN = -1
X_MAX = 60.0
N_MESH = 800

STRONG_PARAMS = [
    {"bg": 0.01, "S":  5.0,  "sigma": 0.5, "label": "S=5"},
    {"bg": 0.01, "S": 10.0,  "sigma": 0.5, "label": "S=10"},
]

# ===================================================================
# 1. RADIAL PROFILE SOLVER
# ===================================================================

def source_gaussian(x, S, sigma):
    return S * np.exp(-x**2 / (2 * sigma**2)) / (2 * np.pi * sigma**2)**1.5

def ode_fun(x, y, bg, S, sigma):
    f = np.maximum(y[0], 1e-15)
    fp = y[1]
    x_safe = np.maximum(x, 1e-10)
    src = SOURCE_SIGN * source_gaussian(x, S, sigma)
    fpp = (-(2.0 / x_safe) * fp
           - ALPHA * fp**2 / f
           - bg * f**2 + bg * f**3
           + src)
    return np.vstack([fp, fpp])

def bc_fun(ya, yb):
    return np.array([ya[1], yb[0] - 1.0])

def solve_profile(bg, S, sigma):
    x_min = 0.01
    x_mesh = np.linspace(x_min, X_MAX, N_MESH)
    peak = S / (4 * np.pi * sigma)
    f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
    fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
    y_guess = np.vstack([f_guess, fp_guess])

    for tol, nodes in [(1e-6, 50000), (1e-4, 100000), (1e-3, 150000)]:
        sol = solve_bvp(
            fun=lambda x, y: ode_fun(x, y, bg, S, sigma),
            bc=bc_fun, x=x_mesh, y=y_guess, tol=tol,
            max_nodes=nodes, verbose=0,
        )
        if sol.success:
            break
        x_mesh = np.linspace(x_min, X_MAX, 2 * len(x_mesh))
        f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
        fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
        y_guess = np.vstack([f_guess, fp_guess])

    if not sol.success:
        return None

    x = sol.x
    f = sol.y[0]
    fp = sol.y[1]
    rhs = ode_fun(x, sol.y, bg, S, sigma)
    fpp = rhs[1]
    return x, f, fp, fpp

# ===================================================================
# 2. POTENTIAL AND TORTOISE
# ===================================================================

def compute_V_eff(x, f, fp, fpp, bg, ell=0):
    f_safe = np.maximum(f, 1e-15)
    x_safe = np.maximum(x, 1e-10)
    chi = f_safe
    bracket = (ell * (ell + 1) / x_safe**2
               + (9.0/4.0) * fpp / f_safe
               + (59.0/16.0) * (fp / f_safe)**2
               + 4.0 * fp / (x_safe * f_safe)
               - 2.0 * bg * chi
               + 3.0 * bg * chi**2)
    V = bracket / f_safe
    return V

def compute_tortoise(x, f):
    integrand = np.sqrt(np.maximum(f, 1e-15))
    r_star = np.zeros_like(x)
    for i in range(1, len(x)):
        r_star[i] = r_star[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (x[i] - x[i-1])
    return r_star

def build_V_spline(x, f, fp, fpp, bg, ell):
    """Build V_eff(r*) as a CubicSpline on tortoise coordinate."""
    V = compute_V_eff(x, f, fp, fpp, bg, ell=ell)
    r_star = compute_tortoise(x, f)

    # Clean data
    mask = np.isfinite(V) & np.isfinite(r_star)
    r_star_c = r_star[mask]
    V_c = V[mask]

    # Monotonic r*
    dr = np.diff(r_star_c)
    mono = np.concatenate([[True], dr > 1e-15])
    r_star_c = r_star_c[mono]
    V_c = V_c[mono]

    # Skip first few noisy points
    skip = min(3, len(r_star_c) // 10)
    r_star_c = r_star_c[skip:]
    V_c = V_c[skip:]

    if len(r_star_c) < 20:
        return None, None, None

    cs = CubicSpline(r_star_c, V_c, extrapolate=True)
    V_inf = V_c[-1]  # V at outer boundary (should be ~ gamma = bg)

    return cs, r_star_c, V_inf

# ===================================================================
# 3. SHOOTING METHOD
# ===================================================================

def shoot(omega, V_spline, r_star_min, r_star_max, V_inf):
    """
    Integrate the perturbation equation from r*_min to r*_max
    for complex frequency omega.

    Returns A_in (amplitude of incoming wave at r*_max).
    QNM condition: A_in = 0.

    The ODE is:
        Phi'' + [omega^2 - V(r*)] Phi = 0

    As a first-order system with complex Phi:
        y[0] = Phi,  y[1] = Phi'
        y[0]' = y[1]
        y[1]' = -(omega^2 - V(r*)) y[0]
    """
    omega2 = omega**2

    def ode(r_star, y):
        V = float(V_spline(r_star))
        return [y[1], -(omega2 - V) * y[0]]

    # Initial condition: Phi(r*_min) = 0, Phi'(r*_min) = 1
    # (evanescent/zero at inner barrier)
    y0 = [0.0 + 0.0j, 1.0 + 0.0j]

    # Integrate
    try:
        sol = solve_ivp(ode, [r_star_min, r_star_max], y0,
                        method='RK45', rtol=1e-8, atol=1e-10,
                        dense_output=False)
        if not sol.success:
            return None, None
    except Exception:
        return None, None

    Phi_end = sol.y[0, -1]
    dPhi_end = sol.y[1, -1]

    # Decompose into outgoing + incoming at r*_max:
    # Phi = A_out e^{+ik r*} + A_in e^{-ik r*}
    # Phi' = ik A_out e^{+ik r*} - ik A_in e^{-ik r*}
    # where k = sqrt(omega^2 - V_inf)

    k_inf = np.sqrt(omega2 - V_inf + 0j)
    # Choose branch with Re(k) > 0 or Im(k) > 0
    if k_inf.real < 0:
        k_inf = -k_inf

    if abs(k_inf) < 1e-15:
        return None, None

    r_end = r_star_max
    exp_plus = np.exp(1j * k_inf * r_end)
    exp_minus = np.exp(-1j * k_inf * r_end)

    # Solve 2x2 system:
    # Phi_end = A_out * exp_plus + A_in * exp_minus
    # dPhi_end = ik * A_out * exp_plus - ik * A_in * exp_minus
    A_out = (Phi_end * 1j * k_inf + dPhi_end) / (2 * 1j * k_inf * exp_plus)
    A_in = (Phi_end * 1j * k_inf - dPhi_end) / (2 * 1j * k_inf * exp_minus)

    return A_in, A_out

def qnm_residual(omega, V_spline, r_star_min, r_star_max, V_inf):
    """
    Returns |A_in(omega)|^2 -- to be minimized.
    QNM has A_in = 0.
    """
    A_in, A_out = shoot(omega, V_spline, r_star_min, r_star_max, V_inf)
    if A_in is None:
        return 1e20
    # Normalize by outgoing amplitude to avoid trivial solution
    if A_out is not None and abs(A_out) > 1e-30:
        return abs(A_in / A_out)**2
    return abs(A_in)**2

def find_qnm_newton(omega_guess, V_spline, r_star_min, r_star_max, V_inf,
                     max_iter=100, tol=1e-8, d_omega=1e-5):
    """
    Find QNM frequency via Newton-Raphson on A_in(omega) = 0
    in the complex omega plane.

    Uses finite-difference Jacobian.
    """
    omega = omega_guess

    for iteration in range(max_iter):
        A_in, A_out = shoot(omega, V_spline, r_star_min, r_star_max, V_inf)
        if A_in is None:
            return None, iteration, "integration failed"

        # Normalize
        if A_out is not None and abs(A_out) > 1e-30:
            F = A_in / A_out
        else:
            F = A_in

        if abs(F) < tol:
            return omega, iteration, "converged"

        # Numerical Jacobian (Wirtinger derivative for complex function)
        # dF/domega_re and dF/domega_im
        A_in_r, A_out_r = shoot(omega + d_omega, V_spline,
                                r_star_min, r_star_max, V_inf)
        A_in_i, A_out_i = shoot(omega + 1j * d_omega, V_spline,
                                r_star_min, r_star_max, V_inf)

        if A_in_r is None or A_in_i is None:
            return None, iteration, "Jacobian integration failed"

        if A_out_r is not None and abs(A_out_r) > 1e-30:
            F_r = A_in_r / A_out_r
        else:
            F_r = A_in_r
        if A_out_i is not None and abs(A_out_i) > 1e-30:
            F_i = A_in_i / A_out_i
        else:
            F_i = A_in_i

        # Complex derivative: dF/domega = (F(omega+d) - F(omega)) / d
        dF_domega = (F_r - F) / d_omega

        if abs(dF_domega) < 1e-30:
            return None, iteration, "zero derivative"

        # Newton step
        delta = F / dF_domega

        # Damped step for stability
        step_size = min(1.0, 0.5 / max(abs(delta) / (abs(omega) + 0.1), 1.0))
        omega_new = omega - step_size * delta

        # Enforce Im(omega) < 0 (damped modes)
        if omega_new.imag > 0:
            omega_new = omega_new.real - abs(omega_new.imag) * 1j

        omega = omega_new

    return omega, max_iter, "max iterations"

def scan_qnm(V_spline, r_star_min, r_star_max, V_inf,
             omega_re_range, omega_im_range, n_re=40, n_im=20):
    """
    Scan complex omega plane for QNM candidates.
    Returns heatmap of |A_in/A_out|^2 and candidate positions.
    """
    omega_re_arr = np.linspace(omega_re_range[0], omega_re_range[1], n_re)
    omega_im_arr = np.linspace(omega_im_range[0], omega_im_range[1], n_im)

    residual_map = np.full((n_im, n_re), np.nan)
    candidates = []

    for i, w_im in enumerate(omega_im_arr):
        for j, w_re in enumerate(omega_re_arr):
            omega = complex(w_re, w_im)
            res = qnm_residual(omega, V_spline, r_star_min, r_star_max, V_inf)
            if res < 1e19:
                residual_map[i, j] = np.log10(max(res, 1e-30))

    # Find local minima in the residual map
    for i in range(1, n_im - 1):
        for j in range(1, n_re - 1):
            val = residual_map[i, j]
            if np.isnan(val):
                continue
            neighbors = [residual_map[i+di, j+dj]
                         for di in [-1, 0, 1] for dj in [-1, 0, 1]
                         if (di, dj) != (0, 0) and not np.isnan(residual_map[i+di, j+dj])]
            if neighbors and val < min(neighbors) and val < -0.5:
                candidates.append((omega_re_arr[j], omega_im_arr[i], val))

    candidates.sort(key=lambda c: c[2])  # sort by residual (lower = better)
    return omega_re_arr, omega_im_arr, residual_map, candidates

# ===================================================================
# 4. MAIN
# ===================================================================

def main():
    print("=" * 70)
    print("  TGP RINGDOWN: QNM VIA SHOOTING METHOD")
    print("  Proper eigenvalue problem in complex omega plane")
    print("=" * 70)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    all_results = []

    for ps in STRONG_PARAMS:
        bg, S, sigma = ps["bg"], ps["S"], ps["sigma"]
        label = ps["label"]
        print(f"\n{'='*60}")
        print(f"  Profile: {label}")
        print(f"{'='*60}")

        out = solve_profile(bg, S, sigma)
        if out is None:
            print("    SKIPPED (solver failed)")
            continue
        x, f, fp, fpp = out
        print(f"    f_max = {f.max():.4f} at x = {x[np.argmax(f)]:.3f}")

        for ell in [0, 2]:
            print(f"\n  --- l = {ell} ---")

            V_spline, r_star_arr, V_inf = build_V_spline(x, f, fp, fpp, bg, ell)
            if V_spline is None:
                print("    SKIPPED (spline construction failed)")
                continue

            r_star_min = r_star_arr[0]
            # Truncate outer domain: once V ~ V_inf, integration is trivial
            # Find where |V - V_inf| < 0.001 * |V_inf| for the first time
            r_star_max_full = r_star_arr[-1]
            test_r = np.linspace(r_star_arr[0], r_star_arr[-1], 500)
            test_V = V_spline(test_r)
            V_inf_est = test_V[-1]
            for k in range(len(test_r) - 1, 0, -1):
                if abs(test_V[k] - V_inf_est) > 0.01:
                    r_star_max = min(test_r[k] + 10.0, r_star_max_full)
                    break
            else:
                r_star_max = r_star_max_full

            V_at_min = float(V_spline(r_star_min))
            V_at_max = float(V_spline(r_star_max))

            print(f"    r* range: [{r_star_min:.2f}, {r_star_max:.2f}]")
            print(f"    V(r*_min) = {V_at_min:.4f}")
            print(f"    V(r*_max) = {V_at_max:.6f} (should be ~ gamma = {bg})")
            print(f"    V_inf = {V_inf:.6f}")

            # For l=0, find where V transitions from negative to positive
            # and start integration from there
            r_start = r_star_min
            if ell == 0:
                # Find outermost V < 0 -> V > 0 transition
                test_r = np.linspace(r_star_min, r_star_max, 1000)
                test_V = V_spline(test_r)
                for k in range(len(test_r) - 1):
                    if test_V[k] < 0 and test_V[k+1] >= 0:
                        r_start = test_r[k+1]
                if r_start > r_star_min:
                    print(f"    l=0: V<0 region truncated, starting at r* = {r_start:.2f}")

            # --- STEP 1: Scan complex omega plane ---
            print(f"    Scanning complex omega plane ...")

            # Scan range: omega_re from sqrt(V_inf) to 2*sqrt(max(V))
            V_test = V_spline(np.linspace(r_start, r_star_max, 500))
            V_max_est = np.max(V_test[np.isfinite(V_test)])
            omega_re_lo = max(0.05, np.sqrt(max(V_inf, 0)))
            omega_re_hi = min(3.0, 2.0 * np.sqrt(max(V_max_est, 0.1)))
            omega_im_lo = -0.8  # upper bound on damping
            omega_im_hi = -0.001  # nearly undamped

            print(f"    Scan range: omega_re in [{omega_re_lo:.3f}, {omega_re_hi:.3f}]")
            print(f"                omega_im in [{omega_im_lo:.3f}, {omega_im_hi:.3f}]")

            omega_re_scan, omega_im_scan, res_map, candidates = scan_qnm(
                V_spline, r_start, r_star_max, V_inf,
                (omega_re_lo, omega_re_hi), (omega_im_lo, omega_im_hi),
                n_re=25, n_im=15
            )

            print(f"    Found {len(candidates)} QNM candidates")
            for k, (w_re, w_im, log_res) in enumerate(candidates[:5]):
                print(f"      candidate {k}: omega ~ {w_re:.4f} {w_im:+.4f}i  "
                      f"(log10|A_in/A_out|^2 = {log_res:.2f})")

            # --- STEP 2: Refine with Newton-Raphson ---
            qnm_found = []
            for w_re, w_im, _ in candidates[:6]:
                omega_guess = complex(w_re, w_im)
                omega_qnm, n_iter, status = find_qnm_newton(
                    omega_guess, V_spline, r_start, r_star_max, V_inf,
                    max_iter=80, tol=1e-10)

                if omega_qnm is not None and status == "converged":
                    # Verify: compute residual at found frequency
                    res_check = qnm_residual(omega_qnm, V_spline,
                                             r_start, r_star_max, V_inf)
                    if res_check < 1e-6:
                        Q = -omega_qnm.real / (2 * omega_qnm.imag) if omega_qnm.imag != 0 else float('inf')
                        print(f"    QNM FOUND: omega = {omega_qnm.real:.6f} "
                              f"{omega_qnm.imag:+.6f}i  Q = {Q:.2f}  "
                              f"(iter={n_iter}, res={res_check:.2e})")
                        qnm_found.append(omega_qnm)

            # Deduplicate (modes within 1% of each other)
            qnm_unique = []
            for w in qnm_found:
                is_dup = False
                for wu in qnm_unique:
                    if abs(w.real - wu.real) / max(abs(w.real), 0.01) < 0.01:
                        is_dup = True
                        break
                if not is_dup and w.real > 0.01:
                    qnm_unique.append(w)

            if not qnm_unique:
                print(f"    No converged QNM found.")
                # Try brute-force initial guesses
                print(f"    Trying brute-force initial guesses ...")
                for w_re_try in np.linspace(omega_re_lo, omega_re_hi, 6):
                    for w_im_try in [-0.05, -0.1, -0.3]:
                        omega_guess = complex(w_re_try, w_im_try)
                        omega_qnm, n_iter, status = find_qnm_newton(
                            omega_guess, V_spline, r_start, r_star_max, V_inf,
                            max_iter=60, tol=1e-8)
                        if omega_qnm is not None and status == "converged":
                            res_check = qnm_residual(omega_qnm, V_spline,
                                                     r_start, r_star_max, V_inf)
                            if res_check < 1e-4 and omega_qnm.real > 0.01:
                                is_dup = any(abs(omega_qnm.real - wu.real) / max(abs(omega_qnm.real), 0.01) < 0.02
                                             for wu in qnm_unique)
                                if not is_dup:
                                    Q = -omega_qnm.real / (2 * omega_qnm.imag) if omega_qnm.imag != 0 else float('inf')
                                    print(f"      QNM: omega = {omega_qnm.real:.6f} "
                                          f"{omega_qnm.imag:+.6f}i  Q = {Q:.2f}  "
                                          f"(res={res_check:.2e})")
                                    qnm_unique.append(omega_qnm)

            # Sort by real part
            qnm_unique.sort(key=lambda w: w.real)

            all_results.append({
                "label": label, "S": S, "ell": ell, "bg": bg,
                "r_star_min": r_start, "r_star_max": r_star_max,
                "V_inf": V_inf,
                "qnm_modes": qnm_unique,
                "scan_re": omega_re_scan, "scan_im": omega_im_scan,
                "scan_map": res_map,
                "candidates": candidates,
            })

    if not all_results:
        print("\nNo results.")
        return

    # =================================================================
    # PLOTS
    # =================================================================
    print(f"\n  Generating plots ...")

    # --- Fig 1: QNM scan heatmaps ---
    n_plots = len(all_results)
    n_cols = 2
    n_rows = (n_plots + 1) // 2
    fig1, axes1 = plt.subplots(n_rows, n_cols, figsize=(7*n_cols, 5*n_rows),
                                squeeze=False)
    fig1.suptitle("TGP QNM: Shooting Method Scan", fontsize=15, y=1.01)

    for idx, res in enumerate(all_results):
        ax = axes1[idx // n_cols, idx % n_cols]

        # Plot heatmap
        extent = [res["scan_re"][0], res["scan_re"][-1],
                  res["scan_im"][0], res["scan_im"][-1]]
        im = ax.imshow(res["scan_map"], origin='lower', aspect='auto',
                       extent=extent, cmap='viridis_r', vmin=-8, vmax=2)
        plt.colorbar(im, ax=ax, label=r"$\log_{10}|A_{\rm in}/A_{\rm out}|^2$")

        # Mark candidates
        for w_re, w_im, _ in res["candidates"][:5]:
            ax.plot(w_re, w_im, "wx", ms=8, mew=2)

        # Mark converged QNMs
        for w in res["qnm_modes"]:
            ax.plot(w.real, w.imag, "r*", ms=15, mew=1.5, zorder=10)

        # Mass gap line
        ax.axvline(np.sqrt(max(res["V_inf"], 0)), color="white", ls="--",
                   lw=0.8, alpha=0.7, label=r"$\omega_{\rm min}$")

        ax.set_xlabel(r"Re($\omega$)", fontsize=12)
        ax.set_ylabel(r"Im($\omega$)", fontsize=12)
        ax.set_title(f"{res['label']}, $\\ell = {res['ell']}$", fontsize=13)
        ax.legend(fontsize=8, loc='upper right')

    # Remove empty subplots
    for idx in range(n_plots, n_rows * n_cols):
        axes1[idx // n_cols, idx % n_cols].set_visible(False)

    fig1.tight_layout()
    path1 = os.path.join(save_dir, "ringdown_qnm_scan.png")
    fig1.savefig(path1, dpi=150, bbox_inches='tight')
    print(f"    Saved {path1}")

    plt.close("all")

    # =================================================================
    # SUMMARY
    # =================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY: QNM SHOOTING METHOD RESULTS")
    print("=" * 70)

    print(f"\n  {'Source':<10} {'l':<4} {'omega (QNM)':<28} {'Q':>8} {'Status'}")
    print("  " + "-" * 60)

    for res in all_results:
        if res["qnm_modes"]:
            for w in res["qnm_modes"]:
                Q = -w.real / (2 * w.imag) if w.imag != 0 else float('inf')
                print(f"  {res['label']:<10} {res['ell']:<4} "
                      f"{w.real:+.6f} {w.imag:+.6f}i      {Q:8.2f} CONVERGED")
        else:
            n_cand = len(res["candidates"])
            print(f"  {res['label']:<10} {res['ell']:<4} "
                  f"{'---':<28} {'---':>8} "
                  f"({n_cand} candidates, not converged)")

    bg0 = STRONG_PARAMS[0]["bg"]
    print(f"\n  Mass gap: omega_min = sqrt(gamma) = {np.sqrt(bg0):.4f}")
    vinf_strs = [f"{r['V_inf']:.6f}" for r in all_results]
    print(f"  V_inf values: {vinf_strs}")

    # l=0 vs l=2 comparison
    print(f"\n  --- l=0 vs l=2 COMPARISON ---")
    for ps in STRONG_PARAMS:
        l0 = [r for r in all_results if r["label"] == ps["label"] and r["ell"] == 0]
        l2 = [r for r in all_results if r["label"] == ps["label"] and r["ell"] == 2]
        if l0 and l2:
            n0 = len(l0[0]["qnm_modes"])
            n2 = len(l2[0]["qnm_modes"])
            print(f"  {ps['label']}: l=0 has {n0} modes, l=2 has {n2} modes")
            if n0 > 0 and n2 > 0:
                w0 = l0[0]["qnm_modes"][0]
                w2 = l2[0]["qnm_modes"][0]
                delta_re = abs(w0.real - w2.real) / max(w0.real, 0.01)
                delta_im = abs(w0.imag - w2.imag) / max(abs(w0.imag), 0.001)
                print(f"    Fundamental: l=0: {w0.real:.4f}{w0.imag:+.4f}i  "
                      f"l=2: {w2.real:.4f}{w2.imag:+.4f}i")
                print(f"    delta_omega_re/omega = {delta_re:.3f}, "
                      f"delta_gamma/gamma = {delta_im:.3f}")

    print()

if __name__ == "__main__":
    main()
