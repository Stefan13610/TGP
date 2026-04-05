#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ringdown_qnm.py -- TGP Quasinormal Mode Spectrum
=================================================
Computes the QNM spectrum for scalar (breathing) perturbations
around a static, spherically symmetric TGP background.

Steps:
  1. Solve strong-field radial profile f(r) via BVP
  2. Compute effective potential V_eff(r) from thm:ringdown (dodatekC)
  3. Find QNM frequencies via WKB approximation (Schutz-Will method)
  4. Compare with GR Schwarzschild QNM

Reference: dodatekC_ringdown.tex, eq:V-eff
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import os
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===================================================================
# PARAMETERS
# ===================================================================
ALPHA = 2.0

# Strong-field parameter sets: increasing source strength S
# beta_hat = gamma_hat = bg (vacuum condition)
# Large S => large f(0) => strong-field regime
STRONG_PARAMS = [
    {"bg": 0.01, "S":  2.0,  "sigma": 0.5, "label": "S=2 (moderate)"},
    {"bg": 0.01, "S":  5.0,  "sigma": 0.5, "label": "S=5 (strong)"},
    {"bg": 0.01, "S": 10.0,  "sigma": 0.5, "label": "S=10 (very strong)"},
    {"bg": 0.01, "S": 20.0,  "sigma": 0.5, "label": "S=20 (ultra strong)"},
]

X_MAX = 60.0
N_MESH = 800
SOURCE_SIGN = -1  # from action (eq:variation-explicit)

# ===================================================================
# 1. RADIAL PROFILE SOLVER (adapted from radial_profile.py)
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
        # densify mesh for retry
        x_mesh = np.linspace(x_min, X_MAX, 2 * len(x_mesh))
        f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
        fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
        y_guess = np.vstack([f_guess, fp_guess])

    if not sol.success:
        print(f"    [FAIL] BVP did not converge: {sol.message}")
        return None

    x = sol.x
    f = sol.y[0]
    fp = sol.y[1]
    rhs = ode_fun(x, sol.y, bg, S, sigma)
    fpp = rhs[1]
    return x, f, fp, fpp

# ===================================================================
# 2. EFFECTIVE POTENTIAL (from thm:ringdown, eq:V-eff)
# ===================================================================

def compute_V_eff(x, f, fp, fpp, bg, ell=0):
    """
    V_eff(r) = (c0^2 Phi0 / f) * [
        ell(ell+1)/r^2
        + (9/4) f''/f
        + (59/16) (f'/f)^2
        + 4 f'/(r f)
        - 2 beta chi + 3 gamma chi^2
    ]

    In dimensionless units (c0=1, Phi0=1, r -> x = r/r0):
    V_hat(x) = (1/f) * [
        ell(ell+1)/x^2
        + (9/4) fpp/f
        + (59/16) (fp/f)^2
        + 4 fp/(x f)
        - 2 bg f + 3 bg f^2
    ]

    (V_hat = V_eff * r0^2 / c0^2)
    """
    f_safe = np.maximum(f, 1e-15)
    x_safe = np.maximum(x, 1e-10)
    chi = f_safe  # since Phi0 = 1 in these units

    bracket = (ell * (ell + 1) / x_safe**2
               + (9.0/4.0) * fpp / f_safe
               + (59.0/16.0) * (fp / f_safe)**2
               + 4.0 * fp / (x_safe * f_safe)
               - 2.0 * bg * chi
               + 3.0 * bg * chi**2)

    V = bracket / f_safe  # the (Phi0/f) prefactor, with Phi0=1
    return V

# ===================================================================
# 3. TORTOISE COORDINATE
# ===================================================================

def compute_tortoise(x, f):
    """
    r* = integral sqrt(f) dx  (with c0=1, Phi0=1)
    dr*/dx = sqrt(f)
    """
    integrand = np.sqrt(np.maximum(f, 1e-15))
    r_star = np.zeros_like(x)
    for i in range(1, len(x)):
        r_star[i] = r_star[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (x[i] - x[i-1])
    return r_star

# ===================================================================
# 4. WKB QUASINORMAL MODE FREQUENCIES (Schutz-Will, 1st order)
# ===================================================================

def find_V_peak(x, V):
    """Find position and value of potential maximum."""
    # Smooth out numerical noise
    idx_valid = np.isfinite(V)
    x_v, V_v = x[idx_valid], V[idx_valid]
    if len(x_v) < 10:
        return None, None, None

    # Find global max (excluding boundaries)
    margin = max(5, len(x_v) // 20)
    idx_peak = margin + np.argmax(V_v[margin:-margin])
    x_peak = x_v[idx_peak]
    V_peak = V_v[idx_peak]

    # Second derivative at peak (for WKB)
    if idx_peak < 2 or idx_peak >= len(x_v) - 2:
        return x_peak, V_peak, None
    dx = x_v[idx_peak+1] - x_v[idx_peak]
    V_pp = (V_v[idx_peak+1] - 2*V_v[idx_peak] + V_v[idx_peak-1]) / dx**2

    return x_peak, V_peak, V_pp

def wkb_qnm(V_peak, V_pp_peak, n=0):
    """
    WKB formula (Schutz-Will, 1st order):
      omega^2 = V_peak - i*(n + 1/2) * sqrt(-2 V_pp_peak)

    Returns complex omega.
    For massive field: omega^2 = V_peak - mass_gap, but V_peak
    already includes the mass term.
    """
    if V_pp_peak is None or V_pp_peak >= 0:
        return None  # no barrier => no QNM

    omega_sq = V_peak - 1j * (n + 0.5) * np.sqrt(-2.0 * V_pp_peak)
    omega = np.sqrt(omega_sq)
    # Convention: positive real part, negative imaginary part
    if omega.real < 0:
        omega = -omega
    return omega

# ===================================================================
# 5. GR SCHWARZSCHILD QNM (for comparison)
# ===================================================================

def schwarzschild_qnm_l0(n=0):
    """
    Schwarzschild QNM for scalar field (l=0).
    Using tabulated values (Berti, Cardoso, Starinets 2009):
    omega * M = 0.1105 - 0.1049 i  (n=0, l=0, scalar)
    In units of 1/M.
    """
    # l=0, n=0 fundamental mode
    return complex(0.1105, -0.1049)

def schwarzschild_qnm_l2(n=0):
    """l=2 tensor mode (GR gravitational wave ringdown)"""
    return complex(0.3737, -0.0890)

# ===================================================================
# 6. MAIN: COMPUTE AND PLOT
# ===================================================================

def main():
    print("=" * 65)
    print("  TGP RINGDOWN: QUASINORMAL MODE COMPUTATION")
    print("=" * 65)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    results = []

    for ps in STRONG_PARAMS:
        bg, S, sigma = ps["bg"], ps["S"], ps["sigma"]
        label = ps["label"]
        print(f"\n  Solving profile: {label} ...")

        out = solve_profile(bg, S, sigma)
        if out is None:
            print("    SKIPPED (solver failed)")
            continue

        x, f, fp, fpp = out
        print(f"    f_max = {f.max():.4f} at x = {x[np.argmax(f)]:.3f}")
        print(f"    f(x_min) = {f[0]:.4f}")
        print(f"    chi_max = {f.max():.4f} (dimensionless)")

        # Compute V_eff for l=0 and l=2
        for ell in [0, 2]:
            V = compute_V_eff(x, f, fp, fpp, bg, ell=ell)
            r_star = compute_tortoise(x, f)

            # Find peak
            x_peak, V_peak, V_pp = find_V_peak(x, V)

            if x_peak is not None and V_peak is not None:
                print(f"    l={ell}: V_peak = {V_peak:.6f} at x = {x_peak:.3f}")
                if V_pp is not None:
                    print(f"           V'' at peak = {V_pp:.6f}")

                # WKB QNM
                omega = wkb_qnm(V_peak, V_pp, n=0)
                if omega is not None:
                    print(f"           omega_0 = {omega.real:.6f} + {omega.imag:.6f}i")
                    print(f"           f_ring = {omega.real/(2*np.pi):.6f} (dimensionless)")
                    print(f"           Q = {-omega.real/(2*omega.imag):.2f} (quality factor)")
                else:
                    print(f"           No QNM (no barrier or V'' >= 0)")

            results.append({
                "x": x, "f": f, "fp": fp, "fpp": fpp,
                "V": V, "r_star": r_star,
                "ell": ell, "bg": bg, "S": S, "label": label,
                "V_peak": V_peak, "V_pp": V_pp,
                "x_peak": x_peak,
            })

    if not results:
        print("\nNo solutions found.")
        return

    # ── MASS GAP ──
    bg0 = STRONG_PARAMS[0]["bg"]
    print(f"\n  Mass gap: omega_min^2 = gamma_hat = {bg0}")
    print(f"           omega_min   = {np.sqrt(bg0):.6f}")
    print(f"  (Modes with omega < omega_min are evanescent)")

    # ── PLOTS ──
    print(f"\n  Generating plots ...")

    # --- Fig 1: Profiles ---
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    seen_labels = set()
    for res in results:
        if res["ell"] != 0:
            continue
        lbl = res["label"]
        if lbl in seen_labels:
            continue
        seen_labels.add(lbl)
        ax1.plot(res["x"], res["f"], lw=1.8, label=lbl)
    ax1.axhline(1, color="k", ls="--", lw=0.8, label=r"$f = 1$ (vacuum)")
    ax1.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax1.set_ylabel(r"$f(x) = \Phi/\Phi_0$", fontsize=13)
    ax1.set_title("TGP strong-field radial profiles", fontsize=14)
    ax1.legend(fontsize=10)
    ax1.set_xlim(0, 30)
    ax1.grid(True, ls=":", alpha=0.4)
    fig1.tight_layout()
    path1 = os.path.join(save_dir, "ringdown_profiles.png")
    fig1.savefig(path1, dpi=150)
    print(f"    Saved {path1}")

    # --- Fig 2: V_eff(x) for different S, l=0 ---
    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

    for res in results:
        if res["ell"] != 0:
            continue
        V_plot = np.clip(res["V"], -5, 5)  # clip for display
        ax2a.plot(res["x"], V_plot, lw=1.5, label=res["label"])
        if res["V_peak"] is not None:
            ax2a.plot(res["x_peak"], res["V_peak"], "o", ms=8, color="red")

    ax2a.axhline(bg0, color="gray", ls="--", lw=0.8,
                 label=r"$\gamma$ (mass gap)")
    ax2a.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax2a.set_ylabel(r"$\hat{V}_{\rm eff}(x)$", fontsize=13)
    ax2a.set_title(r"Effective potential $V_{\rm eff}$ ($\ell=0$)", fontsize=14)
    ax2a.set_xlim(0, 25)
    ax2a.set_ylim(-0.5, 2.0)
    ax2a.legend(fontsize=9)
    ax2a.grid(True, ls=":", alpha=0.4)

    # --- Fig 2b: V_eff in tortoise coordinate ---
    for res in results:
        if res["ell"] != 0:
            continue
        V_plot = np.clip(res["V"], -5, 5)
        ax2b.plot(res["r_star"], V_plot, lw=1.5, label=res["label"])

    ax2b.axhline(bg0, color="gray", ls="--", lw=0.8)
    ax2b.set_xlabel(r"$r_*$ (tortoise)", fontsize=13)
    ax2b.set_ylabel(r"$\hat{V}_{\rm eff}(r_*)$", fontsize=13)
    ax2b.set_title(r"$V_{\rm eff}$ in tortoise coordinate", fontsize=14)
    ax2b.set_ylim(-0.5, 2.0)
    ax2b.legend(fontsize=9)
    ax2b.grid(True, ls=":", alpha=0.4)

    fig2.tight_layout()
    path2 = os.path.join(save_dir, "ringdown_V_eff.png")
    fig2.savefig(path2, dpi=150)
    print(f"    Saved {path2}")

    # --- Fig 3: V_eff for l=0 and l=2 (one S value) ---
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    target_S = 10.0
    for res in results:
        if res["S"] != target_S:
            continue
        V_plot = np.clip(res["V"], -5, 5)
        ax3.plot(res["x"], V_plot, lw=2.0,
                 label=rf"$\ell={res['ell']}$ ({res['label']})")
    ax3.axhline(bg0, color="gray", ls="--", lw=0.8,
                label=r"$\gamma$ (mass gap)")
    ax3.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax3.set_ylabel(r"$\hat{V}_{\rm eff}(x)$", fontsize=13)
    ax3.set_title(rf"$V_{{\rm eff}}$ for $S={target_S}$: multipoles", fontsize=14)
    ax3.set_xlim(0, 25)
    ax3.set_ylim(-0.5, 3.0)
    ax3.legend(fontsize=11)
    ax3.grid(True, ls=":", alpha=0.4)
    fig3.tight_layout()
    path3 = os.path.join(save_dir, "ringdown_V_multipoles.png")
    fig3.savefig(path3, dpi=150)
    print(f"    Saved {path3}")

    plt.close("all")

    # ── SUMMARY ──
    print("\n" + "=" * 65)
    print("  SUMMARY: TGP vs GR QNM")
    print("=" * 65)

    print(f"\n  GR Schwarzschild scalar l=0: omega*M = {schwarzschild_qnm_l0()}")
    print(f"  GR Schwarzschild tensor l=2: omega*M = {schwarzschild_qnm_l2()}")

    print(f"\n  TGP QNM (dimensionless, units of 1/r0):")
    for res in results:
        if res["V_peak"] is None or res["V_pp"] is None:
            continue
        omega = wkb_qnm(res["V_peak"], res["V_pp"], n=0)
        if omega is not None:
            Q = -omega.real / (2 * omega.imag)
            print(f"    {res['label']:30s} l={res['ell']}: "
                  f"omega = {omega.real:+.5f} {omega.imag:+.5f}i  "
                  f"Q = {Q:.1f}")

    print(f"""
  KEY STRUCTURAL DIFFERENCES:

  1. TGP has MASS GAP: omega_min = sqrt(gamma) = {np.sqrt(bg0):.4f}
     GR has NO mass gap (V -> 0 at infinity)

  2. TGP potential RISES toward center (c -> 0 trapping)
     GR potential has peak then DROPS to 0 at horizon

  3. TGP ringdown is SCALAR (breathing mode)
     GR ringdown is TENSOR (h+, hx polarizations)

  4. TGP quality factor Q depends on f_max (compactness)
     GR Q depends only on M and spin
""")

if __name__ == "__main__":
    main()
