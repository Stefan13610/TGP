"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
radial_profile.py  --  Theory of Generated Space (TGP)
======================================================
Numerically solve the static, spherically-symmetric TGP field equation
for a single extended source (Gaussian mass distribution):

    f'' + (2/x)f' + alpha (f')^2 / f + beta_hat f^2 - gamma_hat f^3
        = S_hat * exp(-x^2 / (2 sigma^2)) / (2 pi sigma^2)^{3/2}

where f = Phi/Phi0, x = r/r0, and S_hat = q Phi0 r0^2 M encodes source
strength.

Boundary conditions (regular, no singularity):
    f'(0) = 0       (spherical symmetry)
    f(x_max) = 1    (far-field background)

This avoids the delta-function singularity entirely.  For sigma -> 0
the solution approaches the point-source limit.

Physical regimes visible in the solution:
    I   (large x):   f ~ 1 + C/x           (Newtonian / Yukawa tail)
    II  (medium x):  plateau / bump          (beta-repulsion)
    III (small x):   large f peak            (gamma-confinement)

Vacuum condition: beta_hat = gamma_hat ensures f = 1 is background solution.

Outputs (saved to scripts/plots/):
    radial_phi.png                  -- f(x) profiles
    radial_derivatives.png          -- f'(x) and f''(x)
    radial_force_decomposition.png  -- force term breakdown
"""

import os
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# ── Parameters ──────────────────────────────────────────────────────────────
ALPHA = 2.0

X_MAX = 40.0
N_MESH = 500
SIGMA = 0.5      # source width (dimensionless)

# Source sign: DERIVED from action variation (eq:variation-explicit, sek08).
#
# The action S = ∫d⁴x ψ⁴ [KE - U(ψ) - (q/Φ₀)ρψ] gives, after δS/δψ = 0:
#
#   ∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀² = -qΦ₀ρ
#
# The MINUS sign on the RHS means mass creates EXCESS Φ (space generation):
#   f = Φ/Φ₀ > 1 near mass → c(Φ) < c₀ (slower propagation in denser space).
# This is NOT a convention choice — it follows from the source term
# -(q/Φ₀)ρψ in the Lagrangian (line 1153 of sek08_formalizm.tex).
SOURCE_SIGN = -1   # from action: -qΦ₀ρ on RHS (eq:variation-explicit)

# Parameter sets (beta_hat = gamma_hat for vacuum consistency)
PARAM_SETS = [
    {"bg": 0.005, "S": 1.0,  "sigma": 0.5,
     "label": r"$\hat\beta{=}\hat\gamma{=}0.005,\;S{=}1$"},
    {"bg": 0.01,  "S": 1.0,  "sigma": 0.5,
     "label": r"$\hat\beta{=}\hat\gamma{=}0.01,\;S{=}1$"},
    {"bg": 0.05,  "S": 1.0,  "sigma": 0.5,
     "label": r"$\hat\beta{=}\hat\gamma{=}0.05,\;S{=}1$"},
    {"bg": 0.01,  "S": 5.0,  "sigma": 0.5,
     "label": r"$\hat\beta{=}\hat\gamma{=}0.01,\;S{=}5$"},
    {"bg": 0.01,  "S": 1.0,  "sigma": 1.5,
     "label": r"$\hat\beta{=}\hat\gamma{=}0.01,\;\sigma{=}1.5$"},
]


# ── Source function ─────────────────────────────────────────────────────────
def source(x, S, sigma):
    """Normalized Gaussian source: integral over 4 pi x^2 dx = S."""
    return S * np.exp(-x**2 / (2 * sigma**2)) / (2 * np.pi * sigma**2)**1.5


# ── ODE for solve_bvp ──────────────────────────────────────────────────────
def ode_fun(x, y, bg, S, sigma):
    """
    y[0] = f,  y[1] = f'.
    f'' = -(2/x)*f' - alpha*(f')^2/f - bg*f^2 + bg*f^3 + source
    """
    f = np.maximum(y[0], 1e-15)
    fp = y[1]
    x_safe = np.maximum(x, 1e-10)

    src = SOURCE_SIGN * source(x, S, sigma)
    fpp = -(2.0 / x_safe) * fp \
          - ALPHA * fp**2 / f \
          - bg * f**2 + bg * f**3 \
          + src

    return np.vstack([fp, fpp])


def bc_fun(ya, yb):
    """
    ya = y at x = x_min (near 0): f'(0) = 0
    yb = y at x = x_max:          f(x_max) = 1
    """
    return np.array([ya[1], yb[0] - 1.0])  # f'(0) = 0, f(xmax) = 1


# ── Analytical Jacobian for better convergence ─────────────────────────────
def ode_jac(x, y, bg):
    """df/dy Jacobian (2x2 at each x) for the ODE."""
    f = np.maximum(y[0], 1e-15)
    fp = y[1]
    x_safe = np.maximum(x, 1e-10)

    # d(fpp)/d(f)
    dfpp_df = ALPHA * fp**2 / f**2 - 2 * bg * f + 3 * bg * f**2
    # d(fpp)/d(fp)
    dfpp_dfp = -2.0 / x_safe - 2 * ALPHA * fp / f

    n = x.shape[0]
    jac = np.zeros((2, 2, n))
    jac[0, 1, :] = 1.0
    jac[1, 0, :] = dfpp_df
    jac[1, 1, :] = dfpp_dfp
    return jac


# ── Solver ──────────────────────────────────────────────────────────────────
def solve_profile(bg, S, sigma, x_max=X_MAX, n_mesh=N_MESH):
    """Solve the radial BVP with Gaussian source. Returns x, f, fp, fpp."""
    # Use x_min > 0 to avoid 2/x singularity; f'(x_min) ≈ 0 by symmetry
    x_min = 0.01
    x_mesh = np.linspace(x_min, x_max, n_mesh)

    # Initial guess: 1 + Gaussian bump
    peak = S / (4 * np.pi * sigma)  # rough estimate of central enhancement
    f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
    fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
    y_guess = np.vstack([f_guess, fp_guess])

    sol = solve_bvp(
        fun=lambda x, y: ode_fun(x, y, bg, S, sigma),
        bc=bc_fun,
        x=x_mesh,
        y=y_guess,
        tol=1e-6,
        max_nodes=30000,
        verbose=0,
    )

    if not sol.success:
        print(f"    [warn] BVP: {sol.message}")
        # Retry with denser mesh and relaxed tol
        x_mesh2 = np.linspace(x_min, x_max, 3 * n_mesh)
        f_guess2 = 1.0 + peak * np.exp(-x_mesh2**2 / (2 * sigma**2))
        fp_guess2 = -peak * x_mesh2 / sigma**2 * np.exp(-x_mesh2**2 / (2 * sigma**2))
        y_guess2 = np.vstack([f_guess2, fp_guess2])
        sol = solve_bvp(
            fun=lambda x, y: ode_fun(x, y, bg, S, sigma),
            bc=bc_fun,
            x=x_mesh2,
            y=y_guess2,
            tol=1e-4,
            max_nodes=60000,
            verbose=0,
        )
        if not sol.success:
            print(f"    [FAIL] BVP: {sol.message}")
            return None

    x = sol.x
    f = sol.y[0]
    fp = sol.y[1]

    # Reconstruct f''
    rhs = ode_fun(x, sol.y, bg, S, sigma)
    fpp = rhs[1]

    return x, f, fp, fpp


# ── Inflection points ──────────────────────────────────────────────────────
def find_inflection_points(x, fpp):
    sc = np.where(np.diff(np.sign(fpp)))[0]
    pts = []
    for i in sc:
        w = -fpp[i] / (fpp[i + 1] - fpp[i] + 1e-40)
        pts.append(x[i] + w * (x[i + 1] - x[i]))
    return np.array(pts)


# ── Plotting ────────────────────────────────────────────────────────────────
def plot_all(results, save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(results)))

    # ── Fig 1: profiles ──
    fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(14, 6))

    for (x, f, fp, fpp, meta), col in zip(results, colors):
        ax1a.plot(x, f, color=col, lw=1.8, label=meta["label"])
        ax1b.plot(x, f - 1, color=col, lw=1.8, label=meta["label"])

    ax1a.axhline(1, color="k", ls="--", lw=0.8)
    ax1a.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax1a.set_ylabel(r"$f(x) = \Phi/\Phi_0$", fontsize=13)
    ax1a.set_title("TGP radial field profile", fontsize=14)
    ax1a.legend(fontsize=8)
    ax1a.grid(True, ls=":", alpha=0.4)

    # Log-log of perturbation
    for (x, f, fp, fpp, meta), col in zip(results, colors):
        delta = f - 1
        pos = delta > 0
        if pos.any():
            ax1b.plot(x[pos], delta[pos], color=col, lw=1.5)
    x_ref = np.linspace(2, X_MAX, 200)
    ax1b.plot(x_ref, 0.5 / x_ref, "k:", lw=0.8, label=r"$\propto 1/x$")
    ax1b.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax1b.set_ylabel(r"$\delta(x) = f - 1$", fontsize=13)
    ax1b.set_title(r"Perturbation $\delta\Phi/\Phi_0$ (log-log)", fontsize=14)
    ax1b.set_xscale("log")
    ax1b.set_yscale("log")
    ax1b.legend(fontsize=8)
    ax1b.grid(True, which="both", ls=":", alpha=0.4)

    fig1.tight_layout()
    p = os.path.join(save_dir, "radial_phi.png")
    fig1.savefig(p, dpi=180)
    print(f"  Saved {p}")

    # ── Fig 2: derivatives ──
    fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for (x, f, fp, fpp, meta), col in zip(results, colors):
        ax2a.plot(x, fp, color=col, lw=1.5, label=meta["label"])
        ax2b.plot(x, fpp, color=col, lw=1.5, label=meta["label"])
    ax2a.axhline(0, color="k", lw=0.5)
    ax2a.set_ylabel(r"$f'(x)$", fontsize=13)
    ax2a.set_title("Derivatives of TGP field", fontsize=14)
    ax2a.legend(fontsize=8)
    ax2a.grid(True, ls=":", alpha=0.4)
    ax2b.axhline(0, color="k", lw=0.5)
    ax2b.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax2b.set_ylabel(r"$f''(x)$", fontsize=13)
    ax2b.legend(fontsize=8)
    ax2b.grid(True, ls=":", alpha=0.4)
    fig2.tight_layout()
    p = os.path.join(save_dir, "radial_derivatives.png")
    fig2.savefig(p, dpi=180)
    print(f"  Saved {p}")

    # ── Fig 3: force decomposition (first result) ──
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    x, f, fp, fpp, meta = results[0]
    bg = meta["bg"]
    S, sig = meta["S"], meta["sigma"]

    x_s = np.maximum(x, 1e-10)
    f_s = np.maximum(f, 1e-15)
    F_lapl = -(2.0 / x_s) * fp
    F_nl = -ALPHA * fp**2 / f_s
    F_self = bg * f**2 * (f - 1)   # -bg*f^2 + bg*f^3 combined
    F_src = SOURCE_SIGN * source(x, S, sig)

    ax3.plot(x, F_lapl, "b-", lw=1.2, label=r"$-(2/x)f'$")
    ax3.plot(x, F_nl, "g-", lw=1.2, label=r"$-\alpha(f')^2/f$")
    ax3.plot(x, F_self, "r-", lw=1.2, label=r"$\hat\beta f^2(f{-}1)$")
    ax3.plot(x, F_src, "m--", lw=1.5, label="source")
    ax3.plot(x, fpp, "k:", lw=1.8, label=r"$f''$")
    ax3.axhline(0, color="k", lw=0.5)
    ax3.set_xlabel(r"$x = r/r_0$", fontsize=13)
    ax3.set_ylabel("Force contribution", fontsize=13)
    ax3.set_title(f"Force decomposition ({meta['label']})", fontsize=13)
    ax3.legend(fontsize=10)
    ax3.set_xlim(0, 15)
    ax3.grid(True, ls=":", alpha=0.4)
    fig3.tight_layout()
    p = os.path.join(save_dir, "radial_force_decomposition.png")
    fig3.savefig(p, dpi=180)
    print(f"  Saved {p}")

    plt.close("all")


# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("TGP radial profile solver  (BVP, Gaussian source)")
    print("=" * 60)

    results = []
    for ps in PARAM_SETS:
        bg, S, sig = ps["bg"], ps["S"], ps["sigma"]
        print(f"\n  Solving: bg={bg}, S={S}, sigma={sig} ...")
        out = solve_profile(bg, S, sig)
        if out is None:
            print("    SKIPPED")
            continue
        x, f, fp, fpp = out
        results.append((x, f, fp, fpp, ps))

        # Diagnostics
        print(f"    f_max = {f.max():.6f}  at x = {x[np.argmax(f)]:.3f}")
        print(f"    f(0)  = {f[0]:.6f}")

        x_infl = find_inflection_points(x, fpp)
        if len(x_infl):
            print(f"    Inflection points: x = {np.round(x_infl[:5], 3)}")

        # Far-field Newtonian coefficient
        mask = (x > 0.5 * X_MAX) & np.isfinite(f)
        if mask.sum() > 5:
            C_fit = np.mean((f[mask] - 1.0) * x[mask])
            print(f"    Far-field C ~ {C_fit:.5f}  (f ~ 1 + C/x)")

    if results:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
        plot_all(results, save_dir=save_dir)
    else:
        print("\n  No successful solutions.")

    print("\nDone.")


if __name__ == "__main__":
    main()
