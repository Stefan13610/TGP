#!/usr/bin/env python3
"""
tgp_soliton_bvp.py -- TGP BH soliton via boundary value problem
=================================================================

The TGP BH soliton is a separatrix solution of:
    psi'' + (2/r)psi' + 2(psi')^2/psi + beta*psi^2*(1-psi) = 0

with BCs:
    psi(r_min) >> 1   (frozen interior, psi ~ A/r)
    psi(r_max) = 1    (vacuum)

This is solved using scipy.integrate.solve_bvp (collocation/relaxation).

The substitution u = psi^3 linearizes the kinetic operator:
    D_kin[psi] = psi'' + 2(psi')^2/psi = (1/(3 psi^2)) (u'' where u=psi^3)

Actually that gives: psi'' + 2(psi')^2/psi = (1/(3psi^2)) u''
So: (1/(3psi^2)) u'' + (2/r)psi' + beta*psi^2(1-psi) = 0

In terms of u = psi^3:
    u' = 3 psi^2 psi'
    psi' = u'/(3 psi^2) = u'/(3 u^{2/3})

Hmm, this doesn't simplify nicely for the full equation. Let's stick with (psi, psi').

REGULARIZATION: To handle psi -> inf at r=0, we use the substitution
    w = 1/psi  (w -> 0 at center, w -> 1 at infinity)

Then: psi = 1/w, psi' = -w'/w^2, psi'' = -w''/w^2 + 2(w')^2/w^3

The equation becomes (after substituting):
    [-w''/w^2 + 2(w')^2/w^3] + (2/r)[-w'/w^2] + 2[-w'/w^2]^2/(1/w) + beta/w^2*(1-1/w) = 0

    -w''/w^2 + 2(w')^2/w^3 - 2w'/(r w^2) + 2(w')^2/(w^3) + beta(1-1/w)/w^2 = 0

    -w'' + 2(w')^2/w - 2w'/r + 2(w')^2/w + beta*(w - 1)/w = 0

Wait let me redo. psi = 1/w:
    psi' = -w'/w^2
    psi'' = -w''/w^2 + 2(w')^2/w^3

    Term 1: psi'' = -w''/w^2 + 2(w')^2/w^3
    Term 2: (2/r)psi' = -(2/r) w'/w^2
    Term 3: 2(psi')^2/psi = 2(w')^2/(w^4) * w = 2(w')^2/w^3
    Term 4: beta*psi^2*(1-psi) = beta*(1/w^2)*(1 - 1/w) = beta*(1/w^2 - 1/w^3)

    Sum: -w''/w^2 + 2(w')^2/w^3 - (2/r)w'/w^2 + 2(w')^2/w^3 + beta/w^2 - beta/w^3 = 0

    Multiply by -w^2:
    w'' - 2(w')^2/w + (2/r)w' - 2(w')^2/w - beta + beta/w = 0

    w'' - 4(w')^2/w + (2/r)w' - beta + beta/w = 0

    w'' = 4(w')^2/w - (2/r)w' + beta - beta/w
        = 4(w')^2/w - (2/r)w' + beta*(1 - 1/w)

This is the REGULARIZED equation for w = 1/psi.

BCs:
    w(r_min) = 0+  (frozen interior, psi -> inf)
    Actually w(0) = 0, but we need w'(0) too.

    From psi ~ A/r: w = r/A, w' = 1/A, w'' = 0
    Check: w'' = 0, 4(w')^2/w = 4/(A^2) * (A/r) = 4/(A*r)
    (2/r)w' = 2/(A*r)
    beta*(1-1/w) = beta*(1 - A/r) ~ -beta*A/r for small r

    Sum: 4/(Ar) - 2/(Ar) - beta*A/r = (2/A - beta*A)/r = 0
    => 2/A = beta*A => A^2 = 2/beta ✓

So w ~ r/A with A = sqrt(2/beta) near center, and w'' = 0 at leading order.

BCs for the BVP:
    w(r_min) = r_min / A  (approximately)
    w(r_max) = 1           (vacuum: psi = 1 => w = 1)

Author: Mateusz Serafin (with Claude)
Date: April 2026
"""

import numpy as np
from scipy.integrate import solve_bvp
import sys


# ============================================================
# 1. Regularized equation in w = 1/psi
# ============================================================

def ode_w(r, y, beta):
    """
    ODE system for w = 1/psi:
        w'' = 4(w')^2/w - (2/r)w' + beta*(1 - 1/w)

    y[0] = w, y[1] = w'
    """
    w = y[0]
    dw = y[1]

    # Safety: w must be positive
    w = np.maximum(w, 1e-30)

    ddw = 4.0 * dw**2 / w - (2.0/r) * dw + beta * (1.0 - 1.0/w)

    return np.vstack([dw, ddw])


def ode_w_selfconsistent(r, y, beta):
    """
    Self-consistent equation in w = 1/psi.

    The flat equation has potential beta*psi^2*(1-psi).
    Self-consistent has beta*psi^3*(1-psi).

    In w: beta*psi^3*(1-psi) = beta*(1/w^3)*(1-1/w) = beta*(1/w^3 - 1/w^4)

    Substituting into the full equation (multiply by -w^2):
    w'' = 4(w')^2/w - (2/r)w' + beta/w - beta/w^2

    Compare flat: w'' = 4(w')^2/w - (2/r)w' + beta - beta/w
    """
    w = y[0]
    dw = y[1]
    w = np.maximum(w, 1e-30)

    ddw = 4.0 * dw**2 / w - (2.0/r) * dw + beta / w - beta / w**2

    return np.vstack([dw, ddw])


# ============================================================
# 2. Boundary conditions
# ============================================================

def bc_flat(ya, yb, A, beta):
    """
    BCs for the flat equation:
        w(r_min) = r_min / A  (inner, frozen interior)
        w(r_max) = 1          (outer, vacuum)

    ya = y at r_min, yb = y at r_max
    """
    r_min = 0.01  # will be set by the mesh
    # Actually, ya[0] should be close to r_min/A but we let the solver adjust
    # We impose: w(r_min) ~ r_min/A AND w'(r_min) ~ 1/A
    # But for BVP we can only impose 2 conditions total (one per boundary)
    # So: w(r_min) = r_min/A, w(r_max) = 1

    return np.array([
        ya[0] - ya[0],  # This is always 0 — need proper BC
        yb[0] - 1.0     # w(r_max) = 1
    ])

# Actually, for solve_bvp we need exactly 2 BCs (one per unknown).
# The inner BC should be w(r_min) = r_min/A or w'(r_min) = 1/A.
# The outer BC: w(r_max) = 1.


def solve_soliton_flat(beta=1.0, r_min=0.01, r_max=100.0, n_mesh=500,
                       A_guess=None):
    """
    Solve the flat BH soliton BVP.

    w'' = 4(w')^2/w - (2/r)w' + beta*(1 - 1/w)

    BC: w(r_min) = r_min/A, w(r_max) = 1
    where A = sqrt(2/beta)
    """
    if A_guess is None:
        A_guess = np.sqrt(2.0 / beta)

    A = A_guess

    # Mesh
    r_mesh = np.linspace(r_min, r_max, n_mesh)

    # Initial guess: smooth interpolation from w ~ r/A to w = 1
    # w(r) ~ r/A for small r, transitions to 1 for large r
    # Simple guess: w = 1 - (1 - r/A) * exp(-(r-r_min)/L)
    # Or: w = tanh(r/L) where L adjusts transition
    # Or: w = 1 - exp(-r/r_s) where r_s ~ A (Schwarzschild-like)

    # Try: w = r/(r + A) which gives w(0)=0, w(inf)=1, w'(0)=1/A
    w_guess = r_mesh / (r_mesh + A)
    dw_guess = A / (r_mesh + A)**2

    y_guess = np.vstack([w_guess, dw_guess])

    def ode_func(r, y):
        return ode_w(r, y, beta)

    def bc_func(ya, yb):
        # Inner: w(r_min) = r_min / A (from power-law asymptotics)
        # Outer: w(r_max) = 1
        return np.array([
            ya[0] - r_min / A,
            yb[0] - 1.0
        ])

    print(f"  Solving BVP: beta={beta}, A={A:.4f}, r in [{r_min}, {r_max}], N={n_mesh}")
    print(f"  Inner BC: w({r_min}) = {r_min/A:.6f}")
    print(f"  Outer BC: w({r_max}) = 1.0")

    sol = solve_bvp(ode_func, bc_func, r_mesh, y_guess,
                    tol=1e-6, max_nodes=50000, verbose=0)

    if sol.success:
        print(f"  BVP solved successfully! ({sol.x.size} nodes)")
    else:
        print(f"  BVP failed: {sol.message}")

    return sol


def solve_soliton_sc(beta=1.0, r_min=0.01, r_max=100.0, n_mesh=500,
                     A_guess=None):
    """
    Solve the self-consistent BH soliton BVP.
    """
    if A_guess is None:
        A_guess = np.sqrt(2.0 / beta)

    A = A_guess
    r_mesh = np.linspace(r_min, r_max, n_mesh)

    w_guess = r_mesh / (r_mesh + A)
    dw_guess = A / (r_mesh + A)**2
    y_guess = np.vstack([w_guess, dw_guess])

    def ode_func(r, y):
        return ode_w_selfconsistent(r, y, beta)

    def bc_func(ya, yb):
        return np.array([
            ya[0] - r_min / A,
            yb[0] - 1.0
        ])

    print(f"  Solving SC BVP: beta={beta}, A={A:.4f}")

    sol = solve_bvp(ode_func, bc_func, r_mesh, y_guess,
                    tol=1e-6, max_nodes=50000, verbose=0)

    if sol.success:
        print(f"  SC BVP solved successfully! ({sol.x.size} nodes)")
    else:
        print(f"  SC BVP failed: {sol.message}")

    return sol


# ============================================================
# 3. Post-processing: photon sphere and shadow
# ============================================================

def analyze_solution(sol, beta, label=""):
    """
    Analyze a BVP solution for w(r).

    Convert back to psi = 1/w and compute:
    - Profile psi(r)
    - h(r) = r * psi(r)
    - Photon sphere (minimum of h)
    - Shadow size
    """
    if not sol.success:
        print(f"  [{label}] Solution failed, skipping analysis")
        return None

    r = sol.x
    w = sol.y[0]
    dw = sol.y[1]

    # Convert to psi
    w_safe = np.maximum(w, 1e-30)
    psi = 1.0 / w_safe
    dpsi = -dw / w_safe**2

    print(f"\n  [{label}] SOLUTION ANALYSIS")
    print(f"    r range: [{r[0]:.4f}, {r[-1]:.4f}], N = {len(r)}")
    print(f"    w range: [{w.min():.6e}, {w.max():.6f}]")
    print(f"    psi range: [{psi.min():.6f}, {psi.max():.6e}]")

    # h(r) = r * psi
    h = r * psi
    print(f"    h(r) = r*psi range: [{h.min():.6f}, {h.max():.6f}]")

    # Photon sphere: minimum of h(r)
    idx_min = np.argmin(h)
    h_min = h[idx_min]
    r_ph = r[idx_min]
    psi_ph = psi[idx_min]

    print(f"    h minimum: h = {h_min:.6f} at r = {r_ph:.6f}")
    print(f"    psi(r_ph) = {psi_ph:.6f}")

    # Is this a true minimum (not at boundary)?
    if idx_min > 5 and idx_min < len(h) - 5:
        print(f"    -> Interior minimum: PHOTON SPHERE EXISTS!")

        # Areal radius
        R_ph = r_ph * np.sqrt(psi_ph)
        b_crit = h_min  # = r_ph * psi_ph... wait

        # Actually b_crit in isotropic coords:
        # For ds^2 = -(1/psi)dt^2 + psi(dr^2 + r^2 dOmega^2)
        # b^2 = C/A = psi*r^2 / (1/psi) = psi^2 * r^2
        # b = r * psi
        # Wait no: C = g_{theta theta} = psi * r^2, A = |g_tt| = 1/psi
        # b = sqrt(C/A) = sqrt(psi*r^2*psi) = r*psi

        b_crit = r_ph * psi_ph
        print(f"    b_crit = r_ph * psi_ph = {b_crit:.6f}")

        # Estimate r_s from the Newtonian tail
        # At large r: psi ~ 1/(1 - r_s/r) ~ 1 + r_s/r
        # Wait, w = 1/psi ~ 1 - r_s/r, so psi ~ 1/(1-r_s/r) ~ 1 + r_s/r

        # From w: w ~ 1 - r_s/r => r_s ~ r*(1 - w)
        r_tail = r[-100:]
        w_tail = w[-100:]
        r_s_values = r_tail * (1.0 - w_tail)
        r_s_est = np.median(r_s_values)
        print(f"    Estimated r_s (from tail): {r_s_est:.6f}")

        # GR comparison
        b_gr = (3.0 * np.sqrt(3.0) / 2.0) * r_s_est
        ratio = b_crit / b_gr
        pct = (1.0 - ratio) * 100

        print(f"    b_crit_GR = {b_gr:.6f}")
        print(f"    Shadow ratio = {ratio:.6f}")
        sign = "smaller" if ratio < 1 else "larger"
        print(f"    Shadow: {abs(pct):.2f}% {sign} than Schwarzschild")

        return {
            'r': r, 'psi': psi, 'w': w,
            'r_ph': r_ph, 'psi_ph': psi_ph,
            'b_crit': b_crit, 'r_s': r_s_est,
            'shadow_ratio': ratio,
            'h': h, 'h_min': h_min,
        }
    else:
        print(f"    -> Boundary minimum: no true photon sphere")

        # Still estimate r_s
        r_tail = r[-100:]
        w_tail = w[-100:]
        r_s_values = r_tail * (1.0 - w_tail)
        r_s_est = np.median(r_s_values)
        print(f"    Estimated r_s (from tail): {r_s_est:.6f}")

        return {
            'r': r, 'psi': psi, 'w': w,
            'r_s': r_s_est, 'h': h, 'h_min': h_min,
            'r_ph': None, 'shadow_ratio': None,
        }


# ============================================================
# 4. Profile sampling
# ============================================================

def print_profile(sol, label=""):
    """Print key points of the profile."""
    if not sol.success:
        return

    r = sol.x
    w = sol.y[0]
    psi = 1.0 / np.maximum(w, 1e-30)
    h = r * psi

    print(f"\n  [{label}] Profile sample:")
    r_targets = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0,
                 10.0, 20.0, 50.0, 100.0]
    for rt in r_targets:
        if rt < r[0] or rt > r[-1]:
            continue
        idx = np.argmin(np.abs(r - rt))
        print(f"    r={r[idx]:8.4f}:  w={w[idx]:.6f}  psi={psi[idx]:12.6f}  "
              f"h={h[idx]:10.6f}")


# ============================================================
# 5. Parameter scan
# ============================================================

def scan_beta(betas=None, r_min=0.01, r_max=100.0):
    """Scan over beta values."""
    if betas is None:
        betas = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    print("\n" + "=" * 72)
    print("PARAMETER SCAN (beta)")
    print("=" * 72)

    for beta in betas:
        A = np.sqrt(2.0 / beta)
        print(f"\n  beta = {beta}, A = {A:.4f}")

        sol = solve_soliton_flat(beta, r_min, r_max, n_mesh=500)
        if sol.success:
            result = analyze_solution(sol, beta, f"beta={beta}")
            if result and result.get('shadow_ratio'):
                print(f"    => Shadow ratio = {result['shadow_ratio']:.4f}")
        else:
            print(f"    FAILED")


# ============================================================
# Main
# ============================================================

def main():
    sep = "=" * 72
    print(sep)
    print("TGP BLACK HOLE SOLITON -- BVP SOLVER")
    print(sep)

    beta = 1.0
    A = np.sqrt(2.0 / beta)
    r_min = 0.01
    r_max = 100.0

    # --- A. Flat equation ---
    print("\n" + "-" * 72)
    print("A. FLAT EQUATION")
    print("-" * 72)

    sol_flat = solve_soliton_flat(beta, r_min, r_max, n_mesh=500)
    if sol_flat.success:
        result_flat = analyze_solution(sol_flat, beta, "FLAT")
        print_profile(sol_flat, "FLAT")

    # --- B. Self-consistent equation ---
    print("\n" + "-" * 72)
    print("B. SELF-CONSISTENT EQUATION")
    print("-" * 72)

    sol_sc = solve_soliton_sc(beta, r_min, r_max, n_mesh=500)
    if sol_sc.success:
        result_sc = analyze_solution(sol_sc, beta, "SELF-CONSISTENT")
        print_profile(sol_sc, "SELF-CONSISTENT")

    # --- C. Compare ---
    if sol_flat.success and sol_sc.success:
        print("\n" + "-" * 72)
        print("C. COMPARISON")
        print("-" * 72)

        r_common = np.linspace(r_min, r_max, 5000)
        w_flat = sol_flat.sol(r_common)[0]
        w_sc = sol_sc.sol(r_common)[0]

        psi_flat = 1.0 / np.maximum(w_flat, 1e-30)
        psi_sc = 1.0 / np.maximum(w_sc, 1e-30)

        delta = np.abs(psi_sc - psi_flat) / np.maximum(psi_flat, 1e-10)
        idx_max = np.argmax(delta)

        print(f"  Max relative deviation: {delta[idx_max]:.4e} at r = {r_common[idx_max]:.4f}")
        print(f"  psi_flat = {psi_flat[idx_max]:.6f}, psi_SC = {psi_sc[idx_max]:.6f}")

    # --- D. Parameter scan ---
    if "--scan" in sys.argv:
        scan_beta()

    # --- E. Vary inner BC (different A) ---
    print("\n" + "-" * 72)
    print("D. SENSITIVITY TO INNER BC (varying A)")
    print("-" * 72)

    A_nom = np.sqrt(2.0 / beta)
    A_values = [0.5*A_nom, 0.8*A_nom, A_nom, 1.2*A_nom, 1.5*A_nom, 2.0*A_nom]

    for A_try in A_values:
        sol_try = solve_soliton_flat(beta, r_min, r_max, n_mesh=500, A_guess=A_try)
        if sol_try.success:
            r = sol_try.x
            w = sol_try.y[0]
            psi = 1.0 / np.maximum(w, 1e-30)
            h = r * psi
            idx_hmin = np.argmin(h)

            # Estimate r_s
            r_tail = r[-50:]
            w_tail = w[-50:]
            r_s_est = np.median(r_tail * (1.0 - w_tail))

            ps_info = ""
            if idx_hmin > 5 and idx_hmin < len(h) - 5:
                b_gr = (3*np.sqrt(3)/2)*max(r_s_est, 0.01)
                ratio = (r[idx_hmin]*psi[idx_hmin]) / b_gr
                ps_info = f", shadow ratio={ratio:.4f}"

            print(f"  A={A_try:.4f}: h_min={h[idx_hmin]:.4f} at r={r[idx_hmin]:.4f}"
                  f", r_s~{r_s_est:.4f}{ps_info}")

    print("\nDone.")


if __name__ == "__main__":
    if sys.stdout.encoding != 'utf-8':
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    main()
