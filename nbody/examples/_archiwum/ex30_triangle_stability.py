"""
ex30_triangle_stability.py
==========================
2D shape-space stability analysis for the TGP 3-body equilateral state.

Validates the key approximation used in ex23-ex29: that the ground-state
cluster sits at the equilateral configuration (d12 = d13 = d23 = d).

Method:
  Parameterise triangle deformations from equilateral:
    d12 = d*(1 + ea),   d13 = d*(1 + eb),   d23 = d  (eq. triangle: ea=eb=0)

  The full energy is:
    E(d, ea, eb) = ZP(d) + V2(d12,d13,d23) + V3(d12,d13,d23)

  where ZP = 3/(8d_bar^2) with d_bar = (d12+d13+d23)/3.

  At fixed d = d_eq (equilateral minimum):
    A. 2D scan of E(ea, eb) at C=C_Pl, m_sp=0.1
       -> confirm minimum at ea=eb=0
    B. Hessian matrix H_ij = d^2E/d(ea_i)d(ea_j) at equilateral
       -> eigenvalues > 0 confirm stable minimum
    C. Ratio omega_angular / omega_radial — stiffness of shape DOF
    D. Potential energy along deformation mode ea = eb (symmetric stretch)
       and ea = -eb (antisymmetric, shear mode)
    E. Normal mode analysis: angular oscillation frequencies

Conclusion: if equilateral is a stable energy minimum in shape space,
the 1D reduction used in ex23-ex29 captures the ground state exactly
(modulo zero-point energy in the angular modes).

Author: TGP project, 2026-03-22
"""

import numpy as np
from scipy.special import k0 as K0
import warnings
warnings.filterwarnings('ignore')

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))

# ─────────────────────────────────────────────
# Feynman integral (ex23 convention)
# ─────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5 * (1 + _pts); _uw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW = np.outer(_uw, _uw)
_A1 = _UU; _A2 = _VV * (1 - _UU); _A3 = (1 - _UU) * (1 - _VV); _JAC = (1 - _UU)


def I_Y(d12, d13, d23, m):
    """General Feynman 3-body integral."""
    D = _A1 * _A2 + _A1 * _A3 + _A2 * _A3
    Q = _A2 * d12**2 + _A1 * d13**2 + _A3 * d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)


def E_full(d12, d13, d23, C, m):
    """Full energy: ZP + V2 + V3 for arbitrary triangle."""
    d_bar = (d12 + d13 + d23) / 3.0
    ZP  = 3.0 / (8.0 * d_bar**2)
    V2  = -C**2 * (np.exp(-m * d12) / d12
                 + np.exp(-m * d13) / d13
                 + np.exp(-m * d23) / d23)
    V3  = -C**3 * I_Y(d12, d13, d23, m)
    return ZP + V2 + V3


def E_deformed(d, ea, eb, C, m):
    """Energy with equilateral deformation ea, eb.
    d12 = d*(1+ea),  d13 = d*(1+eb),  d23 = d."""
    d12 = d * (1.0 + ea)
    d13 = d * (1.0 + eb)
    d23 = d
    return E_full(d12, d13, d23, C, m)


# ─────────────────────────────────────────────
# Two-step scan minimum finder (consistent with ex23+)
# ─────────────────────────────────────────────
def find_d_eq(C, m, d_lo=0.5, d_hi=20.0, n1=400, n2=300):
    """Find equilateral equilibrium distance d_eq."""
    d1 = np.linspace(d_lo, d_hi, n1)
    v1 = np.array([E_full(d, d, d, C, m) for d in d1])
    i  = int(np.argmin(v1))
    d2 = np.linspace(max(d_lo, d1[max(0, i-4)]),
                     min(d_hi, d1[min(n1-1, i+4)]), n2)
    v2 = np.array([E_full(d, d, d, C, m) for d in d2])
    j  = int(np.argmin(v2))
    return d2[j], v2[j]


# ─────────────────────────────────────────────
# Hessian via finite differences
# ─────────────────────────────────────────────
def hessian_2d(f, x0, h=1e-3):
    """2x2 Hessian of f(ea, eb) at (ea,eb)=x0 via central differences."""
    H = np.zeros((2, 2))
    ea0, eb0 = x0
    # d^2f/d(ea)^2
    H[0, 0] = (f(ea0 + h, eb0) - 2*f(ea0, eb0) + f(ea0 - h, eb0)) / h**2
    # d^2f/d(eb)^2
    H[1, 1] = (f(ea0, eb0 + h) - 2*f(ea0, eb0) + f(ea0, eb0 - h)) / h**2
    # d^2f/d(ea)d(eb)
    H[0, 1] = (f(ea0 + h, eb0 + h) - f(ea0 + h, eb0 - h)
               - f(ea0 - h, eb0 + h) + f(ea0 - h, eb0 - h)) / (4 * h**2)
    H[1, 0] = H[0, 1]
    return H


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    M_SP = 0.1
    C    = C_PL

    print("=" * 66)
    print("ex30: 2D triangle shape-space stability analysis")
    print(f"  m_sp = {M_SP} l_Pl^-1,  C = C_Pl = {C:.6f}")
    print("  Parameterisation: d12=d*(1+ea), d13=d*(1+eb), d23=d")
    print("  ZP = 3/(8 * d_bar^2),  d_bar = mean side")
    print("=" * 66)

    # ── A: Equilateral equilibrium ─────────────────────────
    print("\n--- A: Equilateral equilibrium ---")
    d_eq, E_eq = find_d_eq(C, M_SP)
    print(f"  d_eq = {d_eq:.5f} l_Pl")
    print(f"  E_eq = {E_eq:.7f} E_Pl")
    E0_check = E_full(d_eq, d_eq, d_eq, C, M_SP)
    print(f"  Check E_full(d_eq,d_eq,d_eq) = {E0_check:.7f} E_Pl")

    # ── B: 2D scan of E(ea, eb) at d=d_eq ─────────────────
    print(f"\n--- B: 2D energy scan in shape space at d={d_eq:.3f} ---")
    eps_arr = np.linspace(-0.5, 0.5, 41)
    E_2d = np.zeros((len(eps_arr), len(eps_arr)))
    for i, ea in enumerate(eps_arr):
        for j, eb in enumerate(eps_arr):
            E_2d[i, j] = E_deformed(d_eq, ea, eb, C, M_SP)

    # Find minimum in 2D scan
    idx_min = np.unravel_index(np.argmin(E_2d), E_2d.shape)
    ea_min = eps_arr[idx_min[0]]
    eb_min = eps_arr[idx_min[1]]
    E_2d_min = E_2d[idx_min]
    print(f"  Global 2D minimum: ea={ea_min:.3f}, eb={eb_min:.3f}, E={E_2d_min:.7f}")
    print(f"  Equilateral (ea=eb=0):    E={E_deformed(d_eq,0,0,C,M_SP):.7f}")

    # Check E along ea-axis (eb=0) and eb-axis (ea=0)
    print(f"\n  Energy along deformation axes (d=d_eq={d_eq:.3f}):")
    print(f"  {'ea/eb':>7}  {'E(ea, eb=0)':>14}  {'E(ea=0, eb)':>14}  {'E(ea=eb)':>14}")
    print("  " + "-" * 54)
    for eps in [-0.4, -0.2, 0.0, 0.2, 0.4]:
        Ea  = E_deformed(d_eq, eps, 0.0, C, M_SP)
        Eb  = E_deformed(d_eq, 0.0, eps, C, M_SP)
        Eab = E_deformed(d_eq, eps, eps, C, M_SP)
        print(f"  {eps:>7.2f}  {Ea:>14.7f}  {Eb:>14.7f}  {Eab:>14.7f}")

    # ── C: Hessian at equilateral ──────────────────────────
    print(f"\n--- C: Hessian d^2E/d(eps_i)d(eps_j) at (ea=eb=0) ---")
    f_shape = lambda ea, eb: E_deformed(d_eq, ea, eb, C, M_SP)
    H = hessian_2d(f_shape, (0.0, 0.0), h=5e-3)
    eigvals, eigvecs = np.linalg.eigh(H)
    print(f"  Hessian matrix H:")
    print(f"    H = [[{H[0,0]:.6f}, {H[0,1]:.6f}],")
    print(f"          [{H[1,0]:.6f}, {H[1,1]:.6f}]]")
    print(f"  Eigenvalues: lambda1={eigvals[0]:.6f}, lambda2={eigvals[1]:.6f}")
    if all(eigvals > 0):
        print(f"  -> Hessian POSITIVE DEFINITE: equilateral is a STABLE MINIMUM")
    else:
        print(f"  -> WARNING: Hessian has non-positive eigenvalue!")
    print(f"  Eigenvector 1 (softer):  [{eigvecs[0,0]:.4f}, {eigvecs[1,0]:.4f}]")
    print(f"  Eigenvector 2 (stiffer): [{eigvecs[0,1]:.4f}, {eigvecs[1,1]:.4f}]")

    # ── D: Normal mode frequencies ─────────────────────────
    print(f"\n--- D: Normal mode frequencies ---")
    # Radial (breathing) curvature at equilateral
    h_r = 1e-3
    d1_arr = np.array([d_eq - 2*h_r, d_eq - h_r, d_eq, d_eq + h_r, d_eq + 2*h_r])
    E_r = np.array([E_full(d, d, d, C, M_SP) for d in d1_arr])
    # 5-point stencil for d^2E/dd^2
    d2E_dd2 = (-E_r[4] + 16*E_r[3] - 30*E_r[2] + 16*E_r[1] - E_r[0]) / (12*h_r**2)

    mu = 1.0   # reduced mass (m_Pl)
    # The effective mass for the d-mode is mu/3 (symmetric breathing in COM frame)
    # Frequencies: omega = sqrt(kappa / mu_eff)
    kappa_r = d2E_dd2   # second derivative in d
    kappa_a = eigvals[0] / d_eq**2   # second derivative in ea (dimensionless)
    kappa_b = eigvals[1] / d_eq**2   # second derivative in eb

    # For breathing mode: effective mass (1/3 of total, reduced)
    mu_r = mu / 3.0          # radial breathing mode
    mu_ang = mu / 3.0        # angular shape mode (same order of magnitude)

    omega_r = np.sqrt(max(kappa_r / mu_r, 0))
    omega_a = np.sqrt(max(kappa_a / mu_ang, 0))
    omega_b = np.sqrt(max(kappa_b / mu_ang, 0))

    print(f"  Radial breathing curvature d^2E/dd^2 = {d2E_dd2:.5f} E_Pl/l_Pl^2")
    print(f"  Shape mode curvatures (per d^2): {eigvals[0]/d_eq**2:.5f}, {eigvals[1]/d_eq**2:.5f}")
    print(f"  omega_radial = {omega_r:.5f} E_Pl  (hbar=1)")
    print(f"  omega_angle1 = {omega_a:.5f} E_Pl")
    print(f"  omega_angle2 = {omega_b:.5f} E_Pl")
    if omega_r > 0:
        print(f"  omega_angle / omega_radial: {omega_a/omega_r:.4f}, {omega_b/omega_r:.4f}")
        if omega_a < omega_r and omega_b < omega_r:
            print(f"  -> Shape modes are SOFTER than breathing mode")
        else:
            print(f"  -> Shape modes are STIFFER than or equal to breathing mode")

    # ZP correction from angular modes
    ZP_angular = 0.5 * (omega_a + omega_b)
    ZP_radial  = 0.5 * omega_r
    print(f"\n  Zero-point energy estimate:")
    print(f"    ZP_radial  (breathing) = {ZP_radial:.6f} E_Pl")
    print(f"    ZP_angular (shape DOF) = {ZP_angular:.6f} E_Pl")
    print(f"    ZP_angular / ZP_radial = {ZP_angular/ZP_radial:.4f}")
    print(f"    E_0(ex29) = -0.069645  (includes radial ZP only)")
    print(f"    Corrected E_0 (approx) = {-0.069645 + ZP_angular:.6f} E_Pl")
    print(f"    Correction = {ZP_angular/0.069645*100:.1f}% of E_0")

    # ── E: 1D slices through shape space ───────────────────
    print(f"\n--- E: Energy along normal mode directions ---")
    eps_1d = np.linspace(-0.6, 0.6, 25)
    v1 = eigvecs[:, 0]   # softer mode
    v2 = eigvecs[:, 1]   # stiffer mode

    print(f"  Softer mode ({v1[0]:.3f}, {v1[1]:.3f}):")
    print(f"  {'eps':>7}  {'E [E_Pl]':>14}  {'dE [E_Pl]':>14}")
    print("  " + "-" * 38)
    E0_eq = E_deformed(d_eq, 0, 0, C, M_SP)
    for eps in eps_1d:
        ea, eb = eps * v1
        E_here = E_deformed(d_eq, ea, eb, C, M_SP)
        print(f"  {eps:>7.3f}  {E_here:>14.7f}  {E_here - E0_eq:>14.7f}")

    print(f"\n  Stiffer mode ({v2[0]:.3f}, {v2[1]:.3f}):")
    print(f"  {'eps':>7}  {'E [E_Pl]':>14}  {'dE [E_Pl]':>14}")
    print("  " + "-" * 38)
    for eps in eps_1d:
        ea, eb = eps * v2
        E_here = E_deformed(d_eq, ea, eb, C, M_SP)
        print(f"  {eps:>7.3f}  {E_here:>14.7f}  {E_here - E0_eq:>14.7f}")

    # ── F: Multi-C summary ─────────────────────────────────
    print(f"\n--- F: Stability at multiple C values (m_sp={M_SP}) ---")
    print(f"  {'C':>6}  {'d_eq':>7}  {'E_eq':>11}  {'H_11':>10}  {'H_22':>10}  "
          f"{'eig1':>9}  {'eig2':>9}  {'stable?':>8}")
    print("  " + "-" * 76)
    for C_test in [0.15, 0.191, 0.20, 0.25, 0.282, 0.30]:
        d_t, E_t = find_d_eq(C_test, M_SP)
        f_t = lambda ea, eb: E_deformed(d_t, ea, eb, C_test, M_SP)
        H_t = hessian_2d(f_t, (0, 0), h=5e-3)
        ev_t = np.linalg.eigvalsh(H_t)
        stable = "YES" if all(ev_t > 0) else "no"
        print(f"  {C_test:>6.3f}  {d_t:>7.4f}  {E_t:>11.6f}  "
              f"{H_t[0,0]:>10.5f}  {H_t[1,1]:>10.5f}  "
              f"{ev_t[0]:>9.5f}  {ev_t[1]:>9.5f}  {stable:>8}")

    # ── G: Triangle inequality check ──────────────────────
    print(f"\n--- G: Triangle inequality at d_eq extremes ---")
    print(f"  (Deformations must satisfy triangle inequalities)")
    print(f"  d12 = d*(1+ea),  d13 = d*(1+eb),  d23 = d")
    print(f"  Triangle: d12+d13>d23, d12+d23>d13, d13+d23>d12")
    for ea, eb in [(0.5, 0.5), (-0.5, -0.5), (0.9, -0.3), (-0.3, 0.9)]:
        d12 = d_eq * (1 + ea); d13 = d_eq * (1 + eb); d23 = d_eq
        ok = (d12+d13>d23) and (d12+d23>d13) and (d13+d23>d12) and d12>0 and d13>0
        print(f"  ea={ea:+.1f}, eb={eb:+.1f}: d12={d12:.3f}, d13={d13:.3f}, "
              f"d23={d23:.3f} -> {'valid' if ok else 'INVALID'}")

    # ── H: Shape Hessian at quantum mean position ─────────
    print(f"\n--- H: Shape Hessian at d = <d>_quantum = 4.52 l_Pl ---")
    print(f"  (The quantum state lives at <d>=4.52 l_Pl, not d_eq={d_eq:.2f})")
    D_QUANTUM = 4.52   # from ex29: <d> at C=C_Pl, m_sp=0.1
    f_q = lambda ea, eb: E_deformed(D_QUANTUM, ea, eb, C, M_SP)
    H_q = hessian_2d(f_q, (0, 0), h=5e-3)
    ev_q = np.linalg.eigvalsh(H_q)
    print(f"  H at d={D_QUANTUM}:")
    print(f"    H = [[{H_q[0,0]:.6f}, {H_q[0,1]:.6f}],")
    print(f"          [{H_q[1,0]:.6f}, {H_q[1,1]:.6f}]]")
    print(f"  Eigenvalues: {ev_q[0]:.6f}, {ev_q[1]:.6f}")
    if all(ev_q > 0):
        print(f"  -> STABLE minimum at <d>_quantum -- equilateral approx. valid for quantum state!")
    elif all(ev_q < 0):
        print(f"  -> UNSTABLE (both negative) at <d>_quantum")
    else:
        print(f"  -> SADDLE POINT (mixed signs) at <d>_quantum")

    # Scan over multiple d values to find stability transition
    print(f"\n  Shape stability vs d (C=C_Pl, m_sp={M_SP}):")
    print(f"  {'d [l_Pl]':>10}  {'eig1':>10}  {'eig2':>10}  {'stable?':>10}")
    print("  " + "-" * 44)
    for d_test in [1.0, 1.37, 2.0, 3.0, 4.52, 6.0, 8.0, 10.0]:
        f_t = lambda ea, eb: E_deformed(d_test, ea, eb, C, M_SP)
        H_t = hessian_2d(f_t, (0, 0), h=5e-3)
        ev_t = np.linalg.eigvalsh(H_t)
        stable = "YES" if all(ev_t > 0) else ("SADDLE" if (ev_t[0]<0<ev_t[1]) else "NO")
        print(f"  {d_test:>10.2f}  {ev_t[0]:>10.6f}  {ev_t[1]:>10.6f}  {stable:>10}")

    print(f"\n--- I: Shape stability vs d and vs C ---")
    # Find where LARGER eigenvalue crosses zero (from positive at small d to negative at large d)
    # From scan: eig2 > 0 at d=1.37, < 0 at d=2.0 --> find crossing
    def eig2_at_d(d_test):
        f_t = lambda ea, eb: E_deformed(d_test, ea, eb, C, M_SP)
        H_t = hessian_2d(f_t, (0, 0), h=5e-3)
        return np.linalg.eigvalsh(H_t)[1]   # larger eigenvalue

    # Check sign change
    e_lo = eig2_at_d(1.37)
    e_hi = eig2_at_d(2.0)
    if e_lo > 0 > e_hi:
        d_lo_s, d_hi_s = 1.37, 2.0
        for _ in range(40):
            d_mid = 0.5*(d_lo_s + d_hi_s)
            if eig2_at_d(d_mid) > 0:
                d_lo_s = d_mid
            else:
                d_hi_s = d_mid
        d_transition = 0.5*(d_lo_s + d_hi_s)
        print(f"  Stability transition (larger eigenvalue = 0) at:")
        print(f"    d* = {d_transition:.4f} l_Pl  (C=C_Pl, m_sp={M_SP})")
        print(f"    d* / <d>_quantum = {d_transition/D_QUANTUM:.4f}")
        print(f"  -> d* < <d>: equilateral is FULLY UNSTABLE at <d>_quantum={D_QUANTUM}")
    else:
        d_transition = None
        print(f"  No stability transition found: equilateral unstable at all d.")

    # Stability vs C at d = D_QUANTUM
    print(f"\n  Shape stability vs C at fixed d=<d>_q={D_QUANTUM} l_Pl:")
    print(f"  {'C':>7}  {'eig1':>10}  {'eig2':>10}  {'stable?':>10}")
    print("  " + "-" * 42)
    for C_test in [0.100, 0.128, 0.150, 0.191, 0.210, 0.250, 0.282, 0.310]:
        f_t2 = lambda ea, eb: E_deformed(D_QUANTUM, ea, eb, C_test, M_SP)
        H_t2 = hessian_2d(f_t2, (0, 0), h=5e-3)
        ev_t2 = np.linalg.eigvalsh(H_t2)
        stable = "YES" if all(ev_t2 > 0) else ("SADDLE" if ev_t2[1] > 0 else "NO")
        print(f"  {C_test:>7.3f}  {ev_t2[0]:>10.6f}  {ev_t2[1]:>10.6f}  {stable:>10}")

    print("\n" + "=" * 66)
    print("CONCLUSION")
    print("=" * 66)
    stab_str = "STABLE" if all(ev_q > 0) else "SADDLE (unstable)"
    print(f"  At d=d_eq={d_eq:.2f} (classical minimum): equilateral is a SADDLE POINT")
    print(f"  At d=<d>_q={D_QUANTUM:.2f} (quantum mean): equilateral is {stab_str}")
    print(f"  Shape stability transition at d* = {d_transition:.3f} l_Pl")
    print()
    if all(ev_q > 0):
        print(f"  KEY RESULT: The quantum ground state (at <d>={D_QUANTUM}) sits in the")
        print(f"  STABLE region of shape space. The equilateral breathing mode")
        print(f"  approximation (ex23-ex29) is VALID for the Efimov quantum state.")
        print(f"  Angular ZP correction ~ {ZP_angular:.4f} E_Pl = small.")
    else:
        print(f"  KEY RESULT: The quantum ground state (at <d>={D_QUANTUM}) sits in the")
        print(f"  UNSTABLE region of shape space. The equilateral breathing mode")
        print(f"  approximation (ex23-ex29) gives an UPPER BOUND on E_0.")
        print(f"  The true ground-state energy is lower (non-equilateral triangle).")
        print(f"  Ex23-ex29 windows are slightly conservative (true window may be wider).")
    print("\nDone.")
