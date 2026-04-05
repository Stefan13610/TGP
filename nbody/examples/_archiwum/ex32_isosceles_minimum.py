"""
ex32_isosceles_minimum.py
=========================
Find the true non-equilateral classical energy minimum for the TGP 3-body
system by scanning the full (a, b) space of isosceles triangles:
    d12 = d13 = a,   d23 = b    (S2 symmetry along the axis d12=d13)

By S3 permutation symmetry, the absolute minimum either lies on the equilateral
line a=b (all three sides equal) or on one of the three isosceles lines.
The isosceles submanifold is the natural complement to the equilateral breathing
mode, and directly tests whether ex30's saddle-point conclusion is correct.

Method:
  1. 2D grid scan over (a, b) with triangle inequalities:
       a > 0,  b > 0,  2a > b  (triangle inequality, b < 2a)
  2. Find grid minimum, refine with scipy.optimize.minimize
  3. Compare isosceles E_min with equilateral E_min at same C, m_sp
  4. Compute energy gain ΔE = E_iso - E_eq (should be negative: iso is lower)
  5. Scan vs m_sp to see how window boundaries shift

Parameters scanned:
  m_sp = 0.1 (nominal), then 0.076, 0.130, 0.198 (window boundaries)
  C = C_Pl = 1/(2sqrt(pi)) ~ 0.2821

Author: TGP project, 2026-03-22
"""

import numpy as np
from scipy.special import k0 as K0
from scipy.optimize import minimize
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
    """Full energy: ZP + V2 + V3 for arbitrary triangle.
    ZP = N/(8*d_bar^2) with N=3."""
    d_bar = (d12 + d13 + d23) / 3.0
    ZP  = 3.0 / (8.0 * d_bar**2)
    V2  = -C**2 * (np.exp(-m * d12) / d12
                 + np.exp(-m * d13) / d13
                 + np.exp(-m * d23) / d23)
    V3  = -C**3 * I_Y(d12, d13, d23, m)
    return ZP + V2 + V3


def E_isosceles(a, b, C, m):
    """Energy for isosceles triangle: d12=d13=a, d23=b.
    Triangle inequality: 2a > b and b > 0 and a > 0."""
    if a <= 0 or b <= 0 or b >= 2.0 * a:
        return 1e10
    return E_full(a, a, b, C, m)


def E_equilateral(d, C, m):
    """Energy for equilateral triangle d12=d13=d23=d."""
    return E_full(d, d, d, C, m)


def find_equilateral_min(C, m, d_lo=0.3, d_hi=25.0, n1=500, n2=400):
    """Two-step scan for equilateral minimum."""
    d1 = np.linspace(d_lo, d_hi, n1)
    v1 = np.array([E_equilateral(d, C, m) for d in d1])
    i  = int(np.argmin(v1))
    d2 = np.linspace(max(d_lo, d1[max(0, i-4)]),
                     min(d_hi, d1[min(n1-1, i+4)]), n2)
    v2 = np.array([E_equilateral(d, C, m) for d in d2])
    j  = int(np.argmin(v2))
    return d2[j], v2[j]


def find_isosceles_min(C, m, na=60, nb=60, a_lo=0.3, a_hi=15.0,
                       b_lo=0.3, b_hi=15.0):
    """2D grid scan + scipy refinement for isosceles minimum.
    Returns (a_min, b_min, E_min)."""
    a_arr = np.linspace(a_lo, a_hi, na)
    b_arr = np.linspace(b_lo, b_hi, nb)
    E_grid = np.full((na, nb), 1e10)
    for i, a in enumerate(a_arr):
        for j, b in enumerate(b_arr):
            if b < 2.0 * a and b > 0:
                E_grid[i, j] = E_isosceles(a, b, C, m)

    # Grid minimum
    ij_min = np.unravel_index(np.argmin(E_grid), E_grid.shape)
    a0 = a_arr[ij_min[0]]
    b0 = b_arr[ij_min[1]]

    # Scipy refinement with bounds to prevent singularities
    b_lo_safe = max(b_lo, 0.1)
    a_lo_safe = max(a_lo, 0.1)
    bounds = [(a_lo_safe, a_hi), (b_lo_safe, a_hi)]  # b < 2a enforced in E_isosceles

    def neg_ok(x):
        a, b = x
        return E_isosceles(a, b, C, m)

    res = minimize(neg_ok, [a0, b0], method='L-BFGS-B', bounds=bounds,
                   options={'ftol': 1e-12, 'gtol': 1e-8, 'maxiter': 5000})
    a_min, b_min = res.x
    E_min = res.fun
    # Fallback to grid if minimizer diverged or hit boundary
    if E_min > E_grid[ij_min] or b_min < b_lo_safe + 1e-6:
        a_min, b_min = a0, b0
        E_min = E_grid[ij_min]
    return a_min, b_min, E_min


def run_analysis(C=C_PL, m=0.1):
    print(f"\n{'='*62}")
    print(f"  TGP ex32: Isosceles minimum   C={C:.4f}  m_sp={m:.4f}")
    print(f"{'='*62}")

    # ── 1. Equilateral minimum ──────────────────────────────────
    d_eq, E_eq = find_equilateral_min(C, m)
    print(f"\n  Equilateral minimum:")
    print(f"    d_eq = {d_eq:.4f} l_Pl")
    print(f"    E_eq = {E_eq:.6f} E_Pl")

    # Components at equilateral
    d_bar_eq = d_eq
    ZP_eq  = 3.0 / (8.0 * d_bar_eq**2)
    V2_eq  = -C**2 * 3 * np.exp(-m * d_eq) / d_eq
    V3_eq  = -C**3 * I_Y(d_eq, d_eq, d_eq, m)
    print(f"    ZP   = {ZP_eq:.6f},  V2 = {V2_eq:.6f},  V3 = {V3_eq:.6f}")

    # ── 2. Isosceles minimum (2D scan) ──────────────────────────
    a_min, b_min, E_iso = find_isosceles_min(C, m)

    # Sanity: ensure it passes triangle inequality
    if b_min >= 2 * a_min or a_min <= 0 or b_min <= 0:
        print(f"\n  WARNING: isosceles optimizer left valid domain. Trying wider grid.")
        a_min, b_min, E_iso = find_isosceles_min(C, m, na=80, nb=80,
                                                  a_lo=0.2, a_hi=20.0,
                                                  b_lo=0.2, b_hi=20.0)

    d_bar_iso = (2 * a_min + b_min) / 3.0
    ZP_iso  = 3.0 / (8.0 * d_bar_iso**2)
    V2_iso  = -C**2 * (2 * np.exp(-m * a_min) / a_min
                      + np.exp(-m * b_min) / b_min)
    V3_iso  = -C**3 * I_Y(a_min, a_min, b_min, m)

    print(f"\n  Isosceles minimum:")
    print(f"    a_min = {a_min:.4f} l_Pl  (d12=d13)")
    print(f"    b_min = {b_min:.4f} l_Pl  (d23)")
    print(f"    b/a   = {b_min/a_min:.4f}  (equilateral: b/a=1)")
    print(f"    d_bar = {d_bar_iso:.4f} l_Pl")
    print(f"    E_iso = {E_iso:.6f} E_Pl")
    print(f"    ZP    = {ZP_iso:.6f},  V2 = {V2_iso:.6f},  V3 = {V3_iso:.6f}")

    # ── 3. Energy gain ───────────────────────────────────────────
    dE = E_iso - E_eq
    print(f"\n  Energy gain (iso - eq):")
    direction = 'LOWER than' if dE < 0 else 'HIGHER than'
    print(f"    dE = {dE:.6f} E_Pl  ({direction} equilateral)")
    if abs(E_eq) > 1e-10:
        print(f"    dE/|E_eq| = {dE/abs(E_eq)*100:.2f}%")

    # ── 4. Deformation parameters ────────────────────────────────
    d_ref = (a_min + b_min) / 2.0   # rough mean bond
    print(f"\n  Geometry vs equilateral d={d_eq:.3f}:")
    print(f"    a/d_eq = {a_min/d_eq:.4f}")
    print(f"    b/d_eq = {b_min/d_eq:.4f}")

    return dict(d_eq=d_eq, E_eq=E_eq, a_min=a_min, b_min=b_min,
                E_iso=E_iso, dE=dE)


def scan_vs_m_sp():
    """Scan the isosceles minimum energy gain vs m_sp."""
    print("\n\n" + "="*62)
    print("  ex32: Scan vs m_sp  (C = C_Pl)")
    print("="*62)
    print(f"\n  {'m_sp':>8}  {'lambda':>8}  {'E_eq':>10}  {'E_iso':>10}"
          f"  {'dE':>10}  {'dE%':>8}  {'b/a':>8}")
    print(f"  {'-'*70}")

    # Window boundaries from ex26/ex27 + nominal
    m_vals = [0.050, 0.076, 0.100, 0.130, 0.150, 0.180, 0.198, 0.250]
    results = []
    for m in m_vals:
        d_eq, E_eq = find_equilateral_min(C_PL, m)
        a_min, b_min, E_iso = find_isosceles_min(C_PL, m)
        dE = E_iso - E_eq
        pct = dE / abs(E_eq) * 100 if abs(E_eq) > 1e-10 else float('nan')
        ba = b_min / a_min
        lam = 1.0 / m
        print(f"  {m:8.3f}  {lam:8.2f}  {E_eq:10.5f}  {E_iso:10.5f}"
              f"  {dE:10.5f}  {pct:8.2f}%  {ba:8.4f}")
        results.append((m, lam, E_eq, E_iso, dE, pct, ba))
    return results


def scan_vs_C(m=0.1):
    """Scan the isosceles minimum energy gain vs C at fixed m_sp."""
    print("\n\n" + "="*62)
    print(f"  ex32: Scan vs C  (m_sp = {m})")
    print("="*62)
    print(f"\n  {'C':>8}  {'E_eq':>10}  {'E_iso':>10}  {'dE':>10}"
          f"  {'dE%':>8}  {'b/a':>8}")
    print(f"  {'-'*60}")

    C_vals = [0.100, 0.128, 0.150, 0.191, 0.220, 0.250, 0.282]
    for C in C_vals:
        d_eq, E_eq = find_equilateral_min(C, m)
        if E_eq >= 0:
            print(f"  {C:8.3f}  (unbound)")
            continue
        a_min, b_min, E_iso = find_isosceles_min(C, m)
        dE = E_iso - E_eq
        pct = dE / abs(E_eq) * 100 if abs(E_eq) > 1e-10 else float('nan')
        ba = b_min / a_min
        print(f"  {C:8.3f}  {E_eq:10.5f}  {E_iso:10.5f}  {dE:10.5f}"
              f"  {pct:8.2f}%  {ba:8.4f}")


if __name__ == '__main__':
    # Main analysis at nominal parameters
    result = run_analysis(C=C_PL, m=0.1)

    # Scan vs m_sp
    scan_results = scan_vs_m_sp()

    # Scan vs C at m_sp=0.1
    scan_vs_C(m=0.1)

    print("\n\n" + "="*62)
    print("  SUMMARY")
    print("="*62)
    r = result
    print(f"""
  At C=C_Pl, m_sp=0.1:
    Equilateral minimum:   d = {r['d_eq']:.3f} l_Pl,  E = {r['E_eq']:.5f} E_Pl
    Isosceles minimum:     a = {r['a_min']:.3f}, b = {r['b_min']:.3f} l_Pl,  E = {r['E_iso']:.5f} E_Pl
    Energy gain dE:        {r['dE']:.5f} E_Pl  ({r['dE']/abs(r['E_eq'])*100:.1f}%)
    b/a ratio:             {r['b_min']/r['a_min']:.4f}

  The equilateral is a SADDLE POINT (confirmed by ex30).
  The isosceles minimum lies {'below' if r['dE']<0 else 'above'} the equilateral minimum.
  => ex23-ex29 give {'UPPER' if r['dE']<0 else 'LOWER'} bounds on E0.
""")
