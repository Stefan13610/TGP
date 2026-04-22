"""
ex25_optimal_geometry.py
========================
Find the optimal 3-body geometry maximising binding energy in TGP.

Convention (consistent with ex22-ex24):
  E_eff(d) = ZP + V2_total + V3_total
  ZP = N/(8*d^2)  where N=3 for triangles  [ex23 convention]
  d  = reference separation (equilateral side; nearest-neighbour otherwise)

  Feynman integral parameterization (ex23 convention):
    alpha1=u, alpha2=v*(1-u), alpha3=(1-u)*(1-v), Jac=(1-u)
    Q = A2*r12^2 + A1*r13^2 + A3*r23^2
    Delta = A1*A2 + A1*A3 + A2*A3

Geometries compared (all parameterised by one scale d):
  1. Equilateral  (r12=r13=r23=d)
  2. Linear       (r12=r23=d, r13=2d)
  3. Isosceles    (r12=r13=d, r23=s*d) with s optimised
  4. Square N=4   (side d, diagonal d*sqrt(2)) — ZP = 4/(8d^2)

Author: TGP project, 2026-03-21
"""

import numpy as np
from scipy.special import k0 as K0
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
# Gauss-Legendre nodes (ex23 convention)
# ─────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_vp = 0.5*(1+_pts); _vw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW = np.outer(_uw, _vw)
_A1 = _UU
_A2 = _VV*(1-_UU)
_A3 = (1-_UU)*(1-_VV)
_JAC = (1-_UU)

C_PL = 0.282094

# ─────────────────────────────────────────────
# Feynman 3-body Yukawa  I_Y(d12, d13, d23; m)
# Matches ex23 exactly (vectorised meshgrid)
# ─────────────────────────────────────────────
def I_Y(d12, d13, d23, m):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12**2 + _A1*d13**2 + _A3*d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u    = np.where(good, m*np.sqrt(Q/D), 1.0)
    val  = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)

def I_Y_equil(d, m):
    return I_Y(d, d, d, m)

# ─────────────────────────────────────────────
# LUT for equilateral  I_Y(d,d,d; m)
# ─────────────────────────────────────────────
def build_V3_lut(m, d_lo=0.2, d_hi=60.0, n=400):
    d_arr  = np.linspace(d_lo, d_hi, n)
    IY_arr = np.array([I_Y_equil(d, m) for d in d_arr])
    lut    = interp1d(d_arr, IY_arr, kind='cubic', bounds_error=False,
                      fill_value=(IY_arr[0], 0.0))
    return lut, d_arr, IY_arr

# ─────────────────────────────────────────────
# Primitives  (ex23 convention)
# ─────────────────────────────────────────────
def V2(d, C, m):         return -3.0*C**2*np.exp(-m*d)/d
def V3_equil(d, C, m, lut): return -C**3*float(lut(d))
def E_ZP(d, N=3):        return N/(8.0*d**2)

# ─────────────────────────────────────────────
# Geometry energies
# ─────────────────────────────────────────────
def E_equilateral(d, C, m, lut):
    """Standard equilateral: r12=r13=r23=d."""
    return E_ZP(d,3) + V2(d,C,m) + V3_equil(d,C,m,lut)

def E_linear(d, C, m):
    """Linear: r12=r23=d, r13=2d.
    V2 = -C^2*(2e^{-md}/d + e^{-2md}/(2d))  [3 pairs]"""
    v2 = -C**2*(2*np.exp(-m*d)/d + np.exp(-m*2*d)/(2*d))
    v3 = -C**3 * I_Y(d, 2*d, d, m)   # (r12, r13, r23) = (d, 2d, d)
    return E_ZP(d,3) + v2 + v3

def E_isosceles(d, s, C, m):
    """Isosceles: r12=r13=d, r23=s*d.
    V2 = -C^2*(2e^{-md}/d + e^{-msd}/(sd))"""
    v2 = -C**2*(2*np.exp(-m*d)/d + np.exp(-m*s*d)/(s*d))
    v3 = -C**3 * I_Y(d, d, s*d, m)   # (r12, r13, r23) = (d, d, s*d)
    return E_ZP(d,3) + v2 + v3

def E_square_4body(d, C, m, lut):
    """Square N=4: side d, diagonal d*sqrt(2).
    6 pairs: 4×(d) + 2×(d√2).  ZP = 4/(8d^2) (4 nearest-neighbour pairs).
    V3: 4 right-triangle triples (d, d, d√2)."""
    r_diag = d * np.sqrt(2.0)
    v2 = -C**2*(4*np.exp(-m*d)/d + 2*np.exp(-m*r_diag)/r_diag)
    IY_right = I_Y(d, d, r_diag, m)   # right-triangle: r12=r13=d, r23=d√2
    v3 = -C**3 * 4.0 * IY_right
    return E_ZP(d,4) + v2 + v3


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    M_SP  = 0.1
    C_VAL = 0.155
    print("=" * 62)
    print("ex25: Optimal 3B geometry in TGP")
    print(f"  m_sp = {M_SP},  C = {C_VAL}  (Efimov-window midpoint)")
    print("  ZP convention: N/(8d^2) [N=3 for triangles, N=4 for square]")
    print("=" * 62)

    print("\nBuilding V3 LUT (equilateral, d in [0.2, 60], n=400)...")
    lut, _, _ = build_V3_lut(M_SP, 0.2, 60.0, 400)
    print("  LUT ready.")

    d_scan  = np.linspace(0.3, 30.0, 600)
    d_coarse = d_scan[::4]

    # ── A: Equilateral ───────────────────────────────────
    print("\n--- A: Equilateral triangle ---")
    E_eq = np.array([E_equilateral(d, C_VAL, M_SP, lut) for d in d_scan])
    idx  = np.argmin(E_eq)
    d_eq_min, E_eq_min = d_scan[idx], E_eq[idx]
    print(f"  E_min = {E_eq_min:.5f} E_Pl   at d_eq = {d_eq_min:.3f} l_Pl")

    # ── B: Linear ────────────────────────────────────────
    print("\n--- B: Linear (r12=r23=d, r13=2d) ---")
    E_lin = np.array([E_linear(d, C_VAL, M_SP) for d in d_coarse])
    idx   = np.argmin(E_lin)
    d_lin_min, E_lin_min = d_coarse[idx], E_lin[idx]
    print(f"  E_min = {E_lin_min:.5f} E_Pl   at d_eq = {d_lin_min:.3f} l_Pl")
    rel_lin = (E_lin_min - E_eq_min)/abs(E_eq_min)*100 if E_eq_min else 0
    print(f"  Relative to equilateral: {rel_lin:+.1f}%")

    # ── C: Isosceles 2D scan ─────────────────────────────
    print("\n--- C: Isosceles 2D scan (d, s=r23/d) ---")
    print("  [Scanning 50x40 = 2000 points...]")
    d_iso = np.linspace(0.5, 20.0, 50)
    s_iso = np.linspace(0.3, 4.0,  40)
    E_best_iso = 0.0; d_best_iso = d_eq_min; s_best_iso = 1.0
    # Also scan along d at each s to find (d_opt, s_opt)
    for d in d_iso:
        for s in s_iso:
            E = E_isosceles(d, s, C_VAL, M_SP)
            if E < E_best_iso:
                E_best_iso, d_best_iso, s_best_iso = E, d, s
    pct_iso = (E_best_iso - E_eq_min)/abs(E_eq_min)*100 if E_eq_min else 0
    print(f"  Global 2D min: E={E_best_iso:.5f} E_Pl  at d={d_best_iso:.3f}, s={s_best_iso:.2f}")
    print(f"    r23 = {s_best_iso*d_best_iso:.3f} l_Pl")
    print(f"  Equilateral (s=1) at same d: E={E_equilateral(d_best_iso,C_VAL,M_SP,lut):.5f}")
    print(f"  Isosceles improvement vs equilateral: {pct_iso:+.1f}%")

    # ── D: Square N=4 ────────────────────────────────────
    print("\n--- D: Square N=4 (side d) ---")
    E_sq = np.array([E_square_4body(d, C_VAL, M_SP, lut) for d in d_coarse])
    idx  = np.argmin(E_sq)
    d_sq_min, E_sq_min = d_coarse[idx], E_sq[idx]
    print(f"  E_min = {E_sq_min:.5f} E_Pl   at d_side = {d_sq_min:.3f} l_Pl")

    # ── E: Summary table ─────────────────────────────────
    print("\n" + "=" * 62)
    print(f"SUMMARY: Binding energy vs geometry  (m_sp={M_SP}, C={C_VAL})")
    print("=" * 62)
    print(f"  {'Geometry':<26} {'E_min [E_Pl]':>13} {'d_ref [l_Pl]':>13} {'Bound?':>7}")
    print("  " + "-" * 60)
    rows = [
        ("Equilateral (N=3)",           E_eq_min,    d_eq_min,    E_eq_min < 0),
        ("Linear (N=3)",                E_lin_min,   d_lin_min,   E_lin_min < 0),
        (f"Isosceles s={s_best_iso:.2f} (N=3)", E_best_iso, d_best_iso, E_best_iso < 0),
        ("Square (N=4)",                E_sq_min,    d_sq_min,    E_sq_min < 0),
    ]
    for name, Emin, deq, bound in rows:
        print(f"  {name:<26} {Emin:>13.5f} {deq:>13.3f} {'YES' if bound else 'NO':>7}")
    print("  " + "-" * 60)
    if E_eq_min < 0:
        print(f"\n  Linear vs equilateral:    {rel_lin:+.1f}% (equilateral stronger)")
        print(f"  Isosceles vs equilateral: {pct_iso:+.1f}%")

    # ── F: 3B-only window — linear vs equilateral ────────
    print("\n--- F: 3B-only window  (equilateral vs linear, m_sp=0.1) ---")
    C_arr = np.linspace(0.08, 0.21, 27)
    print(f"  {'C':>5}  {'E2B(eq)':>10}  {'E3B(eq)':>10}  {'E3B(lin)':>10}  "
          f"{'eq-win?':>8}  {'lin-win?':>9}")
    print("  " + "-" * 62)
    for C in C_arr:
        E_2b     = np.array([E_ZP(d,3) + V2(d,C,M_SP) for d in d_scan])
        E2b_min  = E_2b.min()
        E_3b_eq  = np.array([E_equilateral(d, C, M_SP, lut) for d in d_scan])
        E3b_eq_min = E_3b_eq.min()
        E_3b_lin = np.array([E_linear(d, C, M_SP) for d in d_coarse])
        E3b_lin_min = E_3b_lin.min()
        eq_win   = (E2b_min > 0) and (E3b_eq_min  < 0)
        lin_win  = (E2b_min > 0) and (E3b_lin_min < 0)
        print(f"  {C:>5.3f}  {E2b_min:>10.5f}  {E3b_eq_min:>10.5f}  {E3b_lin_min:>10.5f}  "
              f"{'YES' if eq_win else 'no':>8}  {'YES' if lin_win else 'no':>9}")

    # ── G: Isosceles s-scan at C=0.155 ───────────────────
    print("\n--- G: Isosceles E(s) at d=d_eq_min, C=0.155 ---")
    print(f"  (fixed d={d_eq_min:.3f} l_Pl, vary s)")
    s_arr = np.linspace(0.2, 5.0, 100)
    E_s   = [E_isosceles(d_eq_min, s, C_VAL, M_SP) for s in s_arr]
    idx_s = np.argmin(E_s)
    print(f"  Optimal s = {s_arr[idx_s]:.2f}  =>  E_min = {E_s[idx_s]:.5f} E_Pl")
    print(f"  Equilateral (s=1): E = {E_isosceles(d_eq_min, 1.0, C_VAL, M_SP):.5f}")
    # Print table at selected s values
    print(f"  {'s':>5}  {'r23 [l_Pl]':>11}  {'E [E_Pl]':>10}")
    for s in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]:
        E_here = E_isosceles(d_eq_min, s, C_VAL, M_SP)
        print(f"  {s:>5.2f}  {s*d_eq_min:>11.3f}  {E_here:>10.5f}")

    print("\nDone.")
