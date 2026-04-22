"""
ex28_polygon_geometry.py
========================
Compare binding energy for regular polygons N=3 (equilateral triangle),
N=4 (square), and N=5 (regular pentagon) in TGP.

For each polygon, parameterised by nearest-neighbour side length d:

  N=3 equilateral:
    Pairs: 3 × d
    Triples (C(3,3)=1): (d,d,d)
    ZP = 3/(8d^2)

  N=4 square:
    Pairs: 4×d + 2×d*sqrt(2)
    Triples (C(4,3)=4): 4 × (d,d,d*sqrt(2))  [right-angle type]
    ZP = 4/(8d^2)

  N=5 regular pentagon:
    Pairs: 5×d + 5×phi*d  (phi = golden ratio = 2*cos(pi/5) approx 1.618)
    Triples (C(5,3)=10):
      5 × type A: (d, d, phi*d)         [consecutive vertices]
      5 × type B: (d, phi*d, phi*d)     [skip-one vertices]
    ZP = 5/(8d^2)

  Convention (consistent with ex23-ex27):
    ZP = N/(8d^2)
    V2 = -C^2 * sum_{pairs} e^{-m*r_ij}/r_ij
    V3 = -C^3 * sum_{triples} I_Y(r_ij, r_ik, r_jk; m)

Goals:
  1. Find E_min for each polygon vs C at fixed m_sp=0.1
  2. Determine 3B-only window (E_2B > 0, E_3B < 0)
  3. Compare V3/V2 ratios and equilibrium separations
  4. Pentagon comparison at C=C_Pl=0.282

Author: TGP project, 2026-03-21
"""

import numpy as np
from scipy.special import k0 as K0
import warnings
warnings.filterwarnings('ignore')

PHI = 2.0 * np.cos(np.pi / 5.0)   # golden ratio: ~1.6180339...
C_PL = 1.0 / (2.0 * np.sqrt(np.pi))  # ~0.28209

# ─────────────────────────────────────────────
# Gauss-Legendre quadrature (ex23 convention)
# ─────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5 * (1 + _pts); _uw = 0.5 * _wts
_vp = 0.5 * (1 + _pts); _vw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW = np.outer(_uw, _vw)
_A1 = _UU
_A2 = _VV * (1 - _UU)
_A3 = (1 - _UU) * (1 - _VV)
_JAC = (1 - _UU)


def I_Y(d12, d13, d23, m):
    """Feynman 3-body Yukawa integral I_Y(d12,d13,d23; m).
    Convention: alpha1=u, alpha2=v(1-u), alpha3=(1-u)(1-v), Jac=(1-u).
    Q = A2*d12^2 + A1*d13^2 + A3*d23^2.
    """
    D = _A1 * _A2 + _A1 * _A3 + _A2 * _A3
    Q = _A2 * d12**2 + _A1 * d13**2 + _A3 * d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)


# ─────────────────────────────────────────────
# Energy components
# ─────────────────────────────────────────────
def E_ZP(d, N):
    """Zero-point kinetic energy. Convention: N/(8d^2)."""
    return N / (8.0 * d**2)


def V2_yukawa(r, C, m):
    """Single Yukawa pair energy: -C^2 * exp(-m*r)/r."""
    return -C**2 * np.exp(-m * r) / r


def V3_feynman(r12, r13, r23, C, m):
    """Single 3-body Feynman term: -C^3 * I_Y(...)."""
    return -C**3 * I_Y(r12, r13, r23, m)


# ─────────────────────────────────────────────
# Per-geometry effective energies
# ─────────────────────────────────────────────
def E_equilateral(d, C, m):
    """N=3 equilateral: all sides d."""
    v2 = 3.0 * V2_yukawa(d, C, m)
    v3 = V3_feynman(d, d, d, C, m)
    return E_ZP(d, 3) + v2 + v3


def E_triangle_2B(d, C, m):
    """N=3 equilateral, 2-body only."""
    return E_ZP(d, 3) + 3.0 * V2_yukawa(d, C, m)


def E_square(d, C, m):
    """N=4 square: 4 sides d, 2 diagonals d*sqrt(2).
    4 triples of type (d, d, d*sqrt(2))."""
    r_diag = d * np.sqrt(2.0)
    v2 = 4.0 * V2_yukawa(d, C, m) + 2.0 * V2_yukawa(r_diag, C, m)
    v3 = 4.0 * V3_feynman(d, d, r_diag, C, m)
    return E_ZP(d, 4) + v2 + v3


def E_square_2B(d, C, m):
    """N=4 square, 2-body only."""
    r_diag = d * np.sqrt(2.0)
    return E_ZP(d, 4) + 4.0 * V2_yukawa(d, C, m) + 2.0 * V2_yukawa(r_diag, C, m)


def E_pentagon(d, C, m,
               _IY_A=None, _IY_B=None):
    """N=5 regular pentagon with nearest-neighbour side d.
    phi = 2*cos(pi/5) ~ 1.618 (golden ratio).
    Pairs:  5×d + 5×phi*d
    Triples: 5 × type A (d,d,phi*d), 5 × type B (d,phi*d,phi*d)
    ZP = 5/(8d^2)
    """
    r_phi = PHI * d
    v2 = 5.0 * V2_yukawa(d, C, m) + 5.0 * V2_yukawa(r_phi, C, m)
    IY_A = I_Y(d, d, r_phi, m) if _IY_A is None else _IY_A(d)
    IY_B = I_Y(d, r_phi, r_phi, m) if _IY_B is None else _IY_B(d)
    v3 = -C**3 * (5.0 * IY_A + 5.0 * IY_B)
    return E_ZP(d, 5) + v2 + v3


def E_pentagon_2B(d, C, m):
    """N=5 pentagon, 2-body only."""
    r_phi = PHI * d
    return E_ZP(d, 5) + 5.0 * V2_yukawa(d, C, m) + 5.0 * V2_yukawa(r_phi, C, m)


# ─────────────────────────────────────────────
# Two-step scan minimum finder
# ─────────────────────────────────────────────
def find_Emin_scan(E_func, d_lo=0.3, d_hi=25.0, n_coarse=600, n_fine=500):
    """Find minimum of E_func(d) via two-step scan.
    Avoids Brent's golden-section failure for shallow minima."""
    d1 = np.linspace(d_lo, d_hi, n_coarse)
    v1 = np.array([E_func(d) for d in d1])
    i = int(np.argmin(v1))
    d2 = np.linspace(max(d_lo, d1[max(0, i - 4)]),
                     min(d_hi, d1[min(len(d1) - 1, i + 4)]), n_fine)
    v2 = np.array([E_func(d) for d in d2])
    j = int(np.argmin(v2))
    return v2[j], d2[j]


# ─────────────────────────────────────────────
# Build cached I_Y LUTs for pentagon
# ─────────────────────────────────────────────
def build_pentagon_luts(m, d_lo=0.2, d_hi=30.0, n=300):
    """Build interpolated LUTs for the two pentagon triple types.
    Returns (lut_A, lut_B) callables: d -> I_Y value."""
    from scipy.interpolate import interp1d
    d_arr = np.linspace(d_lo, d_hi, n)
    IY_A = np.array([I_Y(d, d, PHI * d, m) for d in d_arr])
    IY_B = np.array([I_Y(d, PHI * d, PHI * d, m) for d in d_arr])
    lut_A = interp1d(d_arr, IY_A, kind='cubic', bounds_error=False,
                     fill_value=(IY_A[0], 0.0))
    lut_B = interp1d(d_arr, IY_B, kind='cubic', bounds_error=False,
                     fill_value=(IY_B[0], 0.0))
    return lut_A, lut_B


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    M_SP = 0.1
    D_LO, D_HI = 0.3, 25.0

    print("=" * 68)
    print("ex28: Regular polygon binding energies in TGP")
    print(f"  m_sp = {M_SP} l_Pl^-1  (lambda = {1/M_SP:.1f} l_Pl)")
    print(f"  C_Pl = {C_PL:.6f}")
    print(f"  phi  = {PHI:.6f}  (golden ratio = 2*cos(pi/5))")
    print("=" * 68)

    print(f"\nBuilding pentagon I_Y LUTs (n=300, d in [{D_LO},{D_HI}])...")
    lut_A, lut_B = build_pentagon_luts(M_SP, D_LO, D_HI, 300)
    print("  LUTs ready.")

    # ─────────────────────────────────────────
    # Section A: E_min comparison at C=0.155
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print("A: Binding energy at C=0.155 (equilateral 3B-window midpoint)")
    print("=" * 68)
    C_MID = 0.155

    def pent3B(d): return E_pentagon(d, C_MID, M_SP, lut_A, lut_B)
    def pent2B(d): return E_pentagon_2B(d, C_MID, M_SP)
    def tri3B(d):  return E_equilateral(d, C_MID, M_SP)
    def tri2B(d):  return E_triangle_2B(d, C_MID, M_SP)
    def sq3B(d):   return E_square(d, C_MID, M_SP)
    def sq2B(d):   return E_square_2B(d, C_MID, M_SP)

    geoms = [
        ("Triangle (N=3) 2B",   tri2B),
        ("Triangle (N=3) 2B+3B",tri3B),
        ("Square   (N=4) 2B",   sq2B),
        ("Square   (N=4) 2B+3B",sq3B),
        ("Pentagon (N=5) 2B",   pent2B),
        ("Pentagon (N=5) 2B+3B",pent3B),
    ]

    print(f"\n  {'Geometry':<28} {'E_min [E_Pl]':>14} {'d_ref [l_Pl]':>14} {'Bound?':>7}")
    print("  " + "-" * 64)
    results_mid = {}
    for name, func in geoms:
        Emin, dmin = find_Emin_scan(func, D_LO, D_HI)
        bound = "YES" if Emin < 0 else "no"
        print(f"  {name:<28} {Emin:>14.6f} {dmin:>14.3f} {bound:>7}")
        results_mid[name] = (Emin, dmin)
    print()

    # 3B contribution only
    print("  3-body energy contribution at equilibrium:")
    for N, label, e3B_f, e2B_f in [
        (3, "Triangle", tri3B, tri2B),
        (4, "Square",   sq3B,  sq2B),
        (5, "Pentagon", pent3B, pent2B),
    ]:
        Emin3B, d3B = find_Emin_scan(e3B_f, D_LO, D_HI)
        V2_at = e2B_f(d3B) - E_ZP(d3B, N)
        V3_at = Emin3B - E_ZP(d3B, N) - V2_at
        ratio = V3_at / V2_at if abs(V2_at) > 1e-12 else float('nan')
        print(f"  {label:<9}: d_eq={d3B:.3f}, V2={V2_at:+.6f}, V3={V3_at:+.6f}, "
              f"V3/V2={ratio:.3f}")

    # ─────────────────────────────────────────
    # Section B: E_min at C=C_Pl
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print(f"B: Binding energy at C=C_Pl={C_PL:.6f}")
    print("=" * 68)

    def pent3B_cpl(d): return E_pentagon(d, C_PL, M_SP, lut_A, lut_B)
    def pent2B_cpl(d): return E_pentagon_2B(d, C_PL, M_SP)
    def tri3B_cpl(d):  return E_equilateral(d, C_PL, M_SP)
    def tri2B_cpl(d):  return E_triangle_2B(d, C_PL, M_SP)
    def sq3B_cpl(d):   return E_square(d, C_PL, M_SP)
    def sq2B_cpl(d):   return E_square_2B(d, C_PL, M_SP)

    geoms_cpl = [
        ("Triangle (N=3) 2B",    tri2B_cpl),
        ("Triangle (N=3) 2B+3B", tri3B_cpl),
        ("Square   (N=4) 2B",    sq2B_cpl),
        ("Square   (N=4) 2B+3B", sq3B_cpl),
        ("Pentagon (N=5) 2B",    pent2B_cpl),
        ("Pentagon (N=5) 2B+3B", pent3B_cpl),
    ]

    print(f"\n  {'Geometry':<28} {'E_min [E_Pl]':>14} {'d_ref [l_Pl]':>14} {'Bound?':>7}")
    print("  " + "-" * 64)
    for name, func in geoms_cpl:
        Emin, dmin = find_Emin_scan(func, D_LO, D_HI)
        bound = "YES" if Emin < 0 else "no"
        print(f"  {name:<28} {Emin:>14.6f} {dmin:>14.3f} {bound:>7}")

    # ─────────────────────────────────────────
    # Section C: Window scan vs C
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print("C: 3B-only window scan  (E_2B>0 and E_3B<0)  vs C")
    print("=" * 68)
    C_arr = np.linspace(0.08, 0.32, 49)

    print(f"\n  {'C':>5}  {'2B':>8}  {'N=3':>8}  {'N=4':>8}  {'N=5':>8}  "
          f"{'3-only N3':>10}  {'3-only N4':>10}  {'3-only N5':>10}")
    print("  " + "-" * 85)

    for C in C_arr:
        # 2B energies (same for any geometry below? No — different N)
        # Use triangle 2B as the "pair" reference
        E2_tri,  _ = find_Emin_scan(lambda d: E_triangle_2B(d, C, M_SP))
        E3_tri,  _ = find_Emin_scan(lambda d: E_equilateral(d, C, M_SP))
        E2_sq,   _ = find_Emin_scan(lambda d: E_square_2B(d, C, M_SP))
        E3_sq,   _ = find_Emin_scan(lambda d: E_square(d, C, M_SP))
        E2_pent, _ = find_Emin_scan(lambda d: E_pentagon_2B(d, C, M_SP))
        E3_pent, _ = find_Emin_scan(lambda d: E_pentagon(d, C, M_SP, lut_A, lut_B))

        w3 = (E2_tri > 0) and (E3_tri < 0)
        w4 = (E2_sq  > 0) and (E3_sq  < 0)
        w5 = (E2_pent > 0) and (E3_pent < 0)

        print(f"  {C:>5.3f}  {E2_tri:>8.5f}  {E3_tri:>8.5f}  {E3_sq:>8.5f}  "
              f"{E3_pent:>8.5f}  "
              f"{'YES' if w3 else 'no':>10}  {'YES' if w4 else 'no':>10}  "
              f"{'YES' if w5 else 'no':>10}")

    # ─────────────────────────────────────────
    # Section D: I_Y values — how does V3 scale with N?
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print("D: Feynman integral values at equilibrium separations")
    print(f"   (m_sp={M_SP}, C={C_MID})")
    print("=" * 68)

    # Equilateral triangle d_eq
    Etri3, d_tri = find_Emin_scan(lambda d: E_equilateral(d, C_MID, M_SP))
    IY_tri = I_Y(d_tri, d_tri, d_tri, M_SP)

    # Square d_eq
    Esq3, d_sq = find_Emin_scan(lambda d: E_square(d, C_MID, M_SP))
    IY_sq_right = I_Y(d_sq, d_sq, d_sq * np.sqrt(2.0), M_SP)
    IY_sq_total = 4.0 * IY_sq_right

    # Pentagon d_eq
    Epent3, d_pent = find_Emin_scan(lambda d: E_pentagon(d, C_MID, M_SP, lut_A, lut_B))
    IY_pA = float(lut_A(d_pent))
    IY_pB = float(lut_B(d_pent))
    IY_pent_total = 5.0 * IY_pA + 5.0 * IY_pB

    print(f"\n  Triangle (N=3): d_eq={d_tri:.3f}")
    print(f"    I_Y(d,d,d)      = {IY_tri:.6f}")
    print(f"    V3 per C^3      = {-IY_tri:.6f}")
    print(f"    Num triples     = 1")

    print(f"\n  Square (N=4): d_eq={d_sq:.3f}")
    print(f"    I_Y(d,d,d*sqrt2)= {IY_sq_right:.6f}  per right-angle triple")
    print(f"    V3 total / C^3  = {-IY_sq_total:.6f}  (4 triples)")
    print(f"    Num triples     = 4")

    print(f"\n  Pentagon (N=5): d_eq={d_pent:.3f}")
    print(f"    I_Y type A (d,d,phi*d)  = {IY_pA:.6f}  (5 triples)")
    print(f"    I_Y type B (d,phi*d,phi*d)={IY_pB:.6f}  (5 triples)")
    print(f"    V3 total / C^3  = {-IY_pent_total:.6f}  (10 triples)")
    print(f"    Num triples     = 10")

    # ─────────────────────────────────────────
    # Section E: Critical coupling thresholds
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print("E: Critical coupling thresholds (classical, 2-step bisection)")
    print(f"   m_sp={M_SP}")
    print("=" * 68)

    from scipy.optimize import brentq

    def find_C_crit(E_func_of_C, C_lo=0.05, C_hi=0.70, n=40):
        """Find C where E_min(C) crosses zero."""
        def Emin_C(C):
            E, _ = find_Emin_scan(lambda d: E_func_of_C(d, C))
            return E
        # Check if zero crossing exists
        Elo = Emin_C(C_lo)
        Ehi = Emin_C(C_hi)
        if Elo < 0:
            return C_lo   # already bound at C_lo
        if Ehi > 0:
            return None   # never binds
        return brentq(Emin_C, C_lo, C_hi, xtol=1e-6)

    print("\n  N=3 equilateral:")
    C2B_tri = find_C_crit(lambda d, C: E_triangle_2B(d, C, M_SP))
    C3B_tri = find_C_crit(lambda d, C: E_equilateral(d, C, M_SP))
    print(f"    C_crit^(2B) = {C2B_tri:.4f}   (pairs bind)")
    print(f"    C_crit^(3B) = {C3B_tri:.4f}   (triple binds)")
    print(f"    Window: ({C3B_tri:.4f}, {C2B_tri:.4f}),  dC = {C2B_tri - C3B_tri:.4f}")
    print(f"    C_Pl={'inside' if C3B_tri < C_PL < C2B_tri else 'outside'} window")

    print("\n  N=4 square:")
    C2B_sq = find_C_crit(lambda d, C: E_square_2B(d, C, M_SP))
    C3B_sq = find_C_crit(lambda d, C: E_square(d, C, M_SP))
    if C2B_sq and C3B_sq:
        print(f"    C_crit^(2B) = {C2B_sq:.4f}   (pairs bind)")
        print(f"    C_crit^(3B) = {C3B_sq:.4f}   (quartet binds)")
        print(f"    Window: ({C3B_sq:.4f}, {C2B_sq:.4f}),  dC = {C2B_sq - C3B_sq:.4f}")
        print(f"    C_Pl={'inside' if C3B_sq < C_PL < C2B_sq else 'outside'} window")
    else:
        print(f"    C_crit^(2B) = {C2B_sq}")
        print(f"    C_crit^(3B) = {C3B_sq}")

    print("\n  N=5 pentagon:")
    C2B_pent = find_C_crit(lambda d, C: E_pentagon_2B(d, C, M_SP))
    C3B_pent = find_C_crit(lambda d, C: E_pentagon(d, C, M_SP, lut_A, lut_B))
    if C2B_pent and C3B_pent:
        print(f"    C_crit^(2B) = {C2B_pent:.4f}   (pairs bind)")
        print(f"    C_crit^(3B) = {C3B_pent:.4f}   (quintet binds)")
        print(f"    Window: ({C3B_pent:.4f}, {C2B_pent:.4f}),  dC = {C2B_pent - C3B_pent:.4f}")
        print(f"    C_Pl={'inside' if C3B_pent < C_PL < C2B_pent else 'outside'} window")
    else:
        print(f"    C_crit^(2B) = {C2B_pent}")
        print(f"    C_crit^(3B) = {C3B_pent}")

    # ─────────────────────────────────────────
    # Summary table
    # ─────────────────────────────────────────
    print("\n" + "=" * 68)
    print("SUMMARY TABLE")
    print("=" * 68)
    print(f"  m_sp={M_SP}, ZP=N/(8d^2) convention")
    print()
    print(f"  {'Quantity':<35} {'N=3':>10} {'N=4':>10} {'N=5':>10}")
    print("  " + "-" * 66)
    print(f"  {'Num pairs (total)':<35} {'3':>10} {'6':>10} {'10':>10}")
    print(f"  {'Num triples C(N,3)':<35} {'1':>10} {'4':>10} {'10':>10}")
    print(f"  {'Geometry':<35} {'equil.':>10} {'square':>10} {'pentagon':>10}")

    def fmt(x):
        return f"{x:.4f}" if x is not None else "---"

    print(f"  {'C_crit^(2B) (classical)':<35} {fmt(C2B_tri):>10} {fmt(C2B_sq):>10} {fmt(C2B_pent):>10}")
    print(f"  {'C_crit^(3B) (classical)':<35} {fmt(C3B_tri):>10} {fmt(C3B_sq):>10} {fmt(C3B_pent):>10}")

    dC_tri  = C2B_tri - C3B_tri if (C2B_tri and C3B_tri) else None
    dC_sq   = C2B_sq - C3B_sq   if (C2B_sq and C3B_sq)   else None
    dC_pent = C2B_pent - C3B_pent if (C2B_pent and C3B_pent) else None
    print(f"  {'Window width dC':<35} {fmt(dC_tri):>10} {fmt(dC_sq):>10} {fmt(dC_pent):>10}")

    in_tri  = "YES" if (C3B_tri  and C2B_tri  and C3B_tri  < C_PL < C2B_tri)  else "no"
    in_sq   = "YES" if (C3B_sq   and C2B_sq   and C3B_sq   < C_PL < C2B_sq)   else "no"
    in_pent = "YES" if (C3B_pent and C2B_pent and C3B_pent < C_PL < C2B_pent) else "no"
    print(f"  {'C_Pl in classical window?':<35} {in_tri:>10} {in_sq:>10} {in_pent:>10}")

    print(f"  {'d_eq at C=0.155 [l_Pl]':<35} {d_tri:>10.3f} {d_sq:>10.3f} {d_pent:>10.3f}")

    Etri2_mid, _ = find_Emin_scan(lambda d: E_triangle_2B(d, C_MID, M_SP))
    Esq2_mid, _  = find_Emin_scan(lambda d: E_square_2B(d, C_MID, M_SP))
    Epent2_mid,_ = find_Emin_scan(lambda d: E_pentagon_2B(d, C_MID, M_SP))
    print(f"  {'E_min(2B) at C=0.155 [E_Pl]':<35} {Etri2_mid:>10.5f} {Esq2_mid:>10.5f} {Epent2_mid:>10.5f}")
    print(f"  {'E_min(3B) at C=0.155 [E_Pl]':<35} {Etri3:>10.5f} {Esq3:>10.5f} {Epent3:>10.5f}")

    print()
    print(f"  phi (golden ratio, pent. diagonal/side) = {PHI:.6f}")
    print(f"  sqrt(2) (square diagonal/side)          = {np.sqrt(2):.6f}")

    print("\nDone.")
