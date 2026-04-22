"""
ex34_corrected_efimov_window.py
================================
Corrected Efimov-like quantum window scan for TGP 3-body clusters.

MOTIVATION
----------
Ex23–ex29 computed the quantum window using the 1D equilateral breathing
mode approximation:
    E_1D(d) = ZP(d) + 3*V2(d) + V3(d,d,d)
    ZP = 3/(8*d^2)    (adiabatic radial zero-point energy)

Ex30 showed the equilateral configuration is a SADDLE POINT in shape space,
not a minimum. The actual ground state lies in the isosceles submanifold
(d12=d13=a, d23=b). Ex32 computed the classical isosceles minimum. Ex33
solved the 2D Jacobi Schrödinger equation for the isosceles subspace.

Ex33 result (m_sp=0.1, C=C_Pl):
    E0_2D = -0.009 E_Pl   (vs  E0_1D = -0.007 E_Pl,  i.e. ~27% deeper)
    Upper window boundary shifts from m_sp = 0.198 (1D) to ~0.12 (2D).

THIS SCRIPT
-----------
Performs a systematic scan to:
  1. Compute E0_1D(C, m_sp)  — equilateral 1D FD (baseline, ex26 method)
  2. Compute E0_iso(C, m_sp) — isosceles 2D FD (ex33 method, corrected)
  3. Find C_crit^(2B) and C_crit^(3B) for BOTH methods via bisection.
  4. Determine whether C_Pl lies in the 2D-corrected window for each m_sp.
  5. Produce a comparison table and summary.

PHYSICS
-------
The 1D approximation gives an UPPER BOUND on E0 (energy is higher than
true minimum, i.e. less negative). The 2D result gives a better lower bound.
Because the true ground state is deeper (more negative), the binding
threshold C_crit(3B) is lower (easier to bind) in 2D than in 1D.

Key question: does the upper boundary C_crit(2B) also shift, and by how much?

PARAMETERS
----------
m_sp_list  = [0.050, 0.076, 0.080, 0.100, 0.120, 0.150, 0.180, 0.198, 0.250]
C_Pl = 1 / (2*sqrt(pi)) ~ 0.2821

METHOD
------
1D FD: same as ex26 — 1D Schrödinger on d ∈ [0.4, 30], N=2000 points.
       V_eff_1D(d) = 3*V2(d,m) + V3(d,d,d,m) [equilateral breathing mode]
       ZP is automatically included through -∂²/∂d² (not added separately).

2D FD: 2D Schrödinger in Jacobi coords (b=d23, h=height from body 1):
       T = -1/(2*mu_1) * d²/db² - 1/(2*mu_2) * d²/dh²
       mu_1 = 0.5,  mu_2 = 2/3
       a = sqrt(h² + b²/4),  d12=d13=a, d23=b
       V(b,h) = V2(a,a,b,m) + V3(a,a,b,m)
       Grid: 80x60 sparse FD matrix (manageable on CPU)

NOTE: The 2D grid (80x60=4800 points) is chosen for speed vs accuracy
tradeoff. For production results use Nb=120, Nh=80. Here Nb=80, Nh=60.

REFERENCES
----------
  ex26_numerov_bound_state.py  -- 1D FD method (baseline)
  ex33_2d_isosceles_schrodinger.py -- 2D FD method (corrected)
  nbody/_archiwum_docs/WYNIKI_SESJI_2026_03_21.md §4d, §5 — context and prior results
  ANALIZA_NBODY_INTEGRACJA.md §4.2 -- integration summary

Author: TGP project, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os
import numpy as np
from scipy.special import k0 as K0
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
C_PL = 1.0 / (2.0 * np.sqrt(np.pi))

print("=" * 72)
print("EX34: Corrected Efimov Window — 1D vs 2D Jacobi Schrödinger")
print(f"      C_Pl = {C_PL:.6f}  [Planck units]")
print("=" * 72)
print()

# ---------------------------------------------------------------------------
# Feynman integral for triple Yukawa overlap (same as ex23–ex33 convention)
# ---------------------------------------------------------------------------
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5 * (1 + _pts);  _uw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW  = np.outer(_uw, _uw)
_A1  = _UU
_A2  = _VV * (1.0 - _UU)
_A3  = (1.0 - _UU) * (1.0 - _VV)
_JAC = (1.0 - _UU)


def I_Y(d12, d13, d23, m):
    """Feynman 2D integral for triple Yukawa overlap (exact)."""
    D = _A1 * _A2 + _A1 * _A3 + _A2 * _A3
    Q = _A2 * d12**2 + _A1 * d13**2 + _A3 * d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u    = np.where(good, m * np.sqrt(Q / D), 1.0)
    val  = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)


def V2_total(a, b, m):
    """Sum of three Yukawa 2-body potentials for isosceles (d12=d13=a, d23=b).
    Each pair has coupling -C^2 (absorbed into caller).
    Returns sum_pairs exp(-m*d_ij)/d_ij (without C^2 factor).
    """
    return (2.0 * np.exp(-m * a) / a +
                  np.exp(-m * b) / b)


def V2_equil(d, m):
    """3 * Yukawa 2-body potential for equilateral (all sides = d).
    Returns 3*exp(-m*d)/d (without C^2 factor).
    """
    return 3.0 * np.exp(-m * d) / d


# ===========================================================================
# 1D Finite-Difference Schrödinger (equilateral breathing mode)
# ===========================================================================

def solve_1d_equilateral(C, m_sp, N=2000, d_lo=0.3, d_hi=30.0):
    """
    1D equilateral Schrödinger with reduced mass mu_1D = 1/2 (radial mode).

    H = -1/(2*mu) * d²/dd² + V_eff(d)

    V_eff(d) = -C^2 * V2_equil(d, m_sp) - C^3 * I_Y(d,d,d, m_sp)

    Reduced mass for symmetric breathing: mu = m_Pl/2 = 0.5 (standard).

    Returns: lowest eigenvalue E0 (E_Pl), or +inf if no bound state.
    """
    mu = 0.5
    d  = np.linspace(d_lo, d_hi, N)
    dd = d[1] - d[0]

    V = (-C**2 * V2_equil(d, m_sp)
         - C**3 * np.array([I_Y(di, di, di, m_sp) for di in d]))

    # Finite difference Hamiltonian (Dirichlet BC)
    diag  = 1.0 / (mu * dd**2) + V
    off   = -0.5 / (mu * dd**2) * np.ones(N - 1)

    H = (np.diag(diag) +
         np.diag(off, k=1) +
         np.diag(off, k=-1))

    evals = np.linalg.eigvalsh(H)
    E0 = evals[0]
    return E0


def solve_1d_equilateral_sparse(C, m_sp, N=2000, d_lo=0.3, d_hi=30.0):
    """
    Sparse version for speed (no full matrix). Uses scipy.sparse.linalg.eigsh.
    Faster for large N.
    """
    mu = 0.5
    d  = np.linspace(d_lo, d_hi, N)
    dd = d[1] - d[0]

    V = (-C**2 * V2_equil(d, m_sp)
         - C**3 * np.array([I_Y(di, di, di, m_sp) for di in d]))

    diag = 1.0 / (mu * dd**2) + V
    off  = -0.5 / (mu * dd**2)

    H = lil_matrix((N, N))
    H.setdiag(diag, 0)
    H.setdiag(off * np.ones(N - 1), 1)
    H.setdiag(off * np.ones(N - 1), -1)
    H = csr_matrix(H)

    try:
        vals, _ = eigsh(H, k=1, which='SA', tol=1e-8, maxiter=10000)
        return float(vals[0])
    except Exception:
        return np.inf


# ===========================================================================
# 2D Finite-Difference Schrödinger (isosceles Jacobi)
# ===========================================================================

def solve_2d_isosceles(C, m_sp, Nb=80, Nh=60,
                       b_lo=0.3, b_hi=20.0,
                       h_lo=0.15, h_hi=12.0):
    """
    2D Schrödinger in Jacobi coordinates for isosceles triangle (d12=d13=a, d23=b).

    Jacobi coordinates:
        b = d23   (separation of particles 2 & 3)
        h = height from particle 1 to midpoint of 2&3 = sqrt(a² - b²/4)
        a = sqrt(h² + b²/4)

    Kinetic energy:
        T = -1/(2*mu_1) d²/db² - 1/(2*mu_2) d²/dh²
        mu_1 = 0.5,  mu_2 = 2/3

    Potential:
        V(b, h) = -C² * V2_total(a, b, m_sp) - C³ * I_Y(a, a, b, m_sp)
        with a = sqrt(h² + b²/4)

    Boundary conditions: psi = 0 on all boundaries.

    Returns: lowest eigenvalue E0 (E_Pl), or +inf if no bound state.
    """
    mu1 = 0.5
    mu2 = 2.0 / 3.0

    b_arr = np.linspace(b_lo, b_hi, Nb)
    h_arr = np.linspace(h_lo, h_hi, Nh)
    db = b_arr[1] - b_arr[0]
    dh = h_arr[1] - h_arr[0]

    # Build potential on 2D grid
    V_grid = np.zeros((Nb, Nh))
    for ib, bv in enumerate(b_arr):
        for ih, hv in enumerate(h_arr):
            a = np.sqrt(hv**2 + bv**2 / 4.0)
            if a < 1e-6 or bv < 1e-6:
                V_grid[ib, ih] = 1e6  # hard wall
                continue
            v2 = -C**2 * V2_total(a, bv, m_sp)
            v3 = -C**3 * I_Y(a, a, bv, m_sp)
            V_grid[ib, ih] = v2 + v3

    # Flatten index: idx = ib * Nh + ih
    Ntot = Nb * Nh

    def idx(ib, ih):
        return ib * Nh + ih

    H = lil_matrix((Ntot, Ntot))

    coeff_b = 0.5 / (mu1 * db**2)
    coeff_h = 0.5 / (mu2 * dh**2)

    for ib in range(Nb):
        for ih in range(Nh):
            i = idx(ib, ih)
            # Diagonal
            H[i, i] = 2.0 * coeff_b + 2.0 * coeff_h + V_grid[ib, ih]
            # Off-diagonal in b direction
            if ib > 0:
                H[i, idx(ib - 1, ih)] = -coeff_b
            if ib < Nb - 1:
                H[i, idx(ib + 1, ih)] = -coeff_b
            # Off-diagonal in h direction
            if ih > 0:
                H[i, idx(ib, ih - 1)] = -coeff_h
            if ih < Nh - 1:
                H[i, idx(ib, ih + 1)] = -coeff_h

    H = csr_matrix(H)

    try:
        vals, _ = eigsh(H, k=1, which='SA', tol=1e-7, maxiter=20000)
        return float(vals[0])
    except Exception:
        return np.inf


# ===========================================================================
# Bisection to find C_crit (threshold where E0 crosses zero)
# ===========================================================================

def find_C_crit_1d(m_sp, C_lo, C_hi, tol=0.002, max_iter=20,
                   N=1500, d_lo=0.3, d_hi=28.0):
    """
    Find C_crit (2B or 3B, depending on caller) via bisection for 1D equilateral.
    Returns C such that E0_1D(C, m_sp) ≈ 0.
    If E0(C_lo) > 0 and E0(C_hi) < 0 → finds the threshold.
    """
    E_lo = solve_1d_equilateral_sparse(C_lo, m_sp, N=N, d_lo=d_lo, d_hi=d_hi)
    E_hi = solve_1d_equilateral_sparse(C_hi, m_sp, N=N, d_lo=d_lo, d_hi=d_hi)

    if E_lo * E_hi > 0:
        # Both same sign — no crossing in this interval
        return None

    for _ in range(max_iter):
        C_mid = 0.5 * (C_lo + C_hi)
        E_mid = solve_1d_equilateral_sparse(C_mid, m_sp, N=N, d_lo=d_lo, d_hi=d_hi)
        if E_mid * E_lo > 0:
            C_lo = C_mid
            E_lo = E_mid
        else:
            C_hi = C_mid
            E_hi = E_mid
        if C_hi - C_lo < tol:
            break

    return 0.5 * (C_lo + C_hi)


def find_C_crit_2d(m_sp, C_lo, C_hi, tol=0.005, max_iter=15,
                   Nb=60, Nh=45, b_hi=18.0, h_hi=10.0):
    """
    Find C_crit via bisection for 2D isosceles. Coarser grid for speed.
    """
    E_lo = solve_2d_isosceles(C_lo, m_sp, Nb=Nb, Nh=Nh, b_hi=b_hi, h_hi=h_hi)
    E_hi = solve_2d_isosceles(C_hi, m_sp, Nb=Nb, Nh=Nh, b_hi=b_hi, h_hi=h_hi)

    if E_lo * E_hi > 0:
        return None

    for _ in range(max_iter):
        C_mid = 0.5 * (C_lo + C_hi)
        E_mid = solve_2d_isosceles(C_mid, m_sp, Nb=Nb, Nh=Nh, b_hi=b_hi, h_hi=h_hi)
        if E_mid * E_lo > 0:
            C_lo = C_mid
            E_lo = E_mid
        else:
            C_hi = C_mid
            E_hi = E_mid
        if C_hi - C_lo < tol:
            break

    return 0.5 * (C_lo + C_hi)


# ===========================================================================
# Main scan
# ===========================================================================

# Reference result from ex26/ex27 (1D method, reported values)
EX26_1D_REFERENCE = {
    # m_sp: (C_Q_3B, C_Q_2B)  — from ex26/ex27, verify_all.py
    0.10: (0.191, 0.310),
    0.18: (0.267, 0.396),
}

m_sp_list = [0.050, 0.076, 0.100, 0.120, 0.150, 0.180, 0.198, 0.250]

print("─" * 72)
print("PART 1: Quick energy check at C=C_Pl for 1D and 2D methods")
print("─" * 72)
print(f"{'m_sp':>8}  {'E0_1D (E_Pl)':>16}  {'E0_2D (E_Pl)':>16}  "
      f"{'ΔE (2D-1D)':>12}  {'C_Pl bound?':>12}")
print("─" * 72)

results_energy = []
for m_sp in [0.076, 0.100, 0.120, 0.150, 0.198]:
    print(f"  m_sp={m_sp:.3f} computing 1D...", flush=True)
    E1D = solve_1d_equilateral_sparse(C_PL, m_sp, N=1500, d_lo=0.3, d_hi=28.0)
    print(f"  m_sp={m_sp:.3f} computing 2D...", flush=True)
    E2D = solve_2d_isosceles(C_PL, m_sp, Nb=70, Nh=55, b_hi=18.0, h_hi=10.0)

    dE = E2D - E1D
    bound_1D = "1D:YES" if E1D < 0 else "1D:NO "
    bound_2D = "2D:YES" if E2D < 0 else "2D:NO "
    status = f"{bound_1D}/{bound_2D}"

    print(f"  {m_sp:>8.3f}  {E1D:>+16.5f}  {E2D:>+16.5f}  {dE:>+12.5f}  {status:>12}")
    results_energy.append((m_sp, E1D, E2D))

print()
print("─" * 72)
print("PART 2: Window boundary comparison — 1D vs 2D bisection")
print("        (coarser bisection for speed; production: tol=0.001)")
print("─" * 72)
print(f"{'m_sp':>6}  {'C_Q(3B)_1D':>12}  {'C_Q(2B)_1D':>12}  "
      f"{'C_Q(3B)_2D':>12}  {'C_Q(2B)_2D':>12}  "
      f"{'C_Pl 1D?':>9}  {'C_Pl 2D?':>9}")
print("─" * 72)

summary_rows = []
for m_sp in [0.076, 0.100, 0.150, 0.198]:
    print(f"\n  m_sp={m_sp:.3f}:")

    # 1D thresholds
    print(f"    1D C_Q(3B) bisection...", flush=True)
    c3b_1d = find_C_crit_1d(m_sp, C_lo=0.08, C_hi=0.22, tol=0.003,
                              max_iter=18, N=1200, d_lo=0.3, d_hi=25.0)
    print(f"    1D C_Q(2B) bisection...", flush=True)
    c2b_1d = find_C_crit_1d(m_sp, C_lo=0.20, C_hi=0.45, tol=0.003,
                              max_iter=18, N=1200, d_lo=0.3, d_hi=25.0)

    # 2D thresholds (coarser)
    print(f"    2D C_Q(3B) bisection...", flush=True)
    c3b_2d = find_C_crit_2d(m_sp, C_lo=0.08, C_hi=0.22, tol=0.006,
                              max_iter=12, Nb=55, Nh=40, b_hi=16.0, h_hi=9.0)
    print(f"    2D C_Q(2B) bisection...", flush=True)
    c2b_2d = find_C_crit_2d(m_sp, C_lo=0.20, C_hi=0.45, tol=0.006,
                              max_iter=12, Nb=55, Nh=40, b_hi=16.0, h_hi=9.0)

    def fmt(v):
        return f"{v:.3f}" if v is not None else "  N/A "

    in_1d = (c3b_1d is not None and c2b_1d is not None and
             c3b_1d < C_PL < c2b_1d)
    in_2d = (c3b_2d is not None and c2b_2d is not None and
             c3b_2d < C_PL < c2b_2d)

    row = (m_sp,
           c3b_1d, c2b_1d, c3b_2d, c2b_2d,
           in_1d, in_2d)
    summary_rows.append(row)

    print(f"  {m_sp:>6.3f}  {fmt(c3b_1d):>12}  {fmt(c2b_1d):>12}  "
          f"{fmt(c3b_2d):>12}  {fmt(c2b_2d):>12}  "
          f"{'YES' if in_1d else 'no':>9}  {'YES' if in_2d else 'no':>9}")

print()
print("=" * 72)
print("SUMMARY: Efimov window for C_Pl — 1D vs 2D Jacobi")
print("=" * 72)
print()
print(f"  C_Pl = {C_PL:.4f} [Planck units]")
print()
print(f"  {'m_sp':>6}  {'λ=1/m':>8}  "
      f"{'ΔC_1D':>8}  {'ΔC_2D':>8}  {'1D in?':>7}  {'2D in?':>7}")
print("  " + "─" * 58)
for row in summary_rows:
    m_sp, c3b_1d, c2b_1d, c3b_2d, c2b_2d, in_1d, in_2d = row
    lam = 1.0 / m_sp
    dc_1d = (c2b_1d - c3b_1d) if (c3b_1d and c2b_1d) else None
    dc_2d = (c2b_2d - c3b_2d) if (c3b_2d and c2b_2d) else None
    def fmt2(v):
        return f"{v:.3f}" if v is not None else "  — "
    print(f"  {m_sp:>6.3f}  {lam:>8.2f}  "
          f"{fmt2(dc_1d):>8}  {fmt2(dc_2d):>8}  "
          f"{'YES' if in_1d else 'no':>7}  {'YES' if in_2d else 'no':>7}")

print()
print("─" * 72)
print("PHYSICAL INTERPRETATION")
print("─" * 72)
print("""
Key results to extract from this scan:

1. CORRECTION DIRECTION: 2D isosceles gives DEEPER binding (E0_2D < E0_1D),
   consistent with equilateral being a saddle point (ex30).

2. WINDOW SHIFT: The binding threshold C_Q(3B) typically DECREASES in 2D
   (easier to bind — the isosceles configuration finds a better minimum).
   The upper threshold C_Q(2B) may also shift.

3. C_PL MEMBERSHIP: If the 2D window is [C_Q(3B)_2D, C_Q(2B)_2D], check
   whether C_Pl = 0.2821 falls inside.

4. FALSIFIABILITY: The 2D-corrected range of m_sp for which C_Pl is in the
   quantum Efimov window is the sharpened Kill-shot K13 prediction.

5. COMPARISON WITH EX27: The 1D window was (0.076, 0.198) l_Pl^-1.
   If 2D narrows to (0.076, ~0.12), the prediction is more constrained
   but still physically meaningful.

Note: This script uses coarser grids (Nb=55-70, N=1200-1500) for reasonable
runtime. For production accuracy use:
    solve_1d_equilateral_sparse: N=3000, d_hi=30
    solve_2d_isosceles: Nb=120, Nh=80, b_hi=22, h_hi=12
""")

print("=" * 72)
print("EX34 DONE — see ANALIZA_NBODY_INTEGRACJA.md §4.2 for context")
print("=" * 72)
