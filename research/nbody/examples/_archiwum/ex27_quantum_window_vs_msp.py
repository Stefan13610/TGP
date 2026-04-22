"""
ex27_quantum_window_vs_msp.py
==============================
Quantum Efimov window C_Q(3B) and C_Q(2B) as a function of scalar
screening mass m_sp.

For each m_sp the exact FD Schrodinger calculation (ex26) is repeated:
  - build LUT for I_Y(d,d,d; m_sp)
  - bisect to find C_Q(3B): lowest C where E0_FD(V_3B) < 0
  - bisect to find C_Q(2B): lowest C where E0_FD(V_2B) < 0
  - check if C_Pl = 1/(2 sqrt(pi)) ~ 0.2821 lies in (C_Q(3B), C_Q(2B))
  - compare with classical window from ex23

Also finds m_sp values where C_Pl enters / leaves the quantum window,
giving the range of scalar masses for which Planck-mass objects form
quantum 3B-induced bound states.

Physical setup:
  V_eff(d) = 3/(8d^2)  +  V2(d;C,m)  +  V3(d;C,m)  [+V3 for 3B case]
  V2 = -3 C^2 exp(-m d)/d
  V3 = -C^3 I_Y(d,d,d; m)      [equilateral Feynman integral]
  mu = 1 m_Pl,  Dirichlet BC at d_lo=0.4, d_hi=30 l_Pl,  N=2000 pts.

Author: TGP project, 2026-03-21
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import k0 as K0
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────
C_PL  = 1.0 / (2.0 * np.sqrt(np.pi))   # ≈ 0.2821
M_PL  = 2.176e-8                         # kg

# ─────────────────────────────────────────────────────────────
# Feynman integral (ex23 convention): equilateral triangle d=d12=d13=d23
# ─────────────────────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW = np.outer(_uw, _uw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

def I_Y_equil(d, m):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = (_A1+_A2+_A3)*d**2   # all distances equal d → Q = sum(Ai)*d²=d²
    good = (D > 1e-30) & (Q > 1e-30)
    u   = np.where(good, m*np.sqrt(Q/D), 1.0)
    val = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def build_lut(m, d_lo=0.2, d_hi=60.0, n=350):
    d_arr  = np.linspace(d_lo, d_hi, n)
    IY_arr = np.array([I_Y_equil(d, m) for d in d_arr])
    return interp1d(d_arr, IY_arr, kind='cubic', bounds_error=False,
                    fill_value=(IY_arr[0], 0.0))

# ─────────────────────────────────────────────────────────────
# Effective potentials
# ─────────────────────────────────────────────────────────────
def V2B(d, C, m):
    return 3.0/(8.0*d**2) - 3.0*C**2*np.exp(-m*d)/d

def V3B(d, C, m, lut):
    return V2B(d, C, m) - C**3*float(lut(d))

# ─────────────────────────────────────────────────────────────
# FD eigenvalue solver
# ─────────────────────────────────────────────────────────────
def fd_E0(V_arr, d_grid, mu=1.0):
    """Return ground-state FD energy."""
    h = d_grid[1] - d_grid[0]
    diag     =  1.0/(mu*h**2) + V_arr
    off_diag = -0.5/(mu*h**2) * np.ones(len(d_grid)-1)
    vals, _ = eigh_tridiagonal(diag, off_diag, select='i', select_range=(0,0))
    return vals[0]

# ─────────────────────────────────────────────────────────────
# Bisection to find C_Q threshold for one potential type
# ─────────────────────────────────────────────────────────────
def find_C_threshold(m, lut, use_3B,
                     C_lo=0.05, C_hi=0.70,
                     d_lo=0.4, d_hi=30.0, N=2000,
                     n_bisect=30, mu=1.0):
    """
    Find lowest C where E0_FD < 0 (onset of binding).
    use_3B = True  → V3B potential
    use_3B = False → V2B potential
    Returns None if not found in [C_lo, C_hi].
    """
    d_grid = np.linspace(d_lo, d_hi, N)

    def E0(C):
        if use_3B:
            V = np.array([V3B(d, C, m, lut) for d in d_grid])
        else:
            V = np.array([V2B(d, C, m) for d in d_grid])
        return fd_E0(V, d_grid, mu)

    # Check endpoints
    E_lo = E0(C_lo)
    E_hi = E0(C_hi)
    if E_hi > 0:
        return None   # Not bound anywhere in range
    if E_lo < 0:
        return C_lo   # Already bound at lower end

    # Bisect
    for _ in range(n_bisect):
        C_mid = 0.5*(C_lo + C_hi)
        E_mid = E0(C_mid)
        if E_mid > 0:
            C_lo = C_mid
        else:
            C_hi = C_mid
    return 0.5*(C_lo + C_hi)

# ─────────────────────────────────────────────────────────────
# Classical window threshold (from classical min of V_eff)
# ─────────────────────────────────────────────────────────────
def find_C_classical(m, lut, use_3B, C_lo=0.05, C_hi=0.70, n_bisect=30):
    """Find C where classical min of V_eff passes through 0."""
    d_scan = np.linspace(0.5, 20.0, 400)
    def Emin(C):
        if use_3B:
            vals = [V3B(d, C, m, lut) for d in d_scan]
        else:
            vals = [V2B(d, C, m) for d in d_scan]
        return min(vals)
    if Emin(C_hi) > 0:
        return None
    if Emin(C_lo) < 0:
        return C_lo
    for _ in range(n_bisect):
        C_mid = 0.5*(C_lo + C_hi)
        if Emin(C_mid) > 0: C_lo = C_mid
        else:                C_hi = C_mid
    return 0.5*(C_lo + C_hi)

# ─────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────
def _fmt(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return '---'
    return f"{x:.3f}"

# ─────────────────────────────────────────────────────────────
# Main scan
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("="*70)
    print("ex27: Quantum Efimov window C_Q(3B), C_Q(2B) vs m_sp")
    print(f"  C_Pl = 1/(2 sqrt(pi)) = {C_PL:.4f}")
    print("="*70)

    M_SP_LIST = [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20,
                 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00]

    rows = []  # store results

    for m in M_SP_LIST:
        lam = 1.0/m
        print(f"\n  m_sp = {m:.2f}  (lambda = {lam:.1f} l_Pl)")
        print(f"    Building LUT...", end=' ', flush=True)
        lut = build_lut(m, 0.2, 60.0, 350)
        print("done.")

        # Classical thresholds
        Ccl_3B = find_C_classical(m, lut, use_3B=True,  C_lo=0.05, C_hi=0.70)
        Ccl_2B = find_C_classical(m, lut, use_3B=False, C_lo=0.05, C_hi=0.70)

        # Quantum thresholds
        Cq_3B  = find_C_threshold(m, lut, use_3B=True,  C_lo=0.05, C_hi=0.70)
        Cq_2B  = find_C_threshold(m, lut, use_3B=False, C_lo=0.05, C_hi=0.70)

        # Window widths
        dC_cl  = (Ccl_2B - Ccl_3B) if (Ccl_3B is not None and Ccl_2B is not None) else float('nan')
        dC_q   = (Cq_2B  - Cq_3B)  if (Cq_3B  is not None and Cq_2B  is not None) else float('nan')

        # Is C_Pl in quantum window?
        in_window = (Cq_3B is not None and Cq_2B is not None
                     and Cq_3B < C_PL < Cq_2B)

        # Mass range
        def to_mass(C):
            return 2*np.sqrt(np.pi)*C*M_PL if C is not None else float('nan')

        row = dict(
            m=m, lam=lam,
            Ccl_3B=Ccl_3B, Ccl_2B=Ccl_2B, dC_cl=dC_cl,
            Cq_3B=Cq_3B,   Cq_2B=Cq_2B,   dC_q=dC_q,
            in_window=in_window,
            m_q3B=to_mass(Cq_3B)/M_PL, m_q2B=to_mass(Cq_2B)/M_PL,
        )
        rows.append(row)

        flag_cl = 'YES' if (Ccl_3B is not None and Ccl_2B is not None and
                             Ccl_3B < C_PL < Ccl_2B) else 'no'
        flag_q  = '**YES**' if in_window else 'no'
        print(f"    Classical: C_cl(3B)={_fmt(Ccl_3B)}, C_cl(2B)={_fmt(Ccl_2B)}, "
              f"dC={dC_cl:.3f}   C_Pl in? {flag_cl}")
        print(f"    Quantum:   C_q(3B) ={_fmt(Cq_3B)},  C_q(2B) ={_fmt(Cq_2B)},  "
              f"dC={dC_q:.3f}   C_Pl in? {flag_q}")

    # ── Table output ─────────────────────────────────────────
    print("\n\n" + "="*70)
    print("SUMMARY TABLE: Quantum Efimov window vs m_sp")
    print("="*70)
    hdr = (f"{'m_sp':>6} {'lam':>6} "
           f"{'C_cl(3B)':>9} {'C_cl(2B)':>9} {'dC_cl':>7} "
           f"{'C_q(3B)':>8} {'C_q(2B)':>8} {'dC_q':>7} "
           f"{'dC_q/dC_cl':>10} {'C_Pl_in?':>10} {'mass [m_Pl]':>15}")
    print(hdr)
    print("-"*len(hdr))
    for r in rows:
        ratio = r['dC_q']/r['dC_cl'] if (np.isfinite(r['dC_cl']) and r['dC_cl']>0) else float('nan')
        flag = '**YES**' if r['in_window'] else 'no'
        mrange = f"[{r['m_q3B']:.2f},{r['m_q2B']:.2f}]" if np.isfinite(r['m_q3B']) else '---'
        print(f"{r['m']:>6.2f} {r['lam']:>6.1f} "
              f"{_fmt(r['Ccl_3B']):>9} {_fmt(r['Ccl_2B']):>9} {r['dC_cl']:>7.3f} "
              f"{_fmt(r['Cq_3B']):>8} {_fmt(r['Cq_2B']):>8} {r['dC_q']:>7.3f} "
              f"{ratio:>10.2f} {flag:>10} {mrange:>15}")

    # ── Key result: range of m_sp where C_Pl in window ───────
    print("\n")
    win_rows = [r for r in rows if r['in_window']]
    if win_rows:
        m_min = min(r['m'] for r in win_rows)
        m_max = max(r['m'] for r in win_rows)
        print(f"  C_Pl = {C_PL:.4f} lies inside quantum window for:")
        print(f"    m_sp ∈ [{m_min:.2f}, {m_max:.2f}]  l_Pl^-1")
        print(f"    lambda in [{1/m_max:.1f}, {1/m_min:.1f}]  l_Pl")
    else:
        print("  C_Pl not inside quantum window for any m_sp scanned.")

    print("\nDone.")
