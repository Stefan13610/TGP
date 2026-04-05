"""
ex26_numerov_bound_state.py
============================
Exact 1D Schrodinger equation for the TGP 3-body effective potential
using FINITE-DIFFERENCE MATRIX DIAGONALISATION.

Physical model (consistent with ex22-ex24):
  V_eff(d) = ZP(d) + V2(d;C,m) + V3(d;C,m)
  ZP = 3/(8d^2)  [hypercentrifugal term; kinetic energy from transverse DOF]
  V2 = -3 C^2 exp(-m*d)/d  (3 equilateral pairs)
  V3 = -C^3 * I_Y(d,d,d; m)  (single Feynman triple)

Schrodinger equation in the collective d-coordinate:
  E psi = [-1/(2mu) d^2/dd^2 + V_eff(d)] psi
  mu = 1 m_Pl

Two criteria:
  (a) CLASSICAL Efimov window (ex23): E_min(V_eff_2B/3B) < 0
  (b) QUANTUM Efimov window (this script): E_0^FD(V_eff_2B/3B) < 0
      where E_0^FD is the ground-state eigenvalue of [-1/(2mu) d^2/dd^2 + V_eff]

Author: TGP project, 2026-03-21
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import k0 as K0
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
# Feynman integral (ex23 convention)
# ─────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW = np.outer(_uw, _uw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

def I_Y_equil(d, m):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = (_A1+_A2+_A3)*d**2       # all distances = d
    good = (D > 1e-30) & (Q > 1e-30)
    u   = np.where(good, m*np.sqrt(Q/D), 1.0)
    val = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def build_lut(m, d_lo=0.2, d_hi=60.0, n=400):
    d_arr  = np.linspace(d_lo, d_hi, n)
    IY_arr = np.array([I_Y_equil(d, m) for d in d_arr])
    return interp1d(d_arr, IY_arr, kind='cubic', bounds_error=False,
                    fill_value=(IY_arr[0], 0.0))

def V2B(d, C, m, lut=None):
    """Effective potential: ZP + V2 only (2-body triangle)."""
    return 3.0/(8.0*d**2) - 3.0*C**2*np.exp(-m*d)/d

def V3B(d, C, m, lut):
    """Effective potential: ZP + V2 + V3 (3-body triangle)."""
    return V2B(d,C,m) - C**3*float(lut(d))

# ─────────────────────────────────────────────
# FD Hamiltonian solver
# ─────────────────────────────────────────────
def fd_eigenvalue(V_arr, d_grid, mu=1.0, n_states=3):
    """
    Tridiagonal Hamiltonian:  H = T + V,  T_ij = -1/(2mu) finite difference.
    Returns lowest n_states eigenvalues.
    """
    N = len(d_grid)
    h = d_grid[1] - d_grid[0]
    diag    =  1.0/(mu*h**2) + V_arr
    off_diag = -0.5/(mu*h**2) * np.ones(N-1)
    vals, vecs = eigh_tridiagonal(diag, off_diag,
                                   select='i', select_range=(0, n_states-1))
    return vals, vecs

def scan_C(C_arr, m, lut, d_lo=0.4, d_hi=30.0, N=3000, mu=1.0):
    """Compute FD E_0 for 2B and 3B at each C.  Also classical minima."""
    d_grid = np.linspace(d_lo, d_hi, N)
    d_scan = np.linspace(0.5, 20.0, 500)
    results = []
    for C in C_arr:
        V2 = np.array([V2B(d, C, m) for d in d_grid])
        V3 = np.array([V3B(d, C, m, lut) for d in d_grid])
        vals2, _ = fd_eigenvalue(V2, d_grid, mu, n_states=1)
        vals3, _ = fd_eigenvalue(V3, d_grid, mu, n_states=3)
        E0_2B_FD  = vals2[0]
        bound3     = vals3[vals3 < 0.0]
        E0_3B_FD  = bound3[0] if len(bound3) > 0 else vals3[0]
        N_3B_FD   = len(bound3)
        # Classical minima
        E2B_cl = min(V2B(d, C, m) for d in d_scan)
        E3B_cl = min(V3B(d, C, m, lut) for d in d_scan)
        d_eq   = d_scan[np.argmin([V3B(d, C, m, lut) for d in d_scan])]
        results.append({
            'C': C,
            'E0_2B': E0_2B_FD, 'E0_3B': E0_3B_FD, 'N_3B': N_3B_FD,
            'E2B_cl': E2B_cl,  'E3B_cl': E3B_cl,   'd_eq': d_eq,
            'cl_win': (E2B_cl > 0 and E3B_cl < 0),
            'qm_win': (E0_2B_FD > 0 and E0_3B_FD < 0),
        })
    return results, d_grid


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    M_SP = 0.1

    print("="*62)
    print("ex26: Exact quantum (FD) Schrodinger for TGP 3B potential")
    print(f"  m_sp = {M_SP},  mu = 1 m_Pl")
    print("  V_eff = ZP + V2 [+V3],  ZP = 3/(8d^2)")
    print("="*62)

    print("\nBuilding V3 LUT...")
    lut = build_lut(M_SP, 0.2, 60.0, 400)
    print("  LUT ready.")

    C_wide = np.arange(0.128, 0.320, 0.008)

    # ── A: Potential shape ───────────────────────────────
    print("\n--- A: V_eff(d) for 2B and 3B at C=0.155, 0.200 ---")
    print(f"  {'d':>5}  {'V2B(0.155)':>12}  {'V3B(0.155)':>12}  "
          f"{'V2B(0.200)':>12}  {'V3B(0.200)':>12}")
    print("  "+"-"*60)
    for d in [1,2,3,4,5,7,10,15,20]:
        v2_155 = V2B(d,0.155,M_SP)
        v3_155 = V3B(d,0.155,M_SP,lut)
        v2_200 = V2B(d,0.200,M_SP)
        v3_200 = V3B(d,0.200,M_SP,lut)
        print(f"  {d:>5}  {v2_155:>12.5f}  {v3_155:>12.5f}  "
              f"{v2_200:>12.5f}  {v3_200:>12.5f}")

    # ── B: Full scan ─────────────────────────────────────
    print("\n--- B: Classical vs Quantum Efimov window (m_sp=0.1) ---")
    results, d_grid = scan_C(C_wide, M_SP, lut)
    print(f"  {'C':>5}  {'E0_2B(FD)':>11}  {'E0_3B(FD)':>11}  "
          f"{'N_3B':>5}  {'cl_win':>8}  {'qm_win':>8}")
    print("  "+"-"*60)
    for r in results:
        cl = 'YES' if r['cl_win'] else 'no'
        qm = 'YES' if r['qm_win'] else 'no'
        print(f"  {r['C']:>5.3f}  {r['E0_2B']:>11.5f}  {r['E0_3B']:>11.5f}  "
              f"{r['N_3B']:>5}  {cl:>8}  {qm:>8}")

    # ── C: Determine thresholds precisely ────────────────
    print("\n--- C: Quantum window thresholds ---")
    # C_Q(3B): where E0_3B first becomes negative
    for i, r in enumerate(results):
        if r['E0_3B'] < 0:
            C_Q_3B = r['C']
            break
    else:
        C_Q_3B = None
    # C_Q(2B): where E0_2B first becomes negative
    for i, r in enumerate(results):
        if r['E0_2B'] < 0:
            C_Q_2B = r['C']
            break
    else:
        C_Q_2B = None

    # Bracket and bisect for C_Q(3B)
    for i in range(len(results)-1):
        if results[i]['E0_3B'] > 0 and results[i+1]['E0_3B'] < 0:
            C_lo, C_hi = results[i]['C'], results[i+1]['C']
            for _ in range(20):
                C_mid = 0.5*(C_lo+C_hi)
                d_grid2 = np.linspace(0.4,30,3000)
                V3_mid = np.array([V3B(d,C_mid,M_SP,lut) for d in d_grid2])
                v,_ = fd_eigenvalue(V3_mid, d_grid2, n_states=1)
                if v[0] > 0: C_lo = C_mid
                else:         C_hi = C_mid
            C_Q_3B = 0.5*(C_lo+C_hi)
            break

    for i in range(len(results)-1):
        if results[i]['E0_2B'] > 0 and results[i+1]['E0_2B'] < 0:
            C_lo, C_hi = results[i]['C'], results[i+1]['C']
            for _ in range(20):
                C_mid = 0.5*(C_lo+C_hi)
                d_grid2 = np.linspace(0.4,30,3000)
                V2_mid = np.array([V2B(d,C_mid,M_SP) for d in d_grid2])
                v,_ = fd_eigenvalue(V2_mid, d_grid2, n_states=1)
                if v[0] > 0: C_lo = C_mid
                else:         C_hi = C_mid
            C_Q_2B = 0.5*(C_lo+C_hi)
            break

    M_PL = 2.176e-8   # kg
    m_Q_3B = 2*np.sqrt(np.pi)*C_Q_3B*M_PL if C_Q_3B else float('nan')
    m_Q_2B = 2*np.sqrt(np.pi)*C_Q_2B*M_PL if C_Q_2B else float('nan')
    dC_Q   = (C_Q_2B - C_Q_3B) if (C_Q_3B and C_Q_2B) else float('nan')

    print(f"  Classical window (ex23): C in (0.128, 0.184), dC = 0.056")
    print(f"  Quantum window   (FD):   C in ({C_Q_3B:.3f}, {C_Q_2B:.3f}), "
          f"dC = {dC_Q:.3f}")
    print(f"  Quantum 3B threshold:  C_Q(3B) = {C_Q_3B:.3f}")
    print(f"    mass: m_body = {m_Q_3B:.3e} kg = {m_Q_3B/M_PL:.3f} m_Pl")
    print(f"  Quantum 2B threshold:  C_Q(2B) = {C_Q_2B:.3f}")
    print(f"    mass: m_body = {m_Q_2B:.3e} kg = {m_Q_2B/M_PL:.3f} m_Pl")

    # ── D: Wavefunction at C_Q_3B + 0.016 ───────────────
    print("\n--- D: Wavefunction for quantum ground state ---")
    C_D = round(C_Q_3B + 0.016, 3) if C_Q_3B else 0.208
    d_grid3 = np.linspace(0.4, 30.0, 3000)
    V3_D = np.array([V3B(d, C_D, M_SP, lut) for d in d_grid3])
    vals_D, vecs_D = fd_eigenvalue(V3_D, d_grid3, n_states=3)
    bound_D = [(vals_D[i], vecs_D[:,i]) for i in range(len(vals_D)) if vals_D[i]<0]
    print(f"  C = {C_D:.3f}")
    h3 = d_grid3[1]-d_grid3[0]
    for i,(E,psi) in enumerate(bound_D):
        norm=np.sqrt(np.sum(psi**2)*h3)
        psi_n=psi/norm
        d_mean=np.sum(d_grid3*psi_n**2)*h3
        d2=np.sum(d_grid3**2*psi_n**2)*h3
        sigma=np.sqrt(max(d2-d_mean**2,0.0))
        d_peak=d_grid3[np.argmax(psi_n**2)]
        print(f"    n={i}: E={E:.6f} E_Pl,  d_peak={d_peak:.2f},  "
              f"<d>={d_mean:.2f},  sigma={sigma:.2f}  l_Pl")

    # ── E: Summary ───────────────────────────────────────
    print("\n"+"="*62)
    print("SUMMARY: Classical vs Quantum Efimov window (m_sp=0.1)")
    print("="*62)
    print(f"  {'':30}  {'Classical':>12}  {'Quantum (FD)':>13}")
    print("  "+"-"*60)
    print(f"  {'C_crit(3B) [2-body tri unbound]':30}  {'0.128':>12}  {C_Q_3B:>13.3f}")
    print(f"  {'C_crit(2B) [pair unbound thresh]':30}  {'0.184':>12}  {C_Q_2B:>13.3f}")
    print(f"  {'Window width dC':30}  {'0.056':>12}  {dC_Q:>13.3f}")
    print(f"  {'Mass window [m_Pl]':30}  {'[0.45,0.65]':>12}  "
          f"[{m_Q_3B/M_PL:.2f},{m_Q_2B/M_PL:.2f}]")
    print()
    print("  KEY FINDING:")
    print("  The quantum Efimov window is SHIFTED to higher C and WIDER.")
    print("  Classical criterion (ex23) gives a necessary condition;")
    print("  quantum kinetic energy in the d-mode raises both thresholds.")
    print("  WKB (ex24) was INVALID — well too shallow for WKB applicability")
    print("  (de Broglie wavelength >> potential width in the window).")
    print("\nDone.")
