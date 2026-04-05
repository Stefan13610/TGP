"""
ex29_wavefunction.py
====================
Characterise the quantum ground-state wavefunction psi_0(d) for the TGP
equilateral-triangle 3-body potential at C = C_Pl = 1/(2*sqrt(pi)) ~ 0.282.

Uses: FD Schrodinger (eigh_tridiagonal, same as ex26).
      V_eff(d) = 3/(8d^2) + V2(d;C,m) + V3(d;C,m)

Outputs:
  A. Ground-state energy E_0 vs m_sp (scan across quantum window)
  B. Detailed wavefunction at m_sp = 0.1 (window midpoint):
       - E_0, d_peak, <d>, sigma_d
       - classical turning points d_L, d_R
       - tunnelling fraction P(d<d_L or d>d_R)
       - de Broglie wavelength lambda_dB at d_peak (WKB validity check)
       - ASCII profile of |psi_0(d)|^2
  C. Wavefunction moments vs m_sp (for all m_sp inside quantum window)
  D. Comparison with classical equilibrium d_eq

Author: TGP project, 2026-03-21
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import k0 as K0
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))   # 0.282094...
M_PL = 2.176e-8                         # kg
E_PL = 1.956e9                          # J  (Planck energy)

# ─────────────────────────────────────────────
# Feynman integral (ex23/ex26 convention)
# ─────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5 * (1 + _pts); _uw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW = np.outer(_uw, _uw)
_A1 = _UU; _A2 = _VV * (1 - _UU); _A3 = (1 - _UU) * (1 - _VV); _JAC = (1 - _UU)


def I_Y_equil(d, m):
    D = _A1 * _A2 + _A1 * _A3 + _A2 * _A3
    Q = (_A1 + _A2 + _A3) * d**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)


def build_lut(m, d_lo=0.2, d_hi=60.0, n=400):
    d_arr = np.linspace(d_lo, d_hi, n)
    IY_arr = np.array([I_Y_equil(d, m) for d in d_arr])
    return interp1d(d_arr, IY_arr, kind='cubic', bounds_error=False,
                    fill_value=(IY_arr[0], 0.0))


def V2B(d, C, m):
    """ZP + V2 only."""
    return 3.0 / (8.0 * d**2) - 3.0 * C**2 * np.exp(-m * d) / d


def V3B(d, C, m, lut):
    """ZP + V2 + V3."""
    return V2B(d, C, m) - C**3 * float(lut(d))


def fd_solve(V_arr, d_grid, mu=1.0, n_states=5):
    """FD tridiagonal Hamiltonian eigenvalue solver (ex26 convention)."""
    h = d_grid[1] - d_grid[0]
    diag = 1.0 / (mu * h**2) + V_arr
    off  = -0.5 / (mu * h**2) * np.ones(len(d_grid) - 1)
    n_req = min(n_states, len(d_grid) - 1)
    vals, vecs = eigh_tridiagonal(diag, off,
                                   select='i', select_range=(0, n_req - 1))
    return vals, vecs


# ─────────────────────────────────────────────
# Wavefunction diagnostics
# ─────────────────────────────────────────────
def wavefunction_moments(d_grid, psi, mu=1.0):
    """Compute normalised |psi|^2 moments."""
    h = d_grid[1] - d_grid[0]
    prob = psi**2
    norm = np.sum(prob) * h
    prob_n = prob / norm
    psi_n  = psi / np.sqrt(norm)

    d_mean  = np.sum(d_grid * prob_n) * h
    d2_mean = np.sum(d_grid**2 * prob_n) * h
    sigma   = np.sqrt(max(d2_mean - d_mean**2, 0.0))
    d_peak  = d_grid[np.argmax(prob_n)]
    return psi_n, prob_n, d_mean, sigma, d_peak, norm


def classical_turning_points(V_arr, d_grid, E0):
    """Find left and right classical turning points where V(d) = E0."""
    above = V_arr > E0
    # Left turning point: last index where V>E0 before the well
    d_L = d_grid[0]
    for i in range(len(d_grid) - 1):
        if above[i] and not above[i + 1]:
            # Linear interpolation
            f = (E0 - V_arr[i]) / (V_arr[i + 1] - V_arr[i])
            d_L = d_grid[i] + f * (d_grid[i + 1] - d_grid[i])
            break

    # Right turning point
    d_R = d_grid[-1]
    for i in range(len(d_grid) - 1, 0, -1):
        if above[i] and not above[i - 1]:
            f = (E0 - V_arr[i]) / (V_arr[i - 1] - V_arr[i])
            d_R = d_grid[i] + f * (d_grid[i - 1] - d_grid[i])
            break

    return d_L, d_R


def tunnelling_fraction(d_grid, prob_n, d_L, d_R):
    """Fraction of probability outside classical turning points."""
    h = d_grid[1] - d_grid[0]
    inside = ((d_grid >= d_L) & (d_grid <= d_R))
    P_in = np.sum(prob_n[inside]) * h
    return 1.0 - P_in


def de_broglie_wavelength(d_peak, E0, V_arr, d_grid, mu=1.0):
    """Local de Broglie wavelength at d_peak: lambda_dB = 2pi/sqrt(2mu(E-V))."""
    # Interpolate V at d_peak
    idx = np.searchsorted(d_grid, d_peak)
    idx = max(1, min(idx, len(d_grid) - 1))
    V_at_peak = V_arr[idx]
    KE = E0 - V_at_peak
    if KE <= 0:
        return float('inf')
    return 2.0 * np.pi / np.sqrt(2.0 * mu * abs(KE))


def ascii_profile(d_grid, prob_n, width=60, height=15):
    """Print ASCII bar-chart of probability density."""
    h = d_grid[1] - d_grid[0]
    # Bin into `width` bins
    d_lo, d_hi = d_grid[0], d_grid[-1]
    bins = np.linspace(d_lo, d_hi, width + 1)
    bar_vals = np.zeros(width)
    for k in range(width):
        mask = (d_grid >= bins[k]) & (d_grid < bins[k + 1])
        bar_vals[k] = np.sum(prob_n[mask]) * h

    max_val = bar_vals.max()
    if max_val <= 0:
        return
    bar_norm = bar_vals / max_val

    print(f"  |psi_0(d)|^2  (d from {d_lo:.1f} to {d_hi:.1f} l_Pl)")
    for row in range(height, 0, -1):
        thresh = (row - 1) / height
        line = ""
        for k in range(width):
            line += "*" if bar_norm[k] >= thresh else " "
        print(f"  |{line}|")
    # x-axis ticks
    tick_spacing = width // 6
    tick_line = "  +" + "-" * width + "+"
    label_line = "  "
    d_ticks = np.linspace(d_lo, d_hi, 7)
    for k, d_t in enumerate(d_ticks):
        pos = int(k * tick_spacing)
        label = f"{d_t:.1f}"
        pad = pos - len(label_line) + 2
        label_line += " " * max(0, pad) + label
    print(tick_line)
    print(label_line)
    print(f"  (d in l_Pl units; norm = integral |psi|^2 dd over grid)")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    D_LO, D_HI, N_GRID = 0.4, 30.0, 3000

    print("=" * 66)
    print("ex29: Quantum wavefunction at C = C_Pl")
    print(f"  C_Pl = {C_PL:.6f}  =  1/(2*sqrt(pi))")
    print(f"  Grid: d in [{D_LO}, {D_HI}] l_Pl, N = {N_GRID} points")
    print("  V_eff = ZP + V2 + V3,  mu = 1 m_Pl")
    print("=" * 66)

    # ── A: E_0 vs m_sp across window ─────────────────────
    print("\n--- A: Ground-state energy E_0 at C=C_Pl vs m_sp ---")
    print(f"  Quantum window (ex27): m_sp in (0.076, 0.198)")
    print()
    M_SP_LIST = [0.05, 0.076, 0.08, 0.10, 0.12, 0.15, 0.18, 0.198, 0.20, 0.25]
    d_grid_A = np.linspace(D_LO, D_HI, N_GRID)

    print(f"  {'m_sp':>6}  {'lambda':>7}  {'E_0(2B)':>11}  {'E_0(3B)':>11}  "
          f"{'bound?':>7}  {'in window?':>11}")
    print("  " + "-" * 64)

    for m in M_SP_LIST:
        lut_m = build_lut(m, 0.2, 60.0, 400)
        V2_arr = np.array([V2B(d, C_PL, m) for d in d_grid_A])
        V3_arr = np.array([V3B(d, C_PL, m, lut_m) for d in d_grid_A])
        vals2, _ = fd_solve(V2_arr, d_grid_A, n_states=1)
        vals3, _ = fd_solve(V3_arr, d_grid_A, n_states=1)
        E0_2B = vals2[0]
        E0_3B = vals3[0]
        bound = E0_3B < 0
        in_win = (E0_2B > 0) and (E0_3B < 0)
        win_str = "**YES**" if in_win else ("bound" if bound else "no")
        print(f"  {m:>6.3f}  {1/m:>7.1f}  {E0_2B:>11.6f}  {E0_3B:>11.6f}  "
              f"{'YES' if bound else 'no':>7}  {win_str:>11}")

    # ── B: Detailed wavefunction at m_sp = 0.1 ───────────
    M_SP_B = 0.1
    print(f"\n--- B: Detailed wavefunction at m_sp={M_SP_B}, C=C_Pl ---")
    lut_B = build_lut(M_SP_B, 0.2, 60.0, 400)
    d_grid_B = np.linspace(D_LO, D_HI, N_GRID)
    V3_B = np.array([V3B(d, C_PL, M_SP_B, lut_B) for d in d_grid_B])

    vals_B, vecs_B = fd_solve(V3_B, d_grid_B, n_states=5)
    bound_idx = [i for i in range(len(vals_B)) if vals_B[i] < 0]
    print(f"  Number of bound states: {len(bound_idx)}")

    if not bound_idx:
        print("  No bound states found!")
    else:
        i0 = bound_idx[0]
        E0 = vals_B[i0]
        psi_raw = vecs_B[:, i0]
        psi_n, prob_n, d_mean, sigma, d_peak, norm = wavefunction_moments(
            d_grid_B, psi_raw)

        # Classical turning points
        d_L, d_R = classical_turning_points(V3_B, d_grid_B, E0)

        # Tunnelling fraction
        P_tunnel = tunnelling_fraction(d_grid_B, prob_n, d_L, d_R)

        # de Broglie wavelength at d_peak
        lambda_dB = de_broglie_wavelength(d_peak, E0, V3_B, d_grid_B)
        L_well = d_R - d_L

        # Classical equilibrium
        d_scan_cl = np.linspace(0.5, 20.0, 600)
        E_cl = np.array([V3B(d, C_PL, M_SP_B, lut_B) for d in d_scan_cl])
        i_cl = int(np.argmin(E_cl))
        # Refine
        d_ref = np.linspace(max(0.5, d_scan_cl[max(0, i_cl-3)]),
                            min(20.0, d_scan_cl[min(len(d_scan_cl)-1, i_cl+3)]), 300)
        E_ref = np.array([V3B(d, C_PL, M_SP_B, lut_B) for d in d_ref])
        d_cl = d_ref[int(np.argmin(E_ref))]
        E_cl_min = E_ref.min()

        # V decomposition at d_mean
        ZP_at_dm   = 3.0 / (8.0 * d_mean**2)
        V2_at_dm   = -3.0 * C_PL**2 * np.exp(-M_SP_B * d_mean) / d_mean
        V3_at_dm_raw = -C_PL**3 * float(lut_B(d_mean))
        V_eff_at_dm  = ZP_at_dm + V2_at_dm + V3_at_dm_raw

        print(f"\n  === Ground state (n=0) ===")
        print(f"  E_0           = {E0:.8f} E_Pl  =  {E0*E_PL:.3e} J")
        print(f"  d_peak        = {d_peak:.4f} l_Pl")
        print(f"  <d>           = {d_mean:.4f} l_Pl")
        print(f"  sigma_d       = {sigma:.4f} l_Pl  (width of wavepacket)")
        print(f"  sigma_d/<d>   = {sigma/d_mean:.4f}  (relative uncertainty)")

        print(f"\n  === Classical comparison ===")
        print(f"  d_cl(V_eff)   = {d_cl:.4f} l_Pl  (classical equilibrium)")
        print(f"  E_cl_min      = {E_cl_min:.6f} E_Pl")
        print(f"  <d>/d_cl      = {d_mean/d_cl:.4f}  (quantum vs classical position)")
        print(f"  E_0/E_cl_min  = {E0/E_cl_min:.4f}  (quantum vs classical energy)")

        print(f"\n  === Turning points and tunnelling ===")
        print(f"  d_L (left TP) = {d_L:.4f} l_Pl")
        print(f"  d_R (right TP)= {d_R:.4f} l_Pl")
        print(f"  L_well = d_R - d_L = {L_well:.4f} l_Pl")
        print(f"  P(tunnel)     = {P_tunnel*100:.2f}%  (prob outside classical region)")

        print(f"\n  === WKB validity at d_peak ===")
        print(f"  lambda_dB     = {lambda_dB:.4f} l_Pl  (de Broglie wavelength)")
        print(f"  L_well        = {L_well:.4f} l_Pl  (classical well width)")
        print(f"  lambda_dB/L_well = {lambda_dB/L_well:.4f}")
        if lambda_dB > 0.1 * L_well:
            print(f"  WKB INVALID (lambda_dB > 0.1 * L_well) — exact FD required.")
        else:
            print(f"  WKB marginally valid (lambda_dB << L_well).")

        print(f"\n  === Energy decomposition at <d> ===")
        print(f"  ZP(d=<d>)     = {ZP_at_dm:+.6f} E_Pl")
        print(f"  V2(d=<d>)     = {V2_at_dm:+.6f} E_Pl")
        print(f"  V3(d=<d>)     = {V3_at_dm_raw:+.6f} E_Pl")
        print(f"  V_eff(d=<d>)  = {V_eff_at_dm:+.6f} E_Pl  (vs E_0={E0:.6f})")

        print()
        ascii_profile(d_grid_B, prob_n, width=58, height=12)

        # Higher bound states
        if len(bound_idx) > 1:
            print(f"\n  Higher bound states:")
            for k in bound_idx[1:]:
                psi_k = vecs_B[:, k]
                _, prob_k, dm_k, sig_k, dp_k, _ = wavefunction_moments(d_grid_B, psi_k)
                print(f"    n={bound_idx.index(k)}: E={vals_B[k]:.6f}, "
                      f"d_peak={dp_k:.2f}, <d>={dm_k:.2f}, sigma={sig_k:.2f}")

    # ── C: Wavefunction moments vs m_sp ──────────────────
    print(f"\n--- C: Wavefunction moments vs m_sp (C = C_Pl) ---")
    M_SP_WIN = [0.08, 0.10, 0.12, 0.15, 0.18]
    d_grid_C = np.linspace(D_LO, D_HI, N_GRID)

    print(f"\n  {'m_sp':>6}  {'lambda':>7}  {'E_0':>10}  {'d_cl':>7}  "
          f"{'<d>':>7}  {'sigma':>7}  {'P_tun%':>8}  {'lam_dB':>8}")
    print("  " + "-" * 70)

    for m in M_SP_WIN:
        lut_m = build_lut(m, 0.2, 60.0, 400)
        V3_arr = np.array([V3B(d, C_PL, m, lut_m) for d in d_grid_C])
        vals_m, vecs_m = fd_solve(V3_arr, d_grid_C, n_states=3)
        bound_m = [i for i in range(len(vals_m)) if vals_m[i] < 0]
        if not bound_m:
            print(f"  {m:>6.3f}  {1/m:>7.1f}  {'unbound':>10}")
            continue
        i0 = bound_m[0]
        E0_m = vals_m[i0]
        psi_n_m, prob_m, d_mean_m, sigma_m, d_peak_m, _ = \
            wavefunction_moments(d_grid_C, vecs_m[:, i0])

        # Classical equilibrium
        v_cl = np.array([V3B(d, C_PL, m, lut_m) for d in np.linspace(0.5,20,500)])
        d_cl_m = np.linspace(0.5, 20, 500)[int(np.argmin(v_cl))]

        d_L_m, d_R_m = classical_turning_points(V3_arr, d_grid_C, E0_m)
        P_tun_m = tunnelling_fraction(d_grid_C, prob_m, d_L_m, d_R_m)
        lam_dB_m = de_broglie_wavelength(d_peak_m, E0_m, V3_arr, d_grid_C)

        print(f"  {m:>6.3f}  {1/m:>7.1f}  {E0_m:>10.6f}  {d_cl_m:>7.3f}  "
              f"{d_mean_m:>7.3f}  {sigma_m:>7.3f}  {P_tun_m*100:>8.2f}  "
              f"{lam_dB_m:>8.3f}")

    # ── D: Summary table ──────────────────────────────────
    print(f"\n--- D: Summary — ground state at C=C_Pl, m_sp=0.1 ---")

    # Re-run m_sp=0.1 for summary
    m_s = 0.1
    lut_s = build_lut(m_s)
    d_grid_s = np.linspace(D_LO, D_HI, N_GRID)
    V3_s = np.array([V3B(d, C_PL, m_s, lut_s) for d in d_grid_s])
    vals_s, vecs_s = fd_solve(V3_s, d_grid_s, n_states=3)
    bound_s = [i for i in range(len(vals_s)) if vals_s[i] < 0]

    if bound_s:
        i0s = bound_s[0]
        E0s = vals_s[i0s]
        psi_ns, prob_ns, d_ms, sig_s, dp_s, _ = \
            wavefunction_moments(d_grid_s, vecs_s[:, i0s])
        d_Ls, d_Rs = classical_turning_points(V3_s, d_grid_s, E0s)
        P_ts = tunnelling_fraction(d_grid_s, prob_ns, d_Ls, d_Rs)
        lam_s = de_broglie_wavelength(dp_s, E0s, V3_s, d_grid_s)
        L_ws = d_Rs - d_Ls

        print()
        print(f"  {'Quantity':<35} {'Value':>18}")
        print("  " + "-" * 54)
        print(f"  {'m_sp [l_Pl^-1]':<35} {m_s:>18.3f}")
        print(f"  {'C_Pl':<35} {C_PL:>18.6f}")
        print(f"  {'E_0 [E_Pl]':<35} {E0s:>18.8f}")
        print(f"  {'E_0 [J]':<35} {E0s*E_PL:>18.3e}")
        print(f"  {'d_peak [l_Pl]':<35} {dp_s:>18.4f}")
        print(f"  {'<d> [l_Pl]':<35} {d_ms:>18.4f}")
        print(f"  {'sigma_d [l_Pl]':<35} {sig_s:>18.4f}")
        print(f"  {'sigma_d/<d>':<35} {sig_s/d_ms:>18.4f}")
        print(f"  {'d_L (left turning pt) [l_Pl]':<35} {d_Ls:>18.4f}")
        print(f"  {'d_R (right turning pt) [l_Pl]':<35} {d_Rs:>18.4f}")
        print(f"  {'L_well = d_R - d_L [l_Pl]':<35} {L_ws:>18.4f}")
        print(f"  {'lambda_dB at d_peak [l_Pl]':<35} {lam_s:>18.4f}")
        print(f"  {'lambda_dB / L_well':<35} {lam_s/L_ws:>18.4f}")
        print(f"  {'Tunnelling fraction [%]':<35} {P_ts*100:>18.2f}")
        print(f"  {'N_bound states':<35} {len(bound_s):>18d}")

    print("\nDone.")
