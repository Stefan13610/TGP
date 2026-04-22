"""
ex33_2d_isosceles_schrodinger.py
=================================
Full 2D quantum Schrodinger equation for the TGP 3-body isosceles submanifold.

Coordinates:
  Jacobi coordinates for 3 equal-mass particles (mass m = m_Pl = 1):
    xi_1 = r2 - r3       (relative separation of particles 2 and 3)
    xi_2 = r1 - (r2+r3)/2  (particle 1 relative to CM of 2&3)

  Magnitudes:
    b = |xi_1| = d23            (distance between particles 2 and 3)
    h = |xi_2| = sqrt(a^2 - b^2/4)   (height from particle 1 to side 23)
    a = d12 = d13 = sqrt(h^2 + b^2/4)  (distance from particle 1 to either of 2,3)

  Reduced masses:
    mu_1 = m/2 = 0.5    (for xi_1)
    mu_2 = 2m/3 = 2/3   (for xi_2)

Kinetic energy operator (1D particles, hbar=1):
  T = -1/(2*mu_1) * d^2/db^2 - 1/(2*mu_2) * d^2/dh^2
    = -1.0 * d^2/db^2 - 0.75 * d^2/dh^2

This is the correct kinetic energy in Jacobi coordinates WITHOUT
a separate ZP term. The ZP=3/(8d^2) used in ex23-ex29 was an
adiabatic approximation for the angular modes; here it emerges
naturally from the 2D kinetic energy.

The 2D Schrodinger equation:
  [-1.0 * d^2/db^2 - 0.75 * d^2/dh^2 + V(b,h)] psi(b,h) = E * psi(b,h)

Boundary conditions: psi = 0 on all boundaries.

Grid:
  b in [b_lo, b_max],  h in [h_lo, h_max]
  N_b x N_h grid points -> sparse matrix of size (N_b*N_h)^2

Potential (isosceles triangle):
  V(b, h) = V2(a, a, b) + V3(a, a, b, m_sp)
  with a = sqrt(h^2 + b^2/4)

No ZP term — the quantum kinetic energy provides natural support against collapse.

Comparison:
  1D equilateral breathing mode (ex23-ex29): E0_1D
  2D isosceles full FD (this script): E0_2D
  Ratio E0_2D / E0_1D: measures accuracy of equilateral approximation

Parameters: C = C_Pl = 1/(2sqrt(pi)) ~ 0.2821, m_sp = 0.1

Author: TGP project, 2026-03-22
"""

import numpy as np
from scipy.special import k0 as K0
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import warnings
warnings.filterwarnings('ignore')

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))

# ─────────────────────────────────────────────
# Feynman integral (vectorised over scalar inputs)
# ─────────────────────────────────────────────
_pts25, _wts25 = np.polynomial.legendre.leggauss(25)
_up25 = 0.5 * (1 + _pts25); _uw25 = 0.5 * _wts25
_UU25, _VV25 = np.meshgrid(_up25, _up25, indexing='ij')
_WW25 = np.outer(_uw25, _uw25)
_A1_25 = _UU25
_A2_25 = _VV25 * (1 - _UU25)
_A3_25 = (1 - _UU25) * (1 - _VV25)
_JAC25 = (1 - _UU25)


def I_Y(d12, d13, d23, m):
    """Feynman 3-body Yukawa integral (scalar inputs)."""
    A1, A2, A3 = _A1_25, _A2_25, _A3_25
    D = A1 * A2 + A1 * A3 + A2 * A3
    Q = A2 * d12**2 + A1 * d13**2 + A3 * d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW25 * _JAC25 * val)


def V_iso(b, h, C, m_sp):
    """Isosceles potential: V2 + V3 at (b, h).
    a = sqrt(h^2 + b^2/4), b = d23, d12=d13=a."""
    if b <= 0 or h <= 0:
        return 1e10
    a = np.sqrt(h**2 + 0.25 * b**2)
    V2 = -C**2 * (2.0 * np.exp(-m_sp * a) / a + np.exp(-m_sp * b) / b)
    V3 = -C**3 * I_Y(a, a, b, m_sp)
    return V2 + V3


# ─────────────────────────────────────────────
# 1D breathing-mode (equilateral) for comparison
# ─────────────────────────────────────────────
def V_eff_equil(d, C, m_sp):
    """Effective equilateral potential including ZP."""
    ZP = 3.0 / (8.0 * d**2)
    V2 = -C**2 * 3.0 * np.exp(-m_sp * d) / d
    V3 = -C**3 * I_Y(d, d, d, m_sp)
    return ZP + V2 + V3


def solve_1d_equil(C, m_sp, d_lo=0.4, d_hi=30.0, N=2000):
    """1D FD Schrodinger for equilateral breathing mode."""
    d_arr = np.linspace(d_lo, d_hi, N)
    h_step = d_arr[1] - d_arr[0]
    V_arr = np.array([V_eff_equil(d, C, m_sp) for d in d_arr])
    mu = 1.0  # effective mass for breathing mode
    diag = 1.0 / (mu * h_step**2) + V_arr
    off  = -1.0 / (2.0 * mu * h_step**2) * np.ones(N - 1)
    from scipy.linalg import eigh_tridiagonal
    w, v = eigh_tridiagonal(diag, off, select='i', select_range=(0, 1))
    return w, v, d_arr


# ─────────────────────────────────────────────
# 2D FD Schrodinger in (b, h) coordinates
# ─────────────────────────────────────────────
def build_2d_hamiltonian(b_arr, h_arr, C, m_sp):
    """Build sparse 2D FD Hamiltonian in (b, h) Jacobi coords.

    Kinetic energy:
      T = -1.0 * d^2/db^2 - 0.75 * d^2/dh^2
      (reduced masses mu_b = 0.5, mu_h = 2/3 -> coefficients 1/(2*mu) = 1.0, 0.75)

    Boundary conditions: Dirichlet (psi = 0) on all edges.

    Grid: (N_b, N_h), total N = N_b * N_h unknowns (interior points only,
          but we include boundaries implicitly via the FD stencil).
    """
    N_b = len(b_arr)
    N_h = len(h_arr)
    N_tot = N_b * N_h
    db = b_arr[1] - b_arr[0]
    dh = h_arr[1] - h_arr[0]

    # Kinetic prefactors: T = -(hbar^2/2mu) d^2/dx^2
    # mu_b = m/2 = 0.5,  coefficient in T = 1/(2*0.5) = 1.0
    # mu_h = 2m/3 = 2/3, coefficient in T = 1/(2*(2/3)) = 3/4 = 0.75
    kb = 1.0 / db**2    # = 1/(2*mu_b * db^2) with mu_b=0.5
    kh = 0.75 / dh**2   # = 1/(2*mu_h * dh^2) with mu_h=2/3

    print(f"    Grid: {N_b} x {N_h} = {N_tot} points")
    print(f"    db = {db:.4f}, dh = {dh:.4f}")
    print(f"    Precomputing 2D potential V(b,h)...", end='', flush=True)

    # Precompute potential on grid (slow due to I_Y calls)
    V_grid = np.zeros((N_b, N_h))
    n_computed = 0
    for i, b in enumerate(b_arr):
        for j, h in enumerate(h_arr):
            V_grid[i, j] = V_iso(b, h, C, m_sp)
            n_computed += 1
    print(f" done ({n_computed} points)")

    print(f"    Building sparse Hamiltonian...", end='', flush=True)

    # Diagonal: kinetic diagonal + potential
    # Off-diagonals: kinetic FD couplings
    # Index mapping: (i, j) -> i * N_h + j

    diag = np.zeros(N_tot)
    off_b = np.zeros(N_tot - N_h)   # coupling in b direction
    off_h = np.zeros(N_tot - 1)      # coupling in h direction

    for i in range(N_b):
        for j in range(N_h):
            idx = i * N_h + j
            # Diagonal: kinetic sum + potential
            diag[idx] = 2.0 * kb + 2.0 * kh + V_grid[i, j]
            # Boundary: already zero from Dirichlet (no off-diag contribution)

    # b-direction off-diagonals (stride N_h)
    for i in range(N_b - 1):
        for j in range(N_h):
            idx = i * N_h + j
            off_b[idx] = -kb

    # h-direction off-diagonals (stride 1, but beware row boundaries)
    for i in range(N_b):
        for j in range(N_h - 1):
            idx = i * N_h + j
            off_h[idx] = -kh

    print(f" done")
    return diag, off_b, off_h, N_b, N_h, V_grid


def solve_2d_schrodinger(C, m_sp,
                          b_lo=0.15, b_hi=20.0, N_b=60,
                          h_lo=0.15, h_hi=20.0, N_h=60,
                          n_eigvals=3):
    """Solve 2D isosceles Schrodinger equation using sparse eigensolver.

    Returns:
      E0_2D: ground state energy
      psi0:  ground state wavefunction (N_b x N_h array)
      b_arr, h_arr: coordinate grids
    """
    b_arr = np.linspace(b_lo, b_hi, N_b)
    h_arr = np.linspace(h_lo, h_hi, N_h)
    N_tot = N_b * N_h

    diag, off_b, off_h, N_b, N_h, V_grid = build_2d_hamiltonian(
        b_arr, h_arr, C, m_sp)

    # Build sparse tridiagonal-block matrix
    print(f"    Assembling sparse matrix ({N_tot}x{N_tot})...", end='', flush=True)
    H = lil_matrix((N_tot, N_tot))

    for i in range(N_b):
        for j in range(N_h):
            idx = i * N_h + j
            H[idx, idx] = diag[idx]
            # b-direction: couple to (i+1, j)
            if i + 1 < N_b:
                idx2 = (i + 1) * N_h + j
                H[idx, idx2] = -1.0 / (b_arr[1] - b_arr[0])**2
                H[idx2, idx] = -1.0 / (b_arr[1] - b_arr[0])**2
            # h-direction: couple to (i, j+1)
            if j + 1 < N_h:
                idx2 = i * N_h + (j + 1)
                H[idx, idx2] = -0.75 / (h_arr[1] - h_arr[0])**2
                H[idx2, idx] = -0.75 / (h_arr[1] - h_arr[0])**2

    H_csr = csr_matrix(H)
    print(f" done")

    # Solve for lowest n_eigvals eigenvalues
    print(f"    Running eigsh (k={n_eigvals})...", end='', flush=True)
    try:
        eigenvalues, eigenvectors = eigsh(H_csr, k=n_eigvals, which='SA',
                                          tol=1e-8, maxiter=10000)
        eigenvalues = np.sort(eigenvalues)
        print(f" done")
    except Exception as e:
        print(f" FAILED: {e}")
        eigenvalues = np.array([float('nan')] * n_eigvals)
        eigenvectors = None

    psi0 = None
    if eigenvectors is not None:
        psi0 = eigenvectors[:, 0].reshape(N_b, N_h)

    return eigenvalues, psi0, b_arr, h_arr, V_grid


def analyze_wavefunction(psi0, b_arr, h_arr):
    """Compute wavefunction moments in (b, h) coordinates."""
    db = b_arr[1] - b_arr[0]
    dh = h_arr[1] - h_arr[0]
    prob = psi0**2
    norm = np.sum(prob) * db * dh
    prob_n = prob / norm

    B, H = np.meshgrid(b_arr, h_arr, indexing='ij')
    A = np.sqrt(H**2 + 0.25 * B**2)  # a = sqrt(h^2 + b^2/4)
    D_equil = np.sqrt(B**2 + H**2)   # rough "size"

    b_mean = np.sum(B * prob_n) * db * dh
    h_mean = np.sum(H * prob_n) * db * dh
    a_mean = np.sum(A * prob_n) * db * dh
    d_bar_mean = np.sum((2*A + B)/3 * prob_n) * db * dh

    # Peak
    ij_peak = np.unravel_index(np.argmax(prob_n), prob_n.shape)
    b_peak = b_arr[ij_peak[0]]
    h_peak = h_arr[ij_peak[1]]
    a_peak = np.sqrt(h_peak**2 + 0.25 * b_peak**2)

    return dict(b_mean=b_mean, h_mean=h_mean, a_mean=a_mean,
                d_bar_mean=d_bar_mean, b_peak=b_peak, h_peak=h_peak,
                a_peak=a_peak, norm=norm)


def run_comparison(C=C_PL, m_sp=0.1, N_b=50, N_h=50):
    """Compare 1D equilateral vs 2D isosceles ground state energies."""
    print(f"\n{'='*62}")
    print(f"  ex33: 2D Isosceles Schrodinger   C={C:.4f}  m_sp={m_sp:.4f}")
    print(f"{'='*62}")

    # ── 1D equilateral reference ──────────────────────────────────
    print(f"\n  [1] 1D equilateral breathing mode (ex23 reference):")
    w1d, v1d, d_arr = solve_1d_equil(C, m_sp)
    E0_1d = w1d[0]
    print(f"      E0_1D = {E0_1d:.6f} E_Pl")
    if len(w1d) > 1:
        print(f"      E1_1D = {w1d[1]:.6f} E_Pl")

    # ── 2D isosceles full solver ──────────────────────────────────
    print(f"\n  [2] 2D isosceles FD Schrodinger:")
    eigenvalues, psi0, b_arr2, h_arr2, V_grid = solve_2d_schrodinger(
        C, m_sp, N_b=N_b, N_h=N_h, n_eigvals=3)

    E0_2d = eigenvalues[0]
    print(f"\n  Eigenvalues:")
    for k, E in enumerate(eigenvalues):
        print(f"      E{k}_2D = {E:.6f} E_Pl")

    # ── Wavefunction analysis ─────────────────────────────────────
    if psi0 is not None:
        print(f"\n  Ground state wavefunction analysis:")
        mom = analyze_wavefunction(psi0, b_arr2, h_arr2)
        print(f"      <b>    = {mom['b_mean']:.4f} l_Pl  (d23 mean)")
        print(f"      <h>    = {mom['h_mean']:.4f} l_Pl  (height mean)")
        print(f"      <a>    = {mom['a_mean']:.4f} l_Pl  (d12=d13 mean)")
        print(f"      <d_bar>= {mom['d_bar_mean']:.4f} l_Pl")
        print(f"      peak:  b={mom['b_peak']:.4f},  h={mom['h_peak']:.4f},"
              f"  a={mom['a_peak']:.4f}")

        # On equilateral line: h = b*sqrt(3)/2 -> b_equil = 2h/sqrt(3)
        b_peak = mom['b_peak']
        h_peak = mom['h_peak']
        b_equil_at_h = h_peak * 2.0 / np.sqrt(3)
        deform = (b_peak - b_equil_at_h) / b_equil_at_h
        print(f"\n      Equilateral line: b_equil(h_peak) = {b_equil_at_h:.4f}")
        print(f"      Deformation from equilateral: {deform*100:.1f}%")
        print(f"      b/a at peak = {b_peak/mom['a_peak']:.4f} (equilateral: 1.000)")

    # ── Comparison ────────────────────────────────────────────────
    print(f"\n  {'='*40}")
    print(f"  Comparison 1D equilateral vs 2D isosceles:")
    print(f"    E0_1D (equilateral) = {E0_1d:.6f} E_Pl")
    print(f"    E0_2D (isosceles)   = {E0_2d:.6f} E_Pl")
    dE = E0_2d - E0_1d
    print(f"    dE = E0_2D - E0_1D  = {dE:.6f} E_Pl")
    if abs(E0_1d) > 1e-10:
        print(f"    dE / |E0_1D|       = {dE/abs(E0_1d)*100:.2f}%")

    if dE < 0:
        print(f"\n  => 2D ground state is LOWER: equilateral approx is UPPER BOUND")
        print(f"     (as predicted by the saddle-point result of ex30)")
    elif dE > 0:
        print(f"\n  => 2D ground state is HIGHER: boundary effects dominating")
    else:
        print(f"\n  => 1D and 2D give identical result (equilateral exact)")

    return E0_1d, E0_2d


def scan_m_sp_2d():
    """Scan E0_2D vs m_sp and compare with E0_1D."""
    print(f"\n\n{'='*62}")
    print(f"  ex33: Scan vs m_sp (N_b=N_h=40, C=C_Pl)")
    print(f"{'='*62}")
    print(f"\n  {'m_sp':>8}  {'E0_1D':>10}  {'E0_2D':>10}  {'dE':>10}  {'dE%':>8}")
    print(f"  {'-'*55}")

    m_vals = [0.076, 0.100, 0.130, 0.150, 0.180, 0.198]
    for m in m_vals:
        # 1D reference
        w1d, _, _ = solve_1d_equil(C_PL, m)
        E0_1d = w1d[0]
        # 2D full
        evals, _, _, _, _ = solve_2d_schrodinger(
            C_PL, m, N_b=40, N_h=40, n_eigvals=1)
        E0_2d = evals[0]
        dE = E0_2d - E0_1d
        pct = dE / abs(E0_1d) * 100 if abs(E0_1d) > 1e-10 else float('nan')
        print(f"  {m:8.3f}  {E0_1d:10.5f}  {E0_2d:10.5f}  {dE:10.5f}  {pct:8.2f}%")


if __name__ == '__main__':
    # Main comparison at nominal C=C_Pl, m_sp=0.1
    print("Running 2D Schrodinger solver (N=50x50)...")
    E0_1d, E0_2d = run_comparison(C=C_PL, m_sp=0.1, N_b=50, N_h=50)

    # Scan vs m_sp with smaller grid for speed
    print("\n\nRunning m_sp scan (N=40x40, slower)...")
    scan_m_sp_2d()

    print("\n\nDone.")
