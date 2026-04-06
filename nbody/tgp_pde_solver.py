#!/usr/bin/env python3
"""
tgp_pde_solver.py — Full 3D PDE solver for TGP N-body field (Layer 5)
=======================================================================

Solves the static nonlinear TGP field equation on a 3D grid with N
point sources. This is the "Layer 5" approach: no EFT expansion, no
pairwise/3-body decomposition — the full nonlinear field is computed
directly from the PDE.

Field equation for g = Φ/Φ₀ with N point sources (full kinetic mode):

    g²∇²g + g|∇g|² = g²(β g² - γ g³) - Σ Cᵢ δ³(x - xᵢ)

Substitution u = g³ linearizes the kinetic operator:

    ∇²u = 3 u^{2/3} (β u^{1/3} - γ u^{2/3}) - 3 Σ Cᵢ δ_σ(x - xᵢ)

(See tgp_strong_field_solver.py lines 22-28 for 1D derivation.)

In vacuum (β = γ):
    ∇²u = 3 u^{2/3} (1 - u^{1/3})       [β = γ = 1 units]
         = 3 (u^{2/3} - u)

BC: u = 1 on all boundaries (vacuum).

Two-tier solver strategy:
  1. FFT Picard iteration (fast, O(N³ log N) per step)
  2. Newton-Krylov fallback (robust, scipy)

Author: Mateusz Serafin (with Claude)
Date: April 2026
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, List
from scipy.fft import dstn, idstn
import warnings


# ── Grid ──────────────────────────────────────────────────────────────────

@dataclass
class TGPGrid:
    """3D uniform Cartesian grid for PDE solver.

    The grid covers [-Lx/2, Lx/2]³ with Nx interior points per axis.
    Boundary is at the edges of this box (Dirichlet u = 1).

    Interior points are at x = -Lx/2 + (i+1)*dx, i = 0..Nx-1
    where dx = Lx / (Nx + 1).
    """
    Nx: int = 128
    Ny: int = 0  # 0 means same as Nx
    Nz: int = 0
    Lx: float = 20.0
    Ly: float = 0.0  # 0 means same as Lx
    Lz: float = 0.0

    def __post_init__(self):
        if self.Ny == 0:
            self.Ny = self.Nx
        if self.Nz == 0:
            self.Nz = self.Nx
        if self.Ly == 0.0:
            self.Ly = self.Lx
        if self.Lz == 0.0:
            self.Lz = self.Lx

        self.dx = self.Lx / (self.Nx + 1)
        self.dy = self.Ly / (self.Ny + 1)
        self.dz = self.Lz / (self.Nz + 1)

        # 1D coordinate arrays (interior points only)
        self.x1d = -self.Lx / 2 + (np.arange(self.Nx) + 1) * self.dx
        self.y1d = -self.Ly / 2 + (np.arange(self.Ny) + 1) * self.dy
        self.z1d = -self.Lz / 2 + (np.arange(self.Nz) + 1) * self.dz

        # 3D meshgrids (ij indexing: X[i,j,k] varies with first axis)
        self.X, self.Y, self.Z = np.meshgrid(
            self.x1d, self.y1d, self.z1d, indexing='ij'
        )

    @property
    def dV(self) -> float:
        """Volume element."""
        return self.dx * self.dy * self.dz

    @property
    def shape(self) -> Tuple[int, int, int]:
        return (self.Nx, self.Ny, self.Nz)


# ── Source terms ──────────────────────────────────────────────────────────

def gaussian_source(grid: TGPGrid, pos: np.ndarray, C: float,
                    sigma: float) -> np.ndarray:
    """Regularized delta function: Gaussian blob on the grid.

    S(x) = C / (2π σ²)^{3/2} · exp(-|x - pos|² / (2σ²))

    Normalized so that ∫ S d³x ≈ C (to machine precision for σ >> dx).

    Parameters
    ----------
    grid : TGPGrid
    pos : array, shape (3,)
        Source position.
    C : float
        Source strength (integral of the blob).
    sigma : float
        Smoothing width.

    Returns
    -------
    S : ndarray, shape grid.shape
        Source density on the grid.
    """
    r2 = (grid.X - pos[0])**2 + (grid.Y - pos[1])**2 + (grid.Z - pos[2])**2
    norm = (2.0 * np.pi * sigma**2) ** 1.5
    return C * np.exp(-r2 / (2.0 * sigma**2)) / norm


def build_source_field(grid: TGPGrid,
                       positions: np.ndarray,
                       C_values: np.ndarray,
                       sigma: float) -> np.ndarray:
    """Sum of Gaussian-smoothed point sources.

    Parameters
    ----------
    grid : TGPGrid
    positions : array, shape (N, 3)
        Source positions.
    C_values : array, shape (N,)
        Source strengths.
    sigma : float
        Smoothing width (same for all sources).

    Returns
    -------
    S : ndarray, shape grid.shape
        Total source density.
    """
    S = np.zeros(grid.shape)
    for i in range(len(C_values)):
        S += gaussian_source(grid, positions[i], C_values[i], sigma)
    return S


# ── Finite difference operators ──────────────────────────────────────────

def laplacian_3d(u: np.ndarray, dx: float, dy: float = 0.0,
                 dz: float = 0.0, bc_val: float = 1.0) -> np.ndarray:
    """7-point finite difference Laplacian on a 3D grid.

    Dirichlet BC: u = bc_val on all boundaries.

    Parameters
    ----------
    u : ndarray, shape (Nx, Ny, Nz)
        Interior field values.
    dx, dy, dz : float
        Grid spacings. If dy=0 or dz=0, assumed equal to dx.
    bc_val : float
        Boundary value (default 1.0 for vacuum).

    Returns
    -------
    lap : ndarray, shape (Nx, Ny, Nz)
        Discrete Laplacian ∇²u.
    """
    if dy == 0.0:
        dy = dx
    if dz == 0.0:
        dz = dx

    lap = np.zeros_like(u)

    # x-direction
    idx2 = 1.0 / dx**2
    lap[1:, :, :] += u[:-1, :, :] * idx2
    lap[:-1, :, :] += u[1:, :, :] * idx2
    lap[:, :, :] -= 2.0 * u * idx2
    # BC at x boundaries
    lap[0, :, :] += bc_val * idx2
    lap[-1, :, :] += bc_val * idx2

    # y-direction
    idy2 = 1.0 / dy**2
    lap[:, 1:, :] += u[:, :-1, :] * idy2
    lap[:, :-1, :] += u[:, 1:, :] * idy2
    lap[:, :, :] -= 2.0 * u * idy2
    lap[:, 0, :] += bc_val * idy2
    lap[:, -1, :] += bc_val * idy2

    # z-direction
    idz2 = 1.0 / dz**2
    lap[:, :, 1:] += u[:, :, :-1] * idz2
    lap[:, :, :-1] += u[:, :, 1:] * idz2
    lap[:, :, :] -= 2.0 * u * idz2
    lap[:, :, 0] += bc_val * idz2
    lap[:, :, -1] += bc_val * idz2

    return lap


def gradient_sq_3d(f: np.ndarray, dx: float, dy: float = 0.0,
                   dz: float = 0.0, bc_val: float = 1.0) -> np.ndarray:
    """Squared gradient |∇f|² via central differences.

    Uses Dirichlet BC (f = bc_val at boundaries) for ghost cells.

    Parameters
    ----------
    f : ndarray, shape (Nx, Ny, Nz)
    dx, dy, dz : float
    bc_val : float

    Returns
    -------
    grad_sq : ndarray
    """
    if dy == 0.0:
        dy = dx
    if dz == 0.0:
        dz = dx

    Nx, Ny, Nz = f.shape
    grad_sq = np.zeros_like(f)

    # x-component: central diff with BC ghost
    fx = np.empty_like(f)
    fx[1:-1, :, :] = (f[2:, :, :] - f[:-2, :, :]) / (2.0 * dx)
    fx[0, :, :] = (f[1, :, :] - bc_val) / (2.0 * dx)
    fx[-1, :, :] = (bc_val - f[-2, :, :]) / (2.0 * dx)
    grad_sq += fx**2

    # y-component
    fy = np.empty_like(f)
    fy[:, 1:-1, :] = (f[:, 2:, :] - f[:, :-2, :]) / (2.0 * dy)
    fy[:, 0, :] = (f[:, 1, :] - bc_val) / (2.0 * dy)
    fy[:, -1, :] = (bc_val - f[:, -2, :]) / (2.0 * dy)
    grad_sq += fy**2

    # z-component
    fz = np.empty_like(f)
    fz[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
    fz[:, :, 0] = (f[:, :, 1] - bc_val) / (2.0 * dz)
    fz[:, :, -1] = (bc_val - f[:, :, -2]) / (2.0 * dz)
    grad_sq += fz**2

    return grad_sq


def gradient_3d(f: np.ndarray, dx: float, dy: float = 0.0,
                dz: float = 0.0,
                bc_val: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gradient components (df/dx, df/dy, df/dz) via central differences.

    Uses Dirichlet BC (f = bc_val at boundaries) for ghost cells.

    Parameters
    ----------
    f : ndarray, shape (Nx, Ny, Nz)
    dx, dy, dz : float
    bc_val : float

    Returns
    -------
    (fx, fy, fz) : tuple of ndarray
        Gradient components along each axis.
    """
    if dy == 0.0:
        dy = dx
    if dz == 0.0:
        dz = dx

    fx = np.empty_like(f)
    fx[1:-1, :, :] = (f[2:, :, :] - f[:-2, :, :]) / (2.0 * dx)
    fx[0, :, :] = (f[1, :, :] - bc_val) / (2.0 * dx)
    fx[-1, :, :] = (bc_val - f[-2, :, :]) / (2.0 * dx)

    fy = np.empty_like(f)
    fy[:, 1:-1, :] = (f[:, 2:, :] - f[:, :-2, :]) / (2.0 * dy)
    fy[:, 0, :] = (f[:, 1, :] - bc_val) / (2.0 * dy)
    fy[:, -1, :] = (bc_val - f[:, -2, :]) / (2.0 * dy)

    fz = np.empty_like(f)
    fz[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
    fz[:, :, 0] = (f[:, :, 1] - bc_val) / (2.0 * dz)
    fz[:, :, -1] = (bc_val - f[:, :, -2]) / (2.0 * dz)

    return fx, fy, fz


def interpolate_gradient_at(u: np.ndarray, grid: 'TGPGrid',
                            pos: np.ndarray) -> np.ndarray:
    """Interpolate gradient of u at an arbitrary position inside the grid.

    Uses trilinear interpolation of the central-difference gradient
    components to the target position.

    Parameters
    ----------
    u : ndarray, shape grid.shape
        Field array.
    grid : TGPGrid
    pos : array, shape (3,)
        Target position (x, y, z).

    Returns
    -------
    grad_u : ndarray, shape (3,)
        Interpolated gradient (du/dx, du/dy, du/dz) at pos.
    """
    fx, fy, fz = gradient_3d(u, grid.dx, grid.dy, grid.dz, bc_val=1.0)

    # Find fractional grid indices
    # Grid: x1d[i] = -Lx/2 + (i+1)*dx, so i = (x - (-Lx/2))/dx - 1
    ix_f = (pos[0] - grid.x1d[0]) / grid.dx
    iy_f = (pos[1] - grid.y1d[0]) / grid.dy
    iz_f = (pos[2] - grid.z1d[0]) / grid.dz

    # Clamp to valid range [0, N-2] for trilinear
    Nx, Ny, Nz = u.shape
    ix_f = np.clip(ix_f, 0.0, Nx - 1.001)
    iy_f = np.clip(iy_f, 0.0, Ny - 1.001)
    iz_f = np.clip(iz_f, 0.0, Nz - 1.001)

    ix0 = int(ix_f)
    iy0 = int(iy_f)
    iz0 = int(iz_f)
    ix1 = min(ix0 + 1, Nx - 1)
    iy1 = min(iy0 + 1, Ny - 1)
    iz1 = min(iz0 + 1, Nz - 1)

    # Weights
    wx = ix_f - ix0
    wy = iy_f - iy0
    wz = iz_f - iz0

    result = np.zeros(3)
    for c, comp in enumerate([fx, fy, fz]):
        # Trilinear interpolation
        c000 = comp[ix0, iy0, iz0]
        c100 = comp[ix1, iy0, iz0]
        c010 = comp[ix0, iy1, iz0]
        c110 = comp[ix1, iy1, iz0]
        c001 = comp[ix0, iy0, iz1]
        c101 = comp[ix1, iy0, iz1]
        c011 = comp[ix0, iy1, iz1]
        c111 = comp[ix1, iy1, iz1]

        result[c] = (c000 * (1 - wx) * (1 - wy) * (1 - wz)
                     + c100 * wx * (1 - wy) * (1 - wz)
                     + c010 * (1 - wx) * wy * (1 - wz)
                     + c110 * wx * wy * (1 - wz)
                     + c001 * (1 - wx) * (1 - wy) * wz
                     + c101 * wx * (1 - wy) * wz
                     + c011 * (1 - wx) * wy * wz
                     + c111 * wx * wy * wz)

    return result


def force_on_source(u: np.ndarray, grid: 'TGPGrid',
                    source_pos: np.ndarray,
                    C_source: float) -> np.ndarray:
    """Force on a point source from the physical field g = u^{1/3}.

    F = 4*pi * C * grad(g) at the source position.

    Uses grad(g), NOT grad(u). Working with u = g^3 would introduce
    O(C^2) cross-terms from the cube that dominate the O(C^3) three-body
    force in subtractions. Working with g directly ensures linear
    superposition holds at leading order, so 3-body subtraction cleanly
    isolates the nonlinear residual.

    The 4*pi coupling comes from the 3D Yukawa Green's function:
    V2 ~ 4*pi*C1*C2*exp(-md)/d at leading Born order.

    The self-field gradient vanishes at the Gaussian center by symmetry.
    In the 3-body force extraction, the self-field cancels exactly.

    Parameters
    ----------
    u : ndarray, shape grid.shape
    grid : TGPGrid
    source_pos : array, shape (3,)
    C_source : float

    Returns
    -------
    F : ndarray, shape (3,)
        Force vector on the source.
    """
    g = np.maximum(u, 1e-30) ** (1.0 / 3.0)
    grad_g = interpolate_gradient_at(g, grid, source_pos)
    return 4.0 * np.pi * C_source * grad_g


def extract_3body_force(grid: 'TGPGrid',
                        positions: np.ndarray,
                        C_values: np.ndarray,
                        beta: float, gamma: float,
                        sigma: Optional[float] = None,
                        picard_kw: Optional[dict] = None,
                        verbose: bool = False) -> dict:
    """Extract irreducible 3-body force via PDE gradient subtraction.

    Solves 3 PDE problems:
      u_123: all 3 sources
      u_01:  sources 0 and 1
      u_02:  sources 0 and 2

    Returns F_3body = F_0(u_123) - F_0(u_01) - F_0(u_02)

    The self-field of source 0 cancels exactly in the subtraction
    (same Gaussian at same position in all three solves).

    Parameters
    ----------
    grid : TGPGrid
    positions : array, shape (3, 3)
    C_values : array, shape (3,)
    beta, gamma : float
    sigma : float or None
    picard_kw : dict or None
    verbose : bool

    Returns
    -------
    result : dict with keys:
        'F_3body': ndarray(3,) -- irreducible 3-body force on particle 0
        'F_total': ndarray(3,) -- total force on particle 0 from 3-source field
        'F_pair_01': ndarray(3,) -- force on 0 from (0,1) pair field
        'F_pair_02': ndarray(3,) -- force on 0 from (0,2) pair field
        'all_converged': bool
    """
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    assert len(C_values) == 3

    if sigma is None:
        sigma = 1.5 * grid.dx
    kw = dict(sigma=sigma, verbose=False)
    if picard_kw:
        kw['picard_kw'] = picard_kw

    converged_all = True

    # Solve with all 3 sources
    u_123, info_123 = solve_field(grid, positions, C_values, beta, gamma, **kw)
    converged_all &= info_123['converged']

    # Solve with pair (0, 1)
    u_01, info_01 = solve_field(grid, positions[[0, 1]], C_values[[0, 1]],
                                beta, gamma, **kw)
    converged_all &= info_01['converged']

    # Solve with pair (0, 2)
    u_02, info_02 = solve_field(grid, positions[[0, 2]], C_values[[0, 2]],
                                beta, gamma, **kw)
    converged_all &= info_02['converged']

    # Forces on particle 0
    r0 = positions[0]
    C0 = C_values[0]
    F_total = force_on_source(u_123, grid, r0, C0)
    F_pair_01 = force_on_source(u_01, grid, r0, C0)
    F_pair_02 = force_on_source(u_02, grid, r0, C0)

    F_3body = F_total - F_pair_01 - F_pair_02

    if verbose:
        print(f"    F_total  = [{F_total[0]:.4e}, {F_total[1]:.4e}, {F_total[2]:.4e}]")
        print(f"    F_pair01 = [{F_pair_01[0]:.4e}, {F_pair_01[1]:.4e}, {F_pair_01[2]:.4e}]")
        print(f"    F_pair02 = [{F_pair_02[0]:.4e}, {F_pair_02[1]:.4e}, {F_pair_02[2]:.4e}]")
        print(f"    F_3body  = [{F_3body[0]:.4e}, {F_3body[1]:.4e}, {F_3body[2]:.4e}]")

    return {
        'F_3body': F_3body,
        'F_total': F_total,
        'F_pair_01': F_pair_01,
        'F_pair_02': F_pair_02,
        'all_converged': converged_all,
    }


# ── Residual and RHS ─────────────────────────────────────────────────────

def rhs_u(u: np.ndarray, S: np.ndarray,
          beta: float, gamma: float) -> np.ndarray:
    """Right-hand side of nabla^2 u = RHS for the u = g^3 formulation.

    Derivation from the full TGP field equation (ODE form):

      nabla^2 g + 2|grad g|^2/g = gamma g^3 - beta g^2

    The identity (1/3g^2) nabla^2(g^3) = nabla^2 g + 2|grad g|^2/g
    gives:
      (1/3g^2) nabla^2 u = gamma g^3 - beta g^2
      nabla^2 u = 3 gamma g^5 - 3 beta g^4
                = 3 gamma u^{5/3} - 3 beta u^{4/3}

    With point source (smoothed):
      nabla^2 u = 3 gamma u^{5/3} - 3 beta u^{4/3} - 4pi * 3 S
                = 3 gamma u^{5/3} - 3 beta u^{4/3} - 12 pi S

    The 4pi factor comes from the 3D Green's function normalization:
    nabla^2(exp(-r)/r) = exp(-r)/r - 4pi delta^3(x).

    For vacuum (beta = gamma = 1):
      nabla^2 u = 3 u^{4/3}(u^{1/3} - 1) - 12 pi S

    Linearized near u = 1 (w = u - 1):
      nabla^2 w = w - 12 pi S   =>  (nabla^2 - 1)w = -12 pi S  (Yukawa!)

    Parameters
    ----------
    u : ndarray
        Current field u = g^3 (interior points).
    S : ndarray
        Source field (smoothed delta functions, integral = C_i).
    beta, gamma : float
        TGP couplings.

    Returns
    -------
    rhs : ndarray
        Value of RHS = nabla^2 u at convergence.
    """
    u_safe = np.maximum(u, 1e-30)
    u13 = u_safe ** (1.0 / 3.0)
    u43 = u_safe * u13              # u^{4/3}
    u53 = u43 * u13                 # u^{5/3}
    return 3.0 * gamma * u53 - 3.0 * beta * u43 - 12.0 * np.pi * S


def residual_u(u: np.ndarray, grid: TGPGrid, S: np.ndarray,
               beta: float, gamma: float) -> np.ndarray:
    """Residual F(u) = ∇²u - RHS(u).

    At convergence, F(u) = 0.

    Parameters
    ----------
    u : ndarray, shape grid.shape
    grid : TGPGrid
    S : ndarray, shape grid.shape
        Source field.
    beta, gamma : float

    Returns
    -------
    F : ndarray, shape grid.shape
        Residual; ||F||→0 at convergence.
    """
    lap_u = laplacian_3d(u, grid.dx, grid.dy, grid.dz, bc_val=1.0)
    return lap_u - rhs_u(u, S, beta, gamma)


# ── Initial guess ─────────────────────────────────────────────────────────

def initial_guess_u(grid: TGPGrid, positions: np.ndarray,
                    C_values: np.ndarray,
                    beta: float, gamma: float) -> np.ndarray:
    """Initial guess for u = g³ via superposition of Yukawa profiles.

    g(x) ≈ 1 + Σ Cᵢ exp(-m_sp |x - xᵢ|) / |x - xᵢ|

    Then u_guess = g³, clamped to [1e-10, ∞).

    Parameters
    ----------
    grid : TGPGrid
    positions : array, shape (N, 3)
    C_values : array, shape (N,)
    beta, gamma : float

    Returns
    -------
    u0 : ndarray, shape grid.shape
    """
    m_sp = np.sqrt(3.0 * gamma - 2.0 * beta)
    g = np.ones(grid.shape)

    for i in range(len(C_values)):
        r = np.sqrt(
            (grid.X - positions[i, 0])**2 +
            (grid.Y - positions[i, 1])**2 +
            (grid.Z - positions[i, 2])**2
        )
        r_safe = np.maximum(r, 1e-12)
        g += C_values[i] * np.exp(-m_sp * r_safe) / r_safe

    g = np.maximum(g, 1e-10)
    return g ** 3


# ── FFT Picard solver (fast) ─────────────────────────────────────────────

def _dst_laplacian_eigenvalues(grid: TGPGrid) -> np.ndarray:
    """Eigenvalues of the discrete Laplacian for DST-I (Dirichlet BC).

    For a grid with spacing dx and N interior points, the DST-I
    eigenvalues are:
        λ_k = -2(1 - cos(π k / (N+1))) / dx²

    Combined for 3D:
        λ_{ijk} = λ_i^x + λ_j^y + λ_k^z

    Parameters
    ----------
    grid : TGPGrid

    Returns
    -------
    eigenvalues : ndarray, shape grid.shape
        Negative-definite eigenvalues of ∇².
    """
    # Mode indices: 1, 2, ..., N (DST-I convention)
    kx = np.arange(1, grid.Nx + 1)
    ky = np.arange(1, grid.Ny + 1)
    kz = np.arange(1, grid.Nz + 1)

    lx = -2.0 * (1.0 - np.cos(np.pi * kx / (grid.Nx + 1))) / grid.dx**2
    ly = -2.0 * (1.0 - np.cos(np.pi * ky / (grid.Ny + 1))) / grid.dy**2
    lz = -2.0 * (1.0 - np.cos(np.pi * kz / (grid.Nz + 1))) / grid.dz**2

    # Broadcasting: (Nx,1,1) + (1,Ny,1) + (1,1,Nz)
    return lx[:, None, None] + ly[None, :, None] + lz[None, None, :]


def _dst_invert_laplacian(rhs: np.ndarray, eigenvalues: np.ndarray,
                          bc_val: float, grid: TGPGrid) -> np.ndarray:
    """Solve ∇²u = rhs with Dirichlet BC u = bc_val, via DST.

    The DST-I diagonalizes the discrete Laplacian with Dirichlet BC.
    We solve for v = u - bc_val (homogeneous BC), then u = v + bc_val.

    ∇²u = ∇²v + ∇²(bc_val) = ∇²v   (bc_val is constant → ∇² = 0 in interior)

    But the FD Laplacian of a constant field bc_val IS zero in the interior
    only if BC matches. Since our Laplacian uses bc_val at boundaries,
    ∇²(bc_val) = 0 for interior points. So:

    ∇²v = rhs, with v = 0 on boundary → invert via DST → u = v + bc_val.

    Actually, more carefully: the FD laplacian applied to u involves the BC.
    Let's define the problem directly. We want u satisfying:
        L·u = rhs  (interior), u|_∂ = bc_val

    Substituting v = u - bc_val:
        L·v = rhs - L·(bc_val)_interior

    For constant bc_val, L·(bc_val)_interior:
        At each interior point, the 7-point stencil on a constant = 0
        EXCEPT at boundary-adjacent points where it picks up bc_val from BC.
        But the constant field already equals bc_val everywhere, so L·const = 0.

    Therefore L·v = rhs with v|_∂ = 0. Invert via DST.

    Parameters
    ----------
    rhs : ndarray, shape grid.shape
        RHS of Poisson equation ∇²u = rhs.
    eigenvalues : ndarray, shape grid.shape
        DST eigenvalues (negative-definite).
    bc_val : float
        Boundary value for u.
    grid : TGPGrid

    Returns
    -------
    u : ndarray, shape grid.shape
    """
    # Transform RHS to spectral space (DST-I)
    rhs_hat = dstn(rhs, type=1)

    # Solve in spectral space: v_hat = rhs_hat / eigenvalues
    v_hat = rhs_hat / eigenvalues

    # Back to physical space
    # DST-I inverse: the inverse of dstn(type=1) is (2/(N+1))^3 * dstn(type=1)
    # But scipy's idstn handles normalization automatically.
    v = idstn(v_hat, type=1)

    return v + bc_val


def solve_field_picard_fft(grid: TGPGrid, S: np.ndarray,
                           beta: float, gamma: float,
                           u_init: np.ndarray,
                           max_iter: int = 200,
                           tol: float = 1e-6,
                           alpha: float = 0.7,
                           verbose: bool = False) -> Tuple[np.ndarray, dict]:
    """Solve the field equation via linearization-based Picard iteration.

    The equation is: nabla^2 u = N(u) - 12*pi*S
    where N(u) = 3*gamma*u^{5/3} - 3*beta*u^{4/3}.

    At u=1: N'(1) = 5*gamma - 4*beta = mu > 0 (for vacuum).

    We subtract mu*u from both sides:
        (nabla^2 - mu) u = N(u) - mu*u - 12*pi*S

    The key property: d/du[N(u) - mu*u] = 0 at u=1, so the RHS has
    ZERO Jacobian at the fixed point => superlinear convergence.

    The operator (nabla^2 - mu) has eigenvalues (lambda_k - mu) < 0
    for all modes (since lambda_k < 0, mu > 0). Well-conditioned.

    Parameters
    ----------
    grid : TGPGrid
    S : ndarray
        Source field.
    beta, gamma : float
    u_init : ndarray
        Initial guess.
    max_iter : int
        Maximum Picard iterations.
    tol : float
        Convergence tolerance (L-infinity norm of change).
    alpha : float
        Under-relaxation parameter (0 < alpha <= 1).
    verbose : bool

    Returns
    -------
    u : ndarray
        Converged field.
    info : dict
        'converged': bool, 'iterations': int, 'residual_norm': float
    """
    # Shift: mu = N'(1) = 5*gamma - 4*beta
    mu = 5.0 * gamma - 4.0 * beta
    if mu <= 0:
        mu = 1.0

    eigenvalues = _dst_laplacian_eigenvalues(grid)
    # (nabla^2 - mu) has eigenvalues (lambda_k - mu), all strictly negative
    shifted_eig = eigenvalues - mu

    u = u_init.copy()

    for iteration in range(max_iter):
        # N(u) = 3*gamma*u^{5/3} - 3*beta*u^{4/3}  (CORRECT signs)
        u_safe = np.maximum(u, 1e-30)
        u13 = u_safe ** (1.0 / 3.0)
        N_u = 3.0 * gamma * u_safe * u13 * u13 - 3.0 * beta * u_safe * u13

        # (nabla^2 - mu) u_{n+1} = N(u_n) - mu*u_n - 12*pi*S
        # With v = u - 1 (homogeneous BC):
        #   (nabla^2 - mu)(v+1) = N(u_n) - mu*u_n - 12*pi*S
        #   nabla^2 v - mu*v - mu = N(u_n) - mu*u_n - 12*pi*S
        #   (nabla^2 - mu)v = N(u_n) - mu*(u_n - 1) - 12*pi*S
        rhs_v = N_u - mu * (u - 1.0) - 12.0 * np.pi * S

        rhs_hat = dstn(rhs_v, type=1)
        v_hat = rhs_hat / shifted_eig
        v_new = idstn(v_hat, type=1)
        u_new = v_new + 1.0

        # Clamp to positive
        u_new = np.maximum(u_new, 1e-30)

        # Under-relaxation
        u_next = alpha * u_new + (1.0 - alpha) * u

        # Convergence check
        change = np.max(np.abs(u_next - u))
        if verbose and (iteration % 10 == 0 or iteration < 5):
            res_norm = np.max(np.abs(residual_u(u_next, grid, S, beta, gamma)))
            print(f"  Picard iter {iteration:4d}: "
                  f"du_max = {change:.3e}, ||F||_inf = {res_norm:.3e}")

        u = u_next

        if change < tol:
            res_norm = np.max(np.abs(residual_u(u, grid, S, beta, gamma)))
            if verbose:
                print(f"  Picard converged at iter {iteration}: "
                      f"||F||_inf = {res_norm:.3e}")
            return u, {'converged': True, 'iterations': iteration + 1,
                       'residual_norm': float(res_norm)}

    res_norm = np.max(np.abs(residual_u(u, grid, S, beta, gamma)))
    if verbose:
        print(f"  Picard NOT converged after {max_iter} iters: "
              f"||F||_inf = {res_norm:.3e}")
    return u, {'converged': False, 'iterations': max_iter,
               'residual_norm': float(res_norm)}


# ── Newton-Krylov solver (robust fallback) ───────────────────────────────

def solve_field_newton_krylov(grid: TGPGrid, S: np.ndarray,
                              beta: float, gamma: float,
                              u_init: np.ndarray,
                              tol: float = 1e-6,
                              maxiter: int = 200,
                              verbose: bool = False) -> Tuple[np.ndarray, dict]:
    """Solve via scipy's Newton-Krylov (LGMRES inner loop).

    Finds u such that residual_u(u) = 0.

    The f_tol for scipy is the L2 norm of the full residual vector.
    For a grid of N^3 points, ||F||_2 ~ sqrt(N^3) * per_point_residual.
    We scale the tolerance accordingly.

    Parameters
    ----------
    grid : TGPGrid
    S : ndarray
    beta, gamma : float
    u_init : ndarray
    tol : float
        Target L-infinity residual per point.
    maxiter : int
    verbose : bool

    Returns
    -------
    u : ndarray
    info : dict
    """
    from scipy.optimize import newton_krylov, NoConvergence

    n_total = grid.Nx * grid.Ny * grid.Nz
    # Scale tolerance: L2 norm ~ sqrt(N) * L_inf, but the source region
    # has O(1) residual localized in ~(sigma/dx)^3 cells. Be generous.
    f_tol_l2 = tol * np.sqrt(n_total) * 0.1

    def F(u_flat):
        u_3d = u_flat.reshape(grid.shape)
        u_3d = np.maximum(u_3d, 1e-30)
        return residual_u(u_3d, grid, S, beta, gamma).ravel()

    u_sol = None
    try:
        if verbose:
            print(f"  Newton-Krylov starting (f_tol_l2 = {f_tol_l2:.2e})...")
        u_sol_flat = newton_krylov(
            F, u_init.ravel(),
            method='lgmres',
            f_tol=f_tol_l2,
            maxiter=maxiter,
            verbose=verbose,
        )
        u_sol = u_sol_flat.reshape(grid.shape)
    except NoConvergence as e:
        # Extract partial solution from the exception
        if verbose:
            print(f"  Newton-Krylov: maxiter reached, extracting partial solution")
        u_sol_flat = np.asarray(e.args[0])
        if u_sol_flat.size == n_total:
            u_sol = u_sol_flat.reshape(grid.shape)
    except Exception as e:
        if verbose:
            print(f"  Newton-Krylov failed: {e}")

    if u_sol is not None and np.all(np.isfinite(u_sol)):
        u_sol = np.maximum(u_sol, 1e-30)
        res_norm = np.max(np.abs(residual_u(u_sol, grid, S, beta, gamma)))
        converged = res_norm < tol * 100  # relaxed convergence check
        if verbose:
            print(f"  Newton-Krylov result: ||F||_inf = {res_norm:.3e}, "
                  f"converged = {converged}")
        return u_sol, {'converged': converged, 'iterations': -1,
                       'residual_norm': float(res_norm)}
    else:
        return u_init, {'converged': False, 'iterations': -1,
                        'residual_norm': float('inf'),
                        'error': 'NK produced invalid result'}


def solve_field_anderson(grid: TGPGrid, S: np.ndarray,
                         beta: float, gamma: float,
                         u_init: np.ndarray,
                         max_iter: int = 300,
                         tol: float = 1e-6,
                         m_anderson: int = 10,
                         alpha: float = 0.3,
                         verbose: bool = False) -> Tuple[np.ndarray, dict]:
    """Solve via Anderson-accelerated fixed-point iteration.

    Uses the shifted-Laplacian Picard step as the base iteration,
    with Anderson mixing to accelerate convergence and stabilize
    the near-source nonlinear region.

    The idea: compute the Picard update g(u) = (nabla^2 + mu)^{-1}[N(u) + mu*u - 3S],
    then apply Anderson acceleration to the fixed-point problem u = g(u).

    Parameters
    ----------
    grid : TGPGrid
    S : ndarray
    beta, gamma : float
    u_init : ndarray
    max_iter : int
    tol : float
    m_anderson : int
        Anderson mixing depth (number of previous iterates to use).
    alpha : float
        Damping parameter for the base iteration.
    verbose : bool

    Returns
    -------
    u : ndarray
    info : dict
    """
    # Shift parameter: dRHS/du at u=1 is (5*gamma - 4*beta).
    # To zero the Jacobian of the Picard map, use (nabla^2 - mu):
    #   (nabla^2 - mu) u_{n+1} = N(u_n) - mu*u_n - 12*pi*S
    # where N(u) = 3*gamma*u^{5/3} - 3*beta*u^{4/3} and
    # dN/du|_{u=1} = 5*gamma - 4*beta = mu.
    # d/du[N(u) - mu*u] = 0 at u=1 => superlinear convergence.
    mu = 5.0 * gamma - 4.0 * beta
    if mu <= 0:
        mu = 1.0

    eigenvalues = _dst_laplacian_eigenvalues(grid)
    # (nabla^2 - mu) has eigenvalues (lambda_k - mu), all strictly negative
    shifted_eig = eigenvalues - mu

    def picard_step(u):
        """One damped shifted-Picard step.

        Solve: (nabla^2 - mu)(v+1) = N(u_n) - mu*u_n - 12*pi*S
        where v = u - 1 (homogeneous BC).

        Expanding: nabla^2 v - mu*v - mu = N(u_n) - mu*u_n - 12*pi*S
        So: (nabla^2 - mu)v = N(u_n) - mu*u_n - 12*pi*S + mu
                             = N(u_n) - mu*(u_n - 1) - 12*pi*S
        """
        u_safe = np.maximum(u, 1e-30)
        u13 = u_safe ** (1.0 / 3.0)
        # N(u) = 3*gamma*u^{5/3} - 3*beta*u^{4/3}  (CORRECT signs from ODE)
        N_u = 3.0 * gamma * u_safe * u13 * u13 - 3.0 * beta * u_safe * u13
        rhs_v = N_u - mu * (u - 1.0) - 12.0 * np.pi * S
        rhs_hat = dstn(rhs_v, type=1)
        v_hat = rhs_hat / shifted_eig
        v_new = idstn(v_hat, type=1)
        u_new = np.maximum(v_new + 1.0, 1e-30)
        return alpha * u_new + (1.0 - alpha) * u

    u = u_init.copy()

    # Anderson acceleration storage
    # Store residuals r_k = g(u_k) - u_k and iterates
    F_hist = []  # residuals (flat)
    X_hist = []  # iterates (flat)

    for iteration in range(max_iter):
        g_u = picard_step(u)
        r = g_u - u  # residual of fixed-point iteration

        change = np.max(np.abs(r))

        if verbose and (iteration % 10 == 0 or iteration < 5):
            pde_res = np.max(np.abs(residual_u(g_u, grid, S, beta, gamma)))
            print(f"  Anderson iter {iteration:4d}: "
                  f"du_max = {change:.3e}, ||F||_inf = {pde_res:.3e}")

        if change < tol:
            u = g_u
            res_norm = np.max(np.abs(residual_u(u, grid, S, beta, gamma)))
            if verbose:
                print(f"  Anderson converged at iter {iteration}: "
                      f"||F||_inf = {res_norm:.3e}")
            return u, {'converged': True, 'iterations': iteration + 1,
                       'residual_norm': float(res_norm)}

        # Anderson mixing
        r_flat = r.ravel()
        u_flat = u.ravel()

        F_hist.append(r_flat.copy())
        X_hist.append(u_flat.copy())

        if len(F_hist) > m_anderson + 1:
            F_hist.pop(0)
            X_hist.pop(0)

        m_k = len(F_hist)
        if m_k >= 2:
            # Build the matrix of residual differences
            dF = np.column_stack([F_hist[i] - F_hist[-1]
                                  for i in range(m_k - 1)])
            # Solve least squares: min ||dF @ theta - F_hist[-1]||
            theta, _, _, _ = np.linalg.lstsq(dF, -F_hist[-1], rcond=None)

            # Anderson update
            u_new = (1.0 - np.sum(theta)) * (X_hist[-1] + F_hist[-1])
            for i in range(m_k - 1):
                u_new += theta[i] * (X_hist[i] + F_hist[i])

            u = np.maximum(u_new.reshape(grid.shape), 1e-30)
        else:
            u = g_u

    res_norm = np.max(np.abs(residual_u(u, grid, S, beta, gamma)))
    if verbose:
        print(f"  Anderson NOT converged after {max_iter} iters: "
              f"||F||_inf = {res_norm:.3e}")
    return u, {'converged': False, 'iterations': max_iter,
               'residual_norm': float(res_norm)}


# ── Main solver (two-tier) ───────────────────────────────────────────────

def solve_field(grid: TGPGrid,
                positions: np.ndarray,
                C_values: np.ndarray,
                beta: float = 1.0,
                gamma: float = 1.0,
                sigma: Optional[float] = None,
                picard_kw: Optional[dict] = None,
                nk_fallback: bool = True,
                verbose: bool = False) -> Tuple[np.ndarray, dict]:
    """Solve the full nonlinear TGP field equation on a 3D grid.

    Strategy: try Picard FFT first; if it fails, fall back to Newton-Krylov.

    Parameters
    ----------
    grid : TGPGrid
    positions : array, shape (N, 3)
        Source positions.
    C_values : array, shape (N,)
        Source strengths.
    beta, gamma : float
        TGP couplings (default vacuum: β = γ = 1).
    sigma : float or None
        Source smoothing width. If None, uses 1.5 * dx.
    picard_kw : dict or None
        Keyword arguments for solve_field_picard_fft.
    nk_fallback : bool
        Whether to try Newton-Krylov if Picard fails.
    verbose : bool

    Returns
    -------
    u : ndarray, shape grid.shape
        Converged u = g³.
    info : dict
        Solver information.
    """
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)

    if positions.ndim == 1:
        positions = positions[np.newaxis, :]
    if C_values.ndim == 0:
        C_values = C_values[np.newaxis]

    if sigma is None:
        sigma = 1.5 * grid.dx

    # Build source field
    S = build_source_field(grid, positions, C_values, sigma)

    # Check source normalization
    for i, Ci in enumerate(C_values):
        Si = gaussian_source(grid, positions[i], Ci, sigma)
        integral = np.sum(Si) * grid.dV
        if verbose:
            print(f"  Source {i}: C = {Ci:.4f}, int S d3x = {integral:.6f}, "
                  f"error = {abs(integral - Ci) / max(abs(Ci), 1e-30):.2e}")

    # Initial guess
    u_init = initial_guess_u(grid, positions, C_values, beta, gamma)

    # Check boundary values
    bnd_err = max(
        np.max(np.abs(u_init[0, :, :] - 1)),
        np.max(np.abs(u_init[-1, :, :] - 1)),
        np.max(np.abs(u_init[:, 0, :] - 1)),
        np.max(np.abs(u_init[:, -1, :] - 1)),
        np.max(np.abs(u_init[:, :, 0] - 1)),
        np.max(np.abs(u_init[:, :, -1] - 1)),
    )
    if verbose:
        print(f"  Boundary check: max |u_init - 1| at boundary = {bnd_err:.3e}")
    if bnd_err > 0.1:
        warnings.warn(
            f"Domain too small: |g-1| = {bnd_err:.3e} at boundaries (want < 0.01). "
            f"Increase L or decrease C.",
            stacklevel=2
        )

    # --- Anderson-accelerated Picard (primary) ---
    pkw = dict(max_iter=300, tol=1e-6, alpha=0.3, m_anderson=10,
               verbose=verbose)
    if picard_kw:
        pkw.update(picard_kw)

    u, info = solve_field_anderson(grid, S, beta, gamma, u_init, **pkw)
    info['solver'] = 'anderson'

    if info['converged']:
        return u, info

    # --- Newton-Krylov fallback ---
    if nk_fallback:
        if verbose:
            print("  Anderson did not converge; trying Newton-Krylov...")
        # Use the CLEAN initial guess (not the possibly-diverged Anderson result)
        u_start = u if np.all(np.isfinite(u)) else u_init
        u_nk, info_nk = solve_field_newton_krylov(
            grid, S, beta, gamma, u_start, verbose=verbose
        )
        info_nk['solver'] = 'newton_krylov'
        if info_nk['converged']:
            return u_nk, info_nk
        # Return best available
        if info_nk['residual_norm'] < info['residual_norm']:
            return u_nk, info_nk

    return u, info


# ── Energy extraction ────────────────────────────────────────────────────

def total_field_energy(u: np.ndarray, grid: TGPGrid,
                       beta: float, gamma: float) -> float:
    """Total field energy from the converged u = g^3 solution.

    E = int [1/2 g^2 |grad g|^2 + (beta/3) g^3 - (gamma/4) g^4 - e_vac] d^3x

    where K(g) = g^2 (full kinetic mode), e_vac = beta/3 - gamma/4.

    This is the simplified N-body energy functional from the TGP paper.
    The Feynman V2 and V3 references are derived from perturbative
    expansion of THIS functional, so it should be used for comparisons.

    Note: the field equation (ODE) comes from a different action (full
    gravitational). At linear order both give identical Yukawa profiles,
    so evaluating E[g] on the ODE solution is accurate for V2 extraction.
    V3 extraction requires higher grid resolution (see ex216 notes).

    Parameters
    ----------
    u : ndarray, shape grid.shape
        Converged field u = g^3.
    grid : TGPGrid
    beta, gamma : float

    Returns
    -------
    E : float
        Total field energy (vacuum-subtracted).
    """
    u_safe = np.maximum(u, 1e-30)
    g = u_safe ** (1.0 / 3.0)

    # |grad g|^2 via central differences (BC: g = 1)
    grad_g_sq = gradient_sq_3d(g, grid.dx, grid.dy, grid.dz, bc_val=1.0)

    # Kinetic coupling K = g^2 (full mode)
    K = g ** 2

    # Vacuum energy density (constant to subtract)
    e_vac = beta / 3.0 - gamma / 4.0

    # Energy density
    e = 0.5 * K * grad_g_sq + (beta / 3.0) * g**3 - (gamma / 4.0) * g**4 - e_vac

    return float(np.sum(e) * grid.dV)


def total_field_energy_ode(u: np.ndarray, grid: TGPGrid,
                           beta: float, gamma: float,
                           S: Optional[np.ndarray] = None) -> float:
    """Total field energy using ODE-consistent functional.

    E[u] = int [1/2 |grad u|^2 + W(u) - 12*pi*S*u - W_vac] d^3x

    where W(u) = (9*gamma/8)*u^{8/3} - (9*beta/7)*u^{7/3}
    and W_vac = W(1) = 9*gamma/8 - 9*beta/7.

    This is the natural energy for the field equation:
      laplacian(u) = 3*gamma*u^{5/3} - 3*beta*u^{4/3} - 12*pi*S

    Parameters
    ----------
    u : ndarray, shape grid.shape
    grid : TGPGrid
    beta, gamma : float
    S : ndarray or None
        Source density field for coupling term.

    Returns
    -------
    E : float
    """
    u_safe = np.maximum(u, 1e-30)

    grad_u_sq = gradient_sq_3d(u_safe, grid.dx, grid.dy, grid.dz, bc_val=1.0)

    u13 = u_safe ** (1.0 / 3.0)
    u73 = u_safe ** 2 * u13
    u83 = u73 * u13

    W_vac = 9.0 * gamma / 8.0 - 9.0 * beta / 7.0

    e = (0.5 * grad_u_sq
         + (9.0 * gamma / 8.0) * u83
         - (9.0 * beta / 7.0) * u73
         - W_vac)

    if S is not None:
        e -= 12.0 * np.pi * S * u_safe

    return float(np.sum(e) * grid.dV)


# ── Convenience: single-source and interaction energies ──────────────────

_single_energy_cache: dict = {}


def single_source_energy(C: float, grid: TGPGrid,
                         beta: float = 1.0, gamma: float = 1.0,
                         sigma: Optional[float] = None,
                         verbose: bool = False,
                         **solver_kw) -> float:
    """Solve for a single source at grid center and return its energy.

    Results are cached by (C, grid.Nx, grid.Lx, beta, gamma, sigma).
    """
    if sigma is None:
        sigma = 1.5 * grid.dx

    key = (C, grid.Nx, grid.Ny, grid.Nz, grid.Lx, grid.Ly, grid.Lz,
           beta, gamma, sigma)
    if key in _single_energy_cache:
        return _single_energy_cache[key]

    pos = np.array([[0.0, 0.0, 0.0]])
    Cv = np.array([C])
    u, info = solve_field(grid, pos, Cv, beta, gamma, sigma=sigma,
                          verbose=verbose, **solver_kw)
    if not info['converged']:
        warnings.warn(f"Single-source solve did not converge: {info}",
                      stacklevel=2)

    E = total_field_energy(u, grid, beta, gamma)
    _single_energy_cache[key] = E
    return E


def interaction_energy(u: np.ndarray, grid: TGPGrid,
                       C_values: np.ndarray,
                       beta: float = 1.0, gamma: float = 1.0,
                       sigma: Optional[float] = None,
                       verbose: bool = False,
                       **solver_kw) -> float:
    """Interaction energy: E_total - Σ E_single.

    Parameters
    ----------
    u : ndarray
        Converged multi-source field u = g³.
    grid : TGPGrid
    C_values : array
        Source strengths.
    beta, gamma : float
    sigma : float or None
    verbose : bool

    Returns
    -------
    E_int : float
        Interaction energy.
    """
    E_total = total_field_energy(u, grid, beta, gamma)

    E_singles = 0.0
    for C in C_values:
        E_singles += single_source_energy(C, grid, beta, gamma,
                                          sigma=sigma, verbose=verbose,
                                          **solver_kw)

    return E_total - E_singles


def extract_V2_V3(u: np.ndarray, grid: TGPGrid,
                  C_values: np.ndarray,
                  positions: np.ndarray,
                  beta: float = 1.0, gamma: float = 1.0,
                  sigma: Optional[float] = None,
                  verbose: bool = False,
                  **solver_kw) -> dict:
    """Extract interaction energy decomposition from a 3-body PDE solution.

    Computes E_int = E_total - 3*E_single and compares with the
    perturbative prediction 3*V2_Born + V3_Feynman.

    NOTE: Isolating V3 as a PDE residual (E_int - 3*V2) is limited
    by grid resolution -- V3/V2 ~ O(C) ~ 0.2% at typical parameters,
    which is below V2 numerical accuracy at moderate N.

    Parameters
    ----------
    u : ndarray
        Converged 3-source field.
    grid : TGPGrid
    C_values : array, shape (3,)
    positions : array, shape (3, 3)
    beta, gamma : float

    Returns
    -------
    result : dict with keys:
        'E_total', 'E_singles', 'E_int',
        'E_int_pert' (perturbative prediction),
        'E_int_ratio' (PDE/pert),
        'V2_born_total', 'V3_feynman' (if available)
    """
    from .tgp_field import screening_mass

    positions = np.asarray(positions)
    C_values = np.asarray(C_values)
    assert len(C_values) == 3, "extract_V2_V3 requires exactly 3 sources"

    E_total = total_field_energy(u, grid, beta, gamma)

    # Single-source energies
    E_singles = 0.0
    for C in C_values:
        E_singles += single_source_energy(C, grid, beta, gamma,
                                          sigma=sigma, verbose=verbose,
                                          **solver_kw)

    E_int = E_total - E_singles

    # Screened Born V2 for each pair
    m = screening_mass(beta, gamma)
    pairs = [(0, 1), (0, 2), (1, 2)]
    V2_born_total = 0.0
    for i, j in pairs:
        d = np.linalg.norm(positions[i] - positions[j])
        Ci, Cj = C_values[i], C_values[j]
        exp_md = np.exp(-m * d)
        grad_ov = 4.0 * np.pi * Ci * Cj * exp_md / d
        mass_ov = (1.0 + m**2) * 2.0 * np.pi * Ci * Cj * exp_md / m
        V2_born_total += grad_ov - mass_ov

    result = {
        'E_total': float(E_total),
        'E_singles': float(E_singles),
        'E_int': float(E_int),
        'V2_born_total': float(V2_born_total),
    }

    # Feynman V3 if available
    try:
        from .three_body_force_exact import three_body_energy_exact
        d12 = np.linalg.norm(positions[0] - positions[1])
        d13 = np.linalg.norm(positions[0] - positions[2])
        d23 = np.linalg.norm(positions[1] - positions[2])
        V3_feynman = three_body_energy_exact(
            d12, d13, d23,
            C_values[0], C_values[1], C_values[2],
            beta, gamma
        )
        result['V3_feynman'] = float(V3_feynman)
        E_int_pert = V2_born_total + V3_feynman
    except ImportError:
        E_int_pert = V2_born_total

    result['E_int_pert'] = float(E_int_pert)
    if abs(E_int_pert) > 1e-30:
        result['E_int_ratio'] = float(E_int / E_int_pert)

    return result


# ── Radial profile extraction ────────────────────────────────────────────

def extract_radial_profile(u: np.ndarray, grid: TGPGrid,
                           center: np.ndarray = None,
                           n_bins: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Extract azimuthally averaged radial profile g(r) from 3D u field.

    Parameters
    ----------
    u : ndarray
        Field u = g³.
    grid : TGPGrid
    center : array, shape (3,) or None
        Center point (default: grid center = origin).
    n_bins : int
        Number of radial bins.

    Returns
    -------
    r_bins : ndarray, shape (n_bins,)
        Bin centers.
    g_bins : ndarray, shape (n_bins,)
        Averaged g = u^{1/3} in each bin.
    """
    if center is None:
        center = np.array([0.0, 0.0, 0.0])

    g = np.maximum(u, 1e-30) ** (1.0 / 3.0)
    r = np.sqrt(
        (grid.X - center[0])**2 +
        (grid.Y - center[1])**2 +
        (grid.Z - center[2])**2
    )

    r_max = min(grid.Lx, grid.Ly, grid.Lz) / 2.0
    r_edges = np.linspace(0, r_max, n_bins + 1)
    r_bins = 0.5 * (r_edges[:-1] + r_edges[1:])
    g_bins = np.zeros(n_bins)

    for b in range(n_bins):
        mask = (r >= r_edges[b]) & (r < r_edges[b + 1])
        if np.any(mask):
            g_bins[b] = np.mean(g[mask])
        else:
            g_bins[b] = 1.0  # vacuum default

    return r_bins, g_bins
