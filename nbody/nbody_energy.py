"""
nbody_energy.py -- Total n-body energy assembler for TGP
=========================================================
Combines pairwise (2-body) and irreducible 3-body contributions
to compute the total interaction energy, forces, and Hessian matrix
for an arbitrary configuration of N TGP sources.

Energy decomposition
--------------------
    V_total = sum_{i<j} V_2(d_ij) + sum_{i<j<k} V_3(d_ij, d_ik, d_jk)

where V_2 is the pairwise potential from pairwise.py and V_3 is the
irreducible 3-body energy from three_body_terms.py.

Consistency check (CRITICAL)
-----------------------------
For N=2 bodies, V_3 = 0 identically (no triplets exist), so
total_energy must match pairwise.V_eff_total exactly.

Forces and Hessian
-------------------
Computed by finite differences of total_energy.  This is robust
and allows easy inclusion of higher-order multi-body terms in
the future without rederiving analytical gradients.
"""

import numpy as np
from itertools import combinations
from .pairwise import V_eff_total
from .three_body_terms import three_body_energy
from .tgp_field import default_beta_gamma


def total_energy(positions, C_values, beta=None, gamma=None):
    """Total n-body interaction energy.

    V = sum_{i<j} V_2(d_ij) + sum_{i<j<k} V_3(i, j, k)

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Body positions in 3D.
    C_values : array_like, shape (N,)
        Dimensionless source strengths.
    beta, gamma : float or None
        Self-interaction couplings (default: vacuum beta=gamma=1).

    Returns
    -------
    V : float
        Total interaction energy.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)

    V = 0.0

    # Pairwise (2-body) terms
    for i, j in combinations(range(N), 2):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        V += V_eff_total(d_ij, C_values[i], C_values[j], beta, gamma)

    # Irreducible 3-body terms
    for i, j, k in combinations(range(N), 3):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        d_ik = np.linalg.norm(positions[i] - positions[k])
        d_jk = np.linalg.norm(positions[j] - positions[k])
        V += three_body_energy(d_ij, d_ik, d_jk,
                               C_values[i], C_values[j], C_values[k],
                               beta, gamma)

    return V


def forces(positions, C_values, beta=None, gamma=None, eps=1e-6):
    """Forces on all bodies by finite differences.

    F_i = -grad_i V_total

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Body positions.
    C_values : array_like, shape (N,)
        Source strengths.
    beta, gamma : float or None
        Couplings (default: vacuum).
    eps : float
        Finite difference step size.

    Returns
    -------
    F : ndarray, shape (N, 3)
        Force vector on each body.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)

    F = np.zeros((N, 3))

    for i in range(N):
        for dim in range(3):
            pos_plus = positions.copy()
            pos_minus = positions.copy()
            pos_plus[i, dim] += eps
            pos_minus[i, dim] -= eps

            V_plus = total_energy(pos_plus, C_values, beta, gamma)
            V_minus = total_energy(pos_minus, C_values, beta, gamma)

            F[i, dim] = -(V_plus - V_minus) / (2.0 * eps)

    return F


def hessian(positions, C_values, beta=None, gamma=None, eps=1e-5):
    """Full 3N x 3N Hessian matrix of the total energy.

    H_{ia,jb} = d^2 V / (d x_{i,a} d x_{j,b})

    where i,j are body indices and a,b are spatial dimensions.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Body positions.
    C_values : array_like, shape (N,)
        Source strengths.
    beta, gamma : float or None
        Couplings (default: vacuum).
    eps : float
        Finite difference step size.

    Returns
    -------
    H : ndarray, shape (3*N, 3*N)
        Hessian matrix.  Index ordering: [x1,y1,z1, x2,y2,z2, ...].
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)
    dim_total = 3 * N
    H = np.zeros((dim_total, dim_total))

    V0 = total_energy(positions, C_values, beta, gamma)

    for idx_a in range(dim_total):
        i_a, d_a = divmod(idx_a, 3)  # body index, spatial dim
        for idx_b in range(idx_a, dim_total):
            i_b, d_b = divmod(idx_b, 3)

            if idx_a == idx_b:
                # Diagonal: (V(+h) + V(-h) - 2*V0) / h^2
                pos_p = positions.copy()
                pos_m = positions.copy()
                pos_p[i_a, d_a] += eps
                pos_m[i_a, d_a] -= eps
                Vp = total_energy(pos_p, C_values, beta, gamma)
                Vm = total_energy(pos_m, C_values, beta, gamma)
                H[idx_a, idx_a] = (Vp + Vm - 2.0 * V0) / eps**2
            else:
                # Off-diagonal: central difference
                # (V(++) + V(--) - V(+-) - V(-+)) / (4*h^2)
                pos_pp = positions.copy()
                pos_mm = positions.copy()
                pos_pm = positions.copy()
                pos_mp = positions.copy()
                pos_pp[i_a, d_a] += eps; pos_pp[i_b, d_b] += eps
                pos_mm[i_a, d_a] -= eps; pos_mm[i_b, d_b] -= eps
                pos_pm[i_a, d_a] += eps; pos_pm[i_b, d_b] -= eps
                pos_mp[i_a, d_a] -= eps; pos_mp[i_b, d_b] += eps

                Vpp = total_energy(pos_pp, C_values, beta, gamma)
                Vmm = total_energy(pos_mm, C_values, beta, gamma)
                Vpm = total_energy(pos_pm, C_values, beta, gamma)
                Vmp = total_energy(pos_mp, C_values, beta, gamma)

                H[idx_a, idx_b] = (Vpp + Vmm - Vpm - Vmp) / (4.0 * eps**2)
                H[idx_b, idx_a] = H[idx_a, idx_b]  # symmetric

    return H


# ── Decomposed energy (for diagnostics) ──────────────────────────────────

def energy_decomposition(positions, C_values, beta=None, gamma=None):
    """Return pairwise and 3-body contributions separately.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
    C_values : array_like, shape (N,)
    beta, gamma : float or None

    Returns
    -------
    V2_total : float
        Sum of all pairwise interactions.
    V3_total : float
        Sum of all irreducible 3-body interactions.
    V_total : float
        V2_total + V3_total.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)

    V2 = 0.0
    for i, j in combinations(range(N), 2):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        V2 += V_eff_total(d_ij, C_values[i], C_values[j], beta, gamma)

    V3 = 0.0
    for i, j, k in combinations(range(N), 3):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        d_ik = np.linalg.norm(positions[i] - positions[k])
        d_jk = np.linalg.norm(positions[j] - positions[k])
        V3 += three_body_energy(d_ij, d_ik, d_jk,
                                C_values[i], C_values[j], C_values[k],
                                beta, gamma)

    return V2, V3, V2 + V3
