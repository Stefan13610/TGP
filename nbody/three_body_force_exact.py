"""
three_body_force_exact.py — Exact 3-body force from differentiation under
the Feynman integral for the triple Yukawa overlap
==========================================================================

EXACT RESULT (no saddle-point, no approximation):
--------------------------------------------------
From the Fourier/Feynman representation (tgp_yukawa_fourier_feynman_reduction.tex).
Multipole / Legendre series route (semianalytic truncation): tgp_yukawa_IY_multipole_gegenbauer.tex
(\input{} from tgp_yukawa_exact_reduction.tex).

    I_Y(d12, d13, d23; m) = 2 * integral_{Delta_2} Delta^{-3/2}
                             * K_0( m * sqrt(Q / Delta) ) d_alpha

where:
    Q     = alpha_2 * d12^2 + alpha_1 * d13^2 + alpha_3 * d23^2
    Delta = alpha_1*alpha_2 + alpha_1*alpha_3 + alpha_2*alpha_3
    alpha_3 = 1 - alpha_1 - alpha_2

Differentiating under the integral (K_0' = -K_1):

    d I_Y / d d12 = -2*m * d12
                    * integral_{Delta_2} alpha_2 * Delta^{-2} / sqrt(Q)
                                         * K_1(m * sqrt(Q/Delta)) d_alpha

Analogously for d13 (alpha_1 factor) and d23 (alpha_3 factor).

Cartesian force on body 1:
    dI_Y / dx_1 = (dI_Y/dd12) * (x1-x2)/d12
                + (dI_Y/dd13) * (x1-x3)/d13

The 3-body energy for a triplet (1,2,3) is:
    V_3 = -6 * gamma * C1 * C2 * C3 * I_Y

so the force on body 1 is:
    F_1 = -dV_3/dx_1 = 6 * gamma * C1 * C2 * C3 * dI_Y/dx_1

ADVANTAGES over saddle-point approximation:
  - No error from saddle-point: exact to machine precision
  - Correct at all separations (large and small m*d)
  - Correct geometry (not just semi-perimeter dependent)

COST:
  - Each force evaluation requires a 2D numerical integral
    over the standard 2-simplex (Alpha1 in [0,1], Alpha2 in [0,1-Alpha1])
  - With n_quad=50 quadrature points per dimension: ~1250 function evals per
    (i,j,k) triplet per Cartesian direction.  For N bodies: O(N^3) triplets.
  - Use n_quad=30 for dynamics, n_quad=100 for energy/force benchmarking.

COMPARISON WITH NUMERICAL 3D INTEGRATION:
  - triple_overlap_numerical() in three_body_terms.py uses N^3 grid (N=80
    default → 512000 points).
  - This module uses a 2D simplex integral (~O(n_quad^2) ~ 2500 points).
  - This module is BOTH more accurate AND faster by a factor of ~200.

JACOBIAN (exact ``V_3``):
  - ``_d2I_dd_matrix``, ``three_body_force_jacobian_triplet_exact``,
    ``three_body_force_jacobian_exact`` — ``∂F/∂x`` via distance Hessian of
    ``I_Y`` + chain rule (same quadrature grid as forces).
"""

from __future__ import annotations

import numpy as np
from scipy.special import k0, k1
from scipy.integrate import dblquad
from itertools import combinations

from .tgp_field import default_beta_gamma, screening_mass


# ── Core integrals ───────────────────────────────────────────────────────


def _build_quadrature_points(n=40):
    """
    Build a fixed quadrature grid over the standard 2-simplex
    {(a1, a2) : a1>=0, a2>=0, a1+a2<=1}.

    Uses a product Gauss-Legendre grid mapped to the simplex via
    Duffy transformation:  a1 = u,  a2 = v*(1-u).
    Jacobian = (1 - u).

    Parameters
    ----------
    n : int
        Number of quadrature points per dimension (n^2 total).

    Returns
    -------
    a1, a2, a3, weights : ndarray, shape (n^2,)
        Barycentric coordinates and integration weights (including Jacobian).
    """
    pts, wts = np.polynomial.legendre.leggauss(n)
    # Map Gauss-Legendre points from [-1,1] to [0,1]
    u = 0.5 * (pts + 1.0)
    w = 0.5 * wts

    U, V = np.meshgrid(u, u, indexing='ij')
    Wu, Wv = np.meshgrid(w, w, indexing='ij')

    U = U.ravel()
    V = V.ravel()
    Wu = Wu.ravel()
    Wv = Wv.ravel()

    a1 = U
    a2 = V * (1.0 - U)
    a3 = 1.0 - a1 - a2          # guaranteed >= 0 on the simplex

    # Jacobian of Duffy transform: |d(a1,a2)/d(u,v)| = (1-u)
    J = 1.0 - U
    weights = Wu * Wv * J        # shape (n^2,)

    return a1, a2, a3, weights


# Pre-built quadrature grids at different resolutions
_QUAD_CACHE: dict = {}


def _get_quad(n):
    if n not in _QUAD_CACHE:
        _QUAD_CACHE[n] = _build_quadrature_points(n)
    return _QUAD_CACHE[n]


def yukawa_overlap_exact(d12, d13, d23, m, n_quad=40, eps=1e-9):
    """
    Exact triple Yukawa overlap integral via Feynman parametrization.

        I_Y = 2 * int_{Delta_2} Delta^{-3/2} * K_0(m*sqrt(Q/Delta)) d_alpha

    Parameters
    ----------
    d12, d13, d23 : float  (must be > 0)
    m             : float  screening mass (> 0)
    n_quad        : int    quadrature resolution
    eps           : float  small guard against Delta -> 0 singularity

    Returns
    -------
    I_Y : float
    """
    d12, d13, d23 = float(d12), float(d13), float(d23)
    m = float(m)

    a1, a2, a3, weights = _get_quad(n_quad)

    Q     = a2 * d12**2 + a1 * d13**2 + a3 * d23**2
    Delta = a1*a2 + a1*a3 + a2*a3
    Delta = np.maximum(Delta, eps)

    u = m * np.sqrt(Q / Delta)

    # K_0 is integrable (log singularity) near u=0, but Q>0 guarantees u>0
    # for interior simplex points.
    integrand = Delta**(-1.5) * k0(u)
    integrand = np.where(np.isfinite(integrand), integrand, 0.0)

    return 2.0 * np.dot(weights, integrand)


def _dI_dd_components(d12, d13, d23, m, n_quad=40, eps=1e-9):
    """
    Compute the three scalar derivatives dI_Y/dd_ij simultaneously
    using a single quadrature pass.

    Returns (dI_dd12, dI_dd13, dI_dd23).

    Formula (derived in module docstring):
        dI_Y/dd12 = -2*m*d12 * int alpha_2 * Delta^{-2} / sqrt(Q) * K_1(u) d_alpha
        dI_Y/dd13 = -2*m*d13 * int alpha_1 * Delta^{-2} / sqrt(Q) * K_1(u) d_alpha
        dI_Y/dd23 = -2*m*d23 * int alpha_3 * Delta^{-2} / sqrt(Q) * K_1(u) d_alpha

    where u = m*sqrt(Q/Delta).
    """
    a1, a2, a3, weights = _get_quad(n_quad)

    Q     = a2 * d12**2 + a1 * d13**2 + a3 * d23**2
    Delta = a1*a2 + a1*a3 + a2*a3
    Delta = np.maximum(Delta, eps)
    sqrtQ = np.sqrt(np.maximum(Q, eps**2))

    u = m * np.sqrt(Q / Delta)

    # Kernel common to all three derivatives:
    # K_1(u) * Delta^{-2} / sqrt(Q)
    K1u = k1(u)
    kernel = K1u * Delta**(-2.0) / sqrtQ
    kernel = np.where(np.isfinite(kernel), kernel, 0.0)

    J12 = np.dot(weights, a2 * kernel)
    J13 = np.dot(weights, a1 * kernel)
    J23 = np.dot(weights, a3 * kernel)

    dI_dd12 = -2.0 * m * d12 * J12
    dI_dd13 = -2.0 * m * d13 * J13
    dI_dd23 = -2.0 * m * d23 * J23

    return dI_dd12, dI_dd13, dI_dd23


def _d2I_dd_matrix(d12, d13, d23, m, n_quad=40, eps=1e-9):
    """
    Symmetric 3x3 Hessian of ``I_Y`` in pairwise distances
    (indices 0=d12, 1=d13, 2=d23), one quadrature pass.

    Uses the same integrand ``K = Delta^{-2} Q^{-1/2} K_1(u)`` as
    ``_dI_dd_components``, with ``∂K/∂d_ℓ`` from the chain rule on ``Q`` and
    ``u = m\\sqrt{Q/Delta}``.
    """
    d12, d13, d23 = float(d12), float(d13), float(d23)
    m = float(m)

    a1, a2, a3, weights = _get_quad(n_quad)
    w = np.asarray(weights, dtype=float)

    Q = a2 * d12**2 + a1 * d13**2 + a3 * d23**2
    Delta = a1 * a2 + a1 * a3 + a2 * a3
    Delta = np.maximum(Delta, eps)
    sqrtQ = np.sqrt(np.maximum(Q, eps**2))
    u = m * np.sqrt(Q / Delta)

    K1u = k1(u)
    K0u = k0(u)
    u_safe = np.maximum(u, eps)
    K1p = -K0u - K1u / u_safe
    K1p = np.where(np.isfinite(K1p), K1p, 0.0)

    T = Delta ** (-2.0)
    S = Q ** (-0.5)
    K = T * S * K1u
    K = np.where(np.isfinite(K), K, 0.0)

    dQ = np.stack(
        [2.0 * a2 * d12, 2.0 * a1 * d13, 2.0 * a3 * d23], axis=0
    )  # shape (3, npt)

    du_d = (m * m) * dQ / (2.0 * np.maximum(u, eps) * Delta)

    dS_d = -0.5 * (Q ** (-1.5)) * dQ
    dK_d = T * (dS_d * K1u + S * K1p * du_d)  # shape (3, npt)
    dK_d = np.where(np.isfinite(dK_d), dK_d, 0.0)

    alpha_d = np.stack([a2, a1, a3], axis=0)
    dvals = np.array([d12, d13, d23], dtype=float)

    H = np.zeros((3, 3), dtype=float)
    for k in range(3):
        for l in range(3):
            acc = -2.0 * m * dvals[k] * np.dot(w, alpha_d[k] * dK_d[l])
            if k == l:
                acc += -2.0 * m * np.dot(w, alpha_d[k] * K)
            H[k, l] = acc
    return H


def _dd2_single_distance_9x9(e: np.ndarray, d: float, i0: int, j0: int) -> np.ndarray:
    """``∂²d/∂coord²`` for ``d = ||p_i - p_j||``, embedded in 9x9 (3 bodies × 3)."""
    B = (np.eye(3) - np.outer(e, e)) / d
    H = np.zeros((9, 9), dtype=float)
    sl_i = slice(i0, i0 + 3)
    sl_j = slice(j0, j0 + 3)
    H[sl_i, sl_i] = B
    H[sl_j, sl_j] = B
    H[sl_i, sl_j] = -B
    H[sl_j, sl_i] = -B
    return H


def _triplet_I_cartesian_hessian(
    pos1: np.ndarray,
    pos2: np.ndarray,
    pos3: np.ndarray,
    m: float,
    n_quad: int,
    eps: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return ``(h_vec, H9)`` with ``h_vec[k]=∂I/∂d_k`` (length 3) and
    ``H9 = ∂²I/∂X²`` for stacked ``X=(x1,y1,z1,x2,...,x3)``.
    """
    p1 = np.asarray(pos1, dtype=float).ravel()
    p2 = np.asarray(pos2, dtype=float).ravel()
    p3 = np.asarray(pos3, dtype=float).ravel()

    r12 = p1 - p2
    r13 = p1 - p3
    r23 = p2 - p3
    d12 = float(np.linalg.norm(r12))
    d13 = float(np.linalg.norm(r13))
    d23 = float(np.linalg.norm(r23))

    if d12 < 1e-15 or d13 < 1e-15 or d23 < 1e-15:
        return np.zeros(3), np.zeros((9, 9))

    e12 = r12 / d12
    e13 = r13 / d13
    e23 = r23 / d23

    dI12, dI13, dI23 = _dI_dd_components(d12, d13, d23, m, n_quad=n_quad, eps=eps)
    h = np.array([dI12, dI13, dI23], dtype=float)
    H_d = _d2I_dd_matrix(d12, d13, d23, m, n_quad=n_quad, eps=eps)

    A = np.zeros((3, 9), dtype=float)
    A[0, 0:3] = e12
    A[0, 3:6] = -e12
    A[1, 0:3] = e13
    A[1, 6:9] = -e13
    A[2, 3:6] = e23
    A[2, 6:9] = -e23

    H9 = A.T @ H_d @ A
    H9 += h[0] * _dd2_single_distance_9x9(e12, d12, 0, 3)
    H9 += h[1] * _dd2_single_distance_9x9(e13, d13, 0, 6)
    H9 += h[2] * _dd2_single_distance_9x9(e23, d23, 3, 6)

    return h, H9


def three_body_force_jacobian_triplet_exact(
    pos1,
    pos2,
    pos3,
    C1,
    C2,
    C3,
    beta=None,
    gamma=None,
    n_quad=40,
    eps=1e-9,
) -> np.ndarray:
    """
    ``∂(F1,F2,F3)/∂(x1,x2,x3)`` flattened (9×9) for one triplet's exact Yukawa
    ``V_3`` forces (same sign convention as ``three_body_force_triplet_exact``).
    """
    beta, gamma_c = default_beta_gamma(beta, gamma)
    m = screening_mass(beta, gamma_c)
    coeff = 6.0 * gamma_c * C1 * C2 * C3
    _h, H_I = _triplet_I_cartesian_hessian(pos1, pos2, pos3, m, n_quad, eps=eps)
    return coeff * H_I


def three_body_force_jacobian_exact(
    positions: np.ndarray,
    C_values: np.ndarray,
    beta=None,
    gamma=None,
    n_quad=40,
    eps=1e-9,
) -> np.ndarray:
    """
    Global ``∂F/∂x`` (shape ``(3N, 3N)``) from summing exact triplet Yukawa
    ``V_3`` forces over ``combinations(N,3)``.
    """
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)
    m = 3 * N
    dFdq = np.zeros((m, m), dtype=float)

    if N < 3:
        return dFdq

    for i, j, k in combinations(range(N), 3):
        J9 = three_body_force_jacobian_triplet_exact(
            positions[i],
            positions[j],
            positions[k],
            C_values[i],
            C_values[j],
            C_values[k],
            beta=beta,
            gamma=gamma,
            n_quad=n_quad,
            eps=eps,
        )
        ix = [3 * i + a for a in range(3)]
        iy = [3 * j + a for a in range(3)]
        iz = [3 * k + a for a in range(3)]
        idx = ix + iy + iz
        for a in range(9):
            for b in range(9):
                dFdq[idx[a], idx[b]] += J9[a, b]
    return dFdq


# ── 3-body energy (exact) ────────────────────────────────────────────────


def three_body_energy_exact(d12, d13, d23, C1, C2, C3,
                             beta=None, gamma=None, n_quad=40):
    """
    Exact irreducible 3-body energy for a triplet.

        V_3 = -6 * gamma * C1*C2*C3 * I_Y(d12, d13, d23; m)

    Uses the Feynman 2D integral (no saddle-point approximation).

    Parameters
    ----------
    d12, d13, d23 : float   pairwise distances
    C1, C2, C3   : float   source strengths
    beta, gamma  : float   TGP couplings (default: vacuum beta=gamma=1)
    n_quad       : int     quadrature resolution

    Returns
    -------
    V3 : float
    """
    beta, gamma_c = default_beta_gamma(beta, gamma)
    m = screening_mass(beta, gamma_c)
    I_Y = yukawa_overlap_exact(d12, d13, d23, m, n_quad=n_quad)
    # Potential vertex: cubic (+2*beta) + quartic (-6*gamma)
    return (2.0 * beta - 6.0 * gamma_c) * C1 * C2 * C3 * I_Y


def total_three_body_energy_exact(positions, C_values, beta=None, gamma=None,
                                   n_quad=40):
    """
    Total exact 3-body energy for N bodies.

        V3_total = sum_{i<j<k} V_3(i,j,k)

    Parameters
    ----------
    positions : array_like, shape (N, 3)
    C_values  : array_like, shape (N,)
    beta, gamma : float
    n_quad : int

    Returns
    -------
    V3 : float
    """
    beta, gamma_c = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values  = np.asarray(C_values,  dtype=float)
    N = len(positions)
    if N < 3:
        return 0.0

    V3 = 0.0
    for i, j, k in combinations(range(N), 3):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        d_ik = np.linalg.norm(positions[i] - positions[k])
        d_jk = np.linalg.norm(positions[j] - positions[k])
        V3 += three_body_energy_exact(
            d_ij, d_ik, d_jk,
            C_values[i], C_values[j], C_values[k],
            beta=beta, gamma=gamma_c, n_quad=n_quad)
    return V3


# ── 3-body forces (exact) ────────────────────────────────────────────────


def three_body_force_triplet_exact(pos1, pos2, pos3, C1, C2, C3,
                                    beta=None, gamma=None, n_quad=40):
    """
    Exact 3-body forces on all three bodies in a triplet.

    Returns F1, F2, F3 such that F1 + F2 + F3 = 0 (Newton's 3rd law).

    Derivation:
    -----------
    V_3 = (2*beta - 6*gamma)*C1*C2*C3 * I_Y(d12, d13, d23; m)

    F_1 = -dV_3/dx_1
        = (6*gamma - 2*beta)*C1*C2*C3 * [dI_Y/dd12 * (x1-x2)/d12
                                         + dI_Y/dd13 * (x1-x3)/d13]

    (analogous for F_2, F_3)

    Parameters
    ----------
    pos1, pos2, pos3 : array_like, shape (3,)   3D positions
    C1, C2, C3       : float                    source strengths
    beta, gamma      : float                    TGP couplings
    n_quad           : int                      quadrature resolution

    Returns
    -------
    F1, F2, F3 : ndarray, shape (3,)   force vectors
    """
    beta, gamma_c = default_beta_gamma(beta, gamma)
    m = screening_mass(beta, gamma_c)

    p1 = np.asarray(pos1, dtype=float)
    p2 = np.asarray(pos2, dtype=float)
    p3 = np.asarray(pos3, dtype=float)

    r12 = p1 - p2;  d12 = np.linalg.norm(r12)
    r13 = p1 - p3;  d13 = np.linalg.norm(r13)
    r23 = p2 - p3;  d23 = np.linalg.norm(r23)

    if d12 < 1e-15 or d13 < 1e-15 or d23 < 1e-15:
        return np.zeros(3), np.zeros(3), np.zeros(3)

    # Exact derivatives dI_Y/dd_ij
    dI12, dI13, dI23 = _dI_dd_components(d12, d13, d23, m, n_quad=n_quad)

    # F = -dV3/dx; V3 = (2*beta - 6*gamma)*C1C2C3*I_Y
    # → F = (6*gamma - 2*beta)*C1C2C3 * dI/dx
    coeff = (6.0 * gamma_c - 2.0 * beta) * C1 * C2 * C3

    # Unit vectors r_hat_ij = (xi - xj) / d_ij
    e12 = r12 / d12   # points 1 -> away from 2
    e13 = r13 / d13
    e23 = r23 / d23

    # dI_Y / dx_1 = dI12 * e12 + dI13 * e13
    # F_1 = -dV_3/dx_1 = +6*gamma*C1C2C3 * dI/dx_1
    F1 =  coeff * (dI12 * e12 + dI13 * e13)
    F2 =  coeff * (dI12 * (-e12) + dI23 * e23)
    F3 =  coeff * (dI13 * (-e13) + dI23 * (-e23))

    return F1, F2, F3


def three_body_forces_exact(positions, C_values, beta=None, gamma=None,
                             n_quad=40):
    """
    Total exact 3-body forces on all N bodies.

    Sums over all C(N,3) triplets. For each triplet (i,j,k), computes
    the three exact force vectors and adds them to the respective bodies.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
    C_values  : array_like, shape (N,)
    beta, gamma : float
    n_quad : int

    Returns
    -------
    forces : ndarray, shape (N, 3)
        3-body force on each body (summed over all triplets).
    """
    positions = np.asarray(positions, dtype=float)
    C_values  = np.asarray(C_values,  dtype=float)
    N = len(positions)
    forces = np.zeros_like(positions)

    if N < 3:
        return forces

    for i, j, k in combinations(range(N), 3):
        Fi, Fj, Fk = three_body_force_triplet_exact(
            positions[i], positions[j], positions[k],
            C_values[i],  C_values[j],  C_values[k],
            beta=beta, gamma=gamma, n_quad=n_quad)
        forces[i] += Fi
        forces[j] += Fj
        forces[k] += Fk

    return forces


# ── Accuracy diagnostics ─────────────────────────────────────────────────


def benchmark_vs_numerical(d12, d13, d23, m=1.0, n_quad=40):
    """
    Compare exact Feynman integral against brute-force 3D numerical quadrature.

    Returns
    -------
    dict with keys: 'exact', 'numerical', 'relative_error'
    """
    from .three_body_terms import triple_overlap_numerical
    import numpy as np

    I_exact = yukawa_overlap_exact(d12, d13, d23, m, n_quad=n_quad)

    # Build sources at known positions (center at origin, others at d12, d13)
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([d12, 0.0, 0.0])
    # Place p3 so that |p1-p3|=d13, |p2-p3|=d23
    cos_ang = (d12**2 + d13**2 - d23**2) / (2*d12*d13 + 1e-30)
    cos_ang = np.clip(cos_ang, -1, 1)
    p3 = np.array([d13*cos_ang, d13*np.sqrt(1-cos_ang**2), 0.0])

    d_max = max(d12, d13, d23)
    L = max(20.0, 5*d_max)
    sources = [(p1, 1.0, m), (p2, 1.0, m), (p3, 1.0, m)]
    I_numerical = triple_overlap_numerical(sources, L=L, N=100)

    rel_err = abs(I_exact - I_numerical) / (abs(I_numerical) + 1e-30)
    return {'exact': I_exact, 'numerical': I_numerical, 'relative_error': rel_err}


def check_newton_3rd_law(pos1, pos2, pos3, C1=1.0, C2=1.0, C3=1.0,
                          beta=1.0, gamma=1.0, n_quad=50):
    """
    Verify F1 + F2 + F3 = 0 to machine precision.

    Returns residual norm |F1+F2+F3|.
    """
    F1, F2, F3 = three_body_force_triplet_exact(
        pos1, pos2, pos3, C1, C2, C3, beta=beta, gamma=gamma, n_quad=n_quad)
    residual = np.linalg.norm(F1 + F2 + F3)
    return residual, F1, F2, F3
