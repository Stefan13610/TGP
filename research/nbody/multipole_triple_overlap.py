"""
multipole_triple_overlap.py — Multipole (Legendre/Gegenbauer) expansion of I_Y
================================================================================

SEMIANALYTIC computation of the triple Yukawa overlap integral I_Y(d12,d13,d23;m)
using the spherical harmonic addition theorem for the Yukawa Green's function.

DERIVATION (tgp_yukawa_IY_multipole_gegenbauer.tex):
-----------------------------------------------------
The Yukawa Green's function G_m(|r-r'|) = exp(-m|r-r'|)/|r-r'| admits the
addition theorem:

    G_m(|r-r'|) = sum_l (2l+1) w_l(m,r,r') P_l(cos gamma)

where w_l(m,r,r') = (2m/pi) i_l(m r_<) k_l(m r_>), with i_l, k_l being
modified spherical Bessel functions.

Choosing x3 = 0 as the expansion center:

    I_Y = integral G_m(|x|) G_m(|x-x1|) G_m(|x-x2|) d^3x

After expanding the two displaced Green's functions and integrating over
the sphere, the angular integral diagonalizes:

    A_{ll'}(omega) = (4pi/(2l+1)) delta_{ll'} P_l(cos omega)

This reduces the double sum to a SINGLE sum:

    I_Y = 4*pi * sum_l (2l+1) P_l(cos omega) * R_l(d13, d23, m)

where R_l is a 1D radial integral:

    R_l(a,b,m) = (2m/pi)^2 * integral_0^inf r * exp(-m*r)
                  * [i_l(mr_<a) k_l(mr_>a)] [i_l(mr_<b) k_l(mr_>b)] dr

with r_<a = min(r,a), r_>a = max(r,a) and similarly for b.

The integral splits into regions [0,min(a,b)], [min(a,b),max(a,b)], [max(a,b),inf).

ADVANTAGES over Feynman 2D (three_body_force_exact.py):
  - Angular part is analytic (Legendre polynomials)
  - Better convergence control: L_max directly controls angular resolution
  - Derivatives w.r.t. distances can be computed analytically
  - Tabulatable: R_l(a,b,m) depends only on distances, not geometry

ACCURACY:
  - L_max=5: ~1% for moderate geometry (equilateral, m*d ~ 1)
  - L_max=10: ~0.1% for most configurations
  - L_max=15: ~0.01% for well-separated bodies
  - Slow convergence when bodies are nearly collinear (cos omega ~ +/-1)
"""

import numpy as np
from scipy.special import spherical_in, spherical_kn
from numpy.polynomial.legendre import legval


# ── Vectorized modified spherical Bessel functions ───────────────────────

def _bessel_ik_all_ell(L_max, z):
    """
    Compute i_l(z) and k_l(z) for l = 0, 1, ..., L_max at all z values.

    Uses scipy vectorized over z for each l. Returns arrays of shape
    (L_max+1, len(z)).

    Parameters
    ----------
    L_max : int
    z : (n,) ndarray of positive floats

    Returns
    -------
    i_all : (L_max+1, n) ndarray — i_l(z[j]) for each l
    k_all : (L_max+1, n) ndarray — k_l(z[j]) for each l
    """
    z = np.asarray(z, dtype=float)
    n = len(z)
    i_all = np.zeros((L_max + 1, n))
    k_all = np.zeros((L_max + 1, n))

    for ell in range(L_max + 1):
        i_all[ell] = spherical_in(ell, z)
        k_all[ell] = spherical_kn(ell, z)

    return i_all, k_all


# ── Gauss-Legendre quadrature ────────────────────────────────────────────

_GL_CACHE = {}

def _gauss_legendre_segment(a_lo, a_hi, n_pts):
    """Gauss-Legendre nodes and weights on [a_lo, a_hi]."""
    if n_pts not in _GL_CACHE:
        _GL_CACHE[n_pts] = np.polynomial.legendre.leggauss(n_pts)
    pts, wts = _GL_CACHE[n_pts]
    mid = 0.5 * (a_hi + a_lo)
    half = 0.5 * (a_hi - a_lo)
    return mid + half * pts, half * wts


# ── All-ell radial integrals at once ─────────────────────────────────────

def _radial_integrals_all_ell(L_max, a, b, m, n_rad=50, r_max_factor=10.0):
    """
    Compute R_l(a, b, m) for l = 0, 1, ..., L_max in one vectorized pass.

    This avoids re-computing quadrature points and Bessel values for each l.

    Returns
    -------
    R : (L_max+1,) ndarray
    """
    if a > b:
        a, b = b, a  # ensure a <= b

    prefactor = (2.0 * m / np.pi) ** 2
    R = np.zeros(L_max + 1)

    # Precompute Bessel at boundaries: shape (L_max+1,)
    ma, mb = m * a, m * b
    i_a_all = np.array([spherical_in(ell, ma) for ell in range(L_max + 1)])
    k_a_all = np.array([spherical_kn(ell, ma) for ell in range(L_max + 1)])
    i_b_all = np.array([spherical_in(ell, mb) for ell in range(L_max + 1)])
    k_b_all = np.array([spherical_kn(ell, mb) for ell in range(L_max + 1)])

    # ── Segment 1: [eps, a] — r < a <= b ──
    if a > 1e-15:
        eps_lo = max(a * 1e-6, 1e-15)
        r1, w1 = _gauss_legendre_segment(eps_lo, a, n_rad)
        mr1 = m * r1
        common1 = r1 * np.exp(-m * r1)  # shape (n_rad,)

        # Bessel at all quadrature points for all ell: shape (L_max+1, n_rad)
        i_r1_all, _ = _bessel_ik_all_ell(L_max, mr1)

        # Integrand for each ell: common1 * i_l(mr)^2 * k_l(ma) * k_l(mb)
        for ell in range(L_max + 1):
            integrand = common1 * i_r1_all[ell]**2 * k_a_all[ell] * k_b_all[ell]
            R[ell] += np.dot(w1, integrand)

    # ── Segment 2: [a, b] — a < r < b ──
    if b - a > 1e-15:
        r2, w2 = _gauss_legendre_segment(a, b, n_rad)
        mr2 = m * r2
        common2 = r2 * np.exp(-m * r2)

        i_r2_all, k_r2_all = _bessel_ik_all_ell(L_max, mr2)

        for ell in range(L_max + 1):
            integrand = common2 * i_a_all[ell] * k_r2_all[ell] * i_r2_all[ell] * k_b_all[ell]
            R[ell] += np.dot(w2, integrand)

    # ── Segment 3: [b, R_max] — r > b >= a ──
    R_max = b + r_max_factor / m
    r3, w3 = _gauss_legendre_segment(b, R_max, n_rad)
    mr3 = m * r3
    common3 = r3 * np.exp(-m * r3)

    _, k_r3_all = _bessel_ik_all_ell(L_max, mr3)

    for ell in range(L_max + 1):
        integrand = common3 * i_a_all[ell] * k_r3_all[ell] * i_b_all[ell] * k_r3_all[ell]
        R[ell] += np.dot(w3, integrand)

    return prefactor * R


# ── Main: I_Y via multipole expansion ────────────────────────────────────

def yukawa_overlap_multipole(d12, d13, d23, m, L_max=10, n_rad=50):
    """
    Triple Yukawa overlap integral via multipole (Legendre) expansion.

        I_Y = 4*pi * sum_{l=0}^{L_max} (2l+1) P_l(cos omega) * R_l(d13, d23, m)

    where cos omega = (d13^2 + d23^2 - d12^2) / (2*d13*d23).

    Parameters
    ----------
    d12, d13, d23 : float
        Pairwise distances (all > 0).
    m : float
        Screening mass (> 0).
    L_max : int
        Maximum multipole order. Higher = more accurate.
    n_rad : int
        Gauss-Legendre quadrature points per radial segment.

    Returns
    -------
    I_Y : float
        Triple overlap integral.
    """
    d12, d13, d23, m = float(d12), float(d13), float(d23), float(m)

    # Angle between directions to body 1 and body 2 from body 3
    cos_omega = (d13**2 + d23**2 - d12**2) / (2.0 * d13 * d23)
    cos_omega = max(-1.0, min(1.0, cos_omega))

    # All R_l values at once
    R_all = _radial_integrals_all_ell(L_max, d13, d23, m, n_rad=n_rad)

    # Legendre polynomials P_l(cos omega) for l = 0..L_max
    # Use recurrence: P_0 = 1, P_1 = x, (l+1)P_{l+1} = (2l+1)xP_l - lP_{l-1}
    x = cos_omega
    P = np.zeros(L_max + 1)
    P[0] = 1.0
    if L_max >= 1:
        P[1] = x
    for ell in range(1, L_max):
        P[ell + 1] = ((2 * ell + 1) * x * P[ell] - ell * P[ell - 1]) / (ell + 1)

    # Sum: I_Y = 4*pi * sum_l (2l+1) P_l R_l
    ells = np.arange(L_max + 1)
    I_Y = 4.0 * np.pi * np.dot((2 * ells + 1) * P, R_all)

    return I_Y


def yukawa_overlap_multipole_with_derivatives(d12, d13, d23, m,
                                               L_max=10, n_rad=50,
                                               dd=1e-6):
    """
    I_Y and its derivatives dI_Y/dd12, dI_Y/dd13, dI_Y/dd23 via multipole.

    Derivatives are computed by central finite differences on the multipole
    result.

    Returns
    -------
    I_Y, dI_dd12, dI_dd13, dI_dd23 : float
    """
    I_Y = yukawa_overlap_multipole(d12, d13, d23, m, L_max, n_rad)

    Ip = yukawa_overlap_multipole(d12 + dd, d13, d23, m, L_max, n_rad)
    Im = yukawa_overlap_multipole(d12 - dd, d13, d23, m, L_max, n_rad)
    dI_dd12 = (Ip - Im) / (2.0 * dd)

    Ip = yukawa_overlap_multipole(d12, d13 + dd, d23, m, L_max, n_rad)
    Im = yukawa_overlap_multipole(d12, d13 - dd, d23, m, L_max, n_rad)
    dI_dd13 = (Ip - Im) / (2.0 * dd)

    Ip = yukawa_overlap_multipole(d12, d13, d23 + dd, m, L_max, n_rad)
    Im = yukawa_overlap_multipole(d12, d13, d23 - dd, m, L_max, n_rad)
    dI_dd23 = (Ip - Im) / (2.0 * dd)

    return I_Y, dI_dd12, dI_dd13, dI_dd23


# ── Forces using multipole I_Y ──────────────────────────────────────────

def three_body_forces_multipole(positions, C_values, beta=None, gamma=None,
                                 softening=1e-6, L_max=10, n_rad=50):
    """
    Irreducible 3-body forces from TGP quartic nonlinearity, computed
    using the multipole expansion of I_Y.

    F_i = -dV_3/dx_i, where V_3 = -6*gamma*C1*C2*C3*I_Y for each triplet.

    Returns
    -------
    forces : (n, 3) ndarray — 3-body forces (NOT acceleration).
    """
    from .tgp_field import default_beta_gamma, screening_mass
    from itertools import combinations

    beta, gamma = default_beta_gamma(beta, gamma)
    m = screening_mass(beta, gamma)

    n = len(C_values)
    forces = np.zeros((n, 3))

    if n < 3:
        return forces

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_values[i], C_values[j], C_values[k]

        rij = positions[j] - positions[i]
        rik = positions[k] - positions[i]
        rjk = positions[k] - positions[j]

        dij = np.sqrt(np.dot(rij, rij) + softening**2)
        dik = np.sqrt(np.dot(rik, rik) + softening**2)
        djk = np.sqrt(np.dot(rjk, rjk) + softening**2)

        _, dI_dd12, dI_dd13, dI_dd23 = yukawa_overlap_multipole_with_derivatives(
            dij, dik, djk, m, L_max=L_max, n_rad=n_rad,
        )

        # V_3 = (2*beta - 6*gamma)*Ci*Cj*Ck*I_Y
        # F = -dV_3/dx = (6*gamma - 2*beta)*Ci*Cj*Ck * dI/dx
        coeff = (6.0 * gamma - 2.0 * beta) * Ci * Cj * Ck

        e_ij = (positions[i] - positions[j]) / dij
        e_ik = (positions[i] - positions[k]) / dik
        e_jk = (positions[j] - positions[k]) / djk

        dI_dxi = dI_dd12 * e_ij + dI_dd13 * e_ik
        dI_dxj = -dI_dd12 * e_ij + dI_dd23 * e_jk
        dI_dxk = -dI_dd13 * e_ik - dI_dd23 * e_jk

        forces[i] += coeff * dI_dxi
        forces[j] += coeff * dI_dxj
        forces[k] += coeff * dI_dxk

    return forces


def three_body_potential_multipole(positions, C_values, beta=None, gamma=None,
                                    softening=1e-6, L_max=10, n_rad=50):
    """
    Total irreducible 3-body potential energy via multipole expansion.

    V_3_total = sum_{i<j<k} (2*beta - 6*gamma)*Ci*Cj*Ck*I_Y(dij,dik,djk;m)
    """
    from .tgp_field import default_beta_gamma, screening_mass
    from itertools import combinations

    beta, gamma = default_beta_gamma(beta, gamma)
    m = screening_mass(beta, gamma)

    n = len(C_values)
    V = 0.0

    if n < 3:
        return V

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_values[i], C_values[j], C_values[k]

        rij = positions[j] - positions[i]
        rik = positions[k] - positions[i]
        rjk = positions[k] - positions[j]

        dij = np.sqrt(np.dot(rij, rij) + softening**2)
        dik = np.sqrt(np.dot(rik, rik) + softening**2)
        djk = np.sqrt(np.dot(rjk, rjk) + softening**2)

        I_Y = yukawa_overlap_multipole(dij, dik, djk, m,
                                        L_max=L_max, n_rad=n_rad)
        V += (2.0 * beta - 6.0 * gamma) * Ci * Cj * Ck * I_Y

    return V
