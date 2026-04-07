"""
three_body_terms.py -- Irreducible 3-body interactions from TGP nonlinearity
=============================================================================
THIS IS THE NOVEL PHYSICS.

STATUS OF RESULTS:
  APPROXIMATE:  triple_overlap_analytic() — saddle-point estimate of the
                triple Yukawa overlap integral.
                WARNING (revised): actual error is 160-770% for equilateral
                triangles with t=m*d in [0.5, 6].  The formula uses the
                wrong exponential exponent (exp(-3t/2) vs exact exp(-sqrt(3)*t)).
                USE three_body_force_exact.yukawa_overlap_exact() INSTEAD.
  NUMERICAL:    triple_overlap_numerical() — direct 3D grid quadrature.
                Converges as grid resolution increases. PREFERRED for
                quantitative work.
  EXACT:        The Coulomb limit (m_sp → 0) gives I_triple = 8*pi^2/P,
                where P = d12 + d13 + d23 (perimeter). Used in dynamics_v2.py
                for the unscreened 3-body force.

Origin of the 3-body force
---------------------------
The TGP energy functional contains nonlinear self-interaction terms:

    E[Phi] = ... + (beta/3 Phi_0) Phi^3 - (gamma/4 Phi_0^2) Phi^4

When we write the total field as a superposition of single-source profiles
plus the vacuum:

    Phi = Phi_0 * (1 + delta_1 + delta_2 + delta_3)

with delta_i = C_i * exp(-m_sp * r_i) / r_i being the Yukawa profile of
source i, the QUARTIC term -(gamma/4*Phi_0^2) * Phi^4 produces cross-terms.

Expanding (1 + d1 + d2 + d3)^4 via the multinomial theorem, the
irreducible 3-body cross-term (one power of 1, one each of d1, d2, d3)
has multinomial coefficient:

    4! / (1! * 1! * 1! * 1!) = 24

So the 3-body contribution from the quartic term is:

    V_3 = -(gamma/4) * Phi_0^2 * 24 * integral { delta_1 * delta_2 * delta_3 } d^3x
        = -6 * gamma * C1*C2*C3 * I_triple_bare

where I_triple_bare is the triple overlap integral of the SCREENED Yukawa
profiles (not bare 1/r, which would diverge).

Triple overlap integral
------------------------
For Yukawa profiles delta_i = exp(-m*r_i)/r_i:

    I_triple_bare = integral { exp(-m*|x-x1|)/|x-x1| * exp(-m*|x-x2|)/|x-x2|
                              * exp(-m*|x-x3|)/|x-x3| } d^3x

This integral has NO known simple closed form.  We provide:
1. Numerical evaluation on a 3D grid (triple_overlap_numerical) -- PREFERRED
2. An approximate analytic estimate for order-of-magnitude scaling

For comparison, the PAIRWISE Yukawa overlap has the exact result (derived via
Gegenbauer expansion of the angular integral):

    I_pair(d, m) = 2*pi * exp(-m*d) / m

The triple overlap decays faster than exp(-m * semi_perimeter) due to the
geometric constraint that the integration point must be within ~ 1/m of
ALL three sources simultaneously.

Suppression by C
-----------------
The 3-body energy scales as C1*C2*C3 ~ C^3, compared to the pairwise
energy ~ C^2.  For elementary particles (C << 1), 3-body effects are
suppressed by an extra factor of C.  For macroscopic bodies (C >> 1),
3-body effects become significant and could break the Newtonian degeneracy.
"""

import numpy as np
from itertools import combinations
from .tgp_field import default_beta_gamma, yukawa_profile, screening_mass


# ── Triple overlap integral ───────────────────────────────────────────────

def triple_overlap_analytic(d12, d13, d23, C1=1.0, C2=1.0, C3=1.0,
                            m_sp=None, beta=None):
    """Approximate analytic triple overlap integral for Yukawa profiles.

    For three Yukawa profiles exp(-m*r)/r centered at vertices of a
    triangle with sides d12, d13, d23, this returns an ORDER-OF-MAGNITUDE
    estimate based on dimensional analysis and exponential scaling.

    The approximation:

        I_triple ~ (2*pi)^(3/2) / m^2 * exp(-m * s) / s

    where s = (d12 + d13 + d23) / 2 is the semi-perimeter.

    WARNING (revised 2026-03): Comparison against the exact Feynman 2D
    integral shows errors of 160--770% for equilateral triangles with
    t = m*d in [0.5, 6].  The exponential decay rate is WRONG:
    this formula uses exp(-3t/2) while the exact result decays as
    exp(-sqrt(3)*t) for large t (with sqrt(3) > 3/2, so the exact
    result decays FASTER).  The formula overestimates I_Y by a growing
    factor as t increases.

    For quantitative 3-body results, use:
        three_body_force_exact.yukawa_overlap_exact()  (2D integral, EXACT)
    or:
        triple_overlap_numerical()  (3D grid, converges slowly due to 1/r singularity)

    The exponential scaling exp(-m*s) is robust; the polynomial prefactor
    is the uncertain part.

    Parameters
    ----------
    d12, d13, d23 : float
        Pairwise distances between the three sources.
    C1, C2, C3 : float
        Source strengths (included as overall prefactor).
    m_sp : float or None
        Screening mass.  If None, computed from beta via vacuum condition.
    beta : float or None
        Self-interaction coupling (used only if m_sp is None).

    Returns
    -------
    I : float
        Approximate triple overlap integral (with C factors).
    """
    if m_sp is None:
        if beta is None:
            beta = 1.0
        m_sp = screening_mass(beta)

    s = (d12 + d13 + d23) / 2.0  # semi-perimeter
    if s < 1e-15:
        return np.inf

    # Order-of-magnitude estimate.
    # The (2*pi)^(3/2) prefactor and 1/(m^2 * s) polynomial are approximate.
    # The exponential decay exp(-m*s) captures the dominant distance dependence.
    I_bare = (2.0 * np.pi) ** 1.5 / m_sp**2 * np.exp(-m_sp * s) / s

    return C1 * C2 * C3 * I_bare


def triple_overlap_numerical(sources, L=50.0, N=80):
    """Numerical evaluation of the triple overlap integral on a 3D grid.

    Computes:
        I = integral { prod_i [C_i * exp(-m_i*|x-x_i|) / |x-x_i|] } d^3x

    using a midpoint-rule quadrature on a uniform cubic grid.

    This is the PREFERRED method for quantitative 3-body energy calculations.
    The analytic formula (triple_overlap_analytic) provides only an
    order-of-magnitude estimate.

    Grid convergence:  The integral converges as O(1/N) in each dimension.
    For m*d > 1, N=80 with L = max(20, 5*d_max) typically gives < 5% error.
    Increase N for higher accuracy.

    Parameters
    ----------
    sources : list of tuples
        [(pos1, C1, m_sp1), (pos2, C2, m_sp2), (pos3, C3, m_sp3)]
        where pos_i is a 3D position array.  The m_sp values are used
        for Yukawa screening; set m_sp ~ 0 for pure 1/r (will diverge).
    L : float
        Half-size of the integration cube: domain is [-L, L]^3.
    N : int
        Number of grid points per dimension.  Total cost is O(N^3).

    Returns
    -------
    I : float
        Numerical estimate of the overlap integral.
    """
    x_1d = np.linspace(-L, L, N)
    dx = x_1d[1] - x_1d[0]
    dV = dx**3
    X, Y, Z = np.meshgrid(x_1d, x_1d, x_1d, indexing='ij')
    points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    product = np.ones(points.shape[0])
    for pos_i, C_i, m_sp_i in sources:
        pos_i = np.asarray(pos_i, dtype=float)
        r_i = np.linalg.norm(points - pos_i[np.newaxis, :], axis=1)
        if m_sp_i < 1e-10:
            r_safe = np.maximum(r_i, dx * 0.3)
            profile_i = C_i / r_safe
        else:
            profile_i = yukawa_profile(r_i, C_i, m_sp_i)
        product *= profile_i

    return np.sum(product) * dV


# ── 3-body energy ─────────────────────────────────────────────────────────

def three_body_energy(d12, d13, d23, C1, C2, C3, beta=None, gamma=None,
                      use_numerical=False, positions=None):
    """Irreducible 3-body energy from the quartic self-interaction.

    From the expansion of -(gamma/4*Phi_0^2) * Phi^4 with
    Phi = Phi_0*(1 + d1 + d2 + d3), the irreducible 3-body cross-term is:

        V_3 = -(gamma/4) * Phi_0^2 * 24 * C1*C2*C3 * I_triple_bare
            = -6 * gamma * C1*C2*C3 * I_triple_bare

    where I_triple_bare is the triple overlap integral of the Yukawa profiles
    (without C factors), and the coefficient 24 = 4!/(1!1!1!1!) is the
    multinomial coefficient for the term 1^1 * d1^1 * d2^1 * d3^1 in
    (1 + d1 + d2 + d3)^4.

    Physical interpretation: the factor -6*gamma represents the strength
    of the quartic nonlinearity.  The sign is ATTRACTIVE (V_3 < 0 for
    gamma > 0), meaning three bodies attract each other MORE than the
    sum of pairwise attractions.  This is the TGP analogue of the
    Axilrod-Teller-Muto three-body dispersion force.

    Parameters
    ----------
    d12, d13, d23 : float
        Pairwise distances forming the triangle.
    C1, C2, C3 : float
        Source strengths.
    beta, gamma : float or None
        Couplings (default: vacuum beta=gamma=1).
    use_numerical : bool
        If True and positions is provided, use numerical grid integration
        instead of the approximate analytic formula.
    positions : array_like or None
        Required if use_numerical=True; shape (3, 3) giving 3D positions.

    Returns
    -------
    V3 : float
        Irreducible 3-body interaction energy.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    m_sp = screening_mass(beta, gamma)

    if use_numerical and positions is not None:
        sources = [
            (np.asarray(positions[0]), 1.0, m_sp),
            (np.asarray(positions[1]), 1.0, m_sp),
            (np.asarray(positions[2]), 1.0, m_sp),
        ]
        # Adapt grid to source separation
        d_max = max(d12, d13, d23)
        L = max(20.0, 5.0 * d_max)
        I_bare = triple_overlap_numerical(sources, L=L, N=80)
        # Potential vertex: cubic (beta/3)g^3 contributes +2*beta,
        # quartic -(gamma/4)g^4 contributes -6*gamma.
        # Full potential vertex coefficient: (2*beta - 6*gamma).
        V3 = (2.0 * beta - 6.0 * gamma) * C1 * C2 * C3 * I_bare
    else:
        # Use approximate analytic formula (order-of-magnitude).
        # triple_overlap_analytic already includes C factors.
        I_triple = triple_overlap_analytic(d12, d13, d23, C1, C2, C3,
                                           m_sp=m_sp)
        V3 = (2.0 * beta - 6.0 * gamma) * I_triple

    return V3


# ── 3-body force ──────────────────────────────────────────────────────────

def three_body_force_on_1(positions, C_values, beta=None, gamma=None):
    """Force on body 1 from the irreducible 3-body interaction.

    Computed by central finite differences of three_body_energy with
    respect to the position of body 1.

    For N >= 3 bodies, sums over all triplets containing body 1.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Positions of all N bodies.
    C_values : array_like, shape (N,)
        Source strengths.
    beta, gamma : float or None
        Couplings (default: vacuum).

    Returns
    -------
    F1 : ndarray, shape (3,)
        Force vector on body 1 from 3-body terms.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)

    if N < 3:
        return np.zeros(3)

    eps = 1e-6
    F1 = np.zeros(3)

    for dim in range(3):
        pos_plus = positions.copy()
        pos_minus = positions.copy()
        pos_plus[0, dim] += eps
        pos_minus[0, dim] -= eps

        V_plus = 0.0
        V_minus = 0.0

        for j in range(1, N):
            for k in range(j + 1, N):
                d0j_p = np.linalg.norm(pos_plus[0] - positions[j])
                d0k_p = np.linalg.norm(pos_plus[0] - positions[k])
                djk = np.linalg.norm(positions[j] - positions[k])
                V_plus += three_body_energy(d0j_p, d0k_p, djk,
                                            C_values[0], C_values[j],
                                            C_values[k], beta, gamma)

                d0j_m = np.linalg.norm(pos_minus[0] - positions[j])
                d0k_m = np.linalg.norm(pos_minus[0] - positions[k])
                V_minus += three_body_energy(d0j_m, d0k_m, djk,
                                             C_values[0], C_values[j],
                                             C_values[k], beta, gamma)

        F1[dim] = -(V_plus - V_minus) / (2.0 * eps)

    return F1


def total_three_body_energy(positions, C_values, beta=None, gamma=None):
    """Sum of all irreducible 3-body energies for N bodies.

    V_3_total = sum_{i<j<k} V_3(i, j, k)

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Body positions.
    C_values : array_like, shape (N,)
        Source strengths.
    beta, gamma : float or None
        Couplings.

    Returns
    -------
    V3_total : float
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    positions = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    N = len(positions)

    if N < 3:
        return 0.0

    V3 = 0.0
    for i, j, k in combinations(range(N), 3):
        d_ij = np.linalg.norm(positions[i] - positions[j])
        d_ik = np.linalg.norm(positions[i] - positions[k])
        d_jk = np.linalg.norm(positions[j] - positions[k])
        V3 += three_body_energy(d_ij, d_ik, d_jk,
                                C_values[i], C_values[j], C_values[k],
                                beta, gamma)
    return V3
