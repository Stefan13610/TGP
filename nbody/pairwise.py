"""
pairwise.py -- 2-body effective potential and forces in TGP
===========================================================
Analytical interaction energy V_eff(d) between two point sources,
derived from overlap integrals of single-source Yukawa profiles.

The three contributions (Theorem thm:three-regimes, sek08):

    E_grad(d)  = -4*pi*C1*C2 / d          (gradient overlap, attractive)
    E_beta(d)  = +8*pi*beta*C1*C2 / d^2   (quadratic self-interaction, repulsive)
    E_gamma(d) = -24*pi*gamma*C1*C2*(C1+C2)/(2*d^3)  (cubic, confining)

For equal masses (C1=C2=C):
    E_gamma(d) = -24*pi*gamma*C^3 / d^3

Three regimes exist when beta > 9*C/2 (trivially satisfied for C << 1).
"""

import numpy as np
from .tgp_field import default_beta_gamma


def V_eff(d, C1, C2, beta=None, gamma=None):
    """Two-body effective potential in TGP.

    Parameters
    ----------
    d : float or array_like
        Separation distance (dimensionless).
    C1, C2 : float
        Source strengths of the two bodies.
    beta, gamma : float or None
        Self-interaction couplings.  Default: vacuum condition beta=gamma=1.

    Returns
    -------
    V : float or ndarray
        Total interaction energy.
    V_grad, V_beta, V_gamma : float or ndarray
        Individual contributions.
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    d = np.asarray(d, dtype=float)
    d_safe = np.maximum(d, 1e-10)

    # Gradient overlap (attractive, 1/d)
    V_grad = -4.0 * np.pi * C1 * C2 / d_safe

    # Beta term: repulsive, 1/d^2
    V_beta = 8.0 * np.pi * beta * C1 * C2 / d_safe**2

    # Gamma term: confining, 1/d^3
    # For unequal masses: coefficient involves (C1+C2)/2
    V_gamma = -24.0 * np.pi * gamma * C1 * C2 * (C1 + C2) / (2.0 * d_safe**3)

    V_total = V_grad + V_beta + V_gamma
    return V_total, V_grad, V_beta, V_gamma


def V_eff_total(d, C1, C2, beta=None, gamma=None):
    """Return only the total 2-body potential (scalar convenience function)."""
    return V_eff(d, C1, C2, beta, gamma)[0]


def force_2body(d, C1, C2, beta=None, gamma=None):
    """Force between two bodies: F = -dV/dd.

    Positive F means repulsion (force pushes bodies apart).

    Parameters
    ----------
    d : float or array_like
        Separation distance.
    C1, C2 : float
        Source strengths.
    beta, gamma : float or None
        Self-interaction couplings (default: vacuum).

    Returns
    -------
    F : float or ndarray
        Force magnitude (positive = repulsive).
    """
    beta, gamma = default_beta_gamma(beta, gamma)
    d = np.asarray(d, dtype=float)
    d_safe = np.maximum(d, 1e-10)

    # F = -dV/dd
    # V_gamma = -24*pi*gamma*C1*C2*(C1+C2) / (2*d^3)
    # dV_gamma/dd = +36*pi*gamma*C1*C2*(C1+C2) / d^4   (the /2 is absorbed: 24*3/2 = 36)
    # => F_gamma = -dV_gamma/dd = -36*pi*gamma*C1*C2*(C1+C2) / d^4

    F_grad = -4.0 * np.pi * C1 * C2 / d_safe**2          # attractive (negative)
    F_beta = 16.0 * np.pi * beta * C1 * C2 / d_safe**3   # repulsive (positive)
    F_gamma = -36.0 * np.pi * gamma * C1 * C2 * (C1 + C2) / d_safe**4  # attractive (confining)

    return F_grad + F_beta + F_gamma


def force_zeros_2body(C, beta, gamma=None):
    """Find zero-crossings of the 2-body force for equal masses.

    For equal masses (C1=C2=C), F(d) = 0 gives:
        -4*pi*C^2/d^2 + 16*pi*beta*C^2/d^3 - 72*pi*gamma*C^3/d^4 = 0

    Dividing by -4*pi*C^2/d^4:
        d^2 - 4*beta*d + 18*gamma*C = 0

    This is the SAME polynomial as the equilateral equilibrium condition
    (Proposition prop:trzy-rezimy-beta-gamma, sek03).

    Solutions:
        d = 2*beta +/- sqrt(4*beta^2 - 18*gamma*C)

    Three regimes exist when discriminant > 0, i.e., beta > 9*C/2 (for gamma=beta).

    Parameters
    ----------
    C : float
        Source strength (equal masses).
    beta : float
        Self-interaction coupling.
    gamma : float or None
        If None, vacuum condition gamma = beta.

    Returns
    -------
    result : tuple or None
        (d_inner, d_outer) where d_inner is the repulsive barrier
        and d_outer is the well minimum, or None if no real roots.
    """
    if gamma is None:
        gamma = beta

    discriminant = 4.0 * beta**2 - 18.0 * gamma * C
    if discriminant < 0:
        return None

    sqrt_disc = np.sqrt(discriminant)
    d_inner = 2.0 * beta - sqrt_disc  # closer to origin
    d_outer = 2.0 * beta + sqrt_disc  # further from origin

    # Both roots must be positive
    if d_inner <= 0 and d_outer <= 0:
        return None
    if d_inner <= 0:
        return (d_outer,)  # only one physical root

    return (d_inner, d_outer)
