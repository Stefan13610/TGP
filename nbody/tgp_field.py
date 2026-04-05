"""
tgp_field.py -- Core TGP field infrastructure
================================================
Dimensionless units throughout: x = r/r_0, Phi_hat = Phi/Phi_0.

The TGP scalar field Phi satisfies a nonlinear field equation whose
linearized (weak-field) solution around the vacuum Phi = Phi_0 is:

    Phi(r) = Phi_0 * (1 + delta(r))

where delta(r) = C * exp(-m_sp * r) / r is the Yukawa profile, with
C the dimensionless source strength and m_sp the screening mass.

Vacuum condition
----------------
For beta = gamma (the vacuum self-interaction balance):
    m_sp = sqrt(3*gamma - 2*beta) = sqrt(gamma)

Kinetic coupling — two forms
-----------------------------
The kinetic coupling K(g) (where g = Phi/Phi_0) has two representations:

  (a) FULL substrate form (sek10_N0_wyprowadzenie.tex, prop:f-from-substrate):
      K_sub(g) = K_geo * g^2
      Ghost-free (K_sub > 0 for all g > 0). Valid in all regimes.
      ODE solitonu: g^2*g'' + g*(g')^2 + (2/r)*g^2*g' = V'(g)

  (b) LPA (logarithmic) approximation (one-loop ERG truncation):
      f(g) = 1 + 2*alpha*ln(g),  alpha = 2
      First-order truncation of K_sub = g^2 = exp(2*ln(g)).
      Valid only for |ln g| << 1 (weak fields). Develops ghosts for g << 1.

The manuscript (as of 2026-04) uses form (a) as the definitive version.
This module supports both via KINETIC_MODE.

Energy functional (eq:energy-corrected):
    E[Phi] = integral { (1/2) K(Phi) (grad Phi)^2
                       + (beta / 3 Phi_0) Phi^3
                       - (gamma / 4 Phi_0^2) Phi^4 } d^3x

Source sign convention: SOURCE_SIGN = -1 (mass generates excess space,
source term -q*Phi_0*rho on the RHS of the field equation).
"""

import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────
ALPHA = 2.0          # Gradient coupling for LPA mode
SOURCE_SIGN = -1     # Mass -> excess space
K_GEO = 1.0          # Geometric prefactor for K_sub (normalized to 1)

# Kinetic coupling mode: "full" (K_sub = g^2) or "lpa" (f = 1 + 4*ln(g))
# Default: "full" — matches current manuscript (ghost-free, full resummation)
KINETIC_MODE = "full"


def kinetic_coupling(g, mode=None):
    """Kinetic coupling function K(g) where g = Phi/Phi_0.

    Parameters
    ----------
    g : array_like
        Normalized field value Phi/Phi_0.
    mode : str or None
        "full" for K_sub = K_geo * g^2 (ghost-free, manuscript default),
        "lpa" for f = 1 + 2*ALPHA*ln(g) (one-loop ERG truncation).
        If None, uses module-level KINETIC_MODE.

    Returns
    -------
    K : ndarray
        Kinetic coupling at each point.
    """
    if mode is None:
        mode = KINETIC_MODE
    g = np.asarray(g, dtype=float)
    if mode == "full":
        return K_GEO * g**2
    elif mode == "lpa":
        ratio = np.maximum(g, 1e-30)
        return 1.0 + 2.0 * ALPHA * np.log(ratio)
    else:
        raise ValueError(f"Unknown kinetic mode: {mode!r}. Use 'full' or 'lpa'.")


def default_beta_gamma(beta=None, gamma=None):
    """Resolve beta, gamma with vacuum condition beta=gamma as default.

    If both are None, returns (1.0, 1.0).
    If only one is given, the other is set equal (vacuum condition).
    """
    if beta is None and gamma is None:
        return 1.0, 1.0
    if beta is not None and gamma is None:
        return beta, beta
    if gamma is not None and beta is None:
        return gamma, gamma
    return beta, gamma


def screening_mass(beta, gamma=None):
    """Screening mass m_sp = sqrt(3*gamma - 2*beta).

    For vacuum (beta=gamma): m_sp = sqrt(gamma).
    """
    if gamma is None:
        gamma = beta  # vacuum condition
    arg = 3.0 * gamma - 2.0 * beta
    if arg < 0:
        raise ValueError(f"Tachyonic regime: 3*gamma - 2*beta = {arg} < 0")
    return np.sqrt(arg)


# ── Single-source profiles ────────────────────────────────────────────────

def yukawa_profile(r, C, m_sp):
    """Single-source linearized profile.

    delta(r) = C * exp(-m_sp * r) / r

    This is the fractional field excess: delta = (Phi - Phi_0) / Phi_0.
    For a point source of strength C with screening mass m_sp.

    Parameters
    ----------
    r : array_like
        Radial distance (dimensionless, in units of r_0).
    C : float
        Dimensionless source strength.
    m_sp : float
        Screening mass (inverse screening length).

    Returns
    -------
    delta : ndarray
        Profile value delta(r).
    """
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-12)
    return C * np.exp(-m_sp * r_safe) / r_safe


def phi_single(r, C, beta):
    """Full profile Phi/Phi_0 = 1 + delta for a single source in vacuum.

    Uses vacuum condition beta = gamma so m_sp = sqrt(beta).

    Parameters
    ----------
    r : array_like
        Radial distance.
    C : float
        Source strength.
    beta : float
        Self-interaction coupling (= gamma in vacuum).

    Returns
    -------
    phi_hat : ndarray
        Phi / Phi_0 = 1 + C * exp(-m_sp * r) / r
    """
    m_sp = np.sqrt(beta)
    return 1.0 + yukawa_profile(r, C, m_sp)


def grad_phi_single(r, C, m_sp):
    """Radial gradient of the single-source Yukawa profile.

    d(delta)/dr = -C * exp(-m_sp * r) * (m_sp + 1/r) / r

    This is the derivative of delta(r), not of Phi itself.
    The gradient of Phi/Phi_0 equals the gradient of delta (since Phi_0 is constant).

    Parameters
    ----------
    r : array_like
        Radial distance.
    C : float
        Source strength.
    m_sp : float
        Screening mass.

    Returns
    -------
    grad_delta : ndarray
        d(delta)/dr  (negative, pointing inward).
    """
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-12)
    return -C * np.exp(-m_sp * r_safe) * (m_sp + 1.0 / r_safe) / r_safe


# ── Multi-source superposition ────────────────────────────────────────────

def phi_n_sources_linear(x, sources):
    """Linear superposition of Yukawa profiles (weak-field approximation).

    Phi(x)/Phi_0 = 1 + sum_i delta_i(|x - x_i|)

    This is the linearized (Born) approximation: valid when all delta_i << 1.
    Nonlinear corrections (which produce the 3-body forces) are handled
    separately in three_body_terms.py.

    Parameters
    ----------
    x : array_like, shape (3,) or (M, 3)
        Evaluation point(s) in 3D.
    sources : list of tuples
        Each element is (position_3d, C_i, m_sp_i) where
        position_3d is array-like of shape (3,).

    Returns
    -------
    phi_hat : ndarray
        Phi(x) / Phi_0 at each evaluation point.
    """
    x = np.asarray(x, dtype=float)
    single_point = (x.ndim == 1)
    if single_point:
        x = x[np.newaxis, :]  # shape (1, 3)

    phi_hat = np.ones(x.shape[0])
    for pos_i, C_i, m_sp_i in sources:
        pos_i = np.asarray(pos_i, dtype=float)
        r_i = np.linalg.norm(x - pos_i[np.newaxis, :], axis=1)
        phi_hat += yukawa_profile(r_i, C_i, m_sp_i)

    if single_point:
        return phi_hat[0]
    return phi_hat


# ── Energy density ────────────────────────────────────────────────────────

def energy_density(phi, grad_phi_sq, beta, gamma, Phi0=1.0, mode=None):
    """TGP energy density at a point.

    E = (1/2) K(g) |grad phi|^2 + (beta/3) phi^3/Phi0 - (gamma/4) phi^4/Phi0^2

    where K(g) depends on KINETIC_MODE:
      - "full": K = K_geo * g^2  (manuscript default, ghost-free)
      - "lpa":  K = 1 + 2*alpha*ln(g)  (one-loop truncation)

    Parameters
    ----------
    phi : array_like
        Field value Phi (NOT Phi/Phi_0; the actual field).
    grad_phi_sq : array_like
        |grad Phi|^2 at the same points.
    beta, gamma : float
        Self-interaction couplings.
    Phi0 : float
        Vacuum field value (default 1.0 in dimensionless units).
    mode : str or None
        Kinetic coupling mode. If None, uses module-level KINETIC_MODE.

    Returns
    -------
    e : ndarray
        Energy density at each point.
    """
    phi = np.asarray(phi, dtype=float)
    grad_phi_sq = np.asarray(grad_phi_sq, dtype=float)

    g = phi / Phi0
    K = kinetic_coupling(g, mode=mode)

    kinetic = 0.5 * K * grad_phi_sq
    cubic = (beta / 3.0) * phi**3 / Phi0
    quartic = -(gamma / 4.0) * phi**4 / Phi0**2

    return kinetic + cubic + quartic
