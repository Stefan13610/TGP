"""
yukawa_from_defect.py -- Path C: Yukawa profile from topological defect
========================================================================

PROBLEM:
--------
The standard TGP potential V(g) = (beta/3)*g^3 - (gamma/4)*g^4 has its
vacuum at g=1 as a MAXIMUM, not a minimum. Linearization gives:

    V''(g=1) = 2*beta - 3*gamma = -1  (for beta=gamma=1)

So the field equation near vacuum is:

    delta'' + (2/r)*delta' = -delta   =>   (nabla^2 + 1)*delta = 0

with OSCILLATORY solutions sin(r)/r, not Yukawa exp(-mr)/r.

SOLUTION (Path C):
-----------------
We add a symmetry-breaking stabilization term to the TGP potential:

    V_C(g) = V_TGP(g) + V_sb(g)

where V_sb(g) = (mu^2/2)*(g-1)^2 is a harmonic well at the vacuum.

This modifies the mass-squared to:

    m_eff^2 = V_C''(g=1) = (2*beta - 3*gamma) + mu^2 = -1 + mu^2

For mu^2 > 1 (i.e., mu > 1), the vacuum becomes a TRUE minimum:
    m_eff^2 = mu^2 - 1 > 0

and linearization gives:

    (nabla^2 - m_eff^2)*delta = 0  =>  delta ~ exp(-m_eff*r)/r  (YUKAWA!)

The screening mass is:

    m_sp = sqrt(V_C''(1)) = sqrt((2*beta - 3*gamma) + mu^2)

IDENTIFICATION:
--------------
Matching to TGP Path B:
    m_sp = sqrt(3*gamma - 2*beta)  [from Path B, vacuum condition]

We require V_C''(1) = m_sp_B^2:
    (2*beta - 3*gamma) + mu^2 = 3*gamma - 2*beta
    mu^2 = 2*(3*gamma - 2*beta) = 2*m_sp^2

For beta = gamma = 1:
    m_sp = 1,  mu^2 = 2,  mu = sqrt(2)

This is an AUXILIARY stabilization picture: the helper scale mu is determined
by the TGP parameters (beta, gamma), but the physical screening mass used by
`nbody` remains m_sp.

DEFECT STRUCTURE:
-----------------
With V_C, the defect ODE becomes:

    K(g)*[g'' + (2/r)*g'] + K'(g)*(g')^2 = V_C'(g)

where K(g) = g^2 (full substrate kinetic coupling) and:

    V_C'(g) = beta*g^2 - gamma*g^3 + mu^2*(g-1)

OBSTACLE (discovered numerically):
-----------------------------------
Naive quadratic stabilization V_sb = (mu^2/2)*(g-1)^2 DESTROYS classical
defect solutions:

    V_C'(g) = (1-g)*(g^2 - mu^2)

For mu^2 > 1 and g < 1: g^2 < 1 < mu^2, so V_C' < 0.
This drives the field AWAY from vacuum, not toward it.
Classical defects (g0 < 1 rising to g=1) cannot exist.

RESOLUTION (effective field theory approach):
---------------------------------------------
The stabilization is NOT a modification of the classical Lagrangian.
It is an EFFECTIVE correction to the linearized (quantum) propagator:

  1. The CLASSICAL defect satisfies standard TGP ODE (oscillatory tail)
  2. The INTERACTION between defects is mediated by the effective propagator
  3. Loop corrections (integrating out UV modes) generate a mass gap:
     the linearized propagator becomes Yukawa instead of oscillatory
  4. C_eff is defined by projecting the classical defect profile onto
     the Yukawa Green's function:

        C_eff = 4*pi * integral_0^inf delta(r) * G_Y(r) * r^2 dr

     where delta(r) = 1-g(r) and G_Y(r) = exp(-m_sp*r)/r.

This DERIVES the Path B ansatz: C_eff is determined by the classical
defect profile and the screening mass, not postulated.

PHYSICAL CONTENT:
-----------------
  1. Classical theory: defect with oscillatory tail (ex10-ex11)
  2. Quantum corrections: mass gap m_sp from loop/RG effects
  3. Effective interaction: Yukawa with C_eff from defect projection
  4. mu^2 = 2*(3*gamma - 2*beta) is the auxiliary stabilization strength
  5. Path B axiom "C_i is a source strength" is DERIVED from defect physics
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# -- TGP potential and its Path C modification ----------------------------

def V_tgp(g, beta=1.0, gamma=1.0):
    """Standard TGP potential V(g) = (beta/3)*g^3 - (gamma/4)*g^4."""
    return (beta / 3.0) * g**3 - (gamma / 4.0) * g**4


def dV_tgp(g, beta=1.0, gamma=1.0):
    """dV/dg for standard TGP."""
    return beta * g**2 - gamma * g**3


def V_pathC(g, beta=1.0, gamma=1.0, mu2=None):
    """
    Path C potential: V_TGP + symmetry-breaking stabilization.

    V_C(g) = (beta/3)*g^3 - (gamma/4)*g^4 + (mu^2/2)*(g-1)^2

    If mu2 is None, compute from beta, gamma:
        mu^2 = 2*(3*gamma - 2*beta)
    """
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)
    return V_tgp(g, beta, gamma) + 0.5 * mu2 * (g - 1.0)**2


def dV_pathC(g, beta=1.0, gamma=1.0, mu2=None):
    """dV_C/dg for Path C potential."""
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)
    return dV_tgp(g, beta, gamma) + mu2 * (g - 1.0)


def d2V_pathC(g, beta=1.0, gamma=1.0, mu2=None):
    """d^2V_C/dg^2 for Path C potential."""
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)
    return 2.0 * beta * g - 3.0 * gamma * g**2 + mu2


def mu2_self_consistent(beta=1.0, gamma=1.0):
    """
    Self-consistent stabilization strength mu^2.

    From linearized defect ODE at g=1:
      (nabla^2 - V_C''(1))*delta = 0

    We need V_C''(1) = m_sp^2 = 3*gamma - 2*beta  (Path B screening mass).

    V_TGP''(1) = 2*beta - 3*gamma
    V_sb''(1)  = mu^2
    V_C''(1) = (2*beta - 3*gamma) + mu^2 = 3*gamma - 2*beta

    => mu^2 = 2*(3*gamma - 2*beta)

    For beta=gamma=1: mu^2 = 2, m_sp = 1
    For beta=gamma: mu^2 = 2*gamma
    """
    m_sp2 = 3.0 * gamma - 2.0 * beta
    if m_sp2 < 0:
        raise ValueError(f"Tachyonic regime: 3*gamma - 2*beta = {m_sp2} < 0")
    return 2.0 * m_sp2


def effective_mass_squared(beta=1.0, gamma=1.0, mu2=None):
    """
    Effective mass squared at vacuum g=1.

    m_eff^2 = V_C''(1) = (2*beta - 3*gamma) + mu^2

    With self-consistent mu^2 = 2*(3*gamma - 2*beta):
        m_eff^2 = (2*beta - 3*gamma) + 2*(3*gamma - 2*beta)
                = (2*beta - 3*gamma) + (6*gamma - 4*beta)
                = 3*gamma - 2*beta  = m_sp^2 (Path B)
    """
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)
    return (2.0 * beta - 3.0 * gamma) + mu2


def screening_mass_pathC(beta=1.0, gamma=1.0, mu2=None):
    """m_sp = sqrt(V_C''(1)) = sqrt(m_eff^2)."""
    m2 = effective_mass_squared(beta, gamma, mu2)
    if m2 < 0:
        raise ValueError(f"Tachyonic: m_eff^2 = {m2} < 0 (mu^2 too small)")
    return np.sqrt(m2)


# -- Defect ODE solver (Path C) ------------------------------------------

def solve_defect_pathC(g0, beta=1.0, gamma=1.0, mu2=None,
                       r_max=60.0, kinetic="full",
                       n_eval=5000, rtol=1e-10, atol=1e-13):
    """
    Solve the spherically symmetric defect ODE with Path C potential.

    K(g)*[g'' + (2/r)*g'] + K'(g)/2 * (g')^2 = V_C'(g)

    where K(g) depends on kinetic mode:
      "full": K(g) = g^2,  K'(g) = 2*g
      "lpa":  K(g) = 1 + 4*ln(g),  K'(g) = 4/g

    Parameters
    ----------
    g0 : float
        Core value g(r=0) with g'(0)=0.
    beta, gamma : float
        TGP parameters.
    mu2 : float or None
        Stabilization strength. If None, self-consistent.
    r_max : float
        Outer radius of integration.
    kinetic : str
        "full" (K_sub = g^2, ghost-free) or "lpa" (logarithmic, has ghost).

    Returns
    -------
    dict with keys:
        r : radial grid
        g : field profile g(r)
        gp : derivative g'(r)
        C_eff : effective Yukawa coupling (from tail fit)
        m_sp : screening mass
        E_defect : integrated defect energy
        tail_quality : goodness of Yukawa tail fit
    """
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)

    m_eff2 = effective_mass_squared(beta, gamma, mu2)
    if m_eff2 > 0:
        m_sp = np.sqrt(m_eff2)
    else:
        m_sp = 0.0  # oscillatory regime, no Yukawa mass

    ALPHA = 2.0
    if kinetic == "lpa":
        g_ghost = np.exp(-1.0 / (2.0 * ALPHA))
    else:
        g_ghost = 0.0  # Full K_sub = g^2 is ghost-free

    def K(g):
        g = max(g, 1e-15)
        if kinetic == "full":
            return g**2
        else:
            return 1.0 + 2.0 * ALPHA * np.log(g)

    def dK(g):
        g = max(g, 1e-15)
        if kinetic == "full":
            return 2.0 * g
        else:
            return 2.0 * ALPHA / g

    def ode_rhs(r, y):
        g, gp = y[0], y[1]
        g = max(g, 1e-15)

        Kg = K(g)
        dKg = dK(g)
        driving = dV_pathC(g, beta, gamma, mu2)

        if r < 1e-10:
            # L'Hopital: g'' + 2*g'/r -> 3*g'' as r->0
            if abs(Kg) < 1e-12:
                gpp = 0.0
            else:
                gpp = driving / (3.0 * Kg)
            return [gp, gpp]

        if abs(Kg) < 1e-12:
            gpp = 0.0
        else:
            gpp = (driving - 0.5 * dKg * gp**2 - Kg * 2.0 * gp / r) / Kg

        return [gp, gpp]

    # Integration with ghost boundary check (only for LPA)
    events = []
    if kinetic == "lpa":
        def event_ghost(r, y):
            return y[0] - g_ghost - 1e-4
        event_ghost.terminal = True
        event_ghost.direction = -1
        events = [event_ghost]

    r_eval = np.linspace(1e-4, r_max, n_eval)
    sol = solve_ivp(ode_rhs, [1e-4, r_max], [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=rtol, atol=atol,
                    events=events if events else None,
                    max_step=0.1)

    r_arr = sol.t
    g_arr = sol.y[0]
    gp_arr = sol.y[1]

    # Fit Yukawa tail: (1-g)*r ~ C*exp(-m_sp*r) for large r
    C_eff, tail_quality = _fit_yukawa_tail(r_arr, g_arr, m_sp)

    # Compute defect energy
    E_defect = _compute_defect_energy(r_arr, g_arr, gp_arr,
                                       beta, gamma, mu2, kinetic)

    return {
        'r': r_arr,
        'g': g_arr,
        'gp': gp_arr,
        'C_eff': C_eff,
        'm_sp': m_sp,
        'mu2': mu2,
        'E_defect': E_defect,
        'tail_quality': tail_quality,
        'g0': g0,
        'kinetic': kinetic,
        'hit_ghost': len(sol.t_events[0]) > 0 if events and sol.t_events else False,
    }


def _fit_yukawa_tail(r_arr, g_arr, m_sp, r_start_frac=0.5, r_end_frac=0.9):
    """
    Fit the tail of the defect to Yukawa: delta(r) ~ C_eff * exp(-m_sp*r) / r.

    We linearize: (1 - g(r)) * r * exp(m_sp * r) -> C_eff for large r.

    Returns (C_eff, quality) where quality is the relative std of the estimate.
    """
    n = len(r_arr)
    i_start = int(n * r_start_frac)
    i_end = int(n * r_end_frac)
    if i_end - i_start < 10:
        return 0.0, float('inf')

    r_fit = r_arr[i_start:i_end]
    g_fit = g_arr[i_start:i_end]
    delta = 1.0 - g_fit

    # C_eff estimate: delta * r * exp(m_sp * r)
    # For Yukawa tail: delta = C * exp(-m_sp*r) / r
    # => delta * r * exp(m_sp * r) = C
    exponent = m_sp * r_fit
    # Guard against overflow
    max_exp = 500.0
    mask = exponent < max_exp
    if np.sum(mask) < 5:
        # Try with just the closest points
        mask = np.arange(len(r_fit)) < min(20, len(r_fit))

    r_m = r_fit[mask]
    delta_m = delta[mask]
    exp_m = np.exp(m_sp * r_m)

    C_estimates = delta_m * r_m * exp_m
    C_eff = np.median(C_estimates)
    if abs(C_eff) > 1e-30:
        quality = np.std(C_estimates) / abs(C_eff)
    else:
        quality = float('inf')

    return float(C_eff), float(quality)


def _compute_defect_energy(r_arr, g_arr, gp_arr,
                           beta, gamma, mu2, kinetic):
    """
    Compute defect energy by integration:

    E = 4*pi * integral_0^R_max { (1/2)*K(g)*(g')^2 + V_C(g) - V_C(1) } r^2 dr

    We subtract V_C(1) to measure energy relative to vacuum.
    """
    V_vac = V_pathC(1.0, beta, gamma, mu2)
    dr = np.diff(r_arr)

    if kinetic == "full":
        K_arr = g_arr**2
    else:
        K_arr = 1.0 + 4.0 * np.log(np.maximum(g_arr, 1e-15))

    V_arr = np.array([V_pathC(g, beta, gamma, mu2) for g in g_arr])

    # Energy density: (1/2)*K*(g')^2 + V_C(g) - V_C(1)
    rho = 0.5 * K_arr * gp_arr**2 + (V_arr - V_vac)

    # Trapezoidal integration of 4*pi*rho*r^2
    integrand = 4.0 * np.pi * rho * r_arr**2
    # np.trapz was removed in numpy 2.0; use np.trapezoid or fallback
    try:
        E = np.trapezoid(integrand, r_arr)
    except AttributeError:
        E = np.trapz(integrand, r_arr)

    return float(E)


# -- EFT projection: C_eff from classical defect + Yukawa Green's fn ------

def compute_C_eff_projection(r_arr, g_arr, m_sp):
    """
    Compute effective Yukawa coupling by projecting the classical defect
    profile onto the Yukawa Green's function.

    EFT logic:
      1. Classical defect has oscillatory tail (standard TGP)
      2. Loop corrections generate mass gap m_sp (Path B screening mass)
      3. Effective interaction is Yukawa with coupling:

         C_eff = 4*pi * integral_0^inf delta(r) * G_Y(r) * r^2 dr

      where delta(r) = 1 - g(r) is the defect deviation from vacuum,
      and G_Y(r) = exp(-m_sp*r) / (4*pi*r) is the Yukawa Green's function.

    Simplifying:
         C_eff = integral_0^inf delta(r) * exp(-m_sp*r) * r dr

    This DERIVES C_eff from the classical defect shape, connecting
    Path B (phenomenological C_i) to Path C (topological defect).

    Parameters
    ----------
    r_arr : array
        Radial grid from defect ODE solution.
    g_arr : array
        Field profile g(r).
    m_sp : float
        Screening mass (from Path B or self-consistent).

    Returns
    -------
    C_eff : float
        Effective Yukawa coupling strength.
    """
    delta = 1.0 - g_arr
    integrand = delta * np.exp(-m_sp * r_arr) * r_arr

    try:
        C_eff = np.trapezoid(integrand, r_arr)
    except AttributeError:
        C_eff = np.trapz(integrand, r_arr)

    return float(C_eff)


def solve_defect_standard(g0, beta=1.0, gamma=1.0,
                          r_max=60.0, kinetic="full",
                          n_eval=5000, rtol=1e-10, atol=1e-13):
    """
    Solve standard TGP defect ODE (no stabilization, mu2=0) and compute
    C_eff via EFT projection onto Yukawa Green's function.

    This is the physically correct approach:
      - Classical ODE uses standard V_TGP (oscillatory tail)
      - C_eff comes from projecting delta(r) onto exp(-m_sp*r)/r
      - m_sp is the Path B screening mass sqrt(3*gamma - 2*beta)

    Returns dict with all defect data plus:
      C_eff_proj : projection-based coupling
      m_sp_pathB : Path B screening mass used for projection
      tail_type  : 'oscillatory' (always, for standard TGP)
    """
    # Solve with mu2=0 (standard TGP, no stabilization)
    res = solve_defect_pathC(g0, beta, gamma, mu2=0.0,
                             r_max=r_max, kinetic=kinetic,
                             n_eval=n_eval, rtol=rtol, atol=atol)

    # Path B screening mass for projection
    m_sp2 = 3.0 * gamma - 2.0 * beta
    if m_sp2 <= 0:
        raise ValueError(f"Tachyonic: 3*gamma - 2*beta = {m_sp2} <= 0")
    m_sp_B = np.sqrt(m_sp2)

    # EFT projection
    C_eff_proj = compute_C_eff_projection(res['r'], res['g'], m_sp_B)

    # Oscillatory tail fit for comparison
    osc_amp, osc_phase, osc_resid = fit_oscillatory_tail(res['r'], res['g'])

    res['C_eff_proj'] = C_eff_proj
    res['m_sp_pathB'] = m_sp_B
    res['tail_type'] = 'oscillatory'
    res['osc_amplitude'] = osc_amp
    res['osc_phase'] = osc_phase
    res['osc_residual'] = osc_resid

    return res


def demonstrate_stabilization_obstacle(beta=1.0, gamma=1.0, n_points=200):
    """
    Show WHY naive quadratic stabilization destroys classical defects.

    V_C'(g) = beta*g^2 - gamma*g^3 + mu^2*(g-1)
            = (1-g)*(g^2 - mu^2)   [for beta=gamma]

    For mu^2 > 1 and g in (0,1): g^2 < 1 < mu^2, so (g^2 - mu^2) < 0
    and (1-g) > 0, giving V_C'(g) < 0 for ALL g in (0,1).

    This means no turning point exists -- the field is driven away from
    vacuum everywhere, destroying defect solutions.

    Returns dict with analysis data.
    """
    mu2 = mu2_self_consistent(beta, gamma)
    g_grid = np.linspace(0.01, 0.999, n_points)

    dV_standard = np.array([dV_tgp(g, beta, gamma) for g in g_grid])
    dV_stabilized = np.array([dV_pathC(g, beta, gamma, mu2) for g in g_grid])

    # Check: is dV_stabilized ever positive in (0,1)?
    has_positive = np.any(dV_stabilized > 0)
    # For defect to exist, we need dV_C'(g) > 0 somewhere in (0,1)
    # so that the field can be driven TOWARD vacuum

    # Count zero crossings of dV_standard in (0,1)
    sign_changes_std = np.sum(np.diff(np.sign(dV_standard)) != 0)
    sign_changes_stab = np.sum(np.diff(np.sign(dV_stabilized)) != 0)

    return {
        'g_grid': g_grid,
        'dV_standard': dV_standard,
        'dV_stabilized': dV_stabilized,
        'mu2': mu2,
        'has_positive_dV_C': has_positive,
        'sign_changes_standard': sign_changes_std,
        'sign_changes_stabilized': sign_changes_stab,
        'defect_possible_standard': sign_changes_std > 0,
        'defect_impossible_stabilized': not has_positive,
    }


# -- Oscillatory vs Yukawa tail comparison --------------------------------

def fit_oscillatory_tail(r_arr, g_arr, r_start_frac=0.5, r_end_frac=0.9):
    """
    Fit tail to oscillatory form: delta(r) = [A*cos(r) + B*sin(r)] / r.

    Returns (amplitude, phase, residual_rms).
    """
    n = len(r_arr)
    i_start = int(n * r_start_frac)
    i_end = int(n * r_end_frac)
    if i_end - i_start < 10:
        return 0.0, 0.0, float('inf')

    r_fit = r_arr[i_start:i_end]
    g_fit = g_arr[i_start:i_end]
    delta_r = (1.0 - g_fit) * r_fit  # delta * r

    # Fit: delta*r = a*cos(r) + b*sin(r)
    A_mat = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, delta_r, rcond=None)
        a, b = coeffs
        amplitude = np.sqrt(a**2 + b**2)
        phase = np.arctan2(b, a)
        predicted = A_mat @ coeffs
        residual = np.sqrt(np.mean((delta_r - predicted)**2))
        return float(amplitude), float(phase), float(residual)
    except Exception:
        return 0.0, 0.0, float('inf')


def compare_tail_types(r_arr, g_arr, m_sp, r_start_frac=0.5, r_end_frac=0.9):
    """
    Compare oscillatory vs Yukawa tail fits.

    Returns dict with RMSE for each model and which is better.
    """
    n = len(r_arr)
    i_start = int(n * r_start_frac)
    i_end = int(n * r_end_frac)

    r_fit = r_arr[i_start:i_end]
    delta = 1.0 - g_arr[i_start:i_end]

    # Oscillatory: delta = [a*cos(r) + b*sin(r)] / r
    A_osc = np.column_stack([np.cos(r_fit) / r_fit, np.sin(r_fit) / r_fit])
    try:
        c_osc, _, _, _ = np.linalg.lstsq(A_osc, delta, rcond=None)
        pred_osc = A_osc @ c_osc
        rmse_osc = np.sqrt(np.mean((delta - pred_osc)**2))
    except Exception:
        rmse_osc = float('inf')

    # Yukawa: delta = C * exp(-m_sp*r) / r
    yuk_basis = np.exp(-m_sp * r_fit) / r_fit
    c_yuk = np.dot(delta, yuk_basis) / np.dot(yuk_basis, yuk_basis)
    pred_yuk = c_yuk * yuk_basis
    rmse_yuk = np.sqrt(np.mean((delta - pred_yuk)**2))

    better = "yukawa" if rmse_yuk < rmse_osc else "oscillatory"

    return {
        'rmse_oscillatory': rmse_osc,
        'rmse_yukawa': rmse_yuk,
        'C_yukawa': c_yuk,
        'better': better,
        'ratio': rmse_osc / max(rmse_yuk, 1e-30),
    }


# -- Scan: C_eff vs g0 for Path C defects --------------------------------

def scan_defect_coupling(beta=1.0, gamma=1.0, mu2=None,
                         g0_min=0.50, g0_max=0.99, n_g0=20,
                         kinetic="full", r_max=60.0):
    """
    Scan defect core values g0, compute C_eff(g0) and E_defect(g0).

    Returns list of result dicts.
    """
    if mu2 is None:
        mu2 = mu2_self_consistent(beta, gamma)

    g0_vals = np.linspace(g0_min, g0_max, n_g0)
    results = []

    for g0 in g0_vals:
        try:
            res = solve_defect_pathC(g0, beta, gamma, mu2,
                                     r_max=r_max, kinetic=kinetic)
            results.append(res)
        except Exception as e:
            results.append({
                'g0': g0, 'C_eff': None, 'E_defect': None,
                'error': str(e),
            })

    return results


# -- Self-consistency check -----------------------------------------------

def verify_pathC_consistency(beta=1.0, gamma=1.0):
    """
    Verify that Path C gives the same m_sp as Path B.

    Path B: m_sp = sqrt(3*gamma - 2*beta)
    Path C: m_sp = sqrt(V_C''(1)) = sqrt((2b-3g)+mu^2)
            with mu^2 = 2*(3*gamma - 2*beta):
            m_sp = sqrt((2b-3g) + 2*(3g-2b)) = sqrt(3g-2b)  (SAME!)

    Returns dict with comparison.
    """
    from .tgp_field import screening_mass as m_sp_pathB

    m_B = m_sp_pathB(beta, gamma)
    mu2 = mu2_self_consistent(beta, gamma)
    m_C = screening_mass_pathC(beta, gamma, mu2)

    return {
        'm_sp_pathB': m_B,
        'm_sp_pathC': m_C,
        'mu2': mu2,
        'match': abs(m_B - m_C) < 1e-12,
        'relative_error': abs(m_B - m_C) / max(abs(m_B), 1e-30),
    }
