"""
dynamics_v2.py -- TGP n-body dynamics with analytical pairwise forces
=====================================================================

CONVENTION (AXIOM):
  Inertial mass m_i is proportional to source strength C_i:
      m_i = C_i      (equivalence of source and inertia in TGP)

  All functions in this module return FORCES (not accelerations):
      F_i = -nabla_i V_total

  The integrators expect ACCELERATIONS (a_i = F_i / C_i).
  Use dynamics_backends.build_tgp_integration_pair() which wraps
  force functions into acceleration functions (F / C[:, None])
  before passing them to leapfrog_integrate / rk45_integrate.

STATUS OF RESULTS:
  EXACT (pairwise sector):
    - analytical_forces_tgp_pairwise(): closed-form gradient of V_2
    - pairwise_tgp_force_jacobian(): closed-form ∂F/∂x for those forces
    - potential_tgp_pairwise(): closed-form pairwise potential
    - forces_newton(): closed-form Newtonian forces

  APPROXIMATE (irreducible 3-body sector):
    - three_body_forces_approximate(): uses analytic gradient of
      triple overlap integral. EXACT for Coulomb profiles (1/r).
      APPROXIMATE for Yukawa profiles (saddle-point formula with
      WRONG exponential exponent — error 160-770% for m*d in [0.5,6]).
      Use three_body_force_exact.three_body_forces_exact() for exact
      Yukawa 3-body forces (Feynman 2D integral, O(n_quad^2) cost).
    - three_body_potential_approximate(): same approximation level.

  NUMERICAL (full dynamics):
    - leapfrog_integrate(), rk45_integrate(): numerical time integration
      of the equations of motion with softened forces.
    - Softening epsilon is a NUMERICAL REGULATOR, not a physical cutoff.
      Results should be checked for convergence as epsilon -> 0.

TGP PAIRWISE POTENTIAL (exact, between bodies i and j):
    V_2(r) = -4*pi*Ci*Cj/r + 8*pi*beta*Ci*Cj/r^2
             - 12*pi*gamma*Ci*Cj*(Ci+Cj)/r^3

    F_ij = -dV_2/dr * r_hat  (force on i from j, pointing i -> j)

    dV_2/dr = 4*pi*Ci*Cj/r^2 - 16*pi*beta*Ci*Cj/r^3
              + 36*pi*gamma*Ci*Cj*(Ci+Cj)/r^4

Sign convention: positive dV/dr at the attractive side means
F_ij points toward j (attraction).
"""

import numpy as np
from scipy.integrate import solve_ivp
from itertools import combinations


# ── Irreducible 3-body forces (ANALYTICAL) ───────────────────────────

def three_body_forces_approximate(positions, C_values, beta, gamma=None,
                                  softening=1e-6, use_yukawa=False):
    """
    Analytical 3-body forces from the quartic TGP nonlinearity.

    DERIVATION:
    -----------
    The irreducible 3-body energy for triplet (i,j,k) is:

        V_3(i,j,k) = -6*gamma * C_i*C_j*C_k * I_triple

    where I_triple is the triple overlap integral of the field profiles.

    For Coulomb-like (unscreened) profiles:
        I_triple = 8*pi^2 / P,  where P = d_ij + d_ik + d_jk

    The force on mass i is F_i = -dV_3/dx_i:

        dV_3/dx_i = -6*gamma * C_i*C_j*C_k * dI/dx_i

        dI/dx_i = -8*pi^2 / P^2 * dP/dx_i

        dP/dx_i = d(d_ij)/dx_i + d(d_ik)/dx_i
                = (x_i - x_j)/d_ij + (x_i - x_k)/d_ik

    Therefore:
        F_i = -dV_3/dx_i
            = -6*gamma * Ci*Cj*Ck * 8*pi^2 / P^2 * [(x_i-x_j)/d_ij + (x_i-x_k)/d_ik]

    This is a SUM OF TWO UNIT VECTORS, weighted by the coupling strength
    and the inverse square of the perimeter. Completely analytical!

    For Yukawa-screened profiles (m_sp = sqrt(gamma)):
        I_triple ~ (2*pi)^(3/2) / m^2 * exp(-m*s) / s,  s = P/2

        dI/dx_i = I * (-m - 1/s) * (1/2) * dP/dx_i
                = I * (-m*s - 1) / (s * P) * dP/dx_i  (since ds = dP/2)

    Parameters
    ----------
    positions : (n, 3) ndarray
    C_values : (n,) ndarray
    beta, gamma : float
    softening : float
    use_yukawa : bool
        If True, use Yukawa-screened overlap (more accurate for TGP).
        If False, use Coulomb overlap (simpler, valid when m*d << 1).

    Returns
    -------
    forces_3body : (n, 3) ndarray
        Force from 3-body terms on each mass (NOT acceleration).

    STATUS: EXACT for Coulomb profiles (I_triple = 8*pi^2/P).
            APPROXIMATE for Yukawa profiles: saddle-point formula with
            INCORRECT exponential exponent (exp(-1.5*t) used instead of
            exp(-sqrt(3)*t) for equilateral). Error vs Feynman integral:
              t = m*d = 0.5 -> +60%,  t=1.0 -> +86%,  t=2.0 -> +148%,
              t=3.0 -> +208%,  t=5.0 -> +329%.
            Valid only for t = m_sp * d > ~5 (error < 5%).
            For exact Yukawa 3-body forces use:
                three_body_force_exact.three_body_forces_exact()
            (Feynman 2D integral, verified in ex6/ex18; zapis sesji: _archiwum_docs/WYNIKI_SESJI_2026_03_21.md).
    """
    if gamma is None:
        gamma = beta

    n = len(C_values)
    forces = np.zeros((n, 3))

    if n < 3:
        return forces

    if use_yukawa:
        m_sp = np.sqrt(max(3.0 * gamma - 2.0 * beta, 1e-20))

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_values[i], C_values[j], C_values[k]

        rij = positions[j] - positions[i]
        rik = positions[k] - positions[i]
        rjk = positions[k] - positions[j]

        dij = np.sqrt(np.dot(rij, rij) + softening**2)
        dik = np.sqrt(np.dot(rik, rik) + softening**2)
        djk = np.sqrt(np.dot(rjk, rjk) + softening**2)

        P = dij + dik + djk  # perimeter
        # Force coupling: (6*gamma - 2*beta) from corrected potential vertex
        coupling = (6.0 * gamma - 2.0 * beta) * Ci * Cj * Ck

        if use_yukawa:
            s = P / 2.0
            exp_ms = np.exp(-m_sp * s)
            I_val = (2.0 * np.pi) ** 1.5 / m_sp**2 * exp_ms / s
            # dI/dP = I * (-m/2 - 1/(2s)) = I * (-m*s - 1) / (2*s)
            dI_dP = I_val * (-m_sp * s - 1.0) / (2.0 * s)
        else:
            I_val = 8.0 * np.pi**2 / P
            dI_dP = -8.0 * np.pi**2 / P**2

        # Gradient of perimeter w.r.t. each mass position
        # dP/dx_i = (x_i - x_j)/d_ij + (x_i - x_k)/d_ik
        # dP/dx_j = (x_j - x_i)/d_ij + (x_j - x_k)/d_jk
        # dP/dx_k = (x_k - x_i)/d_ik + (x_k - x_j)/d_jk

        dP_dxi = -rij / dij + (-rik) / dik  # (x_i-x_j)/d_ij + (x_i-x_k)/d_ik
        dP_dxj = rij / dij + (-rjk) / djk
        dP_dxk = rik / dik + rjk / djk

        # V_3 = -coupling * I_val
        # F_i = -dV_3/dx_i = coupling * dI/dP * dP/dx_i
        # (the minus from -coupling and the minus from F=-dV cancel,
        #  and dI/dP is already negative for Coulomb)
        # Actually: F_i = -d(-coupling*I)/dx_i = coupling * dI/dx_i
        #         = coupling * dI/dP * dP/dx_i

        forces[i] += coupling * dI_dP * dP_dxi
        forces[j] += coupling * dI_dP * dP_dxj
        forces[k] += coupling * dI_dP * dP_dxk

    return forces


def three_body_potential(positions, C_values, beta, gamma=None,
                          softening=1e-6, use_yukawa=False):
    """Total irreducible 3-body potential energy.

    V_3_total = sum_{i<j<k} V_3(i,j,k)
             = sum_{i<j<k} -6*gamma * Ci*Cj*Ck * I_triple(d_ij, d_ik, d_jk)
    """
    if gamma is None:
        gamma = beta

    n = len(C_values)
    if n < 3:
        return 0.0

    V3 = 0.0

    if use_yukawa:
        m_sp = np.sqrt(max(3.0 * gamma - 2.0 * beta, 1e-20))

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_values[i], C_values[j], C_values[k]

        rij = positions[j] - positions[i]
        rik = positions[k] - positions[i]
        rjk = positions[k] - positions[j]

        dij = np.sqrt(np.dot(rij, rij) + softening**2)
        dik = np.sqrt(np.dot(rik, rik) + softening**2)
        djk = np.sqrt(np.dot(rjk, rjk) + softening**2)

        P = dij + dik + djk

        if use_yukawa:
            s = P / 2.0
            I_val = (2.0 * np.pi) ** 1.5 / m_sp**2 * np.exp(-m_sp * s) / s
        else:
            I_val = 8.0 * np.pi**2 / P

        V3 += -6.0 * gamma * Ci * Cj * Ck * I_val

    return V3


def analytical_forces_tgp_pairwise(positions, C_values, beta, gamma=None, softening=1e-6):
    """
    Exact pairwise TGP forces on all bodies (closed-form gradient of V_2).

    Returns F_i = -nabla_i V_2_total, which is a FORCE (not acceleration).
    To get acceleration, divide by m_i = C_i.

    Parameters
    ----------
    positions : (n, 3) ndarray
    C_values : (n,) ndarray — source strengths (= inertial masses)
    beta : float
    gamma : float or None (defaults to beta for vacuum condition)
    softening : float — numerical regulator (not physical)

    Returns
    -------
    forces : (n, 3) ndarray — force on each body (NOT acceleration)
    """
    if gamma is None:
        gamma = beta

    n = len(C_values)
    forces = np.zeros((n, 3))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            r_vec = positions[j] - positions[i]  # vector from i to j
            r2 = np.dot(r_vec, r_vec) + softening**2  # softened r^2
            r = np.sqrt(r2)
            r_hat = r_vec / r

            Ci, Cj = C_values[i], C_values[j]

            # Force on i from j: F_i = -dV/dx_i = (dV/dr) * r_hat
            # (because dr/dx_i = -r_hat, so -dV/dx_i = +dV/dr * r_hat)
            #
            # dV/dr = 4*pi*Ci*Cj/r^2 - 16*pi*beta*Ci*Cj/r^3
            #         + 36*pi*gamma*Ci*Cj*(Ci+Cj)/r^4
            dVdr = (4.0 * np.pi * Ci * Cj / r2
                    - 16.0 * np.pi * beta * Ci * Cj / (r2 * r)
                    + 36.0 * np.pi * gamma * Ci * Cj * (Ci + Cj) / (r2 * r2))

            # F_i = -dV/dx_i = (dV/dr) * r_hat
            # (because dr/dx_i = -r_hat, so -dV/dx_i = dV/dr * r_hat)
            forces[i] += dVdr * r_hat

    return forces


def pairwise_tgp_force_jacobian(
    positions: np.ndarray,
    C_values: np.ndarray,
    beta: float,
    gamma: float | None = None,
    softening: float = 1e-6,
) -> np.ndarray:
    """
    Jacobian ``∂F/∂x`` for ``analytical_forces_tgp_pairwise`` (forces, not acc.).

    Layout: row index ``3*i+α`` is component ``F_{i,α}``, column ``3*j+β`` is
    ``∂/∂x_{j,β}``, matching ``lyapunov.acceleration_jacobian_newton_softened``.
    """
    if gamma is None:
        gamma = beta

    pos = np.asarray(positions, dtype=float)
    C = np.asarray(C_values, dtype=float)
    n = len(C)
    m = 3 * n
    dFdq = np.zeros((m, m))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            u = pos[j] - pos[i]
            r2 = float(np.dot(u, u) + softening**2)
            r = float(np.sqrt(r2))
            if r < 1e-15:
                continue

            Ci, Cj = C[i], C[j]
            h = (
                4.0 * np.pi * Ci * Cj / r2
                - 16.0 * np.pi * beta * Ci * Cj / (r2 * r)
                + 36.0 * np.pi * gamma * Ci * Cj * (Ci + Cj) / (r2 * r2)
            )
            h2 = (
                -8.0 * np.pi * Ci * Cj / (r2 * r)
                + 48.0 * np.pi * beta * Ci * Cj / (r2 * r2)
                - 144.0 * np.pi * gamma * Ci * Cj * (Ci + Cj) / (r2 * r2 * r)
            )
            coeff_outer = (h2 * r - h) / (r**3)
            coeff_id = h / r
            outer = np.outer(u, u)
            # Ordered pair (i,j): only row F_i gets this pair's contribution
            # (mirror (j,i) iteration handles F_j).
            dFdq[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] += (
                -coeff_id * np.eye(3) - coeff_outer * outer
            )
            dFdq[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] += (
                coeff_id * np.eye(3) + coeff_outer * outer
            )

    return dFdq


def forces_newton(positions, C_values, G=1.0, softening=1e-6):
    """
    Compute Newtonian gravitational forces analytically.

    Parameters
    ----------
    positions : (n, 3) ndarray
    C_values : (n,) ndarray  (mass ~ C)
    G : float
    softening : float

    Returns
    -------
    forces : (n, 3) ndarray
    """
    n = len(C_values)
    forces = np.zeros((n, 3))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r_vec = positions[j] - positions[i]
            r2 = np.dot(r_vec, r_vec) + softening**2
            r = np.sqrt(r2)
            # a_i = G * C_j * r_hat / r^2
            forces[i] += G * C_values[j] * r_vec / (r2 * r)

    return forces


def potential_tgp(positions, C_values, beta, gamma=None, softening=1e-6):
    """Total TGP pairwise potential energy."""
    if gamma is None:
        gamma = beta
    n = len(C_values)
    V = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = positions[j] - positions[i]
            r2 = np.dot(r_vec, r_vec) + softening**2
            r = np.sqrt(r2)
            Ci, Cj = C_values[i], C_values[j]
            V += (-4.0 * np.pi * Ci * Cj / r
                  + 8.0 * np.pi * beta * Ci * Cj / r2
                  - 12.0 * np.pi * gamma * Ci * Cj * (Ci + Cj) / (r2 * r))
    return V


def potential_newton(positions, C_values, G=1.0, softening=1e-6):
    """Total Newtonian potential energy."""
    n = len(C_values)
    V = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = positions[j] - positions[i]
            r2 = np.dot(r_vec, r_vec) + softening**2
            r = np.sqrt(r2)
            V += -G * C_values[i] * C_values[j] / r
    return V


# ── Full TGP forces and potential (2-body + 3-body) ──────────────────────

def full_forces_tgp(positions, C_values, beta, gamma=None,
                    softening=1e-6, include_3body=True, use_yukawa=False):
    """
    TGP forces: exact pairwise + approximate irreducible 3-body.

    Returns F_i = -nabla_i (V_2 + V_3), a FORCE vector.
    The 3-body sector uses the saddle-point approximation for the
    Yukawa triple overlap integral (see three_body_forces_approximate).
    """
    forces = analytical_forces_tgp_pairwise(positions, C_values, beta, gamma, softening)

    if include_3body and len(C_values) >= 3:
        forces += three_body_forces_approximate(positions, C_values, beta, gamma,
                                                softening, use_yukawa)

    return forces


def full_potential_tgp(positions, C_values, beta, gamma=None,
                       softening=1e-6, include_3body=True, use_yukawa=False):
    """
    Complete TGP potential: pairwise + irreducible 3-body.
    """
    V = potential_tgp(positions, C_values, beta, gamma, softening)

    if include_3body and len(C_values) >= 3:
        V += three_body_potential(positions, C_values, beta, gamma,
                                  softening, use_yukawa)

    return V


# ── Integrators ──────────────────────────────────────────────────────────

def leapfrog_integrate(positions, velocities, C_values,
                       force_func, potential_func,
                       t_span=(0, 100), dt=0.01, save_every=10):
    """
    Symplectic leapfrog (Stormer-Verlet) integrator.

    Preserves phase-space volume exactly, giving much better
    long-term energy conservation than Runge-Kutta methods.

    Parameters
    ----------
    positions : (n, 3) ndarray
    velocities : (n, 3) ndarray
    C_values : (n,) ndarray
    force_func : callable
        Returns acceleration (n, 3) array.
        NOTE: if using TGP force functions from this module (which return
        forces, not accelerations), the caller must wrap them to divide
        by C_i: lambda p, c: analytical_forces_tgp_pairwise(p, c, ...) / c[:, None]
    potential_func : callable
        V = potential_func(positions, C_values) -> float
    t_span : (t0, tf)
    dt : float
        Time step.
    save_every : int
        Save state every this many steps.

    Returns
    -------
    result : dict
    """
    t0, tf = t_span
    n_steps = int((tf - t0) / dt)
    n_save = n_steps // save_every + 1

    n = len(C_values)
    pos = positions.copy()
    vel = velocities.copy()

    # Storage
    t_arr = np.zeros(n_save)
    pos_arr = np.zeros((n_save, n, 3))
    vel_arr = np.zeros((n_save, n, 3))
    energy_arr = np.zeros(n_save)

    def compute_energy(p, v):
        # KE = sum 0.5 * m_i * |v_i|^2, with m_i = C_i (TGP axiom)
        T = 0.5 * np.sum(C_values[:, None] * v**2)
        V = potential_func(p, C_values)
        return T + V

    # Save initial state
    save_idx = 0
    t_arr[0] = t0
    pos_arr[0] = pos.copy()
    vel_arr[0] = vel.copy()
    energy_arr[0] = compute_energy(pos, vel)
    save_idx = 1

    # Initial half-step for velocity
    acc = force_func(pos, C_values)

    t = t0
    for step in range(n_steps):
        # Leapfrog: kick-drift-kick
        # Half kick
        vel += 0.5 * dt * acc

        # Full drift
        pos += dt * vel

        # New acceleration
        acc = force_func(pos, C_values)

        # Half kick
        vel += 0.5 * dt * acc

        t += dt

        # Save
        if (step + 1) % save_every == 0 and save_idx < n_save:
            t_arr[save_idx] = t
            pos_arr[save_idx] = pos.copy()
            vel_arr[save_idx] = vel.copy()
            energy_arr[save_idx] = compute_energy(pos, vel)
            save_idx += 1

    # Trim
    t_arr = t_arr[:save_idx]
    pos_arr = pos_arr[:save_idx]
    vel_arr = vel_arr[:save_idx]
    energy_arr = energy_arr[:save_idx]

    E0 = energy_arr[0] if abs(energy_arr[0]) > 1e-20 else 1.0

    return {
        't': t_arr,
        'positions': pos_arr,
        'velocities': vel_arr,
        'energy': energy_arr,
        'energy_error': (energy_arr - energy_arr[0]) / abs(E0),
        'success': True,
        'n_steps': n_steps,
    }


def rk45_integrate(positions, velocities, C_values,
                   force_func, potential_func,
                   t_span=(0, 100), n_output=1000,
                   rtol=1e-10, atol=1e-12,
                   quiet: bool = False):
    """
    High-order adaptive Runge-Kutta integrator (DOP853).

    Parameters
    ----------
    force_func : callable
        Returns acceleration (n, 3) array.
        NOTE: if using TGP force functions from this module (which return
        forces, not accelerations), the caller must wrap them to divide
        by C_i: lambda p, c: analytical_forces_tgp_pairwise(p, c, ...) / c[:, None]
    quiet : bool
        If True, do not print a warning when ``solve_ivp`` reports failure
        (callers may still inspect ``success`` / ``message`` in the result).
    """
    n = len(C_values)

    def ode_rhs(t, state):
        pos = state[:3*n].reshape(n, 3)
        vel = state[3*n:].reshape(n, 3)
        acc = force_func(pos, C_values)
        return np.concatenate([vel.flatten(), acc.flatten()])

    y0 = np.concatenate([positions.flatten(), velocities.flatten()])
    t_eval = np.linspace(t_span[0], t_span[1], n_output)

    sol = solve_ivp(ode_rhs, t_span, y0, t_eval=t_eval,
                    method='DOP853', rtol=rtol, atol=atol)

    if not sol.success and not quiet:
        print("WARNING: Integration failed: %s" % sol.message)

    n_t = len(sol.t)
    pos_traj = sol.y[:3*n].T.reshape(n_t, n, 3)
    vel_traj = sol.y[3*n:].T.reshape(n_t, n, 3)

    energy = np.zeros(n_t)
    for k in range(n_t):
        # KE = sum 0.5 * m_i * |v_i|^2, with m_i = C_i (TGP axiom)
        T = 0.5 * np.sum(C_values[:, None] * vel_traj[k]**2)
        V = potential_func(pos_traj[k], C_values)
        energy[k] = T + V

    E0 = energy[0] if abs(energy[0]) > 1e-20 else 1.0

    return {
        't': sol.t,
        'positions': pos_traj,
        'velocities': vel_traj,
        'energy': energy,
        'energy_error': (energy - energy[0]) / abs(E0),
        'success': sol.success,
        'message': getattr(sol, 'message', ''),
    }
