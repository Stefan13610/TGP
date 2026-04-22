"""
equilibria.py  --  TGP n-body: equilibrium finder
===================================================
Find static equilibrium configurations of n bodies interacting
via the TGP effective potential.

STATUS LEVELS:
  EXACT:        equilateral_pairwise_equilibrium() — closed-form polynomial
                from 2-body sector. Condition: d^2 - 4*beta*d + 18*beta*C = 0.
  APPROXIMATE:  equilateral_equilibrium_corrected() — first-order perturbative
                correction including irreducible 3-body energy V_3.
  NUMERICAL:    collinear_pairwise_equilibrium_numerical() — root scan.
                ngon_pairwise_equilibrium_numerical() — root scan.
                numerical_equilibrium() — general optimizer.

KEY PHYSICS:
In Newtonian gravity, Earnshaw's theorem forbids static equilibria
of point masses in free space. TGP breaks this because the repulsive
regime (1/d^2 term) provides a restoring force.

IMPORTANT: All equilibria computed here are in the PAIRWISE sector
unless explicitly stated otherwise. The irreducible 3-body corrections
shift the equilibrium positions by O(C) (see equilateral_equilibrium_corrected).
"""

import numpy as np
from scipy.optimize import minimize, root_scalar


def equilateral_pairwise_equilibrium(C, beta, gamma=None):
    """
    Pairwise-only equilibrium for equilateral triangle of 3 equal masses.

    STATUS: EXACT (closed-form polynomial from 2-body sector only).
    Does NOT include irreducible 3-body corrections.
    For corrected equilibrium, see equilateral_equilibrium_corrected().

    By symmetry, the total energy depends only on the side length d:
        V(d) = 3 * V_2(d) + V_3(d, d, d)

    The 2-body potential (for equal masses, beta=gamma):
        V_2(d) = -4*pi*C^2/d + 8*pi*beta*C^2/d^2 - 24*pi*beta*C^3/d^3

    The 3-body potential (from quartic coupling, Coulomb limit m_sp->0):
        V_3(d,d,d) = -16*pi^2*gamma*C^3 / d   [skalowanie 1/d, nie 1/d^2!]
        Wyprowadzenie: V_3 = -6*gamma*C^3 * I_triple, I = 8*pi^2/(3d) dla equilateral
        (see three_body_terms.py; TeX: dodatekD_trojcialowe.tex w korzeniu repo)

    Total: V(d) = -12*pi*C^2/d + 24*pi*beta*C^2/d^2 - 72*pi*beta*C^3/d^3
                  + 3-body corrections

    Equilibrium: dV/dd = 0

    For 2-body only (ignoring V_3):
        12*pi*C^2/d^2 - 48*pi*beta*C^2/d^3 + 216*pi*beta*C^3/d^4 = 0
        Multiply by d^4/(12*pi*C^2):
        d^2 - 4*beta*d + 18*beta*C = 0

    This is EXACTLY the same equation as the 2-body force zeros!
    (Proposition prop:trzy-rezimy-beta-gamma, sek03)

    Solutions: d = 2*beta ± sqrt(4*beta^2 - 18*beta*C)
    Exist when beta > 9*C/2.

    Returns
    -------
    dict with keys:
        'd_rep': repulsive equilibrium (unstable in 2-body, but may be
                 stable in 3-body due to additional restoring forces)
        'd_well': confinement equilibrium
        'exists': bool
        'polynomial_coefficients': for the equilibrium equation
    """
    if gamma is None:
        gamma = beta  # vacuum condition

    # 2-body equilibrium polynomial: d^2 - 4*beta*d + 18*gamma*C = 0
    discriminant = 4 * beta**2 - 18 * gamma * C
    if discriminant < 0:
        return {
            'exists': False,
            'd_rep': None,
            'd_well': None,
            'discriminant': discriminant,
            'critical_ratio': beta / C if C > 0 else np.inf,
            'message': f"No equilibrium: beta/C = {beta/C:.2f} < 4.5 (critical)"
        }

    sqrt_disc = np.sqrt(discriminant)
    d_rep = 2 * beta - sqrt_disc
    d_well = 2 * beta + sqrt_disc

    if d_rep <= 0:
        return {
            'exists': False,
            'd_rep': None,
            'd_well': d_well,
            'discriminant': discriminant,
            'message': "d_rep <= 0 (unphysical)"
        }

    return {
        'exists': True,
        'd_rep': d_rep,
        'd_well': d_well,
        'discriminant': discriminant,
        'critical_ratio': beta / C,
        'polynomial_coefficients': [1, -4*beta, 18*beta*C],
        'message': f"Two equilibria at d_rep={d_rep:.4f}, d_well={d_well:.4f}"
    }


def equilateral_equilibrium_corrected(C, beta, gamma=None, use_yukawa=True):
    """
    Perturbatively corrected equilibrium including 3-body effects.

    STATUS: APPROXIMATE (first-order perturbation theory).

    Method:
      1. Start from pairwise equilibrium d0 (exact, from polynomial).
      2. Compute 3-body force at d0: F_3 = -dV_3/dd |_{d=d0}
      3. Compute pairwise Hessian at d0: H_2 = d^2 V_2/dd^2 |_{d=d0}
      4. Corrected equilibrium: d* = d0 - F_3 / H_2

    The expansion parameter is epsilon ~ C (ratio |V_3/V_2|).

    Returns
    -------
    dict with:
        'd0_rep', 'd0_well': pairwise equilibria (exact)
        'd_corrected_rep', 'd_corrected_well': with 3-body correction
        'delta_d_rep', 'delta_d_well': perturbative shifts
        'exists': bool
    """
    if gamma is None:
        gamma = beta

    # Step 1: pairwise equilibrium
    pw = equilateral_pairwise_equilibrium(C, beta, gamma)
    if not pw['exists']:
        return {
            'exists': False,
            'message': pw.get('message', 'No pairwise equilibrium')
        }

    result = {
        'exists': True,
        'd0_rep': pw['d_rep'],
        'd0_well': pw['d_well'],
    }

    # Yukawa mass for 3-body screening
    m_sp = np.sqrt(max(3.0 * gamma - 2.0 * beta, 1e-20)) if use_yukawa else 0.0

    for label, d0 in [('rep', pw['d_rep']), ('well', pw['d_well'])]:
        if d0 is None or d0 <= 0:
            result['d_corrected_' + label] = None
            result['delta_d_' + label] = None
            continue

        # Step 2: 3-body force at d0 (equilateral, all sides = d0)
        # V_3 = -6*gamma*C^3 * I_triple(P), P = 3*d0
        # dV_3/dd = -6*gamma*C^3 * dI/dP * dP/dd, dP/dd = 3
        P = 3.0 * d0
        if use_yukawa:
            s = P / 2.0
            I_val = (2.0 * np.pi)**1.5 / m_sp**2 * np.exp(-m_sp * s) / s
            dI_dP = I_val * (-m_sp * s - 1.0) / (2.0 * s)
        else:
            I_val = 8.0 * np.pi**2 / P
            dI_dP = -8.0 * np.pi**2 / P**2

        dV3_dd = -6.0 * gamma * C**3 * dI_dP * 3.0  # dP/dd = 3 for equilateral
        F_3 = -dV3_dd  # force = -dV/dd (positive = repulsive)

        # Step 3: pairwise Hessian at d0
        # V_2_total = 3 * V_2(d) for equilateral with equal C
        # V_2(d) = -4*pi*C^2/d + 8*pi*beta*C^2/d^2 - 24*pi*beta*C^3/d^3
        # dV_2/dd = 4*pi*C^2/d^2 - 16*pi*beta*C^2/d^3 + 72*pi*beta*C^3/d^4
        # d^2V_2/dd^2 = -8*pi*C^2/d^3 + 48*pi*beta*C^2/d^4 - 288*pi*beta*C^3/d^5
        d2V2 = (-8.0 * np.pi * C**2 / d0**3
                + 48.0 * np.pi * beta * C**2 / d0**4
                - 288.0 * np.pi * beta * C**3 / d0**5)
        H_2 = 3.0 * d2V2  # total Hessian from 3 pairs

        # Step 4: perturbative correction
        if abs(H_2) > 1e-30:
            delta_d = -F_3 / H_2
        else:
            delta_d = 0.0

        d_corrected = d0 + delta_d

        result['d_corrected_' + label] = d_corrected
        result['delta_d_' + label] = delta_d
        result['F_3_' + label] = F_3
        result['H_2_' + label] = H_2
        result['V3_V2_ratio_' + label] = abs(
            6.0 * gamma * C**3 * I_val / (
                3.0 * (-4*np.pi*C**2/d0 + 8*np.pi*beta*C**2/d0**2
                        - 24*np.pi*beta*C**3/d0**3))
        ) if d0 > 0 else None

    return result


def collinear_pairwise_equilibrium_numerical(C1, C2, C3, beta, gamma=None):
    """
    Pairwise-only collinear equilibrium by numerical root scan.

    STATUS: NUMERICAL (sign-change detection + root_scalar).
    Uses pairwise forces only (no 3-body corrections).
    For equal masses, the symmetric ansatz d12=d23=d reduces
    the problem to finding zeros of the total force on mass 1.

    For unequal masses, the equilibrium conditions are two equations:
        dV/d(d12) = 0
        dV/d(d23) = 0

    For equal masses (C1=C2=C3=C), by the mirror symmetry of Euler
    configurations, the symmetric case d12=d23=d is an equilibrium
    candidate, reducing to the same polynomial as equilateral.

    Returns
    -------
    dict with equilibrium separations or None
    """
    if gamma is None:
        gamma = beta

    # For equal masses, symmetric collinear: d12 = d23 = d, d13 = 2d
    if abs(C1 - C2) < 1e-12 and abs(C2 - C3) < 1e-12:
        C = C1
        # V_total = V(d12) + V(d23) + V(d13) = 2*V(d) + V(2d)
        # dV_total/dd = 2*V'(d) + 2*V'(2d) = 0
        # Where V'(d) = 4*pi*C^2/d^2 - 16*pi*beta*C^2/d^3 + 72*pi*beta*C^3/d^4
        # This gives a polynomial in d after substitution

        # Numerical approach: find d where total force on middle mass = 0
        def force_on_end_mass(d):
            """Total pairwise force on mass 1 (at x=0) from masses at x=d and x=2d."""
            if d < 0.01:
                return 1e10
            # Force from mass 1 (at distance d, pulling to the left = negative)
            F1 = -(-4*np.pi*C**2/d**2 + 16*np.pi*beta*C**2/d**3
                    - 72*np.pi*beta*C**3/d**4)
            # Force from mass 3 (at distance d, pulling to the right = positive)
            F3 = -4*np.pi*C**2/d**2 + 16*np.pi*beta*C**2/d**3 \
                 - 72*np.pi*beta*C**3/d**4
            # By symmetry F1 + F3 = 0 for middle mass. Check force on end mass:
            # Force on mass 1 from mass 2 (dist d) and mass 3 (dist 2d)
            F12 = -4*np.pi*C**2/d**2 + 16*np.pi*beta*C**2/d**3 \
                  - 72*np.pi*beta*C**3/d**4
            F13 = -4*np.pi*C**2/(2*d)**2 + 16*np.pi*beta*C**2/(2*d)**3 \
                  - 72*np.pi*beta*C**3/(2*d)**4
            return F12 + F13  # total force on mass 1

        # Scan for sign changes
        d_scan = np.linspace(0.5*C, 10*beta, 2000)
        f_vals = np.array([force_on_end_mass(d) for d in d_scan])

        equilibria = []
        for i in range(len(f_vals) - 1):
            if f_vals[i] * f_vals[i+1] < 0:
                try:
                    result = root_scalar(force_on_end_mass,
                                        bracket=[d_scan[i], d_scan[i+1]])
                    if result.converged:
                        equilibria.append(result.root)
                except Exception:
                    pass

        return {
            'type': 'symmetric_collinear',
            'C': C,
            'equilibria': equilibria,
            'n_equilibria': len(equilibria),
            'message': f"Found {len(equilibria)} symmetric collinear equilibria"
        }

    else:
        # General unequal case: use numerical minimization
        return {
            'type': 'general_collinear',
            'message': "Use numerical_equilibrium() for unequal masses",
            'equilibria': []
        }


def ngon_pairwise_equilibrium_numerical(n, C, beta, gamma=None):
    """
    Find equilibrium radius R for n equal masses on a regular polygon.

    STATUS: NUMERICAL (sign-change detection + root_scalar on radial force).
    Uses pairwise forces only.

    The total energy depends only on R (by symmetry).
    Side length: s = 2*R*sin(pi/n).
    Distance between mass i and mass j: d_ij = 2*R*sin(pi*|i-j|/n).

    Total pairwise energy:
        V = sum_{i<j} V_2(d_ij)

    For large n, the dominant contribution comes from nearest neighbors.
    """
    if gamma is None:
        gamma = beta

    def total_energy(R):
        if R < 0.01:
            return 1e10
        V = 0.0
        for k in range(1, n):
            # Distance from mass 0 to mass k
            d = 2 * R * np.sin(np.pi * k / n)
            # 2-body potential
            V += -4*np.pi*C**2/d + 8*np.pi*beta*C**2/d**2 \
                 - 24*np.pi*beta*C**3/d**3
        return V * n / 2  # count each pair once

    def total_force_radial(R):
        """Radial force on one mass (positive = outward)."""
        if R < 0.01:
            return 1e10
        F_radial = 0.0
        for k in range(1, n):
            theta_k = np.pi * k / n
            d = 2 * R * np.sin(theta_k)
            dd_dR = 2 * np.sin(theta_k)  # d(d_ij)/dR
            # dV_2/dd:
            dV = 4*np.pi*C**2/d**2 - 16*np.pi*beta*C**2/d**3 \
                 + 72*np.pi*beta*C**3/d**4
            # Force = -dV/dR = -dV/dd * dd/dR
            # Radial component: project along radial direction from center
            # For mass 0 at angle 0, mass k at angle 2*pi*k/n,
            # the line connecting them makes angle pi*k/n with the radial
            cos_angle = np.cos(np.pi * k / n)
            F_radial += -dV * dd_dR * cos_angle
        return F_radial

    # Scan for equilibria
    R_scan = np.linspace(0.3*C, 8*beta, 2000)
    f_vals = np.array([total_force_radial(R) for R in R_scan])

    # Find genuine zero crossings, filtering numerical noise.
    # Use central differences of F to confirm the root is a true sign change,
    # not a numerical artifact near F~0.
    f_scale = np.max(np.abs(f_vals)) if np.max(np.abs(f_vals)) > 0 else 1.0
    threshold = f_scale * 1e-5

    equilibria = []
    for i in range(len(f_vals) - 1):
        if f_vals[i] * f_vals[i+1] < 0:
            # Skip tiny sign changes (noise near zero)
            if max(abs(f_vals[i]), abs(f_vals[i+1])) < threshold:
                continue
            try:
                result = root_scalar(total_force_radial,
                                    bracket=[R_scan[i], R_scan[i+1]])
                if result.converged:
                    R_eq = result.root
                    # Verify with a finer check: evaluate F on both sides
                    dR = (R_scan[1] - R_scan[0]) * 0.1
                    F_left = total_force_radial(R_eq - dR)
                    F_right = total_force_radial(R_eq + dR)
                    if F_left * F_right < 0 and max(abs(F_left), abs(F_right)) > threshold:
                        side = 2 * R_eq * np.sin(np.pi / n)
                        equilibria.append({
                            'R': R_eq,
                            'side': side,
                            'energy': total_energy(R_eq),
                        })
            except Exception:
                pass

    return {
        'n': n,
        'C': C,
        'beta': beta,
        'equilibria': equilibria,
        'n_equilibria': len(equilibria),
        'message': f"{n}-gon: found {len(equilibria)} equilibria"
    }


def numerical_equilibrium(positions_init, C_values, beta, gamma=None,
                          energy_func=None):
    """
    Find equilibrium by minimizing total energy numerically.

    Parameters
    ----------
    positions_init : (n, 3) array
        Initial guess for positions.
    C_values : (n,) array
    beta, gamma : float
    energy_func : callable, optional
        total_energy(positions, C_values, beta, gamma) → float.
        If None, uses the pairwise-only approximation.

    Returns
    -------
    dict with optimized positions, energy, convergence info
    """
    if gamma is None:
        gamma = beta

    n = len(C_values)

    if energy_func is None:
        # Default: pairwise only
        def energy_func(pos, C_vals, b, g):
            V = 0.0
            for i in range(len(C_vals)):
                for j in range(i+1, len(C_vals)):
                    d = np.linalg.norm(pos[i] - pos[j])
                    d = max(d, 1e-6)
                    Ci, Cj = C_vals[i], C_vals[j]
                    V += -4*np.pi*Ci*Cj/d \
                         + 8*np.pi*b*Ci*Cj/d**2 \
                         - 24*np.pi*g*Ci*Cj*(Ci+Cj)/(2*d**3)
            return V

    # Remove center-of-mass and rotation degrees of freedom
    # by fixing mass 1 at origin and mass 2 on x-axis
    def energy_reduced(params):
        pos = positions_init.copy()
        # Only vary positions of masses 2..n relative to mass 1
        # (mass 1 stays at origin)
        idx = 0
        for i in range(1, n):
            for k in range(3):
                pos[i, k] = params[idx]
                idx += 1
        return energy_func(pos, C_values, beta, gamma)

    # Initial parameters (flatten positions 1..n)
    x0 = positions_init[1:].flatten()

    result = minimize(energy_reduced, x0, method='Nelder-Mead',
                     options={'maxiter': 50000, 'xatol': 1e-10, 'fatol': 1e-12})

    # Reconstruct positions
    pos_final = positions_init.copy()
    idx = 0
    for i in range(1, n):
        for k in range(3):
            pos_final[i, k] = result.x[idx]
            idx += 1

    return {
        'positions': pos_final,
        'energy': result.fun,
        'converged': result.success,
        'message': result.message,
        'n_iterations': result.nit,
    }
