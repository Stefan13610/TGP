"""
configurations.py  --  TGP n-body: preset configurations
=========================================================
Standard n-body configurations for testing.

Each configuration returns:
  positions: (n, 3) array of 3D positions
  C_values:  (n,) array of dimensionless source strengths
  name:      descriptive string

All in dimensionless units (x = r/r_0).
"""

import numpy as np


def equilateral_triangle(d, C=0.3):
    """
    Three equal masses at vertices of equilateral triangle in xy-plane.

    Parameters
    ----------
    d : float
        Side length.
    C : float
        Dimensionless source strength (same for all three).

    Returns
    -------
    positions : (3, 3) ndarray
    C_values  : (3,) ndarray
    """
    h = d * np.sqrt(3) / 2
    positions = np.array([
        [0.0,       2*h/3,  0.0],   # top vertex
        [-d/2.0,   -h/3,    0.0],   # bottom-left
        [ d/2.0,   -h/3,    0.0],   # bottom-right
    ])
    return positions, np.full(3, C), f"equilateral d={d:.2f} C={C}"


def collinear_equal(d12, d23, C=0.3):
    """
    Three equal masses on the x-axis (Euler-type configuration).

    Parameters
    ----------
    d12, d23 : float
        Separations between masses 1-2 and 2-3.
    C : float
        Source strength (equal).

    Returns
    -------
    positions, C_values
    """
    positions = np.array([
        [0.0,        0.0, 0.0],
        [d12,        0.0, 0.0],
        [d12 + d23,  0.0, 0.0],
    ])
    return positions, np.full(3, C), f"collinear d12={d12:.2f} d23={d23:.2f}"


def collinear_unequal(d12, d23, C1, C2, C3):
    """Three unequal masses on the x-axis."""
    positions = np.array([
        [0.0,        0.0, 0.0],
        [d12,        0.0, 0.0],
        [d12 + d23,  0.0, 0.0],
    ])
    return positions, np.array([C1, C2, C3]), f"collinear unequal"


def regular_ngon(n, R, C=0.3):
    """
    n equal masses at vertices of a regular polygon in the xy-plane.

    Parameters
    ----------
    n : int
        Number of bodies.
    R : float
        Circumradius (distance from center to each vertex).
    C : float
        Source strength.
    """
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    positions = np.column_stack([
        R * np.cos(angles),
        R * np.sin(angles),
        np.zeros(n),
    ])
    side = 2 * R * np.sin(np.pi / n)
    return positions, np.full(n, C), f"{n}-gon R={R:.2f} side={side:.2f}"


def two_body(d, C1=0.3, C2=0.3):
    """Two bodies separated by distance d along x-axis."""
    positions = np.array([
        [-d/2, 0.0, 0.0],
        [ d/2, 0.0, 0.0],
    ])
    return positions, np.array([C1, C2]), f"2-body d={d:.2f}"


def restricted_3body(d_heavy, C_heavy, C_test, test_pos):
    """
    Restricted 3-body problem: two heavy masses + one test mass.

    Parameters
    ----------
    d_heavy : float
        Separation between the two heavy masses (along x-axis).
    C_heavy : float
        Source strength of each heavy mass.
    C_test : float
        Source strength of test mass (C_test << C_heavy).
    test_pos : (3,) array
        Position of test mass.

    Returns
    -------
    positions, C_values
    """
    positions = np.array([
        [-d_heavy/2, 0.0, 0.0],
        [ d_heavy/2, 0.0, 0.0],
        test_pos,
    ])
    C_values = np.array([C_heavy, C_heavy, C_test])
    return positions, C_values, f"restricted 3-body d={d_heavy:.2f}"


def figure_eight_initial(C=0.3, scale=1.0):
    """
    Figure-eight orbit initial conditions (Chenciner-Montgomery).

    In Newtonian gravity, three equal masses can follow a figure-eight
    trajectory. We set up the initial conditions and test whether
    TGP nonlinearity preserves or destroys this orbit.

    Parameters
    ----------
    C : float
        Source strength.
    scale : float
        Overall length scale.

    Returns
    -------
    positions : (3, 3) ndarray
    velocities : (3, 3) ndarray
    C_values : (3,) ndarray
    """
    # Chenciner-Montgomery initial conditions (normalized)
    # x1 = -x2, x3 = 0; v1 = v2, v3 = -2*v1
    x1 = np.array([0.97000436, -0.24308753, 0.0]) * scale
    v1 = np.array([0.93240737/2, 0.86473146/2, 0.0])

    positions = np.array([x1, -x1, [0.0, 0.0, 0.0]])
    velocities = np.array([v1, v1, -2*v1])

    return positions, velocities, np.full(3, C), "figure-8"


def nested_binary(d_inner, d_outer, C_inner=0.3, C_outer=0.3):
    """
    Hierarchical 4-body: inner binary + outer binary (or single).

    Inner binary: two masses at ±d_inner/2 on x-axis.
    Outer mass: at (d_outer, 0, 0).

    Parameters
    ----------
    d_inner : float
        Inner binary separation.
    d_outer : float
        Distance from inner binary center to outer mass.
    """
    positions = np.array([
        [-d_inner/2, 0.0, 0.0],
        [ d_inner/2, 0.0, 0.0],
        [d_outer,    0.0, 0.0],
    ])
    C_values = np.array([C_inner, C_inner, C_outer])
    return positions, C_values, f"nested binary d_in={d_inner:.2f} d_out={d_outer:.2f}"
