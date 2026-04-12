"""
bridge_nbody.py -- Canonical bridge for the N-body layer
========================================================

This module provides the working bridge used inside `nbody`:

    classical defect g(r)
        -> EFT projection
        -> effective Yukawa source (C_eff, m_sp)
        -> N-body inputs for V_2, V_3 and EOM

The goal is not to claim a full substrate-to-metric derivation.  Instead, the
module records the minimal bridge that is already operational in the package
and reusable by higher layers.

Canonical chain
---------------
1. Classical defect:
       g^2 g'' + g (g')^2 + (2/r) g^2 g' = g^2 (1 - g)
   with oscillatory tail ~ sin(r)/r.

2. Effective screening mass:
       m_sp^2 = 3*gamma - 2*beta

3. EFT projection:
       C_eff = integral_0^inf delta(r) exp(-m_sp r) r dr,
       delta(r) = 1 - g(r)

4. Effective N-body source:
       delta_eff(r) = C_eff exp(-m_sp r) / r

5. N-body dynamics:
       C_i x_i'' = -grad_i [V_2 + V_3]

This module wraps the existing defect and projection utilities from
`yukawa_from_defect.py` behind a smaller, canonical API.
"""

from __future__ import annotations

from typing import Any, Dict

from .tgp_field import screening_mass
from .yukawa_from_defect import compute_C_eff_projection, solve_defect_standard


def solve_classical_defect(
    g0: float,
    beta: float = 1.0,
    gamma: float = 1.0,
    *,
    kinetic: str = "full",
    r_max: float = 60.0,
    n_eval: int = 5000,
    rtol: float = 1e-10,
    atol: float = 1e-13,
) -> Dict[str, Any]:
    """Solve the classical defect used as input to the N-body bridge.

    This uses the standard TGP defect equation (no classical stabilization),
    so the resulting tail is oscillatory.  The output can then be projected
    onto the effective Yukawa source sector.
    """
    return solve_defect_standard(
        g0,
        beta=beta,
        gamma=gamma,
        kinetic=kinetic,
        r_max=r_max,
        n_eval=n_eval,
        rtol=rtol,
        atol=atol,
    )


def derive_effective_source_from_defect(
    g0: float,
    beta: float = 1.0,
    gamma: float = 1.0,
    *,
    kinetic: str = "full",
    r_max: float = 60.0,
    n_eval: int = 5000,
    rtol: float = 1e-10,
    atol: float = 1e-13,
) -> Dict[str, Any]:
    """Return the canonical EFT source parameters derived from a defect.

    The returned dictionary is intentionally small and stable, so other layers
    can treat it as the supported hand-off from the soliton sector into the
    N-body sector.
    """
    defect = solve_classical_defect(
        g0,
        beta=beta,
        gamma=gamma,
        kinetic=kinetic,
        r_max=r_max,
        n_eval=n_eval,
        rtol=rtol,
        atol=atol,
    )

    m_sp = screening_mass(beta, gamma)
    C_eff = compute_C_eff_projection(defect["r"], defect["g"], m_sp)
    return {
        "bridge_status": "EFT-DERIVED",
        "g0": float(g0),
        "beta": float(beta),
        "gamma": float(gamma),
        "kinetic_mode": kinetic,
        "m_sp": float(m_sp),
        "C_eff": float(C_eff),
        "tail_type": "oscillatory",
        "osc_amplitude": float(defect["osc_amplitude"]),
        "osc_phase": float(defect["osc_phase"]),
        "osc_residual": float(defect["osc_residual"]),
        "r": defect["r"],
        "g": defect["g"],
        "gp": defect["gp"],
    }


def build_nbody_source_from_defect(
    g0: float,
    beta: float = 1.0,
    gamma: float = 1.0,
    *,
    kinetic: str = "full",
    r_max: float = 60.0,
    n_eval: int = 5000,
    rtol: float = 1e-10,
    atol: float = 1e-13,
) -> Dict[str, Any]:
    """Return the minimal source record needed by the effective N-body layer."""
    derived = derive_effective_source_from_defect(
        g0,
        beta=beta,
        gamma=gamma,
        kinetic=kinetic,
        r_max=r_max,
        n_eval=n_eval,
        rtol=rtol,
        atol=atol,
    )
    return {
        "C": derived["C_eff"],
        "m_sp": derived["m_sp"],
        "beta": derived["beta"],
        "gamma": derived["gamma"],
        "bridge_status": derived["bridge_status"],
        "source_model": "effective_yukawa",
    }
