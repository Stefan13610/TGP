"""
dynamics_backends.py — wybór backendu sił/potencjału pod integratory
=====================================================================

Zwraca parę ``(acc_fn, pot_fn)`` zgodną z ``dynamics_v2.leapfrog_integrate`` /
``rk45_integrate``: obie mają sygnaturę ``(positions, C_values)``.

Warstwa analityczna (wzory) jest w ``eom_tgp.py`` i notatce
``tgp_nbody_lagrangian_eom.tex``; ten moduł to wyłącznie **wybór implementacji**
liczenia $I$ i $V$.
"""

from __future__ import annotations

from typing import Callable, Literal, Tuple

import numpy as np

try:
    from . import dynamics_v2
    from . import eom_tgp
    from . import three_body_force_exact
except ImportError:  # `python ex138_*.py` z sys.path = nbody/
    import dynamics_v2
    import eom_tgp
    import three_body_force_exact

BackendName = Literal[
    "pairwise",
    "coulomb_3b",
    "yukawa_feynman",
    "yukawa_saddle_approx",
]

#: Dostępne identyfikatory (np. CLI, skrypty diagnostyczne).
TGP_INTEGRATION_BACKENDS: tuple[BackendName, ...] = (
    "pairwise",
    "coulomb_3b",
    "yukawa_feynman",
    "yukawa_saddle_approx",
)

ArrayFn = Callable[[np.ndarray, np.ndarray], np.ndarray]
ScalarFn = Callable[[np.ndarray, np.ndarray], float]


def total_force_from_accelerations(
    accelerations: np.ndarray, C_values: np.ndarray
) -> np.ndarray:
    """
    Łączna siła wewnętrzna: ``sum_i F_i = sum_i C_i * a_i`` (przy m_i = C_i).

    Dla potencjału zależnego tylko od względnych odległości oczekujemy
    ``~0`` z dokładnością numeryczną (niezmiennik translacyjny).
    """
    C_values = np.asarray(C_values, dtype=float)
    accelerations = np.asarray(accelerations, dtype=float)
    return np.sum(C_values[:, None] * accelerations, axis=0)


def total_torque_about_origin_from_accelerations(
    positions: np.ndarray,
    accelerations: np.ndarray,
    C_values: np.ndarray,
) -> np.ndarray:
    """
    Moment sił względem początku układu: ``sum_i x_i x (C_i a_i) = sum_i x_i x F_i``.

    Gdy $V$ zależy wyłącznie od odległości $d_{ij}=|\\mathbf{x}_i-\\mathbf{x}_j|$,
    potencjał jest invarianty względem globalnego $SO(3)$, stąd
    $\\mathrm{d}\\mathbf{L}/\\mathrm{d}t=\\mathbf{0}$ dla izolowanego układu
    ($\\mathbf{L}=\\sum_i C_i\\,\\mathbf{x}_i\\times\\dot{\\mathbf{x}}_i$),
    więc $\\sum_i \\mathbf{x}_i\\times\\mathbf{F}_i\\approx\\mathbf{0}$ numerycznie.
    """
    pos = np.asarray(positions, dtype=float)
    acc = np.asarray(accelerations, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    F = C_values[:, None] * acc
    return np.sum(np.cross(pos, F, axis=1), axis=0)


def total_linear_momentum(velocities: np.ndarray, C_values: np.ndarray) -> np.ndarray:
    """``P = sum_i C_i v_i`` (przy m_i = C_i)."""
    v = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    return np.sum(C_values[:, None] * v, axis=0)


def center_of_mass(positions: np.ndarray, C_values: np.ndarray) -> np.ndarray:
    """``R_cm = (sum_i C_i x_i) / M_tot``."""
    pos = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    mtot = float(np.sum(C_values))
    return np.sum(C_values[:, None] * pos, axis=0) / mtot


def center_of_mass_velocity(velocities: np.ndarray, C_values: np.ndarray) -> np.ndarray:
    """``V_cm = P / M_tot``."""
    v = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    mtot = float(np.sum(C_values))
    return total_linear_momentum(v, C_values) / mtot


def conjugate_momenta(velocities: np.ndarray, C_values: np.ndarray) -> np.ndarray:
    """``p_i = C_i v_i`` (jak w Hamiltonianie z ``tgp_nbody_lagrangian_eom.tex``)."""
    v = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    return C_values[:, None] * v


def center_of_mass_acceleration(
    accelerations: np.ndarray, C_values: np.ndarray
) -> np.ndarray:
    """
    ``a_cm = (sum_i C_i a_i) / M_tot = (sum_i F_i) / M_tot``.

    Przy czysto wewnętrznym $V$ od $d_{ij}$: $a_{\\rm cm}\\approx 0$.
    """
    acc = np.asarray(accelerations, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    mtot = float(np.sum(C_values))
    return total_force_from_accelerations(acc, C_values) / mtot


def total_angular_momentum(
    positions: np.ndarray, velocities: np.ndarray, C_values: np.ndarray
) -> np.ndarray:
    """``L = sum_i C_i x_i x v_i`` względem początku układu współrzędnych."""
    pos = np.asarray(positions, dtype=float)
    v = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    return np.sum(C_values[:, None] * np.cross(pos, v, axis=1), axis=0)


def forces_from_potential_central_diff(
    pot_fn: ScalarFn,
    positions: np.ndarray,
    C_values: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    """
    Przybliżenie ``F_i = -nabla_i V`` przez centralne różnice składowych.

    Koszt: ``O(N * 3 * 2)`` ewaluacji ``pot_fn`` — tylko do walidacji, nie do dynamiki.
    """
    pos = np.asarray(positions, dtype=float).copy()
    C_values = np.asarray(C_values, dtype=float)
    n, dim = pos.shape
    assert dim == 3
    F = np.zeros_like(pos)
    for i in range(n):
        for k in range(3):
            pos[i, k] += eps
            vp = float(pot_fn(pos, C_values))
            pos[i, k] -= 2.0 * eps
            vm = float(pot_fn(pos, C_values))
            pos[i, k] += eps
            F[i, k] = -(vp - vm) / (2.0 * eps)
    return F


def build_tgp_integration_pair(
    backend: BackendName = "coulomb_3b",
    *,
    beta: float = 1.0,
    gamma: float | None = None,
    softening: float = 1e-6,
    include_3body: bool = True,
    n_quad_feynman: int = 40,
) -> Tuple[ArrayFn, ScalarFn]:
    """
    Buduje ``(acc_fn, pot_fn)`` dla danego backendu.

    Parameters
    ----------
    backend
        - ``pairwise`` — tylko sektor 2-ciałowy (dokładny).
        - ``coulomb_3b`` — $V_2$ + $I=8\\pi^2/P$ (zgodne z ``eom_tgp`` i
          ``full_potential_tgp(..., use_yukawa=False)``).
        - ``yukawa_feynman`` — $V_2$ (z ``softening``) + dokładne $I_Y$
          (kwadratura 2D; **odległości tripletów bez softeningu**, jak w
          ``three_body_force_exact``).  Do Benettina / stycznej: jawny
          ``∂a/∂x`` — ``lyapunov.acceleration_jacobian_yukawa_feynman_analytic``.
        - ``yukawa_saddle_approx`` — przestarzałe przybliżenie saddle-point z
          ``dynamics_v2.three_body_forces_approximate``; **duży błąd** przy
          $m_{\\rm sp} d \\ll 5$ — tylko do regresji / porównań.

    n_quad_feynman
        Rozdzielczość kwadratury dla ``yukawa_feynman`` (siły i energia $V_3$).
    """
    if gamma is None:
        gamma = beta

    def acc_pairwise(pos: np.ndarray, C: np.ndarray) -> np.ndarray:
        F = dynamics_v2.analytical_forces_tgp_pairwise(
            pos, C, beta, gamma, softening
        )
        return F / C[:, None]

    def pot_pairwise(pos: np.ndarray, C: np.ndarray) -> float:
        return float(dynamics_v2.potential_tgp(pos, C, beta, gamma, softening))

    if backend == "pairwise":
        return acc_pairwise, pot_pairwise

    if backend == "coulomb_3b":

        def acc(pos: np.ndarray, C: np.ndarray) -> np.ndarray:
            return eom_tgp.accelerations_tgp_nbody(
                pos,
                C,
                beta,
                gamma,
                softening,
                include_3body=include_3body,
                three_body_dI=None,
            )

        def pot(pos: np.ndarray, C: np.ndarray) -> float:
            if not include_3body or len(C) < 3:
                return pot_pairwise(pos, C)
            return float(
                dynamics_v2.full_potential_tgp(
                    pos,
                    C,
                    beta,
                    gamma,
                    softening,
                    include_3body=True,
                    use_yukawa=False,
                )
            )

        return acc, pot

    if backend == "yukawa_feynman":

        def acc(pos: np.ndarray, C: np.ndarray) -> np.ndarray:
            F2 = dynamics_v2.analytical_forces_tgp_pairwise(
                pos, C, beta, gamma, softening
            )
            if include_3body and len(C) >= 3:
                F3 = three_body_force_exact.three_body_forces_exact(
                    pos, C, beta=beta, gamma=gamma, n_quad=n_quad_feynman
                )
            else:
                F3 = np.zeros_like(F2)
            return (F2 + F3) / C[:, None]

        def pot(pos: np.ndarray, C: np.ndarray) -> float:
            V2 = dynamics_v2.potential_tgp(pos, C, beta, gamma, softening)
            if include_3body and len(C) >= 3:
                V3 = three_body_force_exact.total_three_body_energy_exact(
                    pos, C, beta=beta, gamma=gamma, n_quad=n_quad_feynman
                )
            else:
                V3 = 0.0
            return float(V2 + V3)

        return acc, pot

    if backend == "yukawa_saddle_approx":

        def acc(pos: np.ndarray, C: np.ndarray) -> np.ndarray:
            F = dynamics_v2.full_forces_tgp(
                pos,
                C,
                beta,
                gamma,
                softening,
                include_3body=include_3body,
                use_yukawa=True,
            )
            return F / C[:, None]

        def pot(pos: np.ndarray, C: np.ndarray) -> float:
            return float(
                dynamics_v2.full_potential_tgp(
                    pos,
                    C,
                    beta,
                    gamma,
                    softening,
                    include_3body=include_3body,
                    use_yukawa=True,
                )
            )

        return acc, pot

    raise ValueError(f"unknown backend: {backend!r}")


def _self_test() -> None:
    rng = np.random.default_rng(1)
    n = 4
    pos = rng.normal(size=(n, 3))
    C = np.abs(rng.normal(size=n)) + 0.15
    beta = gamma = 1.0
    eps = 1e-6
    acc_a, _ = build_tgp_integration_pair(
        "coulomb_3b", beta=beta, gamma=gamma, softening=eps, include_3body=True
    )
    acc_b, _ = build_tgp_integration_pair(
        "pairwise", beta=beta, gamma=gamma, softening=eps, include_3body=True
    )
    a_full = acc_a(pos, C)
    a_pair = acc_b(pos, C)
    assert np.all(np.isfinite(a_full))
    assert np.linalg.norm(a_full - a_pair) > 1e-10  # 3B zmienia dynamikę
    F_ref = dynamics_v2.full_forces_tgp(
        pos, C, beta, gamma, eps, include_3body=True, use_yukawa=False
    )
    assert np.max(np.abs(a_full - F_ref / C[:, None])) < 1e-9


def _self_test_com_acceleration() -> None:
    rng = np.random.default_rng(4)
    for backend in ("coulomb_3b", "yukawa_feynman"):
        for n in (3, 5):
            pos = rng.normal(size=(n, 3))
            C = np.abs(rng.normal(size=n)) + 0.17
            acc_fn, _ = build_tgp_integration_pair(
                backend,
                beta=1.0,
                gamma=1.0,
                softening=1e-6,
                include_3body=True,
                n_quad_feynman=18,
            )
            a = acc_fn(pos, C)
            acm = center_of_mass_acceleration(a, C)
            assert np.linalg.norm(acm) < 1e-9


def _self_test_net_force() -> None:
    rng = np.random.default_rng(2)
    _backends: tuple[BackendName, ...] = (
        "pairwise",
        "coulomb_3b",
        "yukawa_feynman",
    )
    for backend in _backends:
        for n in (2, 3, 4):
            pos = rng.normal(scale=1.2, size=(n, 3))
            C = np.abs(rng.normal(size=n)) + 0.18
            acc_fn, _ = build_tgp_integration_pair(
                backend,
                beta=1.0,
                gamma=1.0,
                softening=1e-6,
                include_3body=True,
                n_quad_feynman=18,
            )
            a = acc_fn(pos, C)
            net = total_force_from_accelerations(a, C)
            fnorm = float(np.max(np.linalg.norm(C[:, None] * a, axis=1)))
            assert np.linalg.norm(net) < 1e-9 * max(1.0, fnorm) + 1e-12


def _self_test_net_torque() -> None:
    rng = np.random.default_rng(3)
    _backends: tuple[BackendName, ...] = (
        "pairwise",
        "coulomb_3b",
        "yukawa_feynman",
        "yukawa_saddle_approx",
    )
    for backend in _backends:
        for n in (2, 3, 4):
            pos = rng.normal(scale=1.0, size=(n, 3))
            C = np.abs(rng.normal(size=n)) + 0.18
            acc_fn, _ = build_tgp_integration_pair(
                backend,
                beta=1.0,
                gamma=1.0,
                softening=1e-6,
                include_3body=True,
                n_quad_feynman=18,
            )
            a = acc_fn(pos, C)
            tau = total_torque_about_origin_from_accelerations(pos, a, C)
            F = C[:, None] * a
            fmag = float(np.max(np.linalg.norm(F, axis=1)))
            rmax = float(np.max(np.linalg.norm(pos, axis=1)))
            scale = rmax * fmag
            assert np.linalg.norm(tau) < 1e-8 * max(1.0, scale) + 1e-12


if __name__ == "__main__":
    _self_test()
    _self_test_net_force()
    _self_test_net_torque()
    _self_test_com_acceleration()
    print("dynamics_backends: self-test OK (force + torque + a_cm ~0)")
