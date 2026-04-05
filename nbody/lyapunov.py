# -*- coding: utf-8 -*-
"""
lyapunov.py
===========
Największy wykładnik Lapunowa (Benettin) dla układu drugiego rzędu

    ddot x_i = a_i(x),   m_i = C_i,   a_i = F_i / C_i.

Stosuje sprzężone RK4 na rozszerzonym stanie Z = [x, v, delta_x, delta_v]
z liniaryzacją

    d/dt delta_x = delta_v,
    d/dt delta_v = (partial a / partial x) delta_x,

gdzie Jacobian ``a`` względem pozycji jest liczony różnicami centralnymi,
lub podawany jawnie (``acceleration_jacobian_newton_softened`` dla Newtona,
``acceleration_jacobian_tgp_pairwise_softened`` dla ``V_2``,
``acceleration_jacobian_yukawa_feynman_split`` — ``V_2`` analitycznie + FD na ``V_3``;
``acceleration_jacobian_yukawa_feynman_analytic`` — ``V_2`` + ``V_3`` bez FD).

Literatura: standardowy algorytm Benettina z okresową renormą ||delta|| -> 1.
Dla ``k>1``: propagacja macierzy stycznej ``U in R^{6n x k}`` i ortogonalizacja
QR co ``renorm_every`` kroków — pierwsze ``k`` wykładników (Benettin–Georgeli).

Dodatkowo: wariant **leapfrog (Velocity Verlet)** na ``(x,v)`` ze sprzężoną
liniaryzacją dyskretną (dwa ewaluacje ``J`` na krok) —
``largest_lyapunov_exponent_benettin_leapfrog``,
``lyapunov_spectrum_benettin_leapfrog`` — lepiej pasuje do symplektycznej
bazy niż RK4 na rozszerzonym ODE.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np

ArrayAccFn = Callable[[np.ndarray, np.ndarray], np.ndarray]
PositionJacobianFn = Callable[[np.ndarray, np.ndarray], np.ndarray]


def acceleration_jacobian_positions(
    positions: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    eps: float = 1e-5,
) -> np.ndarray:
    """
    ``A[i*3+alpha, j*3+beta] = d a_{i,alpha} / d x_{j,beta}`` (układ wierszowy).

    Parametry
    ---------
    positions : (n, 3)
    C_values : (n,)
    acc_fn : (pos, C) -> (n, 3) przyspieszenia
    eps : krok różnic centralnych
    """
    pos = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    n = pos.shape[0]
    m = 3 * n
    A = np.zeros((m, m))
    for j in range(n):
        for beta in range(3):
            dp = np.zeros_like(pos)
            dp[j, beta] = eps
            col = 3 * j + beta
            a_p = acc_fn(pos + dp, C_values).reshape(-1)
            a_m = acc_fn(pos - dp, C_values).reshape(-1)
            A[:, col] = (a_p - a_m) / (2.0 * eps)
    return A


def acceleration_jacobian_newton_softened(
    positions: np.ndarray,
    C_values: np.ndarray,
    *,
    G: float = 1.0,
    softening: float = 1e-6,
) -> np.ndarray:
    """
    Jawny Jacobian ``d a / d x`` dla przyspieszenia Newtona ze ``softening`` jak w
    ``dynamics_v2.forces_newton`` (masa ``C_i``, ``a_i = F_i / C_i``).

    Blok ``[3i:3i+3, 3j:3j+3]`` to ``d a_i / d x_j``.
    """
    pos = np.asarray(positions, dtype=float)
    C = np.asarray(C_values, dtype=float)
    n = pos.shape[0]
    m = 3 * n
    s2 = softening**2
    dFdq = np.zeros((m, m))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            u = pos[j] - pos[i]
            r2 = float(np.dot(u, u) + s2)
            inv_r2 = 1.0 / r2
            inv_r = np.sqrt(inv_r2)
            inv_r3 = inv_r2 * inv_r
            inv_r5 = inv_r2 * inv_r3
            cj = G * C[j]
            outer = np.outer(u, u)
            dTd_xj = cj * (inv_r3 * np.eye(3) - 3.0 * inv_r5 * outer)
            dTd_xi = cj * (-inv_r3 * np.eye(3) + 3.0 * inv_r5 * outer)
            dFdq[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] += dTd_xj
            dFdq[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] += dTd_xi
    J = np.zeros((m, m))
    for i in range(n):
        inv_ci = 1.0 / C[i]
        J[3 * i : 3 * i + 3, :] = inv_ci * dFdq[3 * i : 3 * i + 3, :]
    return J


def acceleration_jacobian_tgp_pairwise_softened(
    positions: np.ndarray,
    C_values: np.ndarray,
    *,
    beta: float,
    gamma: float | None = None,
    softening: float = 1e-6,
) -> np.ndarray:
    """
    Jawny ``∂a/∂x`` dla ``a_i = F_i/C_i`` przy ``F`` z ``dynamics_v2.analytical_forces_tgp_pairwise``.
    """
    from .dynamics_v2 import pairwise_tgp_force_jacobian

    dFdq = pairwise_tgp_force_jacobian(
        positions, C_values, beta, gamma=gamma, softening=softening
    )
    C = np.asarray(C_values, dtype=float)
    n = len(C)
    m = 3 * n
    J = np.zeros((m, m))
    for i in range(n):
        inv_ci = 1.0 / C[i]
        J[3 * i : 3 * i + 3, :] = inv_ci * dFdq[3 * i : 3 * i + 3, :]
    return J


def acceleration_jacobian_yukawa_feynman_split(
    positions: np.ndarray,
    C_values: np.ndarray,
    *,
    beta: float,
    gamma: float | None = None,
    softening: float = 1e-6,
    n_quad_feynman: int = 40,
    jac_eps: float = 1e-5,
) -> np.ndarray:
    """
    ``∂a/∂x`` dla backendu ``yukawa_feynman``: sektor parowy dokładnie
    (``acceleration_jacobian_tgp_pairwise_softened``), irreducible ``V_3``
    przez różnice centralne na ``a_3 = F_3/C`` (tylko siły Feynmana).

    Zgodne z ``dynamics_backends.build_tgp_integration_pair(\"yukawa_feynman\")``
    przy ``include_3body=True``. Dla ``N < 3`` redukuje się do samego ``V_2``.
    """
    from . import three_body_force_exact

    J2 = acceleration_jacobian_tgp_pairwise_softened(
        positions,
        C_values,
        beta=beta,
        gamma=gamma,
        softening=softening,
    )
    n = len(np.asarray(C_values))
    if n < 3:
        return J2

    g = gamma if gamma is not None else beta
    nq = int(n_quad_feynman)

    def acc3_only(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        F3 = three_body_force_exact.three_body_forces_exact(
            p, c, beta=beta, gamma=g, n_quad=nq
        )
        return F3 / c[:, None]

    J3 = acceleration_jacobian_positions(
        positions, C_values, acc3_only, eps=jac_eps
    )
    return J2 + J3


def acceleration_jacobian_three_body_yukawa_exact(
    positions: np.ndarray,
    C_values: np.ndarray,
    *,
    beta: float,
    gamma: float | None = None,
    n_quad_feynman: int = 40,
) -> np.ndarray:
    """
    ``∂a/∂x`` tylko z irreducible Yukawa ``V_3`` (``F_3/C``), bez ``V_2``.
    """
    from . import three_body_force_exact

    dFdq = three_body_force_exact.three_body_force_jacobian_exact(
        positions,
        C_values,
        beta=beta,
        gamma=gamma,
        n_quad=int(n_quad_feynman),
    )
    C = np.asarray(C_values, dtype=float)
    n = len(C)
    m = 3 * n
    J = np.zeros((m, m))
    for i in range(n):
        J[3 * i : 3 * i + 3, :] = dFdq[3 * i : 3 * i + 3, :] / C[i]
    return J


def acceleration_jacobian_yukawa_feynman_analytic(
    positions: np.ndarray,
    C_values: np.ndarray,
    *,
    beta: float,
    gamma: float | None = None,
    softening: float = 1e-6,
    n_quad_feynman: int = 40,
) -> np.ndarray:
    """
    Pełny ``∂a/∂x`` dla ``yukawa_feynman`` (``V_2`` jawny + ``V_3`` przez
    ``three_body_force_jacobian_exact``), bez różnic skończonych.
    """
    J2 = acceleration_jacobian_tgp_pairwise_softened(
        positions,
        C_values,
        beta=beta,
        gamma=gamma,
        softening=softening,
    )
    n = len(np.asarray(C_values))
    if n < 3:
        return J2
    J3 = acceleration_jacobian_three_body_yukawa_exact(
        positions,
        C_values,
        beta=beta,
        gamma=gamma,
        n_quad_feynman=n_quad_feynman,
    )
    return J2 + J3


def _position_jacobian_matrix(
    pos: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
) -> np.ndarray:
    if position_jacobian_fn is not None:
        return position_jacobian_fn(pos, C_values)
    return acceleration_jacobian_positions(pos, C_values, acc_fn, eps=jac_eps)


def _flow_rhs_y(
    y: np.ndarray, C_values: np.ndarray, acc_fn: ArrayAccFn
) -> np.ndarray:
    n = len(C_values)
    pos = y[: 3 * n].reshape(n, 3)
    vel = y[3 * n :].reshape(n, 3)
    acc = acc_fn(pos, C_values)
    return np.concatenate([vel.ravel(), acc.ravel()])


def _extended_rhs(
    Z: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
) -> np.ndarray:
    n = len(C_values)
    y = Z[: 6 * n]
    u = Z[6 * n :]
    pos = y[: 3 * n].reshape(n, 3)
    ky = _flow_rhs_y(y, C_values, acc_fn)
    A = _position_jacobian_matrix(
        pos, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    ux = u[: 3 * n]
    uv = u[3 * n :]
    ku = np.concatenate([uv, A @ ux])
    return np.concatenate([ky, ku])


def _extended_rhs_matrix(
    Z: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
    *,
    n: int,
    k: int,
) -> np.ndarray:
    y = Z[: 6 * n]
    U = Z[6 * n :].reshape(6 * n, k, order="F")
    pos = y[: 3 * n].reshape(n, 3)
    ky = _flow_rhs_y(y, C_values, acc_fn)
    A = _position_jacobian_matrix(
        pos, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    Ux = U[: 3 * n, :]
    Uv = U[3 * n :, :]
    kU = np.vstack([Uv, A @ Ux])
    return np.concatenate([ky, kU.ravel(order="F")])


def _rk4_step_ext(
    Z: np.ndarray,
    dt: float,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
) -> np.ndarray:
    k1 = _extended_rhs(Z, C_values, acc_fn, jac_eps, position_jacobian_fn)
    k2 = _extended_rhs(
        Z + 0.5 * dt * k1, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    k3 = _extended_rhs(
        Z + 0.5 * dt * k2, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    k4 = _extended_rhs(
        Z + dt * k3, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    return Z + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def _rk4_step_ext_matrix(
    Z: np.ndarray,
    dt: float,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
    *,
    n: int,
    k: int,
) -> np.ndarray:
    k1 = _extended_rhs_matrix(
        Z, C_values, acc_fn, jac_eps, position_jacobian_fn, n=n, k=k
    )
    k2 = _extended_rhs_matrix(
        Z + 0.5 * dt * k1,
        C_values,
        acc_fn,
        jac_eps,
        position_jacobian_fn,
        n=n,
        k=k,
    )
    k3 = _extended_rhs_matrix(
        Z + 0.5 * dt * k2,
        C_values,
        acc_fn,
        jac_eps,
        position_jacobian_fn,
        n=n,
        k=k,
    )
    k4 = _extended_rhs_matrix(
        Z + dt * k3,
        C_values,
        acc_fn,
        jac_eps,
        position_jacobian_fn,
        n=n,
        k=k,
    )
    return Z + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def _leapfrog_tangent_step_single(
    x: np.ndarray,
    v: np.ndarray,
    u: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
    dt: float,
    n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Jeden krok Velocity Verlet + sprzężona wariacja (jeden wektor styczny)."""
    pos = x.reshape(n, 3)
    J0 = _position_jacobian_matrix(
        pos, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    a0 = acc_fn(pos, C_values).reshape(-1)
    v_half = v + 0.5 * dt * a0
    ux, uv = u[: 3 * n], u[3 * n :]
    uv_half = uv + 0.5 * dt * (J0 @ ux)
    x_new = x + dt * v_half
    ux_new = ux + dt * uv_half
    pos_new = x_new.reshape(n, 3)
    J1 = _position_jacobian_matrix(
        pos_new, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    a1 = acc_fn(pos_new, C_values).reshape(-1)
    v_new = v_half + 0.5 * dt * a1
    uv_new = uv_half + 0.5 * dt * (J1 @ ux_new)
    u_new = np.concatenate([ux_new, uv_new])
    return x_new, v_new, u_new


def _leapfrog_tangent_step_matrix(
    x: np.ndarray,
    v: np.ndarray,
    U: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    jac_eps: float,
    position_jacobian_fn: Optional[PositionJacobianFn],
    dt: float,
    n: int,
    k: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Velocity Verlet + wariacja na macierzy ``U`` (6n x k, Fortran order w stanie)."""
    pos = x.reshape(n, 3)
    J0 = _position_jacobian_matrix(
        pos, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    a0 = acc_fn(pos, C_values).reshape(-1)
    v_half = v + 0.5 * dt * a0
    Ux = U[: 3 * n, :]
    Uv = U[3 * n :, :]
    Uv_half = Uv + 0.5 * dt * (J0 @ Ux)
    x_new = x + dt * v_half
    Ux_new = Ux + dt * Uv_half
    pos_new = x_new.reshape(n, 3)
    J1 = _position_jacobian_matrix(
        pos_new, C_values, acc_fn, jac_eps, position_jacobian_fn
    )
    a1 = acc_fn(pos_new, C_values).reshape(-1)
    v_new = v_half + 0.5 * dt * a1
    Uv_new = Uv_half + 0.5 * dt * (J1 @ Ux_new)
    U_new = np.vstack([Ux_new, Uv_new])
    return x_new, v_new, U_new


def largest_lyapunov_exponent_benettin_leapfrog(
    positions0: np.ndarray,
    velocities0: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    *,
    t_final: float,
    dt: float,
    renorm_every: int = 25,
    jac_eps: float = 1e-5,
    position_jacobian_fn: Optional[PositionJacobianFn] = None,
    rng: np.random.Generator | None = None,
) -> tuple[float, int]:
    """
    Jak ``largest_lyapunov_exponent_benettin``, lecz orbita i styczna przez
    leapfrog zamiast RK4 na rozszerzonym układzie.
    """
    pos0 = np.asarray(positions0, dtype=float)
    vel0 = np.asarray(velocities0, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    n = len(C_values)
    rng = rng or np.random.default_rng()

    x = pos0.copy().ravel()
    v = vel0.copy().ravel()
    u = rng.standard_normal(6 * n)
    u /= np.linalg.norm(u)

    sum_log = 0.0
    steps = 0
    t = 0.0
    renorm_counter = 0

    while t + dt <= t_final + 1e-12:
        x, v, u = _leapfrog_tangent_step_single(
            x,
            v,
            u,
            C_values,
            acc_fn,
            jac_eps,
            position_jacobian_fn,
            dt,
            n,
        )
        t += dt
        steps += 1
        renorm_counter += 1

        if renorm_counter >= renorm_every:
            nu = float(np.linalg.norm(u))
            if nu < 1e-300:
                return float("nan"), steps
            sum_log += np.log(nu)
            u = u / nu
            renorm_counter = 0

    lam = sum_log / t_final if t_final > 0 else float("nan")
    return lam, steps


def lyapunov_spectrum_benettin_leapfrog(
    positions0: np.ndarray,
    velocities0: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    *,
    n_exponents: int,
    t_final: float,
    dt: float,
    renorm_every: int = 25,
    jac_eps: float = 1e-5,
    position_jacobian_fn: Optional[PositionJacobianFn] = None,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, int]:
    """Pierwsze ``k`` wykładników — Benettin + QR przy integracji leapfrog."""
    pos0 = np.asarray(positions0, dtype=float)
    vel0 = np.asarray(velocities0, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    n = len(C_values)
    k = int(n_exponents)
    dim = 6 * n
    if k < 1 or k > dim:
        raise ValueError(f"n_exponents must be in [1, {dim}], got {k}")
    rng = rng or np.random.default_rng()

    x = pos0.copy().ravel()
    v = vel0.copy().ravel()
    U0 = rng.standard_normal((dim, k))
    U, _ = np.linalg.qr(U0, mode="reduced")

    sum_log = np.zeros(k)
    steps = 0
    t = 0.0
    renorm_counter = 0

    while t + dt <= t_final + 1e-12:
        x, v, U = _leapfrog_tangent_step_matrix(
            x,
            v,
            U,
            C_values,
            acc_fn,
            jac_eps,
            position_jacobian_fn,
            dt,
            n,
            k,
        )
        t += dt
        steps += 1
        renorm_counter += 1

        if renorm_counter >= renorm_every:
            Q, R = np.linalg.qr(U, mode="reduced")
            for j in range(k):
                rjj = float(R[j, j])
                sum_log[j] += np.log(max(abs(rjj), 1e-300))
            U = Q
            renorm_counter = 0

    lambdas = sum_log / t_final if t_final > 0 else np.full(k, np.nan)
    return lambdas, steps


def largest_lyapunov_exponent_benettin(
    positions0: np.ndarray,
    velocities0: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    *,
    t_final: float,
    dt: float,
    renorm_every: int = 25,
    jac_eps: float = 1e-5,
    position_jacobian_fn: Optional[PositionJacobianFn] = None,
    rng: np.random.Generator | None = None,
) -> tuple[float, int]:
    """
    Szacuje największy wykładnik Lapunowa (1/czas) dla podanej trajektorii.

    Jeśli ``position_jacobian_fn`` jest podane: ``(pos, C) -> (3n, 3n)`` zwraca
    ``d a / d x``; w przeciwnym razie Jacobian z różnic centralnych (``jac_eps``).

    Zwraca
    ------
    lambda_max : float
        (1/t_final) * sum log ||delta|| przy renormach (Benettin).
    n_steps : int
        Liczba kroków czasowych.
    """
    pos0 = np.asarray(positions0, dtype=float)
    vel0 = np.asarray(velocities0, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    n = len(C_values)
    rng = rng or np.random.default_rng()

    y0 = np.concatenate([pos0.ravel(), vel0.ravel()])
    u0 = rng.standard_normal(6 * n)
    u0 /= np.linalg.norm(u0)
    Z = np.concatenate([y0, u0])

    sum_log = 0.0
    steps = 0
    t = 0.0
    renorm_counter = 0

    while t + dt <= t_final + 1e-12:
        Z = _rk4_step_ext(
            Z, dt, C_values, acc_fn, jac_eps, position_jacobian_fn
        )
        t += dt
        steps += 1
        renorm_counter += 1

        if renorm_counter >= renorm_every:
            u = Z[6 * n :]
            nu = float(np.linalg.norm(u))
            if nu < 1e-300:
                return float("nan"), steps
            sum_log += np.log(nu)
            Z = Z.copy()
            Z[6 * n :] = u / nu
            renorm_counter = 0

    lam = sum_log / t_final if t_final > 0 else float("nan")
    return lam, steps


def lyapunov_spectrum_benettin(
    positions0: np.ndarray,
    velocities0: np.ndarray,
    C_values: np.ndarray,
    acc_fn: ArrayAccFn,
    *,
    n_exponents: int,
    t_final: float,
    dt: float,
    renorm_every: int = 25,
    jac_eps: float = 1e-5,
    position_jacobian_fn: Optional[PositionJacobianFn] = None,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, int]:
    """
    Pierwsze ``n_exponents`` wykładników Lapunowa (Benettin + QR).

    Opcjonalnie ``position_jacobian_fn(pos, C)`` zamiast FD (patrz
    ``largest_lyapunov_exponent_benettin``).

    Stan styczny: ``U in R^{6n x k}`` (kolumny ortonormalne), układ
    ``dU/dt = J_flow U`` z tym samym ``J_flow`` co w ``largest_lyapunov_exponent_benettin``.
    Po każdym bloku ``renorm_every`` kroków: ``U <- Q`` z rozkładu ``QR(U)``,
    ``sum_log[j] += log|R_{jj}|``.

    Zwraca
    ------
    lambdas : (k,) ndarray
        Szacunki ``lambda_j ~ sum_log[j] / t_final``.
    n_steps : int
    """
    pos0 = np.asarray(positions0, dtype=float)
    vel0 = np.asarray(velocities0, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    n = len(C_values)
    k = int(n_exponents)
    dim = 6 * n
    if k < 1 or k > dim:
        raise ValueError(f"n_exponents must be in [1, {dim}], got {k}")
    rng = rng or np.random.default_rng()

    y0 = np.concatenate([pos0.ravel(), vel0.ravel()])
    U0 = rng.standard_normal((dim, k))
    Q0, _ = np.linalg.qr(U0, mode="reduced")
    Z = np.concatenate([y0, Q0.ravel(order="F")])

    sum_log = np.zeros(k)
    steps = 0
    t = 0.0
    renorm_counter = 0

    while t + dt <= t_final + 1e-12:
        Z = _rk4_step_ext_matrix(
            Z, dt, C_values, acc_fn, jac_eps, position_jacobian_fn, n=n, k=k
        )
        t += dt
        steps += 1
        renorm_counter += 1

        if renorm_counter >= renorm_every:
            y = Z[:dim]
            U = Z[dim:].reshape(dim, k, order="F")
            Q, R = np.linalg.qr(U, mode="reduced")
            for j in range(k):
                rjj = float(R[j, j])
                sum_log[j] += np.log(max(abs(rjj), 1e-300))
            Z = np.concatenate([y, Q.ravel(order="F")])
            renorm_counter = 0

    lambdas = sum_log / t_final if t_final > 0 else np.full(k, np.nan)
    return lambdas, steps


def total_mechanical_energy(
    positions: np.ndarray,
    velocities: np.ndarray,
    C_values: np.ndarray,
    potential_fn: Callable[[np.ndarray, np.ndarray], float],
) -> float:
    """``T + V`` przy ``T = (1/2) sum_i C_i |v_i|^2``."""
    pos = np.asarray(positions, dtype=float)
    vel = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    T = 0.5 * float(np.sum(C_values[:, None] * vel**2))
    return T + float(potential_fn(pos, C_values))


def scale_velocities_match_energy(
    positions: np.ndarray,
    velocities: np.ndarray,
    C_values: np.ndarray,
    potential_fn: Callable[[np.ndarray, np.ndarray], float],
    target_energy: float,
) -> np.ndarray:
    """
    ``v -> s v`` tak, by ``T(s v) + V(x) = target_energy`` (jeśli ``T>0``).

    Jeśli ``||v||`` jest zbyt małe, zwraca ``velocities`` bez zmian — wówczas
    do dopasowania całkowitej energii przy tym samym ``x`` użyj
    ``random_velocities_for_excess_energy``.
    """
    pos = np.asarray(positions, dtype=float)
    vel = np.asarray(velocities, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    T0 = 0.5 * float(np.sum(C_values[:, None] * vel**2))
    if T0 < 1e-30:
        return vel.copy()
    V = float(potential_fn(pos, C_values))
    num = target_energy - V
    if num <= 0:
        return vel.copy()
    s = float(np.sqrt(num / T0))
    return s * vel


def random_velocities_for_excess_energy(
    positions: np.ndarray,
    C_values: np.ndarray,
    potential_fn: Callable[[np.ndarray, np.ndarray], float],
    target_energy: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Losowe ``v`` (izotropowo w przestrzeni 3N) o zadanej energii kinetycznej
    ``T = target_energy - V(x)`` (gdy ta różnica jest dodatnia).

    Przydatne przy ``v=0``: nie da się wtedy dopasować energii skalowaniem ``s v``.
    """
    pos = np.asarray(positions, dtype=float)
    C_values = np.asarray(C_values, dtype=float)
    V = float(potential_fn(pos, C_values))
    T_need = target_energy - V
    if T_need <= 1e-30:
        return np.zeros_like(pos)
    v = rng.standard_normal(pos.shape)
    T0 = 0.5 * float(np.sum(C_values[:, None] * v**2))
    if T0 < 1e-30:
        return np.zeros_like(pos)
    return v * float(np.sqrt(T_need / T0))


def pythagorean_three_body_burrau() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Klasyczne IC problemu pitagorejskiego (2D w z=0), masy 3,4,5, G=1.

    Zwraca ``(positions, velocities, C_values)`` z ``v=0``.
    """
    C = np.array([3.0, 4.0, 5.0], dtype=float)
    pos = np.array(
        [
            [1.0, 3.0, 0.0],
            [-2.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        dtype=float,
    )
    vel = np.zeros_like(pos)
    return pos, vel, C


def _self_test() -> None:
    rng = np.random.default_rng(0)
    n = 2
    C = np.array([1.0, 2.0])
    pos = rng.normal(size=(n, 3))

    def acc_harmonic(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        # a_i = - omega^2 x_i (niezalezne oscylatory — jawny Jacobian)
        omega2 = 0.7
        return -omega2 * p

    A = acceleration_jacobian_positions(pos, C, acc_harmonic, eps=1e-4)
    assert A.shape == (6, 6)
    assert np.max(np.abs(A + 0.7 * np.eye(6))) < 1e-3

    lam, _st = largest_lyapunov_exponent_benettin(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        t_final=30.0,
        dt=0.02,
        renorm_every=30,
        jac_eps=1e-4,
        rng=rng,
    )
    assert abs(lam) < 0.05

    spec4, _ = lyapunov_spectrum_benettin(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        n_exponents=4,
        t_final=35.0,
        dt=0.02,
        renorm_every=35,
        jac_eps=1e-4,
        rng=np.random.default_rng(1),
    )
    assert np.max(np.abs(spec4)) < 0.09

    lam1, _s1 = lyapunov_spectrum_benettin(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        n_exponents=1,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(2),
    )
    lam1b, _s2 = largest_lyapunov_exponent_benettin(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(2),
    )
    assert abs(float(lam1[0]) - lam1b) < 0.04

    from . import dynamics_v2

    n3 = 3
    C3 = np.abs(rng.normal(size=n3)) + 0.4
    p3 = rng.normal(scale=0.8, size=(n3, 3))
    Gn = 1.0
    epsn = 1e-5

    def acc_n(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return dynamics_v2.forces_newton(p, c, G=Gn, softening=epsn) / c[:, None]

    Ja = acceleration_jacobian_newton_softened(p3, C3, G=Gn, softening=epsn)
    Jfd = acceleration_jacobian_positions(p3, C3, acc_n, eps=8e-6)
    assert np.max(np.abs(Ja - Jfd)) < 2e-4

    beta_t = 0.11
    gamma_t = beta_t
    soft_t = 1e-5

    def acc_tgp(p: np.ndarray, c: np.ndarray) -> np.ndarray:
        return dynamics_v2.analytical_forces_tgp_pairwise(
            p, c, beta_t, gamma_t, soft_t
        ) / c[:, None]

    Jt = acceleration_jacobian_tgp_pairwise_softened(
        p3, C3, beta=beta_t, gamma=gamma_t, softening=soft_t
    )
    Jtfd = acceleration_jacobian_positions(p3, C3, acc_tgp, eps=7e-6)
    assert np.max(np.abs(Jt - Jtfd)) < 3e-4

    from . import dynamics_backends

    acc_y, _ = dynamics_backends.build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta_t,
        gamma=gamma_t,
        softening=soft_t,
        include_3body=True,
        n_quad_feynman=14,
    )
    Jy_full = acceleration_jacobian_positions(p3, C3, acc_y, eps=7e-6)
    Jy_split = acceleration_jacobian_yukawa_feynman_split(
        p3,
        C3,
        beta=beta_t,
        gamma=gamma_t,
        softening=soft_t,
        n_quad_feynman=14,
        jac_eps=7e-6,
    )
    assert np.max(np.abs(Jy_split - Jy_full)) < 6e-3

    Jy_an = acceleration_jacobian_yukawa_feynman_analytic(
        p3,
        C3,
        beta=beta_t,
        gamma=gamma_t,
        softening=soft_t,
        n_quad_feynman=14,
    )
    assert np.max(np.abs(Jy_an - Jy_full)) < 5e-3

    lam_lf_h, _ = largest_lyapunov_exponent_benettin_leapfrog(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(11),
    )
    lam_rk_h, _ = largest_lyapunov_exponent_benettin(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(11),
    )
    assert abs(lam_lf_h) < 0.06
    assert abs(lam_rk_h) < 0.06
    assert abs(lam_lf_h - lam_rk_h) < 0.05

    lam_lf_b, _ = largest_lyapunov_exponent_benettin_leapfrog(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(12),
    )
    sp_lf, _ = lyapunov_spectrum_benettin_leapfrog(
        pos,
        np.zeros_like(pos),
        C,
        acc_harmonic,
        n_exponents=1,
        t_final=40.0,
        dt=0.02,
        renorm_every=40,
        jac_eps=1e-4,
        rng=np.random.default_rng(12),
    )
    assert abs(float(sp_lf[0]) - lam_lf_b) < 0.04


if __name__ == "__main__":
    _self_test()
    print("lyapunov: self-test OK")
