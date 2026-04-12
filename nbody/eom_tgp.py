"""
eom_tgp.py — Analityczna postać równań ruchu N ciał TGP (warstwa 1 vs warstwa 2)
================================================================================

Warstwa 1 — ANalityka TGP (bez czasu, bez kwadratur po czasie)
---------------------------------------------------------------
Przy źródłach punktowych i Drodze B (profil Yukawy jako warunek brzegowy) energia
potencjalna w zadanym rzędzie (2-ciałowy + nieredukowalny 3-ciałowy z potencjału)
ma postać:

    V(q) = Σ_{i<j}  V_2(d_ij; C_i,C_j,β,γ)
         + Σ_{i<j<k} V_3^{irr}(d_ij, d_ik, d_jk; C_i,C_j,C_k,β,γ)

    V_3^{irr} = (2β - 6γ) C_i C_j C_k · I(d_ij, d_ik, d_jk)

gdzie I jest całką nakładki potrójnej (Yukawa lub jej granica Coulomba).  Siła:

    F_i = -∇_{x_i} V  =  (2-ciałowy, zamknięty)  +  Σ_{j<k, j,k≠i} F_i^{(ijk)}

Dla pojedynczego tripletu (i,j,k) zależność I tylko od odległości daje:

    ∂I/∂x_i = (∂I/∂d_ij) (x_i-x_j)/d_ij + (∂I/∂d_ik) (x_i-x_k)/d_ik

    F_i^{(ijk)} = -∂V_3/∂x_i
                = (6γ - 2β) C_i C_j C_k · ∂I/∂x_i

(analogicznie dla j, k — indeksy w parach z i).  To jest **pełna postać wektorowa**
niezależna od wyboru backendu dla I.

Aksjomat TGP: m_i = C_i.  Równania ruchu (Newton):

    ẍ_i = (1/C_i) F_i(q).

Warstwa 2 — Przeliczenie na współrzędne (kwadratury / całka po czasie)
----------------------------------------------------------------------
- Wyznaczenie ∂I/∂d_ab **może** wymagać kwadratury 2D (Feynman, moduł
  `three_body_force_exact`) albo siatki 3D (`triple_overlap_numerical`), albo
  być dane zamknięcie (Coulomb: I = 8π²/P, P = d_ij+d_ik+d_jk).
- Całka po czasie (pozycje w t) to osobny problem: `dynamics_v2.leapfrog_integrate`
  itd. — symplektyczny krok wymaga tylko F_i(q) lub a_i = F_i/C_i.

Ten moduł implementuje **warstwę 1** w postaci sumowania tripletów z dostawcą
pochodnych ∂I/∂d (callable).  Nie duplikuje całki Feynmana — użyj
`three_body_force_exact.three_body_forces_exact` gdy I = I_Yukawa dokładne.
"""

from __future__ import annotations

from itertools import combinations
from typing import Callable, Tuple

import numpy as np

try:
    from . import dynamics_v2
except ImportError:  # uruchomienie `python eom_tgp.py` z katalogu nbody/
    import dynamics_v2

# Typ: dla danego tripletu zwraca (∂I/∂d_ij, ∂I/∂d_ik, ∂I/∂d_jk)
TripletIDerivatives = Callable[[float, float, float], Tuple[float, float, float]]


def irreducible_three_body_forces_from_I_derivatives(
    positions: np.ndarray,
    C_values: np.ndarray,
    beta: float,
    gamma: float,
    softening: float,
    dI_triplet: TripletIDerivatives,
) -> np.ndarray:
    """
    Siły 3-ciałowe z
        V_3 = (2β - 6γ) C_i C_j C_k I(d_ij,d_ik,d_jk)
    przy znanych pochodnych I.

    Parameters
    ----------
    positions : (N, 3)
    C_values : (N,)
    beta, gamma : float
    softening : regulator numeryczny na d (jak w dynamics_v2)
    dI_triplet : callable (d_ij, d_ik, d_jk) -> (g_ij, g_ik, g_jk)
        g_ab = ∂I/∂d_ab w punktach aktualnych odległości.

    Returns
    -------
    forces : (N, 3) — wektory siły (nie przyspieszenia).
    """
    n = len(C_values)
    forces = np.zeros((n, 3))
    if n < 3:
        return forces

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_values[i], C_values[j], C_values[k]
        pref = (6.0 * gamma - 2.0 * beta) * Ci * Cj * Ck

        rij = positions[j] - positions[i]
        rik = positions[k] - positions[i]
        rjk = positions[k] - positions[j]

        dij = float(np.sqrt(np.dot(rij, rij) + softening**2))
        dik = float(np.sqrt(np.dot(rik, rik) + softening**2))
        djk = float(np.sqrt(np.dot(rjk, rjk) + softening**2))

        g_ij, g_ik, g_jk = dI_triplet(dij, dik, djk)

        u_ij = -rij / dij  # (x_i - x_j) / d_ij
        u_ik = -rik / dik
        u_ji = rij / dij
        u_jk = -rjk / djk
        u_ki = rik / dik
        u_kj = rjk / djk

        forces[i] += pref * (g_ij * u_ij + g_ik * u_ik)
        forces[j] += pref * (g_ij * u_ji + g_jk * u_jk)
        forces[k] += pref * (g_ik * u_ki + g_jk * u_kj)

    return forces


def dI_coulomb_closed_form(d_ij: float, d_ik: float, d_jk: float) -> Tuple[float, float, float]:
    """
    Granica Coulomba: I = 8π²/P, P = d_ij + d_ik + d_jk.

    ∂I/∂d_ij = ∂I/∂P = -8π²/P² (wszystkie trzy pochodne równe).
    """
    P = d_ij + d_ik + d_jk
    if P <= 0.0:
        return 0.0, 0.0, 0.0
    g = -8.0 * np.pi**2 / (P * P)
    return float(g), float(g), float(g)


def accelerations_tgp_nbody(
    positions: np.ndarray,
    C_values: np.ndarray,
    beta: float,
    gamma: float | None = None,
    softening: float = 1e-6,
    *,
    include_3body: bool = True,
    three_body_dI: TripletIDerivatives | None = None,
) -> np.ndarray:
    """
    Przyspieszenia ẍ_i = F_i / C_i przy rozbiciu 2B (zawsze zamknięte) + 3B z dostawcą ∂I/∂d.

    Jeśli ``three_body_dI is None`` i ``include_3body``, używane jest zamknięcie Coulomba
    na I (spójne z ``dynamics_v2.three_body_forces_approximate(..., use_yukawa=False)``).

    Dla dokładnej Yukawy podaj backend z ``three_body_force_exact`` zamiast tej funkcji,
    albo rozszerz ``three_body_dI`` o kwadraturę pochodnych I_Y.
    """
    if gamma is None:
        gamma = beta

    F2 = dynamics_v2.analytical_forces_tgp_pairwise(
        positions, C_values, beta, gamma, softening
    )
    if not include_3body or len(C_values) < 3:
        return F2 / C_values[:, None]

    dI_fn = three_body_dI if three_body_dI is not None else dI_coulomb_closed_form
    F3 = irreducible_three_body_forces_from_I_derivatives(
        positions, C_values, beta, gamma, softening, dI_fn
    )
    return (F2 + F3) / C_values[:, None]


def _max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)))


def _self_test_coulomb_agreement() -> None:
    """Porównanie z dynamics_v2 przy Coulombie (use_yukawa=False)."""
    rng = np.random.default_rng(0)
    n = 5
    pos = rng.normal(size=(n, 3))
    C = np.abs(rng.normal(size=n)) + 0.2
    beta = gamma = 1.0
    eps = 1e-6

    F_legacy = dynamics_v2.full_forces_tgp(
        pos, C, beta, gamma, eps, include_3body=True, use_yukawa=False
    )
    acc_new = accelerations_tgp_nbody(
        pos, C, beta, gamma, eps, include_3body=True, three_body_dI=None
    )
    F_new = acc_new * C[:, None]
    d = _max_abs_diff(F_legacy, F_new)
    if d > 1e-9:
        raise AssertionError(f"Coulomb 3B force mismatch: max diff {d}")


if __name__ == "__main__":
    _self_test_coulomb_agreement()
    print("eom_tgp: self-test Coulomb vs dynamics_v2 OK  (albo: python -m nbody.eom_tgp z TGP_v1/)")
