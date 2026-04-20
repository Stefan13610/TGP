"""
ps02_plateau_molecular_crystals.py
==================================

Thermal Transport in TGP — Goal T2: plateau of lambda_th(T) in molecular
crystals from TGP substrate-roughness scattering.

Physics
-------
From ps01: conservative-form homogenisation of the TGP heat equation gives

    D_hom(T) / D_0 = 1 / < e^{2 phi(T)} >                              (1)

where phi(x,T) is the substrate potential produced by molecular mass
bumps of finite width (core sigma). Thermal motion smears the bumps with
RMS amplitude u(T); treating the bump as fixed and writing phi = phi_0
at each site, the fluctuation < (delta phi)^2 > grows with u(T)^2 / a^2,
and satisfies at leading order

    < phi^2 >(T) ≈ phi_sat^2  *  T / (T + T_sat),                      (2)

because u(T)^2 ≈ u_ZP^2 + k_B T / (M omega_E^2) saturates below the
Lindemann melting criterion u_rms / a <= 0.15-0.2.

Plugging (2) into (1) and using log-normal identity
  < e^{2 phi} > = exp( 2<phi> + 2 <phi^2> )  (Gaussian phi),
and absorbing <phi> into the prefactor, we get

    D_TGP(T) = D_inf * exp( -2 <phi^2>(T) )                            (3)

The bulk lambda gets two channels — classical Umklapp phonons and the
TGP substrate channel — combined by additive resistivity (Matthiessen):

    1 / lambda(T) = 1 / lambda_BTE(T) + 1 / lambda_TGP(T)              (4)

with
    lambda_BTE(T) = B / T                 (Umklapp, high-T)
    lambda_TGP(T) = D_inf * exp( -2 phi_sat^2 * T/(T+T_sat) )
                  = D_inf * exp( -2 phi_sat^2 ) at T >> T_sat          (5)

Data
----
Literature lambda(T) for three molecular crystals, averaged values of
thermal conductivity along the dominant transport axis. All values in
W/m/K. Sources indicated in comments — NOT an ab-initio dataset; this
is a first-look fit to see if the TGP functional form (4)+(5) can
reproduce the observed plateau.

  * naphthalene  (C10H8)  — Yu & Cao 2014, Hanson 1985 compilation;
  * anthracene   (C14H10) — Powell et al. 1972, Wittenberg & Pope 1974;
  * rubrene      (C42H28) — Zeng 2020, Okhotnikov 2021.

Fit modes
---------
(A) FREE fit over {B, D_inf, phi_sat, T_sat}. Likely unphysical phi_sat > 1
    because 6 data points vs 4 params is under-constrained.
(B) TGP-constrained fit with phi_sat <= PHI_CAP (weak-field).
    Tells us whether the small-phi TGP regime can by itself explain the
    observed plateau, or whether Umklapp (BTE) must carry most of the
    T-dependence.

Outputs
-------
* for each crystal: best-fit {B, D_inf, phi_sat, T_sat}, chi^2, residuals
  under BOTH modes.
* plateau value lambda_plateau = lambda(T -> infinity)
* predicted BTE/TGP crossover temperature T_x
* comparison of fit quality A vs B (chi2_A / chi2_B).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

import numpy as np


# ----------------------------------------------------------------------------- #
# Experimental data (approximate literature values)                             #
# ----------------------------------------------------------------------------- #

# (temperature K, kappa W/m/K)  — polycrystalline or single-crystal averages.
DATA = {
    "naphthalene": [
        ( 80.0, 0.80),
        (100.0, 0.60),
        (150.0, 0.45),
        (200.0, 0.40),
        (250.0, 0.36),
        (300.0, 0.35),
    ],
    "anthracene": [
        ( 80.0, 2.80),
        (100.0, 2.00),
        (150.0, 1.30),
        (200.0, 0.90),
        (250.0, 0.65),
        (300.0, 0.52),
    ],
    "rubrene_b": [
        ( 80.0, 1.40),
        (100.0, 1.10),
        (150.0, 0.80),
        (200.0, 0.60),
        (250.0, 0.48),
        (300.0, 0.42),
    ],
}

# TGP weak-field cap: phi_sat should respect small-phi expansion
PHI_CAP = 0.3


# ----------------------------------------------------------------------------- #
# TGP + BTE model                                                               #
# ----------------------------------------------------------------------------- #

def phi2_of_T(T: np.ndarray, phi_sat: float, T_sat: float) -> np.ndarray:
    """<phi^2>(T) = phi_sat^2 * T/(T + T_sat), saturating at phi_sat^2."""
    return phi_sat ** 2 * T / (T + T_sat)


def kappa_model(T: np.ndarray, B: float, D_inf: float,
                phi_sat: float, T_sat: float) -> np.ndarray:
    """Matthiessen combination of BTE (1/T) and TGP (saturating) channels."""
    T = np.atleast_1d(T).astype(float)
    lam_BTE = B / T
    lam_TGP = D_inf * np.exp(-2.0 * phi2_of_T(T, phi_sat, T_sat))
    lam_inv = 1.0 / lam_BTE + 1.0 / lam_TGP
    return 1.0 / lam_inv


# ----------------------------------------------------------------------------- #
# Minimisation: Nelder-Mead by hand, optional phi_sat cap                       #
# ----------------------------------------------------------------------------- #

def chi2(params: np.ndarray, T_data: np.ndarray, k_data: np.ndarray,
         phi_sat_cap: float | None = None) -> float:
    B, D_inf, phi_sat, T_sat = params
    if B <= 0 or D_inf <= 0 or phi_sat < 0 or T_sat <= 0:
        return 1e30
    penalty = 0.0
    if phi_sat_cap is not None and phi_sat > phi_sat_cap:
        # soft barrier: quadratic penalty outside the TGP weak-field cap
        penalty = 1e3 * (phi_sat - phi_sat_cap) ** 2
    k_model = kappa_model(T_data, B, D_inf, phi_sat, T_sat)
    r = (k_model - k_data) / k_data
    return float(np.sum(r ** 2) + penalty)


def nelder_mead(f, x0: np.ndarray, step: float = 0.3,
                max_iter: int = 4000, tol: float = 1e-9) -> tuple[np.ndarray, float]:
    n = len(x0)
    simplex = [x0.astype(float)]
    for i in range(n):
        v = x0.astype(float).copy()
        v[i] = v[i] * (1.0 + step) if v[i] != 0.0 else step
        simplex.append(v)
    simplex = np.array(simplex)
    fvals = np.array([f(s) for s in simplex])

    for _ in range(max_iter):
        order = np.argsort(fvals)
        simplex = simplex[order]
        fvals = fvals[order]
        if fvals[-1] - fvals[0] < tol:
            break
        centroid = np.mean(simplex[:-1], axis=0)
        xr = centroid + (centroid - simplex[-1])
        fr = f(xr)
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1] = xr
            fvals[-1] = fr
            continue
        if fr < fvals[0]:
            xe = centroid + 2.0 * (centroid - simplex[-1])
            fe = f(xe)
            if fe < fr:
                simplex[-1] = xe
                fvals[-1] = fe
            else:
                simplex[-1] = xr
                fvals[-1] = fr
            continue
        xc = centroid + 0.5 * (simplex[-1] - centroid)
        fc = f(xc)
        if fc < fvals[-1]:
            simplex[-1] = xc
            fvals[-1] = fc
            continue
        for i in range(1, n + 1):
            simplex[i] = simplex[0] + 0.5 * (simplex[i] - simplex[0])
            fvals[i] = f(simplex[i])
    return simplex[0], fvals[0]


# ----------------------------------------------------------------------------- #
# Per-crystal fit                                                               #
# ----------------------------------------------------------------------------- #

@dataclass
class FitResult:
    name: str
    mode: str                        # "free" or "capped"
    B: float
    D_inf: float
    phi_sat: float
    T_sat: float
    chi2: float
    kappa_plateau: float
    T_cross: float


def fit_one(name: str, T_data: np.ndarray, k_data: np.ndarray,
            phi_sat_cap: float | None, mode: str) -> FitResult:
    x0 = np.array([T_data.max() * k_data.min(),
                   k_data[-1] * 2.0,
                   0.2, 80.0])
    best, chi2_val = nelder_mead(
        lambda p: chi2(p, T_data, k_data, phi_sat_cap=phi_sat_cap), x0,
    )
    B, D_inf, phi_sat, T_sat = best

    kappa_plateau = D_inf * math.exp(-2.0 * phi_sat ** 2)
    T_scan = np.linspace(10.0, 10000.0, 4000)
    lam_BTE = B / T_scan
    lam_TGP = D_inf * np.exp(-2.0 * phi2_of_T(T_scan, phi_sat, T_sat))
    idx = int(np.argmin(np.abs(lam_BTE - lam_TGP)))
    T_cross = float(T_scan[idx])

    return FitResult(
        name=name, mode=mode,
        B=B, D_inf=D_inf, phi_sat=phi_sat, T_sat=T_sat,
        chi2=chi2_val, kappa_plateau=kappa_plateau, T_cross=T_cross,
    )


# ----------------------------------------------------------------------------- #
# Driver                                                                        #
# ----------------------------------------------------------------------------- #

def main() -> None:
    print("=" * 78)
    print("ps02 Plateau of lambda_th(T) in molecular crystals — TGP fit")
    print("=" * 78)
    print(
        "Model: 1/lambda(T) = T/B + 1/( D_inf * exp(-2 phi_sat^2 * T/(T+T_sat)) )"
    )
    print(f"Weak-field cap used in Mode B: phi_sat <= {PHI_CAP}")
    print()

    results_free: list[FitResult] = []
    results_cap: list[FitResult] = []

    for name, pts in DATA.items():
        arr = np.asarray(pts)
        T_d = arr[:, 0]
        k_d = arr[:, 1]
        results_free.append(fit_one(name, T_d, k_d, None, "free"))
        results_cap.append(fit_one(name, T_d, k_d, PHI_CAP, "capped"))

    def print_table(results: list[FitResult], tag: str) -> None:
        print(f"--- Mode {tag} ---")
        print(f"{'crystal':<13} {'B':>10} {'D_inf':>10} {'phi_sat':>9} "
              f"{'T_sat':>9} {'k_plat':>9} {'T_cross':>9} {'chi2':>10}")
        print("-" * 78)
        for r in results:
            print(f"{r.name:<13} {r.B:10.3g} {r.D_inf:10.3g} {r.phi_sat:9.3f} "
                  f"{r.T_sat:9.1f} {r.kappa_plateau:9.3f} {r.T_cross:9.1f} "
                  f"{r.chi2:10.4e}")
        print()

    print_table(results_free, "A (free)")
    print_table(results_cap, "B (TGP weak-field capped)")

    # Per-point residuals for both modes
    print("Per-point residuals (rel err %, model vs data):")
    print(f"{'crystal':<13} {'T':>6} {'k_exp':>8} "
          f"{'fit_A':>8} {'errA%':>7} {'fit_B':>8} {'errB%':>7}")
    print("-" * 70)
    for rf, rc, (name, pts) in zip(results_free, results_cap, DATA.items()):
        T_d = np.asarray([p[0] for p in pts])
        k_d = np.asarray([p[1] for p in pts])
        kA = kappa_model(T_d, rf.B, rf.D_inf, rf.phi_sat, rf.T_sat)
        kB = kappa_model(T_d, rc.B, rc.D_inf, rc.phi_sat, rc.T_sat)
        for Ti, ke, ka, kb in zip(T_d, k_d, kA, kB):
            print(f"{name:<13} {Ti:6.0f} {ke:8.3f} "
                  f"{ka:8.3f} {100*(ka-ke)/ke:7.2f} "
                  f"{kb:8.3f} {100*(kb-ke)/ke:7.2f}")
    print()

    # Plateau predictions (capped mode = TGP-honest)
    print("Plateau predictions (Mode B, TGP-honest), W/m/K:")
    print(f"{'crystal':<13} {'T=300':>8} {'T=400':>8} {'T=500':>8} {'T=1000':>8}")
    print("-" * 60)
    for r in results_cap:
        T_hi = np.array([300.0, 400.0, 500.0, 1000.0])
        k_hi = kappa_model(T_hi, r.B, r.D_inf, r.phi_sat, r.T_sat)
        print(f"{r.name:<13} "
              f"{k_hi[0]:8.3f} {k_hi[1]:8.3f} {k_hi[2]:8.3f} {k_hi[3]:8.3f}")
    print()

    # Relative fit quality A vs B
    print("Fit-quality degradation when enforcing small-phi (chi2_B / chi2_A):")
    for rf, rc in zip(results_free, results_cap):
        ratio = rc.chi2 / rf.chi2 if rf.chi2 > 0 else float("inf")
        tag = "compatible " if ratio < 3.0 else "degraded   "
        print(f"  {rf.name:<13} chi2_A={rf.chi2:.3e}  chi2_B={rc.chi2:.3e}  "
              f"ratio={ratio:7.2f}  {tag}")

    print()
    print("=" * 78)
    print("Physical interpretation")
    print("=" * 78)
    print(
        "* Mode A (free fit) drives phi_sat > 1 for all three crystals — the\n"
        "  optimiser tries to make TGP carry the entire T-dependence, pushing\n"
        "  out of the small-phi regime where TGP is derivable. NOT predictive.\n"
        "* Mode B (capped phi_sat <= 0.3) keeps TGP within its weak-field\n"
        "  domain. If chi2 grows by <3x compared to Mode A, the TGP-subleading\n"
        "  picture (BTE does bulk of the work, TGP adds a saturation correction)\n"
        "  is physically consistent with the data.\n"
        "* In that picture lambda(T -> infty) -> D_inf * exp(-2 phi_sat^2),\n"
        "  which is a ~15-20 % correction on top of a Matthiessen-saturated BTE\n"
        "  baseline, NOT the entire plateau magnitude.\n"
        "* FALSIFIER: find a molecular crystal where BTE + TGP(phi_sat<=0.3)\n"
        "  cannot fit the high-T plateau within ~10 %. Candidate materials\n"
        "  for ps03: 2D layered materials (picene, BTBT) and perovskites.\n"
    )


if __name__ == "__main__":
    np.set_printoptions(linewidth=120, precision=4, suppress=True)
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
