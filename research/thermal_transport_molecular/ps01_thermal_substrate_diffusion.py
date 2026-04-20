"""
ps01_thermal_substrate_diffusion.py
===================================

Thermal Transport in TGP — Goal T1: heat diffusion in the exponential metric.

Physics
-------
TGP gives a conformal spatial metric around static matter distributions:

    g_ij(x) = e^{+2 phi(x)} delta_ij,    phi(x) = Phi(x) / Phi_0,

where phi is the dimensionless substrate potential, small (|phi| << 1)
outside of strong-field regions.

The heat equation for an internal-energy density u(x,t) on a curved spatial
slice is (Laplace-Beltrami):

    dt u = (1/sqrt(g)) di ( sqrt(g) g^{ij} dj u ) * D_0

In 3D with g_ij = Omega^2 delta_ij (Omega = e^phi):
    sqrt(g)  = Omega^3
    g^{ij}   = Omega^{-2} delta^{ij}

So
    dt u = D_0 * Omega^{-3} * di (Omega * di u)
         = D_0 * [ Omega^{-2} Laplacian(u) + Omega^{-3} (di Omega)(di u) ]
         = D_0 * e^{-2 phi} [ Laplacian(u) + (di phi)(di u) ]     (1)

In 1D (the case studied here), equation (1) becomes
    dt u = D_0 * e^{-2 phi(x)} [ uxx + phi'(x) ux ]              (2)
         = D_0 * e^{-3 phi(x)} * d/dx [ e^{phi(x)} * du/dx ]     (2')

Form (2') is conservative and what we discretise.

Effective local diffusivity
---------------------------
Reading (2) as d_t u = D_eff(x) * uxx + v(x) * ux with
    D_eff(x) = D_0 * e^{-2 phi(x)}                                (3)
    v(x)     = D_0 * e^{-2 phi(x)} * phi'(x)                      (4)

We have two TGP effects:

  * amplitude effect: regions of *higher* Phi (denser substrate) diffuse
    heat **slower** — exactly the "molecular island" picture;
  * gradient effect: a phi-gradient produces a drift v(x). This is the
    new TGP scattering channel absent in BTE/phonon pictures.

Long-wavelength homogenisation
------------------------------
For a periodic phi(x) with period a, Fick's-law homogenisation of (2')
yields the macroscopic effective diffusivity

    1/D_hom = < e^{2 phi} / D_0 > = (1/D_0) * < e^{2 phi} >

but see the derivation below for the correct weighted average coming from
the conservative form. Concretely we test this numerically by fitting
the MSD < x^2 >(t) of an initial delta-peak evolved in (2').

Numerics
--------
* 1D lattice of N_sites molecular centres at x_n = n * a.
* Substrate potential built as a superposition of Gaussian bumps:
      phi(x) = phi_amp * sum_n exp( -(x - x_n)^2 / (2 sigma^2) )
  so max(phi) ~ phi_amp * (1 + 2 e^{-a^2/(2 sigma^2)}) and the lattice
  is resolved when sigma < a/2.
* Heat equation (2') discretised by conservative finite volumes on a
  uniform grid of spacing dx << sigma. Explicit Euler with
      dt < 0.5 * dx^2 / max D_eff(x).
* Periodic boundary conditions.

Outputs
-------
* D_eff_MSD   — effective D extracted from < x^2 >(t).
* D_hom_pred  — harmonic-mean prediction.
* scattering length l_scat from gradient term.
* scan over phi_amp: D_eff/D_0 vs <phi^2>.

Usage
-----
    python ps01_thermal_substrate_diffusion.py

Outputs a summary table; no external figure files.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

import numpy as np


# ----------------------------------------------------------------------------- #
# Lattice + substrate field construction                                        #
# ----------------------------------------------------------------------------- #

@dataclass
class Lattice1D:
    """Geometry and discretisation of the 1D molecular crystal."""
    N_sites: int          # number of molecular sites
    a: float              # lattice spacing (nm)
    sigma: float          # Gaussian width of each molecular "mass bump" (nm)
    phi_amp: float        # amplitude of phi at one molecular site
    points_per_cell: int  # grid resolution per lattice period

    @property
    def L(self) -> float:
        return self.N_sites * self.a

    @property
    def dx(self) -> float:
        return self.a / self.points_per_cell

    def x_grid(self) -> np.ndarray:
        return np.arange(0.0, self.L, self.dx)

    def site_positions(self) -> np.ndarray:
        return 0.5 * self.a + np.arange(self.N_sites) * self.a


def build_phi(lat: Lattice1D) -> np.ndarray:
    """
    Compose phi(x) as a sum of periodic Gaussian bumps centred on each
    molecular site. Periodic image sum via the k=-1..+1 neighbour cells
    is enough for sigma << L.
    """
    x = lat.x_grid()
    phi = np.zeros_like(x)
    sites = lat.site_positions()
    for x_n in sites:
        for shift in (-lat.L, 0.0, +lat.L):
            phi += lat.phi_amp * np.exp(
                -((x - x_n - shift) ** 2) / (2.0 * lat.sigma ** 2)
            )
    return phi


# ----------------------------------------------------------------------------- #
# Conservative finite-volume solver for eq. (2')                                #
# ----------------------------------------------------------------------------- #

def step_heat_conservative(
    u: np.ndarray,
    exp_phi: np.ndarray,
    exp_neg3phi: np.ndarray,
    D0: float,
    dx: float,
    dt: float,
) -> np.ndarray:
    """
    One explicit step of:

        dt u = D0 * e^{-3 phi} * d/dx ( e^{phi} * du/dx )

    Conservative form, periodic BC.
    Face fluxes computed at i+1/2 using
        flux_{i+1/2} = 0.5*(e^phi_i + e^phi_{i+1}) * (u_{i+1} - u_i)/dx
    then
        u_i^{new} = u_i + dt * D0 * e^{-3 phi_i} * (flux_{i+1/2} - flux_{i-1/2})/dx
    """
    ep_face = 0.5 * (exp_phi + np.roll(exp_phi, -1))   # face i+1/2
    grad = (np.roll(u, -1) - u) / dx                    # face i+1/2
    flux = ep_face * grad                               # face i+1/2
    div = (flux - np.roll(flux, +1)) / dx               # cell i
    return u + dt * D0 * exp_neg3phi * div


def solve_heat(
    lat: Lattice1D,
    D0: float,
    t_total: float,
    CFL: float = 0.25,
    u0_width: float = 0.3,
) -> dict:
    """
    Evolve a narrow Gaussian initial condition under TGP heat equation (2')
    and under classical (phi=0) heat equation for reference. Returns the
    MSD growth and extracted effective diffusivities.
    """
    x = lat.x_grid()
    phi = build_phi(lat)
    exp_phi = np.exp(phi)
    exp_neg3phi = np.exp(-3.0 * phi)

    # CFL step size based on max local D_eff
    D_eff_max = D0 * np.max(np.exp(-2.0 * phi))
    dt = CFL * lat.dx ** 2 / D_eff_max
    n_steps = max(1, int(math.ceil(t_total / dt)))
    dt = t_total / n_steps

    # Initial condition: narrow Gaussian, centre of the box, unit area
    x0 = 0.5 * lat.L
    u_tgp = np.exp(-((x - x0) ** 2) / (2.0 * u0_width ** 2))
    u_tgp /= np.trapezoid(u_tgp, x)
    u_cl = u_tgp.copy()

    # Sample MSD at a few snapshots
    snapshots_t = np.linspace(0.0, t_total, 11)
    snapshots_t = snapshots_t[1:]            # skip t=0
    snap_idx = np.unique(
        np.clip(np.round(snapshots_t / dt).astype(int), 1, n_steps)
    )
    msd_tgp, msd_cl, var0 = [], [], None

    for step in range(1, n_steps + 1):
        u_tgp = step_heat_conservative(u_tgp, exp_phi, exp_neg3phi, D0, lat.dx, dt)
        # Classical reference: phi = 0 ⇒ exp_phi=1, exp_neg3phi=1
        u_cl = step_heat_conservative(
            u_cl, np.ones_like(phi), np.ones_like(phi), D0, lat.dx, dt,
        )
        if step in snap_idx:
            var_tgp = np.trapezoid((x - x0) ** 2 * u_tgp, x) / np.trapezoid(u_tgp, x)
            var_cl = np.trapezoid((x - x0) ** 2 * u_cl, x) / np.trapezoid(u_cl, x)
            if var0 is None:
                var0 = u0_width ** 2
            msd_tgp.append((step * dt, var_tgp))
            msd_cl.append((step * dt, var_cl))

    msd_tgp = np.asarray(msd_tgp)
    msd_cl = np.asarray(msd_cl)

    # Linear fit MSD(t) = 2 * D_eff * t + const
    def _fit_D(arr):
        t = arr[:, 0]
        v = arr[:, 1]
        # last 60 % of points for asymptotic regime
        k = max(3, int(0.4 * len(t)))
        slope, _ = np.polyfit(t[k:], v[k:], 1)
        return 0.5 * slope

    D_eff_MSD = _fit_D(msd_tgp)
    D_eff_cl = _fit_D(msd_cl)

    # Analytical homogenisation predictions for eq (2')
    # Conservative form ⇒ harmonic-mean of e^{2 phi}:
    D_hom_harmonic = D0 / np.mean(np.exp(2.0 * phi))
    # Naive Jensen arithmetic-mean (wrong but instructive):
    D_hom_arith = D0 * np.mean(np.exp(-2.0 * phi))

    return dict(
        dt=dt, n_steps=n_steps,
        D_eff_MSD=D_eff_MSD,
        D_eff_cl=D_eff_cl,
        D_hom_harmonic=D_hom_harmonic,
        D_hom_arith=D_hom_arith,
        phi=phi, x=x,
        msd_tgp=msd_tgp, msd_cl=msd_cl,
    )


# ----------------------------------------------------------------------------- #
# Gradient-scattering length scale                                              #
# ----------------------------------------------------------------------------- #

def gradient_scattering_length(lat: Lattice1D) -> dict:
    """
    Characteristic scattering length from the advection term
        v(x) = D0 * e^{-2 phi} * phi'(x)
    Define l_scat such that v * l_scat / D_eff ~ 1, i.e.
        l_scat = 1 / |phi'(x)|_rms.
    """
    x = lat.x_grid()
    phi = build_phi(lat)
    dphi = np.gradient(phi, lat.dx)
    dphi_rms = float(np.sqrt(np.mean(dphi ** 2)))
    dphi_max = float(np.max(np.abs(dphi)))
    return dict(
        dphi_rms=dphi_rms,
        dphi_max=dphi_max,
        l_scat_rms=1.0 / dphi_rms if dphi_rms > 0 else math.inf,
        l_scat_max=1.0 / dphi_max if dphi_max > 0 else math.inf,
    )


# ----------------------------------------------------------------------------- #
# Main experiment: scan over phi_amp                                            #
# ----------------------------------------------------------------------------- #

def run_scan() -> None:
    D0 = 1.0        # dimensionless (nm^2 / arbitrary-time)
    base = Lattice1D(
        N_sites=32,
        a=1.0,                 # 1 nm lattice spacing
        sigma=0.25,            # molecular "core" radius
        phi_amp=0.0,           # overwritten in the scan
        points_per_cell=32,
    )
    t_total = 12.0             # long enough for diffusion to sample many cells

    print("=" * 78)
    print("ps01 Thermal Substrate Diffusion — TGP heat equation in 1D crystal")
    print("=" * 78)
    print(f"Lattice : N_sites={base.N_sites}, a={base.a} nm, sigma={base.sigma} nm")
    print(f"Grid    : dx={base.dx:.4f} nm,  L={base.L} nm,  D0={D0}")
    print(f"Time    : t_total={t_total}")
    print()

    phi_amp_list = [0.0, 0.01, 0.03, 0.1, 0.3]
    print(f"{'phi_amp':>8} {'<phi>':>9} {'<phi^2>':>10} "
          f"{'D_MSD/D0':>10} {'D_harm/D0':>10} {'D_arith/D0':>11} "
          f"{'l_scat (nm)':>12}")
    print("-" * 78)

    rows = []
    for amp in phi_amp_list:
        lat = Lattice1D(
            N_sites=base.N_sites,
            a=base.a, sigma=base.sigma,
            phi_amp=amp,
            points_per_cell=base.points_per_cell,
        )
        res = solve_heat(lat, D0, t_total)
        phi = res["phi"]
        phi_mean = float(np.mean(phi))
        phi_var = float(np.mean(phi ** 2))
        gs = gradient_scattering_length(lat)

        row = dict(
            phi_amp=amp,
            phi_mean=phi_mean,
            phi_var=phi_var,
            D_MSD=res["D_eff_MSD"],
            D_cl=res["D_eff_cl"],
            D_harm=res["D_hom_harmonic"],
            D_arith=res["D_hom_arith"],
            l_scat=gs["l_scat_rms"],
        )
        rows.append(row)

        print(f"{amp:8.3f} {phi_mean:9.4f} {phi_var:10.5f} "
              f"{row['D_MSD']/D0:10.5f} {row['D_harm']/D0:10.5f} "
              f"{row['D_arith']/D0:11.5f} {row['l_scat']:12.4f}")

    print()
    print("=" * 78)
    print("Notes")
    print("=" * 78)
    print(
        "* D_MSD is the numerical effective diffusivity from MSD(t) = 2 D_eff t.\n"
        "* D_harm = D0 / <e^{2phi}> is the analytic prediction from the\n"
        "  conservative form of (2'). Numerics should track it within 1-2 %.\n"
        "* D_arith = D0 * <e^{-2phi}> is the naive Jensen-arithmetic mean.\n"
        "  It deviates from the true homogenised result and gets worse at\n"
        "  larger phi_amp — a direct test of (2') vs ad-hoc averaging.\n"
        "* l_scat = 1/|phi'|_rms is the gradient-scattering length. At phi~0.3\n"
        "  and sigma=0.25 nm it falls below a (sub-nm), signalling that the\n"
        "  TGP advection term becomes a genuine new scattering channel."
    )

    # --- Physical sanity check: numerics match harmonic-mean homogenisation  -
    print()
    print("Homogenisation residuals (|D_MSD - D_harm|/D_harm):")
    for r in rows:
        rel = abs(r["D_MSD"] - r["D_harm"]) / r["D_harm"]
        tag = "OK " if rel < 0.03 else "!! "
        print(f"  phi_amp={r['phi_amp']:.3f} :  {rel*100:7.3f} %  {tag}")

    # --- Connection to molecular-crystal plateau (T2 teaser) -----------------
    print()
    print("Connection to T2 (plateau at high T):")
    print(
        "  In molecular crystals the relevant 'phi_amp' grows with thermal\n"
        "  rattle amplitude of molecules around lattice sites. Above a\n"
        "  characteristic T, <delta phi^2> saturates because the bumps are\n"
        "  bounded by the molecular size sigma. Hence D_harm saturates too —\n"
        "  this is the TGP explanation of the observed plateau lambda_th(T).\n"
        "  ps02 will fit naphthalene/anthracene/rubren thermal data using\n"
        "  exactly the same D_harm formula used here."
    )


if __name__ == "__main__":
    np.set_printoptions(linewidth=120, precision=6, suppress=True)
    try:
        run_scan()
    except KeyboardInterrupt:
        sys.exit(130)
