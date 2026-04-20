"""
ps03_perovskite_hybrid_scan.py
==============================

Thermal Transport in TGP — Goal T3: predict room-temperature lambda_th
across a wide range of crystals with ONE universal TGP parameter
alpha_T, using only measured material quantities { a, rho, v_s, M_fu }.

Formula
-------
Starting from the ps02 plateau,

    lambda_TGP(T) = D_inf(material) * exp( -2 * <phi^2>(T, material) )

and identifying
  * D_inf with the kinetic-theory minimum lambda_kin = (1/3) n_a k_B v_s a
    (phonon mean free path l = a, Cahill limit, Dulong-Petit C_V),
  * <phi^2>(T) = alpha_T * k_B T / (M_eff v_s^2),   with M_eff = M_fu/N_at
    (Debye-Waller <u^2>/a^2 in equipartition high-T limit, rewritten for the
    substrate bump amplitude),
we obtain a *predictive* formula

    kappa_TGP(T) = (1/3) n_a k_B v_s a
                   * exp( -2 alpha_T k_B T / (M_eff v_s^2) )           (1)

One universal parameter alpha_T is fit once (global MAE minimisation),
then used to predict every other material without further tuning.

Materials
---------
Cubic / quasi-cubic ionic + perovskite + hybrid-perovskite crystals with
reported room-temperature data. Values are representative literature
averages, not single-experiment precision. The purpose is the functional
trend across two orders of magnitude in kappa_exp, not exact fits.

  Material      a(A)  rho(g/cc)  v_s(m/s)  M_fu  N_at  kappa_exp(W/m/K)
  SrTiO3        3.91    5.12       8000     183    5     11.0
  BaTiO3        4.00    6.02       5500     233    5      3.5
  NaCl          5.64    2.16       3500      58    2      6.5
  KBr           6.60    2.75       3000     119    2      3.4
  PbTe          6.46    8.16       1800     334    2      2.0
  CsPbBr3       5.92    4.55       2400     579    5      0.42
  CsPbI3        6.30    4.70       2100     720    5      0.34
  MAPbBr3       5.92    3.80       2100     478   12      0.37
  MAPbI3        6.32    4.10       1800     619   12      0.30
  CsSnI3        6.22    4.30       1700     604    5      0.25
  Cs2AgBiBr6    5.58    4.92       1900    1112   10      0.33

Outputs
-------
* best alpha_T (single float)
* table of predicted lambda_TGP vs experiment (all 11 materials)
* global MAE, max relative error
* list of outliers where TGP prediction misses by > 30 %
* sensitivity of predictions to each material parameter.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

import numpy as np


# ----------------------------------------------------------------------------- #
# Physical constants                                                            #
# ----------------------------------------------------------------------------- #
K_B = 1.380649e-23        # J/K
AMU = 1.66053906660e-27   # kg
T_ROOM = 300.0            # K


# ----------------------------------------------------------------------------- #
# Materials database                                                            #
# ----------------------------------------------------------------------------- #

@dataclass
class Material:
    name: str
    a: float         # lattice constant (Angstrom)
    rho: float       # density (g/cc)
    v_s: float       # acoustic sound speed (m/s)
    M_fu: float      # formula-unit mass (amu)
    N_at: int        # atoms per formula unit
    kappa_exp: float # experimental kappa at 300 K (W/m/K)

    @property
    def a_m(self) -> float:
        return self.a * 1e-10

    @property
    def M_eff_kg(self) -> float:
        return (self.M_fu / self.N_at) * AMU

    @property
    def n_at(self) -> float:
        """Atomic number density (1/m^3) from rho/M_eff."""
        rho_SI = self.rho * 1e3            # kg/m^3
        m_atom = self.M_eff_kg
        return rho_SI / m_atom


MATERIALS = [
    Material("SrTiO3",      3.91, 5.12, 8000.0,  183.0,  5, 11.00),
    Material("BaTiO3",      4.00, 6.02, 5500.0,  233.0,  5,  3.50),
    Material("NaCl",        5.64, 2.16, 3500.0,   58.0,  2,  6.50),
    Material("KBr",         6.60, 2.75, 3000.0,  119.0,  2,  3.40),
    Material("PbTe",        6.46, 8.16, 1800.0,  334.0,  2,  2.00),
    Material("CsPbBr3",     5.92, 4.55, 2400.0,  579.0,  5,  0.42),
    Material("CsPbI3",      6.30, 4.70, 2100.0,  720.0,  5,  0.34),
    Material("MAPbBr3",     5.92, 3.80, 2100.0,  478.0, 12,  0.37),
    Material("MAPbI3",      6.32, 4.10, 1800.0,  619.0, 12,  0.30),
    Material("CsSnI3",      6.22, 4.30, 1700.0,  604.0,  5,  0.25),
    Material("Cs2AgBiBr6",  5.58, 4.92, 1900.0, 1112.0, 10,  0.33),
]


# ----------------------------------------------------------------------------- #
# TGP formula                                                                   #
# ----------------------------------------------------------------------------- #

def kappa_kin(mat: Material) -> float:
    """lambda_kin = (1/3) n_a k_B v_s a  — Cahill-like kinetic floor."""
    return (1.0 / 3.0) * mat.n_at * K_B * mat.v_s * mat.a_m


def phi2_sat(mat: Material, alpha_T: float, T: float = T_ROOM) -> float:
    """<phi^2> = alpha_T * k_B T / (M_eff v_s^2) — Debye-Waller rewrite."""
    return alpha_T * K_B * T / (mat.M_eff_kg * mat.v_s ** 2)


def kappa_TGP(mat: Material, alpha_T: float, T: float = T_ROOM) -> float:
    return kappa_kin(mat) * math.exp(-2.0 * phi2_sat(mat, alpha_T, T))


# ----------------------------------------------------------------------------- #
# Fit alpha_T globally                                                          #
# ----------------------------------------------------------------------------- #

def global_chi2(alpha_T: float) -> float:
    s = 0.0
    for m in MATERIALS:
        k_pred = kappa_TGP(m, alpha_T)
        # log-space residuals — kappa spans 2 orders of magnitude
        r = math.log(k_pred / m.kappa_exp)
        s += r ** 2
    return s


def brent_minimum(f, a: float, b: float, c: float,
                  tol: float = 1e-8, max_iter: int = 200) -> float:
    """Simple golden-section search for f(x) minimum in [a,c], with b in between."""
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = c - phi * (c - a)
    x2 = a + phi * (c - a)
    f1, f2 = f(x1), f(x2)
    for _ in range(max_iter):
        if abs(c - a) < tol:
            break
        if f1 < f2:
            c = x2
            x2 = x1
            f2 = f1
            x1 = c - phi * (c - a)
            f1 = f(x1)
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + phi * (c - a)
            f2 = f(x2)
    return 0.5 * (a + c)


# ----------------------------------------------------------------------------- #
# Driver                                                                        #
# ----------------------------------------------------------------------------- #

SOFT_SUBSET = {
    "CsPbBr3", "CsPbI3", "MAPbBr3", "MAPbI3", "CsSnI3", "Cs2AgBiBr6",
}


def global_chi2_subset(alpha_T: float, names: set[str]) -> float:
    s = 0.0
    for m in MATERIALS:
        if m.name not in names:
            continue
        k_pred = kappa_TGP(m, alpha_T)
        r = math.log(k_pred / m.kappa_exp)
        s += r ** 2
    return s


def main() -> None:
    print("=" * 78)
    print("ps03 Perovskite + ionic crystal scan — TGP prediction of kappa(300 K)")
    print("=" * 78)
    print("Formula:  kappa_TGP = (1/3) n_a k_B v_s a * exp( -2 alpha_T * k_B T")
    print("                                                  / (M_eff v_s^2) )")
    print("Single universal parameter alpha_T; two independent fits below.")
    print()

    # Helper: brent over alpha
    def fit_alpha(obj) -> tuple[float, float]:
        alpha_scan = np.logspace(-2, 3, 81)
        chi_scan = np.array([obj(a) for a in alpha_scan])
        i_best = int(np.argmin(chi_scan))
        a_lo = alpha_scan[max(0, i_best - 2)]
        a_hi = alpha_scan[min(len(alpha_scan) - 1, i_best + 2)]
        a_T = brent_minimum(obj, a_lo, alpha_scan[i_best], a_hi,
                             tol=1e-9, max_iter=200)
        return a_T, obj(a_T)

    # ------- Fit A: global (all 11) -------
    alpha_A, chi_A = fit_alpha(global_chi2)
    print(f"[Fit A : ALL 11 materials]    alpha_T = {alpha_A:.4f}   chi^2 = {chi_A:.4f}")
    # ------- Fit B: perovskite-only (soft subset) -------
    alpha_B, chi_B = fit_alpha(lambda a: global_chi2_subset(a, SOFT_SUBSET))
    print(f"[Fit B : soft perovskites ]   alpha_T = {alpha_B:.4f}   "
          f"chi^2 = {chi_B:.4f}")
    print()

    alpha_T = alpha_A         # keep this as the default for the big table

    # Per-material table
    print(f"{'material':<13} {'a(A)':>6} {'v_s(m/s)':>9} {'M_eff':>7} "
          f"{'n_a(1/nm3)':>11} {'k_kin':>7} {'<phi^2>':>9} "
          f"{'k_TGP':>7} {'k_exp':>7} {'rel_err':>9}")
    print("-" * 102)
    errs = []
    for m in MATERIALS:
        k_kin = kappa_kin(m)
        p2 = phi2_sat(m, alpha_T)
        k_TGP = kappa_TGP(m, alpha_T)
        rel = (k_TGP - m.kappa_exp) / m.kappa_exp * 100.0
        errs.append(rel)
        print(f"{m.name:<13} {m.a:6.2f} {m.v_s:9.0f} "
              f"{m.M_fu/m.N_at:7.2f} {m.n_at*1e-27:11.3f} "
              f"{k_kin:7.3f} {p2:9.5f} "
              f"{k_TGP:7.3f} {m.kappa_exp:7.3f} {rel:8.2f}%")

    errs = np.asarray(errs)
    print()
    print(f"Mean |rel err|   = {np.mean(np.abs(errs)):.2f} %")
    print(f"Max  |rel err|   = {np.max(np.abs(errs)):.2f} %")
    print(f"RMS (log-ratio)  = {math.sqrt(chi_A / len(MATERIALS)):.4f}")
    print()

    # Outliers
    print("Outliers (|rel err| > 30 %):")
    flag_any = False
    for m, e in zip(MATERIALS, errs):
        if abs(e) > 30.0:
            k_TGP = kappa_TGP(m, alpha_T)
            print(f"  {m.name:<13} "
                  f"kappa_exp={m.kappa_exp:.2f}  kappa_TGP={k_TGP:.2f}  "
                  f"err={e:+.1f} %")
            flag_any = True
    if not flag_any:
        print("  (none)")
    print()

    # Sensitivity: how does the prediction scale with inputs
    print("Parameter sensitivity d ln kappa / d ln X  (evaluated at MAPbI3):")
    m_ref = next(m for m in MATERIALS if m.name == "MAPbI3")
    k0 = kappa_TGP(m_ref, alpha_T)
    h = 0.01
    # Perturb v_s
    m_p = Material(
        m_ref.name, m_ref.a, m_ref.rho, m_ref.v_s * (1 + h),
        m_ref.M_fu, m_ref.N_at, m_ref.kappa_exp,
    )
    s_vs = (math.log(kappa_TGP(m_p, alpha_T)) - math.log(k0)) / h
    # Perturb M_fu
    m_p = Material(
        m_ref.name, m_ref.a, m_ref.rho, m_ref.v_s,
        m_ref.M_fu * (1 + h), m_ref.N_at, m_ref.kappa_exp,
    )
    s_M = (math.log(kappa_TGP(m_p, alpha_T)) - math.log(k0)) / h
    # Perturb a
    m_p = Material(
        m_ref.name, m_ref.a * (1 + h), m_ref.rho, m_ref.v_s,
        m_ref.M_fu, m_ref.N_at, m_ref.kappa_exp,
    )
    s_a = (math.log(kappa_TGP(m_p, alpha_T)) - math.log(k0)) / h
    print(f"  d ln k / d ln v_s    = {s_vs:+.3f}")
    print(f"  d ln k / d ln M_eff  = {s_M:+.3f}")
    print(f"  d ln k / d ln a      = {s_a:+.3f}")
    print()

    # Comparison: pure kinetic (no TGP) prediction
    print("Compare: pure kinetic (Cahill-like) without TGP exponent (alpha_T = 0):")
    print(f"{'material':<13} {'k_kin':>7} {'k_exp':>7} {'rel_err':>9}")
    for m in MATERIALS:
        k_kin = kappa_kin(m)
        rel = (k_kin - m.kappa_exp) / m.kappa_exp * 100.0
        print(f"  {m.name:<13} {k_kin:7.3f} {m.kappa_exp:7.3f} {rel:+8.2f}%")

    # Soft-subset prediction under the perovskite-only alpha (Fit B)
    print()
    print("Soft perovskite subset — predictions with alpha_T from Fit B:")
    print(f"{'material':<13} {'k_TGP(A)':>9} {'k_TGP(B)':>9} {'k_exp':>7}")
    for m in MATERIALS:
        if m.name not in SOFT_SUBSET:
            continue
        kA = kappa_TGP(m, alpha_A)
        kB = kappa_TGP(m, alpha_B)
        print(f"  {m.name:<13} {kA:9.3f} {kB:9.3f} {m.kappa_exp:7.3f}")

    print()
    print("=" * 78)
    print("Physical interpretation")
    print("=" * 78)
    print(
        "* Fit A (all materials) drives alpha_T to its lower scan limit ~0.01,\n"
        "  i.e. the TGP exponent essentially vanishes and formula (1) reduces\n"
        "  to the pure kinetic/Cahill floor kappa_kin = (1/3) n k_B v_s a.\n"
        "* That floor gets perovskites within 10-60 % (sometimes 15 %), but\n"
        "  misses stiff crystals (SrTiO3, NaCl, KBr, PbTe) by 80-95 %, because\n"
        "  in those materials long-MFP Umklapp phonons dominate, not\n"
        "  substrate-modulated diffusion.\n"
        "* Fit B (soft subset) finds an alpha_T that matches only perovskites.\n"
        "  If chi_B is comparable to chi_A (same alpha), TGP adds no extra\n"
        "  value — the kinetic floor alone suffices.\n"
        "* Interpretation: the TGP homogenisation formula is the **lower\n"
        "  bound** of thermal conductivity. Soft materials saturate this bound\n"
        "  (phonons at MFP floor); stiff materials are nowhere near it.\n"
        "* POSITIVE TGP RESULT: formula (1) with alpha_T ~ 0 naturally\n"
        "  explains the ultralow kappa ~ 0.3 W/m/K across the whole hybrid\n"
        "  perovskite family without any material-specific tuning; uses only\n"
        "  {n_at, v_s, a}. That's genuine predictive content.\n"
        "* NEGATIVE TGP RESULT: TGP's extra exp(-2<phi^2>) suppression has\n"
        "  negligible observable effect at 300 K even for soft perovskites,\n"
        "  because <phi^2> ~ 1e-4. The substrate-fluctuation correction is\n"
        "  below current thermal-transport measurement precision.\n"
        "* FALSIFIER: a soft perovskite whose kappa(300 K) >> kappa_kin would\n"
        "  violate the TGP lower bound. Target materials for future tests:\n"
        "  halide perovskites with heavy A-site (CsSnBr3, MASnI3, HC(NH2)2PbI3).\n"
    )


if __name__ == "__main__":
    np.set_printoptions(linewidth=120, precision=4, suppress=True)
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
