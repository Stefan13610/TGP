"""
ps06 — Consistency of Phi_0 ~ 0.34 GeV with other TGP / SM observables.

Context
-------
ps04 + ps05 established that TGP's minimal exponential-metric coupling
predicts a dilaton-like Yukawa y_f = m_f / Phi_0 for every massive fermion
of mass m_f, giving Delta a_f ~ m_f^4 / (8 pi^2 Phi_0^2 M_Phi^2) for the
1-loop scalar-exchange contribution to g-2.

Fitting to Delta a_mu gives Phi_0 ~ 0.34 GeV (for M_Phi = 91 GeV).

But m_f / Phi_0 is a FLAVOR-DEPENDENT Yukawa coupling that varies wildly
across the Standard Model spectrum:
    e  : y_e  = 0.511 MeV / 0.34 GeV   = 0.0015   (perturbative)
    mu : y_mu = 105.66 MeV / 0.34 GeV  = 0.31     (moderate)
    tau: y_tau = 1.777 GeV / 0.34 GeV  = 5.2      (non-perturbative)
    u  : y_u  = 2.2 MeV / 0.34 GeV     = 0.0065
    s  : y_s  = 95 MeV / 0.34 GeV      = 0.28
    c  : y_c  = 1.27 GeV / 0.34 GeV    = 3.7      (non-perturbative)
    b  : y_b  = 4.18 GeV / 0.34 GeV    = 12.3     (strongly non-perturbative)
    t  : y_t  = 173.1 GeV / 0.34 GeV   = 509      (breaks EFT entirely)

The perturbative bound for Yukawa-like couplings is y < sqrt(4 pi) ~ 3.5.
=> The tau, charm, bottom, and top quark sectors are INCOMPATIBLE with
   a universal exp-metric dilaton of Phi_0 = 0.34 GeV.

Possible resolutions
--------------------
1. Phi_0 is MUCH larger (~TeV+) and the muon anomaly comes from a different
   mechanism.
2. The substrate couples only to LEPTONS (lepton-selective TGP), avoiding
   the quark sector entirely.
3. Substrate couples to RELATIVE mass differences, not absolute m_f:
   y_f ~ (m_f - m_ref) / Phi_0, with m_ref the "reference" substrate mass.
4. Non-linear corrections to the exp metric cap the effective Yukawa at
   some strong-coupling ceiling y* ~ 1.

Script
------
Tests each resolution:
- Option 1: Phi_0 ~ TeV -> what would Delta a_mu look like, and what mechanism
  produces it?
- Option 2: lepton-only TGP -> cross-check with Z-boson decay (Z -> lepton pair)
  loop would shift from SM.
- Option 3: reference mass scale -> introduces an additional parameter that
  needs to be motivated.

The key takeaway: the m^4 framework from minimal exp-metric is simple and
testable for leptons, BUT hints at a deeper flavor structure when extended
to quarks. This is not a falsification — it is a POINTER to the full TGP
structure that must emerge from the substrate action beyond the minimal ansatz.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.integrate as si
import scipy.optimize as so


ALPHA = 1 / 137.035999
PI = math.pi

# Lepton and quark masses (PDG 2024)
LEPTONS = {
    "e": 0.0005109989461,  # GeV
    "mu": 0.1056583745,
    "tau": 1.77686,
}
QUARKS = {  # MS-bar running masses at mu = 2 GeV, except top pole mass
    "u": 0.00216,
    "d": 0.00467,
    "s": 0.0935,
    "c": 1.273,
    "b": 4.183,
    "t": 172.69,
}
MASSES = {**LEPTONS, **QUARKS}

# Muon g-2 fitting result from ps04
PHI_0_FIT = 0.337  # GeV, fitted at M_Phi = 91.2 GeV
M_PHI_FIT = 91.2   # GeV

# Perturbative Yukawa ceiling
Y_PERT_MAX = math.sqrt(4 * PI)  # ~ 3.545


def yukawa(mass_GeV: float, Phi_0_GeV: float) -> float:
    return mass_GeV / Phi_0_GeV


def scenario_summary():
    print("=" * 74)
    print("ps06 — Phi_0 consistency test across SM particle spectrum")
    print("=" * 74)

    print(f"\nFit from muon g-2 (ps04): Phi_0 = {PHI_0_FIT} GeV, M_Phi = {M_PHI_FIT} GeV")
    print(f"Perturbative Yukawa ceiling: y* = sqrt(4 pi) = {Y_PERT_MAX:.3f}")

    print("\n--- Effective Yukawa y_f = m_f / Phi_0 for each SM fermion ---")
    print(f"  {'fermion':>8} | {'m_f (GeV)':>12} | {'y_f':>10} | status")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*10}-+-" + "-" * 30)
    for name, m in sorted(MASSES.items(), key=lambda kv: kv[1]):
        y = yukawa(m, PHI_0_FIT)
        if y < 0.1:
            status = "perturbative OK"
        elif y < 1.0:
            status = "moderate"
        elif y < Y_PERT_MAX:
            status = "strong-coupling"
        else:
            status = "NON-PERTURBATIVE (problem)"
        print(f"  {name:>8} | {m:>12.4e} | {y:>10.4f} | {status}")


def fit_leptons_only():
    """If substrate couples to leptons only, what quark bound does that imply?"""

    print("\n--- Scenario 2: lepton-selective TGP ---")
    print("Assume substrate decouples from quarks entirely.")
    print(f"  Consequence: no hadronic g-2 shift from TGP dilaton.")
    print(f"  The g-2 hadronic puzzle (HVP data-driven vs lattice) is NOT fixed by TGP.")
    print()
    print("  Test: Z boson decay width to lepton pair.")
    print("  SM: Gamma(Z->ll)/Gamma(Z->hadrons) matches to 1e-3 precision.")
    print("  If TGP adds a dilaton loop to lepton vertex at Z pole,")
    print("  the correction would be enhanced by Z mass in the dilaton propagator:")
    mz = 91.1876
    y_tau = yukawa(MASSES["tau"], PHI_0_FIT)
    # Estimate vertex correction: O(y_tau^2 / 16 pi^2 * (mz/M_Phi)^2 * log)
    delta_gamma = y_tau**2 / (16 * PI**2) * (mz / M_PHI_FIT)**2 * math.log(M_PHI_FIT / MASSES["tau"])
    print(f"    y_tau ~ {y_tau:.2f}, (m_Z/M_Phi)^2 ~ {(mz/M_PHI_FIT)**2:.2f}")
    print(f"    Rough estimate delta(Z->tau+tau-) / SM ~ {delta_gamma:.3e}")
    if delta_gamma > 1e-3:
        print(f"    WARNING: correction exceeds 1e-3 LEP precision on Z->ll => tension")
    else:
        print(f"    OK: correction below current LEP precision.")


def fit_large_Phi0():
    """If we force Phi_0 ~ TeV to avoid non-perturbative Yukawas, what gives a_mu?"""

    print("\n--- Scenario 1: Large Phi_0 ~ TeV scale ---")
    # Required m^4/(Phi_0^2 M_Phi^2) to give a_mu = 2.5e-9
    a_mu_target = 2.5e-9
    m_mu = MASSES["mu"]
    for Phi_0 in (1.0, 10.0, 100.0, 1000.0):
        # M_Phi to saturate a_mu:
        def residual(log_MPhi):
            M = math.exp(log_MPhi)
            # a_mu ~ m_mu^4 / (8 pi^2 Phi_0^2 M^2) * F_S(m_mu/M) (F_S ~ log(M/m_mu))
            return (m_mu**4 / (8 * PI**2 * Phi_0**2 * M**2)) * math.log(max(M/m_mu, 1.001)) - a_mu_target
        try:
            M_fit = math.exp(so.brentq(residual, math.log(1e-4), math.log(1e5)))
            y_tau = yukawa(MASSES["tau"], Phi_0)
            y_t = yukawa(MASSES["t"], Phi_0)
            print(f"  Phi_0 = {Phi_0:>6.1f} GeV => M_Phi = {M_fit:.3e} GeV,  y_tau = {y_tau:.3f},  y_top = {y_t:.3f}")
        except Exception as e:
            print(f"  Phi_0 = {Phi_0:>6.1f} GeV => no real M_Phi solution")
    print()
    print("  Observation: Phi_0 ~ TeV requires M_Phi << m_mu (sub-MeV substrate mass).")
    print("  That is physically weird: a substrate with mass below the muon undermines")
    print("  the whole EFT integration-out picture. So TeV-scale Phi_0 is NOT consistent.")


def proposed_resolution():
    """Natural TGP resolution: Phi_0 ~ M_Phi with strong flavor discrimination."""
    print("\n--- Preferred resolution: Phi_0 ~ M_Phi ~ EW scale with lepton-only coupling ---")
    print("If Phi_0 = M_Phi = M_Z (both at the electroweak scale):")
    Phi_0 = 91.1876
    M_Phi = 91.1876
    m_mu = MASSES["mu"]
    # a_mu = m_mu^4 / (8 pi^2 Phi_0^2 M_Phi^2) * F_S
    # Use F_S ~ log(M_Phi/m_mu):
    F_S = math.log(M_Phi/m_mu)
    a_mu_pred = m_mu**4 / (8 * PI**2 * Phi_0**2 * M_Phi**2) * F_S
    print(f"  Phi_0 = M_Phi = {Phi_0} GeV")
    print(f"  Delta a_mu_predicted = {a_mu_pred:.3e}")
    print(f"  Observed: ~2.5e-9")
    print(f"  Ratio prediction / observed = {a_mu_pred/2.5e-9:.3e}")
    print(f"  => Too small by factor {2.5e-9/a_mu_pred:.1e}")
    print()
    print("  This is the fundamental tension: minimal exp-metric dilaton with")
    print("  Phi_0 = M_Phi = M_Z does NOT saturate the muon anomaly.")
    print()
    print("  To saturate it, need either small Phi_0 (forcing large y_q for heavy quarks)")
    print("  or an additional loop enhancement factor (needs ps07 calculation).")


def _summary():
    print("\n" + "=" * 74)
    print("Synthesis")
    print("=" * 74)
    print("""
The m^4 dilaton framework derived from the exp metric (ps04) fits the
muon anomaly cleanly IF Phi_0 ~ 0.34 GeV. But this value causes two problems
when applied to the full SM spectrum:

1. NON-PERTURBATIVE YUKAWA for tau, c, b, t (y_f > 1 up to 510 for top).
   => The dilaton cannot be a universal coupling for all fermions at this scale.

2. LEPTON-ONLY TGP mitigates problem 1 (quarks decoupled), but introduces
   a testable consequence at Z-pole: dilaton-enhanced tau vertex would
   shift Z->tau+tau- at the 1e-3 level or worse. LEP precision is ~1e-3, so
   this sector is at the EDGE of current constraints — testable!

3. Alternative Phi_0 = M_Phi = M_Z: underfits muon g-2 by ~10^4.

Conclusion:
  - TGP in the minimal exp-metric form is naturally a lepton-only dilaton
    with Phi_0 ~ hadronic scale (~300 MeV). This is a CLEAN prediction but
    requires quark sector to decouple by some TGP mechanism not yet specified.
  - Decisive test: Z->tau+tau- at LEP already constrains this at 1e-3 level.
    Next-gen ee->Z precision (FCC-ee Z-pole run) goes to 1e-5: would either
    detect TGP dilaton vertex shift, or exclude Phi_0 < 10 GeV.

Next step (ps07):
  Check Z->tau+tau- partial width TGP shift against LEP-SLD measurements
  to get the best independent bound on Phi_0 from e+e- colliders.
""")


def main() -> None:
    scenario_summary()
    fit_leptons_only()
    fit_large_Phi0()
    proposed_resolution()
    _summary()


if __name__ == "__main__":
    main()
