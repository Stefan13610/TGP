"""
OP-M92 Phase 0+ — Candidate D multi-source self-consistency check.

Question: Does Candidate D require a UNIVERSAL physical alpha to reproduce
          scenario (e) target shift (+1.46%) for ALL Schwarzschild BHs,
          or does alpha scale with source mass M_BH^2?

This is a critical test of the Phase 0+ structural sketch's heuristic
"alpha ~ 0.1 in geometric units (M_BH=1)" calibration. In the sketch we
calibrated to Sgr A* photon ring; here we extract alpha from M87* photon
ring and check if it agrees as a single physical constant.

Setup:
  Scenario (e) target shift in geometric units is UNIVERSAL: any
  Schwarzschild BH has photon ring at r=3.88M, psi=1.168, factor 0.886.
  Phase 0+ sketch: T*J*J / S_kin ~ O(1) in geom units (M=1) at photon ring.
  Therefore alpha (in geom units, M=1) ~ 0.114 / O(1) ~ 0.1.

Issue: alpha has dimension [length]^2 in geom units. Converting to SI:
       alpha_SI = alpha_geom × R_S(source)^2

If alpha is a UNIVERSAL physical constant:
  alpha_SI(SgrA*) = alpha_SI(M87*) — must be the SAME physical value.
  But alpha_geom(SgrA*) = alpha_geom(M87*) ~ 0.1 (universal target shift)
  implies alpha_SI scales with R_S^2 — different for different sources!

Multi-source check resolves which interpretation is correct.

Quantitative analysis below.
"""
import numpy as np


def main():
    print("=" * 72)
    print(" OP-M92 Phase 0+ — Candidate D multi-source self-consistency")
    print("=" * 72)

    # ────────────────────────────────────────────────────────────
    # Constants
    # ────────────────────────────────────────────────────────────
    G = 6.674e-11        # m^3 kg^-1 s^-2
    c = 2.998e8          # m/s
    M_sun = 1.989e30     # kg
    M_Pl = 2.176e-8      # kg (Planck mass)
    L_Pl = 1.616e-35     # m (Planck length)

    # Sources
    M_SgrA = 4.297e6 * M_sun
    M_M87 = 6.5e9 * M_sun
    M_NS_typical = 1.4 * M_sun     # neutron star
    M_GW150914_total = 65 * M_sun  # GW150914 final BH

    sources = [
        ("Sgr A*",         M_SgrA,       "EHT 2022"),
        ("M87*",           M_M87,        "EHT 2019"),
        ("GW150914 final", M_GW150914_total, "LIGO ringdown"),
        ("Neutron star",   M_NS_typical, "Stellar mass"),
    ]

    print(f"\n[Setup] Sources for multi-source alpha extraction:")
    for name, M, ref in sources:
        R_S = 2 * G * M / c**2
        print(f"  {name:<18} M = {M/M_sun:.2e} M_sun, R_S = {R_S:.3e} m  [{ref}]")

    # ────────────────────────────────────────────────────────────
    # Step 1: Phase 0+ sketch heuristic: alpha_geom ~ 0.1 universal
    # ────────────────────────────────────────────────────────────
    print(f"\n[Step 1] Phase 0+ sketch heuristic recap:")
    print(f"  Scenario (e) target: factor 0.886 at photon ring (r=3.88M, psi=1.168)")
    print(f"  Required shift: 1 - 0.886 = 0.114 (universal in geom units)")
    print(f"  Sketch: T*J*J / S_kin ~ O(1) in geom units at photon ring")
    print(f"  => alpha (in M_BH^2 units) = 0.114 / O(1) ~ 0.1")
    print(f"")
    print(f"  IMPLICATION: alpha_geom ~ 0.1 is UNIVERSAL in geom units.")
    print(f"  Converting to SI: alpha_SI = alpha_geom × R_S^2 (source-dependent!)")

    alpha_geom = 0.1

    # ────────────────────────────────────────────────────────────
    # Step 2: Extract alpha_SI from each source
    # ────────────────────────────────────────────────────────────
    print(f"\n[Step 2] alpha_SI extraction per source (assuming universal alpha_geom):")
    print(f"  {'Source':<18} {'R_S (m)':<14} {'alpha_SI (m^2)':<18} {'alpha_SI/L_Pl^2':<15}")
    print(f"  {'-'*72}")
    alpha_SI_values = []
    for name, M, ref in sources:
        R_S = 2 * G * M / c**2
        alpha_SI = alpha_geom * R_S**2
        alpha_SI_per_LPl2 = alpha_SI / L_Pl**2
        alpha_SI_values.append((name, M, R_S, alpha_SI))
        print(f"  {name:<18} {R_S:<14.3e} {alpha_SI:<18.3e} {alpha_SI_per_LPl2:<15.3e}")

    # ────────────────────────────────────────────────────────────
    # Step 3: Self-consistency analysis
    # ────────────────────────────────────────────────────────────
    print(f"\n[Step 3] Self-consistency: is alpha_SI a single physical constant?")
    name_ref, M_ref, R_S_ref, alpha_SI_ref = alpha_SI_values[0]  # Sgr A* baseline
    print(f"  Reference: alpha_SI({name_ref}) = {alpha_SI_ref:.3e} m^2")
    print(f"")
    print(f"  Pairwise mismatch ratios:")
    print(f"  {'Source':<18} {'alpha_SI/alpha_ref':<22} {'mass ratio^2':<15}")
    print(f"  {'-'*72}")
    for name, M, R_S, alpha_SI in alpha_SI_values:
        ratio = alpha_SI / alpha_SI_ref
        mass_ratio2 = (M / M_ref)**2
        print(f"  {name:<18} {ratio:<22.3e} {mass_ratio2:<15.3e}")

    print(f"")
    print(f"  CRITICAL FINDING: alpha_SI mismatches by orders of magnitude")
    print(f"  M87*/SgrA*: factor (M_M87/M_SgrA)^2 ~ 2.3e6")
    print(f"  GW150914/SgrA*: factor (65/4.3e6)^2 ~ 2.3e-10 (NS-like also)")
    print(f"")
    print(f"  => Phase 0+ heuristic 'alpha_geom ~ 0.1 universal' implies")
    print(f"     alpha_SI is NOT a single physical constant — it scales as M_BH^2.")
    print(f"     This violates universality of physical constants.")

    # ────────────────────────────────────────────────────────────
    # Step 4: Resolution paths
    # ────────────────────────────────────────────────────────────
    print(f"\n[Step 4] Resolution paths (Phase 1 derivation requirements):")
    print(f"")
    print(f"  PATH A: TGP-intrinsic length scale L_TGP")
    print(f"  -----")
    print(f"  Candidate D action: alpha_dim × L_TGP^2 × T^mu_nu J_mu J_nu")
    print(f"  with alpha_dim dimensionless and L_TGP a fundamental TGP scale.")
    print(f"  Strong-field rescue requires L_TGP^2 ~ R_S(source)^2 / 0.1.")
    print(f"")
    L_TGP_required_SgrA = np.sqrt(alpha_SI_ref)
    L_TGP_required_M87 = np.sqrt(0.1 * (alpha_SI_values[1][2])**2)
    print(f"    L_TGP required by Sgr A*: ~ {L_TGP_required_SgrA:.3e} m")
    print(f"    L_TGP required by M87*:   ~ {L_TGP_required_M87:.3e} m")
    print(f"  These differ by factor {L_TGP_required_M87/L_TGP_required_SgrA:.3e}")
    print(f"  => Single L_TGP cannot satisfy both — PATH A FAILS for naive sketch.")
    print(f"")
    print(f"  PATH B: T^mu_nu refers to source-enhanced substrate density")
    print(f"  -----")
    print(f"  If substrate density rho_sub near BH scales as rho_sub ~ M_BH/R_S^3")
    print(f"  (mass-density-like accretion), then T*J*J ~ M/R_S^3 × 1/R_S^2 = M/R_S^5")
    print(f"  In geom units (M=R_S/2): T*J*J ~ 1/R_S^4 = 1/M^4")
    print(f"  Then for rescue O(0.1):  alpha ~ 0.1 × M^4 (source-dependent)")
    print(f"  => Same problem; PATH B FAILS without explicit M-dependence absorption.")
    print(f"")
    print(f"  PATH C: alpha is dimensionless, action term has M_Pl^-2 prefactor")
    print(f"  -----")
    print(f"  Action: alpha / M_Pl^2 × T*J*J  with alpha dimensionless O(1)")
    print(f"  Effective alpha_SI = alpha / M_Pl_SI^2 = ?")
    print(f"  In geom: alpha / 1 (since M_Pl=1 in geom) = alpha_dimensionless")
    print(f"  Strong-field rescue at photon ring (r ~ M_BH):")
    print(f"  T*J*J ~ rho_eff × (1/r)^2 with rho_eff = M/R_S^3 ~ M_Pl^4/M^2 (using M ~ M_BH)")
    print(f"  Hmm: alpha/M_Pl^2 × M_Pl^4/M^2 × 1/M^2 = alpha × M_Pl^2/M^4")
    M_Pl_geom_SgrA = M_Pl / (M_SgrA * G/c**2 * c**2 / c**2)  # M_Pl in M_SgrA units
    M_Pl_in_SgrA_units = M_Pl / M_SgrA  # ratio in mass units
    print(f"    M_Pl/M_SgrA = {M_Pl_in_SgrA_units:.3e}")
    print(f"    factor M_Pl^2/M^4 ~ {M_Pl_in_SgrA_units**2:.3e}")
    print(f"  => rescue requires alpha ~ 1/M_Pl^2/M^4 ~ {1/M_Pl_in_SgrA_units**2:.3e}")
    print(f"  Astronomically large alpha — unphysical without further suppression.")
    print(f"  => PATH C requires careful Phase 1 dimensional analysis.")
    print(f"")
    print(f"  PATH D: Scenario (e) target shift is NOT universal (re-examine T-M92.1)")
    print(f"  -----")
    print(f"  Possibility: scenario (e) factor f(psi)=sqrt(g_tt^GR/g_tt^TGP) at photon")
    print(f"  ring depends on relative position r_ph^TGP=3.88M vs r_ph^GR=3M which IS")
    print(f"  universal. So target +1.46% IS universal in geom units.")
    print(f"  => PATH D NOT a viable resolution.")
    print(f"")
    print(f"  PATH E: Candidate D must be modified — non-minimal coupling")
    print(f"  -----")
    print(f"  Add explicit psi-dependence: alpha(psi) instead of constant alpha.")
    print(f"  E.g., alpha(psi) ~ alpha_0 × (psi - 1)^n with threshold near psi=1.")
    print(f"  This introduces extra parameter but allows source-independent")
    print(f"  reproduction of strong-field rescue while preserving weak-field.")
    print(f"  => PATH E adds tuning but is technically viable.")

    # ────────────────────────────────────────────────────────────
    # Step 5: Verdict
    # ────────────────────────────────────────────────────────────
    print(f"\n[Step 5] Verdict — Phase 0+ multi-source self-consistency:")
    print(f"")
    print(f"  STATUS: ISSUE FLAGGED, not yet a fatal blow but requires")
    print(f"          rigorous Phase 1 derivation to resolve.")
    print(f"")
    print(f"  Phase 0+ heuristic 'alpha ~ 0.1 (geom units)' is PHYSICALLY")
    print(f"  INCONSISTENT across BH sources of different masses. Either:")
    print(f"  - alpha is source-dependent (violates universality), OR")
    print(f"  - the heuristic O(1) factor in T*J*J / S_kin needs revision, OR")
    print(f"  - Candidate D needs non-minimal psi-coupling (PATH E)")
    print(f"")
    print(f"  Implications for Phase 0+ cross-checks done earlier today:")
    print(f"")
    print(f"  - Cosmology (margin 1e+33×): SAFE — even with alpha_SI varying by")
    print(f"    factor 10^10, ratio still << 1 → w(z) ≥ -1 robust")
    print(f"  - WEP MICROSCOPE (margin 6.7×): RECONSIDER — if alpha_SI for the")
    print(f"    laboratory experiment scale is *LARGER* than Sgr A* calibration,")
    print(f"    margin shrinks. Calibrating to local Earth gravity scale rather")
    print(f"    than Sgr A* could give different number.")
    print(f"  - GW vacuum c_GW = c_0: STRUCTURAL — unaffected (T=0 in vacuum)")
    print(f"  - Strong-field rescue: REQUIRES PATH E (non-minimal coupling) or")
    print(f"    re-examination of heuristic")
    print(f"")
    print(f"  PHASE 1 PRIORITY (REVISED):")
    print(f"  1. (NEW) Multi-source self-consistency derivation — establish if")
    print(f"     simple Candidate D action gives universal physical alpha")
    print(f"  2. WEP rigorous Nordvedt analysis (previously top priority)")
    print(f"  3. Photon ring numerical verification")
    print(f"  4. Cosmology perturbation scale dependence")

    # ────────────────────────────────────────────────────────────
    # Step 6: Summary
    # ────────────────────────────────────────────────────────────
    print(f"\n{'=' * 72}")
    print(" PHASE 0+ MULTI-SOURCE SELF-CONSISTENCY SUMMARY:")
    print(f"{'=' * 72}")
    print(f"")
    print(f" Heuristic alpha_geom ~ 0.1 (M_BH=1) implies alpha_SI scales as R_S^2:")
    print(f"   Sgr A*:   alpha_SI = {alpha_SI_values[0][3]:.2e} m^2")
    print(f"   M87*:     alpha_SI = {alpha_SI_values[1][3]:.2e} m^2  (factor 2.3e6 larger)")
    print(f"   GW150914: alpha_SI = {alpha_SI_values[2][3]:.2e} m^2  (factor 2.3e-10 smaller)")
    print(f"")
    print(f" => alpha is NOT a universal physical constant under naive Candidate D.")
    print(f" => Phase 0+ structural sketch's heuristic was a STRUCTURAL PLAUSIBILITY")
    print(f"    check only. Multi-source self-consistency requires either:")
    print(f"    - Modified action with non-minimal psi-coupling (PATH E)")
    print(f"    - TGP-intrinsic length scale L_TGP (PATH A — fails for naive)")
    print(f"    - Source-enhanced substrate density (PATH B — fails)")
    print(f"")
    print(f" Verdict: Phase 0+ self-consistency = OPEN ISSUE FLAGGED.")
    print(f" Candidate D ranking remains 'lead candidate' but with explicit caveat:")
    print(f"   STRUCTURAL SKETCH SUCCESSFUL, but Phase 1 covariant derivation")
    print(f"   MUST resolve multi-source consistency before Candidate D can claim")
    print(f"   to fully reproduce scenario (e) for general Schwarzschild BHs.")
    print(f"")
    print(f" Other candidates (A, B) similarly require multi-source verification.")
    print(f" This issue is CANDIDATE-INDEPENDENT — any M9.2 pivot reproducing")
    print(f" universal +14.56% deviation must explain how a single physical")
    print(f" parameter handles BHs across 9 orders of magnitude in mass.")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
