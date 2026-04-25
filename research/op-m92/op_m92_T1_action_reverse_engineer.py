"""
OP-M92 T-M92.1 — reverse-engineer action principle requirements
                  for scenario (e) coupling.

Scenario (e) coupling (from OP-EHT T3, fitted ansatz):
    f(psi) = sqrt(g_tt^GR / g_tt^TGP)
           = sqrt( (1 - 2M/r) * psi / (4 - 3*psi) )

where g_tt^GR = -(1 - 2M/r) and g_tt^TGP = -(4 - 3 psi)/psi.

The non-locality challenge: f(psi) at the photon ring depends on BOTH:
  - the local TGP psi value (=1.168 at photon ring),
  - the GR Schwarzschild backbone evaluated at the photon ring radius
    (r = 3.88 M for TGP photon ring, gives g_tt^GR = -0.485).

This script analyzes:
  1. Power-law approximation: what alpha makes f(psi) = psi^(-alpha)
     match scenario (e) at photon ring?
  2. Polynomial expansions: minimal coupling form preserving f(1)=1.
  3. Action principle requirements: what L_mat[Phi, psi] reproduces
     the scenario (e) coupling structurally?

Output: structural feasibility map for M9.2 candidates.
"""
import sympy as sp


def main():
    print("=" * 70)
    print(" OP-M92 T-M92.1 — Reverse-engineer scenario (e) action requirements")
    print("=" * 70)

    psi, M_bh, r, U = sp.symbols('psi M r U', positive=True, real=True)
    alpha, beta_p = sp.symbols('alpha beta_p', real=True)

    # ────────────────────────────────────────────────────────
    # Step 1: Scenario (e) coupling at photon ring
    # ────────────────────────────────────────────────────────
    print(f"\n[T-M92.1.1] Scenario (e) coupling at photon ring:")
    psi_ph_val = sp.Rational(1168, 1000)    # M9.1'' photon ring psi (geometric invariant)
    r_ph_TGP = sp.Rational(388, 100)        # M9.1'' photon ring radius (in M)
    r_ph_GR = sp.Rational(3, 1)             # GR Schwarzschild photon ring radius
    M_val = 1

    g_tt_GR = -(1 - 2*M_bh/r)
    g_tt_TGP = -(4 - 3*psi)/psi

    f_e = sp.sqrt(-g_tt_GR / -g_tt_TGP)  # both negative → positive ratio
    f_e_simpl = sp.simplify(f_e)
    print(f"  f_e(psi, r, M) = {f_e_simpl}")

    print(f"\n  Two interpretations of where to evaluate g_tt^GR:")
    f_at_ph_TGP = f_e.subs([(psi, psi_ph_val), (r, r_ph_TGP), (M_bh, M_val)])
    f_at_ph_TGP_val = float(sp.simplify(f_at_ph_TGP))
    print(f"  (i) at TGP photon ring (r=3.88 M):    f_e = {f_at_ph_TGP_val:.4f}")

    f_at_ph_GR = f_e.subs([(psi, psi_ph_val), (r, r_ph_GR), (M_bh, M_val)])
    f_at_ph_GR_val = float(sp.simplify(f_at_ph_GR))
    print(f"  (ii) at GR photon ring (r=3 M):       f_e = {f_at_ph_GR_val:.4f}")
    print(f"  Target from OP-EHT T3 scenario (e):   0.886")
    print(f"  Matches target (interpretation ii): {abs(f_at_ph_GR_val - 0.886) < 0.005}")
    print(f"  ⇒ Scenario (e) ansatz is doubly non-local: compares")
    print(f"    TGP at r_ph^TGP=3.88M to GR at r_ph^GR=3M (DIFFERENT r values)")
    f_at_ph_val = f_at_ph_GR_val

    # ────────────────────────────────────────────────────────
    # Step 2: Power-law approximation f(psi) = psi^(-alpha)
    # ────────────────────────────────────────────────────────
    print(f"\n[T-M92.1.2] Power-law approximation f(psi) = psi^(-alpha):")
    # Auto: f(1) = 1
    # At psi=1.168: 1.168^(-alpha) = f_at_ph_val (using interpretation ii)
    import math
    alpha_fit_val = -math.log(f_at_ph_val) / math.log(float(psi_ph_val))
    print(f"  alpha (fitted at psi_ph, target {f_at_ph_val:.4f}) = {alpha_fit_val:.4f}")
    print(f"  Closest rational: 3/4 = {3/4:.4f}")
    print(f"  Match within 5%: {abs(alpha_fit_val - 0.75) / 0.75 < 0.05}")

    # ────────────────────────────────────────────────────────
    # Step 3: 1PN expansion (psi = 1 + 2U + O(U^2))
    # ────────────────────────────────────────────────────────
    print(f"\n[T-M92.1.3] 1PN expansion of scenario (e) coupling:")
    # We need to expand both g_tt^GR and g_tt^TGP at 1PN.
    # GR: at r ~ M/U → g_tt^GR = -(1 - 2U)
    # TGP: psi = 1 + 2U → g_tt^TGP = -(4 - 3(1+2U))/(1+2U) = -(1 - 6U)/(1+2U)
    #                              = -(1 - 6U)(1 - 2U + 4U^2 - ...) ≈ -(1 - 8U + 16U^2)
    #
    # f²(psi, r) = (1-2U) / [ (1 - 6U)/(1+2U) ]
    #           = (1-2U)(1+2U)/(1-6U) = (1-4U²)/(1-6U)
    # Expand: ≈ (1 - 4U²)(1 + 6U + 36U² + ...) ≈ 1 + 6U + 32U² + O(U^3)
    f2_e_1pn = (1 - 2*U) / ((1 - 6*U)/(1 + 2*U))
    f2_e_1pn_simpl = sp.series(f2_e_1pn, U, 0, 4).removeO()
    print(f"  f^2_e(U) = {sp.simplify(f2_e_1pn_simpl)}")

    f_e_1pn = sp.sqrt(f2_e_1pn)
    f_e_1pn_series = sp.series(f_e_1pn, U, 0, 3).removeO()
    print(f"  f_e(U) ≈ {sp.simplify(f_e_1pn_series)}")
    # If f(U) ≈ 1 + 3U + ... then we have NEGATIVE 1PN deviation in matter coupling

    # Critical check: f_e(U=0) = 1? Yes, since at U=0 both metrics → -1
    f_at_zero = f_e_1pn.subs(U, 0)
    print(f"  f_e(U=0) = {f_at_zero} (must be 1 for 1PN preservation)")

    # ────────────────────────────────────────────────────────
    # Step 4: Action principle minimum requirements
    # ────────────────────────────────────────────────────────
    print(f"\n[T-M92.1.4] Action principle minimum requirements:")
    print(f"  Required structure of L_mat for scenario (e):")
    print(f"    L_mat^(e) = -q × phi × rho × f_e(psi, r) / Phi_0")
    print(f"  ")
    print(f"  Decomposition needed:")
    print(f"    f_e(psi, r) = sqrt(g_tt^GR(r) / g_tt^TGP(psi))")
    print(f"")
    print(f"  Problem: g_tt^GR(r) = -(1 - 2M/r) = -(1 - 2 phi_GR(r))")
    print(f"           where phi_GR is GR Newtonian potential.")
    print(f"  In M9.1'', phi_TGP(r) is the TGP scalar field;")
    print(f"  there is NO independent phi_GR object in single-field substrate.")
    print(f"")
    print(f"  ⇒ M9.2 must introduce a SECOND structure that mimics phi_GR(r)")
    print(f"    in the strong-field while reducing to phi_TGP(r) in weak-field.")
    print(f"")
    print(f"  Three structural options:")
    print(f"    (a) Auxiliary field sigma(x) with non-minimal coupling")
    print(f"    (b) Conformal frame g_mat = Omega^2(psi) g")
    print(f"    (c) ψ-dependent matter charge q(psi) with q(1)=q_0")
    print(f"  Each option needs T-M92.2..T-M92.5 analysis")

    # ────────────────────────────────────────────────────────
    # Step 5: Equivalent power-law derivation as M9.2 baseline
    # ────────────────────────────────────────────────────────
    print(f"\n[T-M92.1.5] Power-law f(psi) = psi^(-3/4) as M9.2 candidate baseline:")
    f_pow = psi**(sp.Rational(-3, 4))
    f_pow_at_ph = f_pow.subs(psi, psi_ph_val)
    f_pow_at_ph_val = float(f_pow_at_ph)
    print(f"  f_pow(psi=1.168) = {f_pow_at_ph_val:.4f}")
    print(f"  vs scenario (e) target: 0.886")
    print(f"  Discrepancy: {abs(f_pow_at_ph_val - 0.886)*100:.1f}%")

    # 1PN of psi^(-3/4)
    f_pow_1pn = sp.series(f_pow.subs(psi, 1 + 2*U), U, 0, 4).removeO()
    print(f"\n  1PN expansion: f_pow(1+2U) = {sp.simplify(f_pow_1pn)}")
    print(f"  Linear-in-U coefficient: -3/2 (vs scenario (e): +3 from above)")
    print(f"  ⇒ Power-law has WRONG sign at 1PN — Mercury γ_PPN shift")

    # ────────────────────────────────────────────────────────
    # Step 6: Scoring summary
    # ────────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(" T-M92.1 SUMMARY — action principle requirements:")
    print(f"{'=' * 70}")
    print(f"")
    print(f" Findings:")
    print(f"  1. Scenario (e) coupling f(psi, r) = sqrt(g_tt^GR(r)/g_tt^TGP(psi))")
    print(f"     is INHERENTLY non-local (depends on r as well as psi).")
    print(f"  2. Single-field M9.1'' cannot reproduce non-locality structurally.")
    print(f"  3. Power-law approx psi^(-3/4) matches at photon ring but breaks")
    print(f"     1PN preservation (linear-U coefficient wrong sign).")
    print(f"  4. M9.2 must introduce a SECOND structure (sigma field, conformal")
    print(f"     factor, q(psi) flow, or back-reaction) to bridge non-locality.")
    print(f"")
    print(f" Verdict: DIAGNOSTIC — establishes that scenario (e) coupling")
    print(f"          requires non-trivial structural extension of M9.1''.")
    print(f"          This justifies the 4-candidate M9.2 program in OP_M92_setup.md.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
