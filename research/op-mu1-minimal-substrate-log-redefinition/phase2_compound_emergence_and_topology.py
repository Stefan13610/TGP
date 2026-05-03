#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
μ.1 Phase 2 — Compound emergence test + topological Σε=2 constraint

P2.1: Multi-soliton ψ-superposition (analytic).
      Pokazać że ψ_total = Σψ_i ⟹ g_total = exp(Σε_i) w linearization.

P2.2: Topological constraint dla Σε = 2. Cztery kandydaci:
      (a) Electric charge: charged-lepton ±1, factor 2 z particle+antiparticle
      (b) Spinor double-cover: 4π identity dla spin-1/2, ratio 2
      (c) R3 winding number n=2 z homotopy structure
      (d) Compound saturation: substrate natural cutoff

P2.3: Verify Σε_topology = 2 ⟹ g_bg = e² ⟹ X = e²/2 dla α=2.

KEY: Compound emergence requires Σε = 2 z first principles, NIE empirical fit.
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


ALPHA = 2.0
G0_TAU = 1.77472
PSI0_TAU = np.log(G0_TAU)


def psi_ode(r, y, alpha):
    psi, psip = y
    rhs = (1 - np.exp(psi)) * np.exp((1 - 2 * alpha) * psi) \
          - (1 + alpha) * psip**2 - (2 / r) * psip
    return [psip, rhs]


def solve_psi_soliton(psi0, alpha=ALPHA, r_max=50.0, n_points=2000):
    sol = solve_ivp(psi_ode, (1e-3, r_max), [psi0, 0.0], args=(alpha,),
                    method='RK45', rtol=1e-10, atol=1e-12,
                    dense_output=True, max_step=0.05)
    r = np.linspace(1e-3, r_max, n_points)
    psi = sol.sol(r)[0]
    return r, psi


def p2_1_multi_soliton_analytic():
    print("=" * 78)
    print("  P2.1 — Multi-soliton ψ-superposition (analytic)")
    print("=" * 78)
    print()

    print("Założenie podstawowe (μ.1 hipoteza):")
    print("  ψ jest pierwotną zmienną substratu, więc multi-soliton background")
    print("  spełnia LINEAR superposition w ψ:")
    print()
    print("    ψ_total(x) = Σ_i ψ_i(x − x_i)        [aksjomat]")
    print()
    print("Konsekwencja dla g (geometric metric):")
    print("    g_total(x) = exp(ψ_total(x)) = exp(Σ_i ψ_i)")
    print("              = ∏_i exp(ψ_i(x))")
    print("              = ∏_i g_i(x)              [multiplicative composition]")
    print()
    print("LINEARYZACJA dla małych perturbacji od vacuum (g_i ≈ 1, ψ_i ≈ 0):")
    print("    g_i = 1 + ε_i  ⟹  ψ_i = log(1 + ε_i) ≈ ε_i − ε_i²/2 + O(ε³)")
    print()
    print("    ψ_total ≈ Σ_i ε_i − (1/2)Σ_i ε_i² + ...")
    print()
    print("    g_total = exp(Σψ_i) ≈ exp(Σε_i)·exp(−Σε_i²/2)")
    print()
    print("    Dla N małych równych perturbacji (ε_i = ε << 1):")
    print("    g_total ≈ exp(N·ε)·exp(−N·ε²/2)")
    print()
    print("    W limit N → ∞, ε → 0, N·ε = const = S:")
    print("    g_total → exp(S)·exp(−ε·S/2)·... → exp(S) [dominant]")
    print()
    print(f"  ⇒ Φ_total → exp(Σε)  (compound formula AUTOMATIC pod ψ-redef.)")
    print()

    # Sympy verification of expansion
    eps, N = sp.symbols('epsilon N', positive=True)
    psi_one = sp.log(1 + eps)
    psi_total_n = N * psi_one  # N solitonów takiej samej amplitudzie
    g_total = sp.exp(psi_total_n)
    g_total_expanded = sp.series(g_total, eps, 0, 4).removeO()
    print(f"  Sympy expansion (N solitonów, ε small):")
    print(f"    g_total ≈ {sp.expand(g_total_expanded)}")
    print()

    # In limit N·ε = S fixed, ε → 0
    print("  W limit N·ε = S fixed, ε → 0:")
    print("    g_total → exp(S) precisely")
    print()
    print("  ✓ P2.1 ANALYTIC PASS: compound emergence trywialna w ψ-substrate")
    print()

    return True


def candidate_a_electric_charge():
    """(a) Electric charge: charged-lepton ±1, factor 2 from particle+antiparticle."""
    print("--- Kandydat (a): Electric charge — particle + antiparticle ---")
    print()
    print("  Argument:")
    print("    Charged lepton (e⁻, μ⁻, τ⁻) ma ładunek q=−1 (elementarny).")
    print("    Vacuum w QFT zawiera quantum fluctuations e⁺e⁻ pairs.")
    print("    ε kontrybucja per particle: 1 (dla particle) + 1 (dla antiparticle)")
    print("    bo |q|² = 1 dla obu, integrating |amplitude|² over both contributions:")
    print()
    print("    Σε = |q_particle|² + |q_antiparticle|² = 1 + 1 = 2")
    print()
    print("  PROBLEMY:")
    print("    1. Czemu QED-like Σ|q|² zamiast Σq (gdzie ładunki by się znosiły)?")
    print("    2. Dla photon-mediated effects, normalna kontrybucja to fine-structure")
    print("       α ≈ 1/137, nie unit Σ.")
    print("    3. Argument wymaga że TGP-substrate liczy charged-lepton self-energy")
    print("       jako particle + antiparticle pair — ale to NIE jest soliton-like,")
    print("       to QFT fluctuation.")
    print()
    print("  SCORE: 2/10 — argument ad hoc, nie z first principles")
    print()
    return False, "Σε=2 z charge particle+antiparticle, ale nieprzekonujący"


def candidate_b_spinor_double_cover():
    """(b) Spinor double-cover: rotation 4π identity."""
    print("--- Kandydat (b): Spinor double-cover — SU(2) → SO(3) ---")
    print()
    print("  Argument:")
    print("    Grupa SO(3) (rotations 3D) ma double-cover SU(2). Spin-1/2 stany")
    print("    transformują się pod SU(2): rotation 2π zwraca −ψ_state, rotation 4π")
    print("    zwraca +ψ_state (identity).")
    print()
    print("    Topologicznie: π₁(SO(3)) = ℤ/2ℤ, więc closed loop w SO(3) ma")
    print("    homotopy class 0 (trivial) lub 1 (nontrivial).")
    print()
    print("    DLA charged-lepton (spin-1/2): substrate background musi 'closure'")
    print("    pod spin-rotation. Background z homotopy class 1 = single windings,")
    print("    z class 0 = double windings. Background ψ closure wymaga 4π:")
    print()
    print("    Σε_substrate = 4π·(1/2π) = 2  [normalizacja przez 2π = unit angle]")
    print()
    print("  PROBLEMY:")
    print("    1. Substrate ψ jest skalarne — czy 'spinor closure' applies?")
    print("    2. Jeśli ψ_background nie nosi spin, to argument ψ ↔ spinor")
    print("       wymaga dodatkowego mappingu.")
    print("    3. Result Σε = 2 dependent on normalizacji jednostek.")
    print()
    print("  SCORE: 4/10 — geometryczny ale skalar/spinor mapping wymagany")
    print()
    return False, "Σε=2 z spinor 4π, ale wymaga mapping skalar→spinor"


def candidate_c_winding_number():
    """(c) R3 winding number z homotopy."""
    print("--- Kandydat (c): R3 winding number n=2 z homotopy ---")
    print()
    print("  Argument:")
    print("    R3 substrate dla α=2 (3D space), soliton jest mapem")
    print("    R³ → ℝ⁺ (dla g) ≅ ℝ (dla ψ). Homotopy:")
    print("      π₃(ℝ⁺) = trivial (kontraktibilny)")
    print("      π₃(target with internal symmetry) = ℤ jeśli target jest S²/S³")
    print()
    print("    Dla TGP charged-lepton z internal U(1) charge: target = S¹ × ℝ⁺,")
    print("    i nontrivial homotopy w π₃(S¹) = 0 — banały, brak topology.")
    print()
    print("    ALE: jeśli rozważyć GROUP STRUCTURE substratu jako ψ ∈ Lie(ℝ⁺),")
    print("    multiplicative composition g₁·g₂ jest natural, i Σε = 2 może")
    print("    być wymuszone przez topology of the bundle R³ → R³/(±1) double-cover.")
    print()
    print("  PROBLEMY:")
    print("    1. R3 substrate w α=2 nie ma natural topological charge (ψ ∈ ℝ).")
    print("    2. Argument requires TGP-specific bundle structure która nie")
    print("       jest jeszcze udokumentowana w obecnym substrate definition.")
    print("    3. Brak konkretnego mechanizm produkujący n=2 z geometrii.")
    print()
    print("  SCORE: 2/10 — wymagałby nowej topology TGP-substratu")
    print()
    return False, "Σε=2 z winding n=2, brak natural topology w ψ ∈ ℝ"


def candidate_d_compound_saturation():
    """(d) Compound saturation: substrate natural cutoff Σε = 2."""
    print("--- Kandydat (d): Compound saturation — natural substrate cutoff ---")
    print()
    print("  Argument:")
    print("    Compound interest formula: (1 + S/N)^N → exp(S) jako N → ∞.")
    print("    Substrate Φ₀ jako produkt nieskończony: Φ₀ = ∏(1 + ε_i)")
    print("    saturates dla różnych S = Σε:")
    print()
    print("      S = 1 ⟹ Φ₀ = e ≈ 2.718")
    print("      S = 2 ⟹ Φ₀ = e² ≈ 7.389")
    print("      S = 3 ⟹ Φ₀ = e³ ≈ 20.09")
    print()
    print("    NATURALNY argument dla S = 2:")
    print("      Substrate ma DWA kanały (amplitude + phase per Section 11.5 audytu).")
    print("      Każdy kanał wnosi unit ε per soliton (background contribution = 1).")
    print("      Total: Σε = 1 (amplitude) + 1 (phase) = 2.")
    print()
    print("    LUB: dla charged-lepton (massive, charged), background musi sumować")
    print("    DWA charges (lepton + reflected antiparticle from vacuum):")
    print("      Σε = ε_self + ε_image = 1 + 1 = 2")
    print()
    print("  PROBLEMY:")
    print("    1. Argument z 'two channels' jest plausible ale wymaga że amplitude")
    print("       i phase contribute równo (po 1 each) — to jest empirical assumption.")
    print("    2. 'Image charge' argument wymaga that vacuum reflects soliton —")
    print("       nie ma direct source w R3 ODE.")
    print("    3. NIE wybiera unikalnie S=2 vs S=1 — czemu nie e zamiast e²?")
    print()
    print("  SCORE: 5/10 — najlepszy z 4 ale nadal motivated, nie derived")
    print()
    return False, "Σε=2 z natural saturation 2 channels, plausible ale empirical"


def p2_2_topology_test():
    print("=" * 78)
    print("  P2.2 — Topological Σε = 2 constraint test (4 kandydaci)")
    print("=" * 78)
    print()
    print("  GATE: AT LEAST 1/4 musi PASS jako 'first-principles' argument.")
    print("        Score >= 7/10 = PASS, < 7/10 = FAIL.")
    print()

    results = []
    pass_a, expl_a = candidate_a_electric_charge()
    results.append(("(a) electric charge", pass_a, 2, expl_a))

    pass_b, expl_b = candidate_b_spinor_double_cover()
    results.append(("(b) spinor double-cover", pass_b, 4, expl_b))

    pass_c, expl_c = candidate_c_winding_number()
    results.append(("(c) winding n=2", pass_c, 2, expl_c))

    pass_d, expl_d = candidate_d_compound_saturation()
    results.append(("(d) natural saturation", pass_d, 5, expl_d))

    print("=" * 78)
    print("  P2.2 summary: topology candidates dla Σε = 2")
    print("=" * 78)
    print()
    print(f"  {'Kandydat':<30} {'Score':>8} {'Pass':>6}")
    print("  " + "-" * 50)
    for name, p, score, _ in results:
        passmark = "PASS" if score >= 7 else "FAIL"
        print(f"  {name:<30} {score}/10{'':<5} {passmark:>6}")
    print()

    n_pass = sum(1 for _, _, score, _ in results if score >= 7)
    print(f"  P2.2 GATE: {n_pass}/4 kandydatów PASS (wymagane >= 1)")
    print()

    if n_pass >= 1:
        print("  ✓ P2.2 PASS")
    else:
        print("  ✗ P2.2 FAIL — żaden kandydat nie daje Σε=2 z first principles")
        print()
        print("  IMPLIKACJA: μ.1 redukuje się do PURE REPARAMETRIZATION.")
        print("  ψ-substrate jest matematycznie identyczny z g-substrate")
        print("  (P1 udowodnił to do floating-point precision), ALE:")
        print()
        print("  - compound emergence formy g_total = exp(Σε) JEST present strukturalnie")
        print("  - VALUE of Σε dla charged-lepton α=2 NIE jest derived z first principles")
        print("  - dla matching mass formula X = e²/2, musimy POSTULATE Σε = 2 — ad hoc")
        print()
        print("  ⇒ μ.1 nie znajduje fizycznego mechanizmu, tylko zmienia język.")

    return n_pass >= 1


def p2_3_compound_value_check():
    """Verify if Σε = 2 (postulated, ad hoc) gives correct e²/2 in mass formula."""
    print("=" * 78)
    print("  P2.3 — Compound value verification (assuming Σε=2 postulated)")
    print("=" * 78)
    print()

    # Jeśli Σε = 2 (postulated), to g_bg = exp(2) = e² = 7.389
    # I X = e²/2 dla α=2 charged-lepton wynika z structural relation
    # m = c·A²·g₀^X = c·A²·exp(X·ψ₀); dla "background-saturated" lepton
    # ψ₀ ≈ Σε/2 (per soliton contribution):

    # Hmm — to nie jest jasne. Sprawdźmy jaki dokładnie jest direct link
    # między Σε = 2 background a X = e²/2 mass formula slope.

    print("  Direct relationship Σε ↔ X niejasny:")
    print()
    print("    Background g_bg = exp(Σε) opisuje multi-soliton vacuum.")
    print("    Mass formula slope X = e²/2 = 3.6945 jest pochodną")
    print("    log(m/A²) względem ψ₀_lepton.")
    print()
    print("    Pytanie: czy slope X = exp(Σε)/2 (jeśli Σε = 2)?")
    print(f"    Test: exp(2)/2 = {np.exp(2)/2:.6f}")
    print(f"          X observed = {np.exp(2)/2:.6f}")
    print(f"          drift: 0%")
    print()
    print("  TO JEST KONSYSTENTNE, ALE TAUTOLOGICZNE:")
    print("    X = e²/2 z λ.1 mass formula (empirical fit)")
    print("    Σε = 2 postulated jako compound saturation")
    print("    → 'derivation' X = exp(2)/2 jest po prostu definicja")
    print()
    print("  KLUCZOWE: bez NIEZALEŻNEGO derivation Σε = 2 (P2.2 FAIL),")
    print("  P2.3 jest cyrkularny.")
    print()
    print("  ✗ P2.3 NEUTRAL: structural form match, ale brak independence.")

    return False  # neutral, not pass


def main():
    print("=" * 78)
    print("  μ.1 Phase 2 — Compound emergence + topological Σε=2")
    print("=" * 78)
    print()

    p2_1 = p2_1_multi_soliton_analytic()
    print()
    p2_2 = p2_2_topology_test()
    print()
    p2_3 = p2_3_compound_value_check()
    print()

    print("=" * 78)
    print("  Phase 2 verdict")
    print("=" * 78)
    print()
    print(f"  P2.1 (multi-soliton ψ analytic):       {'PASS' if p2_1 else 'FAIL'}")
    print(f"  P2.2 (topology Σε=2 candidates):       {'PASS' if p2_2 else 'FAIL'}")
    print(f"  P2.3 (compound→X structural match):    {'NEUTRAL' if not p2_3 else 'PASS'}")
    print()

    if p2_1 and p2_2:
        print("  → Phase 2 PASS — compound mechanism identified z first-principles topology")
    else:
        print("  → Phase 2 PARTIAL/FAIL")
        print()
        print("  KEY FINDING: P1 udowodnił że ψ ≡ log g jest mathematically")
        print("  invariant (drift 1e-13%), ale P2.2 nie znajduje natural source")
        print("  dla Σε = 2 z TGP-substrate first principles.")
        print()
        print("  Best candidate: (d) compound saturation z 2 channels (5/10).")
        print("  Wymaga: explicit definition kanałów amplitude/phase substrate.")
        print()
        print("  μ.1 verdict: REPARAMETRIZATION SUKCES, MECHANISM PARTIAL.")
        print("  Pole g_total = exp(Σε) emerges automatycznie z ψ-substrate,")
        print("  ale value Σε = 2 dla α=2 charged-lepton nie ma derivation.")
    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
