---
title: "Phase EXPLORATION — no-Yukawa hypothetical scope (γ-7)"
type: phase_exploration
status: EXPLORATION_ONLY
phase: EXPLORATION (post-HALT-B)
parent_cycle: op-CE-H-gamma-7-clumping-acceleration-2026-05-24
parent_close: Phase_FINAL_close.md
created_date: 2026-05-24
authorization: "User 2026-05-24: 'warto przeanalizować jaki byłby wynik bez yakuwa, i jaki wymagany rozmiar wszechśiata. Na zasadzie eksploracji nie stwierdzeń'"
halt_b_status: LOCKED_UNCHANGED
anti_lakatos_lock: PRESERVED
modifies_verdicts: NO
modifies_phase_locks: NO
purpose: "Gather scaling info BEFORE any future γ-8; test user hypothesis on mass-dependent screening"
substantive_fp_total: 0  # exploration, not verdict-FPs
script: Phase_EXPLORATION_no_yukawa.py
---

# Phase EXPLORATION — no-Yukawa hypothetical (γ-7 post-HALT-B)

> **STATUS:** EXPLORATION ONLY. **NOT verdict re-derivation.**
> **γ-7 HALT-B** remains LOCKED (user "Zatwierdzam Halt B" 2026-05-24).
> **Anti-Lakatos LOCK** preserved — pre-registered F-γ-7-A/B/C/D verdicts UNCHANGED.
> **F8 four-cycle FAIL pattern** (γ-3 + γ-3' + γ-5 + γ-7) UNCHANGED.
>
> This document is **information-gathering**, *not* a verdict revision. Findings inform future research direction (potential γ-8 if ever opened under fresh pre-registration), but **cannot** be cited as justification to retroactively reopen F8 within current scope.

---

## §1 — Motivation (user query 2026-05-24)

User question (Polish, verbatim):
> "rozumiem, że problemem jest ekranowanie yukawa, pytanie czy na pewno działa ono tak samo w skalach makro, to znaczy mamy mechanizm nakładania się gradientów pól ale może im większa masa, tym granica ekranowania też jest większa.
>
> Wydaje mi się, że to warto przeanalizować, jaki byłby wynik bez yakuwa, i jaki wymagany rozmiar wszechśiata, potem numerycznie sprawdzić proporcje i zobaczyć czy to może być realny mechanizm. **Na zasadzie eksploracji nie stwierdzeń**."

**Translation:** Yukawa screening assumed in Phase 1-2 was linearized KG result. User hypothesizes that at macro scales (massive systems), screening behavior might differ — possibly the screening radius scales with mass concentration. Asks: what would no-Yukawa give, what universe size would be required, can it be a realistic mechanism — **as exploration, not statement**.

This is a **legitimate scaling question**. It does NOT violate anti-Lakatos because:
1. It does not change γ-7 verdicts (HALT-B LOCKED).
2. It does not re-derive Phase 1/2/3 formulas (those remain LOCKED).
3. It is explicit "what-if" labeled hypothetical analysis.
4. Any future γ-8 cycle (if ever opened) would require fresh pre-registration, NOT inheriting this exploration's findings.

---

## §2 — Three hypothetical scenarios analyzed

| Scenario | Field kernel | IR cutoff | Phase 2 formula structure |
|----------|--------------|-----------|---------------------------|
| (A) Yukawa baseline | exp(-μr)/(4πr) | 1/μ_sp (natural) | V_eff/V = (G ρ² V ⟨exp⟩)/(2 μ_sp v²) · ξ |
| (B) Coulomb (no screening) | 1/(4πr) | R_universe (imposed) | V_eff/V = (2π/3) G ρ² R⁴/v² · ξ |
| (C) Hybrid (mass-dep. screening) | exp(-α·r/R)/(4πr) | α·R_universe | between (A) and (B) |

**Setup:**
- v_phi² = 3.03×10⁴⁵ J/m (Appendix E, Phase 2 LOCKED convention)
- ρ_m = Ω_m · ρ_crit = 2.73×10⁻²⁷ kg/m³
- R_obs = c/H₀ = 1.36×10²⁶ m (Hubble radius)
- ξ_clump ∈ {1 (best case), 10⁻⁴ (TGP-empirical from Phase 3)}

---

## §3 — (A) Yukawa baseline numerics

**Formula (Phase 2 LOCKED):** V_eff/V_univ = (G ρ² V_univ ⟨exp⟩_uniform)/(2 μ_sp v²) · ξ_clump

| R_universe | μR | ⟨exp⟩ | V_eff/V (ξ=1) | V_eff/V (ξ_TGP=10⁻⁴) |
|------------|-----|--------|----------------|----------------------|
| R_obs = λ_sp | 1.000 | 0.482 | 5.7×10⁻⁵ | 5.7×10⁻⁹ |
| 5·R_obs | 5.0 | small | grows | — |
| R → ∞ (saturated) | ∞ | → 6/(μR)³ | **7.1×10⁻⁴** (ceiling) | 7.1×10⁻⁸ |

**Key finding:** Yukawa saturation ceiling = 4π G ρ_m² λ_sp⁴/v² ≈ **7.1×10⁻⁴**.

> **No R_universe can break this ceiling.** Factor ~990 short of Ω_DE = 0.7.
> This is an **absolute barrier** for the Yukawa scenario.

---

## §4 — (B) Coulomb (no Yukawa) numerics

**Derivation:** Replace Yukawa pair-overlap integral
$$\int \delta\Phi_i \delta\Phi_j \, dV = \frac{q_i q_j \, e^{-\mu r_{ij}}}{8\pi \mu_{sp}} \quad\to\quad \frac{q_i q_j \, R_{universe}}{4\pi}$$

(Coulomb integral grows linearly with R cutoff — single-source self-pair gives ∫(q/4πr)² · 4πr² dr from 0 to R = q²R/(4π).)

**Formula:**
$$\boxed{\frac{V_{\text{eff}}}{V_{\text{universe}}} \approx \frac{2\pi}{3} \cdot \frac{G \, \rho^2 \, R_{\text{universe}}^4}{v^2} \cdot \xi_{\text{clump}}}$$

| R_universe | V_eff/V (ξ=1) | V_eff/V (ξ_TGP=10⁻⁴) |
|------------|----------------|----------------------|
| R_obs (Hubble) | 1.2×10⁻⁴ | 1.2×10⁻⁸ |
| 5·R_obs | 7.4×10⁻² | 7.4×10⁻⁶ |
| **8.77·R_obs** | **0.70** ✓ (ξ=1) | 7.0×10⁻⁵ |
| **87.7·R_obs** | 7.0×10⁴ (unphysical, theory breaks) | **0.70** ✓ (ξ_TGP) |

**Key finding:** Coulomb V_eff grows as R⁴ — **no saturation barrier**. Required R for Ω_DE = 0.7:
- **R ≈ 8.8 · R_obs** if ξ_clump = 1 (best case, but ξ=1 means full non-linear clumping)
- **R ≈ 88 · R_obs** if ξ_clump = 10⁻⁴ (TGP-empirical from Phase 3)

> Note: V_eff/V_univ > 1 (last row R=50R_obs gives 739) is **unphysical** — linearized perturbation framework breaks when V_eff approaches V_universe (analogous to R1 #17 runaway breakdown). The formula is reliable only in regime V_eff/V << 1.

---

## §5 — (C) Mass-dependent screening hypothesis (sketch)

**Physical motivation:** Nonlinear Mexican-hat Phi-substrate potential
$$V(\Phi) = \frac{\lambda}{4} \left(|\Phi|^2 - \Phi_0^2\right)^2$$

Linearized (Phase 1 used): m_sp² = ∂²V/∂Φ²|_{Φ=Φ_0} = 2λΦ_0² (constant).

**Nonlinear regime:** Near mass concentrations, ⟨|Φ|⟩ deviates from Φ_0. If ⟨|Φ|²⟩ < Φ_0² locally (mass "depletes" the Phi-condensate), then:
$$m_{sp,\text{eff}}^2 = 2\lambda \langle |\Phi|^2 \rangle < 2\lambda \Phi_0^2 \quad\Rightarrow\quad \lambda_{sp,\text{eff}} > \lambda_{sp}$$

If this depletion is significant on cosmological scales (where matter clumps create non-trivial Phi-substrate inhomogeneity), the effective Yukawa range could grow toward R_universe.

**Approximation:** Replace 1/μ_sp with α·R_universe, where α ∈ (0,1].

Required α for V_eff/V_univ = 0.7 (best case ξ=1):

| R_universe | α_required |
|------------|-----------|
| 1·R_obs | 8.77 (IMPOSSIBLE — α capped at 1) |
| 5·R_obs | 1.75 (still impossible) |
| 10·R_obs | 0.88 (extreme screening collapse) |
| 50·R_obs | 0.18 (modest screening relaxation) |

**Interpretation:** If universe is ~10-100× larger than observable (which inflation allows), even modest screening relaxation (α ~ 0.1-0.2) could give Ω_DE = 0.7.

---

## §6 — Comparison vs original γ-7 HALT-B verdict

**γ-7 HALT-B was justified on:**
1. F-γ-7-B FAIL: V_eff/V_univ at observable scale ~10⁻⁸ (factor 10⁷ short)
2. F-γ-7-C MAGNITUDE_FAIL: ä_γ-7 ~ 10⁻¹⁷ vs observed 10⁻¹⁰ (factor 10⁷ short)
3. F-γ-7-D FAIL: timing mismatch (Λ scaling vs structure-formation epoch)
4. R1 #17: TGP naive linear theory unphysical runaway δ growth

**Exploration findings DO NOT change any of these.** What they show:

| Hypothesis | Yukawa-saturated barrier? | Required exotic input | Status |
|------------|----------------------------|------------------------|--------|
| (A) Phase 2 Yukawa baseline | YES (ceiling 7×10⁻⁴) | none | HALT-B confirmed |
| (B) Coulomb + R ≈ 9 R_obs | NO (escapes) | (i) no Yukawa, (ii) ξ_clump=1, (iii) R = inflation-scale | requires three independent exotic assumptions |
| (C) Mass-dep. screening + R = 50 R_obs | NO (escapes) | (i) nonlinear KG, (ii) ξ_clump=1, (iii) α ≈ 0.18 | requires R1 #17 resolution + nonlinear extension |

**Bottom line:** γ-7 HALT-B remains correct **for the scope explicitly pre-registered** (linearized Phi-substrate, Yukawa screening from Appendix E eq. 353, R = c·t_age, ξ_clump from TGP-native R=c·t framework). The exploration identifies that **alternative scenarios require multiple independent exotic assumptions** that are NOT currently in TGP scope.

---

## §7 — What this means for future γ-8 (hypothetical)

If a future γ-8 cycle were ever opened (NOT authorized here), the **cleanest path** would be:

1. **Pre-register before Phase 1:** falsifiers for nonlinear KG predictions specifically (NOT recycle γ-7 falsifiers).
2. **Resolve R1 #17 first:** TGP cosmological perturbation theory must be made consistent (currently gives runaway δ). Without this, no quantitative comparison possible.
3. **Derive (not postulate) m_sp_eff(ρ, Φ):** From Phi-substrate Lagrangian explicitly, including nonlinear |Φ|² coupling. Must reduce to Appendix E eq. 353 in vacuum limit.
4. **Independent test of "large universe" claim:** CMB acoustic peaks + BAO constrain universe geometry. If inflationary R_universe >> R_obs, this is testable independently of γ-8.
5. **Anti-Lakatos discipline:** Must NOT cite this exploration as "γ-7 was incomplete" — γ-7 was complete within its pre-registered scope. γ-8 would need its own justification.

**Note on the four F8 FAIL pattern:** Even if a future γ-8 with mass-dependent screening resolved F8, the four-cycle FAIL pattern (γ-3 R=c·t, γ-3' E_P/ℓ_P, γ-5 quasi-static, γ-7 Yukawa) would still document that **F8 is structurally difficult** within current TGP framework. A single resolved γ-8 would not retroactively change γ-3/γ-3'/γ-5/γ-7 verdicts.

---

## §8 — Summary table (numerical)

| Scenario | V_eff/V_universe | Short of Ω_DE=0.7 |
|----------|-----------------:|------------------:|
| Yukawa, R = R_obs, ξ=1 | 5.7×10⁻⁵ | 1.2×10⁴ |
| Yukawa, R = R_obs, ξ_TGP | 5.7×10⁻⁹ | 1.2×10⁸ |
| Yukawa, R → ∞ (saturated), ξ=1 | **7.1×10⁻⁴** (ceiling) | **990** |
| Yukawa, R → ∞, ξ_TGP | 7.1×10⁻⁸ | 9.9×10⁶ |
| Coulomb, R = R_obs, ξ=1 | 1.2×10⁻⁴ | 5.9×10³ |
| Coulomb, R = R_obs, ξ_TGP | 1.2×10⁻⁸ | 5.9×10⁷ |
| Coulomb, R = 8.8 R_obs, ξ=1 | **0.70 ✓** | 1.0 (target hit) |
| Coulomb, R = 88 R_obs, ξ_TGP | **0.70 ✓** | 1.0 (target hit) |
| OBSERVED Ω_DE | 0.70 | 1.0 |

---

## §9 — Discipline statement (RE-AFFIRMED)

This exploration is provided AS-IS, labeled **information-only**.

**It does NOT:**
- Modify γ-7 HALT-B claim_status (LOCKED 2026-05-24)
- Modify F-γ-7-A/B/C/D pre-registered verdicts (LOCKED)
- Modify Phase 1/2/3/4/5/FINAL LOCKED derivations
- Affect F8 four-cycle FAIL pattern (γ-3 + γ-3' + γ-5 + γ-7)
- Authorize any γ-8 cycle (would require fresh user authorization)

**It DOES:**
- Provide scaling intuition for what assumptions WOULD be required to rescue F8
- Quantify "Yukawa saturation ceiling" (~7×10⁻⁴) as absolute Phase 2 barrier
- Quantify "Coulomb R⁴ growth" as escape path requiring R ≈ 9-88 × R_obs
- Identify R1 #17 (linear theory limitation) and m_sp(ρ) (density-dependent screening) as the two key open questions for any future cycle

**Future research note (R1 #17 extended):**
The TGP perturbation framework would benefit from:
- Nonlinear KG analysis (Mexican-hat potential effects on m_sp)
- Independent constraint on R_universe (CMB + BAO geometry)
- Reformulation of structure-formation framework consistent with R=c·t

These are **research scope items**, not current cycle deliverables.

---

**END EXPLORATION**

**γ-7 closure status:** HALT-B LOCKED (UNCHANGED).
**Cross-cycle propagation:** None required (exploration does not affect STATE.md or PRE_REGISTERED_FALSIFIERS.md verdicts).
**Optional annotation:** R1 #17 in meta could be extended with scaling-hypothesis note (user discretion).
