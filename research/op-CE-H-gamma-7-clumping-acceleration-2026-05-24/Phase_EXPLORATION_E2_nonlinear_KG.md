---
title: "Phase EXPLORATION E2 — nonlinear KG + cosmic-web channeling (γ-7 post-HALT-B)"
type: phase_exploration
status: EXPLORATION_ONLY
phase: EXPLORATION_E2 (post-HALT-B, post-E1)
parent_cycle: op-CE-H-gamma-7-clumping-acceleration-2026-05-24
parent_close: Phase_FINAL_close.md
parent_e1: Phase_EXPLORATION_no_yukawa.md
created_date: 2026-05-24
authorization: "User 2026-05-24: 'tak działaj z Exploration Phase E2, zobaczymy co tutaj wyjdzie ;0'"
halt_b_status: LOCKED_UNCHANGED
anti_lakatos_lock: PRESERVED
modifies_verdicts: NO
modifies_phase_locks: NO
purpose: "Test user hypothesis on density-dependent screening via TGP-native nonlinearity"
substantive_fp_total: 0  # exploration
script: Phase_EXPLORATION_E2_nonlinear_KG.py
---

# Phase EXPLORATION E2 — nonlinear KG + cosmic-web channeling

> **STATUS:** EXPLORATION ONLY. **NOT γ-8 pre-registration.**
> **γ-7 HALT-B** LOCKED unchanged.
> Anti-Lakatos LOCK preserved.
>
> This document tests user hypothesis (2026-05-24) that Yukawa screening behaves differently at macro scales due to nonlinear matter-Phi coupling. Uses TGP-native Lagrangian N[Φ] from `sek02_pole.tex` — does NOT introduce new postulates.

---

## §1 — User hypothesis (verbatim)

> "Rozchodzenie się fal grawitacyjnych to wymuszenie zmiany konfiguracji przestrzeni własnej cząstki względem silnie skonfigurowanej przestrzeni masywnego obiektu.
> [...] taki obiekt zmienia własną konfigurację, to wymusza też konfigurację innych obiektów. Na bardzo dużych odległościach wymuszenie tej konfiguracji faktycznie jest ekranowane (w ramach mocy tego źródła), ale grawitacja w TGP jest silnie skorelowana z tym jak obiekt może 'ustawić swoją własną przestrzeń'. Im dalej od silnego źródła tym to wymuszenie jest słabsze, ale jednocześnie powinno propagować przez inne źródła.
> [...] jeżeli mamy 50 gwiazd ustawionych w linii i losowo wycinamy jakieś gwiazdy — propagacja się zmienia."

**Reformulation in physics language:**
- Yukawa screening in vacuum (Phase 1-3): linearized result, valid only when nonlinear terms in N[Φ] are negligible
- In presence of matter (cosmic web), nonlinear Phi-self-coupling provides "relay paths" through intervening masses
- Effective propagator becomes density-dependent: m_sp_eff(ρ) ≠ const

This is a known physics mechanism in scalar-tensor theories (analog: chameleon mechanism, but with opposite sign for anti-screening).

---

## §2 — TGP-native source (sek02_pole eq. N-preview)

The full TGP field equation includes nonlinear operator (sek02_pole.tex eq:N-preview):

$$\mathcal{N}[\Phi] = \frac{\alpha}{\Phi_0}\frac{(\nabla\Phi)^2}{\Phi} + \beta\frac{\Phi^2}{\Phi_0^2} - \gamma\frac{\Phi^3}{\Phi_0^3}$$

with:
- α = 2 (derived from action variation)
- β, γ > 0
- [β] = [γ] = L⁻²
- vacuum condition (N[Φ_0] = 0): **β = γ**
- m_sp² = γ ≈ H_0²/c² (Appendix E eq. 353)

**Important:** these terms are **already in core TGP**. Phase 1 γ-7 used **linearized** approximation (kept only quadratic term), which deliberately dropped these nonlinearities for tractability. Phase E2 restores them.

---

## §3 — Symbolic linearization (sympy executed)

Expanding N[Φ] around Φ = Φ_0 + δΦ:

| Order in δΦ | Algebraic part (β=γ) | Kinetic α-part (α=2) |
|-------------|-----------------------|----------------------|
| 0 (constant) | **0 ✓** (vacuum) | 0 |
| 1 (linear) | -γ·δΦ/Φ_0 (gives m_sp²=γ) | 0 |
| 2 (quadratic) | -2γ·(δΦ)²/Φ_0² | +(2/Φ_0²)·(∇δΦ)² |
| 3 (cubic) | -γ·(δΦ)³/Φ_0³ | -(2/Φ_0³)·(∇δΦ)²·δΦ |

**Cubic self-coupling** (relevant for relay paths):
$$\sigma_{\text{eff}} = \frac{2\gamma}{\Phi_0}$$

**Quartic self-coupling**:
$$\kappa_{\text{eff}} = \frac{\gamma}{\Phi_0^2}$$

These are **TGP-native**, NOT postulated.

---

## §4 — Effective mass in matter background (KEY RESULT)

Splitting δΦ = ⟨δΦ⟩_bg + δ²Φ (fluctuation on top of matter-induced background), linearizing in δ²Φ:

$$\boxed{m_{sp,\text{eff}}^2(\chi) = \gamma \cdot (1 + 4\chi + 3\chi^2) = \gamma \cdot (1+\chi)(1+3\chi)}$$

where **χ ≡ ⟨δΦ⟩_bg/Φ_0** is the dimensionless field deviation from vacuum.

**Sign behavior:**

| χ range | m_sp_eff² / γ | Behavior | Physical interpretation |
|---------|----------------|----------|--------------------------|
| χ > 0 | > 1 | Stronger screening | Chameleon-like (matter enhances mass) |
| 0 > χ > -1/3 | < 1 | **Anti-screening** | **User hypothesis ✓** |
| χ = -1/3 | 0 | Screening collapse | λ_eff → ∞ |
| -1/3 > χ > -1 | < 0 | **Tachyonic regime** | Unphysical |
| χ < -1 | (1+χ)<0 | Sign-flipped | Non-perturbative regime |

**TGP K2 ontology** (`sek02_pole.tex` §single source): inside matter, Φ < Φ_0 (matter depletes the Phi-condensate). So **χ < 0 in matter regions** — exactly the user's intuition.

**This is the qualitative confirmation of the hypothesis from TGP-native Lagrangian.**

---

## §5 — Numerical: what χ_web is required?

From Phase 2 LOCKED Yukawa-saturated formula:
$$\frac{V_{\text{eff}}}{V_{\text{universe}}}\bigg|_{\text{sat}} = \frac{4\pi G \rho_m^2 \lambda_{sp,\text{eff}}^4}{v_\phi^2}$$

Substituting λ_sp_eff = 1/m_sp_eff = (1+χ)^(-1/2)·(1+3χ)^(-1/2)·λ_sp:

| χ_web | m_eff²/γ | λ_eff/λ_sp | V_eff/V_universe (ξ=1, saturated) |
|-------|----------|-------------|------------------------------------|
| 0 | 1.000 | 1.00 | 7.1×10⁻⁴ |
| -0.05 | 0.808 | 1.11 | 1.1×10⁻³ |
| -0.10 | 0.630 | 1.26 | 1.8×10⁻³ |
| -0.20 | 0.320 | 1.77 | 6.9×10⁻³ |
| **-0.30** | **0.070** | **3.78** | **1.4×10⁻¹** |
| **-0.32** | **0.027** | **6.06** | **0.96 ≈ Ω_DE ✓** |

**Required: χ_web ≈ -0.31 to -0.33** to reach V_eff/V_universe = 0.7.

This requires **substantial Phi-substrate depletion in cosmic web regions**.

Is χ_web ≈ -0.3 physically achievable? **OPEN** — requires nonlinear cosmic-web simulation.

---

## §6 — Lab-scale Newton matching (γ-5 inheritance test)

**CRITICAL TENSION:** If anti-screening operates with χ ≈ -0.3 cosmologically, what happens at lab scale where matter density is much higher?

Naive estimate (using same χ-vs-ρ scaling):
- ρ_Earth/ρ_cosmic ≈ 2×10³⁰
- If χ scales with ρ linearly: χ_lab ≫ 1 → **regime collapse**

**This would destroy lab-scale Newton gravity** unless:
- (a) **Φ_0 calibration** is such that χ_lab ≪ 10⁻⁴ (precision of inverse-square tests) but χ_web ~ 0.3
- (b) **Saturation effect**: χ doesn't scale linearly with ρ; nonlinear amplification activates only above critical density (cosmic-web threshold)
- (c) **Coherence effect**: only **clumping** (not raw density) drives χ; lab is "smooth" so χ_lab ≈ 0 while cosmic web is "clumped" so χ_web ≠ 0

Option (c) is **physically motivated** — V_eff itself is about clumping in γ-7 framework. Anti-screening from clumping (not density) is consistent.

But **quantitatively**, we cannot resolve this without:
1. Explicit Φ_0 SI calibration (open in TGP — `sek02_pole.tex` only gives dimensionless g_0 = 0.87)
2. Coupled nonlinear simulation of cosmic-web Phi-substrate
3. Independent constraint from lab gravity tests

---

## §7 — Relay-path analysis (50 stars in line)

User's example: 50 stars in line, A → ... → Z. Each star = vertex with cubic coupling σ_eff = 2γ/Φ_0.

**Direct amplitude:** ~ q²·exp(-μR_AZ)/(4π·R_AZ)
**Single relay (B midway):** ~ q²·V₃·[exp(-μL/2)/(2πL)]²
**Multi-relay (all 49 hops):** ~ q²·V₃⁴⁸·[exp(-μL)/(4πL)]⁴⁹

Critical observation: **exp factors are SAME** (path length identical: R_AZ = 49L either way).

**Difference is in power-law prefactor:**
$$\frac{\text{relay}}{\text{direct}} = \left(\frac{V_3}{4\pi L^2}\right)^{n-2} \cdot \frac{R_{AZ}}{L}$$

For V_3 = 2γ/Φ_0 = 2μ²/Φ_0:

| Φ_0 calibration | Relay/Direct factor for n=50 stars |
|------------------|-------------------------------------|
| Φ_0 = 1 | ~ 1 (no enhancement) |
| Φ_0 = 0.1 | ~ 10⁴⁸ (relay dominates massively) |
| Φ_0 = 0.01 | ~ 10⁹⁶ (overwhelming) |

**Extreme sensitivity to Φ_0**. Without calibration, this is undetermined.

**User's "wycinamy gwiazdy" test:** if relay paths dominate, removing intermediate stars **drastically reduces** the felt gravity from distant ones. This is **observationally testable** (anisotropy of gravity along cosmic filaments vs voids).

---

## §8 — Honest assessment

### Qualitatively: hypothesis CONFIRMED by TGP-native Lagrangian

The TGP nonlinear operator N[Φ] from `sek02_pole.tex`:
1. **Does give** cubic + quartic self-couplings (σ_eff, κ_eff exist in TGP)
2. **Sign IS anti-screening** in matter regions (χ < 0)
3. **Cosmic-web channeling IS possible** via relay paths
4. **Independent observable predictions** exist (filament-aligned gravity anisotropy)

User's intuition is **physically aligned with TGP core**, not speculative addition.

### Quantitatively: viability is OPEN

Required for F8 acceleration via this mechanism:
- χ_web ≈ -0.3 (cosmologically averaged in cosmic web)
- χ_lab ≪ 10⁻⁴ (preservation of Newton at lab scale)
- These two require **nonlinear amplification mechanism** that selectively kicks in for clumped matter (cosmic web), not uniform matter (lab)

**Three open questions for any future work:**
1. **Φ_0 SI calibration** — currently only dimensionless g_0 = 0.87 is fixed
2. **Cosmic-web χ_RMS computation** — requires nonlinear simulation, not analytical
3. **Lab/cosmic-web bifurcation mechanism** — what makes anti-screening density-clumping-dependent rather than just ρ-dependent

### Verdict for hypothetical γ-8 candidacy

**Status: PROMISING SCAFFOLD but NOT yet enough for γ-8 authorization.**

Would need (minimum) before γ-8 could be cleanly opened:
- (P1) Φ_0 calibration via extended γ-5 Newton matching (lab-scale)
- (P2) Toy-model nonlinear Phi-substrate simulation showing χ_web grows with cosmic web formation
- (P3) Prediction of observable filament-aligned gravity anisotropy (independent test)
- (P4) Solution to R1 #17 (linear theory runaway) via nonlinear corrections

These are **major research items**, not minor extensions.

---

## §9 — What this means for γ-7 closure

**γ-7 HALT-B remains LOCKED.** Phase 1-5 derived linearized framework correctly; the linearized framework genuinely cannot deliver Ω_DE = 0.7. That is the **honest closure** of the linearized pair-overlap approach.

What E2 shows is that **the linearization was the critical assumption**, not the mechanism itself. Restoring TGP-native nonlinearity opens new physics that the pre-registered γ-7 protocol explicitly excluded (Phase 1 §scope was linearized).

**Distinction (anti-Lakatos):**
- γ-7 (linearized pair-overlap in vacuum) — CLOSED, FAIL
- Hypothetical γ-8 (nonlinear cosmic-web channeling) — NEW MECHANISM, **distinct physics**, would require fresh pre-registration

Citing E2 as "γ-7 was incomplete" is **wrong** — γ-7 was complete within its pre-registered scope. Citing E2 as "γ-7 motivates γ-8" is **also wrong** — these are independent investigations of different physical mechanisms.

The correct framing: **E2 is exploratory mapping of nonlinear TGP regime**, useful for any future cosmological investigation in TGP, NOT specifically tied to γ-7's verdict.

---

## §10 — Recommended next steps (NOT proposing γ-8)

If user wants to pursue this thread, the **legitimate path** is:

### Track 1: TGP Lagrangian work (independent of cosmology)
- Calibrate Φ_0 in SI via γ-5 nonlinear Newton matching test
- This is a **standalone TGP foundations** task (improves γ-3'/γ-5 inheritance)
- Output: explicit Φ_0 = X kg/m^? (or whatever units)

### Track 2: Toy-model nonlinear N[Φ] simulation
- 1D or 2D toy: Phi-substrate with point sources + N[Φ] full nonlinearity
- Compute ⟨δΦ⟩ vs density, anti-screening regime mapping
- Standalone numerical TGP investigation

### Track 3: Observable predictions (independent test)
- Compute predicted filament-aligned gravitational anisotropy
- Compare to existing precision tests of weak lensing along/across filaments
- Falsifiable independently of F8 scope

**Only after Tracks 1-3 have substantive findings** would γ-8 pre-registration become well-defined.

---

## §11 — Summary

| Aspect | Finding | Confidence |
|--------|---------|-----------|
| TGP gives cubic/quartic self-couplings | ✓ Yes, in sek02 N[Φ] | HIGH (derived) |
| Sign of correction in matter | ✓ Anti-screening (χ<0) | HIGH (K2 ontology) |
| Magnitude χ ≈ -0.3 achievable | ? Open | LOW (needs calibration) |
| Lab-scale compatibility | ⚠ Tension | OPEN (needs Φ_0 calibration + clumping vs density distinction) |
| Relay-path enhancement | ✓ Possible, ? quantitatively | OPEN (Φ_0 dependent) |
| Justifies γ-8 authorization | ✗ Not yet | (premature) |

**Bottom line:** User's hypothesis is **physically aligned with TGP core** and identifies a **legitimate gap** in Phase 1's linearization. But quantitative viability requires substantial further work (Tracks 1-3) before any γ-8 cycle could be cleanly pre-registered.

γ-7 closure (HALT-B) remains LOCKED.

---

**END EXPLORATION E2**

**Status:** information-gathering complete. Findings logged. No verdict changes. No new authorization implied.
