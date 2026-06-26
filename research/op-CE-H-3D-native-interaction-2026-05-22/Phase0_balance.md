---
title: "Phase 0 — Balance sheet + ANALYTICAL PRE-DERIVATION (op-CE-H-3D-native-interaction-2026-05-22)"
type: phase_balance
status: LOCKED
pre_registration_date: 2026-05-22
phase: 0
parent_cycle: op-CE-H-3D-native-interaction-2026-05-22
methodology_note: "First cycle post CALIBRATION_PROTOCOL §3.6 BINDING (analytical pre-derivation step required dla FP-class falsifiers). §3.6 enforcement visible w §8 below."
---

# Phase 0 — Balance sheet + Analytical pre-derivation

**Status:** LOCKED 2026-05-22.
**Purpose:** Explicit accounting wszystkich inputów (external + structural axioms) + analytical pre-derivation dla F-γ-1 numerical threshold (CALIBRATION §3.6 BINDING). Anti-Lakatos discipline.

---

## §1 — External inputs (co przyjmujemy z zewnątrz)

| ID | Input | Source | Status |
|----|-------|--------|--------|
| EXT-1 | TGP Phi-substrate Lagrangian: $\mathcal{L} = |\partial_\mu \Phi|^2 - (\lambda/4)(|\Phi|^2-v^2)^2$ | meta/TGP_GENERATED_SPACE_COSMOLOGY §3.2 + FFS Phase 3 | DERIVED w previous cycles |
| EXT-2 | Mexican hat potential explicit | Pattern 2.5 §3.5.6; FFS Phase 3 | DERIVED |
| EXT-3 | 3D Euclidean spatial geometry (Minkowski static limit) | Standard QFT; concept paper §3.4 | METHODOLOGICAL CHOICE |
| EXT-4 | S05 U(1) phase as **global** symmetry | TGP_FOUNDATIONS §3.4 (S05 phase mechanism); U(1)_em emergent | THIS CYCLE INTERPRETATION (see §4) |
| EXT-5 | AX-Z2 discrete symmetry preserved (Φ → -Φ) | TGP_FOUNDATIONS §2 | LOCKED |
| EXT-6 | AX-RP2 topology preserved | TGP_FOUNDATIONS §4 | LOCKED |
| EXT-7 | Standard 3D Laplace + Helmholtz Green functions | Standard QFT (Peskin-Schroeder, Weinberg I) | LITERATURE |
| EXT-8 | Sympy 1.12+ symbolic + numerical | Tool dependency | INFRASTRUCTURE |

**Critical note EXT-4:** The interpretation of S05 U(1) (global vs local gauge) directly determines whether Goldstone mode is massless (global) or eaten by gauge field (local). **This cycle pre-registers GLOBAL U(1) interpretation** per minimal axioms (no separate gauge field axiom; gauge field is EFFECTIVE EMERGENT per Foundations §3.4). Local U(1)_em emerges; bare S05 axiom is global. See §4 analytical pre-derivation.

---

## §2 — LOCKED structural axioms (preserved during cycle)

| ID | Axiom | Source | Role |
|----|-------|--------|------|
| AX-S05 | Single scalar Phi z U(1) phase | Foundations §1 | Lagrangian + Goldstone mode source |
| AX-Z2 | Discrete Z₂ symmetry Phi → -Phi | Foundations §2 | Vacuum structure (mexican hat) |
| AX-U1 | U(1) phase symmetry (GLOBAL in basic axiomatization) | Foundations §3 | **Critical: massless Goldstone mode if global** |
| AX-RP2 | RP² topology | Foundations §4 | Hedgehog/orientation field; **used jeśli DEC switch to point defect** |
| AX-DECL-1 | SU(2)_L + SU(3)_c gauge limit | meta/TGP_W_Z_THEORETICAL_LIMIT | PRESERVED (NIE bypassed) |
| AX-DECL-2 | Φ_0_local absolute NOT derivable z minimal axioms | FFS Phase 4 + R3 trigger | PRESERVED |
| AX-CE-STR | CE-H structural feature (R3 3/3 accepted 2026-05-21) | R2 audit 2026-05-22 | TESTED quantitatively w 3D |

---

## §3 — Derived outputs (co claimuje wyprodukować)

| ID | Output | Phase | Pre-registered prediction |
|----|--------|-------|--------------------------|
| OUT-1 | 3D vortex ansatz $\Phi_{vortex}(r,\phi) = \rho(r) e^{in\phi}$ z far-field expansion | Phase 1 | $\rho \to v$ at r→∞; massless Goldstone +massive Higgs decomposition |
| OUT-2 | Mass spectrum: $m_\sigma$ (Higgs) + 0 (Goldstone massless) | Phase 1 | $m_\sigma = \sqrt{2\lambda}v$; Goldstone exactly massless dla global U(1) |
| OUT-3 | Two-defect interaction V_int(L) z propagator analysis | Phase 2 | **Power-law (1/L Coulomb dla point defects) lub log(L)·L_z dla parallel vortices** — see §8 analytical pre-derivation |
| OUT-4 | R²_power vs R²_exp fit comparison | Phase 3 | R²_power > 0.99 at L >> 1/m_σ; R²_exp clearly inferior |

---

## §4 — Tautology test (anti-circular reasoning)

**Question:** Czy zakładamy F-γ-1 PASS przed sympy?

**Pre-registered answer:**

**Co zakładamy (z minimal axioms + structural argument):**
- S05 single Phi field
- Z₂ + mexican hat → SSB
- U(1) phase symmetry **global** w basic axiomatization (justified §4 EXT-4 + §8 analytical)

**Co testujemy (NIE assumed):**
- Czy SSB faktycznie produces massless Goldstone mode w 3D (Phase 1 EL analysis verifies)
- Czy two-defect coupling przez Goldstone faktycznie daje 1/L tail (Phase 2 propagator analysis)
- Czy fit R²_power > R²_exp at large L (Phase 3 differential test)

**Anti-tautology check:**
- Goldstone "expected" z general principle (Goldstone theorem) — BUT explicit verification w TGP-native Lagrangianu Phase 1 wymagana
- 1/L Coulomb "expected" z 3D massless propagator — BUT explicit two-defect calculation Phase 2 wymagana
- Coupling strength + sign NIE pre-registered (compute from sympy)

**Result:** Phase 1-3 są **falsifiable** even with strong analytical expectation. F-γ-1 PASS NIE jest auto-result of pre-registration.

---

## §5 — Falsifiability test

**Question:** Czy F-γ-1 jest falsifiable?

| Output | Falsifiable? | How |
|--------|--------------|-----|
| OUT-1 (vortex ansatz) | YES | Ansatz must satisfy EL equations; jeśli nie, sympy reveals algebraic inconsistency |
| OUT-2 (mass spectrum) | YES | Jeśli Goldstone gets mass z higher-order corrections, F-γ-1 may FAIL |
| OUT-3 (V_int(L)) | YES | Jeśli V_int(L) ~ exp(-mL) bez power-law modulation, F-γ-1 FAIL → HARD HALT scenario |
| OUT-4 (R² comparison) | YES | Jeśli R²_exp > R²_power, F-γ-1 FAIL |

All 4 outputs są falsifiable z explicit threshold.

**HARD HALT scenario explicit pre-registered:** Jeśli sympy reveals że w 3D TGP-native Lagrangianu BRAK massless mode (Goldstone gets mass from some mechanism not anticipated), OR że coupling between defects vanishes, V_int(L) będzie pure exponential. To NIE jest cycle failure — to jest **honest structural finding** wymagający fundamental redesign CE-H framework.

---

## §6 — Anti-BD-drift check (CRITICAL)

**Question:** Czy w jakimkolwiek miejscu fitujemy do Nielsen-Olesen / Vilenkin-Shellard / GR cosmology frameworks?

**Pre-registered answer:** NIE.

- **Vortex ansatz** używana jako native TGP-native S05+Z₂ defect, NIE adopted z Nielsen-Olesen cosmic string framework. Cosmic string theory jest informational anchor.
- **3D propagator analysis** native z Mexican hat V(Φ) + SSB. NIE adopted z Abelian Higgs (which has gauge field — our cycle pre-registers GLOBAL U(1)).
- **Coulomb 1/L expected form** z standard 3D massless scalar propagator (textbook QFT). NIE fitting do observed phenomenology.

**Methodology:** Native equations FIRST. Comparison z literature = post-hoc bonus / informational only.

**User explicit (2026-05-21):**
> "tworzymy natywne równania TGP, sprawdzamy czy wyniki są zgodne z pomiarami, a opcjonalnie robimy mapowanie"

---

## §7 — Independent-path cross-validation

**Pre-registered cross-checks per output:**

- **Path A — Analytical:** Linearized EL equations + Green function method (analytical closed-form expected for far-field)
- **Path B — Numerical:** Sympy nsolve dla ρ(r) profile + numerical integration of two-defect E(L)
- **Path C — Limit checks:** L → ∞ limit (free defects, V_int → 0); L → 0 limit (overlapping defects, breakdown of ansatz)
- **Path D — Mass scaling:** V_int(L)/v² should be dimensionally fixed; check parameter dependence

**Conflict resolution:**
- Path A vs B conflict → flag + investigate (R1 research-tier permissive)
- Path C/D limit fail → structural problem, escalate

---

## §8 — Analytical pre-derivation step (CALIBRATION §3.6 BINDING)

**CRITICAL — first cycle post-§3.6 BINDING 2026-05-22.** This section satisfies new requirement: explicit analytical derivation BEFORE any sympy.

### §8.1 — Setup z TGP-native Lagrangianu

$$\mathcal{L}_{TGP} = |\partial_\mu \Phi|^2 - \frac{\lambda}{4}(|\Phi|^2 - v^2)^2$$

z $\Phi: \mathbb{R}^3 \to \mathbb{C}$ (single Phi field z U(1) phase per AX-S05).

### §8.2 — SSB decomposition

After SSB: $\langle \Phi \rangle = v$ (choosing real positive WLOG, breaking U(1) phase).

Decompose around VEV:
$$\Phi(x) = (v + \sigma(x)) \cdot e^{i\theta(x)/v}$$

(Phase normalized by v dla dimensional consistency.)

Substituting do Lagrangianu:

$$\mathcal{L} = \frac{1}{2}(\partial \sigma)^2 + \frac{1}{2}\left(1 + \frac{\sigma}{v}\right)^2 (\partial \theta)^2 - V(\sigma)$$

z potential expanded:

$$V(\sigma) = \frac{\lambda}{4}\left((v+\sigma)^2 - v^2\right)^2 = \frac{\lambda}{4}(2v\sigma + \sigma^2)^2 = \lambda v^2 \sigma^2 + O(\sigma^3, \sigma^4)$$

### §8.3 — Quadratic Lagrangian (linearized)

Keeping only quadratic terms:

$$\mathcal{L}_{quad} = \frac{1}{2}(\partial\sigma)^2 - \lambda v^2 \sigma^2 + \frac{1}{2}(\partial\theta)^2$$

**Mass spectrum (analytical prediction):**

| Mode | Field | Mass² | Mass |
|------|-------|-------|------|
| Higgs (radial) | σ | $2\lambda v^2$ | $m_\sigma = v\sqrt{2\lambda}$ |
| Goldstone (phase) | θ | 0 | **MASSLESS** |

**Goldstone theorem application:** Continuous global U(1) symmetry broken by ⟨Φ⟩ = v → 1 massless mode (Goldstone). Confirmed analytically z quadratic Lagrangian.

### §8.4 — 3D propagators (analytical)

Static (time-independent) Green functions in 3D Euclidean:

**Goldstone (massless):**

$$\left(-\nabla^2\right) G_\theta(\vec{r}) = \delta^3(\vec{r})$$
$$\Rightarrow G_\theta(\vec{r}) = \frac{1}{4\pi |\vec{r}|}$$

**Higgs (massive):**

$$\left(-\nabla^2 + m_\sigma^2\right) G_\sigma(\vec{r}) = \delta^3(\vec{r})$$
$$\Rightarrow G_\sigma(\vec{r}) = \frac{e^{-m_\sigma |\vec{r}|}}{4\pi |\vec{r}|}$$

### §8.5 — Defect-defect interaction (analytical pre-derivation)

**Setup:** Two static defects at positions $\vec{x}_1, \vec{x}_2$ with separation $L = |\vec{x}_1 - \vec{x}_2|$. Each defect sources field $\theta$ (phase winding) and/or $\sigma$ (radial deviation).

**Point defect (hedgehog/monopole-like) z phase winding $n_i$:**

At far-field (r >> r_core), each defect i contributes to phase $\theta$:
$$\theta_i(\vec{r}) \approx n_i \cdot \Omega(\vec{r} - \vec{x}_i)$$

where Ω is angular function (specific form depends on defect geometry; for monopole-like spherically symmetric, Ω is solid-angle integral).

**Interaction energy at far-field (linear superposition):**

$$E_{int}(L) = \int d^3r \left[\nabla\theta_1 \cdot \nabla\theta_2\right] = -\int d^3r \theta_1 \nabla^2 \theta_2$$

For 3D Coulomb-like point defects coupled to Goldstone:
$$E_{int}(L) \sim n_1 n_2 \cdot v^2 \cdot G_\theta(L) = \frac{n_1 n_2 v^2}{4\pi L}$$

**EXPECTED ANALYTICAL FORM:**

$$\boxed{V_{int}(L) \approx -\frac{C \cdot v^2 \cdot n_1 n_2}{L} + O\left(\frac{e^{-m_\sigma L}}{L}\right)}$$

z C dimensionless coupling depending on defect specific structure (estimated O(1) — sign opposite z anti-aligned windings).

**Higgs (massive) contribution:** Exponentially suppressed at L >> 1/m_σ. Subdominant at large L.

**Goldstone (massless) dominant:** **1/L Coulomb-like power-law** dominant at all L >> r_core.

### §8.6 — Parallel vortex lines (alternative ansatz; DEC reserved)

Jeśli DEC switch to vortex line ansatz (z-translation invariant), parallel vortices at $\vec{x}_{1,\perp} = -L \hat{x}/2$, $\vec{x}_{2,\perp} = +L\hat{x}/2$:

Each vortex sources phase $\theta$ z 2D winding:
$$\theta_i(\vec{r}_\perp) = n_i \cdot \phi_i(\vec{r}_\perp)$$

where $\phi_i$ jest angle around vortex i in perpendicular plane.

**Interaction energy PER UNIT LENGTH (2D problem):**

$$\frac{E_{int}}{L_z}(L) = v^2 \int d^2 r_\perp \nabla\phi_1 \cdot \nabla\phi_2$$

For 2D Laplacian: $\nabla^2 \phi_i = 2\pi n_i \delta^2(\vec{r}_\perp - \vec{x}_{i,\perp})$, and 2D Green function $G_{2D}(\vec{r}_\perp) = -\frac{1}{2\pi}\log(|\vec{r}_\perp|)$.

$$\frac{E_{int}}{L_z}(L) \approx 2\pi v^2 n_1 n_2 \cdot \log(L/r_0)$$

z r_0 IR cutoff (string core size).

**EXPECTED ANALYTICAL FORM (vortex ansatz):**

$$\boxed{\frac{V_{int}(L)}{L_z} \approx 2\pi v^2 n_1 n_2 \log(L/r_0)}$$

**Logarithmic per unit length.** Dla finite vortex length, integration gives log(L) × (vortex length).

### §8.7 — F-γ-1 expected verdict (analytical pre-derivation)

**Pre-registered analytical prediction (LOCKED 2026-05-22):**

1. **Point defect ansatz (DEC-default):** V_int(L) ~ 1/L Coulomb-like (power-law β=1)
2. **Vortex line ansatz (DEC-alternative):** V_int(L)/L_z ~ log(L/r_0) logarithmic

**BOTH satisfy F-γ-1 PASS criterion** (power-law OR logarithmic — clearly distinguishable from pure exponential).

**Subdominant exponential corrections** (from massive Higgs) at L >> 1/m_σ:
$$V_{int}^{Higgs}(L) \sim \frac{e^{-m_\sigma L}}{L}$$

Are present BUT subdominant. Dominant tail jest power-law/log.

### §8.8 — F-γ-1 sanity check vs Poziom β 1D Z2

**1D Z2 result (Poziom β Phase 3):** V_int(L) ~ exp(-m·√2·L), R² = 0.9999 exponential.

**Why 1D differs from 3D:**

In 1D, Z2 is discrete (no continuous symmetry) → no Goldstone mode. Both modes (kink and antikink) are MASSIVE. 1D static propagator for massive field is exp(-m|x|)/(2m) — exponential, no power-law modulation.

In 3D, U(1) continuous → Goldstone massless → 1/L Coulomb.

**This is dimensional + symmetry effect.** Pre-registered analytical prediction:

$$\frac{V_{int}^{3D}(L)}{V_{int}^{1D}(L)} \sim e^{+m\sqrt{2}L}/L \to \infty \text{ at } L \to \infty$$

**3D power-law DOMINATES over 1D exponential** at all sufficient L. Differential test Phase 3 should confirm.

### §8.9 — Caveats explicit (CALIBRATION §3.6.2 forbidden shortcut check)

**Forbidden shortcut check:** Czy używam "heuristic intuition" zamiast analytical derivation?

| Step | Heuristic OR Analytical | Justification |
|------|-------------------------|---------------|
| SSB decomposition | Analytical | Standard QFT method (Goldstone theorem) |
| Quadratic expansion | Analytical | Standard mexican hat linearization |
| 3D Green functions | Analytical | Standard Helmholtz/Laplace Green functions |
| Defect coupling form | **Heuristic** (linear superposition assumption) | Phase 1 sympy must verify linearity |
| 1/L tail | Analytical | Direct consequence of Goldstone + 3D propagator |
| Coupling constant C | **Heuristic** (estimated O(1)) | Phase 2 sympy computes explicitly |

**Honest declaration:** Two steps (defect coupling form + coupling constant) are heuristic — Phase 1/2 sympy MUST verify. Section §8.7 prediction is **EXPECTED form**; Phase 1-3 sympy COMPUTES from Lagrangianu explicitly.

**Anti-Lakatos preservation:** If Phase 1-3 sympy gives different form (e.g., coupling constant zero accidentally, or higher-order suppression), F-γ-1 FAIL outcome IS legitimate, NIE rescued by appealing to §8 analytical prediction. **§8 sets EXPECTATION**; Phase 1-3 sympy is the **VERIFICATION**.

### §8.10 — Pre-registered F-γ-1 numerical threshold (post-§3.6)

Per CALIBRATION §3.6.1 requirement:
- (a) **Analytical derivation:** §8.1-§8.7 explicit derivation z TGP-native Lagrangianu
- (b) **Symbolic computation:** §8.4-§8.6 explicit Green function + interaction formulas
- (c) **Phase 0 documentation:** This entire §8 section

**Numerical threshold F-γ-1:**
- Phase 3 fit R² for V_int(L) at L > 3/m_σ:
  - **R²_power ≥ 0.95 AND R²_power > R²_exp + 0.02** → F-γ-1 PASS
  - R²_power < 0.95 OR R²_power ≤ R²_exp → F-γ-1 PARTIAL or FAIL
- Expected (analytical): R²_power ≈ 0.999 (excellent Coulomb fit), R²_exp ≈ 0.95 (poor for power-law data)

**Bound for expected coupling constant C:**
- Analytical: C ∈ [0.1, 10] (O(1) estimate, depends on defect ansatz specifics)
- Phase 2 sympy result outside this range → R1 flag (NIE FAIL automatic; investigate)

---

## §9 — Choice of ansatz + DEC budget pre-registration

**Default ansatz (Phase 1):** **Vortex line z-translation invariant**

Single vortex profile:
$$\Phi_{vortex}(r_\perp, \phi) = \rho(r_\perp) \cdot e^{in\phi}$$

z BC: $\rho(0) = 0$ (core), $\rho(r_\perp \to \infty) = v$ (vacuum).

EL equation dla ρ(r_⊥):
$$\rho'' + \frac{\rho'}{r_\perp} - \frac{n^2}{r_\perp^2}\rho - \lambda \rho(\rho^2 - v^2) = 0$$

**Reason for vortex default:** U(1) phase symmetry natively gives vortex strings (Nielsen-Olesen analog). σ_ab + RP² hedgehog requires more complex joint config (deferred to alternative if vortex fails).

**DEC budget commit decision rule:**

- IF vortex analysis gives clean 1/L OR log(L) → DEC NOT used (0/1 preserved) → vortex result reported
- IF vortex analysis gives pure exponential (rejecting F-γ-1) → DEC commit (1/1) → switch to point defect (hedgehog + RP²) test for completeness
- IF both vortex AND point defect give pure exponential → F-γ-1 FAIL substantive → HARD HALT scenario

**DEC pre-registered formal: 1/1 budget reserved.**

---

## §10 — Literature checkpoint (informational only)

**Status:** INFORMATIONAL. Literature jest **kontekst**, NIE target.

**Anchors:**
- Nielsen-Olesen 1973 — vortex line in Abelian Higgs (note: gauge model; we have global U(1) — DIFFERENT)
- Vilenkin-Shellard 1994 "Cosmic Strings and Other Topological Defects" — vortex tension + interaction (note: Abelian Higgs framework)
- Manton-Sutcliffe 2004 "Topological Solitons" — solitons w 2D/3D, Skyrme models
- Goldstone-Salam-Weinberg 1962 — Goldstone theorem
- Peskin-Schroeder ch.4, ch.20 — SSB, Goldstone mode

**Key methodological observation:** Most cosmic string literature uses **Abelian Higgs** (local gauge) → all modes massive → exponential interaction at long range. **This differs from TGP S05 pre-registered GLOBAL U(1)** which keeps Goldstone massless → 1/L Coulomb.

**Anti-BD-drift discipline:** Używamy literature dla **methodology** (Goldstone theorem, Green functions). NIE adopting Abelian Higgs framework as analytical baseline.

---

## §11 — Open questions (do rozwiązania w Phase 1+)

1. Explicit ρ(r_⊥) profile dla TGP vortex (analytical lub numerical)?
2. Vortex core size r_0 = ? (z V''(ρ=0) lub similar)
3. Exact coupling constant C w V_int(L) = -C·v²·n_1n_2/L (Phase 2 derivation)?
4. Subleading exp(-m_σL) corrections magnitude vs leading 1/L (crossover scale)?
5. Z₂ symmetry impact na phase winding (n integer, half-integer, or fractional)?
6. RP² topology impact na defect classification (jeśli switch to point defect DEC)?

Wszystkie pre-registered jako open; addressed w odpowiednich fazach.

---

## §12 — Tests planned dla Phase 1/2/3/4/FINAL (counts only)

| Phase | T_xxx tests | DEC budget | LIT/INVENTORY | Substantive FP |
|-------|-------------|------------|----------------|----------------|
| 1 | 5 (ansatz EL, far-field, mass spectrum, energetic stability, RP² check) | 0 | 1 | 4 |
| 2 | 5 (propagator, two-vortex E(L) integration, far-field expansion, scaling, dimensional check) | 0 | 0 | 5 |
| 3 | 3 (fit R²_power, fit R²_exp, differential test) | 0 | 0 | 3 |
| 4 (conditional) | 3 (self-consistency convergence, native bg form, no exogenous D) | ≤1 | 0 | 3 |
| FINAL | 1 (aggregate verdict) | ≤1 | 0 | 1 |
| **Total** | **17** (γ-1) + **3** (γ-2 conditional) | **≤1** | **1** | **16** (γ-1) + **3** (γ-2 conditional) |

**Substantive FP ratio target:** ≥90% pass dla F-γ-1.
**Hardcoded T_pass=True count target:** 0 (strict cycle 1/2/7).

---

## §13 — Status końcowy Phase 0

- ✅ External inputs (EXT-1..EXT-8) inventoried
- ✅ LOCKED structural axioms declared (AX-S05/Z2/U1/RP2 + AX-DECL-1/2 + AX-CE-STR)
- ✅ Derived outputs (OUT-1..OUT-4) z pre-registered predictions
- ✅ Tautology test passed (no circular reasoning)
- ✅ Falsifiability test passed (all outputs falsifiable; HARD HALT scenario explicit)
- ✅ Anti-BD-drift check passed (native equations FIRST; literature informational only)
- ✅ Independent-path cross-validation declared (Paths A/B/C/D)
- ✅ **Analytical pre-derivation step (§3.6 BINDING)** completed §8 — EXPECTED form 1/L Coulomb dla point defect lub log(L) dla vortex line
- ✅ Ansatz declaration locked (vortex line default; DEC reserved for point defect switch)
- ✅ Test counts pre-registered (17 total γ-1, ≤1 DEC)
- ✅ Literature checkpoint informational (NOT validation target)
- ✅ Open questions identified

**Phase 0 LOCKED 2026-05-22. Ready for Phase 1 authorization.**

**§3.6 BINDING compliance check:** ✅ Analytical pre-derivation explicit (§8.1-§8.9). Heuristic steps honestly declared (§8.9). Phase 1-3 sympy will VERIFY analytical expectations, NIE assume them.

---

**END OF PHASE 0 — Balance sheet + Analytical pre-derivation LOCKED 2026-05-22**
