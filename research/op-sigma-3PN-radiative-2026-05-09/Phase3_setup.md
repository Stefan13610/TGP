---
title: "Phase 3 setup — Higher-PN structure of σ-radiation z corrected ξ_eff"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 3
predecessor:
  - "[[./Phase2_results.md]] (UPGRADED do STRUCTURAL DERIVED post-T3.4-amendment)"
  - "[[../op-T34-normalization-amendment-2026-05-09/Phase_FINAL_close.md]] (17/17 PASS, ξ_eff = 4·G·Φ_0²)"
classification: SETUP_PRE_COMPUTATION
critical_decisions:
  - "Scope redefinition: Phase 3 z handoff §3.2 (higher-PN amplitude/phase) NIE z original §2.1 (2nd-order δΦ from σ source)"
  - "Three-channel decomposition of potential h_TT^σ deviation from h_TT^GR at 1PN+ orders"
  - "Yukawa-suppression concern flag dla separate adversarial audit"
  - "Smoking-gun 2PN observable separation: σ-channel (Path A) vs g_eff M9.1'' channel"
---

# Phase 3 setup — Pre-computation locks

## §0 — Phase 3 scope redefinition

### §0.1 — Scope shift z original cycle plan

**Original cycle README §2.1 Phase 3:** "2nd-order δΦ from σ source + far-field expansion"

To było motywowane Phase 1 hipotezą Route A (σ at 3PN+ as resolution dla R5 risk
linearized scalar-mode). Phase 2 Path A direct calculation (24/24 PASS) wykazała,
że σ-radiation pojawia się na **LEADING quadrupole order**, NIE 3PN+. Phase 1
§5.4 "U² ~ 1%" estimate WAS BASED ON WRONG ASSUMPTION (σ as bilinear (∂Φ)²
rather than direct T^TT source).

**HANDOFF §3.2 redefines Phase 3-FINAL post-Phase-2 + post-T3.4-amendment:**
- Phase 3: Higher-order PN corrections (1PN beyond leading, 2PN amplitude)
- Phase 4: Multi-event polarization tests z LIGO O3 catalog
- Phase 5: M9.1'' 2PN deviation specific prediction dla O5+
- Phase FINAL: full closure z amendment integrated

**Phase 3 nowa scope (this setup):** structurally lock PN expansion of TGP h_TT^σ
through 2PN (next-to-next-to-leading) amplitude order in Path A massless
approximation. Identify channels of potential deviation z h_TT^GR. Demonstrate
deviation = 0 channel-by-channel OR identify specific concrete deviation.

### §0.2 — Connection do M9.1'' 2PN deviation smoking-gun

TGP_FOUNDATIONS.md §3.6 + audyt/T01_LIGO3G_falsifier identyfikują **2PN deviation
~0.02 rad at LIGO O5+** jako M9.1''-specific signature. Source mechanism:

- **g_eff M9.1'' channel:** g_tt = -c²·f(ψ) z f(ψ) = (4-3ψ)/ψ form gives
  explicit `|Δg_tt| = (5/6)·U³` deviation od GR Schwarzschild at 2PN level
  (op-newton-momentum/M9_1_pp_setup §2.5)
- **β_ppE^TGP^(b=-1) = -15/4 ≈ -3.75** (post-G_SPA = 48 LOCK; specific (4-3ψ)/ψ
  form RULED OUT 5.02σ by GWTC-3 RE-RUN)
- **Recovery:** post-emergent-metric framework allows alternate f(ψ) forms w
  parametric family β_ppE^new(c_0, ξ_3, a_2, ...) z zero-β region (Phase 4 LOCK)

**Critical separation:** σ-radiative channel (Path A, this Phase 3) i g_eff
M9.1'' channel (op-newton-momentum / op-emergent-metric) are **independent
contributions** to GW phase. Phase 3 isolates σ-channel structure; M9.1''
channel is separate cycle scope.

## §1 — Three channels of potential deviation

Path A z corrected `ξ_eff = 4·G·Φ_0²` daje at LEADING quadrupole:
```
h_TT^σ,(0) = (c_0·ξ_eff/(8π·Φ_0²·c⁴·r))·d²Q^M_TT/dt² = (2G/(c⁴·r))·d²Q^M_TT/dt² = h_TT^GR,(0)
```
EXACTLY (Phase 2 §2 post-amendment).

**At higher PN orders, three potential channels can produce δh = h_TT^σ - h_TT^GR ≠ 0:**

### Channel A — σ self-coupling at higher order

**Hypothesis:** Path A Lagrangian z OP-7 T3.1 contains σ-σ self-interaction terms
beyond -(ξ_eff/2)·σ·T^TT. These would generate nonlinear σ corrections at higher
PN orders.

**Pre-computation expectation:** Path A Lagrangian is **purely linear** w σ:
```
L_σ = -(1/4)(∂σ)² - (1/2)m_σ²·σ² - (ξ_eff/2)·σ·T^TT
```
EOM: `□σ + m_σ²·σ = -ξ_eff·T^TT` jest **linear w σ**. Source T^TT zależy od matter
(NOT od σ self-coupling). Therefore Channel A contributes **zero** deviation at
all PN orders — by structural Lagrangian property.

**Sympy LOCK:** verify EL equations from L_σ are linear (no nonlinear σ coupling).

### Channel B — Massive σ propagation correction

**Hypothesis:** Massless Green function used w Phase 2 derivation; m_σ²σ² term in
EOM produces dispersion correction at higher PN orders.

**EOM full:** `(□ + m_σ²)σ_ab = -ξ_eff·T_ab^TT`

**Massive Green function (frequency domain, retarded):**
```
G_m(k,ω) = 1 / (k² - ω²/c² + m_σ²c²/ℏ²)
G_0(k,ω) = 1 / (k² - ω²/c²)             [massless]

δG/G_0 ≈ -(m_σ c²/(ℏω))² · (1 + O((kc/ω)²))      [for ω-driven radiation]
```

**Phase shift correction structure** (analog Will 1998 for massive graviton):
```
δφ_disp(f) ~ (D_L · m_σ²c⁴) / (2 · ℏ² · (2πf)²)
           = (D_L · π) / (λ_g² · f)            [λ_g = ℏ/(m_σc) Compton wavelength]
```

**Pre-computation expectation z TGP m_σ scale:**
- Phase 2 §1.1 quotes "M_eff/ω_LIGO ~ 10⁹" (massless decoupling regime)
- Path B audit (closure 2026-04-26): M_eff = √2·m_s ≈ 0.71 meV
- For LIGO O5+ band f ~ 100 Hz: ℏω ~ 4·10⁻¹³ eV
- Ratio m_σc²/ℏω ~ 0.71 meV / 4·10⁻¹³ eV ~ 1.8·10⁹ (HEAVY regime, NOT light)
- **Yukawa suppression** at distance D ~ Gpc: e^(-m_σc·D/ℏ) ~ e^(-3·10²¹) — astronomically suppressed

**⚠ HONEST AUDIT FLAG (Phase 3 §1B-AUDIT):** Phase 2 Path A derivation used
massless approximation; with m_σ ≈ 0.71 meV ≫ ℏω_LIGO, the massless propagator
should NOT be valid for σ-mediated radiation at LIGO scales. Either:
- (i) m_σ is effectively renormalized to zero in IR (some EFT mechanism)
- (ii) σ_ab is composite (Path B audit) z light constituents — radiation through different channel
- (iii) "h_TT^σ amplitude at observer" reaches observer through emergent-metric mediation, NOT direct σ propagation
- (iv) Phase 2 Path A formula is structurally correct **as effective coupling** but interpretation as "σ wave radiation z observer" requires re-examination

**This is a separate audit item** (NOT resolved in Phase 3 single-session).
Phase 3 proceeds w Phase 2's effective massless framework as established
methodology, AND flags this concern explicitly w results dla adversarial
follow-up cycle (proposed: `op-sigma-yukawa-audit-2026-05-XX`).

### Channel C — Emergent-metric C(ψ) Taylor expansion at non-vacuum

**Hypothesis:** Phase 4 emergent-metric ansatz `g_eff^ij = δ^ij·B(ψ) + σ^ij·C(ψ)/(Φ_0²·c²)`
uses `C(ψ_0) = c_0` at vacuum. During inspiral, ψ varies near source; observer
in vacuum sees ψ ≈ ψ_0 → C(ψ) ≈ c_0.

**Pre-computation expectation:**
- At observer (vacuum infinity): δψ → 0 → C(ψ_obs) = C(ψ_0) = c_0 EXACT
- Higher PN corrections z C'(ψ_0)·δψ vanish at far-field
- TGP wave amplitude at observer is determined by emitted-field (matter-region)
  σ structure × free vacuum propagation × c_0 factor at observation point
- **No 1PN or 2PN amplitude correction** z this channel structurally

**Sympy LOCK:** verify δC(ψ_obs)/δψ → 0 jest equivalent to vacuum BC `δψ(∞) = 0`,
giving observer-side C = c_0 exact at all orders.

### Channel D — Higher multipoles in T_ab^TT

**Hypothesis:** GR at 2PN amplitude includes mass octupole `d³M_ijk/dt³`,
current quadrupole `d²S_ij/dt²`, and other higher multipoles. Path A z
T_ab^TT linear coupling inherits ALL of these structurally.

**Pre-computation expectation:**
- T_ab^TT for matter point particles contains full multipole hierarchy:
  T_ij^TT = TT[Σ m_a v_a^i v_a^j δ³ - O(v²)·multipoles]
- σ_ab^far via retarded Green function inherits this hierarchy:
  σ_ab^far(observer) = ∫ T_ab^TT(retarded)/r d³y
- **Multipole-by-multipole h_TT^σ matches h_TT^GR exactly** in massless Path A

**Sympy LOCK:** for octupole and current quadrupole moment formulas, verify
the matching condition `c_0·ξ_eff = 16π·G·Φ_0²` (z amendment) reproduces GR
coefficient at each multipole order (NOT just leading mass quadrupole).

## §2 — PN-counting structure dla σ-channel

| PN order | GR contribution | TGP σ-channel (Path A massless) | Deviation source |
|---|---|---|---|
| Leading (Newtonian) | mass quadrupole d²Q^M/dt² | matches GR EXACT post-amendment | none (Phase 2) |
| 1PN amplitude | velocity corrections to Q^M | inherited from T_ab^TT | none (Channel D) |
| 1.5PN amplitude | tail terms (back-scattering) | requires hereditary integral analysis | DEFERRED |
| 2PN amplitude | mass octupole + current quadrupole | inherited from T_ab^TT | none (Channel D) |
| 2.5PN amplitude | back-scattering at 2.5PN | hereditary, deferred | DEFERRED |
| Phase (orbital) | quadrupole formula loss → chirp | identical via Path A T_ab^TT | none |
| Channel B dispersion | n/a (graviton massless) | δφ ~ (m_σ c²)²/(ℏω)² | **AUDIT FLAG** |

**Conclusion z table:** at NON-hereditary amplitude orders (leading, 1PN, 2PN),
TGP σ-channel structurally matches GR via Path A T_ab^TT inheritance — IF
massless approximation valid (Channel B audit flag).

## §3 — Falsifier dla Phase 3

If Phase 3 sympy demonstrates:
1. ✓ Channel A LOCK: σ self-coupling = 0 at all orders
2. ✓ Channel C LOCK: C(ψ_obs) = c_0 exact at vacuum infinity
3. ✓ Channel D LOCK: matching condition `c_0·ξ_eff = 16π·G·Φ_0²` reproduces GR
   coefficient at multipole-by-multipole level (mass octupole + current quadrupole)
4. ⚠ Channel B FLAG: dispersion correction structure derived; magnitude depends
   on resolved m_σ scale (audit flag preserved as separate cycle)

**Then Phase 3 classifies STRUCTURAL DERIVED z honest audit-flag** dla σ-channel
PN structure through 2PN amplitude (NON-hereditary). 1.5PN/2.5PN hereditary tail
+ Channel B Yukawa audit deferred do Phase 4-5 / dedicated audit cycles.

If sympy reveals UNEXPECTED deviation structure (Channel A nonlinear coupling
discovered, OR multipole matching fails, OR Channel C non-trivial vacuum
behavior):
- Phase 3 classifies STRUCTURAL_CONDITIONAL z explicit gap identification
- Trigger sub-cycle audit similar to T3.4 amendment cycle pattern
- Adversarial verification commitment per CALIBRATION_PROTOCOL §4.3

## §4 — Anti-pattern compliance

- ✓ Pre-declared four-channel decomposition (NIE multi-candidate fit)
- ✓ Specific structural expectations PER channel (zero / non-zero / audit-flag)
- ✓ Matching condition reuse (NIE re-derivation; inherits T3.4 amendment LOCK)
- ✓ Honest scope: NON-hereditary amplitude only (1.5PN/2.5PN tail deferred)
- ✓ Honest acknowledgment: Channel B Yukawa concern flagged, NIE swept under rug
- ✓ Adversarial verification commitment dla Channel B audit follow-up

## §5 — Time budget

**Single-session realistic scope:**
- Phase 3 setup (this file)
- Phase 3 sympy: 4 channels × ~4-5 LOCKS each = ~16-20 PASS targeted
- Phase 3 results: structural derivation + audit-flag honest reporting

**Multi-session deferred:**
- Phase 3.5: 1.5PN + 2.5PN hereditary tail terms (technical PN expansion w Path A)
- Phase 4: Multi-event LIGO O3 catalog test (data analysis cycle)
- Phase 5: M9.1''-recovery 2PN signature specific prediction (joint z emergent-metric Phase 4 Path 2)
- Phase FINAL: full closure z amendment + audit flags resolution
- Adversarial: `op-sigma-yukawa-audit` (resolves Channel B m_σ vs ω_LIGO)

## §6 — Cross-references

- [[./README.md]] — cycle setup
- [[./Phase1_results.md]] — Phase 1 foundation (11/11 PASS)
- [[./Phase2_results.md]] — Phase 2 Path A direct (24/24 PASS, UPGRADED post-amendment)
- [[./Phase2_setup.md]] — Phase 2 setup (PN-counting + dimensional analysis)
- [[../op-T34-normalization-amendment-2026-05-09/Phase_FINAL_close.md]] — ξ_eff amendment LOCK
- [[../op7/OP7_T3_results.md]] — T3.4 amendment notice + Path A Lagrangian (T3.1)
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] — Path B audit, M_eff = √2·m_s ≈ 0.71 meV
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — emergent-metric ansatz C(ψ)
- [[../../audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] — 2PN deviation z g_eff M9.1'' (separate channel)
- [[../../TGP_FOUNDATIONS.md]] §3.6.10.6 — post-T3.4-amendment status
