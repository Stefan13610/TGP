---
title: "Phase 2 results — S07 alternative f(ψ): symbolic Bayesian α-mapping + family distinguishability + H1a TENTATIVE verdict draft"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 2
predecessor: "[[./Phase1_results.md]]"
sympy_total: "15/15 PASS (100%); 12 FP (80.0%) + 3 LIT (20.0%) + 2 DEC separate; 0 hardcoded"
phase2_outcome: "α posterior MAP at α_ML ≈ 0 within recovery [-0.832, 0.832] under GWTC-3 1σ; LIGO-O5 A+ ~2027 narrows precision ×3.13 (σ_α^O5 = 80/301 ≈ 0.266); BH5 QNM + ε.1 photon ring families distinguishable via d²f/dψ²(ψ_0) = {0, 2β_q, α²}; cross-cycle Δe_2_native(α) = α/3 z M9.1'' anchor consistency exact; constraint -4ξ_3 + 4 - a_3/8 + 4/3 = α/3 z 1-param {ξ_3, a_3} family per α"
verdict_draft: "H1a TENTATIVE — recovery successful under current GWTC-3 1σ; pending observational LIGO-O5 A+ ~2027 verification + Phase 3 (deferred) BH5/ε.1 numerical for pre-observational family discrimination if needed"
---

# Phase 2 results — S07-reset alternative f(ψ)

## §0 — Headline

```
████████████████████████████████████████████████████████████████████
█  op-S07-reset-alternative-f-psi-2026-05-11  Phase 2              █
█  15/15 PASS — 12 FP (80.0%) / 3 LIT / 0 hardcoded                █
█                                                                  █
█  B.1 Bayesian α posterior:                                       █
█    α_ML(GWTC-3) = 0    (typical ToGR null)                       █
█    σ_α^GWTC-3 = 104/125 ≈ 0.832    (1σ)                          █
█    σ_α^O5     = 80/301  ≈ 0.266    (×3.13 improvement)           █
█                                                                  █
█  B.2 Family distinguishability via d²f/dψ²(ψ_0):                 █
█    polynomial = 0    quadratic = 2β_q    transcendental = α²     █
█    BH5 QNM marker + ε.1 photon ring quad channel STRUCTURALLY    █
█    DISTINCT — Phase 3 numerical evaluation deferred              █
█                                                                  █
█  B.3 Δe_2_native(α) = α/3 derived; M9.1'' α=-4 → -4/3 EXACT      █
█    constraint -4ξ_3 + 4 - a_3/8 + 4/3 = α/3 z c_0·κ_σ=4/3 LOCK   █
█                                                                  █
█  H1a TENTATIVE — pending LIGO-O5 A+ ~2027 + (optional) Phase 3   █
████████████████████████████████████████████████████████████████████
```

## §1 — Test results table

| Test | Klasa | Status | Substance |
|---|---|---|---|
| T1 | FP | PASS | Jacobian α=(16/15)·β_ppE; dα/dβ=16/15 const; product=1 |
| T2 | FP | PASS | α_ML at GWTC-3 β_ML=0: α_ML=0 trivial GR |
| T3 | FP | PASS | α 1σ bound (16/15)·0.78 = 104/125 ≈ 0.832 reproduces Phase 1 T4 |
| T4 | FP | PASS | LIGO-O5 A+ projection σ_α^O5 = 80/301 ≈ 0.266 (×3.13 vs GWTC-3) |
| T5 | FP | PASS | Polynomial d²f/dψ²(ψ_0) = 0 (only α at all PN) |
| T6 | FP | PASS | Quadratic d²f/dψ²(ψ_0) = 2β_quad (NEW family parameter) |
| T7 | FP | PASS | Transcendental d²f/dψ²(ψ_0) = α² (nonlinear coupling) |
| T8 | FP | PASS | BH5 QNM marker δω/ω_GR ∝ d²f/dψ²(ψ_0): families distinguishable |
| T9 | FP | PASS | ε.1 photon ring: linear ch shared (α/3); quad ch poly=0/quad=β_q/9/exp=α²/18 |
| T10 | FP | PASS | Cross-cycle Δe_2_native(α) = α/3; M9.1'' α=-4 → -4/3 EXACT anchor consistency |
| T11 | FP | PASS | Constraint -4ξ_3 + 4 - a_3/8 + 4/3 = α/3 → ξ_3 = -a_3/32 - α/12 + 4/3 (1-param family per α) |
| T12 | FP | PASS | GR-limit α=0: f_poly = f_exp = 1, Δe_2 = 0, β_ppE = 0 (sanity check) |
| T13 | LIT | PASS | GWTC-3 |β_ppE| ≤ 0.78 1σ; M9.1'' rejected at 5.02σ (Abbott+2021) |
| T14 | LIT | PASS | LIGO-O5 A+ ~2027 SNR=15.05σ on M9.1'' (PR-002 inheritance) |
| T15 | LIT | PASS | ngEHT photon ring +14.6% M9.1'' anchor; BH5 Cosmic Explorer ~2030 |

**Declarative (separate count, NIE counted w PASS):**
- T16: Anti-Lakatos LOCKED PR-010 (brak H1c/H1d)
- T17: Three-layer L1/L2/L3 presentation MANDATORY (this section §3)

## §2 — Six P-requirements (Phase 2 final contribution)

| P | Resolution post-Phase-2 |
|---|---|
| P1 | ✅ RESOLVED — 3 families enumerated (Phase 1) + d²f/dψ²(ψ_0) distinguishability marker (Phase 2 T5+T6+T7+T8+T9) |
| P2 | ✅ RESOLVED — β_ppE^poly(α) = (15/16)·α (Phase 1) + α MAP estimate via Jacobian (Phase 2 T1+T2) |
| P3 | ✅ RESOLVED — GR-limit f(ψ_0)=1 dla wszystkich families (Phase 1 + Phase 2 T12 sanity) |
| P4 | ✅ **CLOSED** — symbolic Bayesian α posterior MAP estimation + 1σ/5σ bounds derived (Phase 2 T1-T4); pełny MCMC pipeline = separate cycle, out of substance protocol |
| P5 | ✅ RESOLVED — c_0·κ_σ=4/3 LOCK preserved + Δe_2_native(α) consistency (Phase 2 T10+T11) |
| P6 | ✅ RESOLVED — S05 preserved (Phase 1 T11; Phase 2 T17 declarative) |

**Wszystkie 6 P-requirements RESOLVED.** Phase FINAL closure ready.

## §3 — Three-layer L1/L2/L3 presentation (MANDATORY per PPN_AS_PROJECTION §3.1)

### §3.1 — L1 (Native predictions, PRIMARY)

**Native observable (z LIGO-3G-native A− methodology inherited):**

```
Δφ(f) = -(15/4) · Δe_2_native / (M · (πMf)^(1/3))   [radians]
```

gdzie `Δe_2_native = -4·ξ_3 + 4 - a_3/8 + c_0·κ_σ` (LIGO-3G-native cycle).

**Native Taylor coefs constrained przez Phase 2 (z polynomial f(ψ) family analysis):**

| Native coef | Constraint po Phase 2 | Klasyfikacja per PPN_AS_PROJECTION §3.3 |
|---|---|---|
| `α` (S07 polynomial coef) | α ∈ [-0.832, 0.832] z GWTC-3 1σ; α_ML ≈ 0 typical | **Independent** (Phase 2 directly constrains) |
| `c_0` | inherited LOCK c_0 = 4π z `op-c0-derivation` | Path 2 anchor inheritance (heuristic) |
| `κ_σ` | inherited LOCK κ_σ = 1/(3π) z `op-kappa-sigma-2body-PN` | Path 2 anchor inheritance (heuristic) |
| `c_0 · κ_σ` | **EXACT** = 4/3 (joint product LOCK z emergent-metric Phase 4) | Path 2 anchor LOCK (clean π cancellation) |
| `ξ_3` | constraint ξ_3 = -a_3/32 - α/12 + 4/3 (1-param family per α) | **Coupled** (depends on α + a_3) |
| `a_3` | free (in 1-param family with ξ_3) | **Free** (deferred to higher-PN cycle) |
| `β_quad` | free (quadratic family parameter; constrained by BH5/ε.1 jeśli detection) | **Free** (Phase 3 deferred) |
| `α₁,₂,₃, ζ_i` | ≡ 0 strukturalnie z S05 + Lorentz-invariance | **Forced ≡ 0** (declared, NIE liczone) |

**L1 native parameter audit (per §3.3 mandatory):**

```
Independent Taylor coefs constrained by this cycle: {α (Phase 2 explicit)}
Coupled coefs: {ξ_3, a_3} (constraint -4ξ_3 + 4 - a_3/8 + 4/3 = α/3; 1-param family)
Inherited Path 2 anchors (heuristic): {c_0 = 4π, κ_σ = 1/(3π)} z explicit citation
Free coefs (deferred to Phase 3 / dedicated cycles): {β_quad, a_3, higher-PN coefs}
Forced coefs (substrate symmetry): {ζ_i, α_i, b_ppE-index} ≡ 0 strukturalnie
Total Phase 2 native parameter count: 1 independent (α) + 1 coupled-pair (ξ_3, a_3) + 2 anchored
```

### §3.2 — L2 (PPN/ppE projection consistency map)

| ppE param | Native expression | Phase 2 derived value |
|---|---|---|
| `β_ppE^TGP` (b=-1, 2.5PN phase) | `(45/16) · Δe_2_native` = `(45/16) · (α/3)` = `(15/16) · α` | **Phase 2 T10 verified symbolic** |
| `α_ppE^TGP` (amplitude) | `h_TT^σ / h_TT^GR` ratio (T3.4 amendment cycle) | = 1 EXACT (inherited; Phase 2 nie constraints) |
| `b_ppE^TGP` (PN order index) | `b = -1` (2.5PN phase position) | FORCED z Φ-EOM mass scale (substrate symmetry) |

**Anti-pattern avoidance per §3.1 explicit warning:** Phase 2 results.md NIE prezentuje
"TGP predicts β_ppE = (15/16)·α" jako primary output. Primary L1 = native phase model
Δφ(f) z constraints na {α, ξ_3, a_3} z 1-param coupling. β_ppE = L2 projection
consistency map (this section).

### §3.3 — L3 (Falsification map)

| Observational bound | Constrains native coef | Window | Status |
|---|---|---|---|
| GWTC-3 |β_ppE| ≤ 0.78 (1σ; Abbott+2021 90 BBH) | α via Jacobian | α ∈ [-0.832, 0.832] | **PASSED** if α_ML ≈ 0 (current data consistent z null) |
| GWTC-3 M9.1'' β=-15/4 rejected at 5.02σ | α = -4 specific point | Excluded | **CONSISTENT** — M9.1'' specific point excluded; family neighborhood preserved |
| LIGO-O5 A+ ~2027 SNR=15.05σ on M9.1'' (PR-002) | σ_α improvement ×3.13 | σ_α^O5 ≈ 0.266 | **PENDING OBSERVATIONAL** (~2027) |
| LIGO-O5 single-event 5σ (projection) | |α|_5σ^O5 ≈ 1.33 | Wider than current recovery; CANNOT exclude region single-event | **NEED COMBINED POSTERIOR** lub Phase 3 BH5/ε.1 cross-channel |
| ngEHT photon ring +14.6% M9.1'' (op-eht) | f(ψ_photon) shift family-dependent | Linear ch α/3 + quad ch family-distinct | **DATA POINT** at α=-4 anchor; family-extrapolation deferred Phase 3 |
| BH5 QNM Cosmic Explorer ~2030 | d²f/dψ²(ψ_0) marker | poly=0 / quad=2β_q / exp=α² | **FUTURE TEST** distinguishes families pre-observationally |

**Falsification propagation:** Phase 2 demonstrates że falsyfikacja w L2 (`β_ppE` outside
GWTC-3 window) propaguje do L1 jako constraint na native coef α (NIE jako "TGP
parameter excluded"). Recovery region α ∈ [-0.832, 0.832] = native coef space, NIE PPN
parameter space.

## §4 — Substantive findings z Phase 2 detail

### §4.1 — B.1 Bayesian α posterior (T1-T4)

**Symbolic Jacobian transformation:** dla linear scaling Phase 1 `β_ppE = (15/16)·α`,
inverse `α = (16/15)·β_ppE` jest LINEAR z constant Jacobian `dα/dβ = 16/15`. Posterior
transformation `p_α(α) = p_β((15/16)·α) · (15/16)`. **Bayesian rigor:** linear scaling
ensures posterior shape preserved (gaussian → gaussian, etc.).

**MAP estimate:** dla GWTC-3 ToGR null result `β_ML ≈ 0` (Abbott+2021 reports posterior
consistent z 0), **α_ML ≈ 0** (trivial GR recovery). To jest **central result** Phase 2:
**no detection of TGP S07 polynomial deviation w current GWTC-3 data**, ale finite recovery
region exists `α ∈ [-0.832, 0.832]` — neither H1a confirmed, nor H1b triggered.

**LIGO-O5 A+ projection (PR-002 inheritance):** SNR_M911 = 15.05σ na M9.1'' β = -15/4
implies `σ_β^O5 = (15/4) / 15.05 = 75/301 ≈ 0.249`. Through Jacobian:
`σ_α^O5 = (16/15) · (75/301) = 80/301 ≈ 0.266`. **Improvement factor ×3.13** vs current
GWTC-3 1σ.

**5σ exclusion bound (single-event LIGO-O5):** `|α|_5σ^O5 = 5 · 80/301 ≈ 1.33`. Current
recovery region [-0.832, 0.832] ⊂ [-1.33, 1.33] **single-event O5 5σ window**. Implication:
single-event LIGO-O5 detection nie wyklucza całej recovery region z 5σ — wymaga
**combined posterior** z multiple events lub **cross-channel test** (BH5/ε.1) dla
discrimination.

### §4.2 — B.2 Higher-PN family distinguishability (T5-T9)

**Drugi derivative `d²f/dψ²(ψ_0)` jako family marker:**

| Family | f(ψ) | d²f/dψ²(ψ_0) | Phase 2 marker rola |
|---|---|---|---|
| polynomial | 1 + α(ψ-ψ_0) | **0** | flat curvature; tylko α na wszystkich PN orders |
| quadratic | 1 + α(ψ-ψ_0) + β_q(ψ-ψ_0)² | **2β_q** | new free parameter β_q; distinguishable jeśli β_q ≠ 0 |
| transcendental | exp(α(ψ-ψ_0)) | **α²** | nonlinear w α; distinguishable via α² coupling |

**BH5 QNM marker (T8):** Leading-order δω_QNM/ω_QNM^GR proportional do `d²f/dψ²(ψ_0)`
local kinetic curvature. Trzy struktury (0 / 2β_q / α²) DISTINCT. Cosmic Explorer ~2030
z high-SNR ringdown może pre-observationally discriminate families.

**ε.1 photon ring (T9):** Wszystkie families share linear channel `α · Δψ_ph = α/3`
(Δψ_ph = 1/3 dla M9.1'' canonical ψ_0=1 → ψ_photon=4/3). Quadratic channel distinguishes:

```
poly:  0
quad:  β_q · (Δψ_ph)² = β_q / 9
exp:   α² · (Δψ_ph)² / 2 = α² / 18
```

**Punkt anchor M9.1'' = +14.6% b_crit** w op-eht cycle to data point at α=-4 specific;
family-dependent extrapolation w recovery region deferred do dedicated `op-eht-S07-followup`
cycle (Phase 3 lub separate).

### §4.3 — B.3 Cross-cycle Δe_2_native consistency (T10-T11)

**Symbolic derivation Δe_2_native(α):** Z dwóch independent LOCKs:
- Phase 1 Phase residual cycle: `β_ppE = (45/16) · Δe_2_native`
- Phase 1 polynomial family: `β_ppE^poly = (15/16) · α`

→ `Δe_2_native(α) = (15/16) · α / (45/16) = α/3` **EXACT symbolic**.

**M9.1'' anchor consistency (T10):** Substytucja α=-4 → `Δe_2(α=-4) = -4/3`. **Matches
M9.1'' Path 2 anchor** Δe_2 = -4/3 z LIGO-3G-native A− cycle PR-002. **Cross-cycle
consistency confirmed.**

**Constraint na {ξ_3, a_3} z c_0·κ_σ=4/3 LOCK (T11):**

```
Δe_2_native = -4·ξ_3 + 4 - a_3/8 + c_0·κ_σ
            = -4·ξ_3 + 4 - a_3/8 + 4/3   (Path 2 anchor LOCK)

α/3 = -4·ξ_3(α) + 4 - a_3(α)/8 + 4/3

→ ξ_3(α, a_3) = -a_3/32 - α/12 + 4/3
```

**Sanity checks:**
- α=0 (trivial GR): `ξ_3 = 4/3 - a_3/32` (1-param free)
- α=-4 (M9.1''): `ξ_3 = 5/3 - a_3/32` (M9.1'' specific {ξ_3^M911, a_3^M911} satisfy this)

**Implication:** dla każdego α ∈ recovery region, native parameter freedom = 1 free
(a_3 lub equivalently ξ_3 z constraint coupling). Phase 2 nie constrains a_3 directly;
deferred do higher-PN cycle (Phase 3 lub dedicated).

## §5 — H1a/H1b verdict draft (Phase FINAL setup)

### §5.1 — Anti-Lakatos compliance check

Per PR-010 LOCKED immutable:
- ✅ recovery_scope α ∈ [-0.832, 0.832] preserved (Phase 1 + Phase 2 reproduced)
- ✅ GR-limit f(ψ_0)=1 mandatory (Phase 2 T12 sanity verified)
- ✅ S05 preserved (Phase 2 T17 DEC; Phase 1 T11 sympy)
- ✅ NIE wprowadzono H1c/H1d backstop (anti-Lakatos)
- ✅ NIE post-hoc tuning specific f(ψ) form (3 families pre-declared Phase 1)
- ✅ Brak BD-drift (ASK-RULE Triggers A-D PASS w Phase2_setup §0.1)

### §5.2 — Verdict matrix (Phase 2 evidence applied)

| Scenario | α_ML estimate | Family discriminability | Verdict |
|---|---|---|---|
| α_ML ≈ 0 z GWTC-3 1σ; quad/trans non-degenerate at higher PN | **WITHIN recovery** | **YES (structural)** | **H1a TENTATIVE** |
| α_ML w recovery, family universally degenerate | within recovery | NO | H1a partial |
| α_ML excluded z GWTC-3 5σ | outside recovery | n/a | H1b |

**Phase 2 evidence:**
- α_ML ≈ 0 z current GWTC-3 ToGR null (T2 derived); within recovery [-0.832, 0.832] ✓
- Family discriminability via d²f/dψ²(ψ_0) STRUCTURAL marker (T5+T6+T7+T8+T9) ✓
  - Numerical evaluation per family (BH5 QNM ω shifts, ε.1 photon ring b_crit) deferred Phase 3
- Cross-cycle Δe_2_native(α) = α/3 z c_0·κ_σ=4/3 LOCK consistency ✓

**Phase 2 verdict draft: H1a TENTATIVE.**

### §5.3 — Closure path options

**Opcja A (recommended):** Phase FINAL closure z claim_status **A− (STRUCTURAL_DERIVED_NATIVE
z L2 not-fully-attempted)** pending:
- Observational LIGO-O5 A+ ~2027 single-event verification (PR-002 inheritance LOCKED)
- LIGO-O5 combined posterior z multiple BBH events (precision ×3.13 → potentially tighter
  recovery exclusion)
- Phase 3 numerical BH5/ε.1 evaluation per family (jeśli pre-observational discrimination
  needed; OPTIONAL)

**Opcja B:** Phase 3 BEFORE closure — numerical BH5 QNM ω(α, family) + ε.1 b_crit(α, family)
evaluation per recovery region. Estymata 2-3 sesji additional. Closure z bardziej
substantywnym verdict (e.g., "polynomial family preferred z BH5 ringdown null detection").

**Opcja C:** Phase 2 NIE close cycle — kontynuuj Phase 3 jako multi-session deferred work
deferred do future session decision.

**Recommendation:** **Opcja A** — Phase 2 established H1a TENTATIVE z STRUCTURAL family
distinguishability marker. Numerical refinement Phase 3 nie zmieni binary H1a/H1b decision,
tylko narrows posterior. Anti-Lakatos protocol favors timely close-pending-data over
extended pre-observational refinement.

## §6 — Phase 2 close gate

**Phase 2 GATE: ✅ OPEN — 15/15 PASS z substantive symbolic Bayesian α-mapping +
family distinguishability marker + cross-cycle consistency derivations.**

**Substance metrics:**
- 15/15 sympy PASS (100%)
- 12 FP (80.0%) — exceeds 75% binding target per AUDIT_2026-05-11 substance protocol
- 3 LIT (20.0%) — explicit observational anchors (GWTC-3, PR-002, op-eht)
- 0 hardcoded `T_pass = True` — preserved from Phase 1
- 100% non-trivial — every test verified explicit symbolic statement

**P-requirements final status:** **6/6 RESOLVED** post-Phase-2 (P4 closed via symbolic
Bayesian Jacobian).

**Anti-Lakatos PR-010 compliance:** ✅ all 6 sub-checks PASS (§5.1).

**Three-layer L1/L2/L3:** ✅ explicit (§3.1+§3.2+§3.3).

**M9.1'' Path 2 anchor reframe per M9_RESTRUCTURE §3:** ✅ kategoria (b) inheritance
explicit cited (§4.3 + §3.1 native coef table).

## §7 — Cumulative cycle metrics

| Phase | Sympy | FP | LIT | DEC | Notable |
|---|---|---|---|---|---|
| Phase 1 | 12/12 | 10 (83.3%) | 2 | 2 | β_ppE^poly(α) = (15/16)·α LINEAR SCALING; recovery [-0.832, 0.832] |
| Phase 2 | 15/15 | 12 (80.0%) | 3 | 2 | symbolic Bayesian α-mapping; family distinguishability marker; cross-cycle Δe_2(α)=α/3 |
| **Cumulative** | **27/27** | **22 (81.5%)** | **5 (18.5%)** | **4 separate** | **0 hardcoded; 100% non-trivial** |

## §8 — Status & next session

**Phase 2 closed analytical work.** Phase FINAL `Phase_FINAL_close.md` w **separate session**
(per workflow rule "nie próbuj zamknąć cycle w jednej sesji jeśli wymaga substantywnej
numerical work") z:

1. Closure ceremony per LIGO-3G-native A− template
2. Verdict claim_status = **A−** (Opcja A recommended)
3. Cross-cycle update: PREDICTIONS_REGISTRY entry (M9.1'' family-recovery preserved); PR-010
   status LOCKED-PHASE-2-COMPLETE; future-test annotation LIGO-O5 A+ ~2027 + Cosmic
   Explorer ~2030 BH5
4. STATE.md update: WIP slot 1/5 zwolniony; remaining slot 2/5 = inflation Phase 2 active
5. RETROFIT_SPRINT_2026-05-13_summary.md append jeśli substantive (Phase 2 derivations
   substantive ✓ → append §5b.1.B Phase 2 close)

**Estymata Phase FINAL closure:** 0.5-1 sesja (closure ceremony + cross-cycle propagation).

**Optional Phase 3 (if user authorizes):** numerical BH5/ε.1 per family — 2-3 sesji,
lower-priority post H1a TENTATIVE established.

---

**Phase 2 close.** Symbolic Bayesian α-mapping + higher-PN family distinguishability marker
+ cross-cycle Δe_2_native(α) consistency derived first-principles. **H1a TENTATIVE verdict
ready dla Phase FINAL closure A− pending observational LIGO-O5 A+ ~2027.**
