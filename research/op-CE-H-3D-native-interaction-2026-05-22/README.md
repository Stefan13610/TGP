---
title: "op-CE-H-3D-native-interaction-2026-05-22 — Poziom γ-1: native 3D U(1) interaction F-γ-1 CRUCIAL TEST"
type: research_cycle
status: PRE_REGISTERED_LOCKED
folder_status: active
pre_registration_date: 2026-05-22
parent_concept_paper: meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md
parent_cycle: op-CE-H-two-particle-equilibrium-2026-05-21 (A- conditional Poziom β CLOSED)
parent_audit: op-R2-integration-audit-CE-H-FFS-2026-05-22 (R2_PASS CLOSED 2026-05-22)
hypothesis_codes:
  - CE-H (Cosmic Equilibrium Hypothesis) — structural feature TGP per R3 3/3
  - F-γ-1 CRUCIAL TEST (LOCKED 2026-05-21 PENDING_POZIOM_GAMMA; activated 2026-05-22)
test_level: NATIVE_3D (full S05+Z₂+U(1)+RP² + 3D propagator)
test_scope: defect-defect interaction at large separation w 3D
methodology: native equations FIRST + ANALYTICAL PRE-DERIVATION (CALIBRATION §3.6 BINDING 2026-05-22)
claim_status_target: A | A- (depending on F-γ-1 verdict)
authorization_chain:
  - "2026-05-22: User authorized Path B (Poziom γ-1) per sequence commitment A→B→C"
  - "Post R2_PASS verdict (2026-05-22): methodology BINDING propagated; first cycle post-§3.6 BINDING"
discipline:
  - anti-Lakatos LOCKED
  - strict cycle 1/2/7 (0 hardcoded FP T_pass=True)
  - max 1 DEC budget per cycle
  - R1+R2+R3 two-tier discipline BINDING (CALIBRATION §6)
  - pre-rejestracja PRZED any sympy
  - **§3.6 BINDING analytical pre-derivation step** (Phase 0 MUST include)
falsifiers_pre_registered: F-γ-1 (CRUCIAL TEST, see §3); F-γ-2 (SECONDARY, see §3); F-γ-3/4 NOT activated (cosmological scope = γ-3+)
related:
  - "[[../op-CE-H-two-particle-equilibrium-2026-05-21/]] (Poziom β parent)"
  - "[[../op-R2-integration-audit-CE-H-FFS-2026-05-22/]] (R2 audit parent)"
  - "[[../op-FFS-quark-object-2026-05-20/]] (FFS structural baseline)"
  - "[[../../meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md]] (Poziom α concept paper)"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]] §6 + §3.6 (BINDING 2026-05-22)"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §7 (F-γ-1..4 PENDING entries)"
  - "[[../../meta/TGP_W_Z_THEORETICAL_LIMIT.md]] §6.5 (path η cosmology toy extension)"
---

# op-CE-H-3D-native-interaction — Poziom γ-1 BINDING contract

**Pre-registration date:** 2026-05-22
**Status:** LOCKED — żadnych modyfikacji §3 falsifiers / §4 phase plan / §5 ansatz declarations ex post bez HALT-B.

---

## §0 — Origin i scope

Niniejszy cykl jest **Poziom γ-1** roadmapy zdefiniowanej w:
- Concept paper Poziom α: [[../../meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md]] (LOCKED 2026-05-21)
- CE-H Poziom β Phase FINAL §6: F-γ-1 pre-registered LOCKED 2026-05-21 PENDING_POZIOM_GAMMA
- R2 audit verdict 2026-05-22: PENDING_POZIOM_GAMMA → ACTIVATED 2026-05-22

**Cel cyklu:** Test CRUCIAL F-γ-1 — czy w pełnym 3D TGP (S05+Z₂+U(1)+RP² + 3D propagator) dwa defekty (hedgehog lub vortex) mają **native long-range interaction** (power-law lub logarithmic) jako wynik z TGP-native Lagrangianu, **bez** exogenous D/L^α (jak w 1D Z2 Phase 1b/2 toy)?

**Pre-registered HARD HALT scenario (per concept paper §10.4 + Phase FINAL §6.2):**
Jeśli native 3D też daje czyste exponential exp(-mL) bez power-law/log tail → CE-H bg form NIE jest natively native nawet w 3D → **fundamentalny redesign required**. To **anti-Lakatos honest scenario**, NIE failure to be rescued.

---

## §1 — Pre-registered falsifiers (LOCKED 2026-05-22)

### F-γ-1 — Native 3D long-range interaction (CRUCIAL TEST)

**Pre-registracja:** Dla dwóch defektów (hedgehog lub vortex) w 3D z Lagrangianem S05+Z₂+U(1)+RP², native interaction V_int(L) przy L → ∞ MUSI mieć formę:

- **Power-law:** V_int(L) ~ C/L^β z β finite (NIE β → ∞ exponential limit)
- **Logarithmic:** V_int(L) ~ C·log(L/L_0)
- **Coulomb-like:** V_int(L) ~ -q²/(4πL) (specific case β=1)

**OR exponentially-modified power-law:** V_int(L) ~ C·exp(-mL)/L^β z β finite (Yukawa-like z power-law modulation, NIE pure exponential)

**FAIL criterion:** V_int(L) ~ C·exp(-mL) jako PURE exponential bez any power-law modulation w 3D context.

**Pre-registered tolerance (LOCKED 2026-05-21):**
- Power-law: R² > 0.99 dla power-law/log fit AT large L (L > 3/m_kink)
- Pure exponential clearly distinguishable: R²_exp > R²_power gdyby też było prawdziwie czyste exponential

**Severity:** STRUCTURAL — fail → CE-H bg required exogenously even w 3D → HARD HALT scenario, fundamentalny redesign required.

### F-γ-2 — Self-consistency closure z native bg (SECONDARY)

**Pre-registracja:** Jeśli F-γ-1 PASS, (EQ-1)↔(EQ-2) self-consistency MUSI converge z native 3D bg form (without exogenous D/L^α addition).

**Tolerancja:** Analytical demonstration OR numerical convergence w ≥ factor 10 parameter range.

**Severity:** STRUCTURAL.

**Conditional activation:** F-γ-2 testowane TYLKO jeśli F-γ-1 PASS. Jeśli F-γ-1 FAIL → F-γ-2 NIE applicable.

### F-γ-3, F-γ-4 — NOT activated (cosmological scope)

F-γ-3 (H_0 PRIMARY KILLER) + F-γ-4 (confinement/deconfinement match observed) pozostają **PENDING_POZIOM_GAMMA_3** lub późniejsze.

Niniejszy cykl scope: γ-1 = F-γ-1 + (warunkowo) F-γ-2. Cosmological extension = osobny cykl post-γ-1/γ-2.

---

## §2 — 10 forbidden post-hoc moves (inherited z Poziom β + R2 audit)

Jakikolwiek z poniższych w trakcie γ-1 = **automatyczny HALT-B**:

1. Modyfikacja F-γ-1 power-law/log/Coulomb threshold ex post
2. Renaming "exponential" to "screened power-law" by avoid FAIL
3. Adding "additional fields" by retain failed F-γ-1
4. Cherry-picking ansatz that gives power-law while ignoring alternative ansatzes
5. Re-defining "long-range" by include short-range mid-distance regime
6. Switching defect type mid-cycle (hedgehog → vortex → other) by avoid FAIL
7. Hardcoding FP T_pass=True (strict cycle 1/2/7 violation)
8. Using DEC budget powyżej 1 (max budget exceeded)
9. Introducing new axioms by rescue F-γ-1 (R3 threshold violation)
10. Bypassing CALIBRATION §3.6 analytical pre-derivation step (NEW BINDING 2026-05-22)

---

## §3 — L1 / L2 / L3 falsification map

### L1 (native, cycle-local)
- F-γ-1: native 3D long-range interaction (CRUCIAL)
- F-γ-2: self-consistency closure (SECONDARY, conditional)

### L2 (framework targets, post-cycle)
- TGP-native (EQ-1)-(EQ-6) self-consistency w 3D
- CE-H structural feature confirmed quantitatively w 3D
- C6 PARTIAL FFS → RESOLVED_STRUCTURALLY full closure (conditional)

### L3 (cross-cycle propagation, conditional on success)
- CE-H Poziom β A− → A upgrade trajectory (post γ-1 PASS)
- FFS cycle A− → A upgrade trajectory (post γ-1 + γ-2 PASS)
- Poziom γ-2 (self-consistency native bg) authorization
- Poziom γ-3 (cosmological extension F-γ-3 H_0) authorization gate

---

## §4 — Phase plan (5 substantive faz + Phase 0 + Phase FINAL)

### Phase 0 — Balance sheet + ANALYTICAL PRE-DERIVATION (§3.6 BINDING)
**Scope:** External inputs, LOCKED structural axioms, derived outputs, tautology/falsifiability checks, anti-BD-drift check, **analytical pre-derivation dla F-γ-1 numerical threshold** (NEW BINDING per CALIBRATION §3.6 2026-05-22).
**Critical:** Phase 0 MUST include explicit analytical derivation of expected V_int(L) far-field behavior z S05+Z₂+U(1)+RP² Lagrangianu BEFORE any sympy. Documented symbolically.
**Deliverable:** `Phase0_balance.md`
**Estimated:** 1 dzień (today).

### Phase 1 — 3D hedgehog/vortex ansatz + far-field linearized analysis
**Scope:** Native ansatz dla single defect (decision DEC: hedgehog point defect OR vortex line defect — pre-registered §5). Far-field expansion. Linearized EL equations. Mass spectrum (Higgs massive + Goldstone massless lub gauge field).
**Deliverable:** `Phase1_sympy.py` + `Phase1_results.md`
**Estimated:** 2-3 dni.

### Phase 2 — Two-defect interaction at large separation
**Scope:** Compute V_int(L) between two static defects z full 3D Lagrangianu. Propagator analysis. Identify dominant long-range mode.
**Deliverable:** `Phase2_sympy.py` + `Phase2_results.md`
**Estimated:** 2-3 dni.

### Phase 3 — Differential test vs 1D Z2 baseline
**Scope:** Explicit fit comparison: power-law (R²_power) vs exponential (R²_exp) for 3D V_int(L). Comparison to 1D Z2 baseline (CE-H Poziom β Phase 3 result: pure exp(-m·√2·L)). F-γ-1 verdict.
**Deliverable:** `Phase3_sympy.py` + `Phase3_results.md`
**Estimated:** 1-2 dni.

### Phase 4 (conditional) — Self-consistency closure z native bg (F-γ-2)
**Scope:** ONLY IF F-γ-1 PASS. (EQ-1)↔(EQ-2) self-consistency w 3D z native bg.
**Deliverable:** `Phase4_sympy.py` + `Phase4_results.md`
**Estimated:** 2-3 dni (conditional).

### Phase FINAL — Closure
**Scope:** F-γ-1 verdict, F-γ-2 verdict (if activated), claim_status assignment, cross-cycle propagation map (C6 disposition, parent cycle upgrades, etc.).
**Deliverable:** `Phase_FINAL_close.md`
**Estimated:** 0.5 dnia.

**Total estimated:** 8-12 dni (F-γ-1 only: 6-9 dni; z F-γ-2: +2-3 dni).

---

## §5 — Discipline declarations + pre-registered ansatz choice (DEC budget)

### Strict cycle 1/2/7 pattern
- 0 hardcoded FP T_pass=True
- LIT/INVENTORY informational only
- Substantive FP tests compute-then-compare

### Max DEC budget = 1
- **Reserved dla DEC choice:** hedgehog point defect vs vortex line defect
- Pre-registered preference: **vortex line defect** (analog Nielsen-Olesen + Vilenkin-Shellard standard cosmic string framework)
- Rationale: U(1) phase symmetry natively gives vortex strings; hedgehog requires σ_ab + RP² extension (more complex, deferred)
- Alternative: jeśli vortex analysis fails substantively, hedgehog as fallback (DEC budget commits)

### R1+R2+R3 discipline BINDING (CALIBRATION §6 2026-05-22)
- R1 (research-tier): all phases preserve permissive flagging
- R2 (integration audit): scope = items aggregated z R1, deferred do post-cycle audit
- R3 (multi-line convergence): CE-H structural feature already accepted 3/3 lines; γ-1 verifies quantitatively

### Anti-BD-drift LOCK
- TGP S05+Z₂+U(1)+RP² Lagrangian only
- NO fitting do QCD/Nielsen-Olesen/Vilenkin-Shellard (literature jest informational anchor, NIE target)
- Native equations FIRST; mapping post-hoc bonus only

### Anti-Lakatos LOCK
- Pre-registration LOCKED 2026-05-22 before any sympy
- Each phase reports honestly vs pre-registration
- Any tolerance modification ex post = HALT-B
- §3.6 analytical pre-derivation BINDING (forbidden #10)

### Pre-registered ansatz choice (vortex line)

**Default ansatz (DEC reserved):**

Single vortex line z fractional winding n along z-axis (z-translation invariant):
$$\Phi_{vortex}(r, \phi) = \rho(r) \cdot e^{i n \phi}$$

z $\rho(r) \to 0$ at r=0 (string core) + $\rho(r) \to v$ at r → ∞ (mexican hat vacuum).

**Far-field BC:** $\Phi \to v \cdot e^{i n \phi}$ at r → ∞.

**Two-vortex configuration:** Two parallel z-axis vortices przy x = ±L/2.

**Modified DEC option (point defect / hedgehog):** Jeśli vortex analysis daje pure exponential (potential FAIL), DEC switch do point hedgehog (σ_ab + Φ joint config) z degree-1 winding na S². To DEC commit (1/1 budget) — used tylko jeśli necessary.

---

## §6 — Risk register

| ID | Risk | Mitigation | Severity |
|----|------|-----------|----------|
| R1 | Vortex analysis daje pure exponential at long-range → F-γ-1 FAIL → HARD HALT scenario | Pre-registered honest scenario; CE-H bg required exogenously even w 3D → fundamental redesign | **HIGH** (acknowledged) |
| R2 | Multiple modes coupling complicates analytical solution | Numerical fallback z explicit fit comparison Phase 3 | MEDIUM |
| R3 | Goldstone vs Higgs mode decomposition ambiguity | Phase 0 analytical pre-derivation MUST clarify (§3.6 BINDING) | MEDIUM |
| R4 | DEC budget exhausted by hedgehog/vortex DEC + can't switch later | Pre-registered preference clear (vortex); switch only on substantive FAIL of vortex | LOW |
| R5 | F-γ-1 PASS marginal (R²_power ≈ R²_exp) — interpretation-dependent | Pre-register clear threshold: power-law fit must be at least as good as exponential at large L | MEDIUM |
| R6 | 3D cosmic string theory literature inheritance temptation (BD-drift) | Anti-BD-drift LOCK; native equations FIRST; literature informational only | MEDIUM |
| R7 | Sympy 3D calculations computationally expensive | Phase 0 plan numerical fallback strategy; explicit time budget per phase | LOW |
| R8 | C6 over-claim post γ-1 PASS (FFS Phase 4 resolved fully) | Explicit "RESOLVED_STRUCTURALLY pending γ-2" + Poziom γ-3 cosmological for full A− → A | MEDIUM |
| R9 | Path B → Path C planning premature mid-cycle | NO mid-cycle planning; Path C decision post γ-1 FINAL closure | LOW |
| R10 | Methodology drift (forget §3.6 pre-derivation in later phases) | Phase 0 explicit template; each phase header references §3.6 BINDING | LOW |

---

## §7 — Authorization gate dla Phase 1

Sekwencja autoryzacji:
- ✅ 2026-05-22: Poziom γ-1 scaffold (README + Phase 0) authorized by "Tak działamy ze ścieżką B"
- ⏳ Phase 1 sympy → wymaga osobnego potwierdzenia LUB batch authorization
- ⏳ Phase 2 sympy → wymaga osobnego potwierdzenia LUB batch
- ⏳ Phase 3 sympy → wymaga osobnego potwierdzenia LUB batch
- ⏳ Phase 4 (conditional F-γ-2) → wymaga F-γ-1 PASS
- ⏳ Phase FINAL → wymaga osobnego potwierdzenia

**Alternative:** Batch authorization (analog Poziom β single-session) możliwe jeśli user prefers efficiency.

---

## §8 — Methodology demarcation (post §3.6 BINDING)

| Aspect | Pre-§3.6 cycles | This cycle (post-§3.6 BINDING) |
|--------|-----------------|--------------------------------|
| Pre-registration | Numerical thresholds + heuristic intuition | Numerical thresholds **+ analytical pre-derivation** explicit |
| Phase 0 | External inputs + ansatz + literature | Same **+ analytical derivation slot** dla each FP-class falsifier |
| Risk of T_P3_2-like fail | Possible (FFS-3 + CE-H T_P3_2 patterns) | Pre-empted via §3.6 analytical step |
| Anti-Lakatos discipline | LOCKED | LOCKED **+ §3.6 enforcement** |

**This cycle = first cycle post-§3.6 BINDING.** Methodology rigor explicit visible.

---

## §9 — Conditional outcomes summary

| Scenario | F-γ-1 verdict | Cycle outcome | Parent cycle impact | Path C planning |
|----------|---------------|---------------|---------------------|-----------------|
| **A** | PASS (power-law/log/Coulomb clearly) + F-γ-2 PASS | claim_status A | CE-H Poziom β A−→A; FFS C6 RESOLVED conditional na cosmology γ-3 | Path C = Poziom γ-3 cosmological extension (H_0 PRIMARY KILLER) |
| **B** | PASS + F-γ-2 PARTIAL | claim_status A- | CE-H Poziom β A− → A− preserved | Path C = γ-2 retry or γ-3 with explicit caveat |
| **C** | PASS marginal (R²_power ≈ R²_exp) | claim_status A- conditional | CE-H A− preserved; FFS C6 RESOLVED_STRUCTURALLY pending γ-3 | Path C = decision required |
| **D HALT** | FAIL (pure exponential) | HALT-B scenario | CE-H Poziom β honest caveat strengthened; bg form NIE native | Path C = redesign or abandon CE-H framework cosmological extension |

**Pre-registered honest discipline:** Każdy scenariusz dopuszczalny w kontrakcie. Anti-Lakatos LOCKED — żaden scenariusz NIE rescued ex post.

---

**END OF README — Poziom γ-1 BINDING contract LOCKED 2026-05-22**

**Authorization point:** Phase 1 (Phase 0 → execute now jako part of scaffold + balance sheet zgodnie z R2 audit single-session precedent).
