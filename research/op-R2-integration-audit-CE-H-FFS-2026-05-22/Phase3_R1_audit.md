---
title: "Phase 3 — R1 flag audit (op-R2-integration-audit-CE-H-FFS-2026-05-22)"
type: phase_audit
status: COMPLETE
date_completed: 2026-05-22
phase: 3
parent_cycle: op-R2-integration-audit-CE-H-FFS-2026-05-22
items_audited: 1 (R1-1)
---

# Phase 3 — R1 flag audit

**Status:** COMPLETE 2026-05-22.
**Purpose:** Audit R1-1 item (Phase 3 Poziom β pre-registration analytical pre-derivation gap) + draft methodology addendum.

---

## §0 — Verdict summary

| Item | Verdict | Justification (short) |
|------|---------|------------------------|
| R1-1 | ✅ **CLOSED** | Root cause documented (tanh tail asymptotic); pattern observed twice (CE-H T_P3_2 + FFS Phase 4 q=1 implicit); methodology addendum drafted (§3 below) |

**Aggregate Phase 3:** 1 CLOSED + 0 DEFERRED + 0 ESCALATED (1/1 items).

---

## §1 — R1-1: Pre-registration analytical pre-derivation gap

### §1.1 Status assessment (Path A)

**Source:** CE-H Poziom β Phase 3 T_P3_2 honest fail 2026-05-21.

**Structural facts (from Phase 3 results §2.2-§2.3):**

**Pre-registered threshold (Phase 3):**
> "Native 1D Z2 interaction has form V_int ~ -A·exp(-m·L), fitted decay rate should match m within 10%, expected m_num = 1.0"

**Actual fitted result:**
- V_int ~ -10.04 · exp(-1.40·L), R² = 0.9999
- Fitted decay rate: 1.40
- Pre-registered tolerance: m_fit ∈ [0.9, 1.1] dla m=1.0
- 1.40 > 1.1 → formally FAILS pre-registered threshold

**Analytical correction (Phase 3 §2.3):**
- Kink ansatz $\Phi_K(x) = v \cdot \tanh(m \cdot x / \sqrt{2})$ has tail:
  $$\Phi_K(x) = v - 2v \cdot e^{-m \sqrt{2} \cdot x} + O(e^{-2m\sqrt{2} \cdot x})$$
- So kink-antikink interaction has tail $V_{int}(L) \sim e^{-m \sqrt{2} \cdot L}$
- Analytical decay rate: $m \cdot \sqrt{2} \approx 1.4142$
- Fitted 1.40 vs analytical 1.4142 = **match within 1%**

**Substantive conclusion:** Substance verified at 1% accuracy. ONLY pre-registration threshold was analytically wrong (expected m instead of m·√2).

### §1.2 Root cause analysis (Path A + B)

**Why did pre-registration use m instead of m·√2?**

Pre-registration was based on **heuristic intuition** ("interaction should decay z field mass scale") rather than **explicit analytical derivation** of tail behavior.

Steps that would have caught the error:
1. Write kink ansatz $\Phi_K(x) = v\tanh(mx/\sqrt{2})$
2. Compute asymptotic expansion at large x
3. tanh(y) → 1 - 2e^{-2y} + O(e^{-4y}) for y → ∞
4. With y = mx/√2: $\tanh(mx/\sqrt{2}) \to 1 - 2 e^{-m\sqrt{2}\cdot x}$
5. → tail decay rate is **m·√2**, NIE m
6. Two-kink interaction has same decay rate (exponentially suppressed by tail)

Pre-registration missed step 1-5 — went directly from "interaction with mass scale m" to "exp(-mL)".

### §1.3 Pattern check (Path C cross-cycle)

**Same class of methodology issue observed in:**

**(a) CE-H Phase 3 T_P3_2 (this item)**
- Pre-registered: decay rate m
- Native: m·√2 (factor √2)
- Match w 1%, threshold failed

**(b) FFS pre-screening T7 σ formula (FFS-3 audit Phase 1)**
- Pre-registered: σ = π·v² (implicit q=1 effective)
- Native strict Nielsen-Olesen: σ = π·q²·v² (q² factor explicit)
- For q=1/3: ratio factor ~10 (still within "factor 10" pre-screening claim)

**Pattern:** Pre-registrations occasionally use **convenient approximations** (implicit q=1, decay rate m) instead of **analytically-derived exact values** (q² scaling, m·√2 tail).

**Pattern source:** Heuristic shortcut in pre-registration drafting — pre-registration relies on physical intuition rather than explicit analytical pre-derivation.

### §1.4 Pre-registration check (Path D)

**Question:** Czy methodology framework explicit required analytical pre-derivation?

**Existing methodology (per CALIBRATION_PROTOCOL §3):**
- Pre-registration timestamp PRZED any sympy ✓
- Falsifiers explicit z numerical thresholds ✓
- Anti-Lakatos discipline LOCKED ✓
- **BUT:** explicit requirement dla analytical pre-derivation step **NIE explicit** in current methodology

→ Methodology gap identified. Addendum needed.

---

## §2 — Honest declaration: this audit is NOT Lakatos rescue

**Critical anti-Lakatos check:**

The audit identifies R1-1 as METHODOLOGY ISSUE (pre-registration drafting process), NIE rationalize the FAIL. The fact that fitted 1.40 matches analytical 1.4142 within 1% is **substantive observation**, NIE post-hoc threshold modification.

**What was NOT done (per anti-Lakatos LOCK):**
- ✗ NIE modified threshold ex post
- ✗ NIE re-defined "fitted" to include analytical correction
- ✗ NIE hidden the FAIL
- ✗ NIE re-ran sympy with adjusted tolerance

**What was done (anti-Lakatos preserved):**
- ✓ Reported T_P3_2 jako FAIL per literal pre-registered threshold
- ✓ Documented analytical correction as **methodology lesson**, NIE rescue
- ✓ Substance verified at 1% accuracy noted explicitly
- ✓ R1 flag created dla R2 audit (proper procedure)

**The 1% accuracy verification HAS substantive scientific value** — it confirms native 1D Z2 exponential structure. But the **pre-registration FAIL stands** per anti-Lakatos discipline.

---

## §3 — Methodology addendum draft

### §3.1 Proposed addendum to CALIBRATION_PROTOCOL §3

**Current §3 (Anti-Lakatos pre-registration discipline):**
1. Pre-registration timestamp PRZED any sympy
2. Falsifiers explicit z numerical thresholds
3. Forbidden post-hoc moves enumerated
4. (...)

**Proposed addition §3.6 — Analytical pre-derivation step:**

> **§3.6 — Analytical pre-derivation step (BINDING for falsifiers z numerical thresholds)**
>
> For any pre-registered falsifier specifying numerical threshold for a derived quantity (e.g., "fitted decay rate match m within 10%", "ratio within factor 2"), pre-registration MUST include:
>
> (a) **Analytical derivation** of expected value from underlying ansatz/field theory
> (b) **Symbolic computation** demonstrating expected numerical value
> (c) **Documentation w Phase 0** alongside numerical threshold
>
> **Forbidden shortcut:** Heuristic intuition ("interaction should decay z mass scale m") without explicit analytical derivation.
>
> **Example (R1-1 source):**
> - ❌ "fitted decay rate match m within 10%" (heuristic; expected m=1)
> - ✅ "fitted decay rate match m·√2 within 10%" (analytical: tanh tail decays as e^{-m√2·x}; expected m·√2 ≈ 1.4142)
>
> **Anti-Lakatos preservation:** This step PRE-EMPTS post-hoc threshold modification. Analytical pre-derivation at pre-registration time = honest expectation; failure → real structural finding.
>
> **Permitted relaxation:** Heuristic threshold acceptable IF marked explicitly as "informational only" (LIT/INVENTORY class), NIE FP class.

### §3.2 Cross-cycle propagation requirement

Future cycles should:
1. Include analytical pre-derivation in Phase 0 balance sheet
2. Document expected values symbolically (sympy or LaTeX in Phase 0)
3. Pre-registered thresholds reference analytical derivation explicit
4. Phase X results compare to ANALYTICAL value (not heuristic)

### §3.3 Retrospective application (DEFERRED)

**Question:** Should this addendum be retroactively applied to closed cycles?

**Answer:** NO retroactive re-evaluation (forbidden post-hoc move analog).

**Justification:** Closed cycles preserve their pre-registered LOCKS. Retroactive application would VIOLATE anti-Lakatos discipline. Going forward, new cycles use the addendum.

**Exception:** Cross-cycle documentation MAY annotate previous patterns (FFS-3, CE-H T_P3_2) as "methodology lesson", NIE re-evaluate verdicts.

---

## §4 — Pattern repeat risk assessment

### §4.1 Will pattern repeat?

**Risk factors:**
- Pre-registration often done under time pressure (heuristic shortcut tempting)
- Analytical pre-derivation requires extra Phase 0 work
- Sympy gives "verified" feeling — may compensate for weak pre-registration

**Mitigation (per addendum §3.1):**
- Phase 0 BALANCE sheet template should include analytical pre-derivation explicit slot
- Phase 0 author should verify analytical correctness BEFORE pre-registration LOCK
- R2 audit should be alert dla new instances of this pattern

### §4.2 Methodology robustness check

**Two instances observed (CE-H T_P3_2 + FFS pre-screening q=1):**
- Both honestly reported (FAIL per literal threshold OR explicit honest amendment)
- Anti-Lakatos discipline operationally preserved BOTH times
- Substantive analysis revealed match within reasonable accuracy BOTH times
- → Methodology IS catching these correctly via R1+R2+R3

**Conclusion:** Pattern is detectable + correctable z current methodology + addendum. NO need to escalate dla new structural axiom.

---

## §5 — Verdict R1-1: CLOSED

**Justification:**
- ✅ Root cause documented: heuristic pre-registration shortcut (intuition m instead of analytical m·√2)
- ✅ Pattern observed twice (CE-H T_P3_2 + FFS Phase 4 q=1 implicit)
- ✅ Methodology addendum drafted (§3.1 above) — pre-registration analytical pre-derivation step
- ✅ Pattern repeat risk assessed + mitigation in addendum
- ✅ Anti-Lakatos discipline preserved (no rescue, no threshold modification)
- ✅ Future cycles use addendum; closed cycles preserve their LOCKs

**Decision matrix mapping:** "Root cause documented + methodology addendum drafted" → **CLOSED** per §3 README.

**Cross-cycle propagation implications (deferred do Phase 4):**

| Doc target | Impact from Phase 3 |
|------------|---------------------|
| meta/CALIBRATION_PROTOCOL.md §3 | Add §3.6 "Analytical pre-derivation step" addendum |
| op-CE-H-two-particle-equilibrium-2026-05-21/Phase_FINAL_close.md §5 | Note R1-1 CLOSED status |
| meta/CYCLE_KICKOFF_TEMPLATE.md (jeśli existuje) | Add analytical pre-derivation slot do Phase 0 template |

---

## §6 — Anti-Lakatos discipline check

- ✅ Methodology addendum drafted, NIE retroactively applied to closed cycles
- ✅ No threshold modifications ex post (T_P3_2 FAIL stands per literal threshold)
- ✅ Substantive 1% accuracy noted as scientific observation, NIE rescue
- ✅ Pattern detection (FFS-3 + CE-H T_P3_2) confirms R1+R2+R3 working
- ✅ Cross-cycle inheritance preserved

**Self-audit:** Phase 3 audit clean, anti-Lakatos preserved.

---

**END OF PHASE 3 — R1 flag audit LOCKED 2026-05-22**

**Aggregate Phase 3:** 1 CLOSED. Ready dla Phase 4 (cross-cycle propagation execution).
