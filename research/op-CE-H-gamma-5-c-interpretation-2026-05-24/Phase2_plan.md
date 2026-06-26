---
title: "Phase 2 — c(n_local) derivation z entropy-based crayon box formalism (HANDOFF §3.7)"
type: phase_plan
status: LOCKED
phase: 2
parent_cycle: op-CE-H-gamma-5-c-interpretation-2026-05-24
authorization_date: 2026-05-24
authorization_quote: "działaj Phase 2"
dec_used: "DEC 2 (Configuration space counting method — committed)"
---

# Phase 2 — c(n_local) entropy-driven derivation z crayon box formalism

**Status:** LOCKED 2026-05-24. **In progress.**

**Authorization:** User explicit "działaj Phase 2" (2026-05-24).

---

## §1 — Phase 2 objectives (per Phase 0 §4 + README §4)

1. Configuration space Ω(n) counting (crayon box analogy formal HANDOFF §3.7)
2. Connect Ω(n) → c_eff(n) via substrate reconfiguration rate
3. Identify n_critical where Ω(n) → 1 (event horizon emergence)
4. Compare against pre-registered candidate forms L1-L5 (Phase 0 §5)
5. F-γ-5-D verification (saturating limit at 0 + critical decay)

---

## §2 — DEC 2: Configuration space counting method (COMMITTED)

**Pre-registered choice (per Phase 0 §12):** Configuration counting via **available-slot enumeration** (entropy of unoccupied substrate positions).

**Justification:** HANDOFF §3.7 + §10 user crayon box analogy:
> "im więcej kredek, tym mniej możliwości rekonfiguracji, przy pełnym pudełku, jakiekolwiek potrząsanie nie może wywołać rekonfiguracji"

**MORE sources** in given volume = **FEWER available substrate positions** for reconfiguration = **LOWER propagation rate**.

**Mathematical model (committed PRZED sympy):**

Substrate volume V holds at most N_max = n_critical · V "source slots". With n sources currently occupied (n ≤ N_max), available slots = N_max - n.

**Available reconfiguration moves per source:**
ω(n) = (N_max - n) per source × n sources = n·(N_max - n) total moves

**Normalized reconfiguration capacity:**
Ω_reconfig(n) / Ω_max = (N_max - n) / N_max = 1 - n/n_critical

(Per-source moves; n_critical · V = N_max gives ratio depending only on fraction occupied.)

---

## §3 — Derivation strategy

### §3.1 — Substrate reconfiguration rate

**Definition (committed):** Effective propagation rate c_eff(n_local) = base substrate rate c_0 × "reconfiguration capacity fraction" Ω(n)/Ω_max.

**Per HANDOFF §3.7 entropy interpretation:**
| n_local | Ω (config space) | c_eff |
|---------|------------------|-------|
| n=0 (deep space) | Ω_max | c_0 |
| n << n_critical | Ω_max·(1-n/n_critical) ≈ Ω_max | c_0(1-ε) |
| n = n_critical/2 | Ω_max/2 | c_0/2 |
| n → n_critical | Ω → 0 (single config) | c_eff → 0 |

**Pre-registered formula (DERIVED from slot-count entropy):**

c(n)/c_0 = Ω(n)/Ω_max = 1 - n/n_critical

**Boundary verification (analytical PRZED sympy):**
- c(n=0) = c_0·(1-0) = c_0 ✓ (matches F-γ-5-D property (i): deep space limit)
- c(n=n_critical) = c_0·(1-1) = 0 ✓ (matches F-γ-5-D property (ii): event horizon emergence)
- ∂c/∂n = -c_0/n_critical < 0 ✓ (matches F-γ-5-D property (iii): monotonically decreasing)

### §3.2 — Mapping to pre-registered forms L1-L5

| Form (pre-registered) | Expression | Match z derived? |
|----------------------|------------|------------------|
| **L1** (β=1) | c_0 · (1 - n/n_critical)^β | **MATCH at β=1** ✓ |
| **L2** (γ=1) | c_0 · (Ω(n)/Ω_max)^γ | **MATCH at γ=1** ✓ (since Ω(n)/Ω_max = 1-n/n_c linearly) |
| L3 | c_0 · exp(-α·n/(n_critical-n)) | NIE match (essential singularity, NIE polynomial) |
| L4 | c_0 · S(n)/S_max | Different (linear in entropy, NIE in slot count) |
| L5 (general p) | c_0 · (1 - (n/n_critical)^p) | NIE match unless p=1 (then form differs anyway) |

**Derived form matches L1 (β=1) AND L2 (γ=1) simultaneously** — both reduce to linear form when slot count = entropy proxy.

**Classification:** **CONFIRMED_FORM_L1_LINEAR (β=1)** ≡ CONFIRMED_FORM_L2_LINEAR (γ=1).

### §3.3 — n_critical scaling (dimensional analysis)

**TGP-native candidates dla n_critical:**

| Candidate | Expression | Physical interpretation |
|-----------|------------|--------------------------|
| **A** | n_critical = 1/ℓ_P³ ≈ 2.4×10¹⁰⁵ /m³ | Planck density (each substrate cell ~ ℓ_P size) |
| B | n_critical = (m_σ·c/ℏ)³ | σ-mode Compton density |
| C | n_critical = c²/(G·m_eff) | GR-matching density (assuming m_eff per source) |

**Phase 2 commits to Candidate A:** n_critical = ℓ_P^(-3).

**Justification (Appendix E Thm:natural-cutoff):** Substrate has natural Planck-scale lattice spacing a_sub = ℓ_P; each "source slot" occupies one Planck-volume cell. This is **derivable from TGP first principles** (not external input).

**Honest disposition (Candidate C):** If Phase 5 requires Candidate C (GR-matching n_critical) to satisfy F-γ-5-A, declare honestly. Anti-Lakatos: NIE pre-bias choice of n_critical to save F-γ-5-A.

### §3.4 — Weak-field connection (Phase 5 input preview)

**For mass M at distance r, far-field Yukawa-Phi potential (γ-3' Phase 1 Mechanism C):**

⟨Φ⟩(r) ∝ M/r (assuming massless or far-field Yukawa)

**Local "source density" felt by propagation at distance r:**

n_eff(r) ∝ ⟨Φ⟩(r)/v ∝ M/(r·v)

**Effective c shift:**

c(r)/c_0 = 1 - n_eff(r)/n_critical = 1 - α·M/r

where α encodes coupling and n_critical normalization.

**This gives linear δt/t = α·M/r — GR-LIKE scaling!** Phase 5 detailed analysis.

---

## §4 — Substantive FP tests (Phase 2 sympy execution plan)

**Strict cycle 1/2/7:** 0 hardcoded T_pass=True. Each FP computes-then-compares.

| FP ID | Test | Compute | Compare against |
|-------|------|---------|------------------|
| T_P2_1 | Slot-count entropy: available slots = N_max - n | Symbolic: substract n from N_max | Expected: N_max - n |
| T_P2_2 | Reconfiguration capacity ratio | (N_max - n)/N_max | Expected: 1 - n/n_critical (where n_critical = N_max/V) |
| T_P2_3 | c_eff formula derivation | c_0·Ω(n)/Ω_max | Expected: c_0·(1-n/n_critical) |
| T_P2_4 | F-γ-5-D (i): c(n=0) = c_0 | Substitute n=0 | Expected: c_0 |
| T_P2_5 | F-γ-5-D (ii): c(n=n_critical) = 0 | Substitute n=n_critical | Expected: 0 |
| T_P2_6 | F-γ-5-D (iii): ∂c/∂n < 0 | Symbolic derivative | Expected: -c_0/n_critical (negative) |
| T_P2_7 | F-γ-5-D combined verification | Combined from prior FPs | All 3 properties PASS expected |
| T_P2_8 | n_critical dimensional analysis | ℓ_P^(-3) units check | Expected: [length]^(-3) ≡ [density] |
| T_P2_9 | Weak-field Yukawa connection (linear δt/t in M/r) | c(r)/c_0 = 1 - α·M/r for point source | Expected: GR-like linear scaling |
| T_P2_10 | Linear regime sanity: at n/n_critical=0.001, c/c_0=0.999 | Substitute n/n_c = 10⁻³ | Expected: 0.999 (linear regime) |

**Pre-registered PARTIAL outcomes:**
- PARTIAL_compute: 0 expected
- PARTIAL_concept_mismatch: 1 possible (if n_critical scaling choice not uniquely TGP-derivable, declare honestly)

---

## §5 — Risks acknowledged

1. **β=1 (linear form) IS A CHOICE.** Alternative β values (β=2 percolation, β=3 dimensional cube) physically motivated for different substrate models. Phase 2 commits to β=1 as the cleanest "available slot count" interpretation. Phase 5 may revise IF F-γ-5-A/B require different β. **Anti-Lakatos:** Pre-register β=1; declare honestly if later modification needed.

2. **n_critical = ℓ_P^(-3) is ANSATZ.** Alternative scalings (m_σ Compton, GR-matching c²/G·m_eff) are physically plausible. Phase 2 commits to Planck density as the most TGP-natural choice (per Appendix E Thm:natural-cutoff). Phase 5 will test against F-γ-5-A/B; if FAIL, may need to revise (honest disposition).

3. **Mean-field approximation.** Phase 2 treats sources as uniform density n_local; ignores spatial fluctuations + correlations. Per §3.6.8 BINDING (implicit assumptions enumeration): mean-field approximation declared explicit; fluctuation corrections deferred to future work.

4. **Connection to gravity (Phase 3 + 5) requires additional machinery.** Phase 2 only derives c(n_local) form; Phase 3 ZA gravity-as-configuration-constraint + Phase 5 dla R_s/δt/t numerical evaluation.

---

## §6 — Phase 2 deliverables

- ✅ Phase2_plan.md (this document)
- ☐ Phase2_sympy.py (substantive FP tests)
- ☐ Phase2_results.md (verdicts + cross-reference to candidate forms)

---

**END OF PHASE 2 PLAN — c(n_local) derivation strategy LOCKED 2026-05-24**

**Ready dla sympy execution.**
