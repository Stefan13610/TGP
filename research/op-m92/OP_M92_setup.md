# OP-M92 — M9.2 conditional pivot framework (Track B preparation)

**Data otwarcia:** 2026-04-25
**Status:** **OPEN — DEFERRED RESEARCH PROGRAM**
**Trigger 1:** OP-EHT-A NEGATIVE on naive proper-time path ✓ **SATISFIED 2026-04-25**
**Trigger 2:** ngEHT 2030+ confirms GR Sgr A* shadow at 1% precision — **DEFERRED 2030–2032**
**Cel:** Przygotować formalny scope, koszty strukturalne, decision tree i kandydaty dla pivota M9.1'' → M9.2, w razie potwierdzenia GR shadow przez ngEHT.

**Cross-references:**
- [[research/op-eht-A/OP_EHT_A_final_verdict.md]] — Track A NEGATIVE verdict (parent trigger 1)
- [[research/op-eht/OP_EHT_final_verdict.md]] — OP-EHT closure (CONDITIONAL POSITIVE)
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] — M9.1'' weak-field PASS (must preserve)
- [[research/op7/OP7_T6_results.md]] — c_GW = c_0 (must preserve)
- [[paper/tgp_core.tex]] §applications BH shadows + F4 falsifiability
- [[KNOWN_ISSUES.md]]

---

## 1. Strategiczny kontekst (post OP-EHT-A)

OP-EHT zamknięty CONDITIONAL POSITIVE (13/18 = 72%) z explicit założeniem
że Track A *może* pochodzić z first principles. OP-EHT-A NEGATIVE (8/12 = 67%,
strong-field test T-A3 FAIL) potwierdza że **naive proper-time matter coupling
overshoots scenario (e) target** — magnitude factor 0.652 vs target 0.886,
prowadząc do b_crit deviation **−25.32% vs GR** (wrong direction).

Kluczowe odkrycia z OP-EHT-A:
1. **ψ_ph jest geometric invariant** (1.168) niezależnie od A_eff; self-consistent
   iteration jest trivial single-step.
2. **Weak-field PPN i c_GW = c_0 auto-preserved** pod Track A — nie ma "tanich"
   ścieżek rescue które psują M9.1'' P3 audit lub OP-7 closure.
3. **Scenario (e) coupling f(ψ) = √(g_tt^GR/g_tt^TGP) jest non-local w r** —
   wymaga porównania TGP photon ring (r_ph_TGP = 3.88 M) do GR photon ring
   (r_ph_GR = 3 M) jako benchmark, czego nie ma natural M9.1'' Lagrangian.

**Implikacja:** M9.1'' standalone **nie ma natural rescue path**. Jeśli ngEHT 2030+
potwierdzi GR shadow at 1% precision, M9.1'' falsifies as standalone theory
i M9.2 pivot staje się **unconditional**. OP-M92 przygotowuje to scenariusze
*zanim* nadejdzie verdict, by reagować w 2-4 weeks zamiast 6-12 months.

## 2. Scenario (e) coupling — strukturalna analiza non-lokalności

OP-EHT T3 fitted ansatz:
```
f(ψ) = √(g_tt^GR(r) / g_tt^TGP(r)) = √( (1 - 2M/r) × ψ / (4 - 3ψ) )
```

**Non-local features:**
- `g_tt^GR(r) = -(1 - 2M/r)` — Schwarzschild backbone
- `g_tt^TGP(r) = -(4 - 3ψ(r))/ψ(r)` — M9.1'' metric form
- f(ψ(r)) wymaga *both* metryk **at the same r** — non-trivial matter "vision"

W M9.1'' substrate-jet nie ma GR backbone jako reference object. Aby uzyskać
takie sprzężenie z first principles, M9.2 musi wprowadzić **drugą strukturę**
która mimicuje GR-like response w strong-field bez naruszenia weak-field PPN.

**Cztery strukturalne ścieżki kandydujące** identyfikowane w OP-M92:

### Kandydat A: Dual-field structure (Φ + σ z non-minimal coupling)

```
S = ∫ [ -(∇Φ)² /(8πG) + α (∇σ)² + V(Φ, σ) + L_mat[Φ, σ] ] √-g d⁴x
```

- Drugi skalarny field σ(x) generuje "pseudo-GR backbone" w strong-field
- Matter coupling f(Φ, σ) może mimic √(g_tt^GR/g_tt^TGP)
- **Koszt strukturalny:** narusza single-Φ Z_2 substrate principle
  (TGP_FOUNDATIONS.md §1.2 fundamental ansatz: jeden field Φ)
- **Risk:** Brans-Dicke-style equivalence principle constraints (Cassini γ_PPN)

### Kandydat B: Conformal frame coupling (Jordan-frame matter)

```
g_μν^Jordan = Ω²(ψ) g_μν^Einstein
L_mat[Ω(ψ) × g_μν, ψ_matter]
```

- Conformal factor Ω(ψ) modyfikuje effective metric for matter
- Może dać scenario (e) jeśli Ω²(ψ_ph) = √(g_tt^GR/g_tt^TGP)
- **Koszt strukturalny:** Ω(ψ→1) ≠ 1 → 1PN PPN shift; preservation γ_PPN = 1
  wymaga fine-tuning
- **Risk:** 5th force constraints (composition-dependent acceleration <1e-13)

### Kandydat C: ψ-dependent matter charge q(ψ) (q-flow)

```
L_mat = -q(ψ) × φ × ρ / Φ_0
gdzie q(ψ→1) = q_0 (1PN preserved)
       q(ψ_ph) = q_0 × √(g_tt^GR/g_tt^TGP)|_{ψ_ph}
```

- Variable matter charge skaluje effective coupling w strong-field
- **Koszt strukturalny:** ad-hoc unless q(ψ) wynika z głębszej zasady
- **Risk:** Quantum coherence (q variations affect atom spectroscopy)
- **Otwarte pytanie:** czy q(ψ) może wynikać z β-function flow w renormalization
  group sense? (analogiczne do running coupling QED w polu silnym)

### Kandydat D: Momentum back-reaction (Lenz-law term)

```
S = S_M9.1'' + α ∫ T^μν J_μ J_ν √-g d⁴x
gdzie J_μ = ∂_μ Φ
```

- Self-consistent stress-energy back-reaction na Φ-EOM
- W strong-field T^TT genrowane przez własne Φ tworzy "rezystancję" jak Lenz law
- **Koszt strukturalny:** wprowadza non-linearity O(T²) w action
- **Risk:** unitarity, ghost modes, runaway solutions
- **Atrakcyjne:** mechanizm fizyczny dla self-suppression strong-field
  (zgodne z "back-reaction" intuicja substrate physics)

## 3. Plan testów OP-M92 (T-M92.1..T-M92.6)

| # | Test | Cel | Metoda | Status |
|---|------|-----|--------|--------|
| **T-M92.1** | Reverse-engineer scenario (e) action | Pokazać minimalną formę action principle która generuje f(ψ) = √(g_tt^GR/g_tt^TGP) | sympy: matching constraint + variational | **OPEN** |
| **T-M92.2** | Dual-field cost analysis | Cassini γ_PPN bound → maksymalne α dla σ field | analytic 1PN expansion + observational bounds | **OPEN** |
| **T-M92.3** | Conformal Ω(ψ) constraint | Equivalence principle bounds on Ω' | composition-dependent acceleration analysis | **OPEN** |
| **T-M92.4** | q-flow renormalization group | Czy q(ψ) flow jest consistent z RG fixed-point structure? | β-function sketch | **OPEN** |
| **T-M92.5** | Momentum back-reaction stability | Ghost-mode + unitarity analysis | linear stability + dispersion relation | **OPEN** |
| **T-M92.6** | Decision tree synthesis | Mapa kandydatów → ngEHT 2030/2032 outcomes | meta | **OPEN** |

**Note:** Wszystkie testy są **diagnostic/exploratory**. Nie ma "binary PASS"
verdict do M9.2 selection — wybór kandydata będzie data-driven po ngEHT verdict.
OP-M92 jest **scope document + readiness package**, nie definitive pivot.

## 4. Minimum requirements dla M9.2 (jakikolwiek kandydat)

Każdy M9.2 kandydat musi spełnić:

1. **Strong-field absorption:** generuje f(ψ_ph) → b_crit deviation ≤ 5%
   (within EHT envelope) lub ≤ 1% (within ngEHT 2030+ envelope).
2. **Weak-field PPN preservation:** γ_PPN = β_PPN = 1 at 1PN
   (Cassini |δγ| < 2.3e-5, LLR |δβ| < 8e-5) — zachowane M9.1'' P3 PASS.
3. **GW propagation:** c_GW = c_0 z OP-7 closure validity
   (LIGO/Virgo |c_GW/c_0 - 1| < 1e-15).
4. **Falsifiable signature:** distinct prediction od GR somewhere accessible
   (np. Sgr A* secondary photon ring n=2, gravitational lensing tail,
   pulsar timing signature).
5. **First-principles derivation:** wynika z natural action principle
   z TGP foundations (Z_2 symmetry, substrate physics, variational).

## 5. Decision tree — ngEHT 2030+ outcomes

```
                         ngEHT 2030-2032
                              |
              ┌───────────────┼───────────────┐
              |               |               |
    GR shadow at 1%   Intermediate band   TGP-like deviation
    (matches Schw)    (1-10% deviation)   (≥10% from GR)
              |               |               |
       M9.2 PIVOT      M9.1'' WEAKENED      M9.1'' VINDICATED
       MANDATORY       (refinement only)    (CONDITIONAL POSITIVE
              |               |              upgrades to POSITIVE)
              |               |
    Pick best M9.2    M9.1''' refinement
    candidate from    z modified ψ_ph
    OP-M92 candidates structure
              |
    Re-run all OP closures
    under M9.2 framework
```

**Threshold definitions:**
- "GR shadow at 1%" = b_crit_observed/b_crit_GR - 1 ∈ [-0.01, +0.01]
- "Intermediate" = ∈ [-0.10, -0.01] ∪ [+0.01, +0.10]
- "TGP-like" = ∈ [+0.10, +0.20] (M9.1'' predicts +0.1456)

**Probability assessment** (Bayesian, post OP-EHT-A):
- GR shadow at 1%: ~50% (consistent z całością GR-success story)
- Intermediate: ~30%
- TGP-like: ~20% (M9.1'' standalone success)

## 6. Timeline

| Year | Milestone | OP-M92 action |
|------|-----------|---------------|
| 2026 | OP-EHT-A NEGATIVE → OP-M92 OPEN | Scope document + diagnostic T1 |
| 2027 | OP-M92 candidates analysis | T-M92.2..T-M92.5 (theory readiness) |
| 2028 | OP-M92 decision tree finalized | T-M92.6 + paper draft (deferred section) |
| 2029 | EHT next-gen instrumentation | Monitor ngEHT progress |
| 2030 | ngEHT first results (Sgr A*) | Trigger 2 evaluation |
| 2031–2032 | ngEHT decisive verdict | M9.2 selection if mandated |
| 2033+ | M9.2 paper / M9.1'' vindication | Final closure |

**Critical insight:** OP-M92 *should* take ~2 years of part-time analysis
before ngEHT verdict. Aim is to have **all four candidates pre-analyzed**
so post-verdict response is selection + paper, not 6-month research scramble.

## 7. Pliki OP-M92 (planowane)

- `OP_M92_setup.md` — ten plik (scope + decision tree)
- `op_m92_T1_action_reverse_engineer.py` + `.txt` — T-M92.1
- `OP_M92_T1_results.md` — werdykt T-M92.1
- `op_m92_T2_dual_field_cost.py` — T-M92.2 (Brans-Dicke style bounds)
- `op_m92_T3_conformal_frame.py` — T-M92.3 (equivalence principle)
- `op_m92_T4_q_flow_RG.py` — T-M92.4 (RG analysis)
- `op_m92_T5_momentum_backreaction.py` — T-M92.5 (stability)
- `OP_M92_final_verdict.md` — synthesis (deferred until ngEHT)

## 8. Bottom line — OP-M92 is preparation, not pivot

OP-M92 jest **deferred conditional research program**:

- **Now (2026-04-25):** Open scope document, diagnostic T1.
  M9.1'' status: **alive standalone**, falsifiable in 4-6 lat.
- **2026-2030:** Analiza kandydatów, theory readiness.
- **2030+:** Trigger evaluation z ngEHT verdict.
- **If ngEHT confirms GR:** M9.2 pivot mandatory, candidates pre-analyzed.
- **If ngEHT confirms TGP-like:** M9.1'' vindicated, OP-M92 closes negative.

**Single most important deliverable:** decision tree z pre-analyzed kandydatami,
gotowy do execution w 2-4 tygodnie po ngEHT verdict, nie 6-12 miesięcy.

---

## 9. Cross-program implications

### F4 falsifiability (post-OP-EHT-A, post-OP-M92 OPEN)

OP-EHT closure F2: "TGP M9.1'' standalone predykuje +14.56% strong-field
shadow deviation. ngEHT 2030+ at 1% Sgr A* precision daje >10σ verdict."

OP-M92 status dodaje:
- **If verdict = GR:** M9.2 pivot z pre-analyzed kandydatami
  (2-4 weeks response, nie 6-12 months)
- **If verdict = TGP-like:** M9.1'' vindicated, OP-EHT upgrades to POSITIVE
- **No third option of "we'll figure it out then":** OP-M92 zapewnia
  programatyczny readiness niezależnie od outcome

### M9.1'' status — ostateczna pozycja

M9.1'' pozostaje **active baseline theory** z weakness boundary:
- Weak-field: PASS (Mercury, Cassini, LLR, OP-7 c_GW = c_0)
- Strong-field shadow: NEGATIVE on naive rescue (OP-EHT-A)
- Strong-field shadow: CONDITIONAL on ngEHT 2030+ verdict (OP-M92 prep)
- Cosmology: separate scope (OP-DESI, OP-Hubble — independent)

**Theory health = high** — falsifiable, predictive, locally consistent;
rescue prep w trybie deferred, nie crisis mode.
