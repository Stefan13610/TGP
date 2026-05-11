---
title: "M9 restructure note — co umarło, co przeżyło, jak czytać dziedzictwo M9.1''"
date: 2026-05-10
type: meta-protocol
status: 🟢 ACTIVE — methodological note (lightweight; cycle-by-cycle audit deferred)
binding_scope: "All cycles inheriting M9.1'' as canonical metric, post-2026-05-10"
related:
  - "[[PPN_AS_PROJECTION.md]] (siostrzany doc — PPN as projection, NIE physics)"
  - "[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] (siostrzany doc — anti-BD-drift)"
  - "[[CALIBRATION_PROTOCOL.md]]"
  - "[[../research/op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]] §3.1, §3.6, §5"
  - "[[../TGP_FOUNDATIONS.md]] §3.6 (emergent-metric block)"
parent: "[[README.md]]"
tags:
  - meta
  - methodology
  - M9-restructure
  - post-falsification-recovery
  - tier-hierarchy
  - anchor-not-canonical
---

# M9_RESTRUCTURE_NOTE — co umarło, co przeżyło, jak czytać dziedzictwo M9.1''

## §0 — Po co ten dokument

### §0.1 — Pytanie autora (2026-05-10)

> "Zastanawiam się czy M9 to w takim układzie śmieci, czy wystarczy przeorganizować program
> badawczy."
>
> — autor TGP, post-conversation o PPN-as-projection methodology

Odpowiedź: **głównie reorganizacja**. Ale honest accounting wymaga rozróżnienia, *które
elementy* M9 przeżyły, a które są martwe — bez tego rozróżnienia framework cierpi na
"M9.1''-as-canonical drift" analogiczny do BD-drift z `TGP_NATIVE_COMPUTATIONAL_PATTERNS`.

### §0.2 — Status tego dokumentu

**Lightweight methodology note**, NIE pełen audit cycle-by-cycle. Cykl-by-cykl audit (analog
BD-drift audit) **deferred** do dedykowanego cycle (np. `op-M9-restructure-audit-2026-05-XX`)
— pending decision after current critical-path stabilization (γ-id cascade closed; Branch A
re-asserted).

Ten doc:
- **DOES** klasyfikować elementy M9 (martwe / żywe / re-statusowane)
- **DOES** definiować Tier 0-4 hierarchy
- **DOES** rebrandować M9.1'' z "canonical metric" do "Path 2 anchor"
- **DOES** definiować reading rules dla inherited LOCKs (jak czytać starsze cykle)

- **DOES NOT** auditować poszczególne cykle (lista pending audit cycle)
- **DOES NOT** zmieniać sympy / classification / verdicts żadnego closed cyklu
- **DOES NOT** modyfikować PREDICTIONS_REGISTRY (już zaktualizowane per emergent-metric Phase_FINAL_close §12)

---

## §1 — Honest accounting: 4 buckets

### §1.1 — Bucket A: Martwe (śmieci, nie resuscytować)

| Element | Powód | Reading rule |
|---|---|---|
| **M9.1 original** (β_PPN=4 ansatz) | Falsyfikacja 3·10⁴σ (Mercury perihelion). Wrong functional form, brak recovery path. | Każde inheritance = ⚠️ pre-amendment artefact, do not propagate |
| **M9.1'' jako *jedyna* canonical metryka TGP** | Phase 1 emergent-metric pokazał 3-funkcyjny ansatz `{A(ψ), B(ψ), C(ψ)}`. M9.1'' to specific point z constraintem `A·B=1`, NIE general framework. | Reframing: M9.1'' jest *example point in {A,B,C} space*, NIE the canonical TGP metric |
| **M9.1'' single-source 2.5PN bez σ-coupling** | β_ppE=−15/4, sfalsyfikowane 5σ GWTC-3. Bez `c_0·κ_σ` recovery — point excluded. | Każde inheritance używające M9.1'' params dla ppE *must* explicit cite Path 2 (σ-coupling) recovery |

### §1.2 — Bucket B: Żywe bez zmian (rdzeń)

| Element | Survival reason |
|---|---|
| **M9.2 — Lenz back-reaction (momentum, inertia)** | m_inertial = ∫(∇ε_eq)², m_b=m_g automatic z S05. Phase 5 emergent-metric (10/10 PASS). **Niezależne od specyficznej formy metryki.** Działa dla każdego punktu w {A,B,C}. |
| **Newton I + II structural derivation** | Galilean-translated δΦ_eq + linear back-reaction `F_BR = −m·a`. Niezależne od specyficznej f(ψ). |
| **WEP / equivalence principle** | Z postulatu → automatyczna konsekwencja S05. Wzmocnienie, nie słabość. |
| **PPN constraints w 1PN/2PN** (`b_1=−a_1`, `ξ_2=...`) | Native results dla CAŁEJ rodziny {A,B,C}, NIE specyficznie M9.1''. Per `meta/PPN_AS_PROJECTION.md` reframe — to constraint na native Taylor coefs, projektujący się na PPN basis. |
| **§6 FOUNDATIONS Lenz paradigm** | Foundational, applies regardless of metric specifics. |

### §1.3 — Bucket C: Re-statusowane (z "fundamentalne" do "specyficzny anchor")

| Element | Pre-falsification status | Post-emergent-metric status |
|---|---|---|
| **M9.1'' params** `{a_1^M911, a_2^M911, a_3^M911, b_2^M911, ξ_3^M911}` | "Canonical metryka TGP" | **Path 2 anchor** — specific point in {A,B,C} space, structurally preferred (preserves 3/3 SU(2) paths z Phase 6 cross-consistency), conditional na σ-coupling `c_0·κ_σ = 4/3` żeby pass GWTC-3 |
| **γ_PPN = β_PPN = 1 z M9.1''** | "M9.1'' daje GR limit at 1PN" | **Native result całej rodziny {A,B,C}** z constraintem `b_1=−a_1`, `ξ_2=...`. M9.1'' to jeden punkt spełniający te constraints — nie unique source |
| **Schwarzschild-like 1/r asymptotics z M9.1''** | "Confirms M9.1'' jako correct GR limit" | One sample point's asymptotic; full {A,B,C} family asymptotic landscape pending dedicated cycle |

### §1.4 — Bucket D: Wymaga audytu (cykle inheriting M9.1'')

**Status: list pending**, audit deferred do dedykowanego `op-M9-restructure-audit` cycle.

Każdy cykl pre-2026-05-09 inheriting M9.1'' będzie classified per:
- **(a)** inherits M9.1'' specific form jako *example* w {A,B,C} space → ✅ OK, no change
- **(b)** inherits M9.1'' specific values jako Path 2 anchor → ✅ OK z explicit conditional citation tego doc
- **(c)** inherits M9.1'' jako "the canonical TGP metric" → ⚠️ requires reframe annotation per §3 below
- **(d)** depends on M9.1'' being unique solution → 💀 needs reopening lub closure

Audit subagent template per CALIBRATION_PROTOCOL §4.4 (analog BD-drift audit, dispatched
podczas dedicated audit cycle).

---

## §2 — Tier hierarchy (post-cascade structural picture)

```
┌─────────────────────────────────────────────────────────────────────┐
│ Tier 0 — Axioms (BINDING, untouched by falsification)              │
│   S05 single-field, §5.1 czym TGP NIE jest, §6 Lenz paradigm       │
├─────────────────────────────────────────────────────────────────────┤
│ Tier 1 — Frameworks (BINDING, post-2026-05-10)                     │
│   • Emergent-metric {A(ψ), B(ψ), C(ψ)} (op-emergent-metric Phase 1) │
│   • PPN-as-projection (meta/PPN_AS_PROJECTION.md)                   │
│   • TGP-native computational patterns (meta/TGP_NATIVE_*)           │
│   • Cycle lifecycle + calibration protocol (meta/CALIBRATION_*)     │
├─────────────────────────────────────────────────────────────────────┤
│ Tier 2 — Anchors (CONDITIONAL — specific points/values)             │
│   • M9.1'' params jako Path 2 anchor (conditional na c_0·κ_σ=4/3)   │
│   • c_0 = 4π (op-c0-derivation, heuristic; first-principles pending)│
│   • κ_σ = 1/(3π) (op-kappa-sigma, heuristic)                        │
│   • γ ~ M_Pl² (GF.B-STRUCTURAL, β=γ vacuum stability OPEN)          │
│   • Pattern 2.5: m_Φ_observable env-dependent (BINDING-PRINCIPLE)   │
├─────────────────────────────────────────────────────────────────────┤
│ Tier 3 — Parameters (constrained by observables)                    │
│   Native Taylor coefs: a_n, b_n, ξ_n (n=1,2,3,...)                  │
│   Coupling coefs: c_0, κ_σ                                           │
│   Forced-zero: α_i (i=1,2,3), ζ_i (i=1,2,3,4) — z substrate symmetry│
├─────────────────────────────────────────────────────────────────────┤
│ Tier 4 — Observables (measurement targets)                          │
│   1PN/2PN: deflekcja, Shapiro, perihelion, Nordtvedt, range         │
│   2.5PN+ GW: β_ppE, α_ppE, h_TT^σ, dispersion, c_GW                 │
│   Cosmology: H(z), f_growth, σ_8, w_eff (cross-cycle pending)       │
└─────────────────────────────────────────────────────────────────────┘
```

**Reading direction:** observables (Tier 4) constrain parameters (Tier 3); parameters live
in space defined by frameworks (Tier 1); anchors (Tier 2) są specyficzne wybory w tej
przestrzeni; axioms (Tier 0) gwarantują strukturalne identities (forced-zero PPN params,
EP automatic, etc.).

**M9.1'' lokalizacja w tej hierarchii:** Tier 2 anchor (specific values of Tier 3 params),
NIE Tier 1 framework. To jest istotny demoton statusowy, ale nie deletion.

---

## §3 — Rebranding: M9.1'' → "Path 2 anchor"

### §3.1 — Co się zmienia w naming convention

**Old usage (pre-2026-05-09):**
> "TGP postuluje metrykę M9.1'': `g_tt = -c₀²(4-3ψ)/ψ`"

**New usage (post-emergent-metric closure + this doc):**
> "TGP wybiera Path 2 anchor (M9.1'' params + σ-coupling) jako structurally preferred punkt
> w {A(ψ), B(ψ), C(ψ)} family. M9.1'' label preserved as historical/literature continuity
> reference."

### §3.2 — Kiedy używać starej vs nowej terminologii

- **Historical / literature continuity:** "M9.1''" preserved (old papers, archived cycles,
  cross-references)
- **Active framework presentation:** "Path 2 anchor (M9.1'' params)" — explicit conditional
- **Foundational layer (FOUNDATIONS, registry):** Tier 1 = `{A,B,C}` family; Tier 2 = Path 2
  anchor; "M9.1'' specific form" = exemplar
- **Falsification context:** "M9.1'' specific point in `{A,B,C}` space, observationally
  excluded without σ-coupling" — clear scope

### §3.3 — Co rebranding NIE zmienia

- File names i directory names cykli `op-...-M911-...` preserved (no rename churn)
- Sympy results, classification verdicts, P-requirements status — unchanged
- PREDICTIONS_REGISTRY M911-P1/P2/P3 entries — already updated per emergent-metric Phase_FINAL_close §12 ("FALSIFIED specific + RECOVERY EXISTS parametric family")
- Historical references w archived cycles — preserved as-is

---

## §4 — Reading rules for inherited LOCKs

Każdy cykl post-2026-05-10 czytający LOCK z M9.1'' MUSI explicit citować, w której z
4 kategorii (Bucket D) ten LOCK się znajduje:

```markdown
**M9.1'' inheritance status (per meta/M9_RESTRUCTURE_NOTE.md §1.4):**
- (a) Example point in {A,B,C} space → no change needed
- (b) Path 2 anchor citation → conditional on c_0·κ_σ=4/3
- (c) "Canonical TGP metric" framing → REFRAME per §3.2 above
- (d) M9.1'' as unique solution → ASK-RULE Trigger E (M9-as-canonical drift)
```

**ASK-RULE Trigger E** (proposed extension TGP_NATIVE_COMPUTATIONAL_PATTERNS §1):

> If inherited LOCK uses M9.1'' as the canonical / unique TGP metric (category c or d above),
> pause Phase 1 sympy and ask user: "ten cykl dziedziczy M9.1'' jako canonical metric.
> Per meta/M9_RESTRUCTURE_NOTE.md §3, M9.1'' jest Path 2 anchor, NIE framework. Czy chcesz:
> (i) reframe inheritance jako Path 2 anchor w {A,B,C} space, (ii) audit cyklu pod kątem
> czy dependency na M9.1''-as-canonical jest legitimate, (iii) close cycle if dependency
> jest fatal?"

(Formal addition Trigger E to ASK-RULE deferred to dedicated cycle TGP_NATIVE_COMPUTATIONAL_PATTERNS edit.)

---

## §5 — Implication: framework siła rośnie, nie maleje

Restruktura jest **wzmocnieniem**, nie osłabieniem:

| Aspect | Pre-restructure | Post-restructure |
|---|---|---|
| Foundational metric | M9.1'' jako single canonical form | {A,B,C} family — strukturalnie general |
| Falsification response | "TGP wykluczone 5σ" (catastrophic) | "M9.1'' point excluded; Path 2 (σ-coupling) recovery; neighbourhood otwarte" |
| Parameter freedom | Implicit ~10 PPN params + M9.1'' values | Native ~5-7 Taylor coefs + 4+4 forced-zero PPN params (honest count) |
| Methodology | Postulate metric form, check PPN | Native obserwable z Φ-EOM, PPN/ppE jako L2 projection |
| Cross-level unification | Postulate per level | H6.1 unified mechanism (Phase 6 emergent-metric: g_eff at level 2 ↔ SU(2) at level 3) |

**TGP po restruktrze jest *bardziej* spójny strukturalnie**, nie mniej. Falsyfikacja M9.1''
nie była "TGP fails", tylko "TGP framework jest bardziej ogólny niż początkowo zakładano —
M9.1'' to jeden punkt w przestrzeni, GWTC-3 wyklucza ten punkt bez σ-coupling, recovery
exists strukturalnie".

---

## §6 — Pending items (deferred do dedykowanych cykli)

| # | Item | Rationale | Estimated effort |
|---|---|---|---|
| **1** | `op-M9-restructure-audit-2026-05-XX` cycle | Cycle-by-cycle audit Bucket D classification (a/b/c/d per §1.4) | 3-5 sesji |
| **2** | TGP_FOUNDATIONS § revision: M9.1'' demoted z framework do Tier 2 anchor | Multi-section edit (§3.5, §3.6, §3.7) | 1-2 sesje |
| **3** | TGP_NATIVE_COMPUTATIONAL_PATTERNS §1 ASK-RULE Trigger E formal addition | Per §4 above; small edit | <1 sesja |
| **4** | PREDICTIONS_REGISTRY cross-check post-restructure | Verify M911-P entries still consistent z new tier framing | <1 sesja |
| **5** | audyt/T01_LIGO3G_falsifier reactivation jako native-coefs falsifier (Bucket B/C joint) | Per `meta/PPN_AS_PROJECTION.md §6` pending; complementary work | 2-3 sesji |

**None are blockers** dla active critical path. Można robić incrementally.

---

## §7 — Sign-off

**Doc authored:** 2026-05-10 (post-conversation autor + Claudian o M9 śmieci-vs-reorganizacja).

**Status:** ACTIVE methodological note. Lightweight — cycle-by-cycle audit deferred.

**Insight credit:** autor TGP, 2026-05-10. Doc formalizuje observation jako tier hierarchy
+ rebranding convention.

**Cross-references:**
- `meta/PPN_AS_PROJECTION.md` — siostrzany doc (PPN as projection)
- `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` — siostrzany doc (anti-BD-drift)
- `research/op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md` §3, §5 —
  emergent-metric closure źródło Tier 1 framework
- `TGP_FOUNDATIONS.md` §3.6 — current emergent-metric block (pending revision per item #2)

---

**Three sister methodology docs (2026-05-10 trio):**

1. `TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` — anti-**BD**-drift (BD-form ≠ BD-meaning)
2. `PPN_AS_PROJECTION.md` — anti-**projection**-confusion (PPN/ppE = chart, NIE physics)
3. `M9_RESTRUCTURE_NOTE.md` (this doc) — anti-**M9.1''-as-canonical**-drift (M9.1'' = anchor, NIE framework)

Razem definiują: jak nie pomylić TGP-native physics z (a) standardowym BD/scalar-tensor
toolingiem, (b) PPN/ppE chart parametrami, (c) historical specific anchors. Każdy z tych
trzech typów drift'u ma inny mechanism i inny remediation path.
