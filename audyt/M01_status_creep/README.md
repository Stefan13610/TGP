---
title: "M01 — status creep w PREDICTIONS_REGISTRY"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/M01_status_creep
tags:
  - audit
  - methodology
  - status-creep
  - predictions-registry
  - falsifiability
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../M02_ledger_pollution]]"
  - "[[../S06_circular_anchors]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
tgp_status:
  folder_status: audit
  level: L4
  kind: audit
  core_compatibility: stale
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["registry-audit-status-creep"]
  depends_on: []
  impacts:
    - "[[../M02_ledger_pollution]]"
    - "[[../M03_balance_sheet_missing]]"
  source_of_status:
    - "[[../../PREDICTIONS_REGISTRY.md]] head 80 (multiple promotions)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# M01 — status creep w PREDICTIONS_REGISTRY

## Klasa: CHAOS METODOLOGICZNY • Priorytet: P3

## Diagnoza

[[../../PREDICTIONS_REGISTRY.md]] head 80 demonstruje wzorzec
**multi-step status promotion** bez explicit gate sprawdzającego
falsifiability na każdym kroku.

### Łańcuchy promocji

Z head 80:

> **Status cascade activated:** ξ.1 PARTIALLY DERIVED (refined) → DERIVED (refined²);
> XS.1 PARTIALLY DERIVED (refined) → DERIVED (refined²); UV7 STRUCTURAL-DERIVED → DERIVED.

> **HH5 promoted STRUCTURAL HINT → PARTIALLY DERIVED** (mixing-operator framework);
> η.1 cascade comprehensive **DERIVED (refined²)**

> **HH5/KK5 promoted PARTIALLY DERIVED → DERIVED** (mixing-operator framework
> structurally complete).

> **F6 STRUCTURAL → DERIVED upgrade** (cytuje chi.1 jako źródło)

To są **multi-stage promotions**: `STRUCTURAL HINT → PARTIALLY DERIVED →
DERIVED → DERIVED (refined) → DERIVED (refined²)`.

## Problem strukturalny

Każda promocja powinna mieć dedicated gate:

- **STRUCTURAL HINT** (zauważona koincydencja, brak derywacji)
- **STRUCTURAL** (formalna struktura, ale brak liczb)
- **PARTIALLY DERIVED** (formula z parametrami fitted)
- **DERIVED** (formula bez free params, daje obserwowalne)
- **DERIVED (refined)** (post-correction patches)
- **DERIVED (refined²)** (cross-channel locks)

Ale rzeczywiste promocje w rejestrze:

1. Cykle 5+7+6 sub-tests PASS (lokalne)
2. Bez globalnego algebraicznego sprawdzenia czy structural axioms
   (g*, N_A, K_struct, ...) faktycznie sprzęgają z outputs
3. Gate definiowany jako *score-based*, nie *falsifiability-based*

### Konkretny przykład — ξ.1 promocja

ξ.1 (`research/op-xi1-photon-ring/`):

- 2026-04-29: ξ.1 program END (FULL CONVERGENCE 18/18)
- promocja: STRUCTURAL → DERIVED → DERIVED (refined) → DERIVED (refined²)
- każda promocja triggered przez kolejny zewnętrzny cykl (mini-cycle ε.1,
  ζ.1, θ.1, η.1, α.1, η.2, κ.1, ι.1, μ.1)
- **brak independent test**, że TGP value ξ_grav ≈ 1.06 (lub similar)
  jest konsistent z observed BH shadow

## Późniejsze downgrade'y — przyznanie błędu

Audit § J wymienia downgrade'y:

| Cykl | Original | Audit-aware | Powód |
|------|----------|-------------|-------|
| DE1/DE2 | LOCKED w=-1 | LIVE TENSION | DESI 2024 ~3σ evolving DE |
| ω.2 | LOCKED | LIVE PARTIAL | Δχ²=0.21 << 9 |
| M9 | KOMPLETNY | 1PN-CLOSED | 2PN/strong-field/GW-tensor OPEN |
| TT13-TT18 (ψ.1.v1) | „WYKONALNY DZIŚ" | WITHDRAWN | varying-α not varying-c |
| ZZ2/ZZ3 (ω.3) | LOCKED | LOCKED-CONDITIONAL | inherits ω.2 LIVE PARTIAL |
| F6 (post-χ.1 audit) | DERIVED | (rollback recommended) | tautologia |

To jest **explicit acknowledgment**, że status creep istniał. Liczba
faktycznych downgrade'ów wskazuje, że proces promocji nie miał
falsifiability gate.

## Wzorzec

[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] § 3 identyfikuje
3 incydenty wzorca *systemic over-claiming* (λ.1, χ.1, UV.2). Status
creep to **ten sam wzorzec na poziomie meta-rejestru**:

- Sub-tests patrzą lokalnie; żaden nie sprawdza globalnego
  inputs↔outputs balance
- Promocje są score-based (X/Y PASS rate)
- Przy braku globalnego balance sheet, każdy successful sub-test
  uzasadnia kolejną promocję

## Wpływ na predyktywność

LS-8 audit deklaruje:

> Predictivity ratio = 11/2 = 5.5. Jeśli A=a_Γ/φ potwierdzone → 13/2 = 6.5.

Liczba 5.5 (lub 6.5) zakłada, że wszystkie wpisy DERIVED są faktycznie
DERIVED. Po downgrade'ach (DE1/DE2, ω.2, M9, ψ.1.v1) i potencjalnych
rollback'ach (F6, UV.2, ω.3) faktyczna predyktywność jest **niższa**.

Z punktu widzenia external referee: ratio 5.5 jest *headline number*
w README; rzeczywista wartość wymaga retrospektywnego audytu wszystkich
DERIVED entries.

## Status w audycie

[[../../meta/AUDYT_TGP_2026-05-01.md]] dokumentuje pojedyncze downgrade'y
(§J, §K, §L, §M, §N, §O), ale nie identyfikuje wzorca status creep jako
osobnego problemu. To jest implicit acknowledgment przez liczbę
downgrade'ów.

## Rekomendacja

Otworzyć dedykowany cykl `op-registry-audit-status-creep/`:

### Phase 1 — falsifiability gate retrofit

Dla każdego DERIVED wpisu w PREDICTIONS_REGISTRY (~80 wpisów):

1. Zastosować [[../../meta/CALIBRATION_PROTOCOL.md]] balance sheet:
   - **External inputs**: PDG, CODATA, observational
   - **Structural axioms**: które są niezależnie LOCKED
   - **Derived outputs**: co cykl twierdzi że wyprowadza
   - **Tautology test**: czy outputs wyrażalne wyłącznie przez external
     inputs i axioms?
   - **Falsifiability test**: czy istnieje wartość axiom która
     wykluczyłaby match?

2. Sklasyfikować każdy wpis:
   - **DERIVED FULL** — pass tautology + falsifiability
   - **DERIVED CONDITIONAL** — passes z explicit dependency
   - **NUMEROLOGICAL** — fitted, nie derived (jak χ.1, UV.2)
   - **STRUCTURAL** — formal struktura, brak liczb
   - **STRUCTURAL HINT** — koincydencja

### Phase 2 — registry refactor

Update PREDICTIONS_REGISTRY:

- Każdy entry z explicit gate-passing source
- Counter rozdzielony na DERIVED FULL / CONDITIONAL / NUMEROLOGICAL
- Falsifiability column dla każdego entry

### Phase 3 — re-derive predictivity ratio

Z post-Phase 2 classification, re-derive predictivity. Spodziewane:

- Niższe ratio (np. 4–5 zamiast 5.5)
- Wyższa wiarygodność (każda DERIVED FULL ma falsifiability gate)

**Estymata:** 4–6 tygodni (gate retrofit dla 80 entries) + 1 tydzień
refactor.

## Pliki dotknięte

- `PREDICTIONS_REGISTRY.md` — full audit, retrofit, refactor
- `INDEX.md` — counter update
- `LS-8 prediction taxonomy audit` — re-run z nową klasyfikacją

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §IV.M1
- [[../PRIORITY_MATRIX.md]] klaster C
- [[../../PREDICTIONS_REGISTRY.md]]
- [[../../meta/CALIBRATION_PROTOCOL.md]]
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §J, §K, §L (downgrade history)
- [[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §3 (wzorzec systemic)
- [[../M02_ledger_pollution]] (powiązany — counter pollution)
- [[../M03_balance_sheet_missing]] (powiązany — pre-74394a8 cykle bez balance sheet)
