---
title: "D01 anchor lock + propagacja parametrów (Φ₀, α_s, m_H, Σm_ν, g_0^e)"
date: 2026-05-06
cycle: D01
type: audit-resolution
status: ANCHOR LOCK MANIFEST + Phase 2 propagation (priorytet A: deterministic fixes; priorytet B: annotation; priorytet C: deferred)
parent: "[[../../audyt/D01_drifting_numbers/README.md]]"
predecessors:
  - "[[../op-newton-momentum/B3_v2_alphas_propagation_results.md]]" (B3-v2 closed 2026-05-01, α_s = 0.1184 14 lokacji)
  - "[[../op-L04-ODE-canonicalization-2026-05-04]]" (g_0^e formal lock pending)
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]" §B.3, §B.4, §C.9, §C.10
related:
  - "[[../../audyt/D01_drifting_numbers/README.md]]" (audit-source)
  - "[[../../audyt/D01_drifting_numbers/NEEDS.md]]"
tags:
  - TGP
  - D01
  - anchor-lock
  - propagation
  - parameter-drift
  - audit-resolution
tgp_status:
  folder_status: research
  level: L1
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: 2026-05-06
  may_edit_core: true
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: ["formal-Phi0-derivation", "tooling-Phi0-internal-unification"]
  depends_on:
    - "[[../op-newton-momentum]]" (B3-v2 anchor)
    - "[[../op-L04-ODE-canonicalization-2026-05-04]]" (α=2 canonical decision)
  impacts:
    - "[[../../audyt/D01_drifting_numbers]]"
  source_of_status:
    - "B3-v2 closed 2026-05-01 (α_s = 0.1184 lock + 14 lokacji propagated)"
    - "AUDYT 2026-05-01 §B.3, §B.4, §C.9, §C.10"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# D01 — anchor lock + propagacja parametrów

## Cel

Domknąć status D01 (liczbowy dryft α_s, m_H, Φ₀, Σm_ν, g_0^e) przez:

1. **LOCK MANIFEST** — kanoniczne wartości z formalnym rationale dla każdego parametru
2. **Inwentaryzacja** post-B3-v2 (~B3-v2 zamknął α_s; pozostałe parametry tj. m_H, Σm_ν wymagają propagacji)
3. **Propagacja Priorytet A** — deterministic fixes (clear errors/drift w core/, papers, README)
4. **Annotation Priorytet B** — papers_external z świadomą alt-formułą α_s = 0.1174 (N_f=5)
   otrzymują adnotację z canonical B3-v2 lock 0.1184
5. **Defer Priorytet C** — tooling scripts z internal Φ_0 = 24.65 (B3-v2 explicit
   deferred — łamanie self-consistency `K_geo·m_sp² = π·Φ₀²` etc.)

## Werdykt główny

> **D01 substantialnie ZAMKNIĘTY.**
>
> | Parametr | Status post-D01 |
> |----------|-----------------|
> | **α_s** | LOCKED 0.1184 (B3-v2 2026-05-01, 14 lokacji); +annotation w papers_external (10 wystąpień alt-formuły) |
> | **Φ₀** | LOCKED 24.783 (B3-v2 Brannen, 14 lokacji); 6 tooling scripts internal Φ_0=24.65 deferred |
> | **m_H** | LOCKED 125.31 (TGP-side, formal v×57/112) + 125.20 (PDG 2024 ref); core sync 8 wystąpień |
> | **Σm_ν** | LOCKED 59.01 meV (Z1 anchor B4 2026-05-01); README sync 2 wystąpień |
> | **g_0^e** | LOCKED 0.86941 (substratowa α=1 post-L04 cor:alpha1-preferred); formal derivation: deferred N1 |

## Lock manifest (kanoniczne wartości + rationale)

### α_s(M_Z) = 0.1184

**Formuła kanoniczna:** `α_s = N_c³·g_0^e / (8·Φ_0) = 27·0.86941/(8·24.783) = 0.1184`

**Rationale:**
- B3-v2 lock 2026-05-01 ([[../op-newton-momentum/B3_v2_alphas_propagation_results.md]])
- PDG comparison: 0.1180±0.0009 → 0.4σ (vs 1.16σ w pre-lock)
- 5/5 sympy LOCK PASS

**Alternatywna formuła (papers_external, N_f=5):**
`α_s = N_c³·g_0^e / (8·N_f²) = 27·0.86941/(8·25) = 0.1174`

Companion `tgp_companion.tex:321` explicit dokumentuje obie i preferuje canonical
(0.1184) — papers_external (10 wystąpień) wymagają analogicznej annotation.

### Φ₀ = 24.783 (Brannen vacuum)

**Formuła kanoniczna:** intrinsic z substratowego variational principle (post-G.0 closure 2026-05-02)

**Rationale:**
- B3-v2 anchor decision 2026-05-01: Brannen 24.783 vs Planck-derived 24.65 (Ω_Λ = 0.6847 → Φ₀ = 36·Ω_Λ)
- 24.783/36 = Ω_Λ^TGP = 0.6884 → +0.5σ od Planck (consistent)
- Brannen jest **parameter-free TGP prediction**; Planck-derived imports cosmological measurement (cyrkularność)

**Outstanding:** formal derivation Φ₀ = 24.78296... w sek09 lub sek08a (B3-v2 follow-up #2)

### m_H = 125.31 GeV (TGP-side) / 125.20 GeV (PDG 2024)

**Formuła kanoniczna TGP:** `m_H = v × 57/112 = 125.31 GeV`

z `v = 246.22 GeV` (VEV Higgs) i `57/112` z F11 mass formula.

**Rationale:**
- TGP-side jest **predykcją**, nie kalibracją — wartość derived z formal mass formula
- PDG 2024 reference: 125.20 ± 0.11 GeV
- Diff: 0.11 GeV ~ 1.0σ (PDG band: 0.11)
- Stary CW two-loop derivation (`sek09:1428` = 124, `dodatekU` = 124/125.1) jest **superseded** przez F11 v×57/112

### Σm_ν = 59.01 meV (Z1 anchor)

**Formuła kanoniczna:** Z1 bisection-derived anchor post-B4 (2026-05-01)

**Rationale:**
- B4 cycle [[../op-omicron1-sigmamnu-cosmo]] ustanowił 59.01 jako precyzyjny anchor
- README pre-B4 używa 59.6 (zeroth-order target)
- Diff: 0.6 meV (~1% w skali Σm_ν)

### g_0^e = 0.86941 (substratowa α=1)

**Formuła kanoniczna:** `g_0^e` z φ-fixed-point substratowego ODE (α=1, K_sub = K_geo·g²)

**Rationale:**
- L04 closure ([[../op-L04-ODE-canonicalization-2026-05-04]]) ustaliło że α=1 substratowa daje
  `r_21 = 206.77` (PDG match), `m_τ` z 83 ppm precyzji
- Sek08b `cor:alpha1-preferred` + `rem:formulation-dictionary` (12/12 PASS)
- 0.86941 5 cyfr znaczących — pochodzenie z numerical solver substratowego ODE

**Outstanding:** formal derivation analytical 0.86941 (L04 NEEDS N6, deferred)

## Inwentaryzacja po B3-v2

### α_s = 0.1184 (14 lokacji, B3-v2 2026-05-01 CLOSED)

Już zsynchronizowane:
- `core/sek00_summary:152-154`, `core/sek09_cechowanie:1048-1050,1067-1072,1291-1292`,
  `core/_meta_latex/status_map:1138`
- `tgp_letter.tex:39,143,160-170,290`, `tgp_companion.tex:300-322,1211`
- `tooling/scripts/color_tube_*.py` (3), `sin2thetaW_qcd_correction_tgp.py:310-311`,
  `ls8_prediction_taxonomy_audit.py:121-123,381`,
  `research/tgp_dependency_graph.py:194,240-242`,
  `research/graph_concept_flow.gexf:108,236,239,260`,
  `meta/PLAN_DOMKNIECIA_MASTER.md:207`

**Pozostałe (papers_external, alt-formuła N_f=5):**
- `papers_external/tgp_english_summary/main_en.tex:42`
- `papers_external/tgp_english_summary/sec08_predictions.tex:59,252,310`
- `papers_external/paper_lepton_masses/tgp_lepton_masses.tex:54,108,460,494,503,549,756`

→ **Action:** dodać annotation B3-v2 canonical lock (NON-BREAKING)

### Φ₀ = 24.783 (B3-v2 lock + 6 internal scripts deferred)

Już zsynchronizowane: 14 lokacji (jak α_s).

**Pozostałe (deferred):**
- `tooling/scripts/color_tube_advanced_tgp.py` internal Φ_0 = 24.6492
- `tooling/scripts/color_tube_variational_tgp.py` internal Φ_0 = 24.6492
- `tgp_chain_Phi0_to_masses.py`, `tgp_bridge_substrate_g0e.py`,
  `cosmo_frw_verification_v47.py`, `tgp_unified_predictions_v47.py`,
  `tgp_master_consistency_v47.py`, `tgp_prediction_taxonomy_v47.py` — internal Φ_0=24.65

→ **Action:** DEFERRED (B3-v2 explicit decision: `K_geo·m_sp² = π·Φ_0²` self-consistency
check would break) — tracked in NEEDS.md

### m_H = 125.31 (TGP) / 125.20 (PDG)

Już zsynchronizowane: README.md, companion (z C9 audit re-anchoring).

**Pozostałe (priorytet A):**
- `core/sek09_cechowanie/sek09_cechowanie.tex:1428` — `m_H = 124` (CW two-loop, superseded)
- `core/formalizm/dodatekU_su2_formalizacja.tex:19,197,295,308` — 124/125.1
- `core/formalizm/dodatekV_su3_formalizacja.tex:656` — 125.1
- `core/_meta_latex/status_map.tex:1459` — 125.1
- `tgp_companion.tex:520,1217` — PDG ref 125.25 (stale, should be 125.20 per C9)
- `core/sek00_summary:341` — TGP value ok (125.31), PDG sigma stale (0.3σ z PDG 125.25; should be 1.0σ z PDG 125.20)

→ **Action:** Priorytet A — propagacja edycje (8 wystąpień)

### Σm_ν = 59.01 meV (Z1)

**Pozostałe (priorytet A):**
- `README.md:100` — 59.6 meV (predictions list)
- `README.md:285` — 59.6 meV (predictions table)

→ **Action:** Priorytet A — README sync (2 wystąpień)

### g_0^e = 0.86941

Już zsynchronizowane (B3-v2 + L04). Brak otwartych edycji.

**Outstanding:** formal analytical derivation (L04 NEEDS N6).

## Plan propagacji (Phase 2)

### Phase 2A — Priorytet A (deterministic fixes, ~10 edycji)

| Plik | Linie | Zmiana | Charakter |
|------|-------|--------|-----------|
| `README.md` | 100 | `Σm_ν = 59.6 → 59.01` meV | numeric update |
| `README.md` | 285 | `Σm_ν 59.6 → 59.01` (table) | numeric update |
| `core/sek09_cechowanie/sek09_cechowanie.tex` | 1428 | `m_H = 124` → `m_H = 125.31` (TGP F11 v×57/112) + comment that older CW two-loop superseded | replacement + annotation |
| `core/formalizm/dodatekU_su2_formalizacja.tex` | 19, 197, 295, 308 | m_H 124/125.1 → 125.31 (TGP F11) | replacement w 4 wystąpieniach |
| `core/formalizm/dodatekV_su3_formalizacja.tex` | 656 | m_H 125.1 → 125.31 | replacement |
| `core/_meta_latex/status_map.tex` | 1459 | m_H 125.1 → 125.31 | replacement |
| `tgp_companion.tex` | 520 | PDG `125.25 ± 0.17` → `125.20 ± 0.11` (PDG 2024) + sigma update | replacement |
| `tgp_companion.tex` | 1217 | PDG `125.25 ± 0.17` → `125.20 ± 0.11` + sigma update | replacement |
| `core/sek00_summary/sek00_summary.tex` | 341 | sigma `0.3σ → 1.0σ` (z PDG 125.20 anchor) | sigma update |

### Phase 2B — Priorytet B (annotation papers_external, ~10 edycji)

10 wystąpień `α_s = 0.1174` w 4 plikach `papers_external/`. Strategia:

- **Główny plik per paper** (`main_en.tex` lub `tgp_lepton_masses.tex`): dodaje
  preamble comment block z annotation B3-v2 canonical lock
- **Pozostałe wystąpienia**: pozostają jako alt-formuła N_f=5 (świadoma decyzja
  autora, dokumentowana w companion:321)
- **NIE zmieniam wartości** — to są publikacje z alternatywną formułą

→ **Action:** dodaje **adnotację w preambule** lub przy pierwszym wystąpieniu α_s w każdym z 4 plików, że canonical TGP lock (B3-v2) daje 0.1184 z `N_c³·g_0^e/(8·Φ_0)` zamiast `N_c³·g_0^e/(8·N_f²)`.

### Phase 2C — Priorytet C (deferred)

Tooling scripts internal Φ_0 = 24.6492 — **NIE zmieniam** w tym cyklu:
- B3-v2 explicit decision: scripts używają internal Φ_0 = 24.65 dla
  self-consistency `K_geo·m_sp² = π·Φ_0²` i `κ = 3/(4·Φ_0)`
- Wymaga osobnego cyklu z full re-verification (estymata: 2-3 sesje)

→ **Tracked in NEEDS.md, future cycle `op-tooling-Phi0-unification`**

## Indeks plików

| Plik | Zakres |
|------|--------|
| [[README.md]] | (ten plik) — werdykt + lock manifest + plan |
| [[Phase0_balance.md]] | external inputs + structural axioms (CALIBRATION_PROTOCOL §2) |
| [[FINDINGS.md]] | eksportowalne wyniki |
| [[NEEDS.md]] | otwarte (Phase 2C deferred + formal derivations) |

## Cross-references

- [[../../audyt/D01_drifting_numbers/README.md]] — audit-source
- [[../../audyt/D01_drifting_numbers/POST_ACTION_UPDATE_2026-05-06.md]] — będzie utworzony
- [[../../audyt/PRIORITY_MATRIX.md]] — będzie zaktualizowany
- [[../op-newton-momentum/B3_v2_alphas_propagation_results.md]] — pre-existing α_s lock
- [[../op-L04-ODE-canonicalization-2026-05-04]] — α=1 substratowa decyzja (g_0^e=0.86941)
- [[../../meta/CALIBRATION_PROTOCOL.md]] — anti-overclaim protocol (binding)
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §B.3, §B.4, §C.9, §C.10
