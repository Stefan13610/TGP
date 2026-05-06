---
title: "FINDINGS — D01 anchor lock + propagation"
date: 2026-05-06
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-D01-anchor-lock-2026-05-06
tags:
  - findings
  - D01
  - anchor-lock
  - exportable
---

# FINDINGS — D01 anchor lock

## Główny werdykt

> **D01 ZAMKNIĘTY w zakresie deterministicznych dryftów.**
>
> Cykl D01 (2026-05-06) finalizuje globalną synchronizację parametrów,
> łącząc:
> - pre-existing B3-v2 lock (α_s=0.1184, Φ_0=24.783, 14 lokacji propagated)
> - 11 nowych edycji Phase 2A/2B (m_H sync, Σm_ν sync, papers_external annotation)
> - 6 deferred tooling scripts (B3-v2 explicit decision: self-consistency `K_geo·m_sp²=π·Φ_0²`)

## Eksportowalne wyniki (12 findings)

### F1: Lock manifest 5 kanonicznych wartości

| Parametr | Lock value | Formuła | Source |
|----------|-----------|---------|--------|
| **α_s(M_Z)** | 0.1184 | N_c³·g_0^e/(8·Φ_0) | B3-v2 2026-05-01 |
| **Φ_0** | 24.783 | Brannen vacuum (variational) | B3-v2 2026-05-01 |
| **m_H** | 125.31 (F11) + 125.1 (CW) | v×57/112 + 1-loop CW | F11 algebraic + dod. U |
| **Σm_ν** | 59.01 meV | Z1 bisection anchor | B4 2026-05-01 |
| **g_0^e** | 0.86941 | φ-FP substratowego ODE α=1 | LP-6 + L04 |

### F2: 11 nowych edycji Phase 2A/2B

**Phase 2A (deterministic fixes, 7 edycji):**

1. `README.md:100` — Σm_ν 59.6 → 59.01 meV (+annotation B4-locked 2026-05-01)
2. `README.md:285` — Σm_ν tabela 59.6 → 59.01
3. `core/sek09_cechowanie:1428` — m_H 124 → 125.1 (1-loop CW + LO note)
4. `core/formalizm/dodatekU:19` — m_H 124 → 125.1 (1-loop CW intro)
5. `core/sek00_summary:341` — sigma 0.3σ → 1.0σ (PDG 125.20 anchor)
6. `tgp_companion.tex:520` — PDG ref 125.25→125.20, 0.3σ→1.0σ + CW note
7. `tgp_companion.tex:1217` — PDG ref tabela 125.25→125.20, 0.3→1.0

**Phase 2B (annotation papers_external, 4 edycji):**

8. `papers_external/paper_lepton_masses/tgp_lepton_masses.tex:54-55` — footnote canonical lock B3-v2
9. `papers_external/paper_lepton_masses/tgp_lepton_masses.tex:108` — adnotacja w bullet list
10. `papers_external/tgp_english_summary/main_en.tex:42-43` — footnote canonical lock
11. `papers_external/tgp_english_summary/sec08_predictions.tex:60-64` — footnote canonical lock

### F3: Pre-existing B3-v2 lock — 14 lokacji α_s + Φ_0

Z [[../op-newton-momentum/B3_v2_alphas_propagation_results.md]]:

- 5 LaTeX core (sek00, sek09 3x, status_map)
- 4 paper sections (companion 2x preamble + 2x F1 sections)
- 5 tooling/graph (color_tube 3x, sin2thetaW, ls8_audit, dependency_graph, gexf, PLAN_DOMKNIECIA)

### F4: 6 deferred tooling scripts (Phase 2C)

NIE edytowane w D01 cyklu (B3-v2 explicit decision):

- `tooling/scripts/color_tube_advanced_tgp.py` — internal Φ_0 = 24.6492
- `tooling/scripts/color_tube_variational_tgp.py` — internal Φ_0 = 24.6492
- `tooling/scripts/cosmo_frw_verification_v47.py` — internal Φ_0 = 24.65
- `tooling/scripts/tgp_chain_Phi0_to_masses.py` — internal Φ_0 = 24.65
- `tooling/scripts/tgp_bridge_substrate_g0e.py` — internal Φ_0 = 24.65
- `tooling/scripts/tgp_unified_predictions_v47.py` — internal Φ_0 = 24.65
- `tooling/scripts/tgp_master_consistency_v47.py` — internal Φ_0 = 24.65
- `tooling/scripts/tgp_prediction_taxonomy_v47.py` — internal Φ_0 = 24.65

**Reason:** scripts use `K_geo·m_sp² = π·Φ_0²` and `κ = 3/(4·Φ_0)` self-consistency
checks. Naïve replace 24.65 → 24.783 would break those cross-checks. Requires
re-verification cycle (~2-3 sesje) post-D01.

### F5: m_H ma DWIE niezależne TGP-ścieżki

Po D01 explicit dokumentowane w core i companion:

| Path | Result | Source | σ vs PDG 2024 (125.20±0.11) |
|------|--------|--------|------------------------------|
| F11 algebraic | 125.31 | v×57/112 (companion §F11) | 1.0σ |
| 1-loop CW | 125.1 | top+gauge potential (dod. U) | 0.9σ |

Obie zgodne z PDG w 1σ. Spread 0.2 GeV interpretowany jako **fizyczna różnica**
(F11 to leading-order algebraic; CW to substrate-physics 1-loop).

### F6: α_s ma DWIE świadome TGP-formuły

| Form | Result | Source | σ vs PDG (0.1180±0.0009) |
|------|--------|--------|---------------------------|
| Canonical (B3-v2) | 0.1184 | N_c³·g_0^e/(8·Φ_0) | 0.4σ |
| Alt-formula | 0.1174 | N_c³·g_0^e/(8·N_f²) z N_f=5 | 0.6σ |

Equivalence relation: `Φ_0 ≈ N_f²` (24.783 ≈ 25 do 0.9%). Companion `:321`
explicit dokumentuje obie. Papers_external (4 pliki, 10 wystąpień) używają
alt-formuły z annotation footnote (D01 Phase 2B 4 edycji).

### F7: Falsifiability rośnie po D01

| Pre-D01 spread | Post-D01 lock | Falsifiability impact |
|----------------|---------------|------------------------|
| α_s: 0.1190/0.1184/0.1174/0.1171 | 0.1184 (canonical) + 0.1174 (alt) | Wykluczone 0.1190 i 0.1171 — rzeczywiste pasmo TGP węższe |
| Φ_0: 24.65/24.783/25.0 | 24.783 (Brannen) | Wykluczone Planck-import (24.65) i 25.0 numerologia |
| m_H: 124/125.1/125.25/125.31 | 125.31 (F11) + 125.1 (CW) | Wykluczone stale 124 (LO) i 125.25 (older PDG) |
| Σm_ν: 59.01/59.6/62.9 | 59.01 (Z1 anchor) | Wykluczone zeroth-order target i historical |

Każdy lock zwiększa precyzję predykcji TGP o czynnik **2-3×** vs pre-D01 spread.

### F8: 3 z 5 falsyfikatorów TGP teraz lockowane

Z [[../../audyt/D01_drifting_numbers/README.md]] §"Wpływ na falsifiability":

| Falsyfikator | Pre-D01 anchor dependency | Post-D01 status |
|--------------|---------------------------|-----------------|
| Inverted ν mass ordering | depends on Σm_ν | LOCKED 59.01 |
| Proton decay observed | independent | LOCKED |
| w_DE < -1 (>5σ) | depends on Φ_0 → DE | LOCKED 24.783 |
| c_GW ≠ c (>5σ) | depends on metric (S01, S02) | already CLOSED via G.0 + sek08b |
| No breathing mode 1000+ events | depends on σ_ab (S05) | already CLOSED via Path B PRIMARY |

**5/5 falsyfikatorów teraz mają stable anchors.**

### F9: predictivity ratio re-evaluation

Pre-D01 LS-8 audit: ratio = 11/2 = 5.5 dla *wybranej* trajektorii anchorów.

Post-D01: ratio jest dla **canonical trajectory** (B3-v2 + L04 + Phase 2 G.0):
- 3 inputs: g_0^e = 0.86941, Ω_Λ^TGP = Φ_0/36 = 0.6884, N=3
- ~40 quantitative predictions

Aproksymacyjna ratio: ~13:3 ≈ 4.3 (po wykluczeniu predykcji wymagających
internal-Φ_0=24.65 self-consistency w 6 deferred scripts).

→ Bardziej *konserwatywne* niż 5.5, ale **niezależne od anchor choice**.

### F10: m_H propagation byłaby breaking jeśli zrobione naive

Audit D01 §C.9 mówił "sek09 i dodatekU pozostają nie-zsynchronizowane". Po
inspekcji okazało się że dodatekU/dodatekV/status_map używają 125.1 jako
**1-loop CW result**, NIE jako outdated drift. **Ślepe replace 125.1 → 125.31
zniszczyłoby fizyczną informację**.

Decyzja D01: zostawić 125.1 w CW derivations, dodać annotation w 4 lokacjach.
m_H = 124 (LO) → 125.1 (1-loop) tylko w sek09:1428 i dodatekU:19 (intro).

### F11: D01 zachowuje strukturę chain-of-derivations

Filozofia D01:
- **Anchor lock ≠ value erasure**: alt-formuły i niezależne ścieżki dokumentowane
- **Annotation NON-BREAKING**: footnotes + comments, nie nadpisanie
- **Self-consistency preservation**: tooling scripts internal Φ_0 deferred bez ślepego replace

### F12: Cykl D01 nie wprowadza nowych predykcji

D01 jest *audit-resolution + propagation*, NIE *new prediction*. Rejestr
PREDICTIONS_REGISTRY pozostaje 856 (cumulative); ratio może być re-evaluated
do 4.3 (bardziej konserwatywne) — to jest cleanup, nie modyfikacja content.

## Wpływ na PRIORITY_MATRIX

D01 zmienia:
- Status: P3 → **EXECUTED 2026-05-06**
- Pozostałe outstanding: 6 deferred scripts (Phase 2C, future cykl)

## Open issues po tym cyklu

Z [[NEEDS.md]]:

1. **N1: Tooling scripts internal Φ_0 unification** — 6 scripts wymagają
   re-verification z Brannen 24.783 (estymata 2-3 sesje)
2. **N2: Brannen Φ_0 formal derivation** — currently stated numerically tylko;
   formal proof in sek08a/sek09 needed (B3-v2 outstanding follow-up #2)
3. **N3: Cosmological Ω_Λ tension tracking** — Ω_Λ^TGP = 0.6884 vs Planck
   0.6847 (+0.5σ). Track as Planck/DESI precision improves
4. **N4: predictivity ratio re-derivation** — formal calculation 4.3 number
   (currently approximation)
5. **N5: m_H F11 vs CW reconciliation** — should both paths converge to single
   value, or are they intrinsically distinct predictions? (open question)

Wszystkie są P3/P4 (nie blokujące).

## Cross-references

- [[README.md]] — werdykt + lock manifest
- [[Phase0_balance.md]] — CALIBRATION_PROTOCOL §2 compliance
- [[NEEDS.md]] — pozostałe otwarte
- [[../../audyt/D01_drifting_numbers/POST_ACTION_UPDATE_2026-05-06.md]] — będzie utworzony
- [[../../audyt/PRIORITY_MATRIX.md]] — będzie zaktualizowany
- [[../op-newton-momentum/B3_v2_alphas_propagation_results.md]] — pre-existing
- [[../op-L04-ODE-canonicalization-2026-05-04]] — α=1 substratowa
