---
title: "NEEDS — D01 anchor lock (otwarte problemy)"
date: 2026-05-06
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-D01-anchor-lock-2026-05-06
tags:
  - needs
  - D01
  - tooling-Phi0
  - formal-derivation
  - predictivity-ratio
---

# NEEDS — D01 (otwarte problemy)

## Status

D01 substantialnie ZAMKNIĘTY w zakresie deterministicznych dryftów po cyklu
2026-05-06. Pozostają outstanding follow-ups (P3/P4) — głównie meta i formal.

## Otwarte luki

### N1: Tooling scripts internal Φ_0 = 24.65 unification

**Status:** OPEN — P3 (deferred from B3-v2 2026-05-01).

**Problem:** 6 tooling scripts używają internal `Φ_0 = 24.6492` (Planck-derived
z `36·Ω_Λ`) zamiast canonical Brannen `Φ_0 = 24.783`:

- `tooling/scripts/color_tube_advanced_tgp.py`
- `tooling/scripts/color_tube_variational_tgp.py`
- `tooling/scripts/cosmo_frw_verification_v47.py`
- `tooling/scripts/tgp_chain_Phi0_to_masses.py`
- `tooling/scripts/tgp_bridge_substrate_g0e.py`
- `tooling/scripts/tgp_unified_predictions_v47.py`
- `tooling/scripts/tgp_master_consistency_v47.py`
- `tooling/scripts/tgp_prediction_taxonomy_v47.py`

**Reason for defer:** scripts use `K_geo·m_sp² = π·Φ_0²` and `κ = 3/(4·Φ_0)`
self-consistency checks. Naïve replace 24.65 → 24.783 would break internal
cross-checks bez re-derivation `K_geo` and other coupled parameters.

**Co potrzebne:**

1. Phase 1: explicit identyfikacja wszystkich coupled parameters (`K_geo`,
   `m_sp`, `κ`, `c_src`) per script
2. Phase 2: re-derive coupled values for `Φ_0 = 24.783` (fixing one breaks others)
3. Phase 3: full re-run all 12+ self-consistency checks per script
4. Phase 4: PASS rate verification (currently 5/5 LOCK at 24.65; expected
   5/5 at 24.783 after coupled re-derivation)

**Kandydat dostawcy:** dedicated cycle `op-tooling-Phi0-unification` (estymata:
2-3 sesje, mostly mechanical re-derivation + re-run).

**Typ:** numerical re-verification

**Wpływ:** post-D01 cleanliness; obecnie tylko inconsistency w niezależnych
self-checks. Predykcje TGP nie zależą od tych internal Φ_0 (B3-v2 explicit
documented).

### N2: Brannen Φ_0 formal derivation

**Status:** OPEN — P3 (B3-v2 outstanding follow-up #2).

**Problem:** `Φ_0 = 24.78296...` jest currently *stated numerically* w
ROADMAP_v3:878-892 jako "Brannen variational principle", ale **formal
derivation** w core LaTeX (sek08a, sek09) nie istnieje.

**Co potrzebne:**

1. Identyfikacja variational step w core (`sek08a:???` lub `sek09:???`)
2. Explicit Lagrangian + variational principle dający `Φ_0` jako jedyne
   self-consistent rozwiązanie
3. Sympy LOCK: numerical match z 24.78296 (5+ sig figs)
4. Add as theorem `thm:Phi0-Brannen-derivation` w sek08a lub sek09
   (NON-BREAKING addytywne)

**Kandydat dostawcy:** dedicated cycle `op-Phi0-Brannen-formal-derivation`
(estymata: 2-3 sesje formal physics).

**Typ:** derivation (analytical)

**Wpływ:** krytyczne dla zewnętrznego review — *predykcja Φ_0* musi być
explicit derived, nie tylko *stated numerically*.

### N3: Cosmological Ω_Λ tension tracking

**Status:** OPEN — P4 (low priority, monitoring).

**Problem:** Brannen Φ_0 = 24.783 → Ω_Λ^TGP = Φ_0/36 = 0.6884.
Planck DR3 (2018): `Ω_Λ = 0.6847 ± 0.0073`. Diff: 0.0037 (+0.5σ).

Currently w tolerancji, ale **gdy Planck/DESI precision improves do σ_Ω = 0.005
(spodziewane DESI DR3 + CMB-S4 ~2030)**, tension wyniesie ~0.7σ.

**Co potrzebne:**

1. Track Planck/DESI/CMB-S4 latest precision (cykliczny check raz/rok)
2. If σ < 0.005 i wartość drifts away from 0.6884: critical re-evaluation
   anchor decision (Brannen vs Planck-derived)
3. Possible resolution path: TGP higher-order correction Φ_0 = 24.783 +
   δΦ_0 quantum (1-loop) — może zmniejszyć tension

**Kandydat dostawcy:** monitoring + dedicated cycle when σ_Ω < 0.005
(estymata: dependent on observational schedule, ~2030)

**Typ:** observational tracking + theoretical refinement

**Wpływ:** falsyfikatorowy test for TGP cosmology.

### N4: Predictivity ratio formal re-derivation

**Status:** OPEN — P3.

**Problem:** Pre-D01 LS-8 audit dał `ratio = 11/2 = 5.5` dla *wybranej*
trajektorii anchorów. Post-D01 my approximated ratio jako **~4.3** (bardziej
konserwatywne).

**Co potrzebne:**

1. Formal recount predykcji TGP z locked anchors:
   - Predykcje czysto algebraic (independent of anchor choice): list
   - Predykcje numeryczne z anchor dependency: list
2. Inputs count (3 fundamental: g_0^e, Ω_Λ, N=3)
3. Ratio = (# predictions) / (# free anchors) z explicit accounting
4. Update LS-8 audit z post-D01 number

**Kandydat dostawcy:** krótki audit (1 sesja, mostly bookkeeping).

**Typ:** meta-analytical

**Wpływ:** marketing claim "5.5" → "4.3" jest bardziej konserwatywne ale
**bardziej rzetelne** dla peer review.

### N5: m_H F11 vs CW reconciliation

**Status:** OPEN — P3 (open physics question).

**Problem:** TGP ma DWIE niezależne ścieżki dla m_H:

- **Path A: F11 algebraic** v×57/112 = 125.31 GeV (algebra + VEV)
- **Path B: CW 1-loop** = 125.1 GeV (substrate physics, top + gauge)

Spread: 0.2 GeV (~0.16% w m_H scale). Obie zgodne z PDG 2024 = 125.20±0.11
w 1σ. Ale: **czy te dwie ścieżki powinny się zgadzać do mniejszego niż 0.16%?**

**Co potrzebne:**

1. Analiza: czy F11 i CW są **strukturalnie te same** (np. F11 to leading
   order CW + 57/112 emerging from rep counting), czy **niezależnymi
   prediction paths**?
2. Jeśli te same: derivation 57/112 ze 1-loop CW
3. Jeśli różne: explicit dokumentacja i porównanie z eksperymentem

**Kandydat dostawcy:** krótki dedicated cycle `op-mH-F11-CW-reconciliation`
(estymata: 1-2 sesje).

**Typ:** structural analysis

**Wpływ:** spójność TGP predictions w sektorze elektrosłabym.

### N6: papers_external full unification with B3-v2

**Status:** OPEN — P3.

**Problem:** D01 Phase 2B dodał 4 footnoty z annotation B3-v2 canonical lock,
ale **wartości pozostają 0.1174** (alt-formula N_f=5). Dla pełnej spójności
papers_external z resztą projektu, należy:

**Opcja A:** zachować alt-formuła 0.1174 + annotation (current state, NON-BREAKING)
**Opcja B:** przepisać papers_external na canonical 0.1184 (BREAKING dla
publication, ale jednolicie z core/)

**Co potrzebne:**

1. Decyzja autora: czy publikacja używa alt-formula (więcej intuicyjna z N_f=5)
   czy canonical (Brannen Φ_0)?
2. Jeśli Opcja B: pełny rewrite papers_external sekcji α_s (10 wystąpień
   value + 4 derivation paragraphs)

**Kandydat dostawcy:** decyzja autora + 1 sesja propagacji (jeśli Opcja B).

**Typ:** decision + propagation

**Wpływ:** jednorodność TGP claim w publikacjach.

## Pytania otwarte

- **Q1:** Czy `Φ_0 ≈ N_f²` (24.783 ≈ 25) ma głębszą strukturalną interpretację,
  czy jest numerologią?
- **Q2:** Dla `m_H` — czy F11 v×57/112 może być derived z 1-loop CW + topological
  consideration (rep counting `7×8=56`, +1 dla 1-loop correction = 57)?
- **Q3:** Z 6 deferred scripts internal Φ_0 = 24.65 — czy któryś z nich
  daje predykcję która **istotnie zależy** od Brannen vs Planck-derived
  Φ_0? (testowalność)

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (tooling unification) | dedicated cycle z re-verification | self-consistency preservation |
| B2 | N2 (Brannen formal derivation) | dedicated cycle | core LaTeX edit needed |
| B3 | N3 (Ω_Λ tracking) | observational precision improvement | DESI DR3 / CMB-S4 ~2030 |
| B4 | N4 (predictivity ratio) | krótki audit | bookkeeping only |
| B5 | N5 (m_H reconciliation) | structural analysis | open physics |
| B6 | N6 (papers_external decision) | autor decision | binary choice |

## Closed needs (po D01 cyklu 2026-05-06)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| α_s spread | 0.1184 lock + 0.1174 alt-form annotation | B3-v2 + D01 Phase 2B | 2026-05-01 + 2026-05-06 |
| Φ_0 spread | 24.783 Brannen lock | B3-v2 14 lokacji | 2026-05-01 |
| m_H stale values | 124 → 125.1 (1-loop CW); 125.25 → 125.20 (PDG 2024) | D01 Phase 2A 5 edycji | 2026-05-06 |
| Σm_ν spread | 59.01 Z1 anchor | B4 + D01 Phase 2A 2 edycji | 2026-05-01 + 2026-05-06 |
| g_0^e formulacja | 0.86941 substratowa α=1 | LP-6 + L04 | 2026-05-04 |
| Status m_H z dwóch ścieżek | F11 (125.31) + CW (125.1) explicit dokumentowane | D01 Phase 2A | 2026-05-06 |

## Cross-references

- [[README.md]] — cykl D01 indeks
- [[Phase0_balance.md]] — CALIBRATION_PROTOCOL §2
- [[FINDINGS.md]] — eksportowalne wyniki
- [[../../audyt/D01_drifting_numbers/]] — audit-source
- [[../op-newton-momentum/B3_v2_alphas_propagation_results.md]] — pre-existing α_s lock
- [[../op-tooling-Phi0-unification]] (przyszły cykl, nie istniejący — N1)
- [[../op-Phi0-Brannen-formal-derivation]] (przyszły cykl — N2)
- [[../op-mH-F11-CW-reconciliation]] (przyszły cykl — N5)
