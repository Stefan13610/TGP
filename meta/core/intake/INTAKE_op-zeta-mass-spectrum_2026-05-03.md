---
title: "INTAKE — op-zeta-mass-spectrum → core (ζ.1 neutrino mass spectrum + PMNS)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/op-zeta-mass-spectrum
target_section: core/sek07_predykcje + dodatekF_hierarchia_mas + PREDICTIONS_REGISTRY
hotspot_addressed: H-B4 (Σm_ν anchor — ten cycle dostarcza first-principles derivation)
agent_requesting: Sesja 7 (agent shortlist)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/op-zeta-mass-spectrum/Phase3_results.md]]"
  - "[[research/op-zeta-mass-spectrum/program.md]]"
  - "[[partial_proofs/hierarchia_mas]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
tags:
  - intake
  - core-promotion
  - neutrino-mass
  - PMNS
  - mass-spectrum
  - Z1
  - Sigma-m-nu
---

# INTAKE — `op-zeta-mass-spectrum` → core (ζ.1 neutrino mass spectrum + PMNS)

## 1. Wynik proponowany do promocji

**ζ.1 program END (closed 2026-04-29):** neutrino mass-spectrum + PMNS
first-principles derivation. 6 predictions Z1-Z6 LIVE.

**Cross-link z H-B4:** Audyt 2026-05-01 H-B4 (Σm_ν anchor 59.01 vs 59.6 meV)
zamknięty `RESOLVED-AUDIT-ANNOTATION` w § O.3 — ζ.1 dostarcza
**first-principles derivation** zastępując annotation lock.

Sesja 4 audit: PASS=140, CLOSED=11, DERIVED=11, LOCKED=11 — solidny cycle.

## 2. Source w research/

- `research/op-zeta-mass-spectrum/Phase1_results.md`
- `research/op-zeta-mass-spectrum/Phase2_results.md`
- `research/op-zeta-mass-spectrum/Phase3_results.md` (Z1-Z6 LIVE)
- `research/op-zeta-mass-spectrum/program.md`

## 3. Cross-validation z audytem 2026-05-01

H-B4 (Σm_ν: 59.01 vs 59.6 meV niejednoznaczność) zamknięty annotation-only
w `PREDICTIONS_REGISTRY.md:940`. ζ.1 jest **upgrade** annotation → derivation.

Cross-link z `partial_proofs/hierarchia_mas/dodatekF` (mass hierarchy z R3 ODE)
istnieje — ζ.1 może dostarczyć neutrino part dla full hierarchy picture.

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `core/sek07_predykcje` | Update neutrino mass section — first-principles ζ.1 derivation | YES (Σm_ν = 59.6 meV w predictions table) | replace lub upgrade table entry |
| `partial_proofs/hierarchia_mas/dodatekF` | Extension neutrino mass z R3 ODE | YES (charged lepton + quark) | nowa subsection |
| `PREDICTIONS_REGISTRY.md` | Update Z1 entry + dodać Z2-Z6 | YES (Z1) | upgrade Z1 z `LOCKED` na `DERIVED` z ζ.1 cite |

## 5. Hotspot adresowany

**H-B4** w [[meta/core/CORE_HOTSPOTS.md]] — status `RESOLVED-AUDIT-ANNOTATION`.
Po promotion ζ.1 status zmienia się na `RESOLVED-AUDIT-ANNOTATION + RESOLVED-FULL-PROMOTED`.

## 6. Co MUSI być sprawdzone

- [ ] Σm_ν derivation z ζ.1 zgadza się z 59.01 meV anchor (post-bisection K(m₁)=1/2)
- [ ] PMNS angles z ζ.1 vs NuFIT 5.3 measured values
- [ ] Konsystencja z H-B4 audit-aware annotacjami w PREDICTIONS_REGISTRY
- [ ] Cross-link z `op-nu-majorana-phase-mbb` (neutrino majorana phase, Z3-related)

## 7. Co stanie się po akceptacji

1. **Człowiek** edytuje `core/sek07_predykcje.tex` neutrino mass subsection
2. **Człowiek** dopisuje `partial_proofs/hierarchia_mas/dodatekF` extension
3. **Człowiek** upgradeuje Z1-Z6 entries w `PREDICTIONS_REGISTRY.md`
4. **Człowiek** updateuje H-B4 status w `CORE_HOTSPOTS.md`
5. **Człowiek** dopisuje wpis do `CORE_INVENTORY.md` §1.8

## 8. Decyzja człowieka

- [ ] **ACCEPT** — ζ.1 promoted, H-B4 upgraded (data: ____)
- [ ] **DEFER** — czekamy na JUNO mass ordering verdict (powód: ____)
- [ ] **REJECT** — utrzymane w research (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

## 9. Notatki agent

**PRO:** Bezpośrednio adresuje istniejący hotspot H-B4 (audyt acknowledged jako annotation-only). Promocja jest **upgrade z annotation do derivation** — zwiększa epistemic strength bez zmiany predykcji. PMNS cluster jest spójny z `op-eta-wolfenstein` (CKM analog).

**CONTRA:** Mass ordering jest falsyfikatorem TGP — jeśli JUNO confirms inverted ordering, ζ.1 zostaje obalone. Promotion przed JUNO data point (timeline 2026-2030) ma ryzyko obsoletnessu.

**Rekomendacja:** **ACCEPT z conditional footnote** — promocja OK, ale z explicit caveat "validity contingent on JUNO normal-ordering verdict 2026-2030". To pasuje do TGP falsifiability culture (każda predykcja ma falsifier).
