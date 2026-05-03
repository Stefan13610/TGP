---
title: "INTAKE — op-theta-quark-koide → core (θ.1 Koide K_quark + K-taxonomy completion)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/op-theta-quark-koide
target_section: core/sek07_predykcje + paper_lepton_masses + partial_proofs/koide_fp
hotspot_addressed: NEW (rozszerza istniejące Koide K=2/3 lepton w core)
agent_requesting: Sesja 7 (agent shortlist)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/op-theta-quark-koide/Phase3_results.md]]"
  - "[[research/op-theta-quark-koide/program.md]]"
  - "[[partial_proofs/koide_fp]]"
  - "[[paper_lepton_masses]]"
tags:
  - intake
  - core-promotion
  - Koide
  - quark-sector
  - K-taxonomy
---

# INTAKE — `op-theta-quark-koide` → core (θ.1 Koide K_quark + K-taxonomy)

## 1. Wynik proponowany do promocji

**θ.1 program END (closed 2026-04-29):** quark-sektor Koide first-principles
derivation (`K_quark` analogicznie do `K_lepton = 2/3`) + K-taxonomy completion
(klasyfikacja Koide-like invariants across particle sectors).

6 predictions Q1-Q6 LIVE. Sesja 4 audit: PASS=146, CLOSED=14, DERIVED=94, LOCKED=51,
mtime 2026-04-29 — silny cycle z bogatą substantą.

## 2. Source w research/

- `research/op-theta-quark-koide/Phase1_results.md` (foundation map)
- `research/op-theta-quark-koide/Phase2_results.md` (K_quark derivation tests)
- `research/op-theta-quark-koide/Phase3_results.md` (Q1-Q6 predictions + program END)
- `research/op-theta-quark-koide/program.md` — 3-fazowy plan + Koide-cascade motywacja

## 3. Cross-validation

Audyt 2026-05-01 **nie wymienia** θ.1 explicite (świeży 2026-04-29).

Cross-link z `partial_proofs/koide_fp/` (Koide K_lepton = 2/3 algebraic identity)
jest **silny** — θ.1 jest extension tej samej pracy do quark sector.

`paper_lepton_masses/tgp_lepton_masses.tex` cytuje `K = 2/3` jako exact —
θ.1 dostarcza paralelę dla quarków.

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `core/sek07_predykcje` | Update Koide section — dodać `K_quark` + K-taxonomy | YES (Koide K=2/3 lepton) | nowa subsection ~30 linii |
| `paper_lepton_masses/` | Dodać quark-sector subsection (lub nowy paper) | YES (lepton-only) | dyskusja vs nowy paper |
| `partial_proofs/koide_fp/` | Extension dla K-taxonomy quark side | YES | nowy plik dod-quark.tex |
| `PREDICTIONS_REGISTRY.md` | 6 entries Q1-Q6 | NO | dopisać |

## 5. Hotspot adresowany

Brak istniejącego hotspota. NEW extension Koide K-formula z lepton do quark.

## 6. Co MUSI być sprawdzone

- [ ] K_quark algebraic identity weryfikowalna sympy (folder ma `LOCKED=51` markerów — sprawdzić pliki)
- [ ] Konsystencja z `paper_lepton_masses` Brannen parametrization B = √2
- [ ] Q1-Q6 predictions vs PDG quark masses
- [ ] Cross-link z `partial_proofs/hierarchia_mas/dodatekF`

## 7. Co stanie się po akceptacji

1. **Człowiek** edytuje `core/sek07_predykcje.tex` Koide subsection
2. **Człowiek** decyduje: nowy paper "TGP Quark Masses" czy extension lepton paper
3. **Człowiek** dodaje `partial_proofs/koide_fp/dod-quark.tex` lub równoważne
4. **Człowiek** dopisuje Q1-Q6 do `PREDICTIONS_REGISTRY.md`
5. **Człowiek** dopisuje wpis do `CORE_INVENTORY.md` §1.8

## 8. Decyzja człowieka

- [ ] **ACCEPT** — θ.1 promoted, K-taxonomy complete in core (data: ____)
- [ ] **DEFER** — czekamy na peer-review paper (powód: ____)
- [ ] **REJECT** — utrzymane w research (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

## 9. Notatki agent

**PRO:** Extension istniejącego K=2/3 lepton (już w core) na quarks jest **niskoryzykowna** strukturalnie — używa tej samej Brannen parametrization. Bogaty PASS count (146) + LOCKED markers (51) sugerują dojrzałą pracę.

**CONTRA:** Nie sprawdziłem czy K_quark identyfikacja zachodzi z dokładnością comparable do K_lepton (2/3 EXACT). Wymaga manual review Phase3.

**Rekomendacja:** **ACCEPT z modyfikacją:** najpierw verify Q1-Q6 PDG match, potem promotion. Jeśli match dobre, najnaturalsza ścieżka to extension `paper_lepton_masses` o quark subsection (zachowuje single-paper architekturę).
