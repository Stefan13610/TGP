---
title: "INTAKE — op-xi-photon-ring → core (ξ.1 photon-ring resolution)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/op-xi-photon-ring
target_section: core/sek06_czarne_dziury + paper_bh_shadow + dodatekC_ringdown
hotspot_addressed: NEW (uzupełnia istniejące BH shadow paper)
agent_requesting: Sesja 7 (agent shortlist)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/op-xi-photon-ring/Phase3_results.md]]"
  - "[[research/op-xi-photon-ring/program.md]]"
  - "[[research/op-uv-as-ngfp]]"
  - "[[paper_bh_shadow]]"
tags:
  - intake
  - core-promotion
  - photon-ring
  - EHT
  - BH-shadow
  - xi-factor
---

# INTAKE — `op-xi-photon-ring` → core (ξ.1 photon-ring resolution)

## 1. Wynik proponowany do promocji

**ξ.1 program END (closed 2026-04-29):** ξ-factor structural decomposition
photon-ring scale dla Schwarzschild geometrii, plus UV-route map dla
N_A normalization. ξ-factor RG-invariant verified across UV scales.

**Cross-link z UV.1:** ξ.1 program zostawił "single residual N_A = 8.7719
algebraic provenance" — UV.1 (op-uv-as-ngfp) zamknął tę lukę przez
AS NGFP first-principles. Po UV.1.Phase3 status cascade: ξ.1 PARTIALLY
DERIVED → **DERIVED (refined²)**.

Sesja 4 audit: PASS=135, CLOSED=19, DERIVED=92, LOCKED=86 — silny cycle.

## 2. Source w research/

- `research/op-xi-photon-ring/Phase1_results.md` (frame audit)
- `research/op-xi-photon-ring/Phase2_results.md` (a₂ derivation)
- `research/op-xi-photon-ring/Phase3_results.md` (predictions + UV-route map)
- `research/op-xi-photon-ring/program.md`

## 3. Cross-validation z audytem 2026-05-01

Audyt **nie wymienia** ξ.1 explicite. Closure cascade UV.1 → ξ.1 jest poza
P-list audytu.

**Bridge:** ξ.1 jest źródłem dla `op-eps-photon-ring` (ε.1) — ε.1 jest **inną**
photon-ring decomposition (ε_ph structural decomp), które wraz z ξ.1 tworzy
2-fold consistency check.

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `paper_bh_shadow/tgp_bh_shadow.tex` | Update: dodać ξ-factor decomposition + ngEHT 2030+ falsifier | YES | dyskusja vs new paper |
| `core/sek06_czarne_dziury` | Update photon-sphere subsection z ξ-factor | YES (photon sphere r_ph < 3M) | dopisać 1-2 paragrafy |
| `core/dodatekC_ringdown` | Update: connect ξ-factor z QNM | YES | manualne sprawdzenie |
| `PREDICTIONS_REGISTRY.md` | ξ-cluster predictions | TBD | dopisać entries |

## 5. Hotspot adresowany

Brak istniejącego hotspota. NEW extension `paper_bh_shadow` papieru.

## 6. Co MUSI być sprawdzone

- [ ] ξ-factor RG-invariance verification z UV.1 cascade chain
- [ ] Cross-consistency z ε.1 (op-eps-photon-ring)
- [ ] Connection do EHT M87*/Sgr A* observed shadow radii
- [ ] Falsifier: ngEHT 2030+ photon-ring sharpening band

## 7. Co stanie się po akceptacji

1. **Człowiek** edytuje `paper_bh_shadow/tgp_bh_shadow.tex` dodając ξ-decomposition
2. **Człowiek** updateuje `core/sek06_czarne_dziury` photon-sphere subsection
3. **Człowiek** updateuje `core/dodatekC_ringdown` (jeśli relevant)
4. **Człowiek** dodaje ξ-cluster predictions do `PREDICTIONS_REGISTRY.md`
5. **Człowiek** YAML update `op-xi-photon-ring/README.md` na `level: L4`

## 8. Decyzja człowieka

- [ ] **ACCEPT** — ξ.1 promoted, BH shadow paper extended (data: ____)
- [ ] **DEFER** — czekamy na ngEHT 2030+ data point (powód: ____)
- [ ] **REJECT** — wynik utrzymany w research (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

## 9. Notatki agent

**PRO:** ξ.1 jest tightly cross-validated z UV.1 (cascade) i ε.1 (alternative photon-ring decomposition). Najsilniejszy CLOSED count (19) z 10 L3 candidates. Connection z paper_bh_shadow (już w core jako Paper III) jest naturalne.

**CONTRA:** ξ-factor jest concept relatively niski-poziomu (photon-ring resolution); promocja do core wymaga zdecydowanej koordynacji z istniejącą BH shadow narracją. Może lepiej extension paperu niż core sek06.

**Rekomendacja:** **ACCEPT z preferowaną ścieżką paper extension** — ξ-decomposition trafia do `paper_bh_shadow/tgp_bh_shadow.tex` jako new section (nie core/sek06 directly). Mniej intrusive dla core, zachowuje paper-as-publication semantics.
