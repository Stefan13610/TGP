---
title: "INTAKE — op-eta-wolfenstein → core (η.1 Wolfenstein A, ρ̄, η̄ first-principles)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/op-eta-wolfenstein
target_section: core/sek10_N0_wyprowadzenie + paper_lepton_masses (CKM extension)
hotspot_addressed: NEW (CKM Wolfenstein cycle, poza 43-list)
agent_requesting: Sesja 7 (agent shortlist)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/op-eta-wolfenstein/Phase3_results.md]]"
  - "[[research/op-eta-wolfenstein/program.md]]"
  - "[[research/op-eta2-denom-derivation]]"
  - "[[research/op-mu-pmns-phase-hardening]]"
tags:
  - intake
  - core-promotion
  - Wolfenstein
  - CKM
  - PMNS
  - first-principles-cascade
---

# INTAKE — `op-eta-wolfenstein` → core (η.1 Wolfenstein A, ρ̄, η̄ first-principles)

## 1. Wynik proponowany do promocji

**η.1 program END (closed 2026-04-29):** Wolfenstein parametrization
(A, ρ̄, η̄) first-principles cross-sector cascade. 6 predictions H1-H6 LIVE.

**Cross-sector cascade:** η.1 łączy CKM (kwarki) z analogicznymi PMNS structures
(leptony) — η.2 (op-eta2-denom-derivation) jest extension cycle dla denominator
derivation + α-residual cascade.

Sesja 4 audit: PASS=121, CLOSED=7, DERIVED=85, LOCKED=51.

## 2. Source w research/

- `research/op-eta-wolfenstein/Phase1_results.md` (foundation)
- `research/op-eta-wolfenstein/Phase2_results.md` (derivation tests)
- `research/op-eta-wolfenstein/Phase3_results.md` (H1-H6 LIVE, program END)
- `research/op-eta-wolfenstein/program.md`

## 3. Cross-validation

Audyt 2026-05-01 **nie wymienia** η.1 (świeży 2026-04-29 cycle).

Cross-link z η.2 (op-eta2-denom-derivation) **silny** — η.2 jest follow-up
cycle z α-residual cascade extension. Promotion η.1 powinna być koordynowana
z η.2 jako pakiet "Wolfenstein-cluster".

Cross-link z `paper_lepton_masses` (lepton mass derivation z g_0^e + φ-ladder) —
η.1 dostarcza CKM analog.

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `core/sek10_N0_wyprowadzenie` | Nowa subsection: "Wolfenstein A, ρ̄, η̄ z first-principles" | NO | nowa subsection ~40-60 linii |
| `paper_lepton_masses/` | Update: dodać CKM extension subsection (jeśli paper extends to cross-sector) | YES (lepton-only) | TBD: extend or new paper |
| `PREDICTIONS_REGISTRY.md` | 6 entries H1-H6 + cross-link do η.2 (HH1-HH6) | NO | dodać entries |

## 5. Hotspot adresowany

Brak istniejącego hotspota w `CORE_HOTSPOTS.md`. NEW structural derivation.

## 6. Co MUSI być sprawdzone

- [ ] Wolfenstein values match PDG 2024 (A=0.826±0.012, ρ̄=0.159±0.010, η̄=0.348±0.010)
- [ ] η.2 (denom derivation) consistency check przed cluster promotion
- [ ] Cross-sector cascade z PMNS (μ.1, ν.1, sek10)
- [ ] Konsystencja z `paper_lepton_masses` φ-ladder framework

## 7. Co stanie się po akceptacji

1. **Człowiek** edytuje `core/sek10_N0_wyprowadzenie.tex` Wolfenstein subsection
2. **Człowiek** decyduje cluster promotion (η.1 + η.2) lub seriatim (η.1 first)
3. **Człowiek** updateuje `paper_lepton_masses/` lub tworzy nowy CKM paper
4. **Człowiek** dopisuje H1-H6 + ewentualnie HH1-HH6 do `PREDICTIONS_REGISTRY.md`
5. **Człowiek** dopisuje wpis do `CORE_INVENTORY.md` §1.15

## 8. Decyzja człowieka

- [ ] **ACCEPT** — η.1 promoted (data: ____)
- [ ] **ACCEPT cluster** — η.1 + η.2 jako pakiet Wolfenstein (data: ____)
- [ ] **DEFER** — czekamy na peer-review CKM paper (powód: ____)
- [ ] **REJECT** — utrzymane w research (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

## 9. Notatki agent

**PRO:** Wolfenstein z first-principles to klasyczny goal SM extension. CKM-PMNS cross-sector cascade jest unique features TGP. Strong DERIVED=85 + LOCKED=51 markers.

**CONTRA:** Bardzo świeży cycle (2026-04-29). η.2 follow-up cycle jeszcze nie weryfikowany niezależnie. Promotion bez cluster verification ryzykuje rozjazd η.1 ↔ η.2.

**Rekomendacja:** **DEFER + verify η.2** — najpierw verify η.2 (op-eta2-denom-derivation, też L3 candidate) PASS markers, potem cluster promotion η.1+η.2 jako pakiet. To minimalizuje cross-cycle inconsistency risk.
