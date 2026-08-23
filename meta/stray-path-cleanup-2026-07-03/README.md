---
title: "Kwarantanna — porządek ścieżek podwojonych (2026-07-03, sesja #62)"
date: 2026-07-03
type: housekeeping
tags: [housekeeping, stray-paths, quarantine]
---

# Kwarantanna: sprzątanie podwojonych ścieżek (2026-07-03)

Artefakt wcześniejszych sesji: agenci zapisywali pliki do ścieżek
podwojonych typu `<cykl>/TGP/TGP_v1/<ścieżka-docelowa>` (względna
ścieżka doklejona do bieżącego katalogu zamiast korzenia repo).
Znaleziono 28 plików w 8 lokalizacjach. Dyspozycja:

- **15 plików PRZENIESIONO** do właściwych lokalizacji (istniały TYLKO
  w złej ścieżce — m.in. `Phase_FINAL_close.md` cykli: op-sigma-status-
  propagation-audit-2026-06-20, op-c0-derivation-from-substrate-2026-06-22
  (+README!), op-CE-H-3D-native-interaction-2026-05-22,
  op-Kgeo-from-D-uniqueness-2026-06-26, op-L08-Phase6-Dirac-propagator-
  2026-05-16, op-T34-normalization-amendment-2026-05-09 (+HANDOFF)).
- **9 duplikatów bajt-w-bajt USUNIĘTO**.
- **6 plików RÓŻNIĄCYCH SIĘ od wersji docelowych** → tutaj (starsze/
  częściowe snapshoty; wersje docelowe są nowsze i większe — nic nie
  nadpisano):
  - `research/op-A-derivation-from-CG-2026-06-25/Phase1_chain.txt`
    (317B vs 1409B w celu)
  - `papers/M911_LIGO3G_paper/paper_draft.md` (2026-05-07 vs 05-10)
  - `research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md`
    (05-07 vs 05-09)
  - `research/op-LIGO-3G-deviation/README.md` (05-07 vs 05-09)
  - `research/op-ppE-mapping/Phase1_results.md` (05-07 vs 05-09)
  - `main_build.log` (115B stub latexmk)
- Puste katalogi zagnieżdżone usunięte; `find -path "*/TGP/TGP_v1/*"`
  → 0 trafień.

Po weryfikacji (git diff pokaże przeniesienia) folder można usunąć.
