---
title: "HANDOFF_PROMPT — gotowy prompt dla nowej sesji (agent-implementator Phase 1–3)"
date: 2026-08-23
type: handoff
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: READY
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
---

# Prompt do wklejenia w NOWEJ sesji (skopiuj wszystko poniżej linii)

---

Jesteś agentem-implementatorem cyklu badawczego TGP. Realizujesz Phase 1–3
cyklu `TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/` DOKŁADNIE według
zalockowanego `Phase0_balance.md`. Pracujesz w vaultcie Obsidian; wszystkie
ścieżki względne od korzenia vaultu. Skrypty zapisuj w katalogu cyklu,
outputy do `Phase*_output.txt`.

## DWA PYTANIA BINARNE (szczegóły i kryteria: LOCK §0, §3)
1. **Q1:** czy mod runaway solitonu dostaje ω²>0 przy skończonej gęstości
   sąsiadów n — z ogonami/fazami ZMIERZONYMI (nie z dokumentacji rdzenia)?
2. **Q2:** czy z akcji gałęzi STABILNEJ (eq:field-eq-reproduced) znak
   efektywnego potencjału fluktuacji wokół tła o gęstości n staje się
   tachionowy? (test hipotezy: „różnica znaku = pojedynczy obiekt vs
   obiekt w kąpieli sąsiadów")

## KOLEJNOŚĆ CZYTANIA (obowiązkowa, przed jakimkolwiek kodem)
1. `TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/README.md`
2. `TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/Phase0_balance.md` — CAŁY LOCK; kryteria NIEZMIENNE
3. `TGP/TGP_v1/research/op-lattice-bath-runaway-2026-08-23/PhaseA_report.md` (werdykt bramki: co przeżyło audyt, co padło)
4. `TGP/TGP_v1/research/op-lattice-bath-runaway-2026-08-23/ANALIZA_N1_pochodzenie-faz_2026-08-23.md` + `N1_provenance_check.py` (układy M-P vs M-L, konwencje fitu — kod do reuse w Phase 1)
5. `TGP/TGP_v1/research/op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md` (dwie gałęzie znaku W — fundament Q2)
6. `TGP/TGP_v1/research/op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md` + `RETRO_oscillating_tail_lock.py` (model referencyjny drabiny d*)
7. `TGP/TGP_v1/research/op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py` + `Phase3_output.txt` (#63 V3 — kod do reuse: układ hamiltonowski, gate |ΔE|/E)
8. `TGP/TGP_v1/core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` lin. ~466–520 (prop:field-eq-from-action + Nota kanoniczna — źródło M-S dla Q2)

## REALIZACJA (wg LOCKa §3)
- **Phase 1:** `Phase1_tail_measured.py` — (A, δ) per gatunek {e, μ(, τ)}
  per model {M-P, M-L}; dwa zalockowane okna fitu [50,150]/[120,260]
  + R-kontrola; drabina d* par {ee, eμ, μμ}; kontrola negatywna P1c
  (Yukawa bez cos → 0 minimów, inaczej STOP i debug).
- **Phase 2 (Q1):** `Phase2_bath_runaway.py` — P2a baseline #63 V3
  (gate |ΔE|/E ≤ 1e−6, zbieżność dt; FAIL → STOP Q1), P2b skan 5 gęstości
  z drabiny (spektra ≥2 siatki × ≥2 komórki + ewolucja ±), P2c kontrola
  artefaktu komórki. Werdykt ŚCIŚLE wg Q1-PASS/Q1-FAIL/Q1-INCONCLUSIVE.
- **Phase 3 (Q2, NIEZALEŻNA od Q1):** `Phase3_two_sectors.py` — sektor
  STABILNY M-S: samouzgodnione tła ψ_n dla d ∈ {∞, 8, 6, 4, 2}
  (kontrola d=∞ MUSI odtworzyć Yukawę m²=+γ); spektrum fluktuacji ω²_min(d)
  + statyczna odpowiedź (monotoniczna vs oscylacyjna). Werdykt wg
  Q2-PASS/Q2-FAIL/Q2-INCONCLUSIVE.
- Na koniec: `Phase_FINAL_close.md` (frontmatter w stylu
  `research/op-wall-dynamics-2026-07-03/Phase_FINAL_close.md`; werdykty Q1
  i Q2 wg drzewa LOCKa §6), `NEEDS.md` (user-gated), aktualizacja README.

## FORBIDDEN MOVES (LOCK §4 — bezwzględne)
- Zero zmian kryteriów/progów/punktów/okien po starcie; korekta tylko dla
  udokumentowanego błędu implementacji, opisana PRZED użyciem wyniku.
- Każdy test musi mieć osiągalny FAIL; kontrole P1c/P2c/P3a(d=∞) nieusuwalne.
- (A, δ) w Phase 2 wyłącznie z Phase 1; w Q2 nic z gałęzi tachionowej.
- Wyniki negatywne wprost z liczbami; niezbieżność jako niezbieżność;
  zakaz skanowania do celu.
- NIE dotykać core/*.tex, axioms, papers. NIE edytować STATE.md.
  NIE commitować do gita.
- Q_K=3/2, η_K=181/15, kalibracje g₀ — flagować jako INPUT.
- **Higiena ścieżek:** NIE używać `cd`/`Set-Location`; po każdym zapisie
  pliku zweryfikować materializację (Test-Path); jeśli plik wylądował
  w `<coś>/TGP/TGP_v1/…` — przenieść natychmiast i odnotować.
- Długie przebiegi: uruchamiaj z przechwyceniem do pliku (UTF-8,
  errors='replace'); jeśli skrypt nie kończy się w ~15 min, przerwij
  i odnotuj jako ustalenie (nie czekaj w nieskończoność).

## RAPORT KOŃCOWY (twoja ostatnia wiadomość, po polsku)
(1) gate Phase 1 + zmierzone (A, δ) i d* per para/model; (2) werdykt Q1
z liczbami (ω²_min per n, zbieżność, baseline); (3) werdykt Q2 z liczbami
(ω²_min(d), kontrola d=∞); (4) lista utworzonych plików; (5) odstępstwa/
problemy. Bez upiększania — wyniki negatywne są pełnoprawne.
