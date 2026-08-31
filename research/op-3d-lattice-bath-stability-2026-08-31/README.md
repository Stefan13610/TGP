---
title: "op-3d-lattice-bath-stability — RACHUNEK CENTRALNY w 3D: czy sieć solitonów o skończonej gęstości stabilizuje mod runaway"
date: 2026-08-31
type: research-cycle
status: CLOSED-GATE-FAIL-STOP-P1c
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
authorization: "User 2026-08-31: 'ok 1 Cykl 3D dla stabilności w kąpieli'"
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-symplectic-Jspectrum-2026-08-31/Phase_FINAL_close.md]]"
---

# op-3d-lattice-bath-stability (2026-08-31)

**Status: CLOSED — GATE-FAIL-STOP na P1c (most radialny→kartezjański).
Pytanie Q (ω²_min(3D)) POZOSTAJE OTWARTE — to STOP maszynerii, nie wynik
fizyczny.** Szczegóły: [[Phase_FINAL_close.md]]; opcje kontynuacji
(user-gated): [[NEEDS.md]].

Ostatnia niewykluczona geometria hipotezy „izolacja umiera, kolektyw
żyje": po wyczerpaniu 1D (trzy klasy dynamiki, wszystkie Q-FAIL —
bloch-chain, gradient flow, symplectic-Jspectrum) pytanie ω²(n)
przenosi się do 3D, gdzie negatyw 1D jest jawnie nieprzenośny
(argument węzłowy Hilla nie obowiązuje; ogony zanikają 1/r).

**Pytanie binarne:** czy istnieje separacja d sieci prostej kubicznej
solitonów μ (baseline #63: g₀=2.02117, λ_min=−1.3896, t*=3.62),
przy której ω²_min(d) > 0 na tle samouzgodnionym (Bloch po strefie
Brillouina, po odjęciu modów translacyjnych)?

Konstrukcja fail-fast: tania bramka P1 (dyspersja próżni 3D + kotwica
radialna + most radialny→kartezjański) PRZED ciężką relaksacją 3D;
istnienie sieci = osobna bramka Phase 2 (negatyw = pełnoprawne
CLOSED-GATE-STOP). Kryteria: [[Phase0_balance.md]] (LOCK zacommitowany
przed jakimkolwiek obliczeniem).

## Log faz

- 2026-08-31: Phase 0 LOCK zapisany. Obliczenia: NIEROZPOCZĘTE.
- 2026-08-31: `Phase_method_decisions.md` FROZEN (rozdział modeli:
  bramki P1b/P1c w dosłownym #63 M0-f_ε [tam żyją kotwice −1.3896/3.62],
  rachunek P1a/P2–P4 w akcji kanonicznej; newton_krylov; Bloch ±1;
  eigsh 'SA'; interpolacja splajnowa radialna→3D).
- 2026-08-31: **Phase 1** (`Phase1_gate3d.py` + `Phase1_output.txt`):
  P1a PASS (dyspersja próżni exact, maxerr 6.0e−4/3.4e−4, ratio 3.999);
  P1b PASS (λ_min=−1.38962 vs −1.3896±1e−3; t*=3.62 wszystkie biegi);
  P1d PASS (|ΔE|/E ≤ 7.4e−15, dt-konsystencja 1.3e−13);
  **P1c FAIL** (λ_min(3D)=−8.810/−7.537 przy h=0.3947/0.30 vs kotwica
  −1.3896 ±5%; t*_izo(3D)=0.17 vs 3.62±15%).
- 2026-08-31: Addenda diagnostyczne (2 pliki + outputy): przyczyna FAIL
  strukturalna — wąska sferyczna kieszeń Q₆₃ (−15.5 przy g≈0.75, ściana
  f_ε) niedorozdzielczalna przy h≈0.4; hipoteza błędu konwencji wagi
  OBALONA (T1–T4: RQ modu kotwicznego −1.419 przy N=100 = 2.1% od
  kotwicy; mod −8.81 w 99% na powłoce). Brak legalnej korekty per LOCK
  §4 p.1 (nie ma błędu implementacji).
- 2026-08-31: Korekta 1 (eigsh tol=0; `Phase_correction_note_eigsh.md`)
  udokumentowana przed Phase 3; Phase 3 nie wystartowała.
- 2026-08-31: **STOP wg LOCKa §3/§6 (P1c FAIL ⟹ bez Phase 2–4).**
  Skrypty `Phase2_relax_lattice3d.py`/`Phase3_bloch3d.py` napisane
  w trakcie biegu P1, NIEURUCHOMIONE (zero wyników) — zachowane jako
  artefakty do ewentualnego re-locku. Zamknięcie:
  [[Phase_FINAL_close.md]]; [[NEEDS.md]] (user-gated).
