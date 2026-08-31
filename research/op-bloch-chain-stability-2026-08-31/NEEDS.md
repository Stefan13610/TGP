---
title: "NEEDS — op-bloch-chain-stability (USER-GATED; nic z poniższego nie zostało wykonane w rdzeniu)"
date: 2026-08-31
type: needs
tgp_owner: research/op-bloch-chain-stability-2026-08-31
status: OPEN (user-gate)
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/NEEDS.md]]"
---

# NEEDS (user-gated) — wg drzewa decyzyjnego LOCKa §6, gałąź Q-FAIL

Werdykt cyklu: **Q-FAIL** (odczyt PRIMARY wg [[Phase3_verdict_ruling.md]];
strict-literal: Q-INCONCLUSIVE — oba w [[Phase_FINAL_close.md]]).
Rdzeń `.tex` NIETKNIĘTY; wszystkie wnioski poniżej wymagają zgody usera.

## N1 — dopisek do Limitations (sek08b), dosłownie wg drzewa LOCKa

Proponowany wniosek do sekcji Limitations: **periodyczna kąpiel
samouzgodniona NIE stabilizuje sektora tachionowego w 1D (klasa zbadana:
łańcuch periodyczny EL Noty kanonicznej, tła samouzgodnione d∈{3π,4π,6π},
analiza pasmowa Blocha + test nieliniowy)** — ω²_min(d) = −1.222…−1.230 < 0
zbieżnie, ucieczka nieliniowa ≤ 2t*_ref; dla d ≤ 2π łańcuch w ogóle nie
istnieje (kolaps do próżni; okres orbit ≥ 2π). Pozostaje 3D
i warianty poza-klasowe.

## N2 — kandydat na cykl-następcę: wariant 3D (podniesiony priorytet)

Obserwacja analityczna z zamknięcia (post-hoc, deskryptywna): w 1D wynik
negatywny jest STRUKTURALNY — w formie kanonicznej (u=g³/3) problem
fluktuacji to równanie Hilla, mod translacyjny u₀′ ma węzły, więc stan
podstawowy leży ściśle poniżej zera dla KAŻDEGO niestałego tła
periodycznego 1D. Klasa 1D nie mogła dać Q-PASS inaczej niż numerycznym
artefaktem. **W 3D (sieć solitonów radialnych; czynnik 1/r w ogonach,
inna struktura modów) argument węzłowy w tej formie NIE obowiązuje** —
pytanie o stabilizację gęstością pozostaje tam autentycznie otwarte
i jest jedynym pozostałym wariantem tej linii (zgodnie z zakresem
zapisanym w LOCK §0 tego cyklu i NEEDS N2 poprzednika).

## N3 — nota metodologiczna (opcjonalna, do wiedzy warsztatowej)

Konstrukcja „samouzgodnione tło + Bloch" przeniesiona z Q2 poprzednika
na sektor dynamiczny zadziałała CZYSTO (pełna zbieżność 4/4 teł,
kontrole P3b 5/5 i P3c z kotwicą +1.88310 odtworzoną, cross-check dwóch
reprezentacji operatora do 8.8e−6, zgodność σ_fit z √|ω²_min| do 2.6%):
niezdolność rozstrzygania z Q1 poprzednika była własnością komórki
zero-flux z obciętym profilem, nie sektora. Warta odnotowania przy
ewentualnym projektowaniu wariantu 3D (tam analogiczna higiena:
tło samouzgodnione + BC periodyczne + identyfikacja Goldstone'ów
przed interpretacją).
