---
title: "Phase 1 — FAST AUDYT (F-ND-A topologia + F-ND-B stabilność): value-blind audyt argumentu sek07a za selekcją D=3"
type: phase1_derivation
status: COMPLETE
date: 2026-06-13
cycle: op-nucleation-dimensionality-2026-06-13
authorization: "User 2026-06-13: 'działaj' → Phase 1 FAST AUDYT (Phase0 §5 gate)"
verdicts:
  F-ND-A: "TOPO-NO-SELECTION (+ element GAP: rozmaitość porządku nieustalona; π₂(SO(3)/Z₂)=0 vs π₂(RP²)=ℤ)"
  F-ND-B: "STAB-SELECTS-3-FITTED (D=3 wyróżniony, lecz Δ_d>0 dla pasma d≥3; A,B,C nie-derived; d≥4 wykluczone tylko miękkim Θ(ν⁻¹))"
sympy: "[[Phase1_sympy.py]] / [[Phase1_sympy.txt]] — 13/13 PASS, 0 hardcoded T_pass, werdykty z flag"
anti_lakatos: COMPLIANT
---

# Phase 1 — FAST AUDYT (F-ND-A + F-ND-B)

> **Autoryzacja:** user „działaj" → Phase 1 (Phase0_balance.md §5).
> **Wynik:** [[Phase1_sympy.py]] + [[Phase1_sympy.txt]] — **13/13 PASS**, 0 hardcoded T_pass,
> circularity guard czysty (FP-A6/FP-B6), werdykty WYLICZONE z flag.
> **Obiekt audytu:** `core/sek07a_wymiar_wzmocniony.tex` (`prop:wymiar-quantitative`) — read-only.

## §1 — Co testujemy (przypomnienie ramy)

sek07a wzmacnia argument za D=3 dwiema osiami: **Część I (homotopia)** — wskaźnik
`N_sekt(d)` + teza, że dokładnie w d=3 są „trzy niezależne sektory" (ściana/struna/punkt)
„odpowiadające trzem generacjom"; **Część II (potencjał)** — warunek trzech reżimów
`Δ_d>0` z `Q(d)`→Q(3)=3. Konkluzja sek07a: „d=3 **jedynym realistycznym wyborem**" — przy
jednoczesnej deklaracji „Argument pozostaje **preferencyjny**". Phase 1 audytuje **obie osie
value-blind**: czy selekcja jest DERIVED z mechanizmu, czy artefaktem konstrukcji pod D_obs=3.

## §2 — F-ND-A: oś topologiczna

### §2.1 — Fakty homotopii (standardowa matematyka; zweryfikowane relacjami nakryć — FP-A0)

Reguła klasyfikacji defektów (Kibble/Mermin): w D wymiarach przestrzennych defekt
**punktowy** (0-wymiarowy, kowymiar D) jest stabilny topologicznie ⟺ **π_{D−1}(M_ord) ≠ 0**;
`N_sekt(D) = #{j ∈ 0..D−1 : π_j(M_ord) ≠ 0}` (ta sama reguła dla każdego D — FP-A2).

Grupy homotopii rozmaitości kandydujących (Hatcher §4.1; nakrycie uniwersalne ⇒ π_n iso n≥2):

| M | π₀ | π₁ | π₂ | π₃ | π₄ | π₅ |
|---|----|----|----|----|----|----|
| Z₂ (chiralne, ścianki) | ℤ₂ | 0 | 0 | 0 | 0 | 0 |
| S¹ = U(1) | 0 | ℤ | 0 | 0 | 0 | 0 |
| **S²** | 0 | 0 | **ℤ** | ℤ | ℤ₂ | ℤ₂ |
| **RP² = S²/Z₂** | 0 | ℤ₂ | **ℤ** | ℤ | ℤ₂ | ℤ₂ |
| SO(3) | 0 | ℤ₂ | **0** | ℤ | ℤ₂ | ℤ₂ |
| SO(3)/Z₂ | 0 | (skończona ≠0) | **0** | ℤ | ℤ₂ | ℤ₂ |

(Weryfikacja nakryć FP-A0: π_n(RP²)=π_n(S²) i π_n(SO(3)/Z₂)=π_n(S³) dla n≥2; π₂ dowolnej
grupy Liego i jej ilorazu dyskretnego = 0.)

### §2.2 — Rozstrzygnięcie sporu π₂ (FP-A1) — pierwszy wynik audytu

Phase 0 zarejestrował spór: sek07a pisze `π₂(SO(3)/Z₂)=ℤ`; SCOPING §3a: `π₂(RP²)=ℤ`.
**Rachunek (FP-A1):** π₂ **dowolnej** grupy Liego znika (tw. klasyczne), więc
**π₂(SO(3)/Z₂)=0** — zapis sek07a `π₂(SO(3)/Z₂)=ℤ` jest **niepoprawny jak literalnie
zapisany**. Jedyną rozmaitością w inwentarzu dającą defekty punktowe (π₂=ℤ) jest
**RP²=S²/Z₂** (lub S²). To NIE jest werdykt o teorii — to korekta matematyczna identyfikująca,
że teza „cząstki = defekty punktowe w 3D" wymaga `M_ord ⊇` sektora RP²/S², a nie SO(3)/Z₂
(uzasadnianego w sek07a przez spin-½ z π₁(SO(3))=ℤ₂).

### §2.3 — Dwie gałęzie i dlaczego ŻADNA nie daje czystej selekcji D=3

**Gałąź (α) — rozmaitość jak zapisana, M = chiralZ₂ × SO(3)/Z₂** (uzasadniona spinem-½):
π₂=0 ⇒ **brak stabilnych cząstek punktowych w D=3** (FP-A3/A4). Teza topologiczna
„cząstki w 3D" **upada na własnej rozmaitości TGP**. (Spójne z ustaleniem cyklu GST Phase 1:
π₂(S¹)=0 — uzwojenie U(1) chroni tylko defekty liniowe; punkty wymagają osobnego sektora —
rezyduał W-GST-4 „sektor RP² poza LIVE".)

**Gałąź (β) — naprawa, M = chiralZ₂ × RP²** (π₂=ℤ ⇒ punkty w D=3 ✓): ALE π₃(RP²)=π₃(S²)=ℤ≠0,
więc **π_{D−1}(M)≠0 dla każdego D≥3** ⇒ stabilne defekty punktowe istnieją w D=3 **i** D=4,5,6
(FP-A5). `N_sekt(D)` rośnie monotonicznie (3,4,5,6 dla D=3,4,5,6 — FP-A2): „dokładnie 3
sektory" to **przejście przez 3 przy D=3**, nie pik. D=4 ma **więcej** struktury defektowej
(4 sektory), nie mniej — odwrotnie niż sugeruje narracja „d≥4: średnie pole, brak struktury".

**Wniosek F-ND-A (wyliczony z flag):** **TOPO-NO-SELECTION** z elementem **GAP** (która
rozmaitość jest genuine, nieustalone jednoznacznie z aksjomatów jak podane). Odporna część:
**cząstki punktowe wymagają π₂≠0 ⇒ D≥3** (dolne ograniczenie — autentyczne). Mocna teza
„topologia wycina **dokładnie** D=3" **nie przetrwała** uczciwego testu D>3.

## §3 — F-ND-B: oś stabilności (Δ_d + Derrick)

### §3.1 — Reprodukcja Δ_d (FP-B1, symbolicznie)

Z `V_eff^(d) = −A/r^{d−2} + B/r^{d−1} − C/r^d`, F=−dV/dr; zera siły (u=1/r):
`d·C·u² − (d−1)·B·u + (d−2)·A = 0`; dyskryminanta
**Δ_d = (d−1)²B² − 4d(d−2)AC** — zgodna z sek07a `eq:Delta-d` (FP-B1 PASS, ekwiwalencja zer
zweryfikowana sympy).

### §3.2 — Próg τ_d i pasmo (FP-B2/B3)

Δ_d>0 ⟺ B/√(AC) > **τ_d = 2√(d(d−2))/(d−1)**: τ₃=√3≈1.732, τ₄≈1.886, τ₅≈1.936, τ₆≈1.960
(rosnący, ograniczony →2). sek07a **asercję** (nie derywację w dodatku) B₃/√(A₃C₃)≈3.4 — audyt
DERIVED-vs-FITTED (FP-B3): **A,B,C(d) NIE są wyprowadzone z {β,γ,Φ₀,λ}** w sek07a. Jeśli
stosunek ρ=B/√(AC) jest ~d-niezależny i =3.4, to **Δ_d>0 zachodzi dla d∈{2,3,4,5,6}** — warunek
trzech reżimów **NIE wyklucza d≥4** (sek07a sam to przyznaje: „Δ₄>0 możliwe").

### §3.3 — Derrick / bg-stabilizacja (FP-B4) i nóż Θ(ν⁻¹) (FP-B5)

Skalowanie E[Φ_λ]=λ^{d−2}E_grad+λ^d E_pot: bez stabilizatora brak równowagi L*>0 dla
E_grad,E_pot>0 i d≥2 (Derrick). Z członem tła ∝L^p równowaga możliwa w **paśmie d** (zależnie
od p) — nie wycina pojedynczego d=3. Jedynym czynnikiem wykluczającym d=4 w `Q(d)` jest
**Θ(ν₄⁻¹)=0** („pole średnie dokładne, d_c^Ising=4, fizyka trywialna"). To **fakt RG** (ν=1/2,
η=0 dla d≥4), ale **argument jakościowy** („trywialna" ≠ „niemożliwa/niegenerowalna"), nie
ostry próg derived. ⇒ **d=4 wykluczone wyłącznie miękko** (FP-B5).

**Wniosek F-ND-B:** **STAB-SELECTS-3-FITTED** — D=3 wyróżniony, lecz (i) A,B,C nie-derived,
(ii) Δ_d>0 dla pasma d≥3, (iii) d≥4 wykluczone tylko miękkim Θ(ν⁻¹) ⇒ selekcja
stabilnościowa **PREFERENCYJNA**, nie DERIVED.

## §4 — FP ledger

| FP | Wynik | Treść |
|---|---|---|
| FP-A0 | PASS | relacje nakryć spójne; π₂(grupa Liego)=0 |
| FP-A1 | PASS | spór rozstrzygnięty: π₂(SO(3)/Z₂)=0 ≠ ℤ (zapis sek07a niepoprawny); punkty z RP²/S² |
| FP-A2 | PASS | N_sekt(D) symetrycznie D=1..6 (SO3/Z2: 1,2,2,3,4,5; RP2: 1,2,3,4,5,6) |
| FP-A3 | PASS | punkty stabilne ⟺ π_{D−1}≠0; stablica per-D |
| FP-A4 | PASS | M zapisane (SO3/Z2): brak cząstek punktowych w D=3 (teza upada na własnej rozmaitości) |
| FP-A5 | PASS | uczciwy test D>3: punkty (RP2) dla każdego D≥3, nie tylko 3; N_sekt rośnie |
| FP-A6 | PASS | circularity guard (D_obs=3 tylko w comparison-only) |
| FP-B1 | PASS | Δ_d=(d−1)²B²−4d(d−2)AC symbolicznie (reprodukcja sek07a) |
| FP-B2 | PASS | τ_d=2√(d(d−2))/(d−1) dla d=2..6 |
| FP-B3 | PASS | DERIVED-vs-FITTED: A,B,C nie-derived; ρ=3.4 ⇒ Δ_d>0 dla d∈{2..6} |
| FP-B4 | PASS | Derrick: bez tła brak L*>0 (d≥2); z tłem pasmo d, nie punkt |
| FP-B5 | PASS | Θ(ν⁻¹): d=4 wykluczone tylko miękko (pole średnie) ⇒ PREFERENCYJNE |
| FP-B6 | PASS | circularity guard |

**13/13 PASS · 0 hardcoded T_pass · werdykty z flag · D_obs=3 wyłącznie comparison-only.**

## §5 — Werdykty Phase 1 (klasy CLOSED Phase0 §3)

- **F-ND-A = TOPO-NO-SELECTION** (+ element GAP: rozmaitość porządku nieustalona;
  π₂(SO(3)/Z₂)=0 vs π₂(RP²)=ℤ). Odporna część: **cząstki punktowe ⇒ D≥3**.
- **F-ND-B = STAB-SELECTS-3-FITTED** (D=3 wyróżniony preferencyjnie; nie-derived; pasmo d≥3;
  d≥4 tylko miękko).

**Trajektoria do F-ND-E (INFORMATIONAL; agregat dopiero w Phase FINAL po F-ND-C):**
mocna teza sek07a „**jedyny realistyczny wybór**" **nie przetrwała** value-blind audytu
osi topologicznej i stabilnościowej; przetrwała teza słabsza „D≥3 konieczne dla cząstek
punktowych" + „D=3 jest najmocniejszym kandydatem **preferencyjnym**" — co jest dokładnie
tym, co sek07a deklaruje obok („preferencyjny"), a czego nie obejmuje jego mocniejsze
sformułowanie „jedyny". Kierunek agregatu: **DIM-3-PREFERENTIAL** (na D≥3) z elementem
**SEK07A-CHALLENGED** (na unikalności + korekta π₂). Rozstrzygnięcie końcowe wymaga
**F-ND-C (Phase 2: nukleacja S_E(D) + marginalność grawitacyjna)** — czy któraś z tych osi
dostarcza ostrego (derived) noża wycinającego pojedyncze D=3, czy też potwierdza obraz
„pasmo/preferencja".

## §6 — Comparison-only (po locku; D_obs dozwolone tylko tutaj)

- D_obs=3 leży w paśmie dopuszczonym przez stabilność (d≥3) i spełnia warunek konieczny
  topologii (π₂≠0) — **zgodność, ale nie unikalność**.
- „3 generacje ↔ 3 sektory": N_sekt(3)[RP²]=3 (zgodne), lecz N_sekt(4)=4 — D=4 ma **więcej**
  sektorów. Koincydencja liczbowa, nie selekcja. (Zakaz traktowania jako dowód — Phase0 L3 map.)

## §7 — Anti-Lakatos compliance (Phase 1)

Klasy werdyktów CLOSED użyte literalnie ✓ · werdykty wyliczone z flag (0 hardcoded) ✓ ·
D_obs=3 wyłącznie comparison-only z guardem FP ✓ · uczciwy test D>3 wykonany (D=4,5,6;
punkty istnieją tam też) ✓ · spór π₂ rozstrzygnięty rachunkiem, nie wyborem pod wygodę ✓ ·
korekta sek07a raportowana wbrew interesowi „efektownego wyniku" (teza selekcyjna osłabiona,
nie wzmocniona) ✓ · 0 nowych pól/stałych ✓ · 0 edycji rdzenia (sek07a read-only; rewizja
statusu = user w Phase FINAL) ✓ · rozmaitość porządku jako GAP (nie rozstrzygnięta pod π₂≠0) ✓.

## §8 — Decision menu (po Phase 1)

1. **„działaj"** → **Phase 2 (F-ND-C)**: nukleacja w D wymiarach (S_E(D), tempo Γ) +
   marginalność grawitacyjno-wzrostowa ((1/2)v_c²=G_D M/R^{D−2}; analog γ-3 R=ct; indicial p(D)).
   Pytanie: czy któraś oś dostarcza **ostrego** noża na D=3 (zmiana werdyktu na DERIVED),
   czy potwierdza „pasmo/preferencja" → agregat F-ND-E.
2. **Wcześniejsze domknięcie (F-ND-E częściowy)** — jeśli uznasz, że osie A+B już rozstrzygają
   pytanie cyklu (mocna teza sek07a obalona, słaba D≥3 potwierdzona): przejście do Phase FINAL
   z werdyktem **SEK07A-CHALLENGED / DIM-3-PREFERENTIAL** + dyspozycja statusu sek07a (Twoja decyzja).
3. **Korekta zakresu / dodatkowy audyt** (np. formalne wyprowadzenie A,B,C(d) z {β,γ,Φ₀,λ},
   by przetestować STAB-DERIVED rzetelnie) przed kontynuacją.

**Następny krok wymaga user „działaj" (lub wyboru z menu).** Rdzeń sek07a pozostaje nietknięty;
jakakolwiek jego rewizja (np. korekta π₂, osłabienie „jedyny"→„najmocniejszy preferencyjny")
to wyłącznie Twoja decyzja w Phase FINAL.
