---
title: "Phase 3 — F-ND-D (Q-D2): przegląd asymetrii / sortowanie ND (INFORMATIONAL)"
type: phase3_derivation
status: COMPLETE
date: 2026-06-13
cycle: op-nucleation-dimensionality-2026-06-13
authorization: "User 2026-06-13: 'działaj' → Phase 3 (Q-D2; Phase2_derivation §8 menu poz. 1)"
verdict:
  F-ND-D: "SORT-MONOTONE (genuine oś sortowania malejąca; 'pik w D=3' okna życia = repakowanie F-ND-A Θ(D≥3), nie niezależny selektor)"
sympy: "[[Phase3_sympy.py]] / [[Phase3_sympy.txt]] — 5/5 PASS, 0 hardcoded, INFORMATIONAL"
binding: "H-SORT = mechanizm roboczy (forbidden #12; SCOPING §2) — werdykt NIE podnosi claim_status (ceiling C), NIE zmienia Q-D1"
anti_lakatos: COMPLIANT
---

# Phase 3 — F-ND-D (Q-D2: sortowanie / asymetria ND)

> **Autoryzacja:** „działaj" → Phase 3 (druga połowa pierwotnego pytania user: „wszystkie inne
> możliwości ND asymmetry"). **Wynik:** [[Phase3_sympy.py]] + [[Phase3_sympy.txt]] — **5/5 PASS**.
> **INFORMATIONAL:** H-SORT = mechanizm roboczy (forbidden #12) — nie podnosi claim_status (C).

## §1 — Co i jak (pre-derywacja Phase 0 §4.4)

Ściana = kowymiar 1 (każde D); wymuszony porządek pary jest 1D wzdłuż NORMALNEJ; w D≥2
partnerzy mogą omijać ścianę kanałem bocznym (D−1 kierunków równoległych). Pytanie: czy
wydajność sortowania ma pik w jakimś D — i czy to NIEZALEŻNY selektor.

## §2 — Wyniki

- **FP-D1 (wydajność wzdłuż normalnej):** E_sort(D)=⟨|cosθ|⟩_{S^{D−1}}=Γ(D/2)/(√π·Γ((D+1)/2)):

  | D | 1 | 2 | 3 | 4 | 5 | 6 |
  |---|---|---|---|---|---|---|
  | E_sort | 1.000 | 0.637 | 0.500 | 0.424 | 0.375 | 0.340 |

  **MONOTONICZNIE malejąca** — brak wewnętrznego piku w D=3.
- **FP-D2 (kanały boczne):** = D−1 (rosnące) — spójne: więcej dróg ominięcia ściany z D.
- **FP-D3 (okno życia — test pokusy „iloczynu wskaźników", R-ND-9):**
  W(D)=E_sort(D)·Θ(D≥3): {1:0, 2:0, 3:0.50, 4:0.42, 5:0.38, 6:0.34} ⇒ **pik w D=3**. **ALE**
  czynnik Θ(D≥3) to warunek konieczny TOPOLOGII (F-ND-A: punkty ⟺ π₂≠0 ⟺ D≥3), **nie
  niezależny czynnik** ⇒ pik = **repakowanie A+B** (D≥3 z topologii × spadek geometryczny),
  NIE nowy derived selektor.
- **FP-D4 (guard H-SORT):** mechanizm roboczy — werdykt nie podnosi claim_status, nie zmienia
  Q-D1, nie może być cytowany jako ustalona bariogeneza.
- **FP-D5:** circularity guard (pik w 3 raportowany z jawną atrybucją do topologii; D_obs nieużyte).

## §3 — Werdykt

**F-ND-D = SORT-MONOTONE.** Genuine oś sortowania (wydajność wzdłuż normalnej) **maleje
gładko** z D — nie wyróżnia D=3 jako MECHANIZM. „Pik w D=3" okna życia istnieje, lecz wyłącznie
przez topologiczny czynnik Θ(D≥3) z F-ND-A (nie niezależny) — uczciwie odróżnione od derived
selekcji (odporność na reverse-engineering z R-ND-9). Odpowiedź na Q-D2: **asymetria/sortowanie
ND nie dostarcza dodatkowego, niezależnego wyróżnienia D=3** poza tym, co już daje topologia.

## §4 — Anti-Lakatos compliance (Phase 3)

INFORMATIONAL zadeklarowane ✓ · H-SORT working-mechanism uszanowany (ceiling C nietknięty) ✓ ·
pokusa „iloczynu wskaźników" zidentyfikowana i rozbrojona jawną atrybucją (pik = repakowanie,
nie selektor) ✓ · D_obs=3 nieużyte jako input ✓ · 0 nowych pól/stałych ✓.

## §5 — Decision menu (po Phase 3)

Q-D1 i Q-D2 rozstrzygnięte (4/4 osie). Następny krok:

1. **„działaj" → Phase FINAL (agregat F-ND-E)** — domknięcie cyklu: werdykt
   **DIM-3-PREFERENTIAL + SEK07A-CHALLENGED**, DOUBTS register, dyspozycja statusu sek07a (user).
2. Dodatkowy audyt przed FINAL (np. wyprowadzenie A,B,C(d) z {β,γ,Φ₀,λ}).
