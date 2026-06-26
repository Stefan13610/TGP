---
title: "SCOPING — op-phi-radiative-dof-audit (PR-025 exhaustiveness closure)"
date: 2026-06-13
type: scoping
status: READY (wymaga własnego Phase 0 + autoryzacji „działaj")
parent_verdict: "[[../research/op-PSR-Pdot-energy-balance-2026-06-13/]] PR-025 TRIGGERED (pending ratification)"
origin: "User 2026-06-13: identyfikacja brakującego warunku — «dowód, że Φ jest tylko zmienną pomocniczą dla σ_ab, a nie dodatkowym radiacyjnym stopniem swobody»; falsifier P_total = P_GR vs P_GR + P_φ/6"
tags:
  - scoping
  - PR-025-followup
  - Dirac-constraint-analysis
  - fast-kill
---

# SCOPING — czy Φ może być zmienną pomocniczą (non-radiative) w LIVE TGP?

## §1 — Pytanie (sformułowanie usera, doprecyzowane)

Istniejące elementy:
- C_i = q·m_i ⇒ brak monopolu/dipola skalarnego (PR-025 T0, potwierdzone),
- σ_ab sprzęga się z T_ab^TT i odtwarza GR **w amplitudzie** (κ_E nieprzypięte — PR-025 T5/R1 #23).

Brakujący warunek: jawne usunięcie niezależnego Φ-kwadrupolu, tj. dowód, że
mod oddechowy Φ jest WIĘZEM (constrained/auxiliary), nie radiacyjnym DOF.
Rozstrzyga alternatywę: **P_total = κ_E·P_GR** vs **P_total = κ_E·P_GR + (1/6)P_GR**.

## §2 — Pre-derywacja analityczna (R1 #18 compliance — zapisana PRZED Phase 0)

Oczekiwany wynik: **HONEST_NEGATIVE** („Φ auxiliary" NIE wyprowadzalne w LIVE), z trzech kroków:

1. **Hiperboliczność LOCKED:** zlinearyzowany L (Phase 5 emergent-metric, LOCKED)
   ma π = K₁δΦ̇ odwracalne ⇒ algorytm Diraca: zero więzów pierwotnych ⇒
   δΦ = propagujący DOF. (Test sympy: jawna analiza więzów.)
2. **Statyka ⇔ radiacja:** ten sam wierzchołek qρΦ i ten sam retardowany
   propagator dają −Gm₁m₂/d ORAZ P_φ = (1/6)P_GR. Usunięcie radiacji przy
   zachowaniu 1/r wymaga niekauzalnego propagatora. (Test sympy: rozkład
   bieguna G_ret; część statyczna i falowa z jednej funkcji.)
3. **Brak orbity cechowania:** w GR potencjał newtonowski jest niefizyczny
   radiacyjnie dzięki więzom z dyfeomorfizmów (00, 0i = constraints; tylko TT
   propaguje). Rzeczywisty skalar Z₂ nie ma symetrii ciągłej działającej na Φ
   ⇒ mod oddechowy fizyczny. Wyjścia: (a) nowa symetria cechowania = nowy
   aksjomat (narusza fundament §1 FOUNDATIONS), (b) kinetyka eliptyczna =
   złamanie Lorentza (zawala α_i ≡ 0 strukturalne, sektor PPN L2).

Jeśli rachunek da co innego niż pre-derywacja — rachunek nadrzędny.

## §3 — Falsyfikatory (szkic do LOCK w Phase 0)

| ID | Test | Reguła |
|---|---|---|
| F-AUX-A | Analiza więzów Diraca zlinearyzowanego LIVE L (sektor Φ + σ) — policzyć radiacyjne DOF | DOF_Φ = 0 ⟹ „auxiliary" DERIVED (PR-025 wymaga rewizji); DOF_Φ = 1 ⟹ HONEST_NEGATIVE, PR-025 EXHAUSTIVE-OVER-LIVE |
| F-AUX-B | Rozkład G_ret: czy istnieje modyfikacja LOCKED kinetyki zachowująca 1/r i kasująca część falową przy zachowaniu kauzalności? | NIE ⟹ potwierdzenie kroku 2; TAK ⟹ nowy mechanizm do osobnego PR |
| F-AUX-C | Inwentarz symetrii LIVE działających na Φ (S05 U(1) faza, Z₂) — czy któraś generuje wiąz na mod oddechowy? | NIE ⟹ krok 3 potwierdzony |

## §4 — Format i koszt

Fast-kill wzór GST: 1 faza merytoryczna, 1 sesja, sympy, 0 nowych stałych,
0 dostępu do R_obs (czysto strukturalny — dane pulsarowe nietykane).
Wynik każdego znaku jest informatywny:
- HONEST_NEGATIVE ⟹ PR-025 upgrade: „both branches" → „exhaustive over LIVE";
  sektor radiacyjny domknięty negatywnie analogicznie do galaktycznego (PR-004+GST).
- DERIVED (Φ auxiliary) ⟹ sensacja strukturalna; PR-025 wymaga re-run z
  P_total = κ_E·P_GR i osobnego domknięcia κ_E.

## §5 — Forbidden moves (szkic)

1. Zakaz modyfikacji L poza LOCKED inventory (każda modyfikacja = nowy PR, nie rescue).
2. Zakaz użycia R_obs / danych pulsarowych (cykl strukturalny).
3. Zakaz reinterpretacji werdyktu PR-025 przed zamknięciem F-AUX-A/B/C.
4. Pre-derywacja §2 jest oczekiwaniem, nie progiem — progi definiuje Phase 0.
