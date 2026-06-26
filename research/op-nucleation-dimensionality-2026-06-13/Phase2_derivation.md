---
title: "Phase 2 — F-ND-C: nukleacja S_E(D) + marginalność grawitacyjna w D (value-blind test ostrego selektora)"
type: phase2_derivation
status: COMPLETE
date: 2026-06-13
cycle: op-nucleation-dimensionality-2026-06-13
authorization: "User 2026-06-13: 'działaj' → Phase 2 (Phase1_derivation §8 menu poz. 1)"
verdict:
  F-ND-C: "NUCL-MARG-NO-SELECTION"
sympy: "[[Phase2_sympy.py]] / [[Phase2_sympy.txt]] — 8/8 PASS, 0 hardcoded, werdykt z flag"
reuse_locked: "op-frontier-creation-rate-derivation (FCR): marginalność γ-3 (1/2)v_c²=GM/(ct); ρ̄=3H²/(8πG) EXACT; indicial p=2/3 EdS"
anti_lakatos: COMPLIANT
---

# Phase 2 — F-ND-C: nukleacja + marginalność grawitacyjna

> **Autoryzacja:** „działaj" → Phase 2. **Wynik:** [[Phase2_sympy.py]] + [[Phase2_sympy.txt]]
> — **8/8 PASS**, 0 hardcoded, circularity guard czysty (FP-C8), werdykt z flag.
> **Pytanie:** czy nukleacja w D wymiarach LUB marginalność grawitacyjno-wzrostowa dostarcza
> OSTREGO (derived) noża wycinającego pojedyncze D=3 — czego nie dały F-ND-A/F-ND-B?

## §1 — Nukleacja thin-wall S_E(D)

Bounce w d wymiarach euklidesowych: S_E = σ·Ω_{d−1}·R^{d−1} − ε·(Ω_{d−1}/d)·R^d
(powierzchnia − objętość; Ω_{d−1}=2π^{d/2}/Γ(d/2)).

- **FP-C1:** R_c(d) = (d−1)σ/ε; **B(d) = S_E(R_c) = [Ω_{d−1}/d]·(d−1)^{d−1}·σ^d/ε^{d−1}**
  (wyprowadzone symbolicznie).
- **FP-C2:** czynnik geometryczny g(d)=Ω_{d−1}/d·(d−1)^{d−1} dla σ=ε=1:

  | d | 1 | 2 | 3 | 4 | 5 | 6 |
  |---|---|---|---|---|---|---|
  | g(d) | 2.00 | 3.14 | 16.76 | 133.2 | 1347 | 16149 |

  **B(d) rośnie MONOTONICZNIE** — brak ekstremum wewnętrznego, **brak piku w d=3**. Tempo
  Γ∝exp(−B) jest NAJWYŻSZE dla NISKIEGO d ⇒ nukleacja, wzięta dosłownie jako mechanizm
  selekcji, faworyzowałaby d=1,2, NIE d=3.
- **FP-C3:** B(d) ∝ σ^d/ε^{d−1} (σ,ε wymiarowe) ⇒ porównanie tempa między różnymi D wymaga
  **wstrzykniętej skali** ⇒ brak value-blind selekcji D z nukleacji.

**Część nukleacyjna: NO-SELECTION** (a dosłownie: anty-selekcja D=3, ku niskim d).

## §2 — Marginalność grawitacyjna (uogólnienie FCR γ-3 na D)

Reuse LOCKED (FCR `op-frontier-creation-rate-derivation`): marginalność (1/2)v_c²=GM/(ct);
ρ̄=3H²/(8πG) EXACT; indicial p=2/3 (EdS). Uogólnienie na D = parametr:

- **FP-C4 (zasada):** trychotomia stabilności stacjonarnej kreacji (dE<0 unbounded / dE>0 brak
  kreacji / dE=0 marginal FORCED) **nie odwołuje się do D** ⇒ marginalność jako PRYNCYPIUM
  jest **D-niezależna** ⇒ nie selekcjonuje wymiaru.
- **FP-C5 (domknięcie):** z (1/2)v_c²=G_D M/R^{D−2}, R=ct, M=ρ̄·V_D, H=1/t:
  **ρ̄(D) = [D·v_c²/(2·Ω_{D−1}·c²)]·H²/G_D** — domyka się dla **symbolicznego D** (zero
  specjalnego D). W D=3 (v_c=c, Ω₂=4π): ρ̄=3H²/(8πG) — **odtworzone EXACT zgodnie z FCR LOCKED**.
- **FP-C6 (wzrost):** D-wymiarowy EdS a∝t^{2/D} ⇒ przepływ v_flow=(2/D)c, **indicial p(D)=2/D**
  — gładkie w D, brak wyróżnionego całkowitego D (w D=3: p=2/3, v=2c/3 — zgodne z FCR B-k4).

**Część marginalnościowa: NO-SELECTION** (zasada D-niezależna; księgowość domyka się dla
każdego D z odpowiednim G_D).

## §3 — Audyt niezależności (FP-C7) — dlaczego Bertrand/Ehrenfest NIE liczy się tu

Klasyczny fakt „stabilne orbity związane istnieją tylko w D=3" (Bertrand 1873; Ehrenfest 1917)
**autentycznie** wyróżnia D=3 — ale **nie** jako niezależny derived selektor F-ND-C:
(a) to twierdzenie **mechaniki klasycznej** (import, nie wyprowadzone z aksjomatów TGP);
(b) dotyczy stabilności stanów związanych / studni potencjału — **pokrywa się z osią F-ND-B**
(byłby to double-count). ⇒ traktowane jako **comparison-only**, nie wkład do werdyktu.

## §4 — FP ledger

| FP | Wynik | Treść |
|---|---|---|
| FP-C1 | PASS | R_c(d), B(d) thin-wall symbolicznie |
| FP-C2 | PASS | g(d) monotonicznie rosnąca ⇒ brak piku d=3; Γ faworyzuje niskie d |
| FP-C3 | PASS | B(d)∝σ^d/ε^{d−1} ⇒ selekcja wymaga wstrzykniętej skali |
| FP-C4 | PASS | zasada marginalności D-niezależna |
| FP-C5 | PASS | ρ̄(D) domyka dla symbolicznego D; D=3 ⇒ 3H²/(8πG) EXACT (zgodne FCR) |
| FP-C6 | PASS | indicial p(D)=2/D gładkie; brak wyróżnionego D |
| FP-C7 | PASS | Bertrand/Ehrenfest = comparison-only (import + overlap F-ND-B) |
| FP-C8 | PASS | circularity guard (D=3 tylko comparison-only) |

**8/8 PASS · 0 hardcoded · werdykt z flag · D_obs=3 tylko comparison-only.**

## §5 — Werdykt Phase 2

**F-ND-C = NUCL-MARG-NO-SELECTION.** Ani nukleacja (monotoniczna B(d), faworyzuje niskie d,
wymaga skali), ani marginalność grawitacyjna (zasada D-niezależna, ρ̄ domyka dla każdego D,
indicial gładki) **nie dostarcza ostrego derived noża** wycinającego pojedyncze D=3.
Potwierdza to **anticipated outcome Phase 0 §8** (NUCL-MARG-NO-SELECTION jako najbardziej
prawdopodobny) i obraz „księgowość działa dla każdego D z odpowiednią stałą".

## §6 — Konsolidacja trzech osi (wejście do F-ND-E; agregat finalny w Phase FINAL)

| Oś | Falsyfikator | Werdykt Phase 1/2 | Co wycina |
|---|---|---|---|
| Topologia | **F-ND-A** | TOPO-NO-SELECTION (+GAP) | dolne ograniczenie **D≥3** (punkty ⟺ π₂≠0); brak unikalności (π₃≠0); korekta π₂(SO(3)/Z₂)=0 |
| Stabilność | **F-ND-B** | STAB-SELECTS-3-FITTED | **pasmo d≥3**; d≥4 tylko miękko (Θ(ν⁻¹)); A,B,C nie-derived |
| Nukleacja+marg. | **F-ND-C** | NUCL-MARG-NO-SELECTION | nic ostrego; nukleacja faworyzuje niskie d; marginalność D-parametryczna |

**Spójny obraz:** żadna oś nie dostarcza DERIVED, ostrego selektora pojedynczego D=3.
Wspólny, odporny rdzeń = **D≥3 konieczne** (topologia) + **D=3 najmocniejszy kandydat
preferencyjny** (stabilność: pasmo z miękkim odcięciem d≥4; zgodność marginalności i orbit).
Mocna teza sek07a „d=3 **jedynym realistycznym wyborem**" nie znajduje derived poparcia w
żadnej z trzech osi LIVE; teza słabsza, którą sek07a deklaruje obok („**preferencyjny**"),
jest wspierana.

## §7 — Anti-Lakatos compliance (Phase 2)

Klasy CLOSED literalne ✓ · werdykt z flag (0 hardcoded) ✓ · D_obs=3 tylko comparison-only z
guardem ✓ · reuse FCR jako PRZYPADEK D=3 podstawienia, nie input (uogólnienie value-blind) ✓ ·
Bertrand/Ehrenfest uczciwie wykluczone jako import + double-count (wbrew pokusie „efektownego
D=3") ✓ · 0 nowych pól/stałych ✓ · 0 edycji rdzenia ✓ · wynik = anticipated outcome §8
(negatyw nie łagodzony, pozytyw nie naciągany) ✓.

## §8 — Decision menu (po Phase 2)

Q-D1 (selekcja wymiaru) jest teraz **rozstrzygnięte na poziomie trzech osi LIVE**: brak
derived selekcji; D=3 preferencyjny na D≥3. Pozostają opcje:

1. **„działaj" → Phase 3 (F-ND-D, Q-D2 sortowanie ND)** — INFORMATIONAL; sprawdza, czy
   wydajność sortowania klasy H-SORT ma pik w jakimś D. **Uwaga:** H-SORT = mechanizm roboczy
   (forbidden #12); wynik NIE podniesie claim_status (ceiling C) ani nie zmieni Q-D1.
2. **„działaj" → Phase FINAL (agregat F-ND-E)** — domknięcie z werdyktem
   **DIM-3-PREFERENTIAL + SEK07A-CHALLENGED** (na unikalności + korekta π₂), DOUBTS register,
   i Twoja dyspozycja statusu sek07a (np. „jedyny realistyczny wybór" → „najmocniejszy kandydat
   preferencyjny"; korekta zapisu π₂(SO(3)/Z₂)). Q-D2 można pominąć (INFORMATIONAL).
3. **Dodatkowy audyt** (np. formalne wyprowadzenie A,B,C(d) z {β,γ,Φ₀,λ}) przed FINAL.

**Następny krok wymaga user „działaj" (lub wyboru).** Rdzeń sek07a nietknięty.
