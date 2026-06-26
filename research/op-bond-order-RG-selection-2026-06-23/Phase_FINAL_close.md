---
title: "Phase FINAL — op-bond-order-RG-selection CLOSE. Werdykt value-blind (sympy 5/5): RG-NOT-SELECTED. [g_s]=−s(d−2)−(2s+2)γ; s=5 (α=2) NAJBARDZIEJ irrelevant (d=4: −10 < −4 (s=2) < 0 (free)); RG faworyzuje niższe s; NGFP escape wymaga skrajnego γ≈−0.83. Jedyny selektor α=2 = density-frame klasa konforemna C1–C3 (makro, NIE substrate RG) ⟹ α=2 nieredukowalnie aksjomatyczne (status analog c₀ #37). Ostatnia inżynieryjna droga do 'α=2 derywowane' ZAMKNIĘTA. Rekomendacje ZGŁOSZONE, edycje WSTRZYMANE."
date: 2026-06-23
type: phase_close
phase: FINAL
cycle: op-bond-order-RG-selection-2026-06-23
session: "#39"
status: 🟢 CLOSED — rekomendacje ZASTOSOWANE (user „a" 2026-06-23; sek08 rem:alpha2-pivot-status-pl + rem:amplitude-vs-density-alpha + sek10 rem:K_to_f_amplitude + STATE.md #39 + op-uv-as-ngfp NEEDS; main.tex build exit 0, 553 str.)
anti_lakatos_lock: PRESERVED
predecessor: op-Phi-field-identity-resolution-2026-06-23
verdict: "RG-NOT-SELECTED: power-counting [g_s]=−s(d−2)−(2s+2)γ dyskryminuje s=5 (α=2) jako najbardziej irrelevant; RG faworyzuje niższe s (w tym substratowe s=2/α=½); NGFP escape wymaga skrajnego γ≈−0.83 (niewiarygodny). Jedyny selektor α=2 = density-frame klasa konforemna C1–C3 (makro/geometryczna), NIE substrate-level zasada RG ⟹ α=2 nieredukowalnie aksjomatyczne na gęstości (status analog c₀ #37)."
tags: [alpha2, RG-relevance, power-counting, bond-order, no-go-positive, irreducibly-axiomatic, anti-Lakatos]
---

# Phase FINAL — CLOSE (op-bond-order-RG-selection)

> **Werdykt (value-blind, reguła Phase0 §3 — WYLICZONY, sympy 5/5):** **RG-NOT-SELECTED.**
> Power-counting daje $[g_s]=-s(d-2)-(2s+2)\gamma$. Przy $\gamma=0,\,d>2$: $[g_s]=-s(d-2)$ — **malejący
> w $s$**, więc rząd dający α=2 ($s=5$) jest **najbardziej irrelewantny** (d=4: $[g_5]=-10 < [g_2]=-4 <
> [g_0]=0$). RG-relevance **dyskryminuje** α=2 i faworyzuje niższe $s$, w tym substratowe $s=2$ (α=½).
> Marginalny jest tylko swobodny człon ($s=0$). NGFP/anomalny wymiar nie ratuje: $s=5$ marginalny
> wymaga skrajnego $\gamma\approx-0.83$ (dostrojony, niewiarygodny).

## §1 — Co cykl ROZSTRZYGNĄŁ

#38 zostawił **jedyną** drogę do „α=2 derywowane": znaleźć zasadę selekcjonującą rząd bondu $n=6$
($s=5$). Ten cykl testuje najmocniejszego kandydata — **RG-relevance** — i orzeka **negatywnie**:

| | #38 (REALIZABLE-NONCANONICAL) | **Po tym cyklu (#39)** |
|---|---|---|
| czy α=2 realizowalne? | TAK, przez bond $s=5$ (K∝ŝ¹⁰) | (niezmienione) |
| czy RG selekcjonuje $s=5$? | otwarte | **NIE** — najbardziej irrelevant; RG woli niższe $s$ |
| co RG faworyzuje? | — | $s=0$ (free, marginalny), potem $s=2$ (α=½, substrat) |
| czy NGFP ratuje? | — | tylko przy $\gamma\approx-0.83$ (skrajny, dostrojony) — niewiarygodne |
| status α=2 | aksjomat na gęstości (skwantyfikowany) | **nieredukowalnie aksjomatyczny** (analog c₀ #37) |

**Kluczowa obserwacja:** RG-relevance nie tylko „nie pomaga" — ona **aktywnie idzie w przeciwną stronę**.
Gdyby RG cokolwiek selekcjonowała w sektorze kinetycznym, byłoby to niższe $s$ (bliżej swobodnego/marginalnego),
czyli **bliżej substratowego α=½**, a nie fenomenologicznie wymaganego α=2.

## §2 — Anatomia (precyzyjnie)

- $[\hat s]=(d-2)/2+\gamma$; $O_s=\hat s^{2s}(\nabla\hat s)^2$; $[g_s]=d-[O_s]=-s(d-2)-(2s+2)\gamma$.
- **R1 (sanity):** $s=0$ (K=const, swobodny) → $[g_0]=0$ marginalny ✓. Zwarta postać $[g_s]|_{\gamma=0}=-s(d-2)$.
- **R2:** $\partial_s[g_s]=-(d-2)<0$ dla $d>2$ — monotonicznie malejący ⟹ wyższe $s$ = bardziej irrelevant.
- **R3:** d=4: $[g_0]=0,\,[g_2]=-4,\,[g_5]=-10$. $s=5$ (α=2) najgłębiej irrelevant.
- **R4:** $[g_5]=0$ wymaga $\gamma=-5/6\approx-0.83$ — anomalny wymiar rzędu jedności i ujemny;
  typowe $|\gamma|\ll1$ ⟹ niewiarygodny/dostrojony. (Dla $s=2$: $\gamma=-2/3$ — też duży.)

## §3 — Znaczenie dla programu TGP (uczciwe, dwustronne)

**Co to ZAMYKA:** ostatnią *inżynieryjną* drogę do „α=2 wyprowadzone z substratu". Po #38
(realizowalne, ale nie-kanoniczne) i #39 (RG go nie selekcjonuje, wręcz dyskryminuje), **nie istnieje
naturalna zasada RG** czyniąca α=2 wyróżnionym rzędem bondu. Jedynym selektorem pozostaje **density-frame
klasa konforemna C1–C3** — związana z geometrycznym wymogiem metryki konforemnej $g_{ij}\propto\Phi$ —
która jest **makro/geometryczna**, a NIE substrate-level zasadą RG. To dokładnie status `rem:alpha2-pivot-status`:
**α=2 = aksjomatyczna selekcja na gęstości.** Cykl ją **potwierdza i utwardza** (nieredukowalnie).

**Czego to NIE oznacza (anti-Lakatos, druga strona):**
- TGP **nie jest sfalsyfikowane** — α=2 jest fenomenologicznie wymagane (PPN, masy, Koide, soliton)
  i pozostaje spójnym aksjomatem na gęstości; fundament „jedno pole skalarne Z₂" stoi (#38: bond ŝ¹⁰ legalny).
- Pozostaje **wąska, nie-inżynieryjna** furtka: interagujący NGFP (asymptotic safety, `op-uv-as-ngfp`)
  z dużym ujemnym $\gamma\approx-0.83$ dla operatora $s=5$. To **wieloletni track ontologii UV**, nie
  ścieżka konstrukcyjna — zgłoszony jako honest open question, NIE deferred derivation (analog §4 w #37 c₀).

**Status α=2 jest teraz dokładnie analogiczny do c₀ (#37):** w zasadzie mógłby być wyprowadzony, gdyby
istniała nowa, nietrywialna struktura (Ward/normalizacja dla c₀; NGFP z dużym γ dla α=2), ale **żadna
obecna zasada go nie ustala** ⟹ traktowany jako nieredukowalne wejście aksjomatyczne. Bilans inputów
TGP bez zmian (α=2 było już liczone jako selekcja C1–C3, nie nowy parametr).

## §4 — Rekomendacje dla rdzenia (NIE wykonane — forbidden #3; ZGŁOSZONE, czekają na autoryzację)
- **P1 — `rem:alpha2-pivot-status-pl` + `rem:amplitude-vs-density-alpha`:** dopisać, że selekcja α=2
  na gęstości jest **nieredukowalna do zasady RG** — power-counting ($[g_s]=-s(d-2)$) dyskryminuje
  wymagany rząd $s=5$ (najbardziej irrelevant); jedynym selektorem jest density-frame klasa konforemna
  C1–C3 (geometryczna, $g_{ij}\propto\Phi$). Cross-ref ten cykl + #38.
- **P2 — analogia statusu:** w `rem:alpha2-pivot-status` zaznaczyć **status analog c₀ (#37)** —
  „w zasadzie wyprowadzalne tylko przez nową strukturę (NGFP z γ≈−0.83); brak ścieżki inżynieryjnej".
- **P3 — `op-uv-as-ngfp` / uv-completion:** dopisać konkretny otwarty problem — czy NGFP rodziny
  $(\hat s_i\hat s_j)^n$ generuje $\gamma\approx-5/6$ dla operatora $s=5$ (czyniąc α=2 marginalnym).
  Niski priorytet inżynieryjny, wysoki fundamentalny.
- **P4 — STATE.md:** wpis #39 (po autoryzacji): RG-relevance NIE selekcjonuje α=2; ostatnia inżynieryjna
  droga zamknięta; α=2 nieredukowalnie aksjomatyczne (analog c₀); residual = NGFP track.

## §5 — Track alternatywny (jedyna pozostała droga do α=2 jako predykcji)
Interagujący NGFP (asymptotic safety) rodziny bondów, w którym operator $s=5$ nabywa anomalny wymiar
$\gamma\approx-5/6$ i staje się marginalny/relevant. Wymaga pełnego rachunku FRG na rodzinie
$\hat s^{2s}(\nabla\hat s)^2$ — wieloletni track UV, NIE inżynieria. Zgłoszone jako honest open question.

## §6 — Anti-Lakatos (final)
- ✓ Werdykt **wyliczony** z reguły Phase0 §3 (5× test), nie wybrany — value-blind.
- ✓ Silnik (power-counting) zwalidowany **przed** twierdzeniem (R1: $s=0$ free → marginalny).
- ✓ **Wynik NEGATYWNY zgłoszony wprost** — zamyka ostatnią inżynieryjną drogę do „α=2 derywowane".
- ✓ **Obie strony:** PRO (NGFP z γ≈−0.83 mógłby) vs CONTRA (naive + umiarkowany γ czynią s=5 irrelevant).
- ✓ NIE pomylono density-frame C1–C3 (makro aksjomat) z substrate RG (forbidden #2 chronione).
- ✓ R1–R6 nieruszone (relacja EL, α=(s−1)/2 #38, gęstość kanoniczna, α=2 fenomenologicznie wymagane).
- ✓ Rdzeń NIE edytowany; budżet nowych stałych = 0; rekomendacje §4 osobno.

## §7 — Status końcowy
**🟡 CLOSED-VERDICT (analityczny, sympy 5/5).** **RG-NOT-SELECTED:** RG-relevance dyskryminuje rząd
dający α=2 ($s=5$, najbardziej irrelevant), faworyzując niższe $s$; NGFP escape wymaga skrajnego
$\gamma\approx-0.83$. Jedyny selektor α=2 = density-frame klasa konforemna C1–C3 (makro, NIE substrate RG).
⟹ **α=2 nieredukowalnie aksjomatyczne na gęstości (status analog c₀ #37); ostatnia inżynieryjna droga
do „α=2 derywowane" ZAMKNIĘTA;** residual = NGFP track (wieloletni, fundamentalny). **Edycje rdzenia
(P1–P4) WSTRZYMANE do autoryzacji.**
