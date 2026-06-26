---
title: "Phase 0 — op-bond-order-RG-selection: LOCK. Następczy do op-Phi-field-identity-resolution (#38). Pytanie: czy zasada RG-relevance (skalowanie operatorów na rodzinie bondów (ŝᵢŝⱼ)ⁿ / K∝ŝ^{2s}) SELEKCJONUJE rząd dający α=2 (s=5), czy nie? Jeśli nie — α=2 nieredukowalnie aksjomatyczne (tylko density-frame conformal C1–C3). Value-blind, anti-Lakatos."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-bond-order-RG-selection-2026-06-23
created_date: 2026-06-23
session: "#39"
authorization: "User 2026-06-23: 'b' (otwórz cykl następczy nad zasadą selekcji rzędu bondu n=6)"
origin: "op-Phi-field-identity-resolution-2026-06-23 §5 (#38): α=2 = REALIZABLE-NONCANONICAL; otwarty problem = zasada selekcjonująca rząd bondu n=6"
target: "rodzina bondów Z₂ (ŝᵢŝⱼ)ⁿ → K(ŝ)∝ŝ^{2s}; α_density=(s−1)/2; sek08 rem:alpha2-pivot-status (klasa konforemna C1–C3); op-uv-as-ngfp (NGFP escape hatch)"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-bond-order-RG-selection)

## §1 — Pytanie wiodące

> **Czy istnieje ZASADA RG (relevance/marginality skalowania operatorów) selekcjonująca rząd bondu
> substratu dający α=2 na gęstości (s=5, K∝ŝ¹⁰) — czy RG-relevance go NIE selekcjonuje (a wręcz
> dyskryminuje)?** Jeśli NIE — jedynym selektorem α=2 pozostaje density-frame klasa konforemna C1–C3
> (makro/geometryczna), co czyni α=2 **nieredukowalnie aksjomatycznym** (status analog c₀, #37).

**Kontekst (#38, sympy 5/5):** α_density(s)=(s−1)/2; α=2 ⟺ s=5 ⟺ K(ŝ)∝ŝ¹⁰ (bond rząd n=6). Bond
dopuszczalny (Z₂-skalarny), ale nie-kanoniczny (v2 to s=2→α=½), niezależny od V. #38 zostawił JEDYNĄ
drogę do „α=2 derywowane": zasada selekcjonująca rząd 6. Ten cykl ją testuje — wprost i value-blind.

## §2 — Aparat (LOCKED, standardowy power-counting RG)

Operator kinetyczny rodziny: $O_s = \hat s^{2s}(\nabla\hat s)^2$ (substratowy człon, $K(\hat s)\propto\hat s^{2s}$).
- Wymiar pola z **kanonicznego** członu swobodnego $(\nabla\hat s)^2$ (Gaussian FP, $d$ wymiarów czasoprzestrz.):
  $[\hat s]=(d-2)/2$ (dopuszczamy anomalny wymiar $\gamma$: $[\hat s]=(d-2)/2+\gamma$).
- $[O_s] = (2s+2)[\hat s] + 2$; sprzężenie $g_s$ ma $[g_s] = d - [O_s]$.
- **Relevance:** $[g_s]>0$ relewantny (selekcjonowany w IR), $=0$ marginalny (kandydat na klasę
  konforemną/scale-inv.), $<0$ irrelewantny (tłumiony).
- Relacja do α: $\alpha_{\rm density}=(s-1)/2$ (#38, op-A3). $s=5\Leftrightarrow\alpha=2$; $s=2\Leftrightarrow\alpha=\tfrac12$ (substrat); $s=0\Leftrightarrow\alpha=-\tfrac12$ (swobodny $K$=const).

## §3 — Falsyfikatory / reguła decyzyjna LOCKED (value-blind)

| ID | Test | Reguła |
|---|---|---|
| **R1** | Wyprowadzić $[g_s](d,\gamma)$ z power-counting | sympy; SANITY: $s=0$ (free) → marginalny ($[g_0]=0$ przy γ=0) |
| **R2** | Monotoniczność $[g_s]$ w $s$ (przy γ=0, $d>2$) | rosnący vs malejący — czy wyższe $s$ bardziej czy mniej relewantne? |
| **R3** | Czy $s=5$ (α=2) jest selekcjonowany (relevant/marginal) vs $s=2$ (α=½)? | $[g_5]$ vs $[g_2]$ vs $[g_0]$ value-blind |
| **R4** | NGFP/anomalny wymiar escape: jakie $\gamma$ czyni $s=5$ marginalnym? | rozwiązać $[g_5]=0$ dla γ; ocenić wiarygodność (mały/duży) |

**Reguła agregatu (werdykt WYLICZONY):**
- **RG-SELECTED** — jeśli $s=5$ jest unikalnie relevant/marginal (a $s=2,0$ nie) ⟹ α=2 ma podstawę RG ⟹ upgrade „derywowane".
- **RG-NOT-SELECTED** — jeśli $s=5$ NIE jest relevant/marginal (np. najbardziej irrelewantny; RG faworyzuje niższe $s$) ⟹ α=2 bez podstawy RG ⟹ pozostaje aksjomatem (tylko density-frame C1–C3).
- **NGFP-CONDITIONAL** — jeśli $s=5$ marginalny tylko przy konkretnym γ; werdykt zależny od wiarygodności tego γ (mały γ → wiarygodne; duży/dostrojony → nie).

## §4 — Forbidden moves (anti-Lakatos)
1. **Zakaz przesądzania** SELECTED/NOT — sympy liczy $[g_s]$ value-blind.
2. **Zakaz mylenia** density-frame klasy konforemnej C1–C3 (makro, JUŻ znana selekcja aksjomatyczna) z substrate-level zasadą RG. Jeśli jedyny selektor to C1–C3 — to NIE nowa zasada, lecz przeformułowanie aksjomatu.
3. **Zakaz rewizji rdzenia** — werdykt = rekomendacja.
4. **Obie strony:** za RG-selekcją (NGFP/anomalny wymiar mógłby) i przeciw (naive power-counting faworyzuje niższe s).
5. Relacje #38 (α=(s−1)/2) i op-A3 (EL) — POPRAWNE, nie podważam.
6. Jeśli NOT-SELECTED — uczciwie: to wynik NEGATYWNY, zamykający ostatnią inżynieryjną drogę do „α=2 derywowane".

## §5 — Oczekiwanie (INFORMATIONAL, nie wiąże)
Naive power-counting: $[g_s]=-s(d-2)$ (do weryfikacji sympy) — malejący w $s$ dla $d>2$ ⟹ $s=5$ NAJBARDZIEJ
irrelewantny, $s=0$ (free) marginalny. Prawdopodobny werdykt: **RG-NOT-SELECTED** — RG-relevance dyskryminuje
α=2 (faworyzuje niższe s, w tym substratowe s=2/α=½). α=2 zostaje aksjomatem (density-frame C1–C3). Werdykt rozstrzygnie sympy.

## §6 — Status
**🔒 PHASE 0 LOCKED.** Anti-Lakatos aktywny. Przejście do derywacji sympy (value-blind).
