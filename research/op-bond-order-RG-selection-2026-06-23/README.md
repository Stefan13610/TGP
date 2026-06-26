---
title: "op-bond-order-RG-selection — czy RG-relevance selekcjonuje rząd bondu dający α=2? (następczy do #38)"
date: 2026-06-23
type: research_cycle
status: "🟢 CLOSED — RG-NOT-SELECTED (sympy 5/5); rekomendacje ZASTOSOWANE (sek08 ×2 + sek10 + STATE #39 + op-uv-as-ngfp; build exit 0, 553 str.)"
phase: FINAL
folder_status: active
created_date: 2026-06-23
session: "#39"
authorization: "User 2026-06-23: 'b' (cykl następczy nad zasadą selekcji rzędu bondu)"
predecessor: op-Phi-field-identity-resolution-2026-06-23
verdict: "RG-NOT-SELECTED — power-counting dyskryminuje rząd dający α=2 (s=5, najbardziej irrelevant); jedyny selektor = density-frame klasa konforemna C1–C3 (NIE substrate RG) ⟹ α=2 nieredukowalnie aksjomatyczne (analog c₀ #37)"
anti_lakatos_lock: PRESERVED
---

# op-bond-order-RG-selection (🟡 CLOSED-VERDICT — RG-NOT-SELECTED)

> **Pytanie:** Czy RG-relevance (power-counting na rodzinie bondów Z₂ K(ŝ)∝ŝ^{2s}) selekcjonuje
> rząd dający α=2 (s=5)? **Odpowiedź (sympy, value-blind, 5/5): NIE — RG go dyskryminuje
> (najbardziej irrelevant), faworyzując niższe s. α=2 pozostaje nieredukowalnie aksjomatyczne.**

## Kontekst
Następczy do `op-Phi-field-identity-resolution-2026-06-23` (#38), który orzekł α=2 = REALIZABLE-NONCANONICAL
i zostawił **jedyną** drogę do „α=2 derywowane": zasada selekcjonująca rząd bondu n=6 (s=5). Ten cykl
testuje najmocniejszego kandydata — RG-relevance — value-blind.

## Werdykt
- **[g_s] = −s(d−2) − (2s+2)γ**; przy γ=0: **[g_s] = −s(d−2)** — malejący w s (d>2).
- **d=4, γ=0:** [g₀]=0 (free, marginalny), [g₂]=−4 (α=½, substrat), **[g₅]=−10 (α=2, najbardziej irrelevant)**.
- **RG faworyzuje NIŻSZE s** — bliżej swobodnego/marginalnego, czyli bliżej substratowego α=½, NIE α=2.
- **NGFP escape:** s=5 marginalny wymaga γ≈−0.83 (skrajny, ujemny, dostrojony) — niewiarygodny.
- ⟹ **RG-NOT-SELECTED.** Jedyny selektor α=2 = **density-frame klasa konforemna C1–C3** (makro,
  geometryczna, g_ij∝Φ) — to NIE substrate-level zasada RG, lecz przeformułowanie znanej selekcji aksjomatycznej.

## Znaczenie
- **Zamyka ostatnią inżynieryjną drogę** do „α=2 wyprowadzone z substratu". α=2 = **nieredukowalnie
  aksjomatyczne na gęstości** — status dokładnie **analog c₀ (#37)**: w zasadzie wyprowadzalne tylko przez
  nową, nietrywialną strukturę (NGFP z dużym γ), brak ścieżki konstrukcyjnej.
- **NIE falsyfikacja:** α=2 fenomenologicznie wymagane, fundament „jedno pole skalarne Z₂" stoi (#38).
- **Residual (jedyna furtka):** interagujący NGFP rodziny (ŝᵢŝⱼ)ⁿ z γ≈−5/6 dla operatora s=5 — wieloletni
  track UV (`op-uv-as-ngfp`), nie inżynieria.

## Pliki
- [[./Phase0_lock.md]] — LOCK; aparat power-counting; reguła decyzyjna value-blind.
- [[./Phase1_RG_selection_sympy.py]] + [[./Phase1_output.txt]] + [[./Phase1_results.json]] — derywacja (5/5 PASS).
- [[./Phase_FINAL_close.md]] — werdykt + anatomia + znaczenie dwustronne + rekomendacje rdzenia (§4–5).

## Rekomendacje rdzenia (zgłoszone, NIE wykonane — czekają na autoryzację)
1. `rem:alpha2-pivot-status-pl` + `rem:amplitude-vs-density-alpha` — α=2 nieredukowalne do zasady RG; selektor = density-frame C1–C3.
2. `rem:alpha2-pivot-status` — status analog c₀ (#37): wyprowadzalne tylko przez NGFP (γ≈−0.83), brak ścieżki inżynieryjnej.
3. `op-uv-as-ngfp` / uv-completion — otwarty problem: NGFP rodziny (ŝᵢŝⱼ)ⁿ generujący γ≈−5/6 dla s=5.
4. STATE.md — wpis #39.
