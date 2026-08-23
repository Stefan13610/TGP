---
title: "op-c0-derivation-from-substrate (2026-06-22, #37) — czy c₀ (sprzężenie σ, C(ψ=1)) jest wyprowadzalne z substratu, czy jest wolnym parametrem UV (jak C_σ, #33)?"
date: 2026-06-22
type: research-cycle
status: 🟡 CLOSED-VERDICT — c₀ = WOLNY PARAMETR UV (LOCK §4 wiersz 1, sympy 5/5); rekomendacje rdzenia WSTRZYMANE do autoryzacji
session: "#37"
parent_cycles:
  - "[[../op-CG4-substrate-closure-2026-06-20/Phase_FINAL_close.md]]  (#33)"
  - "[[../op-sigma-ab-pole-residue-2026-06-20/README.md]]  (#34)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]]  (recovery §3.6)"
supersedes_candidate: "[[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]]"
coordination: "[[../../STATE.md]]"
tags: [c0, sigma-coupling, C-sigma, kappa-E, emergent-metric, UV-divergence, predictive-vs-fit, anti-Lakatos]
---

# op-c0-derivation-from-substrate

## Cel
Rozstrzygnąć **value-blind**, czy `c₀ = C(ψ=1)` (wiodące sprzężenie strain-composite $\sigma^{ij}$
w metryce emergentnej $g_{\rm eff}^{ij}$) jest **scheme-independent, substrate-derivable** współczynnikiem
geometrycznym — co przywróciłoby sektorowi grawitacyjnemu TGP falsyfikowalną predykcję amplitudy/fazy GW
po falsyfikacji f(ψ)=(4−3ψ)/ψ (5σ, GWTC-3) — **czy** dziedziczy status WOLNEGO PARAMETRU UV od $C_\sigma$
(#33: rozbieżność liniowa −16/35; #34: brak bieguna spin-2).

## Dlaczego ten cykl (priorytet — wyznaczony 2026-06-22)
Wskazany jako **najważniejsza rzecz do zbadania w TGP_v1** po domknięciu α=2 (#36):
to **jedyny sektor już empirycznie sfalsyfikowany**, a jego „recovery" opiera okno zgodności GWTC-3
na strojeniu `c₀·κ_σ=4/3`. Po tym, jak #33/#34/#35 udowodniły, że κ_E ($\sim C_\sigma$) jest wolnym
parametrem UV, mocny prior mówił, że c₀ — **normalizacja tego samego operatora ∂ŝ∂ŝ** — podzieli ten los.
Dodatkowo istniał **konflikt cykli**: `op-c0-derivation-from-substrate-2026-05-09` ogłosił c₀≈4π
„STRUCTURAL DERIVED", lecz **przed** dowodem #33. Cykl rozstrzyga ten konflikt.

## Pliki cyklu
- `Phase0_lock.md` — fakty F1–F5, pytanie centralne, rozróżnienia R1–R5, reguła werdyktu §4, zakres, forbidden moves, falsyfikator dwustronny, **lean jawny** (pesymistyczny).
- `Phase0_balance.md` — gate (cykl dotyka statusów epistemicznych: downgrade 2026-05-09 + §3.6.8).
- `Phase1_sympy.py` / `Phase1_sympy.txt` — walidacja silnika (V1 −16/35, V2 0) + 3 checki strukturalne (C1–C3).
- `Phase1_results.md` — rekoncyliacja: Argument A (UV) + Argument B (proweniencja matching/kalibracja); werdykt wyliczony.
- `Phase_FINAL_close.md` — werdykt + rekomendacje P1–P5 + track alternatywny + anti-Lakatos.

## Werdykt (value-blind, LOCK §4 wiersz 1, sympy 5/5)
**`c₀` = WOLNY PARAMETR UV** — ten sam operator/normalizacja co $C_\sigma$ ($\kappa_E$). Dwa niezależne,
zbieżne argumenty: (A) współczynnik $p^2$ operatora ∂ŝ∂ŝ w kanale spin-2 jest liniowo rozbieżny w odcięciu
(−16/35, #33), brak tożsamości Warda/bieguna (#34) ⟹ scheme-dependent; (B) „c₀≈4π" z 2026-05-09 pochodzi
z **matchingu** $\xi_{\rm eff}$ (R3, ≠ predykcja) + **kalibracji GW150914**, a iloczyn $4\pi\cdot\frac1{3\pi}=\frac43$
jest **algebraicznie trywialny** (post-hoc-konstruowalny).

## Konsekwencja
Okno recovery „c₀·κ_σ=4/3" to **podwójne strojenie** ⟹ sektor grawitacyjno-GW TGP **nie dostarcza
falsyfikowalnej predykcji amplitudy/fazy**. **Niezmienione (predykcje strukturalne):** 2 TT + breathing
mode (smoking gun 3G), $c_{\rm GW}=c$, $m_\sigma^2=2m_s^2$. Bilans parametrów **3** (c₀ to ta sama stała
radiacyjna co C_σ, nie czwarty parametr).

## Anti-Lakatos
Kryteria + lean zalockowane przed liczbami; werdykt **wyliczony**; silnik zwalidowany (V1/V2); c₀ NIE
sfabrykowane (trywialność 4/3 pokazana); kalibracja wykluczona jako dowód; R1–R5 chronione; budżet stałych 0;
rdzeń/cykl 2026-05-09/FOUNDATIONS **nie edytowane** — rekomendacje P1–P5 czekają na autoryzację.

## Status
🟡 **CLOSED-VERDICT** (analityczny). Rekomendacje rdzenia (P1–P5: FOUNDATIONS §3.6.8/§3.6.10, sek08
`rem:sigma-params`, downgrade 2026-05-09, PREDICTIONS_REGISTRY) **WSTRZYMANE do autoryzacji usera**.
