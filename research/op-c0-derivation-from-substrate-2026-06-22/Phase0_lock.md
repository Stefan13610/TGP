---
title: "Phase 0 LOCK — op-c0-derivation-from-substrate (2026-06-22, #37). Czy c₀ (sprzężenie σ, C(ψ=1)) jest wyprowadzalnym współczynnikiem geometrycznym, czy wolnym parametrem UV — tym samym co C_σ (#33)?"
date: 2026-06-22
type: phase_lock
status: 🔒 LOCKED (przed jakąkolwiek liczbą / edycją rdzenia)
phase: 0
cycle: op-c0-derivation-from-substrate-2026-06-22
session: "#37"
anti_lakatos_lock: ACTIVE
supersedes_candidate: "[[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] (STRUCTURAL DERIVED heuristic — predaty #33)"
parent_cycles:
  - "[[../op-CG4-substrate-closure-2026-06-20/Phase_FINAL_close.md]]  (#33 — C_σ UV-czuły)"
  - "[[../op-sigma-ab-pole-residue-2026-06-20/README.md]]  (#34 — brak bieguna spin-2)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]]  (recovery framework, §3.6)"
coordination: "[[../../STATE.md]]"
tags: [c0, sigma-coupling, C-sigma, kappa-E, emergent-metric, radiation-sector, UV-divergence, predictive-vs-fit, anti-Lakatos]
---

# Phase 0 — LOCK (kryteria zalockowane przed pierwszą liczbą)

> Reguła nadrzędna (anti-Lakatos): kryteria werdyktu, zakres, falsyfikator dwustronny
> oraz **jawny lean** są ustalone **przed** jakąkolwiek liczbą i przed jakąkolwiek edycją
> rdzenia / cykli zależnych. Werdykt jest **wyliczany** z reguły §4, nie wybierany.
> Wynik dowolnego znaku (c₀ DERIVED **lub** c₀ = wolny parametr) jest dopuszczalny i wartościowy.

---

## §1 — Fakty referencyjne (zalockowane, z cykli nadrzędnych)

- `F1` (#33 op-CG4, Phase 3, sympy+num): współczynnik **p²** operatora złożonego
  $O_{ab}=\partial\hat s\,\partial\hat s$ **w kanale TT spin-2** ma **niezerową rozbieżność LINIOWĄ
  UV** — współczynnik kątowy $\int_{-1}^{1}(1-\mu^2)^2(4\mu^2-1)\,d\mu=-16/35\neq0$
  (analit. liniowy $-0.002895$ = num $-0.002891$). ⟹ **nie istnieje scheme-independent
  continuum** gołego operatora; jego normalizacja jest **UV-czuła**.
- `F2` (#34 op-sigma-ab-pole-residue): wolne $\langle\sigma\sigma\rangle$ to **kontinuum**
  (cięcie od $p^2=-4m_s^2$), **brak izolowanego bieguna spin-2** ($\int_{-1}^1 P_2\,dx=0$;
  s-wave nie wiąże d-wave) ⟹ **brak residuum on-shell**, którym można by ustalić normalizację.
- `F3` (#34/#35, zapisane w rdzeniu): $\kappa_E\equiv 8\pi G_0 C_\sigma\sigma_0^2/c^3$
  (≡ stała kinetyczna $C_\sigma$) = **nieredukowalny wolny parametr UV**; bilans parametrów TGP = **3**.
- `F4` (struktura recovery, FOUNDATIONS §3.6.1): metryka emergentna ma
  $g_{\rm eff}^{ij}=\delta^{ij}B(\psi)+\sigma^{ij}\,C(\psi)/(\Phi_0^2c^2)$, gdzie
  $\sigma^{ij}=(\partial^i\Phi)(\partial^j\Phi)-\tfrac13\delta^{ij}(\nabla\Phi)^2$ jest
  **gradient-strain composite** (level-0, OP-7 T2). Definicja: **$c_0\equiv C(\psi=1)$**.
- `F5` (cykl 2026-05-09): c₀≈4π „STRUCTURAL DERIVED (heuristic)" via $\xi_{\rm eff}=4\pi G\Phi_0^2$
  (OP-7 T3.4) + Path A→B; κ_σ≈1/(3π) „plausibility"; iloczyn c₀·κ_σ=4/3 EXACT =
  cel zero-β_ppE Phase 4. **Caveaty własne tego cyklu:** Path A→B „może mieć pominięty O(1)";
  κ_σ „NIE explicit derivation"; wartość „calibrated" używa GW150914; cykl **predaty F1/F2/F3**.

---

## §2 — Pytanie centralne cyklu (jedno zdanie)

> **Czy $c_0=C(\psi=1)$ — wiodące sprzężenie strain-composite $\sigma^{ij}$ w $g_{\rm eff}^{ij}$ —
> jest scheme-independent, substrate-derivable współczynnikiem geometrycznym,
> CZY jest normalizacją tego SAMEGO operatora $\partial\hat s\partial\hat s$, którego współczynnik
> $p^2$ (= $C_\sigma$) udowodniono (F1) jako UV-czuły wolny parametr — a więc $c_0$ dziedziczy
> status wolnego parametru?**

Konsekwencja decyzyjna (do wyliczenia, NIE z góry przesądzona):
- jeśli **c₀ derivable** (chroniony przed rozbieżnością F1, np. jako bezwymiarowy **stosunek**
  / ustalony tożsamością Warda / czynnik czysto drzewiasty) → recovery ma **≥1 genuine predykcję**,
  okno c₀·κ_σ=4/3 staje się testowalne → **pozytyw**;
- jeśli **c₀ = ten sam wolny parametr UV** → okno „c₀·κ_σ=4/3" to **podwójne strojenie** (c₀ wolne
  ∧ κ_E wolne) → sektor grawitacyjno-GW TGP **NIE-PREDYKTYWNY** post-falsyfikacja; cykl 2026-05-09
  do **downgrade** (STRUCTURAL DERIVED → SUPERSEDED/heurystyka-kalibracja), §3.6.8 FOUNDATIONS do korekty.

---

## §3 — Rozróżnienia zalockowane (NIE downgradować, NIE fabrykować)

- `R1`: **m_σ² = 2 m_s²** (GW4, OPE-invariant) = masa, **LOCKED**. NIE jest c₀ ani C_σ. NIE rusza cykl.
- `R2`: **liczba modów / symetria** (2 TT + breathing; c_GW=c; brak modów wektorowych) = niezależne
  od c₀ i κ_E (F3/#35 R2/R4). NIE rusza cykl.
- `R3`: **thm:amplitude-matching** ($\xi_{\rm eff}=4\pi G_0\sigma_0\Phi_0$) = **warunek dopasowania**,
  nie predykcja. „4π" w F5 pochodzi z tego matchingu — sam w sobie NIE jest wyprowadzeniem c₀.
- `R4`: **α=2 / Φ=⟨ŝ²⟩ gęstość** = selekcja na gęstości (#31/#32/#36). NIE reaktywować Opcji B; cykl
  NIE zależy od rozstrzygnięcia α.
- `R5`: **relacja struktury** $\sigma_{ab}$ jako gradient-strain composite (OP-7 T2) = LOCKED definicja;
  cykl bada **normalizację** ($C(\psi)$), nie istnienie struktury.

---

## §4 — Reguła werdyktu (wyliczana, value-blind)

Niech `D` = „rozbieżność liniowa UV współczynnika $p^2$ operatora $\partial\hat s\partial\hat s$
w kanale spin-2" (F1). Niech `c₀` operacyjnie = normalizacja sprzężenia tegoż operatora w $g_{\rm eff}^{ij}$.

| # | Warunek (wyliczany) | Werdykt c₀ | Konsekwencja sektora |
|---|---|---|---|
| 1 | c₀ jest tym samym współczynnikiem $p^2$ co $C_\sigma$ (ta sama subtrakcja) **i** dziedziczy `D` | **WOLNY PARAMETR UV** | NIE-PREDYKTYWNY (podwójne strojenie c₀∧κ_E) |
| 2 | c₀ jest bezwymiarowym **stosunkiem** / ustalony **tożsamością** (Ward/geometria) **anulującą** `D` | **DERIVED (scheme-indep.)** | recovery ma ≥1 genuine predykcję |
| 3 | c₀ = czynnik **drzewiasty** (tree-level), `D` wchodzi dopiero w pętli i jest **renormalizowalny osobno** | **DERIVED-CONDITIONAL** | predyktywny modulo schemat; zgłosić warunki |
| 4 | rozstrzygnięcie 1 vs 2/3 wymaga rachunku niedostarczonego tą drogą | **UNDERDETERMINED** | status zachowany; zgłosić jako gap z dowodem |

**Lean jawny (anti-Lakatos, pesymistyczny — ten, który teorii szkodzi):** wiersz **1**
(c₀ = ten sam wolny parametr UV; sektor NIE-PREDYKTYWNY). Lean zadeklarowany, by werdykt pozytywny
(wiersze 2/3) wymagał **mocniejszego** dowodu niż negatywny.

---

## §5 — Zakres (kotwice; NIE wyczerpujące)

1. `op-c0-derivation-from-substrate-2026-05-09` — Phase1_results + Phase_FINAL (kandydat do downgrade).
2. `op-CG4-substrate-closure-2026-06-20` Phase3_renorm.md — dokładna definicja operatora, którego $C_\sigma$ jest UV-czuły; sprawdzić, czy c₀ to ta sama normalizacja.
3. FOUNDATIONS §3.6.1 (definicja C(ψ)), §3.6.8 (c₀ „framework-derivable, deferred" — kandydat do korekty), §3.6.10.1–3 (c₀=4π, κ_σ, iloczyn 4/3).
4. OP-7 T3.4 (`op7/OP7_T3_results.md`) — źródło $\xi_{\rm eff}=4\pi G\Phi_0^2$; ustalić, czy to matching (R3) czy derivation.
5. sek08 `rem:sigma-params` / `rem:sigma-Csigma-free` — referencja docelowej spójności statusu.
6. PREDICTIONS_REGISTRY GW7 (C_σ FREE-PARAMETER) + ewentualny nowy wiersz dot. c₀.

---

## §6 — Forbidden moves (anti-Lakatos)

1. NIE fabrykować „c₀ przewidziany" — żadnego dopasowania c₀ do 4/3/κ_σ ani do GW150914 jako „derivation".
2. NIE downgradować R1 (m_σ²=2m_s²) ani R2 (modi/symetria) ani R5 (struktura σ_ab).
3. NIE mylić **matchingu** (R3, ξ_eff) z **wyprowadzeniem** normalizacji c₀.
4. Silnik analityczny/sympy ZWALIDOWANY przed użyciem (reprodukcja znanego wyniku F1: −16/35).
5. Edycje rdzenia / cyklu 2026-05-09 / FOUNDATIONS **TYLKO** addytywne, **po autoryzacji usera**, po Phase0_balance gate. Budżet nowych stałych = 0.
6. Wynik negatywny (lub pozytywny) zgłosić **wprost**, wyliczony z §4.
7. NIE wybierać modelu/ścieżki minimalizującej dryf werdyktu; wybór z argumentu strukturalnego.

---

## §7 — Falsyfikator cyklu (dwustronny)

- Cykl **FALSE-NEGATIVE**, jeśli zaklasyfikuje c₀ jako wolny parametr, gdy istnieje jawna
  tożsamość / struktura stosunku anulująca `D` (wiersz 2/3 §4) — pominięta.
- Cykl **FALSE-POSITIVE**, jeśli ogłosi c₀ DERIVED, opierając się na (i) matchingu ξ_eff (R3),
  (ii) kalibracji GW150914, lub (iii) post-hoc iloczynie 4π·1/(3π)=4/3 — wszystkie wykluczone jako dowód.
