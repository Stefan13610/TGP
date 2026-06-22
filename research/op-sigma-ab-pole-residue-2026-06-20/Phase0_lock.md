---
title: "Phase 0 — op-sigma-ab-pole-residue: LOCK. Pytanie: czy framework dostarcza warunek POLE-RESIDUE / kanoniczną normalizację σ_ab, który ustala C_σ (a więc κ_E) jako PREDYKCJĘ, zamiast wolnego parametru UV (op-CG4 Phase 3)? Kryteria PASS/FAIL, falsyfikatory, reguła agregatu, forbidden moves — zalockowane przed pierwszą liczbą."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-sigma-ab-pole-residue-2026-06-20
created_date: 2026-06-20
authorization: "User 2026-06-20: 'działaj z a osobny cykl zobaczmy co z tego wyjdzie'"
parent_cycles:
  - "[[../op-CG4-substrate-closure-2026-06-20/Phase3_renorm.md]] (C_σ UV-czuły = wolny parametr; rozbieżność liniowa −16/35)"
  - "[[../closure_2026-04-26/sigma_ab_pathB/results.md]] (Path B: M²=2m_s² OPE-coeff; spectrum = kontinuum od 2m_s; Bethe-Salpeter OTWARTE §5)"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-sigma-ab-pole-residue)

> **Pytanie cyklu:** op-CG4 Phase 3 dowiódł, że $C_\sigma$ (wsp. $p^2$ kanału TT σ_ab) jest **UV-czuły = wolny
> parametr**. Jedyna droga do uczynienia $\kappa_E$ predykcją: **warunek renormalizacyjny on-shell** —
> izolowany biegun (bound state) σ_ab z **dobrze określonym residuum**, które ustala $C_\sigma$ skończenie.
> **Czy framework go dostarcza?**

## §1 — Stan wejściowy (z parent)

| Fakt | Źródło |
|---|---|
| $C_\sigma$ = wolny parametr UV (rozbieżność liniowa niezerowa) | op-CG4 Phase 3 |
| σ_ab = operator złożony $(\partial\delta\hat s)(\partial\delta\hat s)$, NIE nowy d.o.f. | closure Path B T-PB.4 |
| Heredity EOM: $\Box\sigma_{ab}+2m_s^2\sigma_{ab}=S^{TT}+R^{TT}$; **$M^2{=}2m_s^2$ = coeff OPE** | closure Path B T-PB.1/2 |
| Funkcja spektralna $\rho_O(p^2)$: support $p^2\ge(2m_s)^2$ = **KONTINUUM dwucząstkowe** | closure Path B T-PB.2 |
| **Bethe-Salpeter (bound state) — OTWARTE** | closure Path B §5 |

**Napięcie do rozstrzygnięcia:** „M²=2m_s²" wygląda jak masa cząstki (biegun), ale T-PB.2 mówi, że spektrum to
kontinuum od $(2m_s)^2=4m_s^2$. Biegun bound-state musiałby leżeć PONIŻEJ progu, przy $M^2=2m_s^2<4m_s^2$
($M=\sqrt2\,m_s$). **Czy taki biegun istnieje?** — to decyduje o residuum.

## §2 — Kryteria LOCKED (PASS / FAIL)

### C-POLE — istnienie izolowanego bieguna σ_ab poniżej kontinuum
| PASS | FAIL |
|---|---|
| $\langle\sigma_{ab}\sigma_{cd}\rangle(p)$ ma **izolowany biegun** przy $p^2=-M^2$ z $M<2m_s$ (poniżej progu), wygenerowany przez interakcję substratu w kanale **spin-2** | brak izolowanego bieguna — tylko kontinuum od $4m_s^2$ (kanał spin-2 nie wiąże) |

### C-KERNEL — czy interakcja substratu działa w kanale spin-2
| PASS | FAIL |
|---|---|
| jądro Bethego-Salpetera ma **niezerową** projekcję na falę L=2 (spin-2) i jest atrakcyjne dość, by związać | jądro projektuje się na **zero** w L=2 (np. kontakt φ⁴ = czysta L=0) → brak wiązania |

### C-RESIDUE — czy residuum ustala $C_\sigma$ skończenie i scheme-independent
| PASS | FAIL |
|---|---|
| residuum przy biegunie **skończone, scheme-independent**, daje $C_\sigma$ → $\kappa_E$ liczbowo | brak bieguna ⟹ brak residuum; lub residuum reintrodukuje UV-czułość |

### C-MATCH — czy ewentualny biegun siedzi przy $M^2=2m_s^2$ (zgodność z heredity coeff)
| PASS | FAIL/USTALENIE |
|---|---|
| pozycja bieguna BS = $2m_s^2$ (zgodna z OPE coeff) | pozycja bieguna ≠ $2m_s^2$ ⟹ „M²=2m_s²" jest coeff OPE, NIE masą bieguna (USTALENIE) |

## §3 — Reguła agregatu LOCKED (value-blind)

```
C-KERNEL FAIL (brak L=2 w jądrze)            -> C-POLE FAIL -> C-RESIDUE FAIL
   => WERDYKT: framework NIE dostarcza pole-residue; C_σ pozostaje wolnym parametrem (op-CG4 Phase 3 utrzymane)
   => sektor radiacyjny: kappa_E = WOLNY PARAMETR (opcja b uczciwa)

C-KERNEL PASS i C-POLE PASS i C-RESIDUE PASS:
   => WERDYKT: pole-residue ustala C_σ -> kappa_E PREDYKCJA (liczbowa)
   => wtedy C-MATCH rozstrzyga, czy zgodne z heredity 2m_s^2

C-POLE PASS, C-RESIDUE FAIL (residuum UV-czule): => GAP (residuum nie pomaga)
```
Werdykt **wyliczony** z tej reguły w FINAL, nie wybrany.

## §4 — Forbidden moves (anti-Lakatos)
1. **Zakaz utożsamiania coeff OPE ($2m_s^2$) z fizycznym biegunem** bez dowodu istnienia bieguna w funkcji spektralnej.
2. **Zakaz fabrykowania residuum** — jeśli brak bieguna, residuum NIE istnieje (nie wymyślać normalizacji).
3. **GAP/FAIL jawnie** — wynik negatywny (brak pole-residue) jest dopuszczalny i oczekiwany jako jedno z rozstrzygnięć.
4. **Substrat = M0** (zwalidowany op-CG4: φ⁴ kontakt, klasa Isinga) jako kanoniczny; gradient-bond v2 sprawdzić jako wariant (jawnie).
5. **Zakaz rewizji rdzenia** bez autoryzacji; budżet nowych stałych 0; reuse closure Path B + op-CG4 Phase 3.

## §5 — Hipoteza wejściowa (INFORMATIONAL, nie wiąże)
Kontakt φ⁴ (M0) jest s-wave (L=0); fala spin-2 (L=2) ma zerową projekcję kontaktu (∫P₂=0) ⟹ **przewidywany brak
wiązania spin-2** ⟹ brak bieguna ⟹ C_σ wolny (op-CG4 Phase 3 utrzymane). Werdykt rozstrzygnie rachunek (Phase 1).
Ryzyko obalenia: interakcja pochodna (gradient-bond) mogłaby dać niezerowe L=2 — sprawdzić.

## §6 — Status
**🔒 PHASE 0 LOCKED (2026-06-20).** Kryteria C-POLE/C-KERNEL/C-RESIDUE/C-MATCH, reguła agregatu, forbidden moves —
zalockowane przed pierwszą liczbą. Przejście do Phase 1 (spektrum ⟨σσ⟩ + projekcja jądra na L=2).
