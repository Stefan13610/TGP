---
title: "README — op-adiabatic-source-threshold: czy λ̃_crit jest własnością MODELU, czy PROTOKOŁU nagłego załączenia źródła?"
date: 2026-09-19
type: cycle-readme
tgp_owner: research/op-adiabatic-source-threshold-2026-09-19
folder_status: locked
status: PHASE0-LOCKED
verdict: ""
claim_status: ""
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-matter-induced-creation-2026-09-15/NEEDS.md]]"
  - "[[../op-matter-induced-creation-2026-09-15/Phase_FINAL_close.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_FINAL_close.md]]"
---

# op-adiabatic-source-threshold (2026-09-19)

Realizuje **N4** z [[../op-matter-induced-creation-2026-09-15/NEEDS.md]]
(„próg jest DYNAMICZNY, nie statyczny") oraz — wbudowane w protokół,
nie jako osobne pytania — **N1** (reguła potwierdzenia siatkowego)
i **N3** (martwy człon energetyczny w kryterium RETURN-TO-VACUUM).

## Skąd to pytanie

Dwa cykle zbudowały mapę λ̃ wyłącznie dla rodziny startów „źródło
włączone **skokowo** w t=0" i zmierzyły λ̃_crit = 0.107656 ± 0.000156.
Pre-rejestrowany mechanizm — fold statyczny 0D, λ̃_fold = 0.285770 —
**nie trafił (rozjazd 2.65×)**; kolapsy zachodzą w PIERWSZYM
overshoocie (t_end ∈ [0.60, 0.83]). Reżim λ̃ ∈ (0.1077, 0.2858),
w którym statyczne minimum NADAL ISTNIEJE, a protokół daje wyłącznie
COLLAPSE, jest nietknięty.

## Pytania

- **Q-J1 (centralne):** czy przy adiabatycznym załączaniu (rampa
  smootherstep w górę) istnieją stany osiadłe przy λ̃ > 0.1078?
- **Q-J2:** czy λ̃_crit(Δ_on) saturuje dla Δ_on ∈ {25, 100, 400},
  czy dryfuje bez granicy?
- **Q-J3 (warunkowe na Q-J1-PASS):** czy najgłębszy stan osiągnięty
  adiabatycznie przeżywa wygaszenie źródła? (powtórka Q-I2 na
  stanach głębszych niż dotąd osiągalne)

**Predykcje pre-rejestrowane — dwie, rozbieżne co do znaku:**
Q-J1 **PASS oczekiwany** (pierwsza pozytywna pre-rejestracja w tej
serii; FAIL = mocne zaskoczenie, próg byłby własnością modelu);
Q-J3 **RETURN-TO-VACUUM oczekiwane** (zgodnie z Q-I2-FAIL: zero
histerezy). Odniesienie ilościowe NIE bramkujące: λ̃_fold = 0.285770
jako domniemane ograniczenie górne λ̃_crit(Δ_on→∞).

Zakres: zależność progu od protokołu + trwałość stanów głębokich.
ZAKAZ claimów o masach leptonów, oscylonach i o samouzgodnionym ρ(ψ).

Stawka: przy Q-J1-PASS **cała mapa λ̃ obu poprzedników zostaje
przeklasyfikowana jako własność protokołu, a nie modelu.**

## Log
- 2026-09-19 — LOCK zapisany (sesja główna; autoryzacja: wybór usera
  „rozpisz lock i działaj" po analizie Q-I2-FAIL). Obliczeń zero.
  Prior art sprawdzony (`tooling/similar.py`): najwyższe podobieństwo
  0.152 = cykl bezpośrednio poprzedzający; temat adiabatycznego
  załączania nietknięty w 263 cyklach.
