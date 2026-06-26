---
title: "op-CG34-continuum-closure — domknięcie CG-3 (homogenizacja Φ_B→Φ w H¹) i CG-4 (K_hom=K_TGP)"
type: research_cycle
status: "🟢 CLOSED — CG-3 ZAMKNIĘTY NUM; CG-4 PARTIAL (znaczący postęp) 2026-06-14"
phase: FINAL
folder_status: active
created_date: 2026-06-14
closed_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'działaj z domknięciem CG-3/CG-4'"
target: "dodatekQ CG-3 [OTWARTY] + CG-4 [OTWARTY] (residual cyklu op-Csigma-lattice-MC)"
verdict: "CG-3=ZAMKNIĘTY NUM; CG-4=PARTIAL (c*/β=γ/K_hom-forma zamknięte; α=2↔K(φ) do dopięcia)"
anti_lakatos_lock: PRESERVED
---

# op-CG34-continuum-closure (🟢 CLOSED)

> **Cel:** domknąć CG-3 (zbieżność $\Phi_B\to\Phi$ w $H^1$, homogenizacja) i CG-4 (identyfikacja
> $K_{\rm hom}=K_{\rm TGP}$) — analityczny residual zidentyfikowany w cyklu
> [[../op-Csigma-lattice-MC-2026-06-14/Phase_FINAL_close.md]] (sektor radiacyjny).

## Werdykt
- **CG-3 = ZAMKNIĘTY NUM.** Homogenizacja $\Phi_B\to\Phi$ w $H^1$ potwierdzona (5/5): normy H¹ Cauchy
  $\{1.34,0.38,0.11\}$ monotonicznie maleją (naprawiono bug prior ‖ΔΦ‖=0.000000); rate $L_B^{-1.43}\ge$ bound A5.
- **CG-4 = PARTIAL.** ZAMKNIĘTE: **c*>0 stabilne** (red-flag $c_*\to0$ rozwiązany — był artefaktem
  $\langle|\nabla\Phi|^2\rangle$), **β=γ** (algebraicznie), **K_hom-forma=K_IR** (A2+CG-2).
  DO-PINNIĘCIA: **α=2 ↔ mikroskopowy $K(\phi)$** — samodzielna algebra ($\alpha_{\rm eff}=s-1$) ujawniła
  **niespójność lematu A3** (substrat $Z\propto\phi^2\Rightarrow\alpha=0$, nie 2); zgłoszona, nie sfabrykowana.

## Kluczowe ustalenia (uczciwe)
1. **Naprawa numeryk prior:** H¹ convergence (zera→niezerowe), c* (artefakt→stabilne), CG-4 (R²=0.33→algebra).
2. **Red-flag c*→0 rozwiązany:** błąd metody ($K_1=\langle|\nabla\Phi|^2\rangle$); poprawnie ze structure factor $Z>0$ stabilne.
3. **Substrat $-J(\phi_i\phi_j)^2$ patologiczny dla pure-MC:** λ małe→runaway, λ duże→frozen (brak okna scale-separated).
4. **Niespójność α=2 (lemat A3):** premisa $Z\sim\phi^2$ vs konkluzja $K_1\sim1/\Phi$; krytyczne do rozstrzygnięcia w rdzeniu.

## Pliki cyklu
- [[./Phase0_lock.md]] — LOCK kryteriów C3-A..D, C4-A..D, forbidden moves, diagnoza błędów prior.
- [[./Phase1_engine.md]] + `Phase1_engine.py` — silnik φ⁴, κ_c, **koercywność c*>0** (5/5), diagnoza substratu.
- [[./Phase2_cg3.md]] + `Phase2_cg3.py` — **CG-3 ZAMKNIĘTY NUM** (homogenizacja H¹, 5/5).
- [[./Phase3_cg4.md]] + `Phase3_cg4.py` — **CG-4 PARTIAL** (c*/β=γ/K_hom; α=2 niespójność).
- [[./Phase_FINAL_close.md]] — werdykt agregatowy + rekomendacje rdzenia (§5).

## Rekomendacje dla rdzenia (zgłoszone, NIE wykonane — forbidden move #5)
1. dodatekQ status: CG-3 [OTWARTY]→**[ZAMKNIĘTY NUM]**.
2. dodatekQ status: CG-4 [OTWARTY]→**[CZĘŚCIOWY NUM]**.
3. dodatekQ2 lemat A3: ujednolicić $Z(\phi)$↔α=2 (sek10 vs A3).
4. `substrate_mc_cg3.py`: SUPERSEDED (bug K₁ → artefakt c*→0).
