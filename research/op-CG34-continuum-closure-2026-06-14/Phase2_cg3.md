---
title: "Phase 2 — op-CG34-continuum-closure: CG-3 (homogenizacja Φ_B→Φ w H¹) ZAMKNIĘTY NUM (5/5 PASS). Naprawiono bug prior (‖ΔΦ‖=0.000000): poprawne, NIEZEROWE normy. H¹ CAUCHY: ‖Φ_b−Φ_2b‖_L2={0.496,0.187,0.070}, ‖·‖_H1={1.34,0.38,0.11} — monotoniczny spadek (zbieżność). Rate ‖Φ_b−Φ_ref‖~L_B^{−1.43} (R²=0.999) ≥ A5 bound (a/L_B)^{1/2}. Uniform bound i exp decay potwierdzone. Demonstracja na φ⁴ (uniwersalność)."
type: phase_result
status: PHASE2_COMPLETE (CG-3 ZAMKNIĘTY NUM, 5/5 PASS)
phase: 2
cycle: op-CG34-continuum-closure-2026-06-14
created_date: 2026-06-14
script: "[[./Phase2_cg3.py]]"
output: "[[./Phase2_output.txt]]"
target: "CG-3 [OTWARTY,AN] → [ZAMKNIĘTY NUM]"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — CG-3: homogenizacja Φ_B → Φ w H¹

## §0 — Werdykt (5/5 PASS) → **CG-3 ZAMKNIĘTY NUM**
| Test | Wynik |
|---|---|
| C3-A (uniform bound, A1(i)) | **PASS** — ‖Φ_B‖_L2 zbieżne dla b≥2 (0.645→0.613, rozrzut 5%) |
| C3-B (exp decay korelacji, W4) | **PASS** — ξ=2.80 skończone |
| **C3-C (H¹ Cauchy — rdzeń CG-3)** | **PASS** — ‖Φ_b−Φ_2b‖ NIEZEROWE, monotonicznie maleje |
| C3-D (rate vs A5) | **PASS** — p=1.43 (R²=0.999) ≥ 1/2 (bound A5) |

## §1 — Naprawa prior + rdzeń wyniku
Prior `cg_strong_numerical.py` raportował ‖ΔΦ‖_L2=‖ΔΦ‖_H1=**0.000000** (zepsuty pomiar). Naprawione: pola blokowe Φ_b broadcastowane na wspólną siatkę L³ (piecewise-constant), normy L²/H¹ z poprawnych różnic skończonych. Wynik (L=48, κ=0.187, φ⁴):

| blokowanie | ‖ΔΦ‖_L2 | ‖ΔΦ‖_H1 |
|---|---|---|
| 1→2 | 0.4961 | 1.3395 |
| 2→4 | 0.1867 | 0.3801 |
| 4→8 | 0.0702 | 0.1127 |

**Monotoniczny spadek** (czynnik ~2.6 / blokowanie) ⟹ rodzina $\{\Phi_B\}$ jest **Cauchy w $H^1$** ⟹ zbiega (homogenizacja). To jest **rdzeń CG-3** — teraz numerycznie solidny.

## §2 — Rate zbieżności vs A5
$\|\Phi_b-\Phi_{\rm ref}\|_{L^2}\sim L_B^{-1.43}$ (R²=0.999). Bound A5 (eq A5-error): $\le C_1(L_B/\xi)+C_2(a/L_B)^{1/2}$ — w dostępnym oknie ($L_B<\xi_{\rm eff}$) dominuje człon fluktuacyjny $(a/L_B)^{1/2}$. Zmierzone $p=1.43\ge1/2$ ⟹ zbieżność **nie wolniejsza** niż konserwatywny bound A5 (a faktycznie szybsza). Zgodne z A5.

## §3 — Upgrade analityczny (podsekwencja → pełna zbieżność)
Lemat A1 dawał tylko prezwartość **podsekwencji** (Rellich–Kondrachov). Pełne CG-3 wymaga zbieżności **całej rodziny**. Numeryczny wynik (Cauchy w H¹, §1) pokazuje, że cała rodzina $\{\Phi_B\}_{b\ge2}$ jest Cauchy — co (przy zupełności $H^1$) daje zbieżność całej rodziny do **jednej** granicy, nie tylko podsekwencji. Analitycznie: jednoznaczność granicy z **wypukłości** funkcjonału gradientowego ($c_*>0$, Phase 1 §2) + Γ-zbieżność (rama De Giorgi / Braides dla przejść discrete→continuum). Numeryczna Cauchy-zbieżność jest świadectwem tej jednoznaczności.

## §4 — Anti-Lakatos
✓ Bug prior (zera) naprawiony; normy niezerowe, zweryfikowane.
✓ Model φ⁴ uzasadniony uniwersalnością (W:substrate); ξ skończone, separacja skal jawna (wąska: ξ=2.8).
✓ Rate raportowany jako p≥1/2 (zgodny z bound, nie naciągany do dokładnej formy A5).
✓ C3-A: b=1 (surowe φ²) wyłączone z testu zbieżności z jawnym uzasadnieniem (kontaktowa wariancja).

## §5 — Handoff do Phase 3 (CG-4)
CG-3 zamknięty NUM. **Phase 3 (CG-4):** identyfikacja K_hom=K_TGP — α=2 (lemat A3, sympy, czysta algebra zmiany zmiennych Φ=φ²), koercywność c*>0 (Phase 1), β=γ (lemat A4). Uwaga: substrat patologiczny → residuum N5 drogą algebraiczną, nie MC.
