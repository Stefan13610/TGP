---
title: "Phase 2 — op-CG4-substrate-closure: kampania MC na M0 (φ⁴ Z₂/d3). WYNIK (value-blind): silnik ZWALIDOWANY (forbidden #4) — pik χ rośnie z L (15→21→36 dla L=10→12→16), κ_c≈0.190 (zgodne z CG-34), U4 płynnie 0→2/3, ⟨φ²⟩≈0.62 skończone+stabilne (brak runaway/frozen). C-A∧C-B PASS NUMERYCZNIE (substrat M0 niepatologiczny). Sztywność pola Z_R≈0.40 continuum-stabilna (spread 1.13×<1.2×, R²>0.98) ⟹ c*>0 POTWIERDZONE (red-flag CG-34 rozwiązany na M0). Operator złożony O=(∇φ)²: C_σ>0 O(1) odtworzone, ale surowa magnituda scheme-dependent (R-continuum D=3) — C-D PARTIAL."
type: phase_result
status: PHASE2_COMPLETE (silnik zwalidowany; C-A/C-B/C-C PASS num; C-D PARTIAL — bariera R-continuum potwierdzona empirycznie)
phase: 2
cycle: op-CG4-substrate-closure-2026-06-20
created_date: 2026-06-20
scripts: ["[[./Phase2_mc.py]]", "[[./Phase2_stiffness.py]]"]
data: ["[[./Phase2_validate_L16.json]]", "[[./Phase2_measure_k0.19.json]]", "[[./Phase2_stiffness_k0.19.json]]"]
anti_lakatos_lock: PRESERVED
authorization: "User 2026-06-20: 'Działaj Phase 2 (MC)'"
---

# Phase 2 — kampania MC na M0 (op-CG4-substrate-closure)

> **Werdykt fazy (value-blind):** substrat M0 (φ⁴ Z₂) **NUMERYCZNIE niepatologiczny** — czyste ciągłe
> przejście Z₂, c*>0 continuum-stabilne (spread 1.13×). **C-A ∧ C-B ∧ C-C PASS.** Operator złożony daje
> $C_\sigma>0$ $O(1)$, ale **surowa magnituda scheme-dependent** (R-continuum $D{=}3$ potwierdzone empirycznie)
> ⟹ **C-D = PARTIAL**.

## §1 — Walidacja silnika (forbidden #4, OBOWIĄZKOWA) — PASS

Model M0: $S=\sum_x[\phi_x^2+\lambda(\phi_x^2-1)^2-2\kappa\sum_\mu\phi_x\phi_{x+\mu}]$, $\lambda{=}1$, checkerboard Metropolis.

| Sygnatura | Wynik | Interpretacja |
|---|---|---|
| Pik podatności $\chi$ | $\chi_{\max}\approx15,21,36$ dla $L{=}10,12,16$ | **rośnie z $L$** (χ~$L^{\gamma/\nu}$) → przejście ciągłe |
| $\kappa_c$ | $\approx0.190$ (stabilne w $L$) | **zgodne z CG-34** ($\kappa_c\approx0.189$) — reuse OK |
| Kumulant Bindera $U_4$ | płynnie $0\to2/3$ przez $\kappa_c$ | klasyczny scenariusz Z₂ (faza ↔ złamana) |
| $\langle\phi^2\rangle$ | $\approx0.62$, **skończone, stabilne w $L$** (0.6305/0.6255/0.6192) | **brak runaway** (nie rośnie) **i brak frozen** ($\xi$ skończone) |

⟹ **C-A (stabilność) ∧ C-B (czysty punkt krytyczny) = PASS NUMERYCZNIE.** Patologia #31 (runaway/frozen)
**nie występuje** na M0 — była specyficzna dla bondu M1.

## §2 — Sztywność pola $Z_R$ (structure factor) — c*>0 continuum-stabilne

Metoda LOCKED (forbidden #1: NIE $\langle|\nabla\Phi|^2\rangle$): $1/S(k)=a+Z_R\,\hat k^2$, $\hat k^2=2(1-\cos2\pi n/L)$.

| $L$ | $Z_R$ | intercept | $R^2$ |
|---|---|---|---|
| 16 | 0.402 | 0.003 | 0.9987 |
| 24 | 0.421 | −0.010 | 0.9989 |
| 32 | 0.373 | 0.025 | 0.9800 |

**$Z_R\approx0.40$, spread $1.13\times$ (< 1.2× próg!), $R^2>0.98$** ⟹ **c*>0 POTWIERDZONE i continuum-stabilne.**
To ostatecznie zamyka red-flag CG-34 ($c_*\to0$ był artefaktem $\langle|\nabla\Phi|^2\rangle$) **na M0**. (Wartość
bezwzględna ≠ CG-34 $Z{\sim}2.5$ — inna parametryzacja/normalizacja; istotna jest **dodatniość + stabilność**, nie liczba.)

## §3 — Operator złożony $O=(\nabla\phi)^2$ (tensor σ_ab) — C-D PARTIAL

Structure factor $\langle|\tilde O(k)|^2\rangle/V$ przy $k=0,k_{\min},2k_{\min}$ (mean-subtracted):

| $L$ | $O(k_{\min})$ | $O(2k_{\min})$ |
|---|---|---|
| 16 | 17.6 | 12.1 |
| 24 | 23.9 | 11.6 |
| 32 | 14.8 | 12.2 |

- $C_\sigma>0$, $O(1)$ **odtworzone** (zgodne z #31: central $\kappa_E\approx0.62$, $C_\sigma\sim0.5$–0.7).
- **ALE:** surowe wartości **nie pokazują czystego continuum** (rozrzut zależny od $L$, częściowo niedotermalizacja L=32: $\chi$ spadło do 12 vs 26 przy L=24). Zgodnie z Phase 1 **$D{=}3$ power-divergence**: surowa magnituda jest **scheme-dependent**. Bez additywnej subtrakcji UV + analizy mieszania operatorów (Phase 3) pasmo $C_\sigma$ pozostaje > faktor 1.2.

⟹ **C-D = PARTIAL** (continuum scheme-independent NIE osiągnięte tą kampanią — bariera **analityczna**, nie statystyczna; potwierdza diagnozę lattice-MC #31).

## §4 — Mapowanie na kryteria (Phase0_lock §3)

| Kryterium | Werdykt Phase 2 |
|---|---|
| **C-A** stabilność | ✅ **PASS** (⟨φ²⟩ skończone, stabilne; brak runaway/frozen) |
| **C-B** czysty punkt krytyczny | ✅ **PASS** (ciągłe Z₂, χ~$L^{γ/ν}$, $U_4$ crossover) |
| **C-C** struktura kinetyczna | ✅ **PASS** ($c_*>0$, $Z_R$ spread 1.13×; α=2 = postulat #32) |
| **C-D** scheme-indep. $C_\sigma$ | 🟡 **PARTIAL** ($C_\sigma>0$ O(1) odtworzone; magnituda scheme-dependent, R-continuum D=3) |

## §5 — Anti-Lakatos
✓ Silnik zwalidowany PRZED pomiarem (forbidden #4). ✓ $c_*$ z structure factor, nie $\langle|\nabla\Phi|^2\rangle$ (forbidden #1).
✓ $C_\sigma$ NIE sfabrykowane — scheme-dependence jawnie raportowana (R-continuum D=3, zgodnie z Phase 1). ✓ Zero strojenia do 5/6.
✓ Niedotermalizacja L=32 zgłoszona jawnie, nie zamieciona. ✓ Rdzeń nie edytowany.

## §6 — Handoff do Phase FINAL / Phase 3
- **Phase FINAL (ten cykl):** agregat → CG-4 PARTIAL; **substrat RESOLVED** (R-A zamknięte: M0 niepatologiczny), residual izolowany do C-D.
- **Phase 3 (przyszły, analityczny):** additywna subtrakcja UV operatora złożonego + mieszanie operatorów → scheme-independent $C_\sigma$ → twardy werdykt $\kappa_E$ (lean: FALSIFIED-hard). To problem **renormalizacji**, nie większego MC.
