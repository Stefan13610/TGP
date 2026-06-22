---
title: "Phase 1 — op-CG4-substrate-closure: silnik analityczny. WYNIK (sympy, value-blind): (1) α_eff=s−1 potwierdzone — struktura (∇Φ)²/Φ generowana dla s≠1, ale α=2 NIE wynika z substratu (= postulat na gęstości, #32) ⟹ C-C(obecność) PASS, C-C(wartość) USTALENIE. (2) Klasyfikacja stabilności: M0(φ⁴,u>0) i M3(φ⁶,w>0) bounded-below; M1(−J(s_is_j)²) runaway lub redukcja do M0 (kontrola negatywna potwierdzona); M2(gradient-bond v2) kierunek płaski → frozen. (3) Typ przejścia: M0=ciągły (klasa Isinga), M3=ciągły+trikrytyczny tunable; M1/M2 wykluczone STRUKTURALNIE (forbidden #7 spełniony). (4) C_σ>0 DERIVED (bąbel 3D arctan); R-continuum = power-divergence D=3 (cubic) = bariera scheme-indep. REKOMENDACJA: substrat produkcyjny Phase 2 = M0, cross-check M3."
type: phase_result
status: PHASE1_COMPLETE (silnik analityczny; brama C-A/C-B rozstrzygnięta strukturalnie; C-C ustalone; C-D = Phase 2/3)
phase: 1
cycle: op-CG4-substrate-closure-2026-06-20
created_date: 2026-06-20
script: "[[./Phase1_engine.py]]"
output: "[[./Phase1_output.txt]]"
anti_lakatos_lock: PRESERVED
authorization: "User 2026-06-20: 'działaj z op-CG4-substrate-closure' (Phase 1 = silnik analityczny, bez MC)"
---

# Phase 1 — silnik analityczny (op-CG4-substrate-closure)

> **Werdykt fazy (value-blind, wyliczony sympy):** brama substratu **C-A ∧ C-B rozstrzygnięta strukturalnie** —
> niepatologiczny substrat **istnieje** (M0 φ⁴ klasa Isinga; M3 φ⁶ jako cross-check), a obie patologie
> (M1 runaway, M2 płaski kierunek) są wykluczone **argumentem strukturalnym** (forbidden #7 OK).
> **C-C** (struktura kinetyczna): $(\nabla\Phi)^2/\Phi$ **generowana**, ale α=2 **NIE** z substratu (postulat #32).
> **C-D** (scheme-independent $C_\sigma$): bariera = **R-continuum power-divergence D=3** — pozostaje do Phase 2/3.

## §1 — T1: α_eff = s−1 (reprodukcja CG-34 Phase 3, niezależnie)

Mikroskopowy kinetyk $Z(\phi)=Z_0\phi^{2s}$, zmiana zmiennej $\Phi=\phi^2$:
$$K_1(\Phi)=\frac{Z(\sqrt\Phi)}{4\Phi}=\frac{Z_0}{4}\,\Phi^{\,s-1},\qquad
\frac{K_1'}{2K_1}=\frac{s-1}{2\Phi}\ \Rightarrow\ \boxed{\alpha_{\rm eff}=s-1}$$

| $s$ | $Z(\phi)$ | $K_1(\Phi)$ | $\alpha_{\rm eff}$ |
|---|---|---|---|
| 0 | const (kinetyk standardowy) | $\propto1/\Phi$ | **−1** |
| 1 | $\propto\phi^2$ (substrat sek10) | const | **0** |
| 3 | $\propto\phi^6$ | $\propto\Phi^2$ | **2** |

> **USTALENIE (anti-Lakatos):** człon $(\nabla\Phi)^2/\Phi$ ma **niezerowy** współczynnik dla każdego $s\neq1$ —
> **struktura TGP jest generowana** (C-C obecność = PASS). ALE wartość $\alpha=2$ wymagałaby $s=3$ ($Z\propto\phi^6$),
> czego substrat sek10 ($Z\propto\phi^2$, $s=1$) nie daje. Zgodnie z **sesją #32**, $\alpha=2$ jest **postulatem
> aksjomatycznym na gęstości** (selekcja C1–C3), nie derywacją z substratu (substrat amplituda $\sqrt\Phi$ → α=½).
> **Nie fabrykujemy α=2.** C-C(wartość) = USTALENIE jawne; C-C(obecność) = PASS.

## §2 — T2/T3: klasyfikacja kandydatów (stabilność + typ przejścia)

| Model | bounded-below | typ przejścia Z₂ | werdykt |
|---|---|---|---|
| **M0** φ⁴: $r s^2+u s^4$ | ⟺ $u>0$ | **ciągły** ($b{=}u>0$, klasa Isinga 3D) | ✅ **C-A ∧ C-B PASS** |
| **M3** φ⁶: $+w s^6$ | ⟺ $w>0$ (dopuszcza $u<0$) | ciągły ($u>0$) / **trikrytyczny** ($u=0$) | ✅ **C-A ∧ C-B PASS** (tunable) |
| **M1** $-J(s_is_j)^2$ | ⟺ $u>Jz$ | dla $u<Jz$: runaway; dla $u>Jz$: ≡ M0 | ❌ patologia / brak korzyści (kontrola −) |
| **M2** gradient-bond v2 | człon ≥0, znika dla jednorodnego $s$ | krytyczność z on-site (=M0/M3), nie z bondu | ❌ kierunek płaski → frozen |

**Landau (sympy, T3):** $f(M)=\tfrac{a}{2}M^2+\tfrac{b}{4}M^4+\tfrac{c}{6}M^6$ → przejście **ciągłe** ⟺ $b>0$;
**trikrytyczne** ⟺ $b=0$; **1. rzędu** ⟺ $b<0$ ze skokiem $M^2=-3b/(4c)$ przy $a=3b^2/(16c)$ (zweryfikowane sympy).
M0 ($b=u>0$) leży po stronie ciągłej — to **kanoniczny model klasy Isinga 3D** z rygorystycznie ustalonym
punktem Wilsona–Fishera. **Patologia #31 była specyficzna dla bondu $-J(s_is_j)^2$ (M1), nie dla TGP jako takiego.**

> **Argument strukturalny wyboru (forbidden #7):** M0 wybrany NIE przez minimalizację dryfu, lecz jako
> **minimalny model o uniwersalności Z₂/Ising z udowodnionym ciągłym punktem krytycznym** — ten sam, który
> zwalidował silnik CG-34 Phase 1 ($\kappa_c\approx0.189$). M3 = cross-check robustności (trikrytyczność).

## §3 — T4: C_σ > 0 (bąbel 3D) + bariera R-continuum

Zamknięty bąbel 3D (sympy):
$$\Pi(p)=\frac{\arctan(p/2m)}{4\pi p}=\frac{1}{8\pi m}-\frac{p^2}{96\pi m^3}+\mathcal O(p^4)
\ \Rightarrow\ C_\sigma=-\,\mathrm{coeff}(p^2)=\frac{1}{96\pi m^3}>0$$
Znak $C_\sigma>0$ **DERIVED** (zgodny z op-Csigma-lattice-MC #31).

**Bariera (C-D):** dla operatora złożonego $O_{ab}=\partial_a\hat s\,\partial_b\hat s$ stopień rozbieżności
$$D=d+\#\text{deriv}-2\#\text{prop}=3+4-4=\boxed{3}\ (>0)\ \Rightarrow\ \text{UV power-divergent (cubic)}.$$
To **dokładnie** źródło scheme-dependence prefaktora $C_\sigma$ (R-continuum z #31). Znak i skaling są DERIVED;
**scheme-independent magnituda** wymaga additywnej subtrakcji UV + analizy mieszania operatorów — to **Phase 3**
(na niepatologicznym substracie M0, gdzie czysty punkt krytyczny czyni continuum dobrze postawionym).

## §4 — Mapowanie na kryteria Phase 0

| Kryterium | Stan po Phase 1 |
|---|---|
| **C-A** stabilność | ✅ rozstrzygnięte: M0/M3 bounded-below (analitycznie) |
| **C-B** czysty punkt krytyczny | ✅ rozstrzygnięte strukturalnie: M0 ciągły (Ising), M3 trikrytyczny — **brama otwarta** |
| **C-C** struktura kinetyczna | ☑ obecność PASS; wartość α=2 = USTALENIE (postulat #32, nie fabrykowane) |
| **C-D** scheme-indep. $C_\sigma$ | ⏳ **otwarte** — znak DERIVED; magnituda blokowana przez R-continuum (D=3) → Phase 2/3 |

**Brama agregatu (Phase0_lock §4):** warunek „substrat niepatologiczny = C-A ∧ C-B PASS" — **SPEŁNIONY** (M0/M3).
⟹ cykl **nie** kończy się w GAP-substrat; werdykt sektora zależy teraz **wyłącznie** od C-D (Phase 2/3).

## §5 — Anti-Lakatos
✓ α_eff=s−1 wyliczone samodzielnie (sympy), nie skopiowane z A3. ✓ α=2 **NIE sfabrykowane** (postulat #32, jawnie).
✓ Wybór M0 z **argumentu strukturalnego** (uniwersalność Ising), nie minimalizacji dryfu (forbidden #7).
✓ M1/M2 wykluczone analitycznie (runaway / płaski kierunek), patologia #31 zlokalizowana. ✓ $C_\sigma$ znak DERIVED;
magnituda jawnie otwarta (R-continuum D=3), nie naciągnięta. ✓ Rdzeń nie edytowany; budżet nowych stałych 0.

## §6 — Handoff do Phase 2 (wymaga „działaj Phase 2")
**Plan Phase 2 (numeryczny, ciężki):**
1. Walidacja silnika MC na M0 (faza Z₂, $\langle\phi^2\rangle$, $\xi$, pik $\chi$ przy $\kappa_c\approx0.189$ — reuse CG-34) — **forbidden #4**.
2. Pomiar $C_\sigma$ z **structure factor operatora złożonego** $\langle O_{ab}(k)O_{cd}(-k)\rangle$ na M0, skan $L\in\{16,24,32,48\}$.
3. Additywna subtrakcja UV (R-continuum) + test mieszania operatorów → kandydat na scheme-independent $C_\sigma$.
4. Cross-check M3 (trikrytyczny) — czy $C_\sigma$ stabilne między klasami uniwersalności.
> Phase 2 to realna kampania MC (godziny–dni CPU) + renormalizacja Phase 3. **Tu się rozstrzyga twardy werdykt $\kappa_E$.**
