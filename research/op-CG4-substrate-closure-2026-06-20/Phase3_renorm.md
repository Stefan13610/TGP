---
title: "Phase 3 — op-CG4-substrate-closure: renormalizacja operatora złożonego ∂ŝ∂ŝ (kanał TT spin-2 = σ_ab). WYNIK DECYDUJĄCY (sympy + num, value-blind): współczynnik p² (=C_σ) ma NIEZEROWĄ rozbieżność LINIOWĄ w odcięciu Λ — współczynnik kątowy ∫(1−μ²)²(4μ²−1)dμ=−16/35≠0; predykcja analit. wsp. liniowego −0.002895 vs num −0.002891 (zgodność 4 cyfry). ⟹ NIE istnieje scheme-independent continuum gołego operatora; C_σ jest UV-czułym WOLNYM PARAMETREM (wymaga warunku renormalizacyjnego = residuum bieguna σ_ab), NIEustalany przez bąbel. C-D = GAP (z DOWODEM natury), nie PARTIAL-survival. Tłumaczy pasmo lattice-MC #31 [0.04,11.1] (wartość skaluje się z odcięciem sieci)."
type: phase_result
status: PHASE3_COMPLETE (C-D = GAP z dowodem: C_σ UV-czuły wolny parametr; rozbieżność liniowa niezerowa)
phase: 3
cycle: op-CG4-substrate-closure-2026-06-20
created_date: 2026-06-20
script: "[[./Phase3_renorm.py]]"
output: "[[./Phase3_output.txt]]"
anti_lakatos_lock: PRESERVED
authorization: "User 2026-06-20: 'tak działaj' (Phase 3 renormalizacja)"
---

# Phase 3 — renormalizacja operatora złożonego (op-CG4-substrate-closure)

> **Werdykt fazy (value-blind, wyliczony):** współczynnik $p^2$ w kanale TT spin-2 ($\sigma_{ab}$),
> czyli **$C_\sigma$**, ma **niezerową rozbieżność LINIOWĄ** w odcięciu $\Lambda$. ⟹ **NIE istnieje
> scheme-independent continuum** gołego operatora $\partial\hat s\,\partial\hat s$. **$C_\sigma$ jest
> UV-czułym wolnym parametrem** (ustalany dopiero warunkiem renormalizacyjnym = residuum bieguna $\sigma_{ab}$),
> **nie** przez bąbel. **C-D = GAP** (z dowodem natury obstrukcji).

## §1 — Konstrukcja (free scalar, 3D euklidesowy)

$$\langle O_{ab}(p)\,O_{cd}(-p)\rangle_{\rm conn}=\int\!\frac{d^3k}{(2\pi)^3}\,
\frac{W_{ab}\,W_{cd}}{(k^2+m^2)\big((p-k)^2+m^2\big)},\quad
W_{ab}=\tfrac12(k_aq_b+k_bq_a),\ q=p-k$$

Projekcja TT ($p\parallel z$, kanał σ_ab = symetryczny bezśladowy poprzeczny):
$$|W^{\rm TT}|^2=\tfrac12(W_{xx}-W_{yy})^2+2W_{xy}^2
\;\overset{\rm sympy}{=}\;\tfrac{k_x^4}{2}+k_x^2k_y^2+\tfrac{k_y^4}{2}
\;\xrightarrow{\langle\cdot\rangle_\varphi}\;\tfrac12 k^4\sin^4\theta$$
**Kluczowe:** $|W^{\rm TT}|^2$ **nie zależy od $p$** (p wchodzi tylko przez drugi propagator).

## §2 — Współczynnik $p^2$ = $C_\sigma$ i jego rozbieżność

$$\partial_{p}^2\Big[\tfrac{1}{(p-k)^2+m^2}\Big]_{p=0}=\frac{4k_z^2-k^2-m^2}{(k^2+m^2)^3}
\ \Rightarrow\ C_\sigma=\tfrac12\!\int\!\frac{d^3k}{(2\pi)^3}\,
\frac{\tfrac12k^4\sin^4\theta\,(4k_z^2-k^2-m^2)}{(k^2+m^2)^4}$$

**Analiza wymiarowa:** $[O_{ab}]=3$ ($d{=}3$) ⟹ $[\langle OO\rangle(p)]=3$ ⟹ $[\text{coeff }p^2]=1$ ⟹ **rozbieżność liniowa** ($\sim\Lambda$).

**Współczynnik rozbieżności (kątowy, sympy):**
$$\int_{-1}^{1}(1-\mu^2)^2(4\mu^2-1)\,d\mu=-\frac{16}{35}\neq0\quad(\mu=\cos\theta)$$
**NIEZEROWY** ⟹ rozbieżność liniowa **nie znika** w kanale TT. Wsp. liniowy:
$\tfrac12\cdot\tfrac12\cdot(-\tfrac{16}{35})\cdot\tfrac{2\pi}{(2\pi)^3}=-0.002895$.

## §3 — Weryfikacja numeryczna $C_\sigma(\Lambda)$

| $\Lambda$ | 5 | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|---|
| $C_\sigma(\Lambda)$ | −0.0088 | −0.0228 | −0.0515 | −0.1093 | −0.2251 | −0.4566 |

Dopasowanie $C_\sigma(\Lambda)=a\,\Lambda+b$: **$a=-0.002891$**, $b=0.006$.
**$a_{\rm num}=-0.002891$ vs $a_{\rm analit}=-0.002895$ — zgodność do 4 cyfr znaczących.**
$|a|\gg0$ ⟹ **liniowo rozbieżne**, brak plateau ⟹ brak continuum scheme-independent.

## §4 — Interpretacja (decydująca)

1. **$C_\sigma$ nie jest predykcją.** Goły bąbel operatora złożonego daje wartość rosnącą liniowo z odcięciem; fizyczne $C_\sigma$ wymaga **warunku renormalizacyjnego** (kanonicznej normalizacji σ_ab = residuum bieguna), którego framework **nie dostarcza** niezależnie. ⟹ $C_\sigma$ = **wolny parametr UV** sektora radiacyjnego.
2. **To tłumaczy pasmo lattice-MC #31** $[0.04,11.1]$: sieć dostarcza konkretnego schematu (odcięcie $a=1$), więc daje wartość **skalującą się z odcięciem** — nie ma granicy continuum. Pasmo było sygnaturą **rozbieżności**, nie tylko statystyki.
3. **Spójne z `rem:sigma-params`** („$C_\sigma$ obecnie niezobliczony"): teraz **wiadomo DLACZEGO** — jest UV-czuły, nie tylko trudny. Konwertuje otwarte pytanie obliczeniowe w **zamknięty fakt strukturalny**.

## §5 — Mapowanie na kryteria

| Kryterium | Werdykt |
|---|---|
| **C-D** scheme-indep. $C_\sigma$ | ❌ **GAP** (z dowodem: rozbieżność liniowa niezerowa, wsp. −16/35; $C_\sigma$ = wolny parametr UV) |

**Reguła agregatu (Phase0_lock §4):** brama otwarta ∧ **C-D GAP** ⟹ CG-4 GAP dla pinowania $\kappa_E$;
sektor radiacyjny pozostaje **UNDERDETERMINED**, ale z **dowiedzioną** przyczyną (UV-czułość $C_\sigma$), nie brakiem rachunku.

## §6 — Anti-Lakatos
✓ Werdykt wyliczony (sympy + num, zgodność 4 cyfry), nie wybrany. ✓ $C_\sigma$ **NIE sfabrykowane** — dowiedziono, że jest rozbieżne/wolne (nie podano fałszywej liczby). ✓ Zero strojenia do 5/6. ✓ Wynik **negatywny** (GAP) zgłoszony wprost, nie obejście. ✓ Rachunek na free scalar (kontrolowany); interakcje mogą zmienić współczynnik, NIE fakt UV-czułości (wymiar operatora niezmieniony) — zgłoszone jako założenie. ✓ Rdzeń nie edytowany.

## §7 — Konsekwencja dla werdyktu cyklu
C-D przechodzi z **PARTIAL** (Phase 2: surowa magnituda scheme-dependent) → **GAP z dowodem** (Phase 3: rozbieżność
liniowa niezerowa ⟹ wolny parametr UV). Sektor radiacyjny: **UNDERDETERMINED — i to strukturalnie nieusuwalne
bąblem**. Survival ($\kappa_E=5/6$) możliwe tylko przez **ustawienie wolnego parametru**, nie predykcję
⟹ „fine-tuned" zyskuje **mechanizm**: wolna stała renormalizacyjna UV. Patrz [[Phase_FINAL_close.md]] (zaktualizowany).
