---
title: "Phase 0 — pinowanie κ_E sektora σ_ab z substratu: czy C_σ (stała kinetyczna σ_ab) jest wyprowadzalne z H_Γ, pinując strumień Ṗ_b — i czy daje κ_E=5/6 (survive) czy ≠ (falsified)? Obiekt POPRAWNIE zidentyfikowany (σ_ab, nie skalar/g_eff)."
date: 2026-06-14
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-14
cycle: op-sigma-kinetic-Csigma-2026-06-14
authorization: "User 2026-06-14 (sesja #30): 'ok działaj' (cykl pinowania κ_E z substratu)"
object_id: "κ_E = C_σσ_0² / (wartość GR-lock); C_σ = stała kinetyczna σ_ab (eq:S-sigma); rem:sigma-params: PROBLEM OTWARTY"
parent_verdicts_LOCKED:
  - "op-disformal-hamiltonian-viability: kanał skalarno-disformalny BROKEN; sektor radiacyjny UNDERDETERMINED (po korekcie CL-7)"
  - "op-disformal-radiation-resolution: κ_E unpinned (det J=2)"
  - "PR-025: branże konforemne FALSIFIED (LOCKED)"
tags: [pre-registration, sigma-ab, kappa-E-pinning, C-sigma-open-problem, right-object, anti-Lakatos-LOCKED]
---

# Phase 0 — pinowanie κ_E (sektor σ_ab)

## §1 — Pytanie (i poprawna identyfikacja obiektu)
**Lekcja z 3 korekt: NAJPIERW fizyczny DOF niosący obserwablę.** Obserwabla = Ṗ_b (strata energii
binarnego pulsara). Nośnik GW w TGP = **σ_ab** (`ssec:tensor-substrate`; NIE skalar, NIE induced-TT —
disformalny tłumiony 18 rzędów, `rem:disformal-status`). Strumień σ_ab kontroluje $\kappa_E\propto
C_\sigma\sigma_0^2$ (niezależne od amplitudy $\xi_{\rm eff}$; det J≠0).

> **Czy $C_\sigma$ (stała kinetyczna σ_ab, eq:S-sigma) jest wyprowadzalne z H_Γ (coarse-graining
> korelacji kierunkowej $\sigma_{ab}=\langle\hat s\hat s\rangle^{\rm TF}$, def:sigma-ab) — a jeśli tak,
> czy daje $\kappa_E=5/6$ (przeżycie: σ_ab+skalar=GR) czy $\kappa_E\neq5/6$ (falsyfikacja)?**

Status wejściowy (`rem:sigma-params`): $C_\sigma$ „wyznaczalny w zasadzie z dynamiki substratu, ale
obecnie NIEZOBLICZONY" — **problem otwarty** (redukcja 3→2 parametrów TGP).

## §2 — Wejścia LOCKED
| Input | Forma | Źródło |
|---|---|---|
| Akcja σ_ab | $S_\sigma=\int\sqrt{-g}[\tfrac{C_\sigma}{2}(\partial\sigma)^2-\tfrac{m_\sigma^2}{2}\sigma^2+\tfrac{\xi_{\rm eff}}{\Phi_0^2}\sigma_{ab}\partial^a\Phi\partial^b\Phi]$ | eq:S-sigma |
| Amplituda (pinned) | $\xi_{\rm eff}=4\pi G_0\sigma_0\Phi_0$ (dopasowanie GR) | thm:amplitude-matching |
| Strumień (kontrola) | $\kappa_E\propto C_\sigma\sigma_0^2$; amplituda⊥strumień (det J≠0) | op-disformal-radiation-resolution T5 |
| σ_ab origin | $\sigma_{ab}=K_{ab}-\tfrac13\mathrm{Tr}K\,\delta$, $K_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle$ | def:sigma-ab |
| Skalar konforemny | radiuje $\tfrac16 P_{GR}$ nieuniknienie (nie da się zdrowo wyekranować — viability) | PR-025 + op-disformal-hamiltonian-viability |

## §3 — Falsyfikatory LOCKED
| ID | Test | Reguła |
|---|---|---|
| **F-CS-A** | Czy $C_\sigma$ wyprowadzalne z H_Γ (coarse-graining σ_ab)? | DERIVED / GAP (machinery niewystarczająca) |
| **F-CS-B** | Jeśli DERIVED: wartość $\kappa_E$ vs warunek przeżycia 5/6 | $\kappa_E=5/6\pm$tol → SURVIVE; $\neq$ → FALSIFIED |
| **F-CS-C** | Jeśli GAP: charakter swobody κ_E (miara survival) | survival measure-zero ⟹ UNDERDETERMINED-fine-tuned (lean falsified) |
| **F-CS-D** | Agregat | (A=GAP) ⟹ UNDERDETERMINED-fine-tuned; (A=DERIVED ∧ B≠5/6) ⟹ FALSIFIED; (∧ B=5/6) ⟹ SURVIVE |

## §4 — Pre-derywacja (oczekiwanie; rachunek nadrzędny)
Survival wymaga $\kappa_E=5/6$ EXACT (pojedynczy punkt). Naturalna wartość (σ_ab GR-matched także
w strumieniu, jak w GR gdzie amplituda i strumień zlockowane przez jedno G) = $\kappa_E\approx1$ →
total $7/6$ = **gałąź B PR-025 (2646σ FALSIFIED)**. Ale det J≠0 ⇒ κ_E NIE jest automatycznie zlockowane
(TGP ma 2 parametry tam gdzie GR ma 1) ⇒ κ_E zależy od underived $C_\sigma$. **Oczekiwanie: F-CS-A=GAP**
($C_\sigma$ = uznany problem otwarty), więc werdykt UNDERDETERMINED, ale survival = miara zero.

## §5 — Forbidden moves
1. Zakaz strojenia $C_\sigma$/$\sigma_0$ do κ_E=5/6 (to byłby tuning do danych — survival musi być WYPROWADZONE, nie dobrane).
2. Zakaz mylenia obiektu: κ_E to strumień σ_ab, NIE skalar, NIE g_eff, NIE induced-TT (lekcja 3 korekt).
3. Zakaz rewizji PR-025/parent LOCKED.
4. GAP deklarowany jawnie jeśli C_σ niewyprowadzalne (nie udawać derywacji — anti-Lakatos; lekcja sesji).
5. Budżet nowych stałych 0 (C_σ jest istniejącym parametrem, nie nowym).

## §6 — Anti-Lakatos
✓ Obiekt zidentyfikowany PRZED rachunkiem (§1) — bezpośrednia reakcja na 3 błędy „wrong object".
✓ Werdykt z flag. ✓ GAP dopuszczalny i oczekiwany (C_σ open). ✓ Zakaz strojenia survival. ✓ Liczby poprzedników LOCKED.
