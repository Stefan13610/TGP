---
title: "Phase 1 — pinowanie κ_E: F-CS-A = GAP (C_σ niewyprowadzone — uznany problem otwarty rdzenia), F-CS-C = survival MIARA ZERO (κ_E=5/6). Werdykt: UNDERDETERMINED-fine-tuned (strukturalny lean ku FALSIFIED; naturalna κ_E≈1 → gałąź B PR-025). Brama decydująca precyzyjnie wskazana: C_σ z H_Γ."
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-sigma-kinetic-Csigma-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14: 'ok działaj'"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "EXACT: det J(amp,flux|C_σ,σ_0)=−ξ/C_σ≠0 (amplituda⊥strumień); survival κ_E=5/6; κ_E=1→7/6=gałąź B PR-025."
flag_F_CS_A: "GAP (C_σ niewyprowadzone z H_Γ w tym cyklu — uuznany problem otwarty rem:sigma-params; NIE fabrykowane)"
flag_F_CS_C: "MEASURE-ZERO (survival ⟺ κ_E=5/6 EXACT, pojedynczy punkt)"
flag_F_CS_D: "UNDERDETERMINED-fine-tuned (lean FALSIFIED)"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — pinowanie κ_E (sektor σ_ab)

## §0 — Verdict at a glance
| Element | Wynik |
|---|---|
| **Obiekt** (poprawny) | σ_ab — nośnik GW; κ_E∝C_σσ_0² (strumień), ⊥ ξ_eff (amplituda) |
| **F-CS-A** (C_σ z H_Γ) | **GAP** — $C_\sigma$ niewyprowadzone (uznany problem otwarty rdzenia; NIE fabrykuję) |
| **amplituda⊥strumień** | EXACT: det J=−ξ/C_σ≠0 — κ_E nie zlockowane do amplitudy (TGP: 2 param, GR: 1) |
| **F-CS-C** (survival) | **MIARA ZERO** — przeżycie ⟺ κ_E=5/6 EXACT (skalar 1/6 + σ_ab κ_E = GR) |
| **naturalna κ_E≈1** | total 7/6 = **gałąź B PR-025 (2646σ FALSIFIED)** |
| **F-CS-D** (agregat) | **UNDERDETERMINED-fine-tuned** (lean FALSIFIED) |

## §1 — Co wykazano twardo (EXACT)
1. **κ_E swobodne, NIE zlockowane do amplitudy:** σ_ab daje amplitudę = GR (ξ_eff matched) i NIEZALEŻNIE
   strumień ∝ C_σσ_0² (det J=−ξ/C_σ≠0). To źródło underdeterminacji: GR ma jeden parametr (G) lockujący
   amplitudę↔strumień; TGP-σ ma dwa (ξ_eff, C_σ) ⇒ strumień jest wolny mimo dopasowanej amplitudy.
2. **Warunek przeżycia jest punktem:** Ṗ_b = κ_E·P_GR + (1/6)·P_GR (σ_ab + nieunikniony skalar konforemny;
   skalar 1/6 NIE da się zdrowo wyekranować — viability, op-disformal-hamiltonian-viability). Obserwacja=P_GR
   ⟹ **κ_E = 5/6** dokładnie. Pojedynczy punkt w continuum dopuszczalnych C_σ ⟹ **miara zero**.
3. **Naturalna wartość → falsyfikacja:** jeśli σ_ab jest GR-podobne także w strumieniu (κ_E≈1, oczekiwane
   gdyby wspólne pochodzenie substratowe ⟨ŝŝ⟩ lockowało jak GR), total = 7/6 = **dokładnie gałąź B PR-025
   (2646σ)**. Falsyfikacja jest „naturalnym" wynikiem; przeżycie wymaga niegenerycznego κ_E=5/6.

## §2 — F-CS-A = GAP (uczciwie, NIE fabrykuję C_σ)
Wyprowadzenie $C_\sigma$ wymaga coarse-grainingu H_Γ dla korelacji kierunkowej
$\sigma_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle^{\rm TF}$ (def:sigma-ab) z członu wiązania kierunkowego
$J\sum A_{ij}\hat s_i^2\hat s_j^2(\hat s_j^2-\hat s_i^2)^2$ — analogicznie jak skalarne $K_{\rm geo}=\varphi^4$
z H_Γ→S[Φ], ale dla tensorowego sektora bezśladowego. To jest **uznany problem otwarty** (`rem:sigma-params`:
„obecnie niezobliczony"; `rem:param-counting`: redukcja 3→2 parametrów „problem otwarty"). **W tym cyklu
NIE wyprowadzam $C_\sigma$** — fabrykowanie liczby naruszyłoby anti-Lakatos (i powtórzyłoby błąd sesji).
Deklaruję GAP jawnie.

**Co byłoby potrzebne (specyfikacja dla cyklu domykającego):** pełny coarse-graining tensorowego sektora
substratu (dod. app:substrat machinery) → stała sztywności $C_\sigma$ → porównanie $\kappa_E=C_\sigma\sigma_0^2/(\text{GR-lock})$
z 5/6. To program wielosesyjny (klasa op-gamma-RG-running dla sektora tensorowego).

## §3 — Werdykt F-CS-D
> **UNDERDETERMINED-fine-tuned.** Sektor radiacyjny pozostaje UNDERDETERMINED (κ_E niepinowane — C_σ
> open), ALE: (i) survival to **miara zero** (κ_E=5/6 EXACT), (ii) wartość naturalna κ_E≈1 daje **gałąź B
> PR-025 (falsyfikacja)**. Strukturalny lean: **ku FALSIFIED**. Nie jest to hard-verdykt — bramą jest
> obliczenie C_σ (problem otwarty rdzenia), które rozstrzygnie definitywnie.

**Postęp względem op-disformal-radiation-resolution (UNDERDETERMINED):** tam „κ_E swobodne, status
nierozstrzygnięty". Tu **zaostrzone:** survival = miara zero, naturalna wartość = falsyfikacja, brama
decydująca (C_σ z H_Γ) precyzyjnie wskazana i sklasyfikowana jako otwarty problem. Przestrzeń przeżycia
skurczona z „nieokreślona" do „pojedynczy fine-tuned punkt".

## §4 — Anti-Lakatos
✓ Obiekt zidentyfikowany przed rachunkiem (σ_ab, nie skalar — lekcja 3 korekt). ✓ EXACT (det J, survival).
✓ **GAP zadeklarowany jawnie — NIE sfabrykowano C_σ** (kluczowe wobec błędów sesji). ✓ Survival NIE strojone
(κ_E=5/6 raportowane jako wymóg, nie dobrane). ✓ Liczby poprzedników LOCKED. ✓ Lean ku falsyfikacji
oznaczony jako strukturalny (nie-twardy), bo det J≠0 ⇒ κ_E formalnie wolne. ✓ Budżet stałych 0.

**Następny krok (opcjonalny, user): cykl coarse-grainingu $C_\sigma$ z H_Γ** (wielosesyjny; jedyna droga
od UNDERDETERMINED-fine-tuned do definitywnego werdyktu sektora). Bez niego: sektor pozostaje formalnie
nierozstrzygnięty, ale z przeżyciem o mierze zero.
