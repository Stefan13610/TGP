---
title: "SCOPING — op-Csigma-coarse-graining: wyprowadzenie stałej kinetycznej C_σ (i normalizacji σ_0) sektora tensorowego σ_ab z hamiltonianu substratu H_Γ przez coarse-graining, by przypiąć κ_E i rozstrzygnąć sektor radiacyjny (FALSIFIED / SURVIVE / GAP). Klasa op-gamma-RG-running, wielosesyjny."
date: 2026-06-14
type: meta-scoping
status: "PRE-PHASE-0 NOTE (wymaga własnego Phase 0 + 'działaj'; zero werdyktów)"
origin: "User 2026-06-14 (sesja #30→#31): 'rozpisz ten cykl' — po op-sigma-kinetic-Csigma (UNDERDETERMINED-fine-tuned; brama decydująca = C_σ z H_Γ)"
parent_cycle: "[[../research/op-sigma-kinetic-Csigma-2026-06-14/]] (survival = κ_E=5/6 miara zero; C_σ = problem otwarty)"
target_open_problem: "rem:sigma-params + rem:param-counting (sek08): C_σ 'wyznaczalny w zasadzie z H_Γ, obecnie niezobliczony'; redukcja 3→2 parametrów TGP"
anti_lakatos_note: "Pre-derywacja = oczekiwanie (κ_E≠5/6 generycznie → FALSIFIED), rachunek nadrzędny. GAP (C_σ nieobliczalny obecną machinerią) = pełnoprawny wynik. Zakaz strojenia C_σ/σ_0 do κ_E=5/6."
tags: [scoping, C-sigma, coarse-graining, tensor-sector, kappa-E-pinning, multi-session, anti-Lakatos]
---

# SCOPING — wyprowadzenie C_σ z substratu (domknięcie sektora radiacyjnego)

## §1 — Pytanie i obiekt (precyzyjnie)

op-sigma-kinetic-Csigma ustalił: sektor radiacyjny = UNDERDETERMINED-fine-tuned; przeżycie ⟺
$\kappa_E=5/6$ (miara zero); wartość naturalna $\kappa_E\approx1$ → falsyfikacja. **Jedyna brama do
definitywnego werdyktu: wyprowadzić $C_\sigma$ (i $\sigma_0$) z $H_\Gamma$.**

> **Czy stała kinetyczna $C_\sigma$ w $S_\sigma=\int\sqrt{-g}[\tfrac{C_\sigma}{2}(\partial\sigma)^2-\dots]$
> (eq:S-sigma) jest wyprowadzalna przez coarse-graining korelacji kierunkowej
> $\sigma_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle^{\rm TF}$ (def:sigma-ab) z $H_\Gamma$ —
> a jeśli tak, czy złożona w $\kappa_E=C_\sigma\sigma_0^2/(\text{GR-lock})$ daje $5/6$ (SURVIVE)
> czy $\neq5/6$ (FALSIFIED)?**

**Obiekt LOCKED (lekcja 3 korekt sesji):** $C_\sigma$ = sztywność kinetyczna tensorowego DOF σ_ab
(nośnik GW), NIE skalar, NIE g_eff, NIE induced-TT. Sprzężenie $\xi_{\rm eff}$ JUŻ przypięte
(thm:amplitude-matching); $m_\sigma\to0$ przyjęte; **tylko $C_\sigma$/$\sigma_0$ otwarte**.

## §2 — Trasa wyprowadzenia (fazy; klasa op-gamma-RG-running)

**Phase 1 — Formalne H_Γ dla σ_ab (parallel op-gamma-RG-running Phase 1).** Zidentyfikować człon
$H_\Gamma$ kontrolujący korelacje KIERUNKOWE: $J\sum_{\langle ij\rangle}A_{ij}\hat s_i^2\hat s_j^2(\hat s_j^2-\hat s_i^2)^2$
(eq:B-H, l.6700) + GL-bond v2 ($K_{ij}=J(\phi_i\phi_j)^2$). Operator: $\sigma_{ab}=K_{ab}-\tfrac13\mathrm{Tr}K\,\delta$,
$K_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle$. Ustalić, że σ_ab jest KOMPOZYTEM (jak Φ=⟨ŝ²⟩), więc
jego kinetyka pochodzi z dynamiki korelacji, nie z fundamentalnego pola. Inwentarz parametrów {J, a_sub, μ, m_0², λ_0}.

**Phase 2 — Coarse-graining / gradient expansion → C_σ (RDZEŃ, najtrudniejszy).** Blokowe uśrednianie
$H_\Gamma$ → efektywne działanie σ_ab(x); ekstrakcja współczynnika przy $(\partial_\mu\sigma_{ab})^2$ = $C_\sigma$
w funkcji {J, a_sub, ...}. Analogia: skalarne $K_{\rm geo}\varphi^4$ z metryki konforemnej + block-averaging.
**Wymóg Lorentza:** kinetyka musi być Lorentz-inwariantna (σ_ab propaguje z $c_0$, rem:cGW-tensor) — to
ogranicza formę (czas-przestrzeń symetrycznie, jak emergentny $\Box$ skalara).

**Phase 3 — Normalizacja σ_0 + złożenie κ_E.** $\sigma_0$ z relacji $h^{TT}=2\sigma_{ab}/\sigma_0$
(eq:h-TT-from-sigma). Rozstrzygnąć, czy $C_\sigma$ i $\sigma_0$ to jeden parametr (rem:sigma-params:
„C_σ lub równoważnie σ_0") czy dwa. Złożyć $\kappa_E=C_\sigma\sigma_0^2/(\text{wartość GR-lock})$, gdzie
GR-lock value ~ sztywność grawitonu $\sim c^3/(16\pi G_0)$, $G_0=q^2/(4\pi\Phi_0^2 K_1)$ (sektor skalarny).

**Phase FINAL — κ_E vs 5/6 → werdykt.** Porównanie value-blind (próg pre-LOCKED w Phase 0):
$\kappa_E=5/6\pm$tol → SURVIVE; $\neq$ → FALSIFIED; nieobliczalne → GAP.

## §3 — Falsyfikatory (szkic do LOCK w Phase 0)

| ID | Test | Reguła |
|---|---|---|
| F-CG-A | Czy σ_ab-kinetyka emerguje Lorentz-inwariantnie z H_Γ (C_σ>0, propagacja c_0)? | TAK → przejdź; NIE (brak emergentnej kinetyki) → GAP strukturalny |
| F-CG-B | Forma C_σ({J,a_sub,...}) z gradient expansion | DERIVED (jawna formuła) / PARTIAL (skaling, prefaktor O(1) otwarty) / GAP |
| F-CG-C | σ_0 i złożenie κ_E = C_σσ_0²/(GR-lock) | wartość κ_E lub przedział |
| F-CG-D | κ_E vs survival 5/6 (próg LOCKED Phase 0) | =5/6±tol SURVIVE; ≠ FALSIFIED; tylko skaling → PARTIAL-lean |
| F-CG-E | Agregat | DERIVED+κ_E≠5/6 → FALSIFIED; +5/6 → SURVIVE; GAP/PARTIAL → sektor pozostaje UNDERDETERMINED-fine-tuned (status węższy) |

## §4 — Reuse machinery (LOCKED)
- op-gamma-RG-running-derivation (skalarne H_Γ→S[Φ], Wilsonian, K_geo) — wzorzec metodologiczny.
- dodatekQ_coarse_graining_formal.tex; GL-bond v2 axiom; def:sigma-ab; eq:S-sigma; thm:amplitude-matching.
- op-sigma-kinetic-Csigma (survival κ_E=5/6, det J — strukturalne wejście).

## §5 — Forbidden moves (szkic)
1. **Zakaz strojenia C_σ/σ_0 do κ_E=5/6** (survival musi WYNIKNĄĆ z coarse-grainingu, nie być dobrane — value-blind).
2. Zakaz mylenia obiektu (C_σ = sztywność σ_ab; nie skalar/g_eff/induced-TT — lekcja sesji).
3. Zakaz rewizji PR-025/viability/parent LOCKED.
4. **GAP/PARTIAL deklarowane jawnie** jeśli coarse-graining tensorowy przekracza obecną machinerię — zakaz fabrykowania prefaktora (anti-Lakatos, lekcja #30).
5. Zakaz nowych stałych (C_σ istnieje; {J,a_sub,μ,...} to parametry substratu).
6. Wymóg Lorentza (c_0) na kinetykę σ_ab — nie naruszać.

## §6 — Risk register
- **R-tensor-hard (HIGH):** coarse-graining TENSOROWEGO (bezśladowego) sektora jest trudniejszy niż skalarnego; możliwy PARTIAL (skaling bez prefaktora) lub GAP. To realne — `rem:sigma-params` flaguje jako otwarte.
- **R-Lorentz-emergence (MED):** czy statyczny człon kierunkowy H_Γ w ogóle generuje Lorentz-inwariantną kinetykę $(\partial_t\sigma)^2-(\nabla\sigma)^2$ (a nie tylko przestrzenną sztywność)? Skalar dostał to z GL-bond; dla σ_ab do sprawdzenia (F-CG-A).
- **R-conspiracy (MED, INFORMATIONAL):** κ_E=5/6 wymaga, by tensorowe C_σσ_0² (z ⟨ŝŝ⟩) i skalarne G_0 (z ⟨ŝ²⟩) były powiązane przez wspólne pochodzenie ŝ DOKŁADNIE tak, by dać 5/6. Pre-flag: to spisek dwóch sektorów; generycznie nie. Ale wspólne H_Γ czyni relację możliwą — uczciwy test.
- **R-σ_0-vs-C_σ (MED):** jeśli to jeden parametr (rem:sigma-params), κ_E może zależeć od jednej kombinacji — uprościć w Phase 3.

## §7 — Oczekiwane wyniki (INFORMATIONAL; rachunek nadrzędny)
Najbardziej prawdopodobne: **PARTIAL/GAP** (skaling C_σ ~ J·f(a_sub) bez prefaktora O(1)) — wtedy
κ_E znany co do rzędu, ale nie co do 5/6 ⇒ sektor pozostaje UNDERDETERMINED-fine-tuned, ALE z węższym
przedziałem (jeśli skaling wyklucza 5/6 o rzędy → de facto FALSIFIED). Drugi: **DERIVED → κ_E≠5/6 →
FALSIFIED** (definitywne domknięcie, najmocniejszy wynik). Trzeci (mało prawdopodobny): κ_E=5/6 SURVIVE.

## §8 — Rejestracja
**`op-Csigma-coarse-graining`** — REGISTERED (NOT activated; własny Phase 0 + „działaj"). Wielosesyjny
(klasa op-gamma-RG-running, ~kilka-kilkanaście faz/sesji). Rozstrzyga ostatnią otwartą bramę sektora
grawitacyjnego radiacyjnego. [[../research/op-Csigma-coarse-graining-2026-06-14/]] (README).
