---
title: "Phase 1 — op-Csigma-coarse-graining: formalne H_Γ dla σ_ab. WYNIK: σ_ab = KOMPOZYT bilinowy ⟨ŝ_iŝ_j⟩ (rzut anizotropowy tego samego H_Γ co Φ); kinetyka C_σ pochodzi z dynamiki KORELATORA, nie z pola fundamentalnego. Inwentarz {J,a_sub,μ,m_0²,λ_0}. C_σ NIE wyprowadzone (Phase 2). F-CG-A = PARTIAL-pending; F-CG-B/C/D = OPEN."
type: phase_result
status: PHASE1_COMPLETE (SETUP)
phase: 1
cycle: op-Csigma-coarse-graining-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): autoryzacja fazy 1"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "PASS: σ_ab sym/bezśladowy/5-komp; Z₂-parzysty; [C_σ][σ]²=E²; C_σ>0 z J>0; skaling C_σ~c·J·a_sub^p (c,p=Phase 2)"
flag_F_CG_A: "PARTIAL-pending (c_0 LOCKED rem:cGW-tensor; emergencja członu czasowego (∂_tσ)² ze statycznego H_Γ = test Phase 2)"
flag_F_CG_B: "OPEN → Phase 2 (tu: szkielet skalowania + znak C_σ>0; BEZ prefaktora — anti-Lakatos)"
flag_F_CG_C: "OPEN → Phase 3"
flag_F_CG_D: "OPEN → Phase FINAL (próg 5/6 LOCKED Phase 0)"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — formalne $H_\Gamma$ dla $\sigma_{ab}$ (SETUP)

## §0 — Werdykt fazy w skrócie
| Element | Wynik |
|---|---|
| **Obiekt** | $\sigma_{ab}$ = **kompozyt bilinowy** $\langle\hat s_i\hat s_j\rangle$ — rzut ANIZOTROPOWY tego samego $H_\Gamma$, którego rzut IZOTROPOWY to $\Phi=\langle\hat s^2\rangle$ |
| **Źródło kierunkowych korelacji** | człon wiązania sąsiadów: $-J\sum_{\langle ij\rangle}\hat s_i\hat s_j$ (v1) / GL-bond $J\sum A_{ij}\hat s_i^2\hat s_j^2(\hat s_j^2-\hat s_i^2)^2$ (v2) |
| **Natura kinetyki** | $C_\sigma$ pochodzi z **dynamiki korelatora** (gradient korelatora między blokami), NIE z pola fundamentalnego — jak $K(\Phi)$ dla skalara |
| **Strukturalnie PASS (sympy)** | sym./bezśladowy/5-komp (spin-2); $\mathbb{Z}_2$-parzysty; $[C_\sigma][\sigma]^2=E^2$; $C_\sigma>0\Leftarrow J>0$ |
| **F-CG-A** (Lorentz) | **PARTIAL-pending** — $c_0$ LOCKED (rem:cGW-tensor), ale emergencja $(\partial_t\sigma)^2$ ze STATYCZNEGO $H_\Gamma$ = nietrywialny test Phase 2 |
| **F-CG-B** ($C_\sigma$) | **OPEN** — szkielet $C_\sigma\sim c\,J\,a_{\rm sub}^p$; $(c,p)$ = Phase 2; **prefaktor NIE fabrykowany** |
| **Inwentarz** | $\{J,\,a_{\rm sub},\,\mu,\,m_0^2,\,\lambda_0\}$ |

## §1 — Identyfikacja członu $H_\Gamma$ kontrolującego korelacje KIERUNKOWE

Substrat (eq:B-H): pole $\hat s_i\in\mathbb R$ na sieci $\Gamma$, symetria $\mathbb Z_2$:
$$H_\Gamma=\sum_i\Big[\tfrac{m_0^2}{2}\hat s_i^2+\tfrac{\lambda_0}{4}\hat s_i^4\Big]-J\sum_{\langle ij\rangle}\hat s_i\hat s_j\qquad(\text{v1}),$$
forma kanoniczna v2 z wiązaniem gradientowym GL: $J\sum_{\langle ij\rangle}A_{ij}\,\hat s_i^2\hat s_j^2(\hat s_j^2-\hat s_i^2)^2$.

**Rozdział izotrop./anizotrop. (rdzeń, fbox ssec:tensor-substrate):**
- **Człon on-site** $(\tfrac{m_0^2}{2}\hat s^2+\tfrac{\lambda_0}{4}\hat s^4)$ — kontroluje **amplitudę** $\Phi=\langle\hat s^2\rangle$ (VEV $v^2=|m_0^2|/\lambda_0$, eq:B-v2-vacuum). NIE wnosi informacji kierunkowej.
- **Człon wiązania** $-J\sum_{\langle ij\rangle}\hat s_i\hat s_j$ — to **jedyny** człon, który koduje korelacje **między sąsiadami** $\langle\hat s_i\hat s_{i+\hat a_b}\rangle$, a więc **kierunkowe**. To on źródłuje $K_{ab}$ (eq:K-ab).

> **Konkluzja §1:** $K_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle$ (a więc i $\sigma_{ab}$) jest sterowane **wyłącznie przez sektor wiązaniowy** $H_\Gamma$ (J-term). On-site (m_0², λ_0) wchodzi jedynie pośrednio przez normalizację $\langle\hat s^2\rangle$. To rozdziela rolę parametrów: **$J$ = sztywność kierunkowa (→ $C_\sigma$); $(m_0^2,\lambda_0)$ = amplituda (→ $\sigma_0$ przez VEV).** Ten rozdział jest strukturalnym wejściem do Phase 3 (C_σ vs σ_0).

## §2 — $\sigma_{ab}$ jest KOMPOZYTEM (wynik główny Phase 1)

$\sigma_{ab}$ NIE jest fundamentalnym polem dodanym ad hoc. Z def:sigma-ab:
$$\sigma_{ab}(x)=K_{ab}(x)-\tfrac13\mathrm{Tr}K\,\delta_{ab},\qquad K_{ab}(x)=\tfrac{1}{|\mathcal B|}\sum_{i\in\mathcal B(x)}\langle\hat s_i\hat s_{i+\hat a_b}\rangle.$$
Jest to **dwupunktowy korelator kierunkowy** ŝ, uśredniony blokowo — **operator złożony stopnia 2** w $\hat s$.

**Paralela z $\Phi$ (zweryfikowana sympy, krok 2):**

| | $\Phi$ | $\sigma_{ab}$ |
|---|---|---|
| definicja | $\langle\hat s_i^2\rangle$ (jednopunktowy ŝ²) | $\langle\hat s_i\hat s_{i+\hat a_b}\rangle$ (dwupunktowy, kierunkowy) − ślad |
| stopień w ŝ | 2 (bilinowy) | 2 (bilinowy) |
| $\mathbb Z_2$ | parzysty | parzysty |
| rzut | **izotropowy** | **anizotropowy** (bezśladowy, spin-2) |
| kinetyka | $K(\Phi)(\partial\Phi)^2$ z coarse-grainingu (op-gamma-RG-running) | $C_\sigma(\partial\sigma)^2$ z coarse-grainingu (**ten cykl**) |

> **Konkluzja §2 (LOCKED jako wejście Phase 2):** ponieważ $\sigma_{ab}$ jest kompozytem (jak $\Phi$),
> jego **kinetyka nie jest postulowana — pochodzi z dynamiki samego korelatora**. Współczynnik $C_\sigma$
> mierzy koszt energetyczny **przestrzennej zmienności korelatora kierunkowego** między blokami:
> $\partial_\mu\sigma_{ab}$ to gradient $\langle\hat s_i\hat s_j\rangle$ od bloku do bloku. To definiuje
> precyzyjnie obiekt, który Phase 2 ma wyekstrahować (§4).

## §3 — Weryfikacja strukturalna (sympy, `Phase1_sympy.txt`)
1. **σ_ab sym./bezśladowy/5-komp** — $\mathrm{Tr}\,\sigma=0$ EXACT; $6_{\rm sym}-1_{\rm ślad}=5$ = liczba komponentów spin-2. **PASS**.
2. **$\mathbb Z_2$-parzystość** — $\hat s\to-\hat s$: $\hat s_i\hat s_j$ parzyste; $\sigma_{ab}$ parzyste, spójne z symetrią substratu (rem:własności-σ).
3. **Bilans wymiarowy** $[S_\sigma]=1$ → $[C_\sigma][\sigma]^2=E^2$. Jeśli σ bezwymiarowe: $[C_\sigma]=E^2$ (skala napięcia²). Jeśli $[\sigma]=[\Phi_0]$: $[C_\sigma]=E^2/[\Phi_0]^2$. **Rozstrzygnięcie = Phase 3** (relacja $C_\sigma\leftrightarrow\sigma_0$).
4. **Znak** $C_\sigma>0\Leftarrow J>0$ — druga wariacja kosztu gradientu korelatora ferromagnetycznego $=+J>0$ → ghost-free (zgodne z eq:S-sigma bullet i prob:tensor-modes). **PASS**.

## §4 — Specyfikacja zadania dla Phase 2 (RDZEŃ — co dokładnie liczyć)
Obiekt do ekstrakcji: współczynnik przy $(\partial_\mu\sigma_{ab})^2$ w działaniu efektywnym powstałym z
blokowego uśredniania $H_\Gamma$. Operacyjnie — **gradient korelatora korelatorów**:
$$C_\sigma \;\propto\; \lim_{a_{\rm sub}\to0}\ \frac{\partial^2}{\partial(\partial_\mu\sigma)^2}\,\Big\langle\,H_\Gamma\,\Big\rangle_{\rm blok}\Big|_{\sigma(x)\ \text{wolnozmienne}}.$$
Konkretnie Phase 2 musi:
1. Promować $\sigma_{ab}(x)$ do pola wolnozmiennego; rozwinąć $K_{ab}(x+\delta)=K_{ab}(x)+\delta\cdot\partial K+\tfrac12\delta^2\partial^2K+\dots$ (gradient expansion bilinowego korelatora).
2. Wstawić do energii wiązaniowej i wycałkować po szybkich modach ŝ wewnątrz bloku (analog $\mathcal T_b$ z dodatekQ eq:blocking_operator, ale dla **anizotropowego** sektora $K_{ab}$, nie tylko izotropowego $\Phi_B$).
3. Zebrać współczynnik przy $(\partial\sigma)^2$ → $C_\sigma(J,a_{\rm sub},\dots)$.

**Szkielet skalowania (Phase 1; prefaktor = GAP):** analiza wymiarowa daje $C_\sigma\sim c\cdot J\cdot a_{\rm sub}^{\,p}$,
gdzie wykładnik $p$ (z liczby pochodnych w gradient expansion) i prefaktor $c=O(1)$ są **NIEOKREŚLONE w Phase 1**.
Ich wyznaczenie = Phase 2. **Żadnej liczby nie podaję** (anti-Lakatos; rodzic poprawnie zostawił to jako GAP).

## §5 — Dyspozycja F-CG-A (Lorentz) — najważniejszy test koncepcyjny
$H_\Gamma$ jest **statyczny** (chwilowe wiązanie sąsiadów) — generuje naturalnie sztywność **przestrzenną**
$(\nabla\sigma)^2$. Pełna kinetyka Lorentza wymaga członu **czasowego** $(\partial_t\sigma)^2/c_0^2$.
- **Co LOCKED:** rem:cGW-tensor stwierdza, że σ_ab propaguje z $c_0$ jako perturbacja NA geometrii ustanowionej przez $\Phi$ (Box$=c_0^{-2}\partial_t^2-\nabla^2$). To jest wejście rdzenia.
- **Co OTWARTE (R-Lorentz-emergence MED):** czy ten sam mechanizm emergentnej metryki, który dał skalarowi $\Box\Phi$ (przez GL-bond + dynamikę substratu ustanawiającą $c_0$), generuje człon czasowy dla **bezśladowego** $\sigma_{ab}$ — czy też σ_ab dostaje tylko sztywność przestrzenną (a wtedy $c_T\neq c_0$ i sprzeczność z rem:cGW-tensor / GW170817).
- **Status:** **PARTIAL-pending.** Phase 1 ustala, że pytanie jest dobrze postawione i że odpowiedź NIE jest automatyczna; rozstrzygnięcie (TAK → przejdź do B; NIE → GAP strukturalny) należy do Phase 2. **Nie deklaruję TAK przedwcześnie.**

## §6 — Inwentarz parametrów (LOCKED dla cyklu)
$$\boxed{\{\,J,\ a_{\rm sub},\ \mu,\ m_0^2,\ \lambda_0\,\}}$$
- $J$ — sprzężenie sąsiadów [energia] → **sztywność kierunkowa $C_\sigma$** (rola główna).
- $a_{\rm sub}$ — stała sieci [długość] → wykładnik skalowania $p$.
- $\mu$ — skala RG / wiązania (jak w op-gamma-RG-running, $G_0\leftrightarrow J\mu$ rem:param-counting).
- $m_0^2,\lambda_0$ — sektor on-site → **amplituda / VEV $v^2=|m_0^2|/\lambda_0$ → $\sigma_0$** (przez normalizację, Phase 3).

Budżet nowych stałych: **0** (wszystkie z eq:B-H / rem:param-counting).

## §7 — Anti-Lakatos (checklist Phase 1)
✓ Obiekt zidentyfikowany strukturalnie PRZED rachunkiem wartości (σ_ab = kompozyt, §2).
✓ **$C_\sigma$ NIE wyprowadzone — tylko szkielet skalowania + znak; prefaktor jawnie GAP** (§4) — zgodnie z rodzicem (F-CS-A GAP) i §5 forbidden moves.
✓ F-CG-A NIE zadeklarowane TAK przedwcześnie (§5, PARTIAL-pending) — uczciwe wobec R-Lorentz-emergence.
✓ Liczby poprzedników LOCKED; budżet stałych 0 (§6).
✓ Weryfikacja sympy dla wszystkich twierdzeń strukturalnych (§3); zero strojenia do 5/6.

## §8 — Handoff do Phase 2
**Phase 2 (RDZEŃ): coarse-graining / gradient expansion → $C_\sigma$.** Wejście: §4 (specyfikacja),
§5 (test Lorentza F-CG-A jako pierwsza bramka), §6 (inwentarz). Wzorzec: op-gamma-RG-running Phase 2
(skalarny $H_\Gamma\to S[\Phi]$, block-averaging $\mathcal T_b$) — **rozszerzony na sektor anizotropowy** $K_{ab}$.
Ryzyko: R-tensor-hard (HIGH) — możliwy PARTIAL (skaling bez prefaktora) lub GAP. Wymaga „działaj" usera.
