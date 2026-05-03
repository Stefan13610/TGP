---
title: "μ.1 PLAN — Minimal substrate redefinition (ψ = log g) jako droga do compound emergence"
date: 2026-05-02
cycle: μ.1
status: PLAN — propozycja cyklu, nieuruchomiona
parent: "[[../op-lambda1-e2-amplitude-emergence/README.md]]"
predecessor: "[[../op-lambda1-e2-amplitude-emergence/EXTERNAL_AUDIT_2026-05-02.md]]"
trigger: λ.1 NEGATIVE CLOSURE (6/6 mechanism failures) + użytkownik flagował compound interference jako structurally valid path
related:
  - "[[../op-lambda1-e2-amplitude-emergence/phase2_M6_compound_interference.py]]"
  - "[[../op-lambda1-e2-amplitude-emergence/EXTERNAL_AUDIT_2026-05-02.md]]"
tags:
  - TGP
  - mu1
  - substrate-redefinition
  - compound-interference
  - log-variable
  - plan-only
  - pre-implementation
---

# μ.1 — Minimal substrate redefinition (ψ ≡ log g)

> **Status:** PLAN. Nie uruchomiony. Decyzja go/no-go po analizie ryzyk i korzyści.
>
> **Trigger (z λ.1 closure):** wszystkie 6 testowanych mechanizmów field-theoretic
> dla derywacji X = e²/2 — NEGATIVE. M.6 (compound interference) wykluczył compound
> w obecnym substrate g(r) bo random multi-soliton placement daje liniową
> superpozycję, nie multiplicative.
>
> **Hipoteza μ.1:** TGP-substrate primitive variable nie jest g(r), tylko
> ψ(r) ≡ log g(r). Wtedy multi-soliton composition staje się addytywna w ψ
> (= multiplicative w g) automatycznie, bez dodawania reguł — co naturalnie
> generuje compound exp(Σε_i).

---

## 1. PROBLEM (do rozwiązania)

### 1.1 Co λ.1 udowodniło NEGATIVELY

Wszystkie 6 testowanych mechanizmów dla derywacji X = e²/2 (R3 mass formula
slope dla α=2 charged-leptons) — FAIL:

| # | Mechanism | Verdict | Issue |
|---|-----------|---------|-------|
| P2.1 | 1-loop log det O | NEG | K=−0.97 wrong sign + 4× factor |
| P2.2 | Semiclassical S_sol | NEG | K=−5.92 numerologia |
| P2.3 | Φ_eff stat-mech | NEG | anchor-dependent |
| M.4 | RG flow γ_φ 1-loop | NEG | 0.03 vs 3.69 → 125× miss |
| M.5 | AS NGFP η_φ | NEG | 7-170× miss + \|η\|<2 bound exceeded |
| M.6 | Compound interference | NEG | best I6≈2π z 4.46% drift; brak Σε=2 |

### 1.2 Co M.6 ujawniło strukturalnie

Compound formula `(1+x/N)^N → e^x` jest matematyczną tożsamością — **strukturalnie
zdolna** produkować e^x. ALE w obecnym R3 substrate:

- Random multi-soliton placement N=1..1000 daje **linear superposition** (Φ ~ N),
  nie compound exp
- Brak natural mechanizmu produkującego dokładnie **Σε = 2** dla α=2
- R3 ODE *ma* multiplicative prefactor g^(2α), ale kontroluje pojedynczy soliton,
  nie buildup z N

### 1.3 Diagnoza fundamentalnа

**Źródłem porażki M.6 jest wybór zmiennej superpozycji.** Jeśli pole g(r) jest
primitive, kompozycja N solitonów jest addytywna (g_total = Σg_i, linear).
Compound emergence wymaga multiplicative composition — co znaczy, że *jakaś inna*
zmienna jest pierwotnie liniowa.

**Naturalny kandydat:** ψ(r) ≡ log g(r). Jeśli ψ jest pierwotne, to:
- ψ_total = Σψ_i (liniowa superpozycja)
- g_total = exp(Σψ_i) = ∏g_i (multiplicative composition automatycznie)

**Pytanie μ.1:** Czy TGP-substrate jest naturalniej formułowany w ψ niż w g,
i czy ta redefinicja zachowuje wszystkie wcześniejsze sukcesy TGP-programu?

---

## 2. POTENCJALNA ZMIANA — co dokładnie się zmienia

### 2.1 Substrate primitive (formal)

| Aspekt | Obecny substrate | μ.1 substrate |
|--------|-----------------|--------------|
| Pole pierwotne | g(r) ∈ ℝ⁺ | ψ(r) = log g(r) ∈ ℝ |
| Zakres | g > 0 | ψ ∈ (−∞, ∞) |
| Vacuum | g = 1 | ψ = 0 |
| Kompozycja multi-soliton | g_total = Σg_i (linear) | ψ_total = Σψ_i (linear) |
| Konsekwencja | ε_total = Σε_i | g_total = exp(Σε_i) (compound) |

### 2.2 R3 ODE — transformacja

Obecna ODE (pre-redefinition):
$$g'' + \frac{\alpha}{g}(g')^2 + \frac{2}{r}g' = (1-g)\, g^{2-2\alpha}$$

Pod g = e^ψ, g' = e^ψ · ψ', g'' = e^ψ(ψ'' + (ψ')²):

$$e^{\psi}(\psi'' + (\psi')^2) + \alpha\, e^{\psi}(\psi')^2 + \frac{2}{r}e^{\psi}\psi' = (1 - e^{\psi})\, e^{(2-2\alpha)\psi}$$

Dzieląc przez e^ψ:

$$\boxed{\psi'' + (1+\alpha)(\psi')^2 + \frac{2}{r}\psi' = (e^{-\psi} - 1)\, e^{(1-2\alpha)\psi}}$$

**Uwagi:**
- ODE pozostaje **rozwiązywalna** (smooth, well-posed)
- Vacuum ψ = 0 jest fixed point (RHS znika)
- Dla α=2: RHS = (e^{-ψ} - 1)·e^{-3ψ}
- Kompresja informacji: trzeba sprawdzić czy soliton solution w ψ ma takie samo
  asymptotic decay jak w g

### 2.3 Mass formula — INVARIANT

$$m = c \cdot A^2 \cdot g_0^X = c \cdot A^2 \cdot \exp(X \cdot \psi_0)$$

**Formalnie identyczne.** Zmiana zmiennej tylko relabeluje X-interpretację:
- W g-zmiennej: X = pochodna log(m/A²) względem log(g₀)
- W ψ-zmiennej: X = pochodna log(m/A²) względem ψ₀

Numerycznie: **ten sam X**. **0.0007% μ/e match preserved.**

### 2.4 Compound emergence — analitycznie

Multi-soliton background w ψ-zmiennej:
$$\psi_{\text{bg}}(\vec{x}) = \sum_{i=1}^{N} \psi_i(\vec{x} - \vec{x}_i)$$

Tail każdego solitonu (linearyzacja przy ψ → 0): ψ_i ≈ ε_i.

Wartość tła w punkcie central: ψ_total = Σε_i.
$$g_{\text{bg}} = \exp(\psi_{\text{total}}) = \exp(\Sigma\epsilon_i)$$

**Jeśli** topologiczny constraint daje **Σε = 2** (np. winding number, total charge, spinor structure factor):
$$g_{\text{bg}} = e^2$$

**Compound emergence z e²** wtedy automatyczna **bez nowych reguł**.

### 2.5 Co jest do udowodnienia (kluczowy open)

**Σε = 2 dla α=2 charged-lepton** musi mieć topologiczne/algebraiczne uzasadnienie.
Kandydaci:
- Total electric charge (charged lepton: ±1, ale Σε=2 wymaga factor 2)
- Spinor double-cover: ψ → ψ + 2π pod 2π rotation (factor 2 z spin-1/2)
- Winding number n=2 z R3 topologii dla α=2
- Compound interest natural saturation point (Σε = 2 = "doubling the
  background field" → e²)

**μ.1 nie wybiera apriori — testuje który z tych argumentów ma rację.**

---

## 3. KORZYŚCI (jeśli się powiedzie)

### 3.1 Bezpośrednie

1. **λ.1 zmienia status z NEGATIVE → PARTIAL/POSITIVE CLOSURE** —
   compound mechanism dla X = e²/2 zostaje znaleziony bez dodawania reguł
2. **Naturalne wyjaśnienie Φ² w wielu TGP równaniach** — Φ² = exp(2ψ),
   czyli "doubling background" jest natural compound saturation
3. **Mass formula μ/e 0.0007% match dostaje fizyczne wytłumaczenie**
   (nie empirical fit, ale topological consequence)

### 3.2 Strukturalne dla TGP-programu

1. **Substrate staje się bardziej fundamental** — ψ jest bezdimentional log,
   bardziej "geometric" niż g
2. **Bridge theorem α=1 IFF** może być reinterpreted: dla α=1 ODE w ψ
   ma specjalną właściwość (RHS = (e^{-ψ}-1)·e^{-ψ}) — może być reason
   dlaczego α=1 jest critical case
3. **Amplitude/phase split** (Section 11.5 audytu) zyskuje fundament:
   amplitude = exp(ψ), phase = θ, naturalne że amplitude komponuje multiplicative
4. **Powiązania z RG** — w Wilsonian RG primitive są often log Z; ψ-substrate
   zgodny z tym
5. **Powiązania z QFT** — Born rule \|amplitude\|², gdzie amplitudy mnożą się
   dla composition — naturalnie multiplicative w amplitude domain

### 3.3 Predykcyjne

- **τ/μ ratio** powinien spełnić tę samą compound formułę z innym Σε
- **Inne sektory** (cosmology, BH thresholds) mogą zyskać clean explanations
- **NEW prediction:** ψ-formulation wskazuje gdzie e^n (dla n ≠ 2) może appear
  w innych TGP observables — testable by inspection

---

## 4. RYZYKA (uczciwa lista)

### 4.1 Krytyczne — mogą uniemożliwić cykl

**R1. Σε = 2 nie ma topologicznego źródła.**
   Jeśli żaden z kandydatów (charge, winding, spinor) nie wymusza Σε = 2,
   μ.1 redukuje się do reparametryzacji bez fizycznych konsekwencji.
   **Prawdopodobieństwo:** średnie (40%). **Mitigacja:** P4 jest pierwszym
   testowym etapem; jeśli fail, μ.1 closes cleanly bez kosztu.

**R2. ODE w ψ nie ma physical soliton solutions.**
   Transformowane RHS (e^{-ψ}-1)·e^{(1-2α)ψ} może mieć inne asymptotyki niż
   g-ODE; possibly singularity przy ψ → -∞ (g → 0).
   **Prawdopodobieństwo:** niskie (15%). **Mitigacja:** P1 testuje to numerycznie.

**R3. Bridge theorem α=1 IFF łamie się pod redefinicją.**
   Bridge dotyczy R5↔Phase 2 isomorfizmu w określonej geometrii. Pod ψ-redef.
   ten isomorfizm może wymagać przedefinowania albo padać.
   **Prawdopodobieństwo:** niskie (10%) — bridge jest geometryczny, redef. skalarna
   nie powinna tknąć. **Mitigacja:** P2 sprawdza explicite.

### 4.2 Duże — degradują wartość

**R4. M.4 retest w ψ nadal NEGATIVE.**
   Anomalous dimension γ_φ rescaluje pod field redef., ale niekoniecznie
   o factor 125×. Może być nadal big miss.
   **Prawdopodobieństwo:** średnie (35%). **Konsekwencja:** μ.1 daje kompozycję
   compound, ale nie tłumaczy *wartości* X = e²/2 z RG. Lokalny success.

**R5. M.5 (AS NGFP) zostaje NEGATIVE.**
   Bound \|η\|<2 jest field-redef-invariant. Targety η nie zmienią się pod ψ.
   **Prawdopodobieństwo:** wysokie (70%). **Konsekwencja:** akceptujemy, że
   AS NGFP nie jest właściwym frameworkiem; μ.1 nie żre nic z M.5.

### 4.3 Małe — niepokoją ale nie blokują

**R6. Numerical precision μ/e zostaje invariant ale interpretacja może się zmienić.**
   0.0007% match jest robusto-invariant pod zmianą zmiennej. Ale fizyczna
   interpretacja "dlaczego ten X" jest inna.
   **Prawdopodobieństwo:** pewne (100%) — to nie ryzyko, to feature.

**R7. Praca pełnego retest of P2.1, P2.2 może wymagać Jacobianu pola.**
   1-loop i semiclassical nie są naiwnie field-redef-invariant — Jacobian
   determinant musi być traktowany.
   **Prawdopodobieństwo:** wysokie (80%) **dla retest**. **Mitigacja:** retest
   nie jest core, ale optional w P5.

### 4.4 Hidden — możemy nie zauważyć

**R8. Compound interpretation może być over-fitting.**
   Σε = 2 daje e², ale Σε = 1 daje e, Σε = 3 daje e³. Bez constraint top.
   to jest po prostu wybór. Naukowo: czy to nie "moving target"?
   **Prawdopodobieństwo:** zawsze present. **Mitigacja:** P4 musi mieć FALSIFICATION
   criterion — jeśli żaden topology daje Σε = 2 dla charged-lepton, μ.1 fail.

---

## 5. CYKL μ.1 — proponowana struktura

### Phase 1 (Foundation) — ψ-substrate validity

**P1.1 Sympy derivacja ψ-ODE.** Analitycznie wyprowadzić ODE w ψ z g-ODE.
Sprawdzić smooth, well-posed, fixed point ψ=0 stabilny.

**P1.2 Numeryczne soliton solutions w ψ dla α=2.** Solver ODE w ψ, plot tail
asymptotics, porównać z g-soliton po transformacji ψ = log g. Powinny być
*identical* (to test consistency).

**P1.3 Mass formula in ψ.** Re-fit log(m/A²) = X·ψ₀ + const. Sprawdzić X
identyczny z g-fit. **GATE:** μ/e match musi pozostać 0.0007%.

**Score:** 3/3 PASS minimum. Jeśli fail → μ.1 closes (nie ma sensu kontynuować).

### Phase 2 (Compound emergence) — kluczowy test

**P2.1 Multi-soliton ψ-superposition (analitycznie).** Pokazać formalnie że
ψ_total = Σψ_i ⟹ g_total = exp(Σε_i) w linearization.

**P2.2 Topological constraint dla Σε = 2.** Cztery kandydaci:
   - (a) Electric charge: charged lepton ±1, możliwe że factor 2 z
     "particle + antiparticle in vacuum" = 2
   - (b) Spinor double-cover: rotation 2π → ψ → ψ + 2π zamiast ψ ; dla
     spin-1/2 trzeba 4π dla identity; ratio 2
   - (c) R3 winding number: dla α=2 (3D), homotopy π_2(S²) = ℤ, fundamental
     winding n=1, ale "double winding" n=2 jako bound-state structure
   - (d) Compound saturation: substrate ma natural cutoff Σε = 2 (factor 2 z
     "background + perturbation = doubling")

**Każdy kandydat:** sympy derivation lub numerical demonstration. Co najmniej
JEDEN musi dać Σε = 2 z first principles. **GATE:** at least 1/4 pass.

**P2.3 Verify compound z analytic.** Σε_topology = 2 ⟹ g_bg = e² ⟹ X =
log(m/A²)/ψ₀ = e²/2 dla α=2 charged-lepton. **GATE:** 0.5% drift od 3.69453.

### Phase 3 (Cross-validation) — czy nic nie pęka

**P3.1 Bridge theorem α=1 IFF w ψ-substrate.** Re-derive sympy. **GATE:** musi
przetrwać.

**P3.2 R3 sub-tension τ canonical.** g₀_τ = 1.77472 ⟹ ψ₀_τ = 0.5736. Sprawdzić
czy τ/μ ratio mass formula daje tę samą precision pod compound.

**P3.3 M.4/M.5 retest w ψ.** Optional ale informative. Czy γ_φ rescaluje? Czy
AS NGFP bound nadal binding?

**P3.4 Compound predykcja dla τ/μ.** Σε_τμ powinien być pewną liczbą
naturalną (e.g., 2/3? 3/2? lub Σε_τ - Σε_μ = ...?). **PREDYKCJA NEW:** sprawdzić.

### Phase 4 (Consolidation) — verdict

**P4.1 Status synthesis** dla μ.1 + λ.1.
**P4.2 If POSITIVE:** reopen λ.1, document mechanism, propagacja do innych
TGP-cykli (cosmology, BH, etc.)
**P4.3 If NEGATIVE:** document why, archive μ.1 jako wartościowy negative,
λ.1 zostaje closed.

---

## 6. KOSZT — estimated effort

| Phase | Estimated time | Critical resources |
|-------|---------------|-------------------|
| P1 (Foundation) | 0.5 dnia | sympy, scipy ODE solver |
| P2 (Compound) | 1 dzień | sympy, topology argumenty |
| P3 (Cross-val) | 0.5 dnia | reuse λ.1 numerics |
| P4 (Consolidation) | 0.5 dnia | docs |
| **Total** | **~2-3 dni focused** | |

---

## 7. DECYZJA — kryterium go/no-go

**GO** jeśli:
- Tutaj (PLAN) wszystko strukturalnie się spina ✓ (zakładamy review przed start)
- P1.3 GATE μ/e 0.0007% przeżywa (almost guaranteed by reparametryzation theory)
- P2.2 ma co najmniej 1/4 topological argument PASS

**NO-GO** jeśli:
- Już P1 ujawnia że ψ-ODE ma singularities albo nie ma soliton solutions
- P2.2 fail wszystkie 4 kandydatów dla Σε=2 — wtedy μ.1 = pure reparametryzacja
  bez fizycznych konsekwencji, nie warta cyklu

**SOFT NO-GO** (revisit later):
- Jeśli P3.3 (M.4/M.5 retest) pokazuje że nadal big mismatches w innych
  obszarach — μ.1 lokalny success, ale TGP-program ma większy strukturalny
  problem; wymaga jeszcze innego cyklu

---

## 8. RELACJA do λ.1 i innych cykli

- **λ.1** zostaje **NEGATIVE CLOSURE** *unchanged* przez plan μ.1.
  μ.1 jest **rozszerzeniem ortogonalnym**, nie reopening λ.1 (jeszcze).
- **Jeśli μ.1 POSITIVE** → λ.1 README dostaje postscript Section 13:
  "Mechanism found via μ.1 substrate redefinition; see [[../op-mu1.../README.md]]"
- **G1** (Φ_eff anchor resolution) — pierwotny next-target po λ.1. μ.1 może
  iść **przed** G1 jeśli go-decision, lub **po** G1 jako odkładany.

---

## 9. STATUS PLANU

**Niniejszy dokument:** PLAN tylko. **Nieuruchomiony.**

**Następna akcja:** review przez użytkownika; jeśli plan się "spina" (per user)
→ kick-off Phase 1.

**Pliki μ.1 do wytworzenia podczas implementacji** (nieobecne teraz):
- `phase1_psi_ode_derivation.py`
- `phase1_psi_soliton_numeric.py`
- `phase2_topological_sigma_eps_test.py`
- `phase3_compound_emergence_proof.py`
- `README.md` (after P4)

---

## 10. SAMOOCENA planu (uczciwa)

**Strengths:**
- Najmniej inwazyjna zmiana z 4 rozważanych klas (A: ψ=log g, B: amplitude/phase
  split, C: multiplicative manifold, D: ad hoc rule). Klasa A wybrana jako
  **reparametryzacja** nie nowa reguła
- Wszystkie 6 wcześniejszych negatywów λ.1 są zachowane jako diagnoses;
  μ.1 nie negate ich, tylko proponuje rebase substrate
- GATE-driven structure z falsification criteria w każdej fazie

**Weaknesses:**
- **R1 (Σε=2)** jest serce success-condition; mam 4 kandydatów ale żaden
  nie jest jeszcze przekonujący — to MOŻE być fundamental obstacle
- M.4/M.5 retest mogą zostać NEG nawet w ψ, co znaczy że μ.1 daje kompozycję
  ale nie *wartość* X — partial success only
- "Reparametryzacja" jest filozoficznie cienkie — czy naprawdę zmieniliśmy
  substrate, czy tylko zmienne? User flagował to: "minimalna zmiana substratu" —
  granica pomiędzy reparametryzacją a zmianą substratu jest fuzzy

**Honest verdict on plan itself:** jest **uczciwy**, ma jasne GATEs, kosztowo
2-3 dni. Główne ryzyko: R1 (Σε=2 brak topology) — bez tego μ.1 redukuje się do
relabeling. **Powinien być zrobiony jeśli user chce zamknąć compound path
formalnie**, niezależnie od wyniku.
