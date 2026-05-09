---
title: "λ.1 — e² jako fundamental constant w TGP amplitude sector?"
date: 2026-05-01
cycle: λ.1
status: NEGATIVE CLOSURE (post-audit + M.4 + M.5 + M.6) — wszystkie 6 testowanych mechanizmów FAIL; e² zostaje empirical coincidence
phase1_score: 3.0/5 PASS
phase2_score: 0.5/4 FAIL (specific mechanisms)
phase3_score: 2.5/3 PASS (synthesis + bridge cross-validation + sympy lock)
m4_score: NEGATIVE (RG flow γ_φ 1-loop ~0.03 vs target 3.69; 125x miss)
m5_score: NEGATIVE (AS NGFP η_φ ±0.05-0.5 vs target 3.69; 7-170x miss + perturbativity bound exceeded)
m6_score: NEGATIVE (compound interference Φ_total = ∏(1+ε_i); single soliton best I6≈2π z 4.46% drift, multi-soliton random=linear superposition; brak natural Σε=2 w R3 substrate)
overall_verdict: NEGATIVE CLOSURE z wartościowymi wykluczeniami (6/6 mechanizmów odrzuconych — perturbative + AS NGFP + compound)
external_validation: mass_scaling_k4 R5↔Phase2 bridge dotyczy α=1; R3 używa α=2 — bridge nie waliduje R3
parent: "[[../why_n3/README.md]]"
predecessors:
  - "[[../why_n3/PHASE2_n_alpha_derivation.md]]"
  - "[[../why_n3/PHASE6_alpha_em_connection.md]]"
  - "[[../why_n3/r3_phase7_phi0_screening_e2.py]]"
related:
  - "[[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]"
  - "[[../mass_scaling_k4/RECONCILIATION_R5_vs_phase2_2026-04-30.md]]"
  - "[[../mass_scaling_k4/g0_tau_subtension_diagnostic.py]]"
tags:
  - TGP
  - lambda1
  - e-squared
  - euler-number
  - amplitude-sector
  - NEGATIVE-CLOSURE
  - mechanism-falsified
  - empirical-coincidence
tgp_status:
  folder_status: paused
  level: L1
  kind: derivation
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'λ.1 — e² jako fundamental constant w TGP amplitude sector'"
    - "Phase{1..3}_results.md PASS=155, CLOSED=7, FAIL=26"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

> ## ⛔ FINAL STATUS 2026-05-02 (post-audit + M.4 + M.5 + M.6) — NEGATIVE CLOSURE
>
> **Po external audit + M.4 RG flow + M.5 AS NGFP + M.6 compound interference test (2026-05-02):**
>
> Wszystkie SZEŚĆ testowanych mechanizmów dla derywacji e²/2 — NEGATIVE:
>
> | Mechanizm | Wynik | Target | Status |
> |-----------|-------|--------|--------|
> | P2.1 1-loop log det O | K = -0.97 | +3.69 | wrong sign + factor 4× |
> | P2.2 Semiclassical S_sol | K = -5.92 | ±3.69 / ±7.39 | numerologia |
> | P2.3 Φ_eff stat-mech | anchor-dependent | Brannen 24.78 | nie matche |
> | **M.4 RG flow γ_φ** | **0.03 (1-loop)** | **3.69** | **125× miss** |
> | **M.5 AS NGFP η_φ** | **±[0.05-0.5]** | **3.69** | **7-170× miss + |η|<2 bound exceeded** |
> | **M.6 Compound interference** | **best I6=6.00 ≈ 2π** | **e=2.72, e²=7.39, e²/2=3.69** | **>4% drift, brak natural Σε=2** |
>
> **Systematyczny problem skali (uściślony):**
> - e²/2 ≈ 3.69 jest **rzędu wielkości większe** niż cokolwiek naturalnego z
>   *single-process* perturbative QFT (~1/(16π²) ≈ 0.006) lub AS NGFP η_φ (~0.1-0.5).
> - Compound interference (Φ = ∏(1+ε_i) → exp(Σε)) jest *strukturalnie* zdolne
>   produkować e^x **trywialnie** (to identity, nie hipoteza), ale w R3 substrate
>   random multi-soliton placement daje **linear superposition** (Φ ~ N), nie
>   compound exp. Brak natural TGP-mechanizm produkujący Σε = 2 dla α=2.
>
> **Werdykt audytu:** X = e²/2 zostaje **empirical fit z 0.0007% numerical
> coincidence** z exp(2)/2, **bez derivation** ani z perturbacyjnego, ani z
> NGFP, ani ze strukturalnego compound path. Cross-validation z mass_scaling_k4
> jest częściowo cyrkularna (oba używają tego samego μ/e fit) i bridge theorem
> dotyczy α=1 podczas gdy R3 mass formula używa α=2.
>
> **λ.1 ZAMKNIĘTE OSTATECZNIE: NEGATIVE CLOSURE z 6 mechanizmami wykluczonymi.**
>
> **Pełny audyt + M.4 + M.5 + M.6:** `[[EXTERNAL_AUDIT_2026-05-02.md]]`
> **M.4 test:** `[[phase2_M4_rg_flow_gamma_phi.py]]`
> **M.5 test:** `[[phase2_M5_as_ngfp_eta_phi.py]]`
> **M.6 test:** `[[phase2_M6_compound_interference.py]]`
>
> **Post-script (2026-05-02 noc):** μ.1 substrate redefinition cycle
> uruchomiony (ψ ≡ log g hypothesis) — **NO-GO CLOSURE**: reparametryzacja PASS
> mathematically (1e-13% μ/e drift), ale Σε=2 topology FAIL (0/4 candidates).
> μ.1 ≡ pure relabeling, nie reverses λ.1 NEG. **Total: 7/7 ścieżek wykluczonych.**
> Pełny cykl: [[../op-mu1-minimal-substrate-log-redefinition/README.md]]

---

> ## ✅ FINAL STATUS 2026-05-02 — PARTIAL CLOSURE (post-Phase 3)
>
> **3 fazy zamknięte:**
> - **Phase 1**: 3.0/5 PASS (foundation — hipoteza internally consistent)
> - **Phase 2**: 0.5/4 FAIL (specific field-theoretic mechanisms wykluczone)
> - **Phase 3**: 2.5/3 PASS (synthesis + bridge cross-validation + sympy lock)
> - **Total weighted: 50%** (honest mixed result)
>
> **CO ZOSTAŁO UDOWODNIONE:**
> - ✓ **e² = Euler² strukturalnie zidentyfikowane** w R3 mass formula
>   z **0.0007% match** (cross-validated z mass_scaling_k4)
> - ✓ **Bridge theorem**: R5 K² ≡ Phase 2 IFF α=1 (closed-form sympy proof)
> - ✓ **Mass formula PDG match**: μ/e ±0.000%, τ/e +0.006%, τ/μ +0.015%
> - ✓ **Strukturalna izolacja** amplitude/phase sector (phase wyklucza e_Euler)
> - ✓ **Sub-tensja τ CLOSED**: g₀_τ canonical = 1.77472 (NIE 1.75505)
>
> **CO POZOSTAJE OPEN:**
> - ✗ Konkretny field-theoretic mechanizm produkujący e² z TGP-substrate
>   (Phase 2 wykluczyła log det O, semiclassical, Φ_eff stat-mech)
> - ✗ Wykładnik 2 origin (czemu e², nie e¹/e³/e⁴)
> - ✗ (10/3) w Φ_eff jest anchor-dependent (NIE fundamental)
> - ✗ Cross-cycle universality (neutrina mają osobny mechanizm)
>
> **STATUS PUBLICATION-READY (z phase3_P31_synthesis):**
> R3 mass formula z `e² = exp(2) = Euler²` jest **strukturalnie locked**;
> konkretny derivation mechanizm pozostaje OPEN. Phase 2 negatives są
> **wartościowe** jako wykluczenie specific approaches.
>
> **Patrz `Phase3_results.md` dla pełnego verdict; `MASS_SCALING_K4_CROSS_VALIDATION.md`
> dla bridge implikacji.**

> ## 🔄 STATUS 2026-05-02 — PAUSED (nie zamknięte)
>
> **Phase 1:** 3.0/5 PASS (gate passed)
> **Phase 2:** 0.5/4 FAIL — konkretne mechanizmy testowane (log det O, S_sol,
> Φ_eff substrate stat-mech) **NEGATIVE** dla derywacji e²/4.
>
> **ALE: λ.1 hipoteza NIE jest sfalsyfikowana** — odkrycie zewnętrzne 2026-05-02:
>
> ### Cross-validation z `research/mass_scaling_k4/` (2026-05-02)
>
> Niezależny cykl `mass_scaling_k4` (R5 K² mechanism + analytical bridge)
> dostarczył **silne wsparcie** dla λ.1 hipotezy:
>
> 1. **e² = Euler² = exp(2) = 7.389056** — potwierdzone z μ/e exact fit do
>    **0.0007%** (jeszcze czystsze niż wcześniejsze 0.001%). To NIE jest
>    empirical fit — to **strukturalna identyfikacja** (R5_PHASE2_ANALYTICAL_BRIDGE).
>
> 2. **Sub-tensja τ closed**: pełna formuła Phase 2 dla TGP-canonical α=2 +
>    Koide K=2/3 daje **g₀_τ = 1.77472** i **+0.006% PDG** dla m_τ/m_e
>    (residue -0.085% z r3_alpha2_full_closure.py był artefaktem A³ skrótu).
>
> 3. **Analytical theorem**: R5 K² mass formula ≡ Phase 2 m_obs **IFF α=1**
>    (closed-form proof). Phase 2 jest **fundamental**, R5 K² to derivative.
>
> 4. **Bridge implikacja**: e² wynika strukturalnie z samej Phase 2 universal
>    mass formula, nie wymaga konkretnego mechanizmu typu log det O lub
>    Φ_eff substrate.
>
> ### Re-assessment Phase 2 negative verdict
>
> Moje Phase 2 testy zostają **POPRAWNE jako konkretne wyniki**:
> - P2.1 (log det O fit K=-0.97) — log det NIE jest mechanizm dla e²/4
> - P2.2 (S_sol fit K=-5.92) — semiclassical NIE matche e²-family
> - P2.3 (Φ_eff Brannen mismatch) — Φ_eff·(10/3) jest anchor-dependent
>
> **Te testy stoją.** ALE moja generalizacja "e² nie jest fundamental w TGP"
> była **przedwczesna** — mass_scaling_k4 niezależnie pokazało że e² jest
> fundamental przez **inny** mechanizm (Phase 2 universal mass formula
> structure + R5 bridge).
>
> ### Status λ.1: PAUSED, nie ENDED
>
> - Hipoteza λ.1 ("e² fundamental w TGP amplitude sector") jest **żywa**
> - Konkretne mechanizmy testowane w Phase 2 (log det, semiclassical, Φ_eff)
>   wykluczone — to wartościowe negative results
> - Nowy mechanizm (R5↔Phase 2 bridge analytical theorem) jest **już zrobiony**
>   w mass_scaling_k4 — λ.1 może być wznowiona z tego punktu
> - **Phase 3 może być uruchomiona w przyszłości** jako synteza z mass_scaling_k4
>
> ### Status X = e²/4 (post-cross-validation)
>
> - **Wzmocnione**: e² = Euler² potwierdzone strukturalnie (0.0007% match)
> - **Wzmocnione**: g₀_τ = 1.77472 daje +0.006% PDG (sub-tensja closed)
> - **Otwarte**: konkretny **mechanizm** generujący e² z TGP-substrate (Phase 1+2
>   testy negative wykluczyły log det O, semiclassical, Φ_eff stat-mech)
> - **Ścieżka zamknięcia**: synteza R5 K² ≡ Phase 2 (α=1) bridge + Phase 2
>   universal formula (od which e² jest "natural fit-constant" z μ/e)
>
> **Patrz `Phase2_results.md` dla mechanizm-by-mechanizm assessment**
> i `MASS_SCALING_K4_CROSS_VALIDATION.md` dla bridge implikacji.

---

# λ.1 — e² jako fundamental constant w TGP amplitude sector

> **Cel:** Zbadać czy liczba Eulera (e ≈ 2.71828) jest **fundamentalnie zaszyta**
> w TGP amplitude sector, czy jest tylko empirycznym fitting parameter. Trzy
> niezależne poszlaki numerczne wskazują na realny strukturalny ślad e²;
> ten cykl ma za zadanie albo **zamknąć derywację** z TGP-fundamentu, albo
> jednoznacznie **odrzucić** jako numerologiczną zbieżność.

---

## 1. Streszczenie hipotezy

W trakcie sesji R3 (`research/why_n3/`) 2026-05-01 odkryto **trzy niezależne
numeryczne hint'y** wskazujące że e² (= 7.389) ma realne miejsce w TGP
amplitude sector:

1. **R3 mass formula** dla wszystkich α: `n(α) = (e²/4)·(4-α)` z residuum
   <0.1% dla α ∈ [0.25, 4.0]
2. **PDG match** dla μ/e: m_μ/m_e diff **-0.001%** używając n(2) = e²/2
3. **Φ_eff cosmological**: Φ_eff ≈ (10/3)·e² z diff **0.12%**

Hipoteza λ.1: te trzy hint'y nie są niezależnymi przypadkami — istnieje
**fundamentalna struktura** w TGP-substracie która generuje e² jako
emergent constant w amplitude sector (sektor `J_amp`, NIE `J_phase`/EM).

---

## 2. Background — skąd wzięły się poszlaki

### 2.1 R3 mass formula closure (Faza 2 z why_n3)

R3 cykl ma empirycznie znalezioną mass formula dla TGP-canonical α=2:

```
m_obs(g₀, α) = c_M · A_tail²(g₀, α) · g₀^[n(α)]
gdzie n(α) = e² · (1 - α/4) = (e²/4) · (4-α)
```

Numerycznie zweryfikowane na rozszerzonym zakresie α ∈ [0.25, 4.0]
(14 punktów):

| α | n(α) numerical | (e²/4)·(4-α) | residuum |
|---|----------------|---------------|----------|
| 0.25 | 6.940 | 6.927 | -0.18% |
| **2.00** | **3.6946** | **3.6945** | **-0.001%** |
| 4.00 | -0.006 | 0.000 | numerical zero |

**Source:** `research/why_n3/r3_phase2_n_alpha_derivation.py`

### 2.2 PDG mass ratios match (Faza 5 z why_n3)

Dla α=2 (TGP-canonical, K(φ)=φ⁴):

| Cząstka | TGP prediction | PDG | diff |
|---------|----------------|-----|------|
| m_μ/m_e | 206.77 | 206.7682 | **-0.001%** |
| m_τ/m_e | 3474.28 | 3477.23 | -0.085% |
| m_τ/m_μ | 16.820 | 16.817 | +0.015% |

Wszystkie trzy ratio match <0.1%. Jeśli e²/4 jest tylko empiryczny fit,
match dla μ/e powinien być rozproszony (~0.5%); fakt że jest na poziomie
**0.001%** sugeruje strukturalne pochodzenie.

**Source:** `research/why_n3/r3_phase5_full_propagator.py`

### 2.3 Φ_eff cosmological match (Faza 7 z why_n3)

Z `core/sek00_summary.tex:74-78`:
```
Φ₀ (bare)  = 168·Ω_Λ ≈ 115        [UV scale]
Φ_eff      = Φ₀ · 3/14 ≈ 24.66    [IR, "ekranowany dielektryk"]
```

Numeryczna obserwacja:
```
Φ_eff (TGP, kosmologiczne) = 24.66
(10/3) · e²                = 24.63
                       diff = +0.12%
```

To by oznaczało:
- e² ≈ (3/10) · Φ_eff
- X = e²/4 ≈ (3/40) · Φ_eff

**Source:** `research/why_n3/r3_phase7_phi0_screening_e2.py`

### 2.4 Brakujący link

Te trzy hint'y **nie są skleane przez derywację** w TGP — każdy jest
suggestive samodzielnie, ale brakuje:
- Analitycznej derywacji `e²/4` z R3 ODE lub TGP-action
- Wyjaśnienia czemu Φ_eff = (10/3)·e² (factor 10/3 niewyjaśniony)
- Mechanizmu który tłumaczy że e² żyje w amplitude sector ale **nie**
  w phase/EM sector (gdzie α_em używa Φ₀, J, π — bez e_Euler)

---

## 3. Cel cyklu

### 3.1 Główny cel

**Wyprowadzić e² z fundamentu TGP** w amplitude sector:
1. Pokazać konkretny mechanizm matematyczny dający e² (RG flow, partition
   function, kumulatywny limit)
2. Zamknąć analitycznie wzór n(α) = (e²/4)·(4-α) (lub odrzucić)
3. Zamknąć analitycznie Φ_eff = (10/3)·e² (lub odrzucić)
4. Wyjaśnić selektywność: e² w amplitude, brak w phase

### 3.2 Falsyfikatory

- Jeśli pełna derywacja R3 mass formula daje **inny** współczynnik niż e²/4
  z residuum <0.05% (np. 37/20, (3+e·φ)/4) — λ.1 FALSIFIED jako "e² to
  numerologia"
- Jeśli Φ_eff dokładna wartość okaże się różna od (10/3)·e² o >0.5% po
  precyzyjniejszych pomiarach Ω_Λ — λ.1 FALSIFIED dla cosmological match
- Jeśli explicit mass-formula derivation z R3 ODE daje funkcję NIE-liniową
  w α (kontradiktoryjne z `n(α) = X·(4-α)` z X=const) — λ.1 FALSIFIED

### 3.3 Possible outcomes

- **CASE A — Full closure:** e² wyłania się z TGP-fundamentu, λ.1
  promotuje X od EMPIRICAL do DERIVED, R3 mass formula zamknięta.
- **CASE B — Partial closure:** wyłania się **inna** stała (np. (3+e·φ)/4
  lub 37/20) lepiej dopasowana po precyzyjniejszej derywacji, e²
  wycofane.
- **CASE C — Zarzucone:** brak natural source dla e² w TGP, hint'y
  pozostają numerologicznymi zbieżnościami, R3 X = 1.847 jako empirical
  bez głębszej struktury.

---

## 4. Poszlaki (zebrane wszystkie znalezione)

### 4.1 POZYTYWNE — wskazujące na e² fundamentalność

| # | Poszlaka | Match | Source |
|---|----------|-------|--------|
| P1 | n(α) = (e²/4)·(4-α) fit | <0.1% dla α∈[0.25,4.0] | `why_n3/r3_phase2b_X_constant.py` |
| P2 | m_μ/m_e z n(2)=e²/2 | -0.001% PDG | `why_n3/r3_phase5_full_propagator.py` |
| P3 | Φ_eff ≈ (10/3)·e² | +0.12% (cosmological 36·Ω_Λ) | `why_n3/r3_phase7_phi0_screening_e2.py` |
| P4 | n(4) = 0 EXACT | numerical zero | Hobart-Derrick balance point |
| P5 | e²/4 vs alternativnych X (37/20, (3+eφ)/4) | wygrywa avg residuum | `why_n3/r3_phase2b_X_constant.py` |

### 4.2 NEGATYWNE — wykluczające pewne ścieżki

| # | Negatywne ustalenie | Wniosek |
|---|---------------------|---------|
| N1 | α-em używa 8π, Φ₀, J (BRAK e_Euler) | e² **NIE** pochodzi z phase sector / EM |
| N2 | 5D Kaluza-Klein nie generuje e² trywialnie | R⁵-bridge wycofany (extension nie-TGP) |
| N3 | 3/14 screening = 12/56 (algebraic) | Φ₀ screening pochodzi z power-law algebra, **NIE** e-related |
| N4 | TGP-natural units nie skleają X = α (różnica 4 rzędów) | X **nie jest** α w innych jednostkach |

### 4.3 STRUKTURALNE — TGP-mechanisms które MOGĄ generować e

| # | Mechanizm | Status w TGP |
|---|-----------|--------------|
| M1 | RG flow Z(μ) = exp(∫γ dlog μ) | TGP ma RG (UV.1 NGFP), niewykorzystane do R3 |
| M2 | Partition function Z = Σ exp(-βE) | TGP-substrate ma stat-mech, niezbadane explicit dla R3 |
| M3 | Wave-function renormalization Z_φ | Standardowy QFT, możliwe w R3 amplitude sector |
| M4 | Kumulatywny soliton dressing g(r+dr) = g(r)·(1+Δ/n)^n → exp(Δ) | Hipoteza — niezweryfikowana w R3 ODE |
| M5 | Bare→IR screening Φ₀·(3/14) | TGP MA explicit, ale 3/14 nie e-related (algebraic) |

---

## 5. Miejsca powiązane z badanym zjawiskiem

### 5.1 Source folder (origin)

- **`research/why_n3/`** — origin wszystkich trzech poszlak
  - `PHASE2_n_alpha_derivation.md` (X = e²/4 discovery)
  - `PHASE5_full_propagator.md` (PDG match)
  - `PHASE6_alpha_em_connection.md` (negatywne ustalenie α-em)
  - `r3_phase7_phi0_screening_e2.py` (Φ_eff hint)
  - `tgp_emergent_dirac_propagator.md` (broader context)

### 5.2 Core TGP files (formal definitions)

- **`core/sek00_summary/sek00_summary.tex:74-78`**
  Definicja Φ_eff = Φ₀·(3/14) i bare = 168·Ω_Λ ≈ 115.
  **Kluczowe miejsce dla cosmological hint (P3).**

- **`core/sek00_summary/sek00_summary.tex:60-69`**
  P(φ) = (β/7)φ⁷ - (γ/8)φ⁸ akcja, V(φ) = (γ/3)φ³ - (γ/4)φ⁴ EOM.
  Źródło algebraiczne 3/14 = P(1)/V(1) screening factor.

- **`core/sek10_N0_wyprowadzenie/sek10_N0_wyprowadzenie.tex`**
  Twierdzenie K(φ) = K_geo·φ² na poziomie substratu (vs φ⁴ macroscopic
  TGP_FOUNDATIONS). Niespójność `α dual` (substrate α=1 vs macro α=2).

- **`core/sek08a_akcja_zunifikowana/`**
  K(φ)=φ⁴, β=γ vacuum condition, Φ-EOM. Mass formula bottom-up.

- **`core/sek09_cechowanie/sek09_cechowanie.tex:1050,1077`**
  Brannen lock 2026-05-01: Φ₀ = 24.783 (vs cosmological 24.66).
  **Wpływa na (10/3)·e² match accuracy** — różne lock wartości
  dają różne residua.

- **`core/formalizm/dodatekO_u1_formalizacja.tex:300-403`**
  α_em = Φ₀/(8π·J_phase) z explicit derivation.
  **Dowód że phase sector NIE używa e_Euler** — uzasadnia że
  λ.1 powinien szukać e² **wyłącznie w amplitude sector**.

- **`core/formalizm/dodatekO_u1_formalizacja.tex:405-421`**
  Distinction `J_amp` vs `J_phase`. Klucz dla zakresu λ.1.

- **`TGP_FOUNDATIONS.md:56`**
  K(φ) = K_geo·φ⁴ (macroscopic). Definiuje α=2 dla R3 mass formula.

### 5.3 Powiązane cykle research

- **`research/op-uv-as-ngfp/`** (UV.1)
  AS NGFP fixed point z g* = 0.71, λ* = 0.19, η_N* = -2.
  **Najpoważniejszy kandydat na "γ_φ source"** dla wave-function
  renormalization w R3 (mechanizm M3).

- **`research/op-uv2-mtgp-absolute-scale/`** (UV.2)
  M_TGP scale + dim-less invariants. Może e² łączy się z M_TGP/M_GUT
  ratio przez RG flow integral.

- **`research/op-chi1-newton-constant-derivation/`** (χ.1)
  G_N derivation. Jeśli e² jest cross-sector (amplitude + grav), χ.1
  mogłoby też zawierać e²-trace.

- **`research/op-alpha-fine-structure/`**
  α_em derivation. Negatywny anchor — potwierdza że phase sector NIE
  zawiera e_Euler. **Boundary marker** dla λ.1 scope.

- **`research/op-zeta-mass-spectrum/`** (ζ.1)
  Neutrino mass spectrum + PMNS. Jeśli e² rządzi wszystkimi mass
  ratios w amplitude, neutrino sector też powinien wykazywać hint.
  **Cross-test dla λ.1.**

- **`research/op-phi1-substrate-action-variational/`** (φ.1)
  Substrate action variational + scale-invariance X→λX. Może
  scale-symmetry φ.1 daje natural mechanizm dla e² emergence.

### 5.4 Stałe i wartości używane

- `e ≈ 2.71828` (liczba Eulera, math)
- `e² ≈ 7.38906`
- `e²/4 ≈ 1.84726` (X candidate)
- `e²/2 ≈ 3.69453` (n(2) dla α=2)
- `(10/3)·e² ≈ 24.6302` (Φ_eff candidate)
- `Φ₀ (bare) = 168·Ω_Λ ≈ 115` (sek00:76)
- `Φ_eff (cosmological) = 36·Ω_Λ ≈ 24.66` (sek00:77)
- `Φ_eff (Brannen) = 24.783` (sek09:1077, lock 2026-05-01)
- `Ω_Λ ≈ 0.685` (Planck)
- `3/14` screening factor (P(1)/V(1))
- `8π` w α_em mianowniku (dodatekO:336)

---

## 6. Pytania badawcze

### 6.1 Główne pytania

1. **Czy e²/4 ma derywację analityczną z TGP?**
   Konkretne: pochodzić n(α) z 1-loop (lub non-perturbative) calculation
   w R3 amplitude sector.

2. **Skąd pochodzi (10/3) w Φ_eff = (10/3)·e²?**
   Konkretne: czy 10/3 ma origin w TGP-algebraic ratio (jak 3/14)?

3. **Dlaczego e² w amplitude, nie w phase?**
   Konkretne: jaka strukturalna różnica między J_amp i J_phase wymusza
   że tylko jeden zawiera e_Euler?

4. **Czy e² pojawia się też w innych R3 quantities?**
   Konkretne: cross-test w neutrinach (ζ.1), charged hadrons, BH shadow.

### 6.2 Pomocnicze pytania

- Czy Hobart-Derrick balance point α=4 ma głębszą interpretację?
- Czy n(α) jest dokładnie liniowe, czy ma residual O(α²) korekty?
- Czy istnieje TGP-natural normalizacja jednostek gdzie e²/4 jest
  fundamentalnym coupling?

---

## 7. Ścieżki badawcze (czego próbować, czego NIE)

### 7.1 OBIECUJĄCE — warte eksploracji

1. **Sub-cycle λ.1.RG**: AS NGFP UV.1 + R3 amplitude renormalization
   - Sprawdzić czy γ_φ z anomalous dimension R3 pola amplitude daje
     n(α) = (e²/4)·(4-α)
   - Wymaga: explicit β-function R3 ODE + AS-flow integration

2. **Sub-cycle λ.1.STAT**: Partition function R3 substrate
   - Sprawdzić czy soliton dressing przez fluktuacje substratowe daje
     g₀^(e²/2) factor naturalnie
   - Wymaga: stat-mech derivation gęstości stanów wokół solitonu

3. **Sub-cycle λ.1.CUM**: Iterative dressing g(r+dr) = g(r)·(1+Δ/n)^n
   - User's intuicja "kumulatywne pole" — sprawdzić explicit w R3 ODE
   - Wymaga: continuous limit dyskretnego procesu w r-direction

4. **Sub-cycle λ.1.10/3**: Algebraic search dla 10/3 w TGP
   - Może 10/3 jest related do `5/3·2` (pionowych modes) lub `7/3+1`
     (number-theoretic z TGP integers 7, 8, 12, 14, 56, 36, 168)
   - Wymaga: systematic search w TGP-algebraic ratios

### 7.2 ODRZUCONE — nie wracać do tych ścieżek

1. **R⁵-bridge / 5D Kaluza-Klein** — external import bez ugruntowania
   w TGP-aksjomatach. Patrz `why_n3/PHASE6_Q5_R5_bridge_first_attempt.md`
   (DEPRECATED).

2. **X = α-em w innych jednostkach** — 4 rzędy wielkości różnicy
   (X=1.847 vs α=0.0073), separate sectors w TGP. Patrz
   `why_n3/PHASE6_alpha_em_connection.md` (CLOSED-NEGATIVE).

3. **3/14 screening jako e-related** — pure algebraic ratio P(1)/V(1)
   bez exp/log content. Potwierdzone w sek00:74.

---

## 8. Phase plan (tentative)

### Phase 1 — Structural setup + falsification ranges (~2 tygodnie)

- L1.1 Algebraic search dla 10/3 w TGP integers (number-theoretic)
- L1.2 Explicit β-function dla R3 ODE jako effective theory
- L1.3 Partition function setup dla R3 amplitude fluctuations
- L1.4 Cross-test: czy ζ.1 neutrino spectrum wykazuje e²-hint
- L1.5 Boundary check: dlaczego e² w amplitude, nie w phase (formal proof)

**Score gate:** ≥3/5 PASS → Phase 2 forward

### Phase 2 — Sympy LOCK + analytical derivation (~3-4 tygodnie)

- L2.1 RG flow: Z_φ(μ) z UV.1 AS NGFP, sprawdzić γ_φ = e²/4·(4-α)/4
- L2.2 Cumulative iterative dressing w R3 ODE (user's M4 mechanism)
- L2.3 1-loop scalar self-energy w R3 amplitude sector
- L2.4 Sympy verification: jeśli mechanizm znaleziony, locking exact
- L2.5 Φ_eff = (10/3)·e² pełna derywacja (z Φ₀ bare + screening)

**Score gate:** ≥3/5 PASS → Phase 3 forward

### Phase 3 — Predictions + cross-channel convergence (~2 tygodnie)

- L3.1 R3 mass formula: X promotion EMPIRICAL → DERIVED
- L3.2 Cross-test ζ.1 neutrino: predyktować Σm_ν z e²-mechanism
- L3.3 χ.1 Newton: sprawdzić czy G_N ma e²-trace
- L3.4 Falsyfikacja precyzyjna: PDG α=2 ratio diff <0.001% już matche;
  jeśli e²/4 jest fundamental, drift powinien być structurally limited
- L3.5 Publication anchor: jeśli FULL CONVERGENCE, λ.1 otwiera promotion
  X od `LOCKED EMPIRICAL` do `LOCKED DERIVED` w PREDICTIONS_REGISTRY

**Score gate:** ≥3/5 PASS → λ.1 program END (FULL CONVERGENCE)

---

## 9. Dependencies + prerequisites

### 9.1 Wymagane wcześniejsze cykle (CLOSED)

- ✅ R3 emergent Dirac (`why_n3/PHASE1-5`) — origin poszlak
- ✅ φ.1 substrate-action-variational — scale-invariance setup
- ✅ UV.1 AS NGFP — RG flow infrastructure (jeśli mechanizm M3)
- ✅ R3 mass formula closure (Faza 2 with X = e²/4)

### 9.2 Pomocne ale niewymagane

- UV.2 M_TGP scale (jeśli Φ_eff = (10/3)·e² łączy się z M_TGP)
- χ.1 Newton constant (jeśli e² jest cross-sector)
- ζ.1 neutrino spectrum (cross-test λ.1)

### 9.3 Niezależne (mogą iść równolegle)

- op-alpha-fine-structure (sektor phase, oddzielny)
- op-omega1/2/3 (axion, EM, oddzielne)

---

## 10. Status meta

**Klasyfikacja:** PROPOSED, eksploracyjny.

**Ryzyko:** WYSOKIE — wszystkie trzy poszlaki mają residua 0.001-0.6%,
mogą być statystycznymi coincidences. Brak natural derivation w obecnym
TGP-formalism. Możliwy outcome: zarzucenie po Phase 1 jeśli żadna ścieżka
nie produkuje e² naturally.

**Reward:** WYSOKIE — jeśli e² jest fundamentalna, λ.1 promotuje R3 X
do DERIVED i daje **drugi fundamentalny constant** w TGP po Brannen
K_lep=2/3. To by było **publishable result**: "e² as emergent constant
in TGP amplitude sector".

**Rekomendacja czasowa:** ~2 miesiące Phase 1 + 1-2 dla Phase 2/3 jeśli
Phase 1 daje pozytywny sygnał. Po Phase 1 — **kill criterion**: jeśli
żaden mechanizm M1-M5 nie produkuje e²-trace, zarzucić.

---

## 11. Otwarte pytania meta (dla user / autora TGP)

1. **Czy Φ_eff = 24.66 (cosmological) czy 24.783 (Brannen) jest "the right" wartość?**
   Ten wybór decyduje czy (10/3)·e² match jest 0.12% (cosmological)
   czy 0.62% (Brannen). Te dwa **różne fits** pochodzą z różnych
   strukturalnych anchors.

2. **Czy α (kinetic exponent) = 1 lub 2 (lub coś innego) w TGP-canonical?**
   Niespójność z głównego audytu (sek10 K=φ², sek08a K=φ⁴) wpływa na
   p(α) i derywację e². λ.1 zakłada α=2 (TGP_FOUNDATIONS), ale sektor
   solitonowy preferuje α=1 (sek00:nota o dualizmie α).

3. **Czy istnieje uzasadnienie 10/3 z TGP integers?**
   TGP używa: 3, 7, 8, 12, 14, 36, 56, 168, 24.66, 25, 24.783. Czy
   jakaś kombinacja daje 10/3 = 3.333... naturally?

---

## 12. Końcowa nota

Ten cykl jest **explicit eksploracyjny** — nie ma celu zamykania na siłę.
Jeśli po Phase 1 e² nie pojawi się jako natural emergent constant, lepiej
zarzucić niż forsować numerologiczną pseudo-derywację.

**Najczystszy outcome dla λ.1**: albo **closure** z konkretną mechaniką
(co byłoby przełomem w TGP), albo **honest rejection** z explicit
"e² jest empirical fit, nie fundamental" — bez third option.

---

**Autor proposal:** Sesja R3 + emergent Dirac, 2026-05-01.
**Origin:** `research/why_n3/PHASE2`, `PHASE5`, `PHASE7`.
**Status:** PROPOSED — gotowe do startu Phase 1.
**Następne:** decyzja autora czy uruchomić cykl, lub zostawić jako
otwarte pytanie research.
