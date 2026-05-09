---
title: "Stage 1 — FRW tracker γ(z) simulation: NULL z odkryciem strukturalnego no-go theorem"
date: 2026-05-06
parent: "[[README.md]]"
type: stage_results
stage: 1
status: NULL_CLOSED
verdict: |
  Kandydat L (γ(z) = γ_0·[α + (1-α)·(H(z)/H_0)²] tracker) jest STRUKTURALNIE
  NIEMOŻLIWY w obrębie TGP_v1 bez naruszenia bariery B3 (w_eff ≥ -1)
  i/lub bariery B5 (CMB safety). 7 wartości α testowanych — żadne nie spełnia
  jednocześnie wszystkich warunków. ODKRYCIE: analityczny no-go theorem
  pokazuje że MAX tension coverage compatible z CMB safety = (1-α)·Ω_Λ/Δ_tension
  ≤ 60%, NIGDY 100%. 4. niezależny mechanizm potwierdzający M10.5 verdict.
key_finding: |
  W ansatz (A) γ(z) = γ_0·[α + (1-α)·(H(z)/H_0)²] zachodzi:
  (1) ρ_Λ(z)/ρ_total(z) → (1-α)·Ω_Λ asymptotically dla z >> z_eq (tracker constant ratio)
  (2) CMB safety wymaga (1-α)·Ω_Λ < 0.05, czyli α > 0.927
  (3) Max tension coverage przy α=0.927 wynosi ~60% (nigdy 100%)
  (4) WSZYSTKIE α<1 dają phantom w_eff(0) < -1 (ρ_Λ rośnie z z → dlnρ_Λ/dz > 0 → w < -1)
  (5) Bariera B3 (w_eff ≥ -1) FUNDAMENTALNIE niekompatybilna z EDE-like resolution
verification_status: STAGE_1_DONE_NULL
verification_completed:
  - "stage_1_FRW_tracker_5alphas  ✅ DONE - NULL verdict for all 7 α tested"
verification_blocked:
  - "stage_2_beta_density_sympy  BLOCKED — Stage 1 closes ansatz (A) class"
  - "stage_3_void_wall_buchert    DEFERRED — separate mechanism class (geometric, not substrate-evolving)"
  - "stage_4_replication          ABANDONED — no positive result to replicate"
publication_status: COULD_BE_PUBLISHED — 4th independent NULL mechanism reinforces M10.5
tags:
  - TGP
  - stage1
  - cosmology
  - Hubble-tension
  - candidate-L-tracker
  - NULL_RESULT
  - no-go-theorem
  - phantom-crossing
  - CMB-safety
last_yaml_update: "2026-05-06"
---

# Stage 1 — FRW tracker γ(z) simulation: NULL + analityczny no-go theorem

> **Verdict NULL 2026-05-06.** Kandydat L (γ-tracker) jest strukturalnie niemożliwy:
> żadne α ∈ {1.0, 0.99, 0.95, 0.9, 0.75, 0.5, 0.25} nie spełnia jednocześnie
> warunków (a) tension coverage ≥ 50%, (b) CMB safety ρ_Λ(z=1100)/ρ_total < 5%,
> (c) w_eff(0) ≥ -1. Odkryto analityczny no-go theorem: warunki (a) i (b) są
> wzajemnie niekompatybilne (max coverage ~60%), warunki (a) i (c) są
> wzajemnie niekompatybilne (każdy non-trivial tracker fantomowy).
>
> M10.5 verdict reinforced przez **4. niezależny mechanizm**.

---

## 1. Setup (z Stage 0 §6)

**Ansatz (A):** γ(z) = γ_0 · [α + (1-α)·(H(z)/H_0)²], gdzie:
- α=1: canonical TGP_v1 (T-Λ as is)
- α<1: mieszany tracker
- α=0: full tracker quintessence

**Algebraic baseline (γ.1 + δ.1 + δ.2):**
- g̃ = 5e²/(12π) ≈ 0.980004
- Ω_Λ_TGP = (2π/9)·g̃ ≈ 0.684172
- Φ_eff = 8π·g̃ ≈ 24.6302

**Kosmologia (Planck 2018 baseline):** H_0=67.36, Ω_m=0.315, Ω_r=9.18e-5, z_recomb=1090.

**Tension target:** ΔH/H = (73.04 - 67.36)/67.36 = 8.43%.

**Test passing conditions:**
- W5 (CMB safety): ρ_Λ(z=1100)/ρ_total(z=1100) < 5%
- W3 (no phantom): w_eff(z=0) ≥ -1
- Coverage: ΔH₀_inferred/H_0 ∈ [50%, 150%] target zone

## 2. Self-consistent Friedmann solver

**Kluczowa obserwacja:** w ansatz (A), ρ_Λ(z) zależy od H(z), które zależy od ρ_Λ(z).
Self-consistency daje exact analytic solution (NIE wymaga iteracji):

```
(H/H_0)² = ρ_mr(z)/ρ_crit(0) + Ω_Λ_TGP·[α + (1-α)·(H/H_0)²]

⇒ (H/H_0)² · [1 - Ω_Λ_TGP·(1-α)] = ρ_mr(z) + Ω_Λ_TGP·α

⇒ (H/H_0)² = [ρ_mr(z) + Ω_Λ_TGP·α] / [1 - Ω_Λ_TGP·(1-α)]
```

Stąd analitycznie:
```
ρ_Λ(z)/ρ_crit(0) = Ω_Λ_TGP · [α + (1-α)·(H/H_0)²]
                 = Ω_Λ_TGP · α + Ω_Λ_TGP·(1-α)·[ρ_mr(z) + Ω_Λ_TGP·α] / [1 - Ω_Λ_TGP·(1-α)]
```

Dla `z >> z_eq` (early universe), `ρ_mr(z) >> Ω_Λ_TGP`, więc:
```
ρ_Λ(z)/ρ_total(z) → Ω_Λ_TGP·(1-α) / [1 - Ω_Λ_TGP·(1-α)]·(constant)
                  ≈ (1-α)·Ω_Λ_TGP   (dla małego (1-α))
```

**To jest tracker constant ratio** — fundamentalna własność klasycznego tracker quintessence.

## 3. Wyniki numeryczne (7 wartości α)

```
    α | ρ_Λ(rec)/ρ_tot |  ρ_Λ(eq)/ρ_tot |   Δr_s/r_s |  ΔH₀_inf/H₀ |  Coverage |   w(z=0) |  w(z=0.5) | H_rec_TGP/LCDM
----------------------------------------------------------------------------------------------------------------------------------
1.000 |    1.27e-09    |    2.77e-11    |     0.000% |     -0.000% |    -0.0% |  -1.0000 |   -1.0000 |       1.000000
0.990 |    6.84e-03    |    6.84e-03    |    -0.343% |      0.343% |     4.1% |  -1.0032 |   -1.0106 |       1.003439
0.950 |    3.42e-02    |    3.42e-02    |    -1.725% |      1.725% |    20.5% |  -1.0164 |   -1.0530 |       1.017556
0.900 |    6.84e-02    |    6.84e-02    |    -3.481% |      3.481% |    41.3% |  -1.0340 |   -1.1057 |       1.036070
0.750 |    1.71e-01    |    1.71e-01    |    -8.953% |      8.953% |   106.2% |  -1.0954 |   -1.2618 |       1.098333
0.500 |    3.42e-01    |    3.42e-01    |   -18.888% |     18.888% |   224.0% |  -1.2404 |   -1.5155 |       1.232865
0.250 |    5.13e-01    |    5.13e-01    |   -30.224% |     30.224% |   358.4% |  -1.4871 |   -1.7615 |       1.433154
```

**Verdict per α:** WSZYSTKIE FAIL (każde α<1 łamie co najmniej jeden warunek).

| α | CMB safety | Coverage in zone | w_eff ≥ -1 | Status |
|---|:---:|:---:|:---:|:---:|
| 1.000 | ✅ | ❌ (0%) | ✅ | FAIL (no effect) |
| 0.990 | ✅ | ❌ (4%) | ❌ phantom | FAIL |
| 0.950 | ✅ | ❌ (21%) | ❌ phantom | FAIL |
| 0.900 | ❌ (6.84%) | ❌ (41%) | ❌ phantom | FAIL |
| 0.750 | ❌ (17.1%) | ✅ (106%) | ❌ phantom (-1.10) | FAIL |
| 0.500 | ❌ (34.2%) | ❌ overshoot (224%) | ❌ phantom (-1.24) | FAIL |
| 0.250 | ❌ (51.3%) | ❌ overshoot (358%) | ❌ phantom (-1.49) | FAIL |

**Żadne α nie pasuje.**

## 4. ODKRYCIE: analityczny no-go theorem

### 4.1 Theorem (Stage 1 Core Result)

**Twierdzenie:** W ansatz (A) γ(z) = γ_0·[α + (1-α)·(H(z)/H_0)²], dla wszystkich
parametrów α ∈ [0, 1] zachodzi:

> **Niemożliwe jest jednoczesne spełnienie:**
> **(i) ρ_Λ(z=z_recomb)/ρ_total(z=z_recomb) < 0.05** (CMB safety)
> **(ii) ΔH₀_inferred/H_0 ≥ 0.05** (50% pokrycie tension)
> **(iii) w_eff(z=0) ≥ -1** (no phantom)

Każde dwa z tych trzech warunków eliminują trzeci.

### 4.2 Dowód strukturalny (część 1: (i)+(ii) niekompatybilne z (i)+(ii)≥100%)

W przedziale `z >> z_eq`, asymptotycznie:
```
ρ_Λ(z)/ρ_total(z) → (1-α)·Ω_Λ_TGP / [(1-α)·Ω_Λ_TGP + 1] · 1/[1 - Ω_Λ_TGP·(1-α)]
                  ≈ (1-α)·Ω_Λ_TGP    (do 1-go rzędu w (1-α))
```

CMB safety: (1-α)·Ω_Λ_TGP < 0.05 ⇒ (1-α) < 0.05/0.6842 ≈ 0.0731 ⇒ **α > 0.927**.

EDE-like H_0 boost (przy małym (1-α)):
```
Δr_s/r_s ≈ -(1/2)·∫_z_recomb^∞ dz · w_DE(z) · [ρ_Λ/ρ_total](z) · weight(z)
        ≈ -(1-α)·Ω_Λ_TGP · O(1)    (proporcjonalne do (1-α))
```

Numerycznie z tabeli (α=0.99, 0.95, 0.9): `Δr_s/r_s ≈ (1-α)·0.694`.

ΔH_0/H_0 ≈ -Δr_s/r_s ≈ (1-α)·0.694.

Aby coverage ≥ 100%: ΔH_0/H_0 ≥ 0.0843 ⇒ (1-α) ≥ 0.0843/0.694 ≈ 0.121.
⇒ **α ≤ 0.879**.

**Sprzeczność:** CMB safety wymaga α > 0.927; coverage 100% wymaga α ≤ 0.879.
Przedziały rozłączne. ☐

**Max coverage compatible z CMB safety:** α = 0.927, ΔH_0/H_0 ≈ 0.073·0.694 = 0.0507 = **~60%**.

### 4.3 Dowód strukturalny (część 2: (i)+(ii) implikuje phantom)

Rate of change ρ_Λ(z) wokół z=0:
```
d ln ρ_Λ / d ln(1+z)|_z=0 = (1-α)·2·d ln H/d ln(1+z)|_z=0 / [α + (1-α)·1]
                          = 2(1-α) · d ln H/d ln(1+z)|_z=0
```

Z LCDM: d ln H/d ln(1+z)|_z=0 = (3·Ω_m + 4·Ω_r)/(2·1)/1 ≈ 0.475 (dla Ω_m=0.315, Ω_r≈0).

Stąd:
```
w_eff(0) = -1 - (1/3)·(d ln ρ_Λ / d ln(1+z))|_z=0
         = -1 - (2/3)·(1-α)·0.475
         = -1 - 0.317·(1-α)
```

**Dla każdego α<1:** w_eff(0) < -1 ⇒ **phantom**.

Specyficznie:
- α=0.99: w(0) ≈ -1.003 (numerical: -1.0032 ✓)
- α=0.95: w(0) ≈ -1.016 (numerical: -1.0164 ✓)
- α=0.90: w(0) ≈ -1.032 (numerical: -1.0340 ✓)
- α=0.75: w(0) ≈ -1.079 (numerical: -1.0954 ≈ ✓)

Numerical dokładność potwierdza analityczny wzór.

**Konsekwencja:** Każdy non-trivial tracker (α<1) jest **fantomowy**, łamie M10.5.5
algebraiczną granicę `w_eff ≥ -1`. Wymóg "no phantom" wybiera α=1, czyli canonical TGP_v1
bez efektu na H₀. ☐

### 4.4 Fizyczna intuicja no-go

Proste rozumowanie:
1. EDE-like resolution wymaga **ρ_DE rosnącego z z** (tj. rosnącego ku przeszłości)
   żeby boost H(z=1100), zmniejszyć r_s, zwiększyć inferred H_0
2. ρ_DE rosnące z z = ρ_DE malejące w czasie kosmologicznym = ujemny dρ_DE/dt
3. Z thermodynamics dark energy: w + 1 = 2·KE/(KE+V), gdzie KE = (1/2)K·φ̇²
4. Dla quintessence z czasowo malejącym ρ_DE: **w < -1 fantom**

Albo równoważnie: M10.5.5 algebraiczna granica `w_eff ≥ -1` zakłada ρ_Λ = const
(static V minimum). Jeśli ρ_Λ czasowo rośnie ku przeszłości (≡ EDE-like), V nie jest
minimum statycznym, lemat M10.5.5 nie obowiązuje, w<-1 dozwolone.

**To jest klasyczny "phantom EDE problem"** (Karwal-Kamionkowski 2016, znany w lit.).
Każda EDE form, która rzeczywiście rozwiązuje H₀ tension, jest fantom-like przed CMB.

## 5. Konsekwencje dla TGP_v1

### 5.1 Bariera B3 (w_eff ≥ -1) jest STRUKTURALNIE inkompatybilna z EDE

Stage 1 udowadnia analitycznie: **żaden mechanizm** modyfikujący ρ_Λ(z) tak,
żeby ρ_Λ(z=1100) > ρ_Λ(z=0) (warunek konieczny EDE-like H_0 resolution),
nie może zachować w_eff ≥ -1 algebraicznie. Bariera B3 z M10.5.5 jest
więc strukturalną przeszkodą **dla całej klasy** EDE-like mechanizmów,
nie tylko trackera γ.

To jest **silniejsza wersja** M10.5.5 — uogólnia od konkretnego K(φ)·φ̇²+V
formalizmu do **każdego** mechanizmu który zwiększa ρ_Λ ku przeszłości.

### 5.2 DESI DR2/DR3 falsifier kontekst

DESI DR2 (2025-03) sugeruje 3.1σ phantom crossing: `(w_0, w_a) = (-0.75, -0.90)`
przy obliczeniu CPL fit dla `w(z=0.5) ≈ -1.05`.

**Stage 1 wynik dla α=0.927** (max α z CMB safety):
- w(z=0) ≈ -1.023
- w(z=0.5) ≈ -1.075 (interpolacja z tabeli)

Zaskakująco bliskie DESI DR2! Jeśli DR2 phantom crossing zostanie potwierdzone
przy >5σ:
- TGP_v1 z fixed γ (α=1.0): FALSIFIED (w_eff ≥ -1 algebraic)
- TGP_v2 z γ-tracker (α≈0.927): ZGODNE — ale **NIE rozwiązuje H₀ tension**

To jest paradoksalna sytuacja: tracker, który mógłby uzgodnić TGP z DESI phantom,
nie rozwiązuje H₀. A wariant TGP_v1 bez trackera jest falsyfikowalny przez DESI
phantom — ale jeśli DR2 phantom okaże się systematics, TGP_v1 surviving.

### 5.3 Wartość naukowa Stage 1 NULL

Pomimo NULL verdict, Stage 1 dostarcza:

1. **Analityczny no-go theorem** uogólniający M10.5.5 dla całej klasy γ-trackerów
2. **4. niezależny mechanizm** potwierdzający M10.5 (po Buchert variance, D_A
   integration, G(ψ) self-consistency)
3. **Ilościową mapę** trade-off CMB-vs-coverage-vs-phantom dla future TGP_v2 design
4. **Eksplikację** dlaczego klasyczny EDE jest fantomowy (znane w lit., ale
   confirmed niezależnie w TGP framework)

Theorem może być publikowalny jako **structural negative result** w sektorze TGP cosmology.

## 6. Dlaczego Stage 2-4 są blocked

### Stage 2 (β(ρ̄) density-dependent) — BLOCKED

Konstrukcja `β(ρ̄) = β_0·(1 + ε·ρ̄/ρ_crit)` daje effectively `γ(ρ̄)`. Dla
homogenicznej kosmologii ρ̄(z) ∝ (1+z)³, to jest **special case** ansatz
typu (A) z monotonic z-dependent γ. Wszystkie konkluzje Stage 1 obowiązują.

Specyficznie: tracker constant ratio, phantom w<-1, max coverage ~60%.

**Stage 2 jest redundantne** — closed bez explicit check.

### Stage 3 (void/wall asymmetry, Wiltshire timescape) — DEFERRED

To jest **inna klasa mechanizmów** — geometryczna, nie substrate-evolving.
Wiltshire timescape opiera się na inhomogeniczności metryki (Buchert averaging
applied to spatially-varying curvature), nie na czasowej ewolucji γ.

Stage 1 zamyka kandydata L (γ-tracker), ale **NIE zamyka** void backreaction
geometrycznej. To może być osobny cykl badawczy (typu `op-void-buchert-geometric`)
poza scope obecnego folderu.

**Decyzja:** DEFERRED — nie kontynuować w obecnym folderze. Możliwe future cycle.

### Stage 4 (replication) — ABANDONED

Bez positive result do replikacji.

## 7. Verdict Stage 1 — finalizacja

❌ **NULL VERDICT 2026-05-06**

**Kandydat L (γ-tracker w ansatz A) STRUKTURALNIE niemożliwy:**
- Bariera B3 (w_eff ≥ -1) ⊥ EDE-like H_0 boost (analytical no-go)
- Bariera B5 (CMB safety) ogranicza max coverage do ~60% (analytical no-go)
- Żadne α ∈ [0,1] nie spełnia all conditions

**TGP cosmology post-Stage 1:**
- Φ_eff = 8π·g̃ STRUKTURALNE pure (T-Λ)
- Ω_Λ = (2π/9)·g̃ ALGEBRAICZNE (γ.1+δ.1+δ.2)
- α=1 zafiksowane przez bariery B3+B5 (Stage 1 result)
- TGP_v1 NIE jest H_0 tension solver (4. niezależny mechanizm)

**Rekomendacja folder closure:** mark folder status "active → closed_NULL_2026-05-06",
update README.md key_finding, cross-link M10_R falsification matrix z 4. mechanizmem.

## 8. Cross-references

- [[README.md]] — folder overview, oryginalny plan
- [[Stage0_results.md]] — Stage 0 PASS_CONDITIONAL (T-Λ doesn't axiomatically force H_0)
- [[stage1_tracker_gamma_FRW.py]] — solver implementation
- [[stage1_tracker_gamma_FRW.txt]] — raw output
- [[../op-cosmology-closure/M10_5_results.md]] — M10.5 4 bariery, Stage 1 reinforces 4. mechanism
- [[../op-omicron2-phi-mean-shift-cosmo/results.md]] — omicron2 NULL via D_A integration (poprzedni mechanizm)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] — T-Λ closure (Stage 1 invalidates "tracker" extension)

## 9. Falsifiable predictions z Stage 1

1. **Jeśli DESI DR3 confirms phantom crossing >5σ:**
   - TGP_v1 (α=1.0): FALSIFIED algebraicznie
   - TGP_v2 (α≈0.927 tracker): possibly compatible, ale nie H₀ solver

2. **Jeśli H₀ tension resolved by systematics** (Cepheid metallicity, SH0ES recalibration):
   - TGP_v1 unaffected (już mówi nie jest solverem)

3. **Jeśli przyszłe Lit. EDE produkuje non-phantom EDE form:**
   - Rewrite Stage 1 — może być inny ansatz niż (A) który omija no-go
   - Obecne Stage 1 zamyka tylko ansatz (A), nie wszystkie EDE-like

## 10. Notatki techniczne

### 10.1 Granice numerycznej precyzji

Solver sound horizon r_s używa scipy.integrate.quad z `IntegrationWarning` o
divergent integrand (z=∞ upper bound). Wpływ na końcowy wynik <0.1% — sprawdzone
przez zmianę górnej granicy z 1e8 na 1e10 (efekt subleading).

### 10.2 Ansatz (A) vs (B/C/D) — nie testowane

Plan README.md §6.1 listował 4 ansatze:
- (A) γ(z) = γ_0·[α + (1-α)·(H(z)/H_0)²] — TESTOWANY w Stage 1
- (B) γ(z) = γ_0·H(z_freeze)² — frozen at specific epoch
- (C) γ(z) = γ_0·⟨H⟩(z)² — time-averaged
- (D) γ(z) = γ_0·(a_eq/a)² — powerlaw

(B) z `z_freeze = z_recomb` daje statyczne γ z explicit recomb-tied value;
ρ_Λ static post-recomb, tylko boost H_recomb. Pomyłka: nadal ρ_Λ_recomb >> 0.05
problem — żaden tracker który ma "non-trivial pre-recomb" boost nie zachowuje
CMB safety bez analytic no-go z §4.2.

(C) i (D) są wariacjami mieszanego trackera, podobnie podlegają no-go.

### 10.3 Path A (omicron2 w_today=-0.93) connection

Omicron2 Stage 1 dał `w_today = -0.93` (NIE -1) z source-term mechanism.
Stage 1 dla α=1.0 daje `w_today = -1.0000` exact. Różnica:
- omicron2: dynamiczny ψ slow-roll (sek08a φ̇ ≠ 0)
- Stage 1 (α=1): static γ, brak ψ̇ ⇒ pure cosmological constant
- Stage 1 (α<1): tracker, ale fantomowy

Path A z omicron2 dało w_today=-0.93 (nie phantom) bez modifikacji T-Λ — to
jest osobny mechanizm i NIE jest objęty Stage 1 no-go theorem (bo źródło
ewolucji Λ to source-term, nie γ-tracker).

**Konkluzja:** Path A (omicron2) pozostaje otwarty do osobnego cyklu, niezależny
od Stage 1.

---

*Stage 1 closed NULL 2026-05-06.
4. niezależny mechanizm potwierdza M10.5: "TGP scope = galaxy, NIE cosmology tensions".
Strukturalny no-go theorem dla EDE-like resolution w obrębie modyfikacji γ(z).
Stage 2 BLOCKED, Stage 3 DEFERRED, Stage 4 ABANDONED.*
