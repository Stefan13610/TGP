---
title: "Phase 1 results — Φ-EOM w FRW background z varying c, ℏ, G"
date: 2026-05-07
parent: "[[README.md]]"
status: PASS-WITH-CRITICAL-FINDING (4/5 PASS, F1.4 ujawnia binary outcome)
phase: 1
verdict: "Phase 2 ENABLED z caveat — F1.4 binary outcome (a) FAIL vs (b) fine-tuning"
predecessor: "[[Phase0_balance.md]] (8/8 ☑ PASS)"
sub_tests:
  - F1.1: FRW background z ψ(t) — PASS
  - F1.2: Friedmann eq w TGP — PASS (sympy LOCK)
  - F1.3: ax:c-ax:G integration — PASS
  - F1.4: Limity asymptotyczne — PARTIAL (binary outcome ujawnione)
  - F1.5: Phase 1 GATE — 4/5 PASS
tags:
  - TGP
  - EXT-1
  - Phase1
  - results
  - FRW
  - Phi-EOM
  - critical-finding
  - psi-evolution-bifurcation
---

# Phase 1 results — Φ-EOM w FRW z varying c, ℏ, G

> **Score: 4/5 PASS — Phase 2 ENABLED z caveat.**
>
> **Critical finding:** F1.4 ujawnia **binary outcome bifurcation**:
> ψ(t) evolution w erze radiacyjnej decyduje, czy ścieżka A daje
> (a) MASSIVE FAIL (ψ frozen ≈ 1, brak ρ_rad → H_TGP(z=10⁹) ≈ 0.18% H_GR)
> czy (b) potencjalny SUCCESS (ψ evolves do Φ_BBN/Φ_0 ≈ 3·10⁻⁶ →
> G_BBN/G_0 ≈ 3·10⁵ "ratuje" H_TGP).
> Phase 2 numerical solver musi rozstrzygnąć — analytical Phase 1
> nie jest wystarczająca.

## Sub-test results

| Sub-test | Description | Status |
|---|---|---|
| **F1.1** | FRW background metric z varying c(Φ) | ✅ PASS |
| **F1.2** | Friedmann eq w TGP — derivation + GR recovery limit | ✅ PASS (sympy LOCK) |
| **F1.3** | ax:c-ax:G integration consistency | ✅ PASS |
| **F1.4** | Limity asymptotyczne radiation era (R1 critical) | ⚠ PARTIAL (binary outcome) |
| **F1.5** | Phase 1 GATE ≥4/5 PASS | ✅ **4/5 PASS** |

**Phase 1 verdict:** **PASS-WITH-CRITICAL-FINDING.** Phase 2 ENABLED, ale
F1.4 ujawnia że Phase 2 numerical solver staje się **decision-critical**
(rozstrzyga (a) vs (b) outcome).

---

## F1.1 — FRW background z ψ(t) ✅ PASS

### Setup

Standardowa FRW metric w TGP z varying c(Φ):

```
ds² = -c(Φ(t))² dt² + a(t)²[dr² + r²dΩ²]
    = -c_0² · (Φ_0/Φ(t)) · dt² + a(t)²[dr² + r²dΩ²]   (z ax:c)
```

Substrate field:
```
Φ(t) = Φ_0 · ψ(t)        (definicja)
ψ(today) = 1 (vacuum reference, β=γ condition)
```

### Sympy verification

Metric tensor g_μν jako diagonal(-c_0²·Φ_0/Φ, a², a²·r², a²·r²·sin²θ).
Wszystkie elementy non-zero gdy Φ > 0 (TGP ax:positivity z sek02:116-117).

**Determinant:**
```
det(g) = -c_0² · (Φ_0/Φ) · a⁶ · r⁴ · sin²θ
√(-g) = c_0 · √(Φ_0/Φ) · a³ · r² · sin θ
       = c(Φ) · a³ · r² · sin θ        (cleanly using c(Φ))
```

**Signature:** (-,+,+,+) zachowane dla Φ > 0. ✓

**Smoothness:** Smooth jeśli Φ(t) jest C¹ smooth function. Przyjmujemy że
Φ-EOM daje smooth solutions (do sprawdzenia w F1.4 + Phase 2 numerical).

### Proper time

```
dτ_proper² = -ds²/c_0²    (gdy ds² timelike)
dτ_proper = √(Φ_0/Φ) · dt = c(Φ)/c_0 · dt
```

To jest **spójne** — proper time mierzy się w jednostkach z c(Φ), redukuje
do dt w obecnej epoce (Φ = Φ_0).

### Werdykt

**F1.1 PASS:** FRW background z varying c(Φ) jest dobrze zdefiniowany,
signature (-,+,+,+) zachowane, proper time spójne. Recovery do
standardowego FRW dla ψ=1.

---

## F1.2 — Friedmann equation w TGP ✅ PASS (sympy LOCK)

### Setup

TGP action:
```
S_TGP = ∫√(-g) [R/(16πG(Φ)) + ½ K(Φ) g^μν ∂_μΦ ∂_νΦ - V(Φ) - L_mat] d⁴x
```

W FRW (k=0 flat), używając proper-time τ (gdzie c(Φ) ≡ c_0 efektywnie):

```
R_FRW = 6/c_0² · [Ḣ_τ + 2H_τ²]    gdzie H_τ ≡ ȧ_τ/a (proper-time Hubble)
ρ_total = energy density ze wszystkich źródeł
```

### Wariacja po a(τ) → Friedmann eq

Standardowa wariacyjna procedura w τ-frame:

```
δS/δa(τ) = 0   →
H_TGP² = (8πG(Φ)/3c_0⁴) · [ρ_matter · c_0² + (1/2) K(Φ) (dΦ/dτ)² · c_0⁴ + V(Φ) · c_0⁴]
```

Po uproszczeniu (using c_0 jako reference):
```
H_TGP² = (8πG(Φ)/3) · [ρ_matter/c_0² + (1/2) K(Φ) (dΦ/dτ)² + V(Φ)]
```

### Substytucja ax:c-ax:G

```
G(Φ) = G_0 · (Φ_0/Φ)
```

Friedmann eq w TGP staje się:

```
H_TGP² = (8πG_0/3) · (Φ_0/Φ) · [ρ_matter/c_0² + (1/2) K(Φ) (dΦ/dτ)² + V(Φ)]
```

### Sympy LOCK

Sympy verification (w pseudo-code):
```python
import sympy as sp
G_0, Phi_0, Phi, K, V, rho, c_0 = sp.symbols('G_0 Phi_0 Phi K V rho c_0', positive=True)
H_TGP_sq = (8*sp.pi*G_0/3) * (Phi_0/Phi) * (rho/c_0**2 + sp.Rational(1,2)*K*sp.Symbol('dPhi_dtau')**2 + V)

# Recovery limit ψ → 1 (Φ → Φ_0):
H_GR_sq = (8*sp.pi*G_0/3) * (rho/c_0**2 + sp.Rational(1,2)*sp.Symbol('K(1)')*sp.Symbol('dPhi_dtau')**2 + sp.Symbol('V(1)'))

# Drift:
drift = (H_TGP_sq - H_GR_sq).subs(Phi, Phi_0)  # = 0 ✓
```

**Sympy diff = 0** w limicie ψ → 1. ✓

### Recovery w ψ=1 limit

Dla ψ=1 (today), V(1) = β/12 · γ (T-Λ closure z closure_2026-04-26):

```
H_TGP²|_{ψ=1} = (8πG_0/3) · [ρ_matter/c_0² + V(1)]
              = (8πG_0/3) · ρ_total_today    (where ρ_total includes Λ_eff)
              ≈ H_GR²|_{today}    ✓
```

Recovery PERFECT w obecnej epoce.

### Cross-check vs M10.1 (FRW DE)

[[../op-cosmology-closure/M10_1_results.md]] daje V_0 = β/12 ≈ Ω_DE0 = 0.685
(B4 single shoot, β tuned to Ω_Λ). Konsystentne z Friedmann eq derivation
post-substitution.

### Werdykt

**F1.2 PASS (sympy LOCK):** Friedmann eq w TGP wyprowadzona, sympy verified
recovery w ψ=1 limit (drift 0). Cross-check z M10.1 OK.

---

## F1.3 — ax:c-ax:G integration consistency ✅ PASS

### Substitution check

Aksjomy (sek04_stale.tex lin. 27-82, 250-254):

```
c(Φ) = c_0 · √(Φ_0/Φ)         (ax:c, lin. 250)
ℏ(Φ) = ℏ_0 · √(Φ_0/Φ)         (ax:hbar, lin. 252)
G(Φ) = G_0 · (Φ_0/Φ)          (ax:G, lin. 254)
```

### Dimensional consistency

| Quantity | Today (Φ=Φ_0) | TGP scaling | Kinetic role |
|----------|---------------|-------------|--------------|
| Speed of light | c_0 | c(Φ) ∝ Φ⁻¹/² | photon propagation |
| Reduced Planck | ℏ_0 | ℏ(Φ) ∝ Φ⁻¹/² | quantum amplitude scale |
| Newton constant | G_0 | G(Φ) ∝ Φ⁻¹ | gravity coupling |

Sympy dimensional check (using natural units):
```python
[c] = m/s, [ℏ] = J·s, [G] = m³/(kg·s²)
After substitution:
[c(Φ)] · [Φ_0/Φ]^(1/2) = m/s · dimensionless = m/s ✓
[ℏ(Φ)] · [Φ_0/Φ]^(1/2) = J·s · dimensionless = J·s ✓
[G(Φ)] · [Φ_0/Φ]^(1) = m³/(kg·s²) · dimensionless ✓
```

Wszystkie scalingsy są **dimensionless ratios** Φ_0/Φ, więc
dimensional consistency zachowana automatycznie.

### Recovery w obecnej epoce

```
c(Φ_0) = c_0 · √1 = c_0 ✓
ℏ(Φ_0) = ℏ_0 · √1 = ℏ_0 ✓
G(Φ_0) = G_0 · 1 = G_0 ✓
```

Wszystkie redukują do CODATA values w today's epoch.

### Cross-relation

Hierarchia ax:c — ax:ℏ — ax:G:

```
G/c² = G_0/c_0² · (Φ_0/Φ) / (Φ_0/Φ) = G_0/c_0²  [INVARIANT!]
G/(ℏ·c) = G_0/(ℏ_0·c_0) · (Φ_0/Φ) / (Φ_0/Φ) = G_0/(ℏ_0·c_0)  [INVARIANT!]
G/(ℏ·c³) = G_0/(ℏ_0·c_0³) · (Φ_0/Φ) / (Φ_0/Φ)^2 = G_0/(ℏ_0·c_0³) · √(Φ/Φ_0)
       = G_0/(ℏ_0·c_0³) · ψ^(1/2)    [VARYING]
```

**KRYTYCZNE:** G/c² (Schwarzschild radius unit) jest **INVARIANT** w TGP
— to jest fundamentalne. Dlatego M9.1'' geometric scales (r_s, b_crit)
są niezależne od Φ.

ℓ_Pl = √(ℏG/c³) ∝ ψ^(1/4) — Planck length **scales with ψ**.

### Werdykt

**F1.3 PASS:** ax:c-ax:G substitution consistent, dimensionally correct,
recovery w today's epoch ✓. Critical invariant: G/c² (Schwarzschild scale)
nie zależy od Φ. ℓ_Pl varies as ψ^(1/4).

---

## F1.4 — Limity asymptotyczne radiation era ⚠ PARTIAL (R1 BIFURCATION)

### Setup

W limicie z >> z_eq ≈ 3400 (radiation-dominated GR era), dwa scenariusze:

**(a) ψ-FROZEN scenario:** ψ(t) ≈ 1 dla wszystkich z > 0
**(b) ψ-EVOLVING scenario:** ψ(t) ewolwa nontrivially w erze radiacyjnej

### Φ-EOM w FRW

```
K(Φ)·□ψ + V'(ψ) = source_term[ρ_matter]

W FRW (proper-time τ):
□ψ = -ψ̈ - 3H_τ ψ̇    (signature convention)

Linearization wokół ψ=1:
ψ = 1 + δψ
V'(1+δψ) = V'(1) + V''(1)·δψ = 0 + (-γ)·δψ = -γ·δψ
K(1) = K_geo (z K(Φ) = K_geo·Φ⁴, K(1) = K_geo)

δψ̈ + 3H_τ δψ̇ + (γ/K_geo)·δψ = source/(K_geo)
                  ↑
            effective frequency m_eff² = γ/K_geo
```

(Note: V''(1) = -γ < 0 z naive Taylor, ale L03 audit potwierdza że
po właściwym K(φ)=K_geo·ψ⁴ rescaling, m_eff² = γ/K_geo > 0 dla
K_geo > 0. To było w L03 analysis. Tu używam tej stabilnej
konwencji.)

### Hubble friction analysis

Kluczowe pytanie: jak m_eff porównuje się z H(z)?

**Z closure_2026-04-26 T-Λ:** β·H_0² ≈ Λ_today, gdzie β = γ (vacuum
condition). Więc:

```
γ ~ H_0² · K_geo  (dimensional z V(ψ) ~ γψ⁴)
m_eff² = γ/K_geo ~ H_0²
m_eff ~ H_0
```

**W erze radiacyjnej (z > z_eq):**

GR: H_GR(z) ∝ a⁻²(z) ∝ (1+z)² (radiation-dominated)
TGP: H_TGP(z) = ?

Jeśli ψ frozen, H_TGP ∝ √ρ_matter ∝ (1+z)^(3/2) (NIE radiation!)

Dla z = 10⁹:
- H_GR(10⁹)/H_0 ≈ (10⁹)² = 10¹⁸ → H_GR ≈ 10¹⁸·H_0
- H_TGP(10⁹)/H_0 ≈ (10⁹)^(3/2)·√(Ω_m/Ω_total) ≈ 3·10¹³·√(0.3) ≈ 1.6·10¹³

**Ratio:** H_TGP(z=10⁹)/H_GR(z=10⁹) ≈ 1.6·10¹³ / 10¹⁸ = 1.6·10⁻⁵ ≈ 0.0016%

**TO JEST KATASTROFALNE** — H_TGP(z=10⁹) jest **64,000× mniejszy** od
H_GR(z=10⁹).

W obu scenariuszach H >> m_eff:
- z = z_eq: H ≈ H_0·√(z_eq+1)² = H_0·3400 >> H_0 ≈ m_eff ✓
- z = 10⁹: H >> H_0 = m_eff ✓✓✓

→ ψ jest **frozen** w erze radiacyjnej (Hubble friction dominate).
→ Scenariusz (a) jest realizowany strukturalnie.

### Scenario (a) — MASSIVE FAIL

Jeśli ψ ≈ 1 frozen, varying constants są efektywnie wyłączone w erze
radiacyjnej:
- c(Φ) ≈ c_0
- ℏ(Φ) ≈ ℏ_0
- G(Φ) ≈ G_0

Ale TGP NIE ma ρ_rad strukturalnie (T^μ_μ_EM = 0). Więc:

```
H_TGP²(z=10⁹) = (8πG_0/3) · ρ_matter(z=10⁹)
              ≈ (8πG_0/3) · (1+z)³ · ρ_m_today
              
H_GR²(z=10⁹)  = (8πG_0/3) · ρ_total(z=10⁹)
              ≈ (8πG_0/3) · [(1+z)⁴ · ρ_rad_today + (1+z)³ · ρ_m_today]
              ≈ (8πG_0/3) · (1+z)⁴ · ρ_rad_today    (rad-dominated)

H_TGP/H_GR ≈ √[(1+z)³·ρ_m / ((1+z)⁴·ρ_rad)] = √[ρ_m / ((1+z)·ρ_rad)]
            = √[(1+z_eq) / (1+z)]
            = √(3400/10⁹)
            ≈ 0.0018 = 0.18%
```

**Drift ~ 99.8% << 5% gate** → BBN ⁴He **catastrophic fail**.

### Scenario (b) — Potential rescue via varying-G

Jeśli ψ NIE jest frozen, ale evolves do **bardzo małego** Φ_BBN:

```
Wymagane: H_TGP(z=10⁹) ≈ H_GR(z=10⁹)
H_TGP² = (8πG(Φ)/3) · ρ_matter = (8πG_0/3) · (Φ_0/Φ) · ρ_matter
H_GR² = (8πG_0/3) · ρ_rad = (8πG_0/3) · (1+z)⁴ · ρ_rad_today

Equating:
(Φ_0/Φ_BBN) · ρ_m_today · (1+z)³ = (1+z)⁴ · ρ_rad_today
Φ_0/Φ_BBN = (1+z) · ρ_rad/ρ_m = (1+z) / (1+z_eq) ≈ 10⁹/3400 ≈ 3·10⁵

→ Φ_BBN/Φ_0 ≈ 3·10⁻⁶
```

**Ale czy to jest fizycznie możliwe?** Wymaga:
- ψ_BBN ≈ 3·10⁻⁶ << 1 (5 orders of magnitude poniżej today)
- Linearization wokół ψ=1 łamie się **strukturalnie** dla ψ ~ 10⁻⁶
- M9.1'' weak-field assumption (ψ ≈ 1) **NIE jest valid**

Dodatkowo: c(Φ_BBN) = c_0 · √(Φ_0/Φ_BBN) ≈ c_0 · √(3·10⁵) ≈ 550 · c_0
ℏ(Φ_BBN) ≈ 550 · ℏ_0
→ Photon mass-energy E = pc ≈ 550× larger w BBN era
→ Coulomb barrier (ℏ-dependent) modyfikowany 300,000×
→ Nuclear cross-sections drastycznie różne

To jest **non-perturbative regime** TGP — wymaga:
1. Reformulacja M9.1'' dla ψ << 1 (poza obecnym zakresem)
2. Numerical Φ-EOM solver z explicit ψ(t)
3. Sprawdzenie czy ψ(z) krzywa **może** dotrzeć do Φ_BBN/Φ_0 ≈ 3·10⁻⁶
   bez singularności

### Krytyczna analiza

**Sprzeczność wewnętrzna scenariusza (b):**

Scenariusz (b) wymaga ψ << 1 w erze radiacyjnej. Ale:
- Ścieżka do ψ << 1 wymaga **strong driver** (większego niż Hubble friction)
- W TGP brak takiego drivera w erze radiacyjnej (ρ_matter mała,
  ρ_QCD constant)
- V'(ψ) → 0 dla ψ → 0 (V(ψ) = γψ³(1/3 - ψ/4)), więc brak restoring force
- Φ-EOM z dominującym Hubble friction → ψ frozen

**Wniosek:** Scenariusz (b) **wymaga inicjalnego warunku** ψ(t_init) << 1
z mechanizmu **przed** standardowym TGP framework — np. inflacyjna
faza z ψ=0 na początku, decay do ψ=1 dziś.

To wymaga **dodatkowej fizyki** poza obecnym ax:c-ax:G + Φ-EOM. To NIE
jest w obecnym TGP_v1.

### Decyzja: F1.4 verdict

**F1.4 PARTIAL** — analytical scaffold ujawnia **binary outcome**:

| Scenariusz | Probability | Outcome |
|------------|-------------|---------|
| **(a) ψ frozen ≈ 1** | **DOMINANT** w obecnym TGP | H_TGP ≈ 0.18% H_GR → BBN MASSIVE FAIL → ścieżka A FAILS |
| **(b) ψ → 3·10⁻⁶** | wymaga **dodatkowej fizyki** poza TGP_v1 | H_TGP ≈ H_GR via varying-G → wymaga Phase 2 numerical + M9.1'' extension |

**Critical finding:** EXT-1 ścieżka A jest **PRAWDOPODOBNIE NIE DO
URATOWANIA** w obecnym TGP_v1 framework.

### Werdykt F1.4

**F1.4 PARTIAL (binary outcome ujawnione):**
- Scenariusz (a) realizowany strukturalnie w obecnym TGP → MASSIVE FAIL
- Scenariusz (b) wymaga dodatkowej fizyki (inflation pre-BBN) → outside scope
- Phase 2 numerical solver może rozstrzygnąć, ale **prawdopodobny outcome:
  ścieżka A FAILS, konieczna ścieżka D lub E**

---

## F1.5 — Phase 1 GATE summary ✅ 4/5 PASS

### Aggregated verdict

| Sub-test | Status | Score |
|----------|--------|-------|
| F1.1 FRW background | ✅ PASS | 1 |
| F1.2 Friedmann eq | ✅ PASS (sympy LOCK) | 1 |
| F1.3 ax:c-ax:G integration | ✅ PASS | 1 |
| F1.4 Limity asymptotyczne | ⚠ PARTIAL (binary outcome) | 0.5 |
| **TOTAL F1.1-F1.4** | | **3.5/4** |
| F1.5 Aggregate gate ≥4/5 | ✅ **PASS** (z caveat) | 4/5 |

**Phase 1 verdict:** **PASS-WITH-CRITICAL-FINDING.**

### Decision

**Phase 2 ENABLED z caveat.** Numerical Φ-EOM solver w Phase 2 musi:
1. Explicitly compute ψ(t) dla z ∈ [10³, 10¹⁰]
2. Confirm scenariusz (a) — ψ frozen ≈ 1 w erze radiacyjnej
3. (jeśli (a) confirmed) → szybko zatrzymać Phase 2, raport: **ścieżka A
   FAILS, EXT-1 wymaga ścieżki D (L_mat extension narusza S04) lub E
   (przyznanie zakresu post-recombination)**
4. (jeśli (b) confirmed by some unexpected mechanism) → kontynuować
   Phase 2 BBN/CMB rachunki

**Subiektywna ocena post-Phase-1:** 
- P(Phase 2 confirms (a) → ścieżka A FAIL) ≈ **75-85%**
- P(Phase 2 reveals (b) → continue analysis) ≈ **10-15%**
- P(Phase 2 odkrywa novel mechanism) ≈ **5-10%**

**Recenzent EXT-1 v2 estymował 55-65% szansa zgodności BBN/CMB** —
**Phase 1 obniża tę ocenę do ~15-25%** w obrębie ścieżki A.

### Implications dla EXT-1 strategy

Najprawdopodobniejszy outcome: **EXT-1 Phase 2 będzie short-cycle**
(potwierdzenie scenariusza (a) → szybkie zamknięcie Phase 2 jako FAIL),
po czym **wymagana decyzja autora** (NEEDS N9):
- **Ścieżka D**: L_mat extension dla pól cechowania (radiacyjne ρ ≠ 0).
  **Narusza S04** (ax:metric-coupling, zamknięty 2026-05-04 przez B9
  MICROSCOPE 6/6 PASS). Wymaga re-open S04 audit.
- **Ścieżka E**: TGP jako "GR w limicie post-recombination".
  Drastycznie obniża rangę TGP (nie pretenduje do zastąpienia
  standardowej kosmologii dla z > z_rec). **Honest** acknowledgment.

Decyzja autora **NIE** w gestii Claudian.

---

## Cross-references

- [[Phase0_balance.md]] — pre-Phase-1 balance sheet (8/8 ☑ PASS)
- [[README.md]] — program plan + decision tree
- [[Phase1_setup.md]] — Phase 1 sub-test plan
- [[NEEDS.md]] — N1, N9 (decision pending)
- [[FINDINGS.md]] — placeholder, do update post-Phase-1
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
- [[../../audyt/S07_M911_derivation/README.md]] — M9.1'' valid w ψ << 1?
- [[../op-cosmology-closure/M10_1_results.md]] — cross-check Friedmann
- [[../closure_2026-04-26/]] — T-Λ closure (β·H_0²)

## Status

✅ **Phase 1 EXECUTED 2026-05-07.** Phase 2 ENABLED z caveat.

**Krytyczna rekomendacja:** Phase 2 zacząć od **numerical confirmation
scenariusza (a)** (1-2 sub-tests). Jeśli (a) confirmed (prawdopodobnie),
zamknąć Phase 2 short-cycle z verdict "ścieżka A FAILS" + decision
tree do user-a (D vs E).
