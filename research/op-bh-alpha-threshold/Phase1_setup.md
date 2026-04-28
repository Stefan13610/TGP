---
title: "Phase 1 setup — multi-source dimensional + mass-scaling audit of α"
date: 2026-04-28
cycle: BH.1.Phase1
status: PRE-EXECUTION
predecessor: "[[program.md]]"
tags:
  - TGP
  - BH
  - alpha-psi
  - dimensional-analysis
  - mass-scaling
  - sympy
---

# Phase 1 — Setup: multi-source dimensional + mass-scaling audit of α (Candidate D)

> **Cel:** Sformalizować i zamknąć **Phase 0+ multi-source ISSUE** —
> pokazać że α z Candidate D action `S_int ∝ α · T^μν J_μ J_ν` **nie może**
> być pojedynczą stałą fizyczną pod żadną prostą unit-bridge transformation
> (analog SC.1.Phase1: H₀ unit-cousin α_PB ↔ α_0 → REJECTED).

---

## H₀ vs H₁ recap

| Hipoteza | Twierdzenie | Falsyfikacja |
|---------|-------------|--------------|
| **H₀** (α single physical constant) | Istnieje pojedyncza wartość α_universal (in SI) i naturalny czynnik konwersji f, taki że α_geom ~ 0.1 (Phase 0+ heurystyka) reprezentuje tę samą fizykę dla wszystkich źródeł | M² scaling α_SI = α_geom · R_S² nie jest absorbowany przez żaden strukturalny f bez nielokalnych modyfikacji |
| **H₁** (α(ψ) wymagana) | α z Candidate D **musi** być funkcją scalar-field α(ψ) z progiem (Path E), bo tylko ψ-dependence absorbuje uniwersalność strong-field przy lab suppression | Brak — H₁ jest negatywnym wnioskiem z odrzucenia H₀ + Path A/B/C |

---

## Procedura testów

### T1.1 — Pełna analiza wymiarowa Candidate D action (sympy `dimsys_SI`)

**Cel:** Wyprowadzić wymiar SI α z Candidate D action; sprawdzić co
oznacza "α_geom ~ 0.1 in M_BH=1 units" w SI dla każdego źródła osobno.

**Plan:**
1. Action: `S_int = α · ∫ T^μν J_μ J_ν √-g d⁴x` (Candidate D minimal coupling).
2. Wymiar T^μν: `[T^μν] = J/m³` (gęstość energii). Wymiar J_μ (4-current): `[J_μ] = A·s/m³` lub w geom = M/L³ (substrate flux).
3. Wymiar `T^μν J_μ J_ν √-g d⁴x`: po pełnej redukcji w SI dim[α] musi być takie, żeby S_int [J·s] (akcja).
4. Demonstracja sympy: dim[α] **nie jest** dimensionless w SI — ma wymiar `[powierzchni]` (m²) lub `[powierzchnia × inne]` (zależy od kalibracji T·J·J).
5. To wyjaśnia dlaczego α_SI = α_geom · R_S² ma wymiar m² i skaluje się jak M_BH².

**Wniosek oczekiwany:** dim[α] ≠ dimensionless → α nie ma natury pojedynczej "fundamentalnej stałej". Jakakolwiek "uniwersalna" wartość wymaga reskalowania źródła-zależnego.

### T1.2 — Explicit M² scaling demonstration (4 źródła)

**Cel:** Replikować i sformalizować Phase 0+ wynik w czystej tabelarycznej formie.

**Plan:**
1. Dla 4 źródeł (Sgr A*, M87*, GW150914 final, NS) policzyć:
   - R_S = 2GM/c² [m]
   - α_geom heurystyka = 0.1 (universal w geom units, M_BH=1)
   - α_SI = α_geom · R_S² [m²]
2. Tabela: span od 1.71·10⁶ m² (NS) do 3.69·10²⁵ m² (M87*) — 19 orders.
3. Pearson check: log10(α_SI) vs log10(M_BH/M_⊙) — slope musi być 2.000 (idealne M² scaling).

**Wniosek oczekiwany:** dokładnie M² scaling, span 19 orders dla 9.5 orders w M_BH.

### T1.3 — Closure paths A/B/C explicit fail demonstration

**Cel:** Pokazać że trzy fail-paths z Phase 0+ rzeczywiście zawodzą strukturalnie (nie tylko heurystycznie).

**Plan:**
1. **Path A (L_TGP intrinsic length):** `α_dim · L_TGP² = α_SI` → wymagane `L_TGP(SgrA*) = sqrt(α_SI/α_dim)` dla każdego źródła. Pokaż że L_TGP-source nie jest stałe.
2. **Path B (substrate density ρ_sub ~ M/R_S³):** wprowadź `T·J·J ~ ρ_sub/R_S²` i pokaż że α scaling pozostaje niezmieniony (failure mode identyczny z naive).
3. **Path C (M_Pl⁻² prefactor):** `α dimensionless = (α_SI / M_Pl⁻²) M_Pl⁻²`. Pokaż że wymagane α jest astronomicznie unphysical (~10⁸⁹) dla ratio M_SgrA*/M_Pl.
4. Każda path: jednoznaczny strukturalny fail z numerical demonstration.

**Wniosek oczekiwany:** Path A/B/C wszystkie **strukturalnie** fail. Path D nie-valid (geometria photon ring uniwersalna). Path E (α(ψ)) jako jedyna żywotna kontynuacja.

### T1.4 — Lab-suppression check (Path E sanity)

**Cel:** Sprawdzić że pod hipotezą Path E `α(ψ) = α_0 · (ψ - ψ_th)ⁿ Θ`, lab-scale ψ_lab ≈ 1 daje natural eksponentialne/polynomialne tłumienie eksperymentów lab.

**Plan:**
1. Dla typowych ψ_th ∈ [1.01, 1.05, 1.10] i n ∈ [1, 2, 3]:
   - α(ψ_lab=1) = 0 (Heaviside cuts off)
   - α(ψ_Earth_surface ~ 1 + 7·10⁻¹⁰) — czy tłumi eksperymenty WEP?
   - α(ψ_lab) / α(photon ring ψ=1.168) — ratio.
2. Demonstracja że Path E **automatycznie** rozwiązuje multi-source ISSUE bez fine-tuningu jeśli ψ_th ∈ [1.05, 1.15].

**Wniosek oczekiwany:** Path E z ψ_th ≥ 1.05 daje całkowite lab suppression, aktywuje się tylko strong-field.

### T1.5 — Path E vs Path D/F summary table

**Cel:** Klasyfikacja strukturalnie żywotnych dróg dalej.

**Plan:**
1. Path D (target nie-uniwersalny): pokaż że scenariusz (e) target +14.56% jest **uniwersalny** w geom units (r_ph^TGP/r_ph^GR = 3.88/3.0 niezależnie od M_BH) → nie-valid.
2. Path E (α(ψ)): viable, n + ψ_th to fit-parametry do Phase 2.
3. Path F (composite D + B): możliwa hybrid, ale wymagałaby DWÓCH nowych mechanizmów. Path E preferowana z parsimony.
4. Tabela podsumowująca: co dalej zostaje + jakie predykcje.

**Wniosek oczekiwany:** Path E jako jedyna żywotna minimal-extension; Path F backup jeśli Path E falsifies w Phase 3.

---

## Materiał wykonawczy

- **Skrypt:** `phase1_multi_source_audit.py` (sympy dimensional analysis + numerical multi-source M² check)
- **Output:** `phase1_multi_source_audit.txt` (full execution log)
- **Audit memo:** `Phase1_results.md` (post-execution conclusion)

## Środowisko

Python 3.x + sympy + numpy. Brak zewnętrznych zależności.

```bash
cd TGP/TGP_v1/research/op-bh-alpha-threshold
PYTHONIOENCODING=utf-8 python -X utf8 phase1_multi_source_audit.py 2>&1 | tee phase1_multi_source_audit.txt
```

## Co znajdzie się w wynikach

- **PASS** dla każdego sub-testu T1.1–T1.5 jeśli H₁ poparte (oczekiwane).
- **FAIL** dla T1.1/T1.3 oznaczałoby ukrytą unit-bridge konwersję niewykrytą w Phase 0+ — zaskoczenie, do raportu.

---

## Cross-references

- `program.md` — overall 3-phase plan
- `../op-m92/OP_M92_P0plus_candD_multisource_results.md` — multi-source ISSUE źródłowy
- `../op-m92/OP_M92_Phase1_PLAN.md` — pełny PLAN (świadomie NIE duplikujemy)
- `../op-sc-alpha-origin/Phase1_setup.md` — wzorzec mini-cyklu (analog dimensional audit)

---

## Constants used

```
G  = 6.67430e-11    m³ kg⁻¹ s⁻²
c  = 2.99792458e8   m s⁻¹
M_⊙ = 1.98892e30    kg
M_Pl = 2.176e-8     kg  (Planck mass)

Sources:
  Sgr A*       M = 4.3e6  M_⊙
  M87*         M = 6.5e9  M_⊙
  GW150914 fin M = 65     M_⊙
  NS canonical M = 1.4    M_⊙
```
