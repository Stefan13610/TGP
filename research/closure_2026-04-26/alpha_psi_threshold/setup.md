# T-α: α(ψ) z ψ-threshold — Path E rozwiązanie problemu uniwersalności α dla M9.2-D

**Data:** 2026-04-26
**Status:** OPEN
**Strategic ref:** [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §M9.2-D
**Foundations binding:** [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂)
**Predecessor:** [[research/op-m92/OP_M92_P0plus_candD_multisource_results.md]]
**Goal:** Path E — α(ψ) z thresholdem przy ψ=1 (vacuum)

---

## 1. Cel

OP-M92 Phase 0+ multi-source self-consistency check (2026-04-25) wykrył
**krytyczny problem strukturalny** w naive Candidate D action: stała α ~ 0.1
w geometric units zmusza α_SI do skalowania jak M_BH² (R_S²), czyli α_SI
zmienia się o 19 rzędów wielkości między NS i M87*. **α nie jest single
physical constant pod naive coupling.**

Z 6 zaproponowanych ścieżek rozwiązania (A-F), tylko Path E pozostaje
viable:

> **Path E:** Replace constant α with α(ψ) = α₀ × (ψ - ψ_th)^n × Θ(ψ - ψ_th).
> α(ψ) activates ONLY at strong-field, gdzie ψ przekracza próg.

**Pytanie strategiczne (T-α):**
> Czy struktura α(ψ) z naturalnym thresholdem ψ_th = 1 (vacuum) i
> wykładnikiem n=2 (gładka aktywacja) rozwiązuje multi-source α-universality
> issue, jednocześnie:
> 1. zachowując scenario (e) +14.56% strong-field deviation (universal),
> 2. preserving WEP MICROSCOPE bound η < 10⁻¹⁵,
> 3. NIE wprowadzając nadmiarowych parametrów ponad konieczne minimum?

**Cel T-α:** ustalić falsifiowalny test strukturalny dla Path E.

---

## 2. Setup matematyczny

### 2.1 α(ψ) functional form

Najprostsza forma threshold function spójna z TGP:
```
α(ψ) = α₀ × (ψ - ψ_th)^n × Θ(ψ - ψ_th)         (1)
```

gdzie:
- `α₀` — dimensionless O(1) physical constant (universal)
- `ψ_th` — threshold (postulat: ψ_th = 1 = vacuum)
- `n` — wykładnik aktywacji (postulat: n = 2 dla gładkości)
- `Θ` — Heaviside step

**Postulaty TGP-natural:**

| Parametr | Wartość | Uzasadnienie |
|----------|---------|--------------|
| `ψ_th` | 1 | Vacuum point V'(Φ_eq)=0; α=0 w idealnej próżni |
| `n` | 2 | Quadratic activation: smooth (C¹), strongest WEP suppression |
| `α₀` | TBD | Calibrate z scenario (e): α(ψ_ph) → +14.56% deviation |

### 2.2 Universal scaling argument

**Photon ring uniwersalność:**

W M9.1'' Schwarzschild-like solution: photon ring przy `r_ph/M = 3.88`,
wartość pola: `ψ_ph = 1.168` (UNIVERSAL across all masses).

⇒ `α(ψ_ph) = α₀ × (1.168 - 1)² = α₀ × 0.168² = 0.0282 × α₀`
   jest **universal** (zależy tylko od ψ, nie od M_BH).

⇒ Multi-source consistency restored: M87*, Sgr A*, GW150914 mają
   THE SAME α(ψ_ph) (jako dimensionless number).

### 2.3 WEP suppression przy ψ_lab ≈ 1

**Earth surface gravity:** ψ_Earth - 1 ≈ Φ_Newton/c² ≈ 6.95×10⁻¹⁰

⇒ `α(ψ_Earth) = α₀ × (6.95×10⁻¹⁰)² = α₀ × 4.83×10⁻¹⁹`

Stosunek α_lab / α_ph = 4.83×10⁻¹⁹ / 0.0282 = **1.7×10⁻¹⁷**.

WEP MICROSCOPE marginal η ~ 1.6×10⁻¹⁵ z constant-α calibration.
Z ψ-threshold suppression: η_TGP ~ 1.6×10⁻¹⁵ × 1.7×10⁻¹⁷ ~ **2.7×10⁻³²**.

MICROSCOPE bound η < 1.1×10⁻¹⁵. Margin: **4×10¹⁶** (drastically safe).

### 2.4 Rationale dla n = 2

Trzy alternatywy dla aktywacji:

| n | Activation | C^k smoothness | WEP suppression | Verdict |
|---|------------|----------------|-----------------|---------|
| 1 | Linear | C⁰ (kink przy ψ=1) | × (6.95×10⁻¹⁰) | Marginal |
| **2** | **Quadratic** | **C¹** | **× 4.83×10⁻¹⁹** | **PREFERRED** |
| 3 | Cubic | C² | × 3.36×10⁻²⁸ | Overkill |

n=2 jest **minimal sufficient** dla:
- gładkości (C¹ — first derivative continuous przy ψ_th),
- WEP MICROSCOPE bezpieczeństwa (margin 10¹⁶),
- strong-field aktywacji (α(ψ_ph) = 0.0282 × α₀).

### 2.5 Calibration α₀ z scenario (e)

Phase 0+ structural sketch: O(1) factor `T·J·J / S_kin = ξ` przy photon ring.

Required shift +14.56% ⇔ `α(ψ_ph) × ξ ~ 0.114`.

Z α(ψ_ph) = 0.0282 × α₀:
```
0.0282 × α₀ × ξ = 0.114
α₀ = 4.04 / ξ                  (2)
```

Z ξ ~ O(1) (sketch): **α₀ ≈ 4** (postulat T-α: α₀ = 4 jako natural value).

Alternatywnie z ξ = 1: α₀ = 4.04, czyli przy n=2 i ψ_th=1, α₀ jest
**O(1) dimensionless natural constant** — bez fine-tuning.

---

## 3. Plan testów (T-α.1 ... T-α.5)

| ID | Cel | Metoda | PASS criterion |
|----|-----|--------|----------------|
| **T-α.1** | Define α(ψ) = α₀(ψ-1)² Θ(ψ-1) | sympy: piecewise function | C¹ smoothness przy ψ=1 |
| **T-α.2** | Multi-source universality | numpy: α(ψ_ph) for SgrA*, M87*, GW150914, NS | All identical (dimensionless) |
| **T-α.3** | WEP MICROSCOPE preservation | analytic: α(ψ_Earth)/α(ψ_ph) | ratio < 10⁻¹⁵ (super-safe) |
| **T-α.4** | Calibration α₀ z scenario (e) | numpy: solve α₀ × 0.0282 × ξ = 0.114 | α₀ = O(1) |
| **T-α.5** | Falsifiability — n=2 unique? | compare n=1,2,3 | n=2 minimal sufficient |

**Sukces:** 5/5 PASS → Path E structurally validated jako TGP-natural
resolution multi-source α-universality. Phase 1 covariant derivation
może proceed z α(ψ) parametrization jako baseline.

---

## 4. Falsifiable predictions

**T-α predykcje:**

1. **Multi-source EHT shadow universality:** Wszystkie BHs (Sgr A*, M87*, +5-10 LLAGN by 2030) MUSZĄ pokazać +14.56% deviation z **identyczną** wartością α₀.
   - Jeśli M87* da +12% a SgrA* +16%, T-α z n=2 sfalsyfikowane.

2. **WEP MICROSCOPE-2 (η < 10⁻¹⁷ planned):** TGP z α(ψ) prediction: **strong PASS** (η_TGP ~ 10⁻³²).
   - Jeśli MICROSCOPE-2 wykryje η ~ 10⁻¹⁶, TGP m9.2-D + path E sfalsyfikowane.

3. **NS-NS merger ringdown (LIGO 2030+):** ψ_NS surface ~ 1.4 → strong α activation.
   - Predykcja: α(ψ_NS) × T·J·J → measurable inspiral phase shift.
   - Falsifiable: brak shifu = TGP m9.2-D incorrect.

4. **Solar system PPN:** Sun surface ψ - 1 ~ 2×10⁻⁶ → α(ψ_sun) ~ α₀ × 4×10⁻¹².
   - Cassini/Shapiro time-delay constraint γ - 1 < 2.3×10⁻⁵.
   - TGP α-induced γ deviation ~ α(ψ_sun) ~ 10⁻¹¹ → **PASS by 10⁶ margin**.

---

## 5. Status: structural resolution z postulate parameters

T-α jest **STRUCTURAL POSTULATE z trzema parametrami:** ψ_th = 1, n = 2, α₀ = O(1).

**Co T-α daje:**
- Resolution multi-source α-universality issue (M_BH²-scaling neutralized).
- Drastic WEP MICROSCOPE safety margin (10¹⁶ extra).
- Universal scenario (e) deviation +14.56% preserved.
- Falsifiable predictions na multi-source EHT i WEP-2.

**Co T-α NIE daje (otwarte dla Phase 1):**
- First-principles wyprowadzenie ψ_th = 1 z TGP foundations
  (motywacja: vacuum point, ale brak rigorous derivation).
- First-principles wyprowadzenie n = 2 (motywacja: minimal smooth,
  ale alternatywne wartości n nie są blokowane przez aksjomaty).
- Determinacja α₀ z RG flow (postulat α₀ = 4 jest heuristic).
- Pełna covariant action S[Φ, g, ...] z wyłaniającym się α(ψ).

---

## 6. Pliki

- `setup.md` (ten plik) — design audytu T-α
- `alpha_psi_threshold.py` — sympy + numpy script T-α.1..T-α.5
- `alpha_psi_threshold.txt` — raw output
- `results.md` — synteza po wykonaniu

---

## 7. Cross-references

- [[research/op-m92/OP_M92_P0plus_candD_multisource_results.md]] (problem origin)
- [[research/op-m92/OP_M92_P0plus_candD_wep_results.md]] (WEP baseline)
- [[research/op-m92/op_m92_T1_action_reverse_engineer.py]] (T-M92.1 universality)
- [[research/closure_2026-04-26/f_psi_principle/results.md]] (f(ψ) principle T-FP)
- [[research/closure_2026-04-26/Lambda_from_Phi0/results.md]] (Λ from Φ₀ = T-Λ)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §M9.2-D
