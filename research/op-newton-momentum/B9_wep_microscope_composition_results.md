# B9 — M9.2 WEP MICROSCOPE composition test (two-component)

**Data:** 2026-05-01
**Autor:** Mateusz (zapis: Claudian)
**Status:** ✅ **B9 CLOSED 2026-05-01 — 6/6 PASS**
**Skrypt:** [[B9_wep_microscope_composition.py]] (~430 linii sympy + numpy + scipy)
**Output:** [[B9_wep_microscope_composition.txt]]
**Audit binding:** [[../../meta/AUDYT_TGP_2026-05-01.md]] § B (HIGH, B9), § M
**Predecessor:** [[M9_2_results.md]] (single-source M² universality, NOT two-component)

---

## TL;DR

- ✅ Sympy LOCK: `m_field ~ (qM)²/(4π·σ)` (weak-field) → `m_field/m_g ~ q·U_surface` linear in q oraz σ⁻¹.
- ✅ **Geometryczna kompozycja (Pt vs Ti, same q, M)**: η_geom (natural) = 6.68×10⁻¹ (sigma scaling 1.0 vs 1.68).
- ✅ **Coupling kompozycja (δq/q wariation)**: linear scaling potwierdzone — `eta_AB / (δq/q) → 0.95...0.9995` (5% rel_dev).
- ✅ **Inhomogeneous ρ (core+shell vs pure Gaussian, same M)**: η_inhom = 0.289 (structure-sensitive).
- ✅ **Realistic Pt vs Ti at MICROSCOPE**: η_TGP_lab = **1.32×10⁻²⁶** (without T-α), **6.49×10⁻⁴⁵** (with T-α).
- ✅ **Margins**: MICROSCOPE 2017 (1.1×10⁻¹⁵): 8.3×10¹⁰× safe (no T-α), 1.7×10²⁹× safe (with T-α). MICROSCOPE-2 future (10⁻¹⁷): 7.6×10⁸× and 1.5×10²⁷×.

→ **Audit B9 demand "test kompozycji (dwa źródła o różnym q lub niejednorodnym ρ)" satisfied w OBU TRYBACH.**

---

## 1. Cel testu

M9.2 (`M9_2_results.md`) zamknięte 5/5 PASS, ale verdict M9.2.3 ("WEP intrinsic") oparty był wyłącznie na **single-source** scaling: `m_field ∝ M²` (`rel_dev = 4.65×10⁻¹²` machine precision). Audit critique B9:

> *"M9.2 `m_field/m_g = 3.48×10⁻²` extrapolowane jako WEP MICROSCOPE 10⁻¹⁵ — nie testowane dwuskładnikowo."*
> Akcja: *"Test kompozycji (dwa źródła o różnym q lub niejednorodnym ρ)."*

**Issue:** `m_field ∝ M²` przy fixed `q, σ` to nie jest test WEP. WEP testuje **różnicę przyspieszenia dwóch ciał** w zewnętrznym polu grawitacyjnym. Te ciała mogą się różnić:
1. **q** (efektywny TGP coupling per mass) — bezpośrednia 5-th force violation
2. **σ** (size/density) — geometryczna kompozycja (np. Pt vs Ti przy same M)
3. **ρ structure** (homogeneous vs core+shell) — strukturalna kompozycja

B9 closure musi: (a) sympy-LOCK leading scaling η_AB w obu trybach; (b) numerycznie zweryfikować na BVP (M9.2 Yukawa eq); (c) projektować na realistyczne MICROSCOPE Pt vs Ti; (d) sprawdzić margines vs experimental bound 1.1×10⁻¹⁵.

---

## 2. Metoda

**Skrypt:** `B9_wep_microscope_composition.py` (~430 linii)

| Krok | Test |
|---|---|
| **Step 1** | Sympy LOCK: `m_field ~ (qM)²/σ` w weak-field, `η_AB ~ |δq/q| + |δσ/σ|` |
| **Step 2** | Geometryczna kompozycja: Pt-like (σ=1.0) vs Ti-like (σ=1.68), same q, M |
| **Step 3** | Coupling kompozycja: `δq/q ∈ {0.1, 0.01, 0.001}`, weryfikacja linear scaling |
| **Step 4** | Inhomogeneous ρ: pure Gaussian vs core+shell (60% σ=0.5 + 40% σ=1.5), same M |
| **Step 5** | Realistic MICROSCOPE Pt vs Ti — lab projection (`U_lab ~ 2×10⁻²⁶`) |
| **Step 6** | Combined T-α suppression `α(ψ_lab) ~ (7×10⁻¹⁰)² ~ 5×10⁻¹⁹` + future MICROSCOPE-2 |

BVP: ten sam Yukawa solver co `m9_2_momentum.py` (scipy.solve_bvp na `v(r) = r·ε(r)`).

---

## 3. Wyniki — 6/6 PASS

### 3.1 Step 1 — Sympy LOCK m_field scaling ✅ PASS

**Asymptotyczna analityczna integracja:**
```
∫_σ^∞ 4π r² · (qM/(4π r²))² · exp(-2√β·r) dr  =  M²√β q² Ei(2√β σ e^iπ)/(2π) - i M²√β q²/2 + M²q² exp(-2√β σ)/(4π σ)
limit β→0  =  M² q² / (4π σ)
```

**Wnioski leading scaling:**
- `m_field ~ (qM)² / (4π σ)` — verified
- `m_field/m_g = m_field/(qM) ~ qM/(4π σ) = q·U_surface` (quasi-linear w q i σ⁻¹)
- Dla dwóch ciał (same M, różne q, σ):
  ```
  η_AB ~ |δq/q| + |δσ/σ|   (LEADING)
  ```

### 3.2 Step 2 — Geometryczna kompozycja (Pt vs Ti) ✅ PASS

| Test body | σ (natural) | M | q | m_field | m_field/m_g |
|---|---|---|---|---|---|
| Pt-like | 1.000 | 1.0 | 1.0 | 3.483×10⁻² | 3.483×10⁻² |
| Ti-like | 1.680 | 1.0 | 1.0 | 1.740×10⁻² | 1.740×10⁻² |

**η_geom (natural units)** = `|Δ(m_field/m_g)| / mean = 6.676×10⁻¹` (≈ 50% asymetria w natural).

**Predicted leading**: `|Δσ/σ_avg| = 0.508` → numerical/predicted = **1.32** (O(1) zgodne z `m_field/m_g ~ qM/σ`).

### 3.3 Step 3 — Coupling kompozycja (δq variation) ✅ PASS

| δq/q | q_A | q_B | η_AB | η_AB / (δq/q) |
|---|---|---|---|---|
| 1.00×10⁻¹ | 1.0000 | 1.1000 | 9.524×10⁻² | **0.9524** |
| 1.00×10⁻² | 1.0000 | 1.0100 | 9.950×10⁻³ | **0.9950** |
| 1.00×10⁻³ | 1.0000 | 1.0010 | 9.995×10⁻⁴ | **0.9995** |

**Linearity verification**: `eta_AB / (δq/q) → 1` przy `δq → 0` (rel_dev = 4.80×10⁻²).

→ **Confirms sympy LOCK**: `η_AB ~ (m_field/m_g) × (δq/q)` linearly. Drobne O(δq/q) odchyłki przy `δq/q = 0.1` to nonlinear corrections — w realistycznym lab regimie (`δq/q < 10⁻³`) reżim linearny nadweryfikowany.

### 3.4 Step 4 — Inhomogeneous ρ (core+shell vs pure Gaussian) ✅ PASS

| Distribution | M (total) | m_field | m_field/m_g |
|---|---|---|---|
| Pure Gaussian σ=1.0 | 1.0 | 3.483×10⁻² | 3.483×10⁻² |
| Core+shell (60% σ=0.5 + 40% σ=1.5) | 1.0 | 4.658×10⁻² | 4.658×10⁻² |

**η_inhomogeneous** = 2.886×10⁻¹.

**Interpretacja:** Ten sam całkowity M, ale różny rozkład → non-zero η bo m_field zależy od ⟨1/σ⟩ struktury (sympy LOCK Step 1). Core+shell ma większe `m_field` bo dense core dominuje (`σ_core = 0.5 < σ_pure = 1.0`).

→ **Captures audit demand "inhomogeneous ρ"** — `η ≠ 0` przy fixed M, identyczne q. Realistycznie, materiały o różnym lattice/skład mają różne ρ-struktury → nieistotne korekty O(10⁻²)·U_lab.

### 3.5 Step 5 — Realistic MICROSCOPE Pt vs Ti ✅ PASS

**Real MICROSCOPE setup:**
- Pt: ρ = 21.45 g/cm³, m ~ 0.4 kg, R ~ 1.5 cm
- Ti: ρ = 4.51 g/cm³, m ~ 0.4 kg, R ~ 2.5 cm
- Free fall in Earth's gravity at 700 km altitude.

**Lab-scale `U_surface`:**
- `U_Earth = G·M_E/(c²·R_E) = 6.96×10⁻¹⁰`
- `U_lab_self = G·m_lab/(c²·R_lab) = 1.98×10⁻²⁶`

**TGP eta_AB lab-projection:**
```
η_TGP_geom_lab     = U_lab_self × η_geom_natural          = 1.98×10⁻²⁶ × 0.668 = 1.32×10⁻²⁶
η_TGP_coupling_lab = U_lab_self × |δq/q|_5th_force_bound  = 1.98×10⁻²⁶ × 10⁻³  = 1.98×10⁻²⁹
η_TGP_total_lab    = η_geom + η_coupling                  = 1.32×10⁻²⁶
```

**MICROSCOPE bound 2017:** η < 1.1×10⁻¹⁵ → **safety margin 8.31×10¹⁰× SAFE**.

**MICROSCOPE-2 future:** η < 10⁻¹⁷ → **safety margin 7.55×10⁸× SAFE**.

### 3.6 Step 6 — T-α suppression + MICROSCOPE-2 combined ✅ PASS

**closure_2026-04-26 Phase 4 T-α:** `α(ψ_lab) ~ (ψ_lab)² ~ (7×10⁻¹⁰)² = 4.90×10⁻¹⁹`.

**Combined:**
```
η_TGP_combined = η_TGP_lab × α(ψ_lab) = 1.32×10⁻²⁶ × 4.90×10⁻¹⁹ = 6.49×10⁻⁴⁵
```

**Margins post-T-α:**
- MICROSCOPE 2017 (1.1×10⁻¹⁵): **1.70×10²⁹× SAFE**
- MICROSCOPE-2 (10⁻¹⁷): **1.54×10²⁷× SAFE**

---

## 4. B9 FINAL VERDICT

```
[PASS] Step 1 (sympy scaling LOCK)
[PASS] Step 2 (geometric Pt vs Ti, same q)
[PASS] Step 3 (coupling delta_q variation)
[PASS] Step 4 (inhomogeneous rho core+shell)
[PASS] Step 5 (MICROSCOPE realistic Pt vs Ti)
[PASS] Step 6 (T-alpha + MICROSCOPE-2 combined)
Total: 6/6 PASS
```

> **B9 CLOSURE:** Two-component WEP test PASSES we wszystkich trybach (geometryczny σ, coupling q, inhomogeneous ρ). Realistic Pt vs Ti at MICROSCOPE: `η_TGP ≈ 10⁻²⁶` (no T-α), `≈ 10⁻⁴⁵` (with T-α) — both **wiele rzędów wielkości** poniżej eksperymentalnego bound 1.1×10⁻¹⁵. Audit B9 demand "test kompozycji (dwa źródła o różnym q lub niejednorodnym ρ)" satisfied **eksplicit** w obu modach.

---

## 5. Falsifiable predictions (post-B9)

1. **MICROSCOPE-2 (10⁻¹⁷):** TGP M9.2 prediction `η_TGP ~ 10⁻²⁶` (geom+coupling, no T-α) lub `~ 10⁻⁴⁵` (w T-α). PASS z marginesem ≥ 10⁹× (obowiązkowe nawet bez T-α).

2. **STE-QUEST (10⁻¹⁸ projected):** Margines ≥ 10⁸× (bez T-α). Nadal SAFE.

3. **Falsyfikacja:** Detekcja composition-dependent `η ∈ [10⁻²⁶, 10⁻¹⁷]` byłaby strukturalnie problematyczna dla TGP — wskazywałaby na **brak** Phase-4 T-α suppression i / lub `δq/q > 10⁻³` (sprzeczne z 5-th force constraints).

4. **Inhomogeneous-ρ test:** Material composite (np. nano-composite Pt-Ti core-shell) z `δ⟨1/σ⟩` ~ 30% zmierzony w sub-fm precyzji — nie zaobserwowane.

---

## 6. Outstanding follow-ups (post-B9)

| Item | Status | Priorytet |
|---|---|---|
| Direct sympy LOCK na full BVP (zamiast asymptotic Yukawa) | nice-to-have | LOW |
| Two-component compound MICROSCOPE projection (Pt-Rh + Ti-V alloys) | extension | LOW |
| Strong-field WEP — neutron stars i.e., extended η scaling | dziedziczy z B6 § U.5 | MEDIUM (NS-NS ringdown) |
| Composition-Φ self-coupling explicit derivation z L_mat | dziedziczy z A4 (ax:metric-coupling) | HIGH (A4 dedicated) |

**B9 itself: CLOSED.** Audit credit applied (39/43 → 40/43 = 93% closed).

---

## 7. Pliki

- **Skrypt:** [[B9_wep_microscope_composition.py]]
- **Output:** [[B9_wep_microscope_composition.txt]]
- **Predecessor:** [[M9_2_results.md]] (M9.2.3 single-source baseline; post-B9 closure marker)
- **Audit § W (B9 closure):** [[../../meta/AUDYT_TGP_2026-05-01.md]] § W
- **closure_2026-04-26 T-α:** [[../closure_2026-04-26/alpha_psi_threshold/results.md]] (Phase 4 suppression)
- **B6 closure (precedensowy):** [[B6_m9x_sqrtg_rerun_results.md]]
- **B8 closure (precedensowy):** [[B8_lagrangean_independence_check_results.md]]

---

## 8. Podsumowanie w jednym zdaniu

B9 pod sympy LOCK + numerical BVP pokazuje, że dwukomponentowy test kompozycji WEP (Pt vs Ti via geometryczne `σ`, coupling `q`, inhomogeneous `ρ`) daje `η_TGP ≈ 1.32×10⁻²⁶` w realistycznym lab (6.49×10⁻⁴⁵ z T-α suppression) — wiele rzędów wielkości poniżej MICROSCOPE 2017 (1.1×10⁻¹⁵, margines ≥ 8.3×10¹⁰×) i przyszłej MICROSCOPE-2 (10⁻¹⁷, margines ≥ 7.5×10⁸×); audit B9 demand "test kompozycji (dwa źródła o różnym q lub niejednorodnym ρ)" satisfied explicit w obu modach.
