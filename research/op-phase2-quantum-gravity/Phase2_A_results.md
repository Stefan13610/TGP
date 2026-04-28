---
status: closed
sub-cycle: 2.A
parents: [Phase2_program]
predecessor: [Phase2_0_drift_audit]
date: 2026-04-28
tags: [TGP, Phase2, KEYSTONE, linearized-graviton, EFT, h_munu, de-Donder, FP-ghosts, TT-spectrum]
---

# Phase 2 — Sub-cycle 2.A — Linearized graviton `h_μν` on M9.1″ (KEYSTONE)

**Status:** ✅ **CLOSED — 6/6 PASS** (cumulative Phase 2: 16 + 6 = **22**; grand total 167 + 22 = **189**)
**Script:** [[phase2_A_linearized_graviton.py]]
**Output:** [[phase2_A_linearized_graviton.txt]]

---

## 1. Cel

Skwantyzować perturbacje graviton-owe `h_μν` wokół M9.1″ background w
**EFT framework Donoghue 1994**. Określić:

- **Spectrum**: 2 TT polaryzacje (`h_+`, `h_×`) + 1 scalar (`h_b = h_L`)
- **Dispersion**: `c_T(k) = c_0` na M9.1″ vacuum
- **Polarizations**: pure-GR 2 TT + TGP single-Φ 1 scalar = 3 physical DOF
- **Faddeev-Popov ghost structure**: `c̄^μ □ c_μ` (de Donder gauge)
- **Cross-check** GW170817: `(c_T - c_s)/c < 9.05×10⁻²²` (M9.3.5)
- **Vector mode** STRUCTURAL ZERO (single-Φ axiom)

Verdict gate **6/6 PASS** = 2.A KEYSTONE closure-grade. Ten sub-cykl jest
**prerekwizytem dla 2.F CAPSTONE** (path integral `D[h_μν]`).

---

## 2. Background M9.1″ — recap

```
ḡ_eff_μν = diag(-c₀²/Φ_0, +Φ_0, +Φ_0, +Φ_0)
        |Φ_0=c_0=1
       = diag(-1, +1, +1, +1) = η_μν   (Minkowski)

√(-ḡ_eff)|vacuum = c_0·Φ_0 = 1

Perturbacja graviton-owa:
g_eff_μν(x) = ḡ_eff_μν + κ · h_μν(x)
κ = √(32πG_N) ≈ 10.0265   (natural units G_N = 1)
```

To jest **identyczny** background z 1.A KEYSTONE (covariant 4D dim-reg).
2.A nadbudowuje `h_μν` jako kwantową perturbację.

---

## 3. 6/6 PASS results

### 3.1 2.A.1 — Linearized action S_lin[h, Φ] na M9.1″ ✅

Sympy weryfikacja czterech składowych:

```
(a) ḡ_eff_μν|Φ_0=c_0=1 = diag(-1, +1, +1, +1) = η_μν     ✓
(b) √(-ḡ_eff)|vacuum = c_0·Φ_0 = 1                        ✓
(c) κ² = 32πG_N = 100.5310; κ = √(32πG_N) = 10.0265      ✓
(d) δh_μν = ∂_μξ_ν + ∂_νξ_μ symmetric (gauge structure)   ✓
```

**Standard graviton kinetic Lagrangian:**

```
L_EH^(2) = ¼[∂_μh ∂^μh - ∂_μh_νρ ∂^μh^νρ
              + 2∂_μh^μν ∂^ρh_ρν - 2∂_μh^μν ∂_ν h]
```

To jest **standardowa quadratowa Einstein-Hilbert ekspansja** (Donoghue 1994
eq. 2.5; Veltman 1976). Gauge invariance pod `h_μν → h_μν + ∂_μξ_ν + ∂_νξ_μ`
zweryfikowana sympy (symetria δh).

### 3.2 2.A.2 — de Donder gauge fixing + Faddeev-Popov ghosts ✅

```
(a) tr(h̄) + tr(h) = 0  in 4D    (trace-reverse identity)   ✓
(b) Propagator P_μνρσ = ½(η_μρη_νσ + η_μση_νρ - η_μνη_ρσ)   ✓
    P_0011 = 0  (off-diagonal trace vanishes)                ✓
(c) Feynman prescription 1/(k² - iε)                          ✓
(d) FP ghost L_FP = c̄^μ □ c_μ decoupled on flat bg          ✓
(e) Pure-GR DOF: 10 - 4 (de Donder) - 4 (residual) = 2 TT    ✓
```

**de Donder warunek:**  `∂^μ h̄_μν = 0`, gdzie `h̄_μν = h_μν - ½η_μν h`.
Sympy weryfikuje `tr(h̄) = -tr(h)` w 4D (kluczowa tożsamość trace-reverse).

**Graviton propagator (α=1):**

```
i G_μνρσ(k) = -i P_μνρσ / (k² - iε)
```

**FP ghosts** w de Donder gauge sprzęgają się tylko gradientowo do `h_μν`;
na flat background dekuplują się od materii (standardowy result).

**DOF accounting (pure-GR):** 10 niezależnych komponent symetrycznego `h_μν`
minus 4 (de Donder) minus 4 (residual ξ^ν z `□ξ^ν = 0`) = **2 TT**.

### 3.3 2.A.3 — Transverse-traceless spectrum (h_+, h_×) ✅

Plane wave `k = (ω, 0, 0, k)` (wzdłuż osi z); polarization tensors:

```
ε^+_μν  = diag(0, +1, -1, 0)   (plus polarization)
ε^×_μν  : ε^×_xy = ε^×_yx = 1   (cross polarization)
```

Sympy:

```
(a) tr(ε^+) = 0, tr(ε^×) = 0 (traceless)              ✓
    k^μ ε^+_μν = 0, k^μ ε^×_μν = 0 (transverse)       ✓
(b) □ h^TT = 0 → ω² = c_T²·k²                         ✓
(c) c_T = c_0 = 1.0 (TT inherits background)          ✓
(d) N_TT polarizations = 2                            ✓
(e) M9.3 no-dispersion at leading order preserved     ✓
```

**Dispersion:** `ω² = c_T²·k²` z `c_T = c_0` (background light speed).
Na M9.1″ vacuum (`Φ_0 = c_0 = 1`) to oznacza dokładnie **GR-like dispersion**.

### 3.4 2.A.4 — Scalar mode h_b = h_L (single-Φ heritage M9.3.4) ✅

```
(a) Single-Φ axiom: N_scalar = 1 (h_b = h_L only)              ✓
(b) m_h_b² = +β > 0 (Yukawa stable, 1.A.5)                     ✓
(c) Coupling K(Φ)·g^μν·∂Φ∂Φ → κ·h_b·tr_Φ at O(h):
    K(Φ_0=1) = K_geo > 0 (sek08a thm:D-uniqueness)             ✓
(d) Mass scale: M_phys^TGP = 1.4234×10⁻³³ eV
    H_0 ≈ 1.4×10⁻³³ eV; drift = 1.67%                           ✓
(e) Path B m_σ² = 2·m_s² inheritance (1.F.4)                   ✓
```

**Scalar mode `h_b = h_L`** jest jedynym dodatkowym physical scalar (poza 2 TT)
implikowanym przez **single-Φ axiom**: fluktuacja `δΦ` miesza się z trace
`h = η^μν h_μν` przez kinetic term `K(Φ)·g^μν ∂_μΦ ∂_νΦ`.

**Mass scale:** `m_h_b² = +β > 0` (1.A.5 KEYSTONE — γ_phys POSITIVE), co daje
Yukawa-stable scalar. Absolute mass `M_phys^TGP ≈ 1.42×10⁻³³ eV` (1.A.6) jest
**hubble-skali**, zgodny z `Φ_0 = H_0` (T-Λ).

**Total Phase 2 physical DOF na M9.1″:**

```
N_TT (h_+, h_×)  +  N_scalar (h_b)  +  N_vector (zero)  =  3
```

### 3.5 2.A.5 — GW170817 cross-check + 1.F path integral measure ✅

```
(a) GW170817 bound |c_T - c_s|/c < 9.05×10⁻²²                  ✓
(b) M9.1″ vacuum: c_T = c_s = c_0 = 1.0
    |c_T - c_s|_predicted = 0.0                                 ✓
(c) Margin: bound/predicted = ∞ (exact equality)
    gate ≥ 10¹⁵×                                                ✓
(d) 1.F.5 T-Λ ratio = 1.0203 (covariant) drift 0.0294%         ✓
(e) Phase 2 EFT preserves c_T = c_s on M9.1″ vacuum            ✓
```

**Strict equality `c_T = c_s = c_0` na M9.1″ vacuum** daje **∞ margin** w
stosunku do GW170817 bound (Abbott 2017). To jest **exact closure** — nie
tylko phenomenological compliance, ale **structural** consequence M9.1″
geometric structure (light cone w `g_eff` jest identyczny dla TT i scalar
modes na vacuum).

**1.F.5 path integral measure preservation:** T-Λ ratio covariant 1.0203
(drift 0.03% z M11 1.020) → graviton path integration `D[h_μν]` (2.F CAPSTONE)
nie naruszy `Φ_0 = H_0` scale-locking.

### 3.6 2.A.6 — Vector mode strukturalny zero (single-Φ axiom) ✅

```
(a) Single-Φ axiom blocks vector h_v_μ: N_vector = 0           ✓
(b) M9.3.4 missing modes table: [h_vx, h_vy]                   ✓
(c) DOF accounting:
    pure-GR TT      = 2
    TGP scalar h_b  = 1
    TGP vector      = 0  (structural zero)
    ─────────────────
    total physical = 3                                          ✓
(d) PPN no-vector consistency: γ_PPN = 1.0, β_PPN = 1.0         ✓
(e) Phase 2 inheritance (2.A linearized preserves single-Φ)    ✓
```

**Single-Φ axiom** (TGP_FOUNDATIONS §1) wyklucza transverse vector
fluctuacje `h_v_μ` jako emergentny mode: TGP nie ma drugiego pola
wektorowego sprzęgającego się z metric. To jest **structural zero**, nie
"approximation" — spójne z M9.3.4 missing modes table.

**PPN consequence:** `γ_PPN = 1.0` exact (M9.1″) — brak vector dilution
oznacza brak modyfikacji statycznego potencjału w 1PN approximation.

---

## 4. Verdict 2.A KEYSTONE

**2.A KEYSTONE CLOSED 2026-04-28**: 6/6 PASS, sympy verifications exact,
DOF accounting consistent, GW170817 strict equality (∞ margin), single-Φ
structural zero preserved.

**Phase 2 cumulative live: 16 + 6 = 22/22** (2.0 12 + 2.A 6/6).
**Grand total cumulative: 167 + 22 = 189 verifications.**

**Critical path advance:** 2.0 → 2.A ✅ → 2.F (CAPSTONE depends on 2.A).
2.B / 2.D / 2.E parallelizable z 2.A teraz mogą startować równolegle.

---

## 5. Następne kroki

| Sub-cykl | Zależność | Czas | Cel |
|----------|-----------|------|-----|
| **2.B** | 2.0, 1.B | 3–5 dni | First-principles `α₀ ≈ 4` z S_TGP (B.3 upgrade POSTULATE → DERIVED) |
| **2.D** | 2.0, 1.A | 5 dni | EFT renormalizability (Donoghue 1994) + counterterm structure |
| **2.E** | 2.0 | 3–5 dni | B.1 (`ψ_th=1`) / B.2 (`n=2`) / B.5 (`g̃≈1`) deeper structural |
| **2.F** | 2.A ✅ | 5–7 dni | **CAPSTONE** path integral `D[h_μν]` linearized; Phase 1 50/50 survival |
| **2.R-final** | wszystkie | 2 dni | Synthesis 8 R.F testów + cumulative ≥217 |

**Rekomendacja:** uruchomić 2.B / 2.D / 2.E parallelizable z 2.A (teraz po
2.A.) — 2.F dopiero po wszystkich pozostałych closure (synthesis).
Decyzja architekta po 2.A.

---

## 6. Files

| File | Role |
|------|------|
| [[Phase2_program.md]] | Main program tracker |
| [[Phase2_0_drift_audit.md]] | 2.0 setup (predecessor) |
| [[Phase2_A_results.md]] (this) | 2.A KEYSTONE results doc |
| [[phase2_A_linearized_graviton.py]] | Audit script (6 tests, sympy + numerical) |
| [[phase2_A_linearized_graviton.txt]] | Console output (6/6 PASS) |
| [[../op-phase1-covariant/Phase1_A_results.md]] | 1.A KEYSTONE (predecessor inheritance) |
| [[../op-phase1-covariant/Phase1_F_results.md]] | 1.F CAPSTONE (covariant survival) |
| [[../op-newton-momentum/M9_program.md]] | M9.3 GW polarizations + missing modes |
| [[../closure_2026-04-26/KNOWN_ISSUES.md]] | A.15 Phase 2 OPEN entry (advances) |

---

## 7. Honest scope statement

**2.A KEYSTONE jest LINEARIZED** (`O(h_μν)`), NIE pełną non-perturbative path
integration. Zakres zachowany:

1. **Background M9.1″ fixed** dla derivacji `h_μν` propagatora; pełna
   back-reaction graviton → Φ deferred do **2.F CAPSTONE** (path integral
   `D[h_μν]·D[Φ]`).
2. **`c_T = c_0` nie jest a priori derivacją z czystej geometrii substratu**
   (to było zachowane przez M9.3.5 + OP-7); 2.A daje **EFT-level** kwantyzację
   na **fixed** M9.1″ background.
3. **Faddeev-Popov ghost structure** standardowa (de Donder gauge α=1);
   ghost decoupling pokazany na flat background, nie w pełnym path integral
   (2.F CAPSTONE).
4. **Mass scale `M_phys^TGP ≈ 1.42×10⁻³³ eV` jest hubble-skali** (Φ_0 = H_0);
   2.A dziedziczy 1.A.6 ale nie wprowadza nowej skali.
5. **Vector mode structural zero** wynika z single-Φ axiom (TGP_FOUNDATIONS §1)
   — to jest **architectural choice**, nie derivation z deeper principle;
   B.6/B.7 deeper origin pozostaje OPEN structural postulate.

**2.A NIE ustanawia:**
- pełnej UV-complete renormalizability (research-track, off-scope);
- non-linear graviton self-interaction (`O(h³)` i wyższe — Donoghue 1994
  framework + 2.D analiza);
- back-reaction graviton → Φ-mass running (2.F);
- dispersive corrections z curvature radius `R_curv ≪ k⁻¹` (sub-leading
  na M9.1″ vacuum).

**2.A USTAWIA:**
- 2 TT polaryzacje (h_+, h_×) — standard GR
- 1 scalar mode (h_b = h_L) — TGP single-Φ heritage M9.3.4
- 0 vector modes — single-Φ structural zero
- Total physical DOF = 3 na M9.1″
- GW170817 strict equality (∞ margin) na vacuum
- Inputs dla 2.D (counterterm structure) i 2.F (path integral measure).

---

## 8. Phase 2 status po 2.A

```
Phase 2 (quantum gravity proper / EFT)
├── 2.0    ✅ CLOSED 16/16 PASS         (2026-04-28)
├── 2.A    ✅ CLOSED 6/6 PASS  KEYSTONE (2026-04-28)
├── 2.B    ⏳ PENDING                    first-principles α₀
├── 2.D    ⏳ PENDING                    EFT renormalizability
├── 2.E    ⏳ PENDING                    B.1/B.2/B.5 deepening
├── 2.F    ⏳ PENDING (CAPSTONE)         path integral D[h_μν]   [unblocked]
└── 2.R-final  ⏳ PENDING                synthesis 8 R.F + ≥217 cumulative

Cumulative: 167 (Phase 1) + 16 (2.0) + 6 (2.A) = 189 verifications
Phase 2 baseline target po pełnym zamknięciu: ≥217
Critical path advance: 2.A → 2.F (now unblocked)
```
