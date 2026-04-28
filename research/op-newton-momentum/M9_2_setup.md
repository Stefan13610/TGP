# M9.2 — Pęd, bezwładność, zasada równoważności (setup)

**Cykl:** M9 ("klasyczna dynamika i pęd"), test 2 z 3.
**Data:** 2026-04-26
**Status:** OPEN (kickoff post closure_2026-04-26)
**Solver:** `m9_2_momentum.py` (numpy + scipy ODE/integral)
**Output target:** `m9_2_momentum.txt`, synteza `M9_2_results.md`

**Foundations binding:** [[TGP_FOUNDATIONS.md]] §6 (pęd = Lenz-podobna back-reakcja)
**Predecessor:** [[M9_1_pp_P3_results.md]] (M9.1'' z hiperboliczną metryką, 3/5 PASS + 1 conditional + 1 open)
**Closure_2026-04-26 binding:** [[../closure_2026-04-26/alpha_psi_threshold/results.md]] (T-α α(ψ) baseline dla Candidate D)

---

## 1. Cel M9.2

Po pozytywnym M9.1'' (1PN-PPN PASS dla hiperbolicznej metryki) pozostaje
najgłębsze pytanie programu emergent gravity:

> **Czy klasyczna Φ-EOM TGP odzyskuje pęd, bezwładność i zasadę równoważności
> jako emergentne własności pola, bez postulowania m_inertial niezależnie?**

TGP_FOUNDATIONS §6 postuluje że pęd jest **Lenz-podobnym** efektem
back-reakcji pola Φ na ruchome źródło. M9.2 testuje to obliczeniowo.

**Pytanie strategiczne:**
> Dla statycznego źródła `ρ_0(x)` + małego stałego przyspieszenia `a`,
> czy back-reakcja pola na źródło daje siłę `F_back = -m_b·a` z
> `m_b = m_g` (zasada równoważności emergentna z field theory)?

---

## 2. Setup matematyczny

### 2.1 Φ-EOM w slow-source limit

Z M9.1'' baseline (hiperboliczna metryka):
```
∇²ε + 2(∇ε)²/(1+ε) - β(1+ε)²·ε = -q·ρ(r,t)        (Φ-EOM, ε = φ - 1)
```

Dla ruchomego źródła `ρ(x,t) = ρ_0(x - X(t))` z `X(t) = (1/2) a t²`:
- W spoczynkowej ramie źródła: pole quasi-statyczne `ε_0(x') ≈ ε_eq(|x'|)`.
- W ramie laboratorium: `ε(x,t) = ε_0(x - X(t)) + δε(x,t)` z `δε` retardowane.

### 2.2 Linearyzacja w przyspieszeniu a

Rozwijamy w potędze `a`:
```
ε(x,t) = ε_0(x - X(t)) + a·ε_1(x,t) + O(a²)
```

`ε_0` jest statycznym M9.1'' rozwiązaniem (Yukawa-like z β > 0).
`ε_1` zawiera retardację (skończona prędkość propagacji `c_0 = 1` w units).

Wstawiając do Φ-EOM i biorąc człon liniowy w `a`:
```
□ε_1 - β·ε_1 = -∂_t² ε_0 / a    (źródło z opóźnienia ruchu)
```

`∂_t² ε_0 = (a·∇)² ε_eq + a·∇(∂_t ε_eq)|_{a-source}` → daje "back-flow"
gradient field który odpowiada sile bezwładności.

### 2.3 Back-reakcja: siła pola na źródło

Siła pola na źródło:
```
F_back = -∫ ρ(x,t) ∇φ(x,t) d³x
       = F_static + a·F_inertial + O(a²)
```

W spoczynku (`a = 0`): `F_static = -∫ ρ_0 ∇ε_eq d³x = 0` (równowaga statyczna).

Pierwsze poprawki (linear in `a`):
```
F_inertial = -∫ ρ_0 ∇ε_1 d³x · a + ∫ (a·∇ρ_0) ε_eq d³x
```

**Definicja masy bezwładnościowej:** `F_back = -m_b · a`, więc
```
m_b := ∫ ρ_0 ∇ε_1 d³x · ê_a + ∫ (ê_a·∇ρ_0)·ε_eq d³x        (kanoniczna)
```

### 2.4 Masa grawitacyjna

Z minimalnego sprzężenia `S_int = -(q/Φ_0) ∫ φ ρ d⁴x`:
```
m_g := q · ∫ ρ_0 d³x       (do c_0² factors)
```

W jednostkach gdzie `q = 8πG/c_0²`:
```
m_g = (8πG/c_0²) · M_total       (M_total = ∫ ρ_0 d³x)
```

### 2.5 Test zasady równoważności (WEP)

Predykcja TGP:
```
m_b / m_g = 1 + O(ε²)        (η_TGP_intrinsic ~ U² < 10⁻¹⁵)
```

Dla Earth-scale (U ~ 7×10⁻¹⁰), `O(ε²) ~ 5×10⁻¹⁹` → super-safe pod
MICROSCOPE bound η < 1.1×10⁻¹⁵.

**Connection do closure_2026-04-26 Phase 4:** T-α α(ψ) z thresholdem przy
ψ=1 zachowuje WEP MICROSCOPE-safe **niezależnie** od strong-field
α-coupling przy photon ring. M9.2 testuje **field-theoretical baseline**
(without α(ψ) coupling, czyli α=0 w lab).

### 2.6 Energia radiowana

Larmor-like formula dla scalar field:
```
dE/dt|_radiated = (q²/24π c_0³) · ⟨ä⟩²
```

(scalar analog; dla GW wymaga kwadrupolu, M9.3)

Statyczne `a = const` → brak emisji (bo `ä = 0`); test "Newton I bez radiation".

---

## 3. Plan testów (M9.2.1 ... M9.2.5)

| ID | Cel | Metoda | PASS criterion |
|----|-----|--------|----------------|
| **M9.2.1** | Statyczne BC: `F_back(a=0) = 0` | numerical integration | ‖F_back‖ < 10⁻¹⁰ (Newton's I) |
| **M9.2.2** | Linearna back-reakcja `F = -m_b·a` | scan a ∈ [10⁻⁶, 10⁻³] | linear fit R² > 0.999 |
| **M9.2.3** | WEP: `m_b = m_g` w słabym polu | computed ratio | \|m_b/m_g - 1\| < 10⁻³ |
| **M9.2.4** | Skalowanie z masą źródła | scan M ∈ [10⁻³, 10³] | m_b/M_total = const |
| **M9.2.5** | Stały ruch (a=0) → brak radiation | dE/dt computed | < 10⁻¹⁰ × Larmor unit |

**Sukces:** 5/5 PASS → M9.2 zamyka pęd/bezwładność na poziomie classical
field theory; emergent inertia + WEP **derived**, nie postulowane.

---

## 4. Numeryczna realizacja

### 4.1 Geometria
- Sferycznie symetryczna, `ρ_0(r) = (M/(σ³(2π)^{3/2})) · exp(-r²/(2σ²))`
  (Gaussian source, σ = 1 unit length)
- Dziedzina: `r ∈ [0, R_max]` z `R_max = 50σ` (asymptotyczna płaskość)
- Siatka: 2000 punktów logarytmicznie zagęszczonych

### 4.2 Statyczne pole ε_eq
- BVP solver `scipy.integrate.solve_bvp` z M9.1''
- Hiperboliczna metryka `f(ψ) = (4-3ψ)/ψ` z T-FP closure_2026-04-26
- Yukawa screening dla β > 0

### 4.3 Pierwszy rząd ε_1
- Solver dla `□ε_1 - β·ε_1 = source[ε_0, a]`
- Source z `a·∂_x ε_eq` (dla ruchu w x-direction)
- BC: `ε_1(∞) = 0`, regularność w r=0

### 4.4 Back-reakcja
- `F_inertial = ∫ ρ_0 ∇ε_1 d³x` (numerical Gaussian quadrature)
- Skan `a` → linear fit `F = -m_b · a`

---

## 5. Falsifiable predictions

**M9.2 predykcje (PASS dla TGP):**

1. **WEP MICROSCOPE-2 (η < 10⁻¹⁷):** TGP M9.2 prediction `η_TGP ~ U² ~ 10⁻¹⁹`.
   Combined z closure_2026-04-26 T-α suppression: `η_TGP ~ 10⁻³²` ⇒ strong PASS.

2. **Newton I (constant velocity → no radiation):** TGP `dE/dt|_{a=0} = 0` exactly.
   Falsification: detection of "boost-induced radiation" w precision tests.

3. **m_b = m_g do O(ε²):** TGP composition-independent equivalence.
   Falsification: composition-dependent η at 10⁻¹⁹ level.

---

## 6. Status: structural test po closure_2026-04-26

M9.2 jest **NATURAL FOLLOW-UP** po closure_2026-04-26:
- Phase 1 (Path B σ_ab) — gotowe σ_ab dla M9.3 GW polarizations
- Phase 2 (T-FP f(ψ)) — używamy `f(ψ) = (4-3ψ)/ψ` z M9.1''
- Phase 3 (T-Λ Ω_Λ) — vacuum baseline ε_eq → Φ_eq scale
- Phase 4 (T-α α(ψ)) — α(ψ_lab) ≈ 0 (Earth surface) → M9.2 jest **bezpośrednio**
  weak-field test bez α-induced effects

**Status M9.2:** structural numerical test classical TGP field theory;
zamyka klasyczną fenomenologię gravity (po M9.1'') i otwiera drogę do M9.3 (GW).

---

## 7. Pliki

- `M9_2_setup.md` (ten plik)
- `m9_2_momentum.py` — numpy + scipy script M9.2.1..M9.2.5
- `m9_2_momentum.txt` — raw output
- `M9_2_results.md` — synteza po wykonaniu

---

## 8. Cross-references

- [[M9_program.md]] §4 (M9.2 plan)
- [[M9_1_pp_P3_results.md]] (M9.1'' baseline)
- [[../closure_2026-04-26/alpha_psi_threshold/results.md]] (T-α connection)
- [[../closure_2026-04-26/f_psi_principle/results.md]] (f(ψ) principle)
- [[TGP_FOUNDATIONS.md]] §6 (pęd jako Lenz back-reakcja)
- [[../op-m92/OP_M92_P0plus_candD_results.md]] (Candidate D = momentum back-reaction)
