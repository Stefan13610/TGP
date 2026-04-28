# M9.3 — Fale grawitacyjne, polaryzacje, prędkość propagacji (setup)

**Cykl:** M9 ("klasyczna dynamika i pęd"), test 3 z 3.
**Data:** 2026-04-26
**Status:** OPEN (kickoff post M9.2 closure 5/5 PASS)
**Solver:** `m9_3_gw.py` (numpy + scipy ODE/PDE + linearyzowane δΦ-EOM)
**Output target:** `m9_3_gw.txt`, synteza `M9_3_results.md`

**Foundations binding:** [[TGP_FOUNDATIONS.md]] §5 (gravity = collective Φ fluctuations, brak grawitonu)
**Predecessor:** [[M9_2_results.md]] (5/5 PASS, m_field = ∫(∇ε)² zweryfikowane)
**Path B binding:** [[../closure_2026-04-26/sigma_ab_pathB/results.md]] (σ_ab = ⟨(∂_aŝ)(∂_bŝ)⟩^TF, m_σ² = 2m_s²)
**T-α connection:** [[../closure_2026-04-26/alpha_psi_threshold/results.md]] (α(ψ_NS) ~ 0.65 → measurable phase shift)

---

## 1. Cel M9.3

Po M9.1'' (PPN β=γ=1) i M9.2 (pęd/inercja Lenz-back-reakcja) ostatnim
krokiem klasycznej fenomenologii grawitacyjnej TGP jest **promieniowanie
fal grawitacyjnych** z dynamicznych źródeł.

> **Pytanie centralne:**
> Czy dynamiczna δρ (oscylujący kwadrupol) generuje przez Φ-EOM falowy
> wzór `δΦ + δσ_ab` zgodny z formułą kwadrupolową Peters-Mathews,
> z `c_GW = c₀` (GW170817-safe) i polaryzacjami w granicach LIGO scalar-mode bound?

**Drugorzędne pytanie strukturalne (motywowane intuicją ciśnienie/napięcie z M9.2):**

> **Czy mod skalarny δΦ (analog fali ciśnienia w membranie Φ) i mod tensorowy
> σ_ab (analog fali napięcia z Path B) mają tę samą prędkość propagacji,
> czy różną z powodu m_σ² = 2 m_s²?**

Pierwsza odpowiedź zamyka cykl M9 (compatibility z LIGO/Virgo). Druga
otwiera **nowe testowalne predykcje** TGP odróżniające go od GR.

---

## 2. Intuicja fizyczna: ciśnienie/napięcie membrany Φ

### 2.1 Mapowanie kanoniczne (post-M9.2 dyskusja)

| Obiekt fizyczny | Analog w membranie | TGP form |
|---|---|---|
| **Próżnia** | Naprężona membrana w równowadze | ψ=1, V(1)=β/12 > 0 |
| **Mass M** | Lokalne wybrzuszenie membrany | ε(r) > 0, ψ > 1 |
| **Pęd/inercja** | Opór membrany przeciw deformacji | M9.2 m_field = ∫(∇ε)² ✅ |
| **Mod skalarny δΦ** | Fala ciśnienia (longitudinalna) | δΦ(x,t), m_s² = β |
| **Mod tensorowy σ_ab** | Fala napięcia (transversalna) | σ_ab kompozyt, m_σ² = 2m_s² |
| **GW propagation** | Membrana drgająca jak całość | linearyzacja around Φ_eq |
| **Prędkość fali** | c_s (compression) vs c_T (shear) | **może różne** w TGP |

### 2.2 Strukturalne uzasadnienie różnicy c_s vs c_T

W ośrodku elastycznym (analog mechaniczny) fale podłużne (P-waves)
i poprzeczne (S-waves) mają **różne** prędkości:
```
v_P² = (K + 4G/3)/ρ        (compression)
v_S² = G/ρ                 (shear)
```
gdzie K = bulk modulus, G = shear modulus, ρ = density.

W TGP analog:
- **δΦ (skalarne):** sondują compression-like response Φ — czują **bulk modulus** ~ V''(φ_eq) + K(φ_eq)·k²
- **σ_ab (tensorowe):** sondują shear response substratu przez kompozyt ⟨(∂ŝ)(∂ŝ)⟩^TF — czują **shear modulus** related to ⟨ŝ⟩-dynamics

**Path B closure_2026-04-26 dał kluczową relację:**
```
m_σ² = 2 m_s²        (m_s = effective mass of fundamental ŝ fluctuations)
```

To **bezpośrednio wymusza** różnicę dyspersji:
```
ω_s²(k) = c₀²k² + m_s²c₀⁴/ℏ²       (mod skalarny)
ω_T²(k) = c₀²k² + m_σ²c₀⁴/ℏ² = c₀²k² + 2m_s²c₀⁴/ℏ²   (mod tensorowy)
```

### 2.3 Implikacje obserwacyjne

- **Wysokie k (LIGO ~100Hz):** ω/k → c₀ dla obu, różnica < 10⁻¹⁵ (GW170817-safe)
- **Niskie k (PTA mHz, LISA):** różnica może być mierzalna w fazie inspiral
- **NS-NS ringdown:** scalar i tensor mode mogą "dzwonić" przy różnych częstotliwościach

To jest **nowa falsyfikowalna predykcja TGP**, motivated przez intuicję
ciśnienie/napięcie + Path B m_σ²=2m_s².

---

## 3. Setup matematyczny

### 3.1 Linearyzacja Φ-EOM wokół Φ_eq[ρ_static]

Z Φ-EOM (FOUNDATIONS §3, sek08a `eq:field-eq-reproduced`):
```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀² = -qΦ₀·ρ
```

Pivot:
```
Φ(x,t) = Φ_eq(x) + δΦ(x,t),      |δΦ| ≪ Φ_eq
ρ(x,t) = ρ_static(x) + δρ(x,t)
```

Statyczny baseline `Φ_eq[ρ_static]` z M9.1''-P3 (hiperboliczna metryka).

**Liniowy człon w δΦ:**
```
□_eff δΦ + M²(Φ_eq) · δΦ = -q Φ₀ · δρ                    (δΦ-EOM)
```

gdzie:
- `□_eff = g_eff^μν[Φ_eq] ∂_μ ∂_ν` — d'Alembertian na statycznej metryce hiperbolicznej
- `M²(Φ_eq) = β(2 - 3φ_eq)·φ_eq + nonlinear from kinetic D_kin'`

**W limicie weak-field** (φ_eq → 1): `M²(1) = -β`, ale uwzględniając
kinetic correction (z ∂_μ[K(φ)∂^μφ]):
```
M²_eff(φ_eq=1) = β > 0       (Yukawa, weak-field; potwierdzone M9.1, M9.2)
```

### 3.2 Mod tensorowy z σ_ab (Path B PRIMARY)

Z OP-7 T3 closure_2026-04-26:
```
σ_ab(x,t) = ⟨(∂_a δŝ)(∂_b δŝ)⟩^TF       (kompozyt z poziomu 0)
```

gdzie `δŝ` to fluktuacje substratu wokół ⟨ŝ⟩, ŝ-EOM:
```
□_eff δŝ + m_s² · δŝ = J_ŝ[δρ]            (linearyzowane ŝ-EOM)
```

`σ_ab` jako bilinear `δŝ` ma dyspersję:
```
m_σ² = 2 m_s²        (Path B PRIMARY result)
□_eff σ_ab + m_σ² · σ_ab = T_ab^TF[δρ]   (effective σ_ab-EOM)
```

`T_ab^TF` = trace-free part of stress-energy tensor (kwadrupol moment of δρ).

### 3.3 Metryka efektywna z perturbacji

Z FOUNDATIONS §3 hiperbolicznej metryki:
```
g_tt^eff[Φ] = -c₀²·(4 - 3ψ)/ψ,    g_ii^eff[Φ] = ψ/(4 - 3ψ)
```

Perturbacja:
```
δg_tt = ∂g_tt/∂ψ · (δΦ/Φ₀) = c₀²·(4/ψ²)·(δΦ/Φ₀)|_{ψ=1} = 4c₀²·(δΦ/Φ₀)
δg_ii = ∂g_ii/∂ψ · (δΦ/Φ₀) ≈ -(δΦ/Φ₀)
```

**Tensorowa składowa** wchodzi przez modyfikację metryki pochodzącą z σ_ab
(Path B mechanism):
```
δg_ab^TT = ξ · σ_ab          (ξ = Path B coupling, OP-7 T3 result)
```

### 3.4 Polaryzacje: dekompozycja SO(3) wokół kierunku propagacji

Dla GW propagującej w kierunku ẑ, dekomponujemy `δg_μν` na nieprzywiedne reprezentacje:

| Polaryzacja | Symbol | Źródło w TGP | Status |
|---|---|---|---|
| Plus | h₊ | σ_ab^TT (Path B) | tensorowa |
| Cross | h_× | σ_ab^TT (Path B) | tensorowa |
| Breathing (scalar) | h_b | δΦ (skalarna) | **TGP-specific** |
| Longitudinal (scalar) | h_L | δΦ (subleading) | **TGP-specific** |
| Vector x | h_x | brak | **0 w TGP** |
| Vector y | h_y | brak | **0 w TGP** |

**Klucz strukturalny:** TGP ma **jedno** pole Φ → brak modu wektorowego
(brak struktury 3-wektorowej dla GW). LIGO bound: vector-mode amplitude
< few % — TGP **strukturalnie** spełnia.

Skalarne mody (h_b, h_L) mogą być **ograniczone przez LIGO** (typowo
< 5-10% scalar-tensor ratio) — **falsyfikowalne**.

### 3.5 Energia radiowana: kwadrupol Peters-Mathews + scalar mode

**Mod tensorowy** (analog Einsteina):
```
dE/dt|_tensor = (G/5c⁵) · ⟨Q̈_ij Q̈^ij⟩_{TT}
```
gdzie `Q_ij` = traceless quadrupole moment of δρ.

**Mod skalarny** (analog Larmora dla scalar field):
```
dE/dt|_scalar = (q²/12π c₀³) · ⟨Q̈⟩²       (monopol = 0 z conservation)
              ≈ (q²/60π c₀⁵) · ⟨Q̈_kk⟩²   (kwadrupol scalar trace)
```

W praktyce dla TGP scalar mode jest **subleading** ze względu na
Yukawa screening β > 0 (suppression przez factor `e^(-r√β)` na skalach
> 1/√β).

---

## 4. Plan testów (M9.3.1 ... M9.3.5)

| ID | Cel | Metoda | PASS criterion |
|----|-----|--------|----------------|
| **M9.3.1** | Linearyzacja sanity: □_eff well-defined on hyperbolic Φ_eq | algebraiczna weryfikacja + numerycznie M²(φ=1) = β | M²_eff > 0; m_σ² = 2m_s² confirmed |
| **M9.3.2** | Far-field wave: c_T = c_0, c_s = c_0 (high k) | dispersion relation ω(k) z δΦ-EOM, σ_ab-EOM | \|c_GW - c\|/c < 10⁻¹⁵ at LIGO freq |
| **M9.3.3** | Kwadrupol Peters-Mathews | dE/dt analytycznie + numerycznie dla binary | matching < 10⁻¹ relative to GR |
| **M9.3.4** | Polaryzacje SO(3) decomposition | h₊, h_×, h_b, h_L amplitudes at infinity | h_scalar/h_tensor < LIGO bound (~10%) |
| **M9.3.5** | GW170817 compatibility + dispersion check | check c_T=c_s przy LIGO freq | strong PASS (>10⁻¹³ margin) |

**Sukces:** 5/5 PASS → M9.3 zamyka klasyczną fenomenologię GW TGP;
falsifiable predictions for low-freq (PTA/LISA) and NS-NS ringdown.

**Bonus (jeśli czas):** **M9.3.6** — phase shift NS-NS ringdown z α(ψ_NS)
coupling (T-α z closure_2026-04-26 Phase 4): `α(ψ_NS=1.4) = α₀(0.4)² ≈ 0.65`
→ **measurable signature in LIGO O5 NS-NS template residuals**.

---

## 5. Numeryczna realizacja

### 5.1 Statyczny baseline Φ_eq

- Hiperboliczna metryka z M9.1''-P3, `f(ψ) = (4-3ψ)/ψ` (T-FP closure)
- Spherical source ρ_static = Gaussian (M9.2.1 baseline reused)
- Domain: r ∈ [0, R_max=50σ], log-zagęszczona siatka 2000 pkt
- BVP solver `scipy.integrate.solve_bvp` (M9.2 baseline)

### 5.2 Linearyzacja δΦ-EOM (czasowa propagacja)

- Spherical-harmonic decomposition: `δΦ(r,θ,φ,t) = Σ_lm δΦ_lm(r,t)·Y_lm(θ,φ)`
- Radial PDE: `(∂_t² + ω_l² + Veff_l(r)) δΦ_lm = -qΦ₀·δρ_lm(r,t)`
- Source: oscylujący kwadrupol `δρ ~ ε·sin(ωt)·Y_2m(θ,φ)·ρ_static(r)`
- BC: outgoing wave w `r=R_max` (Sommerfeld); regular w r=0
- Time-stepping: 4th-order Runge-Kutta lub leapfrog

### 5.3 Tensorowa składowa σ_ab (Path B linearyzacja)

- δŝ-EOM linearization (z OP-7 T3 framework)
- Box-of-product algebra: `σ_ab(k) = ∫ d⁴k' (∂_a)(∂_b) ⟨δŝ(k') δŝ(k-k')⟩^TF`
- Effective dispersion: ω² = c₀²k² + 2m_s²c₀⁴/ℏ²
- Source: T_ab^TF[δρ] (trace-free quadrupole stress)

### 5.4 Far-field extraction

- Wave amplitude h_+, h_×, h_b, h_L extracted at r = R_extract = 30σ
- Fourier transform → ω-spectrum
- Polarization decomposition relative to ẑ (propagation direction)
- Energy flux: `dE/dt = ∫_{r=R_extract} S_r dΩ` (radial Poynting-like)

### 5.5 Kwadrupol Peters-Mathews comparison

- Test binary: two point masses M, separation a, ω_orbit = √(GM_tot/a³)
- Predicted GR: `dE/dt|_GR = (32/5)·G⁴M₁²M₂²(M₁+M₂)/(c⁵a⁵)`
- TGP: `dE/dt|_TGP = (32/5)·G⁴M₁²M₂²(M₁+M₂)/(c⁵a⁵) · (1 + δ_TGP)`
- PASS: `|δ_TGP| < 10⁻¹` w nominal regime (testowalne LIGO 3G/ET)

---

## 6. Falsifiable predictions M9.3

### 6.1 Compatibility predictions (PASS = standard GR-zgodność)

1. **GW170817 (`|c_GW - c|/c < 10⁻¹⁵`):** TGP predicts equality strukturalnie
   z konstrukcji g_eff (√-g_eff = c₀φ; metric is algebraic in Φ).

2. **LIGO scalar-mode bound:** TGP scalar mode (h_b) suppressed przez Yukawa
   screening; predicted ratio `h_scalar/h_tensor ~ exp(-r√β)` przy LIGO band.

3. **Peters-Mathews quadrupole:** TGP reproduces GR formula at leading order,
   z ξ-dependent correction `~ (m_σ/M_Pl)² ≪ 1`.

### 6.2 New TGP-specific falsifiable predictions

1. **Different dispersion modes (LISA / PTA / 3G):** TGP predicts
   ```
   ω_s²(k) - ω_T²(k) = (m_s² - m_σ²)c₀⁴/ℏ² = -m_s²c₀⁴/ℏ²
   ```
   At very low k (mHz range), measurable phase difference between scalar and tensor channels in NS-NS / BBH inspirals.
   **Falsification:** equal dispersion within LISA precision falsyfikuje m_σ² = 2m_s² (Path B).

2. **NS-NS ringdown phase shift z α(ψ_NS) (T-α from closure):**
   ```
   α(ψ_NS=1.4) = α₀·(0.4)² ≈ 0.65
   ```
   Predicted phase shift `Δφ ~ α·v²/c²` in LIGO O5 templates relative to GR.
   **Falsification:** absence of phase shift falsifies T-α + M9.2-D combined.

3. **No vector mode (strict zero):** TGP single Φ → vector-mode amplitude
   = 0 strukturalnie. **Falsification:** detection of vector mode in
   LIGO/Virgo/KAGRA falsifies single-Φ aksjomat (FOUNDATIONS §1).

---

## 7. Status: zamknięcie cyklu M9

M9.3 jest **finałem klasycznej fenomenologii grawitacyjnej TGP**:

| Test | M9.x | Status |
|---|---|---|
| Newton II (F = ma) | M9.0 | ✅ |
| 1PN-PPN (β, γ = 1) | M9.1'' | ✅ (3 PASS + 1 cond + 1 open) |
| Pęd / inercja Lenz | M9.2.2 | ✅ |
| WEP MICROSCOPE-safe | M9.2.3 | ✅ (η_TGP ~ 10⁻⁵⁰ ≪ 10⁻¹⁵) |
| Mass scaling m∝M² | M9.2.4 | ✅ |
| Newton I (no rad) | M9.2.5 | ✅ |
| **GW c_GW = c** | **M9.3.5** | **⏳** |
| **GW kwadrupol** | **M9.3.3** | **⏳** |
| **GW polaryzacje** | **M9.3.4** | **⏳** |
| **Dispersion c_s vs c_T** | **M9.3.2** | **⏳** (new prediction) |

**Sukces M9.3 → klasyczna grawitacja TGP ZAMKNIĘTA**, gotowa do publikacji
i testowania observacyjnego w LIGO O5 (2027+) i LISA (2035+).

---

## 8. Pliki

- `M9_3_setup.md` (ten plik)
- `m9_3_gw.py` — numpy + scipy + linearization solver, 5 testów
- `m9_3_gw.txt` — raw output
- `M9_3_results.md` — synteza po wykonaniu

---

## 9. Cross-references

- [[M9_program.md]] §5 (M9.3 plan)
- [[M9_2_results.md]] (M9.2 5/5 PASS, m_field = ∫(∇ε)² baseline)
- [[M9_1_pp_P3_results.md]] (M9.1'' hyperbolic metric baseline)
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] (Path B PRIMARY, m_σ² = 2m_s²)
- [[../closure_2026-04-26/alpha_psi_threshold/results.md]] (T-α α(ψ) → NS-NS phase shift)
- [[../closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]] (4 closure context)
- [[TGP_FOUNDATIONS.md]] §5 (gravity = collective Φ fluctuations, brak grawitonu)
- [[TGP_FOUNDATIONS.md]] §2 (4-poziomowa hierarchia, σ_ab kompozyt poziomu 0)

---

## 10. Hipoteza guidance: ciśnienie/napięcie membrany Φ jako intuicja

Z dyskusji 2026-04-26 (post M9.2 closure):

> **Centralny wgląd:**
> Próżnia TGP ma napięcie strukturalne (V(ψ=1) = β/12 > 0, T-Λ closure).
> Mass tworzy lokalne wybrzuszenie membrany (ε > 0 z Yukawy). Otaczające
> "napięcie" próżni wymusza powrót do równowagi → emergentna grawitacja
> jako manifestation napięcia, nie "zasysania".
>
> Path B m_σ² = 2 m_s² implikuje że mod skalarny (ciśnienie) i mod
> tensorowy (napięcie poprzeczne) propagują z **różnymi dyspersjami**
> przy niskich k. To jest natywna predykcja TGP testowalna w
> LISA/PTA/3G era.

**Ta hipoteza będzie:**
- Strukturalnym zwornikiem M9.3.2 (dispersion test)
- Inspiracją M9.3.4 (polaryzacja jako "P-wave + S-wave membrany Φ")
- Mostem do M9.3.6 NS-NS ringdown (jeśli zdąży)

**Drift check (post FOUNDATIONS):** intuicja w zgodzie z §5.2
("kolektywny efekt fluktuacji"), §6.1 (Lenz-back-reakcja jako
"opór membrany"), Path B σ_ab kompozyt (nie nowe pole). ✅

---

**Status M9.3:** structural test classical GW phenomenology TGP;
zamyka cykl M9 i otwiera observational falsification matrix
(LIGO O5, LISA, PTA, ngEHT, MICROSCOPE-2).
