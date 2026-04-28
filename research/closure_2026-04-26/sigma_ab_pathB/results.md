# T-PB results — σ_ab Path B audit (POSITIVE 11/11)

**Data:** 2026-04-26
**Status:** ✅ POSITIVE
**Plik wykonawczy:** `sigma_ab_pathB_audit.py`
**Raw output:** `sigma_ab_pathB_audit.txt`
**Cross-references:**
- [[research/op7/OP7_T3_results.md]] (T3 STRUCTURAL POSITIVE z Φ₀/m_σ tension)
- [[research/op7/OP7_T2_results.md]] (σ_ab gradient strain definition)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.1 (brainstorm motivating Path B promotion)
- [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂ axiom)
- [[setup.md]] (audyt design)

---

## TL;DR

> Path B (heredity) jest teraz **primary derivation** σ_ab dynamiki.
>
> Wszystkie 5 kryteriów audytu PASS:
> 1. **T-PB.1** — heredity equation `□K_ab + 2 m_s² K_ab = source + grad-coupling`
>    derivowane krok-po-kroku z `□δŝ + m_s² δŝ = J` (sympy exact zero residual).
> 2. **T-PB.2** — `M² = 2 m_s²` derived, nie postulated;
>    spectral threshold `√s_min = 2 m_s` z OPE.
> 3. **T-PB.3** — ghost-free z konstrukcji (Gram-matrix property K_ab ≥ 0).
> 4. **T-PB.4** — żaden nowy d.o.f. (single-Φ axiom *jawnie* zachowany).
> 5. **T-PB.5** — σ_ab = 0 dla static spherical (analitycznie + numerycznie 1e-18).
>
> **Implikacja:** TGP nigdy nie wprowadza Lagrangianu L_σ jako fundamentalny;
> Path A (efektywny opis) używamy tylko jako **convenient computational frame**
> dla quadrupole/PPN, ale ontologicznie wszystko żyje w jednym ŝ-polu.

---

## 1. Co się zmieniło w stosunku do OP-7 T3

OP-7 T3.1 zamknął się ze stwierdzeniem:
> "Path A == Path B structurally; pojedyncze-Φ aksjomat zachowany."

ALE:
- Path A (Lagrangian L_σ z postulowanym m_σ²) służył jako **primary** dla
  T3.2-T3.4 obliczeń.
- Identyfikacja `m_σ² = 2 m_s²` była *post hoc analogia mezonowa*,
  nie wynik bezpośredniej derivacji.

T-PB **odwraca** tę kolejność:
- Path B jest primary derivation: σ_ab dynamika to **operatorowa konsekwencja**
  ŝ-EOM, nie nowa postulat.
- Path A pozostaje jako **effective theory** dla wygody rachunków,
  ale jego współczynniki są wyznaczone (nie wolne parameters):
  - `M_σ² = 2 m_s²` — z box-of-product algebry.
  - `ξ_eff = G·Φ₀² · O(1)` — z T3.4 matching.
  - `Λ(ψ) = 1` — z T6 closure (decoupling regime).

---

## 2. Wyniki testów (11/11 PASS)

### T-PB.1 — Heredity equation z box[bilinear] algebry

**Symboliczna derivacja:**

Start: `□δŝ + m_s² δŝ = J`.

Apply `□` do `K_ab = (∂_a δŝ)(∂_b δŝ)`:
```
□[K_ab] = (□∂_a δŝ)(∂_b δŝ) + 2(∂_μ ∂_a δŝ)(∂^μ ∂_b δŝ) + (∂_a δŝ)(□∂_b δŝ)
```

Substitute `□∂_a δŝ = ∂_a(J − m_s² δŝ) = ∂_a J − m_s² ∂_a δŝ`:
```
□[K_ab] = −2 m_s² K_ab + (∂_a J)(∂_b δŝ) + (∂_a δŝ)(∂_b J) + 2(∂_μ ∂_a δŝ)(∂^μ ∂_b δŝ)
```

Move mass term:
```
□[K_ab] + 2 m_s² K_ab = (∂_a J)(∂_b δŝ) + (∂_a δŝ)(∂_b J) + 2(∂_μ ∂_a δŝ)(∂^μ ∂_b δŝ)
```

Sympy: `LHS_substituted − reconstruction = 0` (exact zero residual).

**Po coarse-grain `⟨...⟩_B` i odjęciu śladu:**
```
□σ_ab + (2 m_s²) σ_ab = S_ab^TT[J,δŝ] + R_ab^TT[higher-OPE]    (heredity EOM)
```

To jest **derived equation**, nie postulat.

| Test | Wynik |
|------|-------|
| T-PB.1a Heredity-form decomposition exact | PASS |
| T-PB.1b M² = 2 m_s² coefficient automatic | PASS |

### T-PB.2 — Composite mass identification

**Heredity EOM coefficient:** `M² = 2 m_s²` (z box-of-product).

**Spectral threshold (OPE):** dla composite operator O(x) = (∂_a δŝ)(∂_b δŝ):
- spectral function `ρ_O(p²)` ma support na `p² ≥ (2 m_s)²`
- dolny kraniec continuum: `√s_min = 2 m_s`

Dwie różne skale — obie *derived*:
- `M² = 2 m_s²` — coefficient w heredity EOM (relevant dla dispersion).
- `s_min = 4 m_s²` — threshold continuum (relevant dla two-particle decay).

**Decoupling regime check** (OP-7 T6, m_s ~ 0.5 meV):
```
M_eff = √2 · m_s ≈ 0.71 meV
ω_LIGO ~ 1e-13 eV
M_eff / ω_LIGO ≈ 7e+9 ≫ 1   → effective masslessness, c_GW = c₀
```

| Test | Wynik |
|------|-------|
| T-PB.2a M² = 2 m_s² derived | PASS |
| T-PB.2b √s_min = 2 m_s | PASS |
| T-PB.2c decoupling M_eff ≫ ω_LIGO | PASS (ratio 7e9) |

### T-PB.3 — Ghost-free z konstrukcji

**Substrate Hamilton density:**
```
H_s = (1/2)π_s² + (1/2)|∇δŝ|² + (1/2) m_s² δŝ²    (≥ 0)
```

Hessian eigenvalues: [1, 1, 1, 1, m_s²] — wszystkie dodatnie.

**K_ab jako Gram matrix:** dla dowolnego pola wektorowego `v_a = ∂_a δŝ`,
correlator `⟨v_a v_b⟩` jest macierzą Grama → eigenvalues ≥ 0.

Konsekwencja: **żaden ghost mode nie może się pojawić** w derived
σ_ab dynamice, bo cała kinetyka pochodzi z dodatnio-określonego H_ŝ.

| Test | Wynik |
|------|-------|
| T-PB.3a H_s ≥ 0 | PASS |
| T-PB.3b K_ab Gram-positivity | PASS |

### T-PB.4 — Degrees-of-freedom counting

**Phase-space substratu:**
- 1 scalar field `δŝ` z conjugate momentum `π_s = ∂_t δŝ`.
- 1 propagating canonical pair → 1 d.o.f.

**σ_ab(x):**
- 5 algebraic components (symmetric traceless in 3D).
- ALE wszystkie określone przez `∂_a δŝ` + averaging procedure.
- Brak nowego canonical pair `(Π_σ, σ)`.
- **σ_ab nie jest niezależnym d.o.f.** — to operator pochodny.

| Test | Wynik |
|------|-------|
| T-PB.4a No new canonical d.o.f. | PASS |
| T-PB.4b Single-Φ axiom preserved | PASS |

### T-PB.5 — Static spherical reduction

**Analitycznie:**
- W static sphericznej config: `δŝ = δŝ(r)`.
- `∂_a δŝ = (dδŝ/dr)(x_a/r)`.
- Z izotropowego averaging: `⟨x_a x_b⟩ = (B²/3) δ_ab`.
- `K_ab = (1/3)(dδŝ/dr)² δ_ab` → `Tr K = (dδŝ/dr)²`.
- `σ_ab = K_ab − (1/3) δ_ab Tr K = 0` **exactly**.

**Numerycznie (N=32 grid, Gaussian profile):**
- `||σ||/Tr K = 2.20e-18` (machine precision).

**Implikacja:** M9.1'' P3 PPN audit (Mercury, Cassini, LLR) **nie jest zaburzony**
przez Path B heredity, bo σ_ab znika w statycznych sferycznych config — zgodnie
z OP-7 T1 no-tensor wynikiem.

| Test | Wynik |
|------|-------|
| T-PB.5a σ_ab = 0 analitycznie | PASS |
| T-PB.5b σ_ab = 0 numerycznie (1e-18) | PASS |

---

## 3. Konsekwencje strategiczne

### 3.1 Spirit-fidelity TGP_FOUNDATIONS §1 jawnie potwierdzona

Single-Φ Z₂ axiom nie jest tylko "structurally OK" (T2 verdict),
ale **operacyjnie redundant** — żadna część TGP dynamiki nie wymaga
niczego więcej poza ŝ-EOM. Path A pozostaje *syntactic sugar* dla
T3.4 quadrupole rachunku, ale nie jest required.

### 3.2 m_σ² staje się prediction, nie założeniem

W Path A: `m_σ²` byłaby free parameter dopasowywaną do danych.
W Path B: `m_σ² = 2 m_s²` — wartość wymuszona.

**Falsyfikowalność wzmocniona:**
- LIGO 3G (Cosmic Explorer ~2030) pomiar dispersion GW.
- Jeśli detekcja `M_GW² ≠ 2 m_s²` (z m_s zmierzone osobno przez ULDM/cosmology),
  TGP single-Φ jest sfalsyfikowany na poziomie tensor sector dynamics.
- Jeśli detekcja zgodna z `M_GW² = 2 m_s²` (lub obie pomijalne) → confirmation.

### 3.3 Decoupling regime jednoznacznie określony

OP-7 T6 closure roboczo przyjął m_σ ~ meV jako "decoupling scale".
Path B precyzuje: `M_eff = √2 · m_s`, więc decoupling regime
jest *automatic* gdy `m_s ≫ ω_observation`.

Dla LIGO band (10 Hz – 10 kHz):
- ω_max ~ 4e-11 eV << m_s ~ 0.5 meV.
- Ratio ~ 1e7 → c_GW = c₀ exact do akceptowalnej precyzji LIGO.

### 3.4 Φ₀/m_σ tension częściowo rozwiązana

OP-7 T3.2 wskazał: m_σ ~ Φ₀ z naturalnej skali → tension z GW170817.
Path B precyzuje: `m_σ = √2 · m_s` (NIE `m_σ ~ Φ₀` bezpośrednio).
- Jeśli `m_s` różny od `Φ₀` (np. m_s ~ 0.5 meV << Φ₀ ~ ?), tension znika.
- Φ₀ to scale aksjonu/normalizacji, m_s to fizyczna masa fluktuacji wokół vacuum.
- Bond-renormalized m_s² > 0 z H_Γ J-coupling (paper §B), niekoniecznie = J·Φ₀.

To jest **częściowa resolution** (R2 z T3.2): `Φ₀ rozłączony od m_s`.
Ostatecznie m_s wyznacza testowalna LIGO 3G dispersion i ULDM bound — nie skala Φ₀.

### 3.5 EHT compatibility (M9.1'' standalone)

T-PB.5 (σ_ab = 0 dla static spherical) jest niezależnym potwierdzeniem
że photon ring prediction +14.56% (M9.1'' standalone) **pozostaje czystą
predykcją single-Φ M9.1''**, nie jest modyfikowane przez tensor sector.

→ Decision tree OP-M92 (Candidate D priority) niezmieniony.

---

## 4. Patch dla OP-7 T3 documentation

W osobnym pliku [[correction_to_OP7_T3.md]] (do utworzenia w step-final):
- T3.1 sub-results table: dodać "T-PB closure 2026-04-26 — Path B promoted to primary".
- T3.2 Φ₀/m_σ tension: zaznaczyć "częściowa resolution z Path B (m_s ≠ Φ₀)".
- T3.3 ghost analysis: cross-link do T-PB.3 (Gram-matrix positivity).

---

## 5. Co T-PB *nie* zamyka

**Otwarte (do następnych faz / closure_2026-04-26):**
- **bond-renormalization m_s skali** — wymaga lattice analizy v2 GL-substrate
  (skala J·Φ₀ vs |V''(s_eq)|). Nie blokuje Path B; precyzuje tylko numerikę.
- **higher-OPE rest R_ab** — T-PB.1 zostawia higher-derivative bilinear correlator
  jako sub-leading w `1/B²·k²`. Pełny rachunek wymaga RG flow w substrate
  (poza T-PB scope; należy do future work).
- **Bethe-Salpeter analysis** dla bound state interpretation `m_σ → 0` Hipoteza C.
  Path B z `M² = 2 m_s²` przenosi pytanie na Bethe-Salpeter: czy coupled
  dwóch ŝ-quanta w channel s tworzy bound state z masą ~ 0?
  Nie wymagane do T-PB closure (m_σ > 0 i tak daje GW-safe decoupling).

---

## 6. Werdykt T-PB

✅ **POSITIVE — 11/11 PASS**

σ_ab Path B (heredity z δŝ-EOM) jest:
1. **Primary derivation** σ_ab dynamiki — algebraicznie z `□δŝ + m_s² δŝ = J`.
2. **Single-substrate-faithful** — żaden nowy d.o.f., żaden nowy postulat.
3. **Ghost-free z konstrukcji** — Gram-matrix positivity.
4. **m_σ² = 2 m_s² derived** — composite-mass identity z box-algebry.
5. **Compatible z M9.1'' P3 PPN** — σ_ab = 0 static spherical (1e-18 numeric).

TGP po T-PB jest **strukturalnie czystszą** teorią niż przed: Path A
zostaje jako efektywny opis, ale ontologia jest *jawnie* single-substrate.

---

## 7. Pliki

- `setup.md` — audyt design + 5 kryteriów PASS
- `sigma_ab_pathB_audit.py` — sympy + numpy script (T-PB.1..T-PB.5)
- `sigma_ab_pathB_audit.txt` — raw output 11/11 PASS
- `results.md` — ten plik (synteza)
- `correction_to_OP7_T3.md` — patch note (do utworzenia w step-final)
