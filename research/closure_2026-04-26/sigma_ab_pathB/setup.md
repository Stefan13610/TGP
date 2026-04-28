# σ_ab Path B audit — single-substrate-faithful heredity (T-PB)

**Data:** 2026-04-26
**Status:** OPEN
**Parent:** [[research/op7/OP7_T3_results.md]] (T3.1: Path A == Path B structurally)
**Strategic ref:** [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.1
**Foundations binding:** [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂ — immovable)

---

## 1. Cel audytu

OP-7 T3.1 pokazał, że **Path A** (efektywny Lagrangian L_σ z postulowanym
masowym członem m_σ²) i **Path B** (heredity z linearizacji ŝ-EOM)
prowadzą strukturalnie do tej samej EOM z identyfikacją `m_σ² = 2 m_s²`.

ALE zamknięcie T3 oparte zostało roboczo na Path A (wygodniejszy
do PPN/quadrupole rachunku). Path B traktowany jako **independent
cross-check**, nie jako primary derivation.

**Audyt T-PB** ma zamknąć tę asymetrię i upewnić się, że:

1. **Path B jest primary** — σ_ab dynamika jest *emergentna* z ŝ-EOM,
   nie wymaga niezależnego Lagrangianu L_σ ani postulowania m_σ².
2. **Identyfikacja `m_σ² = 2 m_s²` jest derived**, nie matched —
   wynika z OPE/composite-operator analizy, nie z analogii mezonowej.
3. **Brak nowych d.o.f.** — σ_ab pozostaje correlator-like,
   nie staje się quasi-fundamental field z własną kinetyką.
4. **Source term** TT-projekcja produkowana przez bilinear coupling
   ŝ-pola do materii (jeśli istnieje) lub przez mode-mode mixing
   w samym substracie.
5. **Ghost-free z konstrukcji** — pozytywno-określona Hamiltonowska
   gęstość ŝ implikuje pozytywno-określoną energię σ_ab.

**Kryterium PASS:** wszystkie 5 punktów weryfikowalne sympy +
analitycznie z `□ŝ + V'(ŝ) = 0` przy V(ŝ) = (β/3)ŝ³ − (γ/4)ŝ⁴ z
M9.1''/OP-7 closure.

**Kryterium FAIL:** którykolwiek punkt wymaga dodatkowego postulatu
poza ŝ-działaniem → TGP straciłby property "single-substrate-faithful"
w sektorze tensora.

---

## 2. Ramowy szkic Path B

### 2.1 Equation of motion ŝ (M9.1''+OP-7 v2)

Z S_TGP[ŝ] z V(ŝ) = (β/3)ŝ³ − (γ/4)ŝ⁴ (paper §B, β=γ vacuum),
EOM dla ŝ wokół próżni s_eq² = Φ₀ = β/γ:

```
□ŝ + V'(ŝ) = J_ŝ[matter]    (1)
```

gdzie V'(ŝ) = β ŝ² − γ ŝ³.

### 2.2 Linearizacja wokół s_eq = √Φ₀

Postaw ŝ = s_eq + δŝ. Z V'(s_eq) = 0 (vacuum) i V''(s_eq) = 2β·s_eq − 3γ·s_eq²:

```
V''(s_eq) = s_eq · (2β − 3γ s_eq) = √Φ₀ · (2β − 3γ √Φ₀)        (2)
```

Z β = γ √Φ₀ (vacuum-condition w nowych jednostkach: M9.1'' P2):
`V''(s_eq) = √Φ₀ · (2γ√Φ₀ − 3γ √Φ₀) = −γ Φ₀`.

(Znak ujemny przy V'' jest cechą bond-renormalized vacuum w M9.1'';
fizyczna masa pojawia się z full H_Γ, paper §B oryginalnie.)

Bond-renormalized fizyczne `m_s² > 0` (z J-coupling, paper §B).
Ozn. `m_s² ≡ J·Φ₀ + |V''(s_eq)|` > 0.

Linearizacja:
```
□(δŝ) + m_s² δŝ = J_ŝ_lin[matter]    (3)
```

### 2.3 σ_ab jako bilinear coarse-grain operator (Path B definicja)

Z OP-7 T2:
```
σ_ab(x) ≡ ⟨(∂_a δŝ)(∂_b δŝ)⟩_B(x) − (1/3) δ_ab Tr ⟨(∂ δŝ)²⟩    (4)
```

gdzie `⟨...⟩_B` to local coarse-grain average na skali B
(np. λ_subst — siatkowa skala substratu).

### 2.4 Heredity equation (kluczowy fragment Path B)

Stosujemy □ do (4). Używamy tożsamości:
```
□[A·B] = (□A)·B + 2 ∂_μ A · ∂^μ B + A · (□B)        (5)
```

z A = ∂_a δŝ, B = ∂_b δŝ:
```
□[(∂_a δŝ)(∂_b δŝ)] =
   (□∂_a δŝ)(∂_b δŝ) + 2 (∂_μ ∂_a δŝ)(∂^μ ∂_b δŝ) + (∂_a δŝ)(□∂_b δŝ)
```

Z (3): `□∂_a δŝ = ∂_a (□ δŝ) = ∂_a (J_ŝ_lin − m_s² δŝ) = ∂_a J_ŝ_lin − m_s² ∂_a δŝ`.

Pierwsze i trzecie człony:
```
(□∂_a δŝ)(∂_b δŝ) + (∂_a δŝ)(□∂_b δŝ) =
   −2 m_s² · (∂_a δŝ)(∂_b δŝ) + (∂_a J_ŝ_lin)(∂_b δŝ) + (∂_a δŝ)(∂_b J_ŝ_lin)
```

Drugi człon `2 (∂_μ ∂_a δŝ)(∂^μ ∂_b δŝ)` to higher-derivative
correlator. W coarse-grain limit (B-skala >> mode wavelength k⁻¹)
ten człon redukuje się do **renormalizacji kinetycznej części**:
po koalescencji punktowej operatorów składowych pozostaje jako
**niezerowa contribution** do gradient norm K_ab; nie wnosi nowych
d.o.f. (wszystko dziedziczy z ŝ).

Po wzięciu coarse-grain ⟨⟩_B i odjęciu śladu:
```
□σ_ab + 2 m_s² σ_ab = S_ab^TT[matter] + R_ab[higher-OPE]    (6)
```

gdzie:
- `M² ≡ 2 m_s²` — **derived composite mass** (jak w mezonach: składową
  jest dwóch masywnych ŝ-quanta przy threshold).
- `S_ab^TT` — TT-projekcja źródła z `(∂_a J_lin)(∂_b δŝ) + ...`,
  symetryzowana i bezśladowa.
- `R_ab` — higher-OPE rest (sub-leading w ekspansji 1/B²·k²).

To jest **derived heredity equation**, nie postulat.

### 2.5 Ghost-free z heredity

ŝ ma:
- Kinetyczna gęstość: `(1/2)(∂_μ δŝ)² ≥ 0` (canonical kinetic).
- Hamilton positive-definite na próżni stable (m_s² > 0 z bond
  renormalization).

σ_ab to bilinear w `∂δŝ` z dodatnimi kowariantnymi metrykami
podstawowego pola → **ghost-free z konstrukcji**.

Ostrogradski: heredity equation pierwszego rzędu w czasie po
linearizacji (∂²σ_ab/∂t² wynika z ∂²δŝ/∂t² w (3)) → brak
higher-derivative ghosts.

### 2.6 Single-substrate axiom zachowany

W całym Path B:
- Tylko jedno pole `δŝ` w działaniu.
- σ_ab to **operator** z ŝ-pola, nie niezależne pole.
- Brak nowego momentu kanonicznego (wszystko z `π_ŝ = ∂L/∂(∂_t δŝ)`).
- Counting d.o.f.: 1 (δŝ scalar) → 1 propagating breathing + emergent
  spin-2 fluktuacje obciążone naturalną korelacją z δŝ-modes.

---

## 3. Plan testów (T-PB.1 ... T-PB.5)

| ID | Cel | Metoda | PASS |
|----|-----|--------|------|
| **T-PB.1** | Sympy: derive (6) z (1) krok-po-kroku | sympy: V(ŝ), linearize, □[bilinear], coarse-grain | EOM (6) literalnie z (1) bez postulatu m_σ² |
| **T-PB.2** | Identyfikacja `M² = 2 m_s²` | OPE coalescence: composite mass z dwóch propagatorów | M² = 2 m_s² **derived** (numeric check) |
| **T-PB.3** | Ghost-free z konstrukcji | sympy: Hamiltonian σ_ab z H_ŝ | H_σ ≥ 0 implikowane przez H_ŝ ≥ 0 |
| **T-PB.4** | Brak nowych d.o.f. (counting) | symbolic phase-space accounting | 1 ŝ-d.o.f. → 1 scalar + emergent σ_ab correlator |
| **T-PB.5** | Reduction do M9.1'' przy ⟨⟩_B-limit | static spherical: σ_ab → 0 | M9.1'' P3 PPN unchanged |

**Sukces:** 5/5 PASS → σ_ab dynamika jest **strukturalnie konsekwentna
z single-Φ Z₂** bez konieczności dodatkowego Lagrangianu L_σ.

---

## 4. Dlaczego to ma znaczenie

T3 zamknął się jako "STRUCTURAL POSITIVE z OPEN TENSION (Φ₀/m_σ)".
Ale *kierunek* domknięcia był taki, że Path A traktowany był jako
primary, a Path B jako consistency check.

Audyt T-PB **odwraca tę kolejność:**

- Path B (heredity) staje się primary derivation, bo jest
  **formalnie redundant w aksjomatach** (nie wprowadza nic poza ŝ-EOM).
- Path A staje się **efektywnym opisem dla obliczeń** (PPN, quadrupole)
  którego używamy gdy heredity equation (6) jest już zamknięta.

Konsekwencje strategiczne:

1. **Spirit-fidelity** TGP_FOUNDATIONS §1 jest *jawnie* potwierdzona —
   żadnego "second field" nie ma w tle nawet operacyjnie.
2. **Ghost-free** staje się **structural property** ŝ-pola, nie
   wynik dobrego dobrania znaków w L_σ.
3. **m_σ² = 2 m_s²** staje się **predykcją**, nie założeniem —
   jeśli dispersion testowi w LIGO 3G da wartość niezgodną z 2 m_s²,
   TGP jest sfalsyfikowany (vs Path A gdzie m_σ było free parameter).
4. **Decoupling regime** (OP-7 T6: 2m_s ~ meV) staje się jednoznacznie
   określony: M_σ = 2 m_s = 2 √Φ₀ (z bond renormalization), nie
   niezależny parameter.

---

## 5. Pliki

- `setup.md` (this file) — audyt design
- `sigma_ab_pathB_audit.py` — sympy script T-PB.1..T-PB.5
- `sigma_ab_pathB_audit.txt` — raw output
- `results.md` — synthesis post-execution
- `correction_to_OP7_T3.md` — patch note dla OP-7 T3.1 (Path B promotion)

---

## 6. Cross-references

- [[research/op7/OP7_T3_results.md]] (T3 closure z Φ₀/m_σ tension)
- [[research/op7/OP7_T6_results.md]] (decoupling regime, Λ(ψ)=1)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.1 (brainstorm)
- [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂ axiom)
- [[../../../tgp-core-paper/paper/tgp_core.tex]] §B (H_Γ, V(ŝ) potential)
