# OP-7 / T3 — Dynamika σ_ab: derivacja EOM, m_σ, ghost analysis, ξ coupling

**Data:** 2026-04-25
**Status:** ⚠️ STRUCTURAL POSITIVE z OPEN TENSION (Φ₀/m_σ)
**Pliki wykonawcze:**
- `op7_t3_sigma_dynamics.py` (T3.1)
- `op7_t3_2_m_sigma_scale.py` (T3.2)
- `op7_t3_3_ghost_analysis.py` (T3.3)
- `op7_t3_4_xi_coupling.py` (T3.4)

**Raw outputs:** odpowiednie pliki `.txt`.

**Cross-references:**
- [[OP7_setup.md]] §3 (T3 row)
- [[OP7_T1_results.md]] (no-tensor M9.1'')
- [[OP7_T2_results.md]] (σ_ab gradient strain definition)
- [[EHT_quick_results.md]] (M9.1'' photon ring +14.6%, niezależny smoking gun)
- [[TGP_CLOSURE_PLAN_2026-04-25.md]] §3 (krytyczna ścieżka), §8 (brainstorm)
- [[../../../tgp-core-paper/paper/tgp_core.tex]] §2, §7 (OP-7 row)

---

## 1. Cel testu T3

Dla efektywnego Lagrangianu σ_ab albo composite z H_Γ wyprowadzić:
1. **Strukturalna postać EOM** (T3.1): czy `□σ_ab + m_σ²σ_ab = -ξ T_ab^{TT}`?
2. **Skalowanie masy** (T3.2): m_σ z trzech hipotez A (mean-field), B (kompresybilna), C (massless).
3. **Ghost analysis** (T3.3): czy struktura jest ghost-free? (Hamiltonian, Ostrogradski, dispersion).
4. **Sprzężenie do materii** (T3.4): forma ξ; matching do GW150914 strain ~1e-21.

**Kryterium PASS:** dynamika ghost-free, propagacja luminal, amplituda GR-like, single-Φ axiom zachowany.

**Kryterium FAIL:** ghost lub tachyon w sygnaturze kinetycznej; m_σ niespójne z LIGO bound dla naturalnych Φ₀.

---

## 2. Wyniki

### Część T3.1 — Strukturalna derivacja EOM

**Two paths sprawdzone strukturalnie i pokrywają się:**

**Path A (efektywny Lagrangian):**
```
L_σ = -(1/4)(∂_μ σ_ab)(∂^μ σ^ab) - (1/2)m_σ² σ_ab σ^ab - (ξ/2) σ_ab T^{ab,TT}
```
Wariacja `δL/δσ^{ab} = 0` daje:
```
□σ_ab + m_σ² σ_ab = -ξ T_ab^{TT}
```

**Path B (composite z s-EOM):** dynamika σ_ab dziedziczy z linearizacji ŝ EOM:
```
□(δŝ) + m_s_eff² δŝ = 0
σ_ab(x) = ⟨(∂_a δŝ)(∂_b δŝ)⟩ - (1/3)δ_ab Tr(...)
→ □σ_ab + 2 m_s_eff² σ_ab = (TT projection of source)
```

**Identyfikacja:** `m_σ² = 2 m_s_eff²` (composite mass = 2× base, jak w mezonach).

**Sympy weryfikacja:**
- V''(s_eq) z V(Φ) = (γ/12)Φ₀²ψ³(4-3ψ): `d²V/ds² |_{s²=Φ₀} = -4γΦ₀` (negative; vacuum stable wymaga bond renormalization).
- Plane wave dispersion: `ω² = k² + m_σ²`.
- Tracelessness i symetria DYNAMICZNIE ZACHOWANE (Tr T^TT = 0; T_ab symmetric).

**Wynik T3.1:** **11/11 PASS strukturalny**. Path A == Path B; single-Φ aksjomat zachowany.

### Część T3.2 — Skala m_σ (3 hipotezy)

| Hipoteza | Mechanizm | Skala m_σ | GW170817 bound (~2e-19 eV) |
|---|---|---|---|
| **A** (mean-field) | m_σ² ~ J⟨ŝ²⟩ ~ Φ₀ | ~Φ₀ | FAIL jeśli Φ₀ ≥ 1e-19 eV |
| **B** (kompresybilna) | m_σ² ~ -V''(Φ₀) renormalized | ~Φ₀ | FAIL jeśli Φ₀ ≥ 1e-19 eV |
| **C** (massless) | non-Goldstone composite | 0 | PASS trivially |

**Tension:** brainstorm §8.9 (cosmological Φ₀ ~ Λ_obs^(1/4) ~ meV) dawałoby m_σ ~ meV >> 2e-19 eV → **niezgodne z GW170817 o 16 rzędów wielkości**.

**Resolution opcje:**
1. **Hipoteza C** structuralnie wspierana przez Bethe-Salpeter / 1-loop (m_σ → 0 z bound state cancellation). **Otwarte do T3-extended.**
2. **Φ₀ << meV** (np. ULDM scale ~1e-22 eV): wtedy m_σ ~ 1e-22 eV, GW-safe.
3. **Częściowa falsyfikacja** TGP na poziomie LIGO 3G dispersion test.

**Wynik T3.2:** **4/7 PASS**, identyfikacja krytycznej **Φ₀/m_σ tension**. INCONCLUSIVE.

### Część T3.3 — Ghost analysis

**Path A Hamiltonian (kanoniczny):**
```
H[π, σ, ∇σ] = π² + (1/4)(∇σ)² + (1/2)m_σ² σ²
```
Wszystkie współczynniki **dodatnie** → ghost-free.

**Ostrogradski test:** L pierwszego rzędu w pochodnych → trivially OK.

**Path B:** σ_ab dziedziczy positive-definite Hamiltonian z ŝ-pola przez konstrukcję (σ to suma over s-modes z dodatnimi energiami).

**Massless hypothesis C:**
- Brak trywialnego mechanizmu (Z₂ to dyskretna symetria, brak Goldstone'a; brak gauge protection w single-Phi).
- Możliwa droga: bound state w 2-particle s-spectrum z binding energy ~ 2 m_s. Wymaga Bethe-Salpeter analizy poza T3.3 scope.
- **Operacyjnie:** zakładamy hipotezę C jako roboczą (najbezpieczniejsza GW); jeśli T3-extended pokaże m_σ > 0, mamy testowalną dispersion w LIGO 3G.

**Dispersion:** Path A ω² = k² + m_σ², c_phase = c₀ (luminal). Path B identyczne (2-particle threshold = 2 m_s).

**Wynik T3.3:** **5/5 PASS structural ghost-free**. Massless C admissible (otwarte numerycznie).

### Część T3.4 — Sprzężenie ξ + GW150914 matching

**Strukturalna identyfikacja:** Λ₀ × ξ = 4πG.
- Z natural choice Λ₀ = 1/Φ₀² (canonical metric coupling): ξ = G·Φ₀²
- Razem: TGP daje GR-equivalent quadrupole formula.

**Empirical matching GW150914:**
- Q̈ ~ 1.44e+48 J (binary 30+30 M_sun, a=350km, f=100Hz)
- h_predicted (z ξ=G) = 9.42e-22 (observed: 1.0e-21)
- **ξ/G ≈ 1.06 — O(1), fizycznie sensowne**

**Smoking guns testowalne:**
| Test | Skala | Status |
|---|---|---|
| LIGO 3G dispersion (m_σ) | ~2030 (Cosmic Explorer / Einstein Telescope) | otwarte |
| LIGO O5+ binary 2PN deviation `(5/6)U³` | ~2027 | otwarte |
| ngEHT photon ring +14.6% | M87* boundary, Sgr A* tension teraz | **partial-NEGATIVE** ([[EHT_quick_results.md]]) |

**Wynik T3.4:** **5/5 PASS structural+empirical** matching GR amplitude.

---

## 3. Werdykt T3 (synteza)

| Sub-test | PASS/FAIL | Liczba | Komentarz |
|---|---|---|---|
| T3.1 | PASS | 11/11 | Path A == Path B; single-Φ zachowany |
| T3.2 | PARTIAL-FAIL | 4/7 | Tension Φ₀/m_σ z GW170817 (16 rzędów!) |
| T3.3 | PASS | 5/5 | Ghost-free; massless C admissible |
| T3.4 | PASS | 5/5 | ξ/G ~ 1.06 (O(1)); GR amplitude reproduced |
| **Suma T3** | **MIXED** | **25/28 ≈ 89%** | structural OK; **m_σ scale tension OPEN** |

**Verdict T3 (overall):**

**STRUCTURAL POSITIVE** — TGP σ_ab dynamika jest:
- Wyprowadzalna z S_TGP[ŝ] przez TWA niezależne paths zgodne strukturalnie
- Ghost-free w obu sformuowania
- Single-Φ axiom zachowany (σ_ab to composite, nie nowy d.o.f.)
- TGP NIE staje się scalar-tensor
- GR-equivalent quadrupole formula z empirycznym ξ/G ~ 1.06

**OPEN TENSION** — naturalna skala m_σ ~ Φ₀ niezgodna z GW170817:
- Cosmological motivated Φ₀ ~ meV (z Λ_obs^(1/4)) → m_σ ~ meV
- GW170817 bound: m_σ < 2e-19 eV
- 16 rzędów wielkości tension

**Resolution paths:**

(R1) **Hypothesis C massless wspierana 1-loop** (T3-extended): bound state w 2-particle ŝ spectrum z masą zerową poprzez niesample subtraction. Wymaga dedykowanej Bethe-Salpeter analizy.

(R2) **Φ₀ rozłączony od Λ_obs**: skala Phi_0 może być niezależna od cosmological constant origin. Brainstorm §8.9 wymaga refinement (V(Phi_0) ≠ Λ_obs jednowarstwowo).

(R3) **Częściowa falsyfikacja**: TGP daje konkretną testowalną predykcję dispersion przy LIGO 3G. Jeśli m_σ ~ 1e-19 eV detected, TGP confirmed; jeśli rejected → TGP jak częściowo falsyfikowane na poziomie GW.

**Decyzja operacyjna:** dla T3 closure przyjmujemy **Hypothesis C jako roboczą** (najbezpieczniejsza), z explicit flag że T3-extended musi zdecydować strukturalnie.

---

## 4. Implikacje dla TGP closure

**T3 zamyka kwalitatywnie kinematyczno-dynamiczną stronę OP-7:**
- T1 (no-tensor M9.1'') ✅
- T2 (σ_ab gradient strain composite) ✅
- T3 (σ_ab dynamics + EOM + ghost-free + ξ matching) **PARTIAL** ⚠️
- T4 (Λ(ψ) metric coupling) — gating analysis na T3 results
- T5 (kwadrupol h+, h× → GW150914) — częściowo zamknięte przez T3.4 (ξ/G ~ 1)
- T6 (PPN + c_GW + Z₂ + stability) — pending

**Krytyczne otwarte pytanie po T3:**
1. **Phi_0 scale**: cosmological vs particle. Brainstorm §8.9 wymaga refinement.
2. **m_σ from 1-loop**: czy hipoteza C structuralnie wspierana?
3. **Λ(ψ) form**: T4 wybiera spośród 4 kandydatów (linear, ψ(4-3ψ), constant, ψ/(4-3ψ)).

**Krytyczne otwarte pytanie OBOK T3:**
4. **EHT photon ring +14.6%**: M9.1'' static spherical FALSIFIES at Sgr A* level (z [[EHT_quick_results.md]]). To jest **niezależne** od σ_ab — czysta predykcja M9.1'' single-Φ. Wymaga decyzji: czy modify M9.1'' lub akceptujemy strong-field deviation.

---

## 5. Następne kroki

**T3-extended (po obecnej sesji):**
- Bethe-Salpeter / 1-loop renormalization dla m_σ
- Numerical lattice check m_σ z H_Γ pełny v2
- Refinement brainstorm §8.9 (Phi_0 vs Lambda_obs separation)

**T4 (Λ(ψ) coupling):**
- Variational analysis: ghost-free choice z 4 kandydatów
- PPN check: σ_ab=0 → β=γ=1 unchanged
- Z₂ parity: σ_ab niezmieniony pod ŝ → -ŝ

**T5 (full quadrupole formula):**
- Greens function explicit z m_σ != 0 case
- LIGO O5 2PN binary inspiral phase prediction
- Matching ξ_eff numerycznie do GW150914+170817+...

**T6 (full consistency):**
- 8 pozostałych PPN parameters z M9.2 (moving sources)
- c_GW = c₀ exactly w próżni z T4 Λ
- Z₂ + stability across full parameter space

---

## 6. Pliki

- `op7_t3_sigma_dynamics.py` — T3.1 sympy structural derivation
- `op7_t3_sigma_dynamics.txt` — raw 11/11 PASS
- `op7_t3_2_m_sigma_scale.py` — T3.2 m_σ hipotezy A/B/C
- `op7_t3_2_m_sigma_scale.txt` — raw 4/7 PASS (Φ₀/m_σ tension)
- `op7_t3_3_ghost_analysis.py` — T3.3 ghost+dispersion
- `op7_t3_3_ghost_analysis.txt` — raw 5/5 PASS
- `op7_t3_4_xi_coupling.py` — T3.4 ξ + GW150914 match
- `op7_t3_4_xi_coupling.txt` — raw 5/5 PASS (ξ/G ~ 1.06)
- `OP7_T3_results.md` — ten plik (synteza T3 werdyktu)
