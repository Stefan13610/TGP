# Plan numeryczny: domknięcie CG-3 i CG-4

## Cel

Operacyjny plan numeryczny wspierający **lematy A1–A5** (dodatekQ2)
i zamknięcie otwartych kroków CG-3 (zbieżność Φ_B → Φ w H¹)
oraz CG-4 (identyfikacja K_hom = K_TGP).

Każdy etap (N1–N5) ma jasno określone:
- **wejście**: dane/skrypty potrzebne do uruchomienia,
- **wyjście**: mierzalna wielkość + kryterium PASS/FAIL,
- **lemat wspierany**: który lemat A1–A5 korzysta z tego wyniku.

## Stan obecny

| Skrypt | CG-krok | Status | Wynik |
|--------|---------|--------|-------|
| `tgp_erg_lpa_prime.py` | CG-2 | ZAMKNIĘTY (8/8) | K_IR/K_UV = 1.000, ν = 0.749 |
| `tgp_cg5_phi0_self_consistency.py` | CG-5 | ZAMKNIĘTY (8/8) | a_Γ·Φ₀ = 1.000000 |
| `tgp_erg_eta_lpa_prime.py` | OP-1 | ZAMKNIĘTY (10/10) | η* = 0.04419 |
| `p121_coarsegraining_convergence.py` | archiwum | ARCHIWALNY (v42) | β/γ ≈ 1 |

Brakuje: dedykowanych skryptów MC dla etapów N1–N5.

---

## Etap N1: Monte Carlo + block averaging

**Cel**: wygenerować konfiguracje substratu i zmierzyć pole blokowe Φ_B
dla wielu rozmiarów bloku b.

**Wspiera**: Lemat A1 (dane wejściowe), Lemat A2 (weryfikacja lokalności)

### Specyfikacja

- **Model**: Z₂-symetryczny, H = -J Σ (φ_i φ_j)², φ_i ∈ ℝ, d = 3
- **Siatka**: L³ z L = 32, 48, 64, 96
- **Temperatura**: T/T_c ∈ {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
- **Bloki**: b ∈ {2, 4, 8, 16} (tj. L_B = b · a_sub)
- **Algorytm**: Wolff cluster + overrelaxation (szybka dekorelacja)
- **Statystyka**: 10⁵ konfiguracji po termalizacji, co 10-ty sweep

### Wyjście

- Konfiguracje Φ_B(x) dla każdego (L, T, b)
- Histogram P(Φ_B) dla każdego zestawu
- Korelator dwupunktowy ⟨Φ_B(x) Φ_B(y)⟩_c → długość korelacji ξ(T, b)

### Kryterium PASS

- ξ(T) rośnie jak |1 - T/T_c|^{-ν} z ν ∈ [0.62, 0.64] (3D Ising)
- Separacja skal: a_sub ≪ L_B ≪ ξ spełniona dla T/T_c ≤ 0.85 i b ≤ 8

### Skrypt

`scripts/substrate/substrate_mc_cg3.py` (DO NAPISANIA)

Bazuje na archiwum: `_archiwum/scripts_exploratory/substrate_mc_z2.py`

---

## Etap N2: Ekstrakcja V_eff(Φ)

**Cel**: z rozkładu histogramowego Φ_B wyznaczyć potencjał efektywny U(Φ).

**Wspiera**: Lemat A4 (identyfikacja współczynników)

### Metoda

1. Z histogramu P(Φ_B) wyznaczyć: U_eff(Φ) = -ln P(Φ_B) + const
2. Dopasować do rozwinięcia:
   ```
   U_eff(Φ) = U₀ + ½ m²_sp (Φ - Φ₀)² + u₃ (Φ - Φ₀)³ + u₄ (Φ - Φ₀)⁴
   ```
3. Zidentyfikować:
   - Φ₀ = argmin U_eff → próżnia
   - m²_sp = U''(Φ₀) → masa skalarna
   - β_eff = m²_sp / (2Φ₀)
   - γ_eff z u₃ → test β = γ

### Wyjście

- Tabela (T, b, Φ₀, m²_sp, β_eff, γ_eff) dla każdego zestawu
- Wykres U_eff(Φ) vs. Φ z dopasowaniem

### Kryterium PASS

- Φ₀ > 0 dla T < T_c (istnienie spontanicznego łamania)
- m²_sp > 0 (stabilność próżni)
- |β_eff/γ_eff - 1| < 0.1 (zgodność z TGP: β = γ)
- Φ₀ zbieżne przy b → ∞ (niezależność od rozdzielczości)

### Skrypt

`scripts/substrate/substrate_veff_extraction.py` (DO NAPISANIA)

---

## Etap N3: Ekstrakcja sektora kinetycznego

**Cel**: zmierzyć efektywną sztywność K₁(Φ) i dolny szacunek c*.

**Wspiera**: Lemat A1 (dolny szacunek gradientowy), Lemat A3 (K₁ ∝ 1/Φ)

### Metoda

1. Obliczyć gradient blokowy: ∇_lat Φ_B(x) = [Φ_B(x+ê) - Φ_B(x)] / L_B
2. Zmierzyć korelator gradientowy:
   ```
   C_kin(Φ) = ⟨|∇_lat Φ_B|² | Φ_B ≈ Φ⟩
   ```
   jako funkcję Φ (binowanie warunkowe)
3. Wyznaczyć K₁(Φ) = C_kin(Φ) / (normalizacja geometryczna)
4. Sprawdzić:
   - K₁(Φ) > 0 dla Φ > 0 (dodatniość → c* > 0)
   - K₁(Φ) ∝ 1/Φ blisko Φ₀ (zgodność z α = 2)

### Wyjście

- Profil K₁(Φ) vs. Φ dla każdego (T, b)
- Wartość c* = min_{Φ > 0} K₁(Φ)
- Test skalowania: K₁(Φ) · Φ ≈ const (α = 2 check)

### Kryterium PASS

- c* > 0 (KRYTYCZNE — bez tego Lemat A1 upada)
- K₁(Φ) · Φ = const ± 10% w zakresie [0.5 Φ₀, 2 Φ₀]
- c* stabilne przy zmianie b (niezależność od rozdzielczości)

### Skrypt

`scripts/substrate/substrate_kinetic_sector.py` (DO NAPISANIA)

---

## Etap N4: Finite-size scaling

**Cel**: ekstrapolacja L → ∞ i kontrola artefaktów skończonego rozmiaru.

**Wspiera**: Lemat A1 (jednostajność ograniczeń), Twierdzenie A5 (kontrola błędu)

### Metoda

1. Dla każdego T i b powtórzyć pomiary z N1–N3 na siatkach L = 32, 48, 64, 96
2. Ekstrapolować do L → ∞:
   - Φ₀(L) → Φ₀(∞) via fit 1/L^d
   - m²_sp(L) → m²_sp(∞)
   - c*(L) → c*(∞)
3. Porównać ξ z L: reżim ξ/L < 0.3 (aby efekty skończonego rozmiaru < 5%)
4. Sprawdzić optymalny rozmiar bloku: b_opt taki, że a_sub ≪ b·a_sub ≪ ξ

### Wyjście

- Tabela ekstrapolowanych wartości (Φ₀, m²_sp, c*, β_eff, γ_eff) w L → ∞
- Wykres Φ₀(L) z fitowaniem FSS

### Kryterium PASS

- Ekstrapolacje stabilne (χ² / ndof < 2)
- Φ₀(L=96) ≈ Φ₀(∞) do 5%
- Reżim b_opt istnieje dla T/T_c ≤ 0.85

### Skrypt

`scripts/substrate/substrate_fss_extrapolation.py` (DO NAPISANIA)

---

## Etap N5: Test zgodności operatora

**Cel**: najsilniejszy test — podstawić zmierzone Φ_B do PDE TGP
i zmierzyć rezydua.

**Wspiera**: Twierdzenie A5 (identyfikacja K_hom = K_TGP), CG-4

### Metoda

1. Dla konfiguracji Φ_B z N1, obliczyć residuum:
   ```
   R(x) = ∇² Φ_B + (∇Φ_B)²/Φ_B + β (Φ_B²/Φ₀ - Φ_B³/Φ₀²)
   ```
   używając β, Φ₀ z N2.
2. Zmierzyć ||R||² / ||∇²Φ_B||² (względne residuum)
3. Sprawdzić, czy residuum maleje z:
   - rosnącym b (lepsza separacja skal)
   - malejącym T/T_c (dłuższa ξ, lepsze continuum)
4. Fitować: ||R|| ~ (L_B/ξ)^p₁ + (a_sub/L_B)^p₂
   i porównać z przewidywaniem A5: p₁ ≥ 1, p₂ ≥ 1/2

### Wyjście

- Tabela (T, b, ||R||_rel) — mapa residuów
- Wykres ||R|| vs. L_B/ξ — skalowanie
- Fitted exponents p₁, p₂

### Kryterium PASS

- ||R||_rel < 0.1 dla co najmniej jednego (T, b) (10% residuum)
- ||R|| maleje monotonicznie z rosnącym b (zbieżność)
- p₁ ≥ 0.8 i p₂ ≥ 0.3 (zgodność z A5 do 50%)

### Skrypt

`scripts/substrate/substrate_operator_test.py` (DO NAPISANIA)

---

## Podsumowanie: mapa zależności

```
N1 (MC + blocks)  ──→  N2 (V_eff)  ──→  N4 (FSS)
       │                    │                 │
       ↓                    ↓                 ↓
  N3 (kinetic)         A4 (coeff)        A5 (error)
       │
       ↓
  A1 (c* > 0)
```

## Kolejność realizacji

| Priorytet | Etap | Szacowany czas | Zależności |
|-----------|------|----------------|------------|
| 1 | N1 (MC + blocks) | 2–3 dni | brak (bazowy) |
| 2 | N3 (kinetic sector) | 1–2 dni | N1 |
| 3 | N2 (V_eff extraction) | 1–2 dni | N1 |
| 4 | N4 (FSS) | 2–3 dni | N1, N2, N3 |
| 5 | N5 (operator test) | 1–2 dni | N1, N2 |

**Czas łączny**: ~8–12 dni roboczych

## Kryterium sukcesu (zamknięcie CG-3 + CG-4 numerycznie)

CG-3 i CG-4 uznajemy za **zamknięte numerycznie** gdy:

1. ✅ c* > 0 (A1 supported) — z N3
2. ✅ β_eff/γ_eff ≈ 1 (β = γ recovered) — z N2
3. ✅ K₁(Φ) · Φ ≈ const (α = 2 recovered) — z N3
4. ✅ ||R||_rel → 0 przy b → ∞ (operator convergence) — z N5
5. ✅ Φ₀ → finite limit przy L → ∞ (continuum exists) — z N4

Spełnienie tych 5 warunków + Lemat A3 (algebraiczny) daje
**pełne numeryczne wsparcie** dla Twierdzenia A5.
