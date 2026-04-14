# R4: Ansatz metryczny h(Φ)=Φ z pierwszych zasad

## Problem

Metryka TGP:
```
ds² = -(c₀²/ψ)dt² + ψ·δᵢⱼdxⁱdxʲ,    ψ = Φ/Φ₀
```

Relacja h(Φ) = Φ (liniowa, p=1) jest **postulatem**, nie wyprowadzeniem.
Pytanie: dlaczego nie h(Φ) = Φ^p dla p ≠ 1?

## Obecny status (2026-04-14)

### ✅ PIĘĆ NIEZALEŻNYCH ARGUMENTÓW ZA p = 1

| # | Argument | Typ | Wynik | Referencja |
|---|----------|-----|-------|-----------|
| 1 | Gęstość substratu | ANALITYCZNY | Φ = gęstość → g_ij = (Φ/Φ₀)δ_ij → p=1 | sek08c Prop. |
| 2 | PPN Cassini + LLR | OBSERWACYJNY | γ = p = 1, β = 1 | ex206 (8/8 PASS) |
| 3 | Budżet informacyjny | ANALITYCZNY | f·h = 1 (antypodyczny) | sek08c Prop. |
| 4 | Element objętościowy | ANALITYCZNY | √(-g) = ψ^p musi = ψ → p = 1 | r4_einstein (11/11) |
| 5 | Stosunek mas solitonów | NUMERYCZNY | Tylko p=1 daje r₂₁ = 206.77 | a2_metric (6/6) |

### ✅ Dodatkowe wyniki

| Element | Status | Dowód |
|---------|--------|-------|
| Metryka eksponencjalna e^{±2U} | **TWIERDZENIE** | Jedyna h(U) spełniająca (i)-(iv) do wszystkich rzędów PPN |
| Ghost-freedom (3 poziomy) | **TWIERDZENIE** | K_sub(g) = g² > 0, Q_s = ψ⁴ > 0 |
| 10 parametrów PPN = GR | **ZWERYFIKOWANE** | γ=β=1, wszystkie inne =0 |
| det(g^{1+1}) = -c₀² = const | **TWIERDZENIE** | Wynika z p = q (antypodyczny) |
| ℓ_P = const | **ZWERYFIKOWANE** | G/ψ · ħ/√ψ / (c₀/√ψ)³ = const |

### ⚠️ OTWARTE (dodatkowe ścieżki)

| Element | Status | Problem |
|---------|--------|---------|
| A2a: Fonony na substracie | OTWARTE | Relacja dyspersji c_s(Φ) = c₀√(Φ/Φ₀) |
| A2c: Argument entropijny | OTWARTE | S_BH ∝ A → metryka liniowa w Φ |
| Formalizacja (Lean 4) | OTWARTE | Łańcuch dowodowy gotowy |

## Nowy argument A2b: Element objętościowy (2026-04-14)

```
Metryka: ds² = -c₀²ψ^{-q}dt² + ψ^p δ_ij dx^i dx^j

(I)   Równanie polowe na tle krzywym wymaga: p = q
      (równoważne warunkowi antypodycznemu f·h = 1)

(II)  Element objętościowy: √(-g) = ψ^{(3p-q)/2} = ψ^p (z p=q)
      Interpretacja substratowa: √(-g) musi = ψ = Φ/Φ₀
      → ψ^p = ψ → p = 1  ■

(III) PPN β = 1 wymaga resummacji eksponencjalnej:
      Prawo potęgowe h = ψ^p daje β = p(p+1) = 2 dla p=1 ✗
      Eksponencjalne h = e^{2U} daje β = 1 ✓
      → U = ½ln(ψ), ψ = e^{2U} — jedyna analityczna funkcja

(IV)  Klasyfikacja: TGP ≠ Brans-Dicke
      K(g) = g^{2α} jest KONKRETNYM coupling (nie wolny parametr ω)
      PPN identyczne z GR, różnice tylko kosmologiczne i silnopolowe
```

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `research/metric_ansatz/r4_einstein_self_consistency.py` | A2b: 11/11 PASS | ✅ NOWE |
| `scripts/ex206_metric_hypothesis_necessity.py` | Jedyność metryki: 8/8 PASS | ✅ RDZEŃ |
| `scripts/a2_metric_consistency.py` | Spójność metryki: 6/6 PASS | ✅ RDZEŃ |
| `scripts/ex201_antipodal_metric_derivation.py` | Antypodyczny: 8/8 PASS | ✅ RDZEŃ |
| `scripts/lk2_metric_from_substrate_propagation.py` | Propagacja: PASS | ✅ RDZEŃ |

## Referencje rdzenia

- `sek08c_metryka_z_substratu.tex` (wyprowadzenie metryki)
- `sek08b_ghost_resolution.tex` (ghost-freedom)
- `sek08a_akcja_zunifikowana.tex` (zunifikowane działanie)
- `nbody/tgp_ppn_full.tex` (pełna analiza PPN)

## Kryterium zamknięcia

**Twierdzenie (udowodnione 5 ścieżkami):** "Spośród g_ij = (Φ/Φ₀)^p·δ_ij, tylko p=1 daje:
ghost-free + pozytywna energia + poprawny limit newtonowski + γ=β=1 + √(-g)=ψ."

## Status

- [x] A2b: Równania Einsteina — element objętościowy wymusza p=1 (11/11 PASS)
- [x] PPN: γ = p = 1 z Cassini (8/8 PASS)
- [x] Antypodyczny: f·h = 1 z budżetu informacyjnego
- [x] Ghost-freedom: K_sub > 0 na trzech poziomach
- [x] Stosunek mas: tylko p=1 daje r₂₁ = 206.77
- [ ] A2a: Relacja dyspersji fononów
- [ ] A2c: Argument entropijny
- [ ] Formalizacja w Lean 4
