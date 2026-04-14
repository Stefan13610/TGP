# R4: Ansatz metryczny h(Φ)=Φ z pierwszych zasad

## Problem

Metryka TGP:
```
ds² = -(c₀²/ψ)dt² + ψ·δᵢⱼdxⁱdxʲ,    ψ = Φ/Φ₀
```

Relacja h(Φ) = Φ (liniowa, p=1) jest **postulatem**, nie wyprowadzeniem.
Pytanie: dlaczego nie h(Φ) = Φ^p dla p ≠ 1?

## Obecne uzasadnienie

| Argument | Typ | Wynik |
|----------|-----|-------|
| Fiber bundle condition | NUMERYCZNY | Wybiera p=1 |
| Mass ratio r₂₁ ≈ 206.77 | NUMERYCZNY | Wymaga p=1 |
| Shapiro delay | ZDEGENEROWANY | Dowolne p daje γ_PPN = 1 |
| PPN γ=1 | ZDEGENEROWANY | Nie rozróżnia p |
| a2_metric_consistency.py | 6/6 PASS | Numeryczne, nie analityczne |

## Plan ataku — trzy niezależne ścieżki

### A2a: Fonony na substracie (2 tygodnie)
- Obliczyć relację dyspersji ω(k) na substracie Γ z H = -J·Σ(φᵢ·φⱼ)²
- Prędkość dźwięku: c_s² = dω/dk|_{k→0}
- Hipoteza: c_s ∝ √Φ wymusza liniowy coupling h=Φ
- **Cel:** c_sound(Φ) = c₀·√(Φ/Φ₀) → ds² = Φ/Φ₀·δᵢⱼ

### A2b: Test równań Einsteina (2–3 tygodnie) ← NAJBARDZIEJ OBIECUJĄCY
- Wstawić g_ij = Φ^p·δ_ij do pełnych równań pola:
  ```
  G_μν = 8πG·T_μν[Φ]
  ```
- Sprawdzić jednocześnie:
  1. Ghost-free: brak duchów w sektorze skalarnym
  2. Pozytywna energia: T₀₀ > 0
  3. Poprawny limit newtonowski: Φ → Φ₀ + δΦ daje -GM/r
  4. Stabilność perturbacji: ω² > 0 dla wszystkich modów
- **Hipoteza:** Tylko p=1 spełnia warunki 1–4 jednocześnie

### A2c: Argument informacyjny (1–2 tygodnie)
- "Bity przestrzenne" ∝ Φ w d=3
- Entropia Bekenstein-Hawking: S_BH ∝ A/4 ∝ r² ∝ Φ (dla sferycznego Φ)
- Wymóg: ds_spatial ∝ Φ^{1/2} → g_ij ∝ Φ → p=1

## Kryterium zamknięcia

**Twierdzenie:** "Spośród g_ij = (Φ/Φ₀)^p·δ_ij, tylko p=1 daje:
ghost-free + pozytywna energia + poprawny limit newtonowski."

## Pliki do scalenia z rdzeniem

- Rozszerzenie `sek08c_metryka_z_substratu.tex`
- Skrypt weryfikacyjny → `scripts/`

## Referencje rdzenia

- `sek08c_metryka_z_substratu.tex`
- `scripts/a2_metric_consistency.py` (6/6 PASS)
- `PLAN_ROZWOJU_v3.md`, linie 53–78

## Status

- [ ] A2b: Tensor Einsteina z g_ij = Φ^p·δ_ij (ogólne p)
- [ ] A2b: Warunek ghost-free → ograniczenie na p
- [ ] A2b: Energia + Newton → jednoznaczność p=1
- [ ] A2a: Relacja dyspersji fononów
- [ ] A2c: Argument entropijny
- [ ] Skrypt weryfikacyjny
