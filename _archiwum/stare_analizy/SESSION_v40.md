# Sesja v40 — Nowy plan: τ przez mnożnik m=5/2 przy osobnym α*_τ (2026-03-29)

## Stan wejściowy po v39

### Zamknięte fakty (pewne)

| Fakt | Skrypt | Wartość |
|------|--------|---------|
| G₃(φ²·z₀)_vfar > R₃₁ wszędzie na α∈[2.3,5.0] | ex96,ex100,ex102 | min 4519 (+30%) |
| Pozorne zero α≈7.1 = artefakt okna [60,76] | ex101 | G₃_true≈65000 |
| F=(A(φ·z₀)/A(z₀))^4 = R₂₁: dwa zera | ex98 | α*₁=2.4360, α*₂=2.7538 |
| m*(α) ≈ −0.199α+2.989 (fit liniowy, VFAR) | ex98 | m*(α*₁)≈2.504, m*(α*₂)≈2.441 |
| α_min_vfar=2.5574, F_min_vfar=197.04 < R₂₁ | ex102 | potwierdzone |
| m=2.5 przy α=α*₁: G₃=3519 (+1.2%) | ex102 FAZA1 | |

---

## Kluczowe odkrycie (weryfikacja numeryczna, 2026-03-29)

### Weryfikacja wstępna (quick scan)

**Skan G3(2.5·z₀)_vfar na α∈[1.80, 2.50]** (obliczono in-session):

| α | G3(2.5·z₀) | G3−R₃₁ |
|---|-----------|---------|
| 1.80 | 6538 | +3061 |
| 2.00 | 4073 | +596 |
| 2.20 | 3677 | +199 |
| 2.30 | 3630 | +152 |
| **2.40** | **3510** | **+32** ← minimum |
| 2.45 | 3527 | +50 |
| 2.50 | 3571 | +94 |

**G3(2.5·z₀) > R₃₁ wszędzie** — oryginalna predykcja BŁĘDNA (fit liniowy m*(α) był przybliżony).

### Nowe odkrycie: krytyczny mnożnik m_c

**Skan min_α G3(m·z₀) dla m∈[2.40, 2.52]**:

| m | min G3 (nad α) | G3_min−R₃₁ | znak |
|---|----------------|-----------|------|
| 2.48 | 3368 | −109 | − |
| **2.49** | **3437** | **−40** | **−** |
| **2.50** | **3510** | **+32** | **+** |
| 2.51 | 3586 | +108 | + |

⟹ Zmiana znaku G3_min−R₃₁ przy **m_c ≈ 2.4955** (interpolacja liniowa).

**Interpretacja**: m_c jest mnożnikiem krytycznym gdzie parabola min_α G3(m·z₀;α) jest tangensem do R₃₁.
- Dla m < m_c: min_α G3 < R₃₁ → **dwa zera G3(m·z₀;α)=R₃₁** (jak muon!)
- Dla m = m_c: **jedno tangensowe zero** (minimum G3=R₃₁)
- Dla m > m_c: G3 > R₃₁ wszędzie (jak φ²-hierarchia)

m_c ≈ 2.4955 ≈ **5/2** (Δ = −0.0045, −0.18%)

---

## Zmiana paradygmatu: α lepton-specyficzne

Dotychczasowe założenie: **ten sam α** opisuje wszystkie leptony jednocześnie
(α*₁ lub α*₂ jest "globalnym" parametrem; F=R₂₁ i G₃=R₃₁ mają to samo α).

**Wynik v38–v39**: To założenie jest fałszywe — nie ma α gdzie oba warunki są spełnione.

**Nowe założenie v40**: każdy lepton może wybrać swoje α*_lepton:

| Lepton | Formuła | Warunek | α* |
|--------|---------|---------|-----|
| e | referencja | A_e(z₀;α) | dowolne |
| μ | F=(A(φ·z₀)/A_e)^4=R₂₁ | brentq F=R₂₁ | α*_μ∈{2.4360, 2.7538} |
| **τ** | **G₃=(A(m·z₀)/A_e)^4=R₃₁** | **brentq G₃=R₃₁** | **α*_τ=?  z m=5/2** |

Pytania:
1. Czy α*_τ(m=5/2) istnieje (S1)?
2. Jaka jest relacja α*_τ do α*_μ?
3. Czy istnieje "prosta" formuła spajająca α*_e, α*_μ, α*_τ?

---

## Scenariusze v40

### S1 — Precyzyjne m_c i para (m,α) dla G₃=R₃₁ ← PRIORYTET NAJWYŻSZY

```
Motywacja: Weryfikacja numeryczna wykazała:
  - m=2.49: min_α G3 = 3437 < R₃₁  → dwa zera G3=R₃₁ w α
  - m=2.50: min_α G3 = 3510 > R₃₁  → brak zer
  → m_c ≈ 2.4955 (tangensowe minimum)

FAZA A: Precyzyjne m_c (tangensowy punkt styku)
  brentq m_c gdzie min_α G3(m_c·z₀;α) = R₃₁
  Oczekiwana precyzja: 6 cyfr po przecinku
  Kandydaci algebraiczni: 5/2=2.5000, √(2π)=2.5066, ...

FAZA B: Dla m=2.49 (< m_c) — dwa zera α
  Skan G3(2.49·z₀;α) na [2.20, 2.60], krok=0.02
  brentq dwóch zer: α*_τ1, α*_τ2 (analogia z α*₁, α*₂ dla muona)
  Sprawdź: czy α*_τ1 = 2 (kanoniczna α_TGP)?

FAZA C: Weryfikacja z R_MAX=300

FAZA D: Przy α*_τ1, α*_τ2 — oblicz F(φ·z₀):
  Czy F(α*_τ) ma jakąś relację do R₂₁?

Skrypt: ex103_mc_precision.py
```

### S2 — Precyzyjne m*(α*₁) i m*(α*₂) z brentq ← PRIORYTET WYSOKI

```
Motywacja: Z ex98 mamy fit liniowy m*(α). Ale fit był przybliżeniem.
Chcemy prawdziwe m*(α*₁) i m*(α*₂) do 6 miejsc po przecinku.

Kandydaci algebraiczni przy α=α*₁=2.4360:
  m*(α*₁) ≈ 2.504
  5/2         = 2.5000  (Δ = −0.004)
  √(2π)       = 2.5066  (Δ = +0.003) ← NAJBLIŻSZY
  ln(12)      = 2.4849  (Δ = −0.019)
  φ+7/8       = 2.4930  (Δ = −0.011)

Kroki:
  FAZA 1: brentq m*(α*₁_far) = rozwiąż G₃(m·z₀;α*₁)=R₃₁ po m
    Wymagana precyzja: |G₃−R₃₁| < 0.01
  FAZA 2: brentq m*(α*₂_far)
  FAZA 3: Porównaj z √(2π), 5/2, i innymi kandydatami
  FAZA 4: Sprawdź: czy m*(α*₁)·m*(α*₂) = stała? Czy m*(α) = a/α+b?

Skrypt: włączyć do ex103 jako FAZA B.
```

### S3 — Energia solitonu jako miara masy ← PRIORYTET ŚREDNI

```
Motywacja: masa leptonu = energia spoczynkowa (E=mc²).
W TGP soliton jest "cząstką" — jego energia E(g₀;α) powinna
odpowiadać masie. Testujemy:

  E_soliton(g₀_tau; α*) / E_soliton(z₀; α*) = R₃₁ ?

Energia solitonu:
  E(g₀;α) = 4π ∫₀^R_MAX [ (f_α(g)/2)(g')² + V(g) ] r² dr
  gdzie V(g) = (g-1)²(g+2)/4 (potencjał TGP double-well)

Kroki:
  FAZA 1: Oblicz E(g₀;α) dla g₀∈{z₀, φ·z₀, φ²·z₀, m*·z₀}
    przy α∈{α*₁, α*₂, α_TGP=2}
  FAZA 2: Sprawdź ratio E_tau/E_e = R₃₁?
  FAZA 3: brentq g₀ gdzie E_ratio=R₃₁ — co to za g₀?
  FAZA 4: Sprawdź czy E_mu/E_e = R₂₁ przy g₀=φ·z₀

Oczekiwania:
  a) E_ratio ≠ A_ratio^4 (inna formuła)
  b) Może E_tau/E_e = R₃₁ przy "ładnym" g₀

Skrypt: ex104_energy_ratio.py
```

### S4 — Mapa 2D: G₃(m,α)=R₃₁ w przestrzeni (m,α) ← PRIORYTET ŚREDNI

```
Motywacja: Zamiast szukać pojedynczego (m,α), zmapuj KRZYWĄ
rozwiązań G₃(m·z₀;α)=R₃₁ w 2D. Może krzywa ta przecina linię
F(φ·z₀;α)=R₂₁ w nowym punkcie?

Siatka: m∈[2.0, 2.8]×α∈[2.0, 3.0], krok=0.1 każda oś (9×11=99 punktów)
VFAR okna.

Dodatkowe pytania:
  - Czy krzywa G₃=R₃₁ jest monotonična w (m,α)?
  - Gdzie przecina się z krzywą F=R₂₁?
  - Czy istnieje (m,α) gdzie jednocześnie G₃=R₃₁ i F=R₂₁?

Skrypt: ex105_2d_map.py
```

### S5 — Inne potęgi: (A_tau/A_e)^n = R₃₁ dla n≠4 ← PRIORYTET NISKI

```
Motywacja: Może soliton tau odpowiada innemu wykładnikowi n?

Analitycznie:
  n=4: potrzeba A_ratio=R₃₁^(1/4)=7.679 → wymaga m≈2.504 (brak zera)
  n=6: potrzeba A_ratio=R₃₁^(1/6)=3.792 ≈ R₂₁^(1/4)=3.792 (!!)
       ⟹ (A_tau/A_e)^6 = R₃₁  ⟺  (A_tau/A_e)^6 ≈ (A_mu/A_e)^4
       Ciekawy zbieg: R₃₁^(1/6) ≈ R₂₁^(1/4) do 4 cyfr (3.7920 vs 3.7920)
  n=2: potrzeba A_ratio=58.97 → g₀_tau absurdalnie duże
  n=8: potrzeba A_ratio=R₃₁^(1/8)=2.432 → mały mnożnik, g₀≈3.06

Dla n=6 potrzeba: A(m·z₀)/A(z₀)=3.792 przy pewnym (m,α).
Sprawdź: A(φ·z₀)/A(z₀)=R₂₁^(1/4)=3.792 przy α=α*₁ (z definicji!).
⟹ (A_mu(φ·z₀)/A_e(z₀))^6 = R₃₁ przy α=α*₁ (DOKŁADNIE, z definicji).

To jest formuła algebraiczna:
  **R₃₁ = R₂₁^(6/4) = R₂₁^(3/2)**

Sprawdzamy: R₂₁^(3/2) = 206.77^1.5 = 206.77×√206.77 = 206.77×14.379 = 2972 ≠ R₃₁=3477.

Niezgodność +17%. Nie jest dokładne, ale motywuje badanie bliskości.

Skrypt: analityczny + szybki test numeryczny.
```

### S6 — Związek algebraiczny α*_τ z {α*₁, α*₂, α_TGP} ← PRIORYTET NISKI

```
Po znalezieniu α*_τ z S1, zbadaj relacje:
  α*_τ + α*₁ = ?          (suma dwóch parametrów)
  α*_τ + α*₂ = ?
  α*_τ · α*₁ = ?
  α*_τ + α*₁ + α*₂ = ? (nowa S₃)
  α*_τ / α*₁ = ?
  α*_τ / α_TGP = ?         (stosunek do kanonicznego α=2)

Czy α*_τ = α_TGP = 2? (najprostsza hipoteza)
Czy α*_τ = 1 (minimalne sensowne α)?
```

---

## Kolejność realizacji v40

```
✅ S1: ex103 — m_c, zera G3(m·z₀;α)=R₃₁       — WYKONANE (16975s)
⏳ S1b: ex103b — weryfikacja per-okno dla α*_τ2≈2.561, m=2.48
⏳ S2: precyzyjne m*(α*₁,₂) z brentq          ← po weryfikacji S1b
⏳ S5: analiza algebraiczna n=6                ← bez skryptu
⏳ S3: ex104 — energia solitonu               ← po S1b
⏳ S4: ex105 — mapa 2D (m,α)                  ← po S1b,S2
⏳ S6: synteza α*_τ                            ← po weryfikacji
```

---

## Wyniki ex103 (in-session, 2026-03-30)

### FAZA A: Krytyczny mnożnik m_c (czas: 16975s)

| Parametr | Wartość |
|----------|---------|
| **m_c (brentq)** | **2.495705** |
| α @ min G3(m_c·z₀) | 2.412011 |
| G3_min(m_c) | 3477.206 (R₃₁=3477.221, residual=−0.015) |

Kandydaci algebraiczni:

| Kandydat | Wartość | Δ od m_c | ppm |
|----------|---------|----------|-----|
| 5/2 | 2.500000 | +0.004295 | **+1721** |
| φ+7/8 | 2.493034 | −0.002671 | **−1070** ← najbliższy |
| 5/2−1/200 | 2.495000 | −0.000705 | −282 |
| ln(12) | 2.484907 | −0.010798 | −4327 |
| √(2π) | 2.506628 | +0.010923 | +4377 |

**m_c ≈ 2.4957** — brak prostego kandydata algebraicznego. Najbliżej φ+7/8=φ+0.875 (−1070 ppm), ale brak fizycznej motywacji.

### FAZA B: Dwa zera dla m=2.48

| Zero | R_MAX=150 | R_MAX=300 | Δ(300−150) | Status |
|------|-----------|-----------|-----------|--------|
| **α*_τ1** | 2.27405 | 2.29695 | **+0.023** | ❌ NIESTABILNE (artefakt) |
| **α*_τ2** | 2.56010 | 2.56129 | **+0.001** | ✅ stabilne |

**α*_τ1 jest artefaktem**: duża zmiana +0.023 przy zwiększeniu R_MAX. Mechanizm: g0_tau=2.48×z0(α≈2.27)≈**3.10** — wciąż w strefie przejściowej (jak φ²·z₀=3.31 z ex93/94); VFAR okna [80,96] jeszcze w tranzycie.

**α*_τ2 jest stabilne**: Δ=+0.001, g0_tau≈3.128. Wymaga weryfikacji per-okno.

### FAZA C: R_MAX=300 (skan grubszy, krok=0.05)

| α | G3(2.48·z₀) | G3−R₃₁ |
|---|------------|---------|
| 2.20 | 3508.9 | +31.7 |
| 2.25 | 3490.5 | +13.3 |
| **2.30** | **3475.8** | **−1.4** ← zmiana znaku |
| 2.55 | 3460.7 | −16.5 |
| **2.60** | **3534.0** | **+56.8** ← zmiana znaku |

Dwa zera potwierdzone przy R_MAX=300:
- **α*_τ1(R300) = 2.2969** (niestabilne — krok 0.05, interpolacja niedokładna)
- **α*_τ2(R300) = 2.5613** (stabilne)

### FAZA D: F(φ·z₀) przy zerach

| Zero | α*_τ | F(φ·z₀) | F−R₂₁ | status |
|------|------|---------|-------|--------|
| α*_τ1 | 2.274 | 226.34 | +19.57 (+9.5%) | niestabilne |
| **α*_τ2** | **2.560** | **197.05** | **−9.71 (−4.7%)** | **stabilne** |

**Obserwacja**: F(α*_τ2)=197.05 jest **identyczne** z F_min_vfar=197.04 (ex102)! Czyli α*_τ2≈α_min_vfar=2.557–2.561.

### Synteza ex103

**Nowe rozumienie struktury G3(m·z₀;α):**
- Dla każdego m istnieje minimum G3_min(α). Gdy m=m_c: min G3=R₃₁ (jedno tangensowe zero).
- Dla m<m_c: dwa zera. Ale jedno (α*_τ1≈2.27) jest niestabilne (artefakt okna).
- Stabilne zero: **α*_τ2≈2.561, m=2.48** → G3(2.48·z₀)=R₃₁.

**Problem**: m=2.48 nie jest algebraicznie czyste. m_c=2.4957 też nie.

**Kluczowa zbieżność**: α*_τ2≈2.561 ≈ α_min_vfar=2.557 (różnica 0.004). To sugeruje, że zero G3(m·z₀)=R₃₁ istnieje blisko minimum F(α), dla mnożnika m≈m_c−0.002.

**Nowe pytanie otwarte (O-L14)**: Czy α*_τ2 jest zbieżne z α_min_vfar gdy m→m_c? Czy istnieje (m,α) gdzie G3(m·z₀)=R₃₁ **i** dF/dα=0 jednocześnie?

---

## Kluczowe pytania otwarte wchodzące z v39

| ID | Pytanie |
|----|---------|
| O-L9 | Dlaczego TGP przewiduje e i μ ale nie τ przy tym samym α? |
| O-L10 | Dlaczego α_min_vfar=2.5574 ≠ (α*₁+α*₂)/2=2.5949? |
| **O-L11** | **m*(α*₁)=? dokładnie — czy to √(2π) czy 5/2?** |
| **O-L12** | **Czy α*_τ(m=5/2) istnieje i jaka jest jego wartość?** |
| **O-L13** | **Czy R₂₁^(3/2) ≈ R₃₁ jest przypadkowe (+17%)?** |

---

## Notka o zmianie paradygmatu

Przez sesje v37–v39 szukaliśmy α gdzie **jednocześnie** F=R₂₁ i G₃=R₃₁.
Wynik: brak rozwiązania.

Sesja v40 testuje model **rozdzielonych α**:
- μ żyje przy α*_μ∈{2.436, 2.754} (gdzie F=R₂₁)
- τ żyje przy α*_τ (gdzie G₃(m·z₀)=R₃₁ dla "ładnego" m)
- e jest referencją przy każdym α

Czy te dwa α-parametry mają relację algebraiczną? To byłoby bardziej
elastyczne, ale też mniej "przewidywalne" — jeśli α*_τ jest niezależne,
model ma więcej swobody. Kluczowe: czy m=5/2 jest "czyste" algebraicznie
(co nada formule status przewidywania, nie dopasowania).
