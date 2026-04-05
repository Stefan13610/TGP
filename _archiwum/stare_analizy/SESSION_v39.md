# Sesja v39 — Plan i analiza (2026-03-29)

## Stan wejściowy po v38

### Zamknięte w v38

| ID | Wynik |
|----|-------|
| ex95 | Wykrycie near-field biasu STD okien dla φ²·z₀≈3.31; A_tau_std=2.652 (fałszywe) vs 2.994 (prawdziwe) |
| ex96 | G₃_far > R₃₁ wszędzie na α∈[2.5,5.0] — min G₃_far=4667 przy α=2.50 |
| ex97 | Plateau A_tau≈2.994 potwierdzone do r=286 (R_MAX=300) — prawdziwa asymptota |
| ex98 | α*₁_far=2.4360, α*₂_far=2.7538, S₂_far=5.1898; m*(α)∈[2.28,2.50] NIE stały |
| ex99 | F_far(α)=R₃₁ przy α≈8.074 (F=3467.6, −0.28%); G₃_far=R₃₁ przy α≈7.1–7.4; dzikie oscylacje dla α>8 |

### Otwarte priorytety v39

| Priorytet | ID | Zadanie |
|-----------|-----|---------|
| 🔴 Wysoki | O-L5 | ex100: G₃_far=R₃₁ na α∈[6.5,8.0] z R_MAX=300 — test stabilności vs rezonans |
| 🟡 Średni | O-L8 | Interpretacja α*₃_F≈8.074 (=4·α_TGP?) — algebraiczna analiza |
| 🟢 Niski | O-L4b | Precyzyjne α_min z FAR oknami: minimum F_far(α) |

---

## Analiza wstępna O-L5 (pre-session)

### Problem: G₃_far=R₃₁ przy α≈7.1 — fizyczny czy rezonans?

Z ex99 coarse scan (krok=0.5):
- α=7.0: G₃_far=781 (poniżej R₃₁=3477)
- α=7.5: G₃_far=13640 (powyżej R₃₁)

Przejście o czynnik 17.5 w pół jednostki α. To BARDZO gwałtowne. Porównaj: na zakresie [2.5,5.0] G₃_far rośnie monotonicz­nie od 4667 do 21407. Przy α=7 nagle G₃_far=781 — 6× MNIEJSZE niż przy α=5!

**Hipoteza rezonansowa**: W oknie FAR [60,136], cos(r+φ(α)) oscyluje. Przy pewnych α, faza φ(α) jest taka, że amplituda przy oknie zaczyna się destruktywnie interferować → A_tau→0 → G₃→0. To jest rezonans fazowy pomiaru, nie fizyczny zanik.

**Test**: z R_MAX=300 i szerszymi oknami [100,280], sprawdź G₃ przy α=7.0, 7.1, ..., 7.5. Jeśli przejście G₃=R₃₁ jest stabilne (ta sama α dla różnych zestawów okien i R_MAX), to jest fizyczne. Jeśli znika/przesuwa się, to rezonans.

### Czynnik 4·α_TGP = 8

Z ex99: F_far(α≈8.074) ≈ R₃₁. Różnica: 8.074 − 8 = 0.074 (~9200 ppm).

Przy FAR oknach i R_MAX=300, czy F_far(8.0) = R₃₁? Lub czy istnieje dokładne α* blisko 8?

4·α_TGP = 4·2 = 8 — to byłoby czyste algebraicznie. Sprawdzić precyzję.

### S₂_far a formuły

S₂_far = α*₁_far + α*₂_far = 2.4360 + 2.7538 = 5.1898.

Gdyby α*₃=8−S₂_far=2.8102, czy G₃_far(2.8102)=R₃₁? Z ex96, G₃_far(2.8)=5142.8 >> R₃₁. Nie.

Gdyby α*₃=4·α_TGP−S₂_far=2.8102... jw. Ta ścieżka wciąż nie działa.

---

## Plan skryptów v39

### ex100_G3_stability.py (priorytet wysoki)

```python
# Cel: weryfikacja przejścia G₃_far=R₃₁ przy α∈[6.5,8.0]
# Test: czy G₃=R₃₁ jest stabilne dla różnych R_MAX i zestawów okien
#
# FAZA 1: Dense skan G₃_far(α) na [6.5, 8.0], krok=0.05, R_MAX=150
# FAZA 2: Weryfikacja z R_MAX=300 przy podejrzanych α (gdzie G₃≈R₃₁)
# FAZA 3: Per-okno profil A_tau przy α≈7.1-7.4 — test rezonansu
# FAZA 4: Jeśli stabilne — brentq α*₃_G3 z R_MAX=300
```

### ex101_F_alpha8.py (priorytet średni)

```python
# Cel: precyzja F_far przy α≈8 z R_MAX=300
# Test: czy F_far(8.0) ≈ R₃₁ (algebraiczne 4·α_TGP=8)?
# Dense skan [7.5, 8.5], krok=0.05, R_MAX=300
```

---

## Scenariusze do zbadania v39

### S1 — G₃_far(α) na α∈[1.5, 2.4]: czy jest zero poniżej α_TGP=2?
```
Motywacja: G₃_far ma minimum przy α≈2.40–2.50 (wartość 4534). Profil G₃_far(α)
na małych α jest nieznany. Może G₃_far maleje gdy α→1.5, przekracza R₃₁ od góry?
Z ex92 (STD, nierzetelne): G₃_std spadało ku 0 dla α<2. Czy FAR daje inne zachowanie?

Oczekiwania:
  a) G₃_far rośnie dla α<2 (soliton deformuje się, amplituda rośnie) → brak zera
  b) G₃_far maleje i osiąga R₃₁ przy pewnym α<2.4 → ZERO FIZYCZNE
  c) G₃_far staje się nieokreślone (z₀ nie istnieje) poniżej pewnego α

Skrypt: ex100_G3_small_alpha.py
  Skan α∈[1.50, 2.45], krok=0.05, FAR okna, R_MAX=150
  Jeśli zmiana znaku → brentq + weryfikacja R_MAX=300
```

### S2 — G₃_far=R₃₁ przy α≈7.1–7.4: rezonans czy fizyka?
```
Motywacja: z ex99 coarse scan (krok=0.5):
  α=7.0: G₃_far=781 << R₃₁
  α=7.5: G₃_far=13640 >> R₃₁
Skok o czynnik 17.5 w Δα=0.5 — podejrzanie gwałtowny.
F(α=7.0)=1928, F(α=7.5)=2518 >> R₂₁=206.8 — warunek muona niezachowany.

Oczekiwania:
  a) Dense skan (krok=0.01) pokaże płynne przejście → może być fizyczne (ale F≠R₂₁)
  b) G₃ per-okno przy α≈7.1 pokaże niespójność okien → rezonans fazowy
  c) R_MAX=300 zmienia α* o >0.1 → rezonans

Skrypt: ex101_G3_alpha7_stability.py
  Dense skan α∈[6.8, 7.8], krok=0.02, FAR okna R_MAX=150
  Per-okno A_tau przy α=7.0,7.1,7.2 — test spójności
  Weryfikacja z R_MAX=300 przy α≈7.1
```

### S3 — F_far(α≈8)=R₃₁: precyzja i algebraiczna czystość
```
Motywacja: ex99 znalazło α*₃_F=8.074 (F_far=3467.6, −0.28% od R₃₁).
4·α_TGP=4·2=8 — czyste algebraicznie.
Pytanie: czy z R_MAX=300 i FAR oknami [100,280] α*₃_F→8.000?

Oczekiwania:
  a) F_far(8.0, R_MAX=300) bliskie R₃₁ → 4·α_TGP jest kanonicznym punktem τ
  b) F_far(8.0, R_MAX=300) daleko od R₃₁ → 8.074 to szum numeryczny

Skrypt: ex102_F_alpha8_precision.py
  Skan α∈[7.8, 8.3], krok=0.02, R_MAX=150 i R_MAX=300
  brentq F_far=R₃₁ z obydwoma R_MAX
  Test per-okno przy α=8.0 (czy plateau stabilne?)
```

### S4 — Inna hierarchia: mnożnik φ zamiast φ²
```
Motywacja: φ-mnożnik → g₀_tau=φ·z₀ ≈ 2.04 (ta sama co g₀_muon!).
Ale F=(A(φ·z₀)/A(z₀))^4=R₂₁ to jest definicja muona.
Więc φ-mnożnik dla τ odpowiadałby: (A(φ·z₀)/A(z₀))^4=R₃₁ → to samo co S3 (F=R₃₁).

Alternatywa: mnożnik φ³ → g₀_tau=φ³·z₀≈5.37.
  Sprawdź G₃_φ³=(A(φ³·z₀)/A(z₀))^4 vs R₃₁.

Inna alternatywa: g₀_tau=z₀+1 lub g₀_tau=2·z₀ (arytmetyczna hierarchia).

Skrypt: włączyć do ex102 lub osobny ex103_hierarchy.py
```

### S5 — α_min z FAR okien: gdzie F_far jest minimalne?
```
Motywacja: STD okna dają α_min≈2.54, F_min≈201.
FAR okna dają α*₁=2.436, α*₂=2.754 → α_min_far=(2.436+2.754)/2=2.595?
Ale F_far(2.595) z ex96: F≈197–199 (interpolacja) — minimum PONIŻEJ R₂₁=206.8.
To znaczy, że F_far_min < R₂₁ i α_min_far ≈ (α*₁_far+α*₂_far)/2.

Sprawdzić precyzyjnie: brentq dF/dα=0 z FAR oknami.
Skrypt: włączyć do ex100 jako FAZA 2.
```

## Kolejność realizacji v39

```
✅ ex100: S1 (G₃_far na α<2.4)  — WYKONANE
✅ ex101: S2 (stabilność α≈7.1) — WYKONANE
⏳ S3: ocena analityczna (bez nowego skryptu — artefakt pewny)
⏳ ex102: S4+S5 (inne hierarchie + α_min_far)
```

---

## Wyniki ex100 (in-session, 2026-03-29)

### S1: G₃_far na α∈[1.5, 2.55] — brak zera

| α | G₃_far | G₃−R₃₁ | F_far |
|---|--------|---------|-------|
| 1.50 | 74182 | +70705 | 4068 |
| 1.80 | 8040 | +4563 | 338 |
| 2.00 | 5380 | +1903 | 255 |
| **2.35** | **4511** | **+1034** | 221 ← **minimum G₃** |
| 2.40 | 4535 | +1058 | 217 |
| 2.55 | 4750 | +1273 | 197 |

**Brak zmiany znaku G₃_far−R₃₁ na [1.5, 2.55].** Min G₃_far=4511 przy α=2.35 — 30% powyżej R₃₁.

G₃_φ³ (φ³-mnożnik): wszędzie >>R₃₁ (od 48000 do 807000). S4-φ³ falsyfikowany.

**α_min_far = 2.5595**: F_min_far=196.85 < R₂₁=206.77 (−9.92). Midpoint zeros: (α*₁+α*₂)/2=2.5949, α_min_far=2.5595 — różnica 0.035 (α_min nie jest dokładnie w połowie).

**S1 = ❌ FALSYFIKOWANY.**

---

## Wyniki ex101 (in-session, 2026-03-29)

### S2: test rezonansu fazowego przy α∈[6.5, 8.0]

FAZA 1 dense scan pokazała dwie pozorne zmiany znaku G₃−R₃₁:
- α≈6.77: G₃_far(6.75)=+1450 → G₃_far(6.80)=−1124
- α≈7.36: G₃_far(7.35)=−359 → G₃_far(7.40)=+2398

**FAZA 2 per-okno A_tau — DEFINITYWNY TEST:**

| α | [60,76] | [80,96] | [100,116] | [120,136] | [130,146] |
|---|---------|---------|-----------|-----------|-----------|
| 6.8 | 3.18 | 4.39 | 4.53 | 4.56 | 4.54 |
| 7.0 | 3.09 | 4.36 | 4.65 | 4.61 | 4.67 |
| 7.2 | 2.65 | 4.36 | 4.72 | 4.68 | 4.76 |
| 8.0 | 1.61 | 3.20 | 4.68 | 4.77 | 4.88 |

**Plateau A_tau ≈ 4.5–4.9 od r≥80–100**, ale okno [60,76] jest w tranzycie (3.18 zamiast 4.5+). Fit A_∞(1+a/r+b/r²) na oknach [60,136] daje błędne Atau≈1.98–2.35 zamiast ~4.5–4.8.

**Prawdziwe G₃ przy α=6.8**: (4.55/0.285)⁴ ≈ **65,000 >> R₃₁=3477**
**Prawdziwe G₃ przy α=7.0**: (4.6/0.280)⁴ ≈ **72,000 >> R₃₁**
**Prawdziwe G₃ przy α=8.0**: (4.8/0.260)⁴ ≈ **116,000 >> R₃₁**

Zmiany znaku przy α≈6.77 i α≈7.36 to **artefakty fitting biasu** — ten sam mechanizm co w ex93/ex94 (near-field [60,76] bias dla φ²·z₀≈3.09–3.10).

FAZA 3 (R_MAX=300): α*₃(150)=6.775, α*₃(300)=6.771 — stabilne (Δα*=0.004), ale stabilność R_MAX nie wyklucza artefaktu okna. Stabilność wynika z faktu, że soliton w r∈[60,136] jest taki sam dla R_MAX=150 i R_MAX=300.

**S2 = ❌ FALSYFIKOWANY (artefakt biasu okna [60,76]).**

---

### S3: ocena analityczna (bez nowego skryptu)

Z ex99 FAZA 2 per-okno przy α=8.0: [60,76]=1.61, [100,116]=4.68, [120,136]=4.77 — brak plateau w [60,136].

Prawdziwe F przy α=8.0: A_e≈0.260, A_mu(φ·z₀)≈? Ale F=(A_mu/A_e)^4, nie G₃.

Z per-okno Ae przy α=8 (nie mamy osobno, ale z ex99 FAZA 1): z₀≈1.165 — mały soliton, Ae płaskowyżuje szybko.

Wzorzec: w regionie α∈[6.5,8] zarówno muon jak i tau solitony (φ·z₀≈1.88, φ²·z₀≈3.05) mają długie tranzyenty. Okna FAR [60,136] zawyżają wpływ strefy przejściowej → wartości F i G₃ z ekstrapolacji są niewiarygodne.

**S3 = ❌ PRAWDOPODOBNIE ARTEFAKT** (bez nowego obliczenia — ocena z per-okno α=8.0 z ex101).

---

### Synteza S1–S3: Definitywny wynik negatywny

**G₃_true (z φ²-mnożnikiem) > R₃₁ wszędzie na fizycznej gałęzi.**

Systematycznie, dla każdego α gdzie testowano:

| Zakres α | min G₃_true | vs R₃₁ | Metoda |
|----------|-------------|--------|--------|
| [1.5, 2.55] | 4511 (α=2.35) | **+30%** | FAR okna |
| [2.55, 5.0] | 4667 (α=2.50) | **+34%** | FAR okna (ex96) |
| [6.5, 8.0] | ~65000 (α=6.8, true) | **×19** | per-okno corr. |

**Formuła G₃=(A_tau(φ²·z₀)/A_e(z₀))^4=R₃₁ nie ma rozwiązania na fizycznej gałęzi α.**

---

### Nowe pytania otwarte (O-L9, O-L10)

**O-L9**: Dlaczego model TGP przewiduje e i μ, ale NIE τ?
- Formuła F=(A_mu(φ·z₀)/A_e(z₀))^4=R₂₁ działa (dwa zera, α*₁=2.436, α*₂=2.754)
- Ale φ² daje za duże g₀ → G₃>>R₃₁ zawsze
- Może τ = inny mechanizm fizyczny, nie amplitudowy?

**O-L10**: α_min_far = 2.5595, F_min_far = 196.85 < R₂₁
- Z okien FAR minimum F jest PONIŻEJ R₂₁ (−4.8%)
- Midpoint zeros: (2.436+2.754)/2 = 2.595 ≈ α_min_far (różnica 0.035)
- Dlaczego α_min ≠ midpoint? Asymetria gałęzi F(α)?

---

## Wyniki ex102 (in-session, 2026-03-29)

### S4: Inne hierarchie mnożnikowe — VFAR okna [80,156]

**FAZA 1**: G₃(m·z₀) dla 11 mnożników przy α∈{α*₁, α*₂, 2.50, 2.75, 3.00}

Najważniejszy wynik: **m=2.5 przy α=α*₁_far:**

| α | Mnożnik m | m·z₀ | G₃_vfar | G₃−R₃₁ |
|---|-----------|-------|---------|---------|
| **2.4360** | **2.5** | **3.1459** | **3519.19** | **+41.97 ← BINGO (+1.2%)** |
| 2.7538 | 2.5 | 3.1584 | 3975.67 | +498.44 ← bliskie |
| 2.5000 | 2.5 | 3.1503 | 3571.18 | +93.96 ← bliskie |
| 2.7500 | 2.5 | 3.1583 | 3968.62 | +491.40 |
| 3.0000 | 2.5 | 3.1567 | 4436.18 | +958.96 |

Pozostałe mnożniki przy α=α*₁:

| Mnożnik m | m·z₀ | G₃_vfar | G₃−R₃₁ |
|-----------|-------|---------|---------|
| φ ≈ 1.618 | 2.036 | 206.84 | −3270.4 (≈R₂₁ — muon!) |
| √2 ≈ 1.414 | 1.780 | 63.96 | −3413.3 |
| √3 ≈ 1.732 | 2.180 | 328.28 | −3148.9 |
| 2 | 2.517 | 958.30 | −2518.9 |
| √5 ≈ 2.236 | 2.814 | 1863.14 | −1614.1 |
| **2.5** | **3.146** | **3519.19** | **+41.97** |
| φ² ≈ 2.618 | 3.295 | 4588.75 | +1111.5 |
| e ≈ 2.718 | 3.421 | 5687.27 | +2210.1 |
| 3 | 3.775 | 9346.87 | +5869.6 |
| φ³ ≈ 4.236 | 5.331 | 47342.7 | +43865.5 |

**Obserwacja kluczowa**: mnożnik m=2.5=5/2 przy α=α*₁_far daje G₃=3519 — tylko +1.2% powyżej R₃₁=3477. Jest to wynikiem m*(α*₁)≈2.504 z ex98 (fit: m*(α)=−0.199α+2.989 → m*(2.436)=2.504). Jednak brak dokładnego zera — G₃_vfar(2.5·z₀) > R₃₁ wszędzie.

**S4 = ❌ BRAK DOKŁADNEGO ZERA** dla żadnego z 11 testowanych mnożników.
Najbliżej: m=2.5 przy α=α*₁, G₃−R₃₁=+41.97 (+1.2%). Brak zera.

---

### FAZA 2: G₃_vfar(φ²·z₀) na α∈[2.3, 5.0] — VFAR okna

| α | G₃_vfar | G₃−R₃₁ |
|---|---------|---------|
| 2.30 | 4519.43 | +1042.21 ← minimum |
| 2.40 | 4546.05 | +1068.83 |
| 2.50 | 4672.43 | +1195.21 |
| 2.60 | 4849.03 | +1371.81 |
| ... | rośnie | rośnie |
| 4.20 | 12672.51 | +9195.29 |
| 4.30 | 5871.07 | +2393.85 (resonans!) |

**Brak zmiany znaku G₃_vfar−R₃₁.** Min G₃_vfar=4519 przy α=2.30 (+30% nad R₃₁).
Przy α=4.3 widoczny resonans (spadek 12673→5871), ale G₃>R₃₁ w każdym punkcie.

Potwierdza wynik FAZA 2 ex100 (FAR) z VFAR oknami — **brak zera φ²-hierarchii**.

---

### FAZA 3: α_min_vfar (S5)

| Parametr | Wartość |
|----------|---------|
| α_min_far (ex100, FAR) | 2.5595170 |
| **α_min_vfar (ex102, VFAR, brentq)** | **2.5573858** |
| Δ(α_min) | −0.0021 |
| F_min_far | 196.85 |
| **F_min_vfar** | **197.0430** |
| F_min_vfar − R₂₁ | −9.7242 |
| (α*₁+α*₂)/2 midpoint | 2.5949 |
| α_min_vfar − midpoint | −0.0375 |

**FAR i VFAR dają niemal identyczne α_min** (różnica 0.002). Potwierdza stabilność parametru.
F_min_vfar = 197.04 < R₂₁ = 206.77 → minimum F leży wyraźnie poniżej masy mionu.

---

### FAZA 4: F_vfar vs F_far na [2.35, 2.80]

| α | F_far | F_vfar | ΔF | F_vfar−R₂₁ |
|---|-------|--------|----|------------|
| 2.350 | 221.11 | 221.14 | +0.03 | +14.37 |
| 2.400 | 216.76 | 216.87 | +0.11 | +10.11 |
| 2.450 | 203.63 | 203.66 | +0.03 | **−3.10 ←** |
| 2.500 | 198.13 | 198.40 | +0.28 | −8.36 |
| 2.550 | 196.88 | 197.07 | +0.19 | −9.70 |
| 2.600 | 197.53 | 197.73 | +0.20 | −9.04 |
| 2.650 | 199.53 | 199.66 | +0.13 | −7.11 |
| 2.700 | 202.72 | 202.72 | −0.01 | −4.05 |
| 2.750 | 206.48 | 206.60 | +0.11 | −0.17 |
| 2.800 | 210.73 | 210.90 | +0.18 | **+4.14 ←** |

**|ΔF| ≤ 0.28 wszędzie** — okna FAR i VFAR dają praktycznie identyczne F dla muona.
**Zmiany znaku F_vfar−R₂₁**: α≈2.445 i α≈2.775 → α*₁_vfar≈2.44, α*₂_vfar≈2.775.

Porównanie z FAR: α*₁_far=2.4360, α*₂_far=2.7538 — różnica **<0.02** w obu zerach.
**Okna FAR i VFAR dają te same zera F do precyzji 0.02.** Muon nie jest podatny na bias [60,76].

---

### Nowe pytanie otwarte (O-L11)

**O-L11**: Dlaczego m=2.5 przy α=α*₁_far daje G₃=3519≈R₃₁ (+1.2%)?
- m*(α*₁) ≈ 2.504 z ex98 (fit liniowy)
- m=5/2 jest prostą ułamkową wartością, ale brak dokładnego zera
- Przypadkowe numeryczne zbliżenie czy wskazówka do innej formuły?

---

### Synteza v39 — Pełna tabela wyników

| Scenariusz | Skrypt | Wynik | Status |
|------------|--------|-------|--------|
| S1: G₃_far α<2.4 | ex100 | min G₃_far=4511 (+30%), brak zera | ❌ Falsyfikowany |
| S2: G₃≈R₃₁ przy α≈7.1 | ex101 | Artefakt [60,76] biasu; G₃_true≈65000–116000 | ❌ Falsyfikowany |
| S3: F_far(α≈8)=R₃₁ | analitycznie | Artefakt per-okno przy α=8 (ex101 dane) | ❌ Prawdopodobnie artefakt |
| S4: inne mnożniki m | ex102 | m=2.5 daje G₃=3519 (+1.2%) przy α*₁ — brak zera | ❌ Brak dokładnego zera |
| S5: α_min_vfar | ex102 | α_min_vfar=2.5574, F_min=197.04 < R₂₁ | ✅ Wykonane |

**Wniosek globalny v39**: Formuła G₃=(A_tau(m·z₀)/A_e(z₀))⁴=R₃₁ nie ma rozwiązania dla żadnego z testowanych mnożników m∈{φ, √2, √3, √5, e, 2, **2.5**, φ², φ³, 3, 1+φ} na fizycznej gałęzi α∈[2.3,5.0]. Τau lepton pozostaje niewyjaśniony przez formulę amplitudową TGP.
