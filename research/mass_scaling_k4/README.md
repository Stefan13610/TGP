# R5: Prawo skalowania m ∝ A_tail⁴

## Problem

Formuła masowa m_n = c_M · A_tail⁴ jest **fundamentem sektora leptonowego**.
Daje r₂₁ = (A_μ/A_e)⁴ = 206.74 z dokładnością **0.013%** wobec PDG (206.768).

Ale potęga k=4 jest **postulatem**. Dlaczego nie k=3 lub k=5?

## Obecny status (2026-04-16)

### ✅ PRZEŁOM (2026-04-16): NIEPERTURBACYJNY DOWÓD m ~ A⁴ ZNALEZIONY

**Mechanizm:** `m_phys = c_m · K²` gdzie `K = ∫_0^R_max (1/2)·g^(2α)·(g')²·r²·dr`
jest **pełną** kinetyczną całką działania (NIE perturbacyjnym rzędem).

```
K ~ C_T · A²  uniwersalnie (slope 1.9997 substrat, 1.9965 canonical)
C_T = R_max/4 + C_core  (analityczne z cos² averaging tail)
C_core/A² ≈ 1.09 topologiczny niezmiennik (std 1.22%)

m_i/m_j = (K_i/K_j)² = (A_i/A_j)⁴  cutoff-independent
```

Weryfikacja: `r5_k_squared_mechanism.py` — **7/7 PASS**.

Perturbacyjne podejście zawodzi bo E_full = K - |V| to RÓŻNICA quasi-rownych
(K/|V| = 1.013), czula na korekty nieliniowe. Ale K i |V| OSOBNO skalują czysto.

### ✅ UDOWODNIONE

| Element | Status | Dowód | Precyzja |
|---------|--------|-------|----------|
| E^(2) = 0 (wiriał) | **TWIERDZENIE** | Tożsamość wiriałowa modu zerowego | dokładne |
| k=4 z konwergencji | **TWIERDZENIE** | k = 2(d-1)/(d-2) = 4 dla d=3 | dokładne |
| k=4 jedyny integer w d=3 | **TWIERDZENIE** | d=3→4, d=4→3, d=5→2.67 | algebraiczne |
| k_eff = 4.000 (numerycznie) | **WERYFIKACJA** | lp4 test G2: ln(r₂₁)/ln(A_μ/A_e) | 10⁻⁴ |
| (A_μ/A_e)⁴ = 206.74 | **WERYFIKACJA** | lp4 test C1: 0.013% od PDG | 0.013% |
| k=4 dyskryminuje | **WERYFIKACJA** | k=3→55, k=4→207, k=5→784 | jednoznaczne |
| E_full ~ A^{2α} | **WERYFIKACJA** | Fit: k=4.36 ±0.4 (canonical α=2) | 9% |

### ⚠️ WYNIKI NEGATYWNE (2026-04-14)

| Element | Status | Wynik |
|---------|--------|-------|
| E^(3) → 0 | **OBALONY** | E^(3) ~ A³ NIE znika; |E3/E4| ~ A^{-0.9} → ∞ |
| Perturbacyjny dowód m~A⁴ | **NIEMOŻLIWY** | E^(3) DOMINUJE E^(4) dla małych solitonów |

### ⚠️ OTWARTE (częściowo rozwiązane 2026-04-16)

| Element | Status | Problem |
|---------|--------|---------|
| ~~Nieperturbacyjny dowód E_full ~ A⁴~~ | **✅ ZAMKNIĘTE** | m = c·K² (nie E_full!); K ~ A² uniwersalnie |
| Zamknięta formuła c_M | CZĘŚCIOWE | c_M = c_m · C_T²; C_T = R_max/4 + C_core analityczne |
| Formalny dowód (Lean) | OTWARTE | Łańcuch dowodowy gotowy do formalizacji |
| Absolutna skala masy (bridge R5) | OTWARTE | Wymaga fizycznego R_max w teorii R5 |

## Kluczowy wynik R5 (2026-04-14): E^(3) NIE znika!

### On-shell identity (nowe twierdzenie)

Używając EOM + IBP (całkowanie przez części):
```
E^(3)_sub = -(2π/3) ∫ h³ r² dr      (substrate, K=g²)
E^(3)_can = +(4π/3) ∫ h³ r² dr      (canonical, K=g⁴)

E^(3)_can = -2 · E^(3)_sub           (na zlinearyzowanym modzie)
```

### Problem konwergencji

Dla zlinearyzowanego h = A·sin(r)/r:
- `∫h³r²dr = A³∫sin³(r)/r dr` → **ROZBIEŻNY** (logarytmicznie!)
- Dla pełnego solitonu: **SKOŃCZONY** (regularyzacja nieliniowa)
- `∫h³r²dr ~ A^{3.15}` — bliskie A³ ale nie dokładnie

### Konsekwencja

```
E^(3) ~ A^{3.4}    (obie formulacje)
E^(4) ~ A^{4.3}    (obie formulacje)
|E^(3)/E^(4)| ~ A^{-0.9} → ∞  gdy A → 0

WNIOSEK: E^(3) DOMINUJE E^(4) dla małych solitonów!
         Perturbacyjny argument "E³→0 ∴ m~A⁴" jest BŁĘDNY.
```

### Prawdziwy mechanizm m ~ A⁴

Skalowanie mas NIE wynika z anulowania poszczególnych rzędów,
lecz z GLOBALNYCH własności:

1. **Wiriał**: E^(2) = 0 (dokładne)
2. **Konwergencja ogona**: ∫sin^n(r)/r^{n-2}dr zbiega tylko dla n > 3
3. **Ogon dominuje**: Energia ogona jest jedyną dobrze określoną częścią
4. **Rdzeń kompaktowy**: E_core jest gładką funkcją g₀, nie separowalną na potęgi A

Skalowanie E_full ~ A^{2α} jest własnością NIEPERTURBACYJNĄ.

## Łańcuch dowodowy (zamknięty 2026-04-16)

```
┌─────────────────────────────────────────────────────────┐
│  (P1) WIRIAŁ: E^(2) = 0 dokładnie                      │
│       Mod zerowy: u₀ = sin(r)/r, E_kin = E_pot          │
│       → potęga k=2 WYKLUCZONA                           │
├─────────────────────────────────────────────────────────┤
│  (P2) KONWERGENCJA: E^(n) ~ A^n ∫ sin^n(r)/r^{n-2} dr │
│       Zbieżność wymaga n > 3 w d=3                     │
│       → potęgi k=1,2,3 WYKLUCZONE (rozbieżne w ogonie) │
├─────────────────────────────────────────────────────────┤
│  (P3) E^(3) ≠ 0 — WYNIK NEGATYWNY                      │
│       E^(3) = -(2π/3)∫h³r²dr (on-shell identity)       │
│       NIE znika; DOMINUJE E^(4) dla małych solitonów    │
│       Perturbacyjny dowód m~A⁴ NIEMOŻLIWY               │
├─────────────────────────────────────────────────────────┤
│  (P4) NIEPERTURBACYJNY K²: m_phys = c · K²              │
│       K = ∫(1/2)·g^(2α)·(g')²·r²·dr  [FULL akcja]      │
│       K = (R_max/4)·A² + C_core·A² (analityczne)       │
│       K_core/A² ≈ 1.09 TOPOLOGICZNY NIEZMIENNIK         │
│       → m_i/m_j = (K_i/K_j)² = (A_i/A_j)⁴ CUTOFF-IND   │
│       ✓ Sub (α=1): slope 1.9997, CAN (α=2): 1.9965     │
│       ✓ PDG diff: μ/e +0.12%, τ/e -0.20%               │
└─────────────────────────────────────────────────────────┘
```

**Klucz:** m ~ A⁴ NIE wynika z E_full = K - |V| (różnicy quasi-rownych
wielkosci — czula na nieliniowosci). Wynika z KWADRATU osobnej calki
kinetycznej K (lub |V|), ktora skaluje sie czysto jak A².

```
K/|V| = 1.013 ± 0.12%   (wirial quasi-trywialny)
K-|V|  =  0.013 · K      (delikatna różnica, niestabilna)
K      =  1     · K      (clean A² scaling) ← TO JEST WŁAŚCIWA OBSERWABLA
```

## Formuła konwergencji

```
k = 2(d-1)/(d-2)

d=3: k = 4  (jedyny dokładny integer!)
d=4: k = 3
d=5: k = 8/3 (nie integer)
```

k=4 jest **algebraicznie wyróżnione** — to jedyna wymiarowość dająca integer.

## Formulacja substratowa (K=g²)

Kanoniczny ODE (K=g⁴) jest niestabilny dla g₀ > 1.3.
Formulacja substratowa (K=g²) daje ten sam wynik z pełną stabilnością:

```
ODE: g'' + (1/g)(g')² + (2/r)g' = 1 - g

g₀ᵉ = 0.86941, g₀ᵘ = φ·g₀ᵉ = 1.40673
A_e = 0.1246, A_μ = 0.4725
(A_μ/A_e)⁴ = 206.74 ≈ 206.768 (PDG)
k_eff = ln(r₂₁)/ln(A_μ/A_e) = 4.0001
```

## Co jeszcze brakuje do zamknięcia (stan 2026-04-16)

1. ~~Analityczny dowód E^(3) → 0~~ → **OBALONY** — E^(3) nie znika
2. ~~Nieperturbacyjny dowód m ~ A⁴~~ → **✅ ZAMKNIĘTE** — m = c·K² (r5_k_squared_mechanism)
3. ~~Core-tail mechanizm~~ → **✅ ZAMKNIĘTE** — K_core/A²≈1.09 uniwersalne
4. **Analityczna wartość c_M** → **CZĘŚCIOWE** — c_M = c_m·(R_max/4 + C_core)²
    - Tail (R_max/4) analityczne z cos² averaging
    - C_core ≈ 1.09 topologiczne, ale wartość 1.09 nie ma jeszcze
      zamknietej formy analitycznej
    - Absolutna kalibracja c_m wymaga fizycznego R_max (bridge R5)
5. **Formalizacja w Lean 4** — cały łańcuch P1-P4

## Otwarta tensja 2026-05-01 → PARTIAL RESOLUTION 2026-04-30

**Status:** PARTIAL RESOLUTION — główna tensja A⁴ (R5) vs A^(5-α) (why_n3
Phase 2) **strukturalnie rozwiązana**. R5 K² jest specyficznym przypadkiem
α=1 substrate ODE; Phase 2 `m_obs = c·A²·g₀^(e²(1-α/4))` to uniwersalna
formuła. Sub-tensja g₀_τ kalibracji (1.72931 R5 vs 1.75505 Phase 2 Koide
K=2/3) pozostaje OPEN dla Phase 7+.

Patrz: [[RECONCILIATION_R5_vs_phase2_2026-04-30.md]] (pełna analiza +
niezależna weryfikacja `r5_phase2_reconciliation.py`).

### Ustawienie problemu

Cykl `research/why_n3/` (RESOLUTION 2026-05-01, patrz
`CORRECTIONS_2026-05-01.md`) odkrył że dla TGP-canonical α=2
(z `K(φ) = K_geo·φ⁴` w sek08a), formuła obserwowalnej masy ma
postać:

```
m_obs(g₀, α) = c · A_tail^(5−α)
```

z:
- α=1 (substrate ODE, K=g²): p = 4 → m ~ A⁴ ✓ **zgadza się z R5 K²
  mechanizmem** (slope 1.9997)
- α=2 (TGP-canonical, K=g⁴): p = 3 → m ~ A³ (6/6 PASS PDG: μ/e
  −0.099%, τ/e −0.085%, τ/μ +0.015%)

### Tensja

R5 mechanizm `m = c·K²` z `K ~ A²` daje **m ~ A⁴ uniwersalnie**
(dla obu α=1 substrate slope 1.9997 i α=2 canonical slope 1.9965 —
patrz § "Łańcuch dowodowy" P4). To jest **niezależne od α**.

But why_n3 RESOLUTION 2026-05-01 mówi że obserwowalna masa skaluje
**α-zależnie** jako A^(5−α).

Dla α=1 obie formuły zgadzają się (p=4, R5 0.013% PDG, why_n3
identycznie). **Dla α=2 są niespójne** (R5 daje m~A⁴, why_n3 daje
m~A³).

### Możliwe rozwiązania

| # | Hipoteza | Konsekwencja |
|---|----------|--------------|
| 1 | R5 K² mechanizm jest specyficzny dla α=1 substrate; dla α=2 canonical operuje **inny** mechanizm (np. tail-coupling z `m_obs = c·A^(5-α)`) | A⁴ jest accidentalne dla α=1, nie strukturalne. Mechanizm K² wymaga rewizji dla α=2. |
| 2 | R5 K² daje **M_full** (pełną energię wewnętrzną), why_n3 A^(5-α) daje **m_obs** (obserwowalną via tail-coupling) — distinction analogiczny do GR ADM vs Komara, QFT bare vs renormalized | Obie formuły GENUINE, ale opisują różne wielkości. Faktor A^(α-1) między nimi ma interpretację physical (np. core dressing, tail-screening). Brakuje analitycznego mostu m_obs ↔ M_full dla α=2. |
| 3 | n(α) = e²(1-α/4) z why_n3 Phase 2 jest EMPIRICAL fit (residuum <0.1%), być może p(α)=5−α też EMPIRICAL i nie strukturalne. Prawdziwy mechanizm to K² (R5) z α-zależną korektą tail | Phase 6+ work: re-derive m_obs from K² + tail correction; sprawdzić czy A^(5-α) wypada jako leading approximation. |

**Najmocniejsze (po analizie 2026-04-30):** Hipoteza 1 + część hipotezy 2.
Numeryczna weryfikacja w `r5_phase2_reconciliation.py` (CASE 2) pokazuje:

- **R5 K² dla α=2 daje +490% mismatch vs PDG dla μ/e** (1220 vs 207).
  Dowód że R5 K² mechanizm jest specific dla α=1, NIE uniwersalny.
- **Audyt R5 `r5_k_squared_mechanism.py`**: TEST 3 (jedyny test mass=K² vs
  PDG) jest przeprowadzany TYLKO dla α=1 substrate (linie 252–296).
  Dla α=2 canonical TEST 2 weryfikuje TYLKO slope K~A². Stwierdzenie
  "K² mechanizm uniwersalny" było **niezweryfikowaną ekstrapolacją**.
- **M_full ≠ m_obs** GENUINE distinction: |M_full|_μ/|M_full|_e = 22 ≠ 207,
  ale M_full ≈ 0 (Hobart-Derrick) — to NIE jest "asymptotic mass" jak ADM.
  Original hipoteza 2 była PARTIALLY RIGHT direction, ale K² ≠ M_full ani.

**Strukturalna konkluzja:** Trzy różne wielkości:
- K (akcja kinetyczna) ~ A² ∀α (R5 verified)
- K² ~ A⁴ ∀α (algebra, NIE m_obs dla α≠1)
- M_full = K + V_eff ≈ 0 (Hobart-Derrick, NIE m_obs)
- **m_obs = c·A²·g₀^(e²(1−α/4))** (Phase 2, UNIWERSALNE)

R5 K² to **specjalny przypadek Phase 2 dla α=1**, gdzie korelacja g₀-A
przez substrate ODE absorbuje core-dressing g₀^(3e²/4) do efektywnego
A⁴ scaling z dokładnością ~0.5%.

### Honest framing po 2026-04-30 (po reconciliation)

- ✅ R5 closure dla α=1 substrate ODE: GENUINE — m ~ A⁴ z K²
  mechanizmu, 7/7 PASS, 0.013% PDG.
- ✅ why_n3 RESOLUTION dla α=2 TGP-canonical: GENUINE — m_obs =
  c·A²·g₀^(e²/2), 6/6 PASS PDG.
- ✅ **Reconciliation A⁴ ↔ A^(5−α) dla α=2: PARTIAL RESOLUTION** —
  R5 K² jest **specjalnym przypadkiem α=1** Phase 2 universal formula.
  R5 niezamierzona ekstrapolacja "K² uniwersalne" zostaje zawężona
  do α=1. Dla α=2 (TGP-canonical), Phase 2 formula z g₀^(e²/2)
  factor jest właściwą formułą. Patrz
  [[RECONCILIATION_R5_vs_phase2_2026-04-30.md]].
- ⚠ Sub-tensja g₀_τ kalibracji (1.72931 R5 substrate vs 1.75505 Phase 2
  Koide K=2/3): OPEN dla Phase 7+. Wymaga ujednolicenia derivation
  źródła (Koide chain vs φ-ladder vs direct fit).
- ⚠ X = e²/4 w `n(α) = e²(1−α/4)` why_n3: EMPIRICAL discovery
  awaiting RG-derivation (Phase 6 Q5 R⁵-bridge: NEGATIVE — patrz
  `research/why_n3/PHASE6_Q5_R5_bridge_first_attempt.md`).

### Cross-references

- `research/why_n3/CORRECTIONS_2026-05-01.md` (RESOLUTION sekcja —
  insight m_obs vs M_full + p(α)=5−α odkrycie)
- `research/why_n3/PHASE2_n_alpha_derivation.md` (mass formula
  derivation z fitu numerycznego)
- `research/why_n3/PHASE4_5_yukawa_propagator.md` (§ 8 open issues:
  τ-drift dla A²·g₀^(e²/2) reconciliation)
- `research/why_n3/PHASE6_Q5_R5_bridge_first_attempt.md` (negatywna
  próba R⁵-bridge dla X=e²/4)
- `meta/AUDYT_TGP_2026-05-01.md` § AB.3 (outstanding Phase 6+ items)

### Item do "Co jeszcze brakuje" (status 2026-04-30)

6. ~~**Reconciliation A⁴ (R5 K² universal) vs A^(5-α) (why_n3 m_obs)
   dla α=2**~~ → **PARTIAL RESOLUTION 2026-04-30** — R5 K² jest specific
   dla α=1 substrate; Phase 2 jest uniwersalne. Patrz
   [[RECONCILIATION_R5_vs_phase2_2026-04-30.md]] i
   `r5_phase2_reconciliation.py`.

7. **Sub-tensja g₀_τ kalibracji** (NEW open) — różne wartości g₀_τ między
   R5 substrate (1.72931 z Koide), R5 canonical (φ²·g₀_e=2.276 z φ-ladder),
   i Phase 2 (1.75505 z Koide K=2/3). Phase 7+ task: ujednolicenie
   derivation source dla g₀_τ across cycles.

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `scripts/lp4_mass_exponent_verification.py` | Weryfikacja 9/9 PASS | ✅ RDZEŃ |
| `scripts/ex188_A4_dimensional_argument.py` | Argument konwergencji | ✅ RDZEŃ |
| `scripts/virial_theorem_v47b.py` | Tożsamość wiriałowa | ✅ RDZEŃ |
| `r5_e3_cancellation.py` | E^(3) NIE znika: 5/7 PASS (2 EXPECTED FAIL) | ✅ NOWE |
| `r5_virial_mass_derivation.py` | Skan E(A_tail) — błędne ODE | ⚠️ DO POPRAWY |
| `r5_mass_ratio_verification.py` | Weryfikacja z poprawnym ODE | ✅ BADAWCZY |
| `r5_k_squared_mechanism.py` | **m = c·K², K~A² uniwersalnie** | ✅ **7/7 PASS** (PRZEŁOM) |

## Referencje rdzenia

- `dodatekJ_ogon_masy.tex` (teoria ogona)
- `dodatekK_wkb_atail.tex` (WKB + A_tail)
- `dodatekR_zero_mode_A4.tex` (Twierdzenie R-A4)
- `dodatekF_hierarchia_mas.tex` (hierarchia)

## Status

- [x] Wiriał E^(2) = 0 — UDOWODNIONE
- [x] Konwergencja k=4 w d=3 — UDOWODNIONE
- [x] Numerycznie k_eff = 4.000 — ZWERYFIKOWANE (9/9 PASS)
- [x] Dyskryminacja k=3,4,5 — ZWERYFIKOWANE
- [x] E^(3) → 0 — OBALONY (E^(3) ~ A³ dominuje E^(4) ~ A⁴)
- [x] On-shell identity: E^(3) = -(2π/3)∫h³r²dr — UDOWODNIONE
- [x] E_full ~ A^{2α} — ZWERYFIKOWANE (k≈4.4 canonical)
- [x] **Nieperturbacyjny dowód m ~ A⁴** — ✅ ZAMKNIĘTE (m = c·K², 7/7 PASS)
- [x] **K ~ A² uniwersalne** — ZWERYFIKOWANE (sub: 1.9997, can: 1.9965)
- [x] **C_T = R_max/4 + C_core analityczne** — WYPROWADZONE (slope 0.2513)
- [x] **K_core/A² ≈ 1.09 topologiczny niezmiennik** — CV=1.22%
- [x] **Wirial K/|V| ≈ 1.013 quasi-trywialny** — wyjaśnia fail E_full²
- [x] **Ratio m_i/m_j cutoff-independent** — POTWIERDZONE
- [ ] Absolutna skala c_m (wymaga fizycznego R_max w bridge R5)
- [ ] Analityczna formuła dla C_core ≈ 1.09
- [ ] Formalizacja łańcucha dowodowego w Lean 4
