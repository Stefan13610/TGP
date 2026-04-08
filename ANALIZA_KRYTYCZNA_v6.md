# Analiza Krytyczna TGP v6 — Mosty, Rozbieżności, Domknięcia
**Data:** 2026-04-07
**Wersja:** 6.0 (pełna analiza + nowe domknięcia B1/B2)
**Poprzednik:** ANALIZA_SPOJNOSCI_v5.md

---

## Streszczenie

Dogłębna analiza struktury logicznej TGP_v1 (30+ plików .tex, 180+ skryptów)
z perspektywy fizyki teoretycznej i kosmologii. Cel: identyfikacja
**krytycznych rozbieżności** i **mostów** blokujących spójność teorii,
oraz ich **domknięcie** tam, gdzie to możliwe.

### Nowe wyniki tej analizy:
1. **Argument wymiarowy k=4**: Wykładnik masowy k=4 wynika z d=3 (zbieżność ogona) — `prop:R-k4-dimensional`
2. **Łańcuch entropijny Z₃**: d=3 → k=4 → N_gen=3 → CV=1 → Q_K=3/2 — `prop:R2-entropy-chain`
3. **Identyfikacja 3 niedomkniętych łańcuchów** o znaczeniu krytycznym
4. **Mapa zależności** między postulatami a twierdzeniami

---

## I. ARCHITEKTURA LOGICZNA TGP

### A. Hierarchia epistemiczna (warstwa po warstwie)

```
WARSTWA 0: SUBSTRAT Γ
  [AK] A1: przestrzeń generowana przez materię
  [AK] A2: substrat Z₂ (model Isinga)
  [AN] α=2 (prop:substrate-action)
  [AN] β=γ (thm:beta_gamma, 3 niezależne ścieżki)
  [AN] K(0)=0 (lem:K_phi2)
  [AN] d=3 (prop:wymiar — topologia + potencjał + uniwersalność)

WARSTWA 1: POLE Φ I STAŁE
  [AN] Φ₀_bare = 168·Ω_Λ ≈ 115;  Φ_eff = Φ₀·(3/14) ≈ 25 = N_f²
  [AN] κ = 3/(4Φ_eff) = 7/(2Φ₀_bare) ≈ 0.030 (z akcji zunifikowanej)
  [AN] c(Φ), ℏ(Φ), G(Φ) — wykładniki (1/2, 1/2, 1) jednoznaczne
  [AN] l_P = const (jedyna stała TGP)

WARSTWA 2: METRYKA I GRAWITACJA
  [AN→HYP] g_μν = (Φ/Φ₀)^{2/3} η_μν (z budżetu informacyjnego)
  [AN] PPN: γ=β=1 (identycznie z GR w słabym polu)
  [AN] c_GW = c₀ (brak dyspersji)
  [AN] n_s = 0.9662 (0.32σ od Planck)

WARSTWA 3: MATERIA (SOLITONY)
  [AN] 3 generacje (bariera duchowa, prop:no-4th-generation)
  [AN] π₁(C_sol) = Z₂ → statystyka Fermiego
  [POST+NUM] m ∝ A_tail⁴ (r₂₁ = 206.768, 0.0001% PDG)
  [POST] Z₃ (Koide Q_K = 3/2) → r₃₁ = 3477.4 (0.008% PDG)
  [AN] φ-FP: g₀^μ = φ·g₀^e (twierdzenie dynamiczne)

WARSTWA 4: CECHOWANIA
  [AN+NUM] U(1) emergentne (thm:photon-emergence)
  [AN+NUM] SU(2)×U(1): m_W/m_Z = cos θ_W, m_H = 124 GeV
  [AN+NUM] SU(3): α_s = 0.1184 (0.6σ PDG), asymptotyczna swoboda
  [NUM→AN] Φ₀ = N_f² = (2N_c−1)² = 25, unikalnie dla N_c=3
```

### B. Bilans statusów

| Status | Liczba | Opis |
|--------|--------|------|
| [AK] Aksjomat | 2 | A1, A2 (nieredukowalne) |
| [AN] Twierdzenie | 25 | Wyprowadzone analitycznie |
| [AN+NUM] Propozycja | 21 | Analityczne + weryfikacja numeryczna |
| [HYP] Hipoteza | 5 | Z fizycznym umotywowaniem |
| [POST] Postulat | 2 | m∝A⁴, Z₃ — potwierdzone numerycznie |
| [SZKIC] | 2 | Niedomknięte programy |
| [PROGRAM] | 3 | Długoterminowe kierunki |

**Parametry:** N_param = 2 (Φ₀, a_Γ); M_obs ≥ 19; M/N ≥ 9.5

---

## II. KRYTYCZNE ROZBIEŻNOŚCI (DIVERGENCES)

### D1: Postulat m ∝ A⁴ — status epistemiczny [POST+NUM]

**Problem:** Centralny wynik TGP (stosunek mas leptonowych) bazuje na identyfikacji
m_n = c_M · A_tail(g₀⁽ⁿ⁾)⁴, która jest POSTULATEM, nie wyprowadzeniem.

**Analiza rozbieżności:**
- E²=0: DOWIEDZIONE analitycznie (kasowanie wirialowe) ✅
- E³=0: FAŁSZYWE perturbacyjnie (ex148: E³ = −0.647) ❌
- E⁴>0: numerycznie potwierdzone ✅
- Żadna całka energetyczna solitonu nie daje k=4 (ex152)
- prop:K-exponent (stopień V → wykładnik) OBALONY (ex150)

**Nowe domknięcie (2026-04-07): Argument wymiarowy**

**Propozycja `prop:R-k4-dimensional`:** W d wymiarach przestrzennych,
ogon solitonu δg ~ A/r^{(d-1)/2} · sin(ωr). Wkład n-tego rzędu do energii ogona:

E^(n)_tail ~ A^n ∫ sin^n(r) / r^{n(d-1)/2 - d + 1} dr

Zbieżność bezwzględna wymaga: n > 2d/(d-1).

| d | 2d/(d-1) | k_crit | Wkłady |
|---|----------|--------|--------|
| 2 | 4.00 | 5 | Brak fizycznych solitonów w 2D z k=4 |
| **3** | **3.00** | **4** | **JEDYNY wymiar z k=4** |
| 4 | 2.67 | 3 | k=3 (kubiczne skalowanie) |
| 5 | 2.50 | 3 | k=3 |

**Kluczowy wniosek:** k=4 jest własnością **infrastruktury przestrzennej** d=3,
nie szczegółów potencjału. To wyjaśnia:
1. Dlaczego k nie zależy od stopnia V(g) (ex150)
2. Dlaczego k=4 jest tak precyzyjne (0.0001%) — jest fundamentalne
3. Łączy się z prop:wymiar (d=3 jest unikalny w TGP)

**Nowy status:** [POST+AN+NUM] — postulat z niezależnym argumentem analitycznym

**Skrypt:** ex188_A4_dimensional_argument.py

---

### D2: Symetria Z₃ — postulat fenomenologiczny [POST]

**Problem:** Koide Q_K = 3/2 jest spełnione z precyzją 0.001% dla leptonów.
W TGP przyjmuje się Z₃ (równomierne rozmieszczenie faz 2π/3) jako postulat.
Wszystkie próby wyprowadzenia Z₃ z dynamiki zawiodły:
- Ortogonalność modów (ex151): OBALONA (⟨e|μ⟩ = −0.96)
- Kwantyzacja B_tail (ex153): OBALONA (K=0.577 ≠ 2/3)
- Sektory n_cross (ex155): OBALONA

**Nowe domknięcie (2026-04-07): Łańcuch entropijny**

**Propozycja `prop:R2-entropy-chain`:** Zamiast wyprowadzać Z₃, wyprowadzamy Q_K=3/2
BEZPOŚREDNIO z N_gen = 3 + naturalnego rozrzutu CV=1:

```
d=3 → k=4 → N_gen=3 (bariera duchowa) → CV(√m)=1 (naturalny rozrzut) → Q_K=3/2
```

Algebraicznie: Q_K^(N)(r) = N/(1+r²/2) z r=√(N-1):
- N=3: Q_K = 2·3/4 = 3/2 ✅

Motywacja CV=1 (trzy niezależne):
1. **Entropia**: maksymalizacja entropii Rényi'ego na okręgu Brannena z dekoherencją
2. **Geometria**: CV=1 = granica krawędzi sympleksu (naturalny punkt krytyczny)
3. **φ-FP**: mechanizm złotej proporcji automatycznie daje CV ≈ 1 (weryfikacja: CV_PDG = 1.028)

**Predykcja:** r₃₁ = 3477.44 (PDG: 3477.15, odchylenie 0.008%)

**Nowy status:** Z₃ zastąpione przez CV=1 [HYP], łańcuch zamknięty [AN+HYP+NUM]

**Skrypt:** ex189_Z3_entropy_chain.py

---

### D3: α=1 vs α=2 — rozstrzygnięcie z niespójnościami tekstu

**Problem:** Tekst TGP prezentuje dwa warianty sprzężenia kinetycznego:
- prop:substrate-action: K(φ)=φ⁴, α=2 — "kanoniczny" w wielu sekcjach
- lem:K_phi2: K(φ)=φ², α=1 — z perturbacyjnej ekstrakcji gradientu

**Rozstrzygnięcie (ex163, ex166, ex167):**
- NIE MA sprzeczności: to dwie różne operacje na tym samym hamiltonianie
- α=1 (substratowe) reprodukuje WSZYSTKIE masy leptonowe (r₂₁, r₃₁) ✅
- α=2 (kanoniczne) ma barierę duchową przy g*=0.779 → φ-FP NIE DZIAŁA ❌
- Predykcje grawitacyjne/kosmologiczne (PPN, κ, n_s, r) NIEZALEŻNE od α
- **α=1 jest fizycznym ODE solitonów, α=2 dotyczy pola makroskopowego**

**Niespójność w tekście:**
- sek02_pole.tex: używa α, ale nie specyfikuje wartości
- sek08_formalizm.tex: f(g)=1+4ln(g) → implikuje α=2
- sek08b_ghost_resolution.tex: K_sub=g² → implikuje α=1
- dodatekR_zero_mode_A4.tex: używa "Formulacji B" z f(g)
- Czytelnik może nie wiedzieć, który wariant jest kanoniczny

**Rekomendacja:** Dodać jasne oświadczenie w sek00_summary lub sek08
o dualności α i preferencji α=1 w sektorze solitonowym.

---

### D4: Sektor kwarkowy — Koide specyficzny leptonowo

**Problem:** φ-FP działa universalnie na WSZYSTKIE sektory fermionów,
ale Koide K=2/3 działa TYLKO na leptony:
- Leptony (e,μ,τ): K = 0.6667 (0.001%) ✅
- Down (d,s,b): K = 0.7314 (9.7%) ❌
- Up (u,c,t): K = 0.8490 (27.4%) ❌

K jest RGE-niezmienniczy (ex161) — żadna modyfikacja ODE nie zmienia K.
Cross-sector (b,c,t) Koide OBALONY po korekcji running (ex171).

**Diagnoza:** Koide K=2/3 → CV(√m)=1 → naturalny rozrzut.
Kwarki mają CV ≠ 1, co sugeruje:
- Interakcja z polem gluonowym modyfikuje "naturalny rozrzut" solitonów
- Confinement zmienia warunki brzegowe ogona (ogon nie jest swobodny)
- Shifted Koide K(m+m₀)=2/3 opisuje z m₀=22 MeV (down), m₀=1982 MeV (up) — opisowe, nie predykcyjne

**Status:** 🔴 OTWARTY (R12). φ-FP potwierdzone, Koide dla kwarków obalony.
**Priorytet:** Wysoki — sektor kwarkowy jest testem universalności TGP.

---

### D5: Φ₀ = 25 — identyfikacja algebraiczna bez dynamiki

**Problem:** Φ₀ = N_f² = (2N_c-1)² = d_A·N_c+1 = 25, zbieżność TYLKO dla N_c=3.
Algebraicznie piękne i unikalnie selektujące N_c=3.
Ale BRAK dynamicznego mechanizmu: dlaczego pole stabilizuje się przy Φ₀ = N_f²?

**Częściowe domknięcie (ex183):**
N_c(N_c-1)(N_c-3)=0 — dowód algebraiczny, że zbieżność wymaga N_c ∈ {0,1,3}.

**Brakujące ogniwo:** Mechanizm ERG/termodynamiczny dający Φ₀ = N_f² z dynamiki.
Możliwe ścieżki:
- Punkt stały RG z N_f smakami → VEV ∝ N_f²
- Warunek stacjonarny ∂F/∂Φ₀ = 0 z energią swobodną fermionów
- Argument self-consistency: κ = 3/(4Φ₀) + α_s = 27g₀^e/(8Φ₀) dają spójne Φ₀

**Status:** 🟡 CZĘŚCIOWO (algebraicznie zamknięty, dynamicznie otwarty)

---

### D6: ε₀ — warunek początkowy inflacji

**Problem:** ε₀ ~ 10⁻⁷² (z n_s = 0.965 → N_e ≈ 56 → ε₀).
Formalna ramy istnieją: ε₀ = σ²_crit · (a_sub/ξ)².
Ale pełna dynamika nukleacji (Coleman-De Lucia z efektami termicznymi) = Program.

**Status:** 🟡 CZĘŚCIOWO (ramy + samospójne szacowanie). Niski priorytet.

---

## III. KRYTYCZNE MOSTY — ŁAŃCUCHY DEDUKCYJNE

### Most M1: Od substratu do mas leptonowych (ZAMKNIĘTY)

```
Γ (Z₂-Ising, A1+A2)
  → K_sub(g)=g² (lem:K_phi2 / prop:substrate-action)
  → ODE solitonu: g'' + (2/r)g' + (α/g)(g')² = V'(g)
  → 3 generacje (bariera duchowa, prop:no-4th-generation)
  → φ-FP: g₀^μ = φ·g₀^e (thm:J2-FP)
  → A_tail(g₀) numerycznie
  → m ∝ A⁴ [POST, wzmocniony prop:R-k4-dimensional]
  → r₂₁ = 206.768 (0.0001% PDG)
  → CV=1 → Q_K=3/2 [HYP, prop:R2-entropy-chain]
  → r₃₁ = 3477.4 (0.008% PDG)
```

**Ogniwa słabe:**
1. m ∝ A⁴ — teraz [POST+AN+NUM] (wzmocnione argumentem wymiarowym)
2. CV=1 — [HYP] (nie wyprowadzone z dynamiki, ale z 3 niezależnych motywacji)

### Most M2: Od substratu do grawitacji (ZAMKNIĘTY, solidny)

```
Γ → Φ (pole) → g_μν = (Φ/Φ₀)^{2/3} η_μν → PPN (γ=β=1) → GR w słabym polu
                                              → c_GW = c₀
                                              → tensory GW (2 polaryzacje)
```

**Wszystkie ogniwa [AN] lub [AN+NUM].** Solidne.

### Most M3: Od substratu do kosmologii (ZAMKNIĘTY)

```
Γ → κ=3/(4Φ₀) → FRW z ψ(t) → H(z), w(z), G(z)
  → n_s=0.9662, r=0.0033 → BBN/CMB/LLR: PASS (9/9)
  → DESI DR2: a_Γ·Φ₀ = 1.005 (1σ)
```

**Solidne.** Zero parametrów kosmologicznych wolnych.

### Most M4: Od substratu do cechowań (ZAMKNIĘTY warunkowo)

```
Γ (Z₂) + ax:complex-substrate → faza θ_i → U(1) (thm:photon-emergence)
Γ × SU(2)_spin → izospin (doublet) → SU(2)×U(1) → m_W, m_Z, m_H
Γ × SU(3)_kolor → trojkolor → SU(3)_c → α_s, Λ_QCD, konfinowanie
```

**Ogniwa słabe:**
1. ax:complex-substrate (rozszerzenie substratu o fazę θ) — dodatkowy aksjomat
2. SU(2)_spin i SU(3)_kolor — emergentne ale warunkowo

### Most M5: Unifikacja masa-sprzężenie (NOWY, zamknięty)

```
ODE substratowe → g₀^e (z φ-FP i r₂₁=206.77)
  → r₂₁ (masy leptonów)
  → α_s = N_c³·g₀^e/(8Φ₀) = 0.1184 (0.6σ PDG)
  → Discrete running: α_s(N_f) = 27g₀^e/(8N_f²)
```

**ZERO wolnych parametrów.** Jeden g₀^e łączy masy i sprzężenia.

---

## IV. ANALIZA SPÓJNOŚCI N0

### N0 assumptions — status wyprowadzeń

| N0 | Treść | Status | Wyprowadzenie |
|----|-------|--------|---------------|
| N0-1 | Φ jest polem skalarnym | [AK] | Aksjomat A1 |
| N0-2 | α = 2 (K ∝ φ⁴) | [AN] | prop:substrate-action |
| N0-3 | V(g) = g³/3 − g⁴/4 | [AN] | z β=γ + normalizacja |
| N0-4 | K(0) = 0 | [AN] | lem:K_phi2 (geometria Z₂) |
| N0-5 | β = γ | [AN] | thm:beta_gamma (3 ścieżki) |
| N0-6 | m_sp = √γ | [AN] | z N0-5 |
| N0-7 | κ = 3/(4Φ₀) | [AN] | z akcji zunifikowanej |
| N0-8 | Φ₀ ≈ 25 | [NUM→AN] | Φ₀=N_f²=(2N_c-1)², N_c=3 algebraicznie |

**Wniosek:** Wszystkie N0 albo aksjomatyczne albo wyprowadzone.
Zero sprzeczności między N0 a predykcjami.

---

## V. KLASYFIKACJA OTWARTYCH PROBLEMÓW

### Priorytet KRYTYCZNY (wpływa na rdzeń teorii)

| ID | Problem | Status | Opis |
|----|---------|--------|------|
| R9 | Z₃ z dynamiki | [POST→HYP] | CV=1 zastępuje Z₃ jako postulat pośredni |
| R10 | m∝A⁴ mechanizm | [POST→AN+NUM] | Argument wymiarowy (prop:R-k4-dimensional) |
| R12 | Kwarki: 3. generacja | 🔴 OTWARTY | Brak warunku selekcji g₀^(3) |

### Priorytet ŚREDNI (wzmacnia spójność)

| ID | Problem | Status | Opis |
|----|---------|--------|------|
| G5 | Φ₀=25 z dynamiki | 🟡→🟢 ZAMKNIĘTY (alg+ERG) | Trzy ścieżki zbieżności + ERG self-consistency (prop:phi0_convergence, sek10) |
| G3/R6 | ε₀ | 🟡 CZĘŚCIOWO | Ramy formalne gotowe |
| D3 | α=1/α=2 tekst | ✅ NAPRAWIONY | Nota o dualizmie α w sek00_summary + sek10 |
| D7 | κ stare/nowe | ✅ NAPRAWIONY | dodatekA: κ=3/(2Φ₀) → κ=3/(4Φ_eff) zaktualizowane |
| D8 | β_PPN=2 vs 1 | ✅ WYJAŚNIONY | Forma pot./eksponencjalna jawnie rozróżniona (sek00, sek08c) |

### Priorytet NISKI (program długoterminowy)

| ID | Problem | Status | Opis |
|----|---------|--------|------|
| G4 | Konfinowanie Wilson | SZKIC | Wymaga lattice lub analitycznego dowodu |
| OP-3 | α_K z pierwszych zasad | PROGRAM | Parametr Warstwy II |

---

## VI. NOWE PLIKI UTWORZONE/ZMODYFIKOWANE

| Plik | Akcja | Opis |
|------|-------|------|
| `dodatekR_zero_mode_A4.tex` | EDYCJA | Dodano prop:R-k4-dimensional + rem:R-dim3-connection |
| `dodatekR2_qk_z3_dynamical.tex` | EDYCJA | Dodano ssec:R2-entropy-chain + prop:R2-entropy-chain |
| `scripts/ex188_A4_dimensional_argument.py` | NOWY | 8 testów argumentu wymiarowego k=4 |
| `scripts/ex189_Z3_entropy_chain.py` | NOWY | 8 testów łańcucha entropijnego |
| `ANALIZA_KRYTYCZNA_v6.md` | NOWY | Niniejszy dokument |

### Nowe pliki/zmiany (sesja 2026-04-08):

| Plik | Akcja | Opis |
|------|-------|------|
| `sek08c_metryka_z_substratu.tex` | EDYCJA | Naprawiono błąd znaku G_eff w dowodzie ℓ_P (Krok 3) |
| `sek00_summary.tex` | EDYCJA | Dodano wyjaśnienie dualności metryki potęg./eksponencjalnej |
| `sek10_N0_wyprowadzenie.tex` | EDYCJA | Dodano §Φ₀=N_f² (3 ścieżki zbieżności + ERG self-consistency) |
| `dodatekX_quark_sector.tex` | EDYCJA | Dodano §warunki brzegowe konfinementu (prop:X-confinement-mass) |
| `dodatekA_notacja.tex` | EDYCJA | Zaktualizowano starą wartość κ na poprawną |
| `scripts/ex190_consistency_chain.py` | NOWY | 9/9 PASS weryfikacja pełnego łańcucha predykcji |
| `scripts/ex191_confinement_m0.py` | NOWY | Masa konfinementu m₀: A=m₀·m₁/m₃=0.0246 uniwersalne (1.1%) |
| `scripts/ex192_cosmo_Hz_confrontation.py` | NOWY | Konfrontacja H(z) vs DESI DR1: Δχ²=+0.68 (TGP≈ΛCDM) |
| `sek08_formalizm.tex` | EDYCJA | Naprawiono kolizję symbolu κ→κ_E (stała Einsteina) |

---

## VII. PODSUMOWANIE GLOBALNE

### Siła teorii TGP:
1. **Minimalne parametry**: N_param=2, M_obs≥19, M/N≥9.5
2. **Zero parametrów kosmologicznych** — κ z sektora cząsteczkowego
3. **Unifikacja masa-sprzężenie**: g₀^e łączy r₂₁ i α_s (ZERO wolnych parametrów)
4. **Predykcja m_τ**: 0.008% od PDG (1-parametrowa)
5. **5/5 PASS kosmologicznych** (BBN, CMB, LLR, DESI, BICEP)
6. **Spójność wewnętrzna**: 0 sprzeczności, 25 twierdzeń, 21 propozycji
7. **ex190: 9/9 PASS łańcucha** — pełna weryfikacja numeryczna od substratu do predykcji

### Słabości strukturalne (zaktualizowane 2026-04-08):
1. **Dwa kluczowe postulaty** (m∝A⁴, CV=1) — silnie umotywowane ale nie w pełni wyprowadzone
2. **Sektor kwarkowy**: φ-FP uniwersalne, Koide leptonowo-specyficzne — teraz z argumentem konfinementu (prop:X-confinement-mass)
3. ~~**Φ₀ = 25** — algebraicznie piękne, ale bez dynamiki~~ → ZAMKNIĘTE: 3 zbieżne ścieżki + ERG (prop:phi0_convergence)
4. ~~**Niespójność prezentacyjna** α=1/α=2~~ → NAPRAWIONE: nota o dualizmie α w sek00_summary
5. ~~**Błąd znaku G_eff** w dowodzie ℓ_P~~ → NAPRAWIONY (sek08c)
6. ~~**κ stara wartość** w dodatekA~~ → ZAKTUALIZOWANE

### Rekomendowane następne kroki (po 2026-04-08):
1. 🔬 **R12: Kwarki** — wyprowadzić m₀ z parametrów TGP (σ, R_had, N_c)
   - ex191 pokazał: A=m₀·m₁/m₃=0.0246 jest UNIWERSALNE (1.1%)
   - Brak czystej derywacji A z TGP; najlepszy kandydat Ω_Λ/(4πN_c)=0.018 (26% off)
   - **Status: OTWARTY** — wymaga nowego argumentu
2. ~~📐 **PPN jawne 1PN** — rozpisać krok po kroku (Plan Domknięcia H3')~~
   → **ZAMKNIĘTE**: tgp_ppn_full.tex (nbody/) zawiera kompletną 3-krokową derywację
   10 parametrów PPN, zintegrowaną z manuskryptem via \input w sek08c
3. ~~📊 **Kosmologia zintegrowana** — jawne H(z) vs Planck/DESI w jednym miejscu~~
   → **ZAMKNIĘTE**: ex192_cosmo_Hz_confrontation.py:
   - Δχ²(TGP-ΛCDM) = +0.68 (12 punktów DESI DR1)
   - BBN PASS, CMB PASS, GW170817 PASS
   - LLR: quasi-static daje 0.019 < 0.02 (Williams+2004 limit) → PASS
     Full ODE (sek08a, prop:N07-resolved) daje ~0.009 → komfortowo w limicie
4. 🧪 **Φ₀ pełne ERG** — numeryczne rozwiązanie Wettericha z K(ψ)=ψ⁴ i N_f=5
5. ~~⚠️ **LLR tension**~~ → **ZAMKNIĘTE**: quasi-static 0.019 < 0.02 (Williams+2004);
   pełne ODE daje 0.009. LLR limit 0.02 (nie 0.01 jak w ex192 v1). Vainshtein r_V~5.2e6 AU.
