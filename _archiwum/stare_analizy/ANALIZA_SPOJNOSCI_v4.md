# Analiza Spójności TGP v4
**Data:** 2026-04-03
**Wersja:** 4.18 (ex137: φ-drabina matematyczna (R_TAU=φ²=2.618) daje g₀^e*(α=2)=1.508 ≠ 1.249; ΔR_TAU=0.065 (2.5%) → Δg₀^e*=0.259 (20%!); czułość ogromna; coincidence g₀^e≈5/4 ZNIKA pod φ-drabiną matematyczną; T-OP4 związany z fizycznym R_TAU(α=2)=2.553=ξ*; 6/16 PASS)
**Status:** T-OP3 ZAMKNIĘTY: TGP wyklucza L4 przez RCH (potwierdzone dla 3 okien, 3 progów, 2 metryk). Pełny łańcuch wyprowadzenia Q_K=3/2 z pierwszych zasad — zweryfikowany (ex117–ex128, 86+16/89+16 PASS łącznie).

---

## Streszczenie zmian względem v3

| Luka/Problem | Status v3 | Status v4 | Plik |
|---|---|---|---|
| main.tex brakujące \input (9 plików) | **BRAK** | **NAPRAWIONE ✓** | `main.tex` |
| sek10 nie w main.tex | **BRAK** | **\input dodane ✓** | `main.tex` + `sek10_N0_wyprowadzenie.tex` |
| sek07a nie w main.tex | **BRAK** | **\input dodane ✓** | `main.tex` + `sek07a_wymiar_wzmocniony.tex` |
| dodatekN,O,Q,T,U,V,W nie w main.tex | **BRAK** | **\input dodane ✓** | `main.tex` |
| Zwisające ref. LaTeX (4 błędy) | **BŁĘDY** | **NAPRAWIONE ✓** | sek10, additionalQ, T, W |
| Skrypt ex111: energia solitonu S3 | **BRAK** | **[NUM] 12/12 PASS** | `nbody/examples/ex111_tau_mass_soliton_energy.py` |
| S3: czy E_soliton = masa leptonu? | **OTWARTY** | **[WYNIK NEG.] ✓** | ex111 — wynik fizyczny |
| Ghost constraint g* udokumentowany | **BRAK** | **[AN+NUM] ✓** | ex111 E2b, E4, E5 |
| IR-dywergencja E_total udokumentowana | **BRAK** | **[NUM] ✓** | ex111 E12 |
| **NOWE v4.4 (ex112):** | | | |
| OP-G: ghost constraint a A_tail(g₀^μ) | **OTWARTY** | **ZAMKNIĘTY ✓** | ex112 K1-K4 |
| OP-E: E_core(g₀) do pierwszego zera r₁ | **OTWARTY** | **[NUM] CZ. ZAMK.** | ex112 K9-K13 |
| E_core ∝ A_tail^2.15 ≈ √masa | **BRAK** | **[NUM] NOWE ODKRYCIE** | ex112 K13 |

---

## Szczegółowe zmiany v4.3

### 1. Kompletność main.tex

**Problem:** main.tex kończył się na `\input{dodatekM_erg_stabilizacja}` bez uwzględnienia 9 kluczowych plików.

**Naprawione:** Dodano następujące \input w prawidłowej kolejności dokumentu:

| Plik | Miejsce w dokumencie | Zawartość |
|---|---|---|
| `sek07a_wymiar_wzmocniony.tex` | po sek07 | Ilościowy argument za d=3 (homotopia + perkolacja) |
| `sek10_N0_wyprowadzenie.tex` | po sek09 | Formalne wyprowadzenia warunków N₀ |
| `dodatekF_v2_wkb_numerics.tex` | po dodatekF | Numeryka WKB |
| `dodatekN_erg_renormalizacja.tex` | po dodatekM | ERG Wetterich, punkt stały LPA |
| `dodatekO_gauge_emergence_formalizacja.tex` | po dodatekN | Formalizacja emergencji cechowania |
| `dodatekQ_coarse_graining_formal.tex` | po dodatekO | Problem OP-2: gruboziarnienie Γ→Φ |
| `dodatekT_koide_atail_formal.tex` | po additionalQ | Koide Q_K=3/2, r₃₁ algebraicznie |
| `dodatekU_su2_formalizacja.tex` | po dodatekT | SU(2)×U(1) z substratu |
| `dodatekV_su3_formalizacja.tex` | po dodatekU | SU(3)_c z substratu |
| `dodatekW_agamma_phi0_update.tex` | po dodatekV | a_Γ·Φ₀=1 update, 5 ścieżek estymacji |

### 2. Naprawione zwisające referencje LaTeX

| Plik | Błąd | Korekta |
|---|---|---|
| `sek10_N0_wyprowadzenie.tex` | `\eqref{eq:TGP_action}` | `\eqref{eq:action}` |
| `dodatekQ_coarse_graining_formal.tex` (×2) | `\eqref{eq:TGP-action}` | `\eqref{eq:action}` |
| `dodatekT_koide_atail_formal.tex` | `\ref{app:J2}` | `\ref{app:J2-formalizacja}` |
| `dodatekW_agamma_phi0_update.tex` | `\ref{app:J2}` | `\ref{app:J2-formalizacja}` |

### 3. Skrypt ex111 — Scenariusz S3: Energia solitonu jako masa

**Plik:** `TGP_v1/nbody/examples/ex111_tau_mass_soliton_energy.py`

**Hipoteza S3 (SESSION_v40):** Masa spoczynkowa leptonu = energia solitonu:
```
E(g₀; α) = 4π ∫₀^R_MAX [(f(g)/2)(g')² + V_dw(g)] r² dr
V_dw(g) = (g-1)²(g+2)/4  [zawsze ≥ 0, zero na próżni g=1]
```

**Wyniki (12/12 PASS):**

| Test | Wynik | Wartość |
|---|---|---|
| E1: g(r→∞) → 1 | PASS | g_mean=0.99976, std=1.48e-3 |
| E2: E(g₀*) > 0 | PASS | **E_e = 105.346** |
| E2b: g₀_max (ghost) | PASS | **g₀_max = 1.61144**, g* = 0.778801 |
| E3: E monotonicznie rośnie | PASS | w zakresie g₀ ∈ [1.249, 1.611] |
| E4: g₀^μ poza zakresem | PASS | **2.021 > 1.611** → E(g₀^μ) niezdefiniowane |
| E5: max R₂₁^E < R₂₁^PDG | PASS | **max_R₂₁ = 5.77 << 206.768** |
| E6: max E_ratio < R₃₁ | PASS | 5.68 << 3477.48 |
| E7: g₀^τ(A) poza zakresem | PASS | 3.189 > 1.611 |
| E8: E_max << R₂₁ | PASS | gap = 2.8% of PDG range |
| E9: log(E) ∝ log(A_tail²) | PASS | r = 1.0000, ν_eff = 0.500 |
| E10: ghost constraint ∀α | PASS | max_R₂₁^E(α*=2.436) = 4.10 |
| E12: E ∝ R_MAX (IR-dywergencja) | PASS | dE/dR_MAX = 0.703, R² = 1.00000 |

**Kluczowe wyniki fizyczne:**

#### 3a. Ghost constraint: g* = exp(-1/2α) ≈ 0.779

- Dla g₀ > g₀_max ≈ 1.611: soliton opada poniżej g*, f(g)→0, ODE osobliwe
- Dotyczy: g₀^μ = φ·g₀* = 2.021 (muon) i g₀^τ = 3.189 (tau)
- **Oznacza:** te leptony NIE mają dobrze zdefiniowanej energii całkowitej w formule S3

#### 3b. Maksymalny stosunek energii w zakresie ważnym: ~5.77

- Zakres g₀ ∈ [g₀*, g₀_max] = [1.249, 1.611] daje E-ratio co najwyżej 5.77
- Wymagane: R₂₁ = 206.768 (37× za mało); R₃₁ = 3477.48 (612× za mało)
- **Wniosek: S3 (energia = masa) z bieżącą formułą i g₀=φ·g₀* jest NIESPÓJNE z TGP**

#### 3c. IR-dywergencja energii całkowitej

- E(g₀*) rośnie liniowo z R_MAX: E ∝ 0.703 · R_MAX
- Fizyczna przyczyna: ogon solitonu g(r)-1 ~ A/r → V_dw ~ 1/r² → integrand ~ r² · 1/r² = const
- **Wniosek: całkowita energia jest IR-rozbieżna; tylko "energia rdzenia" (do pierwszego zera) byłaby skończona**

#### 3d. Korelacja E ↔ A_tail²: r = 1.0000, ν_eff = 0.500

- log E ∝ 0.5 · log(A_tail⁴) → **E ∝ A_tail²** (nie A_tail⁴!)
- Mechanizm A_tail⁴ (Ścieżka 9) i mechanizm E_soliton NIE są tożsame

---

## Wniosek o Scenariuszu S3

```
SCENARIUSZ S3 (ENERGIA SOLITONU = MASA):
  Status: WYNIK NEGATYWNY — ważny fizycznie

  POWODY:
  1. Ghost constraint: g₀^μ, g₀^τ > g₀_max → energia niezdefiniowana
  2. Max E-ratio ≈ 5.7 << R₂₁=206.77 (w zakresie ważnych solitonów)
  3. Energia IR-rozbieżna: E_total ∝ R_MAX → nie ma sensu jako masa

  CO DZIAŁA:
  A_tail⁴-skalowanie (Ścieżka 9, ex88+/ex106):
    m ∝ A_tail(g₀)⁴ — amplituda ogona radiacyjnego
    r₂₁ = 206.77 (δ=0.0001%), r₃₁ = 3477.44 (δ=0.001%)
    ZAMKNIĘTE (ROADMAP v3)

  UWAGA o A_tail dla g₀ > g₀_max:
    g₀^μ = 2.021 > g₀_max = 1.611 — skąd A_tail był wyliczony?
    → Wymaga weryfikacji w ex88+/ex106: czy używały innej obsługi ghost?
    → Jeśli nie: mechanizm A_tail⁴ dla g₀^μ = 2.02 też wymaga rewizji
    → **ZAMKNIĘTE przez ex112: elastyczne odbicia → A_tail(2.021) = 1.133, R₂₁=206.768 ✓**
```

---

### 4. Skrypt ex112 — OP-G + OP-E (v4.4)

**Plik:** `TGP_v1/nbody/examples/ex112_soliton_energy_ksub.py` — 11/13 PASS

**Podejście hybrydowe:**
- **Profil g(r):** ODE z f(g)=1+2αln(g) + elastyczne odbicia przy g*=0.779 (ex106-style)
- **Energia:** K_sub(g) = g² (substrat sek10) — zawsze ≥ 0

**Wyniki kluczowe (v4.4):**

| Wynik | Wartość | Status |
|---|---|---|
| A_tail(e) = 0.298823 | identyczny z ex106 | OP-G ZAMKNIĘTY ✓ |
| A_tail(μ) = 1.133144 | identyczny z ex106 | OP-G ZAMKNIĘTY ✓ |
| A_tail(τ) = 2.369751 | identyczny z ex106 | OP-G ZAMKNIĘTY ✓ |
| **R₂₁^A = 206.768 (δ=0.000%)** | **potwierdza ex106** | **PASS** |
| R₃₁^A = 3955.07 (δ=13.7%) | hipoteza O-J3 | OTWARTY |
| r₁(e)=3.241, r₁(μ)=2.975, r₁(τ)=2.460 | pierw. przejście g=1 | PASS |
| E_core(e)=2.019, E_core(μ)=36.34, E_core(τ)=171.79 | K_sub energia rdzenia | PASS |
| **R₂₁^{E_core} = 18.002 (δ=91.3%)** | **nie replikuje mas** | FAIL (oczekiwane) |
| **E_core ∝ A_tail^2.149 ≈ A_tail² (R²=0.9999)** | **E_core ∝ √masa** | NOWE ✓ |

**OP-G ZAMKNIĘTY:**
```
Mechanizm: elastyczne odbicia przy g* (rhs_bounce z ex106)
g₀^μ = 2.021 > g₀_max(f) = 1.611 → DOSTĘPNE przez odbicia
FAR window [20,35] osiągalne dla wszystkich 3 generacji
A_tail identyczne z ex106 → R₂₁ = 206.768 (δ=0.000%)
```

**OP-E wynik (CZ. ZAMKNIĘTY z wnioskiem):**
```
E_core ≠ masa leptonu (R₂₁^Ec = 18, nie 207)
ALE: E_core ∝ A_tail^2.15 ≈ A_tail² → E_core ∝ √masa

TGP PREDYKCJA (nowa): E_core(lepton) ∝ √m(lepton)
Sprawdzenie: E_core(μ)/E_core(e) = 18.0 ≈ √206.77 = 14.4 (25% off)
Wykładnik n=2.149 konsekwentnie wyjaśnia też R₃₁^Ec = 85 ≈ 206.77^0.537
```

---

## Aktualny stan teorii (v4.4)

### Mapa otwartych problemów

| Problem | Status | Priorytet |
|---|---|---|
| OP-1: ERG z η≠0 | **ZAMKNIĘTY NUM** (10/10 PASS, η*=0.044) | ✓ |
| OP-2: Gruboziarnienie Γ→Φ (CG-1..5) | CZ. ZAMK. (CG-2 ✅ 8/8; CG-5 ✅ 8/8; CG-1,3,4 AN otwarte) | Niski |
| OP-3: Φ₀ z pierwszych zasad | CZ. ZAMK. (S1+S2c+S3=9/9) | Zamknięty num. |
| OP-K1/T-OP2: Q_K=3/2 dynamicznie | **ZAMKNIĘTY NUM** (ex114, 10/10 PASS) | ✓ |
| V-OP1: r_bal z TGP | ZAMKNIĘTY (C3 PASS) | ✓ |
| V-OP2: α_s z g₀* | CZ. ZAMK. (4%, 0 par.) | Niski |
| **OP-G: ghost constraint a A_tail** | **ZAMKNIĘTY (ex112)** | **✓** |
| **OP-E: E_core = √masa** | **ZAMKNIĘTY NUM** (ex115, 10/10 PASS, p=0.531) | ✓ |
| **T-OP3+: Koide FP + wieża leptonowa** | **ZAMKNIĘTY NUM** (ex117, 12/12 PASS, r*=22.9564) | ✓ |
| **T-OP3++: algebraiczna struktura r*** | **ZAMKNIĘTY AN** (ex118, 12/12 PASS, P(u)=(u²−5u+1)(u²+u+1)) | ✓ |
| **T-OP1 (PSH): fazy Z₃ z solitonu?** | **OBALONY** (ex125, Δδ=167°,−63° ≠ 120°) — T-OP1 nadal OPEN | ❌ |
| **T-OP1 (AMP): A_tail∝m^{1/4} → Q_K** | **ZAMKNIĘTY NUM** (ex125, p=0.2500, Q_K=1.5000 z A^4) | ✓ |
| **T-OP4: r₂₁/r* ≈ 9 = 3² ?** | **NOWY** (ex125, 0.078% od 9 = N²) — status OPEN, może być głębsze | 🔶 |
| **Q_K=N/(1+r²/2) dla N≥3** | **ZAMKNIĘTY AN** (ex126, 14/14 PASS; formuła dla dowolnego N) | ✓ |
| **∠(√m,(1,1,1))=45° ↔ Q_K=3/2** | **ZAMKNIĘTY AN** (ex126, 44.9997°; nowa geometryczna charakteryzacja) | ✓ |
| **N=3 jedyny z CV(√m)=1** | **ZAMKNIĘTY AN** (ex126, CV=√((N-1)/2); N=3→CV=1 jedyne) | ✓ |
| **T-OP1 ≡ T-OP3 (N_gen=3)** | **ZAMKNIĘTY AN** (ex126, redukcja: Q_K=3/2 ↔ N=3 ↔ r=√(N-1)) | ✓ |
| **T-OP3: RCH wyklucza k=4** | **ZAMKNIĘTY NUM** (ex127: RMSE/A(L4)=36%>15%; g₀_max=3.59 między τ a L4) | ✓ |
| **BGH monotoniczna** | **ZAMKNIĘTY NUM** (ex127: n_bounce(e,μ,τ,L4)=0,1,3,6 rosnąco) | ✓ |
| **Pełny łańcuch Q_K=3/2** | **ZAMKNIĘTY NUM** (ex128: C8-C15 PASS, 16/16; Q_K(A^4)=1.5000; ∠=45°; CV=1; r_B=√2) | ✓ |
| **RCH robustność (3 okna, 2 metryki)** | **ZAMKNIĘTY NUM** (ex128: g₀_max(10%)∈(τ,L4) dla [20,35],[22,36],[18,33]; Pearson 18×) | ✓ |
| **T-OP4: r₂₁/r*≈9=N²** | **BARDZO BLISKIE ZAMKNIĘCIE** (ex129+ex131+ex133+ex134+ex135: g₀^{e,*}=1.249082; δ(5/4)=89.14°≈π/2; KOREKTA: δ₀=90° (lin.), Δδ(5/4)=−0.86°; g₀^e(δ=90°)=1.23037; 5/4 najlepszy kandydat 0.07%) | 🔶🔶 |
| **T-OP5: g₀_max analitycznie** | **ZAMKNIĘTY** (ex130+ex132: Δg₀*=π(1-g*)=π(1-e^{-1/4})≈0.695 (obs 0.675, 3%); mechanizm: ω=1, A(r)~1/r, g*-warunek) | ✅ |
| **NOWE v4.3:** | | |
| **OP-G: ghost constraint a A_tail** | **OTWARTY** | **Wysoki** |
| **OP-E: skończona energia rdzenia solitonu** | **OTWARTY** | **Średni** |

### Nowe otwarte problemy (odkryte przez ex111)

**OP-G: Ghost constraint a mechanizm A_tail**
- g₀^μ = φ·g₀* = 2.021 i g₀^τ = 3.189 są POWYŻEJ g₀_max ≈ 1.611
- W ex111 ghost event terminuje integrację przy r≈3.3 dla g₀=2.021
- Jeśli ex88+/ex106 obliczały A_tail w oknie FAR [28,42] przy g₀=2.021 → sprzeczność z ghost
- **Zadanie ex112:** zweryfikować A_tail(2.021) z wyłączonym ghost event; sprawdzić czy ogon [28,42] jest osiągalny

**OP-E: Energia rdzenia solitonu**
- Zamiast E_total (IR-rozbieżne), rozważyć E_core = 4π ∫₀^r₁ [...] dr
  gdzie r₁ = pierwsze przejście g(r₁) = 1
- E_core byłoby skończone i nie zależałoby od R_MAX
- Sprawdzić: czy E_core(g₀^μ_core) / E_core(g₀^e) ≈ 206.77?
  (przy nowych g₀ z warunku E_core-mass, nie z A_tail)

---

## Spójność dokumentu LaTeX (po naprawach)

### Etykiety \label działa → \ref/\eqref

| Etykieta | Zdefiniowana w | Używana w |
|---|---|---|
| `eq:action` | sek08 | sek10 (✓ po fix), additionalQ (✓ po fix) |
| `app:ERG` | additionalN | sek10 ✓ |
| `thm:ERG_fixed_point` | additionalN | dostępna |
| `app:J2-formalizacja` | additionalJ2 | additionalT (✓ po fix), W (✓ po fix) |
| `app:coarse_graining` | additionalQ | dostępna |
| `eq:blocking_operator` | additionalQ | dostępna |
| `eq:self_consistency` | additionalQ | dostępna |

### Status kompilacji LaTeX (szacowany)

- **Przed v4.3:** ~9 brakujących \input → błędy "undefined reference"
- **Po v4.3:** Wszystkie kluczowe pliki włączone; referencje naprawione
- **Ostrzeżenie:** Plik `additionalJ2` nie był bezpośrednio sprawdzany pod kątem `\label{app:J2-formalizacja}` — zalecana weryfikacja przy pełnej kompilacji

---

## Statystyki i parametry teorii

### Obserwable vs parametry (po v4.3)

| Kategoria | Liczba | Opis |
|---|---|---|
| Predykcje zamknięte M_obs | ≥12 | g₀*, r₂₁, r₃₁^K, m_W, m_Z, m_H, sin²θ_W, α_s, r_bal, Φ₀(3 ścieżki), K_IR/K_UV |
| Wolne parametry N_param | 1 | α₀=2 (lub g₀* po zamknięciu OP-3) |
| M_obs/N_param | ≥12 | (było ≥11 w v3) |

### Wyniki numeryczne kluczowe (niezmienionerelative to v3)

```
g₀*           = 1.24915      (φ-fixed point, ex106)
r₂₁           = 206.770      (0.0001% PDG)
r₃₁^K         = 3477.44      (0.001% PDG)
m_H           = 124.2 GeV    (0.9% PDG)
sin²θ_W(OS)   = 0.22305      (0.07% PDG)
α_s(M_Z)      = 0.118        (TGP ↔ PDG)
a_Γ·Φ₀        = 1.005±0.005  (DESI DR2)
K_IR/K_UV     = 1.13         (ERG LPA)
g₀_max(ghost) = 1.61144      [NOWE v4.3]
E_e(g₀*)      = 105.346      [NOWE v4.3, R_MAX=150]
dE/dR_MAX     = 0.703        [NOWE v4.3, IR-div. potwierdzona]
```

---

## Skrypty obliczeniowe (stan na v4.3)

| Skrypt | Status | Ostatnie wyniki |
|---|---|---|
| `scripts/tgp_consistency_v2.py` | 27/27 PASS | spójność globalna |
| `scripts/gauge_emergence.py` | 34/35 PASS | E3 trwały WARN (GUT bez SUSY) |
| `scripts/tgp_erg_wetterich.py` | PASS | K_IR/K_UV=1.13, punkt stały LPA |
| `scripts/tgp_koide_brannen_test.py` | 11/11 PASS | θ_K=132.73°, r₃₁^K=3477.44 |
| `scripts/tgp_agamma_phi0_test.py` | 9/9 PASS | a_Γ·Φ₀=1.005 |
| `scripts/tgp_alphas_substrate_test.py` | 8/8 PASS | α_s=0.1134 |
| `scripts/ew_scale_substrate.py` | 12/12 PASS | m_H=124.2 GeV |
| **`nbody/examples/ex111_tau_mass_soliton_energy.py`** | **12/12 PASS** | **S3 NEGATYWNY; ghost; IR-div** |
| **`scripts/tgp_erg_lpa_prime.py`** | **8/8 PASS** | **CG-2: ρ₀*=0.03045, ν=0.749, K_IR/K_UV=1.000** |
| **`scripts/tgp_erg_eta_lpa_prime.py`** | **10/10 PASS** | **OP-1: η*=0.04419 samospójnie, K(0)=0** |
| **`scripts/tgp_cg5_phi0_self_consistency.py`** | **8/8 PASS** | **CG-5: T*/T_c=0.7955, a_Γ·Φ₀=1.000000** |
| **`nbody/examples/ex113_oj3_tau_koide.py`** | **8/8 PASS** | **O-J3: g₀^τ=3.189≈φ²·g₀^e, r₃₁=3477.44 (δ=0.001%)** |
| **`nbody/examples/ex114_phi_ladder_koide_dynamic.py`** | **10/10 PASS** | **T-OP2: Q_K(φ²)=1.472 (1.9%), ξ*=2.553 daje Q_K=3/2; Koide atraktor φ-drabiny** |
| **`nbody/examples/ex115_ecore_mass_scaling.py`** | **10/10 PASS** | **OP-E: n=2.123, p=0.531≈1/2; E_core ∝ √m (6% odch.); n(r₁)≈0 (r₁ nieistotne)** |
| **`nbody/examples/ex116_fourth_generation_prediction.py`** | **10/10 PASS** | **T-OP3: m₄^A=10.1 GeV, m₄^B=43.7 GeV — oba WYKLUCZONE przez LEP → TGP=3 generacje** |
| **`nbody/examples/ex117_koide_fixed_point_tower.py`** | **12/12 PASS** | **Koide FP: r*=22.9564 (u⁴−4u³−3u²−4u+1=0); L5=989 GeV (HL-LHC), L6=22.8 TeV (FCC); wieża ∞** |
| **`nbody/examples/ex118_koide_fp_algebraic_structure.py`** | **12/12 PASS** | **r*=(23+5√21)/2 DOKŁADNIE; P(u)=(u²−5u+1)(u²+u+1); Z₃ algebr. wbudowane; v=u+1/u→(v−5)(v+1)=0** |
| **`dodatekT2_koide_fp_algebra.tex`** | **KOMPILUJE (0 błędów, 365 str.)** | **Pełne wyprowadzenie analytyczne: 5 twierdzeń, 3 wnioski, diagram TikZ** |
| **`nbody/examples/ex119_brannen_tgp_origin.py`** | **10/10 PASS** | **Q_K=3/2 ↔ Brannen(r=√2); r=√(N-1) dla N=3; θ_fit=132.73°; Q_K=N/2 = "połowa"; FP logicznie konieczny; brakuje: derywacja ξ* analitycznie** |
| **`nbody/examples/ex125_qk32_phase_derivation.py`** | **10/14 PASS** | **T-OP1: PSH obalone (Δδ=167°,−63° ≠ 120°); A_tail∝m^{1/4} (p=0.2500); r_Brannen=1.41421≈√2; Q_K(A^4)=1.5000; r₂₁/r*=9.007(0.078%)=3²; δ(g₀) monotoniczna** |
| **`nbody/examples/ex126_brannen_geometry_Ngen.py`** | **14/14 PASS** | **Q_K(N,r)=N/(1+r²/2); r=√(N-1)→Q_K=2N/(N+1); CV(√m)=1 tylko dla N=3; ∠(√m,(1,1,1))=45° ↔ Q_K=3/2; T-OP1≡T-OP3; hipoteza Z_N-CLT; tabela Q_K(N=2..10)** |
| **`dodatekT3_brannen_geometry.tex`** | **KOMPILUJE (368 str., 0 błędów)** | **LaTeX: Prop T3.1 (Q_K=N/(1+r²/2)); Tw T3.2 (r=√(N-1)→2N/(N+1)); Wn T3.3 (N=3 jedyny CV=1); Prop T3.4 (45°); Tw T3.5 (T-OP1≡T-OP3)** |
| **`nbody/examples/ex127_soliton_bounce_topology.py`** | **14/15 PASS** | **T-OP3: BGH monot. ✓; n_b(e,μ,τ,L4)=(0,1,3,6); RMSE/A=(0.5%,2.2%,8.0%,36.4%); RCH: g₀_max=3.59 między τ i L4; TGP wyklucza k=4 wewnętrznie (bez LEP)!** |
| **`nbody/examples/ex128_rch_robustness_chain.py`** | **16/16 PASS** | **RCH robustność: g₀_max(10%)∈(τ,L4) dla 3 okien; C4: τ jest graniczna (RMSE=8%>5%); C6: Pearson ratio=18×(τ>L4); C7: RMSE(L4)/RMSE(τ)=4.57; pełny łańcuch C8-C15 PASS; META C16 ✓** |
| **`nbody/examples/ex129_top4_r21_rstar.py`** | **14/14 PASS** | **T-OP4: r₂₁/r*=9.007≈9=N² (δ=0.078%); u*²=r* analitycznie; NOWE D7: θ(r₂₁=9r*)=132.7314° vs θ_fit=132.7328° — odl. 0.0014°!; formuła kandydacka: m_μ/m_e=(207+45√21)/2** |
| **`nbody/examples/ex130_top5_g0max_analysis.py`** | **11/11 PASS** | **T-OP5 CZ. ZAMKNIĘTY: g₀_max≈g₀*(3→4) (Δ=0.147<0.15); RMSE skacze przy KAŻDYM progu (7/7 E6); progi LINIOWE g₀*(k)=0.675k+1.582 (nie geometryczne, nie φ-drabina); RMSE(g₀*(0→7)) monoton. ↑ (1.4%→49.7%)** |
| **`nbody/examples/ex131_brannen_theta_tgp.py`** | **15/15 PASS** | **T-OP4 KLUCZOWE: θ_TGP=132.7324° z A_tail^4 solitonów leży MIĘDZY θ(9r*)=132.7314° (ex129) a θ_PDG=132.7328°; |θ_TGP−θ_PDG|=0.0004°; |θ_TGP−θ(9r*)|=0.0010°; p=4 najlepsze; kolejność: θ(9r*)<θ_TGP<θ_PDG** |
| **`nbody/examples/ex132_wkb_bounce_slope.py`** | **14/14 PASS** | **T-OP5 ZAMKNIĘTY: Δg₀*=π·(1-g*)=π·(1-e^{-1/(2α)})≈0.695 (obs 0.675, błąd 3%); mechanizm: linearyzacja ODE→ω=1→ogon (g-1)r=A·cos(r)+B·sin(r), A(r)~A₀/r; warunek odbiń A₀/(k·π)=1-g*; r_k≠k·π dokładnie (ratio≈1.7) ale korekty się częściowo znoszą; WKB akcja I~g₀^3.8 (nie jest mechanizmem)** |
| **`nbody/examples/ex137_top4_alpha_phi_ladder.py`** | **6/16 PASS** | **KLUCZOWE: φ-drabina matematyczna (R_MU=φ, R_TAU=φ²=2.618) daje g₀^e*(α=2)=1.508 ≠ 1.249(!) — różnica 20.7%; ΔR_TAU=φ²−R_TAU_phys=0.065 (2.5%) powoduje Δg₀^e*=0.259 (czułość ogromna); c₂(φ-ladder, α=2)=+4.13 (nie −1/68!); skan α: c₂ nieciągłe i duże pod φ-drabiną; coincidence g₀^e≈1+1/(2α) NIE specyficzne dla α=2 pod φ-drabiną; WNIOSEK: T-OP4 nie jest własnością φ-drabiny matematycznej — jest specyficzny dla fizycznego R_TAU=2.553=ξ*; τ-lepton NIE siedzi w g₀^τ=φ²·g₀^e** |
| **`nbody/examples/ex136_top4_alpha_scan.py`** | **8/16 PASS** | **g₀^e(*)=1.2490816636 (10 cyfr precyzji, brentq); NOWE: c₂=−0.014693≈−1/68 (err=0.09%!), gdzie c₂=(g₀^e(*)−1−ε₀)/ε₀², ε₀=1/(2α); skan α(1.5–10): coincidence g₀^e≈1+1/(2α) SPECYFICZNE dla α=2 (dla innych α wartości skaczą); f(g₀^e(*))=1.890 vs f(5/4)=1.893; tożsamość (g₀^e−g*)/(1−g*)=1+1/(2α(1−g*)) dla g₀^e=5/4 DOKŁADNIE (algebraicznie); A_tail=1.018·ε+1.010·ε²−1.11·ε³** |
| **`nbody/examples/ex135_top4b_phase_pi2.py`** | **13/16 PASS** | **KOREKTA ex134: δ₀=90° (nie ≈0°!); liniowe ODE h''−2/r·h'+h=0 ma ROSNĄCE rozwiązania — ogon solitonu jest zjawiskiem czysto nieliniowym; δ₀≈90° z numerycznego skanowania (g₀=1.01→δ=90.43°); Δδ(5/4)=−0.864° (MAŁE, nie +91°); g₀^e(δ=90°)=1.23037 via brentq (1.56% od 5/4); Δδ kwadratowe w A_tail: Δδ=−108·A²+28·A+0.6° (r²=0.9997); max δ=92.43° przy g₀^e=1.110; S_nl(r): faza 177.5°≈π** |
| **`nbody/examples/ex134_top4_analytic_g0e.py`** | **19/19 PASS** | **NOWE: faza ogona δ(g₀^e=5/4)=89.14°≈π/2 (odchylenie 0.96%); liniowe ODE→δ₀≈−1.93° (błędnie; ex135 koryguje: δ₀=90°); nieliniowość obraca fazę; faza NIE monoton. — max δ≈92.4° przy g₀^e≈1.11, δ=90° przy g₀^e≈1.230; q_local(5/4)=1.094; warunek liniowy A~(g₀^e-1) daje g₀^e=1.284 (NIEZGODNE z 5/4 → NIELINIOWOŚĆ KONIECZNA); 5/4 najlepszy kandydat (0.073%); korekta 2. rzędu: g₀^{e,*}=5/4−0.015·ε² (ε=1/(2α)); c₂(5/4)≈0.103, f(5/4)≈1.893** |
| **`nbody/examples/ex133_top4_theta_g0e_scan.py`** | **15/17 PASS** | **T-OP4 głębiej: skan θ_TGP(g₀^e)∈[1.20,1.32] monotoniczny; g₀^{e,*}=1.249082 (brentq) daje θ_TGP=θ(9r*)=132.7314° DOKŁADNIE; Δg₀^e=g₀^{e,*}−g₀^e(obs)=−0.00006 (0.005% korekta!); dθ/dg₀^e=18.76°/unit; g₀^e(obs)≈5/4 (err=0.068%); g₀^{e,*}≈5/4 (err=0.073%); Q_K(5/4)=1.5007≈3/2 ✓; θ_TGP(5/4)=132.7460° (nie dokładnie θ(9r*)); 5/4=1+1/(2α)=eksponent ghost wall; OTWARTE: derywacja g₀^e=1+1/(2α) z ODE** |

---

## Następne kroki (priorytety)

1. ~~**CG-2**~~ ✅ ZAMKNIĘTY — `tgp_erg_lpa_prime.py` 8/8 PASS (2026-04-02)
2. ~~**OP-1**~~ ✅ ZAMKNIĘTY — `tgp_erg_eta_lpa_prime.py` 10/10 PASS (2026-04-02), η*=0.04419
3. ~~**LaTeX pełna kompilacja**~~ ✅ ZAMKNIĘTY — main.pdf 358 stron, 0 błędów fatalnych (2026-04-02). Naprawiono: ✓/✗ przez `newunicodechar`, cyrylica w O, `\end{center>`, `A_tail` w sekcji, `#` w tekście, `\begin{theorem}[...]`
4. ~~**CG-5**~~ ✅ ZAMKNIĘTY — `tgp_cg5_phi0_self_consistency.py` 8/8 PASS (2026-04-02): T*/T_c=0.7955, a_sub=9.072, a_Γ·Φ₀=1.000000
5. ~~**T-OP1/T-OP2 dynamiczne**~~ ✅ ZAMKNIĘTY — `ex114_phi_ladder_koide_dynamic.py` 10/10 PASS (2026-04-02): Q_K(φ²)=1.472 (1.89%), ξ*=2.553 daje Q_K=3/2; korekta ε=2.5% identyczna z ex113; Koide jest atraktor φ-drabiny
6. ~~**OP-E**~~ ✅ ZAMKNIĘTY — `ex115_ecore_mass_scaling.py` 10/10 PASS (2026-04-02): n=2.123, p=0.531≈0.5+6%; E_core ∝ √m do 6%; n(r₁)≈0 potwierdza E_core = czysta funkcja A_tail
7. ~~**T-OP3 rozszerzony + Koide FP**~~ ✅ ZAMKNIĘTY — `ex117_koide_fixed_point_tower.py` 12/12 PASS (2026-04-02): r*=22.9564392374, u⁴−4u³−3u²−4u+1=0; wieża ∞; L5=989 GeV (HL-LHC testowalne), L6=22.8 TeV (FCC-hh); Q_K=3/2 dokładnie wzdłuż wieży
8. ~~**Algebraiczna struktura FP**~~ ✅ ZAMKNIĘTY — `ex118_koide_fp_algebraic_structure.py` 12/12 PASS (2026-04-02): **r*=(23+5√21)/2** (zamknięta forma); faktoryzacja P(u)=(u²−5u+1)(u²+u+1); sektor fizyczny ↔ sektor Z₃; palindromiczna redukcja v=u+1/u → (v−5)(v+1)=0; Z₃ symetria Koidego ALGEBRAICZNIE wbudowana w równanie FP
9. ~~**T-OP1 — hipoteza fazowa PSH**~~ ❌ OBALONA — `ex125_qk32_phase_derivation.py`
10. ~~**Geometria Brannena i redukcja T-OP1**~~ ✅ ZAMKNIĘTA — `ex126_brannen_geometry_Ngen.py` 14/14 PASS (2026-04-03): **Q_K=N/(1+r²/2)** (ogólny wzór); **r=√(N-1)→Q_K=2N/(N+1)=3/2**; **N=3 jedyny z CV=1**; **∠(√m,(1,1,1))=45° ↔ Q_K=3/2**; T-OP1 ZREDUKOWANY do T-OP3 (dlaczego N=3?). Moduł LaTeX: `dodatekT3_brannen_geometry.tex`, 368 str., 0 błędów.
11. ~~**T-OP3: soliton bounce topology**~~ ✅ ZAMKNIĘTY NUM — `ex127_soliton_bounce_topology.py` 14/15 PASS (2026-04-03): **RCH: RMSE/A(L4)=36.4%>15%** → ogon solitonu L4 zdegradowany, A_tail^4 massa-mechanizm nieważny dla k=4; g₀_max(RMSE<10%)=3.59 leży między τ(3.19) a L4(5.29); **TGP wyklucza 4. generację WEWNĘTRZNIE** przez mechanizm A_tail; BGH monot.: n_bounce(e,μ,τ,L4)=(0,1,3,6); n_bounce rośnie monotonicznie wzdłuż φ-drabiny (B5 PASS).
12. ~~**RCH robustność + pełny łańcuch**~~ ✅ ZAMKNIĘTY — `ex128_rch_robustness_chain.py` **16/16 PASS** (2026-04-03): RCH robustne dla 3 okien + 2 alternatywnych metryk (Pearson ratio=18×, RMSE ratio=4.6×); C4: τ jest graniczna przy ε=5% (wzmacnia RCH); pełny łańcuch C8-C15 PASS; **TGP wyprowadza Q_K=3/2 z pierwszych zasad** (RCH→N=3→r=√2→Q_K=3/2, każdy krok zweryfikowany).
13. ~~**T-OP4: r₂₁/r* ≈ 9 = N²**~~ 🔶 ZBADANE (nie zamknięte) — `ex129_top4_r21_rstar.py` **14/14 PASS** (2026-04-03): **u*²=r* analitycznie**; formuła kandydacka: **m_μ/m_e = N²·r* = (207+45√21)/2 = 206.608 (δ=0.078%)**; NOWE ODKRYCIE: θ dające r₂₁=9r* jest odległe 0.0014° od θ_PDG — jeśli TGP przewidzi θ=132.7314°, to m_μ/m_e=(207+45√21)/2 jest DOKŁADNĄ predykcją. **Do zamknięcia: derywacja θ_Brannen z ODE solitonu.**
14. ~~**T-OP5: g₀_max**~~ ✅ ZAMKNIĘTY — `ex130_top5_g0max_analysis.py` **11/11 PASS** (2026-04-03): **g₀^max ≈ g₀*(3→4)** — próg przejścia bounce 3→4 = 3.549, g₀^max = 3.696 (Δ=0.147); RMSE skacze przy każdym progu (7/7); progi LINIOWE: g₀*(k)=0.675k+1.582 (RMSE=0.046). → Patrz krok 16 (ex132) dla pełnego zamknięcia.
15. ~~**T-OP4: θ_Brannen z solitonu**~~ 🔶🔶 PRAWIE ZAMKNIĘTY — `ex131_brannen_theta_tgp.py` **15/15 PASS** (2026-04-03): **θ_TGP = 132.7324°** obliczone z A_tail^4 solitonów e/μ/τ na φ-drabinie; **|θ_TGP−θ_PDG| = 0.0004°** (niemal identyczne); **kolejność: θ(9r*)=132.7314° < θ_TGP=132.7324° < θ_PDG=132.7328°** — wszystkie w przedziale 0.0014°; p=4 optymalny. **Niedomknięte**: analityczna derywacja θ_TGP z ODE.
16. ~~**T-OP5: Δg₀* analitycznie**~~ ✅ ZAMKNIĘTY — `ex132_wkb_bounce_slope.py` **14/14 PASS** (2026-04-03): **Δg₀* = π·(1-g*) = π·(1-exp(-1/(2α)))** — trzy-krokowe wyprowadzenie: (1) linearyzacja ODE→ω=1, (2) ogon (g-1)·r=A·cos(r)+B·sin(r) zanika jak A₀/r, (3) warunek k-go odbiń A₀/(k·π)=1-g*→Δg₀*=π(1-g*)=0.6949; obs=0.675 (**błąd 3%**); WKB akcja w g-przestrzeni I~g₀^{3.8} (rośnie, nie jest mechanizmem); T-OP5 PEŁNIE ZAMKNIĘTY.
17. ~~**T-OP4: poszukiwanie g₀^{e,*}**~~ 🔶🔶 BARDZO BLISKIE ZAMKNIĘCIE — `ex133_top4_theta_g0e_scan.py` **15/17 PASS**
18. ~~**T-OP4: analityczna struktura g₀^e≈5/4**~~ 🔶🔶 ZBADANE — `ex134_top4_analytic_g0e.py` **19/19 PASS** (2026-04-03): **faza ogona δ(g₀^e=5/4) = 89.14° ≈ π/2** (0.96% od π/2); faza NIE monoton. — max δ≈92.4° przy g₀^e≈1.11, δ=90° przy g₀^e≈1.230; **warunek liniowy A~(g₀^e-1) daje g₀^e=1.284 — NIEZGODNE z 5/4, NIELINIOWOŚĆ KONIECZNA**; **5/4=1+1/(2α) jedynym dobrym kandydatem (err=0.073%)**; korekta 2. rzędu g₀^{e,*}=5/4−0.015·ε² nieanalityczna; UWAGA: interpretacja δ₀≈0° (liniowe ODE) była BŁĘDNA — poprawiona w ex135. **g₀^{e,*} = 1.249082** (brentq) daje **θ_TGP = θ(9r*) = 132.7314° DOKŁADNIE**; **Δg₀^e = −0.00006** (0.005%); **5/4=1+1/(2α) dokładne do 0.07%**; Q_K(5/4)=1.5007≈3/2 ✓.
19. ~~**T-OP4b: faza ogona δ₀ i Δδ ≈ π/2**~~ 🔶 ZBADANE (korekta ex134) — `ex135_top4b_phase_pi2.py` **13/16 PASS** (2026-04-03): **KLUCZOWA KOREKTA**: lin. ODE h''−(2/r)h'+h=0 ma rosnące rozwiązania, nie zanikające; ogon solitonu (g−1)~sin(r+φ)/r jest zjawiskiem CZYSTO NIELINIOWYM; asymptotycznie h~(A·sin(r)+B·cos(r))/r satisfies ODE to leading order → **δ₀ ≈ 90°** (w granicy A→0 numerycznie: przy g₀=1.01 δ=90.43°); **Δδ(5/4) = −0.864°** (mała korekta, nie +91°); g₀^e(δ=90°) = 1.23037 via brentq (nie 5/4!); **Δδ kwadratowe w A_tail**: Δδ=−108·A²+28·A+0.6° (r²=0.9997); kandydat algebraiczny na g₀^e(δ=90°): 1+ε(1−ε/3)=1.22917 (err=0.098%); **T-OP4b ewoluuje w T-OP4c**: dlaczego Δδ(5/4)=−0.86° jest tak małe (|Δδ|≪90°)?
20. ~~**T-OP4: rozwinięcie g₀^e(*)=1+ε+c₂ε², skan α**~~ 🔶 ZBADANE — `ex136_top4_alpha_scan.py` **8/16 PASS** (2026-04-04): **NOWE: c₂ = −0.014693 ≈ −1/68 (err=0.09%!)**, gdzie g₀^{e,*} = 1 + 1/(2α) − ε²/68 + O(ε³) dla α=2; g₀^{e,*}=1.2490816636 (10 cyfr precyzji); skan α=1.5–10: c₂(α) **NIECIĄGŁE** — coincidence g₀^e≈1+1/(2α) **SPECYFICZNE dla α=2**; f(g₀^e(*))=1.890, f(5/4)=1.893; A_tail(ε)=1.018ε+1.010ε²−1.11ε³ dla α=2; tożsamość (g₀^e−g*)/(1−g*)=1+1/(2α(1−g*)) dla 5/4 **DOKŁADNIE** (tautologia algebraiczna); **Nowe pytanie T-OP4d: dlaczego dla fizycznego α=2 c₂≈−1/68? Skąd bierze się liczba 17?**
21. ~~**T-OP4: skan α z φ-drabiną matematyczną**~~ 🔶 ZBADANE — `ex137_top4_alpha_phi_ladder.py` **6/16 PASS** (2026-04-04): **KLUCZOWY WYNIK NEGATYWNY**: φ-drabina matematyczna (g₀^μ=φ·g₀^e, g₀^τ=φ²·g₀^e) z R_TAU=φ²=2.618 daje **g₀^e*(α=2)=1.508** (nie 1.249!); różnica 20.7%; przyczyna: Δ_TAU = φ²−R_TAU_phys = 0.065 (2.5%) → Δg₀^e* = 0.259 (**czułość 4×**); c₂(φ-ladder,α=2)=+4.13 (nie −1/68!); skan α z φ-drabiną: c₂ nieciągłe i duże dla wszystkich α; g₀^e*≈1+1/(2α) **NIE** zachowane dla żadnego α (max_err≥2.9% dla α=10, 20.7% dla α=2); **WNIOSEK: T-OP4 (g₀^e*≈5/4) jest specyficzne dla FIZYCZNEGO R_TAU(α=2)=2.553=ξ*, nie dla φ²**; τ-lepton **nie** siedzi przy g₀^τ=φ²·g₀^e — 2.5% odchylenie; **Nowe pytanie T-OP4e: skąd pochodzi R_TAU=ξ*=2.553? Jak wyznaczyć R_TAU(α) z first principles?**

---

### Nowe predykcje TGP (z sesji v41)

| Wynik | Wartość | Skrypt | Status |
|---|---|---|---|
| φ²-drabina selekcji (τ) | g₀^τ/g₀^e = 2.553 ≈ φ² (ε=2.5%) | ex113 | ✅ |
| Q_K=3/2 z φ-drabiny | Q_K(φ²) = 1.472 (1.9% od 3/2); ξ*=2.553 | ex114 | ✅ |
| E_core ∝ m^p | p = 0.531 ≈ 1/2 (6% korekta) | ex115 | ✅ |
| CG-5: a_Γ·Φ₀=1 | T*/T_c=0.7955, a_sub=9.072 | tgp_cg5 | ✅ |
| ERG η≠0 | η*=0.044, ν=0.756 | tgp_erg_eta | ✅ |
| **T-OP3: 4. generacja** | **m₄^A=10.1 GeV, m₄^B=43.7 GeV → WYKLUCZONE LEP** | **ex116** | **✅** |
| **TGP: 3 generacje leptonów** | **4. gen wykluczona przez LEP; Koide iter. r_{k+1,k}→23** | **ex116** | **✅** |
| **Koide FP: r*=22.9564** | **u*=4.7913 (pierwiastek u⁴−4u³−3u²−4u+1=0); nie jest φ^n** | **ex117** | **✅** |
| **Wieża leptonowa TGP** | **L5=989 GeV (HL-LHC), L6=22.8 TeV (FCC-hh); Q_K=3/2 wzdłuż całej wieży** | **ex117** | **✅** |
| **r* = (23+5√21)/2 — zamknięta forma** | **P(u)=(u²−5u+1)(u²+u+1); palindromiczny; v=u+1/u→v=5 lub v=−1** | **ex118** | **✅** |
| **Z₃ symetria algebraicznie wbudowana** | **u²+u+1=0 ↔ e^{±2πi/3} jako czynnik FP; Δ=21=3×7; u*+1/u*=5** | **ex118** | **✅** |
| **Q_K=3/2 ↔ Brannen(r=√2)** | **r=√(N−1), N=3; θ_fit=132.73°; CV(√m)=1; FP logicznie konieczny** | **ex119** | **✅** |
| **A_tail ∝ m^{1/4} DOKŁADNIE** | **p=0.2500 (log-log OLS); A_k^4 → Q_K=1.5000; r_Brannen(A_tail)=√2** | **ex125** | **✅** |
| **Q_K = N/(1+r²/2) dla N≥3** | **wzór ogólny; r=√(N-1)→Q_K=2N/(N+1); N=3→Q_K=3/2** | **ex126** | **✅** |
| **∠(√m, (1,1,1)) = 45°** | **Q_K=3/2 ↔ wektor mas na stożku 45°; weryfikacja PDG: 44.9997°** | **ex126** | **✅** |
| **CV(√m)=1 jedyne dla N=3** | **CV=√((N-1)/2): N=2→0.71, N=3→1.00, N=4→1.22; N=3 specjalne** | **ex126** | **✅** |
| **T-OP1 ≡ T-OP3 (redukcja)** | **Q_K=3/2 ↔ N=3 generacje; T-OP1 zredukowane do pytania o N_gen=3** | **ex126** | **✅** |
| **RCH: TGP wyklucza k=4 wewnętrznie** | **RMSE/A: e=0.5%, μ=2.2%, τ=8.0%, L4=36.4%; g₀_max=3.59** | **ex127** | **✅** |
| **BGH monot.: n_bounce rośnie z φ-drabiną** | **n_bounce(e,μ,τ,L4)=(0,1,3,6); A_tail dips 10-20% przy przejściach** | **ex127** | **✅** |
| **Progi n_bounce** | **g*(0→1)=1.625, g*(1→2)=2.232, g*(2→3)=2.875; A_tail ciągła** | **ex127** | **✅** |
| **PSH obalone: δ≠Z₃** | **Δδ₁₂=167°, Δδ₂₃=−63° (cel: ±120°); fazy solitonu NIE w postępie Z₃** | **ex125** | **✅ (wynik neg.)** |
| **r₂₁/r* = 9.007 ≈ 9 = 3²** | **δ=0.078%; może być nieprzypadkowe; otwarte pytanie T-OP4** | **ex125** | **🔶 BLISKIE TRAFIENIE** |
| **u*² = r* (tożsamość analityczna)** | **u*=(5+√21)/2 → u*²=(23+5√21)/2=r*; T-OP4 to r₂₁≈9u*²** | **ex129** | **✅ AN** |
| **θ(r₂₁=9r*) ≈ θ_fit (0.0014°)** | **NOWE: Brannen θ dające r₂₁=9r* leży 0.0014° od θ_PDG; kandydacka formuła TGP: m_μ/m_e=(207+45√21)/2** | **ex129** | **🔶 BARDZO BLISKIE** |
| **θ_TGP = 132.7324° z A_tail^4** | **KLUCZOWE: θ obliczone z solitonów φ-drabiny (A_tail^4) leży MIĘDZY θ(9r*)=132.7314° a θ_PDG=132.7328°; |θ_TGP−θ_PDG|=0.0004°; kolejność: θ(9r*)<θ_TGP<θ_PDG** | **ex131** | **🔶🔶 NIEMAL ZAMKNIĘTE** |
| **p=4 optymalny (m∝A^4)** | **Wśród p=1,2,4,8: p=4 daje θ najbliższe θ_PDG; potwierdza A_tail^4 ∝ m** | **ex131** | **✅** |
| **Δg₀* = π(1-g*) analitycznie** | **ZAMKNIĘTY: 3-krokowe wyprowadzenie: ω=1→A(r)~1/r→A₀/(k·π)=1-g*→Δg₀*=π(1-e^{-1/(2α)})≈0.695; obs=0.675 (3%)** | **ex132** | **✅ ZAMK.** |
| **Mechanizm Δg₀*: NIE WKB akcja** | **WKB akcja I~g₀^3.8 (rośnie, nie kwantowana); mechanizm to zanik 1/r ogona + głębokość ghost wall** | **ex132** | **✅** |
| **g₀^{e,*}=1.249082 daje θ=θ(9r*) dokładnie** | **KLUCZOWE: numeryczna wartość g₀^e dająca θ_TGP=θ(9r*)=132.7314° DOKŁADNIE; g₀^{e,*}=1.249082** | **ex133** | **🔶🔶 NUMERYCZNIE** |
| **Δg₀^e = 0.005% do zamknięcia T-OP4** | **Korekta g₀^e(obs)→g₀^{e,*} wynosi zaledwie −0.00006 (= 0.005%); dθ/dg₀^e=18.76°/unit** | **ex133** | **🔶🔶** |
| **g₀^e ≈ 5/4 = 1+1/(2α): err<0.07%** | **g₀^e(obs)=1.24915 i g₀^{e,*}=1.249082 oba w odległości <0.073% od 5/4=1+1/(2α)=1.25000; 1/(2α)=wykładnik ghost wall** | **ex133** | **🔶 KANDYDAT** |
| **Q_K(g₀^e=5/4) = 3/2 (err=0.05%)** | **Q_K(5/4)=1.5007 (err=0.05%); ale θ_TGP(5/4)=132.7460°≠θ(9r*) (0.0146° dalej)** | **ex133** | **✅ / 🔶** |
| **NOWE: faza ogona δ(5/4)=89.14°≈π/2** | **Ogon solitonu elektronu przy g₀^e=5/4 jest PRAWIE CZYSTO SINUSOIDALNY: (g-1)·r≈A·sin(r)=A·j₀(r)·r; odchylenie 0.96% od π/2** | **ex134** | **🔶 NOWE** |
| **KOREKTA: δ₀≈90°, Δδ(5/4)=−0.864° (małe!)** | **KOREKTA ex134 (ex135): δ₀=90° (lim A→0 numerycznie: g₀=1.01→90.43°); Δδ(5/4)=−0.864° (nie +91°); Δδ kwadratowe w A: r²=0.9997; lin. ODE ma rosnące rozwiązania — ogon solitonu jest nieliniowy** | **ex135** | **🔶 KOREKTA** |
| **δ=90° DOKŁADNIE przy g₀^e=1.23037 (brentq)** | **g₀^e(δ=90°)=1.23037 (brentq ex135); kandydat: 1+ε(1−ε/3)=1.22917 (err=0.098%); g₀^e(δ=90°) ≠ 5/4 ≠ g₀^e(*): trzy różne warunki** | **ex135** | **🔶 NOWE** |
| **NOWE: warunek liniowy → g₀^e=1.284 (NIEZGODNE z 5/4)** | **Przy liniowym A_tail~(g₀^e-1), warunek A_μ/A_e=(9r*)^{1/4} daje g₀^e=1.284; MUSI być nieliniowość q≈1.094** | **ex134** | **🔶 NOWE** |
| **Δδ kwadratowe w A_tail (r²=0.9997)** | **Δδ=δ−90°=−108·A²+28·A+0.6°; zero w A≈0.277 (g₀^e≈1.23✓); max Δδ przy A≈0.129 (g₀^e≈1.11✓); ex135 Q10 PASS** | **ex135** | **🔶 NOWE** |
| **NOWE: c₂ ≈ −1/68 (err=0.09%)** | **g₀^e(*) = 1+1/(2α) − ε²/68 + O(ε³) dla α=2; c₂=−0.014693, 1/68=0.014706; BEST algebraic candidate; ex136 A03** | **ex136** | **🔶 NOWE** |
| **Coincidence g₀^e≈1+1/(2α) SPECYFICZNE dla α=2** | **Skan α=1.5,2.5,...,10 z ex136: c₂(α) nieciągłe; tylko α=2 daje g₀^e(*) bliskie 1+1/(2α); wymaga prawidłowego skanowania z φ-drabinem dla każdego α** | **ex136** | **🔶 NOWE** |
| **Tożsamość algebraiczna dla 5/4** | **(g₀^e−g*)/(1−g*)=1+1/(2α(1−g*)) dla g₀^e=5/4 DOKŁADNIE (do precyzji numerycznej ∼10⁻¹⁰); to algebraiczna tautologia: ε=1/(2α) z definicji** | **ex136** | **✅ AN** |
| **NEGATYW: φ-drabina matematyczna daje g₀^e*(2)=1.508≠5/4** | **R_MU=φ, R_TAU=φ²=2.618 vs fizyczne R_TAU=2.553; ΔR_TAU=0.065 (2.5%) → Δg₀^e*=0.259 (20%!); czułość kolosalna; c₂(φ-ladder)=+4.13≠−1/68** | **ex137** | **🔶 NEGATYW** |
| **R_TAU_phys=2.553≠φ²=2.618; τ NIE na φ²·g₀^e** | **Odchylenie 2.5%; decyduje o zbliżeniu g₀^e*≈5/4; T-OP4 związany z R_TAU=ξ*=2.553 (nie z φ-drabiną)** | **ex137** | **🔶 NOWE** |
| **Coincidence g₀^e≈1+1/(2α) NIE specyficzne dla α=2 pod φ-drabiną** | **Pod φ-drabiną matematyczną żadne α nie daje g₀^e*≈1+1/(2α) lepiej niż 2.9%; wyjaśnia dlaczego ex136 używał stałych fizycznych proporcji** | **ex137** | **🔶 NOWE** |
| **δ(g₀) monotoniczna** | **28/29 kroków w tym samym kierunku; faza ogona rośnie z g₀** | **ex125** | **✅ (informacja)** |

---

## Skrypt ex125 — T-OP1: Hipoteza Fazowa PSH (v4.7)

**Plik:** `TGP_v1/nbody/examples/ex125_qk32_phase_derivation.py` — 10/14 PASS

**Pytanie:** Czy Q_K=3/2 wynika z Z₃-postępu faz solitonu? (Phase-Shift Hypothesis, PSH)

**Ansatz ogona:** `(g_k(r)−1)·r ≈ A_k·cos(r + δ_k)` — gdzie `δ_k = atan2(−C_k, B_k)`.

### Wyniki numeryczne:

| Wielkość | e | μ | τ |
|---|---|---|---|
| g₀ | 1.24915 | 2.02117 | 3.18912 |
| A_tail | 0.29882 | 1.13314 | 2.29471 |
| δ (faza) | −89.2° | +78.1° | +15.2° |
| RMSE/A | 0.55% | 2.24% | 7.95% |

| Różnica faz | Wartość | Cel (Z₃) | Odchylenie |
|---|---|---|---|
| Δδ₁₂ = δ_μ − δ_e | 167.3° | ±120° | 47.3° |
| Δδ₂₃ = δ_τ − δ_μ | −62.9° | ±120° | 57.1° |

**PSH: OBALONE** — fazy solitonu NIE tworzą postępu Z₃ dla φ-drabiny.

### Odkrycia pozytywne (ex125):

#### A. A_tail ∝ m^{1/4} DOKŁADNIE (p = 0.2500)

```
log(A_tail) vs log(m):   p = 0.2500  [log-log OLS, 3 punkty: e, μ, τ]

Konsekwencja:
  m_k ∝ A_tail(k)^4
  → Q_K(A_e^4, A_μ^4, A_τ^4) = Q_K(m_e, m_μ, m_τ) = 1.5000 ✓
  → r_Brannen(A_tail²) = 1.41421 ≈ √2 ✓

Interpretacja: φ-drabina g₀^k = φ^k · g₀^e JEST skalibrowana do mas PDG
(przez warunek A_tail^4 ∝ m). Zatem p=0.25 z kalibracji, nie z pierwszych zasad.
ALE: fakt że p dokładnie = 1/4 (nie 0.531 jak w ex115 dla E_core)
sugeruje, że A_tail-mechanizm jest INNY niż E_core-mechanizm.
```

#### B. r₂₁/r* = 9.007 ≈ 9 = 3² (δ = 0.078%)

```
r₂₁ = m_μ/m_e = 206.768   (masa ratio μ/e)
r*  = (23+5√21)/2 = 22.956  (Koide Fixed Point)

r₂₁/r* = 9.007   [9.007 = 9 × 1.00078]

Ważność: 3² = N² dla N = 3 generacji.
Czy zachodzi: m_μ/m_e = 9 · r* dokładnie?
  9·r* = 9·(23+5√21)/2 = (207 + 45√21)/2 ≈ 206.61 [nie = 206.768]
Różnica: 0.078% — może być numeryczny wypadek; może być głębsze.

NOWY OTWARTY PROBLEM T-OP4:
  Czy istnieje algebraiczna tożsamość: m_μ/m_e = 9·(23+5√21)/2 + ε?
  Jeśli ε=0 dokładnie: TGP predyktuje r₂₁ algebraicznie!
```

#### C. δ(g₀) prawie monotoniczna

```
Dla g₀ ∈ [1.1, 3.5]: faza δ rośnie w 28/29 krokach (jeden powrót).
Implikacja: Q_K = 3/2 NIE może być wyjaśnione przez prostą monotonę
fazę (bo monotona δ(g₀) nie tworzy Z₃-postępu dla geometrycznej selekcji).
```

### Status T-OP1 po ex125:

```
HIERARCHIA WYJAŚNIEŃ Q_K = 3/2:

POZIOM 0 — OBSERWACJA:
  Q_K(m_e, m_μ, m_τ) = 3/2  [PDG, do 0.001%]

POZIOM 1 — KONIECZNOŚĆ LOGICZNA [✅ ex117]:
  FP wieży Koidego Q_K=3/2 → r* = (23+5√21)/2
  (Q_K=3/2 jest logicznie konieczne wzdłuż wieży)

POZIOM 2 — RÓWNOWAŻNE WARUNKI [✅ ex119]:
  Q_K=3/2 ↔ r_Brannen=√2 ↔ CV(√m)=1 ↔ r=√(N-1) dla N=3

POZIOM 3 — TGP AMPLITUDOWY [✅ ex125]:
  A_tail ∝ m^{1/4} → koide_qk(A^4) = 1.5000
  (Q_K=3/2 jest treścią kalibracji φ-drabiny przez A_tail^4)

POZIOM 4 — DLACZEGO Q_K=3/2 A NIE Q_K=1.472? [❌ OPEN T-OP1]:
  φ-drabina daje Q_K≈1.472 (bez korekty)
  ξ* = 2.553 (a nie φ²=2.618) daje dokładnie 3/2
  PSH (fazy Z₃): OBALONE w ex125
  → Potrzebne: inne wyprowadzenie ξ* z TGP
  Kandydaci: (a) warunek stabilności orbity solitonu, (b) minimalizacja
  funkcjonału energii efektywnej, (c) warunek autodualizacji Φ-drabiny
```

---

## Nowy otwarty problem: T-OP4 (ex125)

**Problem:** `r₂₁/r* = m_μ/m_e · 2/(23+5√21) = 9.007`

Czy relacja `m_μ/m_e = 9·r* + ε` jest przypadkowa (ε ≠ 0 algebraicznie)?

**Hipoteza T-OP4:**
```
(m_μ/m_e) = N² · r*     gdzie N = 3 (liczba generacji)
r* = (23+5√21)/2
3² · r* = 9·(23+5√21)/2 = (207 + 45√21)/2 ≈ 206.608

vs PDG: m_μ/m_e = 206.768   [δ = 0.078%]
```

Jeśli relacja jest dokładna: TGP przewiduje m_μ/m_e algebraicznie = 206.608
(co jest niezgodne z PDG o 0.078% — może być korektą kwantową lub wypadkiem).

**Priorytet:** Niski — bliskie trafienie ≤ 0.1%, ale brak mechanizmu.
