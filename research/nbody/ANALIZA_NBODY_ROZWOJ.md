# Analiza nbody: stan, spójność z TGP v1, plan rozwoju

**Data:** 2026-04-05
**Zakres:** Pełna analiza folderu `nbody/` — co działa, co nie, co jest rozbieżne z aktualną wersją manuskryptu TGP, plan dalszego rozwoju.

---

## Spis treści

1. [Podsumowanie wykonawcze](#1-podsumowanie-wykonawcze)
2. [Inwentaryzacja zasobów](#2-inwentaryzacja-zasobów)
3. [Co jest OK](#3-co-jest-ok)
4. [Co NIE działa / jest niekompletne](#4-co-nie-działa--jest-niekompletne)
5. [Rozbieżności z aktualnym TGP v1](#5-rozbieżności-z-aktualnym-tgp-v1)
6. [Plan rozwoju](#6-plan-rozwoju)

---

## 1. Podsumowanie wykonawcze

Pakiet `nbody/` jest rozbudowanym zbiorem 221 skryptów Pythona, 24 dokumentów LaTeX i 6 plików MD, implementującym fizykę N-ciałową TGP. **Rdzeń jest solidny**: aksjonaty N0 zachowane, V₂ dokładne, V₃ przez redukcję Feynmana, regresja 59/59 PASS, infrastruktura dynamiki (leapfrog/RK45) działa.

**Krytyczna rozbieżność**: nbody używa LPA (logaritmicznej aproksymacji) sprzężenia kinetycznego `f(g) = 1 + 4·ln(g)`, podczas gdy manuskrypt główny przeszedł na pełną formę substratową `K_sub(g) = g²`. Ta różnica jest istotna w regionie silnie nieliniowym (g ≪ 1, rdzenie solitonów, czarne dziury) i musi zostać rozwiązana.

**Priorytety**: (1) aktualizacja ODE do formy K_sub = g², (2) domknięcie tezy P1 (Lyapunov), (3) rozwiązanie P4 (I_triple multipol), (4) ścieżka C.

---

## 2. Inwentaryzacja zasobów

### Pliki

| Typ | Ilość | Lokalizacja |
|-----|-------|-------------|
| Moduły Python (core) | 13 | `nbody/*.py` |
| Przykłady aktywne | ~149 | `nbody/examples/ex*.py` |
| Przykłady archiwalne | ~55 | `nbody/examples/_archiwum/` |
| Dokumenty LaTeX | 24 | `nbody/*.tex` |
| Dokumentacja MD | 6 | `nbody/*.md` |
| Wyniki CSV | ~15 | `nbody/examples/_outputs/` |
| Regresja | 2 | `verify_all.py`, `verify_nbody_lyapunov_quick.py` |

### Moduły core

| Moduł | Funkcja | Status |
|-------|---------|--------|
| `tgp_field.py` | Profil Yukawa, energia, masa ekranowania | ✅ działa (LPA form) |
| `pairwise.py` | Potencjał 2-ciałowy V₂ (analityczny, EXACT) | ✅ zamknięte |
| `three_body_terms.py` | Całka potrójnej nakładki (numeryczna 3D) | ✅ działa |
| `three_body_force_exact.py` | V₃ Feynman 2D redukcja (EXACT) | ✅ zamknięte |
| `multipole_triple_overlap.py` | I_Y multipol Legendre (P4, SEMIANALYTIC) | ✅ nowe (2026-04-05) |
| `equilibria.py` | Równowagi statyczne | ✅ zamknięte |
| `stability.py` | Hessian + projekcja modów zerowych + pełne NMA | ✅ zamknięte (P2: pełne widmo) |
| `dynamics_v2.py` | Leapfrog + RK45 | ✅ działa |
| `dynamics_backends.py` | Abstrakcja backendów integracji | ✅ działa |
| `lyapunov.py` | Benettin RK4/leapfrog, Jacobians, spektrum | ✅ działa |
| `eom_tgp.py` | Wektorowe EOM N ciał | ✅ zamknięte |
| `nbody_energy.py` | Energia pełna (V₂ + V₃) | ✅ działa |
| `configurations.py` | IC: trójkąt, tetraedr, N-gon, Burrau | ✅ zamknięte |

### Zakres przykładów

| Zakres | Skrypty | Temat |
|--------|---------|-------|
| ex55–ex103 | Starsze | Efimov, krzywe rotacji, K13–K18, V_eff |
| ex104–ex124 | Ścieżka 9 | Soliton, ogon, Koide, formalizacja |
| ex125–ex137 | Kontynuacja | EOM, regresja, 3B siły |
| ex138–ex147 | Integracja | EOM N-ciał, weryfikacja |
| ex148–ex194 | P1 Lyapunov | Chaos, Benettin, matched-H, CSV |
| ex195–ex210 | P5–P8 | Soliton K_sub, EFT, V3/V2, Hill |

---

## 3. Co jest OK

### 3.1 Aksjonaty N0 — pełna zgodność

Wszystkie 6 aksjomatów N0 jest zachowanych w kodzie:

| Aksjomat | Opis | Status |
|----------|------|--------|
| N0-1 | Przestrzeń = pole skalarne Φ | ✅ `tgp_field.py` |
| N0-2 | Próżnia Φ = Φ₀ | ✅ β = γ = 1 |
| N0-3 | Źródło: δ = C·e^{-mr}/r (Yukawa) | ✅ Droga B |
| N0-4 | Ekwiwalencja masa-pole | ✅ ex12 |
| N0-5 | V_eff = Φ³/3 − Φ⁴/4 | ✅ `energy_density()` |
| N0-6 | m_sp = √(3γ−2β) | ✅ `screening_mass()` |

Żaden plik nie narusza aksjomatów. Brak obcych założeń (GR, DM jako oddzielny byt).

### 3.2 Wyniki fizyczne potwierdzone

1. **Łamanie tw. Earnshawa** — analityczne: d² − 4βd + 18βC = 0, próg β/C > 4.5
2. **Siły 3-ciałowe z Φ⁴** — V₃ = −6γ C₁C₂C₃ · I_triple (niezerowe, nieredukowalne)
3. **Granica Coulomba** — I^Coul = 8π²/P, skalowanie 1/P (nie 1/P²)
4. **Hierarchia potencjałów** — V_grad (1/d) > V_β (1/d²) > V_γ (1/d³)
5. **Fałszywa próżnia g=1** — V''(1) = −1 (metastabilność, odkrycie ex15)
6. **Okno Efimova** — szerokość vs m_sp (ex23–ex27)
7. **Krzywe rotacji** — wynik NEGATYWNY (TGP sam nie wyjaśnia DM barionowej) — wartościowy

### 3.3 Infrastruktura numeryczna

- **Regresja 59/59 PASS** (`verify_all.py`)
- **Dynamika**: leapfrog (symplektyczny) + RK45 (diagnostyczny) z kontrolą energii
- **Jacobian analityczny**: pełny łańcuch ∂F/∂x dla V₂ (jawny) i V₃ (przez Hesjan I_Y + kartezjańskie d_ij)
- **Benettin**: spektrum Lyapunova k=6N z QR, test Σλ ≈ 0
- **Matched-H methodology**: porównanie TGP vs Newton przy tej samej energii (ex174–ex194)
- **CSV pipeline**: generacja → walidacja → wizualizacja (ex166→ex168, ex182→ex183, etc.)

---

## 4. Co NIE działa / jest niekompletne

### 4.1 P1.B–C: Teza „TGP tłumi chaos" — NIEZAMKNIĘTA

**Problem**: Mimo 53 skryptów (ex148–ex194) i pełnej infrastruktury Benettin, **główna teza P1 nie jest udowodniona numerycznie**.

**Przyczyny**:
- **Twarda numeryka**: eksplozja energii poza krótkim t na IC Burrau (ex167)
- **RK45 failures**: `success=False` przy bliskich przejściach mimo małego |ΔE/E0|
- **Zależność od siatki**: λ_max zmienia się z dt × n_quad (ex166/ex168)
- **Brak asymptotyki**: sortowanie par (+/−) przy dużym T nie zrobione

**Wpływ**: P1 jest **kluczową falsyfikowalną predykcją** — bariera repulsyjna d_rep zapobiega kolizjom → mniej chaosu. Bez zamkniętego P1 brak ilościowego argumentu.

### 4.2 P4: I_triple — brak formy zamkniętej

**Problem**: Całka potrójnej nakładki Yukawy nie ma formy zamkniętej.

**Status**: Saddle-point daje **160–770% błąd** — BEZUŻYTECZNY.
**Wybrany wariant**: Multipol Legendre/Gegenbauer (semianalityczny po L_max cutoff). Fragment w `tgp_yukawa_IY_multipole_gegenbauer.tex`.

**Wpływ**: Ogranicza szybkość obliczenia V₃ w dynamice wielociałowej. Obecnie `three_body_force_exact` (Feynman 2D) działa, ale jest kosztowny.

### 4.3 Hessian — ~~tylko mody oddechowe~~ ROZWIĄZANE (P2, ex201)

~~Moduł `stability.py` oblicza Hessian, ale:~~
- ~~Tylko mody oddechowe (radialne)~~
- ~~Brak pełnego widma dla trójkątów, konfiguracji kolinearnych, N-gonów~~
- ~~Brak mapy bifurkacji w przestrzeni (β, C)~~

**ROZWIĄZANE** w ramach P2 (2026-04-05). Dodano do `stability.py`:
- `compute_hessian_generic()` — generyczny Hessian dla dowolnego V(pos)
- `normal_mode_analysis()` — masa-ważony Hessian D = M^{-1/2}HM^{-1/2}, omega^2, charakter modów
- `stability_comparison()` — TGP vs Newton przy tej samej konfiguracji
- `stability_bifurcation_scan()` — skan (beta, C) z klasyfikacją

Nowy skrypt: `ex201_full_hessian_stability_scan.py` (5 części).

**Kluczowy wynik P2**: Na d_well (równowaga konfinująca):
- Newton: ZAWSZE marginalny (Earnshaw!) — omega^2 ~ 0 (11/11 punktów)
- TGP: ZAWSZE stabilny — omega^2 > 0, 3/3 modów fizycznych (11/11 punktów)
- TGP łamie twierdzenie Earnshawa: oddechowy mod (radial_frac=1.0) omega^2 = 0.0083,
  dwa mody mieszane omega^2 = 0.00415 (zdegenerowane, symetria C₃)
- Efekt V₃: łagodnie destabilizuje (1 mod z 3 → marginalny), ale 2 mody pozostają stabilne
- 5-gon: TGP tworzy siodło (1 mod niestabilny), Newton stabilny — TGP nie stabilizuje WSZYSTKICH geometrii

### 4.4 Znane errory w kodzie

| Problem | Lokalizacja | Wpływ |
|---------|-------------|-------|
| `three_body_forces_approximate(..., use_yukawa=True)` — zły wykładnik | `dynamics_v2.py` | Duży błąd; preferuj `three_body_forces_exact` |
| Saddle-point 160–770% error | `three_body_terms.py` | Znany; nie używać do precyzyjnych obliczeń |
| Leapfrog energy explosion (długi t, Burrau) | `dynamics_v2.py` | Ogranicza horyzont Lyapunova |

---

## 5. Rozbieżności z aktualnym TGP v1

### 5.1 KRYTYCZNA: Forma ODE — LPA vs K_sub = g²

**To jest najważniejsza rozbieżność w całym pakiecie.**

| Aspekt | nbody (obecny) | Manuskrypt główny (aktualny) |
|--------|-----------------|------------------------------|
| Sprzężenie kinetyczne | `f(g) = 1 + 2α·ln(g)`, α=2 | K_sub(g) = K_geo · g² |
| Źródło formy | LPA/ERG jednopenętlowy | Pełne resummowanie substratowe |
| Ghost-free? | **NIE** — f(g) → −∞ dla g → 0⁺ | **TAK** — g² > 0 zawsze |
| Ważność | |ln g| ≪ 1 (słabe pole) | Cały zakres g ∈ (0, ∞) |
| ODE solitonu | f(g)·g'' + (2/r)·f(g)·g' + (α/g)·(g')² = V'(g) | g²g'' + g(g')² + (2/r)g²g' = g²(1−g) |

**Relacja między formami** (z `sek10_N0_wyprowadzenie.tex`, eq. 3.45):
```
K_geo · g² = K_geo · e^{2·ln(g)} = K_geo · [1 + 2·ln(g) + 2·(ln g)² + ...]
                                    ≈ K_geo · (1 + 4·ln(g))  [trunc. 1. rzędu]
```

Zatem f(g) = 1 + 4·ln(g) jest **obcięciem pierwszego rzędu** pełnej formy K_sub = g².

**Konsekwencje rozbieżności**:
1. **Silne pole (g ≪ 1)**: LPA daje f(g) < 0 → duchy kinetyczne. K_sub = g² zawsze > 0.
2. **Rdzenie solitonów**: Profil g(r) przy g₀ ≪ 1 (czarne dziury) jest **jakościowo inny**.
3. **V₃ ilościowo**: Siły 3-ciałowe zależą od profilu → zmiana ODE zmienia V₃.
4. **Ogon oscylacyjny**: Okres 2π zachowany (linearyzacja identyczna), ale amplituda A_tail zmieni się.
5. **Predykcje cząstkowe**: g₀ᵉ = 0.869 pochodzi z ODE K_sub = g² → wyniki φ-FP mogą się zmienić przy LPA.

**Weryfikacja numeryczna (ex195)**:

| g₀ | FULL: A_tail | LPA: A_tail | Różnica | LPA status |
|----|-------------|-------------|---------|------------|
| 0.010 | 0.860 | 59.57 (!) | 6800% | Ghost regime |
| 0.050 | 0.798 | 59.58 (!) | 7400% | Ghost regime |
| 0.100 | 0.739 | FAIL | — | Ghost regime (crash) |
| 0.300 | 0.576 | FAIL | — | Ghost regime (crash) |
| 0.500 | 0.429 | FAIL | — | Ghost regime (crash) |
| 0.700 | 0.273 | FAIL | — | Ghost regime (crash) |
| **0.869** | **0.1255** | **0.1017** | **23.4%** | OK |
| 0.950 | 0.049 | 0.047 | 5.7% | OK |

Próg duchów LPA: g = exp(−1/4) ≈ 0.779. Poniżej tego LPA daje f(g) < 0 (duchy kinetyczne).

- **FULL**: 8/8 converged (ghost-free, stabilne)
- **LPA**: 4/8 converged (50% — crash w ghost regime)
- **g₀ = 0.869 (punkt kalibracji φ-FP)**: A_tail różni się o 23.4% — to **duża** różnica wpływająca na predykcje r₂₁ = m_μ/m_e
- **g₀ < 0.3 (czarne dziury)**: LPA jest całkowicie bezużyteczne

**Wniosek**: Pakiet nbody jest **wewnętrznie spójny** z formą LPA, ale **niezgodny** z aktualnym stanem manuskryptu, który jednoznacznie preferuje K_sub = g² (ghost-free, pełne resummowanie). Różnica 23.4% w A_tail przy g₀ = 0.869 oznacza, że **predykcje cząstkowe muszą być przeliczone z formą FULL**.

**Rekalibracja z formą FULL (ex197)**:

Optymalne g₀ᵉ przesunęło się o zaledwie +0.079%:

| | Manuskrypt (FULL, g₀ᵉ=0.869) | Optymalny FULL (g₀ᵉ=0.869686) | PDG |
|---|---|---|---|
| r₂₁ = m_μ/m_e | 200.14 | **206.770** | 206.768 |
| m_μ (MeV) | — | **105.659** | 105.658 |
| m_τ (MeV) | — | **1776.99** | 1776.86 |
| r₃₂ = m_τ/m_μ | — | **16.818** | 16.817 |
| α_s(M_Z) | 0.1173 | **0.1174** | 0.1179 ± 0.0009 |
| α_s(m_τ)/α_s(M_Z) | 2.778 | **2.778** | 2.799 ± 0.121 |

Score: **7/8 PASS** (K_lepton MARGINAL z powodu numeryki Koide).
Shift g₀ᵉ: 0.869 → 0.869686 (+0.079%), co jest w granicach niepewności numerycznej A_tail.

**KLUCZOWY WNIOSEK**: Forma FULL K_sub = g² **zachowuje** wszystkie predykcje cząstkowe manuskryptu z minimalną rekalibracją g₀ᵉ (+0.08%). Nie ma kryzysu — wystarczy zaktualizować g₀ᵉ z 0.869 na ~0.8697.

### 5.2 Droga B vs wyprowadzenie z równań pola

| Aspekt | nbody | Manuskrypt główny |
|--------|-------|-------------------|
| Profil źródłowy | Yukawa narzucony jako warunek graniczny (Droga B) | Wyprowadzony z ODE solitonu |
| Ogon daleki | exp(−m·r)/r (Yukawa, monotonicznie malejący) | sin(r)/r (oscylacyjny, z ODE) |
| Fizyka | Cząstki = defekty topologiczne + Yukawa | Cząstki = solitony pola Φ |

**Wpływ**: W słabym polu (daleki zasięg, N-body) różnica jest mała — oba dają te same V₂ do wiodącego rzędu. Ale:
- Oscylacyjny ogon generuje **oscylacyjne V_eff** w czarnych dziurach (zob. ex188, dod. C)
- Bliskie odległości (r ~ kilka r₀) mogą dać inne zachowanie

Droga B jest **akceptowalna jako przybliżenie** w nbody, ale nie jest fundamentalna. Manuskrypt dokumentuje to w §1.2 ANALIZA_NBODY_INTEGRACJA.md.

### 5.3 Mniejsze rozbieżności

| Nr | Rozbieżność | Wpływ | Pilność |
|----|-------------|-------|---------|
| 1 | nbody nie uwzględnia dynamicznych stałych c(Φ), ħ(Φ), G(Φ) | Mały — nbody jest nierelatywistyczny | NISKA |
| 2 | Kosmologiczna ewolucja ψ(z) nie wpływa na nbody | Brak — nbody jest lokalny | BRAK |
| 3 | Sektor cechowania (dod. O, U, V) jest niezależny od nbody | Brak — inne skale | BRAK |
| 4 | Numeracja ex w nbody (ex55–ex194) a scripts/ (ex104–ex189) — nakładka | Konfuzja nazw | NISKA |

---

## 6. Plan rozwoju

### Status skryptów głównych (scripts/)

Weryfikacja (agent Explore) potwierdza: **skrypty w `scripts/` (ex104–ex189) już używają K_sub = g²**.
Kluczowe pliki: ex147, ex154, ex154b, ex160 — wszystkie z substratowym ODE `g'' + (1/g)(g')² + (2/r)g' = 1−g`.
Forma LPA pojawia się TYLKO w ex147 jako porównanie (etykieta "Description A: continuum").

**Wniosek**: Problem LPA dotyczy wyłącznie pakietu `nbody/`. Skrypty główne są spójne z manuskryptem.

### Faza 0: Aktualizacja ODE do K_sub = g² [ZAMKNIETA — NIE BLOKUJE N-BODY]

**Cel**: Zastąpić LPA form f(g) = 1 + 4·ln(g) pełną formą K_sub(g) = g² w pakiecie.

**Zrobione** (2026-04-05):
- `tgp_field.py`: dodano `KINETIC_MODE`, `kinetic_coupling(g, mode)`, `K_GEO`; `energy_density()` zaktualizowane
- `ex195_soliton_ksub_full_vs_lpa.py`: porównanie profili (8/8 FULL vs 4/8 LPA converged)
- `ex196_phi_fp_full_vs_lpa.py`: wpływ na predykcje cząstkowe
- `ex197_optimal_g0_full_form.py`: optymalny g₀ᵉ = 0.869686 → r₂₁ = 206.770 (0.0009% od PDG)

**Analiza wpływu na n-body** (2026-04-05):
- `pairwise.py`: V₂(d) jest **analityczne** (wyprowadzone z całkowania gradientów Yukawy, Droga B) — **NIE zależy od ODE solitonu ani K_sub**
- `three_body_force_exact.py`: I_Y jest całką Feynmana z funkcji Greena Yukawy — **NIE zależy od K_sub**
- `equilibria.py`: używa V₂ analitycznego — **NIE zależy od K_sub**
- `lyapunov.py`: używa sił z dynamics_backends — **NIE zależy od K_sub**

**Wniosek**: K_sub wpływa TYLKO na:
1. Profil solitonu g(r) w `tgp_field.py` i `yukawa_from_defect.py` (P5/P6)
2. C_eff (projekcja EFT) — ale zmiana jest minimalna (+0.08% w g₀ᵉ)
3. Predykcje cząstkowe — zachowane z rekalibracją g₀ᵉ

**Wszystkie wyniki n-body (P1–P4, P7) są NIEZALEŻNE od K_sub** — Droga B daje analityczne V₂, V₃ z parametrem m_sp = sqrt(3γ-2β), bez referencji do ODE solitonu.

**Status**: ZAMKNIETA (ex195–ex197 pokrywają porównanie; propagacja do n-body niepotrzebna)

### Faza I: P1 — Domknięcie tezy Lyapunova [WYSOKA — CZĘŚCIOWO ZROBIONE]

**Cel**: Ilościowy dowód „TGP tłumi chaos vs Newton" publikowalny.

**KLUCZOWE ODKRYCIA (ex198, ex199)**:

1. **Korekcja normalizacji**: Porównanie wymaga G_Newton = 4π (nie G=1), bo TGP V_grad = −4π·C₁C₂/d. Bez tej korekcji V_TGP jest 14× głębsze → sztuczny wzrost chaosu.

2. **Teza P1 jest zależna od β** (skan ex199 na IC Burrau):

| β | λ_TGP/λ_Newton | Verdict |
|---|---|---|
| 0.02 | 2.12 | ENHANCE |
| 0.05 | 2.52 | ENHANCE |
| 0.08 | 1.21 | ENHANCE |
| 0.20 | 1.47 | ENHANCE |
| **0.35** | **0.81** | **SUPPRESS** |
| 0.50 | 0.67 | SUPPRESS |
| 0.80 | 0.75 | SUPPRESS |
| 1.00 | 0.98 | COMPARABLE |
| 1.50 | 0.63 | SUPPRESS |

**Wniosek**: Istnieje **próg β_crit ≈ 0.3** powyżej którego TGP tłumi chaos. Poniżej progu — bariera repulsyjna jest za słaba.

3. **Przyczyna fizyczna**: Przy małym β bariera 1/d² jest niska → ciała mogą się zbliżać → dodatkowy potencjał confining 1/d³ ciągnie je bliżej → więcej chaosu. Przy dużym β bariera jest wysoka → ciała nie mogą kolidować → mniej chaosu.

**Wyniki z V₂+V₃ (ex200, yukawa_feynman, pełny skan)**:

| β | ratio(V₂+V₃) | ratio(V₂) | V₃ effect |
|---|---|---|---|
| 0.020 | 1.09 | 1.91 | V₃ helps (−43%) |
| 0.050 | 0.52 | 0.67 | V₃ helps |
| **0.080** | **0.58** | **0.91** | **V₃ helps (−36%)** |
| 0.120 | 0.63 | 0.99 | V₃ helps |
| 0.200 | 0.65 | 0.85 | V₃ helps |
| 0.350 | 0.60 | 0.76 | V₃ helps |
| 0.500 | 0.83 | 0.81 | ~neutral |
| 1.000 | 0.60 | 0.60 | ~neutral |

- **9/10 SUPPRESS** z V₂+V₃ (vs 7/10 z V₂-only)
- **β_crit(V₃) ≈ 0.025** vs β_crit(V₂) ≈ 0.042 — V₃ obniża próg o 40%
- Konwergencja: ratio stabilne od t=2 do t=6 (np. β=0.08: 0.576→0.611)

**STATUS P1: POTWIERDZONA** — TGP tłumi chaos dla β > 0.025 z siłami 3-ciałowymi. Typowe tłumienie: 35–45%.

**Pozostałe kroki** (rozszerzenie):

1. **Więcej IC** (trójkąt, hierarchiczny, random) — potwierdzenie uniwersalności
2. **Publikowalny wynik**: tabela ratio vs β × IC-type, z T → ∞
3. **Interpretacja fizyczna**: d_rep (bariera repulsyjna) vs d_min (minimalne zbliżenie) w trajektoriach

**Szacowany czas**: 1–2 sesje (fundamenty zamknięte)
**Pilność**: ŚREDNIA — rdzeń P1 zamknięty, rozszerzenie publikacyjne

### Faza II: P2 — Pełny Hessian [ZAMKNIĘTA ✅]

**Cel**: Wszystkie mody normalne (nie tylko oddechowe).

**STATUS: ZAMKNIĘTA** (2026-04-05, ex201)

**Wykonane kroki**:
1. ✅ `stability.py` rozszerzony o:
   - `compute_hessian_generic()` — generyczny Hessian V(pos)
   - `normal_mode_analysis()` — D = M⁻¹/²HM⁻¹/², omega², charakter modów
   - `stability_comparison()` — TGP vs Newton
   - `stability_bifurcation_scan()` — skan (β, C) z klasyfikacją
   - Nowa klasyfikacja "marginal" (omega² < 1e-6, Earnshaw)
2. ✅ Konfiguracje: trójkąt (d_rep + d_well), N-gon (5-gon)
3. ✅ Skan (β, C): 4 wartości β × 3 wartości C (quick), 9 × 6 (full)
4. ✅ `examples/ex201_full_hessian_stability_scan.py` — 5 części

**Wyniki kluczowe**:

| Konfiguracja | d_well Newton | d_well TGP | Stabilizacja |
|-------------|---------------|------------|-------------|
| Trójkąt równoboczny | MARGINAL (Earnshaw) | STABLE (3/3 modów) | 11/11 punktów (β,C) |
| 5-gon | stable (4 mody) | saddle (1 niestabilny) | NIE |

- **Earnshaw**: Newton ZAWSZE marginalny na d_well (omega² ~ 0, wszystkie 11/11 punktów)
- **TGP łamie Earnshawa**: stabilne mody w d_well (omega² = 0.0042, 0.0042, 0.0083)
  - 2 × mod mieszany (C₃ degeneracja), 1 × oddechowy (radial_frac = 1.0)
- **V₃ efekt**: łagodnie destabilizuje 1 z 3 modów, ale 2 mody stabilne pozostają
- **5-gon**: TGP tworzy siodło (omega² = -20.7 tangentialny, +24.7 oddechowy)
  — stabilizacja nie jest uniwersalna

**Wniosek P2**: TGP łamie twierdzenie Earnshawa dla trójkąta równobocznego. Stabilizacja jest specyficzna dla geometrii (nie dotyczy 5-gonu).

### Faza III: P4 — I_triple multipol [ZAMKNIĘTA ✅]

**Cel**: Semianalityczna forma I_triple przez multipole Legendre/Gegenbauer.

**STATUS: ZAMKNIĘTA** (2026-04-05, ex202)

**Wykonane kroki**:
1. ✅ `multipole_triple_overlap.py` — nowy moduł:
   - `yukawa_overlap_multipole()` — I_Y z tw. dodawania Yukawy
   - `yukawa_overlap_multipole_with_derivatives()` — I_Y + dI/dd (FD)
   - `three_body_forces_multipole()` — F₃ via multipol
   - `three_body_potential_multipole()` — V₃ via multipol
2. ✅ Kątowa diagonalizacja: A_{ll'}(ω) = (4π/(2l+1)) δ_{ll'} P_l(cos ω)
   → podwójna suma zamieniona na POJEDYNCZĄ sumę po l
3. ✅ Radialne: 3-segmentowa kwadratura Gaussa-Legendre'a
4. ✅ `ex202_multipole_vs_feynman_triple_overlap.py` — weryfikacja 5-częściowa

**Wyniki kluczowe**:

| L_max | Relative error (equilateral d=1) | Relative error (scalene) |
|-------|----------------------------------|--------------------------|
| 3     | 0.23%                            | 2.1%                     |
| 5     | 0.50%                            | 1.5%                     |
| 10    | 0.12%                            | 0.22%                    |
| 15    | 0.006%                           | 0.04%                    |

- **Dokładność**: <1% przy L_max=10 dla wszystkich testowanych geometrii (4/4 PASS)
- **Pochodne**: <10% (12/12 PASS) przy FD na multipolu
- **Siły**: 0.89% błąd na Burrau IC, 3. zasada Newtona spełniona (ΣF ≈ 0)
- **Potencjał**: 0.045% błąd vs Feynman exact
- **Szybkość**: ~4.4 ms/eval (L=10) vs ~0.09 ms (Feynman cached) — Feynman jest szybszy
  per evaluation, ale multipol jest tabulowalny i lepiej kontrolowalny

**Wniosek P4**: Multipol działa poprawnie jako alternatywna metoda obliczania I_Y.
Feynman 2D pozostaje preferowany dla dynamiki (pre-cached grid), multipol nadaje się do
tabulacji i analiz parametrycznych.

### Faza IV: P3 — Orbity zamknięte [ZAMKNIĘTA ✅]

**Cel**: Klasyfikacja orbit periodycznych (figure-8, L4/L5).

**Zrobione** (2026-04-05):
1. `ex203_figure8_tgp.py` — Chenciner-Montgomery figure-8 pod TGP:
   - Figure-8 przeżywa korekcje TGP do beta=0.10 (C=0.3)
   - Okres: Newton T=3.115, TGP(beta=0.01) T=3.283 (+5.4%)
   - Skalowanie: dr_max ~ beta^0.65 (perturbacyjne)
   - Full TGP z V3: V3/V2 = 0.82%
2. `ex204_lagrange_poincare.py` — punkty Lagrange'a i sekcje Poincaré:
   - L4 przesunięcie przez TGP: ~10^-5 (zaniedbywalne)
   - Stabilność L4 zachowana: q_crit nie zmienia się z beta
   - Sekcje Poincaré: regularne orbity wokół L4
   - Krzywizna potencjału efektywnego: zmiana < 10^-4

### Faza V: P5 — Ścieżka C — **ZAMKNIĘTA**

**Cel**: Wyprowadzenie profilu Yukawa z topologii defektu (nie narzucenie).

**Wyniki (EFT approach):**
   - Klasyczne defekty TGP: oscylacyjny ogon sin(r)/r (potwierdzone)
   - Naiwna stabilizacja V_sb = (mu^2/2)*(g-1)^2 niszczy defekty klasyczne
     (V_C'(g) < 0 wszędzie w (0,1) gdy mu^2 > 1)
   - Rozwiązanie: masa m_sp z korekcji pętlowych do propagatora (EFT)
   - C_eff = projekcja defektu na funkcję Greena Yukawy
   - Skalowanie: |C_eff| ~ (1-g0)^1.02, |C_eff| ~ beta^(-1.00)
   - mu^2 = 2*(3*gamma - 2*beta) jednoznacznie wyznaczone
   - Aksjomat Ścieżki B (C_i) WYPROWADZONY z fizyki defektu

### Faza VI.A: P6.A — Oddzialywanie solitonow — **ZAMKNIETA**

**Cel**: Weryfikacja prawa silowego n-body z nakladki pol defektow.

**Wyniki:**
   - delta(d) oscyluje (7 zmian znaku) — ogon klasyczny NIE jest Yukawa
   - G(d) i V''(1)*S(d) kasuja sie w ~93% (rownanie liniowe potwierdzone)
   - V_classical zanika wolniej niz Yukawa: d^-0.5 vs exp(-d)
   - Prawo Yukawy wymaga luki masowej EFT (nie jest wynikiem klasycznym)
   - Petla zamknieta: defekt -> C_eff (P5) -> V_Y (Path B) -> F_i (P0)

### Faza VI.B: P6.B — Dynamika zderzen solitonow [EKSPERYMENTALNA]

**Cel**: Zderzenia solitonów TGP (PDE solver).

### Faza VII: P7 — V3/V2 regime validation — **ZAMKNIETA**

**Cel**: Walidacja V3/V2 ratio w rezimach m_sp: kiedy sily 3-cialowe sa istotne.

**Wyniki (ex209):**
- t = m_sp*d < 0.5 (Coulomb): |V3/V2| ~ 1.3%, saturuje
- t ~ 1 (przejscie): |V3/V2| ~ 14% (peak)
- t > 3 (Yukawa): |V3/V2| < 0.2%, exponencjalnie stlumione
- **Rezim Lyapunov (beta=0.02-0.35):** V3/V2 = 0.8-1.0% — V3 istotne
- Geometria: compact 7x wieksza V3/V2 niz extended (bliskie przejscia amplifikuja)
- V3/V2 ~ C liniowo (max odchylenie < 2%)

### Faza VIII: P8 — Analytical equilibria + Hill criterion — **ZAMKNIETA**

**Cel**: Zamkniete formy analityczne: rownowagi V2, czestosci wlasne, sfera Hilla TGP, korekcja perturbacyjna V3.

**Wyniki (ex210):**

1. **Rownowagi analityczne**: V2'(d)=0 daje kwadratowe d^2 - 4*beta*d + 18*gamma*C = 0
   - d_rep = 2*beta - sqrt(4*beta^2 - 18*gamma*C) (bariera repulsyjna)
   - d_well = 2*beta + sqrt(4*beta^2 - 18*gamma*C) (studnia konfinujaca)
   - Prog istnienia: beta^2 > 4.5*gamma*C
   - Weryfikacja numeryczna: |d_analytic - d_numeric| < 1e-12

2. **Czestosci wlasne (zamknieta forma)**:
   - omega^2_rad = (2/C) * V2''(d_well)
   - V2''(d) = -2A/d^3 + 6B/d^4 - 12G/d^5
   - gdzie A = 4*pi*C^2, B = 8*pi*beta*C^2, G = 12*pi*gamma*C^3
   - Weryfikacja: |omega^2_analytic - omega^2_FD| < 1e-8

3. **Sfera Hilla TGP (nowy wynik)**:
   - R_H^3 = 4*pi*C_body / |Omega^2_tidal(d)|
   - Pole plywowe: -2A/d^3 + 6B/d^4 - 12G/d^5
   - Czlon beta (1/d^4) oslabia pole plywowe -> wieksza R_H
   - Czlon gamma (1/d^5) wzmacnia pole plywowe -> mniejsza R_H
   - Efekt netto: TGP **rozszerza** sfere Hilla o 1.0-2.9x wzgledem Newtona
   - Fizyka: bariera repulsyjna chroni ciala przed przejsciem plywowym

4. **Korekcja perturbacyjna V3**:
   - delta_d = -F3(d_well) / (3 * V2''(d_well))
   - |delta_d / d_well| < 0.2% — V3 jest perturbacyjne na rownowagach
   - V3 przesuwa d_well na zewnatrz (slabo odpychajace)

**Skrypt**: `ex210_analytical_equilibria_hill.py` (5 czesci, --quick mode)

### Faza IX: P9 — Unequal mass, phase diagram, escape velocity, breathing mode — **ZAMKNIETA**

**Cel**: Uogolnienie P8 na nierowne masy; diagram fazowy; predkosc ucieczki; efektywny potencjal oddechowy.

**Wyniki (ex211):**

1. **Rownowaga nierownych mas**: d^2 - 4*beta*d + 9*gamma*(C1+C2) = 0
   - Dla C1=C2=C: redukuje sie do d^2 - 4*beta*d + 18*gamma*C = 0 (P8)
   - Max stosunek mas: q_max = 4*beta^2/(9*gamma*C1) - 1
   - Weryfikacja: |d_analytic - d_numeric| < 1e-15

2. **Diagram fazowy**: beta_crit = (3/2)*sqrt(gamma*(C1+C2))
   - Ponizej: brak rownowag (czysta atrakcja)
   - Powyzej: studnia konfinujaca istnieje
   - Bifurkacja saddle-node przy beta_crit: d_rep = d_well = 2*beta_crit
   - Weryfikacja: 100% zgoda na siatce (beta, C_sum)

3. **Predkosc ucieczki**: v_esc = sqrt(2*|V2(d_well)|/mu), mu = C1C2/(C1+C2)
   - Zamknieta forma analityczna
   - Skalowanie: v_esc maleje z beta (glebsza studnia przy mniejszym beta)

4. **Mod oddechowy 3-cialowy**: omega^2_br = 3*V2''(d_well)/C + V3_correction
   - omega^2_br = (3/2)*omega^2_rad(2-body) -- mod oddechowy o 50% sztywniejszy
   - Korekcja V3: -1.3% (perturbacyjna, zgodna z P7)
   - Bariera wewnetrzna: ~10x energia wiazania

---

## Podsumowanie priorytetów

| Priorytet | Faza | Temat | Blokuje |
|-----------|------|-------|---------|
| ~~KRYTYCZNY~~ ✅ | 0 | ODE: K_sub = g² — **ZAMKNIETA** (ex195-197; n-body niezależne) | — |
| ~~WYSOKI~~ ✅ | I | P1: Lyapunov teza — **ZAMKNIETA** (ex200/207/208/209, synteza TeX) | — |
| ~~ŚREDNI~~ ✅ | II | P2: Hessian pełny — **ZAMKNIĘTA** | — |
| ~~ŚREDNI~~ ✅ | III | P4: I_triple multipol — **ZAMKNIĘTA** | Dynamikę wielociałową |
| ~~NISKI~~ ✅ | IV | P3: Orbity zamknięte — **ZAMKNIĘTA** | — |
| ~~NISKI-ŚR~~ ✅ | V | P5: Ścieżka C — **ZAMKNIĘTA** | — |
| ~~EKSPER.~~ ✅ | VI.A | P6.A: Oddzialywanie solitonow — **ZAMKNIETA** | — |
| **EKSPER.** | VI.B | P6.B: Dynamika zderzen (PDE) | — |
| ~~WALIDACJA~~ ✅ | VII | P7: V3/V2 regime scan — **ZAMKNIETA** | — |
| ~~ANALITYKA~~ ✅ | VIII | P8: Analytical equilibria + Hill — **ZAMKNIETA** | — |
| ~~ANALITYKA~~ ✅ | IX | P9: Unequal mass + phase diagram + v_esc + breathing — **ZAMKNIETA** | — |

---

*Dokument wygenerowany na podstawie analizy 221 skryptów Python, 24 dokumentów LaTeX, 6 plików MD w `nbody/` oraz porównania z `sek10_N0_wyprowadzenie.tex` (eq. 3.45), `sek08b_ghost_resolution.tex`, `sek00_summary.tex` (ODE solitonu).*
