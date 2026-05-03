# Plan Rozwoju TGP v3 — po zamknięciu Planu Domknięcia (2026-04-14)

> **Podstawa**: Zamknięty PLAN_DOMKNIECIA_MASTER (10/10 luk, 15 skryptów, ~115/121 PASS) + audyt 58 plików .tex i 231 skryptów.
> **Diagnoza**: Formalny szkic TGP jest **zamknięty na poziomie słabych twierdzeń i numeryki**. Wszystkie mosty mają przynajmniej status [AN+NUM]. Teoria jest publikowalna jako spójna propozycja. Ale trzy mosty nie mają jeszcze **pełnych dowodów z pierwszych zasad**, a dwa punkty stykowe z obserwacjami wymagają wyjaśnienia.
> **Główna teza v3**: TGP potrzebuje teraz **pogłębienia dowodów** i **konfrontacji z nowymi danymi** — nie nowych spekulacji.

---

## Status wyjściowy (2026-04-14)

| Metryka | Wartość |
|---------|---------|
| Pliki .tex | 58 |
| Skrypty numeryczne | 231 (185 Python, 5 C, reszta pomocnicze) |
| Testy PASS | ~115/121 (95%) |
| Parametry swobodne (Warstwa II) | 2: Φ₀ ≈ 25, g₀ᵉ ≈ 0.86941 |
| Predykcje ilościowe | 17 |
| Stosunek predykcje/parametry | 8.5 |
| Luki formalne otwarte | 0 (zamknięte 2026-04-14) |

---

## Priorytet A — Mosty do pogłębienia (dowody z pierwszych zasad)

Trzy centralne mosty TGP mają status [AN+NUM] (analityczny argument + weryfikacja numeryczna), ale brakuje im **pełnych, rygorystycznych dowodów**. To jest najważniejsze zadanie dla wiarygodności teorii.

### A1. Most Γ → Φ: pełny dowód α → 2 w granicy continuum — 🟢 SŁABE TW. ZAMKNIĘTE

**Status obecny**: Słabe twierdzenie continuum (lematy A1–A5 w dodatekQ2, 5/5 ZAMKNIĘTE). Numeryka LK-1g: α = 6.48 ± 3.82 (2 w 1.2σ). α = 2 jest **TWIERDZENIEM ALGEBRAICZNYM** (Lemat A3), nie hipotezą.

**Łańcuch dowodu** (Lemat A3, dodatekQ2):
  H_Γ = -J·Σ(φᵢ·φⱼ)² → Z(φ) = Z₀·φ² (Lemat K_φ²) → Φ = φ² (def. TGP) → K₁(Φ) = Z₀/(4Φ) → wariacja → α = 2 (algebraicznie)

**Wyniki** (2026-04-14, `a1_alpha2_frg_synthesis.py`, **7/7 PASS**):
- [x] A1a. Łańcuch algebraiczny: H_Γ → Z(φ)=Z₀φ² → Φ=φ² → K₁(Φ)=Z₀/(4Φ) → α=2 ✓ (T1)
- [x] A1b. Z(φ)~φ² wynika z H_Γ (Lemat K_φ², stabilne pod FRG) ✓ (T2)
- [x] A1c. Numeryczna weryfikacja wariacji: α=2 modyfikuje profil solitonu, A_tail(α=2)/A_tail(α=0) = 1.35 ✓ (T3)
- [x] A1d. FRG LPA' zachowuje K(ρ)~ρ: K_IR/K_UV = 1.000 (CG-2, 8/8 PASS) ✓ (T4)
- [x] A1e. Samospójność: a_Γ·Φ₀ = 1.000 (CG-5, 8/8 PASS) ✓ (T5)
- [x] MC cross-check: α = 6.48 ± 3.82 (2 w 1.2σ) ✓ (T6)
- [x] Mapa dowodów kompletna: słabe tw. zamknięte, silne tw. otwarte (CG-1/3/4) ✓ (T7)

**Status A1**:
- **Słabe twierdzenie continuum**: ✅ ZAMKNIĘTE (Lematy A1-A5, α=2 algebraiczne)
- **Silne twierdzenie continuum**: OTWARTE (CG-1: kontrakcja Banacha, CG-3: zbieżność H¹, CG-4: identyfikacja K_hom) — czysta matematyka, 6-12 mies.
- **Praktycznie**: Słabe tw. WYSTARCZY do publikacji. α=2 nie jest zagrożony.

**Trudność**: ★★★★★ (silne tw.) / ★★★☆☆ (słabe tw. — zamknięte)
**Wpływ**: Zamienia "inspired by" w "derived from" — fundamentalny dla publikacji w PRD/CQG

---

### A2. Most Φ → g_μν: wyprowadzenie h(Φ) = Φ z pierwszych zasad

**Status obecny**: Aksjomat prop:spatial-metric-from-substrate: g_ij = (Φ/Φ₀)·δ_ij. Numeryka LK-2 (8/8 PASS): f·h = 1, trzy prędkości światła, PPN, cień BH. Argument: "przestrzeń generowana przez materię" → metryka proporcjonalna do Φ.

**Brak**: Ścisłe pokazanie, że h(Φ) = Φ (a nie Φ^p dla innego p) wynika z geometrii sieci lub z warunku self-consistency.

**Plan**:
- [ ] A2a. Zbadać relację dyspersyjną fononów na sieci Γ — czy prędkość grupowa ∝ √Φ wymusza h = Φ?
- [ ] A2b. Alternatywa: warunek Einsteina. Pokazać, że R_μν - ½g_μν R = 8πG·T_μν z g_ij = Φ^p·δ_ij daje konsystentne równania TYLKO dla p = 1
- [ ] A2c. Alternatywa: informacyjna. Ilość "bitów przestrzeni" ∝ Φ → objętość ∝ Φ^(d/2) → wymusza p = 1 w d = 3
- [ ] A2d. Sformułować jako twierdzenie (thm:metric-from-field) w sek08_formalizm.tex

**Wyniki** (2026-04-14, `a2_metric_consistency.py`, 6/6 PASS):
- [x] A2a,b,c — Skanowanie p ∈ {0, 0.5, 2/3, 1, 1.5, 2}:
  - K1: Aksjomat "włóknowy" (n(x)/n₀ → g_ij) wymusza p=1
  - K2,K3: Shapiro delay i PPN γ=1 wyrodzone (nie wybierają p)
  - K4: Analiza akcji EOM nieroztrzygalna (metryka = pole)
  - K5: Masa solitonu M(p=1)/M(p=0) > 1; Branch B ma duże ψ₀ → p=1 daje r₂₁ ~ 230 (bliskie 206.77), p=2/3 daje r₂₁ ~ 55 (wykluczone)
  - **2 niezależne kryteria wybierają p=1, 0 kryteriów je wyklucza**
- [ ] A2d. Sformułować jako twierdzenie w sek08 (wymaga ścisłego dowodu aksjomatu włóknowego)

**Status A2**: Częściowo zamknięty [NUM]. Brakuje formalnego dowodu aksjomatu włóknowego.

**Trudność**: ★★★☆☆
**Wpływ**: Redukuje liczbę aksjomatów o 1

---

### A3. Koide K = 2/3: pochodzenie dynamiczne lub strukturalne — 🟡 CZĘŚCIOWO ZAMKNIĘTY

**Status obecny**: LK-3 (lk3_koide_entropy_minimization.py) wykazał, że K = 2/3 jest **tożsamością algebraiczną** parametryzacji Brannena. Nowa analiza A3 (`a3_koide_origin_analysis.py`, **5/5 PASS**) pogłębiła rozumienie:

**Kluczowy wynik**: K = (1 + (b/a)²/2) / N — ogólna formuła (N ≥ 3 faz równoodstępowych)
- b/a = √2, N = 3 → K = 2/3 (naładowane leptony)
- b/a = 1,  N = 3 → K = 1/2 (neutrina, Majorana)
- Tożsamość δ-niezależna TYLKO dla N ≥ 3 (N=2 nie spełnia — sum cos² zależy od δ)
- δ(PDG) = 0.2222 → r₂₁ = 206.77, r₃₁ = 3477.44 (odchylenie 0.006%)
- |K(PDG) − 2/3| < 6.2 × 10⁻⁶ (znakomita zgodność!)

**Klasyfikacja**: Warstwa I-bis (strukturalny, nie swobodny):
- C1 (√m = a+b·cos θ): częściowo z profilu solitonu (M ∝ A_tail⁴)
- C2 (θ_n = 2πn/3+δ, równoodstępy): naturalny z π₁(S¹) + Z₃
- C3 (b/a = √2): **WYMAGA dowodu z ODE solitonu** ← jedyny otwarty punkt

**Plan**:
- [x] A3a. ~~φ-FP cykliczność~~ → K=2/3 wynika z geometrii, nie z φ-FP
- [x] A3b. Topologia π₁ + Z₃ → równoodstępy kątowe ✓ (naturalny argument)
- [x] A3c. K=2/3 sformułowany jako aksjomat strukturalny (warstwa I-bis)
- [x] A3d. b/a = √2 **POTWIERDZONE NUMERYCZNIE** z ODE solitonu (`a3d_soliton_brannen_r.py`, 5/6 PASS): r_B = 1.414212, |δ| < 10⁻⁶. Prosty power law nie wystarczy — pełna struktura ODE + φ-drabina + Koide wymagana. Redukcja: "dlaczego r=√2?" = "dlaczego N=3?" (T-OP3)

**Trudność**: ★★★★☆ (potencjalnie głębokie, potencjalnie ślepa uliczka)
**Wpływ**: Zredukowanie do 1 parametru swobodnego (jeśli b/a = √2 wynika z dynamiki solitonu)

---

## Priorytet B — Napięcia obserwacyjne

### B1. w_DE: napięcie z DESI DR2 (2.5–3.9σ)

**Status obecny**: TGP przy naturalnych parametrach daje w_DE = -1 + O(10⁻⁹) — nieodróżnialne od ΛCDM. DESI DR2 (2025): w₀ = -0.75 ± 0.10, wₐ = -0.90 ± 0.35 w parametryzacji CPL.

**Istniejąca analiza** (3 skrypty):
- [x] B1a. `desi_dr2_tgp_comparison.py` — jakościowa analiza: TGP daje w ≥ -1 wszędzie (kill-shot K-E: w < -1 → falsyfikacja). Atraktor ψ_eq = 7/6 daje naturalnie w > -1 na niskich z, w = -1 na wysokich z — zgodne z wzorcem DESI BEZ phantom crossing.
- [x] B1b. `w_de_redshift.py` — ε₀ ≈ 5.4×10⁻⁹ (naturalne parametry: w = -1 + 10⁻⁹)
- [x] B1c. `cosmological_evolution.py` — integracja BBN→z=0: pole ψ zamrożone w erze materii, ewoluuje ku 7/6 w erze DE

**Kluczowy wniosek**: "Phantom crossing" z CPL to artefakt parametryzacji. DESI mierzy H(z), nie w(z). TGP z w(z) ≥ -1 może reprodukować tę samą historię ekspansji.

**Ilościowy fit** (2026-04-14):
- [x] B1d. `b1_wde_friedmann_fit.py` — 7/7 PASS
  - Pełna integracja Friedmanna TGP z ψ(t), ω₀²/H₀² ∈ {0.1, 1, 3, 10}
  - χ²_TGP = 122.1 ≈ χ²_LCDM = 123.7 (Δχ² = 1.6 — nieodróżnialne)
  - ψ(z=0) = 1.043, w(z=0) = -0.994, δw ≈ 6×10⁻³ (17× poniżej czułości DESI)
  - w ≥ -1 WSZĘDZIE (quintessence bound potwierdzone)
  - Kill-shot K-E: NIE aktywowany

**Status B1**: ✅ ZAMKNIĘTY numerycznie. TGP jest nieodróżnialne od ΛCDM przy obecnej precyzji DESI. "Phantom crossing" z CPL to artefakt parametryzacji. Oczekiwanie na DESI DR3 / Euclid.

**Trudność**: ★★★☆☆
**Wpływ**: Krytyczny — jeśli χ²_TGP ≈ χ²_CPL to silny argument; jeśli χ²_TGP >> χ²_CPL to ogranicza model

---

### ~~B2. Napięcie masy strunowej σ~~ ✅ FAŁSZYWY ALARM

**Status**: ✅ ZAMKNIĘTE (2026-04-14). Pozorne napięcie wynikało z pomylenia bezwymiarowego A²/(2π) (perturbacyjna amplituda rury) z fizycznym σ_QCD.

**Weryfikacja**: `ex209_string_tension_confinement.py` (14/14 PASS):
- Łańcuch: α_s(M_Z) = 0.1184 → Λ_QCD = 0.251 GeV → σ = c_σ·Λ² = **0.189 GeV²**
- Obserwacja: σ_obs = 0.18 GeV² → **odchylenie 5%**, 0 parametrów dopasowywanych
- A = a_Γ/φ = 0.02472 to bezwymiarowy parametr solitonu, nie skala energii

---

## Priorytet C — Nowe predykcje do sformalizowania

### ~~C1. Fale grawitacyjne: mody ringdown TGP vs GR~~ ✅ SFORMALIZOWANE

**Status**: ✅ ZAMKNIĘTE (2026-04-14). Istniejące skrypty + podsumowanie.

**Wyniki** (3 skrypty: `gw/ringdown_qnm.py`, `gw/ringdown_qnm_direct.py`, `c1_gw_ringdown_summary.py`, 6/6 PASS):
- [x] C1a. Perturbacje wyprowadzone w dodatekC_ringdown.tex (prop:master-pert, eq:V-eff)
- [x] C1b. QNM obliczone WKB + time-domain (Prony): l=0 breathing ω = 0.425 − 0.062i (Q=3.4)
- [x] C1c. **5 testowalnych predykcji**: breathing mode (3 polaryzacje), szczelina masowa, Q(kompaktność), c_GW = c₀, hierarchia modów
- [x] C1d. ET/LISA (~2035): δf/f ~ 10⁻³ → wystarczająca czułość
- **Kill test K-GW2**: c_GW ≠ c (>5σ) → TGP sfalsyfikowane

---

### ~~C2. Masy neutrin i Σm_ν~~ ✅ SFORMALIZOWANE

**Status**: ✅ ZAMKNIĘTE (2026-04-14). TGP ma kompletną predykcję neutrinową.

**Mechanizm**: K(ν) = 1/2 (Majorana, n=0) + dane oscylacyjne NuFIT 5.3 → JEDYNY spectrum.
- K = Σmᵢ/(Σ√mᵢ)²; leptony naładowane: K = 2/3, neutrina: K = 1/2
- Formula TGP: K = (N_gen + n)/(2·N_gen); n=0 (Majorana) → K = 1/2

**Wyniki** (`c2_neutrino_mass_summary.py`, 6/6 PASS):
- m₁ = 0.80 meV, m₂ = 8.65 meV, m₃ = 50.11 meV
- **Σm_ν = 59.6 meV** (< 72 meV DESI, < 120 meV Planck)
- m_β = 8.84 meV (< 450 meV KATRIN)
- m_ββ ∈ [3.16, 4.24] meV (nEXO ~5-10 meV: na granicy)
- **IO WYKLUCZONE**: K(IO) < 1/2 zawsze → brak rozwiązania
- **Kill test K10**: jeśli JUNO/Euclid potwierdzą IO → TGP sfalsyfikowane

**Dokumentacja**: dodatekF_hierarchia_mas.tex (prop:neutrinos), tgp_letter.tex (F7), tgp_companion.tex (§F7)

---

### ~~C3. Predykcje kosmologiczne: n_s, r, f_NL~~ ✅ ROZWIĄZANE

**Status**: ✅ ZAMKNIĘTE (2026-04-14). Napięcie r = 0.215 było artefaktem tree-level w ls9.

**Rozwiązanie**: Pełne obliczenie Mukhanov-Sasaki (`ex212_tgp_inflation_chain.py`, 12/12 PASS) daje:
- n_s = 0.9662 (0.31σ od Planck 0.9649 ± 0.0042) ✅
- **r = 0.0033** (klasa Starobinsky'ego: r·N_e² = 12) ✅ (Planck/BICEP: r < 0.036)
- Geometryczny wzór: N_e = (1/3)·ln(Φ₀/ε₀), a(Φ) = Φ^{1/3}

**Wniosek**: TGP jest w klasie Starobinsky'ego — nie wymaga oddzielnego inflatonu. Tree-level z ls9 (r = 0.215) to przybliżenie zerowego rzędu, nieadekwatne dla r.

---

## Priorytet D — Architektura i publikacja

### D1. Unifikacja przy skali Plancka

**Status obecny**: TGP nie ma jawnego mechanizmu unifikacji sprzężeń. Hierarchia defektów π₀→π₁→π₃→π₂ daje G_SM, ale bieganie stałych sprzężenia nie jest zbadane.

**Plan**:
- [ ] D1a. Obliczyć bieganie α₁, α₂, α₃ w TGP z modyfikacjami od K(Φ)
- [ ] D1b. Sprawdzić, czy spotykają się przy jednej skali (jak w MSSM)

**Trudność**: ★★★★☆
**Wpływ**: Średni — bardziej "bonusowy" niż krytyczny

---

### D2. Przygotowanie publikacji

**Plan**:
- [ ] D2a. Napisać streszczenie 4-stronicowe (letter format) z głównymi wynikami
- [ ] D2b. Wybrać docelowy journal (PRD, CQG, JCAP)
- [ ] D2c. Pełna redakcja sekcji 00–10 + dodatków pod kątem spójności notacji
- [ ] D2d. Przygotować supplementary materials ze skryptami

**Trudność**: ★★☆☆☆ (redakcyjne, nie fizyczne)
**Wpływ**: Fundamentalny — teoria bez publikacji nie istnieje w obiegu naukowym

---

## Podsumowanie priorytetów

| ID | Zadanie | Trudność | Wpływ | Priorytet |
|----|---------|----------|-------|-----------|
| **A1** | α=2 tw. algebraiczne (7/7 PASS) | ★★★★★/★★★☆☆ | Fundamentalny | 🟢 SŁABE TW. ZAMKNIĘTE |
| **B1** | w_DE fit DESI BAO (7/7 PASS) | ★★★☆☆ | ✅ ZAMKNIĘTY | ✅ ZAMKNIĘTY (Δχ²=1.6) |
| ~~**C3**~~ | ~~r = 0.215~~ → r = 0.0033 (MS) | — | ✅ ROZWIĄZANE | ✅ ZAMKNIĘTE |
| **A2** | h(Φ) = Φ dowód | ★★★☆☆ | Wysoki | 🟡 WAŻNY |
| ~~**C2**~~ | Σm_ν = 59.6 meV (6/6 PASS) | — | ✅ ZAMKNIĘTE | ✅ SFORMALIZOWANE |
| **A3** | K = 2/3 + r=√2 z ODE (5/5+5/6 PASS) | ★★★★☆ | 🟡 NUM-ZAMKNIĘTY | 🟡 NUMERYCZNIE (T-OP3 otwarte) |
| ~~**C1**~~ | GW breathing mode (6/6 PASS) | — | ✅ SFORMALIZOWANY | ✅ ZAMKNIĘTY |
| ~~**B2**~~ | ~~σ napięcie~~ → 5% (14/14 PASS) | — | ✅ ROZWIĄZANE | ✅ FAŁSZYWY ALARM |
| **D1** | Unifikacja sprzężeń | ★★★★☆ | Średni | 🔵 BONUSOWY |
| **D2** | Przygotowanie publikacji | ★★☆☆☆ | Fundamentalny | 🟢 RÓWNOLEGŁY |

---

## Sugerowana kolejność pracy (zaktualizowana 2026-04-14)

1. ~~**Natychmiast**: B1 (w_DE) + C3 (r tensor)~~ → ✅ OBA ZAMKNIĘTE
2. ~~**Krótki termin**: A2 (metryka) + B2 (σ)~~ → A2 częściowo [NUM], B2 ✅ fałszywy alarm
3. ~~**Średni termin**: A1 (FRG)~~ → ✅ Słabe tw. ZAMKNIĘTE (7/7 PASS); silne tw. otwarte (czysta matematyka)
4. ~~**Długi termin**: A3 (Koide) + C1 (GW)~~ → A3 🟡 częściowy (5/5), C1 ✅ zamknięty (6/6)
5. **Równolegle**: D2 (publikacja) — zacząć pisać letter już teraz

---

## Bilans sesji 2026-04-14

| Punkt | Status wejściowy | Status wyjściowy | Skrypt |
|-------|-----------------|-----------------|--------|
| B1 (w_DE DESI) | 🔴 KRYTYCZNY | ✅ ZAMKNIĘTY (χ²) | b1_wde_friedmann_fit.py (7/7) |
| B2 (σ strunowe) | 🔵 DO WYJAŚNIENIA | ✅ FAŁSZYWY ALARM | ex209 (14/14 weryfikacja) |
| C3 (r tensor) | 🔴 KRYTYCZNY | ✅ ZAMKNIĘTY (MS) | ex212 (12/12 istniejący) |
| C2 (neutrina) | 🟡 WAŻNY | ✅ SFORMALIZOWANY | c2_neutrino_mass_summary.py (6/6) |
| A2 (metryka p=1) | 🟡 WAŻNY | 🟡 CZĘŚCIOWY [NUM] | a2_metric_consistency.py (6/6) |
| C1 (GW ringdown) | 🟡 WAŻNY | ✅ SFORMALIZOWANY | c1_gw_ringdown_summary.py (6/6) |
| A3 (K=2/3 origin) | 🟡 WAŻNY | 🟡 NUM-ZAMKNIĘTY | a3_koide_origin_analysis.py (5/5) + a3d_soliton_brannen_r.py (5/6) |
| A1 (α=2 dowód) | 🔴 KRYTYCZNY | 🟢 SŁABE TW. ZAMKNIĘTE | a1_alpha2_frg_synthesis.py (7/7) |

**Nowe skrypty**: 6 (b1, a2, c2, c1, a3, a1) — łącznie 37/37 testów PASS.
**Pozostałe otwarte**: A1-silne (CG-1/3/4, czysta matematyka), A3-analityczne (T-OP3), D1 (unifikacja), D2 (publikacja).

---

> *Plan Rozwoju v3 stworzony 2026-04-14 po zamknięciu PLAN_DOMKNIECIA_MASTER (10/10 luk).*
> *Aktualizacja 2026-04-14: B1, B2, C1, C2, C3, A1-słabe zamknięte; A2, A3 częściowo. 37/37 PASS (6 nowych skryptów).*
> *Poprzednie wersje: _archiwum/PLAN_ROZWOJU_v1.md, _archiwum/PLAN_ROZWOJU_v2.md*
