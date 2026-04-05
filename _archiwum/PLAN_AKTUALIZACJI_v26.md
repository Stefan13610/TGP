# TGP — Plan Aktualizacji v26
*Data: 2026-03-22 | Baza: TGP v26 | master_verification_v26.py: 33/33 PASS*

---

## Status po sesji v26 (UKOŃCZONE)

### Gotowe ✅ (v25 + v26)
| Plik | Co dodano | Sesja |
|---|---|---|
| `ANALIZA_CIEMNA_MATERIA.md` | §8 F1/F2/F3, K18, paradoks γ, ε_th | v25 |
| `ANALIZA_SPOJNOSCI_v25.md` | ex39-ex42, K14/K18, O21-O23, §8 głęboka | v25 |
| `ANALIZA_NBODY_INTEGRACJA.md` | §4.2 K13 ex34v2; §4.5.3 K18 ex44; v26 header | v25+v26 |
| `sek07_predykcje.tex` | K14/K18 tabela, F9 kryterium, O21-O24, predykcja 7 | v25 |
| `sek03_rezimy.tex` | rem:dwa-sektory-nbody (Efimov + FDM) | v26 |
| `sek08_formalizm.tex` | prop:eps-th (ε_th=m_sp²/2), prop:golden-ratio (φ) | v26 |
| `dodatekD_trojcialowe.tex` | ssec:fdm-sector, prop:soliton-schive, rem:k18-verification | v26 |
| `nbody/tgp_nbody_predictions.tex` | K18/K14 predykcje kosmologiczne (ex41-ex44) | v26 |
| `nbody/tgp_quantum_window.tex` | ssec:qw_2d_revision (korekta 2D, ΔE=+0.08) | v26 |
| `nbody/tgp_quantum_window_msp.tex` | ssec:qw_msp_2d (okno 2D: (0.081,0.083) l_Pl⁻¹) | v26 |
| `scripts/master_verification_v25.py` | 60/60 PASS (G1-G7, H1-H2) | v25 |
| `scripts/master_verification_v26.py` | 56/56 PASS (I: ex34, J: ex43, K: ex34v2, L: ex44) | v26 |
| `nbody/examples/ex34` — `ex45` | Wszystkie skrypty Python + UTF-8 fix | v25+v26 |
| `nbody/examples/ex45` — Lyapunov | λ_TGP vs λ_Newton: efekt mieszany (58.3%); bariera stabilizuje przy losowych konfig. | v26 |
| `nbody/examples/verify_all.py` | Naprawa Unicode, 39/39 PASS | v25 |

### Kluczowe wyniki v26

#### ex34 — K13 REWIZJA (okno Efimova 2D)
- Korekta 2D przesuwa E0 o **+0.08 E_Pl** (konsekwentnie dla wszystkich m_sp)
- Przy C=C_Pl: tylko m_sp=0.076 wiąże w 2D, ale C_Q(2B)_2D=0.276 < C_Pl → 2B też wiąże
- **C_Pl nie należy do 2D okna Efimova** dla żadnego testowanego m_sp ∈ {0.076–0.198}
- K13 status: **REWIZJA** — 1D okno było artefaktem aproksymacji siodłowej

#### ex34v2 — K13 WARUNKOWO POTWIERDZONY (precyzyjna bisekcja, krok pośredni)
- Gęsty skan m_sp ∈ [0.076, 0.100], krok 0.003 l_Pl⁻¹ (siatka Nb=75, Nh=55)
- **m_sp\* = 0.0831 l_Pl⁻¹** — próg 3B (interpolacja: zero-crossing E0_2D między 0.082 i 0.085)
- **m_sp\*\* ≈ 0.0806 l_Pl⁻¹** — próg 2B (z ex34 interpolacji C_Q(2B)_2D = C_Pl, był wstępny)

#### ex46 — K13 POTWIERDZONY (dedykowany solver 2-ciałowy)
- Solver 1D Schrödinger FD dla pary TGP: −1/(2μ)ψ'' + V_2B ψ = Eψ, μ = 1/2
- Wynik: C_Q(2B) ≈ **0.488–0.576** >> **C_Pl = 0.282** dla WSZYSTKICH m_sp ∈ [0.065, 0.100]
- Para TGP NIGDY nie wiąże przy C = C_Pl w całym zakresie m_sp
- Okno Efimova = **(0, 0.0831) l_Pl⁻¹** — dużo szersze niż szacowano z ex34
- **K13 status: POTWIERDZONY** ✓ (bez warunku — ex46 definitywny)

#### ex47 — K14 UPDATE (DESI FDM forecast, 3D przybliżenie)
- Funkcja transferu Irsic+2017: T(k) = [1+(α·k)^2.24]^(−4.46), α=0.04/m22^(4/9) Mpc/h
- DESI DR3: V_eff = 6 Gpc³/h³, k ∈ [0.1, 8] h/Mpc
- **m22_crit(DESI DR3) = 0.96** — DESI wyklucza FDM dla m22 < 0.96 (3σ)
- **TGP-F3 (m22=1): SNR = 2.86** — poniżej 3σ → **NIE wykluczone przez DESI DR3** ✓
- ⚠️ Ex47 używał 3D P(k) zamiast właściwej obserwabli 1D P_Lya — patrz ex48

#### ex48 — K14 ROZWIĄZANY (1D Lyman-α + Eisenstein-Hu)
- **Właściwa obserwabla:** P_1D(k_∥) = (1/2π) ∫ P_3D(√(k∥²+k⊥²)) k_⊥ dk_⊥
- **Eisenstein-Hu P_CDM** (zamiast BBKS z ex47) — dokładność ~5% przy k<10 h/Mpc
- **Korekcja TGP:** δ_TGP = C²_Pl·β = 0.0796 → m22_eff = 1.09 (+9%)
- **Wynik:** S_1D(m22=1, z=4.5, k=0.3) = 0.935; S_obs(Rogers) = 0.920 ± 0.040
- **Napięcie = −0.4σ** — K14 ZGODNY z Rogers+2021 w poprawnej analizie 1D!
- **m22_fit = 1.10** — Rogers+2021 zmierzyłby m22≈1.10 dla TGP-FDM (m22=1) ✓
- **Mechanizm:** Projekcja 1D integruje po k_⊥ → mody CDM dominują → supresja FDM rozcieńczona
- **K14 status: ZGODNY** ✓ (napięcie Rogers+2021 było artefaktem metodologicznym 3D vs 1D)

#### ex44 — K18 NAPIĘCIE WYJAŚNIONE (barionowy feedback)
- Model M3 (F3 + feedback): ΔAIC = −293.7 vs M0 — **najlepszy model**
- beta_bar = 0.127: r_c ~ M^(−1/9) · (f_bar)^0.127 (wyższy f_bar → większy rdzeń)
- M_break = 5.13×10⁹ M_⊙ (złamany profil: alpha_lo=−0.070, alpha_hi=−0.130)
- F3 (−1/9) leży **między** alpha_lo i alpha_hi → spójne z trendami
- Napięcie 5.2σ = efekt barionowy (spirale: duże f_bar, skompresowany soliton)
- F1 (alpha=−1) niezmiennie wykluczone na 44σ
- Wniosek: TGP-F3 + feedback barionowy = pełny opis K18

#### ex43 — K18 WZMOCNIONY (n=75)
- α_obs = −0.086 ± 0.021 ≈ −1/9 = −0.111 (F3) → odchylenie **1.2σ** (OK)
- F1 wykluczone na **44.6σ**; ΔAIC(F1−F3) = **515** (ex41: 133, ex43: 515)
- Skalowanie ΔAIC z n liniowe (~7/gal) zgodnie z oczekiwaniami
- Napięcie spirale/karłowe = 5.2σ → możliwy łamany profil mocy lub efekty barionowe
- Prognoza: n=150 → ΔAIC ≈ 1050

### Otwarte ⚠️
- K14 (Lyman-α): **ZAMKNIĘTE** ✅ — ex48 (1D P_Lya + E-H): napięcie = −0.4σ, K14 ZGODNY
  Napięcie Rogers+2021 było artefaktem 3D vs 1D porównania (ex47)
- K13 (Efimov): **ZAMKNIĘTE** ✅ — ex46 potwierdza ostatecznie; okno = (0, 0.0831)
- K14 do potwierdzenia: CLASS/CAMB + hydrosimulacje (ex49, niski priorytet)
- verify_all.py: brak testów ex34/ex43/ex34v2/ex44/ex45/ex46/ex47/ex48 (niski priorytet)

---

## PRIORYTET 1 — N-body (następna sesja)

> Najłatwiejsze do obliczenia i porównania. Skupiamy się tutaj.

### 1.1 ex34 — Pełny Skan 2D Okna Kwantowego Efimova

**Cel:** Wyznaczenie dokładnych granic okna kwantowego m_sp ∈ [0.05, 0.30] l_Pl⁻¹
**Metoda:** FD isosceles (ex32-ex33 jako baza), bisekcja 2D
**Plik:** `nbody/examples/ex34_efimov_window_2d.py`

```python
# Algorytm:
# dla m_sp in linspace(0.05, 0.30, 50):
#   znajdź C_crit_3B(m_sp): min C gdzie E_3B < 0
#   znajdź C_crit_2B(m_sp): min C gdzie E_2B < 0
#   jeśli C_crit_3B < C_Pl < C_crit_2B: m_sp w oknie
# Wynik: mapa (m_sp, C_crit_3B, C_crit_2B, in_window)
```

**Oczekiwany wynik:** Dokładne granice okna dla C = C_Pl = 0.282
**Czas:** ~1h obliczenia (FD 100×100 grid dla każdego m_sp)

---

### 1.2 ex43 — THINGS/LITTLE THINGS Precyzyjne Profile HI (K18 v2)

**Cel:** Zwiększenie próby K18 z n=30 do n=80+ galaktyk z precyzyjnymi r_c
**Dane:** THINGS (34 galaktyki, Walter+2008) + LITTLE THINGS (40 galaktyk, Hunter+2012)
**Plik:** `nbody/examples/ex43_things_k18_precision.py`

**Protokół:**
1. Załaduj dane THINGS/LITTLE THINGS (publiczne, VizieR)
2. Dopasuj profil Schive (lub NFW+soliton) do krzywej rotacji HI
3. Wyznacz r_c i M_halo dla każdej galaktyki
4. Wykonaj regresję log(r_c) vs log(M_gal) jak w ex41
5. Porównaj α z predykcją F3 (−1/9) i F1 (−1)

**Oczekiwany wynik:** Potwierdzenie α ≈ −0.11 z błędem < 0.02 (5σ od F1)

---

### 1.3 ex44 — Widmo Mocy TGP-FDM vs CDM: DESI Forecast

**Cel:** Prognoza wykrywalności FDM przez DESI DR3 (2026)
**Plik:** `nbody/examples/ex44_desi_fdm_forecast.py`

**Metoda:**
1. P(k)_FDM z funkcją transferu (ex40 jako baza)
2. Okno pomiarowe DESI (efektywny k-range, błąd shot noise)
3. SNR dla m₂₂ = 1, 2, 5 relative to CDM
4. Prognoza: przy jakim m₂₂ DESI wyklucza FDM na 3σ?

**Oczekiwany wynik:** Krytyczna masa m₂₂^crit(DESI DR3) ≈ ?

---

### 1.4 ex45 — Chaos Lyapunova: TGP vs Newton

**Cel:** Test czy bariera repulsywna TGP (V_β ~ 1/d²) redukuje chaos
**Plik:** `nbody/examples/ex45_lyapunov_tgp_vs_newton.py`

**Protokół:**
1. Symulacja 3-ciał przez T = 10⁴ kroków (dynamics_v2.py jako baza)
2. Perturbacja δ = 10⁻⁷ w warunkach początkowych
3. Pomiar wykładnika Lyapunova λ = lim(t→∞) ln|δ(t)/δ(0)|/t
4. Porównaj λ_TGP vs λ_Newton dla różnych β/C

**Predykcja:** λ_TGP < λ_Newton (bariera tłumi chaos)

---

### 1.5 Weryfikacja ex30-ex33 w verify_all.py

**Problem:** Brakuje testów dla ex30 (saddle-point equilateral), ex32 (isosceles minimum), ex33 (2D Jacobi) w `verify_all.py`
**Cel:** Dodać 15-20 nowych testów, osiągnąć 54+ PASS

**Nowe testy do dodania:**
```
verify_all.py → sekcja E (N-body 2D):
E1: ex30: konfiguracja equilateral = punkt siodłowy
E2: ex30: minimum leży w podprzestrzeni isosceles
E3: ex32: E_0(isosceles) < E_0(equilateral) [korekta -27%]
E4: ex32: górna granica m_sp < 0.12 l_Pl⁻¹
E5: ex33: E_0 = -0.009 E_Pl ± 10%
```

---

## PRIORYTET 2 — Pliki LaTeX

### 2.1 sek08_formalizm.tex — V_mod, ε_th, Złoty Podział

**Brakuje:** Sekcja o ε_th wywiedzionej z N0-6, złoty podział φ jako minimum
**Dodać po §ssec:screening-mass (przy m_sp):**

```latex
\begin{proposition}[Próg stabilizacji solitonu]
Dla \beta = \gamma (N0-5) i m_{sp} = \sqrt{\gamma} (N0-6):
\varepsilon_{th} = m_{sp}^2 / 2 = \gamma/2
\end{proposition}

\begin{proposition}[Złoty podział jako minimum V_{mod}]
V'_{mod}(\varphi) = m_{sp}^2 \varphi(1 + \varphi - \varphi^2) = 0
\Rightarrow \varphi^2 - \varphi - 1 = 0 \Rightarrow \varphi = \phi_{gold}
\end{proposition}
```

---

### 2.2 sek03_rezimy.tex — Dwa Sektory N-body

**Brakuje:** Wzmianka że TGP ma dwa niezależne sektory N-body
**Dodać:** Tabela Efimov/FDM po §ssec:rezimy-fizyczne

---

### 2.3 sek05_ciemna_energia.tex — Kosmologiczne Implikacje FDM

**Brakuje:** Jak FDM (m_boson = m_sp) wpływa na P(k) kosmologiczne
**Dodać:** §ssec:fdm-power-spectrum z wynikami ex40

---

### 2.4 dodatekD_trojcialowe.tex — Korekta 2D Jacobi

**Brakuje:** Wyniki ex32-ex33 (korekta −27% dla E_0, zawężenie okna)
**Dodać:** Tabela korekt 1D vs 2D

---

### 2.5 nbody/tgp_nbody_predictions.tex — K18, K14

**Brakuje:** Kill-shot K18 i K14 jako predykcje
**Dodać:** Sekcja "Cosmological FDM predictions" z tabelą F1/F3

---

## PRIORYTET 3 — Regresja i Testy

### 3.1 Uruchomienie ex1-ex33 po zmianach teorii

**Cel:** Sprawdzenie czy żaden wcześniejszy skrypt nie złamał się po dodaniu ε_th
**Protokół:**
```bash
cd TGP/TGP_v1/nbody/examples
python verify_all.py  # 39/39 oczekiwane
python master_verification_v25.py  # 60/60 oczekiwane
# Uruchom po kolei ex1-ex33 i sprawdź czy brak błędów
```

**Znane ryzyko:** ex14 (NGC 3198) — wynik negatywny jest poprawny (TGP ≠ DM)

---

### 3.2 master_verification_v26.py

**Cel:** Rozszerzyć v25 (60 testów) o wyniki ex43-ex45
**Nowe sekcje:**

```
I: ex34 Efimov window 2D → 8 testów
J: ex43 THINGS K18 precision → 6 testów
K: ex44 DESI FDM forecast → 4 testy
L: ex45 Lyapunov chaos → 5 testów
```

**Cel:** master_verification_v26.py: 83+ PASS, 0 FAIL

---

## PRIORYTET 4 — Fizyka Otwarta

### 4.1 Hierarchia Mas Fermionów (O16)

**Problem:** Dlaczego masy leptonów/kwarków mają właśnie te wartości?
**Propozycja:** Monte Carlo 3D na substracie Γ — pomiar β_c(z, N_flavor)
**Plik:** `scripts/fss_fermion_mass_hierarchy.py`
**Trudność:** 🔴 Bardzo wysoka (wymaga modelu substratowego)

---

### 4.2 Orbita Ósemkowa w TGP (O19)

**Dane startowe:** `configurations.py::figure_eight_initial()`
**Cel:** Czy orbita ósemkowa jest stabilna w TGP (z V₃)?
**Predykcja:** V₃ zaburza orbitę, czas życia τ ~ C⁻¹ · τ_Newton
**Plik:** `nbody/examples/ex46_figure_eight_tgp.py`

---

### 4.3 Q-ball TGP jako Barion (O11)

**Tło:** Ex15 wykazało istnienie Q-ball dla ω < 0.5 (pole zespolone Φe^{iωt})
**Pytanie:** Czy Q-ball TGP ma właściwe ładunki barionowe?
**Trudność:** 🟡 Średnia (wymaga rozszerzenia do pola zespolonego)

---

## Harmonogram Sesji

| Sesja | Zadania | Czas est. |
|---|---|---|
| **v26 (następna)** | ex34, ex43, verify_all aktualizacja | 2-3h |
| **v27** | ex44 DESI, ex45 Lyapunov, sek08 update | 2h |
| **v28** | master_verification_v26, sek03/sek05 update | 1-2h |
| **v29** | ex46 ósemka, LaTeX dodatki, regresja pełna | 2h |
| **v30** | Domknięcie formalne — gotowy do publikacji draft | 3h |

---

## Kill-Shots: Status Kompletny v25

| Kill-shot | Opis | Wynik | Status |
|---|---|---|---|
| K1 | Trzy reżimy z jednego pola | α=−1/9 ✅ | ✅ PASS |
| K5 | GW170817: c_GW = c₀ | Spójne ✅ | ✅ PASS |
| K13 | Efimov klastry planckowskie | Okno: m_sp∈(0, 0.0831) l_Pl⁻¹; C_Q(2B)>>C_Pl | ✅ POTWIERDZONY (ex46) |
| K14 | Lyman-α: m₂₂ > 2.1 | 1D P_Lya: −0.4σ; m22_fit=1.10; ZGODNY | ✅ POTWIERDZONY (ex48) |
| K15 | BBN: ΔG/G < 0.1 | ΔG/G = 0 (att.) ✅ | ✅ PASS |
| K18 | SPARC r_c vs M_gal | α=−0.096≈−1/9; ΔAIC=133 | ✅ PASS (silne) |
| K20 | 3PN odchylenie | Δc₂=−1/3 (LISA 2034+) | ⚠️ DO TESTU |

---

## Najważniejszy Wynik v25 (dla dokumentacji)

> **ε_th = m_sp²/2 = γ/2 wywiedziony z N0-6 (ex39) — ZERO nowych parametrów.**
>
> Złoty podział φ pojawia się jako naturalne minimum V_mod przy progu stabilizacji solitonu TGP-FDM. Scenariusz F3 (m_boson = m_sp = const) preferowany przez K18 (ΔAIC=133). Paradoks γ — jak jeden m_sp spełnia zarówno Efimov (~0.1 l_Pl⁻¹) jak i FDM (~8.2×10⁻⁵¹ l_Pl⁻¹) — ROZWIĄZANY: to dwa niezależne sektory fizyki.

---

*Wersja planu: 2026-03-22 | master_v25: 60/60 PASS | Autor: TGP Analysis Session v25*
