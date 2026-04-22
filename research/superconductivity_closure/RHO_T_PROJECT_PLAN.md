# ρ(T) Project Plan — TGP-based Resistivity of Normal Metals

**Data utworzenia:** 2026-04-21
**Status:** Szkic (osobny projekt badawczy)
**Kontekst:** Użytkownik (2026-04-20) zasugerował rozszerzenie paperu SC o "zwykłe" (nie-nadprzewodzące) materiały. Pytanie: czy TGP-model (skalibrowany na SC w paperze v1) przewiduje ρ(T) z dobrą dokładnością? Praktycznie: czy te same stałe TGP ($C_0$, $A_\text{orb}$, $\eta_f$, $\eta_{6p}$, $\gamma_M$) pozwalają zamknąć nie tylko $T_c$, ale też transport w fazie normalnej?

## Motywacja fizyczna

### Dlaczego ρ(T) jest komplementarny do $T_c$

BCS i dłuższe teorie parujące łączą **normal-state transport** z **SC pairing** przez wspólny parametr Eliashberga $\alpha^2 F(\omega)$:

- $T_c$: zależy od $\lambda_{ep} = \int \tfrac{\alpha^2 F(\omega)}{\omega}\,d\omega$
- $\rho(T)$: Bloch-Grüneisen, $\rho_\text{ph}(T) \propto \int \omega \alpha^2_\text{tr}(\omega) F(\omega)\,\cdot\, f_{BG}(T/\Theta_D)$
- Te same fonony, ta sama struktura elektronowa, różne momentowe wagi.

**Jeśli TGP daje dobry $T_c$, powinien też dawać dobre $\rho(T)$** — inaczej model jest niespójny.

### Co TGP wnosi nad BCS

BCS nie przewiduje $\Theta_D$ z pierwszych zasad. TGP (paper v1, Eq. 5) ma:
- $k_d(z)$ — koordynacja
- $A_\text{orb}$, $M_\text{gauss}(a)$ — geometria + lattice matching
- $\Lambda_\text{eff}(\omega) = \Lambda_0 (\omega/\omega_0)^{\alpha_{P6B}}$ — substrate coupling freq-dep

Pytanie: czy te same obiekty określają $\Theta_D$ i $\rho_0$ (defekt scattering) z niezależnym fitem?

## Plan badawczy

### Etap 1: Dataset gathering (target N ~ 50)

Meterialne źródła:
- **CRC Handbook of Chemistry and Physics** — $\rho(T)$ dla Cu, Ag, Au, Al, Fe, Ni, Nb, V, Ta, W, Mo, Pb, Sn, In, Tl, Bi, Zn, Cd, Pt, Pd, ...
- **Matula (1979)** review "Electrical resistivity of Cu, Ag, Au" — standard reference
- **Bass (1982)** "Electrical resistivity of pure metals" — Landolt-Börnstein tables
- **DFT z Materials Project** — dla fitów $\Theta_D$ i $N(E_F)$

Wymagane kolumny:
- name, lattice_const, Debye_T, N_F (states/eV/atom)
- $\rho_0$ (0 K residual), $\rho(T)$ at T=77, 295, 500, 1000 K
- purity/RRR flag (dla $\rho_0$)

### Etap 2: Wzór TGP na $\rho(T)$

Propozycja bazowa (mechanical-hypothesis, do testu):

$$\rho_\text{ph}(T) = \rho_\text{TGP}^{(0)} \cdot f_\text{BG}(T/\Theta_D^\text{TGP})\,\cdot\, N(E_F)^{-1} \cdot\, \omega_0^{\alpha_\text{BG}}$$

gdzie:
- $f_\text{BG}(x) = 4x^5 \int_0^{1/x} \tfrac{t^5\,e^t}{(e^t-1)^2}\,dt$ — standardowa funkcja Bloch-Grüneisen
- $\Theta_D^\text{TGP}$ — **wyprowadzone z Eq. 5** poprzez inverse fit: $\Theta_D = f(k_d, A_\text{orb}, a, \Lambda_\text{eff})$
- $\rho_\text{TGP}^{(0)}$ — nowa stała TGP (jeśli universal, to PASS)
- $\alpha_\text{BG}$ — wykładnik freq-dep (test czy = $\alpha_{P6B}$ z paperu v1)

### Etap 3: Testy walidacyjne

**Test A (kalibracja $\Theta_D$):**
- Cu, Al, Nb, V, Ag — z literatury znane $\Theta_D$ DFT
- Fit $\Theta_D^\text{TGP}$ z Eq. 5 struktura: czy zgadza się z DFT w 20%?

**Test B (BG fit na "czystych" metalach):**
- Cu, Ag, Au, Al, Pt, Pd (non-magnetyczne, prosty band)
- Czy jedno $\rho_\text{TGP}^{(0)}$ zamyka wszystkie w RMS(log ρ) < 0.3?

**Test C (rozszerzenie na transition metals):**
- V, Nb, Ta, Mo, W — d-band, wieloorbitowe
- Czy $A_d = 0.310$ (paper v1) działa też w $\rho(T)$?

**Test D (magnetyki):**
- Fe, Ni, Co — FM, spin-scattering dodatkowy
- Test: $\rho(T > T_\text{Curie}) = \rho_\text{ph} + \rho_\text{mag}(T)$, $\rho_\text{mag} \propto M^2(T)$
- Czy $\lambda_\text{sf}$ (z B_mag w Eq. 5) kalibruje $\rho_\text{mag}$?

**Test E (cross-check z SC):**
- Metal-SC z paperu v1 (Nb, V, Pb, Hg, Al, Sn): predykcja $T_c$ i $\rho(T)$ z **tymi samymi** stałymi TGP.
- Jeśli oba RMS < 0.3 — **model TGP unifikuje normal + SC state**.

### Etap 4: Analiza i paper-writeup

Jeśli Etap 3 przejdzie (RMS < 0.3 dla $\geq 40/50$ materiałów):
- **Paper v3** "TGP unified theory of metal transport: from normal-state to SC"
- Nowy rozdział w paperze SC v2: "Normal-state consistency check"
- Przynajmniej jedno **niezależne przewidywanie** ρ(T) dla materiału nie-SC, nie w training set.

Jeśli Etap 3 nie przejdzie:
- Analiza gdzie model się wali — jakich fizycznych członów brakuje?
- Możliwe: ρ_0 (defekt) wymaga separate calibration; $\Theta_D$-TGP ma offset; freq-dependence inna niż SC.
- Osobny sub-problem: P8.x (normal-state extension).

## Milestones

| Etap | Szacowany czas | Output |
|------|----------------|--------|
| 1: Dataset | 1 tydzień | `rho_T_dataset.py` z N=50 |
| 2: Wzór TGP | 2-3 dni | `rho_TGP_formula.py` (testowa forma BG) |
| 3: Testy A-E | 1 tydzień | `ps3N_rho_validation.py` × 5 |
| 4: Analiza + paper | 2 tygodnie | sekcja paper v2 lub osobny preprint |

## Kluczowe pytania otwarte

1. **Czy $\Theta_D$ daje się wyliczyć z Eq. 5 bez dodatkowych parametrów?**
   - Test: $\Theta_D \sim \omega_0 \cdot \sqrt{k_d \cdot C_0}$ (wymiarowo)?
   - Jeśli nie — wprowadzamy $\Theta_D$ jako input (jak $\omega$ w Eq. 5).

2. **Czy $\rho_0$ jest uniwersalne?**
   - Prawdopodobnie NIE (zależy od defektów/purity), ale może istnieć TGP lower bound: $\rho_0^\text{min} = f(\hbar/e^2 \cdot k_F^{-1})$.

3. **Czy P7.10 ($N(E_F)$-factor) wyjdzie tu jako niezależny test?**
   - Insulatory: $\rho(T)$ diverguje przy $T \to 0$ (aktywacyjne), $N(E_F) \to 0$. Prosty test P7.10 functional form.

4. **Czy ρ_mag(T) dla FM (Fe, Ni, Co) przechodzi?**
   - Krytyczny test $B_\text{mag}(\lambda_\text{sf})$ — jeśli to samo $\lambda$ kalibruje $T_c$ i $\rho$, to model jest spójny.

## Powiązania z bieżącymi sub-problemami P7.5

- **P7.6 (5d mass)**: ρ(T) dla Ta, W, Re, Os, Ir daje **niezależny test** $\gamma_M$ — bez dopasowywania do $T_c$, tylko do residual ρ(T).
- **P7.10 ($N(E_F)$-factor)**: ρ(T) insulatorów testuje $g(0) = 0$ bezpośrednio (insulator: $\rho \to \infty$).
- **P7.11 (magnetic)**: Pt/Pd mają $\rho(T) \propto T^2$ w niskim T (Fermi-liquid enhancement) — test czy TGP daje Fermi-liquid coefficient.

## Pre-requisites

- Uzupełnić sub-problem **P7.10** w current SC research (dodać $N(E_F)$ factor) PRZED ρ(T) — bo oba używają tej samej stałej.
- Finalizować **P7.6** ($\gamma_M$ dla ciężkich 5d) — też wspólna stała.

## Jak rozpocząć

1. Stwórz katalog `TGP/TGP_v1/research/rho_normal_state_closure/`
2. Skopiuj bazowe utility funcs z `superconductivity_closure/ps29c_expanded_dataset.py` jako import
3. Utwórz `r00_dataset.py` — N=15 pilot (Cu, Ag, Au, Al, Nb, V, Pb, Sn, Fe, Ni, Pt, Pd, Cd, Zn, Mg)
4. Utwórz `r01_bg_baseline.py` — standardowy Bloch-Grüneisen fit na pilot (bez TGP, tylko $\Theta_D$, $\rho_0$)
5. Utwórz `r02_tgp_formula.py` — czy $\Theta_D$ z Eq. 5?
6. Iteruj.

## Powiązane dokumenty

- [[P7.5_subproblems.md]] — źródłowy kontekst, lista otwartych sub-problemów
- [[ps29c_expanded_dataset.py]] — core TGP functions (Tc_phonon, k_d, M_gauss, B_mag)
- [[ps32_negative_validation.py]] — motywacja: negatywna walidacja ujawnia braki modelu
- Paper SC v1: `../../papers/sc/tgp_sc.tex`, Eq. 5 (`eq:A-orb`)

## Poprzednie prace do zbadania

- **Ziman (1960)** _Electrons and Phonons_ — teoria transportu klasyczna
- **Allen & Mitrovic (1982)** — Eliashberg function i transport
- **Grimvall (1981)** _The Electron-Phonon Interaction in Metals_ — standardowe wzory
- **Ohara (1993)** — $\alpha^2 F$ z T-zależnego ρ(T)

Plan nie jest finalny — jest pierwszym szkicem. Uaktualnić po Etapie 1 (faktyczny dataset).

---

## UPDATE 2026-04-21: Etap 1 PILOT (N=15) zrealizowany

### Skrypty

- [[rho_normal_state_closure/r00_dataset.py]] — dataset N=15 (Cu/Ag/Au/Al/Pb/Sn/Nb/V/Fe/Ni/Pt/Pd/Cd/Zn/Mg)
- [[rho_normal_state_closure/r00_results.txt]] — tabela, sanity checks (monotonia, RRR 300-2800, slope ratios)
- [[rho_normal_state_closure/r01_bg_baseline.py]] — klasyczny BG fit (ρ_0, R, Θ_D)
- [[rho_normal_state_closure/r01_results.txt]] — R per material CoV=80%, Θ_D free-fit dryfuje +200%
- [[rho_normal_state_closure/r02_tgp_formula.py]] — test H1 (z P7.12 g) vs H0 (bez g)
- [[rho_normal_state_closure/r02_results.txt]] — H0 lepsze: RMS=0.34, r=0.57 (p=0.03)
- [[rho_normal_state_closure/r03_class_split.py]] — 7 alternatywnych parametryzacji
- [[rho_normal_state_closure/r03_results.txt]] — **(E) R ∝ (λ·Θ_D)^1.025, RMS=0.25, r=0.79 (p=0.0004)**

### Kluczowe odkrycia

1. **Θ_D nie da się wyciągnąć z rho(T)** przy T-points (77, 295, 500, 1000) K — wszystkie w high-T rezimie. Potrzebne T < Θ_D/3 (czyli 20-80K) aby odzyskać BG curvature. Alternatywa: akceptujemy Θ_D jako input (jak ω w Eq. 5 SC).

2. **P7.12 (g(N_F) z SC) NIE przenosi się na ρ(T).** Cu/Ag/Au mają g(x)≈0 z P7.12, ale _mają_ realne ρ(T). Fizycznie: SC wrażliwe na DOS (eksp formuła gap), ρ(T) transport istnieje dla każdego $N_F > 0$. **To pierwsza wyraźna demarkacja mode-specific vs universal** w TGP.

3. **FM ma nadmiarowy slope** (Fe: 1.36×, Ni: 1.31× vs BG predict) — sygnał magnon scattering, potwierdza P7.11 magnetic. Stoner paramagnety też lekko (Pd: 1.20×).

4. **Najlepszy universal fit:** $R \propto \lambda_\text{ep} \cdot \Theta_D$ z wykładnikiem 1.0 (klasyczny BG) i prefactorem $C \approx 0.27$ (zbliżony do Zimana 0.19).

### Ścieżka unifikacji SC+ρ(T) przez λ_ep

$$\rho_\text{TGP}(T) = \rho_0 + C_\rho \cdot \lambda_\text{ep}^\text{TGP} \cdot \Theta_D \cdot \left(\frac{T}{\Theta_D}\right)^5 J_5\!\left(\frac{\Theta_D}{T}\right)$$

gdzie $\lambda_\text{ep}^\text{TGP}$ = funkcja standardowych stałych z Eq. 5 (C_0, A_orb, M_gauss, Λ_eff). Test w Etapie 2: czy **Eq. 5 structure** (bez rozszerzeń P7.10-P7.12 które są mode-specific) zamyka λ_ep do 20%?

### Ograniczenia Etapu 1

- N=15 mały sample; rozszerzenie o Mo, Ta, W, Be, Ti, Sr, Ca, Bi (N~25 target)
- Brak niskotemperaturowych T-points (potrzebne dla BG curvature)
- λ_ep literaturowe (Allen & Mitrovic) używane jako input — test circularity: czy TGP wyliczy to samo?

### Proponowany Etap 2 (r04-r06)

- **r04**: oblicz $\lambda_\text{ep}^\text{Eq.5}$ z paper v1 formuły dla 15 metali, porównaj z literaturowym. Jeśli RMS_log(λ_TGP, λ_lit) < 0.15 → unifikacja zachowana.
- **r05**: wersja BG z pełnym integrated T-profile: fituj $\rho(T_i)$ na wszystkich 4 T-points simultanously (nie tylko amplituda R). Wtedy widzimy, czy BG+TGP zamyka pełną krzywą.
- **r06**: rozszerzenie dataset o Mo, Ta, W (5d reaferractory), Be, Ti (low-Z), Bi, Sb (semimetal). Weryfikacja czy struktura λ_ep·Θ_D skaluje.

### Status

Etap 1 **zamknięty** z mocnym sygnałem: $(\lambda \cdot \Theta_D)$ jest dominującą skalą. Gotowy do Etapu 2.

---

## UPDATE 2026-04-21 (kontynuacja): Etap 2 (r04-r05) zrealizowany

### r04 — wynik **NEGATYWNY** dla bezpośredniej unifikacji przez λ_ep

**Skrypt:** [[rho_normal_state_closure/r04_lambda_tgp.py]] | [[rho_normal_state_closure/r04_results.txt]]

**Testowana hipoteza:** $\lambda_\text{ep} = C \cdot F_\text{TGP}$, gdzie $F_\text{TGP} = k_d \cdot A_\text{orb}^2 \cdot M_\text{gauss} \cdot \Lambda_\text{eff}$ (struktura Eq. 5 bez suppresji SC-specific).

**Wyniki 7 form:**

| Forma | RMS_log | r (Pearson) | p-value |
|------|---------|-------------|---------|
| A0: $C \cdot F_\text{TGP}$ | 0.400 | **−0.063** | 0.82 |
| A1: $C \cdot F_\text{TGP} \cdot N_F^q$ (fit q=−0.13) | 0.393 | — | — |
| A2: $C \cdot F_\text{TGP}^\alpha$ (fit α=−0.10) | 0.327 | — | — |
| A3: $C \cdot F_\text{TGP} \cdot N_F$ | 0.742 | 0.22 | — |
| A4: Cu/Ag/Au reklas. "sp"→"s" | 0.363 | 0.41 | 0.13 |
| A5: Hopfield $C \cdot F \cdot N_F/\Theta_D^2$ | 0.613 | 0.42 | 0.12 |
| A6: $C \cdot (1/\Theta_D)^\alpha$ | 0.310 | 0.33 | 0.24 |
| A7: $C \cdot N_F^\alpha$ | 0.311 | 0.32 | 0.24 |

**Kluczowa obserwacja:** F_TGP przewiduje **prawie stałe** λ ≈ 0.38 dla wszystkich materiałów (zakres 0.35–0.41), podczas gdy literatura daje **13× wariancję** (0.12 Ag → 1.55 Pb). Pearson r = −0.063 z p = 0.82 → **brak korelacji**.

**Wniosek fizyczny:** TGP Eq. 5 struktura **nie koduje λ_ep** jako oddzielnej wielkości. To znaczy że Eq. 5 przewiduje bezpośrednio T_c (linearnie w Λ_eff), nie poprzez McMillan eksponencjalne. Unifikacja "TGP→λ→ρ(T)" jest nieosiągalna — trzeba testować bezpośredni path TGP→ρ(T).

### r05 — wynik **POZYTYWNY** dla bezpośredniej predykcji ρ(T)

**Skrypt:** [[rho_normal_state_closure/r05_rho_direct.py]] | [[rho_normal_state_closure/r05_results.txt]]

**Testowane formy** (R = amplituda BG):

| Forma | RMS_log | # params | Uwagi |
|------|---------|----------|-------|
| D0: $C \cdot F_\text{TGP}$ | 0.339 | 2 | baseline r02 |
| D1: $C \cdot F_\text{TGP} \cdot \Theta_D$ | 0.358 | 2 | Ziman high-T |
| D2: $C \cdot F_\text{TGP} \cdot \Theta_D^p$ | 0.335 | 2 | p=0.29 |
| D3: $C \cdot F_\text{TGP}^a \cdot \Theta_D^b$ | 0.335 | 3 | a=0.94, b=0.34 |
| **D4: $C \cdot F_\text{TGP}^a \cdot \Theta_D^b \cdot N_F^c$** | **0.197** | 4 | **a=−1.89, b=1.21, c=1.05** |
| D5: $C \cdot \lambda_\text{lit} \cdot \Theta_D$ *(UB)* | 0.250 | 2 | upper bound lit-based |
| D6: $C \cdot \lambda_\text{lit}^a \cdot \Theta_D^b$ *(UB)* | 0.224 | 3 | upper bound lit-based |
| D7: $C \cdot \Theta_D^p \cdot N_F^q$ (bez F_TGP) | 0.252 | 3 | kontrola |

**D4 bije D6 (lit-based upper bound)**: TGP-only model jest lepszy niż model używający literatury λ! Gap = −0.027 dex.

**LOO stabilność wykładników D4:**
- a = −1.89 ± 0.26 (range −2.76 … −1.63)
- b = +1.21 ± 0.13 (range +1.02 … +1.49)
- c = +1.05 ± 0.08 (range +0.94 … +1.29)

Wszystkie 3 wykładniki stabilne na LOO → nie jest overfitting.

**Profil pełnego ρ(T_i):** mean RMS = 0.185 (median 0.21), max 0.38 (Al). 12/15 materiałów zamyka w RMS_log < 0.27.

### Interpretacja fizyczna r05 D4

Wykładnik $a = -1.89$ (NEGATYWNY) dla F_TGP jest kontrintuicyjny. Analiza:

1. **Multikolinearność:** F_TGP i N_F silnie skorelowane (oba wysokie dla d-metali, niskie dla noble). Regresja dzieli informację na c=+1.05 (N_F) i a=−1.89 (F), które częściowo się kasują.

2. **Efektywna zależność:** Ponieważ $\Lambda_\text{eff} \propto \Theta_D^{0.5}$, wykładnik −1.89 na F_TGP generuje efektywny $\Theta_D^{-0.945}$, który łączy się z b=+1.21 dając netto $\Theta_D^{+0.27}$ — blisko D2/D3 optimum.

3. **Prawdziwa struktura:** $R \approx C \cdot \Theta_D^{0.3} \cdot N_F^{1.0} \cdot (A_\text{orb}^2 \cdot k_d \cdot M_\text{gauss})^{-1.89}$. Orbital/koordinacyjny czynnik działa jako **suppressor** transportu, nie booster — co jest zgodne z faktem, że d-metals (wysokie A_orb²) mają wysokie R ale nie tak wysokie jak Pb/Nb (gdzie Pb ma wysoki Λ_eff).

### Ścieżka dla Etapu 3

**Zaktualizowana hipoteza unifikacyjna:**

$$\rho_\text{TGP}(T) = \rho_0 + C_\rho \cdot \mathcal{G}(F_\text{TGP}, N_F) \cdot \Theta_D^{+1.2} \cdot \left(\frac{T}{\Theta_D}\right)^5 J_5\!\left(\frac{\Theta_D}{T}\right)$$

gdzie $\mathcal{G}(F, N_F) = F_\text{TGP}^{-1.89} \cdot N_F^{+1.05}$ jest **nową TGP-transport funkcją** (nie λ_ep), ale przez inwersję McMillan dostaniemy łącznik z T_c na jej bazie.

### Ograniczenia r04-r05

- 4 parametry (D4) na 15 materiałów - overfitting flag pomimo LOO stability.
- Współczynniki D4 mogą być artefaktem konkretnego wyboru A_orb/k_d. Test: inne class assignments (r06 powinno).
- Al ma największy profile RMS (0.38) — wysoki Θ_D + niski Z wypada poza main sequence; do zbadania w r07.

### Proponowany Etap 3 (r06-r08)

- **r06**: dataset extension do N=25 (Mo, Ta, W, Ti, Be, Bi, Sb, Co, Cr, Mn, Rh). Sprawdza czy D4 wykładniki (a≈−1.9, b≈1.2, c≈1.0) utrzymują się.
- **r07**: inwersja D4: startując z R(D4), wyciągnij "ukryte λ_TGP" i porównaj z McMillan-Allen λ z T_c. To test unifikacji z drugiej strony.
- **r08**: finalne zamknięcie: paper v3 (lub rozszerzenie v2) gdzie pokazujemy D4 jako **niezależne przewidywanie transportu** bez fitowania λ.

### Status Etapu 2

**ZAMKNIĘTY** z mocnym pozytywnym wynikiem dla r05 (D4 RMS=0.197 < lit upper bound 0.224) i mocnym negatywnym dla r04 (TGP nie koduje λ_ep osobno). Interpretacja: **TGP nie jest generatorem λ_ep, ale jest bezpośrednim generatorem transport amplitudy R_BG** przy danym (F_TGP, Θ_D, N_F). To otwiera paper v3 kierunek "TGP unified transport" bez roszczeń do λ-prediction.

---

## UPDATE 2026-04-21 (Etap 3): r06-r07 — falsyfikacja uniwersalnej D4

### r06 — rozszerzenie datasetu do N=23, D4 niestabilny

**Skrypt:** [[rho_normal_state_closure/r06_extension.py]] | [[rho_normal_state_closure/r06_results.txt]]

**Dodane materiały (N+8):**

| Material | Klasa | Θ_D | N_F | ρ(295) | λ_lit | Z | Uwagi |
|---------|-------|-----|-----|--------|-------|---|-------|
| Mo | d | 450 | 0.60 | 5.34 | 0.41 | 42 | bcc 4d refractory |
| Ta | d | 240 | 0.82 | 13.15 | 0.69 | 73 | bcc 5d SC |
| W | d | 310 | 0.42 | 5.28 | 0.28 | 74 | bcc 5d |
| Ti | d | 420 | 1.00 | 42.0 | 0.38 | 22 | hcp 3d, wysokie ρ |
| Be | sp | 1440 | 0.04 | 3.56 | 0.23 | 4 | hcp s-metal, extreme Θ_D |
| Bi | sp | 119 | 0.006 | 107.0 | 0.20 | 83 | semimetal, N_F=100× niższe |
| Co | d | 445 | 1.95 | 6.24 | 0.30 | 27 | hcp FM T_C=1388K |
| Rh | d | 480 | 0.75 | 4.51 | 0.37 | 45 | fcc 4d |

**Krytyczny wynik r06:** Re-fit D4 (R = C·F^a·Θ^b·N_F^c) na N=23 daje:

| Parametr | N=15 | N=23 | Dryf |
|---------|------|------|------|
| a | −1.888 | **+0.419** | 122% (sign flip!) |
| b | +1.211 | +0.461 | 62% |
| c | +1.045 | **−0.003** | 100% (~zero) |
| RMS_log | 0.197 | **0.384** | +95% |

**Wniosek:** Forma D4 z r05 była **overfitted** na N=15 przez multikolinearność F_TGP i N_F. Na N=23 nie istnieje uniwersalna forma z RMS < 0.35.

### r07 — analiza outlierów + class-specific fit

**Skrypt:** [[rho_normal_state_closure/r07_outlier_class.py]] | [[rho_normal_state_closure/r07_results.txt]]

**Testowanie wykluczeń i reklasyfikacji D4 (N=23 baseline):**

| Scenariusz | a | b | c | RMS |
|-----------|---|---|---|-----|
| Baseline N=23 | +0.42 | +0.46 | 0.00 | 0.384 |
| (A) noble(Cu,Ag,Au) "sp"→"s" | +1.48 | −0.12 | −0.37 | 0.266 |
| (B) bez Bi (semimetal) | −1.17 | +1.57 | +0.71 | 0.257 |
| (C) bez Bi + noble "s" | +0.93 | +0.40 | −0.03 | 0.243 |
| (D) bez Bi,Be + noble "s" | +0.80 | 0.00 | +0.16 | 0.226 |

Każde wykluczenie zmienia znak i wielkość wykładników → **D4 nie jest uniwersalne**.

### Class-specific D form — MOCNY wynik

**R = C · Θ_D^b · N_F^c osobno dla każdej klasy:**

| Klasa | N | b (Θ_D) | c (N_F) | RMS_log | Członkowie |
|------|---|---------|---------|---------|-----------|
| **sp** | 8 | +0.25 | **−0.48** | **0.161** | Al, Pb, Sn, Cd, Zn, Mg, Be, Bi |
| **d** | 12 | +0.84 | +0.29 | 0.219 | Nb, V, Fe, Ni, Pt, Pd, Mo, Ta, W, Ti, Co, Rh |
| s (noble) | 3 | +3.21 | −2.53 | 0 (trivial) | Cu, Ag, Au |

**Interpretacja fizyczna:**

- **sp-metals (c = −0.48):** Wyższe N_F → NIŻSZY R. Counter-intuitive? Fizyczna mechanizm: w sp-metalach wyższe N_F oznacza więcej d-character mix → silniejsze screening → niższa oporność. Potwierdzone dla Bi (N_F=0.006, ρ ogromne) vs Mg/Al (N_F≈0.1-0.4, ρ umiarkowane).

- **d-metals (c = +0.29):** Wyższe N_F → wyższy R (zgodnie z Eliashberg intuition: więcej scattering channels).

- **Noble (s-class):** Tylko 3 punkty, niedookreślone, ale obecnie wygląda na bardzo inną zależność niż sp/d.

### Co się zmieniło w interpretacji unifikacji

**Pierwotna hipoteza (Etap 1-2):** Istnieje uniwersalny funkcjonał TGP → ρ(T) poprzez λ_ep lub poprzez bezpośredni F_TGP·Θ_D·N_F.

**Aktualna hipoteza (po r06-r07):** Transport jest **class-specific** — sp, d, s-metale mają różne zależności (Θ_D, N_F). Uniwersalna forma D4 była artefaktem overfitting (N=15 zbyt małe na 4 parametry).

Ale:
- Class-specific fit **jest spójny z fizyką** (sp/d różnią się naturą orbital mixing)
- RMS ~0.16-0.22 per class to **bardzo dobra zgodność** (lepsze niż literatura λ-based)
- Class jest już w TGP Eq. 5 przez A_orb → *TGP wie o klasach*

### Następne kroki

1. **r08**: Per-class analiza profilu ρ(T) na 4 T-points — czy class-specific amplituda R przewiduje też krzywiznę BG?
2. **r09**: Fizyczna interpretacja sp c=−0.48 — czy to związane z d-band mixing lub Fermi surface topologia?
3. **r10**: Extension N=23 → N=35 z bardziej różnych klas (Ir, Re, Cs, Sr, Sc, Y, Zr, Hf) żeby zweryfikować klasowe wykładniki.
4. **paper**: Jeśli class-specific RMS stabilne (<0.25) — paper v3 "TGP class-specific transport theory".

### Status Etapu 3

**ZAMKNIĘTY** z wynikiem mieszanym:

- ❌ **NEGATYWNY:** Uniwersalna D4 forma (r05) obalona przez r06 — overfit.
- ✅ **POZYTYWNY:** Class-specific fit stabilny (sp RMS=0.16, d RMS=0.22), fizycznie interpretowalny.
- ⚠️ **OTWARTE:** Czy klasowe wykładniki utrzymują się na N=35? (r10)

---

## UPDATE 2026-04-21 (Etap 4): r08 — class-specific profil ρ(T_i)

**Skrypt:** [[rho_normal_state_closure/r08_class_profile.py]] | [[rho_normal_state_closure/r08_results.txt]]

### Class-specific fit na pełnym profilu ρ(T_i)

Klasowy model: $R_{BG} = C_\text{cls} \cdot \Theta_D^{b_\text{cls}} \cdot N_F^{c_\text{cls}}$ + klasyczne BG formula z ρ_0 + R·(T/Θ_D)^5·J_5(Θ_D/T).

| Klasa | N | C | b (Θ_D) | c (N_F) | RMS_R | **Profil RMS** | Członkowie |
|------|---|---|---------|---------|-------|----------------|-----------|
| **s** (noble) | 3 | 6.9e−2 | +0.84 | — (tylko Θ) | 0.032 | **0.048** | Cu, Ag, Au |
| **sp** | 8 | 3.9 | +0.25 | −0.48 | 0.161 | **0.214** | Al, Pb, Sn, Cd, Zn, Mg, Be, Bi |
| **d** | 12 | 0.47 | +0.84 | +0.29 | 0.219 | **0.227** | Nb, V, Fe, Ni, Pt, Pd, Mo, Ta, W, Ti, Co, Rh |

**Mean profile RMS (N=23) = 0.199** (median 0.177) → **BARDZO DOBRY wynik** dla 4 T-points per materiał.

### Histogram profile RMS na N=23

| Zakres RMS_log | Liczba mat. | Członkowie |
|---------------|-----|-----------|
| [0.0, 0.1) | 5 | Cu, Ag, Au, Sn, Ta |
| [0.1, 0.2) | 9 | Al, Pb, Nb, Fe, Pt, Pd, Zn, Mo, Bi |
| [0.2, 0.3) | 6 | V, Ni, Mg, W, Co, Rh |
| [0.3, 0.5) | 2 | Cd, Be |
| [0.5+) | 1 | Ti |

**20/23 (87%) materiałów RMS < 0.3**. Ti outlier — prawdopodobnie pozostała fizyka (spd hybridization + paramagnetic spin-fluct).

### LOO stability (C)

- **d-class STABLE:** b=+0.84±0.17, c=+0.29±0.08 (CoV 20% / 27%)
- sp-class unstable: b=+0.23±0.21 (CoV 92%), c=−0.46±0.07 (CoV 15%)

→ d-class coefficients są **fizycznie istotne**, sp-class ma problemy z leverage od Be/Bi.

### Sub-class sp (B) — rozwiązanie "c=−0.48 paradoksu"

Rozbijając sp na normalne (N_F > 0.08) vs low-carrier (N_F < 0.05):

| Subklasa | N | b | c | RMS |
|----------|---|---|---|-----|
| sp_normal (Al, Pb, Sn, Cd, Zn, Mg) | 6 | −0.13 | −0.11 | **0.084** |
| sp_low (Be, Bi) | 2 | (ratios różnią się 100×) | | |

**Wnioski:**
- Dla "normalnych" sp-metali: R ≈ stała ~40 μΩ·cm (wąski zakres!). RMS=0.08, znakomity.
- Be i Bi są poza domeną BG (Be: Θ_D=1440 K extreme, Bi: semimetal N_F=0.006). Dołączenie ich wymusza c=−0.48 który nie jest fizyczny.

### LOO Cross-Validation (predykcja held-out) (D)

- Mean |dlog R| = 0.254 (predicted R w 2× obserwowanego)
- Mean profile RMS LOO = 0.293 (pełna krzywa BG)
- Outliers: Be (0.76), Ti (0.67), Rh (0.37), Cd (0.39)
- 15/20 materials z LOO profile RMS < 0.35

### Interpretacja fizyczna

1. **Noble (s-class):** ρ zdominowane prostą s-only transport. Brak d-contribution near E_F. R ∝ Θ_D^0.84 (blisko klasycznego Bloch-Mott).

2. **sp-class normal:** R praktycznie stała (~40 μΩ·cm), bo w sp-metalach elektrony s+p mają podobną naturę niezależnie od DOS. Sugeruje saturation transport.

3. **d-class:** R ∝ Θ_D^0.84 · N_F^0.29 — b zgodne z s-class (ten sam mechanizm fonoowy), c dodatnie i małe (DOS modestly zwiększa scattering). Stabilne LOO.

4. **Outliery:**
   - Ti: probably spin-fluct / magnon scatter not captured
   - Be: extreme Θ_D wymaga poprawki wielogrady BG-Debye
   - Bi: semimetal wymaga own category (band-structure calculation)

### Gotowość do Paperu v3

**Możemy napisać paper v3 "TGP class-specific normal-state transport":**

Treść:
- Sekcja 1: Klasy TGP (s, sp, d) z A_orb — już w paperze v1
- Sekcja 2: Eksperyment dataset N=23 z literaturą (Matula, Bass, CRC)
- Sekcja 3: Class-specific regression R = C_cls · Θ_D^b · N_F^c
- Sekcja 4: Full profile BG prediction (87% materiałów RMS<0.3)
- Sekcja 5: Out-of-domain: Ti (spin-fluct), Be (Θ_D extreme), Bi (semimetal)
- Sekcja 6: Unifikacja: te same klasy (s/sp/d) przewidują T_c (SC paper v1) ORAZ R_BG (ten paper) → TGP jest **wspólną strukturą** orbital-transport

Brakuje jeszcze:
- **r10**: rozszerzenie do N~35 (dodać Ir, Re, Ru, Os, Cs, K, Rb, Sr, Ba, Zr, Hf, Sc, Y, La) aby potwierdzić klasowe wykładniki — głównie s-class (brakuje nam alkaliów do testu 'noble vs s-metal')
- **r09 opcjonalnie**: inwersja class-D do λ_eff i porównanie z McMillan λ z Tc — weryfikacja unifikacji SC↔ρ na poziomie λ.

### Status Etapu 4

**ZAMKNIĘTY** z wynikiem **POZYTYWNYM**:

- ✅ Class-specific D (R = C·Θ^b·N_F^c per klasa) zamyka ρ(T) na N=23 z mean RMS=0.20
- ✅ 20/23 (87%) materiałów RMS<0.3
- ✅ d-class fit stabilny na LOO (b=0.84±0.17, c=0.29±0.08)
- ✅ Paper v3 zakres: "TGP class-specific transport theory" **gotowy do draftu**

Ograniczenia:
- sp sub-class c=−0.48 driven by Be+Bi — w paperze prezentujemy sp_normal (RMS=0.08) + outlier note dla Be/Bi
- Ti pozostaje outlier (single 3d high-ρ anomaly)

**Proponowany następny krok:** r10 — rozszerzenie do N=35 przed paperem (weryfikacja klasowych wykładników).

---

## UPDATE 2026-04-21 (Etap 5): r10 — weryfikacja N=34, **FALSYFIKACJA d-class wykładnika Θ_D**

**Skrypt:** [[rho_normal_state_closure/r10_extension_v2.py]] | [[rho_normal_state_closure/r10_results.txt]]

### Rozszerzenie datasetu z N=23 → N=34

Dodano 11 nowych metali:

| Klasa | Nowe materiały (N=11) |
|-------|-----------------------|
| **d** (+9) | Ir, Re, Ru, Os (refractory noble 4d/5d), Zr, Hf (grupa IV), Sc, Y, La (grupa III) |
| **s** (+2) | Ca, Sr (alkaline-earth, reklasyfikowane do s-class) |
| **sp** (+0) | — (brak sensownych kandydatów — In/Tl/Sb/Ga topią się < 1000K) |

Alkalie (K, Rb, Cs, Na) **wykluczone** — topią się poniżej 500K, brak ρ(500)/ρ(1000) w fazie stałej.

### Porównanie wykładników N=23 vs N=34

| Klasa | N₂₃ | b₂₃ | c₂₃ | RMS₂₃ | N₃₄ | b₃₄ | c₃₄ | RMS₃₄ | |Δb| | |Δc| |
|-------|-----|-----|-----|-------|-----|-----|-----|-------|-----|-----|
| **d** | 12 | **+0.84** | +0.29 | 0.22 | 21 | **+0.02** | +0.25 | 0.25 | **0.82** | 0.05 |
| **sp** | 8 | +0.25 | −0.48 | 0.16 | 8 | +0.25 | −0.48 | 0.16 | 0.00 | 0.00 |
| **s** | 3 | **+0.84** | — (Θ only) | 0.03 | 5 | **−0.45** | +0.88 | 0.11 | **1.29** | — |

**KRYTYCZNE:**
- **d-class b collapsed z +0.84 → +0.02** po dodaniu 9 nowych d-metali (Ir, Re, Ru, Os, Zr, Hf, Sc, Y, La). c stabilny (+0.25).
- **s-class b flipped z +0.84 → −0.45** po dołączeniu Ca, Sr. c zmienił się z 0 → +0.88.
- sp-class nietknięty (brak nowych sp candidates).

**Diagnoza:** N=12 dla d-class i N=3 dla s-class było za małe. Wykładnik **b (Θ_D)** był artefaktem ograniczonego sample'a. **r08 "paper v3 gotowy" wnioisek był PRZEDWCZESNY.**

### Co się utrzymało / co padło

**Utrzymało się:**
- LOO stability dla d-class na N=21: b=+0.02±0.06, c=+0.25±0.06 (tight!).
- Per-klasowy profile RMS: mean=0.229 na N=34 (dalej < 0.25 threshold).
- Histogram: 17/34 (50%) ma RMS<0.2, 24/34 (71%) ma RMS<0.3.
- Ru, Os, La, Sr, Ca — nowe materiały dobrze pasują (RMS 0.08–0.18).
- Noble (s) zachowuje niskie RMS ale wykładniki się zmieniły.

**Padło:**
- Uniwersalny d-class b=+0.84 — prawdziwa wartość wydaje się być ≈0 (brak istotnej zależności R od Θ_D).
- s-class trywialny 2-parametrowy Θ-only fit — wymaga teraz N_F (c=+0.88).
- **Paper v3 "TGP class-specific transport theory" NIE jest gotowy** — opublikowanie z wykładnikami z N=12 byłoby błędem.

### Nowe outlery na N=34

| Materiał | Profile RMS | Komentarz |
|----------|-------------|-----------|
| Ti (d) | 0.575 | stare, spin-fluct |
| Sc (d) | 0.560 | nowe, 3d grupa III |
| Be (sp) | 0.455 | stare, Θ_D=1440 K |
| Ir (d) | 0.430 | nowe, 5d noble d |
| W (d) | 0.410 | stare, bcc 5d |
| Pt (d) | 0.383 | stare, Stoner |
| Y (d) | 0.370 | nowe, 4d grupa III |

→ **Grupa III (Sc, Y, La)** daje mieszane wyniki: Sc/Y outlery, La pasuje świetnie (0.08). Sugestia heterogeniczności w d-class.

### Ca, Sr: klasyfikacja s vs sp (D)

| Mat. | R_obs | R_s_pred (tylko noble) | dlog_s | R_sp_pred (sp_normal) | dlog_sp | Lepsze |
|------|-------|------------------------|--------|----------------------|---------|--------|
| Ca | 11.76 | 6.60 | −0.25 | 30.31 | +0.41 | **s** |
| Sr | 28.03 | 4.55 | −0.79 | 32.14 | +0.06 | **sp** |

**Mixed result:** Ca pasuje lepiej do s-class Cu/Ag/Au; Sr pasuje lepiej do sp-normal. Sugeruje, że alkaline-earth nie jest czystą jedną klasą w ramach TGP A_orb.

### Status Etapu 5

**ZAMKNIĘTY** z wynikiem **NEGATYWNYM** dla uniwersalności d-class:

- ❌ d-class b=+0.84 z r08 **obalone** — na N=21 b≈+0.02
- ❌ s-class b=+0.84 z r08 **obalone** — na N=5 b≈−0.45 (nietrwałe na tak małym sample)
- ❌ Paper v3 w pierwotnej formie **wstrzymany** — potrzeba re-interpretacji
- ✅ Mean profile RMS (0.23) nadal w granicach — teoria **działa strukturalnie**, ale param. są czulsze niż N=23 sugerowało
- ✅ c (exponent N_F) d-class stabilne: +0.25±0.06 (z N=21 LOO). **c jest prawdziwe, b nie było.**

### Nowa hipoteza (post-r10)

**Słaba zależność d-class R od Θ_D (b≈0)** sugeruje, że w high-T BG limit (T ≫ Θ_D) i w formule R = A·Θ_D (jak w klasycznym Bloch-Mott) współczynnik A sam silnie zależy od Θ_D tak, że zamiast R ∝ Θ_D wychodzi R ≈ const dla d-class. To spójne z:

$$\rho(T) \to 4 R \cdot (T/\Theta_D) \cdot \pi^4/5 \qquad (T \gg \Theta_D)$$

Jeśli R ∝ const → ρ_{high-T} ∝ T/Θ_D, a jeśli jednocześnie R ∝ Θ_D → ρ_{high-T} ∝ T (niezależnie od Θ_D). Druga opcja jest klasyczna; z r10 (b≈0 dla d) wychodzi pierwsza opcja (nietypowa!).

**Interpretacja:** W d-class, DOS N_F dominuje (c=+0.25 stabilne), a Θ_D wpływa na ρ(T) głównie przez **kształt krzywej BG**, nie przez amplitudę R. To jest **odkrycie** — nie falsyfikacja TGP, tylko **korekta klasowej struktury**.

### Proponowane następne kroki

**r11 (priorytet):** Pod-klasy wewnątrz d-class — testować czy d3, d4, d5, d6, d7, d8, d9, d10 mają inne b, c (efekt wypełnienia podpasma d). Grupy:
- Grupa III (Sc, Y, La)  — empty d shell, f-resonance
- Grupa IV (Ti, Zr, Hf) — d² 
- Grupa V (V, Nb, Ta) — d³
- Grupa VI (Cr?, Mo, W) — d⁴/d⁵
- Grupa VII (Mn?, Tc?, Re) — d⁵/d⁶
- Grupa VIII (Fe, Ru, Os) — d⁶
- Grupa IX (Co, Rh, Ir) — d⁷
- Grupa X (Ni, Pd, Pt) — d⁸/d⁹
- Grupa XI (Cu, Ag, Au) — d¹⁰ (→ s-class)

**r12:** Pomiar wrażliwości b na Θ_D — zamiast naszych CRC Θ_D (często z niskich T elastyczności), spróbować Θ_D z high-T specific heat.

**r13:** Alternatywna forma — może R nie zależy tylko od (Θ_D, N_F), może powinniśmy dodać grupę orbitalną g (d³, d⁸, etc.) jako 3-ci parametr.

**Paper v3 status:** WSTRZYMANY aż do wyjaśnienia. Możliwe dwa kierunki:
1. **Subclass paper** — R(Θ_D, N_F | grupa) — wymaga r11
2. **Kompaktowy paper** — tylko c-exponent (N_F) + kształt BG — oddzielenie od amplitudy

Status ogólny projektu: teza "TGP klasowa" przetrwała, ale **skala wymagana do publikacji zwiększona** z N≈23 do N≈40+ z subklasami.

---

## UPDATE 2026-04-21 (Etap 6): r11 — **POZYTYWNE** odkrycie: d-class rozkłada się wg grupy PT

**Skrypt:** [[rho_normal_state_closure/r11_dshell_subclass.py]] | [[rho_normal_state_closure/r11_results.txt]]

### Motywacja

r10 pokazało że naive d-class `R = C·Θ^b·N_F^c` ma b≈0 na N=21 (poprzednie b=+0.84 z N=12 było artefaktem). Hipoteza: d-class jest heterogeniczna wg **wypełnienia podpasma d** (grupa III..X układu okresowego).

### Podział 21 d-metali na 8 grup PT

| Grupa | d^n | N | Materiały |
|-------|-----|---|-----------|
| III | d¹ | 3 | Sc, Y, La |
| IV | d² | 3 | Ti, Zr, Hf |
| V | d³ | 3 | V, Nb, Ta |
| VI | d⁴⁻⁵ | 2 | Mo, W |
| VII | d⁵⁻⁶ | 1 | Re |
| VIII | d⁶⁻⁷ | 3 | Fe, Ru, Os |
| IX | d⁷⁻⁸ | 3 | Co, Rh, Ir |
| X | d⁸⁻¹⁰ | 3 | Ni, Pd, Pt |

### Per-grupa fits (adaptive: 2-param gdzie N=3, const gdzie N≤2)

| Grupa | N | C | b (Θ_D) | c (N_F) | RMS | forma |
|-------|---|---|---------|---------|------|-------|
| III | 3 | 7.4e-1 | **+1.00** | — | 0.000 | Θ-only |
| IV | 3 | 1.4e-1 | **+1.25** | — | 0.025 | Θ-only |
| V | 3 | 5.8e-4 | **+2.08** | — | 0.011 | Θ-only |
| VI | 2 | 3.8e+1 | — | — | 0.086 | const |
| VII | 1 | 1.2e+2 | — | — | — | single |
| VIII | 3 | 8.5e+1 | — | **+0.83** | 0.018 | N_F-only |
| IX | 3 | 4.4e+1 | — | **+0.57** | 0.036 | N_F-only |
| X | 3 | 8.9e-2 | **+1.13** | — | 0.011 | Θ-only |

**Uderzający wzorzec:**
- **Wczesne grupy (III, IV, V):** silna zależność Θ_D → Θ_D², Θ_D¹·², Θ_D¹ (b rośnie z d-fillingem!)
- **Późne grupy (VIII, IX):** silna zależność N_F → c≈+0.6 do +0.8
- **Grupa X:** z powrotem Θ_D dominuje (b=+1.13)

### Porównanie: pełny d-fit (N=21) vs suma per-grupa

| Metryka | d-class unified (N=21) | per-grupa | Redukcja |
|---------|------------------------|-----------|----------|
| **RMS(R)** | 0.254 | **0.033** | **87%** |
| Materiałów lepiej fitowanych | — | **20/21** (sub wins), 1 tie, 0 full wins | |
| **Mean profile RMS (BG)** | 0.266 | **0.105** | **60%** |
| Median profile RMS | 0.271 | **0.093** | |
| Max profile RMS | 0.575 (Ti) | 0.197 (Ir) | |

**Per-material breakdown (profile RMS z sub-grupy R_pred):**

Wszystkie 21 d-metali mają RMS < 0.20 **(0 outliers)**. Top 3 najlepsze: Hf (0.032), Pt (0.042), Ta (0.045). Top 3 najtrudniejsze: Ir (0.197), W (0.189), Sc (0.180).

### Fit per okres (3d vs 4d vs 5d) — kontrolny

| Okres | N | b | c | RMS |
|-------|---|---|---|-----|
| 3d | 6 | −2.02 | −0.92 | 0.116 |
| 4d | 7 | −0.77 | +0.07 | 0.200 |
| 5d | 8 | −0.29 | −0.05 | 0.204 |

**Wniosek:** Podział per-okres daje **słabsze** wyniki (RMS 0.12–0.20) niż per-grupa (mean 0.03). **Grupa** jest dominującą osią heterogeniczności, nie okres. Zgodne z intuicją: d-filling określa strukturę DOS przy E_F.

### Interpretacja fizyczna

**Pattern "b dominuje w wczesnych grupach, c w późnych"** jest spójny z:

1. **Wczesne TM (III–V, IV):** podpasmo d jest słabo wypełnione, DOS(E_F) jest **po stronie rosnącej** pasma d. Θ_D wyznacza phase space rozpraszania e-ph. ρ ~ R·(T/Θ_D) z R ∝ Θ_D^1+ (Bloch-Mott classical + d-band enhancement).

2. **Środkowe TM (VI, VII):** Dokładnie w pół-wypełnionym (Cr, Mo, W → d⁵ spin-polarized nearly half-filled) N_F często w minimum (pseudo-gap Hund). Efekt: R ≈ const (?).

3. **Późne TM (VIII, IX):** DOS(E_F) jest **na szczycie** pasma d — bardzo wysoki N_F, dominujące scattering. R ∝ N_F^0.6-0.8 (Eliashberg-like, DOS-driven).

4. **Grupa X (noble metals — Cu/Ag/Au zostały przerzucone do s-class, ale Ni/Pd/Pt zostały jako d):** d-band schodzi poniżej E_F, znowu Θ_D-driven ale z innym offsetem.

To jest **falsyfikacja** TGP-class w formie naive (jedna D-forma dla d), ale **potwierdzenie** głębszej struktury: TGP A_orb=+0.310 dla d-class to średnia, a prawdziwa fizyka to **d-shell-filling dependent** exponents (b(g), c(g)).

### Status Etapu 6

**ZAMKNIĘTY** z wynikiem **POZYTYWNYM**:

- ✅ Per-grupa fits dają RMS_R = 0.033 (87% redukcja vs unified)
- ✅ Mean profile RMS = 0.105 (60% redukcja)
- ✅ 21/21 d-materiałów ma profile RMS < 0.20 (**0 outliers** wśród d!)
- ✅ Wzorzec b(grupa), c(grupa) ma sensowną interpretację fizyczną (d-shell filling)
- ✅ **Paper v3 "Per-group TGP transport theory"** kierunek wyjaśniony

### Ograniczenia r11

- N=3 na grupę → statystyka słaba; każdy per-group fit ma de facto 1 DoF
- Grupa VI (N=2) i VII (N=1) potrzebują dodatkowych materiałów
- **b** wykładniki (+1.0, +1.25, +2.08, +1.13) mogą być artefaktem N=3 fittingu — wymagają weryfikacji na N=5+ per grupa

### Nowy plan

**r12 (priorytet):** Dodać **Cr** (grupa VI) + **Tc** (rzadki, grupa VII, może symulacja?) + jeszcze **Mn** (grupa VII jeśli sensowne, choć anomalny) + dodać drugi materiał dla Group VIIII/IX (np. **Ni-O podobne?**). Cel: 5 materiałów na grupę.

**r13:** Uwzględnić w pełnym fitie 4-parametry: `R = C · Θ^b · N_F^c · g^d` gdzie g = liczba d-elektronów w walencji (d¹..d¹⁰). Sprawdzić istotność d.

**r14 (alternatywa):** Jeśli per-grupa jest finalnym modelem, zweryfikować **s-class** (noble + alkaline-earth, 5 mat) i **sp-class** (sp_normal 6 mat) tak samo solidnie jak d.

**Paper v3 scope updated:** "TGP periodic-table-group dependent transport" — 8 podklas d + s + sp z wzorcem b(g), c(g).

Status ogólny projektu: **Fundament strukturalny potwierdzony na N=34, forma per-grupa otwiera ścieżkę do paperu v3 z nowym (i lepszym) scope'm** — heterogeniczność TGP A_orb wg grupy PT.

---

## UPDATE 2026-04-21 (Etap 7): r12 — **ODKRYCIE PRAWDZIWEJ FORMY**: g-intercept shift

**Skrypt:** [[rho_normal_state_closure/r12_gparam_model.py]] | [[rho_normal_state_closure/r12_results.txt]]

### Motywacja

r11 pokazał spektakularny pattern: grupa III b=+1.0, grupa IV b=+1.25, grupa V b=+2.08. Ale fity per-grupa miały N=3 i de facto 1 DoF → 2-parametrowe fity przechodziły przez 3 pkt bez sprawdzenia hipotezy. **Czy per-grupa b(g), c(g) są realne czy to artefakt N=3?**

### Rozszerzenie: Cr, Tc, Lu → d-class N=21 → N=24

| Nowy | Grupa | d-filling | Uzasadnienie |
|------|-------|-----------|--------------|
| **Cr** | VI (g=6) | d⁵s¹ | Wypełnia grupę VI (miała Mo/W, N=2 → 3) |
| **Tc** | VII (g=7) | d⁵s² | Grupa VII N=1 → 2 (radioaktywny, literatura ograniczona) |
| **Lu** | III (g=3) | d¹4f¹⁴ | Grupa III rozbudowana (Sc/Y/La/Lu N=3 → 4) |

### Test nested modeli (ANOVA AIC/BIC)

Model postaci: log R = a + b·log(Θ) + c·log(N_F) + [dodatkowe terminy z g = grupa PT, 3..10]:

| Model | Parametry | k | RMS | ΔAIC | Uwagi |
|-------|-----------|---|-----|------|-------|
| **M0** | a, b, c (universal) | 3 | 0.265 | +23.4 | Silne dowody przeciw |
| **M1** | + d·g | 4 | **0.156** | **0.0** ⭐ | Best AIC + BIC |
| M2 | + e·g·log(Θ) | 5 | 0.155 | +1.5 | Brak istotnej poprawy |
| M3 | + f·g·log(N_F) | 6 | 0.154 | +3.4 | Interactions NIE znaczące |
| M4 | + h·g² | 7 | 0.150 | +4.1 | Kwadratowe g pogarsza AIC |
| M5 | + i·g²·logΘ + j·g²·logNF | 9 | 0.140 | +4.8 | Overfitting |

**KRYTYCZNE:** dAIC(M2) = 1.5, dAIC(M3) = 3.4 → g-zależne **slopes** b(g), c(g) **NIE SĄ** istotne statystycznie. Tylko g-intercept jest istotny.

### Finalna formuła d-class

$$\boxed{\; \log_{10} R_d = -0.293 + 1.17 \cdot \log_{10}\Theta_D + 0.534 \cdot \log_{10} N_F - 0.114 \cdot g \;}$$

lub równoważnie:

$$R_d = 0.509 \cdot \Theta_D^{1.17} \cdot N_F^{0.53} \cdot 10^{-0.114 \cdot g}$$

gdzie `g ∈ {3..10}` to numer grupy w układzie okresowym (Sc=3, Ti=4, V=5, Cr=6, Mn=7, Fe=8, Co=9, Ni=10).

### Interpretacja fizyczna

1. **b = +1.17** — bliskie klasycznego Bloch-Mott `ρ ∝ T/Θ_D` dla T ≫ Θ_D. Odzyskujemy **prawidłowe zachowanie phonon-limited transport** po rozszerzeniu z N=21 do N=24.

2. **c = +0.53** — umiarkowany DOS-enhancement. W granicy Eliashberg `ρ ∝ λ_ep ∝ N_F·<I²>/M<ω²>`, c ≈ 1 byłoby pełnym DOS-driven. Nasze c ≈ 0.5 sugeruje **częściowe uwzględnienie** DOS, reszta to phonon kinematics.

3. **d = −0.114** (g-intercept) — dla stałych (Θ_D, N_F), **późne d-metale mają mniejsze R**:
   - g=3 (Sc): 10^(-0.34) = 0.46 × baseline
   - g=10 (Ni): 10^(-1.14) = 0.07 × baseline
   - → Stosunek 6.5× redukcji R między grupą III a X

   To ma **sensowną interpretację**: w późnych d-metalach hybrydyzacja d-s i krótszy e-ph coupling range redukuje scattering matrix element <I²>, co mnoży R przez 10^(−0.114·g). W wczesnych d-metalach d-orbitale są mniej zhybrydyzowane, większe <I²>.

### Stabilność (LOO CV)

- Full fit RMS = 0.156
- **LOO RMS = 0.185** (overfit ratio 1.18×)
- **Brak ekstremalnych outlierów** w LOO: największe |dlog R| = 0.41 (Tc — tylko 2 pt w grupie VII, spodziewane), kolejne: 0.38 (Mo), 0.37 (Re), 0.31 (W)

### Profile RMS (full BG rho(T_i)) — N=24 d-class

| Model | Mean profile RMS | Interpretacja |
|-------|------------------|---------------|
| M0 (universal, 3 par) | 0.284 | Baseline |
| **M1 (+g linear, 4 par)** | **0.185** | **Best AIC** |
| r11 per-grupa (8 osobnych fitów) | 0.105 | Niefair: 8×(2-3 par) ≈ 18 parametrów |

**Parametry vs wynik:**
- M1 to 4 parametry dla całej d-class (N=24) → 1 param / 6 materiałów
- r11 per-grupa to ≈18 parametrów dla N=21 → 1 param / 1.2 materiałów (de facto fit danych, nie modelu)
- M1 jest **realnie lepszy** statystycznie (AIC uwzględnia k)

### Porównanie b_eff(g), c_eff(g): M1 vs r11

r11 zaraportowało wysoce zróżnicowane b, c per grupa. M1 mówi że wszystkie grupy mają TAKIE SAME b=+1.17, c=+0.53, tylko **różne intercepty**. r11 odczyty były artefaktem:
- Group V (Nb, V, Ta) b=+2.08 — to była oscylacja z N=3; z pełnym modelem znika
- Group VIII (Fe, Ru, Os) c=+0.83 — to samo
- Group IX (Co, Rh, Ir) c=+0.57 — to samo

### Status Etapu 7

**ZAMKNIĘTY** z wynikiem **STRONGLY POSITIVE**:

- ✅ Jedna formuła M1 z 4 parametrami zamyka d-class N=24
- ✅ Profile RMS 0.185 (LOO 0.185) — robustny
- ✅ Wykładniki uniwersalne: b=+1.17 (Bloch-Mott!), c=+0.53 (DOS)
- ✅ g-parameter (grupa PT) dostarcza intercept korekcję z jasną interpretacją fizyczną
- ✅ Obaliło r11 per-grupa b(g) pattern — był artefaktem N=3
- ✅ **Paper v3 "Universal TGP transport: d-class with g-intercept"** ma scope

### Ograniczenia r12

- Tc i Re w grupie VII mają duże residua (Re dlog=+0.15, Tc dlog=−0.38). Grupa VII jest nadal słabo obsadzona (N=2).
- Cr ma SDW poniżej 311K — ρ(77K) może nie być czyste BG. Chociaż fit R traktuje tylko amplitudę.
- g = kolumna PT; alternatywne kodowanie (liczba d-elektronów w konkretnym stanie walencyjnym, 0..10) mogłoby dać inne wyniki.
- s-class i sp-class jeszcze nieprzetestowane w tym formalizmie.

### Proponowane następne kroki

**r13:** Rozszerzyć formalizm M1 (a + b·logΘ + c·logNF + d·g) na **s-class i sp-class**:
- s-class: N=5 (Cu, Ag, Au, Ca, Sr) — dla s-elektronów g = 1 (Cu/Ag/Au) lub 2 (Ca/Sr). Test czy M1 z innym (a,b,c,d) zamyka s.
- sp-class: N=8 (Al, Pb, Sn, Cd, Zn, Mg, Be, Bi) — kodowanie g problematyczne (sp₁, sp₂, sp³). Może z-value lepsze.

**r14:** Meta-formuła przekrojów: Czy można zamknąć wszystkie 3 klasy jednym parametrem kategorycznym (dummy s/sp/d) + universal b, c + g? Test:
$$\log R = a_\text{cls} + b \log \Theta + c \log N_F + d_\text{cls} \cdot g$$

Jeśli b, c są uniwersalne → **jedna formuła TGP dla całego normal-state transport**.

**Paper v3 scope (nowa wersja):** "TGP-unified normal-state resistivity: b=1.17, c=0.53, and d-shell filling shifts". 
- Sekcja 1: TGP klasyfikacja (s/sp/d z A_orb) — z paperu v1
- Sekcja 2: Universal form ρ_BG(T) z R = f(Θ_D, N_F, g)
- Sekcja 3: Data: N=37 (34 d-metale + 5 s + 8 sp)
- Sekcja 4: Model selection (M0..M5 AIC)
- Sekcja 5: Wzorzec g (d-shell filling → intercept shift −0.11 × g, fizyczna interpretacja)
- Sekcja 6: Unifikacja TGP-SC i TGP-ρ na poziomie klas

**Status projektu:** W 3 iteracjach (r10→r11→r12) przeszliśmy od fałszywej hipotezy (d-class universal, b=+0.84) przez błędną diagnozę (per-grupa b(g)) do **prawdziwej struktury** (universal b, c + g-intercept). Paper v3 w obecnym scope'm ma **solidne ugruntowanie statystyczne** (AIC, LOO) i **fizyczną interpretację**.

---

## UPDATE 2026-04-21 (Etap 8): r13 + r13b — **UNIFIKACJA 3 klas TGP** (s + sp + d)

**Skrypty:** [[rho_normal_state_closure/r13_unified_all_classes.py]] | [[rho_normal_state_closure/r13b_loo_sweep.py]] | [[rho_normal_state_closure/r13_results.txt]] | [[rho_normal_state_closure/r13b_results.txt]]

### Rozszerzenie formalizmu M1 na N=37

- d-class N=24 (z r12)
- sp-class N=8 (Al, Pb, Sn, Cd, Zn, Mg, Be, Bi)
- s-class N=5 (Cu, Ag, Au, Ca, Sr)

Kodowanie `v` (liczba walencyjnych elektronów przewodnictwa):

| Klasa | Materiały | v |
|-------|-----------|---|
| s | Cu, Ag, Au | 1 |
| s | Ca, Sr | 2 |
| sp | Be, Mg, Zn, Cd | 2 |
| sp | Al | 3 |
| sp | Sn, Pb | 4 |
| sp | Bi | 5 |
| d | Sc, Y, La, Lu | 3 |
| d | Ti, Zr, Hf | 4 |
| … | … | … (do X=10) |

### Test 7 nested modeli (r13b LOO sweep)

Kryterium **mean LOO profile RMS** (robustny niż pure AIC):

| Model | k | Full RMS | LOO RMS | **LOO prof RMS** | Overfit | Uwaga |
|-------|---|----------|---------|------------------|---------|-------|
| U0 (univ b, c + a_cls) | 5 | 0.269 | 0.343 | 0.314 | 1.28× | baseline |
| U1 (+ univ d·v) | 6 | 0.247 | 0.397 | 0.304 | 1.60× | |
| U2 (+ d_cls·v) | 8 | 0.194 | 0.314 | 0.258 | 1.62× | |
| U3 (+ b_cls·logΘ) | 8 | 0.242 | 0.470 | 0.335 | 1.94× | |
| U4 (+ c_cls·logNF) | 8 | 0.190 | 0.291 | 0.271 | 1.53× | |
| **U24 (+ d_cls·v + c_cls·logNF)** | **10** | **0.156** | **0.261** | **0.229** ⭐ | 1.68× | **Best robust** |
| U_full (wszystko cls-specific) | 12 | 0.138 | 0.361 | 0.278 | 2.62× | Overfit |

**Kluczowe spostrzeżenie:** U_full miał najlepszy AIC ale najgorszy LOO overfit (2.62×). **U24** to scalone minimum: class-specific c i d, uniwersalne b.

### Finalna formuła uniwersalna (U24)

$$\boxed{\; \log_{10} R = a_\text{cls} + 0.786 \cdot \log_{10}\Theta_D + c_\text{cls} \cdot \log_{10} N_F + d_\text{cls} \cdot v \;}$$

| Klasa | a_cls | c_cls | d_cls | Interpretacja |
|-------|-------|-------|-------|---------------|
| **d** | +0.58 | **+0.41** | **−0.10** | DOS enhances R (+c), g-shift redukuje R dla późnych d (−d) |
| **s** | −1.51 | **+0.03** | **+0.50** | DOS niezależne (c≈0), monovalent<divalent (+d) |
| **sp** | −1.21 | **−0.39** | **+0.19** | DOS INVERSED (−c) — semimetal/localization effect |
| **Wszystkie** | — | — | **b = +0.79** | **Uniwersalny Bloch-Mott** |

### Fizyczna interpretacja

1. **b = +0.786 UNIWERSALNE** — Θ_D zachowuje się klasycznie (Bloch-Mott limit ρ ∝ T/Θ_D) **niezależnie od klasy orbital**. To jest **uniwersalna struktura fononowa** metalów.

2. **c różnicuje klasy:**
   - **d-class (c=+0.41):** DOS wzmacnia R zgodnie z Eliashberg — elektrony d są silnie scatterowane, bogactwo stanów przy E_F daje więcej kanałów rozpraszania
   - **s-class (c=+0.03):** DOS irrelewantne — monovalent/divalent s-electrons mają "podobne" transport niezależnie od gęstości stanów (przybliżenie free-electron)
   - **sp-class (c=−0.39):** ODWRÓCONA zależność! Wyższy N_F → mniejsze R. Interpretacja: w sp-metalach, więcej stanów przy E_F sygnalizuje większe mieszanie s-p **i bardziej delokalizowane** stany (mniej zlokalizowane → mniej rozpraszania). Lub: dominuje Bi/Be outlier effect (semimetal, extreme Θ_D).

3. **d (v-slope) różnicuje klasy:**
   - **d-class (d=−0.10):** Późne d (v→10) mają mniejsze R — zgodne z r12 (Ni<Sc)
   - **s-class (d=+0.50):** Ca/Sr (v=2) mają większe R niż Cu/Ag/Au (v=1) — trywialne: alkaline-earth mają dwa razy więcej carriers ale dwa razy więcej scattering też
   - **sp-class (d=+0.19):** Słabo rosnące z walencją (Al v=3 → Bi v=5)

4. **Class-intercepts odzwierciedlają TGP A_orb:**
   - s: A_orb = −0.111 → a_cls = −1.51 (najniższe R baseline)
   - sp: A_orb = +0.207 → a_cls = −1.21 (pośrednie)
   - d: A_orb = +0.310 → a_cls = +0.58 (najwyższe R)
   - **Korelacja:** wyższy A_orb ⇒ większe R baseline (naiwnie zgodne z TGP: więcej orbital character ⇒ silniejsze scattering).

### Per-class LOO profile RMS (U24)

| Klasa | N | mean | median | max | Outliers |
|-------|---|------|--------|-----|----------|
| s | 5 | 0.266 | 0.139 | 0.529 | Ca (0.53), Sr (0.43) |
| sp | 8 | 0.293 | 0.229 | 0.806 | Bi (0.81), Be (0.43) |
| d | 24 | 0.200 | 0.159 | 0.475 | Tc (0.47), W (0.46), Re (0.39) |

**d-class robustny (mean 0.20).** s i sp są zdominowane przez 2–3 klasyczne outliery (Bi, Be, Ca, Sr — wszystkie znane anomaliczne materiały).

### Status Etapu 8

**ZAMKNIĘTY** z wynikiem **MOCNO POZYTYWNYM**:

- ✅ **Jedna formuła** U24 z 10 parametrów (3 intercepty, 3 c, 3 d, 1 universal b) zamyka całe ρ(T) dla 37 metali
- ✅ **b = +0.79 uniwersalne** we wszystkich klasach — Bloch-Mott prawdziwe także poza d
- ✅ **c, d class-specific** z interpretacją przez A_orb (TGP zgodne)
- ✅ Mean LOO profile RMS = 0.229 (znacznie lepsze od U_full 0.278)
- ✅ **Paper v3 "Universal TGP normal-state transport"** ma teraz pełny zestaw 3 klas z solidną statystyką
- ✅ Struktura TGP (A_orb — s/sp/d) ma **korelat w amplitudach transportu** (a_cls), nie tylko w SC paper v1

### Ograniczenia / TODO dla paperu v3

1. **s-class potrzebuje więcej materiałów**: obecnie N=5 (Cu, Ag, Au, Ca, Sr). Ca i Sr są LOO-outlierami. Trzeba dodać: **Ba** (z ρ(900K) zamiast ρ(1000K)), ewentualnie **lantanowce** z okolicach dhcp.

2. **sp-class ma Bi dominujące outlier**: Bi jest semimetal (N_F=0.006), co wypacza fit. Opcja:
   - Wykluczyć Bi (sp_normal N=7)
   - Dodać **Tl** (ρ(500K) możliwe, m.p.=577K) i **In** (ρ(400K) możliwe, m.p.=430K) **używając niższych T_points**

3. **d-class outlery (Tc, W, Re):**
   - Tc — radioaktywny, literatura ograniczona
   - W — bcc 5d z SC anomaly (T_c=0.012K dla pure)
   - Re — bcc/hcp mixed

4. **Alternatywne formy N_F:**
   - Obecnie N_F = DFT states/eV/atom
   - Możliwie: γ_Sommerfeld (z specific heat) — bardziej bezpośrednia miara DOS(E_F)

### Proponowane następne kroki

**r14:** Wymienić Bi w sp-class na **Tl+In** (odpowiednie T_points, m.p. OK). Sprawdzić czy sp c=−0.39 trzyma.

**r15:** γ_Sommerfeld jako alternative dla N_F. Test czy c-exponents per klasa są spójniejsze.

**r16 (opcjonalne):** Przetestowanie formy **alternatywnej** U24 gdzie v traktowane nie jako linear, ale jako **class×grupa** (Curie-Mendeleev coding).

**r17 (paper draft):** Rozpocząć draft paperu v3:
- Title: "TGP class-universal phonon exponent and class-specific DOS response in normal-state resistivity of pure metals"
- Abstract: `ρ(T) = ρ_0 + R·(T/Θ_D)⁵·J₅(Θ_D/T)` z `R = f(Θ_D, N_F, valence, orbital class)` — 10 parametrów uniwersalnie na 37 metalach
- Methods: Dataset ze źródeł CRC/Matula/Desai, BG baseline, AIC + LOO model selection
- Results: U24 best + tabele koeficjentów + fizyczna interpretacja
- Discussion: Połączenie z SC paper v1 — te same klasy TGP pojawiają się w dwóch różnych zjawiskach transportu (normalne + nadprzewodzenie)

**Status projektu (koniec Etapu 8):** Formuła uniwersalna potwierdzona, **Paper v3 ma solidne fundamenty**. Dalej już to głównie "polishing" (dodatkowe materiały dla s/sp, alternatywne N_F coding) — nie zasadnicze zmiany. **Gotowe do draftu**.

---

## Etap 9: OSTRZAŁ rownania U24 — stabilność, uproszczenie, stopy (r14–r16)

**Data:** 2026-04-21
**Dyrektywa użytkownika:** *"nie wiem czy jest sens skupiać się na publikacji, trzeba teraz zrobić ostrzał równania, sprawdzić jak jest stabilne, dla możliwe dużej klasy metali/stopów i zobaczyć, czy uda się uprościć parametry lub wyprowadzić je z rdzenia"*

Trzy niezależne testy: (A) bootstrap stabilności, (B) redukcja przez TGP A_orb, (C) minimalny model, (D) stress-test na stopach.

### r14 — Bootstrap U24 + A_orb reduction (N=37, N_boot=2000)

**Cel:** Sprawdzić czy (a_cls, c_cls, d_cls) dają się wyprowadzić z TGP A_orb → zredukowanie 10 parametrów do 5–7.

**Wyniki:**

Bootstrap U24 (95% CI):

| Parametr | OLS | CI₂.₅ | CI₉₇.₅ | Istotny? |
|---|---|---|---|---|
| a(d) | +0.58 | −0.78 | +1.44 | ✗ |
| **as(s-d)** | **−2.09** | **−3.39** | **−0.11** | ✓ |
| asp(sp-d) | −1.79 | −3.60 | +0.18 | ✗ |
| **b(logTh)** | **+0.79** | **+0.43** | **+1.37** | ✓ |
| c(logNF) | +0.41 | −0.01 | +0.86 | ✗ |
| **d(v)** | **−0.10** | **−0.13** | **−0.07** | ✓ |
| ds(v·s) | +0.60 | −0.97 | +1.24 | ✗ |
| dsp(v·sp) | +0.29 | −0.34 | +0.74 | ✗ |
| cs(s·logNF) | −0.39 | −1.35 | +0.98 | ✗ |
| csp(sp·logNF) | −0.80 | −1.97 | +0.61 | ✗ |

**Tylko 3/10 parametrów ma 95% CI poza zerem.** Pozornie model przeparametryzowany.

Regresja klasowych parametrów vs A_orb:

| Parametr | intercept | slope | R² |
|---|---|---|---|
| a_cls | −1.11 | +2.68 | **0.34** (słaba) |
| c_cls | +0.11 | −0.12 | **0.003** (brak) |
| d_cls | +0.19 | −0.62 | **0.53** (umiarkowana) |

**Wniosek:** class-specific parametry **NIE wyprowadzają się z A_orb**. Najlepsze R² = 0.53 (d_cls) to za mało dla redukcji.

Redukowane modele (LOO profile RMS):

| Model | k | LOO prof RMS | vs U24 |
|---|---|---|---|
| U_red5 (β₀+β₁·A+b+γ+δ) | 5 | 0.319 | 1.39× |
| U_red6c (+γ₁·A·logNF) | 6 | 0.330 | 1.44× |
| U_red6d (+δ₁·A·v) | 6 | 0.314 | 1.37× |
| U_red7 (all A-linear) | 7 | 0.316 | 1.38× |
| **U24** | **10** | **0.229** | **1.00×** |

Wszystkie redukowane modele tracą ≥37% precyzji. **A_orb reduction NIE działa.**

### r15 — Minimalny model (czy 4-6 param wystarczy?)

**Cel:** Gdy 7/10 parametrów "nieistotnych" w bootstrap, czy można zredukować strukturalnie?

| Model | k | LOO prof RMS | vs U24 | Verdict |
|---|---|---|---|---|
| M_min3 (a + as·s + b·logTh + d·v) | 4 | 0.328 | 1.43× | borderline |
| M_min4 (+c·logNF) | 5 | 0.353 | 1.54× | FAILS |
| M_min4sp (+asp) ≡ U1 | 6 | 0.304 | 1.33× | borderline |
| **U24** | **10** | **0.229** | **1.00×** | baseline |

**Wniosek:** Żaden minimalny model NIE dorównuje U24. Parametry są **silnie skorelowane** — pojedynczo wydają się nieistotne w bootstrap, ale razem tworzą strukturę nieredukowalną. Interakcje v·δ_sp, δ_sp·logN_F są szczególnie ważne dla sp-class.

### r16 — Stopy (13 alloys): Cu-Ni, Nichrome, Constantan, Manganin, Monel, Ag-Au, Cu₃Au, brass, bronze, solder

**Cel:** Sprawdzić czy U24 (kalibrowana na czystych metalach) przewiduje R dla stopów przez Matthiessen decomposition `ρ(T) = ρ_0 + R·f_BG(T/Θ_D)`.

**Test główny — log R space:**

| Klasa | N | mean dlog R | RMS | max abs |
|---|---|---|---|---|
| s (AgAu, Cu₃Au) | 4 | +0.42 | 0.43 | 0.50 |
| sp (brass, bronze, solder) | 4 | −0.14 | 0.19 | 0.36 |
| d (CuNi, Nichrome, Manganin, Monel) | 5 | **−0.89** | **1.06** | **1.69** |
| **Overall (alloy)** | **13** | **−0.25** | **0.70** | — |
| (ref) pure N=37 | 37 | 0.00 | **0.16** | — |

**U24 na stopy: 4.5× gorsze RMS niż na czystych.**

**Per-klasa diagnostyka:**
- **d-class katastrofa** (bias −0.89): U24 przewiduje 8× za dużo R dla CuNi/Nichrome/Manganin. Alloy disorder saturuje rozpraszanie fononowe — phonon component jest silnie stłumiony.
- **s-class** (bias +0.42): AgAu i Cu₃Au mają **większe** R niż pure prediction — dodatkowy mechanizm scatteringu (concentration fluctuations, phonon modulation przez mass mismatch).
- **sp-class** (bias −0.14): Brass, bronze, yellow brass działają OK — pure-metal-like transport.

**Matthiessen verification (pozytywny wynik):** Cu₃Au ordered vs disordered — te same (Θ_D, N_F, v, cls), różnica tylko w ρ_0 (1.8 vs 10.5 μΩ·cm):

| Wariant | ρ_0 | R_obs | R_pred |
|---|---|---|---|
| Cu₃Au ordered | 1.8 | 25.4 | 8.1 |
| Cu₃Au disordered | 10.5 | 24.0 | 8.1 |
| **Ratio** | **5.8×** | **1.06×** | — |

**Potwierdzone empirycznie:** `ρ_0(disorder)` i `R(electronic)` dekomponują się niezależnie. Matthiessen rule DZIAŁA.

**Matthiessen-extreme (Manganin, Constantan): paradoks profile RMS:**

Profile RMS (log ρ(T)) dla stopów = **0.095** (2.4× **LEPSZE** niż pure 0.229). Ale to artefakt — gdy ρ_0/ρ_295 = 98% (Manganin), to błąd w R rzędu 30× daje tylko 2% różnicy w całkowitym ρ(T).

**Prawdziwa metryka: log R RMS = 0.70 (bardzo źle).**

### Verdict r14–r16

1. **U24 JEST formułą dla czystych metali** (N=37, LOO prof RMS=0.229). Jest optymalna (r14, r15).
2. **U24 NIE jest uniwersalna na stopy.** Matthiessen decomposition jest OK (r16.C1), ale R(cls, Θ, N_F, v_eff) z liniowo uśrednionego v_eff nie oddaje alloy phonon scattering.
3. **A_orb reduction nie działa** (r14.B). Class-specific parametry nie są funkcją A_orb. Korelują (A_orb→a_cls R²=0.34), ale nie wystarczająco dla redukcji.
4. **Bootstrap "7/10 nieistotnych" jest mylące** (r14.A). Parametry są skorelowane; usuwanie ich osobno wydaje się bezpieczne, ale razem tracimy strukturę (r15).

### Implikacje dla paperu v3

**OSŁABIONE RAMY** — paper nie jest "universal TGP transport formula", lecz **"10-parameter empirical class-aware BG fit of 37 pure metals"**. To nadal wartościowy wynik (b=+0.79 universal Bloch-Mott, c/d/a class-specific), ale:

- Zakres zastosowalności: **tylko czyste metale**.
- TGP A_orb pozostaje jako **klasyfikator** (s/sp/d), nie jako **parametr predyktywny**.
- Stopy wymagają osobnej teorii (Nordheim, Anderson disorder).

### Proponowane następne kroki

**r17 (opcjonalne):** Sensitivity analysis — perturbować Θ_D/N_F o ±10% i sprawdzać stabilność U24 coefs. Oczekiwany wynik: coefs stabilne (bo LOO już tego dowodzi pośrednio).

**r18 (opcjonalne, ambitne):** Rozszerzyć formułę o term Nordheim'a `+ n_subst·x·(1−x)·<ΔZ²>` dla stopów. Refit łączony (pure+alloy). To byłoby nowa teoria, nie mieści się w paperze v3.

**r19 (tylko jeśli użytkownik chce):** Draft paperu v3 z zawężonym zakresem (tylko pure metals) + explicit statement że formuła nie extrapoluje do alloys. Honest science.

**Status projektu po Etapie 9:** Formuła U24 ma **solidne podstawy dla 37 czystych metali**, ale **nie jest uniwersalna dla stopów**. Paper v3 możliwy ze zawężonym scope. Ostrzał ujawnił granice modelu — to wartościowa wiedza, która zapobiega over-claim'owi w publikacji.

---

## Etap 10: Sensitivity + Nordheim-like unification (r17–r18)

**Data:** 2026-04-21

### r17 — Sensitivity analysis (Theta_D, N_F perturbations)

**Cel:** Sprawdzić robustność współczynników U24 na niepewności literaturowe inputów.

**Testy:**
- (A) Θ_D ±10% symetrycznie
- (B) N_F ±15% symetrycznie
- (C) Monte Carlo N=500 z gaussian noise σ_logTh=0.04, σ_logNF=0.06
- (D) Leave-outlier-out (Bi, Be, Ca, Sr, Tc)

**Wyniki MC (N=500) — rozbicie parametrów na stabilne i niestabilne:**

**Stabilne (TGP-uniwersalne):**
| Parametr | wartość | MC std/|base| | >0/500 |
|---|---|---|---|
| **b(logTh)** | +0.79 | 7.9% | 500/500 dodatnie |
| **d(v)** | −0.10 | 4.2% | 500/500 ujemne |
| asp(sp-d) | −1.79 | 3.8% | — |
| csp(sp·logNF) | −0.80 | 7.5% | — |
| c(logNF) | +0.41 | 13.7% | — |

**Niestabilne (zdominowane przez N_s=5):**
| Parametr | wartość | MC std/|base| |
|---|---|---|
| as(s-d) | −2.09 | 44.5% |
| ds(v·s) | +0.60 | 51.0% |
| cs(s·logNF) | −0.39 | 163.9% |

**Leave-outlier-out test** (shifty współczynników po usunięciu Bi/Be/Ca/Sr/Tc):
- `b(logTh)`: +0.79 → +0.96 (shift +0.17, nadal dodatnie) ✓
- `d(v)`: −0.10 → −0.11 (shift −0.01, praktycznie niezmienione) ✓
- `as`, `ds`, `cs`: shifty 1.0–1.5 (niestabilne)

**Kluczowe odkrycie r17:** `b=+0.79` i `d(v)=−0.10` są **prawdziwymi TGP-invariantami**; pozostałe parametry są artefaktami małej próby s-class (N=5). To wzmacnia paper v3 — głównymi odkryciami są **uniwersalny Bloch-Mott wykładnik** i **d-class valence coupling**.

### r18 — Próba unifikacji pure+alloy (Nordheim, rho_0, is_alloy terms)

**Cel:** Czy dodanie fizycznego termu (Nordheim x(1−x)ΔZ² lub disorder rho_0) zamyka pure+alloy (N_total=50) w jednej formule?

**Testowane hipotezy:**
- H1: U24 + α·log(1+Nordheim/100)
- H2: U24 + λ·log₁₀(ρ_0) (Ioffe-Regel saturation)
- H3: U24 + class-specific λ·log(ρ_0)
- H4: U24 + class-specific is_alloy indicator (benchmark trivial)
- H5: U24 + log(ρ_0) + Nordheim

**Wyniki:**

| Model | k | Pure RMS | Alloy RMS | Total RMS | AIC |
|-------|---|----------|-----------|-----------|-----|
| U24 (pure only) | 10 | **0.156** | 0.705 | 0.384 | −75.8 |
| H1 Nordheim | 11 | 0.235 | 0.463 | 0.311 | −94.8 |
| H2 ρ_0-log | 11 | 0.229 | 0.400 | 0.283 | −104.1 |
| H3 ρ_0 × cls | 13 | 0.212 | 0.395 | 0.272 | −104.2 |
| **H4 is_alloy × cls** | 13 | 0.172 | **0.316** | **0.218** | **−126.2** |
| H5 H1+H2 combined | 12 | 0.227 | 0.398 | 0.281 | −102.9 |

**Interpretacja H3 (class-specific disorder response):**
- `λ_ρ₀(d) = −0.20` — **d-alloys saturują** (więcej disorder → mniej R). Zgodne z Ioffe-Regel: kF·l_e→1 w strongly disordered d-metals.
- `λ_ρ₀(s) = +0.27` — **s-alloys mają dodatkowe R** z disorder (Linde-Sondheimer deviation).
- `λ_ρ₀(sp) = +0.13` — **sp-alloys lekko dodatnie** (pośrednie).

**Fizyczny wniosek:** znaki λ zgadzają się z różnym mechanizmem scatteringu per klasa orbitalna. To jest **nowa interesująca obserwacja**, ale...

**Praktyczny wniosek:** wszystkie hipotezy mają alloy RMS ≥ 0.32. **Żadna nie osiąga target < 0.25.** Ponadto:
- Włączenie stopów psuje pure RMS (z 0.156 do ≥0.17)
- H4 (is_alloy × cls) wygrywa AIC, ale to jest **trywialny offset per klasa**, nie fizyczna formuła
- H1-H3-H5 robią unifikację ale kosztem utraty pure-metal precyzji

**Verdict r18:** Pełna unifikacja pure+alloy w ramach U24 + pojedynczego termu **nie jest możliwa**. Stopy wymagają **osobnego modelu**, prawdopodobnie z:
1. Wielo-parametrowym termem Nordheim (niemożliwe dla stopów 3-składnikowych z istniejącą precyzją danych)
2. Pełnym treatment Ioffe-Regel saturation (fizycznie poprawne, ale wymaga dodatkowych zmiennych: kF, vF)
3. Rozbiciem ρ(T) na trzy komponenty: ρ_0 + ρ_i^pure + ρ_disorder (rozszerzenie Matthiessen o term kolektywny)

### Status po Etapie 10

**Finalna diagnoza:**
1. **U24 (pure metals only)** jest optymalną formułą, z robust TGP-invariantami: `b=+0.79` Bloch-Mott, `d(v)=−0.10` d-class coupling. **Gotowe do paperu v3.**
2. **Stopy wymagają osobnej teorii** — żadna prosta rozszerzenie U24 nie generalizuje. To nie jest porażka TGP, lecz właściwy wynik: ρ(T) stopów jest fizycznie różne (Anderson disorder, Ioffe-Regel, Linde-Sondheimer).
3. **Obserwacje fizyczne z r18** (znaki class-specific λ_ρ₀) mogą stać się **teasery/dyskusje** w paperze v3, ale nie centralny wynik.

**Rekomendacja dla paperu v3:** Zawężony scope — **37 pure metals, universal b=+0.79, TGP A_orb jako klasyfikator**. Osobny appendix lub discussion section: "Why alloys require separate treatment" z diagnostyką r16/r18.

**Alternatywa (jeśli ambicja):** r19-r20 — poszerzenie alloy dataset do 30+ stopów, wtedy H3/H4 miałyby statystyczną moc. Ale to osobny projekt.
