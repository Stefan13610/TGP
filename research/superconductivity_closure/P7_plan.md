# P7 — Plan domknięcia nadprzewodnictwa (post-P6)

**Data:** 2026-04-19
**Status P6:** complete (A+B+C+D), r=0.875 na 29 materiałach (global), r=0.930 na 16 kalibracyjnych.

## Punkt startu P7

Po P6 pozostają trzy kategorie otwartych kwestii:

### A. Outliers z ps17 master validation

| Material | dlog | Kierunek | Hipoteza |
|----------|------|----------|----------|
| Y_amb | +1.15 | OVER 14× | λ_sf ambient Y zaniżone |
| Th_amb | +0.88 | OVER 7.5× | λ_sf dla actinides |
| La_amb | +0.38 | OVER 2.4× | j.w. |
| Ce_5GPa | +0.55 | OVER 3.5× | η(5GPa) < 0.8 może być |
| H3S | -0.43 | UNDER 2.7× | S 3d hybrydyzacja |

### B. Fenomenologiczne parametry do pierwszych zasad

- **λ_sf(material)** — obecnie z literatury (Stoner enhancement). P7 powinien derywowac z TGP: λ_sf ~ I·N(E_F), gdzie I i N wychodzą z substratu Φ.
- **P_scale(element)** — obecnie fitted dla Ce (5.8 GPa). Przewidzieć dla Yb, Eu, Sm z 4f binding energy.
- **ω_phonon(material)** — obecnie z eksperymentu. Predykcja z TGP substratu potrzebuje P6.B → P7 bridge.

### C. Room-temperature SC ambient (niezrealizowane)

Wymagania P6:
- Metaliczny
- ω_ph ≥ 200 meV
- a ≈ 4.088 Å (harmonika TGP)
- d-orbital (η=1)
- λ_sf ≈ 0

Żaden znany materiał nie łączy wszystkich. Kandydaci spekulatywni:
- Dopowany diament metaliczny (B- lub N-dopowany): ω~165 meV, a=3.57 Å
- BC₃ / BC₄ metaliczne warstwowe: syntetyczny wyzwał
- Li-H ambient-metaliczny: wymaga stabilizacji chemicznej

## Strategia P7

### P7.1 — Derivacja λ_sf z TGP substratu

Hipoteza: w TGP, paramagnon to lokalna fluktuacja $\Phi$ na energii $E_F$. Stoner I w TGP wychodzi z:
$$I = \left(\frac{\partial^2 E[\Phi]}{\partial m^2}\right)_{\!m=0}$$
gdzie $m$ to magnetyzacja substratu. Powinna być wyrazzalna przez $\Lambda_E$ i gęstość stanów.

Cel numeryczny (ps18 przyszły): reprodukować λ_sf dla V=0.6, Nb=0.2, FeSe_bulk=0.9 z pierwszych zasad (błąd ≤ 30%).

### P7.2 — Lantanowe pair-breaking (ZAMKNIĘTE ps21, 2026-04-19)

**Wnioski po analizie eksperymentalnej (PrH9 Tc=5K, NdH9 Tc=4.5K vs naiwne 200K):**

Główny efekt NIE jest η(P) ani P_scale — lecz **Abrikosov-Gorkov pair breaking**
od lokalnych 4f momentów magnetycznych. Formuła:
$$B_\text{PB}(\mu) = \exp(-\alpha_\text{PB} \cdot \mu_\text{eff}^2)$$

gdzie $\mu_\text{eff} = g_J\sqrt{J(J+1)}$ (Hund), a **α_PB = 0.2887 μ_B⁻²** (fit na PrH9+NdH9).

Zamiast P_scale(4f^n) z E_4f-binding, dostaliśmy **czystszy mechanizm μ_eff(Hund)**,
który przewiduje T_c dla wszystkich 15 lantanowców z jednej uniwersalnej stałej.

**Kluczowa predykcja P7.2:** LuH10 @ ~170 GPa → T_c ≈ 200-250 K 
(Lu³⁺ 4f¹⁴ closed, μ=0, brak pair-breaking).

Szczegóły: [[P7B_summary.md]], [[ps21_p7b_lanthanide_pair_breaking.py]].

### P7.3 — Room-temp candidates deep dive

**Kandydat 1: B-doped metalic diamond**
- Struktura: diament sp³ z domieszką B (akceptor)
- a(C-C) = 1.54 Å → kubiczna komórka a = 3.57 Å
- ω_phonon(C-C) = 165 meV
- Metaliczność: B > 3% dopant (realizowalne, obserwowane T_c=11K)
- TGP prediction: ω=165, orb=sp, a=3.57 → nieharmoniczne z a*, ale w ogonie
- Wymaga: fit do a*_TGP (diament ma a=3.57, a*=4.088 — off by 0.5Å)

**Kandydat 2: BC₃ / BC₄ layered**
- Syntetyczny boron carbide
- ω_phonon~180 meV (B-C bond stiffness)
- Warstwowa struktura — można dostroić a
- Wymaga synthesis (nie jest jeszcze eksperymentalnie stabilna)

**Kandydat 3: Yb superhydride pod P** *(ZREWIDOWANY ps18, 2026-04-19)*
- ~~YbH₉ @ 300 GPa: 215 K~~ → po refit z Yb4H23 (Sharps 2025, 11.5 K @ 180 GPa): **38 K**
- ~~YbH₁₀ @ 400 GPa: 267 K~~ → **71 K**
- P_scale_Yb = 552 GPa (nie 10 GPa jak wcześniej zakładano)
- **Wniosek:** Yb superhydrydy NIE są ścieżką do RT-SC. 4f¹⁴ zbyt zlokalizowane.

**Kandydat 4 (NOWY po weryfikacji 2026): pressure-quench Hg-cuprates**
- Hg1223 ambient po P-quench = 151 K ([Deng/Chu UH 2026](https://physics.aps.org/articles/v19/37))
- Hg1245 ambient TGP-pred = 178 K (nasz P6.A Tier 1)
- Potencjalnie Hg1245 + P-quench → bliżej 200 K bez wysokiego ciśnienia
- Wymaga: MBE + P-quench toolchain (np. 20 GPa → quench → ambient lock)

### P7.4 — Unified Tc master plot

Zintegrować wszystkie ps# w jeden skrypt produkujący:
1. Tabela 50+ materiałów z każdą klasą
2. Log-log scatter plot T_obs vs T_pred
3. Residuals histogram
4. Per-class RMS breakdown

### P7.5 — Relatywistyczna rewizja Eq. 5 (class-split η)  (OPEN, 2026-04-20)

**Źródła:** ps23 (korelacja), ps24 (OLS / 2-bin test), **ps25** (test propozycji $(1+\eta Z^2/137^2)$).

**Proponowana forma rozszerzenia Eq. 5 [eq:A-orb] w paperze SC:**
$$A_{\mathrm{orb}}(Z) = A_{\mathrm{orb}}^{(0)} \cdot \sqrt{1 + \eta_{\mathrm{orb}} \cdot Z^2/137^2}$$

gdzie $\eta_{\mathrm{orb}}$ jest **class-indexed** (nie uniwersalna — to kluczowy wniosek ps25).

**Wyniki ps25 (CORE N=25, F-test):**

| Wariant | # par. | RMS_log | Redukcja | F-test p |
|---------|-------|---------|----------|----------|
| Baseline (intercept) | 1 | 0.2073 | — | — |
| Globalne η (Z=Z_pair) | 1 | 0.2105 | −1.6% | 0.66 (n.s.) |
| Globalne η (Z=Z_heavy) | 1 | 0.2111 | −1.9% | 0.80 (n.s.) |
| **T6: η_f + η_6p (split)** | **2** | **0.1839** | **+11.3%** | **0.020*** |

**Wniosek:** uniwersalna stała η NIE działa (kompromis między dodatnimi i ujemnymi znakami na różnych klasach). Split na dwie klasy orbitalne rozbija P7.5 na **P7.5a** i **P7.5b**.

#### P7.5a — A_f reduction dla 4f (Kondo/Hund screening)

Wynik ps25-T4: **η_f = −1.308** (fit na N=4: LaH10, CeH9, CeH10, La_amb).

Konsekwencja dla Eq. 5:
$$A_f(Z) = 2.034 \cdot \sqrt{1 - 1.308 \cdot Z^2/137^2}$$

Przykładowe efektywne wartości:
- Z=57 (La): $A_f^{\mathrm{eff}} = 2.034 \cdot 0.774 = \mathbf{1.574}$
- Z=58 (Ce): $A_f^{\mathrm{eff}} = 2.034 \cdot 0.770 = 1.566$
- Z=90 (Th): $A_f^{\mathrm{eff}} = 2.034 \cdot 0.542 = 1.103$  ← wyjaśnia **Th_amb outlier**

**Efekty uboczne (korzystne):**
- Th_amb: dlog $+0.877 \to +0.156$ (outlier praktycznie zamknięty)
- Ce_5GPa: dlog $+0.549 \to +0.317$ (redukcja 42%)
- F-test standalone: F(1,24) = 5.46, **p = 0.028***

**Fizyka:** 4f elektrony podlegają Hund/Russell-Saunders splittingowi (~1 eV) i Kondo screeningowi (niska gęstość stanów na Fermim w lantanowcach). Paper v1 używa A_f = 2.034 skalibrowane na Y2BaCuO (faktycznie pairing na Cu), ale ekstrapoluje to na właściwe 4f-hydridki — co je **nad-szacuje**. Poprawka Z² jest naturalnym nośnikiem Hund/SOC wagi.

**Rewizja do paperu SC v2:** zastąpić w Eq. 5 stałą $A_f = 2.034$ postacią $A_f(Z) = 2.034\sqrt{1 - 1.308\,Z^2/137^2}$, przy czym fit $\eta_f$ jest outcomem P7.5a i powinien być przedstawiony jako nowa stała TGP (ósma obok $C_0$, $a^{\ast}$, $\sigma$, $\eta$, czterech $A$'s).

#### P7.5b — SOC correction dla 6p-heavy spectators (Pb, Hg, Tl, Bi)

Wynik ps25-T5: **η_6p = +0.331** (fit na N=7: Pb, Hg_elem, BSCCO2212, Tl2212, Tl2223, Hg1223, Bi2201), aplikowane via **Z_heavy** (ciężki kation, nie Cu).

Konsekwencja: dla cupratów z ciężkim kationem (Bi₂, Tl, Hg):
$$A_d^{\mathrm{eff}}(Z_{\mathrm{heavy}}) = 0.310 \cdot \sqrt{1 + 0.331 \cdot Z_{\mathrm{heavy}}^2/137^2}$$

Przykładowe wartości boostu na T_c $\sim A^2$:
- Z_heavy=80 (Hg1223): boost $(1 + 0.113)^2 = 1.24\times$
- Z_heavy=83 (BSCCO, Bi2201): boost $1.25\times$
- Z_heavy=82 (Pb, Tl2212): boost $1.25\times$

**Status:** P7.5b standalone jest borderline (p = 0.24), ale w T6 **łącznie z P7.5a** osiąga znaczność p = 0.020. Mniej "statystycznie dojrzałe" niż P7.5a — wymaga więcej danych (np. Hg1245 MBE, Tl2223 oczyszczone próbki).

**Fizyka:** 6p elementów ciężkich (Bi, Tl, Hg, Pb) mają SOC > 1 eV; pojawia się jako dodatkowy kanał mieszania stanów parujących Cu 3d z kationowym 6p. W paperze v1 brak — Tl2212, Hg1223 są niedoszacowane o ~factor 1.25×.

**Rewizja do paperu SC v2:** rozszerzenie A_orb o Z_heavy-zależny czynnik dla d-orbital cupratów, z fitem η_6p wykonanym na zunifikowanym zbiorze (potrzebne ≥ 10 materiałów).

**Status P7.5 ogólny:**
- **P7.5a (4f): zalecana inkluzja do SC v2** (p = 0.028, wyjaśnia Th_amb)
- **P7.5b (6p spectator): kandydat do inkluzji** (borderline, potrzebuje więcej danych)
- **η_d, η_s, η_sp ≈ 0:** 3d/4d/s/sp nie wymagają relatywistycznej poprawki (brak sygnału w ps25-T3)

**Cross-refs:** [[OSTRZAL_ROWNANIE_Tc_SERIA2_2026-04-20.md]], [[ps23_residual_correlation_Z_Mmol.py]], [[ps24_relativistic_Z_split.py]], [[ps25_eta_relativistic_correction.py]], [[ps25_results.txt]], [[P7.5_subproblems.md]].

## Predykcje experimentalne do falsyfikacji P6

1. **Hg1245 z SrTiO3 podłożem** (Tier 1, MBE): T_c_pred = 178 K
   - Falsyfikuje P6.A (n-dep) + P6.B (interface boost)
   - Jeśli obserwowane T_c > 200 K: P6 pełne validation
   - Jeśli T_c < 100 K: P6.A formula wymaga rewizji

2. **FeSe/BaTiO₃** (alternate substrate): T_c_pred = 120 K
   - Falsyfikuje hipotezy o specjalności SrTiO3 Fuchs-Kliewer

3. **Pd-H ambient** (ps15 scenariusz P6.D SF damping): T_c_pred = 86 K
   - Falsyfikuje P6.D (obserwowane PdH ma T_c~9K przy małej zawartości H)
   - Wymaga wyższej zawartości H i ω > 60 meV

4. **YbH₉/YbH₁₀ high-P**: T_c_pred **38-71 K** (po rewizji ps18)
   - Test orbital switching dla Yb z P_scale=552 GPa
   - Synteza nie jest jeszcze wykonana
   - Weryfikator: Yb4H23 @ 180 GPa = 11.5 K (Sharps 2025) już zwalidowany

5. **Hg1223 pressure-quench Tc=151K** (Deng/Chu 2026) — POTWIERDZONE
   - Nasz P6.A undershootuje (pred 93K), ale kierunek hipotezy OK
   - Hg1245 + quench to najbardziej realistyczny kandydat RT-SC dla P7

## Pliki P7 (status)

- **ps18**: ✅ weryfikacja + P_scale_Yb refit (552 GPa)
- **ps19**: ✅ λ_sf pierwsze zasady TGP (κ_TGP = 2.012)
- **ps20**: ✅ master plot, N=31 materiałów, r=0.875
- **ps21**: ✅ lantanowe pair-breaking (α_PB = 0.2887, predykcje 15 Ln)
- **ps22**: ✅ open candidates predictions (LuH10, ScH9, Li2AuH6)
- **ps23**: ✅ residual correlation test (P7.5 źródłowy: r=−0.41 dla Z_heavy)
- **ps24**: ✅ relativistic Z-split / Z_heavy correction test (P7.5)
- **ps25**: ✅ η(1+ηZ²/137²) proposal test → P7.5a (η_f=−1.308) + P7.5b (η_6p=+0.331)
- **ps26+**: (przyszłe) P7.3 Ce Kondo screening; konkretne synthesis proposals dla LuH10
- **dodatekP7**: (przyszłe) room-temp ambient summary + eksperymentalne predykcje P7.1+P7.2

## Co już zamknięte (nie dotykamy)

- P4: sektor cząstek → dodatekP4
- P5: IFF + T_c-skala → dodatekP5
- P6: cuprates + phonon + orbital + magnetyzm → dodatekP6
- Wszystkie ps1-ps17

## Logistic note

Gdyby empiryk zrobił Hg1245/SrTiO₃ i uzyskał T_c > 150 K, P6 ma strong validation.
Jeśli T_c < 30 K, trzeba przepisać P6.A (layer factor √n prawdopodobnie wrong).
