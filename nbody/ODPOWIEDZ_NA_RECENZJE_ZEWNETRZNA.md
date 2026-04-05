# Odpowiedz na krytykę zewnętrznego agenta (2026-04-05)

**Kontekst:** Agent zewnętrzny (prawdopodobnie GPT-class) przeanalizował PDF manuskryptu TGP i zaproponował roadmapę H1-H7 (wysoki), M1-M3 (średni), L1 (niski) do "domknięcia pod recenzję" w czasopiśmie grawitacji/kosmologii.

---

## Diagnoza ogólna

Krytyka jest **metodologicznie poprawna** w identyfikacji standardów recenzji (Will 2014, Planck 2018, DESI 2024, BIPM/CODATA). Jednak:

1. **Agent nie miał dostępu do pełnego drzewa źródeł** — widział tylko PDF. Wiele zidentyfikowanych "luk" jest adresowanych w plikach .tex i .py, których nie mógł przeczytać.
2. **Krytyka dotyczy głównie formatu/widoczności**, nie braku substancji.
3. **Kilka punktów jest naprawdę trafnych** — te wymagają pracy.

---

## Mapowanie punkt po punkcie

### H1: Operacjonalizacja "czasu substratu" — CZESCIOWO ADRESOWANE

**Co jest:**
- `sek08_formalizm.tex` (rem:substrate-time, linie ~3012-3024): dt_sub zdefiniowane
- `sek08c_metryka_z_substratu.tex`: wyprowadzenie g_tt = -c₀²·Φ₀/Φ z budżetu informacyjnego
- `sek04_stale.tex` (rem:c-interpretacja, rem:SI-constants): relacja c_lok → dt_sub

**Co brakuje:**
- ❌ Jawny słownik {dt_sub, dτ, zegar Cs} w jednym miejscu
- ❌ Analiza EEP/SEP: czy TGP łamie EEP (tak — c(Φ) jest nieuniwersalne) i jakie testy to wykrywają
- ❌ Mapowanie na bezwymiarowe obserwable zegarowe (porównania częstotliwości)

**Wysiłek:** 3-5 dni — głównie redakcja i formalizacja istniejącego materiału
**Blokuje:** H2, H5

### H2: Most Φ → g_eff^μν — W DUŻEJ MIERZE ZAMKNIĘTY

**Co jest:**
- `sek08c_metryka_z_substratu.tex`: pełne wyprowadzenie (propozycja, nie hipoteza)
- `sek08_formalizm.tex`: metryka wykładnicza (Stwierdzenie 208), konformalna jednoznaczność (prop:conformal-unique)
- `sek08_formalizm.tex` (linie ~3887-6656): dysformalność — kiedy potrzebna, kiedy nie
- Wynik: γ_PPN = β_PPN = 1 dokładnie w formie wykładniczej

**Co brakuje:**
- ❌ Jednoznaczne rozdzielenie "model minimalny" vs "rozszerzenia" w jednej tabeli
- ❌ Jawna deklaracja: w jakiej klasie transformacji jest metryka jednoznaczna

**Wysiłek:** 2-3 dni — porządkowanie istniejącego materiału
**Status:** Substancja jest, brak nawigacji/podsumowania

### H3: PPN/1PN w standardzie Will 2014 — CZESCIOWO

**Co jest:**
- γ_PPN = β_PPN = 1 (prop:PPN w sek08c, ex167 weryfikacja)
- Metryka wykładnicza automatycznie daje poprawne wyższe rzędy

**Co brakuje:**
- ❌ Komplet parametrów PPN (α₁, α₂, ξ, ζ₁-ζ₄) — nie tylko γ, β
- ❌ Jawne mapowanie na tabelę Will 2014
- ❌ Screening mechanism (czy potrzebny? Prawdopodobnie nie: m_sp ~ H₀ ~ 10⁻²⁶ m⁻¹ → ekranowanie na skalach kosmologicznych, ale nie w Układzie Słonecznym → potencjalny problem!)
- ❌ Preferred-frame effects (α₁, α₂)

**Wysiłek:** 2-4 tygodnie — wymaga obliczeń algebraicznych
**Ważność:** KRYTYCZNA — recenzent z grawitacji doświadczalnej to sprawdzi
**Ryzyko:** Jeśli substrat łamie izotropię Lorentza → α₁,α₂ ≠ 0 → potencjalny problem

### H4: Fale grawitacyjne — ZAMKNIĘTE W DUŻEJ MIERZE

**Co jest:**
- `sek07_predykcje.tex`: c_GW = c₀ dokładnie (stw. prop:cT, klasa: tożsamość algebraiczna)
- `scripts/gw/disformal_waveform.py`: pełne obliczenia tensorowe + SVT
  - GW170817 constraint: |c_GW/c_EM - 1| < 3×10⁻¹⁵ spełniony
- `scripts/gw/gw_breathing_mode.py`: mod oddechowy (spin-0, ekranowany masą m_sp)
- `dodatekC_ringdown.tex`: ringdown spektroskopia

**Co brakuje:**
- ❌ Jawna tabela: try/polaryzacja/prędkość → constraint → TGP prediction → status
- ❌ Dyspersja: ω² = c₀²(k² + m_sp²) → dla LIGO band to efektywnie c_GW = c₀

**Wysiłek:** 3-5 dni — głównie formatowanie istniejących wyników
**Kill-shot status:** PRZEŻYTY (c_GW = c₀ jest tożsamością w TGP)

### H5: Metrologia dynamicznych stałych — CZESCIOWO

**Co jest:**
- `sek04_stale.tex`: pełna sekcja z c(Φ), ħ(Φ), G(Φ)
- rem:SI-constants (linie 619-664): dyskusja SI i BIPM 2019
- rem:c-interpretacja: c₀ jest zawsze lokalne, zmienność widoczna w porównaniach

**Co brakuje:**
- ❌ Jawne przejście na bezwymiarowe kombinacje (α_EM, μ = m_p/m_e, etc.)
- ❌ Odniesienie do CODATA 2018/2022 procedur
- ❌ Protokoły: jak porównać zegary w różnych Φ (ACES, GPS, LLR)

**Wysiłek:** 5-10 dni
**Uwaga:** To jest bardziej "prezentacja" niż "brak fizyki"

### H6: Kosmologia Planck/DESI — CZESCIOWO

**Co jest:**
- `sek07_predykcje.tex` (linie 1068-1142): DESI CPL tension, w₀ ≈ -0.82
- `scripts/ex187_cosmo_data_confrontation.py`: H(z), G_eff(z), w_DE(z), porównanie z Planck/BBN/DESI
- `scripts/cosmology/tgp_cosmo.py`: solver kosmologiczny
- `scripts/cosmology/growth_factor_tgp.py`: wzrost struktur

**Co brakuje:**
- ❌ Formalne likelihood (χ², MCMC) z publicznymi danymi Planck PLA i DESI DR1 VAC
- ❌ Łańcuchy MCMC / contour plots (standard w kosmologii)
- ❌ Tabela bestfit z niepewnościami
- ❌ Porównanie modeli: ΔAIC, BIC vs ΛCDM

**Wysiłek:** 4-8 tygodni — to jest duży projekt
**Uwaga:** To jest **realnie najważniejsza luka**. Bez formalnych fits kosmologicznych recenzent odrzuci paper jako "qualitative"

### H7: Paczka reprodukowalności — CZESCIOWO

**Co jest:**
- `nbody/examples/verify_all.py`: 59/59 PASS
- `nbody/examples/verify_nbody_lyapunov_quick.py`: regresja Lyapunov
- `scripts/`: ex104-ex189 + ex195-ex197 + gw/ + cosmology/
- PLAN_DOMKNIECIA_v1.md: identyfikuje lukę

**Co brakuje:**
- ❌ README z instrukcją "one-command reproduce"
- ❌ requirements.txt / environment.yml
- ❌ Manifest: {skrypt → tabela/wykres w publikacji}

**Wysiłek:** 3-5 dni

---

## Priorytety realne (nasza ocena)

| Rang | ID | Zadanie | Wysiłek | Wpływ na recenzję |
|------|----|---------|---------|--------------------|
| 1 | H6 | Formalne likelihood Planck+DESI | 4-8 tyg | **KRYTYCZNY** — bez tego brak ilościowej kosmologii |
| 2 | H3 | Pełne PPN (komplet parametrów) | 2-4 tyg | **WYSOKI** — standard w grawitacji doświadczalnej |
| 3 | H1 | Słownik czasu substratu + EEP | 3-5 dni | ŚREDNI — konceptualnie ważne |
| 4 | H7 | Paczka reprodukowalności | 3-5 dni | ŚREDNI — standard współczesny |
| 5 | H5 | Metrologia + CODATA | 5-10 dni | NISKI-ŚR — ważne ale nie blokujące |
| 6 | H2 | Podsumowanie mostu metrycznego | 2-3 dni | NISKI — substancja jest |
| 7 | H4 | Tabela GW | 3-5 dni | NISKI — substancja jest |

### Co agent zewnętrzny przeoczył

1. **Sektor cząsteczkowy** — najsilniejszy argument TGP (r₂₁ = 206.77, α_s = 0.1174 z jednego parametru). Recenzent z grawitacji tego nie zna, ale to czyni TGP unikalnym. Agent nie uwzględnił tego w roadmapie.

2. **N-body jako wzorzec jakości** (M1 agenta) — to jest trafne. N-body ma najlepszą strukturę "definicje → równania → kod → testy" w całym TGP. Powinno być modelem dla PPN i kosmologii.

3. **Screening** — agent nie poruszył ryzyka: jeśli m_sp ~ √γ jest kosmologicznie mały, to pole skalarne TGP jest dalekozasięgowe w Układzie Słonecznym → potencjalne napięcie z testem Cassini (γ-1 < 2.3×10⁻⁵). TGP twierdzi c_GW = c₀ (tożsamość) i γ_PPN = 1 (twierdzenie), ale komplet parametrów PPN wymaga jawnego rachunku.

---

## Rekomendacja: co robić dalej

### Ścieżka A: "Paper N-body" (szybka)
Opublikować osobno notę N-body (tgp_nbody_results_clean.tex) w journalu typu CQG/PRD/EPJC. Wymaga:
- Tłumaczenie na angielski (TeX jest już po angielsku)
- Dodanie bibliography (.bib)
- Porównanie z literaturą (modyfikowane grawitacje 3-ciałowe)
- Reprodukowalność (H7 — szybkie)

**Wysiłek:** 2-3 tygodnie. **Wpływ:** pier wsza publikacja z TGP.

### Ścieżka B: "Full TGP paper" (długa)
Wymaga H1-H7 + M1-M2. Minimum 3-4 miesiące.
Sens: dopiero po N-body paper jako "proof of concept".

### Rekomendacja: **Ścieżka A → B**
1. Domknij N-body paper (2-3 tyg)
2. Równolegle: H3 (PPN komplet) + H7 (reproduce)
3. Potem: H6 (kosmologia formalna)
4. Na końcu: full paper

---

*Nota wygenerowana 2026-04-05 na podstawie analizy krytyki zewnętrznej + pełnego przeglądu drzewa TGP_v1.*
