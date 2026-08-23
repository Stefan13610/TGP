---
title: "op-wall-dynamics — interpretacja ściany duchowej + stabilność solitonów μ/τ z więzem budżetu tworzonej przestrzeni"
type: research_cycle
status: CLOSED-EXECUTED
phase: FINAL
folder_status: closed-mixed-verdict
claim_status: "CLOSED-EXECUTED 2026-07-03 — W1 NEGATYWNE (wersja liniowa): mody μ/τ przetrwały wszystkie więzy K1–K3 i pre-deklarowane kombinacje; minimum osiągalne μ:1/τ:1 (pary z K4), rezydualny mod NIE jest kierunkiem rodziny (overlap ≈0,005 ≪ 0,9). Struktura: K4 (budżet rdzeniowy) usuwa dokładnie mody GŁĘBOKIE (−1,282/−4,216). W2 NEGATYWNE dla gładkich regularyzacji: soft-wall f_ε — τ KOLABUJE dla każdego ε (odbicie funkcjonalnie konieczne także wśród gładkich modeli EL); λ_min(ε→0) nie zbiega; dryf r₂₁ (μ-only) +1,9%…+23% ≫ 0,1%. W3a: g₀_wall=1,6114 ≈ 8/5 (0,71%) — górny próg g_crit (H7) pokrywa się z progiem aktywacji ściany dolnej g*: dwa progi = jeden mechanizm ścienny; B_core bez ekstremum przy progach (max ≈3,06); E_core nierozstrzygalne (szum kinków). Hipoteza budżetowa autora: obalona w wersji LINIOWEJ, wzmocniona strukturalnie (K4/g_crit⟺g*) — wymaga więzu nieliniowego/ładunkowego."
created_date: 2026-07-03
closed_date: 2026-07-03
authorization: "User 2026-07-03: 'rozpisz pierwszą fazę 0 cyklu op-wall-dynamics, zajmiemy się tym w następnej sesji (osobny agent)'; realizacja Phase 1–3: sesja #62 (osobny agent, zgodnie z autoryzacją)"
anti_lakatos_lock: PRESERVED
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]"
  - "[[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-03.md]]"
---

# op-wall-dynamics (realizacja Phase 1–3)

## Cel (LOCK)

Rozstrzygnąć, czy mody ujemne solitonów μ/τ z CP-7 (l=0: μ: 2, τ: 3
zlokalizowane) znikają po nałożeniu fizycznego więzu budżetu tworzonej
przestrzeni (hipoteza autora) lub poprawnego modelu ściany — tj. czy
stabilność korony jest twierdzeniem warunkowym, a nie założeniem.
Kryteria W1–W3, więzy K1–K4, kombinacje i rodzina regularyzacji f_ε
zalockowane w [[Phase0_balance.md]] PRZED jakimkolwiek obliczeniem.

## Metoda

- **Phase 1** ([[Phase1_constrained_spectrum.py]] → [[Phase1_output.txt]]):
  dokładne spektrum P L̂ P na podprzestrzeni więzów (rachuba inercji
  Haynswortha na tridiagonalu CP-7 + bisekcja; więz w współrzędnych
  symetryzowanych c_j = w(r_j)·r_j — waga B=r² uwzględniona).
  Walidacja: (a) unconstrained = CP-7 (μ: −1,2822/−1,0566;
  τ: −4,2164/−1,1142/−1,0098) ✓; (b) gęsta diagonalizacja P L̂ P przy
  N=2000: 4/4 PASS, max|Δλ| ≈ 5e−10.
- **Phase 2** ([[Phase2_wall_models.py]] → [[Phase2_output.txt]]):
  W2a Dirichlet na zbiorze kontaktu (dokładne odsprzęgnięcie segmentów)
  + suplement stożkowy skończonej amplitudy; W2b rodzina
  f_ε=½[f+√(f²+ε²)], re-solve profili EL i spektra, ε∈{0,2;0,1;0,05;0,02}.
- **Phase 3** ([[Phase3_budget_thresholds.py]] → [[Phase3_output.txt]]):
  sympy exact (f(g*)=0; tożsamość dH/dr=−(2/r)f g′²; dolne ograniczenie
  kontaktu W(g₀)≤W(g*)) + skan numeryczny g₀∈[1,04; 3,40] (120 punktów,
  B_core/E_core/odbicia/N_loc) + bisekcja progu kontaktu g₀_wall.

## Werdykty względem zalockowanych kryteriów

| Kryt. | Wynik | Szczegół |
|---|---|---|
| W1 | **WYNIK NEGATYWNY (wersja liniowa; zaraportowany wprost)** | Mody zlokalizowane przetrwały WSZYSTKIE pojedyncze więzy K1–K3 (μ: 2→2, τ: 3→2) oraz wszystkie pre-deklarowane rozszerzenia: K4 (μ: 2→1, τ: 3→2), pary K_i∧K_j (najlepiej K_i∧K4: μ→1, τ→1). Warunek „N(μ)=0" nigdy nieosiągnięty; rezydualny mod τ/μ NIE jest kierunkiem rodziny profili: overlap z ∂g/∂g₀ = 0,004–0,008 ≪ 0,9 (δ∈{1e−4,1e−3}). Zbieżność: N∈{2k,4k,8k} — identyczne liczby modów; R-kontrola: mody przy krawędzi (≈−1,006…−1,010) R-zależne (krawędź kontinuum), mody głębokie stabilne. Wg LOCK: „stabilizacja budżetem obalona w wersji liniowej; hipoteza autora wymaga więzu nieliniowego/innego ładunku". |
| W2a | **ROZSTRZYGNIĘTE: warunek jednostronny nie usuwa modów** | Zbiór kontaktu odbitego profilu jest niemal pusty (min g = 0,7876/0,7863 vs próg g*+0,01=0,7888; 0–2 punkty siatki) — Dirichlet na kontakcie zostawia μ: 2, τ: 3 (zmiany <0,3%). LCP: w ścisłej linearyzacji stożek przeszkody nieaktywny (przeszkoda w skończonej odległości g_eq−g* ≥ 0,005) — ograniczenie ODNOTOWANE per LOCK; suplement skończonej amplitudy (rzut na stożek, amp∈{0,01;0,05}): minimum Rayleigha = λ_min bez więzu (0 aktywnych węzłów). |
| W2b | **ROZSTRZYGNIĘTE: spektrum NIE jest zbieżne przy ε→0; mechanizm wrażliwy na model ściany (wynik negatywny do Limitations)** | (i) τ: profil EL z f_ε **KOLABUJE dla każdego ε** (urywa się przy r≈2,6–2,7, jak substrat α=1) — soliton τ istnieje TYLKO z ad-hoc odbiciem; wśród gładkich modeli f_ε mechanizm gen-3 nie ma reprezentanta EL. (ii) λ_min(ε) nie zbiega (μ: −1,39/−1,38/−1,29/−1,32; τ: −4,3/−226/−190/−4,8) i przy ustalonym ε spektra są NIEZBIEŻNE w N dla ε≤0,1 (μ) i dla wszystkich ε (τ); jedynie μ przy ε=0,2 zbiega (−1,389) — min f_ε ~ ε²/4 poniżej zdolności rozdzielczej siatek — kwantyfikacja regularization-dependence z CP-7. (iii) Tabela dryfu r₂₁/r₃₁(ε) formalnie nieobliczalna (brak ogona τ); dryf μ-only: r₂₁ +1,9% (ε=0,02) … +23% (ε=0,2) ≫ dopuszczalne 0,1% — **mechanizm korony jest wrażliwy na model ściany: raportowane wprost**. Baseline hard-wall odtworzony: r₂₁=206,73, r₃₁=3479,6 (PDG 206,77/3477,4). |
| W3a | **ROZSTRZYGNIĘTE (tak — dla pary progów; nie — dla B_core/E_core jako nośnika)** | **g₀_wall = 1,6114 ≈ 8/5 = g_crit (zgodność 0,71%)**: górny ogranicznik rdzenia (H7) pokrywa się numerycznie z progiem PIERWSZEGO kontaktu profilu ze ścianą dolną g* — dwa progi są dwiema projekcjami jednego mechanizmu ściennego (w polu: g*; w amplitudzie centralnej: g_crit). Sympy exact: f(g*)=0; dH/dr=−(2/r)f g′² (PASS); warunek konieczny kontaktu (bez tarcia) g₀ ≥ 1,1696 — tarcie 2/r przesuwa realny próg do 1,61. Natomiast wielkość budżetowa B_core NIE ma ekstremum przy żadnym progu (monotoniczna do maksimum przy g₀≈3,06±0,05, potem łagodny spadek); E_core: pochodna zdominowana szumem kinków odbić → **nierozstrzygalne z powodu szumu numerycznego** (raportowane per LOCK: „tak/nie/nierozstrzygalne z powodu X"). |
| W3b | **DESKRYPTYWNE (SPECULATIVE, zero claimów)** | N_loc(g₀) rośnie globalnie (0 przy g₀≤1,3 → 4 przy 3,4), ale NIE jest prostą funkcją liczby odbić: tuż po każdym skoku liczby odbić (g₀≈1,61/2,25/2,89) N_loc chwilowo SPADA (1,7→0; 2,3→2; 2,9→2); mody zlokalizowane istnieją już PRZED kontaktem ze ścianą (g₀≥1,4; λ_min pogłębia się −1,01→−3,53 przy zbliżaniu do progu). Punkty e/μ/τ (0/2/3) leżą na tej krzywej. Tylko dokumentacja korelacji. |

## Wynik centralny i struktura (bez nadinterpretacji)

1. **Hipoteza budżetowa autora w wersji LINIOWEJ jest obalona**: żaden
   z zalockowanych więzów liniowych (K1–K4, pary) nie zeruje indeksu μ
   (minimum osiągalne: 1) ani nie zostawia w τ modu rodziny profili.
   Zgodnie z zapisem LOCK-a hipoteza wymaga teraz więzu
   NIELINIOWEGO/ładunkowego (właściwa analogia Q-ball: ładunek z
   symetrii, nie liniowa deflacja) — do ewentualnego osobnego cyklu.
2. **Struktura godna odnotowania (fakt numeryczny, nie claim):** K4
   (lokalny budżet rdzenia, r<r_core) usuwa dokładnie mody GŁĘBOKIE
   (μ: −1,282; τ: −4,216), których nie ruszają więzy globalne K1–K3
   (te są niemal współliniowe: cos(c_i,c_j) > 0,995, zdominowane ogonem).
   Kierunek najgłębszej niestabilności μ/τ leży w podprzestrzeni zmian
   budżetu RDZENIA — jakościowo zgodne z intuicją hipotezy autora,
   ilościowo niewystarczające w wersji liniowej.
3. **Ściana jest funkcjonalnie konieczna i nie ma gładkiego zamiennika
   EL w rodzinie f_ε**: kolaps τ przy każdym ε (W2b) + niemal pusty
   zbiór kontaktu (W2a) + niezbieżność λ_min(ε→0). Napięcie z
   `cor:ghost-artifact` (sek08b) pogłębione: regularyzacja gładka NIE
   przywraca solitona gen-3.
4. **Dwa progi = jeden mechanizm**: g_crit = 8/5 (H7, górny) ≈ próg
   aktywacji ściany dolnej g* (0,71%) — pierwsze bezpośrednie numeryczne
   powiązanie obu progów. UWAGA: to wiąże progi ze SOBĄ (przez dynamikę
   ODE), ale nie wskazuje B_core/E_core jako wspólnej wielkości
   budżetowej (brak ekstremów przy progach).

## Czego wyniki NIE unieważniają

- Dopasowań mas korony r₂₁/r₃₁ w modelu hard-wall (baseline odtworzony:
  206,73/3479,6) — unieważniają natomiast ich NIEZALEŻNOŚĆ od modelu
  ściany (dryf ≫ 0,1% w rodzinie f_ε; τ bez odbicia nie istnieje).
- Werdyktów CP-7 (odtworzone co do 1e−4), mechanizmu Koide N=3 (#56),
  H7/H8 (wynik W3a je WZMACNIA strukturalnie).

## Pliki

- [[Phase0_balance.md]] — LOCK (kryteria, forbidden moves) — NIETKNIĘTY
- [[Phase1_constrained_spectrum.py]] / [[Phase1_output.txt]] — W1
- [[Phase2_wall_models.py]] / [[Phase2_output.txt]] — W2a/W2b
- [[Phase3_budget_thresholds.py]] / [[Phase3_output.txt]] — W3a/W3b
- [[NEEDS.md]] — propozycje edycji core (user-gated)

## Anti-Lakatos

✓ Kryteria, więzy, kombinacje i rodzina f_ε zalockowane PRZED kodem;
ZERO zmian po uruchomieniu (kombinacje K_i∧K_j i K4 były pre-deklarowane
w LOCK-u, nie dobrane post-hoc).
✓ Trzy wyniki negatywne zgłoszone wprost (W1 liniowy; W2b niezbieżność
+ wrażliwość r₂₁ na model ściany; W3a brak ekstremów B_core/E_core przy
progach), z liczbami i zbieżnością siatek.
✓ Mody raportowane wyłącznie z testem zbieżności (3 siatki + R-kontrola);
niezbieżności W2b raportowane JAKO niezbieżności, nie jako wartości.
✓ Metoda zwalidowana niezależnie (gęsta projekcja vs inercja: 4/4 PASS).
✓ Rdzeń .tex NIETKNIĘTY — propozycje w [[NEEDS.md]] (user-gated).
✓ Zero nowych claimów pozytywnych; W3b i obserwacja K4 oznaczone jako
deskryptywne/strukturalne; zgodność g₀_wall≈8/5 podana z niepewnością
(0,71%, konwencja guard g*+0,005), bez promocji do „exact".
