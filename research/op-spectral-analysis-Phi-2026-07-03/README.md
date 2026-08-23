---
title: "op-spectral-analysis-Phi — CP-7/L03: numeryczna diagonalizacja operatora fluktuacji na tłach vacuum/Yukawa/soliton + ghost wall"
type: research_cycle
status: CLOSED-EXECUTED
phase: FINAL
folder_status: closed-mixed-verdict
claim_status: "CLOSED-EXECUTED 2026-07-03 — sektor grawitacyjny F-A: spektralnie CZYSTY (C1/C2 PASS; C3 PASS na tłach EL-konsystentnych do amp 1,28; linear-profile artefakty zdiagnozowane). Sektor solitonowy F-S (forma korony): wynik NEGATYWNY — próżnia g=1 funkcjonału solitonowego jest tachioniczna (kontinuum od −1, potwierdzone box-count = floor(R/π)); mody zlokalizowane: e=0, μ=2, τ=3 (l=0). Ghost wall: AKTYWNY dla μ/τ (min f(g)≈0,04); forma substratowa α=1 traci solitona τ (kolaps g→0,16). Twierdzenie syntezy 2026-05-06 (σ⊂[0,∞) dla wszystkich tł) OBALONE dla sektora solitonowego — obowiązuje TYLKO w F-A."
created_date: 2026-07-03
closed_date: 2026-07-03
authorization: "User 2026-07-03: 'Ok działaj z op-spectral-analysis-Phi' (CP-7 z AUDYT_GLEBOKI_2026-06-28 §3 Tier 2)"
anti_lakatos_lock: PRESERVED
related:
  - "[[../../audyt/L03_K_phi_stability/README.md]]"
  - "[[../op-L03-spectral-stability-2026-05-06/spectral_synthesis.md]]"
  - "[[../../meta/AUDYT_GLEBOKI_2026-06-28.md]]"
---

# op-spectral-analysis-Phi (CP-7 / L03)

## Cel (AUDYT_GLEBOKI CP-7)

Diagonalizacja `K(φ)∂² + V''(φ)` na tłach próżnia/Yukawa/soliton
(sympy + BVP numeryczny), `m² > 0` dla `φ > 0`, rozstrzygnięcie
ghost-wall przy `φ → 0`. Poprzedni cykl (2026-05-06) był analityczny;
audyt: „tylko punkt próżni domknięty".

## Metoda

- **Phase 0** ([[Phase0_balance.md]]): LOCK operatora, teł, kryteriów
  C1–C6 i forbidden moves PRZED napisaniem kodu.
- **Phase 1** ([[Phase1_operator_sympy.py]] → [[Phase1_output.txt]]):
  dokładna druga wariacja `E = ∫ ½F(u)|∇u|² + W(u)` (sympy, 10/11 PASS):
  `L̂[v] = −(1/r²)(r²F v′)′ + Q v`,
  `Q = W″(u₀) − ½F″(u₀)u₀′² − F′(u₀)[u₀″ + (2/r)u₀′]`.
  Zweryfikowano tożsamości EL↔EOM: akcja F-A ⇔ `thm:field-eq` (α=2),
  funkcjonał F-S ⇔ ODE korony (a3d/ls10), F-S′ ⇔ ODE słownikowe α=1.
- **Phase 1b** ([[Phase1b_endpoint_classification.py]]): klasyfikacja
  końca ψ→0 (C6b).
- **Phase 2** ([[Phase2_bvp_spectrum.py]] → [[Phase2_output.txt]],
  [[Phase2b_supplement.py]] → [[Phase2b_output.txt]]): samosprzężona
  dyskretyzacja S-L (siatka staggered, Dirichlet w R), zbieżność
  N ∈ {2000,4000,8000}, R ∈ {40,60,80}; tła: próżnie, Yukawa liniowa
  i **nieliniowa (Newton, residuum <3e−12)**, solitony e/μ/τ w formie
  korony (α=2, log, odbicia od ściany jak a3d) i substratowej (α=1).

## Werdykty względem zalockowanych kryteriów

| Kryt. | Wynik | Szczegół |
|---|---|---|
| C1 | **PASS** | sympy exact: m_sp² = γ na próżni F-A (U_A″(1)/K(1) = γ) |
| C2 | **PASS** | próżnia F-A: N_neg=0, krawędź 1,0027 (błąd 0,27%), zbieżne |
| C3 | **FAIL (raw) / PASS (EL-konsystentne)** | profile liniowe A=0,01/0,1: N_neg=1–2, ale λ_min **dywerguje** z N (−1,6e4→−3,6e7→…): artefakt tła nie-EL (1/r core poza reżimem liniowym) — nie spełnia zalockowanego wymogu potwierdzenia zbieżnością. Tła nieliniowe EL-konsystentne (Newton, amp = 0,08/0,34/0,66/0,93/**1,28**): **N_neg=0 wszędzie** |
| C4 | **WYNIK NEGATYWNY (zaraportowany wprost)** | patrz niżej |
| C5 | **CONFIRMED** | próżnia F-S: krawędź −0,9973 ≈ −1 = W″(1)/f(1); dualizm L04 zmierzony u źródła |
| C6a | **ROZSTRZYGNIĘTE** | e: min g=0,932 (ściana nietknięta); μ: 1 odbicie, min f(g)=0,045; τ: 3 odbicia, min f=0,038 — **ściana AKTYWNA dla gen 2–3**; substrat α=1: τ **kolabuje** do g≈0,16 (brak ściany), profil kończy się przy r≈3 |
| C6b | **ROZSTRZYGNIĘTE** | koniec ψ→0 (F-A): **miękki** (U″(χ)→0⁻, wiodąco −4·3^{1/3}γK·χ^{1/3}); hipoteza bariery odpychającej OBALONA (Phase1 T7b FAIL, uczciwie); wykluczenie ψ→0 jest **aksjomat-warunkowe** (łańcuch no-absolute-vacuum W(1)=γ/3≠0 + ρ≥0 ⇒ δψ>0), nie dynamiczne |

## C4 — wynik centralny (negatywny dla claimu syntezy 2026-05-06)

1. **Próżnia funkcjonału solitonowego F-S jest tachioniczna:** kontinuum
   σ(L̂) zaczyna się od **−1** (jedn. γ=1), bo W″(1) = −1 < 0. Liczba
   modów ujemnych w pudle = floor(R/π) **dokładnie** (12/19/25 dla
   R=40/60/80) — to stany kontinuum, nie niestabilności solitonu.
2. **Mody zlokalizowane** (λ < −1−10⁻³, zbieżne na siatkach):

   | Tło | l=0 | l=1 |
   |---|---|---|
   | soliton e (0 odbić) | **0** | 0 |
   | soliton μ (1 odbicie) | **2** (−1,282; −1,057) | 1 (−1,031) |
   | soliton τ (3 odbicia) | **3** (−4,216; −1,114; −1,010) | 2 |
   | substrat: e / μ / τ | 0 / 2 / 1* | 0 / 0 / 0* |

   *τ w substracie: profil skolabowany (nie-EL do r_max) — niemiarodajne.
3. **Interpretacja (bez nadinterpretacji):** elektron jest konfiguracją
   spektralnie czystą WZGLĘDEM swojego (tachionicznego) kontinuum;
   μ i τ są **punktami siodłowymi** funkcjonału E_S (indeks 2 i 3).
   Twierdzenie `thm:spectral-synthesis-L03` („σ(L̂) ⊂ [0,∞) dla
   wszystkich fizycznych tł") jest **prawdziwe tylko w formulacji
   grawitacyjnej F-A**; w formulacji solitonowej (tej, w której
   istnieją tła korony) jest **obalone**. Przyczyna u źródła: synteza
   2026-05-06 założyła Q→+γ asymptotycznie — to własność F-A, nie F-S
   (tam Q→−1). Konflacja formulacji = dokładnie dualizm L04.
4. **Czego ten wynik NIE unieważnia:** dopasowań mas (r₂₁, r₃₁ przez
   Koide) — to własności profili ODE, nie spektrum fluktuacji. Unieważnia
   claim „pełna stabilność spektralna wspiera koronę" i przenosi ciężar
   stabilności μ/τ na dynamikę NIE-wariacyjną (odbicia od ściany) lub
   przyszłą stabilizację więzem (typ Q-ball, ładunek/masa ustalona) — OPEN.
5. **Obserwacja (SPECULATIVE, zero claimów):** liczba modów
   zlokalizowanych l=0 rośnie z generacją (0/2/3) i koreluje z liczbą
   odbić (0/1/3). Czy indeks siodłowy = wewnętrzna etykieta generacji —
   pytanie otwarte (wymaga interpretacji dynamicznej ściany).

## Ghost wall — rozstrzygnięcie (cel CP-7)

- **F-A (grawitacja):** ściana nie istnieje (K=ψ⁴>0); koniec ψ=0 jest
  miękki, wykluczony aksjomatycznie (C6b) — status `cor:ghost-artifact`
  dla sektora grawitacyjnego POTWIERDZONY, ale z doprecyzowaniem:
  wykluczenie NIE jest dynamiczne.
- **F-S (korona, log α=2):** dla μ/τ ściana jest **aktywnym składnikiem
  dynamiki** (odbicia; min f(g) ≈ 0,04 — operator na granicy utraty
  eliptyczności). Odbicie nie jest rozwiązaniem EL (regularyzacja
  ad-hoc a3d: g′→−g′ przy g*+0,005) ⇒ analiza spektralna μ/τ jest
  **zależna od regularyzacji ściany**. `cor:ghost-artifact` („ściana =
  artefakt log-aproksymacji") jest w NAPIĘCIU z faktycznym użyciem
  ściany przez mechanizm selekcji N=3 (#56): jeśli ściana to artefakt,
  mechanizm odbić wymaga niezależnego uzasadnienia; jeśli ściana jest
  fizyczna, `cor:ghost-artifact` wymaga zawężenia zakresu.
- **F-S′ (substrat α=1, preferowany przez sek08b):** brak ściany ⇒
  soliton τ (g₀=3,189) **kolabuje** (min g=0,158, profil urywa się
  przy r≈3) — substratowa forma NIE reprodukuje mechanizmu gen-3
  korony w obecnej postaci. (Uwaga: `cor:alpha1-preferred` deklaruje
  r₂₁=206,8 w substracie — zgodne: e/μ przeżywają; problem dotyczy
  wyłącznie τ.)

## Dyspozycja L03 (aktualizacja po CP-7)

- **Sektor grawitacyjny (F-A): L03 CLOSED-RESOLVED numerycznie**
  (diagonalizacja wykonana; próżnia + Yukawa EL-konsystentna do
  amp 1,28: σ ≥ 0; koniec ψ→0 sklasyfikowany).
- **Sektor solitonowy (F-S): L03 → OPEN-RECLASSIFIED** — nie „brak
  analizy", lecz **zmierzony wynik negatywny**: tachioniczne kontinuum
  + siodłowość μ/τ. Wymaga interpretacji dynamicznej (ściana/więz) —
  łączy się z L04 (dualizm) i T-OP3/T-OP4 (korona, Limitations).

## Pliki

- [[Phase0_balance.md]] — LOCK (kryteria, forbidden moves)
- [[Phase1_operator_sympy.py]] / [[Phase1_output.txt]] — operator exact (10/11)
- [[Phase1b_endpoint_classification.py]] / [[Phase1b_output.txt]] — C6b
- [[Phase2_bvp_spectrum.py]] / [[Phase2_output.txt]] — diagonalizacja główna
- [[Phase2b_supplement.py]] / [[Phase2b_output.txt]] — amp→1,28; mody zlokalizowane; box-count
- [[NEEDS.md]] — proponowane edycje core (user-gated)

## Anti-Lakatos

✓ Kryteria i tła zalockowane przed kodem; ZERO zmian po uruchomieniu.
✓ Wyniki negatywne (T7b FAIL; C3 raw FAIL; C4 siodłowość μ/τ; obalenie
twierdzenia syntezy dla F-S) zgłoszone wprost, z liczbami i zbieżnością.
✓ Artefakty odróżnione od fizyki testem zbieżności (λ_min dywergujące ⇒
artefakt; λ zbieżne ⇒ fizyka operatora).
✓ Rdzeń .tex NIETKNIĘTY (propozycje w NEEDS.md, user-gated).
✓ Zero nowych claimów pozytywnych; obserwacja generacyjna oznaczona
SPECULATIVE.
