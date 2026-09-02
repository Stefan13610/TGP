---
title: "Phase0_balance — LOCK: relaksacja z OBUSTRONNYM domknięciem — podłoga QB-2 + GRANICA METRYCZNA ψ=4/3 z korpusu (M9.1''); większe pudło dla genezy; nukleacja = pre-rejestrowany pozytyw (dziedziczone)"
date: 2026-09-02
type: phase0-lock
tgp_owner: research/op-metric-closure-relaxation-2026-09-02
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-02: „ok, działaj" na sekwencję N1→re-lock Q2 z NEEDS op-metametric-boundary (N1: górne domknięcie z korpusu — ZNALEZIONE dokumentacyjnie: biegun metryki efektywnej ψ=4/3, sek08a M9.1''; N2: większe pudło). Kryterium nukleacyjne autora dziedziczone."
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../op-metametric-boundary-2026-09-01/NEEDS.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase 0 — LOCK cyklu `op-metric-closure-relaxation`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

---

## 0. Rozstrzygnięcie N1 (dokumentacyjne, wejście do locka)

Górne domknięcie pola ISTNIEJE w korpusie i jest METRYCZNE (nie ad-hoc):
M9.1'' canonical (G.0 closure LOCK 2026-05-02; sek08a):
- element objętości: √(−g_eff) = c₀ψ/(4−3ψ) (eq:vol-element-M911; linie ~315, 358, 366),
- metryka: ds² = −c₀²(4−3ψ)/ψ dt² + ψ/(4−3ψ) δᵢⱼ dxⁱdxʲ (~linia 355),
- V_M9.1''(ψ) = −γψ²(4−3ψ)²/12 (podwójne zero w ψ=4/3).

**⟹ ψ_max = 4/3; w zmiennej g (ψ=g², bo φ=(Φ/Φ₀)^{1/2}): g_ceil = √(4/3)
= 1.1547005.** Przy ψ→4/3: g_tt→0 (nieskończona dylatacja czasu — GRANICA
METRYKI dosłownie), √−g→∞ (nieskończony koszt objętościowy — dynamiczna
bariera). To jest dokładnie „granica metryczna" z hipotezy autora
(op-blocked-soliton-bang: „wielki wybuch trwa na granicy metryki").

**Obserwacja pre-rejestrowana (do weryfikacji sanity, nie kryterium):**
tła łańcucha 1D z cyklu bloch miały g_max = 1.1406/1.1427/1.1429 —
tuż POD g_ceil=1.15470; zrelaksowane struktury zdają się naturalnie
żyć przy granicy metrycznej.

## 1. Pytanie binarne (Q — reponowane Q2 poprzednika)

**Czy z OBUSTRONNYM domknięciem (podłoga QB-2 + granica metryczna 4/3)
relaksacja sektora tachionowego osiąga stan metametryczny: (a) stacjonarny
NIESTAŁY stan strukturalny, LUB (b) kaskadę spontanicznej nukleacji
obiektów (pre-rejestrowany POZYTYW autora, detektor dziedziczony)?**

## 2. Model ZAMKNIĘTY

- **Sektor:** zunifikowany z czynnikiem metrycznym (wymóg korpusu:
  pytania metryczne → M9.1''): funkcjonał relaksowany
  E[g] = ∫ w(ψ)·[½K(g)|∇g|² + U(g)] dx, ψ=g², **w(ψ)=ψ/(4−3ψ)**
  (PRIMARY — domknięcie górne pochodzi z FIZYKI: w→∞ przy g→g_ceil);
  U=g⁷/7−g⁸/8, K=g⁴.
- **Wariant kontrolny C-BAR (nieusuwalny):** w≡1 + bariera C²
  (κ/3)(g−g_ceil)³ dla g>g_ceil — oddziela „fizykę metryki" od
  „jakiegokolwiek sufitu". Rozjazd PRIMARY↔C-BAR raportować.
- **Podłoga (dziedziczona):** kara C² z g_floor=0.5459 (próg QB-2
  Φ_c/Φ_vac=0.298; wrażliwość: {0.4438, 0.5753} tylko dla biegu
  genezowego).
- **Starty (zalockowane):** (i) geneza: próżnia+szum, **3D pudło L=4π**
  (N∈{48,64}), seed=20260902, amp=1e−3 — adresuje N2 poprzednika
  (L=2π tłumiło k≥1; teraz k_min=0.5<1 — pasmo tachionowe próżni
  częściowo w pudle, struktura MOŻE powstać); (ii) soliton radialny μ
  (R=60, h∈{0.025,0.0125}); (iii) sieć 2π z npz (READ-ONLY, N∈{32,48}).
- **Detektor nukleacji:** DZIEDZICZONY dosłownie z locka poprzednika
  (spójne regiony g<(1+g_floor)/2 ponad zasiane, ≥10 j.cz., zbieżność
  siatka×dt, ±1 obiekt). DODATKOWO (zalockowane TERAZ) detektor
  obiektów GÓRNYCH: spójne regiony g>(1+g_ceil)/2=1.0774 — kreacja
  przy granicy metrycznej liczy się tak samo (kierunek ucieczki
  poprzednika był W GÓRĘ).
- **Dynamika:** tłumiony gradient flow E (schemat FROZEN);
  do ‖ġ‖∞≤1e−8 LUB nukleacji LUB t_max=200.
- **Rejestr WEJŚĆ:** g_ceil=√(4/3) [z M9.1'' — CYTAT], g_floor (QB-2),
  seed, g₀_μ=φ·0.90548, β=γ=1.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — bramka maszynerii
- P1a: próżnia g≡1 bez zaburzeń zostaje w g≡1 (oba warianty domknięcia);
  gate dryfu ≤1e−10.
- P1b: pojedynczy bieg diagnostyczny startu solitonowego BEZ domknięć
  (reprodukcja BREAKDOWN g→+∞ poprzednika, t≈2.7–3.2) — osiągalny FAIL
  ciągłości z poprzednikiem.
- P1c (kontrola sektora stabilnego m²=+γ, nieusuwalna): relaksacja do
  próżni, ZERO nukleacji obu detektorów.
- Dowolny FAIL ⟹ STOP.

### Phase 2 — RACHUNEK CENTRALNY: relaksacja z obustronnym domknięciem
- 3 starty × wariant PRIMARY (+ C-BAR dla startu solitonowego i sieci;
  geneza: PRIMARY + wrażliwość podłogi ×3) × 2 siatki.
- **Werdykty (litera):**
  - **Q-PASS-NUCLEATION:** nukleacja (dolna lub górna) wg detektora,
    ZBIEŻNA — reżim kreacji potwierdzony.
  - **Q-PASS-STATIC:** stan stacjonarny NIESTAŁY (‖g−const‖∞≥0.05,
    zbieżny ≤5e−3). Sanity pre-rejestrowane: czy g_max stanu → okolice
    g_ceil (obserwacja z bloch).
  - **Q-FAIL:** wszystko relaksuje do jednorodnego.
  - **Q-INCONCLUSIVE:** załamanie nie-nukleacyjne / niezbieżność —
    wprost, nie pozytyw.

### Phase 3 (warunkowe przy Q-PASS-STATIC) — widmo na stanie metametrycznym
Jak u poprzednika (operator radialny / Bloch wg typu tła; mody zerowe
przed interpretacją; zbieżność ≤0.05·max(|ω²_min|,0.1); Q3-PASS:
ω²_min≥−1e−3). Przy Q-PASS-NUCLEATION: deskryptywna charakterystyka
kaskady bez progów.

## 4. Forbidden moves
Dziedziczone dosłownie z poprzednika (§4) + (8) czynnik metryczny
w(ψ) WYŁĄCZNIE w formie z eq:vol-element-M911 (zakaz modyfikacji
wykładników); (9) detektory (oba) niezmienialne po pierwszym biegu;
(10) wariant C-BAR nieusuwalny; (11) rdzeń .tex/STATE/git — jak zawsze
(git = sesja główna).

## 5. Deliverables
`Phase_method_decisions.md` (FROZEN; w tym schemat i implementacja w(ψ)),
`Phase1_gate.py`+output, `Phase2_two_sided_relax.py`+output+npz,
`Phase3_spectrum.py` (warunkowy)+output, `Phase_FINAL_close.md`,
`NEEDS.md`, `README.md`.

## 6. Drzewo decyzyjne
```text
P1 FAIL → STOP
Q-PASS-NUCLEATION → NEEDS: hipoteza „kreacja na granicy metryki" z nośnikiem
   numerycznym; PILNE: geneza Γ+s_i (poziom 0) + reinterpretacja Q-FAIL-i
   kanonicznych jako relaksacji (user-gate, dopiski core)
Q-PASS-STATIC → Phase 3; przy Q3-PASS: pierwszy STABILNY obiekt sektora
   tachionowego względem właściwej granicy (konsekwencje: znak W, sek08b —
   user-gate); sanity g_max vs g_ceil raportować
Q-FAIL → NEEDS: nawet obustronne domknięcie nie generuje struktury —
   hipoteza granicy osłabiona w klasie kontinuum; zostaje geneza poziomu 0
Q-INCONCLUSIVE → NEEDS: wniosek metodologiczny
PRIMARY↔C-BAR rozjazd → raportować deskryptywnie (fizyka metryki vs sufit)
```

---

**LOCK ZAMKNIĘTY 2026-09-02. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move.**
