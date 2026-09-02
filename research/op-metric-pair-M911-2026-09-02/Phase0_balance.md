---
title: "Phase0_balance — LOCK: właściwa para metryczna (w, V_M9.1'') — krajobraz sektora grawitacyjnego, relaksacja/geneza i kreacja przy granicy metryki (N3 poprzednika)"
date: 2026-09-02
type: phase0-lock
tgp_owner: research/op-metric-pair-M911-2026-09-02
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-02: „ok, rozpisz cykl dla nowego agenta" (N3 z NEEDS op-metric-closure-relaxation: re-lock z parą (w, V_M9.1'') zamiast hybrydy z kanonicznym U — diagnoza rozjazdu PRIMARY↔C-BAR)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-metric-closure-relaxation-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metric-closure-relaxation-2026-09-02/NEEDS.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase 0 — LOCK cyklu `op-metric-pair-M911`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

---

## 0. Pytanie i dziedziczona diagnoza

Poprzednik (op-metric-closure-relaxation, Q-PASS-NUCLEATION) wykazał:
hybryda „czynnik metryczny w(ψ)=ψ/(4−3ψ) × kanoniczne U materii" jest
niekompatybilna (U(g_ceil)−U(1)<0 ⟹ biegun PRZYCIĄGA). Korpus ma
własną, kompatybilną parę sektora grawitacyjnego (G.0 closure LOCK
2026-05-02, sek08a ~340–380): **w(ψ)=ψ/(4−3ψ)** (eq:vol-element-M911)
+ **V_M9.1''(ψ) = −γψ²(4−3ψ)²/12** (podwójne zero w ψ=4/3)
+ **K(ψ)=K_geo·ψ⁴**. Reguła korpusu (sek08a ~215): cykle
gravity-related MUSZĄ używać V_M9.1'' — ten cykl jest pierwszym, który
robi to w programie granicy metametrycznej.

**Pytania binarne:**
- **Q-A (krajobraz, analitycznie):** czy sektor grawitacyjny
  (w, V_M9.1'', K) jest SAMODOMKNIĘTY — tj. efektywna gęstość
  ρ_eff(ψ)=w·V ma dobrze określone minimum próżniowe, granica ψ=4/3
  jest „pod górkę", a E jest ograniczone z dołu na dziedzinie ψ∈(0,4/3)
  BEZ podłóg/barier ad-hoc?
- **Q-B (RACHUNEK CENTRALNY):** czy relaksacja w tej parze wytwarza
  STRUKTURĘ: nukleację obiektów (pre-rejestrowany POZYTYW autora,
  detektory dziedziczone) lub nietrywialny stan stacjonarny — czy
  wszystko relaksuje do próżni?
- **Q-C (warunkowe przy PASS-STATIC):** widmo stabilności stanu.

Uwaga interpretacyjna zalockowana: Q-B-FAIL (wszystko do próżni) też
jest ważnym wynikiem — mówi, że kreacja przy granicy metryki NIE jest
własnością czystego sektora grawitacyjnego i kieruje program do genezy
Γ+s_i (poziom 0) lub sprzężenia z materią.

## 1. Model ZAMKNIĘTY

- **Funkcjonał:** E[ψ] = ∫ w(ψ)·[½K_geo ψ⁴|∇ψ|² + V_M911(ψ)] dx,
  dziedzina ψ∈(0, 4/3); K_geo=1, γ=1 (bezwymiarowo). Implementator
  MUSI odtworzyć dokładne formy z sek08a (~340–380) i zamrozić
  z cytatem w method_decisions; ewentualny rozjazd form (np. czy w
  mnoży też człon kinetyczny) rozstrzygnąć CYTATEM z eq:S-unified-psi /
  eq:vol-element-M911 PRZED startem, oba odczyty odnotować.
- **Zmienna:** ψ (nie g); mapowanie ψ=g² przy imporcie startów.
- **Starty (zalockowane 3):** (i) geneza: próżnia ψ=1 + szum, 3D pudło
  L=4π, N∈{48,64}, seed=20260903, amp=1e−3; (ii) bump radialny
  (gauss, ψ_max=1.30 — poniżej granicy; R=60, h∈{0.025,0.0125});
  (iii) sieć 2π: ψ=g² z `../op-3d-canonical-lattice-2026-08-31/
  Phase2_backgrounds3d.npz` (READ-ONLY; UWAGA: ψ_max=g_max²≈2.17 > 4/3
  ⟹ przed użyciem przeskalować amplitudę tak, by ψ_max=1.30 —
  procedura FROZEN; oryginalny kształt raportować).
- **Detektory (dziedziczone reguły, progi w ψ zalockowane TERAZ):**
  dolny: spójne regiony ψ < 5/6 (średnia progu między próżnią 1
  a zerem V′ w 2/3); górny: ψ > 7/6 (średnia między 1 a 4/3);
  reguły zbieżności jak u poprzednika (≥10 j.cz., siatka×dt, ±1);
  kontrola z zasianym obiektem (osiągalny FAIL) nieusuwalna.
- **Obsługa granicy:** brak barier ad-hoc (to jest sedno Q-A);
  numeryczna regularizacja WYŁĄCZNIE w pasie ψ>4/3−1e−6, udokumentowana;
  wejście pola w pas = BREAKDOWN-BOUNDARY (klasyfikacja, nie werdykt).
- **Rejestr WEJŚĆ:** seed, ψ_max startów=1.30, progi detektorów,
  K_geo=γ=1, procedura skalowania sieci.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — Q-A: krajobraz (sympy + wykresy kontrolne, zero relaksacji)
- P1a: punkty krytyczne i krzywizna ρ_eff=w·V na (0,4/3): położenie
  minimum, wartość w próżni, zachowanie przy 0⁺ i 4/3⁻; analogicznie
  pełna gęstość z członem gradientowym (znak współczynnika kinetycznego
  na dziedzinie). **Q-A-PASS:** minimum w ψ*∈(0,4/3) z dodatnią
  krzywizną ORAZ ρ_eff(4/3)>ρ_eff(ψ*) ORAZ brak kierunku ucieczki
  (E ograniczone z dołu). **Q-A-FAIL:** przeciwnie (raport wprost —
  wtedy sektor też wymaga domknięć i program wraca do NEEDS).
- P1b (gate zgodności): wartości ρ_eff w {0.5, 1, 7/6, 1.3} policzone
  sympy vs float — zgodność 1e−12 (osiągalny FAIL implementacji).

### Phase 2 — bramka maszynerii relaksacyjnej
- P2a: próżnia ψ≡1 bez zaburzeń zostaje (dryf ≤1e−10; jeśli P1a wykaże,
  że ψ=1 NIE jest punktem krytycznym pełnego funkcjonału — użyć ψ*
  z P1a jako próżni odniesienia i odnotować).
- P2b: kontrola detektorów z zasianym obiektem dolnym i górnym
  (wykrycie 1±0) oraz na czystej próżni (zero alarmów).
- FAIL ⟹ STOP.

### Phase 3 — Q-B: RACHUNEK CENTRALNY (relaksacja)
- 3 starty × 2 siatki (+dt/2 przy zdarzeniach); gradient flow
  do ‖ψ̇‖∞≤1e−8 LUB nukleacji LUB t_max=200.
- **Werdykty (litera):** Q-B-PASS-NUCLEATION (detektor, zbieżnie) /
  Q-B-PASS-STATIC (stacjonarny NIESTAŁY ‖ψ−const‖∞≥0.05, zbieżny
  ≤5e−3; sanity: ψ_max vs 4/3) / Q-B-FAIL (wszystko do jednorodnej
  próżni) / Q-B-INCONCLUSIVE (w tym BREAKDOWN-BOUNDARY zbieżny —
  raportować jako osobną kategorię deskryptywną: „pole wybiera
  granicę", bez nadinterpretacji).
- Deskryptywnie obowiązkowo: los startu (iii) — czy struktura sieci
  przeżywa w sektorze metrycznym.

### Phase 4 — Q-C (tylko przy Q-B-PASS-STATIC)
Widmo wokół stanu (Hessian pełnego E z czynnikiem w — konsekwentnie);
mody zerowe przed interpretacją; zbieżność ≤0.05·max(|ω²_min|,0.1);
Q-C-PASS: ω²_min≥−1e−3. Przy PASS-NUCLEATION: charakterystyka kaskady
bez progów.

## 3. Forbidden moves
Dziedziczone z poprzednika (§4) + (a) formy (w, V, K) wyłącznie
z sek08a z cytatem — zakaz modyfikacji wykładników/współczynników;
(b) zakaz DODAWANIA podłóg/barier (sedno Q-A); (c) detektory
niezmienialne po pierwszym biegu; (d) rdzeń .tex/STATE/git nietykane
(commit robi sesja główna); (e) katalogi innych cykli tylko odczyt;
(f) INCONCLUSIVE/BREAKDOWN-BOUNDARY ≠ pozytyw.

## 4. Deliverables
`Phase_method_decisions.md` (FROZEN, z cytatami form), `Phase1_landscape.py`
+output, `Phase2_gate.py`+output, `Phase3_relax_M911.py`+output+npz,
`Phase4_spectrum.py` (warunkowy)+output, `Phase_FINAL_close.md`,
`NEEDS.md`, `README.md`.

## 5. Drzewo decyzyjne
```text
Q-A-FAIL → raport + Phase 2–3 NADAL wykonać (krajobraz zły ≠ brak struktur),
   ale z flagą; NEEDS: domknięcie sektora grawitacyjnego
P2 FAIL → STOP (maszyneria)
Q-B-PASS-NUCLEATION → NEEDS: kreacja przy granicy metryki potwierdzona
   we WŁAŚCIWYM sektorze (mocna wersja hipotezy autora); PILNE: geneza
   Γ+s_i + user-gate dopiski core (retrospektywa Q-FAIL-i)
Q-B-PASS-STATIC → Phase 4; przy Q-C-PASS: pierwszy STABILNY obiekt
   programu w sektorze metrycznym (user-gate: konsekwencje dla znaku W)
Q-B-FAIL → NEEDS: czysty sektor grawitacyjny nie kreuje — hipoteza
   kierowana do genezy poziomu 0 / sprzężenia z materią (osobny LOCK)
Q-B-INCONCLUSIVE/BOUNDARY → NEEDS: wniosek metodologiczny + kategoria
   „pole wybiera granicę" deskryptywnie
```

---

**LOCK ZAMKNIĘTY 2026-09-02. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move.**
