---
title: "Phase0_balance — LOCK: następca Q1 metodą zwalidowaną w Q2 — samouzgodniony łańcuch periodyczny + analiza pasmowa Blocha w sektorze dynamicznym (ω²(n) modu runaway)"
date: 2026-08-31
type: phase0-lock
tgp_owner: research/op-bloch-chain-stability-2026-08-31
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-31: 'ok działaj z krokami 1,2,3' (krok 3 = cykl-następca Q1 wg NEEDS N2 cyklu op-bath-two-sectors: tło SAMOUZGODNIONE + periodyczne BC + analiza pasmowa Blocha — przeniesienie konstrukcji, która w Q2 dała czystą zbieżność, na sektor dynamiczny)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/NEEDS.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Phase 0 — LOCK cyklu `op-bloch-chain-stability`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria,
modele, progi, punkty skanu i forbidden moves zamknięte PRZED kodem.**

---

## 0. Pytanie binarne (Q)

**Czy w sektorze tachionowym (maszyneria 2 / EL Noty kanonicznej) istnieje
skończona gęstość struktur (separacja d łańcucha periodycznego), przy
której najniższa gałąź pasmowa fluktuacji wokół SAMOUZGODNIONEGO tła jest
dodatnia: ω²_min(d) = min_k λ_min(L̂_d(k)) > 0?**

- TAK ⟹ stabilność jest własnością konfiguracji o skończonej gęstości
  (ślepa plamka z retrospektywy 2026-08-16 zamknięta pozytywnie; ontologia
  „kolektyw żyje" dostaje pierwszy policzalny nośnik w sektorze dynamicznym).
- NIE ⟹ kąpiel periodyczna nie stabilizuje w klasie zbadanej — ostatnia
  wersja hipotezy stabilizacji gęstością w tym sektorze zamknięta negatywnie.
- Niezbieżność ⟹ INCONCLUSIVE, raportować jako niezbieżność.

**Zakres jawnie ograniczony (zapisany PRZED startem):** łańcuch 1D
periodyczny NIE jest solitonem radialnym 3D z #63. Cykl odpowiada na
pytanie STRUKTURALNE (czy periodyczna gęstość może odwrócić znak ω²
w sektorze tachionowym przy tle samouzgodnionym), nie na pytanie o μ
dosłownie. Przeniesienie na 3D = ewentualny osobny cykl (NEEDS).

## 1. Kontekst dziedziczony (NIE powtarzać tych rachunków)

- **#63 V3:** runaway solitonu μ w izolacji potwierdzony nieliniowo
  (M0-f_ε; λ_min(w1)=−1.389, t*=3.62).
- **op-bath-two-sectors Q1-INCONCLUSIVE:** komórka WS z obciętym profilem
  + zero-flux BC NIEZDOLNA rozstrzygnąć (0/15 zbieżnych; znak ω² podąża
  za rozmiarem komórki; komparator bez kąpieli identyczny).
- **op-bath-two-sectors Q2 (metoda zwalidowana):** samouzgodnione tło
  + komórka periodyczna dały pełną zbieżność (|dv|≤5.7e−6) i czystą
  kontrolę d=∞ — tę konstrukcję przenosimy na sektor dynamiczny.
- **ANALIZA_N2:** sektor tachionowy = EL Noty kanonicznej 2026-04-07
  (kinetyka K=g⁴, źródło g²(1−g), W″(1)=−γ); w 1D próżnia g=1 ma
  ogony oscylacyjne NIEZANIKAJĄCE (brak czynnika 1/r) — łańcuch
  periodyczny jest w 1D obiektem naturalnym, izolowany soliton nie.

## 2. Modele ZAMKNIĘTE

- **M-C (łańcuch, PRIMARY):** statyka 1D sektora tachionowego:
  `g″ + (2/g)g′² = g²(1−g)` (redukcja 1D EL Noty kanonicznej; forma
  z podzieleniem przez K jak w bramce A cyklu lattice-bath — implementator
  MUSI zweryfikować zgodność formy z `Phase3_nonlinear_evolution.py` #63
  i rozjazd udokumentować PRZED startem w `Phase_method_decisions.md`).
  Komórka [0, d) z periodycznymi BC, jeden obiekt na komórkę.
- **M-D (dynamika/linearyzacja):** M0-f_ε jak #63 Phase 3 (ta sama
  regularyzacja f_ε, ε=0.2; kontrola ε=0.1), zredukowana do 1D;
  linearyzacja wokół tła z M-C; warunek Blocha δg(x+d)=e^{ikd}δg(x).
- **Tła startowe relaksacji (zalockowane 2):** bump „głęboki"
  (g_min≈0.3) i „płytki" (g_min≈0.7) — gauss wycentrowany w komórce;
  wolno dodać starty, NIE usuwać.
- **Skan d (zalockowany, 5 punktów):** d ∈ {π, 2π, 3π, 4π, 6π}
  (długość fali ogona oscylacyjnego 2π przy |m²|=1); wolno dodać,
  NIE usuwać.
- **Siatka k (zalockowana):** ≥16 punktów równomiernie w [0, π/d]
  (strefa zredukowana; kontrola ×2 = 32 punkty).
- **Rejestr WEJŚĆ:** β=γ=1 (jednostki bezwymiarowe, jak Q2); ε=0.2/0.1
  [INPUT z #63]; brak kalibracji gatunkowych g₀ w 1D (amplituda tła
  wyznaczona przez relaksację — raportować, nie stroić).

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — bramka maszynerii (osiągalne FAIL-e analityczne)

- **P1a (dyspersja próżni, exact):** Bloch na próżni g≡1 w tej samej
  maszynerii: sektor tachionowy MUSI dać ω²(k)=k²−1, sektor stabilny
  (kontrola: źródło z m²=+1) MUSI dać ω²(k)=k²+1. Gate: błąd względny
  ≤1e−3 na całej siatce k przy siatce bazowej ORAZ redukcja błędu ~×4
  przy h/2 (zbieżność rzędu 2).
- **P1b (gate energii ewolucji):** kod ewolucji nieliniowej na próżni
  + mała perturbacja: |ΔE|/E ≤ 1e−6, zbieżność dt (×2).
- **Gate Phase 1 FAIL ⟹ STOP (kod nieważny).**

### Phase 2 — samouzgodnione tło łańcucha (istnienie)

- **P2a:** relaksacja (Newton lub tłumiony gradient flow — wybór zapisać
  w `Phase_method_decisions.md` PRZED startem) do ‖residuum EL‖_∞ ≤ 1e−10,
  dla 5 zalockowanych d × 2 starty; ≥2 siatki (N i 2N).
- **Kryterium istnienia:** rozwiązanie NIESTAŁE (odróżnialne od próżni:
  ‖g−1‖_∞ ≥ 0.05) i zbieżne siatkowo (‖g_N−g_2N‖_∞ ≤ 1e−4).
- **Jeśli dla ŻADNEGO d żaden start nie daje rozwiązania niestałego
  (kolaps do próżni lub g→0) ⟹ CLOSED-GATE-STOP z raportem** (wynik
  negatywny istnienia = pełnoprawne rozstrzygnięcie; zawęża hipotezę
  drabiny do 3D).

### Phase 3 — spektrum pasmowe (RACHUNEK CENTRALNY)

- **P3a:** dla każdego d z tłem z Phase 2: ω²_min(d) = min po siatce k
  z λ_min(L̂_d(k)); ≥2 siatki przestrzenne × 2 siatki k.
  **Kryterium zbieżności: |Δω²_min| ≤ 0.01·max(|ω²_min|, 0.1)**
  (identyczne jak LOCK poprzednika — porównywalność).
- **P3b (kontrola artefaktu, nieusuwalna):** to samo spektrum dla próżni
  g≡1 w superkomórce m·d (m=4): wszystkie gałęzie MUSZĄ odtworzyć
  ω²(k)=k²−1 zwinięte do strefy (identyfikacja modów próżni PRZED
  interpretacją ω²_min tła — lekcja P2c poprzednika).
- **P3c (kontrola sektora stabilnego, nieusuwalna):** ta sama procedura
  (relaksacja + Bloch) w sektorze stabilnym m²=+1 ze źródłem gaussowskim
  (konstrukcja Q2): ω²_min MUSI być > 0 — osiągalny FAIL maszynerii.

### Phase 4 (warunkowa — tylko przy zbieżnym Phase 3) — test nieliniowy

- Superkomórka 4·d, perturbacja ± wzdłuż modu minimalnego (amplituda
  0.01·‖tło‖); ewolucja M-D do t = 3·t*_ref, gdzie t*_ref = 3.62 (#63).
- Kryteria końcowe:
  - **Q-PASS:** istnieje d z ω²_min(d) > 0 (zbieżne wg P3a) ORAZ brak
    ucieczki (g→0 lub ‖δg‖ > 100% tła) do 3·t*_ref przy obu znakach.
  - **Q-FAIL:** ω²_min(d) < 0 dla WSZYSTKICH zalockowanych d (zbieżnie)
    ORAZ ucieczka w t ≤ 2·t*_ref.
  - **Q-INCONCLUSIVE:** pozostałe przypadki (w tym niezbieżność) —
    raportować wprost.

## 4. Forbidden moves

1. Zero zmian kryteriów/progów/punktów d/siatek k po starcie obliczeń;
   korekta wyłącznie dla udokumentowanego błędu implementacji wykrytego
   kontrolą (opis w `Phase*_correction_note.md` PRZED użyciem wyniku).
2. Każdy test musi mieć osiągalny FAIL; kontrole P1a, P3b, P3c
   nieusuwalne.
3. Zakaz wyboru „lepszej" siatki/startu post-hoc; oba starty raportowane.
4. Wyniki negatywne wprost; niezbieżność jako niezbieżność; zakaz
   skanowania do celu (dodatkowe d wolno dodać wyłącznie PRZED Phase 3).
5. Rdzeń `.tex` NIETKNIĘTY; wnioski → NEEDS (user-gate). NIE edytować
   STATE.md. NIE commitować (commit robi sesja główna za zgodą usera).
6. Decyzje implementacyjne (forma EL vs #63, schemat relaksacji,
   dyskretyzacja Blocha) zapisane w `Phase_method_decisions.md` PRZED
   startem obliczeń danej fazy.
7. Skrypty przed uruchomieniem; outputy do `Phase*_output.txt`.
8. **Higiena ścieżek (lekcja 2026-08-23):** NIE zmieniać cwd; po KAŻDYM
   zapisie pliku zweryfikować materializację (`ls`/`Test-Path`) —
   artefakt `<x>/TGP/TGP_v1/…` powstaje przy zapisie vault-relative
   spod cwd repo.

## 5. Deliverables

- `Phase_method_decisions.md` (przed startem), `Phase1_machinery_gate.py`,
  `Phase2_relax_chain.py`, `Phase3_bloch_spectrum.py`,
  `Phase4_nonlinear.py` (warunkowy) + outputy `Phase*_output.txt`.
- `Phase_FINAL_close.md` (werdykt wg §3), `NEEDS.md` (user-gated),
  `README.md` aktualizowany po fazach.

## 6. Drzewo decyzyjne

```text
P1 gate FAIL → STOP (kod nieważny; raport)
P2 brak istnienia dla wszystkich d → CLOSED-GATE-STOP (negatyw istnienia;
      NEEDS: hipoteza drabiny wymaga 3D — kandydat na osobny cykl)
Q-PASS → NEEDS: „stabilność = własność konfiguracji o skończonej gęstości"
      (user-gate, dopisek do retrospektywy + sek08b Limitations)
      + charakterystyka ω²(d) deskryptywnie + kandydat 3D
Q-FAIL → NEEDS: Limitations — periodyczna kąpiel samouzgodniona nie
      stabilizuje w 1D (klasa zbadana); pozostaje 3D/poza-klasowe
Q-INCONCLUSIVE → NEEDS: wniosek metodologiczny, bez ontologii
```

---

**LOCK ZAMKNIĘTY 2026-08-31. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
