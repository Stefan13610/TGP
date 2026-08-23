---
title: "Phase0_balance — LOCK: runaway w kąpieli na fazach ZMIERZONYCH (Q1) + derywacyjny test hipotezy dwóch sektorów: znak efektywnego W z gęstości (Q2)"
date: 2026-08-23
type: phase0-lock
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-23: 'ok, przygotuj prompt' (N3 z NEEDS cyklu op-lattice-bath-runaway + hipoteza autora z N2: 'różnica znaku może wynikać z tego czy mierzymy pojedynczy obiekt czy obiekt w kąpieli sąsiadów'); realizacja Phase 1–3: nowa sesja, osobny agent"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/PhaseA_report.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N1_pochodzenie-faz_2026-08-23.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
---

# Phase 0 — LOCK cyklu `op-bath-two-sectors`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria, modele,
progi i forbidden moves zamknięte PRZED kodem.**

---

## 0. Dwa pytania binarne

- **Q1 (rachunek centralny, N3):** czy mod runaway solitonu (potwierdzony
  w izolacji: #63 V3, g→0, t*≈1.7–3.6) dostaje **ω² > 0 przy jakiejś
  skończonej gęstości sąsiadów n** — z ogonami i fazami **ZMIERZONYMI**
  z faktycznych układów (nie z dokumentacji rdzenia, która padła w A3)?
- **Q2 (hipoteza dwóch sektorów, autor 2026-08-23):** czy z akcji gałęzi
  STABILNEJ (eq:field-eq-reproduced, m²=+γ, W″(1)>0) efektywna krzywizna
  potencjału fluktuacji **wokół tła o skończonej gęstości źródeł n** zmienia
  znak na tachionowy? Jeśli TAK — rozdwojenie znaku W z N2 nie jest
  konfliktem dwóch akcji, tylko **dwoma sektorami jednej akcji**
  (izolowany obiekt w próżni vs obiekt w kąpieli), a maszyneria 2 dostaje
  derywacyjną kotwicę. Jeśli NIE — wybór znaku pozostaje otwartym
  problemem aksjomatycznym (do decyzji ontologicznej autora).

Q1 i Q2 są niezależne logicznie; oba wyniki negatywne = pełnoprawne
rozstrzygnięcia (raportować wprost).

## 1. Kontekst (dziedziczony — NIE powtarzać tych rachunków)

- Bramka audytu maszynerii 2 (cykl [[../op-lattice-bath-runaway-2026-08-23/README.md]],
  CLOSED-GATE-STOP): **A1 PASS** (próg 8/5, 3.7e−11), **A2 PASS** (ω=1
  uniwersalne), **A3 FAIL** (fazy dokumentacji niereprodukowalne z układu
  potęgowego), **A5** (rozdwojenie = znak W).
- **N1:** fazy dodatekH żyją w układzie LOGARYTMICZNYM eq:J-ode
  (F=1+2α_eff·ln g; Δ(e→μ)=120.01° reprodukowane co do setnej) — nie
  w potęgowym.
- **N2:** rdzeń ma OBA znaki W (sek08a: prop:field-eq-from-action vs Nota
  kanoniczna 2026-04-07, której EL = układ potęgowy M2).
- Model referencyjny locka oscylacyjnego: `RETRO_oscillating_tail_lock.py`
  (4/4 PASS; drabina minimów co 2π, kontrola negatywna Yukawa).

## 2. Modele ZAMKNIĘTE

- **M-P (potęgowy, PRIMARY dla Q1):**
  `g″ + (2/r)g′ + (2/g)g′² = g²(1−g)` — EL Noty kanonicznej (K=g⁴);
  audytowany (A1/A2 PASS). Obiekt: solitony g₀<8/5; gatunki: g₀_e=0.90548,
  g₀_μ=φ·g₀_e (kalibracja jak w bramce A); τ: brak regularnego
  reprezentanta blisko progu — dopuszczone g₀_τ=1.5696 (z Q_K=3/2,
  **flagować INPUT**) z kontrolą wrażliwości ±2% progu.
- **M-L (logarytmiczny, CROSS-CHECK dla Q1):** eq:J-ode z runningiem:
  `F(g)g″+(2/r)g′=g²(1−g)`, F=1+2α_eff·ln g, α_eff=2/(1+η_K(g−1)²),
  η_K=181/15; g₀_e=0.90548, g₀_μ=φ·g₀_e (nośnik faz Δ=120.01° — N1).
- **M-B (kąpiel dla Q1):** komórka Wignera–Seitza / klatka periodyczna;
  tło = soliton centralny + suma ogonów sąsiadów A·cos(r−δ)/r
  z (A, δ) per gatunek WYŁĄCZNIE z Phase 1 (κ≈0 dopuszczone — próżnia
  tachionowa nie tłumi wykładniczo; jeśli fit wykryje tłumienie, użyć
  zmierzonego). Konfiguracja jednorodna jednogatunkowa (ee, μμ) +
  jedna mieszana (eμ) — 3 konfiguracje zalockowane.
- **M-S (sektor stabilny dla Q2):** eq:field-eq-reproduced (sek08a):
  `∇²ψ + 2(∇ψ)²/ψ + βψ² − γψ³ = −qρ̃` (β=γ=1, jednostki bezwymiarowe);
  źródła ρ̃ = periodyczna sieć profili gaussowskich o szerokości σ_s
  i separacji d (gęstość n~1/d³ w 3D radialnie / 1/d w komórce 1D —
  wariant wymiarowy rozstrzygnąć i zapisać w raporcie PRZED startem
  obliczeń Q2, potem NIE zmieniać). Pytanie: znak minimalnej wartości
  własnej operatora fluktuacji wokół samouzgodnionego tła ψ_n(x)
  jako funkcja d.
- **Rejestr WEJŚĆ:** Q_K=3/2 (g₀_τ), r₂₁/φ-FP (kalibracja g₀_e/g₀_μ),
  η_K=181/15 — wszystkie flagować INPUT w każdym zależnym wyniku.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — ekstrakcja zmierzonych (A, δ) i drabina d*

- **P1a:** dla M-P i M-L, gatunki e/μ(/τ w M-P): fit ogona
  (g−1)·r ~ B·cos r + C·sin r; **dwa zalockowane okna [50,150]
  i [120,260]** (rozjazd okien raportować; >3° = flaga WINDOW-SENSITIVE),
  R-kontrola (R=300 vs 450, |Δδ|≤1% zakresu), konwencja δ=atan2(C,B).
- **P1b:** E_int(d) par {ee, eμ, μμ} z sumy parowej ogonów; pierwsze
  minimum d* i odstęp drabiny per para per model. Predykcja referencyjna
  (zalockowana TERAZ): odstęp = 2π ± 5% (korekta 1/d malejąca).
- **P1c (kontrola negatywna, nieusuwalna):** profil Yukawa bez cos
  (te same A): **0 minimów**. Minimum w kontroli ⟹ STOP i debug PRZED
  interpretacją.
- **Gate Phase 1:** fity zbieżne (R² ≥ 0.999) i R-stabilne; kontrola czysta.

### Phase 2 — Q1: runaway w kąpieli (RACHUNEK CENTRALNY)

- **P2a (baseline, gate):** reprodukcja runaway w izolacji (#63 V3) na tym
  samym kodzie/siatce co kąpiel; **gate |ΔE|/E ≤ 1e−6**, zbieżność dt (×2).
  FAIL ⟹ STOP (kod nieważny).
- **P2b (skan):** n z drabiny Phase 1: d ∈ {d*₁, d*₂, d*₃} ∪ {0.5·d*₁,
  1.5·d*₁} (5 punktów minimum per konfiguracja; wolno dodać, NIE usuwać).
  Dla każdego: (i) ω²_min spektrum linearyzacji wokół konfiguracji kąpieli
  (≥2 siatki × ≥2 rozmiary komórki); (ii) ewolucja nieliniowa
  z perturbacją ±wzdłuż modu runaway.
- **Kryteria Q1:**
  - **Q1-PASS:** istnieje n z ω²_min>0 (zbieżne do 2 cyfr, stabilne na
    rozmiarach komórki) ORAZ brak runaway do t=3·t*_izolacji przy ±.
  - **Q1-FAIL:** ω²_min<0 dla wszystkich zalockowanych n ORAZ ucieczka
    w t ≤ 2·t*_izolacji.
  - **Q1-INCONCLUSIVE:** niezbieżność siatek/komórek — raportować JAKO
    niezbieżność; zakaz wyboru „lepszej" siatki post-hoc.
- **P2c (kontrola artefaktu komórki, nieusuwalna):** spektrum próżni
  jednorodnej w tej samej komórce; mody komórki zidentyfikowane PRZED
  interpretacją ω²_min (lekcja N_neg=floor(R/π)).

### Phase 3 — Q2: derywacyjny test dwóch sektorów

- **P3a:** sektor M-S: samouzgodnione tło ψ_n dla sieci źródeł przy
  d ∈ {∞ (kontrola: pojedyncze źródło), 8, 6, 4, 2} (5 punktów
  zalockowanych; skala σ_s zapisana przed startem). Kontrola d=∞ MUSI
  odtworzyć Yukawę m²=+γ (to jest osiągalny FAIL procedury).
- **P3b:** operator fluktuacji L̂_n wokół ψ_n: znak krzywizny efektywnej —
  ω²_min(d) na ≥2 siatkach. Dodatkowo statyczna odpowiedź na punktowe
  zaburzenie tła ψ_n: monotoniczna (Yukawa) czy oscylacyjna.
- **Kryteria Q2:**
  - **Q2-PASS (hipoteza autora potwierdzona):** istnieje d z ω²_min<0
    lub odpowiedzią oscylacyjną, zbieżnie, przy czystej kontroli d=∞ —
    znak tachionowy EMERGUJE z gęstości w akcji stabilnej.
  - **Q2-FAIL:** ω²_min>0 i odpowiedź monotoniczna dla wszystkich d —
    gęstość nie zmienia znaku; rozdwojenie W pozostaje problemem
    aksjomatycznym.
  - **Q2-INCONCLUSIVE:** jak wyżej, raportować.
- Q2 wykonywać NIEZALEŻNIE od werdyktu Q1 (nie jest warunkowa).

## 4. Forbidden moves

1. Zero zmian kryteriów/progów/list punktów/okien po starcie obliczeń;
   korekta wyłącznie dla udokumentowanego błędu implementacji wykrytego
   kontrolą (opis PRZED użyciem wyniku).
2. Każdy test musi mieć osiągalny FAIL; kontrole P1c, P2c, P3a(d=∞)
   nieusuwalne.
3. (A, δ) w Phase 2 wyłącznie z Phase 1; parametry Q2 wyłącznie z M-S
   (zakaz przenoszenia czegokolwiek z gałęzi tachionowej do Q2 —
   Q2 ma sens tylko jako derywacja z akcji STABILNEJ).
4. Wyniki negatywne wprost; niezbieżność jako niezbieżność; zakaz
   skanowania do celu.
5. Rdzeń `.tex` NIETKNIĘTY; wnioski → NEEDS (user-gate). NIE edytować
   STATE.md. NIE commitować bez zgody.
6. Q_K=3/2, η_K, kalibracje — flagować INPUT.
7. Skrypty przed uruchomieniem; outputy do `Phase*_output.txt`.
8. **Higiena ścieżek (lekcja 2026-08-23):** NIE zmieniać cwd (`cd`/
   `Set-Location`); po KAŻDYM zapisie pliku zweryfikować materializację
   (Test-Path/ls) — artefakt `<x>/TGP/TGP_v1/…` powstaje przy zapisie
   vault-relative spod cwd repo.

## 5. Deliverables

- `Phase1_tail_measured.py` (+ output) — (A, δ) per gatunek per model,
  d* per para; `Phase2_bath_runaway.py` (+ output) — baseline + skan Q1;
  `Phase3_two_sectors.py` (+ output) — tła ψ_n + spektra Q2.
- `Phase_FINAL_close.md` (werdykty Q1 i Q2 wg §3), `NEEDS.md` (user-gated),
  README aktualizowany po fazach.

## 6. Drzewo decyzyjne

```text
P1 gate FAIL → STOP (raport; bez Phase 2; Phase 3/Q2 WOLNO wykonać —
               nie zależy od ogonów Phase 1)
P2a FAIL → STOP Q1 (kod nieważny); Q2 kontynuować
Q1-PASS → NEEDS: dopisek „stabilność = własność konfiguracji o skończonej
          gęstości" (user-gate) + charakterystyka ω²(n) deskryptywnie
Q1-FAIL → NEEDS: Limitations — kąpiel oscylacyjna nie stabilizuje
          w klasie zbadanej (ostatnia kryjówka ontologii policzona)
Q2-PASS → NEEDS: rozdwojenie znaku W = dwa sektory jednej akcji
          (derywacyjna kotwica maszynerii 2; wpis do sek08a — user-gate)
Q2-FAIL → NEEDS: znak W = otwarty problem aksjomatyczny (decyzja
          ontologiczna autora, nie numeryka)
```

---

**LOCK ZAMKNIĘTY 2026-08-23. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
