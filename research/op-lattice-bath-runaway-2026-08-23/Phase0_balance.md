---
title: "Phase0_balance — LOCK: audyt maszynerii 2 (bramka) + test runaway solitonu w kąpieli sąsiadów (rachunek centralny retrospektywy 2026-08-16)"
date: 2026-08-23
type: phase0-lock
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-23: 'tak zrób też audyt i start z fazą 0' (audyt maszynerii 2 WEWNĄTRZ cyklu jako bramka; realizacja Phase A–3: następna sesja / osobny agent)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/AUDYT_trzy-rezimy_ZAMKNIECIE_2026-08-10.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
  - "[[../../core/_meta_latex/status_map.tex]]"
  - "[[../../core/formalizm/dodatekH_lancuch_wyprowadzen.tex]]"
---

# Phase 0 — LOCK cyklu `op-lattice-bath-runaway`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Wszystkie kryteria,
progi, modele i forbidden moves zamknięte PRZED napisaniem jakiegokolwiek kodu.**

---

## 0. Cel i pytanie binarne

**Pytanie centralne (retrospektywa 2026-08-16, §6 pkt 1):**

> Czy mod runaway solitonu (potwierdzony w izolacji: #63 V3, g→0, t*≈1.7–3.6)
> **dostaje ω² > 0 przy jakiejś skończonej gęstości sąsiadów n** —
> tj. czy kąpiel sąsiadów z ogonami oscylacyjnymi (z faktycznego ODE rdzenia,
> nie z modelu) stabilizuje konfigurację, którą ontologia TGP uznaje za fizyczną?

To jest pierwszy w korpusie test stabilności konfiguracji o skończonej gęstości
źródeł (ślepa plamka udokumentowana: wszystkie testy #60–#63 + audyty = pojedynczy
obiekt/para w próżni). Wynik pozytywny = mechanizm stabilności ontologii balonów
policzony; wynik negatywny = ontologia „ciśnienia sąsiadów" traci ostatnią
niepoliczoną kryjówkę. **Oba wyniki są pełnoprawnymi rozstrzygnięciami.**

**Warunek konieczny (bramka):** maszyneria 2 (ODE / O-L5 / dodatekH), z której
bierzemy ogony, NIE przeszła audytu takiego jak maszyneria 1 (EFT Φ — padła
2026-08-10). Retrospektywa §5: „Nie budować na niej bez audytu." Dlatego
**Phase A = audyt maszynerii 2 jest bramką całego cyklu**.

## 1. Wejście dla agenta (kolejność czytania)

1. [[README.md]] tego cyklu (skrót) → ten LOCK w całości.
2. [[../op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]] + `RETRO_oscillating_tail_lock.py` (model referencyjny locka).
3. `core/_meta_latex/status_map.tex` lin. 1283–1316 (O-L5: próg kolapsu, ogon oscylacyjny, status Q_K).
4. `core/formalizm/dodatekH_lancuch_wyprowadzen.tex` lin. 1110–1161 + p131–p146 (fazy ogona, KOREKTA p134e–g, K-samodzielność).
5. [[../op-native-pressure-lepton-stability-2026-07-27/AUDYT_KRYTYCZNY_2026-07-28.md]] Dodatek A (runaway F-A kanonicznej — niespójność do rozstrzygnięcia w A5).
6. [[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]] + `Phase3_output.txt` (#63 V3 — kod do reuse: układ hamiltonowski semi-dyskretny, gate energii).
7. [[../op-native-pressure-lepton-stability-2026-07-27/AUDYT_trzy-rezimy_ZAMKNIECIE_2026-08-10.md]] (anty-wzorce audytowe: „17 PASS bez możliwego FAIL").
8. Skrypty rdzenia maszynerii 2 (do audytu w A4): `tooling/scripts/ngen_collapse_proof_v47b.py`, `gcrit_energy_proof_v47b.py`, `gcrit_pohozaev_v47b.py`, `atail_asymptotic_v47b.py`, `atail_functional_v47b.py`, `ode_koide_formA_exact_v47b.py`, `collapse_exponent_v47b.py`, `a3d_soliton_brannen_r.py`.

## 2. Modele ZAMKNIĘTE

- **M2 (maszyneria 2, obiekt audytu):** radialne ODE kanonicznej Form A
  (α=2, źródło g²(1−g)), zmienna kanoniczna f=g^(α+1); twierdzenia rdzenia
  do weryfikacji: g₀,crit = 2(α+2)/[2(α+2)−d] = 8/5 (d=3, α=2);
  ogon: linearyzacja → r. Bessela sferycznego, ω=1 uniwersalne;
  A_g ≈ |g₀−1| przy próżni, A_g ~ (8/5−g₀)^(−0.23) przy kolapsie;
  fazy δ_e=−81.4°, δ_μ=+38.6°, δ_τ(g₀=4)=−27.3°, Δ(e→μ)=120.01°≈2π/3.
- **MB (kąpiel):** komórka periodyczna (1D-radialna klatka Wignera–Seitza
  z warunkiem periodycznym LUB sieć 3D zredukowana do komórki z symetrią) —
  wybór wariantu wolno rozstrzygnąć w Phase A na podstawie audytu, ale
  MUSI być zapisany w raporcie A PRZED startem Phase 2 i potem NIE zmieniany.
  Tło = superpozycja solitonu centralnego + ogonów sąsiadów z parametrami
  (κ, φ, A) wyłącznie z Phase 1 (zakaz parametrów „z ręki").
- **MD (dynamika):** dokładny układ hamiltonowski semi-dyskretny wg wzorca
  #63 V3 (energia zachowana exact w ODE ⟹ gate mierzy czysto błąd całkowania).
- **Rejestr WEJŚĆ (jawny):** Q_K = 3/2 jest parametrem wejściowym TGP
  (status_map: 12 ścieżek wyprowadzenia zawiodło) — każdy wynik cyklu,
  który od niego zależy, MUSI to flagować w raporcie.

## 3. Fazy i kryteria PASS/FAIL (zalockowane)

### Phase A — AUDYT MASZYNERII 2 (BRAMKA CYKLU)

Wzorzec: audyt maszynerii 1 (2026-08-10). Niezależna re-implementacja,
nie odpalenie skryptów rdzenia „na wiarę".

- **A1 (próg kolapsu):** niezależna implementacja ODE (inna metoda całkowania
  niż w skryptach rdzenia) reprodukuje g₀,crit = 8/5 z |Δ| ≤ 1e−6 oraz
  formułę 2(α+2)/[2(α+2)−d] na ≥3 parach (α,d) ≠ (2,3).
  **FAIL:** rozbieżność > 1e−3 lub formuła nie generalizuje.
- **A2 (ogon oscylacyjny):** z niezależnego rozwiązania: ω_tail = 1.000 ± 0.01
  dla ≥3 wartości α; ogon oscyluje (nie zanika wykładniczo bez oscylacji)
  w oknie r gdzie |g−1| ≤ 0.2 (reżim, w którym funkcjonał wg audytu
  dwa-zrodla JEST wiarygodny). **FAIL:** ω zależne od α poza tolerancją
  lub brak oscylacji w oknie ważności.
- **A3 (fazy i amplitudy):** niezależny fit ogona reprodukuje
  δ_e=−81.4°, δ_μ=+38.6° (tolerancja ±2°), Δ(e→μ)=120° ± 1°;
  A_g ≈ |g₀−1| (współczynnik 1 ± 10%) w reżimie przy próżni.
  Okno fitu ODLEGŁE od ściany siatki: koniec okna ≤ 0.75·R; kontrola
  R∈{R₀, 1.5·R₀}: parametry stabilne do 3 cyfr. **FAIL:** fazy poza
  tolerancją lub R-zależność parametrów > 1%.
- **A4 (audyt skryptów rdzenia):** dla 8 skryptów z §1 pkt 8: (i) czy każdy
  ma test zdolny dać FAIL (kontrola negatywna / próg, który realne dane mogą
  złamać); (ii) czy SUMMARY zgadza się z faktycznym outputem; (iii) rejestr
  wejść vs wyjść (co jest INPUT: Q_K=3/2, K=g⁴, φ-FP; co faktycznie DERIVED).
  Wynik = tabela per skrypt. **To jest raport deskryptywny — nie ma PASS/FAIL,
  ma OBOWIĄZEK kompletności** (wszystkie 8, bez pomijania).
- **A5 (rozstrzygnięcie niespójności istnienia):** AUDYT_KRYTYCZNY 2026-07-28
  Dodatek A: F-A kanoniczna (u″+2u′/r = u³−u²−2u′²/u) — runaway dla wszystkich
  g₀. DodatekH/O-L5: solitony F-A kanonicznej istnieją, τ przy g₀=1.5696.
  Obowiązek: wskazać PRECYZYJNIE (równanie po równaniu) czym różnią się te
  dwa układy (zmienna? źródło? warunek brzegowy? odbicie ad-hoc?) i który
  opisuje obiekt, którego ogon bierzemy do Phase 1. **FAIL:** różnica nie
  daje się zlokalizować (⟹ maszyneria 2 wewnętrznie niespójna).
- **A6 (korekta p134e–g, uczciwość):** raport MUSI powtórzyć jawnie:
  A_τ monotonicznie rosnąca, r₃₁=3477.5 to kres skanu (g₀≤4), nie fizyczne
  maksimum; Δ(μ→τ) nie osiąga ±120°. Zakaz cytowania r₃₁ jako „wyprowadzonego".

**Reguła bramki:** A1–A3 PASS ∧ A5 rozstrzygnięte → Phase 1.
Którykolwiek z A1/A2/A3/A5 FAIL → **STOP CAŁEGO CYKLU** przed Phase 1;
raport negatywny wprost (maszyneria 2 pada jak maszyneria 1); NEEDS
z propozycją dalszej drogi (user-gate). A4/A6 są obowiązkowe niezależnie
od werdyktu bramki.

### Phase 1 — ekstrakcja (κ, φ, A) i drabina d* dla realnych par

- **P1a:** parametry ogona (κ, φ, A) per gatunek {e, μ, τ} z FAKTYCZNEGO
  rozwiązania ODE (audytowanego w A). Kontrole jak w A3 (okno, R-kontrola).
- **P1b:** E_int(d) dla par {ee, eμ, μτ, eτ} z sumy parowej ogonów;
  wyznaczyć pierwsze minimum d* i odstęp drabiny per para.
  Predykcja modelu referencyjnego (RETRO, zalockowana TERAZ): odstęp = 2π·r_core
  (tolerancja ±5%, korekta 1/d malejąca z d); d*(model, κ=0–0.5) ≈ 6.0–6.1.
  **Rozbieżność ≠ FAIL** — różne fazy gatunków mogą przesuwać d*
  (pre-rejestrowane w retrospektywie §5); rozbieżność RAPORTOWAĆ z liczbami.
- **P1c (kontrola negatywna, obowiązkowa):** ogon Yukawa bez cos (te same κ, A):
  0 minimów. Jeśli kontrola daje minima ⟹ błąd procedury ⟹ STOP i debug
  PRZED interpretacją.

### Phase 2 — RACHUNEK CENTRALNY: runaway w kąpieli

- **P2a (baseline):** reprodukcja wyniku #63 V3 w izolacji (runaway, wyjście
  z dziedziny g→0) na tym samym kodzie, którym potem liczymy kąpiel —
  ta sama siatka, ten sam gate. **Gate energii: |ΔE|/E ≤ 1e−6** (wzorzec #63);
  zbieżność dt (dwa kroki, stosunek 2). Baseline FAIL ⟹ STOP (kod nieważny).
- **P2b (skan gęstości):** kąpiel MB z n odpowiadającym odległościom
  d ∈ {d*₁, d*₂, d*₃ z Phase 1} ∪ {0.5·d*₁, 1.5·d*₁} (5 punktów minimum;
  wolno dodać punkty, NIE wolno usuwać zalockowanych). Dla każdego n:
  (i) spektrum linearyzacji wokół konfiguracji kąpieli (mod najniższy ω²)
  na ≥2 siatkach i ≥2 rozmiarach komórki; (ii) ewolucja nieliniowa
  z perturbacją wzdłuż modu runaway (obie amplitudy ±, jak #63 V3).
- **Kryteria (zalockowane):**
  - **V-PASS (stabilizacja):** istnieje n z ω²_min > 0 (zbieżne na siatkach
    do 2 cyfr, stabilne na rozmiarach komórki) ORAZ ewolucja nieliniowa
    bez runaway do t = 3·t*_izolacji przy obu znakach perturbacji.
  - **V-FAIL (brak stabilizacji):** ω²_min < 0 dla WSZYSTKICH zalockowanych n
    ORAZ ewolucja ucieka z dziedziny w t ≤ 2·t*_izolacji.
  - **V-INCONCLUSIVE:** wyniki siatek/komórek niezbieżne — raportować JAKO
    niezbieżność (zakaz wyboru „lepszej" siatki post-hoc).
- **P2c (kontrola artefaktu komórki):** próżnia jednorodna w tej samej komórce
  periodycznej — policzone spektrum; mody komórki odjęte/zidentyfikowane
  PRZED interpretacją ω²_min (lekcja N_neg=floor(R/π) z AUDYT_KRYTYCZNY).

### Phase 3 (warunkowa — wyłącznie po V-PASS)

- ω²(n) deskryptywnie; porównanie n_stab z drabiną d* z Phase 1; zależność
  od pary gatunków. Zero claimów interpretacyjnych poza deskrypcją
  (wzorzec W3b z #62: SPECULATIVE, deskryptywnie).

## 4. Forbidden moves (zalockowane)

1. Zakaz zmiany kryteriów, progów, listy n, okien fitu i konwencji po starcie
   obliczeń. Jedyna dopuszczalna korekta: udokumentowany błąd IMPLEMENTACJI
   wykryty przez kontrolę negatywną/baseline — korekta opisana PRZED użyciem
   wyniku (wzorzec T3 z RETRO: pierwotny FAIL udokumentowany w skrypcie).
2. **Każdy test MUSI mieć osiągalny FAIL** (anty-wzorzec „17 PASS bez
   możliwego FAIL" z audytu maszynerii 1). Kontrole negatywne: P1c, P2c —
   obowiązkowe, nieusuwalne.
3. Zakaz parametrów ogona „z ręki" lub z modelu RETRO w Phase 2 —
   wyłącznie z Phase 1 (z faktycznego ODE).
4. Wynik negatywny (V-FAIL, FAIL bramki A) zgłaszany wprost z liczbami
   i zbieżnością — bez reinterpretacji w tym samym raporcie.
5. Zakaz skanowania czegokolwiek „do trafienia w cel" (anty-wzorzec
   „pressure+loops=111%"). Skan n jest zalockowany w P2b.
6. Rdzeń `.tex` NIETKNIĘTY w tym cyklu; wnioski → NEEDS (user-gate).
7. NIE commitować bez zgody użytkownika.
8. Niezbieżności raportowane JAKO niezbieżności (V-INCONCLUSIVE istnieje
   po to, żeby nie było pokusy rozstrzygania na siłę).
9. Skrypty pisane PRZED uruchomieniem; outputy zapisywane do plików
   `Phase*_output.txt` w cyklu.

## 5. Deliverables

- `PhaseA_audit_machinery2.py` (+ `PhaseA_output.txt`) + `PhaseA_report.md`
  (w tym tabela A4 per skrypt i rozstrzygnięcie A5).
- `Phase1_tail_extraction.py` (+ output) — (κ, φ, A) per gatunek + d* per para.
- `Phase2_bath_runaway.py` (+ output) — baseline, skan n, spektra, ewolucje.
- `Phase_FINAL_close.md` — werdykt względem reguły bramki i V-PASS/V-FAIL.
- `NEEDS.md` — user-gated (dopiski do rdzenia zależnie od werdyktu).
- README cyklu aktualizowany po każdej fazie.

## 6. Reguła decyzyjna (drzewo, zalockowane)

```text
A1–A3, A5:  FAIL któregokolwiek → STOP; raport „maszyneria 2 pada";
            NEEDS: decyzja użytkownika (audyt głębszy / porzucenie gałęzi ODE)
A PASS →
P2a baseline FAIL → STOP (kod nieważny; debug ≠ obliczenia cyklu)
P2a PASS →
  V-PASS  → Phase 3 (deskryptywnie) + NEEDS (dopisek sek08b/status_map:
            „stabilność = własność konfiguracji o skończonej gęstości" — user-gate)
  V-FAIL  → raport wprost + NEEDS (dopisek Limitations: ciśnienie sąsiadów
            z ogonem oscylacyjnym NIE stabilizuje — ostatnia kryjówka
            ontologii policzona negatywnie — user-gate)
  V-INCONCLUSIVE → raport niezbieżności + propozycja następnego kroku
            metodycznego (bez wykonywania go w tym cyklu)
```

---

**LOCK ZAMKNIĘTY 2026-08-23. Od tego momentu zmiany w tym pliku są
forbidden move (poza dopisaniem daty realizacji faz).**
