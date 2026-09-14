---
title: "Phase0_balance — LOCK: kieszeń przedmetryczna na poziomie 0 (N1 op-dynamics-class-M911) — los kieszeni fazy gołej Φ=0 w próżni metrycznej, kieszeń z przestrzałem znaku (kropla antyfazowa Z2), ściana domenowa jako trwała struktura przedmetryczna, oraz czy kolaps kieszeni zostawia zlokalizowany obiekt na ścianie"
date: 2026-09-02
type: phase0-lock
tgp_owner: research/op-premetric-pocket-level0-2026-09-02
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-06: „ok, działaj z N1" (N1 z NEEDS op-dynamics-class-M911: fizyka przejścia ψ→0 — LOCK poziomu 0 dla kieszeni przedmetrycznej; ścieżka (a): substrat, bez zmiany modelu efektywnego)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-dynamics-class-M911-2026-09-02/NEEDS.md]]"
  - "[[../op-dynamics-class-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase0_balance.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
---

# Phase 0 — LOCK cyklu `op-premetric-pocket-level0`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

## 0. Pytanie i dziedziczona diagnoza

Poziom 1 (op-dynamics-class-M911, wariant bezwładny pary M9.1''):
dip słabopolowy kolapsuje do ψ→0 w skończonym, zbieżnym czasie
(t≈5.55), z przestrzałem przez zero w nierozwiązywalnym finale —
model efektywny tam się KOŃCZY (B→0, w→0). N1: rozstrzygnięcie należy
do poziomu 0.

Model poziomu 0 dziedziczony VERBATIM z `op-bare-substrate-genesis`
(CLOSED, G1–G6 PASS — „fundament mostu do grawitacji STOI"):
substrat s(x)∈ℝ na siatce 2D, Φ=s²; **Φ=0 = faza niemetryczna
(przedmetryczna), metastabilna**; prawdziwa próżnia metryczna
|s|=s*≈1.174 (Φ*≈1.3787); bariera |s_bar|≈0.426 (Φ_bar≈0.181);
dynamika: gradient flow w parametrze selekcji τ (nie czasie fizycznym).

**Fakt strukturalny odnotowany PRZED obliczeniami (analitycznie):**
substrat ma symetrię Z2 (s→−s), obserwablem jest Φ — poziom 1 (ψ) jest
ŚLEPY na znak s. Każda ciągła ścieżka między −s* a +s* przechodzi przez
s=0 ⟹ ściana domenowa Z2 ma z konieczności rdzeń Φ=0. Kolaps poziomu 1
z przestrzałem (ψ→0, potem „eksplozja") ma naturalny odpowiednik
poziomu 0: lokalne przejście s do PRZECIWNEJ studni. Stąd trzy pytania.

**Pytania binarne:**
- **Q-PKT-A (kieszeń goła, ten sam znak):** kropla fazy gołej (s≈0)
  o promieniu R wewnątrz próżni +s* — czy ZAWSZE zasklepia się
  (HEAL: front próżni zamyka kieszeń; prawo τ_close(R) deskryptywnie),
  czy istnieje R z PERSIST/FRAGMENT?
- **Q-PKT-B (kieszeń z przestrzałem znaku — RACHUNEK CENTRALNY):**
  kropla antyfazowa (rdzeń s=−s*) o promieniu R w próżni +s* — jej
  brzeg to zamknięta ściana Z2 = powłoka przedmetryczna (Φ→0 na
  brzegu). Czy (a) kurczy się i znika (lifetime τ_life(R), prawo
  skalowania — kandydat R²), (b) TRWA (pierwszy stabilny obiekt
  przedmetryczny), (c) fragmentuje?
- **Q-PKT-C (ściana jako trwała struktura przedmetryczna):** czy
  płaska ściana Z2 (geometria pasków na torusie) jest stacjonarna
  i trwała (persistence_tail ≥ 0.99), z pasem przedmetrycznym Φ<ε
  o stabilnej szerokości — tj. czy poziom 0 MA trwałe struktury ψ≈0
  (topologiczne), niewidoczne w opisie poziomu 1?
- **Q-PKT-D (pozostałość po kieszeni na ścianie):** kieszeń goła
  osadzona NA ścianie Z2 — czy po zasklepieniu zostaje zlokalizowany
  nadmiar fazy gołej ponad bazową ścianę („koralik", pre-rejestrowany
  POZYTYW autora — nośnik kreacji), czy ściana wraca do czystego
  profilu bazowego?

Uwaga interpretacyjna zalockowana: HEAL-ALL w A + zanik w B + PASS w C
to spójna odpowiedź N1: „kieszeń przedmetryczna re-metryzuje się
w skończonym τ; trwałość przedmetryczna żyje wyłącznie na topologii Z2
(ściany), której poziom 1 nie widzi". Każdy inny wynik raportować wg
litery. ZERO claimów grawitacyjnych/obserwacyjnych (poziom 0).

## 1. Model ZAMKNIĘTY (dziedziczony VERBATIM, cytat)

Z LOCKa `op-bare-substrate-genesis` §1–2 (CLOSED; cytat):
> „Phi(x) = s(x)^2 ≥ 0; ds/dtau = kappa·Lap(s) − V′(s);
> V(s) = 0.5·a·s² − (b/3)|s|³ + 0.25·c·s⁴; a=0.50, b=1.60, c=1.00;
> kappa=0.50; grid N=128, L=64.0 (dx=0.5); dt=0.02; brzegi periodyczne
> (torus), Laplasjan 5-punktowy; Phi_metric_threshold eps=0.30;
> A_min = 4 węzły/N²; bare noise U(−0.05,+0.05)."

Konsekwencje analityczne (zweryfikowane w tamtym cyklu, Phase0_output):
s*=1.174 (Φ*=1.3787), s_bar=0.426 (Φ_bar=0.181), V(s*)=−0.044,
szerokość ścianki ~√(2κ/V″(s*))≈1.07. ZAKAZ zmiany a,b,c,κ,ε,A_min.

**Rejestr WEJŚĆ tego cyklu (nowe, zalockowane TERAZ):**
- seed szumu = 20260905; szum wnętrza kieszeni U(−0.05,+0.05)
  (amplituda dziedziczona — ≪ bariera);
- promienie R ∈ {4, 8, 16} (w jednostkach L; R=16 = 0.25·L — łagodnie
  poniżej strefy boundary_contact 0.4·L);
- profil brzegowy kieszeni/kropli: tanh o szerokości δ=1.07
  (= szerokość ścianki; A: s = s*·½(1+tanh((r−R)/δ)) + szum·½(1−tanh(…));
  B: s = s*·tanh((r−R)/δ) — kropla −s* w +s*);
- geometria C/D: paski ±s* (dwie proste ściany pionowe w x=L/4, 3L/4,
  profil tanh δ=1.07); D: kieszeń (jak A, R∈{4,8}) centrowana na
  ścianie x=L/4;
- kroki: steps_max = 30000 (τ_max=600); próbkowanie co 50 kroków (Δτ=1);
- **detektor przedmetryczny:** bare_area(τ) = frakcja węzłów z Φ<ε
  (ε=0.30 dziedziczone); metric_area = frakcja Φ>ε; persistence_tail
  = bare_area(krok 30000)/bare_area(krok 27000) (ostatnie 10%,
  reguła dziedziczona);
- **mapowanie na poziom 1 (deskryptywne, obowiązkowe):** ψ ≡ Φ/Φ*
  (próżnia ψ=1); pas przedmetryczny Φ<ε ⟺ ψ<0.218;
- kontrola siatki (dziedziczona reguła pinningu G5): przypadki
  POZYTYWNE/graniczne (PERSIST, RESIDUAL, oraz C zawsze) re-run
  N=256, dx=0.25, dt=0.005 — klasyfikacja musi się jakościowo zgadzać;
- flaga boundary_contact (Chebyshev > 0.4·L od centrum) — dziedziczona;
  dotknięcie = dalsze kroki poza ocenami;
- 1D profil ściany (pomocniczy do C): N=4096, dx=0.05, relaksacja
  z BC ±s* do ‖ds/dτ‖∞ ≤ 1e−10; wielkości: Φ_min, szerokości pasów
  Φ<ε i Φ<Φ_bar, napięcie σ_w = H[ściana]−H[próżnia]; kontrola
  dx/2 — zmiana Φ_min i szerokości < 1%.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — bramka ciągłości z poprzednikiem (FAIL ⟹ STOP)
- P1a: `bare` (szum U(−0.05,0.05), seed 20260905, 6000 kroków):
  metric_area < A_min (reprodukcja G2 klasy);
- P1b: `single(A0=1.4, w=1.5)` po 6000 kroków: metric_area < A_min
  (reprodukcja G3: podkrytyczny single zanika; A0=1.4 = największy
  podkrytyczny z tamtego skanu);
- P1c: H_Γ(τ) nierosnące we wszystkich biegach P1 (guard przepływu);
- P1d: 1D ściana: zbieżność dx→dx/2 (<1% w Φ_min i szerokościach).

### Phase 2 — rachunki (A, B, C, D)
- **Q-PKT-A:** HEAL(R) ⟺ bare_area < A_min przed τ_max i pozostaje
  (wszystkie późniejsze próbki); τ_close(R) = pierwsze τ z bare_area
  < A_min; deskryptywnie fit τ_close vs R (liniowy vs kwadratowy — R²
  raportować, bez progu). PERSIST(R) ⟺ bare_area ≥ A_min z
  persistence_tail ≥ 0.99 (wtedy kontrola N=256 obowiązkowa).
  **Werdykt:** Q-PKT-A-HEAL-ALL / Q-PKT-A-PERSIST(R…) / INCONCLUSIVE
  (boundary_contact / niezbieżność kontroli).
- **Q-PKT-B (centralny):** τ_life(R) = pierwsze τ, w którym znika
  region s<−s_bar (rdzeń antyfazowy; równoważnie max(−s) < s_bar).
  ZANIK ⟺ τ_life(R) < τ_max dla wszystkich R (prawo skalowania
  τ_life vs R deskryptywnie: liniowe vs R² — R² fitów raportować);
  TRWA ⟺ dla pewnego R rdzeń antyfazowy przetrwa τ_max
  z persistence_tail(bare_area) ≥ 0.99 — wtedy kontrola N=256
  obowiązkowa; PASS-OBJECT tylko przy zgodnej kontroli. Deskryptywnie:
  bare_area(τ) powłoki, Φ_min(τ), kształt końcowy.
- **Q-PKT-C:** paski bez kieszeni, 30000 kroków: PASS-SHEET ⟺
  bare_area pasa ścian stabilna (persistence_tail ≥ 0.99) i profil
  poprzeczny zgodny z 1D (Φ_min < Φ_bar w rdzeniu obu ścian)
  + kontrola N=256 zgodna. FAIL ⟺ ściany znikają/dryfują do anihilacji
  przed τ_max (raportować mechanizm).
- **Q-PKT-D:** stan końcowy `pocketD(R)` vs bazowa ściana (C):
  excess_area = bare_area_final − bare_area_bazowej_ściany.
  RESIDUAL-OBJECT ⟺ excess_area ≥ A_min ORAZ persistence_tail
  nadmiaru ≥ 0.99 ORAZ kontrola N=256 zgodna (pre-rejestrowany POZYTYW
  autora); CLEAN-WALL ⟺ excess_area < A_min; inne = INCONCLUSIVE.

### Phase 3 — zamknięcie
`Phase_FINAL_close.md` (werdykty + mapowanie ψ=Φ/Φ* + odpowiedź N1
wprost), `NEEDS.md` (user-gated), `README.md`.

## 3. Forbidden moves
1. Zmiana a,b,c,κ,ε,A_min,dx,dt po starcie obliczeń — zakaz (formy
   i progi dziedziczone z CLOSED cyklu, cytat §1).
2. Zmiana kryteriów/promieni/seedów/definicji obszarów po pierwszym
   biegu — zakaz.
3. Zero reguł spawn/przeżycia; los struktur wyłącznie z pola s
   (dziedziczone).
4. Rdzeń .tex/STATE/git nietykane; katalogi innych cykli tylko odczyt.
5. Zero claimów grawitacyjnych/obserwacyjnych; τ ≠ czas fizyczny;
   wynik negatywny raportowany wprost.
6. INCONCLUSIVE ≠ pozytyw; kontrola N=256 nieusuwalna dla pozytywów.

## 4. Deliverables
`Phase1_gate.py`+output, `Phase2_pocket.py`+output(+npz stanów),
`Phase_FINAL_close.md`, `NEEDS.md`, `README.md`.

## 5. Drzewo decyzyjne
```text
P1 FAIL → STOP (brak ciągłości z zamkniętym cyklem — problem maszynerii)
A-HEAL-ALL ∧ B-ZANIK ∧ C-PASS-SHEET → odpowiedź N1: kieszeń re-metryzuje
   się w skończonym τ (bounce), przestrzał znaku daje POWŁOKĘ przedmetryczną
   o życiu τ_life(R) (długowieczną przy dużym R), trwałość przedmetryczna
   = topologia Z2 (ściany) niewidoczna dla poziomu 1; NEEDS: dopisek core
   (user-gate) + pytanie o 3D/sieć ścian
B-TRWA (obiekt antyfazowy stabilny) → PIERWSZY trwały obiekt przedmetryczny
   programu; NEEDS PILNE: charakterystyka + odpowiednik 3D + user-gate core
D-RESIDUAL-OBJECT → nośnik hipotezy usera (kolaps kieszeni ZOSTAWIA obiekt);
   NEEDS: mechanizm + wersja 3D
C-FAIL (ściany nietrwałe) → trwałość przedmetryczna nie istnieje nawet
   topologicznie w 2D — kieszenie ψ→0 są zawsze przejściowe; raport wprost
INCONCLUSIVE → wniosek metodologiczny (pinning/boundary), ewentualny re-lock
```

---

**LOCK ZAMKNIĘTY. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move.**
