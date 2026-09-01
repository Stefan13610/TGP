---
title: "Phase_FINAL_close — zamknięcie: Q1-POS (μ(sol|próżnia)=−0.179<0 — próżnia już w reżimie opłacalnej kreacji; znaki mieszane w policzonym zbiorze), Q2-INCONCLUSIVE (18/18 załamanie NIE-nukleacyjne: ucieczka g→+∞, podłoga QB-2 nigdy nieaktywowana; P2a 12/12 PASS, P2c 18/18 zero fałszywych alarmów), Q3 NIE WYKONANE (warunek Q2-PASS-STATIC niespełniony)"
date: 2026-09-01
type: phase-final-close
tgp_owner: research/op-metametric-boundary-2026-09-01
status: CLOSED
verdict: "Q1: POS wg litery — w policzonym zbiorze konfiguracja z μ=0 oddziela reżimy: ΔE_create(sol μ|próżnia) = −0.178999 (R=60, h=0.0125; zbieżność 4.1e−5) i ε(2π) = −0.000519 (N48; N32: −0.000555) UJEMNE vs ΔE_create(sol|stan pusty) = +16156.583 DODATNIE (zdominowane objętościowym U(1)·V(R)≈16156.8); P1a sympy 4/4 (U(0)=0<U(1)=1/56, U″(1)=−1, U′=g⁶(1−g)>0 na (0,1) — brak minimum). Q2: INCONCLUSIVE wg litery — WSZYSTKIE 18 biegów tachionowych (3 starty × 3 podłogi × 2 siatki) = BREAKDOWN przez skończono-czasową ucieczkę g→+∞ (U=g⁷/7−g⁸/8 → −∞): sol t≈2.75/3.13, sieć 2π t≈1.25/1.49, próżnia+szum t=11.06 (identyczny na obu siatkach — wspólne pole pasmowe); min g NIE osiągnęło g_floor w ŻADNYM biegu (0.6065/0.7773 bez spadku; szum: całe pole poszło W GÓRĘ, min g→1.85); ZERO detekcji nukleacji — „załamanie przez generowanie obiektów" NIE wystąpiło; P2a gate 12/12 PASS (próżnia zostaje w g≡1, |ΔE|/E ≤ 6e−16); P2c 18/18 relaksacja do próżni, zero alarmów detektora. Q3: NIE WYKONANE (LOCK: tylko przy Q2-PASS-STATIC; wariant kaskadowy też nie — brak nukleacji). g_floor z derywacji: √{0.197,0.298,0.331} = {0.4438468,0.5458938,0.5753260}."
anti_lakatos_lock: PRESERVED
tags: [metametric-boundary, creation-cost, substrate-floor, nucleation-detector, inconclusive-q2, upward-runaway, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_breakdown_handling.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../op-substrate-fluctuation-channel-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-metametric-boundary

**Status: CLOSED-EXECUTED (2026-09-01, jedna sesja: LOCK → method_decisions
→ Phase 1 → Phase 2; Phase 3 warunkowa — warunek niespełniony).**
Kryteria LOCKa (`Phase0_balance.md` §3, §6) stosowane DOSŁOWNIE; zero zmian
progów/detektora/seeda po starcie (jedna korekta implementacyjna — §4).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q1** (granica w kontinuum) | **Q1-POS** (litera §3 Phase 1) | znaki kosztów kreacji w policzonym zbiorze MIESZANE — istnieje konfiguracja z μ=0 oddzielająca reżimy; uderzające: **próżnia jest już PO stronie opłacalnej kreacji** (μ(sol|vac) < 0) |
| P2a (gate podłogi) | **PASS 12/12** | próżnia bez zaburzeń zostaje w g≡1 (∥g−1∥∞=0); gate energii hamiltonowskiej ≤ 6.0e−16 (próg 1e−6) tam, gdzie podłoga nieaktywna |
| **Q2** (podłoga substratowa) | **Q2-INCONCLUSIVE** (litera §3 Phase 2: „załamanie NIE-nukleacyjne = INCONCLUSIVE, nie pozytyw") | 18/18 biegów tachionowych załamuje się przez ucieczkę **g→+∞** (nie przez generowanie obiektów); podłoga QB-2 nigdy nieaktywowana |
| P2c (kontrola detektora) | **PASS 18/18** | sektor m²=+γ: wszystkie starty relaksują do próżni, ZERO fałszywych alarmów — detektor ważny |
| **Q3** (stabilność stanu zrelaksowanego) | **NIE WYKONANE** | LOCK dopuszcza Phase 3 wyłącznie przy Q2-PASS-STATIC; wariant kaskadowy wymagał Q2-PASS-NUCLEATION — żaden nie zaszedł |

## 1. Derywacja g_floor (FROZEN przed obliczeniami, z cytatem — MD §1)

Pole g akcji kanonicznej ≡ φ z prop:substrate-action dodatku B
(człon kinetyczny φ⁴(∇φ)² = K(g)=g⁴; EL = eq:EL-alpha2 z potencjałem);
eq:K-geometric: „φ_i = (Φ_i/Φ₀)^{1/2}"; Φ₀=Φ_vac z Drogi 1
thm:beta-eq-gamma-triple („Definicja Φ₀ jako wartości próżniowej…
φ_vac=1"). Stąd **g_floor = √(Φ_c/Φ_vac)**:

| Φ_c/Φ_vac (QB-2) [INPUT] | g_floor | g_thr detektora |
|---|---|---|
| 0.197 | 0.4438468 | 0.7219234 |
| 0.298 | 0.5458938 | 0.7729469 |
| 0.331 | 0.5753260 | 0.7876630 |

Implementacja: kara C² V_fl=(κ/3)(g_floor−g)³ dla g<g_floor, κ=100
(FROZEN, bez strojenia); zachowuje strukturę wariacyjną (gate energii
osiągalny i zdany).

## 2. Phase 1 — Q1: równanie stanu kreacji (liczby)

- **P1a (sympy, 4/4 PASS):** U(0)=0; U(1)=1/56; U″(1)=−1;
  solveset(U′=0,(0,1))=∅, U′=g⁶(1−g)>0 na (0,1) — U ściśle rosnące,
  próżnia na maksimum lokalnym, relaksacja w dół bez pośredniego minimum.
- **P1b (kwadratura, R=60):** ΔE_create(soliton μ | próżnia):
  h=0.05: −0.178794; h=0.025: −0.178958; **h=0.0125: −0.178999**
  (zbieżność 4.1e−5). ΔE_create(soliton μ | stan pusty):
  **+16156.583** (h=0.0125; zdominowane objętościowym kosztem próżni
  U(1)·(4π/3)R³ ≈ 16156.8 — rośnie z R³).
- **P1c (sieć 2π z npz, bez relaksowania):** ε(2π) = **−0.00051944**
  (A1.0, N=48; N=32: −0.00055526; A0.7 identyczne — duplikat tła;
  |Δ32→48|=3.6e−5). Sieć leży PONIŻEJ próżni energetycznie.
- **Werdykt wg litery:** znaki {−, +, −} w policzonym zbiorze ⟹
  istnieje konfiguracja z μ=0 oddzielająca reżim μ<0 od μ>0 ⟹
  **Q1-POS**. Charakteryzacja: separacja μ=0 przebiega między stanem
  pustym a konfiguracjami metrycznymi — **próżnia g=1 i sieć 2π są już
  w reżimie opłacalnej kreacji** (μ<0). LOCK spodziewał się negatywu;
  wynik pozytywny jest mocniejszy: nie trzeba szukać granicy „za
  próżnią" — próżnia sektora tachionowego już ją przekroczyła.

## 3. Phase 2 — Q2: relaksacja z podłogą (rachunek centralny)

**P2a (gate):** 12/12 PASS — (i) gradient flow z podłogą na dokładnej
próżni: ∥g−1∥∞ = 0.0 przez t=10 (radialnie h=0.0125 i 3D N=32, 3 podłogi);
(ii) hamiltonowska ewolucja K_ε=0.2 z U_tot tam, gdzie podłoga
nieaktywna: |ΔE|/E ≤ 1.77e−16 (radialnie) / 6.02e−16 (3D) — próg 1e−6.

**P2b (18 biegów tachionowych): wyniki per start × podłoga × siatka**
(status; t_break; ostatnia zdrowa próbka [min g, max g]):

| Start | Siatka | f1 (0.4438) | f2 (0.5459) | f3 (0.5753) |
|---|---|---|---|---|
| soliton μ (radial) | h=0.025 | BRKDWN t=2.75 [0.777,1.48] | BRKDWN t=2.75 | BRKDWN t=2.75 |
| soliton μ (radial) | h=0.0125 | BRKDWN t=3.13 [0.778,1.93] | BRKDWN t=3.13 | BRKDWN t=3.13 |
| sieć 2π (npz) | N=32 | BRKDWN t=1.25 [0.607,1.63] | BRKDWN t=1.25 | BRKDWN t=1.25 |
| sieć 2π (npz) | N=48 | BRKDWN t=1.49 [0.610,1.52] | BRKDWN t=1.49 | BRKDWN t=1.49 |
| próżnia+szum | N=32 | BRKDWN t=11.06 [1.908,1.912] | BRKDWN t=11.06 | BRKDWN t=11.06 |
| próżnia+szum | N=48 | BRKDWN t=11.06 [1.854,1.857] | BRKDWN t=11.06 | BRKDWN t=11.06 |

Wszystkie załamania = **skończono-czasowa ucieczka g→+∞** (U(g)→−∞ dla
g→∞; energia spada monotonicznie — np. szum N=32: E 4.43 → −2205 przed
załamaniem). Fakty kluczowe:
1. **Podłoga nigdy nieaktywowana:** min g w żadnym biegu nie osiągnęło
   g_floor (soliton: 0.777 bez spadku; sieć: 0.607 bez spadku; szum:
   min g ROSŁO do 1.85 — całe pole poszło w górę). Wynik jest
   IDENTYCZNY dla wszystkich trzech podłóg — zerowa czułość na g_floor.
2. **Zero nukleacji:** detektor (ważny — P2c) nie zarejestrował ani
   jednego wzrostu liczby obiektów; pre-rejestrowany pozytyw autora
   („załamanie przez GENEROWANIE OBIEKTÓW") NIE wystąpił — załamanie
   idzie przez WZROST jednego istniejącego maksimum, nie mnożenie.
3. **Spójność między siatkami:** kierunek i charakter załamania zgodne
   na obu siatkach każdej pary; dla szumu t_break identyczny co do
   0.01 (wspólne pasmowe pole startowe); trajektoria zweryfikowana
   niezależnie jawnym Eulerem dt=2e−4 (ten sam runaway, szybszy —
   semi-implicit nie generuje artefaktu).
4. **Start genezowy (szum, seed=20260901):** mody k≥1 w pudle L=2π są
   tłumione (tempo 1−k²K(g) ≤ 0 dla k≥1) — szum homogenizuje się w
   ~1 j.cz., a jedyny rosnący stopień swobody (k=0, znak = znak
   średniej seeda: +2e−4) prowadzi całe pole w górę. W zalockowanym
   pudle 2π start genezowy NIE MÓGŁ wytworzyć struktury przestrzennej
   (ograniczenie klasy zbadanej, raportowane wprost).

**Klasyfikacja wg litery LOCKa:** nie Q2-PASS-NUCLEATION (brak
nukleacji), nie Q2-PASS-STATIC (żaden stan stacjonarny niestały), nie
Q2-FAIL (biegi NIE relaksują do stanu jednorodnego g≡1/g≡g_floor —
załamują się). Pozostaje: **„Załamanie NIE-nukleacyjne (NaN…) =
INCONCLUSIVE, nie pozytyw"** — Q2-INCONCLUSIVE, bez reinterpretacji.

**P2c (kontrola nieusuwalna): PASS 18/18** — sektor m²=+γ
(U_stab=−g⁷/7+g⁸/8): soliton → próżnia (t=16), sieć 2π → próżnia (t=19),
szum → próżnia (t=6); wszystkie STATIONARY ∥ġ∥∞<1e−8, g≡1; ZERO alarmów
detektora (w tym dla startów z zasianymi obiektami N_seed=1) — detektor
zdolny do FAIL i czysty.

## 4. Anti-Lakatos / korekty / higiena

- ✓ LOCK + method_decisions zamknięte PRZED kodem; kryteria/progi/
  detektor/seed nietknięte po starcie; zakaz strojenia g_floor i κ
  dotrzymany (wyniki identyczne dla 3 podłóg — strojenie nie miałoby
  zresztą czego uratować).
- ✓ **Korekta 1** (`Phase_correction_note_breakdown_handling.md`,
  zapisana PRZED użyciem wyników): silnik nie przechwytywał
  niefinityczności (ValueError/OverflowError solvera) — załamanie
  crashowało proces zamiast klasyfikacji BREAKDOWN. Czysta obsługa
  błędu; pierwotny output zachowany
  (`Phase2_results/job_tach_sol_f2_h025_pre_correction.txt`).
- ✓ Kosmetyka: jawne `axes` w irfftn (DeprecationWarning numpy 2.x),
  zero wpływu na wartości; P2a wykonane przed tą poprawką pozostaje
  ważne (wartości deterministyczne).
- ✓ Tła npz READ-ONLY: mtime 2026-08-31 21:41:07 niezmieniony po
  wszystkich odczytach (weryfikowany w P1c i w każdym starcie `lat`).
- ✓ Wynik negatywny/nierozstrzygnięty wprost; INCONCLUSIVE NIE
  reinterpretowane jako „prawie nukleacja" (forbidden move 5).
- ✓ Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowane; git nieużywany;
  katalogi innych cykli tylko odczyt; pełne ścieżki bez `cd`;
  `ls` po zapisach (zero artefaktów zagnieżdżonych).
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0.

## 5. Odczyt (deskryptywnie, bez claimów poza klasą zbadaną)

1. **Q1-POS przesuwa obraz:** koszt kreacji solitonu μ z próżni jest
   UJEMNY (−0.179) — „granica metametryczna" nie leży przed próżnią;
   sektor tachionowy w próżni już jest w reżimie opłacalnej kreacji.
   To spójne z Q2: relaksacja natychmiast konsumuje ten zysk — ale
   kanałem WZROSTU istniejącego obiektu (g→+∞, Φ→∞), nie mnożenia
   obiektów.
2. **Podłoga QB-2 adresuje niewłaściwy brzeg:** spinodala rozrzedzenia
   ogranicza pole od DOŁU (g_floor<1), a dominujący kanał ucieczki
   zalockowanego modelu jest od GÓRY (U nieograniczone z dołu dla
   g→∞ — własność znana, MD op-3d-canonical §1). W klasie zbadanej
   hipoteza „układ relaksuje do granicy metametrycznej i tam się
   zatrzymuje" nie dostała nośnika.
3. Rozstrzygnięcie wymaga któregoś z: (a) mechanizmu górnego
   ograniczenia z korpusu (nie ad-hoc), (b) pudła nadkrytycznego
   L>2π dla struktury przestrzennej, (c) programu genezy Γ+s_i
   poziomu 0 — wszystkie user-gated w NEEDS.md.

## 6. Deliverables

`Phase0_balance.md` (LOCK) · `Phase_method_decisions.md` (FROZEN) ·
`Phase1_creation_cost.py` → `Phase1_output.txt` ·
`Phase2_floor_relax.py` → `Phase2_output.txt` +
`Phase2_relaxed_states.npz` (36 stanów) + `Phase2_results/` (logi
per-bieg, json serii, npz stanów, output pierwotny pre-correction) ·
`Phase_correction_note_breakdown_handling.md` ·
`Phase3_metametric_spectrum.py` — NIE UTWORZONY (warunek LOCKa
niespełniony) · `NEEDS.md` · `README.md` (log).
