---
title: "Phase_FINAL_close — zamknięcie: GATE-FAIL-STOP na P1c (most radialny→kartezjański): dyskretyzacja 3D przy zalockowanym h≈0.4 NIE reprodukuje kotwicy modelu M0-f_ε (λ_min(3D)=−8.81 vs −1.3896; t*_izo(3D)=0.17 vs 3.62±15%) — przyczyna strukturalna: wąska sferyczna kieszeń Q₆₃ (−15.5 przy g≈0.75, strefa ściany f_ε) wymaga h≲0.025 (≈1500³ w 3D, niewykonalne); P1a/P1b/P1d czyste; Phase 2–4 NIEURUCHOMIONE wg LOCKa; pytanie Q (ω²_min(3D)) pozostaje OTWARTE"
date: 2026-08-31
type: phase-final-close
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
status: CLOSED
verdict: "GATE Phase 1: FAIL na P1c ⟹ STOP-MASZYNERIA (LOCK §3 Phase 1, dosłownie; drzewo §6 gałąź 'P1b/P1c FAIL → STOP'). Składowe: P1a PASS (dyspersja próżni 3D exact: maxerr najniższej gałęzi 6.02e−4 tach / 3.44e−4 stab przy N=32, ratio N32/N64 = 3.999/3.999 — rząd 2 czysty); P1b PASS (kotwica radialna #63 verbatim: λ_min(w1)=−1.38962 vs −1.3896±1e−3; t*=3.62 we wszystkich biegach a=±0.01, dt=0.004/0.002; gate energii ≤1.64e−8); P1c FAIL (λ_min(3D, h=0.3947)=−8.80982, λ_min(3D, h=0.30)=−7.53713 vs kotwica −1.3896 — odchył 534%/442% przy gate ±5%; 'poprawa z zagęszczeniem' formalnie zachodzi, ale bez zbieżności do kotwicy; t*_izo(3D)=0.17–0.18 vs 3.62±15% FAIL; gate energii ewolucji 3D PASS ≤5.2e−7); P1d PASS (|ΔE|/E ≤ 7.4e−15 dla K_ε=0.2/0.1 i F_ε=0.2, dt-konsystencja ≤1.3e−13). Diagnoza (addendum diagnostyczne, zero zmian kryteriów): operator fluktuacji modelu M0-f_ε ma wąską SFERYCZNĄ kieszeń potencjału Q₆₃: min pointwise −15.51 przy r=3.382 (g=0.7509 — strefa przejścia f_log przez 0 przy g=e^{−1/4}≈0.779, gdzie F′,F″ regularyzacji ε spike'ują), F_ε(g_min=0.4107)=3.9e−3; kod RADIALNY przy TYM SAMYM h=0.3947 daje −1.652 (nie −8.8; gruba siatka radialna aliasuje powłokę, gęsta kątowo siatka 3D zawsze ją trafia setkami węzłów → klaster modów pasożytniczych −8.81). Hipoteza błędu konwencji wagi (3D=F-ważone zamiast weight-1; podejrzenie z bliskości −7.54 do F-ważonej kotwicy −7.8579) PRZETESTOWANA I OBALONA (addendum 2): radialne F-ważone przy h=0.3947/0.30 daje −114.5/−117.6 ≠ −8.8/−7.5; iloraz Rayleigha DOKŁADNEGO modu w1 interpolowanego do operatora 3D = −1.518 (N=76; 9.3% od kotwicy) i −1.419 (N=100; 2.1%) — operator 3D weight-1 reprezentuje mod kotwiczny poprawnie, a λ_min jest przejęte przez mody pasożytnicze powłoki (mod −8.81: 99% masy w |r−3.38|<0.3, ⟨r⟩=3.380, 0% w rdzeniu; min diag(A) = −8.43/−6.90 śledzi λ_min). Tłumienie artefaktu wymaga 6F_powłoki/h² > głębokość kieszeni ⟹ h ≲ 0.15 (N≳200³ przy L=30), a zbieżność weight-1 do ±5% dodatkowo h ≲ 0.1–0.2 (radialnie: 4.5% przy h=0.2, 3.1% przy h=0.1) — czyli bramka w formie zalockowanej (h≈0.4) jest NIEOSIĄGALNA dla poprawnej implementacji. Maszyneria 3D per se zwalidowana (P1a exact, P1d 1e−15, RQ modu kotwicznego 2.1% przy N=100) — FAIL jest barierą ROZDZIELCZOŚCI zalockowanego mostu w modelu f_ε, nie błędem dyskretyzacji kartezjańskiej ani konwencji. WERDYKT Q: NIEROZSTRZYGNIĘTE (STOP przed Phase 2–4; kategorie Q-PASS/Q-FAIL/Q-INCONCLUSIVE Phase 4 nie były ewaluowane) — pytanie ω²(n) w 3D pozostaje OTWARTE, ostatnia droga rachunkowa hipotezy stabilizacji gęstością NIE została ani potwierdzona, ani obalona."
anti_lakatos_lock: PRESERVED
tags: [3d-lattice, tachyonic-sector, machinery-gate-fail, bridge-gate, soft-wall-resolution, stop-report, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_eigsh.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-3d-lattice-bath-stability

Obliczenia i zamknięcie 2026-08-31 (jedna sesja implementatora). Kryteria
LOCKa (Phase0_balance.md §3, §6) stosowane DOSŁOWNIE; zero zmian
progów/punktów/siatek po starcie. STOP wykonany dokładnie tam, gdzie LOCK
go pre-rejestrował: P1c FAIL ⟹ „STOP (dyskretyzacja 3D nieadekwatna —
raport bez Phase 2–4)".

**Rejestr WEJŚĆ:** g₀=2.02117 (kalibracja μ #63) [użyte w P1b/P1c],
ε=0.2, kontrola 0.1 [P1b/P1c/P1d], d*₁=3.0790 [NIEUŻYTE — Phase 2 nie
wystartowała], β=γ=1 [P1a/P1d].

---

## 0. Werdykt

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| Gate Phase 1 | **FAIL (P1c)** | most radialny→kartezjański nie reprodukuje kotwicy −1.3896/t*=3.62 przy zalockowanym h≈0.4 (λ_min(3D)=−8.81; t*=0.17) |
| Phase 2–4 | **NIEURUCHOMIONE** | dosłowna litera LOCKa §3/§6: P1c FAIL ⟹ STOP bez Phase 2–4 |
| **Q** (∃d: ω²_min(d)>0 w 3D?) | **NIEROZSTRZYGNIĘTE — OTWARTE** | żadna liczba ω²_min(3D) nie powstała; negatyw 1D (−1.22) nadal jawnie nieprzenośny (LOCK §0) |

**Sens metodologiczny (deskryptywnie):** to nie jest negatyw hipotezy —
to negatyw WYKONALNOŚCI zalockowanego mostu P1c w dosłownym modelu
M0-f_ε przy siatkach osiągalnych w 3D. Wąska sferyczna struktura ściany
f_ε (skala ~0.02–0.1) jest tania w 1D radialnym (N=4000), a zaporowa
w 3D kartezjańskim (~1500³). Maszyneria 3D sama w sobie przeszła
wszystkie pozostałe bramki (P1a exact rząd 2; P1d gate energii 1e−15;
P1b kotwica radialna verbatim).

## 1. Phase 1 — wyniki bramek (Phase1_output.txt)

- **P1a (dyspersja próżni 3D exact): PASS.** Pudło L=2π, operator
  kanoniczny (K=g⁴, spektra bez regularyzacji): najniższa gałąź
  w Γ/X/M/R: maxerr = 6.022e−4 (tach, N=32) i 3.441e−4 (stab, N=32),
  gate ≤1e−3 ✓; ratio N32/N64 = 3.999/3.999 ∈ [3,5] ✓ (rząd 2).
  Γ: błąd 6e−14 (mod stały reprezentowany dokładnie). Deskryptywna
  tabela 10 gałęzi nosi artefakt ARPACK (gubione krotności przy
  tol=1e−6) — wyjaśnione i skorygowane w
  [[Phase_correction_note_eigsh.md]] (bramka NIEDOTKNIĘTA: jej metryka
  FROZEN to gałąź najniższa, niezdegenerowana).
- **P1b (kotwica radialna, model #63 verbatim): PASS.**
  λ_min(w1) = −1.38962 (kotwica −1.3896 ± 1e−3 ✓); t* = 3.62 we
  wszystkich trzech biegach (a=±0.01; dt=0.004 i 0.002; detektor
  g≤0.01 co krok — korekta 1 poprzednika); gate energii ≤ 1.64e−8 ✓.
- **P1c (most radialny→kartezjański, model #63): FAIL.**
  Pudło L=30 periodyczne, soliton w centrum, interpolacja splajnowa
  profilu radialnego (g, g′, g″ analitycznie):

  | siatka | h | λ_min(3D) | vs kotwica −1.3896 |
  |---|---|---|---|
  | N=76³ | 0.3947 | **−8.80982** | odchył 534% (gate ±5%) — FAIL |
  | N=100³ | 0.3000 | **−7.53713** | odchył 442%; „poprawa" formalna, bez zbieżności |

  10 najniższych przy N=76: klaster −8.810…−8.804 (mody pasożytnicze
  zlokalizowane na powłoce ściany). t*_izo(3D) = 0.17–0.18 (a=±0.01,
  dt=0.02/0.01) vs 3.62±15% = [3.08, 4.16] — FAIL (natychmiastowy
  breakdown w nierozdzielczonej ścianie; gate energii samej ewolucji
  PASS ≤ 5.2e−7 — integrator czysty, pole niefizyczne).
- **P1d (gate energii ewolucji 3D): PASS.** L=2π, N=32, a=1e−3,
  t_end=4: |ΔE|/E ≤ 7.42e−15 dla wszystkich wariantów (kanoniczna
  K_ε=0.2, K_ε=0.1, #63-owa F_ε=0.2; dt=0.004 i 0.002);
  dt-konsystencja ‖Δg‖∞ ≤ 1.29e−13 (gate 1e−6) ✓.

## 2. Diagnoza FAIL P1c (addendum diagnostyczne; zero zmian kryteriów)

`Phase1_addendum_diag_p1c.py` + `Phase1_output_addendum_diag.txt`
(dokumentacja przyczyny, policzona PO orzeczeniu FAIL, nie wpływa na
werdykt):

1. **Profil μ:** g ∈ [0.4107, 2.0212]; F_ε(g_min) = 3.90e−3 (ściana
   miękka: kinetyka niemal znika), F_ε(g_max) = 3.82.
2. **Potencjał Q₆₃ pointwise:** min = **−15.51 przy r=3.382**
   (g=0.7509) — wąska kieszeń w strefie przejścia f_log przez zero
   (g = e^{−1/4} ≈ 0.779), gdzie pochodne regularyzacji F′, F″
   spike'ują; w rdzeniu (g≈2.02) Q ≈ −6.0 szeroko.
3. **Kod radialny przy równym h** (R=30, ten sam profil):
   h=0.3947 → λ_min=−1.652; h=0.30 → −1.650; h=0.20 → −1.327;
   h=0.10 → −1.347; h=0.05 → −1.427; h=0.025 → −1.3904;
   h=0.015 → −1.3896. Zbieżność do kotwicy dopiero przy h ≲ 0.025;
   przebieg niemonotoniczny (aliasing wąskiej powłoki).
4. **Dlaczego 3D ≠ radial przy równym h:** gruba siatka radialna
   trafia powłokę r≈3.38 jednym–dwoma węzłami (albo wcale — aliasing);
   siatka kartezjańska 3D próbkuje CAŁĄ sferę setkami węzłów —
   zawsze zawiera węzły z Q ≈ −9…−15 przy kinetyce sprzęgła zbyt
   słabej, by je ukarać ⟹ klaster modów pasożytniczych ≈ −8.8.

**Addendum 2 — rozstrzygnięcie hipotezy błędu konwencji wagi**
(`Phase1_addendum_diag2_weight.py` + `Phase1_output_addendum_diag2.txt`;
test wykonany na żądanie sesji koordynującej PRZED utrzymaniem STOP;
wynik: BRAK błędu implementacji ⟹ korekta per LOCK §4 p.1 NIELEGALNA):

- **T1 (test rozstrzygający koordynatora):** radialne widmo F-WAŻONE
  przy h=0.3947/0.30 (R=30): λ_min(F) = −114.47/−117.56 — NIE pokrywa
  się z 3D (−8.81/−7.54); przy h=0.015: −7.8579 = kotwica F-ważona #63
  odtworzona. Bliskość −7.54 do −7.8579 była koincydencją. Operator 3D
  NIE rozwiązuje problemu F-ważonego — konwencja weight-1 potwierdzona.
- **T2 (dowód poprawności operatora):** iloraz Rayleigha dokładnego
  modu w1 (u₁/r z N=4000, interpolowany) w operatorze 3D:
  RQ = −1.51846 (N=76; 9.3% od −1.38962) i −1.41924 (N=100; 2.1%) —
  operator 3D reprezentuje mod kotwiczny z błędem malejącym z h;
  problemem nie jest reprezentacja modu, lecz ISTNIENIE głębszych
  modów pasożytniczych na niedorozdzielczonej powłoce.
- **T3 (lokalizacja):** mod λ=−8.8098 ma ⟨r⟩=3.380, 99.0% masy
  w powłoce |r−3.38|<0.3, 0.0% w rdzeniu r<2 — czysty artefakt kieszeni.
- **T4 (skala):** min diag(A) = −8.433 (N=76) / −6.902 (N=100) przy
  r=3.81 — pointwise minima diagonali śledzą λ_min(3D): mody
  pasożytnicze żyją na pojedynczych węzłach powłoki.

5. **Skala wykonalności (ilościowo):** artefakt kieszeni jest tłumiony
   dopiero, gdy diagonalna kara kinetyczna 6F_powłoki/h² (F≈0.05)
   przewyższy głębokość kieszeni (~14) ⟹ h ≲ 0.15 ⟹ N ≳ 200 na wymiar
   przy L=30 (≈8×10⁶ węzłów — ciężkie, ale nie absurdalne); niezależnie
   błąd dyskretyzacji weight-1 ±5% wymaga h ≲ 0.1–0.2 (tabela radialna:
   4.5% przy h=0.2; 3.1% przy h=0.1; pełna zbieżność h≲0.025). Bramka
   ZALOCKOWANA żąda ±5% PRZY h≈0.4 — nieosiągalne dla poprawnej
   implementacji; kryterium NIE jest zmieniane (anti-Lakatos), STOP
   utrzymany. Bramka P1c zadziałała dokładnie tak, jak ją LOCK
   pre-rejestrował: wykryła nieadekwatność dyskretyzacji 3D (w sensie:
   zalockowanej rozdzielczości) dla modelu mostu f_ε.

## 3. Obowiązkowy punkt LOCKa: trend ω²_min(3D) vs 1D (−1.22)

**NIEDOSTĘPNY.** STOP na bramce Phase 1 ⟹ żadne tło 3D nie zostało
zrelaksowane i żadne ω²_min(3D) nie istnieje. Trend względem wyniku 1D
(bloch-chain: −1.222…−1.230) NIE może być raportowany; pytanie, czy 3D
podnosi ω²_min względem −1.22, pozostaje w całości otwarte. Argument
przenośności negatywu 1D (LOCK §0: twierdzenie węzłowe Hilla nie
obowiązuje w 3D; ogony 1/r) pozostaje w mocy — kierunek NIEZNANY.

## 4. Korekty, dodatki, incydenty (pełna lista)

1. **Rozdział modeli P1 vs P2–4** zamrożony PRZED startem
   (Phase_method_decisions.md §0): bramki P1b/P1c w dosłownym modelu
   #63 M0-f_ε (bo tam żyją kotwice −1.3896/3.62), rachunek centralny
   w akcji kanonicznej (przyjęcie LOCKa za bloch-chain). FAIL P1c
   dotyczy mostu w modelu f_ε.
2. **Korekta 1** ([[Phase_correction_note_eigsh.md]], zapisana przed
   Phase 3): eigsh tol=1e−6 gubi krotności zdegenerowane (dowód:
   diagnostyka Γ + tabela P1a); korekta tol=0 + kontrola pokrycia
   translacji. Phase 3 nigdy nie wystartowała — korekta pozostaje
   udokumentowana, użyta wyłącznie w niewykonanym skrypcie.
   Pierwotny `Phase1_output.txt` zachowany bez zmian.
3. **Poprawka freeze przed pierwszym biegiem** (bez correction_note —
   zero obliczeń istniało): v0 Lanczosa zmieniony ze stałego na
   deterministyczny pseudolosowy (stały leżałby w podprzestrzeni
   symetrycznej i mógłby ominąć mody antysymetryczne).
4. **Skrypty Phase 2/Phase 3 istnieją, NIEURUCHOMIONE** (pisane
   równolegle z biegiem P1, przed poznaniem wyniku bramki):
   `Phase2_relax_lattice3d.py`, `Phase3_bloch3d.py` — zero wyników,
   zero plików wyjściowych (brak `Phase2_backgrounds3d.npz`,
   `Phase3_results3d.json`). Zachowane jako artefakty dla ewentualnego
   przyszłego cyklu (user-gate w NEEDS).
5. **Addenda diagnostyczne** (dozwolona dokumentacja przyczyny FAIL,
   po orzeczeniu FAIL; zero zmian kryteriów):
   `Phase1_addendum_diag_p1c.py` + `Phase1_output_addendum_diag.txt`
   (kieszeń Q, radial przy równych h) oraz
   `Phase1_addendum_diag2_weight.py` + `Phase1_output_addendum_diag2.txt`
   (T1–T4: obalenie hipotezy błędu konwencji wagi, RQ modu kotwicznego,
   lokalizacja modu pasożytniczego). Ścieżka korekty per LOCK §4 p.1
   rozważona i ZAMKNIĘTA: brak udokumentowanego błędu implementacji.
6. Incydent środowiskowy bez wpływu na wyniki: sandbox odrzuca
   niektóre heredoc-i i ścieżki poza vaultem — wszystkie rachunki
   przez pliki skryptów w katalogu cyklu, pełne ścieżki, bez `cd`.

## 5. Pliki cyklu

Obliczeniowe: `Phase1_gate3d.py` + `Phase1_output.txt`;
`Phase1_addendum_diag_p1c.py` + `Phase1_output_addendum_diag.txt`;
`Phase1_addendum_diag2_weight.py` + `Phase1_output_addendum_diag2.txt`.
Nieuruchomione (artefakty): `Phase2_relax_lattice3d.py`,
`Phase3_bloch3d.py`. Metodologiczne: `Phase_method_decisions.md`,
`Phase_correction_note_eigsh.md`. Zamykające: ten plik, `NEEDS.md`
(user-gated), `README.md` (zaktualizowany). Rdzeń `.tex`, STATE.md,
git, katalogi innych cykli — NIETKNIĘTE.

## 6. Mapowanie na drzewo decyzyjne LOCKa §6

Gałąź: **„P1b/P1c FAIL → STOP (maszyneria; raport, bez Phase 2–4)"** —
wykonana dosłownie. Do NEEDS trafia wniosek METODOLOGICZNY (bariera
rozdzielczości mostu f_ε w 3D + opcje user-gated: most w modelu
kanonicznym / współrzędne adaptowane / AMR), NIE ontologiczny.
Hipoteza stabilizacji gęstością w 3D pozostaje OTWARTA i nadal jest
jedyną niewykluczoną geometrią tej linii (negatyw 1D nieprzenośny —
pre-rejestracja LOCKa §0 zachowuje ważność).
