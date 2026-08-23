---
title: "Phase_FINAL_close — zamkniecie: geneza samopodtrzymujacych struktur z golego substratu"
date: 2026-07-04
type: phase-final-close
tgp_owner: research/op-bare-substrate-genesis-2026-07-04
status: CLOSED
verdict: PASS (G1-G6)
anti_lakatos_lock: PRESERVED
tags: [substrate-genesis, soliton, lock, self-sustaining, nucleation, gravity-bridge, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../op-lock-interaction-gravity-2026-07-04/Phase0_balance.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS, G3 PASS, G4 PASS, G5 PASS, G6 PASS (wersja mocna).**

Wg zalockowanej reguly decyzyjnej (Phase0_balance.md, sekcja 6): substrat
POTRAFI wygenerowac samopodtrzymujace, lockowane struktury `Phi > 0` przez
kolektywny lock. Fundament mostu do grawitacji STOI. Nastepny cykl
autoryzowany do zaprojektowania.

## 1. Co testowano

Hipoteza autora (po zamknieciu `op-blocked-soliton-bang`): gola faza
substratu (`Phi = 0`) jest metastabilna; pojedyncze zaburzenie relaksuje
z powrotem, ale dostatecznie gesta konfiguracja nakladajacych sie zaburzen
podtrzymuje sie NAWZAJEM (lock), z jednej dynamiki pola, bez osobnej
reguly przezycia. Model: siatka 2D, `s(x)`, `Phi = s^2`, przeplyw
gradientowy `ds/dtau = kappa*Lap(s) - V'(s)`, potencjal z metastabilnym
`s = 0` i glebsza prawdziwa proznia `s* ~ 1.174` (a=0.5, b=1.6, c=1.0).

## 2. Przebieg

- **Phase 0**: LOCK modelu + rachunki analityczne; 6 poprawek pre-code
  (sekcja 11 Phase0_balance.md), wszystkie PRZED pierwsza linijka kodu.
  Weryfikacja numeryczna rachunkow: `Phase0_output.txt` (PASS).
- **Phase 1**: scenariusze `bare`, `single(A0)` skan, `cluster(A0, D)`.
  Wynik: G1-G4 PASS. `Phase1_output.txt`.
- **Phase 2**: skan progu (bisekcja po A0 klastra przy D in {2,3,4}),
  rozszerzenie skanu singla, kontrola pinningu N=256, `DeltaE_create`
  FAR/CORE/FRONT. Wynik: G5-G6 PASS. `Phase2_output.txt`.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Progi

```text
bare:                relaksuje do Phi = 0 w calosci (G2)
single, w=1.5:       A_c(single) w (4.0, 5.0)   -> single potrzebuje
                     amplitudy ~4x powyzej prawdziwej prozni s*=1.174
cluster K=7, D=2:    A0_c ~ 0.4578  (przedzial 0.45312-0.46250)
cluster K=7, D=3:    A0_c ~ 0.4672  (przedzial 0.46250-0.47188)
cluster K=7, D=4:    A0_c ~ 0.5703  (przedzial 0.56562-0.57500)
przewaga kolektywna: >= 7.0x - 8.7x (dolne oszacowanie)
```

Klaster lockuje przy amplitudzie ledwie ~7-34% powyzej amplitudy bariery
(|s_bar| = 0.426) — pojedynczy seed przy tej amplitudzie zanika
bezwarunkowo. "Single decays / many lock" potwierdzone wprost.

### 3.2 Charakter progu (G5)

Ostra bifurkacja: szerokosc przedzialu progowego 0.0094 po A0, zero stanow
posrednich (po jednej stronie pelny zanik do Phi=0, po drugiej lock
dokladnie w Phi* = 1.3787). Pelna relaksacja przez caly run (H_Gamma
nierosnace, zero clampow). Kontrola pinningu: przypadki graniczne D=3
powtorzone na N=256/dx=0.25/dt=0.005 — klasyfikacja identyczna. Lock NIE
jest artefaktem siatki.

### 3.3 Koszt kreacji (G6, kluczowe wg autora)

Stan: cluster(1.4, D=3) po steps_probe=3000. Bump A=0.10, w=1.0:

```text
FAR   (Phi=0, ~29 od struktury):  DeltaE = +0.014390
CORE  (wnetrze locku, s ~ s*):    DeltaE = +0.022790
FRONT (granica metryki):          DeltaE = -0.014581   < 0
```

Kreacja kolejnego zaburzenia NA GRANICY regionu metrycznego jest
energetycznie samoczynna; w martwym substracie kosztuje; we wnetrzu
kosztuje najwiecej. Predykcja pre-code (CORE > FAR, bo V''(s*)=0.879 >
V''(0)=0.5) potwierdzona. To daje energetyczna podstawe intuicji
"kreacja trwa na granicy metryki".

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- minimalny substrat (jedno pole, jedna dynamika, zero regul pomocniczych)
  generuje trwale struktury metryczne wylacznie przez kolektywnosc,
- prog jest prawdziwa bifurkacja typu nukleacyjnego (krytyczny rozmiar
  obszaru nadbarierowego, nie krytyczna amplituda pojedynczego seedu),
- granica metryki jest miejscem ujemnego kosztu kreacji.

CZEGO WYNIK NIE ZNACZY:
- lockowana struktura to ROSNACA domena prawdziwej prozni (front
  ekspanduje), nie skonczony soliton — stabilizacja rozmiaru pozostaje
  otwarta (por. op-wall-dynamics),
- `tau` to parametr relaksacji, nie czas fizyczny,
- zaden element nie jest dowodem GR/MOND/kosmologii; mechanizm jest
  klasyczna teoria nukleacji (Allen-Cahn z przechylonym podwojnym dolkiem)
  zrealizowana w ontologii substratu — cykl potwierdzil, ze TA klasa
  modeli robi to, czego wymagala hipoteza, nie wiecej.

## 5. Errata / transparencja

Dwa bledy RAPORTOWANIA (nie dynamiki, nie progow) wykryte i naprawione
w dniu runow, przed interpretacja:

1. Phase 1: linia podsumowania oceniala G3 ostrzej niz brzmienie LOCKa
   (wymagala progu z obu stron skanu). Poprawiono do brzmienia LOCKa;
   brak nadkrytycznego singla w zalockowanym skanie raportowany jako
   obserwacja poza kryterium.
2. Phase 2: podsumowanie przewagi kolektywnej zakladalo "A_c(single)>5.0"
   mimo ze sekcja [A] znalazla lock singla przy A0=5.0. Poprawiono na
   dolne oszacowanie z bracketa (4.0, 5.0).

Zadne progi G1-G6, wspolczynniki ani procedury nie byly zmieniane po
pierwszym runie. Jedyne dozwolone wyjatki techniczne (obnizenie dt,
kontrolny run na drobnej siatce) — dt NIE musialo byc obnizane; kontrolny
run wykonano zgodnie z LOCKiem.

## 6. Decyzja i nastepny cykl

Regula decyzyjna: G2+G3+G4+G5+G6 PASS -> przejdz do mostu do grawitacji.

Nastepny cykl: [[../op-lock-interaction-gravity-2026-07-04/Phase0_balance.md]]
— hipoteza autora: grawitacja jako oddzialywanie pol miedzy zlockowanymi
obiektami, policzalne dla prostych przypadkow, z proba wyprowadzenia
algebraicznych rownan oddzialywania. Status: DRAFT-FOR-LOCK (wymaga
autoryzacji przed kodem).

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK + poprawki pre-code (sekcja 11)
- `Phase0_check.py` / `Phase0_output.txt` — weryfikacja rachunkow
- `Phase1_substrate_genesis.py` / `Phase1_output.txt` — G1-G4
- `Phase2_threshold_and_create_cost.py` / `Phase2_output.txt` — G5-G6
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
