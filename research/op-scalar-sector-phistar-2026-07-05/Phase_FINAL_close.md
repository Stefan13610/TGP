---
title: "Phase_FINAL_close — zamkniecie: kanal dylatacyjny na tle Phi* (galaz C)"
date: 2026-07-05
type: phase-final-close
tgp_owner: research/op-scalar-sector-phistar-2026-07-05
status: CLOSED
verdict: "G1-G2 PASS; G3 FAIL -> STOP przed interpretacja G4/G5 (wg LOCKa); kanal dylatacyjny NIEROZSTRZYGNIETY (nie sfalsyfikowany) — procedura statyczna S/O lamie wlasny warunek symetrii"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-C, dilaton, amplitude-channel, quadrupole, measurement-failure, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS, G2 PASS, G3 FAIL -> STOP przed interpretacja G4/G5/G6
(dokladnie wg zalockowanej reguly: "G2/G3 FAIL -> stop przed
interpretacja (fundament pomiaru)").**

Kanal dylatacyjny NIE zostal rozstrzygniety — ani potwierdzony, ani
sfalsyfikowany. Cykl zdiagnozowal precyzyjnie, DLACZEGO pomiar
w zalockowanej procedurze jest niewazny (sekcja 3.3) — to jest
produkt tego opa, wraz z potwierdzonym fundamentem (ogon algebraiczny
wiru) i zalockowana predykcja B_pred na przyszlosc.

## 1. Co testowano

Po B1 (dlugi zasieg, zly znak) i B2 (geometria niemierzalna na tle
zamrozonym): czy kanal amplitudowy (dylatacyjny) na tle Phi* jest
UNIWERSALNIE przyciagajacy i policzalny. Konstrukcja: kwadrupole S/O
o identycznych pozycjach 4 wirow rozniace sie tylko znakami — sektor
fazowy kasuje sie w E_sym = (E_S+E_O)/2 dokladnie (algebra w LOCKu),
zostaje kanal slepy na znak. Model identyczny z B1 (zero strojenia).

## 2. Przebieg

- **Phase 0**: LOCK przy utworzeniu (autoryzacja: wybor drogi
  "sektor skalarny" przez autora po zamknieciu B1/B2) + jedno
  doprecyzowanie pre-code (floor w P2) przed pierwsza linijka kodu.
  `Phase0_output.txt`: ogon ODE vs 1 - xi^2/r^2 (zgodnosc 0.07%
  przy r=8); E_sym_ansatz(d) rosnace (przewidywane przyciaganie);
  **B_pred = 13.546 (glowna) / 13.486 (kontrola)** — spojnosc miedzy
  siatkami; sygnal przewidywany 0.76. UWAGA pre-code, ktora okazala
  sie prorocza: E_anty/teoria-Greena = 1.34-1.61 przy d in [5,12]
  (rdzenie modyfikuja kanal fazowy przy tych separacjach).
- **Phase 1**: ogon amplitudowy na siatce 2D (G2). PASS (0.26%).
- **Phase 2**: skan S/O + fit laczny + kontrola dx=0.25 + obserwacja
  P5. G3 FAIL -> stop.

## 3. Wyniki liczbowe (kluczowe)

### 3.1 Fundament (G1, G2)

```text
G1: determinizm bitowy, zero NaN/runaway, H_Gamma nierosnace,
    kret zachowany (4x |n|=1, n_tot=0) na WSZYSTKICH runach: PASS
G2: |f_2D(r) / [s*(1 - xi^2/r^2)] - 1| <= 0.256% dla r in [6,10]
    (prog 1%): PASS — ogon amplitudowy wiru jest ALGEBRAICZNY
    (~1/r^2), zrodlowany faza; to trwaly wynik strukturalny
```

### 3.2 Pomiar glowny (G3) — FAIL

```text
fit laczny (10 runow, 200 krokow): Cphi/(2piJ) = 1.373  (pasmo
  [0.8, 1.2])  -> G3 FAIL;  RSS = 109 (sygnal ~0.8!)
czulosc (poza kryterium, bez skolapsowanych d_meas<3):
  Cphi/(2piJ) = 1.454, B/B_pred = 3.89, Delta_AIC(Y-P) = +1.7
  — dekompozycja nie dziala rowniez po wykluczeniu kolapsu
floor (P2, stabilnosc 150 vs 200) = 4.43 >> rozstep sygnalu 0.76
```

### 3.3 Diagnoza niewaznosci (produkt cyklu)

Trzy sprzezone przyczyny, wszystkie widoczne wprost w danych:

1. **Relaksacja lamie warunek symetrii konstrukcji.** Kasowanie fazy
   w E_sym wymaga IDENTYCZNYCH pozycji w S i O. Ale w S pary
   rownoimienne sie odpychaja (d_meas rosnie: 5 -> 6.19), a w O
   przyciagaja (d_meas maleje: 5 -> 0.85). Po 200 krokach konfiguracje
   nie sa juz "te same pozycje, inne znaki" — dokladna algebra
   kasowania przestaje obowiazywac, a gph_sym nabiera d-zaleznosci,
   ktora przecieka do kanalu parzystego.
2. **Kolaps par O przy malych d0.** d0=5 -> d_meas=0.85, d0=6 -> 2.90:
   pary schodza ponizej skali rdzenia (2xi ~ 1.5) w trakcie relaksacji
   — punkt opuszcza okno waznosci bazy punktowej Greena.
3. **Kanal fazowy przy d in [5,12] jest silnie modyfikowany rdzeniami**
   (zapowiedziane w Phase 0: E_anty/teoria = 1.34-1.61, dryf z d) —
   baza gph z punktowej funkcji Greena nie modeluje kanalu
   nieparzystego wystarczajaco dobrze, a jego rezyduum jest wieksze
   od sygnalu dylatacyjnego.

Wniosek metodologiczny: procedura "statyka po ustalonym tau_relax",
ktora wystarczala w B1 (sygnal ~6 jednostek), jest za gruba dla
kanalu dylatacyjnego (sygnal ~0.8) przy zrodlach, ktore sie PORUSZAJA
roznicowo.

### 3.4 Raport pozostalych pozycji (BEZ interpretacji, wg reguly stop)

```text
G4/G5/G6: wydruki w Phase2_output.txt (G5 formalnie "PASS (Y)" na
skorumpowanym ficie) sa NIEWAZNE na mocy zalockowanej reguly
"G3 FAIL -> stop przed interpretacja" — zapisane tu dla transparencji,
nie wnosza tresci fizycznej.
P5 (obserwacja): v*d pary rownoimiennej spada 1.38 -> 0.91 na d 4->15;
nieodroznialne od poprawek logarytmicznych mobilnosci (B1 widzialo
podobny dryf) — niekonkluzywne, raport wprost.
```

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- fundament teoretyczny kanalu dylatacyjnego stoi: ogon algebraiczny
  ~ 1 - xi^2/r^2 potwierdzony na 0.26%, predykcja B_pred = 13.55
  spojna miedzy siatkami — HIPOTEZA POZOSTAJE OTWARTA I POLICZALNA,
- pomiar wymaga procedury utrzymujacej identycznosc pozycji S/O:
  (a) ekstrapolacja tau_relax -> 0 z serii krotkich relaksacji,
  (b) wiezy symetrii (pomiar na ansatzu z relaksacja TYLKO amplitudy
  przy zamrozonej fazie — osobna, deklarowalna procedura), albo
  (c) wieksze d0 z gestszym skanem i dluzsza baza Greena — kazda
  z tych opcji to NOWY lock, nie poprawka tego.

CZEGO WYNIK NIE ZNACZY:
- NIE znaczy, ze kanal dylatacyjny nie istnieje albo ze odpycha —
  zaden z kierunkow nie zostal zmierzony (stop przed interpretacja),
- NIE podwaza wynikow B1 (tam sygnal fazowy >> floor i procedura
  byla adekwatna).

## 5. Errata / transparencja

1. Jedno doprecyzowanie pre-code (definicja floor w P2) dopisane
   PRZED pierwsza linijka kodu — udokumentowane w samym LOCKu.
2. Zadne progi, procedury ani parametry nie byly zmieniane po
   pierwszym runie. Fit czulosci (bez d_meas<3) raportowany wylacznie
   poza kryterium (precedens erraty #2 z B1); NIE zmienil werdyktu
   (G3 FAIL takze w czulosci).
3. Werdykt "STOP" pochodzi wprost z zalockowanej reguly decyzyjnej;
   G5-wydruk na skorumpowanym ficie oznaczony jako niewazny (3.4).

## 6. Decyzja i nastepny krok

Wg zalockowanej reguly ("G2/G3 FAIL -> stop przed interpretacja"):
cykl zamkniety bez rozstrzygniecia hipotezy. Nastepne kroki (kazdy
jako osobny op z nowym LOCKiem, do autoryzacji autora):

1. **Powtorka kanalu dylatacyjnego z procedura zachowujaca symetrie**
   (ekstrapolacja tau->0 albo relaksacja amplitudy przy zamrozonej
   fazie) — B_pred = 13.546 zostaje jako zalockowany zaklad.
2. **Stabilizacja obiektu** (op-wall-dynamics) — odblokowuje B2
   i daje stacjonarne tla takze dla innych pomiarow.
3. Dynamiczna wersja pomiaru dylatacyjnego (precyzyjny pomiar
   odchylki v*d od czystej fazy przy malych d) — wymaga najpierw
   ilosciowej teorii mobilnosci (osobny fundament).

Bilans mostu do grawitacji po pieciu opach: substrat ma kanal
dlugozasiegowy (fazowy, ladunkowy) i policzalna geometrie eikonalna
(czeka na tlo stacjonarne); uniwersalnosc znaku pozostaje otwartym
frontem programu — dwie drogi pomiarowe zdefiniowane, zadna jeszcze
nie rozstrzygnieta.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK (autoryzacja 2026-07-05)
- `Phase0_check.py` / `Phase0_output.txt` — algebra, ogon ODE,
  B_pred, walidacja algebry kasowania na ansatzu
- `Phase1_vortex_tail.py` / `Phase1_output.txt` — G2 (+G1) PASS
- `Phase2_channel_decomposition.py` / `Phase2_output.txt` — G3 FAIL,
  diagnoza; czulosc w tym pliku (sekcja 3.2)
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
