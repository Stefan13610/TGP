---
title: "Phase_FINAL_close — zamkniecie: B2-prime proba #3 — refrakcja krotkofalowa na siatce L=256 (domkniecie P1 w domenie eikonalu)"
date: 2026-07-05
type: phase-final-close
tgp_owner: research/op-shortwave-lattice-2026-07-05
status: CLOSED
verdict: "G1 PASS; G2 BRAMKA PASS (0.0000 przemieszczenia przez tau=600, dryft energii 0.00e+00); G3 PASS (1.73% w pasmie rozszerzonym [0.5,1.7]); G4 PASS; G5 PASS 6/6 — WSZYSTKIE punkty kryterialne w pasmie [0.5,2.0] ze znakiem KU wirowi -> P1 ROZSTRZYGNIETE POZYTYWNIE (lacznie z 5/6 proby #2); G6a PASS (4.9e-11 deg / bitowo); G6b PASS (roznica 0.000 deg)"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, shortwave, L256, closed, P1-resolved]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G5 PASS 6/6 — PYTANIE P1 (POLICZALNOSC) ROZSTRZYGNIETE
POZYTYWNIE.** Wszystkie szesc punktow kryterialnych,
zaprojektowanych z gory w domenie stosowalnosci eikonalu
(b_nom >= 1.5 lambda + 0.9), daje ugiecie KU wirowi o wielkosci
zgodnej z eikonalem substratu w pasmie [0.5, 2.0]:
ratio = 0.788-1.573. Lacznie z proba #2 (5/6; wypadl wylacznie
punkt dyfrakcyjny b_eff ~ lambda) program ma **11 waznych punktow
kryterialnych w pasmie na dwoch geometriach siatki (L=128
i L=256), w zakresie omega 1.1-1.5 i b_eff/lambda 1.5-4.2** —
policzalna, uniwersalnie skupiajaca geometria propagacji wokol
materii-defektow, z jawnie zmierzona granica stosowalnosci
eikonalu (b_eff ~ lambda; obserwacja 1.1/8: ratio 0.398 przy
L=256 vs 0.445 przy L=128 — wlasnosc rezimu falowego, nie
geometrii siatki). Tlo szachownicowe skaluje sie do L=256 bez
degradacji (G2: 0.0000 przez tau=600). Superpozycja soczewek
zaobserwowana bezposrednio (D16 wykonana: pakiet za druga
soczewka, kat surowy P ujemny = ugiecie ku drugiemu wirowi rzedu).

Glowne USTALENIA cyklu:

1. **P1 domkniete pomiarem, nie deklaracja**: eikonal liczony
   z samego substratu (n^2 z V'' na zrelaksowanym tle) przewiduje
   kat ugiecia pakietu amplitudowego w calej swojej domenie
   stosowalnosci do czynnika <= 1.6, z poprawnym znakiem
   (KU defektowi) w 100% runow (10/10 waznych w tym cyklu,
   19/19 lacznie z proba #2).
2. **Granica stosowalnosci eikonalu POTWIERDZONA niezaleznie od
   siatki**: punkt (1.1, b=8) na L=256 daje ratio 0.398 przy
   b_eff/lambda = 0.91 (L=128: 0.445 przy 0.93) — rozjazd
   dyfrakcyjny jest wlasnoscia fizyki falowej (fala usrednia
   soczewke po swojej szerokosci), nie artefaktem geometrii.
   Wraz z (1.3, 8): 0.718 (L=128: 0.741) i (1.7, 8): 0.804 —
   przejscie do pasma jest monotoniczne w b_eff/lambda.
3. **Tlo szachownicowe jest skalowalna konstrukcja pomiarowa**:
   przy L=256 (rozstaw soczewek 128) rdzenie nie przemieszczaja
   sie W OGOLE pod pelna dynamika II rzedu przez tau=600
   (kryterium: 450 >= 3x przelot); residuum relaksacji 9.6e-07;
   symetrie tla dokladne do 4.9e-11 deg (S1) / bitowo (S2).
4. **Dyspersja sieciowa pod kontrola przy k0(1.5) = 1.66**:
   zmierzony blad -1.73% (analitycznie: -1.72%) — siatka dx=0.5
   niesie rzetelnie kanal krotkofalowy; omega=1.7 (k0*dx ~ 1.0,
   blad ~ -2.9%) pozostaje poza kryteriami, ale obserwacyjnie
   tez daje ratio w pasmie (0.804).

## 1. Co testowano

Sciezka reguly decyzyjnej proby #2 ("G5 FAIL wielkosc przy dobrym
znaku -> powtorka krotkofalowa jako osobny op"): to samo pytanie
P1, ale KAZDY punkt kryterialny zaprojektowany z gory w domenie
stosowalnosci eikonalu (b_nom >= 1.5 lambda + 0.9; margines 0.9
z pomiaru dipu b_eff w probie #2), dwiema dzwigniami: seria
omega=1.5 (lambda=3.79) i L=256 (wieksze b dla omega=1.1 bez
nakladania soczewek). Model bez strojenia (kappa=0.5, a=0.5,
b=1.6, c=1.0; N=512, L=256, dx=0.5 bez zmian); tlo: szachownica
4 wirow (+1,-1/-1,+1) w rozstawie L/2=128; bramka G2
przeskalowana do nowego przelotu (tau<=450, okno 600); pasmo
kalibracji G3 rozszerzone do [0.5, 1.7] z pre-rejestrowanym
oszacowaniem bledu sieciowego.

## 2. Przebieg

- **LOCK**: autoryzacja autora (polecenie realizacji cyklu,
  2026-07-05); doprecyzowania pre-code D1-D16 dopisane do sekcji
  10 PRZED pierwsza linijka kodu (adaptacja geometryczna wzorca
  proby #2; D16 doprecyzowana: wynik kryterialny punktu (1.5,8)
  z runu standardowego D9, superpozycja z osobnego runu
  dedykowanego).
- **Phase 0** (`Phase0_output.txt`): tabela projektowa — 6/6
  punktow w domenie (zapas +1.41 do +9.41 nad progiem
  1.5 lambda + 0.9); bilans sil L=256: kazda powloka D4 = 0.0
  DOKLADNIE (obciecie pudelkowe: 9.6e-06 ~K^-2, artefakt);
  dyspersja sieciowa analitycznie: -1.72% przy k0(1.5), -2.86%
  przy k0(1.7); tlo zrelaksowane: residuum 9.6e-07, szew 0.124
  rad; zaklad eikonalny WIAZACY na tle L=256 + czulosc plaszczyzn;
  budzet czasu (przelot 140 przy omega=1.1; superpozycja ~315
  -> limit 500 wystarcza).
- **Phase 1** (`Phase1_output.txt`): **BRAMKA G2 PASS**: max
  przemieszczenie rdzeni 0.0000 (tau<=450 i tau<=600), n_tot=0,
  4 rdzenie przez cale okno, dryft energii 0.00e+00, gamma_s:
  brak widocznej ucieczki; **G3 PASS**: max |blad| 1.73%
  (k=1.669; prog 5%); determinizm bitowy (relaksacja + fala).
  G1(czesc) PASS.
- **Phase 2** (`Phase2_output.txt`): 11 runow, WSZYSTKIE WAZNE
  (5 probek kazdy; zadnych flag; E_dryft <= 2.2e-14; gamma_core
  kryterialnych <= -0.038; kret zachowany). G4 PASS; **G5 PASS
  6/6**; G6b PASS (roznica 0.000 deg przy u0/2); obserwacje
  1.1/8, 1.3/8, 1.7/8 + superpozycja D16 WYKONANA.
- **Phase 3** (`Phase3_output.txt`): G6a PASS — b=-14 zgodne do
  4.9e-11 deg; tlo lustrzane: pole po relaksacji ROWNE sprzezeniu
  (max roznica 0.0e+00), run bitowo identyczny (D11 potwierdzona
  co do bitu). alpha_odd = +2.5e-11 deg = 0 z konstrukcji
  (NIE raportowane jako pomiar P2, forbidden move #5).

## 3. Wyniki liczbowe (kluczowe)

```text
BRAMKA G2 (tlo bez pakietu, dynamika II rzedu, tau=600, L=256):
  max przemieszczenie rdzenia: 0.0000 (wszystkie 4, cale okno)
  dryft sekularny energii: 0.00e+00; gamma_s: brak ucieczki
G3 (dyspersja, pasmo rozszerzone [0.5, 1.7]):
  max |blad| = 1.73% przy k=1.669 (analitycznie -1.72%)  PASS
skan glowny (alpha zmierzone vs alpha_pred(b_eff) [deg]; ratio;
b_eff/lambda) — zestaw kryterialny ZAMKNIETY, 6 punktow:
  omega=1.1 b=14: +5.555+/-0.196  vs +3.765  ratio 1.475  (1.72)  OK
  omega=1.1 b=16: +4.177+/-0.160  vs +2.656  ratio 1.573  (2.00)  OK
  omega=1.1 b=20: +2.168+/-0.064  vs +1.510  ratio 1.435  (2.54)  OK
  omega=1.5 b= 8: +2.668+/-0.035  vs +3.385  ratio 0.788  (2.01)  OK
  omega=1.5 b=12: +1.765+/-0.048  vs +1.215  ratio 1.453  (3.08)  OK
  omega=1.5 b=16: +0.835+/-0.014  vs +0.617  ratio 1.354  (4.15)  OK
G4: |alpha(14,1.1)| = 5.55 > 3 sigma = 0.59; monotonia w b: TAK;
    |alpha(1.1,16)| > |alpha(1.5,16)|: TAK
G6a: |alpha(+14)-alpha(-14)| = 4.9e-11 deg; lustro: 0.0 (bitowo)
G6b: |alpha(u0)-alpha(u0/2)| = 0.000 deg (prog 0.556)
obserwacje poza kryteriami (D12):
  1.1/8: alpha=+8.62, pred=+21.66, ratio 0.398, b_eff/lambda=0.91
         (L=128: 0.445 przy 0.93 -> wlasnosc rezimu, nie siatki)
  1.3/8: alpha=+4.49, pred=+6.24,  ratio 0.718, b_eff/lambda=1.51
         (L=128: 0.741 — sprzeg miedzy probami spojny)
  1.7/8: alpha=+1.80, pred=+2.24,  ratio 0.804, b_eff/lambda=2.45
         (ultra-krotkofalowo; dyspersja sieciowa ~ -2.9%)
  superpozycja (D16, WYKONANA — proba #2 nie osiagnela okna):
    tau=347, frakcja 0.48 za druga soczewka; centroid y = -53.3;
    szerokosc y = 37.9; kat surowy P = -2.91 deg (znak ku drugiemu
    wirowi rzedu — jakosciowo zgodny z geometria dwusoczewkowa)
```

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- **P1 ROZSTRZYGNIETE POZYTYWNIE** (w brzmieniu LOCKa, lacznie
  z 5/6 proby #2): kat ugiecia pakietu amplitudowego na defekcie
  jest POLICZALNY z samego substratu (eikonal n^2 z V'') w calej
  domenie stosowalnosci eikonalu, na stabilnym tle
  wielosoczewkowym, w dwoch geometriach siatki,
- **kanal geometryczny ma zmierzony znak grawitacji
  i policzalnosc**: 19/19 waznych runow ugiecia w programie
  (proba #2 + #3) daje ugiecie KU defektowi; wielkosc w pasmie
  [0.5, 2.0] wszedzie tam, gdzie b_eff >= 1.5 lambda,
- **granica stosowalnosci eikonalu jest zmierzona i uniwersalna**:
  ratio spada ponizej pasma wylacznie dla b_eff ~ lambda
  (0.398-0.445 niezaleznie od L) — spojne z usrednianiem soczewki
  przez fale o szerokosci ~lambda,
- **systematyka w pasmie**: punkty dalekie (b_eff >= 1.7 lambda)
  daja ratio 1.35-1.57 (pomiar > eikonal), punkt najkrotszy
  w domenie (1.5, 8; 2.0 lambda) daje 0.79 — pasmo [0.5, 2.0]
  pozostaje wlasciwa miara policzalnosci rzedu 1; struktura
  subtelniejsza (skad nadwyzka ~1.4-1.6 przy omega=1.1) = materia
  na przyszle opy, nie do reinterpretacji tutaj.

CZEGO WYNIK NIE ZNACZY:
- NIE testuje uniwersalnosci w znaku kretu (P2) nietrywialnie:
  alpha_odd = 0 z konstrukcji tla (pre-rejestrowane; pelny P2 =
  osobny op na tle asymetrycznym),
- NIE mierzy skalowania ugiecia z |n|; nie wyprowadza metryki;
- zadnych deklaracji GR/soczewkowania obserwacyjnego; tau != czas;
  bez promocji do core/ (osobna autoryzacja).

## 5. Errata / transparencja

1. **b_eff w wariancie S1 (b=-14)**: 13.785 vs 13.285 (baseline) —
   roznica 0.50 = 2 x 0.25 przesuniecia detektora plakietowego
   wzgledem osi odbicia (jak errata #4 proby #2); alpha zgodne do
   4.9e-11 deg — wlasnosc detektora, nie dynamiki.
2. **gamma_core runu obserwacyjnego D16**: +0.0102 (prog
   kryterialny 0.01) — run POZA kryteriami (obserwacja
   superpozycji, dluga kontynuacja tau=347 z pakietem
   przechodzacym przy rdzeniach); wszystkie runy KRYTERIALNE
   maja gamma_core ujemne (<= -0.038). Raportowane dla pelnej
   przejrzystosci; bez wplywu na G1 (ewaluacja G1 obejmuje runy
   kryterialne + G6b, zgodnie z LOCKiem).
3. **Probki alpha G6b bitowo identyczne z baseline** (roznica
   0.000 deg do 4+ miejsc): liniowosc sektora delta_psi przy
   u0 = 1e-3 jest dokladna do precyzji pomiaru — spojne z proba
   #2; raportowane jako wlasnosc, nie anomalia.
4. Zadnych zmian progow, parametrow ani procedur po pierwszym
   runie; zadnych powtorek runow kryterialnych; zestaw G5
   pozostal zamkniety (obserwacje 1.1/8, 1.3/8, 1.7/8
   niedolaczone).

## 6. Decyzja i nastepny krok

Wg reguly LOCKa — sciezka "G5 PASS (+ G4, G6 PASS)":
**P1 ROZSTRZYGNIETE POZYTYWNIE** (lacznie z 5/6 proby #2):
policzalna, uniwersalnie skupiajaca geometria propagacji wokol
materii-defektow, z jawnie zmierzona granica stosowalnosci
eikonalu (b_eff ~ lambda). Nastepne opy (kazdy za OSOBNA
autoryzacja i osobnym lockiem):

1. **pelny P2 na tle asymetrycznym** — uniwersalnosc w znaku
   kretu; czesc cyrkulacyjna nieskasowana symetria (priorytet
   wg reguly),
2. **skalowanie z |n|** — czy ugiecie ~ n^2 (defekty wyzszego
   kretu),
3. **granica newtonowska 2D i G_eff z parametrow substratu** —
   pierwszy krok ilosciowy od refrakcji do "grawitacji" programu.

Bilans mostu do grawitacji po dziewieciu opach: kanal
geometryczny wokol materii-defektow ma (a) stabilne, skalowalne
tlo pomiarowe (L=128 i L=256), (b) zmierzony uniwersalnie
przyciagajacy znak (19/19), (c) POLICZALNOSC POTWIERDZONA
w calej domenie eikonalu (11 punktow kryterialnych w pasmie,
0 wyjatkow w domenie), (d) zmierzona granice stosowalnosci
(b_eff ~ lambda) i (e) bezposrednia obserwacje przejscia pakietu
za druga soczewke. P1 ZAMKNIETE.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK (sekcja 10: D1-D16)
- `Phase0_check.py` / `Phase0_output.txt` — tabela projektowa
  (6/6 w domenie), bilans sil L=256 (0.0 per powloka), dyspersja
  sieciowa analityczna, tlo (residuum 9.6e-07), eikonal wiazacy
  + czulosc plaszczyzn, budzet czasu
- `Phase1_gate_and_calibration.py` / `Phase1_output.txt` —
  **BRAMKA G2 PASS** (0.0000 / tau=600), G3 PASS (1.73%
  w [0.5, 1.7]), background_check, determinizm
- `Phase2_deflection.py` / `Phase2_output.txt` — 11/11 runow
  waznych; G4 PASS, **G5 PASS 6/6**, G6b PASS + obserwacje
  (w tym superpozycja D16 wykonana)
- `Phase3_symmetry.py` / `Phase3_output.txt` — G6a PASS
  (S1: 4.9e-11 deg; S2: bitowo 0)
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
