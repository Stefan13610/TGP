---
title: "Phase_FINAL_close — zamkniecie: tlo stacjonarne dla sektora falowego (babel vs wir)"
date: 2026-07-05
type: phase-final-close
tgp_owner: research/op-stationary-background-2026-07-05
status: CLOSED
verdict: "G1 PASS; G2 FAIL (noga lambda, R_c potwierdzone); G3 FAIL strukturalne (para na torusie zjezdza z maksimum); G4 formalnie FAIL (noga globalna przez ruch tla) — ALE sektor falowy wokol wiru BEZ TACHIONU (lambda_full=+0.0008, lambda_toy=+0.0198, gamma_rdzenie=0.0003); STOP wg reguly"
anti_lakatos_lock: PRESERVED
tags: [gravity-bridge, stationary-background, critical-bubble, vortex, tachyon-resolution, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[../op-phi-metric-refraction-2026-07-04/Phase_FINAL_close.md]]"
  - "[[../op-goldstone-mediator-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamkniecie cyklu

## 0. Werdykt

**G1 PASS. G2 FAIL (formalnie: noga lambda_min; R_c potwierdzone
w pasmie). G3 FAIL — strukturalne, utrzymane po pre-rejestrowanym
wariancie x4 -> STOP wg reguly. G4 formalnie FAIL (noga globalna G4b
przez ruch tla), przy czym obie miary TACHIONOWE sa czyste.
G5 (klasyfikator): pakiet ugina sie KU wirowi.**

Dwa glowne USTALENIA cyklu (przetrwaja werdykt formalny):

1. **Sektor falowy wokol wiru NIE MA tachionu** — w pelnej
   linearyzacji zespolonej `lambda_min = +0.0008` (>= -1e-3),
   empirycznie `gamma(rdzenie) = 0.0003` (B2 na obiekcie mialo
   0.125!), a WBREW pre-code oczekiwaniu takze toy-sektor
   rzeczywisty jest czysty: `lambda_min(L_toy) = +0.0198 > 0`
   (mod rdzeniowy DODATNI). Wniosek strukturalny: tachion B2 byl
   wlasnoscia SCIANY duzego rosnacego obiektu (R=12: dlugi, niemal
   prosty pierscien V''<0 wiaze stan), a NIE sektora falowego ani
   substratu: pierscien wokol rdzenia wiru (r ~ 0.5-2.3) jest za
   krotki i za mocno zakrzywiony, by wiazac. Falsyfikacja wlasnego
   oczekiwania pre-code zaraportowana wprost.
2. **Scisla stacjonarnosc jest przeszkoda topologiczno-geometryczna,
   nie falowa:** na torusie kret calkowity musi byc 0, wiec "wir"
   = para (+1,-1); konfiguracja maksymalnie symetryczna
   (szachownicowa) jest MAKSIMUM energii wzglednej — kazda
   asymetria inicjuje przyspieszajacy zjazd ku anihilacji
   (dryf 0.42 -> 1.21 na kolejnych 4000 krokow; wariant x4
   pre-rejestrowany: FAIL). Dryf jest jednak WOLNY: v_tla ~ 0.015
   na jednostke tau vs v_g pakietu ~ 0.37 (separacja skal ~25x).

## 1. Co testowano

Po B2 (tachion sciany zamrozonego obiektu): czy substrat posiada
tlo stacjonarne i falowo stabilne. Kandydaci: (K1) babel krytyczny,
(K2) wir/para (stacjonarnosc topologia — zgodnie z reinterpretacja
B1). Sektor falowy = II pochodna (jedyna struktura B2) na polu
zespolonym B1. Model bez strojenia, N=256, L=128.

## 2. Przebieg

- **Phase 0**: LOCK (autoryzacja: opcja 2 po zamknieciu galezi C).
  `Phase0_output.txt`: R_c_pred = 3.0626 (z CYTOWANEGO prawa frontu
  op-lock-interaction); n(r) wokol wiru > 1 wszedzie (deficyt =
  osrodek gestszy), ogon algebraiczny n^2-1 ~ 2.569/r^2/(omega^2-V'');
  **predykcja eikonalna: refrakcja KU wirowi, UNIWERSALNA (n slepe
  na znak kretu)**: alpha_pred(omega=1.1) = +43/+14/+5.3 deg dla
  b = 6/8/12 (b=4, omega=1.0: strefa przechwytu).
- **Phase 1**: bisekcja babla: **R_c = 3.227 +/- 0.04** (5.4% od
  predykcji, w pasmie) — dychotomia SHRINK/GROW = empiryczny dowod
  siodlowosci K1. Noga lambda_min: +0.0066 (patrz errata #3).
- **Phase 2**: tlo wirowe: widma L_toy/L_full, noise-run, demo
  pakietu; wariant x4 (Phase 2b). Wyniki w sekcji 0 i 3.

## 3. Wyniki liczbowe (kluczowe)

```text
K1 (babel):  R_c = 3.227 +/- 0.039  (pred. 3.063; pasmo OK)
             ponizej: kurczy sie; powyzej: rosnie (siodlo dynamiczne)
K2 (para wirow, ansatz zlozony — errata #1):
  |psi| daleko od rdzeni: min 1.018 (s* = 1.174; wadliwy ansatz: 0.10)
  lambda_min(L_toy)  = +0.01979  (mod 99.8% w rdzeniach, DODATNI)
  lambda_min(L_full) = +0.00079  (prog -1e-3: PASS; dno ~ Goldstone)
  noise-run tau=150: gamma(rdzenie) = +0.00028 (PASS),
    gamma(globalnie) = +0.0115 (prog 0.01: FAIL — ale to ruch tla
    /zjazd pary i relaksacja artefaktu C, nie wzrost modu falowego;
    dryft sekularny energii 6.9e-07, zero flag)
  stacjonarnosc: residuum 8.1e-03 (prog 1e-4: FAIL), dryf rdzeni
    0.42/1000 krokow -> wariant x4: dryf 1.21/4000 (przyspiesza) FAIL
demo pakietu (omega=1.1, b=6): kierunek KU wirowi (probki katowe
  -135 deg spojne; wielkosc >> eikonal +43 deg — pakiet szeroki
  wzgledem b, strefa silnego ugiecia; TYLKO klasyfikator kierunku)
```

## 4. Interpretacja w granicach modelu

CO WYNIK ZNACZY:
- **problem tachionowy B2 jest ROZWIAZANY zrozumieniem**: to sciana
  duzego rosnacego obiektu, nie sektor falowy; wokol defektow
  (wirow) sektor falowy — nawet minimalny toy — jest stabilny,
- **defekt jest soczewka skupiajaca**: n > 1 z ogonem 1/r^2,
  uniwersalna (slepa na znak kretu); demo potwierdza kierunek KU;
  to pierwszy kanal w programie o ZNAKU grawitacji w geometrii
  propagacji — czeka na pomiar ilosciowy,
- pozostala przeszkoda jest czysto techniczna: scisla stacjonarnosc
  pary na skonczonym torusie; separacja skal (25x) czyni tlo
  QUASI-stacjonarnym w oknie pomiaru pakietu.

CZEGO WYNIK NIE ZNACZY:
- nie deklarujemy "soczewkowania grawitacyjnego" (forbidden move 3;
  ilosciowy pomiar alpha(b, omega) = nastepny op),
- G4 formalnie FAIL — nie oglaszamy "odblokowania B2-prime" na mocy
  reguly; oglaszamy zmierzone przeslanki do nowego locka.

## 5. Errata / transparencja

1. **Wadliwy ansatz pierwszego runu Phase 2 (naprawa w dniu runow,
   przed interpretacja):** iloczyn sinh (D1 z B1) dla pary na
   PRZEKATNEJ zostawia skok fazy pi na szwie x (272 linki; dla par
   wspolliniowych z B1 asymptotyki kasuja sie — dla roznych y NIE).
   Relaksacja zamienila szew w sztuczna ciemna sciane (|psi| -> 0.10)
   z falszywym tachionem (lambda = -0.061..-0.079) SPRZECZNYM
   z empirycznym noise-runem (gamma ~ 0.003) — sprzecznosc wykryta,
   zdiagnozowana i usunieta. Poprawka: kompozycja dokladna dwoch par
   wspolliniowych (pozioma A->C + transponowana pionowa C->B; ladunek
   w C sie kasuje). Wadliwy output zachowany:
   `Phase2_output_run1_defective_ansatz.txt`. Pozycje wirow z LOCKa
   bez zmian; zero zmian progow.
2. **Rezydualny artefakt punktu C:** kompozycja zostawia lokalnie
   duzy gradient fazy przy C=(32,-32) (max 1.58 rad na linku, zero
   kretu netto) — relaksuje sie do plytkiego sladu (|psi|min = 1.02);
   wklad do residuum G3 i do gamma(globalnie).
3. **G2, noga lambda_min:** Hessian mierzony na ZAMROZONYM ansatzu
   tanh przy R_c z bisekcji (+0.0066) — do prawdziwego siodla nie da
   sie dorelaksowac przeplywem gradientowym, wiec noga jest
   niekonkluzywna metodologicznie; klasyfikacje K1 = siodlo
   rozstrzyga dychotomia bisekcji (kurczy/rosnie). Raport wprost;
   G2 formalnie FAIL zgodnie z brzmieniem LOCKa.
4. Pre-code oczekiwanie "L_toy ma tachion na tle wirowym"
   SFALSYFIKOWANE przez pomiar (+0.0198) — zaraportowane jako
   pelnoprawny wynik (tak dziala protokol: oczekiwania sa zapisane,
   dane rozstrzygaja).

## 6. Decyzja i nastepny krok

Wg reguly LOCKa (G3 FAIL po wariancie x4 -> STOP): cykl zamkniety
bez formalnego "odblokowania". Zmierzone przeslanki definiuja
nastepny op (za autoryzacja):

**B2-prime na tle quasi-stacjonarnym** — nowy lock z jawnie
pre-rejestrowanymi elementami: (i) tlo = para wirow, dryf tla
zmierzony TUTAJ (0.015/tau, 25x wolniej od pakietu) i skorygowany
w obserwablach; (ii) sektor falowy pelny zespolony (czysty:
lambda = +0.0008); (iii) zaklad eikonalny alpha_pred(b, omega)
z Phase0 tego cyklu; (iv) alternatywy geometryczne: L=256 (dryf
~4x mniejszy) albo relaksacja z wiezem symetrii obrotu o pi
(do rozstrzygniecia w locku, czy to nie narusza zakazu pinningu).

Bilans mostu do grawitacji po szesciu opach:
- kanal fazowy: dlugozasiegowy, policzalny, LADUNKOWY (B1),
- kanal dylatacyjny: otwarty (galaz C, procedura do wymiany),
- kanal geometryczny: tachion B2 zrozumiany i usuniety (sciana, nie
  sektor); wokol DEFEKTOW propagacja stabilna, defekt = soczewka
  skupiajaca (n > 1, ogon 1/r^2, uniwersalna) — najblizszy pomiaru
  kanal o znaku grawitacji.

## 7. Pliki cyklu

- [[Phase0_balance.md]] — LOCK
- `Phase0_check.py` / `Phase0_output.txt` — R_c_pred, n^2(r),
  alpha_pred, progi
- `Phase1_critical_bubble.py` / `Phase1_output.txt` — G2
- `Phase2_vortex_background.py` / `Phase2_output.txt` — G3, G4, G5
  (+ `Phase2_output_run1_defective_ansatz.txt` — errata #1)
- `Phase2b_relax_x4_variant.py` / `Phase2b_output.txt` — wariant x4
- [[README.md]] — przeglad opa
- ten plik — zamkniecie
