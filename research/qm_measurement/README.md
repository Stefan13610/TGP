# Q1: Niepewnosc pomiarowa z samozwrotnosci Phi

## Problem

W TGP czastka jest solitonem pola Phi, ktore JEST przestrzenia.
Pomiar pozycji czastki wymaga innego pola Phi (detektora).
Ale pola Phi interferuja — pomiar ZMIENIA to, co mierzymy.

To nie jest ograniczenie techniczne. To jest ONTOLOGICZNE:
nie istnieje zewnetrzny uklad odniesienia niezalezny od Phi.

## Obecny stan w TGP

Propozycja 3.3 (sek03_rezimy.tex): schematyczny dowod order-of-magnitude:
```
E_kin = hbar^2 / (2m * Dx^2)       <-- lokalizacja wymaga energii
E_field = beta*qm^2 / (4pi*Phi_0^2) * Dx   <-- pole Phi rosnie z Dx
Minimalizacja E_total => Dx* ~ lambda_C
=> Dx * Dp >= hbar(Phi)
```

Problem: to daje SKALE niepewnosci, ale nie MECHANIZM.
Nie wynika z tego rozklad prawdopodobienstwa.

## Nowe podejscie: pomiar jako interakcja soliton-soliton

### Fizyka

Detektor = duzy soliton o amplitudzie g0_det >> 1
Czastka = maly soliton o amplitudzie g0_part ~ 0.87

Kazdy soliton ma ogon: g(r) ~ 1 + A*sin(r+delta)/r

Gdy sa blisko, ogony sie nakladaja:
```
g_total(x) = g_det(x - x_det) + g_part(x - x_part) - 1
                                                    ^-- odejmij vacuum
```

Energia interakcji zalezy od fazy wzglednej ogonow:
```
E_int ~ integral A_det * A_part * cos(|x_det - x_part| + fazy) / r^2 dr
```

To OSCYLUJE z odlegloscia! => wynik pomiaru zalezy od fazy

### Kluczowy mechanizm

Faza ogona zalezy od g0 (amplituda solitonu).
Ale g0 jest zmienione przez obecnosc drugiego solitonu.
=> samozwrotnosc: pomiar zmienia stan

Niepewnosc wynika z tego, ze nie mozna znac fazy z dokladnoscia
lepsza niz ~1 (ogon oscyluje z okresem 2pi w jednostkach naturalnych).
To daje: Dx ~ 2pi / k ~ lambda_Compton

## Plan ataku

### Krok 1: Dwa solitony w 1D (q1_self_referential.py)
- Rozwiaz ODE dla pojedynczego solitonu
- Naloz dwa solitony w roznych odleglosciach
- Oblicz energie interakcji E_int(d) vs odleglosc d
- Sprawdz: czy E_int oscyluje? Z jakim okresem?

### Krok 2: Back-reaction (q1_back_reaction.py)
- Rozwiaz PELNE rownanie z dwoma solitonami
- Sprawdz jak obecnosc detektora zmienia profil czastki
- Oblicz Delta_g0 jako funkcje odleglosci
- => to jest "perturbacja pomiaru"

### Krok 3: Rozklad wynikow (q1_measurement_distribution.py)
- Detektor w losowej fazie (termiczny szum substratu)
- Zbierz histogram "wynikow pomiaru" (E_int lub Delta_g0)
- Sprawdz: czy rozklad ~ |A_tail|^2?
- Jesli tak => REGULA BORNA z TGP!

### Krok 4: Analityczne ograniczenie (q1_uncertainty_bound.py)
- Z rozwiazania perturbacyjnego ogonow
- Wyprowadz Dx * Dp >= hbar(Phi) FORMALNIE
- Porowna z Propozycja 3.3

## Kryterium zamkniecia

- [ ] E_int(d) oscyluje z okresem lambda_C
- [ ] Back-reaction Delta_g0 ~ A_det * A_part
- [ ] Rozklad wynikow ~ |A_tail|^2
- [ ] Formalna nierownosc Dx*Dp >= hbar z samozwrotnosci

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q1_self_referential.py | Interakcja dwoch solitonow | W TRAKCIE |
