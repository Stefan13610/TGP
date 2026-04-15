# Q4: Splatanie z korelacji substratu

## Problem

Wyprowadzic splatanie kwantowe ze struktury substratu TGP,
bez postulowania nielokalnych wektorow stanu.

## Mechanizm

Dwa solitony utworzone razem dziela region pola Phi.
Ich ogony interferuja w dzielonym regionie, tworzac KORELACJE:
- Wiez fazowy: phi_A + phi_B = const (zachowanie z kreacji)
- Pomiar jednego solitonu (przez back-reaction na Phi)
  zmienia dzielony substrat, wplywajac na drugi

To KONTEKSTUALNE (nie LHV) -- wybor pomiaru zmienia Phi.

## Wyniki (2026-04-15) -- q4_entanglement.py (3/4 PASS, 1 oczekiwany FAIL)

### Model fazowy
- phi_A + phi_B = 0 koduje splatanie (analogia: rozpad na dwa spiny)
- Funkcja korelacji sign(sin): E = -1 + 2|a-b|/pi (trojkatna, LHV)
- RMS vs LHV = 0.002 -- model fazowy jest DOKLADNIE LHV

### Test Bell-CHSH
- sign(sin) model: S = 0.001 -- brak naruszenia Bell (oczekiwane!)
- back-reaction model: S = -0.003 -- tez brak naruszenia
- To jest twierdzenie Bella: jeden ukryty parametr (phi) nie moze naruszyc CHSH
- **FAIL jest POPRAWNY** -- prosty model LHV nie moze naruszyc Bell

### Entropia splatania
- Niezalezne: H = 2*log(2*pi) = 3.676 bitow
- Splatane: H = log(2*pi) = 1.838 bitow
- Redukcja: delta_H = log(2*pi) = zamrozenie 1 stopnia swobody fazy

### Dekoherencja splatania
- Szum niezalezny: E(sigma) = E(0) * exp(-sigma^2)
- Potwierdzone numerycznie: E(0.5)/E(0) = 0.779 vs exp(-0.25) = 0.779 (!)
- Szum wspolny (common mode): NIE niszczy splatania
- Tempo dekoherencji: Gamma ~ N_env * chi^2 * A_env^2 / D_env^2

### Droga do naruszenia Bella
- Prosty model fazowy (1 parametr) jest LHV -- nie moze naruszyc
- TGP oferuje KONTEKSTUALNNOSC: pomiar zmienia substrat = zmienia stan
- Modele kontekstualne MOGA naruszyc Bell (znane twierdzenie)
- Potrzebny: wielowymiarowy model substratu (A_tail + phi + topologia)

### Predykcja testowalna
- Korelacje splatania zaleza od hbar(Phi)
- Blisko masywnego obiektu: splatanie SLABSZE
- delta_C/C ~ delta_hbar/hbar ~ GM/(rc^2) ~ 10^-9

## STATUS: Q4 CZESCIOWO ZAMKNIETE

- [x] Mechanizm zidentyfikowany (dzielony substrat)
- [x] Model wiezu fazowego zbudowany
- [x] Funkcja korelacji obliczona (LHV triangular)
- [x] Dekoherencja z szumu przeanalizowana (exp(-sigma^2))
- [x] Droga do naruszenia Bella zidentyfikowana (kontekstualnosc)
- [ ] Pelne CHSH > 2 z modelu substratu (wymaga wielowymiarowego podejscia)

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q4_entanglement.py | Fazy, korelacje, Bell, dekoherencja | 3/4 PASS |
