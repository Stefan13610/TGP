---
title: "op-asymmetric-lattice (B2-prime, proba #4) — pelny P2: uniwersalnosc w znaku kretu na siatce skosnej (eksperyment zerowy kanalu cyrkulacyjnego)"
date: 2026-07-05
type: research-op
status: CLOSED (STOP na bramce G2, 2026-07-12)
tags: [gravity-bridge, branch-B2prime, oblique-lattice, shear-instability, gate-stop, circulation-channel, universality, P2, L256]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-shortwave-lattice-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
---

# op-asymmetric-lattice-2026-07-05 (B2-prime, proba #4)

**Status: CLOSED — STOP na bramce G2 (2026-07-12), wg reguly
decyzyjnej LOCKa.** Autoryzacja autora 2026-07-12; LOCK + sekcja 10
(D1-D19) przed kodem; wykonawca: osobny agent, nowa sesja.

**Wynik:** siatka skosna (scinanie 64/rzad, L=256) jest dynamicznie
NIEZDATNA jako tlo: niestabilny mod scinajacy podsieci (+1 sunie
w +x, -1 w -x; gamma_s = +0.284/tau), utrata rdzeni (anihilacje
w rzedach) przy tau = 227.0 << wymagane 450. Sila netto na kazdym
wirze = 0 DOKLADNIE (pary inwersyjne v/-v, zweryfikowane), ale
konfiguracja jest SIODLEM — zero sily nie implikuje stabilnosci.
Scenariusz pre-rejestrowany jako RYZYKO #1; bramka zadzialala PRZED
jakimkolwiek runem pomiarowym. **Zero runow ugiecia; G3-G6 nie
mierzone; P2 pozostaje OTWARTE.** Deliverable: czasy zycia (skosna
227.0 vs szachownica L=256: 0.0000 przez tau=600 vs para: 107.5).
Kandydaci na nastepny op (osobne locki): scinanie s=32; L=512.
Szczegoly, erraty, interpretacja: [[Phase_FINAL_close.md]].

**Kontynuacja (autoryzowana 2026-07-12): proba #5 lamie lustro
GEOMETRIA WIAZKI (atan(1/2)) na zwalidowanej szachownicy zamiast
modyfikowac tlo** —
[[../op-oblique-beam-2026-07-12/Phase0_balance.md]].
**Werdykt proby #5 (2026-07-14): KANAL CYRKULACYJNY ZMIERZONY**
(G5b FAIL 3/4 par, spojna ujemna nadwyzka -0.98..-3.25 deg
ponad 3 sigma_Delta i rozrzut plaszczyzn; G6a PASS) — cel tego
opa (pomiar P2 bez oslony symetrii) osiagniety trzecia proba,
z wynikiem NEGATYWNYM dla "slepoty na kret"; eksperyment zerowy
przewidziany tutaj dal nie-zero na geometrii wiazki skosnej.
Szczegoly:
[[../op-oblique-beam-2026-07-12/Phase_FINAL_close.md]].

Pytanie: **czy kanal geometryczny jest slepy na znak kretu tam,
gdzie zadna symetria tego nie wymusza?** W probach #2/#3
alpha_odd = 0 bylo wlasnoscia KONSTRUKCJI tla (lustro+sprzezenie).
Ten cykl buduje tlo, ktore rozdziela kanaly: **siatka skosna**
(szachownica scinana o 64/rzad) ma pozycje soczewek lustrzane
wzgledem osi wiazki (eikonal, slepy na kret: Delta_eik ~ 0), ale
wzor LADUNKOW lamie lustro calkowicie — sprzezenie pakietu
z cyrkulacja tla nie jest niczym skasowane. Pomiar
Delta = alpha(+b) - alpha(-b) na 5 parach to **eksperyment zerowy
typu Eotvosa**: Delta zgodna z eikonalem = uniwersalnosc
ZMIERZONA; nadwyzka spojna w >= 2 parach = ODKRYCIE kanalu
cyrkulacyjnego (analogia grawitomagnetyzmu; bez deklaracji GR).
Sila netto na kazdym wirze = 0 DOKLADNIE z symetrii inwersyjnej
(nie D4); czas zycia rozstrzyga bramka G2.

Szczegoly, kryteria G1-G6, regula decyzyjna: [[Phase0_balance.md]].

## Notatki wykonawcze dla nowego agenta (doswiadczenie sesji #B2prime-1/2/3)

**Protokol (bezwzgledny):**
1. NAJPIERW autoryzacja autora -> status LOCKED -> doprecyzowania
   pre-code do sekcji 10 -> DOPIERO kod. Wzorzec: sekcja 10 LOCKa
   `../op-shortwave-lattice-2026-07-05/` (D1-D16); zmiany wskazane
   w sekcji 9 draftu (bilans sil PARAMI INWERSYJNYMI, detektor
   4+4, pakiet lewobiezny dla S-inv, definicje Delta per para).
2. Zadnych zmian progow/procedur po pierwszym runie; erraty
   w dniu runow, przed interpretacja. Bramki G2 nie wolno oslabic;
   przy G2 FAIL zadnych runow ugiecia.
3. Kazda faza: skrypt + output.txt w katalogu opa; na koncu
   Phase_FINAL_close.md + aktualizacja tego README.
4. KOLEJNOSC FAZ = CZESC PROTOKOLU: bramka G2 PRZED jakimkolwiek
   runem z pakietem.

**Co juz dziala (zwalidowane w probach #2/#3 — przenosic wprost):**
- caly aparat pomiarowy: `make_lattice` z par poziomych
  `theta_h_pair` (tu: 4 pary — sekcja 2 draftu), lockstep
  z psi_ref, wyzwolenie centroidem + okno zamrazane, pomiar
  wektorowy (Px,Py), b_eff, ray tracing z plaszczyzna centroidu
  3. probki. Kod wzorcowy: `../op-shortwave-lattice-2026-07-05/
  Phase2_deflection.py` (N=512, L=256 — ten sam grid).
- detektor rdzeni: klastrowanie windingu (D14) — rozszerzyc
  z 2+2 na 4+4 (DOKLADNIE 4 klastry na znak; nominaly w sekcji 0
  draftu).
- G3 na tym gridzie zmierzone w probie #3: bledy do -1.73%
  przy k=1.669 — powtorzyc pomiar (tani), nie cytowac.
- dip b_eff: b_nom - b_eff in [0.27, 0.84] (obie proby); margines
  0.9 w projekcie punktow.

**Znane pulapki specyficzne dla TEGO cyklu:**
- **bilans sil**: na siatce skosnej powloki-kwadraty (D4) NIE
  kasuja sie osobno — kasowanie jest PARAMI v/-v (inwersja).
  Sumowac parami; nie kopiowac petli powlokowej z prob #2/#3
  bez zmiany (bledny "brak kasowania" to bylby artefakt metody,
  analogia erraty #1 proby #2).
- **mody scinajace**: siatka skosna moze byc dynamicznie miekksza
  niz szachownica — G2 FAIL jest realnym scenariuszem; wtedy STOP
  i deliverable czasow zycia (kandydaci: scinanie 32, L=512).
- **artefakt detektora ~0.5 w b_eff miedzy +b/-b** (erraty #4/#1
  prob #2/#3): dlatego Delta_eik liczy sie z predykcji przy
  WLASNYM b_eff kazdego runu — artefakt sie odejmuje. Nie
  raportowac surowej Delta_meas bez Delta_eik (forbidden move #5).
- **pakiet lewobiezny (S-inv)**: nowy wariant kodu (start
  (-28, -64-b), kierunek -x, okno lustrzane, psi(-dt) z faza
  -omega dt po odbiciu — doprecyzowac w sekcji 10 PRZED kodem);
  uzywany WYLACZNIE w G6a, nie w pomiarze.
- **sigma_Delta**: najczystsze zero daje seria omega=1.5
  (sigma ~ 0.02-0.05 deg); seria 1.1 ma sigma do 0.2 deg —
  G5b tam jest slabszym testem (pre-rejestrowane w 1.3).
- **Narzedzia**: nie uzywac `cd`; sciezki od korzenia vaulta;
  outputy przez `> plik.txt`; python + numpy; dluzsze runy w tle.

**Kluczowe liczby (cytowane):**
- substrat: s* = 1.174166, V''(s*) = 0.878665, m_TV = 0.9374,
  A_t = 2.5690
- k0/lambda/v_g: (1.1): 0.814/7.72/0.370; (1.5): 1.656/3.79/0.552
- formula ogona (kalibracja: trace #3, blad 3%):
  t_deg(d) = 698/d^2 (1.1); 169/d^2 (1.5)
- proba #3, sigma pomiaru: 0.196 (1.1,14) ... 0.014 (1.5,16)
- G2 kontrast: szachownica L=256: 0.0000 przez tau=600;
  para: anihilacja tau=107.5
- P1 (zamkniete): 11/12 punktow kryterialnych w pasmie [0.5,2.0];
  19/19 runow znak KU wirowi
