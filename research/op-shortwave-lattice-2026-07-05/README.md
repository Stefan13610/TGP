---
title: "op-shortwave-lattice (B2-prime, proba #3) — refrakcja krotkofalowa na siatce L=256: domkniecie P1 w domenie eikonalu"
date: 2026-07-05
type: research-op
status: CLOSED
verdict: "G1-G6 WSZYSTKIE PASS; G5 6/6 -> P1 ROZSTRZYGNIETE POZYTYWNIE (lacznie z 5/6 proby #2)"
tags: [gravity-bridge, branch-B2prime, lattice-background, vortex-lensing, refraction, eikonal, shortwave, L256, closed, P1-resolved]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-lattice-background-2026-07-05/Phase_FINAL_close.md]]"
  - "[[../op-vortex-refraction-2026-07-05/Phase_FINAL_close.md]]"
---

# op-shortwave-lattice-2026-07-05 (B2-prime, proba #3)

**Status: CLOSED (2026-07-05). Werdykt: G1 PASS; G2 BRAMKA PASS
(0.0000 przemieszczenia przez tau=600); G3 PASS (1.73% w pasmie
rozszerzonym [0.5, 1.7]); G4 PASS; G5 PASS 6/6; G6a PASS
(4.9e-11 deg / bitowo); G6b PASS (0.000 deg).
PYTANIE P1 (POLICZALNOSC) ROZSTRZYGNIETE POZYTYWNIE.**

Pytanie bylo: czy eikonal substratu przewiduje ugiecie ilosciowo
(pasmo [0.5, 2.0]) we WSZYSTKICH punktach zaprojektowanych w jego
domenie stosowalnosci (b_eff >= 1.5 lambda)? Odpowiedz: **TAK** —
6/6 punktow kryterialnych (omega=1.1 x b {14,16,20}; omega=1.5 x
b {8,12,16}) daje znak KU wirowi i ratio 0.788-1.573; lacznie
z proba #2: 11 punktow kryterialnych w pasmie na dwoch
geometriach siatki, w zakresie b_eff/lambda 1.5-4.2. Granica
stosowalnosci eikonalu potwierdzona niezaleznie od siatki
(obserwacja 1.1/8: ratio 0.398 przy L=256 vs 0.445 przy L=128,
b_eff ~ lambda). Superpozycja soczewek zaobserwowana bezposrednio
(D16 wykonana: tau=347, frakcja 0.48 za druga soczewka).

Szczegoly: [[Phase_FINAL_close.md]] (werdykt, liczby, errata,
nastepne kroki) i [[Phase0_balance.md]] (LOCK, D1-D16).

## Wyniki kluczowe (skrot)

```text
G5 (zestaw ZAMKNIETY, 6 punktow):
  omega=1.1 b=14: +5.555 vs +3.765  ratio 1.475  b_eff/lam 1.72
  omega=1.1 b=16: +4.177 vs +2.656  ratio 1.573  b_eff/lam 2.00
  omega=1.1 b=20: +2.168 vs +1.510  ratio 1.435  b_eff/lam 2.54
  omega=1.5 b= 8: +2.668 vs +3.385  ratio 0.788  b_eff/lam 2.01
  omega=1.5 b=12: +1.765 vs +1.215  ratio 1.453  b_eff/lam 3.08
  omega=1.5 b=16: +0.835 vs +0.617  ratio 1.354  b_eff/lam 4.15
obserwacje: 1.1/8: 0.398 (0.91 lam); 1.3/8: 0.718 (1.51 lam);
            1.7/8: 0.804 (2.45 lam)
```

## Nastepne opy (kazdy za osobna autoryzacja, osobny lock)

1. pelny P2 na tle asymetrycznym (uniwersalnosc w znaku kretu)
   — **draft gotowy**:
   [[../op-asymmetric-lattice-2026-07-05/Phase0_balance.md]]
   (siatka skosna: pozycje lustrzane, kret lamie lustro —
   eksperyment zerowy kanalu cyrkulacyjnego).
   **WERDYKT proby #4 (2026-07-12): STOP na bramce G2** — siatka
   skosna (scinanie 64/rzad) dynamicznie niezdatna: niestabilny
   mod scinajacy podsieci (gamma_s = +0.284/tau), utrata rdzeni
   przy tau = 227.0 << 450; zero runow ugiecia, P2 pozostaje
   OTWARTE; kandydaci na kolejna probe (osobne locki): scinanie
   s=32, L=512
   ([[../op-asymmetric-lattice-2026-07-05/Phase_FINAL_close.md]]).
   **Proba #5 (autoryzowana 2026-07-12): wiazka skosna atan(1/2)
   na zwalidowanej szachownicy** —
   [[../op-oblique-beam-2026-07-12/Phase0_balance.md]].
   **WERDYKT proby #5 (2026-07-14): KANAL CYRKULACYJNY
   ZMIERZONY** — G5b FAIL w 3/4 par ze spojna ujemna nadwyzka
   Delta_meas - Delta_eik (-0.98 do -3.25 deg vs 3 sigma_Delta
   0.24-0.43 i rozrzut plaszczyzn 0.09-0.37), G6a PASS (S-inv
   0.056 deg, S-conj bitowo — sygnal jest wlasnoscia pola, nie
   ramki); sygnatura w b_eff: +b odpychane, -b przyciagane.
   P2 w brzmieniu "slepota na kret" rozstrzygniete NEGATYWNIE;
   status P1 (ten op) BEZ ZMIAN — czesc parzysta ugiecia
   pozostaje w pasmie eikonalu we wszystkich 4 parach (post hoc);
   G4/G5a FAIL = skutek tego samego kanalu (tlumienie alpha na
   galezi +b serii 1.5 az do odwrocenia znaku). Nastepny op:
   charakteryzacja kanalu cyrkulacyjnego (osobny lock)
   ([[../op-oblique-beam-2026-07-12/Phase_FINAL_close.md]]).
2. skalowanie ugiecia z |n| (defekty wyzszego kretu),
3. granica newtonowska 2D i G_eff z parametrow substratu.

## Notatki wykonawcze (dla przyszlych opow; doswiadczenie sesji #B2prime-3)

**Protokol utrzymany w calosci**: autoryzacja -> LOCK -> sekcja 10
(D1-D16) -> kod; bramka G2 przed pierwszym pakietem; zero zmian
progow po pierwszym runie; zestaw G5 zamkniety (obserwacje
niedolaczone); erraty spisane w dniu runow (S1-offset detektora,
gamma_core runu D16 poza kryteriami).

**Co dziala i przenosi sie dalej (zwalidowane na L=128 i L=256):**
- ansatz siatki = suma dwoch par poziomych `theta_h_pair`
  (szew 0.124 rad niezaleznie od L); relaksacja 2000 krokow
  -> residuum ~1e-06,
- bramka G2: szachownica NIE przemieszcza rdzeni w ogole
  (0.0000) — konstrukcja skaluje sie z L,
- detektor 2+2, referencja blizniacza, wyzwolenie centroidem
  z oknem zamrazanym, pomiar wektorowy (Px, Py), eikonal
  z plaszczyzna centroidu 3. probki — wszystko bez modyfikacji
  poza nominalami geometrii,
- dyspersja sieciowa: analityka (omega_lat) przewiduje pomiar
  do 0.01 p.p. — mozna projektowac pasma G3 z gory,
- dip b_eff: b_nom - b_eff in [0.27, 0.72] tutaj (L=128:
  [0.32, 0.84]) — margines 0.9 w kryterium projektowym wystarcza.

**Pulapki na przyszlosc:**
- gamma_core runow z dluga kontynuacja (D16) potrafi przekroczyc
  prog 0.01 — projektowac progi metrologiczne osobno dla runow
  obserwacyjnych z kontynuacja,
- superpozycja: centroid delta_psi jest wolniejszy od v_g (ogon
  w oknie); limit krokow liczyc z zapasem ~10% nad (dystans/v_g),
- koszt N=512: caly cykl ~30-60 min na fazę pomiarowa; dluzsze
  runy w tle, interpretacja dopiero po komplecie wynikow.

**Kluczowe liczby (cytowane dalej):**
- substrat: s* = 1.174166, V''(s*) = 0.878665, m_TV = 0.9374,
  A_t = 2.5690
- k0/lambda/v_g: (1.1): 0.814/7.72/0.370; (1.3): 1.274/4.93/0.490;
  (1.5): 1.656/3.79/0.552; (1.7): 2.006/3.13/0.590
- G3 na L=256: bledy {-0.03, -0.12, -0.37, -0.93, -1.73}% dla
  k {0.540, 0.736, 1.006, 1.350, 1.669}
- bilans P1 programu: 19/19 waznych runow ze znakiem KU wirowi;
  11/12 punktow kryterialnych obu prob w pasmie [0.5, 2.0]
  (jedyny poza: rezim dyfrakcyjny b_eff ~ lambda, zastrzezony
  z gory w D12 proby #2)
