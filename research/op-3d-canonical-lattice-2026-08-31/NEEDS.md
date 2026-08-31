---
title: "NEEDS — user-gated po Q-FAIL: hipoteza stabilizacji gęstością zamknięta także w 3D sc (model kanoniczny); pierwszy pełny negatyw rachunkowy ω²_min(3D)=−1.674; konsekwencje dla znaku W i drabiny mas wymagają decyzji użytkownika"
date: 2026-08-31
type: needs
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/NEEDS.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/NEEDS.md]]"
---

# NEEDS (wszystko poniżej wymaga decyzji użytkownika; zero automatycznych zmian)

Gałąź drzewa LOCKa §6 wykonana: **„Q-FAIL → hipoteza stabilizacji
gęstością zamknięta także w 3D sc (model kanoniczny); decyzja
aksjomatyczna o znaku W bez podpory"**.

## 1. Stan pytania Q (fakt rachunkowy, do odnotowania)

- Po raz pierwszy w linii ω²(n) istnieje LICZBA 3D:
  **ω²_min(3D, sc, d=2π) = −1.674350** (zbieżnie, po odjęciu translacji,
  potwierdzona nieliniowo: ucieczka t_esc≈4.1 ≤ 2·t*_izo przy obu
  znakach, dt-konsystentnie, ε-konsystentnie).
- Trend vs 1D: −1.674 < −1.222 — **trzeci wymiar pogłębia niestabilność**
  (argument nieprzenośności negatywu 1D skonsumowany rachunkiem).
- Zakres negatywu: PRIMARY po tłach istniejących (jedyne: d=2π);
  strict: dla d∈{π, d*₁=3.0790, 3π, 4π} sieć sc kanoniczna NIE ISTNIEJE
  w klasie startów superpozycyjnych (negatyw istnienia, nie stabilności).
- ZAKRES MODELOWY: akcja kanoniczna (K=g⁴), geometria sc, jeden soliton
  μ na komórkę, starty superpozycyjne — poza tą klasą wynik nie orzeka.

## 2. Konsekwencje wymagające decyzji (user-gate)

a) **Znak W / sektor tachionowy:** ostatnia niewykluczona droga
   rachunkowa „stabilizacji gęstością" w geometriach zbadanych (1D
   łańcuch, 3D sc) jest zamknięta. Decyzja aksjomatyczna o znaku W
   pozostaje bez podpory rachunkowej — do rozstrzygnięcia na poziomie
   aksjomatyki (retrospektywa + ewentualny zapis w rdzeniu, WYŁĄCZNIE
   decyzją użytkownika).
b) **Drabina mas / sek08b:** jeżeli linia ω²(n) ma pozostać w programie,
   wymaga nowej hipotezy (inna geometria, inna klasa teł, inny mechanizm)
   — nie kosmetyki obecnej.

## 3. Opcje kontynuacji rachunkowej (KAŻDA wymaga nowego Phase 0 LOCK)

a) **Inne sieci: bcc/fcc** (drzewo LOCKa przy braku istnienia wskazywało
   bcc/fcc; tu istnienie było częściowe, ale geometria sc jest
   najrzadsza — gęstsze upakowania zmieniają i istnienie, i widmo).
b) **Niepewność jednostronna d=3π:** N=48 samodzielnie dał kandydata
   tła (‖R‖∞=6.7e−11, ‖g−1‖∞=0.586), N=32 nie zbiega (NO-CONV) —
   punkt formalnie niezbieżny siatkowo. Ewentualny re-lock z N∈{48,64}
   (LOCK pozwalał dodać 64, nie zamienić) mógłby rozstrzygnąć istnienie
   przy 3π; kierunek widma przy 2π sugeruje jednak ten sam znak.
c) **Klasy startów:** starty superpozycyjne (1.0×, 0.7×) nie wyczerpują
   przestrzeni teł; tła wielosolitonowe/nietrywialne topologicznie —
   osobna hipoteza, osobny LOCK.
d) **Nic nie robić:** przyjąć negatyw jako domknięcie linii ω²(n)
   w geometriach kryształowych i wrócić do aksjomatyki znaku W.

## 4. Higiena cytowania

Wolno cytować: „sieć sc solitonów kanonicznych μ istnieje przy d=2π
i jest tachionowa (ω²_min=−1.674, potwierdzone nieliniowo)". NIE wolno
cytować jako: „wszystkie sieci 3D są niestabilne" (zbadano sc, jedną
kalibrację, jedną klasę startów) ani jako negatyw istnienia przy d=2π.
