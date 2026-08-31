---
title: "NEEDS — user-gated po GATE-FAIL-STOP (P1c): pytanie ω²(n) w 3D pozostaje OTWARTE; wniosek metodologiczny + opcje nowego LOCKa mostu"
date: 2026-08-31
type: needs
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/NEEDS.md]]"
---

# NEEDS (wszystko poniżej wymaga decyzji użytkownika; zero automatycznych zmian)

Gałąź drzewa LOCKa §6: **„P1b/P1c FAIL → STOP (maszyneria; raport, bez
Phase 2–4)"**. Zgodnie z konstrukcją fail-fast cykl zakończył się na taniej
bramce, PRZED ciężkim rachunkiem — bez kosztu i bez fałszywych liczb.

## 1. Stan pytania Q (do odnotowania, nie do zmiany)

- **ω²_min(3D) NIE ISTNIEJE jako liczba** — hipoteza stabilizacji
  gęstością w 3D sc NIE została ani potwierdzona, ani obalona.
- Pre-rejestracja LOCKa §0 (negatyw 1D nieprzenośny: twierdzenie węzłowe
  Hilla nie działa w 3D, ogony 1/r) pozostaje w mocy: 3D to nadal jedyna
  niewykluczona geometria linii ω²(n).
- ZAKAZ cytowania tego cyklu jako „3D też negatywne" — to STOP maszynerii,
  nie wynik fizyczny.

## 2. Wniosek metodologiczny (udokumentowany w close, §2)

Mod runaway w dosłownym modelu #63 M0-f_ε wymaga rozdzielczości h≈0.015–0.1
(wąska sferyczna kieszeń Q₆₃ ≈ −15.5 na powłoce g≈0.75, gdzie pochodne
regularyzacji f_ε spike'ują). Wykonalne siatki 3D zalockowane w P1c (h≈0.4,
N=76–100 na wymiar) mają błąd dyskretyzacji 19% (radialnie) i generują mody
pasożytnicze −8.8 (99% masy na powłoce) — **kartezjańska maszyneria 3D
w tej klasie i przy tym budżecie jest NIEZDOLNA certyfikować mostu
radialny→3D**. Operator sam w sobie jest poprawny (RQ modu kotwicznego
2.1% przy N=100; konwencja wagi zweryfikowana T1–T4).

## 3. Opcje kontynuacji (KAŻDA wymaga nowego Phase 0 LOCK; user-gate)

a) **Re-lock mostu przy wykonalnej rozdzielczości:** ilościowa granica
   z addendum 2: artefakt kieszeni tłumiony dla h ≲ 0.15 (N ≳ 200/wymiar
   przy L=30, ~8e6 węzłów — ciężkie, wykonalne w tle); błąd weight-1
   3.1% przy h=0.1. Kandydat: P1c przy h≈0.10–0.15 (N=200–300), z gate
   ±5% i ewolucją t* na tej samej siatce (koszt RK4 ~godziny/bieg).
b) **Most w modelu KANONICZNYM (K=g⁴):** rachunek centralny (Phase 2–4)
   i tak jest w akcji kanonicznej (K wielomianowe — bez wąskiej powłoki ε);
   most można przedefiniować jako porównanie radialnego i kartezjańskiego
   widma operatora kanonicznego wokół TEGO SAMEGO profilu interpolowanego
   (kotwica radialna kanoniczna do policzenia w 1D, tanio). Wymaga decyzji,
   że kotwice #63 (−1.3896/3.62) przestają być bramką 3D (zmiana litery
   LOCKa — wyłącznie user-gate).
c) **Metody niesiatkowe / adaptowane:** bazy sferyczne sprzężone
   (rozwinięcie w harmoniki Y_lm × siatka radialna h=0.015 — łączy
   dokładność radialną z anizotropią sieci), metody spektralne,
   AMR wokół powłoki, wykorzystanie pełnej symetrii kubicznej O_h
   (redukcja 48×).
d) **Nic nie robić:** pozostawić ω²(n) w 3D jako jawnie otwarte
   (decyzja aksjomatyczna o znaku W pozostaje bez podpory rachunkowej —
   jak w NEEDS bloch-chain).

## 4. Artefakty gotowe do reużycia (bez ponownego pisania)

`Phase2_relax_lattice3d.py` i `Phase3_bloch3d.py` — kompletne,
skompilowane, NIEURUCHOMIONE (zero wyników); operator 3D, Bloch ±1,
newton_krylov z preconditionerem, kontrole P3b/P3c i identyfikacja
translacji zaimplementowane. Korekta eigsh (tol=0, pokrycie translacji):
`Phase_correction_note_eigsh.md`. Do ewentualnego nowego cyklu wystarczy
nowy LOCK + ewentualna zmiana bramki mostu.
