---
title: "Phase_FINAL_close — zamknięcie: Q-NET-A-PASS-MOSAIC (10/10 genez lockuje w mozaikę OBU znaków, metric 0.80–0.94 w τ=100) + Q-NET-B-INCONCLUSIVE wg litery (sieci grubieją przez skalę pudła przed τ=1000 — 4/5 okien nieważnych; jedyne ważne: slope −0.249 spłaszczony plateau pasków) + Q-NET-C-PASS-TOPOLOGICAL: geneza s=20260908 (A=1.0) kończy w TRWAŁYCH paskach nawiniętych na torus (L_wall=128.0 = 2 proste ściany, tails=1.0000; kontrola N=256 IDENTYCZNA: +0.516/−0.414, Lw=128.0) — trwała sieć przedmetryczna z LOSOWEJ genezy potwierdzona; łącznie 3/10 genez → stan nawinięty (~1/3, zgodne z literaturą quenchów Z2), 7/10 → jedna domena"
date: 2026-09-02
type: phase-final-close
tgp_owner: research/op-z2-wall-network-2026-09-02
status: CLOSED
verdict: "P1 PASS (po korekcie 1 komparatora: pola ifft2 rejestrowane węzłowo — współpołożenie DOKŁADNE 0.0000; L_wall paski=128.0 exact, kropla=64.0=8R (A1), ciągłość τ_life ✓). Q-NET-A: PASS-MOSAIC — wszystkie 10 genez (A∈{0.6,1.0}×5 seedów) lockuje do τ=100 w mozaikę obu znaków (A=1.0: plus/minus ≥0.107/≥0.342, metric 0.869–0.922, L_wall 154–253); ani FAIL-ONESIGN, ani FAIL-BARE nie wystąpił. Q-NET-B: INCONCLUSIVE wg litery — 4/5 sieci (A=1.0) znika przed końcem okna [100,1000] (grubienie przechodzi przez skalę pudła L=64; okna nieważne 80–316/451 próbek), jedyne ważne okno (s=20260908) daje slope −0.249 (R²=0.898) spłaszczony przez plateau stanu paskowego — wyznaczenie wykładnika wymaga większego pudła (mechanizm raportowany, nie strojono okna). Q-NET-C: PASS-TOPOLOGICAL — s=20260908 A=1.0 kończy w WOUND-STRIPES: dwie proste ściany nawinięte na torus, L_wall=128.0, tail_Lw=1.0000, tail_bare=1.0000, Ndom=1/1; kontrola pinningu N=256 (τ_max ten sam) ZGODNA i niemal identyczna liczbowo (+0.516/−0.414, Lw=128.0, tails 1.0000). Deskryptywnie: 3/10 wszystkich genez (1/5 przy A=1.0, 2/5 przy A=0.6) kończy w stanie nawiniętym (~1/3 — zgodne ze statystyką quenchów Z2 w 2D), 7/10 grubieje do jednej domeny. WIDOK POZIOMU 1: stan nawinięty niesie trwałą frakcję przedmetryczną 7.81% objętości (pasy ψ<0.218) przy „próżni" ψ≈1 w domenach; skala mozaiki 29.5 j.; etykieta znaku (±s*) — czyli CO oddzielają ściany — jest dla ψ niewidoczna."
anti_lakatos_lock: PRESERVED
tags: [z2-wall-network, genesis-quench, wound-stripes, topological-premetric, coarsening, level0, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_correction_note_p1c_colocation.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-premetric-pocket-level0-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-z2-wall-network

**Status: CLOSED-EXECUTED (jedna sesja: LOCK+A1 pre-code → Phase 1
(+korekta 1 komparatora) → Phase 2 (10 biegów + kontrola N=256)
→ zamknięcie).** Kryteria stosowane DOSŁOWNIE; zero zmian po starcie
obliczeń; jedna korekta implementacji TESTU (udokumentowana przed
użyciem).

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| P1 (estymatory+ciągłość) | **PASS** (po korekcie 1) | L_wall exact (paski 128.0, kropla 64.0=8R wg A1); GRF współpołożony 0.0000; τ_life ciągłe; H monotone |
| **Q-NET-A** (mozaika) | **PASS-MOSAIC** | 10/10 genez lockuje w mozaikę OBU znaków (metric 0.80–0.94, L_wall 119–286 w τ=100) |
| **Q-NET-B** (grubienie) | **INCONCLUSIVE** (litera) | 4/5 sieci grubieje przez skalę pudła przed τ=1000 — wykładnik wymaga większego L; ważne okno: −0.249 (plateau pasków) |
| **Q-NET-C** (trwałość topologiczna) | **PASS-TOPOLOGICAL** | s=20260908: paski nawinięte na torus TRWAŁE (tails=1.0000), kontrola N=256 niemal identyczna (Lw=128.0, +0.516/−0.414) |

## 1. Model i wejścia

Model substratu VERBATIM (op-bare-substrate-genesis; łańcuch cytowań
w LOCKu §1). Start: GRF pasmowy grid-niezależny (|n|≤16, obwiednia
exp(−(|n|/8)²), std z budowy N=256); A∈{0.6,1.0}; seedy
{20260906–20260910}; τ_max=2000 (100k kroków); estymator L_wall
manhattański (Amendment A1 pre-code); okno fitu τ∈[100,1000] FROZEN.

## 2. Wyniki

- **Mozaika (τ=100):** wszystkie geneza lockują z obiema fazami
  (przykład A=1.0: plus 0.107–0.549, minus 0.342–0.815, bare
  0.078–0.131 — sieć ścian + resztkowe jeziora gołe); żaden bieg nie
  został goły (kontrast informacyjny z G2: struktura nadbarierowa
  lockuje zawsze).
- **Grubienie:** sieci L_wall(τ) maleją i w 7/10 biegów znikają
  całkowicie (jedna domena) przed τ_max; w pudle L=64 okno [100,1000]
  jest za długie względem czasu życia sieci — wykładnik −1/2
  nieweryfikowalny w tej skali (raport wprost; okna NIE przesuwano).
  Jedyny bieg z pełnym oknem to ten, który wszedł w plateau paskowe
  (slope −0.249, R²=0.898 — mieszanka grubienia i plateau).
- **Stan nawinięty (topologiczny):** s=20260908 A=1.0 → dwie proste,
  równoległe ściany zamknięte przez torus: L_wall = 128.0 (dokładnie
  2L), plus=0.516/minus=0.406, tails 1.0000/1.0000, Ndom=1/1;
  **kontrola N=256 (dx/2, dt/4, ten sam GRF ciągły): ZGODNA
  i liczbowo niemal identyczna** (+0.516/−0.414, Lw=128.0, tails
  1.0000) — nie pinning. Deskryptywnie: nawinięte kończą też
  s=20260907 i s=20260910 przy A=0.6 ⟹ 3/10 (~1/3, zgodnie ze
  znaną statystyką stanów paskowych w quenchach Z2 2D — obserwacja
  porównawcza, zero claimów).
- **Widok poziomu 1 (obowiązkowy):** stan nawinięty = trwałe pasy
  przedmetryczne (ψ<0.218) o łącznej frakcji 7.81% objętości,
  oddzielające domeny, w których ψ≈1 jest nieodróżnialne od próżni;
  skala mozaiki L²·metric/L_wall = 29.5 j. Etykieta znaku (to, ŻE
  domeny są różne) jest dla ψ niewidoczna — sieć jest strukturą
  substratu, nie pola efektywnego.

## 3. Odpowiedź N1 (wprost; klasa: 2D, kwench relaksacyjny, tor)

**Losowa geneza substratu generyczne produkuje mozaikę domen Z2
z siecią arkuszy przedmetrycznych; sieć jest przejściowa (grubieje),
ALE skończona frakcja genez (~1/3 w zbadanej klasie) kończy w stanie
topologicznie TRWAŁYM — nawinięte ściany, które nie znikają nigdy.**
Trwała wielkoskalowa struktura przedmetryczna z czystej genezy istnieje
i jest wybierana topologią przypadku, nie strojeniem. W pudle L=64
przeżywa forma minimalna (2 proste ściany); statystyka form i prawo
grubienia w większych pudłach — NEEDS.

## 4. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ Amendment A1 (czynnik stereologiczny 4/π estymatora) — PRE-CODE,
  zapisany w LOCKu §6 przed napisaniem kodu.
- ✓ **Korekta 1** (`Phase_correction_note_p1c_colocation.md`, PRZED
  ponownym biegiem bramki i przed Phase 2): komparator GRF uśredniał
  2×2 przy rejestrze węzłowym pól ifft2 — błąd TESTU; po korekcie
  współpołożenie DOKŁADNE (0.0000). Pierwotny FAIL zachowany
  (`Phase1_output_pre_correction.txt`). Progi nietknięte.
- ✓ Zero zmian modelu/okien/seedów po starcie; wszystkie 10 biegów
  raportowane; kontrola N=256 wykonana dla pozytywu; H_Γ nierosnące
  we wszystkich biegach; okno fitu NIE przesuwane mimo nieważności
  (INCONCLUSIVE wg litery).
- ✓ Katalogi innych cykli tylko odczyt; rdzeń .tex/STATE/git
  nietykane; τ ≠ czas fizyczny; zero claimów kosmologicznych.
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1.

## 5. Pliki cyklu

`Phase0_balance.md` (LOCK + A1) · `Phase1_gate.py` →
`Phase1_output.txt` (+ `Phase1_output_pre_correction.txt`) ·
`Phase_correction_note_p1c_colocation.md` · `Phase2_network.py` →
`Phase2_output.txt` + `Phase2_results/` (json + npz stanów s i ψ)
+ `Phase2_run.log` + `Phase2_ctrl.log` · `NEEDS.md` (user-gated) ·
`README.md`.

## 6. Mapowanie na drzewo LOCKa §5

Gałąź: **A-PASS ∧ C-PASS-TOPOLOGICAL** (B INCONCLUSIVE z mechanizmem
skali pudła) ⟹ „geneza produkuje trwałą sieć przedmetryczną (ukrytą
dla ψ): NEEDS PILNE: wersja 3D (arkusze), statystyka konfiguracji
nawiniętych, user-gate dopisek core" — [[NEEDS.md]].
