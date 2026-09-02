---
title: "NEEDS — op-metametric-boundary: pytania user-gated po Q1-POS / Q2-INCONCLUSIVE"
date: 2026-09-01
type: needs
tgp_owner: research/op-metametric-boundary-2026-09-01
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
---

# NEEDS (wszystkie pozycje wymagają decyzji usera; zero samowolnych zmian korpusu)

Drzewo decyzyjne LOCKa §6 nie ma jawnej gałęzi dla Q2-INCONCLUSIVE —
poniższe pozycje raportują potrzeby wynikające z litery wyników, bez
przesądzania.

> **LOG WYKONANIA (2026-09-02, user-gate „ok, działaj"):**
> - **N1 ✅ ROZSTRZYGNIĘTE (dokumentacyjnie):** korpus MA naturalne górne
>   domknięcie — biegun metryki efektywnej M9.1'' (sek08a,
>   eq:vol-element-M911: √−g_eff=c₀ψ/(4−3ψ); ds²: g_tt→0 przy ψ=4/3)
>   ⟹ g_ceil=√(4/3)=1.15470. Obserwacja: tła bloch g_max=1.141–1.143
>   tuż pod granicą. Domknięcie METRYCZNE — zbieżne z hipotezą
>   „granicy metryki" z op-blocked-soliton-bang.
> - **N2 ✅ ADRESOWANE:** pudło L=4π (k_min=0.5<1) w cyklu-następcy.
> - **N1+N2 → realizowane jako
>   [[../op-metric-closure-relaxation-2026-09-02/Phase0_balance.md]]**
>   (obustronne domknięcie + większe pudło + dziedziczony detektor
>   nukleacji + nowy detektor obiektów górnych).
> - N3 (dopisek Q1-POS), N4 (geneza Γ+s_i), N5 (kanał g→+∞
>   w dynamice hamiltonowskiej) — OPEN, bez zmian.

## N1 — Kanał ucieczki g→+∞: czy korpus ma górne domknięcie? (kluczowe)

18/18 biegów relaksacyjnych sektora tachionowego załamuje się przez
skończono-czasową ucieczkę g→+∞ (U=g⁷/7−g⁸/8 nieograniczone z dołu),
ZANIM podłoga dolna (QB-2) w ogóle się aktywuje. Pytanie do korpusu:
czy istnieje NIE-ad-hoc mechanizm górnego ograniczenia Φ (analog górnej
spinodali; entropia substratu S_Γ = k_B N_B γ(φ−lnφ−1)
z cor:entropy-potential rośnie z φ, ale przy T_Γ≪1 nie domyka U)?
Bez tego każdy przyszły test „relaksacji do granicy metametrycznej"
w tym modelu skończy się tym samym kanałem.

## N2 — Pudło nadkrytyczne dla struktury przestrzennej

W zalockowanym pudle L=2π mody k≥1 od próżni są tłumione (jedyny
rosnący stopień swobody: k=0). Start genezowy nie mógł więc z zasady
wytworzyć struktury/nukleacji przestrzennej. Powtórka Q2 wymaga L>2π
(np. L∈{4π,8π}) — zmiana zalockowanej siatki = nowy cykl (user-gate).

## N3 — Q1-POS: dopisek interpretacyjny (user-gate)

μ(soliton μ|próżnia) = −0.179 < 0 i ε(2π) = −5.2e−4 < 0: próżnia
sektora tachionowego jest już PO stronie opłacalnej kreacji. Kandydat
na dopisek do materiałów programu „granica metametryczna": separacja
μ=0 leży między stanem pustym a konfiguracjami metrycznymi, nie „za"
próżnią; dotychczasowe Q-FAIL-e stabilności mierzone względem próżni
mierzyły stan, który sam jest źródłem zysku kreacyjnego.

## N4 — Program genezy Γ+s_i (dziedziczone z op-blocked-soliton-bang)

Gałąź drzewa dla Q2-FAIL (a fortiori dla INCONCLUSIVE): granica
metametryczna może wymagać poziomu 0 wprost (goły substrat relacyjny,
POST_CLOSE_BRAINSTORM_UPDATE op-blocked-soliton-bang §3-4) — osobny
cykl `op-bare-substrate-genesis` (propozycja tam sformułowana),
zamiast kolejnych prób w kontinuum kanonicznym.

## N5 — Charakter kanału g→+∞ a dynamika hamiltonowska

Gradient flow (przetłumiony) wybiera g→+∞; statyka newton_krylov
poprzedników znajdowała ucieczki g→0; ewolucje hamiltonowskie miały
detektor obustronny (min≤0.01 LUB max≥2g₀). Pytanie porządkujące
(deskryptywne, tani rachunek): który kanał dominuje w pełnej dynamice
z tłumieniem skończonym (nie ∞)? — user-gate, bo wymaga nowego locka.
