---
title: "Phase_correction_note — korekta macierzy kontroli P1c: start sieciowy 2π leży częściowo POZA domeną metryki PRIMARY (g_max=1.4734 > g_ceil, 6.6% objętości), więc jego włączenie do KONTROLI PRIMARY naruszało zamrożoną zasadę domenową MD §6; pierwotny output zachowany"
date: 2026-09-02
type: correction-note
tgp_owner: research/op-metric-closure-relaxation-2026-09-02
status: RECORDED-BEFORE-USE
related:
  - "[[Phase_method_decisions.md]]"
  - "[[Phase1_output_pre_correction.txt]]"
  - "[[Phase0_balance.md]]"
---

# Correction note (zapisana PRZED użyciem skorygowanego wyniku P1c)

## Błąd (udokumentowany)

Zamrożona macierz P1c (MD §6) włączyła bieg `lat × PRIMARY` do KONTROLI
sektora stabilnego, mimo że ta sama sekcja MD zamroziła zasadę: kontrola
wymaga startu w DOMENIE WAŻNOŚCI modelu (na tej podstawie z kontroli
PRIMARY wyłączono soliton, którego rdzeń g₀=1.465 > g_ceil).

Przesłanka faktograficzna włączenia `lat × PRIMARY` była FAŁSZYWA:
obserwacja LOCKa §0 („tła bloch g_max = 1.1406/1.1427/1.1429 — tuż POD
g_ceil") dotyczy tła ŁAŃCUCHA 1D z cyklu bloch, a nie tła 3D. Weryfikacja
w danych wejściowych (READ-ONLY, mtime niezmieniony):

- `2pi__A1.0__N32`: g ∈ [0.6064, **1.4734**], frakcja g > g_ceil = 6.59%
- `2pi__A1.0__N48`: g ∈ [0.6100, **1.4688**], frakcja g > g_ceil = 6.59%

Start sieciowy leży zatem częściowo POZA domeną metryki PRIMARY
(ψ > 4/3 — g_tt zmienia znak, metryka M9.1'' nie istnieje) — dokładnie
tak samo jak soliton. Włączenie go do kontroli PRIMARY było błędem
implementacji zamrożonej zasady MD §6 (niesprawdzenie zakresu g tła 3D),
nie nową decyzją merytoryczną.

## Pierwotny wynik (zachowany: `Phase1_output_pre_correction.txt`)

`p1c_lat_PRIM`: BREAKDOWN t=0.04 (start w strefie regularizacji bieguna:
‖ġ‖∞(t=0) = 1.07e8 — sztywność sił W̃_clip≈1.9e5 przy dziedziczonym
dt=0.01). Pozostałe 4 biegi kontroli: PASS (lat CBAR, gen PRIM,
gen CBAR, sol CBAR — wszystkie STATIONARY w próżni, ZERO alarmów obu
detektorów, w tym przy zasianych obiektach N_seed_dn=1, N_seed_up=1).

## Korekta (jedyna zmiana)

Macierz KONTROLI P1c zredukowana zgodnie z literą zamrożonej zasady
domenowej MD §6: `lat × PRIMARY` wyłączony z kontroli (jak sol × PRIMARY).
Kontrola PRIMARY na starcie w domenie: `gen × PRIMARY` (PASS).
P1c po korekcie: 4 biegi.

BEZ ZMIAN pozostają: kryteria/progi/detektory/seed/siatki; macierz P2
(LOCK §3 nakazuje `lat × PRIMARY` i `sol × PRIMARY` w sektorze
TACHIONOWYM — biegną wg litery, ich zachowanie przy starcie poza domeną
jest WYNIKIEM Phase 2, raportowanym wprost); zachowanie
`lat × PRIMARY` w sektorze stabilnym raportowane deskryptywnie
w Phase_FINAL (BREAKDOWN t=0.04, pre-correction output).

Ważność detektorów (cel P1c wg LOCKa: „ZERO alarmów OBU detektorów")
jest ustalona przez 4 biegi PASS — detektory zdolne do FAIL i czyste.
