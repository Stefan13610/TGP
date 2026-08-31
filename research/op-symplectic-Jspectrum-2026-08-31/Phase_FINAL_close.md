---
title: "Phase_FINAL_close — zamknięcie: Q-FAIL (odczyt PRIMARY tol wg rulingu FROZEN; odczyt literalny tol zdegenerowany — flagowany): tła sektora tachionowego NIE są spektralnie stabilne w dynamice symplektycznej I rzędu — σ(JL̂) ma zbieżnie Re λ = +0.1396…+0.1434 ≫ tol dla WSZYSTKICH 4 teł i wszystkich k; sama próżnia u=1 jest niestabilna symplektycznie (analitycznie: max Re λ = γ/4); kontrole Gate-A/C1/C2/C3 czyste"
date: 2026-08-31
type: phase-final-close
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
status: CLOSED
verdict: "Phase 1 (sympy, 16/16 exact): P1a inwariantność widmowa zanurzeń II rzędu — ω²=k²−γ dla DOWOLNEJ wagi F (V′=F·S) oraz dla obu nazwanych konwencji (energia K=g⁴/U i źródło #63 f=1+4ln g/W′); P1b gradient flow rządzi się tym samym Hessianem L₊; P1c: L₊/L₋ wyprowadzone z E[u] (czynnik ½ Wirtingera; JL̂=[[0,½L₋],[−½L₊,0]], λ²=−ν/4), tożsamość L₋g_d=0 on-shell EXACT, Q₊ ≡ operator Phase 3 bloch, dyspersja próżni λ²(κ)=−¼κ²(κ²−γ) ⟹ PRÓŻNIA TACHIONOWA NIESTABILNA TAKŻE SYMPLEKTYCZNIE (max Re λ=γ/4 przy κ²=γ/2). Phase 2 (bramka): Gate A kotwica λ_min(3π)=−1.222191 odtworzona (|δ|=1.19e−7 ≤ 1e−4); C1 NLS kubiczny STABILNY: max Re λ=1.5e−7 ≤ tol (mimo ujemnego kierunku L₊; ‖L₋φ‖=fp-0); C2 NLS |u|⁶u NIESTABILNY wykryty: λ=2.9152/2.9076 zbieżnie; C3 próżnia vs analityka: maxerr 6.8e−4 ≤ 1e−3 (ratio N400/N800 = 4.00, rząd 2) — PASS 4/4. Phase 3 (centralna): max Re λ (N=800) ZBIEŻNE we wszystkich 4 tłach × 9 k (Δ ≤ 4.3e−5 przy progu ~1.4e−3): 3π +0.143396; 4π +0.139810; 6π 2-garb +0.143413; 6π 1-garb +0.139571 — WSZYSTKIE ≫ tol PRIMARY ≈ 1.2e−2, niemal płasko w k (mod zlokalizowany); cross-check produktowy λ²=−ν/4 zgodny do ≤3.4e−8; ‖L₋g‖∞ ≤ 3.6e−12 (fp-0), symetria czwórek ≤ 8.2e−12; artefakty modów symetrii zidentyfikowane osobno (λ ~ 3e−4…7e−4, halving z N — O(h)). WERDYKT: Q-FAIL (PRIMARY; wszystkie tła — bez MIXED). Odczyt literalny-nieograniczony tol (1.5–6.1, rośnie ∝N²) dawałby formalnie 'Q-PASS' — zdegenerowany (przy N=2048 sklasyfikowałby ZNANY niestabilny C2 jako stabilny: λ=2.91 < tol_lit=10.5); ruling FROZEN przed obliczeniami (Phase_method_decisions §4), oba odczyty raportowane. Phase 4: NIE uruchomiona (LOCK: TYLKO przy Q-PASS). Progi/kryteria/tła LOCKa nietknięte."
anti_lakatos_lock: PRESERVED
tags: [symplectic-dynamics, J-spectrum, NLS, tachyonic-sector, orbital-stability, negative-result, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase1_correction_note.md]]"
  - "[[Phase2_correction_note.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Phase FINAL — zamknięcie cyklu op-symplectic-Jspectrum

Obliczenia i zamknięcie 2026-08-31 (jedna sesja implementatora).
Kryteria LOCKa (Phase0_balance.md §0, §2–3, §6) stosowane DOSŁOWNIE;
zero zmian progów/tol/list teł po starcie. Jedyna niedomkniętość litery
LOCKa (zbiór, po którym biegnie max|Im λ| w definicji tol) rozstrzygnięta
rulingiem zamrożonym PRZED obliczeniami (Phase_method_decisions.md §4) —
oba odczyty raportowane w każdym punkcie.

---

## 0. Werdykt

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| Gate Phase 1 (P1c) | **PASS** | L₋g_d = 0 on-shell sympy-EXACT; L₊ ≡ operator bloch; czynnik ½ wyprowadzony |
| Gate Phase 2 | **PASS 4/4** | kotwica |δ|=1.19e−7; C1 stabilny przechodzi; C2 niestabilny wykryty; C3 ≤ 6.8e−4 |
| **Q** (σ(JL̂) ⊂ iℝ na tłach łańcucha?) | **Q-FAIL** (PRIMARY; wszystkie 4 tła) | Re λ = +0.1396…+0.1434 zbieżnie dla WSZYSTKICH teł i wszystkich k; ≫ tol ≈ 0.012 |

**Sens rachunkowy (deskryptywnie, bez interpretacji ontologicznej):**
przejście z dynamiki II rzędu do symplektycznej I rzędu NIE stabilizuje
sektora tachionowego w 1D. Mechanizm C1 (soliton NLS stabilny mimo
ujemnego kierunku L₊) tu NIE działa, bo niestabilność nie siedzi
w dyskretnym modzie ujemnym neutralizowanym więzami — już SAMA PRÓŻNIA
u=1 jest symplektycznie niestabilna (analitycznie λ² = −¼κ²(κ²−γ) > 0
dla 0<κ<√γ, max Re λ = γ/4 = 0.25). Tła łańcucha spowalniają wzrost
(~0.14 vs 0.25 próżni, czyli ~1.8×), ale go nie usuwają. Hipoteza
stabilizacji przez zmianę klasy dynamiki zamknięta negatywnie w klasie
zbadanej (1D, akcja kanoniczna) — drzewo LOCKa §6, gałąź Q-FAIL.

## 1. Phase 1 — formalizacja analityczna (Phase1_output.txt): PASS 16/16

- **P1a:** dyspersja próżni ω² = k² − γ wychodzi (i) dla OGÓLNEGO
  zanurzenia II rzędu V′=F·S z DOWOLNĄ wagą F (sympy, Function F —
  silniejsze niż żądane), (ii) dla zanurzenia energetycznego (K=g⁴,
  U′=K·S), (iii) dla źródłowego #63 (f=1+4ln g, W′=S). Relabeling
  W_source/V_energy widmowo inwariantny w całej klasie II rzędu.
- **P1b:** linearyzacja gradient flow = −L₊ (ten sam Hessian; sympy);
  problem II rzędu to forma ważona L₊φ=ω²Kφ, K>0 (Sylvester — te same
  znaki widma).
- **P1c (centralna):** z E[u] = ∫[½K(|u|)|u_x|² + U(|u|)]dx przy
  u = g₀+a+ib: E₂ = ½⟨a,L₊a⟩+½⟨b,L₋b⟩, Q₊ = U″−K′g″−½K″g′²
  (on-shell ≡ K(2βg−3γg²+2g′²/g²) — operator bloch), Q₋ = ½K′g′²/g+U′/g
  (on-shell ≡ (Kg′)′/g); Wirtinger: i u_t = δE/δu* ⟹ ∂_t(a,b) =
  (½L₋b, −½L₊a), λ² = −ν/4, ν = spec(L₋L₊). **Tożsamość L₋g_d = 0
  on-shell EXACT — PASS** (warunek kontynuacji cyklu).
- **Analityka próżni (do C3):** λ²(κ) = −¼κ²(κ²−γ); niestabilność
  modulacyjna próżni tachionowej w klasie symplektycznej z max
  Re λ = γ/4 przy κ² = γ/2.

## 2. Phase 2 — bramka maszynerii (Phase2_output.txt): PASS 4/4

| Kontrola | Wynik | Kryterium |
|---|---|---|
| Gate A (kotwica bloch) | λ_min(3π) = −1.222191, |δ|=1.19e−7 | ≤ 1e−4 **PASS** |
| C1 (NLS kubiczny, znany stabilny) | max Re λ = 1.50e−7 (N=1024) / 5.5e−7 (N=2048) | ≤ tol PRIMARY 1.17e−2 (i literalnego) **PASS** |
| C2 (NLS |u|⁶u, znany niestabilny) | λ = 2.9152 / 2.9076 (zbieżnie, |Δ|=7.6e−3 ≤ 2.9e−2) | Re λ > 0 **PASS**; separacja od C1: 5.3e6× |
| C3 (próżnia vs analityka P1) | maxerr = 6.80e−4 (N=400), 1.70e−4 (N=800), ratio 4.00 | ≤ 1e−3 **PASS** (rząd 2) |

- C1 przechodzi MIMO ujemnego kierunku L₊ solitonu — dokładnie mechanizm,
  o który pyta cykl; maszyneria umie odróżnić oba przypadki (FAIL-e
  osiągalne w obie strony, zrealizowane przez C2).
- Dyskretyzacja on-shell Q₋ = (K·Dg)′/g daje ‖L₋φ‖, ‖L₋g‖ = fp-0
  (3.6e−12) — mod fazowy dokładny; solitony bramek Newton-polish do
  residuum ≤ 1e−12 na pół-komórce parzystej (pinning translacji).
- **Empiryczna korroboracja rulingu tol (§4 method decisions):** przy
  N=2048 tol literalny-nieograniczony = 10.49 > λ_C2 = 2.91 — odczyt
  literalny sklasyfikowałby ZNANY niestabilny przypadek jako „stabilny";
  jest grid-rozbieżny (∝N²) i sprzeczny z wymogiem osiągalnego FAIL.

## 3. Phase 3 — RACHUNEK CENTRALNY (Phase3_output.txt + Phase3_results.json)

**Tabela max Re λ (wartość główna N=800; wszystkie punkty ZBIEŻNE:
Δ_siatka ≤ 4.3e−5 przy progu ~1.4e−3; tol PRIMARY ≈ 1.2e−2,
tol literalny 1.5–6.1):** [INPUT-ONTO] u₀ = g_d (moduł pola = g).

| tło | d | max Re λ (k=0) | max po k (argmax) | werdykt PRIMARY | (literalny) |
|---|---|---|---|---|---|
| 3π (1-garb) | 9.4248 | **+0.143396** | +0.143405 (k=0.208) | FAIL | (PASS) |
| 4π (1-garb) | 12.5664 | **+0.139810** | +0.139810 (k=0) | FAIL | (PASS) |
| 6π (2-garb) | 18.8496 | **+0.143413** | +0.143415 (k=0.125) | FAIL | (PASS) |
| 6π (1-garb, S3single) | 18.8496 | **+0.139571** | +0.139571 (k=0.167) | FAIL | (PASS) |

- **Zbieżność:** wszystkie 4 tła × 9 punktów k × 2 siatki — ZBIEŻNE
  (kryterium LOCKa dosłownie). Wzrost niemal niezależny od k (rozstęp
  ≤ 1e−5 wzdłuż pasma — mod zlokalizowany na garbie); pary teł:
  6π 2-garb ≈ 3π (powielona orbita), 6π 1-garb ≈ 4π (głębsza amplituda)
  — spójne strukturalnie z tabelą ω²_min cyklu bloch.
- **Cross-check produktowy (LOCK):** max Re λ z JL̂ vs z λ²=−ν/4,
  ν=spec(L₋L₊): |δ| ≤ 3.4e−8 (wszystkie tła, obie siatki).
- **Tożsamości dyskretne:** ‖L₋g‖∞ ≤ 3.6e−12 (fp-0, mod fazowy);
  ‖L₊g′‖∞ = 8e−6…1.4e−4 (O(h²), halving z N ✓); symetria czwórek
  ±λ ≤ 8.2e−12.
- **Mody symetrii (diagnostyka §5 method decisions):** artefakty
  rozszczepienia Jordana widoczne jako λ ≈ +3e−4…+7e−4 (halving z N —
  O(h)); leżą 3 rzędy wielkości POD sygnałem 0.14. Uwaga: heurystyka
  overlapu ≥0.9 dla tła 3π flagowała także sam mod niestabilny
  (składowa b silnie przekrywa g — cecha prawdziwego modu wzrostu,
  nie artefaktu); flaga jest wyłącznie deskryptywna, kryterium używa
  pełnego max Re λ, a autentyczność modu 0.1434 potwierdzają: zbieżność
  siatkowa, płaskość w k, cross-check produktowy i brak flagi na
  pozostałych 3 tłach.
- **Werdykt per tło: FAIL ×4 ⟹ Q-FAIL (bez MIXED).** Odczyt literalny
  tol: „PASS" ×4 ⟹ formalnie „Q-PASS" — zdegenerowany (patrz §2 punkt
  ostatni i ruling §4 method decisions zamrożony PRZED obliczeniami);
  oba odczyty raportowane, PRIMARY rządzi werdyktem nagłówkowym.

## 4. Phase 4 — NIE uruchomiona (zgodnie z LOCK)

LOCK §3: Phase 4 „warunkowa, TYLKO przy Q-PASS". Werdykt PRIMARY =
Q-FAIL ⟹ Phase 4 nie przysługuje. Rozbieżność odczytów tol NIE
uruchamia Phase 4 (ruling FROZEN: PRIMARY rządzi; uruchomienie ewolucji
poza literą LOCKa byłoby forbidden move). Incydentów P3↔P4 brak
(Phase 4 nie było).

## 5. Złożenie werdyktu (litera LOCKa §3 Phase 3)

- **Q-PASS** wymaga max Re λ ≤ tol zbieżnie dla WSZYSTKICH teł i k —
  nie zaszło w odczycie PRIMARY (zaszłoby tylko w zdegenerowanym
  odczycie literalnym tol, który przy N=2048 nie wykrywa znanego
  niestabilnego C2 — flagowane, nieużyte).
- **Q-FAIL** — „istnieje tło z Re λ > próg zbieżnie": zaszło dla
  WSZYSTKICH 4 teł (silniejsze niż egzystencjalne; brak MIXED).
- Kontrole C1/C2/C3 + Gate A czyste ⟹ warunek jakości spełniony.

**Obserwacja analityczna (nie-kryterialna; korroboruje wynik):**
niestabilność symplektyczna jest tu strukturalna jak w cyklu bloch:
L₋ ≥ 0 (stan podstawowy g_d bezwęzłowy, wartość 0), a L₊ ma kierunek
ujemny ⟹ ν_min = min spec(L₋L₊) < 0 generycznie, λ = √(−ν/4) > 0.
W NLS ratuje soliton struktura więzów (ładunek; kryterium VK) działająca
na DYSKRETNY mod ujemny nad stabilną próżnią — tu ujemna jest już
PRÓŻNIA (kontinuum 0<κ<√γ), więc nie ma czego neutralizować więzami.
Zmiana klasy dynamiki nie omija tachionu próżni w 1D.

## 6. Korekty, dodatki, incydenty (pełna lista)

1. **Phase1_correction_note.md** (PRZED użyciem wyniku; pierwotny output
   zachowany: `Phase1_output_run1_original.txt`): (a) brak `positive=True`
   pól sympy ⟹ nieredukowalne (g₀²)^{7/2}; (b) zbędny znak − przy
   konwencji euler_equations w P1b. Czyste błędy implementacyjne
   skryptu symbolicznego; zero zmian testowanych tożsamości.
2. **Phase2_correction_note.md** (PRZED użyciem wyniku; pierwotny output
   zachowany: `Phase2_output_run1_original.txt`): błędny wzór WYŁĄCZNIE
   deskryptywnej diagnostyki „symetria czwórek" (suma zamiast różnicy
   posortowanych widm); bramki nietknięte i identyczne w obu biegach.
3. **Ruling tol** (Phase_method_decisions.md §4, FROZEN przed
   obliczeniami): LOCK nie specyfikuje zbioru dla max|Im λ|; odczyt
   PRIMARY pasmowy (|λ| ≤ 12) vs literalny-nieograniczony (∝N²,
   zdegenerowany — empirycznie nie wykrywa C2 przy N=2048). Oba
   raportowane w każdym punkcie; w Phase 3 odczyty ROZJEŻDŻAJĄ werdykt
   (FAIL vs formalny PASS) — raportowane jawnie, PRIMARY rządzi.
4. Heurystyka overlapu modów symetrii dała fałszywy-dodatni flag na
   prawdziwym modzie niestabilnym tła 3π (deskryptywna; §3 wyżej).
5. Środowiskowy: sandbox blokuje `cd` poza vaultem — wszystkie wywołania
   z pełnymi ścieżkami (i tak wymagane higieną LOCKa §4 p.6).

## 7. Pliki cyklu

Obliczeniowe: `Phase1_analytic.py` + `Phase1_output.txt`
(+ `Phase1_output_run1_original.txt`); `Phase2_gate_nls.py`
+ `Phase2_output.txt` (+ `Phase2_output_run1_original.txt`);
`Phase3_Jspectrum.py` + `Phase3_output.txt` + `Phase3_results.json`.
Metodologiczne: `Phase_method_decisions.md`, `Phase1_correction_note.md`,
`Phase2_correction_note.md`. Zamykające: ten plik, `NEEDS.md`
(user-gated), `README.md` (zaktualizowany). Wejście zewnętrzne
(READ-ONLY): `../op-bloch-chain-stability-2026-08-31/
Phase2_backgrounds.npz`. Rdzeń `.tex`, STATE.md, git, katalogi innych
cykli — NIETKNIĘTE.

## 8. Mapowanie na drzewo decyzyjne LOCKa §6

Q-FAIL → NEEDS: „klasa symplektyczna nie ratuje 1D; decyzja o znaku W
wraca do autora bez tej podpory; pozostaje 3D" — dosłownie wg drzewa.
Materiał dodatkowy (deskryptywny, do decyzji autora, nie decyzja):
niestabilna jest już próżnia klasy symplektycznej (γ/4), więc
ewentualna kontynuacja tej linii wymagałaby teł niehomogenicznych
w 3D lub innej klasy dynamiki — poza zakresem zbadanym.
