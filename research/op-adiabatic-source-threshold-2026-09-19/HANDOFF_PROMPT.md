# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-adiabatic-source-threshold-2026-09-19`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-adiabatic-source-threshold-2026-09-19`. Pracujesz w vaulcie
Obsidian; repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu
vaulta. ZAKAZ zmiany cwd (`cd`); po KAŻDYM zapisie zweryfikuj `ls`,
że plik wylądował we właściwym miejscu (znany artefakt
`TGP/TGP_v1/TGP/...`). Konwencje: przeczytaj `TGP/TGP_v1/CLAUDE.md`.
NIE czytaj STATE.md ani INDEX.md.

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-adiabatic-source-threshold-2026-09-19/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów / progów / detektorów / listy λ̃ / siatki Δ_on /
rampy / okien / sponge / reguły potwierdzenia siatkowego po pierwszym
biegu produkcyjnym; korekty tylko dla udokumentowanego błędu
implementacji — correction_note PRZED użyciem wyniku, pierwotne
outputy zachowane). Wykonaj WSZYSTKIE deliverables (FINAL + NEEDS +
dopis README).

**KRYTYCZNE (trzykrotna lekcja poprzedników): NIE kończ tury między
fazami ani czekając na obliczenia w tle — czekaj AKTYWNIE (skrypt
czekający w foreground, pętla z time.sleep sprawdzająca plik
wynikowy). Tura może się skończyć dopiero, gdy wszystkie deliverables
istnieją na dysku.**

## Pytanie cyklu w jednym zdaniu
Dwa poprzednie cykle zmierzyły λ̃_crit = 0.1077 **wyłącznie dla
nagłego załączenia źródła w t=0**; fold statyczny leży 2.65× wyżej
(0.2858) — sprawdzasz, czy przy ADIABATYCZNYM załączaniu stany
osiadłe istnieją powyżej 0.1078, czy próg saturuje z Δ_on, i czy
najgłębszy taki stan przeżywa wygaszenie źródła.

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/research/op-matter-induced-creation-2026-09-15/`
   — **cykl bezpośrednio poprzedzający, główne źródło**:
   `Phase_method_decisions.md` (§6 klasyfikatory i definicje
   pomiarowe — dziedziczysz 1:1; §3 rampa w dół — dziedziczysz dla
   Q-J3), `engine_core.py` (**SKOPIUJ**; ma już λ̃(t) — jedyna zmiana
   merytoryczna to profil rampy W GÓRĘ wg LOCK §3),
   `Phase1_output.txt` (tabela ψ_min(λ̃) i λ̃_fold=0.285769734239
   — CYTAT + niezależna kontrola do 1e−9),
   `Phase3_qi1_output.txt` (λ̃_crit=0.107656±0.000156; kotwice
   regresyjne), `Phase3_qi2_output.txt` (progi wygaszania),
   `NEEDS.md` (**N1, N3, N4 — ten cykl je realizuje**).
2. `TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14/Phase1_output.txt`
   — wyprowadzenie 𝒰_mat (CYTAT, nie wyprowadzaj ponownie).

## Kolejność (LOCK §2–§3)
0. `Phase_method_decisions.md` FROZEN + `integrity_snapshot.txt`
   (SHA256 LOCKa + MD + engine_core po FROZEN). Zapisz JAWNIE:
   λ̃(t)=λ̃·S(t/Δ_on) dla rampy w górę (S smootherstep),
   t_on=Δ_on+600, okno klasyfikacji [t_on−100, t_on]; klasyfikatory
   dziedziczone; **regułę potwierdzenia siatkowego (a)+(b) i regułę
   rozjazdu — decyduje h=0.025, etykieta GRID-DIVERGENT**; profil
   zjazdu dla Q-J3.
1. **Phase 1** (`Phase1_analytic.py`+output): P1-J1 gate form 1e−12
   + ciągłość λ̃(t) i pochodnej na końcach rampy; P1-J2 niezależne
   przeliczenie λ̃_fold (gate zgodności 1e−9 z 0.285769734239)
   + tabela ψ_min(λ̃) dla listy z LOCK §3 **zapisana PRZED Phase 3**;
   P1-J3 kryterium adiabatyczności (T₀=2π; Δ_on w jednostkach T₀);
   P1-J4 cytaty faktów.
2. **Phase 2** (`Phase2_gate.py`+output): P2a próżnia z pełną rampą
   ≤1e−10; P2b regresje dziedziczone (ψ̄(0)=0.865982 ±1%;
   COLLAPSE t=0.375 ±5%); **P2c NOWY: λ̃=0.05 z rampą Δ_on=100 vs
   bez rampy, |Δψ̄(0)| ≤ 1e−3**; P2d dryf ≤1e−6/100T₀.
   **FAIL ⟹ STOP** (zapisz co zawiodło i zakończ cykl raportem).
3. **Phase 3 — Q-J1** (`Phase3_qj1_settle.py`+output): Δ_on=100,
   pełna lista λ̃ {0.10,0.12,0.15,0.18,0.21,0.24,0.27,0.30}
   + kotwica 0.05; klasyfikacja w oknie [700,800]; potwierdzenia
   siatkowe (a) najgłębszy SETTLED-* ORAZ (b) największe λ̃ ≤
   0.8·λ̃_crit(Δ_on=100). Werdykt Q-J1 wg litery LOCK §4.
   Deskryptywnie: ψ̄(0) zmierzone vs tabela P1-J2.
4. **Phase 3 — Q-J2** (`Phase3_qj2_bisect.py`+output): bisekcja
   λ̃_crit dla Δ_on ∈ {25, 100, 400}, 6 kroków każda; dla
   Δ_on∈{25,400} przedział startowy z 4 biegów sondujących
   {0.10,0.15,0.21,0.27}. Werdykt Q-J2 wg litery (I₁, I₂, próg 0.25).
5. **Phase 3 — Q-J3** (`Phase3_qj3_rampoff.py`+output) —
   **TYLKO jeśli Q-J1-PASS**; inaczej zapisz `Q-J3-NIEURUCHOMIONE
   (warunkowe)` i przejdź dalej. Najgłębszy ZBIEŻNY SETTLED-*
   z Δ_on=100: kontynuacja z checkpointu t=t_on, rampa w dół
   (Δ=100, t_off=t_on), ewolucja swobodna do t_on+1100; klasyfikacja
   od t_off+Δ; **kryterium RETURN-TO-VACUUM WYŁĄCZNIE amplitudowe
   (max|ψ−1|(r≤40)<1e−3 w oknie końcowym) — człon energetyczny
   raportuj deskryptywnie, NIE bramkuj** (realizuje N3 poprzednika);
   PERSISTENT-OBJECT → potwierdzenie h=0.025 + dt/2. Konfrontacja
   z predykcją pre-rejestrowaną — bez reinterpretacji.
6. `Phase_FINAL_close.md` (wzorzec frontmattera: poprzednik
   `op-matter-induced-creation-2026-09-15/Phase_FINAL_close.md`),
   `NEEDS.md` (drzewo LOCK §5), dopis logu `README.md`
   (folder_status: closed, verdict w frontmatterze) + regeneracja
   `python TGP/TGP_v1/tooling/build_cycles_index.py`.

## Wskazówki techniczne (nie zmieniają LOCKa)
- **Budżet:** Q-J1 9 biegów do t=700 + 2 potwierdzenia h/2;
  Q-J2 3×(4 sondujące + 6 bisekcji), przy czym Δ_on=400 ma
  t_on=1000 (najdroższe); Q-J3 do 3 biegów długich. Batch w tle
  + checkpointy npz co 100 j.cz. + AKTYWNE czekanie. Etapy ≤50 min.
- **Uwaga na COLLAPSE przy rampie:** przy nagłym załączeniu kolaps
  zachodził w t<1; przy rampie może zajść DOPIERO w okolicy
  t≈Δ_on (gdy λ̃ dochodzi do plateau). Detektor COLLAPSE ma
  działać przez CAŁY bieg, nie tylko na początku.
- Energia: przy λ̃=const zachowana (gate P2d); podczas rampy NIE
  (praca źródła) — raportuj bilans deskryptywnie, nie bramkuj.
- **Konsola Windows to cp1250** — ustaw `PYTHONIOENCODING=utf-8`
  albo pisz ASCII do stdout; outputy i tak idą do
  `Phase*_output.txt` (UTF-8), nie do konsoli.
- Python: `python` (numpy/scipy/sympy). Sandbox: bez `/dev/null`,
  bez heredoc, **bez zapisu poza vaultem** (żadnego `/tmp`)
  — wszystko przez pliki skryptów w katalogu cyklu.

## Forbidden moves (egzekwowane)
Rdzeń `.tex` / STATE.md / git NIETYKANE; katalogi innych cykli tylko
odczyt; formy bez modyfikacji; lista λ̃ / siatka Δ_on / rampa / progi
/ detektory / okna / sponge / reguła potwierdzenia i rozjazdu
niezmienialne po pierwszym biegu; ZAKAZ podłóg/barier poza pasem;
ZAKAZ pól dynamicznych i ρ(ψ) (S05); ZAKAZ claimów o masach
leptonów i oscylonach; **predykcje P1-J2, Q-J1 (PASS oczekiwany)
i Q-J3 (RETURN oczekiwany) nienaruszalne — Q-J1-FAIL NIE wolno
przedstawić jako częściowego potwierdzenia**; INCONCLUSIVE ≠
pozytyw; rejestr WEJŚĆ flagowany [INPUT].

## Raport końcowy
Werdykty Q-J1 / Q-J2 / Q-J3 z literą; tabela Δ_on=100 (λ̃: klasa,
ψ̄(0) zmierzone vs P1-J2); **tabela λ̃_crit(Δ_on) dla {25,100,400}
z I₁, I₂ i testem saturacji**; wynik obu potwierdzeń siatkowych
(a) i (b) z jawnym stwierdzeniem, czy zaszedł GRID-DIVERGENT;
tabela wygaszania (jeśli uruchomione); **konfrontacja z λ̃_fold
= 0.2858 i z λ̃_crit^nagłe = 0.1077**; konfrontacja z predykcjami;
pliki; korekty/incydenty; higiena.
Pracuj do skutku.
