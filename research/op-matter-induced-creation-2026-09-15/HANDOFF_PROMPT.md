# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-matter-induced-creation-2026-09-15`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-matter-induced-creation-2026-09-15`. Pracujesz w vaulcie Obsidian;
repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu vaulta. ZAKAZ
zmiany cwd (`cd`); po KAŻDYM zapisie zweryfikuj `ls`, że plik wylądował
we właściwym miejscu (artefakt `TGP/TGP_v1/TGP/...`). Konwencje:
przeczytaj `TGP/TGP_v1/CLAUDE.md`. NIE czytaj STATE.md ani INDEX.md.

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-matter-induced-creation-2026-09-15/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/detektorów/listy λ̃/rampy/okien/sponge po
pierwszym biegu produkcyjnym; korekty tylko dla udokumentowanego błędu
implementacji — correction_note PRZED użyciem wyniku, pierwotne
outputy zachowane). Wykonaj WSZYSTKIE deliverables (FINAL + NEEDS +
dopis README). **KRYTYCZNE (dwukrotna lekcja poprzedników): NIE kończ
tury między fazami ani czekając na obliczenia w tle — czekaj AKTYWNIE
(skrypt czekający w foreground, pętla z time.sleep sprawdzająca plik
wynikowy). Tura może się skończyć dopiero, gdy wszystkie deliverables
istnieją na dysku.**

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14/`:
   `Phase_method_decisions.md` + `engine_core.py` (SKOPIUJ z cytatem;
   jedyna zmiana merytoryczna: λ̃(t) wg zamrożonej rampy z LOCKa §3),
   `Phase1_output.txt` (formy + 𝒰_mat + δψ_lin — CYTATY),
   `Phase3_qh1_output.txt` (kotwice regresji: λ̃=0.05 → δψ(0)=−0.134018;
   λ̃=0.5 → COLLAPSE t=0.375), `Phase_FINAL_close.md` (kontekst
   λ̃_crit∈(0.05,0.2], stan podprogowy ψ(0)→0.61).
2. `TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/NEEDS.md`
   (N2 — geneza wątku; detektor progowy 5/6).

## Kolejność (LOCK §2–§3)
0. `Phase_method_decisions.md` FROZEN: cytaty form; JAWNY zapis
   λ̃(t)=λ̃·S((t_off+Δ−t)/Δ) (smootherstep S, Δ=100, t_off=600,
   λ̃=0 od t=700); klasyfikatory SETTLED-SUB/SETTLED-DEF/COLLAPSE/
   UNSETTLED (faza włączona, okno [500,600]) i PERSISTENT-OBJECT/
   RETURN-TO-VACUUM/COLLAPSE/INCONCLUSIVE-RUN (po wygaszeniu, od
   t=700); reguła bisekcji λ̃_crit (6 kroków); okna pomiarowe.
1. **Phase 1** (`Phase1_analytic.py`+output): P1-I1 gate form 1e−12;
   P1-I2 krzywa ψ_min(λ̃) i λ̃_fold z minimum algebraicznego
   𝒰+λ̃ψ²/(4−3ψ) (sympy/numeryka 0D) — TABELA dla λ̃∈{0.06…0.18}
   zapisana PRZED Phase 3; P1-I3 cytat kontekstu (ψ≡1 jedyny stan
   trwały bez źródła).
2. **Phase 2** (`Phase2_gate.py`+output): P2a próżnia ‖ψ−1‖∞≤1e−10;
   P2b regresje (λ̃=0.05: ψ̄(0)=1−0.134018 ±1%; λ̃=0.5: COLLAPSE
   t=0.375±5%); P2c dryf ≤1e−6/100T₀ przy λ̃=const=0.05. FAIL ⟹ STOP.
3. **Phase 3 — Q-I1** (`Phase3_qi1_settle.py`+output): 7 biegów λ̃
   z listy + kotwice {0.05,0.20}, h=0.05, t_on=600; klasyfikacja
   w oknie [500,600]; potwierdzenie h=0.025 dla najgłębszego
   SETTLED-SUB (lub najgłębszego SETTLED-DEF, jeśli SUB brak);
   bisekcja λ̃_crit (6 kroków, FROZEN). Werdykt Q-I1 wg litery.
   Deskryptywnie: ψ̄(0) vs tabela P1-I2.
4. **Phase 3 — Q-I2** (`Phase3_qi2_rampoff.py`+output): dla KAŻDEGO
   SETTLED-*: rampa (checkpoint z t=600) + ewolucja swobodna do
   t=1700; klasyfikacja od t=700; potwierdzenia: PERSISTENT-OBJECT →
   h=0.025+dt/2; jeden RETURN-TO-VACUUM (najniższe λ̃ SETTLED) →
   h=0.025; kontrola czystości wygaszania λ̃=0.01 (RETURN oczekiwany,
   max|ψ−1|<1e−3). Werdykt Q-I2 wg litery. Konfrontacja z predykcją
   pre-rejestrowaną (RETURN/COLLAPSE oczekiwane) — bez reinterpretacji.
5. `Phase_FINAL_close.md` (wzorzec frontmattera: poprzednik
   `op-collapse-matter-source-2026-09-14/Phase_FINAL_close.md`),
   `NEEDS.md` (drzewo LOCK §5), dopis logu `README.md`
   (folder_status: closed, verdict: w frontmatterze) + regeneracja
   `python TGP/TGP_v1/tooling/build_cycles_index.py`.

## Wskazówki techniczne (nie zmieniają LOCKa)
- Budżet umiarkowany: faza włączona 9 biegów × 120k kroków × 4000 pkt;
  bisekcja ~6 krótkich; wygaszanie: do 9 biegów × 220k kroków.
  Batch w tle + checkpointy npz co 100 j.cz. + AKTYWNE czekanie.
- Energia: przy λ̃=const zachowana (gate P2c); podczas rampy NIE —
  raportuj bilans deskryptywnie (praca źródła), nie bramkuj.
- `integrity_snapshot.txt`: SHA256 LOCKa+MD po FROZEN, weryfikacja
  przy zamknięciu.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
  Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko
odczyt; formy bez modyfikacji; lista λ̃/rampa/progi/detektory/okna/
sponge niezmienialne po pierwszym biegu; ZAKAZ podłóg/barier poza
pasem; ZAKAZ pól dynamicznych i ρ(ψ) (S05); ZAKAZ claimów o masach/
oscylonach; predykcje P1-I2 i Q-I2 nienaruszalne; INCONCLUSIVE ≠
pozytyw; rejestr WEJŚĆ flagowany [INPUT].

## Raport końcowy
Werdykty Q-I1/Q-I2 z literą; tabela fazy włączonej (λ̃: klasa, ψ̄(0)
zmierzone vs P1-I2); λ̃_crit±okno (bisekcja); tabela wygaszania
(λ̃: klasa po rampie, E_core(1700)/E_core(700), max|ψ−1| końcowe);
konfrontacja z predykcją Q-I2; pliki; korekty/incydenty; higiena.
Pracuj do skutku.
