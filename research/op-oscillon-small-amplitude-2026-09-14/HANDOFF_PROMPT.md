# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-oscillon-small-amplitude-2026-09-14`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-oscillon-small-amplitude-2026-09-14`. Pracujesz w vaulcie Obsidian;
repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu vaulta. ZAKAZ
zmiany cwd (`cd`); po KAŻDYM zapisie zweryfikuj `ls`, że plik wylądował
we właściwym miejscu (znany artefakt `TGP/TGP_v1/TGP/...`). Konwencje
projektu: przeczytaj `TGP/TGP_v1/CLAUDE.md` (krótki). NIE czytaj
STATE.md ani INDEX.md.

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-oscillon-small-amplitude-2026-09-14/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/detektora/rodzin startów/triage/sponge po
pierwszym biegu produkcyjnym; korekty tylko dla udokumentowanego błędu
implementacji — correction_note PRZED użyciem wyniku, pierwotne outputy
zachowane). Wykonaj WSZYSTKIE deliverables (Phase_FINAL_close.md +
NEEDS.md + dopis README) — cykl bez nich NIE jest zakończony; czekaj
aktywnie na batche, nie kończ tury między fazami.

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/research/op-r3-stationary-states-2026-09-14/`:
   `Phase_method_decisions.md` (schemat integratora, detektor, sponge —
   DZIEDZICZYSZ z cytatem), `engine_core.py` (SKOPIUJ do swojego
   katalogu z cytatem w method_decisions; zmiany tylko parametryczne
   R/t_max/sponge-okno + staging), `Phase1_output.txt` (formy M,𝒦,𝒰 +
   ω₂=−139/24 — CYTAT, nie wyprowadzaj ponownie),
   `Phase_correction_note_2_energy_eval.md` (tożsamość energii bez
   kancelacji — dziedziczysz), `Phase3_output.txt` (baseline τ=206.9
   dla regresji P2b).
2. `TGP/TGP_v1/research/op-r3-stationary-states-2026-09-14/NEEDS.md`
   (N1 — geneza tego cyklu).

## Kolejność (LOCK §2–§4)
0. `Phase_method_decisions.md` FROZEN: cytaty form i wartości ω₂;
   parametry R=400, t_max=10000, sponge [320,400]; staging (Etap A
   t=2000 → triage E_core≥0.2·E_ref → Etap B t=10000); detektor
   OSCILLON (3 warunki + potwierdzenia), kategorie OSCILLON-WEAK/
   RADIATED/COLLAPSE(nadkategoria)/INCONCLUSIVE-RUN; tabela
   przewidywanych ω(a) (miękka) ZAPISANA PRZED Phase 3.
1. **Phase 1** (`Phase1_analytic.py`+output): P1a′ gate cytatów form
   sympy vs float 1e−12 w ψ∈{0.9,1,1.1}; P1b′ wpis ω₂ (cytat) + tabela
   ω(a)≈1+ω₂a² dla a∈{0.02,0.05,0.08,0.10}.
2. **Phase 2** (`Phase2_gate.py`+output): P2a próżnia ‖ψ−1‖∞≤1e−10;
   P2b regresja (a=+0.15,σ=3): τ=206.9 ±5%; P2c dryf ≤1e−6/100T₀
   + odbicie sponge ≤1e−3 (różnicowo vs R=800, t=350). FAIL ⟹ STOP.
3. **Phase 3 — Q-G** (`Phase3_evolve.py`+output+json/npz per bieg):
   Etap A: 12 startów gauss (a∈{0.02,0.05,0.08,0.10}×σ∈{3,6,10}) + vac,
   h=0.05, do t=2000; triage FROZEN; Etap B: żywe do t_max=10000
   (checkpointy npz co 500 j.cz., wznawiaj z checkpointów). FFT ψ(0,t)
   (okno Hanna, segment stabilny ≥1000 j.cz.) dla żywych. Potwierdzenia:
   kandydaci OSCILLON/OSCILLON-WEAK → h=0.025 + dt/2; obowiązkowo
   (a=0.05,σ=6) → h=0.025 niezależnie od klasy; COLLAPSE → dt/2.
   Werdykt Q-G wg litery LOCKa §5. Obowiązkowo deskryptywnie: tabela
   ω_peak vs przewidywane ω(a) (bez progu PASS/FAIL — miękka).
4. `Phase_FINAL_close.md` (wzorzec frontmattera: poprzednik
   `op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md`),
   `NEEDS.md` (drzewo LOCK §6), dopis logu w `README.md`
   (+ `folder_status: closed`, `verdict:` w frontmatterze README).

## Wskazówki techniczne (nie zmieniają LOCKa)
- Budżet: Etap A ≈ 13 biegów × 400k kroków × 8000 pkt — batchuj w tle,
  etapy ≤50 min, checkpointy npz + log postępu; zapisuj ψ(0,t)
  i E_core(t) co dt_out=0.1; pełny profil co 100 j.cz.
- Δω FFT: okno ≥1000 j.cz. ⟹ Δω≈2π/1000·(bin) — raportuj Δω obok
  ω_peak; interpolacja paraboliczna piku dozwolona (zapisz w MD).
- CFL: przy małych amplitudach c≈1, dt=0.005 bezpieczne; monitoruj
  max|ψ−1| i lokalne c=(4−3ψ)/ψ przy każdym checkpoincie.
- `integrity_snapshot.txt`: SHA256 LOCKa + MD po FROZEN, weryfikacja
  przy zamknięciu.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
  Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko odczyt;
formy M,𝒦,𝒰 bez modyfikacji (cytaty); ZAKAZ podłóg/barier; starty/
detektor/progi/triage/sponge niezmienialne po pierwszym biegu; mapa
ω(a) pozostaje miękka; INCONCLUSIVE ≠ pozytyw; ZAKAZ claimów o masach
leptonów i dyskretności rodzin; rejestr WEJŚĆ flagowany [INPUT].

## Raport końcowy
Werdykt Q-G z literą; tabela biegów (start × siatka: klasa, τ/czas
zdarzenia, ω_peak±Δω, E_end/E_ref); konfrontacja ω_peak vs mapa LP
(deskryptywnie); los triage (kto odpadł w Etapie A); pliki; korekty/
incydenty; higiena. Pracuj do skutku.
