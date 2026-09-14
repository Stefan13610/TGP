# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-r3-stationary-states-2026-09-14`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-r3-stationary-states-2026-09-14`. Pracujesz w vaulcie Obsidian;
repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu vaulta.
ZAKAZ zmiany cwd (`cd`) — pełne ścieżki; po KAŻDYM zapisie pliku
zweryfikuj `ls`, że zmaterializował się we właściwym miejscu (znany
artefakt zagnieżdżonych `TGP/TGP_v1/TGP/...`).

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-r3-stationary-states-2026-09-14/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/detektora/rodzin startów/sponge po
starcie; korekty wyłącznie dla udokumentowanego błędu implementacji —
correction_note PRZED użyciem wyniku, pierwotne outputy zachowane).
UWAGA: wykonaj WSZYSTKIE deliverables do końca (Phase_FINAL_close.md
+ NEEDS.md + dopis README) — cykl bez nich NIE jest zakończony;
czekaj aktywnie na batche, nie kończ tury między fazami.

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/research/op-action-audit-spectrum-insert-2026-09-13/`:
   `Phase1_output.txt` (jawne M(ψ)=ψ⁶/(4−3ψ)², 𝒦=ψ⁴, 𝒰=γ(ψ⁴/4−ψ³/3),
   π=Mψ̇, EOM — DZIEDZICZYSZ formy z cytatem w method_decisions),
   `Phase2_output.txt` (ω²=k²+1 — do porównania w P2b),
   `Phase_FINAL_close.md` (kontekst decyzji |g^tt|).
2. `TGP/TGP_v1/core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`:
   rem:W-sign-axiomatic(iv) + rem:psi-EOM-R3-branch-status (dopiski
   2026-09-14 — kontekst konwencji; TYLKO ODCZYT).
3. `TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/Phase3_relax_M911.py`
   (siatka radialna, energia 4π∫r²dr — do adaptacji; UWAGA: tam był
   gradient flow, tu dynamika 2. rzędu — silnik piszesz od nowa).

## Kolejność (LOCK §2)
0. `Phase_method_decisions.md` FROZEN: cytaty form; schemat
   integratora 2. rzędu na parze (ψ, π=Mψ̇) (leapfrog/Verlet z jawną
   obsługą M(ψ): scheme FROZEN z uzasadnieniem symetrii); forma
   sponge γ_sp(r) gładka na r∈[160,200]; detektor oscylonu; E_core
   (r≤80); klasyfikacja BREAKDOWN-BOUNDARY (ψ>4/3−1e−6 lub ψ<1e−6).
1. **Phase 1** (`Phase1_stationary.py`+output): P1a κ²=ω²−1
   (wyprowadzenie sympy, klasy ω<1/ω>1, mapowanie na linearyzację R3
   κ=1⟺ω²=2); P1b znak przesunięcia częstości Lindstedt–Poincaré
   O(a²) z 𝒰‴(1),𝒰⁗(1),M′(1) — PREDYKCJA pre-rejestrowana; P1c gate
   sympy vs float 1e−12 w {0.9, 1, 1.1}.
2. **Phase 2** (`Phase2_gate_dynamics.py`+output): P2a próżnia 100 T₀
   ‖ψ−1‖∞≤1e−10; P2b test dyspersji (puls a=1e−3, FFT (r,t), ≥3 mody,
   zgodność z k²+1 ≤1%); P2c dryf energii ≤1e−6/100 T₀ + odbicie
   sponge ≤1e−3. FAIL ⟹ STOP.
3. **Phase 3 — Q-E** (`Phase3_evolve.py`+output+npz/json per bieg):
   10 startów × 2 siatki (h∈{0.05,0.025}), R=200, dt=0.005,
   t_max=1000; klasyfikacja OSCILLON/RADIATED/BREAKDOWN-BOUNDARY/
   INCONCLUSIVE wg detektora FROZEN (E_core≥0.5·E_core(50) przez
   ≥100·T₀ + ≥50 przejść ψ(0,t) przez 1; potwierdzenie: h i dt/2,
   czas życia ±10%); FFT ω dla kandydatów. Werdykt Q-E wg litery.
   Deskryptywnie: los startów quasi-R3.
4. **Phase 4 — Q-F** (tylko przy Q-E-PASS; `Phase4_families.py`
   +output): (n, ω, czas życia) per oscylon; rodziny; werdykt
   Q-F-PASS/PARTIAL/FAIL wg litery.
5. `Phase_FINAL_close.md` (wzorzec frontmattera:
   `TGP/TGP_v1/research/op-action-audit-spectrum-insert-2026-09-13/Phase_FINAL_close.md`),
   `NEEDS.md` (drzewo §5), dopis logu w `README.md`.

## Wskazówki techniczne (nie zmieniają LOCKa)
- Ewolucje są dłuższe niż u poprzedników (t_max=1000, N~8000 przy
  h=0.025): batch w tle z checkpointami npz co 100 j.cz. i logiem;
  dziel na etapy ≤50 min; zapisuj ψ(0,t) i E_core(t) co dt_out=0.1
  (do FFT), pełny profil co 50 j.cz.
- Laplasjan radialny: (1/r²)(r²𝒦ψ′)′ przez różnice centralne
  z regularnością w r=0 (ψ′(0)=0, l'Hôpital); π-formulacja: ψ̇=π/M,
  π̇=RHS−(wkład ½M′ψ̇² przez π²M′/2M²) — zapisz JAWNIE w MD.
- FFT: scipy.fft, okno Hanna na stabilnym odcinku.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
- Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko
odczyt; formy M,𝒦,𝒰 bez modyfikacji (cytaty); ZAKAZ podłóg/barier;
sponge/detektor/progi/rodziny startów niezmienialne po pierwszym
biegu; INCONCLUSIVE ≠ pozytyw; ZAKAZ claimów o masach leptonów;
P1b pre-rejestracja nienaruszalna; rejestr WEJŚĆ flagowany.

## Raport końcowy
Werdykty Q-E/Q-F z literą; predykcja P1b vs wynik Phase 3; tabela
biegów (start × siatka: klasyfikacja, czas życia, ω, węzły);
los quasi-R3; pliki; korekty/incydenty; higiena. Pracuj do skutku.
