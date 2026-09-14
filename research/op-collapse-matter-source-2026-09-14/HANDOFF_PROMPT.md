# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-collapse-matter-source-2026-09-14`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-collapse-matter-source-2026-09-14`. Pracujesz w vaulcie Obsidian;
repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu vaulta. ZAKAZ
zmiany cwd (`cd`); po KAŻDYM zapisie zweryfikuj `ls`, że plik wylądował
we właściwym miejscu (znany artefakt `TGP/TGP_v1/TGP/...`). Konwencje
projektu: przeczytaj `TGP/TGP_v1/CLAUDE.md` (krótki). NIE czytaj
STATE.md ani INDEX.md.

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/detektorów/startów/listy λ̃/ρ̂/sponge po
pierwszym biegu produkcyjnym; korekty tylko dla udokumentowanego błędu
implementacji — correction_note PRZED użyciem wyniku, pierwotne outputy
zachowane). Wykonaj WSZYSTKIE deliverables (Phase_FINAL_close.md +
NEEDS.md + dopis README) — cykl bez nich NIE jest zakończony; czekaj
aktywnie na batche, nie kończ tury między fazami.

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/research/op-r3-stationary-states-2026-09-14/`:
   `Phase_method_decisions.md` + `engine_core.py` (SKOPIUJ do swojego
   katalogu z cytatem; jedyna zmiana merytoryczna: człon materii
   −∂𝒰_mat/∂ψ w RHS + 𝒰_mat w ewaluatorze energii),
   `Phase_correction_note_2_energy_eval.md` (tożsamość energii bez
   kancelacji — dziedziczysz), `Phase1_output.txt` (formy — CYTAT),
   `Phase3_output.txt` + `Phase3_results/` (baseline λ̃=0: kategorie
   i czasy zdarzeń 4 startów-reprezentantów; PROFILE STARTOWE odtwórz
   IDENTYCZNIE wg definicji z tamtego LOCKa/MD).
2. `TGP/TGP_v1/core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`:
   eq:L-mat-unified (~linia 229) + audit note A4 2026-05-01 + L01 formal
   definition ρ (~linie 239–265) + eq:vol-element-M911 — TYLKO ODCZYT
   (Grep po etykiecie, potem Read z offsetem).
3. `TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/NEEDS.md` (N2 —
   geneza pytania Q-H1; detektory progowe 5/6, 7/6).

## Kolejność (LOCK §2–§3)
0. `Phase_method_decisions.md` FROZEN: cytaty form; wyprowadzona forma
   𝒰_mat (po P1-H1 — zapisz placeholder „wg P1-H1" i uzupełnij WYNIKIEM
   przed pierwszym biegiem numerycznym); lista λ̃; ρ̂(r)=exp(−r²/18);
   klasyfikatory DEFORMATION/THRESHOLD-PULL/COLLAPSE (nadkategoria)/
   RADIATED/STABILIZED/INCONCLUSIVE-RUN; okna pomiarowe.
1. **Phase 1** (`Phase1_matter_analytic.py`+output): P1-H1 wyprowadzenie
   𝒰_mat=λ̃ρ̂ψ²/(4−3ψ) z literalnego √−g·(q/Φ₀)ψρ (sympy, gate
   tożsamości + pochodna ψ(8−3ψ)/(4−3ψ)², 1e−12); P1-H2 odpowiedź
   zlinearyzowana δψ=−5λ̃(G_Yuk∗ρ̂) — wyprowadzenie + kwadratura splotu
   (wzorzec do gate'u Q-H1); P1-H3 fakty brzegowe (granice 𝒰_mat przy
   ψ→4/3⁻ i ψ→0⁺, sympy limit).
2. **Phase 2** (`Phase2_gate.py`+output): P2a próżnia λ̃=0 ‖ψ−1‖∞≤1e−10;
   P2b regresja qR3 a=−0.20: COLLAPSE t=4.74±2%; P2c dryf energii ze
   źródłem λ̃=0.05 ≤1e−6/100T₀. FAIL ⟹ STOP.
3. **Phase 3 — Q-H1** (`Phase3_qh1_response.py`+output): ψ≡1 + źródło,
   λ̃∈{0.01,0.05,0.2,0.5}, t=300, h=0.05 (λ̃=0.01 i 0.5 też h=0.025);
   klasyfikacja DEFORMATION/THRESHOLD-PULL/COLLAPSE; gate liniowy przy
   λ̃=0.01 (≤5% vs splot na r≤40); deskryptywnie δψ(0) vs liniowa dla
   wszystkich λ̃.
4. **Phase 3 — Q-H2** (`Phase3_qh2_stabilize.py`+output): 4 starty
   (qR3 ±0.20, gauss a=−0.30 σ=3, gauss a=+0.15 σ=6) × λ̃∈{0.05,0.2,0.5},
   t_max=1000, h=0.05; potwierdzenia wg reguły LOCKa (zmiana kategorii
   vs baseline → h=0.025+dt/2; kontrola negatywu qR3+0.20 przy
   najniższym λ̃ → h=0.025). Werdykty Q-H1/Q-H2 wg litery LOCKa §4.
5. `Phase_FINAL_close.md` (wzorzec frontmattera: poprzednik
   `op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md`),
   `NEEDS.md` (drzewo LOCK §5), dopis logu w `README.md`
   (+ `folder_status: closed`, `verdict:` w frontmatterze README).

## Wskazówki techniczne (nie zmieniają LOCKa)
- Splot Yukawy sferycznie symetryczny: δψ(r) = −(5λ̃/r)∫₀^∞ r′ρ̂(r′)
  [e^{−|r−r′|}−e^{−(r+r′)}]/2 dr′ (kwadratura scipy; wyprowadź w P1-H2).
- Budżet mały: Q-H1 ~10 krótkich biegów, Q-H2 ~12+potwierdzenia
  (kolapsy kończą się w t≈2–20; STABILIZED biegnie pełne t=1000).
  Batch w tle, checkpointy npz co 100 j.cz.
- `integrity_snapshot.txt`: SHA256 LOCKa + MD po FROZEN, weryfikacja
  przy zamknięciu.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
  Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko odczyt;
formy M,𝒦,𝒰,𝒰_mat bez modyfikacji po P1-H1; lista λ̃/ρ̂/starty/progi/
detektory/sponge niezmienialne po pierwszym biegu; ZAKAZ podłóg/barier
poza pasem klasyfikacyjnym; ZAKAZ pól dynamicznych (S05 — ρ to źródło
statyczne); λ̃<0 poza zakresem; INCONCLUSIVE ≠ pozytyw; ZAKAZ claimów
o masach leptonów/oscylonach; rejestr WEJŚĆ flagowany [INPUT].

## Raport końcowy
Werdykty Q-H1/Q-H2 z literą; tabela Q-H1 (λ̃ × siatka: klasa, δψ(0)
zmierzone vs liniowe); tabela Q-H2 (start × λ̃: klasa, czas zdarzenia/τ,
vs baseline λ̃=0); konfrontacja z P1-H2/P1-H3 (znak δψ, asymetria
sufit/podłoga); pliki; korekty/incydenty; higiena. Pracuj do skutku.
