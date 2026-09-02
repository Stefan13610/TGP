# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-metric-pair-M911-2026-09-02`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP `op-metric-pair-M911-2026-09-02`.
Pracujesz w vaulcie Obsidian; repo = TGP/TGP_v1. WSZYSTKIE ścieżki
RELATYWNIE od rootu vaulta. ZAKAZ zmiany cwd (`cd`) — pełne ścieżki;
po KAŻDYM zapisie pliku zweryfikuj `ls`, że zmaterializował się we
właściwym miejscu (znany artefakt zagnieżdżonych `TGP/TGP_v1/TGP/...`).

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/detektorów/seedów/form po starcie; korekty
wyłącznie dla udokumentowanego błędu implementacji — correction_note
PRZED użyciem wyniku, pierwotne outputy zachowane).

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`
   (linie ~170–200 i ~310–380): DOKŁADNE formy w(ψ)=ψ/(4−3ψ)
   (eq:vol-element-M911), V_M9.1''(ψ)=−γψ²(4−3ψ)²/12, K(ψ)=K_geo·ψ⁴,
   metryka ds². Zamroź formy w method_decisions Z CYTATEM; rozjazd
   odczytów (czy w mnoży człon kinetyczny) rozstrzygnij cytatem
   z eq:S-unified-psi PRZED startem, oba odczyty odnotuj.
2. `TGP/TGP_v1/research/op-metric-closure-relaxation-2026-09-02/`:
   `Phase_FINAL_close.md` (diagnoza niekompatybilnej hybrydy — kontekst),
   `Phase_method_decisions.md` + `Phase2_two_sided_relax.py` (silnik
   gradient flow i detektory do adaptacji; detektory przenosisz
   z progami W ψ wg LOCKa: dolny ψ<5/6, górny ψ>7/6).
3. `TGP/TGP_v1/research/op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz`
   — READ-ONLY (weryfikuj mtime po odczycie); start (iii): ψ=g²,
   przeskalowanie do ψ_max=1.30 wg procedury FROZEN.

## Kolejność (LOCK §2)
0. `Phase_method_decisions.md` FROZEN (formy z cytatami; wariacja
   δE/δψ z członami w′ i K′; schemat gradient flow; obsługa pasa
   granicznego ψ>4/3−1e−6 jako BREAKDOWN-BOUNDARY; detektory).
1. **Phase 1 — Q-A** (`Phase1_landscape.py`+output): sympy — punkty
   krytyczne i krzywizna ρ_eff=w·V na (0,4/3), zachowanie przy 0⁺/4/3⁻,
   ograniczoność E; gate P1b: sympy vs float 1e−12 w {0.5, 1, 7/6, 1.3}.
   Werdykt Q-A wg litery.
2. **Phase 2 — bramka** (`Phase2_gate.py`+output): P2a próżnia zostaje
   (dryf ≤1e−10; jeśli ψ=1 nie jest punktem krytycznym pełnego E —
   użyj ψ* z P1a i odnotuj); P2b detektory: zasiany obiekt dolny
   i górny wykryte 1±0, czysta próżnia zero alarmów. FAIL ⟹ STOP.
3. **Phase 3 — Q-B RACHUNEK CENTRALNY** (`Phase3_relax_M911.py`
   +output+npz): 3 starty × 2 siatki (+dt/2 przy zdarzeniach);
   do ‖ψ̇‖∞≤1e−8 / nukleacji / t_max=200. Werdykty:
   Q-B-PASS-NUCLEATION (zbieżna: siatka×dt, ±1) / Q-B-PASS-STATIC
   (‖ψ−const‖∞≥0.05, zbieżność ≤5e−3; raportuj ψ_max vs 4/3) /
   Q-B-FAIL (wszystko do próżni) / Q-B-INCONCLUSIVE
   (BREAKDOWN-BOUNDARY = osobna kategoria deskryptywna, nie pozytyw).
   Deskryptywnie: los startu sieciowego.
4. **Phase 4 — Q-C** (tylko przy Q-B-PASS-STATIC;
   `Phase4_spectrum.py`+output): Hessian pełnego E (konsekwentnie
   z czynnikiem w), mody zerowe przed interpretacją, zbieżność
   ≤0.05·max(|ω²_min|,0.1), Q-C-PASS: ω²_min≥−1e−3. Przy
   PASS-NUCLEATION: charakterystyka kaskady bez progów.
5. `Phase_FINAL_close.md` (wzorzec frontmattera:
   `TGP/TGP_v1/research/op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md`),
   `NEEDS.md` (user-gated, drzewo §5 LOCKa), `README.md` (status+log).

## Wskazówki techniczne (nie zmieniają LOCKa)
- Runy >10 min w tle z checkpointami npz (limit ~55 min/proces —
  dziel na etapy); czekaj aktywnie, nie kończ tury między fazami.
- Detektory: scipy.ndimage.label z periodycznym sklejaniem.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
- Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE (commit robi sesja główna);
katalogi innych cykli tylko odczyt; formy (w,V,K) bez modyfikacji;
ZAKAZ dodawania podłóg/barier (sedno Q-A); detektory niezmienialne po
pierwszym biegu; INCONCLUSIVE ≠ pozytyw; rejestr WEJŚĆ flagowany.

## Raport końcowy
Werdykty Q-A/Q-B(/Q-C) z literą; krajobraz (punkty krytyczne, krzywizny,
wartości ρ_eff); tabela relaksacji per start×siatka (stan końcowy /
nukleacje dn-up / BOUNDARY / INCOMPLETE; ψ_max vs 4/3); los sieci;
pliki; korekty/incydenty; higiena. Pracuj do skutku.
