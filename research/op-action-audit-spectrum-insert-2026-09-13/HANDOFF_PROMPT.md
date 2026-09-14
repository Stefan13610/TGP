# HANDOFF_PROMPT — dla agenta-implementatora cyklu `op-action-audit-spectrum-insert-2026-09-13`

(Do wklejenia w całości jako zlecenie dla nowego agenta.)

---

Jesteś agentem-implementatorem cyklu TGP
`op-action-audit-spectrum-insert-2026-09-13`. Pracujesz w vaulcie
Obsidian; repo = TGP/TGP_v1. WSZYSTKIE ścieżki RELATYWNIE od rootu
vaulta. ZAKAZ zmiany cwd (`cd`) — pełne ścieżki; po KAŻDYM zapisie
pliku zweryfikuj `ls`, że zmaterializował się we właściwym miejscu
(znany artefakt zagnieżdżonych `TGP/TGP_v1/TGP/...`).

## Zadanie
Wykonaj cykl WYŁĄCZNIE wg LOCKa:
`TGP/TGP_v1/research/op-action-audit-spectrum-insert-2026-09-13/Phase0_balance.md`
— przeczytaj W CAŁOŚCI jako pierwszy, stosuj DOSŁOWNIE (anti-Lakatos:
zero zmian kryteriów/progów/form/rodziny A po starcie; korekty
wyłącznie dla udokumentowanego błędu implementacji — correction_note
PRZED użyciem wyniku, pierwotne outputy zachowane).

## Kontekst dziedziczony (przeczytaj PRZED kodem)
1. `TGP/TGP_v1/core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`
   (~148–380: eq:S-TGP-unified, eq:K-coupling-unified,
   eq:S-TGP-unified-M911-canonical, eq:V-M911 ~977) oraz
   `TGP/TGP_v1/core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex`
   (~511–520: rem/eq:vol-element-M911, pełna metryka ds²). Zamroź
   formy (w, V, K, metryka, √−g) w method_decisions Z CYTATAMI.
2. `TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/`:
   `Phase_method_decisions.md` §2 (rozstrzygnięcie odczytu B —
   DZIEDZICZYSZ, nie rozstrzygasz ponownie), `Phase_FINAL_close.md`
   (funkcjonał PRIMARY do gate'u P1b), `Phase3_relax_M911.py`
   (silnik gradient flow radialny do adaptacji pod więz pinu).
3. `TGP/TGP_v1/research/op-metametric-boundary-2026-09-01/Phase_FINAL_close.md`
   (Q1-POS: liczby −0.179 / +16156.6 — do odniesienia deskryptywnego
   w Phase 3, NIE do reprodukcji).

## Kolejność (LOCK §2)
0. `Phase_method_decisions.md` FROZEN: cytaty form; jawny zapis
   członu czasowego akcji i wariacji; schemat relaksacji z więzem
   ψ(0)=A (projekcja pinu po każdym kroku, udokumentowana); obsługa
   pasa ψ>4/3−1e−6 (BREAKDOWN-BOUNDARY, klasyfikacja).
1. **Phase 1 — kanonika** (`Phase1_canonical.py`+output): sympy —
   L=½Mψ̇²−½𝒦|∇ψ|²−𝒰 z literalnej akcji (M jest WYNIKIEM); π, H,
   równanie ruchu; P1b: tożsamość δH/δψ|_{π=0} ≡ δE_PRIMARY/δψ
   (simplify=0) — FAIL ⟹ STOP; P1c: sympy vs float 1e−12
   w {0.5, 1, 7/6, 1.3}.
2. **Phase 2 — Q-D1** (`Phase2_dispersion.py`+output): linearyzacja
   wokół ψ*=1, ω²(k)=[𝒦(1)k²+𝒰″(1)]/M(1); m², c_s²; znaki M,𝒦 na
   (0,4/3); odczyt A równolegle (odnotowanie). Werdykt wg litery.
3. **Phase 3 — Q-D2** (`Phase3_insert_cost.py`+output+json):
   P3a bramka (próżnia dryf ≤1e−10; pin A=1 ⟹ ΔE=0±1e−10; FAIL⟹STOP);
   P3b macierz A∈{0.50,0.70,5/6,7/6,1.25,1.30} × R∈{60,120} ×
   h∈{0.025,0.0125}; ΔE_insert(A;R,h)=E[ψ_A^relax]−E[ψ≡1] na
   IDENTYCZNYM pudle/siatce/brzegu; zbieżność h ≤5e−3 rel.; werdykt
   Q-D2-COST / Q-D2-CHANNEL / Q-D2-INCONCLUSIVE wg litery; tabela
   ΔE_insert(A) + profile (rozciągłość, ψ_max/min) + odniesienie
   definicyjne do Q1-POS.
4. `Phase_FINAL_close.md` (wzorzec frontmattera:
   `TGP/TGP_v1/research/op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md`),
   `NEEDS.md` (user-gated, drzewo §5 LOCKa), `README.md` (status+log).

## Wskazówki techniczne (nie zmieniają LOCKa)
- Rachunki radialne są tanie (1D) — batch w tle z logiem, czekaj
  aktywnie, nie kończ tury między fazami.
- E radialne: 4π∫r²[...]dr (trapez); Neumann przez ghost points;
  więz pinu: nadpisanie ψ[0]=A po każdym kroku + raport residuum.
- Python: `python` (numpy/scipy/sympy). Outputy do `Phase*_output.txt`.
- Sandbox: bez /dev/null, bez heredoc — wszystko przez pliki skryptów.

## Forbidden moves (egzekwowane)
Rdzeń `.tex`/STATE.md/git NIETYKANE (commit/STATE robi sesja główna);
katalogi innych cykli tylko odczyt; formy bez modyfikacji; M(ψ)
wyprowadzone, nie założone; ZAKAZ podłóg/barier; rodzina A/progi/więz
niezmienialne po pierwszym biegu; INCONCLUSIVE ≠ pozytyw; ZAKAZ
wnioskowania o barierze kreacji z ΔE_insert; rejestr WEJŚĆ flagowany.

## Raport końcowy
Werdykty Q-D1/Q-D2 z literą; jawne M(ψ),𝒦(ψ),𝒰(ψ), π, H, równanie
ruchu; ω²(k), m², c_s² (oba odczyty); tabela ΔE_insert(A;R,h)
ze zbieżnościami; pliki; korekty/incydenty; higiena. Pracuj do skutku.
