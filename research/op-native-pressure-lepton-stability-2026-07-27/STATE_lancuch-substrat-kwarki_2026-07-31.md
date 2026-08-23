# STATE — łańcuch substrat→…→kwarki (sesja 2026-07-30/31)

**Typ:** dokument-stan / punkt wznowienia. Uczciwy, po audycie adwersarialnym.
**Zasada:** rozdzielamy DERIVED / TAUTOLOGICZNE-CYRKULARNE / KONCEPT / OTWARTE.

---

## 0. Misja (niezmieniona)

Łańcuch od dołu: **substrat → wstążki → [pośrednie] → kwarki → solitony**.
Twarda kolejność. Reguły RAMY §5. Oś §0: **kolor i ładunek z JEDNEJ topologii — lock ma WYPAŚĆ.**
Filar (DERIVED, wcześniejszy): spin-½ z nawinięcia (π₃(S³)=ℤ + FR; niezależnie RP²+Berry).
Balony (e–p) = finalna bramka, NIE dotykana w tej sesji.

## 1. Mapa strzałek — status po sesji

| Strzałka | Status | Uwaga |
|---|---|---|
| substrat → wstążki | **domknięta strukturalnie** | 2T wymuszone (DERIVED); background-free = częściowo przemianowanie (audyt) |
| wstążki → kwarki | **rozpoczęta, §0 NIE domknięte** | ℤ₃-lock tautologiczny; nowe tylko: spin bezbarwny |
| kwarki → solitony | nietknięta | — |

## 2. Co jest GENUINNIE DERIVED (przetrwało niezależny audyt SL(2,3))

1. **Uniqueness 2T.** Wśród WSZYSTKICH skończonych podgrup SU(2) (π₁ dla SO(3)/H), **2T
   (binarna tetraedralna, rz.24) jest JEDYNĄ nieabelową z ℤ₃ w abelianizacji.** Dic_n→|ab|=4,
   2O→ℤ₂, 2I→triv, ℤ_n abelowe. Odporne na cutoff (n=1..15). → kolor=ℤ₃ wymusza target.
2. **Spin bezbarwny.** −1 ∈ komutatorze Q₈=[2T,2T], więc każdy charakter χ(−1)=1 → spin nie
   niesie ładunku ułamkowego. Wymuszony, nietrywialny, checkowalny fakt fizyczny.
3. **Struktura 2T=Q₈⋊ℤ₃:** spin=centralne −1 (w Q₈), kolor=ℤ₃ (iloraz). Spin+kolor z jednej grupy.
4. **Topologia M=S³/2T:** π₁=2T, π₂=0, π₃=ℤ. (π₂=0 → brak ucieczki monopolowej.)

## 3. Co jest TAUTOLOGICZNE / CYRKULARNE (NIE budować na tym jako osiągnięciu)

- **„Lock ładunek↔kolor w trzecich"** — automat: Hom(G,U(1))=Hom(G^ab,U(1)) dla każdej grupy;
  każdy płaski ładunek U(1) MUSI faktoryzować przez abelianizację. „Ładunek=funkcja koloru" to
  jedyne, co ładunek może robić. Wymaga DWÓCH wejść (kolor=ℤ₃, ładunek=holonomia), nie jednego.
- **„3" ładunku** = to samo ℤ₃ co kolor przez dualność Pontriagina; wejście, nie wyjście.
- **Symetria T substratu** „wymuszona trialnością" — cyrkularna (T bo żądamy koloru=ℤ₃).

## 4. KONCEPT (nieudowodnione identyfikacje, nie budować jako fakt)

- **Hom(π₁,U(1)) = ładunek elektryczny** — to fazy Aharonova-Bohma, NIE ładunek cechowania
  (prawo Gaussa). Skok. Cała strona ładunkowa wisi na tym.
- **„M=S³/2T jako reper samo-generowanej metryki" = background-free** — identyczna matematyka
  z zewnętrznym targetem; ontologia background-free jest interpretacją, nie derywacją.

## 5. OTWARTE (główne fronty)

1. **§0 lock realny** — wymaga ŁADUNKU Z NIEZALEŻNEGO ŹRÓDŁA (nie płaskiej holonomii, bo ta jest
   tautologiczna). Dwaj kandydaci:
   - **(EM)** Emergencja EM (RAMY §3): holonomia → realny dynamiczny ładunek (Gauss, prąd topol.).
   - **(GMN)** Ładunek z π₃-winding B/2 (Gell-Mann–Nishijima), niezależny od koloru; lock jako
     nietrywialna kwantyzacja kombinacji.
2. **Lokalizacja cząstki** (skończony rozmiar) — otwarta. Kanoniczny F-A jest runaway (brak
   solitonu). Status jak filar (liczba kwantowa ≠ zlokalizowana cząstka).
3. **Dokładne wartości kwarków** (+2/3, −1/3) — NIE wolno fitować do PDG.
4. **Generacje** (3 pokolenia ≠ 3 kolory) — nietknięte, otwarte.

## 6. ZDYSKWALIFIKOWANE w tej sesji (lekcje)

- **Płaski skyrmion + człon Skyrme (krok 2 dynamiki):** nie-natywny (płaskie tło, Skyrme wsadzony)
  → łamie §5.4. Derrick, który wymusił Skyrme, jest artefaktem płaskiego tła. `RIBBONS_stepB2_*`.
- **Pierwszy przebieg charge_lock:** bug etykietowania (fałszywe „NIE Z3") — zdiagnozowany, poprawiony.
- (wcześniejsze, z AUDYTU: bounce-hierarchy, lock z gołych warkoczy — patrz `AUDYT_KRYTYCZNY`).

## 7. Punkt wznowienia (rekomendacja)

Następna sesja: zaatakować **§0 przez niezależne źródło ładunku**. Preferencja: najpierw **(EM)
emergencja** — bez niej „ładunek" to tylko przemianowany kolor. Dopiero gdy holonomia jest realnym
ładunkiem, lock kolor↔ładunek przestaje być tautologią. Alternatywnie **(GMN)** z B/2.
NIE iść do kwarki→solitony ani do balonów, dopóki §0 nie ma nietautologicznej treści.

## 8. Indeks plików sesji (folder op-native-pressure-lepton-stability-2026-07-27)

**Wyniki (md):**
- `WYNIKI_substrat-wstazki_2T_2026-07-30.md` — uniqueness 2T
- `WYNIKI_stepB_dynamika-topologia_2026-07-30.md` — krok 1 (π₃/framing)
- `WYNIKI_stepB2_dynamika-rozmiar_2026-07-30.md` — krok 2 (ZDYSKWALIFIKOWANY, nie-natywny)
- `WYNIKI_s54_background-free_reconciliation_2026-07-30.md` — catch §5.4
- `WYNIKI_s54_domkniecie_frame-native_2026-07-30.md` — domknięcie §5.4 (reper), 2T=Q₈⋊ℤ₃
- `WYNIKI_charge-lock_Z3_2026-07-30.md` — lock ℤ₃ + KOREKTA PO AUDYCIE (czytać nagłówek)

**Skrypty (py):** `RIBBONS_group_compare.py`, `RIBBONS_uniqueness_2T.py`,
`RIBBONS_stepB_topology.py`, `RIBBONS_stepB2_skyrme_v3.py` (zdyskw.),
`RIBBONS_s54_frame_structure.py`, `RIBBONS_charge_lock.py`;
audyt: `AUDIT2_group_verify.py`, `AUDIT2_uniqueness_probe.py`.

**Źródła rdzenia:** `core/sek08a_akcja_zunifikowana.tex` (akcja: Φ skalarne, K=ψ⁴, √-g),
`core/sek08c_metryka_z_substratu.tex` (metryka algebraiczna z Φ = background-free),
`research/qm_spin/README.md` (filar spinu z reperu).

## 9. Jednozdaniowo

> Sesja wymusiła grupę koloru (2T, jedyna) i wyprowadziła bezbarwność spinu, ale audyt pokazał,
> że „lock kolor↔ładunek" via płaska holonomia jest tautologiczny; **oś §0 pozostaje otwarta i
> wymaga ładunku z niezależnego źródła (emergencja EM lub winding B/2)** — to jest właściwy
> punkt startu następnej sesji.
