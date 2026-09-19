# AGENTS.md — konwencje TGP_v1 dla agentów

Repo = `TGP/TGP_v1` wewnątrz vaulta Obsidian. Projekt badawczy (fizyka teoretyczna, TGP),
prowadzony w **cyklach badawczych** z pre-rejestracją. Język roboczy: polski.

> **Ten plik opisuje JAK pracować.** Co jest aktualnie robione → [[STATE.md]] (czytaj w całości, jest krótki).

---

## 1. Budżet czytania — przeczytaj to zanim otworzysz jakikolwiek plik

Repo ma ~5600 plików / 142 MB. Kilka plików referencyjnych jest tak dużych, że przeczytanie
jednego zjada ćwierć okna kontekstu. **Reguła twarda:**

| Plik | Rozmiar | Jak czytać |
|---|---|---|
| `STATE.md` | ~22 KB | **CAŁOŚĆ** — to jest punkt startu |
| `AGENTS.md` (ten) | ~10 KB | **CAŁOŚĆ** |
| `CYCLES.tsv` | ~44 KB | **Grep lub head** — mapa wszystkich 263 cykli, 1 linia = 1 cykl |
| `research/<cykl>/README.md` | 2–4 KB | **CAŁOŚĆ** — karta cyklu z werdyktami |
| `research/<cykl>/HANDOFF_PROMPT.md` | ~5 KB | **CAŁOŚĆ** jeśli realizujesz ten cykl |
| `research/<cykl>/Phase0_balance.md` | 5–15 KB | **CAŁOŚĆ** jeśli realizujesz ten cykl (to LOCK) |
| `PREDICTIONS_REGISTRY.md` | 360 KB | **TYLKO Grep** |
| `INDEX.md` | 336 KB | **TYLKO Grep** (i tak stale od 2026-05-10) |
| `DEPENDENCIES.md` / `_REVERSE.md` | 172 KB | **TYLKO Grep** (auto-gen) |
| `TGP_FOUNDATIONS.md` | 74 KB | **TYLKO Grep** po numerze sekcji (§3.5.6 itd.) |
| `meta/STATE_ARCHIVE_*.md` | 80–400 KB | **TYLKO Grep** — historia sesji, nie czytać |
| `core/**/*.tex` | 10–200 KB | Grep po etykiecie (`rem:`, `prop:`, `eq:`), potem Read z `offset` |
| `main*.log`, `*_build*.log`, `main.aux` | do 724 KB | **NIGDY** — artefakty LaTeX-a, gitignored |

Nie czytaj pliku „na wszelki wypadek". Jeśli nie wiesz gdzie coś jest — Grep, nie Read.

**Zanim zaproponujesz nowy cykl — sprawdź, czy to już nie było liczone.** Grep odpowiada tylko
wtedy, gdy znasz słowo klucz; do pytania „czy ktoś już to badał" służy:

```
python tooling/similar.py "opis zagadnienia wlasnymi slowami"
python tooling/similar.py --like op-metric-pair-M911-2026-09-02   # cykle pokrewne
```

Zwraca ścieżki i werdykty (nie zdania), więc wynik jest weryfikowalny w źródle.

## 2. Warstwy repo — co czym jest i co wolno dotknąć

| Warstwa | Zawartość | Prawo zapisu |
|---|---|---|
| `axioms/` | Aksjomaty (notacja, substrat, różnica N0) | **user-gate** |
| `core/sek*/` | Korpus `.tex` — twierdzenia, propozycje, równania | **user-gate** |
| `TGP_FOUNDATIONS.md` | Referencja aksjomatyczna (W/E/P/H, dual-V §3.5) | **user-gate** |
| `research/op-*/` | Cykle badawcze — tu toczy się praca | swobodnie **we własnym cyklu** |
| `meta/` | Polityki, plany, handoffy, audyty, archiwa STATE | wg zadania |
| `audyt/` | Strukturalne długi (S01–S07, L01–L08, M, D, T) | wg zadania |
| `papers/`, `*.tex` w roocie | Publikacje (main / letter / companion) | wg zadania |
| `tooling/` | Skrypty pomocnicze i numeryczne | wg zadania |
| `partial_proofs/` | Dowody częściowe per temat | wg zadania |

**user-gate** = zmiana wymaga jawnej autoryzacji autora w tej sesji. Agent-implementator cyklu
**nigdy** nie dotyka `core/`, `STATE.md` ani gita — zgłasza to w `NEEDS.md` i kończy.

## 3. Protokół cyklu badawczego

Cykl = folder `research/op-<slug>-<YYYY-MM-DD>/`. Sekwencja:

```
Phase0_balance.md   LOCK — zero obliczeń. Pytania Q-*, kryteria, progi, detektory,
                    rodziny startów, pre-rejestrowane predykcje, drzewo decyzyjne.
HANDOFF_PROMPT.md   Samowystarczalne zlecenie dla agenta-implementatora.
Phase_method_decisions.md   FROZEN przed pierwszym biegiem: cytaty form, schemat, detektor.
Phase1..N_*.py + Phase*_output.txt    Realizacja. Output do pliku, nie do konsoli.
Phase_FINAL_close.md    Werdykty wg litery.
NEEDS.md            Co zostało otwarte / co wymaga decyzji autora (N1, N2, …).
README.md           Frontmatter + status + log cyklu.
```

**Cykl bez `Phase_FINAL_close.md` + `NEEDS.md` + dopisu w `README.md` NIE jest zakończony.**

### Anti-Lakatos (reguła nadrzędna)

Po starcie obliczeń **nie wolno** zmieniać: kryteriów, progów, detektorów, rodzin startów,
definicji werdyktu. Korekta dopuszczalna **wyłącznie** dla udokumentowanego błędu implementacji:
`Phase_correction_note_*.md` **przed** użyciem poprawionego wyniku, pierwotne outputy zachowane.

### Werdykty

- `Q-*-PASS` / `Q-*-FAIL` / `Q-*-INCONCLUSIVE` — zawsze **wg litery** zapisanej w LOCK-u.
- **`INCONCLUSIVE` ≠ pozytyw.** Nie przeformułowuj FAIL-a na „częściowy sukces".
- Wynik negatywny jest pełnoprawnym wynikiem. Nie ratuj hipotezy dopiskiem.
- `CONDITIONAL-ON-BRANCH` — twierdzenie ważne tylko przy jednej z konkurencyjnych konwencji.

### Native-first (obserwable, nie mimikra algebraiczna)

Reprodukcja wyniku MUSI być w **skali obserwabli** (arcsec, Hz, strain), **nie** w analogii
algebraicznej (parametry PPN/ppE). PPN to chart Willa, nie fizyka: γ jest natywne, β induced.
Dla cykli grawitacyjnych obowiązuje trójwarstwa L1 (natywne) / L2 (projekcja) / L3 (mapa
falsyfikacji) — [[meta/PPN_AS_PROJECTION.md]], [[meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]].

### S05 — jedno pole

TGP jest jednopolowe. Zanim zaproponujesz drugie pole / dodatkowy nośnik — sprawdź, czy
zjawisko nie jest sygnaturą tego samego pola. Naruszenie S05 wymaga user-gate.

### Dwa rejestry cykli

Nie myl ich w ocenie: **structural-emergence** (czy struktura w ogóle się wyłania, spójność)
vs **empirical-novelty** (czy przewiduje nową liczbę). Cykl spójnościowy nie przegrywa przez
brak nowej predykcji.

## 4. Statusy i WIP

`folder_status` w README cyklu to pole **maszynowe**, znormalizowane 2026-09-14 do **6 wartości**:

| wartość | znaczenie |
|---|---|
| `wip` | realnie w ruchu w tej / następnej sesji |
| `locked` | `Phase0_balance.md` zapisany (LOCK), realizacja nie ruszyła |
| `paused` | zaczęty, świadomie zamrożony |
| `parking` | pomysł zarejestrowany, Phase 0 nie istnieje |
| `closed` | zamknięty — **dowolny** werdykt |
| `legacy` | folder generacji sprzed 2026-05, status nigdy nie zapisany; referencja, nie WIP |

Treść opisowa żyje w **osobnych** polach: `verdict:` (werdykt słowny), `claim_status:` (A+/A/A−/B/C/D),
`legacy_status:` (poprzednia wartość `folder_status`, zachowana przy normalizacji).

**WIP-limit: max 5 cykli `wip`** (slot krytycznej ścieżki osobno). Historyczna pułapka, przed którą
chroni ta normalizacja: `active` dryfowało do „nigdy formalnie nie zamknięty" — w 2026-05 było 80,
w 2026-09 nadal 27, przy realnym WIP = 1. **Realny WIP jest wyłącznie w [[STATE.md]] §Active WIP.**

Mapa wszystkich cykli: `CYCLES.tsv` (generowana, 1 linia = 1 cykl). Regeneracja po zmianie statusu:

```
python tooling/build_cycles_index.py           # regeneruje CYCLES.tsv
python tooling/build_cycles_index.py --check   # exit 1 jesli nieaktualny
python tooling/normalize_folder_status.py      # dry-run normalizacji (--apply zapisuje)
```

Pełna polityka cyklu życia (WIP-limit, warunki przejścia): [[meta/CYCLE_LIFECYCLE.md]] — słownik
statusów w tamtym dokumencie pochodzi z 2026-05 i jest mapowany na powyższe 6 wartości.

## 5. Higiena techniczna (Windows + Obsidian)

- **Ścieżki relatywnie od rootu vaulta.** ZAKAZ `cd` — używaj pełnych ścieżek.
- Po każdym zapisie **zweryfikuj, że plik wylądował we właściwym miejscu** — znany artefakt
  zagnieżdżonych ścieżek `TGP/TGP_v1/TGP/TGP_v1/...` (czyszczony już dwukrotnie, patrz `.gitignore`).
- Python: `python` (numpy / scipy / sympy). Wersje: `tooling/requirements.lock`.
- Outputy numeryczne do `Phase*_output.txt` — nie do konsoli, nie do STATE.
- Długie ewolucje: batch w tle + checkpointy `.npz` + log; etapy ≤50 min.
- Linki Obsidian: `[[ścieżka/plik.md]]`. Naprawa zbiorcza: `tooling/fix_wikilinks.py`.
- Weryfikacja spójności: `tooling/check_status_drift.py`, `tooling/check_stale_cycles.py`,
  `tooling/validate_kickoff.py`, `tooling/build_deps_graph.py` (regen DEPENDENCIES).

## 6. Co robi sesja główna, a co agent-implementator

| | Sesja główna (z autorem) | Agent-implementator cyklu |
|---|---|---|
| Pisze `Phase0_balance.md` (LOCK) | TAK | nie |
| Wykonuje fazy obliczeniowe | zwykle nie | TAK |
| Dotyka `core/` | tylko za user-gate | **NIGDY** |
| Aktualizuje `STATE.md` | TAK | **NIGDY** |
| Operacje git | za zgodą autora | **NIGDY** |
| Czyta cudze cykle | TAK | tylko odczyt, wg listy w handoffie |

**Zasada handoffu:** `HANDOFF_PROMPT.md` musi być samowystarczalny — agent realizujący cykl
czyta handoff + `Phase0_balance.md` + wskazane 2–3 pliki źródłowe. **Nie czyta STATE ani INDEX.**

## 7. Start sesji — kolejność

1. `STATE.md` w całości → krytyczna ścieżka + WIP + 3 ostatnie sesje.
2. README cyklu, którego dotyczy zadanie.
3. Dopiero potem Grep po konkretach.

Jeśli zadanie jest badawcze i nie ma jeszcze LOCK-a — **najpierw LOCK (zero obliczeń)**,
potem realizacja. Nie licz przed pre-rejestracją.
