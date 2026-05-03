---
title: "Plan wdrożenia nowego workflow `research/` (multi-agent + warstwa statusowa)"
date: 2026-05-02
type: plan
status: APPROVED v1 (decyzje człowieka 2026-05-02 zatwierdzone, sekcja 9)
author: Claudian (planner)
scope: organizacja struktury folderu `research/` + warstwa `meta/research/` + `meta/core/`; NIE rozwija fizyki
related:
  - "[[README.md]]"
  - "[[INDEX.md]]"
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
  - "[[meta/PLAN_DOMKNIECIA_MASTER.md]]"
  - "[[meta/PLAN_PUBLIKACJI_MASTER.md]]"
  - "[[meta/ROADMAP_v3.md]]"
  - "[[_archiwum/research_program_docs_2026-04/REDIRECT_PROGRAM_2026-04-19.md]]"
  - "[[_archiwum/research_program_docs_2026-04/TGP_STATUS_2026-04-19.md]]"
  - "[[research/closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]]"
  - "[[research/closure_2026-04-26/KNOWN_ISSUES.md]]"
tags:
  - meta
  - workflow
  - multi-agent
  - research-organization
  - status-layer
  - planning
---

# Plan wdrożenia nowego workflow `research/` (multi-agent + warstwa statusowa)

> **Cel planu:** Przygotować bezpieczne, krokowe wdrożenie nowej warstwy
> statusowej i meta-organizacyjnej dla folderu `research/`, tak aby wielu
> agentów mogło pracować równolegle bez deptania sobie po nogach, a wyniki
> jednego folderu były wykrywane jako wejście dla innego.
>
> **NIE jest celem:** rozwijanie fizyki, dopisywanie wyprowadzeń, ocena czy
> teoria jest prawdziwa. Plan jest **organizacyjny**.
>
> **NIE wykonuj:** żadnych `mv`, masowych edycji rdzenia ani zmian w
> `core/`, `axioms/`, `partial_proofs/`, `papers_external/` bez akceptacji.

---

## 0. Najważniejsze ostrzeżenia, przeczytać przed czymkolwiek

1. **Kanoniczny `research/` = `TGP/TGP_v1/research/`** (decyzja człowieka,
   2026-05-02). Vault-rootowy `./research/` z 11 folderami `op-*` jest
   traktowany jako **out-of-scope tego planu**. W Sesji 1 agent dodatkowo
   raportuje stan vault-rootowego `./research/` (diff z canonical), żeby
   człowiek mógł później zdecydować, czy go usunąć / zsynchronizować.
   Plan **nie modyfikuje** `./research/` (vault root) w żadnej sesji 1–9.

2. **Linki są kruche.** `INDEX.md` w `TGP/TGP_v1/INDEX.md` zawiera dziesiątki
   wpisów typu `[[research/op-xi-photon-ring/Phase3_results.md]]`.
   Każde fizyczne przeniesienie folderu rozwala te linki. Z tego powodu
   plan **preferuje warstwę meta + frontmatter `tgp_status`** zamiast
   fizycznych nadfolderów statusowych. Fizyczne przenoszenie jest opcją,
   nie domyślnym ruchem.

3. **Już istnieje `meta/` i `_archiwum/`** na poziomie `TGP/TGP_v1/`.
   - `TGP/TGP_v1/meta/` — istniejące audyty, plany, propozycje.
   - `TGP/TGP_v1/_archiwum/` — historyczne plany, stare analizy.
   Nowa warstwa meta dla `research/` musi być **kompatybilna** z tymi
   istniejącymi konwencjami, a nie z nimi konkurować.

4. **Subagenty mają historię over-claimingu.** Patrz
   [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] (commit `74394a8` z czterema
   "FULL CONVERGENCE" cykli, które okazały się polluted).
   Z tego wynika, że plan **nie wolno** budować na ślepym zaufaniu do
   werdyktów agenta — każda klasyfikacja statusowa folderu wymaga
   weryfikowalnego źródła (PASS/FAIL ze skryptem + plik wynikowy + audyt
   dat) oraz, w przypadku promocji do core, decyzji człowieka.

5. **Awans do core ≠ przeniesienie pliku.** "Core-promoted" w tym planie
   znaczy: wynik został uznany za rdzeń teorii (występuje w `main.tex` /
   `core/` / `axioms/` / `partial_proofs/`), ale folder źródłowy nadal
   żyje w `research/` jako historia/uzasadnienie. **Nie wolno bez decyzji
   człowieka** dotykać `core/`, `axioms/`, `papers_external/`, `main.tex`,
   `partial_proofs/`.

---

## 1. Audyt — co już wiemy o `research/` (wkład w plan)

Skrócony obraz, który spinam, żeby kolejne sesje miały punkt zaczepienia.
**Ten audyt jest informacyjny, nie klasyfikacyjny.** Pełna klasyfikacja
przychodzi w Sesji 4.

### 1.1 Skala

- `TGP/TGP_v1/research/`: ~95 podfolderów + kilkanaście top-level plików
  (`graph_*.gexf/png`, `tgp_dependency_graph.py`, `TGP_STATUS_*.md`,
  `REDIRECT_PROGRAM_*.md`, `NEW_DIRECTIONS_*.md`).
- Z tego ~56 folderów ma prefiks `op-` (kanoniczna konwencja "open problem").
- Reszta to: `qm_*` (kwantowe), `closure_*` (sesje domknięć),
  `external_review_*` (wzajemne audyty), tematyczne (`hubble_tension`,
  `galaxy_scaling`, `em_from_substrate`, `nbody`, `cosmo_tensions`,
  `metric_ansatz`, `continuum_limit`, …) i historyczne (`op6`, `op7`,
  `op1-op2-op4`, `uv_completion`).

### 1.2 Wzorce wewnątrz folderów (heurystyki, nie reguły)

- **3-fazowy `op-*` z `program.md`:** np. `op-uv-as-ngfp/`,
  `op-alpha-fine-structure/`. Pliki: `program.md`,
  `Phase{1,2,3}_setup.md`, `Phase{1,2,3}_results.md`, skrypty + `.txt`.
  Te foldery są zwykle **active** lub **core-ready**.
- **Topic + README:** np. `qm_born_rule/`, `hubble_tension/`, `s8_tension/`,
  `continuum_limit/`. Często jeden README + 1–2 skrypty + `.txt`.
  Te są zwykle **needs-bridge** lub **sandbox**.
- **Closure-aggregator:** `closure_2026-04-26/` z czterema podzamknięciami
  (`sigma_ab_pathB/`, `f_psi_principle/`, `Lambda_from_Phi0/`,
  `alpha_psi_threshold/`) + `CLOSURE_*_SUMMARY.md` + `KNOWN_ISSUES.md`.
  To jest **specjalny typ** — meta-folder wieloośkowych domknięć.
- **External review:** `external_review_2026-04-25/` — wzajemne audyty
  gotowości materiału do publikacji. **Należy do osobnej kategorii
  `audit`/`review`**, niemiałczone z badaniem pierwotnym.
- **Legacy `opN`:** `op6/`, `op7/`, `op1-op2-op4/`. Mają mieszany status —
  częściowo **core-promoted** (op-7 zamknięte, weszło do `main.tex`),
  częściowo nadal **active** podzamknięcia.
- **Top-level pliki w `research/`:** `TGP_STATUS_*.md`,
  `REDIRECT_PROGRAM_*.md`, `NEW_DIRECTIONS_*.md` — to są dokumenty
  programowe, nie foldery. **Powinny zostać w `research/` w roli
  index/program**, ale mogą zostać dowiązane z `meta/`.

### 1.3 Sygnały statusowe już obecne w plikach

W treści plików widać już bogatą semantykę statusów (z których plan czerpie):

- `LOCKED`, `CLOSED`, `OPEN`, `STRUCTURAL`, `DERIVED`, `PARTIALLY DERIVED`,
  `STRUCTURAL HINT`, `BLOCKED`, `DEFERRED`, `WITHDRAWN`, `FALSIFIED`.
- "Phase X CLOSED N/N PASS" — werdykt agentowy z liczbami.
- "Promotion → DERIVED", "scenario A RATIFIED", "Path B PRIMARY".
- Ostrzeżenia: "⚠ B3-locked 2026-05-01", "⚠ C9-PDG-update".

Plan **nie tworzy nowego słownika** — wprowadza tylko cienką warstwę
mapującą tę istniejącą semantykę na 9 znormalizowanych pól w
frontmatterze + poziomy L0–L4.

### 1.4 Co już istnieje — nie duplikujemy

| Istniejący artefakt | Lokalizacja | Czy zastępujemy? |
|---|---|---|
| Audyty całości | `meta/AUDYT_TGP_*.md` | NIE — uzupełniamy o `meta/research/`. |
| Master plan domknięć | `meta/PLAN_DOMKNIECIA_MASTER.md` | NIE — `IMPACT_MATRIX.md` go uzupełnia. |
| Plan publikacji | `meta/PLAN_PUBLIKACJI_MASTER.md` | NIE. |
| Roadmapa | `meta/ROADMAP_v3.md` | NIE. |
| Słownik formalizmu | `axioms/notacja/slownik_formalizmu.tex` | TAK, **delegujemy** — `meta/research/GLOSSARY.md` tylko linkuje + dodaje pojęcia organizacyjne. |
| Indeks repo | `INDEX.md` | NIE — `FOLDER_STATUS_INDEX.md` uzupełnia. |
| Graph zależności | `tgp_dependency_graph.py`, `graph_*.gexf` | NIE — `RESEARCH_BUS.md` to ludzkojęzyczny komplement. |
| Dependencies forward/reverse | `DEPENDENCIES.md`, `DEPENDENCIES_REVERSE.md` | NIE — input dla `IMPACT_MATRIX.md`. |
| Archiwum | `_archiwum/` | NIE — pozostaje docelową lokalizacją historycznych analiz. Plan **nie tworzy** drugiego archiwum. |

---

## 2. Docelowa filozofia

### 2.1 Każdy folder w `research/` ma:

1. Jasny **status** (`folder_status`).
2. Jasny **poziom dojrzałości** L0–L4 lub status pomocniczy.
3. `README.md` — opis zawartości i intencji.
4. `NEEDS.md` — czego brakuje (otwarte luki, pytania, dane).
5. `FINDINGS.md` — eksportowalne wyniki (numeryczne i analityczne).
6. Pole `core_compatibility` — czy zgodny z aktualnym core.
7. Pole `may_edit_core` — czy folder może wpływać na core.
8. Pole `last_reviewed_against_core` — kiedy ostatni raz porównano.
9. Pola `depends_on` / `impacts` / `open_bridges` — relacje międzyfolderowe.

### 2.2 Warstwa meta jako single source of truth

Plan **preferuje warstwę meta** (frontmatter + `meta/research/*.md`) nad
fizyczną reorganizacją. Dlaczego:

- Mniej zerwanych wikilinks i `\input{}`/`\ref{}` (pamiętaj — w `INDEX.md`
  jest 152 wikilinków i 1436 `\ref{}` edges; każdy przeniesiony folder
  to potencjalna eksplozja sieroctw).
- Status folderu zmienia się częściej niż jego fizyczna pozycja.
- Człowiek może otrzymać "view" przez raporty (`FOLDER_STATUS_INDEX.md`),
  bez konieczności renamowania ścieżek.

Fizyczne przeniesienia są zarezerwowane dla kategorii **archive** /
**sandbox** / **needs-migration**, gdzie rozróżnienie ma realny sens
(człowiek ma "nie szukać tu aktualnej teorii").

---

## 3. Proponowana struktura `research/` — krytyczna ocena strawmana

Strawman z briefu (active / core-ready / core-promoted / needs-bridge /
needs-migration / phenomenology / numerical / archive / sandbox) jest dobry,
ale po audycie repo wymaga modyfikacji.

### 3.1 Co chcę **zostawić** ze strawmana

- `archive/` — sensowne **tylko** dla rzeczy, które nie pasują do
  istniejącego `_archiwum/`. Generalnie wolałbym przesyłać historyczne
  rzeczy do `_archiwum/`. Wartością `research/archive/` jest co najwyżej
  podtrzymanie izolacji starych podejść z poziomu badawczego (np. `op6`
  pre-pivot v1).
- `sandbox/` — TAK, potrzebne. Obecnie nie ma jasnego miejsca na
  eksperymenty agentów ("a może by tak…").
- `needs-bridge/` — TAK. Pasuje do tego, co audyt 2026-05-01 nazywa
  "STRUCTURAL HINT" / "PARTIALLY DERIVED".
- `needs-migration/` — TAK. Z `AUDYT_TGP_2026-05-01.md` wynika, że
  kilka folderów (M9.x, sek08c) wymaga aktualizacji po pivot M9.1''.

### 3.2 Co chcę **odrzucić** lub zmienić

- **Nie tworzę** `research/core-promoted/` i `research/core-ready/` jako
  fizycznych nadfolderów. Zamiast tego — `tgp_status.folder_status` w
  YAMLu + `meta/research/CORE_CANDIDATES.md`. Powód: 56 folderów `op-*`
  z linkami w `INDEX.md` nie zniesie hurtowego renamowania, a
  "core-promoted" jest właściwością relacyjną folderu wobec `main.tex`,
  nie jego fizyczną lokalizacją.
- **Nie tworzę** dichotomii `phenomenology/` vs `numerical/`. W tym repo
  to nie jest naturalna oś — większość folderów `op-*` ma JEDNO i drugie
  (np. `op-eht`, `op-rho1-71Ge-cross-section`). Zamiast tego dodaję
  pomocnicze tagi `kind: derivation | numerical | phenomenology |
  audit | review | bridge | program-doc` w YAMLu.
- **Audyty zewnętrzne** (`external_review_*`) zasługują na osobny
  kind/folder pomocniczy, niezmieszany z aktywnym researchem.

### 3.3 Wybrana struktura — wariant minimalny (decyzja człowieka, 2026-05-02)

```
TGP/TGP_v1/research/
├── (95 folderów badawczych, NIEPRZENIESIONYCH; każdy z YAML status block)
├── _sandbox/      ← NOWY: bezpieczna piaskownica dla agentów (NIE part teorii)
├── _archive/      ← NOWY: research-level archiwum dla obsoletnych podejść,
│                    FALSIFIED, WITHDRAWN, superseded prób badawczych
│                    (różne od `_archiwum/` na poziomie TGP_v1, który
│                    trzyma stare plany / sesje / analizy meta)
├── TGP_STATUS_*.md, REDIRECT_PROGRAM_*.md, NEW_DIRECTIONS_*.md (zostają)
└── graph_*.png, graph_*.gexf, tgp_dependency_graph.py (zostają)

TGP/TGP_v1/meta/research/   ← NOWY (warstwa meta dla badań)
├── RESEARCH_BUS.md
├── CANDIDATE_BRIDGES.md
├── CORE_CANDIDATES.md
├── IMPACT_MATRIX.md
├── AGENT_PROTOCOL.md
├── FOLDER_STATUS_INDEX.md
├── MIGRATION_LOG.md
├── GLOSSARY.md            (cienki adapter do axioms/notacja/slownik_formalizmu.tex)
├── _examples/             (3 podglądy z Sesji 2)
└── templates/             (4 szablony z Sesji 2)

TGP/TGP_v1/meta/core/       ← NOWY (warstwa meta dla rdzenia, sekcja 5A)
├── CORE_INVENTORY.md
├── CORE_HOTSPOTS.md
└── CORE_INTAKE.md
```

**Status folderu jest w YAML**, nie w ścieżce. To jest **kontrakt** tego
planu. Jedyny wyjątek to fizyczne `_sandbox/` i `_archive/` — te
istnieją jako wyraźny sygnał "tu nie szukasz aktualnej teorii".

### 3.4 Reguły wpisu do `research/_archive/`

`_archive/` w `research/` jest **rygorystyczny**. Folder trafia tam
**wyłącznie** gdy spełnia jeden z warunków:

1. Treść jest oznaczona w samym folderze jako `WITHDRAWN`, `FALSIFIED`,
   `OBSOLETE`, `SUPERSEDED` w sposób potwierdzony co najmniej dwoma
   źródłami (np. notka w `README.md` + wpis w `KNOWN_ISSUES.md`
   lub w jakimś `AUDYT_*.md`).
2. **Lub** człowiek explicite zdecydował o archiwizacji.

Zasadniczo do `_archive/` **NIE** trafiają:
- foldery aktualnie używane jako historia/uzasadnienie (zostają in place
  z `folder_status: core-promoted`),
- foldery z otwartym researchem (zostają in place z `active`/`needs-bridge`),
- foldery typu `closure_2026-04-26/` (closure-aggregator, in place),
- legacy `op6/`, `op7/`, `op1-op2-op4/` (zbyt gęsto linkowane, in place
  z odpowiednim `folder_status`).

Reguła **różnicy** wobec `_archiwum/` (TGP_v1 level):
- `_archiwum/` (TGP_v1) = stare **plany / sesje / analizy meta** (typ:
  `meta` historical).
- `research/_archive/` = stare **ścieżki badawcze / próby derivacyjne**
  (typ: `research` historical).

Każde przeniesienie do `_archive/` jest fizycznym `git mv` + patchem
linków (zob. Sesja 8.5 niżej, "selektywna archiwizacja").

---

## 4. YAML status block — schemat do każdego `README.md`

Dodawany **na samym początku** istniejącego README. Jeśli README ma już
frontmatter, sekcja `tgp_status:` jest doklejana jako pole; nie zastępuje
istniejących pól.

```yaml
---
# (istniejący frontmatter, jeśli był)
tgp_status:
  folder_status: active        # active | core-ready | core-promoted
                               # | needs-bridge | needs-migration
                               # | audit | review | program-doc
                               # | sandbox | archive
  level: L2                    # L0 | L1 | L2 | L3 | L4 | mixed | unknown
  kind: derivation             # derivation | numerical | phenomenology
                               # | bridge | audit | review | program-doc
                               # | closure-aggregator
  core_compatibility: current  # current | partial | stale | broken | unknown
  last_reviewed_against_core: 2026-05-02   # ISO date or "unknown"
  may_edit_core: false         # bool — czy ten folder ma prawo dotykać core
  exports_findings: true       # bool — czy plik FINDINGS.md istnieje i jest świeży
  has_needs_file: true         # bool
  has_findings_file: true      # bool
  open_bridges:                # listy lemmatów / luk wymagających mostu
    - "M9.1'' volume element re-derivation"
  depends_on:                  # foldery, których wyników ten folder potrzebuje
    - "research/op-phase3-uv-completion"
    - "research/closure_2026-04-26/sigma_ab_pathB"
  impacts:                     # foldery, które mogą przyjąć wyniki tego folderu
    - "research/op-uv-as-ngfp"
  source_of_status:            # WERYFIKOWALNE źródła statusu (anty-overclaim)
    - "Phase3_results.md (94/97 PASS, 2026-04-25)"
    - "meta/AUDYT_TGP_2026-05-01.md §A1"
  promoted_to_core: null       # null | "main.tex §X.Y" | "core/sek0Z_*" | "axioms/*"
---
```

### 4.1 Definicje poziomów L0–L4

- **L0 — pomysł / sandbox**: hipoteza, brak skryptów, brak weryfikacji.
- **L1 — szkic numeryczny**: skrypt biegnie, wyniki istnieją, brak
  analitycznej interpretacji.
- **L2 — analityczny szkic**: jest derivation albo bridge, ale bez pełnej
  weryfikacji, lub testy częściowe.
- **L3 — domknięty wynik**: derivation + numeryka zgodne, brak otwartych
  bridge'y, audyt dat OK; gotowe do `core-ready`.
- **L4 — promowany do core**: wynik jest w `main.tex` / `core/` /
  `axioms/`; folder źródłowy żyje jako historia.
- **mixed** — folder zawiera podzbiory na różnych poziomach (np.
  `closure_2026-04-26/` ma L3 podfoldery).
- **unknown** — przed klasyfikacją lub po nieudanym audycie.

### 4.2 Reguła "anty-overclaim"

`folder_status: core-promoted` lub `level: L4` jest dopuszczalne **tylko**
gdy `promoted_to_core` zawiera konkretną referencję do pliku core
(`main.tex §X.Y` itp.) **i** ta referencja jest weryfikowalnie obecna w
repozytorium. W przeciwnym razie agent **musi** użyć `core-ready` + pole
`source_of_status` z konkretnymi PASS-ami.

To jest bezpośredni wniosek z [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]].

---

## 5. Warstwa meta — pliki w `meta/research/`

Kompatybilna z istniejącym `meta/`, ale **nie miesza się** z nim
(nie psuje istniejących planów / audytów). Wszystkie ścieżki względne do
`TGP/TGP_v1/`.

### 5.1 `meta/research/RESEARCH_BUS.md`
Tablica ogłoszeń międzyfolderowych. Każdy wpis to "ten folder
wyprodukował wynik X w dacie D, który jest potencjalnie wejściem dla
folderów {A, B, C}". Format tabeli + wikilinks. Aktualizowany przez
agentów po Sesji 6 jako standard. **Nie** powtarza zawartości
`FINDINGS.md` — tylko linkuje i dodaje broadcast metadata.

### 5.2 `meta/research/CANDIDATE_BRIDGES.md`
Lista luk analitycznych ("jest numerycznie, brak bridge'a"). Pochodzi
z agregowania `NEEDS.md` ze wszystkich folderów + ręcznych wstawek z
audytu. Format:

```
| ID    | Skąd brak (folder) | Czego brakuje (lemma)        | Kto może dostarczyć | Status |
| BR-01 | op-eht             | analityczny rdzeń photon ring| op-xi-photon-ring   | OPEN   |
```

### 5.3 `meta/research/CORE_CANDIDATES.md`
Lista folderów oznaczonych jako `core-ready` z notatką, **co** ma trafić
do core, **gdzie** w `main.tex`/`core/` i **kiedy** decyzja człowieka jest
wymagana. Bez tego pliku **nie wolno** modyfikować core.

### 5.4 `meta/research/IMPACT_MATRIX.md`
Macierz "zmiana w folderze X może wpłynąć na foldery {…}". Bazuje na:
- `tgp_status.depends_on` / `impacts` w YAMLach,
- istniejącym `DEPENDENCIES.md` / `DEPENDENCIES_REVERSE.md`,
- ręcznych anotacjach z audytu (`AUDYT_TGP_2026-05-01.md` daje gotową
  listę propagacji M9.1'' → sek08c, M9.2, M9.3 itd.).

Plik jest punktem wejścia dla "regression review" w Sesji 8.

### 5.5 `meta/research/AGENT_PROTOCOL.md`
Reguły dla agentów (multi-agent rules of engagement):
- jak otwierać sesję na folderze (lock plik / odczyt + zapis YAML),
- co wolno, czego nie wolno,
- co znaczy `may_edit_core: false`,
- jak nadawać status (z weryfikowalnym źródłem),
- jak zgłaszać wynik na `RESEARCH_BUS`,
- jak używać `_sandbox/`,
- jak postępować przy konflikcie ze stanem core (zawsze: stop →
  `needs-migration` → human),
- ban na "FULL CONVERGENCE" bez weryfikowalnego ledger entry.

To jest **bezpośrednia odpowiedź** na incydent z commit `74394a8`.

### 5.6 `meta/research/FOLDER_STATUS_INDEX.md`
Wygenerowany (i potem manualnie poprawialny) widok wszystkich folderów
z ich `tgp_status`. Format:

```
| Folder | folder_status | level | kind | core_compatibility | last_reviewed |
| op-alpha-fine-structure | core-ready | L3 | derivation | current | 2026-05-02 |
| qm_born_rule            | needs-bridge | L1 | derivation | unknown | -- |
```

Generowany przez prosty skrypt, który czyta wszystkie YAMLe. Skrypt
zostanie napisany w Sesji 4. Plik jest **stale aktualizowany**, więc
ma być deterministyczny (sortowanie alfabetyczne).

### 5.7 `meta/research/MIGRATION_LOG.md`
Każda fizyczna zmiana lokalizacji folderu (do `_sandbox`, do
`needs-migration`, do `_archiwum/`) musi być zalogowana z datą,
ID-iem decyzji człowieka i listą zaktualizowanych linków.

### 5.8 `meta/research/GLOSSARY.md`
**Cienki adapter** do istniejącego `axioms/notacja/slownik_formalizmu.tex`.
Zawiera tylko terminy organizacyjne (status, level, kind, bridge,
research bus, …) + linki do słownika fizycznego. Nie powtarzamy
fizyki.

---

## 5A. Warstwa `meta/core/` (decyzja człowieka, 2026-05-02)

`meta/research/` opisuje **źródło** (ścieżki badawcze). `meta/core/`
opisuje **cel** (rdzeń teorii). Bez warstwy `meta/core/` agenci nie mają
gdzie zapisywać "co już jest w core" w sposób wyszukiwalny — muszą za
każdym razem grepować `main.tex` i `core/sek*.tex`. To była realna
dziura, która ułatwiła incydent `74394a8` (agenci awansowali do `DERIVED`
bez weryfikowalnej odpowiedzi "tak, ale **gdzie** w core to jest?").

`meta/core/` zawiera 3 pliki:

### 5A.1 `meta/core/CORE_INVENTORY.md`

**Statyczna mapa** rdzenia: dla każdej sekcji `core/sek0X_*` /
`axioms/*` / `partial_proofs/*` jednolinijkowy opis, jakie wyniki tam
żyją, plus link "promoted from research/...". Format:

```
| Core path                          | Topic                          | Promoted from           | Last verified |
| core/sek08c_metryka_z_substratu/   | Effective metric M9.1''        | research/op-newton-momentum | 2026-04-25 |
| axioms/M1A_substrate_v2/           | Substrate axiom v2             | research/op6/v2_pivot_summary.md | 2026-04-24 |
| partial_proofs/koide_fp/           | Koide K=2/3 algebraic identity | research/cabibbo_correction/ + paper_lepton_masses/ | 2026-04-15 |
```

To jest **jedyne** miejsce, w którym agent szuka odpowiedzi "czy X jest
już w core". Wpisy są dodawane w Sesji 7 (gdy folder otrzymuje
`folder_status: core-promoted`).

`CORE_INVENTORY.md` jest **read-mostly** dla agentów — modyfikacja
wymaga tej samej weryfikacji co reguła "anty-overclaim" (sekcja 4.2):
referencja do `core/...` lub `main.tex §X.Y` musi być sprawdzalna.

### 5A.2 `meta/core/CORE_HOTSPOTS.md`

**Living checklist** otwartych sprzeczności / luk migracyjnych w samym
core. Pochodzi bezpośrednio z [[meta/AUDYT_TGP_2026-05-01.md]]
(8 CRITICAL + 12 HIGH) i jest na bieżąco redukowana, gdy hotspot zostanie
zaadresowany. Format:

```
| ID  | Hotspot                                              | Source         | Status   | Linked needs-migration folders |
| H1  | sek08c — 4 wzajemnie sprzeczne formy metryki         | AUDYT 2026-05-01 §A1 | OPEN | research/op-newton-momentum, research/metric_ansatz |
| H2  | √(-g) dla M9.1'' — używana stara forma c·ψ           | AUDYT 2026-05-01 §A2 | OPEN | research/op-newton-momentum (M9.2, M9.3) |
| ... |                                                      |                |          |                                |
```

Każdy hotspot z `OPEN` automatycznie generuje wpis "needs-migration"
w `IMPACT_MATRIX.md` dla powiązanych folderów badawczych. Sesja 7
**konsumuje** `CORE_HOTSPOTS.md` przy klasyfikacji `needs-migration`.

### 5A.3 `meta/core/CORE_INTAKE.md`

**Protokół promocji** research → core. Twardy gate dla człowieka.
Definiuje:

- co znaczy "core-ready" (kandydat do core, ale jeszcze nie w core),
- jak agent przygotowuje **wniosek o promocję** (PR-like artefakt),
- co musi zawierać wniosek (referencja do `FINDINGS.md` źródła +
  proponowana lokacja w core + diff w core, jeśli zmiana wymaga
  modyfikacji istniejących sek*),
- **człowiek jest jedynym aktorem**, który może zaakceptować promocję
  i wykonać edycję core,
- po akceptacji: aktualizacja `CORE_INVENTORY.md` + zmiana
  `folder_status` źródłowego folderu na `core-promoted` + wpis
  `promoted_to_core` z konkretną referencją.

`CORE_INTAKE.md` jest tworzony w Sesji 3 jako protokół (statyczny);
poszczególne wnioski o promocję trafiają jako osobne pliki
`meta/core/intake/INTAKE_<folder>_<date>.md` w Sesji 7+.

### 5A.4 Czego `meta/core/` świadomie NIE robi

- Nie powtarza `INDEX.md` (który jest hubem nawigacyjnym całego repo).
- Nie powtarza `DEPENDENCIES.md` (który jest auto-generowanym grafem).
- Nie zastępuje `meta/AUDYT_*` (które są snapshot'ami punkt-w-czasie).
- Nie modyfikuje plików w `core/` / `axioms/` / `partial_proofs/`
  bezpośrednio. Wszystkie zmiany w core są **wyłącznie** przez gate
  `CORE_INTAKE.md` z akceptacją człowieka.

---

## 6. Multi-session execution plan

Domyślnie wszystkie sesje są **read-only + write w `meta/research/` +
write `README.md`/`NEEDS.md`/`FINDINGS.md` w folderach badawczych**.
Edycja core, `axioms/`, `partial_proofs/`, `papers_external/`,
`main.tex`, `tgp_letter.tex`, `tgp_companion.tex` — **zakazana** we
wszystkich sesjach poza wyraźnie oznaczonymi gate'ami z aprobatą człowieka.

### Sesja 1 — Audyt struktury i decyzja o wariancie

**Cel:** zebrać twarde dane o `research/` i poprosić człowieka o decyzję
między wariantem minimalnym (pkt 3.3) a rozszerzonym (3.4).

**Wejście:**
- `TGP/TGP_v1/research/` (cały)
- `./research/` (vault root) — porównanie z TGP_v1/research/
- `INDEX.md`, `DEPENDENCIES.md`, `DEPENDENCIES_REVERSE.md`
- `meta/AUDYT_TGP_2026-05-01.md`, `meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md`

**Co agent ma zrobić:**
1. Wygenerować płaską listę folderów `research/` z 5 metrykami:
   - czy ma `program.md`,
   - czy ma `Phase{1,2,3}_*.md`,
   - czy ma `README.md` z YAMLem,
   - czy ma jakieś `PASS/CLOSED/DERIVED` w treści,
   - data ostatniej modyfikacji (`mtime` najnowszego pliku w folderze).
2. Zidentyfikować **duplikaty** między `./research/` a
   `TGP/TGP_v1/research/`. Podać diff (które pliki różnią się treścią).
3. Wypisać **kandydatów do różnych klas statusowych** WSTĘPNIE
   (heurystycznie, bez nadawania).
4. Sprawdzić, ile linków w `INDEX.md` / `DEPENDENCIES.md` zostanie
   złamanych dla każdego z dwóch wariantów (3.3 vs 3.4).
5. Wygenerować raport `meta/research/AUDIT_RESEARCH_S1.md`.

**Czego agentowi nie wolno:**
- Przenosić plików.
- Modyfikować plików w folderach badawczych.
- Edytować `core/`, `axioms/`, `partial_proofs/`, `main.tex`.
- Klasyfikować folderów ostatecznie.

**Rezultat:**
- `meta/research/AUDIT_RESEARCH_S1.md`
- Raport diff `research/` vs `TGP/TGP_v1/research/`
- Lista pytań decyzyjnych dla człowieka

**Pliki, które powstają / aktualizują się:**
- `meta/research/AUDIT_RESEARCH_S1.md` (nowy)

**Kryterium zakończenia:**
- Człowiek odpowiedział na 3 pytania:
  (1) który `research/` jest kanoniczny,
  (2) wariant minimalny czy rozszerzony,
  (3) czy `_archiwum/` (TGP_v1) jest jedynym archiwum, czy potrzebny też
      `research/_archive/`.

**Ryzyka:**
- Agent może źle zliczyć linki w `INDEX.md`. Mitygacja: skrypt + manual
  ręczny sample.

---

### Sesja 2 — Zaprojektowanie statusów, szablonów, schematu YAML

**Cel:** sfinalizować konwencję statusu (sekcja 4) jako szablony i
przykłady.

**Wejście:**
- Niniejszy plan (sekcja 4)
- `meta/research/AUDIT_RESEARCH_S1.md` z Sesji 1
- 3 reprezentatywne foldery (np. `op-uv-as-ngfp/`, `qm_born_rule/`,
  `closure_2026-04-26/`) jako case studies

**Co agent ma zrobić:**
1. Stworzyć szablony:
   - `meta/research/templates/README.template.md`
   - `meta/research/templates/NEEDS.template.md`
   - `meta/research/templates/FINDINGS.template.md`
   - `meta/research/templates/STATUS_BLOCK.yaml`
2. Dla 3 case studies — **prywatnie** w sandbox sketcha — pokazać, jak
   te szablony wyglądałyby zaaplikowane (output do
   `meta/research/_examples/op-uv-as-ngfp_README.preview.md` itd.).
   **Nie modyfikuje oryginalnych folderów.**
3. Dopisać do `AGENT_PROTOCOL.md` (na razie jako stub) regułę "jak
   tworzyć NEEDS / FINDINGS bez fabrykowania PASS".
4. Sformalizować definicje L0–L4 + reguły "anty-overclaim".

**Czego agentowi nie wolno:**
- Dotykać oryginalnych README w folderach badawczych.
- Klasyfikować folderów (klasyfikacja jest w Sesji 4).

**Rezultat:**
- 4 szablony w `meta/research/templates/`
- 3 podglądy w `meta/research/_examples/`
- Stub `AGENT_PROTOCOL.md`

**Pliki, które powstają:**
- `meta/research/templates/*` (nowe)
- `meta/research/_examples/*` (nowe)
- `meta/research/AGENT_PROTOCOL.md` (stub)

**Kryterium zakończenia:**
- Człowiek zaakceptował szablony lub wskazał poprawki.
- 3 podglądy nie zawierają wymyślonych statusów.

**Ryzyka:**
- Pole `source_of_status` w YAMLu może być za rygorystyczne i zablokować
  agenta nawet w jasnych przypadkach. Mitygacja: dopuścić wartość
  `"heuristic-by-folder-content (Sesja 1 audit)"` jako tymczasową, z
  oznaczeniem `level: unknown`.

---

### Sesja 3 — Utworzenie warstwy `meta/research/` + `meta/core/`

**Cel:** stworzyć wszystkie pliki z sekcji 5 i 5A jako stuby z poprawnym
frontmatterem i krzyżowymi linkami, ALE bez treści wymagającej
klasyfikacji folderów (to przyjdzie później). Plus utworzyć fizycznie
`research/_sandbox/` i `research/_archive/` z README opisującymi reguły.

**Wejście:**
- Niniejszy plan, sekcje 5 i 5A
- Szablony z Sesji 2

**Co agent ma zrobić:**
1. Utworzyć w `meta/research/`:
   - `RESEARCH_BUS.md` (pusta tablica + nagłówki)
   - `CANDIDATE_BRIDGES.md` (pusta tabela + reguły wpisu)
   - `CORE_CANDIDATES.md` (pusta tabela + sekcja
     "Required human approval gate")
   - `IMPACT_MATRIX.md` (kościec macierzy + import nagłówków
     z `DEPENDENCIES.md`)
   - `AGENT_PROTOCOL.md` (pełna wersja)
   - `FOLDER_STATUS_INDEX.md` (placeholder)
   - `MIGRATION_LOG.md` (pusty + sekcja "Reguły wpisu")
   - `GLOSSARY.md` (linki do `axioms/notacja/slownik_formalizmu.tex`
     + organizacyjne pojęcia)
2. Utworzyć w `meta/core/`:
   - `CORE_INVENTORY.md` (pusta mapa + reguły wpisu)
   - `CORE_HOTSPOTS.md` — **wstępnie wypełniony 8 CRITICAL + 12 HIGH** z
     [[meta/AUDYT_TGP_2026-05-01.md]] (transkrypcja, nie nowa fizyka)
   - `CORE_INTAKE.md` (pełny protokół promocji, statyczny)
   - `meta/core/intake/.gitkeep` (puste, na przyszłe wnioski)
3. Utworzyć fizycznie:
   - `research/_sandbox/README.md` z regułami "co to jest, czego tu nie
     trzymamy, kiedy folder z _sandbox awansuje do research/"
   - `research/_archive/README.md` z regułami z sekcji 3.4
     (kiedy coś tu trafia, jak różni się od `_archiwum/`)
4. Dodać każdemu z plików `related:` w frontmatterze (nawigacja
   krzyżowa).
5. Dodać do `meta/README.md` sekcję "Pliki dla `research/` workflow:"
   z linkami do `meta/research/*` i `meta/core/*`.

**Czego agentowi nie wolno:**
- Klasyfikować jakiegokolwiek folderu w `research/`.
- Edytować `INDEX.md` / głównego `README.md` / `core/` / `axioms/` /
  `main.tex` / `partial_proofs/` / `papers_external/`.
- Wpisywać do `CORE_INVENTORY.md` żadnego mappingu — plik jest pusty
  do Sesji 7.
- Modyfikować `CORE_HOTSPOTS.md` poza początkową transkrypcją z audytu
  (nie wolno samodzielnie dodawać nowych hotspotów ani zmieniać statusu
  "OPEN" na cokolwiek innego).

**Rezultat:**
- 8 plików w `meta/research/` + 3 pliki + folder `intake/` w `meta/core/`
- 2 nowe foldery fizyczne: `research/_sandbox/`, `research/_archive/`
  (każdy z `README.md` opisującym reguły)
- Aktualizacja `meta/README.md` (mała sekcja).

**Kryterium zakończenia:**
- Wszystkie pliki istnieją i mają zgodny frontmatter.
- `CORE_HOTSPOTS.md` zawiera 8 CRITICAL + 12 HIGH z audytu 2026-05-01,
  każdy ze źródłem.
- Człowiek zaakceptował treść `AGENT_PROTOCOL.md` **i** `CORE_INTAKE.md`
  — to są dwa twarde gate'y.

**Ryzyka:**
- `AGENT_PROTOCOL.md` może być za luźny i nie zapobiec przyszłym
  over-claimom. Mitygacja: w treści protokołu **zacytować** incydent
  `74394a8` jako case study, żeby agenci kolejnych sesji przeczytali
  i zrozumieli powód reguł.

---

### Sesja 4 — Klasyfikacja folderów `research/` (warstwa YAML)

**Cel:** dla każdego folderu w `research/` ustawić wstępne `tgp_status`
w jego `README.md` oraz wpisać go do `FOLDER_STATUS_INDEX.md`.

**Wejście:**
- Wszystkie foldery `research/`
- `meta/research/AUDIT_RESEARCH_S1.md`
- Szablony z Sesji 2
- `meta/AUDYT_TGP_2026-05-01.md` (źródło dla `core_compatibility`)
- `meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md` (zatrute foldery do
  ostrożnej klasyfikacji)

**Co agent ma zrobić:**
1. Dla każdego folderu (w batchach po ~10):
   - Odczytać istniejący README + program.md + Phase*_results.md (jeśli są).
   - Zaproponować `tgp_status` zgodnie z regułami sekcji 4.
   - **Z zastrzeżeniem reguły "anty-overclaim":** jeśli brak
     weryfikowalnego źródła statusu — `level: unknown`.
   - Wstawić frontmatter do `README.md` (jeśli folder nie ma README,
     stworzyć minimalny z szablonu).
2. Foldery z incydentu `74394a8` (`op-chi1-newton-constant-derivation`,
   `op-uv2-mtgp-absolute-scale`, `op-omega2-axion-coupling-lock`,
   `op-omega3-axion-decay-constant`) — automatycznie ustawić
   `level: unknown`, `core_compatibility: unknown`, `folder_status:
   needs-bridge` z notatką w `source_of_status` o audycie 2026-05-02.
3. Foldery z `closure_2026-04-26/` — oznaczyć **rodzic** jako
   `kind: closure-aggregator`, `level: mixed`. Każde podzamknięcie
   dostaje własny YAML w swoim README.
4. Wygenerować `FOLDER_STATUS_INDEX.md` skryptem (Python),
   deterministycznie sortowany alfabetycznie.

**Czego agentowi nie wolno:**
- Nadawać `level: L4` / `folder_status: core-promoted` bez weryfikowalnej
  referencji do `main.tex` / `core/` / `axioms/`.
- Modyfikować części README poza wstawieniem frontmattera (chyba że
  README nie istnieje — wtedy minimalny szablon).
- Przenosić folderów.

**Rezultat:**
- ~95 zaktualizowanych README.md w `research/*/`
- `meta/research/FOLDER_STATUS_INDEX.md` z pełną mapą.

**Kryterium zakończenia:**
- 100% folderów ma `tgp_status` w README.
- 0 folderów z `level: L4` bez `promoted_to_core` field.
- `FOLDER_STATUS_INDEX.md` zgadza się z YAMLami (skrypt-checker).

**Ryzyka:**
- Agent zacznie syntetyzować PASSe z treści README. Mitygacja: regułą
  protokołu — `source_of_status` musi cytować plik + sekcję, nie sam
  przymiotnik.
- 95 folderów = ryzyko zmęczenia kontekstu. Mitygacja: batche po 10 +
  checkpoint w `MIGRATION_LOG.md`.

---

### Sesja 5 — Dodanie `NEEDS.md` i `FINDINGS.md` do każdego folderu

**Cel:** wymusić eksportowalność wyników i jawność luk.

**Wejście:**
- Foldery z `tgp_status` ustawionym w Sesji 4
- Szablony

**Co agent ma zrobić:**
1. Dla każdego folderu, którego `has_findings_file: true` jeszcze nie
   spełnione — utworzyć `FINDINGS.md` na podstawie: README, program.md,
   Phase*_results.md, plików `.txt` z wynikami. **Tylko cytaty** —
   żadnych nowych derivation.
2. Dla każdego folderu — utworzyć `NEEDS.md` na podstawie: sekcji
   "Open issues" / "Phase X pending" / "TODO" / "BLOCKED" z istniejących
   plików, plus pól `open_bridges` z YAML.
3. Po utworzeniu — zaktualizować flagi `has_needs_file`, `has_findings_file`
   w YAMLu README.
4. Zsynchronizować z `RESEARCH_BUS.md`: każdy nowy `FINDINGS.md`
   generuje wpis broadcast.

**Czego agentowi nie wolno:**
- Wymyślać wyników, których nie ma w plikach folderu.
- Wpisywać do `NEEDS.md` luk, których nikt nigdzie nie zidentyfikował
  (każdy item musi mieć `source:` cytujący istniejący plik).
- Dotykać core.

**Rezultat:**
- ~95 plików `FINDINGS.md` (po jednym na folder)
- ~95 plików `NEEDS.md`
- Aktualizacja YAMLów (flagi)
- Zaktualizowany `RESEARCH_BUS.md`

**Kryterium zakończenia:**
- 100% pokrycia. Skrypt-checker potwierdza spójność flag.

**Ryzyka:**
- Agent stworzy puste pliki "zaliczeniowe". Mitygacja: skrypt-checker
  weryfikuje minimalną długość + obecność co najmniej 1 cytatu w każdym
  pliku, w przeciwnym razie ustawia flagę `has_findings_file: false`
  z notatką "EMPTY".

---

### Sesja 6 — Research Bus i Candidate Bridges

**Cel:** pierwsze realne broadcasty międzyfolderowe i mapowanie luk.

**Wejście:**
- Wszystkie `FINDINGS.md`, `NEEDS.md` z Sesji 5
- `IMPACT_MATRIX.md` (skeleton z Sesji 3)
- `DEPENDENCIES.md`, `DEPENDENCIES_REVERSE.md`

**Co agent ma zrobić:**
1. Dla każdego `NEEDS.md` — sprawdzić, czy któryś `FINDINGS.md` w innym
   folderze mówi o tej luce. Jeśli tak — dopisać do `CANDIDATE_BRIDGES.md`
   wpis "BR-NN: folder A potrzebuje X; folder B ma X (refer)".
2. Każdy match ⇒ wpis w `RESEARCH_BUS.md` typu "PROPOSED BRIDGE".
3. Zaktualizować `IMPACT_MATRIX.md` z relacjami (nadpisać `depends_on` /
   `impacts` w YAMLach, jeśli ujawniono nowe powiązania).
4. Wypisać **3 najbardziej obiecujące mosty** do uwagi człowieka
   (sekcja "Top bridges for human review" w `CANDIDATE_BRIDGES.md`).

**Czego agentowi nie wolno:**
- Wykonywać samego mostu (to jest nowy research, nie jest częścią
  organizacji).
- Modyfikować plików innych folderów poza dodawaniem wpisów do
  `RESEARCH_BUS.md`.

**Rezultat:**
- Wypełniony `CANDIDATE_BRIDGES.md` z N kandydatami
- Wypełniony `RESEARCH_BUS.md` z broadcastami
- Pełna `IMPACT_MATRIX.md`

**Kryterium zakończenia:**
- Każdy `NEEDS.md` ma w `CANDIDATE_BRIDGES.md` przynajmniej wpis
  "scanned, no match" lub konkretnego kandydata.
- Człowiek zaakceptował top 3 mosty (gate dla Sesji 7).

**Ryzyka:**
- Agent dopasuje na siłę luźne podobieństwa terminów (np. dwa różne
  "α"). Mitygacja: regułą protokołu — match musi cytować dokładny lemma
  / formułę / numeryczną wartość, nie sam keyword.

---

### Sesja 7 — Klasyfikacja core-ready / core-promoted / needs-migration

**Cel:** zidentyfikować foldery, które są (a) gotowe do awansu do core,
(b) już są de facto w core, (c) wymagają migracji po pivotcie M9.1''.

**Wejście:**
- `FOLDER_STATUS_INDEX.md` z Sesji 4
- `CANDIDATE_BRIDGES.md` z Sesji 6
- `meta/AUDYT_TGP_2026-05-01.md` (lista 8 CRITICAL)
- `meta/core/CORE_HOTSPOTS.md` (Sesja 3, transkrypcja audytu)
- `meta/core/CORE_INTAKE.md` (Sesja 3, protokół promocji)
- `meta/core/CORE_INVENTORY.md` (pusta przed sesją; po sesji wypełniona
  dla folderów L4)
- `main.tex`, `core/sek*.tex`, `axioms/*`, `partial_proofs/*` —
  **read-only**, żeby sprawdzić, co już jest w core.

**Co agent ma zrobić:**
1. Dla każdego folderu z `level: L3` — sprawdzić w `core/`, czy jego
   wynik **już** jest w sekcji `sek0X` lub w `axioms/` lub w
   `partial_proofs/`. Jeśli tak ⇒ podnieść do `level: L4`,
   `folder_status: core-promoted`, wypełnić `promoted_to_core` z
   konkretną referencją (`core/sek08c_metryka_z_substratu/m911.tex
   §3.2` itp.). Dopisać wiersz do `meta/core/CORE_INVENTORY.md`
   z linkiem zwrotnym `Promoted from`.
2. Dla pozostałych `level: L3` ⇒ `folder_status: core-ready`. Wpis w
   `CORE_CANDIDATES.md` z polem **"required human decision"**. Agent
   **może proponować** awans do L4 przez utworzenie wniosku
   `meta/core/intake/INTAKE_<folder>_2026-05-XX.md` zgodnie z protokołem
   `CORE_INTAKE.md` — wniosek to artefakt do oceny człowieka, NIE jest
   wykonaniem promocji.
3. Foldery z `core_compatibility: stale | broken` (typowo M9.x post-pivot,
   sek08c-related) ⇒ `folder_status: needs-migration`. Wpis w
   `IMPACT_MATRIX.md` z listą wymaganych zmian (cytując
   `meta/core/CORE_HOTSPOTS.md §HN` zamiast bezpośrednio audytu, żeby
   ślad był live, nie snapshot).
4. Wygenerować "shortlist" 5 najbardziej dojrzałych core-ready do
   ręcznego review (osobno od wniosków `INTAKE_*.md`).
5. Dla każdego rozwiązanego hotspota w `CORE_HOTSPOTS.md` (jeśli
   znajdzie się folder, który już go adresuje) — proponować zmianę
   statusu hotspota z `OPEN` na `READY-FOR-INTAKE` (decyzja człowieka).

**Czego agentowi nie wolno:**
- Modyfikować `core/`, `axioms/`, `partial_proofs/`, `main.tex`,
  `tgp_letter.tex`, `tgp_companion.tex`, `arxiv_submission/`,
  `paper_lepton_masses/`, `paper_bh_shadow/`.
- Samodzielnie awansować folderu do `level: L4` bez konkretnej
  referencji do core (tylko gdy plik core **już** zawiera ten wynik —
  zgodnie z regułą "core-promoted = już w core, nie 'powinno być'").
- Wykonywać migracji "needs-migration" (to jest physics work, nie
  organizacja).
- Zmieniać statusu hotspota w `CORE_HOTSPOTS.md` na `RESOLVED`
  (to wymaga decyzji człowieka po faktycznej migracji core).

**Rezultat:**
- Zaktualizowane YAMLe w ~5–15 folderach (zależnie od stanu)
- Wypełniony `CORE_CANDIDATES.md`
- Wypełniony `CORE_INVENTORY.md` (dla folderów L4)
- Lista needs-migration w `IMPACT_MATRIX.md`
- 0–N wniosków `meta/core/intake/INTAKE_*.md` do oceny człowieka

**Kryterium zakończenia:**
- 0 folderów z niespójnym `promoted_to_core`.
- Każdy wpis w `CORE_INVENTORY.md` ma weryfikowalny link do core.
- Człowiek zaakceptował shortlist core-ready.

**Ryzyka:**
- Agent źle przeczyta `main.tex` i błędnie awansuje folder. Mitygacja:
  wymaganie zacytowania konkretnej linii / sekcji `main.tex` w
  `promoted_to_core`.
- Agent pominie folder, bo `core/` używa innej terminologii niż
  `research/`. Mitygacja: użycie `GLOSSARY.md` jako mapy synonimów.

---

### Sesja 8 — Regression / impact review (read-only kontrola jakości)

**Cel:** sprawdzić spójność całej warstwy meta + odporność na
over-claiming.

**Wejście:**
- Cała warstwa `meta/research/`
- Wszystkie YAMLe w folderach badawczych
- `INDEX.md`, `DEPENDENCIES.md`, `DEPENDENCIES_REVERSE.md`

**Co agent ma zrobić:**
1. Skrypt-checker:
   - Każdy folder ma README + YAML + NEEDS + FINDINGS.
   - Każde `level: L4` ma `promoted_to_core` z weryfikowalną referencją.
   - Każdy wpis w `RESEARCH_BUS.md` ma swoje źródło w jakimś `FINDINGS.md`.
   - Każdy bridge w `CANDIDATE_BRIDGES.md` ma source/target istniejące
     w repo.
   - Każda `depends_on` / `impacts` referencja prowadzi do realnego
     folderu.
   - Każdy `last_reviewed_against_core` jest ≤ dziś.
2. Audyt anty-overclaim:
   - Wyszukać w wszystkich `FINDINGS.md` i `NEEDS.md` fraz "FULL
     CONVERGENCE", "100%", "DERIVED" bez `source_of_status`.
   - Wyszukać `source_of_status` z formą "heuristic" tam, gdzie YAML
     deklaruje `level: L3` lub wyżej (sprzeczność).
3. Raport `meta/research/REGRESSION_S8.md` z listą problemów + rekomendacji.

**Czego agentowi nie wolno:**
- Naprawiać znalezionych błędów automatycznie. Tylko raport.
- Edytować plików innych niż `meta/research/REGRESSION_S8.md`.

**Rezultat:**
- `meta/research/REGRESSION_S8.md`
- Lista poprawek dla Sesji 4–7 (jeśli potrzebne)

**Kryterium zakończenia:**
- Raport istnieje i ma jawną sekcję "0 critical findings" lub listę
  z konkretnymi PR-like poprawkami.

**Ryzyka:**
- Skrypt-checker będzie miał false positives (np. uznać legalny `L3`
  za podejrzany). Mitygacja: każdy false positive musi być w raporcie
  oznaczony i można go ręcznie zaakceptować w Sesji 4-7 retrospektywnie.

---

### Sesja 8.5 — Selektywna archiwizacja do `research/_archive/`

**Cel:** wykonać fizyczne `git mv` **tylko** dla folderów, które
spełniają rygorystyczne kryteria archiwizacji z sekcji 3.4. To jest
jedyne fizyczne przenoszenie w wariancie minimalnym.

**Wejście:**
- `FOLDER_STATUS_INDEX.md` (po Sesji 8)
- `meta/research/REGRESSION_S8.md`
- Ręczna shortlista kandydatów do archiwizacji od człowieka **lub**
  agentowa propozycja (każdy item z 2 źródłami zgodnie z 3.4)

**Co agent ma zrobić:**
1. Wygenerować propozycję `meta/research/ARCHIVE_PROPOSAL_S8_5.md`
   z listą folderów spełniających kryteria 3.4 (każdy z dwoma cytatami
   źródłowymi: `WITHDRAWN`/`FALSIFIED`/`SUPERSEDED`).
2. **Czekać na akceptację człowieka** dla każdej pozycji osobno (gate).
3. Dla każdego zaakceptowanego folderu:
   a. `git mv research/<folder> research/_archive/<folder>`
   b. Wyszukać wszystkie referencje (`[[research/<folder>/...]]`,
      `\input{research/<folder>/...}`, `\ref{research/<folder>/...}`,
      `import` w `*.py`) w całym repo.
   c. Zaktualizować je do `research/_archive/<folder>/`.
   d. W `README.md` przeniesionego folderu zmienić
      `folder_status: archive`, dodać `archived_date`, `archived_reason`,
      `archived_by_human_decision: <data>`.
   e. Wpis w `meta/research/MIGRATION_LOG.md` (data, ID akceptacji,
      lista zaktualizowanych ścieżek, hash commita).
4. Po wszystkich `mv` — uruchomić `tooling/build_deps_graph.py` (jeśli
   istnieje) i porównać orphans z poprzednim runem. Orphans dla
   przenoszonych folderów muszą być 0.
5. Skompilować `main.tex` (sanity check) — jeśli nie kompiluje się,
   przerwać i `git revert`.

**Czego agentowi nie wolno:**
- Wykonywać `mv` bez explicite akceptacji człowieka dla każdej pozycji.
- Przenosić folderów spełniających któryś z warunków sekcji 3.4
  ("zasadniczo do `_archive/` NIE trafiają …"):
  - foldery `core-promoted` (zostają in place jako historia),
  - `closure_*` (closure-aggregator, in place),
  - `op6/`, `op7/`, `op1-op2-op4/` (legacy z gęstym graphem linków),
  - `external_review_*` (wzajemne audyty, in place).
- Łączyć wielu `mv` w jeden commit. Każdy folder = osobny commit.
- Dotykać `core/`, `axioms/`, `partial_proofs/`, `papers_external/`,
  `main.tex` poza poprawą ścieżek `\input`/`\ref`/wikilink dla
  przenoszonych folderów.

**Rezultat:**
- N folderów fizycznie w `research/_archive/<folder>/`
- 0 broken links
- `main.tex` kompiluje się
- N wpisów w `MIGRATION_LOG.md`
- Aktualizacja `FOLDER_STATUS_INDEX.md`

**Kryterium zakończenia:**
- `pdflatex main && pdflatex main` przechodzi.
- Re-run `tooling/build_deps_graph.py` — orphans = 0.
- Każdy przeniesiony folder ma wpis w `MIGRATION_LOG.md` z datą
  akceptacji człowieka.

**Ryzyka:**
- Linki w plikach `.txt` z wynikami skryptów — zwykle tylko ścieżki
  string-owe, ale `git mv` ich nie patchuje. Mitygacja: scan `*.txt`
  + `*.py` na hardcoded paths przed `mv`.
- Human zaakceptuje archiwizację folderu, który ma żywego dependenta.
  Mitygacja: pre-flight check `IMPACT_MATRIX.md` — jeśli folder X jest
  `depends_on:` jakiegoś żywego folderu, ostrzec i wymagać dodatkowej
  akceptacji.

---

### Sesja 9 (opcjonalna, prawdopodobnie SKIP w wariancie minimalnym) — Fizyczna migracja folderów

**Cel:** **TYLKO jeśli** człowiek wybrał wariant rozszerzony w Sesji 1
i zaakceptował shortlist. Wykonuje fizyczne `mv`.

**Wejście:**
- Akceptacja człowieka.
- Lista folderów do przeniesienia (z `FOLDER_STATUS_INDEX.md`).
- Lista linków do naprawy (z `INDEX.md`, `DEPENDENCIES.md`,
  `program.md` w przenoszonych folderach, `meta/research/*`).

**Co agent ma zrobić:**
1. Dla każdego folderu z listy:
   a. `git mv research/<folder> research/<podfolder>/<folder>`
   b. Wyszukać wszystkie `[[research/<folder>/...]]` i
      `\input{research/<folder>/...}` i `\ref{research/<folder>/...}`
      w całym repo.
   c. Zaktualizować je do nowej ścieżki.
   d. Wpisać do `MIGRATION_LOG.md`.
2. Po wszystkich `mv` — uruchomić istniejący `tooling/build_deps_graph.py`
   (jeśli istnieje), wygenerować nowe `DEPENDENCIES.md` i porównać z
   poprzednim — orphans muszą być 0.
3. Skompilować `main.tex` (sanity check) — jeśli się nie kompiluje,
   przerwać i zwrócić raport.

**Czego agentowi nie wolno:**
- Wykonywać `mv` bez wcześniejszego pełnego patcha linków w jednej
  atomic operation (per folder).
- Robić squash commitów. Każdy folder = osobny commit z czytelną
  wiadomością.
- Modyfikować zawartości plików `core/` / `axioms/` poza poprawą
  ścieżek `\input`/`\ref`/wikilink.

**Rezultat:**
- Fizycznie zreorganizowane `research/`
- 0 broken links
- `main.tex` kompiluje się
- Pełny `MIGRATION_LOG.md`

**Kryterium zakończenia:**
- `pdflatex main && pdflatex main` przechodzi.
- `tooling/build_deps_graph.py` (lub odpowiednik) zwraca 0 orphans.

**Ryzyka:**
- Linki w plikach `.txt`, `.py` (np. `program.md` referuje
  `../op-foo/Phase3_results.md`) zostaną złamane. Mitygacja: scan na
  plikach `*.md`, `*.tex`, `*.py` (raw string w komentarzach też).

---

## 7. Ryzyka globalne i safeguardy

| Ryzyko | Mitygacja |
|---|---|
| Agent over-claim (powtórka `74394a8`) | `AGENT_PROTOCOL.md` + reguła `source_of_status` + audit Sesja 8 |
| Złamanie linków w `INDEX.md` | Domyślnie wariant minimalny (bez `mv`); Sesja 9 opcjonalna |
| Konflikt z `_archiwum/` | Pierwsze pytanie do człowieka w Sesji 1 |
| Niezgodność `./research/` (vault root) z `TGP/TGP_v1/research/` | Decyzja w Sesji 1, plan działa głównie na canonical |
| Ekspansja kontekstu (95 folderów) | Batche po 10 + checkpoint w `MIGRATION_LOG.md` |
| Niezamierzone dotknięcie core | Twardy ban na edycję core/axioms/main.tex/papers w sesjach 1–8 |
| `level: L4` bez podstawy | Wymaganie `promoted_to_core` z linią/sekcją `main.tex` |
| `FINDINGS.md` jako pusta zaliczeniówka | Skrypt-checker minimal-length + cytat |
| `CANDIDATE_BRIDGES` mache słów | Reguła: cytat z formułą / wartością, nie sam keyword |

---

## 8. Proponowana migracja (wariant minimalny — wybrany)

### 8.1 Setup w Sesji 3 (mechaniczny, bez decyzji)

```bash
# Sandbox dla agentów (nowy)
mkdir -p TGP/TGP_v1/research/_sandbox
# (README.md tworzony przez Sesja 3 z regułami)

# Research-level archiwum (nowy)
mkdir -p TGP/TGP_v1/research/_archive
# (README.md tworzony przez Sesja 3 z regułami sekcji 3.4)

# Warstwa meta dla badań
mkdir -p TGP/TGP_v1/meta/research/templates
mkdir -p TGP/TGP_v1/meta/research/_examples

# Warstwa meta dla rdzenia
mkdir -p TGP/TGP_v1/meta/core/intake
```

### 8.2 Selektywna archiwizacja w Sesji 8.5 (z akceptacją człowieka)

```bash
# Tylko dla folderów spełniających kryteria sekcji 3.4
# i z explicite akceptacją per folder:
git mv TGP/TGP_v1/research/<obsolete-folder> \
       TGP/TGP_v1/research/_archive/<obsolete-folder>

# + patch linków + osobny commit per folder
# (zob. Sesja 8.5)
```

### 8.3 Co NIGDY nie powinno być przenoszone bez świadomej decyzji

- `research/op6/`, `research/op7/`, `research/op1-op2-op4/` — legacy z
  bardzo gęstym graph linkowym.
- `research/closure_2026-04-26/` — closure-aggregator z 4 podzamknięciami,
  gęsto linkowany z `KNOWN_ISSUES.md`.
- `research/external_review_2026-04-25/` — review historyczne,
  bezpieczniej zostawić in place.
- ~~`research/TGP_STATUS_*.md` / `REDIRECT_PROGRAM_*.md` /
  `NEW_DIRECTIONS_*.md`~~ — **przeniesione 2026-05-03** do
  `_archiwum/research_program_docs_2026-04/` (decyzja Q3, sekcja 9.2).
- Foldery z `folder_status: core-promoted` (są in place jako historia
  i uzasadnienie tego, co poszło do core).

### 8.4 Co NIE jest częścią wariantu minimalnego

- `research/active/`, `research/core-ready/`, `research/core-promoted/`,
  `research/needs-bridge/`, `research/needs-migration/`,
  `research/phenomenology/`, `research/numerical/` — **nie tworzymy**.
  Te statusy żyją w YAMLu i `FOLDER_STATUS_INDEX.md`.
- Sesja 9 jest **opcjonalna** i domyślnie pomijana.

---

## 9. Decyzje człowieka (zatwierdzone 2026-05-02)

| # | Pytanie | Decyzja |
|---|---|---|
| 1 | Który `research/` jest kanoniczny? | `TGP/TGP_v1/research/` (canonical). `./research/` (vault root) jest **out-of-scope** tego planu — Sesja 1 raportuje stan, ale nie modyfikuje. |
| 2 | Wariant minimalny czy rozszerzony? | **Minimalny.** Bez fizycznych nadfolderów statusowych. Status w YAMLu. |
| 3 | Drugie archiwum w `research/`? | **Tak**, `research/_archive/` z rygorystycznymi regułami (sekcja 3.4). Różne od `_archiwum/` (TGP_v1 level). |
| 4 | Lokalizacja warstwy meta? | `meta/research/` **+ `meta/core/`** (3 pliki, sekcja 5A). |
| 5 | Limit batcha w Sesji 4? | **10 folderów na batch** (potwierdzone). |
| 6 | Czy agent może proponować awans do L4? | **Tak.** Wniosek `meta/core/intake/INTAKE_*.md` z weryfikowalnym cytatem. Decyzja człowieka jako gate. Stan obecny = "śmietnik do posprzątania", więc agent ma aktywnie sugerować promocje. |

**Wniosek z decyzji 6:** Sesja 7 jest istotnym kierunkiem clean-upu —
agent celowo szuka folderów, które są de facto już w core, oraz tych,
które są L3 i tylko czekają na promocję. Nie ma neutralnej "obserwacji" —
plan dostarcza shortlistę i wnioski.

### 9.1 Pytania post-Sesja 1 (zatwierdzone 2026-05-03)

Po raporcie [[meta/research/AUDIT_RESEARCH_S1.md]]:

| # | Pytanie | Decyzja |
|---|---|---|
| Q1 | Co z vault-rootowym `./research/` (11 stub-folderów, 14 KB, content = stderr z nieudanych runów)? | **Skasować.** Wykonane 2026-05-03 (`shutil.rmtree`). Vault root nie zawiera już duplikatu. |
| Q2 | Foldery z subfolderami (np. `closure_2026-04-26/sigma_ab_pathB/`) — pełna warstwa per dziecko, czy tylko rodzic? | **Tylko rodzic.** Dzieci są lekkie — pojawiają się w `subfolder_names` rodzica, ale nie dostają osobnych YAML/NEEDS/FINDINGS. Wyjątek: jeśli dziecko jest realnie samodzielnym workstreamem (decyduje człowiek per case), agent w Sesji 4 zostawia notatkę w README rodzica i czeka. |
| Q3 | Top-level pliki w `research/` (`TGP_STATUS_*`, `REDIRECT_PROGRAM_*`, `NEW_DIRECTIONS_*`) — co z nimi? | **Przeniesione do `_archiwum/research_program_docs_2026-04/`** (2026-05-03). Wikilinks `[[FILENAME.md]]` (8 referencji w research) działają nadal — Obsidian rozwiązuje po nazwie pliku. Eksplicytne `[[research/...]]` paths zaktualizowane (tylko niniejszy plan). |
| Q4 | `FINDINGS.md`/`NEEDS.md` — nadpisywać istniejące, czy tylko uzupełniać? | **Tylko brakujące.** Domyślne zachowanie: agent nie nadpisuje istniejących plików. Jeśli plik istnieje, agent dopisuje pole `pre_existing: true` w YAML rodzica. (W Sesji 1 audyt wykazał 0/86 — ale reguła zostaje na przyszłość.) |

### 9.2 Aktualizacja stanu repo (2026-05-03)

- ✅ `./research/` (vault root) — **usunięte** (11 stub-folderów po `shutil.rmtree`).
- ✅ `research/TGP_STATUS_2026-04-18.md` → `_archiwum/research_program_docs_2026-04/`
- ✅ `research/TGP_STATUS_2026-04-19.md` → `_archiwum/research_program_docs_2026-04/`
- ✅ `research/REDIRECT_PROGRAM_2026-04-19.md` → `_archiwum/research_program_docs_2026-04/`
- ✅ `research/NEW_DIRECTIONS_2026-04-20.md` → `_archiwum/research_program_docs_2026-04/`
- Top-level pliki w `research/` po przeniesieniu: tylko grafy
  (`graph_*.gexf`, `graph_*.png`, `tgp_dependency_graph.py`).

### 9.3 Sesja 3.5 — Hotspot Reality Check (insert post-Sesji 3, 2026-05-03)

**Trigger:** człowiek 2026-05-03: "sporo z tych rzeczy jest już pozamykane,
zastanawiam się czy warto się zatrzymać i zrobić audyt tych 8 krytycznych
problemów przed dalszą sesją".

**Co odkryto:** `AUDYT_TGP_2026-05-01.md` zawiera **31 paragrafów § A–§ AB**.
Sesja 3 (transkrypcja do `CORE_HOTSPOTS.md`) czytała tylko § A (CRITICAL)
i § B (HIGH) z sekcji diagnozy — ignorując § J–§ AB ("Naprawa…" sekcje
zamykające 43/43 = 100% items).

**Wynik:** Mój `CORE_HOTSPOTS.md` z Sesji 3 fałszywie oznaczył **15 z 20
hotspotów jako OPEN**. Faktycznie wszystkie 20 są zamknięte
w samym audycie (różne typy closure):

| Closure type | Liczba hotspotów A+B |
|--------------|---------------------:|
| RESOLVED-AUDIT-FULL | 11 |
| RESOLVED-AUDIT-ANNOTATION | 6 |
| RESOLVED-AUDIT-CONVENTION | 1 (A3) |
| RESOLVED-AUDIT-ARCHITECTURAL | 1 (A4) |
| RESOLVED-AUDIT-FULL (przez wcześniejszy ω.3 cycle) | 1 (A7) |
| **Razem** | **20** |

**Naprawione pliki post-Sesja 3.5:**

- ✅ Nowy raport: `meta/research/HOTSPOT_AUDIT_S3_5.md` (per-hotspot werdykt + closure-level taxonomy)
- ✅ Przepisany: `meta/core/CORE_HOTSPOTS.md` (z `15 OPEN + 4 READY-FOR-INTAKE + 1 RESOLVED` → `0 OPEN + 20 RESOLVED-AUDIT-*`)
- ✅ Zaktualizowany: `meta/research/AGENT_PROTOCOL.md` §3 — dopisany **case study #2** (under-claim Sesji 3) jako symetryczny do case study #1 (over-claim `74394a8`); reguła 3.0 wymaga przeczytania całego dokumentu źródłowego, włącznie z final closure summary
- ✅ Zaktualizowany: `meta/research/IMPACT_MATRIX.md` §3 — krawędzie z hotspotów oznaczone jako `HIST` (historical paths of closure), `CROSS-VAL` (alternative formal tracks), `NEW` (structural contributions poza 43-list)
- ✅ Zaktualizowany: `meta/research/CORE_CANDIDATES.md` — reframing kandydatów `closure_2026-04-26/*` jako cross-validation (audyt zamknął te same hotspoty niezależnie przez `op-newton-momentum/B6/B8/B9_*.py`)

**Lesson learned (do AGENT_PROTOCOL §3, case #2):** Audyty TGP zawierają
zarówno diagnozę jak i naprawę w jednym pliku. Skanowanie nazw paragrafów
to pierwszy krok, nie ostatni — trzeba wyłuskać `final tally` /
`closure aggregate` (w `AUDYT_TGP_2026-05-01.md` sekcja § AB.5 daje
TOTAL: 43/43 = 100%).

**Implikacje dla Sesji 4:**

- `needs-migration` jako kategoria może być **niemal pusta** — większość
  M9.x folderów ma już audit-aware markery (`B6-CLOSED`, `B8-CLOSED`,
  `B9-CLOSED`). Folder zostaje `needs-migration` **tylko** gdy ma
  dowód, że żywa luka istnieje **poza** sekcjami § A–§ AB audytu.
- `core_compatibility: stale | broken` powinno być rzadkie (większość:
  `current` lub `partial`).
- `op-tau3-substrate-clock-acceleration` ma bogaty zestaw audit-aware
  patches (multiplicative m_e_eff, B7 Greens, B7-v2 numerical, Phase 4
  Adams positivity) — przykład folderu, który zaadresował 3 hotspoty
  (H-A5, H-B1, H-B7) i jest w stanie `current` po patches.

---

## 10. Konwencje dla agentów (skrót — pełna wersja w `AGENT_PROTOCOL.md` po Sesji 3)

- **Nigdy** nie edytuj `core/`, `axioms/`, `partial_proofs/`,
  `papers_external/`, `main.tex`, `tgp_letter.tex`, `tgp_companion.tex`
  bez explicite rozkazu człowieka.
- **Nigdy** nie wymyślaj PASSów, „FULL CONVERGENCE" ani „DERIVED"
  bez `source_of_status` cytującego konkretny plik + linijka/sekcja.
- **Zawsze** zapisuj do `RESEARCH_BUS.md` po znalezieniu wyniku, który
  inny folder mógłby wykorzystać.
- **Zawsze** trzymaj `last_reviewed_against_core` świeży po dotknięciu
  YAMLa.
- **Zawsze** używaj relatywnych ścieżek od `TGP/TGP_v1/` w wikilinkach
  i frontmatterze.
- **W przypadku konfliktu** ze stanem core lub niepewności, ustaw
  `folder_status: needs-migration`, dopisz do `MIGRATION_LOG.md` i
  zatrzymaj sesję — czekaj na człowieka.

---

## 11. Co plan **świadomie pomija**

- Nie projektuje formatu skryptów weryfikacyjnych (to jest physics
  research; plan jest organizacyjny).
- Nie definiuje, co konkretnie znaczy "PASS" w Phase results — szanuje
  istniejące konwencje (`PASS/FAIL`, `LOCKED`, `STRUCTURAL`).
- Nie zmienia systemu wersji predykcji (`PREDICTIONS_REGISTRY.md` jest
  poza scope tego planu, ale jego kompatybilność jest zachowana).
- Nie podejmuje decyzji o `./research/` (vault root) vs
  `TGP/TGP_v1/research/` — to pytanie do człowieka.
- Nie tworzy nowego archiwum, dopóki człowiek nie zdecyduje o relacji
  z `_archiwum/`.

---

## 12. Pliki, które ten plan **przewiduje** stworzyć (sumarycznie)

```
TGP/TGP_v1/meta/research/
├── AGENT_PROTOCOL.md                    (Sesja 3)
├── ARCHIVE_PROPOSAL_S8_5.md             (Sesja 8.5, jeśli są kandydaci)
├── AUDIT_RESEARCH_S1.md                 (Sesja 1)
├── CANDIDATE_BRIDGES.md                 (Sesja 3 stub, Sesja 6 wypełnienie)
├── CORE_CANDIDATES.md                   (Sesja 3 stub, Sesja 7 wypełnienie)
├── FOLDER_STATUS_INDEX.md               (Sesja 3 stub, Sesja 4 wypełnienie)
├── GLOSSARY.md                          (Sesja 3)
├── IMPACT_MATRIX.md                     (Sesja 3 stub, Sesja 6+7 wypełnienie)
├── MIGRATION_LOG.md                     (Sesja 3 stub, Sesja 8.5 wpisy)
├── REGRESSION_S8.md                     (Sesja 8)
├── RESEARCH_BUS.md                      (Sesja 3 stub, Sesja 5+6 wpisy)
├── _examples/                           (Sesja 2)
│   ├── op-uv-as-ngfp_README.preview.md
│   ├── qm_born_rule_README.preview.md
│   └── closure_2026-04-26_README.preview.md
└── templates/                           (Sesja 2)
    ├── FINDINGS.template.md
    ├── NEEDS.template.md
    ├── README.template.md
    └── STATUS_BLOCK.yaml

TGP/TGP_v1/meta/core/                    (NOWA warstwa, Sesja 3)
├── CORE_INVENTORY.md                    (Sesja 3 stub, Sesja 7 wypełnienie)
├── CORE_HOTSPOTS.md                     (Sesja 3, transkrypcja audytu 2026-05-01)
├── CORE_INTAKE.md                       (Sesja 3, statyczny protokół)
└── intake/                              (Sesja 3, .gitkeep)
    └── INTAKE_<folder>_<date>.md        (Sesja 7+, per wniosek)

TGP/TGP_v1/research/
├── _sandbox/                            (Sesja 3, fizyczny folder)
│   └── README.md                        (reguły piaskownicy)
├── _archive/                            (Sesja 3, fizyczny folder)
│   ├── README.md                        (reguły archiwizacji 3.4)
│   └── <archived-folder>/               (Sesja 8.5, per akceptacja)
└── (~95 folderów, każdy z dopisanym `tgp_status:` w README + NEEDS.md + FINDINGS.md)

TGP/TGP_v1/meta/README.md                (mała aktualizacja, Sesja 3)
```

### 12.1 Skrót: ile plików powstaje per sesja

| Sesja | Nowe pliki | Edytowane pliki | Fizyczne `mv` |
|---|---|---|---|
| 1   | 1 (AUDIT_RESEARCH_S1.md) | 0 | 0 |
| 2   | ~7 (4 templates + 3 examples) | 0 | 0 |
| 3   | ~14 (8 meta/research + 3 meta/core + 2 sandbox/archive READMEs + meta/README update) | 1 (meta/README.md) | 0 |
| 4   | ~95 (NEEDS+FINDINGS są w Sesji 5; tu tylko aktualizacja README YAMLów) | ~95 (READMEs) | 0 |
| 5   | ~190 (NEEDS+FINDINGS, 2 per folder) | ~95 (READMEs flag updates) | 0 |
| 6   | 0 | 4 (RESEARCH_BUS, CANDIDATE_BRIDGES, IMPACT_MATRIX, YAMLs) | 0 |
| 7   | 0–N (INTAKE_*.md per propozycja) | ~5–15 (YAMLs w core-promoted) + CORE_INVENTORY + CORE_CANDIDATES | 0 |
| 8   | 1 (REGRESSION_S8.md) | 0 | 0 |
| 8.5 | 1 (ARCHIVE_PROPOSAL_S8_5.md) | per-folder linki + MIGRATION_LOG | 0–N (z akceptacją) |
| 9   | (opcjonalna, prawdopodobnie SKIP) | — | — |

---

## 13. Kolejność wykonania (TL;DR)

1. ✅ **Człowiek zdecydował** (sekcja 9, 6 pytań — 2026-05-02).
2. **Sesja 1** → audyt + raport `AUDIT_RESEARCH_S1.md`.
3. **Sesja 2** → szablony + 3 podglądy.
4. **Sesja 3** → warstwa `meta/research/` + **`meta/core/`** + fizyczne
   `research/_sandbox/` i `research/_archive/`.
5. **Sesja 4** → klasyfikacja YAML w ~95 folderach (batch po 10).
6. **Sesja 5** → `NEEDS.md` + `FINDINGS.md` w ~95 folderach.
7. **Sesja 6** → Research Bus + Candidate Bridges.
8. **Sesja 7** → core-ready / core-promoted / needs-migration; agent
   **proponuje** awanse przez `INTAKE_*.md` (decyzja 6).
9. **Sesja 8** → regression + audyt anty-overclaim.
10. **Sesja 8.5** → selektywna archiwizacja do `research/_archive/`
    (jedyne fizyczne `mv` w wariancie minimalnym, każde z akceptacją
    człowieka).
11. **Sesja 9** (opcjonalna, domyślnie SKIP) → masowe `mv`, gdyby
    człowiek zmienił zdanie i chciał wariant rozszerzony.

Każda sesja kończy się **gate'em człowieka** zgodnie z polem
"Kryterium zakończenia". Brak gate'a = brak przejścia do następnej sesji.

### 13.1 Gate'y krytyczne (twarde)

- Po Sesji 3: akceptacja `AGENT_PROTOCOL.md` **i** `CORE_INTAKE.md`.
- Po Sesji 4: 0 folderów z `level: L4` bez `promoted_to_core`.
- Po Sesji 7: akceptacja shortlisty core-ready **i** każdego
  `INTAKE_*.md` osobno (jeśli powstały).
- Po Sesji 8: `REGRESSION_S8.md` ma "0 critical findings" lub powrót
  do Sesji 4–7 z poprawkami.
- Sesja 8.5: per-folder akceptacja człowieka, każdy `mv` to osobny
  commit z zaktualizowanymi linkami.
