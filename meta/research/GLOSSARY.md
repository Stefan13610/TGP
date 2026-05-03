---
title: "GLOSSARY — terminy organizacyjne workflow research/ + meta/"
date: 2026-05-03
type: glossary
status: STATIC v1.0 (Sesja 3)
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[axioms/notacja/slownik_formalizmu.tex]]"
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/templates/STATUS_BLOCK.yaml]]"
tags:
  - glossary
  - terminology
  - organizational
---

# GLOSSARY — terminy organizacyjne

> **Cel pliku:** ujednolicić **organizacyjną** terminologię workflow.
> **NIE powtarza** słownika fizycznego (formuły, oznaczenia zmiennych,
> formulacje TGP) — to żyje w
> [[axioms/notacja/slownik_formalizmu.tex]].
>
> Zakres: terminy z `tgp_status` YAMLa, statusy hotspotów, role plików
> meta, terminy bridge / bus / intake.

---

## A. Słownik fizyczny TGP — delegacja

Pojęcia fizyczne (Φ, ψ, Φ₀, K(φ), V(φ), σ_ab, M9.1'', PPN, AS NGFP,
itd.) — patrz:

- **`axioms/notacja/slownik_formalizmu.tex`** — kanoniczny słownik
  formalizmu (Form. A / Sub / B, K, α, ghost-point, ODE solitonu)
- **`TGP_FOUNDATIONS.md`** — fundamenty teorii
- **`core/sek00_summary/`** — executive summary
- **`PREDICTIONS_REGISTRY.md`** — rejestr predykcji

**Reguła:** nowy termin fizyczny pojawia się w
`slownik_formalizmu.tex`, NIE w tym pliku.

---

## B. Statusy folderu (`tgp_status.folder_status`)

| Termin | Definicja | Kryterium twarde |
|--------|-----------|------------------|
| **active** | Aktywny research, na bieżąco modyfikowany. | Domyślny dla nowych folderów. |
| **core-ready** | L3, gotowy do awansu, ale jeszcze nie w core. | `level: L3` + `source_of_status` z PASS-ami + brak otwartych blokujących bridge'y. |
| **core-promoted** | Wynik **JEST** w `main.tex` / `core/` / `axioms/`. Folder żyje jako historia. | `promoted_to_core` z weryfikowalną referencją. |
| **needs-bridge** | Numerka istnieje, brak analitycznego mostu LUB folder oczekuje na rezultat z innego. | `open_bridges` niepuste. |
| **needs-migration** | Folder wymaga aktualizacji po zmianie core. | `core_compatibility ∈ {stale, broken}` + powiązany hotspot w `CORE_HOTSPOTS.md`. |
| **audit** | Audyt wewnętrzny / cross-check. | Folder nie produkuje fizyki, tylko analizuje. |
| **review** | Audyt zewnętrzny / external review. | Np. `external_review_2026-04-25/`. |
| **program-doc** | Dokument programowy (status / redirect). | Po decyzji Q3 z 2026-05-03 outdated trafia do `_archiwum/`. |
| **sandbox** | Luźna eksploracja, NIE część teorii. | Fizycznie w `research/_sandbox/`. |
| **archive** | Stara ścieżka badawcza, obsoletna. | Fizycznie w `research/_archive/`. WITHDRAWN/OBSOLETE/SUPERSEDED z 2 źródłami + akceptacja człowieka. |

---

## C. Poziomy dojrzałości (`tgp_status.level`)

| Termin | Definicja | Required evidence |
|--------|-----------|-------------------|
| **L0** | Pomysł / sandbox; brak skryptów, brak weryfikacji. | Tekst hipotezy. |
| **L1** | Szkic numeryczny; skrypt biegnie, są pliki `.txt` z wynikami. Brak interpretacji analitycznej. | ≥1 skrypt z output. |
| **L2** | Analityczny szkic; derivation lub bridge istnieje, ale weryfikacja częściowa lub ≥ 1 luka znana. | ≥1 PASS lub strukturalne testy częściowe. |
| **L3** | Domknięty wynik; `Phase{1..N}_results.md` zamknięte z **≥ 90% PASS**. Brak otwartych *blokujących* bridge'y. | `core_compatibility ∈ {current, partial}` z explicite wyjaśnieniem. |
| **L4** | Promowany do core; wynik **JEST** w `main.tex` / `core/`. | `promoted_to_core` ścieżka + sekcja. |
| **mixed** | Folder zawiera podzbiory na różnych poziomach. | Typowo closure-aggregator. |
| **unknown** | Przed klasyfikacją lub po nieudanym audycie. | Zawsze legal. |

**Anty-overclaim:** L4 wymaga `promoted_to_core` weryfikowalnego;
L3 wymaga `source_of_status` cytującego konkretny plik+sekcję.

---

## D. Charakter pracy (`tgp_status.kind`)

| Termin | Definicja |
|--------|-----------|
| **derivation** | Analityczne wyprowadzenie. |
| **numerical** | Głównie skrypty, wyniki numeryczne. |
| **phenomenology** | Dopasowania do danych, predykcje obserwacyjne. |
| **bridge** | Folder dostarcza brakujący lemat / most analityczny. |
| **audit** | Audyt / cross-check. |
| **review** | Review zewnętrzny. |
| **program-doc** | Dokument programowy. |
| **closure-aggregator** | Folder-rodzic z wieloma podzamknięciami. |

---

## E. Zgodność z core (`tgp_status.core_compatibility`)

| Termin | Definicja |
|--------|-----------|
| **current** | W pełni zgodny z aktualnym `main.tex` / `core/`. |
| **partial** | Zgodny w części; istnieją **explicite** udokumentowane lokalne tensje. |
| **stale** | Folder poprawny pod swoje założenia, ale core się ruszył (typowo M9.x post-pivot M9.1''). Wymaga re-runu. |
| **broken** | Strukturalna sprzeczność z core — należy traktować jako `needs-migration`. |
| **unknown** | Niesprawdzane lub niemożliwe do określenia bez human review. |

---

## F. Statusy hotspotów (`CORE_HOTSPOTS.md`)

| Termin | Kto może ustawić |
|--------|------------------|
| **OPEN** | Domyślny — z audytu lub stwierdzenia agenta. |
| **READY-FOR-INTAKE** | Agent — gdy znajdzie folder, który adresuje hotspot. |
| **IN-INTAKE** | Agent — gdy stworzy `INTAKE_*.md`. |
| **PROMOTING** | **Tylko człowiek** — gdy zaakceptuje INTAKE i rozpocznie edycję core. |
| **RESOLVED** | **Tylko człowiek** — po fizycznej zmianie w core. |

---

## G. Statusy INTAKE (`meta/core/intake/INTAKE_*.md`)

| Termin | Kto może ustawić |
|--------|------------------|
| **PENDING-HUMAN-REVIEW** | Agent — stan początkowy. |
| **ACCEPTED-AND-PROMOTED** | Człowiek. |
| **DEFERRED** | Człowiek. |
| **REJECTED** | Człowiek. |
| **NEEDS-MORE-EVIDENCE** | Człowiek. |

---

## H. Statusy bridge (`CANDIDATE_BRIDGES.md`)

| Termin | Definicja |
|--------|-----------|
| **PROPOSED** | Wpis utworzony przez agenta (Sesja 6). |
| **HUMAN-CONFIRMED** | Człowiek zatwierdził match. |
| **EXECUTED** | Bridge zrealizowany — wynik skopiowany do source FINDINGS z notatką pochodzenia. |
| **REJECTED** | Nie pasuje, z uzasadnieniem. |

| Match strength | Definicja |
|----------------|-----------|
| **EXACT** | Formuła / wartość liczbowa identyczna z dokładnością do tolerancji. |
| **PARTIAL** | Wynik w odpowiedniej formie, wymaga adaptacji (jednostek, frame'a). |
| **HEURISTIC** | Tylko podobieństwo strukturalne. Wymaga ręcznego przeglądu. |

---

## I. Statusy bus (`RESEARCH_BUS.md`)

| Termin | Definicja |
|--------|-----------|
| **BROADCAST** | Wpis świeży, nieprzeczytany przez konsumentów. |
| **CONSUMED** | ≥1 konsument potwierdził użycie (dopisał do `depends_on`). |
| **STALE** | ≥ 90 dni bez akcji. Trafia do § 5 historii. |

---

## J. Anty-overclaim — terminy zakazane

Z `AGENT_PROTOCOL.md` §3 + audytu `74394a8`:

| Termin | Kontekst zakazu |
|--------|-----------------|
| **FULL CONVERGENCE** | Zarezerwowane. NIGDY w `FINDINGS.md`, `source_of_status`, README. |
| **DERIVED** (bez referencji do core) | Tylko jeśli `level: L4` + weryfikowalne `promoted_to_core`. W innych przypadkach: `STRUCTURAL`, `PARTIALLY DERIVED`, `STRUCTURAL HINT`. |
| **LOCKED** (bez liczby PASS i pliku) | Wymaga konkretnego `Phase{N}_results.md` z PASS X/Y. |
| **CLOSED** (bez pliku) | Wymaga `Phase*_results.md` lub `KNOWN_ISSUES.md` cytatu. |

---

## K. Role plików meta/

| Plik | Rola | Kto edytuje |
|------|------|-------------|
| `meta/PLAN_RESEARCH_WORKFLOW_v1.md` | Master plan workflow | Człowiek + agent (sekcje statusowe) |
| `meta/research/AGENT_PROTOCOL.md` | Reguły dla agentów | Człowiek (gate akceptacji) |
| `meta/research/AUDIT_RESEARCH_S1.md` | Snapshot audytu Sesji 1 | Generowany skryptem (`_audit_s1*.py`) |
| `meta/research/RESEARCH_BUS.md` | Tablica ogłoszeń międzyfolderowych | Agent (każdy nowy `FINDINGS.md` item) |
| `meta/research/CANDIDATE_BRIDGES.md` | Match `NEEDS` × `FINDINGS` | Agent (Sesja 6) + człowiek (HUMAN-CONFIRMED) |
| `meta/research/CORE_CANDIDATES.md` | Lista core-ready / core-promoted | Sesja 7 + INTAKE workflow |
| `meta/research/IMPACT_MATRIX.md` | Macierz wpływu zmian | Agent (auto-import) + Sesja 6/7 |
| `meta/research/FOLDER_STATUS_INDEX.md` | Globalna mapa statusu | Skrypt `_build_status_index.py` (Sesja 4) |
| `meta/research/MIGRATION_LOG.md` | Log fizycznych ruchów | Agent (każde mv/delete) + akceptacja człowieka |
| `meta/research/GLOSSARY.md` | Ten plik | Człowiek (statyczny) |
| `meta/research/templates/*` | Szablony README/NEEDS/FINDINGS/STATUS_BLOCK | Człowiek (statyczne) |
| `meta/research/_examples/*` | Podglądy szablonów | Człowiek (Sesja 2) |
| `meta/core/CORE_INVENTORY.md` | Mapa rdzenia (co w core, skąd) | Człowiek (przez INTAKE) |
| `meta/core/CORE_HOTSPOTS.md` | Living checklist core sprzeczności | Agent (OPEN→READY/IN-INTAKE) + człowiek (PROMOTING/RESOLVED) |
| `meta/core/CORE_INTAKE.md` | Protokół promocji | Człowiek (statyczny) |
| `meta/core/intake/INTAKE_*.md` | Wnioski o promocję | Agent (utworzenie) + człowiek (decyzja) |

---

## L. Decyzje człowieka — referencje

| Decyzja | Plik | Sekcja |
|---------|------|--------|
| Q-2026-05-02 #1–6 | `meta/PLAN_RESEARCH_WORKFLOW_v1.md` | §9 |
| Q-2026-05-03 Q1–Q4 | `meta/PLAN_RESEARCH_WORKFLOW_v1.md` | §9.1 |
| Stan repo po decyzjach | `meta/PLAN_RESEARCH_WORKFLOW_v1.md` | §9.2 |
| Migracje fizyczne | `meta/research/MIGRATION_LOG.md` | § 2 |
