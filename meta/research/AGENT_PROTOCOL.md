---
title: "AGENT_PROTOCOL — reguły dla agentów pracujących w research/"
date: 2026-05-03
type: protocol
status: STUB v0.2 (Sesja 2) — pełna wersja w Sesji 3
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/templates/STATUS_BLOCK.yaml]]"
  - "[[meta/research/templates/README.template.md]]"
  - "[[meta/research/templates/NEEDS.template.md]]"
  - "[[meta/research/templates/FINDINGS.template.md]]"
  - "[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
tags:
  - protocol
  - multi-agent
  - anti-overclaim
  - workflow
---

# AGENT_PROTOCOL — reguły dla agentów pracujących w `research/`

> **Status:** STUB Sesji 2. Sesja 3 rozszerzy o pełny opis multi-agent
> orchestracji (locki, sesje, RESEARCH_BUS workflow). Tu są spisane
> reguły, które są już TWARDE i obowiązują od momentu rozpoczęcia
> Sesji 4.

---

## 0. Najważniejsze ostrzeżenia

1. **NIGDY** nie edytuj `core/`, `axioms/`, `partial_proofs/`,
   `papers_external/`, `main.tex`, `tgp_letter.tex`, `tgp_companion.tex`,
   `arxiv_submission/`, `paper_lepton_masses/`, `paper_bh_shadow/`
   bez **explicite rozkazu człowieka** w bieżącej sesji.
2. **NIGDY** nie nadawaj `level: L4` ani `folder_status: core-promoted`
   bez weryfikowalnego `promoted_to_core` cytującego konkretny plik
   core (zob. § 3 niżej).
3. **NIGDY** nie wpisuj „FULL CONVERGENCE", „DERIVED", „LOCKED",
   „CLOSED" do `FINDINGS.md` lub `source_of_status` bez cytowania pliku
   + sekcji + liczby PASS-ów.
4. **NIGDY** nie usuwaj plików pod `research/` bez wpisu w
   [[MIGRATION_LOG.md]] i akceptacji człowieka.
5. **NIGDY** nie pomijaj zatrutych folderów z incydentu `74394a8`
   (lista poniżej, § 6) — automatycznie stosuj klauzulę kwarantanny.

---

## 1. Definicje poziomów L0–L4

| Level | Co znaczy | Wymagane evidence |
|------:|-----------|--------------------|
| **L0** | Pomysł / sandbox | Tekst hipotezy. Brak skryptu. |
| **L1** | Szkic numeryczny | Skrypt biegnie, są pliki `.txt` z wynikami. Brak interpretacji analitycznej. |
| **L2** | Analityczny szkic | Derivation lub bridge istnieje. Weryfikacja częściowa (≥ 1 PASS) lub testy strukturalne wykonane, ale ≥ 1 luka znana. |
| **L3** | Domknięty wynik | Phase{1..N}_results.md zamknięte z **≥ 90% PASS**. Brak otwartych bridge'y *blokujących*. `core_compatibility ∈ {current, partial}` z explicite wyjaśnioną częściową niezgodnością. |
| **L4** | Promowany do core | Wynik **JEST** w `main.tex` / `core/` / `axioms/` / `partial_proofs/`. Pole `promoted_to_core` zawiera ścieżkę + sekcję. |

**Mixed**: folder zawiera podzbiory na różnych poziomach
(typowo closure-aggregator, np. `closure_2026-04-26/`).
**Unknown**: przed klasyfikacją lub po nieudanym audycie.

### 1.1 Zakaz „aspiracyjnego" L4

Jeśli wynik *powinien* być w core, ale jeszcze go tam nie ma:
to jest **L3 + folder_status: core-ready** plus **wniosek**
`meta/core/intake/INTAKE_<folder>_<date>.md` (decyzja Q6 z 2026-05-02:
agent ma proponować awanse, ale nie wykonywać). Ustawienie L4 bez
referencji w core jest **falszywym wpisem ledger** w sensie audytu
2026-05-02 (zob. [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]).

---

## 2. Reguła „source_of_status"

Każda wartość statusu, która jest „mocniejsza" niż domyślna
(`active` + `L0/L1/unknown`), wymaga wpisu w `tgp_status.source_of_status`.
Wpis ma jeden z następujących formatów:

| Format | Przykład |
|--------|----------|
| Plik wewnątrz folderu + sekcja / nazwa testu | `Phase3_results.md (94/97 PASS, 2026-04-25)` |
| Audyt globalny + paragraf | `meta/AUDYT_TGP_2026-05-01.md §A1 (CRITICAL)` |
| Sąsiedni folder + plik | `closure_2026-04-26/sigma_ab_pathB/results.md (11/11 PASS)` |
| Wpis w core / main.tex | `main.tex §8.3 (eq:m911-metric)` (WYŁĄCZNIE dla L4) |

**Niedozwolone formaty:**
- `general intuition`
- `obvious from setup`
- `the script clearly shows...`
- `as discussed earlier`
- pojedyncze keywordy bez pliku (np. tylko `PASS` lub `CLOSED`)

Reguła twarda: brak `source_of_status` ⇒ status nie może być wyższy niż
`{folder_status: active, level: L1}` (lub `sandbox/L0`).

---

## 3. Reguła „anty-overclaim"

Reguła ma **dwa case studies symetryczne** — błędy w przeciwnych kierunkach,
oba popełnione w obrębie tego workflow:

### Case study #1 — over-claim (`74394a8`, subagent)

Subagent wprowadził 4 fałszywe „FULL CONVERGENCE" cykli, modyfikując
globalny ledger (`INDEX.md` 784 → 856 master count) bez weryfikowalnej
podstawy (zob. [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §1).

**Klucz:** dopisanie PASS-ów / DERIVED / LOCKED bez source z konkretnego
pliku.

### Case study #2 — under-claim (Sesja 3, ten agent)

Ten agent w Sesji 3 transkrybował 20 hotspotów z `AUDYT_TGP_2026-05-01.md`
(sekcje § A "CRITICAL" + § B "HIGH") do `meta/core/CORE_HOTSPOTS.md`,
**wszystkie oznaczając jako OPEN**. Faktycznie ten sam dokument audytu
zawiera sekcje § J–§ AB ("Naprawa…") z **43/43 = 100% closure** (zob.
[[meta/research/HOTSPOT_AUDIT_S3_5.md]]).

**Klucz:** czytanie tylko sekcji diagnozy i ignorowanie sekcji naprawy
w tym samym pliku.

**Symetria:** case #1 = "twierdzę zamknięte gdy jest otwarte"; case #2 =
"twierdzę otwarte gdy jest zamknięte". Oba są naruszeniem reguły
`source_of_status` — pierwszy przez fabrykację, drugi przez niedopełnienie
przeczytania źródła do końca.

### 3.0 Twarda reguła wynikowa

**Każdy "audit-driven" wpis (hotspot / intake / candidate) wymaga
przeczytania całego dokumentu źródłowego, włącznie z final closure
summary, zanim jakikolwiek item zostanie wpisany jako OPEN.**

Audyty TGP zawierają zarówno diagnozę jak i naprawę w jednym pliku.
Skanowanie nazw paragrafów (`grep '^## '`) jest **pierwszym krokiem**,
nie ostatnim — zawsze trzeba wyłuskać `final tally` / `closure aggregate`
/ `total closed/open` (typowo w sekcji "Final" / "Status final" / "§ S").

### 3.0.1 Cite-only ekstrakcja MUSI filtrować banned phrases

Reguła wynikowa z Sesji 8 phrase-cleanup (2026-05-03):

> Cite-only extractor (np. `_extract_s5.py`) **NIE WOLNO** propagować
> banned phrases (lista w § 3.1) do auto-generated artefaktów
> (`FINDINGS.md`, `NEEDS.md`, broadcast w `RESEARCH_BUS.md`),
> nawet jeśli źródło ich używa.

**Implementacja:** `sanitize_banned_phrases()` w `_extract_s5.py`
zastępuje banned phrases placeholderem `[phrase redacted per
AGENT_PROTOCOL §3]` przed zapisem do filesystem, zachowując PASS-count
gdzie to możliwe (np. "6/6 PASS — FULL CONVERGENCE" → "6/6 PASS [phrase
redacted per AGENT_PROTOCOL §3]").

**Powód:** Sesja 8 regression check (2026-05-03) wykrył 3 MEDIUM
findings z "FULL CONVERGENCE" propagated z Phase3_results.md frontmatter
do auto-generated FINDINGS.md przez cite-only extractor. Heurystyka
"cytujemy bez modyfikacji" generuje konflikt z anti-overclaim regułą.
Cite-with-sanitize jest właściwym kompromisem.

### 3.1 Twarde zasady (oba case studies adresowane)

> ⚠ **2026-05-04: Phase 0 balance sheet wymagany dla nowych claimów DERIVED.**
> Pełny protokół: [[../CALIBRATION_PROTOCOL.md]] — Phase 0 balance sheet
> + tautology test + falsifiability test + independent-path cross-validation.
> **Brak balance sheet → max status STRUCTURAL.**

### 3.1 Twarde zasady

1. **„FULL CONVERGENCE" jest zarezerwowane** — agent nie używa tej frazy
   nigdzie w `FINDINGS.md`, `source_of_status` ani w README. Jeśli folder
   zasługuje na taki werdykt, wynika to z kombinacji `level: L3` +
   `Phase*_results.md` z konkretną liczbą PASS-ów + Phase 0 balance sheet
   per [[../CALIBRATION_PROTOCOL.md]] §2.
2. **„DERIVED" wymaga referencji do core + Phase 0 balance sheet** — jeśli
   wynik nie jest w core lub balance sheet brak, używaj `STRUCTURAL`,
   `PARTIALLY DERIVED`, `STRUCTURAL HINT`, `ANSATZ`, `NUMEROLOGICAL OBSERVATION`
   (zgodnie z taksonomią [[../CALIBRATION_PROTOCOL.md]] §1).
3. **Nie modyfikujemy globalnych ledgerów (`INDEX.md`,
   `PREDICTIONS_REGISTRY.md`, `DEPENDENCIES.md`)** w sesjach 1–8 bez
   explicite rozkazu człowieka. Te pliki są częścią rdzenia.
4. **Cumulative test counts zostają invariantem** — agent nie zwiększa
   liczników typu „281 cumulative verifications" bez fizycznych testów
   z plikiem wynikowym + akceptacji człowieka.
5. **Każdy wpis statusu jest dodawany WRAZ z `last_yaml_update`** w polu
   YAML, żeby audyt dat działał. Skrypt-checker w Sesji 8 sprawdza:
   `last_yaml_update <= today` i `>= mtime_najnowszy_w_folderze`.

### 3.2 Lista zatrutych folderów (klauzula kwarantanny)

Z [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]:

- `research/op-chi1-newton-constant-derivation/`
- `research/op-uv2-mtgp-absolute-scale/`
- `research/op-omega2-axion-coupling-lock/`
- `research/op-omega3-axion-decay-constant/`

Dla TYCH folderów do czasu forward-patch (decyzja człowieka):
- `level: unknown`
- `core_compatibility: unknown`
- `folder_status: needs-bridge`
- `polluted_74394a8: true`
- `source_of_status` zawiera notatkę:
  `"Polluted 74394a8 — frozen pending forward-patch (zob. meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md)"`

Agent **NIE WOLNO** awansować ich do `core-ready` lub `core-promoted`,
**NIE WOLNO** dodawać `INTAKE_*.md` dla nich, **NIE WOLNO** zmieniać
ich statusu na `archive` (to zatuszowałoby ślad).

---

## 4. Reguły wpisu do `RESEARCH_BUS.md`

Każdy nowy wynik (każdy nowy item w `FINDINGS.md`) wymaga wpisu w
[[RESEARCH_BUS.md]] z polem `consumers:` (lista folderów, które
mogą skonsumować wynik). Jeśli `consumers: []`, wpis jest dopuszczalny,
ale agent musi explicite napisać `# scanned, no consumers found
2026-MM-DD`.

Format wpisu (do uszczegółowienia w Sesji 3):

```
| Date | Source folder | Finding ID | Headline | Consumers |
| 2026-05-03 | op-uv-as-ngfp | F1 | "N_A = 8.7719 ± 0.068% (PARTIALLY DERIVED)" | op-bh-alpha-threshold, op-eps-photon-ring |
```

---

## 5. Reguły wpisu do `CANDIDATE_BRIDGES.md`

Match między `NEEDS.md` jednego folderu a `FINDINGS.md` innego musi
spełniać:

1. Cytat **dokładnej formuły / wartości / lemma** z obydwu stron.
2. Sam keyword nie wystarczy — np. dwa różne „α" w różnych kontekstach
   (α_em vs α_QCD vs α dim) NIE tworzą bridge'a.
3. Match musi być przeglądalny przez człowieka — wpis ma format
   „Folder A potrzebuje X (cytat); folder B ma X (cytat)".

Patrz Sesja 6 planu.

---

## 6. Konflikt ze stanem core

Jeśli agent w trakcie pracy w folderze odkryje sprzeczność z aktualnym
`core/` / `main.tex` / `axioms/`:

1. **STOP.** Nie modyfikuj core.
2. Ustaw `tgp_status.folder_status: needs-migration` w README folderu
   źródłowego.
3. Ustaw `core_compatibility: stale | broken` (zależnie od skali).
4. Dopisz wpis do [[CORE_HOTSPOTS.md]] (nowy hotspot z `OPEN`) —
   chyba że hotspot już tam jest, wtedy dopisz folder do
   `Linked needs-migration folders`.
5. Wpisz do [[MIGRATION_LOG.md]] z datą + opisem.
6. **Zatrzymaj sesję** — czekaj na human review.

Niedozwolone: „naprawienie" sprzeczności w `research/` przez retroaktywne
zmienianie wyników, żeby pasowały do core. Jeśli core zmienił reguły
i folder jest stale, to **folder** jest stale, nie odwrotnie.

---

## 7. Konwencje ścieżek

- Wszystkie ścieżki wewnątrz `meta/research/*` używają **względnych**
  ścieżek od `TGP/TGP_v1/` (np. `research/op-uv-as-ngfp/Phase3_results.md`).
- Wikilinks: domyślnie używaj formy z pełną ścieżką
  (`[[research/op-uv-as-ngfp/Phase3_results.md]]`) — Obsidian rozwiąże,
  a jednocześnie audyt linków (Sesja 8) ma jasny target.
- Wikilinks w `research/<folder>/...` mogą używać formy bez ścieżki
  (`[[Phase3_results.md]]` w obrębie tego samego folderu) lub
  `[[../<inny>/<plik>]]` dla bliskich sąsiadów.

---

## 8. Decyzje człowieka odzwierciedlone w protokole

| Decyzja | Implikacja w protokole |
|---|---|
| Q-2026-05-02 #1 (canonical) | Ścieżki względne od `TGP/TGP_v1/`. Vault root `./research/` nie istnieje. |
| Q-2026-05-02 #2 (wariant minimalny) | Brak fizycznych nadfolderów statusowych. Status w YAMLu. |
| Q-2026-05-02 #3 (research/_archive/) | Sesja 8.5 jest jedyną sesją z `git mv` (bez decyzji ad-hoc). |
| Q-2026-05-02 #6 (proponuj L4) | Agent **może** tworzyć `meta/core/intake/INTAKE_*.md` z weryfikowalnym cytatem. **Nie może** modyfikować core. |
| Q-2026-05-03 #2 (tylko rodzic) | Subfoldery z podzamknięciami nie mają osobnego YAML/NEEDS/FINDINGS. Zsumowane w `subfolder_summary` rodzica. |
| Q-2026-05-03 #4 (nie nadpisywać) | Jeśli `FINDINGS.md`/`NEEDS.md` istnieje przed Sesją 5 — flag `pre_existing_findings: true` / `pre_existing_needs: true`, agent **NIE nadpisuje**, dopisuje tylko strukturalny header. |

---

## 9. Self-check przed zakończeniem zadania (per folder)

Każdy agent kończąc pracę z folderem **musi** odpowiedzieć sobie
na te pytania (i opcjonalnie zalogować w `MIGRATION_LOG.md` przy
pracach większych niż edycja YAML):

1. Czy `tgp_status.last_yaml_update` jest dzisiejszy?
2. Czy `source_of_status` cytuje konkretny plik + sekcję dla każdego
   pola statusu „mocniejszego" niż domyślny?
3. Czy `level` jest spójny z `folder_status` (zob. § 1)?
4. Czy `promoted_to_core` jest null **lub** zawiera weryfikowalną
   ścieżkę do core?
5. Czy nowe wyniki w `FINDINGS.md` są zsynchronizowane z
   `RESEARCH_BUS.md`?
6. Czy nowe luki w `NEEDS.md` są zsynchronizowane z
   `tgp_status.open_bridges`?
7. Jeśli zmienił się `core_compatibility` — czy `CORE_HOTSPOTS.md` /
   `IMPACT_MATRIX.md` zostały zaktualizowane?

---

## 10. Co protokół świadomie pomija (do Sesji 3)

- Mechanika multi-agent locków (kto edytuje co, jak unikać konfliktów).
- Pełna mechanika `RESEARCH_BUS.md` (subscribers, ack, diff broadcast).
- Algorytm batchowy dla Sesji 4 (10 folderów / batch — kolejność,
  resume, partial state).
- Procedura `INTAKE_*.md` (forma wniosku, co dokładnie human gate
  zatwierdza).

To zostanie dopisane w Sesji 3, gdy `meta/research/` zostanie ustanowione
jako żywa warstwa.
