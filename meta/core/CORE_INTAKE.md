---
title: "CORE_INTAKE — protokół promocji research/ → core/"
date: 2026-05-03
type: protocol
status: STATIC v1.0 (Sesja 3)
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/core/CORE_INVENTORY.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
  - "[[meta/research/CORE_CANDIDATES.md]]"
tags:
  - protocol
  - core-promotion
  - human-gate
---

# CORE_INTAKE — protokół promocji `research/` → `core/`

> **Twardy gate.** Każda promocja folderu badawczego do statusu
> `core-promoted` (`level: L4`) wymaga:
>
> 1. Wniosku `meta/core/intake/INTAKE_<folder>_<date>.md` od agenta.
> 2. Akceptacji człowieka.
> 3. Manualnego edytu core przez **człowieka** (NIGDY agenta).
> 4. Aktualizacji `CORE_INVENTORY.md`.
> 5. Aktualizacji YAML źródłowego folderu (`promoted_to_core: <ref>`).
>
> **Agent jest tylko proponentem.** Nie wykonuje promocji.

---

## 1. Kiedy agent powinien stworzyć INTAKE

Agent tworzy `INTAKE_<folder>_<date>.md`, gdy:

1. Folder ma `tgp_status.level >= L3` z weryfikowalnym `source_of_status`.
2. Folder ma `folder_status: core-ready` (po klasyfikacji w Sesji 7).
3. Wynik folderu adresuje hotspot z [[CORE_HOTSPOTS.md]] (idealnie
   `READY-FOR-INTAKE`).
4. Folder NIE jest na liście kwarantanny (`polluted_74394a8: false`).

Agent NIE tworzy INTAKE, gdy:

- Folder ma `level: L1/L2` lub `unknown`.
- `source_of_status` jest pusta lub zawiera tylko keywordy.
- Folder ma `polluted_74394a8: true`.
- Folder ma `core_compatibility: stale | broken` (najpierw migracja).

---

## 2. Format wniosku INTAKE

Plik: `meta/core/intake/INTAKE_<folder-slug>_<YYYY-MM-DD>.md`

**Slug**: nazwa folderu badawczego z myślnikami (np.
`INTAKE_op-uv-as-ngfp_2026-05-15.md`).

```markdown
---
title: "INTAKE — <folder> → core (<topic>)"
date: <YYYY-MM-DD>
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/<folder>
target_section: <main.tex §X.Y | core/sek0Z_* | axioms/* | partial_proofs/*>
hotspot_addressed: <H-A1 | H-B6 | none>
agent_requesting: <session ID, np. "Sesja 7 batch N">
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/<folder>/README.md]]"
  - "[[research/<folder>/FINDINGS.md]]"
  - (link do hotspota, jeśli adresuje)
tags:
  - intake
  - core-promotion
---

# INTAKE — research/<folder> → core (<topic>)

## 1. Wynik proponowany do promocji

<Jednolinijkowe stwierdzenie wyniku, np.: "M9.1'' f(ψ) = (4-3ψ)/ψ z
n = deg(V) = 4 unique exponent (T-FP principle)".>

## 2. Source w research/

<Pliki w research/<folder>, które zawierają wynik. Każdy z dokładnym
cytowaniem (plik:linia / sekcja / PASS-count).>

- `research/<folder>/Phase3_results.md §X.Y` — <co tu jest>
- `research/<folder>/FINDINGS.md F<id>` — <eksportowany wynik>
- `research/<folder>/<plik>.txt` — <numeryczny PASS X/Y>

## 3. Proponowana lokacja w core

| Target | Co dokładnie ma się tam pojawić | Czy istnieje teraz | Diff (dla istniejących sekcji) |
|--------|----------------------------------|---------------------|---------------------------------|
| `<core/sek0Z_*.tex>` | <statement / formuła / dowód> | NO | (nowa sekcja) |
| `<main.tex §X.Y>` | <reference / tabela predykcji> | YES (lin. <N>) | <konkretny diff> |

## 4. Hotspot, który ten wynik zamyka (jeśli)

<Link do `meta/core/CORE_HOTSPOTS.md §H-XN`. Cytat hotspota.>

## 5. Co MUSI być sprawdzone przed promocją

- [ ] Wynik kompiluje się (sympy/numerical) na bieżącej wersji core
- [ ] Brak konfliktu z aktualnym `main.tex` (read-only sprawdzenie)
- [ ] `core_compatibility: current` w YAML folderu źródłowego
- [ ] Wszystkie testy w `Phase*_results.md` PASS (X/Y zawarte w wniosku)
- [ ] `source_of_status` w YAMLu zawiera linki do plików source

## 6. Co się STANIE po akceptacji człowieka

1. **Człowiek** edytuje `core/<target>` ręcznie.
2. **Człowiek** dopisuje wpis do [[meta/core/CORE_INVENTORY.md]].
3. **Człowiek** zmienia w `research/<folder>/README.md`:
   - `folder_status: core-promoted`
   - `level: L4`
   - `promoted_to_core: "<konkretna ścieżka>"`
   - `last_yaml_update: <data promocji>`
4. **Człowiek** zmienia status hotspota w `CORE_HOTSPOTS.md` na `RESOLVED`
   (jeśli był).
5. **Człowiek** zamyka ten wniosek INTAKE: ustawia `status: ACCEPTED-AND-PROMOTED`
   + dopisuje datę.

## 7. Decyzja człowieka

- [ ] **ACCEPT** — wynik promowany do core (data: ____)
- [ ] **DEFER** — wynik OK, ale promocja przesunięta (powód: ____)
- [ ] **REJECT** — wynik nie nadaje się do core (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — wymaga dopełnienia w `research/<folder>` (lista: ____)

---

(po decyzji człowieka — wpisuje się notatkę o akceptacji + ewentualne
poprawki, zostaje archiwum decyzji)
```

---

## 3. Stany cyklu życia INTAKE

```
PENDING-HUMAN-REVIEW   ← stan początkowy (agent stworzył wniosek)
        ↓
        ├──→ ACCEPTED-AND-PROMOTED  (człowiek wykonał promocję)
        ├──→ DEFERRED               (człowiek odracza)
        ├──→ REJECTED               (człowiek odrzuca)
        └──→ NEEDS-MORE-EVIDENCE    (człowiek żąda dopełnienia)
```

- **Agent** może ustawić tylko `PENDING-HUMAN-REVIEW`.
- **Człowiek** może ustawić każdy z 4 końcowych stanów.
- `NEEDS-MORE-EVIDENCE` → agent może (po dopełnieniu w research) stworzyć
  nowy wniosek `INTAKE_<folder>_<new-date>.md` z linkiem zwrotnym
  do poprzedniego.

---

## 4. Anty-overclaim w INTAKE

Wniosek INTAKE jest **traktowany jak research-level finding** w sensie
reguł `AGENT_PROTOCOL.md` §3:

- "Wynik proponowany do promocji" musi cytować plik + sekcję + liczbę
  PASS-ów. **Bez cytatu wniosek jest nieważny.**
- "Source w research/" musi mieć ≥ 2 niezależne pliki source (chyba że
  folder jest aggregator-em z jednym SUMMARY).
- "Diff (dla istniejących sekcji)" musi być **konkretny** — line numbers
  + przed/po (cytat). Nie wystarczy "update metryka".

---

## 5. Czego INTAKE NIE robi

- **Nie modyfikuje core.** Modyfikacja jest aktem człowieka.
- **Nie generuje predykcji.** Predykcje pochodzą z `research/` i są w
  `FINDINGS.md` źródłowego folderu.
- **Nie aktualizuje `INDEX.md`, `DEPENDENCIES.md`, `PREDICTIONS_REGISTRY.md`**
  ani innych globalnych ledgerów. Te aktualizacje są częścią aktu
  promocji i wykonuje je człowiek.
- **Nie obejmuje migracji** — jeśli folder ma `core_compatibility: stale`
  lub `broken`, najpierw idzie ścieżką `needs-migration` w
  `IMPACT_MATRIX.md`, a INTAKE jest wstrzymany.

---

## 6. Obecny stan intake/

Folder: [[meta/core/intake/]]
Plik: `meta/core/intake/.gitkeep` (placeholder).

Wnioski INTAKE pojawią się w Sesji 7 (klasyfikacja core-ready /
core-promoted) lub później (gdy agent w pracy nad folderem stwierdzi
gotowość).

**Pierwsze realne kandydaty do INTAKE** (na bazie audytu Sesji 1 +
[[CORE_HOTSPOTS.md]] §B):

- `research/closure_2026-04-26/sigma_ab_pathB` → `core/sek08_formalizm`
  (adresuje H-B6: tensor σ_ab Path B PRIMARY)
- `research/closure_2026-04-26/f_psi_principle` → `core/sek08c_metryka_z_substratu`
  (adresuje H-B8: f(ψ) deeper principle T-FP)
- `research/closure_2026-04-26/Lambda_from_Phi0` → `core/sek05_ciemna_energia`
  (adresuje vacuum catastrophe; Ω_Λ input → prediction)
- `research/closure_2026-04-26/alpha_psi_threshold` → `core/sek06` lub
  nowa sekcja (adresuje H-B9: WEP margin 4×10¹⁶×)

Wszystkie 4 są w stanie `READY-FOR-INTAKE` w `CORE_HOTSPOTS.md`. Sesja 7
sklasyfikuje je formalnie i (jeśli `core_compatibility: current`) stworzy
wnioski INTAKE.
