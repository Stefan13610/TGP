---
title: "IMPACT_MATRIX — macierz wpływu zmian między folderami research/ i core"
date: 2026-05-03
type: matrix
status: SKELETON (Sesja 3) — wypełniany w Sesji 6+7
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[DEPENDENCIES.md]]"
  - "[[DEPENDENCIES_REVERSE.md]]"
  - "[[meta/research/RESEARCH_BUS.md]]"
  - "[[meta/research/CANDIDATE_BRIDGES.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
tags:
  - impact
  - matrix
  - dependencies
  - regression
---

# IMPACT_MATRIX — macierz wpływu zmian

> **Cel pliku:** odpowiedzieć na pytanie „jeśli zmienię X, co jeszcze
> wymaga uwagi?". Wzbogacenie auto-generowanego `DEPENDENCIES.md`
> o **ludzkojęzyczne** anotacje: rodzaj wpływu, severity, hotspot
> linkowy. Punkt wejścia dla regression review w Sesji 8.
>
> **Źródła wpisów:**
> 1. `tgp_status.depends_on` / `impacts` w YAMLach folderów (Sesja 4).
> 2. Auto-import z `DEPENDENCIES.md` / `DEPENDENCIES_REVERSE.md`
>    (wykonywany skryptowo).
> 3. Ręczne anotacje z audytów (`AUDYT_TGP_*.md`,
>    `meta/core/CORE_HOTSPOTS.md`).
> 4. Match z `CANDIDATE_BRIDGES.md` (Sesja 6).

---

## 0. Reguły wpisu

1. **Wiersz** = pojedyncza krawędź wpływu `Source → Target`.
2. **Kierunek**: Source = folder, którego zmiana wpływa na Target.
3. **Severity**:
   - **CRITICAL** — zmiana w Source unieważnia wynik w Target
     (np. M9.1'' pivot unieważnia M9.2 m_field — H-A2)
   - **HIGH** — zmiana w Source wymaga re-runu w Target, ale wynik może
     przetrwać
   - **MEDIUM** — Source dostarcza wejścia dla Target, ale Target ma
     własne testy spójności
   - **LOW** — luźne powiązanie tematyczne
4. **Edge type**:
   - `data` — Source produkuje liczbę / formułę używaną w Target
   - `lemma` — Source dostarcza lemat / dowód używany w Target
   - `axiom` — Source jest aksjomatem / postulatem dla Target
   - `methodology` — Source dostarcza metody / formalizmu
   - `hotspot` — krawędź pochodzi z hotspota core
5. **Linked hotspot** (opcjonalnie): ID z `CORE_HOTSPOTS.md` (np. `H-A2`).

---

## 1. Statystyka

| Severity | Liczba krawędzi |
|----------|----------------:|
| CRITICAL | 0 |
| HIGH | 0 |
| MEDIUM | 0 |
| LOW | 0 |
| **Razem** | **0** |

(stub — Sesja 3 initial seed)

**Z `DEPENDENCIES.md`** (auto, generated 2026-04-22):
- 1657 krawędzi total: `\input` 69, `\ref` 1436, `[[wiki]]` 152.
- 91 folderów (fine), 10 top-level (coarse).
- 0 orphans.

Tę liczbę traktujemy jako **dolną granicę** — `IMPACT_MATRIX.md`
dodaje semantykę (severity + edge type + hotspot link), której
auto-generated graph nie ma.

---

## 2. Coarse view (top-level)

(import z `DEPENDENCIES.md` § "Coarse view")

| Source | Targets (depends-on) | Notes |
|--------|----------------------|-------|
| `<root>` (main.tex) | `axioms`, `core`, `core/_meta_latex`, `core/formalizm`, `papers_external`, `partial_proofs` | Manuscript root |
| `axioms` | `core`, `core/formalizm`, `partial_proofs` | A1A axiom v2 + slownik |
| `core` | `axioms`, `core/_meta_latex`, `core/formalizm`, `partial_proofs`, `research` | Sek01–Sek10 + dodatki |
| `core/_meta_latex` | `axioms`, `core`, `core/formalizm`, `partial_proofs` | Mapa statusu |
| `core/formalizm` | `axioms`, `core`, `partial_proofs` | Dodatki formalizujące |
| `meta` | -- | Audyty / plany (no LaTeX deps) |
| `papers_external` | `core`, `core/formalizm` | External response papers |
| `partial_proofs` | `axioms`, `core`, `core/formalizm` | Częściowe dowody |
| `research` | varies (~85 folderów) | Active research workstreams |

---

## 3. Fine view — historyczne ścieżki closure z audytu 2026-05-01

> **Aktualizacja Sesji 3.5 (2026-05-03):** wszystkie 20 hotspotów z
> `CORE_HOTSPOTS.md` zostały zamknięte w samym audycie 2026-05-01
> (zob. [[meta/research/HOTSPOT_AUDIT_S3_5.md]]). Krawędzie poniżej są
> więc **historical paths of closure**, nie aktywne `needs-migration`.
> Severity zmieniona na `HIST` (historical) — folder już ma audit-aware
> markery (`B6-CLOSED`, `B8-CLOSED`, `B9-CLOSED` w wynikach M9.x).

| Source folder | → Target (core file lub research) | Severity | Edge type | Hotspot (zamknięty) | Closure source |
|---------------|------------------------------------|---------:|-----------|----------------------|-----------------|
| `research/op-newton-momentum` (B6_m9x_sqrtg_rerun.py 6/6 PASS) | → `core/sek08c_metryka_z_substratu.tex` | HIST | data | H-A2, H-B6, H-B10 | § U |
| `research/op-newton-momentum` (B8_lagrangean_independence_check.py 5/5 PASS) | → `core/sek08c` (M9.1'' Lagrangean) | HIST | lemma | H-B8 | § V |
| `research/op-newton-momentum` (B9 WEP MICROSCOPE 6/6 PASS) | → `core/sek06_czarne_dziury` | HIST | data | H-B9 | § W |
| `research/op-tau3-substrate-clock-acceleration` (Phase 4 Adams-positivity 5/5 PASS) | → `core/sek10_N0_wyprowadzenie` (Adams) | HIST | methodology | H-B1 | § Z |
| `research/op-tau3-substrate-clock-acceleration` (B7_greens_function 2-regime sympy LOCK) | → ψ.1.v2 / Phase 5 corrections | HIST | data | H-B7 | § O.6 + § T |
| `research/op-tau3-substrate-clock-acceleration` (multiplicative m_e_eff patch) | → registry annotations TT7-TT12 | HIST | methodology | H-A5 | § L |
| `research/op-psi1-substrate-light-acceleration` (Phase1-3 INVALIDATED markers) | → registry [WITHDRAWN] TT13-TT18 | HIST | annotation | H-A6, H-A8 | § K |
| (B5 q-dimension reconciliation across LaTeX core) | → multiple core files | HIST | data | H-B5 | § X |
| (B3-v2 α_s 0.1184 propagation 14+ files) | → README/sek00/companion/scripts | HIST | data | H-B3 | § Y |
| `research/closure_2026-04-26/sigma_ab_pathB` (cross-validation 11/11 PASS) | → `core/sek08_formalizm` (alternative) | CROSS-VAL | lemma | H-B6 (alt track) | closure_2026-04-26 |
| `research/closure_2026-04-26/f_psi_principle` (cross-validation 12/12 PASS) | → `core/sek08c` (alternative) | CROSS-VAL | lemma | H-B8 (alt track) | closure_2026-04-26 |
| `research/closure_2026-04-26/Lambda_from_Phi0` (7/7 PASS) | → `core/sek05_ciemna_energia` (Ω_Λ vacuum cat.) | NEW | data | (nowy strukturalny wkład; nie w 43-list) | closure_2026-04-26 |
| `research/closure_2026-04-26/alpha_psi_threshold` (cross-validation 5/5 PASS) | → `core/sek06_czarne_dziury` (alternative) | CROSS-VAL | lemma | H-B9 (alt track) | closure_2026-04-26 |
| `research/why_n3` (5-phase emergent Dirac propagator) | → cross-deepen H-A1+H-A2 (R3↔M9.1'' horizon) | CROSS-VAL | structural | H-A1, H-A2 (deepening) | § AB |

---

## 4. Fine view — krawędzie inter-research (Sesja 6 wypełnione 2026-05-03)

> **4 PROPOSED bridges** z heurystycznego matchingu w Sesji 6 (zob.
> [[meta/research/CANDIDATE_BRIDGES.md]] §3). Konserwatywna heurystyka —
> 86 folderów × 86 = ~7400 par, z czego tylko 4 przeszły anti-spam threshold.

| ID | Source folder | → Target folder | Severity | Edge type | Strength | Common signals |
|----|---------------|------------------|---------:|-----------|----------|----------------|
| BR-001 | `op-delta1-g-tilde-derivation` (NEEDS) | ← `op-uv3-phi0-renormalization` (FINDINGS) | LOW | data | HEURISTIC (score=10) | tokens: n_c, sek00, α_s, φ_eff, ω_λ |
| BR-002 | `op-gamma1-phi-eff-anchor-resolution` (NEEDS) | ← `op-uv3-phi0-renormalization` (FINDINGS) | LOW | data | HEURISTIC (score=6) | tokens: sek00, α_s, φ_eff |
| BR-003 | `atomic_shells_closure` (NEEDS) | ← `qm_superposition` (FINDINGS) | LOW | methodology | PARTIAL (score=5) | explicit folder cross-reference |
| BR-004 | `mass_scaling_k4` (NEEDS) | ← `op-lambda1-e2-amplitude-emergence` (FINDINGS) | LOW | data | PARTIAL (score=5) | explicit folder cross-reference |

**Implikacja:** UV.3 (Z_Φ = 14/3 derivation) jest source dla 2 folderów
(`δ.1`, `γ.1`) — to wzmacnia kandydaturę UV.3 na core (zob.
INTAKE candidate dla UV.3 w przyszłej Sesji 7+).

**45 orphan NEEDS** (folderów bez bridge candidate) i **67 orphan FINDINGS**
(eksportujących bez konsumenta) → szczegóły w `CANDIDATE_BRIDGES.md` §5/§6.

---

## 5. Regression triggers (Sesja 8)

> **Severity legenda zaktualizowana po Sesji 3.5:**
> - `CRITICAL` / `HIGH` / `MEDIUM` / `LOW` — *aktywne* krawędzie (audit-driven jeśli odnoszą się do nowego hotspota)
> - `HIST` — historical path of closure (hotspot już zamknięty audytem 2026-05-01); folder ma audit-aware markery
> - `CROSS-VAL` — alternative formal track (np. closure_2026-04-26 vs § U/V/W); wymaga decyzji człowieka czy core ma cytować obie ścieżki
> - `NEW` — strukturalny wkład poza 43-item list (np. Λ_from_Phi0, why_n3 cross-deepening)

Skrypt-checker w Sesji 8 sprawdza dla każdej krawędzi:

- Source `mtime_newest > Target last_reviewed` ⇒ stale-dependency (severity zachowana)
- Source `folder_status: needs-migration` && Target `core_compatibility: current` ⇒ forward-incompatibility (CRITICAL)
- Source w kwarantannie 74394a8 && jakikolwiek Target ⇒ **CRITICAL: poisoned-edge**
- `HIST` krawędzie ⇒ INFO only (nie CRITICAL); sprawdzane czy folder źródłowy ma `B*-CLOSED` markery zgodne z audytem
- `CROSS-VAL` krawędzie ⇒ INFO only; sprawdzane czy `INTAKE_*.md` istnieje w `meta/core/intake/` (jeśli jest taka decyzja człowieka)
- `NEW` krawędzie ⇒ HIGH (nowe strukturalne wkłady wymagają jawnej oceny czy zostają jako stand-alone czy promowane do core)

Wynik trafia do `meta/research/REGRESSION_S8.md`.

---

## 6. Synchronizacja z innymi plikami

- Każdy folder z `tgp_status.depends_on: [...]` ⇒ wiersze w § 4.
- Każdy folder z `tgp_status.impacts: [...]` ⇒ wiersze w § 4 (odwrócona
  perspektywa — sprawdzamy spójność).
- Każdy hotspot z `CORE_HOTSPOTS.md` z `Linked needs-migration folders`
  ⇒ wiersze w § 3.
- `EXECUTED` bridges z `CANDIDATE_BRIDGES.md` ⇒ wiersze w § 4.

---

## 7. Co Sesja 6+7 doda

- **Sesja 6** (research bus + bridges): krawędzie `inter-research` (§ 4),
  ~30–80 wpisów.
- **Sesja 7** (core-ready/promoted/migration): krawędzie `research → core`
  (§ 3), ~15–25 wpisów na bazie `CORE_HOTSPOTS.md`.
- Po Sesji 7: pełna macierz ~50–100 krawędzi z severity i hotspot link.
