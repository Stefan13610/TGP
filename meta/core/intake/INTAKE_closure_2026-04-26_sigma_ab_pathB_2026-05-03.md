---
title: "INTAKE — closure_2026-04-26/sigma_ab_pathB → core (alternative track for σ_ab dynamics)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/closure_2026-04-26/sigma_ab_pathB
target_section: core/sek08_formalizm (option: dodać alternative derivation block)
hotspot_addressed: H-B6 (alternative track; main closure: § U)
agent_requesting: Sesja 3.5 Q6 (agent)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/closure_2026-04-26/sigma_ab_pathB/results.md]]"
  - "[[research/closure_2026-04-26/sigma_ab_pathB/setup.md]]"
  - "[[research/closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
tags:
  - intake
  - core-promotion
  - cross-validation
  - sigma-ab
  - path-B-primary
---

# INTAKE — `closure_2026-04-26/sigma_ab_pathB` → core (alternative track for σ_ab dynamics)

## 1. Wynik proponowany do promocji

**σ_ab Path B PRIMARY**: tensor σ_ab dynamika jest **operatorową konsekwencją**
ŝ-EOM (heredity construction), NIE niezależny Lagrangian L_σ. Trzy
strukturalne wnioski:

1. **`□σ_ab + 2 m_s² σ_ab = source + grad-coupling`** — heredity equation
   derivowana krok-po-kroku z `□δŝ + m_s² δŝ = J` (sympy exact zero residual)
2. **`m_σ² = 2 m_s²`** — derived (NIE postulated) z box-of-product algebra;
   spectral threshold `√s_min = 2 m_s` z OPE
3. **Ghost-free przez Gram-positivity** — strukturalnie z konstrukcji
   (`K_ab ≥ 0` jako Gram matrix); single-Φ Z₂ aksjomat zachowany

## 2. Source w research/

- `research/closure_2026-04-26/sigma_ab_pathB/results.md` (10.4 KB) —
  TL;DR + 5 sub-tests T-PB.1 do T-PB.5 (11/11 PASS)
- `research/closure_2026-04-26/sigma_ab_pathB/sigma_ab_pathB_audit.py`
  (16 KB) + `.txt` (5.6 KB) — sympy + numerical verification
- `research/closure_2026-04-26/sigma_ab_pathB/setup.md` (8.3 KB) —
  audit design (5 kryteriów)
- `research/closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md` — Phase 1
  4-phase closure aggregator (35/35 PASS total)

## 3. Cross-validation z audytem 2026-05-01

**KLUCZOWY KONTEKST:** Audyt 2026-05-01 § U **niezależnie** zamknął ten
sam strukturalny problem (H-B6) przez:

- `research/op-newton-momentum/B6_m9x_sqrtg_rerun.py` (6/6 PASS)
- M9.x full re-run z poprawnym `√(-g) = c·ψ/(4-3ψ)`
- Path: `√(-g)` correction + scalar/tensor mode ratio z δΦ

**Dwie niezależne ścieżki zamknięcia tego samego hotspota** (H-B6):

| Ścieżka | Folder | PASS | Filozofia |
|---------|--------|-----:|-----------|
| Audyt main (§ U) | `op-newton-momentum/B6_m9x_sqrtg_rerun.py` | 6/6 | √(-g) re-run z scalar/tensor mode ratio |
| Cross-validation | `closure_2026-04-26/sigma_ab_pathB/sigma_ab_pathB_audit.py` | 11/11 | Path B heredity z ŝ-EOM (single-Φ axiom preserved) |

To jest **zdrowy znak** matematyczny — dwa niezależne tracki dochodzą
do tego samego strukturalnego rezultatu (`m_σ² = 2 m_s²` ghost-free).

## 4. Proponowana lokacja w core

| Target | Co dokładnie ma się tam pojawić | Czy istnieje teraz | Diff (dla istniejących sekcji) |
|--------|----------------------------------|---------------------|---------------------------------|
| `core/sek08_formalizm.tex` | Nowa subsection: "σ_ab as Path B heredity construction" — alternative derivation z `□δŝ + m_s² δŝ = J` | NO (tylko § U cytat) | (nowa subsection ~30 linii LaTeX) |
| `core/sek08_formalizm.tex` przy `ax:metric-coupling` | Cross-link do Path B preserving single-Φ axiom | YES (lin. 11132+) | dopisać 1-2 linie cross-reference |

## 5. Hotspot, który ten wynik adresuje

**H-B6** w [[meta/core/CORE_HOTSPOTS.md]] §B (status: `RESOLVED-AUDIT-FULL`
przez § U). Ten INTAKE proponuje **dodanie alternative track** —
core cytowałby obie ścieżki jako mutual cross-validation.

Cytat hotspota (oryginalny audyt § B.6): "No-graviton claim vs M9.3 daje
TYLKO scalar mode `h_b=h_L=4δψ` z fundamentalnej linearyzacji δΦ; tensor
h_+, h_× wymagają osobnego pola σ_ab (ŁAMIE single-Φ axiom)".

Path B (sigma_ab_pathB) **explicite preserves** single-Φ axiom — σ_ab to
operator pochodny, nie nowe pole. To jest **dokładnie** ta charakterystyka,
której hotspot wymagał.

## 6. Co MUSI być sprawdzone przed promocją

- [ ] Wynik kompiluje się (sympy/numerical) na bieżącej wersji core
- [ ] Brak konfliktu z aktualnym `main.tex` (read-only sprawdzenie)
- [ ] `core_compatibility: current` w YAML folderu źródłowego (do uzupełnienia w Sesji 4)
- [ ] Wszystkie 11 testów T-PB.1–T-PB.5 PASS (potwierdzone w `results.md`)
- [ ] `source_of_status` w YAMLu zawiera linki do plików source (do uzupełnienia w Sesji 4)
- [ ] Ścieżki § U (`B6_m9x_sqrtg_rerun.py`) i Path B (`sigma_ab_pathB_audit.py`)
      dają zgodne predykcje dla `m_σ²`, scalar/tensor ratio (cross-validation
      sanity check)

## 7. Co się STANIE po akceptacji człowieka

1. **Człowiek** edytuje `core/sek08_formalizm.tex` ręcznie, dodając
   subsection "σ_ab as Path B heredity construction"
2. **Człowiek** dopisuje wpis do [[meta/core/CORE_INVENTORY.md]] §1.10
3. **Człowiek** zmienia w `research/closure_2026-04-26/README.md`:
   - `tgp_status.subfolder_summary[0].promoted: true` (lub równoważne)
   - `last_yaml_update: <data promocji>`
   - W `closure_2026-04-26/sigma_ab_pathB/` dodaje notatkę post-promotion
4. **Człowiek** zmienia status hotspota H-B6 w `CORE_HOTSPOTS.md`
   z `RESOLVED-AUDIT-FULL` na `RESOLVED-AUDIT-FULL + RESOLVED-FULL-PROMOTED`
   (dual closure z dwoma ścieżkami w core)
5. **Człowiek** zamyka ten wniosek INTAKE: ustawia
   `status: ACCEPTED-AND-PROMOTED` + dopisuje datę

## 8. Decyzja człowieka

- [ ] **ACCEPT** — Path B promowany jako alternative track w core (data: ____)
- [ ] **DEFER** — wynik OK, ale promocja przesunięta (powód: ____)
- [ ] **REJECT** — § U closure jest wystarczający; alternative track nie wnosi (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — wymaga dopełnienia (lista: ____)

---

## 9. Notatki dla decyzji

### 9.1 Argumenty PRO

1. **Cross-validation matematyczna** — dwa niezależne formalne tracki
   dające zgodne wyniki to silny sygnał poprawności
2. **Path B explicite zachowuje single-Φ axiom** — hotspot H-B6 oryginalnie
   zwracał uwagę na "ŁAMIE single-Φ axiom"; Path B to fix-it-by-construction
3. **Algebraiczna derywacja `m_σ² = 2 m_s²`** zamiast post-hoc analogii
   mezonowej (improvement over OP-7 T3.1)

### 9.2 Argumenty CONTRA / risk factors

1. **Audyt main § U już closurował hotspot z markerami w core** — dodanie
   drugiej ścieżki może być **redundancja**, nie value-add
2. **Path B introduces grad-coupling term** którego § U nie ma —
   verification, że oba dają **identyczne** predykcje obserwacyjne, jest
   non-trivial
3. **`closure_2026-04-26/` jest sub-folder closure-aggregator** — promocja
   jednego dziecka może być inkoherentna z resztą rodzeństwa

### 9.3 Rekomendacja agent

**ACCEPT z modyfikacją:** dodać subsection w `sek08_formalizm.tex`,
ale jako **note/remark** ("Path B alternative derivation, preserving
single-Φ axiom"), nie jako pełny replacement § U logiki. Daje to
mathematical depth bez konfliktu z istniejącym closurem audytu.

To jest tylko opinia agenta. Decyzja należy do człowieka.
