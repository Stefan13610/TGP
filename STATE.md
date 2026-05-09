---
title: "STATE.md — TGP_v1 single-source coordination point"
date: 2026-05-09
type: state
status: ACTIVE
purpose: "Jedyny plik aktualizowany po każdej sesji. Inne warstwy (INDEX, audyt/PRIORITY_MATRIX, meta/PLAN_*) są referencyjne."
update_policy: "Aktualizować po: (a) closure cyklu, (b) zmianie krytycznej ścieżki, (c) zmianie WIP."
---

# STATE.md — current state of TGP_v1 framework

> **Po co ten plik?** Single-source-of-truth dla "co się dzieje teraz".
> Diagnoza 2026-05-09: 80 cykli z `folder_status: active` w README ≠ realnie WIP.
> Bez WIP-limitu i centralnego entry-point każda sesja zaczyna się od audytu stanu.
>
> **Reguła:** ten plik aktualizować po każdej sesji. INDEX.md, audyt/PRIORITY_MATRIX,
> meta/PLAN_* zostają, ale są referencyjne — nie są źródłem prawdy o aktualnym WIP.

## 🔴 Critical path

| Cykl | Faza | Blocker dla | Owner |
|---|---|---|---|
| **[[research/op-S07-alternative-f-psi-derivation-2026-05-09/]]** | Phase 0 (setup, ongoing) | **Cały sektor grawitacji** TGP — M911-P1/P2/P3 w limbo | next session |

**Dlaczego krytyczne:** M9.1'' (4−3ψ)/ψ FALSIFIED 5.02σ przez GWTC-3 RERUN
([[research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]).
Wszystkie predykcje sektora grawitacji (Newton limit, PPN, GW phase, BH ringdown)
zależą od specific f(ψ). Dopóki S07 nie znajdzie alternative spełniającej C1-C10
z first-principles derivation — TGP nie produkuje falsyfikowalnych predykcji
grawitacyjnych.

**Honest probability assessment** (z [[research/op-S07-alternative-f-psi-derivation-2026-05-09/README.md]]):
- 25–35% alternative f(ψ) found + first-principles derived
- 30–40% alternative empirical only (constraint-satisfying ale nie first-principles)
- 20–30% NO alternative → M9.1'' approach FALSIFIED structurally → radical gravity rethink

## 🟡 Active WIP (limit: 5 równolegle)

Cykle które realnie poruszają się w tej i następnej sesji. **Critical path
(S07) ma zarezerwowany slot, nie liczy się do limitu 5** (per
[[meta/CYCLE_LIFECYCLE.md]] §WIP-limit).

| # | Cykl | Faza / status | Następny krok |
|---|---|---|---|
| ★ | **[[research/op-S07-alternative-f-psi-derivation-2026-05-09/]]** | Phase 0+ (critical path) | C1-C10 sympy formalization + Phase 1 candidate enumeration |
| 1 | [[research/op-FRW-radiation-era-varying-c-2026-05-06/]] | Phase 2 PASS, ścieżka A FAILS | decyzja D/E/F (pivot L_mat?) |
| 2 | [[research/op-Phi-decomposition-photon-2026-05-07/]] | aktywny | kontynuacja dekompozycji Φ → fotony (V-independent) |
| 3 | (slot wolny) | — | — |
| 4 | (slot wolny) | — | — |
| 5 | (slot wolny) | — | — |

> **Korekta WIP z 2026-05-09 wieczór:** `op-MAG-anomalous-moment-2026-05-09` był początkowo na liście WIP-5, ale jego YAML ma `status: EARLY_HALT_2026-05-09` (sympy 2/2 PASS, classification `EARLY_HALT_HONEST`) — czyli już zamknięty z honest acknowledgment. Reklasyfikowany na `closed-NULL`, zwolnił WIP slot. Nie ma silnego kandydata na zastępcę z reszty Bucket A — uczciwiej zostawić 2 wolne sloty niż wpychać słabego kandydata.
>
> **Korekta WIP z 2026-05-09 noc:** `op-emergent-metric-from-interaction-2026-05-09` zamknięty przez parallel agent (Phase 1-6 complete, **57/57 sympy PASS, STRUCTURAL_DERIVED**). 6/6 wymagań P1-P6 RESOLVED, 13/14 NEEDS resolved. Reklasyfikowany na `closed-resolved`, zwolnił kolejny WIP slot. Wynik **bezpośrednio relevantny dla S07**: g_eff = G[{Φ_i}] proposal może być fundamentem alternative f(ψ) (interaction-emergent zamiast postulate-functional).

**Co poszło do `paused`** (z poprzedniej listy / Bucket A):

- `op-D01-anchor-lock-2026-05-06` — strukturalny audit, można wznowić po S07
- `audyt/T01_LIGO3G_falsifier/` — FALSIFIER_STATEMENT_DRAFT czeka na M9.1'' resolution z S07; trzymany jako meta-debt #3 (osobny handoff w środku — do skonsolidowania)
- `papers/M911_LIGO3G_paper/` — paper draft na pause do post-S07 (M9.1'' status pending)

**Reguła WIP:** maksymalnie 5 cykli `active` (poza critical-path slot) w jednym
czasie. Wszystkie inne oznaczone w `folder_status` jako jeden z: `paused`,
`needs-bridge`, `parking`, `closed-resolved`, `closed-NULL`,
`closed-superseded`. Pełna polityka: [[meta/CYCLE_LIFECYCLE.md]].

## ✅ Recent closures (last 5–7)

Wszystkie 2026-05-09:

| Cykl | Sympy | Verdict |
|---|---|---|
| [[research/op-Phi-vacuum-scale-2026-05-09/]] | 84/88 (95.5%) | STRUCTURAL_DERIVED_CONDITIONAL_HALT |
| [[research/op-V-canonical-consistency-audit-2026-05-09/]] | 10/10 | dual-V framework confirmed |
| [[research/op-MAG-Phase5-V-reference-clarification-2026-05-09/]] | 10/10 | erratum applied |
| [[research/op-dual-V-structure-clarification-2026-05-09/]] | 10/10 | TGP_FOUNDATIONS §3.5 added |
| [[research/op-Phase5-MAG-erratum-2026-05-09/]] | 5/5 | γ = m_C² correction |
| [[research/op-Phi0-spatial-variation-predictions-2026-05-09/]] | 6/6 | atomic clocks + EP predictions logged |

**Cumulative:** 103/107 PASS (96.3%) z dual-V audit chain.

## ⚠ Outstanding meta-debt

Sygnał że framework wymaga porządków obok pracy badawczej.

### Załatwione 2026-05-09 (post-cleanup)

| # | Dług | Status | Co zrobiono |
|---|---|---|---|
| ~~1~~ | INDEX.md stale (2026-04-28) | ✅ **DONE 2026-05-09** | Dodano banner critical-blocker S07 + STATE.md jako primary entry-point + audyt/, CYCLE_LIFECYCLE, CALIBRATION_PROTOCOL w Top-level entry points; date 2026-04-28 → 2026-05-09 |
| ~~2~~ | DEPENDENCIES.md stale (2026-04-22) | ✅ **DONE 2026-05-09** | Regenerated via `tooling/build_deps_graph.py`: 117 tex / 1098 md / 70 inputs / 1469 refs / 5891 wikilinks (z ~1657 dependencies poprzednio — ×4 wzrost) |
| ~~3~~ | Drugi handoff w audyt/T01 | ✅ **DONE 2026-05-09** | Zarchiwizowany jako stub: [[audyt/T01_LIGO3G_falsifier/HANDOFF_PROMPT_NEXT_SESSION.md]] (treść była pre-falsification, β=−5/64 ; faktycznie po RERUN β=−15/4 → TGP RULED OUT 5σ). T01 paused do post-S07 |
| ~~4~~ | 80 cykli z `folder_status: active` (realnie ~5) | ✅ **DONE 2026-05-09** | Mass-triage: 85 → `paused` (auto), 9 → `closed-resolved` (cascade), 1 → `closed-NULL` (MAG-anomalous), 4 → manual fix (M03/L01/L04/void-flat-modes), 2 → `parking` (SPIN-MAG-leakage, tensor-modes-FUTURE). Patrz commit `67e0677` |
| ~~5~~ | Brak cycle-lifecycle policy | ✅ **DONE 2026-05-09** | Spisane: [[meta/CYCLE_LIFECYCLE.md]] (9 statusów, WIP-limit, anti-patterns, mapping legacy) |

### Załatwione 2026-05-09 (post-cleanup, runda 6-10)

| # | Dług | Status | Co zrobiono |
|---|---|---|---|
| ~~6~~ | LaTeX cruft committed historycznie | ✅ **FALSE ALARM 2026-05-09** | `git ls-files \| grep -E '\.(aux\|log\|bbl\|...)$'` zwrócił 0 wyników. Pliki NIGDY nie były tracked — .gitignore działa od początku. Lokalne build artifacts pozostają tylko na dysku |
| ~~7~~ | 3 PDF kanoniczne? | ✅ **DOCUMENTED 2026-05-09** | Spisane w [[PAPER_LAYOUT.md]]: main.pdf=full PL thesis (autorska), tgp_letter.pdf=PRL English (krótki submission), tgp_companion.pdf=PRD English (długi technical). Trójdzielny layout standardowy. Decyzja "który kanoniczny" zależy od kontekstu — patrz tabela w PAPER_LAYOUT.md |
| ~~8~~ | Documentation drift `status` ↔ `folder_status` | ✅ **TOOLING + 2 manual fixes 2026-05-09** | Skrypt detekcji: [[tooling/check_status_drift.py]] (read-only). Zastosowane 2 oczywiste fixy: op-g0-r3-from-canonical-projection (paused → closed-resolved, text "PHASE 4 CLOSED-POSITIVE"), op-omicron2-phi-mean-shift-cosmo (paused → closed-NULL, text "STAGE_1_NULL_CLOSED_2026-05-03"). Pozostałe drifty pozostają — `folder_status` jest source of truth, text status — manual fix per cykl |
| ~~9~~ | Brak skryptu auto-pause stale cycles | ✅ **DONE 2026-05-09** | Spisane: [[tooling/check_stale_cycles.py]] (read-only weekly report). Domyślny próg 30 dni, `--strict` daje 14 dni. Exit code 1 jeśli znaleziono stale-active (do CI/cron) |
| ~~10~~ | DEPENDENCIES_REVERSE.md duplikat | ✅ **NO ACTION 2026-05-09** | Świadoma decyzja: zostawić (`tooling/build_deps_graph.py` generuje oba). Niskoryzyko duplicate, czasem przydatny dla "kto cytuje X". Można usunąć w przyszłości jeśli nigdy się nie używa |

### Otwarte (do osobnych sesji)

| # | Dług | Lokalizacja | Notatka |
|---|---|---|---|
| 11 | Text status drift w ~15 cyklach (post #8 partial fix) | `research/op-*/README.md` | Skrypt `check_status_drift.py` raportuje pozostałe; manual review każdego rzędu (np. dual-V cascade dalej "PHASE0_PHASE1_IN_PROGRESS") |
| 12 | `tgp_companionNotes.bib` + `tgp_letterNotes.bib` placeholders | root | ~104 bytes każdy. Letter+companion używają `tgp_main.bib` — te pliki to noise. Niskoryzyko cleanup |
| 13 | INDEX.md cycle-list nieaktualne | INDEX.md | Zaktualizowane top-level entry points + critical-blocker banner, ale long Phase ledger + per-cykl listy mogą wymagać regen po quartet of 2026-05-09 closures |

## 🗂 Coordination layers — co czym jest

Żeby uniknąć duplikatów i drift'u:

| Plik | Rola | Aktualizacja |
|---|---|---|
| **STATE.md** (TEN) | Critical path + WIP + recent closures + meta-debt | Po każdej sesji |
| [[INDEX.md]] | Indeks plików / głęboka nawigacja | Co kilka tygodni; obecnie stale |
| [[README.md]] | Entry point dla nowych — filozofia + high-level | Rzadko; stabilny |
| [[TGP_FOUNDATIONS.md]] | Aksjomatyczna referencja (W/E/P/H, dual-V §3.5) | Przy zmianach strukturalnych |
| [[PREDICTIONS_REGISTRY.md]] | Wszystkie predykcje (FALSIFIED/PASS/PENDING) | Po każdym Phase 4-5 closure |
| [[DEPENDENCIES.md]] | Auto-generated graph zależności | `tooling/build_deps_graph.py` |
| [[audyt/README.md]] + [[audyt/PRIORITY_MATRIX.md]] | Strukturalne długi (S/L/D/M/T/EXT) | Po każdym audit closure |
| `meta/PLAN_*`, `meta/CALIBRATION_PROTOCOL.md` | Procedury i meta-zasady | Rzadko; stabilne |

**Zasada:** STATE.md wskazuje JEDNĄ rzecz krytyczną + max 5 WIP. Reszta to zasoby
referencyjne. Nie kopiować ich treści tutaj.

## 📋 WIP lifecycle (proposal — nie wdrożone strukturalnie)

Reguła kiedy cykl wchodzi w jaki status (do przepisania w `meta/CYCLE_LIFECYCLE.md`
w odpowiedniej sesji):

| Status | Warunek wejścia | Warunek wyjścia |
|---|---|---|
| `active` | Wybrany na critical path lub WIP slot wolny | Phase FINAL closed lub pivot do `paused` |
| `paused` | Świadomie zamrożony; blocker udokumentowany w README | Blocker rozwiązany → `active` |
| `needs-bridge` | Czeka na poprzednika (op-X CLOSED dependency) | Poprzednik CLOSED → `active` |
| `parking` | Pomysł zarejestrowany, niegotowy do startu | User decyzja → `active` |
| `closed-resolved` | Phase FINAL z verdict DERIVED/STRUCTURAL_CONDITIONAL | — |
| `closed-NULL` | Phase FINAL z verdict EARLY_HALT honest | — |
| `closed-superseded` | Inny cykl objął zakres | Link do następcy w README |
| (auto-pause) | Brak commita >30 dni | — (wymaga skryptu `tooling/auto_pause_stale.py`) |

## 📜 Migration log

| Data | Zmiana |
|---|---|
| 2026-05-09 | STATE.md utworzony jako single-source coordination point |
| 2026-05-09 | Handoff `HANDOFF_NEXT_SESSION_S07_alternative_f_psi.md` (root) → migrated do `op-S07-alternative-f-psi-derivation-2026-05-09/`; root file zamieniony na stub |
| 2026-05-09 | Cycle `op-S07-alternative-f-psi-derivation-2026-05-09` otwarty (Phase 0) |
| 2026-05-09 | `meta/CYCLE_LIFECYCLE.md` policy spisana (dwa poziomy statusu, WIP-limit, słownik 9 statusów, anti-patterns) |
| 2026-05-09 | Inwentaryzacja 116 cykli `research/`: A=19 active-recent, B=3 mislabeled-closed, C=91 stale-active, D=6 needs-bridge, E=10 unknown |
| 2026-05-09 | WIP-5 selected: S07 (★) + FRW + emergent-metric + MAG-anomalous + Phi-decomposition-photon. D01 + audyt-T01 + M911-paper → paused/meta-debt |
| 2026-05-09 | `tooling/reclassify_cycles_2026-05-09.py` script (mass-triage Bucket A+B+C, dry-run domyślnie) |
| 2026-05-09 | Mass-triage applied: 85 cykli `active`/`research` → `paused` (auto via skrypt) |
| 2026-05-09 | Manual fix 4: M03/L01/L04 → `closed-resolved`; void-flat-modes naming `closed_NULL` → `closed-NULL` |
| 2026-05-09 | 15 edge cases bez `folder_status` field — dodane top-level: 3× `active` (S07, emergent-metric, Phi-decomposition-photon), 9× `closed-resolved` (Phi-vacuum + dual-V cascade + MAG-Lorentz/resonance, SPIN-SU2), 1× `closed-NULL` (MAG-anomalous EARLY_HALT odkryte przy edycji), 2× `parking` (SPIN-MAG-leakage informal, tensor-modes-FUTURE placeholder) |
| 2026-05-09 | **Documentation drift wykryty:** 5 cykli z dual-V cascade ma w README `status: PHASE0_PHASE1_IN_PROGRESS` mimo że parent `op-Phi-vacuum-scale/Phase_FINAL_close.md` dokumentuje je jako zamknięte. Tekstowy `status:` field nie został zaktualizowany przy cascade closure 2026-05-09. `folder_status: closed-resolved` dodane na podstawie parent's claim — text status do osobnego cleanupu |
| 2026-05-09 | **Outstanding-debt #1-#5 załatwione:** INDEX.md update (banner S07 + STATE.md primary entry-point + audyt/CYCLE/CALIBRATION w entry points), DEPENDENCIES.md regenerated (×4 wzrost dependencies), audyt/T01 HANDOFF zarchiwizowany jako stub (pre-falsification, β=−5/64 stale), #4+#5 oznaczone DONE (mass-triage + CYCLE_LIFECYCLE policy z poprzednich rund) |
| 2026-05-09 | **Outstanding-debt #6-#10 załatwione:** #6 false alarm (LaTeX cruft nigdy nie tracked), #7 PAPER_LAYOUT.md (3 PDF role spisane), #8 check_status_drift.py + 2 manual fixes (g0-r3 → closed-resolved, omicron2 → closed-NULL), #9 check_stale_cycles.py, #10 no action (świadomie) |
| 2026-05-09 | **op-emergent-metric-from-interaction CLOSED:** parallel agent zamknął cykl (Phase 1-6 complete, 57/57 sympy PASS, **STRUCTURAL_DERIVED**). Bezpośrednio relevantny dla S07 — g_eff = G[{Φ_i}] może być fundamentem alternative f(ψ) (interaction-emergent zamiast postulate-functional). WIP-5 zwolniło 2 sloty (z poprzedniego MAG-anomalous EARLY_HALT discovery + emergent-metric closure) |
