---
title: "Sesja 1 — Audyt struktury research/"
date: 2026-05-03
type: audit
status: COMPLETED
session: S1
generated_at: "2026-05-03T12:56:28"
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
  - "[[INDEX.md]]"
  - "[[DEPENDENCIES.md]]"
tags:
  - audit
  - research-workflow
  - session-1
---

# Sesja 1 — Audyt struktury `research/`

> **Read-only audyt.** Nie nadaje statusów (to Sesja 4), nie modyfikuje folderów,
> nie przenosi plików. Wynik jest wejściem dla Sesji 2 (szablony) i Sesji 4 (klasyfikacja).

- Canonical: `research`
- Liczba folderów badawczych: **86**
- Generated: `2026-05-03T12:56:28`
- Skrypty audytu: `meta/research/_audit_s1.py` + `meta/research/_audit_s1_summarize.py`
- Surowe dane: `meta/research/_audit_s1_raw.json` (176 007 B)

## 1. Agregaty

| Metryka | Liczba | % z 86 |
|---|---:|---:|
| `has_program_md` | 32 | 37% |
| `has_readme` | 33 | 38% |
| `has_findings` | 0 | 0% |
| `has_needs` | 0 | 0% |
| `readme_has_frontmatter` | 8 | 9% |
| `readme_has_tgp_status` | 0 | 0% |
| `is_polluted_74394a8` | 4 | 4% |
| `is_legacy_stay` | 3 | 3% |
| folders with subfolders | 13 | 15% |
| `FULL CONVERGENCE` occurrences (anti-overclaim flag) | 142 | — |

**Kluczowe obserwacje agregatów:**

- 0/86 folderów ma już `FINDINGS.md` (cel Sesji 5: 100%).
- 0/86 folderów ma już `NEEDS.md` (cel Sesji 5: 100%).
- 0/86 folderów ma już `tgp_status:` w README (cel Sesji 4: 100%).
- 32/86 folderów ma `program.md` (kanoniczna konwencja 3-fazowa).
- `FULL CONVERGENCE`: 142 wystąpień — kandydaty do weryfikacji anty-overclaim w Sesji 8.

> **Nota o liczbie folderów:** plan szacował "~95". Audyt liczy
> faktycznie **86** podfolderów (różnica = top-level pliki w `research/`,
> sekcja 10: programowe dokumenty + grafy). Liczba 86 jest precyzyjna.

## 2. Heurystyczna pre-klasyfikacja (NIE jest werdyktem)

> **WAŻNE:** Heurystyka skanuje keywordy (`PASS`, `CLOSED`, `FALSIFIED`, …) i nazwy.
> Nie jest decyzją statusową. Sesja 4 nadaje status z weryfikowalnym źródłem.
>
> **Znane ograniczenia heurystyki (kalibracja 2026-05-03):**
>
> - Słowo `FALSIFIED` w aktywnym tekście oznacza zazwyczaj **regułę**
>   falsyfikacji predykcji (np. "if observed → TGP falsified"), NIE status
>   archiwalny folderu. Dlatego archive-flag wymaga `WITHDRAWN`/`OBSOLETE`/
>   `SUPERSEDED` (te są używane chirurgicznie) lub `FALSIFIED` przy bardzo
>   niskim PASS.
> - `candidate-active (heuristic)` jest klasą "szumową" — pokrywa wszystko,
>   co nie pasuje do bardziej specyficznych klas. Sesja 4 musi je faktycznie
>   sklasyfikować.
> - `candidate-needs-bridge-or-sandbox` to topic-folders z minimalnym
>   scaffoldingiem (1–2 .md, ≤3 .py). Część z nich jest realnie active,
>   tylko jeszcze nie ma `program.md`.

| Wstępna klasa heurystyczna | Liczba |
|---|---:|
| candidate-active (heuristic) | 25 |
| candidate-active (3-phase op-*) | 25 |
| candidate-needs-bridge-or-sandbox (heuristic) | 19 |
| candidate-archive (heuristic, verify) | 6 |
| needs-bridge (polluted-74394a8) | 4 |
| legacy-stay-in-place | 3 |
| closure-aggregator | 1 |
| audit/review | 1 |
| candidate-needs-bridge (qm cluster) | 1 |
| candidate-promoted-or-active (legacy uv) | 1 |

### 2.x. `audit/review` (1)

- `research/external_review_2026-04-25`

### 2.x. `candidate-active (3-phase op-*)` (25)

- `research/op-alpha-fine-structure`
- `research/op-bh-alpha-threshold`
- `research/op-cross-sector-charge`
- `research/op-eps-photon-ring`
- `research/op-eta-wolfenstein`
- `research/op-eta2-denom-derivation`
- `research/op-iota-charge-pmns-unification`
- `research/op-kappa-mixing-numerator`
- `research/op-mu-pmns-phase-hardening`
- `research/op-nu-majorana-phase-mbb`
- `research/op-omega1-substrate-em-coupling`
- `research/op-omicron1-sigmamnu-cosmo`
- `research/op-phi1-substrate-action-variational`
- `research/op-pi1-bb0nu-nme-isotope`
- `research/op-rho1-71Ge-cross-section`
- `research/op-sc-alpha-origin`
- `research/op-sigma1-substrate-light-dispersion`
- `research/op-tau1-closure-overlap-coulomb`
- `research/op-tau2-substrate-time-coupling`
- `research/op-theta-quark-koide`
- `research/op-upsilon1-closure-cross-family`
- `research/op-uv-as-ngfp`
- `research/op-xi-photon-ring`
- `research/op-xi2-sterile-nu-5sector`
- `research/op-zeta-mass-spectrum`

### 2.x. `candidate-active (heuristic)` (25)

- `research/atom_from_soliton`
- `research/brannen_sqrt2`
- `research/casimir_mof`
- `research/desi_dark_energy`
- `research/em_from_substrate`
- `research/liquid_viscosity`
- `research/mass_scaling_k4`
- `research/muon_g_minus_2`
- `research/nbody`
- `research/op-eht`
- `research/op-eht-A`
- `research/op-g0-r3-from-canonical-projection`
- `research/op-lambda1-e2-amplitude-emergence`
- `research/op-m92`
- `research/op-mu1-minimal-substrate-log-redefinition`
- `research/op-phase1-covariant`
- `research/op-phase2-quantum-gravity`
- `research/op-phase3-uv-completion`
- `research/op-quantum-closure`
- `research/op-uv3-phi0-renormalization`
- `research/particle_sector_closure`
- `research/rho_normal_state_closure`
- `research/superconductivity_closure`
- `research/thermal_transport_molecular`
- `research/why_n3`

### 2.x. `candidate-archive (heuristic, verify)` (6)

- `research/cosmo_tensions`
- `research/galaxy_scaling`
- `research/op-cosmology-closure`
- `research/op-newton-momentum`
- `research/op-psi1-substrate-light-acceleration`
- `research/op-tau3-substrate-clock-acceleration`

### 2.x. `candidate-needs-bridge (qm cluster)` (1)

- `research/qm_measurement`

### 2.x. `candidate-needs-bridge-or-sandbox (heuristic)` (19)

- `research/atomic_shells_closure`
- `research/cabibbo_correction`
- `research/cohesion_closure`
- `research/continuum_limit`
- `research/hubble_tension`
- `research/metric_ansatz`
- `research/neutrino_msw`
- `research/op-delta1-g-tilde-derivation`
- `research/op-delta2-Nf-derivation`
- `research/op-gamma1-phi-eff-anchor-resolution`
- `research/op-uv-renormalizability-research`
- `research/qm_born_rule`
- `research/qm_decoherence`
- `research/qm_entanglement`
- `research/qm_foundations`
- `research/qm_spin`
- `research/qm_statistics`
- `research/qm_superposition`
- `research/s8_tension`

### 2.x. `candidate-promoted-or-active (legacy uv)` (1)

- `research/uv_completion`

### 2.x. `closure-aggregator` (1)

- `research/closure_2026-04-26`

### 2.x. `legacy-stay-in-place` (3)

- `research/op1-op2-op4`
- `research/op6`
- `research/op7`

### 2.x. `needs-bridge (polluted-74394a8)` (4)

- `research/op-chi1-newton-constant-derivation`
- `research/op-omega2-axion-coupling-lock`
- `research/op-omega3-axion-decay-constant`
- `research/op-uv2-mtgp-absolute-scale`

## 3. Top 15 najgęściej zalinkowanych folderów

> Liczba referencji do ścieżki `research/<folder>` w `INDEX.md` + `DEPENDENCIES.md` +
> `DEPENDENCIES_REVERSE.md`. Wysoka liczba = drogie do fizycznego przeniesienia.
> W wariancie minimalnym **niczego nie ruszamy**, ale to lista folderów, których w
> Sesji 8.5 NIE ARCHIWIZUJEMY (zbyt drogie).

| Folder | INDEX.md | DEPENDENCIES.md | DEP_REV.md | Total |
|---|---:|---:|---:|---:|
| `nbody` | 0 | 15 | 13 | **28** |
| `op-psi1-substrate-light-acceleration` | 18 | 0 | 0 | **18** |
| `em_from_substrate` | 0 | 7 | 9 | **16** |
| `atom_from_soliton` | 0 | 3 | 7 | **10** |
| `op-phase3-uv-completion` | 10 | 0 | 0 | **10** |
| `atomic_shells_closure` | 0 | 3 | 6 | **9** |
| `op-alpha-fine-structure` | 9 | 0 | 0 | **9** |
| `op-cross-sector-charge` | 9 | 0 | 0 | **9** |
| `op-eps-photon-ring` | 9 | 0 | 0 | **9** |
| `op-eta-wolfenstein` | 9 | 0 | 0 | **9** |
| `op-eta2-denom-derivation` | 9 | 0 | 0 | **9** |
| `op-iota-charge-pmns-unification` | 9 | 0 | 0 | **9** |
| `op-kappa-mixing-numerator` | 9 | 0 | 0 | **9** |
| `op-omega1-substrate-em-coupling` | 9 | 0 | 0 | **9** |
| `op-phi1-substrate-action-variational` | 9 | 0 | 0 | **9** |

**Foldery z 0 linkami w INDEX/DEPENDENCIES:** 17 sztuk
(bezpieczne kandydaty do Sesji 8.5 archiwizacji, jeśli spełnią kryteria 3.4).

<details><summary>Lista folderów z 0 linkami</summary>

- `research/external_review_2026-04-25` — audit/review
- `research/op-cosmology-closure` — candidate-archive (heuristic, verify)
- `research/op-delta1-g-tilde-derivation` — candidate-needs-bridge-or-sandbox (heuristic)
- `research/op-delta2-Nf-derivation` — candidate-needs-bridge-or-sandbox (heuristic)
- `research/op-eht` — candidate-active (heuristic)
- `research/op-eht-A` — candidate-active (heuristic)
- `research/op-g0-r3-from-canonical-projection` — candidate-active (heuristic)
- `research/op-gamma1-phi-eff-anchor-resolution` — candidate-needs-bridge-or-sandbox (heuristic)
- `research/op-lambda1-e2-amplitude-emergence` — candidate-active (heuristic)
- `research/op-m92` — candidate-active (heuristic)
- `research/op-mu1-minimal-substrate-log-redefinition` — candidate-active (heuristic)
- `research/op-newton-momentum` — candidate-archive (heuristic, verify)
- `research/op-omega3-axion-decay-constant` — needs-bridge (polluted-74394a8)
- `research/op-uv3-phi0-renormalization` — candidate-active (heuristic)
- `research/op1-op2-op4` — legacy-stay-in-place
- `research/op6` — legacy-stay-in-place
- `research/op7` — legacy-stay-in-place

</details>

## 4. Aktywność (mtime ostatniego pliku w folderze)

### 4.1 Top 10 najnowszych (najprawdopodobniej `active`)

| Folder | mtime | klasa | PASS | CLOSED |
|---|---|---|---:|---:|
| `op-uv2-mtgp-absolute-scale` | 2026-05-02 | needs-bridge (polluted-74394a8) | 103 | 2 |
| `op-uv3-phi0-renormalization` | 2026-05-02 | candidate-active (heuristic) | 108 | 1 |
| `op-delta1-g-tilde-derivation` | 2026-05-02 | candidate-needs-bridge-or-sandbox (heuristic) | 2 | 1 |
| `closure_2026-04-26` | 2026-05-02 | closure-aggregator | 251 | 113 |
| `op-lambda1-e2-amplitude-emergence` | 2026-05-02 | candidate-active (heuristic) | 155 | 7 |
| `op-delta2-Nf-derivation` | 2026-05-02 | candidate-needs-bridge-or-sandbox (heuristic) | 4 | 0 |
| `op-g0-r3-from-canonical-projection` | 2026-05-02 | candidate-active (heuristic) | 206 | 51 |
| `op-gamma1-phi-eff-anchor-resolution` | 2026-05-02 | candidate-needs-bridge-or-sandbox (heuristic) | 9 | 0 |
| `op-mu1-minimal-substrate-log-redefinition` | 2026-05-02 | candidate-active (heuristic) | 42 | 1 |
| `op-chi1-newton-constant-derivation` | 2026-05-02 | needs-bridge (polluted-74394a8) | 104 | 9 |

### 4.2 Top 10 najstarszych (kandydaci do `archive` / `core-promoted` / legacy)

| Folder | mtime | klasa | PASS | FALSIFIED |
|---|---|---|---:|---:|
| `uv_completion` | 2026-04-14 | candidate-promoted-or-active (legacy uv) | 0 | 0 |
| `cabibbo_correction` | 2026-04-14 | candidate-needs-bridge-or-sandbox (heuristic) | 4 | 0 |
| `metric_ansatz` | 2026-04-14 | candidate-needs-bridge-or-sandbox (heuristic) | 11 | 0 |
| `qm_born_rule` | 2026-04-15 | candidate-needs-bridge-or-sandbox (heuristic) | 0 | 0 |
| `qm_measurement` | 2026-04-15 | candidate-needs-bridge (qm cluster) | 36 | 0 |
| `qm_superposition` | 2026-04-15 | candidate-needs-bridge-or-sandbox (heuristic) | 9 | 0 |
| `qm_statistics` | 2026-04-15 | candidate-needs-bridge-or-sandbox (heuristic) | 9 | 0 |
| `qm_decoherence` | 2026-04-15 | candidate-needs-bridge-or-sandbox (heuristic) | 9 | 0 |
| `qm_entanglement` | 2026-04-16 | candidate-needs-bridge-or-sandbox (heuristic) | 21 | 0 |
| `qm_foundations` | 2026-04-16 | candidate-needs-bridge-or-sandbox (heuristic) | 22 | 0 |

## 5. Top 15 folderów z najwyższym sygnałem PASS+CLOSED+DERIVED

> Heurystyka kandydatów do `core-ready` lub już `core-promoted`.
> **Anty-overclaim**: te liczby same nie wystarczą do nadania statusu — Sesja 4
> wymaga `source_of_status` cytującego konkretny plik. Polluted folders z 74394a8
> są wyłączone z dalszej promocji do czasu forward-patch (zob. Sesja 7).

| Folder | PASS | CLOSED | DERIVED | LOCKED | mtime | links | flagi |
|---|---:|---:|---:|---:|---|---:|---|
| `nbody` | 682 | 0 | 34 | 0 | 2026-04-22 | 28 | — |
| `op-cosmology-closure` | 592 | 42 | 0 | 0 | 2026-04-26 | 0 | — |
| `op-quantum-closure` | 486 | 139 | 0 | 0 | 2026-04-28 | 1 | — |
| `op-phase2-quantum-gravity` | 188 | 175 | 208 | 0 | 2026-04-28 | 2 | — |
| `op-phase3-uv-completion` | 206 | 154 | 122 | 0 | 2026-04-28 | 10 | — |
| `op-newton-momentum` | 358 | 49 | 1 | 4 | 2026-05-01 | 0 | — |
| `closure_2026-04-26` | 251 | 113 | 34 | 0 | 2026-05-02 | 4 | — |
| `op7` | 361 | 13 | 0 | 0 | 2026-04-25 | 0 | legacy-stay |
| `op-uv-as-ngfp` | 164 | 15 | 133 | 77 | 2026-04-29 | 9 | — |
| `op-phase1-covariant` | 187 | 90 | 1 | 0 | 2026-04-27 | 2 | — |
| `op-g0-r3-from-canonical-projection` | 206 | 51 | 1 | 0 | 2026-05-02 | 0 | — |
| `op-theta-quark-koide` | 146 | 14 | 94 | 51 | 2026-04-29 | 9 | — |
| `op-kappa-mixing-numerator` | 112 | 5 | 134 | 9 | 2026-04-30 | 9 | — |
| `op-xi-photon-ring` | 135 | 19 | 92 | 86 | 2026-04-29 | 9 | — |
| `op-eta2-denom-derivation` | 138 | 6 | 101 | 23 | 2026-04-30 | 9 | — |

## 6. Sygnał archiwalny (FALSIFIED/WITHDRAWN/OBSOLETE/SUPERSEDED)

> Top 10 folderów z najwyższym sygnałem archiwalnym.
> **NIE** kwalifikuje automatycznie do `_archive/` — Sesja 8.5 wymaga 2 źródeł
> per folder + akceptację człowieka per folder.

| Folder | FALSIFIED | WITHDRAWN | OBSOLETE | SUPERSEDED | links | mtime |
|---|---:|---:|---:|---:|---:|---|
| `op-cosmology-closure` | 30 | 0 | 1 | 79 | 0 | 2026-04-26 |
| `op-omicron1-sigmamnu-cosmo` | 36 | 0 | 0 | 0 | 6 | 2026-04-30 |
| `why_n3` | 35 | 0 | 0 | 0 | 2 | 2026-05-01 |
| `op-g0-r3-from-canonical-projection` | 33 | 0 | 1 | 0 | 0 | 2026-05-02 |
| `op-eps-photon-ring` | 29 | 0 | 0 | 0 | 9 | 2026-04-29 |
| `op-pi1-bb0nu-nme-isotope` | 25 | 0 | 0 | 0 | 6 | 2026-04-30 |
| `op-nu-majorana-phase-mbb` | 24 | 0 | 0 | 0 | 6 | 2026-04-30 |
| `op-theta-quark-koide` | 24 | 0 | 0 | 0 | 9 | 2026-04-29 |
| `op-iota-charge-pmns-unification` | 23 | 0 | 0 | 0 | 9 | 2026-04-30 |
| `op-kappa-mixing-numerator` | 22 | 0 | 0 | 0 | 9 | 2026-04-30 |

## 7. Foldery zatrute incydentem 74394a8

> Z [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]. **Plan**: w Sesji 4 te foldery
> dostają `level: unknown`, `core_compatibility: unknown`,
> `folder_status: needs-bridge` z notatką. **Nie wolno** awansować ich w Sesji 7
> bez forward-patch.

| Folder | PASS | FULL_CONVERGENCE | mtime | links |
|---|---:|---:|---|---:|
| `op-chi1-newton-constant-derivation` | 104 | 15 | 2026-05-02 | 5 |
| `op-omega2-axion-coupling-lock` | 108 | 1 | 2026-05-01 | 3 |
| `op-omega3-axion-decay-constant` | 109 | 6 | 2026-05-01 | 0 |
| `op-uv2-mtgp-absolute-scale` | 103 | 14 | 2026-05-02 | 3 |

## 8. Closure-aggregator i foldery z podfolderami

| Folder | subfolders | nazwy podfolderów |
|---|---:|---|
| `brannen_sqrt2` | 2 | `TGP`, `research` |
| `closure_2026-04-26` | 4 | `Lambda_from_Phi0`, `alpha_psi_threshold`, `f_psi_principle`, `sigma_ab_pathB` |
| `em_from_substrate` | 1 | `TGP` |
| `galaxy_scaling` | 1 | `Rotmod_LTG` |
| `nbody` | 6 | `__pycache__`, `_archiwum_docs`, `docs`, `examples`, `paper`, `plots` |
| `op-g0-r3-from-canonical-projection` | 1 | `.claude` |
| `op-newton-momentum` | 1 | `__pycache__` |
| `op-quantum-closure` | 1 | `__pycache__` |
| `op-uv2-mtgp-absolute-scale` | 1 | `.claude` |
| `op1-op2-op4` | 1 | `__pycache__` |
| `op6` | 1 | `__pycache__` |
| `rho_normal_state_closure` | 1 | `__pycache__` |
| `superconductivity_closure` | 2 | `TGP`, `__pycache__` |

## 9. Diff vault-rootowego `./research/` vs canonical

Vault-rootowy `./research/` nie istnieje.
## 10. Top-level pliki w `research/` (programowe dokumenty, NIE foldery)

Te pliki **zostają in place** — to programowe dokumenty (status / redirect / new-directions),
nie foldery badawcze. Plan ich nie dotyka.

- `graph_concept_flow.gexf`
- `graph_concept_flow.png`
- `graph_concept_flow_v2.gexf`
- `graph_concept_flow_v2.png`
- `graph_python_modules.gexf`
- `graph_python_modules.png`
- `graph_research_problems.gexf`
- `graph_research_problems.png`
- `graph_tex_includes.gexf`
- `graph_tex_includes.png`
- `tgp_dependency_graph.py`

## 11. Pełna lista folderów (kompaktowa tabela)

| # | Folder | klasa heurystyczna | program | phases | YAML | PASS | CLOSED | FALS | mtime | links |
|---:|---|---|:-:|---:|:-:|---:|---:|---:|---|---:|
| 1 | `atom_from_soliton` | candidate-active (heuristic) | — | 0 | — | 38 | 0 | 0 | 2026-04-22 | 10 |
| 2 | `atomic_shells_closure` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 5 | 1 | 0 | 2026-04-22 | 9 |
| 3 | `brannen_sqrt2` | candidate-active (heuristic) | — | 0 | — | 70 | 9 | 0 | 2026-04-17 | 3 |
| 4 | `cabibbo_correction` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 4 | 0 | 0 | 2026-04-14 | 2 |
| 5 | `casimir_mof` | candidate-active (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-22 | 7 |
| 6 | `closure_2026-04-26` | closure-aggregator | — | 0 | — | 251 | 113 | 4 | 2026-05-02 | 4 |
| 7 | `cohesion_closure` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 20 | 0 | 0 | 2026-04-22 | 7 |
| 8 | `continuum_limit` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 5 | 2 | 0 | 2026-04-18 | 2 |
| 9 | `cosmo_tensions` | candidate-archive (heuristic, verify) | — | 0 | — | 32 | 0 | 0 | 2026-04-26 | 3 |
| 10 | `desi_dark_energy` | candidate-active (heuristic) | — | 0 | — | 8 | 0 | 21 | 2026-04-22 | 4 |
| 11 | `em_from_substrate` | candidate-active (heuristic) | — | 0 | — | 42 | 0 | 0 | 2026-04-22 | 16 |
| 12 | `external_review_2026-04-25` | audit/review | — | 0 | — | 23 | 1 | 0 | 2026-04-26 | 0 |
| 13 | `galaxy_scaling` | candidate-archive (heuristic, verify) | — | 0 | — | 72 | 4 | 13 | 2026-04-26 | 6 |
| 14 | `hubble_tension` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-18 | 2 |
| 15 | `liquid_viscosity` | candidate-active (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-22 | 6 |
| 16 | `mass_scaling_k4` | candidate-active (heuristic) | — | 0 | — | 32 | 8 | 0 | 2026-05-02 | 2 |
| 17 | `metric_ansatz` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 11 | 0 | 0 | 2026-04-14 | 2 |
| 18 | `muon_g_minus_2` | candidate-active (heuristic) | — | 0 | — | 0 | 0 | 1 | 2026-04-22 | 7 |
| 19 | `nbody` | candidate-active (heuristic) | — | 0 | — | 682 | 0 | 11 | 2026-04-22 | 28 |
| 20 | `neutrino_msw` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-22 | 3 |
| 21 | `op-alpha-fine-structure` | candidate-active (3-phase op-*) | ✓ | 6 | — | 127 | 6 | 16 | 2026-04-30 | 9 |
| 22 | `op-bh-alpha-threshold` | candidate-active (3-phase op-*) | ✓ | 6 | — | 144 | 11 | 0 | 2026-04-28 | 7 |
| 23 | `op-chi1-newton-constant-derivation` | needs-bridge (polluted-74394a8) | ✓ | 6 | — | 104 | 9 | 1 | 2026-05-02 | 5 |
| 24 | `op-cosmology-closure` | candidate-archive (heuristic, verify) | — | 0 | — | 592 | 42 | 30 | 2026-04-26 | 0 |
| 25 | `op-cross-sector-charge` | candidate-active (3-phase op-*) | ✓ | 6 | — | 138 | 10 | 0 | 2026-04-28 | 9 |
| 26 | `op-delta1-g-tilde-derivation` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | ✓ | 2 | 1 | 0 | 2026-05-02 | 0 |
| 27 | `op-delta2-Nf-derivation` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | ✓ | 4 | 0 | 0 | 2026-05-02 | 0 |
| 28 | `op-eht` | candidate-active (heuristic) | — | 0 | — | 112 | 3 | 0 | 2026-04-25 | 0 |
| 29 | `op-eht-A` | candidate-active (heuristic) | — | 0 | — | 50 | 4 | 0 | 2026-04-25 | 0 |
| 30 | `op-eps-photon-ring` | candidate-active (3-phase op-*) | ✓ | 6 | — | 166 | 16 | 29 | 2026-04-29 | 9 |
| 31 | `op-eta-wolfenstein` | candidate-active (3-phase op-*) | ✓ | 6 | — | 121 | 7 | 20 | 2026-04-29 | 9 |
| 32 | `op-eta2-denom-derivation` | candidate-active (3-phase op-*) | ✓ | 6 | — | 138 | 6 | 8 | 2026-04-30 | 9 |
| 33 | `op-g0-r3-from-canonical-projection` | candidate-active (heuristic) | — | 6 | ✓ | 206 | 51 | 33 | 2026-05-02 | 0 |
| 34 | `op-gamma1-phi-eff-anchor-resolution` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | ✓ | 9 | 0 | 0 | 2026-05-02 | 0 |
| 35 | `op-iota-charge-pmns-unification` | candidate-active (3-phase op-*) | ✓ | 6 | — | 106 | 3 | 23 | 2026-04-30 | 9 |
| 36 | `op-kappa-mixing-numerator` | candidate-active (3-phase op-*) | ✓ | 6 | — | 112 | 5 | 22 | 2026-04-30 | 9 |
| 37 | `op-lambda1-e2-amplitude-emergence` | candidate-active (heuristic) | — | 8 | ✓ | 155 | 7 | 9 | 2026-05-02 | 0 |
| 38 | `op-m92` | candidate-active (heuristic) | — | 0 | — | 31 | 0 | 0 | 2026-04-25 | 0 |
| 39 | `op-mu-pmns-phase-hardening` | candidate-active (3-phase op-*) | ✓ | 6 | — | 99 | 3 | 7 | 2026-04-30 | 6 |
| 40 | `op-mu1-minimal-substrate-log-redefinition` | candidate-active (heuristic) | — | 0 | ✓ | 42 | 1 | 0 | 2026-05-02 | 0 |
| 41 | `op-newton-momentum` | candidate-archive (heuristic, verify) | — | 0 | — | 358 | 49 | 11 | 2026-05-01 | 0 |
| 42 | `op-nu-majorana-phase-mbb` | candidate-active (3-phase op-*) | ✓ | 6 | — | 102 | 0 | 24 | 2026-04-30 | 6 |
| 43 | `op-omega1-substrate-em-coupling` | candidate-active (3-phase op-*) | ✓ | 6 | — | 63 | 3 | 8 | 2026-05-01 | 9 |
| 44 | `op-omega2-axion-coupling-lock` | needs-bridge (polluted-74394a8) | ✓ | 6 | — | 108 | 0 | 2 | 2026-05-01 | 3 |
| 45 | `op-omega3-axion-decay-constant` | needs-bridge (polluted-74394a8) | ✓ | 6 | — | 109 | 1 | 1 | 2026-05-01 | 0 |
| 46 | `op-omicron1-sigmamnu-cosmo` | candidate-active (3-phase op-*) | ✓ | 6 | — | 129 | 0 | 36 | 2026-04-30 | 6 |
| 47 | `op-phase1-covariant` | candidate-active (heuristic) | — | 8 | — | 187 | 90 | 0 | 2026-04-27 | 2 |
| 48 | `op-phase2-quantum-gravity` | candidate-active (heuristic) | — | 8 | — | 188 | 175 | 0 | 2026-04-28 | 2 |
| 49 | `op-phase3-uv-completion` | candidate-active (heuristic) | — | 9 | — | 206 | 154 | 0 | 2026-04-28 | 10 |
| 50 | `op-phi1-substrate-action-variational` | candidate-active (3-phase op-*) | ✓ | 6 | — | 60 | 0 | 7 | 2026-04-30 | 9 |
| 51 | `op-pi1-bb0nu-nme-isotope` | candidate-active (3-phase op-*) | ✓ | 6 | — | 111 | 0 | 25 | 2026-04-30 | 6 |
| 52 | `op-psi1-substrate-light-acceleration` | candidate-archive (heuristic, verify) | ✓ | 13 | — | 169 | 23 | 0 | 2026-05-01 | 18 |
| 53 | `op-quantum-closure` | candidate-active (heuristic) | — | 0 | — | 486 | 139 | 0 | 2026-04-28 | 1 |
| 54 | `op-rho1-71Ge-cross-section` | candidate-active (3-phase op-*) | ✓ | 6 | — | 107 | 0 | 0 | 2026-04-30 | 9 |
| 55 | `op-sc-alpha-origin` | candidate-active (3-phase op-*) | ✓ | 6 | — | 82 | 10 | 0 | 2026-04-28 | 7 |
| 56 | `op-sigma1-substrate-light-dispersion` | candidate-active (3-phase op-*) | ✓ | 6 | — | 60 | 1 | 6 | 2026-05-01 | 9 |
| 57 | `op-tau1-closure-overlap-coulomb` | candidate-active (3-phase op-*) | ✓ | 6 | — | 105 | 0 | 6 | 2026-04-30 | 9 |
| 58 | `op-tau2-substrate-time-coupling` | candidate-active (3-phase op-*) | ✓ | 6 | — | 66 | 0 | 19 | 2026-04-30 | 9 |
| 59 | `op-tau3-substrate-clock-acceleration` | candidate-archive (heuristic, verify) | ✓ | 7 | — | 88 | 20 | 2 | 2026-05-01 | 9 |
| 60 | `op-theta-quark-koide` | candidate-active (3-phase op-*) | ✓ | 6 | — | 146 | 14 | 24 | 2026-04-29 | 9 |
| 61 | `op-upsilon1-closure-cross-family` | candidate-active (3-phase op-*) | ✓ | 6 | — | 61 | 0 | 2 | 2026-04-30 | 9 |
| 62 | `op-uv-as-ngfp` | candidate-active (3-phase op-*) | ✓ | 6 | — | 164 | 15 | 0 | 2026-04-29 | 9 |
| 63 | `op-uv-renormalizability-research` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | ✓ | 3 | 7 | 0 | 2026-04-28 | 1 |
| 64 | `op-uv2-mtgp-absolute-scale` | needs-bridge (polluted-74394a8) | ✓ | 6 | — | 103 | 2 | 0 | 2026-05-02 | 3 |
| 65 | `op-uv3-phi0-renormalization` | candidate-active (heuristic) | ✓ | 3 | ✓ | 108 | 1 | 0 | 2026-05-02 | 0 |
| 66 | `op-xi-photon-ring` | candidate-active (3-phase op-*) | ✓ | 6 | — | 135 | 19 | 4 | 2026-04-29 | 9 |
| 67 | `op-xi2-sterile-nu-5sector` | candidate-active (3-phase op-*) | ✓ | 6 | — | 94 | 0 | 14 | 2026-04-30 | 6 |
| 68 | `op-zeta-mass-spectrum` | candidate-active (3-phase op-*) | ✓ | 6 | — | 140 | 11 | 18 | 2026-04-29 | 9 |
| 69 | `op1-op2-op4` | legacy-stay-in-place | — | 0 | — | 0 | 4 | 0 | 2026-04-27 | 0 |
| 70 | `op6` | legacy-stay-in-place | — | 0 | — | 1 | 3 | 0 | 2026-04-24 | 0 |
| 71 | `op7` | legacy-stay-in-place | — | 0 | — | 361 | 13 | 1 | 2026-04-25 | 0 |
| 72 | `particle_sector_closure` | candidate-active (heuristic) | — | 0 | — | 0 | 5 | 0 | 2026-05-01 | 5 |
| 73 | `qm_born_rule` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-15 | 2 |
| 74 | `qm_decoherence` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 9 | 0 | 0 | 2026-04-15 | 2 |
| 75 | `qm_entanglement` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 21 | 1 | 0 | 2026-04-16 | 2 |
| 76 | `qm_foundations` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 22 | 0 | 0 | 2026-04-16 | 2 |
| 77 | `qm_measurement` | candidate-needs-bridge (qm cluster) | — | 0 | — | 36 | 1 | 0 | 2026-04-15 | 2 |
| 78 | `qm_spin` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 10 | 0 | 0 | 2026-05-01 | 2 |
| 79 | `qm_statistics` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 9 | 0 | 0 | 2026-04-15 | 3 |
| 80 | `qm_superposition` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 9 | 1 | 0 | 2026-04-15 | 2 |
| 81 | `rho_normal_state_closure` | candidate-active (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-21 | 3 |
| 82 | `s8_tension` | candidate-needs-bridge-or-sandbox (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-18 | 2 |
| 83 | `superconductivity_closure` | candidate-active (heuristic) | — | 0 | — | 137 | 1 | 0 | 2026-04-22 | 8 |
| 84 | `thermal_transport_molecular` | candidate-active (heuristic) | — | 0 | — | 0 | 0 | 0 | 2026-04-22 | 5 |
| 85 | `uv_completion` | candidate-promoted-or-active (legacy uv) | — | 0 | — | 0 | 0 | 0 | 2026-04-14 | 2 |
| 86 | `why_n3` | candidate-active (heuristic) | — | 6 | — | 146 | 25 | 35 | 2026-05-01 | 2 |

## 12. Pytania pozostające otwarte po Sesji 1

Decyzje 1–6 z [[meta/PLAN_RESEARCH_WORKFLOW_v1.md#9-decyzje-człowieka-zatwierdzone-2026-05-02]]
już są zatwierdzone. Pytania, które ten audyt **dodatkowo** odsłonił:

**Q1.** Vault-root `./research/` (11 folderów) ma 100% pokrycie z canonical (zob. sekcja 9).
Czy człowiek chce:
- (a) zachować bez zmian (out-of-scope na zawsze),
- (b) zaplanować kasowanie po zakończeniu workflow (Sesja 9.x dodatkowo),
- (c) zsynchronizować z canonical jako mirror (np. via skrypt po commitach)?

**Q2.** Foldery z subfolderami (sekcja 8) — czy każde podzamknięcie
(np. `closure_2026-04-26/sigma_ab_pathB/`) traktujemy jako osobny **byt**
z własnym YAMLem, NEEDS, FINDINGS, czy tylko rodzic dostaje pełną warstwę,
a dzieci są lekkie?

**Q3.** Top-level pliki w `research/` (sekcja 10): czy program-doc dostaje też
YAML `tgp_status` (`folder_status: program-doc`, `kind: program-doc`),
czy zostaje bez frontmattera?

**Q4.** *(nieaktualne)* — żaden folder nie ma jeszcze `FINDINGS.md`/`NEEDS.md`, więc
Sesja 5 zaczyna od zera. Brak konfliktu z istniejącymi plikami.

## 13. Następne kroki

Po akceptacji odpowiedzi Q1–Q4:

- **Sesja 2**: szablony README/NEEDS/FINDINGS + 3 case-study previews.
- **Sesja 3**: warstwa `meta/research/` + `meta/core/` + fizyczne `_sandbox`/`_archive`.
- **Sesja 4**: klasyfikacja YAML w 86 folderach (batch po 10 → ~9 batchy).

---

*Raport wygenerowany automatycznie. Surowe dane: [[_audit_s1_raw.json]].*