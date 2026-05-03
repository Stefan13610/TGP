---
title: "CORE_CANDIDATES — foldery gotowe do awansu (lub już awansowane) do core"
date: 2026-05-03
type: candidates
status: STUB (Sesja 3) — wypełniany w Sesji 7
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/core/CORE_INTAKE.md]]"
  - "[[meta/core/CORE_INVENTORY.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
tags:
  - core-ready
  - core-promoted
  - human-gate
---

# CORE_CANDIDATES — kandydaci do awansu do core

> **Cel pliku:** lista folderów z `tgp_status.folder_status: core-ready`
> lub `core-promoted` — z linkiem do INTAKE i decyzji człowieka.
>
> **Twardy gate:** żaden folder nie awansuje do core bez wpisu tutaj
> + akceptacji człowieka w `INTAKE_*.md`.

---

## 0. Reguły wpisu

1. Wpis powstaje **automatycznie** w Sesji 7, gdy folder otrzymuje
   `tgp_status.folder_status: core-ready` lub `core-promoted`.
2. **`core-ready`** = `level: L3`, brak otwartych blokujących bridge'y,
   `core_compatibility ∈ {current, partial}`, NIE w kwarantannie
   `polluted_74394a8`.
3. **`core-promoted`** = wynik **już** jest w core; folder żyje jako
   historia. `promoted_to_core` ma weryfikowalną referencję.
4. **Required human approval gate**: `core-ready → core-promoted`
   wymaga zaakceptowanego `INTAKE_*.md` przez człowieka.

---

## 1. Statystyka

| Status | Liczba folderów |
|--------|----------------:|
| `core-ready` (oczekuje INTAKE) | **5** (Sesja 7 stay-as-L3 — niższy priorytet) |
| `core-ready` (INTAKE w `PENDING-HUMAN-REVIEW`) | **9** |
| `core-promoted` (już w core) | 0 |
| Łącznie | **14** kandydatów |

(Sesja 3.5 Q6: 4 INTAKE dla closure_2026-04-26 children. Sesja 7: 5 dodatkowych
INTAKE dla top L3 op-* + 5 stay-as-L3 z niższym priorytetem.)

---

## 2. core-ready (oczekują INTAKE)

| Folder | Level | Zamykany hotspot | Proponowany target | Source of readiness | Action |
|--------|-------|------------------|--------------------|---------------------|--------|
| (puste w Sesji 3) |  |  |  |  |  |

---

## 3. core-ready (INTAKE w PENDING-HUMAN-REVIEW)

### 3.1 Sesja 3.5 Q6 — closure_2026-04-26 children (4)

| Folder | INTAKE file | Hotspot | Data wniosku | Decyzja człowieka |
|--------|-------------|---------|---------------|--------------------|
| `closure_2026-04-26/sigma_ab_pathB` | [[meta/core/intake/INTAKE_closure_2026-04-26_sigma_ab_pathB_2026-05-03.md]] | H-B6 (alt track) | 2026-05-03 | PENDING |
| `closure_2026-04-26/f_psi_principle` | [[meta/core/intake/INTAKE_closure_2026-04-26_f_psi_principle_2026-05-03.md]] | H-B8 (alt track) | 2026-05-03 | PENDING |
| `closure_2026-04-26/alpha_psi_threshold` | [[meta/core/intake/INTAKE_closure_2026-04-26_alpha_psi_threshold_2026-05-03.md]] | H-B9 (alt track) | 2026-05-03 | PENDING |
| `closure_2026-04-26/Lambda_from_Phi0` | [[meta/core/intake/INTAKE_closure_2026-04-26_Lambda_from_Phi0_2026-05-03.md]] | NEW (vacuum catastrophe, poza 43-list) | 2026-05-03 | PENDING |

### 3.2 Sesja 7 — top 5 L3 op-* shortlist (5)

| Folder | INTAKE file | Hotspot / Cross-link | Sektor | Decyzja |
|--------|-------------|----------------------|--------|---------|
| `op-uv-as-ngfp` | [[meta/core/intake/INTAKE_op-uv-as-ngfp_2026-05-03.md]] | NEW (UV completion, status cascade ξ.1+XS.1+UV7→DERIVED) | UV completion | PENDING |
| `op-theta-quark-koide` | [[meta/core/intake/INTAKE_op-theta-quark-koide_2026-05-03.md]] | NEW (extension Koide K=2/3 lepton → quark) | quark sector | PENDING |
| `op-xi-photon-ring` | [[meta/core/intake/INTAKE_op-xi-photon-ring_2026-05-03.md]] | NEW (photon-ring decomposition; cross-validation z UV.1) | BH shadow | PENDING |
| `op-zeta-mass-spectrum` | [[meta/core/intake/INTAKE_op-zeta-mass-spectrum_2026-05-03.md]] | H-B4 (Σm_ν upgrade z annotation → derivation) | neutrino mass | PENDING |
| `op-eta-wolfenstein` | [[meta/core/intake/INTAKE_op-eta-wolfenstein_2026-05-03.md]] | NEW (CKM Wolfenstein A, ρ̄, η̄ first-principles) | CKM/PMNS | PENDING |

### 3.3 Sesja 7 — pozostałe 5 L3 op-* (stay-as-L3, niższy priorytet)

Te 5 folderów ma `level: L3` ale nie wpadły w shortlist Sesji 7. Pozostają
jako `core-ready` bez INTAKE — kandydaci do INTAKE w przyszłej sesji
(zwłaszcza po review pierwszych 9 wniosków przez człowieka).

| Folder | Cycle | Powód odłożenia |
|--------|-------|------------------|
| `op-alpha-fine-structure` | α.1 | Topical overlap z `dodatekT_alpha` (już w core); może być duplicate path |
| `op-cross-sector-charge` | XS.1 | Cross-validated z UV.1 cascade; sensible jako bundle z UV.1 (osobny INTAKE bundle) |
| `op-eps-photon-ring` | ε.1 | Topical overlap z ξ.1 (op-xi-photon-ring); 2-fold consistency check, możliwy bundle |
| `op-eta2-denom-derivation` | η.2 | Cluster z η.1 (op-eta-wolfenstein); cluster INTAKE proponowany w INTAKE_op-eta-wolfenstein §8 |
| `op-sc-alpha-origin` | SC.1 | Bundle candidate z ξ.1 (oba photon-ring) lub UV.1 (cascade) |

---

## 4. core-promoted (już w core)

| Folder | promoted_to_core | Data promocji | INTAKE wniosek | INVENTORY ID |
|--------|-------------------|---------------|------------------|---------------|
| (puste w Sesji 3 — Sesja 7 znajdzie pre-existing) |  |  |  |  |

---

## 5. Shortlist kandydatów (zaktualizowana po Sesji 3.5)

> **Reframing post-Sesja 3.5:** `closure_2026-04-26/*` foldery zamykają
> hotspoty H-B6, H-B8, H-B9 — ale audyt 2026-05-01 **niezależnie** zamknął
> te same hotspoty przez `op-newton-momentum/B6/B8/B9_*.py` skrypty
> (§ U/V/W). To są **alternative formal tracks** (cross-validation), nie
> primary closure paths.
>
> Pełen reality-check: [[meta/research/HOTSPOT_AUDIT_S3_5.md]] §6.3.

### 5.1 Cross-validation INTAKE candidates (closure_2026-04-26)

| # | Folder | Adresuje (cross-val) | Audyt main closure | Decyzja człowieka |
|---|--------|----------------------|---------------------|--------------------|
| 1 | `closure_2026-04-26/sigma_ab_pathB` (11/11 PASS) | tensor σ_ab Path B PRIMARY (H-B6 alt) | § U (B6 6/6 PASS via M9.x re-run) | czy core ma cytować obie ścieżki? |
| 2 | `closure_2026-04-26/f_psi_principle` (12/12 PASS) | n=deg(V)=4 unique exponent (H-B8 alt) | § V (B8 5/5 PASS via independence check) | czy core ma cytować obie? |
| 3 | `closure_2026-04-26/alpha_psi_threshold` (5/5 PASS) | α(ψ) threshold for OP-M92 (H-B9 alt) | § W (B9 6/6 PASS via composition test) | czy core ma cytować obie? |

### 5.2 NEW structural contribution candidates (poza 43-item list)

| # | Folder | Strukturalny wkład | Audyt status |
|---|--------|--------------------|--------------|
| 4 | `closure_2026-04-26/Lambda_from_Phi0` (7/7 PASS) | Ω_Λ input → prediction (vacuum catastrophe avoided structurally) | NIE w 43-item list audytu; niezależny strukturalny kandydat |
| 5 | `research/why_n3` (5-phase emergent Dirac propagator, 6/6 PASS m_obs=cA^(5-α)) | R3↔M9.1'' horizon coincidence; m_μ/m_e diff −0.001% PDG | § AB cross-cycle deepening; warstwa 3c szkic strukturalny zamknięty |

### 5.3 INTAKE workflow rekomendacja

**Decyzja proponowana w Sesji 7:**

- **§ 5.1 (cross-validation, 3 itemy):** hybryda INTAKE — jeden wniosek
  od rodzica `closure_2026-04-26` z 3 sekcjami (każda dziecko = osobna
  sekcja decyzyjna). Człowiek decyduje per case czy core cytuje:
  (a) tylko § U/V/W path (audyt main),
  (b) tylko closure_2026-04-26 path (alternative),
  (c) obie ścieżki jako cross-validation.
- **§ 5.2 (NEW, 2 itemy):** osobne INTAKE wnioski — to są nowe
  strukturalne wkłady, nie zamknięcia hotspotów. `Lambda_from_Phi0`
  zasługuje na osobny rozdział w `core/sek05_ciemna_energia`;
  `why_n3` jest długoterminowy (Phase 6+), na razie INTAKE
  preliminary.

**Brak innych core-ready kandydatów** w obecnym przeglądzie. Foldery
z wysokim PASS-count w audycie Sesji 1 (`op-uv-as-ngfp` PASS=164 LOCKED=77,
`op-phase3-uv-completion` PASS=206 CLOSED=154 DERIVED=122) wymagają
ręcznej oceny w Sesji 4 — ich wyniki mogą już być częścią core (UV completion
sekcji `sek08*`, gauge sektor `sek09`) lub kandydatami na Sesję 7 INTAKE.

---

## 6. Synchronizacja z innymi plikami

- Zmiana `tgp_status.folder_status` w `research/<X>/README.md` ⇒
  obowiązkowy update tabel § 2/3/4.
- Akceptacja `INTAKE_*.md` ⇒ § 3 → § 4 (zostaje wpis history).
- `core-promoted` w § 4 ⇒ wpis w `CORE_INVENTORY.md` § 1.X (ten sam ID).

---

## 7. Anty-overclaim w tym pliku

- `core-promoted` w § 4 wymaga `promoted_to_core` z weryfikowalną
  ścieżką do core (`main.tex §X.Y` lub `core/sek0Z_*` lub `axioms/*`).
- Skrypt-checker (Sesja 8) ⇒ jeśli plik core wskazany w `promoted_to_core`
  nie zawiera odniesienia do source folderu (np. komentarza
  `% promoted from research/<X>/Phase3_results.md`), wpis jest
  CRITICAL finding.
