---
title: "G.0 Phase 3 — Integration audit specifikacja (sek08 v2.0 plan, NIE actual mods)"
date: 2026-05-02
phase: 3
parent: "[[README.md]]"
predecessor: "[[Phase2_results.md]]"
status: ACTIVE — start 2026-05-02
score_gate: "≥3/4 PASS dla recommendation actual core mods"
scope_constraint: "Phase 3 produkuje TYLKO specifications/drafts. Actual core
                  modifications sek08a/sek08c są SEPARATE task po user approval."
---

# G.0 Phase 3 — Integration audit (specifications, NIE actual mods)

## 0. Scope constraint (KRYTYCZNE)

Per user instruction z 2026-05-02: *"przy aktualizacji tej sekcji niestety
trzeba zrobić pełny audyt wszystkich powiazanych z nią elementów a ona niestety
jest fundamentalna więc na poczatku osobny program "badawczy" z analizą dopiero
potem zcalenie tego z sekcją 8 jeżeli wszystko będzie się spinało"*.

**Phase 3 produkuje TYLKO:**
- Specification documents w `research/op-g0-r3-from-canonical-projection/`
- Cross-reference impact analysis
- Draft text dla proposed changes

**Phase 3 NIE robi:**
- Modyfikacji `core/sek08a/`
- Modyfikacji `core/sek08c/`
- Modyfikacji `TGP_FOUNDATIONS.md`

Actual core modifications = **separate later phase** (Phase 4 "Integration"),
po user review wynikow Phase 3.

---

## 1. Cel Phase 3

Po Phase 1 + Phase 2 CLOSED-POSITIVE z fundingiem V_M911, Phase 3 odpowiada na:

1. **Czy q·Φ_0 = 4πG_0/c² nadal trzyma w G.0 framework?** (P32)
2. **Co dokladnie wymaga zmiany w sek08a / sek08c?** (P31, P33, P34)
3. **Jak nowy κ wplynie na observational constraints?** (P32)
4. **Czy istnieja ukryte zaleznosci w core ktore zmiana wskaze?** (P33)

**Output Phase 3:** complete specification ready dla Phase 4 integration.

---

## 2. Sub-tasks (4)

### P32 — Newton limit re-derivation z V_M911 (PHYSICS CHECK)

**Cel:** Sprawdz, czy R3 ODE w static spherical weak-field limicie z matter
source (na akcji G.0) reprodukuje Newton's law z odpowiednią wartością G_0.

**Metoda:**
1. Linearize R3 ODE z matter source: ψ = 1 + δψ, |δψ| << 1, ρ_matter > 0
2. Identify Newtonian potential U z PPN: U = 2δψ (z P23)
3. Reduce do Poisson-like: ∇²δψ = m_sp²·δψ + (5q/Φ_0)·ρ_matter (z P24 source coupling)
4. W r << 1/m_sp (Solar System): ∇²δψ ≈ (5q/Φ_0)·ρ_matter
5. Demand ∇²U = 4πG_0·ρ_matter (Newton)
6. Solve: q·Φ_0 nowe vs sek08a stary
7. Numerical: integration R3 ODE z matter source (rho ~ r^-3 w Solar System), extract Newton constant
8. Compare obu: stary (q·Φ_0 = 4πG_0/c²) vs nowy (q·Φ_0 = ?)

**PASS criteria:**
- Sympy LOCK new q·Φ_0 relation
- Numerical Newton G_0 reproduction <1% error z R3 ODE solver
- Φ_0 nowe consistent z BBN/LLR/CMB constraints po update

**Plik:** `phase3_P32_newton_limit_rederivation.py`

---

### P33 — Cross-reference audit (impact analysis)

**Cel:** Identyfikuj WSZYSTKIE miejsca w core/ i research/ ktore wymagaja
update z V_orig → V_M911 + sqrt(-g) update + κ update.

**Metoda:**
1. Grep dla `V_orig`, `V_TGP`, `\beta/3 \psi^3`, `psi**3/3 - psi**4/4`
2. Grep dla `\kappa = 3/(4\Phi_0)`, `kappa = 3/(4*Phi_0)`, `3/(4Phi_0)`
3. Grep dla `m_sp = γ`, `m_sp^2 = gamma`
4. Grep dla `\sqrt{-g} = c \psi`, `sqrt(-g) = c*psi`
5. Grep dla `M9.1 (FALSIFIED)`, `M9.1` references w core
6. Sklasyfikuj: PHYSICAL (wymaga re-derivation), ANNOTATION (wymaga update tekstu),
   COMMENT (wymaga rewording)
7. Sortuj impact: HIGH (foundational), MEDIUM (derivation), LOW (commentary)

**PASS criteria:**
- Pelna lista plikow z linenumbers (w research/op-g0-..../P33_audit_results.md)
- Klasyfikacja impact dla kazdego match
- Estimate czasu Phase 4 dla kazdego cluster

**Plik:** `phase3_P33_cross_reference_audit.py` (skrypt) + `P33_audit_results.md`

---

### P31 — Sek08a v2.0 specification draft

**Cel:** Napisz `sek08a_v2_specification.md` ktory specyfikuje DOKLADNIE jakie
zmiany sa potrzebne w sek08a, z gotowym draft tekstem dla kluczowych
sub-sections.

**Sections do update (wstepna lista):**

| Section | Update type | Source |
|---|---|---|
| prop:K_psi-uniqueness | unchanged ✓ | Phase 1 |
| prop:V-tgp-canonical | NEW (V_M911) | Phase 1 G0a |
| prop:psi-EOM | replaced (R3 ODE) | Phase 1 G0a |
| prop:vacuum-stability | FIXED (m_sp²=+γ vs old -γ) | Phase 2 P21 |
| prop:kappa-corrected | UPDATED (15/(8Φ_0) vs 3/(4Φ_0)) | Phase 2 P24 |
| eq:cosmo-linearized-unified | UPDATED (5q vs 2q source) | Phase 2 P24 |
| eq:newton-limit | UPDATED (q·Φ_0 = ? c² G_0) | Phase 3 P32 |
| All Φ-EOM derivations | re-done with √(-g) = c·ψ/(4-3ψ) | Phase 3 P33 |
| hyp:vacuum-mass | annotated (G.0 closure) | Phase 2 P21 |

**PASS criteria:**
- Spec document ready z draft text dla wszystkich changed sections
- Backwards compatibility statement (jak relacja stary -> nowy)
- Versioning notes (sek08a v1.x → v2.0)

**Plik:** `sek08a_v2_specification.md`

---

### P34 — Sek08c A1/A2/A3 closure draft

**Cel:** Audit annotations A1, A2, A3 w sek08c (current status: OPEN) closeable
po G.0. Napisz `sek08c_A1_A2_A3_closure_draft.md` z proposed text.

**Annotations do close:**

| Annot | Issue | G.0 resolution |
|---|---|---|
| A1 | Φ-EOM mismatch z R3 ODE | RESOLVED (Phase 1 G0a, V_M911) |
| A2 | √(-g) = c·ψ vs c·ψ/(4-3ψ) inconsistency | RESOLVED (M9.1'' canonical adopted) |
| A3 | 4 different metric forms | RESOLVED (M9.1'' UNIQUE canonical) |

**PASS criteria:**
- Closure text gotowy dla kazdej annotation
- Linkage do Phase 1 + Phase 2 results
- Cross-reference z sek08a v2.0 specification (P31)

**Plik:** `sek08c_A1_A2_A3_closure_draft.md`

---

## 3. Order

**Order:** P32 → P33 → P31 → P34 → Phase3_results.md

| Sub-task | Trudność | Czas | Wymaga |
|---|---|---|---|
| P32 | Średnia-wysoka | ~2 dni | sympy + numerical R3 z source |
| P33 | Niska | ~0.5 dnia | grep + classification |
| P31 | Wysoka | ~2 dni | draft text dla sek08a |
| P34 | Niska | ~0.5 dnia | annotation closure text |
| Phase3_results.md | Średnia | ~0.5 dnia | synthesis |

**Total Phase 3: ~5.5 dni** (~1 tydzien kalendarzowy).

---

## 4. Score gate

```
Score = sum (P3x PASS), x ∈ {1, 2, 3, 4}

≥ 3/4 PASS → Phase 4 forward APPROVED (rekomendacja core mods)
2/4 PASS    → partial; analiza failed sub-tasks; ewentualnie pivot
< 2/4 PASS  → G.0 closure needs revision; back to Phase 1/2
```

**Kill criterion:** Jeśli P32 (Newton limit) FAIL → q·Φ_0 nie reprodukuje
G_0 z V_M911 — to oznacza ze G.0 update jest fundamentalnie incompatible z
empirical Newton's law. To bylby najbardziej powazny finding (G.0 NEGATIVE).

**Najprawdopodobniej PASS** bo z P24 wiemy ze source coupling form jest
preserved (5q/Φ_0 vs 2q/Φ_0), wiec q·Φ_0 = (factor)·4πG_0/c² powinno
zachodzic z odpowiednim factor.

---

## 5. Hard anchors po Phase 3 closure

Po complete Phase 3 mamy gotowe:

| Element | Source | Status |
|---|---|---|
| V_M911 = -γψ²(4-3ψ)²/12 | P1 G0a | LOCK |
| K(ψ) = ψ⁴ | sek08a (zachowane) | LOCK |
| √(-g) = c·ψ/(4-3ψ) | sek08c M9.1'' (zachowane) | LOCK |
| Vacuum ψ=1, m_sp²=+γ | P21 | LOCK |
| Mass spectrum lepton | P22 | LOCK |
| PPN γ=β=1 | P23 | LOCK |
| FRW κ=15/(8Φ_0) | P24 | LOCK |
| **q·Φ_0 NEW value** | **P32** | **TBD** |
| **Φ_0 NEW value** | **P32** | **TBD** |
| Sek08a v2.0 spec | P31 | TBD |
| Sek08c A1/A2/A3 closure | P34 | TBD |
| Cross-reference impact | P33 | TBD |

---

## 6. Plik scaffolding Phase 3

```
research/op-g0-r3-from-canonical-projection/
├── README.md                                       ✓ updated z Phase 2
├── Phase1_setup.md                                 ✓
├── Phase1_results.md                               ✓ CLOSED-POSITIVE
├── phase1_G0a_volume_integration.py                ✓ PASS 4/4
├── phase1_G0b_field_redefinition.py                ✓ NEGATIVE-INFORMATIVE
├── phase1_G0c_einstein_frame_projection.py         ✓ PASS 3/3
├── Phase2_setup.md                                 ✓
├── Phase2_results.md                               ✓ CLOSED-POSITIVE 4/4
├── phase2_P21_vacuum_uniqueness.py                 ✓ PASS 4/5
├── phase2_P22_mass_spectrum_verification.py        ✓ PASS 5/5
├── phase2_P23_PPN_verification.py                  ✓ PASS 5/5
├── phase2_P24_FRW_cosmology.py                     ✓ PASS 5/5
├── Phase3_setup.md                                 ✓ ten dokument
├── phase3_P32_newton_limit_rederivation.py         ← następne
├── phase3_P33_cross_reference_audit.py             ← potem
├── P33_audit_results.md                            ← potem
├── sek08a_v2_specification.md                      ← potem
├── sek08c_A1_A2_A3_closure_draft.md                ← potem
└── Phase3_results.md                               ← finally
```

---

## 7. Open questions z Phase 2 (resolution targets w Phase 3)

1. **Czy q·Φ_0 = 4πG_0/c² nadal trzyma?** — P32
2. **Jakie powinno byc Φ_0_new?** — P32 (z fit Newton G_0)
3. **Ile miejsc w core/ wymaga update?** — P33
4. **Czy A1/A2/A3 sek08c closeable?** — P34 (yes, expected)
5. **Czy sa side effects w innych sekcjach (sek08b, sek09, ...)?** — P33

---

**Status:** Phase 3 ACTIVE — start 2026-05-02 (kontynuacja sesji).

**Order startu:** P32 (Newton limit, MUST-PASS) → P33 → P31 → P34 → P3 results.
