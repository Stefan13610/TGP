---
title: "L01 — operacyjna definicja ρ (warstwa 3b)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L01_rho_operational
tags:
  - audit
  - ontology
  - rho
  - L_mat
  - matter-hierarchy
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../S04_metric_coupling_axiom]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["rho-stress-energy-bridge"]
  depends_on: []
  impacts:
    - "[[../S04_metric_coupling_axiom]]"
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §C.3, §N.2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L01 — operacyjna definicja ρ (warstwa 3b)

## Klasa: LUKA ONTOLOGICZNA • Priorytet: P2

## Diagnoza

[[../../TGP_FOUNDATIONS.md]] §4 rozróżnia trzy warstwy materii:

| Warstwa | Co to jest | Sprzęga z Φ przez | Status |
|---------|-----------|-------------------|--------|
| 3a — Pola materii ψ_m | Klasyczne i kwantowe pola SM | g_eff^μν[Φ] w L_mat, **nie** Φ bezpośrednio | (P) |
| **3b — Gęstość ρ** | Skalarna gęstość w `L_mat = -(q/Φ_0)·φ·ρ` | Bezpośrednie sprzężenie minimalne (czynnik φ) | (W) — to jest **źródło** dla Φ-EOM |
| 3c — Kinki / defekty | Hipoteza fermionów jako struktur w Φ | (otwarty problem) | (Hipoteza) |

**Problem:** `ρ` w warstwie 3b jest *operacyjnie nieokreślone*.

Pojawia się jako „skalarna gęstość", ale brakuje formalnej definicji
kowariantnej:

- Czy `ρ = T^μ_μ/c²` (śląd stress-energy tensora)?
- Czy `ρ = ρ_rest` (gęstość masy spoczynkowej)?
- Czy `ρ = |ψ_D|²·m` (norma pola Diraca × masa)?
- Czy `ρ = ε/c²` (gęstość energii)?

Dla pól nierelatywistycznych: `T^μ_μ = -ρ_rest·c²` (z minusem od mostly-plus
metryki) ⇒ obie definicje zgadzają się do znaku. **Ale dla pól
relatywistycznych:**

- Foton (EM): `T^μ_μ = 0` ⇒ `ρ_EM = 0` ⇒ **fotony nie sprzęgają z TGP**?
- Fermion bez masy: `T^μ_μ = 0` ⇒ analogicznie?
- Pole skalarne ze złamaną symetrią: `T^μ_μ ≠ 0` ale wymaga explicit
  reguł.

## „Rozwiązanie" audytu

§ C.3 zamyka *via* § N (Option-2 dla S04):

> ρ ≡ T^μ_μ / c_0² (relativistic invariant, scalar density)
>
> Dla materii nierelatywistycznej T^μ_μ = -ρ_rest·c² → ρ = -ρ_rest
> (z minusem od η^00 mostly-plus). Definicja kompatybilna z polem ψ_D
> w QFT (T^μ_μ = m·|ψ|²) oraz z stress-energy continuum.

**To jest decyzja konwencyjna**, nie wyprowadzenie. Brakuje:

1. Formalnej spójności ze wszystkimi sektorami SM (foton — czy `ρ=0`
   znaczy „brak sprzęgania" czy „minimal coupling przez g_eff^μν tylko"?).
2. Sprawdzenia, czy ta definicja `ρ` jest stabilna w pętlach kwantowych
   (renormalizacja `T^μ_μ` w QFT na curved background — anomalia śladu).
3. Zgodności z N-body symulacjami (SPARC galaxy fits, gdzie `ρ_baryon`
   jest *ad hoc* identifikowane bez explicit T^μ_μ derivation).

## Wpływ na predykcje

- **GW170817** (`c_GW = c`): jeśli fotony sprzęgają tylko przez `g_eff^μν`,
  to dyspersja sygnału EM jest jak GR. Jeśli też przez `φ·ρ_EM`, to
  `ρ_EM = 0` (T^μ_μ = 0) i nic się nie zmienia. To OK, ale wymaga
  formal acknowledgment.
- **N-body / SPARC**: aktualne fits używają `ρ = ρ_baryon` bez explicit
  derivation z `T^μ_μ`. Jeśli to nie jest spójne, fits są conditionalne.
- **Test piątej siły** (S04): `dφ/dr · ρ` zależy od jednoznacznej
  definicji `ρ`. Bez niej, MICROSCOPE 10⁻¹⁵ test jest formalnie nie
  dobrze postawiony.

## Status w audycie

§ C.3 — MEDIUM, oznaczony „CLOSED via A4" (przez Option-2). Faktyczny
stan: **CLOSED-by-decision-only**, formalnie nie wykazane.

## Rekomendacja

Krótki dedykowany cykl `op-rho-stress-energy-bridge/`:

### Phase 1 — formal definition

Lock `ρ = -T^μ_μ/c²` (mostly-plus signature) z explicit:

- mapping na kanoniczne sektory SM (Dirac, EM, scalar, gauge)
- treatment edge cases (massless fields, conformal coupling)
- spójność z N-body (`ρ → ρ_baryon` w ścisle określonym limicie)

### Phase 2 — quantum corrections

Trace anomaly w QFT na curved background — czy `ρ` jest stabilne pod
renormalizacją? Standard QFT result: `T^μ_μ ≠ 0` jest klasycznie zero
dla conformal fields, ale anomalia trace daje `~R²` corrections.

### Phase 3 — N-body consistency

Re-derive SPARC fits z explicit `T^μ_μ → ρ` mapping. Sprawdzić, czy
zgodność z obserwacjami się utrzymuje.

**Estymata:** 3–4 tygodnie.

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `sek08a_akcja_zunifikowana.tex` lin. 79–83 | dodać explicit ρ = -T^μ_μ/c² |
| `TGP_FOUNDATIONS.md` §4 (warstwa 3b) | dodać formal definition |
| `core/sek08_formalizm.tex` (rem:materia-hierarchia) | sync |
| `nbody/` skrypty SPARC | dodać explicit T^μ_μ derivation |
| nowy: `research/op-rho-stress-energy-bridge/` | nowy cykl |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L1
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §C.3, §N.2
- [[../S04_metric_coupling_axiom]] (powiązany — wymaga ρ definition)
