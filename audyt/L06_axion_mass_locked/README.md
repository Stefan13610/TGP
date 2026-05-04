---
title: "L06 — m_X 'locked' 100 MeV jako fenomenologia"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L06_axion_mass_locked
tags:
  - audit
  - ontology
  - axion
  - m_X
  - phenomenology
  - omega-2
  - omega-3
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../S06_circular_anchors]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["omega4-axion-mass-derivation"]
  depends_on:
    - "[[../S06_circular_anchors]]"
  impacts: []
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7, §D.4, §J"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L06 — m_X 'locked' 100 MeV — fenomenologia

## Klasa: LUKA ONTOLOGICZNA / NIEJAWNY FREE PARAMETER • Priorytet: P2

## Diagnoza

`m_X` (mass of TGP axion-like X-particle) jest centralny dla:

- ω.2 (axion coupling lock) — TT13/TT14 Sagnac SNR predictions
- ω.3 (axion decay constant) — `f_a = N_A · 2π² · M_GUT / E_TGP`
- ψ.1.v2 (substrate-light-acceleration) — Yukawa range L_eff ~ 1 µm
- τ.3 (clock-acceleration) — Λ-scan

W większości plików występuje wartość `m_X = 100 MeV`, oznaczona jako
„locked" i prezentowana jako parametr derived.

## Faktyczny status (z audytu)

Audit § A.7 — **CRITICAL**:

> ω.2 nie locks m_X — locks tylko coupling g, degenerowane z α_em.
> Nazwa cyklu sugeruje "axion mass lock", ale faktycznie ω.2 lockuje
> bezwymiarowy g_anomaly = α·E_TGP/(2π) ≈ 8.30·10⁻³.
>
> Δχ² (winner vs plain α_em) = **0.21** (znacznie poniżej 9 = 3σ).
> PASS verdykt warunkowy na **forward 2027+ data**.
>
> Wartość m_X ≈ 100 MeV pojawia się **tylko** w
> `op-psi1.../Phase1_setup.md:50-56` jako phenomenological input
> (Yukawa range L_eff ~ 1 µm) bez derywacji z UV matching, BBN czy
> RG fixed point.

Audit § D.4 — LOW:

> m_X = 100 MeV jest "design choice" dla SNR experimental, NIE derived;
> default τ.3 lab parameter z WW8 anchor jest m_X = 0.83 MeV.

To znaczy: **w plikach żywych są dwa różne `m_X`**:

- `m_X = 100 MeV` — phenomenological w ψ.1, ω.2 dla TT13 SNR
- `m_X = 0.83 MeV` — faktyczny default τ.3 z g·f_X (g_ω.1=8.3·10⁻³ × f_X=100 MeV)

B7 (`research/op-tau3-substrate-clock-acceleration/B7_greens_function_results.md`)
KEY PHYSICS finding pokazuje, że dla `m_X = 0.83 MeV` w lab przy `L ~ 1 mm`,
**m_X·L ~ 4·10⁹ → heavy regime universal w lab**. Bulk signal ZERO
(lnX uniform inside source) — TT7-TT12 numerical predictions sub-detectable.

## ω.3 closure: A7 option-2 zamknięte tautologicznie

Audit § J (post-ω.3 update) deklaruje A7 `option-2 CLOSED`:

> ω.3 program END (FULL CONVERGENCE 18/18) explicit spełnia A7 option 2
> ("jawnie przyznać że m_X jest wolnym parametrem"):
> - TGP axion classified **ALP** (E-only anomaly, no QCD N anomaly)
> - m_a (≡ m_X) is **FREE PARAMETER** post-ω.3
> - forward-gate ω.4+ structural m_a derivation noted

Ale ω.3 jest sam w 74394a8-polluted set ([[../S06_circular_anchors]]).
`f_a = (N_A·2π²·M_GUT)/E_TGP` używa K_struct fitted z UV.2 — **kaskada
do cyrkularności**. Klasyfikacja TGP axion jako ALP jest poprawnym
acknowledgment, ale ledger entries ZZ1-ZZ6 zostały dodane do counter
856 jako jeśli były niezależnymi predykcjami.

## Strukturalna sprzeczność

Status `m_X` w plikach:

| Lokalizacja | Wartość | Status w pliku | Faktyczny status |
|-------------|---------|----------------|------------------|
| ψ.1 Phase1_setup | 100 MeV | „input" | phenomenological choice |
| ω.2 program.md | implicit (g_axion locked) | „LOCKED" | LIVE PARTIAL Δχ²=0.21 |
| τ.3 default | 0.83 MeV (z g·f_X) | „derived" | inherited from g_ω.1 |
| ω.3 ZZ1-ZZ6 | f_a = 4.85·10¹⁷ GeV | „LOCKED" | LOCKED-CONDITIONAL na UV.2 |
| Audit D.4 | „design choice" | „CLOSED annotation" | acknowledged as free |

## Wpływ na predykcje

- **TT7-TT12** (Sagnac SNR, frontier reach, magnetar polar shift): wszystkie
  inherit `m_X` value. A5-patched + B7 KEY PHYSICS: post-multiplicative
  formula + heavy regime reality, predictions wymagają B7-v2 dedicated
  re-derivation.
- **TT13** (lab Sagnac „WYKONALNY DZIŚ"): WITHDRAWN 2026-05-01 (audit § K)
  — bo opierał się na fitted m_X=100 MeV i niefizycznej parametrize
  E=10¹⁵ V/m + B=1000T.
- **f_a = 4.85·10¹⁷ GeV** (ω.3 ZZ1): LOCKED-CONDITIONAL na UV.2 K_struct.
  Jeśli S06 rollback UV.2 → f_a tez.

## Status w audycie

§ A.7 — **CRITICAL**, oznaczony „CLOSED via ω.3 option 2" (§ J.1).
Faktyczny status: **CLOSED-by-redefinition** — m_X przeklasyfikowane
jako free, ale ω.3 i ω.2 wciąż w rejestrze jako DERIVED.

§ D.4 — LOW, „CLOSED annotation". Faktyczny status: dwa różne `m_X`
nadal koegzystują w plikach (100 MeV phenomenological vs 0.83 MeV default).

## Rekomendacja

Otworzyć dedykowany cykl `op-omega4-axion-mass/` z trzema fazami:

### Phase 1 — explicit m_X = FREE w wszystkich plikach

1. Przejść przez ψ.1, ω.1, ω.2, ω.3, τ.3 i jednoznacznie oznaczyć:
   - `m_X = m_a (FREE PARAMETER post-ω.3)` — nie „locked"
   - Każda predykcja zależna od `m_X` powinna mieć explicit „inherits
     m_X from phenomenological choice X MeV" caveat
2. WW7-WW12, ZZ1-ZZ6 → annotation „LIVE PARTIAL — m_X free, locked
   conditional na ω.4 derivation"

### Phase 2 — structural derivation attempt

Próba wyprowadzenia `m_X` z:

- substrate-action (pełna analiza vacuum solitonu w v2 GL)
- instanton-analog (PQ-like breaking w substrate)
- explicit PQ-breaking term jeśli relevant

Jeśli derywacja nie działa — explicit acknowledgment „m_X jest jednym
z fundamentalnych free parameters TGP", takiego samego rzędu jak g₀^e
czy Φ₀.

### Phase 3 — re-fit predictions

Z m_X jako free param:

- TT7-TT12 prawdopodobnie zostają jako forecasts conditional na m_X
- f_a (ω.3) prawdopodobnie wymaga downgrade do PARTIAL
- ψ.1 Yukawa range L_eff jest free w `1/m_X`

**Estymata:** 4–6 tygodni (Phase 2 jest realnym problemem fizycznym).

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `research/op-psi1-substrate-light-acceleration/Phase1_setup.md` | explicit free m_X |
| `research/op-omega2-axion-coupling-lock/` | downgrade do LIVE PARTIAL |
| `research/op-omega3-axion-decay-constant/` | downgrade do LOCKED-CONDITIONAL |
| `research/op-tau3-substrate-clock-acceleration/` | sync m_X = 0.83 MeV vs 100 MeV |
| `PREDICTIONS_REGISTRY.md` | TT7-TT12 + WW7-WW12 + ZZ1-ZZ6 annotacje |
| nowy: `research/op-omega4-axion-mass/` | nowy cykl |

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L6
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7, §D.4, §J
- [[../S06_circular_anchors]] (ω.3 cyrkularność)
