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
  open_bridges: ["~~omega4-axion-mass-derivation~~ → op-L06-axion-mass-derivation-2026-05-16 (PARTIAL B+ closure)"]
  depends_on:
    - "[[../S06_circular_anchors]]"
  impacts:
    - "ω.3 m_a FREE classification — STRENGTHENED 2026-05-16 z explicit obstruction proofs"
    - "ψ.1 m_X = 100 MeV — confirmed phenomenological SNR choice (NIE derived)"
    - "τ.3 m_X = g·f_X = 0.83 MeV — confirmed algebraic relation z phenomenological f_X"
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7, §D.4, §J"
    - "[[../../research/op-L06-axion-mass-derivation-2026-05-16/Phase_FINAL_close.md]] (B+ partial closure 2026-05-16)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-16
---

# L06 — m_X 'locked' 100 MeV — fenomenologia

## 🟡 STATUS UPDATE 2026-05-16 — **PATH 2 PARTIAL CLOSURE B+** via op-L06-axion-mass-derivation

**Closure cycle:** [[../../research/op-L06-axion-mass-derivation-2026-05-16/Phase_FINAL_close.md]]
(STRUCTURAL_PARTIAL — Paths A-D obstructed; Path E FREE PARAMETER strukturalnie verified;
1 numerical anchor documented; 11/11 sympy PASS, 10 FP + 1 LIT, 6/6 P-requirements RESOLVED, 0 hardcoded).

**Disposition update:**

| L06 component | Pre-cycle status | Post-cycle status (2026-05-16) |
|---|---|---|
| **m_X status** | "locked 100 MeV" / FREE post-ω.3 (mixed) | ✅ **FREE PARAMETER strukturalnie verified** |
| **Path 2** (forward derivation) | unattempted | ✅ **partially successful** (4 paths tested z obstruction proofs) |
| **Path A** (substrate breathing mode V''(1)) | hypothetical | ❌ obstructed (tachyonic + OOM mismatch 10) |
| **Path B** (m_X = g·f_X) | τ.3 derivation chain | 🟡 algebraic z phenomenological f_X |
| **Path C** (dimensional enumeration) | hypothetical | ❌ obstructed (0 derivation hits); ⚠ 1 numerical anchor |
| **Path D** (Coleman-Weinberg) | hypothetical | ❌ obstructed (all 3 cutoff scenarios miss) |
| **Path E** (FREE PARAMETER ack) | audit § A.7 option 2 | ✅ **strukturalnie verified** (Goldstone + S05 + ω.1 emergent) |
| **ψ.1/τ.3 cross-cycle inconsistency** | open | ✅ **dispositioned** as phenomenological choice diversity |
| **ω.4 forward-gate** (from ω.3) | open | **partially closed** by this cycle |

**Key structural argument (Path E confirmation):**
```
L07 cycle (today): H_Γ[φ] = H_Γ[-φ] Z₂-exact substrate symmetry DERIVED
↓
Goldstone theorem: pure-Z₂ symmetric Lagrangian + spontaneous breaking
                   → axion-like excitation = Goldstone (m = 0 strukturalnie)
↓
S05 single-Φ: NO explicit Z₂-breaking term w fundamental TGP (Φ-only Lagrangian)
↓
ω.1 g·φ·F·F̃ coupling: Z₂-EVEN (φ Z₂-odd × F·F̃ Z₂-odd = even)
                       → NIE wprowadza Z₂-breaking sam z siebie
                       → m_X² ~ ⟨F·F̃⟩²·loop BACKGROUND-DEPENDENT
↓
CONCLUSION: m_X NIE constant TGP property; m_X = FREE PARAMETER
            observed values (ψ.1: 100 MeV; τ.3: 0.83 MeV) są SNR-optimization CHOICES
```

**Notable finding — NUMERICAL ANCHOR:**

```
(M_Pl² · H_0)^(1/3) = 6.07·10⁷ eV ≈ 60 MeV

Distance from phenomenological target (100 MeV):
  - factor 1.7 (40% off)
  - within ±0.5 OOM anchor tolerance
  - OUTSIDE ±0.041 OOM derivation tolerance (10% precision required)

Status: NUMERICAL ANCHOR / NUMERICAL COINCIDENCE
  - NIE structural derivation (40% off target)
  - NO known TGP structural mechanism
  - ANALOG L08 e_Euler² classification (PHASE6 §11)
  - Documented dla future investigation (if mechanism found, extension cycle)
```

**Cross-cycle harmonization:**
- ψ.1 m_X = 100 MeV: confirmed phenomenological choice (NIE derived; SNR optimization for 1.97 μm Yukawa range)
- τ.3 m_X = g·f_X = 0.83 MeV: confirmed algebraic relation (NIE structural derivation; f_X phenomenological)
- BOTH values są free choices dla different applications; NIE structural sprzeczność

**Recommendation:** **NO immediate action needed.** Both phenomenological m_X values remain
valid w their respective SNR optimization scenarios. Future housekeeping (low priority):
add annotations to ψ.1/τ.3/ω.2 z reference do L06 closure ("m_X confirmed FREE per L06 2026-05-16").

**Audit problem disposition:** **PARTIAL CLOSURE B+** — pre-registered acceptable outcome.
m_X = FREE PARAMETER **strukturalnie verified**; Paths A-D **obstructed** z explicit proofs;
numerical anchor `(M_Pl²·H_0)^(1/3)` **documented** but NOT derivation. Open: pursuit of
structural mechanism dla numerical anchor (extension cycle territory, low priority).

---

## Klasa: LUKA ONTOLOGICZNA / NIEJAWNY FREE PARAMETER • Priorytet: **P2** → **P2 PARTIAL B+ 2026-05-16**

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

### Phase 2 — structural derivation attempt — ✅ **EXECUTED 2026-05-16 (B+ partial)**

**Realized as:** [[../../research/op-L06-axion-mass-derivation-2026-05-16/]]
(STRUCTURAL_PARTIAL_DERIVED, 11/11 sympy PASS, claim B+).

Próba wyprowadzenia `m_X` z:

- ❌ **Path A (substrate breathing mode V''(1)):** V''(1) = -γ < 0 tachyonic;
  √(M_Pl·H_0) ~ 4 meV ≠ 100 MeV (OOM mismatch 10)
- 🟡 **Path B (m_X = g·f_X = 0.83 MeV):** algebraic, f_X phenomenological
- ❌ **Path C (dimensional enumeration TGP scales):** 0 derivation-level hits;
  ⚠ 1 NUMERICAL ANCHOR documented: `(M_Pl²·H_0)^(1/3) ≈ 60 MeV` (factor 1.7 z target,
  NO known structural mechanism — analog L08 e_Euler² classification)
- ❌ **Path D (Coleman-Weinberg radiative):** Planck cutoff TOO BIG; QCD TOO SMALL;
  f_X circular

**Path E (FREE PARAMETER acknowledgment):** ✅ **CONFIRMED strukturalnie**
- Goldstone theorem applied to L07-derived Z₂-exact substrate → pure axion = Goldstone (m=0 strukturalnie)
- S05 single-Φ: NO explicit Z₂-breaking term in fundamental TGP
- ω.1 g·φ·F·F̃ coupling is Z₂-EVEN → emergent mass background-dependent (NIE constant)
- ψ.1 (100 MeV) i τ.3 (0.83 MeV) confirmed jako PHENOMENOLOGICAL CHOICES for SNR scenarios

**Verdict:** `m_X jest jednym z fundamentalnych free parameters TGP`, takiego samego rzędu
jak inne phenomenological inputs. **Audit § A.7 option 2 endorsed strukturalnie.**

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
| `research/op-tau3-substrate-clock-acceleration/` | sync m_X = 0.83 MeV vs 100 MeV — DISPOSITIONED 2026-05-16 (phenomenological choice diversity, NIE conflict) |
| `PREDICTIONS_REGISTRY.md` | TT7-TT12 + WW7-WW12 + ZZ1-ZZ6 annotacje — deferred housekeeping |
| ~~nowy: `research/op-omega4-axion-mass/`~~ → **EXECUTED 2026-05-16** as [[../../research/op-L06-axion-mass-derivation-2026-05-16/]] |

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L6
- [[../PRIORITY_MATRIX.md]] klaster D — L06 P2 status (post-2026-05-16: PARTIAL B+)
- [[../../research/op-L06-axion-mass-derivation-2026-05-16/]] — **closure cycle Path 2 (B+ partial)**
- [[../../research/op-L06-axion-mass-derivation-2026-05-16/Phase_FINAL_close.md]] — closure ceremony details
- [[../../research/op-L06-axion-mass-derivation-2026-05-16/Phase1_results.md]] — derivation results (11/11 PASS + numerical anchor finding)
- [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/]] — Z₂ structure inherited for Goldstone application (T7)
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7, §D.4, §J
- [[../S06_circular_anchors]] (ω.3 cyrkularność)
- [[../../STATE.md]] — 7th cycle of sesja 2026-05-16 entry
