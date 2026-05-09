---
title: "NEEDS — op-Phi-vacuum-scale-2026-05-09"
date: 2026-05-09
type: needs-list
status: 🔒 CLOSED — final status assigned per P-item
parent: "[[./README.md]]"
final_close: "[[./Phase_FINAL_close.md]]"
tags:
  - needs
  - phase0
  - phi-vacuum-scale
  - cycle-closed
---

> **CYCLE CLOSED 2026-05-09.** Final P1-P13 status w priority matrix poniżej +
> §7 [[./Phase_FINAL_close.md]] dla pełnego retrospektywy.

# NEEDS — explicit wymagania dla cyklu op-Phi-vacuum-scale-2026-05-09

> **Note:** ta lista jest **post-inventory** — odzwierciedla że duża część
> oryginalnego scope (sub-problem B, P3, P6) jest **pre-closed** w prior
> cycles. Realny scope cyklu zawęża się do P1+P4 + dokumentacja.

## CRITICAL needs (P1-P2)

### P1: Φ_0 absolute eV scale dla MAG Phase 5 m_e predictivity
**Priority:** CRITICAL — **MAIN open problem**
**Phase:** 1 (reconnaissance), 2-3 (jeśli CONTINUE)
**Status:** OPEN

**Description:** Phase 5 MAG formula:
```
m_Mach = (3γq²)/(16π Φ_0² m_C) · ⟨δΦ²_bg⟩
```
wymaga **absolute eV value** Φ_0 dla quantitative m_e prediction.
Phase 5 quantitative test pokazał:

| Scenariusz | Φ_0 (eV) | Self-consistency |
|---|---|---|
| (a) H₀ cosmological | 1.5×10⁻³³ | FAILS perturbative (δ_bg ~ 10²² · Φ_0) |
| (b) **v_EW** | 2.46×10¹¹ | **REASONABLE** (δ_bg ~ 0.4% Φ_0) |
| (c) M_Pl | 1.22×10²⁸ | OK ale Φ_0 ~ M_Pl unusual |
| (d) atomic 1 eV | 1.0 | UNPHYSICAL |
| (e) m_e Compton | 5.11×10⁵ | UNPHYSICAL |

**Resolution paths (Phase 1 candidates A1-A6):**
- A1 [√(Λ/G)] ≡ A1 → H₀ algebraicznie (z T-Λ M_Pl²·H₀² ~ ρ_obs)
- A2 [v_EW = 246 GeV] — **scenariusz (b)** Phase 5 working
- A3 [m_e Compton] — definicyjnie circular dla m_e prediction
- A4 [M_Pl] — scenariusz (c) Phase 5
- A5 [M9.1'' geometric ψ=4/3] — dimensionless, nie absolute eV
- A6 [FRG NGFP] — UV.1 daje structural ale nie absolute scale

**Phase 1 quick verdict:** structural compatibility każdego, NIE pełne derivation.

### P2: Sub-problem B (UV/IR Φ_0 reconciliation)
**Priority:** DOCUMENTED — **PRE-CLOSED w UV.3**
**Phase:** 1 (documentation only)
**Status:** ✅ CLOSED via Z_Φ = 14/3 (op-uv3-phi0-renormalization, 16/16 PASS)

**Description:** User flagged "z UV normalization mamy 2 różne wartości Φ_0".
Inventory wskazuje że to literalnie odnosi się do UV.3 result:

```
Φ_0^bare ≈ 115     (UV side, dimensionless = 168·Ω_Λ)
Φ_eff   ≈ 24.65   (IR side, dimensionless = 36·Ω_Λ)
Z_Φ     = 14/3   (algebraic = V(1)/P(1) wave-function renormalization)
```

**Resolution path:** Phase 1 dokumentuje cross-reference do UV.3, NIE re-derivation.

**Open frontiers post-UV.3** (NIE w scope tego cyklu):
- Mechanizm dynamiczny eksponentów (7,8,3,4) w P, V (UV.3 §"NIE rozwiązuje" #1)
- Z_Φ vs UV.1 η_N* = -2 unification (UV.3 §"NIE rozwiązuje" #2)
- γ.1 trade-off 0.88% Ω_Λ ↔ α_s first-principles (UV.3 §"NIE rozwiązuje" #3)

## HIGH needs (P3-P4)

### P3: Connection do Λ_CC (cosmological constant)
**Priority:** DOCUMENTED — **PRE-CLOSED w T-Λ**
**Phase:** 1 (documentation only)
**Status:** ✅ CLOSED (closure_2026-04-26/Lambda_from_Phi0, 7/7 PASS)

**Description:** T-Λ pokazał:
```
ρ_vac,TGP = V(Φ_eq) = γ·Φ_eq²/12
         z Φ_eq = H₀, γ = M_Pl²·g̃ (g̃ ≈ 1)
         = M_Pl²·H₀²/12 ≈ ρ_vac,obs (2% match)
```

**To znaczy że Φ_eq (cosmological vacuum scale) = H₀ JEST TGP-natywnie ustalone**.

**TENSION zidentyfikowana:** Phase 5 MAG scenariusz (a) [Φ_0=H₀] FAILS.
To sugeruje że **mass-generation Φ_0 ≠ cosmological Φ_eq**, lub że MAG
Phase 5 formula wymaga subleading correction (per Phase5b interpretation).

**Resolution path Phase 1:** dokumentacja existing T-Λ closure i explicit
notowanie tej tension jako open question (NIE rozwiązanie w niniejszym
reconnaissance cycle).

### P4: Connection do v_EW
**Priority:** HIGH — **MAIN structural question**
**Phase:** 1-2
**Status:** PARTIAL (δ.2 sketch; pełen mechanism open)

**Description:** Phase 5 MAG scenariusz (b) [Φ_0 = v_EW = 246 GeV] daje
self-consistent reproduction m_e=511 keV. Pytanie: **dlaczego v_EW jest
TGP-natywną Φ_0**?

**Existing TGP-native v_EW closure (δ.2):**
```
v_W = ℓ_P · exp(-4π²/(3·J_EW²)) ≈ 246.2 GeV  (Coleman-Weinberg w sek09 §O14)
```

**To jest TGP-natywne** (TGP EWSB via Coleman-Weinberg), ale formula
wymaga J_EW jako input. δ.2 nie domyka pełnego J_EW derivation.

**Resolution path Phase 1:** dokumentacja δ.2 status; verdict: **A2 (Φ_0 ~ v_EW)
jest najobiecujący kandydat** o ile δ.2 zostaje utrzymany.

## MEDIUM needs (P5-P6)

### P5: Hierarchia α/G_N — out of scope
**Priority:** OUT_OF_SCOPE
**Phase:** N/A
**Status:** Falsified separately w op-MAG-Lorentz Phase 3

**Description:** Original README listed P5 jako requirement ale to MISNOMER.
Hierarchia α/G_N (40 dekad EM ≫ G_N) jest **osobnym** problem:
- op-MAG-Lorentz Phase 3 wykazał że TGP-natywne mechanism dla F_EM/F_grav
  nie wynika z direct A_μ coupling
- To **nie** jest naturalny scope op-Phi-vacuum-scale

**Resolution path Phase 1:** acknowledgment że P5 jest **out of scope**,
listed dla completeness ale NIE addressed.

### P6: Particle spectrum compatibility (m_μ, m_τ, kwarki)
**Priority:** DOCUMENTED — **PRE-CLOSED w P4 + Phase5b**
**Phase:** 1 (documentation only)
**Status:** ✅ CLOSED

**Description:**
- P4 (particle_sector_closure) zamknął lepton tower **niezależnie** od Φ_0
  (Koide K=2/3 z A_tail solitonowej amplitudy + PDG m_τ/m_e)
- Phase5b consistency check (sympy 2/2 PASS) pokazał że MAG Mach formula
  jest **subleading** (m ∝ A²) vs P4 intrinsic energy (m ∝ A⁴ lub A³)
- Lepton tower jest TGP-natywnie domknięty bez Φ_0 fixing

**Resolution path Phase 1:** dokumentacja pre-closure, NIE re-derivation.

## SUPPORTING needs (P7-P9)

### P7: Anti-tautology gate dla Phase 1 candidate scan
**Priority:** CRITICAL (CALIBRATION_PROTOCOL binding)
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** Phase 1 sympy verification MUSI explicit prevent:
- Multi-candidate fit z minimum drift selection
- Constructed criterion to select winner
- Algebraic re-arrangement masquerading as second path
- Definitional tautology
- Sympy-rationalization "DERIVED" without first-principles

**Resolution path:** sympy script Phase 1 zawiera explicit check() helper
z PASS/FAIL counter, każdy A1-A6 candidate oceniany przez **wyprzedzająco
zdefiniowane** kryteria (struktura V(Φ), zgodność z T-Λ, zgodność z UV.3,
zgodność z δ.2 EWSB).

### P8: Honest verdict CONTINUE / STRUCTURAL CONDITIONAL / EARLY_HALT
**Priority:** HIGH
**Phase:** 1 (results)
**Status:** OPEN

**Description:** Po Phase 1 reconnaissance, cykl REKOMENDUJE jedną z opcji:
- **CONTINUE Phase 2-N:** jeśli specific candidate (najpewniej A2 ~ v_EW)
  ma reasonable prospect closure poprzez δ.2 follow-up
- **STRUCTURAL CONDITIONAL halt:** jeśli landscape mapping clarifies framework,
  ale full Φ_0 derivation nie achievable bez deep RG/EWSB infrastructure
- **EARLY_HALT:** jeśli problem strukturalnie nie wykonuje się na poziomie
  reconnaissance, defer until more infrastructure (FRG tools, etc.)

**Resolution path:** Phase 1 results z explicit decyzją + uzasadnieniem.

**User explicit guidance:** "moja decyzja" o starcie pełnego cyklu;
reconnaissance pomaga TĘ decyzję informować.

### P9: Cross-references map
**Priority:** MEDIUM
**Phase:** 1 (documentation)
**Status:** OPEN

**Description:** Cykl ma cross-link do co najmniej:
- UV.3 (sub-problem B closure)
- T-Λ (P3 closure)
- P4 (P6 closure)
- Phase5b (consistency MAG vs P4)
- δ.2 (v_W EWSB sketch)
- δ.1 (g̃ N_f QCD)
- γ.1 (Φ_eff algebraic)
- S07 audyt (M9.1'' postulate)
- M9.1'' canonical V_M911

**Resolution path:** Phase 1 results + balance sheet zawierają explicit
cross-references; wikilinks używane.

## CRITICAL needs added post user iteration (P11-P13)

### P11: Canonical V_M9.1'' vs deprecated V_orig — framework consistency audit
**Priority:** CRITICAL — **BLOCKER dla całego cyklu**
**Phase:** 1.6 → audit chain → Path C confirmed
**Status:** ✅ **FULLY RESOLVED** (post Path C confirmation 2026-05-09)
**Final verdict:** **DUAL-V STRUCTURE** — V_M9.1'' (gravity) + V_orig (matter)

**Audit chain:**
1. [[../op-V-canonical-consistency-audit-2026-05-09/]] — 2 residual gaps zlokalizowane
2. [[../op-MAG-Phase5-V-reference-clarification-2026-05-09/]] — Path C suggested
3. [[../op-dual-V-structure-clarification-2026-05-09/]] — **Path C CONFIRMED** (sympy 10/10 PASS)

**Gravity-only deprecation:** G.0 closure A4 marker (Phase1_results.md linia 266)
explicit reserved matter coupling separate verification. Niniejsze audit chain
formalnie zrealizowało A4 marker.

**Framework status post-resolution:**
- V_M9.1'' canonical = **gravitational sektor** (G.0 P21 LOCKED)
- V_orig canonical = **matter sektor** (A4 verification 2026-05-09)
- T-Λ + Phase 5 V_orig usage: ✅ legitimate (matter sektor)
- sek08a annotation: ⚠️ wymaga update (gravity-only deprecation clarification)

**Description:** Phase 1.6 user iteration (sympy 15/15 PASS) ujawniła że:
- V_orig(Φ) = (β/3)Φ³ - (γ/4)Φ⁴ jest **DEPRECATED 2026-05-02** (sek08a linie 95-110)
- Canonical TGP używa: V_M9.1''(ψ) = -γψ²(4-3ψ)²/12 (quartic w ψ, brak β param)
- Phase 1 reconnaissance (T13) cytowała **deprecated formułę**

**Implikacje:**
- Phase 1 sympy T13 ("β=γ tension") jest **artefaktem deprecated V_orig**
- W canonical V_M9.1'' nie ma β/γ ratio do tuningowania
- Multi-vacuum strukture: ψ ∈ {0, 2/3, 4/3} **zamiast** unique vacuum w V_orig

**Resolution path:**
- Audit każdego cyklu (T-Λ, Phase 5 MAG, UV.3, particle_sector) dla użycia
  V_orig vs V_M9.1''
- Spawn `op-V-canonical-consistency-audit-YYYY-MM-DD`
- Cross-reference: [[./Phase1_6_strong_field_canonical_sympy.py]]

### P12: Multi-vacuum identification — który ψ odpowiada Phase 5 v_EW?
**Priority:** CRITICAL — **kluczowy dla zrozumienia hierarchy v_EW/H_0 ~ 10⁴⁴**
**Phase:** 1.6 discovery, wymaga follow-up cycle
**Status:** OPEN

**Description:** Canonical V_M9.1'' ma TRZY punkty krytyczne:

| ψ | V(ψ) | Interpretacja |
|---|---|---|
| 0 | 0 | Trivial vacuum (asymptotic free?) |
| **2/3** | **-4γ/27** | **GLOBAL MINIMUM** (cosmological vacuum?) |
| 4/3 | 0 | Horyzont M9.1'' (degenerate boundary, "silne pole") |

**Open questions:**
- (a) Czy Phi_eq = (2/3)·Phi_0 = H_0 (z T-Λ) → Phi_0 = (3/2)·H_0?
- (b) Czy v_EW odpowiada **innemu** punktowi krytycznemu (np. ψ=4/3 horyzontowi
  jako "EW vacuum" w innym sensie)?
- (c) Czy istnieją **dwie fizyczne fazy** vacuum TGP (cosmo + EW) oddzielone
  przez running coupling lub phase transition?

**Resolution path:**
- Spawn `op-multi-vacuum-identification-YYYY-MM-DD`
- Pre-requisite: P11 audit (clarify który framework component używa którego V)
- Connection do A6 (FRG NGFP) — multi-vacuum może być natural FRG result

### P13: Phase 5 MAG Phi_0 explicit definition audit
**Priority:** HIGH — wymaga przed jakimkolwiek "spawn op-EWSB-from-substrate"
**Phase:** 1.5 discovery, audit-style
**Status:** OPEN

**Description:** Phase 5 MAG formula:
```
m_Mach = (3γq²)/(16π Phi_0² m_C) · ⟨δΦ²_bg⟩
```

uses "Phi_0" jako parameter, ale **nie jest jasne** czy to:
- (a) Lagrangian parameter z V(Phi) (UV/bare)
- (b) Vacuum value Phi_eq (IR/renormalized)
- (c) Effective scale w EW regime
- (d) Coś innego (np. mass scale m_C)

**Phase 5 sympy** wykazał scale-invariance: delta_bg/Phi_0 ratio jest
**niezależny od Phi_0** (formula struktur. skanije Phi_0² i delta_bg²
proporcjonalnie). To znaczy "scenariusz b reproducuje m_e" jest prawdziwe
**dla każdego Phi_0** o ile delta_bg² odpowiednio dobrane.

**Open question:** czy Phase 5 MAG **faktycznie wymaga** Phi_0 = v_EW,
czy to jest tylko **wybór parametryzacji** delta_bg²?

**Resolution path:**
- Audit [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]]
- Sprawdź derivation Phase 5 z V_M9.1'' (po P11 audit)
- Cross-reference: [[./Phase1_5_user_iteration_collective_sympy.py]] T15

## Priority matrix (updated post user iterations)

| ID | Title | Priority | Phase | Status | Resolution |
|---|---|---|---|---|---|
| P1 | Φ_0 absolute eV scale | CRITICAL | 1+ | OPEN | A1-A6 candidate scan (Phase 1 done) |
| P2 | UV/IR reconciliation | DOCUMENTED | 1 | ✅ CLOSED (UV.3) | Cross-ref |
| P3 | Λ_CC connection | DOCUMENTED | 1 | ✅ CLOSED (T-Λ) — needs P11 audit | Cross-ref + audit |
| P4 | v_EW connection | HIGH | 1-2 | PARTIAL (δ.2) | Document + verdict |
| P5 | α/G_N hierarchy | OUT_OF_SCOPE | N/A | Falsified separately | Acknowledgment |
| P6 | Particle spectrum | DOCUMENTED | 1 | ✅ CLOSED (P4+Phase5b) | Cross-ref |
| P7 | Anti-tautology gate | CRITICAL | 1 | ✅ CLOSED | Sympy explicit checks (37/40 PASS) |
| P8 | Honest verdict | HIGH | 1 | ✅ CLOSED | STRUCTURAL CONDITIONAL HALT |
| P9 | Cross-references | MEDIUM | 1 | ✅ CLOSED | Wikilinks complete |
| **P11** | **V_M9.1'' canonical audit** | **CRITICAL (BLOCKER)** | **1.6** | **OPEN** | Spawn `op-V-canonical-consistency-audit` |
| **P12** | **Multi-vacuum identification** | **CRITICAL** | **1.6** | **OPEN** | Spawn `op-multi-vacuum-identification` |
| **P13** | **Phase 5 Phi_0 audit** | **HIGH** | **1.5** | **OPEN** | Audit Phase 5 MAG derivation |

## Reduced scope post-inventory

**Original ambitious scope (per README):** 6 candidate paths, 2 sub-problems.

**Realistic Phase 1 scope (post-inventory):**
1. **Document** P2, P3, P6 jako pre-closed (UV.3, T-Λ, P4+Phase5b)
2. **Quick scan** A1-A6 z anti-tautology gates (P1, P4, P7)
3. **Honest verdict** CONTINUE / STRUCTURAL CONDITIONAL / EARLY_HALT (P8)
4. **Cross-reference map** (P9)

**To jest 1-session deliverable**, NIE multi-week ambitious cycle.

## Cross-references

- [[./README.md]] — cycle scoping
- [[./Phase0_balance.md]] — balance sheet
- [[./Phase1_reconnaissance_sympy.py]] — sympy verification
- [[./Phase1_reconnaissance_results.md]] — verdict
- [[../op-uv3-phi0-renormalization/README.md]] — P2 closure
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] — P3 closure
- [[../particle_sector_closure/README.md]] — P6 closure
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] — Phase 5 MAG
- [[../op-delta2-Nf-derivation/]] — δ.2 v_W EWSB
- [[../../meta/CALIBRATION_PROTOCOL.md]] — anti-tautology binding
