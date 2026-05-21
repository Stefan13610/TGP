---
title: "L07 — Warunek zerowej sumy: aksjomat vs derywacja"
date: 2026-05-06
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L07_zero_sum_axiom
tags:
  - audit
  - ontology
  - dark-energy
  - zero-sum
  - axiom
  - L07
  - EXT-2
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-L07-zero-sum-Z2-derivation-2026-05-16 (Path A partial closure)", "op-nonlocal-foundations (Path D, deferred multi-session)"]
  depends_on: []
  impacts: ["sek05_ciemna_energia (Λ_eff > 0) — foundation STRENGTHENED 2026-05-16", "ax:zero w sek01_ontologia — ZS1 derived as Z₂-tożsamość 2026-05-16"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2"
    - "[[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/Phase_FINAL_close.md]] (B+ partial closure 2026-05-16)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-16
---

# L07 — Warunek zerowej sumy: aksjomat vs derywacja

## 🟡 STATUS UPDATE 2026-05-16 — **PATH A PARTIAL CLOSURE B+** via op-L07-zero-sum-Z2-derivation

**Closure cycle:** [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/Phase_FINAL_close.md]]
(STRUCTURAL_PARTIAL_DERIVED, 11/11 sympy PASS, 10 FP + 1 LIT, 6/6 P-requirements RESOLVED, 0 hardcoded).

**Disposition update:**

| L07 component | Pre-cycle status | Post-cycle status (2026-05-16) |
|---|---|---|
| **ZS1** (chiralna, ∫Δ√h = 0) | aksjomat (ax:zero w sek01 §417) | ✅ **DERIVED AS Z₂-TOŻSAMOŚĆ** (clean A−) |
| **ZS2** linear part (∝ δφ) | aksjomat | ✅ **DERIVED** (Z₂-orbit balance, parallel ZS1) |
| **ZS2** quadratic part (∝ (δφ)²) | aksjomat | 🟡 **GAUGE FIXING** (Φ₀ ≡ ⟨Φ⟩_Σ boundary condition, NIE raw axiom) |
| **prop:Lambda-positive** foundation | wisi na raw ZS2 aksjomacie | ✅ **STRENGTHENED** (ZS1 + boundary + ⟨δφ²⟩>0) |
| **Cosmological constant problem** disposition | open w TGP | ✅ **Foundations clarified** strukturalnie |

**Key derivation mechanism (ZS1):**
```
H_Γ[φ] = H_Γ[-φ] (Z₂-invariant);  Δ(x) Z₂-odd;  P_Z₂|Ψ⟩ = |Ψ⟩
⇒ ⟨Ψ|Δ(x)|Ψ⟩ = -⟨Ψ|Δ(x)|Ψ⟩ ⇒ ⟨Δ(x)⟩ = 0 pointwise
⇒ ZS1: ∫_Σ ⟨Δ⟩_Ψ √h d³x = 0   ✅ DERIVED AS Z₂-TOŻSAMOŚĆ
Analog: QCD ⟨q̄γ⁵q⟩=0 (Goldstone-Nambu 1960-61)
```

**Key insight (ZS2 character):**
```
Φ(φ) = (φ/v)²·Φ₀ jest Z₂-EVEN — czysta Z₂-tożsamość NIE działa dla ZS2.
δΦ = (2Φ₀/v)·δφ + (Φ₀/v²)·(δφ)²
     |________|    |__________|
     Z₂-orbit       gauge fixing
     vanishes       (Φ₀ ≡ ⟨Φ⟩_Σ definitional)
```

**Path enumeration (post-cycle):**
| Path | Pre-cycle | Post-cycle |
|---|---|---|
| Path A (Z₂-tożsamość) | hypothetical | ✅ **partially successful** (ZS1 clean; ZS2 partial+gauge) |
| Path B (Lagrange multiplier) | alternative | NIE attempted (B+ achieved without) |
| Path C (φ_eff redefinition) | alternative | partially overlapping (T9 boundary condition is C-like) |
| Path D (nonlokalność cosmological) | alternative | **reserved for ZS2 full pure-Z₂-tożsamość** (deferred multi-session) |

**Proposed core update (review-only — `may_edit_core: false`):**

W cycle's `Phase_FINAL_close.md` §5.1 zaproponowane annotations dla:
- `core/sek01_ontologia/sek01_ontologia.tex` ax:zero (§417): cross-link "ZS1 derived as
  Z₂-tożsamość; ZS2 gauge fixing character (op-L07-zero-sum-Z2-derivation-2026-05-16)"
- `core/sek05_ciemna_energia/sek05_ciemna_energia.tex` prop:Lambda-positive (§240-293):
  "Foundation STRENGTHENED — no longer hangs on raw ZS2 axiom; Λ_eff > 0 emerges from
  ZS1 + boundary condition + ⟨δφ²⟩>0"

Core updates pending separate cycle z `may_edit_core: true` authorization.

**Audit problem disposition:** **PARTIAL CLOSURE B+** — pre-registered acceptable outcome.
Path D (nonlokalność full structural derivation of ZS2 quadratic) **REMAINS OPEN** for
multi-session extension.

---

## Klasa: LUKA ONTOLOGICZNA (zewnętrzna recenzja EXT-2) • Priorytet: **P2** → **P2 PARTIAL B+ post-2026-05-16**

## Diagnoza (z EXT-2)

`sek05_ciemna_energia.tex` lin. 240–293 (`prop:Lambda-positive`)
opiera **całość mechanizmu Λ_eff** na warunku zerowej sumy:

```
∫_Σ φ √h d³x = 0     (eq:zero-sum-phi, lin. 250)
```

zwany **"zasadą zerowej sumy"** (`ax:zero` w sek01).

Bez tego warunku ⟨φ²_min⟩ niekoniecznie > 0, a Λ_eff = (8πG_0/c_0⁴)
⟨U(φ_min)⟩ niekoniecznie > 0. **Cała sek05 (ciemna energia z fluktuacji)
zawisa na tym jednym wymaganiu.**

Warunek jest **globalny** — więz na całej hypersurface przestrzennej Σ.
Nie wynika z lokalnej dynamiki (Φ-EOM), nie wynika z H_Γ jawnie.

## Pliki dotknięte

- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]
  lin. 146–161 (subsection "Związek z zasadą zerowej sumy"),
  240–293 (prop:Lambda-positive)
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] aksjomat ax:zero
  (linia do zlokalizowania w pełnym pliku 1435 lin.)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] lin. 736+
  (H_Γ definition v2 GL-bond)
- `axioms/substrat/` (substrat directory — referencje do roli Z₂)

## Problemy jakie rodzi

1. **Lokalność vs globalność.** Globalny więz bez Lagrange'a multiplier
   z dynamiki to **klasyczna czerwona flaga w fizyce teoretycznej**.
   Albo jest konsekwencją gauge/redundancji (jak constraint Hamiltoniana
   w GR z dyfeomorfizmu), albo wprowadza nielokalność w postaci
   "system jako całość zna globalny rozkład materii".

2. **Operacyjna nieobserwowalność.** Warunek dotyczy całej hypersurface
   Σ jednocześnie. Lokalny eksperymentator nie może go zweryfikować ani
   złamać. To czyni go **nie-falsyfikowalnym** — niezgodne z duchem
   TGP (`TGP_FOUNDATIONS.md` § 5.3: "GR jako numeryczny analog, nie
   izomorfizm; otwiera drzwi falsyfikacji").

3. **Mechaniczna ad-hoc-owość.** Wprowadzenie ax:zero jako *aksjomat*
   przyznaje, że jest on potrzebny, ale nie wynika z głębszej
   struktury. TGP ma 3 fundamenty:
   - jedno pole Φ z Z₂ (§1 FOUNDATIONS)
   - hamiltonian H_Γ z GL-bond
   - **warunek zerowej sumy ∫ φ √h = 0**

   Trzeci wygląda na *montowanie* zamiast *wyprowadzenia*.

4. **Kompatybilność z lokalnymi rozwiązaniami.** Jeśli rozważam
   izolowaną cząstkę w TGP (radialny kink Φ(r)), to lokalnie φ > 0
   (rośnie wokół źródła). ∫ φ √h > 0 dla każdej skończonej kuli. Więz
   ∫_Σ φ = 0 wymaga, by **gdzieś indziej** φ był ujemny — ale
   `sek02_pole` lin. 116–117 explicit zakazuje φ < 0 ("Φ(r) > 0 wszędzie,
   przestrzeń jest generowana"). **Sprzeczność?** Albo "deficyt φ < 0"
   oznacza coś innego niż dosłowna ujemność — ale wtedy potrzebne jest
   doprecyzowanie operacyjne.

5. **Cosmological constant problem zwija się tylko częściowo.**
   Λ_eff ~ γ/12 ~ Φ_0 H_0²/(12 c_0²). Wymaga γ w skali H_0² — ale
   *dlaczego* γ jest w skali H_0², a nie M_Pl²? Standardowy CC problem
   (10¹²⁰ rozbieżność QFT vacuum vs obs) został przesunięty z "dlaczego
   Λ jest takie małe" na "dlaczego γ jest takie małe". **Bez
   wyprowadzenia γ z UV, problem wraca tylnymi drzwiami.**

## Potencjalne ścieżki domknięcia

### Ścieżka A — derywacja zerowej sumy z Z₂-symetrii substratu — ✅ **PARTIAL CLOSURE 2026-05-16**

**Hipoteza:** ∫_Σ φ √h = 0 może być konsekwencją tożsamości typu
Bianchi dla pola φ na grafie Γ. Konkretnie: jeśli H_Γ = −J Σ (φ_i φ_j)²
ma chiralną Z₂ (`sek01_ontologia.tex:83`), to operator parzystości
P: φ_i ↔ −φ_i pozostawia H_Γ invariantnym. Suma ∫ φ √h nad
hypersurface zachowującą P powinna być zerem z konstrukcji.

**Cykl wykonany:** [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/]]
(STRUCTURAL_PARTIAL_DERIVED, 11/11 sympy PASS, claim B+).

**Outcome (per Phase_FINAL_close.md):**
- ✅ **ZS1** (chiralna, ∫⟨Δ⟩√h = 0): DERIVED AS Z₂-tożsamość — clean operator-identity
  argument, analog do QCD ⟨q̄γ⁵q⟩=0 (Goldstone-Nambu 1960-61)
- ✅ **ZS2 linear part** (∝ δφ): DERIVED — Z₂-orbit balance (parallel ZS1)
- 🟡 **ZS2 quadratic part** (∝ (δφ)²): GAUGE FIXING (Φ₀ ≡ ⟨Φ⟩_Σ boundary condition);
  NIE pure Z₂-tożsamość (bo Φ = (φ/v)²·Φ₀ jest Z₂-EVEN)
- ✅ **prop:Lambda-positive** foundation STRENGTHENED — nie wisi na raw ZS2 axiom

**Status:** Path A **PARTIALLY SUCCESSFUL** dla ZS1 + ZS2-linear; ZS2-quadratic
zachowuje **gauge fixing character** (NIE raw axiom). Path D nonlokalność rezerwowana
dla full pure-Z₂-tożsamość derivation ZS2-quadratic.

### Ścieżka B — Lagrange'a multiplier w akcji zunifikowanej

Dodać do S_TGP człon λ ∫ φ √h. Wariacja po λ daje warunek; po φ daje
EOM z λ jako efektywnym Λ. To by *zlikwidowało* ax:zero jako aksjomat,
redukując go do dynamicznego więzu z explicit mnożnikiem. Cena: λ
jest dodatkowym parametrem, którego skala wymaga wyjaśnienia.

### Ścieżka C — Wprowadzenie φ_eff = φ − ⟨φ⟩_Σ

Reinterpretacja: "zerowa suma" jest *definicją* φ_eff jako odchylenia
od średniej. Wówczas ∫ φ_eff = 0 jest tożsamością. Λ_eff przestaje
wynikać z warunku, a zaczyna wynikać z ⟨φ⟩_Σ — co jest mierzalne
(tło kosmologiczne). Wymaga przekształcenia całej sek05.

### Ścieżka D — domknięcie przez nielokalność — ✅ **PARTIAL CLOSURE B+ 2026-05-16** (5 sub-paths investigated)

Przyznanie, że TGP jest **nielokalna na skali kosmologicznej** (nie
sprzeczne z relativity jeśli nielokalność jest *spacelike* na
hypersurface), z explicit dyskusją horyzontów i kompatybilności
z FRW. Cykl `op-nonlocal-foundations/`.

**Cykl wykonany:** [[../../research/op-L07-Path-D-nonlocal-foundations-2026-05-16/]]
(STRUCTURAL_PARTIAL, 11/11 sympy PASS, claim B+).

**5 sub-paths investigated (Phase 1 outcomes):**

| Sub-path | Mechanism | Status |
|---|---|---|
| **D1** FRW horizon truncation | Mode cutoff k_max = H_0/c | ❌ OBSTRUCTED — ⟨φ²⟩_(k<k_max) ≈ (H_0)²/4π² > 0 |
| **D2** dS SO(4,1) symmetry | Conformal Killing constraints | 🟡 **PARTIAL** — homogeneity ⟨φ²(x)⟩ = const, NIE = 0 |
| **D3** Bunch-Davies vacuum | (H_0/(2π))²·log scaling | ❌ OBSTRUCTED — explicit ~10⁻⁶⁵ eV² > 0 |
| **D4** Wheeler-DeWitt mini-superspace | H_Ψ\|Ψ⟩ = 0 global constraint | ❌ OBSTRUCTED — equivalent do L07 gauge fixing |
| **D5** Closed-FRW S³ topology | π₃(S³) = ℤ winding | ❌ OBSTRUCTED — π₃ trivial dla real scalar φ ∈ ℝ |

**Outcome:**
- ✅ D2 provides REAL partial structural constraint (homogeneity dla Bunch-Davies vacuum)
- ❌ 4 of 5 sub-paths explicitly OBSTRUCTED z first-principles calculations
- ✅ **ZS2 gauge-fixing character SOLIDIFIED jako canonical** z 4 obstruction proofs
- ⚠ A− NIE achieved — full pure-Z₂-tożsamość derivation of ZS2 quadratic = 0 NIE z standard
  cosmology + QFT tools
- ⚠ HALT-B NIE — D2 partial constraint jest real structural contribution

**Status:** Path D **PARTIALLY SUCCESSFUL** w B+ sensie — provides additional structural
clarity on ZS2 quadratic character. Deeper "full structural" derivation deferred do
multi-session/multi-year efforts: full quantum gravity (beyond mini-superspace), holographic
principle, entropic gravity reformulation.

**Post-2026-05-16 status:** ZS2 quadratic remainder (`(Φ₀/v²)·V_Σ·⟨(δφ)²⟩_Σ > 0`)
**CANONICAL DISPOSITION** jest **gauge fixing** (Φ₀ ≡ ⟨Φ⟩_Σ boundary condition).
D2 partial homogeneity constraint dodaje strukturalnie informację, że dla Bunch-Davies
vacuum w dS, ⟨φ²(x)⟩ jest position-independent, co consistent z gauge-fixing
interpretation. Path D extension cycles (full QG, holographic) **DEFERRED**.

## Rekomendowany priorytet

**P2 — wysoki** (pre-2026-05-16) → **P2 PARTIAL B+** (post-2026-05-16).

**Pre-cycle:** Bez ścieżki A lub B, sek05 (główna predykcja kosmologiczna TGP) opierała
się na aksjomacie, którego status nie był ani derived, ani falsyfikowalny.

**Post-cycle (2026-05-16):** Path A osiągnięty PARTIAL B+ —
- ZS1 derived clean (Z₂-tożsamość);
- ZS2 partial (linear Z₂-derived; quadratic gauge fixing);
- prop:Lambda-positive foundation strengthened (no longer wisi na raw axiom).

Pozostałe open: ZS2 quadratic full pure-Z₂-tożsamość (Path D), lub acceptance
of gauge fixing character jako definitive disposition.

## Powiązanie z istniejącym audytem

**Nowa klasa L07** (proponowana w EXT-2). Brak odpowiednika w obecnym
S/L/D/M. Skala konsekwencji: mechanizm ciemnej energii w TGP.

## Cross-references

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2 — recenzja źródłowa
- [[../README.md]] — indeks audytu
- [[../PRIORITY_MATRIX.md]] — L07 P2 status (post-2026-05-16: PARTIAL B+)
- [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/]] — **closure cycle Path A (B+ partial)**
- [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/Phase_FINAL_close.md]] — closure ceremony details
- [[../../research/op-L07-zero-sum-Z2-derivation-2026-05-16/Phase1_results.md]] — derivation results (11/11 PASS)
- [[../../research/op-L07-Path-D-nonlocal-foundations-2026-05-16/]] — **closure cycle Path D (B+ partial; 5 sub-paths)**
- [[../../research/op-L07-Path-D-nonlocal-foundations-2026-05-16/Phase_FINAL_close.md]] — Path D closure details
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]
  prop:Lambda-positive (foundation strengthened 2026-05-16; proposed annotation
  pending separate core update cycle)
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] ax:zero
  (ZS1 derived as Z₂-tożsamość 2026-05-16; proposed annotation pending separate core update)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] H_Γ definition
- [[../../TGP_FOUNDATIONS.md]] § 5.3 falsyfikacja
- [[../../STATE.md]] — 6th cycle of sesja 2026-05-16 entry
