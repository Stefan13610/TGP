---
title: "ADDENDUM 2026-05-10 — Native-observables-first reframing of Phase 2-4 results"
date: 2026-05-10
parent: "[[./README.md]]"
type: addendum
addendum_kind: methodological-reframing
modifies_classification: NO
modifies_sympy_count: NO
status: 🟢 ACTIVE — interpretive overlay (no derivation change)
related:
  - "[[../../meta/PPN_AS_PROJECTION.md]] (parent methodology doc)"
  - "[[./Phase_FINAL_close.md]] §3.2, §3.3 (results being reframed)"
  - "[[../../TGP_FOUNDATIONS.md]] §3.6.2 (target of FOUNDATIONS reframe)"
tags:
  - addendum
  - methodological-reframing
  - native-observables-first
  - PPN-as-projection
  - interpretive-overlay
---

# ADDENDUM 2026-05-10 — Native-observables-first reframing

## §0 — Status uwagi

Ten addendum **NIE zmienia** STRUCTURAL_DERIVED classification cyklu, **NIE zmienia**
sympy 57/57 PASS count, **NIE zmienia** P1-P6 resolution status. Jest **interpretive
overlay** — reframe sposobu prezentacji wyników Phase 2-4, motywowany insightem autora
2026-05-10:

> "γ jest natywne dla TGP, β wydaje się szukaniem na siłę parametru, którego nie ma w
> teorii. Mapowanie na siłę parametru β, który jest tylko zmienną ukrytą TGP zależną od
> samego γ i konfiguracji pola Φ, nie wydaje się mieć sensu."

Insight jest poprawny. Formalizowane w `meta/PPN_AS_PROJECTION.md` (parent doc tego
addendum'a). Niniejszy addendum aplikuje methodology do wyników tego cyklu.

## §1 — Reframing Phase 2 results

### §1.1 — Original presentation (Phase_FINAL_close §3.2)

```
γ_PPN = 1 ⇔ b_1 = -a_1               (1PN constraint)
β_PPN = 1 ⇔ ξ_2 = ξ - a_2·ξ³/2       (2PN constraint)
σ-coupling C(ψ) FREE w 1PN/2PN regime  (enters at 2.5PN)
```

Format: PPN parameter as function of constraint on native Taylor coefs.

### §1.2 — Native-first reframing

**L1 — Native predictions (primary):**

W 1PN/2PN regime (1 dominujące źródło, weak-field), TGP daje obserwable bezpośrednio z
geodezyk w `g_eff[Φ_⊙]`:

| Observable | TGP-native formula | Constraint na native coefs |
|---|---|---|
| Light deflection | θ = 2(1+a_1+b_1)·GM/(bc²) | a_1, b_1 (Cassini constraint) |
| Shapiro delay | Δt ∝ (1+a_1+b_1)·log(...) | a_1, b_1 (Cassini constraint) |
| Perihelion shift | Δϕ = f(a_1, a_2, b_1, b_2, ξ_2, ξ_3) | a_2, b_2, ξ_2 (Mercury, LLR) |
| Nordtvedt effect | η_N = g(a_1, a_2, b_1, ξ_2) | LLR test of EP |

**L2 — PPN projection (consistency check):**

```
γ - 1 = (b_1 + a_1) / 1                      ⟹ γ = 1 ⟺ b_1 = -a_1
β - 1 = (ξ_2 - ξ + a_2·ξ³/2) / 1             ⟹ β = 1 ⟺ ξ_2 = ξ - a_2·ξ³/2
ξ_PPN = 0     (Whitehead)                    ⟹ TGP preserves at 2PN
α_i = 0  (preferred frame, i=1,2,3)         ⟹ FORCED z S05 + §5.1 Lorentz-invariance
ζ_i = 0  (energy/momentum non-conservation) ⟹ FORCED z covariant Φ-EOM
```

**L3 — Falsification map:**

| Bound | Constrains native coefs | Window |
|---|---|---|
| Cassini (γ) | a_1, b_1 (constraint b_1 + a_1 = 0 ± 2.3·10⁻⁵) | preserved EXACT w obecnym ansatz |
| LLR (β) | a_2, b_2, ξ_2 (constraint ξ_2 - ξ + a_2·ξ³/2 = 0 ± 8·10⁻⁵) | preserved EXACT w obecnym ansatz |
| Mercury perihelion | a_2, ξ_2 | preserved (subset of LLR constraints) |

### §1.3 — Co się zmienia

**Logically:** nic (constraint system identyczny).

**Methodologically:** dwie różnice istotne:

1. **Parameter count.** Native: cykl constrainuje ~5 coefs (`a_1, a_2, b_1, b_2, ξ_2`); ξ
   reszta (ξ_3, c_0, κ_σ) deferred. PPN: wygląda jakby było 10 swobodnych parametrów
   (β, γ, ξ, α₁₂₃, ζ₁₂₃₄), z których 4 są forced-zero. Native count ujawnia, że TGP ma
   znacznie mniej swobody niż PPN-language sugeruje.

2. **Falsyfikacja propaguje czysto.** Jeśli kiedyś LLR zmierzy `|β−1| > 10⁻⁵`, native-first
   reframe natychmiast mówi *które coefs* są constrained — bez dwuznaczności "is β a free
   parameter or derived?".

## §2 — Reframing Phase 3-4 results

### §2.1 — Original (Phase_FINAL_close §3.3)

```
β_ppE^new(c_0) at η=1/4:
   β_ppE = (45/16)·Δe_2 + (45/16)·c_0·κ_σ(η)
Δe_2 = -a_1·ξ_3 - 3 - 4·a_2/a_1² + 4·b_2/a_1² - 8·a_3/a_1³ + 16·a_2²/a_1⁴
```

### §2.2 — Native-first reframing

**L1 — Native prediction (primary):**

W 2.5PN binary inspiral regime (η = m₁m₂/(m₁+m₂)² = 1/4 dla equal-mass), TGP daje phase
modification GW przez:

```
Δφ(f) = native function of (a_1, a_2, a_3, b_2, ξ_3, c_0, κ_σ; η, M_chirp)
```

To liczy się **bezpośrednio z g_eff[Φ_1, Φ_2] + Φ-EOM dla 2-body inspiral**, używając
gradient cross-terms `∂_μΦ_1·∂_νΦ_2` jako structurally-new contribution względem
single-source.

**L2 — ppE projection (consistency check):**

```
β_ppE^TGP = (45/16)·[-a_1·ξ_3 - 3 - 4·a_2/a_1² + 4·b_2/a_1² - 8·a_3/a_1³ + 16·a_2²/a_1⁴] 
           + (45/16)·c_0·κ_σ
```

**L3 — Falsification map:**

| Bound | Constrains native coefs | Status |
|---|---|---|
| GWTC-3 \|β_ppE\| ≤ 0.78 (1σ) | Linear combination of `{a_1, a_2, a_3, b_2, ξ_3, c_0·κ_σ}` | Window EXISTS w przestrzeni |
| GWTC-3 + M9.1'' specific values | M9.1'' point `{a_n^M911, ξ_n^M911}` | **FALSIFIED 5σ** (5 sigma na specific point, NIE on space) |
| Recovery Path 1: c_0=0 | ξ_3 = (32 - a_3)/32 (3PN tuning) | preserved (zero-β solution) |
| Recovery Path 2: M9.1'' params | c_0·κ_σ = 4/3 (σ-coupling addition) | preserved (canonical, structurally preferred per §5 Phase_FINAL_close) |

### §2.3 — Kluczowa konsekwencja

Falsyfikacja `β_ppE^TGP = -15/4` przez GWTC-3 (2026-05-09) **NIE jest falsyfikacją TGP**
ani "falsyfikacją β_ppE jako parametru". Jest **falsyfikacją punktu** `{a_n^M911,
ξ_n^M911}` w przestrzeni Taylor coefs `g_eff[Φ]`. Native-first reframe tej obserwacji jest
clean: `M9.1'' point excluded; neighbourhood open`.

W PPN/ppE-only language to wyglądało jak "TGP β predicts -15/4, observation rules out, TGP
fails". Native-first language: "TGP space of `g_eff[Φ]` Taylor expansions has a region with
β_ppE ∈ [-0.78, 0.78]; M9.1'' specific point is outside that region; recovery exists by
moving to Path 1 or Path 2 within the same space."

To jest **bardziej honest ontologically**, bo TGP nie był nigdy "metric theory with β as
free parameter" — był substrate theory with `g_eff[Φ]` functional, którego Taylor expansion
projektuje się na ppE β.

## §3 — Implication dla downstream cykli

### §3.1 — Cykle już zaplanowane post-emergent-metric

| Cykl | Implication z native-first reframing |
|---|---|
| `op-c0-derivation-from-substrate-2026-05-09` | c_0 jest native coef (nie ppE artifact). Derivation z substrate jest *natywna*, NIE "kalibracja parametru ppE" |
| `op-kappa-sigma-2body-PN-2026-05-09` | Tożsamie: κ_σ jest native, derivation z 2-body PN jest natywna |
| `op-scalar-mode-LIGO-bound-2026-05-09` | h_TT^σ = h_TT^GR EXACTLY jest native result (Phase 2 σ-3PN); ppE α-projection jest L2 |
| `op-recovery-V-LIGO-regime-2026-05-10` | Recovery cycles **mogą jasno** mówić "ten coef constrainowany" zamiast "ten parameter passes" |

### §3.2 — Cykle gravity-related po 2026-05-10

Per `meta/PPN_AS_PROJECTION.md §3.1`: **mandatory three-layer presentation** (L1 native, L2
projection, L3 falsification map). Każdy nowy gravity-sector cycle Phase FINAL musi mieć te
trzy warstwy.

### §3.3 — TGP_FOUNDATIONS §3.6.2 reframe

Obecna sekcja `TGP_FOUNDATIONS.md §3.6.2 PPN constraints jako derivation (NIE postulat)`
prezentuje constraint w *inverse* formie (PPN → native). Native-first wersja powinna być
*forward* (native coefs first, projection second).

To **nie blocker** — multi-session work; addendum `meta/PPN_AS_PROJECTION.md §5` notes
implication. Pending dedicated edit cycle.

### §3.4 — audyt/T01_LIGO3G_falsifier reactivation

T01 falsifier framework (currently HANDOFF stub, pre-falsification, β=−5/64 stale) — jeśli
reactivowany, powinien być **native-coefs falsifier**:

```
Native form: "GWTC-X bounds rule out region in (a_1, a_2, a_3, b_2, ξ_3, c_0·κ_σ) space
              of size W_native"
ppE projection: "GWTC-X bounds give |β_ppE| ≤ B"
Map: B-bound translates to W_native via β_ppE = (45/16)·Δe_2 + (45/16)·c_0·κ_σ
```

Pending decision na reactivation.

## §4 — Co addendum NIE zmienia

Explicite, dla unikania niejasności:

- **Sympy 57/57 PASS count:** unchanged. Addendum jest interpretive overlay, nie nowy
  derivation.
- **STRUCTURAL DERIVED classification:** unchanged. Cykl pozostaje closed STRUCTURAL_DERIVED.
- **6/6 P-requirements RESOLVED:** unchanged. P2 "1PN reproduction γ=β=1 z derivation"
  jest preserved — addendum tylko wyjaśnia że "γ=β=1" to L2 form rezultatu, którego L1 form
  to constraints na native Taylor coefs.
- **Phase_FINAL_close.md content:** unchanged. Addendum cross-references je, ale ich nie
  modyfikuje. Jeśli FUTURE wersja Phase_FINAL byłaby przepisana w native-first form,
  byłoby to *separate edit*, nie ten addendum.
- **Six requirements P1-P6 status:** unchanged.

Addendum służy jako *interpretive lens* dla downstream cykli i FOUNDATIONS reframe, nie
jako zmiana derivation.

## §5 — Sign-off

**Addendum authored:** 2026-05-10 (post-conversation autor + Claudian o γ-natywność vs
β-induced).

**Status:** ACTIVE interpretive overlay. Cycle classification preserved.

**Cross-references:**
- `meta/PPN_AS_PROJECTION.md` — parent methodology doc (binding 2026-05-10+)
- `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md §4` — siostrzany form-meaning protocol
  (analogous structure)
- `Phase_FINAL_close.md §3.2, §3.3` — original presentation that this addendum reframes
- `TGP_FOUNDATIONS.md §3.6.2` — target of future native-first reframe edit

**Insight credit:** autor TGP, 2026-05-10. Addendum formalizuje observation jako binding
methodology dla cyklu.
