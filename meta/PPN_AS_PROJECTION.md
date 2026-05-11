---
title: "PPN/ppE jako projekcja, nie fizyka — native-observables-first methodology"
date: 2026-05-10
type: meta-protocol
status: 🟢 ACTIVE — methodological note
binding_scope: "All TGP cycles touching gravity tests (PPN, ppE, cosmology), 2026-05-10+"
related:
  - "[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §4 (form-meaning mapping — analogous protocol for BD-form)"
  - "[[CALIBRATION_PROTOCOL.md]]"
  - "[[../TGP_FOUNDATIONS.md]] §3.6.2 (PPN constraints — to be reframed per this doc)"
  - "[[../research/op-emergent-metric-from-interaction-2026-05-09/]] §3.2 (γ=β=1 derivation — primary application)"
parent: "[[README.md]]"
tags:
  - meta
  - methodology
  - PPN
  - ppE
  - native-observables-first
  - anti-projection-confusion
---

# PPN_AS_PROJECTION — native-observables-first methodology

## §0 — Po co ten dokument

### §0.1 — Insight (autor, 2026-05-10)

> "γ jest natywne dla TGP, β wydaje się szukaniem na siłę parametru, którego nie ma w teorii.
> Mapowanie na siłę parametru β, który jest tylko zmienną ukrytą TGP zależną od samego γ
> i konfiguracji pola Φ, nie wydaje się mieć sensu."
>
> — autor TGP, 2026-05-10

Insight jest poprawny i ma istotne konsekwencje metodologiczne dla framework'a.

### §0.2 — Krótka diagnoza

PPN (Will, Nordtvedt) został zaprojektowany **dla teorii metrycznych** — gdzie `g_μν` jest
obiektem fundamentalnym, a parametry {β, γ, ξ, α₁, α₂, α₃, ζ_i} to współczynniki Taylor
expansion `g_μν` w bazie potencjału `U`. W tym frameworku β i γ są równoprawne.

W TGP **fundamentem nie jest g_μν**, tylko Φ + emergent functional `g_eff[Φ]`. Wtedy:

- **γ** ↔ leading-order odpowiedź `g_ij[Φ]` na perturbację Φ — jedna pochodna funkcjonału.
  **Natywny** obiekt teorii.
- **β** ↔ kombinacja drugich pochodnych Taylor expansion `g_eff[Φ]`. **Induced quantity** —
  projekcja na bazę PPN, która ma własną nazwę u Willa, ale nie ma natywnego sensu w TGP.
- **α_i (preferred-frame)** ↔ parametry Lorentz-violation. **Forced ≡ 0** strukturalnie
  z Lorentz-invariance substratu. Nie są swobodne — nawet nie są obliczalne, są
  tożsamością.

### §0.3 — Konsekwencja: PPN to język tłumaczenia, nie fizyka

**Prawdziwymi przewidywaniami TGP są obserwable** (deflekcja, Shapiro, perihelion, Nordtvedt,
range residual, GW dispersion, ...), liczone bezpośrednio z `g_eff[Φ]` + Φ-EOM. PPN
parameters to *projekcja* tych obserwabli na chart Willa — wygodne dla porównania z
literaturą bound'ów, ale **nie są pierwotnymi obiektami teorii**.

Analogia z `TGP_NATIVE_COMPUTATIONAL_PATTERNS.md §4`: tak jak BD-form formuły mogą mieć
TGP-meaning (form-meaning split), tak PPN basis może mapować się na TGP native bez
implikacji że TGP *jest* metric theory.

---

## §1 — Klasyfikacja PPN parametrów względem TGP native objects

| PPN | TGP-native counterpart | Klasa | Status |
|---|---|---|---|
| **γ** | `δg_ij/δΦ` leading-order | **NATYWNY** (1-st derivative funkcjonału) | 1-1 mapping; γ=1 ⟺ `b_1 = −a_1` (Phase 2 emergent-metric) |
| **β** | Kombinacja `{a_2, ξ_2, ξ_3, b_2, a_1²}` | **INDUCED** (mix 2nd derivatives) | β=1 ⟺ `ξ_2 = ξ − a_2·ξ³/2` (Phase 2 emergent-metric) |
| **ξ (Whitehead)** | Higher-order Taylor coef of `g_tt[Φ]` | INDUCED | TGP w 2PN preserves ξ=0 |
| **α₁, α₂, α₃** | Preferred-frame breaking term | **FORCED ≡ 0** | Strukturalnie, z Lorentz-invariance Φ-substratu |
| **ζ_i (energy/momentum non-conservation)** | Non-conservation w T_μν | **FORCED ≡ 0** | Strukturalnie, z covariant Φ-EOM |

**Obserwacja:** TGP nie ma 10 swobodnych PPN parametrów. Ma kilka Taylor coefs `{a_n, ξ_n,
b_n, c_0, κ_σ}` emergent functional — zazwyczaj 3-5 niezależnych. Reszta PPN parametrów
albo *projektuje się* na te coefs (γ, β, ξ), albo jest *strukturalnie zerowa* (α_i, ζ_i).

---

## §2 — Tabela analogiczna dla ppE (post-2.5PN regime GW)

ppE basis (Yunes & Pretorius) jest *innym* projection chart'em — dla phase modification
inspiral GW. Ten sam paradygmat się stosuje:

| ppE | TGP-native counterpart | Klasa | Status |
|---|---|---|---|
| **β_ppE (phase deviation)** | Kombinacja `{a_1, a_2, a_3, b_2, ξ_3, c_0, κ_σ}` z 2.5PN binary inspiral | **INDUCED** | β_ppE^new = (45/16)·Δe_2 + (45/16)·c_0·κ_σ (Phase 4 emergent-metric) |
| **α_ppE (amplitude deviation)** | Kombinacja TT-projection σ-coupling C(ψ) z gradient-strain | INDUCED | h_TT^σ = h_TT^GR EXACTLY (T3.4 amendment + σ-3PN Phase 2) |
| **b_ppE (PN order index)** | Wynika z Φ-EOM mass scale (m_Φ env-dependent) | FORCED | Determinowane przez którego rzędu coupling kontroluje fazę |

**Implikacja:** falsyfikacja `β_ppE^TGP = −15/4` przez GWTC-3 (2026-05-09) NIE była
"falsyfikacją β jako parametru", tylko **falsyfikacją specific Taylor expansion** M9.1''
`(4-3ψ)/ψ` — czyli konkretnych wartości `{a_1, a_2, a_3, b_2, ξ_3}`. Cykl emergent-metric
recovery (Phase 3-4) pokazał, że **przestrzeń `g_eff[Φ]` zawiera zero-β region** — falsyfikacja
dotyczyła punktu w przestrzeni, nie samej przestrzeni.

To jest exactly to, co native-first methodology pozwala zobaczyć clean: w PPN-/ppE-only
języku falsyfikacja "β_ppE=−15/4" wygląda jak "TGP wykluczone"; w native języku to "punkt
{a_n^M911} wykluczony, sąsiedztwo otwarte".

---

## §3 — Operacyjne reguły (binding dla cykli 2026-05-10+)

### §3.1 — Three-layer presentation MANDATORY

Każdy cykl gravity-related od 2026-05-10 raportuje wyniki w **trzech warstwach**:

| Layer | Co zawiera | Status |
|---|---|---|
| **L1 — Native predictions** (primary) | Obserwable liczone bezpośrednio z `g_eff[Φ]` + Φ-EOM (deflekcja, Shapiro, perihelion, range residual, GW phase, GW amplitude, dispersion). Każda obserwabla z constraint na native Taylor coefs `{a_n, ξ_n, b_n, c_0, κ_σ}`. | **PRIMARY** |
| **L2 — PPN/ppE projection** (consistency map) | Tabela: dla każdego PPN/ppE parametru, *jaka kombinacja* native coefs się na niego mapuje. Format: `β = f(a_2, ξ_2, ξ_3, b_2, a_1²)`. | **CONSISTENCY CHECK** |
| **L3 — Falsification map** | Dla każdego observational bound (Cassini, LLR, GWTC, BBN), *który* native coef jest constrained, w jakim window. Falsyfikacja propaguje do struktury teorii (`{a_n}` constraints), NIE do "parameter X excluded". | **DIAGNOSTIC** |

**Anti-pattern (do avoidance):** prezentacja "TGP predicts β_PPN = 1 EXACT" jako rezultat
Phase 2 — to L2-only, brakuje L1 (jakie obserwable?) + L3 (jakie native coefs constrained?).
Poprawnie: "L1: deflekcja, Shapiro, perihelion shift TGP-formula = GR-formula at 1PN order;
L2: γ=β=1 EXACT; L3: Cassini constrains `a_1`, `b_1` z `b_1=−a_1`; LLR constrains `ξ_2,
a_2`."

### §3.2 — Forced-zero parameters: nie liczyć, deklarować

α₁, α₂, α₃, ζ_i: NIE wykonywać sympy verification "α_i = 0". To jest tożsamość z
Lorentz-invariance substratu, nie predykcja. Format: jedna linijka deklaracji + cite
S05 + §5.1 FOUNDATIONS.

### §3.3 — Native parameter count audit

W każdym closing Phase, sekcja `## §X — Native parameter audit`:

```
Independent Taylor coefs constrained by this cycle: {a_1, a_2, b_1, ξ_2}
Free coefs (deferred to other cycles): {a_3, ξ_3, c_0, κ_σ}
Forced coefs (from substrate symmetry): {ζ_i, α_i} ≡ 0
Total native parameter count for gravity sector: ~5-7 (vs PPN-language 10)
```

To honest accounting parameter freedom — przeciwdziała sztucznej inflacji liczby
swobodnych parametrów.

### §3.4 — Connection do TGP_NATIVE_COMPUTATIONAL_PATTERNS

Ten dokument jest **siostrzany** względem `TGP_NATIVE_COMPUTATIONAL_PATTERNS.md`:

- TGP_NATIVE_COMPUTATIONAL_PATTERNS: **anti-BD-drift** — formuły wyglądające jak BD ale
  TGP-meaning (form-meaning mapping)
- PPN_AS_PROJECTION (ten doc): **anti-projection-confusion** — parametry phenomenological
  framework'ów (PPN, ppE, …) wyglądające jak fundamental ale TGP-projection (chart-meaning
  mapping)

Oba dokumenty dotyczą tego samego głębszego problemu: **język standardowej fizyki ma chart
zaprojektowane dla metric theories i Lagrangian-fundamental theories. TGP, jako substrate
theory, używa tych chartów dla porównania z observational bounds, ale nie jest tymi
chartami.**

---

## §4 — Generalizacja na cosmology (FRW, perturbations)

PPN nie działa w cosmology (post-Newtonian expansion zakłada slowly varying U; FRW jest
homogenic). Tam analogiczne projekcje to:

- **w_eff (effective dark energy EoS)** — projekcja `g_eff[Φ̄(t)]` na FRW chart
- **f_growth, σ_8 deviations** — projekcja `δg_ij[δΦ]` na perturbed FRW
- **Σ, μ (modified gravity post-Friedmann)** — projekcje analogiczne do γ, β ale
  cosmological-scale

Wszystkie te są **induced** w sensie tego dokumentu — TGP daje je z Φ-EOM przy odpowiednim
choice of Φ̄ background. Native-first methodology w cosmology cycle: liczyć `H(z)`,
`f_growth(z, k)`, `δT/T(l)` bezpośrednio, mapować na `{w_eff, Σ, μ}` jako consistency map.

---

## §5 — Implikacje dla TGP_FOUNDATIONS §3.6.2

Obecna sekcja `TGP_FOUNDATIONS.md §3.6.2 PPN constraints jako derivation (NIE postulat)`
prezentuje:

```
γ_PPN = 1 ⟺ b_1 = -a_1
β_PPN = 1 ⟺ ξ_2 = ξ - a_2·ξ³/2
```

**Zgodne z native-first methodology**, ale prezentacja inwersyjna: *PPN parameter as
function of constraint on native coefs*. Native-first wersja: *native coefs constraint set
by observational bound on PPN projection*. Te są równoważne logicznie, ale różnią się
metodologicznie:

- Inverse (current): TGP musi *passować* PPN test, więc *wymuszamy* `b_1=−a_1`
- Native-first: native coefs `{a_1, b_1}` mają obserwowalny constraint z Cassini +
  Shapiro; observational projection na γ-basis jest 1.0 ± 2.3·10⁻⁵

Native-first jest **structurally cleaner**: TGP nie "passuje testu" jako post-hoc tuning;
TGP daje obserwable, observational bound mapuje się na native coef bound. Brak tuning
language — jest physics + measurement.

**Recommendation:** §3.6.2 dostać reframe w next FOUNDATIONS edit (multi-session work, nie
blocker), gdzie:
- Najpierw lista obserwabli z native derivation
- Potem mapowanie na PPN basis jako tabela
- Constraint na native coefs explicitly z observational bounds

---

## §6 — Status (2026-05-10) i next steps

**Status:** 🟢 ACTIVE — methodological note. Binding dla wszystkich cykli post-2026-05-10.

**Cycle integrations (already aligned):**
- `op-emergent-metric-from-interaction-2026-05-09` Phase 2-4: results ALREADY in native
  form (`b_1=−a_1`, `ξ_2=...`, β_ppE^new(c_0)). **Addendum 2026-05-10** dodany do
  cyklu jasno mówi że to native results, PPN/ppE są L2 projection.

**Pending integrations:**
- `TGP_FOUNDATIONS.md §3.6.2`: reframe inversion direction (per §5 above) — multi-session
- `audyt/T01_LIGO3G_falsifier`: aktualizacja jako native-coefs falsifier zamiast β_ppE
  parameter falsifier — pending audit cycle
- Future gravity-sector cycles: MANDATORY three-layer presentation (§3.1)

**Cross-link:**
- `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` — siostrzany doc (anti-BD-drift)
- `meta/CALIBRATION_PROTOCOL.md` — może wymagać dodania anti-pattern "PPN-only
  presentation without native layer" (pending iterative review)

---

**Cycle authored:** 2026-05-10 (post-conversation z autorem o γ vs β natywności).

**Insight credit:** autor (TGP), 2026-05-10. Doc formalizuje obserwację jako binding
methodology.
