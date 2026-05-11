---
title: "TGP Native Computational Patterns — anti-BD-drift binding protocol"
date: 2026-05-10
type: meta-protocol
status: 🟢 ACTIVE — All 7 patterns written (Patterns 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7); pending T1.B + T1.C cascade
binding_scope: "All research cycles touching gravity / inertia / momentum / mass / GW sektor (post-2026-05-10)"
supersedes: null
related:
  - "[[CALIBRATION_PROTOCOL.md]]"
  - "[[CYCLE_LIFECYCLE.md]]"
  - "[[../TGP_FOUNDATIONS.md]] §3.5 (dual-V), §5 (czym TGP NIE jest), §6 (Lenz back-reaction)"
  - "[[TGP_TECH_DEBT_REGISTRY.md]] (planned, not yet created — early-stage BD-bridges audit)"
parent: "[[README.md]]"
tags:
  - meta
  - anti-BD-drift
  - TGP-native
  - computational-patterns
  - binding-protocol
  - ASK-RULE
---

# TGP_NATIVE_COMPUTATIONAL_PATTERNS

## §0 — Po co ten dokument

### §0.1 — Zdiagnozowany problem

Po sesji burzy mózgów 2026-05-09/2026-05-10 (post-cykl `op-recovery-V-mPhi-parametric-analysis`)
zidentyfikowano **systematic BD-drift** w pracach agentów na frameworku TGP. Przy braku
TGP-native computational tools, agenci sięgają po standardowe narzędzia fizyki (Brans-Dicke,
Horndeski, scalar-tensor QFT) i wstawiają je w obliczenia TGP, co prowadzi do:

- Kategorialnych błędów interpretacyjnych (np. "δΦ-quantum jako graviton-substitute")
- Fałszywych downgrade verdiktów (np. "mechanism iii FAILS bo m_Φ ~ M_Pl" — interpretacja
  oparta o std-physics; w TGP m_Φ jest environment-dependent observable, NIE universal mass)
- Cascade propagacji błędnych framingów do downstream cykli
- Wzajemnie sprzecznych LOCKs między cyklami (BD-form formuły z TGP-meaning bez explicit annotation)

### §0.2 — Trzy diagnozy źródłowe (ze burzy 2026-05-10)

| Diagnoza | Co znaczy | Solution w tym dokumencie |
|---|---|---|
| **A** | Brakuje TGP-native matematyki (computational tools) | §2 Worked-example patterns |
| **B** | Form-meaning mismatch nieoznaczony w predecessor cycles | §4 Form-meaning mapping table |
| **C** | LLM training bias na standard physics (BD/Horndeski/QFT) | §1 ASK-RULE (mandatory pause-and-ask) + §3 Red flags + §5 Pre-flight checklist |

### §0.3 — Tech-debt warning

Framework TGP_v1 ma znaczny **dług technologiczny**: wczesne cykle (2026-Q1, pre-foundations)
budowały "mosty do standardowej fizyki" jako narzędzia selekcji w czasie gdy formalizmu
TGP-native jeszcze nie było. **Te mosty nie są ground truth** — wiele finalnie nie przetrwa
i zostanie zastąpione TGP-native derivations.

**Implikacja dla agenta:** inheriting LOCK z predecessor cycle ≠ "to jest TGP truth".
Predecessor LOCK może być:
- (a) TGP-native LIVE — można dziedziczyć bez restrykcji
- (b) BD-form / TGP-meaning — można dziedziczyć ALE z explicit annotation z §4
- (c) Tech-debt bridge — NIE można dziedziczyć bez ASK-RULE (§1 Trigger B)

Lista (c) cykli — `meta/TGP_TECH_DEBT_REGISTRY.md` (planowany, nie istnieje jeszcze; pending
dedicated audit cycle). **Do czasu utworzenia listy: każdy LOCK z formułą BD-form wymaga
ASK-RULE.**

### §0.4 — Binding scope

**Wszystkie nowe cykle z 2026-05-10+ MUSZĄ:**

1. **Phase 0 README:** zawierać explicit sekcję `## §X — TGP-native check (mandatory)` z
   pre-flight checklist (§5 ten dokument)
2. **Pre-Phase-1:** wykonać ASK-RULE (§1) gdy spełniony jakikolwiek z 4 trigerów
3. **Phase 1+ sympy:** explicit cite z §4 form-meaning mapping dla każdej formuły BD-form
4. **Phase FINAL:** spawn `bd-drift-audit` subagent (per `CALIBRATION_PROTOCOL.md` §4.4 — to
   be added) który read Phase 1+ outputs i flag potential drifts
5. **Cycle classification:** `STRUCTURAL DERIVED` ZAREZERWOWANE dla cykli z bd-drift-audit PASS;
   cykle bez audyt → automatically `STRUCTURAL_CONDITIONAL z BD-drift-audit-pending`

**Cykle pre-2026-05-10:** zachowane w current state (`folder_status` preserved); audit BD-drift
post-hoc w dedicated cycles (np. T2.A `op-mPhi-verification-fluid-analog-audit` jako pierwszy).

---

## §1 — MANDATORY ASK-RULE (binding for all agents)

### §1.1 — Cztery triggery STOP-and-ASK

Agent MUSI przerwać pracę i zapytać użytkownika gdy spełniony JAKIKOLWIEK z poniższych warunków:

#### Trigger A — No TGP analogy visible

**Definicja:** próbujesz policzyć obserwowalne X (Newton force, GW amplitude, observable
mass, decay rate, scattering cross-section, ...) i NIE WIDZISZ w:
- `TGP_FOUNDATIONS.md` §3-§6
- Ten dokument §2 patterns
- Worked-example cykle z §6 cross-references

explicit TGP-native sposobu jak to liczyć.

**Co masz NIE robić:** NIE sięgać domyślnie po standardową fizykę (BD/QFT/GR derivation)
"żeby coś policzyć". To jest known anti-pattern (§3) który prowadzi do BD-drift.

**Co masz robić:** **STOP and ASK USER** (template §1.2).

#### Trigger B — Predecessor analogy seems outdated / tech-debt

**Definicja:** widzisz w predecessor cycle LOCK / formułę / interpretation która:
- (i) Wygląda jak BD/Yukawa/scalar-tensor (BD-form)
- (ii) NIE ma explicit TGP-native annotation w §4 form-meaning mapping
- (iii) Mogła być "most do standardowej fizyki" z wczesnego etapu (pre-foundations)

**Co masz NIE robić:** NIE inheriting LOCK assuming że jest TGP-native LIVE. Tech-debt LOCKs
mogą być semantically incorrect mimo że algebraically valid.

**Co masz robić:** **STOP and ASK USER** czy LOCK jest:
- (a) TGP-native LIVE (proceed z inheritance)
- (b) BD-form / TGP-meaning (proceed z explicit annotation, dodaj do §4)
- (c) Tech-debt bridge (NIE inherit, znaleźć alternative TGP-native derivation)

#### Trigger C — Just reproducing standard physics without TGP mechanism

**Definicja:** wpadłeś w obliczenie które mechanicznie reproduktuje standardową fizykę
(np. Yukawa screening, BD ω parameter, scalar-tensor PPN deviation, Feynman scattering)
**bez explicit TGP mechanism** wyjaśniającego skąd te liczby się biorą W TGP STRUCTURE.

**Diagnoza:** jeśli Twoja derywacja nie używa NIGDZIE:
- Pełnego nieliniowego Φ-EOM (z `D_kin = (1/3φ²)∇²(φ³)` operator)
- T^μν momentum flux integral
- σ_ab gradient-strain composite
- g_eff[Φ] funkcjonału (NIE g_eff jako independent variable)
- Lenz back-reaction picture (Phase 5 emergent-metric / M9.2)
- Environment-dependent observable mass (Pattern 2.5)

...to prawdopodobnie robisz BD-translation, NIE TGP-native derivation.

**Co masz robić:** **STOP and ASK USER** czy:
- (a) TGP genuinely recovers std physics by design w tym regime (proceed z explicit annotation)
- (b) To jest BD-drift i wymaga TGP-native re-derivation

#### Trigger D — Inheriting LOCK suspect

**Definicja:** Twój cycle Phase 1+ MUSI inheritować LOCK z predecessor cycle. Sprawdź:
- Czy LOCK ma TGP-native interpretację (z §4 mapping)?
- Czy źródłowy cycle ma `bd-drift-audit` PASS?
- Czy formuła używa standard-physics terminology (Yukawa, BD, propagator, exchange)?

**Jeśli ANY answer jest "no" lub "unsure":** **STOP and ASK USER**.

### §1.2 — Template ASK message

Gdy któryś trigger spełniony, agent musi napisać do user'a wiadomość w formacie:

```
🚨 TGP-NATIVE ASK-RULE TRIGGERED

**Trigger:** [A / B / C / D]
**Cycle:** [name of current cycle]
**Phase:** [Phase 0/1/2/3/...]

**Context (1-2 zdania):**
[Co próbujesz policzyć, krótko]

**Standard-physics default I would use:**
[Konkretnie jaką BD/QFT/GR derivation byś użył, np. "Yukawa propagator (1/(4π·r))·exp(-m·r)
dla matter-matter scalar exchange"]

**What I don't see (TGP-native gap):**
[Konkretnie czego brak w foundations / patterns / cross-references, np. "Brak Pattern 2.X
dla [observable]; predecessor cycle Y uses formula Z bez annotation"]

**Question for user (1-3 opcje, jeśli możliwe):**
1. [opcja A]
2. [opcja B]
3. [opcja C — np. "spawn dedicated TGP-native derivation cycle"]

**Impact if I guessed wrong (BD-drift propagation):**
[np. "Cycle outcome interpretive-only; downstream cycles inherit BD-styled framing;
 amendment needed; framework cascade re-evaluation"]
```

### §1.3 — Konsekwencje pominięcia ASK-RULE

Jeśli agent NIE wykona ASK-RULE i wstawi BD-drift do cyklu:

| Stage | Co się stanie |
|---|---|
| Phase FINAL | `bd-drift-audit` subagent (mandatory post-2026-05-10 per `CALIBRATION_PROTOCOL §4.4`) wyłapie drift |
| Cycle verdict | Auto-downgrade: planned `STRUCTURAL DERIVED` → `STRUCTURAL_CONDITIONAL z BD-drift flag` |
| Downstream | Wymagany amendment cycle z TGP-native re-derivation (analog T3.4 amendment chain) |
| Net effect | **Więcej pracy niż gdyby zapytano up-front.** ASK-RULE oszczędza czas. |

**Konkretny precedens:** Phase 1 cyklu `op-recovery-V-mPhi-parametric-analysis-2026-05-09`
(38/38 sympy PASS) — algebraicznie poprawny, ale interpretive framing BD-styled. Wymaga
amendment + re-frame. Gdyby ASK-RULE istniał wcześniej i zostal wykonany pre-Phase-1,
zaoszczędziłby ~2 sesje pracy.

### §1.4 — Wyjątki (kiedy NIE pytać)

✅ **Pure mathematical sympy proof** (algebraic manipulation, symbol disjointness, identity verification)
   — np. moja Phase 1 §2.1-2.6 "free symbols of β_ppE^new disjoint od V params" jest
   pure algebra, nie wymaga ASK-RULE; ALE interpretive jump "więc V is free parameter"
   wymaga ASK-RULE.

✅ **Numerical simulation z explicit code** gdzie TGP-native equation już zdefiniowane
   w predecessor cycle (np. M9.2 BVP solver dla Φ-EOM z β-mass term).

✅ **Documentation / cycle setup / cross-reference work** (no new derivation, no interpretive
   claims).

✅ **Closing cycle z verdict CONDITIONAL/HALT** (no need to derive anything new — tylko
   verdict statement).

### §1.5 — Wyjątki NIE applies gdy

❌ **Interpretacja sympy proof → wnioski fizyczne.** Każdy interpretive jump z algebra do
   physics wymaga ASK-RULE check. Sympy nie kłamie; interpretacja może.

❌ **"Naturalne" mapping na obserwacje.** Jeśli widzisz formułę i myślisz "to wygląda jak
   Yukawa/BD/standard X" — to red flag, ASK.

❌ **Inherited LOCK z BD-form bez §4 annotation.** Bez annotation: ASK.

❌ **Cycle predecessor nie ma `bd-drift-audit` PASS** (czyli cykl pre-2026-05-10): treat
   suspect, ASK przy inheritance.

### §1.6 — Self-check (przed every Phase 1 sympy)

Agent przed odpaleniem Phase 1 sympy MUSI mentally run:

1. ❓ "Czy mój Phase 1 plan używa TGP-native Pattern (z §2)? Który? __"
2. ❓ "Jeśli używa BD/std-physics derivation: czy mam justification (§4 mapping lub explicit reason)?"
3. ❓ "Czy m_Φ używam jako universal stałą czy environment-dependent (Pattern 2.5)?"
4. ❓ "Czy interpretacja outcome zakłada Φ-quantum carrier picture? (anti-pattern z §3)"
5. ❓ "Czy inheritance LOCKs z predecessor wszystkie mają §4 entry lub `bd-drift-audit` PASS?"

**Jeśli ANY odpowiedź sugeruje BD-drift risk: WYKONAJ ASK-RULE.**

---

## §2 — TGP-native computational patterns

> **Status patterns:** 2.2 + 2.5 napisane (2026-05-10); reszta pending iterative review.

### Pattern 2.1: Static Φ_eq[ρ] from arbitrary source (foundation) ⭐ **NAPISANY**

#### 2.1.1 — Setup

**Problem:** zadana statyczna gęstość źródła `ρ(x)`. Pytanie: jaka jest equilibrium configuration
`Φ_eq(x)` pola Φ wokół tego źródła?

**To jest fundament wszystkich innych patterns.** Bez poprawnego Φ_eq[ρ] żadna analiza
Newton'a (Pattern 2.2), bezwładności (2.3), GW emission (2.4), pressure-dome (2.6) ani
environment-dependent observable mass (2.5) nie ma sensownego punktu wyjścia.

**Standard physics default (anti-pattern, §3 red flag #1):**

W BD/scalar-tensor: linearize Φ = Φ_0 + δΦ, dropuje nonlinear terms, rozwiązuje Klein-Gordon
Helmholtz equation:

```
(∇² - m²) δΦ = -q·ρ/Φ_0          ← LINEARIZED ANTI-PATTERN
δΦ_eq(x) = -∫ G_Yukawa(x,x') · q·ρ(x')/Φ_0 d³x'
G_Yukawa(x,x') = exp(-m|x-x'|) / (4π|x-x'|)
```

To jest BD-mode framing. Działa **tylko w weak field daleko od source** gdzie nieliniowość
istotnie znikoma.

**TGP-native podejście (Pattern 2.1):**

Pełne nieliniowe równanie pola TGP (foundations §3 eq 134):

```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ_0 - γΦ³/Φ_0² = -qΦ_0·ρ              [STATIC FORM]

w canonical postaci:  D_kin[Φ] + V'(Φ) = -qΦ_0·ρ
gdzie  D_kin[Φ] = (1/3φ²)·∇²(φ³)         [foundations §3 eq 137]
       φ = Φ/Φ_0
```

**KLUCZOWE:** `D_kin` jest NIE Klein-Gordon `∇²`. Zawiera nieliniowość `2(∇Φ)²/Φ` która
fundamentalnie zmienia asymptotykę solution. **NIE upraszczać do ∇² bez explicit weak-field
justification.**

#### 2.1.2 — Method (step-by-step)

**Step 1:** Identify regime — czy weak field, czy strong field?

| Regime | Kryterium | Method |
|---|---|---|
| Asymptotically weak | r → ∞, φ → 1 (i.e. Φ → Φ_0), Δφ ≪ 1 | Linearization OK; ALE explicit annotate jako weak-field approximation, NIE general truth |
| Strong field | Bliskie source (φ daleko od 1), φ ~ 2/3 lub ~ 4/3 (M9.1'' multi-vacuum) | Pełna nieliniowość MUSI być zachowana |
| Mixed | Source z gradientem, gdzie różne regiony mają różne φ | BVP solver z proper boundary conditions; np. M9.2 scipy.solve_bvp template |

**Step 2:** Wybór reprezentacji — `Φ`, `φ = Φ/Φ_0`, lub `ψ = Φ/Φ_0` (M9.1'' notation).

W M9.1'' canonical: `ψ ∈ [0, 4/3]` z multi-vacuum (`ψ=0` trivial unstable, `ψ=2/3` cosmological
stable, `ψ=4/3` BH horyzont degenerate). Pełna nieliniowa solution może wędrować po tej
range w zależności od source strength i geometry.

**Step 3:** Setup BVP (boundary value problem) lub initial value problem:

- **BC at infinity:** `Φ(r → ∞) → Φ_0_cosmological` (asymptotic vacuum)
- **BC at origin** (lub source location): regular, no singularity
- **BC at finite extent source:** `Φ` continuous + first derivative continuous across surface

Dla spherically symmetric source, redukcja do ODE w `r`. Reference template: M9.2
`m9_2_momentum.py` z `scipy.integrate.solve_bvp` na `v(r) = r·ε(r)`.

**Step 4:** Solve nieliniową ODE/PDE numerically (analitycznie tylko w specjalnych przypadkach
typu pure M9.1'' Schwarzschild-like).

**Step 5 (KRYTYCZNE):** Verify asymptotic behavior:

| Asymptotyka | Interpretacja | Implication |
|---|---|---|
| `Φ - Φ_0 ~ A/r` (1/r tail) | M9.1''-Schwarzschild-like; Newton range INFINITE; m_Φ irrelevant for static | TGP-native; `m_Φ_intrinsic` heavy NIE blokuje Newton@AU |
| `Φ - Φ_0 ~ A·exp(-m·r)/r` (Yukawa tail) | Weak-field linearization regime; m_Φ controls range | BD-form result; mark §4 mapping F2 PARTIAL |
| Other (M9.1'' multi-vacuum, kink/defekt) | TGP-specific, dependent na full ψ-EOM nonlinearity | Document case-by-case |

**Step 6:** Compute observables (Pattern 2.2 momentum flux, 2.3 inertia, 2.5 local m_Φ_observable).

#### 2.1.3 — Anti-pattern warnings (specific to Pattern 2.1)

🚨 **Warning A:** NIE linearize before identifying regime. "Domyślnie linearize" jest
BD-drift habit. W M9.1'' (i likely jego recovery V analogues) nieliniowość D_kin daje
fundamentalnie inną asymptotykę.

🚨 **Warning B:** NIE używaj Helmholtz Green function exp(-m·r)/(4π·r) jako "exact static
field". To jest weak-field linearization result. W TGP-native Pattern 2.1 jest *jednym z
możliwych asymptotic forms*, NIE general truth.

🚨 **Warning C:** Multi-vacuum awareness — w V_M9.1'' są 3 critical points (ψ=0, 2/3, 4/3).
Solution może asymptotycznie reach inny vacuum niż "obvious" cosmological ψ=2/3 w pewnych
geometries. Sprawdź który vacuum solution naprawdę reach.

🚨 **Warning D:** Newton range w TGP NIE jest bezpośrednio `1/m_Φ` (jak w BD). Określa go
asymptotyka pełnego nieliniowego solution. Heavy m_Φ_intrinsic ≠ short Newton range.

#### 2.1.4 — Worked-example references

- **op-newton-momentum / M9.2:** spherical source z β-mass, BVP solver z Yukawa decay BC.
  ⚠️ ALE: użyto `β = 0.01` (relatively small mass²); dla heavy m_Φ trzeba sprawdzić czy
  asymptotyka się nie zmienia.
- **op-emergent-metric Phase 1-3:** definicja g_eff = G[{Φ_i}] funkcjonału (foundations §3.6.1
  refined ansatz {A, B, C}). Daje strukturalny framework dla Φ_eq solution wraz z g_eff
  dependence.
- **Open work:** explicit Pattern 2.1 derivation dla M9.1'' V_M9.1''(ψ) = -γ·ψ²·(4-3ψ)²/12
  z heavy γ ~ M_Pl² — **czy asymptotyka jest 1/r (Schwarzschild) czy exp(-m·r)/r (Yukawa)?**
  Spawned cycle candidate: `op-V-M911-static-asymptotics-audit` (z burza §13).

#### 2.1.5 — Open work

- 🔲 Explicit Pattern 2.1 dla M9.1''-recovery V z heavy m_Φ_intrinsic — sprawdź asymptotykę
- 🔲 Two-source BVP solver template (potrzebne dla Pattern 2.2 two-body Newton)
- 🔲 Multi-vacuum boundary handling (ψ=2/3 vs ψ=4/3 transitions w geometriach z BH)
- 🔲 Connection do Pattern 2.5 — local m_Φ_observable wyłania się z `V''(Φ_eq(x))` z Pattern 2.1

---

### Pattern 2.2: Newton force from momentum-flux integral ⭐ **NAPISANY**

#### 2.2.1 — Setup (TGP-native, NIE BD)

**Problem:** dwa źródła `M_1, M_2` w odległości `r`. Pytanie: jaka siła grawitacyjna `F^i_{1←2}`
działa na M_1 od M_2?

**Standard physics default (anti-pattern, §3 red flag #1):**

W BD/scalar-tensor: matter sources ρ_i couplają się bezpośrednio do scalar field Φ przez
vertex `q·ρ·Φ/Φ_0`. Dwie cząstki testowe wymieniają Φ-quanta (analog Yukawa exchange).
Force computed via tree-level Φ exchange propagator:

```
F_{BD-Yukawa} = -G_eff·M_1·M_2/r² · exp(-m·r)·(1 + m·r)        ← ANTI-PATTERN
```

To jest fundamentalnie BD-mode framing (matter ↔ Φ ↔ matter exchange channel).

**TGP-native podejście (Pattern 2.2):**

W TGP (foundations §5.1): nie ma δΦ-quantum jako exchange particle, NIE ma graviton-substitute.
Matter NIE sprzęga się bezpośrednio do Φ-quanta przez vertex; matter porusza się po geodezykach
`g_eff[Φ_total]`, gdzie `Φ_total = Φ_eq[ρ_1 + ρ_2]` jest **pełnym nieliniowym rozwiązaniem**
Φ-EOM z combined source.

**Klucz:** "siła Newton'a" jest manifestacją nierównowagi pola Φ (foundations §6 Lenz-podobny
obraz). Computowana przez **integral momentum flux** stress-energy tensora T^μν[Φ_total]
przez powierzchnię otaczającą source 1.

#### 2.2.2 — Method (step-by-step)

**Step 1:** Solve full nonlinear Φ-EOM dla combined source:

```
D_kin[Φ_eq] + V'(Φ_eq) = -q · (ρ_1 + ρ_2)            [foundations eq 134]
D_kin[Φ] = (1/3φ²)·∇²(φ³)         [canonical TGP form, NIE Klein-Gordon]
```

**NIE linearize prematurely.** Nieliniowość `D_kin` jest essential — w slabosym polu
linearization daje Klein-Gordon, ale w realistycznym Φ_eq z M9.1''-style metric
asymptotyka jest 1/r tail (Schwarzschild), NIE Yukawa.

**Step 2:** Compute scalar field stress-energy tensor:

```
T^μν[Φ] = (∂^μΦ)(∂^νΦ) - g_eff^μν · L_TGP[Φ]
       = (∂^μΦ)(∂^νΦ) - g_eff^μν · [(1/2)·K(φ)·g_eff^αβ·∂_αφ·∂_βφ - V(φ)]
```

z `K(φ) = K_geo·φ⁴` (foundations §3 z α=2 selection).

**Step 3:** Choose closed surface S_1 around source 1 (sphere of radius R, gdzie R ≪ r distance).

**Step 4:** Compute force jako integral momentum flux:

```
F^i_{1←2} = -∮_{S_1} T^ij[Φ_eq(x)] · n_j dA              [Newton force z field momentum]
```

**Step 5:** Verify: w weak-field limit (`Φ ≈ Φ_0 + δΦ_1 + δΦ_2`, gdzie `δΦ_i ≪ Φ_0`):
- Pole δΦ_1 spada jak 1/r (M9.1''-Schwarzschild) — dominanty contribution z linear part
- Cross-term `(∂δΦ_1)·(∂δΦ_2)` w T^ij daje `~ δΦ_1·δΦ_2 / r²` momentum flux
- Integral over S_1 daje `F ~ G·M_1·M_2/r²` (Newton)

**Step 6:** Identify Newton constant z TGP framework (Phase 5 emergent-metric LOCK):

```
G_eff = q² / (4π·Φ_0²·K_1)                  [BD-form, TGP-meaning per §4 mapping]
```

To jest TGP-native expression for Newton G — NIE jest to "Yukawa exchange coupling", ale
**coefficient w T^μν momentum flux integral** dla weak-field limit.

#### 2.2.3 — Anti-pattern warnings (specific to Pattern 2.2)

🚨 **Warning A:** NIE używaj `δΦ propagator` jako "matter-matter exchange channel". W TGP
matter NIE sprzęga się bezpośrednio do Φ-quanta — matter follows geodezyks of g_eff[Φ_total].
Propagator picture jest BD-translation.

🚨 **Warning B:** NIE używaj `exp(-m·r)/r` jako "Yukawa screening of Newton". W weak-field
linearization to się pojawia, ALE w pełnym nieliniowym rozwiązaniu (M9.1''-style) asymptotyka
może być inna (1/r Schwarzschild). Sprawdź pełne rozwiązanie before claiming Yukawa.

🚨 **Warning C:** NIE używaj G_eff = q²/(4π·Φ_0²·K_1) jako "Newton constant from coupling
vertex". To jest formula Phase 5 emergent-metric **z momentum-flux derivation**; BD-interpretation
("vertex coupling") jest tech-debt remnant. Cite §4 mapping.

#### 2.2.4 — Worked-example references

- **M9.2 / op-newton-momentum:** wykonana derivacja `F_back^i = -∫ρ·(∂ε_eq/∂x^i) d³x` dla
  spherical source — proof-of-concept Pattern 2.2 dla single-source self-force.
- **op-emergent-metric Phase 5:** wykonana derivacja G_eff = q²/(4π·Φ_0²·K_1) z linearized
  Φ-EOM (BD-form ALE TGP-meaning per §4).

#### 2.2.5 — Open work (Pattern 2.2 not yet fully developed)

- 🔲 Two-source Pattern 2.2 explicit derivation w pełnej nieliniowości (M9.1''-style)
- 🔲 Verification że full nonlinear Φ-EOM daje 1/r tail (NIE Yukawa) dla ψ=2/3 vacuum
- 🔲 Integration z Pattern 2.5 (environment-dependent m_Φ — czy m_Φ wpływa na momentum flux
  via local V''(⟨Φ⟩))
- 🔲 Spawn dedicated cycle `op-Newton-from-momentum-flux-TGP` (na razie placeholder, conditional
  na T2.A audit outcome)

---

### Pattern 2.5: Environment-dependent observable m_Φ (fluid analog) ⭐ **NAPISANY**

#### 2.5.1 — Setup: dlaczego standard "m_Φ jako fixed parameter" jest BD-drift

W standardowej fizyce skalarnej (BD, scalar-tensor, Higgs): masa pola jest **intrinsic
parameter** określony jednoznacznie przez Lagrangian. Computed once at vacuum:

```
m_Φ² ≡ V''(Φ_vacuum)         ← STANDARD: jeden uniwersalny parametr
```

Następnie używany wszędzie: w Yukawa propagator, w Cassini PPN deviation, w GW dispersion,
w particle physics scattering — to samo `m_Φ`.

**To jest BD-drift gdy zastosowane do TGP.** TGP single-Φ Lagrangian (foundations §3):
- Φ jest **jedno fundamentalne pole substratu**
- Cała materia, przestrzeń, geometria — emergent z dynamiki Φ
- Particle masses **EMERGENT** z structure pola (kinki, defekty per warstwa 3c §4 foundations)
- "Mass" jako observable jest **interakcja** particle z otaczającym Φ-tłem, NIE intrinsic
  parameter Lagrangianu

#### 2.5.2 — TGP-native definicja (fluid analog)

**Postulat (foundations §3.5.6 DRAFT, pending T2.A audit):**

```
m_Φ_observable(x) = V''(⟨Φ⟩_local(x))         ← TGP: environment-dependent
```

gdzie `⟨Φ⟩_local(x)` jest **lokalnie uśrednioną wartością tła Φ** w obszarze x.

**Fluid analog (sformułowany przez user'a 2026-05-10):**

> "Pole generowane przez pojedyńczą cząstkę, jej gęstość jest 'stała', ale zaburzenia pola
> i globalna wartość masy własnej cząstki **zależy od otoczenia**, nie tylko od jej własnych
> parametrów. Na Marsie ważę mniej niż na ziemi."

> "Tutaj wchodzimy w coś jak fizyka cieczy, im mniejsza gęstość tym mniej energii potrzeba
> do ruchu cząstki."

**Interpretacja:** cząstka ma **własny wkład** do Φ-pola (jej "static field configuration"),
ale **observable mass** = interakcja particle z TOTAL Φ-tłem (własne + global average).
Total Φ zależy od environment (other sources, cosmological background).

#### 2.5.3 — Kluczowe rozróżnienie

| Kategoria | Symbol | Co to jest | Przykład |
|---|---|---|---|
| **Intrinsic gap** | `m_Φ_intrinsic ≡ V''(Φ_0_vacuum)` | Mass-gap dla fluctuations wokół DEEP cosmological vacuum | M9.1'' V_M9.1''(ψ=2/3) → m_Φ_intrinsic = (2/√3)·M_Pl |
| **Local observable** | `m_Φ_observable(x) = V''(⟨Φ⟩_local(x))` | Mass-gap dla fluctuations w LOKALNYM environment | Near Sun: m_Φ_observable wartość ZNACZNIE INNA niż intrinsic |
| **Particle observable mass** | `m_particle(x)` | Pełna observable mass particle's z env interaction | Mars vs Ziemia: różne `m_particle` dla tego samego "własnego" wkładu |

#### 2.5.4 — Implikacje dla pre-existing cykli (KRYTYCZNE)

**op-mPhi-level0-verification Phase 1 (24/24 PASS):**
- Computed: `m_ψ = (2/√3)·M_Pl ≈ 1.41·10²⁸ eV`
- **Co to JEST (TGP-native interpretation):** `m_Φ_intrinsic` przy DEEP cosmological vacuum
  ψ=2/3 (V_M9.1'' specific form, falsified ALE local interpretation OK).
- **Co to NIE JEST:** `m_Φ_observable` w środowisku LIGO source (binary BH coalescence).
- **Co to NIE JEST:** `m_Φ_observable` przy AU od Sun (Cassini test environment).
- **Wnioski mPhi-verification "mechanism iii FAILS" mogą być BD-drift:** zastosowano
  `m_Φ_intrinsic` (universal heavy) tam gdzie powinno być `m_Φ_observable(LIGO source)`
  (potentially much lighter w środowisku z dużym ρ).

**T2.A audit (planowany, pending) sprawdzi:** czy `m_Φ_observable(LIGO BBH source)` jest
naprawdę M_Pl, czy znacznie lżejsze w środowisku z high local Φ deviation z vacuum.

#### 2.5.5 — Connection do Vainshtein-style screening (Pattern 2.7 emerges)

Variable `m_Φ_observable` z fluid-analog **automatycznie daje Vainshtein-style screening**:

- Dense matter region (high local Φ_eq deviation z vacuum) → `m_Φ_observable(local)` shifted
- Could be much heavier than `m_Φ_intrinsic` w sparse regions
- Heavy local m_Φ → short Yukawa range w tym obszarze → 5th force suppressed
- **Screening EMERGES naturalnie z TGP framework, NIE postulated jak w chameleon/symmetron**

To jest **strukturalna przewaga TGP** nad standardowymi modified gravity models.

#### 2.5.6 — Anti-pattern warnings (specific to Pattern 2.5)

🚨 **Warning A:** NIE używaj computed `V''(Φ_vacuum)` jako "m_Φ for all calculations".
Sprawdź jaki environment jest relevant dla observable, użyj `V''(⟨Φ⟩_local)`.

🚨 **Warning B:** NIE inheritować "m_Φ ~ M_Pl" z mPhi-verification verdict bez kontekstu.
To jest `m_Φ_intrinsic`; downstream cycles often need `m_Φ_observable(local environment)`.

🚨 **Warning C:** NIE assumować że "m_Φ ≪ ℏω_LIGO" requirement applies do `m_Φ_intrinsic`.
Requirement applies do `m_Φ_observable(LIGO source environment)`. Te mogą być znacznie
różne wartości.

🚨 **Warning D (subtelność):** "Particle mass in TGP" jest **trzecia** kategoria — nie m_Φ_intrinsic
ani m_Φ_observable, tylko pełna observable mass particle's w danym environment. To jest
domena warstwa 3c kinki/defekty (foundations §4) — currently strukturalny szkic, nie pełen
formalism. NIE conflate trzy kategorie.

#### 2.5.7 — Worked-example references

- **TGP_FOUNDATIONS §3.5.6 DRAFT** (pending T2.A confirmation, w przygotowaniu) — formal
  postulate "m_Φ_observable jako environment-dependent"
- **Op-newton-momentum / M9.2** — implicitly uses environment-dependent m_field via Lenz
  back-reaction (m_field zależy od `ε_eq` distribution, czyli env)
- **T2.A audit (planowany)** — explicit verification dla mPhi-verification environments

#### 2.5.8 — Open work (Pattern 2.5 not yet fully developed)

- 🔲 Sympy verification: dla concrete V form (np. V_orig dual-V z foundations §3.5), compute
  V''(⟨Φ⟩_local) dla typowych environments (cosmological, near-Sun, near-binary)
- 🔲 Quantitative comparison m_Φ_observable values across environments — czy 12-orders-of-magnitude
  variation jest plausible (analog do "Mars vs Ziemia" weight scaling)
- 🔲 Cross-check z Phase 5 emergent-metric Lenz back-reaction — gdzie m_field implicit env-dependent
- 🔲 Spawn `op-mPhi-environment-dependence-foundations` cycle dla formal derivation z H_Γ
  substrate dynamics (multi-session, deferred do post-T2.A)

---

### Pattern 2.3: Inertia from Lenz back-reaction ⭐ **NAPISANY**

#### 2.3.1 — Setup (TGP-native foundational picture)

**Problem:** zadane source o masie M (z Pattern 2.1 Φ_eq[ρ] solution). Pytanie: jaka jest
**bezwładność** tego source — masa którą "czuje" przy próbie przyspieszenia?

**Standardowy obraz (anti-pattern):**

W standardowej fizyce: m_inertial jest postulowany jako parametr Lagrangianu (Newton I axiom),
NIEZALEŻNY od m_grav. Zasada równoważności (m_b = m_g) jest dodatkowym postulatem (Einstein
weak equivalence principle).

**TGP-native obraz (foundations §6, kluczowy):**

> "Spoczynek = lokalna równowaga pola Φ wokół źródła ρ. Przyspieszenie źródła łamie tę
> równowagę → pole reaguje zwrotnie → siła back-reakcji opiera się zmianie konfiguracji
> → manifestuje się jako bezwładność. Stały ruch = pole 'ślizgające się' wraz ze źródłem,
> brak back-reakcji, zachowanie pędu (analog reguły Lenza)."

**Kluczowa różnica od BD/std:**

- m_inertial NIE jest postulatem, jest **EMERGENT** z TGP field theory
- m_inertial = m_grav AUTOMATYCZNIE (jeden q determinuje obie) — WEP **wynika strukturalnie**
- Newton I (constant velocity → no force) wynika z translational invariance Φ_eq
- Newton II (F = -m·a) wynika z linear back-reaction wokół Galilean-translated equilibrium

#### 2.3.2 — Method (step-by-step)

**Step 1:** Solve Pattern 2.1 dla statycznego source — get `Φ_eq(x)` lub `ε_eq(x) = Φ_eq - Φ_0`.

**Step 2:** Constant velocity case (Newton I check):

```
Source moving z v = const:  ρ(x, t) = ρ_0(x - v·t)
Test: czy Φ(x, t) = Φ_eq(x - v·t) jest exact solution Φ-EOM?
```

**Verification w non-relativistic limit (v ≪ c):**

```
∂_t² Φ ~ (v/c)² · ∇² Φ ≈ 0   (negligible w v ≪ c)
```

⟹ Galilean-translated `Φ_eq(x - v·t)` JEST exact solution Φ-EOM dla v=const.
⟹ **Brak back-reaction force** — Newton I structurally satisfied.

**Step 3:** Acceleration case (Newton II + m_inertial computation):

```
Source accelerating:  X(t) z ẍ = a (linear order)
Linear order ansatz:  Φ(x, t) = Φ_eq(x - X(t)) + δΦ_BR(x, t; a) + O(a²)
```

`δΦ_BR` jest back-reaction perturbation z **acceleration** (NIE z velocity).

**Step 4:** Substitute do Φ-EOM, linearize w a:

```
(∂_t² - c²∇² + m_eff²) δΦ_BR = -ä·∂_X(δΦ_eq) + (other O(a) terms)
```

(Tu m_eff² = V''(Φ_eq) jest **environment-dependent local m_Φ** per Pattern 2.5; NIE
universal m_Φ_intrinsic.)

**Step 5:** Compute back-reaction force on source:

```
F_BR^i = -q·M · ∂_i δΦ_BR(at source position)
       = -m_inertial · a^i           [linear w a, structural Newton II]
```

**Step 6:** Compute m_inertial jako field momentum integral:

```
m_inertial = ∫ (∇ε_eq)² d³x          [field gradient squared]
           = E_static / c²             [self-energy of static field config / c²]
```

To jest **canonical TGP definicja inertii** (M9.2.2 reference).

#### 2.3.3 — Equivalence Principle automatic from S05

**Klucz:** TGP single-field axiom (S05) — jedno pole Φ, jeden charge q.

| Mass | Source | Value |
|---|---|---|
| m_grav | S_mat = -q·ρ·Φ/Φ_0 coupling | M (q·M = single coupling) |
| m_inertial | E_static/c² + bare renormalization | M (after renormalization) |

**Oba derivują z TYCH SAMYCH parametrów** (q, ρ_source, K_1, Φ_0).

S05 lock: ratio m_inertial / m_grav = 1 EXACTLY, **NIE postulat, NIE fit**.

```
m_inertial / m_grav = 1   ← AUTOMATIC consequence of single-field S05
```

Weak Equivalence Principle (m_b = m_g for same source) jest **fundamentalną konsekwencją
TGP S05**, nie aksjomatem dodatkowym.

#### 2.3.4 — Anti-pattern warnings

🚨 **Warning A:** NIE postulate m_inertial niezależnie od m_grav. To jest std-physics axiom,
NIE TGP. W TGP m_inertial **derivuje** z field theory back-reaction.

🚨 **Warning B:** NIE używaj m_eff = m_Φ_intrinsic w back-reaction equation. Local m_Φ
zależy od Φ_eq w obszarze wokół source (Pattern 2.5). W realistycznych geometriach
m_eff_local może być znacznie inne niż universal V''(Φ_0_cosmological).

🚨 **Warning C:** NIE conflate m_inertial (bezwładność) z m_Φ (mass-gap fluctuations).
To są DWIE RÓŻNE OBSERWABLE:
- **m_inertial** = response of source ρ to external perturbation, ~ q²·M²
- **m_Φ_observable** = mass-gap dla fluctuations δΦ wokół local ⟨Φ⟩

#### 2.3.5 — Worked-example references

- **`op-newton-momentum/M9_2_results.md` §M9.2.2:** explicit derivation
  `m_field = ∫(∇ε_eq)² d³x = 3.483·10⁻²` dla Gaussian source z β-mass term. **Foundational
  template** — wszystkie inne Pattern 2.3 derivations powinny follow this structure.
- **`op-emergent-metric Phase 5`:** structural derivation m_inertial = m_grav AUTOMATIC z S05.
  Cykl 10/10 PASS verifies WEP automatic property.
- **`op-newton-momentum/B9_wep_microscope_composition_results.md`:** two-component WEP test
  (geometric, coupling, inhomogeneous ρ); η_TGP_lab = 1.32·10⁻²⁶ ≪ MICROSCOPE 1.1·10⁻¹⁵.

#### 2.3.6 — Open work

- 🔲 Pattern 2.3 dla strong-field source (M ≥ 0.5 M_BH-natural) — gdzie nieliniowość D_kin
  istotnie modyfikuje formula
- 🔲 Connection z Pattern 2.5: m_inertial computed at LOCAL m_Φ_observable, NIE intrinsic
- 🔲 Two-source m_inertial — czy bezwładność jednego ciała zmienia się w obecności drugiego
  (analog Mach principle: distant masses affecting local inertia)

---

### Pattern 2.4: GW emission from time-varying T^μν ⭐ **NAPISANY**

#### 2.4.1 — Setup (TGP-native vs std physics)

**Problem:** binary system (np. dwie BH/NS o czasowo zmiennej konfiguracji `ρ(x,t)`).
Pytanie: jaka jest amplituda GW signal wyemitowanego z tego systemu, w tym jaki jest
TT-component obserwowany przez LIGO?

**Standardowy obraz (anti-pattern, §3 red flag #4 + #9):**

W BD/scalar-tensor: matter sources ρ tworzą δΦ-quanta jako exchange particle. δΦ propaguje
z mass m_Φ (Goldstone-like dla light m_Φ). LIGO observable jest amplituda wave packet δΦ
z effective propagator (1/r)·exp(-m_Φ·r) far field. Dla m_Φ ≪ ℏω_LIGO: light scalar mode
"carries" GW signal.

To jest fundamentalna BD-mode picture. **Konflikt z foundations §5.1:**
> "Nie ma grawitonu — ani fundamentalnego, ani kompozytowego. Nie ma propagującej cząstki spin-2."

**TGP-native obraz (Pattern 2.4):**

GW jest **kolektywnym wzorcem** czasowo zmieniających się Φ-fluktuacji wokół dynamic
equilibrium `Φ_eq[ρ(t)]`, propagowanym przez `g_eff[Φ]` non-trivial structure.

NIE ma "Φ-quantum carrier"; jest **time-varying T^μν[Φ_total]** stress-energy field, którego
outgoing momentum flux at far field IS gravitational wave (analog momentum flux Pattern 2.2,
ale dla time-varying source).

**Kluczowa rola σ_ab gradient-strain composite (foundations §3.6.1):**

```
g_eff^ij = δ^ij·B(ψ) + σ^ij·C(ψ)/(Φ_0²·c²)
σ^ij = (∂^iΦ)(∂^jΦ) - (1/3)δ^ij·(∇Φ)²            ← TRACE-LESS gradient strain
```

`σ_ab` jest **TRACE-LESS** po konstrukcji → automatycznie ma TT-component przy odpowiedniej
projekcji. **TT modes h_+, h_× emergują z σ-coupling, NIE z δΦ-as-graviton.**

To rozwiązuje strukturalny "scalar-only at linearized" problem (foundations §3.6.10.4):
- Linearized δg_eff^ij = δ^ij·b_1·δΦ jest isotropic spatial scalar → TT-projekcja IDENTYCZNIE = 0
- Ale σ-channel (full nonlinear, z `(∂Φ)²` composite) daje nie-zerowy TT — w sposób TGP-native

#### 2.4.2 — Method (step-by-step)

**Step 1:** Quasi-static Φ_eq for time-varying source:

```
Source:  ρ(x, t) (binary system z orbital motion)
         Quasi-statyczny limit: Φ_eq[ρ(t)] ≈ instantaneous Φ_eq snapshot dla każdego t
```

Quasi-static valid gdy orbital frequency ω_orb ≪ ω_Phi_relax (Φ relaxation time-scale,
related do local m_Φ_observable per Pattern 2.5).

**Step 2:** Compute time-varying T^μν[Φ_eq(t)]:

```
T^μν[Φ_eq(t)] = (∂^μΦ_eq)·(∂^νΦ_eq) - g_eff^μν · L_TGP[Φ_eq]
```

Czasowa zależność wynika z motion ρ(x,t) → Φ_eq depends on t.

**Step 3:** Identify σ_ab dynamic component:

```
σ^ij(x, t) = (∂^iΦ_eq(x, t))·(∂^jΦ_eq(x, t)) - (1/3)δ^ij·(∇Φ_eq)²
```

W binary system z time-varying quadrupole moment `Q^TT_ij(t)`:

```
σ^ij_far-field(x, t) ~ Q̈^TT_ij(t - r/c) / r        (analog quadrupole formula)
```

(Per σ-3PN Phase 2 result, z BD-form/TGP-meaning per §4 mapping F3.)

**Step 4:** Compute outgoing momentum flux at far-field surface S_∞:

```
dE/dt = ∮_{S_∞} T^0i[Φ_eq(t)] · n_i dA
```

To jest **energia GW radiated outward** — w TGP-native picture, **emergent from collective
σ_ab pattern**, NIE z "Φ-quanta z mass m_Φ propagating".

**Step 5:** Compute h_TT amplitude observed by LIGO at distance D:

```
h_TT^ij(D, t) = (TT projection of g_eff^ij perturbation)
              = (C(ψ)/(Φ_0²·c²)) · σ^ij(x_observer, t - D/c)        [σ-channel]
              ~ (c_0 · ξ_eff/(8π·c⁴·D)) · Q̈^TT_ij(t - D/c)            [Phase 2 LOCK]
```

z `c_0 = 4π`, `ξ_eff = 4·G·Φ_0²` (post-T3.4-amendment LOCKs):

```
h_TT^σ = h_TT^GR EXACTLY at leading PN                     [§4 mapping F3]
```

**Step 6:** Verify role of m_Φ — KLUCZOWE rozróżnienie BD vs TGP:

| Picture | m_Φ rola |
|---|---|
| **BD anti-pattern** | m_Φ ≪ ℏω_LIGO required dla "carrier propagation" |
| **TGP Pattern 2.4** | m_Φ_observable(local) controls dispersion of collective wave packet w obszarze propagacji; dla LIGO: needs `1/m_Φ_observable > Gpc` (path length); ALE TO JEST environment-dependent (Pattern 2.5) — m_Φ_observable może być MUCH lighter w intergalactic medium niż w deep cosmological vacuum |

**Implication:** "mechanism iii FAILS bo m_Φ_intrinsic ~ M_Pl" verdict z mPhi-verification
może być BD-drift artifact. Real question: czy `m_Φ_observable(LIGO propagation environment)`
jest light enough — to T2.A audit scope.

#### 2.4.3 — Honest caveats (Yukawa-resolution audit)

Per `op-sigma-yukawa-audit Phase 1` verdict (35/35 PASS), σ-3PN Phase 2 + T3.4 amendment
calculations użyły **massless retarded Green function explicitly** (formal m → 0 limit).

**To jest acceptable w Pattern 2.4 picture jeśli:**
- m_Φ_observable(LIGO propagation environment) is light enough → m → 0 limit valid
- LUB: σ_ab gradient-strain composite provides effective mediation niezależnie od m_Φ

**Open verification scope (Pattern 2.4 incomplete):**
- 🔲 Explicit derivation outgoing momentum flux z time-varying Φ_eq[ρ(t)] z heavy m_Φ_intrinsic
  ALE light m_Φ_observable(intergalactic)
- 🔲 Compare quasi-static derivation z dynamical (full retarded) — gdzie quasi-static breaks
- 🔲 Verify że σ-channel TT-projection daje GR amplitude in TGP framework (per Phase 2 LOCK
  ale z explicit TGP-native annotation)

#### 2.4.4 — Anti-pattern warnings

🚨 **Warning A:** NIE używaj "δΦ propagator z mass m_Φ" jako derivation of GW emission.
W TGP-native, GW jest collective σ_ab pattern, NIE Φ-quantum exchange.

🚨 **Warning B:** NIE assume "m_Φ ≪ ℏω_LIGO" jako universal requirement bez spec which m_Φ
(intrinsic vs observable). Per Pattern 2.5: requirement applies do local m_Φ_observable
w propagation environment, NIE do m_Φ_intrinsic deep cosmological.

🚨 **Warning C:** NIE używaj light scalar Goldstone framing. TGP ma Z₂ discrete symmetry
(foundations §3.5), no continuous symmetry breaking, no Goldstone mode. σ-channel emerges
INACZEJ — przez gradient strain composite, NIE przez spontaneous symmetry breaking.

🚨 **Warning D:** Quasi-static limit MUST be justified. Jeśli source motion timescale
porównywalna do Φ relaxation timescale, quasi-static breaks i potrzebny full dynamic Φ-EOM
solution. (Otwarty problem dla compact binaries near merger.)

#### 2.4.5 — Worked-example references

- **`op-sigma-3PN-radiative-2026-05-09/Phase2_results.md` (24/24 PASS post-T3.4-amendment):**
  derivation `h_TT^σ = h_TT^GR EXACTLY at leading PN order`. ⚠️ z BD-form / TGP-meaning
  annotation per §4 mapping F3.
- **`op-T34-normalization-amendment-2026-05-09/Phase1_sympy.py` (17/17 PASS):**
  matching condition `c_0·ξ_eff = 16π·G·Φ_0²` clean re-derivation.
- **`op-sigma-yukawa-audit-2026-05-09/Phase1_results.md` (35/35 PASS):**
  Honest verdict — calculations są formal m → 0 limit, NIE direct LIGO observable bez
  Pattern 2.5 environment-dependent m_Φ verification.
- **`op-mPhi-level0-verification-2026-05-09/Phase1_results.md` (24/24 PASS):**
  Computed `m_ψ_intrinsic = 1.41·10²⁸ eV`. ⚠️ Verdict "mechanism iii FAILS" jest BD-drift
  IF m_Φ_observable ≠ m_Φ_intrinsic (T2.A audit scope).

#### 2.4.6 — Open work

- 🔲 Pattern 2.4 explicit derivation z Pattern 2.5 environment-dependent m_Φ
- 🔲 Verification że σ-channel TT-projection w nieliniowej strukturze daje GR amplitude
  bez relying na "Φ-quantum carrier" picture
- 🔲 Quasi-static vs dynamic comparison dla compact binaries
- 🔲 Connection z `op-recovery-V-mPhi` cycle (re-frame post T2.A audit)

---

### Pattern 2.6: Pressure-dome interaction (geometric) ⭐ **NAPISANY**

#### 2.6.1 — User's foundational intuition (burza 2026-05-10)

> "Załóżmy, że mamy 2 źródła pola Φ, tworzące wspólną przestrzeń. Im bliżej siebie się
> znajdują (uwzględniając to, że metryka sama jest zależna od Φ) tym więcej przestrzeni
> powstaje poza obiektami a mniej między nimi. To tworzy naturalny gradient ciśnienia
> obserwowalny jako przyspieszenie grawitacyjne."

> "Czyli przyspieszenie / pęd wynika z nierównowagi pola. Patrz na równania pędu — one
> powinny być corem dla tego."

To jest **najgłębsza TGP-native intuicja** dla grawitacji — fundamentalnie różna od:
- BD/scalar-tensor (matter-matter exchange via δΦ)
- GR (geodezyks of independent g_μν)

W TGP grawitacja jest **geometric pressure-gradient** wynikający z nieliniowej zależności
volume `√(-g_eff[Φ])` od combined Φ_total. Ekwivalentny do Pattern 2.2 momentum-flux
**w wyniku liczbowym** ale z explicit GEOMETRIC interpretation.

#### 2.6.2 — Setup (geometric framing)

**Problem:** N źródeł `ρ_1, ρ_2, ..., ρ_N` w pewnej przestrzennej konfiguracji. Pytanie:
jak rozkład "ilości przestrzeni" (volume element) zależy od konfiguracji, i jak to daje
acceleration na każde source?

**Standardowe podejście (anti-pattern):**

W BD/Newton: total potential U(x) = -Σ_i G·M_i / |x - x_i| (linear superposition).
Acceleration na test particle = -∇U. Brak geometric interpretation; pole jest "abstrakcyjny
potencjał".

**TGP-native (Pattern 2.6):**

Combined Φ_eq[ρ_total] (z Pattern 2.1, **NIELINIOWE** rozwiązanie) determinuje volume distribution:

```
dV(x) = √(-g_eff[Φ_eq(x)]) · d³x
```

Per foundations §3 eq przed §3.5: `√(-g_eff) = c_0 · φ` (M9.1'' specific) — czyli volume
element jest **bezpośrednio** funkcją local φ = Φ/Φ_0.

**Kluczowa obserwacja (user's geometric picture):**

W obszarze między dwoma sources φ jest "wciągnięte" w ich kierunku (oba sources mają
attractive coupling) → φ_between ma niższą wartość niż gdyby był jeden source → mniej
volume between sources.

W obszarze poza dwoma sources φ jest niezakłócone (lub mniej zakłócone w tej kierunce) →
więcej volume outside.

**To jest obraz geometric pressure** — i może być explicitly computed:

```
"Pressure" P(x) ≡ d(volume density) / d(some Φ-config parameter)
                ~ ∂(√(-g_eff))/∂Φ · ∂Φ/∂x   (in some sense)
```

Acceleration na source 1 wynika z gradient tej "pressure":

```
a^i_{source 1} ∝ -∇^i (volume density) at source 1 location
              = -∇^i (√(-g_eff[Φ_eq(x)]))    at x = position(source 1)
```

#### 2.6.3 — Ekvivalencja z Pattern 2.2 (rygorystyczne uzgodnienie)

**Claim (do verification w future cycle):** Pattern 2.6 (geometric pressure gradient) i
Pattern 2.2 (T^μν momentum flux) **DAJĄ TĘ SAMĄ NUMERYCZNĄ wartość** dla acceleration
of source.

**Heurystyczny argument:**
- T^μν w TGP zawiera kinetic part (∂Φ)(∂Φ) AND potential part V(Φ)
- Volume element √(-g_eff) zależy od Φ
- Energy-momentum conservation `∇_μ T^μν = 0` w obecności g_eff[Φ] łączy te dwa
- Momentum flux integral (Pattern 2.2) i volume gradient (Pattern 2.6) są **dwie strony
  tego samego coin** — jak hydrostatic pressure i momentum flux w fluid mechanics

**Rygorystyczna verification: open work (cycle candidate `op-pattern-26-pattern-22-equivalence`).**

#### 2.6.4 — Fluid analog (rozszerzenie Pattern 2.5)

Per user (burza 2026-05-10):

> "Tutaj wchodzimy w coś jak fizyka cieczy, im mniejsza gęstość tym mniej energii potrzeba
> do ruchu cząstki, nie mam tutaj lepszego analogu."

Pattern 2.6 + Pattern 2.5 razem dają picture analog do **incompressible fluid z lokalną
zmianą gęstości**:

| Fluid concept | TGP analog |
|---|---|
| Fluid density ρ_fluid(x) | `√(-g_eff[Φ_eq(x)])` (volume density) |
| Fluid pressure P(x) | "Geometric pressure" — gradient √(-g_eff) (Pattern 2.6) |
| Particle in fluid (effective mass) | Particle z observable mass `m_particle(x)` zależnym od local Φ-density (Pattern 2.5) |
| Buoyancy / pressure-driven force | Gravitational acceleration |
| Density variation | Multi-source Φ_eq[ρ_total] geometry |

**KLUCZOWE: TO JEST ANALOGIA, NIE FORMALNA REDUKCJA.** Φ-EOM ma fundamentalnie inne struktury
niż Navier-Stokes. ALE intuicyjnie pomocna.

#### 2.6.5 — Connection do "Mars vs Ziemia" weight

User's example (Pattern 2.5 §2.5.2): "Na Marsie ważę mniej niż na ziemi."

**Standardowa interpretation:** acceleration g_planet = G·M_planet/r²; Mars ma mniejszą
M, więc mniejsze g.

**Pattern 2.6 + 2.5 interpretation:**
- Earth jest source → tworzy local Φ_eq z determined volume distribution
- Background ⟨Φ⟩ near Earth surface ≠ ⟨Φ⟩ near Mars surface (Pattern 2.5)
- m_observable człowieka zależy od local ⟨Φ⟩ (Pattern 2.5)
- "Weight" = m_observable · acceleration_geometric_pressure
- **OBA effects** (mass change AND acceleration) wynikają z **jednej** local Φ-structure

To może implikować że "Mars vs Ziemia" weight differ by **slightly different ratio** niż
samo standardowe `g_Mars/g_Earth = 3.71/9.81 ≈ 0.378` — jeśli m_observable też się zmienia
(o tinych O((Φ_planet - Φ_cosmological)/Φ_cosmological)² corrections).

**Open prediction:** environment-dependent observable mass DOES give corrections do "weight"
formula at very high precision. Może być testable z atomic clock comparisons (Earth orbit
vs deep space) — connection z `op-Phi0-spatial-variation-predictions-2026-05-09` cycle
(6/6 PASS, atomic clocks + EP predictions logged).

#### 2.6.6 — Anti-pattern warnings

🚨 **Warning A:** NIE używaj linear superposition `U_total = Σ U_i` jako derivation TGP gravity.
W Pattern 2.6 combined Φ_eq[ρ_total] jest **NIELINIOWE** — inny niż sum of Φ_eq[ρ_i].

🚨 **Warning B:** NIE conflate "volume density" z energy density. W TGP `√(-g_eff[Φ])` jest
geometric volume element, NIE energy density. Choć they're related przez T^μν, są DWIE
RÓŻNE rzeczy.

🚨 **Warning C:** Fluid analog (§2.6.4) jest INTUICYJNY pomocniczy, NIE formal redukcja.
NIE dosłownie używać Navier-Stokes equations dla TGP.

🚨 **Warning D:** "Pressure gradient" jest interpretive language — actual computational tool
to T^μν momentum flux (Pattern 2.2). Pattern 2.6 dostarcza intuicji, NIE alternative computational
formalism.

#### 2.6.7 — Worked-example references

- **None yet** — Pattern 2.6 jest formalization user's intuition, dotychczas nie
  explicitly derived w żadnym cyklu.
- **Closest existing:** `op-emergent-metric Phase 6` (H6.1 structural unification) —
  "TGP ma JEDNĄ ZASADĘ generowania tensor structure z interactions"; Pattern 2.6 jest
  narrative interpretation tej zasady dla case Newton'a.
- **Cross-reference (open):** Verlinde 2010 entropic gravity, Padmanabhan emergent gravity,
  BEC analog gravity — analogous geometric pressure pictures w innych frameworkach (NIE
  identyczne, ale spokrewnione intuicją).

#### 2.6.8 — Open work

- 🔲 Explicit derivation Pattern 2.6 = Pattern 2.2 ekwivalencji (cycle candidate)
- 🔲 Two-source explicit volume distribution computation dla M9.1'' z `√(-g_eff) = c_0·φ`
- 🔲 Connection z Pattern 2.5 dla "Mars vs Ziemia" precise prediction
- 🔲 Spawn cycle `op-pressure-dome-newton-derivation-TGP` (multi-session)

---

### Pattern 2.7: Vainshtein-style screening as natural emergence ⭐ **NAPISANY**

#### 2.7.1 — Problem standardowych modified gravity

W standardowych modyfikacjach grawitacji (Brans-Dicke, scalar-tensor, Horndeski) zachodzi
**konflikt observational vs cosmological**:

- **Cassini bound:** \|γ_PPN − 1\| ≤ 2.3·10⁻⁵ → wymusza `ω_BD > 4·10⁴` (silne sprzężenie)
- **Cosmological dark energy phenomenology:** wymaga light scalar (m ~ H_0 ~ 10⁻³³ eV)

Te dwa są NIEKOMPATYBILNE w naive scalar-tensor — light scalar generuje 5th force violating
Cassini. Standard solution: **postulate Vainshtein-style screening** (lub chameleon, symmetron)
— mechanizmy AD-HOC dodawane do Lagrangianu by ratować phenomenology:

- **Chameleon:** m_φ(ρ) zależy od local matter density; w dense regions (Solar System) heavy → no 5th force; w sparse regions (cosmological) light → DE phenomenology
- **Symmetron:** Z₂-symmetric potential z spontaneous breaking dependent na ρ
- **Vainshtein:** kinetic term nonlinearity blokuje 5th force w high-curvature regions

**Wszystkie trzy są POSTULATAMI** — explicitly added do Lagrangian żeby pasowało do obserwacji.
Brak głębszego "why".

#### 2.7.2 — TGP framework: screening EMERGES naturally

**Klucz: Pattern 2.5 (environment-dependent observable m_Φ) jest already enough.**

W TGP `m_Φ_observable(x) = V''(⟨Φ⟩_local(x))`. Variable ⟨Φ⟩_local automatycznie daje:

- **Dense matter region:** local Φ_eq mocno odchylone od cosmological vacuum → V''(⟨Φ⟩_local)
  może być znacznie inne (potentially heavier) niż V''(Φ_0_cosmological)
- **Sparse region (cosmological background):** ⟨Φ⟩_local ≈ Φ_0_cosmological → m_Φ_observable
  ≈ m_Φ_intrinsic (potentially light)

**To IS Vainshtein-style screening** — dense regions have effective shift do heavier m_Φ,
suppressing long-range scalar effects, podczas gdy cosmological regions retain light m_Φ.

**Different od chameleon/symmetron:**
- Chameleon: m_φ(ρ) jest postulated funkcją w Lagrangianie
- TGP Pattern 2.7: m_Φ(x) emerges z **substrate dynamics** (V'' evaluated at local Φ_eq solution
  determined by ρ via Pattern 2.1)

**Tactical advantage:** TGP nie potrzebuje DODATKOWYCH parametrów dla screening. Wszystko
wynika z **jednej** struktury V(Φ) + nieliniowej Φ-EOM.

#### 2.7.3 — Method (analyzing screening dla danej geometrii)

**Step 1:** Identify environment regions:
- Dense (e.g., near Sun, near binary BH, lab Earth surface): high local |Φ_eq - Φ_0_cosmological|
- Sparse (e.g., intergalactic medium, deep cosmic vacuum): low local |Φ_eq - Φ_0_cosmological|

**Step 2:** Solve Pattern 2.1 dla każdego regionu — get local Φ_eq(x).

**Step 3:** Compute Pattern 2.5 — get m_Φ_observable(x) = V''(Φ_eq(x)).

**Step 4:** Identify screening pattern:

| Region | m_Φ_observable | Yukawa range 1/m_Φ | Effective screening |
|---|---|---|---|
| Dense (e.g., Solar System) | Potentially heavy (V''(Φ_eq_near_Sun)) | Short | 5th force suppressed; γ_PPN deviations small |
| Sparse (cosmological vacuum) | Potentially light (V''(Φ_0_cosmological)) | Long (Hubble-scale) | DE-like phenomenology emerges |

**Step 5 (KRYTYCZNE):** verify że screening pattern jest **automatic** z framework, NIE
requires fine-tuning V parameters.

#### 2.7.4 — Connection do mPhi-verification cycle (T2.A audit scope)

**op-mPhi-level0-verification verdict** ("m_ψ ~ M_Pl → mechanism iii FAILS"):

W BD-drift interpretation: m_ψ_universal = M_Pl, więc Yukawa range = Planck length, więc
GW Yukawa-suppressed everywhere. Mechanism iii fails uniwersalnie.

W TGP Pattern 2.7 interpretation:
- m_ψ_intrinsic = (2/√3)·M_Pl jest computed at deep cosmological vacuum ψ=2/3 → relevant
  dla DEEP COSMIC VACUUM environment
- W LIGO source environment (binary BH coalescence, high local Φ deviation z vacuum):
  m_Φ_observable może być znacznie INNE (lighter w sparse propagation paths, heavier w merger
  region)
- Mechanism iii może **STILL realize** jeśli effective m_Φ_observable along propagation
  path jest light enough

**T2.A audit scope:** quantitative computation `m_Φ_observable(LIGO propagation environment)`
dla V_M9.1'' (i recovery candidates) — verify czy fluid-analog interpretation odwraca
mPhi-verification verdict.

#### 2.7.5 — Predykcje dla testowalnych observables (z Pattern 2.7)

Pattern 2.7 implikuje (testable):

| Observable | TGP prediction | vs Standard |
|---|---|---|
| Atomic clock at Earth orbit vs Mars orbit | Slight Δω (z m_observable variation) | Standard GR: brak Δω, tylko gravitational time dilation |
| Lab tests of WEP w high vs low Φ-density | Maybe slight η_WEP variation | Std: WEP exact w lab |
| GW signal z mergers w voids vs dense regions | Slightly different effective mass / amplitude | Std: identyczne |
| 5th force search w lab vs cosmological | Suppressed lab, allowed cosmological | Std: same suppression everywhere |

**Cross-reference:** `op-Phi0-spatial-variation-predictions-2026-05-09/` (6/6 PASS) —
atomic clocks + EP predictions już logged. Pattern 2.7 dostarcza framework dla tych predictions.

#### 2.7.6 — Anti-pattern warnings

🚨 **Warning A:** NIE postulate Vainshtein/chameleon/symmetron jako fix dla TGP. Per Pattern 2.7,
screening emerges naturally — dodawanie ad-hoc mechanism jest BD-drift.

🚨 **Warning B:** NIE assume "TGP ma single universal m_Φ" co potem trzeba "screened". Pattern
2.5 already zapewnia variable m_Φ_observable.

🚨 **Warning C:** NIE używaj chameleon literature formuł "as is" w TGP. Te formuły zawierają
postulated `m_φ(ρ)` funkcje; TGP m_Φ(x) wynika z **substrate dynamics**, NIE postulated function.

🚨 **Warning D:** Pattern 2.7 jest **PROPOSAL**, nie verified TGP feature yet. Pełna verification
wymaga T2.A audit + dedicated cycle (op-Pattern-27-screening-derivation-TGP). Treat jako
working hypothesis, NIE proven LOCK.

#### 2.7.7 — Worked-example references

- **None yet** — Pattern 2.7 jest derivative Pattern 2.5; pełna verification wymaga
  T2.A audit + post-T2.A cycle.
- **Cross-reference:** chameleon/symmetron literature (Khoury-Weltman 2003; Hinterbichler-Khoury
  2010) — analogous emergent screening proposals w innych frameworkach. TGP picture jest
  **strukturalnie głębszy** (emerging z substrate dynamics, NIE postulated potential
  modification).

#### 2.7.8 — Open work

- 🔲 T2.A audit verification że m_Φ_observable variation jest sufficient dla screening
- 🔲 Spawn cycle `op-Pattern-27-screening-derivation-TGP` (post-T2.A, multi-session)
- 🔲 Cross-validation z `op-Phi0-spatial-variation-predictions` (already PASS — verify
  consistency)
- 🔲 Quantitative comparison z chameleon predictions (czy TGP daje similar phenomenology
  ALE z fewer free parameters)
- 🔲 Test predictions §2.7.5 — lab atomic clocks (most achievable near-term)

---

---

## §3 — Anti-pattern red flags table

Lista std-physics red flags które agent musi detect i flag (trigger ASK-RULE §1):

| # | 🚨 Red Flag | Co to wygląda | TGP-native alternatywa | Detect w sympy / kodzie |
|---|---|---|---|---|
| 1 | Yukawa exchange Newton | "G_eff = q²/(4π·Φ²·K)·exp(-m·r)/r ; matter-matter via δΦ propagator" | Pattern 2.2 momentum-flux integral | search: `exp(-m*r)`, `Yukawa`, `propagator`, `vertex coupling` |
| 2 | Universal m_Φ | "m_Φ² = V''(Φ_0) jako stała w wszystkich obliczeniach" | Pattern 2.5 environment-dependent observable | search: `m_Phi` w globalnym scope; brak `_local` qualifier |
| 3 | BD γ_PPN deviation | "(m·r)² Yukawa correction to γ_PPN" | γ_PPN = 1 EXACT structural identity (b_1 = -a_1) | search: `gamma_PPN` deviation z formula z `m_Phi` |
| 4 | Φ-quantum carrier | "δΦ z mass m_Φ propagates GW signal" | Pattern 2.4 collective T^μν patterns; brak grawitonu (foundations §5.1) | search: `Phi propagator` + `GW emission` ; `Phi-quantum` |
| 5 | g_eff as independent | "δg_eff dynamic mode, h_μν EOM" | g_eff = G[{Φ_i}] funkcjonał (foundations §3.6) | search: `g_eff EOM`, `metric perturbation EOM` |
| 6 | Scalar-tensor language | "ω_BD parameter", "Horndeski class", "Brans-Dicke" | TGP single-Φ S05; foundations §5.1 explicit "NIE BD/Horndeski" | search literally tych terms |
| 7 | "GR limit" obscured | "Match GR by tuning" bez TGP mechanism | "GR jako numerical analog z TGP-native structure (Pattern 2.X)" | search: `GR match`, `GR limit` bez explicit `via Pattern X.Y` cite |
| 8 | δΦ exchange between matter | "Two test masses interact via Φ-exchange" | Pattern 2.2 — matter follows geodezyks of g_eff[Φ_total] | search: "exchange between", "matter-matter via Phi" |
| 9 | Goldstone-like m → 0 limit | "Light scalar mode dla GW signal" | Foundations §5.1: brak grawitonu/spin-2 propagating; σ_ab gradient strain | search: `Goldstone`, `massless scalar mode for GW` |
| 10 | Postulated Vainshtein | "Add Vainshtein screening fix dla Cassini" | Pattern 2.7 — emerges naturalnie z Pattern 2.5 | search: `Vainshtein`, `chameleon`, `symmetron` |
| 11 | m_inertial postulated | "m_inertial as independent parameter" | Pattern 2.3 Lenz back-reaction emergent | search: `m_inertial` as input parameter, NIE derived |
| 12 | Universal Newton's G | "G_N as fundamental constant" | G_eff jako observable z Pattern 2.2 momentum flux; environment-dependent? (open) | search: `G_N` w deep theory derivations |

---

## §4 — Form-meaning mapping (BD-form ≠ BD-meaning)

Tabela formuł które WYGLĄDAJĄ jak BD/std physics ale mają TGP-native interpretację. Każdy
inherited LOCK z BD-form MUSI mieć entry tutaj — inaczej Trigger D ASK-RULE.

| # | Formula | BD interpretation (anti) | TGP-native interpretation | Cycle source | Status |
|---|---|---|---|---|---|
| F1 | `G_eff = q²/(4π·Φ_0²·K_1)` | "Effective Newton coupling z scalar exchange vertex" | Coefficient w T^μν momentum-flux integral dla weak-field limit; **derivation Phase 5 z linearized Φ-EOM jest BD-form, ale meaning jest momentum-flux**. Pattern 2.2 worked example. | Phase 5 emergent-metric LOCK | ✅ TGP-meaning verified |
| F2 | `δΦ_eq(r) = -q·M·exp(-m_eff·r)/(4π·K_1·Φ_0·r)` | "Yukawa propagator z scalar mass m_eff" | Static field response z linearized weak-field; **w pełnym nieliniowym M9.1''-style może mieć inną asymptotykę (1/r Schwarzschild)**; m_eff parametrizes weak-field linearization tail, NIE intrinsic propagator pole | Phase 5 emergent-metric + M9.2 | ⚠️ **PARTIAL** — pełna nieliniowa derivation pending |
| F3 | `h_TT^σ = h_TT^GR EXACTLY` | "Scalar mode reproduces tensor mode by tuning" | σ_ab gradient-strain composite ma TT-projection nie-zerową; matching condition c_0·ξ_eff=16π·G·Φ_0² jest **formal m → 0 limit dla TGP collective wave packet** (σ-3PN Phase 2 + T3.4 amendment); Yukawa-audit Phase 1 explicitly flagged "massless retarded Green used explicitly" | T3.4 amendment + σ-3PN Phase 2 | ⚠️ **CONDITIONAL** — Yukawa concern open |
| F4 | `β_ppE^new = (45/16)·Δe_2 + (45/16)·c_0·κ_σ` | "ppE deviation parameter z scalar-tensor 2.5PN" | 2.5PN binary inspiral phase coefficient z g_eff Taylor coefs (a_1, a_2, a_3, b_2, ξ_3) + σ-coupling (c_0, κ_σ); Phase 4 emergent-metric LOCK | Phase 4 emergent-metric | ✅ TGP-meaning verified |
| F5 | `m_inertial = ∫(∇ε_eq)² d³x` | "Field mass z gradient energy" | Lenz back-reaction integral: bezwładność jako field momentum w response do source acceleration. Pattern 2.3. | M9.2.2, Phase 5 emergent-metric | ✅ TGP-native LIVE |
| F6 | `c_0 = 4π` | "Magic number coupling" | Geometric coefficient z OP-7 T3.4 chain (`ξ_eff = 4·G·Φ_0²` z LOCK matching condition `c_0·ξ_eff = 16π·G·Φ_0²`) | op-c0-derivation cycle | ✅ TGP-meaning verified |
| F7 | `κ_σ = 1/(3π)` | "Magic number κ" | 1/π factor z orbital phase averaging × 1/3 z σ-trace structure | op-kappa-sigma cycle | ✅ TGP-meaning verified |
| F8 | `c_0 · κ_σ = 4/3 EXACT` | "Numerical match by coincidence" | Clean π cancellation z dwóch niezależnych source — strong consistency, NIE post-hoc fit | joint LOCK | ✅ TGP-meaning verified |

**Rule:** dodawać per cykl gdy nowy LOCK używa BD-form ale TGP-meaning. Status:
- ✅ TGP-meaning verified — można inheritować bez restrykcji
- ⚠️ PARTIAL / CONDITIONAL — można inheritować z explicit caveat citation
- ❌ Tech-debt confirmed — NIE inheritować bez ASK-RULE Trigger B (use TGP_TECH_DEBT_REGISTRY)

---

## §5 — Pre-flight checklist (Phase 0 mandatory)

Każdy nowy cycle Phase 0 README MUSI zawierać sekcję:

```markdown
## §X — TGP-native check (mandatory, pre-Phase-1)

[ ] **Q1 (Pattern coverage):** Czy reviewuję §2 patterns relevant dla mojego problemu?
        Lista relevantnych patterns: __________________
[ ] **Q2 (Red flags):** Czy zidentyfikowałem §3 red flags w moim Phase 1 plan?
        Lista znalezionych red flags: __________________
        Jeśli brak: explicit explanation czemu plan jest red-flag-free: __________________
[ ] **Q3 (Inherited LOCKs §4 mapping):** Lista inherited LOCKs i ich §4 entry status:
        - LOCK 1: ____________  →  §4 entry F#: ____  status: ✅/⚠️/❌
        - LOCK 2: ____________  →  §4 entry F#: ____  status: ✅/⚠️/❌
        Jeśli ANY ❌ lub brak entry: WYKONAJ ASK-RULE Trigger B.
[ ] **Q4 (Standard-physics tools):** Czy mój Phase 1 plan używa standardowych tools 
        (Yukawa, BD ω, GR derivation, propagator, exchange)?
        Jeśli TAK: explicit justify czemu NOT BD-translation (link do §2 pattern lub §4 mapping):
        __________________
[ ] **Q5 (m_Φ usage):** Czy m_Φ używam jako universal stałą czy environment-dependent (Pattern 2.5)?
        Jeśli universal: justify or use Pattern 2.5 environment-dependent: __________________
[ ] **Q6 (GR limit framing):** Czy plan distinguishes "TGP gives same number as GR (via TGP mechanism)" 
        vs "TGP IS GR by translation"?
        Mechanism cite: __________________
[ ] **Q7 (ASK-RULE self-check):** Czy są niewyjaśnione gaps gdzie sięgam po std physics?
        Jeśli TAK: WYKONAJ §1 ASK-RULE before Phase 1 sympy.
        Triggers fired: __________________
        Asked user: yes / no  (jeśli no: explain czemu wyjątki §1.4 applies)
```

---

## §6 — Worked-example cross-references

| Pattern | Cycle / file | Status w cycle | Notes |
|---|---|---|---|
| 2.1 (static Φ_eq) | `op-newton-momentum/M9_2_results.md` | implicit, partial | BVP solver dla single source z β-mass term |
| 2.2 (Newton momentum flux) | `op-newton-momentum/M9_2.2` | partial — single-source self-force | Two-source nieliniowy: open work |
| 2.3 (Lenz back-reaction) | `op-newton-momentum/M9_2.2` + `op-emergent-metric Phase 5` | ✅ done | Two cycle confirmation |
| 2.4 (GW emission) | `op-sigma-3PN-radiative-2026-05-09 Phase 2` | ⚠️ partial — z BD-form/TGP-meaning annotation per §4 F3 | Yukawa-audit flag pending |
| 2.5 (env-dependent m_Φ) | None yet — first explicit derivation w T2.A audit | 🔲 open | Pattern THIS document jest pierwszą explicit formulation |
| 2.6 (pressure-dome) | None yet | 🔲 open | User's intuition — needs formal cycle |
| 2.7 (Vainshtein emergent) | None yet — derivative of 2.5 | 🔲 open | Conditional na 2.5 confirmation |

**Related documents:**

- `meta/TGP_TECH_DEBT_REGISTRY.md` — **planned, not yet created.** Lista cykli oznaczonych
  jako BD-bridges z early stages. Pending dedicated audit cycle.
- `TGP_FOUNDATIONS.md §3.5.6 DRAFT` — pending T2.A confirmation; formal postulate
  environment-dependent m_Φ (Pattern 2.5 source).
- `meta/CALIBRATION_PROTOCOL.md §4.4` — pending addition of `bd-drift-audit` subagent
  protocol (T1.C deliverable).
- `meta/CYCLE_LIFECYCLE.md` — pending addition of pre-flight checklist (§5 ten dokument)
  to Phase 0 README template.

---

## §7 — Update log

| Date | Change | Author |
|---|---|---|
| 2026-05-10 | Document created post-burza-mózgów op-recovery-V-mPhi cycle. Outline §§0-7 + §1 ASK-RULE complete + Patterns 2.5 (env-dependent m_Φ) + 2.2 (Newton momentum flux) written. Patterns 2.1, 2.3, 2.4, 2.6, 2.7 placeholders pending iterative review. §3 red flags table complete (12 entries). §4 form-meaning mapping initial (8 entries: F1-F8). §5 pre-flight checklist complete. §6 cross-references draft. | Claudian (initial draft) + user iterative review (planned) |
| 2026-05-10 (later) | All 7 patterns now fully written. Pattern 2.1 (static Φ_eq foundation, with regime identification + asymptotic verification + multi-vacuum awareness). Pattern 2.3 (Lenz back-reaction inertia, m_inertial = ∫(∇ε_eq)² + EP automatic from S05). Pattern 2.4 (GW emission z time-varying T^μν z σ_ab gradient strain composite, NIE Φ-quantum carrier; honest Yukawa-resolution caveats). Pattern 2.6 (geometric pressure-dome, user's foundational intuition formalized; ekvivalencja z Pattern 2.2 jako open work; fluid analog explicit). Pattern 2.7 (Vainshtein-style screening as natural emergence z Pattern 2.5; predictions for atomic clocks + WEP testing). Document promoted z DRAFT do ACTIVE binding. | Claudian (post user "działaj" approval 2026-05-10) |

---

## §8 — Critical review checklist (każdy update tego dokumentu)

Przed merge update to TGP_NATIVE_COMPUTATIONAL_PATTERNS.md, sprawdzić:

- [ ] Czy nowy pattern ma worked-example reference (§6)?
- [ ] Czy nowy pattern ma anti-pattern warnings (§3 red flags update)?
- [ ] Czy zmiana wpływa na §4 form-meaning mapping (any LOCK status change)?
- [ ] Czy ASK-RULE (§1) wymaga update (np. nowy trigger)?
- [ ] Czy pre-flight checklist (§5) wymaga update?
- [ ] Czy update affects istniejące cykle (cascade implications)?
- [ ] Czy `STATE.md` wymaga update (jeśli ten dokument introduces new active rules)?
