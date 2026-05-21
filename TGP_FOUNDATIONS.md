# TGP — fundamenty interpretacyjne

**Cel dokumentu:** kanoniczny, jednostronicowy zapis ontologicznych
i interpretacyjnych zasad TGP, do którego sięga każda sesja agenta
przed rozpoczęciem pracy obliczeniowej. Ten dokument **nie zastępuje**
formalizmu w `core/sek08*.tex`, **definiuje natomiast jego znaczenie**.

**Data:** 2026-04-25.
**Autor:** ustalenia z dyskusji z autorem TGP (Mateusz, environmenstefan@gmail.com).
**Status:** binding dla wszystkich przyszłych prac w `research/`.

---

## 1. Ontologia TGP — jedno zdanie

> **TGP postuluje jedno fundamentalne pole skalarne `Φ` z symetrią Z₂.
> Wszystko inne — przestrzeń, czas, materia, grawitacja, interakcje,
> efektywne stopnie swobody — jest emergentne z dynamiki tego jednego
> pola.**

To jest **fundament nieruchomy**. Każda propozycja modyfikacji:
- dodania drugiego pola fundamentalnego (skalar, wektor, tensor),
- zmiany symetrii (poza Z₂),
- rezygnacji ze skalarnego charakteru,

**narusza fundament programu TGP** i jest a priori odrzucona, niezależnie
od tego, jak elegancko rozwiązuje konkretny problem techniczny. Pivot
substratu **w obrębie skalarnego Z₂** (np. zmiana potencjału, kinetyki,
struktury bondu) jest dozwolony jako narzędzie inżynieryjne.

## 2. Hierarchia formalizmu (cztery poziomy)

Patrz `rem:hierarchia-sektorow` w `core/sek08_formalizm.tex` (lin. 15–47).

| Poziom | Obiekt | Zawartość | Status |
|---|---|---|---|
| **0** | Substrat dyskretny `Γ = (V, E)` | Hamilton `H_Γ` (GL-bond, v2 2026-04-24), symetria Z₂, coarse-graining; `Φ = ⟨ŝ²⟩`, **`σ_ab = K_ab − (1/3)δ_ab Tr(K)`, `K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩`** (gradient strain composite, OP-7 T2 2026-04-25) | (W) |
| **1** | Równanie pola Φ | Operator `D_kin`, samointerferencja `N`, akcja TGP, α=2, β=γ | (W) |
| **2** | Metryka efektywna `g_eff^μν[Φ, σ_ab]` | **Metryka hiperboliczna (M9.1'', 2026-04-25)** `g_tt = -c₀²(4-3ψ)/ψ`, emergencja w limicie GR, FRW, PPN (γ=β=1 exact). **⚠ STATUS 2026-05-09: specific (4-3ψ)/ψ form FALSIFIED-OBSERVATIONAL** at 5σ przez GWTC-3 (zob. § 3 CRITICAL UPDATE poniżej + [[research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]]) | (E) but specific form falsified |
| **3** | Materia i pola cechowania | Sprzężenie metryczne, Dirac, **U(1)_em DERIVED native (Abelian z S05 phase), SU(2)_L + SU(3)_c DECLARED LIMIT (non-Abelian; per [[meta/TGP_W_Z_THEORETICAL_LIMIT.md]] post-2026-05-18 audit)**, generacje | (P) — split per audit 2026-05-18: Abelian native ✅ / non-Abelian declared limit 🔴 |
| **3a** | Hadron composition rule (color singlet $N_q - N_{\bar q} \equiv 0 \pmod 3$) | Compact U(1) winding mechanism; **warunkowy na input ±1/3, ±2/3 charges** | (D) cycle 2026-05-16 hadron-topology A− conditional |
| **3b** | Quark mass spectrum | Universal Φ-kink formula tested **HALT-B 2026-05-16** (structural ceiling 2.68× vs required 80,000×) | (R) requires extension |
| **3c** | Non-Abelian gauge dynamics (W/Z + gluons + Yang-Mills + asymptotic freedom + confinement σ) | **DECLARED LIMIT** — TGP minimal axioms structurally cannot derive | (L) per [[meta/TGP_W_Z_THEORETICAL_LIMIT.md]] |

**Materia (poziom 3) sprzęga się z `Φ` wyłącznie przez `g_eff`**, nie
bezpośrednio przez `Φ` (aksjomat `ax:metric-coupling`,
`sek08_formalizm.tex` lin. 11115–11132).

## 3. Akcja zunifikowana

Z `sek08a_akcja_zunifikowana.tex` lin. 29–97:

```
S_TGP[Φ, ψ_m] = ∫ d⁴x √(-g_eff) [ L_field(Φ) + L_mat(Φ, ψ_m) ]
```

z:

- **L_field** = (1/2) K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ),
  - K(φ) = K_geo φ⁴ (z α=2 selection),
  - V(φ) = (β/3) φ³ - (γ/4) φ⁴, **β = γ** (warunek próżni),
    - **⚠ POST-2026-05-09 dual-V clarification:** ta formuła V_orig odnosi się
      teraz do **MATTER SECTOR V** (V_orig, używanego w Phase 5 Mach inertia,
      T-Λ ρ_vac, particle masses przez field expansion wokół Φ_0 vacuum).
      W **gravity sektorze** używa się V_M9.1''(ψ) = -γψ²(4-3ψ)²/12 (G.0 P21).
      Patrz §3.5 dual-V structure poniżej.
  - φ = Φ/Φ_0 (zmienna bezwymiarowa).

- **L_mat** = -(q/Φ_0) φ ρ
  (sprzężenie minimalne źródła ρ z Φ przez czynnik φ).

- **Metryka efektywna** (forma hiperboliczna, M9.1'' 2026-04-25, izotropowa diagonalna):
  ```
  ds² = -c_0² (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j,    ψ = Φ/Φ_0
  ```
  Daje **γ_PPN = β_PPN = 1 exact** w 1PN; współczynniki PPN wyższych rzędów
  c₂=-1, c₃=+5/3, c₄=-10/3 (Schwarzschild reproduced, M9.1'' P1).

  **Forma potęgowa** `g_tt = -c²/ψ` (poprzednie M9.1) **została sfalsyfikowana**
  2026-04-25 (M9.1 T3: β_PPN=4 vs obs 1, 3·10⁴σ).

  **Forma eksponencjalna** `g_tt = -c² e^(-2U)` (pre-pivot 2026-04-24, weak-field)
  jest **równoważna M9.1'' do O(U)** ale różni się przy 2PN+:
  M9.1'' przewiduje explicit `|Δg_tt| = (5/6)U³` deviation od GR (testowalne LIGO 3G).

  > ## ⚠ CRITICAL UPDATE 2026-05-09 — M9.1'' (4-3ψ)/ψ form OBSERVATIONALLY FALSIFIED
  >
  > **[[research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] Phase 1.5** wyprowadza
  > G_SPA = 48 sympy-exact (test-particle, **4-level verified**: sympy LOCK 5/5
  > + independent hand-calculation E_TGP(U³)/m = 49/48 + numerical sanity
  > at U=0.1 6-decimal match + alternative SPA derivation orthogonal route).
  > Phase 1's heuristic G_SPA ≈ 1 (cytat Sampson-Yunes-Cornish 2013 dla
  > "metric-only deviations") **NIE applies** do TGP M9.1'' — SYC 2013
  > derived dla SMALL-PERTURBATION regimes (BD z 1/ω_BD, dCS z ζ_dCS),
  > a TGP M9.1'' ma STRUCTURAL O(1) modification.
  >
  > **Konsekwencja:** β_ppE^TGP^(b=-1) = -15/4 ≈ -3.75 (η=1/4, test-p exact),
  > **factor 48× LARGER** niż Phase 1 estymowało (-5/64 ≈ -0.078).
  >
  > **GWTC-3 RE-RUN** ([[research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]):
  > z corrected β, **TGP M9.1'' specific (4-3ψ)/ψ ansatz RULED OUT at 5.02σ**
  > by GWTC-3 combined ~90 BBH posterior (BF_TGP/GR = 3.5·10⁻⁶, log10(BF) =
  > -5.45 → "OVERWHELMING GR preference"). δφ̂_4_TGP = -0.865 vs observed
  > +0.05 ± 0.18 (1σ).
  >
  > **Status (4-3ψ)/ψ specific form:**
  > - **SFALSYFIKOWANE OBSERWACYJNIE** (5σ, GWTC-3, TIGER framework analysis).
  > - PREDICTIONS_REGISTRY M911-P1: LIVE → **FALSIFIED-OBSERVATIONAL**.
  > - PREDICTIONS_REGISTRY M911-P2: ratios {-23/10, -38/23, +337/228}
  >   (Phase 1 heurystyka) **WRONG** — alternative SPA daje -11161/504 ≠ -23/10
  >   → **WITHDRAWN-NEEDS-REDERIVATION**.
  > - PREDICTIONS_REGISTRY M911-P3: 2/4 channels invalid → **PARTIAL-FALSIFIED**.
  >
  > **Ocena programu TGP:**
  > - Falsyfikacja jest **wąska**: dotyczy SPECYFICZNEJ formy hyperbolic
  >   f(ψ) = (4-3ψ)/ψ. Sam program TGP (Φ z Z₂ symetrią, emergencja, etc.)
  >   **nie jest sfalsyfikowany**.
  > - Path forward: **S07 audyt** (M9.1'' derivation z fundamentu) ma
  >   eksplorować alternative f(ψ) ansatze konsystentne z α=2 vacuum Φ-EOM
  >   + 1PN (β=γ=1) + GWTC-3 constraints (|δφ̂_4| ≲ 0.18 1σ → β_ppE ≲ 0.78
  >   absolute, |G_SPA·Δα_3| ≲ 0.78 / (3/(128·1/4)) = 8.32 → constraint
  >   na alternative ansatze).
  > - Methodological finding: **G_SPA = 48 dla structural-modification
  >   theories** jest publishable result na własną. SYC 2013 framework
  >   nie applies poza small-coupling regime.
  >
  > Dla operational consistency programu, ten paragraph nadrzędny nad
  > zdaniem "M9.1'' przewiduje (5/6)U³ deviation testowalne LIGO 3G"
  > (które jest historyczne; teraz znaczy: "M9.1'' przewidywało, zostało
  > przetestowane, sfalsyfikowane").

- **Element objętościowy:** √(-g_eff) = c_0 · φ (dokładnie).

- **Wariacja `δS_TGP/δΦ = 0`** daje równanie pola
  (`prop:field-eq-from-action`, `eq:field-eq-reproduced`):
  ```
  ∇²Φ + 2 (∇Φ)²/Φ + β Φ²/Φ_0 - γ Φ³/Φ_0² = -q Φ_0 ρ.
  ```
  Operator kinetyczny `D_kin[Φ] = ∇²Φ + 2(∇Φ)²/Φ = (1/3φ²) ∇²(φ³)`
  jest **kanoniczną postacią TGP**.

## 3.5 Dual-V structure: gravity vs matter sektora (2026-05-09)

> **Framework clarification ustalona 2026-05-09** w cyklu
> [[research/op-Phi-vacuum-scale-2026-05-09/]] (sympy 84/88 PASS) +
> 4 spawned audit cykle, finalised w [[research/op-dual-V-structure-clarification-2026-05-09/]].

TGP framework posiada **dwa potencjały** dla różnych sektorów dynamiki —
**NIE konflikt**, lecz **complementary structure** wynikająca z dual-aspect
single-Φ field:

| Potencjał | Sektor | Status | Source |
|---|---|---|---|
| **V_M9.1''(ψ) = -γψ²(4-3ψ)²/12** | **Gravitational** | G.0 P21 LOCKED 2026-05-02 (specific form FALSIFIED 2026-05-09 — pending S07 alternative f(ψ)) | derived z R3 ODE + M9.1'' metric + K=ψ⁴ constraint |
| **V_orig(Φ) = -βΦ³/(3Φ_0) + γΦ⁴/(4Φ_0²)** z β=γ | **Matter** | CANONICAL maintained (A4 marker realization 2026-05-09) | field theory expansion wokół vacuum Φ_eq=Φ_0 |

### 3.5.1 Geneza dual-V

G.0 closure 2026-05-02 ([[research/op-g0-r3-from-canonical-projection/]])
explicit zostawił matter sektor jako "separate verification" (Phase1_results.md
linia 266: *"A4 (matter coupling) — wymaga osobnego sprawdzenia (G.0 nie dotyka L_mat)"*).
Niniejsza A4 verifikacja zrealizowana 2026-05-09 w audit chain:
- [[research/op-V-canonical-consistency-audit-2026-05-09/]] (10/10 PASS)
- [[research/op-MAG-Phase5-V-reference-clarification-2026-05-09/]] (10/10 PASS)
- [[research/op-dual-V-structure-clarification-2026-05-09/]] (10/10 PASS)

### 3.5.2 Reguła użycia (binding post-2026-05-09)

Każdy nowy cykl używający V(Φ) MUSI explicit cytować sektor:
- **Cykle gravity-related** (M9.1'' metryka, R3 ODE, Newton limit, PPN, GW)
  → MUSZĄ używać **V_M9.1''**
- **Cykle matter-related** (Mach inertia, ρ_vac cosmological, particle masses,
  field fluctuations) → MUSZĄ explicit cytować **V_orig** + sektor

V_M9.1'' i V_orig **NIE są Taylor expansions** of each other (sympy verified
[[research/op-dual-V-structure-clarification-2026-05-09/Phase1_dual_V_sympy.py]]
T3 PASS) — są niezależne funkcjonalne formy.

### 3.5.3 Φ_0 jako EFT scale-dependent free parameter

**Kluczowa konstatacja 2026-05-09** (post audit chain):

> Φ_0 jest **EFT scale-dependent free parameter** (analogicznie do Higgs VEV
> w Standardowym Modelu), **NIE single first-principles value**.

Wynika z niedualnej możliwości skalowania: wszystkie scenariusze (Φ_0=H_0
cosmological, Φ_0=v_EW EW, Φ_0=M_Pl Planck) działają **perturbative** w Phase 5
Mach inertia po **erratum 2026-05-09** ([[research/op-Phase5-MAG-erratum-2026-05-09/]]).

| Domena | Preferred Phi_0 | Source |
|---|---|---|
| Cosmological | Φ_0 ~ H_0 ~ 10⁻³³ eV | T-Λ closure (ρ_vac match 1.020) |
| EW (jeśli δ.2 derives) | Φ_0 ~ v_EW ~ 246 GeV | jeśli future EWSB derivation |
| Phase 5 Mach inertia | freely choosable post-erratum | wszystkie scenariusze działają |

#### 3.5.3.1 Quantitative framework update (2026-05-10)

**Status update post-Cycles 1+3:** Original §3.5.3 declaration (Φ_0 jest EFT scale-dependent
free parameter) **QUANTITATIVELY SUBSTANTIATED** przez:
- [[research/op-gamma-RG-running-derivation-2026-05-10/]] (88/88 sympy PASS, **GF.B-STRUCTURAL** verdict)
- [[research/op-EFT-Phi0-multi-scale-2026-05-10/]] (10/10 sympy PASS, formal EFT framework)

**γ_eff(μ) one-loop running** (standard Coleman-Weinberg ϕ⁴, Peskin-Schroeder Ch. 12):
$$\gamma_{\text{eff}}(\mu) = \frac{\gamma_0}{1 - \frac{3\gamma_0}{16\pi^2}\ln(\mu/\mu_0)}$$

**Φ_0(μ) one-loop running** (z mass anomalous dimension γ_m² = γ/(16π²),
sign per minimal subtraction; magnitude robust):
$$\frac{\Phi_0(\mu)}{\Phi_0(\mu_0)} \approx 1 - \frac{2\gamma}{16\pi^2}\ln(\mu/\mu_0)$$

**Numerical multi-scale** (z γ(M_Pl) = 0.1, Φ_0(M_Pl) = 1.0 natural units):

| Scale | γ_eff(μ) | Φ_0(μ)/Φ_0(M_Pl) | γ·Φ_0² (relative) |
|---|---|---|---|
| M_Pl (1.22·10²⁸ eV) | 0.1000 | 1.000 | 1.000 |
| M_Z (9.12·10¹⁰ eV) | 0.0930 | 1.050 | 1.025 |
| ω_LIGO (4·10⁻¹³ eV) | 0.0850 | 1.118 | 1.063 |
| H_0 (1.44·10⁻³³ eV) | 0.0790 | 1.178 | 1.096 |

**KLUCZOWE:** across **61 orders of magnitude w μ**, mild log running only — γ varies
factor 0.79, Φ_0 varies factor 1.18. **NIE 10⁸² scale separation needed dla Branch D
quantitative substantiation.**

**T-Λ closure (cosmological anchor):** γ(μ_cosm)·Φ_0²(μ_cosm) = 12·ρ_vac.
Pod Branch A (γ ~ M_Pl²·g̃, Φ_0 ~ H_0): g̃ = 12·ρ_vac/(M_Pl²·H_0²) ≈ 0.98 — Λ-CDM
cosmological coincidence (NIE first-principles derivation; STRUCTURAL POSTULATE z
T-Λ closure consistency).

**Branch identification (post-Cycle-1 GF.B-STRUCTURAL):**
- ✅ Branch A re-asserted (single-scale γ z mild log running)
- ❌ Branch D quantitative substantiation FAILS (one-loop daje O(log) only, NIE multi-scale separation)
- ❌ Branch B (γ ~ ω_LIGO²) UNREACHABLE z one-loop ϕ⁴ flow
- Pattern 2.5 (§3.5.6) preserved jako BINDING-PRINCIPLE z PHYSICAL APPLICATION CONDITIONAL na extreme environments

**Honest open questions:**
1. **β=γ vacuum condition** — generic fine-tuning lub TGP-specific RG fixed-point (Cycle 1 §3.1)
2. **β-cubic coupling running** — deferred
3. **Non-perturbative corrections** to γ-flow — deferred

**OP-1 M2** (M-derivation U(φ) z H_Γ) status: **PARTIALLY RESOLVED** — β-function derivable
i RG flow analytical, ALE specific γ ~ M_Pl² value remains STRUCTURAL POSTULATE z T-Λ closure
consistency, NIE first-principles derived.

### 3.5.4 Phase 5 erratum (2026-05-09)

Phase 5 MAG Mach inertia derivation
([[research/op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]])
miało **internal inconsistency**: zakładało simultaneously β=γ (vacuum) i
β<<γ (m_C²~3γ approximation), które są wzajemnie sprzeczne.

**Correct identification** (sympy [[research/op-Phase5-MAG-erratum-2026-05-09/Phase1_erratum_sympy.py]] 5/5 PASS):
$$V''(\Phi_0)\big|_{\beta=\gamma} = \gamma \quad \Rightarrow \quad m_C^2 = \gamma$$

(NIE m_C² = γ/3). Z m_C = M_Pl (T-Λ consistency): wszystkie Φ_0 scenariusze
działają. **44-rzędowa hierarchia v_EW/H_0** była **artefaktem** internal
inconsistency, NIE physics.

### 3.5.5 Wpływ M9.1'' falsification (2026-05-09)

Inny agent (2026-05-09) zaktualizował M9.1'' specific (4-3ψ)/ψ form jako
**OBSERWACYJNIE SFALSYFIKOWANY at 5σ** przez GWTC-3
([[research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]).

**Konsekwencje dla dual-V:**
- **Gravity sector V_M9.1'' specific form:** ✅ recovery path identified
  ([[research/op-emergent-metric-from-interaction-2026-05-09/]] CLOSED 57/57 PASS)
- **Matter sector V_orig:** **NIE AFFECTED** — niezalezne od gravity specific form
- **Dual-V framework structure:** **CONFIRMED** przez sektor separation

Cross-cycle consistency check
([[research/op-Phi-vacuum-scale-2026-05-09/CONSISTENCY_REPORT_post_M911_falsification.md]],
sympy 7/7 PASS): essential matter-sector findings (Phase 5 erratum, EFT Phi_0,
sek08a annotation) survive falsification; multi-vacuum gravity specific values
(ψ=0, 2/3, 4/3) require post-S07 re-analysis but methodology preserved.

### 3.5.6 m_Φ jako environment-dependent observable

> **🟢 STATUS UPDATE 2026-05-10:** Pattern 2.5 **BINDING-PRINCIPLE-CONFIRMED-ALGEBRAIC**
> z **PHYSICAL APPLICATION CONDITIONAL** na extreme environments.
>
> **Verification chain:**
> - T2.A audit ([[research/op-mPhi-verification-fluid-analog-audit-2026-05-10/]]):
>   identified ψ_± = (6 ± 2√3)/9 inflection points (CONDITIONAL initial finding)
> - T3 Phase 1 ([[research/op-V-M911-psi-profile-near-degenerate-2026-05-10/Phase1_sympy.py]],
>   23/23 PASS): V''(ψ_±) = 0 EXACT z ψ_+ ≈ 1.052 inflection above vacuum ψ=2/3;
>   linearization "fixed m_Φ" valid TYLKO dla \|δψ\| ≪ 0.385 — **BINDING-PRINCIPLE-CONFIRMED-ALGEBRAIC**
> - T3 Phase 2 (14/14 PASS): numerical BVP confirmed M_critical=15.80 for δψ → ψ_+
> - T3 Phase 3 (13/13 PASS): δψ_LIGO ≈ 10⁻¹⁰⁴ pod Branch A (typical LIGO sources);
>   mechanism iii FAILS dla typical sources
> - Cycle 1 Phase 4 ([[research/op-gamma-RG-running-derivation-2026-05-10/Phase4_matching.md]]):
>   GF.B-STRUCTURAL verdict (Branch A re-asserted); Pattern 2.5 NEGATIVE for typical LIGO
>
> **Final status (post-Cycle-1 GF.B):**
> - **BINDING-PRINCIPLE: CONFIRMED-ALGEBRAIC** — m_Φ field-dep mechanism strukturalnie real
> - **PHYSICAL APPLICATION: CONDITIONAL** — active TYLKO w extreme environments (δψ ~ 0.3+)
> - **NEGATIVE for typical LIGO** (δψ_LIGO ≈ 10⁻¹⁰⁴ → V''(ψ_local) ≈ V''(ψ_vac))
> - **NEGATIVE for Solar System** (Cassini bounds preserved; γ_PPN=β_PPN=1 EXACT 1PN)
> - **POTENTIALLY ACTIVE** dla binary BH near-horizon (study deferred do recovery V cycle reactivation)
>
> **Pattern 2.5 jest TGP-NATIVE composite-field mechanism** (Φ=⟨ŝ²⟩ field-dep), DISTINCT
> z RG-scale running γ(μ) (§3.5.3.1). Combined picture:
> $$m_{\Phi,\text{observable}}^2(\psi, \mu) = V''(\psi_{\text{local}}) \cdot \gamma_{\text{RG}}(\mu)$$

**Framework clarification ustalona 2026-05-10** w sesji burzy mózgów post-cycle
[[research/op-recovery-V-mPhi-parametric-analysis-2026-05-09/Phase1_results.md]].
Identifikuje LIVE TGP feature dotychczas nie spisane explicit, prowadzące do systemic
BD-drift w cyklach używających `m_Φ`.

#### 3.5.6.1 — Fundamentalna konstatacja

**Standard physics framing (BD-drift anti-pattern):**

W standardowej fizyce (BD, scalar-tensor, Higgs, particle physics): masa pola jest
**intrinsic parameter** określony jednoznacznie przez Lagrangian, computed once at vacuum:

```
m_Φ² ≡ V''(Φ_vacuum)             ← STANDARD: jeden uniwersalny parametr
```

Następnie używany wszędzie: w Yukawa propagator, Cassini PPN deviation, GW dispersion,
particle scattering — to samo m_Φ.

**TGP-native framing:**

W TGP single-Φ Lagrangian: **Φ jest jedno fundamentalne pole substratu**, cała materia
i geometria emergentne. Particle "mass" jako observable jest **interakcja** particle z
otaczającym Φ-tłem — NIE intrinsic parameter Lagrangianu.

**Postulate (DRAFT):**

```
m_Φ_observable(x) = V''(⟨Φ⟩_local(x))            [environment-dependent observable]
```

gdzie `⟨Φ⟩_local(x)` jest **lokalnie uśrednioną wartością tła Φ** w obszarze x.

#### 3.5.6.2 — Fluid analog (sformułowany przez autora TGP)

> "Pole generowane przez pojedyńczą cząstkę, jej gęstość jest 'stała', ale zaburzenia pola
> i globalna wartość masy własnej cząstki **zależy od otoczenia**, nie tylko od jej
> własnych parametrów. Na Marsie ważę mniej niż na ziemi."

> "Tutaj wchodzimy w coś jak fizyka cieczy, im mniejsza gęstość tym mniej energii potrzeba
> do ruchu cząstki."

**Interpretacja:** cząstka ma **własny wkład** do Φ-pola (jej "static field configuration"),
ale **observable mass** = interakcja particle z TOTAL Φ-tłem (własne + global average).
Total Φ zależy od environment (other sources, cosmological background).

#### 3.5.6.3 — Trzy kategorie m_Φ (rozróżnienie BINDING)

| Kategoria | Symbol | Co to jest | Przykład |
|---|---|---|---|
| **Intrinsic gap** | `m_Φ_intrinsic ≡ V''(Φ_0_vacuum)` | Mass-gap dla fluctuations wokół DEEP cosmological vacuum | M9.1'' V''(ψ=2/3) → m_Φ_intrinsic = (2/√3)·M_Pl |
| **Local observable** | `m_Φ_observable(x) = V''(⟨Φ⟩_local(x))` | Mass-gap dla fluctuations w LOKALNYM environment | Near Sun / near LIGO source: m_Φ_observable ≠ m_Φ_intrinsic |
| **Particle observable mass** | `m_particle(x)` | Pełna observable mass particle's z env interaction | Mars vs Ziemia: różne `m_particle` |

**Reguła użycia (binding post-2026-05-10):**

Każdy nowy cycle używający m_Φ w obliczeniach MUSI explicit cytować KTÓRA z trzech
kategorii jest used:
- ✅ "Used m_Φ_intrinsic from V''(Φ_0_cosmological)" — applies tylko dla deep cosmological vacuum
- ✅ "Used m_Φ_observable(local environment X) from V''(⟨Φ⟩_local(X))" — applies dla
   observables w specific environment
- ⚠️ "Used m_Φ" bez qualifier — **AMBIGUOUS, wymaga ASK-RULE** per
   [[meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1 Trigger A

#### 3.5.6.4 — Implikacje dla pre-existing cykli

**[[research/op-mPhi-level0-verification-2026-05-09/Phase1_results.md]] verdict (24/24 PASS):**
- Computed: `m_ψ = (2/√3)·M_Pl ≈ 1.41·10²⁸ eV`
- **Co to JEST (TGP-native interpretation):** `m_Φ_intrinsic` przy DEEP cosmological vacuum
  ψ=2/3 (V_M9.1'' specific form, falsified ALE local interpretation OK)
- **Co to NIE JEST:** `m_Φ_observable` w środowisku LIGO source lub Cassini test
- **Verdict "mechanism iii FAILS"** może być BD-drift artifact: zastosowano `m_Φ_intrinsic`
  (universal heavy) tam gdzie powinno być `m_Φ_observable(LIGO source environment)`

**T2.A audit (planowany):** quantitative computation `m_Φ_observable(LIGO BBH source)` —
verify czy fluid-analog interpretation odwraca mPhi-verification verdict.

**[[research/op-recovery-V-mPhi-parametric-analysis-2026-05-09/Phase1_results.md]]:**
- 38/38 sympy PASS algebraic claims preserved
- Interpretive claims ("joint window m_Φ ≪ 6·10⁻²¹ eV Cassini-dominated") FLAGGED jako
  conditional pending T2.A — assumed universal m_Φ (BD-drift); actual constraint applies
  do environment-dependent m_Φ_observable

#### 3.5.6.5 — Vainshtein-style screening jako automatic emergence

Variable `m_Φ_observable` automatycznie daje **Vainshtein-style screening** (Pattern 2.7
w [[meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]]):

- Dense matter region (high local Φ_eq deviation z vacuum) → `m_Φ_observable(local)` shifted
  (potentially heavier) → short Yukawa range → 5th force suppressed w dense regions
- Sparse region (cosmological vacuum) → `m_Φ_observable` ≈ `m_Φ_intrinsic` (potentially
  lighter) → DE-like phenomenology

**Strukturalna przewaga TGP:** screening EMERGES naturalnie z framework, NIE postulated jak
w chameleon/symmetron/Vainshtein models. Wyjaśnia automatycznie:
- Cassini-compatible (Solar System: dense → heavy m_Φ_observable → no 5th force)
- DE-compatible (cosmological: sparse → light m_Φ_observable → long-range effects)
- Bez konfliktu observation vs cosmology charakterystycznego dla naive scalar-tensor

#### 3.5.6.6 — Status DRAFT — co T2.A audit musi rozstrzygnąć

T2.A audit ([[research/op-mPhi-verification-fluid-analog-audit-2026-05-10/]] planowany)
verifies czy postulate §3.5.6.1-§3.5.6.5 jest **quantitatively correct** w TGP framework:

**Claims do verification:**

| Claim | Test |
|---|---|
| C1 | V''(⟨Φ⟩_local) jest dobrze zdefiniowana (NIE pathological dla typical environments) |
| C2 | m_Φ_observable(LIGO BBH source) ≠ m_Φ_intrinsic — quantitative magnitude |
| C3 | m_Φ_observable(Earth lab) ≠ m_Φ_intrinsic — quantitative magnitude |
| C4 | Variation jest enough do explain Cassini compliance + DE phenomenology bez fine-tuning |
| C5 | Mechanism iii w `op-mPhi-verification` verdict "FAILS" jest BD-drift artifact (verdict reverses) |

**Possible outcomes:**
- (i) All C1-C5 PASS → §3.5.6 promoted to BINDING; mPhi-verification cascade DOWNGRADE-REVERSAL;
  recovery V cycle re-frame (see post-T2.A roadmap)
- (ii) Partial PASS → §3.5.6 BINDING z explicit caveats, downstream cascade z mixed implications
- (iii) C1-C5 FAIL → §3.5.6 WITHDRAWN; need alternative TGP-native explanation; foundations
  remain pre-2026-05-10 state z `m_Φ` ambiguity unresolved

#### 3.5.6.7 — Cross-references

- [[meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §2.5 (Pattern 2.5 explicit derivation),
  §2.7 (Pattern 2.7 Vainshtein-style screening), §1 (ASK-RULE Trigger A applies przy m_Φ)
- [[research/op-mPhi-level0-verification-2026-05-09/Phase1_results.md]] (verdict pending re-interpretation)
- [[research/op-recovery-V-mPhi-parametric-analysis-2026-05-09/Phase1_results.md]] (BD-drift detected)
- [[research/op-Phi0-spatial-variation-predictions-2026-05-09/]] (atomic clocks + EP predictions —
  testowalne consequences §3.5.6)
- [[research/op-Phi-vacuum-scale-2026-05-09/Phase_FINAL_close.md]] §3.5.3 (Φ_0 jako EFT
  scale-dependent free parameter — already established that Φ_0 nie jest universal; §3.5.6
  rozszerza to do m_Φ)

## 3.6 Emergent-metric framework: post-falsification recovery (2026-05-09)

> **Framework recovery established 2026-05-09** w cyklu
> [[research/op-emergent-metric-from-interaction-2026-05-09/]]
> (CLOSED, **STRUCTURAL DERIVED**, sympy 57/57 PASS, six requirements P1-P6: 6/6 RESOLVED).

Post-2026-05-09 GWTC-3 falsifikacji M9.1'' specific f(ψ) = (4−3ψ)/ψ formy,
TGP gravity sektor został odbudowany strukturalnie poprzez **interaction-driven
emergent metric framework**. Zamiast postulować lokalną funkcję f(ψ), g_eff
emerguje jako **funkcjonał konfiguracji wielu Φ-źródeł**.

### 3.6.1 Refined ansatz (3 niezależne funkcje)

```
g_eff^00 = -A(ψ)
g_eff^ij = δ^ij · B(ψ) + σ^ij · C(ψ) / (Φ_0² c²)
g_eff^0i = 0  (statyczny limit; gravitomagnetic 0i deferred)
```

gdzie σ^ij = (∂^iΦ)(∂^jΦ) − (1/3)δ^ij (∇Φ)² jest **gradient strain composite**
(level 0 obiekt z OP-7 T2 2026-04-25), aktywowany jako tensor source dla g_eff.

**Trzy niezależne funkcje** A(ψ), B(ψ), C(ψ) — w przeciwieństwie do
M9.1'' canonical, gdzie A·B = 1 było dodatkowym constraintem. M9.1'' jest
**konkretnym przypadkiem** w nowej rodzinie, NIE jedynym.

### 3.6.2 1PN/2PN constraints (native-first, PPN jako projekcja)

> **Methodology note (2026-05-10):** ta sekcja prezentuje wyniki w **native-first** form
> per [[meta/PPN_AS_PROJECTION.md]] binding methodology. PPN parametry {β, γ, ξ, α_i, ζ_i}
> są **projekcją** native Taylor coefs `g_eff[Φ]` na chart Willa, NIE fundamental
> obiektami TGP. Patrz `meta/PPN_AS_PROJECTION.md §1` dla klasyfikacji.

#### L1 — Native predictions (primary)

Phase 2 cyklu `op-emergent-metric-from-interaction` (sympy 7/7 PASS) wyprowadziło z
action variation, że obserwable 1PN/2PN regime liczą się bezpośrednio z geodezyk w
`g_eff[Φ_⊙]` jako funkcje native Taylor coefs:

| Observable | Native formula | Constrained coefs |
|---|---|---|
| Light deflection | θ ∝ (1+a_1+b_1)·GM/(bc²) | a_1, b_1 |
| Shapiro delay | Δt ∝ (1+a_1+b_1)·log(...) | a_1, b_1 |
| Perihelion shift | Δϕ = f(a_1, a_2, b_1, b_2, ξ_2, ξ_3) | a_2, b_2, ξ_2, ξ_3 |
| Nordtvedt (LLR) | η_N = g(a_1, a_2, b_1, ξ_2) | a_2, ξ_2 |

#### L2 — PPN projection (consistency map)

Mapowanie native → PPN basis:

```
γ_PPN − 1 = b_1 + a_1                          ⟹ γ = 1 ⟺ b_1 = −a_1
β_PPN − 1 = ξ_2 − ξ + a_2·ξ³/2                 ⟹ β = 1 ⟺ ξ_2 = ξ − a_2·ξ³/2
ξ_PPN (Whitehead) = 0                          ⟹ TGP preserves at 2PN
α_1 = α_2 = α_3 = 0 (preferred-frame breaking) ⟹ FORCED ≡ 0 z S05 + §5.1 Lorentz-invariance Φ-substratu
ζ_1 = ζ_2 = ζ_3 = ζ_4 = 0 (energy/momentum non-conservation) ⟹ FORCED ≡ 0 z covariant Φ-EOM
```

α_i, ζ_i są **strukturalnymi tożsamościami**, NIE predykcjami — wynikają z symetrii
substratu, nie wymagają sympy verification.

#### L3 — Falsification map (observational bounds → native coef constraints)

| Bound | Constrains native coefs | Window |
|---|---|---|
| Cassini (|γ−1| ≤ 2.3·10⁻⁵) | a_1, b_1 (przez `b_1 + a_1`) | preserved EXACT w obecnym ansatz |
| Mercury perihelion + LLR (|β−1| ≤ 8·10⁻⁵) | a_2, b_2, ξ_2 | preserved EXACT w obecnym ansatz |

**Solar system tests** są strukturalnie satysfakcjonowane przez `b_1 = −a_1` i `ξ_2 = ξ − a_2·ξ³/2`.
Nie jest to empiryczny fit specyficznej formy — to **constraint relacyjny** na native coefs,
który observational bounds mapują na 1.0 ± machine-precision-tail.

#### Native parameter freedom audit

```
1PN/2PN constraints: a_1+b_1 = 0, ξ_2 = ξ − a_2·ξ³/2     (2 native relations)
Free w 1PN/2PN regime: a_2, a_3, b_2, ξ_3, σ-coupling C(ψ)  (enters at 2.5PN)
Forced ≡ 0 from substrate symmetry: α_i, ζ_i               (4+4 = 8 PPN params, structurally zero)
```

**σ-coupling C(ψ) jest WOLNA w 1PN/2PN** — wnosi parametryczną swobodę dopiero przy
2.5PN binary inspiral (zob. §3.6.3).

### 3.6.3 β_ppE^new(c_0) parametric family (LOCK)

Phase 3 cyklu (sympy 5/5 PASS) wyprowadziło uogólniony SPA chain:

$$\beta_{\rm ppE}^{\rm new}\big|_{\eta=1/4} = \frac{45}{16} \cdot \Delta e_2 + \frac{45}{16} \cdot c_0 \cdot \kappa_\sigma(\eta)$$

gdzie:

$$\Delta e_2 = -a_1\xi_3 - 3 - \frac{4 a_2}{a_1^2} + \frac{4 b_2}{a_1^2} - \frac{8 a_3}{a_1^3} + \frac{16 a_2^2}{a_1^4}$$

**Recovery rodzina:** dla każdego (a_3, ξ_3) z **ξ_3 = (32 − a_3)/32** i c_0 = 0
otrzymujemy β_ppE = 0 EXACT (post-falsification recovery z 3PN parameter tuning).
Alternatywnie, z M9.1''-canonical params i **c_0·κ_σ = 4/3** otrzymujemy ten
sam zero-β rezultat (recovery z σ-coupling addition).

### 3.6.4 Phase 4 GWTC-3 compliance window (LOCK)

Phase 4 cyklu (sympy 8/8 PASS) zidentyfikowało:
- |β_ppE^new| ≤ 0.78 (1σ): **window EXISTS**
- Width on (a_3, ξ_3) parametric region: ~0.144 w ξ_3 space
- Width on c_0·κ_σ: [1.056, 1.611] (centered at 4/3)
- **c_GW = c structurally** (no Lorentz-violation, GW170817 constraint trivial)

### 3.6.5 Equivalence principle automatic z S05 (Phase 5)

Phase 5 cyklu (sympy 10/10 PASS) wyprowadziło:

$$\frac{m_{\rm inertial}}{m_{\rm grav}} = 1 \quad \text{AUTOMATYCZNIE z S05}$$

Ten sam q (single-field charge) determinuje OBIE masy:
- m_grav z S_mat coupling q·ρ·Φ/Φ_0
- m_inertial z back-reaction integral E_static/c² (analog Lenz)

**Newton I + II strukturalnie derived:**
- Newton I: Galilean-translated δΦ_eq = exact solution dla v=const (no back-reaction)
- Newton II: F_BR = -m_inertial · a (linear back-reaction structural form)

WEP (m_b = m_g) jest **automatyczną konsekwencją S05**, NIE postulatem.

### 3.6.6 H6.1 structural unification (Phase 6)

Phase 6 cyklu (sympy 11/11 PASS) potwierdziło **H6.1: structural unification**:

> TGP ma **JEDNĄ ZASADĘ** generowania tensor structure z interactions,
> applied at multiple levels:
> - **Level 2 (g_eff metryka):** functional G[{Φ_i}, σ_ab, Φ̄] z gradient cross-terms
> - **Level 3 (SU(2) spin):** dynamic-equilibrium soliton-Φ̄ interaction

Cross-consistency z [[research/op-SPIN-SU2-substrate-derivation-2026-05-08/]] (47/47 PASS):
- **Path A** (V_matter bifurcation): c_0-INDEPENDENT (dual-V lock) ✓
- **Path B** (M9.1'' horizon multipole): preserved by Phase 4 Path 2 ✓
- **Path C** (external embedding): c_0-INDEPENDENT (geometric) ✓

⟹ SU(2) emergence ROBUST regardless of Phase 4 family choice; ≥2 of 3 paths
preserved dla każdego punktu w GWTC-3 compliance window.

### 3.6.7 Phase 4 Path 2 strukturalnie preferowana

| Phase 4 path | Path A SU(2) | Path B SU(2) | Path C SU(2) | Sum |
|---|---|---|---|---|
| Path 1 (zmiana 3PN params) | ✅ | ❌ may break | ✅ | 2/3 OK |
| **Path 2 (keep params + add c_0)** | **✅** | **✅** | **✅** | **3/3 OK** |

**Canonical TGP recovery uses σ-coupling (c_0·κ_σ ≈ 4/3), NIE 3PN tuning.**
Argument strukturalny (preserve all 3 SU(2) paths), NIE empiryczny fit.

### 3.6.8 c_0 status — derivable, deferred multi-session

c_0 (leading σ-coupling coefficient C(ψ=1)) jest **framework-derivable**:
1. σ_ab coarse-graining z H_Γ substrate (~5-10 sesji)
2. Dynamic-equilibrium balance analog SPIN N16 (~3-5 sesji) ← preferred
3. SU(2) Path B exact preservation (~2-4 sesji)

Heuristic target: c_0·κ_σ ≈ 4/3 (z Phase 4 GWTC-3 zero-β analysis). Numerical
pinning HONESTLY DEFERRED do dedicated cycle `op-c0-derivation-from-substrate`.

### 3.6.9 Six requirements P1-P6 (6/6 RESOLVED)

| # | Requirement | Resolution |
|---|---|---|
| P1 | Formal definition g_eff = G[{Φ_i}] | ✅ Phase 1 (16/16) |
| P2 | 1PN reproduction γ=β=1 z derivation | ✅ Phase 2 (7/7) |
| P3 | 2.5PN β_ppE alternative do -15/4 | ✅ Phase 3 (parametric LOCK) |
| P4 | M9.2 Lenz back-reaction → m_inertial | ✅ Phase 5 (m_b=m_g AUTOMATIC) |
| P5 | Cross-consistency z 3 SU(2) paths | ✅ Phase 6 (H6.1 CONFIRMED) |
| P6 | Falsifiability w GWTC-3 | ✅ Phase 4 (compliance window) |

**Post-cycle deferred items RESOLVED (2026-05-09 follow-up cycles):**

### 3.6.10 Joint follow-up cycles closure (2026-05-09)

Trzy follow-up cykle uruchomione + zamknięte 2026-05-09:

#### 3.6.10.1 — c_0 derivation z OP-7 T3.4 chain
[[research/op-c0-derivation-from-substrate-2026-05-09/]] STRUCTURAL DERIVED
(heuristic numerical, 5/5 sympy):
```
c_0 ≈ 4π (geometric clean)  lub  c_0 ≈ 4π·1.06 ≈ 13.32 (z GW150914 ξ/G calibration)
```
Source: OP-7 T3.4 LOCK ξ_eff = 4π·G·Φ_0², Path A → Path B conversion.

#### 3.6.10.2 — κ_σ z 2-body PN orbital averaging
[[research/op-kappa-sigma-2body-PN-2026-05-09/]] STRUCTURAL DERIVED
(heuristic numerical, 7/7 sympy):
```
κ_σ ≈ 1/(3π) ≈ 0.106
```
Source: 1/π factor z orbital phase averaging × 1/3 z σ-trace structure.

#### 3.6.10.3 — JOINT RESULT (cykle #1 + #2)
**π factory cancelają cleanly** z dwóch niezależnych źródeł:
```
c_0 · κ_σ = 4π · 1/(3π) = 4/3 EXACT (clean π cancellation)
```
**Reproduces Phase 4 zero-β_ppE target structurally** — strong consistency,
NIE post-hoc fit. GW150914 6% deviation = real-world calibration artifact.

#### 3.6.10.4 — N14 LIGO scalar mode AMENDED 2026-05-09 (R5 risk RESTORED → **RESOLVED** post-T3.4 amendment, see §3.6.10.6)
[[research/op-scalar-mode-LIGO-bound-2026-05-09/]] **DOWNGRADED 2026-05-09**
do STRUCTURAL_CONDITIONAL po rigorous verification w
[[research/op-h-TT-calibration-2026-05-09/]] (8/8 PASS Phase 2 amendment).

**Phase journey:**
- Phase 1 (5/5): R5 risk identified (h_S/h_T ≈ 2√π ≈ 3.5 naive)
- Phase 2 (6/6): STRUCTURAL_NO_GO at linearized
- Phase 3 (9/9): "RE-EXAMINATION via §6.4" sphere-average multipole argument **— INCORRECT**
- **Calibration cycle Phase 1+2 (8/8 amendment):** **Phase 3 error rigorously identified**

**Klucz błędu Phase 3:** sphere-averaged ⟨δΦ⟩ = 0 ≠ h_S(observer) at given LIGO detector.

**RIGOROUS VERIFICATION (calibration Phase 2 sympy 8/8 PASS):**
TT-projection of δ^ij·X (isotropic spatial scalar metric perturbation):
```
Λ^ij_kl·δ^kl = P^ij·(1 - tr(P)/2)·X
              = P^ij·X·(1 - 1) [for unit n: tr(P) = 3 - n² = 2]
              = 0 IDENTICALLY
```

⟹ TGP linearized z δg_eff^ij = δ^ij·b_1·δΦ ansatz **NIE może** produkować
   h_+, h_× modes at observer.

| Status | Original Phase 3 (incorrect) | Amended (rigorously verified) |
|---|---|---|
| h_S at observer | 0 EXACT | b_1·δΦ(obs) ≠ 0 generically |
| h_+, h_× | "z l=2 multipole non-zero" | **0 IDENTICALLY** at linearized |
| LIGO bound | trivially satisfied | **VIOLATED ~70×** at linearized |
| R5 risk | mitigated | **RESTORED active** |

**Implikacje dla framework:**
- emergent-metric Phase 4 Path 2: STILL VALID dla 1PN/2PN (γ=β=1, Newton, EP)
- emergent-metric Phase 4 Path 2: **THREATENED** at GW polarization level
- Multi-session escape routes needed (σ at 3PN+, nonlinear δΦ self-coupling,
  velocity-dependent g_eff coupling, multi-Φ extension)

**Strukturalny issue:** w TGP single-field linearized, source coupling q·ρ·δΦ/Φ_0
jest scalar-only. NIE ma tensor coupling at LINEAR order. To FUNDAMENTAL
różnica z GR (gdzie T_ij ma tensor structure dziedziczona przez h_ij).

**Honest reporting (CALIBRATION_PROTOCOL):** original Phase 3 verdict was
incorrect; this amendment based on rigorous TT-projection identity. NO
framework-protection bias — adversarial calibration analysis correctly
caught error before propagation.

#### 3.6.10.5 — Cumulative status post-follow-up (AMENDED 2026-05-09)

| Component | Status |
|---|---|
| c_0 (heuristic) | ≈ 4π z OP-7 T3.4 chain |
| κ_σ (heuristic) | ≈ 1/(3π) z orbital averaging |
| c_0·κ_σ | **4/3 EXACT** (Phase 4 target reproduced) |
| **N14 LIGO scalar mode** | **AMENDED 2026-05-09 (morning):** R5 risk RESTORED — DOWNGRADED to STRUCTURAL_CONDITIONAL.<br/>**RE-AMENDED 2026-05-09 (evening, see §3.6.10.6):** R5 risk **RESOLVED** post-T3.4 ξ_eff factor-4 correction; cycle UPGRADED to **STRUCTURAL DERIVED** |
| h_+, h_× TT modes (linearized) | **0 IDENTICALLY** via δ^ij·b₁·δΦ ansatz alone (rigorously: TT projection of δ^ij·X = 0); **non-zero via σ-channel** post-T3.4 amendment (h_TT^σ = h_TT^GR EXACTLY at leading PN, see §3.6.10.6) |
| **6/6 P-requirements** | morning state: **5/6 RESOLVED** (P6 z R5 risk active).<br/>**evening state post-T3.4 amendment: 6/6 RESOLVED** (P6 RESOLVED via Path A direct + ξ_eff = 4·G·Φ_0², see §3.6.10.6) |

**Cumulative sympy: 105/105 PASS** (57 emergent-metric + 5+7+20 follow-up + 8+8 calibration cycle amendment).
**Updated 2026-05-09 (evening) post-T3.4 amendment cascade: 157/157 PASS** (105 prior + 24 σ-3PN Phase 2 + 11 σ-3PN Phase 1 foundation + 17 T3.4 amendment Phase 1 = +52). See §3.6.10.6.

**Honest caveats (multi-session deferred):**
- Precise c_0 numerical pinning (OP-7 T3.4 to higher precision)
- Precise κ_σ z explicit Hadamard regularization 2-body PN
- **R5 risk RESTORED at linearized level — NOT calibration question, deeper structural issue**
- Multi-session escape routes (σ at 3PN+, nonlinear δΦ, framework extension)

**Strukturalnie (AMENDED):** post-falsification recovery dla TGP gravity sektora
**STRUCTURAL_CONDITIONAL** at qualitative level. Recovery framework solves:
- ✅ 1PN/2PN observables (γ=β=1, Cassini, Mercury, Newton)
- ✅ 2.5PN orbital binding (β_ppE^new compliance window)
- ✅ Equivalence principle (m_b=m_g automatic)
- ⚠️ **NIE solves GW polarization mode content** (linearized predicts scalar-only,
  observation TT-dominant) — multi-session escape route needed

**Conservative position:** TGP framework PASSES all 1PN/2PN/2.5PN observable
tests, ALE wymaga multi-session work dla pełnego matching z GW polarization
content observed by LIGO. Status STRUCTURAL_CONDITIONAL accurately reflects this.

> ⚠ **Post-T3.4-amendment update (2026-05-09 evening):** Status powyżej
> ("STRUCTURAL_CONDITIONAL", "5/6 RESOLVED", "R5 risk active") opisuje
> stan PRZED amendment cascade z §3.6.10.6. Po T3.4 ξ_eff factor-4
> korekcie framework jest **STRUCTURAL DERIVED 6/6 RESOLVED** (R5 RESOLVED).
> Sekcje §3.6.10.4 i §3.6.10.5 pozostają jako audit-trail morning state;
> aktualne LIVE status w §3.6.10.6.

#### 3.6.10.6 — T3.4 ξ_eff normalization amendment 2026-05-09 (evening) — R5 RESOLVED

[[research/op-T34-normalization-amendment-2026-05-09/]] **STRUCTURAL DERIVED**
(17/17 PASS Phase 1) — clean first-principles re-derivation z explicit factor
tracking (Misner-Thorne-Wheeler 1973 §36, Maggiore 2008 §3.1-3.4, Wald 1984 §11.2)
identifies factor-4 algebraic gap w OP-7 T3.4 source.

**Trigger sequence (multi-cycle resolution path):**

1. [[research/op-h-TT-calibration-2026-05-09/]] (16/16 PASS) — adversarial
   calibration cycle wymusiła rigorous audit § 3.6.10.4 (caught Phase 3 cycle #3
   sphere-average error)
2. [[research/op-sigma-3PN-radiative-2026-05-09/]] Phase 2 (24/24 PASS) — Path A
   direct h_TT^σ amplitude derivation, ujawniła h_TT^σ / h_TT^GR ≈ 0.265 z
   literal LOCKS (factor 1/4 gap)
3. [[research/op-sigma-3PN-radiative-2026-05-09/Phase2_adversarial_verification.md]]
   — independent agent confirmed compound factor-4 gap (Gap 1 line 132 + Gap 2
   line 140 w `op7_t3_4_xi_coupling.py`)
4. [[research/op-T34-normalization-amendment-2026-05-09/Phase1_sympy.py]]
   (17/17 PASS) — clean re-derivation z first principles, **NIE inheriting**
   żadnej z trzech istniejących wartości ξ_eff w cycle chain

**Three-source disambiguation (z chain):**

| Source | Pre-amendment value | Post-amendment value |
|---|---|---|
| OP-7 T3.4 text | `ξ_eff = G·Φ_0²` ✗ (factor 1/4 low) | `ξ_eff = 4·G·Φ_0²` ✓ |
| c_0 cycle Phase 1 sympy line 65 | `ξ_eff = 4π·G·Φ_0²` ✗ (factor π high) | `ξ_eff = 4·G·Φ_0²` ✓ |
| Clean re-derivation (this) | — | **`ξ_eff = 4·G·Φ_0²` (correct)** |

**Matching condition (derived from MTW + Maggiore textbook chain):**
```
c_0 · ξ_eff = 16π · G · Φ_0²        [first-principles, h_TGP^σ = h_TT^GR EXACTLY]
```
Z `c_0 = 4π` LOCK (cycle #1, preserved): `ξ_eff = 4·G·Φ_0²`.

**Identified gaps w `op7_t3_4_xi_coupling.py` (Gap 1 + Gap 2 = factor 4):**

- **Gap 1 (line ~132):** missing PN-(1/2) factor z standard identity
  `∫T^ij d³y = (1/2)·d²Q^M_ij/dt²` (Maggiore Eq. 3.81). Stated:
  `σ_ab(r,t) ~ -(ξ/(4π·c⁴))·Q̈_ab^TT/r`; correct:
  `σ_ab = -(ξ/(8π·c⁴·r))·Q̈_ab^TT`. **Factor 2 low.**
- **Gap 2 (line ~140):** algebraic mismatch z line 139 (h_GR ma explicit factor 2).
  Stated: `Λ_0·ξ = 4πG`; correct (z line 139 equation): `Λ_0·ξ = 8πG`. **Factor 2 low.**
- **Compound:** `2 × 2 = 4`. Konsystentne z Phase 2 σ-3PN finding.

**Preserved LOCKS (NIE affected — single-coefficient correction):**

| LOCK | Value | Status |
|---|---|---|
| c_0 | 4π | ✓ preserved (cycle #1) |
| κ_σ | 1/(3π) | ✓ preserved (cycle #2) |
| c_0·κ_σ | 4/3 EXACT | ✓ preserved (joint LOCK) |
| β_ppE | 0 (z 2.5PN GR phase match) | ✓ preserved (Phase 4) |
| γ_PPN = β_PPN | 1 EXACT | ✓ preserved (1PN/2PN tests) |
| m_inertial = m_grav | automatic z S05 | ✓ preserved (WEP) |
| **ξ_eff** | **4·G·Φ_0² (was G·Φ_0²)** | **AMENDED** |

**Post-amendment framework status:**

| # | P-requirement | Pre-amendment | Post-amendment |
|---|---|---|---|
| P1 | Newton I, II structural | ✓ | ✓ |
| P2 | 1PN/2PN tests (Cassini, Mercury, LLR) | ✓ | ✓ |
| P3 | β_ppE GWTC-3 1σ window | ✓ | ✓ |
| P4 | m_b = m_g (S05) | ✓ | ✓ |
| P5 | ngEHT consistency (Sgr A*, M87*) | ✓ | ✓ |
| **P6** | **LIGO O3 polarization + amplitude (R5)** | ⚠️ active | **✓ RESOLVED** |

**6/6 P-requirements RESOLVED post-T3.4-amendment.**

**Cycle status upgrades (cascade):**

| Cycle | Pre-amendment | Post-amendment |
|---|---|---|
| op-sigma-3PN-radiative Phase 2 | STRUCTURAL_CONDITIONAL (norm gap) | **STRUCTURAL DERIVED** |
| op-scalar-mode-LIGO-bound (#3) | STRUCTURAL_CONDITIONAL (R5 active) | **STRUCTURAL DERIVED** (R5 RESOLVED) |
| op-h-TT-calibration | STRUCTURAL_CONDITIONAL_HALT | preserved (adversarial-historical artifact) |
| op-emergent-metric Phase 4 (Path 2) | STRUCTURAL DERIVED qualitatively | **STRUCTURAL DERIVED quantitatively** (full GR amplitude reproduction) |
| op-T34-normalization-amendment | — (this cycle) | **STRUCTURAL DERIVED** (17/17 PASS) |

**Smoking-gun predictions post-amendment (testable across regimes):**

- `h_TT^σ = h_TT^GR` EXACTLY at leading PN order (verified Phase 1 + Phase 2 σ-3PN)
- `β_ppE = 0` at 2.5PN phase (Phase 4 emergent-metric LOCK preserved)
- 2PN deviation amplitude/phase **~0.02 rad** at LIGO O5+ (M9.1''-specific)
- `m_σ² = 2m_s²` dispersion ≈ (0.71 meV)² testable via Cosmic Explorer (~2030)
- ngEHT photon ring **+14.6%** at M87*/Sgr A* (M9.1'' static spherical, σ-sector independent)
- GW170817 `c_GW = c` structural (Lorentz-violation absent)

Każda predykcja jest **niezależną obserwacją** i **niezależnym falsyfikatorem**.

**Adversarial verification protocol — value DEMONSTRATED 2× this day:**

- 1× w `op-h-TT-calibration` cycle: caught Phase 3 cycle #3 "h_+ ≠ 0 z l=2 multipole"
  sphere-average error
- 1× w `op-T34-normalization-amendment` cycle: caught compound factor-4 ξ_eff gap
  w T3.4 algebraic chain

Bez protokołu framework byłby opublikowany z normalization inconsistency.
Maintain CALIBRATION_PROTOCOL §4.3 adversarial commitment jako default w **wszystkich**
quantitative cycles (binding rule, [[meta/CALIBRATION_PROTOCOL.md]]).

**Cumulative sympy post-T3.4-amendment: 157/157 PASS**
(105 prior + 11 σ-3PN Phase 1 + 24 σ-3PN Phase 2 + 17 T3.4 amendment Phase 1).

**Updated 2026-05-09 (post-σ-3PN-Phase-3): 176/176 PASS** (157 prior + 19 σ-3PN
Phase 3 higher-PN four-channel structural derivation). σ-channel of TGP h_TT^σ
locked z GR through 2PN amplitude (non-hereditary multipoles). **Channel B
audit flag preserved** (m_σ ≈ 0.71 meV ≫ ℏω_LIGO ~ 4·10⁻¹³ eV Yukawa concern;
4 resolution mechanisms identified; sub-cycle pending `op-sigma-yukawa-audit`).
[[research/op-sigma-3PN-radiative-2026-05-09/Phase3_results.md]] dla detail.

**Updated 2026-05-09 (post-Yukawa-audit Phase 1): 211/211 PASS** (176 prior +
35 Yukawa audit Phase 1). [[research/op-sigma-yukawa-audit-2026-05-09/Phase1_results.md]]
Adversarial cycle ujawnił: Phase 2 + T3.4 amendment wykorzystały **massless retarded
Green explicitly**; przy m_σ ≈ 0.71 meV ≫ ℏω_LIGO calculation jest formal m → 0
limit, NIE direct LIGO observable. **4-mechanism verdict:**
- (i) Goldstone: NO realization (Z₂ discrete, no continuous symmetry breaks)
- (ii) composite: δŝ itself heavy (m_s ≈ 0.5 meV) → does NOT resolve
- (iii) emergent-metric δΦ-mediation: **PLAUSIBLE** pending m_Φ at level 0
  verification (m_Φ ~ Λ_cosm ~ 10⁻³³ eV would give λ_C ~ Hubble scale, no
  Yukawa suppression in observable universe; nonlinear (∂Φ)² → σ effective
  composite would mediate h_TT)
- (iv) Path A as effective contact: INTERPRETIVE (Phase 2 formula = formal
  matching condition; combines z iii dla full physical picture)

**Conservative classification recommendation adopted:** framework status preserved
STRUCTURAL DERIVED z **explicit Yukawa-resolution-pending caveat** (calculations
remain mathematically valid; classification refined NIE invalidated). Aggressive
alternative (DOWNGRADE) deferred pending mechanism (iii) m_Φ verification work.

**Updated 2026-05-09 (post-mPhi-verification Phase 1): 235/235 PASS** (211 prior +
24 mPhi-verification Phase 1). [[research/op-mPhi-level0-verification-2026-05-09/Phase1_results.md]]
Phase 1 sympy zweryfikował V''(ψ=2/3) = (4/3)·γ EXACT dla
V_M9.1''(ψ) = -γ·ψ²·(4-3ψ)²/12 (G.0 closure 2026-05-02 LOCK form).
**Result:** m_ψ² = (4/3)·M_Pl²·g̃ → m_ψ = (2/√3)·√g̃·M_Pl ≈ 1.155·M_Pl
≈ **1.41·10²⁸ eV** (at g̃=1). **Verifies op-Phi-vacuum-scale Phase_FINAL §2.1 line
99 'm_ψ ~ M_Pl' claim** via clean derivation.

**Numerical scale comparison:** m_ψ/ℏω_LIGO ≈ 3.5·10⁴⁰ (mechanism iii prerequisite
VIOLATED by 41 orders of magnitude); λ_C(m_ψ) ≈ Planck length; D/λ_C at LIGO 1 Gpc
distance ≈ 10⁶⁰ — Yukawa suppression exp(-D/λ_C) ≈ exp(-10⁶⁰+).

**Verdict on mechanism (iii):**
- At **FALSIFIED V_M9.1''** (specific (4-3ψ)/ψ form, RULED OUT 5σ by GWTC-3):
  mechanism (iii) **FAILS**
- At **RECOVERY V parametric family** (post-emergent-metric Phase 4 zero-β region):
  **OPEN question** (multi-session continuation, P1 priority)

**Framework cascade DOWNGRADE applied** (analog T3.4 amendment cycle pattern but
reversed direction):
- σ-3PN Phase 2 + amendment + Phase 3: STRUCTURAL DERIVED z caveat → **STRUCTURAL_CONDITIONAL pending recovery V**
- op-scalar-mode-LIGO-bound (#3): R5 RESOLVED conditional → **R5 RESTORED** at LIGO amplitude level
- **6/6 → 5/6 P-requirements RESOLVED** (P6 z R5 risk active)
- **Calculations preserved:** 235/235 sympy PASS (calculations remain mathematically valid w stated framework; classification refined honestly per CALIBRATION_PROTOCOL §4.3 pattern)

**P1 OPEN PATH (multi-session):** explicit recovery V form analysis in zero-β
region of emergent-metric Phase 4 parametric family. If ANY zero-β-compatible
V has near-degenerate minimum (V''(Φ_0) ≪ ℏω_LIGO ~ 4·10⁻¹³ eV) → mechanism (iii)
realizes for that V → framework status restorable do STRUCTURAL DERIVED.
If ruled out → mechanism v (framework extension: additional massless tensor mode,
nonlinear δΦ products beyond level 0) — multi-session deep theoretical work.

**Honest scientific outcome.** Adversarial verification protocol value
DEMONSTRATED **5× this day**. Pattern: each audit step identifies hidden
assumption before publication-grade claims propagate. Maintain CALIBRATION_PROTOCOL
§4.3 default w wszystkich quantitative cycles jako binding rule.

**Strukturalnie post-amendment:** TGP gravity sector post-falsification recovery
**STRUCTURAL DERIVED** (qualitative + quantitative). Framework reproduces
GR-equivalent quadrupole formula z explicit factor calibration; smoking-gun
deviations testable at LIGO O5+, ngEHT, Cosmic Explorer scales.

**Honest reporting (CALIBRATION_PROTOCOL):** original T3.4 (2026-04-25) closure
preserved jako historical artifact; amendment notice issued w
[[research/op7/OP7_T3_results.md]] §0 i `op7_t3_4_xi_coupling.py` top-of-file.
Single algebraic correction sufficient — żadnej framework-protection retreat,
żadnego multi-candidate fit. Pattern: "Full honest amendment" preferred over
framework protection (analog Errata + Amendment w published papers).

## 4. Co to jest "materia" i co to jest "źródło"

Trzy warstwy operacyjne, **każdą z innym statusem**:

| Warstwa | Co to jest | Sprzęga z Φ przez | Status |
|---|---|---|---|
| **3a — Pola materii** ψ_m (fermiony, cechowania, płyn) | Klasyczne i kwantowe pola standardowego modelu, żyjące na `g_eff` | `g_eff^μν[Φ]` w `L_mat`, **nie** `Φ` bezpośrednio | (P) |
| **3b — Gęstość ρ** | Skalarna gęstość w `L_mat = -(q/Φ_0) φ ρ` | Bezpośrednie sprzężenie minimalne (czynnik φ) | (W) — to jest źródło dla Φ-EOM |
| **3c — Kinki / defekty** (cząstka = radialny kink Φ + topologia chiralna) | **Hipoteza/roadmap** alternatywnego opisu fermionów jako struktur w samym Φ; emergent Dirac propagator. **Strukturalny szkic 5-fazowy zamknięty 2026-05-01** w `research/why_n3/` (Phase 1: ψ↔g₀ liniowa identyfikacja, bariera R3 ≡ M9.1'' Lorentzian horizon ψ=4/3; Phase 2: n(α)=e²(1−α/4); Phase 3: spin-1/2 z RP² + Berry phase π; Phase 4: Yukawa + m₀=0 vacuum; Phase 5: pełen S_TGP(p;ψ)). **POST-2026-05-16 STATUS UPDATE — Phase 6+ partial closure:** L08 audit problems #1 (spin-statistics theorem) i #4 (Dirac algebra Clifford) **OPERATIONALLY CLOSED A−** via [[research/op-L08-Phase6-FR-antisymmetry-2026-05-16/]] (Finkelstein-Rubinstein 2-particle exchange π₁(C_2-defect)=Z₂ → Berry phase π → antisymmetric wavefunction) + [[research/op-L08-Phase6-Clifford-emergence-2026-05-16/]] (Cl(1,3) algebra inherited z M9.1'' Lorentz tetrad geometry). Problem #2 (e²/4, three generations) **PARTIAL B+** z algebraic reconciliation F1↔F2 + e_Euler² classification as NUMERICAL ANCHOR (PHASE6 §11 reinforced 2026-05-16). L05 mass-exponent k(α,d) z m_obs vs M_full distinction **CLOSED A−** via [[research/op-L05-mass-exponent-k-alpha-d-2026-05-16/]]. **POST-2026-05-17 NEUTRINO LINE:** problem #3 sub-component **neutrinos A− REINFORCED** via sesja 7-cycle z PR-016 dual-scenario + 7-bound astrofizyczny survey ([[research/op-neutrino-mu-nu-astrophysical-discrimination-2026-05-17/]]). **POST-2026-05-18 W/Z DECLARED LIMIT:** problem #3 sub-component **bosons (W/Z, SU(2)) — DECLARED THEORETICAL LIMIT** post 5-path exhaustion (cycle 6 paths α/β/γ/δ + cycle ε composite Higgs sesja-1-of-N 2026-05-18 [[research/op-composite-higgs-substrate-attempt-2026-05-18/]]). Combined disposition Option A + Option C adopted w [[meta/TGP_W_Z_THEORETICAL_LIMIT.md]] — W/Z + Higgs doublet structure jako input phenomenology (analog SM treating Higgs); future structural-extension Option B preserved optional. **Status warstwy 3c: (H) hipoteza → partial-(D) z DECLARED LIMIT dla SU(2) post-2026-05-18** (quarks A− + neutrinos A− REINFORCED + bosons DECLARED LIMIT + spin-statistics + Dirac algebra CLOSED A−) | (problem #3 boson: declared limit; quarks + neutrinos closed) — `rem:materia-hierarchia` | **(Hipoteza → partial-(D) z DECLARED LIMIT post-2026-05-18)** |

**Operacyjne źródło dla Φ-EOM** to **ρ** (warstwa 3b). Punkt do
których wraca każdy rachunek typu Newton/PPN/GW.

Warstwa 3c (kinki) jest **długoterminową hipotezą** — TGP ma w
przyszłości pokazać, że ψ_m emerguje z topologii Φ — ale dziś **nie
zastępuje** ρ w bieżących rachunkach.

## 5. Grawitacja w TGP — co to jest, czym **nie** jest

### 5.1 Czym **nie** jest

- **Nie ma grawitonu** — ani fundamentalnego, ani kompozytowego.
- **Nie ma propagującej cząstki spin-2.**
- **Nie ma teorii skalarno-tensorowej** w stylu Branstr-Dicke / Horndeski
  (`rem:not-scalar-tensor`, `sek08a` lin. 129–151).
- **Metryka `g_eff` nie jest niezależną zmienną dynamiczną** — jest
  wynikowo zdefiniowana przez Φ.

### 5.2 Czym **jest**

> Grawitacja w TGP jest **kolektywnym efektem fluktuacji jednego pola
> substratu Φ**. Wszystkie zjawiska grawitacyjne (Newton, GR weak field,
> GW, soczewkowanie, dylatacja czasu, redshift) są **wzorcami statystycznymi
> długo-falowych zaburzeń pola Φ** wokół klasycznego rozwiązania
> równowagi `Φ_eq[ρ]`.

To jest tradycja **emergentnej grawitacji** (Sakharov 1968, Verlinde 2010,
Padmanabhan, Volovik). Brak grawitonu **nie jest** brakiem — jest
**cechą strukturalną** programu, eliminującą problem renormalizacji
grawitonu i automatycznie respektującą twierdzenie Weinberga-Wittena
(które dotyczy fundamentalnych cząstek, nie kompozytów ani efektów
kolektywnych).

### 5.3 Cel teoretyczny

> **GR w limicie**, jako analog numeryczny — nie izomorfizm analityczny.

TGP **nie musi** wyprowadzać równań Einsteina jako równanie Lagrange'a
dla `g_eff`. **Musi** dawać zgodne z GR predykcje liczbowe w odpowiednim
regimie (system słoneczny + LIGO/Virgo + EHT, w precyzji obserwacji).
**Odchylenia** od GR są **dozwolone** poza tym regimem (np. UV substrat,
kosmologiczne ekstremalne skale).

To jest **mocniejsza** pozycja teoretyczna niż "TGP = GR", nie słabsza:
otwiera drzwi falsyfikacji w skalach UV i daje konkretne predykcje
różnicujące.

## 6. Pęd i bezwładność — obraz Lenz-podobny

### 6.1 Fundamentalny obraz fizyczny (kanoniczny dla TGP)

> **Spoczynek = lokalna równowaga pola Φ wokół źródła ρ.** Przyspieszenie
> źródła łamie tę równowagę → pole reaguje zwrotnie → siła back-reakcji
> opiera się zmianie konfiguracji → manifestuje się jako bezwładność.
> Stały ruch = pole "ślizgające się" wraz ze źródłem, brak back-reakcji,
> zachowanie pędu (analog reguły Lenza).

### 6.2 Konsekwencje strukturalne

| Zjawisko | Standardowy obraz | Obraz TGP |
|---|---|---|
| **Bezwładność** | Aksjomat (Newton I) | Back-reakcja Φ-pola na zmianę konfiguracji równowagi |
| **Masa bezwładnościowa** | Niewyjaśniony parametr | Współczynnik back-reakcji, obliczalny z linearyzacji Φ-EOM wokół `Φ_eq` |
| **Masa grawitacyjna** | Sprzężenie z polem grawitacyjnym | Siła sprzężenia źródła z Φ (`q/Φ_0`) |
| **Zasada równoważności** (m_b = m_g) | Postulat / dyfeo. niezmienniczość | **Automatyczna** z konstrukcji jednopolowej (oba pochodzą z tej samej stałej q) |
| **Newton I (brak tarcia)** | Postulat | Translacyjna niezmienniczość minimum lokalnej równowagi → brak preferowanego układu odniesienia → brak back-reakcji w stanie ustalonym |
| **Newton II (F = ma)** | Postulat | Liniowa odpowiedź pola na zmianę konfiguracji równowagi |
| **Promieniowanie z przyspieszającego źródła** | Larmor (EM), kwadrupol (GR) | Energia uciekająca jako fluktuacje Φ na nieskończoność (analog) |

### 6.3 Skala mikro vs makro (interpretacja kwantowa — robocza)

W skali mikro fluktuacje pola Φ są porównywalne z lokalizacją źródła
(cząstki). "Pomiar pozycji cząstki" = pomiar centroidu fluktuującej
konfiguracji pola → emergentna nieoznaczoność `Δx · Δp ≥ ℏ/2`. To
jest stochastyczna interpretacja QM (Nelson; SED Marshalla/Boyera) —
**spekulatywna, ale spójna z architekturą TGP, do badania w późniejszych
cyklach**.

W skali makro wiele źródeł generuje uśrednione pole efektywne →
fluktuacje uśredniają się → klasyczny limit.

## 7. Status M3–M8 (archiwum)

W cyklach M3–M8 (kwiecień 2026) zbadano własności **wykładników
krytycznych** ratio `β/γ` w punkcie stałym Wilsona-Fishera 3D Isinga
metodami MK-RG (M3) i NPRG/Wetterich-LPA (M8), wraz z trzema kanałami
zamknięcia (M4 H-S Jacobian, M5 Z_Φ, M6/M7 wiązanie GL).

**Wnioski techniczne (zachowują wartość):**
- W jednoskładnikowym skalarnym Z₂ na poziomie LPA β/γ w punkcie stałym
  WF jest ujemne (≈ −0.3 do −0.5) w obu schematach RG.
- Wszystkie trzy single-channel kanały zamknięcia zostały sfalsyfikowane.

**Wnioski reinterpretacyjne (pod nową ramą "GR w limicie"):**
- M3–M8 odpowiadało na pytanie o **uniwersalność klasy krytycznej**
  (wykładniki w punkcie stałym WF), nie na pytanie o **fenomenologię
  TGP** (klasyczna dynamika pola średniego w fazie złamanej, na
  skończonej skali korelacji).
- TGP fizycznie nie żyje na T = T_c (faza krytyczna), tylko w **fazie
  złamanej** (`m₀² < 0`, `v² = |m₀²|/λ₀`), z **klasyczną dynamiką
  Φ-EOM** plus poprawki z fluktuacji.
- Wartość `β/γ` w punkcie stałym WF **nie jest** liczbą wiążącą dla
  fenomenologii TGP. **Wiążące jest:** liczbowa zgodność rozwiązań
  Φ-EOM z fenomenologią GR w odpowiednim regimie (Newton + γ/β_PPN +
  GW kwadrupol + c_GW=c).

**Status:** M3–M8 zachowane jako "dobrze wykonane na poprzednim torze",
ich wynik nie jest używany jako kryterium otwartego/zamkniętego dla
fenomenologii TGP. Sekwencja M9+ rozpoczyna nowy tor.

## 8. Aktualny program badawczy — M9+ ("klasyczna dynamika i pęd")

**Punkt wyjścia:** Φ-EOM (`eq:field-eq-reproduced`) jako **klasyczne
równanie pola w fazie złamanej**, ze źródłem `ρ` (warstwa 3b).

**Cel:** sprawdzenie liczbowe, że rozwiązania Φ-EOM odtwarzają
fenomenologię GR w precyzji obserwacji, w trzech sytuacjach:

- **M9.1 — Statyka i potencjał Newtona/PPN.** Sferyczne źródło,
  słabe pole, dopasowanie `q ↔ G`, weryfikacja `γ_PPN = 1` (dokładnie
  z konstrukcji metryki) i `β_PPN = 1` (z formy eksponencjalnej).

- **M9.2 — Pęd i bezwładność (Lenz-podobnie).** Statyczne i poruszające
  się źródło, back-reakcja pola na przyspieszenie, identyfikacja masy
  bezwładnościowej z back-reakcją, weryfikacja zasady równoważności.

- **M9.3 — Promieniowanie GW z dynamicznego źródła.** Linearyzacja
  Φ-EOM wokół `Φ_eq`, wzór kwadrupolowy, polaryzacje fluktuacji,
  porównanie z LIGO scalar-mode bound (< few %) i GW170817
  (`c_GW = c` do 10⁻¹⁵).

Pliki wykonawcze: `research/op-newton-momentum/M9_*.md` i `*.py`.

## 9. Reguły dla agentów

1. **Przed jakąkolwiek pracą obliczeniową w `research/`**: przeczytać
   ten dokument oraz odnośne `sek08*.tex`.
2. **Każda nowa hipoteza/test** musi się odnosić do hierarchii (poziom
   0/1/2/3) i deklarować, na którym poziomie operuje.
3. **Pivot substratu jest dozwolony** w obrębie aksjomatu §1 (jedno pole
   skalarne Z₂). Pivot poza ten aksjomat jest **odrzucony bez
   dyskusji** — wymaga osobnej rozmowy z autorem.
4. **GR jako numeryczny analog, nie izomorfizm.** Każdy test grawitacyjny
   formułować jako zgodność liczbową w precyzji obserwacji, nie jako
   próbę wyprowadzenia analitycznego równań Einsteina.
5. **Pęd jako Lenz-podobny back-reakcja Φ-pola.** Każda interpretacja
   inercji/masy/pędu musi być zgodna z §6 — w razie wątpliwości,
   poprosić autora o klaryfikację.
6. **Nie wprowadzać grawitonu** — ani fundamentalnego, ani
   kompozytowego (np. `σ_μν = ∂φ ⊗ ∂φ` jako "emergentny grawiton" jest
   **niezgodne** z §5: TGP nie używa cząstkowej ramy dla grawitacji).

## 10. Pliki referencyjne (canonical sources of truth)

- `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` —
  akcja zunifikowana, wariacja → Φ-EOM, β=γ jako warunek próżni.
- `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` —
  metryka efektywna z budżetu substratu, weryfikacja PPN.
- `core/sek08_formalizm/sek08_formalizm.tex` — pełny formalizm,
  hierarchia poziomów, sprzężenie metryczne, hipoteza materii-hierarchii.
- `KNOWN_ISSUES.md` — historia decyzji, statusy testów,
  sfalsyfikowane warianty.
- `research/op1-op2-op4/M{3..8}_*.md` — archiwum poprzedniego cyklu
  (krytyczne RG, jednoskładnikowy skalar, β/γ przy WF FP).
- `research/op-newton-momentum/M9_*.md` — bieżący cykl (klasyczna
  dynamika, pęd, GW jako analog numeryczny).

---

**Ostatnia aktualizacja:** 2026-04-25, po rozmowie ustalającej obraz
"grawitacja = efekt fluktuacji pola, brak grawitonu, pęd = Lenz-podobny
back-reakcja". Następne aktualizacje: po istotnych decyzjach
ontologicznych lub interpretacyjnych. **Nie** aktualizować po
pojedynczym wyniku liczbowym (te idą do `KNOWN_ISSUES.md` i
odnośnych `M*_results.md`).
