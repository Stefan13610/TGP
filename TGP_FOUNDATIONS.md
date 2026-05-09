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
| **3** | Materia i pola cechowania | Sprzężenie metryczne, Dirac, U(1)×SU(2)×SU(3), generacje | (P) |

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

### 3.6.2 PPN constraints jako derivation (NIE postulat)

Phase 2 cyklu (sympy 7/7 PASS) ustaliło z action variation:

```
γ_PPN = 1   ⟺   b_1 = -a_1     (1PN constraint na Taylor coefs)
β_PPN = 1   ⟺   ξ_2 = ξ - a_2·ξ³/2     (2PN consistency)
```

**Solar system constraints** (Cassini |γ−1| ≤ 2.3·10⁻⁵, Mercury |β−1| ≤ 8·10⁻⁵)
strukturalnie satysfakcjonowane: γ = β = 1 EXACT z constraintów strukturalnych,
NIE z empirycznego fitu specyficznej formy.

**σ-coupling C(ψ) jest WOLNA w 1PN/2PN** — wnosi parametryczną swobodę
dopiero przy 2.5PN binary inspiral.

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

**Open (post-cycle, deferred):** N14 LIGO scalar mode amplitude (R5 risk, multi-session).

## 4. Co to jest "materia" i co to jest "źródło"

Trzy warstwy operacyjne, **każdą z innym statusem**:

| Warstwa | Co to jest | Sprzęga z Φ przez | Status |
|---|---|---|---|
| **3a — Pola materii** ψ_m (fermiony, cechowania, płyn) | Klasyczne i kwantowe pola standardowego modelu, żyjące na `g_eff` | `g_eff^μν[Φ]` w `L_mat`, **nie** `Φ` bezpośrednio | (P) |
| **3b — Gęstość ρ** | Skalarna gęstość w `L_mat = -(q/Φ_0) φ ρ` | Bezpośrednie sprzężenie minimalne (czynnik φ) | (W) — to jest źródło dla Φ-EOM |
| **3c — Kinki / defekty** (cząstka = radialny kink Φ + topologia chiralna) | **Hipoteza/roadmap** alternatywnego opisu fermionów jako struktur w samym Φ; emergent Dirac propagator. **Strukturalny szkic 5-fazowy zamknięty 2026-05-01** w `research/why_n3/` (Phase 1: ψ↔g₀ liniowa identyfikacja, bariera R3 ≡ M9.1'' Lorentzian horizon ψ=4/3; Phase 2: n(α)=e²(1−α/4); Phase 3: spin-1/2 z RP² + Berry phase π; Phase 4: Yukawa + m₀=0 vacuum; Phase 5: pełen S_TGP(p;ψ)). Otwarte (Phase 6+): X=e²/4 analytic derivation, A^(5−α) vs A²·g₀^(e²/2) reconciliation dla τ/e | (otwarty problem) — `rem:materia-hierarchia` | (Hipoteza, strukturalny szkic CLOSED) |

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
