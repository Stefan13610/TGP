# TGP — Plan Domknięcia Teorii (Master Closure Plan)

**Data:** 2026-04-25
**Status:** STRATEGY DOC (post M9.1''+ OP-7 T1+T2 closure)
**Autor:** Synteza 5 niezależnych sond audytowych (2026-04-25)
**Następne zatwierdzenie:** wymaga akceptacji użytkownika przed egzekucją.

---

## 0. Streszczenie wykonawcze (TL;DR)

TGP po pivotach 2026-04-24 (v2 GL-substrate) → 2026-04-25 (M9.1'' hyperbolic metric) → 2026-04-25 (OP-7 T1+T2 kinematics) **stoi w punkcie przełomowym**:

- **Strukturalnie zdrowa**: aksjomaty single-Φ Z₂ + emergentna grawitacja są spójne; M9.1'' i OP-7 T1+T2 wzmocniły je, nie zburzyły.
- **Falsyfikowalna**: explicit 2PN deviation `|Δg_tt| = (5/6)U³`, scalar-only GW (jeszcze) niezgodne z GW170817 czeka na OP-7 T3-T6.
- **W długu technicznym**: ~18 dokumentowych sprzeczności + ~20 skryptów do depekacji/rewrite po pivotach. Paper jest 90% dopasowany, ale §2 σ_ab formula stale.

**Główna oś krytyczna do domknięcia:**
```
OP-7 T3 (dynamika σ_ab)
  → OP-7 T4 (Λ(ψ) coupling)
    → OP-7 T5 (kwadrupol h+, h×)
      → OP-7 T6 (PPN/Z₂/c_GW konsystencja)
        → GW170817 unconditional + GW150914 quantitative
          → tgp_core.tex v3 (Zenodo) + EHT path
```

Sześć tygodni przy obecnym tempie, **jeśli T3 daje ghost-free EOM**. Jeśli T3 odkrywa ghost lub patologię — TGP wymaga drugiego pivotu (Λ(ψ) ad hoc albo dodatkowe pole). To jest "moment prawdy".

---

## 1. Mapa stanu teorii (2026-04-25)

### 1.1 Co zamknięte pozytywnie ✅

| ID | Co | Data | Komentarz |
|---|---|---|---|
| **M1-M8** | Wilson-Fisher RG cycle | pre-2026-04 | ARCHIVED jako pre-M9 program |
| **M9.1''** | Hyperbolic metric `g_tt = -c₀²(4-3ψ)/ψ` | 2026-04-25 | β_PPN=γ_PPN=1 exact at 1PN |
| **M9.1'' P1** | Higher PPN orders (Schwarzschild) | 2026-04-25 | c₂=-1, c₃=+5/3, c₄=-10/3 exact |
| **M9.1'' P2** | Variational principle V(Φ)/Φ⁴ | 2026-04-25 | Triple convergence (P2-C/D/E) |
| **M9.1'' P3** | Mercury, Cassini, BBN, LLR | 2026-04-25 | Weak-field PASS; GW170817 conditional |
| **M9.1'' P4** | Paper rewrite (tgp_core.tex) | 2026-04-25 | Hyperbolic form integrated |
| **OP-7 T1** | No-tensor for single-Φ M9.1'' | 2026-04-25 | 7/7 PASS; only breathing mode |
| **OP-7 T2** | σ_ab from H_Γ as composite | 2026-04-25 | 12/12 PASS; gradient strain `K_ab - (1/3)δ_ab Tr(K)` |
| **OP-6 (M1-A')** | H_GL effective axiom adopted | 2026-04-24 | Pivot from minimal bilinear (M1-B falsified) |
| **QM emergence** | 79/79 testów | pre-2026-04 | Substrate-independent |
| **Particle sector** | Lepton masses, α₃ | pre-2026-04 | Decoupled from gravity pivots |

### 1.2 Co krytycznie otwarte 🔴 (gating closure)

| ID | Co | Status | Blocker |
|---|---|---|---|
| **OP-7 T3** | Dynamika σ_ab (EOM) | OPEN | wariacyjnie z S_TGP[ŝ]; m_σ origin; T_ab^TT source |
| **OP-7 T4** | Metric coupling Λ(ψ) | OPEN | ghost-free; c_GW=c₀; redukcja do M9.1'' przy σ=0 |
| **OP-7 T5** | Kwadrupol h+, h× ∝ Q̈/r | OPEN | matching ξ_eff do amplitudy GW150914 (~1.0e-21) |
| **OP-7 T6** | PPN/c_GW/Z₂ konsystencja | OPEN | full PPN scope, ghost analysis, parity |
| **GW170817** | Unconditional pass | conditional | gated by T3-T6 |
| **GW150914** | Amplitude match | open | gated by T5 |
| **EHT** | Black hole shadow | open | gated by OP-7 + strong-field nonlinear Φ-EOM |

### 1.3 Co umiarkowanie otwarte 🟡 (programmes)

| ID | Co | Priorytet | Notatka |
|---|---|---|---|
| **OP-1** | Uniqueness U(φ) | MED | M2 derivation z H_Γ blocked on structural analysis |
| **OP-2** | β=γ all scales | LOW | CLOSED-NEG; archived |
| **OP-3** | a_Γ = 1/Φ₀ identification | LOW | Postulate, needs QFT renormalization argument |
| **OP-4** | γ = g₀^ε derivation | LOW | Blocked by OP-1 M2 |
| **OP-5** | N-body nonlinear | LOW | Not yet formalized |
| **OP-8** | U(1) gauge emergence | LOW | Deferred post-OP-7 |
| **OP-9** | Strong-field BH | MED | Gated by OP-7 + EHT path |
| **OP-10** | Quantum theory of Φ | LOW | Deferred |
| **M9.2** | Momentum/Lenz back-reaction | MED | Gated by M9.1'' (DONE) |
| **M9.3** | GW radiation | HIGH | Gated by OP-7 T5 |

### 1.4 Co sfalsyfikowane ❌

| ID | Co | Data | Konsekwencja |
|---|---|---|---|
| **M9.1 (power)** | `g_tt = -c²/ψ` | 2026-04-25 | β_PPN = 4 vs 1 (3·10⁴σ); pivot to M9.1' |
| **M9.1' (no-rescue)** | Power-law p≠0 nie ratuje β=1 | 2026-04-25 | Pivot to M9.1'' hyperbolic |
| **OP-6 M1-B** | Minimal bilinear → α=2 | 2026-04-24 | Pivot to H_GL axiom |
| **OP-6 M2-c** | "H₃ bilinear" α=1 | 2026-04-24 | Real source = H_GL not minimal |
| **disformal_waveform.py** | Disformal mechanism | pre-2026-04 | 18 rzędów falsyfikacji |
| **Galaxy scaling derivation** | TGP-native ν(y) function | 2026-04-19 | Phenomenological only |

---

## 2. Krytyczne sprzeczności do natychmiastowego naprawienia

Z 5 sond zebrałem 30+ unikalnych sprzeczności. Poniżej priorytety:

### 2.1 P0 (CRITICAL — must fix in 24h)

**[P0-1] §2 paperu — stale σ_ab formula**
- Plik: `tgp-core-paper/paper/tgp_core.tex` linie 295-311
- Stare: `σ_ab(x) ∝ ⟨ŝ_i · ŝ_{i+ê_a}⟩_B^TF`
- Nowe: `σ_ab = K_ab − (1/3)δ_ab Tr(K), K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩_B`
- Footnote: "(OP-7 T2, 2026-04-25: validated as composite from H_Γ kinetic term, replacing nearest-neighbor postulate.)"

**[P0-2] OP7_setup.md linia 7 — stary σ_ab w "punkcie wyjścia ontologicznym"**
- Plik: `TGP_v1/research/op7/OP7_setup.md`
- Stare: `σ_ab ∝ ⟨ŝ·ŝ_{+â}⟩^TF`
- Nowe: `σ_ab = K_ab − (1/3)δ_ab Tr(K)` (post-T2 supersession)

**[P0-3] TGP_FOUNDATIONS.md linie 64-71 — exponential metric jako "równoważna"**
- Stare: "Forma eksponencjalna (równoważna do O(U)): g_tt = -c₀² e^(-2U)"
- Nowe: "Forma hiperboliczna (M9.1'', 2026-04-25): g_tt = -c₀² (4-3ψ)/ψ"
- Add: "Power-form `g_tt = -c²/ψ` was falsified in M9.1 T3 (2026-04-25; β_PPN=4 vs obs 1)."

**[P0-4] TGP_v1/README.md linia 28 — stale exponential**
- Stare: "algebraic metric g_ij = e^{2U/c₀²}δ_ij, ... All 40 quantitative predictions"
- Nowe: "hyperbolic metric (M9.1'') ... 21 confirmed *(GW170817 conditional on OP-7 T3-T5)*"

**[P0-5] Skrypty hardcoded exponential metric**
- `tooling/scripts/cosmology/tgp_lensing_formal.py` (multiple lines)
- `tooling/scripts/gravity/ppn_3pn_corrections.py`
- `tooling/scripts/ex201_antipodal_metric_derivation.py`, `ex206_metric_hypothesis_necessity.py`
- Akcja: header DEPRECATED + reference to M9.1'' hyperbolic form

### 2.2 P1 (IMPORTANT — fix in 1 week)

- **[P1-1]** Paper §7 OP-7 description should reflect T1+T2 closure: add note "Kinematic half closed 2026-04-25; dynamical (T3-T5) open."
- **[P1-2]** Paper Abstract line 76-77 "luminal GW consistent with GW170817" — add footnote redirecting to cor:cGW (scalar sector caveat).
- **[P1-3]** README "GW prediction: scalar breathing mode (3rd polarization)" → "(1 polarization, scalar-only until OP-7 closure)".
- **[P1-4]** "All 10 PPN = GR" wszędzie → "γ=β=1 (1PN); 8 pozostałych pending OP-7".
- **[P1-5]** TGP_STATUS_2026-04-18.md → utworzyć TGP_STATUS_2026-04-25.md uwzględniający M9.1'' P1-P4 + OP-7 T1+T2.
- **[P1-6]** Globalny grep `c_GW = c` → wszędzie dodać "[scalar sector only — tensor pending OP-7]".
- **[P1-7]** Date stamps causality: README "v2 pivot 2026-04-24 → M9.1'' 2026-04-25 → OP-7 T1+T2 2026-04-25".

### 2.3 P2 (HOUSEKEEPING — fix in 2 weeks)

- **[P2-1]** Deprecate scripts: 6 plików (ex201, ex206, a2_metric_consistency, disformal_waveform, cmb_tensor_perturbations, ppn_3pn_corrections).
- **[P2-2]** Rewrite scripts: 12-15 plików (two_body_Veff, bh_*, profiles/*, dirac_tgp, tensor_from_substrate, gw_breathing_mode, p73_perturbations_CMB...).
- **[P2-3]** Patch tracking table w KNOWN_ISSUES.md: P1.1-P1.6 review patches z commit hashes.

---

## 3. Plan domknięcia OP-7 (krytyczna ścieżka)

To jest **rdzeń wszystkiego** — od domknięcia T3-T6 zależy czy TGP zostanie sfalsyfikowany na poziomie LIGO/Virgo czy zamknie się jako pełna teoria.

### 3.1 OP-7 T3 — Dynamika σ_ab (najtrudniejsze)

**Cel:** Wyprowadzić wariacyjnie z S_TGP[ŝ]:
```
□σ_ab + m_σ² σ_ab = -ξ T_ab^TT[matter]
```

**Główne wyzwanie:** Skąd bierze się masa `m_σ`? Trzy hipotezy do przetestowania:

**Hipoteza A (mean-field):** `m_σ² ~ J ⟨ŝ²⟩ = J·v₀²`
- Plus: naturalne z H_Γ self-consistent (J-term w GL-bond)
- Minus: m_σ ~ Φ₀ (Planck-scale), tłumi GW poniżej obserwabilnej amplitudy
- Test: numerycznie sprawdzić skalę `m_σ/Φ₀` z lattice MC

**Hipoteza B (kompresybilna fluktuacja):** `m_σ² ~ U''(Φ₀) - σ-mode-renormalization`
- Plus: skala odpowiada "potential mass" (μ²=2 dla TGP), fenomenologicznie OK
- Minus: wymaga renormalizacji σ-pętli, nietrywialna
- Test: 1-loop calc w sympy + numerycznie z H_Γ

**Hipoteza C (massless)**: `m_σ = 0`
- Plus: 1/r tail, długi zasięg, GR-like radiation pattern
- Minus: konstrukt graviton-like, ale skoro σ_ab jest composite — może być OK
- Test: sprawdzić Z₂ + ghost-free constraint

**Plan T3:**
1. **T3.1** (1 tydzień): wariacja S_TGP[ŝ] z odpowiednim coupling do macierzowego σ_ab; otrzymać EOM symbolicznie (sympy)
2. **T3.2** (1 tydzień): wyznaczyć m_σ z self-consistent MF (lattice + 1-loop)
3. **T3.3** (3 dni): ghost analysis (kinetic term sign, dispersion relation)
4. **T3.4** (3 dni): T_ab^TT source z stress-energy; hipoteza ξ ~ G_eff·Φ₀⁻²

**Werdykt T3:**
- POSITIVE: dynamika ghost-free, propagująca z c_σ (=c₀ jeśli m_σ=0) → idziemy do T4
- NEGATIVE: ghost lub tachyon → TGP wymaga dodatkowego mechanizmu (poważny pivot)

### 3.2 OP-7 T4 — Sprzężenie metryki Λ(ψ)

**Cel:** Postulat:
```
g_ij = h(ψ)·δ_ij + Λ(ψ)·σ_ij
```
zachowujący:
- `Λ(ψ_eq=1) = const` (definiuje skalę GW)
- `g_eff` ghost-free (no kinetic instabilities)
- `c_GW = c₀` w próżni (=> Λ niezależne od ψ jeśli σ_ab ma c₀ dispersion z T3)
- M9.1'' redukcja przy σ=0

**Naturalne kandydaty Λ(ψ):**
1. `Λ(ψ) = ψ` (linear, konformalny do trace)
2. `Λ(ψ) = ψ·(4-3ψ)` (matching M9.1'' algebra)
3. `Λ(ψ) = 1` (constant, najprostsze)
4. `Λ(ψ) = ψ/(4-3ψ)` (matching h(ψ))

**Plan T4:**
1. **T4.1**: wyznaczyć Λ(ψ) z requirement c_GW=c₀ + ghost-free (sympy linearization)
2. **T4.2**: PPN check: σ_ab = 0 → β=γ=1 unchanged
3. **T4.3**: Z₂ parity: σ_ab → σ_ab pod ŝ → -ŝ; metric niezmieniona

### 3.3 OP-7 T5 — Formuła kwadrupolowa

**Cel:** Obliczyć amplitudę `h+, h× ∝ Q̈_ij/r` dla binary inspiral; matching do GW150914.

**Plan T5:**
1. **T5.1**: T3 EOM + plane-wave w próżni → 2 d.o.f. propagujące (z OP-7 T1 wiemy że są)
2. **T5.2**: Source term T_ab^TT[binary] (point-mass approximation)
3. **T5.3**: Greens function dla □σ_ab + m_σ² σ_ab; far-field amplitude
4. **T5.4**: matching ξ_eff (jednowarstwowy parametr) do GW150914 strain ~1e-21

**Werdykt T5:**
- POSITIVE jeśli ξ_eff jest **physically reasonable** (= G_N·Φ₀² lub podobne); ALSO jeśli amplituda spłaszcza obserwowane wartości GW150914 (±20%)
- NEGATIVE jeśli ξ_eff jest absurdalnie różny od G — TGP fenomenologicznie nie pasuje

### 3.4 OP-7 T6 — Konsystencja

**Cel:** Pełne audyty:
- PPN: γ=β=1 unchanged + 8 pozostałych (α₁₋₃, ξ, ζ₁₋₄) — w M9.1'' static σ_ab=0 więc trivially zachowane, ALE moving matter (M9.2) może zaburzyć
- c_GW = c₀ exactly w próżni (z T4 Λ)
- Ghost-free w pełnym kinetic sector
- Z₂ parity: cała struktura σ_ab niezmieniona pod `ŝ → -ŝ`
- Stability under perturbations (linearized solutions bounded)

**Werdykt OP-7 (cały):**
- T1+T2+T3+T4+T5+T6 wszystko PASS → OP-7 CLOSED-POSITIVE → tgp_core.tex v3 → submission to journal
- Jakikolwiek FAIL → identyfikujemy konkretny mechanizm i decydujemy: pivot (akceptowalny jeśli zachowuje single-Φ Z₂) czy falsyfikacja TGP

---

## 4. Plan paralelnych programów (niezależne od OP-7)

Dopóki OP-7 idzie, **inne programy mogą iść równolegle**:

### 4.1 M9.2 — Momentum / Lenz back-reaction (HIGH PRIORITY)

**Cel:** Sprawdzić czy ruchome źródła w M9.1'' produkują 8 pozostałych PPN parameters zgodnie z GR.

**Plan:**
- **M9.2 P1** (1 tydzień): Lenz mechanism (substrate inertia from ŝ-redistribution); sympy + numerical
- **M9.2 P2** (1 tydzień): PPN α₁, α₂, α₃ z moving mass distribution
- **M9.2 P3** (3 dni): preferred-frame, preferred-location effects (currently obs limit ~10⁻⁴)

**Krytyczne:** M9.2 wymaga spójnego z M9.1'' static rozwiązania; może wymusić Λ(ψ) z OP-7 T4 niezależnym kanałem.

### 4.2 OP-EHT — Strong-field (lower priority, but visibility)

**Cel:** Numerical solution full nonlinear Φ-EOM dla static spherical (no weak-field expansion).

**Plan:**
- **EHT-P1**: scipy.solve_bvp na M9.1'' Φ-EOM, photon ring computation
- **EHT-P2**: shadow size dla M87*, Sgr A* + comparison z observation

**Status:** Może iść równolegle do OP-7 T3-T6, ale M87* shadow ma rozmiary ±5% co jest miękki test.

### 4.3 OP-1 / OP-4 — Potential uniqueness (LOW priority, theoretical)

**Plan:**
- **OP-1 M2** (gdy OP-7 closed): wariacyjna derivation U(φ) = β/3 φ³ - γ/4 φ⁴ z H_Γ poprzez HS-transform i 1-loop integration
- **OP-4 M4**: matching γ = g₀^ε do RG fixed point (M3-c sygnalizował obstrukcję, ale H_GL pivot może uratować)

**Status:** Niech czeka — jeśli OP-7 nie zamknie się, OP-1 jest nieistotne.

### 4.4 OP-8 — U(1) gauge emergence (DEFERRED)

**Status:** Zarezerwowany do v2 papera (v1 zamknięty na grawitacji).

---

## 5. Plan housekeepingu kodu (równolegle z teorią)

### 5.1 Faza A — Deprecation (1 tydzień, niski wysiłek)

Stworzyć `tooling/scripts/deprecated/` i przenieść:
- `ex201_antipodal_metric_derivation.py`
- `ex206_metric_hypothesis_necessity.py`
- `a2_metric_consistency.py`
- `gw/disformal_waveform.py`
- `gw/cmb_tensor_perturbations.py`
- `gravity/ppn_3pn_corrections.py`
- `cosmology/p73_perturbations_CMB.py`

Każdy plik dostaje header:
```python
# DEPRECATED 2026-04-25
# Pre-pivot script. Replaced by:
#   - M9.1'' hyperbolic metric (research/op-newton-momentum/M9_1_pp_*.md)
#   - OP-7 T1+T2 (research/op7/OP7_T{1,2}_results.md)
# Kept for reference; do NOT use for new analyses.
```

### 5.2 Faza B — Rewrite (3-4 tygodnie)

Priorytet (każdy = 2-3 dni pracy):
1. `gravity/two_body_Veff.py` (M9.1'' hyperbolic V_eff)
2. `bh_isco_raytrace_tgp.py`, `bh_shadow_tgp_vs_gr.py` (geodesics on hyperbolic metric)
3. `profiles/single_source_profile.py` + 3 inne (numerical hyperbolic ODE)
4. `substrate/tensor_from_substrate.py` → przepisać jako `gradient_strain_substrate.py` z OP-7 T2 formula
5. `gw/gw_breathing_mode.py` (m_sp(ψ) dla hyperbolic; rozszerzyć o predictions T3 σ_ab)
6. `cosmology/tgp_perturbations_formal.py` (dodać σ_ab perturbation section)

### 5.3 Faza C — Update docstrings tylko (kilka godzin)

- `tooling/scripts/substrate/tensor_from_substrate.py`: docstring → reference OP-7 T2 + flag pre-pivot
- Wszystkie scripty wspomniane w §7 paperu: dodać citation block z M9.1'' / OP-7 setup

---

## 6. Plan integracji do paperu (tgp_core.tex)

### 6.1 Faza P1 (immediate, ≤1 dzień)

- §2 σ_ab formula update (P0-1)
- §7 OP-7 box: dodać note T1+T2 closed
- Abstract footnote: tensor caveat (P1-2)

### 6.2 Faza P2 (po OP-7 T3-T4)

- §6 (NEW): "Tensor sector dynamics" — rezultat T3 (EOM + m_σ + ghost analysis)
- §6.1: σ_ab definition (T2)
- §6.2: dynamics (T3)
- §6.3: metric coupling Λ(ψ) (T4)

### 6.3 Faza P3 (po OP-7 T5-T6, ostateczna submission)

- §7 update F3 (luminal GW): tensor sector closed
- §8 conclusion: GW170817 unconditional pass
- New corollary `cor:GW-amplitude`: kwadrupol formula z ξ_eff
- Przepisać abstract: drop "OP-7 pending" caveats

### 6.4 Wersja v3 do Zenodo

Kiedy:
- OP-7 T3-T6 PASS
- M9.2 PPN spójne z α₁=α₂=α₃=0 (lub explicit deviations)
- Critical sprzeczności (P0/P1) wszystkie patched

Co:
- Bump tgp_core.tex do v3
- New DOI
- Update KNOWN_ISSUES.md (zamknąć większość P-patches)

---

## 7. Plan TGP_FOUNDATIONS update

### 7.1 Critical edits (24h)

1. Zaktualizować §3 metryka: hyperbolic form jako kanonical, exponential = pre-pivot
2. Dodać §3.1 (NEW): "M9.1'' Pivot — od boxed do hyperbolicznej"
3. Update §2 (Hierarchia): dodać row σ_ab z OP-7 T1/T2 status
4. Update §7 (M3-M8): explicit "ARCHIVED, pre-M9 program"

### 7.2 New axiom proposal (po brainstorm)

Rozważyć wprowadzenie `ax:metric-from-potential`:
```
Effektywna metryka jest określona przez znormalizowaną gęstość potencjału:
  f(ψ) = V(Φ)/Φ⁴ z f(1) = 1
Dla V(Φ) = (γ/12)·Φ₀²·ψ³·(4-3ψ):
  f(ψ) = (4-3ψ)/ψ                    [hyperbolic form]
```
**Status:** postulat z potrójną substratową motywacją (P2-C/D/E), bez jednostopniowego wyprowadzenia z S_TGP. **To jest naturalny kandydat do "deeper principle" derivation.**

---

## 8. BRAINSTORM — co naprawdę wymaga myślenia teoretycznego

Synteza nie wystarczy. Tu są **kreatywne otwarte pytania** które gating closure nie z powodu pracy mechanicznej, ale dlatego że odkrywamy gdzie TGP może wymagać prawdziwych nowych insightów.

### 8.1 σ_ab dynamika — czy mamy ducha?

**Centralna obawa:** σ_ab jest gradient strain `⟨(∂_a ŝ)(∂_b ŝ)⟩`. To jest **bilinearny** w `ŝ`, więc kwadratowy w fluktuacjach. Wariacyjna pochodna względem σ_ab daje EOM, ale czy ta EOM jest **kanoniczna** czy **non-canonical (drugiego rzędu w czasie z wielokrotnymi rozwiązaniami)**?

**Hipoteza-kontra:** Jeśli σ_ab jest composite, jego "EOM" to faktycznie **łańcuchowa równość operatorowa** (correlator equation), nie classical PDE. Wtedy dynamiczne zachowanie σ_ab bierze się **automatycznie** z dynamiki bazowego ŝ — nie trzeba T3 wyprowadzać niezależnej PDE!

**Sugestia:** w T3 sprawdzić **dwie ścieżki**:
- A: traktować σ_ab jak quasi-fundamental (wprowadzić Lagrangian L[σ_ab] explicitly) — może produkować ghost
- B: dziedziczyć dynamikę z `□ŝ + V'(ŝ) = 0` i uśredniać `⟨(∂_a ŝ)(∂_b ŝ)⟩` (effective Boltzmann-like equation)

**Path B jest bardziej "single-substrate-faithful"** — żaden nowy d.o.f., żadnych ghostów z konstrukcji, bo wszystko ma kinetyczny term ŝ-pochodny.

### 8.2 M9.1'' — dlaczego (4-3ψ)/ψ a nie inna funkcja?

**Jeszcze nie znamy "głębszej zasady"** wyboru hyperbolic form. Triple convergence (P2-C/D/E) jest podejrzanie zbieżne — to często sygnał że **istnieje jedno głębokie principle** z którego wszystkie trzy wynikają.

**Kandydaty:**
1. **Information geometry**: f(ψ) = V(ψ)/ψ⁴ to jest naturalna metryka Fishera dla rozkładu p(ŝ) ~ exp(-V(ŝ))?
2. **Holographic**: ψ⁴ to volume-form dla niekomutatywnych współrzędnych (ŝ-substrate), V to energia, V/ψ⁴ to densja energii per node?
3. **Substrate hydrodynamics**: f(ψ) jest "speed of sound squared" dla substrate density wave; (4-3ψ)/ψ to algebraiczna postać "compression of substrate"?

**Test:** wyprowadzić f(ψ) = V/Φ⁴ z **jednego** axiomatu (nie trzech). Naturalny kandydat: minimization of free energy on coarse-grained substrate.

### 8.3 m_σ — Yukawa hint?

Częstotliwość masy substratu w fizyce TGP to **`μ = √2`** w units Φ₀ (kink mass, paper §B). Czy `m_σ ~ μ`?

**Pomysł:** σ_ab to **kompozyt dwóch derivatives ŝ** → kompozytowa masa to **suma mas składników** (analogia do mezonów: m_π² ~ m_q + Λ_QCD²).
- Jeśli `m_ŝ ~ μ = √2`, to `m_σ ~ 2μ = 2√2` (units Φ₀).
- W konwencjach SI: `m_σ ~ √2·m_ŝ ~ √2·(Planck mass)/(Φ₀-norm)`
- Dla TGP-galaktycznego skalarna `Φ₀ ~ Hubble-1` to daje `m_σ ~ µeV` — w obszarze ULDM!

**Spekulatywne:** σ_ab z masą µeV-scale? Wtedy GW poniżej kHz są massless-like, ale powyżej kHz pojawiają się dispersyjne efekty. **Predykcja testowalna LIGO future.**

### 8.4 Λ(ψ) — czy ją znamy z first principles?

**Pytanie:** czy Λ(ψ) jest wyznaczane przez wymóg `c_GW = c₀ exactly` w próżni? Sprawdzenie dimensional analysis:
- σ_ab ma [Φ²/L²] (gradient strain)
- g_ij to bezwymiarowe → Λ(ψ)·σ_ij musi być bezwymiarowe → Λ ~ [L²/Φ²] = 1/Φ²
- Naturalna skala: `Λ(ψ) = 1/Φ₀² · λ̃(ψ)` gdzie λ̃ bezwymiarowe

**Hipoteza:** `Λ(ψ) = Λ₀/ψ²` (matching `1/Φ²`); wtedy g_ij = h(ψ)δ_ij + Λ₀·σ_ij/ψ²

**Test:** sprawdzić czy ta postać daje c_GW=c₀ ghost-free w T4.

### 8.5 Czy single-substrate przetrwa T3?

**Najgłębsza obawa.** Jeśli T3 wymaga:
- niezależnego d.o.f. dla σ_ab (= TGP staje się scalar-tensor)
- dodatkowego pola źródłowego (= TGP nie jest single-Φ)
- ghosta którego nie umiemy usunąć

— to TGP jako "single-Φ Z₂ scalar emergent gravity" jest **zfalszyfikowany** w obecnej formie.

**Path forward jeśli to się stanie:**
1. **Przyjąć multi-component substrate** (ŝ₁, ŝ₂, ...) — to jest porażka aksjomatu, ale teoria może być uratowana.
2. **Przyjąć non-trivial spacetime topology** (foliated substrate) — dramatyczny pivot, ale potencjalnie elegancki.
3. **Zaakceptować TGP jako "approximate" theory** — nie obejmuje tensor sector, działa tylko dla scalars (cosmology, particle, Φ-static).

**Moja ocena:** prawdopodobieństwo że T3 PASS ze single-Φ ~ 60-70%, z drugim pivotem (Λ(ψ) form-fitted) ~ 20%, falsyfikacja ~ 10-15%. To jest moment prawdy.

### 8.6 Co jeśli σ_ab nie ma masy?

Jeśli `m_σ = 0`, σ_ab propaguje na c₀, jest spin-2 composite, **dokładnie jak emergent graviton** w teoriach typu Volovik / Verlinde.

**Implikacja:** mamy emergent spin-2 mode z properties identicznych z gravitonem — co prawda jako collective mode, ale fenomenologicznie to to samo. **TGP wtedy zbliża się do "induced gravity" Sakharova** w sensie precyzyjnym.

**Pomysł:** w T6 sprawdzić czy emergent spin-2 mode ma 2 polarizations (z OP-7 T2 wiemy że 5 d.o.f. → 2 po TT projection w propagacji), a kwadrupol formula z T5 daje GR amplitude do leading order. Jeśli tak — TGP jest **MIKROSKOPICZNĄ INSTANCJĄ** scenariusza Sakharov-Verlinde, czyli realizuje filozoficzną wizję która od dekad pozostaje "in principle" ale bez konkretnej teorii.

To byłby **ogromny rezultat**.

### 8.7 GW170817 reasoning

**Obecna sytuacja:**
- M9.1'' single-Φ → tylko breathing mode (T1)
- LIGO mierzy h+, h× (TT modes) z 5% bound
- Bez σ_ab dynamics: TGP scalar-only signal jest **zerowy** w TT (różne polarisations!), więc nie ma "5% mismatch" — jest **brak signal**.

**Niuans:** w obecności matter binary inspiral, scalar breathing mode też produkuje detectable strain (różne template). LIGO **też nie wykryło** tego — co znaczy że amplituda breathing mode jest <5% TT amplitude.

**Predykcja TGP single-Φ (bez σ_ab):**
- Breathing mode z amplitude ~ G·M/r (kwadrupol w Φ-projection)
- LIGO antenna pattern dla breathing = `(1/2) sin²θ` (vs. F+ = `(1+cos²θ)/2 cos2φ`)
- Unrejected → daje upper bound na TGP scalar coupling

**Z σ_ab dynamics (T3 PASS):** TT modes pojawiają się jako **second-order** w substrate fluctuations (σ_ab ~ ⟨(∂ŝ)²⟩). Amplituda TT jest **suppressed** względem GR przez ~ξ_eff/G_N. Czy ξ_eff = G_N daje normalne amplitudy GR-like? Test w T5.

**Rough estimate:** jeśli ξ_eff ~ G_N (coupling normalization), TGP TT amplitude ~ GR amplitude do leading order, second-order corrections ~ (v/c)² (już w GR znanej Skrócić jako 2.5PN). Spójność pełna. **Test:** T5 numerycznie.

### 8.8 Czy mamy "smoking gun" dla TGP vs GR?

Najbardziej krytyczna prediction TGP **niezależna od OP-7**:
- 2PN deviation `|Δg_tt| = (5/6)U³` w solar system

**Eksperymentalna realność:**
- Cassini: γ-1 < 2.3·10⁻⁵ przy U_Sun ~ 10⁻⁸ → 2PN signal `~ U³ ~ 10⁻²⁴` — **nieobserwabilne**
- Mercury precesja: 2PN ~ 10⁻¹² rad/orbit — nieobserwabilne
- BNS inspiral GW170817: 2PN phase ~ ±0.5 rad observed, TGP add `5/6·U³ ~ 0.01 rad` — **w tle** GR predykcji
- LIGO third-generation (Cosmic Explorer ~2030): może rozróżnić ±0.001 rad → **może wykryć TGP**

**Strong-field smoking gun:**
- Photon ring around BH: M9.1'' przewiduje `r_ring/r_g` różne od Kerr o ~10% — **EHT obserwable**
- Ringdown spectrum: M9.1'' nieliniowy → odmienne quasi-normal modes vs Kerr — **LIGO BBH obserwable**
- Both wymagają OP-EHT + nonlinear M9.1'' P3 strong-field closure

### 8.9 Cosmologia — zgłoszony "TGP-native dark energy"?

W M9.1'' V(Φ) = (γ/12)Φ₀²ψ³(4-3ψ); w próżni V(1) = (γ/12)Φ₀² ≠ 0 → **vacuum energy**.

**Hipoteza:** TGP automatycznie predicts cosmological constant Λ ~ V(1) = (γ/12)Φ₀². 
- Jeśli Φ₀ ~ Planck, to Λ ~ M_P² → 10¹²² × obserwowane (vacuum catastrophe).
- Jeśli Φ₀ ~ µeV (galaktyczna), to Λ ~ µeV² ~ obserwowane Λ_obs!

**To jest bardzo ciekawe.** Skala Φ₀ wymaga wyznaczenia z OP-3 (`a_Γ = 1/Φ₀`). Jeśli Φ₀ ~ Hubble^(-1) (µeV), to Λ_TGP ~ Λ_obs **automatycznie** → cosmological constant problem solved by TGP!

**Plan:** dodać do OP-3 milestone `Λ_TGP from Φ₀ scale` — to może być **najgłośniejsze prediction TGP** poza GR closure.

### 8.10 EHT — quick win?

Z OP-7 T1 wiemy że σ_ab=0 dla statycznego sferycznego BH (izotropia → ⟨(∂_a ŝ)(∂_b ŝ)⟩ ∝ δ_ab → traceless = 0). Więc **photon ring jest identyczny w M9.1'' single-Φ i M9.1''+σ_ab dla static sphere**.

**Implikacja:** OP-EHT można zrobić bez czekania na T3-T6! Tylko M9.1'' Φ-EOM nonlinear jest potrzebne.

**Plan EHT-quick:**
1. scipy.solve_bvp na M9.1'' Φ-EOM dla static sphere
2. Photon geodesic equation w hyperbolic metric → ring radius `r_ph`
3. Compare z Kerr `r_ph = 3M_BH` → predict M9.1'' `r_ph(M_BH, Φ₀)`
4. Match z M87* observed shadow size

**Czas:** 2-3 dni przy obecnej infrastrukturze. **Może dać wczesny smoking gun.**

---

## 9. Roadmap — kwartał po kwartale

### 9.1 2026 Q2 (kwiecień-czerwiec) — DOMKNIĘCIE OP-7 + EHT-quick

| Tydzień | Cel | Output |
|---|---|---|
| **TYG 17 (24-30 IV)** | P0/P1 patches, OP-7 T3.1 setup | clean repo + sympy EOM draft |
| **TYG 18 (1-7 V)** | OP-7 T3.1-T3.4 (full T3) | T3 verdict (POSITIVE/NEG); ghost analysis |
| **TYG 19 (8-14 V)** | EHT-quick (photon ring static); M9.2 P1 | M87* prediction; PPN α₁₋₃ |
| **TYG 20 (15-21 V)** | OP-7 T4 Λ(ψ) | T4 verdict; ghost-free metric |
| **TYG 21 (22-28 V)** | OP-7 T5 quadrupole | h+, h× formula; ξ_eff matching |
| **TYG 22-26 (czerwiec)** | OP-7 T6, M9.2 P2-P3, paper §6 draft | OP-7 closed; M9.2 closed; v3 paper draft |

**Q2 Output:** OP-7 fully closed (positive lub identified path forward), v3 paper ready do circulation.

### 9.2 2026 Q3 (lipiec-wrzesień) — OBSERVATIONAL CLOSURE + v3 SUBMISSION

| Miesiąc | Cel |
|---|---|
| **VII** | M9.3 GW radiation + GW170817 unconditional + GW150914 quantitative |
| **VIII** | EHT full nonlinear; OP-EHT closure |
| **IX** | tgp_core.tex v3 → arXiv + Zenodo; community circulation |

**Q3 Output:** TGP fully closed at first-order observational predictions. Paper submitted.

### 9.3 2026 Q4 (październik-grudzień) — RESPONSE + EXTENSIONS

- Q4 dla peer review responses + corrections
- Start OP-1, OP-3, OP-4 (potential uniqueness, scale identification)
- Start OP-8 (gauge emergence) — zaczyna się v2 paper

### 9.4 2027+ — v2 PROGRAM

- OP-8 U(1) gauge emergence
- OP-10 quantum theory of Φ
- Extended particle sector (quark/Yukawa derivations)
- Cosmology paper (Λ from Φ₀, BBN, CMB)

---

## 10. Decyzje wymagające zatwierdzenia użytkownika

Zanim ruszę z czymkolwiek, potrzebuję twojego sign-off na:

### 10.1 Architektura — kolejność egzekucji

**Opcja A (SEKWENCYJNIE)**: P0 patches → OP-7 T3 → ... → T6 → paper v3
- Plus: minimalne ryzyko regresji
- Minus: ~6-8 tygodni do v3

**Opcja B (RÓWNOLEGLE)**: P0 patches + OP-7 T3 + EHT-quick + M9.2 P1 razem
- Plus: maksymalna prędkość, do 4 tygodni do v3
- Minus: trudniej zarządzać contextem; jeśli OP-7 T3 FAIL, trzeba undo M9.2/EHT

**Moja rekomendacja:** Opcja B z safety hatch — startujemy P0/P1 (1 dzień) + EHT-quick (3 dni) + OP-7 T3 (1 tydzień) **równolegle**. Jeśli OP-7 T3 NEGATIVE, zatrzymujemy M9.2 i wracamy do brainstorm 8.5.

### 10.2 Backup theoretyczny

**Co jeśli OP-7 T3 odkryje ghost?** Plan B:
- (B1) Λ(ψ) ad-hoc form-fit (akceptujemy że nie deriwujemy z first principles)
- (B2) Restrict TGP to scalar-only sector (emergent gravity opisuje tylko Φ-projection, GW polarisations są **outside** TGP scope)
- (B3) Multi-substrate pivot (TGP staje się scalar-tensor) — utrata aksjomatu

**Pytanie do ciebie:** czy w razie B1 chciałbyś jechać dalej z form-fitted Λ, czy raczej zatrzymać i opublikować v3 jako "scalar sector closure" (bez tensor)?

### 10.3 Priorytety brainstorm

Z 10 punktów brainstorm w §8, które chcesz teraz pogłębiać?
- **§8.5 + §8.9** (single-substrate survival + cosmological constant) → najgłębsze theoretical work
- **§8.10** (EHT quick) → najszybszy observational win
- **§8.3** (m_σ ~ µeV) → potencjalnie ULDM connection
- **§8.6** (TGP jako mikroskopowy Sakharov) → potencjalnie najgłośniejsze claim

---

## 11. Bottom line

TGP jest **6-8 tygodni od pełnego domknięcia** lub **identyfikacji konkretnego pivotu**. Główne ryzyka:

1. **OP-7 T3 ghost** — 10-20% prawdopodobieństwo; mamy plan B
2. **Λ(ψ) trudne form-fit** — 30% prawdopodobieństwo; zarządzamy poprzez T4 ghost analysis
3. **GW170817 amplitude mismatch** — 20% prawdopodobieństwo w T5; oznaczałby częściową falsyfikację

Reszta to **pracowite kontynuowanie** istniejącej trajektorii. Wszystkie krytyczne prywatne sprzeczności zostały zidentyfikowane (~30) i mają plan naprawy.

**Naprzód:** zatwierdź §10.1 + §10.2 + §10.3, a startuję natychmiast.

---

## 12. Cross-references do tej syntezy

- [[OP7_setup.md]] (master plan T1-T6)
- [[OP7_T1_results.md]] (no-tensor 7/7)
- [[OP7_T2_results.md]] (σ_ab gradient strain 12/12)
- [[../op-newton-momentum/M9_1_pp_P3_results.md]] (P3 GW170817 conditional)
- [[../../../tgp-core-paper/paper/tgp_core.tex]] (paper)
- [[../../../tgp-core-paper/KNOWN_ISSUES.md]] (gap inventory)
- [[../../TGP_FOUNDATIONS.md]] (axioms)

---

**Przeznaczenie tego dokumentu:** living strategy doc, aktualizowany po każdym milestone. Następna rewizja: po OP-7 T3 verdict.
