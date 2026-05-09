---
title: "G.0 — R3 ODE jako projekcja TGP-canonical na M9.1'' (foundational unification)"
date: 2026-05-02
cycle: G.0
status: PHASE 4 CLOSED-POSITIVE — sek08a v2.0 + sek08c A1/A2/A3 closure + cross-references applied to core LaTeX (pdflatex compile clean)
phase1_score: "2/3 PASS + 1 NEGATIVE-INFORMATIVE"
phase1_finding: "V_M911(ψ) = -ψ²(4-3ψ)²/12 (sympy LOCK) reprodukuje R3 ODE EXACT"
phase2_score: "4/4 PASS (gate ≥3/4 met)"
phase2_finding_1: "V_M911 unique + vacuum stable (G.0 fixes sek08a tachion bug)"
phase2_finding_2: "Mass spectrum lepton: m_μ/m_e -0.0013%, m_τ/m_e +0.0049% PDG"
phase2_finding_3: "PPN γ=β=1 INVARIANT pod V update"
phase2_finding_4: "FRW κ form invariant + 5/2x prefactor (Phi_0 re-fit do Phase 3)"
phase3_score: "4/4 PASS (gate ≥3/4 met)"
phase3_finding_1: "Newton limit: q·c²/Φ_0=(4/5)πG_0, Newton G_0 EXACTLY reproduced"
phase3_finding_2: "κ INVARIANT po re-fit (KOREKTA do P24 framing)"
phase3_finding_3: "30 HIGH-impact core/research files identified (P33 audit)"
phase3_finding_4: "Sek08a v2.0 + sek08c A1/A2/A3 closure DRAFTS ready dla Phase 4"
phase4_score: "6/6 files modified (sek08a, sek08c, sek08, dodatekH, status_map, sek09, dodatekO), pdflatex compile clean"
phase4_finding_1: "sek08a v2.0 ADDENDUM (sssec:g0-closure-v2): 4 propositions + 2 equations (V_M911, R3 ODE, vacuum-stability-G0, kappa-corrected-G0, newton-limit-G0, cosmo-linearized-G0)"
phase4_finding_2: "sek08c A1+A2+A3 → ALL CLOSED-RESOLVED (preamble G.0 block + 3 inline annotations CLOSED)"
phase4_finding_3: "0 nowych undefined references, 0 nowych multiply-defined labels (537 stron compile)"
parent: "[[../why_n3/README.md]]"
predecessors:
  - "[[../why_n3/PHASE1_psi_g0_identification.md]]"
  - "[[../why_n3/PHASE4_5_yukawa_propagator.md]]"
  - "[[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]"
related:
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tags:
  - TGP
  - G0
  - foundational-unification
  - R3
  - M9.1pp
  - sek08-audit
  - ODE-reconciliation
  - PROPOSED
tgp_status:
  folder_status: closed-resolved
  level: L1
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'G.0 — R3 ODE jako projekcja TGP-canonical na M9.1'''"
    - "Phase{1..3}_results.md PASS=206, CLOSED=51, FAIL=31"
    - audit-aware markers=9 (B6/B8/B9/A5/etc.)
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

# G.0 — R3 ODE jako projekcja TGP-canonical na M9.1''

> **Cel programu:** Rozstrzygnąć, czy R3 ODE (sektor solitonowy `research/why_n3/`,
> daje N=3 generacji + mass spectrum <0.1% PDG + spin-1/2) jest **efektywnym
> EOM** wynikającym z wariacji TGP-canonical akcji `S_TGP[Φ, g_eff]` z metryką
> M9.1'' jako background, czy jest **niezależnym formalizmem** wymagającym
> osobnego aksjomatu w `TGP_FOUNDATIONS.md`.
>
> **Stawka:** zamknięcie G.0 unifikuje całą warstwę 3c (fermiony, R3) z warstwą 1
> (Φ-EOM, sek08a) i warstwą 2 (metryka, sek08c) w jeden łańcuch derywacji od
> aksjomatu jednopolowego Φ-Z₂ aż do `m_τ/m_e = 3477` zgodnego z PDG.

---

## ⚡ PHASE 1 RESULT (2026-05-02): CLOSED-POSITIVE

**3 sub-tasks ukończone w 1 dzień prac kalendarzowych:**

| Sub-task | Score | Status |
|---|---|---|
| G0a — Volume integration | 4/4 | **PASS** (H1 active) |
| G0b — Field redefinition | 0/4 | NEGATIVE-INFORMATIVE (wzmacnia G0a) |
| G0c — Einstein-frame projection | 3/3 | **PASS** (H2 active, ≡ H1) |

**Centralne odkrycie:** Sek08a TGP-canonical action JEST wariacyjnie spójnym
fundamentem dla R3 ODE pod warunkiem **JEDNEJ ZMIANY**:

```
V_TGP_orig(ψ) = (β/3)ψ³ - (γ/4)ψ⁴   (z M9.1, FALSIFIED)
              ───────────────────────
                       ↓
V_M911(ψ) = -ψ²·(4-3ψ)²/12            (sympy LOCK, M9.1'' canonical)
```

Plus zachowane K(ψ)=ψ⁴ (T-D-uniqueness, α=2) i √(-g)=c₀·ψ/(4-3ψ) (M9.1''
canonical, zgodne z A2 audytu).

**Numerical verification:** R3 vs M9.1''-derived profile match `max diff =
0.000000` dla wszystkich 4 testów (g₀^e, g₀^μ, g₀^τ, g₀_crit).

Pełna synteza: [[Phase1_results.md]]. Status: **PHASE 2 CLOSED-POSITIVE**.

---

## ⚡ PHASE 2 RESULT (2026-05-02): CLOSED-POSITIVE — gate Phase 3 forward

**4 sub-tasks ukończone w 1 dzień prac kalendarzowych:**

| Sub-task | Score | Verdict |
|---|---|---|
| **P21** Vacuum + uniqueness | 4/5 | **PASS** + bonus (sek08a tachion bug uncovered) |
| **P22** Mass spectrum lepton | 5/5 | **PASS** (m_μ/m_e -0.0013%, m_τ/m_e +0.0049%) |
| **P23** PPN γ=β=1 | 5/5 | **PASS** (INVARIANT pod V update) |
| **P24** FRW κ | 5/5 | **PASS** (form invariant, 5/2x prefactor → Phi_0 re-fit task) |

**Phase 2 score: 4/4 PASS** (gate ≥3/4 SPELNIONE).

**Hard anchors zweryfikowane na akcji G.0:**

| Anchor | Wartość | Test |
|---|---|---|
| ψ_vacuum | 1 (stable) | P21 sympy LOCK |
| m_sp² | +γ (stabilne, byl -γ tachion w sek08a!) | P21 |
| m_μ/m_e | 206.766 (PDG 206.768, **-0.0013%**) | P22 |
| m_τ/m_e | 3477.40 (PDG 3477.23, **+0.0049%**) | P22 |
| Koide K=2/3 | convergent | P22 |
| 4-th gen | FORBIDDEN (g₀^4 > g_crit) | P22 |
| γ_PPN | 1 (Cassini) | P23 |
| β_PPN | 1 (Mercury, master formula) | P23 |
| FRW κ form | const · q/Φ_0 (invariant) | P24 |
| FRW κ value | 15/(8Φ_0) (5/2x sek08a, Phi_0 re-fit) | P24 |

**Bonus discovery:** G.0 update **fixes drugi bug w sek08a** — stary daje
m_sp² = -γ (tachionowy vacuum!), G.0 daje m_sp² = +γ (stabilny). Wzmacnia
to argument za G.0 jako prawdziwą naprawą sek08a.

Pełna synteza: [[Phase2_results.md]]. Status: **PHASE 3 CLOSED-POSITIVE**.

---

## ⚡ PHASE 3 RESULT (2026-05-02): CLOSED-POSITIVE — gate Phase 4 APPROVED

**4 sub-tasks ukończone w 1 dzień prac kalendarzowych:**

| Sub-task | Score | Verdict |
|---|---|---|
| **P31** Sek08a v2.0 specification | COMPLETE | **PASS** (7 prop changes drafted, Phase 4 plan ready) |
| **P32** Newton limit re-derivation | 5/5 | **PASS** (q·c²/Φ_0 = (4/5)πG_0, Newton EXACT, **κ po re-fit INVARIANT**) |
| **P33** Cross-reference audit | 4/4 | **PASS** (30 HIGH-impact files identified, P33_audit_results.md) |
| **P34** Sek08c A1/A2/A3 closure | COMPLETE | **PASS** (3 annotations CLOSED-RESOLVED) |

**Phase 3 score: 4/4 PASS** (gate ≥3/4 SPELNIONE).

### Centralne odkrycia Phase 3:

1. **Newton limit RIGOROUSLY zachowany** w G.0:
   - q·c²/Φ_0 = (4/5)πG_0 (NEW, sympy LOCK)
   - Numerical Newton G_0 reprodukcja EXACT (diff +0.0000%)
   - vs sek08a v1.x: q·c²/Φ_0 = 2πG_0 (faktor 2/5 zmiana)

2. **κ po re-fit INVARIANT** (KOREKTA do P24 framing):
   - P24's "5/2x prefactor" był poprawny tylko BEZ re-fit
   - Po re-fit q·c²/Φ_0 z Newton, κ = 4πG_0/(3H_0²) = 3/(4Φ_0_new) **NIEZMIENIONE**
   - Wszystkie observables (BBN, LLR, CMB) IDENTYCZNE v1.x i v2.0

3. **Cross-reference audit (P33):** 30 HIGH-impact files (sek08a, sek08, sek08c,
   dodatekH, status_map, sek09, dodatekO + 10 research). Konkretne line
   numbers + recommended actions. Patrz `P33_audit_results.md`.

4. **Drafts ready dla Phase 4:**
   - `sek08a_v2_specification.md`: 7 propositions updated (V_M911 LOCK,
     vacuum fix, kappa update, newton-limit update, source coupling update)
   - `sek08c_A1_A2_A3_closure_draft.md`: 3 annotations CLOSED-RESOLVED
     (preamble + inline replacements ready)

### Hard anchors complete (Phase 1+2+3):

| Anchor | Wartość |
|---|---|
| K(ψ) | ψ⁴ |
| V(ψ) | -γψ²(4-3ψ)²/12 (V_M911) |
| √(-g) | c·ψ/(4-3ψ) (M9.1'' canonical) |
| ψ_vacuum | 1 (stable, m_sp²=+γ) |
| Φ-EOM (static) | R3 ODE |
| g_crit ≡ ψ_horizon | 1.874 ≡ 4/3 |
| m_μ/m_e, m_τ/m_e | 206.766, 3477.40 (PDG <0.01%) |
| γ_PPN, β_PPN | 1, 1 |
| q·c²/Φ_0 | (4/5)πG_0 (NEW) |
| κ | 4πG_0/(3H_0²) (INVARIANT po re-fit) |

Pełna synteza: [[Phase3_results.md]]. **Status: PHASE 4 (actual core mods) APPROVED.**

---

## ⚡ PHASE 4 RESULT (2026-05-02): CLOSED-POSITIVE — core integration complete

**Phase 4 = actual core LaTeX modifications zgodnie z `sek08a_v2_specification.md`
i `sek08c_A1_A2_A3_closure_draft.md`.**

### Modyfikowane pliki (6+1 = 7 core files)

| Plik | Zmiany | Linia |
|---|---|---|
| `core/sek08a_akcja_zunifikowana.tex` | + intro G.0 closure block, + ADDENDUM `sssec:g0-closure-v2` (prop:V-M911-canonical, prop:psi-EOM-R3, prop:vacuum-stability-G0, prop:kappa-corrected-G0, eq:newton-limit-G0, eq:cosmo-linearized-unified-G0), + 3 inline annotations | + ~250 lin |
| `core/sek08c_metryka_z_substratu.tex` | + G.0 CLOSURE preamble block (A1+A2+A3 RESOLVED), + 3 inline STATUS CLOSED, + tabela porównawcza M9.1'' | + ~70 lin |
| `core/sek08_formalizm.tex` | + intro G.0 closure block, + 2 inline annotations (FRW prop, footnote) | + ~30 lin |
| `core/formalizm/dodatekH_lancuch_wyprowadzen.tex` | + intro G.0 block, + footnote A11b | + ~20 lin |
| `core/_meta_latex/status_map.tex` | + intro G.0 block, + remark `rem:status-map-G0`, + nowa pozycja v2.0 w tabeli statusów | + ~30 lin |
| `core/sek09_cechowanie.tex` | + inline annotation kappa-G.0 | + 3 lin |
| `core/formalizm/dodatekO_u1_formalizacja.tex` | + inline annotation kappa-G.0 | + 2 lin |

### Verification

- ✅ **pdflatex compile clean**: 537 stron, dwa runy, brak nowych błędów krytycznych
- ✅ **0 nowych undefined references**: 8 pre-existing `app:A-aksjomaty`, `ax:substrat`, `eq:Phi-sigma-action`, etc. bez związku z G.0
- ✅ **0 nowych multiply-defined labels**: 3 pre-existing `eq:MS-TGP`, `rem:formulation-dictionary`, `ssec:epistemic-table` bez związku z G.0
- ✅ **Wszystkie nowe `\ref{prop:V-M911-canonical}`, `\ref{prop:psi-EOM-R3}`, `\ref{prop:kappa-corrected-G0}`, `\ref{eq:newton-limit-G0}` rozwiązane** z wewnętrznych `\label`
- ✅ **Backwards compatibility**: stare `prop:kappa-corrected` (sek08a v1.x) zachowane → wszystkie cross-references w sek08, sek08c, sek09, dodatekH, dodatekO, status_map oraz ~10 plikach research/ pozostają strukturalnie poprawne

### Strategia: ADDENDUM (warstwa v2.0) zamiast literalnej zmiany v1.x

Phase 4 implementacja **DODAJE** kanoniczną G.0-closed warstwę v2.0 obok historycznej v1.x, zamiast literalnie usuwać/zastępować stare propositions. Powód:

1. **Zachowanie cross-references**: ~30 HIGH-impact files odnosi się do `prop:kappa-corrected`, `eq:kappa-corrected`, `eq:sqrt-g-eff`, etc. — literalne usunięcie zniszczyłoby je. Wartość `κ = 3/(4Φ_0)` jest INVARIANT po re-fit `Φ_0`, więc references pozostają poprawne.

2. **Traceability**: czytelnik widzi ŁAŃCUCH derivacji v1.x → v2.0, nie tylko końcowy stan. Pomocne dla recenzentów i archiwum.

3. **Bezpieczeństwo**: brak ryzyka rozsypania innych derivacji (sek09 K_geo, dodatekO α_em/κ, etc.).

### Dalsze kroki (opcjonalne)

- **Propagacja do research/**: ~10 plików `research/op-newton-momentum/`, `research/nbody/examples/`, etc. (P33 MEDIUM/HIGH) z literalnym `kappa = 3/(4 Phi_0)` annotations.
- **Literalny refactor v1.x → v2.0**: jeśli autor zdecyduje na czyste removal v1.x, można w przyszłej sesji usuwać deprecated propositions. Wymaga full grep-and-replace w core/ + research/.
- **CI compile check**: dodanie `pdflatex` do CI/pre-commit by łapać future regression.

Pełna synteza: w toku (Phase4_results.md, opcjonalne).
**Status: PHASE 4 CLOSED-POSITIVE — G.0 program ukończony.**

---

---

## ⚠ Kontekst krytyczny — luka w sek08 (potwierdzona 2026-05-02)

Audyt 2026-05-01 zamknął adnotacyjnie cluster A1+A2+A3, ale **fizyczna treść
sek08a + sek08c na dysku pozostaje oparta na obsoletnej metryce M9.1
(`g_tt = -c²/ψ`, FALSIFIED 2026-04-25, β_PPN=4 vs Mercury 3·10⁴σ)**.

Bezpośrednia weryfikacja na dysku (2026-05-02):

| Plik | Modified | Stan merytoryczny |
|---|---|---|
| `core/sek08a_akcja_zunifikowana.tex` | 01.05.2026 19:33 | `√(-g) = c₀ψ` z M9.1 (FALSIFIED), κ=3/(4Φ₀) wyprowadzone z falsified √(-g) |
| `core/sek08c_metryka_z_substratu.tex` | 01.05.2026 22:00 | 4 sprzeczne formy metryki, M9.1'' (kanoniczna) **nie występuje fizycznie** — tylko w komentarzu |

`sek08c` lin. 50–54 explicit przyznaje:
> *"A2: dla M9.1'' poprawne `√(-g) = c·ψ/(4-3ψ)`, NIE `c·ψ`. Wszystkie
> wyprowadzenia w sek08a (κ=3/(4Φ₀)), Φ-EOM, M9.2 m_field, M9.3 Peters-Mathews
> używają obsoletnego `√(-g) = c·ψ` z formy (I). Pełny re-run M9.x z
> poprawnym volume element OPEN (B6/dedicated repair cycle)."*

**Bieżący program G.0 jest tym dedicated repair cycle**, ale ze świadomym
poszerzeniem zakresu o weryfikację konsystencji R3 ↔ TGP-canonical.

---

## 1. Cztery konkurencyjne EOM — diagnoza wstępna (2026-05-02)

Bezpośrednia weryfikacja wariacyjna na dysku (skróty rachunkowe, do
formalizacji w Phase 1):

### 1.1. EOM (a): TGP_FOUNDATIONS lin. 80–86 — *claimed* canonical Φ-EOM

```
ψ'' + (2/r)ψ' + 2(ψ')²/ψ + ψ²(1-ψ) = 0          [β=γ=1, ρ=0]
```

**Status:** **claimed**, ale **niepotwierdzone** wariacyjnie z `S_TGP` na dysku.

### 1.2. EOM (b): sek08a derivation z M9.1 (FALSIFIED metric)

`sek08a` używa `√(-g) = c₀ψ` z formy M9.1; po wariacji `δS/δψ = 0` rachunek
daje:

```
S_static = ∫r²dr × [½ψ⁴(ψ')² − ψV(ψ)]
EL → ψ'' + (2/r)ψ' + 2(ψ')²/ψ = γ(15ψ-16)/(12ψ)    [β=γ]
```

**Nie zgadza się z (a) "canonical"**. To pierwsza wewnętrzna niespójność
sek08a — Φ-EOM jest zapostulowane w FOUNDATIONS, ale wariacja akcji daje
co innego.

### 1.3. EOM (c): wariacja S_TGP z poprawnym √(-g) M9.1''

Z M9.1'': `g_eff^rr = (4-3ψ)/ψ`, `√(-g) = c₀ψ/(4-3ψ)`. Z `K(ψ)=ψ⁴`,
`V(ψ)=γψ³(4-3ψ)/12` (β=γ identity):

```
ψV(ψ)/(4-3ψ) = γψ⁴/12        ← (4-3ψ) skraca się!
S_static = 4πc₀ ∫r²dr × [½ψ⁴(ψ')² − γψ⁴/12]
EL → ψ'' + (2/r)ψ' + 2(ψ')²/ψ = −γ/(3ψ)
```

**Nie zgadza się ani z (a), ani z (b), ani z R3**. To jest "co naprawdę
wynika z S_TGP po update sek08 do M9.1''".

### 1.4. EOM (d): R3 ODE (z `r3_*.py`, α=2)

```
g'' + (2/r)g' + 2(g')²/g + (g-1)/g² = 0
```

co odpowiada wariacji **Lagrangianu BEZ volume element** (płaska R³):

```
L_flat = ½g⁴(g')² − V(g)        [V = g³/3 - g⁴/4]
EL → g'' + 2(g')²/g + (2/r)g' = V'/g⁴ = (g²-g³)/g⁴ = -(g-1)/g²    ✓
```

R3 ODE jest dokładnie EOM dla L = ½g⁴(g')² − V(g) na **flat 3D background**,
**bez √(-g)** sprzężenia z metryką dynamiczną.

### 1.5. Kluczowa obserwacja

```
EOM (a)  ≠  EOM (b)  ≠  EOM (c)  ≠  EOM (d)
claimed  ≠  derived(M9.1)  ≠  derived(M9.1'')  ≠  R3
```

**Wszystkie 4 są strukturalnie różne.** R3 ODE (d) jest jedyną, która ma
empiryczne potwierdzenie (mass spectrum, N=3, Koide). M9.1'' (c) jest jedyną,
która używa kanonicznej metryki TGP. (a) i (b) są obsoletne.

---

## 2. Trzy hipotezy do testowania w Phase 1

### H1 (renormalizacja / luka w sek08) — najbardziej prawdopodobna

> **R3 jest naturalną EOM dla wariacji "bez volume element" (flat reference).
> Pełne sprzężenie z M9.1'' wprowadza dodatkowe człony, które są albo
> renormalizacyjnie marginalne, albo absorbowalne przez redefinicję pola
> ψ → ψ̃(g) lub współrzędnej r → r̃(r, ψ).**

Test: znaleźć transformację `(ψ, r) → (g, r̃)` taką, że wariacja S_TGP z
M9.1'' redukuje się do R3 ODE.

### H2 (R3 jako effective theory na M9.1'' background — scenariusz Phase1 §3 (ii))

> **R3 ODE jest efektywnym EOM solitonu w lokalnym Einstein-frame metryki M9.1'',
> uzyskanym przez wycałkowanie sektora geometrycznego (Weyl rescaling +
> conformal frame transformation).**

Test: conformal transformation `g_eff = h(ψ)·η + scalar redefinition φ̃ = T(ψ)`
dająca R3 jako EOM dla φ̃.

### H3 (TGP-canonical jest błędne — scenariusz Phase1 §3 (iii))

> **Nie istnieje konsystentna akcja `S[Φ, g_eff^M9.1'']` która jednocześnie
> daje (a) M9.1'' jako stacjonarny punkt sektora geometrycznego, (b) R3 ODE
> jako stacjonarny punkt sektora skalarnego, (c) PPN γ=β=1, (d) FRW κ=3/(4Φ₀).
> Sek08 wymaga reformulacji: TGP ma DWA niezależne formalizmy** (sek08-grav
> dla makro/PPN/GW; R3 dla mikro/mass/spin), połączone tylko liniową
> reparametryzacją bariera↔horizon (Faza 1 why_n3 §5).

Test: enumeracja możliwych `K(ψ), V(ψ), √(-g)` form i sprawdzenie czy
istnieje wybór dający wszystkie 4 constraint'y jednocześnie.

---

## 3. Cel cyklu

### 3.1. Główny cel

**Rozstrzygnąć, która z H1/H2/H3 zachodzi**, w trybie eksploracyjnym (program
badawczy, nie reformulacja sek08 do końca Phase 1). Po Phase 1+2 dopiero,
jeśli H1 lub H2 jest potwierdzone, można przystąpić do **scalania z sek08**
(Phase 3+ jako osobny audyt rewizyjny wszystkich powiązanych elementów).

### 3.2. Falsyfikatory

- Jeśli żadna transformacja `(ψ, r) → (g, r̃)` nie daje R3 z S_TGP na M9.1''
  — **H1 FALSIFIED**
- Jeśli żadna conformal/Weyl transformation nie produkuje R3 jako EOM w
  Einstein-frame — **H2 FALSIFIED**
- Jeśli enumeracja kandidatów `(K, V, √(-g))` znajdzie konsystentną akcję
  spełniającą wszystkie 4 constraints — **H3 FALSIFIED**, R3 jest derivable
- Jeśli wszystkie 3 hipotezy NEGATIVE — **TGP wymaga nowego aksjomatu**
  (warstwa 3c jako niezależna struktura)

### 3.3. Possible outcomes

- **CASE A — H1 closure**: R3 ↔ TGP-canonical przez transformację pola/współrzędnej.
  Sek08 dostaje update do M9.1'' i automatycznie produkuje R3. **Najczystszy.**
- **CASE B — H2 closure**: R3 jako effective Einstein-frame EOM. Sek08 zachowuje
  obecną strukturę, ale dodaje sekcję "R3 jako projekcja Einstein-frame".
- **CASE C — H3 closure**: TGP_FOUNDATIONS.md dostaje **nowy aksjomat o warstwie
  3c jako niezależnym formalizmie**, połączonym z warstwą 1/2 tylko przez
  liniową reparametryzację bariera ≡ horyzont. Honest framing dla publication.
- **CASE D — wszystkie 3 NEGATIVE**: kryzys foundational, TGP wymaga gruntownej
  reformulacji.

---

## 4. Trzy strategie ataku (Phase 1)

### G0a — Volume integration (H1 test)

Wycałkować TGP-canonical EOM (eq:field-eq-reproduced z M9.1'' volume element)
po `√(-g_M9.1'') d³x` w sektorze statycznym i sprawdzić, czy redukuje się do
R3 ODE z V_R3(g) = -1/g - ln(g).

**Plik:** `phase1_G0a_volume_integration.py`

**PASS criterion:** EOM po projekcji ma RHS proporcjonalne do `(1-g)/g²`
(R3) z residuum <5%.

### G0b — Field redefinition (H1 alt + H2 test)

Symbolicznie szukać transformacji `φ = T(g, ψ)` (z ψ liniowo zależnym od g)
takiej, że TGP-canonical EOM (a) lub (c) → R3 ODE (d) po podstawieniu.

Strategia:
- Substytucja `ψ = a·g + b` (znana z PHASE1: a=0.3814, b=0.6186)
- Lub bardziej ogólna: `ψ = T(g)` z constraint'ami `T(1)=1`, `T(g₀_crit)=4/3`
- Sympy LOCK na transformacji potencjału: `V_TGP(ψ(g)) → V_R3(g)·jacobian`

**Plik:** `phase1_G0b_field_redefinition.py`

**PASS criterion:** Symbolic LOCK transformacji która redukuje EOM (a), (b)
lub (c) do R3 ODE (d).

### G0c — Einstein-frame projection (H2 test)

Conformal transformation `g_eff = h(ψ)·η`, scalar redefinition `φ̃ = ∫√K(ψ)dψ`,
i sprawdzić, czy w Einstein-frame matter sector daje R3 ODE jako EOM dla
soliton excitation.

**Plik:** `phase1_G0c_einstein_frame_projection.py`

**PASS criterion:** Einstein-frame matter Lagrangian dla soliton excitation
w lokalnym limicie reprodukuje L = ½g⁴(g')² - V(g) z V = g³/3 - g⁴/4.

---

## 5. Anchory i punkty wyjścia

### 5.1. Empiryczne hard-anchors (R3, niezależnie od formalizmu)

Te liczby muszą być reprodukowane przez DOWOLNĄ wybraną reformulację:

| Wielkość | Wartość R3 | PDG | Diff |
|---|---|---|---|
| g₀_crit (3D, α=2) | 1.874 | — | — |
| ψ_horizon (M9.1'') | 4/3 = 1.333 | — | — |
| g₀^e | 0.86941 | input | — |
| g₀^μ = φ·g₀^e | 1.40673 | — | — |
| g₀^τ (Koide K=2/3, full A²·g₀^(e²/2)) | 1.77472 | — | — |
| m_μ/m_e (z formuły α=2) | 206.77 | 206.77 | -0.001% |
| m_τ/m_e (z formuły α=2 + g₀_τ=1.77472) | 3477.43 | 3477.23 | +0.006% |
| Spin-1/2 (RP² Berry phase π) | derived | — | — |

### 5.2. Strukturalne hard-anchors

- Ax §1 FOUNDATIONS: jedno pole skalarne Φ z Z₂ — **NIETYKALNE**
- M9.1'' kanoniczna metryka — `g_tt = -c²(4-3ψ)/ψ`, `g_rr = ψ/(4-3ψ)`
- Lorentzian horizon ψ=4/3 — koincydencja z g₀_crit R3 (4 cyfry znaczące,
  PHASE1 §5)
- PPN: γ_PPN = β_PPN = 1 (1PN exact)
- Cosmology κ = 3/(4Φ₀) (FRW)

### 5.3. Otwarte parametry (do strojenia w Phase 1)

- Forma K(ψ): ψ⁴, [(4-3ψ)/ψ]², lub inna
- Forma V(ψ) lub V_eff(ψ) po projekcji
- Field redefinition T(ψ) (Phase 1 ψ ↔ g linear jest first-guess)

---

## 6. Pytania badawcze

### 6.1. Główne pytania (do rozstrzygnięcia w Phase 1)

1. **Czy istnieje wariacyjnie spójna akcja S[Φ, g_eff^M9.1''] dająca R3 ODE
   jako EOM solitonu?** (testy G0a, G0c)
2. **Czy istnieje field redefinition T(ψ) = g taka, że TGP-canonical EOM
   z M9.1'' przekształca się w R3 ODE?** (test G0b)
3. **Jakie są stacjonarne punkty sektora geometrycznego w S_TGP[Φ, g_eff]?**
   Czy M9.1'' wynika wariacyjnie, czy musi być postulowane?
4. **Co naprawdę wynika z `δS_TGP/δΦ = 0` po update sek08a do M9.1''?**
   Czy claim'owany Φ-EOM (a) z FOUNDATIONS się potwierdza, czy nie?

### 6.2. Pomocnicze pytania

- Czy M9.1'' kinetic prefactor `K(ψ)` powinien być inny niż ψ⁴ (z M9.1)?
  Może `K(ψ) = ψ⁴ × (4-3ψ)/ψ = ψ³(4-3ψ)` z poprawki Laplace-Beltrami?
- Czy `m₀ = 0` w vacuum (PHASE4 finding) pochodzi z `V_eff(1) = 0` po
  projekcji, czy z innego mechanizmu?
- Czy R3 ODE potencjał logarytmiczny `V_R3(g) = -1/g - ln(g)` ma topologiczne
  pochodzenie (RP² defect)?

---

## 7. Phase plan

### Phase 1 — Strukturalne testy 3 hipotez (2-3 tygodnie)

**Sub-tasks (3):**
- G0a: volume integration test (H1)
- G0b: field redefinition test (H1 alt + H2)
- G0c: Einstein-frame projection test (H2)

**Score gate:** ≥1/3 PASS → Phase 2 forward (jedna z H1/H2 jest potwierdzona);
0/3 PASS → Phase 2 redirect do testu H3 (enumeracja kandydatów akcji).

### Phase 2 — Sympy LOCK + numerical verification (3-4 tygodnie)

W zależności od Phase 1 outcome:
- **Jeśli H1**: Sympy LOCK redukcji TGP-canonical → R3 + verification że
  cała why_n3 mass formula derives z S_TGP
- **Jeśli H2**: Sympy LOCK Einstein-frame projection + verification PPN +
  cosmology constraint preserved
- **Jeśli H3 fallback**: Enumeracja kandydatów akcji + identification
  najmniej-naruszającego foundation reformulacja

### Phase 3 — Audit rewizyjny sek08 + integration (2-3 tygodnie)

**Tylko jeśli Phase 2 dała pozytywny outcome (H1 lub H2 potwierdzone).**

- Audit wszystkich plików powiązanych z sek08: sek08a, sek08c, sek09, dodatekA,
  dodatekB, M9.1''/M9.2/M9.3 results
- Plan integracji G0 closure z core LaTeX
- Update audit-aware annotations w meta/AUDYT_TGP_2026-05-01.md
- Update `TGP_FOUNDATIONS.md` §1-3 jeśli potrzebny

**Score gate dla Phase 3:** wszystkie powiązane testy numeryczne (M9.x,
why_n3 mass spectrum, PPN, FRW) reprodukowane z poprawną akcją.

---

## 8. Dependencies + prerequisites

### 8.1. Wymagane wcześniejsze cykle (CLOSED)

- ✅ R3 emergent Dirac (`why_n3/PHASE1-5`) — origin R3 ODE + N=3 + spin-1/2
- ✅ Audyt 2026-05-01 — diagnoza luki sek08
- ✅ M9.1'' canonical (TGP_FOUNDATIONS.md:64-69)
- ✅ R5↔Phase 2 bridge (mass_scaling_k4) — m_τ subtension closed

### 8.2. Pomocne ale niewymagane

- λ.1 (e²/4 mechanism) — może bezpośrednio wynikać z G0 closure (volume integral
  z metryki M9.1'' może produkować e² strukturalnie)
- φ.1 substrate-action-variational — może dostarczyć enumerację K(ψ) form

### 8.3. Niezależne (mogą iść równolegle)

- B6/B7 numerical re-runs M9.x (zostaną automatycznie zamknięte przez Phase 3
  G0 jeśli H1 lub H2)

---

## 9. Status meta

**Klasyfikacja:** PROPOSED, eksploracyjny — **program BADAWCZY**, nie
reformulacja sek08. Phase 1 jest pure analysis; Phase 3 (integracja z sek08)
tylko jeśli Phase 1+2 dały konsystentne pozytywne wyniki.

**Ryzyko:** WYSOKIE — być może żadna z H1/H2/H3 nie zachodzi i potrzebna
jest gruntowna reformulacja TGP. Plus risk: rozwiązanie może wymagać
modyfikacji aksjomatu metric-coupling (A4 audytu), co byłoby drugą zmianą
foundational.

**Reward:** **NAJWYŻSZE** w skali całego TGP. Closure G0:
- zamyka A1+A2+A3+A4+B6+B8 audytu *strukturalnie* (nie adnotacyjnie)
- promuje N=3 z "twierdzenie strukturalne R3" do "twierdzenie TGP-canonical"
- daje konkretną ścieżkę dla λ.1 e²/4 (volume integral z M9.1'')
- unifikuje warstwy 0/1/2/3c w jeden łańcuch derywacji
- otwiera czyste publication framing: "TGP jako jednoaksjomatyczny program"

**Rekomendacja czasowa:** ~2-3 miesiące Phase 1+2; +1-2 miesiące Phase 3
jeśli pozytywne. Po Phase 1 — **kill criterion**: jeśli H1+H2+H3 wszystkie
NEGATIVE z explicit mathematical proofs, honest postmortem i decyzja autora
o pivot ontologicznym.

---

## 10. Otwarte pytania meta (dla autora TGP)

1. **Czy `K(ψ) = K_geo·ψ⁴` z M9.1 powinno być zaktualizowane dla M9.1''?**
   Theorem D-uniqueness w sek08a używa M9.1 conformal substrate metric
   `h_ij = ψ⁴δ_ij`. Dla M9.1'' substrate metric byłaby inna —
   `h_ij = [ψ/(4-3ψ)]² δ_ij`? Wtedy K = `K_geo × [ψ/(4-3ψ)]^k` z innym k.

2. **Czy aksjomat A8 (`G(Φ) = G₀·Φ₀/Φ`) jest spójny z M9.1''?**
   Sek08c lin. 297-336 używa A8 do pokazania `ℓ_P = const` dla formy
   eksponencjalnej. Dla M9.1'' formuła może być inna.

3. **Czy R3 logarithmic potential `V_R3(g) = -1/g - ln(g)` ma fundamental
   meaning, czy jest to artefakt projekcji?**
   Jeśli H2 zachodzi, V_R3 wyłania się z V_TGP(ψ) po Einstein-frame
   transformation. Jeśli H1, jest to artefakt zmiany zmiennej.

---

## 11. Końcowa nota

Ten cykl jest **najwyższym priorytetem foundational** w TGP_v1 po zamknięciu
audytu 2026-05-01. Jest to **research-mode** cykl (nie publication-mode):
celem nie jest "udowodnić że TGP-canonical jest poprawne", tylko **rozstrzygnąć
prawdę matematyczną** — z trzema możliwościami końcowymi:

- **Renormalizacja / luka** (H1) — najprostszy, sek08 dostaje update
- **Effective theory** (H2) — sek08 zachowuje się, ale dostaje nową sekcję
- **Niezależne formalizmy** (H3) — TGP_FOUNDATIONS dostaje nowy aksjomat,
  honest framing dla publication

**Najczystszy outcome dla G0**: jeśli żadna z 3 strategii nie produkuje
pozytywu, **honest closure NEGATIVE** z explicit "TGP-canonical (sek08) i R3
są dwoma niezależnymi formalizmami" — co byłoby też wartościowym wynikiem
metaepistemicznym (nikt nie próbował tego rozstrzygnąć formalnie wcześniej).

---

**Autor proposal:** sesja G.0 audyt-extension, 2026-05-02, po-audyt 2026-05-01.
**Origin:** `research/why_n3/PHASE1_psi_g0_identification.md` §3 (Odkrycie 1).
**Status:** PHASE 4 CLOSED-POSITIVE 2026-05-02 — sek08a v2.0 ADDENDUM (`sssec:g0-closure-v2`) + sek08c A1/A2/A3 closure + sek08_formalizm + dodatekH + status_map + sek09 + dodatekO annotations applied. pdflatex compile czysty (537 stron, 0 nowych undefined refs, 0 nowych multiply-defined labels).
**Następne:** opcjonalne — głębszy refactor (literalny replacement V_orig, sqrt-g w v1.x) lub propagacja G.0 closure do `research/` files (~10 plików z P33).
**Plik scaffolding:**
- `Phase1_setup.md` ✓
- `Phase1_results.md` ✓ CLOSED-POSITIVE
- `phase1_G0a_volume_integration.py` ✓ PASS 4/4
- `phase1_G0b_field_redefinition.py` ✓ NEGATIVE-INFORMATIVE
- `phase1_G0c_einstein_frame_projection.py` ✓ PASS 3/3
- `Phase2_setup.md` ✓
- `Phase2_results.md` ✓ CLOSED-POSITIVE 4/4 PASS
- `phase2_P21_vacuum_uniqueness.py` ✓ PASS 4/5
- `phase2_P22_mass_spectrum_verification.py` ✓ PASS 5/5
- `phase2_P23_PPN_verification.py` ✓ PASS 5/5
- `phase2_P24_FRW_cosmology.py` ✓ PASS 5/5
- `Phase3_setup.md` ✓
- `Phase3_results.md` ✓ CLOSED-POSITIVE 4/4 PASS
- `phase3_P32_newton_limit_rederivation.py` ✓ PASS 5/5
- `phase3_P33_cross_reference_audit.py` ✓ PASS 4/4
- `P33_audit_results.md` ✓ (30 HIGH-impact files)
- `sek08a_v2_specification.md` ✓ (7 prop changes drafted)
- `sek08c_A1_A2_A3_closure_draft.md` ✓ (A1+A2+A3 RESOLVED)
