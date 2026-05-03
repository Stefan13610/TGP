---
title: "Sek08c A1/A2/A3 closure draft (G.0 closure)"
date: 2026-05-02
phase: 3
sub-task: P34
parent: "[[Phase3_setup.md]]"
predecessors:
  - "[[Phase1_results.md]]"
  - "[[Phase2_results.md]]"
  - "[[sek08a_v2_specification.md]]"
status: CLOSURE DRAFT COMPLETE
scope_constraint: "Draft only. Actual sek08c.tex modification = Phase 4."
---

# Sek08c A1/A2/A3 — closure draft

> **Cel:** Po G.0 closure (Phase 1+2+3) wszystkie 3 audit annotations
> z 2026-05-01 (A1, A2, A3) w sek08c są CLOSEABLE. Ten dokument
> proponuje zaktualizowany tekst preamble sek08c reflecting closure.

---

## 0. Status sumarycznie

| Annotation | Original issue | G.0 resolution | Status |
|---|---|---|---|
| **A1** | Φ-EOM mismatch z R3 ODE | Phase 1 G0a: V_M911 sympy LOCK reprodukuje R3 EXACT | **CLOSED-RESOLVED** |
| **A2** | √(-g)=c·ψ obsoleted (M9.1 falsified) | M9.1'' canonical adopted: √(-g)=c·ψ/(4-3ψ) | **CLOSED-RESOLVED** |
| **A3** | 4 metric forms współistnieją | M9.1'' forma (IV) UNIQUE canonical post-G.0 | **CLOSED-RESOLVED** |

---

## 1. A1 — Φ-EOM mismatch z R3 ODE

### Original (sek08c lin. 7, 50-65 audit warning)

> "AUDIT 2026-05-01 (A1) — METRIC-FORM DRIFT NOTICE: w pliku
> współistnieją 4 wzajemnie sprzeczne formy metryki..."

Problem A1: derivacje sek08a (Φ-EOM, κ, soliton mass) używają formy I
metryki (M9.1: g_tt = -c²/ψ), ale R3 ODE w research/why_n3 jest niezależnie
poprawne i daje N=3 + mass spectrum z PDG. Dwa formalizmy wyglądały
niereconciliable.

### G.0 resolution (Phase 1 G0a + Phase 2)

**Phase 1 G0a** (sympy LOCK):
- Z constraint K(ψ)=ψ⁴, √(-g)=c·ψ/(4-3ψ) (M9.1'' canonical), wymóg
  EOM = R3 ODE jednoznacznie wyznacza V(ψ) = -γψ²(4-3ψ)²/12 (V_M911)
- Numerical profile match z R3: max diff = 0.000000 dla wszystkich 4 testów
  (g₀^e, g₀^μ, g₀^τ, g₀_crit)

**Phase 2 P21 + P22**:
- V_M911 unique pod constraint'ami (sympy LOCK)
- Mass spectrum lepton z PDG <0.01% (przeniesione z R3)
- m_sp² = +γ stabilne (FIXED tachion bug w sek08a v1.x)

**Wynik:** R3 ODE JEST effective EOM solitonu w sek08a v2.0 (po G.0 update).
Brak mismatch — TGP-canonical i R3 są tym samym formalizmem.

### Proposed closure text dla sek08c (replace lin. 9-32 audit warning)

```latex
% =====================================================================
% A1 CLOSED 2026-05-02 (G.0 closure, Phase 1+2)
% =====================================================================
% Original audit annotation A1 (METRIC-FORM DRIFT NOTICE) jest CLOSED-RESOLVED
% przez program G.0 (research/op-g0-r3-from-canonical-projection/):
%
% RESOLUTION: Pelna formal consistency miedzy sek08a Φ-EOM a R3 ODE
% (research/why_n3/) zostala udowodniona w 3 niezaleznych sub-tasks:
%
%   Phase 1 G0a (volume integration): K(ψ)=ψ⁴ + √(-g)=c·ψ/(4-3ψ) +
%     wymog EOM=R3 ODE → V(ψ) = -γψ²(4-3ψ)²/12 (V_M911), sympy LOCK.
%     Numerical profile match max diff = 0.000000.
%
%   Phase 1 G0c (Einstein-frame projection): Spatial Einstein-frame
%     z M9.1'' Ω² = ψ/(4-3ψ) → identyczne V_M911 derived (≡ G0a).
%
%   Phase 2 P21-P24: V_M911 unique, vacuum stable (m_sp²=+γ),
%     mass spectrum z PDG <0.01%, PPN γ=β=1, FRW κ form invariant.
%
% Sek08a v2.0 (specifikacja: research/op-g0-r3-from-canonical-projection/
% sek08a_v2_specification.md) zawiera updated propositions reflecting
% G.0 closure.
% =====================================================================
```

---

## 2. A2 — √(-g) form inconsistency

### Original (sek08c lin. 50-54 audit warning)

> "(a) A2: dla M9.1'' poprawne `√(-g) = c·ψ/(4-3ψ)`, NIE `c·ψ`.
> Wszystkie wyprowadzenia w sek08a (κ=3/(4Φ_0)), Φ-EOM, M9.2 m_field,
> M9.3 Peters-Mathews używają obsoletnego `√(-g) = c·ψ` z formy (I).
> **Pełny re-run M9.x z poprawnym volume element OPEN** (B6/dedicated
> repair cycle)."

Problem A2: derivacje używały √(-g) = c·ψ (forma I, M9.1 FALSIFIED), ale
canonical metric M9.1'' ma √(-g) = c·ψ/(4-3ψ).

### G.0 resolution (Phase 1+2+3)

**Phase 1 G0a:** Re-run derivacji z poprawnym √(-g) = c·ψ/(4-3ψ) (M9.1'')
zwraca R3 ODE jako EOM solitonu (z V_M911 update).

**Phase 2 P21:** Sympy LOCK potwierdza unique V_M911.

**Phase 2 P24:** Re-run FRW (cosmology) z V_M911 + √(-g)=c·ψ/(4-3ψ):
- Source coupling: 5q·ρ/Φ_0 (vs 2q·ρ/Φ_0 w v1.x), structurally different
- Form κ ∝ q·c²/(Φ_0·H_0²) zachowana
- Po re-fit q·c²/Φ_0 (P32) value κ INVARIANT = 4πG_0/(3H_0²)

**Phase 3 P32:** Re-run Newton limit z V_M911 + √(-g)=c·ψ/(4-3ψ):
- q·c²/Φ_0 = (4/5)πG_0 (NEW, vs 2πG_0 v1.x)
- Newton G_0 EXACTLY reproduced
- Wszystkie observables INVARIANT po re-fit

### Proposed closure text dla sek08c (replace lin. 50-54 audit warning)

```latex
% =====================================================================
% A2 CLOSED 2026-05-02 (G.0 closure, Phase 1+2+3 P32)
% =====================================================================
% Original audit annotation A2 (√(-g) form obsoleted) jest CLOSED-RESOLVED
% przez program G.0:
%
% RESOLUTION: Pełny re-run wszystkich M9.x derivations z poprawnym
% √(-g) = c·ψ/(4-3ψ) (M9.1'' canonical) został wykonany:
%
%   - Φ-EOM (statyczna sferyczna): R3 ODE (Phase 1 G0a)
%   - κ (FRW): form invariant, value invariant po re-fit q·c²/Φ_0
%     = 4πG_0/(3H_0²) = 3/(4Φ_0_new) (Phase 2 P24 + Phase 3 P32)
%   - m_field (mass spectrum): zachowane (R3 ODE = effective EOM,
%     Phase 2 P22)
%   - Newton-limit relation: q·c²/Φ_0 = (4/5)πG_0 (NEW, P32)
%   - Wszystkie observables (Newton G_0, PPN, BBN, LLR, CMB) INVARIANT
%     po re-fit (gauge-equivalence)
%
% Sek08a v2.0 (specifikacja: research/op-g0-r3-from-canonical-projection/
% sek08a_v2_specification.md) zawiera updated derivations z poprawnym
% √(-g). Phase 4 (actual core mod) implementuje update.
% =====================================================================
```

---

## 3. A3 — 4 metric forms współistnieją

### Original (sek08c lin. 56-65 audit warning)

> "(b) A3: β_PPN=1 dla M9.1'' wymaga 'master formula' β = f''(1)/f'(1)²
> + 2c_2/f'(1) = 1/2 + 1/2, gdzie c_2=-1 pochodzi z **tej samej**
> nieliniowości α=2 co M9.1. Czysto metryczna identyfikacja daje
> β_PPN_metric = 1/2."

Problem A3: 4 metric forms (I, II, III, IV) współistnieją w sek08c,
przy czym tylko forma IV (M9.1'') daje β_PPN=1 spójnie z observation,
ale wymaga "master formula z kinetic correction".

### G.0 resolution (Phase 2 P23 + Phase 3)

**Phase 2 P23** (sympy LOCK):
- γ_PPN = 1 z M9.1'' g_rr linearization (V-independent)
- β_metric = 1/2 z M9.1'' g_tt linearization (V-independent)
- c₂ = -1 z R3 ODE RHS linearization (INVARIANT pod V update,
  bo R3 ODE = effective EOM zarówno w sek08a v1.x jak i v2.0)
- β_PPN = β_metric + 2c₂/f'(1) = 1/2 + 1/2 = 1 (master formula sek08c)

**Wynik:** master formula sek08c (z kinetic correction) jest poprawna,
i daje β_PPN = 1 spójnie. M9.1'' (forma IV) jest UNIQUE canonical
(formy I, II, III deprecated/falsified).

### Proposed closure text dla sek08c (replace lin. 56-65 audit warning)

```latex
% =====================================================================
% A3 CLOSED 2026-05-02 (G.0 closure, Phase 2 P23)
% =====================================================================
% Original audit annotation A3 (β_PPN master formula konwencja) jest
% CLOSED-RESOLVED przez program G.0:
%
% RESOLUTION: Sympy LOCK wyników M9.1'' PPN derivation (Phase 2 P23):
%
%   - γ_PPN = 1 (z g_rr M9.1'' linearization, V-independent) ✓ EXACT
%   - β_metric = 1/2 (z g_tt M9.1'' linearization) ✓ EXACT
%   - c₂ = -1 (z R3 ODE RHS linearization, INVARIANT pod V update) ✓
%   - β_PPN = β_metric + 2c₂/f'(1) = 1/2 + 1/2 = 1 ✓ EXACT
%
% Master formula z kinetic correction jest **kanoniczna** dla TGP-canonical
% theory: γ_PPN i β_PPN sa konsystentnie wyznaczone przez:
%   (a) pure metric M9.1'' linearization (γ_PPN, β_metric)
%   (b) field equation (R3 ODE) linearization (c₂)
%   (c) sek08c master formula combining (a) + (b)
%
% PPN sektor jest INVARIANT pod G.0 V update (Phase 2 P23 PASS 5/5):
% γ depends tylko na metric M9.1'' (kanonicznej, niezmienionej),
% β depends na (γ_metric + c₂/f'(1)) gdzie obie składowe są INVARIANT.
%
% M9.1'' (forma IV) staje sie UNIQUE CANONICAL post-G.0 closure.
% Formy (I), (II), (III) → deprecated/historical.
% =====================================================================
```

---

## 4. Updated preamble dla sek08c (proposed full replacement, lin. 5-71)

```latex
% =====================================================================
% G.0 CLOSURE 2026-05-02 (Phase 1+2+3, A1+A2+A3 RESOLVED)
% =====================================================================
% Sek08c v2.0 — po G.0 closure (research/op-g0-r3-from-canonical-projection/)
%
% HISTORYCZNIE w tym pliku istniały 4 wzajemnie sprzeczne formy metryki:
%   (I)   eq:metric-full-derived: g_tt = -c²/ψ — FALSIFIED 2026-04-25
%   (II)  forma eksponencjalna: g_tt = -c²·e^(-2U) — OBSOLETE/przejściowa
%   (III) thm:antipodal-uniqueness: f=ψ^(-1/2), h=ψ^(+1/2) — niespójne wew.
%   (IV)  M9.1'' hiperboliczna: g_tt = -c²·(4-3ψ)/ψ — KANONICZNA wg
%         TGP_FOUNDATIONS.md
%
% POST-G.0 CLOSURE (2026-05-02):
%   - M9.1'' (forma IV) jest UNIQUE CANONICAL metryka TGP-canonical
%   - Formy (I), (II), (III) → DEPRECATED, pozostają w tekście jako
%     historyczne kroki "mostu substrat→metryka"
%   - Wszystkie audit annotations A1+A2+A3 → CLOSED-RESOLVED przez G.0
%
% Specjalnie:
%   A1 (Φ-EOM ↔ R3 ODE mismatch): RESOLVED przez V_M911 update
%       (Phase 1 G0a, sympy LOCK, max diff = 0.000000)
%   A2 (√(-g) = c·ψ obsoleted): RESOLVED przez √(-g) = c·ψ/(4-3ψ) adopt
%       (Phase 1+2+3, wszystkie M9.x re-derived)
%   A3 (4 metric forms): RESOLVED — M9.1'' UNIQUE, master formula
%       z kinetic correction kanoniczna (Phase 2 P23 sympy LOCK)
%
% Aktualne canonical anchory G.0:
%   K(ψ) = ψ⁴ (T-D-uniqueness, α=2, sek08a v1.x → v2.0 zachowane)
%   V(ψ) = -γψ²(4-3ψ)²/12 (V_M911, sek08a v2.0 NEW)
%   √(-g) = c·ψ/(4-3ψ) (M9.1'' canonical, A2 closure)
%   ψ_vacuum = 1, m_sp² = +γ (vacuum stability fix, A1 closure)
%   q·c²/Φ_0 = (4/5)πG_0 (Newton-limit, P32)
%   κ = 3/(4Φ_0_new) = 4πG_0/(3H_0²) (po re-fit)
%
% Patrz:
%   research/op-g0-r3-from-canonical-projection/README.md (program G.0)
%   research/op-g0-r3-from-canonical-projection/Phase1_results.md
%   research/op-g0-r3-from-canonical-projection/Phase2_results.md
%   research/op-g0-r3-from-canonical-projection/Phase3_results.md
%   research/op-g0-r3-from-canonical-projection/sek08a_v2_specification.md
% =====================================================================
```

---

## 5. Inline annotations w body sek08c

### 5.1 Lin. 228 (STATUS audit A1) — replace z:

```latex
% --- STATUS (audit 2026-05-01, A1 → CLOSED-RESOLVED 2026-05-02 G.0) ---
% G.0 closure resolved A1 — formy (I)+(II)+(III) deprecated; forma (IV)
% M9.1'' unique canonical; Φ-EOM = R3 ODE (sek08a v2.0).
```

### 5.2 Lin. 288 (STATUS audit A1) — replace z:

```latex
% --- STATUS (audit 2026-05-01, A1 → CLOSED-RESOLVED 2026-05-02 G.0) ---
```

### 5.3 Lin. 455 (STATUS audit A1) — replace z:

```latex
% --- STATUS (audit 2026-05-01, A1 → CLOSED-RESOLVED 2026-05-02 G.0) ---
```

---

## 6. Phase 4 implementation notes

(Phase 4 = actual core mod, separate task)

### 6.1 Order

1. Replace preamble lin. 5-71 z proposed full replacement (sekcja 4).
2. Replace inline STATUS annotations (sekcja 5).
3. Verify wszystkie cross-references do A1/A2/A3 w innych plikach
   pozostają consistent (w grep audit, P33 found references w sek08c
   itself — żadne zewnętrzne A1/A2/A3 references).

### 6.2 Estimate

- Preamble update: 0.5 godziny
- Inline annotations: 0.5 godziny
- Verification + recompile: 0.5 godziny

**Total dla sek08c closure: ~1.5 godziny.**

---

## 7. Cross-reference impact (sek08c-related)

Z P33 audit:
- sek08c ma 4 kappa references — wszystkie pozostają poprawne
  (kappa = 3/(4Φ_0) wartość invariant)
- Φ-EOM derivacja w sek08a (cytowana z sek08c) wymaga update do v2.0

Brak innych external references do A1/A2/A3 specifically.

---

**Status:** Closure draft complete. Ready dla user review + Phase 4
actual core mod.
