---
title: "M10.2 setup — Inflation audit (ex261, n_s, r, GL(3,F₂) origin)"
date: 2026-04-26
cycle: M10
phase: M10.2
status: SETUP
predecessor: "[[M10_1_results.md]]"
audit_target: "[[../nbody/examples/ex261_inflation_tgp.py]]"
parent: "[[M10_program.md]]"
tags:
  - TGP
  - M10
  - inflation
  - audit-setup
---

# M10.2 — Inflation audit setup (ex261)

> **Cel:** weryfikacja czy istniejący draft [[../nbody/examples/ex261_inflation_tgp.py]]
> jest spójny z **aktualną sek08a akcją** (V=(β/3)φ³-(γ/4)φ⁴, K=K_geo·φ⁴), oraz czy
> jego predykcje (`n_s = 1-2/N_e ≈ 0.967`, `r ≪ 0.036`, GL(3,F₂) p=2N-3=3 origin)
> wynikają strukturalnie z sek08a — czy też wymagają **rebuild from scratch**.
>
> **CRITICAL drift do verify:** ex261 używa V = `(β/7)g⁷ - (γ/8)g⁸` (potęgi 7-8),
> sek08a ma V = `(β/3)φ³ - (γ/4)φ⁴` (potęgi 3-4). Symbolic pre-check pokazuje, że
> **ŻADNE field redefinition g = φ^p nie odwzorowuje tych potencjałów na siebie**
> (potęgi V_7 → V_3 wymagają p=3/7, potęgi V_8 → V_4 wymagają p=1/2 — niespójne).

---

## 1. Cel & predykcje

### 1.1 Cel sub-cyklu

Z drift audit M10.0: ex261 jest YELLOW (V powers mismatch). Audit weryfikuje:

1. **Strukturalnie:** czy ex261 V form jest related do sek08a przez ŻADNĄ field redefinition? (sympy proof)
2. **Hilltop slow-roll w sek08a:** jakie są n_s, r prediction wprost z sek08a `V=(β/3)φ³-(γ/4)φ⁴`, `K=K_geo·φ⁴`?
3. **GL(3,F₂) origin:** czy `p=2N-3=3` z N=3 generations jest derivable z sek08a, czy tylko z ex261 specific action?
4. **Falsifiability:** które predykcje TGP inflation utrzymują się?

### 1.2 Predykcje (TGP — pre-audit)

Z ex261 (10/10 PASS):
- `n_s = 1 - 2/N_e ≈ 0.967` (Starobinsky-class, 0.4σ od Planck 0.9649)
- `r ~ 12/N_e² ≈ 0.0033` (BICEP/Keck-safe, < 0.036)
- `dn_s/dlnk = -2/N_e² ≈ -5.6×10⁻⁴`
- `T_reh ~ 10¹¹ GeV` (sufficient for leptogenesis)
- No monopoles (`π_2(GL(3,F₂))` trivial)

**Pytanie M10.2:** czy te utrzymują się **strukturalnie** dla sek08a action?

### 1.3 Cross-checks z closures

| Closure | Use w M10.2 |
|---------|-------------|
| sek08a (V=(β/3)φ³-(γ/4)φ⁴, β=γ) | główna referencja |
| sek08a (K=K_geo·φ⁴) | non-canonical kinetic w canonical frame |
| M10.1 (V''(1)=-β slow-roll max) | relevant: φ=1 jest hilltopem inflacji? |
| T-Λ (Φ_eq=H_0, V(1)=β/12) | residual Λ source post-inflation |
| T-α (α(ψ)=α_0(ψ-1)²Θ(ψ-1)) | post-inflation (ψ_NS regime) — n/a |

---

## 2. Setup matematyczny

### 2.1 Akcja TGP (sek08a)

```
S_TGP = ∫ d⁴x √(-g_eff) [ ½ K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ) - (q/Φ_0) φ ρ ]
K(φ)  = K_geo φ⁴
V(φ)  = (β/3)φ³ - (γ/4)φ⁴,    β = γ
```

### 2.2 Akcja ex261 (drift)

```
L_ex261 = ½ g⁴ (∂g)² + (β/7) g⁷ - (γ/8) g⁸
```

W zwykłej Minkowski (no g_eff). Powers w V: 7 i 8.

### 2.3 Field redefinition test (g = φ^p)

Jeśli g = φ^p:
- Kinetic: g⁴(∂g)² = p²·φ^(6p-2)(∂φ)²
- V power 7: (β/7)φ^(7p)
- V power 8: (γ/8)φ^(8p)

Dla matching sek08a:
- Kinetic K_geo·φ⁴: 6p-2 = 4 → **p=1** (g = φ rename only)
- V power 3: 7p=3 → p=3/7
- V power 4: 8p=4 → p=1/2

**Niespójne!** Żadne single-power g = φ^p nie odwzorowuje ex261 → sek08a.

**Konkluzja:** ex261 i sek08a to **structurally different actions** w sensie field potential.

### 2.4 Sek08a slow-roll w canonical χ frame

Canonicalize: dχ/dφ = √K = √K_geo · φ²
```
χ(φ) = √K_geo · φ³ / 3
φ(χ) = (3χ/√K_geo)^(1/3)
```

W kanonicznym χ frame (K_geo=1, β=γ=1):
```
V(χ) = (1/3)·(3χ)^(3·1/3) - (1/4)·(3χ)^(4·1/3)
     = χ - (3^(4/3)/4)·χ^(4/3)
     = χ − 1.082·χ^(4/3)
```

Hilltop:
```
V'(χ) = 1 − 3^(1/3)·χ^(1/3) = 0 → χ_max = 1/3
V(χ_max) = β/12   (→ T-Λ residual)
V''(χ_max) = −β
```

**Slow-roll na hilltopie sek08a:**
```
η = V''/V = −β / (β/12) = −12
```

**|η| = 12 ≫ 1 — slow-roll FAILS** w sek08a hilltop bez additional plateau V₀.

**Wniosek:** sek08a w canonical χ frame **nie wspiera bezpośrednio** slow-roll inflation. Slow-roll wymaga (a) dodatkowego skalu V₀ ≫ V_max (plateau-extended hilltop), albo (b) różnego frame transformation (e.g., conformal Einstein), które ex261 zastosował dla swojej akcji `g⁷-g⁸`.

### 2.5 Hilltop with plateau V₀

Jeśli sub-leading V₀ ≫ V_max (η-suppression):

```
V_infl(χ) = V₀ - (β/2)·12·(χ-1/3)²·V_max + ...
          ≈ V₀ · [1 - (V_max/V₀)·12·(χ-1/3)² + ...]
```

Effective hilltop with η_infl = −12·V_max/V₀.

Dla V₀/V_max ~ 100: |η_infl| ~ 0.12. Still ≫ 0.01 needed for n_s ≈ 0.96.

Konkluzja: sek08a hilltop wymaga BARDZO dużego plateau V₀ (V₀/V_max ~ 1000+) aby η_infl ≪ 1.

**To nie jest "natural inflation"** w sensie minimal sek08a — wymaga dodatkowego parametru.

### 2.6 Slow-roll near origin (φ → 0)

Alternative: inflation near φ=0 (degenerate critical point V(0)=V'(0)=V''(0)=0, V'''(0)=2β).

W canonical χ frame, χ → 0 odpowiada φ → 0:
```
V(χ) ≈ χ      (linear leading dla małych χ)
V'(χ) ≈ 1
ε(χ) = (V'/V)²/2 = 1/(2χ²) → ∞ as χ → 0
```

ε diverges — slow-roll fails na drugim końcu. Jedyny region z ε~1 jest pomiędzy χ=0 i χ=1/3, ale tam |η|→∞.

**Wniosek:** sek08a sam w sobie nie ma "good slow-roll region". Inflation TGP wymaga dodatkowych elementów (e.g., conformal couplings, plateau, multi-field).

### 2.7 ex261 vs sek08a — structural comparison

| Aspect | ex261 | sek08a |
|--------|-------|--------|
| V powers | 7, 8 | 3, 4 |
| Kinetic | g⁴(∂g)² | K_geo·φ⁴(∂φ)² |
| Background | Conformal g_μν=g²η_μν | Hyperbolic g_eff (M9.1'') |
| Hilltop pos | g→1 (P max) | χ=1/3 (V max) in canonical |
| Slow-roll near hilltop | OK (with V₀ from GL(3,F₂)/168) | FAILS (|η|=12) |
| n_s | 1-2/N_e ≈ 0.967 | undefined w/o V₀ plateau |
| r | ~10⁻³ | undefined w/o V₀ plateau |

**Strukturalnie:** ex261 i sek08a są DIFFERENT theories. ex261's predictions NIE są derived z sek08a directly.

---

## 3. Plan testów (6 sub-testów)

### M10.2.1 — Action structure drift confirmation (sympy)

**Cel:** sympy proof, że ex261 action ≠ sek08a action.

**Metoda:** test field redefinitions g = φ^p.

**Sub-tests:**
- (a) p=1 (rename): kinetic OK, ale V powers 7-8 vs 3-4 differ.
- (b) Match V_7→V_3 demands p=3/7; match V_8→V_4 demands p=1/2 → INCONSISTENT.
- (c) No single field redefinition links the two.
- (d) Conclude: ex261 V form is OBSOLETE/UNRELATED to sek08a.

### M10.2.2 — Sek08a canonical χ frame & hilltop (sympy)

**Cel:** derive V(χ) w canonical frame i pokazać hilltop properties.

**Metoda:** sympy:
- χ = √K_geo · φ³/3
- V(χ) = symbolic
- V'(χ_max) = 0 → χ_max = 1/3
- V_max = V(1/3) = β/12
- V''(χ_max) = −β

**Sub-tests:**
- (a) χ_max = 1/3 — sympy
- (b) V_max = β/12 — matches T-Λ residual ✓
- (c) V''(χ_max) = −β — slow-roll max
- (d) η_hilltop = V''/V_max = −12 → |η|=12 ≫ 1 → **slow-roll FAILS bez plateau**

### M10.2.3 — Slow-roll near origin in sek08a (numerical)

**Cel:** pokazać że slow-roll w canonical sek08a fails wszedzie.

**Metoda:** numerical:
- ε(χ), η(χ) na siatce χ ∈ (0, 1.5)
- check ε<1 i |η|<1 simultaneously
- Find optimal slow-roll window (jeśli istnieje)

**Sub-tests:**
- (a) ε(χ) → ∞ as χ → 0 (linear V)
- (b) |η(χ)| → ∞ as χ → 0 (V'' diverges)
- (c) χ_max=1/3 daje η=−12
- (d) **No good slow-roll window** w sek08a w canonical frame

### M10.2.4 — ex261 with plateau V₀ (consistency check)

**Cel:** verify że ex261's `n_s=1-2/N_e` jest internally consistent dla GL(3,F₂)/168 plateau.

**Metoda:** mimic ex261 calculation:
- V_infl(g) = V₀(1 - (g/μ)³ + ...)
- p=3 hilltop, N_e=60
- Compute n_s, r

**Sub-tests:**
- (a) p=3 hilltop → n_s = 1 - 5/(3N_e) ≈ 0.972 (Boubekeur-Lyth)
- (b) ex261 simplification n_s = 1 - 2/N_e ≈ 0.967 (Starobinsky scaling)
- (c) Both within 1-2σ of Planck 0.9649 ± 0.0042
- (d) r ~ 12·9/(4·N_e²) ≈ 7.5×10⁻³ < 0.036 — BICEP-safe

### M10.2.5 — GL(3,F₂) origin claim verification (structural)

**Cel:** verify ex261's claim że `p = 2N-3 = 3` z N=3 generations ma structural origin.

**Metoda:** sympy:
- ex261: V = g^(2N+1)/(2N+1) - g^(2N+2)/(2(N+1))
- For N=3: V = g⁷/7 - g⁸/8 ✓ (matches ex261 formula)
- For N=2: V = g⁵/5 - g⁶/6 (different theory, p=2N-3=1 — not hilltop)
- Check: n_s sensitivity to N

**Sub-tests:**
- (a) ex261's V = g^(2N+1)/(2N+1) - g^(2N+2)/(2N+2) for general N
- (b) For N=3: matches g⁷-g⁸ ✓
- (c) Conformal frame V_eff = -g^N/(2N+1) + g^(N+1)/(2N+2) (leading power N=3)
- (d) p_eff = 2N-3 → for N=3: p=3 ✓ (cubic hilltop)
- (e) N=3 from generations is structural input, not derived

### M10.2.6 — Honest verdict on TGP inflation predictions

**Cel:** synthesize honest verdict.

**Metoda:** documentation only.

**Sub-tests:**
- (a) ex261 action ≠ sek08a action — confirmed (M10.2.1)
- (b) sek08a does NOT directly support slow-roll inflation w canonical frame (M10.2.2-3)
- (c) ex261's predictions assume specific (older) action; structurally separate from sek08a
- (d) **HONEST VERDICT:** ex261 jest YELLOW (drift confirmed) lub RED (rebuild from sek08a needed)
- (e) Open question: czy TGP inflation derivable z full hyperbolic metric M9.1'' z odpowiednim V₀? (M10.2 doesn't resolve — flag for future cycle)

---

## 4. Numerical realization (m10_2_inflation.py outline)

```python
# Top: sympy
def test_M10_2_1():  # action drift
def test_M10_2_2():  # canonical frame hilltop
def test_M10_2_5():  # GL(3,F2) origin

# Middle: numerical
def epsilon_eta(V_func, Vp_func, Vpp_func, chi_grid):
    """Compute slow-roll params on grid"""

def test_M10_2_3():  # sek08a slow-roll fails

# ex261-style hilltop
def test_M10_2_4():  # plateau V0 with p=3
def test_M10_2_6():  # honest verdict
```

**Bezpieczne defaults:**
- `K_geo = 1`, `β = γ = 1` (normalized)
- N_e = 60 (canonical)
- `Planck 2018: n_s = 0.9649 ± 0.0042`
- `BICEP/Keck 2021: r < 0.036`

**Output:** `m10_2_inflation.txt`. Expected verdict: 4-5/6 PASS (honest documentation).

---

## 5. Falsifiable predictions (post-M10.2)

1. **TGP inflation ≠ minimal sek08a:** sek08a w canonical frame nie produkuje natural slow-roll. TGP inflation requires dodatkowy element (V₀ plateau, conformal frame, modified action).
2. **n_s = 1-2/N_e (heuristic):** matches Planck within 1σ jeśli sek08a ma efektywny p=3 hilltop. Ale to NIE jest derived z minimal sek08a.
3. **r ~ 10⁻³:** BICEP-safe; LiteBIRD detectable. Predyktywne dla TGP plateau scenarios.
4. **GL(3,F₂)/N=3 origin:** structural input from generation count; nie derived first-principles z sek08a.

## 6. Drift-check matrix (post-execution expected)

| Constraint | Test | Expected status |
|------------|------|-----------------|
| sek08a `V=(β/3)φ³-(γ/4)φ⁴` | M10.2.1 vs ex261's `(β/7)g⁷-(γ/8)g⁸` | DRIFT confirmed |
| sek08a `K=K_geo·φ⁴` | M10.2.2 canonical frame | OK in principle |
| Hilltop near vacuum | M10.2.2 sek08a vs ex261 | sek08a needs plateau |
| n_s=1-2/N_e (Starobinsky scaling) | M10.2.4 with V₀ | consistent with ex261 |
| r ≪ 0.036 | M10.2.4 hilltop p=3 | BICEP-safe |
| GL(3,F₂)/N=3 origin | M10.2.5 structural | structural input |

## 7. Następne (post-M10.2)

Jeśli verdict YELLOW (4-5/6 PASS):
- ex261 → "structural reference" status (nie closure-grade)
- M10.2 documents drift; ex261's predictions are HEURISTIC consistency check
- Future M11 task: derive TGP inflation from full hyperbolic metric M9.1'' + sek08a action

Jeśli verdict RED (≤3/6 PASS):
- ex261 → REBUILD required (next cycle)
- ex261 archived as historical reference

Most likely: **5/6 PASS** with M10.2.6 documenting honest drift, ex261 keeps qualitative status.

---

*M10.2 setup completed 2026-04-26. Ready for execution.*
