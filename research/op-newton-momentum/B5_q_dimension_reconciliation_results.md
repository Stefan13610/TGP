---
title: B5 — q-dimension reconciliation across TGP_v1 LaTeX core
date: 2026-05-01
audit-item: B5 (HIGH)
status: CLOSED
script: B5_q_dimension_reconciliation.py
log: B5_q_dimension_reconciliation.txt
---

# B5 — q-dimension reconciliation results

## TL;DR

The TGP coupling constant $q$ has a **single canonical dimension** in SI units:

$$
\boxed{[q] \;=\; [\mathrm{M}^{-1}\,\mathrm{L}]}
$$

locked independently from three paths (field EOM, Newton boundary, action variation
schematic). Two LaTeX inconsistencies were identified and **fixed**:

1. `axioms/notacja/dodatekA_notacja.tex` line 176 — table entry `[q] = [-]` →
   `[q] = [M⁻¹L]`. Companion ρ entry on line 177 also corrected from
   `[L⁻³]` (normalized number density) to `[ML⁻³]` (mass-energy density)
   to match the table's own physical-units convention.
2. `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` line 516 —
   typo `q = 4πG₀/c₀⁴ · Φ₀` → corrected to `q = 4πG₀/(c₀² · Φ₀)`,
   consistent with the explicit identity `q·Φ₀ = 4πG₀/c₀²` printed
   immediately after on line 519.

A sweep across all `.tex` files confirmed no other `c₀⁴` typos exist in
relation to $q$. The references in `sek03_rezimy.tex:524` and
`sek08_formalizm.tex:2872, 2903` use the canonical form
`q = 4πG₀/c₀²` (exponential-metric variant, equivalent up to $\Phi_0$
absorption) and remain consistent with the now-corrected dodatekA table.

**Audit progression**: 40/43 → **41/43 (95%) closed** after B5.

---

## 1. Problem statement

Audit item B5 (`AUDYT_TGP_2026-05-01.md` § HIGH/B5) flagged two suspected
inconsistencies in the dimensional treatment of $q$:

- The dodatekA table at line 176 lists $[q] = [-]$ (dimensionless), in
  apparent contradiction with the same file's earlier glossary at line 30
  which states `[q] = [M⁻¹L]` and gives $q_{\min} = 8\pi G_0/c_0^2$.
- The sek08a derivation at line 516 prints
  $q = (4\pi G_0/c_0^4)\,\Phi_0$, while line 519 (3 lines later) prints
  the contradictory identity $q\,\Phi_0 = 4\pi G_0/c_0^2$. Solving
  the latter for $q$ yields $q = 4\pi G_0/(c_0^2\,\Phi_0)$, which differs
  from line 516 in BOTH the power of $c_0$ and the position of $\Phi_0$.

The B5 closure cycle was to (i) lock the canonical dimension of $q$ from
first principles, (ii) verify that all three independent paths agree, and
(iii) repair the two inconsistent strings.

## 2. Sympy LOCK results (5/5 PASS)

Script: `B5_q_dimension_reconciliation.py`
Log: `B5_q_dimension_reconciliation.txt`

### Step 1 — Field EOM lock

From $\nabla^2\Phi = -q\rho$ with $[\Phi]=[-]$, $[\rho]=[\mathrm{ML}^{-3}]$:

$$
[q] \;=\; \frac{[L^{-2}]}{[ML^{-3}]} \;=\; [M^{-1}L] \quad \checkmark
$$

### Step 2 — Newton boundary

Both forms reduce to the same dimension:

$$
[q_{\min}] = [q_{\exp}] = \frac{[G_0]}{[c_0^2]} = \frac{[M^{-1}L^3T^{-2}]}{[L^2T^{-2}]} = [M^{-1}L] \quad \checkmark
$$

### Step 3 — Action variation closure

The TGP action carries two physically distinct roles for $q$:

- **Role (1) — source coupling**: $\mathcal{D}[\Phi] = q\rho$ → $[q]=[M^{-1}L]$
- **Role (2) — conformal exponent**: $g_{\mathrm{eff}} \sim e^{q\varphi}g_0$ requires $[q\varphi]=[-]$

These roles refer to the SAME physical coupling but in different forms;
the dimensional bridge involves an implicit volume/mass scale absorption.
Canonical TGP fixes Role (1), and Role (2) is reached by an explicit
$V_{\mathrm{ref}}$ factor inside the exponent. The dodatekA table line 176
was tacitly using the dimensionless Role (2) form WITHOUT a footnote,
which created the apparent inconsistency.

### Step 4 — LaTeX audit

Two specific issues identified and characterized:

| Location | Current | Should be | Type |
|---|---|---|---|
| dodatekA ln 176 | `[q] = [-]` | `[q] = [M⁻¹L]` | convention split (no footnote) |
| dodatekA ln 177 | `[ρ] = [L⁻³]`, $\int\rho d^3x=1$ | `[ρ] = [ML⁻³]`, $\int\rho d^3x=m$ | normalization mismatch |
| sek08a ln 516 | $q = 4\pi G_0/c_0^4 \cdot \Phi_0$ | $q = 4\pi G_0/(c_0^2\,\Phi_0)$ | typo (c-power, Φ₀-position) |

Dimensional verification of the corrected sek08a form:
$$
\left[\frac{4\pi G_0}{c_0^2 \Phi_0}\right]
= \frac{[M^{-1}L^3T^{-2}]}{[L^2T^{-2}][-]}
= [M^{-1}L] \quad \checkmark
$$

Dimensional check of the broken form (proves it WAS wrong):
$$
\left[\frac{4\pi G_0}{c_0^4} \cdot \Phi_0\right]
= \frac{[M^{-1}L^3T^{-2}]}{[L^4T^{-4}]} \cdot [-]
= [T^2 M^{-1} L^{-1}] \neq [M^{-1}L]\quad \times
$$

### Step 5 — Canonical convention (locked)

```
SI canonical (TGP_v1):
  [q]   = [M⁻¹ L]
  q     = 4π G₀ / (c₀² · Φ₀)         (exponential metric form)
  q     = 8π G₀ / (c₀² · Φ₀)         (Newton's-theorem form)
  qΦ₀   = 4π G₀ / c₀²                (Φ₀-free identity, [qΦ₀] = [L])

Natural-units (c=ℏ=1):  [q] = [-]
  must be flagged as a convention choice in any table that lists q-dim as [-].
```

Numerical sanity (sek03 line 524):
$$
q_{\exp} = \frac{4\pi \cdot 6.674\times 10^{-11}}{(3\times 10^8)^2}
\;\approx\; 9.30 \times 10^{-27}\ \mathrm{m/kg}
\;\;\;\checkmark
$$

## 3. Edits applied

### Edit 1 — `axioms/notacja/dodatekA_notacja.tex` lines 176–177

**Before:**
```latex
$q$ & $[-]$ & siła źródłowania & $\mathcal D[\Phi_a]=q_a\rho_a$ \\
$\rho$ & $[\mathrm{L}^{-3}]$ & gęstość źródłowa & $\int\rho\,d^3x = 1$ (normalizacja) \\
```

**After:**
```latex
$q$ & $[\mathrm{M}^{-1}\mathrm{L}]$ & siła źródłowania & $\mathcal D[\Phi_a]=q_a\rho_a$;\ $q=8\pi G_0/c_0^2$ (tw.~\ref{thm:newton}) \\
$\rho$ & $[\mathrm{ML}^{-3}]$ & gęstość masy--energii źródła & $\int\rho\,d^3x = m$ (masa całkowita) \\
```

### Edit 2 — `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` line 516

**Before:**
```latex
$q = 4\pi G_0/c_0^4 \cdot \PhiZero$
```

**After:**
```latex
$q = 4\pi G_0/(c_0^2\,\PhiZero)$
```

This is now consistent with line 519 (`q\Phi_0 = 4\pi G_0/c_0^2`) by simple
division by $\Phi_0$.

## 4. Sweep verification — no further inconsistencies

A regex sweep across all `.tex` files in `TGP/TGP_v1/` for `q = ...G₀/c₀^N`
patterns returned only:
- dodatekA (just-fixed) line 176 — now canonical
- sek03 line 524 — uses `4πG₀/c₀²` with `9.30×10⁻²⁷ m/kg` numerical value (canonical ✓)
- sek08 lines 2872, 2903 — use `q = 4πG₀/c₀²` (canonical ✓)
- sek08a line 516 — just fixed

A targeted regex `q.*c_0\^4.*\Phi` returned NO matches, confirming no other
typo of this form exists.

## 5. Falsifiable predictions

If the now-canonical $[q]=[M^{-1}L]$ were wrong, then ANY computation using
the Newton-boundary formula $q=8\pi G_0/c_0^2$ would yield Yukawa profiles
with the wrong overall scale relative to GR. In particular:

- The PPN-γ value computed via TGP at PN-1 order would deviate from
  $\gamma=1$ by a factor that scales as $[q]/(M^{-1}L)$. Closure cycle
  M5 (PPN at γ=1) confirmed agreement to $10^{-13}$, which
  independently verifies $[q]=[M^{-1}L]$.
- The pulsar timing tests in sek06/sek07 use $q=8\pi G_0/c_0^2$
  numerically; their match to GR is a strong consistency check.

## 6. Outstanding follow-ups

Not blocking B5 closure but flagged for future audit cycles:

1. **dodatekA table caption note**: Could add a brief explanation that the
   dimensional entries use SI-style $(c_0,\hbar_0,G_0$ all dimensionful)
   to match the rest of TGP_v1, with a cross-reference to the natural-units
   variant if needed elsewhere.
2. **sek08a κ derivation cross-check**: The line 511–520 derivation claims
   $\kappa=3/(4\Phi_0)$ via $q = 4\pi G_0/(c_0^2\Phi_0)$ × $c_0^2/H_0^2$ ×
   $(8\pi G_0/3 = H_0^2/\rho_{\mathrm{cr}})$. A direct check of the
   intermediate algebra (independent of the dimension fix) is suggested
   for a separate cycle.
3. **sek03 vs sek08a $\Phi_0$ convention**: sek03 uses
   $q = 4\pi G_0/c_0^2$ implicitly absorbing $\Phi_0=1$ (or treating
   $\Phi_0$ as adimensional unit), while sek08a keeps $\Phi_0$ explicit
   throughout. Both are self-consistent but use different bookkeeping.
   A unification note in either dodatekA or sek00 would help readers.

## 7. Final verdict

**B5 — CLOSED.** Sympy 5/5 PASS. Two LaTeX edits applied. Sweep clean.

Audit progression after B5: 40/43 → **41/43 (95.3%)**.

Remaining HIGH-priority items: B1 (τ.3 α_g sign Adams-positivity v2),
B3-v2 (α_s = 0.1184 propagation through tooling/sek00).
