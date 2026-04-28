---
title: "M10.1 setup — FRW dark energy w(z) audit (de2 audit)"
date: 2026-04-26
cycle: M10
phase: M10.1
status: SETUP
predecessor: "[[M10_0_drift_audit.md]]"
audit_target: "[[../desi_dark_energy/de2_tgp_frw_evolution.py]]"
parent: "[[M10_program.md]]"
tags:
  - TGP
  - M10
  - dark-energy
  - audit-setup
---

# M10.1 — FRW dark energy w(z) audit setup

> **Cel:** weryfikacja że istniejący draft [[../desi_dark_energy/de2_tgp_frw_evolution.py]] daje `w(z) ≥ -1` strukturalnie z **full sek08a kinetic** `K(φ)=K_geo·φ⁴` (NIE tylko canonical K=1), oraz że CPL fit `(w_0, w_a)` near-ΛCDM utrzymuje się.
>
> **Drift do verify:** de2 używa canonical K=1 (line 100), ale sek08a ma K(φ)=K_geo·φ⁴ non-canonical. Near vacuum ψ≈1 to sub-leading correction (K≈K_geo·1=const), ale powinno być sprawdzone.

---

## 1. Cel & predykcje

### 1.1 Cel sub-cyklu

Z drift audit M10.0: de2 jest YELLOW (canonical K=1 vs sek08a K=K_geo·φ⁴). Audit weryfikuje:
1. Strukturalnie: **bound `w ≥ -1`** — czy zachowany przy non-canonical K?
2. Numerycznie: **w(z)** profile near-ΛCDM (`w_0 ≈ -1`, `w_a ≈ 0`)?
3. Comparison: **CPL fit do DESI DR1** — TGP near-ΛCDM, DESI w_0=-0.45, w_a=-1.79 → falsifiability?

### 1.2 Predykcje TGP

Z sek08a + de2:
- `V(ψ) = (β/3)ψ³ - (γ/4)ψ⁴`, `β=γ`, `V''(1)=-β` slow-roll max
- Slow-roll inflation-like dynamics dla ψ ≠ 1 z Hubble damping
- **Strukturalnie:** `w + 1 = (½K(φ)φ̇²) / ρ_ψ ≥ 0` — bo K(φ)>0 zawsze (positive definite)
- **TGP unique:** brak phantom crossing (`w < -1`) bez non-canonical kinetic, którego SEK08A NIE MA w sensie ghost (K=K_geo·φ⁴ > 0)

### 1.3 Cross-checks z closures

| Closure | Use w M10.1 |
|---------|-------------|
| T-Λ (Φ_eq=H_0, ρ_vac=M_Pl²H_0²/12) | `V_0` shooting parameter z `Ω_DE0=0.685` |
| sek08a (V(ψ)=(β/3)ψ³-(γ/4)ψ⁴) | potential form |
| sek08a (K(φ)=K_geo·φ⁴) | non-canonical kinetic — TEST DRIFT |
| M9.1'' (hyperbolic metric) | not direct (FRW limit) |

---

## 2. Setup matematyczny

### 2.1 Akcja TGP (sek08a)

```
S_TGP = ∫ d⁴x √(-g_eff) [ ½ K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ) - (q/Φ_0) φ ρ ]
K(φ)  = K_geo φ⁴
V(φ)  = (β/3)φ³ - (γ/4)φ⁴,    β = γ
```

### 2.2 FRW EOM (homogeneous ψ(t))

W FRW (-,+,+,+) signature, dla homogeneous `ψ(t)`:

**Canonical kinetic K=1 (de2 default):**
```
ψ̈ + 3H ψ̇ + V'(ψ) = 0
ρ_ψ = ½ ψ̇² + V(ψ)
p_ψ = ½ ψ̇² - V(ψ)
w_ψ = p_ψ / ρ_ψ
```

**Non-canonical K(φ)=K_geo·φ⁴ (sek08a):**

Z action density `L = √(-g_eff)[½K(φ)g^μν∂_μφ∂_νφ - V(φ)]`:
- Kinetic energy density: `ρ_K = ½ K(ψ) ψ̇² = ½ K_geo ψ⁴ ψ̇²`
- Total: `ρ_ψ = ½ K_geo ψ⁴ ψ̇² + V(ψ)`
- Pressure: `p_ψ = ½ K_geo ψ⁴ ψ̇² - V(ψ)`
- EOM: `K(ψ) ψ̈ + ½ K'(ψ) ψ̇² + 3H K(ψ) ψ̇ + V'(ψ) = 0`
- Z `K'(ψ) = 4 K_geo ψ³`: `K_geo ψ⁴ ψ̈ + 2 K_geo ψ³ ψ̇² + 3H K_geo ψ⁴ ψ̇ + V'(ψ) = 0`
- Dzieląc przez `K_geo ψ⁴`: `ψ̈ + 2 ψ̇²/ψ + 3H ψ̇ + V'(ψ)/(K_geo ψ⁴) = 0`

**Kluczowa różnica (vs canonical):**
- Dodatkowy term `+2ψ̇²/ψ` (geometric friction)
- Effective force `V'(ψ)/(K_geo ψ⁴)` skalowana przez `ψ⁴`

**Bound `w ≥ -1`:**
- `w + 1 = 2·(½ K(φ)ψ̇²) / ρ_ψ = K(φ)ψ̇² / ρ_ψ`
- Dla `K(φ) = K_geo φ⁴ > 0`: `w + 1 ≥ 0` ZACHOWANE dla `ψ > 0`.
- **Strukturalna konkluzja:** bound `w≥-1` IS robust z non-canonical kinetic, **provided K(ψ) > 0**.

### 2.3 Friedmann equation

```
H² = (8πG/3)(ρ_m + ρ_r + ρ_ψ)
```

W jednostkach H_0=1, 8πG/3=1:
```
H²(a) = Ω_m0/a³ + Ω_r0/a⁴ + ρ_ψ
ρ_ψ(a) = ½ K(ψ) ψ̇² + V_0·V(ψ)   (K(ψ)=K_geo·ψ⁴ for non-canonical)
```

z `V_0` shooting parameter tak żeby `ρ_ψ(a=1) = Ω_DE0 = 0.685`.

### 2.4 CPL parametrization

DESI DR1 raportuje `(w_0, w_a)` w CPL form:
```
w(a) = w_0 + w_a (1 - a)
w(z) = w_0 + w_a · z/(1+z)
```

DESI DR1 (2024) constraints:
- `w_0 = -0.45 ± 0.21`
- `w_a = -1.79 ± 0.65`

ΛCDM: `w_0=-1, w_a=0`.

---

## 3. Plan testów (5 sub-testów)

### M10.1.1 — Action structure verification (symboliczny)

**Cel:** weryfikacja sek08a V form, V'(1)=0, V''(1)=-β slow-roll max, K(φ)=K_geo·φ⁴ poprawna.

**Metoda:** sympy:
- `V_sym = (β/3)ψ³ - (γ/4)ψ⁴`
- Substitute β=γ
- Compute V_sym, V'_sym, V''_sym at ψ=1
- Verify: V(1)=β/12, V'(1)=0, V''(1)=-β

**Sub-tests:**
- (a) V'(1) = 0 (vacuum cond) — symbolic
- (b) V''(1) = -β (slow-roll max) — symbolic
- (c) V(1) = β/12 > 0 (residual vacuum energy = Λ source via T-Λ) — symbolic
- (d) K_geo·φ⁴ > 0 dla ψ>0 (positivity guarantees w≥-1) — symbolic

### M10.1.2 — Bound `w ≥ -1` z non-canonical kinetic

**Cel:** weryfikacja że bound `w ≥ -1` jest robust przy `K(φ)=K_geo·φ⁴`.

**Metoda:** sympy:
- `w + 1 = K(ψ)ψ̇² / ρ_ψ`
- ρ_ψ = ½K(ψ)ψ̇² + V(ψ)
- Show: `w + 1 ≥ 0` iff `K(ψ) > 0` (independent of V sign)

**Sub-tests:**
- (a) Symbolic `w + 1 = K(ψ)ψ̇²/ρ_ψ` (sympy)
- (b) `K(ψ) = K_geo·ψ⁴ > 0` for ψ>0 — proven
- (c) Implication: `w + 1 ≥ 0` strukturalnie — proven
- (d) Compare with canonical K=1: same form, same bound

### M10.1.3 — Numerical FRW evolution: canonical vs non-canonical

**Cel:** porównanie w(z) z canonical (K=1) i non-canonical (K=K_geo·ψ⁴).

**Metoda:** dwie integracje:
- IC: `ψ_i = 1 + δ`, `ψ̇_i = 0`, z `δ ∈ {1e-4, 1e-3, 1e-2}` at `a_i = 0.01`
- ODE solve: scipy.integrate.solve_ivp, rtol=1e-9
- Extract w(z) at z={0, 0.5, 1, 2, 5}
- For non-canonical: include `K_geo·ψ⁴` factors in EOM and ρ_ψ, p_ψ
- Use K_geo = 1 (sub-leading shift; absolute value irrelevant for w(z))

**Sub-tests:**
- (a) Canonical w(z=0) ≈ -1 (frozen) — quantitative
- (b) Non-canonical w(z=0) ≈ -1 (frozen) — same
- (c) Difference |w_canonical - w_non-canonical| < 0.01 (sub-leading) — verify
- (d) Both: w(z) ≥ -1 for all z, all δ — verified
- (e) Hubble damping: ψ slow-rolls but ψ remains close to 1 today

### M10.1.4 — CPL projection (w_0, w_a)

**Cel:** fit TGP w(z) do CPL parametryzacji.

**Metoda:**
- Generate w(z) na siatce z∈[0, 2] (40 punktów)
- Linear lstsq fit: `w(a) = w_0 + w_a (1-a)` w bazie `[1, 1-a]`
- Compare: TGP `(w_0, w_a)` vs ΛCDM `(-1, 0)` vs DESI DR1 `(-0.45, -1.79)`

**Sub-tests:**
- (a) TGP `|w_0 + 1| < 0.1` (near-ΛCDM)
- (b) TGP `|w_a| < 0.5` (small evolution)
- (c) Distance in (w_0, w_a) plane: |TGP - ΛCDM| < |TGP - DESI|
- (d) Non-canonical and canonical CPL fits agree within 1%

### M10.1.5 — DESI DR1 falsifiability

**Cel:** quantify falsifiability of TGP vs DESI DR1.

**Metoda:**
- Compute χ² distance: TGP (w_0_TGP, w_a_TGP) vs DESI mean (-0.45, -1.79) z DESI covariance (assumed σ_w0=0.21, σ_wa=0.65, ρ=-0.5 typical)
- Convert χ² → significance σ
- Statement: "TGP excluded at >Nσ jeśli DESI DR1 values hold"

**Sub-tests:**
- (a) χ²_TGP_vs_DESI > χ²_3sigma → TGP falsifiable now (3σ)
- (b) Cross-check vs ΛCDM: same conclusion (ΛCDM też 3σ od DESI DR1)
- (c) Falsifiability vs DESI DR2/DR3 (oczekiwane improved precision)
- (d) Honest: DESI DR1 may shift; final verdict needs DR2

### M10.1.6 — T-Λ closure consistency check (cross-validation)

**Cel:** spójność z [[../closure_2026-04-26/Lambda_from_Phi0/results.md]].

**Metoda:**
- T-Λ closure: `ρ_vac = M_Pl² H_0² / 12 = β·Φ_0² / 12` (analog)
- Z M10.1: `V_0` shoot parameter ≈ Ω_DE0 = 0.685
- Verify: V_0 ≈ V(1) ratio matches T-Λ prediction

**Sub-tests:**
- (a) `V_0 / V(1) ≈ Ω_DE0 / Ω_total` — check ratio
- (b) `V(1) = β/12` — matches T-Λ ρ_vac form
- (c) Cross-validation z [[../closure_2026-04-26/Lambda_from_Phi0/results.md]]

---

## 4. Numerical realization (m10_1_de.py outline)

```python
# Top: sympy symbolic checks
def test_M10_1_1():  # V form
def test_M10_1_2():  # w >= -1 bound

# Middle: numerical FRW
def integrate_TGP(K_form='canonical' | 'non_canonical', delta=1e-3):
    """ODE: y = [a, psi, u]
       du/dt = -3H u - V'/(K(psi))  (canonical: K=1; non: K=K_geo psi^4)
       Add geometric friction term for non-canonical
    """

def test_M10_1_3(...):  # canonical vs non-canonical
def test_M10_1_4(...):  # CPL fit
def test_M10_1_5(...):  # DESI DR1

# Bottom:
def test_M10_1_6(...):  # T-Lambda cross-check
def main():  # 6 tests aggregate, expect 6/6 PASS
```

**Bezpieczne defaults:**
- `K_geo = 1` (units; absolute value irrelevant)
- `Φ_0 = 1` (normalized)
- `β = γ = 1` (vacuum cond, normalized so V(1)=1/12 absorbed in V_0 shoot)
- DESI DR1: `(w_0, w_a) = (-0.45, -1.79) ± (0.21, 0.65)`

**Output:** `m10_1_de.txt` z zapisem testów + werdyktem 6/6 PASS expected.

---

## 5. Falsifiable predictions

1. **TGP DE structurally `w(z) ≥ -1`:** any phantom crossing at >3σ falsifies TGP.
2. **TGP CPL fit:** `(w_0, w_a) ≈ (-1, 0)` (near-ΛCDM). DESI DR2/DR3 confirming `(w_0, w_a) = (-0.45, -1.79)` at >5σ falsifies TGP.
3. **Non-canonical kinetic K=φ⁴:** sub-leading correction to canonical near vacuum ψ=1 (% level), but should be detectable in precision DESI DR3 if K_geo·φ⁴ different from 1 by O(0.1%).

---

## 6. Drift-check matrix (post-execution)

| Constraint | Test |
|------------|------|
| sek08a `V(ψ)=(β/3)ψ³-(γ/4)ψ⁴` | M10.1.1 verifies V form |
| sek08a `K(φ)=K_geo·φ⁴` non-canonical | M10.1.3 tests both K=1 and K=φ⁴ |
| `V''(1)=-β` slow-roll max | M10.1.1 (c), M10.1.3 (slow-roll dynamics) |
| Bound `w≥-1` | M10.1.2 symbolic, M10.1.3 numerical |
| T-Λ closure consistency | M10.1.6 cross-validation |

---

## 7. Następne (post M10.1)

Jeśli 6/6 PASS:
- M10.1 ZAMKNIĘTE
- Przejście do M10.2 (ex261 inflation audit)

Jeśli FAIL (np. K=φ⁴ daje istotną różnicę):
- Investigate; może wymagać M10.1.b (rework) z full hyperbolic metric M9.1''.

---

*M10.1 setup completed 2026-04-26. Ready for execution.*
