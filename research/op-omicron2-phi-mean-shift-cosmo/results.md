---
title: "omicron.2 results — TGP cosmological Phi_0(t) tracking — NULL Hubble tension"
date: 2026-05-03
parent: "[[README.md]]"
type: results
status: STAGE_1_NULL — Stage 0 magnitude was incorrect; tension impact actually ~0.6%
verdict: "TGP does NOT solve Hubble tension via Phi_0 evolution mechanism (Stage 1 confirms)"
key_finding: |
  Stage 1 (re-shoot V_0 + D_A integration) ujawnia że Stage 0 conclusion był strukturalnie
  WRONG. Stage 0 użył formula |dH/H| = 0.5·|dL/L|·Omega_L (LDE-style impact today). Hubble
  tension wymaga EDE-style mechanism: H(recomb) higher → smaller r_s → higher H_0 inferred.
  W TGP, pomimo że Lambda_recomb/Lambda_today = 35,319, total energy at recomb dominated by
  matter (Lambda_recomb / rho_total_recomb = 4.4e-5), więc H(recomb) shift to 0.0022%, dając
  finalnie dH/H ~ 0.6% via D_A integration — 14× za mało.
verification_status: STAGE_1_DONE_NULL
verification_required:
  - stage_1_re_shoot_V0  ✅ DONE - NULL verdict
  - stage_2_DESI_w_z_compatibility  ABANDONED (Stage 1 nullified premise)
  - stage_3_CMB_and_BBN_check  ABANDONED
  - stage_4_independent_replication  ABANDONED
publication_status: ABANDONED — TGP NOT Hubble tension solver via this mechanism
tags:
  - TGP
  - omicron2
  - cosmology
  - Hubble-tension
  - NULL_RESULT
  - Phi_0-cosmological-evolution
  - sek01-sek08a-reconciliation
  - Stage_0_magnitude_error
---

# omicron.2 results — TGP cosmological Phi_0(t) tracking — NULL Hubble tension

> **STATUS Stage 1 NULL 2026-05-03.** Stage 0 conclusion (107% Hubble tension match)
> okazał się **strukturalnie błędny** w magnitude. Stage 1 ujawnił właściwy
> mechanism (D_A integration) i dał ~0.6% — 14× za mało.
>
> **Bottom line:** TGP NIE rozwiązuje Hubble tension przez Φ_0(t) tracking matter.
> M10.5 verdict efektywnie reinforced (różnymi tools, but same conclusion).

## 1. STAGE 1 NULL VERDICT — w jednym zdaniu

> **Stage 0 magnitude formula był WRONG.** Stage 1 D_A integration shows actual
> dH/H ≈ 0.6% (5/5 ICs), nie 8.9%. TGP **NIE** rozwiązuje Hubble tension.

## 2. Pełna historia (Stage 0 → Stage 1)

### Stage 0 (zapisany 2026-05-03 wczesny popołudnie):

```
s_natural (precise) = 18·Ω_DE0 / (Φ_eff²·g̃) = 2.07×10⁻²
|ΔH/H| (Stage 0 formula) = 0.5·|ΔΛ/Λ|·Ω_Λ = 8.93%
Coverage: 107% (overshoot 7%)
```

**EXCITING — but ultimately misleading.**

### Stage 1 (2026-05-03 wieczór) — RE-SHOOT V_0 + proper D_A integration:

```
5 initial conditions tested z proper V_0 shooting (Ω_DE0(today) = 0.685 enforced)

Direction: ALL EDE-like ✓ (Λ_recomb > Λ_today)
Magnitude: ALL ~0.6% via D_A ratio (NOT 8.9%)
Best match: 7.3% of required 8.37%
```

**NULL VERDICT — TGP nie pokrywa Hubble tension.**

## 2. Sequence wydarzeń (jak doszło do wyniku)

### 2.1 Background — initial drift audit (2026-05-03)

User zapytał o stan kosmologii TGP wobec DESI. Drift-check 4 plików kosmologicznych
ujawnił sek05/sek08a inkonsystencja co do Φ_0:
- **sek01**: "Φ_0 jest generowane przez kosmologiczną materię"
- **sek08a**: Φ_0 = vacuum minimum, niezależne od matter

### 2.2 Quick check (Opcja A) — confirms M10.5 ~10⁻⁸

Linear variance + Buchert backreaction → **2×10⁻⁸**. Gap 8 dekad do tension.
M10.5 verdict POŚREDNIO confirmed.

### 2.3 Audit (13 candidates) — identifies kandydat K

Centralna sprzeczność sek01 ↔ sek08a sugeruje że **matter source term w Φ-EOM**
mógł być pominięty. M10.5 obliczyło variance, NIE mean shift z source.

### 2.4 K-test FRW (s = 2.5e-3 naive) — surprise

Naive s gave **2.0×10⁻³** (NOT 10⁻⁸). 6 rzędów wielkości większe niż M10.5.
Gap do tension: ~factor 40.

### 2.5 Z-test precise (8.38× K_geo amplification) — KEY RESULT

Pełna derywacja z sek08a + closure cascade → s_precise = 2.07×10⁻². FRW: **8.93%**.

## 3. Quantitative result

### 3.1 Closed-form s_natural

```
s_de2 = 18 · Ω_DE0 / (Φ_eff² · g̃)
      = 1 / (2 · Φ_eff · g̃)              [alt form, using Ω_L = Φ_0/36]
      = 18 · Ω_DE0 / (64π² · g̃³)         [full cascade form, Φ_eff = 8π·g̃]
```

| Quantity | Value | Source |
|---|---:|---|
| Ω_DE0 | 0.6847 | Planck 2018 |
| Φ_eff | 24.6302 | γ.1 closure: (10/3)·e² |
| g̃ | 0.98003 | δ.1+δ.2: 5·e²/(12π), N_f=5 |
| **s_de2 (precise)** | **2.073×10⁻²** | this calc |

### 3.2 K_geo derivation chain

Crucial: K_geo MUSI być 0.119 (nie 1!) by V_0 = Ω_DE0 spełnione pod T-Λ:

```
T-Λ closure:  γ = g̃·M_Pl²·H_0²
de2 V_0:      V_0 = γ/(12·K_geo·H_0²) = Ω_DE0
⇒  K_geo = g̃·M_Pl² / (12·Ω_DE0) = 0.1193  (Planck units, M_Pl=1)
```

**To jest non-trivial: previous M10.5 + de2 used K=1 (canonical), ale
sek08a z T-Λ closure forces K_geo = 0.119.**

### 3.3 FRW solver result

| s | ψ_today | ψ_recomb | Δψ | ΔΛ/Λ | \|ΔH/H\| |
|---:|---:|---:|---:|---:|---:|
| 0 (no source) | 0.997 | 0.999 | -1.6e-3 | -3.5e-5 | 1.2e-5 |
| 2.5e-3 (naive) | 0.968 | 0.997 | -2.9e-2 | -6e-3 | 2.0e-3 |
| **2.07e-2 (PRECISE)** | **0.783** | **0.983** | **-0.20** | **-0.26** | **🔴 0.0893** |
| 4.1e-2 (2× precise) | 0.619 | 0.968 | -0.35 | -0.96 | 0.328 |
| 1e-2 (0.5× precise) | 0.881 | 0.991 | -0.11 | -0.077 | 0.026 |

Required: 0.0837 (8.37%). **Precise s daje 0.0893 = 107% pokrycie.**

## 4. Caveats — to NIE jest jeszcze final

### 4.1 ⚠️ w(z) prediction shifts dramatically

Pre-cascade prediction (M10.1): w₀ = -1.000, w_a = 0.000 (essentially ΛCDM).

**Post-cascade z precise s:**
- ψ_today = 0.78 (significant slow-roll, NIE ψ ≈ 1 frozen)
- φ̇ ≈ 0.20 (znaczące, NIE ~0)
- KE_today = 0.02; PE_today = 0.54
- **w_DE_today ≈ -0.93** (NIE -1.0)

To jest **sprzeczne z M10.1 numerical verification** (w(z) ≈ -1 w 10⁻⁴).
Z precise s, M10.1 conclusion nie zachodzi.

**Zgodność z DESI DR2:** w₀ = -0.93 jest **bliżej -1 niż DESI** (-0.75 ± 0.10).
Może być w 1-2σ gap. Nie sfalsyfikowane.

### 4.2 ⚠️ Ω_DE0 today ≠ 0.685 z precise s

W skrypcie używaliśmy V_0 = 0.685 (jako de2 baseline). Ale z source pushing,
prawdziwe Ω_DE0 today wynosi ~0.56 (KE+PE), nie 0.685.

**Need to re-shoot V_0** by dać Ω_DE0 = 0.685 today AFTER 26% Λ-decay.

### 4.3 ⚠️ Recombination Λ-eff = 1.27× today's

ψ_recomb = 0.983 vs ψ_today = 0.78. V(0.983) = 0.998, V(0.78) = 0.788.
Λ_eff(recomb)/Λ_eff(today) = 0.998/0.788 = 1.27.

**Implications dla CMB:**
- Sound horizon at recomb depends on H(z=1100), które zależy od Λ(recomb)
- ISW effect zmieniony
- Acoustic peaks mogą shiftować

**Need to verify Planck CMB fit nie psuje się.**

### 4.4 ⚠️ Single-pass calculation, no cross-check

Z_precise_derivation.py wyniki SĄ z **one calculation pass**. Niepewność:
- 8.38× amplification z K_geo — czy poprawnie wyprowadzone z T-Λ?
- 4π factor w Newton limit q·Φ_0 = 4π·G/c² — czy to jest TGP convention?
- Phi_eff = 8π·g̃ algebraic identity — czy spójne ze sek00?

**Independent replication required przed publikacją.**

### 4.5 ⚠️ Other M10.x verdicts shift

Predykcje M10.x (M10.1 w(z), M10.2 inflation, M10.3 propagator, M10.4 CMB) były
oparte na de2-style slow-roll z ψ≈1 frozen. **Wszystkie wymagają re-evaluation
w post-cascade scenario** z precise s.

## 5. Strategic implications

### 5.1 Pre-cascade TGP (M10.5 era, < 2026-04-26)

- Λ ≈ const cosmologicznie
- TGP-DE indistinguishable from ΛCDM
- NIE rozwiązuje H₀ tension
- Honest scope: galaxy-scale + structural DE
- Φ_0 = parametr dopasowania

### 5.2 Post-cascade TGP (Z-test era, ≥ 2026-05-03)

- Λ ewoluuje cosmologicznie (decay 26% od recomb)
- TGP-DE = quintessence-like (w(z) = -0.93)
- **POTENTIALLY rozwiązuje H₀ tension strukturalnie** (107% match)
- Scope rozszerza się do **DE evolution + H₀ tension**
- Φ_0 = algebraic prediction, cosmologically evolving

To jest **fundamentalny shift** w statusie kosmologicznym TGP.

### 5.3 Co potwierdza, co reverse-uje

| M10.5 / M10.R / M10.1 claim | Status post-Z |
|---|---|
| "B_ψ/H_0² ~ 10⁻⁸" (M10.5.3) | ✅ correct dla Buchert variance kanału |
| "Gap 7 rzędów do H_0 tension" (M10.5.6) | 🔴 **WRONG** — gap zamknięty z source term |
| "TGP nie jest H_0 tension solver" (M10.R.6) | 🔴 **REVISED** — może być solver |
| "w(z) ≈ -1 numerycznie" (M10.1.3) | 🔴 **WRONG** — w(z) = -0.93 z source term |
| "CPL (w₀, w_a) ≈ (-1, 0)" (M10.1.4) | 🔴 **WRONG** — (-0.93, ?) potencjalnie matchuje DESI |
| "Honest scope: galaxy NIE cosmology" (M10.R.6) | 🔴 **REVISED** — TGP pokrywa cosmology |
| "Path B σ_ab" (closure_2026-04-26) | ✅ niezmienione |
| "T-Λ Φ_eq = ℏH_0" (closure_2026-04-26) | ✅ niezmienione (ZASILA Z-test) |
| "γ.1 + δ.1 + δ.2 closures" (2026-05-02) | ✅ niezmienione (ZASILAJĄ Z-test) |

## 6. Verification roadmap (4 stages)

### Stage 1 — Re-shoot V_0 (~1-2h)

**Cel:** Sprawdzić czy z source term ON, Ω_DE0 = 0.685 wciąż osiągalny przez
shooting V_0_initial, **oraz** czy ΔH/H pozostaje ~8% z properly shooted V_0.

**Skrypt:** `stage1_reshoot_V0.py`

**Test passing condition:** Ω_DE0(today) = 0.685 ± 0.001 ze shooting, ΔH/H
pozostaje w przedziale [0.05, 0.12] (50-150% Hubble tension).

### Stage 2 — DESI w(z) compatibility (~2-3h)

**Cel:** Compute (w₀, w_a) CPL fit z full FRW history po Stage 1. Compare
z DESI DR2 `(-0.75, -0.90)` 3.1σ.

**Skrypt:** `stage2_DESI_compatibility.py`

**Test passing condition:** TGP CPL prediction w 1-2σ od DESI DR2 best-fit.
Akceptowalny: w₀ ∈ [-0.95, -0.7], w_a ∈ [-1.1, -0.5].

### Stage 3 — CMB + BBN check (~3-4h)

**Cel:** Sprawdzić czy Λ_eff(recomb) = 1.27× today's psuje:
- Planck CMB acoustic peaks (sound horizon, ISW)
- BBN He⁴ abundance (G_eff drift)

**Skrypt:** `stage3_CMB_BBN_check.py`

**Test passing condition:** CMB χ² nie pogarsza się więcej niż 5 vs LCDM.
BBN He⁴ w pre-existing M10.4 limits.

### Stage 4 — Independent replication (~1-2 dni)

**Cel:** Wykorzystując subagenta lub niezależną sesję, **niezależna derywacja
s_natural** z sek08a action, niezależna FRW integration. Cross-verify s = 2.07e-2.

**Akcja:** spawn fresh agent w nowej sesji, dac mu task "derive s_natural from
sek08a + T-Λ + UV.3", compare z naszym wynikiem.

**Test passing condition:** Independent reasults s ∈ [1.5e-2, 2.7e-2] (within ~30%).

## 7. Cross-references

### 7.1 Source code
- [[Z_precise_derivation.py]] — main derivation script
- [[Z_output.txt]] — full output
- [[omicron2_phi_mean_shift_check.py]] — Opcja A (linear/Buchert variance)
- [[omicron2_K_test.py]] — Opcja A naive coupling test
- [[quick_test.py]] — fast numerical scan

### 7.2 Background docs
- [[README.md]] — overview tego cyklu
- [[AUDIT_blind_spots.md]] — 13 candidates analysis
- [[ROADMAP.md]] — staged verification plan

### 7.3 Closure cascade (2026-04-26 → 2026-05-02)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] — T-Λ (γ = M_Pl² H_0² g̃)
- [[../op-uv3-phi0-renormalization/Phase2_results.md]] — Z_Φ = 14/3
- [[../op-gamma1-phi-eff-anchor-resolution/]] — Φ_eff = 8π g̃ = (10/3)·e²
- [[../op-delta1-g-tilde-derivation/]] — g̃ = N_f·e²/(12π)
- [[../op-delta2-Nf-derivation/]] — N_f = 5

### 7.4 Affected M10 cycle results (need re-evaluation)
- [[../op-cosmology-closure/M10_R_results.md]] — M10 synthesis (12 falsifiers)
- [[../op-cosmology-closure/M10_5_results.md]] — H_0/S_8 audit
- [[../op-cosmology-closure/M10_1_results.md]] — DE w(z) audit
- [[../op-cosmology-closure/M10_4_results.md]] — CMB safety
- [[../desi_dark_energy/README.md]] — Program P2 redirect

### 7.5 Sek01 ↔ sek08a reconciliation pending
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] §rem:tablica-stanow (linia 700-703): "Φ_0 generowane przez kosmologiczną materię"
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]: V'(Φ_0) = 0 vacuum

Z-test pokazuje że obie wykładni są **zgodne** strukturalnie — Φ_0 IS vacuum
minimum lokalnie, ALE matter source term forces cosmological evolution. **Brak
sprzeczności po dokładnym rachunku.**

## 8. Bottom line — ZAKTUALIZOWANY post Stage 1

> **TGP cosmology post-cascade NIE jest Hubble tension solver.**
>
> Stage 0 wstępny wynik (107% pokrycie) okazał się strukturalnie błędny —
> użyta formula `dH/H = 0.5·|dL/L|·Ω_L` mierzy fractional H change today,
> NIE EDE-like H boost at recombination który jest właściwym mechanism dla tension.
>
> Stage 1 z proper D_A integration: TGP daje **0.6% pokrycia tension** (vs
> required 8.4%). Wszystkie 5 testowanych initial conditions dają ten sam wynik
> w obrębie 0.1%. **Magnitude error w Stage 0 = 14×.**
>
> **Strukturalne wyjaśnienie:** TGP source term w Φ-EOM **istnieje** (potwierdza
> sek01 ontology), Φ rzeczywiście tracks ρ̄(t), ALE coupling jest "active" tylko
> w late universe gdzie ρ_m ~ ρ_DE. Przy recombinacji ρ_m_recomb >> Λ_recomb_TGP,
> więc TGP source nie boostuje H(recomb), nie zmniejsza r_s, nie shiftuje
> H_0_inferred via CMB.
>
> M10.5 verdict ("TGP NIE jest H_0 tension solver, gap structural") — efektywnie
> **REINFORCED** w Stage 1, choć przez inne mechanism niż M10.5 use.
> 
> User intuition (Φ_0 should track matter): **strukturalnie POTWIERDZONA**
> (sek01 ontology + matter source term faktycznie istnieje), ale **kwantytatywnie
> NIE rozwiązuje Hubble tension** w ramach minimalnego TGP coupling.

## 9. Lessons learned (this is also a result)

### 9.1 Stage 0 magnitude error postmortem

**Where Stage 0 went wrong:**
1. ✅ Identified matter source term w Φ-EOM (correct)
2. ✅ Derived precise s_natural z UV.3+γ.1+δ.1+δ.2+K_geo (correct algebra)
3. ✅ FRW solver gave large ψ shifts (real)
4. 🔴 **Used wrong tension impact formula:** `dH/H = 0.5·|dL/L|·Ω_L` is for
     LDE-style today's Lambda shift, NOT Hubble tension's EDE-style mechanism

**Right formula (Stage 1):**
```
dH/H_inferred ~ +(D_A_LCDM - D_A_TGP) / D_A_TGP
              = (1/D_A_ratio - 1)
```
where D_A integration carries actual H_recomb effect through to inferred H_0.

### 9.2 What Stage 0 missed

- **Ratio** `Λ_recomb/Λ_today = 35,319` looked huge BUT in absolute units
  `Λ_recomb / ρ_total_recomb = 4.4e-5` (negligible)
- **r_s** integral dominated by matter+radiation at recomb, where Λ-shift is irrelevant
- **EDE-like solutions** require Λ-like component **comparable to ρ_total** at recomb
  — TGP physically cannot do this without modifying ax:metric-coupling

### 9.3 What's NOT changed by Stage 1

- **Sek01 ontology:** Φ_0 generated by matter — still TRUE (qualitatively)
- **Matter source term in Φ-EOM:** structurally exists (sek08a)
- **ψ tracks ρ̄(t):** real (Stage 1 confirms ψ shifts ~20% cosmologically)
- **closure_2026-04-26 + UV.3 + γ.1 + δ.1 + δ.2:** all valid
- **w(z) prediction shift:** real (M10.1 conclusion w(z)≈-1 was based on de2 without
  matter source; with proper source, w_today ≈ -0.93). Independent issue.

## 10. Implications for M10 results revision

| M10 result | Pre-Stage_1 status | Post-Stage_1 status |
|---|---|---|
| M10.1: w(z) ≈ -1 | unchanged (de2 without source) | 🟡 **needs revision z proper source** (separate issue from Hubble tension) |
| M10.5: TGP NIE H_0 solver | 🟡 challenged by Stage 0 | ✅ **REINFORCED** by Stage 1 |
| M10.R: 12 falsifiers | unchanged | unchanged |
| sek01 ↔ sek08a reconciliation | sprzeczność | ✅ **resolved** — matter source IS in formalism, just doesn't solve tension |

## 11. Final verdict

🔴 **TGP NIE jest Hubble tension solver przez Φ_0(t) tracking mechanism.**

Stage 1 verification confirms M10.5 conclusion through independent route:
- M10.5: B_ψ/H_0² ~ 10⁻⁸ via Buchert variance
- Stage 1: dH/H_inferred ~ 0.6% via D_A integration
- Both confirm: TGP coupling structurally insufficient

🟢 **Sek01 ontology vs sek08a formalism:** RECONCILED. Matter source term
istnieje w sek08a, daje ψ tracking ρ̄(t), strukturalnie zgodne z sek01
"Φ_0 generated by matter". Ale magnitude (przez minimal coupling) jest
za mała by rozwiązać kosmologiczne napięcia.

🟡 **Affected M10 predictions:** Glob check needed (Stage 5 abbreviated).
Specifically:
- M10.1 w(z) prediction may shift from -1 to -0.93 (independent issue)
- DESI compatibility actually IMPROVES with -0.93 (closer to DR2 -0.75)

## 12. Co dalej

### Honest path forward:

**(A) Close omicron2 with NULL verdict** — record this as honest negative result.
Update sek05/sek08a/M10.5 z explicit nota: "matter source confirmed exists,
but cosmologically negligible for Hubble tension; M10.5 conclusion stands".

**(B) Pursue M10.1 revision** — separate cycle to compute w(z) with proper
source term. Result: w_today ≈ -0.93. Compare with DESI DR2 (-0.75 ± 0.10).
This is a more productive question than Hubble tension.

**(C) Look for non-minimal coupling** — would TGP with extra ρ-Φ coupling
beyond ax:metric-coupling solve tension? Speculative, requires modifying
TGP_FOUNDATIONS §4.

**Recommendation:** **(A) + (B)**. Close NULL on Hubble tension, pursue
w(z) revision separately as it's well-motivated by Stage 1 finding.

---

*results.md updated 2026-05-03 wieczór po Stage 1 NULL verdict.
Stage 0 conclusion REVERSED. M10.5 effectively reinforced.
Publication blocked permanently for Hubble tension claim.
w(z) revision opens as separate research direction.*
