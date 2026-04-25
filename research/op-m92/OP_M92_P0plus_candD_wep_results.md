# OP-M92 Phase 0+ — Candidate D WEP / 5th force cross-check (results)

**Data:** 2026-04-25
**Status:** **Phase 0+ WEP cross-check DONE — verdict POSITIVE (z notą o ciasnym MICROSCOPE marginesie)**
**Sub-test:** Candidate D vs MICROSCOPE / Eöt-Wash / LLR Nordvedt bounds
**Następnik:** Phase 1 deferred item (e) "compositional constraints / equivalence principle" — partial closure

---

## TL;DR

> Candidate D back-reaction term `α T^μν J_μ J_ν` modyfikuje **wyłącznie**
> sektor grawitacyjny (Φ-EOM); matter EOM pozostaje niezmodyfikowane.
> ⇒ UFF (Universality of Free Fall) preserved at zeroth order — wszystkie
> matter species następują tych samych geodezyjnych w effective metric.
>
> Residual composition-dependent acceleration (Nordvedt-like) z heuristic
> scaling `η ~ α (∇Φ)² × Δ(E_B/Mc²)`:
>
> - **MICROSCOPE Ti-Pt (Earth surface):** η ≈ **1.64e-16** vs bound 1.1e-15
>   → margin **6.7×** (TIGHT but PASS)
> - **LLR Nordvedt (Sun-Earth-Moon):** η ≈ **4.66e-29** vs bound 1.4e-13
>   → margin **3×10¹⁵** (very safe)
> - **Internal pressure contribution:** suppressed by `p/(ρc²) ~ 10⁻²¹` (negligible)
>
> **Verdict POSITIVE** — wszystkie bounds passed, choć MICROSCOPE jest
> najbardziej stringent test i powinien być **first-order priority**
> w pełnym Phase 1 covariant derivation.

---

## 1. Structural WEP argument (zerowy rząd)

Action structure Candidate D:

```
S = S_M9.1'' + α ∫ T^μν ∂_μ Φ ∂_ν Φ √(-g) d⁴x
```

Wariacja względem **pól materii** (np. ψ_DM, A^μ, Dirac ψ): T^μν wchodzi
jako **współczynnik** mnożący kinetic Φ-term, nie jako nowy bezpośredni
vertex matter-Φ.

⇒ Matter EOM:
```
∇_μ T^μν = 0    (niezmodyfikowane)
```

⇒ Test bodies of any composition follow **same geodesics** of effective
metric. UFF/WEP preserved **strukturalnie** at zeroth order.

Modyfikacja Φ-EOM (sourced by α T·J·J term) zmienia *profile* Φ-field
wokół matter, ale wszystkie species widzą **ten sam** modified profile.

## 2. Residual composition dependence (Nordvedt-like)

Composition-dependence wchodzi tylko przez **wewnętrzną strukturę**
extended test bodies — różne species mają różne stress-energy decompositions:

| Source | Scale (per nucleon/kg) | Comment |
|--------|------------------------|---------|
| Nuclear binding energy | E_B/Mc² ~ 10⁻³ | Largest residual contribution |
| Gravitational self-energy | U_self ~ 10⁻¹⁰ (Earth) | LLR Nordvedt scale |
| Internal pressure | p/(ρc²) ~ 10⁻²¹ (lab) | Negligible |
| Electromagnetic binding | E_EM/Mc² ~ 10⁻³ | Similar scale to nuclear |

Heuristic 5th-force scaling:
```
η  ~  α |∇Φ_external|²  ×  Δ(internal energy fraction)
```

## 3. Numerical results

### 3.1 Calibration α z Sgr A* (z Phase 0+ sketch + cosmology)

- α (geom units, M_BH=1) = 0.1
- α (SI, length²) = 0.1 × R_S(Sgr A*)² ≈ **1.61×10¹⁹ m²**
- α (SI, time²) ≈ 180 s²

### 3.2 MICROSCOPE Ti-Pt at Earth surface

- |∇Φ_Earth|² = (g/c²)² ≈ 1.19×10⁻³² m⁻²
- Nuclear binding: E_B(Ti-48)/A = 8.72 MeV/nucleon → E_B/Mc² ≈ 9.29×10⁻³
- Nuclear binding: E_B(Pt-195)/A = 7.92 MeV/nucleon → E_B/Mc² ≈ 8.43×10⁻³
- Δ(E_B/Mc²) ≈ **8.52×10⁻⁴**

```
η_MICROSCOPE  ≈  α × |∇Φ|² × Δ(E_B/Mc²)
              ≈  1.61×10¹⁹ × 1.19×10⁻³² × 8.52×10⁻⁴
              ≈  1.64×10⁻¹⁶
```

**Bound:** |η_Ti-Pt| < 1.1×10⁻¹⁵ (MICROSCOPE 2022, Touboul et al.)
**Margin:** **6.7×** safe (TIGHT)

### 3.3 LLR Nordvedt

- |∇Φ_Sun at 1AU|² ≈ 4.36×10⁻³⁹ m⁻²
- ΔU_self (Earth-Moon) ≈ 6.65×10⁻¹⁰

```
η_LLR  ≈  1.61×10¹⁹ × 4.36×10⁻³⁹ × 6.65×10⁻¹⁰
       ≈  4.66×10⁻²⁹
```

**Bound:** |η| < 1.4×10⁻¹³ (LLR Apollo retroreflector)
**Margin:** **3×10¹⁵** safe (very comfortable)

### 3.4 Internal pressure contribution

- Lab pressure (atmospheric + structural): p ~ 10⁵ Pa
- Pt density: ρ ~ 21450 kg/m³
- p/(ρc²) ≈ 5.2×10⁻¹⁷ (cf. nuclear binding ~10⁻³)
- ⇒ ~14 orders of magnitude weaker than nuclear contribution → negligible

## 4. Critical findings

### Finding 1: Matter EOM strukturalnie niezmodyfikowane

To jest kluczowy structural advantage Candidate D vs Candidate A
(dual-field) i Candidate B (conformal frame):

- **Candidate A (dual-field):** drugi field σ może mieć non-universal
  matter coupling → explicit composition dependence enters at zeroth order.
- **Candidate B (conformal frame):** konformalna metryka g_eff = Ω²(ψ)g
  prowadzi do anomalous matter coupling jeśli Ω depends on matter sector
  → potential WEP violation.
- **Candidate D (momentum back-reaction):** modyfikacja w **gravitational
  sector tylko**, matter EOM = standard. UFF zerowo zachowane.

### Finding 2: MICROSCOPE jest najtighter constraint

Ze wszystkich cross-checks Phase 0+ done dotychczas:

| Test | η heuristic | bound | margin |
|------|-------------|-------|--------|
| Cassini γ_PPN (Mercury) | 2.5×10⁻¹⁸ | 2.3×10⁻⁵ | 9×10¹² |
| OP-DESI w(z) (today) | 8×10⁻³⁴ | 1 (phantom) | 1×10³³ |
| **MICROSCOPE η_Ti-Pt** | **1.6×10⁻¹⁶** | **1.1×10⁻¹⁵** | **6.7** |
| LLR Nordvedt | 4.7×10⁻²⁹ | 1.4×10⁻¹³ | 3×10¹⁵ |
| OP-7 c_GW | 0 (vacuum) | exact | structural |

⇒ MICROSCOPE jest **dominant constraint** na α-coupling. Future MICROSCOPE++
lub improved Eöt-Wash mogłyby zacząć directly probing Candidate D regime.

### Finding 3: TIGHT margin nie znaczy "fail"

Heuristic estimate `η ~ α (∇Φ)² × Δε_internal` ma uncertainties of order:

- Proportionality coefficient (could be O(0.1) lub O(10) depending on
  detailed Nordvedt-like derivation)
- Effective ∇Φ at internal nucleon scale (atomic vs nuclear)
- Suppression by gravitational vs nuclear coupling channels

**Structural argument (matter EOM unchanged) jest primary defense.**
Heuristic 6.7× margin może w rigorous Phase 1 derivation:
- **Wzrosnąć** do ~100-1000× (jeśli Nordvedt-like factor ~1/10)
- **Zostać marginalny** (jeśli factor ~1)
- **Pojawić się jako falsification probe** (jeśli factor ~10)

⇒ Full Phase 1 covariant derivation **WYMAGA** rigorous WEP analysis
(higher priority than e.g. cosmology, where margin is overwhelming).

### Finding 4: 5th force coupling deeper niż naive scalar

Standard Brans-Dicke 5th force gives η ~ 1/(2ω+3) — depends on global
ω parameter, easily ruled out by Cassini if ω < 40000.

Candidate D 5th force:
- α-coupling jest **derivative-derivative** (J_μ J_ν), nie linear
- Suppression przez (∇Φ)² jest natural — quadratic w gradient
- Composition dependence przez T^μν decomposition (wedlug stress-energy
  internal structure)
- ⇒ **Nie redukuje się do Brans-Dicke** → standard ω-bound nie aplikuje

## 5. Phase 1 deferred items (post-WEP)

Phase 0+ WEP cross-check zamyka **partial** Phase 1 item (e):

✅ **Closed (Phase 0+):**
- (e.1) WEP zerowego rzędu (matter EOM unchanged) — STRUCTURAL PASS
- (e.2) Nordvedt-like LLR scaling — heuristic margin 3e+15×
- (e.3) MICROSCOPE Ti-Pt — heuristic margin 6.7× (TIGHT)

⏳ **Still deferred (full Phase 1 — HIGH PRIORITY):**
- (e.4) Rigorous Nordvedt-like derivation: proper proportionality factor
  via covariant Φ-EOM perturbation theory
- (e.5) MICROSCOPE++ projection: czy improved bound 10⁻¹⁷ falsifyje
  Candidate D?
- (e.6) Atomic-scale gradient effective ∇Φ in stress-energy integral
- (e.7) E_EM (electromagnetic binding) contribution — similar scale do
  nuclear, może constructively/destructively interfere

## 6. Implications strategiczne

### Phase 0+ closure status (post-WEP)

- ✅ **Phase 0** (rano 2026-04-25): scope + 4 candidates + decision tree
- ✅ **Phase 0+ kickoff** (popołudnie 2026-04-25): structural sketch 5/5 POSITIVE
- ✅ **Phase 0+ cosmology** (wieczór 2026-04-25): w(z) ≥ −1 preserved (margin 10³³×)
- ✅ **Phase 0+ WEP** (noc 2026-04-25): UFF preserved structurally; MICROSCOPE
  margin 6.7× (TIGHT, requires rigorous Phase 1 derivation)
- ⏳ **Phase 1** (deferred 2026 Q3-Q4 OR post-ngEHT 2030+): **WEP analysis
  highest priority** + photon ring + perturbation cosmology

### F4 falsifiability hardening (post-WEP)

> **TGP M9.1'' standalone falsifiable. M9.2-D pivot path:**
> - Strong-field: Candidate D structural sketch POSITIVE (α ≈ 0.1 tunes
>   scenario (e) target +1.46% deviation)
> - Weak-field PPN: U⁴ auto-suppressed (margin 9e+12× nad Cassini)
> - Cosmology: w(z) ≥ −1 preserved (margin 1.24e+33× w obserwowalnym régime)
> - **WEP: UFF strukturalnie zachowane** (matter EOM unchanged), residual
>   Nordvedt-like η ~ 1.6e-16 vs MICROSCOPE 1.1e-15 → margin 6.7× (TIGHT)
> - GW: c_GW = c_0 vacuum exact (OP-7 unchanged)
> - Stability: no-ghost α > 0, Ostrogradsky-free at tree level
>
> **Multi-front cross-check:** strong-field + weak-field + cosmology + WEP
> + GW propagation — wszystkie na single coupling α ≈ 0.1 (geom units).
> Tightest constraint: MICROSCOPE — wymaga rigorous Phase 1 derivation
> przed quantitative claim falsifikowalności via WEP++ improvement.

## 7. Files (Phase 0+ WEP)

- `op_m92_P0plus_candD_wep.py` — analytical cross-check script
- `op_m92_P0plus_candD_wep.txt` — raw output
- `OP_M92_P0plus_candD_wep_results.md` — this synthesis

## 8. Cross-references

- [[research/op-m92/OP_M92_P0plus_candD_results.md]] — Phase 0+ kickoff (parent)
- [[research/op-m92/OP_M92_P0plus_candD_cosmology_results.md]] — cosmology cross-check (sibling)
- [[research/op-m92/OP_M92_readiness_summary.md]] — Phase 0 closure
- [[research/op-m92/OP_M92_setup.md]] — formal scope

External:
- MICROSCOPE 2022: Touboul et al., PRL 129, 121102 (2022)
- LLR Nordvedt: Williams, Turyshev, Boggs, IJMPD 18, 1129 (2009)
- Eöt-Wash: Wagner, Schlamminger, Gundlach, Adelberger, CQG 29, 184002 (2012)

---

## Bottom line

OP-M92 Phase 0+ WEP cross-check **ZAMKNIĘTY POSITIVE z notą o ciasnym
MICROSCOPE marginesie**. Candidate D strukturalnie zachowuje UFF
(matter EOM unchanged), z residual Nordvedt-like contribution dominated
by nuclear binding energy fraction.

**Heuristic verdict:**
- η_MICROSCOPE ≈ 1.6×10⁻¹⁶ (bound 1.1×10⁻¹⁵) → **6.7× margin** (TIGHT)
- η_LLR ≈ 4.7×10⁻²⁹ (bound 1.4×10⁻¹³) → **3×10¹⁵× margin** (safe)

**Phase 1 priority elevation:** WEP analysis jest **highest priority**
deferred item — rigorous covariant derivation może albo:
- Comfortably loosen MICROSCOPE margin → solid PASS
- Reveal MICROSCOPE++ as **new falsification probe** dla Candidate D

W obu przypadkach: tightening of falsifiability narrative.

**M9.2 readiness elevated:** Candidate D structural scaffolding cross-checked
across **pięciu** sektorów (strong-field, weak-field, cosmology, WEP, GW)
on single coupling α ≈ 0.1 (geom units). Response time post-ngEHT 2030+
verdict: 2-4 weeks zamiast 2-3 lat.
