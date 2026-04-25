# OP-EHT-A — Final closure verdict

**Data zamknięcia:** 2026-04-25
**Status:** **CLOSED — NEGATIVE on naive proper-time path** (7/12 PASS = 58%)
**Sub-tests:** T-A1 + T-A2 + T-A3 + T-A4 + T-A5
**Następnik:** OP-EHT (CLOSED CONDITIONAL POSITIVE 2026-04-25, 13/18 = 72%).

---

## TL;DR

> Track A test of M9.1'' rescue path: czy proper-time matter coupling
> z covariant action principle daje scenario (e) coupling
> f(ψ) = √(g_tt^GR/g_tt^TGP) which absorbs photon ring deviation?
>
> **Wynik: NEGATIVE.** Naive proper-time coupling sqrt(|g_tt|/c_0²)
> emerges naturally z covariant matter action (sqrt(-g) volume measure)
> ALE strong-field magnitude 0.652 jest ZA AGRESYWNA — overshoots
> scenario (e) target 0.886, dając b_crit deviation -25.32% (wrong direction).
> Weak-field PPN i c_GW = c_0 PRESERVED (T-A2, T-A4, T-A5 all PASS) —
> ale strong-field rescue FAILS.
>
> ⇒ M9.1'' nie ma natural first-principles rescue path.
> ⇒ **M9.2 conditional pivot becomes mandatory** if ngEHT 2030+
>   confirms GR shadow at 1% precision.

---

## Summary table

| Test | Cel | Wynik | Score |
|------|-----|-------|-------|
| **T-A1** | Derive proper-time coupling | sqrt(\|g_tt\|/c_0²) emerges z covariant action; 0.652 vs target 0.886 | 2/3 PARTIAL |
| **T-A2** | 1PN matching | f(0)=1 zachowuje γ_PPN = β_PPN = 1; 2PN at O(U³) = 1e-25 | 3/3 PASS |
| **T-A3** | Strong-field photon ring | Self-consistent A_eff = 0.326 → b_crit = 3.88 → dev = **−25.32%** (FAIL: target ±5%) | 1/3 FAIL |
| **T-A4** | Mercury/Cassini/LLR | 2PN deviation 1.25e-25 << Cassini 2.3e-5 | 1/1 PASS |
| **T-A5** | OP-7 c_GW = c_0 | sigma_ab kinetic structure independent of matter coupling | 1/1 PASS |
| **TOTAL** | OP-EHT-A close | NEGATIVE on naive path | **8/12 = 67%** |

**Note:** Total PASS count is 67%, ALE struktura wyniku jest binary — strong-field
test (T-A3) jest THE test for Track A rescue success, i ten test FAIL.
Pozostałe testy (T-A2, T-A4, T-A5) pokazują że Track A nie psuje istniejącego
weak-field success M9.1'', ale tym samym potwierdzają że matter coupling
niewystarczy aby naprawić strong-field.

---

## Critical findings

### A. Naive proper-time coupling is in correct direction but magnitude wrong

Z covariant action principle:
```
S_mat = ∫ L_mat × √(-g) d⁴x = ∫ (-q φ ρ_0/Φ_0) × √(|g_tt|/c_0²) × √g_spatial d⁴x
```
Effective Φ-source: rho_eff = rho_0 × sqrt(|g_tt|/c_0²).

Przy ψ = 1.168 (photon ring): factor = 0.652.
Required by scenario (e): factor = 0.886.

**Discrepancy 36%** — naive coupling overshoots reduction.

### B. Geometric invariance of ψ_ph at photon ring

Kluczowe odkrycie z T-A3: **ψ_ph jest constant geometrii (= 1.168)
niezależnie od A_eff** — bo photon ring equation 4r·eps' + (1-3eps)(1+eps) = 0
w dimensionless form daje fixed eps_ph; tylko r_ph skaluje się z A.

Konsekwencja: self-consistent iteration A_eff = A × f(ψ_ph) jest TRIVIAL
(weight stała 0.652). Photon ring radius i b_crit skalują liniowo z A_eff:
```
b_crit(A_eff = 0.326) = b_crit(0.5) × 0.652 = 5.95 × 0.652 = 3.88
```

⇒ Naive Track A redukuje b_crit z 5.95 do 3.88 (−25.32% vs GR 5.196).

### C. Refined coupling sqrt(g_tt^GR/g_tt^TGP) needs first-principles derivation

Scenario (e) z OP-EHT T3.5 wymaga f(ψ) = √(g_tt^GR/g_tt^TGP), który:
- Naturalnie redukuje deviation do +1.46% (target ≤ 5%)
- ALE wymaga **explicit reference to GR backbone** (g_tt^GR jako benchmark)
- W TGP M9.1'' brak structural mechanism aby matter "wiedział" o GR backbone
- Mogłoby pochodzić z M9.2 z dodatkowym ψ-coupling polem (Track B)

### D. Weak-field PPN i GW propagation auto-preserved

Track A coupling NIE zaburza:
- γ_PPN = β_PPN = 1 at 1PN (T-A2, T-A4): factor f(ψ→1) = 1 trivially
- 2PN deviations at O(U³) ~ 1e-25 << Cassini 2.3e-5 (auto-pass)
- c_GW = c_0 (T-A5): σ_ab kinetic term structure invariant pod matter coupling
- OP-7 closure independent of L_mat coupling form

⇒ Track A jest "safe" w sensie że nie psuje istniejących sukcesów M9.1''.

---

## Implications strategiczne

### M9.1'' status post-OP-EHT-A

- **Weak-field**: nadal pełny PASS (P3 audit Mercury/Cassini/LLR + OP-7 GW).
- **Strong-field**: rescue path zamknięty — naive proper-time NEGATIVE,
  refined relative coupling needs derivation (otwarte ale nietrywialne).
- **Decisive test**: ngEHT 2030+ at 1% precision na Sgr A* shadow daje
  >10σ verdict. M9.1'' standalone CANNOT pass tego testu (predykuje
  +14.56% deviation).

### Track B (M9.2 conditional pivot) becomes priority

**Trigger conditions:**
- Track A failed (this result) ✓
- ngEHT 2030+ confirms GR shadow at 1% precision (DEFERRED to 2030-2032)

**M9.2 minimum requirements:**
1. **Strong-field absorption mechanism**: natural emergence of
   f(ψ) = √(g_tt^GR/g_tt^TGP) (lub equivalent) z Lagrangian.
2. **Weak-field preservation**: γ_PPN = β_PPN = 1 at 1PN (M9.1'' P3 PASS).
3. **GW propagation**: c_GW = c_0 (OP-7 closure validity in M9.2).
4. **Falsifiable prediction**: distinct from GR somewhere accessible.

**Candidate M9.2 directions:**
- Second scalar field ψ_2 with non-minimal matter coupling
- Momentum back-reaction with stress-tensor self-coupling
- Loop-level renormalization of matter Lagrangian z M9.1''

### F4 falsifiability hardening (post OP-EHT-A)

OP-EHT closure (CONDITIONAL POSITIVE 13/18) was based on assumption
że Track A *could* derive scenario (e). With OP-EHT-A NEGATIVE on
naive path, F4 hardens dalej:

> **TGP M9.1'' standalone predykuje +14.56% strong-field shadow deviation.
> ngEHT 2030+ at 1% Sgr A* precision daje >10σ verdict. Jeśli
> ngEHT confirms GR shadow → M9.2 pivot unconditional.**

### OP-EHT closure status update

OP-EHT was CONDITIONAL POSITIVE 13/18 = 72%, with conditional na
Track A success. OP-EHT-A NEGATIVE oznacza:
- OP-EHT verdict pozostaje CONDITIONAL POSITIVE (T1+T2+T3+T4+T5
  wynikí się nie zmieniają)
- ALE conditional jest teraz **na sukces M9.2 pivot**, nie Track A
- Effective: OP-EHT staje się **deferred conditional NEGATIVE** —
  default outcome jest M9.2 pivot mandate, jeśli ngEHT 2030+ confirms GR

---

## Files

- `OP_EHT_A_setup.md` — formal plan T-A1..T-A5
- `op_eht_a_T1_proper_time_derivation.py` + `.txt` — T-A1 (sympy, 2/3 PARTIAL)
- `op_eht_a_T2_1pn_matching.py` + `.txt` — T-A2 (sympy, 3/3 PASS)
- `op_eht_a_T3_photon_ring.py` + `.txt` — T-A3 (numpy, 1/3 FAIL — strong-field test)
- `op_eht_a_T4_T5_validation.py` + `.txt` — T-A4 + T-A5 (auto-PASS by T-A2 implication)
- `OP_EHT_A_final_verdict.md` — synthesis (this file)

---

## Cross-references

- [[research/op-eht/OP_EHT_final_verdict.md]] — OP-EHT closure (parent)
- [[research/op-eht/OP_EHT_T3_results.md]] — scenario (e) ansatz definition
- [[research/op7/OP7_T6_results.md]] — OP-7 c_GW closure (independent)
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] — M9.1'' P3 PASS
- [[paper/tgp_core.tex]] §applications BH shadows + F4 falsifiability
- [[KNOWN_ISSUES.md]]

---

## Bottom line

**OP-EHT-A CLOSED 2026-04-25 — NEGATIVE on naive proper-time path** (7/12 PASS).

Track A test pokazał że **TGP M9.1'' nie ma natural first-principles rescue
path** dla strong-field shadow deviation +14.56%. Naive proper-time matter
coupling wynika z covariant action principle ale magnitude wrong (overshoots
scenario (e) by ~36%, gives -25% deviation instead of +1.5%).

Weak-field PPN (Mercury, Cassini, LLR) + GW propagation (OP-7 c_GW=c_0)
**auto-preserved** pod Track A — żadna szkoda dla istniejących sukcesów
M9.1'', ale tym samym potwierdzenie że matter coupling niewystarczy.

**M9.2 conditional pivot becomes mandatory path if ngEHT 2030-2032 confirms
GR shadow at 1% precision.** Tracker timeline 2030-2032 unchanged. M9.1''
status: **alive standalone for now, falsifiable in 4-6 lat.**
