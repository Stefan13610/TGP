---
title: "M11 — Branch strategy (soliton-based vs FRG-based quantization)"
date: 2026-04-26
cycle: M11
sub-cycle: M11.0+ (strategic refinement, pre M11.S/I/G launch)
status: STRATEGY-DRAFT
predecessor: "[[M11_0_drift_audit.md]]"
related:
  - "[[M11_program.md]]"
  - "[[../op-cosmology-closure/M10_R_results.md]]"
  - "[[../op-newton-momentum/M9_3_results.md]]"
  - "[[../continuum_limit/cg_results.txt]]"
tags:
  - TGP
  - M11
  - quantization
  - soliton
  - FRG
  - strategy
---

# M11 — Branch strategy: soliton-based ↔ FRG-based quantization

> **Cel dokumentu:** sformalizować dwutorową strategię kwantyzacji TGP (Branch I = soliton; Branch II = FRG/perturbatywny), wyspecyfikować math relację między branches, zdefiniować warunki konsystencji w M11.R przed jakimkolwiek liczeniem.
> **Powód:** standardowa kwantyzacja perturbatywna wokół Φ=0 lub Φ_0 jest **technicznie podejrzana** dla TGP z trzech strukturalnych powodów (single-Φ axiom, K(φ)=K_geo·φ⁴ degeneracja, V''(Φ_0)=-β tachyon w czasie). Solitonowa kwantyzacja jest fizycznie naturalniejsza (matter = soliton, single-Φ axiom respected); FRG jest universal cross-check (3D Ising class). Trzeba je obie + M11.R weryfikuje konsystencję.

---

## 1. Motywacja: dlaczego dwie gałęzie

### 1.1 Trzy obstacles dla czystej QFT

**(O1) K(φ) = K_geo·φ⁴ degeneruje przy Φ→0**

Kanoniczna kwantyzacja wymaga normalizowalnego kinetycznego propagatora: ½(∂φ)². W TGP:
```
K(0) = 0   ⇒   propagator ill-defined wokół Φ=0
K(Φ_0) = K_geo · 1 ≠ 0   ⇒   OK wokół Φ_0
```
Perturbacja wokół Φ_0 wymaga rescalingu φ → ξ = ∫ √K(φ) dφ żeby ξ było canonical. To wprowadza nieliniową map między fizycznym Φ a quantum field ξ — interpretacyjnie nietrywialne.

**(O2) V''(Φ_0=1) = -β < 0 (slow-roll MAXIMUM w czasie kosmologicznym)**

Z sek08a (β=γ vacuum):
```
V(φ) = (β/3)φ³ - (γ/4)φ⁴ = (β/3)φ³(1 - 3φ/4)|_{β=γ}
V'(φ) = β·φ²(1 - φ),       V'(1) = 0  ✓
V''(φ) = β·φ(2 - 3φ),      V''(1) = -β < 0
```
Tachyonic w czasie ⇒ slow-roll cosmologiczny (M9, M10 confirmed: kosmologia w czasie wyjaśniona). Ale dla **kwantowania w przestrzeni** (M9.3.1 mass term Yukawa) M_eff² = +β > 0 — przeciwny znak. Cel: kwantyzacja musi się zgodzić z **oboma** signami w odpowiednich kanałach.

**(O3) Source coupling (q/Φ_0)·φ·ρ — solitony naturalniejsze niż quanta**

Z sek08a action:
```
S_int = -∫d⁴x √(-g_eff) (q/Φ_0) · φ · ρ
```
ρ to klasyczna density (matter) — lokalnie skupiona wokół point sources. EOM:
```
-∂_μ[K g^μν ∂_νφ] + V'(φ) + (q/Φ_0)ρ = 0
```
Dla ρ = M·δ³(r) (point source) rozwiązaniem jest **klasyczny soliton** Φ_sol(r). Standardowa QFT trakuje ρ jako external, ale TGP ma single-Φ axiom — matter (ρ) i pole (Φ) są na tym samym ontologicznym poziomie. Naturalniejsze: ρ = ρ[Φ_sol], czyli sources są klasycznymi konfiguracjami pola Φ samego.

**Wniosek:** TGP może być teorią tipo Skyrme / sine-Gordon / 't Hooft-Polyakov monopole, gdzie cząstki = solitony, NIE Fock states.

### 1.2 Dlaczego nie wystarczy sam soliton (Branch I)

- Universalność klasy Ising 3D (CG-2 η=0.044) jest własnością FRG fixed-point, NIE soliton spectrum.
- RG running γ(k) (M10.5.4 prediction γ(k_LSS)/γ(k_CMB)=1.244) musi być derived z β-functions → potrzeba Branch II.
- KNOWN_ISSUES C.3 (γ=M_Pl²·g̃) wymaga RG fixed-point analysis.

### 1.3 Dlaczego nie wystarczy sam FRG (Branch II)

- Nie respektuje single-Φ + source axiom (matter wciąż external).
- (O1)+(O2) niewyrugowane perturbatywnie wokół Φ=0/Φ_0.
- Most do M9 klasycznej grawitacji nie jest naturalny (Newton const wynika pośrednio z RG).

⇒ **Obie potrzebne, każda na innym terenie.**

---

## 2. Branch I — Soliton-based (bottom-up)

### 2.1 M11.S — single soliton

**Setup matematyczny:**

Klasyczne EOM w background M9.1'' (hyperbolic):
```
-∂_μ[K(Φ) g_eff^μν ∂_νΦ] + V'(Φ) + (q/Φ_0)ρ = 0
```
Dla statycznej, sferycznie symetrycznej konfiguracji wokół point source M at origin:
```
ρ(r) = M·δ³(r)
Φ_sol = Φ_sol(r),     g_eff = M9.1''  hyperbolic
```
Asymptotyka:
```
Φ_sol(r → ∞) → Φ_0 = 1   (próżnia kosmologiczna)
Φ_sol(r → 0) → Φ_core(M)   (regularna; M9.3.1 ψ→0 limit, no naked singularity)
```

**Linearization (kwantowe fluktuacje):**
```
Φ(x,t) = Φ_sol(r) + δΦ(x,t)
```
Action expanded do O(δΦ²):
```
S_quad[δΦ] = ½ ∫d⁴x √(-g_eff) · K(Φ_sol) · g_eff^μν ∂_μδΦ ∂_νδΦ
            - ½ ∫d⁴x √(-g_eff) · M_sol²(r) · δΦ²
gdzie:
M_sol²(r) = V''(Φ_sol(r)) + K-corrections (od derivatives K w Φ_sol)
```
Operator fluktuacyjny:
```
D̂_sol = -∂_μ[K(Φ_sol) g_eff^μν ∂_ν · ] + M_sol²(r)
```

**Mode expansion + spectrum:**
```
δΦ(x,t) = Σ_n [a_n u_n(x) e^{-iω_n t} + h.c.]
D̂_sol u_n = ω_n² K(Φ_sol) u_n   (eigenvalue problem)
```
**Zero modes** (3 expected, translacyjne):
```
u_0^(i) = ∂_i Φ_sol(r),  i=1,2,3   (Goldstone'y od broken translational sym)
```
**Continuous spectrum** (radiation): ω_n² → β + k² jako asymptotyczna postać (M9.3.1 Yukawa).

**Collective coordinate quantization:**

Standard Christ-Lee / Tomboulis:
```
Φ(x,t) = Φ_sol(x - R(t)) + δΦ_⊥(x - R(t), t)
```
- R(t) → quantum DOF, conjugate P = M·dR/dt
- δΦ_⊥ orthogonal do zero modes (radiation only)
- Hamiltonian: H = P²/(2M_phys) + H_rad
- M_phys = M_class + ½ Σ_n ω_n - counterterms (1-loop renorm)

**Pre-test hypotheses:**
- **H1:** Φ_sol(r) istnieje jako nontrivial classical solution (numerical shooting from r=0 to r→∞)
- **H2:** Φ_sol(r) regularna w r=0 (no naked singularity); Φ_sol(r) → Φ_0 dla r → ∞
- **H3:** {ω_n²} ≥ 0 dla wszystkich n except 3 zero modes (kwantowa stabilność wokół Φ_sol)
- **H4:** δM_phys/M_class finite po regularization (UV cutoff or zeta-function)

**Deliverables:**
- `m11_S_soliton.py`: numerical solver Φ_sol(r), eigenmode spectrum
- `M11_S_results.md`: closure-grade audit ≥6/6 PASS

### 2.2 M11.I — multi-soliton interference

**Setup:**

Two-soliton ansatz (asymptotic):
```
Φ_2sol(r; r_1, r_2) ≈ Φ_sol(r-r_1) + Φ_sol(r-r_2) - Φ_0   (background subtraction)
```
Lepszy ansatz: full nonlinear EOM z dwoma sources, solved numerically.

**Interakcja:**
```
V_int(r_12) = S[Φ_2sol] - S[Φ_sol(r-r_1)] - S[Φ_sol(r-r_2)]
```
Asymptotyka dla r_12 ≫ λ_C = 1/√β (Yukawa range ~ Hubble radius dla β~H_0²):
```
V_int(r_12) ≃ -G_TGP · M² · exp(-r_12 √β) / r_12   (Yukawa)
```
gdzie G_TGP wynika z structure sek08a (M9 already gives leading Newton const).

**Quantum correction:**

Funkcjonalny determinant od fluktuacji wokół Φ_2sol vs around Φ_sol osobno:
```
δV_quant(r_12) = ½ Tr ln D̂_2sol - ½ Tr ln D̂_sol(r_1) - ½ Tr ln D̂_sol(r_2)
```
For r_12 ≫ λ_C: Casimir-like contribution ~ exp(-2 r_12 √β)/r_12² (sub-leading).

**Pre-test hypotheses:**
- **H5:** V_int(r_12 → ∞) → -G_TGP·M²·exp(-r_12/λ_C)/r_12 (Yukawa, M9 reproduced)
- **H6:** V_int(r_12 → 0) finite (regularny merge)
- **H7:** N-soliton gas → mean-field = M9 klasyczna grawitacja (consistency)

**Deliverables:**
- `m11_I_interference.py`: 2-soliton solver + V_int(r_12) numerically
- `M11_I_results.md`: closure-grade ≥6/6 PASS

### 2.3 M11.G — global field z source extraction

**Setup:**

Decompozycja:
```
Φ(x,t) = Φ_cl[{r_i(t)}](x) + δΦ_rad(x,t)
```
gdzie:
- Φ_cl klasycznie spełnia EOM dla source ρ = Σ_i M_i·δ³(x-r_i(t))
- δΦ_rad = quantum radiation, orthogonal do zero modes każdego solitona
- {r_i(t)} = klasyczne (c-number) trajektorie; ich kwantyzacja jest collective coord problem (M11.S)

**Effective action:**
```
S_eff[r_i, δΦ_rad] = S_cl[Φ_cl[{r_i}]] + ½∫δΦ_rad D̂_global δΦ_rad + S_int(higher order)
```
Integrate out δΦ_rad at 1-loop:
```
W[r_i] = -i·ln Z[r_i] = S_cl[r_i] - (i/2) Tr ln D̂_global[r_i] + O(ℏ²)
V_eff({r_i}) ← real part of W
```

**Match z M9 + 1-loop quantum corrections:**

Leading order (mean-field, classical):
```
V_eff^(0)({r_i}) = -Σ_{i<j} G_TGP M_i M_j exp(-|r_i-r_j|/λ_C) / |r_i-r_j|
                                                                  ↑
                                                          M9 Newton + Yukawa
```

1-loop correction (quantum vacuum polarization):
```
V_eff^(1)({r_i}) ~ G_TGP² Σ_{i<j} M_i M_j · F(|r_i-r_j|, β, η) / |r_i-r_j|³
```
gdzie F zawiera anomalous dim η = 0.044 (jeżeli Branch I ↔ Branch II self-consistent).

**Pre-test hypotheses:**
- **H8:** Mean-field over δΦ_rad → exactly M9 Φ_0(r) profile (consistency check)
- **H9:** 1-loop V_eff(r) corrections sub-leading przy classical Newtonian regime (post-Newtonian quantum)
- **H10:** η wynikające z anomalous dim δΦ_rad operator wokół background = 0.044 (= CG-2; Branch I ↔ Branch II convergence)

**Deliverables:**
- `m11_G_extraction.py`: 1-loop V_eff({r_i}) computation + match z M9
- `M11_G_results.md`: closure-grade ≥6/6 PASS, **kluczowe** dla branch consistency

---

## 3. Branch II — FRG-based (top-down, oryginalny plan)

Bez zmian względem M11_program.md sekcja 3.

### 3.1 M11.1 — 1-loop V_eff (m2b audit)

Background field method:
```
V_eff(Φ̄) = V_tree(Φ̄) + ½ Tr ln(D̂² + M²(Φ̄))
M²(Φ̄) = V''(Φ̄) = β·Φ̄·(2 - 3Φ̄)
```
Audit `m2b_loop.py` against canonical sek08a.

### 3.2 M11.2 — β-functions FRG (Wetterich)

Wetterich equation:
```
∂_t Γ_k = ½ STr [(Γ_k^(2) + R_k)^{-1} ∂_t R_k]
```
LPA': K(φ), V(φ), η = -∂_t ln Z_k running. Reproduce CG-2 η=0.044.

### 3.3 M11.3 — γ(k) running

Solve flow eq. Predict γ(k_LSS)/γ(k_CMB) = 1.244 (M10.5.4).

### 3.4 M11.4 — RG-driven structural derivation

C.3, B.3, B.5, B.2 closure z fixed-point analysis.

---

## 4. Branch consistency conditions (M11.R)

**To są warunki które MUSZĄ być spełnione żeby M11 cycle się zamknął:**

### 4.1 Anomalous dimension match

```
η_BI := anomalous dim δΦ_rad operator wokół Φ_sol/Φ_2sol/Φ_cl background  (Branch I, M11.G H10)
η_BII := LPA' running of Z_k          (Branch II, M11.2)
η_CG2 := 0.044                         (CG-2 closed reference)

Wymóg M11.R:
    |η_BI - η_BII| < 0.01
    |η_BI - η_CG2| < 0.01
    |η_BII - η_CG2| < 0.01
```

Jeśli NIE → MAJOR INCONSISTENCY, M11 cycle wraca do drawing board.

### 4.2 Effective Newton constant agreement

```
G_TGP^{Branch I}  := from V_int(r_12) Yukawa coefficient (M11.I, H5)
G_TGP^{Branch II} := from γ at IR fixed point (M11.4: γ = M_Pl² · g̃)

Wymóg M11.R:
    |G_TGP^{BI}/G_TGP^{BII} - 1| < 0.01
```
Plus oba muszą reproducować Newton G_N = (8πG)/(c⁴) standardową kalibrację (T-Λ + T-α).

### 4.3 Mass scale self-consistency

```
λ_C^{Branch I}  := soliton Yukawa range = 1/√β (z M11.I H5)
λ_C^{Branch II} := = 1/√β z RG-flowed β at IR (M11.2 z initial β/H_0² = 1)

Wymóg:
    λ_C^{BI} = λ_C^{BII}   (analytically, both should equal 1/H_0 = Hubble radius w β/H_0²=1)
```

### 4.4 Universality class agreement

```
Branch I: soliton phase transition uniwersalność (kondensacja solitonów); class TBD przez M11.S/I
Branch II: 3D Ising (CG-2 verified)

Wymóg M11.R:
    Both → 3D Ising, η = 0.044, ν ≈ 0.649
```

Jeśli Branch I daje INNĄ klasę → problem fundamentalny.

### 4.5 KNOWN_ISSUES closure agreement

```
C.3 (γ = M_Pl² · g̃):     Branch II M11.4 derivation MUST agree with Branch I extracted G_N
B.3 (α₀ ≈ 4):              Branch II M11.4 RG running MUST agree with Branch I scattering off solitons (jeśli computable)
B.5 (g̃ ≈ 0.98 vs 1.0):    Branch II RG-corrected = Branch I extracted within tolerance
B.2 (n = 2):               Branch I quantum stability (H3) + Branch II RG smoothness — both n=2 stable
```

### 4.6 Most do M9 (klasyczna grawitacja)

```
Branch I (M11.G mean-field, leading order) MUST exactly reproduce M9 Φ_0(r) profile
Branch II (M11.4 g̃ at IR + γ) MUST give Newton G_N matching M9 calibration
```

---

## 5. Risks & mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| **R1:** Φ_sol(r) nontrivial nie istnieje (e.g. tylko Φ=Φ_0 trivial sol) | LOW | HIGH | M11.S H1 numerical shooting; jeśli FAIL → revert to perturbative-only |
| **R2:** Spectrum {ω_n²} ma negative modes (instability) | MED | HIGH | M11.S H3; jeśli FAIL → soliton metastable not stable, requires re-thinking |
| **R3:** η_BI ≠ η_BII | MED | HIGH | M11.R 4.1 explicit check; jeśli FAIL → fundamental incompatibility |
| **R4:** G_TGP mismatch | LOW | HIGH | M11.R 4.2; possibly różne schemes — match w common scheme |
| **R5:** Soliton CoM quantization ordering ambiguity | LOW | LOW | Use established Christ-Lee 1975 / Gervais-Sakita 1975 prescription |
| **R6:** UV divergences w M11.S 1-loop nie kasują się standardowymi counterterms | MED | MED | Zeta regularization (Hawking-style), or dimensional reg adapted to non-canonical K |
| **R7:** Soliton Compton ~ Hubble radius (β~H_0²) → "cząstka rozmiaru obserwowalnego Wszechświata" — interpretacyjny problem | HIGH | LOW | To **feature** TGP, nie bug — soliton = cosmological mass distribution; dla local objects effective mass scale jest wyższy z renormalization |

---

## 6. Sub-cycle execution order

```
Phase 1 (parallel possible):
  M11.S (single soliton)         ⊕   M11.1 (m2b 1-loop audit)
        ↓                                ↓
Phase 2:
  M11.I (multi-soliton)              M11.2 (FRG β-functions)
        ↓                                ↓
Phase 3:
  M11.G (global extraction)          M11.3 (γ(k) running)
                                        ↓
                                     M11.4 (RG-driven C.3/B.3/B.5)

Phase 4 (final):
  M11.R — branch consistency synthesis
        verifies §4.1–4.6 conditions
        aggregates 7 sub-cycle PASS counts (target: ≥42/42)
```

**Krytyczna ścieżka:** M11.S → M11.G → M11.R (Branch I core).
**Niezależny tor:** M11.1 → M11.2 → M11.4 (Branch II).

Iteracja: M11.S najpierw, bo bez Φ_sol(r) classical solution nie ma jak liczyć M11.I/G.

---

## 7. Pre-test hypothesis matrix (consolidated)

| ID | Sub-cycle | Hipoteza | Verification |
|---|---|---|---|
| H1 | M11.S | Φ_sol(r) nontrivial istnieje | Numerical shooting |
| H2 | M11.S | Φ_sol(0) regularne, Φ_sol(∞) → Φ_0 | Asymptotic matching |
| H3 | M11.S | Spectrum stabilny (3 zero modes only) | Eigenvalue solver |
| H4 | M11.S | δM 1-loop finite | Zeta or dim reg |
| H5 | M11.I | V_int(∞) → Yukawa, M9 reproduced | Asymptotic match |
| H6 | M11.I | V_int(0) finite | Numerical small-r |
| H7 | M11.I | N-soliton mean-field = M9 | Statistical mechanics |
| H8 | M11.G | <Φ> = M9 Φ_0(r) profile | M9 cross-check |
| H9 | M11.G | 1-loop sub-leading | Power-counting |
| H10 | M11.G | η_BI = 0.044 | M11.R §4.1 critical |
| H11 | M11.1 | m2b ↔ canonical sek08a equivalent | Audit |
| H12 | M11.2 | LPA' reproduces η = 0.044 | CG-2 cross-check |
| H13 | M11.3 | γ(k_LSS)/γ(k_CMB) = 1.244 | M10.5.4 prediction |
| H14 | M11.4 | g̃, α₀ from RG fixed-point | KNOWN_ISSUES C.3/B.3/B.5 |

---

## 8. Deliverables map

```
M11.S  → m11_S_soliton.py        + M11_S_results.md
M11.I  → m11_I_interference.py   + M11_I_results.md
M11.G  → m11_G_extraction.py     + M11_G_results.md
M11.1  → m11_1_audit.py          + M11_1_audit_results.md
M11.2  → m11_2_betafn.py         + M11_2_results.md
M11.3  → m11_3_gamma_k.py        + M11_3_results.md
M11.4  → m11_4_structural.py     + M11_4_results.md
M11.R  → m11_R_consistency.py    + M11_R_results.md  (final; aggregates 7+R)
```

Plus update:
- `M11_program.md` (sub-cycle table → 8 rows: 0/S/I/G/1/2/3/4/R)
- `[[../closure_2026-04-26/KNOWN_ISSUES.md]]` (C.3/B.3/B.5/B.2 closure, A.13 entry)
- Cross-link w `[[../op-cosmology-closure/M10_program.md]]` jeśli M10.5.4 prediction (γ(k_LSS)/γ(k_CMB)=1.244) zweryfikowana w M11.3.

---

## 9. References

- **TGP foundations:**
  - sek08a action: [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  - M9.1'' hyperbolic metric: [[../op-newton-momentum/M9_3_results.md]]
  - M9.3.1 stable Yukawa M_eff²=+β: ditto
  - T-Λ closure (β~H_0²): [[../closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]]
  - CG-2 LPA' η=0.044: [[../continuum_limit/cg_results.txt]]
  - M10 cycle CLOSED: [[../op-cosmology-closure/M10_R_results.md]]

- **Soliton quantization (literature):**
  - Christ & Lee 1975 (Phys. Rev. D 12, 1606) — CoM separation
  - Gervais & Sakita 1975 (Phys. Rev. D 11, 2943) — collective coords
  - Tomboulis 1976 (Phys. Rev. D 12, 1678) — quantization w soliton sectors
  - Rajaraman "Solitons and Instantons" 1982 — book ref
  - Coleman "Aspects of Symmetry" 1985 — semiclassical methods

- **FRG / Wetterich:**
  - Wetterich 1993 (Phys. Lett. B 301, 90) — exact RG eq
  - Berges, Tetradis, Wetterich 2002 (Phys. Rep. 363, 223) — review
  - LPA': Tetradis & Wetterich 1994 (Nucl. Phys. B 422, 541)

---

## 10. Verdict M11.0+

**Strategia OK do execution.** Branch I (M11.S → M11.I → M11.G) jest fizycznie naturalniejsza dla TGP single-Φ + source coupling axiom. Branch II (M11.1-4) zachowana jako universal cross-check (Ising class + RG running).

Konsystencja branches w M11.R (§4) jest kryterium zamknięcia M11 cyklu. Jeśli §4.1/4.2 FAIL → reset.

**Następny krok:** mini-PoC dla M11.S — sprawdzić czy Φ_sol(r) nontrivial existuje numerycznie (H1) ZANIM commitujemy się formalnie do całej restruktury M11_program.md.

---

*M11.0+ strategy draft 2026-04-26. Pre-execution scaffolding before M11.S launch.*
