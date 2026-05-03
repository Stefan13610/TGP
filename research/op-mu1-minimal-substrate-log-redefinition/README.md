---
title: "μ.1 — Minimal substrate redefinition (ψ ≡ log g)"
date: 2026-05-02
cycle: μ.1
status: NO-GO CLOSURE — reparametryzacja PASS (mathematically invariant), ale Σε=2 topology NIE-DERIVED
phase1_score: 3/3 PASS (sympy derivation + numerical consistency 4.6e-11 + μ/e gate 1.2e-13%)
phase2_score: 1/3 (P2.1 PASS, P2.2 FAIL 0/4 candidates, P2.3 NEUTRAL tautological)
overall_verdict: REPARAMETRIZATION SUCCESS, MECHANISM NOT FOUND
parent: "[[../op-lambda1-e2-amplitude-emergence/README.md]]"
predecessor_plan: "[[PLAN.md]]"
related:
  - "[[../op-lambda1-e2-amplitude-emergence/EXTERNAL_AUDIT_2026-05-02.md]]"
  - "[[../op-lambda1-e2-amplitude-emergence/phase2_M6_compound_interference.py]]"
tags:
  - TGP
  - mu1
  - substrate-redefinition
  - log-variable
  - NO-GO
  - reparametrization
  - valuable-negative
tgp_status:
  folder_status: active
  level: L1
  kind: derivation
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'μ.1 — Minimal substrate redefinition (ψ ≡ log g)'"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

# μ.1 — Minimal substrate redefinition (ψ ≡ log g)

> **Status: NO-GO CLOSURE** zgodnie z explicit kryterium z [[PLAN.md]] §7.
>
> **TL;DR:** Reparametryzacja `ψ = log g` JEST mathematycznie identyczna z g-substrate
> (P1: drift 1e-13% w μ/e ratio). Compound formula `g_total = exp(Σε_i)` emerges
> automatycznie z linear superposition w ψ. **ALE:** żaden z 4 testowanych
> topologicznych kandydatów nie daje Σε = 2 z first principles (P2.2 0/4 PASS).
>
> ⇒ μ.1 ≡ **PURE RELABELING** bez fizycznych konsekwencji.
>
> **λ.1 status: zostaje NEGATIVE CLOSURE.** μ.1 nie reopen.

---

## 1. Rezultaty Phase 1 — reparametryzacja działa perfekcyjnie

### P1.1 Sympy derivation ψ-ODE — PASS

Pod g(r) = exp(ψ(r)), R3 ODE:

$$g'' + \frac{\alpha}{g}(g')^2 + \frac{2}{r}g' = (1-g)\, g^{2-2\alpha}$$

przekształca się w:

$$\boxed{\psi'' + (1+\alpha)(\psi')^2 + \frac{2}{r}\psi' = (1 - e^\psi)\, e^{(1-2\alpha)\psi}}$$

Sympy potwierdza: `diff = 0` po `powsimp(force=True)`.

**Specjalność α=1:** RHS = e^(−ψ) − 1 (najsimplest exp form, BEZ r-zależności w prefactorze)
— zgodne z bridge theorem (α=1 IFF jako critical case, possibly integrable).

**Dla α=2 (charged-lepton):** RHS = e^(−3ψ) − e^(−2ψ), vacuum stable, linear restoring
−ψ przy ψ → 0.

### P1.2 Numerical consistency — PASS spectacularly

ODE solved równolegle w g i ψ dla α=2, g₀=1.77472:

| Metric | Wartość |
|--------|---------|
| Max drift \|ψ(r) − log g(r)\| | **4.6 × 10⁻¹¹** |
| RMS drift | 2.2 × 10⁻¹¹ |
| Range r | [0.001, 100.0], 2000 points |

ψ(r) ≡ log g(r) do floating-point precision. **Identyczne soliton solutions.**

### P1.3 Mass formula μ/e gate — PASS spectacularly

Mass formula `m = c·A²·exp(X·ψ₀)` z X = e²/2 = 3.694528:

| Ratio | Predicted | Observed | Drift |
|-------|-----------|----------|-------|
| μ/e | 206.768283 | 206.768283 | **1.24 × 10⁻¹³%** |
| τ/e | 3477.228280 | 3477.228280 | **1.96 × 10⁻¹³%** |

Wszystkie λ.1 sukcesy mass formula **preserved** pod ψ-redefinicją. To jest **invariance test** — zmiana zmiennej nie psuje numeryki (jak oczekiwano dla reparametryzacji).

**Phase 1 GATE: 3/3 PASS.** Reparametryzacja jest exact.

---

## 2. Rezultaty Phase 2 — kluczowy test KEKLOSURE

### P2.1 Multi-soliton ψ-superposition — PASS (analytic)

Pod założeniem (μ.1 hipoteza): ψ_total(x) = Σ_i ψ_i(x − x_i) [linear w ψ]

Konsekwencja:

$$g_{\text{total}}(x) = \exp(\psi_{\text{total}}) = \prod_i g_i(x) \quad\text{(multiplicative)}$$

Linearyzacja dla N małych perturbacji ε_i = ε:

$$g_{\text{total}} \to \exp(N\cdot\epsilon)\cdot\exp(-N\epsilon^2/2) \stackrel{N\epsilon=S\text{ fixed}, \epsilon\to 0}{\longrightarrow} \exp(S)$$

**Compound emergence trywialna w ψ-substrate** — to nie hipoteza, to identity.

### P2.2 Topology Σε = 2 — **FAIL 0/4 (KRYTYCZNY)**

Cztery kandydaci dla derivation Σε = 2 z TGP first principles:

| # | Kandydat | Argument | Score | Verdict |
|---|----------|----------|-------|---------|
| (a) | Electric charge | charged-lepton + antiparticle = \|q\|²·2 | **2/10** | ad hoc |
| (b) | Spinor double-cover | SU(2) → SO(3), 4π identity | **4/10** | wymaga skalar↔spinor mapping |
| (c) | R3 winding n=2 | π_3 homotopy + bundle ℝ³/(±1) | **2/10** | brak natural topology |
| (d) | Compound saturation | 2 channels (amplitude+phase) | **5/10** | plausible ale empirical |

GATE: ≥7/10 = PASS. **Wszystkie FAIL.** Best score: (d) z 5/10.

**Konsekwencja:** μ.1 znajduje strukturalną formę `g_total = exp(Σε)`, ale **wartość**
Σε = 2 dla α=2 charged-lepton **nie jest derived** — tylko empirically postulated.

### P2.3 Compound value match — NEUTRAL (tautological)

Postulując Σε = 2, otrzymujemy g_bg = exp(2) = e² ≈ 7.389. Mass formula slope:

$$X = e^2/2 = 3.694528 \quad\text{(observed in λ.1)}$$
$$X = \exp(\Sigma\epsilon)/2 = \exp(2)/2 = 3.694528 \quad\text{(predicted with Σε=2)}$$

Match exact, ale **cyrkularny**: bez niezależnego derivation Σε=2 (P2.2 FAIL),
"derivation" X = exp(2)/2 jest po prostu definicja.

---

## 3. Verdict ostateczny

### 3.1 Co μ.1 udowodniło POSITIVELY

1. **ψ ≡ log g jest exact reparametryzacją** — matematyczna równoważność do floating-point precision
2. **Compound formula `g_total = exp(Σε)` emerges automatycznie** w ψ-substrate
   z linear superposition (analytic, sympy verified)
3. **Mass formula μ/e i τ/e match preserved** pod redefinicją (1e-13% drift)
4. **α=1 specjalność** w ψ-formulation (czysty exp without r-dep) zgodna z bridge theorem

### 3.2 Co μ.1 udowodniło NEGATIVELY

1. **Σε = 2 nie ma natural source** w TGP-substrate (0/4 candidates pass)
2. **Compound interpretation jest tautologiczna** bez topology argumentu
3. **μ.1 ≡ reparametryzacja**, nie zmiana fizyki

### 3.3 Implikacje dla λ.1

**λ.1 NEGATIVE CLOSURE NIE jest reopened.** μ.1 nie znajduje mechanizmu derivacji
X = e²/2; tylko zmienia zmienną w której opisujemy ten sam empirical fit.

Lambda.1 zostaje:
- 6/6 mechanisms tested NEGATIVE (P2.1, P2.2, P2.3, M.4, M.5, M.6)
- + μ.1 reparametrization NO-GO (substrate redefinition explored, mathematically valid, fizycznie pusty)

**Total uczciwie wykluczonych ścieżek: 7.** μ.1 dodaje wartościowy negative,
nie reverses λ.1.

### 3.4 Co BY potrzebowało żeby mechanism znaleźć

Aby μ.1 mogło dać POSITIVE closure, potrzeba **niezależnego argumentu** dla Σε = 2:
- Albo eksplicit definition amplitude/phase channels w substrate (kandydat d, najlepszy)
- Albo nowa topology TGP-substratu (kandydat c, wymaga rozszerzenia definicji)
- Albo wyprowadzenie Σε = 2 z gauge/QFT structure (kandydat a, wymaga nowych założeń)

Żaden z tych nie jest "minimal change" — wszystkie wymagają **dodania struktury** do
substratu, czego user explicitly chciał uniknąć ("nie chcę dodawać reguł na pałę").

**μ.1 dojrzewa jako conscious negative:** próbowaliśmy, znaleźliśmy że minimalna
zmiana (reparametryzacja) jest mathematically pewna ale fizycznie pusta. Większa
zmiana wymagałaby "dodania reguł" — odłożone.

---

## 4. Pliki

- `PLAN.md` — pre-implementation plan (10 sekcji)
- `phase1_psi_ode_derivation.py` + `.txt` — P1.1 sympy proof
- `phase1_psi_soliton_and_mass_formula.py` + `.txt` — P1.2 numerical + P1.3 μ/e gate
- `phase2_compound_emergence_and_topology.py` + `.txt` — P2.1 multi-soliton + P2.2 topology + P2.3 verify
- `README.md` (this file) — synthesis + verdict

## 5. Co dalej (TGP-program)

Per [[PLAN.md]] §8 ("Relacja do λ.1 i innych cykli"):

**μ.1 closes cleanly bez kosztów dla λ.1.** Następne naturalne cele:

1. **G1** (Φ_eff anchor inconsistency resolution) — pierwotny next-target,
   nie blokowany przez μ.1 outcome
2. **Inne sektory portfolio TGP** (cosmology, BH, mixing) — gdzie compound
   formulas mogą mieć inne zastosowanie z natural Σε derivable

**Lessons learned z μ.1 dla program:**

- Reparametryzacja ≠ zmiana fizyki — różnica między relabeling a substrate
  modification jest ostra
- Compound emergence wymaga mechanizmu wymuszającego konkretną Σε wartość;
  sama struktura ψ-superposition jest niewystarczająca
- "Minimalność" zmiany substratu jest zwodnicza — najmniejsza zmiana często
  oznacza brak fizycznego efektu

---

## 6. Status syntactic

**μ.1 ZAMKNIĘTE: NO-GO CLOSURE z reparametryzacją PASS i mechanism FAIL.**

Per user "Go" decision (2026-05-02 wieczór) — uruchomiono cykl, P1+P2 wykonane,
NO-GO trigger z PLAN.md §7 zadziałał (P2.2 fail wszystkich 4 kandydatów dla Σε=2).

**μ.1 produkuje wartościowe negatives:**
- Reparametryzacja jest mathematycznie czysta (P1 PASS gate)
- Compound emergence formy `exp(Σε)` jest structurally achievable
- ALE wartość Σε = 2 dla α=2 wymaga external mechanism, nie z substratu

**To closes compound interference path equally do M.6:** już 7/7 mechanisms
NEGATIVE for X = e²/2 fundamental derivation.
