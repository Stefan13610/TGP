---
title: "L05 — wykładnik masy: k=4 (LP-4) vs p=5−α (R3 2026-05-01)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L05_mass_exponent_drift
tags:
  - audit
  - ontology
  - mass-formula
  - lepton-masses
  - exponent-drift
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../L04_ODE_dualism_alpha]]"
  - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]]"
  - "[[../../research/why_n3/CORRECTIONS_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["mass-formula-reconcile"]
  depends_on:
    - "[[../L04_ODE_dualism_alpha]]"
  impacts: []
  source_of_status:
    - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4"
    - "[[../../research/why_n3/CORRECTIONS_2026-05-01.md]]"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L05 — wykładnik masy: k=4 vs p=5−α

## 🟢 STATUS UPDATE 2026-05-16 — **CLOSED-RESOLVED A−**

**Closure source:** [[../../research/op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase_FINAL_close.md]]

**Verdict:** STRUCTURAL_DERIVED_NATIVE (A−), 12/12 sympy PASS, 11 FP (91.7%), 6/6 P-requirements
RESOLVED. Możliwość A constructively confirmed; Możliwości B i C eliminated.

**Key analytical results:**
- `k_full(α, d) = 4 + d(α-2)/2` — volumetric M_full z Derrick virial scaling
- `σ_match(α, d) = 1 + (d-1)(α-2)/4` — core-tail matching A_tail ∝ A^σ_match
- `k_obs(α, d=3) = 5 − α = p_crit_Sobolev(d=3) − α` — structural d=3 conformal critical

**Reconciliation theorem:**
- LP-4 "M ∝ A^4" = m_obs(α=1, d=3) = 4 (NOT k_full=5/2)
- R3 "m_obs ∝ A_tail^(5-α)" = k_obs(α, d=3) z Sobolev p_crit structure
- m_obs ≠ M_full operationally formalized (ADM-vs-Komara analog)

**Audit Możliwości dispositioned:**
- Możliwość A (LP-4 dotyczy M_full, R3 dotyczy m_obs; LP-4 niedoprecyzowane) — ✅ **CONFIRMED constructively**
- Możliwość B (R3 fitting artifact) — ❌ ELIMINATED (Sobolev structural origin)
- Możliwość C (LP-4 wrong) — ❌ ELIMINATED (LP-4 correct dla m_obs at α=1)

**Downstream impact:**
- audyt/PRIORITY_MATRIX klaster D L05 P2 OPEN → CLOSED-RESOLVED
- audyt/L08 (kink fermion closure) — m_obs vs M_full distinction available dla emergent Dirac pole-mass
- research/why_n3/CORRECTIONS_2026-05-01.md — analytical backbone added (m_obs ≠ M_full now derived)
- research/mass_scaling_k4 — reinterpreted (LP-4 = m_obs at α=1, NOT M_full)
- core/sek08b_ghost_resolution thm:B1'' — reinterpretation note recommended (specific to M_full)

---

## Klasa: LUKA ONTOLOGICZNA / NIEJAWNY PIVOT • Priorytet: ~~P2~~ **CLOSED-RESOLVED 2026-05-16**

## Diagnoza

Wykładnik wzoru masy `M ∝ A_tail^k` był pierwotnie **k=4**, zamknięty
w LP-4 (PLAN_DOMKNIECIA_MASTER) jako **9/9 PASS** z trzema niezależnymi
argumentami:

1. **Zero-mode** (E_kin = E_pot, virial theorem): E₂ = 0, eliminuje O(A²)
2. **Wymiarowa konwergencja**: `k = 2(d-1)/(d-2) = 4`, jedyne całkowite k w d=3
3. **Dyskryminacja**: k=3 → r₂₁ = 55, k=4 → 207, k=5 → 784 — tylko k=4
   trafia w PDG [200, 210]

Dodatkowo „twierdzenie B1''" w core: *„k=4 jedyne całkowite k dające
skończoną masę solitonu w d=3"*.

## Pivot R3 (2026-05-01) — niejawne unieważnienie

[[../../research/why_n3/CORRECTIONS_2026-05-01.md]] wprowadza:

```
m_obs = c · A_tail^(5−α)
```

gdzie:
- α=1 (substratowa K=g²) ⇒ p=4 (zgodne z LP-4)
- α=2 (kanoniczna K=g⁴) ⇒ p=3

Numeryczna weryfikacja w `r3_observable_vs_full_mass.py`:

| α | p (numer.) | 5−α | diff |
|---|-----------|------|------|
| **1.00** | **4.001** | **4.000** | **0.0%** |
| **2.00** | **3.001** | **3.000** | **0.0%** |

Strukturalna interpretacja w R3:

> Twoja distinction między **obserwowalną a pełną masą** jest analogią do:
> - GR: masa ADM (asymptotyczna) ≠ masa Komara (wewnętrzna)
> - QFT: bare mass (UV) ≠ renormalized mass (IR)
>
> W R3:
> - **m_obs = c·A_tail^(5−α)** = "wagę z dystans" przez tail-coupling
> - **M_full** (= K + V_eff) = pełna energia struktury wewnętrznej

R3 deklaruje pełną zgodność dla TGP-canonical α=2:

| Wielkość | TGP α=2 (z A³) | PDG | Diff |
|----------|----------------|-----|------|
| m_μ/m_e | 206.56 | 206.77 | −0.099% ✓ |
| m_τ/m_e (z Koide) | 3474.28 | 3477.23 | −0.085% ✓ |

## Strukturalna sprzeczność

LP-4 argument konwergencyjny `k = 2(d-1)/(d-2) = 4` jest **niezależny
od α** (od formy K(φ)). Wynika z analizy zachowania się asymptotyk
solitonu w d=3, niezależnie od kinetic prefactor.

R3 twierdzi `p = 5 − α`, co dla α=2 daje p=3. **Argument konwergencyjny
mówi, że jedyne całkowite k w d=3 to k=4** — więc `p=3` jest *niedopuszczalne*
zgodnie z LP-4.

Możliwości:

### Możliwość A: argument konwergencyjny LP-4 dotyczy `M_full`, nie `m_obs`

R3 explicit: `m_obs ≠ M_full`. Może argument konwergencyjny `k=4` dotyczy
*pełnej masy* (M_full ∝ A^4), a `m_obs ∝ A^3` jest *projekcją* (tail coupling).

→ Wówczas LP-4 jest **niedoprecyzowane**: nie różnicuje `m_obs` od `M_full`.

### Możliwość B: R3 jest błędny dla α=2

Empiryczna zgodność `p=3` dla α=2 (z dokładnością 0.099%) może być
artefaktem dopasowania przez Koide K=2/3 + zerowy parametr g₀^τ. R3 nie
wyprowadza p=3 — *dopasowuje*.

→ Wówczas LP-4 (k=4) zostaje, R3 jest *niezamierzonym pivotem* przez fitting.

### Możliwość C: argument konwergencyjny LP-4 jest błędny

`k = 2(d-1)/(d-2)` może wymagać formy K(φ)·∂² operatora kanonicznego
(K=const). Dla niestandardowego K(φ)=φ⁴ ten argument może nie obowiązywać.

→ Wówczas LP-4 jest **specyficzne dla K=g²**, a dla K=g⁴ wzór jest inny.

**Aktualny stan: niejasny.** Każda możliwość wymaga formalnego
rozstrzygnięcia.

## Status w plikach

- LP-4 9/9 PASS z `lp4_mass_exponent_verification.py` używającym
  K=g² (substratowa) — *więc faktycznie sprawdza k=4 dla α=1*
- R3 `r3_observable_vs_full_mass.py` używa parametryzacji `p` jako
  numerical fit dla danego α — *więc nie jest niezależną pochodną
  wzoru masy*
- README explicit: „K(φ)=φ⁴ (z α=2 selection)" + „Mass formula M ∝
  A_tail^4" — *te dwie deklaracje są wzajemnie sprzeczne*

## Status w audycie

[[../../meta/AUDYT_TGP_2026-05-01.md]] nie wymienia wprost. R3
CORRECTIONS sam deklaruje: *„Sprzeczność audytu A1+A2 jest rozwiązana"*
— ale to oznacza tylko, że *jakieś* sprzeczność strukturalna jest
adresowana, nie że L05 jest zamknięty. Audit § A.1 dotyczył metryki, nie
wzoru masy.

## Rekomendacja

### Phase 1 — formal rozstrzygnięcie argumentu konwergencyjnego LP-4

Re-derive `k = 2(d-1)/(d-2)` dla **niestandardowego K(φ)**. Krok po kroku:

1. Dla `K(φ)·∂²φ + V'(φ) = ρ`, asymptotyka `φ → 1 + δ`:
   linearyzowane EOM ma rozwiązanie tail `δ ∝ exp(-m·r)/r`.
2. Energia całkowita `E = ∫(½K(φ)·(∂φ)² + V(φ))d³x` musi być skończona.
3. Pokazać, że dla `K(φ) = K_geo·φ^α` warunek skończonej energii daje
   konkretną relację `k(α, d)`.

Jeśli wyjdzie `k = (4-α)·(d-1)/(d-2)` lub podobne — to ujednolica LP-4
i R3.

### Phase 2 — distinction `m_obs` vs `M_full`

Sformalizować w sek08 czym dokładnie jest `M_full` (pełna energia
solitonu) vs `m_obs` (asymptotyczne tail coupling). Pokazać, że:

```
M_full = c₁ · A^k_full
m_obs = c₂ · A^k_obs
k_full ≠ k_obs (chyba że α ma specific value)
```

I konkretnie wyznaczyć `k_full(α, d)` i `k_obs(α, d)`.

### Phase 3 — sync z decyzją L04

Po decyzji L04 (α=1 czy α=2):

- Jeśli α=1: k_obs = 4 (LP-4 oryginalne) potwierdzone, R3 cykl `p=4` weryfikujący
- Jeśli α=2: k_obs = 3 (R3 nowy) — LP-4 wymaga revisji jako specyficzne dla K=g²

**Estymata:** 2–3 tygodnie (po L04 decyzji).

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `core/sek08b_ghost_resolution/sek08b.tex` (twierdzenie B1'') | re-derive k(α, d) |
| `core/sek10_N0_wyprowadzenie/` | sync z formal mass formula |
| `research/why_n3/CORRECTIONS_2026-05-01.md` | sync z formal derivation |
| `research/mass_scaling_k4/` | rename / split na k_obs vs k_full |
| `scripts/lp4_mass_exponent_verification.py` | sync z formal derivation |
| README mass formula deklaracja | reformulacja |

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L5
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4
- [[../../research/why_n3/CORRECTIONS_2026-05-01.md]]
- [[../L04_ODE_dualism_alpha]] (zależność wstępna)
