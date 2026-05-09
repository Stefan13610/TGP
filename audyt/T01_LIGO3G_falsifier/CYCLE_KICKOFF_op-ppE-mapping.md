---
title: "Cycle kickoff brief — op-ppE-mapping (Path B; closes T01 N2)"
date: 2026-05-07
parent: "[[README.md]]"
type: cycle-kickoff
tgp_owner: audyt/T01_LIGO3G_falsifier
target_cycle: research/op-ppE-mapping/
tags:
  - kickoff
  - cycle-brief
  - ppE
  - 2PN-phase
  - M911
  - SPA
  - inspiral
  - T01
  - EXT-5
  - path-B
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[CONVENTION_DECISION.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
---

# Cycle kickoff — `research/op-ppE-mapping/`

> **Cel pliku.** Dostarczyć **gotowy do skopiowania** brief
> uruchomieniowy dla cyklu `research/op-ppE-mapping/`, który ma
> zamknąć ścieżkę B z [[README.md]] (mapowanie (5/6) U³ → ppE
> phase parameters). Brief jest *gotowy do startu* — autor (lub
> nowy subagent) może uruchomić cykl bez dodatkowej preparacji
> kontekstu.
>
> **Polityka:** ten plik **przygotowuje** cykl, ale **nie tworzy go**
> — utworzenie folderu `research/op-ppE-mapping/` z plikami inicjalnymi
> należy do autora (T01 ma `may_edit_core: false`, w tym brak
> uprawnień do tworzenia nowych folderów `research/op-…/`).

## §1 — Charakterystyka cyklu

| Pole | Wartość |
|------|---------|
| **Folder docelowy** | `research/op-ppE-mapping/` |
| **Klasa** | analytical-derivation |
| **Domyka** | T01 N1, N2 (z [[NEEDS.md]]) |
| **Otwiera bridge** | tak — `op-LIGO-3G-deviation/` używa output Phase 1.4 |
| **Estymata pracy** | ~2–3 sesje analityczne (1 sesja = 4–6 godzin pracy + sympy) |
| **Output główny** | `Phase1_results.md` — β_ppE^TGP^(2PN-phase, b=−1) liczbowo zamknięty + multi-coefficient pattern |
| **Output secondary** | `Phase2_validation.md` — sympy LOCK 5/5 (style M9_1_pp_P1) + literature cross-check |

## §2 — Cele i kryteria sukcesu

### 2.1 Phase 1 (analytical core)

**Cel:** Pełne wyprowadzenie β_ppE^TGP^(b=−1) z M9.1''
two-body Lagrangian + dE/dt luminosity + SPA inversion.

**Kryteria PASS (sympy LOCK):**
- [P1.1] Two-body Lagrangian w M9.1'' policzony do v⁸ (4PN orbital).
- [P1.2] E_orb(v) reprodukuje GR przy v² (Newton), v⁴ (1PN) — zgodne
  z β_PPN = γ_PPN = 1 EXACT (z M9_1_pp_verify.py).
- [P1.3] E_orb(v) odbiega od GR przy v⁶ (2PN orbital), zgodnie z
  metric (5/6) U³ — **kluczowy test spójności**. Sympy weryfikacja:
  współczynnik deviation w E_orb(v⁶) musi być pochodną
  z (5/6) przez chain rule (oczekiwana skala O(1)·(5/6)).
- [P1.4] dE/dt do v¹⁸ (4PN luminosity); jeśli M9.1'' modyfikuje
  quadrupole formula (założenie A2/A3 z [[PPN_TO_PPE_MAPPING.md]] §5),
  to explicit policzyć.
- [P1.5] SPA inversion → Ψ_TGP(f) z explicit β_(N-PN_phase)
  współczynnikami dla N = 2, 3, 4 (phase convention).
- [P1.6] β_ppE^TGP^(2PN-phase, b=−1) = liczba ± 30%, ±5% target po
  Phase 2 numerical cross-check.

**Kryterium FAIL** (wymaga rewizji założeń A1–A6 z
[[PPN_TO_PPE_MAPPING.md]] §5):
- E_orb(v⁶) NIE odbiega od GR pomimo (5/6) U³ w metric → A1 (two-body
  Lagrangian) wymaga rewizji.
- E_orb(v⁶) odbiega o O(2)·(5/6) zamiast O(1)·(5/6) → A2 (quadrupole
  formula) wymaga modyfikacji w M9.1''.
- Multi-coefficient pattern (β_2PN, β_3PN, β_4PN) NIE jest pochodną
  z {(5/6), (23/12), (19/6)} → konwencja PN counting wymaga rewizji
  (CONVENTION_DECISION może wymagać aktualizacji).

### 2.2 Phase 2 (literature cross-check & uniqueness)

**Cel:** Pokazać, że TGP-pattern (β_2PN, β_3PN, β_4PN) jest
*specyficzny* dla M9.1'' i nieproducowalny z innymi modyfikacjami GR
bez fitting.

**Kryteria PASS:**
- [P2.1] Catalog Yunes–Yagi–Pretorius 2016 + dCS, sGB, EÆ, Brans-Dicke
  — żaden NIE reprodukuje TGP pattern z bare coefficients.
- [P2.2] Najlepsze fit z 2-parameter modyfikacji GR ma residual >> 5σ
  od TGP pattern w ET-D + CE detection regime.
- [P2.3] Sympy LOCK 5/5 (struktura testów jak M9_1_pp_P1) na
  cross-coefficients ratios.

### 2.3 Phase 3 (paper-ready output)

**Cel:** Output w formie `Phase3_paper-ready.md` z:
- Tabela TGP β_(N-PN_phase) z N = 2, 3, 4, 5, 6.
- Comparison plot z catalog (10–15 modyfikacji GR).
- Falsifier statement liczbowy (zaktualizowany [[FALSIFIER_STATEMENT_DRAFT.md]]).

## §3 — Setup analityczny (gotowy do skopiowania jako Phase 1.0)

### 3.1 M9.1'' two-body Lagrangian — schemat

Punkt wyjścia: M9.1'' static metric (M9_1_pp_setup.md §1):

```
ds² = -c₀² · (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j
ψ = 1 + ε,        ε(r) = U₁/r + U₂/r² + …  (multi-source decomposition)
```

Two-body uogólnienie: ψ → ψ(x; m₁, m₂) z liniową superpozycją (do
order 1/r²). Linear superposition dla statycznej weak-field jest
dobre — non-linear corrections (DJS Hamiltonian style) wchodzą od
2PN orbital.

```
ψ(x) = 1 + ε(x)
ε(x) = ε₁(|x − x₁|) + ε₂(|x − x₂|) + ε_int(x; x₁, x₂)
ε_i(r) = (G m_i)/(r c²) · [1 + c_2 (Gm_i/(r c²)) + c_3 (Gm_i/(r c²))² + …]
ε_int = O((G m₁/r₁₂)(G m₂/r₂)) — non-linear interaction term
```

**Zadanie Phase 1.1:** policzyć ε_int explicit z Φ-EOM α=2 dla
two-source ansatz. Z (5/6) U³ wynika: ε_int kontrybuuje deviation
do orbital energy E_orb przy v⁶.

### 3.2 Lagrangian formulation

Body Lagrangian z M9.1'' metric (relativistic point particle):

```
L_body_a = -m_a c² √(-g_μν(x_a) v_a^μ v_a^ν / c²)
        = -m_a c² √(f(ψ_a) - h(ψ_a) v_a²/c²)
```

dla a ∈ {1, 2}. Po expansion w v/c i ε:

```
L_body_a ≈ ½ m_a v_a² + (m_a v_a⁴)/(8 c²) − m_a c² · ε(x_a)·(...)
        + O(v⁶/c⁴, ε², ε·v²/c²)
```

(należy careful policzyć prefactors; dla GR Schwarzschild izotropowych
to standardowy DJS — dla M9.1'' identical structure ale z hyperbolic
f).

**Zadanie Phase 1.2:** wyprowadzić full L_2body do v⁸ używając:
- f(1+ε) = 1 - 4ε + 4ε² - 4ε³ + … (z M9_1_pp_setup §2.3 eq:f-Taylor)
- h(1+ε) = 1 + ... (compute from f·h=1 expansion)

Sympy assist: `sympy.series(f(1+e), e, 0, 8) * sympy.series(h(1+e),
e, 0, 8)` z post-processing.

### 3.3 dE/dt luminosity

GR formula (Peters 1964): `dE/dt = -(32/5) G⁴ m₁² m₂² (m₁+m₂)/(c⁵ a⁵)`
dla circular orbit. M9.1'' uogólnienie:

**Założenie A2 (do walidacji):** quadrupole formula struktura zachowana
w M9.1'' przy ψ → 1, deviation O((5/6)U³) wchodzi przez modyfikację
metric perturbation w retarded Green's function.

**Zadanie Phase 1.3:** wyprowadzić dE/dt^TGP do v¹⁸. Reference:
Blanchet *Living Reviews* 2014 §6 — pełna 3.5PN luminosity GR.
TGP modyfikacja: substitute `f^(GR)(ψ) → f^(M9.1'')(ψ)` w retarded
field calc, evaluate by perturbation w (5/6) U³.

### 3.4 SPA inversion

Standard chain (Cutler–Flanagan 1994):
```
v(f) = (πM f)^(1/3)         M = total mass in geometric units
dE/df = (dE/dr)(dr/df)
df/dt = (df/dE)(dE/dt) = -dE/dt / (dE/df)        — sekulárny chirp
t(f) = ∫ df / (df/dt)
Ψ(f) = ∫ 2π t(f) df + … (stationary phase)
```

W M9.1'' modyfikacje:
1. dE/dr → dE^TGP/dr (z Phase 1.2 Lagrangian)
2. dE/dt → dE^TGP/dt (z Phase 1.3)
3. v(f) relation może być modyfikowana przez Kepler-modify w M9.1''
   (do sprawdzenia; założenie A4 — likely zachowana adiabaticity).

**Zadanie Phase 1.4:** wykonać full SPA inversion w M9.1''. Output:
explicit współczynniki β_(N-PN_phase) w fazie waveformu.

## §4 — Wymagania techniczne

### 4.1 Tools

- **sympy** — algebraic manipulation (jak M9_1_pp_P1).
- **Optional: Mathematica** — niektóre PN expansions z DJS są łatwiejsze
  w Math (ale sympy też wystarczy).
- **Reference codes:** TaylorF2 implementation w `lalsimulation`
  (LIGO public software) — do cross-check.

### 4.2 Skill set (subagent)

- Post-Newtonian theory (Blanchet *Living Reviews* 2014 fluently).
- ppE framework (Yunes–Pretorius 2009 + Sampson–Yunes–Cornish 2013).
- sympy proficiency (>200 lines symbolic algebra).
- M9.1'' baseline (M9_1_pp_setup, M9_1_pp_P1 — already in context).

### 4.3 Estymowane sesje

| Phase | Praca | Szacunkowo godzin |
|-------|-------|-------------------|
| Phase 1.1 (two-body Lagrangian setup) | DJS-like derivation | 2–3 h |
| Phase 1.2 (E_orb to v⁸) | sympy + cross-check GR | 3–4 h |
| Phase 1.3 (dE/dt to v¹⁸) | most complex | 4–6 h |
| Phase 1.4 (SPA inversion → β_(N-PN_phase)) | mechanical | 1–2 h |
| Phase 1.5 (multi-coefficient pattern verification) | sympy LOCK | 1–2 h |
| Phase 2 (literature cross-check) | catalog comparison | 2–3 h |
| Phase 3 (paper-ready output) | formatting | 1–2 h |
| **Razem** | | **~14–22 h, ~2–4 sesje** |

## §5 — Walidacja i gates

### 5.1 Pre-commit gate (Phase 0 balance check)

Zgodnie z **CALIBRATION_PROTOCOL § ABSOLUTE BINDING gate**
(M03 Phase 6 enforcement post-2026-05-06):

Cykl `op-ppE-mapping/` MUSI utworzyć `Phase0_balance.md` zawierający:
- inputs: lista plików read (M9_1_pp_setup, M9_1_pp_P1, sek08c,
  sek08a, ten kickoff, [[CONVENTION_DECISION.md]])
- outputs: lista plików/predykcji którymi cykl będzie kontrybuował
- predyktywność ratio (jeśli dotyczy): N_outputs / N_locked_inputs
- epistemic class tag (per-row): DERIVED, ANSATZ, or NUMERICAL

### 5.2 Output validation gates

| Gate | Test | Pass criterion |
|------|------|----------------|
| G1 | Sympy LOCK on c_n consistency | 5/5 PASS (M9_1_pp_P1 baseline reproduced) |
| G2 | E_orb(v⁴) = E_orb^GR(v⁴) | EXACT match (β_PPN=1 verification) |
| G3 | E_orb(v⁶) deviation from GR | proportional to (5/6) within numerical |
| G4 | β_(2PN-phase) numerical | within OOM of [[PPN_TO_PPE_MAPPING.md]] §2.2 estimate |
| G5 | Multi-coefficient (β_2PN, β_3PN, β_4PN) | all three computed, ratios reproducible |
| G6 | TaylorF2 cross-check | M9.1'' reduces to GR when (5/6)→0 (sanity) |

## §6 — Output format (proponowany)

```
research/op-ppE-mapping/
├── README.md                    — overview, status, cycle log
├── Phase0_balance.md            — M03 gate
├── Phase1_setup.md              — analytical setup, this kickoff adapted
├── Phase1_two_body_lagrangian.md       — §3.1, §3.2
├── Phase1_orbital_energy.md            — §3.3 with sympy outputs
├── Phase1_luminosity.md                — §3.4
├── Phase1_SPA_inversion.md             — §3.5
├── Phase1_results.md            — β_ppE^TGP locked, multi-coefficient pattern
├── Phase2_literature_crosscheck.md
├── Phase3_paper_ready.md        — sec1: setup, sec2: results, sec3: comparison
├── *.py                         — sympy scripts (jak M9_1_pp_p1_higher_pn.py)
├── *.txt                        — sympy outputs
└── NEEDS.md                     — open follow-ups (Q4, Q5 z T01 NEEDS)
```

## §7 — Output → T01 closure path

Po zamknięciu `op-ppE-mapping/` Phase 1:

1. **Update** [[FALSIFIER_STATEMENT_DRAFT.md]] §1: zastąpić `[β_th]`
   liczbą z `Phase1_results.md`.
2. **Update** [[SENSITIVITY_BACK_OF_ENVELOPE.md]] §4.2: zaktualizować
   "OOM" do "locked-from-Phase1" w kolumnie β_ppE^TGP.
3. **Update** [[NEEDS.md]] N1, N2: zamknąć status "OPEN" → "EXECUTED via op-ppE-mapping".
4. **Update** [[README.md]] sekcja "Postęp domknięcia" tabela:
   "PPN_TO_PPE_MAPPING preview" → "Phase B EXECUTED, locked numerically".
5. **(Opcjonalnie autor)** Wkleić zaktualizowany falsifier do
   `PREDICTIONS_REGISTRY.md` jako wpis `M911-P1` z liczbowym threshold.

T01 status awansuje wtedy z **PREPARED** do **CLOSED-PARTIAL
(B+C executed; A pending)**.

## §8 — Powiązanie z innymi cyklami

- **`op-LIGO-3G-deviation/`** (Path A) — **używa output Phase 1.4 jako input.**
  Po zamknięciu Phase 1 op-ppE-mapping, Path A jest gotowa do startu.
- **`op-newton-momentum/`** (existing) — **input source** (M9_1_pp_P1
  daje c_n=2..7, M9_1_pp_setup daje metric ansatz).
- **S07 (M9.1'' derivation)** — **meta-bloker B1 z [[NEEDS.md]]**.
  S07 closure NIE blokuje op-ppE-mapping (test ansatzu wystarczy
  by liczyć (5/6) → β_ppE^TGP), ale po S07 closure cały M911-P1
  awansuje z (P) do (W).

## §9 — Risk register

| Ryzyko | Prawdopodobieństwo | Impact | Mitigation |
|--------|--------------------|--------| -----------|
| A2/A3 fail (quadrupole formula non-trivially modyfikowana) | medium (30%) | high (rewrite Phase 1.3) | iterate, dopuścić extended timeline |
| TGP β_ppE^TGP exceeds LIGO-O3 bound (deviation już sfalsyfikowana) | low (10%, OOM suggests safe margin) | very high (TGP M9.1'' upada) | proszę autor; honest reporting |
| Multi-coefficient pattern NIE TGP-specific (zdegenerowane z dCS) | medium (20%) | medium (paper weaker) | wzbogacić cross-check Phase 2 |
| Konwencja PN counting wymaga rewizji | low (5%, decision już made) | medium (artifacts patches) | [[CONVENTION_DECISION.md]] §5 OPCJA B fallback |

## Cross-references

- [[README.md]] — diagnoza T01 problemu #3 (Konkurencja modeli)
- [[NEEDS.md]] — N1, N2, Q4, Q5 (closed by this cycle)
- [[CONVENTION_DECISION.md]] — adopted phase-PN convention
- [[PPN_TO_PPE_MAPPING.md]] — preview document, basis dla Phase 1
- [[FALSIFIER_STATEMENT_DRAFT.md]] — target update po Phase 1.4
- [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] — następujący cykl (Path A)
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — input source
- [[../../research/op-newton-momentum/M9_1_pp_setup.md]] — M9.1'' setup
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 0 gate enforcement
- [[../../PREDICTIONS_REGISTRY.md]] — target rejestru after closure

## Bibliografia operacyjna

- Blanchet, *Living Rev. Relativ.* **17**, 2 (2014) — primary reference
- Damour, Jaranowski, Schäfer, Phys. Rev. D **89**, 064058 (2014) — 4PN ADM
- Buonanno, Iyer, Ochsner, Pan, Sathyaprakash, Phys. Rev. D **80**, 084043 (2009) — TaylorF2
- Mishra, Iyer, Sundararajan, Phys. Rev. D **93**, 084054 (2016) — 3PN waveform tail
- Sampson, Yunes, Cornish, Phys. Rev. D **88**, 064056 (2013) — ppE consistency
- Yunes, Pretorius, Phys. Rev. D **80**, 122003 (2009) — ppE original
- Will, Wiseman, Phys. Rev. D **54**, 4813 (1996) — original 2PN
