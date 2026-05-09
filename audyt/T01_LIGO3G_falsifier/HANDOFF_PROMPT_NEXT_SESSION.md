---
title: "Handoff prompt dla następnej sesji — T01 follow-up (#2 G_SPA + #3 Production Fisher)"
date: 2026-05-07
parent: "[[README.md]]"
type: handoff-prompt
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - handoff
  - next-session
  - prompt
  - T01
  - G_SPA
  - production-fisher
  - workflow-A
related:
  - "[[README.md]]"
  - "[[FINDINGS.md]]"
  - "[[SESSION_REPORT_2026-05-07_C-B-A-D.md]]"
  - "[[../../research/op-ppE-mapping/Phase1_results.md]]"
  - "[[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]"
  - "[[../../research/op-GWTC3-reanalysis/Phase3_verdict.md]]"
---

# Handoff prompt — następna sesja

> **Cel pliku.** Self-contained prompt dla nowego Claude agenta
> kontynuującego pracę nad T01 LIGO 3G falsifier, gdzie poprzednie
> sesje zamknęły C-light + B + A + D + Tier 5 (GWTC-3). Pozostają:
> #2 G_SPA tighter lock i #3 Production Fisher rerun (+ opcjonalne
> Workflow A pełna implementacja, Workflow C multi-coef ratio test,
> paper revision per Phase3_verdict §3.3).
>
> **Użycie:** wkleić sekcję §1 (PROMPT) jako pierwszy message w
> nowej sesji Claude.

---

## §1 — PROMPT (do wklejenia w nowej sesji)

```
Jesteś ekspertem w fizyce teoretycznej (post-Newtonian gravitational
wave theory, ppE framework, Fisher matrix forecasting). Twoje
zadanie: kontynuować pracę nad T01 LIGO 3G falsifier w vault TGP.

## Kontekst (z poprzednich sesji 2026-05-07)

TGP (Theory of Geometry of Pulsation) to scalar-tensor program
z emergentną metryką M9.1'' = -c²(4-3ψ)/ψ dt² + h(ψ) δ_ij dx^i dx^j.
M9.1'' reprodukuje GR exactly do 1PN (β_PPN = γ_PPN = 1) ale
deviates przy 3PN-energy / 2PN-phase: |Δg_tt|/c² = -(5/6) U³.

W poprzednich sesjach wykonano (zacznij od read-once tych plików
w kolejności):

1. **[[audyt/T01_LIGO3G_falsifier/SESSION_REPORT_2026-05-07_C-B-A-D.md]]**
   — pełna execution C-light → B → A → D (sympy LOCK 14/14, β_ppE^TGP^(b=-1)
   = -5/64 ≈ -7.81·10⁻² LOCKED, Fisher matrix forecasts, paper draft).

2. **[[research/op-ppE-mapping/Phase1_results.md]]** — sympy LOCK
   β_ppE^TGP + multi-coefficient ratios {-23/10, -38/23, +337/228}.

3. **[[research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]**
   — calibrated Fisher thresholds (CE single-event decisive, ET-D
   stack ~10 events, LIGO-O5 stack ~250 events).

4. **[[research/op-GWTC3-reanalysis/Phase3_verdict.md]]** — Tier 5
   GWTC-3 reanalysis EXECUTED. **Verdict: TGP CONSISTENT z GWTC-3
   within 1σ; BF ≈ 0.97 INCONCLUSIVE.** KEY DISCOVERY: detection
   wymaga **dedicated TGP-specific Bayes pipeline** (single-coef
   prior); generic LIGO ToGR multi-coef marginalized factor ~50×
   weaker.

5. **[[audyt/T01_LIGO3G_falsifier/FINDINGS.md]]** — top-level
   audit synthesis (+ §2.5 Tier 5 EXECUTED note).

## Stan locked (READ-ONLY)

- α=2 vacuum Φ-EOM kinetic exponent (LOCKED, F2 w PREDICTIONS_REGISTRY)
- f(ψ) = (4-3ψ)/ψ M9.1'' canonical form (LOCKED z M9_1_pp_setup §2)
- c_n = a_n/a_1^n współczynniki (sympy LOCK 5/5: -1, +5/3, -10/3,
  +22/3, -154/9, +374/9 dla n=2..7)
- α_n^TGP w g_tt: 1, -2, +2, -7/3, +35/12, -91/24, +91/18 (LOCKED 7/7)
- Δα_n = α^TGP - α^GR: 0, 0, 0, **-5/6**, +23/12, -19/6, +337/72 (LOCKED 7/7)
- β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81·10⁻² (η=1/4, G_SPA=1 central)
- OOM window: |β_ppE^TGP| ∈ [5.5·10⁻², 1.2·10⁻¹]
- Multi-coefficient ratios (M911-P2): {-23/10, -38/23, +337/228}
- M911-P1, P2, P3 wpisy w PREDICTIONS_REGISTRY.md Sector 2b LIVE

## Twoje zadanie

Wykonaj kolejno (priorytet w tej kolejności):

### Zadanie #2 — G_SPA tighter lock (op-ppE-mapping Phase 1.5)

**Cel:** zmniejszyć OOM uncertainty G_SPA z ~30% do ~5% precision.

**Background:** β_ppE^TGP = -(3/(128 η))·(5/6)·G_SPA, gdzie G_SPA
jest "SPA chain prefactor" z dE/dr → dE/dt → SPA chain. Phase 1
założyło G_SPA = 1.0 jako central (OOM justified przez "metric-only
deviation" argument w analogicznych modyfikacjach GR — Sampson-Yunes-
Cornish 2013). Phase 1.5 wymaga *explicit* derivation.

**Wymaga:**
1. Two-body Lagrangian w M9.1'' do v⁸ (4PN orbital). Reference: DJS
   (Damour-Jaranowski-Schäfer 2014, Phys. Rev. D 89:064058) dla GR;
   uogólnienie na hyperbolic f(ψ) wymaga substitution `c² → c²·(4-3ψ)²/ψ²`.
2. Modified quadrupole formula w M9.1'' (retarded Green's function w
   hyperbolic metric). Reference: Blanchet *Living Reviews* 17:2 (2014)
   §6 dla GR; M9.1'' uogólnienie wymaga retarded scalar field equation
   solution.
3. Pełen SPA chain z modified E_orb(v) i dE/dt(v) → Ψ(f) → β_ppE^TGP
   liczbowy do 5% precision.

**Output:**
- `research/op-ppE-mapping/Phase1.5_G_SPA_lock.md` — pełen derivation
- `research/op-ppE-mapping/scripts/phase1_5_G_SPA_derivation.py` — sympy
  + integrals
- Update [[research/op-ppE-mapping/Phase1_results.md]] §6 z locked
  G_SPA value (zastąpić "G_SPA = 1.0 (PRELIMINARY)" liczbowym lockiem)
- Update downstream: FALSIFIER_STATEMENT_DRAFT, PREDICTIONS_REGISTRY
  M911-P1, paper draft

**Estymata:** ~1 sesja analytic (4-6 godz pracy + sympy)

**Założenia (z [[audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]] §5):**
- A1 (linear superposition two-body): VALID weak-field
- A2 (quadrupole formula structure): **MUSI walidować** w Phase 1.5
- A3 (dE/dt 2PN+ propagation): **MUSI walidować**
- A4 (SPA adiabatic): VALID
- A5 (no new radiative DOF): VALID via GW1, GW2 w rejestrze
- A6 (PN counting convention): VALID PHASE primary

**Pass criteria:**
- G_SPA central locked z ±5% precision (vs current ±30% OOM)
- Założenia A2, A3 explicitly validated (lub modified quadrupole
  formula explicit)
- Sympy LOCK 5/5 PASS

### Zadanie #3 — Production Fisher rerun (bilby/pycbc)

**Cel:** zastąpić analityczne ASD fits w op-LIGO-3G-deviation z
production-grade noise curves + bilby Fisher matrix; calibrowane
absolute β bounds z 5% precision (vs current factor ~3-5x off).

**Background:** Path A Phase 2 użyła analytical ASD fits (LIGO-O5
A+, ET-D, CE) które są factor ~3-5× off od production-grade.
SNR estimates underestimate ~3-5× literature values (np. moje
ET-D GW150914-like SNR ~50, literatura ~250).

**Wymaga:**
1. Pobrać oficjalne noise curves:
   - LIGO A+: DCC T1800042 (`aLIGOAPlusDesignSensitivityT1800042.txt`)
   - ET-D: ETDSensitivityCurveT0900288 (Hild 2010 design study)
   - CE: cosmicexplorer.org/data/ce1_strain.txt
2. `pip install bilby pycbc lalsimulation` (jeśli nie zainstalowane)
3. Implementować TaylorF2 + ppE deviation w bilby:
   - Use `bilby.gw.source.lal_binary_black_hole` + custom phase
     correction `δΨ = β_ppE · u^(-1)` przy b_ppE = -1
   - Lub zmodyfikować `pycbc.waveform` z ppE inflation
4. Rerun Fisher matrix dla 5 scenariuszy z Phase 2 + dodatkowo:
   - GW170817-like BNS (Λ_tidal = 0)
   - Heavy BBH M_tot = 100 M_⊙
   - Asymmetric q = 4 BBH
5. Validation gates (z CYCLE_KICKOFF_op-LIGO-3G-deviation §5.2):
   - G1: GR Fisher LIGO-O3 GW150914 SNR ~24 (within 5% literature)
   - G2: LIGO O3 ToGR ppE bound (within 30%)
   - G3: ET-D Fisher GW150914 single SNR ~250 (within 5%)
   - G4: β_TGP detection scenario clean 5σ
   - G5: degeneracy with χ_eff < 80%
   - G6: stack √N internal consistency

**Output:**
- `research/op-LIGO-3G-deviation/Phase2.5_production_fisher.md`
- `research/op-LIGO-3G-deviation/scripts/phase2_5_bilby_fisher.py`
- Update [[research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]
  §2 tabela z production numbers (zamiast calibrated OOM)
- Update downstream artifacts

**Estymata:** ~2-3 sesje numerical

**Caveat:** wymaga internet access dla download noise curves. Jeśli
brak, użyć cached lub fallback do Hild 2010 published values
+ calibration jak w Phase 2.

### Zadania opcjonalne (dependent na czasie)

#### R1 — Workflow A (TGP-specific Bayes z public GWTC-3)

**Cel:** pełna implementacja Workflow A z [[research/op-GWTC3-reanalysis/Phase3_verdict.md]]
§3.1: TGP-specific single-coefficient prior (delta function at β_TGP =
-5/64) Bayes inference na public GWTC-3 events.

**Estymata:** ~3-5 sesji (najwyższy ROI dla pre-publication signal).

**Expected output:** log10(BF_TGP/GR) z 90 BBH events; ~1-2σ
tentative signal jeśli TGP correct.

#### Workflow C — multi-coefficient ratio test (M911-P2)

**Cel:** test consistency ratios {-23/10, -38/23, +337/228} z
published GWTC-3 multi-coefficient posteriors. Stronger niż Workflow
B, weaker niż A.

**Estymata:** ~2-3 sesje.

#### Paper revision (Path D §4-§5)

**Cel:** zaktualizować [[papers/M911_LIGO3G_paper/paper_draft.md]]
sekcje 4-5 zgodnie z [[research/op-GWTC3-reanalysis/Phase3_verdict.md]]
§3.3 — clarify że detection thresholds wymagają TGP-specific Bayes
pipeline, nie generic LIGO ToGR.

**Estymata:** ~1 sesja.

## Polityka i zakres

- **Scope:** możesz edytować w `audyt/T01_LIGO3G_falsifier/`,
  `research/op-ppE-mapping/`, `research/op-LIGO-3G-deviation/`,
  `research/op-GWTC3-reanalysis/`, `papers/M911_LIGO3G_paper/`,
  `PREDICTIONS_REGISTRY.md`. NIE edytować core/ bez explicit
  autoryzacji.
- **M03 gate enforcement:** każdy nowy cykl wymaga Phase 0 balance
  (zob. CALIBRATION_PROTOCOL).
- **Honest reporting:** jeśli verdict negatywny lub niepewny, raportuj
  honestnie (nie udawaj że wynik jest pozytywny).
- **TodoWrite:** używać dla multi-step tasks (>2 kroki).

## Recommended start

1. Read SESSION_REPORT_2026-05-07_C-B-A-D.md (overview).
2. Read FINDINGS.md (audit synthesis).
3. Read Phase 3 verdict files (3 cycles): op-ppE-mapping, op-LIGO-3G-
   deviation, op-GWTC3-reanalysis.
4. Decide priorytet: zwykle G_SPA tighter lock (#2) **przed**
   production Fisher (#3), bo β_TGP central value enters in Fisher.
5. Spawn subagents jeśli zadania są parallelizable (np. literature
   search dla DJS-style two-body Lagrangian może iść w tle).

Zaczynaj od Zadania #2 chyba że user explicitly powiedział inaczej.
```

---

## §2 — Notes dla autora przy uruchomieniu nowej sesji

### Kontekst którego prompt nie zawiera

- **Vault root path:** `C:\Users\Mateusz\Documents\ObsydnianMain\TGP\TGP_v1\`
  (z system reminder w nowej sesji)
- **Email:** mateusz.serafin@traveltech.pl (z system reminder)
- **Date:** w sesji 2026-05-07 ten audit-folder + cykle są na świeżo

### Jeśli nowa sesja idzie w niewłaściwym kierunku

Reset z explicit instrukcją:
> "Zacznij od read [[SESSION_REPORT_2026-05-07_C-B-A-D.md]] +
> [[FINDINGS.md]]. Pracuj WYŁĄCZNIE w T01-related folders i scope
> zdefiniowanym w prompt §1. Zaczynaj od Zadania #2."

### Jak walidować że nowa sesja działa correctly

Po pierwszej response sprawdź:
1. Czy agent przeczytał kluczowe pliki (linki w sesji powinny być
   wikilinks).
2. Czy używa TodoWrite dla multi-step tracking.
3. Czy plan agentowy zaczyna od Zadania #2 (G_SPA), nie #3 (Production).
4. Czy honest reporting jest aktywny (nie ukrywa niepewności).

### Optional: rozpoczęcie z subagent kickoff

Jeśli chcesz zacząć z parallel research, prompt może być rozszerzony:
> "Spawn subagent (Explore typu) dla literature search:
> 'Find references for Damour-Jaranowski-Schäfer 4PN ADM Hamiltonian
> in Schwarzschild metric (Phys. Rev. D 89:064058 2014); identify
> the key derivation chain dE/dr → dE/dt → SPA in eq. (X.Y); collect
> cross-references to similar derivations in modified-gravity
> theories (Brans-Dicke, dCS, sGB)'."

To może iść równolegle z głównym agentem rozpoczynającym Zadanie #2
analytic setup.

## §3 — Cross-references

- [[README.md]] — T01 audit
- [[SESSION_REPORT_2026-05-07_C-B-A-D.md]] — pełna execution
- [[FINDINGS.md]] — audit synthesis
- [[../../research/op-ppE-mapping/Phase1_results.md]] — β_ppE^TGP locked (Phase 1.5 będzie tighter)
- [[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] — calibrated thresholds (Production Fisher będzie precise)
- [[../../research/op-GWTC3-reanalysis/Phase3_verdict.md]] — Tier 5 verdict + workflow caveat
- [[../../papers/M911_LIGO3G_paper/paper_draft.md]] — paper draft (revision pending per Phase3_verdict §3.3)
- [[CYCLE_KICKOFF_op-ppE-mapping.md]] — original Path B kickoff (rozszerza Phase 1.5)
- [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] — original Path A kickoff (rozszerza Phase 2.5)
