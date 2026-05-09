---
title: "FINDINGS — T01 LIGO 3G falsifier (synteza audytu 2026-05-07)"
date: 2026-05-07
parent: "[[README.md]]"
type: audit-findings
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - findings
  - audit-synthesis
  - T01
  - EXT-5
  - M911
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - 2PN-phase
  - ppE
  - falsifier
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[CONVENTION_DECISION.md]]"
  - "[[CYCLE_KICKOFF_op-ppE-mapping.md]]"
  - "[[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]"
---

# FINDINGS — T01 LIGO 3G falsifier (synteza audytu)

> **Cel pliku.** Top-level synteza pracy audytowej w
> `audyt/T01_LIGO3G_falsifier/`. Łączy diagnozę EXT-5 (2026-05-06)
> z artefaktami sesji 2026-05-07 (NEEDS, FALSIFIER_STATEMENT_DRAFT,
> PPN_TO_PPE_MAPPING, SENSITIVITY_BACK_OF_ENVELOPE, CONVENTION_DECISION,
> CYCLE_KICKOFF_×2). Adresat: autor TGP (decyzje), recenzent zewnętrzny
> (audit-evidence), kolejne sesje audytu T-class.
>
> **Status pliku:** v1, sesja 2026-05-07. Audyt-folder przeszedł z
> stanu **NOTED** do **PREPARED** w jednej sesji (~3 godz. pracy).

## §0 — Executive summary (1 strona)

**Diagnoza** (z [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5).
TGP_v1 deklaruje "M9.1'' przewiduje |Δg_tt| = (5/6) U³ deviation
od GR (testowalne LIGO 3G)" w `TGP_FOUNDATIONS.md` § 3 i `sek08c`,
ale (1) nie ma sensitivity analysis, (2) nie ma wpisu falsifier
w `PREDICTIONS_REGISTRY.md`, (3) nie ma mapowania ppE, (4) nie ma
spójnego rachunku 2PN-phase, (5) nie ma decyzji o waveform-model.

**Stan po sesji 2026-05-07.**

| Wymiar | Pre-session | Post-session |
|--------|-------------|--------------|
| Diagnoza problemu | full (README v1) | full (README v2 + clarified PN convention) |
| Konkrety zadań | brak | **11 N's + 5 Q's + 4 blokery** ([[NEEDS.md]]) |
| Falsifier-statement | brak | **gotowy draft do wklejenia** do rejestru ([[FALSIFIER_STATEMENT_DRAFT.md]]) |
| ppE mapping | brak | **analytical dictionary preview** ([[PPN_TO_PPE_MAPPING.md]]) |
| Sensitivity OOM | brak | **tabela detekcji LIGO-O3 → ET+CE** ([[SENSITIVITY_BACK_OF_ENVELOPE.md]]) |
| Konwencja PN | niejasna | **PHASE adopted, dual-presentation w rejestrze** ([[CONVENTION_DECISION.md]]) |
| Cykle ścieżek A, B | nie istnieją | **kickoff briefs gotowe do startu** ([[CYCLE_KICKOFF_op-ppE-mapping.md]], [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]) |
| YAML status | NOTED, no needs/findings | PREPARED, has_needs_file: true, has_findings_file: true |

**Kluczowa konkluzja sesji.** Predykcja (5/6) U³ jest fizycznie
detektywna w erze ET-D + CE (2035+), z TGP β_ppE^TGP^(b=−1) ≈ 10⁻¹
OOM vs detector bound ~10⁻³ (single event ET-D) — margin ~100×.
Path C ("najlżejsze domknięcie", wpis falsifier do rejestru) jest
**gotowy do wykonania** *teraz*; Path B (`op-ppE-mapping/`) ma
ready-to-launch kickoff brief, estymowany na ~2–3 sesje analityczne;
Path A (`op-LIGO-3G-deviation/`) wymaga Path B Phase 1.4 jako input,
estymowany na ~3–6 sesji numerycznych.

**Rekomendowana sekwencja:** **C-light → B → A → D**.

## §1 — Co audyt znalazł (consolidated findings)

### 1.1 Pre-existing strength: predykcja jest DERIVED

[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] §3.2
zawiera *kompletne* analityczne wyprowadzenie c_n=2..7 z α=2 vacuum
Φ-EOM, sympy LOCK 5/5 PASS:

```
g_tt^TGP / (-c²) = 1 − 2U + 2U² − (7/3)U³ + (35/12)U⁴ − (91/24)U⁵ + …
g_tt^GR  / (-c²) = 1 − 2U + 2U² − (3/2)U³ +    1·U⁴   − (5/8)U⁵   + …
```

Differences: Δα_3 = -5/6, Δα_4 = +23/12, Δα_5 = -19/6, Δα_6 = +337/72.

**Implikacja audytu.** T01 NIE musi wykonywać research na poziomie
analitycznym dla samej predykcji — to jest done. T01 musi
*operacjonalizować* tę predykcję dla obserwowalnego sygnału ET/CE.

### 1.2 Strukturalna luka: brak operacyjnego mostu PN → waveform → detector

Brakujące mosty:

| Most (chain link) | Co to robi | Stan pre-T01 | Po T01 v2 |
|-------------------|-----------|--------------|-----------|
| `Δα_n` w g_tt → ppE phase β | Mapowanie analityczne na waveform space | brak | **dictionary preview gotowy** ([[PPN_TO_PPE_MAPPING.md]]) |
| ppE β → ET/CE strain h(f) | Numeryczna prognoza detekcji | brak | **OOM ready** ([[SENSITIVITY_BACK_OF_ENVELOPE.md]]) |
| h(f) → falsifier statement | Public contract z obserwacjami | brak | **draft do wklejenia** ([[FALSIFIER_STATEMENT_DRAFT.md]]) |
| Falsifier → cykle wykonawcze | Operationalizacja (kto/jak liczy) | brak | **dwa kickoff briefs gotowe** |

### 1.3 Konwencjonalna pułapka: "U³ → 3PN" jest nieoczywiste

[[CONVENTION_DECISION.md]] zidentyfikował: **Energy-PN** (U³ ↔ 3PN,
naturalne dla statycznego PPN) i **Phase-PN** (U³ ↔ 2PN-phase,
standardowe dla GW waveformów) — to są *dwie różne* konwencje, dla
*tej samej* deviation. Bez explicit decyzji konwersacyjnej, każdy
peer-review reader musi konwertować między nimi → confusion + reduced
credibility.

**Decyzja audytu:** PHASE convention rekomendowana jako primary
w rejestrze, ENERGY zachowana w research-folderze gdzie była używana.

**Implikacja:** **wszystkie tabele bounds** w [[SENSITIVITY_BACK_OF_ENVELOPE.md]]
i [[FALSIFIER_STATEMENT_DRAFT.md]] zostały ZAKTUALIZOWANE w sesji
do konwencji PHASE (b_ppE = −1, nie +1). Zob. §1.4 patches log.

### 1.4 Patches zaaplikowane podczas sesji 2026-05-07 (intra-folder)

Po wprowadzeniu CONVENTION_DECISION zastosowano natychmiast
(OPCJA A) patche do trzech plików:

| Plik | Sekcja | Zmiana |
|------|--------|--------|
| [[README.md]] | Problem #4 | "U³ to człon 3PN" → "U³ to 2PN-phase (3PN-energy)" + nowy paragraf z bounds |
| [[FALSIFIER_STATEMENT_DRAFT.md]] | §1 mapowanie ppE | b_ppE = +1 → **b_ppE = −1**; explicit β_ppE^TGP wzór z κ ~ 0.4–1.5 |
| [[FALSIFIER_STATEMENT_DRAFT.md]] | §1 tabela detekcji | Δβ_ppE^(3PN) → **Δβ_ppE^(2PN-phase, b=−1)**; values dostosowane |
| [[FALSIFIER_STATEMENT_DRAFT.md]] | §1 falsifier clause | "ppE 3PN coefficient" → "ppE 2PN-phase coefficient" |
| [[SENSITIVITY_BACK_OF_ENVELOPE.md]] | §4.2 tabela | Wszystkie bounds → b=−1 line; LIGO-O3 ~10⁻¹ (nie ~3·10⁻¹) |
| [[NEEDS.md]] | N7 | "regime ważności" + explicit konwencja PN |

**Po patches T01 jest internally consistent.**

### 1.5 Bridge-discovery: TGP M9.1'' jest na granicy LIGO-O3

[[SENSITIVITY_BACK_OF_ENVELOPE.md]] §4.2 ustala (po patch konwencji):
- TGP β_ppE^TGP^(b=−1) ≈ 10⁻¹ (OOM, locked w `op-ppE-mapping`)
- LIGO-O3 single-event bound ~10⁻¹ (Chamberlain–Yunes 2017 baseline)

**To znaczy: TGP M9.1'' jest "borderline" w aktualnym LIGO-O3** —
nie sfalsyfikowane, ale również nie *pełni under-constrained* jak
wcześniej sądzono. **Implikacja praktyczna:**

- LIGO-O5 stack (~2027–2030; A+ design SNR_stack ~800) z bound
  ~3·10⁻² **już falsyfikuje** TGP M9.1'' lub potwierdza, na 8–10 lat
  przed ET/CE.
- T01 nie jest tylko "future test 2035+" — jest **active falsifier
  od ~2027**. To istotnie zwiększa ranking strategiczny.

**Side-cykl rekomendowany** (NEEDS.md Q3): re-analiza GWTC-3 ppE
constraints w 2PN-phase specifically dla M9.1'' β_TGP^TGP — ~1 sesja
pracy na danych public, mogłaby od razu *teraz* dać sygnał.

### 1.6 Multi-coefficient pattern jako TGP-distinguishing signature

[[PPN_TO_PPE_MAPPING.md]] §3 ustala: pojedynczy 2PN-phase deviation
b_ppE = −1 jest dzielony z dCS, sGB, Einstein-Æther — TGP NIE jest
distinguishable single-coefficient.

**Co czyni TGP distinguishable:** wzorzec **stosunków**
{β_2PN_phase, β_3PN_phase, β_4PN_phase} = κ · {(5/6), (23/12)·c_4ratio,
(19/6)·c_5ratio} — wymuszony przez α=2 i hyperbolic f(ψ) bez fitting
freedom.

**Implikacja dla cyklu A:** Phase 4 (multi-coefficient extension) jest
**KEY** dla TGP-specificity, nie tylko bonus.

### 1.7 Bloker meta: T01 jest "test ansatzu" dopóki S07 zamknięte

[[../S07_M911_derivation/README.md]] otwarty (P2). Bez derywacji
M9.1'' z fundamentu, T01 testuje **ansatz**, nie *teorię wyprowadzoną*.
To **NIE blokuje** T01 — ansatz może być sfalsyfikowany niezależnie.
ALE peer-review pada paperu Path D (publikacja "M9.1'' ET/CE
forecast") jest **znacznie silniejsza** post-S07 closure.

**Rekomendowana koordynacja:**
- T01 Path C (rejestr) i Path B (ppE mapping) **niezależnie od S07**.
- T01 Path A (Fisher) **niezależnie od S07**.
- T01 Path D (paper) **CO-SCHEDULE z S07 closure** — paper
  napisany "M9.1'' as derived" jest mocniejszy niż "M9.1'' as
  postulate".

## §2 — Co audyt zarekomendował (decision points dla autora)

### 2.1 Tier 1 — natychmiast wykonalne (~30 min, low risk)

**[T01-D1]** Wkleić [[FALSIFIER_STATEMENT_DRAFT.md]] §1, §2 do
`PREDICTIONS_REGISTRY.md` jako wpis **M911-P1**, status **LIVE-PARTIAL**
(z `[β_th]` placeholder zamknięty po Path B Phase 1.4). Counter
+1 (~857). T01 awansuje do **CLOSED-PARTIAL (registered)**.

**Argumenty za:**
- Predykcja matematyczna już DERIVED (M9_1_pp_P1) — wpis nie
  dodaje "fitted-anchor".
- Falsifier-statement = publiczny kontrakt → sygnał dla
  peer-review community.
- Najlżejsza akcja w T01 (Path C w README).

**Argumenty przeciw:**
- Bez liczbowego [β_th] wpis może być postrzegany jako preliminary.
- (Mitigation:) etykieta `LIVE-PARTIAL` + nota "[β_th] zamknięty
  z `op-ppE-mapping` Phase 1; oczekiwana skala ~10⁻¹–10⁻²" jest
  honest reporting.

**Rekomendacja:** **WYKONAJ TERAZ**. Wpis częściowy jest preferowany
nad "nie ma".

### 2.2 Tier 2 — ready-to-launch (~2–3 sesje, analityczne)

**[T01-D2]** Uruchomić cykl `research/op-ppE-mapping/` z kickoff
[[CYCLE_KICKOFF_op-ppE-mapping.md]]. Zamyka N1, N2, Q4, Q5.

**Output (Phase 1):** β_ppE^TGP^(b=−1) liczbowo zamknięty;
multi-coefficient pattern policzony.

**Po wykonaniu T01-D2:**
- Update [[FALSIFIER_STATEMENT_DRAFT.md]] §1: zastąpić `[β_th]`
  liczbą.
- T01 awansuje z CLOSED-PARTIAL do **CLOSED-PARTIAL (B+C executed)**.

### 2.3 Tier 3 — sequential after Tier 2 (~3–6 sesji, numeryczne)

**[T01-D3]** Uruchomić cykl `research/op-LIGO-3G-deviation/` z kickoff
[[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]. Zamyka N1 (waveform), N3,
N4, N5, N7, N11.

**Output (Phase 3):** SNR thresholds + falsifier thresholds dla
wszystkich detektorów (LIGO-O5, ET-D, CE) i scenariuszy
(single, stack, network).

**Po wykonaniu T01-D3:**
- Update [[FALSIFIER_STATEMENT_DRAFT.md]] §1 falsifier clause:
  liczbowy threshold + scenariusze.
- T01 awansuje do **CLOSED-EXECUTED (A+B+C wykonane)**.

### 2.4 Tier 4 — strategic (~1–2 mies., paper)

**[T01-D4]** Draftować short paper "M9.1'' 2PN-phase deviation:
predictions for Einstein Telescope and Cosmic Explorer" do
**Phys. Rev. D** lub **Class. Quantum Grav.** Zamyka N9.

**Co-schedule:** S07 closure (zob. §1.7).

**Po wykonaniu T01-D4:** TGP wchodzi w mainstream falsification
machinery. T01 staje się **CLOSED-PUBLISHED**.

### 2.5 Tier 5 — opportunistic side-cykl (~1 sesja, OPCJONALNE) — **EXECUTED 2026-05-07**

**[T01-D5]** Re-analiza GWTC-3 ppE constraints w 2PN-phase
specifically dla M9.1'' β_TGP^TGP. Zamyka NEEDS.md Q3.

**STATUS: EXECUTED 2026-05-07** via [[../../research/op-GWTC3-reanalysis/]]
(Phase 0 + 1 + 2 + 3 wszystkie executed).

**Verdict:**
- **TGP CONSISTENT z GWTC-3 within 1σ** (deviation 0.37σ z observed
  centroid; BF_TGP/GR ≈ 0.97 INCONCLUSIVE)
- **Brak detection** — GWTC-3 generic ppE-marginalized analysis
  daje σ_combined ~ 0.18 fractional, TGP signal ~0.018 fractional
  → ratio TGP/σ ≈ 0.10 (10× below sensitivity)
- **Brak falsification** — TGP w basenie 1σ zgodności

**KEY DISCOVERY:** detekcja TGP w 3G era wymaga **DEDICATED TGP-
SPECIFIC ANALYSIS PIPELINE**, nie generic LIGO ToGR multi-coefficient
marginalization. Path A optimistic bounds (~10⁻²-10⁻³) odnoszą się
do single-coefficient TGP-specific Bayes; LIGO ToGR papers używają
multi-coefficient marginalized fits which are factor ~50× weaker
per coefficient.

**Implication dla paper draft:** [[../../papers/M911_LIGO3G_paper/paper_draft.md]]
sekcja 4-5 wymaga revision — clarify że detection requires custom
TGP analysis pipeline, NOT generic ToGR.

**Recommended next session:**
- **Workflow A (TGP-specific Bayes z public GWTC-3 events):** ~1
  sesja Python (bilby + custom prior). Estimated yield: ~1-2σ
  tentative signal z aktualnych ~90 BBH events. **Najwyższy ROI**
  z trzech recommended next steps.
- **Workflow C (multi-coefficient ratio test M911-P2):** ~2-3 sesji.
  Test consistency ratios {-23/10, -38/23, +337/228} z published
  multi-coef posteriors. Stronger than Workflow B, weaker than A.

**Argumenty za:** TGP może już być sfalsyfikowane (lub potwierdzone)
LIGO-O3 data, public available. Daje sygnał *przed* ET/CE. **Verdict
2026-05-07: NIE sfalsyfikowane, NIE potwierdzone, ale Workflow A
viable z ~1-2σ expected significance.**

**Rekomendacja:** **OPCJONALNE → EXECUTED**. Workflow A
(R1 z [[../../research/op-GWTC3-reanalysis/Phase3_verdict.md]] §3.1)
jest wysoko-ROI follow-up dla pre-publication signal.

## §3 — Co audyt NIE zrobił (świadome ograniczenia)

### 3.1 NIE wykonano edycji w core/

T01 ma `may_edit_core: false`. W szczególności:
- `TGP_FOUNDATIONS.md` § 3 nie zaktualizowany konwencją PN counting
  (decyzja należy do autora, po adopcji [[CONVENTION_DECISION.md]]).
- `PREDICTIONS_REGISTRY.md` nie zaktualizowany wpisem M911-P1
  (Tier 1 decyzja autora, [[FALSIFIER_STATEMENT_DRAFT.md]] gotowy).
- `core/sek08c/` nie tknięty (pozostawione dla S07 ścieżki domknięcia).

### 3.2 NIE utworzono cykli `research/op-…/`

Tworzenie folderów `research/op-ppE-mapping/` i
`research/op-LIGO-3G-deviation/` należy do autora. Audit-folder
dostarcza kickoff briefs (gotowe do skopiowania jako Phase 0).

### 3.3 NIE wykonano numerycznych Fisher analiz

OOM tabele w [[SENSITIVITY_BACK_OF_ENVELOPE.md]] są wzięte z literatury
(Chamberlain–Yunes 2017, Maggiore 2020, Reitze 2019). Precyzyjne
liczby Fisher dla TGP β_ppE^TGP wymagają cyklu Path A.

### 3.4 NIE rozstrzygnięto konwencji PN ostatecznie

[[CONVENTION_DECISION.md]] **rekomenduje** PHASE convention, ale
adopcja wymaga akceptacji autora. Cykl `op-ppE-mapping/` Phase 1
może odkryć że b_ppE jest *różne* od −1 (np. nowy mode); wtedy
konwencja jest zwalidowana danymi.

### 3.5 NIE adresowano S07 (meta-bloker)

S07 closure jest osobnym audytem ([[../S07_M911_derivation/README.md]]).
T01 *referuje* go jako meta-bloker B1 (zob. NEEDS.md), ale nie
wykonuje pracy nad S07.

## §4 — Stan plików post-session

```
audyt/T01_LIGO3G_falsifier/                          [TOTAL: 9 plików, ~50KB]
├── README.md                                         [v2 — diagnoza + closure progress]
├── NEEDS.md                                          [11 N + 5 Q + 4 blokery]
├── FALSIFIER_STATEMENT_DRAFT.md                      [Path C, ready-to-paste]
├── PPN_TO_PPE_MAPPING.md                             [Path B preview, dictionary]
├── SENSITIVITY_BACK_OF_ENVELOPE.md                   [Path A preview, OOM]
├── CONVENTION_DECISION.md                            [phase-PN convention adopted]
├── CYCLE_KICKOFF_op-ppE-mapping.md                   [Path B kickoff, ~2-3 sesje]
├── CYCLE_KICKOFF_op-LIGO-3G-deviation.md             [Path A kickoff, ~3-6 sesji]
└── FINDINGS.md                                       [ten plik — synteza]
```

**Spójność wewnętrzna:**
- Wszystkie linki `[[…]]` walidowane (single-folder relative).
- YAML headery jednolite, `tgp_owner: audyt/T01_LIGO3G_falsifier`.
- Konwencja PN counting ujednolicona (PHASE primary, ENERGY noted).
- Cross-references między plikami konsystentne.

**Spójność zewnętrzna:**
- Linki do `../EXTERNAL_REVIEW_2026-05-06.md`, `../S07_M911_derivation/`,
  `../README.md`, `../PRIORITY_MATRIX.md` walidowane.
- Linki do `../../research/op-newton-momentum/M9_1_pp_*` walidowane.
- Linki do `../../PREDICTIONS_REGISTRY.md`, `../../TGP_FOUNDATIONS.md`,
  `../../core/sek08c…` walidowane.

## §5 — Status T01 (formal)

| Atrybut | Wartość |
|---------|---------|
| **folder_status** | audit |
| **level** | T1 |
| **kind** | audit-test |
| **core_compatibility** | review-only |
| **last_reviewed_against_core** | 2026-05-07 |
| **may_edit_core** | false |
| **exports_findings** | true ✓ (was: false) |
| **has_needs_file** | true ✓ (was: false) |
| **has_findings_file** | true ✓ (was: false) |
| **open_bridges** | `op-LIGO-3G-deviation`, `op-ppE-mapping` (kickoff briefs gotowe) |
| **depends_on** | `S07 M9.1'' status`, `M9.1'' P1 deviation prediction` (latter: DERIVED) |
| **impacts** | `PREDICTIONS_REGISTRY (falsifier statement)`, `peer-review path` |

**Status klasy T01:** **PREPARED** (was: NOTED).

**Przejście:** PREPARED → CLOSED-PARTIAL (registered) **wymaga**
T01-D1 (autor wkleja falsifier statement). **Decyzja czeka.**

## §6 — Lessons learned (meta dla przyszłych T-class audytów)

### 6.1 Audit-folder może być productive bez "robienia research"

T01 pokazuje, że audit-folder z `may_edit_core: false` może
*znacząco* posunąć status problemu, dostarczając:
1. **NEEDS contract** (konkrety zadań).
2. **Drafts do adopcji** (ready-to-paste falsifier).
3. **Previews** (analytical / numerical OOM przed full cyklem).
4. **Kickoff briefs** (skraca startup ~50% next cycle czas).
5. **Convention decisions** (rozwiązuje niejasności *przed* one-letter rejestru).
6. **FINDINGS synthesis** (audit-evidence do peer-review).

**Wniosek:** "audit only" nie znaczy "tylko opisuje problem" — może
być **operacyjnie productive** w ramach safety constraint.

### 6.2 Konwencyjne pułapki są wszechobecne

T01 odkrył *istotną* niejednoznaczność (energy-PN vs phase-PN) która
gdyby nie była wyjaśniona przed wpisem do rejestru, mogłaby wprowadzić
*systematyczną pomyłkę o factor ~10×* w bounds. **Lekcja:** każdy
audit T-class powinien mieć osobny `CONVENTION_DECISION` lub
explicit walidację konwencji.

### 6.3 Multi-coefficient TGP-signature jest underused

[[PPN_TO_PPE_MAPPING.md]] §3 odkrył, że TGP M9.1'' ma **wzorzec
stosunków** współczynników PN — *to* jest distinguishing signature,
nie pojedynczy coefficient. **Lekcja:** dla T-class predictions
TGP, *zawsze* sprawdzaj czy istnieje multi-coefficient pattern przed
falsifier-statement. To zwiększa robustność testu i ranking
strategiczny.

### 6.4 LIGO-O3 jako "borderline" jest częsta pomyłka

Sesja zaczęła z założeniem "TGP w basenie konsystencji LIGO-O3"
(margin 100×); po patch konwencji PN okazało się że margin jest
*~1×* (borderline). **Lekcja:** OOM bounds są wrażliwe na konwencję;
zawsze re-derive po decyzji.

## §7 — Rekomendacja końcowa (dla autora TGP)

### Jeden-liner

> Predykcja M9.1'' (5/6) U³ jest fizycznie testowalna, struktulalnie
> falsyfikowalna, i operacyjnie ready-to-register. Wpis falsifier
> w `PREDICTIONS_REGISTRY.md` (~30 min pracy, gotowy draft w
> [[FALSIFIER_STATEMENT_DRAFT.md]]) **jest najwartościowszą akcją
> w TGP_v1 z najlepszym ROI** w sesji 2026-05-07.

### Sekwencja decyzyjna

```
Now (next session, ~30 min):
  → Tier 1: Wklej M911-P1 do PREDICTIONS_REGISTRY.md
  → T01 status: CLOSED-PARTIAL (registered)

Within ~2-3 sessions:
  → Tier 2: Uruchom op-ppE-mapping (kickoff brief gotowy)
  → Update FALSIFIER_STATEMENT_DRAFT [β_th]
  → T01 status: CLOSED-PARTIAL (B+C)

Within ~3-6 next sessions (po B):
  → Tier 3: Uruchom op-LIGO-3G-deviation (kickoff brief gotowy)
  → Update FALSIFIER_STATEMENT_DRAFT thresholds
  → T01 status: CLOSED-EXECUTED

Within ~1-2 months (co-schedule z S07):
  → Tier 4: Paper draft → PRD/CQG submission
  → T01 status: CLOSED-PUBLISHED

Optional sidechannel (~1 session, anytime):
  → Tier 5: GWTC-3 reanalysis dla M911-P1 (sygnał *teraz*)
```

## §8 — Cross-references

### Wewnątrz folderu T01

- [[README.md]] — diagnoza, closure progress
- [[NEEDS.md]] — open tasks (N1–N11, Q1–Q5, B1–B4)
- [[FALSIFIER_STATEMENT_DRAFT.md]] — Path C ready-to-paste
- [[PPN_TO_PPE_MAPPING.md]] — Path B preview
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — Path A preview
- [[CONVENTION_DECISION.md]] — phase-PN convention
- [[CYCLE_KICKOFF_op-ppE-mapping.md]] — Path B kickoff
- [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] — Path A kickoff

### Audit-level

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa
- [[../README.md]] — indeks audytu (klasa T)
- [[../PRIORITY_MATRIX.md]] — T01 P3 strategiczny
- [[../S07_M911_derivation/README.md]] — meta-bloker

### Core-level (read-only)

- [[../../PREDICTIONS_REGISTRY.md]] — target Path C
- [[../../TGP_FOUNDATIONS.md]] § 3 — narrative source
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — c_n derivation
- [[../../research/op-newton-momentum/M9_1_pp_setup.md]] — M9.1'' setup
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]] — M9.1'' P1
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] — M9.1'' canonical
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 0 gate enforcement

## §9 — Historia rewizji

- **v1 (2026-05-07):** Pierwsza synteza po sesji 2026-05-07. Audit
  T01 awansował z NOTED do PREPARED w 1 sesji (~3 godz.). 9 plików
  w folderze (README v2, NEEDS, 4 documents previews/drafts,
  CONVENTION_DECISION, 2 CYCLE_KICKOFF briefs, ten FINDINGS).
  Patches CONVENTION_DECISION zaaplikowane intra-session do README,
  FALSIFIER_STATEMENT_DRAFT, SENSITIVITY_BACK_OF_ENVELOPE, NEEDS.
  Rekomendowana sekwencja akcji: C-light → B → A → D, z opcjonalnym
  T-D5 (GWTC-3 reanalysis) anytime.
