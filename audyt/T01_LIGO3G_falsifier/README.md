---
title: "T01 — Test różnicujący: |Δg_tt| = (5/6) U³ vs LIGO 3G (Einstein Telescope / Cosmic Explorer)"
date: 2026-05-06
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - audit
  - test-falsification
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - M911
  - 3PN
  - ppE
  - T01
  - EXT-5
  - new-class-T
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../S07_M911_derivation/README.md]]"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[CONVENTION_DECISION.md]]"
  - "[[CYCLE_KICKOFF_op-ppE-mapping.md]]"
  - "[[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]]"
  - "[[FINDINGS.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
tgp_status:
  folder_status: audit
  level: T1
  kind: audit-test
  core_compatibility: review-only
  last_reviewed_against_core: 2026-05-07
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: ["op-LIGO-3G-deviation", "op-ppE-mapping"]
  depends_on: ["S07 M9.1'' status", "M9.1'' P1 deviation prediction"]
  impacts: ["PREDICTIONS_REGISTRY (falsifier statement)", "peer-review path"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-07
revisions:
  - v1 (2026-05-06): Pierwsza wersja — diagnoza EXT-5, ścieżki A/B/C/D
  - v2 (2026-05-07): Dodano NEEDS.md, FALSIFIER_STATEMENT_DRAFT.md
    (Path C draft gotowy do adopcji), PPN_TO_PPE_MAPPING.md (Path B
    preview — analytical dictionary), SENSITIVITY_BACK_OF_ENVELOPE.md
    (Path A preview — OOM detektywności ET-D / CE). Audit-folder
    przygotowany do uruchomienia cykli `op-LIGO-3G-deviation/` i
    `op-ppE-mapping/`. T01 status: PREPARED (NEEDS contract gotowy,
    falsifier draft gotowy do wklejenia do PREDICTIONS_REGISTRY).
  - v2.1 (2026-05-07, kontynuacja sesji): Dodano CONVENTION_DECISION.md
    (PN counting: energy vs phase, PHASE primary adoptowana),
    CYCLE_KICKOFF_op-ppE-mapping.md (Path B brief uruchomieniowy,
    ~2-3 sesje), CYCLE_KICKOFF_op-LIGO-3G-deviation.md (Path A brief,
    ~3-6 sesji), FINDINGS.md (synteza audytu). Patches CONVENTION
    zaaplikowane do README problem #4, FALSIFIER_STATEMENT_DRAFT §1
    (b_ppE = -1, nie +1), SENSITIVITY_BACK_OF_ENVELOPE §4.2 (LIGO-O3
    bound ~10⁻¹), NEEDS N7. has_findings_file: true. Folder z 9
    plików, ~50KB, internally consistent. Status pozostaje PREPARED;
    przejście do CLOSED-PARTIAL (registered) wymaga T01-D1 (autor
    wkleja falsifier do PREDICTIONS_REGISTRY).
  - v3 (2026-05-07, sesja C-B-A-D execution): autor autoryzował
    pełną sekwencję C-light → B → A → D. **Wszystkie 4 ścieżki
    EXECUTED**:
    (C) M911-P1/P2/P3 wpisany do PREDICTIONS_REGISTRY.md Sector 2b +
    3 roadmap rows;
    (B) cykl `research/op-ppE-mapping/` utworzony i wykonany
    (Phase 0+1+2+3, sympy LOCK 14/14, β_ppE^TGP^(b=-1) = -5/64
    ≈ -7.81·10⁻² LOCKED, multi-coefficient ratios {-23/10, -38/23,
    +337/228} LOCKED);
    (A) cykl `research/op-LIGO-3G-deviation/` utworzony i wykonany
    (Phase 0+1+2+3, Fisher matrix Python script, calibrated
    thresholds: CE single-event decisive >5σ, ET-D stack ~10 events,
    LIGO-O5 stack ~16-250 events);
    (D) paper draft `papers/M911_LIGO3G_paper/paper_draft.md`
    gotowy do LaTeX porting (~700 lines, paper-ready).
    Patches do FALSIFIER_STATEMENT_DRAFT, NEEDS, SENSITIVITY,
    rejestru zaaplikowane. Side-cykl GWTC-3 reanalysis recommended
    (~1 sesja Python; pre-publication confirmation channel).
    **T01 status: CLOSED-EXECUTED.** SESSION_REPORT_2026-05-07_C-B-A-D.md
    dokumentuje pełną execution.
  - v3.1 (2026-05-07, Tier 5 GWTC-3 reanalysis EXECUTED): cykl
    `research/op-GWTC3-reanalysis/` utworzony i wykonany (Phase
    0+1+2+3, Bayes factor Python script). **Verdict: TGP M911-P1
    CONSISTENT z GWTC-3 within 1σ** (BF_TGP/GR ≈ 0.97 INCONCLUSIVE;
    deviation 0.37σ z observed centroid). Q3 z NEEDS CLOSED.
    **KEY DISCOVERY:** detekcja TGP wymaga **dedicated TGP-specific
    Bayes pipeline** (single-coefficient β_TGP prior); generic
    LIGO ToGR multi-coefficient marginalized analysis jest factor
    ~50× weaker. M911-P1 entry rejestru zaktualizowany z workflow
    caveat. Paper draft (Path D) wymaga §4-§5 revision per
    Phase3_verdict.md §3.3. **Recommended next:** Workflow A
    (TGP-specific Bayes z public GWTC-3 events, ~1 sesja Python,
    expected ~1-2σ tentative signal).
---

# T01 — Test różnicujący: |Δg_tt| = (5/6) U³ vs LIGO 3G

## Klasa: NOWA KLASA T (testy falsyfikujące, EXT-5) • Priorytet: **P3 (strategiczny)**

> **Pierwsza pozycja w nowej klasie T** (testy zewnętrzne) —
> rekomendowana przez EXT-5 jako uzupełnienie obecnych klas
> S/L/D/M (które są o spójności wewnętrznej).

## Diagnoza (z EXT-5)

`TGP_FOUNDATIONS.md` § 3 sygnalizuje konkretną predykcję różnicującą:

> M9.1'' przewiduje explicit |Δg_tt| = (5/6) U³ deviation od GR
> (testowalne LIGO 3G).

Ten jeden konkret odróżnia TGP od GR w obserwowalnym reżimie — i jako
taki jest **najmocniejszym istniejącym kandydatem na test
falsyfikujący**.

## Pliki dotknięte

- [[../../TGP_FOUNDATIONS.md]] § 3 (eksplicit deviation)
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
  (M9.1'' P1, derywacja deviation)
- [[../../research/op-newton-momentum/]] (M9.x PPN higher-order coefficients)
- [[../../PREDICTIONS_REGISTRY.md]] (status predykcji — bez falsifier statement)

## Problemy jakie rodzi

### 1. Brak explicit sensitivity analysis

TGP_v1 deklaruje "testowalne LIGO 3G", ale **nie ma w `research/`**
policzenia:
- jaki SNR w Einstein Telescope / Cosmic Explorer wymagany
- przy jakich masach BBH/BNS deviation jest detektowalna
- jaki strain h(f) odpowiada (5/6) U³ deviation
- parameter degeneracy z innymi modyfikacjami GR (ppE framework)

### 2. Brak wpisania do PREDICTIONS_REGISTRY z statusem FALSIFIABLE

Predykcja (5/6) U³ powinna być (W)/(P) z **explicit warunkiem
falsyfikacji**:
> "jeśli LIGO 3G nie zobaczy deviation o określonej amplitudzie,
> M9.1'' upada"

Bez tego jest *opisem konsekwencji*, nie *kontraktem z obserwacjami*.

### 3. Konkurencja modeli

ppE framework (parametrized post-Einstein, Yunes-Pretorius) ma
**5 parametrów modyfikujących GR**. TGP musi pokazać, że (5/6) U³
jest **specyficznym wyborem** w tej parametryzacji (które α_ppE,
β_ppE odpowiada), żeby być rozróżnialna od innych modyfikacji.

### 4. Reżim ważności

U³ deviation w `g_tt` to **2PN-phase** correction w stationary phase
approximation waveformu (równoważnie: 3PN-energy w konwencji PPN/Will;
zob. [[CONVENTION_DECISION.md]] dla pełnego dictionary i
rekomendacji konwencji).

LIGO obecnie ma 3.5PN-phase complete (Mishra et al. 2016), LIGO 3G
celuje w 4PN-phase (ET science case 2020). Dla 2PN-phase coefficient
(M911-P1 b_ppE = −1), aktualne LIGO-O3 bounds są ~10⁻¹; ET-D + CE
osiągną ~10⁻³ na single event, ~10⁻⁴ na stack ([[SENSITIVITY_BACK_OF_ENVELOPE.md]]).

**Wymaga spójnego rachunku 2PN-phase w TGP** — co znaczy: pełny
dE/dr + dE/dt + SPA chain w M9.1''. Pełne 2PN-phase matching to
dużo pracy ale wykonalne (cykl `op-ppE-mapping/` Phase 1; zob.
[[CYCLE_KICKOFF_op-ppE-mapping.md]]).

### 5. Dualizm modeli waveform

Czy w TGP używamy GR waveform z deviation (małe perturbacja) czy
*pełnego* TGP waveform (numerical relativity w TGP)? **Brak wyboru
jest blokujące.**

## Potencjalne ścieżki domknięcia

### Ścieżka A — explicit calc Δh(f) dla GW150914-like event

**Cykl proponowany:** `op-LIGO-3G-deviation/`. Użyć GW150914 jako
szablon (M_1+M_2, spinów, distance), policzyć GR waveform + (5/6) U³
deviation w fazie inspiral. Określić, jaki SNR potrzebny do detekcji.

**Wynik:** konkretna liczba dla Einstein Telescope (ET) i Cosmic
Explorer (CE).

### Ścieżka B — mapowanie na ppE

**Cykl proponowany:** `op-ppE-mapping/`. Pokazać, że (5/6) U³
deviation odpowiada konkretnym (α_ppE, β_ppE, b) z Yunes-Pretorius
parameterization. To umieszcza TGP w obecnym benchmarkingu modyfikacji
GR.

### Ścieżka C — falsifier statement w PREDICTIONS_REGISTRY

**Najlżejsza akcja:** dodać w `PREDICTIONS_REGISTRY.md` explicit wpis:

> **M9.1'' P1 — (5/6) U³ deviation.** Falsyfikacja: jeśli ET/CE z 5σ
> nie widzi deviation > X·10^(−Y) w fazie inspiral dla M_BBH > Z M_⊙,
> M9.1'' jest sfalsyfikowana. Odpowiedzialna grupa:
> `research/op-LIGO-3G/`.

### Ścieżka D — peer-review submit

**Najmocniejszy ruch:** napisać paper:

> "Strong-field test of M9.1'' metric: predictions for Einstein
> Telescope"

i wysłać do **PRD** lub **Class. Quantum Grav.** To wprowadza TGP
do mainstream falsification machinery.

## Rekomendowany priorytet

**P3 — średni, ale strategicznie ważny.** Pozostałe 4 punkty
(EXT-1..4) są **strukturalne** (długi teoretyczne). Ten jest
**strategiczny** — dotyczy *sposobu pokazania światu*, że TGP jest
falsyfikowalna.

## Powiązanie z istniejącym audytem

**Brak odpowiednika w S/L/D/M** (te są o spójności wewnętrznej, nie
o testach zewnętrznych).

EXT-5 proponuje **klasę T01** — "test różnicujący" jako pierwsza
pozycja w **nowej klasie T** (testy falsyfikujące).

**Powiązane ze [[../S07_M911_derivation/README.md]]** — bez derywacji
M9.1'' z fundamentu (S07), test (5/6) U³ jest *testem ansatzu*, a nie
*testem fundamentu*. Po zamknięciu S07 (M9.1'' DERIVED), T01 staje się
testem **derywowanej** predykcji.

## Postęp domknięcia (2026-05-07)

Status T01 awansował z **NOTED** (pierwsze v1 README) do **PREPARED**:
audit-folder dostarcza teraz wszystkie *przygotowawcze artifacts*
potrzebne, by uruchomić cykle `research/op-LIGO-3G-deviation/` i
`research/op-ppE-mapping/`. Sam wpis falsifier do
`PREDICTIONS_REGISTRY.md` (Path C, "najlżejsze domknięcie") jest
gotowy do *wklejenia* przez autora bez modyfikacji.

### Mapa artifacts → ścieżek z README

| Artifact | Plik | Path | Stan domknięcia |
|----------|------|------|-----------------|
| **NEEDS contract** (11 luk + 5 pytań Q + 4 blokery) | [[NEEDS.md]] | wszystkie A/B/C/D | **N1, N2, N8, Q4 CLOSED 2026-05-07; pozostałe pending Path A** |
| **Falsifier statement draft** | [[FALSIFIER_STATEMENT_DRAFT.md]] | C | **EXECUTED 2026-05-07** — wpisany do PREDICTIONS_REGISTRY.md jako M911-P1/P2/P3; lock z B Phase 1 zaaplikowany (β_ppE^TGP^(b=-1) = -5/64 central) |
| **ppE mapping preview** | [[PPN_TO_PPE_MAPPING.md]] | B | **SUPERSEDED przez** [[../../research/op-ppE-mapping/Phase1_results.md]] (locked 2026-05-07) |
| **Sensitivity OOM → LOCKED** | [[SENSITIVITY_BACK_OF_ENVELOPE.md]] | A | **§4 zaktualizowany z Phase 1 lock** (TGP/bound ratio 0.78–780); pełen Fisher pending Path A |
| **PN convention decision** | [[CONVENTION_DECISION.md]] | meta-konwencyjny | EXECUTED — patches zaaplikowane |
| **Path B kickoff brief** | [[CYCLE_KICKOFF_op-ppE-mapping.md]] | B | **EXECUTED 2026-05-07** → cykl [[../../research/op-ppE-mapping/]] uruchomiony i zamknięty Phase 1+2+3 |
| **Path A kickoff brief** | [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] | A | KICKOFF-READY (cykl dependent na Path B Phase 1.4 — teraz odblokowany) |
| **Audit synthesis** | [[FINDINGS.md]] | meta-syntetyczny | EXECUTED 2026-05-07 v1 |
| **Tier 5 GWTC-3 reanalysis** | [[../../research/op-GWTC3-reanalysis/]] | side-channel pre-publication | **EXECUTED 2026-05-07** — TGP CONSISTENT z GWTC-3 within 1σ; BF ≈ 0.97 INCONCLUSIVE; key discovery: detection requires TGP-specific Bayes pipeline |

### Kluczowe konkluzje preparatorne

1. **Predykcja jest już DERIVED** w `research/op-newton-momentum/`
   ([[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] §3.2):
   wszystkie współczynniki deviation (5/6, 23/12, 19/6, 337/72) są
   policzone analitycznie z α=2 vacuum Φ-EOM, sympy LOCK 5/5 PASS.
   T01 NIE musi wyprowadzać deviation — musi **skonkretyzować ją
   do strain h(f) i zarejestrować jako falsifier**.

2. **Detekcja jest fizycznie wykonalna** ([[SENSITIVITY_BACK_OF_ENVELOPE.md]]
   §4.2; po patch konwencji PHASE z [[CONVENTION_DECISION.md]]):
   TGP M9.1'' β_ppE^TGP_(b=−1) ≈ 10⁻¹ (OOM) jest **borderline
   w aktualnym LIGO-O3** (bound ~10⁻¹, na granicy wykrywalności,
   nie w basenie 100×) i staje się **decisive** w ET-D + CE
   (bound ~10⁻³–10⁻⁴). Window of testability: **już od ~2027**
   (LIGO-O5 stack ~3·10⁻²) → 2035 (ET-D + CE single-event ~10⁻³).
   Side-cykl GWTC-3 reanalysis (~1 sesja, NEEDS Q3 / FINDINGS Tier 5)
   może dać sygnał *teraz* na public data.

3. **Multi-coefficient pattern jako TGP-specific signature**
   ([[PPN_TO_PPE_MAPPING.md]] §3): pojedynczy 3PN deviation b_ppE = −1
   jest dzielony z dCS, sGB, Einstein–Æther. Co czyni TGP
   distinguishable: **wzorzec stosunków (β_3PN, β_4PN, β_5PN)**
   wymuszony przez α=2 i hyperbolic f(ψ) — niereprodukowalny
   z innych modyfikacji bez fitting.

4. **Bloker meta = S07.** T01 jest *testem ansatzu* dopóki S07 nie
   zamknie M9.1'' z fundamentu (ścieżki A/B/D w
   [[../S07_M911_derivation/README.md]]). To **NIE blokuje** T01 —
   ansatz może być sfalsyfikowany niezależnie. Ale po zamknięciu
   S07 awans T01 do "test derywowanej predykcji" znacząco wzmacnia
   peer-review path D.

5. **Konwencja PN counting wymaga decyzji** ([[PPN_TO_PPE_MAPPING.md]]
   §1 KRYTYCZNA UWAGA): "U³ → 3PN" w README używa konwencji
   energetycznej; standardowa literatura inspiral nazywa tę samą
   deviation "2PN phase correction". Cykl `op-ppE-mapping/` musi
   zatwierdzić jedną konwencję dla rejestru.

### Następne kroki (kandydatne, dla autora)

| Akcja | Wymagane | Zwrot |
|-------|----------|-------|
| **(C-light)** Wkleić [[FALSIFIER_STATEMENT_DRAFT.md]] §1, §2 do `PREDICTIONS_REGISTRY.md` z `[β_th]` placeholder | autorska decyzja, ~30 minut | T01 status: **CLOSED-PARTIAL (registered)**. Public falsifier-contract dla TGP. Counter rejestru +1 (~857). |
| **(B)** Uruchomić cykl `research/op-ppE-mapping/` | new cycle, ~2–3 sesje analytic | β_ppE^TGP zamknięty liczbowo. T01 N2 closed. |
| **(A)** Uruchomić cykl `research/op-LIGO-3G-deviation/` | new cycle, ~3–6 sesji numerical (Fisher) | SNR thresholds + degeneracy. T01 N1, N3, N4, N5, N7, N11 closed. |
| **(D)** Draftować short paper "M9.1'' 3PN deviation: predictions for ET/CE" | new research output, ~1–2 mies. | Peer-review path. T01 N9 closed. |
| **(meta)** Zamknąć S07 ścieżką A lub B | EXT-3 cykl | T01 awansuje z "test ansatzu" do "test derywowanej predykcji". B1 closed. |

**Rekomendacja sekwencji:** **C-light → B → A → D**. C-light jest
darmowy; B jest analitycznie tani; A wymaga more work ale daje
*decisive* contract; D pójdzie po A naturalnie.

## Cross-references

### Wewnątrz folderu T01

- [[NEEDS.md]] — concrete open tasks (11 N + 5 Q + 4 blokery)
- [[FALSIFIER_STATEMENT_DRAFT.md]] — Path C, gotowy draft do `PREDICTIONS_REGISTRY.md`
- [[PPN_TO_PPE_MAPPING.md]] — Path B preview (analytical dictionary)
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — Path A preview (OOM detectability)
- [[CONVENTION_DECISION.md]] — PN counting convention (PHASE adopted) + patches log
- [[CYCLE_KICKOFF_op-ppE-mapping.md]] — Path B cycle kickoff brief
- [[CYCLE_KICKOFF_op-LIGO-3G-deviation.md]] — Path A cycle kickoff brief
- [[FINDINGS.md]] — top-level audit synthesis + decision tiers dla autora

### Audit-level

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa
- [[../README.md]] — indeks audytu (klasa T dodana)
- [[../PRIORITY_MATRIX.md]] — T01 P3 strategiczny
- [[../S07_M911_derivation/README.md]] — meta-bloker (test ansatzu vs test fundamentu)

### Core-level (read-only z perspektywy T01)

- [[../../PREDICTIONS_REGISTRY.md]] — target dla falsifier statement (Path C)
- [[../../TGP_FOUNDATIONS.md]] § 3 — narrative source predykcji (5/6) U³
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
  M9.1'' P1 derivation
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — analytical c_n derivation, sympy 5/5
- [[../../research/op-newton-momentum/M9_1_pp_setup.md]] — M9.1'' setup, c²(ψ) modyfikacja
- [[../../research/op-newton-momentum/]] — pełen folder M9.x PPN
