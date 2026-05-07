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
  - T01
  - EXT-5
  - new-class-T
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../S07_M911_derivation/README.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
tgp_status:
  folder_status: audit
  level: T1
  kind: audit-test
  core_compatibility: review-only
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-LIGO-3G-deviation", "op-ppE-mapping"]
  depends_on: ["S07 M9.1'' status", "M9.1'' P1 deviation prediction"]
  impacts: ["PREDICTIONS_REGISTRY (falsifier statement)", "peer-review path"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
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

U³ to człon **3PN** (post-Newtonian rzędu 3). LIGO obecnie ma 3.5PN
waveformy. LIGO 3G będzie miał ~4PN+. **Wymaga spójnego rachunku
3PN w TGP** (czy TGP daje *wszystkie* 3PN współczynniki, czy tylko
jeden? Pełne 3PN matching to dużo pracy).

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

## Cross-references

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa
- [[../README.md]] — indeks audytu (do dodania klasa T)
- [[../PRIORITY_MATRIX.md]] — do update z T01 P3
- [[../S07_M911_derivation/README.md]] — powiązane (test ansatzu vs test fundamentu)
- [[../../PREDICTIONS_REGISTRY.md]] — target dla falsifier statement
- [[../../TGP_FOUNDATIONS.md]] § 3 explicit deviation
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
  M9.1'' P1
- [[../../research/op-newton-momentum/]] M9.x PPN
