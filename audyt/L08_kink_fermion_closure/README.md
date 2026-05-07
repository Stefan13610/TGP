---
title: "L08 — Phase 6+ why_n3 (warstwa 3c: kinki jako fermiony)"
date: 2026-05-06
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L08_kink_fermion_closure
tags:
  - audit
  - ontology
  - kink-fermions
  - emergent-dirac
  - spin-statistics
  - L08
  - EXT-4
  - why_n3
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../L05_mass_exponent_drift/README.md]]"
  - "[[../../research/why_n3/]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-why_n3-Phase6-dirac", "op-why_n3-Phase6-mass-exponent", "op-why_n3-Phase6-quarks"]
  depends_on: ["why_n3 Phase 1-5 closed 2026-05-01"]
  impacts: ["TGP_FOUNDATIONS § 4 warstwa 3c", "L05 mass exponent drift"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-4"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# L08 — Phase 6+ why_n3 (warstwa 3c: kinki jako fermiony)

## Klasa: LUKA ONTOLOGICZNA (zewnętrzna recenzja EXT-4) • Priorytet: **P2**

## Diagnoza (z EXT-4)

`TGP_FOUNDATIONS.md` § 4 deklaruje warstwę 3c jako **hipotezę**:

> **3c — Kinki / defekty** (cząstka = radialny kink Φ + topologia
> chiralna). Hipoteza/roadmap alternatywnego opisu fermionów jako
> struktur w samym Φ; emergent Dirac propagator. **Strukturalny szkic
> 5-fazowy zamknięty 2026-05-01** w `research/why_n3/`. Otwarte
> (Phase 6+): X = e²/4 analytic derivation, A^(5−α) vs A²·g_0^(e²/2)
> reconciliation dla τ/e.

Phase 1–5 są zamknięte strukturalnie. **Phase 6+ pozostaje otwarta**,
i to ona zawiera **kluczowe domknięcia analytyczne**.

## Pliki dotknięte

- [[../../research/why_n3/]] (Phase 1–5 closed; Phase 6+ open)
- [[../../TGP_FOUNDATIONS.md]] § 4 (3c hipoteza)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] lin. 9658+
  (energie kinku K_n(α))
- [[../L05_mass_exponent_drift/README.md]] (k=4 LP-4 vs p=5−α R3
  2026-05-01, P2 OPEN)

## Problemy jakie rodzi

### 1. Spin-statistics theorem

Poprawny opis fermionów wymaga, by stany dwucząstkowe były
antysymetryczne (exchange phase = −1). Phase 3 `why_n3` używa
"RP² + Berry phase π" — to znana ścieżka (Finkelstein-Rubinstein 1968),
ale wymaga, by Berry phase **spójnie indukowała Lorentzowską strukturę
spinową**.

**Bez explicit konstrukcji emergentnego propagatora Diraca z odpowiednimi
antykomutacyjnymi własnościami, "kink jako fermion" pozostaje
*roszczeniem strukturalnym*, nie *konstrukcją operacyjną*.**

### 2. Trzy generacje (e/μ/τ)

Phase 5 closed strukturalnie ma "uniwersalną formułę masy":

```
m_obs(g_0, α) = c_M · A_tail² · g_0^(e²(1−α/4))
```

(`research/why_n3/Phase5_*.md`). Dla α=2 daje g_0^(e²/2) ≈ g_0^3.69.
Reprodukuje m_μ/m_e i m_τ/m_e z PDG <0.01%.

**Ale e² w wykładniku jest empirycznym dopasowaniem** — bez derywacji
wykładnika z głębszej struktury, formuła jest *spektakularnym
numerologicznym sukcesem*, nie wyprowadzeniem. L05 (mass exponent
drift) jasno to sygnalizuje (k=4 LP-4 vs p=5−α R3 — niezgodność
wykładników między różnymi sektorami).

### 3. Kwarki, neutrina, bozony cechowania

TGP_v1 koncentruje się **na leptonach**. Kwarki (g_0 ∈ [0.817, 0.891])
są "uniwersalne" via ten sam ODE substratowy (`sek08b:529`), ale
**explicit predykcje mas kwarków NIE są w PREDICTIONS_REGISTRY**.
Neutrina (Σm_ν) są w D01 jako anchor lock, nie jako derywacja.
Bozony cechowania (W, Z, gluon) **nie mają realizacji w warstwie 3c**.

### 4. Algebra Diraca

Pola fermionowe w SM mają strukturę Cliffordowską {γ^μ, γ^ν} = 2g^μν.
Z kinka skalarnego Φ z Z₂ wyprowadzić tę algebrę jest nietrywialne.
Skyrme model (1962) jest klasycznym precedensem (skyrmion = baryon
w QCD-like chiralnej teorii), ale wymaga grupy chiralnej
SU(N)_L × SU(N)_R, **nie samej Z₂**. TGP ma tylko Z₂ — to *za mało*
dla pełnej algebry spinowej.

### 5. Zamiast Diraca: czy emergentna SUSY?

Niektóre programy (np. Wen-Levin string-net) generują emergentne
fermiony przez string-net condensation. To wymaga dyskretnej
geometrii i deconfinement. TGP ma graf Γ, więc droga w zasadzie
otwarta — ale niewystarczalność Z₂ pozostaje.

## Potencjalne ścieżki domknięcia

### Ścieżka A — explicit emergent Dirac propagator

**Cykl proponowany:** `op-why_n3-Phase6-dirac/` z konstrukcją
2-cząstkowego stanu fermionowego z dwóch kinków:
- weryfikacja antysymetrii pod exchange
- identyfikacja γ^μ z transportu Berry'ego po RP²
- spin-1/2 ma być *konsekwencją*, nie *postulatem*

### Ścieżka B — derywacja e² w wykładniku masy

**Cykl proponowany:** `op-why_n3-Phase6-mass-exponent/`.

**Hipoteza:** e² w g_0^(e²/2) wynika z ilości stopni swobody
w ogonie oscylacyjnym (2 polaryzacje × Berry phase 2π). Wymaga explicit
obliczenia: dlaczego e² (ładunek), a nie g (sprzężenie grawitacyjne)?

### Ścieżka C — przyznanie statusu hipotezy w prognozach

**Najszczersze rozwiązanie:** warstwa 3c pozostaje (H) *na zawsze*
w `TGP_FOUNDATIONS.md`, dopóki Phase 6+ nie zamknie się analitycznie.
TGP funkcjonuje jako **teoria grawitacji + materia jako warstwa 3a/3b**
(Dirac fields + g_eff coupling), bez roszczenia, że wszystko emerguje
z Φ.

To **odróżnia** TGP od programu unifikacyjnego — czyniąc ją **teorią
grawitacji emergentnej**, jak Verlinde 2010.

### Ścieżka D — rozszerzenie symetrii substratu

Z₂ jest niewystarczająca dla algebry spinowej. Rozszerzenie do SU(2)
lub Spin(3) na poziomie substratu **narusza § 1 FOUNDATIONS**
("dyskretna Z₂ chiralna na substracie") — to jest poważna decyzja
**pivot substratu** (zakazana bez rozmowy z autorem). Ale jeśli
Phase 6+ nie zamknie się z samą Z₂, rozważenie tej ścieżki będzie
konieczne.

## Rekomendowany priorytet

**P2 — wysoki.** Bez Phase 6+ TGP jest **teorią grawitacji emergentnej
+ container na SM**, nie *unifikacją*.

**Decyzja strategiczna:** czy TGP chce być kompletną teorią materii
(wymaga zamknięcia 3c), czy ograniczyć się do grawitacji
(3a/3b wystarczają)?

## Powiązanie z istniejącym audytem

Powiązane z [[../L05_mass_exponent_drift/README.md]] (k=4 vs p=5−α
drift między sektorami). L05 jest P2 OPEN.

EXT-4 proponuje **skonsolidować pod nową klasą L08** — "warstwa 3c
(kinki jako fermiony) — domknięcie analityczne". L05 + L08 razem
= status warstwy 3c kompletny audit.

## Cross-references

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-4 — recenzja źródłowa
- [[../README.md]] — indeks audytu
- [[../PRIORITY_MATRIX.md]] — do update z L08 P2
- [[../L05_mass_exponent_drift/README.md]] — powiązane (k vs p drift)
- [[../../research/why_n3/]] — Phase 1-5 closed, Phase 6+ open
- [[../../TGP_FOUNDATIONS.md]] § 4 warstwa 3c hipoteza
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] lin. 9658+ kinki
