---
title: "Sesja 3.5 — Hotspot Reality Check (audyt mojego audytu)"
date: 2026-05-03
type: audit
status: COMPLETED
session: S3.5
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
  - "[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
tags:
  - audit
  - self-correction
  - anti-overclaim
  - lesson-learned
  - session-3-5
---

# Sesja 3.5 — Hotspot Reality Check

> **Cel:** zweryfikować, czy 20 hotspotów wpisanych w `CORE_HOTSPOTS.md`
> w Sesji 3 (transkrypcja sekcji § A i § B audytu) **rzeczywiście są
> otwarte**. Trigger: pytanie człowieka 2026-05-03 — "sporo z tych rzeczy
> jest już pozamykane".
>
> **Wynik:** **20/20 hotspotów już zamkniętych** w samym audycie
> 2026-05-01 przez sekcje § J–§ AB. Mój `CORE_HOTSPOTS.md` z Sesji 3 był
> globalnie przekłamany jako OPEN.
>
> **Jest to dokładnie ten sam typ błędu, którego kategorycznie zakazuje
> [[AGENT_PROTOCOL.md]] §3** — "anti-overclaim". Tylko że tym razem
> dotyczy mnie, nie subagenta z incydentu `74394a8`. Lesson-learned
> wpisuję do AGENT_PROTOCOL §3 jako drugi case study.

---

## 1. Szok poznawczy

`AUDYT_TGP_2026-05-01.md` ma **31 paragrafów § A–§ AB** (linia 32–2509),
z których:

- **§ A–§ E** = diagnoza problemów (8 CRITICAL + 12 HIGH + 13 MEDIUM + 5 LOW + 5 cross-cycle = 43 items)
- **§ F–§ I** = rekomendacje + pliki referencyjne
- **§ J–§ AB** = **NAPRAWA i DOMKNIĘCIE każdego z 43 items, w tym samym dokumencie**

Sekcja § AB.5 (linia 2496) podaje **final tally**:

| Klaster | Items | Closed | OPEN |
|---|---:|---:|---:|
| A (CRITICAL) | 8 | **8/8** ✓ | 0 |
| B (HIGH) | 12 | **12/12** ✓ | 0 |
| C (MEDIUM) | 13 | 12/13 | 1 (C6) |
| D (LOW) | 5 | **5/5** ✓ | 0 |
| E (cross-cycle) | 5 | **5/5** ✓ | 0 |
| **TOTAL** | **43** | **43/43 = 100%** | (C6 only) |

C6 (τ.2-v2 line-of-sight integral CMB θ derivation) to jedyny zostawiony
"open" item, ale jest **MEDIUM**, nie hotspot. Plus 100% closure jest
osiągnięte przez markery `CLOSED 2026-05-01` w samym audycie.

**Mój CORE_HOTSPOTS.md w Sesji 3** miał:
- 8 CRITICAL — wszystkie OPEN
- 12 HIGH — 7 OPEN, 4 READY-FOR-INTAKE, 1 RESOLVED

**Faktyczny stan:**
- 8 CRITICAL — wszystkie zamknięte (różne typy closure)
- 12 HIGH — wszystkie zamknięte (różne typy closure)

---

## 2. Taksonomia poziomów domknięcia (nowa)

Poprzednia mapa statusów (`OPEN | READY-FOR-INTAKE | IN-INTAKE | PROMOTING | RESOLVED`)
nie wystarcza. Audyt zamknął hotspoty na **5 różnych poziomach**:

| Poziom closure | Co zostało zrobione | Czy wymaga jeszcze pracy nad core? |
|----------------|----------------------|--------------------------------------|
| **AUDIT-FULL** | Pełny code change w core LaTeX + sympy LOCK + propagacja przez wiele plików | **NIE** — zamknięte trwale |
| **AUDIT-STRUCTURAL** | Sympy LOCK / derivation w research/, częściowa integracja core. Praca formalna OK, finalna integracja core może być follow-up. | Może wymagać follow-up cycle (wskazane w audycie) |
| **AUDIT-ANNOTATION** | Markery inline (`STATUS`, `[WITHDRAWN]`, audit-aware blockquote) bez przepisania logiki core | Logiki core nie zmieniono; future revision opcjonalna |
| **AUDIT-CONVENTION** | Wybór konwencji acknowledged (np. A3 master formula) | Convention obowiązuje; further work opcjonalny |
| **AUDIT-ARCHITECTURAL** | Decyzja architektoniczna o zachowaniu axiom (A4 Option-2) | Architektura zachowana; future Option-1 ALT możliwe |

Status `RESOLVED` w starym sensie ("wynik wszedł do core jako primary content")
**nie pasuje** do większości zamknięć audytu. Audyt zrobił to, co dało się
zrobić w ramach jednej sesji (markery, sympy LOCK, registry annotacje), ale
**pełny rewrite core LaTeX** dla większości hotspotów A nadal jest
**legitymacyjnym follow-up**, tyle że już nie jako sprzeczność, tylko jako
"polishing pełnej integracji".

---

## 3. Per-hotspot werdykt (8 CRITICAL)

### H-A1 — 4 sprzeczne formy metryki w sek08c.tex

- **Domknięte przez:** § M (linie 651–737)
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:** Inline `STATUS` markery w sek08c.tex przy każdej z 4 form:
  - Forma (I) boxed power → `FALSIFIED 2026-04-25` (M9.1: β_PPN=4)
  - Forma (II) eksponencjalna → `OBSOLETE`
  - Forma (III) antypodyczna → `NIESPÓJNA wewnętrznie`
  - Forma (IV) M9.1'' hiperboliczna → `KANONICZNA wg FOUNDATIONS`
- **Cytat:** "→ A1 **CLOSED** (markery inline w sek08c.tex + status każdej z 4 form explicit)." (M.5)
- **Cross-deepening:** § AB (research/why_n3 5-phase Dirac propagator) — niezależny test M9.1'' przez R3↔M9.1'' horizon coincidence (`g₀_crit=1.874 ↔ ψ=4/3`)
- **Co JESZCZE może być zrobione (nie blocking):** pełny rewrite sek08c.tex usuwający 3 obsoletne formy (zostawiając tylko Form-IV jako primary derivation). Decyzja: deferred — annotation acknowledgment wystarcza dla audytu.

### H-A2 — `√(-g)` wrong volume element

- **Domknięte przez:** § M (initial structural) + **§ U (B6 full re-run, 6/6 PASS)**
- **Closure type:** `AUDIT-FULL` (po B6 § U)
- **Co zrobiono:** Pełny re-run M9.x z poprawnym `√(-g) = c₀·ψ/(4-3ψ)`:
  - Step 1: sympy LOCK Form-IV `det(g) = -c₀²·ψ²/(3ψ-4)²`
  - Step 2: PPN re-derivation — `β_PPN = 1 EXACT`, `γ_PPN = 1 EXACT`
  - Step 3: M9.2 `m_field` central value `3.48×10⁻²` → **`3.98×10⁻²`** (+14.2% shift)
  - Step 4: M9.3 dispersion + Peters-Mathews — vacuum UNCHANGED, NS interior +43% correction
  - Step 5: scalar/tensor mode ratio — `1.27×10⁻⁶` (LIGO-safe)
  - Step 6: discrepancy table — pełna mapa Form-I vs Form-IV
- **Cytat:** "→ **B6 FULL CLOSURE**: M9.x re-run z Form-IV `√(-g) = c·ψ/(4-3ψ)` zakończony, wszystkie inheritance markers updated" (U.6)
- **PASS-count:** 6/6 PASS (sympy + numpy + scipy)
- **Pliki:** `research/op-newton-momentum/B6_m9x_sqrtg_rerun.py` + `.txt` + `_results.md`
- **Note:** M9.1''-P3 PPN status upgrade: "conditional" → **"EXACT"**

### H-A3 — β_PPN convention dla M9.1''

- **Domknięte przez:** § M (convention choice)
- **Closure type:** `AUDIT-CONVENTION`
- **Co zrobiono:** OBOWIĄZUJĄCA KONWENCJA — "master formula z kinetic correction":
  `β = f''(1)/f'(1)² + 2c_2/f'(1) = 1/2 + 1/2 = 1` (z `c₂=-1` z α=2)
- **Cytat:** "→ A3 **CLOSED konwencyjnie** (master formula obowiązująca; algebraic dopasowanie acknowledged)." (M.5)
- **Wzmocnione przez:** § U Step 2 — Form-IV + c₂=-1 daje `β=1 EXACT`, **niezależnie** od convention choice; konfirmuje że master formula nie jest tylko algebraic fit
- **Honest acknowledgment:** convention is "algebraic dopasowanie" — niezależna derivacja `c₂=-1` z first principles to long-term cycle (closure_2026-04-26 c₂ derivation pending)

### H-A4 — `ax:metric-coupling` vs `L_mat = -(q/Φ₀)·φ·ρ`

- **Domknięte przez:** § N (linie 739–842)
- **Closure type:** `AUDIT-ARCHITECTURAL`
- **Co zrobiono:** Decyzja **OPTION 2 (preserve axiom)**:
  - Aksjomat `ax:metric-coupling` zachowany
  - `φ` w `eq:L-mat-unified` reinterpretowane jako derived consequence struktury `√(-g_eff) ∝ φ` + ontologicznej roli `Φ` (`ax:zrodlo`)
  - `ρ = T^μ_μ/c₀²` z stress-energy tensor w kanonicznym sprzęganiu metrycznym
  - **NIE** independent dilaton coupling Brans-Dicke type
- **Pliki:** `core/sek08_formalizm.tex` + `core/sek08a_akcja_zunifikowana.tex` audit-aware komentarze
- **Cytat:** "→ A4 **CLOSED architecturally** (Option 2 preserves axiom)" (N.5)
- **Konsekwencje:** Eöt-Wash / MICROSCOPE consistency intact; no fifth force; equivalence principle preserved
- **Future:** Option-1 ALT (modify axiom) explicit odrzucony — long-term reanalysis A4-v2 możliwy ale nie blocking

### H-A5 — `m_e_eff` τ.3 dimensional inconsistency

- **Domknięte przez:** § L (linie 540–648)
- **Closure type:** `AUDIT-STRUCTURAL` (numerical TT7-TT12 follow-up zaadresowany przez § O.6 + § T)
- **Co zrobiono:**
  - Multiplicative form: `m_e_eff = m_e^(0)·[1 + (α_g/Λ²)(∂lnX)²]`
  - `δω/ω = (α_g/Λ²)(∂lnX)²` **bez** podziału przez `m_e^(0)`
  - sympy LOCK confirmed: `delta omega/omega = alpha_g*dlnX**2/Lambda**2`, diff=0
  - Λ-scan detection gates: `Λ ≲ 100 MeV` → **~GeV scale** (~3 OOM upward)
  - 6 plików τ.3 patched (program.md, Phase1-3 results, Python skrypty, PREDICTIONS_REGISTRY)
- **Cytat:** "→ A5 **CLOSED structural** (formula multiplicative LOCKED, propagacja do wszystkich τ.3 plików + cross-channel registry annotacje)" (L.5)
- **Numerical follow-up:** § O.6 (B7 sympy LOCK two-regime) + **§ T (B7-v2 numerical re-derivation z edge geometry)** — TT7-TT12 numerical thresholds re-derived

### H-A6 — ψ.1.v1 "Z(x)F² → Δc/c" wave-function renormalization

- **Domknięte przez:** § K (linie 469–537)
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:**
  - YAML `PASS` → `INVALIDATED` w Phase1/2/3_results.md
  - Header `⛔ INVALIDATED` block w każdym z 3 plików
  - Replacement reference do Phase4/5/6 v2
  - Cross-link do AUDYT_TGP_2026-05-01 A6/A8
- **Cytat:** "→ A6 **CLOSED** (header markers + program.md note)" (K.4)
- **Note:** ψ.1.v2 (Phase 4+5+6) jest aktywna i poprawna — v1 jest historią

### H-A7 — ω.2 nie locks `m_X`

- **Domknięte przez:** wcześniejsza naprawa ω.3 (mentioned w K.4)
- **Closure type:** `AUDIT-FULL` (przez ω.3 cycle)
- **Cytat z K.4:** "Po naprawie A6+A8 i wcześniej A7 (przez ω.3): **3 z 8 CRITICAL CLOSED**"
- **Note:** `op-omega3-axion-decay-constant/` jest w kwarantannie 74394a8 — ale to jest sprawa **innego incydentu** (subagent over-claim na ω.2/ω.3 commit `74394a8`); audyt 2026-05-01 A7 closure przez ω.3 jest legitymny i predates incydent

### H-A8 — Sagnac SNR ≈ 3×10⁴ artifact

- **Domknięte przez:** § K (linie 505–516)
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:**
  - TT13 entry w `PREDICTIONS_REGISTRY.md:1085` → `[WITHDRAWN 2026-05-01]`
  - TT13-TT18 (6 entries) wszystkie `[WITHDRAWN]` z replacement pointers (TT19-TT23)
  - Header note w `Phase1_results.md` o niefizycznych konfiguracjach E+B
- **Cytat:** "→ A8 **CLOSED** (registry [WITHDRAWN] + header markers)" (K.4)

---

## 4. Per-hotspot werdykt (12 HIGH)

### H-B1 — β_g, α_g sign Adams positivity

- **Domknięte przez:** § Z (linie 2185+)
- **Closure type:** `AUDIT-FULL`
- **PASS-count:** ψ.1.v2 (Phase 4+6.T6.5) + τ.3 Phase 4 `phase4_tau3_adams_positivity.py` 5/5 PASS
- **Wynik:** `α_g > 0 STRICT, UV-independent`
- **Cytat:** "B1 (β_g, α_g sign — Adams positivity v2 robust) | **FULL CLOSED 2026-05-01 (§ Z)**" (O.7)

### H-B2 — DESI 2024 evolving DE tension

- **Domknięte przez:** § O.1
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:** DE1/DE2 `LOCKED` → `LIVE TENSION` + `⚠ B2-tension 2026-05-01` markery w `PREDICTIONS_REGISTRY.md:859-860`
- **Cytat:** "→ B2 **CLOSED annotation-only**" (O.1)

### H-B3 — α_s(M_Z) lock 0.1184

- **Domknięte przez:** § O.2 (initial) + **§ Y (B3-v2 full propagation)**
- **Closure type:** `AUDIT-FULL` (po § Y)
- **Co zrobiono:** Lock 0.1184 propagated through 14+ aktywnych lokalizacji w LaTeX core, papers, tooling
- **Cytat:** "B3 — α_s = 0.1184 propagation through LaTeX core, papers and tooling (CLOSED 2026-05-01)" (Y header)

### H-B4 — Σm_ν anchor 59.01 vs 59.6 meV

- **Domknięte przez:** § O.3
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:** Z1 entry `PREDICTIONS_REGISTRY.md:940` → anchor 59.01 meV LOCKED + explicit "59.6 meV był zeroth-order pre-bisection target"
- **Cytat:** "→ B4 **CLOSED annotation-only**" (O.3)

### H-B5 — wymiar `q` niespójny

- **Domknięte przez:** § X
- **Closure type:** `AUDIT-FULL`
- **Cytat:** "B5 — q-dimension reconciliation across LaTeX core (CLOSED 2026-05-01)" (X header)
- **Co zrobiono:** Pełna reconciliation across LaTeX core — single dimension propagated

### H-B6 — M9.x re-run z poprawnym √(-g)

- **Domknięte przez:** § U (już opisane w H-A2)
- **Closure type:** `AUDIT-FULL`
- **PASS-count:** 6/6 PASS

### H-B7 — (∂lnX)² Greens function explicit

- **Domknięte przez:** § O.6 (structural sympy LOCK) + **§ T (B7-v2 numerical)**
- **Closure type:** `AUDIT-FULL` (po § T)
- **Co zrobiono:**
  - sympy LOCK two-regime: heavy regime `(∂lnX)²_edge = J²/m_X²`; light regime `(∂lnX)²_bulk = J²L²/(16π²)`
  - KEY PHYSICS: default τ.3 lab params → heavy regime universal → bulk signal ZERO; tylko edge shell `1/m_X` kontrybuuje
  - § T: B7-v2 numerical re-derivation TT7-TT12 z geometric edge analysis
- **Konsekwencja:** TT7-TT12 numerical thresholds patched with realistic detector volume integration

### H-B8 — M9.1'' Lagrangean derivation

- **Domknięte przez:** § V
- **Closure type:** `AUDIT-FULL` (closure via "honesty acknowledgment")
- **PASS-count:** 5/5 PASS
- **Cytat:** "✅ **CLOSED 2026-05-01 — 5/5 PASS** (closure via honesty acknowledgment)" (V header)
- **Pliki:** `research/op-newton-momentum/B8_lagrangean_independence_check.py` + `.txt`

### H-B9 — M9.2 WEP MICROSCOPE composition test

- **Domknięte przez:** § W
- **Closure type:** `AUDIT-FULL`
- **PASS-count:** 6/6 PASS (sympy + numerical)
- **Cytat:** "✅ **CLOSED 2026-05-01 — 6/6 PASS**" (W header)

### H-B10 — forma minimalna PPN β=2

- **Domknięte przez:** § O.7 (closed via A1 — Form-IV M9.1'' kanoniczna)
- **Closure type:** `AUDIT-FULL` (inheritance)
- **Cytat:** "B10 (β=2 reconciliation) | **closed via A1** (M9.1'' kanoniczna form (IV))" (O.7)

### H-B11 — DESI w₀w_a hint propagation

- **Domknięte przez:** § O.4
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:** WW8/WW11 entries w `PREDICTIONS_REGISTRY.md` → `⚠ B11-DESI-caveat 2026-05-01` markery (PARTIAL candidate)
- **Cytat:** "→ B11 **CLOSED annotation-only**" (O.4)

### H-B12 — niefizyczne pole konfiguracje E+B

- **Domknięte przez:** § O.5
- **Closure type:** `AUDIT-ANNOTATION`
- **Co zrobiono:** TT7 entry `⚠ B12-field-caveat 2026-05-01` + realistic ELI-NP schedule documented (E=10¹³ V/m + B=30T)
- **Cytat:** "→ B12 **CLOSED annotation-only** (realistic schedule documented, SNR re-derivation B7-pending)" (O.5)
- **Note:** B7 SNR re-derivation już zamknięte w § T

---

## 5. Mapa "co naprawdę pozostaje do roboty" w core

Po zamknięciach audytu, **realne** OPEN punkty wobec core to:

| ID | Co | Severity | Source |
|----|-----|---------:|--------|
| **C6** | τ.2-v2 line-of-sight integral CMB θ derivation | MEDIUM | § AB.5; facilitated by L₅'_b basis from § AA |
| (long-term) | A4-v2 ax:metric-coupling vs L_mat (Option-1 alt) | optional | § N (Option-2 wybrany; Option-1 ALT not blocking) |
| (long-term) | NS-NS ringdown form-IV numerical | follow-up | U.5 |
| (long-term) | Closure_2026-04-26 c₂=-1 derivation z first principles | follow-up | U.5 |
| (long-term) | Strong-field α(ψ) regulator dla ψ → 4/3 | follow-up | U.5 (T-α może to dostarczyć) |
| (long-term) | ψ.1-v3 Phase 8 parity-odd β̃_g Adams | follow-up | § AB closing remarks |
| (long-term) | ω.3 super-light substrate | follow-up | § AB closing remarks |
| (long-term) | why_n3 Phase 6+ A^(5−α) reconciliation + analytic X=e²/4 derivation | follow-up | § AB.3 |

**To jest kompletnie inny rozkład** niż mój CORE_HOTSPOTS.md sugerował.
Zamiast 15 OPEN + 4 READY-FOR-INTAKE + 1 RESOLVED (=20 hotspotów):

- **0 hotspotów blocking core** (wszystkie 8A + 12B zamknięte)
- **1 MEDIUM** poza scope hotspotów (C6)
- **~7 long-term follow-ups** poza scope audytu (next-audit triggers)

---

## 6. Implikacje dla `meta/research/` workflow

### 6.1 `CORE_HOTSPOTS.md` — wymaga pełnej rewizji

Z 20 hotspotów `OPEN` należy zmienić wszystkie na odpowiednie `RESOLVED-AUDIT-{TYPE}`.
Plik staje się **historią closure**, nie listą todo.

### 6.2 `IMPACT_MATRIX.md` — krawędzie hotspot→folder są nieaktualne

§ 3 mojego IMPACT_MATRIX.md zakładała, że hotspoty są aktywne i propagują
do `needs-migration` folderów. Faktycznie:
- Większość krawędzi to **historyczne ścieżki audytu** (gdzie zaszła naprawa)
- Foldery (np. `research/op-newton-momentum/`) **nie są** `needs-migration` w sensie planu — ich M9.x results mają już `B6-CLOSED`/`B8-CLOSED`/`B9-CLOSED` markery

### 6.3 `CORE_CANDIDATES.md` — initial seed kandydatów wymaga przeglądu

Sugerowani kandydaci:
1. `closure_2026-04-26/sigma_ab_pathB` → core (H-B6) — **H-B6 ZAMKNIĘTE w § U niezależnie**
2. `closure_2026-04-26/f_psi_principle` → core (H-B8) — **H-B8 ZAMKNIĘTE w § V niezależnie**
3. `closure_2026-04-26/Lambda_from_Phi0` → vacuum catastrophe — **NOT in 43-item list**, niezależny strukturalny kandydat
4. `closure_2026-04-26/alpha_psi_threshold` → core (H-B9) — **H-B9 ZAMKNIĘTE w § W niezależnie**

Czyli `closure_2026-04-26/` 4 podzamknięcia są **niezależnie** od audytu
2026-05-01 — zamykają luki, ale audyt też je domknął **innymi ścieżkami**
(B6 przez `op-newton-momentum/B6_m9x_sqrtg_rerun.py`; B8 przez
`op-newton-momentum/B8_lagrangean_independence_check.py`; B9 przez
analogiczny skrypt). Dwa niezależne sources zamykające ten sam hotspot.

To jest **ZDROWY znak** (cross-validation), ale oznacza, że `CORE_CANDIDATES.md`
shortlist wymaga **innego framing'u** — nie "ten folder zamyka hotspot",
ale "ten folder dostarcza alternative formal track dla wyniku, który już
ma markery audytu".

### 6.4 Sesja 4 (klasyfikacja YAML) — zmiana priorytetów

Pierwotnie planowałem:
- Sesja 4 batch 1 = aktywne `op-*` 3-fazowe
- Sesja 4 batch 2 = topic folders
- ...
- needs-migration folders dostają `core_compatibility: stale | broken`

**Po Sesji 3.5:** mało folderów ma realną `stale | broken` compatibility.
Większość M9.x folderów ma już `B6-CLOSED`/`B8-CLOSED` markery → ich
`core_compatibility` to `current` lub `partial` (co najwyżej).

`needs-migration` jako kategoria może być niemal pusta — to byłaby **dobra
wiadomość** dla TGP, ale wymaga ręcznej weryfikacji per folder w Sesji 4.

---

## 7. Lesson-learned (do AGENT_PROTOCOL §3)

**Co poszło nie tak w Sesji 3:**

1. **Czytałem tylko § A i § B audytu**, ignorując § J–§ AB (paragrafy "Naprawa…").
   30-sekundowy skan nazw paragrafów wystarczył, by zidentyfikować problem
   (co właśnie zrobił człowiek, sprawdzając mój zapis).

2. **Założyłem, że audyt = lista TODO.** W rzeczywistości audyt 2026-05-01
   to **audit + closure trail** w jednym dokumencie. To nie jest "snapshot
   problemów", to jest "snapshot problemów + akt naprawy".

3. **Pominąłem self-reported `CLOSED` markery** w § A/§ B nawet tam, gdzie
   były explicite wpisane (B1, B3 — te akurat złapałem; B5, B9 — te
   miały `pending` w § O.7, ale potem CLOSED w § X/§ W których pominąłem).

4. **Nie zweryfikowałem `final tally`** w § S i § AB.5 — gdyby tak,
   pierwszą rzeczą, którą bym zobaczył, byłoby "TOTAL: 43/43 = 100%".

**Reguła wynikowa (do AGENT_PROTOCOL §3):**

> Każdy "audit-driven" wpis (hotspot, intake, candidate) wymaga **przeczytania
> całego dokumentu źródłowego, włącznie z final closure summary**, zanim
> jakikolwiek item zostanie wpisany jako OPEN. Nie wystarczy przeczytać
> sekcji diagnozy. **Audyty TGP zawierają zarówno diagnozę jak i naprawę
> w jednym pliku.**

Plus: **niech to będzie drugi case study** w AGENT_PROTOCOL §3, obok
incydentu `74394a8`. Pierwsze case study jest "subagent fabrykował PASSes";
drugie case study (Sesja 3 mnie) jest "agent wpisał OPEN tam gdzie audyt
sam zamknął". Symetryczne błędy w przeciwnych kierunkach.

---

## 8. Następne kroki

1. **Aktualizacja `CORE_HOTSPOTS.md`** — wszystkie 20 hotspotów → odpowiednie
   `RESOLVED-AUDIT-{TYPE}` z cytatem source.
2. **Aktualizacja `AGENT_PROTOCOL.md` §3** — dopisanie reguły z §7 + drugi
   case study.
3. **Aktualizacja `IMPACT_MATRIX.md` § 3** — krawędzie z hotspotów do
   folderów oznaczone jako "historical paths of closure", nie active
   `needs-migration`.
4. **Aktualizacja `CORE_CANDIDATES.md`** — initial seed kandydatów dostaje
   notatkę o cross-validation z audytem § U/V/W.
5. **Aktualizacja `meta/PLAN_RESEARCH_WORKFLOW_v1.md` §9.1** — dopisanie
   Sesji 3.5 jako gate'a + uzasadnienie post-faktum.

Po tych aktualizacjach: ruszamy Sesję 4 (klasyfikacja YAML) **z czystym
stanem hotspotów**, czyli z realistycznym oczekiwaniem, że bardzo mało
folderów ma `core_compatibility ∈ {stale, broken}`.
