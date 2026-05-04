---
title: "Priority Matrix — audyt TGP_v1 2026-05-04"
date: 2026-05-04
parent: "[[README.md]]"
type: audit-matrix
tgp_owner: audyt
tags:
  - audit
  - priority
  - matrix
  - planning
related:
  - "[[README.md]]"
  - "[[SUMMARY_2026-05-04.md]]"
---

# Priority Matrix — audyt TGP_v1

## Mapowanie problemów → priorytety → cykle naprawcze

| ID | Problem (skrót) | Klasa | Priorytet | Cykl naprawczy | Dotknięte LOCKED |
|----|-----------------|-------|-----------|----------------|-------------------|
| S01 | 4 sprzeczne formy metryki w sek08c | Strukturalny | ~~P1~~ **CLOSED-RESOLVED via G.0 (2026-05-02)** | ~~sek08c-rewrite~~ → cosmetic literal cleanup (Phase 5 G.0, optional) | M9 PPN, M9.x, sek08a κ |
| S02 | √(-g) niespójny z M9.1'' | Strukturalny | ~~P1~~ **CLOSED-RESOLVED via G.0 P32+P24** | ~~B6 M9.x re-run~~ (wykonane) → propagacja research/M9.x (~10 plików, opcjonalna) | M9.2 m_field, M9.3 quadrupole |
| S03 | β_PPN: 1/2 (metric) vs 1 (master) | Strukturalny | ~~P1~~ **CLOSED-RESOLVED via G.0 P23** | ~~PPN-convention-lock~~ (wykonane) | M9.1'' γ=β=1 |
| S04 | ax:metric-coupling ⊥ L_mat | Strukturalny | ~~P1~~ **CLOSED-DERIVED 2026-05-04** (B9 phenomenological + L01 formal cykl: Option-2 audytu promowane do DERIVED przez perturbation theory wokół ψ=1) | piąta siła, MICROSCOPE, L_mat formal |
| S05 | σ_ab łamie single-Φ axiom | Strukturalny | ~~P1~~ **CLOSED-RESOLVED via Path B PRIMARY (2026-04-26)** | ~~OP-7 axiom-decision~~ (wykonane) | OP-7 T2-T6, GW prediction |
| S06 | χ.1/UV.2 cyrkularność | Strukturalny | ~~P1~~ **CLOSED substantially via CRITIQUE 2026-05-02 + AUDIT_omega2/3 2026-05-04** | F6 rollback decyzja autora; counter reconciliation | F6, M_TGP, f_a (ω.3) |
| L01 | ρ operacyjna definicja | Ontologiczny | ~~P2~~ **EXECUTED 2026-05-04 via cykl L01** | research/op-L01-rho-stress-energy-bridge-2026-05-04 (formal `ρ ≡ -T^μ_μ/c_0²` + mapping SM 5 sektorów + photon treatment + sek08a addytywne) | L_mat, MICROSCOPE, S04 N1 closure |
| L02 | β/γ WF FP vs faza złamana | Ontologiczny | ~~P2~~ **EXECUTED 2026-05-04** | dodatekA `app:A-beta-gamma-distinction` dodana | M3-M8, LK-1d |
| L03 | V''(1)<0 + K(φ)=K_geo·φ⁴ | Ontologiczny | **P2** | full-spectral-analysis | wszystkie soliton ODE |
| L04 | α=1 vs α=2, K=g² vs K=g⁴ | Ontologiczny | ~~P2~~ **RESOLVED 2026-05-04 via cykl L04** | research/op-L04-ODE-canonicalization-2026-05-04 (analytical decision: α=2 canonical via thm:D-uniqueness + Phase 2 universal mass formula + R5 ↔ Phase 2 IFF α=1) | LP-6 reframing pending |
| L05 | k=4 vs p=5-α | Ontologiczny | **P2** | mass-formula-reconcile | LP-4, R3, masa leptonów |
| L06 | m_X "locked" 100 MeV | Ontologiczny | **P2** | omega4-axion-mass | ω.2, ω.3, TT7-TT12 |
| D01 | Liczbowy dryft (α_s, m_H, Φ₀, …) | Liczbowy | **P3** | C10-v2 + B3-v2 | wszystkie predykcje liczbowe |
| M01 | Status creep w PREDICTIONS_REGISTRY | Metodologiczny | **P3** | registry-audit | counter 856 |
| M02 | 74394a8 forward-patch | Metodologiczny | **P3** | retrospect-rollback | counter 856, F6, M_TGP |
| M03 | Brak balance sheet pre-74394a8 | Metodologiczny | **P4** | full-balance-sheet | 27+ cykli DERIVED |

## Klastry naprawcze (1 cykl = 1 problem zamknięty strukturalnie)

### ⚠ UPDATE 2026-05-04: Klaster A strukturalnie zamknięty przez G.0

**Re-priorytetyzacja po inspekcji aktualnego stanu plików rdzenia
2026-05-04:** cykl
[[../research/op-g0-r3-from-canonical-projection]] zamknął PHASE 4
dnia 2026-05-02 i **strukturalnie rozwiązał A1+A2+A3** (= S01+S02+S03):

- sek08a v2.0 ADDENDUM (sssec:g0-closure-v2): V_M911 LOCK + R3 ODE +
  vacuum stability fix + κ correction + Newton-limit (4/5)πG_0
- sek08c preamble G.0 CLOSURE block + 3 inline STATUS CLOSED markery
- 7 plików core LaTeX zmodyfikowanych, pdflatex compile clean (537 stron)

| Test G.0 | Result |
|----------|--------|
| Phase 1 G0a (volume integration) | 4/4 PASS — V_M911 LOCK reprodukuje R3 ODE |
| Phase 2 P21 (vacuum uniqueness) | 4/5 PASS + bonus (sek08a tachion bug) |
| Phase 2 P22 (mass spectrum) | 5/5 PASS — m_μ/m_e -0.0013%, m_τ/m_e +0.0049% |
| Phase 2 P23 (PPN γ=β=1) | 5/5 PASS — INVARIANT pod V update |
| Phase 2 P24 (FRW κ) | 5/5 PASS — form invariant po re-fit |
| Phase 3 P32 (Newton limit) | 5/5 PASS — q·c²/Φ_0 = (4/5)πG_0 |

Patrz POST_ACTION_UPDATE w S01/S02/S03 folderach.

**Co opcjonalnie pozostaje (Phase 5 G.0, niska priorytet):**
- Cosmetic literal cleanup form (I)-(III) w body sek08c
- Propagacja G.0 do ~10 plików `research/op-newton-momentum/`,
  `research/nbody/examples/`

### Klaster A (deprecated po G.0)

~~Trzy problemy o wspólnym podłożu — wszystkie biorą się z M9.1 → M9.1''
pivotu, który nie został doprowadzony do końca w plikach źródłowych.~~

~~**Rekomendowana sekwencja:** sek08c-rewrite → M9.x re-run → β_PPN
formal derivation. Estymata: 2–3 tygodnie pracy fizyka teoretyka +
review.~~

**Status post-G.0:** wszystkie 3 zamknięte strukturalnie + numerycznie.
Klaster A nie jest już aktywny dla rozwoju.

### ⚠ UPDATE 2026-05-04 (drugi pass): Klaster B również zamknięty

Druga inspekcja w sesji 2026-05-04 wykryła:

- **S04** zamknięte fenomenologicznie przez **B9 WEP MICROSCOPE composition test**
  (2026-05-01, 6/6 PASS): η_TGP_lab = 1.32×10⁻²⁶ vs MICROSCOPE 1.1×10⁻¹⁵
  → **8.3×10¹⁰× safe**. Plus closure_2026-04-26 T-α threshold (5/5 PASS)
  daje dalsze 4×10¹⁶× margin upgrade. Patrz S04 POST_ACTION_UPDATE.
- **S05** zamknięte strukturalnie przez **closure_2026-04-26 Path B PRIMARY**
  (11/11 PASS): σ_ab = composite operator z heredity equation,
  M² = 2m_s² *derived* z box-of-product, ghost-free przez Gram-positivity,
  single-Φ axiom *strictly preserved*. FOUNDATIONS §2 warstwa 0
  zsynchronizowane. Patrz S05 POST_ACTION_UPDATE.

**Wniosek:** klastry A i B są strukturalnie zamknięte przez wcześniejsze
cykle (G.0 2026-05-02, B9 2026-05-01, closure_2026-04-26). Mój audit
z 2026-05-04 nie zaktualizował tego automatycznie — POST_ACTION_UPDATE
files dodane retrospektywnie.

### Klaster B (deprecated po POST_ACTION analysis)

~~Po zamknięciu klastra A przez G.0, **klaster B staje się najpoważniejszą
otwartą luką w TGP_v1**.~~

**Status post-2026-05-04 second pass:** klaster B zamknięty.

### Faktyczny najpoważniejszy P1 otwarty: S06

Pozostaje **S06 (cyrkularność χ.1/UV.2)** jako jedyny faktycznie otwarty
P1. SUBAGENT_AUDIT_74394a8 explicit acknowledges, decyzja użytkownika:
„brak rollbacku, forward-patch w przyszłej sesji" (2026-05-02).
**Sesja 2026-05-04 jest tą przyszłą sesją.**

Sesja przekierowuje fokus na S06, patrz
`research/op-S06-omega2-omega3-circularity-audit-2026-05-04/`.

### Klaster C: rejestr (P2-P3)

- **S06** + **M01** + **M02** → wszystkie dotyczą tego samego: rejestr
  ma 856 wpisów, z których część jest tautologiczna lub fitted.
  Naprawcze: rollback χ.1 F6, downgrade UV.2/ω.3, audyt status creep.
- **D01** → wymaga jednorazowego globalnego locka i propagacji przez
  15+ plików.

### Klaster D: ontologia (P2-P3)

- **L01** (ρ definition) — **EXECUTED 2026-05-04** via cykl L01 (formal `ρ ≡ -T^μ_μ/c_0²` + mapping SM 5 sektorów + photon treatment + sek08a addytywne edit). Bonus: zamyka S04 N1 (formal kowariantna derywacja `L_mat`).
- **L02** (β/γ semantyka) — **EXECUTED 2026-05-04** (renotacja w dodatekA).
- **L03** (V'' stability) — wymaga formalnej analizy spektralnej.
- **L04** (α dualism) — **RESOLVED 2026-05-04** (cykl op-L04-ODE-canonicalization). α=2 jest jedyną kanoniczną formulacją; α=1 to specjalny case Phase 2 universal mass formula. Distinction `m_obs ≠ M_full` (insight użytkownika 2026-05-01) wyjaśnia pozorny dualizm.
- **L05** (mass exponent) — **częściowo zaadresowane** przez L04: k=4 (LP-4) i p=3 (R3 α=2) oba są specjalnymi przypadkami Phase 2 dwuwykładnikowej formuły `m=c·A²·g₀^[e²(1−α/4)]`. Pełna formal derivation k(α, d) generalization wciąż OPEN.
- **L06** (m_X) — czeka na ω.4 cycle.

## Czas naprawy (estymata)

| Klaster | Estymata | Charakter |
|---------|----------|-----------|
| A (rdzeń metryki) | **2–3 tygodnie** | LaTeX rewrite + numerical re-run |
| B (aksjomatyka) | **4–6 tygodni** | formal physics work |
| C (rejestr) | **1–2 tygodnie** | audit + propagacja |
| D (ontologia) | **3–4 tygodnie** | mix formal + decyzje |

**Razem (P1+P2+P3):** ok. 10–15 tygodni pracy fizyka teoretyka, jeśli
P1 i P2 idą sekwencyjnie. Klastry A, B, C mogą iść równolegle (różne
domeny), klaster D zależy od decyzji o L04 i L05 (autor).

## Czego ten audit NIE rozwiązuje

- **Nie wprowadza** nowej fizyki — tylko inwentaryzuje.
- **Nie modyfikuje** plików rdzenia (tylko meta-rejestr w `audyt/`).
- **Nie waliduje** numerycznie predykcji TGP — to robi `scripts/`.
- **Nie rozstrzyga** decyzji autorskich (L04 dualism, S05 σ_ab kierunek).

## Cross-references

- [[SUMMARY_2026-05-04.md]] — pełny raport
- [[README.md]] — indeks folderu
- [[../meta/AUDYT_TGP_2026-05-01.md]] — meta-audit z 2026-05-01
- [[../meta/PLAN_DOMKNIECIA_MASTER.md]] — plan z kwietnia (zamknięty)
- [[../meta/CALIBRATION_PROTOCOL.md]] — protokół anchor lock
