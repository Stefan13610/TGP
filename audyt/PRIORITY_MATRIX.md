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
| L03 | V''(1)<0 + K(φ)=K_geo·φ⁴ | Ontologiczny | ~~P2~~ **EXECUTED 2026-05-06 via cykl L03** | research/op-L03-spectral-stability-2026-05-06 (synteza pre-existing 70% + mode counting Z₂ + tachyon check 4 profile + S-L thm:spectral-synthesis-L03; sek08b ssec:spectral-synthesis-L03 addytywne) | wszystkie soliton ODE |
| L04 | α=1 vs α=2, K=g² vs K=g⁴ | Ontologiczny | ~~P2~~ **RESOLVED 2026-05-04 via cykl L04** | research/op-L04-ODE-canonicalization-2026-05-04 (analytical decision: α=2 canonical via thm:D-uniqueness + Phase 2 universal mass formula + R5 ↔ Phase 2 IFF α=1) | LP-6 reframing pending |
| L05 | k=4 vs p=5-α | Ontologiczny | **P2** | mass-formula-reconcile | LP-4, R3, masa leptonów |
| L06 | m_X "locked" 100 MeV | Ontologiczny | **P2** | omega4-axion-mass | ω.2, ω.3, TT7-TT12 |
| D01 | Liczbowy dryft (α_s, m_H, Φ₀, …) | Liczbowy | ~~P3~~ **EXECUTED 2026-05-06 via cykl D01** | research/op-D01-anchor-lock-2026-05-06 (Lock manifest 5 anchors + Phase 2A 7 edycji + Phase 2B 4 papers_external footnotes + Phase 2C 6 deferred scripts; pre-existing B3-v2 14 lokacji α_s=0.1184/Φ_0=24.783) | wszystkie predykcje liczbowe (5/5 falsyfikatorów lockowanych) |
| M01 | Status creep w PREDICTIONS_REGISTRY | Metodologiczny | **P3** | registry-audit | counter 856 |
| M02 | 74394a8 forward-patch | Metodologiczny | **P3** | retrospect-rollback | counter 856, F6, M_TGP |
| **EXT-1** | **Kosmologia radiacyjna z varying c, ℏ, G → BBN/CMB** | rozszerzenie L01 | **P1 OTWARTE RYZYKO** (EXT-1 v2 2026-05-06) | research/op-FRW-radiation-era-varying-c (PROPOSED, najpilniejsze); plus op-BBN-TGP, op-CMB-TGP. Subiektywna ocena recenzenta: 55-65% szansa na zgodność BBN/CMB w 5%/0.5% tolerance przez varying c, ℏ, G. T^μ_μ_EM=0 trzyma się strukturalnie (Weyl-niezmienniczość 4D), ALE H_TGP(z) z varying constants może dać własną strukturę kinematyczną ekspansji. **Patrz [[L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]** | L01 NEEDS N1, N2 |
| **EXT-2 / L07** | **Warunek zerowej sumy ∫_Σ φ √h = 0 — aksjomat vs derywacja** | Ontologiczny (NEW class L07) | **P2** | op-zero-sum-derivation lub op-nonlocal-foundations. Ścieżki: A (Z₂-symetria substratu), B (Lagrange'a multiplier), C (φ_eff = φ − ⟨φ⟩_Σ), D (nielokalność spacelike). **Patrz [[L07_zero_sum_axiom/README.md]]** | sek05 ciemna energia (Λ_eff > 0), ax:zero w sek01 |
| **EXT-3 / S07** | **Metryka M9.1'' jako postulat vs derywacja z pierwszych zasad** | Strukturalny (NEW class S07) | **P2** | op-metric-from-substrate, op-metric-uniqueness, op-entropic-metric. M9.1'' = trzecia iteracja po falsyfikacji form I, II — empirycznie wybrana, NIE strukturalnie wyprowadzona. Ścieżki: A (coarse-graining H_Γ), B (wariacyjne kryterium minimalności), C (przyznanie statusu (P)), D (Verlinde/Padmanabhan entropic). **Patrz [[S07_M911_derivation/README.md]]** | sek08c, sek08a M9.1'' canonical, TGP_FOUNDATIONS § 2 status (E) |
| **EXT-4 / L08** | **Phase 6+ why_n3 (warstwa 3c: kinki jako fermiony) — domknięcie analityczne** | Ontologiczny (NEW class L08) | **P2** | op-why_n3-Phase6-dirac, op-why_n3-Phase6-mass-exponent, op-why_n3-Phase6-quarks. Brak: emergent Dirac propagator z antysymetrią, derywacja e² w wykładniku g_0^(e²/2), kwarki/neutrina/bozony w warstwie 3c, algebra Cliffordowska z Z₂ (niewystarczające, vs SU(N)_L×SU(N)_R Skyrme). **Patrz [[L08_kink_fermion_closure/README.md]]** | TGP_FOUNDATIONS § 4, L05 mass exponent drift |
| **EXT-5 / T01** | **\|Δg_tt\|=(5/6)U³ vs LIGO 3G (Einstein Telescope/Cosmic Explorer)** | NEW class T (testy falsyfikujące) | **P3 (strategiczny)** | op-LIGO-3G-deviation, op-ppE-mapping. Ścieżki: A (explicit Δh(f) calc), B (mapowanie ppE), C (falsifier statement w PREDICTIONS_REGISTRY), D (peer-review submit do PRD/CQG). Pierwsza pozycja w nowej klasie T. **Patrz [[T01_LIGO3G_falsifier/README.md]]** | M9.1'' P1, peer-review path |
| M03 | Brak balance sheet pre-74394a8 | Metodologiczny | ~~P4~~ ~~IN_PROGRESS~~ **🎉 EXECUTED ALL PHASES — 100% COMPLETE 2026-05-06** | research/op-M03-balance-sheet-retrofit-2026-05-06 (Framework 8 dokumentów + **40 retrofitów** + Phase5_registry_refactor_draft 12 cykli + Phase 6 gate enforcement w meta/). **Phase 1-4 100% COMPLETE** w 1 sesji 2026-05-06 (10 mini-sessions A-J): Phase 2 high-risk **11/11**, Phase 3 medium-risk **15/15**, Phase 4 low-risk **13/13** + 5 pre-M03 (chi.1, UV.2, λ.1, ω.2, ω.3). **Honest reporting baseline 33/40 = 82.5%** ★. **Mixing-operator family ISOLATED** do **4 cykli** (κ.1+ι.1+μ.1+**ν.1** sesja C — single coherent NUMEROLOGICAL/ANSATZ cascade). 9 pre-74394a8 systemic patterns CONFIRMED. **3 NULL closure pattern** (μ.1' pre-imp + ο.2 post-imp ABANDONED + op-eht-A negative-on-rescue). **3 SPLIT verdict pattern** (ν.1 + π.1 + ο.1). **4 multi-stage self-correction** (★★★ ψ.1 3-version canonical Phase 6 §4 + ο.2 + μ.1' + τ.3). **5 aggregate cycles** (M10 42/42 + M11 62/62 + Phase 1 50/50 + Phase 2 54/54 + Phase 3 60/60 GRAND TOTAL **281**). **5-cycle 2026-05-01 self-correction cluster** (ω.1+σ.1+ψ.1+τ.3+τ.2). **2 corrective upgrades** (ψ.1 v1→v2 + UV.3 DEPRECATES UV.2 K_struct TAUTOLOGY). **8 ★★ exemplary** + 1 ★★★ canonical (ψ.1). **1 research-track wieloletni** (op-uv-research OPEN no deadline). γ.1 resolves D01 NEEDS N2. Estimated post-M03: counter **856 → ~712 effective**, ratio **5.5 → 3.5-4.5**. Phase 6 ABSOLUTE BINDING gate enforced (CALIBRATION_PROTOCOL + CALIBRATION_GATE_ENFORCEMENT + AGENT_PROTOCOL §0). Future cycles MUST create Phase0_balance.md przed registry commit. **Patrz [[M03_balance_sheet_missing/POST_ACTION_FINAL_2026-05-06.md]]**. Phase 5 full implementation (per-row epistemic class tagging) deferred ~2-3 sesje. | 40 cykli audited (33 ★ + 7 ⚠) |

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
- **D01** → **EXECUTED 2026-05-06 via cykl D01** (research/op-D01-anchor-lock-2026-05-06):
  Lock manifest 5 anchors (α_s=0.1184, Φ_0=24.783, m_H=125.31/125.1 dwie ścieżki,
  Σm_ν=59.01, g_0^e=0.86941). Phase 2A: 7 edycji deterministicznych (README Σm_ν,
  sek09:1428, dodatekU:19, sek00:341, companion 520+1217). Phase 2B: 4 footnotes
  papers_external z B3-v2 canonical lock. Phase 2C deferred: 6 tooling scripts
  internal Φ_0=24.65 (self-consistency `K_geo·m_sp²=π·Φ_0²`). Pre-existing B3-v2
  14 lokacji α_s + Φ_0. Razem ~25 lokacji propagated, 5/5 falsyfikatorów TGP
  lockowanych. Outstanding (P3/P4): N1 tooling re-verification, N2 Brannen
  formal derivation, N3 Ω_Λ tracking, N4 predictivity ratio re-derivation,
  N5 m_H F11 vs CW reconciliation.

### Klaster D: ontologia (P2-P3)

- **L01** (ρ definition) — **EXECUTED 2026-05-04** via cykl L01 (formal `ρ ≡ -T^μ_μ/c_0²` + mapping SM 5 sektorów + photon treatment + sek08a addytywne edit). Bonus: zamyka S04 N1 (formal kowariantna derywacja `L_mat`).
- **L02** (β/γ semantyka) — **EXECUTED 2026-05-04** (renotacja w dodatekA).
- **L03** (V'' stability) — **EXECUTED 2026-05-06 via cykl L03** (research/op-L03-spectral-stability-2026-05-06): synteza pre-existing pokrycia (~70% już w sek08b 3-tier ghost-freedom + sek08_formalizm prop:vacuum-stability + sek08a prop:vacuum-stability-G0) z 3 nowymi elementami (mode counting Z₂, tachyon check 4 profile, S-L thm:spectral-synthesis-L03). NON-BREAKING addytywne sek08b ssec:spectral-synthesis-L03 (~140 linii LaTeX). Pozostałe N1-N5 z NEEDS.md są P3/P4 (rigorous bounds, 1-loop, extreme regimes — nie blokują predykcji).
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
